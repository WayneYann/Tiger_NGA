module dns_box
  use precision
  use param
  implicit none
  
  include 'fftw3.f'
  
  ! Length of the domain
  real(WP) :: L
  real(WP) :: dx
  character(len=str_medium) :: data_type
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  
  ! Chemistry parameters
  real(WP) :: density_ratio
  
end module dns_box

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine dns_box_grid
  use dns_box
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('Domain size',L)
  ny = nx
  nz = nx
  
  ! Set the periodicity
  xper = 1
  yper = 1
  zper = 1
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  dx = L/nx
  do i=1,nx+1
     x(i) = (i-1)*L/nx
  end do
  do j=1,ny+1
     y(j) = (j-1)*L/ny - 0.5_WP*L
  end do
  do k=1,nz+1
     z(k) = (k-1)*L/nz - 0.5_WP*L
  end do
  
  ! Create the mid points
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do

  ! Create the masks
  mask = 0
  
  return
end subroutine dns_box_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine dns_box_data
  use dns_box
  use random
  use parser
  use fileio
  implicit none
  
  ! Turbulent velocity
  real(WP) :: Ut
  
  ! Spectrum type
  character(len=str_short) :: spectrum
  real(WP) :: le,ld
  real(WP) :: epsilon
  
  ! Spectrum computation
  real(WP) :: psr,ps1,ps2
  real(WP) :: ke,kd,ks,dk,ksk0ratio,kc,kcksratio,kk,kx,ky,kz,kk2
  real(WP) :: alpha,spec_amp,eps,amp_disc
  integer  :: nk
  real(WP) :: e_total,energy_spec
  real(WP), dimension(:,:), pointer :: spect
  complex(WP), dimension(:,:,:), pointer :: ak,bk
  
  ! Other
  integer :: i,j,k,ik,iunit,dim
  complex(WP) :: ii=(0.0_WP,1.0_WP)
  real(WP) :: rand,pi
  
  ! Fourier coefficients
  integer(KIND=8) :: plan_r2c,plan_c2r
  complex(WP), dimension(:,:,:), pointer :: Uk,Vk,Wk
  complex(WP), dimension(:,:,:), pointer :: Cbuf
  real(WP), dimension(:,:,:), pointer :: Rbuf
  real(WP) :: f_phi
  
  ! Initialize the random number generator
  call random_init
  
  ! Decide how to generate the input file
  call parser_read('Data type',data_type)
  
  ! Allocate the array data
  select case (trim(data_type))
  case ('cold')
     nvar = 4
     allocate(data(nx,nx,nx,nvar))
     allocate(names(nvar))
     names = ''
     data = 0.0_WP
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
  case ('passive mixing')
     nvar = 5
     allocate(data(nx,nx,nx,nvar))
     allocate(names(nvar))
     names = ''
     data = 0.0_WP
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ZMIX => data(:,:,:,5);  names(5)  = 'ZMIX'
  case ('hot')
     nvar = 7
     allocate(data(nx,nx,nx,nvar))
     allocate(names(nvar))
     names = ''
     data = 0.0_WP
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ZMIX => data(:,:,:,5);  names(5)  = 'ZMIX'
     RHO  => data(:,:,:,6);  names(6)  = 'RHO'
     dRHO => data(:,:,:,7);  names(7)  = 'dRHO'
  case ('level set')
     nvar = 5
     allocate(data(nx,nx,nx,nvar))
     allocate(names(nvar))
     names = ''
     data = 0.0_WP
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ZMIX => data(:,:,:,5);  names(5)  = 'LVLSET'
  case default
     stop "dns_box_data: unknown Data type"
  end select
  
  ! Read spectrum parameters
  call parser_read('Fluctuations',Ut)
  call parser_read('Spectrum form',spectrum)
  call parser_read('Energetic scale',le)
  if (trim(spectrum).eq.'VKP') then
     call parser_read('Dissipative scale',ld)
     call parser_read('Dissipation',epsilon)
  end if
  
  ! Create pi
  pi = acos(-1.0_WP)
  
  ! =================
  ! Velocity Spectrum
  
  ! Spectrum computation
  ke = 2.0_WP*pi/le
  if (trim(spectrum).eq.'VKP') kd = 2.0_WP*pi/ld
  dk = 2.0_WP*pi/L
  kc = real(nx/2,WP)*dk
  eps=ke/1000000.0_WP
  if (trim(spectrum).eq.'PP') then
     spec_amp = 16.0_WP*sqrt(2.0_WP/pi)*Ut**2/ke
  else if (trim(spectrum).eq.'VKP') then
     alpha=1.5_WP
     spec_amp = 1.5_WP*Ut**5/epsilon
  else
     stop "dns_box_data: unknown Spectrum form"
  end if
  amp_disc = sqrt(dk)**3
  
  ! Compute spectrum
  nk = nx/2+1
  allocate(ak(nk,ny,nz),bk(nk,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nk
           ! Random numbers
           call random_number(rand)
           psr=2.0_WP*pi*(rand-0.5_WP)
           call random_number(rand)
           ps1=2.0_WP*pi*(rand-0.5_WP)
           call random_number(rand)
           ps2=2.0_WP*pi*(rand-0.5_WP)
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
           kk=sqrt(kx**2+ky**2+kz**2)
           ! Spectrums
           if (trim(spectrum).eq.'PP') then
              energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
           else if (trim(spectrum).eq.'VKP') then
              energy_spec=spec_amp*(kk/ke)**4/(1.0_WP+(kk/ke)**2)**(17.0_WP/6.0_WP)*exp(-1.5_WP*alpha*(kk/kd)**(4.0_WP/3.0_WP))
           end if
           ! Coeff
           ak(i,j,k)=0.0_WP
           bk(i,j,k)=0.0_WP
           if ((kk.gt.eps).and.(kk.le.kc)) then
              ak(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*exp(ii*ps1)*cos(psr)
              bk(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*exp(ii*ps2)*sin(psr)
           end if
        end do
     end do
  end do
  
  ! Output spectrum for comparison
  allocate(spect(2,nx+1))
  e_total=0.0_WP
  do ik=1,nx+1
     kk = dk*(ik-1)
     spect(1,ik)=kk
     if (trim(spectrum).eq.'PP') then
        energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
     else if (trim(spectrum).eq.'VKP') then
        energy_spec=spec_amp*(kk/ke)**4/(1.0_WP+(kk/ke)**2)**(17.0_WP/6.0_WP)*exp(-1.5_WP*alpha*(kk/kd)**(4.0_WP/3.0_WP))
     end if
     if ((kk.gt.eps).and.(kk.le.kc)) then
        spect(2,ik) = energy_spec
        e_total=e_total+dk*energy_spec
     else
        spect(2,ik) = 0.0_WP
     end if
  end do
  iunit=iopen()
  open(iunit,file='spectrum.analytic',form='formatted')
  do i=1,nx+1
     if ((spect(1,i).ne.0.0_WP).and.(spect(2,i).ne.0.0_WP)) then
        write(11,*) spect(1,i),spect(2,i)  
     end if
  end do
  close(iclose(iunit))
  
  ! Compute 3D velocity field
  allocate(Uk(nk,ny,nz))
  allocate(Vk(nk,ny,nz))
  allocate(Wk(nk,ny,nz))
  Uk=(0.0_WP,0.0_WP)
  Vk=(0.0_WP,0.0_WP)
  Wk=(0.0_WP,0.0_WP)
  
  do dim=1,3
     
     ! Compute the Fourier coefficients
     do k=1,nz
        do j=1,ny
           do i=1,nk
              
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
              kk =sqrt(kx**2+ky**2+kz**2)
              kk2=sqrt(kx**2+ky**2)
              
              ! Compute the Fourier coefficients
              if ((kk.gt.eps).and.(kk.le.kc)) then
                 if (dim.eq.1) then
                    if (kk2.lt.eps) then
                       Uk(i,j,k)=(ak(i,j,k)+bk(i,j,k))/sqrt(2.0_WP)
                    else
                       Uk(i,j,k)=(ak(i,j,k)*kk*ky+bk(i,j,k)*kx*kz)/(kk*kk2)
                    end if
                 end if
                 if (dim.eq.2) then
                    if (kk2.lt.eps) then
                       Vk(i,j,k)=(bk(i,j,k)-ak(i,j,k))/sqrt(2.0_WP)
                    else
                       Vk(i,j,k)=(bk(i,j,k)*ky*kz-ak(i,j,k)*kk*kx)/(kk*kk2)
                    end if
                 end if
                 if (dim.eq.3) Wk(i,j,k)=-bk(i,j,k)*kk2/kk
              end if
              
           end do
        end do
     end do
     
  end do
  
  ! Oddball
  do k=2,nz
     do j=nk+1,ny
        Uk(1,j,k)=conjg(Uk(1,ny+2-j,nz+2-k))
        Vk(1,j,k)=conjg(Vk(1,ny+2-j,nz+2-k))
        Wk(1,j,k)=conjg(Wk(1,ny+2-j,nz+2-k))
     end do
  end do
  do k=nk+1,nz
     Uk(1,1,k)=conjg(Uk(1,1,nz+2-k))
     Vk(1,1,k)=conjg(Vk(1,1,nz+2-k))
     Wk(1,1,k)=conjg(Wk(1,1,nz+2-k))
  end do
  
  ! Inverse Fourier transform
  allocate(Cbuf(nk,ny,nz))
  allocate(Rbuf(nx,ny,nz))
  call dfftw_plan_dft_c2r_3d(plan_c2r,nx,ny,nz,Cbuf,Rbuf,FFTW_ESTIMATE)
  call dfftw_plan_dft_r2c_3d(plan_r2c,nx,ny,nz,Rbuf,Cbuf,FFTW_ESTIMATE)
  
  ! Execute the plans
  do k=1,nz
     do j=1,ny
        do i=1,nk
           Cbuf(i,j,k) = Uk(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan_c2r)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           U(i,j,k) = Rbuf(i,j,k)
        end do
     end do
  end do
  do k=1,nz
     do j=1,ny
        do i=1,nk
           Cbuf(i,j,k) = Vk(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan_c2r)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           V(i,j,k) = Rbuf(i,j,k)
        end do
     end do
  end do
  do k=1,nz
     do j=1,ny
        do i=1,nk
           Cbuf(i,j,k) = Wk(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan_c2r)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           W(i,j,k) = Rbuf(i,j,k)
        end do
     end do
  end do
    
  ! Clean up
  deallocate(Uk)
  deallocate(Vk)
  deallocate(Wk)
  deallocate(ak)
  deallocate(bk)
  
  ! Check the results
  call dns_box_check(spectrum,Ut,ke)
  
  ! ===============
  ! Scalar Spectrum
  select case (trim(data_type))
  case ('cold')
     ! Nothing to do, no scalar
  case ('passive mixing','hot')
     
     ! Initialize in similar manner to Eswaran and Pope 1988
     call parser_read('ks/ko',ksk0ratio)
     ks = ksk0ratio*dk
     call parser_read('kc/ks',kcksratio)
     kc = kcksratio*ks
     
     ! Compute the Fourier coefficients
     do k=1,nz
        do j=1,ny
           do i=1,nk
              
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
              kk =sqrt(kx**2+ky**2+kz**2)
              kk2=sqrt(kx**2+ky**2)
              
              ! Compute the Fourier coefficients
              if ((ks-dk/2.0_WP.le.kk).and.(kk.le.ks+dk/2.0_WP)) then
                 f_phi = 1.0_WP
              else
                 f_phi = 0.0_WP
              end if
              call random_number(rand)
              if (kk.lt.1e-10) then
                 Cbuf(i,j,k) = 0.0_WP
              else
                 Cbuf(i,j,k) = sqrt(f_phi/(4.0_WP*pi*kk**2))*exp(ii*2.0_WP*pi*rand)
              end if
              
           end do
        end do
     end do
     
     ! Oddball
     do k=2,nz
        do j=nk+1,ny
           Cbuf(1,j,k)=conjg(Cbuf(1,ny+2-j,nz+2-k))
        end do
     end do
     do k=nk+1,nz
        Cbuf(1,1,k)=conjg(Cbuf(1,1,nz+2-k))
     end do
     
     ! Inverse Fourier transform
     call dfftw_execute(plan_c2r)
     
     ! Force 'double-delta' pdf on scalar field
     do k=1,nz
        do j=1,ny
           do i=1,nx
              if (Rbuf(i,j,k).le.0.0_WP) then
                 Rbuf(i,j,k) = 0.0_WP
              else
                 Rbuf(i,j,k) = 1.0_WP
              end if
           end do
        end do
     end do
     
     ! Fourier Transform and filter to smooth
     call dfftw_execute(plan_r2c)
     
     do k=1,nz
        do j=1,ny
           do i=1,nk
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
              kk =sqrt(kx**2+ky**2+kz**2)
              kk2=sqrt(kx**2+ky**2)
              
              ! Filter to remove high wavenumber components
              if (kk.le.kc) then
                 Cbuf(i,j,k) = Cbuf(i,j,k) * 1.0_WP
              else
                 Cbuf(i,j,k) = Cbuf(i,j,k) * (kc/kk)**2
              end if
              
           end do
        end do
     end do
     
     ! Oddball
     do k=2,nz
        do j=nk+1,ny
           Cbuf(1,j,k)=conjg(Cbuf(1,ny+2-j,nz+2-k))
        end do
     end do
     do k=nk+1,nz
        Cbuf(1,1,k)=conjg(Cbuf(1,1,nz+2-k))
     end do
     
     ! Fourier Transform back to real
     call dfftw_execute(plan_c2r)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = Rbuf(i,j,k)/real(nx*ny*nz,WP)
           end do
        end do
     end do
     
     ! Destroy the plans
     call dfftw_destroy_plan(plan_c2r)
     call dfftw_destroy_plan(plan_r2c)
     
     ! Clean up
     deallocate(Cbuf)
     deallocate(Rbuf)

  case ('level set')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = xm(i) - xm(floor(real(nx,WP)/4.0_WP))
           end do
        end do
     end do
  end select
  
  ! Density field if necessary
  if (trim(data_type).eq.'hot') then
     ! Density
     call parser_read('Density ratio',density_ratio)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = min(1.0_WP,max(0.0_WP,ZMIX(i,j,k)))
              RHO(i,j,k) = 1.0_WP/(1.0_WP+(density_ratio-1.0_WP)*ZMIX(i,j,k))
              U(i,j,k)  = U(i,j,k) / RHO(i,j,k)
              V(i,j,k)  = V(i,j,k) / RHO(i,j,k)
              W(i,j,k)  = W(i,j,k) / RHO(i,j,k)
           end do
        end do
     end do
     dRHO = 0.0_WP
  end if
  
  return
end subroutine dns_box_data

! ========================= !
! Create the chemtable data !
! ========================= !
subroutine dns_box_chemtable
  use dns_box
  use parser
  implicit none

  integer  :: i
  real(WP) :: visc,diff
  
  ! Return if chemtable not needed
  if (trim(data_type).ne.'hot') return
  
  ! Setup the coordinates
  n1 = 100
  n2 = 2
  n3 = 2

  ! Allocate arrays
  allocate(x1(n1))
  allocate(x2(n2))
  allocate(x3(n3))

  ! Setup the arrays
  do i=1,n1
     x1(i) = real(i-1,WP)/real(n1-1,WP)
  end do
  x2 = 0.0_WP
  x3 = 0.0_WP

  ! Chemtable model
  combModel = 'Steady Flamelet'

  ! Table itself
  nvar_chem = 4
  allocate(table(n1,n2,n3,nvar_chem))
  allocate(names_chem(nvar_chem))
  allocate(chem_mask(n1,n2,n3))

  ! Masked values in the chemtable
  chem_mask = 0

  ! Density
  names_chem(1)  = 'RHO'
  do i=1,n1
     table(i,:,:,1) = 1.0_WP/(1.0_WP+(density_ratio-1.0_WP)*x1(i))
  end do
  ! Temperature
  names_chem(2)  = 'T'
  table(:,:,:,2) = 1.0_WP/table(:,:,:,1)
  ! Viscosity
  call parser_read('Kinematic viscosity',visc)
  names_chem(3)  = 'VISC'
  table(:,:,:,3) = visc*table(:,:,:,1)
  ! Diffusivity
  call parser_read('Kinematic diffusivity',diff)
  names_chem(4)  = 'DIFF'
  table(:,:,:,4) = diff*table(:,:,:,1)

  return
end subroutine dns_box_chemtable


subroutine dns_box_check(spectrum,Ut,ke)
  use dns_box
  use fileio
  use parser
  implicit none
  integer :: iunit,i
  real(WP) :: pi
  real(WP) :: Ut,ke,mu_cin,rho_cin
  character(len=*) :: spectrum
  real(WP) :: kcin,epsi,l11,lt,re_turb,tau_epsi,l_k,tau_k
  real(WP), dimension(:,:), pointer :: spect
  
  pi = acos(-1.0_WP)
  call parser_read('Viscosity',mu_cin)
  call parser_read('Density',rho_cin)
  kcin     = 3.0_WP*0.5_WP*Ut*Ut
  if (trim(spectrum).eq.'PP') then
     epsi  = 15.0_WP*Ut*Ut*ke*ke*mu_cin/rho_cin/4.0_WP
     l11   = sqrt(2.0_WP*pi)/ke
  else if (trim(spectrum).eq.'VKP') then
     call parser_read('Dissipation',epsi)
  end if
  lt       = Ut**3.0_WP/epsi
  re_turb  = rho_cin*Ut**4.0_WP/mu_cin/epsi
  tau_epsi = kcin/epsi
  l_k      = ((mu_cin/rho_cin)**0.75_WP)/(epsi**0.25_WP)
  tau_k    = sqrt(mu_cin/rho_cin/epsi)
  
  write(*,*)  
  write(*,*)' ======================================== '
  write(*,*)' Debugging turbulent values               '
  write(*,*)' ---------------------------------------- '
  if (trim(spectrum).eq.'PP') then
     write(*,*)' -Spectrum type ----------> Passot-Pouquet'
  else if (trim(spectrum).eq.'VKP') then
     write(*,*)' -Spectrum type ----------> Von Karman-Pao'
  end if
  write(*,*)' -Turbulent Reynolds '
  write(*,*)'   Re_t --------> ', re_turb
  write(*,*)' -Turbulent kinetic energy '
  write(*,*)'   k -----------> ', kcin
  write(*,*)' -Turbulent dissipation rate '
  write(*,*)'   epsilon -----> ', epsi
  write(*,*)' -Size of the biggest eddy '
  write(*,*)'   l_t ---------> ', lt
  write(*,*)'   C1 ----------> ', L/lt
  write(*,*)' -Eddy turn over time '
  write(*,*)'   tau ---------> ', tau_epsi
  write(*,*)' -Kolmogorov scale '
  write(*,*)'   l_k ---------> ', l_k
  write(*,*)'   C2 ----------> ', l_k/dx
  write(*,*)' -Kolmogorov time scale '
  write(*,*)'   tau_k -------> ', tau_k
  write(*,*)' -Reynolds lambda '
  write(*,*)'   Re_lambda ---> ', sqrt(re_turb*15.0_WP)
  if (trim(spectrum).eq.'PP') then
     write(*,*)' -Integral length scale '
     write(*,*)'   Lii --------> ', l11
     write(*,*)' -Turbulent Reynolds based on Lii '
     write(*,*)'   Re_Lii -----> ', rho_cin*Ut*l11/mu_cin
  end if
  write(*,*)' ======================================== '
  write(*,*)
  write(*,*)
  
  ! Need to change that
  call dns_box_stat
  
  ! Recompute the spectrum
  allocate(spect(2,nx+1))
  spect = 0.0_WP
  call get_spectrum(nx,L,U,spect)
  call get_spectrum(nx,L,V,spect)
  call get_spectrum(nx,L,W,spect)
  
  ! Output it
  iunit=iopen()
  open(iunit,file='spectrum.numeric',form='formatted')
  do i=1,nx+1
     if ((spect(1,i).ne.0.0_WP).and.(spect(2,i).ne.0.0_WP)) then
        write(11,*) spect(1,i),spect(2,i)
     end if
  end do
  close(iclose(iunit))
  
  return
end subroutine dns_box_check



subroutine get_spectrum(nx,L,A,S)
  use precision
  implicit none
  include 'fftw3.f'  
  integer, intent(in) :: nx
  real(WP), intent(in) :: L
  real(WP), dimension(nx,nx,nx), intent(in) :: A
  real(WP), dimension(2,nx+1), intent(out) :: S
  complex(WP), dimension(:,:,:), pointer :: in,out
  real(WP), dimension(:,:,:), pointer :: B
  real(WP) :: pi,dk,kc,eps,kx,ky,kz,kk
  integer(KIND=8) :: plan
  integer :: i,j,k,ik
  
  ! Create the plan
  allocate(in(nx,nx,nx))
  allocate(out(nx,nx,nx))
  call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  
  ! FFT
  allocate(B(nx,nx,nx))
  do k=1,nx
     do j=1,nx
        do i=1,nx
           in(i,j,k) = A(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           out(i,j,k) = out(i,j,k)/(nx**3)
           B(i,j,k) = sqrt(real(out(i,j,k)*conjg(out(i,j,k))))
        end do
     end do
  end do
  
  ! Compute the spectrum
  pi=acos(-1.0_WP)
  dk=2.0_WP*pi/L
  kc=pi*nx/L
  eps=kc/1000000.0_WP
  do k=1,nx
     do j=1,nx
        do i=1,nx
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.(nx/2+1)) ky=-real(nx-j+1,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.(nx/2+1)) kz=-real(nx-k+1,WP)*dk
           kk=sqrt(kx**2+ky**2+kz**2)
           ! Spectrum
           ik=1+idint(kk/dk+0.5_WP)
           if ((kk.gt.eps).and.(kk.le.kc)) then
              S(2,ik)=S(2,ik)+0.5_WP*(B(i,j,k)**2)/dk
           end if
        end do
     end do
  end do
  do ik=1,nx+1
     S(1,ik)=dk*(ik-1)
  end do
  
  ! Clean up
  call dfftw_destroy_plan(plan)
  deallocate(B)
  deallocate(out)
  deallocate(in)
  
  return
end subroutine get_spectrum


subroutine dns_box_stat
  use dns_box
  use parser
  use fileio
  implicit none
  integer :: i,j,k
  real(WP) :: upvp,upwp,vpwp,Ut_final,kcin_final,div,div_max,div_b
  real(WP), dimension(3) :: um,umax,umin,uminval
  
  ! Computing mean velocity
  um=0.0_WP
  umax=-1.0e+70_WP
  umin=99.0e+70_WP
  uminval=99.0e+70_WP   
  do k=1,nz
     do j=1,ny
        do i=1,nx
           um(1)=um(1)+U(i,j,k)
           umax(1)=max(umax(1),(U(i,j,k)))
           umin(1)=min(umin(1),(U(i,j,k)))
           uminval(1)=min(uminval(1),abs(U(i,j,k)))
           um(2)=um(2)+V(i,j,k)
           umax(2)=max(umax(2),(V(i,j,k)))
           umin(2)=min(umin(2),(V(i,j,k)))
           uminval(2)=min(uminval(2),abs(V(i,j,k)))
           um(3)=um(3)+W(i,j,k)
           umax(3)=max(umax(3),(W(i,j,k)))
           umin(3)=min(umin(3),(W(i,j,k)))
           uminval(3)=min(uminval(3),abs(W(i,j,k)))
        end do
     end do
  end do
  um(1)=um(1)/(nx*ny*nz)
  um(2)=um(2)/(nx*ny*nz)
  um(3)=um(3)/(nx*ny*nz)
  
  ! Computing cross correlations to confirm isotropy
  upvp=0.0_WP
  upwp=0.0_WP
  vpwp=0.0_WP
  do k=1,nz 
     do j=1,ny
        do i=1,nx
           upvp=upvp+(U(i,j,k)-um(1))*(V(i,j,k)-um(2))
           upwp=upwp+(U(i,j,k)-um(1))*(W(i,j,k)-um(3))
           vpwp=vpwp+(V(i,j,k)-um(2))*(W(i,j,k)-um(3))
        end do
     end do
  end do
  upvp = upvp / (nx*ny*nz)
  upwp = upwp / (nx*ny*nz)
  vpwp = vpwp / (nx*ny*nz)
  
  ! Computing the final Ut
  Ut_final = 0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           Ut_final = Ut_final + (U(i,j,k)-um(1))**2+(V(i,j,k)-um(2))**2+(W(i,j,k)-um(3))**2
        end do
     end do
  end do
  Ut_final = sqrt(Ut_final/nx/ny/nz/3.0_WP)
  kcin_final = 3.0_WP*0.5_WP*Ut_final*Ut_final
  
  ! Computing divergence values
  div_b=0.0_WP
  div_max=0.0_WP
  do k=2,nz-1
     do j=2,ny-1
        do i=2,nx-1
           div=abs(U(i+1,j,k)-U(i-1,j,k)+V(i,j+1,k)-V(i,j-1,k)+W(i,j,k+1)-W(i,j,k-1))
           div_max=max(div_max,div)
           div_b=div_b+div
        end do
     end do
  end do
  div_b=div_b/(nx-2)/(ny-2)/(nz-2)
  
  ! Output it ------------------------------------------
  write(*,*)  
  write(*,*)' ======================================== '
  write(*,*)' Debugging turbulent values given by      '
  write(*,*)' the generated field                      '
  write(*,*)' ---------------------------------------- '
  write(*,*)' -Turbulent velocity '
  write(*,*)'   up ---------> ', Ut_final
  write(*,*)' -Turbulent kinetic energy '
  write(*,*)'   k ----------> ', kcin_final
  write(*,*)' -Cross correlations'
  write(*,*)'   <upvp> -----> ', upvp
  write(*,*)'   <upwp> -----> ', upwp
  write(*,*)'   <vpwp> -----> ', vpwp
  write(*,*)' -Velocity min, max and mean'
  write(*,*)'   u_min ------> ', umin(1)
  write(*,*)'   u_mean -----> ', um(1)
  write(*,*)'   u_max ------> ', umax(1)
  write(*,*)'   |u|_min ----> ', uminval(1)
  write(*,*)'   v_min ------> ', umin(2)
  write(*,*)'   v_mean -----> ', um(2)
  write(*,*)'   v_max ------> ', umax(2)
  write(*,*)'   |v|_min ----> ', uminval(2)
  write(*,*)'   w_min ------> ', umin(3)
  write(*,*)'   w_mean -----> ', um(3)
  write(*,*)'   w_max ------> ', umax(3)
  write(*,*)'   |w|_min ----> ', uminval(3)
  write(*,*)' -Average divergence -1 order'
  write(*,*)'   Div_b ------> ', div_b
  write(*,*)' -Maximum divergence'
  write(*,*)'   Div_max ----> ', div_max
  write(*,*)' ======================================== '
  write(*,*)
  write(*,*)
  
  return
end subroutine dns_box_stat
