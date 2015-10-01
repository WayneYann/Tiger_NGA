program spectrum
  use precision
  use string
  use fileio
  implicit none
  
  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz
  integer :: xper,yper,zper
  integer :: icyl
  real(WP), dimension(:), pointer :: x,y,z
  integer, dimension(:,:), pointer :: mask
  !integer, dimension(:,:,:), pointer :: iblank
  integer :: iunit,ierr
  character(len=str_medium) :: filename1,filename2,filename3
  character(len=str_medium) :: config
  integer :: i
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data
  integer :: var,nvar
  real(WP) :: L
  real(WP) :: dt,time
  real(WP), dimension(:,:), pointer :: spect
  real(WP) :: epsilon
  
  ! Read file name from standard input
  print*,'=============================='
  print*,'| ARTS - Spectrum calculator |'
  print*,'=============================='
  print*
  print "(a,$)", " config file : "
  read "(a)", filename1
  print "(a,$)", " data file : "
  read "(a)", filename2
  print "(a,$)", " spectrum file : "
  read "(a)", filename3
  
  ! ** Open the config file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename1),"r",ierr)
  
  ! Read sizes
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  print*,'Config : ',config
  print*,'Grid :',nx,'x',ny,'x',nz
  
  ! Read grid field
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Get domain size
  L = x(nx+1) - x(1)
  print*,'Domain size:',L
  
  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"r",ierr)
  
  ! Read sizes
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
  end do
  
  ! Read data field
  allocate(data(nx,ny,nz,nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,data(:,:,:,var),nx*ny*nz,kind(data),ierr)
  end do
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Compute the spectrum
  allocate(spect(2,nx+1))
  spect   = 0.0_WP
  epsilon = 0.0_WP
  call get_spectrum(nx,L,data(:,:,:,1),spect,epsilon)
  call get_spectrum(nx,L,data(:,:,:,2),spect,epsilon)
  call get_spectrum(nx,L,data(:,:,:,3),spect,epsilon)
  
  ! Output it
  iunit=iopen()
  open(iunit,file=filename3,form='formatted')
  do i=1,nx+1
     if ((spect(1,i).ne.0.0_WP).and.(spect(2,i).ne.0.0_WP)) then
        write(11,*) spect(1,i),',',spect(2,i)
     end if
  end do
  close(iclose(iunit))
  
  ! Get and print some mean quantities
  print*,""
  print*,"TKE : ",sum(spect(2,:))
  print*,"epsilon (integral) : ",2.0_WP*sum(spect(1,:)**2*spect(2,:))
  print*,"epsilon (averaged) : ",epsilon/real(nx*nx*nx,WP)
  
  ! Deallocation
  deallocate(spect)
  deallocate(x,y,z)
  deallocate(mask)
  deallocate(data)
  
end program spectrum


! ============================================================ !
! Compute the spectrum of a velocity component and add it to S !
! ============================================================ !
subroutine get_spectrum(nx,L,A,S,epsilon)
  use precision
  implicit none
  include 'fftw3.f'  
  integer, intent(in) :: nx
  real(WP), intent(inout) :: L,epsilon
  real(WP), dimension(nx,nx,nx), intent(in) :: A
  real(WP), dimension(2,nx+1), intent(out) :: S
  complex(WP), dimension(:,:,:), pointer :: in,out,vel
  real(WP), dimension(:,:,:), pointer :: B
  real(WP) :: pi,dk,kc,eps,kx,ky,kz,kk
  integer(KIND=8) :: plan
  integer :: i,j,k,ik
  complex(WP), parameter :: ii =(0.0_WP,1.0_WP)
  
  ! Create the plan
  allocate(in (nx,nx,nx))
  allocate(out(nx,nx,nx))
  allocate(vel(nx,nx,nx))
  call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  
  ! FFT
  in = A
  call dfftw_execute(plan)
  out = out/(nx**3)
  allocate(B(nx,nx,nx))
  do k=1,nx
     do j=1,nx
        do i=1,nx
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
  
  ! Create the plan
  call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
  
  ! Save velocity
  do k=1,nx
     do j=1,nx
        do i=1,nx
           vel(i,j,k) = out(i,j,k)
        end do
     end do
  end do
  
  ! Compute the derivatives - X
  do k=1,nx
     do j=1,nx
        do i=1,nx
           kx=real(i-1,WP)*dk
           if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
           in(i,j,k) = ii*kx*vel(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           epsilon = epsilon + real(out(i,j,k),WP)**2
        end do
     end do
  end do
  
  ! Compute the derivatives - Y
  do k=1,nx
     do j=1,nx
        do i=1,nx
           ky=real(j-1,WP)*dk
           if (j.gt.(nx/2+1)) ky=-real(nx-j+1,WP)*dk
           in(i,j,k) = ii*ky*vel(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           epsilon = epsilon + real(out(i,j,k),WP)**2
        end do
     end do
  end do
  
  ! Compute the derivatives - Z
  do k=1,nx
     do j=1,nx
        do i=1,nx
           kz=real(k-1,WP)*dk
           if (k.gt.(nx/2+1)) kz=-real(nx-k+1,WP)*dk
           in(i,j,k) = ii*kz*vel(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           epsilon = epsilon + real(out(i,j,k),WP)**2
        end do
     end do
  end do
  
  ! Clean up
  call dfftw_destroy_plan(plan)
  deallocate(B)
  deallocate(vel)
  deallocate(out)
  deallocate(in)
  
  return
end subroutine get_spectrum

