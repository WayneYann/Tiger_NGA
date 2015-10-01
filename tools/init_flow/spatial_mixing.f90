module spatial_mixing
  use precision
  use param
  implicit none

  ! Length and diameter of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: out_a
  integer  :: out_ny,core_ny

  ! Inflow profile
  character(len=str_medium) :: inflow_type
  ! Convective velocity and velocity difference
  real(WP) :: U_convective,delta_U
  real(WP) :: U1,U2
  ! Half thickness of the layer
  real(WP) :: thick
  ! Reynolds number of the mean flow
  real(WP) :: Re
    
  ! Blasius solution
  real(WP) :: sqrt_nu_x,thick_gauss,yp
  integer, parameter :: no = 10
  real(WP), dimension(-no:+no) :: sten1,sten2,gauss,xg

  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U,Ur,Ui
  real(WP), dimension(:,:,:), pointer :: V,Vr,Vi
  real(WP), dimension(:,:,:), pointer :: W,Wr,Wi
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: rhoU,rhoV,rhoW
  real(WP), dimension(:,:,:), pointer :: RHO,dRHO
  
  ! Flow Field stats
  real(WP), dimension(:), pointer :: Ubase   ! [m/s]
  real(WP), dimension(:), pointer :: Up      ! [1/s]
  real(WP), dimension(:), pointer :: Upp     ! [1/m*s]
  
  ! Spatial Orr Sommerfeld
  integer :: nos
  real(WP), dimension(:), pointer :: xos,yos
  real(WP), parameter :: stretching = 3.0_WP
  real(WP), parameter ::  h = 1.0E-4_WP
  real(WP),    dimension(:,:), pointer :: Id,D1,D2,D3,D4,T,G,Tp
  complex(WP), dimension(:,:), pointer :: A,B,C
  complex(WP), dimension(:,:), pointer :: A11,A12,A13,A14,A21,A32,A43,B11,B22,B33,B44
  complex(WP), dimension(:),   pointer :: eta_tilde,u_tilde,v_tilde,w_tilde
  
  ! Orr-Sommerfeld - lapack solver
  integer :: ierr,lwork
  integer,     dimension(:),   pointer :: iwork
  complex(WP), dimension(:),   pointer :: work
  real(WP),    dimension(:),   pointer :: rwork
  complex(WP), dimension(:),   pointer :: eigval_a,eigval_b  ! a/b = lambda = eigenvalue
  complex(WP), dimension(:,:), pointer :: eigvec_l           ! left eigenvector
  complex(WP), dimension(:,:), pointer :: eigvec_r           ! right eigenvector
  
  ! Complex constants
  complex(WP), parameter :: zero = (0.0_WP,0.0_WP)
  complex(WP), parameter :: one  = (1.0_WP,0.0_WP)
  complex(WP), parameter :: ii   = (0.0_WP,1.0_WP)

end module spatial_mixing

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine spatial_mixing_grid
  use spatial_mixing
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',core_ny)
  call parser_read('nz',nz)

  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Lz',Lz)
  
  call parser_read('Outer ny',out_ny)
  call parser_read('Outer ratio',out_a)

  ny = core_ny + out_ny

  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = (i-1)*Lx/real(nx,WP)
  end do
  y(ny/2+1) = 0.0_WP
  do j=1,core_ny/2
     y(ny/2+1+j) = +real(j,WP)*Ly/real(core_ny,WP)
     y(ny/2+1-j) = -real(j,WP)*Ly/real(core_ny,WP)
  end do
  do j=1,out_ny/2
     y(ny/2+1+core_ny/2+j) = y(ny/2+1+core_ny/2+j-1) + out_a*(y(ny/2+1+core_ny/2+j-1)-y(ny/2+1+core_ny/2+j-2))
     y(ny/2+1-core_ny/2-j) = y(ny/2+1-core_ny/2-j+1) + out_a*(y(ny/2+1-core_ny/2-j+1)-y(ny/2+1-core_ny/2-j+2))
  end do
  do k=1,nz+1
     z(k) = (k-1)*Lz/real(nz,WP) - 0.5_WP*Lz
  end do

  print*,maxval(abs(y))

  ! Create the mid points 
  do i=1,nx
     xm(i)= 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)= 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k)= 0.5_WP*(z(k)+z(k+1))
  end do
  
  ! Create the masks
  mask = 0  

  return
end subroutine spatial_mixing_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine spatial_mixing_data
  use spatial_mixing
  implicit none

  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  ZMIX = 0.0_WP

  ! Preset the mixing layer for faster convergence
  !U(:,1:ny/2,:)    = U1
  !U(:,ny/2+1:ny,:) = U2
  !ZMIX(:,1:ny/2,:)    = 0.0_WP
  !ZMIX(:,ny/2+1:ny,:) = 1.0_WP

  return
end subroutine spatial_mixing_data

! ===================================================== !
! Create an inflow profile for the spatial mixing layer !
! ===================================================== !
subroutine spatial_mixing_inflow
  use spatial_mixing
  use parser
  use fileio
  use math
  implicit none

  integer :: i,j,k,n
  
  ! Velocity generation
  real(WP) :: omega(9),beta(9),amp(9)
  complex(WP) :: alpha
  complex(WP), dimension(:), pointer :: u_fluct,v_fluct,w_fluct
  
  ! Max growth loop
  real(WP), dimension(3) :: growth
  real(WP) :: err_rel,omega_opt
  integer :: iter,iunit

  ! Init random module
  call random_init
  
  ! Read some parameters: sizes
  call parser_read('ntime',ntime)
  nvar_inflow = 8

  ! Allocate some arrays
  allocate(t_inflow(ntime))
  allocate(inflow(ntime,ny,nz,nvar_inflow))
  allocate(names_inflow(nvar_inflow))
  
  ! Link the pointers
  U  => inflow(:,:,:,1); names_inflow(1) = 'U'
  Ur => inflow(:,:,:,2); names_inflow(2) = 'Ur'
  Ui => inflow(:,:,:,3); names_inflow(3) = 'Ui'
  Vr => inflow(:,:,:,4); names_inflow(4) = 'Vr'
  Vi => inflow(:,:,:,5); names_inflow(5) = 'Vi'
  Wr => inflow(:,:,:,6); names_inflow(6) = 'Wr'
  Wi => inflow(:,:,:,7); names_inflow(7) = 'Wi'
  ZMIX => inflow(:,:,:,8); names_inflow(8) = 'ZMIX'

  ! Initialize the Orr-Sommerfeld solver
  call orr_sommerfeld_init

  ! Allocate perturbation arrays
  allocate(u_fluct(ny))
  allocate(v_fluct(ny))
  allocate(w_fluct(ny))

  ! *** Get most unstable mode ***

  ! Initial guess for omega (Valid for high Reynolds number)
  select case (trim(inflow_type))
  case ('tanh')
     omega_opt = 0.40_WP*U_convective/thick
  case ('blasius')
     omega_opt = 0.88_WP*U_convective/thick
  end select

  growth = 0.0_WP
  err_rel = huge(1.0_WP)
  iter = 1

  iunit = iopen()
  open(iunit,file='newton.txt')
  write(*,*)
  write(*,'(a)') "# Getting the most unstable frequency"
  write(iunit,'(a)') "# Getting the most unstable frequency"

  do while ((err_rel.gt.h).and.(iter.lt.10))
     call orr_sommerfeld_get_mode(omega_opt-h,0.0_WP,alpha,u_fluct,v_fluct,w_fluct)
     growth(1) = -aimag(alpha)
     write(iunit,*) omega_opt-h,alpha
     call orr_sommerfeld_get_mode(omega_opt,0.0_WP,alpha,u_fluct,v_fluct,w_fluct)
     growth(2) = -aimag(alpha)
     write(iunit,*) omega_opt,alpha
     call orr_sommerfeld_get_mode(omega_opt+h,0.0_WP,alpha,u_fluct,v_fluct,w_fluct)
     growth(3) = -aimag(alpha)
     write(iunit,*) omega_opt+h,alpha
     err_rel = 0.5_WP*abs(h*(growth(3)-growth(1))/(growth(1)+growth(3)-2.0_WP*growth(2)))
     omega_opt = omega_opt - 0.5_WP*h*(growth(3)-growth(1))/(growth(1)+growth(3)-2.0_WP*growth(2))
     write(*,'(a12,e16.6,a16,e12.4)') "omega_opt =", omega_opt, "rel error =", err_rel
     write(iunit,'(a12,e16.6,a16,e12.4)') "omega_opt =", omega_opt, "rel error =", err_rel
  end do

  close(iclose(iunit))
  write(*,*)

  ! Compute time grid
  time_inflow = 4.0_WP*twoPi/omega_opt
  dt_inflow   = time_inflow/real(ntime,WP)  
  do i=1,ntime
     t_inflow(i) = dt_inflow * real(i-1,WP)
  end do

  ! *** Inflow initialization ***
  
  ! Mean flow
  Ur = 0.0_WP
  Vr = 0.0_WP
  Wr = 0.0_WP
  Ui = 0.0_WP
  Vi = 0.0_WP
  Wi = 0.0_WP
  select case (trim(inflow_type))
  case ('tanh')
     do j=1,ny
        U(:,j,:) = U_convective + 0.5_WP*delta_U * tanh(ym(j)/thick)
        if (ym(j).ge.0.0_WP) ZMIX(:,j,:) = 1.0_WP
     end do
  case ('blasius')
     do j=1,ny
        do i=-no,no
           yp = ym(j) + 0.25_WP*xg(i)*thick
           if (yp.ge.0.0_WP) then
              sten1(i) = U2 * blasius1(yp*sqrt(U2)/sqrt_nu_x)
              sten2(i) = 1.0_WP
           else
              sten1(i) = U1 * blasius1(-yp*sqrt(U1)/sqrt_nu_x)
              sten2(i) = 0.0_WP
           end if
        end do
        U   (:,j,:) = sum(gauss*sten1)
        ZMIX(:,j,:) = sum(gauss*sten2)
     end do
  end select

  ! 3D - Fluctuations
  beta = 0.0_WP

  omega(1) = omega_opt
  omega(2) = omega_opt/2.0_WP
  omega(3) = omega_opt/4.0_WP
  omega(4) = omega_opt
  omega(5) = omega_opt/2.0_WP
  omega(6) = omega_opt/4.0_WP
  omega(7) = omega_opt
  omega(8) = omega_opt/2.0_WP
  omega(9) = omega_opt/4.0_WP
  
  call parser_read('Fluct. amplitude',amp)
  write(*,'(a)') "# Setting the modes"
  write(*,'(a6,4a16)') "mode","beta","omega","alpha_r","alpha_i"

  ! Get purturbations
  do n=1,9
     call orr_sommerfeld_get_mode(omega(n),beta(n),alpha,u_fluct,v_fluct,w_fluct)
     write(*,'(i6,4e16.6)') n,beta(n),omega(n),alpha

     if (n.le.3) then
        beta(n+3) = +alpha
        beta(n+6) = -alpha
     end if
     
     do k=1,nz
        do j=1,ny
           do i=1,ntime
              Ur(i,j,k) = Ur(i,j,k) + amp(n)*real (u_fluct(j)*exp(ii*(beta(n)*zm(k)-omega(n)*t_inflow(i))))
              Vr(i,j,k) = Vr(i,j,k) + amp(n)*real (v_fluct(j)*exp(ii*(beta(n)*zm(k)-omega(n)*t_inflow(i))))
              Wr(i,j,k) = Wr(i,j,k) + amp(n)*real (w_fluct(j)*exp(ii*(beta(n)*z(k) -omega(n)*t_inflow(i))))
              Ui(i,j,k) = Ui(i,j,k) + amp(n)*aimag(u_fluct(j)*exp(ii*(beta(n)*zm(k)-omega(n)*t_inflow(i))))
              Vi(i,j,k) = Vi(i,j,k) + amp(n)*aimag(v_fluct(j)*exp(ii*(beta(n)*zm(k)-omega(n)*t_inflow(i))))
              Wi(i,j,k) = Wi(i,j,k) + amp(n)*aimag(w_fluct(j)*exp(ii*(beta(n)*z(k) -omega(n)*t_inflow(i))))
           end do
        end do
     end do
  end do
  
  return
end subroutine spatial_mixing_inflow


! ====================================== !
! Initialize the Orr_Sommerfeld solver   !
! Allocate and pre-compute some matrices !
! ====================================== !
subroutine orr_sommerfeld_init
  use spatial_mixing
  use parser
  use fileio
  use math
  implicit none

  integer :: i,j,iunit

  ! Read number of grid points
  call parser_read('Orr-Sommerfeld points',nos)

  ! Read some parameters: mean flow
  call parser_read('Inflow type', inflow_type)
  call parser_read('Convective velocity',U_convective)
  call parser_read('Velocity difference',delta_U)
  call parser_read('Half layer thickness',thick)
  U1 = U_convective - 0.5_WP*delta_U
  U2 = U_convective + 0.5_WP*delta_U
  
  ! Get the operating Reynolds number
  call parser_read('Reynolds number',Re)

  ! Create grid
  allocate(xos(nos))
  allocate(yos(nos))
  xos(1)   = +1.0_WP
  yos(1)   = +huge(1.0_WP)
  xos(nos) = -1.0_WP
  yos(nos) = -huge(1.0_WP)
  do j=2,nos-1
     xos(j) = cos(Pi*real(j-1,WP)/real(nos-1,WP))
     yos(j) = stretching * 0.5_WP * thick * log((1.0_WP+xos(j))/(1.0_WP-xos(j)))
  end do
  
  ! Orr-Sommerfeld - lapack solver
  allocate(work(8*nos))
  allocate(rwork(32*nos))
  allocate(eigval_a(4*nos))
  allocate(eigval_b(4*nos))
  allocate(eigvec_l(1,4*nos))
  allocate(eigvec_r(4*nos,4*nos))
  allocate(iwork(nos))

  ! Temporary solution from Orr-Sommerfeld
  allocate(u_tilde(nos))
  allocate(v_tilde(nos))
  allocate(w_tilde(nos))
  allocate(eta_tilde(nos))

  ! Setup eigenvalue problem
  allocate(A(4*nos,4*nos))
  allocate(B(4*nos,4*nos))
  allocate(C(1*nos,1*nos))
  
  ! Define the matrix by its blocs
  A11 => A(0*nos+1 :1*nos ,0*nos+1 :1*nos )
  A12 => A(0*nos+1 :1*nos ,1*nos+1 :2*nos )
  A13 => A(0*nos+1 :1*nos ,2*nos+1 :3*nos )
  A14 => A(0*nos+1 :1*nos ,3*nos+1 :4*nos )

  A21 => A(1*nos+1 :2*nos ,0*nos+1 :1*nos )
  A32 => A(2*nos+1 :3*nos ,1*nos+1 :2*nos )
  A43 => A(3*nos+1 :4*nos ,2*nos+1 :3*nos )

  B11 => B(0*nos+1 :1*nos ,0*nos+1 :1*nos )
  B22 => B(1*nos+1 :2*nos ,1*nos+1 :2*nos )
  B33 => B(2*nos+1 :3*nos ,2*nos+1 :3*nos )
  B44 => B(3*nos+1 :4*nos ,3*nos+1 :4*nos )
  
  ! Precompute the stencil for derivations
  allocate(Id(nos,nos))
  allocate(D1(nos,nos))
  allocate(D2(nos,nos))
  allocate(D3(nos,nos))
  allocate(D4(nos,nos))
  
  ! Chebyshev matrices
  allocate(T (nos,nos))
  allocate(G (nos,nos))
  allocate(Tp(nos,nos))
  
  Id = 0.0_WP
  D1 = 0.0_WP
  D2 = 0.0_WP
  D3 = 0.0_WP
  D4 = 0.0_WP

  do i=1,nos

     Id(i,i) = 1.0_WP

     do j=1,nos
        ! Chebyshev to real
        T(i,j) = cos(Pi*real((i-1)*(j-1))/real(nos-1))

        ! Real to Chebyshev
        Tp(i,j) = 2.0_WP * cos(Pi*real((i-1)*(j-1))/real(nos-1)) / real(nos-1)
        if (i.eq.1 .or. i.eq.nos) Tp(i,j) = 0.5_WP * Tp(i,j)
        if (j.eq.1 .or. j.eq.nos) Tp(i,j) = 0.5_WP * Tp(i,j)
        
        ! Derivation in Chebyshev space
        if (i.ge.j .or. mod(i+j,2).eq.0) then
           G(i,j) = 0.0_WP
        else
           if (i.eq.1 .or. i.eq.nos) then
              G(i,j) = real(j-1,WP)
           else
              G(i,j) = 2.0_WP*real(j-1,WP)
           end if
        end if
        
     end do
  end do
  
  ! Derivation in real space
  call dgemm('N','N',nos,nos,nos,1.0_WP,T,nos,G,nos,0.0_WP,D1,nos)
  G = D1
  call dgemm('N','N',nos,nos,nos,1.0_WP,G,nos,Tp,nos,0.0_WP,D1,nos)
  
  ! Account for the mapping
  do i=1,nos
     D1(i,:) = D1(i,:) * (1.0_WP-xos(i)**2) / (thick*stretching)
  end do
     
  ! Higher order derivatives
  call dgemm('N','N',nos,nos,nos,1.0_WP,D1,nos,D1,nos,0.0_WP,D2,nos)
  call dgemm('N','N',nos,nos,nos,1.0_WP,D1,nos,D2,nos,0.0_WP,D3,nos)
  call dgemm('N','N',nos,nos,nos,1.0_WP,D1,nos,D3,nos,0.0_WP,D4,nos)
  
  ! Create base flow
  allocate(Ubase(nos))
  allocate(Up   (nos))
  allocate(Upp  (nos))

  select case (trim(inflow_type))
  case ('tanh')
     do j=1,nos
        Ubase(j) = U_convective + 0.5_WP*delta_U * tanh(yos(j)/thick)
        Up(j)    = 0.5_WP*delta_U/thick / cosh(yos(j)/thick)**2
        Upp(j)   = -delta_U/thick**2 * tanh(yos(j)/thick) / cosh(yos(j)/thick)**2
     end do

  case ('blasius')
     sqrt_nu_x = thick / (1.5_WP*(U1**(-0.5_WP)+U2**(-0.5_WP)))
     do i=-no,no
        xg(i) = 3.0_WP * real(i,WP) / real(no,WP)
     end do
     gauss = exp(-0.5_WP*xg**2)
     gauss = gauss / sum(gauss)
     do j=1,nos
        do i=-no,no
           yp = yos(j) + 0.25_WP*xg(i)*thick
           if (yp.ge.0.0_WP) then
              sten1(i) = U2 * blasius1(yp*sqrt(U2)/sqrt_nu_x)
           else
              sten1(i) = U1 * blasius1(-yp*sqrt(U1)/sqrt_nu_x)
           end if
        end do
        Ubase(j) = sum(gauss*sten1)
     end do
     do j=1,nos
        Up(j) = sum(D1(j,:)*Ubase)
     end do
     do j=1,nos
        Upp(j) = sum(D1(j,:)*Up)
     end do
  end select

  ! Save to a file the mean flow
  iunit = iopen()
  open(iunit,file='mean.txt')
  do j=1,nos
     write(iunit,'(i3,4e16.6E3)') j,yos(j),Ubase(j),Up(j),Upp(j)
  end do
  close(iclose(iunit))

  return
end subroutine orr_sommerfeld_init


! ======================================================= !
! Compute the solution of the Orr-Sommerfeld equation for !
! -> Re    : Reynolds number                              !
! -> omega : time frequency                               !
! -> beta  : wave number in z                             !
! ======================================================= !
subroutine orr_sommerfeld_get_mode(omega,beta,alpha,u_fluct,v_fluct,w_fluct)
  use spatial_mixing
  use fileio
  implicit none
  
  real(WP), intent(in) :: omega,beta
  complex(WP), intent(out) :: alpha
  complex(WP), dimension(ny) :: u_fluct,v_fluct,w_fluct

  complex(WP), dimension(:), pointer :: coeff_u,coeff_v,coeff_w
  real(WP) :: speed,Umin,Umax,norm,xx,sumBC
  integer  :: j,jj,mode_index,iunit
  complex(WP) :: alpha_

  A = zero
  B = zero

  ! Interior points
  do j=2,nos

     A11(j,:) = -ii*Ubase(j)*Id(j,:)
     A12(j,:) = ii*omega*Id(j,:) + 2.0_WP*(D2(j,:)-beta**2*Id(j,:))/Re
     A13(j,:) = -ii*Upp(j)*Id(j,:) + ii*Ubase(j)*(D2(j,:)-beta**2*Id(j,:))
     A14(j,:) = -(D4(j,:)-2.0_WP*beta**2*D2(j,:)+beta**4*Id(j,:))/Re - ii*omega*(D2(j,:)-beta**2*Id(j,:))
     
     B11(j,:) = Id(j,:)/Re
  end do
  
  A21 = Id
  A32 = Id
  A43 = Id
  
  B22 = Id
  B33 = Id
  B44 = Id

  ! Boundary Conditions
  A14(1,1)     = one
  A14(nos,nos) = one
  
  ! Get eigenvalues and eigenvectors
  lwork = 8*nos
  call ZGGEV('N','V',4*nos,A,4*nos,B,4*nos,eigval_a,eigval_b,eigvec_l,1,eigvec_r,4*nos,work,lwork,rwork,ierr)
  if (ierr.ne.0) stop "ZGGEV: unable to compute the eigenvalues/vectors"
  
  ! Get the most unstable mode
  ! -> phase speed (real part) between min and max of Ubase
  ! -> negative imaginary part
  ! -> BCs at zero
  Umin = minval(Ubase)
  Umax = maxval(Ubase)

  mode_index = -1
  alpha = (0.0_WP,1.0_WP)*huge(1.0_WP)

  do j=1,4*nos
     alpha_ = eigval_a(j)/(eigval_b(j)+epsilon(abs(eigval_b(j))))
     speed = omega/real(alpha_)
     sumBC = abs(eigvec_r(3*nos+1,j)) + abs(eigvec_r(4*nos,j))
     
!!$     if (speed.ge.Umin .and. speed.le.Umax .and. aimag(alpha_).lt.0.0_WP) & 
!!$          write(*,'(A,i3,1f12.4,3f12.4,2e12.4)') '####  ',j,omega,alpha_,speed,&
!!$          abs(eigvec_r(3*nos+1,j)),abs(eigvec_r(4*nos,j))
     
     
     if (speed.ge.Umin .and. speed.le.Umax .and. sumBC.lt.1E-10_WP .and. &
          aimag(alpha_).lt.0.0_WP .and. aimag(alpha_).gt.-10.0_WP) then
        
        if (aimag(alpha_).lt.aimag(alpha)) then
           mode_index = j
           alpha = alpha_
        end if
     end if
  end do
  if (mode_index.eq.-1) stop "orr_sommerfeld_get_mode: Could not find a correct mode"

  ! Get v tilde
  v_tilde = eigvec_r(3*nos+1:4*nos,mode_index)

  ! Create the linear system for eta
  C = 0.0_WP
  do j=2,nos-1
     C(j,:) = (alpha*Ubase(j)-omega)*Id(j,:) + ii*(D2(j,:)-(alpha**2+beta**2)*Id(j,:))/Re
     eta_tilde(j) = -ii*beta*Up(j)*v_tilde(j)
  end do
  
  ! Boundary conditons
  C(1,1)     = one
  C(nos,nos) = one
  eta_tilde(1)   = zero
  eta_tilde(nos) = zero
  
  ! Solve it (linear solver)
  call ZGESV(nos,1,C,nos,iwork,eta_tilde,nos,ierr)
  
  ! Get back the two remaining components
  do j=1,nos
     u_tilde(j) = ii*sum(alpha*D1(j,:)*v_tilde-beta*eta_tilde) / (alpha**2+beta**2)
  end do
  w_tilde = (ii*eta_tilde+beta*u_tilde)/alpha
  
  ! Compute the Chebyshev coefficients
  coeff_u => work(1:nos)
  coeff_v => work(nos+1:2*nos)
  coeff_w => work(2*nos+1:3*nos)
  do j=1,nos
     coeff_u(j) = sum(Tp(j,:)*u_tilde)
     coeff_v(j) = sum(Tp(j,:)*v_tilde)
     coeff_w(j) = sum(Tp(j,:)*w_tilde)
  end do
  
  ! Compute the exact staggered fluctuations
  do jj=1,ny
     xx = acos(tanh(ym(jj)/(thick*stretching)))
     u_fluct(jj) = zero
     w_fluct(jj) = zero
     do j=1,nos
        u_fluct(jj) = u_fluct(jj) + coeff_u(j)*cos(real(j-1,WP)*xx)
        w_fluct(jj) = w_fluct(jj) + coeff_w(j)*cos(real(j-1,WP)*xx)
     end do
     xx = acos(tanh(y(jj)/(thick*stretching)))
     v_fluct(jj) = zero
     do j=1,nos
        v_fluct(jj) = v_fluct(jj) + coeff_v(j)*cos(real(j-1,WP)*xx)
     end do
  end do

  ! Normalize everything
  norm = maxval(abs(u_fluct)**2+abs(v_fluct)**2+abs(w_fluct)**2)
  u_fluct = u_fluct / norm
  v_fluct = v_fluct / norm
  w_fluct = w_fluct / norm

  ! Write the solution to a file
  iunit = iopen()
  open(iunit,file='mode.txt',form="formatted",iostat=ierr)
  do jj=1,ny
     write(iunit,'(i3,8e14.6)') jj,y(jj),ym(jj),u_fluct(jj),v_fluct(jj),w_fluct(jj)
  end do
  close(iclose(iunit))     

  return
end subroutine orr_sommerfeld_get_mode
