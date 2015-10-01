module boundary_layer
  use precision
  use param
  implicit none

  ! Length of the domain
  real(WP) :: Lx,Ly,Lz

  ! Stretching
  real(WP) :: r
  
  ! U velocity at infinity
  real(WP) :: Uinft

  ! Momentum thickness
  real(WP) :: theta,delta
  
  ! Parameters for 'Law of the wall'
  real(WP), parameter :: kappa = 0.41_WP
  real(WP), parameter :: BigPI = 0.5_WP
  real(WP), parameter :: SiPi = 1.852_WP

  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P

end module boundary_layer

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine boundary_layer_grid
  use boundary_layer
  use parser
  implicit none

  integer :: i,j,k
  real(WP) :: ytilde

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Lz',Lz)

  call parser_read('Stretching',r)

  ! Set the periodicity
  xper = 1
  yper = 0
  zper = 1

  ! Cartesian
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1)*Lx/real(nx)
  end do
  do j=2,ny+1
     ytilde = real(ny+1-j)/real(ny-1)
     y(j) = Ly * (1.0_WP-tanh(r*ytilde)/tanh(r))
  end do
  y(1) = 2.0_WP*y(2)-y(3)
  do k=1,nz+1
     z(k) = real(k-1)/real(nz)*Lz - 0.5_WP*Lz
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
  mask(:,1)  = 1
  
  return
end subroutine boundary_layer_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine boundary_layer_data
  use boundary_layer
  use parser
  use random
  use math
  implicit none

  integer  :: i,j,k
  real(WP) :: amp,rnd,kl,Cf,eta

  ! Initialize the random number generator
  call random_init
  
  ! Read in the amplitude of random numbers
  call parser_read('Fluctuation magnitude',amp)
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Read the mean velocities
  call parser_read('Velocity at infinity',Uinft)

  ! Momentum thickness
  call parser_read('Momentum thickness',theta)
  call parser_read('Skin friction',Cf)
  kl = kappa*sqrt(2.0_WP/Cf)
  delta = theta / ( &
       + (1.0_WP+BigPi)/kl & 
       - (2.0_WP+2.0_WP*BigPi*(SiPi/Pi+1.0_WP)+1.5_WP*BigPi**2)/kl**2)
  print*,delta

  ! Create the velocities
  do j=2,ny
     eta = ym(j)/delta
     if (eta.le.1.0_WP) then
        U(:,j,:) = Uinft * (1.0_WP + &
             (log(eta)+2.0_WP*BigPi*((sin(0.5_WP*Pi*eta))**2-1.0_WP))/kl)
        do i=1,nx
           do k=1,nz
              call random_number(rnd)
              rnd = amp * (rnd-0.5_WP)
              U(i,j,k) = U(i,j,k) + rnd
           end do
        end do
     else
        U(:,j,:) = Uinft
     end if
  end do
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  U(:,1,:)  = 0.0_WP
  
  return
end subroutine boundary_layer_data

