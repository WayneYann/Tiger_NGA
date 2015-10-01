module wave
  use string
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module wave

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine wave_grid
  use wave
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  ny = nx
  nz = 1
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
     x(i) = real(i-1,WP)*twoPi/real(nx,WP)-Pi
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*twoPi/real(nx,WP)-Pi
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nx,WP)
  end do
  
  ! Create the cell centered grid
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
end subroutine wave_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine wave_data
  use wave
  use parser
  use math
  implicit none
  real(WP) :: amp
  integer :: i,j,k
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  
  ! Initialize
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  G = huge(1.0_WP)
  
  ! Setup the interface
  call parser_read('Initial amplitude',amp,0.01_WP*twoPi)
  print*,'amp=',amp
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = -ym(j)+amp*cos(xm(i))
        end do
     end do
  end do
  
  return
end subroutine wave_data
