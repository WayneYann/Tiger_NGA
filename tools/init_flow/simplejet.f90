module simplejet
  use precision
  use param
  implicit none
  
  ! Length and diameter of the combustor
  real(WP) :: length, diameter
  ! Diameter of the pipe
  real(WP) :: jet_diameter
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in radial direction
  integer :: inflow_points, jet_points, wall_points
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX

end module simplejet

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine simplejet_grid
  use simplejet
  use parser
  implicit none

  integer :: i,j,k
  real(WP), parameter :: twoPi = 2.0_WP*acos(-1.0_WP)

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Domain length',length)
  call parser_read('Domain diameter',diameter)

  call parser_read('Jet diameter',jet_diameter)

  call parser_read('Inflow length',inflow_length)

  call parser_read('Inflow points',inflow_points)
  call parser_read('Jet points',jet_points)
  call parser_read('Wall points',wall_points)

  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid

  ! Inflow
  do i=1,inflow_points+1
     x(i) = real(i-1,WP)*inflow_length/real(inflow_points,WP) - inflow_length
  end do
  ! Remainder of domain
  do i=inflow_points+1,nx+1
     x(i) = 0.0_WP
  end do

  ! Central jet
  do j=1,jet_points+1
     y(j) = 0.0_WP
  end do
  ! Wall
  do j=jet_points+1,jet_points+wall_points+1
     y(j) = 0.0_WP
  end do
  ! Coflow
  do j=jet_points+wall_points+1,ny+1
     y(j) = 0.0_WP
  end do

  ! Circumferential
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do

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
  mask(1:inflow_points,jet_points+1:jet_points+wall_points) = 1

  return
end subroutine simplejet_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine simplejet_data
  use simplejet
  use parser
  implicit none

  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U    => data(:,:,:,1); names(1) = 'U'
  V    => data(:,:,:,2); names(2) = 'V'
  W    => data(:,:,:,3); names(3) = 'W'
  P    => data(:,:,:,4); names(4) = 'P'
  ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  
  ! Create them
  U    = 0.0_WP
  V    = 0.0_WP
  W    = 0.0_WP
  P    = 0.0_WP
  ZMIX = 0.0_WP
  
  return
end subroutine simplejet_data

