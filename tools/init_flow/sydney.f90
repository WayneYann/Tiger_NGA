module sydney
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P

end module sydney

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine sydney_grid
  use sydney
  use parser
  use math
  implicit none

  integer :: i,j,k
  real :: jet_diameter, pilot_diameter
  real :: inner_wall, outer_wall

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  ! Geometrical parameters
  jet_diameter   = 0.0075_WP
  pilot_diameter = 0.018_WP
  inner_wall     = 0.00025_WP
  outer_wall     = 0.00020_WP

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

  ! Create the x-grid
  if (nx.eq.256) then
     ! Inflow (1D)
     do i=1,9
        x(i) = jet_diameter/8.0_WP*real(i-1,WP)-jet_diameter
     end do
     ! Flame
     do i=9,nx+1
        x(i) = 0.05678388_WP*1.00818261_WP**i-0.06110516_WP
     end do
  else
     write(*,'(a)') 'nx must be 256'
     call kill
  end if

  do i=1,nx+1
     print*, i,x(i)
  end do

  ! Create the y-grid
  if (ny.eq.128) then
     ! Fuel Jet
     y(1) = 0.0_WP
     do j=2,17
        y(j) = -0.00479441_WP*0.80993658_WP**j+0.00388317_WP
        !y(j) = -0.00255106_WP*0.85800087_WP**j+0.00218883_WP
     end do
     ! Inner Wall
     do j=17,25
        y(j) = inner_wall/8.0_WP*real(j-17,WP) + jet_diameter/2.0_WP
     end do
     ! Pilot
     do j=25,57
        y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
             tanh(2.0_WP*(2.0_WP*real(j-25,WP)/real(32,WP)-1.0_WP))/tanh(2.0_WP) + &
             inner_wall+jet_diameter/2.0_WP + &
             0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
     end do
     ! Outer Wall
     do j=57,63
        y(j) = outer_wall/6.0_WP*real(j-57,WP) + pilot_diameter/2.0_WP
     end do
     ! Coflow
     do j=63,129
        y(j) = 0.000000588244_WP*1.10501068_WP**j+0.00888257_WP
     end do
  else
     write(*,'(a)') 'ny must be 128'
     call kill
  end if

  do j=1,ny+1
     print*, j, y(j)
  end do

  ! Create the z-grid
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do

  ! Create the mask
  mask = 0
  mask(1:8,17:24) = 1
  mask(1:8,57:62) = 1
  
  return
end subroutine sydney_grid


! =============================== !
! Create the variable array: Data !
! =============================== !
subroutine sydney_data
  use sydney
  implicit none
  
  ! Allocate the array data
  nvar = 11
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Initialize the variables and names: Gas-phase
  data(:,:,:, 1) = 0.0_WP;      names( 1) = 'U'
  data(:,:,:, 2) = 0.0_WP;      names( 2) = 'V'
  data(:,:,:, 3) = 0.0_WP;      names( 3) = 'W'
  data(:,:,:, 4) = 0.0_WP;      names( 4) = 'P'
  data(:,:,:, 5) = 1.179640_WP; names( 5) = 'RHO'
  data(:,:,:, 6) = 0.0_WP;      names( 6) = 'dRHO'
  data(:,:,:, 7) = 0.0_WP;      names( 7) = 'VISC'
  data(:,:,:, 8) = 0.0_WP;      names( 8) = 'DIFF'
  data(:,:,:, 9) = 0.0_WP;      names( 9) = 'ZMIX'
  data(:,:,:,10) = 0.0_WP;      names(10) = 'ZMIX2'
  data(:,:,:,11) = 0.0_WP;      names(11) = 'PROG'
  
  return
end subroutine sydney_data

