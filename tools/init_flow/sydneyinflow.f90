module sydneyinflow
  use precision
  use param
  implicit none
  
end module sydneyinflow

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine sydneyinflow_grid
  use sydneyinflow
  use parser
  use math
  implicit none

  integer :: i,j,k
  integer :: nx_int, nx_inf
  real(WP) :: jet_diameter, ann_diameter, wall_thickness
  real(WP) :: recess
  real(WP) :: height

  ! Read in the size of the domain
  call parser_read('nx',nx_int)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  ! Recess
  call parser_read('Recess',recess)
  if (recess.ne.0.1_WP .and. (recess.ne.0.3_WP .and. recess.ne.0.075_WP)) then
     write(*,'(a)') 'Recess must be 100mm or 300mm'
  end if
  
  ! Add points for inflow
  if (recess.eq.0.1_WP) then
     nx_inf = 3*nx_int/40
     nx = nx_int + nx_inf
  elseif (recess.eq.0.3_WP) then
     nx_inf = nx_int/40
     nx = nx_int + nx_inf
  elseif (recess.eq.0.075_WP) then
     nx_inf = 4*nx_int/40
     nx = nx_int + nx_inf
  end if

  ! Account for outer wall
  ny = ny + 1

  ! Geometrical parameters
  jet_diameter   = 0.0040_WP
  ann_diameter   = 0.0075_WP
  wall_thickness = 0.00025_WP

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
  print*, nx_inf, nx_int
  do i=1,nx+1
     x(i) = recess/real(nx_int,WP)*real(i-1,WP) - real(nx_inf,WP)/real(nx_int,WP)*recess
  end do

  do i=1,nx+1
     print*, i,x(i)
  end do

  ! Create the y-grid
  if (ny-1.eq.160) then
     ! Fuel Jet
     do j=1,65
        y(j) = jet_diameter/2.0_WP*tanh(2.0_WP*real(j-1,WP)/64.0_WP)/tanh(2.0_WP)
     end do
     ! Wall
     do j=65,97
        y(j) = wall_thickness/32.0_WP*real(j-65,WP) + jet_diameter/2.0_WP
     end do
     ! Air annulus
     do j=97,161
        height = ann_diameter/2.0_WP - (jet_diameter/2.0_WP+wall_thickness)
        y(j) = 0.5_WP*height*tanh(2.0_WP*(2.0_WP*real(j-97,WP)/64.0_WP-1.0_WP))/tanh(2.0_WP) + &
            0.5_WP*jet_diameter+wall_thickness + 0.5_WP*height
     end do
     y(ny+1) = 2.0_WP*y(ny)-y(ny-1)
  else
     write(*,'(a)') 'ny must be 160'
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
  mask(1:nx_inf,65:96) = 1
  mask(:,ny) = 1
  
  return
end subroutine sydneyinflow_grid


! =============================== !
! Create the variable array: Data !
! =============================== !
subroutine sydneyinflow_data
  use sydneyinflow
  implicit none
  
  ! Allocate the array data
  nvar = 10
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
  
  return
end subroutine sydneyinflow_data

