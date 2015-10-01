module jetcyl
  use precision
  use param
  implicit none
  
  ! Length and diameter of the combustor
  real(WP) :: length, diameter
  ! Enclosed or not
  logical :: enclosed
  ! Diameter of the pipe
  real(WP) :: pipe_diameter
  ! Inner/outer diameter of the annulus
  real(WP) :: ann_inner, ann_outer
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in radial direction
  integer :: pipe_points, wall_points, ann_points
  ! Radial stretching
  real(WP) :: pipe_stretch,ann_stretch,outer_stretch
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX

end module jetcyl

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine jetcyl_grid
  use jetcyl
  use parser
  implicit none

  integer :: i,j,k
  real(WP), parameter :: twoPi = 2.0_WP*acos(-1.0_WP)

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor diameter',diameter)

  call parser_read('Enclosed',enclosed)

  call parser_read('Pipe diameter',pipe_diameter)
  call parser_read('Ann. inner diameter',ann_inner)
  call parser_read('Ann. outer diameter',ann_outer)

  call parser_read('Inflow length',inflow_length)

  call parser_read('Pipe points',pipe_points)
  call parser_read('Pipe stretching',pipe_stretch)
  call parser_read('Wall points',wall_points)
  call parser_read('Ann. points',ann_points)
  call parser_read('Ann. stretching',ann_stretch)
  call parser_read('Outer stretching',outer_stretch)

  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 1

  ! If enclosed, add the extra point
  if (enclosed) ny = ny + 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid  -- SOLVE FOR STRETCHING FACTORS?
  do i=1,nx+1
     x(i) = real(i-1,WP)*(length+inflow_length)/real(nx,WP) - inflow_length
  end do
  ! Central jet
  do j=1,pipe_points+1
     y(j) = 0.5_WP*pipe_diameter*tanh(pipe_stretch*real(j-1,WP)/real(pipe_points,WP))/tanh(pipe_stretch)
!!$     print*, 'pipe', j, y(j)
  end do
  ! Wall
  do j=pipe_points+2,pipe_points+wall_points
     y(j) = y(j-1) + 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)
!!$     print*, 'wall', j, y(j)
  end do
  ! Annulus
  do j=pipe_points+wall_points+1,pipe_points+wall_points+ann_points+1
     y(j) = 0.25_WP*(ann_outer-ann_inner)*tanh(ann_stretch*(2.0_WP*real(j-pipe_points-wall_points-1,WP)/real(ann_points,WP)-1.0_WP))/tanh(ann_stretch) + 0.25_WP*(ann_outer+ann_inner)
!!$     print*, 'annulus', j, y(j)
  end do
  ! Remainder
  ! SD
  if (pipe_points+wall_points+ann_points+1.lt.ny+1) then
     if (enclosed) then
        do j=pipe_points+wall_points+ann_points+1,ny
           y(j) = 0.5_WP*(diameter-ann_outer)*tanh(outer_stretch*(real(j-pipe_points-wall_points-ann_points-1,WP)/real(ny-1-ann_points-wall_points-pipe_points,WP)-1.0_WP))/tanh(outer_stretch) + 0.5_WP*diameter
!!$        print*, 'outside', j, y(j)
        end do
        y(ny+1) = 2.0_WP*y(ny) - y(ny-1)
     else
        do j=pipe_points+wall_points+ann_points+1,ny+1
           y(j) = 0.5_WP*(diameter-ann_outer)*tanh(outer_stretch*(real(j-pipe_points-wall_points-ann_points-1,WP)/real(ny-ann_points-wall_points-pipe_points,WP)-1.0_WP))/tanh(outer_stretch) + 0.5_WP*diameter
!!$        print*, 'outside', j, y(j)
        end do
     end if
  end if
  ! SD END

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
  do i=1,nx
     do j=1,ny
        if (ym(j).gt.0.5_WP*pipe_diameter .and. &
            ym(j).lt.0.5_WP*ann_inner .and. &
            xm(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
        if (ym(j).gt.0.5_WP*ann_outer .and. &
            xm(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
  end do

  ! If enclosed, add walls
  if (enclosed) mask(:,ny) = 1
 
  return
end subroutine jetcyl_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine jetcyl_data
  use jetcyl
  use parser
  implicit none

  character(str_medium) :: data_type
  real(WP) :: rho_init, T_init, U_init, V_init, W_init, P_init
  real(WP), parameter :: R_cst   = 8314.34_WP   ! J/[kmol K]

  ! Data type
  call parser_read('Data type',data_type,'cold')

  select case(trim(adjustl(data_type)))
  case ('cold')
     ! Allocate the array data
     nvar = 4
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U => data(:,:,:,1); names(1) = 'U'
     V => data(:,:,:,2); names(2) = 'V'
     W => data(:,:,:,3); names(3) = 'W'
     P => data(:,:,:,4); names(4) = 'P'
  
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
  case ('passive mixing')
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
  case ('cold mixing')
     ! Allocate the array data
     nvar = 7
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,7); names(7) = 'ZMIX'
  
     ! Create them
     U    = 0.0_WP
     V    = 0.0_WP
     W    = 0.0_WP
     P    = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
     ZMIX = 0.0_WP
  case ('finite chem')
     ! Allocate the array data
     nvar = 49 ! 7 + N species + T + two scalars to watch the flow
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,47); names(47) = 'ZMIX'
     names(7) = 'N2'
     names(8) = 'H'
     names(9) = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'HO2'
     names(15) = 'H2O2'
     names(16) = 'CO'
     names(17) = 'CO2'
     names(18) = 'HCO'
     names(19) = 'CH3'
     names(20) = 'CH4'
     names(21) = 'CH2O'
     names(22) = 'CH3O'
     names(23) = 'C2H6'
     names(24) = 'CH2OH'
     names(25) = 'C2H5'
     names(26) = 'CH2'
     names(27) = 'CH2X'
     names(28) = 'C2H4'
     names(29) = 'CH3HCO'
     names(30) = 'C2H2'
     names(31) = 'C2H3'
     names(32) = 'CH3OCH3'
     names(33) = 'CH3OCH2'
     names(34) = 'CH3OCH2O'
     names(35) = 'CH3OCHO'
     names(36) = 'CH3OCO'
     names(37) = 'RO2'
     names(38) = 'ROH'
     names(39) = 'QOOH'
     names(40) = 'O2QOOH'
     names(41) = 'HO2QHO'
     names(42) = 'OCH2OCHO'
     names(43) = 'HOCH2OCO'
     names(44) = 'HOCH2O'
     names(45) = 'HCOOH'
     names(46) = 'T'
     names(48) = 'Z_inner'
     names(49) = 'Z_outer'
  
     ! Create them
     call parser_read('U',U_init)
     call parser_read('V',V_init)
     call parser_read('W',W_init)
     call parser_read('Pressure',P_init)
     U = U_init
     V = V_init
     W = W_init
     P = P_init
     data(:,:,:,7) = 0.768_WP
     data(:,:,:,8) = 0.0_WP
     data(:,:,:,9) = 0.232_WP
     data(:,:,:,10) = 0.0_WP
     data(:,:,:,11) = 0.0_WP
     data(:,:,:,12) = 0.0_WP
     data(:,:,:,13) = 0.0_WP
     data(:,:,:,14) = 0.0_WP
     data(:,:,:,15) = 0.0_WP
     data(:,:,:,16) = 0.0_WP
     data(:,:,:,17) = 0.0_WP
     data(:,:,:,18) = 0.0_WP
     data(:,:,:,19) = 0.0_WP
     data(:,:,:,20) = 0.0_WP
     data(:,:,:,21) = 0.0_WP
     data(:,:,:,22) = 0.0_WP
     data(:,:,:,23) = 0.0_WP
     data(:,:,:,24) = 0.0_WP
     data(:,:,:,25) = 0.0_WP
     data(:,:,:,26) = 0.0_WP
     data(:,:,:,27) = 0.0_WP
     data(:,:,:,28) = 0.0_WP
     data(:,:,:,29) = 0.0_WP
     data(:,:,:,30) = 0.0_WP
     data(:,:,:,31) = 0.0_WP
     data(:,:,:,32) = 0.0_WP
     data(:,:,:,33) = 0.0_WP
     data(:,:,:,34) = 0.0_WP
     data(:,:,:,35) = 0.0_WP
     data(:,:,:,36) = 0.0_WP
     data(:,:,:,37) = 0.0_WP
     data(:,:,:,38) = 0.0_WP
     data(:,:,:,39) = 0.0_WP
     data(:,:,:,40) = 0.0_WP
     data(:,:,:,41) = 0.0_WP
     data(:,:,:,42) = 0.0_WP
     data(:,:,:,43) = 0.0_WP
     data(:,:,:,44) = 0.0_WP
     data(:,:,:,45) = 0.0_WP
     data(:,:,:,48) = 0.0_WP
     data(:,:,:,49) = 0.0_WP

     call parser_read('Initial temperature',T_init)
     data(:,:,:,46) = T_init
     RHO  = P * 28.84_WP / (T_init * R_cst)
     dRHO = 0.0_WP
     ZMIX = 0.0_WP 
  case default
     print*, "Data type not recognized..."
  end select
  
  return
end subroutine jetcyl_data

