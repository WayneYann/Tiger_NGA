module jetcart2
  use precision
  use param
  implicit none

  ! Length, height, and width of the combustor
  real(WP) :: length,height,width
  ! Enclosed or not
  logical :: enclosed
  ! Height of the channel
  real(WP) :: channel_height_1,channel_height_2
  ! Thickness of the wall/number of points in wall
  real(WP) :: wall_thickness_1,wall_thickness_2
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in each part
  integer :: channel_points_1,channel_points_2
  integer :: wall_points_1,wall_points_2
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX
  
end module jetcart2

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine jetcart2_grid
  use jetcart2
  use parser
  implicit none
  
  integer :: i,j,k
  real(WP), parameter :: piovertwo = 2.0_WP*atan(1.0_WP)
  real(WP) :: s,dytmp,sx,i0
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor full height',height)
  call parser_read('Combustor full width',width)
  
  call parser_read('Channel 1 full height',channel_height_1)
  call parser_read('Channel 2 full height',channel_height_2)
  call parser_read('Wall 1 thickness',wall_thickness_1)
  call parser_read('Wall 2 thickness',wall_thickness_2)
  call parser_read('Inflow length',inflow_length)

  call parser_read('Channel 1 points',channel_points_1)
  call parser_read('Channel 2 points',channel_points_2)
  call parser_read('Wall 1 points',wall_points_1)
  call parser_read('Wall 2 points',wall_points_2)

  call parser_read('y-Stretching',s)
  call parser_read('x-Stretching',sx)
  
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
  i0 = nx*inflow_length**(1/s)/(length**(1/s)+inflow_length**(1/s))
  do i=0,nx
     if (i.le.i0) then
        x(i+1) = -(i0-i)**s*inflow_length/i0**s
     else
        x(i+1) = (i-i0)**s*length/(nx-i0)**s
     end if
!!$  do i=1,nx+1
!!$     x(i) = (i-1)*(length+inflow_length)/nx - inflow_length
  end do

  ! Central jet
  do j=ny/2+2,ny/2+channel_points_1/2+1
     y(j) = channel_height_1/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1-channel_points_1/2),WP)/real(channel_points_1,WP)-1.0_WP))/tanh(s)
  end do
  ! Wall around central jet
  do j=ny/2+channel_points_1/2+2,ny/2+channel_points_1/2+wall_points_1+1
     y(j) = y(j-1) + wall_thickness_1/real(wall_points_1,WP)
  end do
  ! Auxiliary jet
  do j=ny/2+channel_points_1/2+wall_points_1+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+1
     y(j) = channel_height_2/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1+channel_points_1/2+wall_points_1+channel_points_2/2))/real(channel_points_2,WP)))/tanh(s)+channel_height_1/2.0_WP+wall_thickness_1+channel_height_2/2.0_WP
  end do
  ! Wall around aux jet
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1
     y(j) = y(j-1) + wall_thickness_2/real(wall_points_2,WP)
  end do
  ! Coflow
  dytmp = y(ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1)
  s = log((wall_thickness_2/real(wall_points_2,WP))/(height/2.0_WP-channel_height_1/2.0_WP-wall_thickness_1-channel_height_2-wall_thickness_2))/log(1.0_WP/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+2,ny+1
     y(j) = dytmp+(height/2.0_WP-dytmp)*(real(j-ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2-1,WP)/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))**s
  end do
  ! Set mid-point to zero
  y(ny/2+1) = 0.0_WP
  ! Mirror the y-grid
  do j=ny/2+2,ny+1
     y(ny+2-j) = -y(j)
  end do
  do k=1,nz+1
     z(k) = (k-1)*width/nz - 0.5_WP*width
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
     do j=1,ny/2
        if (y(j).lt.-0.5_WP*channel_height_1 .and. &
            y(j+wall_points_1).ge.-0.5_WP*channel_height_1 .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j) = 1
        end if
        if (y(j).lt.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             y(j+wall_points_2).ge.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             x(i+1).le.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
     do j=ny/2+1,ny
        if (y(j).gt.0.5_WP*channel_height_1 .and. &
            y(j-wall_points_1).le.0.5_WP*channel_height_1 .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j-1) = 1
        end if
        if (y(j).gt.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             y(j-wall_points_2).le.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             x(i+1).le.0.0_WP) then
           mask(i,j-1) = 1
        end if
     end do
  end do
  
  return
end subroutine jetcart2_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine jetcart2_data
  use jetcart2
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
!!$     nvar = 47 ! 7 + N species + T 
     nvar = 17
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
!!$     ZMIX => data(:,:,:,47); names(47) = 'ZMIX'
     ZMIX => data(:,:,:,17); names(17) = 'ZMIX'
     names(7)  = 'N2'
     names(8)  = 'H'
     names(9)  = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'HO2'
     names(15) = 'H2O2'
!!$     names(16) = 'CO'
!!$     names(17) = 'CO2'
!!$     names(18) = 'HCO'
!!$     names(19) = 'CH3'
!!$     names(20) = 'CH4'
!!$     names(21) = 'CH2O'
!!$     names(22) = 'CH3O'
!!$     names(23) = 'C2H6'
!!$     names(24) = 'CH2OH'
!!$     names(25) = 'C2H5'
!!$     names(26) = 'CH2'
!!$     names(27) = 'CH2X'
!!$     names(28) = 'C2H4'
!!$     names(29) = 'CH3HCO'
!!$     names(30) = 'C2H2'
!!$     names(31) = 'C2H3'
!!$     names(32) = 'CH3OCH3'
!!$     names(33) = 'CH3OCH2'
!!$     names(34) = 'CH3OCH2O'
!!$     names(35) = 'CH3OCHO'
!!$     names(36) = 'CH3OCO'
!!$     names(37) = 'RO2'
!!$     names(38) = 'ROH'
!!$     names(39) = 'QOOH'
!!$     names(40) = 'O2QOOH'
!!$     names(41) = 'HO2QHO'
!!$     names(42) = 'OCH2OCHO'
!!$     names(43) = 'HOCH2OCO'
!!$     names(44) = 'HOCH2O'
!!$     names(45) = 'HCOOH'
!!$     names(46) = 'T'
     names(16) = 'T'
       
     ! Create them
     call parser_read('U',U_init)
     call parser_read('V',V_init)
     call parser_read('W',W_init)
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     U = U_init
     V = V_init
     W = W_init
     P = 0.0_WP !Hydrodynamic pressure
     data(:,:,:,7)  = 0.768_WP
     data(:,:,:,8)  = 0.0_WP
     data(:,:,:,9)  = 0.232_WP
     data(:,:,:,10) = 0.0_WP
     data(:,:,:,11) = 0.0_WP
     data(:,:,:,12) = 0.0_WP
     data(:,:,:,13) = 0.0_WP
     data(:,:,:,14) = 0.0_WP
     data(:,:,:,15) = 0.0_WP
!!$     data(:,:,:,16) = 0.0_WP
!!$     data(:,:,:,17) = 0.0_WP
!!$     data(:,:,:,18) = 0.0_WP
!!$     data(:,:,:,19) = 0.0_WP
!!$     data(:,:,:,20) = 0.0_WP
!!$     data(:,:,:,21) = 0.0_WP
!!$     data(:,:,:,22) = 0.0_WP
!!$     data(:,:,:,23) = 0.0_WP
!!$     data(:,:,:,24) = 0.0_WP
!!$     data(:,:,:,25) = 0.0_WP
!!$     data(:,:,:,26) = 0.0_WP
!!$     data(:,:,:,27) = 0.0_WP
!!$     data(:,:,:,28) = 0.0_WP
!!$     data(:,:,:,29) = 0.0_WP
!!$     data(:,:,:,30) = 0.0_WP
!!$     data(:,:,:,31) = 0.0_WP
!!$     data(:,:,:,32) = 0.0_WP
!!$     data(:,:,:,33) = 0.0_WP
!!$     data(:,:,:,34) = 0.0_WP
!!$     data(:,:,:,35) = 0.0_WP
!!$     data(:,:,:,36) = 0.0_WP
!!$     data(:,:,:,37) = 0.0_WP
!!$     data(:,:,:,38) = 0.0_WP
!!$     data(:,:,:,39) = 0.0_WP
!!$     data(:,:,:,40) = 0.0_WP
!!$     data(:,:,:,41) = 0.0_WP
!!$     data(:,:,:,42) = 0.0_WP
!!$     data(:,:,:,43) = 0.0_WP
!!$     data(:,:,:,44) = 0.0_WP
!!$     data(:,:,:,45) = 0.0_WP
    
     call parser_read('Initial temperature',T_init)
!!$     data(:,:,:,46) = T_init
     data(:,:,:,16) = T_init
     RHO  = P_init * 28.84_WP / (T_init * R_cst)
     dRHO = 0.0_WP
     ZMIX = 0.0_WP 
  case default
     print*, "Data type not recognized..."
  end select

  return
end subroutine jetcart2_data

