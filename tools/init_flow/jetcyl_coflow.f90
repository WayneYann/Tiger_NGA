module jetcyl_coflow
  use precision
  use param
  implicit none
  
  ! Length and diameter of the combustor
  real(WP) :: length, diameter
  ! Enclosed or not
  logical :: enclosed, xuniform
  ! Diameter of the pipe
  real(WP) :: pipe_diameter
  ! Inner diameter of the annulus
  real(WP) :: ann_inner
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in radial direction
  integer :: pipe_points, wall_points, refine_times
  ! Radial stretching
  real(WP) :: s
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX

end module jetcyl_coflow

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine jetcyl_coflow_grid
  use jetcyl_coflow
  use parser
  implicit none

  integer :: i,j,k,xz_1,zone_points
  real(WP), parameter :: twoPi = 2.0_WP*acos(-1.0_WP)
  real(WP) :: dytmp,zone_begin,zone_length

  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor diameter',diameter)

  call parser_read('Enclosed',enclosed)

  call parser_read('Pipe diameter',pipe_diameter)
  call parser_read('Ann. inner diameter',ann_inner)
  
  call parser_read('Inflow length',inflow_length)

  call parser_read('Pipe points',pipe_points)
  call parser_read('Wall points',wall_points)

  call parser_read('X uniform',xuniform,.false.)
  call parser_read('Refined mesh',refine_times,1)
  

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

!!$  ! Create the grid  -- SOLVE FOR STRETCHING FACTORS?
!!$  do i=1,nx+1
!!$     x(i) = real(i-1,WP)*(length+inflow_length)/real(nx,WP) - inflow_length
!!$  end do
!!$  ! Central jet
!!$  s = log(0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP) / (0.5_WP*pipe_diameter)) / log(1.0_WP/real(pipe_points,WP))
!!$  dytmp = 0.5_wp * pipe_diameter
!!$  do j=1,pipe_points
!!$     y(j) = dytmp - 0.5_WP*pipe_diameter * (real(pipe_points+1-j,WP) / real(pipe_points,WP))**s
!!$     print*, 'pipe', j, y(j)
!!$  end do
!!$  ! Wall
!!$  do j=pipe_points+1,pipe_points+wall_points+1
!!$     y(j) = y(j-1) + 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)
!!$     print*, 'wall', j, y(j)
!!$  end do
!!$  ! Annulus
!!$  dytmp = y(pipe_points+wall_points+1)
!!$  s = log(0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP) / (0.5_WP*(diameter-ann_inner))) / log(1.0_WP/real(ny-wall_points-pipe_points,WP))
!!$  do j=pipe_points+wall_points+2,ny+1
!!$     y(j) = dytmp + 0.5_WP * (diameter-ann_inner) * (real(j-pipe_points-wall_points-1,WP) / real(ny-pipe_points-wall_points,WP))**s
!!$     print*, 'annulus', j, y(j)
!!$  end do
!!$  do k=1,nz+1
!!$     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
!!$  end do
!!$  if (y(2)-y(1) .lt. 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)) print*, 'Use fewer points in the pipe.'
!!$  if (y(ny+1)-y(ny) .lt. 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)) print*, 'Use fewer points in the coflow.'
!!$  print*, 'Stretch rates: pipe and annulus:',( y(pipe_points)-y(pipe_points-1))/(y(pipe_points+1)-y(pipe_points)), (y(pipe_points+wall_points+3)-y(pipe_points+wall_points+2))/(y(pipe_points+1)-y(pipe_points))

! Create the grid  -- SOLVE FOR STRETCHING FACTORS?
  if (.not.xuniform) then
!!$     call parser_read('zone_begin',zone_begin)
!!$     call parser_read('zone_length',zone_length)
!!$     call parser_read('xz_1',xz_1)
!!$     call parser_read('zone_points',zone_points)
!!$     ! Zone_1 including the inflow default: 4 points in the inflow
!!$     s = 0.97
!!$     dytmp = zone_begin / (s**(xz_1) - s**5)
!!$     do i=1,xz_1
!!$        if(i.le.4) then
!!$           x(i) = -(5-i) * (- dytmp * s**5 + dytmp * s**6)
!!$        else
!!$           x(i) = - dytmp * s**5 + dytmp * s**i
!!$        end if
!!$        print*, 'Zone_1', i, x(i)
!!$     end do
!!$     ! Zone_2 of interest--Uniform
!!$     do i=xz_1+1,xz_1+zone_points-1
!!$        x(i) = x(i-1) + (zone_length)/real(zone_points,WP)
!!$        print*, 'Zone of interest', i, x(i)
!!$     end do
!!$     ! Zone_3 to the end of the domain
!!$     s = 1.03
!!$     dytmp = (length - (zone_begin + zone_length)) / (s**(nx+1) - s**(xz_1+zone_points))
!!$     do i=xz_1+zone_points,nx+1
!!$        x(i) = length - dytmp * s**(nx+1) + dytmp * s**i
!!$        print*, 'Zone_3', i, x(i)
!!$     end do
!!$     print*,'Zone_1 (0.97-1.03):', (x(xz_1)-x(xz_1-1))/(zone_length/real(zone_points,WP))
!!$     print*,'Zone_3 (0.97-1.03):', (x(xz_1+zone_points+2)-x(xz_1+zone_points+1))/(zone_length/real(zone_points,WP))
     call parser_read('zone_length',zone_length)
     call parser_read('zone_points',zone_points)
     ! Zone_1 of interest--Uniform
     do i=1,zone_points+1
        x(i) = real(i-refine_times-1,WP) * (zone_length)/real(zone_points-refine_times,WP) ! Double grid points
        !x(i) = real(i-2,WP) * (zone_length)/real(zone_points,WP)
     end do
     
     ! Zone_2 to the end of the domain
     s = 1.03
     dytmp = (length - x(zone_points+1)) / (s**(nx+1) - s**(zone_points+1))
     do i=zone_points+2,nx+1
        x(i) = length - dytmp * s**(nx+1) + dytmp * s**i
        print*, 'Zone_2', i, x(i)
     end do
     print*, 'Zone of interest', x(1), x(zone_points+1)
     print*,'Zone_2 (0.97-1.03):', (x(zone_points+2)-x(zone_points+1))/(zone_length/real(zone_points,WP))
  else
     print*, 'Uniform X is used.'
     do i=1,nx+1
        !x(i) = real(i-1,WP)*(length+inflow_length)/real(nx,WP) - inflow_length
        x(i) = real(i-refine_times-1,WP)*length/real(nx-refine_times,WP) ! for double grid points
        !x(i) = real(i-2,WP)*length/real(nx,WP)
        print*,'x', i, x(i)
     end do
     print*,'X_left and X_right:',x(1),x(nx+1)
  end if
  ! Central jet
  s = 0.97   ! stretch rate y(j) = A + B * S**j dytmp = B  DEFAULT: 0.97
  dytmp = (0.5_WP * pipe_diameter - 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)) / (s**(pipe_points) - s)
  do j=1,pipe_points+1
     y(j) = -dytmp * s + dytmp * s**j
     if (j.eq.pipe_points+1) y(j) = y(j-1) + 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)
     print*, 'pipe', j, y(j)
  end do
  ! Wall
  do j=pipe_points+2,pipe_points+wall_points
     y(j) = y(j-1) + 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)
     print*, 'wall', j, y(j)
  end do
  ! Annulus
  s = 1.03   ! DEFAULT: 1.03
  dytmp = (0.5_WP * (diameter - ann_inner) - 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)) / (s**(ny+1) - s**(pipe_points+wall_points+2))
  do j=pipe_points+wall_points+1,ny+1
     !y(j) = 0.5_WP * ann_inner - dytmp * s**(pipe_points+wall_points+1) + dytmp * s**j
     y(j) = 0.5_WP * diameter - dytmp * s**(ny+1) + dytmp * s**j
     if (j.eq.pipe_points+wall_points+1) y(j) = y(j-1) + 0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP)
     print*, 'annulus', j, y(j)
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do
  print*, 'The stretch rate in the pipe is 0.99'
  print*, 'In the pipe (0.97-1.03):', (y(pipe_points)-y(pipe_points-1))/(0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP))
  print*, 'In the coflow (0.97-1.03):', (y(pipe_points+wall_points+3)-y(pipe_points+wall_points+2))/(0.5_WP*(ann_inner-pipe_diameter)/real(wall_points,WP))

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

!!$  ! Create the masks
!!$  mask = 0
!!$  do i=1,nx
!!$     do j=1,ny
!!$        if (ym(j).gt.0.5_WP*pipe_diameter .and. &
!!$            ym(j).lt.0.5_WP*ann_inner .and. &
!!$            xm(i).lt.0.0_WP) then
!!$           mask(i,j) = 1
!!$        end if
!!$     end do
!!$  end do

  ! Create the masks
  mask = 0
!!$  do i=1,nx
!!$     do j=1,ny
!!$        if (j.gt.pipe_points .and. &
!!$             j.lt.pipe_points+wall_points+1 .and. &
!!$             xm(i).lt.0.0_WP) then
!!$           mask(i,j) = 1
!!$           print*, 'Mask', '(',i,',',j,')'
!!$        end if
!!$     end do
!!$  end do

! SD if double the mesh, make sure that it has two points in the inlet
     do j=1,ny
        if (j.gt.pipe_points .and. &
             j.lt.pipe_points+wall_points+1) then
           mask(1:refine_times,j) = 1
           print*, 'Mask', '(1:',refine_times,',',j,')'  ! For double grid points
        end if
     end do

  ! If enclosed, add walls
  if (enclosed) mask(:,ny) = 1
 
  return
end subroutine jetcyl_coflow_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine jetcyl_coflow_data
  use jetcyl_coflow
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
     nvar = 7 + 158 + 1 !47 ! 7 + N species + T 
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,nvar); names(nvar) = 'ZMIX'
     ! 39 Species DME mech
!!$     names(7) = 'N2'
!!$     names(8) = 'H'
!!$     names(9) = 'O2'
!!$     names(10) = 'O'
!!$     names(11) = 'OH'
!!$     names(12) = 'H2'
!!$     names(13) = 'H2O'
!!$     names(14) = 'HO2'
!!$     names(15) = 'H2O2'
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

     ! 160 Species TheSoot mech
     names(7) = 'N2'
     names(8) = 'H'
     names(9) = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'CO2'
     names(15) = 'HO2'
     names(16) = 'H2O2'
     names(17) = 'CO'
     names(18) = 'HCO'
     names(19) = 'C'
     names(20) = 'CH'
     names(21) = 'TXCH2'
     names(22) = 'CH3'
     names(23) = 'CH2O'
     names(24) = 'HCCO'
     names(25) = 'C2H'
     names(26) = 'CH2CO'
     names(27) = 'C2H2'
     names(28) = 'SXCH2'
     names(29) = 'AR'
     names(30) = 'CH3OH'
     names(31) = 'CH2OH'
     names(32) = 'CH3O'
     names(33) = 'CH4'
     names(34) = 'CH3O2'
     names(35) = 'C2H3'
     names(36) = 'C2H4'
     names(37) = 'C2H5'
     names(38) = 'HCCOH'
     names(39) = 'CH2CHO'
     names(40) = 'CH3CHO'
     names(41) = 'H2C2'
     names(42) = 'C2H5O'
     names(43) = 'NXC3H7'
     names(44) = 'C2H6'
     names(45) = 'C3H8'
     names(46) = 'C3H6'
     names(47) = 'C3H3'
     names(48) = 'PXC3H4'
     names(49) = 'AXC3H4'
     names(50) = 'SXC3H5'
     names(51) = 'NXC4H3'
     names(52) = 'C2H3CHO'
     names(53) = 'AXC3H5'
     names(54) = 'C2O'
     names(55) = 'C4H4'
     names(56) = 'C3H2'
     names(57) = 'C3H2O'
     names(58) = 'C4H2'
     names(59) = 'IXC4H3'
     names(60) = 'TXC3H5'
     names(61) = 'C3H5O'
     names(62) = 'C4H'
     names(63) = 'C8H2'
     names(64) = 'C6H2'
     names(65) = 'C4H6'
     names(66) = 'NXC4H5'
     names(67) = 'IXC4H5'
     names(68) = 'A1XC6H6'
     names(69) = 'NXC7H16'
     names(70) = 'C5H11'
     names(71) = 'PXC4H9'
     names(72) = 'C7H15'
     names(73) = 'PXC4H8'
     names(74) = 'C5H10'
     names(75) = 'C7H14'
     names(76) = 'C7H15O'
     names(77) = 'C3H7CHO'
     names(78) = 'C4H7'
     names(79) = 'C7H13'
     names(80) = 'C5H9'
     names(81) = 'C4H7O'
     names(82) = 'NXC3H7O'
     names(83) = 'IXC8H18'
     names(84) = 'YXC7H15'
     names(85) = 'IXC4H8'
     names(86) = 'IXC3H7'
     names(87) = 'TXC4H9'
     names(88) = 'CXC8H17'
     names(89) = 'YXC7H14'
     names(90) = 'DXC8H17O'
     names(91) = 'CH3COCH3'
     names(92) = 'IXC4H7'
     names(93) = 'XXC7H13'
     names(94) = 'IXC3H5CH'
     names(95) = 'TXC4H9O'
     names(96) = 'IXC4H7O'
     names(97) = 'C5H4CH2'
     names(98) = 'A1XXC6H5'
     names(99) = 'A1C2H2XC'
     names(100) = 'A1C2H3XC'
     names(101) = 'A1C2HXC8'
     names(102) = 'A1C2HYXC'
     names(103) = 'A1C2H3YX'
     names(104) = 'A2XXC10H'
     names(105) = 'A2XC10H8'
     names(106) = 'A2YXC10H'
     names(107) = 'A2C2H2AX'
     names(108) = 'A2C2H2BX'
     names(109) = 'A2C2HAXC'
     names(110) = 'A2C2HBXC'
     names(111) = 'A2C2HAYX'
     names(112) = 'A2C2HBYX'
     names(113) = 'A2R5XC12'
     names(114) = 'A2R5XXC1'
     names(115) = 'A2R5C2H2'
     names(116) = 'A2R5C2HX'
     names(117) = 'A2R5C2HY'
     names(118) = 'P2XC12H1'
     names(119) = 'P2XXC12H'
     names(120) = 'A3XXC14H'
     names(121) = 'A3XC14H1'
     names(122) = 'A3YXC14H'
     names(123) = 'A3R5XXC1'
     names(124) = 'A3R5XC16'
     names(125) = 'A4XC16H1'
     names(126) = 'A4XXC16H'
     names(127) = 'A4R5XC18'
     names(128) = 'FLTNXC16'
     names(129) = 'C5H6'
     names(130) = 'C5H5'
     names(131) = 'TXC5H5O'
     names(132) = 'C5H4O'
     names(133) = 'SXC5H5O'
     names(134) = 'C9H8'
     names(135) = 'C9H7'
     names(136) = 'A1CH2XC7'
     names(137) = 'C9H6O'
     names(138) = 'OXC6H4'
     names(139) = 'A1CH3XC7'
     names(140) = 'A1OHXC6H'
     names(141) = 'HOA1CH3X'
     names(142) = 'OA1CH3XC'
     names(143) = 'A1CH2OXC'
     names(144) = 'A1CH2OHX'
     names(145) = 'A1CHOXC7'
     names(146) = 'A1OXC6H5'
     names(147) = 'A1CH3YXC'
     names(148) = 'A1C2H4XC'
     names(149) = 'A1C2H5XC'
     names(150) = 'C8H9O2'
     names(151) = 'C8H8OOH'
     names(152) = 'OC8H7OOH'
     names(153) = 'A1CH3CH3'
     names(154) = 'A1CH3CH2'
     names(155) = 'A1CH3CHO'
     names(156) = 'A2CH3XC1'
     names(157) = 'A1CHOCH2'
     names(158) = 'A1CHOCHO'
     names(159) = 'A2OHXC10'
     names(160) = 'A2CH2XC1'
     names(161) = 'A2CH2OXC'
     names(162) = 'A2CHOXC1'
     names(163) = 'A2OXC10H'
     names(164) = 'OC6H4O'
     names(nvar-1) = 'T'
       
     ! Create them
     call parser_read('U',U_init)
     call parser_read('V',V_init)
     call parser_read('W',W_init)
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     U = U_init
     V = V_init
     W = W_init
     P = 0.0_WP !Dydrodynamic pressure
     ! 39 Species DME mech
!!$     data(:,:,:,7) = 0.768_WP
!!$     data(:,:,:,8) = 0.0_WP
!!$     data(:,:,:,9) = 0.232_WP
!!$     data(:,:,:,10) = 0.0_WP
!!$     data(:,:,:,11) = 0.0_WP
!!$     data(:,:,:,12) = 0.0_WP
!!$     data(:,:,:,13) = 0.0_WP
!!$     data(:,:,:,14) = 0.0_WP
!!$     data(:,:,:,15) = 0.0_WP
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

     ! 160 Species TheSoot mech
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
     data(:,:,:,46) = 0.0_WP
     data(:,:,:,47) = 0.0_WP
     data(:,:,:,48) = 0.0_WP
     data(:,:,:,49) = 0.0_WP
     data(:,:,:,50) = 0.0_WP
     data(:,:,:,51) = 0.0_WP
     data(:,:,:,52) = 0.0_WP
     data(:,:,:,53) = 0.0_WP
     data(:,:,:,54) = 0.0_WP
     data(:,:,:,55) = 0.0_WP
     data(:,:,:,56) = 0.0_WP
     data(:,:,:,57) = 0.0_WP
     data(:,:,:,58) = 0.0_WP
     data(:,:,:,59) = 0.0_WP
     data(:,:,:,60) = 0.0_WP
     data(:,:,:,61) = 0.0_WP
     data(:,:,:,62) = 0.0_WP
     data(:,:,:,63) = 0.0_WP
     data(:,:,:,64) = 0.0_WP
     data(:,:,:,65) = 0.0_WP
     data(:,:,:,66) = 0.0_WP
     data(:,:,:,67) = 0.0_WP
     data(:,:,:,68) = 0.0_WP
     data(:,:,:,69) = 0.0_WP
     data(:,:,:,70) = 0.0_WP
     data(:,:,:,71) = 0.0_WP
     data(:,:,:,72) = 0.0_WP
     data(:,:,:,73) = 0.0_WP
     data(:,:,:,74) = 0.0_WP
     data(:,:,:,75) = 0.0_WP
     data(:,:,:,76) = 0.0_WP
     data(:,:,:,77) = 0.0_WP
     data(:,:,:,78) = 0.0_WP
     data(:,:,:,79) = 0.0_WP
     data(:,:,:,80) = 0.0_WP
     data(:,:,:,81) = 0.0_WP
     data(:,:,:,82) = 0.0_WP
     data(:,:,:,83) = 0.0_WP
     data(:,:,:,84) = 0.0_WP
     data(:,:,:,85) = 0.0_WP
     data(:,:,:,86) = 0.0_WP
     data(:,:,:,87) = 0.0_WP
     data(:,:,:,88) = 0.0_WP
     data(:,:,:,89) = 0.0_WP
     data(:,:,:,90) = 0.0_WP
     data(:,:,:,91) = 0.0_WP
     data(:,:,:,92) = 0.0_WP
     data(:,:,:,93) = 0.0_WP
     data(:,:,:,94) = 0.0_WP
     data(:,:,:,95) = 0.0_WP
     data(:,:,:,96) = 0.0_WP
     data(:,:,:,97) = 0.0_WP
     data(:,:,:,98) = 0.0_WP
     data(:,:,:,99) = 0.0_WP
     data(:,:,:,100) = 0.0_WP
     data(:,:,:,101) = 0.0_WP
     data(:,:,:,102) = 0.0_WP
     data(:,:,:,103) = 0.0_WP
     data(:,:,:,104) = 0.0_WP
     data(:,:,:,105) = 0.0_WP
     data(:,:,:,106) = 0.0_WP
     data(:,:,:,107) = 0.0_WP
     data(:,:,:,108) = 0.0_WP
     data(:,:,:,109) = 0.0_WP
     data(:,:,:,110) = 0.0_WP
     data(:,:,:,111) = 0.0_WP
     data(:,:,:,112) = 0.0_WP
     data(:,:,:,113) = 0.0_WP
     data(:,:,:,114) = 0.0_WP
     data(:,:,:,115) = 0.0_WP
     data(:,:,:,116) = 0.0_WP
     data(:,:,:,117) = 0.0_WP
     data(:,:,:,118) = 0.0_WP
     data(:,:,:,119) = 0.0_WP
     data(:,:,:,120) = 0.0_WP
     data(:,:,:,121) = 0.0_WP
     data(:,:,:,122) = 0.0_WP
     data(:,:,:,123) = 0.0_WP
     data(:,:,:,124) = 0.0_WP
     data(:,:,:,125) = 0.0_WP
     data(:,:,:,126) = 0.0_WP
     data(:,:,:,127) = 0.0_WP
     data(:,:,:,128) = 0.0_WP
     data(:,:,:,129) = 0.0_WP
     data(:,:,:,130) = 0.0_WP
     data(:,:,:,131) = 0.0_WP
     data(:,:,:,132) = 0.0_WP
     data(:,:,:,133) = 0.0_WP
     data(:,:,:,134) = 0.0_WP
     data(:,:,:,135) = 0.0_WP
     data(:,:,:,136) = 0.0_WP
     data(:,:,:,137) = 0.0_WP
     data(:,:,:,138) = 0.0_WP
     data(:,:,:,139) = 0.0_WP
     data(:,:,:,140) = 0.0_WP
     data(:,:,:,141) = 0.0_WP
     data(:,:,:,142) = 0.0_WP
     data(:,:,:,143) = 0.0_WP
     data(:,:,:,144) = 0.0_WP
     data(:,:,:,145) = 0.0_WP
     data(:,:,:,146) = 0.0_WP
     data(:,:,:,147) = 0.0_WP
     data(:,:,:,148) = 0.0_WP
     data(:,:,:,149) = 0.0_WP
     data(:,:,:,150) = 0.0_WP
     data(:,:,:,151) = 0.0_WP
     data(:,:,:,152) = 0.0_WP
     data(:,:,:,153) = 0.0_WP
     data(:,:,:,154) = 0.0_WP
     data(:,:,:,155) = 0.0_WP
     data(:,:,:,156) = 0.0_WP
     data(:,:,:,157) = 0.0_WP
     data(:,:,:,158) = 0.0_WP
     data(:,:,:,159) = 0.0_WP
     data(:,:,:,160) = 0.0_WP
     data(:,:,:,161) = 0.0_WP
     data(:,:,:,162) = 0.0_WP
     data(:,:,:,163) = 0.0_WP
     data(:,:,:,164) = 0.0_WP


   
     call parser_read('Initial temperature',T_init)
     data(:,:,:,nvar-1) = T_init
     RHO  = P_init * 28.84_WP / (T_init * R_cst)
     dRHO = 0.0_WP
     ZMIX = 0.0_WP 
  case default
     print*, "Data type not recognized..."
  end select
  
  return
end subroutine jetcyl_coflow_data

