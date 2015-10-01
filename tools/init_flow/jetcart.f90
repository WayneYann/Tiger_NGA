module jetcart
  use precision
  use param
  implicit none

  ! Length, height, and width of the combustor
  real(WP) :: length,height,width
  ! Enclosed or not
  logical :: enclosed
  ! Height of the channel
  real(WP) :: channel_height
  ! Thickness of the wall/number of points in wall
  real(WP) :: wall_thickness
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in each part
  integer :: channel_points,wall_points
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX
  
end module jetcart

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine jetcart_grid
  use jetcart
  use parser
  implicit none
  
  integer :: i,j,k
  real(WP), parameter :: piovertwo = 2.0_WP*atan(1.0_WP)
  real(WP) :: s,dytmp
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor full height',height)
  call parser_read('Combustor full width',width)
  
  call parser_read('Channel full height',channel_height)
  call parser_read('Wall thickness',wall_thickness)
  call parser_read('Inflow length',inflow_length)

  call parser_read('Channel points',channel_points)
  call parser_read('Stretching',s)
  call parser_read('Wall points',wall_points)
  
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
     x(i) = (i-1)*(length+inflow_length)/nx - inflow_length
  end do
  ! Central jet
  do j=ny/2+2,ny/2+channel_points/2+1
     y(j) = channel_height/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1-channel_points/2),WP)/real(channel_points,WP)-1.0_WP))/tanh(s)
  end do
  ! Wall
  do j=ny/2+channel_points/2+2,ny/2+channel_points/2+wall_points+1
     y(j) = y(j-1) + wall_thickness/real(wall_points,WP)
  end do
  ! Coflow
  dytmp = y(ny/2+channel_points/2+wall_points+1)
  s = log((wall_thickness/real(wall_points,WP))/(height/2.0_WP-channel_height/2.0_WP-wall_thickness))/log(1.0_WP/real(ny/2-channel_points/2-wall_points,WP))
  do j=ny/2+channel_points/2+wall_points+2,ny+1
     y(j) = dytmp+(height/2.0_WP-dytmp)*(real(j-ny/2-channel_points/2-wall_points-1,WP)/real(ny/2-channel_points/2-wall_points,WP))**s
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
        if (y(j).lt.-0.5_WP*channel_height .and. &
            y(j+wall_points).gt.-0.5_WP*channel_height .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
     do j=ny/2+1,ny
        if (y(j).gt.0.5_WP*channel_height .and. &
            y(j-wall_points).lt.0.5_WP*channel_height .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j-1) = 1
        end if
     end do
  end do
  
  return
end subroutine jetcart_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine jetcart_data
  use jetcart
  use parser
  implicit none

  character(str_medium) :: data_type
  real(WP) :: rho_init

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
  case default
     print*, "Data type not recognized..."
  end select

  return
end subroutine jetcart_data

