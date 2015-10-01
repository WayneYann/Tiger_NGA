module taylor
  use string
  use precision
  use param
  implicit none
  
  ! Mesh
  real(WP) :: dx,dy,dz
  real(WP) :: sx,sy
  
  ! Orientation
  character(len=str_medium) :: orientation
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  
end module taylor


! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine taylor_grid
  use taylor
  use math
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read the number of points
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  
  ! Read the stretching
  call parser_read('Stretching x',sx)
  call parser_read('Stretching y',sy)
  
  ! Set the periodicity
  xper = 1
  yper = 1
  zper = 1
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  dx = twoPi/real(nx,WP)
  dy = twoPi/real(ny,WP)
  dz = twoPi/real(nz,WP)
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*dy
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*dz
  end do
  
  ! Apply stretching
  x = x + sx*sin(x)
  y = y + sy*sin(y)
  
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
end subroutine taylor_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine taylor_data
  use taylor
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Read in the orientation
  call parser_read('Orientation',orientation)
  
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
  
  ! Create the fields
  select case (trim(orientation))
  case('x')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) = 0.0_WP
              V(i,j,k) =  cos(y (j))*sin(zm(k))
              W(i,j,k) = -sin(ym(j))*cos(z (k))
              P(i,j,k) = 0.5_WP*(sin(ym(j))**2+sin(zm(k))**2)
              ZMIX(i,j,k) = 0.5_WP*(1.0_WP+sin(ym(j)))
           end do
        end do
     end do
  case('y')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) = -sin(zm(k))*cos(x (i))
              V(i,j,k) = 0.0_WP
              W(i,j,k) =  cos(z (k))*sin(xm(i))
              P(i,j,k) = 0.5_WP*(sin(zm(k))**2+sin(xm(i))**2)
              ZMIX(i,j,k) = 0.5_WP*(1.0_WP+sin(zm(k)))
           end do
        end do
     end do
  case('z')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) =  cos(x (i))*sin(ym(j))
              V(i,j,k) = -sin(xm(i))*cos(y (j))
              W(i,j,k) = 0.0_WP
              P(i,j,k) = 0.5_WP*(sin(xm(i))**2+sin(ym(j))**2)
              ZMIX(i,j,k) = 0.5_WP*(1.0_WP+sin(xm(i)))
           end do
        end do
     end do
  end select
  
  return
end subroutine taylor_data

