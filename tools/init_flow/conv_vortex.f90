module conv_vortex
  use string
  use precision
  use param
  implicit none
  
  ! Size of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  real(WP) :: sx,sy
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module conv_vortex

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine conv_vortex_grid
  use conv_vortex
  use parser
  use math
  use random
  implicit none
  
  integer  :: i,j,k
  real(WP) :: rand
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Lz',Lz)
  
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
  dx = Lx/real(nx,WP)
  dy = Ly/real(ny,WP)
  dz = Lz/real(nz,WP)
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx-0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*dy-0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*dz-0.5_WP*Lz
  end do
  
  ! Apply stretching
  if (sx.gt.0.0_WP) then
     x = x + sx*sin(1.0_WP*twoPi*x/Lx)
  else if (sx.lt.0.0_WP) then
     call random_init
     do i=2,nx+1
        call random_number(rand)
        rand=1.0_WP+sx-2.0_WP*sx*rand
        x(i) = x(i-1)+rand*(x(i)-x(i-1))
     end do
     x=x+0.5_WP*Lx
     rand=x(nx+1)
     x=x*Lx/rand
     x=x-0.5_WP*Lx
  end if
  if (sy.gt.0.0_WP) then
     y = y + sy*sin(1.0_WP*twoPi*y/Ly)
  else if (sy.lt.0.0_WP) then
     call random_init
     do j=2,ny+1
        call random_number(rand)
        rand=1.0_WP+sy-2.0_WP*sy*rand
        y(j) = y(j-1)+rand*(y(j)-y(j-1))
     end do
     y=y+0.5_WP*Ly
     rand=y(ny+1)
     y=y*Ly/rand
     y=y-0.5_WP*Ly
  end if
  
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
end subroutine conv_vortex_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine conv_vortex_data
  use conv_vortex
  use parser
  implicit none
  
  integer :: i,j,k,dim
  real(WP) :: tau,Umax,A
  character(len=str_short) :: direction
  real(WP), dimension(3) :: Uc
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Set be convective velocity
  call parser_getsize('Convective velocity',dim)
  if (dim.ne.3) stop 'Convective velocity should be of size 3'
  call parser_read('Convective velocity',Uc)
  U = Uc(1)
  V = Uc(2)
  W = Uc(3)
  P = 0.0_WP
  
  ! Create the vortex
  call parser_read('Vortex orientation',direction)
  call parser_read('Vortex size',tau)
  call parser_read('Vortex intensity',Umax)
  A = sqrt(2.0_WP)*Umax*exp(0.5_WP)/tau
  select case (trim(direction))
  case('x')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              V(i,j,k) = V(i,j,k) + A*exp(-(y(j)**2+zm(k)**2)/tau**2)*(-zm(k))
              W(i,j,k) = W(i,j,k) + A*exp(-(ym(j)**2+z(k)**2)/tau**2)*(ym(j))
           end do
        end do
     end do
  case('y')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) = U(i,j,k) + A*exp(-(x(i)**2+zm(k)**2)/tau**2)*(-zm(k))
              W(i,j,k) = W(i,j,k) + A*exp(-(xm(i)**2+z(k)**2)/tau**2)*(xm(i))
           end do
        end do
     end do
  case('z')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) = U(i,j,k) + A*exp(-(x(i)**2+ym(j)**2)/tau**2)*(-ym(j))
              V(i,j,k) = V(i,j,k) + A*exp(-(xm(i)**2+y(j)**2)/tau**2)*(xm(i))
           end do
        end do
     end do
  end select
  
  return
end subroutine conv_vortex_data

