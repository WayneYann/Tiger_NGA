module conv_scalar
  use string
  use precision
  use param
  implicit none
  
  ! Size of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX

end module conv_scalar

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine conv_scalar_grid
  use conv_scalar
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Lz',Lz)
  
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
  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  do i=1,nx+1
     x(i) = (i-nx/2-1)*dx
  end do
  do j=1,ny+1
     y(j) = (j-ny/2-1)*dy
  end do
  do k=1,nz+1
     z(k) = (k-nz/2-1)*dz
  end do
  
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
end subroutine conv_scalar_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine conv_scalar_data
  use conv_scalar
  use parser
  implicit none
  
  integer :: i,j,k,dim
  real(WP) :: tau
  character(len=str_short) :: direction
  real(WP), dimension(3) :: Uc
  
  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  
  ! Set be convective velocity
  call parser_getsize('Convective velocity',dim)
  if (dim.ne.3) stop 'Convective velocity should be of size 3'
  call parser_read('Convective velocity',Uc)
  U = Uc(1)
  V = Uc(2)
  W = Uc(3)
  P = 0.0_WP
  ZMIX = 0.0_WP
  
  ! Create the scalar
  call parser_read('Scalar orientation',direction)
  call parser_read('Scalar size',tau)
  select case (trim(direction))
  case('x')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = exp(-(ym(j)**2+zm(k)**2)/tau**2)
           end do
        end do
     end do
  case('y')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = exp(-(xm(i)**2+zm(k)**2)/tau**2)
           end do
        end do
     end do
  case('z')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ZMIX(i,j,k) = exp(-(xm(i)**2+ym(j)**2)/tau**2)
           end do
        end do
     end do
  end select
  
  return
end subroutine conv_scalar_data

