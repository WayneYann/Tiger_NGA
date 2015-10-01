module uniform
  use precision
  use param
  implicit none

  ! Length of the domain
  real(WP) :: Lx,Ly,Lz

  ! Mean velocities
  real(WP) :: Umean,Vmean,Wmean

  ! Pointers to variables in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: T

end module uniform


subroutine uniform_grid
  use uniform
  use parser
  implicit none

  integer :: i,j,k

  ! Read the size of the domain
  call parser_read('nx',nx)
  call parser_read('Lx',Lx)
  call parser_read('ny',ny)
  call parser_read('Ly',Ly)
  call parser_read('nz',nz)
  call parser_read('Lz',Lz)
  
  ! Set the periodicity
!!$  xper = 1
!!$  yper = 1
!!$  zper = 1
  
  call parser_read('xper',xper)
  call parser_read('yper',yper)
  call parser_read('zper',zper)

  ! Cartesian
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)
  end do

  ! Create the mid points
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do

  ! Create the masks
  mask = 0

  return
end subroutine uniform_grid


subroutine uniform_data
  use uniform
  use parser
  implicit none

  integer :: i,j,k

  ! Read the means
  call parser_read('Mean U Velocity',Umean)
  call parser_read('Mean V Velocity',Vmean)
  call parser_read('Mean W Velocity',Wmean)

  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))

  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  T => data(:,:,:,4); names(4) = 'T'

  ! Create them
  U = Umean
  V = Vmean
  W = Wmean
  T = 298.0_WP

  return
end subroutine uniform_data
