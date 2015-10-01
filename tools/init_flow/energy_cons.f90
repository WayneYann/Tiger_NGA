module energy_cons
  use precision
  use param
  implicit none
  
  ! Length of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  real(WP) :: sx,sy
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module energy_cons


! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine energy_cons_grid
  use energy_cons
  use math
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Cylindrical or cartesian
  call parser_read('icyl',icyl)
  
  ! Establish periodicity
  xper = 1
  if (icyl.eq.1) then
     yper = 0
  else
     yper = 1
  end if
  !yper = 0
  zper = 1
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  if (icyl.eq.1) then
     Lz = twoPi
  else
     call parser_read('Lz',Lz)
  end if
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  dx = Lx/real(nx,WP)
  dy = Ly/real(ny,WP)
  dz = Lz/real(nz,WP)
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
  call parser_read('Stretching x',sx)
  call parser_read('Stretching y',sy)
  x = x + sx*sin(twoPi*x/Lx)
  y = y + sy*sin(twoPi*y/Ly)
  
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
  !mask(:,1) = 1
  !mask(:,ny) = 1
  
  return
end subroutine energy_cons_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine energy_cons_data
  use energy_cons
  use random
  use parser
  implicit none
  
  integer  :: i,j,k
  real(WP) :: Uamp,rand
  
  ! Initialize random numbers
  call random_init
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Initialize with random numbers
  P = 0.0_WP
  call parser_read('Uamp',Uamp)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           call random_number(rand)
           rand = (2.0_WP*rand-1.0_WP)*Uamp
           U(i,j,k) = rand
           call random_number(rand)
           rand = (2.0_WP*rand-1.0_WP)*Uamp
           V(i,j,k) = rand
           call random_number(rand)
           rand = (2.0_WP*rand-1.0_WP)*Uamp
           W(i,j,k) = rand
        end do
     end do
  end do
    
  return
end subroutine energy_cons_data
