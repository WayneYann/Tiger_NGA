module round_cart_jet
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module round_cart_jet


! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine round_cart_jet_grid
  use round_cart_jet
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Domain definition
  icyl=0
  xper=0
  yper=1
  zper=1
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('Lx',Lx)
  dx = Lx/real(nx,WP)
  call parser_read('ny',ny)
  call parser_read('Ly',Ly)
  dy = Ly/real(ny,WP)
  call parser_read('nz',nz)
  call parser_read('Lz',Lz)
  dz = Lz/real(nz,WP)
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*dy-0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*dz-0.5_WP*Lz
  end do
  
  ! Create the cell centered grid
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
end subroutine round_cart_jet_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine round_cart_jet_data
  use round_cart_jet
  use parser
  implicit none
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Set the initial fields
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  
  return
end subroutine round_cart_jet_data
