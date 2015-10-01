module electrospray
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
  real(WP), dimension(:,:,:), pointer :: G
  
end module electrospray

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine electrospray_grid
  use electrospray
  use parser
  use math
  implicit none
  integer :: i,j,k
  real(WP) :: NL,NH
  
  ! Cylindrical
  icyl=0
  
  ! Periodicity
  xper=0
  yper=0
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
  
  ! Get dimensions
  call parser_read('Nozzle length',NL)
  call parser_read('Nozzle height',NH)

  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx-NL
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
  do i=1,nx-1
     do j=1,ny
        if ((abs(y(j)).ge.0.5_WP*NH).and.(x(i+1).le.0.0_WP)) then
           mask(i,j) = 1
        end if
     end do
  end do
  mask(:,1 ) = 1
  mask(:,ny) = 1
  mask(nx,:) = mask(nx-1,:)
  
  return
end subroutine electrospray_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine electrospray_data
  use electrospray
  use parser
  use math
  use string
  implicit none
   
  integer :: i,j,k
  real(WP) :: NL

  call parser_read('Nozzle length',NL)

  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  
  ! Set the initial fields
  U=0.0_WP
  V=0.0_WP
  W=0.0_WP
  G=0.0_WP
  do i=1,nx
     do j=1,ny
        do k=1,nz
           G(i,j,k) = -xm(i)-0.5_WP*NL
        end do
     end do
  end do
  
  return
end subroutine electrospray_data
