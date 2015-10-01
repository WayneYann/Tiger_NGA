module glass
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  integer :: dim
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module glass

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine glass_grid
  use glass
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  Lx=0.3_WP;Ly=0.3_WP;Lz=0.3_WP

  ! Find dimensionality
  dim=3
  if (nz.eq.1) then
     dim=2
     Lz=Lx/real(nx,WP)
  end if
  
  ! Cartesian
  icyl = 0
  
  ! Periodicity
  xper=0
  yper=0
  zper=1

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
     z(k) = real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
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
  mask(:, 1) = 1
  mask(:,ny) = 1
  mask(1 ,:) = 1
  mask(nx,:) = 1
  ! Jet opening
  mask(1 ,6*ny/8:7*ny/8)=0
  mask(nx,3*ny/8:ny-1)=0
  ! Glass
  mask(nx/2:nx/2+1,ny/8:4*ny/8)=1
  mask(3*nx/4:3*nx/4+1,ny/8:4*ny/8)=1
  mask(nx/2:3*nx/4,ny/8:ny/8+1)=1

  return
end subroutine glass_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine glass_data
  use glass
  use parser
  implicit none
  
  integer :: i,j,k

  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  
  ! Initialize
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  G = -huge(1.0_WP)

  ! Generate small chunk of liquid around the jet
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k)=0.02_WP-sqrt((ym(j)-0.242_WP)**2+xm(i)**2)
        end do
     end do
  end do

  return
end subroutine glass_data
