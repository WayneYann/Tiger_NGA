module milk_crown
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  integer :: nsect
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module milk_crown

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine milk_crown_grid
  use milk_crown
  use parser
  use math
  implicit none
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  Lz=twoPi
  call parser_read('nsect',nsect,1)
  Lz=Lz/real(nsect,WP)
  xper=0;yper=0;zper=1
  icyl=1
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-2,WP)*Lx/real(nx,WP)
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)
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
  mask(1,:) = 1
  !mask(nx,:) = 1
  
  return
end subroutine milk_crown_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine milk_crown_data
  use milk_crown
  use parser
  use math
  implicit none
  real(WP) :: dradius,bradius,thick,pheight,Udrop
  real(WP) :: dcenter,bcenter
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
  G = huge(1.0_WP)
  
  ! Set the initial distance field
  call parser_read('Drop radius',dradius,0.0_WP)
  if (dradius.gt.0.0_WP) call parser_read('Drop center',dcenter)
  call parser_read('Bubble radius',bradius,0.0_WP)
  if (bradius.gt.0.0_WP) call parser_read('Bubble center',bcenter)
  call parser_read('Initial G thickness',thick)
  call parser_read('Pool height',pheight,0.0_WP)
  call parser_read('Drop velocity',Udrop,0.0_WP)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! Compute distance function for the sphere
           if (dradius.gt.0.0_WP) G(i,j,k) = dradius - sqrt((xm(i)-dcenter)**2+ym(j)**2)
           ! Compute distance function for the sphere
           if (bradius.gt.0.0_WP) G(i,j,k) = - bradius + sqrt((xm(i)-bcenter)**2+ym(j)**2)
           ! Add water velocity
           if (G(i,j,k).gt.0.0_WP) U(i,j,k) = -Udrop
           ! Add the water pool
           if (pheight.gt.0.0_WP) G(i,j,k) = max(G(i,j,k),pheight-xm(i))
        end do
     end do
  end do
  
  ! Convert to tanh profile
  if (thick.gt.0.0_WP) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              G(i,j,k) = 0.5_WP*(tanh(0.5_WP*G(i,j,k)/thick)+1.0_WP)
           end do
        end do
     end do
  end if
  
  return
end subroutine milk_crown_data
