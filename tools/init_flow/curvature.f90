module curvature
  use string
  use precision
  use param
  implicit none
  
  ! Domain size
  real(WP) :: Lx,Ly,Lz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module curvature

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine curvature_grid
  use curvature
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  xper = 1
  yper = 1
  zper = 1
  Lx = 2.0_WP
  Ly = Lx/real(nx,WP)*real(ny,WP)
  Lz = Lx/real(nx,WP)*real(nz,WP)
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = Lx*real(i-1,WP)/real(nx,WP)-0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j) = Ly*real(j-1,WP)/real(ny,WP)-0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = Lz*real(k-1,WP)/real(nz,WP)-0.5_WP*Lz
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
end subroutine curvature_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine curvature_data
  use curvature
  use parser
  use math
  implicit none
  real(WP) :: thick,amp
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
  
  ! Setup the interface
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = 0.5_WP-sqrt(xm(i)**2+ym(j)**2+zm(k)**2)
        end do
     end do
  end do
  
  ! Convert to tanh profile
  call parser_read('Initial G thickness',thick,0.75_WP*0.5_WP*Lx/real(nx,WP))
  print*,'eps=',thick
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
end subroutine curvature_data
