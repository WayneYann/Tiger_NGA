module rt_vortex
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  integer  :: nyc
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module rt_vortex

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine rt_vortex_grid
  use rt_vortex
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  nx = 1
  call parser_read('ny core',nyc)
  call parser_read('ny total',ny)
  call parser_read('nz',nz)
  call parser_read('Ly core',Ly)
  Lx = Ly/nyc
  Lz = twoPi
  xper = 1
  yper = 0
  zper = 1
  
  ! Cylindrical
  icyl = 1
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)-0.5_WP*Lx
  end do
  do j=1,nyc+1
     y(j) = real(j-1,WP)*Ly/real(nyc,WP)
  end do
  do j=nyc+2,ny+1
     y(j) = y(j-1)+(y(j-1)-y(j-2))*1.2_WP
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
  
  return
end subroutine rt_vortex_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine rt_vortex_data
  use rt_vortex
  use parser
  use math
  implicit none
  real(WP) :: drad,amp,mode,thick
  real(WP) :: vrad,vcirc
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
  
  ! Setup the Oseen vortex
  call parser_read('Vortex radius',vrad)
  call parser_read('Vortex circulation',vcirc)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           W(i,j,k) = vcirc*(1.0_WP-exp(-ym(j)**2/vrad**2))/(twoPi*ym(j))
        end do
     end do
  end do
  
  ! Setup the initial distance field
  call parser_read('Drop radius',drad)
  call parser_read('Amplitude',amp)
  call parser_read('Mode',mode)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = drad-ym(j)-amp*cos(mode*zm(k))
        end do
     end do
  end do
  
  ! Convert to tanh profile
  call parser_read('Initial G thickness',thick)
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
end subroutine rt_vortex_data
