module vortex_ring
  use string
  use precision
  use param
  implicit none
  
  ! Size of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  logical  :: bwall
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module vortex_ring


! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine vortex_ring_grid
  use vortex_ring
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Cylindrical
  call parser_read('icyl',icyl)
  
  ! Bottom wall
  call parser_read('Bottom wall',bwall)
  
  ! Periodicity
  if (bwall) then
     xper = 0
  else
     xper = 1
  end if
  if (icyl.eq.1) then
     yper = 0
  else
     yper = 1
  end if
  zper = 1
  
  ! Domain size
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  if (icyl.eq.1) then
     if (nz.eq.1) then
        Lz=twoPi/real(max(10*nx,10*ny),WP)
     else
        Lz=twoPi
     end if
  else
     call parser_read('Lz',Lz)
  end if
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  if (bwall) then
     do i=1,nx+1
        x(i) = real(i-2,WP)*Lx/real(nx-1,WP)
     end do
  else
     do i=1,nx+1
        x(i) = real(i-1,WP)*Lx/real(nx,WP)
     end do
  end if
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)+real(icyl-1,WP)*0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)+real(icyl-1,WP)*0.5_WP*Lz
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
  if (bwall) mask(1,:) = 1
  
  return
end subroutine vortex_ring_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine vortex_ring_data
  use vortex_ring
  use parser
  use math
  implicit none
  
  integer  :: i,j,k
  real(WP) :: vheight,vint,vradius
  real(WP) :: s
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Initialize
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  
  ! Generate the vortex
  call parser_read('Vortex height',vheight)
  call parser_read('Vortex radius',vradius)
  call parser_read('Vortex intensity',vint)
  if (icyl.eq.1) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! Main vortex
              s = sqrt((x(i)-vheight)**2+(ym(j)-vradius)**2)+tiny(1.0_WP)
              U(i,j,k) = +1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(ym(j)-vradius)/s
              s = sqrt((xm(i)-vheight)**2+(y(j)-vradius)**2)+tiny(1.0_WP)
              V(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(xm(i)-vheight)/s
              ! Secondary vortex to ensure zero vorticity at the axis
              s = sqrt((x(i)-vheight)**2+(ym(j)+vradius)**2)+tiny(1.0_WP)
              U(i,j,k) = U(i,j,k) -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(ym(j)+vradius)/s
              s = sqrt((xm(i)-vheight)**2+(y(j)+vradius)**2)+tiny(1.0_WP)
              V(i,j,k) = V(i,j,k) +1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(xm(i)-vheight)/s
           end do
        end do
     end do
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! Axial
              s = sqrt((x(i)-vheight)**2+(sqrt(ym(j)**2+zm(k)**2)-vradius)**2)+tiny(1.0_WP)
              U(i,j,k) = +1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(sqrt(ym(j)**2+zm(k)**2)-vradius)/s
              ! V
              s = sqrt((xm(i)-vheight)**2+(sqrt(y(j)**2+zm(k)**2)-vradius)**2)+tiny(1.0_WP)
              V(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(xm(i)-vheight)/s*y(j)/sqrt(y(j)**2+zm(k)**2)
              ! W
              s = sqrt((xm(i)-vheight)**2+(sqrt(ym(j)**2+z(k)**2)-vradius)**2)+tiny(1.0_WP)
              W(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(xm(i)-vheight)/s*z(k)/sqrt(ym(j)**2+z(k)**2)
           end do
        end do
     end do
  end if
  
  return
end subroutine vortex_ring_data
