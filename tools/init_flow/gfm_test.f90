module gfm_test
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  character(len=str_medium) :: dim
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: LVLSET
  
end module gfm_test

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine gfm_test_grid
  use gfm_test
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
  call parser_read('xper',xper)
  call parser_read('yper',yper)
  call parser_read('zper',zper)
  
  ! Find dimensionality
  dim='xyz'
  if      (nx.le.2) then
     dim='yz'
  else if (ny.le.2) then
     dim='xz'
  else if (nz.le.2) then
     dim='xy'
  end if
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)-0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)-0.5_WP*Ly
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
  
  return
end subroutine gfm_test_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine gfm_test_data
  use gfm_test
  use parser
  use math
  implicit none
  real(WP) :: dradius,bradius,thick,pheight,pwidth,dist,amp,mode,angle
  real(WP), dimension(3) :: dcenter,bcenter
  integer :: i,j,k
  
  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  LVLSET => data(:,:,:,5); names(5) = 'LVLSET'
  
  ! Initialize
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  LVLSET = huge(1.0_WP)
  
  ! Set the initial distance field
  call parser_read('Drop radius',dradius,0.0_WP)
  if (dradius.gt.0.0_WP) call parser_read('Drop center',dcenter)
  call parser_read('Bubble radius',bradius,0.0_WP)
  if (bradius.gt.0.0_WP) call parser_read('Bubble center',bcenter)
  call parser_read('Amplitude',amp,0.0_WP)
  call parser_read('Mode',mode,0.0_WP)
  !call parser_read('Initial G thickness',thick)
  call parser_read('Pool height',pheight,0.0_WP)
  call parser_read('Pool width',pwidth,0.0_WP)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! Compute distance function for the sphere
           if (dradius.gt.0.0_WP) then
              select case (trim(dim))
              case ('xyz')
                 LVLSET(i,j,k) = dradius - sqrt((xm(i)-dcenter(1))**2+(ym(j)-dcenter(2))**2+(zm(k)-dcenter(3))**2)
              case ('yz')
                 angle = arctan(ym(j),zm(k))
                 LVLSET(i,j,k) = dradius-amp*cos(mode*angle) - sqrt((ym(j)-dcenter(2))**2+(zm(k)-dcenter(3))**2)
              case ('xz')
                 LVLSET(i,j,k) = dradius - sqrt((xm(1)-dcenter(1))**2+(zm(k)-dcenter(3))**2)
              case ('xy')
                 angle = arctan(xm(i),ym(j))
                 LVLSET(i,j,k) = dradius-amp*cos(mode*angle) - sqrt((xm(i)-dcenter(1))**2+(ym(j)-dcenter(2))**2)
              end select
           end if
           ! Compute distance function for the sphere
           if (bradius.gt.0.0_WP) then
              select case (trim(dim))
              case ('xyz')
                 LVLSET(i,j,k) = - bradius + sqrt((xm(i)-bcenter(1))**2+(ym(j)-bcenter(2))**2+(zm(k)-bcenter(3))**2)
              case ('yz')
                 LVLSET(i,j,k) = - bradius + sqrt((ym(j)-bcenter(2))**2+(zm(k)-bcenter(3))**2)
              case ('xz')
                 LVLSET(i,j,k) = - bradius + sqrt((xm(1)-bcenter(1))**2+(zm(k)-bcenter(3))**2)
              case ('xy')
                 LVLSET(i,j,k) = - bradius + sqrt((xm(i)-bcenter(1))**2+(ym(j)-bcenter(2))**2)
              end select
           end if
           ! Add the water pool
           if (pheight.gt.0.0_WP .and. pwidth .gt.0.0_WP) then
              ! Compute the distance from the pool
              if     (ym(j).le.-0.5_WP*Ly+pheight .and. zm(k).le.0.5_WP*pwidth) then
                 dist = min(-ym(j)+(-0.5_WP*Ly+pheight),-abs(zm(k))+0.5_WP*pwidth)
              elseif (ym(j).le.-0.5_WP*Ly+pheight .and. abs(zm(k)).gt.0.5_WP*pwidth) then
                 dist = -abs(zm(k))+0.5_WP*pwidth
              elseif (ym(j).gt.-0.5_WP*Ly+pheight .and. abs(zm(k)).le.0.5_WP*pwidth) then
                 dist = -ym(j)+(-0.5_WP*Ly+pheight)
              else
                 dist = -sqrt((-ym(j)+(-0.5_WP*Ly+pheight))**2+(-abs(zm(k))+0.5_WP*pwidth)**2)
              end if
              ! If closer, replace G by dist
              if (abs(LVLSET(i,j,k)).gt.abs(dist)) then
                 LVLSET(i,j,k) = dist
              end if
           end if
        end do
     end do
  end do
  
  ! Convert to tanh profile
  !if (thick.gt.0.0_WP) then
  !   do k=1,nz
  !      do j=1,ny
  !         do i=1,nx
  !            G(i,j,k) = 0.5_WP*(tanh(0.5_WP*G(i,j,k)/thick)+1.0_WP)
  !         end do
  !      end do
  !   end do
  !end if
  
  return
end subroutine gfm_test_data
