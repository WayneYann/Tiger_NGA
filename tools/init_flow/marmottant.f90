module marmottant
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP), parameter :: D =  7.6e-3_WP
  real(WP), parameter :: R1 = 3.8e-3_WP
  real(WP), parameter :: R2 = 4.0e-3_WP
  real(WP), parameter :: R3 = 5.7e-3_WP
  real(WP) :: Lx,Ly,Lz
  real(WP) :: Lp
  real(WP) :: dx
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: G
  
end module marmottant

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine marmottant_grid
  use marmottant
  use parser
  use math
  implicit none
  integer :: i,j,k
  
  ! Read in the size of the domain
  !call parser_read('nx',nx)
  !call parser_read('ny',ny)
  !call parser_read('nz',nz)
  call parser_read('dx',dx)
  call parser_read('Lp',Lp)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Cylindrical coordinates',icyl)
  if (icyl.eq.0) then
     call parser_read('Lz',Lz)
  else
     call parser_read('nz',nz)
     if (nz.gt.1) then
        Lz = twoPi
     else
        Lz = 2.0_WP*dx/D
        print*,'Lz for same equivalent dx:',Lz
     end if
  end if
  nx = int((Lp+Lx)*D/dx)
  ny = int(    Ly *D/dx)
  if (icyl.eq.0) then
     nz = int( Lz*D/dx)
  else
     !nz = int(0.5_WP*Lz*D/dx)
     print*,'nz for same equivalent dx:',int(Lz*D/dx)!,int(0.5_WP*Lz*D/dx)
  end if
  print*,nx,ny,nz
  
  ! Periodicity
  xper = 0
  yper = 0
  zper = 1
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx-Lp*D
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly*D/real(ny,WP)
  end do
  if (icyl.eq.0) then
     do k=1,nz+1
        z(k) = real(k-1,WP)*Lz*D/real(nz,WP)
     end do
  else
     do k=1,nz+1
        z(k) = real(k-1,WP)*Lz/real(nz,WP)
     end do
  end if
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
  do i=1,nx
     do j=1,ny
        if (x(i).lt.0.0_WP) then
           if (abs(y(j)).ge.R1 .and. abs(y(j)).lt.R2) mask(i,j) = 1
           if (abs(y(j)).ge.R3                      ) mask(i,j) = 1
        end if
     end do
  end do
  
  return
end subroutine marmottant_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine marmottant_data
  use marmottant
  use parser
  implicit none
  
  real(WP) :: thick,glocationx,glocationy,Ul,Ug,Wg
  integer :: i,j,k
  logical :: use_G
  
  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  W => data(:,:,:,4); names(4) = 'W'
  G => data(:,:,:,5); names(5) = 'LVLSET'
  
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  G = 0.0_WP
  
  ! Initialize the velocity
  call parser_read('U liquid',Ul,0.0_WP)
  call parser_read('U gas',Ug,0.0_WP)
  call parser_read('W gas',Wg,0.0_WP)
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if (abs(ym(j)).lt.R1                       ) U(i,j,k) = Ul
           if (abs(ym(j)).ge.R2 .and. abs(ym(j)).le.R3) U(i,j,k) = Ug
           if (abs(ym(j)).ge.R2 .and. abs(ym(j)).le.R3) W(i,j,k) = Wg
        end do
     end do
  end do
  
  ! Initialize G as distance
  call parser_is_defined('Initial G x',use_G)
  if (.not.use_G) return
  call parser_read('Initial G x',glocationx)
  glocationy = 0.5_WP*(R1+R2)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if      ( xm(i).le.glocationx .and. abs(ym(j)).le.glocationy ) then
              G(i,j,k) = min(-xm(i)+glocationx,-abs(ym(j))+glocationy)
              if (xm(i).lt.0.0_WP) G(i,j,k) = min(sqrt((-xm(i))**2+(-abs(ym(j))+glocationy)**2),-xm(i)+glocationx)
           else if ( xm(i).le.glocationx .and. abs(ym(j)).gt.glocationy ) then
              G(i,j,k) = -abs(ym(j))+glocationy
              if (xm(i).lt.0.0_WP) G(i,j,k) = min(-sqrt((-xm(i))**2+(-abs(ym(j))+glocationy)**2),-abs(ym(j))+glocationy)
           else if ( xm(i).gt.glocationx .and. abs(ym(j)).le.glocationy ) then
              G(i,j,k) = -xm(i)+glocationx
           else
              G(i,j,k) = -sqrt((-xm(i)+glocationx)**2+(-abs(ym(j))+glocationy)**2)
           end if
        end do
     end do
  end do
  
  ! Convert to tanh profile
  !call parser_read('Initial G thickness',thick)
  !do k=1,nz
  !   do j=1,ny
  !      do i=1,nx
  !         G(i,j,k) = 0.5_WP*(tanh(0.5_WP*G(i,j,k)/thick)+1.0_WP)
  !      end do
  !   end do
  !end do
  
  return
end subroutine marmottant_data
