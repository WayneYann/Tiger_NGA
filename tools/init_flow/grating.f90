module grating
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module grating

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine grating_grid
  use grating
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  nz=1
  call parser_read('Wt',Lx)
  call parser_read('Ht',Ly)
  Lz=Lx/real(nx,WP)
  xper=1
  yper=0
  zper=1
  
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
     y(j) = real(j-2,WP)*Ly/real(ny,WP)
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
  mask(:,1) = 1
  
  return
end subroutine grating_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine grating_data
  use grating
  use parser
  use math
  implicit none
  real(WP) :: Hg,Wg,Hr
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
  call parser_read('Hg',Hg)
  call parser_read('Hr',Hr)
  call parser_read('Wg',Wg)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if     (ym(j).le.Hg+Hr .and. ym(j).gt.Hr .and. xm(i).le.0.5_WP*Wg) then
              G(i,j,k) = min(-ym(j)+Hg+Hr,-abs(xm(i))+0.5_WP*Wg)
           elseif (ym(j).le.Hg+Hr .and. ym(j).gt.Hr .and. abs(xm(i)).gt.0.5_WP*Wg) then
              G(i,j,k) = min(-ym(j)+Hr   ,-abs(xm(i))+0.5_WP*Wg)
           elseif (ym(j).gt.Hg+Hr .and. abs(xm(i)).le.0.5_WP*Wg) then
              G(i,j,k) = -ym(j)+Hg+Hr
           elseif (ym(j).gt.Hg+Hr .and. abs(xm(i)).gt.0.5_WP*Wg) then
              G(i,j,k) = -sqrt((-ym(j)+Hg+Hr)**2+(-abs(xm(i))+0.5_WP*Wg)**2)
           elseif (ym(j).le. Hr .and. xm(i).le.0.5_WP*Wg) then
              G(i,j,k) = sqrt((-ym(j)+Hr)**2+(-abs(xm(i))+0.5_WP*Wg)**2)
           elseif (ym(j).le. Hr .and. xm(i).gt.0.5_WP*Wg) then
              G(i,j,k) = -ym(j)+Hr
           end if
        end do
     end do
  end do
  
  return
end subroutine grating_data
