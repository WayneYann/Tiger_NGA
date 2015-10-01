module bubble
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  logical  :: two_dim
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module bubble

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine bubble_grid
  use bubble
  use parser
  implicit none
  
  integer :: i,j,k
  logical :: btm_wall
  
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
  if (nz.eq.1) then
     Lz=Lx/real(nx,WP)
     two_dim=.true.
  else
     two_dim=.false.
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
  
  ! Masks
  mask=0
  
  return
end subroutine bubble_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine bubble_data
  use bubble
  use parser
  use math
  implicit none
  real(WP) :: R1,R2,pheight
  real(WP), dimension(3) :: bc1,bc2
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
  call parser_read('Bubble 1 radius',R1,0.0_WP)
  if (R1.gt.0.0_WP) call parser_read('Bubble 1 center',bc1)
  call parser_read('Bubble 2 radius',R2,0.0_WP)
  if (R2.gt.0.0_WP) call parser_read('Bubble 2 center',bc2)
  call parser_read('Pool height',pheight,0.0_WP)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! Compute distance from the first bubble
           if (R1.gt.0.0_WP) then
              if (two_dim) then
                 G(i,j,k) = sqrt((xm(i)-bc1(1))**2+(ym(j)-bc1(2))**2)-R1
              else
                 G(i,j,k) = sqrt((xm(i)-bc1(1))**2+(ym(j)-bc1(2))**2+(zm(k)-bc1(3))**2)-R1
              end if
           end if
           ! Compute distance from the second bubble
           if (R2.gt.0.0_WP) then
              if (two_dim) then
                 G(i,j,k) = min(G(i,j,k),sqrt((xm(i)-bc2(1))**2+(ym(j)-bc2(2))**2)-R2)
              else
                 G(i,j,k) = min(G(i,j,k),sqrt((xm(i)-bc2(1))**2+(ym(j)-bc2(2))**2+(zm(k)-bc2(3))**2)-R2)
              end if
           end if
           ! Compute distance from pool
           if (pheight.gt.0.0_WP) then
              G(i,j,k) = min(G(i,j,k),pheight-ym(j))
           end if
        end do
     end do
  end do
  
  return
end subroutine bubble_data

