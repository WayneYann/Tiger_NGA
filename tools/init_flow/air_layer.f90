module air_layer
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
 
  ! Dimension
  real(WP) :: Lx,Ly,Lz
  real(WP) :: aH,lH,lL
 
end module air_layer

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine air_layer_grid
  use air_layer
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
  call parser_read('Air height',aH)
  call parser_read('Lip height',lH)
  call parser_read('Lip length',lL)
  
  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)-Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
  end do

  ! Create the mid points 
  do i=1,nx
     xm(i)= 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)= 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k)= 0.5_WP*(z(k)+z(k+1))
  end do
  
  ! Create the masks
  mask = 0
  do i=1,nx-1
     do j=1,ny
        if ((abs(y(j)).ge.aH).and.(abs(y(j)).le.aH+lH).and.(x(i+1).le.lL)) then
           mask(i,j) = 1
        end if
     end do
  end do
  mask(:, 1) = 0
  mask(:,ny) = 1
  mask(nx,:) = mask(nx-1,:)
  
  return
end subroutine air_layer_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine air_layer_data
  use air_layer
  use parser
  implicit none

  integer :: i,j,k
  real(WP) :: Ugas,Uliq,thick

  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  G = 0.0_WP

  ! Initialize the flow field
  call parser_read('Initial G thickness',thick,-1.0_WP)
  call parser_read('Initial Uliq',Uliq,0.0_WP)
  call parser_read('Initial Ugas',Ugas,0.0_WP)  
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k)=abs(ym(j))-(aH+lH)
           if (G(i,j,k).gt.0.0_WP) then
              U(i,j,k)=Uliq
           else
              U(i,j,k)=Ugas
           end if
        end do
     end do
  end do

  ! Transform to tanh
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
end subroutine air_layer_data

