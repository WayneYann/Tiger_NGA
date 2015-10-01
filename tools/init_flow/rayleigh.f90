module rayleigh
  use string
  use precision
  use param
  implicit none
  
  ! Problem description
  real(WP) :: Lx,Ly,Lz
  real(WP) :: lambda,r0
  real(WP) :: dx,dy,dz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module rayleigh

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine rayleigh_grid
  use rayleigh
  use parser
  use math
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('Ly',Ly)
  call parser_read('ny',ny)
  call parser_read('nx',nx)
  dy = Ly/real(ny,WP)
  dx = dy
  Lx = real(nx,WP)*dx
  lambda = Lx
  r0 = Ly/3.0_WP
  dz = dx
  nz = 1
  Lz = twoPi/real(2048,WP)
  xper = 1
  yper = 0
  zper = 1
  
  ! Print out parameters
  print*,'Lx =',Lx
  print*,'nx =',nx
  print*,'dx =',dx
  print*,'Ly =',Ly
  print*,'ny =',ny
  print*,'dy =',dy
  print*,'l_ =',lambda
  print*,'k_o=',twoPi*r0/lambda
  
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
  
  return
end subroutine rayleigh_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine rayleigh_data
  use rayleigh
  use parser
  use math
  implicit none
  real(WP) :: amp,U0,thick
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
  call parser_read('Amplitude',amp,1.0e-4_WP*dx)
  print*,'eta_0 =',amp
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = r0-ym(j)+amp*cos(twoPi*xm(i)/lambda)
        end do
     end do
  end do
  
  ! Setup th initial velocity
  call parser_read('Cylinder velocity',U0,0.0_WP)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if (G(i,j,k).ge.0.0_WP) U(i,j,k)=U0
        end do
     end do
  end do
  
  ! Convert to tanh profile
  call parser_read('Initial G thickness',thick,0.75_WP*0.5_WP*dx)
  print*,'eps =',thick
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
end subroutine rayleigh_data
