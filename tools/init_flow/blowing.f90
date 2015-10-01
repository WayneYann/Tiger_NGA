module blowing
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module blowing

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine blowing_grid
  use blowing
  use parser
  implicit none
  
  integer  :: i,j,k
  real(WP) :: Lx,Ly,Lz
  real(WP) :: L,H

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  nz=1
  Lx=0.01_WP
  Ly=0.01_WP
  Lz=Lx/real(nx,WP)
  H=0.001_WP
  L=0.002_WP

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
     y(j) = real(j-1,WP)*Ly/real(ny,WP)-0.5_WP*Ly
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
        if ((abs(y(j)).ge.0.5_WP*H).and.(x(i+1).le.L)) then
           mask(i,j) = 1
        end if
     end do
  end do
  mask(:,1:2) = 1
  mask(:,ny) = 1
  mask(nx,:) = mask(nx-1,:)
  
  return
end subroutine blowing_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine blowing_data
  use blowing
  use parser
  implicit none

  integer :: i,j,k
  real(WP) :: xloc,liquid_lvl

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
  
  ! Place an interface
  call parser_read("Front location",xloc,0.001_WP)
  call parser_read("Liquid level",liquid_lvl,0.005_WP)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = xm(i)-xloc
           G(i,j,k) = min(liquid_lvl-xm(i),G(i,j,k))
        end do
     end do
  end do
  
  return
end subroutine blowing_data

