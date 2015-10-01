module drop_coll
  use string
  use precision
  use param
  implicit none

  ! Length of the domain
  real(WP) :: Lx,Lz,Ly

  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G

end module drop_coll

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine drop_coll_grid
  use drop_coll
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

  ! Cartesian and periodicity
  icyl = 0
  xper = 1
  yper = 1
  zper = 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create a uniform grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)-0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)-0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
  end do

  ! Create the mid points
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
end subroutine drop_coll_grid


subroutine drop_coll_data
  use drop_coll
  use parser
  use random
  use math
  implicit none

  integer :: i,j,k

  real(WP) :: d1,d2,dist1,dist2
  real(WP), dimension(3) :: cent1,cent2
  real(WP), dimension(3) :: U1,U2

  ! Allocate the array inflow
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))

  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  data = 0.0_WP

  ! Read parameters
  call parser_read('Drop 1 diameter',d1)
  call parser_read('Drop 1 center',cent1(1:3))
  call parser_read('Drop 1 velocity',U1(1:3))
  call parser_read('Drop 2 diameter',d2)
  call parser_read('Drop 2 center',cent2(1:3))
  call parser_read('Drop 2 velocity',U2(1:3))
  
  ! Set level set and velocity
  do k=1,nz
    do j=1,ny
       do i=1,nx
          dist1 = sqrt((xm(i)-cent1(1))**2 + (ym(j)-cent1(2))**2 + (zm(k)-cent1(3))**2 )
          if (dist1.le.0.5_WP*d1) then
             U(i,j,k) = U1(1)
             V(i,j,k) = U1(2)
             W(i,j,k) = U1(3)
          end if
          dist2 = sqrt((xm(i)-cent2(1))**2 + (ym(j)-cent2(2))**2 + (zm(k)-cent2(3))**2 )
          if (dist2.le.0.5_WP*d2) then
             U(i,j,k) = U2(1)
             V(i,j,k) = U2(2)
             W(i,j,k) = U2(3)
          end if
          G(i,j,k) = max(0.5_WP*d1-dist1,0.5_WP*d2-dist2)
       end do
    end do
  end do

  return
end subroutine drop_coll_data
