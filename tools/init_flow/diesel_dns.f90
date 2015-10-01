module diesel_dns
  use string
  use precision
  use param
  implicit none
  
  ! Domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: G
  
end module diesel_dns

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine diesel_dns_grid
  use diesel_dns
  use parser
  use math
  implicit none
  integer :: i,j,k
  
  ! Cartesian
  call parser_read('icyl',icyl)
  
  ! Periodicity
  call parser_read('xper',xper)
  if (icyl.eq.0) then
     yper = 1
  else
     yper = 0
  end if
  zper = 1
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('Lx',Lx)
  dx = Lx/real(nx,WP)
  call parser_read('ny',ny)
  call parser_read('Ly',Ly)
  if (icyl.eq.1) Ly = 0.5_WP*Ly
  dy = Ly/real(ny,WP)
  call parser_read('nz',nz)
  if (icyl.eq.1) then
     Lz = twoPi
  else
     call parser_read('Lz',Lz)
  end if
  dz = Lz/real(nz,WP)
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*dx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*dy-0.5_WP*Ly*(1.0_WP-real(icyl,WP))
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*dz-0.5_WP*Lz*(1.0_WP-real(icyl,WP))
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
end subroutine diesel_dns_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine diesel_dns_data
  use diesel_dns
  use parser
  use math
  use string
  implicit none
  
  real(WP) :: thick,rand,Djet,Ujet,d
  integer  :: i,j,k
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  G => data(:,:,:,4); names(4) = 'G'
  
  ! Level set
  call parser_read('Djet',Djet)
  G=0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if (icyl.eq.0) then
              d = sqrt(ym(j)**2+zm(k)**2+xm(i)**2)
           else
              d = sqrt(ym(j)**2+xm(i)**2)
           end if
           G(i,j,k) = 0.5_WP*Djet-d
        end do
     end do
  end do
  
  ! Velocity
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  !do j=1,ny
  !   do k=1,nz
  !      if (icyl.eq.0) then
  !         d = sqrt(ym(j)**2+zm(k)**2)
  !      else
  !         d = ym(j)
  !      end if
  !      if (d.le.0.5_WP*Djet) U(:,j,k) = Ujet
  !   end do
  !end do
  
  return
end subroutine diesel_dns_data
