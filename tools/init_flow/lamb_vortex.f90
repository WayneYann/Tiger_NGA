module lamb_vortex
  use precision
  use math
  use param
  implicit none

  ! Length and diameter of the domain
  real(WP) :: L,radius
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  
end module lamb_vortex

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine lamb_vortex_grid
  use lamb_vortex
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)
  
  call parser_read('Full length',L)
  call parser_read('Full radius',radius)
  
  ! Set the periodicity
  xper = 1
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = (i-1)*L/nx
  end do
  do j=1,ny+1
     y(j) = (j-1)*radius/ny
  end do
  do k=1,nz+1
     z(k) = (k-1)*twoPi/nz
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

  return
end subroutine lamb_vortex_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine lamb_vortex_data
  use lamb_vortex
  use parser
  implicit none

  integer :: j,k
  real(WP) :: a,C,vel
  real(WP) :: r,cos_theta,sin_theta
  real(WP) :: x0,y0,ux,uy,ur,ut
  
  ! Allocate the array data
  nvar = 4
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  
  ! Lamb Dipole function
  call parser_read('Dipole radius',a)
  call parser_read('Dipole velocity',vel)
  call parser_read('Dipole loc. x',x0)
  call parser_read('Dipole loc. y',y0)
  C = 2.0_WP/bessj0(bessj1_zero)

  do j=1,ny
     do k=1,nz
        ! Radial velocity
        r = sqrt((x0-y(j)*cos(zm(k)))**2+(y0-y(j)*sin(zm(k)))**2)
        !r = r + epsilon(r)
        if (r.eq.0.0_WP) then
           cos_theta = cos(zm(k))
           sin_theta = sin(zm(k))
           r = r+epsilon(r)
        else
           cos_theta = (y(j)*cos(zm(k))-x0)/r
           sin_theta = (y(j)*sin(zm(k))-y0)/r
        end if
        if (r.le.a) then
           r = bessj1_zero*r/a
           ur = vel * ( C*bessj1(r)/r - 1.0_WP ) * cos_theta
           ut = vel * ( 1.0_WP - C*(bessj0(r)-bessj1(r)/r) ) * sin_theta
        else
           r = r/a
           ur = - vel / r**2 * cos_theta
           ut = - vel / r**2 * sin_theta
        end if
        ux = cos_theta*ur - sin_theta*ut
        uy = sin_theta*ur + cos_theta*ut
        V(:,j,k) = cos(zm(k))*ux + sin(zm(k))*uy
        
        ! Orthoradial velocity
        r = sqrt((x0-ym(j)*cos(z(k)))**2+(y0-ym(j)*sin(z(k)))**2)
        !r = r + epsilon(r)
        cos_theta = (ym(j)*cos(z(k))-x0)/r
        sin_theta = (ym(j)*sin(z(k))-y0)/r
        if (r.le.a) then
           r = bessj1_zero*r/a
           ur = vel * ( C*bessj1(r)/r - 1.0_WP ) * cos_theta
           ut = vel * ( 1.0_WP - C*(bessj0(r)-bessj1(r)/r) ) * sin_theta
        else
           r = r/a
           ur = - vel / r**2 * cos_theta
           ut = - vel / r**2 * sin_theta
        end if
        ux = cos_theta*ur - sin_theta*ut
        uy = sin_theta*ur + cos_theta*ut
        W(:,j,k) = -sin(z(k))*ux + cos(z(k))*uy

        ! Pressure
        r = sqrt((x0-ym(j)*cos(zm(k)))**2+(y0-ym(j)*sin(zm(k)))**2)
        !r = r + epsilon(r)
        cos_theta = (ym(j)*cos(zm(k))-x0)/r
        sin_theta = (ym(j)*sin(zm(k))-y0)/r
        if (r.le.a) then
           r = bessj1_zero*r/a
           P(:,j,k) = - 0.5_WP * vel**2 * C**2 * ( &
                ((bessj0(r)-bessj1(r)/r)*sin_theta)**2 + &
                (bessj1(r)/r*cos_theta) **2            + &
                (bessj1(r)*sin_theta)**2 )
        else
           r = r/a
           P(:,j,k) = - 0.5_WP * vel**2 * ( &
                1.0_WP + 1.0_WP/r**4 - 2.0_WP/r**2*(cos_theta**2-sin_theta**2))
        end if
     end do
  end do
  
  U = 0.0_WP

  return
end subroutine lamb_vortex_data

