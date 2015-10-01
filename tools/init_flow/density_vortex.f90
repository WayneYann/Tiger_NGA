module density_vortex
  use string
  use precision
  use param
  implicit none
  
  ! Size of the domain
  real(WP) :: Lx,Ly,Lz
  real(WP) :: dx,dy,dz
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: rhoU
  real(WP), dimension(:,:,:), pointer :: rhoV
  real(WP), dimension(:,:,:), pointer :: rhoW
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  
  ! Density info
  real(WP) :: density_ratio

end module density_vortex

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine density_vortex_grid
  use density_vortex
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
  
  ! Set the periodicity
  xper = 1
  yper = 1
  zper = 1
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz
  do i=1,nx+1
     x(i) = (i-1)*dx
  end do
  do j=1,ny+1
     y(j) = (j-1)*dy
  end do
  do k=1,nz+1
     z(k) = (k-1)*dz
  end do
  
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
end subroutine density_vortex_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine density_vortex_data
  use density_vortex
  use parser
  implicit none
  
  integer :: i,j,k,dim
  real(WP) :: radius,Umax,r
  character(len=str_short) :: direction
  real(WP), dimension(3) :: Uc,center
  
 ! Allocate the array data
  nvar = 7
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U    => data(:,:,:,1); names(1) = 'U'
  V    => data(:,:,:,2); names(2) = 'V'
  W    => data(:,:,:,3); names(3) = 'W'
  P    => data(:,:,:,4); names(4) = 'P'
  ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  RHO  => data(:,:,:,6); names(6) = 'RHO'
  dRHO => data(:,:,:,7); names(7) = 'dRHO'

  ! General values
  dRHO = 0.0_WP
  P = 0.0_WP
  
  ! Set be convective velocity
  call parser_getsize('Convective velocity',dim)
  if (dim.ne.3) stop 'Convective velocity should be of size 3'
  call parser_read('Convective velocity',Uc)
  U = Uc(1)
  V = Uc(2)
  W = Uc(3)
  Umax = sqrt(sum(Uc**2))
  
  ! Get the center of the vortex
  call parser_getsize('Vortex center',dim)
  if (dim.ne.3) stop 'Vortex center should be of size 3'
  call parser_read('Vortex center',center)

  ! Create the vortex
  call parser_read('Vortex orientation',direction)
  call parser_read('Vortex radius',radius)

  ! Density stuff
  call parser_read('Density ratio',density_ratio)

  select case (trim(direction))
  case('x')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! V
              r = sqrt((y(j)-center(2))**2+(zm(k)-center(3))**2)
              if (r.lt.0.5_WP*radius) then
                 V(i,j,k) = V(i,j,k) - 2.0_WP*Umax/radius * (zm(k)-center(3))
              else if (r.lt.radius) then
                 V(i,j,k) = V(i,j,k) - 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (zm(k)-center(3))
              end if
              ! W
              r = sqrt((ym(j)-center(2))**2+(z(k)-center(3))**2)
              if (r.lt.0.5_WP*radius) then
                 W(i,j,k) = W(i,j,k) + 2.0_WP*Umax/radius * (ym(j)-center(2))
              else if (r.lt.radius) then
                 W(i,j,k) = W(i,j,k) + 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (ym(j)-center(2))
              end if
              ! rhoV
              r = sqrt((y(j)-center(2))**2+(zm(k)-center(3))**2)
              if (r.lt.radius) then
                 rhoV(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * V(i,j,k)
              else
                 rhoV(i,j,k) = V(i,j,k)
              end if
              ! rhoW
              r = sqrt((ym(j)-center(2))**2+(z(k)-center(3))**2)
              if (r.lt.radius) then
                 rhoW(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * W(i,j,k)
              else
                 rhoW(i,j,k) = W(i,j,k)
              end if
              ! ZMIX
              r = sqrt((ym(j)-center(2))**2+(zm(k)-center(3))**2)
              if (r.lt.radius) then
                 ZMIX(i,j,k) = 1.0_WP - (r/radius)**2
              else
                 ZMIX(i,j,k) = 0.0_WP
              end if
              ! RHO
              RHO(i,j,k)  = 1.0_WP + ZMIX(i,j,k)*(density_ratio-1.0_WP)
           end do
        end do
     end do
  case('y')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! W
              r = sqrt((z(k)-center(3))**2+(xm(i)-center(1))**2)
              if (r.lt.0.5_WP*radius) then
                 W(i,j,k) = W(i,j,k) - 2.0_WP*Umax/radius * (xm(i)-center(1))
              else if (r.lt.radius) then
                 W(i,j,k) = W(i,j,k) - 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (xm(i)-center(1))
              end if
              ! U
              r = sqrt((zm(k)-center(3))**2+(x(i)-center(1))**2)
              if (r.lt.0.5_WP*radius) then
                 U(i,j,k) = U(i,j,k) + 2.0_WP*Umax/radius * (zm(k)-center(3))
              else if (r.lt.radius) then
                 U(i,j,k) = U(i,j,k) + 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (zm(k)-center(3))
              end if
              ! rhoW
              r = sqrt((z(k)-center(3))**2+(xm(i)-center(1))**2)
              if (r.lt.radius) then
                 rhoW(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * W(i,j,k)
              else
                 rhoW(i,j,k) = W(i,j,k)
              end if
              ! rhoU
              r = sqrt((zm(k)-center(3))**2+(x(i)-center(1))**2)
              if (r.lt.radius) then
                 rhoU(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * U(i,j,k)
              else
                 rhoU(i,j,k) = U(i,j,k)
              end if
              ! ZMIX
              r = sqrt((zm(k)-center(3))**2+(xm(i)-center(1))**2)
              if (r.lt.radius) then
                 ZMIX(i,j,k) = 1.0_WP - (r/radius)**2
              else
                 ZMIX(i,j,k) = 0.0_WP
              end if
              ! RHO
              RHO(i,j,k)  = 1.0_WP + ZMIX(i,j,k)*(density_ratio-1.0_WP)
           end do
        end do
     end do
  case('z')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! U
              r = sqrt((x(i)-center(1))**2+(ym(j)-center(2))**2)
              if (r.lt.0.5_WP*radius) then
                 U(i,j,k) = U(i,j,k) - 2.0_WP*Umax/radius * (ym(j)-center(2))
              else if (r.lt.radius) then
                 U(i,j,k) = U(i,j,k) - 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (ym(j)-center(2))
              end if
              ! V
              r = sqrt((xm(i)-center(1))**2+(y(j)-center(2))**2)
              if (r.lt.0.5_WP*radius) then
                 V(i,j,k) = V(i,j,k) + 2.0_WP*Umax/radius * (xm(i)-center(1))
              else if (r.lt.radius) then
                 V(i,j,k) = V(i,j,k) + 2.0_WP*Umax*(1.0_WP/r-1.0_WP/radius) * (xm(i)-center(1))
              end if
              ! rhoU
              r = sqrt((x(i)-center(1))**2+(ym(j)-center(2))**2)
              if (r.lt.radius) then
                 rhoU(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * U(i,j,k)
              else
                 rhoU(i,j,k) = U(i,j,k)
              end if
              ! rhoV
              r = sqrt((xm(i)-center(1))**2+(y(j)-center(2))**2)
              if (r.lt.radius) then
                 rhoV(i,j,k) = (density_ratio-(density_ratio-1.0_WP)*r**2/radius**2) * V(i,j,k)
              else
                 rhoV(i,j,k) = V(i,j,k)
              end if
              ! ZMIX
              r = sqrt((xm(i)-center(1))**2+(ym(j)-center(2))**2)
              if (r.lt.radius) then
                 ZMIX(i,j,k) = 1.0_WP - (r/radius)**2
              else
                 ZMIX(i,j,k) = 0.0_WP
              end if
              ! RHO
              RHO(i,j,k)  = 1.0_WP + ZMIX(i,j,k)*(density_ratio-1.0_WP)
           end do
        end do
     end do
  end select

  ! Density

  return
end subroutine density_vortex_data


! ========================= !
! Create the chemtable data !
! ========================= !
subroutine density_vortex_chemtable
  use density_vortex
  implicit none

  ! Setup the coordinates
  n1 = 2
  n2 = 2
  n3 = 2

  ! Allocate arrays
  allocate(x1(n1))
  allocate(x2(n2))
  allocate(x3(n3))

  ! Setup the arrays
  x1(1) = 0.0_WP
  x1(2) = 1.0_WP
  x2 = 0.0_WP
  x3 = 0.0_WP

  ! Chemtable model
  combModel = 'Steady Flamelet'

  ! Table itself
  nvar_chem = 4
  allocate(table(n1,n2,n3,nvar_chem))
  allocate(names_chem(nvar_chem))
  allocate(chem_mask(n1,n2,n3))

  ! Masked values in the chemtable
  chem_mask = 0

  ! Density
  names_chem(1)  = 'RHO'
  table(1,:,:,1) = 1.0_WP
  table(2,:,:,1) = density_ratio
  ! Temperature
  names_chem(2)  = 'T'
  table(1,:,:,2) = 1.0_WP
  table(2,:,:,2) = 1.0_WP/density_ratio
  ! Viscosity
  names_chem(3)  = 'VISC'
  table(1,:,:,3) = 0.0_WP
  table(2,:,:,3) = 0.0_WP
  ! Diffusivity
  names_chem(4)  = 'DIFF'
  table(1,:,:,4) = 0.0_WP
  table(2,:,:,4) = 0.0_WP

  return
end subroutine density_vortex_chemtable
