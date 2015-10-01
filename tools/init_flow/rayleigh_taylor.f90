module rayleigh_taylor
  use string
  use precision
  use param
  implicit none
  
  ! Length of the domain
  real(WP) :: Lx,Ly,Lz
  
  ! Scalar model
  character(len=str_medium) :: scalar_model
  real(WP) :: density_ratio

  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: G
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  
end module rayleigh_taylor

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine rayleigh_taylor_grid
  use rayleigh_taylor
  use parser
  implicit none
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  nz = 1
  Lz = Lx/real(nx,WP)
  
  ! Set the periodicity
  xper = 1
  yper = 0
  zper = 1
  
  ! Cartesian
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1) ,y(ny+1) ,z(nz+1))
  allocate(xm(nx+1),ym(ny+1),zm(nz+1))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP) - 0.5_WP*Lx
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP) - 0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP) - 0.5_WP*Lz
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
  mask(:,1)  = 1
  mask(:,ny) = 1
  
  return
end subroutine rayleigh_taylor_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine rayleigh_taylor_data
  use rayleigh_taylor
  use parser
  use math
  use random
  implicit none
  integer :: i,j,k,nmodes
  real(WP) :: thick,amp
  real(WP), dimension(:), pointer :: modes
  
  ! Type of density model
  call parser_read('Scalar model',scalar_model)
  
  select case(trim(scalar_model))
  case ('levelset')
      ! Allocate the array data
     nvar = 4
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U => data(:,:,:,1); names(1) = 'U'
     V => data(:,:,:,2); names(2) = 'V'
     W => data(:,:,:,3); names(3) = 'W'
     G => data(:,:,:,4); names(4) = 'G'
  case ('mixture fraction')
     ! Allocate the array data
     nvar = 7
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U => data(:,:,:,1); names(1) = 'U'
     V => data(:,:,:,2); names(2) = 'V'
     W => data(:,:,:,3); names(3) = 'W'
     P => data(:,:,:,4); names(4) = 'P'
     ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
     RHO  => data(:,:,:,6); names(6) = 'RHO'
     dRHO => data(:,:,:,7); names(7) = 'dRHO'
     G => ZMIX
  case default
     stop "rayleigh_taylor_data: unknown Scalar model"
  end select
  
  ! Read input
  call parser_read('Initial amplitude',amp)

  ! Set velocity
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  
  ! Get modes
  call parser_getsize('Modes',nmodes)
  allocate(modes(nmodes))
  call parser_read('Modes',modes)
  
  ! Set G
  do k=1,nz
     do j=1,ny
        do i=1,nx
           G(i,j,k) = ym(j) + amp*sum(cos(modes(:)*Pi*xm(i)/Lx))
        end do
     end do
  end do
  
  ! Other variables
  select case(trim(scalar_model))
  case ('levelset')
     
     ! All done
     
  case ('mixture fraction')
     
     ! Convert to tanh profile
     call parser_read('Initial thickness',thick)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              G(i,j,k) = 0.5_WP*(tanh(-0.5_WP*G(i,j,k)/thick)+1.0_WP)
           end do
        end do
     end do
     
     ! Density
     call parser_read('Density ratio',density_ratio)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              RHO(i,j,k) = 1.0_WP/(1.0_WP+(density_ratio-1.0_WP)*G(i,j,k))
           end do
        end do
     end do
     dRHO = 0.0_WP
     P    = 0.0_WP
     
  case default
     stop "rayleigh_taylor_data: unknown Scalar model"
  end select
  
  return
end subroutine rayleigh_taylor_data




! ========================= !
! Create the chemtable data !
! ========================= !
subroutine rayleigh_taylor_chemtable
  use rayleigh_taylor
  use parser
  implicit none

  integer  :: i
  real(WP) :: visc,diff
  
  ! Return if chemtable not needed
  if (trim(scalar_model).ne.'mixture fraction') return
  
  ! Setup the coordinates
  n1 = 100
  n2 = 2
  n3 = 2

  ! Allocate arrays
  allocate(x1(n1))
  allocate(x2(n2))
  allocate(x3(n3))

  ! Setup the arrays
  do i=1,n1
     x1(i) = real(i-1,WP)/real(n1-1,WP)
  end do
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
  do i=1,n1
     table(i,:,:,1) = 1.0_WP/(1.0_WP+(density_ratio-1.0_WP)*x1(i))
  end do
  ! Temperature
  names_chem(2)  = 'T'
  table(:,:,:,2) = 1.0_WP/table(:,:,:,1)
  ! Viscosity
  call parser_read('Kinematic viscosity',visc)
  names_chem(3)  = 'VISC'
  table(:,:,:,3) = visc*table(:,:,:,1)
  ! Diffusivity
  call parser_read('Kinematic diffusivity',diff)
  names_chem(4)  = 'DIFF'
  table(:,:,:,4) = diff*table(:,:,:,1)

  return
end subroutine rayleigh_taylor_chemtable
