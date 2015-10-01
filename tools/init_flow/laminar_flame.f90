module laminar_flame
  use precision
  use param
  implicit none

  ! Length and diameter of the domain
  real(WP) :: L,height,width
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: Yf
  real(WP), dimension(:,:,:), pointer :: Yo
  real(WP), dimension(:,:,:), pointer :: Yp
  real(WP), dimension(:,:,:), pointer :: T
  
  ! Stoichiometric coefficients
  real(WP) :: nu_f
  real(WP) :: nu_o
  real(WP) :: nu_p
  
  ! Molecular weight
  real(WP) :: W_f   ! [kg/mol]
  real(WP) :: W_o   ! [kg/mol]
  real(WP) :: W_p   ! [kg/mol]
  real(WP) :: W_d   ! [kg/mol]
  
  ! Composition
  real(WP) :: Yf_1,Yf_2,Yo_1,Yo_2
  real(WP) :: Yf_u, Yo_u, Yd_u
  real(WP) :: nu
  real(WP) :: phi,ZMIX

end module laminar_flame

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine laminar_flame_grid
  use laminar_flame
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  ny = 1
  nz = 1
  
  call parser_read('Length',L)
  height = L/real(nx,WP)
  width  = L/real(nx,WP)
  
  ! Set the periodicity
  xper = 0
  yper = 1
  zper = 1

  ! Cartesian
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = (i-1)*L/nx
  end do
  do j=1,ny+1
     y(j) = (j-1)*height/ny
  end do
  do k=1,nz+1
     z(k) = (k-1)*width/nz
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
  mask(1,:) = 1

  return
end subroutine laminar_flame_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine laminar_flame_data
  use laminar_flame
  use parser
  implicit none
  real(WP) :: dx,xloc
  integer  :: i
  logical  :: isdef
  
  ! Read the stoichiometric coefficients
  call parser_read('nu_f',nu_f)
  call parser_read('nu_o',nu_o)
  call parser_read('nu_p',nu_p)
  nu_f = -abs(nu_f)
  nu_o = -abs(nu_o)
  nu_p = +abs(nu_p)
  
  ! Read the molecular masses
  call parser_read('W_f',W_f)
  call parser_read('W_o',W_o)
  call parser_read('W_d',W_d)
  W_p = (abs(nu_f)*W_f + abs(nu_o)*W_o)/nu_p
  
  ! Get the composition of both streams
  call parser_read('Yf_1',Yf_1,1.0_WP)
  call parser_read('Yf_2',Yf_2,0.0_WP)
  call parser_read('Yo_1',Yo_1,0.0_WP)
  call parser_read('Yo_2',Yo_2,0.233_WP)
  nu  = nu_o*W_o/(nu_f*W_f)
  
  ! Compute composition
  call parser_is_defined('phi',isdef)
  if (isdef) then
     call parser_read('phi',phi)
     Yd_u = 1.0_WP/(1.0_WP+Yo_2/(1.0_WP-Yo_2)*(1.0_WP+phi/nu))
     Yo_u = Yd_u * Yo_2/(1.0_WP-Yo_2)
     Yf_u = Yo_u*phi/nu
  else
     call parser_read('zmix', ZMIX)
     Yf_u = Yf_1*ZMIX + Yf_2*(1.0_WP-ZMIX)
     Yo_u = Yo_1*ZMIX + Yo_2*(1.0_WP-ZMIX)
     phi  = nu*Yf_u/(Yo_u+epsilon(Yo_u))
  end if
  
  ! Allocate the array data
  nvar = 10
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U    => data(:,:,:,1);  names(1)  = 'U'
  V    => data(:,:,:,2);  names(2)  = 'V'
  W    => data(:,:,:,3);  names(3)  = 'W'
  P    => data(:,:,:,4);  names(4)  = 'P'
  RHO  => data(:,:,:,5);  names(5)  = 'RHO'
  dRHO => data(:,:,:,6);  names(6)  = 'dRHO'
  Yf   => data(:,:,:,7);  names(7)  = 'Yf'
  Yo   => data(:,:,:,8);  names(8)  = 'Yo'
  Yp   => data(:,:,:,9);  names(9)  = 'Yp'
  T    => data(:,:,:,10); names(10) = 'T'
  
  ! Create them
  U    = 0.0_WP
  V    = 0.0_WP
  W    = 0.0_WP
  P    = 0.0_WP
  RHO  = 1.0_WP
  dRHO = 0.0_WP
  Yp   = 0.0_WP
  
  ! Some fully burning conditions
  xloc = xm(7*nx/8)
  dx   = 4.0_WP*L/real(nx,WP)
  do i=1,nx
     Yf(i,:,:) = Yf_u * 0.5_WP*(1.0_WP-tanh((xm(i)-xloc)/dx))
     Yo(i,:,:) = Yo_u * 0.5_WP*(1.0_WP-tanh((xm(i)-xloc)/dx))
     T (i,:,:) = 298.0_WP + 2702.0_WP * 0.5_WP*(1.0_WP+tanh((xm(i)-xloc)/dx))
  end do

  return
end subroutine laminar_flame_data


