module gambit
  use precision
  use param
  use string
  
  ! Size of the domain along z
  real(WP) :: zL
  
  ! Open boundaries/walls
  integer :: open_boundary,twall,bwall
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: rhoU
  real(WP), dimension(:,:,:), pointer :: rhoV
  real(WP), dimension(:,:,:), pointer :: rhoW
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: ZMIX2
  real(WP), dimension(:,:,:), pointer :: ZVAR
  real(WP), dimension(:,:,:), pointer :: T
  real(WP), dimension(:,:,:), pointer :: H
  real(WP), dimension(:,:,:), pointer :: G
  real(WP), dimension(:,:,:), pointer :: PROG
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  
end module gambit

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine gambit_grid
  use gambit
  use parser
  implicit none
  integer :: i,j,k
  real(WP) :: twopi
  
  ! Create 2*pi
  twopi = 2.0_WP * acos(-1.0_WP)
  
  ! Assumes the periodicity
  xper = 0
  yper = 0
  zper = 1
  
  ! Cartesian / cylindrical
  call parser_read('Cylindrical coords.',icyl)
  
  ! Read in the size of the domain
  call parser_read('nz',nz)
  
  ! Read mesh file
  if (icyl.eq.1) then
     call parser_read('Open boundaries',open_boundary)
  else
     call parser_read('Top wall',twall)
     call parser_read('Bottom wall',bwall)
     if (twall.eq.0 .and. bwall.eq.0) then
        open_boundary=1
     else if (twall.eq.1 .and. bwall.eq.1) then
        open_boundary=0
     else if (twall.eq.0 .and. bwall.eq.1) then
        open_boundary=2
     else if (twall.eq.1 .and. bwall.eq.0) then
        open_boundary=3
     end if
  end if
  call gambit_readmesh(open_boundary)

  write(*,'(a,i4,a,i4,a6)') " Mesh created :",nx," x ",ny," cells"
  
  ! Create the z direction mesh
  allocate(z(nz+1))
  if (icyl.eq.1) then
     call parser_read('Size along z',zL,twopi)
     do k=1,nz+1
        z(k) = real(k-1,WP)*zL/real(nz,WP)
     end do
  else
     call parser_read('Size along z',zL,1.0_WP)
     do k=1,nz+1
        z(k) = real(k-1,WP)*zL/real(nz,WP) - 0.5_WP*zL
     end do
  end if
  
  ! Cell centered mesh
  allocate(xm(nx),ym(ny),zm(nz))
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do
  
  return
end subroutine gambit_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine gambit_data
  use gambit
  use parser
  use random
  implicit none
  character(str_medium) :: data_type
  real(WP) :: rho_init,trad,thick,dist,xloc
  integer :: i,j,k
  logical :: isdef
  
  call parser_read('Data type',data_type)
  
  select case (trim(adjustl(data_type)))
  case ('hot')
     ! Allocate the array data
     nvar = 8
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
     PROG => data(:,:,:,6); names(6) = 'PROG'
     RHO  => data(:,:,:,7); names(7) = 'RHO'
     dRHO => data(:,:,:,8); names(8) = 'dRHO'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P    = 0.0_WP
     ZMIX = 0.0_WP
     PROG = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
  case ('ZMIX^2')
     ! Allocate the array data
     nvar = 9
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U     => data(:,:,:,1); names(1) = 'U'
     V     => data(:,:,:,2); names(2) = 'V'
     W     => data(:,:,:,3); names(3) = 'W'
     P     => data(:,:,:,4); names(4) = 'P'
     ZMIX  => data(:,:,:,5); names(5) = 'ZMIX'
     ZMIX2 => data(:,:,:,6); names(6) = 'ZMIX^2'
     PROG  => data(:,:,:,7); names(7) = 'PROG'
     RHO   => data(:,:,:,8); names(8) = 'RHO'
     dRHO  => data(:,:,:,9); names(9) = 'dRHO'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
     ZMIX  = 0.0_WP
     ZMIX2 = 0.0_WP
     PROG  = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
  case ('ZVAR')
     ! Allocate the array data
     nvar = 9
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
     ZVAR => data(:,:,:,6); names(6) = 'ZVAR'
     PROG => data(:,:,:,7); names(7) = 'PROG'
     RHO  => data(:,:,:,8); names(8) = 'RHO'
     dRHO => data(:,:,:,9); names(9) = 'dRHO'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
     ZMIX = 0.0_WP
     ZVAR = 0.0_WP
     PROG = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
  case ('premixed')
     ! Allocate the array data
     nvar = 10
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ZMIX => data(:,:,:,5);  names(5)  = 'ZMIX'
     T    => data(:,:,:,6);  names(6)  = 'T'
     H    => data(:,:,:,7);  names(7)  = 'H'
     G    => data(:,:,:,8);  names(8)  = 'LVLSET'
     RHO  => data(:,:,:,9);  names(9)  = 'RHO'
     dRHO => data(:,:,:,10); names(10) = 'dRHO'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
     ZMIX = 0.0_WP
     T = 300.0_WP
     H = 0.0_WP
     G = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
  case ('cold')
     ! Allocate the array data
     nvar = 4
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
  case ('cold mixing')
     ! Allocate the array data
     nvar = 5
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     ! Link the pointers
     U    => data(:,:,:,1);  names(1)  = 'U'
     V    => data(:,:,:,2);  names(2)  = 'V'
     W    => data(:,:,:,3);  names(3)  = 'W'
     P    => data(:,:,:,4);  names(4)  = 'P'
     ZMIX => data(:,:,:,5);  names(5)  = 'ZMIX'
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
     ZMIX = 0.0_WP
  end select
  
  return
end subroutine gambit_data

! ========================= !
! Create the variable array !
! ========================= !
subroutine gambit_optdata
  use gambit
  use parser
  use random
  implicit none
  integer :: isoot

  ! Soot Opt Data
  call parser_read('Opdata with soot',isoot,0)
  if (isoot.eq.1) then
     nod = 4
     allocate(OD(nx,ny,nz,nod))
     allocate(OD_names(nod))
     OD_names(1) = 'SOOT_W1'
     OD_names(2) = 'SOOT_C1'
     OD_names(3) = 'SOOT_S1'
     OD_names(4) = 'SOOT_H1'
     OD = 0.0_WP
  end if
  
  return
end subroutine gambit_optdata
  
