module helium_plume
  use precision
  use param
  use string
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  
end module helium_plume

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine helium_plume_grid
  use helium_plume
  use gambit
  use parser
  use math
  implicit none
  integer :: i,j,k
  
  ! Assumes the periodicity
  xper = 0
  yper = 0
  zper = 1
  
  ! Cylindrical
  icyl = 1
  
  ! Read in the size of the domain
  call parser_read('nz',nz)
  
  ! Read mesh file
  call gambit_readmesh(0)
  
  write(*,'(a,i4,a,i4,a6)') " Mesh created :",nx," x ",ny," cells"
  
  ! Create the z direction mesh
  allocate(z(nz+1))
  zL = twopi
  if (nz.eq.1) zL = twopi/128.0_WP
  do k=1,nz+1
     z(k) = real(k-1,WP)*zL/real(nz,WP) - 0.5_WP*zL
  end do
  
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
  
  ! Force a thin wall for the air flow
  do j=1,ny-1
     do i=1,nx
        if ((xm(i).lt.-1.74_WP) .and. (ym(j).lt.2.30_WP) .and. (ym(j+1).gt.2.30_WP)) then
           mask(i,j) = 1
        end if
     end do
  end do
  
  ! Force some mask values
  do j=1,ny
     do i=1,nx
        ! Helium
        if ((xm(i).lt.0.0_WP) .and.(ym(j).lt.0.5_WP)) then
           mask(i,j) = 3
        end if
        ! Air
        !if ((xm(i).lt.-1.74_WP) .and. (ym(j).gt.2.30_WP) .and. (ym(j).lt.2.91_WP)) then
        !   mask(i,j) = 3
        !end if
     end do
  end do
  
  ! Force stair stepping at the exit
  do j=1,ny
     do i=1,nx
        if ( (y(j+1).gt.(x(i+1)-3.55_WP)/(4.56_WP-3.55_WP)*(1.30_WP-2.91_WP)+2.91_WP) &
             .and. (ym(j).gt.1.30_WP)) &
             mask(i,j) = 1
     end do
  end do
  
  return
end subroutine helium_plume_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine helium_plume_data
  use helium_plume
  use parser
  use random
  implicit none
  
  integer  :: i,j
  real(WP) :: U_air,U_he,rho_air,rho_he
  
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
  
  ! Create them
  U = 0.0_WP
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  ZMIX = 0.0_WP
  
  ! Read values
  call parser_read('Air velocity',U_air)
  call parser_read('Helium velocity',U_he)  
  call parser_read('Air density',rho_air)
  call parser_read('Helium density',rho_he)
  RHO  = rho_air
  dRHO = 0.0_WP
  
  ! Preset the values in the inlets
  do j=1,ny
     do i=1,nx
        ! Helium
        if ((xm(i).lt.0.0_WP) .and.(ym(j).lt.0.5_WP)) then
           U   (i:i+1,j,:) = U_he
           RHO (i,j,:) = rho_he
           ZMIX(i,j,:) = 1.0_WP
        end if
        ! Air
        if ((xm(i).lt.-1.74_WP) .and. (ym(j).gt.2.30_WP) .and. (ym(j).lt.2.91_WP)) then
           U   (i:i+1,j,:) = U_air
           RHO (i,j,:) = rho_air
           ZMIX(i,j,:) = 0.0_WP
        end if
     end do
  end do
  
  return
end subroutine helium_plume_data

! ========================= !
! Create the variable array !
! ========================= !
subroutine helium_plume_optdata
  use helium_plume
  use parser
  use random
  implicit none
  integer :: iopt

  ! Soot Opt Data
  call parser_read('Optdata with SGS',iopt,0)
  if (iopt.eq.1) then
     nod = 4
     allocate(OD(nx,ny,nz,nod))
     allocate(OD_names(nod))
     OD_names(1) = 'LM_VEL'
     OD_names(2) = 'MM_VEL'
     OD_names(3) = 'LM_ZMIX'
     OD_names(4) = 'MM_ZMIX'
     OD = 0.0_WP
  end if
  
  return
end subroutine helium_plume_optdata
