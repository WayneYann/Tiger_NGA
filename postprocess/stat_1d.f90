module stat_1d
  use stat
  implicit none
  
  ! Sampled variables
  integer :: stat_nvar
  character(len=str_medium), dimension(:), pointer :: stat_name

  ! Sampled variables - conditional
  integer :: stat_nvar_cond
  character(len=str_medium), dimension(:), pointer :: stat_name_cond
  
  ! ----------------------------------------------------
  ! 1D-Y STATISTICS 
  !  -> as a function of y
  !  -> at a given station in x
  !
  ! nxloc = -1  for no station (x periodic)
  ! xloc belongs to [ x(iloc)  ,x(iloc+1)   [
  ! xloc belongs to [ xm(imloc),xm(imloc+1) [
  ! xloc = cx(iloc) *x(iloc)  + (1-cx(iloc)) *x(iloc+1)
  ! xmloc = cxm(iloc)*xm(iloc) + (1-cxm)iloc))*xm(iloc+1)
  ! -----------------------------------------------------
  logical :: allx
  integer :: nxloc
  real(WP), dimension(:),     allocatable :: xloc
  integer,  dimension(:),     allocatable :: iloc
  integer,  dimension(:),     allocatable :: iloc_ifwrite
  integer,  dimension(:),     allocatable :: imloc
  real(WP), dimension(:,:),   allocatable :: cxloc
  real(WP), dimension(:,:),   allocatable :: cxmloc
  real(WP), dimension(:,:,:), allocatable :: stat_y
  real(WP), dimension(:,:,:), allocatable :: buf_y
  real(WP) :: Delta_ty
  
  ! Conditional
  real(WP), dimension(:,:,:), allocatable :: stat_y_cond
  real(WP), dimension(:,:,:), allocatable :: buf_y_cond
  real(WP), dimension(:,:),   allocatable :: nSamples_cond
  real(WP), dimension(:,:),   allocatable :: buf_nSamples_cond
  integer :: nbins_cond = 200

  ! -----------------------------------------------------
  ! 1D-X STATISTICS 
  !  -> as a function of x
  !  -> at a given station in y
  !
  ! nyloc = -1  for no station (y periodic)
  ! yloc belongs to [ y(jloc)  ,y(jloc+1)   [
  ! yloc belongs to [ ym(jmloc),ym(jmloc+1) [
  ! yloc = cy(jloc) *y(jloc)  + (1-cy(jloc)) *y(jloc+1)
  ! ymloc = cym(jloc)*ym(jloc) + (1-cym(jloc))*ym(jloc+1)
  ! -----------------------------------------------------
  integer :: nyloc
  real(WP), dimension(:),     allocatable :: yloc
  integer,  dimension(:),     allocatable :: jloc
  integer,  dimension(:),     allocatable :: jloc_ifwrite
  integer,  dimension(:),     allocatable :: jmloc
  real(WP), dimension(:,:),   allocatable :: cyloc
  real(WP), dimension(:,:),   allocatable :: cymloc
  real(WP), dimension(:,:,:), allocatable :: stat_x
  real(WP), dimension(:,:,:), allocatable :: buf_x
  real(WP) :: Delta_tx
  
contains
  
  ! ================================= !
  ! General initialization of STAT-1D !
  ! ================================= !
  subroutine stat_1d_init
    use data
    use combustion
    use string
    implicit none
    
    integer :: isc,ns
    character(len=str_short) :: name
    
    ! Return if init already done
    if (associated(stat_name)) return
    
    ! Create a timer
    call timing_create('stat_1d')
    
    stat_nvar = 0
    ! Density 
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+3
    ! Velocity
    stat_nvar = stat_nvar+9
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+6
    ! Scalars
    do isc=1,nscalar
       name = SC_name(isc)
       if (name(1:2).ne.'S_') then
          stat_nvar = stat_nvar+2
          if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
       end if
    end do
    ! Combustion
    select case (trim(chemistry))
    case ('chemtable')
       stat_nvar = stat_nvar+32
    end select
    ! Soot
    if (use_soot) then
       stat_nvar = stat_nvar + 18
       if (.not.use_pah) stat_nvar = stat_nvar + 4
    end if
    
    ! Allocate
    allocate(stat_name(stat_nvar))

    ! Conditional statistics
    if (cond_stat) then
       if (trim(chemistry).ne.'chemtable') then
          call die('Conditional statistics only for nonpremixed chemtable')
       end if
       stat_nvar_cond = 0
       ! Density
       stat_nvar_cond = stat_nvar_cond + 1
       ! Scalars
       do isc=1,nscalar
          name = SC_name(isc)
          if (name(1:2).ne.'S_') then
             stat_nvar_cond = stat_nvar_cond + 4
          end if
       end do
       ! Scalar Dissipation Rate
       if (isc_ZMIX.ne.0) stat_nvar_cond = stat_nvar_cond + 4
       ! Combustion
       stat_nvar_cond = stat_nvar_cond + 32
       ! Soot
       if (use_soot) then
          stat_nvar_cond = stat_nvar_cond + 18
          if (.not.use_pah) then
             stat_nvar_cond = stat_nvar_cond + 4
          end if
       end if
       
       ! Allocate
       allocate(stat_name_cond(stat_nvar_cond))
    end if
       
    ns = 0

    ! Density statistics
    if (trim(chemistry).ne.'none') then
       stat_name(ns+1) = 'RHO'
       stat_name(ns+2) = 'RHO^2'
       stat_name(ns+3) = 'RHO_F'
       ns = ns+3
    end if

    ! Velocity statistics
    stat_name(ns+1) = 'U'
    stat_name(ns+2) = 'V'
    stat_name(ns+3) = 'W'
    stat_name(ns+4) = 'U^2'
    stat_name(ns+5) = 'V^2'
    stat_name(ns+6) = 'W^2'
    stat_name(ns+7) = 'U^3'
    stat_name(ns+8) = 'V^3'
    stat_name(ns+9) = 'W^3'
    ns = ns+9
    if (trim(chemistry).ne.'none') then
       stat_name(ns+1) = 'rhoU'
       stat_name(ns+2) = 'rhoV'
       stat_name(ns+3) = 'rhoW'
       stat_name(ns+4) = 'rhoU^2'
       stat_name(ns+5) = 'rhoV^2'
       stat_name(ns+6) = 'rhoW^2'
       ns = ns+6
    end if
    
    ! Scalars
    do isc=1,nscalar
       name = SC_name(isc)
       if (name(1:2).ne.'S_') then
          stat_name(ns+1) = 'SC-'  // trim(adjustl(sc_name(isc)))
          stat_name(ns+2) = 'SC^2-'// trim(adjustl(sc_name(isc)))
          ns = ns+2
          if (trim(chemistry).ne.'none') then
             stat_name(ns+1) = 'rhoS-'  // trim(adjustl(sc_name(isc)))
             stat_name(ns+2) = 'rhoS^2-'// trim(adjustl(sc_name(isc)))
             ns = ns+2
          end if
       end if
    end do
    
    ! Combustion
    if (trim(chemistry).eq.'chemtable') then
       stat_name(ns+ 1) = 'T'
       stat_name(ns+ 2) = 'rhoT'
       stat_name(ns+ 3) = 'T^2'
       stat_name(ns+ 4) = 'rhoT^2'
       stat_name(ns+ 5) = 'Y_F'
       stat_name(ns+ 6) = 'rhoY_F'
       stat_name(ns+ 7) = 'Y_F^2'
       stat_name(ns+ 8) = 'rhoY_F^2'
       stat_name(ns+ 9) = 'Y_O2'
       stat_name(ns+10) = 'rhoY_O2'
       stat_name(ns+11) = 'Y_O2^2'
       stat_name(ns+12) = 'rhoY_O2^2'
       stat_name(ns+13) = 'Y_CO'
       stat_name(ns+14) = 'rhoY_CO'
       stat_name(ns+15) = 'Y_CO^2'
       stat_name(ns+16) = 'rhoY_CO^2'
       stat_name(ns+17) = 'Y_CO2'
       stat_name(ns+18) = 'rhoY_CO2'
       stat_name(ns+19) = 'Y_CO2^2'
       stat_name(ns+20) = 'rhoY_CO2^2'
       stat_name(ns+21) = 'Y_H2'
       stat_name(ns+22) = 'rhoY_H2'
       stat_name(ns+23) = 'Y_H2^2'
       stat_name(ns+24) = 'rhoY_H2^2'
       stat_name(ns+25) = 'Y_H2O'
       stat_name(ns+26) = 'rhoY_H2O'
       stat_name(ns+27) = 'Y_H2O^2'
       stat_name(ns+28) = 'rhoY_H2O^2'
       stat_name(ns+29) = 'Y_OH'
       stat_name(ns+30) = 'rhoY_OH'
       stat_name(ns+31) = 'Y_OH^2'
       stat_name(ns+32) = 'rhoY_OH^2'
       ns = ns+32
    end if
    
    ! Soot
    if (use_soot) then
       stat_name(ns+ 1) = 'fV'
       stat_name(ns+ 2) = 'fV^2'
       stat_name(ns+ 3) = 'N'
       stat_name(ns+ 4) = 'N^2'
       stat_name(ns+ 5) = 'dp'
       stat_name(ns+ 6) = 'dp^2'
       stat_name(ns+ 7) = 'np'
       stat_name(ns+ 8) = 'np^2'
       stat_name(ns+ 9) = 'intermit'
       stat_name(ns+10) = 'intermit^2'
       stat_name(ns+11) = 'dNdt_nucl'
       stat_name(ns+12) = 'dNdt_coag'
       stat_name(ns+13) = 'dNdt_ox'
       stat_name(ns+14) = 'dNdt_frag'
       stat_name(ns+15) = 'dfvdt_nucl'
       stat_name(ns+16) = 'dfvdt_cond'
       stat_name(ns+17) = 'dfvdt_sg'
       stat_name(ns+18) = 'dfvdt_ox'
       ns = ns+18
       if (.not.use_pah) then
          stat_name(ns+1) = 'Y_PAH'
          stat_name(ns+2) = 'rhoY_PAH'
          stat_name(ns+3) = 'Y_PAH^2'
          stat_name(ns+4) = 'rhoY_PAH^2'
          ns = ns+4
       end if
    end if

    ! Conditional statistics
    if (cond_stat) then
       stat_nvar_cond = 0
       ! Density
       stat_name_cond(stat_nvar_cond+1) = 'RHO|Z'
       stat_nvar_cond = stat_nvar_cond + 1
       ! Scalars
       do isc=1,nscalar
          name = SC_name(isc)
          if (name(1:2).ne.'S_') then
             stat_name_cond(stat_nvar_cond+1) = 'SC-'     // trim(adjustl(sc_name(isc))) // '|Z'
             stat_name_cond(stat_nvar_cond+2) = 'SC^2-'   // trim(adjustl(sc_name(isc))) // '|Z'
             stat_name_cond(stat_nvar_cond+3) = 'rhoS-'   // trim(adjustl(sc_name(isc))) // '|Z'
             stat_name_cond(stat_nvar_cond+4) = 'rhoS^2-' // trim(adjustl(sc_name(isc))) // '|Z'
             stat_nvar_cond = stat_nvar_cond + 4
          end if
       end do
       ! Scalar Dissipation Rate
       if (isc_ZMIX.ne.0) then
          stat_name_cond(stat_nvar_cond+1) = 'CHI'
          stat_name_cond(stat_nvar_cond+2) = 'CHI^2'
          stat_name_cond(stat_nvar_cond+3) = 'rhoCHI'
          stat_name_cond(stat_nvar_cond+4) = 'rhoCHI^2'
          stat_nvar_cond = stat_nvar_cond + 4
       end if
       ! Combustion
       stat_name_cond(stat_nvar_cond+ 1) = 'T|Z'
       stat_name_cond(stat_nvar_cond+ 2) = 'rhoT|Z'
       stat_name_cond(stat_nvar_cond+ 3) = 'T^2|Z'
       stat_name_cond(stat_nvar_cond+ 4) = 'rhoT^2|Z'
       stat_name_cond(stat_nvar_cond+ 5) = 'Y_F|Z'
       stat_name_cond(stat_nvar_cond+ 6) = 'rhoY_F|Z'
       stat_name_cond(stat_nvar_cond+ 7) = 'Y_F^2|Z'
       stat_name_cond(stat_nvar_cond+ 8) = 'rhoY_F^2|Z'
       stat_name_cond(stat_nvar_cond+ 9) = 'Y_O2|Z'
       stat_name_cond(stat_nvar_cond+10) = 'rhoY_O2|Z'
       stat_name_cond(stat_nvar_cond+11) = 'Y_O2^2|Z'
       stat_name_cond(stat_nvar_cond+12) = 'rhoY_O2^2|Z'
       stat_name_cond(stat_nvar_cond+13) = 'Y_CO|Z'
       stat_name_cond(stat_nvar_cond+14) = 'rhoY_CO|Z'
       stat_name_cond(stat_nvar_cond+15) = 'Y_CO^2|Z'
       stat_name_cond(stat_nvar_cond+16) = 'rhoY_CO^2|Z'
       stat_name_cond(stat_nvar_cond+17) = 'Y_CO2|Z'
       stat_name_cond(stat_nvar_cond+18) = 'rhoY_CO2|Z'
       stat_name_cond(stat_nvar_cond+19) = 'Y_CO2^2|Z'
       stat_name_cond(stat_nvar_cond+20) = 'rhoY_CO2^2|Z'
       stat_name_cond(stat_nvar_cond+21) = 'Y_H2|Z'
       stat_name_cond(stat_nvar_cond+22) = 'rhoY_H2|Z'
       stat_name_cond(stat_nvar_cond+23) = 'Y_H2^2|Z'
       stat_name_cond(stat_nvar_cond+24) = 'rhoY_H2^2|Z'
       stat_name_cond(stat_nvar_cond+25) = 'Y_H2O|Z'
       stat_name_cond(stat_nvar_cond+26) = 'rhoY_H2O|Z'
       stat_name_cond(stat_nvar_cond+27) = 'Y_H2O^2|Z'
       stat_name_cond(stat_nvar_cond+28) = 'rhoY_H2O^2|Z'
       stat_name_cond(stat_nvar_cond+29) = 'Y_OH|Z'
       stat_name_cond(stat_nvar_cond+30) = 'rhoY_OH|Z'
       stat_name_cond(stat_nvar_cond+31) = 'Y_OH^2|Z'
       stat_name_cond(stat_nvar_cond+32) = 'rhoY_OH^2|Z'
       stat_nvar_cond = stat_nvar_cond + 32
       ! Soot
       if (use_soot) then
          stat_name_cond(stat_nvar_cond+ 1) = 'fV|Z'
          stat_name_cond(stat_nvar_cond+ 2) = 'fV^2|Z'
          stat_name_cond(stat_nvar_cond+ 3) = 'N|Z'
          stat_name_cond(stat_nvar_cond+ 4) = 'N^2|Z'
          stat_name_cond(stat_nvar_cond+ 5) = 'dp|Z'
          stat_name_cond(stat_nvar_cond+ 6) = 'dp^2|Z'
          stat_name_cond(stat_nvar_cond+ 7) = 'np|Z'
          stat_name_cond(stat_nvar_cond+ 8) = 'np^2|Z'
          stat_name_cond(stat_nvar_cond+ 9) = 'intermit|Z'
          stat_name_cond(stat_nvar_cond+10) = 'intermit^2|Z'
          stat_name_cond(stat_nvar_cond+11) = 'dNdt_nucl|Z'
          stat_name_cond(stat_nvar_cond+12) = 'dNdt_coag|Z'
          stat_name_cond(stat_nvar_cond+13) = 'dNdt_ox|Z'
          stat_name_cond(stat_nvar_cond+14) = 'dNdt_frag|Z'
          stat_name_cond(stat_nvar_cond+15) = 'dfvdt_nucl|Z'
          stat_name_cond(stat_nvar_cond+16) = 'dfvdt_cond|Z'
          stat_name_cond(stat_nvar_cond+17) = 'dfvdt_sg|Z'
          stat_name_cond(stat_nvar_cond+18) = 'dfvdt_ox|Z'
          stat_nvar_cond = stat_nvar_cond + 18
          if (.not.use_pah) then
             stat_name_cond(stat_nvar_cond+1) = 'Y_PAH|Z'
             stat_name_cond(stat_nvar_cond+2) = 'rhoY_PAH|Z'
             stat_name_cond(stat_nvar_cond+3) = 'Y_PAH^2|Z'
             stat_name_cond(stat_nvar_cond+4) = 'rhoY_PAH^2|Z'
             stat_nvar_cond = stat_nvar_cond + 4
          end if
       end if
    end if
    
    return
  end subroutine stat_1d_init
  
  ! ============================== !
  ! Setup the 1D metric            !
  !  -> indices                    !
  !  -> interpolation coefficients !
  ! ============================== !
  subroutine stat_1dx_metric
    implicit none
    integer :: i,loc
    
    ! Allocate the arrays
    allocate(iloc(nxloc))
    allocate(iloc_ifwrite(nxloc))
    allocate(imloc(nxloc))
    allocate(cxloc(nxloc,2))
    allocate(cxmloc(nxloc,2))
    
    ! Locate the first i index to the left
    do loc=1,nxloc
       if (x(imin_).gt.xloc(loc) .or. x(imax_+1).le.xloc(loc)) then
          iloc(loc) = imin_
          imloc(loc) = imin_
          cxloc(loc,:) = 0.0_WP
          cxmloc(loc,:) = 0.0_WP
          iloc_ifwrite(loc) = 0
       else
          i = imin_
          do while(i.le.imax_ .and. x(i+1).le.xloc(loc))
             i = i+1
          end do
          iloc(loc) = i
          cxloc(loc,1) = (x(iloc(loc)+1)-xloc(loc))/(x(iloc(loc)+1)-x(iloc(loc)))
          cxloc(loc,2) = 1.0_WP-cxloc(loc,1)
          
          if (xloc(loc).ge.xm(i)) then
             imloc(loc) = i
          else
             imloc(loc) = i-1
          end if
          cxmloc(loc,1) = (xm(imloc(loc)+1)-xloc(loc))/(xm(imloc(loc)+1)-xm(imloc(loc)))
          cxmloc(loc,2) = 1.0_WP-cxmloc(loc,1)
          iloc_ifwrite(loc) = 1
       end if
    end do
    
    return
  end subroutine stat_1dx_metric
  
  subroutine stat_1dy_metric
    implicit none
    integer :: j,loc
    
    ! Allocate the arrays
    allocate(jloc(nyloc))
    allocate(jloc_ifwrite(nyloc))
    allocate(jmloc(nyloc))
    allocate(cyloc(nyloc,2))
    allocate(cymloc(nyloc,2))
    
    ! Locate the first j index to the bottom - cell face
    do loc=1,nyloc
       if (y(jmin_).gt.yloc(loc) .or. y(jmax_+1).le.yloc(loc)) then
          jloc(loc) = jmin_
          jmloc(loc) = jmin_
          cyloc(loc,:) = 0.0_WP
          cymloc(loc,:) = 0.0_WP
          jloc_ifwrite(loc) = 0
       else
          j = jmin_
          do while(j.le.jmax_ .and. y(j+1).le.yloc(loc))
             j = j+1
          end do
          jloc(loc) = j
          cyloc(loc,1) = (y(jloc(loc)+1)-yloc(loc))/(y(jloc(loc)+1)-y(jloc(loc)))
          cyloc(loc,2) = 1.0_WP-cyloc(loc,1)
          
          if (yloc(loc).ge.ym(j)) then
             jmloc(loc) = j
          else
             jmloc(loc) = j-1
          end if
          cymloc(loc,1) = (ym(jmloc(loc)+1)-yloc(loc))/(ym(jmloc(loc)+1)-ym(jmloc(loc)))
          cymloc(loc,2) = 1.0_WP-cymloc(loc,1)
          jloc_ifwrite(loc) = 1
       end if
    end do
    
    return
  end subroutine stat_1dy_metric
  
end module stat_1d


! ================================== !
! Initialize the 1D statistic module !
!   -> at a given x                  !
!   -> as a function of y            !
! ================================== !
subroutine stat_1dy_init
  use stat_1d
  use parallel
  use parser
  implicit none
  logical :: isdef,file_is_there
  character(len=str_medium) :: text
  
  ! General 1D non specific initialization
  call stat_1d_init
  
  ! Start the timer
  call timing_start('stat_1d')
  
  ! No stat at a given x (as a function of y) if y periodic
  if (yper.eq.1) call die('stat_1dy_init: 1D-y statistics impossible (y is periodic)')
  
  ! Any x locations ?
  call parser_is_defined('Statistics locations x',isdef)
  if (isdef) then
     call parser_read('Statistics locations x',text)
     if (trim(text).eq.'all') then
        if (xper.ne.1) call die('stat_1dy_init: all incompatible with xper')
        ! Create all locations
        allx = .true.
        nxloc = nx
        allocate(xloc(nxloc))
        xloc = xm(imin:imax)
     else
        ! Read the locations
        allx = .false.
        call parser_getsize('Statistics locations x',nxloc)
        allocate(xloc(nxloc))
        call parser_read('Statistics locations x',xloc)
     end if
  else
     call die('stat_1dy_init: Statistics locations x not specified')
  end if
  
  ! Create metric
  call stat_1dx_metric

  ! Allocate stat arrays
  allocate(stat_y(nxloc,jmin:jmax,stat_nvar))
  allocate(buf_y (nxloc,jmin:jmax,stat_nvar))
  stat_y = 0.0_WP
  
  ! Open file if it exists
  inquire(file='stat/stat-1Dy',exist=file_is_there)
  if (file_is_there) then
     call stat_1dy_read
  else
     Delta_ty = 0.0_WP
  end if
  
  ! Conditional statistics
  if (cond_stat) then

     ! Allocate stat arrays
     allocate(stat_y_cond(nxloc,nbins_cond,stat_nvar_cond))
     allocate(buf_y_cond(nxloc,nbins_cond,stat_nvar_cond))
     allocate(nSamples_cond(nxloc,nbins_cond))
     allocate(buf_nSamples_cond(nxloc,nbins_cond))
     stat_y_cond = 0.0_WP
     nSamples_cond = 0.0_WP

     ! Open file if it exists
     inquire(file='stat/stat-1Dy-cond',exist=file_is_there)
     if (file_is_there) then
        call stat_1dy_cond_read
     end if
  end if

  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dy_init


! ================================== !
! Initialize the 1D statistic module !
!   -> at a given y                  !
!   -> as a function of x            !
! ================================== !
subroutine stat_1dx_init
  use stat_1d
  use parallel
  use parser
  implicit none
  logical :: isdef,file_is_there
  
  ! General 1D non specific initialization
  call stat_1d_init
  
  ! Start the timer
  call timing_start('stat_1d')
  
  ! No stat at a given y (as a function of x) if x periodic
  if (xper.eq.1) call die('stat_1dx_init: 1D-x statistics impossible (x is periodic)')
  
  ! Any y locations ?
  call parser_is_defined('Statistics locations y',isdef)
  if (isdef) then
     call parser_getsize('Statistics locations y',nyloc)
     allocate(yloc(nyloc))
     call parser_read('Statistics locations y',yloc)
  else 
     call die('stat_1dy_init: Statistics locations y not specified')
  end if
  
  ! Create metric
  call stat_1dy_metric
  
  ! Allocate stat arrays
  allocate(stat_x(nyloc,imin:imax,stat_nvar))
  allocate(buf_x(nyloc,imin:imax,stat_nvar))
  stat_x = 0.0_WP
  
  ! Open file if it exists
  inquire(file='stat/stat-1Dx',exist=file_is_there)
  if (file_is_there) then
     call stat_1dx_read
  else
     Delta_tx = 0.0_WP
  end if
  
  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dx_init


! ======================= !
! Sample the statistics   !
! at the center: xm or ym !
! ======================= !
subroutine stat_1dy_sample
  use stat_1d
  use data
  use combustion
  use time_info
  use memory
  use soot
  use pah
  use parallel  
  use metric_generic
  implicit none
  
  integer  :: i,im,j,k,isc,loc,ns,n,bin
  logical  :: found_loc
  real(WP) :: xcyl,ycyl,zcyl,ucyl,vcyl,wcyl
  character(len=str_short) :: name
  
  ! Start the timer
  call timing_start('stat_1d')
  
  ! Prepare the chemical variable if necessary
  if     (trim(chemistry).eq.'chemtable') then
     call chemtable_lookup('T'    ,tmp1)
     call chemtable_lookup('Y_F'  ,tmp2)
     call chemtable_lookup('Y_O2' ,tmp3)
     call chemtable_lookup('Y_CO' ,tmp4)
     call chemtable_lookup('Y_CO2',tmp5)
     call chemtable_lookup('Y_H2' ,tmp6)
     call chemtable_lookup('Y_H2O',tmp7)
     call chemtable_lookup('Y_OH' ,tmp8)
     if (use_soot .and. .not.use_pah) call chemtable_lookup('Y_PAH',tmp9)
  end if

  ! Prepare the density at the y-faces
  if (trim(chemistry).ne.'none') then
     do i=imin_,imax_
        do j=jmin_,jmax_
           do k=kmin_,kmax_
              tmp9(i,j,k) = sum(interp_sc_y(i,j,:)*RHO(i,j-st2:j+st1,k))
           end do
        end do
     end do
  end if
  
  ! Gather the stats
  do loc=1,nxloc
     i  = iloc (loc)
     im = imloc(loc)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           
           ns = 0

           ! Density
           if (trim(chemistry).ne.'none') then
              stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxmloc(loc,:)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*RHO(im:im+1,j,k)**2)
              stat_y(loc,j,ns+3) = stat_y(loc,j,ns+3) + dt*sum(cxloc (loc,:)*tmp9(i:i+1,j,k))
              ns = ns+3
           end if
           
           ! Velocity
           stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxloc (loc,:)*U(i:i+1  ,j,k))
           stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*V(im:im+1,j,k))
           stat_y(loc,j,ns+3) = stat_y(loc,j,ns+3) + dt*sum(cxmloc(loc,:)*W(im:im+1,j,k))
           stat_y(loc,j,ns+4) = stat_y(loc,j,ns+4) + dt*sum(cxloc (loc,:)*U(i:i+1  ,j,k)**2)
           stat_y(loc,j,ns+5) = stat_y(loc,j,ns+5) + dt*sum(cxmloc(loc,:)*V(im:im+1,j,k)**2)
           stat_y(loc,j,ns+6) = stat_y(loc,j,ns+6) + dt*sum(cxmloc(loc,:)*W(im:im+1,j,k)**2)
           stat_y(loc,j,ns+7) = stat_y(loc,j,ns+7) + dt*sum(cxloc (loc,:)*U(i:i+1  ,j,k)**3)
           stat_y(loc,j,ns+8) = stat_y(loc,j,ns+8) + dt*sum(cxmloc(loc,:)*V(im:im+1,j,k)**3)
           stat_y(loc,j,ns+9) = stat_y(loc,j,ns+9) + dt*sum(cxmloc(loc,:)*W(im:im+1,j,k)**3)
           ns = ns+9
           
           if (trim(chemistry).ne.'none') then
              stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxloc (loc,:)*rhoU(i:i+1  ,j,k))
              stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*rhoV(im:im+1,j,k))
              stat_y(loc,j,ns+3) = stat_y(loc,j,ns+3) + dt*sum(cxmloc(loc,:)*rhoW(im:im+1,j,k))
              stat_y(loc,j,ns+4) = stat_y(loc,j,ns+4) + dt*sum(cxloc (loc,:)*rhoU(i:i+1  ,j,k)*U(i:i+1  ,j,k))
              stat_y(loc,j,ns+5) = stat_y(loc,j,ns+5) + dt*sum(cxmloc(loc,:)*rhoV(im:im+1,j,k)*V(im:im+1,j,k))
              stat_y(loc,j,ns+6) = stat_y(loc,j,ns+6) + dt*sum(cxmloc(loc,:)*rhoW(im:im+1,j,k)*W(im:im+1,j,k))
              ns = ns+6
           end if
           
           ! Scalars
           do isc=1,nscalar
              name = SC_name(isc)
              if (name(1:2).ne.'S_') then
                 stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc))
                 stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)**2)
                 ns = ns+2
                 if (trim(chemistry).ne.'none') then
                    stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)*RHO(im:im+1,j,k))
                    stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)**2*RHO(im:im+1,j,k))
                    ns = ns+2
                 end if
              end if
           end do
           
           ! Combustion
           if (trim(chemistry).eq.'chemtable') then
              stat_y(loc,j,ns+ 1) = stat_y(loc,j,ns+ 1) + dt*sum(cxmloc(loc,:)*tmp1(im:im+1,j,k))
              stat_y(loc,j,ns+ 2) = stat_y(loc,j,ns+ 2) + dt*sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+ 3) = stat_y(loc,j,ns+ 3) + dt*sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 4) = stat_y(loc,j,ns+ 4) + dt*sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+ 5) = stat_y(loc,j,ns+ 5) + dt*sum(cxmloc(loc,:)*tmp2(im:im+1,j,k))
              stat_y(loc,j,ns+ 6) = stat_y(loc,j,ns+ 6) + dt*sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+ 7) = stat_y(loc,j,ns+ 7) + dt*sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 8) = stat_y(loc,j,ns+ 8) + dt*sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+ 9) = stat_y(loc,j,ns+ 9) + dt*sum(cxmloc(loc,:)*tmp3(im:im+1,j,k))
              stat_y(loc,j,ns+10) = stat_y(loc,j,ns+10) + dt*sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+11) = stat_y(loc,j,ns+11) + dt*sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)**2)
              stat_y(loc,j,ns+12) = stat_y(loc,j,ns+12) + dt*sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+13) = stat_y(loc,j,ns+13) + dt*sum(cxmloc(loc,:)*tmp4(im:im+1,j,k))
              stat_y(loc,j,ns+14) = stat_y(loc,j,ns+14) + dt*sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+15) = stat_y(loc,j,ns+15) + dt*sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)**2)
              stat_y(loc,j,ns+16) = stat_y(loc,j,ns+16) + dt*sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+17) = stat_y(loc,j,ns+17) + dt*sum(cxmloc(loc,:)*tmp5(im:im+1,j,k))
              stat_y(loc,j,ns+18) = stat_y(loc,j,ns+18) + dt*sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+19) = stat_y(loc,j,ns+19) + dt*sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)**2)
              stat_y(loc,j,ns+20) = stat_y(loc,j,ns+20) + dt*sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+21) = stat_y(loc,j,ns+21) + dt*sum(cxmloc(loc,:)*tmp6(im:im+1,j,k))
              stat_y(loc,j,ns+22) = stat_y(loc,j,ns+22) + dt*sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+23) = stat_y(loc,j,ns+23) + dt*sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)**2)
              stat_y(loc,j,ns+24) = stat_y(loc,j,ns+24) + dt*sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+25) = stat_y(loc,j,ns+25) + dt*sum(cxmloc(loc,:)*tmp7(im:im+1,j,k))
              stat_y(loc,j,ns+26) = stat_y(loc,j,ns+26) + dt*sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+27) = stat_y(loc,j,ns+27) + dt*sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)**2)
              stat_y(loc,j,ns+28) = stat_y(loc,j,ns+28) + dt*sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+29) = stat_y(loc,j,ns+29) + dt*sum(cxmloc(loc,:)*tmp8(im:im+1,j,k))
              stat_y(loc,j,ns+30) = stat_y(loc,j,ns+30) + dt*sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y(loc,j,ns+31) = stat_y(loc,j,ns+31) + dt*sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)**2)
              stat_y(loc,j,ns+32) = stat_y(loc,j,ns+32) + dt*sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              ns = ns+32
           end if
           
           ! Soot
           if (use_soot) then
              stat_y(loc,j,ns+ 1) = stat_y(loc,j,ns+ 1) + dt*sum(cxmloc(loc,:)*volfrac(im:im+1,j,k))
              stat_y(loc,j,ns+ 2) = stat_y(loc,j,ns+ 2) + dt*sum(cxmloc(loc,:)*volfrac(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 3) = stat_y(loc,j,ns+ 3) + dt*sum(cxmloc(loc,:)*numdens(im:im+1,j,k))
              stat_y(loc,j,ns+ 4) = stat_y(loc,j,ns+ 4) + dt*sum(cxmloc(loc,:)*numdens(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 5) = stat_y(loc,j,ns+ 5) + dt*sum(cxmloc(loc,:)*partdiam(im:im+1,j,k))
              stat_y(loc,j,ns+ 6) = stat_y(loc,j,ns+ 6) + dt*sum(cxmloc(loc,:)*partdiam(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 7) = stat_y(loc,j,ns+ 7) + dt*sum(cxmloc(loc,:)*partaggr(im:im+1,j,k))
              stat_y(loc,j,ns+ 8) = stat_y(loc,j,ns+ 8) + dt*sum(cxmloc(loc,:)*partaggr(im:im+1,j,k)**2)
              stat_y(loc,j,ns+ 9) = stat_y(loc,j,ns+ 9) + dt*sum(cxmloc(loc,:)*intermit(im:im+1,j,k))
              stat_y(loc,j,ns+10) = stat_y(loc,j,ns+10) + dt*sum(cxmloc(loc,:)*intermit(im:im+1,j,k)**2)
              stat_y(loc,j,ns+11) = stat_y(loc,j,ns+11) + dt*sum(cxmloc(loc,:)*Nsrc_nucl(im:im+1,j,k))
              stat_y(loc,j,ns+12) = stat_y(loc,j,ns+12) + dt*sum(cxmloc(loc,:)*Nsrc_coag(im:im+1,j,k))
              stat_y(loc,j,ns+13) = stat_y(loc,j,ns+13) + dt*sum(cxmloc(loc,:)*Nsrc_ox  (im:im+1,j,k))
              stat_y(loc,j,ns+14) = stat_y(loc,j,ns+14) + dt*sum(cxmloc(loc,:)*Nsrc_frag(im:im+1,j,k))
              stat_y(loc,j,ns+15) = stat_y(loc,j,ns+15) + dt*sum(cxmloc(loc,:)*FVsrc_nucl(im:im+1,j,k))
              stat_y(loc,j,ns+16) = stat_y(loc,j,ns+16) + dt*sum(cxmloc(loc,:)*FVsrc_cond(im:im+1,j,k))
              stat_y(loc,j,ns+17) = stat_y(loc,j,ns+17) + dt*sum(cxmloc(loc,:)*FVsrc_sg  (im:im+1,j,k))
              stat_y(loc,j,ns+18) = stat_y(loc,j,ns+18) + dt*sum(cxmloc(loc,:)*FVsrc_ox  (im:im+1,j,k))
              ns = ns+18
              if (.not.use_pah) then
                 stat_y(loc,j,ns+1) = stat_y(loc,j,ns+1) + dt*sum(cxmloc(loc,:)*tmp9(im:im+1,j,k))
                 stat_y(loc,j,ns+2) = stat_y(loc,j,ns+2) + dt*sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)*RHO(im:im+1,j,k))
                 stat_y(loc,j,ns+3) = stat_y(loc,j,ns+3) + dt*sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)**2)
                 stat_y(loc,j,ns+4) = stat_y(loc,j,ns+4) + dt*sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)**2*RHO(im:im+1,j,k))
                 ns = ns + 4
              end if
           end if


           ! Conditional statistics
           if (cond_stat) then
              ! Find correct bin
              bin = min(nbins_cond,max(1,floor(sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc_ZMIX))*nbins_cond)+1))
              ! Increase number of samples by one
              nSamples_cond(loc,bin) = nSamples_cond(loc,bin) + 1.0_WP
              ns = 0
              ! Density
              stat_y_cond(loc,bin,ns+1) = stat_y_cond(loc,bin,ns+1) + sum(cxmloc(loc,:)*RHO(im:im+1,j,k))
              ns = ns + 1
              ! Scalars
              do isc=1,nscalar
                 name = SC_name(isc)
                 if (name(1:2).ne.'S_') then
                    stat_y_cond(loc,bin,ns+1) = stat_y_cond(loc,bin,ns+1) + sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc))
                    stat_y_cond(loc,bin,ns+2) = stat_y_cond(loc,bin,ns+2) + sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)**2)
                    stat_y_cond(loc,bin,ns+3) = stat_y_cond(loc,bin,ns+3) + sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)*RHO(im:im+1,j,k))
                    stat_y_cond(loc,bin,ns+4) = stat_y_cond(loc,bin,ns+4) + sum(cxmloc(loc,:)*SC(im:im+1,j,k,isc)**2*RHO(im:im+1,j,k))
                    ns = ns + 4
                 end if
              end do
              ! Scalar Dissipation Rate
              if (isc_ZMIX.ne.0) then
                 stat_y_cond(loc,bin,ns+1) = stat_y_cond(loc,bin,ns+1) + sum(cxmloc(loc,:)*CHI(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+2) = stat_y_cond(loc,bin,ns+2) + sum(cxmloc(loc,:)*CHI(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+3) = stat_y_cond(loc,bin,ns+3) + sum(cxmloc(loc,:)*CHI(im:im+1,j,k)*RHO(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+4) = stat_y_cond(loc,bin,ns+4) + sum(cxmloc(loc,:)*CHI(im:im+1,j,k)**2*RHO(im:im+1,j,k))
                 ns = ns + 4
              end if
              ! Combustion
              stat_y_cond(loc,bin,ns+ 1) = stat_y_cond(loc,bin,ns+ 1) + sum(cxmloc(loc,:)*tmp1(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 2) = stat_y_cond(loc,bin,ns+ 2) + sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 3) = stat_y_cond(loc,bin,ns+ 3) + sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+ 4) = stat_y_cond(loc,bin,ns+ 4) + sum(cxmloc(loc,:)*tmp1(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 5) = stat_y_cond(loc,bin,ns+ 5) + sum(cxmloc(loc,:)*tmp2(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 6) = stat_y_cond(loc,bin,ns+ 6) + sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 7) = stat_y_cond(loc,bin,ns+ 7) + sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+ 8) = stat_y_cond(loc,bin,ns+ 8) + sum(cxmloc(loc,:)*tmp2(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+ 9) = stat_y_cond(loc,bin,ns+ 9) + sum(cxmloc(loc,:)*tmp3(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+10) = stat_y_cond(loc,bin,ns+10) + sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+11) = stat_y_cond(loc,bin,ns+11) + sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+12) = stat_y_cond(loc,bin,ns+12) + sum(cxmloc(loc,:)*tmp3(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+13) = stat_y_cond(loc,bin,ns+13) + sum(cxmloc(loc,:)*tmp4(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+14) = stat_y_cond(loc,bin,ns+14) + sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+15) = stat_y_cond(loc,bin,ns+15) + sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+16) = stat_y_cond(loc,bin,ns+16) + sum(cxmloc(loc,:)*tmp4(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+17) = stat_y_cond(loc,bin,ns+17) + sum(cxmloc(loc,:)*tmp5(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+18) = stat_y_cond(loc,bin,ns+18) + sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+19) = stat_y_cond(loc,bin,ns+19) + sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+20) = stat_y_cond(loc,bin,ns+20) + sum(cxmloc(loc,:)*tmp5(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+21) = stat_y_cond(loc,bin,ns+21) + sum(cxmloc(loc,:)*tmp6(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+22) = stat_y_cond(loc,bin,ns+22) + sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+23) = stat_y_cond(loc,bin,ns+23) + sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+24) = stat_y_cond(loc,bin,ns+24) + sum(cxmloc(loc,:)*tmp6(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+25) = stat_y_cond(loc,bin,ns+25) + sum(cxmloc(loc,:)*tmp7(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+26) = stat_y_cond(loc,bin,ns+26) + sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+27) = stat_y_cond(loc,bin,ns+27) + sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+28) = stat_y_cond(loc,bin,ns+28) + sum(cxmloc(loc,:)*tmp7(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+29) = stat_y_cond(loc,bin,ns+29) + sum(cxmloc(loc,:)*tmp8(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+30) = stat_y_cond(loc,bin,ns+30) + sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)*RHO(im:im+1,j,k))
              stat_y_cond(loc,bin,ns+31) = stat_y_cond(loc,bin,ns+31) + sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)**2)
              stat_y_cond(loc,bin,ns+32) = stat_y_cond(loc,bin,ns+32) + sum(cxmloc(loc,:)*tmp8(im:im+1,j,k)**2*RHO(im:im+1,j,k))
              ns = ns + 32
              ! Soot
              if (use_soot) then
                 stat_y_cond(loc,bin,ns+ 1) = stat_y_cond(loc,bin,ns+ 1) + sum(cxmloc(loc,:)*volfrac(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+ 2) = stat_y_cond(loc,bin,ns+ 2) + sum(cxmloc(loc,:)*volfrac(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+ 3) = stat_y_cond(loc,bin,ns+ 3) + sum(cxmloc(loc,:)*numdens(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+ 4) = stat_y_cond(loc,bin,ns+ 4) + sum(cxmloc(loc,:)*numdens(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+ 5) = stat_y_cond(loc,bin,ns+ 5) + sum(cxmloc(loc,:)*partdiam(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+ 6) = stat_y_cond(loc,bin,ns+ 6) + sum(cxmloc(loc,:)*partdiam(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+ 7) = stat_y_cond(loc,bin,ns+ 7) + sum(cxmloc(loc,:)*partaggr(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+ 8) = stat_y_cond(loc,bin,ns+ 8) + sum(cxmloc(loc,:)*partaggr(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+ 9) = stat_y_cond(loc,bin,ns+ 9) + sum(cxmloc(loc,:)*intermit(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+10) = stat_y_cond(loc,bin,ns+10) + sum(cxmloc(loc,:)*intermit(im:im+1,j,k)**2)
                 stat_y_cond(loc,bin,ns+11) = stat_y_cond(loc,bin,ns+11) + sum(cxmloc(loc,:)*Nsrc_nucl(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+12) = stat_y_cond(loc,bin,ns+12) + sum(cxmloc(loc,:)*Nsrc_coag(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+13) = stat_y_cond(loc,bin,ns+13) + sum(cxmloc(loc,:)*Nsrc_ox  (im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+14) = stat_y_cond(loc,bin,ns+14) + sum(cxmloc(loc,:)*Nsrc_frag(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+15) = stat_y_cond(loc,bin,ns+15) + sum(cxmloc(loc,:)*FVsrc_nucl(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+16) = stat_y_cond(loc,bin,ns+16) + sum(cxmloc(loc,:)*FVsrc_cond(im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+17) = stat_y_cond(loc,bin,ns+17) + sum(cxmloc(loc,:)*FVsrc_sg  (im:im+1,j,k))
                 stat_y_cond(loc,bin,ns+18) = stat_y_cond(loc,bin,ns+18) + sum(cxmloc(loc,:)*FVsrc_ox  (im:im+1,j,k))
                 ns = ns + 18
                 if (.not.use_pah) then
                    stat_y_cond(loc,bin,ns+1) = stat_y_cond(loc,bin,ns+1) + sum(cxmloc(loc,:)*tmp9(im:im+1,j,k))
                    stat_y_cond(loc,bin,ns+2) = stat_y_cond(loc,bin,ns+2) + sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)*RHO(im:im+1,j,k))
                    stat_y_cond(loc,bin,ns+3) = stat_y_cond(loc,bin,ns+3) + sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)**2)
                    stat_y_cond(loc,bin,ns+4) = stat_y_cond(loc,bin,ns+4) + sum(cxmloc(loc,:)*tmp9(im:im+1,j,k)**2*RHO(im:im+1,j,k))
                    ns = ns + 4
                 end if
              end if
           end if

           
        end do
     end do
  end do
  
  Delta_ty = Delta_ty + dt
  
  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dy_sample

subroutine stat_1dx_sample
  use stat_1d
  use data
  use combustion
  use time_info
  use memory
  use soot
  use pah
  use parallel
  use metric_generic
  implicit none
  
  integer  :: i,j,jm,k,isc,loc,ns,n
  logical  :: found_loc
  real(WP) :: xcyl,ycyl,zcyl,ucyl,vcyl,wcyl
  character(len=str_short) :: name

  ! Start the timer
  call timing_start('stat_1d')
  
  ! Prepare the chemical variable if necessary
  if (trim(chemistry).eq.'chemtable') then
     call chemtable_lookup('T'    ,tmp1)
     call chemtable_lookup('Y_F'  ,tmp2)
     call chemtable_lookup('Y_O2' ,tmp3)
     call chemtable_lookup('Y_CO' ,tmp4)
     call chemtable_lookup('Y_CO2',tmp5)
     call chemtable_lookup('Y_H2' ,tmp6)
     call chemtable_lookup('Y_H2O',tmp7)
     call chemtable_lookup('Y_OH' ,tmp8)
     if (use_soot .and. .not.use_pah) call chemtable_lookup('Y_PAH',tmp9)
  end if

  ! Prepare the density at the x-faces
  if (trim(chemistry).ne.'none') then
     do i=imin_,imax_
        do j=jmin_,jmax_
           do k=kmin_,kmax_
              tmp9(i,j,k) = sum(interp_sc_x(i,j,:)*RHO(i-st2:i+st1,j,k))
           end do
        end do
     end do
  end if
  
  ! Gather the stats
  do loc=1,nyloc
     j  = jloc (loc)
     jm = jmloc(loc)
     do k=kmin_,kmax_
        do i=imin_,imax_
           
           ns = 0
           
           ! Density
           if (trim(chemistry).ne.'none') then
              stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cymloc(loc,:)*RHO(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+3) = stat_x(loc,i,ns+3) + dt*sum(cyloc (loc,:)*tmp9(i,j:j+1,k))
              ns = ns+3
           end if

           ! Velocity statistics
           stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*U(i,jm:jm+1,k))
           stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cyloc (loc,:)*V(i,j:j+1  ,k))
           stat_x(loc,i,ns+3) = stat_x(loc,i,ns+3) + dt*sum(cymloc(loc,:)*W(i,jm:jm+1,k))
           stat_x(loc,i,ns+4) = stat_x(loc,i,ns+4) + dt*sum(cymloc(loc,:)*U(i,jm:jm+1,k)**2)
           stat_x(loc,i,ns+5) = stat_x(loc,i,ns+5) + dt*sum(cyloc (loc,:)*V(i,j:j+1  ,k)**2)
           stat_x(loc,i,ns+6) = stat_x(loc,i,ns+6) + dt*sum(cymloc(loc,:)*W(i,jm:jm+1,k)**2)
           stat_x(loc,i,ns+7) = stat_x(loc,i,ns+7) + dt*sum(cymloc(loc,:)*U(i,jm:jm+1,k)**3)
           stat_x(loc,i,ns+8) = stat_x(loc,i,ns+8) + dt*sum(cyloc (loc,:)*V(i,j:j+1  ,k)**3)
           stat_x(loc,i,ns+9) = stat_x(loc,i,ns+9) + dt*sum(cymloc(loc,:)*W(i,jm:jm+1,k)**3)
           ns = ns+9
           
           if (trim(chemistry).ne.'none') then
              stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*rhoU(i,jm:jm+1,k))
              stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cyloc (loc,:)*rhoV(i,j:j+1  ,k))
              stat_x(loc,i,ns+3) = stat_x(loc,i,ns+3) + dt*sum(cymloc(loc,:)*rhoW(i,jm:jm+1,k))
              stat_x(loc,i,ns+4) = stat_x(loc,i,ns+4) + dt*sum(cymloc(loc,:)*rhoU(i,jm:jm+1,k)*U(i,jm:jm+1,k))
              stat_x(loc,i,ns+5) = stat_x(loc,i,ns+5) + dt*sum(cyloc (loc,:)*rhoV(i,j:j+1  ,k)*V(i,j:j+1  ,k))
              stat_x(loc,i,ns+6) = stat_x(loc,i,ns+6) + dt*sum(cymloc(loc,:)*rhoW(i,jm:jm+1,k)*W(i,jm:jm+1,k))
              ns = ns+6
           end if
           
           ! Scalars
           do isc=1,nscalar
              name = SC_name(isc)
              if (name(1:2).ne.'S_') then
                 stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*SC(i,jm:jm+1,k,isc))
                 stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cymloc(loc,:)*SC(i,jm:jm+1,k,isc)**2)
                 ns = ns+2
                 if (trim(chemistry).ne.'none') then
                    stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*SC(i,jm:jm+1,k,isc)*RHO(i,jm:jm+1,k))
                    stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cymloc(loc,:)*SC(i,jm:jm+1,k,isc)**2*RHO(i,jm:jm+1,k))
                    ns = ns+2
                 end if
              end if
           end do
           
           ! Combustion
           if (trim(chemistry).eq.'chemtable') then
              stat_x(loc,i,ns+ 1) = stat_x(loc,i,ns+ 1) + dt*sum(cymloc(loc,:)*tmp1(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 2) = stat_x(loc,i,ns+ 2) + dt*sum(cymloc(loc,:)*tmp1(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 3) = stat_x(loc,i,ns+ 3) + dt*sum(cymloc(loc,:)*tmp1(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 4) = stat_x(loc,i,ns+ 4) + dt*sum(cymloc(loc,:)*tmp1(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 5) = stat_x(loc,i,ns+ 5) + dt*sum(cymloc(loc,:)*tmp2(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 6) = stat_x(loc,i,ns+ 6) + dt*sum(cymloc(loc,:)*tmp2(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 7) = stat_x(loc,i,ns+ 7) + dt*sum(cymloc(loc,:)*tmp2(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 8) = stat_x(loc,i,ns+ 8) + dt*sum(cymloc(loc,:)*tmp2(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 9) = stat_x(loc,i,ns+ 9) + dt*sum(cymloc(loc,:)*tmp3(i,jm:jm+1,k))
              stat_x(loc,i,ns+10) = stat_x(loc,i,ns+10) + dt*sum(cymloc(loc,:)*tmp3(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+11) = stat_x(loc,i,ns+11) + dt*sum(cymloc(loc,:)*tmp3(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+12) = stat_x(loc,i,ns+12) + dt*sum(cymloc(loc,:)*tmp3(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+13) = stat_x(loc,i,ns+13) + dt*sum(cymloc(loc,:)*tmp4(i,jm:jm+1,k))
              stat_x(loc,i,ns+14) = stat_x(loc,i,ns+14) + dt*sum(cymloc(loc,:)*tmp4(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+15) = stat_x(loc,i,ns+15) + dt*sum(cymloc(loc,:)*tmp4(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+16) = stat_x(loc,i,ns+16) + dt*sum(cymloc(loc,:)*tmp4(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+17) = stat_x(loc,i,ns+17) + dt*sum(cymloc(loc,:)*tmp5(i,jm:jm+1,k))
              stat_x(loc,i,ns+18) = stat_x(loc,i,ns+18) + dt*sum(cymloc(loc,:)*tmp5(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+19) = stat_x(loc,i,ns+19) + dt*sum(cymloc(loc,:)*tmp5(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+20) = stat_x(loc,i,ns+20) + dt*sum(cymloc(loc,:)*tmp5(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+21) = stat_x(loc,i,ns+21) + dt*sum(cymloc(loc,:)*tmp6(i,jm:jm+1,k))
              stat_x(loc,i,ns+22) = stat_x(loc,i,ns+22) + dt*sum(cymloc(loc,:)*tmp6(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+23) = stat_x(loc,i,ns+23) + dt*sum(cymloc(loc,:)*tmp6(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+24) = stat_x(loc,i,ns+24) + dt*sum(cymloc(loc,:)*tmp6(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+25) = stat_x(loc,i,ns+25) + dt*sum(cymloc(loc,:)*tmp7(i,jm:jm+1,k))
              stat_x(loc,i,ns+26) = stat_x(loc,i,ns+26) + dt*sum(cymloc(loc,:)*tmp7(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+27) = stat_x(loc,i,ns+27) + dt*sum(cymloc(loc,:)*tmp7(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+28) = stat_x(loc,i,ns+28) + dt*sum(cymloc(loc,:)*tmp7(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+29) = stat_x(loc,i,ns+29) + dt*sum(cymloc(loc,:)*tmp8(i,jm:jm+1,k))
              stat_x(loc,i,ns+30) = stat_x(loc,i,ns+30) + dt*sum(cymloc(loc,:)*tmp8(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
              stat_x(loc,i,ns+31) = stat_x(loc,i,ns+31) + dt*sum(cymloc(loc,:)*tmp8(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+32) = stat_x(loc,i,ns+32) + dt*sum(cymloc(loc,:)*tmp8(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
              ns = ns+32
           end if
           
           ! Soot
           if (use_soot) then
              stat_x(loc,i,ns+ 1) = stat_x(loc,i,ns+ 1) + dt*sum(cymloc(loc,:)*volfrac(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 2) = stat_x(loc,i,ns+ 2) + dt*sum(cymloc(loc,:)*volfrac(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 3) = stat_x(loc,i,ns+ 3) + dt*sum(cymloc(loc,:)*numdens(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 4) = stat_x(loc,i,ns+ 4) + dt*sum(cymloc(loc,:)*numdens(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 5) = stat_x(loc,i,ns+ 5) + dt*sum(cymloc(loc,:)*partdiam(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 6) = stat_x(loc,i,ns+ 6) + dt*sum(cymloc(loc,:)*partdiam(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 7) = stat_x(loc,i,ns+ 7) + dt*sum(cymloc(loc,:)*partaggr(i,jm:jm+1,k))
              stat_x(loc,i,ns+ 8) = stat_x(loc,i,ns+ 8) + dt*sum(cymloc(loc,:)*partaggr(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+ 9) = stat_x(loc,i,ns+ 9) + dt*sum(cymloc(loc,:)*intermit(i,jm:jm+1,k))
              stat_x(loc,i,ns+10) = stat_x(loc,i,ns+10) + dt*sum(cymloc(loc,:)*intermit(i,jm:jm+1,k)**2)
              stat_x(loc,i,ns+11) = stat_x(loc,i,ns+11) + dt*sum(cymloc(loc,:)*Nsrc_nucl(i,jm:jm+1,k))
              stat_x(loc,i,ns+12) = stat_x(loc,i,ns+12) + dt*sum(cymloc(loc,:)*Nsrc_coag(i,jm:jm+1,k))
              stat_x(loc,i,ns+13) = stat_x(loc,i,ns+13) + dt*sum(cymloc(loc,:)*Nsrc_ox  (i,jm:jm+1,k))
              stat_x(loc,i,ns+14) = stat_x(loc,i,ns+14) + dt*sum(cymloc(loc,:)*Nsrc_frag(i,jm:jm+1,k))
              stat_x(loc,i,ns+15) = stat_x(loc,i,ns+15) + dt*sum(cymloc(loc,:)*FVsrc_nucl(i,jm:jm+1,k))
              stat_x(loc,i,ns+16) = stat_x(loc,i,ns+16) + dt*sum(cymloc(loc,:)*FVsrc_cond(i,jm:jm+1,k))
              stat_x(loc,i,ns+17) = stat_x(loc,i,ns+17) + dt*sum(cymloc(loc,:)*FVsrc_sg  (i,jm:jm+1,k))
              stat_x(loc,i,ns+18) = stat_x(loc,i,ns+18) + dt*sum(cymloc(loc,:)*FVsrc_ox  (i,jm:jm+1,k))
              ns = ns+18
              if (.not.use_pah) then
                 stat_x(loc,i,ns+1) = stat_x(loc,i,ns+1) + dt*sum(cymloc(loc,:)*tmp9(i,jm:jm+1,k))
                 stat_x(loc,i,ns+2) = stat_x(loc,i,ns+2) + dt*sum(cymloc(loc,:)*tmp9(i,jm:jm+1,k)*RHO(i,jm:jm+1,k))
                 stat_x(loc,i,ns+3) = stat_x(loc,i,ns+3) + dt*sum(cymloc(loc,:)*tmp9(i,jm:jm+1,k)**2)
                 stat_x(loc,i,ns+4) = stat_x(loc,i,ns+4) + dt*sum(cymloc(loc,:)*tmp9(i,jm:jm+1,k)**2*RHO(i,jm:jm+1,k))
                 ns = ns+4
              end if
           end if
           
        end do
     end do
  end do
  
  Delta_tx = Delta_tx + dt
  
  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dx_sample


! ================================= !
! Read the statistics from the disk !
! ================================= !
subroutine stat_1dy_read
  use stat_1d
  use parallel
  implicit none
  
  real(WP) :: time
  integer, dimension(3) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  character(len=str_medium)  :: name
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,loc
  
  ! Open the file to write
  filename = trim(mpiiofs)//":stat/stat-1Dy"
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
  
  ! Read dimensions from header
  call MPI_FILE_READ_ALL(ifile,dims,3,MPI_INTEGER,status,ierr)
  if (dims(1).ne.nxloc .or. dims(2).ne.ny) then
     print*, 'expected = ',nxloc,ny
     print*, 'stat = ',dims(1),dims(2)
     call die('stat_1dy_read: The size of the stat file is incorrect')
  end if
  if (dims(3).ne.stat_nvar) call die('stat_1dy_read: Wrong number of variables in stat file')
  
  ! Read some headers
  call MPI_FILE_READ_ALL(ifile,Delta_ty,1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(ifile,time,1,MPI_REAL_WP,status,ierr)
  
  ! Read variable names
  do var=1,stat_nvar
     call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierr)
     if (name.ne.stat_name(var)) then
        print*,irank,name,stat_name(var)
        call die('stat_1dy_read: Variables names in stat and data files do not correspond')
     end if
  end do
  
  ! Read
  call MPI_FILE_READ_ALL(ifile,buf_y,nxloc*ny*stat_nvar,MPI_REAL_WP,status,ierr)
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  ! Restart the stats
  do loc=1,nxloc
     if (iloc_ifwrite(loc).eq.1) stat_y(loc,jmin_:jmax_,:) = buf_y(loc,jmin_:jmax_,:) * Delta_ty*real(nz_,WP)
  end do
  
  return
end subroutine stat_1dy_read

subroutine stat_1dx_read
  use stat_1d
  use parallel
  implicit none
  
  real(WP) :: time
  integer, dimension(3) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  character(len=str_medium)  :: name
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,loc
  
  ! Open the file to write
  filename = trim(mpiiofs)//":stat/stat-1Dx"
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
  
  ! Read dimensions from header
  call MPI_FILE_READ_ALL(ifile,dims,3,MPI_INTEGER,status,ierr)
  if ((dims(1).ne.nyloc) .or. (dims(2).ne.nx)) then
     print*, 'expected = ',nyloc,nx
     print*, 'stat = ',dims(1),dims(2)
     call die('stat_1dx_read: The size of the stat file is incorrect')
  end if
  if (dims(3).ne.stat_nvar) call die('stat_1dx_read: Wrong number of variables in stat file')
  
  ! Read some headers
  call MPI_FILE_READ_ALL(ifile,Delta_tx,1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(ifile,time,1,MPI_REAL_WP,status,ierr)
  
  ! Read variable names
  do var=1,stat_nvar
     call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierr)
     if (name.ne.stat_name(var)) then
        call die('stat_1dx_read: Variables names in stat and data files do not correspond')
     end if
  end do
  
  ! Read
  call MPI_FILE_READ_ALL(ifile,buf_x,nyloc*nx*stat_nvar,MPI_REAL_WP,status,ierr)

  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  ! Restart the stats
  do loc=1,nyloc
     if (jloc_ifwrite(loc).eq.1) stat_x(loc,imin_:imax_,:) = buf_x(loc,imin_:imax_,:) * Delta_tx*real(nz_,WP)
  end do
  
  return
end subroutine stat_1dx_read


! =============================================== !
! Read the statistics from the disk - CONDITIONAL !
! =============================================== !
subroutine stat_1dy_cond_read
  use stat_1d
  use parallel
  implicit none

  integer, dimension(3) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  character(len=str_medium) :: name
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,loc,bin

  ! Open the file to write
  filename = trim(mpiiofs)//":stat/stat-1Dy-cond"
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)

  ! Read dimension from header
  call MPI_FILE_READ_ALL(ifile,dims,3,MPI_INTEGER,status,ierr)
  if (dims(1).ne.nxloc .and. dims(2).ne.nbins_cond) then
     print*, 'expected = ',nxloc,nbins_cond
     print*, 'stat= ',dims(1),dims(2)
     call die('stat_1dy_cond_read: The size of the stat file is incorrect')
  end if
  if (dims(3).ne.stat_nvar_cond) call die('stat_1dy_cond_read: Wrong number of variables in stat file')

  ! Read variable names
  do var=1,stat_nvar_cond
     call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierr)
     if (name.ne.stat_name_cond(var)) then
        print*,irank,name,stat_name_cond(var)
        call die('stat_1dy_cond_read: Variable names in stat and data files do no correspond')
     end if
  end do

  ! Read
  call MPI_FILE_READ_ALL(ifile,buf_nSamples_cond,nxloc*nbins_cond,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(ifile,buf_y_cond,nxloc*nbins_cond*stat_nvar_cond,MPI_REAL_WP,status,ierr)

  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)

  ! Restart the stat
  do loc=1,nxloc
     if (iloc_ifwrite(loc).eq.1) then
        do bin=1,nbins_cond
           nSamples_cond(loc,bin) = buf_nSamples_cond(loc,bin)
           stat_y_cond(loc,bin,:) = buf_y_cond(loc,bin,:)*nSamples_cond(loc,bin)
        end do
     end if
  end do

  return
end subroutine stat_1dy_cond_read


! ================================ !
! Write the statistics to the disk !
! ================================ !
subroutine stat_1dy_write
  use stat_1d
  use parallel
  use time_info
  use fileio
  implicit none
  
  character(len=str_medium) :: filename
  integer  :: j,loc,iunit,ierr,var

  ! Start the timer
  call timing_start('stat_1d')
  
  ! Nothing to write
  if (Delta_ty.eq.0.0_WP) return

  ! Gather the data
  call parallel_sum(stat_y,buf_y)
  buf_y = buf_y / (Delta_ty*real(nz,WP))

  ! Return if not the root process
  if (irank.ne.iroot) then
     if (cond_stat) call stat_1dy_cond_write
     call timing_stop('stat_1d')  
     return
  end if

  ! ------- WRITE THE BINARY FILE ------
  ! Open the file
  filename = "stat/stat-1Dy"
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

  ! Write dimensions
  call BINARY_FILE_WRITE(iunit,nxloc,1,kind(nxloc),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,stat_nvar,1,kind(stat_nvar),ierr)
  call BINARY_FILE_WRITE(iunit,Delta_ty,1,kind(Delta_ty),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)

  ! Write variable names
  do var=1,stat_nvar
     call BINARY_FILE_WRITE(iunit,stat_name(var),str_medium,kind(stat_name),ierr)
  end do
  
  ! Write the stats
  call BINARY_FILE_WRITE(iunit,buf_y,ny*stat_nvar*nxloc,kind(buf_y),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! ------- WRITE THE ASCII FILE ------
  if (allx) then
     iunit = iopen()
     write(filename,'(a18)') "stat/stat-1Dy.txt"
     open (iunit, file=filename, form="formatted", iostat=ierr)
     write(iunit,'(10000a20)') 'y', 'ym', (trim(adjustl(stat_name(var))),var=1,stat_nvar)
     do j=jmin,jmax
        write(iunit,'(10000ES20.12)') y(j)+epsilon(y(j)),ym(j)+epsilon(ym(j)),sum(buf_y(:,j,:),dim=1)/real(nx,WP)
     end do
     close(iclose(iunit))
  else
     do loc=1,nxloc
        iunit = iopen()
        write(filename,'(a14,e12.6e2,a4)') "stat/stat-1Dy-",xloc(loc),".txt"
        open (iunit, file=filename, form="formatted", iostat=ierr)
        write(iunit,'(10000a20)') 'y', 'ym', (trim(adjustl(stat_name(var))),var=1,stat_nvar)
        do j=jmin,jmax
           write(iunit,'(10000ES20.12)') y(j)+epsilon(y(j)),ym(j)+epsilon(ym(j)),buf_y(loc,j,:)
        end do
        close(iclose(iunit))
     end do
  end if
  
  ! Log
  call monitor_log("1D-Y STATISTICS FILE WRITTEN")

  ! Conditional statistics
  if (cond_stat) call stat_1dy_cond_write
  
  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dy_write


subroutine stat_1dx_write
  use stat_1d
  use parallel
  use time_info
  use fileio
  implicit none
  
  character(len=str_medium) :: filename
  integer  :: i,loc,iunit,ierr,var

  ! Start the timer
  call timing_start('stat_1d')
  
  ! Nothing to write
  if (Delta_tx.eq.0.0_WP) return

  ! Gather the data
  call parallel_sum(stat_x,buf_x)
  buf_x = buf_x / (Delta_tx*real(nz,WP))

  ! Return if not the root process
  if (irank.ne.iroot) then
     call timing_stop('stat_1d')  
     return
  end if

  ! ------- WRITE THE BINARY FILE ------
  ! Open the file
  filename = "stat/stat-1Dx"
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

  ! Write dimensions
  call BINARY_FILE_WRITE(iunit,nyloc,1,kind(nyloc),ierr)
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,stat_nvar,1,kind(stat_nvar),ierr)
  call BINARY_FILE_WRITE(iunit,Delta_tx,1,kind(Delta_tx),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)

  ! Write variable names
  do var=1,stat_nvar
     call BINARY_FILE_WRITE(iunit,stat_name(var),str_medium,kind(stat_name),ierr)
  end do
  
  ! Write the stats
  call BINARY_FILE_WRITE(iunit,buf_x,nx*stat_nvar*nyloc,kind(buf_x),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! ------- WRITE THE ASCII FILE ------
  do loc=1,nyloc
     iunit = iopen()
     write(filename,'(a14,e12.6e2,a4)') "stat/stat-1Dx-",yloc(loc),".txt"
     open (iunit, file=filename, form="formatted", iostat=ierr)
     write(iunit,'(10000a20)') 'x', 'xm', (trim(adjustl(stat_name(var))),var=1,stat_nvar)
     do i=imin,imax
        write(iunit,'(10000ES20.12)') x(i)+epsilon(x(i)),xm(i)+epsilon(xm(i)),buf_x(loc,i,:)
     end do
     close(iclose(iunit))
  end do
  
  ! Log
  call monitor_log("1D-X STATISTICS FILE WRITTEN")
  
  ! Stop the timer
  call timing_stop('stat_1d')
  
  return
end subroutine stat_1dx_write


! ============================================== !
! Write the statistics to the disk - CONDITIONAL !
! ============================================== !
subroutine stat_1dy_cond_write
  use stat_1d
  use parallel
  use time_info
  use fileio
  implicit none

  character(len=str_medium) :: filename
  integer :: bin,loc,iunit,ierr,var

  ! Nothing to write
  if (maxval(nSamples_cond).lt.0.1_WP) return

  ! Gather the data
  call parallel_sum(stat_y_cond,buf_y_cond)
  call parallel_sum(nSamples_cond,buf_nSamples_cond)

  ! Only the root continues
  if (irank.ne.iroot) return

  ! Compute conditional data
  do loc=1,nxloc
     do bin=1,nbins_cond
        if (buf_nSamples_cond(loc,bin).gt.0.1_WP) then
           buf_y_cond(loc,bin,:) = buf_y_cond(loc,bin,:) / buf_nSamples_cond(loc,bin)
        else
           buf_y_cond(loc,bin,:) = 0.0_WP
        end if
     end do
  end do

  ! ------- WRITE THE BINARY FILE -------
  ! Open the file
  filename = "stat/stat-1Dy-cond"
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

  ! Write dimensions
  call BINARY_FILE_WRITE(iunit,nxloc,1,kind(nxloc),ierr)
  call BINARY_FILE_WRITE(iunit,nbins_cond,1,kind(nbins_cond),ierr)
  call BINARY_FILE_WRITE(iunit,stat_nvar_cond,1,kind(stat_nvar_cond),ierr)

  ! Write variable names
  do var=1,stat_nvar_cond
     call BINARY_FILE_WRITE(iunit,stat_name_cond(var),str_medium,kind(stat_name_cond),ierr)
  end do

  ! Write the stats
  call BINARY_FILE_WRITE(iunit,buf_nSamples_cond,nbins_cond*nxloc,kind(buf_nSamples_cond),ierr)
  call BINARY_FILE_WRITE(iunit,buf_y_cond,nbins_cond*stat_nvar_cond*nxloc,kind(buf_y_cond),ierr)

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! ------- WRITE THE ASCII FILE -------
  if (allx) then
     iunit = iopen()
     write(filename,'(a22)') "stat/stat-1Dy-cond.txt"
     open (iunit, file=filename, form="formatted", iostat=ierr)
     write(iunit,'(10000a20)') 'ZMIX', (trim(adjustl(stat_name_cond(var))),var=1,stat_nvar_cond)
     do bin =1,nbins_cond
        write(iunit,'(10000ES20.12)') (real(bin,WP)-0.5_WP)/real(nbins_cond,WP),sum(buf_y_cond(:,bin,:),dim=1)/real(nx,WP)
     end do
     close(iclose(iunit))
  else
     do loc=1,nxloc
        iunit = iopen()
        write(filename,'(a19,e12.6e2,a4)') "stat/stat-1Dy-cond-",xloc(loc),".txt"
        open (iunit, file=filename, form="formatted", iostat=ierr)
        write(iunit,'(10000a20)') 'ZMIX', (trim(adjustl(stat_name_cond(var))),var=1,stat_nvar_cond)
        do bin=1,nbins_cond
           write(iunit,'(10000ES20.12)') (real(bin,WP)-0.5_WP)/real(nbins_cond,WP),buf_y_cond(loc,bin,:)
        end do
        close(iclose(iunit))
     end do
  end if

  ! Log
  call monitor_log("1D-Y COND STATISTICS FILE WRITTEN")

  return
end subroutine stat_1dy_cond_write
