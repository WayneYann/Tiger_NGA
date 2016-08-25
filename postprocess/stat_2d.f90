module stat_2d
  use stat
  implicit none

  ! Sampled variables
  integer :: stat_nvar
  character(len=str_medium), dimension(:), pointer :: stat_name

  ! 2D Statistics
  real(WP), dimension(:,:,:), allocatable :: stat_xy
  real(WP), dimension(:,:,:), allocatable :: buf
  real(WP) :: Delta_t
  
  ! Fileview
  integer :: fileview
  
contains
  
  subroutine stat_2d_init_names
    use data
    use combustion
    use string
    implicit none
    
    integer :: isc,ns
    character(len=str_short) :: name
    
    ! ================================================================================================ !
    ! Count the variables
    stat_nvar = 0
    ! Density 
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    ! Velocity
    stat_nvar = stat_nvar+6
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+6
    ! Scalars
    do isc=1,nscalar
       name = SC_name(isc)
       if (name(1:2).ne.'S_') then
          stat_nvar = stat_nvar+2
          if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
       end if
    end do
    ! Scalar Variance
    if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
       stat_nvar = stat_nvar+2
       if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    end if
    ! Scalar Dissipation Rate
    if (isc_ZMIX.ne.0) then
       stat_nvar = stat_nvar+2
       if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    end if
    ! Combustion
    select case (trim(chemistry))
    case ('chemtable')
       stat_nvar = stat_nvar+32
    ! SD store DIFF and HR_gp
    case ('finite chem')
       stat_nvar = stat_nvar+2+2 ! for DIFF
       stat_nvar = stat_nvar+2+2 ! for HR_gp
       stat_nvar = stat_nvar+3+3 ! for three species' chemsrc and diff
    end select
    ! Soot
    if (use_soot) then
       stat_nvar = stat_nvar + 18
       if (.not.use_pah) then
          stat_nvar = stat_nvar + 4
       end if
    end if
    ! Allocate
    allocate(stat_name(stat_nvar))
    
    ns = 0
    ! Density statistics
    if (trim(chemistry).ne.'none') then
       stat_name(ns+1) = 'RHO'
       stat_name(ns+2) = 'RHO^2'
       ns = ns+2
    end if
    ! Velocity statistics
    stat_name(ns+1) = 'U'
    stat_name(ns+2) = 'V'
    stat_name(ns+3) = 'W'
    stat_name(ns+4) = 'U^2'
    stat_name(ns+5) = 'V^2'
    stat_name(ns+6) = 'W^2'
    ns = ns+6
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

    ! Scalar Variance
    if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
       stat_name(ns+1) = 'ZVAR'
       stat_name(ns+2) = 'ZVAR^2'
       ns = ns+2
       if (trim(chemistry).ne.'none') then
          stat_name(ns+1) = 'rhoZVAR'
          stat_name(ns+2) = 'rhoZVAR^2'
          ns = ns+2
       end if
    end if

    ! Scalar Dissipation Rate
    if (isc_ZMIX.ne.0) then
       stat_name(ns+1) = 'CHI'
       stat_name(ns+2) = 'CHI^2'
       ns = ns+2
       if (trim(chemistry).ne.'none') then
          stat_name(ns+1) = 'rhoCHI'
          stat_name(ns+2) = 'rhoCHI^2'
          ns = ns+2
       end if
    end if
    
    ! Combustion
    select case (trim(chemistry))
    case ('chemtable')
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
    ! SD store DIFF and HR_gp
    case ('finite chem')
       stat_name(ns+1) = 'DIFF'
       stat_name(ns+2) = 'DIFF^2'
       stat_name(ns+3) = 'rhoDIFF'
       stat_name(ns+4) = 'rhoDIFF^2'
       stat_name(ns+5) = 'HRR'
       stat_name(ns+6) = 'HRR^2'
       stat_name(ns+7) = 'rhoHRR'
       stat_name(ns+8) = 'rhoHRR^2'
       ns = ns+8
       ! SD store three species chemsrc and diff
       stat_name(ns+1) = 'SRC_H2O'
       stat_name(ns+2) = 'SRC_CO'
       stat_name(ns+3) = 'SRC_CO2'
       stat_name(ns+4) = 'D_H2O'
       stat_name(ns+5) = 'D_CO'
       stat_name(ns+6) = 'D_CO2'
       ns = ns+6
    end select
    
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

    return
  end subroutine stat_2d_init_names
  
end module stat_2d


! ================================== !
! Initialize the 2D statistic module !
! ================================== !
subroutine stat_2d_init
  use stat_2d
  use parallel
  use parser
  implicit none
  
  integer, dimension(2) :: gsizes,lsizes,start
  integer :: ierr
  
  ! Test if we can gather stats
  if (xper.eq.1) call die('stat_2d_init: 2D statistics in x impossible (x is periodic)')
  if (yper.eq.1) call die('stat_2d_init: 2D statistics in y impossible (y is periodic)')
  
  ! Get the number of variables and names
  call stat_2d_init_names
  
  ! Allocate the storage space
  allocate(stat_xy(imin_:imax_,jmin_:jmax_,stat_nvar))
  allocate(buf(imin_:imax_,jmin_:jmax_,stat_nvar))
  stat_xy = 0.0_WP
  
  ! Generate the fileview
  gsizes(1) = nx
  lsizes(1) = nx_
  start(1)  = imin_-imin
  gsizes(2) = ny
  lsizes(2) = ny_
  start(2)  = jmin_-jmin
  call MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview,ierr)
  call MPI_TYPE_COMMIT(fileview,ierr)
  
  ! Read the stat file
  call stat_2d_read
  
  return
end subroutine stat_2d_init


! ===================== !
! Sample the statistics !
! ===================== !
subroutine stat_2d_sample
  use stat_2d
  use data
  use combustion
  use finitechem
  use time_info
  use memory
  use soot
  use pah
  implicit none

  integer :: i,j,k,isc,ns,n
  character(len=str_short) :: name

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
     
  ! Gather the stats
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_

           ns = 0

           ! Density
           if (trim(chemistry).ne.'none') then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)**2
              ns = ns+2
           end if
           
           ! Velocity
           stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*U(i,j,k)
           stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*V(i,j,k)
           stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*W(i,j,k)
           stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*U(i,j,k)**2
           stat_xy(i,j,ns+5) = stat_xy(i,j,ns+5) + dt*V(i,j,k)**2
           stat_xy(i,j,ns+6) = stat_xy(i,j,ns+6) + dt*W(i,j,k)**2
           ns = ns+6
           
           if (trim(chemistry).ne.'none') then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*rhoU(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*rhoV(i,j,k)
              stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*rhoW(i,j,k)
              stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*rhoU(i,j,k)*U(i,j,k)
              stat_xy(i,j,ns+5) = stat_xy(i,j,ns+5) + dt*rhoV(i,j,k)*V(i,j,k)
              stat_xy(i,j,ns+6) = stat_xy(i,j,ns+6) + dt*rhoW(i,j,k)*W(i,j,k)
              ns = ns+6
           end if
           
           ! Scalars
           do isc=1,nscalar
              name = SC_name(isc)
              if (name(1:2).ne.'S_') then  ! SD
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*SC(i,j,k,isc)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*SC(i,j,k,isc)**2
                 ns = ns+2
                 if (trim(chemistry).ne.'none') then
                    stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*SC(i,j,k,isc)
                    stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*SC(i,j,k,isc)**2
                    ns = ns+2
                 end if
              end if
           end do

           ! Scalar Variance
           if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*ZVAR(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*ZVAR(i,j,k)**2
              ns = ns+2
              if (trim(chemistry).ne.'none') then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*ZVAR(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*ZVAR(i,j,k)**2
                 ns = ns+2
              end if
           end if

           ! Scalar Dissipation Rate
           if (isc_ZMIX.ne.0) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*CHI(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*CHI(i,j,k)**2
              ns = ns+2
              if (trim(chemistry).ne.'none') then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*CHI(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*CHI(i,j,k)**2
                 ns = ns+2
              end if
           end if
           
           ! Combustion
           select case (trim(chemistry))
           case ('chemtable')
              stat_xy(i,j,ns+ 1) = stat_xy(i,j,ns+ 1) + dt*tmp1(i,j,k)
              stat_xy(i,j,ns+ 2) = stat_xy(i,j,ns+ 2) + dt*RHO(i,j,k)*tmp1(i,j,k)
              stat_xy(i,j,ns+ 3) = stat_xy(i,j,ns+ 3) + dt*tmp1(i,j,k)**2
              stat_xy(i,j,ns+ 4) = stat_xy(i,j,ns+ 4) + dt*RHO(i,j,k)*tmp1(i,j,k)**2
              stat_xy(i,j,ns+ 5) = stat_xy(i,j,ns+ 5) + dt*tmp2(i,j,k)
              stat_xy(i,j,ns+ 6) = stat_xy(i,j,ns+ 6) + dt*RHO(i,j,k)*tmp2(i,j,k)
              stat_xy(i,j,ns+ 7) = stat_xy(i,j,ns+ 7) + dt*tmp2(i,j,k)**2
              stat_xy(i,j,ns+ 8) = stat_xy(i,j,ns+ 8) + dt*RHO(i,j,k)*tmp2(i,j,k)**2
              stat_xy(i,j,ns+ 9) = stat_xy(i,j,ns+ 9) + dt*tmp3(i,j,k)
              stat_xy(i,j,ns+10) = stat_xy(i,j,ns+10) + dt*RHO(i,j,k)*tmp3(i,j,k)
              stat_xy(i,j,ns+11) = stat_xy(i,j,ns+11) + dt*tmp3(i,j,k)**2
              stat_xy(i,j,ns+12) = stat_xy(i,j,ns+12) + dt*RHO(i,j,k)*tmp3(i,j,k)**2
              stat_xy(i,j,ns+13) = stat_xy(i,j,ns+13) + dt*tmp4(i,j,k)
              stat_xy(i,j,ns+14) = stat_xy(i,j,ns+14) + dt*RHO(i,j,k)*tmp4(i,j,k)
              stat_xy(i,j,ns+15) = stat_xy(i,j,ns+15) + dt*tmp4(i,j,k)**2
              stat_xy(i,j,ns+16) = stat_xy(i,j,ns+16) + dt*RHO(i,j,k)*tmp4(i,j,k)**2
              stat_xy(i,j,ns+17) = stat_xy(i,j,ns+17) + dt*tmp5(i,j,k)
              stat_xy(i,j,ns+18) = stat_xy(i,j,ns+18) + dt*RHO(i,j,k)*tmp5(i,j,k)
              stat_xy(i,j,ns+19) = stat_xy(i,j,ns+19) + dt*tmp5(i,j,k)**2
              stat_xy(i,j,ns+20) = stat_xy(i,j,ns+20) + dt*RHO(i,j,k)*tmp5(i,j,k)**2
              stat_xy(i,j,ns+21) = stat_xy(i,j,ns+21) + dt*tmp6(i,j,k)
              stat_xy(i,j,ns+22) = stat_xy(i,j,ns+22) + dt*RHO(i,j,k)*tmp6(i,j,k)
              stat_xy(i,j,ns+23) = stat_xy(i,j,ns+23) + dt*tmp6(i,j,k)**2
              stat_xy(i,j,ns+24) = stat_xy(i,j,ns+24) + dt*RHO(i,j,k)*tmp6(i,j,k)**2
              stat_xy(i,j,ns+25) = stat_xy(i,j,ns+25) + dt*tmp7(i,j,k)
              stat_xy(i,j,ns+26) = stat_xy(i,j,ns+26) + dt*RHO(i,j,k)*tmp7(i,j,k)
              stat_xy(i,j,ns+27) = stat_xy(i,j,ns+27) + dt*tmp7(i,j,k)**2
              stat_xy(i,j,ns+28) = stat_xy(i,j,ns+28) + dt*RHO(i,j,k)*tmp7(i,j,k)**2
              stat_xy(i,j,ns+29) = stat_xy(i,j,ns+29) + dt*tmp8(i,j,k)
              stat_xy(i,j,ns+30) = stat_xy(i,j,ns+30) + dt*RHO(i,j,k)*tmp8(i,j,k)
              stat_xy(i,j,ns+31) = stat_xy(i,j,ns+31) + dt*tmp8(i,j,k)**2
              stat_xy(i,j,ns+32) = stat_xy(i,j,ns+32) + dt*RHO(i,j,k)*tmp8(i,j,k)**2
              ns = ns+32
           ! SD store DIFF and HR_gp
           case ('finite chem')
              stat_xy(i,j,ns+ 1) = stat_xy(i,j,ns+ 1) + dt*DIFF(i,j,k,isc_ZMIX)
              stat_xy(i,j,ns+ 2) = stat_xy(i,j,ns+ 2) + dt*DIFF(i,j,k,isc_ZMIX)**2
              stat_xy(i,j,ns+ 3) = stat_xy(i,j,ns+ 3) + dt*RHO(i,j,k)*DIFF(i,j,k,isc_ZMIX)
              stat_xy(i,j,ns+ 4) = stat_xy(i,j,ns+ 4) + dt*RHO(i,j,k)*DIFF(i,j,k,isc_ZMIX)**2
              stat_xy(i,j,ns+ 5) = stat_xy(i,j,ns+ 5) + dt*HR_gp(i,j,k)
              stat_xy(i,j,ns+ 6) = stat_xy(i,j,ns+ 6) + dt*HR_gp(i,j,k)**2
              stat_xy(i,j,ns+ 7) = stat_xy(i,j,ns+ 7) + dt*RHO(i,j,k)*HR_gp(i,j,k)
              stat_xy(i,j,ns+ 8) = stat_xy(i,j,ns+ 8) + dt*RHO(i,j,k)*HR_gp(i,j,k)**2
              ns = ns+8
              ! SD store three species' chemsrc and diff
              stat_xy(i,j,ns+ 1) = stat_xy(i,j,ns+ 1) + dt*chemsrc_gp(i,j,k,1)
              stat_xy(i,j,ns+ 2) = stat_xy(i,j,ns+ 2) + dt*chemsrc_gp(i,j,k,2)
              stat_xy(i,j,ns+ 3) = stat_xy(i,j,ns+ 3) + dt*chemsrc_gp(i,j,k,3)
              stat_xy(i,j,ns+ 4) = stat_xy(i,j,ns+ 4) + dt*diff_gp(i,j,k,1)
              stat_xy(i,j,ns+ 5) = stat_xy(i,j,ns+ 5) + dt*diff_gp(i,j,k,2)
              stat_xy(i,j,ns+ 6) = stat_xy(i,j,ns+ 6) + dt*diff_gp(i,j,k,3)
              ns = ns+6
           end select
           
           ! Soot
           if (use_soot) then
              stat_xy(i,j,ns+ 1) = stat_xy(i,j,ns+ 1) + dt*volfrac(i,j,k)
              stat_xy(i,j,ns+ 2) = stat_xy(i,j,ns+ 2) + dt*volfrac(i,j,k)**2
              stat_xy(i,j,ns+ 3) = stat_xy(i,j,ns+ 3) + dt*numdens(i,j,k)
              stat_xy(i,j,ns+ 4) = stat_xy(i,j,ns+ 4) + dt*numdens(i,j,k)**2
              stat_xy(i,j,ns+ 5) = stat_xy(i,j,ns+ 5) + dt*partdiam(i,j,k)
              stat_xy(i,j,ns+ 6) = stat_xy(i,j,ns+ 6) + dt*partdiam(i,j,k)**2
              stat_xy(i,j,ns+ 7) = stat_xy(i,j,ns+ 7) + dt*partaggr(i,j,k)
              stat_xy(i,j,ns+ 8) = stat_xy(i,j,ns+ 8) + dt*partaggr(i,j,k)**2
              stat_xy(i,j,ns+ 9) = stat_xy(i,j,ns+ 9) + dt*intermit(i,j,k)
              stat_xy(i,j,ns+10) = stat_xy(i,j,ns+10) + dt*intermit(i,j,k)**2
              stat_xy(i,j,ns+11) = stat_xy(i,j,ns+11) + dt*Nsrc_nucl(i,j,k)
              stat_xy(i,j,ns+12) = stat_xy(i,j,ns+12) + dt*Nsrc_coag(i,j,k)
              stat_xy(i,j,ns+13) = stat_xy(i,j,ns+13) + dt*Nsrc_ox  (i,j,k)
              stat_xy(i,j,ns+14) = stat_xy(i,j,ns+14) + dt*Nsrc_frag(i,j,k)
              stat_xy(i,j,ns+15) = stat_xy(i,j,ns+15) + dt*FVsrc_nucl(i,j,k)
              stat_xy(i,j,ns+16) = stat_xy(i,j,ns+16) + dt*FVsrc_cond(i,j,k)
              stat_xy(i,j,ns+17) = stat_xy(i,j,ns+17) + dt*FVsrc_sg  (i,j,k)
              stat_xy(i,j,ns+18) = stat_xy(i,j,ns+18) + dt*FVsrc_ox  (i,j,k)
              ns = ns+18
              if (.not.use_pah) then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp9(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*tmp9(i,j,k)
                 stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*tmp9(i,j,k)**2
                 stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*RHO(i,j,k)*tmp9(i,j,k)**2
                 ns = ns+4
              end if
           end if
           
        end do
     end do
  end do
  !print*,'FVsrc:',maxval(FVsrc_nucl),maxval(FVsrc_cond),maxval(FVsrc_sg),maxval(FVsrc_ox)
  Delta_t = Delta_t+dt
  
  return
end subroutine stat_2d_sample


! ================================= !
! Read the statistics from the disk !
! ================================= !
subroutine stat_2d_read
  use stat_2d
  use parallel
  implicit none
  
  real(WP) :: time
  integer, dimension(4) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  character(len=str_medium)  :: name
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,data_size
  logical :: file_is_there
  
  ! -- GAS PHASE ---------------------------------------------------------------------
  inquire(file='stat/stat-2D',exist=file_is_there)
  if (file_is_there) then
     
     ! Open the file to write
     filename = trim(mpiiofs) // ":stat/stat-2D"
     call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,ifile,ierr)
     
     ! Read dimensions from header
     call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
     if ((dims(1).ne.nx) .or. (dims(2).ne.ny) .or. (dims(3).ne.1)) then
        print*, 'expected = ',nx,ny,1
        print*, 'stat = ',dims(1),dims(2),dims(3)
        call die('stat_2d_read: The size of the stat file is incorrect')
     end if
     if (dims(4).ne.stat_nvar) call die('stat_2d_read: Wrong number of variables in stat file')
     
     ! Read some headers
     call MPI_FILE_READ_ALL(ifile,Delta_t,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_READ_ALL(ifile,time,1,MPI_REAL_WP,status,ierr)
     
     ! Read variable names
     do var=1,stat_nvar
        call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierr)
        if (name.ne.stat_name(var)) then
           call die('stat_2d_read: Variables names in stat are incorrect')
        end if
     end do
          
     ! Read each variables
     data_size = nx_*ny_
     do var=1,stat_nvar
        disp = 4*4 + str_medium*stat_nvar + 2*WP + real(var-1,WP)*nx*ny*WP
        call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(ifile,buf(:,:,var),data_size,MPI_REAL_WP,status,ierr)
     end do
     
     ! Close the file
     call MPI_FILE_CLOSE(ifile,ierr)
     
     ! Recompute the stats
     stat_xy = buf*Delta_t*real(nz_)
     
  else
     
     ! Start from scratch
     Delta_t = 0.0_WP
     
  end if
  
  return
end subroutine stat_2d_read


! ================================ !
! Write the statistics to the disk !
! ================================ !
subroutine stat_2d_write
  use stat_2d
  use parallel
  use time_info
  implicit none
  
  integer, dimension(4) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,data_size,i,j
  
  ! -- GAS PHASE ---------------------------------------------------------------------
  
  ! Gather the data
  call parallel_sum_dir(stat_xy,buf,'z')
  buf = buf / (Delta_t*real(nz))
  
  ! Open the file to write
  filename = trim(mpiiofs) // ":stat/stat-2D"
  if (irank.eq.iroot) call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,ifile,ierr)
  
  ! Write the headers
  if (irank.eq.iroot) then
     ! Write dimensions
     dims(1) = nx
     dims(2) = ny
     dims(3) = 1
     dims(4) = stat_nvar
     call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
     call MPI_FILE_WRITE(ifile,Delta_t,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     ! Write variable names
     do var=1,stat_nvar
        call MPI_FILE_WRITE(ifile,stat_name(var),str_medium,MPI_CHARACTER,status,ierr)
     end do
  end if
  
  ! Write each variables
  data_size = nx_*ny_
  if (kproc.ne.1) data_size = 0
  do var=1,stat_nvar
     disp = 4*4 + str_medium*stat_nvar + 2*WP + real(var-1,WP)*nx*ny*WP
     call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(ifile,buf(:,:,var),data_size,MPI_REAL_WP,status,ierr)
  end do
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  return
end subroutine stat_2d_write
