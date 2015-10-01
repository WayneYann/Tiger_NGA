module spray
  use combustion
  use chemtable
  implicit none

  ! Indices
  integer :: isc_SP_M0, isc_SP_M1

  ! Derived mean quantities
  real(WP), dimension(:,:,:), pointer :: numdens
  real(WP), dimension(:,:,:), pointer :: volfrac
  real(WP), dimension(:,:,:), pointer :: dropdiam

  ! Arrays of source terms (for computing statistics)

  ! Values to monitor
  real(WP) :: min_N, max_N
  real(WP) :: min_fv, max_fv
  real(WP) :: min_smd, max_smd

  ! Liquid properties
  real(WP) :: density_liquid

  ! Local data
  real(WP) :: massrate
  real(WP) :: heatrate
  
  !$OMP THREADPRIVATE(massrate,heatrate)

contains

  ! Dispersed-phase source terms
  ! ----------------------------
  subroutine spray_source_liquid(SC_,srcSC_)
    use time_info
    implicit none
    
    real(WP), dimension(nScalar), intent(in) :: SC_
    real(WP), dimension(nScalar), intent(inout) :: srcSC_    

    ! Mass Transfer Rate
    massrate = density_liquid*SC_(isc_SP_M1)/1.0e-3_WP
    srcSC_(isc_SP_M1) = srcSC_(isc_SP_M1) - massrate*dt

    ! Heat Transfer Rate
    heatrate = 0.0_WP

    return
  end subroutine spray_source_liquid

  ! Gas-phase source terms
  ! ----------------------
  subroutine spray_source_gas(SC_,srcSC_,srcP_)
    use time_info
    implicit none

    real(WP), dimension(nScalar), intent(in) :: SC_
    real(WP), dimension(nScalar), intent(inout) :: srcSC_
    real(WP), intent(inout) :: srcP_

    ! Continuity
    srcP_ = srcP_ + massrate*dt

    ! Mixture Fraction
    ! CHECK THIS
    srcSC_(isc_ZMIX) = srcSC_(isc_ZMIX) + massrate

    ! Mixture Fraction Squared
    ! Nothing yet; need to implement saturation model
    ! WHAT ABOUT THE SECOND SOURCE TERM????
    srcSC_(isc_ZMIX2) = srcSC_(isc_ZMIX2) + massrate*SC_(isc_ZMIX)

    ! Mixture Fraction Variance
    ! Nothing yet; need to implement saturation model

    ! Enthalpy
    ! NEED TO IMPLEMENT

    ! Temperature
    ! NEED TO IMPLEMENT

    ! Species Mass Fractions

    return
  end subroutine spray_source_gas

  ! Spray momentum
  ! --------------
  subroutine spray_momentum
    implicit none

    ! For now, same as gas momentum/velocity

    return
  end subroutine spray_momentum

  ! Spray diffusivity
  ! -----------------
  subroutine spray_diffusivity
    implicit none

    integer :: i,j,k

    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             DIFF(i,j,k,isc_SP_M0) = 0.0_WP
             DIFFmol(i,j,k,isc_SP_M0) = 0.0_WP
             
             DIFF(i,j,k,isc_SP_M1) = 0.0_WP
             DIFFmol(i,j,k,isc_SP_M1) = 0.0_WP
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    return
  end subroutine spray_diffusivity

end module spray


! =========================== !
! Initialize the spray module !
! =========================== !
subroutine spray_init
  use spray
  implicit none

  integer :: isc

  ! If not using spray model, then return
  call parser_read('Use spray model',use_spray,.false.)
  if (.not.use_spray) return

  ! Create and start the time
  call timing_create('spray')
  call timing_start ('spray')

  ! Detect variable name
  isc_SP_M0 = 0
  isc_SP_M1 = 0
  do isc=1,nScalar
     select case(trim(SC_name(isc)))
     case ('SP_M0')
        isc_SP_M0 = isc
     case ('SP_M1')
        isc_SP_M1 = isc
     end select
  end do

  ! Check for missing variables
  if (isc_SP_M0.eq.0 .or. isc_SP_M1.eq.0) &
       call die('spray_init: missing spray variable')

  ! Allocate/point arrays for derived quantities
  numdens => SC(:,:,:,isc_SP_M0)
  volfrac => SC(:,:,:,isc_SP_M1)
  allocate(dropdiam(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Initialize spray diffusivity
  call spray_diffusivity

  ! Get liquid density
  call parser_read('Liquid density',density_liquid)

  ! Compute derived spray quantities
  call spray_poststep

  ! Create a new file to monitor at each iteration
  call monitor_create_file_step('spray',6)
  call monitor_set_header(1,'min_N','r')
  call monitor_set_header(2,'max_N','r')
  call monitor_set_header(3,'min_fv','r')
  call monitor_set_header(4,'max_fv','r')
  call monitor_set_header(5,'min_smd','r')
  call monitor_set_header(6,'max_smd','r')

  ! Stop the timer
  call timing_stop('spray')

  return
end subroutine spray_init


! =============================================== !
! Pre-timestep routine for spray                  !
! - Compute scalar momentum (not yet implemented) !
! - Compute scalar diffusivity                    !
! =============================================== !
subroutine spray_prestep
  use spray
  implicit none

  ! If not spray, then exit
  if (.not.use_spray) return

  ! Start the timer
  call timing_start('spray')

  ! Compute spray momentum
  call spray_momentum

  ! Compute spray diffusivity
  call spray_diffusivity

  ! Stop the timer
  call timing_stop('spray')

  return
end subroutine spray_prestep


! =============================================== !
! Compute the spray source terms for the momentum !
! =============================================== !
subroutine spray_source_momentum
  use spray
  implicit none

  integer :: i,j,k

  ! Momentum source due to mass transfer
!!$  !$OMP PARALLEL DO
!!$  do k=kmin_,kmax_
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           srcU(i,j,k) = srcU(i,j,k) + dt_uvw*0.5_WP*(Uold(i,j,k)+U(i,j,k))*sum(interp_sc_x(i,j,:)*rhodot(i-st2:i+st1,j,k))
!!$           srcV(i,j,k) = srcV(i,j,k) + dt_uvw*0.5_WP*(Vold(i,j,k)+V(i,j,k))*sum(interp_sc_y(i,j,:)*rhodot(i,j-st2:j+st1,k))
!!$           srcW(i,j,k) = srcW(i,j,k) + dt*uvw*0.5_WP*(Wold(i,j,k)+W(i,j,k))*sum(interp_sc_z(i,j,:)*rhodot(i,j,k-st2:k+st1))
!!$        end do
!!$     end do
!!$  end do
!!$  !$OMP END PARALLEL DO

  ! Momentum source due to drag
  ! NOT YET IMPLEMENTED
!!$  !$OMP PARALLEL DO
!!$  do k=kmin_,kmax_
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           srcU(i,j,k) = srcU(i,j,k) + 0.0_WP
!!$           srcV(i,j,k) = srcV(i,j,k) + 0.0_WP
!!$           srcW(i,j,k) = srcW(i,j,k) + 0.0_WP
!!$        end do
!!$     end do
!!$  end do
!!$  !$OMP END PARALLEL DO

  return
end subroutine spray_source_momentum


! ============================================== !
! Compute the spray source terms for the scalars !
! ============================================== !
subroutine spray_source_scalar(SC_,srcSC_,srcP_)
  use spray
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(in) :: SC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: srcP_

  integer :: i,j,k

  ! If not spray, then return
  if (.not.use_spray) return

  ! Start the timer
  call timing_start('spray')

  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! Dispersed-phase source terms
           call spray_source_liquid(SC_(i,j,k,:),srcSC_(i,j,k,:))

           ! Gas-phase source terms
           call spray_source_gas(SC_(i,j,k,:),srcSC_(i,j,k,:),srcP_(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Stop the timer
  call timing_stop('spray')

  return
end subroutine spray_source_scalar


! =============================== !
! Post-timestep routine for spray !
! - Compute derived quantities    !
! =============================== !
subroutine spray_poststep
  use spray
  implicit none

  integer :: i,j,k

  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_

           ! Droplet diameter (in um)
           dropdiam(i,j,k) = (SC(i,j,k,isc_SP_M1)/SC(i,j,k,isc_SP_M0))**(1.0_WP/3.0_WP)*1.0e6_WP

        end do
     end do
  end do

  return
end subroutine spray_poststep


! ======================== !
! Monitor the spray module !
! ======================== !
subroutine spray_monitor
  use spray
  use time_info
  use parallel
  implicit none

  ! If not spray, then exit
  if (.not.use_nox) return

  ! Compute min/max
  call parallel_max( numdens,max_N)
  call parallel_min(-numdens,min_N)
  min_N = -min_N

  call parallel_max( volfrac,max_fv)
  call parallel_max(-volfrac,min_fv)
  min_fv = -min_fv

  call parallel_max( dropdiam,max_smd)
  call parallel_max(-dropdiam,min_smd)
  min_smd = -min_smd

  ! Transfer values to monitor
  call monitor_select_file('spray')
  call monitor_set_single_value(1,min_N)
  call monitor_set_single_value(2,max_N)
  call monitor_set_single_value(3,min_fv)
  call monitor_set_single_value(4,max_fv)
  call monitor_set_single_value(5,min_smd)
  call monitor_set_single_value(6,max_smd)

  return
end subroutine spray_monitor
