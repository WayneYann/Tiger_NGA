module nox
  use combustion
  use chemtable
  implicit none

  ! Index
  integer :: isc_NOX

  ! Values to monitor
  real(WP) :: min_nox, max_nox

  ! Arrays of derived quantities (for computing statistics)
  real(WP), dimension(:,:,:), pointer :: NOXsteady
  real(WP), dimension(:,:,:), pointer :: NOXsrc_pos
  real(WP), dimension(:,:,:), pointer :: NOXsrc_neg

contains

  ! NOX momentum
  ! ------------
  subroutine nox_momentum
    implicit none

    ! Nothing to do

    return
  end subroutine nox_momentum

  ! NOX diffusivity
  ! ---------------
  subroutine nox_diffusivity
    implicit none

    integer :: i,j,k

    ! Same as mixture fraction
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             DIFF(i,j,k,isc_NOX) = DIFF(i,j,k,isc_ZMIX)
             DIFFmol(i,j,k,isc_NOX) = DIFFmol(i,j,k,isc_ZMIX)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    return
  end subroutine nox_diffusivity

end module nox


! ========================= !
! Initialize the NOX module !
! ========================= !
subroutine nox_init
  use nox
  implicit none

  integer :: isc

  ! If not using nox model, then return
  call parser_read('Use nox model',use_nox,.false.)

  if (.not.use_nox) return

  ! NOx model only for chemtable
  if (trim(chemistry).ne.'chemtable') &
       call die ('nox_init: Only use NOX transport equation model with chemtable')

  ! Create and start the timer
  call timing_create('nox')
  call timing_start ('nox')

  ! Allocate arrays for source terms
  allocate(NOXsteady(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(NOXsrc_pos(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(NOXsrc_neg(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Detect variable name
  isc_NOX = 0
  do isc=1,nScalar
     select case(trim(SC_name(isc)))
     case ('NOX')
        isc_NOX = isc
     end select
  end do

  ! Initialize NOX diffusivity
  call nox_diffusivity

  ! Compute derived NOX quantities
  call nox_poststep

  ! Create a new file to monitor at each iteration
  call monitor_create_file_step('nox',2)
  call monitor_set_header(1,'min_NOX','r')
  call monitor_set_header(2,'max_NOX','r')

  ! Stop the timer
  call timing_stop('nox')

  return
end subroutine nox_init


! =============================================== !
! Pre-timestep routine for NOX                    !
! - Compute scalar momentum (not yet implemented) !
! =============================================== !
subroutine nox_prestep
  use nox
  implicit none

  ! If not NOX, then exit
  if (.not.use_nox) return

  ! Start the timer
  call timing_start('nox')

!!$  ! Compute NOX momentum
!!$  call nox_momentum

  ! Compute NOX diffusivity
  call nox_diffusivity

  ! Stop the timer
  call timing_stop('nox')

  return
end subroutine nox_prestep


! ============================ !
! Compute the NOX source terms !
! ============================ !
subroutine nox_source_scalar(SC_,srcSC_)
  use nox
  use memory
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(in) :: SC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC_

  integer :: i,j,k

  ! If not NOX, then return
  if (.not.use_nox) return

  ! Start the timer
  call timing_start('nox')

  ! Chemical production and consumption
  call chemtable_lookup('Y_NOX',tmp1)
  call chemtable_lookup('SRC_NOX_POS',tmp2)
  call chemtable_lookup('SRC_NOX_NEG',tmp3)

  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! Rescale consumption with transported value
           srcSC_(i,j,k,isc_NOX) = tmp2(i,j,k) - tmp3(i,j,k) * (SC_(i,j,k,isc_NOX)/tmp1(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Stop the timer
  call timing_stop('nox')

  return
end subroutine nox_source_scalar


! ============================= !
! Post-timestep routine for NOX !
! - Compute derived quantities  !
! ============================= !
subroutine nox_poststep
  use nox
  implicit none

  ! If not NOX, then return
  if (.not.use_nox) return

  ! Start the timer
  call timing_start('nox')

  ! Save derived for post-processing
  call chemtable_lookup('Y_NOX',NOXsteady)
  call chemtable_lookup('SRC_NOX_POS',NOXsrc_pos)
  call chemtable_lookup('SRC_NOX_NEG',NOXsrc_neg)

  ! Stop the timer
  call timing_stop('nox')

  return
end subroutine nox_poststep


! ====================== !
! Monitor the NOX module !
! ====================== !
subroutine nox_monitor
  use nox
  use time_info
  use parallel
  implicit none

  ! If not NOX, then exit
  if (.not.use_nox) return

  ! Compute min/max
  call parallel_max( SC(:,:,:,isc_NOX),max_nox)
  call parallel_max(-SC(:,:,:,isc_NOX),min_nox)
  min_nox = -min_nox

  ! Transfer values to monitor
  call monitor_select_file('nox')
  call monitor_set_single_value(1,min_nox)
  call monitor_set_single_value(2,max_nox)

  return
end subroutine nox_monitor
