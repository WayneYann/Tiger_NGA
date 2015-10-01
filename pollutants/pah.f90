module pah
  use combustion
  use chemtable
  implicit none

  ! Index
  integer :: isc_PAH

  ! Values to monitor
  real(WP) :: min_pah, max_pah

  ! Arrays of derived quantites (for computing statistics)
  real(WP), dimension(:,:,:), pointer :: PAHsteady
  real(WP), dimension(:,:,:), pointer :: PAHsrc_pos
  real(WP), dimension(:,:,:), pointer :: PAHsrc_neg

contains

  ! PAH momentum
  ! ------------
  subroutine pah_momentum
    implicit none

    ! Nothing to do

    return
  end subroutine pah_momentum

  ! PAH diffusivity
  ! ---------------
  subroutine pah_diffusivity
    implicit none

    integer :: i,j,k

    ! Same as mixture fraction
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             DIFF(i,j,k,isc_PAH) = DIFF(i,j,k,isc_ZMIX)
             DIFFmol(i,j,k,isc_PAH) = DIFFmol(i,j,k,isc_ZMIX)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    return
  end subroutine pah_diffusivity

end module pah


! ========================= !
! Initialize the PAH module !
! ========================= !
subroutine pah_init
  use pah
  implicit none

  integer :: isc

  ! If not using PAH  model, then return
  call parser_read('Use PAH model',use_pah,.false.)

  if (.not.use_pah) return

  ! PAH model only for chemtable
  if (trim(chemistry).ne.'chemtable') &
       call die('pah_init: Only use PAH transport equation model with chemtable')

  ! Create and start the timer
  call timing_create('pah')
  call timing_start('pah')

  ! Allocate arrays for derived quantities
  allocate(PAHsteady(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(PAHsrc_pos(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(PAHsrc_neg(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Detect variable name
  isc_PAH = -1
  do isc=1,nScalar
     select case(trim(SC_name(isc)))
     case ('PAH')
        isc_PAH = isc
     end select
  end do
  if (isc_PAH.eq.-1) call die('pah_init: missing pah variable')

  ! Initialize PAH diffusivity
  call pah_diffusivity

  ! Compute derived PAH quantities
  call pah_poststep

  ! Create a new file to monitor at each iteration
  call monitor_create_file_step('pah',2)
  call monitor_set_header(1,'min_PAH','r')
  call monitor_set_header(2,'max_PAH','r')

  ! Stop the timer
  call timing_stop('pah')

  return
end subroutine pah_init


! =============================================== !
! Pre-timestep routine for PAH                    !
! - Compute scalar momentum (not yet implemented) !
! - Compute scalar diffusivity                    !
! =============================================== !
subroutine pah_prestep
  use pah
  implicit none

  ! If not PAH, then exit
  if (.not.use_pah) return

  ! Start the timer
  call timing_start('pah')

!!$  ! Compute PAH momentum
!!$  call pah_momentum

  ! Compute PAH diffusivity
  call pah_diffusivity

  ! Stop the timer
  call timing_stop('pah')

  return
end subroutine pah_prestep


! ============================ !
! Compute the PAH source terms !
! ============================ !
subroutine pah_source_scalar(SC_,srcSC_)
  use pah
  use memory
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(in) :: SC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC_

  integer :: i,j,k

  ! If not PAH, then return
  if (.not.use_pah) return

  ! Start the timer
  call timing_start('pah')

  ! Chemical production and consumption
  call chemtable_lookup('Y_PAH',tmp1)
  call chemtable_lookup('SRC_PAH_POS',tmp2)
  call chemtable_lookup('SRC_PAH_NEG',tmp3)

  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! Production: Not rescaled
           srcSC_(i,j,k,isc_PAH) = srcSC_(i,j,k,isc_PAH) + dt_ * max(tmp2(i,j,k),0.0_WP)
           ! Consumption: Rescaled with transported value
           srcSC_(i,j,k,isc_PAH) = srcSC_(i,j,k,isc_PAH) + dt_ * min(tmp3(i,j,k),0.0_WP) * min(SC_(i,j,k,isc_PAH)/max(tmp1(i,j,k),1.0e-60_WP),1.0e2_WP) ! M coded it to be 1.0e2_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Stop the timer
  call timing_stop('pah')

  return
end subroutine pah_source_scalar


! ============================= !
! Post-timestep routine for PAH !
! - Compute derived quantities  !
! ============================= !
subroutine pah_poststep
  use pah
  implicit none

  integer :: i,j,k

  ! If not PAH, then return
  if (.not.use_pah) return

  ! Clip PAH
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           SC(i,j,k,isc_PAH) = max(SC(i,j,k,isc_PAH),1.0e-60_WP)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Save derived quantities for post-processing
  call chemtable_lookup('Y_PAH',PAHsteady)
  call chemtable_lookup('SRC_PAH_POS',PAHsrc_pos)
  call chemtable_lookup('SRC_PAH_NEG',PAHsrc_neg)  

  return
end subroutine pah_poststep


! ====================== !
! Monitor the PAH module !
! ====================== !
subroutine pah_monitor
  use pah
  use time_info
  use parallel
  implicit none

  ! If not PAH, then exit
  if (.not.use_pah) return

  call pah_poststep

  ! Compute min/max
  call parallel_max( SC(:,:,:,isc_PAH),max_pah)
  call parallel_max(-SC(:,:,:,isc_PAH),min_pah)
  min_pah = -min_pah

  ! Transfer values to monitor
  call monitor_select_file('pah')
  call monitor_set_single_value(1,min_pah)
  call monitor_set_single_value(2,max_pah)

  return
end subroutine pah_monitor
