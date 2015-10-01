! =================================================================== !
! Pollutants                                                          !
! - Unsteady evolution of NOX, PAH, SOOT, and other future quantities !
! - I need to come up with a way to interchangably switch between     !       
!       this and computed chemistry for NOX, PAH, etc.                !
! =================================================================== !
module pollutants
  use combustion
  use chemtable
  implicit none

end module pollutants


! ================================ !
! Initialize the pollutants module !
! ================================ !
subroutine pollutants_init
  use pollutants
  implicit none

  ! Initialize the NOX module
  call nox_init

  ! Initialize the PAH module
  call pah_init

  ! Initialize the soot module
  call soot_init

  return
end subroutine pollutants_init


! =============================================== !
! Pre-timestep routine for pollutants             !
! - Compute scalar momentum (not yet implemented) !
! - Compute scalar diffusivity                    !
! =============================================== !
subroutine pollutants_prestep
  use pollutants
  implicit none

  ! Pre-timestep routine for NOX
  call nox_prestep

  ! Pre-timestep routine for PAH
  call pah_prestep

  ! Pre-timestep routine for soot
  call soot_prestep  

  return
end subroutine pollutants_prestep


! =================================================== !
! Compute the pollutant source terms for the momentum !
! =================================================== !
subroutine pollutants_source_momentum
  use pollutants
  implicit none

  ! Nothing for NOX or PAH

  ! Soot source terms
  call soot_source_momentum

  return
end subroutine pollutants_source_momentum


! ================================================== !
! Compute the pollutant source terms for the scalars !
! ================================================== !
subroutine pollutants_source_scalar(SC_,srcSC_,srcP_)
  use pollutants
  implicit none

  real(WP), dimension(nxo_,nyo_,nzo_,nScalar), intent(in) :: SC_
  real(WP), dimension(nxo_,nyo_,nzo_,nScalar), intent(inout) :: srcSC_
  real(WP), dimension(nxo_,nyo_,nzo_), intent(inout) :: srcP_

  ! NOX source terms
  call nox_source_scalar(SC_,srcSC_)

  ! PAH source terms
  call pah_source_scalar(SC_,srcSC_)

  ! Soot source terms
  call soot_source_scalar(SC_,srcSC_,srcP_)

  return
end subroutine pollutants_source_scalar


! ==================================== !
! Post-timestep routine for pollutants !
! - Compute derived quantities         !
! ==================================== !
subroutine pollutants_poststep
  use pollutants
  implicit none

  ! Post-timestep routine for NOX
  call nox_poststep
  
  ! Post-timestep routine for PAH
  call pah_poststep

  ! Post-timestep routine for soot
  call soot_poststep

  return
end subroutine pollutants_poststep


! ============================= !
! Monitor the pollutants module !
! ============================= !
subroutine pollutants_monitor
  use pollutants
  implicit none

  ! Monitor the NOX module
  call nox_monitor
  
  ! Monitor the PAH module
  call pah_monitor
  
  ! Monitor the soot module
  call soot_monitor

  return
end subroutine pollutants_monitor
