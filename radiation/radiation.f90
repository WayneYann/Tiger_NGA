module radiation
  use precision

  ! Radiation model
  logical :: use_radiation

  ! Number of ordinates
  ! Ordinates
  ! Weights

  ! Radiation intensity
  real(WP), dimension(:,:,:,:), pointer :: Intensity
  real(WP), dimension(:,:,:,:), pointer :: IntensityOLD

  contains

    ! Set optimal ordinates
    subroutine radiation_ordinates
      implicit none

!!$      if (nordinates.eq.1) then
!!$         ordinate(1,:) = 5
!!$      elseif (nordinates.eq.1) then
!!$         ordinate(1,:) = 4
!!$         ordinate(2,:) = 5
!!$      end if

      return
    end subroutine radiation_ordinates

    ! Compute weights for source term
    subroutine radiation_weights
      implicit none

      return
    end subroutine radiation_weights

end module radiation


! =============================== !
! Initialize the radiation module !
! =============================== !
subroutine radiation_init
  use radiation
  implicit none

  ! Read parser to set logical

  ! Set ordinates
  ! Set weights
  
  return
end subroutine radiation_init


! ======= !
! Prestep !
! ======= !
subroutine radation_prestep
  use radiation
  implicit none

  IntensityOLD = Intensity

  return
end subroutine radation_prestep


! ========= !
! Solve RTE !
! ========= !
subroutine radiation_step
  use radiation
  implicit none

  ! Solve the RTE

  return
end subroutine radiation_step


! ======================= !
! Compute the source term !
! ======================= !
subroutine radiation_compute_source
  use radiation
  implicit none

  ! Compute qrad for enthalpy equation

  return
end subroutine
