module config
  use precision
  use string
  implicit none

  ! Cylindrical or Cartesian geometry
  integer :: icyl
  ! Sector (1) or full 360 (0)
  integer :: isect
  ! Directions of periodicity
  integer :: xper,yper,zper

  ! Name of simulation
  character(len=str_medium) :: simu_name
  ! Type of simulation
  character(len=str_medium) :: simu_type
  
  ! Do we use combustion
  logical :: combust
  ! Do we use SGS
  logical :: use_sgs
  ! Do we use levelset
  logical :: use_lvlset
  ! Dow we have spray
  logical :: use_spray
  ! Do we have pollutants
  logical :: use_nox
  logical :: use_soot
  logical :: use_soot_sgs
  logical :: use_pah
  
  ! Velocity/Pressure and Scalar scheme
  integer :: vel_conv_order
  integer :: vel_visc_order
  character(len=str_medium) :: scalar_scheme
  
end module config


! ======================================================== !
! Read in the schemes for the velocity/pressure and scalar !
! Return the number of overlap required for the schemes    !
! ======================================================== !
subroutine config_get_schemes(nover)
  use config
  use parser
  implicit none
  integer, intent(out) :: nover
  integer :: loc,nh
  character(len=str_medium) :: buffer
  
  ! Velocity scheme
  call parser_read('Velocity conv scheme',vel_conv_order,2)
  call parser_read('Velocity visc scheme',vel_visc_order,2)
  if (mod(vel_conv_order,2).ne.0 .or. mod(vel_visc_order,2).ne.0) &
       call die('config_get_schemes: velocity schemes must be even.')
  nover = max(vel_conv_order,vel_visc_order)-1
  
  ! Scalar scheme
  call parser_read('Scalar scheme',scalar_scheme,'none')
  loc = index(scalar_scheme,' ')
  if (loc.ne.0) then
     buffer=trim(scalar_scheme(loc+1:))
     scalar_scheme=trim(scalar_scheme(1:loc))
  end if
  select case (trim(scalar_scheme))
  case ('up')
     nover = max(nover,1)
  case ('quick')
     nover = max(nover,2)
  case ('bquick')
     nover = max(nover,2)
  case ('weno3')
     nover = max(nover,2)
  case ('weno5')
     nover = max(nover,3)
  case ('houc')
     read(buffer,'(i10)') nh
     nover = max(nover,(nh+1)/2)
  case ('none')
     ! No scalar has been defined
  end select
  
  return
end subroutine config_get_schemes
