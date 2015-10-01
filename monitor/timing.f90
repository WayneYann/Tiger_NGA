module timing
  use precision
  use string
  use parallel
  implicit none

  ! Number of timer
  integer :: ntimers

  ! Maximum of timers to be used
  integer, parameter :: ntimers_max = 32
  
  ! Definition of a timer
  type :: timer
     character(len=str_medium) :: name
     real(WP) :: time
     real(WP) :: time_in
     logical  :: started
  end type timer
  
  ! Array of timers
  type(timer), dimension(:), pointer :: timers

  ! Number and List of timers started
  integer :: nstarted
  integer, dimension(:), pointer :: started

  ! Timing for full timestep
  real(WP) :: full_time_in

  ! Array of values to transfer to monitor
  real(WP), dimension(:), pointer :: mval

contains
  
  ! Find a timer of a given name
  ! ----------------------------
  subroutine timing_find_timer(name,itimer)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(out) :: itimer

    loop: do itimer=1,ntimers
       if (trim(timers(itimer)%name).eq.trim(name)) exit loop
    end do loop
    if (itimer.gt.ntimers) then
       print*,'Timer : ', name
       call die('timing_find_timer: unknown timer')
    end if
    
    return
  end subroutine timing_find_timer

end module timing


! ============================ !
! Initialize the timing module !
! ============================ !
subroutine timing_pre_init
  use timing
  implicit none
  
  ! Start with no timers
  ntimers  = 0
  nstarted = 0

  ! Get time for full time step
  full_time_in = MPI_WTIME()

  ! Allocate the timers
  allocate(timers (ntimers_max))
  allocate(started(ntimers_max))

  return
end subroutine timing_pre_init
  

! ============================================= !
! Get all the timers and return that to monitor !
! ============================================= !
subroutine timing_post_init
  use timing
  implicit none
  integer :: itimer
  
  ! Allocate the array of data to transfer to monitor
  allocate(mval(2*ntimers+4))
  
  ! Create the file to monitor
  call monitor_create_file_step('timing',2*ntimers+4)
  call monitor_set_header(1,'Total [s]','r')
  call monitor_set_header(2,'Time/points [us]','r')
  do itimer=1,ntimers
     call monitor_set_header(itimer*2+1,trim(timers(itimer)%name)//' [s]','r')
     call monitor_set_header(itimer*2+2,trim(timers(itimer)%name)//' [%]','r')
  end do
  call monitor_set_header(2*ntimers+3,'Rest [s]','r')
  call monitor_set_header(2*ntimers+4,'Rest [%]','r')

  return
end subroutine timing_post_init


! ========================= !
! Create a new timer object !
! ========================= !
subroutine timing_create(name)
  use timing
  implicit none
  character(len=*), intent(in) :: name

  ! Add a timer
  ntimers = ntimers+1
  if (ntimers.gt.ntimers_max) call die ('timing_create: too many timers')

  ! Preset the values
  timers(ntimers)%name = name
  timers(ntimers)%time = 0.0_WP
  timers(ntimers)%time_in = 0.0_WP
  timers(ntimers)%started = .false.
  
  return
end subroutine timing_create


! =============== !
! Start the timer !
! =============== !
subroutine timing_start(name)
  use timing
  implicit none
  character(len=*), intent(in) :: name
  integer ::itimer
  real(WP) :: time_in

  call timing_find_timer(name,itimer)

  if (timers(itimer)%started) then
     print*,'Timer : ', name
     call die('timing_start: timer already started')
  end if
  
  time_in = MPI_WTIME()

  timers(itimer)%time_in = time_in
  timers(itimer)%started = .true.

  if (nstarted.ne.0) then
     timers(started(nstarted))%time = timers(started(nstarted))%time &
          + time_in - timers(started(nstarted))%time_in
     timers(started(nstarted))%time_in = 0.0_WP
  end if
  nstarted = nstarted+1
  started(nstarted) = itimer
  
  return
end subroutine timing_start


! ============== !
! Stop the timer !
! ============== !
subroutine timing_stop(name)
  use timing
  implicit none
  character(len=*), intent(in) :: name
  integer ::itimer
  real(WP) :: time_out

  call timing_find_timer(name,itimer)

  if (.not.timers(itimer)%started) then
     print*,'Timer : ', name
     call die('timing_stop: timer already stopped')
  end if
  
  time_out = MPI_WTIME()

  timers(itimer)%time = timers(itimer)%time + time_out-timers(itimer)%time_in
  timers(itimer)%time_in = 0.0_WP
  timers(itimer)%started = .false.

  nstarted = nstarted-1  
  if (nstarted.ne.0) then
     timers(started(nstarted))%time_in = time_out
  end if
  
  return
end subroutine timing_stop


! ================== !
! Monitor the timers !
! ================== !
subroutine timing_monitor
  use timing
  use geometry
  implicit none
  integer :: itimer
  real(WP) :: full_time_out,rest_time

  ! Get time for full time step
  full_time_out = MPI_WTIME()

  ! Create the array of values
  rest_time = 0.0_WP
  mval(1) = full_time_out-full_time_in + epsilon(1.0_WP)
  mval(2) = 1e6_WP * mval(1) * real(nproc,WP) / real(nx*ny*nz,WP)
  do itimer=1,ntimers
     rest_time = rest_time+timers(itimer)%time
     mval(2*itimer+1) = timers(itimer)%time
     mval(2*itimer+2) = 100.0_WP * timers(itimer)%time / mval(1)
  end do
  mval(2*ntimers+3) = mval(1) - rest_time
  mval(2*ntimers+4) = 100.0_WP * (mval(1) - rest_time) / mval(1)
  
  ! Transfer data to monitor
  call monitor_select_file('timing')
  call monitor_set_array_values(mval)
  
  ! Reset the timers
  full_time_in = full_time_out
  do itimer=1,ntimers
     timers(itimer)%time = 0.0_WP
     timers(itimer)%time_in = 0.0_WP
  end do
  
  return
end subroutine timing_monitor

