module time_info
  use precision
  implicit none
  
  ! Time description --------------------------------- !
  ! time : scalar/density time                         !
  ! dt : time between 2 scalar fields                  !
  ! dt_uvw : time between 2 velocity fields            !
  ! -------------------------------------------------- !
  !                       RHOmid                       !
  !                       SCmid                        !
  !             Uold        U                          !
  !               |         |                          !
  !               | dt_uvw  |                          !
  !               |<------->|                          !
  !               |         |                          !
  !    n-1        n        n+1       n+2               !
  ! ----|----|----|----|----|----|----|----|----> time !
  !        n-1/2     n+1/2     n+3/2                   !
  !                    |         |                     !
  !                    |<------->|                     !
  !                    |    dt   |                     !
  !                    |         |                     !
  !                  Umid                              !
  !                 RHOold      RHO                    !
  !                 SCold        SC                    !
  !                 P , DP                             !
  ! -------------------------------------------------- !
  
  
  ! Current iteration/subiteration/time/wall-time/dt/CFL
  integer  :: ntime
  integer  :: niter
  real(WP) :: time
  real(WP) :: wtime
  real(WP) :: dt, dt_
  real(WP) :: CFL
  
  ! Maximum iteration/subiteration/time/wall-time/dt/CFL
  integer  :: max_ntime
  integer  :: max_niter
  real(WP) :: max_time
  real(WP) :: max_wtime
  real(WP) :: max_dt
  real(WP) :: max_CFL
  
  ! Additional time info
  real(WP) :: dt_old
  real(WP) :: dt_uvw
  real(WP) :: init_wtime
  
contains
  
  ! ------------------------------------------- !
  ! Check whether the simulation is done or not !
  ! ------------------------------------------- !
  function done()
    use parallel
    use fileio
    use string
    implicit none
    
    logical :: done,file_is_there
    integer :: ifile,ierr,size_int,i
    character(len=str_short) :: command
    integer(kind=MPI_Offset_kind) :: size_off
    integer, dimension(MPI_STATUS_SIZE) :: status
    character(len=str_medium) :: filename
    
    done = .false.
    
    ! Normal ending tests
    if (ntime.ge.max_ntime) then
       done = .true.
       call monitor_log("MAXIMUM NUMBER OF ITERATION REACHED")
    end if
    if ( time.ge.max_time) then
       done = .true.
       call monitor_log("MAXIMUM SIMULATION TIME REACHED")
    end if
    call parallel_time(wtime)
    wtime = wtime - init_wtime
    if ((wtime+600.0_WP).ge.max_wtime .and. max_wtime.gt.0.0_WP) then
       done = .true.
       call monitor_log("MAXIMUM WALL TIME REACHED")
    end if

    ! Check the control file
    inquire(file='control',exist=file_is_there)
    if (file_is_there) then
       
       ! Open the file
       filename = trim(mpiiofs)//":control"
       call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
       
       ! Read in the command
       call MPI_FILE_GET_SIZE(ifile,size_off,ierr)
       size_int = size_off
       call MPI_FILE_READ_ALL(ifile,command,size_int,MPI_CHARACTER,status,ierr)
       
       ! Close and remove the file
       call MPI_FILE_CLOSE(ifile,ierr)
       if (irank.eq.iroot) call MPI_FILE_DELETE(filename,mpi_info,ierr)
       
       ! Keep only letters
       do i=1,size_int
          if ( ichar(command(i:i)).lt.ichar('A') .or. &
               ichar(command(i:i)).gt.ichar('Z')) then 
             command(i:str_short) = ''
          end if
       end do
       
       ! Obey the command if recognized
       select case (trim(adjustl(command)))
       case ('KILL')
          call monitor_log("RECEIVED KILL SIGNAL BY USER")
          call die("Killed by user")
       case ('SAVE')
          call simulation_write(.true.)
       case ('END')
          call monitor_log("RECEIVED END SIGNAL BY USER")
          done = .true.
       case ('IGNITE')
          call monitor_log("RECEIVED IGNITE SIGNAL BY USER")
          call combustion_ignite
       end select
    end if
    
    return
  end function done
  
  ! ------------------------------------ !
  ! Predict the new time step size based !
  ! on dt, max_dt, CFL and max_CFL       !
  ! ------------------------------------ !
  subroutine predict_timestepsize
    implicit none
    real(WP), parameter :: alpha = 0.7_WP
    
    dt_old = dt
    dt = min(max_CFL/CFL*dt_old,max_dt)
    if (dt.gt.dt_old) dt = alpha*dt + (1.0_WP-alpha)*dt_old
    dt_uvw = 0.5_WP*(dt+dt_old)
    
    return
  end subroutine predict_timestepsize
  
end module time_info


! ================================== !
! Initialization of time information !
! ================================== !
subroutine time_info_init
  use parser
  use time_info
  implicit none
  
  ! Read in the time parameters
  call parser_read('Timestep size',max_dt)
  call parser_read('CFL number',max_CFL)
  
  ! Setup initial time step size
  ntime = 0
  if (time.eq.0.0_WP) then
     dt = max_dt
     !call velocity_CFL
     call predict_timestepsize
  end if
  dt_uvw = dt
  
  ! Read in the time limits
  call parser_read('Subiterations',max_niter)
  call parser_read('Maximum time',max_time)
  call parser_read('Maximum iterations',max_ntime)
  call parser_read('Maximum wall time',max_wtime)
  max_wtime = max_wtime * 3600.0_WP
  call parallel_time(init_wtime)
  
  return
end subroutine time_info_init
