module monitor
  use precision
  use string
  implicit none

  ! Log file
  integer :: file_log

  ! Maximum number of files to monitor
  integer, parameter :: nfiles_max = 32

  ! Preset some length for the column length
  ! -> Present inside the code at 2 different locations
  integer, parameter :: col_len = 14
  character(len=*), parameter :: f1 = '(a12)'
  character(len=*), parameter :: f2 = '(i12)'
  character(len=*), parameter :: f3 = '(ES12.5)'
  ! -> Parameters for double-precision output (in monitor_dump_values_timestep ONLY)
  integer, parameter :: col_len_d = 22
  character(len=*), parameter :: f4 = '(a20)'
  character(len=*), parameter :: f5 = '(ES20.13)'

  ! Number of files to monitor
  integer :: nfiles

  ! Pointer to current file
  integer :: ifile

  ! Definition of type 'mfile'
  type mfile
     integer :: iunit
     character(len=str_medium) :: filename
     integer :: freq
     integer :: ncols
     character(len=str_medium), dimension(:), pointer   :: header
     character(len=1),  dimension(:), pointer :: col_type
     real(WP),  dimension(:), pointer :: val
  end type mfile

  ! Array of mfiles
  type(mfile), dimension(:), pointer :: mfiles

contains
  
  
  ! Dump the values to the files at each iteration
  ! ----------------------------------------------
  subroutine monitor_dump_values_timestep
    use parallel
    use fileio
    use time_info
    implicit none
    
    integer :: icol,offset
    character(len=str_long) :: line
    
    ! Only the root process does something
    if (irank.ne.iroot) return
    
    do ifile=1,nfiles
       
       ! Test if we need to dump the values
       if (mfiles(ifile)%freq.eq.1) then
          
          ! Create the line to dump
          write (line(1+0*col_len:),f2) ntime
          write (line(1+1*col_len:),f3) time
          offset = 2*col_len
          
          do icol=1,mfiles(ifile)%ncols
             
             select case(mfiles(ifile)%col_type(icol))
             case ('i')
                write (line(offset+1+(icol-1)*col_len:),f2) int(mfiles(ifile)%val(icol))
             case ('r')
                write (line(offset+1+(icol-1)*col_len:),f3) mfiles(ifile)%val(icol) 
             case ('d')
                write (line(offset+1+(icol-1)*col_len_d:),f5) mfiles(ifile)%val(icol) 
             end select
             
          end do
          
          ! Dump the header
          write(mfiles(ifile)%iunit,'(a)') trim(line)
       end if
       
    end do
    
    return
  end subroutine monitor_dump_values_timestep
  
  
  ! Dump the values to the files at each iteration
  ! -------------------------------------------------
  subroutine monitor_dump_values_iteration
    use parallel
    use fileio
    use time_info
    implicit none
    
    integer :: icol,offset
    character(len=str_long) :: line
    
    ! Only the root process does something
    if (irank.ne.iroot) return
    
    do ifile=1,nfiles
       
       ! Test if we need to dump the values
       if (mfiles(ifile)%freq.eq.2) then
          
          ! Create the line to dump
          write (line(1+0*col_len:),f2) ntime
          write (line(1+1*col_len:),f3) time
          write (line(1+2*col_len:),f2) niter
          offset = 3*col_len
          
          do icol=1,mfiles(ifile)%ncols
             
             select case(mfiles(ifile)%col_type(icol))
             case ('i')
                write (line(offset+1+(icol-1)*col_len:),f2) int(mfiles(ifile)%val(icol))
             case ('r')
                write (line(offset+1+(icol-1)*col_len:),f3) mfiles(ifile)%val(icol) 
             end select
             
          end do
          
          ! Dump the header
          write(mfiles(ifile)%iunit,'(a)') trim(line)
       end if
       
    end do
    
    return
  end subroutine monitor_dump_values_iteration
  

  ! Print some values to the screen
  ! -------------------------------
  subroutine monitor_to_screen
    use parallel
    use time_info
    implicit none

    real(WP) :: max_u,max_v,max_w,max_divg
    integer  :: icol
    
    ! Only the root process does something
    if (irank.ne.iroot) return
    
    ! Find the file
    loop: do ifile=1,nfiles
       if (trim(mfiles(ifile)%filename).eq.'velocity') exit loop
    end do loop
    
    ! Get the values
    do icol=1,mfiles(ifile)%ncols
       select case (trim(mfiles(ifile)%header(icol)))
       case ('max_u') 
          max_u = mfiles(ifile)%val(icol)
       case ('max_v') 
          max_v = mfiles(ifile)%val(icol)
       case ('max_w') 
          max_w = mfiles(ifile)%val(icol)
       case ('max_divg') 
          max_divg = mfiles(ifile)%val(icol)
       end select
    end do

    ! Print to the screen
    write(*,'(i12,a2,1ES12.5,1F12.4,4ES12.3)') &
         ntime,'  ',time,CFL,max_u,max_v,max_w,max_divg
    
    return
  end subroutine monitor_to_screen
  
end module monitor

! ============================= !
! Initialize the monitor module !
! ============================= !
subroutine monitor_pre_init
  use monitor
  use fileio
  use config
  use parallel
  implicit none
  integer :: ierr

  ! Set number of files to zero
  nfiles = 0

  ! But preallocate the number of files to some value
  allocate(mfiles(nfiles_max))

  ! Only the root process does the following
  if (irank.eq.iroot) then
    
     ! Create the directory
     !call CREATE_FOLDER("monitor")
     call system("mkdir -p monitor")

     ! Open log file
     file_log = iopen()
     open(file_log,file="monitor/log", form="formatted",iostat=ierr,status="REPLACE")
     
     ! Start to print on the screen
     write(*,'(a12,a2,6a12)') 'Step','  ','Time  ','CFLmax','Umax  ','Vmax  ','Wmax  ','Divergence'

  end if

  ! Initialize the timer
  call timing_pre_init

  ! Initialize the probes
  call probe_init

  ! Init the sub modules
  call monitor_conservation_init
  select case(trim(simu_type))
  case ("mixing layer")
     call monitor_mixinglayer_init
  case("hit")
     call monitor_hit_init
  case ("boundary layer")
     call monitor_boundarylayer_init
  case ("channel")
     call monitor_pipechannel_init
  case ("pipe")
     call monitor_pipechannel_init
  end select
  
  return
end subroutine monitor_pre_init


! ============================= !
! Dump the headers of the files !
! ============================= !
subroutine monitor_post_init
  use monitor
  use parallel
  use fileio
  implicit none
  
  integer :: ierr,icol,offset
  character(len=col_len) :: col
  character(len=str_long) :: header
  character(len=str_medium) :: filename,buffer
  logical :: twoLines
  integer :: index1
  
  ! Get the info about the timers
  call timing_post_init

  ! Only the root process does something
  if (irank.ne.iroot) return
  
  do ifile=1,nfiles

     twoLines = .false.
     filename = 'monitor/'//trim(mfiles(ifile)%filename)

     ! Open the file
     mfiles(ifile)%iunit = iopen()
     open(mfiles(ifile)%iunit,file=filename,form="formatted",iostat=ierr,status="REPLACE")

     ! Create the header
     select case (mfiles(ifile)%freq)
     case (1)
        write(header(1+0*col_len:),f1) 'Step'
        write(header(1+1*col_len:),f1) 'Time'
        offset = 2*col_len
     case (2)
        write(header(1+0*col_len:),f1) 'Step'
        write(header(1+1*col_len:),f1) 'Time'
        write(header(1+2*col_len:),f1) 'Niter'
        offset = 3*col_len
     end select

     ! Extract the first line
     do icol=1,mfiles(ifile)%ncols
        read (mfiles(ifile)%header(icol),'(a)') buffer
        index1 = index(trim(buffer),' ')

        if (index1.ne.0 .and. index1.lt.col_len-1) then
           twoLines = .true.
           read(buffer(1:index1),f1) col
        else
           read(buffer,f1) col
        end if
        if (mfiles(ifile)%col_type(icol).eq.'d') then
           write(header(offset+1+(icol-1)*col_len_d:),f4) trim(col)
        else
           write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
        end if
     end do
     
     ! Dump the header
     write(mfiles(ifile)%iunit,'(a)') trim(header)
     
     ! Dump second line if necessary
     if (twoLines) then
        header = ''
        do icol=1,mfiles(ifile)%ncols
           read (mfiles(ifile)%header(icol),'(a)') buffer
           index1 = index(trim(buffer),' ')
           if (index1.ne.0.and. index1.lt.col_len-1) then
              read(buffer(index1:),f1) col
           else
              col = ''
           end if
           if (mfiles(ifile)%col_type(icol).eq.'d') then
              write(header(offset+1+(icol-1)*col_len_d:),f4) trim(col)
           else
              write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
           end if
        end do
        write(mfiles(ifile)%iunit,'(a)') trim(header)
     end if

  end do

  return
end subroutine monitor_post_init


! ============================================= !
! Create a new file to monitor at each timestep !
! ============================================= !
subroutine monitor_create_file_step(filename,ncols)
  use monitor
  implicit none
  
  character(len=*), intent(in) :: filename
  integer, intent(in) :: ncols

  ! Add a file
  nfiles = nfiles+1
  if (nfiles.gt.nfiles_max) call die ('monitor_create_file_step: too many files to monitor')
  ifile  = nfiles
  
  ! Preset the values
  mfiles(ifile)%filename = filename
  mfiles(ifile)%freq     = 1
  mfiles(ifile)%ncols    = ncols

  ! Allocate the arrays
  allocate(mfiles(ifile)%header(ncols))
  allocate(mfiles(ifile)%col_type(ncols))
  allocate(mfiles(ifile)%val(ncols))

  return
end subroutine monitor_create_file_step


! ============================================== !
! Create a new file to monitor at each iteration !
! ============================================== !
subroutine monitor_create_file_iter(filename,ncols)
  use monitor
  implicit none
  
  character(len=*), intent(in) :: filename
  integer, intent(in) :: ncols

  ! Add a file
  nfiles = nfiles+1
  if (nfiles.gt.nfiles_max) call die ('monitor_create_file: too many files to monitor')
  ifile  = nfiles
  
  ! Preset the values
  mfiles(ifile)%filename = filename
  mfiles(ifile)%freq     = 2
  mfiles(ifile)%ncols    = ncols

  ! Allocate the arrays
  allocate(mfiles(ifile)%header(ncols))
  allocate(mfiles(ifile)%col_type(ncols))
  allocate(mfiles(ifile)%val(ncols))

  return
end subroutine monitor_create_file_iter


! ======================================= !
! Search for the file and set the pointer !
! ======================================= !
subroutine monitor_select_file(filename)
  use monitor
  implicit none
  
  character(len=*), intent(in) :: filename

  ! Find the file
  loop: do ifile=1,nfiles
     if (trim(mfiles(ifile)%filename).eq.trim(filename)) exit loop
  end do loop
  if (ifile.gt.nfiles) then
     print*,'filename : ', filename
     call die('monitor_select_file: unknown file to monitor')
  end if
  
  return
end subroutine monitor_select_file


! ========================== !
! Set the header of the file !
! ========================== !
subroutine monitor_set_header(icol,header,col_type)
  use monitor
  implicit none
  
  integer, intent(in) :: icol
  character(len=*), intent(in) :: header
  character(len=1), intent(in) :: col_type

  ! Test if enough columns
  if (icol.gt.mfiles(ifile)%ncols) then
     print*,'filename : ', mfiles(ifile)%filename
     print*,'column index', icol
     call die('monitor_set_header: column index too large')
  end if
  
  ! Set the header
  mfiles(ifile)%header(icol)   = trim(header)
  mfiles(ifile)%col_type(icol) = col_type

  return
end subroutine monitor_set_header


! ================================== !
! Set the values to dump in the file !
! ================================== !
subroutine monitor_set_single_value(int,val)
  use monitor
  implicit none
  
  integer,  intent(in) :: int
  real(WP), intent(in) :: val
  
  if (int.gt.mfiles(ifile)%ncols) then
     call die('monitor_set: too many values')
  else
     mfiles(ifile)%val(int) = val
  end if
  
  return
end subroutine monitor_set_single_value


! ================================== !
! Set the values to dump in the file !
! ================================== !
subroutine monitor_set_array_values(val)
  use monitor
  implicit none
  
  real(WP), dimension(mfiles(ifile)%ncols) :: val
  
  mfiles(ifile)%val = val
  
  return
end subroutine monitor_set_array_values


! ============================= !
! Monitor some additional stuff !
! ============================= !
subroutine monitor_log(text)
  use monitor
  use parallel
  implicit none

  character(len=*), intent(in) :: text
  integer, dimension(8) :: date
  
  if (irank.eq.iroot) then
     call date_and_time(VALUES=date)
     write(file_log,100) date(2),date(3),date(1),date(5),date(6),date(7),text
  end if
  
100 format(i2.2,'/',i2.2,'/',i4.4,' @ ',i2.2,':',i2.2,':',i2.2,a59)
  
  return
end subroutine monitor_log


! ======================= !
! Finalize the monitoring !
! ======================= !
subroutine monitor_finalize
  use monitor
  use parallel
  use fileio
  implicit none

  if (irank.eq.iroot) then
     ! Close the files
     do ifile=1,nfiles
        close(iclose(mfiles(ifile)%iunit))
     end do
     ! Close the log file
     close(iclose(file_log))
  end if
  
  return
end subroutine monitor_finalize


! ====================== !
! Monitor a subiteration !
! ====================== !
subroutine monitor_iteration
  use monitor
  implicit none
  
  ! Dump the data
  call monitor_dump_values_iteration
  
  return
end subroutine monitor_iteration


! ================== !
! Monitor a timestep !
! ================== !
subroutine monitor_timestep
  use monitor
  use config
  implicit none
  
  ! Compute the quantities relavant for each config to monitor
  call monitor_conservation_compute
  select case(trim(simu_type))
  case ("mixing layer")
     call monitor_mixinglayer_compute
  case("hit")
     call monitor_hit_compute
  case ("boundary layer")
     call monitor_boundarylayer_compute
  case ("channel")
     call monitor_pipechannel_compute
  case ("pipe")
     call monitor_pipechannel_compute
  end select
  
  ! Call the different monitor subroutines
  call velocity_monitor
  call scalar_monitor
  call combustion_monitor
  call implicit_monitor
  call pollutants_monitor
  call timing_monitor
  call probe_monitor
  
  ! Dump the data
  call monitor_dump_values_timestep
  
  ! Print data to screen
  call monitor_to_screen

  return
end subroutine monitor_timestep


