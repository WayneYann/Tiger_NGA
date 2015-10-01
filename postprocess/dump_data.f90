module dump_data
  use precision
  use string
  implicit none
  
  ! Number of outputs
  integer :: nout
  
  ! Type of outputs
  character(str_medium), dimension(:), allocatable :: output_type
  
  ! Frequency of output
  real(WP), dimension(:), allocatable :: output_freq
  
  ! Output info
  integer, dimension(:), allocatable :: output_save
  
end module dump_data

! ================================== !
! Initialize data dumping procedures !
! ================================== !
subroutine dump_init
  use dump_data
  use parser
  use time_info
  implicit none

  logical :: output_isdef
  integer :: i,nfreq
  real(WP) :: buffer
  
  ! Read the input file
  call parser_is_defined('Output type',output_isdef)
  nout = 0
  if (output_isdef) then
     
     call parser_getsize('Output type',nout)
     allocate(output_type(nout))
     call parser_read('Output type',output_type)
     
     call parser_getsize('Output frequency',nfreq)
     allocate(output_freq(nfreq))
     call parser_read('Output frequency',output_freq)
     if (nfreq.EQ.1) then
        buffer = output_freq(1)
        deallocate(output_freq)
        allocate(output_freq(nout))
        output_freq = buffer
     else if (nfreq.NE.nout) then
        call die('dump_data: not enough output frequencies specified')
     end if
     
     allocate(output_save(nout))
     output_save(:) = int(time/output_freq(:))
     
  end if
  
  ! Initialize
  do i=1,nout
     select case (trim(output_type(i)))
     case ('ensight-3D')
        call dump_ensight_3D_init
     case ('ensight-str-3D')
        call dump_ensight_str_3D_init
     case ('plot3D-2D')
        call dump_plot3D_2D_init
     case ('plot3D-3D')
        call dump_plot3D_3D_init
     case default
     end select
  end do
  
  return
end subroutine dump_init

! ======================================== !
! Dump_data at each time step if necessary !
! ======================================== !
subroutine dump_result
  use dump_data
  use time_info
  implicit none

  integer :: i
  
  ! Dump data if needed
  do i=1,nout
     select case (trim(output_type(i)))
     case ('ensight-3D')
        if (int(time/output_freq(i)).NE.output_save(i)) then
           output_save(i) = int(time/output_freq(i))
           call dump_ensight_3D
           call monitor_log("3D ENSIGHT FILE WRITTEN")
        end if
        
     case ('ensight-str-3D')
        if (int(time/output_freq(i)).NE.output_save(i)) then
           output_save(i) = int(time/output_freq(i))
           call dump_ensight_str_3D
           call monitor_log("3D ENSIGHT-STR FILE WRITTEN")
        end if
        
     case ('plot3D-2D')
        if (int(time/output_freq(i)).NE.output_save(i)) then
           output_save(i) = int(time/output_freq(i))
           call dump_plot3D_2D
           call monitor_log("2D PLOT3D FILE WRITTEN")
        end if
        
     case ('plot3D-3D')
        if (int(time/output_freq(i)).NE.output_save(i)) then
           output_save(i) = int(time/output_freq(i))
           call dump_plot3D_3D
           call monitor_log("3D PLOT3D FILE WRITTEN")
        end if
        
     case default
        call die('dump_data: unknown output type ('//trim(output_type(i))//')')
     end select
  end do
  
  return
end subroutine dump_result

