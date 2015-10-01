program arts
  
  call main_init
  call simulation_run
  call main_stop
  
end program arts




subroutine main_init
  use string
  implicit none
  character(len=str_medium) :: input_name
  
  ! Initialize parallel environment
  call parallel_init
  
  ! Initialize the random number generator
  call random_init
  
  ! Parse the command line
  call parallel_get_inputname(input_name)
  
  ! Parse the input file
  call parser_init
  call parser_parsefile(input_name)
  
  ! Geometry initialization
  call geometry_init
 
  ! Data initialization
  call data_init
  call optdata_init
  
  ! Simulation initialization
  call simulation_init
  
  return
end subroutine main_init




subroutine main_stop
  implicit none
  
  ! Stop the simulation
  call simulation_finalize
  
  ! Finalize the parallel environment
  call parallel_final
  
  return
end subroutine main_stop





! -----------------------------------------
subroutine die(error_text)
  use parallel
  implicit none
  character(len=*), intent(in) :: error_text
  
  call monitor_log("KILLED")
  call parallel_kill(error_text)
  
  return
end subroutine die
