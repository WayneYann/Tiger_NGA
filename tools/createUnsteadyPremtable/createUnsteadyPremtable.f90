program createUnsteadyPremtable
  use unsteady_prem_table
  use parser
  use precision
  use string
  implicit none

  integer :: n
  
  ! Filename of input file
  character(len=str_medium) :: name

  ! -------------------------------------------------------

  ! Parse the input file
  call commandline_args(name)
  call parser_init
  call parser_parsefile(name)
  
  ! Read the list of files
  call parser_getsize("List of Flamelets", nfiles)
  allocate(files(nfiles))
  call parser_read("List of Flamelets", files)
  
  ! Initialize the modules
  call unsteady_prem_flamelet_init
  call unsteady_prem_table_init

  ! -------------------------------------------------------

  print*,''
  print*,'** Files in the table **'
  do n=1,nfiles
     write(*,'(a)') trim(files(n))
     ! Read the file
     call unsteady_prem_flamelet_readfile(files(n))
     ! Convolute with PDF
     call unsteady_prem_table_convolute(n)
     ! Deallocate the data array
     call unsteady_prem_flamelet_cleanup
  end do
  
  ! Change the variables names
  call unsteady_prem_table_convert_names
  ! Compute the table
  call unsteady_prem_table_setup
  ! Extent the bounds of the chemtable
  !call unsteady_table_extent
  ! Print some statistics
  call unsteady_prem_table_stats
  ! Write the table
  call unsteady_prem_table_write

end program createUnsteadyPremtable



! -----------------------------------------
subroutine commandline_args(input_name)
  implicit none
  integer, external :: iargc
  character(len=*) :: input_name
  external :: getarg
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
