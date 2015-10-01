program create1DChemtable
  use table_1D
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
  call parser_read("Flamelet file", flameletfile)
  
  ! Initialize the modules
  call flamelet_1D_init
  call table_1D_init

  ! -------------------------------------------------------

  print*,''
  print*,'** Convoluting flamelet **'
  write(*,'(a)') trim(flameletfile)

  ! Read the file
  call flamelet_1D_readfile(flameletfile)
  ! Convolute with PDF
  call table_1D_convolute
  ! Deallocate the data array
  call flamelet_1D_cleanup
  
  ! Change the variables names
  call table_1D_convert_names
  ! Compute the table
  call table_1D_setup
  ! Print some statistics
  call table_1D_stats
  ! Write the table
  call table_1D_write

end program create1DChemtable



! -----------------------------------------
subroutine commandline_args(input_name)
  implicit none
  integer, external :: iargc
  character(len=*) :: input_name
  external :: getarg
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
