program createChemtable
  use table
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
  call flamelet_init
  call table_init

  ! -------------------------------------------------------

  print*,''
  print*,'** Files in the table **'
  do n=1,nfiles
     write(*,'(a)') trim(files(n))
     ! Read the file
     call flamelet_readfile(files(n))
     ! Convolute with PDF
     call table_convolute(n)
     ! Deallocate the data array
     call flamelet_cleanup
  end do
  
  ! Change the variables names
  call table_convert_names
  ! Compute the table
  call table_setup
  ! Extent the bounds of the chemtable
  !call table_extent
  ! Print some statistics
  call table_stats
  ! Write the table
  call table_write

end program createChemtable



! -----------------------------------------
subroutine commandline_args(input_name)
  implicit none
  integer, external :: iargc
  character(len=*) :: input_name
  external :: getarg
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
