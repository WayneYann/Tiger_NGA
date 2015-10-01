module utilities
  implicit none

end module utilities

! ======================================= !
! Directory creation routine              !
!   - checks for use_create_command flag  !
!   and uses right command                !
! ======================================= !
subroutine make_directory(name)
  use utilities
  use string
  use parser
  implicit none
  
  character(len=*),intent(in) :: name
  character(len=str_medium) :: temp
  ! Are we on a cluster needing Create Folder command
  logical :: use_create
  
  ! Read if cluster type is specified
  ! Set to false by default
  call parser_read('Use Create Folder',use_create,.false.)
  if (use_create) then
     call CREATE_FOLDER(trim(name))
  else
     temp = "mkdir -p "//trim(name)
     call system(trim(temp))
  end if
  return
end subroutine make_directory
