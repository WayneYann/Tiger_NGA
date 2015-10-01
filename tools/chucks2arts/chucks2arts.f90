program chucks2arts
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer(4) :: nx4,ny4,nz4,nvar4
  integer :: nx,ny,nz,nvar
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data8
  integer :: iunit,ierr
  integer :: var
  character(len=str_medium) :: filename1,filename2
  real(WP) :: dt,time

  ! Read file name from standard input
  print*,'==================================='
  print*,'| ARTS - Chucks to ARTS converter |'
  print*,'==================================='
  print*
  print "(a32,$)", " Unformatted Chucks data file : "
  read "(a)", filename1
  print "(a25,$)", " Binary ARTS data file : "
  read "(a)", filename2

  ! ** Open the data file to read **
  iunit = iopen()
  open (iunit, file=filename1, form="unformatted", &
       status="old", iostat=ierr)
  if (ierr .ne. 0) stop "File not found"

  ! Read header
  read (iunit) nx4, ny4, nz4, nvar4
  print*,'Grid :',nx4,'x',ny4,'x',nz4

  ! Allocate
  allocate(data8(nx4,ny4,nz4,nvar4))

  ! Read data
  do var=1,nvar4
     read (iunit) data8(:,:,:,var)
     print*,maxval(abs(data8(:,:,:,var))), ' at ', maxloc(abs(data8(:,:,:,var)))
  end do

  ! Close the file
  close (iclose(iunit))

  ! ** Open the data file to write **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"w",ierr)

  ! Write sizes
  nx = nx4-2
  ny = ny4-2
  nz = nz4-2
  nvar = 4!nvar4
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
  
  ! Write additional stuff
  dt = 0.0_WP
  time = 0.0_WP
  call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
  
  ! Write variable names
  allocate(names(nvar))
  names(1) = 'U'
  names(2) = 'V'
  names(3) = 'W'
  names(4) = 'P'
  !names(5) = 'VISC'
  do var=6,nvar
     names(var) = 'Z' // achar(var)
  end do
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
  end do

  ! Write data field
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,data8(2:nx+1,2:nx+1,2:nx+1,var),nx*ny*nz,kind(data8),ierr)
  end do

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  

end program chucks2arts
