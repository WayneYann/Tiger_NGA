program mergeData
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar1,nvar2,nvar3
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data
  integer :: iunit1,iunit2,iunit3,ierr,var
  character(len=str_medium) :: filename1,filename2,filename3
  real(WP) :: dt,time

  ! Read file name from standard input
  print*,'======================'
  print*,'| ARTS - data Mergeor |'
  print*,'======================'
  print*
  print "(a15,$)", " data file 1 : "
  read "(a)", filename1
  print "(a15,$)", " data file 2 : "
  read "(a)", filename2
  print "(a20,$)", " data file merged : "
  read "(a)", filename3

  ! ** Open the data file 1 to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! ** Open the data file 2 to read **
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit2,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit2,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit2,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit2,nvar2,1,kind(nvar2),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit2,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit2,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  

  ! Read variable names and merge them
  nvar3 = nvar1+nvar2
  allocate(names(nvar3))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  do var=nvar1+1,nvar3
     call BINARY_FILE_READ(iunit2,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  
  ! Allocate arrays
  allocate(data(nx,ny,nz))
  
  ! Write header
  time = 0.0_WP
  call BINARY_FILE_OPEN(iunit3,trim(filename3),"w",ierr)
  call BINARY_FILE_WRITE(iunit3,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit3,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit3,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit3,nvar3,1,kind(nvar3),ierr)
  call BINARY_FILE_WRITE(iunit3,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit3,time,1,kind(time),ierr)
  do var=1,nvar3
     call BINARY_FILE_WRITE(iunit3,names(var),str_short,kind(names),ierr)
  end do
  
  ! Read data and write them
  do var=1,nvar1
     call BINARY_FILE_READ (iunit1,data,nx*ny*nz,kind(data),ierr)
     call BINARY_FILE_WRITE(iunit3,data,nx*ny*nz,kind(data),ierr)
  end do
  do var=nvar1+1,nvar3
     call BINARY_FILE_READ (iunit2,data,nx*ny*nz,kind(data),ierr)
     call BINARY_FILE_WRITE(iunit3,data,nx*ny*nz,kind(data),ierr)
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  call BINARY_FILE_CLOSE(iunit3,ierr)
  
end program mergeData
