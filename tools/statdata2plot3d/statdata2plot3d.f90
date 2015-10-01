program statdata2plot3d
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar
  character(len=str_medium), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data8
  real(SP), dimension(:,:), pointer :: data4
  integer :: iunit,iunit8,iunit4,ierr
  integer :: var
  character(len=str_medium) :: filename1,filename2,directory
  real(WP) :: dt,time

  ! Read file name from standard input
  print*,'===================================='
  print*,'| ARTS - stats to PLOT3D converter |'
  print*,'===================================='
  print*
  print "(a,$)", " stat file : "
  read "(a)", filename1
  print "(a20,$)", " PLOT3D directory : "
  read "(a)", directory

  !call CREATE_FOLDER(trim(directory))
  call system("mkdir -p "//trim(directory))
  filename2 = trim(directory) // '/data'
  
  ! ** Open the stat file to read **
  call BINARY_FILE_OPEN(iunit8,trim(filename1),"r",ierr)
  ! ** Open the data file to write **
  call BINARY_FILE_OPEN(iunit4,trim(filename2),"w",ierr)
  
  ! Read sizes
  call BINARY_FILE_READ(iunit8,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit8,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit8,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit8,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
  if (nz.ne.1) stop "stat should be 2D"
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit4,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit4,ny,1,kind(ny),ierr)
  !call BINARY_FILE_WRITE(iunit4,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit4,nvar,1,kind(nvar),ierr)
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit8,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit8,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit8,names(var),str_medium,kind(names),ierr)
  end do
  print*,'Variables : ',names
  
  ! Allocate arrays
  allocate(data8(nx,ny,nz))
  allocate(data4(nx,ny))
  
  do var=1,nvar
     ! Read data field
     call BINARY_FILE_READ(iunit8,data8,nx*ny*nz,kind(data8),ierr)
     print*,maxval(abs(data8)), ' at ', maxloc(abs(data8))
     
     ! Convert the data
     data4(:,:) = data8(:,:,1)
     
     ! Write data field
     call BINARY_FILE_WRITE(iunit4,data4,nx*ny,kind(data4),ierr)
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit8,ierr)
  call BINARY_FILE_CLOSE(iunit4,ierr)
  
  ! ** Open the name file to write **
  filename2 = trim(directory) // '/data.nam'
  iunit = iopen()
  open (iunit, file=filename2, form="formatted", iostat=ierr)

  ! Write variable names
  do var=1,nvar
     write(iunit,'(a)') trim(names(var)) 
  end do
  close(iclose(iunit))
  
  ! ** Open the ensight result file to write **
  filename2 = trim(directory) // '/data.res'
  iunit = iopen()
  open (iunit, file=filename2, form="formatted", iostat=ierr)

  ! Write variable names
  write(iunit,'(i3,x,a3)') nvar,'0 0'
  write(iunit,'(a1)') '1'
  write(iunit,'(a3)') '0.0'
  do var=1,nvar
     write(iunit,'(a,x,a1,x,i3,x,a)') 'data','F',var,trim(names(var))
  end do
  close(iclose(iunit))
  
end program statdata2plot3d
