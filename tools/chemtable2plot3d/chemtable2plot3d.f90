program chemtable2plot3d
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: n1,n2,n3,nvar
  integer :: i,j,k
  character(len=str_medium), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data8
  real(SP), dimension(:,:,:), pointer :: data4
  integer :: iunit,iunit8,iunit4,ierr,var
  character(len=str_medium) :: filename1,filename2,filename3,buffer
  character(len=str_medium) :: combModel
  real(WP), dimension(:), pointer :: x1,x2,x3
  real(SP), dimension(:,:,:), pointer :: grid
  integer, dimension(:,:,:), pointer :: iblank

  ! Read file name from standard input
  print*,'========================================'
  print*,'| ARTS - chemtable to PLOT3D converter |'
  print*,'========================================'
  print*
  print "(a18,$)", " chemtable file : "
  read "(a)", filename1
  print "(a20,$)", " PLOT3D grid file : "
  read "(a)", filename2
  print "(a20,$)", " PLOT3D data file : "
  read "(a)", filename3
  buffer = filename3

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit8,trim(filename1),"r",ierr)
  ! ** Open the grid file to write **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"w",ierr)
  ! ** Open the data file to write **
  call BINARY_FILE_OPEN(iunit4,trim(filename3),"w",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit8,n1,1,kind(n1),ierr)
  call BINARY_FILE_READ(iunit8,n2,1,kind(n2),ierr)
  call BINARY_FILE_READ(iunit8,n3,1,kind(n3),ierr)
  call BINARY_FILE_READ(iunit8,nvar,1,kind(nvar),ierr)
  print*,'Grid :',n1,'x',n2,'x',n3
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit4,n1,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit4,n2,1,kind(n2),ierr)
  call BINARY_FILE_WRITE(iunit4,n3,1,kind(n1),ierr)
  call BINARY_FILE_WRITE(iunit4,nvar,1,kind(nvar),ierr)

  call BINARY_FILE_WRITE(iunit,n1,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit,n2,1,kind(n2),ierr)
  call BINARY_FILE_WRITE(iunit,n3,1,kind(n1),ierr)

  ! Read the axis coordinates
  allocate(x1(n1),x2(n2),x3(n3))
  call BINARY_FILE_READ(iunit8,x1,n1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit8,x2,n2,kind(x2),ierr)
  call BINARY_FILE_READ(iunit8,x3,n3,kind(x3),ierr)

  ! Write the axis coordinates
  allocate(grid(n1,n2,n3))
  do i=1,n1
     grid(i,:,:) = x1(i)/maxval(abs(x1))
  end do
  call BINARY_FILE_WRITE(iunit,grid,n1*n2*n3,kind(grid),ierr)
  do j=1,n2
     grid(:,j,:) = x2(j)/maxval(abs(x2))
  end do
  call BINARY_FILE_WRITE(iunit,grid,n1*n2*n3,kind(grid),ierr)
  do k=1,n3
     grid(:,:,k) = x3(k)/maxval(abs(x3))
  end do
  call BINARY_FILE_WRITE(iunit,grid,n1*n2*n3,kind(grid),ierr)
  
  ! Read the mask
  allocate(iblank(n1,n2,n3))
  call BINARY_FILE_READ(iunit8,iblank,n1*n2*n3,kind(iblank),ierr)
!!$  iblank = 0
  
  ! Write the Iblanks and close the file
  iblank = 1-iblank
  call BINARY_FILE_WRITE(iunit,iblank,n1*n2*n3,kind(iblank),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Read additional stuff
  call BINARY_FILE_READ(iunit8,combModel,str_medium,kind(combModel),ierr)
  print*,'Combustion Model : ',trim(combModel)

  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit8,names(var),str_medium,kind(names),ierr)
  end do
  
  ! Allocate arrays
  allocate(data8(n1,n2,n3))
  allocate(data4(n1,n2,n3))
  
  do var=1,nvar
     ! Read data field
     call BINARY_FILE_READ(iunit8,data8,n1*n2*n3,kind(data8),ierr)
     print"(a8,a9,e14.6,a8,e14.6)",trim(names(var)),' -> min: ',minval(data8),' - max: ',maxval(data8)

     ! Convert
     data4 = data8

     ! Write data field
     call BINARY_FILE_WRITE(iunit4,data4,n1*n2*n3,kind(data4),ierr)
  end do

  ! Close the files
  call BINARY_FILE_CLOSE(iunit8,ierr)
  call BINARY_FILE_CLOSE(iunit4,ierr)
  
  ! ** Open the name file to write **
  filename3 = trim(filename3) // '.nam'
  iunit = iopen()
  open (iunit, file=filename3, form="formatted", iostat=ierr)

  ! Write variable names
  do var=1,nvar
     write(iunit,'(a)') trim(names(var)) 
  end do
  close(iclose(iunit))

  ! ** Open the ensight result file to write **
  filename3 = trim(buffer) // '.res'
  iunit = iopen()
  open (iunit, file=filename3, form="formatted", iostat=ierr)
  
  ! Write variable names
  write(iunit,'(i3,x,a3)') nvar,'0 0'
  write(iunit,'(a1)') '1'
  write(iunit,'(a3)') '0.0'
  do var=1,nvar
     write(iunit,'(a,x,a1,x,i3,x,a)') trim(buffer),'F',var,trim(names(var))
  end do
  close(iclose(iunit))
  
end program chemtable2plot3d
