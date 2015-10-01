program dns2les
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar,size
  integer :: i,j,k,var
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data1
  real(WP), dimension(:,:,:), pointer :: data2
  integer :: iunit1,iunit2,ierr
  character(len=str_medium) :: filename1,filename2
  real(WP) :: dt,time

  ! Read file name from standard input
  print*,'================================'
  print*,'| ARTS - DNS to LES conversion !'
  print*,'================================'
  print*
  print "(a17,$)", " DNS data file : "
  read "(a)", filename1
  print "(a17,$)", " LES data file : "
  read "(a)", filename2
  print "(a30,$)", " Number of cells to average : "
  read "(i3)", size

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  ! ** Open the data file to write **
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz

  ! Some checks
  if (size.eq.0) stop "Number of cells to average cannot be zero."
  if (int(real(nx)/real(size)).ne.nx/size) stop "nx cannot be divided by size."
  if (int(real(ny)/real(size)).ne.ny/size) stop "ny cannot be divided by size."
  if (int(real(nz)/real(size)).ne.nz/size) stop "nz cannot be divided by size."
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit2,nx/size,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit2,ny/size,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit2,nz/size,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)

  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Write additional stuff
  call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)

  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
     call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names

  ! Allocate arrays
  allocate(data1(nx,ny,nz))
  allocate(data2(nx/size,ny/size,nz/size))

  do var=1,nvar
     ! Read data field
     call BINARY_FILE_READ(iunit1,data1,nx*ny*nz,kind(data1),ierr)
     print*,maxval(abs(data1(:,:,:))), ' at ', maxloc(abs(data1(:,:,:)))

     ! Average the data
     do k=1,nz/size
        do j=1,ny/size
           do i=1,nx/size
              data2(i,j,k) = sum(data1((i-1)*size+1:i*size,(j-1)*size+1:j*size,(k-1)*size+1:k*size)) / size**3
           end do
        end do
     end do
     
     ! Write data field
     call BINARY_FILE_WRITE(iunit2,data2,nx*ny*nz/size**3,kind(data2),ierr)
  end do

  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
end program dns2les
