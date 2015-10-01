program statconfig2plot3d
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz
  integer :: xper,yper,zper
  integer :: icyl
  real(WP), dimension(:), pointer :: x,y,z
  real(SP), dimension(:,:), pointer :: grid
  integer, dimension(:,:), pointer :: mask
  integer, dimension(:,:), pointer :: iblank
  integer :: ierr,iunit
  character(len=str_medium) :: filename1,filename2,directory
  character(len=str_medium) :: config
  integer :: i,j

  ! Read file name from standard input
  print*,'===================================='
  print*,'| ARTS - config to PLOT3D converter |'
  print*,'===================================='
  print*
  print "(a,$)", " config file : "
  read "(a)", filename1
  print "(a20,$)", " PLOT3D directory : "
  read "(a)", directory

  !call CREATE_FOLDER(trim(directory))
  call system("mkdir -p "//trim(directory))
  filename2 = trim(directory) // '/grid'

  ! ** Open the config file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  print*,'Config : ',config
  print*,'Grid :',nx,'x',ny,'x',nz
  
  ! Read grid field
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! ** Open the grid file to write **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"w",ierr)
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  
  ! Write the grid
  allocate(grid(nx,ny))
  allocate(iblank(nx,ny))
  do i=1,nx
     grid(i,:) = 0.5_WP*(x(i)+x(i+1))
  end do
  call BINARY_FILE_WRITE(iunit,grid,nx*ny,kind(grid),ierr)
  do j=1,ny
     grid(:,j) = 0.5_WP*(y(j)+y(j+1))
  end do
  call BINARY_FILE_WRITE(iunit,grid,nx*ny,kind(grid),ierr)
  iblank = 1-mask
  call BINARY_FILE_WRITE(iunit,iblank,nx*ny,kind(iblank),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
end program statconfig2plot3d
