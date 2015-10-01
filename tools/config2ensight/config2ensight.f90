program config2ensight
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz
  integer :: xper,yper,zper
  integer :: icyl
  real(WP), dimension(:), pointer :: xWP,yWP,zWP
  integer, dimension(:,:), pointer :: mask
  integer :: ierr,iunit
  character(len=str_medium) :: filename1,filename2, directory
  character(len=str_medium) :: config
  integer :: i,j,k

  integer, dimension(:,:,:), pointer :: iblank
  character(len=80) :: buffer
  real(SP), dimension(:), pointer :: x,y,z
  real(SP), dimension(:), pointer :: xm,ym,zm
  real(SP) :: max_x,max_y,max_z
  real(SP) :: min_x,min_y,min_z
  integer :: ipart

  ! Read file name from standard input
  print*,'==========================================='
  print*,'| ARTS - config to ENSIGHT GOLD converter |'
  print*,'==========================================='
  print*
  print "(a15,$)", " config file : "
  read "(a)", filename1
  print "(a21,$)", " ensight directory : "
  read "(a)", directory

  !call CREATE_FOLDER(trim(directory))
  call system("mkdir -p "//trim(directory))
  filename2 = trim(directory) // '/geometry'

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
  allocate(xWP(nx+1),yWP(ny+1),zWP(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,xWP,nx+1,kind(xWP),ierr)
  call BINARY_FILE_READ(iunit,yWP,ny+1,kind(yWP),ierr)
  call BINARY_FILE_READ(iunit,zWP,nz+1,kind(zWP),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Create arrays necessary for Ensight format
  allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
  allocate(xm(1:nx),ym(1:ny),zm(1:nz))
  allocate(iblank(nx,ny,nz))
  do k=1,nz
     iblank(:,:,k) = 1-mask
  end do

  x = xWP
  y = yWP
  z = zWP
  do i=1,nx
     xm(i) = 0.5_SP*(x(i)+x(i+1))
     !print*,'x',i,x(i)
  end do
  do j=1,ny
     ym(j) = 0.5_SP*(y(j)+y(j+1))
     !print*,'y',j,y(j)
  end do
  do k=1,nz
     zm(k) = 0.5_SP*(z(k)+z(k+1))
  end do

  max_x = x(nx+1)
  max_y = y(ny+1)
  max_z = z(nz+1)
  min_x = x(1)
  min_y = y(1)
  min_z = z(1)

  ! ** Open the grid file to write **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"w",ierr)

  ! Write the geometry
  buffer = 'C Binary'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'Ensight Gold Geometry File'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'Structured Geometry'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'node id off'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'element id off'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)

!!$  buffer = 'extents'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  call BINARY_FILE_WRITE(min_x,1,kind(min_x),ierr)
!!$  call BINARY_FILE_WRITE(max_x,1,kind(max_x),ierr)
!!$  call BINARY_FILE_WRITE(min_y,1,kind(min_y),ierr)
!!$  call BINARY_FILE_WRITE(max_y,1,kind(max_y),ierr)
!!$  call BINARY_FILE_WRITE(min_z,1,kind(min_z),ierr)
!!$  call BINARY_FILE_WRITE(max_z,1,kind(max_z),ierr)

  ! Cell centers
  buffer = 'part'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  ipart = 1
  call BINARY_FILE_WRITE(iunit,ipart,1,kind(ipart),ierr)

  buffer = 'Complete geometry'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'block rectilinear iblanked'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)

  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)

  call BINARY_FILE_WRITE(iunit,xm,nx,kind(xm),ierr)
  call BINARY_FILE_WRITE(iunit,ym,ny,kind(ym),ierr)
  call BINARY_FILE_WRITE(iunit,zm,nz,kind(zm),ierr)
  call BINARY_FILE_WRITE(iunit,iblank,nx*ny*nz,kind(iblank),ierr)

!!$  ! Face in x-direction
!!$  buffer = 'part'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  ipart = 2
!!$  call BINARY_FILE_WRITE(ipart,1,kind(ipart),ierr)
!!$
!!$  buffer = 'Complete geometry'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  buffer = 'block rectilinear'! iblanked'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$
!!$  call BINARY_FILE_WRITE(nx,1,kind(nx),ierr)
!!$  call BINARY_FILE_WRITE(ny,1,kind(ny),ierr)
!!$  call BINARY_FILE_WRITE(nz,1,kind(nz),ierr)
!!$
!!$  call BINARY_FILE_WRITE(x(2:nx+1),nx,kind(x),ierr)
!!$  call BINARY_FILE_WRITE(ym,ny,kind(ym),ierr)
!!$  call BINARY_FILE_WRITE(zm,nz,kind(zm),ierr)
!!$  !call BINARY_FILE_WRITE(iblank,nx*ny*nz,kind(iblank),ierr)
!!$
!!$  ! Face in y-direction
!!$  buffer = 'part'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  ipart = 3
!!$  call BINARY_FILE_WRITE(ipart,1,kind(ipart),ierr)
!!$
!!$  buffer = 'Complete geometry'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  buffer = 'block rectilinear'! iblanked'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$
!!$  call BINARY_FILE_WRITE(nx,1,kind(nx),ierr)
!!$  call BINARY_FILE_WRITE(ny,1,kind(ny),ierr)
!!$  call BINARY_FILE_WRITE(nz,1,kind(nz),ierr)
!!$
!!$  call BINARY_FILE_WRITE(xm,nx,kind(xm),ierr)
!!$  call BINARY_FILE_WRITE(y(2:ny+1),ny,kind(y),ierr)
!!$  call BINARY_FILE_WRITE(zm,nz,kind(zm),ierr)
!!$  !call BINARY_FILE_WRITE(iblank,nx*ny*nz,kind(iblank),ierr)
!!$
!!$  ! Face in z-direction
!!$  buffer = 'part'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  ipart = 4
!!$  call BINARY_FILE_WRITE(ipart,1,kind(ipart),ierr)
!!$
!!$  buffer = 'Complete geometry'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$  buffer = 'block rectilinear'! iblanked'
!!$  call BINARY_FILE_WRITE(buffer,80,kind(buffer),ierr)
!!$
!!$  call BINARY_FILE_WRITE(nx,1,kind(nx),ierr)
!!$  call BINARY_FILE_WRITE(ny,1,kind(ny),ierr)
!!$  call BINARY_FILE_WRITE(nz,1,kind(nz),ierr)
!!$
!!$  call BINARY_FILE_WRITE(xm,nx,kind(xm),ierr)
!!$  call BINARY_FILE_WRITE(ym,ny,kind(ym),ierr)
!!$  call BINARY_FILE_WRITE(z(2:nz+1),nz,kind(z),ierr)
!!$  !call BINARY_FILE_WRITE(iblank,nx*ny*nz,kind(iblank),ierr)

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

end program config2ensight
