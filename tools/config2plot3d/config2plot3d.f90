program config2plot3d
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz
  integer :: xper,yper,zper
  integer :: icyl
  real(WP), dimension(:), pointer :: x,y,z
  real(SP), dimension(:,:,:), pointer :: grid
  integer, dimension(:,:), pointer :: mask
  integer, dimension(:,:,:), pointer :: iblank
  integer :: ierr,iunit
  character(len=str_medium) :: filename1,filename2
  character(len=str_medium) :: config
  integer :: i,j,k

  ! Read file name from standard input
  print*,'===================================='
  print*,'| ARTS - config to PLOT3D converter |'
  print*,'===================================='
  print*
  print "(a15,$)", " config file : "
  read "(a)", filename1
  print "(a15,$)", " PLOT3D file : "
  read "(a)", filename2

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
  allocate(x(nx+1),y(ny+1),z(nz+2))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Add the "periodic point"
  z(nz+2) = 2.0_WP*z(nz+1)-z(nz)
  nz = nz+1

  ! ** Open the grid file to write **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"w",ierr)

  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)

  ! Write the grid
  allocate(grid(nx,ny,nz))
  allocate(iblank(nx,ny,nz))
  
  ! Cartesian
  if (icyl.eq.0) then
     do i=1,nx
        grid(i,:,:) = 0.5_WP*(x(i)+x(i+1))
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     do j=1,ny
        grid(:,j,:) = 0.5_WP*(y(j)+y(j+1))
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     do k=1,nz
        grid(:,:,k) = 0.5_WP*(z(k)+z(k+1))
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)

  else ! Cylindrical
     do i=1,nx
        grid(i,:,:) = 0.5_WP*(x(i)+x(i+1))
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     do k=1,nz
        do j=1,ny
           grid(:,j,k) = 0.5_WP*(y(j)+y(j+1))*cos(0.5_WP*(z(k)+z(k+1)))
        end do
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     do k=1,nz
        do j=1,ny
           grid(:,j,k) = 0.5_WP*(y(j)+y(j+1))*sin(0.5_WP*(z(k)+z(k+1)))
        end do
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
  end if
  
  do k=1,nz
     iblank(:,:,k) = 1-mask(:,:)
  end do
  call BINARY_FILE_WRITE(iunit,iblank,nx*ny*nz,kind(iblank),ierr)

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

end program config2plot3d
