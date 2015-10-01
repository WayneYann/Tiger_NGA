program bl2ml
  use precision
  use string
  use fileio
  use cli_reader
  implicit none

  character(len=str_medium) :: filename
  integer :: iunit,ierr
  
  ! Boundary Layer
  integer  :: nx1,ny1,nz1,nvar1
  real(WP), dimension(:),       pointer :: xBL,yBL,zBL
  integer,  dimension(:,:),     pointer :: maskBL
  real(WP), dimension(:,:,:,:), pointer :: dataBL
  real(WP) :: Uinf
  
  ! Mixing Layer
  integer  :: nx,ny,nz,nvar
  integer  :: xper,yper,zper
  integer  :: icyl
  real(WP) :: dt,time
  real(WP), dimension(:),   pointer :: xML,yML,zML
  integer,  dimension(:,:), pointer :: maskML
  character(len=str_medium) :: config
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: dataML
  
  ! Variables
  integer :: i,j,var
  
  
  ! Read file names from standard input
  print*,'==========================================='
  print*,'| ARTS - Mixing Layer from Boundary Layer |'
  print*,'==========================================='

  
  ! *** BL CONFIG ***
  
  ! Open the config file of the Boundary Layer
  print*
  print "(a18,$)", " BL config file : "
  read "(a)", filename
  call BINARY_FILE_OPEN(iunit,trim(filename),"r",ierr)
  
  ! Read sizes
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  print*,'BL Grid :',nx1,'x',ny1,'x',nz1
  
  ! Read grid field
  allocate(xBL(nx1+1),yBL(ny1+1),zBL(nz1+1))
  allocate(maskBL(nx1,ny1))
  call BINARY_FILE_READ(iunit,xBL,nx1+1,kind(xBL),ierr)
  call BINARY_FILE_READ(iunit,yBL,ny1+1,kind(yBL),ierr)
  call BINARY_FILE_READ(iunit,zBL,nz1+1,kind(zBL),ierr)
  call BINARY_FILE_READ(iunit,maskBL,nx1*ny1,kind(maskBL),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! *** BL DATA ***
  
  ! Open the data file of the Boundary Layer
  print "(a16,$)", " BL data file : "
  read "(a)", filename
  call BINARY_FILE_OPEN(iunit,trim(filename),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  call BINARY_FILE_READ(iunit,nvar1,1,kind(nvar1),ierr)
  print*,'BL data :',nx1,'x',ny1,'x',nz1
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
  
  ! Read variable names
  allocate(names(nvar1+1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names

  ! Read data field
  allocate(dataBL(nx1,ny1,nz1,nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit,dataBL(:,:,:,var),nx1*ny1*nz1,kind(dataBL),ierr)
     print*,maxval(abs(dataBL(:,:,:,var))), ' at ', maxloc(abs(dataBL(:,:,:,var)))
  end do
  call BINARY_FILE_CLOSE(iunit,ierr)
  

  ! *** BL TO ML CONVERSION ***
  
  ! Config
  config = 'Mixing Layer'
  ! Grid
  nx = nx1
  ny = 2*(ny1-1)
  nz = nz1
  allocate(xML(nx+1),yML(-ny1+1:ny1-1),zML(nz+1))
  xML = xBL
  zML = zBL
  do j=2,ny1+1
     yML(j-2) = +yBL(j)
     yML(2-j) = -yBL(j)
  end do
  ! Mask
  allocate(maskML(nx,ny))
  maskML = 0
  
  ! Data
  dt   = 0.0_WP
  time = 0.0_WP
  nvar = nvar1+1
  allocate(dataML(nx,-ny1+2:ny1-1,nz,nvar))
  ! Convert
  do var=1,nvar1
     select case(trim(names(var)))
     case ('U')
        Uinf = sum(dataBL(:,ny1,:,var))/real(nx1*nz1,WP)
        dataBL(:,ny1,:,var) = Uinf
        Print*,"Uinf = ",Uinf
        do j=2,ny1
           do i=1,nx1
              dataML(i,j-1,:,var) = +dataBL(i,j,:,var)
              dataML(i,2-j,:,var) = -dataBL(nx1-i+1,j,:,var)
           end do
        end do
     case ('V')
        dataBL(:,ny1,:,var) = 0.0_WP
        do j=2,ny1
           do i=1,nx1
              dataML(i,j-2,:,var) = +dataBL(i,j,:,var)
              dataML(i,2-j,:,var) = -dataBL(nx1-i+1,j,:,var)
           end do
        end do
     case ('W')
        dataBL(:,ny1,:,var) = 0.0_WP
        do j=2,ny1
           do i=1,nx1
              dataML(i,j-1,:,var) = dataBL(i,j,:,var)
              dataML(i,2-j,:,var) = dataBL(nx1-i+1,j,:,var)
           end do
        end do
     case default
        do j=2,ny1
           do i=1,nx1
              dataML(i,j-1,:,var) = dataBL(i,j,:,var)
              dataML(i,2-j,:,var) = dataBL(nx1-i+1,j,:,var)
           end do
        end do
     end select
  end do
  ! Scalar
  names(nvar) = 'ZMIX'
  dataML(:,-ny1+2: 0,:,nvar) = 0.0_WP
  dataML(:,+1:+ny1-1,:,nvar) = 1.0_WP
  
  
  ! *** ML CONFIG ***
  
  ! Open the config file of the Boundary Layer
  print*
  print "(a18,$)", " ML config file : "
  read "(a)", filename
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_WRITE(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_WRITE(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
  print*,'ML Grid :',nx,'x',ny,'x',nz
  
  ! Write grid field
  call BINARY_FILE_WRITE(iunit,xML,nx+1,kind(xML),ierr)
  call BINARY_FILE_WRITE(iunit,yML,ny+1,kind(yML),ierr)
  call BINARY_FILE_WRITE(iunit,zML,nz+1,kind(zML),ierr)
  call BINARY_FILE_WRITE(iunit,maskML,nx*ny,kind(maskML),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! *** ML DATA ***
  
  ! Open the data file of the Boundary Layer
  print "(a16,$)", " ML data file : "
  read "(a)", filename
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
  
  ! Write additional stuff
  call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
  
  ! Read variable names
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names

  ! Read data field
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,dataML(:,:,:,var),nx*ny*nz,kind(dataML),ierr)
     print*,maxval(abs(dataML(:,:,:,var))), ' at ', maxloc(abs(dataML(:,:,:,var)))
  end do
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  
end program bl2ml
