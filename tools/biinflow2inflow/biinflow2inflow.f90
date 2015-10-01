program biinflow2inflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none

  character(len=str_medium) :: input_name
  character(len=str_medium) :: filename,filename1,filename2
  integer  :: ierr,var,n,j,nscalar
  integer  :: iunit,iunit1,iunit2
  integer  :: ntime,ntime1,ntime2
  integer  :: ny,ny1,ny2
  integer  :: nz,nz1,nz2
  integer  :: nvar,nvar1,nvar2
  real(WP) :: dt,dt1,dt2
  real(WP) :: time,time1,time2
  character(len=str_short), dimension(:), pointer :: names,names1,names2
  real(WP), dimension(:,:), pointer :: inflow,inflow1,inflow2
  real(WP), dimension(:), pointer :: y,y1,y2
  real(WP), dimension(:), pointer :: z,z1,z2
  character(len=str_short), dimension(:), pointer :: list1,list2
  real(WP), dimension(:), pointer :: sc1,sc2
  real(WP)  :: U1,U2
     
  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)

  ! Open upper inflow file
  call parser_read('Upper TBL filename',filename1)
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  if (ierr.ne.0) stop "biinflow2inflow: cannot open file: upper"
  call BINARY_FILE_READ(iunit1,ntime1,1,kind(ntime1),ierr)
  call BINARY_FILE_READ(iunit1,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit1,nz1,1,kind(nz1),ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  call BINARY_FILE_READ(iunit1,dt1,1,kind(dt1),ierr)
  call BINARY_FILE_READ(iunit1,time1,1,kind(time1),ierr)
  allocate(names1(nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names1(var),str_short,kind(names1),ierr)
  end do
  
  ! Open lower inflow file
  call parser_read('Lower TBL filename',filename2)
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"r",ierr)
  if (ierr.ne.0) stop "biinflow2inflow: cannot open file: lower"
  call BINARY_FILE_READ(iunit2,ntime2,1,kind(ntime2),ierr)
  call BINARY_FILE_READ(iunit2,ny2,1,kind(ny2),ierr)
  call BINARY_FILE_READ(iunit2,nz2,1,kind(nz2),ierr)
  call BINARY_FILE_READ(iunit2,nvar2,1,kind(nvar2),ierr)
  call BINARY_FILE_READ(iunit2,dt2,1,kind(dt2),ierr)
  call BINARY_FILE_READ(iunit2,time2,1,kind(time2),ierr)
  allocate(names2(nvar2))
  do var=1,nvar2
     call BINARY_FILE_READ(iunit2,names2(var),str_short,kind(names2),ierr)
  end do
  
  ! Check the consistency of the two inflows
  if (dt1.ne.dt2) &
       stop "biinflow2inflow: two inflows should have same dt"
  if (nz1.ne.nz2) &
       stop "biinflow2inflow: two inflows should have same nz"
  if (nvar1.ne.nvar2) &
       stop "biinflow2inflow: two inflows should have same variables"
  do var=1,min(nvar1,nvar2)
     if (trim(names1(var)).ne.trim(names2(var))) &
          stop "biinflow2inflow: two inflows should have same variables"
  end do
  ntime = min(ntime1,ntime2)
  ny = ny1+ny2-2
  nz = nz1
  dt = dt1
  time = 0.0_WP

  ! Account for additional scalars
  call parser_getsize('Upper scalar values',nscalar)
  nscalar = nscalar/2
  allocate(list1(2*nscalar))
  allocate(list2(2*nscalar))
  allocate(sc1(nscalar))
  allocate(sc2(nscalar))
  call parser_read('Upper scalar values',list1)
  call parser_read('Lower scalar values',list2)
  
  nvar = nvar1 + nscalar
  allocate(names(nvar))
  names(1:nvar1) = names1
  do var=1,nscalar
     read (list1(2*var-1),*) names(nvar1+var)
     read (list1(2*var),*) sc1(var)
     read (list2(2*var),*) sc2(var)
  end do
  
  ! Read the free stream values
  call parser_read('Upper U velocity',U1)
  call parser_read('Lower U velocity',U2)
  
  ! Print some stats
  print*
  print*,'Common parameters:'
  print*,'-> nz',nz
  print*,'-> dt',dt
  print*
  print*,'Updated parameters:'
  write(*,'(a9,3i6)') ' -> ntime',ntime1,ntime2,ntime
  write(*,'(a9,3i6)') ' -> ny   ',ny1,ny2,ny
  print*
  print*,'Variables '
  print*,'-> ',names
  print*

  ! Read the grids
  allocate(y1(ny1+1))
  allocate(z1(nz1+1))
  allocate(y2(ny2+1))
  allocate(z2(nz2+1))
  call BINARY_FILE_READ(iunit1,y1,ny1+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit1,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit2,y2,ny2+1,kind(y2),ierr)
  call BINARY_FILE_READ(iunit2,z2,nz2+1,kind(z2),ierr)

  ! Create the new grid
  allocate(y(-ny2+1:ny1-1))
  allocate(z(nz+1))
  z = z1
  do j=2,ny1+1
     y(j-2) = y1(j)
  end do
  do j=2,ny2+1
     y(2-j) = -y2(j)
  end do

  ! Write the new header
  call parser_read('Inflow filename',filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  call BINARY_FILE_WRITE(iunit,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
  end do
  call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)

  ! Read and write data field
  allocate(inflow (-ny2+2:ny1-1,nz))
  allocate(inflow1(ny1,nz1))
  allocate(inflow2(ny2,nz2))
  do n=1,ntime
     do var=1,nvar1
        ! Read 
        call BINARY_FILE_READ(iunit1,inflow1(:,:),ny1*nz1,kind(inflow1),ierr)
        call BINARY_FILE_READ(iunit2,inflow2(:,:),ny2*nz2,kind(inflow2),ierr)
        ! Convert
        select case(trim(names(var)))
        case ('U')
           do j=2,ny1
              inflow(j-1,:) = inflow1(j,:)
           end do
           do j=2,ny2
              inflow(2-j,:) = inflow2(j,:)
           end do
           inflow(2-ny2,:) = U2
           inflow(ny1-1,:) = U1
        case ('V')
           do j=2,ny1
              inflow(j-2,:) = inflow1(j,:)
           end do
           do j=2,ny2
              inflow(2-j,:) = -inflow2(j,:)
           end do
           inflow(2-ny2,:) = 0.0_WP
           inflow(ny1-2,:) = 0.0_WP
           inflow(ny1-1,:) = 0.0_WP
        case ('W')
           do j=2,ny1
              inflow(j-1,:) = inflow1(j,:)
           end do
           do j=2,ny2
              inflow(2-j,:) = -inflow2(j,:)
           end do
           inflow(2-ny2,:) = 0.0_WP
           inflow(ny1-1,:) = 0.0_WP
        case default
           do j=2,ny1
              inflow(j-1,:) = inflow1(j,:)
           end do
           do j=2,ny2
              inflow(2-j,:) = inflow2(j,:)
           end do
        end select
        ! Write
        call BINARY_FILE_WRITE(iunit,inflow(:,:),ny*nz,kind(inflow),ierr)
     end do
     do var=1,nscalar
        ! Convert
        do j=2,ny1
           inflow(j-1,:) = sc1(var)
        end do
        do j=2,ny2
           inflow(2-j,:) = sc2(var)
        end do
        ! Write
        call BINARY_FILE_WRITE(iunit,inflow(:,:),ny*nz,kind(inflow),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  

end program biinflow2inflow
