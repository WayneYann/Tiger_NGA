program combineInflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none

  character(len=str_medium) :: input_name
  character(len=str_medium) :: filename,filename_pipe,filename_coflow
  integer  :: ierr,var,n,j,nscalar
  integer  :: iunit,iunit_pipe,iunit_coflow
  integer  :: ntime,ntime_pipe,ntime_coflow
  integer  :: ny,ny_pipe,ny_coflow
  integer  :: nz,nz_pipe,nz_coflow
  integer  :: nvar,nvar_pipe,nvar_coflow
  integer  :: icyl
  real(WP) :: dt,dt_pipe,dt_coflow
  real(WP) :: time,time_pipe,time_coflow
  character(len=str_short), dimension(:), pointer :: names,names_pipe,names_coflow
  real(WP), dimension(:,:), pointer :: inflow,inflow_pipe,inflow_coflow
  real(WP), dimension(:), pointer :: y,y_pipe,y_coflow
  real(WP), dimension(:), pointer :: z,z_pipe,z_coflow

  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)

  ! Open the pipe inflow
  call parser_read('Pipe inflow file',filename_pipe)
  call BINARY_FILE_OPEN(iunit_pipe,trim(filename_pipe),"r",ierr)
  if (ierr.ne.0) stop "combineInflow: cannot open pipe file"
  call BINARY_FILE_READ(iunit_pipe,ntime_pipe,1,kind(ntime_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,ny_pipe,1,kind(ny_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,nz_pipe,1,kind(nz_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,nvar_pipe,1,kind(nvar_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,dt_pipe,1,kind(dt_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,time_pipe,1,kind(time_pipe),ierr)
  allocate(names_pipe(nvar_pipe))
  do var=1,nvar_pipe
     call BINARY_FILE_READ(iunit_pipe,names_pipe(var),str_short,kind(names_pipe),ierr)
  end do

  ! Open the coflow inflow
  call parser_read('Coflow inflow file',filename_coflow)
  call BINARY_FILE_OPEN(iunit_coflow,trim(filename_coflow),"r",ierr)
  if (ierr.ne.0) stop "combineInflow: cannot open coflow file"
  call BINARY_FILE_READ(iunit_coflow,ntime_coflow,1,kind(ntime_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,ny_coflow,1,kind(ny_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,nz_coflow,1,kind(nz_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,nvar_coflow,1,kind(nvar_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,dt_coflow,1,kind(dt_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,time_coflow,1,kind(time_coflow),ierr)
  allocate(names_coflow(nvar_coflow))
  do var=1,nvar_coflow
     call BINARY_FILE_READ(iunit_coflow,names_coflow(var),str_short,kind(names_coflow),ierr)
  end do

  ! Check the consistency of the two inflows
  if (nz_pipe.ne.nz_coflow) stop "combineInflow: two inflows should have same nz"
  if (nvar_pipe.ne.nvar_coflow) stop "combineInflow: two inflows should have same variables"
  if (dt_pipe.ne.dt_coflow) stop "combineInflow: two inflows should have same dt"
  do var=1,nvar_pipe
     if (trim(names_pipe(var)).ne.trim(names_coflow(var))) &
          stop "combineInflow: two inflows should have same variables"
  end do
  ntime = min(ntime_pipe,ntime_coflow)
  ny = ny_pipe + ny_coflow + 1
  nz = nz_pipe
  nvar = nvar_pipe
  dt = dt_pipe
  time = 0.0_WP
  allocate(names(nvar))
  names = names_pipe

  ! Print some stats
  print*
  print*,'Common parameters:'
  print*,'-> nz',nz
  print*,'-> dt',dt
  print*
  print*,'Updated parameters:'
  write(*,'(a9,3i6)') ' -> ntime',ntime_pipe,ntime_coflow,ntime
  write(*,'(a9,3i6)') ' -> ny   ',ny_pipe,ny_coflow,ny
  print*
  print*,'Variables '
  print*,'-> ',names
  print*

  ! Read the grids
  allocate(y_pipe(ny_pipe+1))
  allocate(z_pipe(nz_pipe+1))
  call BINARY_FILE_READ(iunit_pipe,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit_pipe,y_pipe,ny_pipe+1,kind(y_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,z_pipe,nz_pipe+1,kind(z_pipe),ierr)
  allocate(y_coflow(ny_coflow+1))
  allocate(z_coflow(nz_coflow+1))
  call BINARY_FILE_READ(iunit_coflow,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit_coflow,y_coflow,ny_coflow+1,kind(y_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,z_coflow,nz_coflow+1,kind(z_coflow),ierr)

  ! Create the new grid
  allocate(y(ny+1))
  y(1:ny_pipe+1) = y_pipe
  y(ny_pipe+2:ny_pipe+ny_coflow+2) = y_coflow
  allocate(z(nz+1))
  z = z_pipe

  ! Write the new header
  call parser_read('Combined inflow file',filename)
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
  call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)

  ! Read and write data field
  allocate(inflow(ny,nz))
  allocate(inflow_pipe(ny_pipe,nz_pipe))
  allocate(inflow_coflow(ny_coflow,nz_coflow))
  do n=1,ntime
     do var=1,nvar
        inflow = 0.0_WP
        ! Read
        call BINARY_FILE_READ(iunit_pipe  ,inflow_pipe  (:,:),ny_pipe  *nz_pipe  ,kind(inflow_pipe  ),ierr)
        call BINARY_FILE_READ(iunit_coflow,inflow_coflow(:,:),ny_coflow*nz_coflow,kind(inflow_coflow),ierr)
        ! Convert
        inflow(1:ny_pipe,:) = inflow_pipe
        inflow(ny_pipe+2:ny_pipe+ny_coflow+1,:) = inflow_coflow
        ! Write
        call BINARY_FILE_WRITE(iunit,inflow(:,:),ny*nz,kind(inflow),ierr)
     end do
  end do

  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)
  call BINARY_FILE_CLOSE(iunit_pipe,ierr)
  call BINARY_FILE_CLOSE(iunit_coflow,ierr)

end program combineInflow
