program rescaler
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none
  
  character(len=str_medium) :: input_name
  character(len=str_medium) :: file_type
  
  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)
  
  ! Check what needs to be done
  call parser_read('File type',file_type)
  select case (trim(file_type))
  case ('inflow')
     call rescale_inflow
  case ('G inflow')
     call G_inflow
  case ('rotate inflow')
     call rotate_inflow
  case default
     stop 'Not implemented yet: DIY.'
  end select
  
end program rescaler


subroutine rescale_inflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none
  
  character(len=str_medium) :: filename1,filename2
  integer  :: n,var
  integer  :: iunit1,iunit2,ierr
  integer  :: ntime
  integer  :: ny
  integer  :: nz
  integer  :: nvar
  integer  :: icyl
  real(WP) :: dt1,dt2
  real(WP) :: time1,time2
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:), pointer :: inflow1,inflow2
  real(WP), dimension(:), pointer :: y1,y2
  real(WP), dimension(:), pointer :: z1,z2
  real(WP) :: u_scale,d_scale,t_scale
  
  ! Filenames
  call parser_read('Inflow to read',filename1)
  call parser_read('Inflow to write',filename2)
  
  ! Read inflow file
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  if (ierr.ne.0) stop "Inflow rescaler: cannot open inflow file to rescale."
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_READ(iunit1,dt1,1,kind(dt1),ierr)
  call BINARY_FILE_READ(iunit1,time1,1,kind(time1),ierr)
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  allocate(y1(ny+1)); call BINARY_FILE_READ(iunit1,y1,ny+1,kind(y1),ierr)
  allocate(z1(nz+1)); call BINARY_FILE_READ(iunit1,z1,nz+1,kind(z1),ierr)
  
  ! Scaling values
  call parser_read('Scaling velocity',u_scale)
  call parser_read('Scaling distance',d_scale)
  t_scale=d_scale/u_scale
  
  ! Rescale and write the new inflow
  time2=time1/t_scale
  dt2=dt1/t_scale
  allocate(y2(ny+1)); y2=y1/d_scale
  if (icyl.eq.0) then
     allocate(z2(nz+1)); z2=z1/d_scale
  else
     allocate(z2(nz+1)); z2=z1
  end if
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit2,dt2,1,kind(dt2),ierr)
  call BINARY_FILE_WRITE(iunit2,time2,1,kind(time2),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
  end do
  call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit2,y2,ny+1,kind(y2),ierr)
  call BINARY_FILE_WRITE(iunit2,z2,nz+1,kind(z2),ierr)
  allocate(inflow1(ny,nz))
  allocate(inflow2(ny,nz))
  do n=1,ntime
     do var=1,nvar
        ! Read
        call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
        ! Convert
        select case(trim(names(var)))
        case ('U')
           inflow2 = inflow1/u_scale
        case ('V')
           inflow2 = inflow1/u_scale
        case ('W')
           inflow2 = inflow1/u_scale
        case default
           stop "Unknown unit for a variable in inflow."
        end select
        ! Write
        call BINARY_FILE_WRITE(iunit2,inflow2,ny*nz,kind(inflow2),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
  return
end subroutine rescale_inflow


subroutine G_inflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none
  
  character(len=str_medium) :: filename1,filename2
  integer  :: n,var
  integer  :: iunit1,iunit2,ierr
  integer  :: ntime
  integer  :: ny
  integer  :: nz
  integer  :: nvar1,nvar2
  integer  :: icyl
  real(WP) :: dt
  real(WP) :: time
  character(len=str_short), dimension(:), pointer :: names1,names2
  real(WP), dimension(:,:), pointer :: inflow1,inflow2
  real(WP), dimension(:), pointer :: y1,y2
  real(WP), dimension(:), pointer :: z1,z2
  
  ! Filenames
  call parser_read('Inflow to read',filename1)
  call parser_read('Inflow to write',filename2)
  
  ! Read inflow file
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  if (ierr.ne.0) stop "Inflow rescaler: cannot open inflow file to rescale."
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  allocate(names1(nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names1(var),str_short,kind(names1),ierr)
  end do
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  allocate(y1(ny+1)); call BINARY_FILE_READ(iunit1,y1,ny+1,kind(y1),ierr)
  allocate(z1(nz+1)); call BINARY_FILE_READ(iunit1,z1,nz+1,kind(z1),ierr)
  
  ! Generate the new inflow
  allocate(y2(ny+1)); y2=y1
  allocate(z2(nz+1)); z2=z1
  nvar2=nvar1+1
  allocate(names2(nvar2))
  names2(1:nvar1)=names1
  names2(nvar2)='G'
  
  ! Write the new inflow
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit2,nvar2,1,kind(nvar2),ierr)
  call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
  do var=1,nvar2
     call BINARY_FILE_WRITE(iunit2,names2(var),str_short,kind(names2),ierr)
  end do
  call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit2,y2,ny+1,kind(y2),ierr)
  call BINARY_FILE_WRITE(iunit2,z2,nz+1,kind(z2),ierr)
  allocate(inflow1(ny,nz))
  allocate(inflow2(ny,nz))
  do n=1,ntime
     do var=1,nvar2
        select case(trim(names2(var)))
        case ('U')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        case ('V')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        case ('W')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        case ('G')
           ! Set
           inflow2(1:ny-1,:) = 1.0_WP
           inflow2(ny,:) = 0.0_WP
        case default
           stop "Unknown unit for a variable in inflow."
        end select
        ! Write
        call BINARY_FILE_WRITE(iunit2,inflow2,ny*nz,kind(inflow2),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
  return
end subroutine G_inflow


subroutine rotate_inflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none
  
  character(len=str_medium) :: filename1,filename2
  integer  :: n,var
  integer  :: iunit1,iunit2,ierr
  integer  :: ntime
  integer  :: ny
  integer  :: nz
  integer  :: nvar1,nvar2
  integer  :: icyl
  real(WP) :: dt
  real(WP) :: time
  character(len=str_short), dimension(:), pointer :: names1,names2
  real(WP), dimension(:,:), pointer :: inflow1,inflow2
  real(WP), dimension(:), pointer :: y1,y2
  real(WP), dimension(:), pointer :: z1,z2
  
  ! Filenames
  call parser_read('Inflow to read',filename1)
  call parser_read('Inflow to write',filename2)
  
  ! Read inflow file
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  if (ierr.ne.0) stop "Inflow rescaler: cannot open inflow file to rescale."
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  allocate(names1(nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names1(var),str_short,kind(names1),ierr)
  end do
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  allocate(y1(ny+1)); call BINARY_FILE_READ(iunit1,y1,ny+1,kind(y1),ierr)
  allocate(z1(nz+1)); call BINARY_FILE_READ(iunit1,z1,nz+1,kind(z1),ierr)
  
  ! Generate the new inflow
  allocate(y2(ny+1)); y2=y1
  allocate(z2(nz+1)); z2=z1
  nvar2=nvar1
  allocate(names2(nvar2))
  names2(1)='V' ! Should be U
  names2(2)='U' ! Should be V
  names2(3)='W' ! Should be W
  
  ! Write the new inflow
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit2,nvar2,1,kind(nvar2),ierr)
  call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
  do var=1,nvar2
     call BINARY_FILE_WRITE(iunit2,names2(var),str_short,kind(names2),ierr)
  end do
  call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit2,y2,ny+1,kind(y2),ierr)
  call BINARY_FILE_WRITE(iunit2,z2,nz+1,kind(z2),ierr)
  allocate(inflow1(ny,nz))
  allocate(inflow2(ny,nz))
  do n=1,ntime
     do var=1,nvar2
        select case(trim(names2(var)))
        case ('U')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        case ('V')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        case ('W')
           ! Read
           call BINARY_FILE_READ(iunit1,inflow1,ny*nz,kind(inflow1),ierr)
           ! Convert
           inflow2 = inflow1
        !case ('G')
        !   ! Set
        !   inflow2(1:ny-1,:) = 1.0_WP
        !   inflow2(ny,:) = 0.0_WP
        case default
           stop "Unknown unit for a variable in inflow."
        end select
        ! Write
        call BINARY_FILE_WRITE(iunit2,inflow2,ny*nz,kind(inflow2),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
  return
end subroutine rotate_inflow
