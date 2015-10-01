module optdata
  use data
  implicit none
  
  ! =========================================== !
  ! Additional storage space for any other data !
  ! =========================================== !
  
  ! Do we use this additional file or not
  logical :: optdata_present
  
  ! Everything in one big array
  integer :: nod
  real(WP), dimension(:,:,:,:), pointer :: OD
  character(len=str_short), dimension(:), pointer :: OD_name
  
  ! Write frequency
  integer :: optdata_last
  
  ! Optional data time info
  real(WP) :: optdata_dt,optdata_time
  
  ! Array with each variables
  integer :: MPI_IO_NVARS_OPT
  type(MPI_IO_VAR), dimension(:), pointer :: MPI_IO_DATA_OPT
  
contains
  
  ! =============================================================== !
  ! Define the variables to read/write and the views to the file    !
  !                                                                 !
  ! To read/write a new variable:                                   !
  !   - make a link to the new data array                           !
  ! =============================================================== !
  subroutine optdata_mpi_init
    implicit none
    
    integer, dimension(4) :: gsizes,lsizes,start
    integer :: ierr,var
    
    ! Put everything in one array
    nod = MPI_IO_NVARS_OPT
    
    ! Allocate the memory
    if (nod.ne.0) then
       allocate(OD(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nod))
       allocate(OD_name(1:nOD))
    end if
    
    ! Link the variables
    do var=1,MPI_IO_NVARS_OPT
       MPI_IO_DATA_OPT(var)%var => OD(imin_:imax_,jmin_:jmax_,kmin_:kmax_,var)
       OD_name(var) = trim(MPI_IO_DATA_OPT(var)%name)
    end do
    
    ! Define global(g) and local(l) sizes
    gsizes(1) = nx
    gsizes(2) = ny
    gsizes(3) = nz
    gsizes(4) = MPI_IO_NVARS_OPT
    
    lsizes(1) = nx_
    lsizes(2) = ny_
    lsizes(3) = nz_
    lsizes(4) = 1
    
    ! Define starting points
    start(1) = imin_-imin
    start(2) = jmin_-jmin
    start(3) = kmin_-kmin
    
    ! Define the view for each variable
    do var=1,MPI_IO_NVARS_OPT
       
       start(4) = var - 1
       
       call MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,start,&
            MPI_ORDER_FORTRAN,MPI_REAL_WP,MPI_IO_DATA_OPT(var)%view,ierr)
       call MPI_TYPE_COMMIT(MPI_IO_DATA_OPT(var)%view,ierr)
       
    end do
    
    return
  end subroutine optdata_mpi_init
  
  
  ! ======================================== !
  ! Read the full 3D PLOT3D file in parallel !
  ! ======================================== !
  subroutine optdata_read
    implicit none
    
    integer :: ifile,ierr,var,data_size
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(4) :: dims
    integer(kind=MPI_Offset_kind) :: disp
    character(len=str_medium) :: filename
    
    ! Get the name of the file to read in
    call parser_is_defined('Optional data to read',optdata_present)
    if (.not.optdata_present) return
    call parser_read('Optional data to read',filename)
    filename = trim(mpiiofs)//":" // trim(filename)
    
    ! Open the file
    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
    
    ! Read dimensions from header
    call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
    if ((dims(1).ne.nx) .or. (dims(2).ne.ny) .or. (dims(3).ne.nz)) then
       print*, 'grid = ',nx,ny,nz
       print*, 'data = ',dims(1),dims(2),dims(3)
       call die('The size of the optional data file does not correspond to the grid file')
    end if
    
    ! Read additional stuff
    call MPI_FILE_READ_ALL(ifile,optdata_dt,1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_READ_ALL(ifile,optdata_time,1,MPI_REAL_WP,status,ierr)
    
    ! Read the names and set the views to the file
    MPI_IO_NVARS_OPT = dims(4)
    allocate(MPI_IO_DATA_OPT(MPI_IO_NVARS_OPT))
    do var=1,MPI_IO_NVARS_OPT
       call MPI_FILE_READ_ALL(ifile,MPI_IO_DATA_OPT(var)%name,str_short,MPI_CHARACTER,status,ierr)
    end do
    call optdata_mpi_init()
    
    ! Displacement and size of arrays
    disp = 4*4 + str_short*MPI_IO_NVARS_OPT + 2*WP
    data_size = nx_*ny_*nz_
    
    ! Read the data for each variable
    do var=1,MPI_IO_NVARS_OPT
       call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_IO_DATA_OPT(var)%view, &
            "native",mpi_info,ierr)
       call MPI_FILE_READ_ALL(ifile,MPI_IO_DATA_OPT(var)%var,data_size, &
            MPI_REAL_WP,status,ierr)
    end do
    
    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)
    
    return
  end subroutine optdata_read
  
end module optdata


! ======================== !
! Data Initialization      !
!   - array allocation     !
!   - MPI structure for IO !
! ======================== !
subroutine optdata_init
  use optdata
  use simulation
  implicit none
  
  ! Frequency to write data file
  optdata_last = int(time/data_freq)
  
  ! Read in the data
  call optdata_read
  
  return
end subroutine optdata_init


! ===================================== !
! Test if we need to save the data file !
! ===================================== !
subroutine optdata_write(flag)
  use optdata
  use time_info
  implicit none
  logical, intent(in) :: flag
  
  ! Return if not used
  if (.not.optdata_present) return
  
  if (int(time/data_freq).ne.optdata_last .or. flag) then
     optdata_last = int(time/data_freq)
     call optdata_write_full3D
  end if
  
  return
end subroutine optdata_write


! ========================================= !
! Write the full 3D PLOT3D file in parallel !
! ========================================= !
subroutine optdata_write_full3D
  use optdata
  use time_info
  implicit none
  
  integer :: ifile,ierr,var,data_size
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  character(len=str_medium) :: filename,buffer
  integer :: overwrite
  logical :: file_is_there
  
  ! Get the name of the file to write to
  call parser_read('Optional data to write',filename)
  filename = trim(mpiiofs)//":" // trim(filename)
  
  ! Add time info to the file name
  call parser_read('Data overwrite',overwrite,1)
  if (overwrite.eq.0) then
     write(buffer,'(ES12.3)') time
     filename = trim(adjustl(filename))//'_'//trim(adjustl(buffer))
  end if
     
  ! Open the file to write
  inquire(file=filename,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(filename,mpi_info,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,ifile,ierr)
  
  ! Write header
  if (irank.eq.iroot) then
     ! Write dimensions
     dims(1) = nx
     dims(2) = ny
     dims(3) = nz
     dims(4) = MPI_IO_NVARS_OPT
     call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
     ! Write additional stuff
     call MPI_FILE_WRITE(ifile,dt,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     ! Write variable names
     do var=1,MPI_IO_NVARS_OPT
        call MPI_FILE_WRITE(ifile,MPI_IO_DATA_OPT(var)%name,str_short,MPI_CHARACTER,status,ierr)
     end do
  end if
  
  ! Displacement and size of arrays
  disp = 4*4 + str_short*MPI_IO_NVARS_OPT + 2*WP
  data_size = nx_*ny_*nz_
  
  ! Write the data for each variable
  do var=1,MPI_IO_NVARS_OPT
     call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_IO_DATA_OPT(var)%view, &
          "native",mpi_info,ierr)
     call MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA_OPT(var)%var,data_size, &
          MPI_REAL_WP,status,ierr)
  end do
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  ! Log
  call monitor_log("3D OPTIONAL DATA WRITTEN")
  
  return
end subroutine optdata_write_full3D


! ====================================== !
! Operations when simulation is finished !
!   - write out the data field           !
! ====================================== !
subroutine optdata_finalize
  use optdata
  implicit none
  
  ! Return if not used
  if (.not.optdata_present) return
  
  ! Write the data file
  call optdata_write_full3D
  
  return
end subroutine optdata_finalize
