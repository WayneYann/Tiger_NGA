module inflow_generation
  use precision
  use geometry
  use partition
  use string
  implicit none
  
  ! Flag to indicate if we generate inflow profiles
  integer :: iinflow
  
  ! Name of output and file id
  character(len=str_medium) :: filename
  integer :: ifile
  
  ! Frequency of output
  real(WP) :: inflow_freq
  
  ! Location of inflow generation
  real(WP) :: inflow_loc
  integer  :: iloc_inf
  
  ! Number of time steps and last save
  integer  :: inflow_ntime
  real(WP) :: inflow_save
  
  ! Inflow values at the last plane
  integer :: nvar
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: buf
  integer, dimension(:), pointer :: inf_ind

  ! View in the file
  integer :: default_view
  integer :: data_size
  integer, dimension(:), pointer :: fileview
  
end module inflow_generation


! =================================== !
! Initialize Inflow Generation module !
! =================================== !
subroutine inflow_generation_init
  use inflow_generation
  use parallel
  use parser
  use time_info
  use data
  implicit none
  
  logical :: file_is_there,isdefined
  integer, dimension(3) :: gsizes, lsizes, start
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer :: ierr,var
  character(len=str_short), dimension(:), pointer :: names_tmp
  integer :: isc,n
  
  ! Do we output an inflow profile?
  call parser_is_defined('Inflow file to write',isdefined)
  
  ! If not return
  iinflow = 0
  if (isdefined) iinflow = 1
  if (iinflow.eq.0) return
  
  ! Read filename and frequency
  call parser_read('Inflow file to write',filename)
  call parser_read('Inflow frequency',inflow_freq)
  call parser_read('Inflow location',inflow_loc,0.5_WP*xL)
  
  ! Handle inflow generation location
  if (xper.eq.1) then
     ! Generate the inflow at the end of the domain
     iloc_inf=imax
  else
     ! Generate the inflow elsewhere because of convective outflow
     if (inflow_loc.gt.xm(imax)) then
        iloc_inf=imax
     else
        iloc_inf=imin
        do while (xm(iloc_inf).lt.inflow_loc)
           iloc_inf=iloc_inf+1
        end do
     end if
     ! Do not interpolate, just take the closest point => minimize errors
  end if
  
  ! Additional inflow variables
  call parser_is_defined('Additional inflow variables',isdefined)
  if (isdefined) then
     call parser_getsize('Additional inflow variables',nvar)
     allocate(names_tmp(nvar))
     allocate(inf_ind(nvar))
     call parser_read('Additional inflow variables',names_tmp)
  end if

  ! Setup the names
  nvar = nvar + 3
  allocate(names(nvar))
  names(1) = 'U'
  names(2) = 'V'
  names(3) = 'W'
  do n=4,nvar
     names(n) = names_tmp(n-3)
  end do

  ! Find the scalar index
  do n=4,nvar
     loop:do isc=1,nscalar
        if (trim(SC_name(isc)).eq.trim(names(n))) then
           inf_ind(n-3) = isc
           exit loop
        end if
     end do loop
     if (isc.eq.nscalar+1) &
          call die('inflow_generation_init: unknown scalar name in input')
  end do     
  
  ! Allocate arrays
  allocate(buf(jmin_:jmax_,kmin_:kmax_,nvar))
  allocate(fileview(nvar))
  
  ! Open the file to write
  inquire(file=filename,exist=file_is_there)
  filename = trim(mpiiofs) // ":" // trim(filename)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
       MPI_INFO_NULL,ifile,ierr)
  
  ! Set headers
  inflow_save = time
  inflow_ntime = 0
  if (irank.eq.iroot) then
     ! Write dimensions
     dims(1) = inflow_ntime
     dims(2) = ny
     dims(3) = nz
     dims(4) = nvar
     call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
     call MPI_FILE_WRITE(ifile,inflow_freq,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     ! Write variable names
     do var=1,nvar
        call MPI_FILE_WRITE(ifile,names(var),str_short,MPI_CHARACTER,status,ierr)
     end do
     ! Write the grid
     call MPI_FILE_WRITE(ifile,icyl,1,MPI_INTEGER,status,ierr)
     call MPI_FILE_WRITE(ifile,y(jmin:jmax+1),ny+1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,z(kmin:kmax+1),nz+1,MPI_REAL_WP,status,ierr)
  end if
  
  ! Global sizes
  gsizes(1) = ny
  gsizes(2) = nz
  gsizes(3) = nvar
  ! Local sizes
  lsizes(1) = ny_
  lsizes(2) = nz_
  lsizes(3) = 1
  ! Starting index
  start(1) = jmin_-jmin
  start(2) = kmin_-kmin
  ! Define the size
  data_size = lsizes(1)*lsizes(2)*lsizes(3)
  if (.not.(iloc_inf.ge.imin_ .and. iloc_inf.lt.imax_+1)) data_size = 0
  ! Predefine the view  
  do var=1,nvar
     start(3) = var - 1
     call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,&
          MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview(var),ierr)
     call MPI_TYPE_COMMIT(fileview(var),ierr)
  end do
  
  return
end subroutine inflow_generation_init


! ================================ !
! Dump inflow profile if necessary !
! ================================ !
subroutine inflow_generation_save
  use inflow_generation
  use parallel
  use metric_generic
  use data
  use time_info
  implicit none
  
  integer :: var,ierr, tmpint1, tmpint2
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: j,k,n
  
  ! If not return
  if (iinflow.eq.0) return
  
  ! New time step to store
  if (int(time/inflow_freq).eq.inflow_save) return
  inflow_save = int(time/inflow_freq)
  inflow_ntime = inflow_ntime + 1
  
  ! Update the number of time steps stored
  disp = 0
  call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,MPI_INTEGER, &
       "native",MPI_INFO_NULL,ierr)
  if (irank.eq.iroot) then
     call MPI_FILE_WRITE(ifile,inflow_ntime,1,MPI_INTEGER,status,ierr)
  end if
  
  ! Compute all the quantities...
  if (iloc_inf.ge.imin_ .and. iloc_inf.lt.imax_+1) then
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           buf(j,k,1) = sum(interp_u_xm(iloc_inf,j,:)*U(iloc_inf-st1:iloc_inf+st2,j,k))
           buf(j,k,2) = V(iloc_inf,j,k)
           buf(j,k,3) = W(iloc_inf,j,k)
           do n=4,nvar
              buf(j,k,n) = SC(iloc_inf,j,k,inf_ind(n-3))
           end do
        end do
     end do
  end if
  
  ! Displacement
  disp = int(4,8)*int(4,8) + int(2,8)*int(WP,8) + int(nvar,8)*int(str_short,8) + &
            int(4,8) + (int(ny,8)+int(nz,8)+int(2,8))*int(WP,8) + &
            (int(inflow_ntime,8)-int(1,8))*int(nvar,8)*int(ny,8)*int(nz,8)*int(WP,8)

  ! Write the data
  do var=1,nvar
     call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,fileview(var), &
          "native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(ifile,buf(:,:,var),data_size, &
          MPI_REAL_WP,status,ierr)
  end do

  ! Log
  call monitor_log("INFLOW PROFILE FILE WRITTEN")

  return
end subroutine inflow_generation_save


! ============== !
! Close the file !
! ============== !
subroutine inflow_generation_close
  use inflow_generation
  use parallel
  implicit none
  integer :: ierr
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  return
end subroutine inflow_generation_close
