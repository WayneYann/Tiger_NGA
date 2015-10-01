module dump_plot3D
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  implicit none
  
  ! Time info
  integer :: nout_time
  real(WP), dimension(:), allocatable :: out_time
  
  ! File numbers, fileviews and size of data to dump
  integer :: nfile_2D, nfile_3D
  integer,dimension(:), pointer :: fileview_2D1, fileview_2D2, fileview_3D
  integer :: data_size_2D1, data_size_2D2, data_size_3D
  
  ! SP buffers
  real(SP), dimension(:,:,:), pointer :: buffer1_SP
  real(SP), dimension(:,:,:,:), pointer :: buffer3_SP

  ! Number of variables
  integer :: nvar_2D,nvar_3D
  
contains

  ! =============================================== !
  ! Dump a scalar to a plot3D 2D data file          !
  ! =============================================== !
  subroutine write_2D_scalar(iunit,fileview,data_size,scalar_SP)
    use data
    use string
    implicit none
    
    integer, intent(inout) :: iunit
    integer, intent(in) :: fileview, data_size
    real(SP), dimension(imin_:imax_,jmin_:jmax_,1), intent(in) :: scalar_SP
    
    integer(kind=MPI_OFFSET_KIND), parameter :: disp = 12
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: ierr
    
    ! Return if no data to write
    if (data_size.eq.0) return

    call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",MPI_INFO_NULL,ierr)
    call MPI_FILE_WRITE_ALL(iunit,scalar_SP,data_size,MPI_REAL_SP,status,ierr)
    
    return 
  end subroutine write_2D_scalar

  ! =============================================== !
  ! Dump a scalar to a plot3D 3D data file          !
  ! =============================================== !
  subroutine write_3D_scalar(iunit,fileview,data_size,scalar_SP)
    use data
    use string
    implicit none
    
    integer, intent(inout) :: iunit
    integer, intent(in) :: fileview, data_size
    real(SP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: scalar_SP
    
    integer(kind=MPI_OFFSET_KIND), parameter :: disp = 16
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: ierr
    
    ! Return if no data to write
    if (data_size.eq.0) return

    call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",MPI_INFO_NULL,ierr)
    call MPI_FILE_WRITE_ALL(iunit,scalar_SP,data_size,MPI_REAL_SP,status,ierr)
    
    return 
  end subroutine write_3D_scalar

end module dump_plot3D


! ----------------------------------- 2D -------------------------------------

! ================================================ !
! Dump 2D binary plot3D data - Initialization      !
! ================================================ !
subroutine dump_plot3D_2D_init
  use dump_plot3D
  use parallel
  use combustion
  use fileio
  implicit none

  integer :: iunit,isc,var,ierr
  integer,  dimension(3) :: gsizes,lsizes,start,start2
  
  ! Allocate buffers
  if (associated(buffer1_SP)) call die('dump_plot3D_2D_init : Cannot have both 3D and 2D plot3D output')
  allocate(buffer1_SP(imin_:imax_,jmin_:jmax_,1  ))
  allocate(buffer3_SP(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3))
  
  ! CHANGE HERE
  ! For now restart the files from 0
  nfile_2D = 0

  ! Number of variables
  nvar_2D = 5
  if (trim(chemistry).ne.'none') nvar_2D = nvar_2D + 1       ! Density
  if (nscalar.ge.1)              nvar_2D = nvar_2D + 1       ! Diffusivity of first scalar
  if (nscalar.ge.1)              nvar_2D = nvar_2D + nscalar ! The scalars
  
  if (irank.eq.iroot) then
     ! Create directory
     !call CREATE_FOLDER("plot3D-2D")
     call system("mkdir -p plot3D-2D")
     ! Write the geometry
     call dump_plot3D_2D_geometry
     
     ! Write a names file for fieldview
     iunit = iopen()
     open(iunit,file='plot3D-2D/data.nam',form='formatted',status='replace')
     if (trim(chemistry).ne.'none') write(iunit,'(A)') 'RHO'
     write(iunit,'(A)') 'U'
     write(iunit,'(A)') 'V'
     write(iunit,'(A)') 'W'
     write(iunit,'(A)') 'P'
     write(iunit,'(A)') 'VISC'
     if (nscalar.ge.1)  write(iunit,'(A)') 'DIFF'
     do isc=1,nscalar
        write(iunit,'(A)') SC_name(isc)
     end do
     close(iclose(iunit))
  end if

  ! Allocate the views
  allocate(fileview_2D1(nvar_2D))
  allocate(fileview_2D2(nvar_2D))

  if (icyl.eq.0) then
     ! Global sizes
     gsizes(1) = nx
     gsizes(2) = ny
     gsizes(3) = nvar_2D
     ! Local sizes
     lsizes(1) = nx_
     lsizes(2) = ny_
     lsizes(3) = 1
     ! Starting index
     start(1) = imin_-imin
     start(2) = jmin_-jmin
     start(3) = 0
     ! Size of arrays
     if (kproc.eq.1) then
        data_size_2D1 = lsizes(1)*lsizes(2)*lsizes(3)
     else
        data_size_2D1 = 0
     end if
     ! Commit the type
     do var=1,nvar_2D
        start(3) = var-1
        call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview_2D1(var),ierr)
        call MPI_TYPE_COMMIT(fileview_2D1(var),ierr)
     end do
  else
     ! Global sizes
     gsizes(1) = nx
     gsizes(2) = 2*ny
     gsizes(3) = 1
     ! Local sizes
     lsizes(1) = nx_
     lsizes(2) = ny_
     lsizes(3) = 1
     ! Starting index (top)
     start(1)  = imin_-imin
     start(2)  = ny + jmin_-jmin
     start(3)  = 0
     ! Starting index (bottom)
     start2(1)  = imin_-imin
     start2(2)  = ny - (jmin_-jmin) - ny_
     start2(3)  = 0
     ! Size of arrays
     data_size_2D1 = lsizes(1)*lsizes(2)*lsizes(3)
     data_size_2D2 = lsizes(1)*lsizes(2)*lsizes(3)
     ! Commit the type
     do var=1,nvar_2D
        start(3) = var-1
        call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview_2D1(var),ierr)
        call MPI_TYPE_COMMIT(fileview_2D1(var),ierr)
        start2(3) = var-1
        call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start2,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview_2D2(var),ierr)
        call MPI_TYPE_COMMIT(fileview_2D2(var),ierr)
     end do
  end if

  call dump_plot3D_2D

  return
end subroutine dump_plot3D_2D_init


! =========================================== !
! Dump 2D binary plot3D data - geometry       !
! => one processor only - test before         !
! =========================================== !
subroutine dump_plot3D_2D_geometry
  use dump_plot3D
  implicit none

  integer,  dimension(:,:), pointer :: iblank
  real(SP), dimension(:,:), pointer :: grid
  integer  :: iunit,ierr,i,j
  real(SP) :: max_x,max_y,max_z
  real(SP) :: min_x,min_y,min_z
  
  ! Get single precision mesh
  max_x = x(imax+1)
  max_y = y(jmax+1)
  max_z = z(kmax+1)
  min_x = x(imin)
  min_y = y(jmin)
  min_z = z(kmin)
  
  ! Open the grid file to write
  call BINARY_FILE_OPEN(iunit,"plot3D-2D/grid","w",ierr)
     
  ! Cartesian
  if (icyl.eq.0) then
     ! Allocate
     allocate(iblank(nx,ny))
     allocate(grid(nx,ny))

     ! Write the sizes
     call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
 
     ! Write the geometry
     do i=1,nx
        grid(i,:) = xm(i+imin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny,kind(grid),ierr)
     
     do j=1,ny
        grid(:,j) = ym(j+jmin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny,kind(grid),ierr)
     
     iblank(1:nx,1:ny) = 1-mask(imin:imax,jmin:jmax)
     call BINARY_FILE_WRITE(iunit,iblank,nx*ny,kind(iblank),ierr)

  ! Cylindrical
  else
     ! Allocate
     allocate(iblank(nx,2*ny))
     allocate(grid(nx,2*ny))

     ! Write the sizes
     call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit,2*ny,1,kind(ny),ierr)
 
     ! Write the geometry
     do i=1,nx
        grid(i,:) = xm(i+imin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,2*nx*ny,kind(grid),ierr)
     
     do j=1,ny
        grid(:,ny+j)     =  ym(j+jmin-1)
        grid(:,ny-(j-1)) = -ym(j+jmin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,2*nx*ny,kind(grid),ierr)

     ! Get the iblank in 3D
     iblank(1:nx,ny+1:2*ny) = 1-mask(imin:imax,jmin:jmax)
     iblank(1:nx,1:ny)      = 1-mask(imin:imax,jmax:jmin:-1)
     call BINARY_FILE_WRITE(iunit,iblank,2*nx*ny,kind(iblank),ierr)
  end if
     
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  deallocate(grid)
  deallocate(iblank)
  
  return
end subroutine dump_plot3D_2D_geometry


! ========================================= !
! Dump 2D binary data in plot3D format      !
! ========================================= !
subroutine dump_plot3D_2D
  use dump_plot3D
  use data
  use partition
  use string
  use combustion
  use interpolate
  implicit none

  character(len=str_medium) :: filename
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: i,j,k,isc,var
  integer :: ierr, iunit
  integer, dimension(3) :: dims
  logical :: file_is_there
  
  ! Get the file to write - might need to change indexing here
  nfile_2D = nfile_2D + 1
  filename = trim(mpiiofs) // ":plot3D-2D/data" 
  write(filename(len_trim(filename)+1:len_trim(filename)+5),'(i5.5)') nfile_2D
  
  ! Open the file to write
  inquire(file=filename,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                     MPI_INFO_NULL,iunit,ierr)
  
  ! Start wiuth the first variable
  var = 0

  if (icyl.eq.0) then     
     ! Write header
     if (irank.eq.iroot) then          
        ! Write dimensions
        dims(1) = nx
        dims(2) = ny 
        dims(3) = nvar_2D
        call MPI_FILE_WRITE(iunit,dims,3,MPI_INTEGER,status,ierr)
     end if
     
     ! Density
     if (trim(chemistry).ne.'none') then
        var = var + 1
        buffer1_SP(:,:,1) = RHO(imin_:imax_,jmin_:jmax_,kmin)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     end if

     ! U velocity
     var = var + 1
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer1_SP(i,j,1) = Ui(i,j,1)
        end do
     end do
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     
     ! V velocity
     var = var + 1
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer1_SP(i,j,1) = Vi(i,j,1)
        end do
     end do
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     
     ! W velocity
     var = var + 1
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer1_SP(i,j,1) = Wi(i,j,1)
        end do
     end do
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     
     ! Pressure
     var = var + 1
     buffer1_SP(:,:,1) = P(imin_:imax_,jmin_:jmax_,kmin)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     
     ! Viscosity
     var = var + 1
     buffer1_SP(:,:,1) = VISC(imin_:imax_,jmin_:jmax_,kmin)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     
     ! Diffusivity
     if (nscalar.ge.1) then
        var = var + 1
        buffer1_SP(:,:,1) = DIFF(imin_:imax_,jmin_:jmax_,kmin,1)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     end if

     ! Scalars
     do isc=1,nscalar
        var = var + 1
        buffer1_SP(:,:,1) = SC(imin_:imax_,jmin_:jmax_,kmin,isc)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     end do
     
  else
     
     ! Write header
     if (irank.eq.iroot) then          
        ! Write dimensions
        dims(1) = nx
        dims(2) = 2*ny 
        dims(3) = nvar_2D
        call MPI_FILE_WRITE(iunit,dims,3,MPI_INTEGER,status,ierr)
     end if
  
     ! Velocities
     k=kmin
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer3_SP(i,j,k,1) = Ui(i,j,k)
           buffer3_SP(i,j,k,2) = Vi(i,j,k)
           buffer3_SP(i,j,k,3) = Wi(i,j,k)
        end do
     end do
     k=kmin+nz/2
     do j=jmin_,jmax_
        do i=imin_,imax_
           buffer3_SP(i,j,k,1) = + Ui(i,j,k)
           buffer3_SP(i,j,k,2) = - Vi(i,j,k)
           buffer3_SP(i,j,k,3) = - Wi(i,j,k)
        end do
     end do
     
     ! Density
     if (trim(chemistry).ne.'none') then
        var = var + 1
        buffer1_SP(:,:,1) = RHO(imin_:imax_,jmin_:jmax_,kmin)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
        buffer1_SP(:,:,1) = RHO(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2)
        call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
     end if

     ! U velocity
     var = var + 1
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmin_:jmax_,kmin,1)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2,1)
     call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
    
     ! V velocity
     var = var + 1
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmin_:jmax_,kmin,2)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2,2)
     call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)

     ! W velocity
     var = var + 1
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmin_:jmax_,kmin,3)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     buffer1_SP(:,:,1) = buffer3_SP(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2,3)
     call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)

     ! Pressure 
     var = var + 1
     buffer1_SP(:,:,1) = P(imin_:imax_,jmin_:jmax_,kmin)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     buffer1_SP(:,:,1) = P(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2)
     call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
  
     ! Viscosity
     var = var + 1
     buffer1_SP(:,:,1) = VISC(imin_:imax_,jmin_:jmax_,kmin)
     call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
     buffer1_SP(:,:,1) = VISC(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2)
     call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
     
     ! Diffusivity
     if (nscalar.ge.1) then
        var = var + 1
        buffer1_SP(:,:,1) = DIFF(imin_:imax_,jmin_:jmax_,kmin,1)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
        buffer1_SP(:,:,1) = DIFF(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2,1)
        call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
     end if

     ! Scalars
     do isc=1,nscalar
        var = var + 1
        buffer1_SP(:,:,1) = SC(imin_:imax_,jmin_:jmax_,kmin,isc)
        call write_2D_scalar(iunit,fileview_2D1(var),data_size_2D1,buffer1_SP)
        buffer1_SP(:,:,1) = SC(imin_:imax_,jmax_:jmin_:-1,kmin+nz/2,isc)
        call write_2D_scalar(iunit,fileview_2D2(var),data_size_2D2,buffer1_SP)
     end do

  end if
     
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine dump_plot3D_2D


! --------------------------------------- 3D ------------------------------------------


! ================================================ !
! Dump 3D binary plot3D data - Initialization      !
! ================================================ !
subroutine dump_plot3D_3D_init
  use dump_plot3D
  use parallel
  use combustion
  use fileio
  implicit none

  integer :: iunit,isc,var,ierr
  integer,  dimension(4) :: gsizes,lsizes,start
  
  ! Allocate buffers
  if (associated(buffer1_SP)) call die('dump_plot3D_2D_init : Cannot have both 3D and 2D plot3D output')
  if (jproc.eq.1) then
     allocate(buffer1_SP(imin_:imax_,jmin_-icyl:jmax_,kmin_:kmax_+icyl  ))
     allocate(buffer3_SP(imin_:imax_,jmin_-icyl:jmax_,kmin_:kmax_+icyl,3))
  else
     allocate(buffer1_SP(imin_:imax_,jmin_:jmax_,kmin_:kmax_+icyl  ))
     allocate(buffer3_SP(imin_:imax_,jmin_:jmax_,kmin_:kmax_+icyl,3))
  end if
  
  ! CHANGE HERE
  ! For now restart the files from 0
  nfile_3D = 0

  ! Number of variables
  nvar_3D = 5
  if (trim(chemistry).ne.'none') nvar_3D = nvar_3D + 1       ! Density
  if (nscalar.ge.1)              nvar_3D = nvar_3D + 1       ! Diffusivity of first scalar
  if (nscalar.ge.1)              nvar_3D = nvar_3D + nscalar ! The scalars

  if (irank.eq.iroot) then
     ! Create directory
     !call CREATE_FOLDER("plot3D-3D")
     call system("mkdir -p plot3D-3D")
     ! Write the geometry
     call dump_plot3D_3D_geometry
  
     ! Write a name file for fieldview
     iunit = iopen()
     open(iunit,file='plot3D-3D/data.nam',form='formatted',status='replace')
     if (trim(chemistry).ne.'none') write(iunit,'(A)') 'RHO'
     write(iunit,'(A)') 'U'
     write(iunit,'(A)') 'V'
     write(iunit,'(A)') 'W'
     write(iunit,'(A)') 'P'
     write(iunit,'(A)') 'VISC'
     if (nscalar.ge.1)  write(iunit,'(A)') 'DIFF'
     do isc=1,nscalar
        write(iunit,'(A)') SC_name(isc)
     end do
     close(iclose(iunit))
  end if

  ! Allocate fileview
  allocate(fileview_3D(nvar_3D))

  ! Global sizes
  gsizes(1) = nx
  gsizes(2) = ny+icyl
  gsizes(3) = nz+icyl
  gsizes(4) = nvar_3D
  ! Local sizes
  lsizes(1) = nx_
  lsizes(2) = ny_
  if (jproc.eq.1) lsizes(2) = lsizes(2) + icyl
  lsizes(3) = nz_+icyl
  lsizes(4) = 1
  ! Starting index
  start(1) = imin_-imin
  start(2) = jmin_-jmin
  if (jproc.gt.1) start(2) = start(2) + icyl
  start(3) = kmin_-kmin
  start(4) = 0
  ! Size of arrays
  data_size_3D = lsizes(1)*lsizes(2)*lsizes(3)
  ! Commit the type
  do var=1,nvar_3D
     start(4) = var-1
     call MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview_3D(var),ierr)
     call MPI_TYPE_COMMIT(fileview_3D(var),ierr)
  end do
 
  call dump_plot3D_3D

  return
end subroutine dump_plot3D_3D_init


! =========================================== !
! Dump 3D binary plot3D data - geometry       !
! => one processor only - test before         !
! =========================================== !
subroutine dump_plot3D_3D_geometry
  use dump_plot3D
  implicit none

  integer,  dimension(:,:,:), pointer :: iblank
  real(SP), dimension(:,:,:), pointer :: grid
  integer  :: iunit,ierr,i,j,k
  
  ! Open the grid file to write
  call BINARY_FILE_OPEN(iunit,"plot3D-3D/grid","w",ierr)
     
  ! Cartesian
  if (icyl.eq.0) then
     ! Allocate
     allocate(iblank(nx,ny,nz))
     allocate(grid(nx,ny,nz))

     ! Write the sizes
     call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
 
     ! Write the geometry
     do i=1,nx
        grid(i,:,:) = xm(i+imin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     
     do j=1,ny
        grid(:,j,:) = ym(j+jmin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     
     do k=1,nz
        grid(:,:,k) = zm(k+kmin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*ny*nz,kind(grid),ierr)
     
     do k=1,nz
        iblank(1:nx,1:ny,k) = 1-mask(imin:imax,jmin:jmax)
     end do
     call BINARY_FILE_WRITE(iunit,iblank,nx*ny*nz,kind(iblank),ierr)

  ! Cylindrical
  else
     ! Allocate
     allocate(iblank(nx,ny+1,nz+1))
     allocate(grid(nx,ny+1,nz+1))

     ! Write the sizes
     call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit,ny+1,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit,nz+1,1,kind(nz),ierr)
 
     ! Write the geometry
     do i=1,nx
        grid(i,:,:) = xm(i+imin-1)
     end do
     call BINARY_FILE_WRITE(iunit,grid,nx*(ny+1)*(nz+1),kind(grid),ierr)
     
     do k=1,nz
        do j=2,ny+1
           grid(:,j,k) = ym(j+jmin-2) * cos(zm(k+kmin-1))
        end do
     end do
     grid(:,:,nz+1) = grid(:,:,1)
     grid(:,1,:) = 0.0_SP
     call BINARY_FILE_WRITE(iunit,grid,nx*(ny+1)*(nz+1),kind(grid),ierr)
     
     do k=1,nz
        do j=2,ny+1
           grid(:,j,k) = ym(j+jmin-2) * sin(zm(k+kmin-1))
        end do
     end do
     grid(:,:,nz+1) = grid(:,:,1)
     grid(:,1,:) = 0.0_SP
     call BINARY_FILE_WRITE(iunit,grid,nx*(ny+1)*(nz+1),kind(grid),ierr)

     do k=1,nz+1
        iblank(:,2:ny+1,k) = 1-mask(imin:imax,jmin:jmax)
     end do
     iblank(:,1,:) = iblank(:,2,:)
     call BINARY_FILE_WRITE(iunit,iblank,nx*(ny+1)*(nz+1),kind(iblank),ierr)
  end if
     
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  deallocate(grid)
  deallocate(iblank)

  return
end subroutine dump_plot3D_3D_geometry


! ========================================= !
! Dump 3D binary data in plot3D format      !
! ========================================= !
subroutine dump_plot3D_3D
  use dump_plot3D
  use data
  use partition
  use string
  use combustion
  use interpolate
  implicit none

  character(len=str_medium) :: filename
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: i,j,k,isc,var,j0
  integer :: ierr, iunit
  integer, dimension(4) :: dims
  logical :: file_is_there
  
  ! Get the file to write - might need to change indexing here
  nfile_3D = nfile_3D + 1
  filename = trim(mpiiofs) // ":plot3D-3D/data" 
  write(filename(len_trim(filename)+1:len_trim(filename)+5),'(i5.5)') nfile_3D

  ! Open the file to write
  inquire(file=filename,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                     MPI_INFO_NULL,iunit,ierr)

  ! Write header
  if (irank.eq.iroot) then          
     ! Write dimensions
     dims(1) = nx
     dims(2) = ny+icyl
     dims(3) = nz+icyl
     dims(4) = nvar_3D
     call MPI_FILE_WRITE(iunit,dims,4,MPI_INTEGER,status,ierr)
  end if
  var = 0
  
  ! Handling the additional point on the axis
  if (icyl.eq.1 .and. jproc.eq.1) then
     j0 = jmin_-1
  else
     j0 = jmin_
  end if

  ! Velocities
  if (icyl.eq.0) then
     do k=kmin_,kmax_
        do j=j0,jmax_
           do i=imin_,imax_
              buffer3_SP(i,j,k,1) = Ui(i,j,k)
              buffer3_SP(i,j,k,2) = Vi(i,j,k)
              buffer3_SP(i,j,k,3) = Wi(i,j,k)
           end do
        end do
     end do
  else
     do k=kmin_,kmax_+1
        do j=j0,jmax_
           do i=imin_,imax_
              buffer3_SP(i,j,k,1) = Ui(i,j,k)
              buffer3_SP(i,j,k,2) = Vi(i,j,k) * cos(zm(k)) - Wi(i,j,k) * sin(zm(k))
              buffer3_SP(i,j,k,3) = Vi(i,j,k) * sin(zm(k)) + Wi(i,j,k) * cos(zm(k))
           end do
        end do
     end do
  end if
  
  ! Density
  if (trim(chemistry).ne.'none') then
     var = var + 1
     buffer1_SP = RHO(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl)
     call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
  end if

  ! U velocity
  var = var + 1
  buffer1_SP = buffer3_SP(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl,1)
  call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
     
  ! V velocity
  var = var + 1
  buffer1_SP = buffer3_SP(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl,2)
  call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
     
  ! W velocity
  var = var + 1
  buffer1_SP = buffer3_SP(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl,3)
  call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
     
  ! Pressure
  var = var + 1
  buffer1_SP = P(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl)
  call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
     
  ! Viscosity
  var = var + 1
  buffer1_SP = VISC(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl)
  call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
  
  ! Diffusivity
  if (nscalar.ge.1) then
     var = var + 1
     buffer1_SP = DIFF(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl,1)
     call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
  end if
  
  ! Scalars
  do isc=1,nscalar
     var = var + 1
     buffer1_SP = SC(imin_:imax_,j0:jmax_,kmin_:kmax_+icyl,isc)
     call write_3D_scalar(iunit,fileview_3D(var),data_size_3D,buffer1_SP)
  end do

  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine dump_plot3D_3D


