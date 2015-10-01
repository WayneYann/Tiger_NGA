module data
  use simulation
  use geometry
  use partition
  use parser
  use parallel
  use string
  use precision
  implicit none
  
  ! =========================================== !
  ! STAGGERED representation of variables       !
  !                                             !
  ! U(i,j,k)  -> x(i),ym(j),zm(k)  at t=n       !
  ! V(i,j,k)  -> xm(i),y(j),zm(k)  at t=n       !
  ! W(i,j,k)  -> xm(i),ym(j),z(k)  at t=n       !
  ! RHO(i,j,k)-> xm(i),ym(j),zm(k) at t=n+1/2   !
  ! P(i,j,k)  -> xm(i),ym(j),zm(k) at t=n-1/2   !
  ! VISC(i,j,k) -> xm(i),ym(j),zm(k) at t=n+1/2 !
  ! =========================================== !
  
  ! Velocity
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: rhoU
  real(WP), dimension(:,:,:), pointer :: rhoV
  real(WP), dimension(:,:,:), pointer :: rhoW
  
  ! Pressure variation
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: Pold

  ! Viscosity (mu=rho*nu) - molecular
  real(WP), dimension(:,:,:), pointer :: VISCmol
  
  ! Viscosity (mu=rho*nu) - molecular + turbulent
  real(WP), dimension(:,:,:), pointer :: VISC
  
  ! Density
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: RHOold
  real(WP), dimension(:,:,:), pointer :: dRHO

  ! Velocity -- For scalar transport (thermophoresis, inertial effects, etc.)
  real(WP), dimension(:,:,:,:), pointer :: Ut
  real(WP), dimension(:,:,:,:), pointer :: Vt
  real(WP), dimension(:,:,:,:), pointer :: Wt
  real(WP), dimension(:,:,:,:), pointer :: rhoUt
  real(WP), dimension(:,:,:,:), pointer :: rhoVt
  real(WP), dimension(:,:,:,:), pointer :: rhoWt

  ! Scalars (transported only)
  integer :: nscalar
  real(WP), dimension(:,:,:,:), pointer :: SC
  character(len=str_short), dimension(:), pointer :: SC_name
  
  ! Diffusivity (rho*D) - molecular
  real(WP), dimension(:,:,:,:), pointer :: DIFFmol
  
  ! Diffusivity (rhod*D) - molecular + turbulent
  real(WP), dimension(:,:,:,:), pointer :: DIFF
  
  ! Logical : was the array in the data file?
  logical ::        mom_present   ! momentum : rhoU, rhoV, rhoW
  logical ::        vel_present   ! velocity : U, V, W
  logical ::        rho_present
  logical ::       drho_present
  logical ::       visc_present
  logical ::       diff_present
  
  ! ================================== !
  ! Type definitions for MPI IO        !
  ! for a given variable :             !
  !   - name : name of the variable    !
  !   - var  : pointer to data array   !
  !   - view : MPI IO view on the file !
  ! ================================== !
  type MPI_IO_VAR
     character(len=str_short) :: name
     real(WP), dimension(:,:,:), pointer :: var
     integer :: view
  end type MPI_IO_VAR
  
  ! Array with each variables
  integer :: MPI_IO_NVARS
  type(MPI_IO_VAR), dimension(:), pointer :: MPI_IO_DATA
  
  ! Frequency to write data file
  real(WP) :: data_freq
  integer  :: data_last
  
contains
  
  ! =============================================================== !
  ! Define the variables to read/write and the views to the file    !
  !                                                                 !
  ! To read/write a new variable:                                   !
  !   - make a link to the new data array                           !
  ! =============================================================== !
  subroutine data_mpi_init
    implicit none
    
    integer, dimension(3) :: gsizes, lsizes, start
    integer :: var, ierr, isc
    
    nscalar = 0
    mom_present        = .false.
    vel_present        = .false.
    rho_present        = .false.
    drho_present       = .false.
    visc_present       = .false.
    diff_present       = .false.
    
    ! Count the number of scalars and allocate scalar array
    do var=1,MPI_IO_NVARS
       select case(trim(MPI_IO_DATA(var)%name))
       case ('rhoU')
          mom_present = .true.
       case ('rhoV')
          mom_present = .true.
       case ('rhoW')
          mom_present = .true.
       case ('U')
          vel_present = .true.
       case ('V')
          vel_present = .true.
       case ('W')
          vel_present = .true.
       case ('P')
       case ('RHO')
          rho_present = .true.
       case ('dRHO')
          drho_present = .true.
       case ('VISC')
          visc_present = .true.
       case ('DIFF')
          diff_present = .true.
       case default
          nscalar = nscalar + 1
       end select
    end do
    if (nscalar.ne.0) then
       allocate(SC  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
       allocate(DIFF   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
       allocate(DIFFmol(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
       allocate(SC_name(1:nscalar))
    end if
    
    ! Link the variables
    isc = 1
    do var=1,MPI_IO_NVARS
       select case(trim(MPI_IO_DATA(var)%name))
          ! Momentum
       case ('rhoU')
          MPI_IO_DATA(var)%var => rhoU(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('rhoV')
          MPI_IO_DATA(var)%var => rhoV(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('rhoW')
          MPI_IO_DATA(var)%var => rhoW(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
          ! Velocity
       case ('U')
          MPI_IO_DATA(var)%var => U(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('V')
          MPI_IO_DATA(var)%var => V(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('W')
          MPI_IO_DATA(var)%var => W(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('P')
          MPI_IO_DATA(var)%var => P(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
          ! Density
       case ('RHO')
          MPI_IO_DATA(var)%var => RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('dRHO')
          MPI_IO_DATA(var)%var => dRHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
          ! Viscosity / Diffusivity
       case ('VISC')
          MPI_IO_DATA(var)%var => VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
       case ('DIFF')
          MPI_IO_DATA(var)%var => DIFF(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1)
          ! Scalars
       case default
          MPI_IO_DATA(var)%var => SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
          SC_name(isc) = trim(MPI_IO_DATA(var)%name)
          isc = isc + 1
       end select
    end do
    
    ! Define global(g) and local(l) sizes
    gsizes(1) = nx
    gsizes(2) = ny
    gsizes(3) = nz
    
    lsizes(1) = nx_
    lsizes(2) = ny_
    lsizes(3) = nz_
    
    ! Define starting points
    start(1) = imin_-imin
    start(2) = jmin_-jmin
    start(3) = kmin_-kmin
    
    ! Define the view for each variable
    do var=1,MPI_IO_NVARS
       
       call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,&
            MPI_ORDER_FORTRAN,MPI_REAL_WP,MPI_IO_DATA(var)%view,ierr)
       call MPI_TYPE_COMMIT(MPI_IO_DATA(var)%view,ierr)
       
    end do
    
    return
  end subroutine data_mpi_init
  
  
  ! ================================= !
  ! Read the full 3D file in parallel !
  ! ================================= !
  subroutine data_read
    use simulation
    implicit none
    
    integer :: ifile,ierr,var,data_size
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(4) :: dims
    integer(kind=MPI_Offset_kind) :: disp
    integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
    integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
    integer(kind=MPI_Offset_kind) :: NVARS_MOK
    character(len=str_medium) :: filename

    ! Get the name of the file to read in
    call parser_read('Data file to read',filename)
    filename = trim(mpiiofs)//":" // trim(filename)
    
    ! Open the file
    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
    
    ! Read dimensions from header
    call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
    if ((dims(1).ne.nx) .or. (dims(2).ne.ny) .or. (dims(3).ne.nz)) then
       print*, 'grid = ',nx,ny,nz
       print*, 'data = ',dims(1),dims(2),dims(3)
       call die('The size of the data file does not correspond to the grid file')
    end if
    
    ! Read additional stuff
    call MPI_FILE_READ_ALL(ifile,dt,1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_READ_ALL(ifile,time,1,MPI_REAL_WP,status,ierr)
    
    ! Read the names and set the views to the file
    MPI_IO_NVARS = dims(4)
    allocate(MPI_IO_DATA(MPI_IO_NVARS))
    do var=1,MPI_IO_NVARS
       call MPI_FILE_READ_ALL(ifile,MPI_IO_DATA(var)%name,str_short,MPI_CHARACTER,status,ierr)
    end do
    call data_mpi_init()
    
    ! Size of local arrays
    data_size = nx_*ny_*nz_

    ! Resize some integers so MPI can read even the biggest files
    nx_MOK    = int(nx,          MPI_Offset_kind)
    ny_MOK    = int(ny,          MPI_Offset_kind)
    nz_MOK    = int(nz,          MPI_Offset_kind)
    WP_MOK    = int(WP,          MPI_Offset_kind)
    str_MOK   = int(str_short,   MPI_Offset_kind)
    NVARS_MOK = int(MPI_IO_NVARS,MPI_Offset_kind)
    
    ! Read the data for each variable
    do var=1,MPI_IO_NVARS
       var_MOK = int(var,MPI_Offset_kind)
       disp = 4*4 + str_MOK*NVARS_MOK + 2*WP_MOK + & 
                    nx_MOK*ny_MOK*nz_MOK*WP_MOK*(var_MOK-1)
       call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_IO_DATA(var)%view, &
            "native",mpi_info,ierr)
       call MPI_FILE_READ_ALL(ifile,MPI_IO_DATA(var)%var,data_size, &
            MPI_REAL_WP,status,ierr)
    end do
    
    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)
    
    return
  end subroutine data_read
  
end module data


! ======================== !
! Data Initialization      !
!   - array allocation     !
!   - MPI structure for IO !
! ======================== !
subroutine data_init
  use data
  use simulation
  implicit none

  integer :: i,j,k
  
  ! Allocate arrays
  allocate(U(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(V(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(W(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  allocate(rhoU(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoV(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoW(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  allocate(P(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Pold(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  allocate(VISC   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(VISCmol(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  allocate(RHO   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(RHOold(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(dRHO  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))  
  
  ! Default values
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           U      (i,j,k) = 0.0_WP
           V      (i,j,k) = 0.0_WP
           W      (i,j,k) = 0.0_WP
           rhoU   (i,j,k) = 0.0_WP
           rhoV   (i,j,k) = 0.0_WP
           rhoW   (i,j,k) = 0.0_WP
           P      (i,j,k) = 0.0_WP
           RHO    (i,j,k) = 1.0_WP ! For rho_divide
           RHOold (i,j,k) = 1.0_WP ! For rho_divide
           dRHO   (i,j,k) = 0.0_WP
           VISC   (i,j,k) = 0.0_WP
           VISCmol(i,j,k) = 0.0_WP
           if (nscalar.ne.0) then
              SC     (i,j,k,:) = 0.0_WP
              DIFF   (i,j,k,:) = 0.0_WP
              DIFFmol(i,j,k,:) = 0.0_WP
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Read in the data
  call data_read
  
  ! Frequency to write data file
  call parser_read('Data frequency',data_freq)
  data_last = int(time/data_freq)
  
  return
end subroutine data_init


! ===================================== !
! Test if we need to save the data file !
! ===================================== !
subroutine data_write(flag)
  use data
  use time_info
  implicit none
  logical, intent(in) :: flag
  
  if (int(time/data_freq).ne.data_last .or. flag) then
     data_last = int(time/data_freq)
     call data_write_full3D
  end if
  
  return
end subroutine data_write


! ================================== !
! Write the full 3D file in parallel !
! ================================== !
subroutine data_write_full3D
  use data
  use time_info
  implicit none
  
  integer :: ifile,ierr,var,data_size
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK
  character(len=str_medium) :: filename,buffer
  integer :: overwrite
  logical :: file_is_there
  
  ! Get the name of the file to write to
  call parser_read('Data file to write',filename)
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
     dims(4) = MPI_IO_NVARS
     call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
     ! Write additional stuff
     call MPI_FILE_WRITE(ifile,dt,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     ! Write variable names
     do var=1,MPI_IO_NVARS
        call MPI_FILE_WRITE(ifile,MPI_IO_DATA(var)%name,str_short,MPI_CHARACTER,status,ierr)
     end do
  end if
  
  ! Size of local arrays
  data_size = nx_*ny_*nz_
    
  ! Resize some integers so MPI can read even the biggest files
  nx_MOK    = int(nx,          MPI_Offset_kind)
  ny_MOK    = int(ny,          MPI_Offset_kind)
  nz_MOK    = int(nz,          MPI_Offset_kind)
  WP_MOK    = int(WP,          MPI_Offset_kind)
  str_MOK   = int(str_short,   MPI_Offset_kind)
  NVARS_MOK = int(MPI_IO_NVARS,MPI_Offset_kind)
  
  ! Write the data for each variable
  do var=1,MPI_IO_NVARS
     var_MOK = int(var,MPI_Offset_kind)
     disp = 4*4 + str_MOK*NVARS_MOK + 2*WP_MOK + & 
                  nx_MOK*ny_MOK*nz_MOK*WP_MOK*(var_MOK-1)
     call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_IO_DATA(var)%view, &
          "native",mpi_info,ierr)
     call MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA(var)%var,data_size, &
          MPI_REAL_WP,status,ierr)
  end do
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  ! Log
  call monitor_log("3D DATA WRITTEN")
  
  return
end subroutine data_write_full3D


! ====================================== !
! Operations when simulation is finished !
!   - write out the data field           !
! ====================================== !
subroutine data_finalize
  use data
  implicit none
  
  ! Write the data file
  call data_write_full3D
  
  return
end subroutine data_finalize
