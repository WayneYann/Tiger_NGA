module inflow_file
  use inflow
  use parallel
  implicit none
  
  ! Communicator for reading inflow from file
  integer :: inflow_comm
  integer :: file_inflow
  
  ! Name of inflow file and file id
  character(len=str_medium) :: filename
  integer :: ifile
  
  ! Variables
  integer :: inflow_nvar
  character(len=str_short), dimension(:), pointer :: inflow_name
  real(WP), dimension(:,:,:), pointer :: inflow_prev,inflow_next
  
  ! Axis variables
  real(WP), dimension(:), pointer :: Uaxis,Vaxis,Waxis
  real(WP), dimension(:), pointer :: Gaxis
  real(WP), dimension(:,:), pointer :: SCaxis
  
  ! Mapping between inflow variables and real variables
  integer :: index_u,index_v,index_w
  integer :: index_ur,index_vr,index_wr
  integer :: index_ui,index_vi,index_wi
  integer :: index_G
  integer, dimension(:), pointer :: index_sc
  
  ! Grid from the inflow file
  integer :: inflow_icyl
  integer :: inflow_ny,inflow_nz
  real(WP), dimension(:), pointer :: grid_y,grid_z
  real(WP), dimension(:), pointer :: grid_ym,grid_zm
  real(WP) :: inflow_Lz
  
  ! Size of random phase step for mixing layer
  real(WP) :: walk_step,phase
  
  ! Number and size of timesteps
  integer :: inflow_ntime
  real(WP) :: inflow_freq
  real(WP) :: inflow_time
  
  ! Current timestep index and time interpolation coeff
  integer :: itime
  real(WP) :: interp_time
  
  ! Interpolation coefficients in space
  integer,  dimension(:,:), pointer :: index_cc_y,index_cc_z
  integer,  dimension(:,:), pointer :: index_fy_y,index_fy_z
  integer,  dimension(:,:), pointer :: index_fz_y,index_fz_z
  real(WP), dimension(:,:), pointer :: interp_cc_y,interp_cc_z
  real(WP), dimension(:,:), pointer :: interp_fy_y,interp_fy_z
  real(WP), dimension(:,:), pointer :: interp_fz_y,interp_fz_z
  
contains
  
  ! Read the corresponding time from the file
  subroutine inflow_file_read
    use time_info
    implicit none
    
    integer :: itime_new,var,ierr,offset
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(WP) :: tmp
    
    ! Get the corresponding time in the inflow modulo...
    tmp = mod(time/inflow_freq,real(inflow_ntime,WP))
    if (tmp.lt.0.0_WP) tmp = tmp + real(inflow_ntime,WP)
    
    ! Get the index of the timestep right after time
    itime_new = int(tmp) + 1
    
    ! Get the coeff of interpolation in time
    interp_time = real(itime_new,WP) - tmp
    
    ! Nothing to do if still same timestep index
    if (itime_new.eq.itime) return
    
    ! Read one or two time steps?
    if (itime_new.ne.mod(itime,inflow_ntime)+1) then
       ! Compute displacement
       offset = itime_new-2
       if (offset.lt.0) offset = inflow_ntime-1
       disp = int(4,MPI_OFFSET_KIND)*int(4,MPI_OFFSET_KIND) + int(2,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND) + &
            int(inflow_nvar,MPI_OFFSET_KIND)*int(str_short,MPI_OFFSET_KIND) + int(4,MPI_OFFSET_KIND) + &
            (int(inflow_ny,MPI_OFFSET_KIND)+int(inflow_nz,MPI_OFFSET_KIND)+int(2,MPI_OFFSET_KIND))*int(WP,MPI_OFFSET_KIND) + &
            int(offset,MPI_OFFSET_KIND)*int(inflow_nvar,MPI_OFFSET_KIND)*int(inflow_ny,MPI_OFFSET_KIND)*int(inflow_nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
       call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_REAL_WP,"native",MPI_INFO_NULL,ierr)
       
       ! Read the data
       do var=1,inflow_nvar
          call MPI_FILE_READ_ALL(ifile,inflow_next(:,:,var),inflow_ny*inflow_nz,MPI_REAL_WP,status,ierr)
       end do
    end if
    
    ! Compute next displacement
    if (itime_new.ne.itime+1) then
       offset = itime_new-1
       disp = int(4,MPI_OFFSET_KIND)*int(4,MPI_OFFSET_KIND) + int(2,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND) + &
            int(inflow_nvar,MPI_OFFSET_KIND)*int(str_short,MPI_OFFSET_KIND) + int(4,MPI_OFFSET_KIND) + &
            (int(inflow_ny,MPI_OFFSET_KIND)+int(inflow_nz,MPI_OFFSET_KIND)+int(2,MPI_OFFSET_KIND))*int(WP,MPI_OFFSET_KIND) + &
            int(offset,MPI_OFFSET_KIND)*int(inflow_nvar,MPI_OFFSET_KIND)*int(inflow_ny,MPI_OFFSET_KIND)*int(inflow_nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
       call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_REAL_WP,"native",MPI_INFO_NULL,ierr)
    end if
    
    ! Read the data
    itime = itime_new
    inflow_prev = inflow_next
    do var=1,inflow_nvar
       call MPI_FILE_READ_ALL(ifile,inflow_next(:,:,var),inflow_ny*inflow_nz,MPI_REAL_WP,status,ierr)
    end do
    
    return
  end subroutine inflow_file_read
  
  
  ! Compute the interpolation coefficient for one given point
  ! ---------------------------------------------------------
  subroutine inflow_file_get_interp(yin,zin,yloc,zloc,j,k,wj,wk)
    use math
    implicit none
    
    real(WP),         intent(in)  :: yin,zin
    character(len=*), intent(in)  :: yloc,zloc
    real(WP),         intent(out) :: wj,wk
    integer,          intent(out) :: j,k
    
    real(WP), dimension(:), pointer :: mesh_y
    real(WP), dimension(:), pointer :: mesh_z
    real(WP) :: yp,zp,zi
    
    ! Prepare the mesh
    if (trim(yloc).eq.'ym') then
       mesh_y=>grid_ym
    else
       mesh_y=>grid_y
    end if
    if (trim(zloc).eq.'zm') then
       mesh_z=>grid_zm
    else
       mesh_z=>grid_z
    end if
    
    ! Get the points in the inflow config
    if (inflow_icyl.eq.1) then
       if (icyl.eq.1) then
          if (yin.lt.0.0_WP) then
             yp = -yin
             zp = zin+Pi
          else
             yp = yin
             zp = zin
          end if
       else
          yp = sqrt(yin**2+zin**2)
          if (zin/(yp+epsilon(1.0_WP)).ge.0.0_WP) then
             zp = acos(yin/(yp+epsilon(1.0_WP)))
          else
             zp = twoPi-acos(yin/(yp+epsilon(1.0_WP)))
          end if
       end if
    else
       if (icyl.eq.1) then
          yp = yin*cos(zin)
          zp = yin*sin(zin)
       else
          yp = yin
          zp = zin
       end if
    end if
    
    ! Take care of periodicity in inflow
    ! z always periodic
    zp = mod(zp-grid_z(1),inflow_Lz)+grid_z(1)
    !if (zp.lt.grid_z(1)) zp = zp+inflow_Lz
    
    ! Find the point in y
    if (yp.lt.mesh_y(1)) then
       j = 1
       wj = 1.0_WP
    else if (yp.ge.mesh_y(inflow_ny)) then
       j = inflow_ny-1
       wj = 0.0_WP
    else
       loop1:do j=1,inflow_ny-1
          if (yp.lt.mesh_y(j+1)) exit loop1
       end do loop1
       wj = (mesh_y(j+1)-yp)/(mesh_y(j+1)-mesh_y(j))
    end if
    
    ! Find the point in z
    if (zp.lt.mesh_z(1)) then
       k  = inflow_nz
       zi = mesh_z(inflow_nz)-inflow_Lz
       wk = (mesh_z(1)-zp)/(mesh_z(1)-zi)
    else if (zp.ge.mesh_z(inflow_nz)) then
       k  = inflow_nz
       zi = mesh_z(1)+inflow_Lz
       wk = (zi-zp)/(zi-mesh_z(inflow_nz))
    else
       loop3:do k=1,inflow_nz-1
          if (zp.lt.mesh_z(k+1)) exit loop3
       end do loop3
       wk = (mesh_z(k+1)-zp)/(mesh_z(k+1)-mesh_z(k))
    end if
    
    return
  end subroutine inflow_file_get_interp
  
  
  ! Precompute the interpolation coefficients
  ! Include the ghost cells
  ! HERE axis is WRONG. Done later in "inflow_file_precompute"
  subroutine inflow_file_interp
    implicit none
    
    integer  :: j,k,nflow
    
    ! Different types of inflows
    do nflow=1,ninlet
       if (trim(inlet_type(nflow)).eq.'file') then
          
          ! Get the index and coeff
          do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
             do k=kmino_,kmaxo_
                
                ! Cell centers
                call inflow_file_get_interp(ym(j),zm(k),'ym','zm', &
                     index_cc_y(j,k),index_cc_z(j,k),interp_cc_y(j,k),interp_cc_z(j,k))
                
                ! Face in y
                call inflow_file_get_interp(y(j) ,zm(k),'y' ,'zm', &
                     index_fy_y(j,k),index_fy_z(j,k),interp_fy_y(j,k),interp_fy_z(j,k))
                
                ! Face in z
                call inflow_file_get_interp(ym(j),z(k) ,'ym','z' , &
                     index_fz_y(j,k),index_fz_z(j,k),interp_fz_y(j,k),interp_fz_z(j,k))
                
             end do
          end do
       end if
    end do
    
    return
  end subroutine inflow_file_interp

end module inflow_file


! ===================== !
! Initialize the module !
! ===================== !
subroutine inflow_file_init
  use inflow_file
  use data
  use parser
  use math
  implicit none
  
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer :: ierr,var,j,k,isc,nflow
  
  ! Split the communicator
  if (iproc.ne.1) then
     file_inflow = 0
     call MPI_COMM_SPLIT(comm,file_inflow,MPI_UNDEFINED,inflow_comm,ierr)
     return
  end if
  
  ! Does the processor have an type 'file' inlet?
  ! Split the communicator
  file_inflow = 0
  do nflow=1,ninlet
     if (trim(inlet_type(nflow)).eq.'file') then
        if (min(jmax_,inlet(nflow)%jmax)-max(jmin_,inlet(nflow)%jmin).gt.0) file_inflow = 1
     end if
  end do
  call MPI_COMM_SPLIT(comm,file_inflow,MPI_UNDEFINED,inflow_comm,ierr)
  
  ! Nothing to be done if no 'file' infow type
  if (file_inflow.ne.1) return

  ! Open the file to read
  call parser_read('Inflow file to read',filename)
  filename = trim(mpiiofs) // ":" // trim(filename)
  call MPI_FILE_OPEN(inflow_comm,filename,MPI_MODE_RDONLY,mpi_info,ifile,ierr)
  
  ! Read dimensions
  call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(ifile,inflow_freq,1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(ifile,inflow_time,1,MPI_REAL_WP,status,ierr)
  inflow_ntime = dims(1)
  inflow_ny    = dims(2)
  inflow_nz    = dims(3)
  inflow_nvar  = dims(4)
  
  ! Read variable names
  allocate(inflow_name(inflow_nvar))
  allocate(index_sc(inflow_nvar))
  index_u = 0
  index_v = 0
  index_w = 0
  index_G = 0
  index_sc = 0
  do var=1,inflow_nvar
     call MPI_FILE_READ_ALL(ifile,inflow_name(var),str_short,MPI_CHARACTER,status,ierr)
     ! Detect the name of the variable
     select case(trim(inflow_name(var)))
     case('U')
        index_u = var
     case('V')
        index_v = var
     case('W')
        index_w = var
     case('Ur')
        index_ur = var
     case('Vr')
        index_vr = var
     case('Wr')
        index_wr = var
     case('Ui')
        index_ui = var
     case('Vi')
        index_vi = var
     case('Wi')
        index_wi = var
     case('G')
        index_G  = var
     case default
        loop:do isc=1,nscalar
           if (trim(sc_name(isc)).eq.trim(inflow_name(var))) exit loop
        end do loop
        if (isc.eq.nscalar+1) then
           call die('inflow_file_init: Wrong scalar name in inflow file')
        else
           index_sc(var) = isc
        end if
     end select
  end do
  
  ! Check the variables and set the initial phase
  phase = 0.0_WP
  call parser_read('Random phase step',walk_step,0.0_WP)
  walk_step = walk_step * twoPi/180.0_WP
  if (walk_step.eq.0.0_WP) then
     if (index_u.eq.0) call die('inflow_file_init: U not present in inflow file')
     if (index_v.eq.0 .and. index_vr.eq.0) call die('inflow_file_init: V or Vr not present in inflow file')
     if (index_w.eq.0 .and. index_wr.eq.0) call die('inflow_file_init: W or Wr not present in inflow file')
  else
     if (index_u.eq.0 .or. index_ur.eq.0 .or. index_ui.eq.0) &
          call die('inflow_file_init: U or Ur or Ui not present in inflow file')
     if (index_vr.eq.0 .or. index_vi.eq.0) call die('inflow_file_init: Vr or Vi not present in inflow file')
     if (index_wr.eq.0 .or. index_wi.eq.0) call die('inflow_file_init: Wr or Wi not present in inflow file')     
  end if
  
  ! Read the grid
  allocate(grid_y (inflow_ny+1))
  allocate(grid_ym(inflow_ny))
  allocate(grid_z (inflow_nz+1))
  allocate(grid_zm(inflow_nz))
  call MPI_FILE_READ_ALL(ifile,inflow_icyl,1,MPI_INTEGER,status,ierr)
  call MPI_FILE_READ_ALL(ifile,grid_y,inflow_ny+1,MPI_REAL_WP,status,ierr)
  call MPI_FILE_READ_ALL(ifile,grid_z,inflow_nz+1,MPI_REAL_WP,status,ierr)
  
  ! Take care of periodicity
  ! Shift to have first point at 0
  ! GB removed - grid_z = grid_z - grid_z(1)
  inflow_Lz = grid_z(inflow_nz+1)-grid_z(1)
  
  ! Create mid points
  do j=1,inflow_ny
     grid_ym(j) = 0.5*(grid_y(j)+grid_y(j+1))
  end do
  do k=1,inflow_nz
     grid_zm(k) = 0.5*(grid_z(k)+grid_z(k+1))
  end do
  
  ! Allocate arrays
  allocate(inflow_prev(inflow_ny,inflow_nz,inflow_nvar))
  allocate(inflow_next(inflow_ny,inflow_nz,inflow_nvar))
  
  allocate(index_cc_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(index_cc_z(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(index_fy_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(index_fy_z(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(index_fz_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(index_fz_z(jmino_:jmaxo_,kmino_:kmaxo_))
  
  allocate(interp_cc_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(interp_cc_z(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(interp_fy_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(interp_fy_z(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(interp_fz_y(jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(interp_fz_z(jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Allocate the axis variables
  allocate(Uaxis(kmin:kmax),Vaxis(kmin:kmax),Waxis(kmin:kmax))
  allocate(Gaxis(kmin:kmax))
  allocate(SCaxis(kmin:kmax,nscalar))
  
  ! Compute the interpolation coefficients
  call inflow_file_interp
  
  ! Read the first two time steps
  itime = -1
  call inflow_file_read
  
  return
end subroutine inflow_file_init


! =========================================================== !
! Precompute the U,V,W fields at the inlet for the given time !
! =========================================================== !
subroutine inflow_file_velocity
  use inflow_file
  use random
  implicit none
  
  real(WP) :: prev,next
  real(WP) :: wj,wk,rnd
  integer  :: j,k,jj,k1,k2
  real(WP) :: umean,wmean,utmp,wtmp
  integer  :: nflow,var,k0,kPi
  real(WP) :: theta
  
  ! Nothing to be done if no 'file' infow type
  ! Read if necessary the inflow file
  if (file_inflow.eq.1) then
     call inflow_file_read
  end if
  
  ! Random number for a random walk
  call random_number(rnd)
  if (rnd.ge.0.5_WP) then
     phase = phase + walk_step
  else
     phase = phase - walk_step
  end if
  
  ! Different types of inflows
  do nflow=1,ninlet
     if (trim(inlet_type(nflow)).eq.'file') then
        !$OMP PARALLEL PRIVATE(theta,jj,k1,k2,wj,wk,prev,next)
        do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
           !$OMP DO
           do k=kmino_,kmaxo_
              
              ! First set to zero
              Uin(j,k) = 0.0_WP
              Vin(j,k) = 0.0_WP
              Win(j,k) = 0.0_WP
              
              do var=1,inflow_nvar
                 
                 ! Set angle to default
                 theta = 0.0_WP
                 
                 select case(trim(inflow_name(var)))
                 case('U')
                    ! Interpolation in space and time for: U
                    jj = index_cc_y(j,k)
                    k1 = index_cc_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_cc_y(j,k)
                    wk = interp_cc_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Uin(j,k) = Uin(j,k) + interp_time*prev + (1.0_WP-interp_time)*next
                    
                 case('Ur')
                    ! Interpolation in space and time for: Ur (real part of fluctuation)
                    jj = index_cc_y(j,k)
                    k1 = index_cc_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_cc_y(j,k)
                    wk = interp_cc_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Uin(j,k) = Uin(j,k) + cos(phase) * (interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case('Ui')
                    ! Interpolation in space and time for: Ui (imag part of fluctuation)
                    jj = index_cc_y(j,k)
                    k1 = index_cc_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_cc_y(j,k)
                    wk = interp_cc_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Uin(j,k) = Uin(j,k) - sin(phase) * (interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case('V')
                    ! Interpolation in space and time for: V
                    jj = index_fy_y(j,k)
                    k1 = index_fy_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fy_y(j,k)
                    wk = interp_fy_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    if ((inflow_icyl.eq.1) .and. (icyl.eq.0)) theta = grid_zm(k1)+(1.0_WP-wk)*inflow_Lz/real(inflow_nz,WP)
                    Vin(j,k) = Vin(j,k) + (interp_time*prev + (1.0_WP-interp_time)*next)*cos(theta)
                    Win(j,k) = Win(j,k) + (interp_time*prev + (1.0_WP-interp_time)*next)*sin(theta)
                    
                 case('Vr')
                    ! Interpolation in space and time for: Vr (real part of fluctuation)
                    jj = index_fy_y(j,k)
                    k1 = index_fy_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fy_y(j,k)
                    wk = interp_fy_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Vin(j,k) = Vin(j,k) + cos(phase) * (interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case('Vi')
                    ! Interpolation in space and time for: Vi (imag part of fluctuation)
                    jj = index_fy_y(j,k)
                    k1 = index_fy_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fy_y(j,k)
                    wk = interp_fy_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Vin(j,k) = Vin(j,k) - sin(phase) * ( interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case('W')
                    ! Interpolation in space and time for: W
                    jj = index_fz_y(j,k)
                    k1 = index_fz_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fz_y(j,k)
                    wk = interp_fz_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    if ((inflow_icyl.eq.1) .and. (icyl.eq.0)) theta = grid_z(k1)+(1.0_WP-wk)*inflow_Lz/real(inflow_nz,WP)
                    Vin(j,k) = Vin(j,k) - (interp_time*prev + (1.0_WP-interp_time)*next)*sin(theta)
                    Win(j,k) = Win(j,k) + (interp_time*prev + (1.0_WP-interp_time)*next)*cos(theta)
                    
                 case('Wr')
                    ! Interpolation in space and time for: Wr (real part of fluctuation)
                    jj = index_fz_y(j,k)
                    k1 = index_fz_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fz_y(j,k)
                    wk = interp_fz_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Win(j,k) = Win(j,k) + cos(phase) * (interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case('Wi')
                    ! Interpolation in space and time for: Wi (imag part of fluctuation)
                    jj = index_fz_y(j,k)
                    k1 = index_fz_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_fz_y(j,k)
                    wk = interp_fz_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    Win(j,k) = Win(j,k) - sin(phase) * ( interp_time*prev + (1.0_WP-interp_time)*next)
                    
                 case default
                    ! Interpolation in space and time for: scalars
                    ! Done in inflow_file_scalar
                 end select
              end do
           end do
           !$OMP END DO
        end do

        !$OMP END PARALLEL
        
        ! Cylindrical case
        if (icyl.eq.1 .and. max(inlet(nflow)%jmino,jmino_).eq.jmino) then
           
           ! Get values around the axis
           call parallel_gather_dir(Uin(jmin,  kmin_:kmax_),Uaxis(kmin:kmax),'z')
           call parallel_gather_dir(Vin(jmin+1,kmin_:kmax_),Vaxis(kmin:kmax),'z')
           call parallel_gather_dir(Win(jmin,  kmin_:kmax_),Waxis(kmin:kmax),'z')
           
           !$OMP PARALLEL DO PRIVATE(k0,kPi)
           do k=kmino_,kmaxo_
              k0 = kmin+modulo(k-kmin,nz)
              kPi = k0+nz/2
              if (kPi.gt.kmax) kPi = kPi-nz  
              
              Vin(jmin,k) = 0.5_WP*(Vaxis(k0)+Vaxis(kPi))
              do j=1,nover
                 Uin (jmin-j,k)   = + Uaxis(kPi)
                 Vin (jmin-j,k)   = - Vaxis(kPi)
                 Win (jmin-j,k)   = - Waxis(kPi)
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     
     end if

     if (additive_inflow) then
        !$OMP PARALLEL
        do j = max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
           !$OMP DO
           do k=kmino_,kmaxo_
              Uin(j,k) = Uin(j,k) + uadd(nflow)
           end do
           !$OMP END DO
        end do
        !$OMP END PARALLEL
     end if
        
     if (rescale_inflow) then
        ! Compute the mean velocities
        utmp = 0.0_WP
        wtmp = 0.0_WP
        !$OMP PARALLEL DO REDUCTION(+:utmp,wtmp)
        do j = max(inlet(nflow)%jmin,jmin_),min(inlet(nflow)%jmax,jmax_)
           do k=kmin_,kmax_
              utmp = utmp + Uin(j,k)*dA(j)
              wtmp = wtmp + Win(j,k)*dA(j)
           end do
        end do
        !$OMP END PARALLEL DO
        call parallel_sum_dir(utmp,umean,'yz')
        call parallel_sum_dir(wtmp,wmean,'yz')
        umean = umean/inlet(nflow)%A
        wmean = wmean/inlet(nflow)%A
           
        ! Rescale the velocities by keeping constant the turbulent intensity
        !$OMP PARALLEL
        do j = max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
           !$OMP DO
           do k=kmino_,kmaxo_
              Uin(j,k) = Uin(j,k) * ubulk(nflow)/(umean+epsilon(umean))
              Win(j,k) = Win(j,k) * wbulk(nflow)/(wmean+epsilon(wmean))
           end do
           !$OMP END DO
        end do
        !$OMP END PARALLEL
     end if
        
  end do
  
  return
end subroutine inflow_file_velocity


! ============================================================ !
! Precompute the scalar fields at the inlet for the given time !
! ============================================================ !
subroutine inflow_file_scalar
  use inflow_file
  implicit none
  
  real(WP) :: prev,next
  real(WP) :: wj,wk
  integer  :: j,k,jj,k1,k2
  integer  :: nflow,var,isc,k0,kPi
  
  ! Nothing to be done if no 'file' infow type
  if (file_inflow.ne.1) return
  
  ! Read if necessary the inflow file
  call inflow_file_read
  
  ! Different types of inflows
  do nflow=1,ninlet
     if (trim(inlet_type(nflow)).eq.'file') then
        !$OMP PARALLEL PRIVATE(isc,jj,k1,k2,wj,wk,prev,next)
        do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
           !$OMP DO
           do k=kmino_,kmaxo_
              
              do var=1,inflow_nvar
                 select case(trim(inflow_name(var)))
                 case('U')
                    ! Interpolation in space and time for: U
                    ! Done in inflow_file_velocity
                 case('Ur')
                    ! Interpolation in space and time for: Ur
                    ! Done in inflow_file_velocity
                 case('Ui')
                    ! Interpolation in space and time for: Ui
                    ! Done in inflow_file_velocity
                 case('V')
                    ! Interpolation in space and time for: V
                    ! Done in inflow_file_velocity
                 case('Vr')
                    ! Interpolation in space and time for: Vr
                    ! Done in inflow_file_velocity
                 case('Vi')
                    ! Interpolation in space and time for: Vi
                    ! Done in inflow_file_velocity
                 case('W')
                    ! Interpolation in space and time for: W
                    ! Done in inflow_file_velocity
                 case('Wr')
                    ! Interpolation in space and time for: Wr
                    ! Done in inflow_file_velocity
                 case('Wi')
                    ! Interpolation in space and time for: Wi
                    ! Done in inflow_file_velocity
                 case default
                    ! Interpolation in space and time for: scalars
                    isc = index_sc(var)
                    jj = index_cc_y(j,k)
                    k1 = index_cc_z(j,k); k2=mod(k1,inflow_nz)+1
                    wj = interp_cc_y(j,k)
                    wk = interp_cc_z(j,k)
                    prev =       wk *(wj*inflow_prev(jj,k1,var) + (1.0_WP-wj)*inflow_prev(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_prev(jj,k2,var) + (1.0_WP-wj)*inflow_prev(jj+1,k2,var))
                    next =       wk *(wj*inflow_next(jj,k1,var) + (1.0_WP-wj)*inflow_next(jj+1,k1,var)) + &
                         (1.0_WP-wk)*(wj*inflow_next(jj,k2,var) + (1.0_WP-wj)*inflow_next(jj+1,k2,var))
                    SCin(j,k,isc) = interp_time*prev + (1.0_WP-interp_time)*next
                 end select
              end do
           end do
           !$OMP END DO
        end do

        !$OMP END PARALLEL
        
        ! Cylindrical case
        if (icyl.eq.1 .and. max(inlet(nflow)%jmino,jmino_).eq.jmino) then
           
           do var=1,inflow_nvar
              isc = index_sc(var)
              if (isc.eq.0) cycle
              
              ! Get values around the axis
              call parallel_gather_dir(SCin(jmin,kmin_:kmax_,isc),SCaxis(kmin:kmax,isc),'z')
              
              !$OMP PARALLEL DO PRIVATE(k0,kPi)
              do k=kmino_,kmaxo_
                 k0 = k
                 if (k.lt.kmin) k0 = k+nz
                 if (k.gt.kmax) k0 = k-nz
                 kPi = k0+nz/2
                 if (kPi.gt.kmax) kPi = kPi-nz  
                 
                 do j=1,nover
                    SCin(jmin-j,k,isc) = + SCaxis(kPi,isc)
                 end do
              end do
              !$OMP END PARALLEL DO
              
           end do
        end if
        
     end if
  end do
  
  return
end subroutine inflow_file_scalar
