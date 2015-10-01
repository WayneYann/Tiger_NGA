module stat_3d
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use unstruct
  implicit none
  
  ! SP buffers
  real(SP), dimension(:), pointer :: buffer_hexa
  real(SP), dimension(:), pointer :: buffer_wedge
  real(SP), dimension(:,:), pointer :: buffer3_hexa
  real(SP), dimension(:,:), pointer :: buffer3_wedge
  
  ! Fileviews
  integer :: fileview_hexa
  integer :: fileview_wedge
  integer :: fileview_node
  integer :: fileview_hexa_conn
  integer :: fileview_wedge_conn
  integer, dimension(:), pointer :: map
  
  ! Time of stat
  real(WP) :: stat_time
  
  ! Velocity statistics
  real(WP), dimension(:,:,:),   pointer :: stat_U,stat_V,stat_W
  real(WP), dimension(:,:,:),   pointer :: stat_U2,stat_V2,stat_W2
  real(WP), dimension(:,:,:,:), pointer :: stat_SC
  
  ! Scalar statistics
  integer :: stat_nscalar
  character(len=str_short), dimension(:), pointer :: stat_sc_name
  
end module stat_3d


! =============================================================== !
! Dump Statistics in 3D binary ensight gold data - Initialization !
! =============================================================== !
subroutine stat_3d_init
  use stat_3d
  use parser
  use data
  implicit none
  integer :: ierr,isc
  integer, dimension(:), pointer :: blocklength
  
  ! Create & Start the timer
  call timing_create('stat-3D')
  call timing_start ('stat-3D')
  
  ! Initialize stat time
  stat_time = 0.0_WP
  
  ! Initialize the arrays for Velocity stat
  allocate(stat_U (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(stat_V (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(stat_W (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(stat_U2(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(stat_V2(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(stat_W2(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  stat_U  = 0.0_WP
  stat_V  = 0.0_WP
  stat_W  = 0.0_WP
  stat_U2 = 0.0_WP
  stat_V2 = 0.0_WP
  stat_W2 = 0.0_WP
  
  ! Initialize the arrays for Scalar stat
  stat_nscalar = 2*nscalar
  allocate(stat_sc_name(stat_nscalar))
  allocate(stat_SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,stat_nscalar))
  do isc=1,nscalar
     stat_sc_name(2*isc-1) = trim(SC_name(isc))
     stat_sc_name(2*isc+0) = trim(SC_name(isc))//"2"
  end do
  if (stat_nscalar.gt.0) stat_sc = 0.0_WP
  
  ! Allocate buffers
  allocate(buffer_hexa(ncells_hexa_))
  allocate(buffer_wedge(ncells_wedge_))
  allocate(buffer3_hexa(ncells_hexa_,3))
  allocate(buffer3_wedge(ncells_wedge_,3))

  ! For MPI-1 derived datatypes
  allocate(blocklength(max(nnodes_,ncells_hexa_,ncells_wedge_)))
  
  ! Create the views  - hex
  allocate(map(ncells_hexa_))
  map = hexa_-1
  blocklength = 1
  call MPI_TYPE_INDEXED(ncells_hexa_,blocklength,map,MPI_REAL_SP,fileview_hexa,ierr)
  call MPI_TYPE_COMMIT(fileview_hexa,ierr)
  map = map * 8
  blocklength = 8
  call MPI_TYPE_INDEXED(ncells_hexa_,blocklength,map,MPI_INTEGER,fileview_hexa_conn,ierr)
  call MPI_TYPE_COMMIT(fileview_hexa_conn,ierr)
  deallocate(map)
  ! Create the views  - wedge
  allocate(map(ncells_wedge_))
  map = wedge_-1
  blocklength = 1
  call MPI_TYPE_INDEXED(ncells_wedge_,blocklength,map,MPI_REAL_SP,fileview_wedge,ierr)
  call MPI_TYPE_COMMIT(fileview_wedge,ierr)
  map = map * 6
  blocklength = 6
  call MPI_TYPE_INDEXED(ncells_wedge_,blocklength,map,MPI_INTEGER,fileview_wedge_conn,ierr)
  call MPI_TYPE_COMMIT(fileview_wedge_conn,ierr)
  deallocate(map)
  ! Create the views  - nodes
  allocate(map(nnodes_))
  map = nodes_-1
  blocklength = 1
  call MPI_TYPE_INDEXED(nnodes_,blocklength,map,MPI_INTEGER,fileview_node,ierr)
  call MPI_TYPE_COMMIT(fileview_node,ierr)
  deallocate(blocklength)

  ! Stop the timer
  call timing_stop ('stat-3D')

  ! Write the geometry
  call stat_3d_geometry
  
  ! Write the case file
  call stat_3d_case
  
  return
end subroutine stat_3d_init


! ============================ !
! Compute the statistics in 3D !
! ============================ !
subroutine stat_3d_sample
  use stat_3d
  use data
  use interpolate
  use time_info
  implicit none
  integer :: i,j,k,isc
  
  ! Compute the sum
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           ! Velocity
           stat_U(i,j,k)  = stat_U(i,j,k) + dt*Ui(i,j,k)
           stat_V(i,j,k)  = stat_V(i,j,k) + dt*Vi(i,j,k)
           stat_W(i,j,k)  = stat_W(i,j,k) + dt*Wi(i,j,k)
           
           ! Velocity^2
           stat_U2(i,j,k)  = stat_U2(i,j,k) + dt*Ui(i,j,k)**2
           stat_V2(i,j,k)  = stat_V2(i,j,k) + dt*Vi(i,j,k)**2
           stat_W2(i,j,k)  = stat_W2(i,j,k) + dt*Wi(i,j,k)**2
           
           ! Scalars
           do isc=1,nscalar
              stat_SC(i,j,k,isc*2-1)  = stat_SC(i,j,k,2*isc-1) + dt*SC(i,j,k,isc)
              stat_SC(i,j,k,isc*2+0)  = stat_SC(i,j,k,2*isc+0) + dt*SC(i,j,k,isc)**2
           end do
           
        end do
     end do
  end do
  
  ! Add to the time
  stat_time = stat_time + dt
  
  return
end subroutine stat_3d_sample


! ============================================== !
! Dump 3D Statistics in binary ensight gold data !
! ============================================== !
subroutine stat_3d_write
  use stat_3d
  use data
  use partition
  use string
  use interpolate
  use memory
  implicit none
  character(len=str_short) :: name
  integer :: i,j,k,isc
  
  ! Start the timer
  call timing_start('stat-3D')
  
  ! Velocity
  name = "V"
  if (icyl.eq.0) then
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = stat_U(i,j,k) / real(stat_time,WP)
              tmp2(i,j,k) = stat_V(i,j,k) / real(stat_time,WP)
              tmp3(i,j,k) = stat_W(i,j,k) / real(stat_time,WP)
           end do
        end do
     end do
  else
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = stat_U(i,j,k) / real(stat_time,WP)
              tmp2(i,j,k) = (stat_V(i,j,k) * cos(zm(k)) - stat_W(i,j,k) * sin(zm(k))) / real(stat_time,WP)
              tmp3(i,j,k) = (stat_V(i,j,k) * sin(zm(k)) + stat_W(i,j,k) * cos(zm(k))) / real(stat_time,WP)
           end do
        end do
     end do
  end if
  call stat_3d_vector(tmp1,tmp2,tmp3,name)
  
  ! Velocity^2
  name = "V2"
  if (icyl.eq.0) then
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = stat_U2(i,j,k) / real(stat_time,WP)
              tmp2(i,j,k) = stat_V2(i,j,k) / real(stat_time,WP)
              tmp3(i,j,k) = stat_W2(i,j,k) / real(stat_time,WP)
           end do
        end do
     end do
  else
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = stat_U2(i,j,k) / real(stat_time,WP)
              tmp2(i,j,k) = (stat_V2(i,j,k) * cos(zm(k)) - stat_W2(i,j,k) * sin(zm(k))) / real(stat_time,WP)
              tmp3(i,j,k) = (stat_V2(i,j,k) * sin(zm(k)) + stat_W2(i,j,k) * cos(zm(k))) / real(stat_time,WP)
           end do
        end do
     end do
  end if
  call stat_3d_vector(tmp1,tmp2,tmp3,name)
  
  ! Scalars
  do isc=1,stat_nscalar
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = stat_SC(i,j,k,isc) / real(stat_time,WP)
           end do
        end do
     end do
     call stat_3d_scalar(tmp1,stat_sc_name(isc))
  end do
  
  ! Stop the timer
  call timing_stop('stat-3D')
  
  return
end subroutine stat_3d_write


! ========================================================= !
! Dump 3D statistics in binary ensight gold data - geometry !
! => one processor only - test before                       !
! ========================================================= !
subroutine stat_3d_geometry
  use stat_3d
  implicit none
  integer :: ierr,ibuffer,iunit,i
  character(len=80) :: buffer
  character(len=str_medium) :: file
  real(SP), dimension(:), pointer :: xbuf,ybuf,zbuf
  real(SP) :: max_x,max_y,max_z
  real(SP) :: min_x,min_y,min_z
  logical  :: file_is_there
  integer  :: size
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  ! Get single precision mesh
  allocate(xbuf(nnodes_),ybuf(nnodes_),zbuf(nnodes_))
  if (icyl.eq.0) then
     do i=1,nnodes_
        xbuf(i)=x(unstr2str(i,1))
        ybuf(i)=y(unstr2str(i,2))
        zbuf(i)=z(unstr2str(i,3))
     end do
  else
     do i=1,nnodes_
        xbuf(i)=x(unstr2str(i,1))
        ybuf(i)=y(unstr2str(i,2))*cos(z(unstr2str(i,3)))
        zbuf(i)=y(unstr2str(i,2))*sin(z(unstr2str(i,3)))
     end do
  end if
  min_x=minval(xbuf)
  min_y=minval(ybuf)
  min_z=minval(zbuf)
  max_x=maxval(xbuf)
  max_y=maxval(ybuf)
  max_z=maxval(zbuf)
  
  ! Generate the geometry file in parallel
  file="stat/geometry"
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,iunit,ierr)
  
  ! Write header (only root)
  if (irank.eq.iroot) then
     ! Global header
     buffer = 'C Binary'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'Ensight Gold Geometry File'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'Unstructured Geometry from ARTS'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'node id given'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'element id given'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ! Extents
     buffer = 'extents'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,min_x,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_x,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,min_y,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_y,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,min_z,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_z,size,MPI_REAL_SP,status,ierr)
     ! Part header
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'ARTS 3D geometry'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ! Nodes list
     buffer = 'coordinates'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,nnodes,size,MPI_INTEGER,status,ierr)
  end if
  
  ! Write the node positions
  disp = 752
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,nodes_,nnodes_,MPI_INTEGER,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,xbuf,nnodes_,MPI_REAL_SP,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,ybuf,nnodes_,MPI_REAL_SP,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,zbuf,nnodes_,MPI_REAL_SP,status,ierr)
  
  ! Write the hexa connectivity
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
  if (irank.eq.iroot) then
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,ncells_hexa,size,MPI_INTEGER,status,ierr)
     !do i=1,ncells_hexa
     !   call MPI_FILE_WRITE(iunit,i,size,MPI_INTEGER,status,ierr)
     !end do
  end if
  disp = disp+84
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_hexa_conn,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,hexa_,ncells_hexa_,MPI_INTEGER,status,ierr)
  disp = disp+4*ncells_hexa
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_hexa_conn,"native",MPI_INFO_NULL,ierr)
  size = 8*ncells_hexa_
  call MPI_FILE_WRITE_ALL(iunit,conn_hexa,size,MPI_INTEGER,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = disp+8*4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
        size = 1
        call MPI_FILE_WRITE(iunit,ncells_wedge,size,MPI_INTEGER,status,ierr)
        !do i=1,ncells_wedge
        !   call MPI_FILE_WRITE(iunit,i,size,MPI_INTEGER,status,ierr)
        !end do
     end if
     disp = disp+84
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_wedge_conn,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,wedge_,ncells_wedge_,MPI_INTEGER,status,ierr)
     disp = disp+4*ncells_wedge
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_wedge_conn,"native",MPI_INFO_NULL,ierr)
     size = 6*ncells_wedge_
     call MPI_FILE_WRITE_ALL(iunit,conn_wedge,size,MPI_INTEGER,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)

  ! Deallocate arrays
  deallocate(xbuf,ybuf,zbuf)
  
  return
end subroutine stat_3d_geometry


! ========================================================== !
! Dump 3D statistics in binary ensight gold data - case file !
! ========================================================== !
subroutine stat_3d_case
  use stat_3d
  use data
  use time_info
  implicit none
  integer :: iunit,ierr,i
  character(len=80) :: str
  
  ! Write - Single proc & parallel => only root writes (in ASCII)
  if (irank.eq.iroot) then
     ! Open the file
     iunit = iopen()
     open(iunit,file="stat/case",form="formatted",iostat=ierr,status="REPLACE")
     ! Write the case
     str='FORMAT'
     write(iunit,'(a80)') str
     str='type: ensight gold'
     write(iunit,'(a80)') str
     str='GEOMETRY'
     write(iunit,'(a80)') str
     str='model: geometry'
     write(iunit,'(a80)') str
     str='VARIABLE'
     write(iunit,'(a80)') str
     ! Scalars
     do i=1,nscalar
        str = 'scalar per element: 1 '//trim(SC_name(i))//' '//trim(SC_name(i))//'/'//trim(SC_name(i))
        write(iunit,'(a80)') str
     end do
     ! Velocity
     str='vector per element: 1 V V'
     write(iunit,'(a80)') str
     ! Velocity
     str='vector per element: 1 V2 V2'
     write(iunit,'(a80)') str
     ! Time section
     str='TIME'
     write(iunit,'(a80)') str
     str='time set: 1'
     write(iunit,'(a80)') str
     str='number of steps:'
     write(iunit,'(a16,x,i12)') str,1
     str='filename start number: 0'
     write(iunit,'(a80)') str
     str='filename increment: 1'
     write(iunit,'(a80)') str
     str='time values:'
     write(iunit,'(a12,x,ES12.5,x)') str,time
     ! Close the file
     close(iclose(iunit))
  end if
  
  return
end subroutine stat_3d_case


! ======================================================= !
! Dump 3D statistics in binary ensight gold data - scalar !
! ======================================================= !
subroutine stat_3d_scalar(scalar,name)
  use stat_3d
  use string
  use data
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer,i,j,k,count_wedge_,count_hexa_
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(WP),dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: scalar
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  
  ! Extract the data
  count_wedge_ = 0
  count_hexa_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 count_wedge_ = count_wedge_ + 1
                 buffer_wedge(count_wedge_) = real(scalar(i,j,k),SP)
              else
                 count_hexa_ = count_hexa_ + 1
                 buffer_hexa(count_hexa_) = real(scalar(i,j,k),SP)
              end if
           end if
        end do
     end do
  end do
  
  ! Generate the file
  file = "stat/" // trim(adjustl(name))
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,iunit,ierr)
  
  ! Write header (only root)
  if (irank.eq.iroot) then
     buffer = trim(adjustl(name))
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the file - hexa
  disp = 3*80+4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer_hexa,ncells_hexa_,MPI_REAL_SP,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = 3*80+4+4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     end if
     disp = 3*80+4+80+4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer_wedge,ncells_wedge_,MPI_REAL_SP,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine stat_3d_scalar


! ========================================= !
! Dump 3D binary ensight gold data - vector !
! ========================================= !
subroutine stat_3d_vector(vec1,vec2,vec3,name)
  use stat_3d
  use string
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer,i,j,k,count_wedge_,count_hexa_
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec1
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec2
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec3
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  
  ! Extract the data
  count_wedge_ = 0
  count_hexa_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 count_wedge_ = count_wedge_ + 1
                 buffer3_wedge(count_wedge_,1) = real(vec1(i,j,k),SP)
                 buffer3_wedge(count_wedge_,2) = real(vec2(i,j,k),SP)
                 buffer3_wedge(count_wedge_,3) = real(vec3(i,j,k),SP)
              else
                 count_hexa_ = count_hexa_ + 1
                 buffer3_hexa(count_hexa_,1) = real(vec1(i,j,k),SP)
                 buffer3_hexa(count_hexa_,2) = real(vec2(i,j,k),SP)
                 buffer3_hexa(count_hexa_,3) = real(vec3(i,j,k),SP)
              end if
           end if
        end do
     end do
  end do
  
  ! Generate the file
  file = "stat/" // trim(adjustl(name))
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,iunit,ierr)
  
  ! Write header (only root)
  if (irank.eq.iroot) then
     buffer = trim(adjustl(name))
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the data
  disp = 3*80+4+0*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,1),ncells_hexa_,MPI_REAL_SP,status,ierr)
  disp = 3*80+4+1*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,2),ncells_hexa_,MPI_REAL_SP,status,ierr)
  disp = 3*80+4+2*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,3),ncells_hexa_,MPI_REAL_SP,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = 3*80+4+3*ncells_hexa*4+0*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     end if
     disp = 3*80+4+3*ncells_hexa*4+80+0*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,1),ncells_wedge_,MPI_REAL_SP,status,ierr)
     disp = 3*80+4+3*ncells_hexa*4+80+1*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,2),ncells_wedge_,MPI_REAL_SP,status,ierr)
     disp = 3*80+4+3*ncells_hexa*4+80+2*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,3),ncells_wedge_,MPI_REAL_SP,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine stat_3d_vector

