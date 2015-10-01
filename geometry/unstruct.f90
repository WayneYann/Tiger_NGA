! ------------------------------ !
! Unstructured mesh module       !
! Generates an unstructured mesh !
! equivalent to our mesh         !
! ------------------------------ !
module unstruct
  use partition
  use geometry
  use masks
  use parallel
  implicit none
  
  ! Number of active cells
  integer :: ncells
  integer :: ncells_hexa
  integer :: ncells_wedge
  
  ! Number of active nodes
  integer :: nnodes
  
  ! Conversion structured - unstructured for nodes
  integer, dimension(:,:,:), pointer :: str2unstr ! => global node numbering
  integer, dimension(:,:),   pointer :: unstr2str ! => local node numbering
  
  ! Connectivity
  integer, dimension(:,:), pointer :: conn_hexa
  integer, dimension(:,:), pointer :: conn_wedge
  
  ! Local list of active nodes
  integer :: nnodes_
  integer, dimension(:), pointer :: nodes_
  
  ! Local list of cells
  integer :: ncells_
  integer :: ncells_hexa_
  integer, dimension(:), pointer :: hexa_
  integer :: ncells_wedge_
  integer, dimension(:), pointer :: wedge_
  
end module unstruct

subroutine unstruct_init
  use unstruct
  use math
  implicit none
  integer :: i,j,k,count
  integer :: count_hexa_,count_wedge_
  integer :: i1,i2,j1,j2,k1,k2
  integer :: ierr
  integer, dimension(:), pointer :: nnodes_proc,nhexas_proc,nwedges_proc
  
  ! Count the active cells first
  ncells_hexa_  = 0
  ncells_wedge_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 ncells_wedge_ = ncells_wedge_ + 1
              else
                 ncells_hexa_ = ncells_hexa_ + 1
              end if
           end if
        end do
     end do
  end do
  call parallel_sum(ncells_hexa_ ,ncells_hexa )
  call parallel_sum(ncells_wedge_,ncells_wedge)
  ncells  = ncells_hexa  + ncells_wedge
  ncells_ = ncells_hexa_ + ncells_wedge_
  
  ! Create mapping local -> global
  allocate(nhexas_proc (nproc))
  allocate(nwedges_proc(nproc))
  allocate(hexa_ (ncells_hexa_ ))
  allocate(wedge_(ncells_wedge_))
  call MPI_allgather(ncells_hexa_, 1,MPI_INTEGER,nhexas_proc, 1,MPI_INTEGER,comm,ierr)
  call MPI_allgather(ncells_wedge_,1,MPI_INTEGER,nwedges_proc,1,MPI_INTEGER,comm,ierr)
  do i=1,ncells_hexa_
     hexa_(i)  = i + sum(nhexas_proc (1:irank-1))
  end do
  do i=1,ncells_wedge_
     wedge_(i) = i + sum(nwedges_proc(1:irank-1))
  end do
  
  ! Count the active local nodes
  nnodes_ = 0
  k1 = kmin_
  k2 = kmax_
  if (kproc.eq.npz) then
     if (icyl.eq.0 .or. abs(zL-twopi).ge.dz) then
        k2 = k2 + 1
     end if
  end if
  j1 = jmin_
  if (icyl.eq.1 .and. jproc.eq.1) j1 = j1 + 1 ! Count the centerline once if cylindrical
  j2 = jmax_
  if (jproc.eq.npy) j2 = j2 + 1  ! Add the last point to close the domain if last proc
  i1 = imin_
  i2 = imax_
  if (iproc.eq.npx) i2 = i2 + 1  ! Add the last point to close the domain if last proc
  do k=k1,k2
     do j=j1,j2
        do i=i1,i2
           if (mask_node(i,j).eq.0) then
              nnodes_ = nnodes_ + 1
           end if
        end do
     end do
  end do
  if (icyl.eq.1 .and. jproc.eq.1) then ! Now count only once the centerline if cylindrical
     do i=i1,i2
        if (mask_node(i,jmin).eq.0) then
           nnodes_ = nnodes_ + 1
        end if
     end do
  end if
  call parallel_sum(nnodes_,nnodes)
  
  ! Create mapping local -> global
  allocate(nnodes_proc(nproc))
  allocate(nodes_(nnodes_))
  call MPI_allgather(nnodes_,1,MPI_INTEGER,nnodes_proc,1,MPI_INTEGER,comm,ierr)
  do i=1,nnodes_
     nodes_(i) = i + sum(nnodes_proc(1:irank-1))
  end do
  
  ! Allocate the arrays
  allocate(str2unstr(imin_:imax_+1,jmin_:jmax_+1,kmin_:kmax_+1))
  allocate(unstr2str(nnodes_,3))
  allocate(conn_hexa (8,ncells_hexa_ ))
  allocate(conn_wedge(6,ncells_wedge_))
  
  ! Compute str2unstr and unstr2str
  count = 0
  do k=k1,k2
     do j=j1,j2
        do i=i1,i2
           if (mask_node(i,j).eq.0) then
              count=count+1
              str2unstr(i,j,k)=nodes_(count)
              unstr2str(count,1)=i
              unstr2str(count,2)=j
              unstr2str(count,3)=k
           end if
        end do
     end do
  end do
  if (icyl.eq.1 .and. jproc.eq.1) then ! Now count only once the centerline if cylindrical
     do i=i1,i2
        if (mask_node(i,jmin).eq.0) then
           count=count+1
           str2unstr(i,jmin,:)=nodes_(count)
           unstr2str(count,1)=i
           unstr2str(count,2)=jmin
           unstr2str(count,3)=kmin
        end if
     end do
  end if
  call unstruct_comm
  
  ! Compute conn
  count_hexa_  = 0
  count_wedge_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 count_wedge_ = count_wedge_ + 1
                 conn_wedge(1,count_wedge_) = str2unstr(i  ,j+1,k  )
                 conn_wedge(2,count_wedge_) = str2unstr(i  ,j  ,k  )
                 conn_wedge(3,count_wedge_) = str2unstr(i  ,j+1,k+1)
                 conn_wedge(4,count_wedge_) = str2unstr(i+1,j+1,k  )
                 conn_wedge(5,count_wedge_) = str2unstr(i+1,j  ,k  )
                 conn_wedge(6,count_wedge_) = str2unstr(i+1,j+1,k+1)
              else
                 count_hexa_ = count_hexa_ + 1
                 conn_hexa(1,count_hexa_) = str2unstr(i  ,j  ,k  )
                 conn_hexa(2,count_hexa_) = str2unstr(i  ,j  ,k+1)
                 conn_hexa(3,count_hexa_) = str2unstr(i  ,j+1,k+1)
                 conn_hexa(4,count_hexa_) = str2unstr(i  ,j+1,k  )
                 conn_hexa(5,count_hexa_) = str2unstr(i+1,j  ,k  )
                 conn_hexa(6,count_hexa_) = str2unstr(i+1,j  ,k+1)
                 conn_hexa(7,count_hexa_) = str2unstr(i+1,j+1,k+1)
                 conn_hexa(8,count_hexa_) = str2unstr(i+1,j+1,k  )
              end if
           end if
        end do
     end do
  end do
  
  return
end subroutine unstruct_init


subroutine unstruct_comm
  use unstruct
  use math
  implicit none
  
  ! Local variables
  integer, dimension(:,:), pointer :: buf1
  integer, dimension(:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr,k
  
  ! Centerline ---------------------------------------------
  if (icyl.eq.1 .and. jproc.eq.1) then
     do k=kmin_+1,kmax_+1
        str2unstr(:,jmin,k) = str2unstr(:,jmin,kmin)
     end do
  end if
  
  ! X direction --------------------------------------------
  ! Initialize buffer
  allocate(buf1(ny_+1,nz_+1))
  allocate(buf2(ny_+1,nz_+1))
  icount = (ny_+1)*(nz_+1)
  ! Copy left buffer
  buf1(:,:) = str2unstr(imin_,:,:)
  ! Send left buffer to left neighbour
  call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL .and. iproc.ne.npx) str2unstr(imax_+1,:,:) = buf2(:,:)
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  ! Y direction --------------------------------------------
  ! Initialize buffer
  allocate(buf1(nx_+1,nz_+1))
  allocate(buf2(nx_+1,nz_+1))
  icount = (nx_+1)*(nz_+1)
  ! Copy left buffer
  buf1(:,:) = str2unstr(:,jmin_,:)
  ! Send lower buffer to lower neighbour
  call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  ! Paste lower buffer to top
  if (isource.NE.MPI_PROC_NULL .and. jproc.ne.npy) str2unstr(:,jmax_+1,:) = buf2(:,:)
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  ! Z direction --------------------------------------------
  ! Initialize buffer
  allocate(buf1(nx_+1,ny_+1))
  allocate(buf2(nx_+1,ny_+1))
  icount = (nx_+1)*(ny_+1)
  ! Copy left buffer
  buf1(:,:) = str2unstr(:,:,kmin_)
  ! Send left buffer to lower neighbour
  call MPI_CART_SHIFT(comm,2,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) then
     if (kproc.ne.npz) then
        str2unstr(:,:,kmax_+1) = buf2(:,:)
     else
        if (icyl.eq.1 .and. abs(zL-twopi).lt.dz) then
           str2unstr(:,:,kmax_+1) = buf2(:,:)
        end if
     end if
  end if
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine unstruct_comm
