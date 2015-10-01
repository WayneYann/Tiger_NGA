module structure
  use precision
  implicit none
  
  ! Cartesian  : x-direction
  ! Cylindrical: x-direction
  integer :: ibmino
  integer :: ibmin
  integer :: ibmaxo
  integer :: ibmax
  integer :: nbx
  integer :: nbxo
  integer :: ibmin_
  integer :: ibmax_

  ! Cartesian  : y-direction
  ! Cylindrical: r-direction
  integer :: jbmino
  integer :: jbmin
  integer :: jbmaxo
  integer :: jbmax
  integer :: nby
  integer :: nbyo
  integer :: jbmin_
  integer :: jbmax_

  ! Cartesian  : z-direction
  ! Cylindrical: theta-direction
  integer :: kbmino
  integer :: kbmin
  integer :: kbmaxo
  integer :: kbmax
  integer :: nbz
  integer :: nbzo
  integer :: kbmin_
  integer :: kbmax_
  
  ! Node locations (lower left of the cell)
  integer, dimension(:), pointer :: ib1,ib2
  integer, dimension(:), pointer :: jb1,jb2
  integer, dimension(:), pointer :: kb1,kb2
  
  ! Masks
  ! 0: interior point
  ! 1: wall
  ! 2: boundary
  integer, dimension(:,:,:), pointer :: bmask
  
  ! CPU
  integer, dimension(:,:,:), pointer :: bcpu
  
end module structure

! ======================== !
! Initialize the structure !
! ======================== !
subroutine structure_init()
  use structure
  use masks
  use geometry
  use partition
  use parallel
  implicit none
  
  integer :: i,j,k,ierr,count
  integer, dimension(:), pointer :: border_i,border_i2
  integer, dimension(:), pointer :: border_j,border_j2
  integer, dimension(:), pointer :: border_k,border_k2
  integer, dimension(:,:,:), pointer :: bcpu2
  
  ! Initialize
  ibmino = 1
  ibmin  = ibmino + 1
  jbmino = 1
  jbmin  = jbmino + 1
  kbmino = 1
  kbmin  = kbmino + 1
  
  ! Allocate and initialize
  allocate(border_i(imin:imax),border_i2(imin:imax))
  allocate(border_j(jmin:jmax),border_j2(jmin:jmax))
  allocate(border_k(kmin:kmax),border_k2(kmin:kmax))
  border_i=0
  border_j=0
  border_k=0
  border_i2=0
  border_j2=0
  border_k2=0
  
  ! Identify the walls
  do j=jmin_,jmax_
     do i=imin_+1,imax_
        if (mask(i,j).ne.mask(i-1,j)) then
           border_i2(i) = 1
        end if
     end do
  end do
  do j=jmin_+1,jmax_
     do i=imin_,imax_
        if (mask(i,j).ne.mask(i,j-1)) then
           border_j2(j) = 1
        end if
     end do
  end do
  
  ! Add the interprocessor boundaries
  border_i2(imin_) = 1
  border_j2(jmin_) = 1
  border_k2(kmin_) = 1
  
  ! Communicate that
  call MPI_Allreduce(border_i2,border_i,nx,MPI_INTEGER,MPI_MAX,comm,ierr)
  call MPI_Allreduce(border_j2,border_j,ny,MPI_INTEGER,MPI_MAX,comm,ierr)
  call MPI_Allreduce(border_k2,border_k,nz,MPI_INTEGER,MPI_MAX,comm,ierr)
  
  ! Count the block
  nbx    = sum(border_i)
  nby    = sum(border_j)
  nbz    = sum(border_k)
  nbxo   = nbx + 2
  nbyo   = nby + 2
  nbzo   = nbz + 2
  ibmax  = ibmin + nbx - 1
  jbmax  = jbmin + nby - 1
  kbmax  = kbmin + nbz - 1
  ibmaxo = ibmax + 1
  jbmaxo = jbmax + 1
  kbmaxo = kbmax + 1
  
  ! Allocate the arrays
  allocate(bmask(ibmino:ibmaxo,jbmino:jbmaxo,kbmino:kbmaxo))
  allocate(bcpu(ibmin:ibmax,jbmin:jbmax,kbmin:kbmax))
  allocate(ib1(ibmino:ibmaxo),ib2(ibmino:ibmaxo))
  allocate(jb1(jbmino:jbmaxo),jb2(jbmino:jbmaxo))
  allocate(kb1(kbmino:kbmaxo),kb2(kbmino:kbmaxo))
  bmask=0
  bcpu=0
  ib1=0
  jb1=0
  kb1=0
  ib2=0
  jb2=0
  kb2=0
  
  ! Define the blocks
  ib1(ibmino) = imino
  jb1(jbmino) = jmino
  kb1(kbmino) = kmino
  ib2(ibmino) = imin-1
  jb2(jbmino) = jmin-1
  kb2(kbmino) = kmin-1
  
  ib1(ibmaxo) = imax+1
  jb1(jbmaxo) = jmax+1
  kb1(kbmaxo) = kmax+1
  ib2(ibmaxo) = imaxo
  jb2(jbmaxo) = jmaxo
  kb2(kbmaxo) = kmaxo
  
  count = ibmin
  do i=imin,imax
     if (border_i(i).eq.1) then
        ib1(count) = i
        ib2(count-1) = i-1
        count = count + 1
     end if
  end do
  ib2(ibmax) = imax
  
  count = jbmin
  do j=jmin,jmax
     if (border_j(j).eq.1) then
        jb1(count) = j
        jb2(count-1) = j-1
        count = count + 1
     end if
  end do
  jb2(jbmax) = jmax
  
  count = kbmin
  do k=kmin,kmax
     if (border_k(k).eq.1) then
        kb1(count) = k
        kb2(count-1) = k-1
        count = count + 1
     end if
  end do
  kb2(kbmax) = kmax
  
  ! Get the mask and cpu
  allocate(bcpu2(ibmin:ibmax,jbmin:jbmax,kbmin:kbmax))
  bcpu2 = 0
  bcpu  = 0
  bmask(ibmino,:,:) = 2
  bmask(ibmaxo,:,:) = 2
  bmask(:,jbmino,:) = 2
  bmask(:,jbmaxo,:) = 2
  bmask(:,:,kbmino) = 2
  bmask(:,:,kbmaxo) = 2
  do k=kbmin,kbmax
     do j=jbmin,jbmax
        do i=ibmin,ibmax
           bmask(i,j,k) = mask(ib1(i),jb1(j))
           if (ib1(i).GE.imin_ .AND. ib1(i).LE.imax_) then
              if (jb1(j).GE.jmin_ .AND. jb1(j).LE.jmax_) then
                 if (kb1(k).GE.kmin_ .AND. kb1(k).LE.kmax_) then
                    bcpu2(i,j,k)=irank
                 end if
              end if
           end if
        end do
     end do
  end do
  
  ! Communicate the cpu
  call MPI_Allreduce(bcpu2,bcpu,nbx*nby*nbz,MPI_INTEGER,MPI_MAX,comm,ierr)
  
  ! Get the local structure info
  do i=ibmin,ibmax
     if (ib1(i).EQ.imin_) ibmin_=i
     if (ib2(i).EQ.imax_) ibmax_=i
  end do
  do j=jbmin,jbmax
     if (jb1(j).EQ.jmin_) jbmin_=j
     if (jb2(j).EQ.jmax_) jbmax_=j
  end do
  do k=kbmin,kbmax
     if (kb1(k).EQ.kmin_) kbmin_=k
     if (kb2(k).EQ.kmax_) kbmax_=k
  end do
  
  ! Test
  !do k=kbmin,kbmax
  !   do j=jbmin,jbmax
  !      do i=ibmin,ibmax
  !         ! Print
  !         print*,i,j,k
  !         print*,"----->",bmask(i,j,k),bcpu(i,j,k)
  !         print*,"----->",ib1(i),ib2(i)
  !         print*,"----->",jb1(j),jb2(j)
  !         print*,"----->",kb1(k),kb2(k)
  !      end do
  !   end do
  !end do
  
  return
end subroutine structure_init
