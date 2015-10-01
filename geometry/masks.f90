module masks
  use precision
  implicit none
  
  ! Masks to identify interior/exterior/boundary grid points
  ! Mask values:
  !  0  normal interior point
  !  1  boundary point: wall
  !  2  boundary point: non-wall
  !  3  Dirichlet condition
  ! Component mask values:
  !  0  normal point
  !  1  boundary point: wall
  !  2  boundary point: non-wall
  !  3  Dirichlet condition
  ! Node mask values:
  !  0  normal point
  !  1  node surrounded by ONLY masked cells
  ! Directional distances to wall
  !  0  cell is a wall
  !  n  distance to the first wall
  
  
  ! Masks
  integer, dimension(:,:), pointer :: mask

  ! Component masks
  integer, dimension(:,:), pointer :: mask_u
  integer, dimension(:,:), pointer :: mask_v
  integer, dimension(:,:), pointer :: mask_w

  ! Node masks
  integer, dimension(:,:), pointer :: mask_node

  ! Directional distances to wall
  integer, dimension(:,:), pointer :: dist_xp
  integer, dimension(:,:), pointer :: dist_xm
  integer, dimension(:,:), pointer :: dist_yp
  integer, dimension(:,:), pointer :: dist_ym

end module masks


subroutine masks_init(iunit)
  use masks
  use fileio
  use geometry
  use parallel
  implicit none
  integer :: i,j,n,ierr
  integer, intent(inout) :: iunit
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  ! Allocate the arrays
  allocate(mask(imino:imaxo,jmino:jmaxo))
  
  allocate(mask_u(imino:imaxo,jmino:jmaxo))
  allocate(mask_v(imino:imaxo,jmino:jmaxo))
  allocate(mask_w(imino:imaxo,jmino:jmaxo))
  
  allocate(mask_node(imin:imax+1,jmin:jmax+1))
  
  allocate(dist_xp(imino:imaxo,jmino:jmaxo))
  allocate(dist_xm(imino:imaxo,jmino:jmaxo))
  allocate(dist_yp(imino:imaxo,jmino:jmaxo))
  allocate(dist_ym(imino:imaxo,jmino:jmaxo))
  
  ! Read the masks
  call MPI_FILE_READ_ALL(iunit,mask(imin:imax,jmin:jmax),nx*ny,MPI_INTEGER,status,ierr)
  call MPI_FILE_CLOSE(iunit,ierr)
  
  ! Domain x and y boundaries
  do j=jmino,jmaxo
     if (mask(imin,j).eq.0 .or. mask(imin,j).eq.3) then
        mask(imino:imin-1,j) = 2
     else
        mask(imino:imin-1,j) = mask(imin,j)
     end if
     if (mask(imax,j).eq.0) then
        mask(imax+1:imaxo,j) = 2
     else
        mask(imax+1:imaxo,j) = mask(imax,j)
     end if
  end do
  do i=imino,imaxo
     if (mask(i,jmin).eq.0) then
        mask(i,jmino:jmin-1) = 2
     else
        mask(i,jmino:jmin-1) = mask(i,jmin)
     end if
     if (mask(i,jmax).eq.0) then
        mask(i,jmax+1:jmaxo) = 2
     else
        mask(i,jmax+1:jmaxo) = mask(i,jmax)
     end if
  end do
  
  ! Component masks
  mask_u = 0
  do j=jmino,jmaxo
     ! Walls first
     do i=imino,imaxo-1
        if (mask(i,j).eq.1) mask_u(i:i+1,j) = 1
     end do
     if (mask(imaxo,j).eq.1) mask_u(imaxo,j) = 1
     ! Then the rest if not a wall
     do i=imino+1,imaxo
        if (mask_u(i,j).ne.1) then
           if (mask(i-1,j).eq.3 .or. mask(i,j).eq.3) mask_u(i,j) = 3
           if (mask(i-1,j).eq.2 .or. mask(i,j).eq.2) mask_u(i,j) = 2
        end if
     end do
  end do
  
  mask_v = 0
  do i=imino,imaxo
     ! Walls first
     do j=jmino,jmaxo-1
        if (mask(i,j).eq.1) mask_v(i,j:j+1) = 1
     end do
     if (mask(i,jmaxo).eq.1) mask_v(i,jmaxo) = 1
     ! Then the rest if not a wall
     do j=jmino+1,jmaxo
        if (mask_v(i,j).ne.1) then
           if (mask(i,j-1).eq.3 .or. mask(i,j).eq.3) mask_v(i,j) = 3
           if (mask(i,j-1).eq.2 .or. mask(i,j).eq.2) mask_v(i,j) = 2
        end if
     end do
  end do
  
  mask_w = 0
  do j=jmino,jmaxo
     do i=imino,imaxo
        if (mask(i,j).eq.3) mask_w(i,j) = 3
        if (mask(i,j).eq.2) mask_w(i,j) = 2
        if (mask(i,j).eq.1) mask_w(i,j) = 1
     end do
  end do
  
  ! Node masks
  mask_node = 0
  do j=jmin-1,jmax
     do i=imin-1,imax
        if ( (mask(i,j).eq.1     .or. mask(i,j).eq.3) .and. &
             (mask(i+1,j).eq.1   .or. mask(i+1,j).eq.3) .and. &
             (mask(i,j+1).eq.1   .or. mask(i,j+1).eq.3) .and. &
             (mask(i+1,j+1).eq.1 .or. mask(i+1,j+1).eq.3)) then
           mask_node(i+1,j+1) = 1
        end if
     end do
  end do
  
  ! Directional distances to wall
  dist_xm = nxo
  do j=jmino,jmaxo
     do i=imino,imaxo
        loopxm:do n=0,i-1
           if (mask(i-n,j).eq.1) exit loopxm
        end do loopxm
        if (n.ne.i) dist_xm(i,j) = n
     end do
  end do

  dist_xp = nxo
  do j=jmino,jmaxo
     do i=imino,imaxo
        loopxp:do n=0,nxo-i-1
           if (mask(i+n,j).eq.1) exit loopxp
        end do loopxp
        if (n.ne.nxo-i) dist_xp(i,j) = n
     end do
  end do
  
  dist_ym = nyo
  do j=jmino,jmaxo
     do i=imino,imaxo
        loopym:do n=0,j-1
           if (mask(i,j-n).eq.1) exit loopym
        end do loopym
        if (n.ne.j) dist_ym(i,j) = n
     end do
  end do

  dist_yp = nyo
  do j=jmino,jmaxo
     do i=imino,imaxo
        loopyp:do n=0,nyo-j-1
           if (mask(i,j+n).eq.1) exit loopyp
        end do loopyp
        if (n.ne.nyo-j) dist_yp(i,j) = n
     end do
  end do

  ! Determine the configuration we run
  call masks_get_config
  
  return
end subroutine masks_init


! ======================================= !
! Determine if we run a particular config !
! ======================================= !
subroutine masks_get_config
  use masks
  use geometry
  implicit none
  integer :: upper_wall,lower_wall
  
  ! First set the type of simulation to null
  simu_type = ""
  
  ! Mixing Layer
  if ((icyl.eq.0) .and. (xper.eq.1) .and. (yper.eq.0) &
       .and. maxval(mask(imin:imax,jmin:jmax)).eq.0) then
     simu_type = "mixing layer"
  end if
  
  ! HIT
  if ((xper.eq.1) .and. (yper.eq.1)) then
     simu_type = "hit"
  end if
  
  ! Wall detection
  upper_wall = maxval(abs(1-mask(imin:imax,jmax)))
  lower_wall = maxval(abs(1-mask(imin:imax,jmin)))
  
  ! Pipes and channels
  if ((xper.eq.1) .and. (yper.eq.0) .and. &
       (upper_wall.eq.0) .and. (lower_wall.eq.0) .and. icyl.eq.0) then
     simu_type = "channel"
  end if
  if ((xper.eq.1) .and. (yper.eq.0) .and. &
       (upper_wall.eq.0) .and. icyl.eq.1) then
     simu_type = "pipe"
  end if
  
  ! Boundary layer
  if ( (xper.eq.1) .and. (yper.eq.0) .and. (icyl.eq.0) .and. &
       (lower_wall.eq.0) .and. (upper_wall.eq.1)) then
     simu_type = "boundary layer"
  end if
  
  return
end subroutine masks_get_config
