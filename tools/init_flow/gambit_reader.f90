module gambit_reader
  use precision
  use string
  use param
  implicit none
  
  integer :: numnp
  integer :: nxpts, nypts
  integer,  dimension(:,:), pointer :: got_points
  real(WP), dimension(:)  , pointer :: xval,yval
  real(WP), parameter :: thres=1.0e-10_WP
  integer :: open_boundary
  
end module gambit_reader


subroutine gambit_readmesh(bc)
  use gambit_reader
  implicit none
  
  integer, intent(in) :: bc
  open_boundary = bc
  
  ! Read the mesh file
  call read_gambit
  
  ! Read the points
  call read_points
  
  ! Create the mask
  call mesh_mask
  
  return
end subroutine gambit_readmesh


! ============================================= !
! Read the gambit file and store all the points !
! ============================================= !
subroutine read_gambit
  use gambit_reader
  use fileio
  use parser
  implicit none
  
  character(str_medium) :: filename  
  integer :: meshfile,ierr
  integer :: nelem,ngprs,nbsets,ndfcd,ndfvl
  integer :: np
  integer :: dummy
  character*120 :: buffer

  ! Get name of file
  call parser_read('Gambit neutral file',filename)

  ! Open mesh file
  meshfile = iopen()
  open (meshfile, file=filename, status="old", iostat=ierr)
  if (ierr.ne.0) stop "Could not open the mesh file"
  
  ! Read header
  read(meshfile, "(a)", iostat=ierr) buffer ! Header Record Descriptor
  read(meshfile, "(a)", iostat=ierr) buffer ! Neutral-File Header
  if (trim(buffer).ne."** GAMBIT NEUTRAL FILE") stop "Requires a Gambit Neutral File"
  read(meshfile, "(a)", iostat=ierr) buffer ! User-Defined Title
  read(meshfile, "(a)", iostat=ierr) buffer ! Data Source and Revision Level
  read(meshfile, "(a)", iostat=ierr) buffer ! Date and Time
  read(meshfile, "(a)", iostat=ierr) buffer ! Problem Size Parameters Headings
  read(meshfile, "(a)", iostat=ierr) buffer ! Problem Size Parameters 
  read(buffer,"(6I10)") numnp,nelem,ngprs,nbsets,ndfcd,ndfvl
  if (ndfcd.ne.2) stop "Requires a 2D Gambit Neutral File"
  read(meshfile, "(a)", iostat=ierr) buffer ! End of Section
  if (trim(buffer).ne."ENDOFSECTION") stop "ENDOFSECTION expected (header)"

  ! Allocate arrays
  allocate(xval(numnp))
  allocate(yval(numnp))

  ! Read Nodal Coordinates
  read(meshfile, "(a)", iostat=ierr) buffer
  if (trim(buffer(1:20)).ne."   NODAL COORDINATES") stop "NODAL COORDINATES expected"
  do np=1,numnp
     read(meshfile, "(a)", iostat=ierr) buffer
     read(buffer,"(i10,1x,2e20.12)") dummy, xval(np), yval(np)
  end do
  read(meshfile, "(a)", iostat=ierr) buffer
  if (trim(buffer).ne."ENDOFSECTION") stop "ENDOFSECTION expected (nodal coordinates)"

  ! Close file
  close (iclose(meshfile))
  
  return
end subroutine read_gambit


! ============================================ !
! Trim and Sort the points by insertion method !
! ============================================ !
subroutine read_points
  use gambit_reader
  implicit none
  
  integer  :: np,i,j
  real(WP), dimension(:), pointer :: xtmp, ytmp
  integer,  dimension(:), pointer :: xnbr, ynbr
  integer  :: index,itmp
  real(WP) :: alpha,rtmp

  ! Allocate temporary arrays
  allocate(xtmp(numnp))
  allocate(ytmp(numnp))
  allocate(xnbr(numnp))
  allocate(ynbr(numnp))
  nxpts = 0
  nypts = 0
  
  ! Trim
  do np=1,numnp
     ! x - direction
     loopi:do i=1,nxpts
        if (abs(xtmp(i)-xval(np)).le.thres*abs(xval(np))) exit loopi
     end do loopi
     if (i>nxpts) then
        nxpts = nxpts + 1
        xtmp(nxpts) = xval(np)
        xnbr(nxpts) = 1
     else
        xnbr(i) = xnbr(i) + 1
     end if
     
     ! y - direction
     loopj:do j=1,nypts
        if (abs(ytmp(j)-yval(np)).le.thres*abs(yval(np))) exit loopj
     end do loopj
     if (j>nypts) then
        nypts = nypts + 1
        ytmp(nypts) = yval(np)
        ynbr(nypts) = 1
     else
        ynbr(j) = ynbr(j) + 1
     end if
  end do
  
  ! Sort
  do np=1,nxpts-1
     index = np
     do i=np+1,nxpts
        if (xtmp(i)<xtmp(index)) index = i
     end do
     rtmp = xtmp(index)
     xtmp(index) = xtmp(np)
     xtmp(np) = rtmp
     itmp = xnbr(index)
     xnbr(index) = xnbr(np)
     xnbr(np) = itmp
  end do
  do np=1,nypts-1
     index = np
     do j=np+1,nypts
        if (ytmp(j)<ytmp(index)) index = j
     end do
     rtmp = ytmp(index)
     ytmp(index) = ytmp(np)
     ytmp(np) = rtmp
     itmp = ynbr(index)
     ynbr(index) = ynbr(np)
     ynbr(np) = itmp
  end do
  
  ! Check the points
  i = 2
  do while(i<=nxpts-1)
     alpha = (xtmp(i+1)-xtmp(i))/(xtmp(i+1)-xtmp(i-1))
     if (alpha.lt.0.25_WP) then
        write(*,"(a24,ES14.7,a2,i3,a6,ES14.7,a2,i3,a1)"), "Replacing point in x : ", &
             xtmp(i)," (",xnbr(i),") by ",xtmp(i+1)," (",xnbr(i+1),")"
        do np=1,numnp
           if (abs(xval(np)-xtmp(i)).lt.thres*abs(xtmp(i))) xval(np) = xtmp(i+1)
        end do
        xtmp(i:nxpts-1) = xtmp(i+1:nxpts)
        itmp = xnbr(i)
        xnbr(i:nxpts-1) = xnbr(i+1:nxpts)
        xnbr(i) = xnbr(i) + itmp
        nxpts = nxpts-1
     else
        i = i+1
     end if
  end do
  j = 2
  do while (j<=nypts-1)
     alpha = (ytmp(j+1)-ytmp(j))/(ytmp(j+1)-ytmp(j-1))
     if (alpha.lt.0.25_WP) then
        write(*,"(a24,ES14.7,a2,i3,a6,ES14.7,a2,i3,a1)"), "Replacing point in y : ", &
             ytmp(j)," (",ynbr(j),") by ",ytmp(j+1)," (",ynbr(j+1),")"
        do np=1,numnp
           if (abs(yval(np)-ytmp(j)).lt.thres*abs(ytmp(j))) yval(np) = ytmp(j+1)
        end do
        ytmp(j:nypts-1) = ytmp(j+1:nypts)
        itmp = ynbr(j)
        ynbr(j:nypts-1) = ynbr(j+1:nypts)
        ynbr(j) = ynbr(j) + itmp
        nypts = nypts-1
     else
        j = j+1
     end if
  end do
  
  ! Recheck just for fun
  do i=1,nxpts
     if (xnbr(i).ne.nypts) print*,"Warning: uncorrect occurrence for point x ",i," with",xnbr(i)
  end do
  do j=1,nypts
     if (ynbr(j).ne.nxpts) print*,"Warning: uncorrect occurrence for point y ",j," with",ynbr(j)
  end do
  
  ! Allocate the x,y arrays
  nx = nxpts-1
  allocate(x(nx+1))
  x = xtmp(1:nxpts)
  
  if (icyl.eq.1) then
     if (open_boundary.eq.1) then
        ny = nypts-1
        allocate(y(ny+1))
        y = ytmp(1:nypts)
     else
        ny = nypts
        allocate(y(ny+1))
        y(1:ny) = ytmp(1:nypts)
        y(ny+1) = 2.0_WP*y(ny) - y(ny-1)
     end if
  else
     if (open_boundary.eq.1) then
        ny = nypts-1
        allocate(y(ny+1))
        y = ytmp(1:nypts)
     else if (open_boundary.eq.2) then
        ny = nypts
        allocate(y(ny+1))
        y(2:ny+1) = ytmp(1:nypts)
        y(1) = 2.0_WP*y(2) - y(3)
     else if (open_boundary.eq.3) then
        ny = nypts
        allocate(y(ny+1))
        y(1:ny) = ytmp(1:nypts)
        y(ny+1) = 2.0_WP*y(ny) - y(ny-1)
     else
        ny = nypts+1
        allocate(y(ny+1))
        y(2:ny) = ytmp(1:nypts)
        y(1)    = 2.0_WP*y(2)  - y(3)
        y(ny+1) = 2.0_WP*y(ny) - y(ny-1)
     end if
  end if
  
  write(*,'(a,i4,a,i4)') " Number of points detected :",nxpts," x ",nypts

  ! Deallocate
  deallocate(xtmp)
  deallocate(ytmp)
  deallocate(xnbr)
  deallocate(ynbr)

  return
end subroutine read_points


! ======================== !
! Detect the masked points !
! ======================== !
subroutine mesh_mask
  use gambit_reader
  implicit none
  
  real(WP) :: dist
  integer  :: i,j,k,kk,np
  
  ! Create the mapping structured -> unstructured
  allocate(got_points(nxpts,nypts))
  got_points = 0
  i = -1
  j = -1
  do np=1,numnp
     ! Find x index
     dist = huge(1.0_WP)
     do k=1,nxpts
        if (abs(x(k)-xval(np)).le.dist) then
           dist = abs(x(k)-xval(np))
           i = k
        end if
     end do
     ! Find y index
     dist = huge(1.0_WP)
     do k=1,nypts
        if (icyl.eq.0 .and. open_boundary.eq.0) then
           ! Shift by one cell
           kk = k+1
        else if (icyl.eq.0 .and. open_boundary.eq.2) then
           ! Shift by one cell
           kk = k+1
        else
           kk = k
        end if
        if (abs(y(kk)-yval(np)).le.dist) then
           dist = abs(y(kk)-yval(np))
           j = k
        end if
     end do
     if (i.eq.-1 .or. j.eq.-1) stop "mesh_mask: unknown point"
     got_points(i,j) = 1
  end do
  
  ! Masked boundaries (0=interior; 1=wall)
  allocate (mask(nx,ny))
  mask = 1
  
  ! Check if cell has 4 points
  do j=1,nypts-1
     do i=1,nxpts-1
        ! Has 4 points => not wall
        if (got_points(i,j)  .eq.1 .and. got_points(i,j+1)  .eq.1 .and. &
             got_points(i+1,j).eq.1 .and. got_points(i+1,j+1).eq.1 ) then
           ! Depends on config
           if (icyl.eq.0 .and. open_boundary.eq.0) then
              ! Shift by one cell
              mask(i,j+1) = 0
           else if (icyl.eq.0 .and. open_boundary.eq.2) then
              ! Shift by one cell
              mask(i,j+1) = 0
           else
              mask(i,j) = 0
           end if
        end if
     end do
  end do
  
  ! Add walls if not open_boundary 
  if (open_boundary.eq.0) then
     ! Walls at top and bottom
     mask(:,ny) = 1
     if (icyl.eq.0) then
        mask(:,1) = 1
     end if
  else if (open_boundary.eq.2) then
     ! Walls at bottom
     if (icyl.eq.0) then
        mask(:,1) = 1
     end if
  else if (open_boundary.eq.3) then
     ! Walls at top
     mask(:,ny) = 1
  end if
  
  return
end subroutine mesh_mask


