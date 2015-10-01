! ============================================ !
! Update the borders:                          !
!   - periodic boundaries                      !
!   - inside boundaries (domain decomposition) !
! ============================================ ! 

!========================================!
! BASIC version, suitable for most cases !
!========================================!
subroutine communication_border(A)
  use parallel
  use geometry
  use partition
  implicit none
  real(WP), dimension(nxo_,nyo_,nzo_) :: A
  
  call communication_border_x(A,nxo_,nyo_,nzo_,nover)
  call communication_border_y(A,nxo_,nyo_,nzo_,nover)
  call communication_border_z(A,nxo_,nyo_,nzo_,nover)
  
  return
end subroutine communication_border

!========================================!
! BASIC version, suitable for most cases !
! FOR INTEGERS                           !
!========================================!
subroutine communication_border_int(A)
  use parallel
  use geometry
  use partition
  implicit none
  integer, dimension(nxo_,nyo_,nzo_) :: A
  
  call communication_border_x_int(A,nxo_,nyo_,nzo_,nover)
  call communication_border_y_int(A,nxo_,nyo_,nzo_,nover)
  call communication_border_z_int(A,nxo_,nyo_,nzo_,nover)
  
  return
end subroutine communication_border_int

!=====================================!
! ADVANCED version, for SGS typically !
!=====================================!
subroutine communication_border_adv(A,n1,n2,n3,no)
  use parallel
  use partition
  implicit none
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  
  call communication_border_x(A,n1,n2,n3,no)
  call communication_border_y(A,n1,n2,n3,no)
  call communication_border_z(A,n1,n2,n3,no)
  
  return
end subroutine communication_border_adv

!=============================================================================!

subroutine communication_border_x(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: i
  
  ! If nx=1 then it is 2D
  if (nx.eq.1) then
     do i=1,no
        A(i,:,:) = A(1+no,:,:)
        A(n1-i+1,:,:) = A(1+no,:,:)
     end do
     return
  end if
  
  ! Initialize buffer
  allocate(buf1(no,n2,n3))
  allocate(buf2(no,n2,n3))
  icount = no*n2*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(1+no:1+no+no-1,:,:)
  
  ! Send left buffer to left neighbour
  call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(n1-no+1:n1,:,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(n1-no-no+1:n1-no,:,:)
  
  ! Send right buffer to right neighbour
  call MPI_CART_SHIFT(comm,0,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(1:no,:,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_x

!=============================================================================!

subroutine communication_border_y(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: j
  
  ! If ny=1 then it is 2D
  if (ny.eq.1) then
     do j=1,no
        A(:,j,:) = A(:,1+no,:)
        A(:,n2-j+1,:) = A(:,1+no,:)
     end do
     return
  end if

  ! Initialize buffer
  allocate(buf1(n1,no,n3))
  allocate(buf2(n1,no,n3))
  icount = n1*no*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,1+no:1+no+no-1,:)
  
  ! Send lower buffer to lower neighbour
  call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste lower buffer to top
  if (isource.NE.MPI_PROC_NULL) A(:,n2-no+1:n2,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,n2-no-no+1:n2-no,:)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,1,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,1:no,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_y

!=============================================================================!

subroutine communication_border_z(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: k
  
  ! If nz=1 then it is 2D
  if (nz.eq.1) then
     do k=1,no
        A(:,:,k) = A(:,:,1+no)
        A(:,:,n3-k+1) = A(:,:,1+no)
     end do
     return
  end if

  ! Initialize buffer
  allocate(buf1(n1,n2,no))
  allocate(buf2(n1,n2,no))
  icount = n1*n2*no
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,:,1+no:1+no+no-1)
  
  ! Send left buffer to lower neighbour
  call MPI_CART_SHIFT(comm,2,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(:,:,n3-no+1:n3) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,:,n3-no-no+1:n3-no)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,2,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,:,1:no) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_z

!=============================================================================!

subroutine communication_border_x_int(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  integer, dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  integer, dimension(:,:,:), pointer :: buf1
  integer, dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: i

  ! If ny=1 then it is 2D
  if (nx.eq.1) then
     do i=1,no
        A(i,:,:) = A(1+no,:,:)
        A(n1-i+1,:,:) = A(1+no,:,:)
     end do
     return
  end if
  
  ! Initialize buffer
  allocate(buf1(no,n2,n3))
  allocate(buf2(no,n2,n3))
  icount = no*n2*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(1+no:1+no+no-1,:,:)
  
  ! Send left buffer to left neighbour
  call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(n1-no+1:n1,:,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(n1-no-no+1:n1-no,:,:)
  
  ! Send right buffer to right neighbour
  call MPI_CART_SHIFT(comm,0,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(1:no,:,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_x_int

!=============================================================================!

subroutine communication_border_y_int(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  integer, dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  integer, dimension(:,:,:), pointer :: buf1
  integer, dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: j
  
  ! If ny=1 then it is 2D
  if (ny.eq.1) then
     do j=1,no
        A(:,j,:) = A(:,1+no,:)
        A(:,n2-j+1,:) = A(:,1+no,:)
     end do
     return
  end if
  
  ! Initialize buffer
  allocate(buf1(n1,no,n3))
  allocate(buf2(n1,no,n3))
  icount = n1*no*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,1+no:1+no+no-1,:)
  
  ! Send lower buffer to lower neighbour
  call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste lower buffer to top
  if (isource.NE.MPI_PROC_NULL) A(:,n2-no+1:n2,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,n2-no-no+1:n2-no,:)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,1,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,1:no,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_y_int

!=============================================================================!

subroutine communication_border_z_int(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  integer, dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  integer, dimension(:,:,:), pointer :: buf1
  integer, dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: k
  
  ! If nz=1 then it is 2D
  if (nz.eq.1) then
     do k=1,no
        A(:,:,k) = A(:,:,1+no)
        A(:,:,n3-k+1) = A(:,:,1+no)
     end do
     return
  end if

  ! Initialize buffer
  allocate(buf1(n1,n2,no))
  allocate(buf2(n1,n2,no))
  icount = n1*n2*no
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,:,1+no:1+no+no-1)
  
  ! Send left buffer to lower neighbour
  call MPI_CART_SHIFT(comm,2,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(:,:,n3-no+1:n3) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,:,n3-no-no+1:n3-no)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,2,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_INTEGER,idest,0,buf2,icount,MPI_INTEGER,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,:,1:no) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine communication_border_z_int

!========================================!
! BASIC version, suitable for most cases !
!========================================!
subroutine summation_border(A)
  use parallel
  use geometry
  use partition
  implicit none
  real(WP), dimension(nxo_,nyo_,nzo_) :: A
  
  call summation_border_x(A,nxo_,nyo_,nzo_,nover)
  call summation_border_y(A,nxo_,nyo_,nzo_,nover)
  call summation_border_z(A,nxo_,nyo_,nzo_,nover)
  
  return
end subroutine summation_border

!=============================================================================!

subroutine summation_border_x(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: i

  ! If nx=1 then it is 2D
  if (nx.eq.1) then
     do i=1,no
        A(1+no,:,:) = A(1+no,:,:) + A(i,:,:) + A(n1-i+1,:,:)
     end do
     return
  end if
  
  ! Initialize buffer
  allocate(buf1(no,n2,n3))
  allocate(buf2(no,n2,n3))
  icount = no*n2*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(1:no,:,:)
  !buf1(:,:,:) = A(1+no:1+no+no-1,:,:)
  
  ! Send left buffer to left neighbour
  call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(n1-no-no+1:n1-no,:,:) = A(n1-no-no+1:n1-no,:,:) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(n1-no+1:n1,:,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(n1-no+1:n1,:,:)
  !buf1(:,:,:) = A(n1-no-no+1:n1-no,:,:)
  
  ! Send right buffer to right neighbour
  call MPI_CART_SHIFT(comm,0,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(1+no:1+no+no-1,:,:) = A(1+no:1+no+no-1,:,:) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(1:no,:,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine summation_border_x

!=============================================================================!

subroutine summation_border_y(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: j

  ! If ny=1 then it is 2D
  if (ny.eq.1) then
     do j=1,no
        A(:,1+no,:) = A(:,1+no,:) + A(:,j,:) + A(:,n2-j+1,:)
     end do
     return
  end if
  
  ! Initialize buffer
  allocate(buf1(n1,no,n3))
  allocate(buf2(n1,no,n3))
  icount = n1*no*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,1:no,:)
  !buf1(:,:,:) = A(:,1+no:1+no+no-1,:)
  
  ! Send lower buffer to lower neighbour
  call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste lower buffer to top
  if (isource.NE.MPI_PROC_NULL) A(:,n2-no-no+1:n2-no,:) = A(:,n2-no-no+1:n2-no,:) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(:,n2-no+1:n2,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,n2-no+1:n2,:)
  !buf1(:,:,:) = A(:,n2-no-no+1:n2-no,:)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,1,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,1+no:1+no+no-1,:) = A(:,1+no:1+no+no-1,:) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(:,1:no,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine summation_border_y

!=============================================================================!

subroutine summation_border_z(A,n1,n2,n3,no)
  use parallel
  use geometry
  use partition
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  integer :: k
  
  ! If nz=1 then it is 2D
  if (nz.eq.1) then
     do k=1,no
        A(:,:,1+no) = A(:,:,1+no) + A(:,:,n3-k+1) + A(:,:,k)
     end do
     return
  end if

  ! Initialize buffer
  allocate(buf1(n1,n2,no))
  allocate(buf2(n1,n2,no))
  icount = n1*n2*no
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,:,1:no)
  !buf1(:,:,:) = A(:,:,1+no:1+no+no-1)
  
  ! Send left buffer to lower neighbour
  call MPI_CART_SHIFT(comm,2,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(:,:,n3-no-no+1:n3-no) = A(:,:,n3-no-no+1:n3-no) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(:,:,n3-no+1:n3) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,:,n3-no+1:n3)
  !buf1(:,:,:) = A(:,:,n3-no-no+1:n3-no)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,2,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,:,1+no:1+no+no-1) = A(:,:,1+no:1+no+no-1) + buf2(:,:,:)
  !if (isource.NE.MPI_PROC_NULL) A(:,:,1:no) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine summation_border_z
