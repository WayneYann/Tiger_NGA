module parallel
  use precision
  use string
  use mpi
  !$ use omp_lib
  implicit none
  
  ! FORTRAN notation (1->nproc)
  integer :: nproc,npx,npy,npz
  integer :: irank,iroot
  integer :: iproc,jproc,kproc
  integer, dimension(3) :: periodicity
  integer :: comm,comm_x,comm_y,comm_z
  integer :: comm_xy,comm_yz,comm_xz
  integer :: irank_x,irank_y,irank_z
  integer :: MPI_REAL_WP,MPI_REAL_SP
  !include 'mpif.h'

  ! OpenMP Threads
  !$ integer :: nthreads
  
  ! MPI Info for panfs
  character(len=str_short) :: mpiiofs
  integer :: mpi_info
  
  ! ====================================== !
  ! Compute the global maximum of anything !
  ! ====================================== !
  interface parallel_max
     module procedure parallel_max_real_3d
     module procedure parallel_max_real_0d
  end interface parallel_max
  ! =========================================== !
  ! Compute the directional maximum of anything !
  ! =========================================== !
  interface parallel_max_dir
     module procedure parallel_max_real_1d
     module procedure parallel_max_real_2d
  end interface
  ! ====================================== !
  ! Compute the global minimum of anything !
  ! ====================================== !
  interface parallel_min
     module procedure parallel_min_real_3d
     module procedure parallel_min_real_0d
  end interface parallel_min
  ! ================================== !
  ! Compute the global sum of anything !
  ! ================================== !
  interface parallel_sum
     module procedure parallel_sum_int_0d
     module procedure parallel_sum_int_2d
     module procedure parallel_sum_real_0d
     module procedure parallel_sum_real_1d
     module procedure parallel_sum_real_2d
     module procedure parallel_sum_real_3d
  end interface parallel_sum
  ! ================= !
  ! Perform broadcast !
  ! ================= !
  interface parallel_bc
     module procedure parallel_bc_char
  end interface parallel_bc
  ! ============================= !
  ! Perform directional summation !
  ! ============================= !
  interface parallel_sum_dir
     module procedure parallel_sum_dir_real_0d
     module procedure parallel_sum_dir_real_1d
     module procedure parallel_sum_dir_real_2d
     module procedure parallel_sum_dir_real_3d
  end interface parallel_sum_dir
  ! ============================= !
  ! Perform directional gathering !
  ! ============================= !
  interface parallel_gather_dir
     module procedure parallel_gather_dir_real_1d
  end interface parallel_gather_dir
  
contains

  ! MPI ALLGATHER
  subroutine parallel_gather_dir_real_1d(A,B,dir)
    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: ierr,n

    n = size(A)

    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLGATHER(A,n,MPI_REAL_WP,B,n,MPI_REAL_WP,comm_x,ierr)
    case ('y')
       call MPI_ALLGATHER(A,n,MPI_REAL_WP,B,n,MPI_REAL_WP,comm_y,ierr)
    case ('z')
       call MPI_ALLGATHER(A,n,MPI_REAL_WP,B,n,MPI_REAL_WP,comm_z,ierr)
    end select
    
    return
  end subroutine parallel_gather_dir_real_1d

  ! MPI SUM-DIR
  subroutine parallel_sum_dir_real_0d(A,B,dir)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: ierr,n
    
    n = 1
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_sum_dir_real_0d
  
  ! MPI SUM-DIR
  subroutine parallel_sum_dir_real_1d(A,B,dir)
    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: ierr,n
    
    n = size(A)
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_sum_dir_real_1d
  
  subroutine parallel_sum_dir_real_2d(A,B,dir)
    implicit none
    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: ierr,n
    
    n = size(A)
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_sum_dir_real_2d
  
  subroutine parallel_sum_dir_real_3d(A,B,dir)
    implicit none
    real(WP), dimension(:,:,:), intent(in)  :: A
    real(WP), dimension(:,:,:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: ierr,n
    
    n = size(A)
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_SUM,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_sum_dir_real_3d
  
  ! MPI MAX
  subroutine parallel_max_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierr
    C = maxval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
    return
  end subroutine parallel_max_real_3d
  
  subroutine parallel_max_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
    return
  end subroutine parallel_max_real_0d
  
  ! MPI DIRECTIONAL MAX
  subroutine parallel_max_real_1d(A,B,dir)
    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: n,ierr
    
    n = size(A)
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_max_real_1d
  
  subroutine parallel_max_real_2d(A,B,dir)
    implicit none
    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    character(len=*), intent(in) :: dir
    integer :: n,ierr
    
    n = size(A)
    
    select case(adjustl(trim(dir)))
    case ('x')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_x,ierr)
    case ('y')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_y,ierr)
    case ('z')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_z,ierr)
    case ('xy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xy,ierr)
    case ('yx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xy,ierr)
    case ('xz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xz,ierr)
    case ('zx')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_xz,ierr)
    case ('yz')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_yz,ierr)
    case ('zy')
       call MPI_ALLREDUCE(A,B,n,MPI_REAL_WP,MPI_MAX,comm_yz,ierr)
    end select
    
    return
  end subroutine parallel_max_real_2d
  
  ! MPI MIN
  subroutine parallel_min_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierr
    C = minval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
    return
  end subroutine parallel_min_real_3d
  
  subroutine parallel_min_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
    return
  end subroutine parallel_min_real_0d
  
  ! MPI SUM
  subroutine parallel_sum_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_real_0d
  
  subroutine parallel_sum_int_0d(A,B)
    implicit none
    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_int_0d

  subroutine parallel_sum_int_2d(A,B)
    implicit none
    integer, dimension(:,:), intent(in)  :: A
    integer, dimension(:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_int_2d
  
  subroutine parallel_sum_real_1d(A,B)
    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_real_1d
  
  subroutine parallel_sum_real_2d(A,B)
    implicit none
    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_real_2d

  subroutine parallel_sum_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in)  :: A
    real(WP), dimension(:,:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
    return
  end subroutine parallel_sum_real_3d
  
  ! MPI KILL
  subroutine parallel_kill(error_text)
    implicit none
    integer :: ierr
    character(len=*), intent(in), optional :: error_text
    
    ! Specify who sends the abort signal and what it means
    if (present(error_text)) then
       write(*,'(a,i3,a,/,a)') '[',irank,'] initiated general abort signal due to the following error:',trim(error_text)
    else
       write(*,'(a,i3,a)') '[',irank,'] initiated general abort signal due to an unknown error'
    end if
    
    ! Call general abort
    call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
    
    return
  end subroutine parallel_kill
  
  subroutine parallel_bc_char(text)
    implicit none
    integer :: ierr
    character(len=*) :: text
    
    call MPI_BCAST(text,len(text),MPI_CHARACTER,iroot-1,comm,ierr)
    
    return
  end subroutine parallel_bc_char
  
  subroutine parallel_bcast(r0,r1,r2,r3,i0,i1,i2,i3)
    implicit none
    integer :: ierr
    real(WP),                   optional, intent(inout) :: r0
    real(WP), dimension(:),     optional, intent(inout) :: r1
    real(WP), dimension(:,:),   optional, intent(inout) :: r2
    real(WP), dimension(:,:,:), optional, intent(inout) :: r3
    integer,                    optional, intent(inout) :: i0
    integer,  dimension(:),     optional, intent(inout) :: i1
    integer,  dimension(:,:),   optional, intent(inout) :: i2
    integer,  dimension(:,:,:), optional, intent(inout) :: i3

    if (present(r0)) call MPI_BCAST(r0,1,MPI_REAL_WP,iroot-1,comm,ierr)
    if (present(r1)) call MPI_BCAST(r1,size(r1),MPI_REAL_WP,iroot-1,comm,ierr)
    if (present(r2)) call MPI_BCAST(r2,size(r2),MPI_REAL_WP,iroot-1,comm,ierr)
    if (present(r3)) call MPI_BCAST(r3,size(r3),MPI_REAL_WP,iroot-1,comm,ierr)
    if (present(i0)) call MPI_BCAST(i0,1,MPI_INTEGER,iroot-1,comm,ierr)
    if (present(i1)) call MPI_BCAST(i1,size(i1),MPI_INTEGER,iroot-1,comm,ierr)
    if (present(i2)) call MPI_BCAST(i2,size(i2),MPI_INTEGER,iroot-1,comm,ierr)
    if (present(i3)) call MPI_BCAST(i3,size(i3),MPI_INTEGER,iroot-1,comm,ierr)
    
    return
  end subroutine parallel_bcast
  
end module parallel

  
subroutine parallel_init
  use parallel
  use parser
  implicit none
  integer :: ierr
  integer :: size_real,size_dp
  !$ integer :: thread_supp
  !$ logical :: use_threads
  
  ! Initialize a first basic MPI environment
  !$ use_threads = .true.
  !$ if (use_threads) then
  !$   call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,thread_supp,ierr)
  !!$   if (thread_supp.lt.MPI_THREAD_FUNNELED) then
  !!$      call parallel_kill('Error in parallel_init: no threads supported in MPI')
  !!$   end if
  !$ else
       call MPI_INIT(ierr)
  !$ end if
  call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr) 
  irank = irank+1
  iroot = 1
  
  ! Set MPI working precision - WP
  call MPI_TYPE_SIZE(MPI_REAL,size_real,ierr)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_dp,ierr)
  if (WP .eq. size_real) then
     MPI_REAL_WP = MPI_REAL
  else if (WP .eq. size_dp) then
     MPI_REAL_WP = MPI_DOUBLE_PRECISION
  else
     call parallel_kill('Error in parallel_init: no WP equivalent in MPI')
  end if
  
  ! Set MPI single precision
  call MPI_TYPE_SIZE(MPI_REAL,size_real,ierr)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_dp,ierr)
  if (SP .eq. size_real) then
     MPI_REAL_SP = MPI_REAL
  else if (SP .eq. size_dp) then
     MPI_REAL_SP = MPI_DOUBLE_PRECISION
  else
     call parallel_kill('Error in parallel_init: no SP equivalent in MPI')
  end if
  
  return
end subroutine parallel_init


subroutine parallel_init_topology(xper,yper,zper)
  use parallel
  use parser
  implicit none
  integer :: ierr
  integer, dimension(3) :: dims
  logical, dimension(3) :: isper
  logical, dimension(3) :: dir
  logical :: reorder
  integer :: ndims
  integer, dimension(3) :: coords
  integer, intent(in) :: xper,yper,zper
  
  ! Save periodicity
  periodicity(1) = xper
  periodicity(2) = yper
  periodicity(3) = zper
  
  ! Read topology from input file
  call parser_read('Processors along X',npx)
  call parser_read('Processors along Y',npy)
  call parser_read('Processors along Z',npz)
  !$ call parser_read('OpenMP threads',nthreads)
  !$ call OMP_SET_NUM_THREADS(nthreads)
  
  ! Test if nproc is correct
  if (nproc .ne. npx*npy*npz) call parallel_kill('Wrong number of cpus specified in input file')
  
  ! Set MPI topology
  ndims = 3
  dims(1) = npx
  dims(2) = npy
  dims(3) = npz
  if (xper.EQ.1) then
     isper(1) = .true.
  else
     isper(1) = .false.
  end if
  if (yper.EQ.1) then
     isper(2) = .true.
  else
     isper(2) = .false.
  end if
  if (zper.EQ.1) then
     isper(3) = .true.
  else
     isper(3) = .false.
  end if
  reorder = .true.
  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,isper,reorder,comm,ierr)
  call MPI_COMM_RANK(comm,irank,ierr)
  call MPI_CART_COORDS(comm,irank,ndims,coords,ierr)
  irank = irank + 1
  iproc = coords(1) + 1
  jproc = coords(2) + 1
  kproc = coords(3) + 1
  
  ! Define a root processor at coordinates (1,1,1)
  dims = 0
  call MPI_CART_RANK(comm,dims,iroot,ierr)
  iroot = iroot + 1
  
  ! Create line communicators
  ! Along x
  dir(1) = .true.
  dir(2) = .false.
  dir(3) = .false.
  call MPI_CART_SUB(comm,dir,comm_x,ierr)
  call MPI_COMM_RANK(comm_x,irank_x,ierr)
  ! Along y
  dir(1) = .false.
  dir(2) = .true.
  dir(3) = .false.
  call MPI_CART_SUB(comm,dir,comm_y,ierr)
  call MPI_COMM_RANK(comm_y,irank_y,ierr)
  ! Along z
  dir(1) = .false.
  dir(2) = .false.
  dir(3) = .true.
  call MPI_CART_SUB(comm,dir,comm_z,ierr)
  call MPI_COMM_RANK(comm_z,irank_z,ierr)
  
  ! Create planar communicators
  ! Along xy
  dir(1) = .true.
  dir(2) = .true.
  dir(3) = .false.
  call MPI_CART_SUB(comm,dir,comm_xy,ierr)
  ! Along yz
  dir(1) = .false.
  dir(2) = .true.
  dir(3) = .true.
  call MPI_CART_SUB(comm,dir,comm_yz,ierr)
  ! Along xz
  dir(1) = .true.
  dir(2) = .false.
  dir(3) = .true.
  call MPI_CART_SUB(comm,dir,comm_xz,ierr)
  
  ! Test
  !call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'-x'
  !call MPI_CART_SHIFT(comm,0,+1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'+x'
  !call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'-y'
  !call MPI_CART_SHIFT(comm,1,+1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'+y'
  !call MPI_CART_SHIFT(comm,2,-1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'-z'
  !call MPI_CART_SHIFT(comm,2,+1,isource,idest,ierr)
  !print*,'[',irank,']',isource,idest,'+z'
  !stop
   
  
  return
end subroutine parallel_init_topology


subroutine parallel_init_filesystem
  use parallel
  use parser
  implicit none
  integer :: ierr

  ! Create info
  call parser_read('MPIIO fs',mpiiofs,'ufs')
  select case(trim(mpiiofs))
  case('panfs')
     call MPI_INFO_CREATE(mpi_info,ierr)
     call MPI_INFO_SET(mpi_info,"panfs_concurrent_write","1",ierr)
  case('ufs')
     mpi_info = MPI_INFO_NULL
  case('lustre')
     call MPI_INFO_CREATE(mpi_info,ierr)
     call MPI_INFO_SET(mpi_info,"romio_ds_write","disable",ierr)     
  case default
     call parallel_kill('parallel_init_topology: unknown mpiio fs')
  end select  

  return
end subroutine parallel_init_filesystem

subroutine parallel_final
  use parallel
  implicit none
  integer :: ierr
  
  ! Finalize MPI
  call MPI_FINALIZE(ierr)
  
  return
end subroutine parallel_final


subroutine parallel_time(wtime)
  use parallel
  implicit none
  real(WP), intent(out) :: wtime
  
  ! Get time
  wtime = MPI_WTIME()
  
  return
end subroutine parallel_time


subroutine parallel_get_inputname(input_name)
  use parallel
  use cli_reader
  implicit none
  character(len=str_medium), intent(inout) :: input_name
  integer :: ierr
  
  if (irank.eq.iroot) then
     if (command_argument_count().ge.1) then
        call get_command_argument(1,input_name)
        if (input_name(1:1).eq.'-' .or. len_trim(input_name).eq.0) then
           print*,'No input file name was detected, using "input".'
           input_name='input'
        end if
     else
        print*,'No input file name was detected, using "input".'
        input_name='input'
     end if
  end if
  
  call MPI_BCAST(input_name,str_medium,MPI_CHARACTER,iroot-1,MPI_COMM_WORLD,ierr)
  
  return
end subroutine parallel_get_inputname
