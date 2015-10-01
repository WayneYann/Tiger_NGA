module stat
  use precision
  use geometry
  use partition
  use string
  implicit none

  ! Number and type of statistics
  integer :: nstat
  character(str_medium), dimension(:), allocatable :: stat_type
  logical :: cond_stat

  ! Frequency of writing
  integer :: nfreq_writing
  real(WP), dimension(:), allocatable :: writing_freq
  integer,  dimension(:), allocatable :: writing_iter
  
end module stat

! ======================================== !
! Initialize statistics dumping procedures !
! ======================================== !
subroutine stat_init
  use stat
  use parser
  use time_info
  implicit none

  logical :: isdef
  integer :: i
  real(WP) :: buffer
  
  ! Read the input file
  call parser_is_defined('Statistics type',isdef)
  nstat = 0

  ! If nothing to do returns
  if (.not. isdef) return
  
  ! Create the directory
  call system("mkdir -p stat")
  
  call parser_getsize('Statistics type',nstat)
  allocate(stat_type(nstat))
  call parser_read('Statistics type',stat_type)
  call parser_read('Conditional statistics',cond_stat,.false.)
  
  call parser_getsize('Statistics frequency',nfreq_writing)
  allocate(writing_freq(nfreq_writing))
  call parser_read('Statistics frequency',writing_freq)
  if (nfreq_writing.EQ.1) then
     buffer = writing_freq(1)
     deallocate(writing_freq)
     allocate(writing_freq(nstat))
     writing_freq = buffer
  else if (nfreq_writing.NE.nstat) then
     call die('stat_init: wrong number of statistics writing frequencies specified')
  end if
  allocate(writing_iter(nstat))
  writing_iter(:) = int(time/writing_freq(:))
  
  ! Initialize
  do i=1,nstat
     select case (trim(stat_type(i)))
     case ('1D-x')
        call stat_1Dx_init
     case ('1D-y')
        call stat_1Dy_init
     case ('2D')
        call stat_2D_init
     case ('3D')
        call stat_3D_init
     case default
        call die('stat_init: unknown output type ('//trim(stat_type(i))//')')
     end select
  end do
  
  return
end subroutine stat_init


! ======================================== !
! Dump_stat at each time step if necessary !
! ======================================== !
subroutine dump_statistics
  use stat
  use time_info
  implicit none
  integer :: i
  
  ! Compute stat
  do i=1,nstat
     select case(trim(stat_type(i)))
     case('1D-x')
        call stat_1Dx_sample
     case('1D-y')
        call stat_1Dy_sample
    case('2D')
        call stat_2D_sample
     case('3D')
        call stat_3D_sample
     end select
  end do
  
  ! Dump stat if needed
  do i=1,nstat
     if (int(time/writing_freq(i)).NE.writing_iter(i)) then
        writing_iter(i) = int(time/writing_freq(i))
        select case(trim(stat_type(i)))
        case('1D-x')
           call stat_1Dx_write
        case('1D-y')
           call stat_1Dy_write
       case('2D')
           call stat_2D_write
           call monitor_log("2D STATISTICS FILE WRITTEN")
        case('3D')
           call stat_3D_write
        end select
     end if
  end do
  
  return
end subroutine dump_statistics


! ======================================== !
! Dump_stat at each time step if necessary !
! ======================================== !
subroutine stat_save
  use stat
  implicit none
  integer :: i
  
  do i=1,nstat
     select case(trim(stat_type(i)))
     case('1D-x')
        call stat_1Dx_write
     case('1D-y')
        call stat_1Dy_write
     case('2D')
        call stat_2D_write
        call monitor_log("2D STATISTICS FILE WRITTEN")
     case('3D')
        call stat_3D_write
     end select
  end do
  
  return
end subroutine stat_save



