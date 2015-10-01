module quicksort
  use precision
  implicit none
  
  public :: quick_sort
  private :: qs_partition
  
contains
  
  ! ------------------------------------------------ !
  ! Recursive QuickSort algorithm : this routine     !
  ! recursively sorts a real(WP) array by increasing !
  ! order and sorts a integer array at the same time !
  ! ------------------------------------------------ !
  recursive subroutine quick_sort(A,B)
    use precision
    implicit none
    
    real(WP), dimension(:) :: A
    integer, dimension(:)  :: B
    integer :: imark
    
    if (size(A).gt.1) then
       call qs_partition(A,B,imark)
       call quick_sort(A(:imark-1),B(:imark-1))
       call quick_sort(A(imark:),B(imark:))
    end if
    
    return
  end subroutine quick_sort
  
  ! -------------------------------------- !
  ! Private routine belonging to quicksort !
  ! -------------------------------------- !
  subroutine qs_partition(A,B,marker)
    use precision
    implicit none
    
    real(WP), dimension(:) :: A
    integer,  dimension(:) :: B
    integer, intent(out) :: marker
    integer :: i,j
    real(WP) :: dtmp
    integer  :: itmp
    real(WP) :: x
    
    x = A(1)
    i = 0
    j = size(A) + 1
    
    do
       j = j-1
       do
          if (A(j).le.x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i).ge.x) exit
          i = i+1
       end do
       if (i.lt.j) then
          ! Exchange A(i) and A(j)
          dtmp = A(i)
          A(i) = A(j)
          A(j) = dtmp
          ! Also exchange B(i) and B(j)
          itmp = B(i)
          B(i) = B(j)
          B(j) = itmp
       else if (i.eq.j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
    
    return
  end subroutine qs_partition
  
end module quicksort
