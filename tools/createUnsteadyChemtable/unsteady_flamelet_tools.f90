module unsteady_flamelet_tools
  use precision
  implicit none

  integer :: n_triangles
  integer, dimension(:,:), pointer :: tri_verts  !vertices of the triangles
  real(WP), dimension(:), pointer :: Z3_nrmlzd, Z4_nrmlzd  !normalized C & H 

  integer, dimension(:), allocatable :: file_mapping

end module unsteady_flamelet_tools


!==================================================================!
! Function to determine whether a location is inside of a triangle !
!==================================================================!
function chk_pt_in_triangle(a,b,c,p,lambda,strict)
  use precision 
  implicit none

  logical :: chk_pt_in_triangle    ! the point is in the triangle (or not)
  integer, intent(in) :: strict
  real(WP), dimension(2), intent(in) :: a,b,c,p
  real(WP), dimension(3), intent(out) :: lambda
  real(WP) :: inv_deter, deter
  real(WP) :: u_bound, l_bound
  integer :: i

  ! Idea: Transform to barycentric coordinates lambda, so that any point p 
  ! can be written as p = lambda_1 * r_1 + lambda_2 * r_2 + lambda_3 * r_3,
  ! where r_1, r_2, and r_3 are the (x,y) values describing the three 
  ! triangular vertices.  Then, a point is inside the triangle if all the 
  ! lambda's are between zero and one.  

  if (strict.eq.1) then
     u_bound = 1.0_WP+1.0E-8_WP; l_bound = 0.0_WP-1.0E-8_WP
  else
     u_bound = 1.0_WP+1.0E-3_WP; l_bound = 0.0_WP-1.0E-3_WP
  end if

  ! Compute 1/determinant associated with the coordinate transformation
  deter = (a(1)-c(1))*(b(2)-c(2))-(b(1)-c(1))*(a(2)-c(2))
  !write(*,*) 'determinant = ',deter
  !write(*,*) 'determinant 1 = ',(a(1)-c(1))*(b(2)-c(2))
  !write(*,*) 'determinant 2 = ',(b(1)-c(1))*(a(2)-c(2))
  if (abs(deter).le.1.0E-10_WP) then
  !if (deter.eq.0.0_WP) then
     chk_pt_in_triangle = .false.
     lambda(:) = 0.0_WP
     return
  end if
  inv_deter = 1.0_WP / deter 

  lambda(1) = ( (b(2)-c(2))*(p(1)-c(1)) - (b(1)-c(1))*(p(2)-c(2)))*inv_deter
  lambda(2) = (-(a(2)-c(2))*(p(1)-c(1)) + (a(1)-c(1))*(p(2)-c(2)))*inv_deter
  lambda(3) = 1.0_WP - lambda(1) - lambda(2)
  
  chk_pt_in_triangle = .true.
  chkptloop: do i=1,3
     if (lambda(i).gt.u_bound) then
        chk_pt_in_triangle = .false.
        exit chkptloop
     end if
     if (lambda(i).lt.l_bound) then
        chk_pt_in_triangle = .false.
        exit chkptloop
     end if
  end do chkptloop

  return
end function


function finddet(matrix,n)
  use precision
  implicit none

  real(WP) :: finddet
  real(WP), dimension(n,n) :: matrix
  integer, intent(in) :: n
  real(WP) :: m, temp
  integer :: i,j,k,l
  logical :: detexists
  
  l = 1
 
  detexists = .true.
  ! Convert to upper triangular form
  do k=1,n-1
     if (matrix(k,k).eq.0) then
        detexists = .false. 
        do i=k+1,n
           if (matrix(i,k).ne.0) then
              do j=1,n
                 temp = matrix(i,j)
                 matrix(i,j) = matrix(k,j)
                 matrix(k,j) = temp
              end do
              detexists = .true.
              l=-l
              exit
           end if
        end do
        if (detexists.eq..false.) then
           finddet = 0.0_WP
           return
        end if
     end if
     do j=k+1,n
        m=matrix(j,k)/matrix(k,k)
        do i=k+1,n
           matrix(j,i) = matrix(j,i)-m*matrix(k,i)
        end do
     end do
  end do

  ! Calculate determinant by finding product of diagonal elements
  finddet = real(l,WP)
  do i=1,n
     finddet = finddet * matrix(i,i)
  end do

  return 
end function finddet

! Inverse the linear system : AX=B
subroutine solve_linear_system(A,B,X,n)
  use precision
  implicit none

  integer, intent(in) :: n
  real(WP), dimension(n,n) :: A
  real(WP), dimension(n) :: B
  real(WP), dimension(n) :: X
  integer :: i,l

  ! Forward elimination
  do i=1,n
     B(i) = B(i) / A(i,i)
     A(i,:) = A(i,:) / A(i,i)
     do l=i+1,n
        B(l) = B(l) - A(l,i)*B(i)
        A(l,:) = A(l,:) - A(l,i)*A(i,:)
     end do
  end do

  ! Backward substitution
  do i=n,1,-1
     X(i) = B(i)
     do l=i+1,n
        X(i) = X(i) - A(i,l)*X(l)
     end do
  end do

  return
end subroutine solve_linear_system



!*******************************************************************************
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2-D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!    On abnormal return, IERROR is set to 8, 224, or 225.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer NPT, the number of vertices.
!
!    Input, integer MAXST, the maximum size available for the STACK array; 
!    should be about NPT to be safe, but MAX(10,2*LOG2(NPT)) is usually enough.
!
!    Input, double precision VCL(2,NPT), the coordinates of the vertices.
!
!    Input/output, integer IND(NPT), the indices in VCL of the vertices 
!    to be triangulated.  On output, IND has been permuted by the sort.
!
!    Output, integer NTRI, the number of triangles in the triangulation; 
!    NTRI is equal to 2*NPT - NB - 2, where NB is the number of boundary 
!    vertices.
!
!    Output, integer TIL(3,NTRI), the nodes that make up each triangle.
!    The elements are indices of VCL.  The vertices of the triangles are 
!    in counter clockwise order.
!
!    Output, integer TNBR(3,NTRI), the triangle neighbor list.
!    Positive elements are indices of TIL; negative elements are used for links
!    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(MAXST), used for a stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer IERROR, an error flag, nonzero if an error occurred.
!
subroutine dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )
  use precision
  implicit none
  
  integer :: maxst
  integer :: npt

  real(WP) :: cmax
  integer :: e
  integer :: i
  integer :: ierror
  integer :: ind(npt)
  integer :: j
  integer :: k
  integer :: l
  integer :: ledg
  integer :: lr
  integer :: lrline
  integer :: ltri
  integer :: m
  integer :: m1
  integer :: m2
  integer, parameter :: msglvl = 0
  integer :: n
  integer :: ntri
  integer :: redg
  integer :: rtri
  integer, dimension(maxst) :: stack
  integer :: t
  integer :: til(3,npt*2)
  integer :: tnbr(3,npt*2)
  integer :: top
  real(WP) :: tol
  real(WP), dimension(2,npt) :: vcl

  ierror = 0
  tol = 100.0_WP * epsilon ( tol )

  ierror = 0

  ! Sort the vertices.
  call dhpsrt ( 2, npt, 2, vcl, ind )

  ! Ensure that no two consecutive points are too close.
  m1 = ind(1)

  do i = 2, npt
     m = m1
     m1 = ind(i)

     k = 0
     do j = 1, 2
        cmax = max ( abs ( vcl(j,m) ), abs ( vcl(j,m1) ) )
        if ( abs ( vcl(j,m) - vcl(j,m1) ) > tol * cmax .and. cmax > tol ) then
           k = j
           exit
        end if
     end do

     if ( k == 0 ) then
        ierror = 224
        return
     end if
  end do

  !  Take the first two points, M1 and M2, and find a suitable non-collinear
  !  third, M.  All points between M2 and M are very close to collinear
  !  with M1 and M2.

  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do
     if ( j > npt ) then
        ierror = 225
        return
     end if

     m = ind(j)
     lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
       vcl(2,m2), 0.0_WP )

     if ( lr /= 0 ) then
        exit
     end if

     j = j + 1
  
  end do

  ntri = j - 2
  
  ! Depending on the orientation of M1, M2, and M, set up the initial
  ! triangle data.

  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri

      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3 * i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1

    end do

    tnbr(1,ntri) = -3 * ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3 * i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3 * ntri
    tnbr(2,1) = -3 * ntri - 2
    ledg = 2
    ltri = 1

  end if

  if ( msglvl == 4 ) then

    m2 = ind(1)
    write ( *, 600 ) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)

    do i = 2, j-1
      m1 = m2
      m2 = ind(i)
      write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
      write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    end do

  end if
!
!  Insert vertices one at a time from outside the convex hull, determine
!  the visible boundary edges, and apply diagonal edge swaps until
!  the Delaunay triangulation of the vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    if ( msglvl == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, 600 ) i
    end if

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0_WP )

    if ( lr > 0 ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )

    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( top > maxst ) then
        ierror = 8
        return
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
      end if

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    if ( msglvl == 4 ) then
      write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
    end if

    tnbr(ledg,ltri) = -3 * n - 1
    tnbr(2,n) = -3 * ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, maxst, ltri, ledg, vcl, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      return
    end if

  end do

  if ( msglvl == 4 ) then
    write ( *, '(i7)' ) npt + 1
  end if

  600 format (1x,i7,4f15.7)

  return
end subroutine dtris2



!*******************************************************************************
!! DHPSRT sorts points into lexicographic order using heap sort
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    points so that the points are in lexicographic increasing order.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Modified:
!
!    19 February 2001
!
!  Parameters:
!
!    Input, integer K, the dimension of the points (for instance, 2
!    for points in the plane).
!
!    Input, integer N, the number of points.
!
!    Input, integer LDA, the leading dimension of array A in the calling
!    routine; LDA should be at least K.
!
!    Input, double precision A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
subroutine dhpsrt ( k, n, lda, a, map )
  use precision
  implicit none

  integer, intent(in) :: lda
  integer, intent(in) :: n
  
  real(WP), dimension(lda,*) :: a
  integer, dimension(n) :: map
  integer :: i, k
  
  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    call i_swap ( map(1), map(i) )
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end subroutine dhpsrt



!*******************************************************************************
!! DSFTDW sifts A(*,MAP(L)) down a heap of size U.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer L, U, the lower and upper indices of part of the heap.
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, integer LDA, the leading dimension of A in the calling routine.
!
!    Input, double precision A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted.
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
subroutine dsftdw ( l, u, k, lda, a, map )
  use precision
  implicit none

  integer, intent(in) :: lda

  real(WP), dimension(lda,*) :: a
  logical :: dless
  integer :: i, j, k
  integer :: t, u, l
  integer, dimension(*) :: map
  
  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end subroutine dsftdw



!*******************************************************************************
!! DLESS determine whether P is lexicographically less than Q.
!
!  Discussion:
!
!    P and Q are K-dimensional points.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, double precision P(K), Q(K), the points to be compared.
!
!    Output, logical RLESS, is TRUE if P < Q, FALSE otherwise.
!
function dless ( k, p, q )
  use precision
  implicit none

  integer, intent(in) :: k

  logical dless
  integer :: i
  real(WP), dimension(k) :: p, q
  real(WP) :: tol
  real(WP) :: cmax

  tol = 100.0_WP * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) > tol * cmax .and. cmax > tol ) then

      if ( p(i) < q(i) ) then
        dless = .true.
      else
        dless = .false.
      end if

      return
    end if

  end do

  dless = .false.

  return
end function dless



!*******************************************************************************
!
!! I_SWAP swaps two integer values.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
subroutine i_swap ( i, j )
  use precision
  implicit none

  integer :: i,j,k
  
  k = i
  i = j
  j = k

  return
end subroutine i_swap



!*******************************************************************************
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is paralled to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, double precision XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, double precision DV, the signed distance of the directed line 
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).  
!    DV is positive for a line to the left of the base line.
! 
!    Output, integer LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )
  use precision
  implicit none

  integer :: lrline
  real(WP) :: dv, dx, dxu, dy, dyu
  real(WP) :: t, tol, tolabs
  real(WP) :: xu, xv1, xv2
  real(WP) :: yu, yv1, yv2

  tol = 100.0_WP * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end function lrline



!*******************************************************************************
!! VBEDG determines visible boundary edges of a 2D triangulation.
!
!  Purpose: 
!
!    Determine boundary edges of 2-D triangulation which are
!    visible from point (X,Y) outside convex hull.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a 2-D point outside
!    the convex hull.
!
!    Input, double precision VCL(1:2,1:*), the coordinates of 2-D vertices.
!
!    Input, integer TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer TNBR(1:3,1:*), the triangle neighbor list; negative 
!    values are used for links of counter clockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer LTRI, LEDG.  On input, if LTRI /= 0 then they 
!    are assumed to be as defined below and are not changed, else they are 
!    updated.  On output, LTRI is the index of the boundary triangle to the
!    left of leftmost boundary triangle visible from (X,Y), and LEDG is the
!    boundary edge of triangle LTRI to left of leftmost
!    boundary edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer RTRI, on input, the index of boundary triangle 
!    to begin search at.  On output, the index of rightmost boundary triangle 
!    visible from (X,Y).
!
!    Input/output, integer REDG.  On input, the edge of triangle RTRI that 
!    is visible from (X,Y).  On output, REDG has been updated so that this
!    is still true. 1 <= REDG <= 3.
!
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )
  use precision
  implicit none

  integer, intent(in), dimension(3,*) :: til, tnbr
  integer, intent(inout) :: ledg, ltri, rtri, redg
  integer :: a, b, e, l
  logical :: ldone
  integer :: lr
  integer :: t
  integer, external :: i_wrap, lrline
  real(WP), dimension(2,*) :: vcl
  real(WP) :: x, y

  !  Find rightmost visible boundary edge using links, then possibly
  !  leftmost visible boundary edge using triangle neighbor information.

  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

10 continue

  l = -tnbr(redg,rtri)
  t = l / 3
  e = mod ( l, 3 ) + 1
  a = til(e,t)

  if ( e <= 2 ) then
    b = til(e+1,t)
  else
    b = til(1,t)
  end if

  lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0_WP )

  if ( lr > 0 ) then
    rtri = t
    redg = e
    go to 10
  end if

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)
    e = i_wrap ( e-1, 1, 3 )

    do while ( tnbr(e,t) > 0 )

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), 0.0_WP)

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end subroutine vbedg



!*******************************************************************************
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2-D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to triangulation.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer I, the index in VCL of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input, integer MAXST, the maximum size available for the STACK array.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, double precision VCL(2,*), the coordinates of the vertices.
!
!    Input/output, integer TIL(3,*), the triangle incidence list.  May be updated
!    on output because of swaps.
!
!    Input/output, integer TNBR(3,*), the triangle neighbor list; negative 
!    values are used for links of the counter-clockwise linked list of boundary 
!    edges; May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(1:MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERROR is set to 8 for abnormal return.
!
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierror )
  use precision
  implicit none 

  integer, intent(in) :: maxst
  integer, intent(inout) :: btri, bedg
  integer, external :: diaedg, i_wrap
  integer, dimension(maxst) :: stack
  integer, dimension(3,*), intent(inout) :: til, tnbr

  integer :: a, b, c, e
  integer :: ee, em1, ep1, f, fm1, fp1
  integer :: i, l, r, s, t
  integer :: ierror
  integer :: swap
  integer :: top
  integer :: tt, u
  
  real(WP), dimension(2,*) :: vcl
  real(WP) :: x, y

  !  Determine whether triangles in stack are Delaunay, and swap
  !  diagonal edge of convex quadrilateral if not.

  ierror = 0
  x = vcl(1,i)
  y = vcl(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( til(1,t) == i ) then
      e = 2
      b = til(3,t)
    else if ( til(2,t) == i ) then
      e = 3
      b = til(1,t)
    else
      e = 1
      b = til(2,t)
    end if

    a = til(e,t)
    u = tnbr(e,t)

    if ( tnbr(1,u) == t ) then
      f = 1
      c = til(3,u)
    else if ( tnbr(2,u) == t ) then
      f = 2
      c = til(1,u)
    else
      f = 3
      c = til(2,u)
    end if

    swap = diaedg ( x, y, vcl(1,a), vcl(2,a), vcl(1,c), vcl(2,c), &
      vcl(1,b), vcl(2,b) )

    if ( swap == 1 ) then

      em1 = i_wrap ( e - 1, 1, 3 )
      ep1 = i_wrap ( e + 1, 1, 3 )
      fm1 = i_wrap ( f - 1, 1, 3 )
      fp1 = i_wrap ( f + 1, 1, 3 )
      
      til(ep1,t) = c
      til(fp1,u) = i
      r = tnbr(ep1,t)
      s = tnbr(fp1,u)
      tnbr(ep1,t) = u
      tnbr(fp1,u) = t
      tnbr(e,t) = s
      tnbr(f,u) = r

      if ( tnbr(fm1,u) > 0 ) then
        top = top + 1
        stack(top) = u
      end if

      if ( s > 0 ) then

        if ( tnbr(1,s) == u ) then
          tnbr(1,s) = t
        else if ( tnbr(2,s) == u ) then
          tnbr(2,s) = t
        else
          tnbr(3,s) = t
        end if

        top = top + 1

        if ( top > maxst ) then
          ierror = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == a ) then
            ee = 3
          else if ( til(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

      if ( r > 0 ) then

        if ( tnbr(1,r) == t ) then
          tnbr(1,r) = u
        else if ( tnbr(2,r) == t ) then
          tnbr(2,r) = u
        else
          tnbr(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == b ) then
            ee = 3
          else if ( til(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

    end if

  end do

  return
end subroutine swapec



!*******************************************************************************
!! I_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I_WRAP, a "wrapped" version of IVAL.
!
function i_wrap ( ival, ilo, ihi )
  use precision
  implicit none

  integer, intent(in) :: ival, ihi, ilo
  integer :: i_modp
  integer :: i_wrap
  integer :: wide
  
  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i_wrap = ilo
  else
    i_wrap = ilo + i_modp ( ival-ilo, wide )
  end if

  return
end function i_wrap




!*******************************************************************************
!! DIAEDG chooses one of the diagonals of a quadrilateral.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge 
!    that should be chosen, based on the circumcircle criterion, where 
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple 
!    quadrilateral in counterclockwise order.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, double precision X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )
  use precision
  implicit none

  integer diaedg
  real(WP) :: ca
  real(WP) :: cb
  real(WP) :: dx10
  real(WP) :: dx12
  real(WP) :: dx30
  real(WP) :: dx32
  real(WP) :: dy10
  real(WP) :: dy12
  real(WP) :: dy30
  real(WP) :: dy32
  real(WP) :: s
  real(WP) :: tol
  real(WP) :: tola
  real(WP) :: tolb
  real(WP) :: x0
  real(WP) :: x1
  real(WP) :: x2
  real(WP) :: x3
  real(WP) :: y0
  real(WP) :: y1
  real(WP) :: y2
  real(WP) :: y3

  tol = 100.0_WP * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( ca > tola .and. cb  >  tolb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( s > tola ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end function diaedg




!*******************************************************************************
!! I_MODP returns the nonnegative remainder of integer division.
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  I_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I_MODP, the nonnegative remainder when I is
!    divided by J.
!
function i_modp ( i, j )
  use precision
  implicit none

  integer, intent(in) :: i, j
  integer :: i_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
  end if

  return
end function i_modp
