module math
  use precision

  ! Trigonometric parameters
  real(WP), parameter :: Pi    = 3.1415926535897932385_WP
  real(WP), parameter :: twoPi = 6.2831853071795864770_WP
  
  ! Bessel first zero
  real(WP), parameter :: bessj1_zero = 3.8317059702075123115_WP

  ! Blasius data points
  real(WP), dimension(0:9) :: by0 = (/ &
       0.0_WP, 0.165571818583440_WP, 0.650024518764203_WP, 1.39680822972500_WP, &
       2.30574664618049_WP, 3.28327391871370_WP, 4.27962110517696_WP, &
       5.27923901129384_WP, 6.27921363832835_WP, 7.27921257797747_WP /)
  real(WP), dimension(0:9) :: by1 = (/ &
       0.0_WP, 0.329780063306651_WP, 0.629765721178679_WP, 0.84604458266019_WP, &
       0.95551827831671_WP, 0.99154183259084_WP, 0.99897290050990_WP, &
        0.9999216098795_WP, 0.99999627301467_WP, 0.99999989265063_WP /)
  real(WP), dimension(0:9) :: by2 = (/ &
       0.332057384255589_WP, 0.323007152241930_WP, 0.266751564401387_WP, 0.161360240845588_WP, &
       0.06423404047594_WP, 0.01590689966410_WP, 0.00240199722109_WP, &
       0.00022016340923_WP, 0.00001224984692_WP, 0.00000041090325_WP /)

contains
  
  ! Returns the gamma function
  function gamma(xx)
    implicit none
    
    real(WP) :: gamma
    real(WP), intent(in) :: xx
    
    gamma=exp(gammaln(xx))
    
    return
  end function gamma
  
  
  ! Returns the log of the gamma function
  function gammaln(xx)
    implicit none
  
    real(WP) :: gammaln
    real(WP), intent(in) :: xx
    
    real(WP), parameter :: stp = 2.5066282746310005_WP
    real(WP), dimension(6), parameter :: cof = (/ 76.18009172947146_WP, &
         -86.50532032941677_WP, 24.01409824083091_WP,-1.231739572450155_WP, &
         .1208650973866179E-2_WP, -.5395239384953E-5_WP /)
    
    real(WP) :: ser,tmp,x,y
    integer :: j

    x = xx
    y = x
    tmp = x + 5.5_WP
    tmp = (x+0.5_WP)*log(tmp)-tmp
    ser = 1.000000000190015_WP
    do j=1,6
       y = y + 1.0_WP
       ser = ser+cof(j)/y
    end do
    gammaln = tmp + log(stp*ser/x)
    
    return
  end function gammaln
  

  ! Returns the Bessel function J0(x) for any real x.
  ! [Numerical Recipes in Fortran]
  function bessj0(x)
    implicit none

    real(WP) :: bessj0
    real(WP), intent(in) :: x
    
    real(WP), dimension(6), parameter :: r = (/57568490574.0_WP,-13362590354.0_WP,651619640.7_WP, &
         -11214424.18_WP,77392.33017_WP,-184.9052456_WP/)
    real(WP), dimension(6), parameter :: s = (/57568490411.0_WP,1029532985.0_WP,9494680.718_WP, &
         59272.64853_WP,267.8532712_WP,1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP,-.1098628627E-2_WP,.2734510407E-4_WP, &
         -.2073370639E-5_WP,.2093887211E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/-.1562499995E-1_WP,.1430488765E-3_WP,-.6911147651E-5_WP, &
         .7621095161E-6_WP,-.934945152E-7_WP/)

    real(WP) :: ax,xx,z,y
    
    if(abs(x).lt.8.0_WP)then
       y=x**2
       bessj0=(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) / &
            (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*s(6))))))
    else
       ax=abs(x)
       z=8.0_WP/ax
       y=z**2
       xx=ax-.785398164_WP
       bessj0=sqrt(.636619772_WP/ax)*(cos(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) - &
            z*sin(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))))
    endif
    return
  end function bessj0

  ! Returns the Bessel function J1(x) for any real x.
  ! [Numerical Recipes in Fortran]
  function bessj1(x)
    implicit none

    real(WP) :: bessj1
    real(WP), intent(in) :: x
    
    real(WP), dimension(6), parameter :: r = (/72362614232.0_WP,-7895059235.0_WP,242396853.1_WP, &
         -2972611.439_WP,15704.48260_WP,-30.16036606_WP/)
    real(WP), dimension(6), parameter :: s = (/144725228442.0_WP,2300535178.0_WP,18583304.74_WP, &
         99447.43394_WP,376.9991397_WP,1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP,.183105E-2_WP,-.3516396496E-4_WP,.2457520174E-5_WP, &
         -.240337019E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/.04687499995_WP,-.2002690873E-3_WP,.8449199096E-5_WP, &
         -.88228987E-6_WP,.105787412E-6_WP/)

    real(WP) :: ax,xx,z,y
    
    if(abs(x).lt.8.0_WP)then
       y=x**2
       bessj1=x*(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) / &
            (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*s(6))))))
    else
       ax=abs(x)
       z=8.0_WP/ax
       y=z**2
       xx=ax-2.356194491_WP
       bessj1=sqrt(.636619772_WP/ax)*(cos(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) - &
            z*sin(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))))*sign(1.0_WP,x)
    endif
    return
  end function bessj1
  

  ! Return the Blasius function
  function blasius0(x)
    implicit none

    real(WP) :: blasius0
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius0 = 0.0_WP
    else if (x.ge.9.0_WP) then
       blasius0 = by0(9) + (x-9.0_WP)
    else
       i = int(x)
       f0l = by0(i)
       f0r = by0(i+1)
       f2l = by2(i)
       f2r = by2(i+1)
       blasius0 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius0
  
  ! Return the first derivative of the Blasius function
  function blasius1(x)
    implicit none

    real(WP) :: blasius1
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius1 = 0.0_WP
    else if (x.ge.9.0_WP) then
       blasius1 = 1.0_WP
    else
       i = int(x)
       f0l = by1(i)
       f0r = by1(i+1)
       f2l = -0.5_WP*by0(i)*by2(i)
       f2r = -0.5_WP*by0(i+1)*by2(i+1)
       blasius1 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius1
  
  ! Return the second derivative of the Blasius function
  function blasius2(x)
    implicit none

    real(WP) :: blasius2
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius2 = by2(0)
    else if (x.ge.9.0_WP) then
       blasius2 = 0.0_WP
    else
       i = int(x)
       f0l = by2(i)
       f0r = by2(i+1)
       f2l = 0.25_WP*by0(i)**2*by2(i) - 0.5_WP*by1(i)*by2(i)
       f2r = 0.25_WP*by0(i+1)**2*by2(i+1) - 0.5_WP*by1(i+1)*by2(i+1)
       blasius2 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius2
  
  ! Return the third derivative of the Blasius function
  function blasius3(x)
    implicit none

    real(WP) :: blasius3
    real(WP), intent(in) :: x

    blasius3 = -0.5_WP*blasius0(x)*blasius2(x)
    
    return
  end function blasius3
  
  function arctan(dx,dy)
    implicit none
    
    real(WP), intent(in) :: dx,dy
    real(WP) :: arctan
    
    ! Evaluate atan
    if (abs(dx)+abs(dy).lt.1.0e-9_WP) then
       arctan = 0.0_WP
    else
       arctan = atan(dy/dx)
    end if
    
    ! Quadrant correction
    if (dx.le.0.0_WP) then
       arctan = Pi+arctan
    elseif (dy.le.0.0_WP .and. dx.gt.0.0_WP) then
       arctan = twoPi+arctan
    end if
    
    return
  end function arctan
  
  ! Inverse the linear system : AX=B
  subroutine solve_linear_system(A,B,X,n)
    implicit none

    integer, intent(in) :: n
    real(WP), dimension(n,n) :: A
    real(WP), dimension(n) :: B
    real(WP), dimension(n) :: X
    integer :: i,j,l
    
    ! Forward elimination
    do i=1,n
       ! Pivoting
       j = i-1+maxloc(abs(A(i:n,i)),1)
       X = A(i,:)
       A(i,:) = A(j,:)
       A(j,:) = X
       X(1) = B(i)
       B(i) = B(j)
       B(j) = X(1)
       ! Elimination
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
  
  
  ! Inverse the matrix
  subroutine inverse_matrix(A,B,n)
    implicit none

    integer, intent(in) :: n
    real(WP), dimension(n,n) :: A,B
    integer :: i,l

    B=0.0_WP

    ! Forward elimination
    do i=1,n
       B(i,i) = 1.0_WP
       B(i,:) = B(i,:) / A(i,i)
       A(i,:) = A(i,:) / A(i,i)
       do l=i+1,n
          B(l,:) = B(l,:) - A(l,i)*B(i,:)
          A(l,:) = A(l,:) - A(l,i)*A(i,:)
       end do
    end do

    ! Backward substitution
    do i=n,1,-1
       do l=i+1,n
          B(i,:) = B(i,:) - A(i,l)*B(l,:)
       end do
    end do
    
    return
  end subroutine inverse_matrix
  
  ! Third order finite volume interpolation
  subroutine poly3_coeff(xx,xp,coeff)
    implicit none
    
    real(WP), intent(in),  dimension(4) :: xx
    real(WP), intent(out), dimension(3) :: coeff
    real(WP), intent(in)     :: xp
    real(WP), dimension(3,3) :: A,B
    
    ! Form the matrix
    A(1,1) = (xx(2)**3-xx(1)**3)/(3.0_WP*(xx(2)-xx(1)))
    A(2,1) = (xx(3)**3-xx(2)**3)/(3.0_WP*(xx(3)-xx(2)))
    A(3,1) = (xx(4)**3-xx(3)**3)/(3.0_WP*(xx(4)-xx(3)))
    
    A(1,2) = (xx(2)**2-xx(1)**2)/(2.0_WP*(xx(2)-xx(1)))
    A(2,2) = (xx(3)**2-xx(2)**2)/(2.0_WP*(xx(3)-xx(2)))
    A(3,2) = (xx(4)**2-xx(3)**2)/(2.0_WP*(xx(4)-xx(3)))
    
    A(1,3) = (xx(2)**1-xx(1)**1)/(1.0_WP*(xx(2)-xx(1)))
    A(2,3) = (xx(3)**1-xx(2)**1)/(1.0_WP*(xx(3)-xx(2)))
    A(3,3) = (xx(4)**1-xx(3)**1)/(1.0_WP*(xx(4)-xx(3)))
    
    ! Invert it
    call inverse_matrix(A,B,3)
    
    ! Compute metrics
    coeff(1) = B(1,1)*xp**2+B(2,1)*xp**1+B(3,1)*xp**0
    coeff(2) = B(1,2)*xp**2+B(2,2)*xp**1+B(3,2)*xp**0
    coeff(3) = B(1,3)*xp**2+B(2,3)*xp**1+B(3,3)*xp**0
    
    return
  end subroutine poly3_coeff
  
  ! HO finite volume interpolation
  subroutine hofvi(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n+1) :: xx
    real(WP), intent(out), dimension(n)   :: coeff
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i+1)**(n+1-j)-xx(i)**(n+1-j))/(real(n+1-j,WP)*(xx(i+1)-xx(i)))
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + B(j,i)*xp**(n-j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofvi
  
  ! HO finite difference interpolation
  subroutine hofdi(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xp)**(j-1)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = B(1,:)
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi
  
  ! HO finite difference interpolation
  ! with Dirichlet condition
  subroutine hofdi_d(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + B(j,i)*(xp-xbc)**(j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi_d
  
  ! HO finite difference interpolation
  ! with Neumann condition
  subroutine hofdi_n(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       A(i,1) = 1.0_WP
       do j=2,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       coeff(i) = B(1,i)
       do j=2,n
          coeff(i) = coeff(i) + B(j,i)*(xp-xbc)**(j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi_n
  
  ! HO finite difference differenciation
  subroutine hofdd(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xp)**(j-1)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = B(2,:)
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd
  
  ! HO finite difference differenciation
  ! with Dirichlet condition
  subroutine hofdd_d(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + real(j,WP)*B(j,i)*(xp-xbc)**(j-1)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd_d
  
  ! HO finite difference differenciation
  ! with Neumann condition
  subroutine hofdd_n(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       A(i,1) = 1.0_WP
       do j=2,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=2,n
          coeff(i) = coeff(i) + real(j,WP)*B(j,i)*(xp-xbc)**(j-1)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd_n
  
  ! Eigenvalue/vector solver for a real symmetric NxN matrix A
  ! A is destroyed, d contains the eigenvalues, v contains the eigenvectors in columns
  subroutine eigensolver(A,n,d,v,nrot)
    use precision
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: nrot
    real(WP), dimension(:), intent(out) :: d
    real(WP), dimension(:,:), intent(inout) :: A
    real(WP), dimension(:,:), intent(out) :: v
    integer :: i,ip,iq
    real(WP) :: c,g,h,s,sm,t,tau,theta,tresh
    real(WP), dimension(size(d)) :: b,z
    
    ! Initialize v to be the identity matrix
    v = 0.0_WP
    do ip=1,n
       do iq=1,n
          v(ip,iq) = 0.0_WP
       end do
       v(ip,ip) = 1.0_WP
    end do
    
    ! Initialize b and d to the diagonal of A and z to 0
    do ip=1,n
       b(ip) = A(ip,ip)
       d(ip) = b(ip)
       z(ip) = 0.0_WP
    end do
    
    ! Eigenvalue solver
    nrot = 0
    do i=1,50
       sm=0.0_WP
       do ip=1,n-1
          do iq=ip+1,n
             sm=sm+abs(A(ip,iq))
          end do
       end do
       if (sm.eq.0.0_WP) return
       if (i.lt.4) then
          tresh=0.2_WP*sm/n**2
       else
          tresh=0.0_WP
       end if
       do ip=1,n-1
          do iq=ip+1,n
             g=100.0_WP*abs(A(ip,iq))
             if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
                A(ip,iq)=0.0_WP
             else if (abs(A(ip,iq)).gt.tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g.eq.abs(h)) then
                   t=A(ip,iq)/h
                else
                   theta=0.5_WP*h/A(ip,iq)
                   t=1.0_WP/(abs(theta)+sqrt(1.0_WP+theta**2))
                   if (theta.lt.0.0_WP) t=-t
                end if
                c=1.0_WP/sqrt(1.0_WP+t**2)
                s=t*c
                tau=s/(1.0_WP+c)
                h=t*A(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                A(ip,iq)=0.0_WP
                call jrotate(A(1:ip-1,ip),   A(1:ip-1,iq)   )
                call jrotate(A(ip,ip+1:iq-1),A(ip+1:iq-1,iq))
                call jrotate(A(ip,iq+1:n),   A(iq,iq+1:n)   )
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:)
       z(:)=0.0_WP
    end do
    print*,'Too many iterations in eigensolver'
    
  contains
    
    subroutine jrotate(a1,a2)
      implicit none
      
      real(WP), dimension(:), intent(inout) :: a1,a2
      real(WP), dimension(size(a1)) :: wk1
      
      wk1(:)=a1(:)
      a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
      a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
      
      return
    end subroutine jrotate
    
  end subroutine eigensolver
  
  ! Eigenvalue solver for tri-diag matrices
  ! d contains the diagonal on input, the eigenvalues on output
  ! e contains the lower diag, e(1) is not used. e is destroyed
  subroutine tqli(d,n,e,z)
    use precision
    implicit none
    
    integer, intent(in) :: n
    real(WP), dimension(n), intent(inout) :: d,e
    real(WP), dimension(n,n), intent(inout) :: z
    integer :: i,iter,l,m
    real(WP) :: b,c,dd,f,g,p,r,s
    real(WP), dimension(n) :: ff
    
    ! Renumber e
    do i=2,n
       e(i-1)=e(i)
    end do
    
    ! Do it
    do l=1,n
       iter=0
       iterate: do
          do m=l,n-1
             dd=abs(d(m))+abs(d(m+1))
             if (abs(e(m))+dd.eq.dd) exit
          end do
          if (m.eq.l) exit iterate
          if (iter.eq.30) stop 'tqli: too many iterations.'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0_WP*e(l))
          r=pythag(g,1.0_WP)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0_WP
          c=1.0_WP
          p=0.0_WP
          do i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             if (r.eq.0.0_WP) then
                d(i+1)=d(i+1)-p
                e(m)=0.0_WP
                cycle iterate
             end if
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.0_WP*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             ! Eigenvector
             ff(1:n)=z(1:n,i+1)
             z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
             z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
          end do
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0_WP
       end do iterate
    end do
    
  contains
    
    function pythag(a,b)
      implicit none
      real(WP) :: a,b
      real(WP) :: pythag
      real(WP) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if (absa.gt.absb) then
         pythag=absa*sqrt(1.0_WP+(absb/absa)**2)
      else
         if (absb.eq.0.0_WP) then
            pythag=0.0_WP
         else
            pythag=absb*sqrt(1.0_WP+(absa/absb)**2)
         end if
      end if
      return
    end function pythag
    
  end subroutine tqli
  
end module math
