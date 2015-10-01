module metric_velocity_conv
  use precision
  use geometry
  use partition
  implicit none
  
  ! Stencil lengths
  ! ---------------
  integer :: stc1,stc2
  integer :: stcp
  
  ! Interpolation and derivation stencil
  ! ------------------------------------
  real(WP), dimension(:,:), pointer :: coeffc_deriv
  real(WP), dimension(:,:), pointer :: coeffc_interp  
  real(WP), dimension(:),   pointer :: coeffc_deriv_full
  real(WP), dimension(:),   pointer :: coeffc_interp_full
  
  ! Interpolation operators
  ! -----------------------
  ! Velocity interpolators
  real(WP), dimension(:,:,:), pointer :: interp_Ju_xm,interp_Ju_y, interp_Ju_z
  real(WP), dimension(:,:,:), pointer :: interp_Jv_x, interp_Jv_ym,interp_Jv_z
  real(WP), dimension(:,:,:), pointer :: interp_Jw_x, interp_Jw_y, interp_Jw_zm
  ! Velocity interpolators for cylindrical coordinates
  real(WP), dimension(:,:,:), pointer :: interp_cyl_v_ym,interp_cyl_w_zm
  ! Flux interpolators for cylindrical coordinates
  real(WP), dimension(:,:,:), pointer :: interp_cyl_F_y, interp_cyl_F_z
  
  ! Divergence operators
  ! --------------------
  ! Divergence for continuity
  real(WP), dimension(:,:,:), pointer :: divc_u,divc_v,divc_w
  ! Divergence of the convective fluxes
  real(WP), dimension(:,:,:), pointer :: divc_xx,divc_xy,divc_xz
  real(WP), dimension(:,:,:), pointer :: divc_yx,divc_yy,divc_yz
  real(WP), dimension(:,:,:), pointer :: divc_zx,divc_zy,divc_zz
  ! Length of interpolation used in div
  integer, dimension(:,:,:), pointer :: interp_xx,interp_xy,interp_xz
  integer, dimension(:,:,:), pointer :: interp_yx,interp_yy,interp_yz
  integer, dimension(:,:,:), pointer :: interp_zx,interp_zy,interp_zz
  
  ! Gradient operators
  ! ------------------
  ! Pressure gradient
  real(WP), dimension(:,:,:), pointer :: grad_Px,grad_Py,grad_Pz
  
  ! Laplacian operator
  ! ------------------
  ! Poisson equation
  real(WP), dimension(:,:,:,:,:), pointer :: lap
  
  ! Coefficient for BCs
  ! -------------------
  real(WP), dimension(:,:), pointer :: coeffc_bc
  
contains
  
  ! Allocate all the arrays
  ! -----------------------
  subroutine metric_velocity_conv_allocate
    implicit none
    
    ! Larger because of momentum fluxes
    allocate(interp_Ju_xm(imin_-stc2:imax_+stc1,jmin_-stc2:jmax_+stc1,-stc1:stc2))
    allocate(interp_Ju_y (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    allocate(interp_Ju_z (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    
    allocate(interp_Jv_ym(imin_-stc2:imax_+stc1,jmin_-stc2:jmax_+stc1,-stc1:stc2))
    allocate(interp_Jv_x (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    allocate(interp_Jv_z (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    
    allocate(interp_Jw_zm(imin_-stc2:imax_+stc1,jmin_-stc2:jmax_+stc1,-stc1:stc2))
    allocate(interp_Jw_x (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    allocate(interp_Jw_y (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    
    ! Interior only for fluxes
    allocate(interp_cyl_F_y(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    allocate(interp_cyl_F_z(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    
    ! Larger because of momentum fluxes
    allocate(interp_cyl_v_ym(imin_-stc2:imax_+stc1,jmin_-stc2:jmax_+stc1,-stc1:stc2))
    allocate(interp_cyl_w_zm(imin_-stc2:imax_+stc1,jmin_-stc2:jmax_+stc1,-stc1:stc2))
    
    ! Larger because of momentum fluxes
    allocate(divc_u(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_v(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_w(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    
    ! Only inside the domain
    allocate(divc_xx(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    allocate(divc_xy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_xz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_yx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_yy(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    allocate(divc_yz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_zx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_zy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(divc_zz(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    
    ! Only inside the domain
    allocate(interp_xx(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    allocate(interp_xy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_xz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_yx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_yy(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    allocate(interp_yz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_zx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_zy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
    allocate(interp_zz(imin_:imax_,jmin_:jmax_,-stc2:stc1))
    
    ! Larger for diffusion term in scalar equation
    allocate(grad_Px(imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    allocate(grad_Py(imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    allocate(grad_Pz(imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
    
    ! Global for icc preconditioner
    allocate(lap(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,3,-stcp:stcp))
    
    return
  end subroutine metric_velocity_conv_allocate
  
  
  ! Compute the interpolation coefficients
  ! --------------------------------------
  subroutine metric_velocity_conv_compute_coeffs
    use math
    implicit none

    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:),   pointer :: B
    real(WP), dimension(:),   pointer :: alpha
    integer :: i,l,n
    
    ! Create the pointers to the full stencil
    ! ---------------------------------------
    allocate(coeffc_deriv_full (-stc1:stc2))
    allocate(coeffc_interp_full(-stc1:stc2))
    n = stc2
    
    ! Allocate the matrix/vector
    allocate(A(n,n))
    allocate(B(n))
    allocate(alpha(n))
    
    ! Set up the linear system
    B = 0.0_WP
    B(1) = 1.0_WP
    do i=1,n
       do l=1,n
          A(i,l) = real((2*l-1)**(2*(i-1)),WP)
       end do
    end do
    
    ! Solve the linear system
    call solve_linear_system(A,B,alpha,n)
    
    ! Save the BCs
    do i=1,n
       coeffc_interp_full(1-i) = 0.5_WP * alpha(i)
       coeffc_interp_full(i)   = 0.5_WP * alpha(i)
       coeffc_deriv_full (1-i) = - alpha(i) / real(2*i-1,WP)
       coeffc_deriv_full (i)   = + alpha(i) / real(2*i-1,WP)
    end do
    
    ! Clean up memory
    deallocate(A);nullify(A)
    deallocate(B);nullify(B)
    deallocate(alpha);nullify(alpha)
    
    ! Create the pointers for BCs at walls
    ! ------------------------------------
    allocate(coeffc_deriv (-stc1:stc2,-stc1:stc2))
    allocate(coeffc_interp(-stc1:stc2,-stc1:stc2))
    coeffc_deriv  = 0.0_WP
    coeffc_interp = 0.0_WP
    
    do n=1,stc1
       ! Allocate the matrix/vector
       allocate(A(n,n))
       allocate(B(n))
       allocate(alpha(n))
       
       ! Set up the linear system
       B = 0.0_WP
       B(1) = 1.0_WP
       do i=1,n
          do l=1,n
             A(i,l) = real((2*l-1)**(2*(i-1)),WP)
          end do
       end do
       
       ! Solve the linear system
       call solve_linear_system(A,B,alpha,n)
       
       ! Save the BCs
       do i=1,n
          coeffc_interp(n+1,1-i) = 0.5_WP * alpha(i)
          coeffc_interp(n+1,i)   = 0.5_WP * alpha(i)
          coeffc_deriv (n+1,1-i) = - alpha(i) / real(2*i-1,WP)
          coeffc_deriv (n+1,i)   = + alpha(i) / real(2*i-1,WP)
       end do
    
       ! Clean up memory
       deallocate(A);nullify(A)
       deallocate(B);nullify(B)
       deallocate(alpha);nullify(alpha)
       
    end do
    
    ! Create first order coeffs
    coeffc_interp(1,0) =  0.5_WP
    coeffc_interp(1,1) =  0.5_WP
    coeffc_deriv (1,0) = -1.0_WP
    coeffc_deriv (1,1) = +1.0_WP
    
    ! Symmetrize
    do i=-stc1,0
       do l=-stc1,stc2
          coeffc_deriv (i,l) = - coeffc_deriv (-i+1,-l+1)
          coeffc_interp(i,l) = + coeffc_interp(-i+1,-l+1)
       end do
    end do
    
    ! New coefficients for conservation of mass/momentum
    ! --------------------------------------------------
    allocate(coeffc_bc(-stc1:stc2,-stc1:stc2))    
    
    ! New coefficients for conservation of mass
    select case(vel_conv_order)
    case(2)
       ! Nothing to do for the derivative
    case(4)
       ! One stencil to change
       coeffc_bc(0,-1) = 0.0_WP
       coeffc_bc(0, 0) = coeffc_deriv_full( 0) + 2.0_WP*coeffc_deriv_full(-1)
       coeffc_bc(0,+1) = coeffc_deriv_full(+1) - 1.0_WP*coeffc_deriv_full(-1)
       coeffc_bc(0,+2) = coeffc_deriv_full(+2)
    case(6)
       ! First stencil to change
       coeffc_bc(0,-2) = 0.0_WP
       coeffc_bc(0,-1) = 0.0_WP
       coeffc_bc(0, 0) = coeffc_deriv_full( 0) + 2.0_WP*coeffc_deriv_full(-1) + 1.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(0,+1) = coeffc_deriv_full(+1) - 1.0_WP*coeffc_deriv_full(-1) + 2.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(0,+2) = coeffc_deriv_full(+2)                                - 2.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(0,+3) = coeffc_deriv_full(+3)
       ! Second stencil to change
       coeffc_bc(-1,-2) = 0.0_WP
       coeffc_bc(-1,-1) = coeffc_deriv_full(-1) + 3.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(-1, 0) = coeffc_deriv_full( 0) - 3.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(-1,+1) = coeffc_deriv_full(+1) + 1.0_WP*coeffc_deriv_full(-2)
       coeffc_bc(-1,+2) = coeffc_deriv_full(+2)
       coeffc_bc(-1,+3) = coeffc_deriv_full(+3)
    case default
       !call die('metric_velocity_conv_compute_coeffs: order not coded yet')
    end select
    
    ! For a point on the left
    coeffc_bc(-stc1,:) = coeffc_deriv_full
    
    ! Symmetrize
    do i=1,stc2
       do l=-stc1,stc2
          coeffc_bc(i,l) = - coeffc_bc(-i+1,-l+1)
       end do
    end do
    
    return
  end subroutine metric_velocity_conv_compute_coeffs
  
  
  ! Compute the metric coeffs in physical space
  ! -------------------------------------------
  subroutine metric_velocity_conv_compute_primary
    implicit none
    
    integer :: i,j,st
    
    ! Will be set later in metric_cyl if necessary
    interp_cyl_F_y  = 0.0_WP
    interp_cyl_F_z  = 0.0_WP
    interp_cyl_v_ym = 0.0_WP
    interp_cyl_w_zm = 0.0_WP
    
    ! Interpolation of the velocities at the faces
    do j=jmin_-stc2,jmax_+stc1
       do i=imin_-stc2,imax_+stc1
          interp_Ju_xm(i,j,:) = coeffc_interp_full
          interp_Jv_ym(i,j,:) = coeffc_interp_full
          interp_Jw_zm(i,j,:) = coeffc_interp_full
       end do
    end do
    
    ! Interpolation of the velocities at the corners
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          interp_Ju_y(i,j,:) = coeffc_interp_full * dy(j-stc2:j+stc1) * dymi(j-1)
          interp_Ju_z(i,j,:) = coeffc_interp_full
          
          interp_Jv_x(i,j,:) = coeffc_interp_full * dx(i-stc2:i+stc1) * dxmi(i-1)
          interp_Jv_z(i,j,:) = coeffc_interp_full
          
          interp_Jw_x(i,j,:) = coeffc_interp_full * dx(i-stc2:i+stc1) * dxmi(i-1)
          interp_Jw_y(i,j,:) = coeffc_interp_full * dy(j-stc2:j+stc1) * dymi(j-1)
       end do
    end do
    
    ! Divergence of a vector
    do j=jmin_,jmax_
       do i=imin_,imax_
          divc_u(i,j,:) = coeffc_deriv_full * dxi(i)
          divc_v(i,j,:) = coeffc_deriv_full * dyi(j)
          divc_w(i,j,:) = coeffc_deriv_full * dzi
       end do
    end do
    
    ! Divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          divc_xx(i,j,:) = coeffc_deriv_full * dxmi(i-1)
          divc_xy(i,j,:) = coeffc_deriv_full * dyi(j)
          divc_xz(i,j,:) = coeffc_deriv_full * dzi
          
          divc_yx(i,j,:) = coeffc_deriv_full * dxi(i)
          divc_yy(i,j,:) = coeffc_deriv_full * dymi(j-1)
          divc_yz(i,j,:) = coeffc_deriv_full * dzi
          
          divc_zx(i,j,:) = coeffc_deriv_full * dxi(i)
          divc_zy(i,j,:) = coeffc_deriv_full * dyi(j)
          divc_zz(i,j,:) = coeffc_deriv_full * dzi
       end do
    end do
    
    ! Length of interpolation for divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          do st=-stc2,stc1
             interp_xx(i,j,st) = st
             interp_yy(i,j,st) = st
             interp_zz(i,j,st) = st
          end do
          do st=-stc1,stc2
             interp_xy(i,j,st) = st
             interp_xz(i,j,st) = st
             interp_yx(i,j,st) = st
             interp_yz(i,j,st) = st           
             interp_zx(i,j,st) = st
             interp_zy(i,j,st) = st
          end do
       end do
    end do
    
    ! Gradient of pressure
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          grad_Px(i,j,:) = coeffc_deriv_full * dxmi(i-1)
          grad_Py(i,j,:) = coeffc_deriv_full * dymi(j-1)
          grad_Pz(i,j,:) = coeffc_deriv_full * dzi
       end do
    end do
    
    ! Interpolation of the velocities at the faces
    do j=jmin_-stc2,jmax_+stc1
       do i=imin_-stc2,imax_+stc1
          interp_Jv_ym   (i,j,:) = coeffc_interp_full
          interp_cyl_v_ym(i,j,:) = coeffc_interp_full
          interp_cyl_w_zm(i,j,:) = coeffc_interp_full
       end do
    end do
    
    ! Interpolation of fluxes at the corners
    do j=jmin_,jmax_
       do i=imin_,imax_
          interp_cyl_F_y(i,j,:) = coeffc_interp_full * dy(j-stc2:j+stc1) * dymi(j-1)
          interp_cyl_F_z(i,j,:) = coeffc_interp_full
       end do
    end do
    
    return
  end subroutine metric_velocity_conv_compute_primary
  
  
  ! Compute the metric coeffs in physical space
  ! Force second order in y for cylindrical cases
  ! ---------------------------------------------
  subroutine metric_velocity_conv_enforce_cyl
    use parallel
    implicit none
    
    integer  :: i,j,st
    
    ! If not cylindrical set all interpolations to zero and exit
    if (icyl.ne.1) then
       interp_cyl_v_ym = 0.0_WP
       interp_cyl_w_zm = 0.0_WP
       interp_cyl_F_y  = 0.0_WP
       interp_cyl_F_z  = 0.0_WP
       return
    end if
    
    ! Interpolation of the velocities at the faces
    do j=jmin_-stc2,jmax_+stc1
       do i=imin_-stc2,imax_+stc1
          !interp_Jv_ym(i,j,:) = interp_Jv_ym(i,j,:) * abs(y(j-stc1:j+stc2)) * ymi(j)
          interp_Jv_ym(i,j,:) = interp_Jv_ym(i,j,:) * y(j-stc1:j+stc2) * ymi(j)
       end do
    end do
    
    ! Interpolation of the velocities at the corners
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          !interp_Ju_y(i,j,:) = interp_Ju_y(i,j,:) * abs(ym(j-stc2:j+stc1)) * ymmi(j)
          !interp_Jw_y(i,j,:) = interp_Jw_y(i,j,:) * sign(1.0_WP,ym(j-stc2:j+stc1))
          interp_Ju_y(i,j,:) = interp_Ju_y(i,j,:) * ym(j-stc2:j+stc1) * ymmi(j)
       end do
    end do
    
    ! Interpolation of fluxes at the corners
    do j=jmin_,jmax_
       do i=imin_,imax_
          !interp_cyl_F_y(i,j,:) = interp_cyl_F_y(i,j,:) * sign(1.0_WP,ym(j-stc2:j+stc1))
       end do
    end do
    
    ! Divergence of a vector
    do j=jmin_,jmax_
       do i=imin_,imax_
          !divc_v(i,j,:) = divc_v(i,j,:) * abs(y(j-stc1:j+stc2)) * ymi(j)
          divc_v(i,j,:) = divc_v(i,j,:) * y(j-stc1:j+stc2) * ymi(j)
          divc_w(i,j,:) = divc_w(i,j,:) * ymi(j)
       end do
    end do
    
    ! Divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          divc_xy(i,j,:) = divc_xy(i,j,:) * y(j-stc1:j+stc2)  * ymi(j)
          divc_yy(i,j,:) = divc_yy(i,j,:) * ym(j-stc2:j+stc1) * ymmi(j)
          divc_zy(i,j,:) = divc_zy(i,j,:) * y(j-stc1:j+stc2)  * ymi(j)
          
          divc_xz(i,j,:) = divc_xz(i,j,:) * ymi(j)
          divc_yz(i,j,:) = divc_yz(i,j,:) * ymmi(j)
          divc_zz(i,j,:) = divc_zz(i,j,:) * ymi(j)
       end do
    end do
    
    ! Gradient of pressure
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          grad_Py(i,j,:) = grad_Py(i,j,:) * y(j) * ymmi(j)
          grad_Pz(i,j,:) = grad_Pz(i,j,:) * ymi(j)
       end do
    end do
    
    ! For full energy conservation
    !if (jproc.eq.1) grad_Py(:,jmin,:) = 0.0_WP
    
    ! Mass conservation
    if (jproc.eq.1) then
       do st=-stc1,0
          j = jmin-st
          do i=imin_,imax_
             divc_v(i,j,:) = coeffc_bc(st,:) * dyi(j) * y(j-stc1:j+stc2) * ymi(j)
          end do
       end do
    end if
    
    return
  end subroutine metric_velocity_conv_enforce_cyl
  
  
  ! Change the values of the metric to enforce wall conditions
  ! ----------------------------------------------------------
  subroutine metric_velocity_conv_enforce_walls
    use masks
    implicit none
    
    integer  :: i,j,st
    
    ! Interpolation of the velocities at the faces
    ! -> Set to zero in walls
    do j=jmin_-stc2,jmax_+stc1
       do i=imin_-stc2,imax_+stc1
          do st=-stc1,stc2
             if (mask_u(i+st,j).eq.1) interp_Ju_xm   (i,j,st) = 0.0_WP
             if (mask_v(i,j+st).eq.1) interp_Jv_ym   (i,j,st) = 0.0_WP
             if (mask_v(i,j+st).eq.1) interp_cyl_v_ym(i,j,st) = 0.0_WP
          end do
          if (mask(i,j).eq.1) then
             interp_Ju_xm   (i,j,:) = 0.0_WP
             interp_Jv_ym   (i,j,:) = 0.0_WP
             interp_cyl_v_ym(i,j,:) = 0.0_WP
             interp_cyl_w_zm(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Interpolation of fluxes - cylindrical
    do j=jmin_,jmax_
       do i=imin_,imax_
          if (mask_v(i,j).eq.1) interp_cyl_F_y(i,j,:) = 0.0_WP
          if (mask_w(i,j).eq.1) interp_cyl_F_z(i,j,:) = 0.0_WP
       end do
    end do
    
    ! Interpolation of the velocities at the corners
    ! -> Set to zero in walls
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          do st=-stc2,stc1
             if (mask_u(i,j+st).eq.1) interp_Ju_y(i,j,st) = 0.0_WP
             if (mask_w(i,j+st).eq.1) interp_Jw_y(i,j,st) = 0.0_WP
             if (mask_v(i+st,j).eq.1) interp_Jv_x(i,j,st) = 0.0_WP
             if (mask_w(i+st,j).eq.1) interp_Jw_x(i,j,st) = 0.0_WP
          end do
       end do
    end do
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          if (mask_u(i,j).eq.1) interp_Ju_z(i,j,:) = 0.0_WP
          if (mask_v(i,j).eq.1) interp_Jv_z(i,j,:) = 0.0_WP
       end do
    end do
    do j=jmin_-stc1,jmax_+stc2-1
       do i=imin_-stc1,imax_+stc2-1
          if (mask(i,j).eq.1) then
             interp_Jv_x(i:i+1,j:j+1,:) = 0.0_WP
             interp_Ju_y(i:i+1,j:j+1,:) = 0.0_WP
             interp_Jw_x(i:i+1,j    ,:) = 0.0_WP
             interp_Jw_y(i    ,j:j+1,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Divergence of a vector
    ! -> Conservation of mass
    do j=jmin_,jmax_
       do i=imin_,imax_
          
          do st=-stc1,0
             if (mask_u(i+st,j).eq.1) then
                divc_u(i,j,:)  = coeffc_bc(st,:) * dxi(i)
                divc_u(i,j,st) = 0.0_WP
             end if
             if (mask_v(i,j+st).eq.1) then
                divc_v(i,j,:)  = coeffc_bc(st,:) * dyi(j)
                divc_v(i,j,st) = 0.0_WP
             end if
          end do
          
          do st=stc2,1,-1
             if (mask_u(i+st,j).eq.1) then
                divc_u(i,j,:)  = coeffc_bc(st,:) * dxi(i)
                divc_u(i,j,st) = 0.0_WP
             end if
             if (mask_v(i,j+st).eq.1) then
                divc_v(i,j,:)  = coeffc_bc(st,:) * dyi(j)
                divc_v(i,j,st) = 0.0_WP
             end if
          end do
          
          if (mask(i,j).eq.1) then
             divc_u(i,j,:) = 0.0_WP
             divc_v(i,j,:) = 0.0_WP
             divc_w(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Divergence of a matrix
    ! -> Conservation of momentum    
    do j=jmin_,jmax_
       do i=imin_,imax_
          if (mask_u(i,j).eq.1) then
             divc_xx(i,j,:) = 0.0_WP
             divc_xy(i,j,:) = 0.0_WP
             divc_xz(i,j,:) = 0.0_WP
          else
             do st=dist_xp(i,j),stc1
                divc_xx(i,j,2*dist_xp(i,j)-st-2) = divc_xx(i,j,2*dist_xp(i,j)-st-2) - divc_xx(i,j,st)
                divc_xx(i,j,dist_xp(i,j)-1) = divc_xx(i,j,dist_xp(i,j)-1) + 2.0_WP*divc_xx(i,j,st)
                divc_xx(i,j,st) = 0.0_WP
             end do
             do st=-stc2,-dist_xm(i,j)
                divc_xx(i,j,-2*dist_xm(i,j)-st+2) = divc_xx(i,j,-2*dist_xm(i,j)-st+2) - divc_xx(i,j,st)
                divc_xx(i,j,-dist_xm(i,j)+1) = divc_xx(i,j,-dist_xm(i,j)+1) + 2.0_WP*divc_xx(i,j,st)
                divc_xx(i,j,st) = 0.0_WP
             end do
             do st=dist_yp(i,j)+1,stc2
                divc_xy(i,j,2*dist_yp(i,j)-st) = divc_xy(i,j,2*dist_yp(i,j)-st) - divc_xy(i,j,st)
                divc_xy(i,j,dist_yp(i,j)) = divc_xy(i,j,dist_yp(i,j)) + 2.0_WP*divc_xy(i,j,st)
                divc_xy(i,j,st) = 0.0_WP
             end do
             do st=-stc1,-dist_ym(i,j)
                divc_xy(i,j,-2*dist_ym(i,j)-st+2) = divc_xy(i,j,-2*dist_ym(i,j)-st+2) - divc_xy(i,j,st)
                divc_xy(i,j,-dist_ym(i,j)+1) = divc_xy(i,j,-dist_ym(i,j)+1) + 2.0_WP*divc_xy(i,j,st)
                divc_xy(i,j,st) = 0.0_WP
             end do
          end if

          if (mask_v(i,j).eq.1) then
             divc_yx(i,j,:) = 0.0_WP
             divc_yy(i,j,:) = 0.0_WP
             divc_yz(i,j,:) = 0.0_WP
          else
             do st=dist_xp(i,j)+1,stc2
                divc_yx(i,j,2*dist_xp(i,j)-st) = divc_yx(i,j,2*dist_xp(i,j)-st) - divc_yx(i,j,st)
                divc_yx(i,j,dist_xp(i,j)) = divc_yx(i,j,dist_xp(i,j)) + 2.0_WP*divc_yx(i,j,st)
                divc_yx(i,j,st) = 0.0_WP
             end do
             do st=-stc1,-dist_xm(i,j)
                divc_yx(i,j,-2*dist_xm(i,j)-st+2) = divc_yx(i,j,-2*dist_xm(i,j)-st+2) - divc_yx(i,j,st)
                divc_yx(i,j,-dist_xm(i,j)+1) = divc_yx(i,j,-dist_xm(i,j)+1) + 2.0_WP*divc_yx(i,j,st)
                divc_yx(i,j,st) = 0.0_WP
             end do
             do st=dist_yp(i,j),stc1
                divc_yy(i,j,2*dist_yp(i,j)-st-2) = divc_yy(i,j,2*dist_yp(i,j)-st-2) - divc_yy(i,j,st)
                divc_yy(i,j,dist_yp(i,j)-1) = divc_yy(i,j,dist_yp(i,j)-1) + 2.0_WP*divc_yy(i,j,st)
                divc_yy(i,j,st) = 0.0_WP
             end do
             do st=-stc2,-dist_ym(i,j)
                divc_yy(i,j,-2*dist_ym(i,j)-st+2) = divc_yy(i,j,-2*dist_ym(i,j)-st+2) - divc_yy(i,j,st)
                divc_yy(i,j,-dist_ym(i,j)+1) = divc_yy(i,j,-dist_ym(i,j)+1) + 2.0_WP*divc_yy(i,j,st)
                divc_yy(i,j,st) = 0.0_WP
             end do
          end if

          if (mask_w(i,j).eq.1) then
             divc_zx(i,j,:) = 0.0_WP
             divc_zy(i,j,:) = 0.0_WP
             divc_zz(i,j,:) = 0.0_WP
          else
             do st=dist_xp(i,j)+1,stc2
                divc_zx(i,j,2*dist_xp(i,j)-st) = divc_zx(i,j,2*dist_xp(i,j)-st) - divc_zx(i,j,st)
                divc_zx(i,j,dist_xp(i,j)) = divc_zx(i,j,dist_xp(i,j)) + 2.0_WP*divc_zx(i,j,st)
                divc_zx(i,j,st) = 0.0_WP
             end do
             do st=-stc1,-dist_xm(i,j)
                divc_zx(i,j,-2*dist_xm(i,j)-st+2) = divc_zx(i,j,-2*dist_xm(i,j)-st+2) - divc_zx(i,j,st)
                divc_zx(i,j,-dist_xm(i,j)+1) = divc_zx(i,j,-dist_xm(i,j)+1) + 2.0_WP*divc_zx(i,j,st)
                divc_zx(i,j,st) = 0.0_WP
             end do
             do st=dist_yp(i,j)+1,stc2
                divc_zy(i,j,2*dist_yp(i,j)-st) = divc_zy(i,j,2*dist_yp(i,j)-st) - divc_zy(i,j,st)
                divc_zy(i,j,dist_yp(i,j)) = divc_zy(i,j,dist_yp(i,j)) + 2.0_WP*divc_zy(i,j,st)
                divc_zy(i,j,st) = 0.0_WP
             end do
             do st=-stc1,-dist_ym(i,j)
                divc_zy(i,j,-2*dist_ym(i,j)-st+2) = divc_zy(i,j,-2*dist_ym(i,j)-st+2) - divc_zy(i,j,st)
                divc_zy(i,j,-dist_ym(i,j)+1) = divc_zy(i,j,-dist_ym(i,j)+1) + 2.0_WP*divc_zy(i,j,st)
                divc_zy(i,j,st) = 0.0_WP
             end do
          end if

       end do
    end do

    ! Length of interpolation for divergence of a matrix
    ! -> Conservation of momentum    
    do j=jmin_,jmax_
       do i=imin_,imax_

          do st=-stc2,stc1
             if (st.le.-dist_xm(i+st,j))   interp_xx(i,j,st) = -dist_xm(i+st,j)
             if (st.ge. dist_xp(i+st,j)-1) interp_xx(i,j,st) =  dist_xp(i+st,j)-1
             
             if (st.le.-dist_ym(i,j+st))   interp_yy(i,j,st) = -dist_ym(i,j+st)
             if (st.ge. dist_yp(i,j+st)-1) interp_yy(i,j,st) =  dist_yp(i,j+st)-1
          end do
          
          do st=-stc1,stc2
             if (st.lt.-dist_ym(i,j+st-1)+1) then
                interp_xy(i,j,st) = -dist_ym(i,j+st-1)+1
                interp_zy(i,j,st) = -dist_ym(i,j+st-1)+1
             end if
             if (st.gt. dist_yp(i,j+st)) then
                interp_xy(i,j,st) =  dist_yp(i,j+st)
                interp_zy(i,j,st) =  dist_yp(i,j+st)
             end if
             
             if (st.lt.-dist_xm(i+st-1,j)+1) then
                interp_yx(i,j,st) = -dist_xm(i+st-1,j)+1
                interp_zx(i,j,st) = -dist_xm(i+st-1,j)+1
             end if
             if (st.gt. dist_xp(i+st,j)) then
                interp_yx(i,j,st) =  dist_xp(i+st,j)
                interp_zx(i,j,st) =  dist_xp(i+st,j)
             end if
          end do
       end do
    end do
    
    ! Gradient of pressure
    ! -> Degenerate the order
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          if (mask_u(i,j).eq.1) then
             grad_Px(i,j,:) = 0.0_WP
          else
             do st=dist_xp(i,j),stc1
                grad_Px(i,j,2*dist_xp(i,j)-st-2) = grad_Px(i,j,2*dist_xp(i,j)-st-2) - grad_Px(i,j,st)
                grad_Px(i,j,dist_xp(i,j)-1) = grad_Px(i,j,dist_xp(i,j)-1) + 2.0_WP*grad_Px(i,j,st)
                grad_Px(i,j,st) = 0.0_WP
             end do
             do st=-stc2,-dist_xm(i,j)
                grad_Px(i,j,-2*dist_xm(i,j)-st+2) = grad_Px(i,j,-2*dist_xm(i,j)-st+2) - grad_Px(i,j,st)
                grad_Px(i,j,-dist_xm(i,j)+1) = grad_Px(i,j,-dist_xm(i,j)+1) + 2.0_WP*grad_Px(i,j,st)
                grad_Px(i,j,st) = 0.0_WP
             end do
          end if
          if (mask_v(i,j).eq.1) then
             grad_Py(i,j,:) = 0.0_WP
          else
             do st=dist_yp(i,j),stc1
                grad_Py(i,j,2*dist_yp(i,j)-st-2) = grad_Py(i,j,2*dist_yp(i,j)-st-2) - grad_Py(i,j,st)
                grad_Py(i,j,dist_yp(i,j)-1) = grad_Py(i,j,dist_yp(i,j)-1) + 2.0_WP*grad_Py(i,j,st)
                grad_Py(i,j,st) = 0.0_WP
             end do
             do st=-stc2,-dist_ym(i,j)
                grad_Py(i,j,-2*dist_ym(i,j)-st+2) = grad_Py(i,j,-2*dist_ym(i,j)-st+2) - grad_Py(i,j,st)
                grad_Py(i,j,-dist_ym(i,j)+1) = grad_Py(i,j,-dist_ym(i,j)+1) + 2.0_WP*grad_Py(i,j,st)
                grad_Py(i,j,st) = 0.0_WP
             end do
          end if
          if (mask_w(i,j).eq.1) grad_Pz(i,j,:) = 0.0_WP
       end do
    end do
    
    return
  end subroutine metric_velocity_conv_enforce_walls


  ! Enforce the physical boundary conditions of the domain
  ! ------------------------------------------------------
  subroutine metric_velocity_conv_enforce_bc
    use parallel
    implicit none
    integer :: i,j,st
    
    if (xper.ne.1) then
       ! Left boundary
       ! -> Newmann on P
       ! -> Dirichlet on U,V,W
       if (iproc.eq.1) then
          divc_xx(imin,:,:) = 0.0_WP
          divc_xy(imin,:,:) = 0.0_WP
          divc_xz(imin,:,:) = 0.0_WP
          
          ! Mass conservation
          do st=-stc1,0
             i = imin-st
             do j=jmin_,jmax_
                divc_u(i,j,:) = coeffc_bc(st,:) * dxi(i)
             end do
          end do
          
          ! Momentum conservation
          do i=imin-stc1,imin
             grad_Px(i,:,:) = 0.0_WP
          end do
          do i=imin+1,imax_
             do st=-stc2,imin-1-i
                divc_xx(i,:,-2*(i-imin)-st) = divc_xx(i,:,-2*(i-imin)-st) - divc_xx(i,:,st)
                divc_xx(i,:,imin-i) = divc_xx(i,:,imin-i) + 2.0_WP*divc_xx(i,:,st)
                divc_xx(i,:,st) = 0.0_WP
                grad_Px(i,:,-2*(i-imin)-st) = grad_Px(i,:,-2*(i-imin)-st) - grad_Px(i,:,st)
                grad_Px(i,:,imin-i) = grad_Px(i,:,imin-i) + 2.0_WP*grad_Px(i,:,st)
                grad_Px(i,:,st) = 0.0_WP
             end do
          end do
          do i=imin,imax_
             do st=-stc1,imin-1-i
                divc_yx(i,:,-2*(i-imin)-st) = divc_yx(i,:,-2*(i-imin)-st) - divc_yx(i,:,st)
                divc_yx(i,:,imin-i) = divc_yx(i,:,imin-i) + 2.0_WP*divc_yx(i,:,st)
                divc_yx(i,:,st) = 0.0_WP
                divc_zx(i,:,-2*(i-imin)-st) = divc_zx(i,:,-2*(i-imin)-st) - divc_zx(i,:,st)
                divc_zx(i,:,imin-i) = divc_zx(i,:,imin-i) + 2.0_WP*divc_zx(i,:,st)
                divc_zx(i,:,st) = 0.0_WP
             end do
             do st=-stc2,stc1
                if (2*st.le.imin-i-1) interp_xx(i,:,st) = -i+imin-st-1
             end do
             do st=-stc1,stc2
                if (2*st.lt.imin-i+1) then
                   interp_yx(i,:,st) = -i+imin-st+1
                   interp_zx(i,:,st) = -i+imin-st+1
                end if
             end do
          end do

          ! Imposed U velocity
          do i=imin-stc2,imin-1
             interp_Ju_xm(i,:,:) = 0.0_WP
          end do
          do i=imin,imin+stc1-1
             do st=-stc1,imin-i-1
                interp_Ju_xm(i,:,imin-i) = interp_Ju_xm(i,:,imin-i) + interp_Ju_xm(i,:,st)
                interp_Ju_xm(i,:,st) = 0.0_WP
             end do
          end do
          
          ! Imposed V & W velocity
          do i=imin-stc1,imin-1
             interp_Jv_x(i,:,:) = 0.0_WP
             interp_Jw_x(i,:,:) = 0.0_WP
          end do
          do i=imin,imin+stc1-1
             do st=-stc2,imin-2-i
                interp_Jv_x(i,:,imin-1-i) = interp_Jv_x(i,:,imin-1-i) + interp_Jv_x(i,:,st)
                interp_Jw_x(i,:,imin-1-i) = interp_Jw_x(i,:,imin-1-i) + interp_Jw_x(i,:,st)
                interp_Jv_x(i,:,st) = 0.0_WP
                interp_Jw_x(i,:,st) = 0.0_WP
             end do
          end do
       end if

       ! Right boundary
       ! -> Newmann on P
       ! -> Dirichlet on U,V,W
       if (iproc.eq.npx) then
          
          ! Mass conservation
          do st=1,stc2
             i = imax+1-st
             do j=jmin_,jmax_
                divc_u(i,j,:) = coeffc_bc(st,:) * dxi(i)
             end do
          end do
          
          ! Momentum conservation
          do i=imax+1,imax+stc2
             grad_Px(i,:,:) = 0.0_WP
          end do
          do i=imin_,imax
             do st=imax-i+1,stc1
                divc_xx(i,:,2*(imax-i)-st) = divc_xx(i,:,2*(imax-i)-st) - divc_xx(i,:,st)
                divc_xx(i,:,imax-i) = divc_xx(i,:,imax-i) + 2.0_WP*divc_xx(i,:,st)
                divc_xx(i,:,st) = 0.0_WP
                grad_Px(i,:,2*(imax-i)-st) = grad_Px(i,:,2*(imax-i)-st) - grad_Px(i,:,st)
                grad_Px(i,:,imax-i) = grad_Px(i,:,imax-i) + 2.0_WP*grad_Px(i,:,st)
                grad_Px(i,:,st) = 0.0_WP
             end do          
             do st=imax-i+2,stc2
                divc_yx(i,:,2*(imax-i)-st+2) = divc_yx(i,:,2*(imax-i)-st+2) - divc_yx(i,:,st)
                divc_yx(i,:,imax-i+1) = divc_yx(i,:,imax-i+1) + 2.0_WP*divc_yx(i,:,st)
                divc_yx(i,:,st) = 0.0_WP
                divc_zx(i,:,2*(imax-i)-st+2) = divc_zx(i,:,2*(imax-i)-st+2) - divc_zx(i,:,st)
                divc_zx(i,:,imax-i+1) = divc_zx(i,:,imax-i+1) + 2.0_WP*divc_zx(i,:,st)
                divc_zx(i,:,st) = 0.0_WP
             end do
             do st=-stc2,stc1
                if (2*st.ge.imax-i) interp_xx(i,:,st) = imax-i-st
             end do
             do st=-stc1,stc2
                if (2*st.gt.imax-i+1) then
                   interp_yx(i,:,st) = imax-i-st+1
                   interp_zx(i,:,st) = imax-i-st+1
                end if
             end do
          end do

          ! Imposed U velocity
          do i=imax+1,imax+stc1
             interp_Ju_xm(i,:,:) = 0.0_WP
          end do
          do i=imax+1-stc1,imax
             do st=imax+2-i,stc2
                interp_Ju_xm(i,:,imax+1-i) = interp_Ju_xm(i,:,imax+1-i) + interp_Ju_xm(i,:,st)
                interp_Ju_xm(i,:,st) = 0.0_WP
             end do
          end do
          
          ! Imposed V & W velocity
          do i=imax+2,imax+stc2
             interp_Jv_x(i,:,:) = 0.0_WP
             interp_Jw_x(i,:,:) = 0.0_WP
          end do
          do i=imax+2-stc1,imax+1
             do st=imax+2-i,stc1
                interp_Jv_x(i,:,imax+1-i) = interp_Jv_x(i,:,imax+1-i) + interp_Jv_x(i,:,st)
                interp_Jw_x(i,:,imax+1-i) = interp_Jw_x(i,:,imax+1-i) + interp_Jw_x(i,:,st)
                interp_Jv_x(i,:,st) = 0.0_WP
                interp_Jw_x(i,:,st) = 0.0_WP
             end do
          end do
       end if
    end if
    
    if (yper.ne.1) then
       ! Lower boundary
       ! -> Newmann on P
       ! -> Newmann on U/W
       ! -> Dirichlet on V
       if (jproc.eq.1 .and. icyl.eq.0) then
          divc_yx(:,jmin,:) = 0.0_WP
          divc_yy(:,jmin,:) = 0.0_WP
          divc_yz(:,jmin,:) = 0.0_WP
          
          ! Mass conservation
          do st=-stc1,0
             j = jmin-st
             do i=imin_,imax_
                divc_v(i,j,:) = coeffc_bc(st,:) * dyi(j)
                divc_v(i,j,st) = 0.0_WP
             end do
          end do
          
          ! Momentum conservation
          do j=jmin-stc1,jmin
             grad_Py(:,j,:) = 0.0_WP
          end do
          do j=jmin+1,jmax_
             do st=-stc2,jmin-1-j
                divc_yy(:,j,-2*(j-jmin)-st) = divc_yy(:,j,-2*(j-jmin)-st) - divc_yy(:,j,st)
                divc_yy(:,j,jmin-j) = divc_yy(:,j,jmin-j) + 2.0_WP*divc_yy(:,j,st)
                divc_yy(:,j,st) = 0.0_WP
                grad_Py(:,j,2*(jmin-j)-st) = grad_Py(:,j,2*(jmin-j)-st) - grad_Py(:,j,st)
                grad_Py(:,j,jmin-j) = grad_Py(:,j,jmin-j) + 2.0_WP*grad_Py(:,j,st)
                grad_Py(:,j,st) = 0.0_WP
             end do
          end do
          do j=jmin,jmax_
             do st=-stc1,jmin-1-j
                divc_xy(:,j,-2*(j-jmin)-st) = divc_xy(:,j,-2*(j-jmin)-st) - divc_xy(:,j,st)
                divc_xy(:,j,jmin-j) = divc_xy(:,j,jmin-j) + 2.0_WP*divc_xy(:,j,st)
                divc_xy(:,j,st) = 0.0_WP
                divc_zy(:,j,-2*(j-jmin)-st) = divc_zy(:,j,-2*(j-jmin)-st) - divc_zy(:,j,st)
                divc_zy(:,j,jmin-j) = divc_zy(:,j,jmin-j) + 2.0_WP*divc_zy(:,j,st)
                divc_zy(:,j,st) = 0.0_WP
             end do
             do st=-stc2,stc1
                if (2*st.le.jmin-j-1) interp_yy(:,j,st) = -j+jmin-st-1
             end do
             do st=-stc1,stc2
                if (2*st.lt.jmin-j+1) then
                   interp_xy(:,j,st) = -j+jmin-st+1
                   interp_zy(:,j,st) = -j+jmin-st+1
                end if
             end do
          end do
          
          ! Dirichlet condition on V
          do j=jmin-stc2,jmin-1
             interp_Jv_ym(:,j,:) = 0.0_WP
          end do
          do j=jmin,jmin+stc1
             do st=-stc1,jmin-j
!!$                interp_Jv_ym(:,j,jmin-j) = interp_Jv_ym(:,j,jmin-j) + interp_Jv_ym(:,j,st) ! Neumann
                interp_Jv_ym(:,j,st) = 0.0_WP
             end do
          end do
          
          ! Newman condition on U,W
          do j=jmin-stc1,jmin-1
             interp_Ju_y(:,j,:) = 0.0_WP
             interp_Jw_y(:,j,:) = 0.0_WP
          end do
          do j=jmin,jmin+stc1
             do st=-stc2,jmin-1-j
                interp_Ju_y(:,j,jmin-j) = interp_Ju_y(:,j,jmin-j) + interp_Ju_y(:,j,st)
                interp_Jw_y(:,j,jmin-j) = interp_Jw_y(:,j,jmin-j) + interp_Jw_y(:,j,st)
                interp_Ju_y(:,j,st) = 0.0_WP
                interp_Jw_y(:,j,st) = 0.0_WP
             end do
          end do
       end if

       ! Upper boundary
       ! -> Newmann on P
       ! -> Newmann on U/W
       ! -> Dirichlet on V
       if (jproc.eq.npy) then
          
          ! Mass conservation
          do st=1,stc2
             j = jmax+1-st
             do i=imin_,imax_
                divc_v(i,j,:) = coeffc_bc(st,:) * dyi(j)
                divc_v(i,j,st) = 0.0_WP
             end do
          end do
          
          ! Momentum conservation
          do j=jmax+1,jmax+stc2
             grad_Py(:,j,:) = 0.0_WP
          end do
          do j=jmin_,jmax
             do st=jmax-j+1,stc1
                divc_yy(:,j,2*(jmax-j)-st) = divc_yy(:,j,2*(jmax-j)-st) - divc_yy(:,j,st)
                divc_yy(:,j,jmax-j) = divc_yy(:,j,jmax-j) + 2.0_WP*divc_yy(:,j,st)
                divc_yy(:,j,st) = 0.0_WP
                grad_Py(:,j,2*(jmax-j)-st) = grad_Py(:,j,2*(jmax-j)-st) - grad_Py(:,j,st)
                grad_Py(:,j,jmax-j) = grad_Py(:,j,jmax-j) + 2.0_WP*grad_Py(:,j,st)
                grad_Py(:,j,st) = 0.0_WP
             end do          
             do st=jmax-j+2,stc2
                divc_xy(:,j,2*(jmax-j)-st+2) = divc_xy(:,j,2*(jmax-j)-st+2) - divc_xy(:,j,st)
                divc_xy(:,j,jmax-j+1) = divc_xy(:,j,jmax-j+1) + 2.0_WP*divc_xy(:,j,st)
                divc_xy(:,j,st) = 0.0_WP
                divc_zy(:,j,2*(jmax-j)-st+2) = divc_zy(:,j,2*(jmax-j)-st+2) - divc_zy(:,j,st)
                divc_zy(:,j,jmax-j+1) = divc_zy(:,j,jmax-j+1) + 2.0_WP*divc_zy(:,j,st)
                divc_zy(:,j,st) = 0.0_WP
             end do
             do st=-stc2,stc1
                if (2*st.ge.jmax-j) interp_yy(:,j,st) = jmax-j-st
             end do
             do st=-stc1,stc2
                if (2*st.gt.jmax-j+1) then
                   interp_xy(:,j,st) = jmax-j-st+1
                   interp_zy(:,j,st) = jmax-j-st+1
                end if
             end do
          end do

          ! Dirichlet condition on V
          do j=jmax+1,jmax+stc1
             interp_Jv_ym   (:,j,:) = 0.0_WP
             interp_cyl_v_ym(:,j,:) = 0.0_WP
          end do
          do j=jmax-stc1,jmax
             do st=jmax+1-j,stc2
!!$                interp_Jv_ym   (:,j,jmax-j) = interp_Jv_ym   (:,j,jmax-j) + interp_Jv_ym   (:,j,st) ! Neumann 
!!$                interp_cyl_v_ym(:,j,jmax-j) = interp_cyl_v_ym(:,j,jmax-j) + interp_cyl_v_ym(:,j,st) ! Neumann
                interp_Jv_ym   (:,j,st) = 0.0_WP
                interp_cyl_v_ym(:,j,st) = 0.0_WP
             end do
          end do
          
          ! Newman condition on U,W
          do j=jmax+1,jmax+stc2
             interp_Ju_y(:,j,:) = 0.0_WP
             interp_Jw_y(:,j,:) = 0.0_WP
          end do
          do j=jmax+1-stc1,jmax
             do st=jmax+1-j,stc1
                interp_Ju_y(:,j,jmax-j) = interp_Ju_y(:,j,jmax-j) + interp_Ju_y(:,j,st)
                interp_Jw_y(:,j,jmax-j) = interp_Jw_y(:,j,jmax-j) + interp_Jw_y(:,j,st)
                interp_Ju_y(:,j,st) = 0.0_WP
                interp_Jw_y(:,j,st) = 0.0_WP
             end do
          end do
       end if
    end if
    
    return
  end subroutine metric_velocity_conv_enforce_bc
  
  
  ! Enforce the dimension of the problem
  ! ------------------------------------
  subroutine metric_velocity_conv_enforce_dim
    use parallel
    implicit none

    ! Remove dependance in x and U
    if (nx.eq.1) then
       interp_Ju_xm = 0.0_WP
       interp_Ju_y  = 0.0_WP
       interp_Ju_z  = 0.0_WP
       
       divc_u = 0.0_WP
       
       divc_xx = 0.0_WP
       divc_yx = 0.0_WP
       divc_zx = 0.0_WP
       divc_xy = 0.0_WP
       divc_xz = 0.0_WP
       
       grad_Px = 0.0_WP
    end if

    ! Remove dependance in y and V
    if (ny.eq.1) then
       interp_Jv_ym = 0.0_WP
       interp_Jv_x  = 0.0_WP
       interp_Jv_z  = 0.0_WP
       
       interp_cyl_F_y  = 0.0_WP
       interp_cyl_v_ym = 0.0_WP
       
       divc_v = 0.0_WP
       
       divc_xy = 0.0_WP
       divc_yy = 0.0_WP
       divc_zy = 0.0_WP
       divc_yx = 0.0_WP
       divc_yz = 0.0_WP
       
       grad_Py = 0.0_WP
    end if

    ! Remove dependance in z and W
    if (nz.eq.1) then
       interp_Jw_zm = 0.0_WP
       interp_Jw_y  = 0.0_WP
       interp_Jw_x  = 0.0_WP
       
       interp_cyl_F_z  = 0.0_WP
       interp_cyl_w_zm = 0.0_WP
       
       divc_w = 0.0_WP
       
       divc_xz = 0.0_WP
       divc_yz = 0.0_WP
       divc_zx = 0.0_WP
       divc_zy = 0.0_WP
       divc_zz = 0.0_WP
       
       grad_Pz = 0.0_WP
       
       if (icyl.eq.1 .and. jproc.eq.1) then
          divc_yx(:,jmin,:) = 0.0_WP
          divc_yy(:,jmin,:) = 0.0_WP
          divc_yz(:,jmin,:) = 0.0_WP
          interp_cyl_F_y(:,jmin,:) = 0.0_WP
       end if
    end if
    
    ! Handle sector computation
    if (icyl.eq.1 .and. jproc.eq.1 .and. isect.eq.1) then
       divc_yx(:,jmin,:) = 0.0_WP
       divc_yy(:,jmin,:) = 0.0_WP
       divc_yz(:,jmin,:) = 0.0_WP
       interp_cyl_F_y(:,jmin,:) = 0.0_WP
    end if

    return
  end subroutine metric_velocity_conv_enforce_dim

  ! Enforce the Dirichlet conditions
  ! --------------------------------
  subroutine metric_velocity_conv_dirichlet
    use parallel
    use masks
    implicit none
    
    integer :: i,j,st
    
    do j=jmin_,jmax_
       do i=imin_,imax_

          ! Mass conservation
          if (mask(i,j).eq.0) then
             do st=-stc1,0
                if (mask_u(i+st,j).eq.3) divc_u(i,j,:)  = coeffc_bc(st,:) * dxi(i)
                if (mask_v(i,j+st).eq.3) divc_v(i,j,:)  = coeffc_bc(st,:) * dyi(j)
             end do
             do st=stc2,1,-1
                if (mask_u(i+st,j).eq.3) divc_u(i,j,:)  = coeffc_bc(st,:) * dxi(i)
                if (mask_v(i,j+st).eq.3) divc_v(i,j,:)  = coeffc_bc(st,:) * dyi(j)
             end do
          end if
          
          ! Remove dependance in U
          if (mask_u(i,j).eq.3) then
             divc_xx(i,j,:) = 0.0_WP
             divc_xy(i,j,:) = 0.0_WP
             divc_xz(i,j,:) = 0.0_WP
             
             grad_Px(i,j,:) = 0.0_WP
          end if

          ! Remove dependance in V
          if (mask_v(i,j).eq.3) then
             interp_cyl_F_y(i,j,:) = 0.0_WP

             divc_yx(i,j,:) = 0.0_WP
             divc_yy(i,j,:) = 0.0_WP
             divc_yz(i,j,:) = 0.0_WP
             
             grad_Py(i,j,:) = 0.0_WP
          end if
          
          ! Remove dependance in W
          if (mask_w(i,j).eq.3) then
             interp_cyl_F_z(i,j,:) = 0.0_WP
             
             divc_zx(i,j,:) = 0.0_WP
             divc_zy(i,j,:) = 0.0_WP
             divc_zz(i,j,:) = 0.0_WP
             
             grad_Pz(i,j,:) = 0.0_WP
          end if

          ! Take care of the full cell
          if (mask(i,j).eq.3) then
             divc_u(i,j,:) = 0.0_WP
             divc_v(i,j,:) = 0.0_WP
             divc_w(i,j,:) = 0.0_WP
          end if
         
       end do
    end do
    
    ! Gradient of pressure
    ! -> Degenerate the order
    do j=jmin_-stc1,jmax_+stc2
       do i=imin_-stc1,imax_+stc2
          do st=-stc2,-1
             if (mask(i+st,j).eq.3) then
                grad_Px(i,j,st+1) = grad_Px(i,j,st+1) + grad_Px(i,j,st)
                grad_Px(i,j,st) = 0.0_WP
             end if
             if (mask(i,j+st).eq.3) then
                grad_Py(i,j,st+1) = grad_Py(i,j,st+1) + grad_Py(i,j,st)
                grad_Py(i,j,st) = 0.0_WP
             end if
          end do
          do st=stc1,0,-1
             if (mask(i+st,j).eq.3) then
                grad_Px(i,j,st-1) = grad_Px(i,j,st-1) + grad_Px(i,j,st)
                grad_Px(i,j,st) = 0.0_WP
             end if
             if (mask(i,j+st).eq.3) then
                grad_Py(i,j,st-1) = grad_Py(i,j,st-1) + grad_Py(i,j,st)
                grad_Py(i,j,st) = 0.0_WP
             end if
          end do
       end do
    end do
    
    return
  end subroutine metric_velocity_conv_dirichlet

  ! Compute the secondary metric dependant on primary metric
  ! --------------------------------------------------------
  subroutine metric_velocity_conv_compute_secondary
    implicit none
    integer :: i,j,k
    integer :: s1,s2
    
    ! Initialize
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             lap(i,j,k,1,:) = 0.0_WP
             lap(i,j,k,2,:) = 0.0_WP
             lap(i,j,k,3,:) = 0.0_WP
          end do
       end do
    end do
    
    ! Laplacian operator: lap = div( grad )
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             ! Convolute
             do s1=-stc1,stc2
                do s2=-stc2,stc1
                   lap(i,j,k,1,s1+s2) = lap(i,j,k,1,s1+s2) + &
                        divc_u(i,j,s1)*grad_Px(i+s1,j,s2)
                   lap(i,j,k,2,s1+s2) = lap(i,j,k,2,s1+s2) + &
                        divc_v(i,j,s1)*grad_Py(i,j+s1,s2)
                   lap(i,j,k,3,s1+s2) = lap(i,j,k,3,s1+s2) + &
                        divc_w(i,j,s1)*grad_Pz(i,j,s2)
                end do
             end do
          end do
       end do
    end do
    
    do s1=-stcp,stcp
       call boundary_update_border(lap(:,:,:,1,s1),'+','ym')
       call boundary_update_border(lap(:,:,:,2,s1),'+','ym')
       call boundary_update_border(lap(:,:,:,3,s1),'+','ym')
    end do
    
    return
  end subroutine metric_velocity_conv_compute_secondary
  
end module metric_velocity_conv


subroutine metric_velocity_conv_init
  use metric_velocity_conv
  implicit none
  
  ! Set the stencil lengths
  stc2 = vel_conv_order/2
  stc1 = stc2 - 1
  stcp = stc1 + stc2
  
  ! Allocate all the arrays
  call metric_velocity_conv_allocate

  ! Compute the interpolation weights
  call metric_velocity_conv_compute_coeffs

  ! Compute the metric
  call metric_velocity_conv_compute_primary

  ! Impose BCs
  call metric_velocity_conv_enforce_bc
  call metric_velocity_conv_enforce_walls

  ! Enforce some properties
  call metric_velocity_conv_enforce_cyl
  call metric_velocity_conv_enforce_dim
  call metric_velocity_conv_dirichlet
  
  ! Compute the secondary metrics
  call metric_velocity_conv_compute_secondary

  return
end subroutine metric_velocity_conv_init
