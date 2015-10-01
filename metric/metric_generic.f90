module metric_generic
  use precision
  use geometry
  use partition
  implicit none
  
  ! Stencil lengths
  ! ---------------
  integer, parameter :: generic_order=2
  integer :: st1,st2
  integer :: stp
  
  ! Interpolation and derivation stencil
  ! ------------------------------------
  real(WP), dimension(:),   pointer :: coeff_deriv
  real(WP), dimension(:),   pointer :: coeff_interp  
  real(WP), dimension(:,:), pointer :: coeff_interp_2D
  
  ! Interpolation operators
  ! -----------------------
  ! Velocity interpolation
  real(WP), dimension(:,:,:), pointer :: interp_u_xm,interp_v_ym,interp_w_zm
  real(WP), dimension(:,:,:), pointer :: interp_uvw_x,interp_uvw_y,interp_uvw_z
  ! Scalar interpolation
  real(WP), dimension(:,:,:),   pointer :: interp_sc_x,interp_sc_y,interp_sc_z
  real(WP), dimension(:,:,:,:), pointer :: interp_sc_xy,interp_sc_yz,interp_sc_xz
  
  ! Divergence operator
  ! -------------------
  real(WP), dimension(:,:,:), pointer :: div_u,div_v,div_w
  
  ! Gradient operators
  ! ------------------
  real(WP), dimension(:,:,:), pointer :: grad_x,grad_y,grad_z
  real(WP), dimension(:,:,:), pointer :: grad_xm,grad_ym,grad_zm
  real(WP), dimension(:,:,:), pointer :: secder_xm,secder_ym,secder_zm
  
contains

  ! Allocate all the arrays
  ! -----------------------
  subroutine metric_generic_allocate
    implicit none
    
    ! Inside only
    allocate(interp_u_xm(imin_:imax_,jmin_:jmax_,-st1:st2))
    allocate(interp_v_ym(imin_:imax_,jmin_:jmax_,-st1:st2))
    allocate(interp_w_zm(imin_:imax_,jmin_:jmax_,-st1:st2))
    
    ! Larger for strain rate calculation
    allocate(interp_uvw_x(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    allocate(interp_uvw_y(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    allocate(interp_uvw_z(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    
    ! Larger for diffusion term in scalar equation
    ! Larger for inflow/outflow
    allocate(interp_sc_x(imin_-st1:imax_+st2,jmino_   :jmaxo_   ,-st2:st1))
    allocate(interp_sc_y(imino_   :imaxo_   ,jmin_-st1:jmax_+st2,-st2:st1))
    allocate(interp_sc_z(imino_   :imaxo_   ,jmin_-st1:jmax_+st2,-st2:st1))
    
    ! Larger because of momentum fluxes
    allocate(interp_sc_xy(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
    allocate(interp_sc_yz(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
    allocate(interp_sc_xz(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
    
    ! Inside only
    allocate(div_u(imin_:imax_,jmin_:jmax_,-st1:st2))
    allocate(div_v(imin_:imax_,jmin_:jmax_,-st1:st2))
    allocate(div_w(imin_:imax_,jmin_:jmax_,-st1:st2))
    
    ! Larger for diffusion term in scalar equation
    allocate(grad_x(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    allocate(grad_y(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    allocate(grad_z(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-st2:st1))
    
    ! Only inside the domain
    allocate(grad_xm(imin_:imax_,jmin_:jmax_,-stp:stp))
    allocate(grad_ym(imin_:imax_,jmin_:jmax_,-stp:stp))
    allocate(grad_zm(imin_:imax_,jmin_:jmax_,-stp:stp))
    
    ! Only inside the domain
    allocate(secder_xm(imin_:imax_,jmin_:jmax_,-stp:stp))
    allocate(secder_ym(imin_:imax_,jmin_:jmax_,-stp:stp))
    allocate(secder_zm(imin_:imax_,jmin_:jmax_,-stp:stp))
    
    return
  end subroutine metric_generic_allocate
  
  
  ! Compute the interpolation coefficients
  ! --------------------------------------
  subroutine metric_generic_compute_coeffs
    use math
    implicit none
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:),   pointer :: B
    real(WP), dimension(:),   pointer :: alpha
    integer :: i,l,i1,i2,n
    
    allocate(coeff_deriv(generic_order))
    allocate(coeff_interp(generic_order))
    allocate(coeff_interp_2D(generic_order,generic_order))
    
    ! Allocate the matrix/vector
    n = generic_order/2
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
    
    ! Compute the coefficients for the 'full' order
    do i=1,n
       i1 = n+1-i
       i2 = n+i
       coeff_interp(i1) = 0.5_WP * alpha(i)
       coeff_interp(i2) = 0.5_WP * alpha(i)
       coeff_deriv (i1) = - alpha(i) / real(2*i-1,WP)
       coeff_deriv (i2) = + alpha(i) / real(2*i-1,WP)
    end do
    
    do i1=1,generic_order
       do i2=1,generic_order
          coeff_interp_2D(i1,i2) = coeff_interp(i1)*coeff_interp(i2)
       end do
    end do
    
    ! Clean up memory
    deallocate(A);nullify(A)
    deallocate(B);nullify(B)
    deallocate(alpha);nullify(alpha)
    
    return
  end subroutine metric_generic_compute_coeffs


  ! Compute the metric coeffs in physical space
  ! -------------------------------------------
  subroutine metric_generic_compute_primary
    use math
    implicit none
    
    integer  :: i,j,n,st
    real(WP), dimension(:,:), pointer :: A,B
    
    ! Allocate necessary arrays
    allocate(A(-st2:st1,generic_order))
    allocate(B(generic_order,-st2:st1))
    
    ! Interpolation of a scalar at the faces
    do i=imin_-st1,imax_+st2
       do st=-st2,st1
          do n=1,generic_order
             A(st,n) = (xm(i+st)-x(i))**(n-1)
          end do
       end do
       call inverse_matrix(A,B,generic_order)
       do j=jmino_,jmaxo_
          interp_sc_x(i,j,:) = B(1,:)
       end do
    end do
    do j=jmin_-st1,jmax_+st2
       do st=-st2,st1
          do n=1,generic_order
             A(st,n) = (ym(j+st)-y(j))**(n-1)
          end do
       end do
       call inverse_matrix(A,B,generic_order)
       do i=imino_,imaxo_
          interp_sc_y(i,j,:) = B(1,:)
       end do
    end do
    do j=jmin_-st1,jmax_+st2
       do i=imino_,imaxo_
          interp_sc_z(i,j,:) = coeff_interp
       end do
    end do
    
    ! Interpolation of the velocities at the center
    do j=jmin_,jmax_
       do i=imin_,imax_
          interp_u_xm(i,j,:) = coeff_interp
          interp_v_ym(i,j,:) = coeff_interp
          interp_w_zm(i,j,:) = coeff_interp
       end do
    end do
    
    ! Interpolation of the scalars at the corners
    do j=jmino_+st2,jmaxo_-st1
       do i=imino_+st2,imaxo_-st1
          interp_sc_xy(i,j,:,:) = coeff_interp_2D
          interp_sc_yz(i,j,:,:) = coeff_interp_2D
          interp_sc_xz(i,j,:,:) = coeff_interp_2D
       end do
    end do
    
    ! Interpolate centered velocities at the faces
    do j=jmin_-st1,jmax_+st2
       do i=imin_-st1,imax_+st2
          interp_uvw_x(i,j,:) = coeff_interp
          interp_uvw_y(i,j,:) = coeff_interp
          interp_uvw_z(i,j,:) = coeff_interp
       end do
    end do
    
    ! Divergence of a vector
    do j=jmin_,jmax_
       do i=imin_,imax_
          div_u(i,j,:) = coeff_deriv * dxi(i)
          div_v(i,j,:) = coeff_deriv * dyi(j)
          div_w(i,j,:) = coeff_deriv * dzi
       end do
    end do
    
    ! Gradient of a scalar
    do j=jmin_-st1,jmax_+st2
       do i=imin_-st1,imax_+st2
          grad_x(i,j,:) = coeff_deriv * dxmi(i-1)
          grad_y(i,j,:) = coeff_deriv * dymi(j-1)
          grad_z(i,j,:) = coeff_deriv * dzi
       end do
    end do
    
    return
  end subroutine metric_generic_compute_primary
  
  
  ! Compute the metric coeffs in physical space
  ! Force second order in y for cylindrical cases
  ! ---------------------------------------------
  subroutine metric_generic_enforce_cyl
    implicit none
    
    integer  :: i,j
    
    ! Do nothing more if not cylindrical
    if (icyl.ne.1) return
    
    ! Divergence of a vector
    do j=jmin_,jmax_
       do i=imin_,imax_
          div_v(i,j,:) = coeff_deriv * dyi(j) * y(j-st1:j+st2) * ymi(j)
          div_w(i,j,:) = coeff_deriv * dzi * ymi(j)
       end do
    end do
    
    ! Gradient of a scalar
    do j=jmin_-st1,jmax_+st2
       do i=imin_,imax_
          grad_y(i,j,:) = coeff_deriv * dymi(j-1)
       end do
    end do
    do j=jmin_,jmax_
       do i=imin_,imax_
          grad_z(i,j,:) = coeff_deriv * dzi * ymi(j)
       end do
    end do
    
    return
  end subroutine metric_generic_enforce_cyl
  
  
  ! Change the values of the metric to enforce wall conditions
  ! ----------------------------------------------------------
  subroutine metric_generic_enforce_walls
    use masks
    implicit none
    
    integer  :: i,j,st,n
    
    ! Interpolation of a scalar at the faces
    do j=jmino_,jmaxo_
       do i=imin_-st1,imax_+st2
          do st=-st2,-1
             if (mask(i+st,j).eq.1) then
                interp_sc_x(i,j,st+1) = interp_sc_x(i,j,st+1) + interp_sc_x(i,j,st)
                interp_sc_x(i,j,st)   = 0.0_WP
             end if
          end do
          do st=st1,0,-1
             if (mask(i+st,j).eq.1) then
                interp_sc_x(i,j,st-1) = interp_sc_x(i,j,st-1) + interp_sc_x(i,j,st)
                interp_sc_x(i,j,st)   = 0.0_WP
             end if
          end do
       end do
    end do
    do j=jmin_-st1,jmax_+st2
       do i=imino_,imaxo_
          do st=-st2,-1
             if (mask(i,j+st).eq.1) then
                interp_sc_y(i,j,st+1) = interp_sc_y(i,j,st+1) + interp_sc_y(i,j,st)
                interp_sc_y(i,j,st)   = 0.0_WP
             end if
          end do
          do st=st1,0,-1
             if (mask(i,j+st).eq.1) then
                interp_sc_y(i,j,st-1) = interp_sc_y(i,j,st-1) + interp_sc_y(i,j,st)
                interp_sc_y(i,j,st)   = 0.0_WP
             end if
          end do
       end do
    end do
    
    ! Interpolation of a scalar at the corners
    do j=jmino_+st2,jmaxo_-st1
       do i=imino_+st2,imaxo_-st1
          do st=-st2,st1
             do n=-st2,st1
                if (mask(i+st,j+n).eq.1) interp_sc_xy(i,j,st,n) = 0.0_WP
             end do
             if (mask(i+st,j).eq.1) interp_sc_xz(i,j,st,:) = 0.0_WP
             if (mask(i,j+st).eq.1) interp_sc_yz(i,j,st,:) = 0.0_WP
          end do
          interp_sc_xy(i,j,:,:) = interp_sc_xy(i,j,:,:) / (sum(interp_sc_xy(i,j,:,:))+epsilon(1.0_WP))
          interp_sc_xz(i,j,:,:) = interp_sc_xz(i,j,:,:) / (sum(interp_sc_xz(i,j,:,:))+epsilon(1.0_WP))
          interp_sc_yz(i,j,:,:) = interp_sc_yz(i,j,:,:) / (sum(interp_sc_yz(i,j,:,:))+epsilon(1.0_WP))
       end do
    end do
    
    ! Interpolation of the velocities at the center
    do j=jmin_,jmax_
       do i=imin_,imax_
          do st=-st1,st2
             if (mask_u(i+st,j).eq.1) interp_u_xm(i,j,st) = 0.0_WP
             if (mask_v(i,j+st).eq.1) interp_v_ym(i,j,st) = 0.0_WP
          end do
          if (mask(i,j).eq.1) then
             interp_u_xm(i,j,:) = 0.0_WP
             interp_v_ym(i,j,:) = 0.0_WP
             interp_w_zm(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Interpolate centered velocities at the faces
    do j=jmin_-st1,jmax_+st2
       do i=imin_-st1,imax_+st2
          do st=-st2,st1
             if (mask(i+st,j).eq.1) interp_uvw_x(i,j,st) = 0.0_WP
             if (mask(i,j+st).eq.1) interp_uvw_y(i,j,st) = 0.0_WP
          end do
          if (mask_u(i,j).eq.1) interp_uvw_x(i,j,:) = 0.0_WP
          if (mask_v(i,j).eq.1) interp_uvw_y(i,j,:) = 0.0_WP
          if (mask_w(i,j).eq.1) interp_uvw_z(i,j,:) = 0.0_WP
       end do
    end do
    
    ! Divergence of a vector
    do j=jmin_,jmax_
       do i=imin_,imax_
          do st=-st1,st2
             if (mask_u(i+st,j).eq.1) div_u(i,j,st) = 0.0_WP
             if (mask_v(i,j+st).eq.1) div_v(i,j,st) = 0.0_WP
          end do
          if (mask(i,j).eq.1) then
             div_u(i,j,:) = 0.0_WP
             div_v(i,j,:) = 0.0_WP
             div_w(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Gradient of a scalar
    do j=jmin_-st1,jmax_+st2
       do i=imin_-st1,imax_+st2
          do st=-st2,-1
             if (mask(i+st,j).eq.1) then
                grad_x(i,j,st+1) = grad_x(i,j,st+1) + grad_x(i,j,st)
                grad_x(i,j,st)   = 0.0_WP
             end if
             if (mask(i,j+st).eq.1) then
                grad_y(i,j,st+1) = grad_y(i,j,st+1) + grad_y(i,j,st)
                grad_y(i,j,st)   = 0.0_WP
             end if
          end do
          do st=st1,0,-1
             if (mask(i+st,j).eq.1) then
                grad_x(i,j,st-1) = grad_x(i,j,st-1) + grad_x(i,j,st)
                grad_x(i,j,st)   = 0.0_WP
             end if
             if (mask(i,j+st).eq.1) then
                grad_y(i,j,st-1) = grad_y(i,j,st-1) + grad_y(i,j,st)
                grad_y(i,j,st)   = 0.0_WP
             end if
          end do
          if (mask_u(i,j).eq.1) grad_x(i,j,:) = 0.0_WP
          if (mask_v(i,j).eq.1) grad_y(i,j,:) = 0.0_WP
          if (mask_w(i,j).eq.1) grad_z(i,j,:) = 0.0_WP
       end do
    end do
    
    return
  end subroutine metric_generic_enforce_walls
  
  
  ! Enforce the physical boundary conditions of the domain
  ! ------------------------------------------------------
  subroutine metric_generic_enforce_bc
    use parallel
    implicit none
    integer :: i,j,st
    
    if (xper.ne.1) then
       ! Left boundary
       ! -> Force grad_x to zero for sgs
       ! -> Nothing for interp_sc as physical values in ghost cells
       if (iproc.eq.1) then
          
          do i=imin-st1,imax_+st2
             do st=-st2,-1
                if (i+st.lt.imin) then
                   grad_x(i,:,st+1) = grad_x(i,:,st+1) + grad_x(i,:,st)
                   grad_x(i,:,st)   = 0.0_WP
                end if
             end do
             if (i.le.imin) grad_x(i,:,:) = 0.0_WP
          end do
          
          do i=imin,imax_
             do st=-st1,0
                if (i+st.lt.imin) then
                   div_u(i,:,st+1) = div_u(i,:,st+1) + div_u(i,:,st)
                   div_u(i,:,st)   = 0.0_WP
                end if
             end do
             if (i.lt.imin) then
                div_u(i,:,:) = 0.0_WP
                div_v(i,:,:) = 0.0_WP
                div_w(i,:,:) = 0.0_WP
             end if
          end do
       end if
       ! Right boundary
       ! -> Force grad_x to zero for sgs
       ! -> Nothing for interp_sc as physical values in ghost cells
       if (iproc.eq.npx) then
          do i=imin_-st1,imax+st2
             do st=st1,0,-1
                if (i+st.ge.imax+1) then
                   grad_x(i,:,st-1) = grad_x(i,:,st-1) + grad_x(i,:,st)
                   grad_x(i,:,st)   = 0.0_WP
                end if
             end do
             if (i.ge.imax+1) grad_x(i,:,:) = 0.0_WP
          end do

          do i=imin_,imax_
             do st=st2,1,-1
                if (i+st.gt.imax+1) then
                   div_u(i,:,st-1) = div_u(i,:,st-1) + div_u(i,:,st)
                   div_u(i,:,st)   = 0.0_WP
                end if
             end do
             if (i.gt.imax) then
                div_u(i,:,:) = 0.0_WP
                div_v(i,:,:) = 0.0_WP
                div_w(i,:,:) = 0.0_WP
             end if
          end do
       end if
    end if
    
    if (yper.ne.1) then
       ! Lower boundary
       ! -> Newmann on scalars
       ! -> Dirichlet on V
       if (jproc.eq.1 .and. icyl.eq.0) then
          
          do j=jmin-st1,jmax_+st2
             do st=-st2,-1
                if (j+st.lt.jmin) then
                   grad_y(:,j,st+1) = grad_y(:,j,st+1) + grad_y(:,j,st)
                   grad_y(:,j,st)   = 0.0_WP
                   interp_sc_y(:,j,st+1) = interp_sc_y(:,j,st+1) + interp_sc_y(:,j,st)
                   interp_sc_y(:,j,st)   = 0.0_WP
                end if
             end do
             if (j.le.jmin) grad_y(:,j,:) = 0.0_WP
          end do
          
          do j=jmin_,jmax_
             do st=-st1,0
                if (j+st.le.jmin) then
!!$                   div_v(:,j,st+1) = div_v(:,j,st+1) + div_v(:,j,st) ! Neumann
                   div_v(:,j,st) = 0.0_WP
                end if
             end do
             if (j.lt.jmin) then
                div_u(:,j,:) = 0.0_WP
                div_v(:,j,:) = 0.0_WP
                div_w(:,j,:) = 0.0_WP
             end if
          end do
       end if
       ! Upper boundary
       ! -> Newmann on scalars
       ! -> Dirichlet on V
       if (jproc.eq.npy) then
          do j=jmin_-st1,jmax+st2
             do st=st1,0,-1
                if (j+st.gt.jmax+1) then
                   grad_y(:,j,st-1) = grad_y(:,j,st-1) + grad_y(:,j,st)
                   grad_y(:,j,st)   = 0.0_WP
                   interp_sc_y(:,j,st-1) = interp_sc_y(:,j,st-1) + interp_sc_y(:,j,st)
                   interp_sc_y(:,j,st)   = 0.0_WP
                end if
             end do
             if (j.ge.jmax+1) grad_y(:,j,:) = 0.0_WP
          end do
          
          do j=jmin_,jmax_
             do st=st2,1,-1
                if (j+st.ge.jmax+1) then
!!$                   div_v(:,j,st-1) = div_v(:,j,st-1) + div_v(:,j,st) ! Neumann
                   div_v(:,j,st) = 0.0_WP
                end if
             end do
             if (j.gt.jmax) then
                div_u(:,j,:) = 0.0_WP
                div_v(:,j,:) = 0.0_WP
                div_w(:,j,:) = 0.0_WP
             end if
          end do
       end if
    end if
    
    return
  end subroutine metric_generic_enforce_bc
  
  
  ! Enforce the dimension of the problem
  ! ------------------------------------
  subroutine metric_generic_enforce_dim
    implicit none
    
    ! Remove dependance in x and U
    if (nx.eq.1) then
       interp_u_xm = 0.0_WP
       div_u = 0.0_WP
       grad_x = 0.0_WP
    end if
    
    ! Remove dependance in y and V
    if (ny.eq.1) then
       interp_v_ym = 0.0_WP
       div_v = 0.0_WP
       grad_y = 0.0_WP
    end if
    
    ! Remove dependance in z and W
    if (nz.eq.1) then
       interp_w_zm = 0.0_WP
       div_w = 0.0_WP
       grad_z = 0.0_WP
    end if
    
    return
  end subroutine metric_generic_enforce_dim
  
  
  ! Enforce Dirichlet condition
  ! ---------------------------
  subroutine metric_generic_dirichlet
    use masks
    implicit none
    
    integer :: i,j
    
    do j=jmin_,jmax_
       do i=imin_,imax_
          ! Force zero RHS
          if (mask(i,j).eq.3) then
             div_u(i,j,:) = 0.0_WP
             div_v(i,j,:) = 0.0_WP
             div_w(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    do j=jmin_-st1,jmax_+st2
       do i=imin_-st1,imax_+st2
          ! Force correct inflow flux of density
          if (mask_u(i,j).eq.3) then
             interp_sc_x(i,j,:)  = 0.0_WP
             interp_sc_x(i,j,-1) = 1.0_WP
          end if
          ! No scalar flux
          if (mask_u(i,j).eq.3) grad_x(i,j,:) = 0.0_WP
          if (mask_v(i,j).eq.3) grad_y(i,j,:) = 0.0_WP
          if (mask_w(i,j).eq.3) grad_z(i,j,:) = 0.0_WP
       end do
    end do
    
    return
  end subroutine metric_generic_dirichlet
  
  
  ! Compute the secondary metric dependant on primary metric
  ! --------------------------------------------------------
  subroutine metric_generic_compute_secondary
    implicit none
    integer :: i,j
    integer :: s1,s2
    
    ! Centered Gradient operator
    do j=jmin_,jmax_
       do i=imin_,imax_
          ! Initialize
          grad_xm(i,j,:) = 0.0_WP
          grad_ym(i,j,:) = 0.0_WP
          grad_zm(i,j,:) = 0.0_WP
          ! Convolute
          do s1=-st1,st2
             do s2=-st2,st1
                grad_xm(i,j,s1+s2) = grad_xm(i,j,s1+s2) + &
                     interp_u_xm(i,j,s1)*grad_x(i+s1,j,s2)
                grad_ym(i,j,s1+s2) = grad_ym(i,j,s1+s2) + &
                     interp_v_ym(i,j,s1)*grad_y(i,j+s1,s2)
                grad_zm(i,j,s1+s2) = grad_zm(i,j,s1+s2) + &
                     interp_w_zm(i,j,s1)*grad_z(i,j,s2)
             end do
          end do
       end do
    end do
    
    ! Centered second derivative operators
    do j=jmin_,jmax_
       do i=imin_,imax_
          ! Initialize
          secder_xm(i,j,:) = 0.0_WP
          secder_ym(i,j,:) = 0.0_WP
          secder_zm(i,j,:) = 0.0_WP
          ! Convolute
          do s1=-st1,st2
             do s2=-st2,st1
                secder_xm(i,j,s1+s2) = secder_xm(i,j,s1+s2) + &
                     div_u(i,j,s1)*grad_x(i+s1,j,s2)
                secder_ym(i,j,s1+s2) = secder_ym(i,j,s1+s2) + &
                     div_v(i,j,s1)*grad_y(i,j+s1,s2)
                secder_zm(i,j,s1+s2) = secder_zm(i,j,s1+s2) + &
                     div_w(i,j,s1)*grad_z(i,j,s2)
             end do
          end do
       end do
    end do
    
    return
  end subroutine metric_generic_compute_secondary
  
end module metric_generic


subroutine metric_generic_init
  use metric_generic
  use masks
  implicit none
  
  integer :: i,j
  
  ! Compute the cells volume
  do j=jmino,jmaxo
     do i=imino,imaxo
        vol(i,j) = dx(i)*dy(j)*dz
        if (icyl.eq.1) vol(i,j) = vol(i,j) * abs(ym(j))
        if (mask(i,j).eq.1) vol(i,j) = 0.0_WP
     end do
  end do
  vol_total = sum(vol(imin:imax,jmin:jmax))*nz
  
  ! Set the stencil lengths
  st2 = generic_order/2
  st1 = st2 - 1
  stp = st1 + st2
  
  ! Allocate all the arrays
  call metric_generic_allocate
  
  ! Compute the interpolation weights
  call metric_generic_compute_coeffs
  
  ! Compute the metric
  call metric_generic_compute_primary
  
  ! Enforce some properties
  call metric_generic_enforce_cyl
  call metric_generic_enforce_dim
  call metric_generic_enforce_bc
  call metric_generic_enforce_walls
  call metric_generic_dirichlet
  
  ! Compute the secondary metrics
  call metric_generic_compute_secondary
  
  return
end subroutine metric_generic_init
