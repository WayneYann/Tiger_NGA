module metric_velocity_visc
  use precision
  use geometry
  use partition
  implicit none
  
  ! Stencil lengths
  ! ---------------
  integer :: stv1,stv2
  integer :: stvp
  
  ! Interpolation operators
  ! -----------------------
  ! Additional cylindrical terms
  real(WP), dimension(:,:,:), pointer :: interpv_cyl_F_ym,interpv_cyl_F_y
  real(WP), dimension(:,:,:), pointer :: interpv_cyl_v_ym,interpv_cyl_w_y
  
  ! Divergence operators
  ! --------------------
  ! Trace of velocity gradient tensor
  real(WP), dimension(:,:,:), pointer :: divv_u,divv_v,divv_w
  ! Divergence of the viscous fluxes
  real(WP), dimension(:,:,:), pointer :: divv_xx,divv_xy,divv_xz
  real(WP), dimension(:,:,:), pointer :: divv_yx,divv_yy,divv_yz
  real(WP), dimension(:,:,:), pointer :: divv_zx,divv_zy,divv_zz
  
  ! Gradient operators
  ! ------------------
  ! Velocity gradients
  real(WP), dimension(:,:,:), pointer :: grad_u_x,grad_u_y,grad_u_z
  real(WP), dimension(:,:,:), pointer :: grad_v_x,grad_v_y,grad_v_z
  real(WP), dimension(:,:,:), pointer :: grad_w_x,grad_w_y,grad_w_z
  
contains
  
  ! Allocate all the arrays
  ! -----------------------
  subroutine metric_velocity_visc_allocate
    implicit none
    
    ! Interior only for fluxes
    allocate(interpv_cyl_F_ym(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(interpv_cyl_F_y (imin_:imax_,jmin_:jmax_,-stv2:stv1))
    
    ! Larger because of momentum fluxes
    allocate(interpv_cyl_v_ym(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    allocate(interpv_cyl_w_y (imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    
    ! Larger because of momentum fluxes
    allocate(divv_u(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    allocate(divv_v(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    allocate(divv_w(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    
    ! Only inside the domain
    allocate(divv_xx(imin_:imax_,jmin_:jmax_,-stv2:stv1))
    allocate(divv_xy(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_xz(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_yx(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_yy(imin_:imax_,jmin_:jmax_,-stv2:stv1))
    allocate(divv_yz(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_zx(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_zy(imin_:imax_,jmin_:jmax_,-stv1:stv2))
    allocate(divv_zz(imin_:imax_,jmin_:jmax_,-stv2:stv1))
    
    ! Larger because of momentum fluxes
    allocate(grad_u_x(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    allocate(grad_u_y(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    allocate(grad_u_z(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    
    ! Larger because of momentum fluxes
    allocate(grad_v_x(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    allocate(grad_v_y(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    allocate(grad_v_z(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    
    ! Larger because of momentum fluxes
    allocate(grad_w_x(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    allocate(grad_w_y(imin_-stv1:imax_+stv2,jmin_-stv1:jmax_+stv2,-stv2:stv1))
    allocate(grad_w_z(imin_-stv2:imax_+stv1,jmin_-stv2:jmax_+stv1,-stv1:stv2))
    
    return
  end subroutine metric_velocity_visc_allocate
  
  
  ! Compute the metric coeffs in physical space
  ! -------------------------------------------
  subroutine metric_velocity_visc_compute_primary
    use math
    implicit none
    
    integer  :: i,j,k
    
    ! Use first plane
    k = kmin_
    
    ! Divergence of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          call hofdd(vel_visc_order,x(i-stv1:i+stv2),xm(i),divv_u(i,j,:))
          call hofdd(vel_visc_order,y(j-stv1:j+stv2),ym(j),divv_v(i,j,:))
          call hofdd(vel_visc_order,z(k-stv1:k+stv2),zm(k),divv_w(i,j,:))
       end do
    end do
    
    ! Divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          call hofdd(vel_visc_order,xm(i-stv2:i+stv1),x (i),divv_xx(i,j,:))
          call hofdd(vel_visc_order,y (j-stv1:j+stv2),ym(j),divv_xy(i,j,:))
          call hofdd(vel_visc_order,z (k-stv1:k+stv2),zm(k),divv_xz(i,j,:))
          
          call hofdd(vel_visc_order,x (i-stv1:i+stv2),xm(i),divv_yx(i,j,:))
          call hofdd(vel_visc_order,ym(j-stv2:j+stv1),y (j),divv_yy(i,j,:))
          call hofdd(vel_visc_order,z (k-stv1:k+stv2),zm(k),divv_yz(i,j,:))
          
          call hofdd(vel_visc_order,x (i-stv1:i+stv2),xm(i),divv_zx(i,j,:))
          call hofdd(vel_visc_order,y (j-stv1:j+stv2),ym(j),divv_zy(i,j,:))
          call hofdd(vel_visc_order,zm(k-stv2:k+stv1),z (k),divv_zz(i,j,:))
       end do
    end do
    
    ! Gradient of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          call hofdd(vel_visc_order,x(i-stv1:i+stv2),xm(i),grad_u_x(i,j,:))
          call hofdd(vel_visc_order,y(j-stv1:j+stv2),ym(j),grad_v_y(i,j,:))
          call hofdd(vel_visc_order,z(k-stv1:k+stv2),zm(k),grad_w_z(i,j,:))
       end do
    end do
    do j=jmin_-stv1,jmax_+stv2
       do i=imin_-stv1,imax_+stv2
          call hofdd(vel_visc_order,xm(i-stv2:i+stv1),x(i),grad_v_x(i,j,:))
          call hofdd(vel_visc_order,xm(i-stv2:i+stv1),x(i),grad_w_x(i,j,:))
          call hofdd(vel_visc_order,ym(j-stv2:j+stv1),y(j),grad_u_y(i,j,:))
          call hofdd(vel_visc_order,ym(j-stv2:j+stv1),y(j),grad_w_y(i,j,:))
          call hofdd(vel_visc_order,zm(k-stv2:k+stv1),z(k),grad_u_z(i,j,:))
          call hofdd(vel_visc_order,zm(k-stv2:k+stv1),z(k),grad_v_z(i,j,:))
       end do
    end do
    
    ! Interpolation of the velocities at the faces
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          call hofdi(vel_visc_order,y (j-stv1:j+stv2),ym(j),interpv_cyl_v_ym(i,j,:))
       end do
    end do
    
    ! Interpolation of the velocities at the corners
    do j=jmin_-stv1,jmax_+stv2
       do i=imin_-stv1,imax_+stv2
          call hofdi(vel_visc_order,ym(j-stv2:j+stv1),y (j),interpv_cyl_w_y (i,j,:))
       end do
    end do
    
    ! Interpolation of fluxes at the faces
    do j=jmin_,jmax_
       do i=imin_,imax_
          call hofdi(vel_visc_order,y (j-stv1:j+stv2),ym(j),interpv_cyl_F_ym(i,j,:))
          call hofdi(vel_visc_order,ym(j-stv2:j+stv1),y (j),interpv_cyl_F_y (i,j,:))
       end do
    end do
    
    return
  end subroutine metric_velocity_visc_compute_primary
  
  
  ! Compute the metric coeffs in physical space
  ! Force second order in y for cylindrical cases
  ! ---------------------------------------------
  subroutine metric_velocity_visc_enforce_cyl
    use math
    use parallel
    implicit none
    
    integer  :: i,j
    
    ! If not cylindrical set all interpolations to zero and exit
    if (icyl.ne.1) then
       interpv_cyl_v_ym = 0.0_WP
       interpv_cyl_w_y  = 0.0_WP
       interpv_cyl_F_ym = 0.0_WP
       interpv_cyl_F_y  = 0.0_WP
       return
    end if
    
    ! Divergence of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          divv_v(i,j,:) = divv_v(i,j,:) * y(j-stv1:j+stv2) * ymi(j)
          divv_w(i,j,:) = divv_w(i,j,:) * ymi(j)
       end do
    end do
    
    ! Divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          divv_xy(i,j,:) = divv_xy(i,j,:) * y (j-stv1:j+stv2) * ymi(j)
          divv_yy(i,j,:) = divv_yy(i,j,:) * ym(j-stv2:j+stv1) * yi (j)
          divv_zy(i,j,:) = divv_zy(i,j,:) * y (j-stv1:j+stv2) * ymi(j)
          
          divv_xz(i,j,:) = divv_xz(i,j,:) * ymi(j)
          divv_yz(i,j,:) = divv_yz(i,j,:) * yi (j)
          divv_zz(i,j,:) = divv_zz(i,j,:) * ymi(j)
       end do
    end do
    
    ! Gradient of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          grad_w_z(i,j,:) = grad_w_z(i,j,:) * ymi(j)
       end do
    end do
    do j=jmin_-stv1,jmax_+stv2
       do i=imin_-stv1,imax_+stv2
          grad_u_z(i,j,:) = grad_u_z(i,j,:) * ymi(j)
          grad_v_z(i,j,:) = grad_v_z(i,j,:) * yi (j)
       end do
    end do
    
    ! Force exact (zero) flux at the axis
    if (jproc.eq.1) grad_w_y(:,jmin_,:) = 0.0_WP
    
    return
  end subroutine metric_velocity_visc_enforce_cyl
  
  
  ! Change the values of the metric to enforce wall conditions
  ! ----------------------------------------------------------
  subroutine metric_velocity_visc_enforce_walls
    use masks
    use math
    use parallel
    implicit none
    
    integer  :: i,j,k,st
    
    ! Use first plane
    k = kmin_
    
    ! Divergence of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          do st=-stv1,0
             if (mask_u(i+st,j).eq.1) then
                divv_u(i,j,:) = 0.0_WP
                call hofdd_d(stv2-st,x(i+st+1:i+stv2),xm(i),divv_u(i,j,st+1:stv2),x(i+st))
             end if
             if (mask_v(i,j+st).eq.1) then
                divv_v(i,j,:) = 0.0_WP
                call hofdd_d(stv2-st,y(j+st+1:j+stv2),ym(j),divv_v(i,j,st+1:stv2),y(j+st))
             end if
          end do
          do st=stv2,1,-1
             if (mask_u(i+st,j).eq.1) then
                divv_u(i,j,:) = 0.0_WP
                call hofdd_d(st+stv1,x(i-stv1:i+st-1),xm(i),divv_u(i,j,-stv1:st-1),x(i+st))
             end if
             if (mask_v(i,j+st).eq.1) then
                divv_v(i,j,:) = 0.0_WP
                call hofdd_d(st+stv1,y(j-stv1:j+st-1),ym(j),divv_v(i,j,-stv1:st-1),y(j+st))
             end if
          end do
          if (mask(i,j).eq.1) then
             divv_u(i,j,:) = 0.0_WP
             divv_v(i,j,:) = 0.0_WP
             divv_w(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Divergence of a matrix
    do j=jmin_,jmax_
       do i=imin_,imax_
          
          if (mask_u(i,j).eq.1) then
             divv_xx(i,j,:) = 0.0_WP
             divv_xy(i,j,:) = 0.0_WP
             divv_xz(i,j,:) = 0.0_WP
          else
             do st=-stv2,-1
                if (mask(i+st,j).eq.1) then
                   divv_xx(i,j,:) = 0.0_WP
                   call hofdd(stv1-st,xm(i+st+1:i+stv1),x(i),divv_xx(i,j,st+1:stv1))
                end if
             end do
             do st=stv1,0,-1
                if (mask(i+st,j).eq.1) then
                   divv_xx(i,j,:) = 0.0_WP
                   call hofdd(st+stv2,xm(i-stv2:i+st-1),x(i),divv_xx(i,j,-stv2:st-1))
                end if
             end do
             do st=-stv1,0
                if (mask_u(i,j+st-1).eq.1) then
                   divv_xy(i,j,:) = 0.0_WP
                   call hofdd(stv2-st+1,y(j+st:j+stv2),ym(j),divv_xy(i,j,st:stv2))
                end if
             end do
             do st=stv2,1,-1
                if (mask_u(i,j+st).eq.1) then
                   divv_xy(i,j,:) = 0.0_WP
                   call hofdd(st+stv1+1,y(j-stv1:j+st),ym(j),divv_xy(i,j,-stv1:st))
                end if
             end do
          end if
          
          if (mask_v(i,j).eq.1) then
             divv_yx(i,j,:) = 0.0_WP
             divv_yy(i,j,:) = 0.0_WP
             divv_yz(i,j,:) = 0.0_WP
          else
             do st=-stv2,-1
                if (mask(i,j+st).eq.1) then
                   divv_yy(i,j,:) = 0.0_WP
                   call hofdd(stv1-st,ym(j+st+1:j+stv1),y(j),divv_yy(i,j,st+1:stv1))
                end if
             end do
             do st=stv1,0,-1
                if (mask(i,j+st).eq.1) then
                   divv_yy(i,j,:) = 0.0_WP
                   call hofdd(st+stv2,ym(j-stv2:j+st-1),y(j),divv_yy(i,j,-stv2:st-1))
                end if
             end do
             do st=-stv1,0
                if (mask_v(i+st-1,j).eq.1) then
                   divv_yx(i,j,:) = 0.0_WP
                   call hofdd(stv2-st+1,x(i+st:i+stv2),xm(i),divv_yx(i,j,st:stv2))
                end if
             end do
             do st=stv2,1,-1
                if (mask_v(i+st,j).eq.1) then
                   divv_yx(i,j,:) = 0.0_WP
                   call hofdd(st+stv1+1,x(i-stv1:i+st),xm(i),divv_yx(i,j,-stv1:st))
                end if
             end do
          end if
          
          if (mask_w(i,j).eq.1) then
             divv_zx(i,j,:) = 0.0_WP
             divv_zy(i,j,:) = 0.0_WP
             divv_zz(i,j,:) = 0.0_WP
          else
             do st=-stv1,0
                if (mask_w(i+st-1,j).eq.1) then
                   divv_zx(i,j,:) = 0.0_WP
                   call hofdd(stv2-st+1,x(i+st:i+stv2),xm(i),divv_zx(i,j,st:stv2))
                end if
                if (mask_w(i,j+st-1).eq.1) then
                   divv_zy(i,j,:) = 0.0_WP
                   call hofdd(stv2-st+1,y(j+st:j+stv2),ym(j),divv_zy(i,j,st:stv2))
                end if
             end do
             do st=stv2,1,-1
                if (mask_w(i+st,j).eq.1) then
                   divv_zx(i,j,:) = 0.0_WP
                   call hofdd(st+stv1+1,x(i-stv1:i+st),xm(i),divv_zx(i,j,-stv1:st))
                end if
                if (mask_w(i,j+st).eq.1) then
                   divv_zy(i,j,:) = 0.0_WP
                   call hofdd(st+stv1+1,y(j-stv1:j+st),ym(j),divv_zy(i,j,-stv1:st))
                end if
             end do
          end if

       end do
    end do
    
    ! Gradient of a vector
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          do st=-stv1,0
             if (mask_u(i+st,j).eq.1) then
                grad_u_x(i,j,:) = 0.0_WP
                call hofdd_d(stv2-st,x(i+st+1:i+stv2),xm(i),grad_u_x(i,j,st+1:stv2),x(i+st))
             end if
             if (mask_v(i,j+st).eq.1) then
                grad_v_y(i,j,:) = 0.0_WP
                call hofdd_d(stv2-st,y(j+st+1:j+stv2),ym(j),grad_v_y(i,j,st+1:stv2),y(j+st))
             end if
          end do
          do st=stv2,1,-1
             if (mask_u(i+st,j).eq.1) then
                grad_u_x(i,j,:) = 0.0_WP
                call hofdd_d(st+stv1,x(i-stv1:i+st-1),xm(i),grad_u_x(i,j,-stv1:st-1),x(i+st))
             end if
             if (mask_v(i,j+st).eq.1) then
                grad_v_y(i,j,:) = 0.0_WP
                call hofdd_d(st+stv1,y(j-stv1:j+st-1),ym(j),grad_v_y(i,j,-stv1:st-1),y(j+st))
             end if
          end do
          if (mask(i,j).eq.1) then
             grad_u_x(i,j,:) = 0.0_WP
             grad_v_y(i,j,:) = 0.0_WP
             grad_w_z(i,j,:) = 0.0_WP
          end if
       end do
    end do
    do j=jmin_-stv1,jmax_+stv2
       do i=imin_-stv1,imax_+stv2
          do st=-stv2,-1
             if (mask_v(i+st,j).eq.1) then
                grad_v_x(i,j,:) = 0.0_WP
                call hofdd_d(stv1-st,xm(i+st+1:i+stv1),x(i),grad_v_x(i,j,st+1:stv1),x(i+st+1))
             end if
             if (mask_w(i+st,j).eq.1) then
                grad_w_x(i,j,:) = 0.0_WP
                call hofdd_d(stv1-st,xm(i+st+1:i+stv1),x(i),grad_w_x(i,j,st+1:stv1),x(i+st+1))
             end if
             if (mask_u(i,j+st).eq.1) then
                grad_u_y(i,j,:) = 0.0_WP
                call hofdd_d(stv1-st,ym(j+st+1:j+stv1),y(j),grad_u_y(i,j,st+1:stv1),y(j+st+1))
             end if
             if (mask_w(i,j+st).eq.1) then
                grad_w_y(i,j,:) = 0.0_WP
                call hofdd_d(stv1-st,ym(j+st+1:j+stv1),y(j),grad_w_y(i,j,st+1:stv1),y(j+st+1))
             end if
          end do
          do st=stv1,0,-1
             if (mask_v(i+st,j).eq.1) then
                grad_v_x(i,j,:) = 0.0_WP
                call hofdd_d(st+stv2,xm(i-stv2:i+st-1),x(i),grad_v_x(i,j,-stv2:st-1),x(i+st))
             end if
             if (mask_w(i+st,j).eq.1) then
                grad_w_x(i,j,:) = 0.0_WP
                call hofdd_d(st+stv2,xm(i-stv2:i+st-1),x(i),grad_w_x(i,j,-stv2:st-1),x(i+st))
             end if
             if (mask_u(i,j+st).eq.1) then
                grad_u_y(i,j,:) = 0.0_WP
                call hofdd_d(st+stv2,ym(j-stv2:j+st-1),y(j),grad_u_y(i,j,-stv2:st-1),y(j+st))
             end if
             if (mask_w(i,j+st).eq.1) then
                grad_w_y(i,j,:) = 0.0_WP
                call hofdd_d(st+stv2,ym(j-stv2:j+st-1),y(j),grad_w_y(i,j,-stv2:st-1),y(j+st))
             end if
          end do
          if (mask(i,j).eq.1) then
             grad_u_z(i,j,:) = 0.0_WP
             grad_v_z(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    ! Interpolation of the velocities at the faces
    do j=jmin_-stv2,jmax_+stv1
       do i=imin_-stv2,imax_+stv1
          do st=-stv1,0
             if (mask_v(i,j+st).eq.1) then
                interpv_cyl_v_ym(i,j,:) = 0.0_WP
                call hofdi_d(stv2-st,y(j+st+1:j+stv2),ym(j),interpv_cyl_v_ym(i,j,st+1:stv2),y(j+st))
             end if
          end do
          do st=stv2,1,-1
             if (mask_v(i,j+st).eq.1) then
                interpv_cyl_v_ym(i,j,:) = 0.0_WP
                call hofdi_d(st+stv1,y(j-stv1:j+st-1),ym(j),interpv_cyl_v_ym(i,j,-stv1:st-1),y(j+st))
             end if
          end do
          if (mask(i,j).eq.1) interpv_cyl_v_ym(i,j,:) = 0.0_WP
       end do
    end do
    
    ! Interpolation of the velocities at the corners
    do j=jmin_-stv1,jmax_+stv2
       do i=imin_-stv1,imax_+stv2
          do st=-stv2,-1
             if (mask_w(i,j+st).eq.1) then
                interpv_cyl_w_y(i,j,:) = 0.0_WP
                call hofdi_d(stv1-st,ym(j+st+1:j+stv1),y(j),interpv_cyl_w_y(i,j,st+1:stv1),y(j+st+1))
             end if
          end do
          do st=stv1,0,-1
             if (mask_w(i,j+st).eq.1) then
                interpv_cyl_w_y(i,j,:) = 0.0_WP
                call hofdi_d(st+stv2,ym(j-stv2:j+st-1),y(j),interpv_cyl_w_y(i,j,-stv2:st-1),y(j+st))
             end if
          end do
       end do
    end do
    
    ! Interpolation of fluxes at the faces
    do j=jmin_,jmax_
       do i=imin_,imax_
          if (mask_v(i,j).eq.1) then
             interpv_cyl_F_y(i,j,:) = 0.0_WP
          else
             do st=-stv2,-1
                if (mask(i,j+st).eq.1) then
                   interpv_cyl_F_y(i,j,:) = 0.0_WP
                   call hofdi_d(stv1-st,ym(j+st+1:j+stv1),y(j),interpv_cyl_F_y(i,j,st+1:stv1),y(j+st+1))
                end if
             end do
             do st=stv1,0,-1
                if (mask(i,j+st).eq.1) then
                   interpv_cyl_F_y(i,j,:) = 0.0_WP
                   call hofdi_d(st+stv2,ym(j-stv2:j+st-1),y(j),interpv_cyl_F_y(i,j,-stv2:st-1),y(j+st))
                end if
             end do
          end if
          if (mask_w(i,j).eq.1) then
             interpv_cyl_F_ym(i,j,:) = 0.0_WP
          else
             do st=-stv1,0
                if (mask_w(i,j+st-1).eq.1) then
                   interpv_cyl_F_ym(i,j,:) = 0.0_WP
                   call hofdi_d(stv2-st,y(j+st+1:j+stv2),ym(j),interpv_cyl_F_ym(i,j,st+1:stv2),y(j+st))
                end if
             end do
             do st=stv2,1,-1
                if (mask_w(i,j+st).eq.1) then
                   interpv_cyl_F_ym(i,j,:) = 0.0_WP
                   call hofdi_d(st+stv1,y(j-stv1:j+st-1),ym(j),interpv_cyl_F_ym(i,j,-stv1:st-1),y(j+st))
                end if
             end do
          end if
       end do
    end do
    
    if (xper.ne.1) then
       ! Left boundary
       ! -> Newmann on P
       ! -> Dirichlet on U,V,W
       if (iproc.eq.1) then
          divv_xx(imin,:,:) = 0.0_WP
          divv_xy(imin,:,:) = 0.0_WP
          divv_xz(imin,:,:) = 0.0_WP
       end if
    end if
    
    if (yper.ne.1) then
       ! Lower boundary - Cartesian
       ! -> Newmann on P
       ! -> Newmann on U/W
       ! -> Dirichlet on V
       if (jproc.eq.1 .and. icyl.eq.0) then
          divv_yx(:,jmin,:) = 0.0_WP
          divv_yy(:,jmin,:) = 0.0_WP
          divv_yz(:,jmin,:) = 0.0_WP
       end if
    end if

    return
  end subroutine metric_velocity_visc_enforce_walls
  
  
  ! Enforce the physical boundary conditions of the domain
  ! ------------------------------------------------------
  subroutine metric_velocity_visc_enforce_bc
    use parallel
    use math
    implicit none
    integer :: i,j,st
    
    if (xper.ne.1) then
       ! Left boundary
       ! -> Newmann on P
       ! -> Dirichlet on U,V,W
       if (iproc.eq.1) then
          divv_xx(imin,:,:) = 0.0_WP
          divv_xy(imin,:,:) = 0.0_WP
          divv_xz(imin,:,:) = 0.0_WP
          
          ! Nothing outside
          do i=imin-stv2,imin-1
             divv_u(i,:,:) = 0.0_WP
             divv_v(i,:,:) = 0.0_WP
             divv_w(i,:,:) = 0.0_WP
          end do
          
          ! Dirichlet on U - NOT ZERO
          ! but here we prefer Newmann on U
          do i=imin,imin+stv1-1
             st = imin-i
             divv_u  (i,:,:) = 0.0_WP
             grad_u_x(i,:,:) = 0.0_WP
             do j=jmin_-stv2,jmax_+stv1
                call hofdd_n(stv2-st+1,x(imin:i+stv2),xm(i),divv_u  (i,j,st:stv2),x(imin))
                call hofdd_n(stv2-st+1,x(imin:i+stv2),xm(i),grad_u_x(i,j,st:stv2),x(imin))
             end do
          end do
          
          ! Newmann on V & W
          do i=imin-stv1,imin
             grad_v_x(i,:,:) = 0.0_WP
             grad_w_x(i,:,:) = 0.0_WP
          end do
          do i=imin+1,imin+stv2
             st = imin-1-i
             grad_v_x(i,:,:) = 0.0_WP
             grad_w_x(i,:,:) = 0.0_WP
             do j=jmin_-stv1,jmax_+stv2
                call hofdd_n(stv1-st,xm(imin:i+stv1),x(i),grad_v_x(i,j,st+1:stv1),x(imin))
                call hofdd_n(stv1-st,xm(imin:i+stv1),x(i),grad_w_x(i,j,st+1:stv1),x(imin))
             end do
          end do
          
          ! Reduced stencil
          do i=imin+1,imin+stv1
             st = imin-1-i
             divv_xx(i,:,:) = 0.0_WP
             do j=jmin_,jmax_
                call hofdd(stv1-st,xm(imin:i+stv1),x(i),divv_xx(i,j,st+1:stv1))
             end do
          end do
          do i=imin,imin+stv1
             st = imin-i
             divv_yx(i,:,:) = 0.0_WP
             divv_zx(i,:,:) = 0.0_WP
             do j=jmin_,jmax_
                call hofdd(stv2-st+1,x(imin:i+stv2),xm(i),divv_yx(i,j,st:stv2))
                call hofdd(stv2-st+1,x(imin:i+stv2),xm(i),divv_zx(i,j,st:stv2))
             end do
          end do
       end if

       ! Right boundary
       ! -> Newmann on P
       ! -> Dirichlet on U,V,W
       if (iproc.eq.npx) then
          
          ! Nothing outside
          do i=imax+1,imax+stv1
             divv_u(i,:,:) = 0.0_WP
             divv_v(i,:,:) = 0.0_WP
             divv_w(i,:,:) = 0.0_WP
          end do
          
          ! Dirichlet on U - NOT ZERO
          ! but here we prefer Newmann on U
          do i=imax+1-stv1,imax
             st = imax+1-i
             divv_u  (i,:,:) = 0.0_WP
             grad_u_x(i,:,:) = 0.0_WP
             do j=jmin_-stv2,jmax_+stv1
                call hofdd_n(st+stv1+1,x(i-stv1:imax+1),xm(i),divv_u  (i,j,-stv1:st),x(imax+1))
                call hofdd_n(st+stv1+1,x(i-stv1:imax+1),xm(i),grad_u_x(i,j,-stv1:st),x(imax+1))
             end do
          end do
          
          ! Newmann on V & W
          do i=imax+1,imax+stv2
             grad_v_x(i,:,:) = 0.0_WP
             grad_w_x(i,:,:) = 0.0_WP
          end do
          do i=imax-stv1,imax
             st = imax+1-i
             grad_v_x(i,:,:) = 0.0_WP
             grad_w_x(i,:,:) = 0.0_WP
             do j=jmin_-stv1,jmax_+stv2
                call hofdd_n(st+stv2,xm(i-stv2:imax),x(i),grad_v_x(i,j,-stv2:st-1),x(imax+1))
                call hofdd_n(st+stv2,xm(i-stv2:imax),x(i),grad_w_x(i,j,-stv2:st-1),x(imax+1))
             end do
          end do
          
          ! Reduced stencil
          do i=imax+1-stv1,imax
             st = imax+1-i
             divv_xx(i,:,:) = 0.0_WP
             do j=jmin_,jmax_
                call hofdd(st+stv2,xm(i-stv2:imax),x(i),divv_xx(i,j,-stv2:st-1))
             end do
          end do
          do i=imax-stv1,imax
             st = imax+1-i
             divv_yx(i,:,:) = 0.0_WP
             divv_zx(i,:,:) = 0.0_WP
             do j=jmin_,jmax_
                call hofdd(st+stv1+1,x(i-stv1:imax+1),xm(i),divv_yx(i,j,-stv1:st))
                call hofdd(st+stv1+1,x(i-stv1:imax+1),xm(i),divv_zx(i,j,-stv1:st))
             end do
          end do
       end if
    end if
    
    if (yper.ne.1) then
       ! Lower boundary - Cartesian
       ! -> Newmann on P
       ! -> Newmann on U/W
       ! -> Dirichlet on V
       if (jproc.eq.1 .and. icyl.eq.0) then
          divv_yx(:,jmin,:) = 0.0_WP
          divv_yy(:,jmin,:) = 0.0_WP
          divv_yz(:,jmin,:) = 0.0_WP
          
          ! Nothing outside
          do j=jmin-stv2,jmin-1
             divv_u(:,j,:) = 0.0_WP
             divv_v(:,j,:) = 0.0_WP
             divv_w(:,j,:) = 0.0_WP
          end do
          
          ! Dirichlet on V
          do j=jmin,jmin+stv1
             st = jmin-j
             divv_v  (:,j,:) = 0.0_WP
             grad_v_y(:,j,:) = 0.0_WP
             do i=imin_-stv2,imax_+stv1
                call hofdd_d(stv2-st,y(jmin+1:j+stv2),ym(j),divv_v  (i,j,st+1:stv2),y(jmin))
                call hofdd_d(stv2-st,y(jmin+1:j+stv2),ym(j),grad_v_y(i,j,st+1:stv2),y(jmin))
!!$                call hofdd_n(stv2-st,y(jmin+1:j+stv2),ym(j),divv_v  (i,j,st+1:stv2),y(jmin))
!!$                call hofdd_n(stv2-st,y(jmin+1:j+stv2),ym(j),grad_v_y(i,j,st+1:stv2),y(jmin))
             end do
          end do
          
          ! Newmann on U & W
          do j=jmin-stv1,jmin
             grad_u_y(:,j,:) = 0.0_WP
             grad_w_y(:,j,:) = 0.0_WP
             interpv_cyl_w_y(:,j,:) = 0.0_WP
          end do
          do j=jmin+1,jmin+stv2
             st = jmin-1-j
             grad_u_y(:,j,:) = 0.0_WP
             grad_w_y(:,j,:) = 0.0_WP
             interpv_cyl_w_y(:,j,:) = 0.0_WP
             do i=imin_-stv1,imax_+stv2
                call hofdd_n(stv1-st,ym(jmin:j+stv1),y(j),grad_u_y(i,j,st+1:stv1),y(jmin))
                call hofdd_n(stv1-st,ym(jmin:j+stv1),y(j),grad_w_y(i,j,st+1:stv1),y(jmin))
                call hofdi_n(stv1-st,ym(jmin:j+stv1),y(j),interpv_cyl_w_y(i,j,st+1:stv1),y(jmin))
             end do
          end do
          
          ! Reduced stencil
          do j=jmin+1,jmin+stv1
             st = jmin-1-j
             divv_yy(:,j,:) = 0.0_WP
             interpv_cyl_F_y(:,j,:) = 0.0_WP
             do i=imin_,imax_
                call hofdd(stv1-st,ym(jmin:j+stv1),y(j),divv_yy(i,j,st+1:stv1))
                call hofdi(stv1-st,ym(jmin:j+stv1),y(j),interpv_cyl_F_y(i,j,st+1:stv1))
             end do
          end do
          do j=jmin,jmin+stv1
             st = jmin-j
             divv_xy(:,j,:) = 0.0_WP
             divv_zy(:,j,:) = 0.0_WP
             interpv_cyl_F_ym(:,j,:) = 0.0_WP
             do i=imin_,imax_
                call hofdd(stv2-st+1,y(jmin:j+stv2),ym(j),divv_xy(i,j,st:stv2))
                call hofdd(stv2-st+1,y(jmin:j+stv2),ym(j),divv_zy(i,j,st:stv2))
                call hofdi(stv2-st+1,y(jmin:j+stv2),ym(j),interpv_cyl_F_ym(i,j,st:stv2))
             end do
          end do
       end if
       
       ! Upper boundary
       ! -> Newmann on P
       ! -> Newmann on U/W
       ! -> Dirichlet on V
       if (jproc.eq.npy) then
          
          ! Nothing outside
          do j=jmax+1,jmax+stv1
             divv_u(:,j,:) = 0.0_WP
             divv_v(:,j,:) = 0.0_WP
             divv_w(:,j,:) = 0.0_WP
          end do
          
          ! Dirichlet on V
          do j=jmax-stv2,jmax
             st = jmax+1-j
             divv_v  (:,j,:) = 0.0_WP
             grad_v_y(:,j,:) = 0.0_WP
             do i=imin_-stv2,imax_+stv1
                call hofdd_d(st+stv1,y(j-stv1:jmax),ym(j),divv_v  (i,j,-stv1:st-1),y(jmax+1))
                call hofdd_d(st+stv1,y(j-stv1:jmax),ym(j),grad_v_y(i,j,-stv1:st-1),y(jmax+1))
!!$                call hofdd_n(st+stv1,y(j-stv1:jmax),ym(j),divv_v  (i,j,-stv1:st-1),y(jmax+1))
!!$                call hofdd_n(st+stv1,y(j-stv1:jmax),ym(j),grad_v_y(i,j,-stv1:st-1),y(jmax+1))
             end do
          end do
          
          ! Newmann on U & W
          do j=jmax+1,jmax+stv2
             grad_u_y(:,j,:) = 0.0_WP
             grad_w_y(:,j,:) = 0.0_WP
             interpv_cyl_w_y(:,j,:) = 0.0_WP
          end do
          do j=jmax-stv1,jmax
             st = jmax+1-j
             grad_u_y(:,j,:) = 0.0_WP
             grad_w_y(:,j,:) = 0.0_WP
             interpv_cyl_w_y(:,j,:) = 0.0_WP
             do i=imin_-stv1,imax_+stv2
                call hofdd_n(st+stv2,ym(j-stv2:jmax),y(j),grad_u_y(i,j,-stv2:st-1),y(jmax+1))
                call hofdd_n(st+stv2,ym(j-stv2:jmax),y(j),grad_w_y(i,j,-stv2:st-1),y(jmax+1))
                call hofdi_n(st+stv2,ym(j-stv2:jmax),y(j),interpv_cyl_w_y(i,j,-stv2:st-1),y(jmax+1))
             end do
          end do
          
          ! Reduced stencil
          do j=jmax+1-stv1,jmax
             st = jmax+1-j
             divv_yy(:,j,:) = 0.0_WP
             interpv_cyl_F_y(:,j,:) = 0.0_WP
             do i=imin_,imax_
                call hofdd(st+stv2,ym(j-stv2:jmax),y(j),divv_yy(i,j,-stv2:st-1))
                call hofdi(st+stv2,ym(j-stv2:jmax),y(j),interpv_cyl_F_y(i,j,-stv2:st-1))
             end do
          end do
          do j=jmax-stv1,jmax
             st = jmax+1-j
             divv_xy(:,j,:) = 0.0_WP
             divv_zy(:,j,:) = 0.0_WP
             interpv_cyl_F_ym(:,j,:) = 0.0_WP
             do i=imin_,imax_
                call hofdd(st+stv1+1,y(j-stv1:jmax+1),ym(j),divv_xy(i,j,-stv1:st))
                call hofdd(st+stv1+1,y(j-stv1:jmax+1),ym(j),divv_zy(i,j,-stv1:st))
                call hofdi(st+stv1+1,y(j-stv1:jmax+1),ym(j),interpv_cyl_F_ym(i,j,-stv1:st))
             end do
          end do
       end if
    end if
    
    return
  end subroutine metric_velocity_visc_enforce_bc
  
  
  ! Enforce the dimension of the problem
  ! ------------------------------------
  subroutine metric_velocity_visc_enforce_dim
    use parallel
    implicit none
    
    ! Remove dependance in x and U
    if (nx.eq.1) then
       divv_u = 0.0_WP
       
       divv_xx = 0.0_WP
       divv_yx = 0.0_WP
       divv_zx = 0.0_WP
       divv_xy = 0.0_WP
       divv_xz = 0.0_WP
       
       grad_u_x = 0.0_WP
       grad_u_z = 0.0_WP
       grad_u_y = 0.0_WP
       grad_v_x = 0.0_WP
       grad_w_x = 0.0_WP
    end if
    
    ! Remove dependance in y and V
    if (ny.eq.1) then
       interpv_cyl_v_ym = 0.0_WP
       interpv_cyl_F_y  = 0.0_WP
       
       divv_v = 0.0_WP
       
       divv_xy = 0.0_WP
       divv_yy = 0.0_WP
       divv_zy = 0.0_WP
       divv_yx = 0.0_WP
       divv_yz = 0.0_WP
       
       grad_v_x = 0.0_WP
       grad_v_y = 0.0_WP
       grad_v_z = 0.0_WP
       grad_u_y = 0.0_WP
       grad_w_y = 0.0_WP
    end if
    
    ! Remove dependance in z and W
    if (nz.eq.1) then
       interpv_cyl_w_y  = 0.0_WP
       interpv_cyl_F_ym = 0.0_WP
       
       divv_w = 0.0_WP
       
       divv_xz = 0.0_WP
       divv_yz = 0.0_WP
       divv_zx = 0.0_WP
       divv_zy = 0.0_WP
       divv_zz = 0.0_WP
       
       grad_w_x = 0.0_WP
       grad_w_y = 0.0_WP
       grad_w_z = 0.0_WP
       grad_u_z = 0.0_WP
       grad_v_z = 0.0_WP
       
       if (icyl.eq.1 .and. jproc.eq.1) then
          divv_yx(:,jmin,:) = 0.0_WP
          divv_yy(:,jmin,:) = 0.0_WP
          divv_yz(:,jmin,:) = 0.0_WP
          interpv_cyl_F_y(:,jmin,:) = 0.0_WP
          grad_u_y(:,jmin,:) = 0.0_WP
       end if
    end if
    
    ! Handle sector computation
    if (icyl.eq.1 .and. jproc.eq.1 .and. isect.eq.1) then
       divv_yx(:,jmin,:) = 0.0_WP
       divv_yy(:,jmin,:) = 0.0_WP
       divv_yz(:,jmin,:) = 0.0_WP
       interpv_cyl_F_y(:,jmin,:) = 0.0_WP
       grad_u_y(:,jmin,:) = 0.0_WP
       grad_w_y(:,jmin,:) = 0.0_WP
    end if
    
    return
  end subroutine metric_velocity_visc_enforce_dim
  
  ! Enforce Dirichlet condition
  ! ---------------------------
  subroutine metric_velocity_visc_dirichlet
    use parallel
    use masks
    implicit none
    
    integer :: i,j
    
    do j=jmin_,jmax_
       do i=imin_,imax_

          ! Remove dependance in U
          if (mask_u(i,j).eq.3) then
             divv_xx(i,j,:) = 0.0_WP
             divv_xy(i,j,:) = 0.0_WP
             divv_xz(i,j,:) = 0.0_WP
          end if
          
          ! Remove dependance in V
          if (mask_v(i,j).eq.3) then
             interpv_cyl_F_y (i,j,:) = 0.0_WP
             divv_yx(i,j,:) = 0.0_WP
             divv_yy(i,j,:) = 0.0_WP
             divv_yz(i,j,:) = 0.0_WP
          end if
          
          ! Remove dependance in W
          if (mask_w(i,j).eq.3) then
             interpv_cyl_F_ym(i,j,:) = 0.0_WP
             divv_zx(i,j,:) = 0.0_WP
             divv_zy(i,j,:) = 0.0_WP
             divv_zz(i,j,:) = 0.0_WP
          end if
       end do
    end do
    
    return
  end subroutine metric_velocity_visc_dirichlet
  
end module metric_velocity_visc


subroutine metric_velocity_visc_init
  use metric_velocity_visc
  implicit none
  
  ! Set the stencil lengths
  stv2 = vel_visc_order/2
  stv1 = stv2 - 1
  stvp = stv1 + stv2
  
  ! Allocate all the arrays
  call metric_velocity_visc_allocate
  
  ! Compute the metric
  call metric_velocity_visc_compute_primary
  
  ! Impose BCs
  call metric_velocity_visc_enforce_bc
  call metric_velocity_visc_enforce_walls
  
  ! Enforce some properties
  call metric_velocity_visc_enforce_cyl
  call metric_velocity_visc_enforce_dim
  call metric_velocity_visc_dirichlet
  
  return
end subroutine metric_velocity_visc_init
