module bodyforce
  use precision
  use geometry
  use borders
  use partition
  implicit none

  ! Body force for channels and pipes
  real(WP), dimension(:), pointer :: Fx
  real(WP), dimension(:), pointer :: Fz
  real(WP), dimension(:), pointer :: umean
  real(WP), dimension(:), pointer :: wmean
  real(WP), dimension(:), allocatable :: ubulk
  real(WP), dimension(:), allocatable :: wbulk
  
  ! Forcing coefficient for HIT
  real(WP) :: Ahit

  ! Body force for boundary layers
  real(WP), dimension(:), pointer :: Uy
  real(WP) :: deltas,deltas0
  
  ! Gravity term
  logical :: use_gravity
  real(WP), dimension(3) :: gravity
  
contains
  
  ! Compute friction force on the walls
  subroutine compute_friction
    use data
    use metric_generic
    use metric_velocity_visc
    use parallel
    use masks
    implicit none
    
    integer :: i,j,k,nflow
    real(WP) :: tmp_x,tmp_z
    
    do nflow=1,ninlet
       
       ! Lower wall
       tmp_x = 0.0_WP
       tmp_z = 0.0_WP
       if (inlet(nflow)%jmin.ge.jmin_ .and. inlet(nflow)%jmin.le.jmax_) then 
          j = max(inlet(nflow)%jmin,jmin_)
          if (mask(imin_,j-1).eq.1) then
             !$OMP PARALLEL DO REDUCTION(+:tmp_x,tmp_z)
             do k=kmin_,kmax_
                do i=imin_,imax_
                   tmp_x = tmp_x + dz_v(j) * dxm(i-1) * &
                        sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k))* &
                        sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k))

                   tmp_z = tmp_z + dz_v(j) * dxm(i-1) * &
                        sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k-st2:k+st1)) * &
                        ( sum(grad_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) &
                        - ymmi(j)*sum(interpv_cyl_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) )
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       end if

       ! Upper wall
       if (inlet(nflow)%jmax+1.ge.jmin_ .and. inlet(nflow)%jmax+1.le.jmax_) then
          j = min(inlet(nflow)%jmax+1,jmax_+1)
          if (mask(imin_,j).eq.1) then
             !$OMP PARALLEL DO REDUCTION(+:tmp_x,tmp_z)
             do k=kmin_,kmax_
                do i=imin_,imax_
                   tmp_x = tmp_x - dz_v(j) * dxm(i-1) * &
                        sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k))* &
                        sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k))

                   tmp_z = tmp_z - dz_v(j) * dxm(i-1) * &
                        sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k-st2:k+st1)) * &
                        ( sum(grad_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) &
                        - ymmi(j)*sum(interpv_cyl_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) )
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       end if

       call parallel_sum(tmp_x,Fx(nflow))
       call parallel_sum(tmp_z,Fz(nflow))

       Fx(nflow) = Fx(nflow)/inlet(nflow)%A
       Fz(nflow) = Fz(nflow)/inlet(nflow)%A
    end do
    
    Fx  = Fx /(xm(imax)-xm(imin-1))
    Fz  = Fz /(xm(imax)-xm(imin-1))

    return
  end subroutine compute_friction
  
  ! Compute mass flow rate between the walls
  subroutine compute_massflowrate
    use parallel
    use data
    implicit none
    
    integer :: i,j,k,nflow
    real(WP) :: tmp_x
    real(WP) :: tmp_z
    
    do nflow=1,ninlet
       tmp_x = 0.0_WP
       tmp_z = 0.0_WP
       
       !$OMP PARALLEL DO REDUCTION(+:tmp_x,tmp_z)
       do j = max(inlet(nflow)%jmin,jmin_),min(inlet(nflow)%jmax,jmax_)
          do k=kmin_,kmax_
             do i=imin_,imax_
                tmp_x = tmp_x + rhoU(i,j,k)*dA(j)*dxm(i-1)
                tmp_z = tmp_z + rhoW(i,j,k)*dA(j)*dxm(i-1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       
       call parallel_sum(tmp_x,umean(nflow))
       call parallel_sum(tmp_z,wmean(nflow))

       umean(nflow) = umean(nflow)/inlet(nflow)%A
       wmean(nflow) = wmean(nflow)/inlet(nflow)%A
    end do
    
    umean = umean/(xm(imax)-xm(imin-1))
    wmean = wmean/(xm(imax)-xm(imin-1))
    
    return
  end subroutine compute_massflowrate
  
  ! Compute momentum thickness
  subroutine compute_disp_thickness
    use parallel
    use data
    implicit none
    
    integer  :: i,j,k
    real(WP) :: tmp,Uinft
    
    ! Compute mean profile
    Uy = 0.0_WP
    do j = jmin_,jmax_
       do k=kmin_,kmax_
          do i=imin_,imax_
             Uy(j) = Uy(j) + U(i,j,k)
          end do
       end do
    end do
    do j=jmin,jmax
       call parallel_sum(Uy(j),tmp)
       Uy(j) = tmp/real(nx*nz)
    end do
    
    ! Normalize
    Uinft = Uy(jmax) + epsilon(1.0_WP)
    Uy = Uy / Uinft

    ! Compute momentum thickness and gamma terms
    deltas = 0.0_WP
    do j=jmin,jmax
       deltas = deltas + (1.0_WP-Uy(j))*dy(j)
    end do

    return
  end subroutine compute_disp_thickness
  

  ! ================================================== !
  ! Source term for pipes and channels                 !
  !   -> relaxation towards mean bluk velocity (alpha) !
  !   -> counter pressure drop due to friction         !
  ! ================================================== !
  subroutine bodyforce_pipe_channel
    use velocity
    implicit none
    
    integer :: j,nflow
    real(WP), parameter :: alpha = 1.0_WP

    call compute_friction
    call compute_massflowrate
    
    do nflow=1,ninlet
       !$OMP PARALLEL DO
       do j = max(inlet(nflow)%jmin,jmin_),min(inlet(nflow)%jmax,jmax_)
          srcU(:,j,:) = alpha*(ubulk(nflow)-umean(nflow)) + dt_uvw*Fx(nflow)
          srcW(:,j,:) = alpha*(wbulk(nflow)-wmean(nflow)) + dt_uvw*Fz(nflow)
       end do
       !$OMP END PARALLEL DO
    end do
    
    return
  end subroutine bodyforce_pipe_channel
  

  ! ================================================== !
  ! Source term for boundary layers                    !
  !   -> relaxation towards mean momentum thickness    !
  !   -> counter pressure drop due to friction         !
  ! ================================================== !
  subroutine bodyforce_boundary_layer
    use velocity
    implicit none
    
    real(WP), parameter :: alpha = 1.0_WP
    integer :: j

    call compute_friction
    call compute_massflowrate
    !call compute_disp_thickness
    
    do j = max(inlet(1)%jmin,jmin_),min(inlet(1)%jmax,jmax_)
       srcU(:,j,:) = alpha*(ubulk(1)-umean(1)) + dt_uvw*Fx(1)
    end do
    
    return
  end subroutine bodyforce_boundary_layer
  

  ! ================================================== !
  ! Source term for forced isotropic turbulence        !
  !  -> Meneveau, with "improvements" by Blanquart     !                                  
  ! ================================================== !
  subroutine bodyforce_hit
    use data
    use velocity
    use metric_generic
    use parser
    implicit none
    
    real(WP) :: srcUtmp,srcVtmp,srcWtmp
    integer :: i,j,k
    real(WP) :: RHOmean, buf1, buf2, tke, tke_0

    call parser_read("Initial TKE", tke_0, 0.0_WP)

    ! Favre(urms^2) = <rhoU.U/3> / <RHO>
    buf1 = 0.0_WP
    buf2 = 0.0_WP
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             buf1 = buf1 + 0.5_WP * vol(i,j) * ( &
                  sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)*rhoU(i-st1:i+st2,j,k)) + &
                  sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)*rhoV(i,j-st1:j+st2,k)) + &
                  sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)*rhoW(i,j,k-st1:k+st2)) )
             buf2 = buf2 + RHO(i,j,k) * vol(i,j)
          end do
       end do
    end do
    call parallel_sum(buf1,tke)
    tke = tke/vol_total
    call parallel_sum(buf2,RHOmean)
    RHOmean = RHOmean/vol_total
    tke = tke/RHOmean
    
    ! Following Meneveau PoF 17 (2005)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             srcU(i,j,k) = dt_uvw*Ahit*rhoU(i,j,k)*(tke_0/tke)
             srcV(i,j,k) = dt_uvw*Ahit*rhoV(i,j,k)*(tke_0/tke)
             srcW(i,j,k) = dt_uvw*Ahit*rhoW(i,j,k)*(tke_0/tke)
          end do
       end do
    end do
    
    ! Ensure that <srcU/V/W> = 0
    call parallel_sum(sum(srcU(imin_:imax_,jmin_:jmax_,kmin_:kmax_)),srcUtmp)
    call parallel_sum(sum(srcV(imin_:imax_,jmin_:jmax_,kmin_:kmax_)),srcVtmp)
    call parallel_sum(sum(srcW(imin_:imax_,jmin_:jmax_,kmin_:kmax_)),srcWtmp)
    srcU = srcU - srcUtmp/real(nx*ny*nz,WP)
    srcV = srcV - srcVtmp/real(nx*ny*nz,WP)
    srcW = srcW - srcWtmp/real(nx*ny*nz,WP)
    
    return
  end subroutine bodyforce_hit
  
  
  ! =================== !
  ! Gravity source term !
  ! =================== !
  subroutine bodyforce_gravity
    use data
    use velocity
    use masks
    use metric_generic
    
    integer :: i,j,k
    real(WP) :: RHOx,RHOy,RHOz
    
    !$OMP PARALLEL DO PRIVATE(RHOx,RHOy,RHOz)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             RHOx = sum(interp_sc_x(i,j,:)*RHO(i-st2:i+st1,j,k))
             RHOy = sum(interp_sc_y(i,j,:)*RHO(i,j-st2:j+st1,k))
             RHOz = sum(interp_sc_z(i,j,:)*RHO(i,j,k-st2:k+st1))
             if (mask_u(i,j).eq.0) srcU(i,j,k) = srcU(i,j,k)+dt_uvw*RHOx*gravity(1)
             if (mask_v(i,j).eq.0) srcV(i,j,k) = srcV(i,j,k)+dt_uvw*RHOy*gravity(2)
             if (mask_w(i,j).eq.0) srcW(i,j,k) = srcW(i,j,k)+dt_uvw*RHOz*gravity(3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine bodyforce_gravity
  
end module bodyforce


! ===================== !
! INITIALIZE the module !
! ===================== !
subroutine bodyforce_init
  use bodyforce
  use parser
  use config
  implicit none
  
  ! Channels and Pipes and boundary layers
  if ( trim(simu_type).eq."channel" .or. trim(simu_type).eq."pipe") then
     
     ! Allocate the arrays
     allocate(Fx(ninlet))
     allocate(Fz(ninlet))
     allocate(umean(ninlet))
     allocate(wmean(ninlet))
     allocate(ubulk(ninlet))
     allocate(wbulk(ninlet))
     
     ! Compute initial mass flow rates and impose them 
     call compute_massflowrate
     ubulk = umean
     wbulk = wmean
     
  end if
  
  ! HIT
  if (trim(simu_type).eq."hit") then
     call parser_read('Forcing coefficient',Ahit,0.0_WP)
  end if
  
  ! Boundary Layer
  if (trim(simu_type).eq."boundary layer") then

     ! Allocate the arrays
     allocate(Fx(ninlet))
     allocate(Fz(ninlet))
     allocate(Uy(jmin:jmax))

     ! Compute inital displacement thickness
     call compute_disp_thickness
     deltas0 = deltas

     allocate(umean(ninlet))
     allocate(wmean(ninlet))
     allocate(ubulk(ninlet))
     allocate(wbulk(ninlet))
     
     ! Compute initial mass flow rates and impose them 
     call compute_massflowrate
     ubulk = umean
     wbulk = wmean

  end if
  
  ! Gravity
  gravity = 0.0_WP
  call parser_is_defined('Gravity',use_gravity)
  if (use_gravity) call parser_read('Gravity',gravity)
  if (icyl.eq.1) then
     if (gravity(2).ne.0.0_WP .or. gravity(3).ne.0.0_WP) then
        call die("bodyforce_init: in cylindrical coordinates, gravity must be aligned with x.")
     end if
  end if
  
  return
end subroutine bodyforce_init


! ================================================= !
! Compute the SOURCE TERM for the momentum equation !
! ================================================= !
subroutine bodyforce_src
  use bodyforce
  use config
  implicit none
  
  ! Channels and Pipes
  if ( trim(simu_type).eq."channel" .or. trim(simu_type).eq."pipe") then
     call bodyforce_pipe_channel
  end if
  
  ! HIT - Meneveau forcing
  if (trim(simu_type).eq."hit") then
     call bodyforce_hit
  end if
  
  ! Boundary Layer
  if (trim(simu_type).eq."boundary layer") then
     call bodyforce_boundary_layer
  end if
  
  ! Gravity
  if (use_gravity) call bodyforce_gravity
  
  return
end subroutine bodyforce_src
