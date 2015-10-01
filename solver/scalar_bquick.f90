module scalar_bquick
  use scalar_quick
  use data
  use metric_generic
  implicit none

  ! ================================================= !
  ! Definition of the representation of the metrics   !
  ! ================================================= !
  !                                                   !
  ! 1D interpolation or differentiation               !
  ! [point 1] [point 2] [point 3]            case U>0 !
  !           [point 4] [point 5] [point 6]  case U<0 !
  ! ---------------------------------> x,y,z          !
  !                                                   !
  ! ================================================= !
  
  ! Coefficient for the Bounded QUICK
  real(WP), dimension(:,:,:,:), pointer :: bquick_xp,bquick_yp,bquick_zp
  real(WP), dimension(:,:,:,:), pointer :: bquick_xm,bquick_ym,bquick_zm

  ! Withint physical bounds
  integer, dimension(:,:,:,:), pointer :: bounded
  
  ! Temporary array
  real(WP), dimension(:), pointer :: SC_

  !$OMP THREADPRIVATE(SC_)
  
  ! Values to monitor
  real(WP), dimension(:), pointer :: not_bounded
  
end module scalar_bquick


! ============================ !
! Initialize the BQUICK module !
! ============================ !
subroutine scalar_bquick_init
  use scalar_bquick
  implicit none
  
  integer :: isc
  
  ! Allocate arrays
  allocate(bquick_xp(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-2: 0))
  allocate(bquick_yp(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-2: 0))
  allocate(bquick_zp(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-2: 0))
  allocate(bquick_xm(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-1:+1))
  allocate(bquick_ym(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-1:+1))
  allocate(bquick_zm(imin_-st1-1:imax_+st2+1,jmin_-st1-1:jmax_+st2+1,kmin_-st1-1:kmax_+st2+1,-1:+1))

  allocate(bounded (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  
  ! Initialize the quick scheme
  call scalar_quick_init
  
  ! Allocate temp array
  !$OMP PARALLEL
  allocate(SC_(nscalar))
  !$OMP END PARALLEL
  
  ! Create new file to monitor at each subiterations
  call monitor_create_file_iter('bquick',nscalar)
  allocate(not_bounded(nscalar))
  do isc=1,nscalar
     call monitor_set_header(isc,'Not_Bounded '//SC_name(isc),'i')
  end do
  
  return
end subroutine scalar_bquick_init


! ==================================================== !
! Compute the BQUICK coefficients for the given scalar !
! ==================================================== !
subroutine scalar_bquick_check_bounds
  use scalar_bquick
  use time_info
  use masks
  implicit none
  
  integer  :: i,j,k,isc,tmpi
  real(WP) :: tmp
  
  ! Starting at the third iteration keep the non-bounded points
  if (niter.le.2) then

     !$OMP PARALLEL

     do isc=1,nscalar
        !$OMP DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 bounded(i,j,k,isc) = 1
              end do
           end do
        end do
        !$OMP END DO
     end do

     !$OMP END PARALLEL
     
  end if

  ! Check if scalar is bounded
  do isc=1,nscalar
     !$OMP PARALLEL DO PRIVATE(tmpi)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              SC_ = ResSC(i,j,k,:)
              call combustion_check_bounds(SC_,i,j,k,isc,tmpi)
              !call spray_check_bounds(SC_,i,j,k,isc,tmpi)
              !call pollutants_check_bounds(SC_,i,j,k,isc,tmpi)
              bounded(i,j,k,isc) = min(bounded(i,j,k,isc),tmpi)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call boundary_update_border_int(bounded(:,:,:,isc),'+','ym')     
  end do
  
  ! Force bounded in walls
  !$OMP PARALLEL DO
  do j=jmin_,jmax_
     do i=imin_,imax_
        if (mask(i,j).ne.0) bounded(i,j,:,:) = 1
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Count the points not bounded
  do isc=1,nscalar
     tmp = real(nx_*ny_*nz_-sum(bounded(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)),WP)
     call parallel_sum(tmp,not_bounded(isc))
  end do
  
  ! Transfer values to monitor
  call monitor_select_file('bquick')
  call monitor_set_array_values(not_bounded)
  
  return
end subroutine scalar_bquick_check_bounds


! ==================================================== !
! Compute the BQUICK coefficients for the given scalar !
! ==================================================== !
subroutine scalar_bquick_coeff(isc)
  use scalar_bquick
  implicit none
  
  integer , intent(in) :: isc
  integer :: i,j,k
  
  ! By default, set the coefficients to quick
  !$OMP DO
  do k=kmin_-st1,kmax_+st2
     do j=jmin_-st1,jmax_+st2
        do i=imin_-st1,imax_+st2
           bquick_xp(i,j,k,:) = quick_xp(i,j,:)
           bquick_yp(i,j,k,:) = quick_yp(i,j,:)
           bquick_zp(i,j,k,:) = quick_zp(i,j,:)
           bquick_xm(i,j,k,:) = quick_xm(i,j,:)
           bquick_ym(i,j,k,:) = quick_ym(i,j,:)
           bquick_zm(i,j,k,:) = quick_zm(i,j,:)
        end do
     end do
  end do
  !$OMP END DO
  
  ! Switch to first order if locally out of bounds
  !$OMP DO
  do k=kmin_-1,kmax_+1
     do j=jmin_-1,jmax_+1
        do i=imin_-1,imax_+1
           if (bounded(i,j,k,isc).eq.0) then
              ! Direction x - U>0
              bquick_xp(i:i+1,j,k,-2) = 0.0_WP
              bquick_xp(i:i+1,j,k,-1) = 1.0_WP
              bquick_xp(i:i+1,j,k, 0) = 0.0_WP
              ! Direction x - U<0
              bquick_xm(i:i+1,j,k,-1) = 0.0_WP
              bquick_xm(i:i+1,j,k, 0) = 1.0_WP
              bquick_xm(i:i+1,j,k,+1) = 0.0_WP
              ! Direction y - V>0
              bquick_yp(i,j:j+1,k,-2) = 0.0_WP
              bquick_yp(i,j:j+1,k,-1) = 1.0_WP
              bquick_yp(i,j:j+1,k, 0) = 0.0_WP
              ! Direction y - V<0
              bquick_ym(i,j:j+1,k,-1) = 0.0_WP
              bquick_ym(i,j:j+1,k, 0) = 1.0_WP
              bquick_ym(i,j:j+1,k,+1) = 0.0_WP
              ! Direction z - W>0
              bquick_zp(i,j,k:k+1,-2) = 0.0_WP
              bquick_zp(i,j,k:k+1,-1) = 1.0_WP
              bquick_zp(i,j,k:k+1, 0) = 0.0_WP
              ! Direction z - W<0
              bquick_zm(i,j,k:k+1,-1) = 0.0_WP
              bquick_zm(i,j,k:k+1, 0) = 1.0_WP
              bquick_zm(i,j,k:k+1,+1) = 0.0_WP
           end if
        end do
     end do
  end do
  !$OMP END DO
  
  return
end subroutine scalar_bquick_coeff



! =========================================================== !
! Compute the residuals of the scalar equations               !
!                                                             !
! - velocity field n+1 stored in U/rhoU                       !
! - scalar field n+1 stored in SCmid                          !
!                                                             !
! 3 working arrays of size (at least) (nx_+1)*(ny_+1)*(nz_+1) !
! =========================================================== !
subroutine scalar_bquick_residual
  use scalar_bquick
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: rhs

  !$OMP PARALLEL

  ! First, Get an estimation with the quick scheme
  do isc=1,nscalar
     !$OMP DO
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              
              FX(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoUt(i,j,k,isc)+abs(rhoUt(i,j,k,isc))) * sum(quick_xp(i,j,:)*SC(i-2:i  ,j,k,isc)) &
                   - 0.5_WP*(rhoUt(i,j,k,isc)-abs(rhoUt(i,j,k,isc))) * sum(quick_xm(i,j,:)*SC(i-1:i+1,j,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_x(i,j,:)*DIFF(i-st2:i+st1,j,k,isc)) * &
                     sum(grad_x(i,j,:)*SC(i-st2:i+st1,j,k,isc))
              
              FY(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoVt(i,j,k,isc)+abs(rhoVt(i,j,k,isc))) * sum(quick_yp(i,j,:)*SC(i,j-2:j  ,k,isc)) &
                   - 0.5_WP*(rhoVt(i,j,k,isc)-abs(rhoVt(i,j,k,isc))) * sum(quick_ym(i,j,:)*SC(i,j-1:j+1,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_y(i,j,:)*DIFF(i,j-st2:j+st1,k,isc)) * &
                     sum(grad_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
              
              FZ(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoWt(i,j,k,isc)+abs(rhoWt(i,j,k,isc))) * sum(quick_zp(i,j,:)*SC(i,j,k-2:k  ,isc)) &
                   - 0.5_WP*(rhoWt(i,j,k,isc)-abs(rhoWt(i,j,k,isc))) * sum(quick_zm(i,j,:)*SC(i,j,k-1:k+1,isc)) &
                   ! Viscous term
                   + sum(interp_sc_z(i,j,:)*DIFF(i,j,k-st2:k+st1,isc)) * &
                     sum(grad_z(i,j,:)*SC(i,j,k-st2:k+st1,isc))
              
           end do
        end do
     end do
     !$OMP END DO
     
     !$OMP DO PRIVATE(rhs)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              
              rhs =+sum(div_u(i,j,:)*FX(i-st1:i+st2,j,k)) &
                   +sum(div_v(i,j,:)*FY(i,j-st1:j+st2,k)) &
                   +sum(div_w(i,j,:)*FZ(i,j,k-st1:k+st2)) 

              ResSC(i,j,k,isc) =  -2.0_WP*SC(i,j,k,isc)+SCold(i,j,k,isc) &
                   + ( RHOold(i,j,k)*SCold(i,j,k,isc) + dt*rhs &
                   + srcSCmid(i,j,k,isc) + srcSCfull(i,j,k,isc) ) /RHO(i,j,k)
              
              ResSC(i,j,k,isc) = 2.0_WP*SC(i,j,k,isc)-SCold(i,j,k,isc) + ResSC(i,j,k,isc)

           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL
  
  ! Check the bounds
  call scalar_bquick_check_bounds

  !$OMP PARALLEL

  do isc=1,nscalar   
     ! Compute the coefficients
     call scalar_bquick_coeff(isc)
     
     ! Second, compute the real residual with the corrected coefficients
     !$OMP DO
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              
              FX(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoUt(i,j,k,isc)+abs(rhoUt(i,j,k,isc))) * sum(bquick_xp(i,j,k,:)*SC(i-2:i  ,j,k,isc)) &
                   - 0.5_WP*(rhoUt(i,j,k,isc)-abs(rhoUt(i,j,k,isc))) * sum(bquick_xm(i,j,k,:)*SC(i-1:i+1,j,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_x(i,j,:)*DIFF(i-st2:i+st1,j,k,isc)) * &
                     sum(grad_x(i,j,:)*SC(i-st2:i+st1,j,k,isc))
              
              FY(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoVt(i,j,k,isc)+abs(rhoVt(i,j,k,isc))) * sum(bquick_yp(i,j,k,:)*SC(i,j-2:j  ,k,isc)) &
                   - 0.5_WP*(rhoVt(i,j,k,isc)-abs(rhoVt(i,j,k,isc))) * sum(bquick_ym(i,j,k,:)*SC(i,j-1:j+1,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_y(i,j,:)*DIFF(i,j-st2:j+st1,k,isc)) * &
                     sum(grad_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
              
              FZ(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoWt(i,j,k,isc)+abs(rhoWt(i,j,k,isc))) * sum(bquick_zp(i,j,k,:)*SC(i,j,k-2:k  ,isc)) &
                   - 0.5_WP*(rhoWt(i,j,k,isc)-abs(rhoWt(i,j,k,isc))) * sum(bquick_zm(i,j,k,:)*SC(i,j,k-1:k+1,isc)) &
                   ! Viscous term
                   + sum(interp_sc_z(i,j,:)*DIFF(i,j,k-st2:k+st1,isc)) * &
                     sum(grad_z(i,j,:)*SC(i,j,k-st2:k+st1,isc))
              
           end do
        end do
     end do
     !$OMP END DO
     
     !$OMP DO PRIVATE(rhs)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              
              rhs =+sum(div_u(i,j,:)*FX(i-st1:i+st2,j,k)) &
                   +sum(div_v(i,j,:)*FY(i,j-st1:j+st2,k)) &
                   +sum(div_w(i,j,:)*FZ(i,j,k-st1:k+st2)) 

              ResSC(i,j,k,isc) =  -2.0_WP*SC(i,j,k,isc)+SCold(i,j,k,isc) &
                   + ( RHOold(i,j,k)*SCold(i,j,k,isc) + dt*rhs &
                   + srcSCmid(i,j,k,isc) + srcSCfull(i,j,k,isc) ) /RHO(i,j,k)
              
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL

  return
end subroutine scalar_bquick_residual


! =========================================================== !
! Inverse the linear system obtained from the implicit scalar !
! transport equation                                          !
! =========================================================== !
subroutine scalar_bquick_inverse
  use scalar_bquick
  use parallel
  use memory
  use implicit
  use time_info
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: conv1,conv2,conv3,conv4
  real(WP) :: visc1,visc2
  real(WP) :: dt2
  
  dt2 = dt/2.0_WP
  
  ! If purely explicit return
  if (.not.implicit_any) return
  
  do isc=1,nscalar
     
     ! Get coefficients first
     !$OMP PARALLEL
     call scalar_bquick_coeff(isc)
     !$OMP END PARALLEL
     
     ! X-direction
     if (implicit_x) then

        !$OMP PARALLEL DO PRIVATE(conv1,conv2,conv3,conv4,visc1,visc2)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoUt(i  ,j,k,isc) + abs(rhoUt(i  ,j,k,isc)))
                 conv2 = 0.5_WP*(rhoUt(i  ,j,k,isc) - abs(rhoUt(i  ,j,k,isc)))
                 conv3 = 0.5_WP*(rhoUt(i+1,j,k,isc) + abs(rhoUt(i+1,j,k,isc)))
                 conv4 = 0.5_WP*(rhoUt(i+1,j,k,isc) - abs(rhoUt(i+1,j,k,isc)))
                 
                 ! No convective implicit correction at the axis => tight coupling x-z
                 if (j.eq.jmin) then
                    conv1 = 0.0_WP
                    conv2 = 0.0_WP
                    conv3 = 0.0_WP
                    conv4 = 0.0_WP
                 end if
                 
                 visc1 = sum(interp_sc_x(i  ,j,:)*DIFF(i  -st2:i  +st1,j,k,isc))
                 visc2 = sum(interp_sc_x(i+1,j,:)*DIFF(i+1-st2:i+1+st1,j,k,isc))
                 
                 Ax(j,k,i,-2) = + dt2 * div_u(i,j,0)*conv1*bquick_xp(i,j,k,-2)
                 
                 Ax(j,k,i,-1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*bquick_xp(i  ,j,k,-1) - conv2*bquick_xm(i  ,j,k,-1) + visc1*grad_x(i,j,-1)) &
                      + div_u(i,j,1)*( - conv3*bquick_xp(i+1,j,k,-2)))
                 
                 Ax(j,k,i, 0) = RHO(i,j,k) - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*bquick_xp(i  ,j,k, 0) - conv2*bquick_xm(i  ,j,k, 0) + visc1*grad_x(i,j,0)) &
                      + div_u(i,j,1)*( - conv3*bquick_xp(i+1,j,k,-1) - conv4*bquick_xm(i+1,j,k,-1) + visc2*grad_x(i+1,j,-1)))
                 
                 Ax(j,k,i,+1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv2*bquick_xm(i  ,j,k,+1)) &
                      + div_u(i,j,1)*( - conv3*bquick_xp(i+1,j,k, 0) - conv4*bquick_xm(i+1,j,k, 0) + visc2*grad_x(i+1,j,0)))
                 
                 Ax(j,k,i,+2) = + dt2 * div_u(i,j,1)*conv4*bquick_xm(i+1,j,k,+1)
                  
                 Rx(j,k,i) = RHO(i,j,k) * ResSC(i,j,k,isc)
                 
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        call implicit_solve_x(5)
     else
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 Rx(j,k,i) = ResSC(i,j,k,isc)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     
     ! Y-direction
     if (implicit_y) then

        !$OMP PARALLEL PRIVATE(conv1,conv2,conv3,conv4,visc1,visc2)

        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoVt(i,j  ,k,isc) + abs(rhoVt(i,j  ,k,isc)))
                 conv2 = 0.5_WP*(rhoVt(i,j  ,k,isc) - abs(rhoVt(i,j  ,k,isc)))
                 conv3 = 0.5_WP*(rhoVt(i,j+1,k,isc) + abs(rhoVt(i,j+1,k,isc)))
                 conv4 = 0.5_WP*(rhoVt(i,j+1,k,isc) - abs(rhoVt(i,j+1,k,isc)))
                 
                 visc1 = sum(interp_sc_y(i,j  ,:)*DIFF(i,j  -st2:j  +st1,k,isc))
                 visc2 = sum(interp_sc_y(i,j+1,:)*DIFF(i,j+1-st2:j+1+st1,k,isc))
                 
                 Ay(i,k,j,-2) = + dt2 * div_v(i,j,0)*conv1*bquick_yp(i,j,k,-2)
                 
                 Ay(i,k,j,-1) = - dt2*( &
                      + div_v(i,j,0)*( - conv1*bquick_yp(i,j  ,k,-1) - conv2*bquick_ym(i,j  ,k,-1) + visc1*grad_y(i,j,-1)) &
                      + div_v(i,j,1)*( - conv3*bquick_yp(i,j+1,k,-2)))
                           
                 Ay(i,k,j, 0) = RHO(i,j,k) - dt2 * ( &
                      + div_v(i,j,0)*( - conv1*bquick_yp(i,j  ,k, 0) - conv2*bquick_ym(i,j  ,k, 0) + visc1*grad_y(i,j,0)) &
                      + div_v(i,j,1)*( - conv3*bquick_yp(i,j+1,k,-1) - conv4*bquick_ym(i,j+1,k,-1) + visc2*grad_y(i,j+1,-1)))
                 
                 Ay(i,k,j,+1) = - dt2*( &
                      + div_v(i,j,0)*( - conv2*bquick_ym(i,j  ,k,+1)) &
                      + div_v(i,j,1)*( - conv3*bquick_yp(i,j+1,k, 0) - conv4*bquick_ym(i,j+1,k, 0) + visc2*grad_y(i,j+1,0)))
                 
                 Ay(i,k,j,+2) = dt2 * div_v(i,j,1)*conv4*bquick_ym(i,j+1,k,+1)
                 
                 Ry(i,k,j) = RHO(i,j,k) * Rx(j,k,i)
                 
              end do
           end do
        end do
        !$OMP END DO
        if (icyl.eq.1 .and. jproc.eq.1) then
           !$OMP DO
           do k=kmin_,kmax_
              Ay(:,k,jmin+1,-1) = Ay(:,k,jmin+1,-1) + Ay(:,k,jmin+1,-2)
              Ay(:,k,jmin+1,-2) = 0.0_WP
              Ay(:,k,jmin, 0) = Ay(:,k,jmin, 0) + Ay(:,k,jmin,-1) + Ay(:,k,jmin,-2)
              Ay(:,k,jmin,-1) = 0.0_WP
              Ay(:,k,jmin,-2) = 0.0_WP
           end do
           !$OMP END DO
        end if

        !$OMP END PARALLEL

        call implicit_solve_y(5)
     else 
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 Ry(i,k,j) = Rx(j,k,i)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     
     ! Z-direction
     if (implicit_z) then
        !$OMP PARALLEL DO PRIVATE(conv1,conv2,conv3,conv4,visc1,visc2)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoWt(i,j,k  ,isc) + abs(rhoWt(i,j,k  ,isc)))
                 conv2 = 0.5_WP*(rhoWt(i,j,k  ,isc) - abs(rhoWt(i,j,k  ,isc)))
                 conv3 = 0.5_WP*(rhoWt(i,j,k+1,isc) + abs(rhoWt(i,j,k+1,isc)))
                 conv4 = 0.5_WP*(rhoWt(i,j,k+1,isc) - abs(rhoWt(i,j,k+1,isc)))
                 
                 visc1 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k  -st2:k  +st1,isc))
                 visc2 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k+1-st2:k+1+st1,isc))
                 
                 Az(i,j,k,-2) = + dt2 * div_w(i,j,0)*conv1*bquick_zp(i,j,k,-2)
                 
                 Az(i,j,k,-1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*bquick_zp(i,j,k  ,-1) - conv2*bquick_zm(i,j,k,-1) + visc1*grad_z(i,j,-1)) &
                      + div_w(i,j,1)*( - conv3*bquick_zp(i,j,k+1,-2)))
                 
                 Az(i,j,k, 0) = RHO(i,j,k) - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*bquick_zp(i,j,k  , 0) - conv2*bquick_zm(i,j,k  , 0) + visc1*grad_z(i,j,0)) &
                      + div_w(i,j,1)*( - conv3*bquick_zp(i,j,k+1,-1) - conv4*bquick_zm(i,j,k+1,-1) + visc2*grad_z(i,j,-1)))
                 
                 Az(i,j,k,+1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv2*bquick_zm(i,j,k  ,+1)) &
                      + div_w(i,j,1)*( - conv3*bquick_zp(i,j,k+1, 0) - conv4*bquick_zm(i,j,k+1, 0) + visc2*grad_z(i,j,0)))
                 
                 Az(i,j,k,+2) = + dt2 * div_w(i,j,1)*conv4*bquick_zm(i,j,k+1,+1)
                 
                 Rz(i,j,k) = RHO(i,j,k) * Ry(i,k,j)
                 
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        call implicit_solve_z(5)
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 ResSC(i,j,k,isc) = Rz(i,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     else
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 ResSC(i,j,k,isc) = Ry(i,k,j)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     
  end do
  
  return
end subroutine scalar_bquick_inverse
