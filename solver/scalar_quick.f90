module scalar_quick
  use scalar
  use data
  use metric_generic
  implicit none

  real(WP), dimension(:,:,:), pointer :: quick_xp,quick_yp,quick_zp
  real(WP), dimension(:,:,:), pointer :: quick_xm,quick_ym,quick_zm
  
end module scalar_quick


! =========================================== !
! Initialize the metrics for the QUICK scheme !
! =========================================== !
subroutine scalar_quick_init
  use scalar_quick
  use masks
  use math
  implicit none
  integer :: i,j
  
  ! Allocate the arrays
  allocate(quick_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+0))
  allocate(quick_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+0))
  allocate(quick_zp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+0))
  allocate(quick_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(quick_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(quick_zm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  
  !$OMP PARALLEL

  ! Interpolation at the faces
  !$OMP DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        
        call poly3_coeff(x(i-2:i+1),x(i),quick_xp(i,j,:))
        call poly3_coeff(x(i-1:i+2),x(i),quick_xm(i,j,:))
        
        if (mask(i-2,j).eq.1) then
           quick_xp(i,j,-1) = quick_xp(i,j,-1) + quick_xp(i,j,-2)
           quick_xp(i,j,-2) = 0.0_WP
        end if
        if (mask(i-1,j).eq.1) then
           quick_xp(i,j,-2) = 0.0_WP
           quick_xp(i,j,-1) = 0.0_WP
           quick_xp(i,j, 0) = 1.0_WP
           quick_xm(i,j, 0) = quick_xm(i,j, 0) + quick_xm(i,j,-1)
           quick_xm(i,j,-1) = 0.0_WP
        end if
        if (mask(i,j).eq.1) then
           quick_xp(i,j,-1) = quick_xp(i,j,-1) + quick_xp(i,j, 0)
           quick_xp(i,j, 0) = 0.0_WP
           quick_xm(i,j,-1) = 1.0_WP
           quick_xm(i,j, 0) = 0.0_WP
           quick_xm(i,j,+1) = 0.0_WP
        end if
        if (mask(i+1,j).eq.1) then
           quick_xm(i,j, 0) = quick_xm(i,j, 0) + quick_xm(i,j,+1)
           quick_xm(i,j,+1) = 0.0_WP
        end if
     end do
  end do
  !$OMP END DO
  
  !$OMP DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        
        call poly3_coeff(y(j-2:j+1),y(j),quick_yp(i,j,:))
        call poly3_coeff(y(j-1:j+2),y(j),quick_ym(i,j,:))
        
        if (mask(i,j-2).eq.1) then
           quick_yp(i,j,-1) = quick_yp(i,j,-1) + quick_yp(i,j,-2)
           quick_yp(i,j,-2) = 0.0_WP
        end if
        if (mask(i,j-1).eq.1) then
           quick_yp(i,j,-2) = 0.0_WP
           quick_yp(i,j,-1) = 0.0_WP
           quick_yp(i,j, 0) = 1.0_WP
           quick_ym(i,j, 0) = quick_ym(i,j, 0) + quick_ym(i,j,-1)
           quick_ym(i,j,-1) = 0.0_WP
        end if
        if (mask(i,j).eq.1) then
           quick_yp(i,j,-1) = quick_yp(i,j,-1) + quick_yp(i,j, 0)
           quick_yp(i,j, 0) = 0.0_WP
           quick_ym(i,j,-1) = 1.0_WP
           quick_ym(i,j, 0) = 0.0_WP
           quick_ym(i,j,+1) = 0.0_WP
        end if
        if (mask(i,j+1).eq.1) then
           quick_ym(i,j, 0) = quick_ym(i,j, 0) + quick_ym(i,j,+1)
           quick_ym(i,j,+1) = 0.0_WP
        end if
     end do
  end do
  !$OMP END DO
  
  !$OMP DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        quick_zp(i,j,-2) = -1.0_WP/6.0_WP
        quick_zp(i,j,-1) =  5.0_WP/6.0_WP
        quick_zp(i,j, 0) =  2.0_WP/6.0_WP
        quick_zm(i,j,-1) =  2.0_WP/6.0_WP
        quick_zm(i,j, 0) =  5.0_WP/6.0_WP
        quick_zm(i,j,+1) = -1.0_WP/6.0_WP
     end do
  end do
  !$OMP END DO
  
  ! Dirichlet conditions
  ! --------------------
  !$OMP DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        if (mask_u(i,j).eq.3) then
           quick_xp(i,j,:)  = 0.0_WP
           quick_xp(i,j,-1) = 1.0_WP
           quick_xm(i,j,:)  = 0.0_WP
           quick_xm(i,j,-1) = 1.0_WP
        end if
     end do
  end do
  !$OMP END DO

  ! Boundary Conditions
  ! -------------------
  
  ! In x
  ! - Periodic
  !   -> nothing to be done
  ! - Not Periodic
  !   -> left  : Dirichlet : nothing to be done
  !   -> right : Convective outflow: nothing to be done
  
  ! In y
  ! - Periodic
  !   -> nothing to be done
  ! - Not Periodic
  !   -> up/down : Neumann
  if (yper.ne.1 .and. jproc.eq.1 .and. icyl.eq.0) then
     !$OMP DO
     do i=imin_-st1,imax_+st2
        quick_yp(i,jmin,  -2) = 0.0_WP
        quick_yp(i,jmin,  -1) = 0.0_WP
        quick_yp(i,jmin,   0) = 1.0_WP
        quick_ym(i,jmin,   0) = quick_ym(i,jmin,0) + quick_ym(i,jmin,-1)
        quick_ym(i,jmin,  -1) = 0.0_WP
        quick_yp(i,jmin+1,-1) = quick_yp(i,jmin+1,-1) + quick_yp(i,jmin+1,-2)
        quick_yp(i,jmin+1,-2) = 0.0_WP
     end do
     !$OMP END DO
  end if
  if (yper.ne.1 .and. jproc.eq.npy) then
     !$OMP DO
     do i=imin_-st1,imax_+st2
        quick_ym(i,jmax+1,-1) = 1.0_WP
        quick_ym(i,jmax+1, 0) = 0.0_WP
        quick_ym(i,jmax+1,+1) = 0.0_WP
        quick_yp(i,jmax+1,-1) = quick_yp(i,jmax+1,-1) + quick_yp(i,jmax+1,0)
        quick_yp(i,jmax+1, 0) = 0.0_WP
        quick_ym(i,jmax,   0) = quick_ym(i,jmax,0) + quick_ym(i,jmax,+1)
        quick_ym(i,jmax,  +1) = 0.0_WP
     end do
     !$OMP END DO
  end if
  
  ! In z
  ! - Periodic
  !   -> nothing to be done

  !$OMP END PARALLEL
  
  return
end subroutine scalar_quick_init


! =========================================================== !
! Compute the residuals of the scalar equations               !
!                                                             !
! - velocity field n+1 stored in U/rhoU                       !
! - scalar field n+1 stored in SCmid                          !
!                                                             !
! 3 working arrays of size (at least) (nx_+1)*(ny_+1)*(nz_+1) !
! =========================================================== !
subroutine scalar_quick_residual
  use scalar_quick
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: rhs

  do isc=1,nscalar
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
     
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              
              rhs =+sum(div_u(i,j,:)*FX(i-st1:i+st2,j,k)) &
                   +sum(div_v(i,j,:)*FY(i,j-st1:j+st2,k)) &
                   +sum(div_w(i,j,:)*FZ(i,j,k-st1:k+st2)) 

              ResSC(i,j,k,isc) =  -2.0_WP*SC(i,j,k,isc)+SCold(i,j,k,isc) &
                   + ( RHO_(i,j,k)*SC_(i,j,k,isc) + dt_*rhs  &
                   + srcSCmid(i,j,k,isc) + srcSCfull(i,j,k,isc) ) /RHO_s(i,j,k)
              
           end do
        end do
     end do
  end do
  
  return
end subroutine scalar_quick_residual


! =========================================================== !
! Inverse the linear system obtained from the implicit scalar !
! transport equation                                          !
! =========================================================== !
subroutine scalar_quick_inverse
  use scalar_quick
  use parallel
  use memory
  use implicit
  use time_info
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: conv1,conv2,conv3,conv4
  real(WP) :: visc1,visc2
  real(WP) :: dt2
  
  dt2 = dt_/2.0_WP
  
  ! If purely explicit return
  if ((.not.implicit_any)) return
  
  do isc=1,nscalar
     
     ! X-direction
     if (implicit_x) then
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoUt(i  ,j,k,isc) + abs(rhoUt(i  ,j,k,isc)))
                 conv2 = 0.5_WP*(rhoUt(i  ,j,k,isc) - abs(rhoUt(i  ,j,k,isc)))
                 conv3 = 0.5_WP*(rhoUt(i+1,j,k,isc) + abs(rhoUt(i+1,j,k,isc)))
                 conv4 = 0.5_WP*(rhoUt(i+1,j,k,isc) - abs(rhoUt(i+1,j,k,isc)))
                 
                 visc1 = sum(interp_sc_x(i  ,j,:)*DIFF(i  -st2:i  +st1,j,k,isc))
                 visc2 = sum(interp_sc_x(i+1,j,:)*DIFF(i+1-st2:i+1+st1,j,k,isc))
                 
                 Ax(j,k,i,-2) = + dt2 * div_u(i,j,0)*conv1*quick_xp(i,j,-2)
                 
                 Ax(j,k,i,-1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*quick_xp(i  ,j,-1) - conv2*quick_xm(i  ,j,-1) + visc1*grad_x(i,j,-1)) &
                      + div_u(i,j,1)*( - conv3*quick_xp(i+1,j,-2)))
                 
                 Ax(j,k,i, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*quick_xp(i  ,j, 0) - conv2*quick_xm(i  ,j, 0) + visc1*grad_x(i,j,0)) &
                      + div_u(i,j,1)*( - conv3*quick_xp(i+1,j,-1) - conv4*quick_xm(i+1,j,-1) + visc2*grad_x(i+1,j,-1)))
                 
                 Ax(j,k,i,+1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv2*quick_xm(i  ,j,+1)) &
                      + div_u(i,j,1)*( - conv3*quick_xp(i+1,j, 0) - conv4*quick_xm(i+1,j, 0) + visc2*grad_x(i+1,j, 0)))
                 
                 Ax(j,k,i,+2) = + dt2 * div_u(i,j,1)*conv4*quick_xm(i+1,j,+1)
                 
                 ! JFM 6/13/14
                 ! Use RHO_s, since that's the density at the time to which we're updating
                 Rx(j,k,i) = RHO_s(i,j,k) * ResSC(i,j,k,isc)
                 
              end do
           end do
        end do
        call implicit_solve_x(5)
     else
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 Rx(j,k,i) = ResSC(i,j,k,isc)
              end do
           end do
        end do
     end if
     
     ! Y-direction
     if (implicit_y) then
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoVt(i,j  ,k,isc) + abs(rhoVt(i,j  ,k,isc)))
                 conv2 = 0.5_WP*(rhoVt(i,j  ,k,isc) - abs(rhoVt(i,j  ,k,isc)))
                 conv3 = 0.5_WP*(rhoVt(i,j+1,k,isc) + abs(rhoVt(i,j+1,k,isc)))
                 conv4 = 0.5_WP*(rhoVt(i,j+1,k,isc) - abs(rhoVt(i,j+1,k,isc)))
                 
                 visc1 = sum(interp_sc_y(i,j  ,:)*DIFF(i,j  -st2:j  +st1,k,isc))
                 visc2 = sum(interp_sc_y(i,j+1,:)*DIFF(i,j+1-st2:j+1+st1,k,isc))
                 
                 Ay(i,k,j,-2) = + dt2 * div_v(i,j,0)*conv1*quick_yp(i,j,-2)
                 
                 Ay(i,k,j,-1) = - dt2*( &
                      + div_v(i,j,0)*( - conv1*quick_yp(i,j  ,-1) - conv2*quick_ym(i,j  ,-1) + visc1*grad_y(i,j,-1)) &
                      + div_v(i,j,1)*( - conv3*quick_yp(i,j+1,-2)))
                           
                 Ay(i,k,j, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_v(i,j,0)*( - conv1*quick_yp(i,j  , 0) - conv2*quick_ym(i,j  , 0) + visc1*grad_y(i,j,0)) &
                      + div_v(i,j,1)*( - conv3*quick_yp(i,j+1,-1) - conv4*quick_ym(i,j+1,-1) + visc2*grad_y(i,j+1,-1)))
                 
                 Ay(i,k,j,+1) = - dt2*( &
                      + div_v(i,j,0)*( - conv2*quick_ym(i,j  ,+1)) &
                      + div_v(i,j,1)*( - conv3*quick_yp(i,j+1, 0) - conv4*quick_ym(i,j+1, 0) + visc2*grad_y(i,j+1,0)))
                 
                 Ay(i,k,j,+2) = dt2 * div_v(i,j,1)*conv4*quick_ym(i,j+1,+1)
                 
                 Ry(i,k,j) = RHO_s(i,j,k) * Rx(j,k,i)
                 
              end do
           end do
        end do
        if (icyl.eq.1 .and. jproc.eq.1) then
           do k=kmin_,kmax_
              Ay(:,k,jmin+1,-1) = Ay(:,k,jmin+1,-1) + Ay(:,k,jmin+1,-2)
              Ay(:,k,jmin+1,-2) = 0.0_WP
              Ay(:,k,jmin, 0) = Ay(:,k,jmin, 0) + Ay(:,k,jmin,-1) + Ay(:,k,jmin,-2)
              Ay(:,k,jmin,-1) = 0.0_WP
              Ay(:,k,jmin,-2) = 0.0_WP
           end do
        end if
        call implicit_solve_y(5)
     else 
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 Ry(i,k,j) = Rx(j,k,i)
              end do
           end do
        end do
     end if
     
     ! Z-direction
     if (implicit_z) then
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 conv1 = 0.5_WP*(rhoWt(i,j,k  ,isc) + abs(rhoWt(i,j,k  ,isc)))
                 conv2 = 0.5_WP*(rhoWt(i,j,k  ,isc) - abs(rhoWt(i,j,k  ,isc)))
                 conv3 = 0.5_WP*(rhoWt(i,j,k+1,isc) + abs(rhoWt(i,j,k+1,isc)))
                 conv4 = 0.5_WP*(rhoWt(i,j,k+1,isc) - abs(rhoWt(i,j,k+1,isc)))
                 
                 visc1 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k  -st2:k  +st1,isc))
                 visc2 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k+1-st2:k+1+st1,isc))
                 
                 Az(i,j,k,-2) = + dt2 * div_w(i,j,0)*conv1*quick_zp(i,j, -2)
                 
                 Az(i,j,k,-1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*quick_zp(i,j,-1) - conv2*quick_zm(i,j,-1) + visc1*grad_z(i,j,-1)) &
                      + div_w(i,j,1)*( - conv3*quick_zp(i,j,-2)))
                 
                 Az(i,j,k, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*quick_zp(i,j, 0) - conv2*quick_zm(i,j, 0) + visc1*grad_z(i,j,0)) &
                      + div_w(i,j,1)*( - conv3*quick_zp(i,j,-1) - conv4*quick_zm(i,j,-1) + visc2*grad_z(i,j,-1)))
                 
                 Az(i,j,k,+1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv2*quick_zm(i,j,+1)) &
                      + div_w(i,j,1)*( - conv3*quick_zp(i,j, 0) - conv4*quick_zm(i,j, 0) + visc2*grad_z(i,j,0)))
                 
                 Az(i,j,k,+2) = + dt2 * div_w(i,j,1)*conv4*quick_zm(i,j,+1)
                 
                 Rz(i,j,k) = RHO_s(i,j,k) * Ry(i,k,j)
                 
              end do
           end do
        end do
        call implicit_solve_z(5)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 ResSC(i,j,k,isc) = Rz(i,j,k)
              end do
           end do
        end do
     else
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 ResSC(i,j,k,isc) = Ry(i,k,j)
              end do
           end do
        end do
     end if
     
  end do
  
  return
end subroutine scalar_quick_inverse
