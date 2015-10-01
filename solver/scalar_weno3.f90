module scalar_weno3
  use scalar
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
  
  ! Coefficient for the weno 3rd
  real(WP), dimension(:,:,:,:), pointer :: weno3_xp,weno3_yp,weno3_zp
  real(WP), dimension(:,:,:,:), pointer :: weno3_xm,weno3_ym,weno3_zm

  ! Coefficient for the
  ! -> CD : Central difference
  ! -> UP : Upwind biased
  real(WP), dimension(:,:,:), pointer :: cd_xp,cd_yp,cd_xm,cd_ym
  real(WP), dimension(:,:,:), pointer :: up_xp,up_yp,up_xm,up_ym
  real(WP), dimension(:),     pointer :: cd_zp,cd_zm,up_zp,up_zm
  
end module scalar_weno3


! ============================ !
! Initialize the weno 3 module !
! ============================ !
subroutine scalar_weno3_init
  use scalar_weno3
  use masks
  implicit none
  
  integer :: i,j
  
  ! Allocate arrays
  allocate(weno3_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2: 0))
  allocate(weno3_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2: 0))
  allocate(weno3_zp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2: 0))
  allocate(weno3_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-1:+1))
  allocate(weno3_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-1:+1))
  allocate(weno3_zm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-1:+1))
  
  allocate(cd_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2: 0))
  allocate(cd_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2: 0))
  allocate(cd_zp(-2: 0))
  allocate(cd_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(cd_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(cd_zm(-1:+1))

  allocate(up_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2: 0))
  allocate(up_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2: 0))
  allocate(up_zp(-2: 0))
  allocate(up_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(up_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-1:+1))
  allocate(up_zm(-1:+1))
  
  ! Interpolation at the faces
  !$OMP PARALLEL DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        cd_xp(i,j,:) = 0.0_WP; cd_xm(i,j,:) = 0.0_WP
        cd_yp(i,j,:) = 0.0_WP; cd_ym(i,j,:) = 0.0_WP
     end do
  end do
  !$OMP END PARALLEL DO
  cd_zp = 0.0_WP; cd_zm = 0.0_WP
  cd_zp(-1) = 0.5_WP
  cd_zp( 0) = 0.5_WP
  cd_zm(-1) = 0.5_WP
  cd_zm( 0) = 0.5_WP
  
  !$OMP PARALLEL DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        up_xp(i,j,:) = 0.0_WP; up_xm(i,j,:) = 0.0_WP
        up_yp(i,j,:) = 0.0_WP; up_ym(i,j,:) = 0.0_WP
     end do
  end do
  !$OMP END PARALLEL DO
  up_zp = 0.0_WP; up_zm = 0.0_WP
  up_zp(-2) = -0.5_WP
  up_zp(-1) =  1.5_WP
  up_zm( 0) =  1.5_WP
  up_zm(+1) = -0.5_WP
  
  !$OMP PARALLEL

  !$OMP DO
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        
        ! CENTRAL DIFFERENCE Scheme
        ! Along x
        cd_xp(i,j,-1) = (xm(i)-x (i  ))/(xm(i)-xm(i-1))
        cd_xp(i,j, 0) = (x (i)-xm(i-1))/(xm(i)-xm(i-1))
        cd_xm(i,j,-1) = (xm(i)-x (i  ))/(xm(i)-xm(i-1))
        cd_xm(i,j, 0) = (x (i)-xm(i-1))/(xm(i)-xm(i-1))
        ! Along y
        cd_yp(i,j,-1) = (ym(j)-y (j  ))/(ym(j)-ym(j-1))
        cd_yp(i,j, 0) = (y (j)-ym(j-1))/(ym(j)-ym(j-1))
        cd_ym(i,j,-1) = (ym(j)-y (j  ))/(ym(j)-ym(j-1))
        cd_ym(i,j, 0) = (y (j)-ym(j-1))/(ym(j)-ym(j-1))
        
        ! UPWIND Scheme
        ! Along x
        up_xp(i,j,-2) = (xm(i-1)-x (i  ))/(xm(i-1)-xm(i-2))
        up_xp(i,j,-1) = (x (i  )-xm(i-2))/(xm(i-1)-xm(i-2))
        up_xm(i,j, 0) = (xm(i+1)-x (i  ))/(xm(i+1)-xm(i  ))
        up_xm(i,j,+1) = (x (i  )-xm(i  ))/(xm(i+1)-xm(i  ))
        ! Along y
        up_yp(i,j,-2) = (ym(j-1)-y (j  ))/(ym(j-1)-ym(j-2))
        up_yp(i,j,-1) = (y (j  )-ym(j-2))/(ym(j-1)-ym(j-2))
        up_ym(i,j, 0) = (ym(j+1)-y (j  ))/(ym(j+1)-ym(j  ))
        up_ym(i,j,+1) = (y (j  )-ym(j  ))/(ym(j+1)-ym(j  ))
        
        ! Wall BC
        if (mask(i-2,j).eq.1) then
           cd_xp(i,j,: ) = 0.0_WP
           cd_xp(i,j,-1) = 1.0_WP
           up_xp(i,j,: ) = 0.0_WP
           up_xp(i,j,-1) = 1.0_WP
        end if
        if (mask(i-1,j).eq.1) then
           cd_xp(i,j,: ) = 0.0_WP
           cd_xp(i,j, 0) = 1.0_WP
           up_xp(i,j,: ) = 0.0_WP
           up_xp(i,j, 0) = 1.0_WP
        end if
        if (mask(i,j).eq.1) then
           cd_xm(i,j,: ) = 0.0_WP
           cd_xm(i,j,-1) = 1.0_WP
           up_xm(i,j,: ) = 0.0_WP
           up_xm(i,j,-1) = 1.0_WP
        end if
        if (mask(i+1,j).eq.1) then
           cd_xm(i,j,: ) = 0.0_WP
           cd_xm(i,j, 0) = 1.0_WP
           up_xm(i,j,: ) = 0.0_WP
           up_xm(i,j, 0) = 1.0_WP
        end if
        if (mask(i,j-2).eq.1) then
           cd_yp(i,j,: ) = 0.0_WP
           cd_yp(i,j,-1) = 1.0_WP
           up_yp(i,j,: ) = 0.0_WP
           up_yp(i,j,-1) = 1.0_WP
        end if
        if (mask(i,j-1).eq.1) then
           cd_yp(i,j,: ) = 0.0_WP
           cd_yp(i,j, 0) = 1.0_WP
           up_yp(i,j,: ) = 0.0_WP
           up_yp(i,j, 0) = 1.0_WP
        end if
        if (mask(i,j).eq.1) then
           cd_ym(i,j,: ) = 0.0_WP
           cd_ym(i,j,-1) = 1.0_WP
           up_ym(i,j,: ) = 0.0_WP
           up_ym(i,j,-1) = 1.0_WP
        end if
        if (mask(i,j+1).eq.1) then
           cd_ym(i,j,: ) = 0.0_WP
           cd_ym(i,j, 0) = 1.0_WP
           up_ym(i,j,: ) = 0.0_WP
           up_ym(i,j, 0) = 1.0_WP
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
        cd_yp(i,jmin,:) = 0.0_WP
        cd_yp(i,jmin,0) = 1.0_WP
        up_yp(i,jmin,:) = 0.0_WP
        up_yp(i,jmin,0) = 1.0_WP
        cd_ym(i,jmin,:) = 0.0_WP
        cd_ym(i,jmin,0) = 1.0_WP
        up_ym(i,jmin,:) = 0.0_WP
        up_ym(i,jmin,0) = 1.0_WP
        
        cd_yp(i,jmin+1,: ) = 0.0_WP
        cd_yp(i,jmin+1,-1) = 1.0_WP
        up_yp(i,jmin+1,: ) = 0.0_WP
        up_yp(i,jmin+1,-1) = 1.0_WP
     end do
     !$OMP END DO
  end if
  if (yper.ne.1 .and. jproc.eq.npy) then
     !$OMP DO
     do i=imin_-st1,imax_+st2
        cd_yp(i,jmax+1,: ) = 0.0_WP
        cd_yp(i,jmax+1,-1) = 1.0_WP
        up_yp(i,jmax+1,: ) = 0.0_WP
        up_yp(i,jmax+1,-1) = 1.0_WP
        cd_ym(i,jmax+1,: ) = 0.0_WP
        cd_ym(i,jmax+1,-1) = 1.0_WP
        up_ym(i,jmax+1,: ) = 0.0_WP
        up_ym(i,jmax+1,-1) = 1.0_WP
     
        cd_ym(i,jmax,:) = 0.0_WP
        cd_ym(i,jmax,0) = 1.0_WP
        up_ym(i,jmax,:) = 0.0_WP
        up_ym(i,jmax,0) = 1.0_WP
     end do
     !$OMP END DO
  end if
  
  ! In z
  ! - Periodic
  !   -> nothing to be done
  
  !$OMP END PARALLEL

  return
end subroutine scalar_weno3_init


! ==================================================== !
! Compute the weno 3 coefficients for the given scalar !
! ==================================================== !
subroutine scalar_weno3_coeff(isc)
  use scalar_weno3
  implicit none
  
  integer, intent(in) :: isc
  integer :: i,j,k
  real(WP) :: alpha,r
  real(WP), parameter :: epsilon = 1.0e-6_WP
  
  !$OMP PARALLEL DO PRIVATE(r,alpha)
  do k=kmin_-st1,kmax_+st2
     do j=jmin_-st1,jmax_+st2
        do i=imin_-st1,imax_+st2
           ! Direction x - U>0
           r = (epsilon+((SC(i-1,j,k,isc)-SC(i-2,j,k,isc))*dxmi(i-2))**2)/(epsilon+((SC(i,j,k,isc)-SC(i-1,j,k,isc))*dxmi(i-1))**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_xp(i,j,k,:) = (1.0_WP-alpha) * cd_xp(i,j,:) + alpha * up_xp(i,j,:)
           
           ! Direction x - U<0
           r = (epsilon+((SC(i+1,j,k,isc)-SC(i,j,k,isc))*dxmi(i))**2)/(epsilon+((SC(i,j,k,isc)-SC(i-1,j,k,isc))*dxmi(i-1))**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_xm(i,j,k,:) = (1.0_WP-alpha) * cd_xm(i,j,:) + alpha * up_xm(i,j,:)
           
           ! Direction y - V>0
           r = (epsilon+((SC(i,j-1,k,isc)-SC(i,j-2,k,isc))*dymi(j-2))**2)/(epsilon+((SC(i,j,k,isc)-SC(i,j-1,k,isc))*dymi(j-1))**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_yp(i,j,k,:) = (1.0_WP-alpha) * cd_yp(i,j,:) + alpha * up_yp(i,j,:)
           
           ! Direction y - V<0
           r = (epsilon+((SC(i,j+1,k,isc)-SC(i,j,k,isc))*dymi(j))**2)/(epsilon+((SC(i,j,k,isc)-SC(i,j-1,k,isc))*dymi(j-1))**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_ym(i,j,k,:) = (1.0_WP-alpha) * cd_ym(i,j,:) + alpha * up_ym(i,j,:)
           
           ! Direction z - W>0
           r = (epsilon+((SC(i,j,k-1,isc)-SC(i,j,k-2,isc))*dzi)**2)/(epsilon+((SC(i,j,k,isc)-SC(i,j,k-1,isc))*dzi)**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_zp(i,j,k,:) = (1.0_WP-alpha) * cd_zp(:) + alpha * up_zp(:)
           
           ! Direction z - W<0
           r = (epsilon+((SC(i,j,k+1,isc)-SC(i,j,k,isc))*dzi)**2)/(epsilon+((SC(i,j,k,isc)-SC(i,j,k-1,isc))*dzi)**2)
           alpha = 1.0_WP/(1.0_WP+2.0_WP*r**2)
           weno3_zm(i,j,k,:) = (1.0_WP-alpha) * cd_zm(:) + alpha * up_zm(:)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine scalar_weno3_coeff


! =========================================================== !
! Compute the residuals of the scalar equations               !
!                                                             !
! - velocity field n+1 stored in U/rhoU                       !
! - scalar field n+1 stored in SCmid                          !
!                                                             !
! 3 working arrays of size (at least) (nx_+1)*(ny_+1)*(nz_+1) !
! =========================================================== !
subroutine scalar_weno3_residual
  use scalar_weno3
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: rhs
  
  do isc=1,nscalar
     ! Get coefficients first
     call scalar_weno3_coeff(isc)
     
     !$OMP PARALLEL PRIVATE(rhs)

     ! Now compute residuals
     !$OMP DO
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              
              FX(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoUt(i,j,k,isc)+abs(rhoUt(i,j,k,isc))) * sum(weno3_xp(i,j,k,:)*SC(i-2:i  ,j,k,isc)) &
                   - 0.5_WP*(rhoUt(i,j,k,isc)-abs(rhoUt(i,j,k,isc))) * sum(weno3_xm(i,j,k,:)*SC(i-1:i+1,j,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_x(i,j,:)*DIFF(i-st2:i+st1,j,k,isc)) * &
                     sum(grad_x(i,j,:)*SC(i-st2:i+st1,j,k,isc))

              FY(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoVt(i,j,k,isc)+abs(rhoVt(i,j,k,isc))) * sum(weno3_yp(i,j,k,:)*SC(i,j-2:j  ,k,isc)) &
                   - 0.5_WP*(rhoVt(i,j,k,isc)-abs(rhoVt(i,j,k,isc))) * sum(weno3_ym(i,j,k,:)*SC(i,j-1:j+1,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_y(i,j,:)*DIFF(i,j-st2:j+st1,k,isc)) * &
                     sum(grad_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
              
              FZ(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoWt(i,j,k,isc)+abs(rhoWt(i,j,k,isc))) * sum(weno3_zp(i,j,k,:)*SC(i,j,k-2:k  ,isc)) &
                   - 0.5_WP*(rhoWt(i,j,k,isc)-abs(rhoWt(i,j,k,isc))) * sum(weno3_zm(i,j,k,:)*SC(i,j,k-1:k+1,isc)) &
                   ! Viscous term
                   + sum(interp_sc_z(i,j,:)*DIFF(i,j,k-st2:k+st1,isc)) * &
                     sum(grad_z(i,j,:)*SC(i,j,k-st2:k+st1,isc))
           end do
        end do
     end do
     !$OMP END DO
     
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              
              rhs =+sum(div_u(i,j,:)*FX(i-st1:i+st2,j,k)) &
                   +sum(div_v(i,j,:)*FY(i,j-st1:j+st2,k)) &
                   +sum(div_w(i,j,:)*FZ(i,j,k-st1:k+st2)) 
              
              ResSC(i,j,k,isc) =  -2.0_WP*SC(i,j,k,isc)+SC_(i,j,k,isc) &
                   + ( RHO_(i,j,k)*SC_(i,j,k,isc) + dt_*rhs &
                   + srcSCmid(i,j,k,isc) +srcSCfull(i,j,k,isc) ) /RHO_s(i,j,k)

           end do
        end do
     end do
     !$OMP END DO

     !$OMP END PARALLEL

  end do

  return
end subroutine scalar_weno3_residual


! =========================================================== !
! Inverse the linear system obtained from the implicit scalar !
! transport equation                                          !
! =========================================================== !
subroutine scalar_weno3_inverse
  use scalar_weno3
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
     
     ! Get coefficients first
     call scalar_weno3_coeff(isc)
     
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
                 
                 visc1 = sum(interp_sc_x(i  ,j,:)*DIFF(i  -st2:i  +st1,j,k,isc))
                 visc2 = sum(interp_sc_x(i+1,j,:)*DIFF(i+1-st2:i+1+st1,j,k,isc))
                 
                 Ax(j,k,i,-2) = + dt2 * div_u(i,j,0)*conv1*weno3_xp(i,j,k,-2)
                 
                 Ax(j,k,i,-1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*weno3_xp(i  ,j,k,-1) - conv2*weno3_xm(i  ,j,k,-1) + visc1*grad_x(i,j,-1)) &
                      + div_u(i,j,1)*( - conv3*weno3_xp(i+1,j,k,-2)))
                 
                 Ax(j,k,i, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_u(i,j,0)*( - conv1*weno3_xp(i  ,j,k, 0) - conv2*weno3_xm(i  ,j,k, 0) + visc1*grad_x(i,j,0)) &
                      + div_u(i,j,1)*( - conv3*weno3_xp(i+1,j,k,-1) - conv4*weno3_xm(i+1,j,k,-1) + visc2*grad_x(i+1,j,-1)))
                 
                 Ax(j,k,i,+1) = - dt2 * ( &
                      + div_u(i,j,0)*( - conv2*weno3_xm(i  ,j,k,+1)) &
                      + div_u(i,j,1)*( - conv3*weno3_xp(i+1,j,k, 0) - conv4*weno3_xm(i+1,j,k, 0) + visc2*grad_x(i+1,j,0)))
                 
                 Ax(j,k,i,+2) = + dt2 * div_u(i,j,1)*conv4*weno3_xm(i+1,j,k,+1)
                 
                 ! JFM 6/13/14
                 ! Use RHO_s, since that's the density at the time to which we're updating
                 Rx(j,k,i) = RHO_s(i,j,k) * ResSC(i,j,k,isc)
                 
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
                 
                 Ay(i,k,j,-2) = + dt2 * div_v(i,j,0)*conv1*weno3_yp(i,j,k,-2)
                 
                 Ay(i,k,j,-1) = - dt2*( &
                      + div_v(i,j,0)*( - conv1*weno3_yp(i,j  ,k,-1) - conv2*weno3_ym(i,j  ,k,-1) + visc1*grad_y(i,j,-1)) &
                      + div_v(i,j,1)*( - conv3*weno3_yp(i,j+1,k,-2)))
                           
                 Ay(i,k,j, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_v(i,j,0)*( - conv1*weno3_yp(i,j  ,k, 0) - conv2*weno3_ym(i,j  ,k, 0) + visc1*grad_y(i,j,0)) &
                      + div_v(i,j,1)*( - conv3*weno3_yp(i,j+1,k,-1) - conv4*weno3_ym(i,j+1,k,-1) + visc2*grad_y(i,j+1,-1)))
                 
                 Ay(i,k,j,+1) = - dt2*( &
                      + div_v(i,j,0)*( - conv2*weno3_ym(i,j  ,k,+1)) &
                      + div_v(i,j,1)*( - conv3*weno3_yp(i,j+1,k, 0) - conv4*weno3_ym(i,j+1,k, 0) + visc2*grad_y(i,j+1,0)))
                 
                 Ay(i,k,j,+2) = dt2 * div_v(i,j,1)*conv4*weno3_ym(i,j+1,k,+1)
                 
                 Ry(i,k,j) = RHO_s(i,j,k) * Rx(j,k,i)
                 
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
                 
                 Az(i,j,k,-2) = + dt2 * div_w(i,j,0)*conv1*weno3_zp(i,j,k,-2)
                 
                 Az(i,j,k,-1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*weno3_zp(i,j,k  ,-1) - conv2*weno3_zm(i,j,k,-1) + visc1*grad_z(i,j,-1)) &
                      + div_w(i,j,1)*( - conv3*weno3_zp(i,j,k+1,-2)))
                 
                 Az(i,j,k, 0) = RHO_s(i,j,k) - dt2 * ( &
                      + div_w(i,j,0)*( - conv1*weno3_zp(i,j,k  , 0) - conv2*weno3_zm(i,j,k  , 0) + visc1*grad_z(i,j,0)) &
                      + div_w(i,j,1)*( - conv3*weno3_zp(i,j,k+1,-1) - conv4*weno3_zm(i,j,k+1,-1) + visc2*grad_z(i,j,-1)))
                 
                 Az(i,j,k,+1) = - dt2 * ( &
                      + div_w(i,j,0)*( - conv2*weno3_zm(i,j,k  ,+1)) &
                      + div_w(i,j,1)*( - conv3*weno3_zp(i,j,k+1, 0) - conv4*weno3_zm(i,j,k+1, 0) + visc2*grad_z(i,j,0)))
                 
                 Az(i,j,k,+2) = + dt2 * div_w(i,j,1)*conv4*weno3_zm(i,j,k+1,+1)
                 
                 Rz(i,j,k) = RHO_s(i,j,k) * Ry(i,k,j)
                 
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
end subroutine scalar_weno3_inverse
