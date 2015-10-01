module scalar_houc
  use scalar
  use data
  use metric_generic
  implicit none
  
  ! Order
  integer :: nhouc,sthp1,sthp2,sthm1,sthm2
  
  ! HOUC coefficients
  real(WP), dimension(:,:,:), pointer :: houc_xp,houc_yp
  real(WP), dimension(:,:,:), pointer :: houc_xm,houc_ym
  real(WP), dimension(:),     pointer :: houc_zp,houc_zm
  
end module scalar_houc


! ========================== !
! Initialize the HOUC module !
! ========================== !
subroutine scalar_houc_init
  use scalar_houc
  use masks
  use math
  use parser
  use string
  implicit none
  
  integer  :: i,j,st,n
  character(len=str_medium) :: scheme
  
  ! Order of HOUC
  call parser_read('Scalar scheme',scheme)
  scheme=trim(scheme(5:))
  read(scheme,'(i10)') nhouc
  if (mod(nhouc,2).eq.0) call die('HOUC must be of odd order')
  sthp1=(nhouc+1)/2
  sthp2=nhouc-sthp1-1
  sthm1=(nhouc-1)/2
  sthm2=nhouc-sthm1-1
  
  ! Allocate the metrics
  allocate(houc_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-sthp1:+sthp2))
  allocate(houc_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-sthp1:+sthp2))
  allocate(houc_zp(-sthp1:+sthp2))
  allocate(houc_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-sthm1:+sthm2))
  allocate(houc_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-sthm1:+sthm2))
  allocate(houc_zm(-sthm1:+sthm2))
  
  ! Computation of the metrics in Z
  call hofvi(nhouc,z(kmin-sthp1:kmin+sthp2+1),z(kmin),houc_zp)
  call hofvi(nhouc,z(kmin-sthm1:kmin+sthm2+1),z(kmin),houc_zm)
  
  ! Computation of the metrics in X and Y
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        call hofvi(nhouc,x(i-sthp1:i+sthp2+1),x(i),houc_xp(i,j,:))
        call hofvi(nhouc,x(i-sthm1:i+sthm2+1),x(i),houc_xm(i,j,:))
        call hofvi(nhouc,y(j-sthp1:j+sthp2+1),y(j),houc_yp(i,j,:))
        call hofvi(nhouc,y(j-sthm1:j+sthm2+1),y(j),houc_ym(i,j,:))
     end do
  end do
  
  ! Boundary Conditions
  ! -------------------
  
  ! ======================================================================== !
  ! Definition of the stencil                                                !
  !                                                                          !
  !                             SCface                                       !
  !                               |                                          !
  !   pppp      pppp      pppp    |   pppp      pppp                         !
  !                               |                                          !
  !             mmmm      mmmm    |   mmmm      mmmm     mmmmm               !
  !                               |                                          !
  !    i-3       i-2       i-1    |    i         i+1      i+2                !
  ! [point 1] [point 2] [point 3] | [point 4] [point 5]             case U>0 !
  !           [point 6] [point 7] | [point 8] [point 9] [point 10]  case U<0 !
  !                               |                                          !
  ! -----------------------------------------------------------------> x,y,z !
  ! ======================================================================== !
  
  ! Walls
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        ! x+
        do st=-sthp1,-1
           if (mask(i+st,j).eq.1) then
              houc_xp(i,j,st+1) = houc_xp(i,j,st+1) + houc_xp(i,j,st)
              houc_xp(i,j,st)   = 0.0_WP
           end if
        end do
        do st=+sthp2,0,-1
           if (mask(i+st,j).eq.1) then
              houc_xp(i,j,st-1) = houc_xp(i,j,st-1) + houc_xp(i,j,st)
              houc_xp(i,j,st)   = 0.0_WP
           end if
        end do
        ! x-
        do st=-sthm1,-1
           if (mask(i+st,j).eq.1) then
              houc_xm(i,j,st+1) = houc_xm(i,j,st+1) + houc_xm(i,j,st)
              houc_xm(i,j,st)   = 0.0_WP
           end if
        end do
        do st=+sthm2,0,-1
           if (mask(i+st,j).eq.1) then
              houc_xm(i,j,st-1) = houc_xm(i,j,st-1) + houc_xm(i,j,st)
              houc_xm(i,j,st)   = 0.0_WP
           end if
        end do
        ! y+
        do st=-sthp1,-1
           if (mask(i,j+st).eq.1) then
              houc_yp(i,j,st+1) = houc_yp(i,j,st+1) + houc_yp(i,j,st)
              houc_yp(i,j,st)   = 0.0_WP
           end if
        end do
        do st=+sthp2,0,-1
           if (mask(i,j+st).eq.1) then
              houc_yp(i,j,st-1) = houc_yp(i,j,st-1) + houc_yp(i,j,st)
              houc_yp(i,j,st)   = 0.0_WP
           end if
        end do
        ! y-
        do st=-sthm1,-1
           if (mask(i,j+st).eq.1) then
              houc_ym(i,j,st+1) = houc_ym(i,j,st+1) + houc_ym(i,j,st)
              houc_ym(i,j,st)   = 0.0_WP
           end if
        end do
        do st=+sthm2,0,-1
           if (mask(i,j+st).eq.1) then
              houc_ym(i,j,st-1) = houc_ym(i,j,st-1) + houc_ym(i,j,st)
              houc_ym(i,j,st)   = 0.0_WP
           end if
        end do
     end do
  end do
  
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
     
     ! y+
     do n=0,sthp1-1
        do st=-sthp1,-n-1
           houc_yp(:,jmin+n,-n) = houc_yp(:,jmin+n,-n) + houc_yp(:,jmin+n,st)
           houc_yp(:,jmin+n,st) = 0.0_WP
        end do
     end do
     ! y-
     do n=0,sthm1-1
        do st=-sthm1,-n-1
           houc_ym(:,jmin+n,-n) = houc_ym(:,jmin+n,-n) + houc_ym(:,jmin+n,st)
           houc_ym(:,jmin+n,st) = 0.0_WP
        end do
     end do
     
  end if
  
  if (yper.ne.1 .and. jproc.eq.npy) then
     
     ! y+
     do n=0,-sthp2+1,-1
        do st=-n,sthp2
           houc_yp(:,jmax+1+n,-n) = houc_yp(:,jmax+1+n,-n) + houc_yp(:,jmax+1+n,st)
           houc_yp(:,jmax+1+n,st) = 0.0_WP
        end do
     end do
     ! y-
     do n=0,-sthm2+1,-1
        do st=-n,sthm2
           houc_ym(:,jmax+1+n,-n) = houc_ym(:,jmax+1+n,-n) + houc_ym(:,jmax+1+n,st)
           houc_ym(:,jmax+1+n,st) = 0.0_WP
        end do
     end do
     
  end if
  
  ! In z
  ! - Periodic
  !   -> nothing to be done
  
  return
end subroutine scalar_houc_init


! =========================================================== !
! Compute the residuals of the scalar equations               !
!                                                             !
! - velocity field n+1 stored in U/rhoU                       !
! - scalar field n+1 stored in SCmid                          !
!                                                             !
! 3 working arrays of size (at least) (nx_+1)*(ny_+1)*(nz_+1) !
! =========================================================== !
subroutine scalar_houc_residual
  use scalar_houc
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: rhs
  
  do isc=1,nscalar
     
     ! Now compute residuals
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              
              FX(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoUt(i,j,k,isc)+abs(rhoUt(i,j,k,isc)))*sum(houc_xp(i,j,:)*SC(i-sthp1:i+sthp2,j,k,isc)) &
                   - 0.5_WP*(rhoUt(i,j,k,isc)-abs(rhoUt(i,j,k,isc)))*sum(houc_xm(i,j,:)*SC(i-sthm1:i+sthm2,j,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_x(i,j,:)*DIFF(i-st2:i+st1,j,k,isc)) * &
                     sum(grad_x(i,j,:)*SC(i-st2:i+st1,j,k,isc))
              
              FY(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoVt(i,j,k,isc)+abs(rhoVt(i,j,k,isc)))*sum(houc_yp(i,j,:)*SC(i,j-sthp1:j+sthp2,k,isc)) &
                   - 0.5_WP*(rhoVt(i,j,k,isc)-abs(rhoVt(i,j,k,isc)))*sum(houc_ym(i,j,:)*SC(i,j-sthm1:j+sthm2,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_y(i,j,:)*DIFF(i,j-st2:j+st1,k,isc)) * &
                     sum(grad_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
              
              FZ(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoWt(i,j,k,isc)+abs(rhoWt(i,j,k,isc)))*sum(houc_zp(:)    *SC(i,j,k-sthp1:k+sthp2,isc)) &
                   - 0.5_WP*(rhoWt(i,j,k,isc)-abs(rhoWt(i,j,k,isc)))*sum(houc_zm(:)    *SC(i,j,k-sthm1:k+sthm2,isc)) &
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
                   + ( RHO_(i,j,k)*SC_(i,j,k,isc) + dt_*rhs &
                   + srcSCmid(i,j,k,isc) +srcSCfull(i,j,k,isc) ) /RHO_s(i,j,k)
              
           end do
        end do
     end do
     
  end do
  
  return
end subroutine scalar_houc_residual


! =========================================================== !
! Inverse the linear system obtained from the implicit scalar !
! transport equation                                          !
! =========================================================== !
subroutine scalar_houc_inverse
  use scalar_houc
  use parallel
  use memory
  use implicit
  use time_info
  implicit none
  
  integer  :: i,j,k,isc,st,n
  real(WP) :: visc1,visc2
  real(WP) :: dt2
  
  dt2 = dt/2.0_WP

  ! If purely explicit return
  if ((.not.implicit_any) .or. (strang_splitting)) return
  
  do isc=1,nscalar
     
     ! X-direction
     if (implicit_x) then
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 visc1 = sum(interp_sc_x(i  ,j,:)*DIFF(i  -st2:i  +st1,j,k,isc))
                 visc2 = sum(interp_sc_x(i+1,j,:)*DIFF(i+1-st2:i+1+st1,j,k,isc))
                 
                 Ax(j,k,i,:) = 0.0_WP
                 
                 Ax(j,k,i,-1) = - dt2 * div_u(i,j,0)*visc1*grad_x(i,j,-1)
                 
                 Ax(j,k,i, 0) = RHO(i,j,k) - dt2*(div_u(i,j,0)*visc1*grad_x(i,j,0) + div_u(i,j,1)*visc2*grad_x(i+1,j,-1))
                 
                 Ax(j,k,i,+1) = - dt2 * div_u(i,j,1)*visc2*grad_x(i+1,j,0)
                 
                 Rx(j,k,i) = RHO(i,j,k) * ResSC(i,j,k,isc)
                 
              end do
           end do
        end do
        
        if (npx.eq.1 .or. nhouc.le.3) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    do st=-st1,st2
                       do n=-sthp1,sthp2
                          Ax(j,k,i,st+n) = Ax(j,k,i,st+n)+dt2*div_u(i,j,st)*0.5_WP*(rhoUt(i+st,j,k,isc)+abs(rhoUt(i+st,j,k,isc)))*houc_xp(i+st,j,n)
                       end do
                       do n=-sthm1,sthm2
                          Ax(j,k,i,st+n) = Ax(j,k,i,st+n)+dt2*div_u(i,j,st)*0.5_WP*(rhoUt(i+st,j,k,isc)-abs(rhoUt(i+st,j,k,isc)))*houc_xm(i+st,j,n)
                       end do
                    end do
                 end do
              end do
           end do
           call implicit_solve_x(nhouc+2)
        else
           call implicit_solve_x(3)
        end if
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
                 
                 visc1 = sum(interp_sc_y(i,j  ,:)*DIFF(i,j  -st2:j  +st1,k,isc))
                 visc2 = sum(interp_sc_y(i,j+1,:)*DIFF(i,j+1-st2:j+1+st1,k,isc))
                 
                 Ay(i,k,j,:) = 0.0_WP

                 Ay(i,k,j,-1) = - dt2 * div_v(i,j,0)*visc1*grad_y(i,j,-1)
                           
                 Ay(i,k,j, 0) = RHO(i,j,k) - dt2*(div_v(i,j,0)*visc1*grad_y(i,j,0) + div_v(i,j,1)*visc2*grad_y(i,j+1,-1))
                 
                 Ay(i,k,j,+1) = - dt2 * div_v(i,j,1)*visc2*grad_y(i,j+1,0)
                 
                 Ry(i,k,j) = RHO(i,j,k) * Rx(j,k,i)
                 
              end do
           end do
        end do

        if (npy.eq.1 .or. nhouc.le.3) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    do st=-st1,st2
                       do n=-sthp1,sthp2
                          Ay(i,k,j,st+n) = Ay(i,k,j,st+n)+dt2*div_v(i,j,st)*0.5_WP*(rhoVt(i,j+st,k,isc)+abs(rhoVt(i,j+st,k,isc)))*houc_yp(i,j+st,n)
                       end do
                       do n=-sthm1,sthm2
                          Ay(i,k,j,st+n) = Ay(i,k,j,st+n)+dt2*div_v(i,j,st)*0.5_WP*(rhoVt(i,j+st,k,isc)-abs(rhoVt(i,j+st,k,isc)))*houc_ym(i,j+st,n)
                       end do
                    end do
                 end do
              end do
           end do
           if (icyl.eq.1 .and. jproc.eq.1) then
              do k=kmin_,kmax_
                 do n=0,sthp1-1
                    do st=-sthp1,-n-1
                       Ay(:,k,jmin+n,-n) = Ay(:,k,jmin+n,-n) + Ay(:,k,jmin+n,st)
                       Ay(:,k,jmin+n,st) = 0.0_WP
                    end do
                 end do
              end do
           end if
           call implicit_solve_y(nhouc+2)
        else
           if (icyl.eq.1 .and. jproc.eq.1) then
              do k=kmin_,kmax_
                 Ay(:,k,jmin, 0) = Ay(:,k,jmin, 0) + Ay(:,k,jmin,-1)
                 Ay(:,k,jmin,-1) = 0.0_WP
              end do
           end if
           call implicit_solve_y(3)
        end if
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
                 
                 visc1 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k  -st2:k  +st1,isc))
                 visc2 = sum(interp_sc_z(i,j,:)*DIFF(i,j,k+1-st2:k+1+st1,isc))
                 
                 Az(i,j,k,:) = 0.0_WP

                 Az(i,j,k,-1) = - dt2 * div_w(i,j,0)*visc1*grad_z(i,j,-1)
                 
                 Az(i,j,k, 0) = RHO(i,j,k) - dt2*(div_w(i,j,0)*visc1*grad_z(i,j,0) + div_w(i,j,1)*visc2*grad_z(i,j,-1))
                 
                 Az(i,j,k,+1) = - dt2 * div_w(i,j,1)*visc2*grad_z(i,j,0)
                 
                 Rz(i,j,k) = RHO(i,j,k) * Ry(i,k,j)
                 
              end do
           end do
        end do
        
        if (npz.eq.1 .or. nhouc.le.3) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    do st=-st1,st2
                       do n=-sthp1,sthp2
                          Az(i,j,k,st+n) = Az(i,j,k,st+n)+dt2*div_w(i,j,st)*0.5_WP*(rhoWt(i,j,k+st,isc)+abs(rhoWt(i,j,k+st,isc)))*houc_zp(n)
                       end do
                       do n=-sthm1,sthm2
                          Az(i,j,k,st+n) = Az(i,j,k,st+n)+dt2*div_w(i,j,st)*0.5_WP*(rhoWt(i,j,k+st,isc)-abs(rhoWt(i,j,k+st,isc)))*houc_zm(n)
                       end do
                    end do
                 end do
              end do
           end do
           call implicit_solve_z(nhouc+2)
        else
           call implicit_solve_z(3)
        end if
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
end subroutine scalar_houc_inverse
