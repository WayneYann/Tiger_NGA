module outflow
  use data
  use parallel
  use boundary
  use borders
  implicit none
  
end module outflow


! ====================== !
! Initialize the Outflow !
! ====================== !
subroutine outflow_init
  use outflow
  use parser
  implicit none

  ! Nothing to do if:
  ! -> no outlet
  ! -> not last proc in x
  ! -> periodic in x
  if (noutlet.eq.0 .or. iproc.ne.npx .or. xper.eq.1) return

  if (compressible) call nscbc_init
  
  return     
end subroutine outflow_init


! ============================================================ !
! Convective outflow condition for velocity/momentum           !
! Update the exit plane                                        !
!                                                              !
! DOES NOT update the ghost cells                              !
! boundary_update_border REQUIRED after it (except if ntime=0) !
! ============================================================ !
subroutine outflow_velocity
  use velocity
  use outflow
  use metric_generic
  use time_info
  implicit none
  
  real(WP) :: u_a,tmp,alpha,sigma
  real(WP) :: rhox,rhoy,rhoz
  real(WP) :: newU,newV,newW
  integer :: i,j,k
  
  ! Nothing to do if:
  ! -> no outlet
  ! -> periodic in x
  if (noutlet.eq.0 .or. xper.eq.1) return

  ! Take u_a to be the maximum velocity at the exit plane
  if (iproc.eq.npx) then
     tmp = maxval(abs(U(imax,:,:)))
  else
     tmp = 0.0_WP
  end if
  call parallel_max(tmp,u_a)
  
  ! Outlow condition only if last cpu in x
  if (iproc.eq.npx) then
     if (.not.compressible) then
        if (ntime.ne.0) then
           ! Convective condition: u,t = -c u,x
           ! Implicit formulation for stability reasons
           
           !$OMP PARALLEL PRIVATE(sigma,alpha,newU,newV,newW,rhox,rhoy,rhoz)
           
           sigma = dt_uvw*u_a*dxi(imax)
           alpha = (1.0_WP-0.5_WP*sigma)/(1.0_WP+0.5_WP*sigma)
           !$OMP DO
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 newU = alpha * rhoUold(imax+1,j,k) &
                      + (1.0_WP-alpha) * 0.5_WP*(rhoU(imax,j,k)+rhoUold(imax,j,k))
                 rhox = sum(interp_sc_x(imax+1,j,:)*RHOmid(imax+1-st2:imax+1+st1,j,k))
                 do i=imax+1,imaxo
                    rhoU(i,j,k) = newU
                    U(i,j,k) = newU / rhox
                 end do
              end do
           end do
           !$OMP END DO
           
           sigma = dt_uvw*u_a*dxmi(imax)
           alpha = (1.0_WP-0.5_WP*sigma)/(1.0_WP+0.5_WP*sigma)
           !$OMP DO
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 newV = alpha * Vold(imax+1,j,k) &
                      + (1.0_WP-alpha) * 0.5_WP*(V(imax,j,k)+Vold(imax,j,k))
                 newW = alpha * Wold(imax+1,j,k) &
                      + (1.0_WP-alpha) * 0.5_WP*(W(imax,j,k)+Wold(imax,j,k))
                 rhoy = sum(interp_sc_y(imax+1,j,:)*RHOmid(imax+1,j-st2:j+st1,k))
                 rhoz = sum(interp_sc_z(imax+1,j,:)*RHOmid(imax+1,j,k-st2:k+st1))
                 do i=imax+1,imaxo
                    V(i,j,k) = newV
                    rhoV(i,j,k) = newV * rhoy
                    W(i,j,k) = newW
                    rhoW(i,j,k) = newW * rhoz
                 end do
              end do
           end do
           !$OMP END DO
           
           !$OMP END PARALLEL
        else
           ! Neumann condition
           
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imax+1,imaxo
                    U(i,j,k) = U(imax,j,k)
                    V(i,j,k) = V(imax,j,k)
                    W(i,j,k) = W(imax,j,k)
                    rhoU(i,j,k) = rhoU(imax,j,k)
                    rhoV(i,j,k) = rhoV(imax,j,k)
                    rhoW(i,j,k) = rhoW(imax,j,k)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     else
        call nscbc_outflow_velocity
     end if
  end if
  
  return
end subroutine outflow_velocity


! =========================================== !
! Convective outflow condition for scalars    !
! Update the exit plane                       !
!                                             !
! DOES update the ghost cells                 !
! NO boundary_update_border required after it !
! =========================================== !
subroutine outflow_scalar
  use scalar
  use outflow
  use time_info
  implicit none
  
  real(WP) :: u_a,tmp,alpha,sigma,newSC
  integer :: i,j,k,isc
  
  ! Nothing to do if:
  ! -> no outlet
  ! -> periodic in x
  if (noutlet.eq.0 .or. xper.eq.1) return

  ! Take u_a to be the maximum velocity at the exit plane
  if (iproc.eq.npx) then
     tmp = maxval(abs(U(imax,:,:)))
  else
     tmp = 0.0_WP
  end if
  call parallel_max(tmp,u_a)
  
  ! Outlow condition only if last cpu in x
  if (iproc.eq.npx) then
     if (.not.compressible) then
        if (ntime.ne.0) then
           ! Convective condition: u,t = -c u,x
           ! Implicit formulation for stability reasons
           
           !$OMP PARALLEL PRIVATE(sigma,alpha,newSC)
           
           sigma = dt*u_a*dxmi(imax)
           alpha = (1.0_WP-0.5_WP*sigma)/(1.0_WP+0.5_WP*sigma)
           do isc=1,nscalar
              !$OMP DO
              do k=kmin_,kmax_
                 do j=jmin_,jmax_
                    newSC = alpha * SCold(imax+1,j,k,isc) &
                         + (1.0_WP-alpha) * 0.5_WP*(SC(imax,j,k,isc)+SCold(imax,j,k,isc))
                    
                    do i=imax+1,imaxo
                       SC(i,j,k,isc) = newSC
                    end do
                 end do
              end do
              !$OMP END DO
           end do
           
           !$OMP END PARALLEL
        else
           ! Neumann condition
           
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imax+1,imaxo
                    SC(i,j,k,:) = SC(imax,j,k,:)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     else
        call nscbc_outflow_scalar
     end if
  end if
  
  return
end subroutine outflow_scalar


! ============================================================ !
! Apply a correction to the outflow to ensure mass consistency !
! Update the exit plane                                        !
!                                                              !
! DOES update the ghost cells                                  !
! NO boundary_update_border required after it                  !
! ============================================================ !
subroutine outflow_correction
  use outflow
  use velocity
  use metric_generic
  implicit none
  integer :: i,j,k,nflow
  real(WP) :: newU,rhoi,alpha

  ! Nothing to do if:
  ! -> no outlet
  ! -> not last proc in x
  ! -> periodic in x
  if (noutlet.eq.0 .or. iproc.ne.npx .or. xper.eq.1) return
  
  if (massflux_exit.eq.0.0_WP) then
     !$OMP PARALLEL DO PRIVATE(newU,rhoi)
     do k=kmino_,kmaxo_
        do nflow=1,noutlet
           do j=max(jmino_,outlet(nflow)%jmino),min(jmaxo_,outlet(nflow)%jmaxo)
              newU = rhoU(imax+1,j,k) + masscorrection/outlet_area
              rhoi = sum(interp_sc_x(imax+1,j,:)*RHOmid(imax+1-st2:imax+1+st1,j,k))
              do i=imax+1,imaxo
                 rhoU(i,j,k) = newU
                 U(i,j,k) = newU / rhoi
              end do
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     alpha = 1.0_WP + masscorrection/massflux_exit
     !$OMP PARALLEL DO PRIVATE(newU,rhoi)
     do k=kmino_,kmaxo_
        do nflow=1,noutlet
           do j=max(jmino_,outlet(nflow)%jmino),min(jmaxo_,outlet(nflow)%jmaxo)
              newU = rhoU(imax+1,j,k) * alpha
              rhoi = sum(interp_sc_x(imax+1,j,:)*RHOmid(imax+1-st2:imax+1+st1,j,k))
              do i=imax+1,imaxo
                 rhoU(i,j,k) = newU
                 U(i,j,k) = newU / rhoi
              end do
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine outflow_correction


! ============================================================ !
! Apply a correction to the outflow to ensure mass consistency !
! Update the exit plane                                        !
!                                                              !
! DOES update the ghost cells                                  !
! NO boundary_update_border required after it                  !
! ============================================================ !
subroutine outflow_clip_negative
  use outflow
  use velocity
  implicit none
  integer :: i,j,k,nflow
  
  ! Nothing to do if:
  ! -> no outlet
  ! -> not last proc in x
  ! -> periodic in x
  if (noutlet.eq.0 .or. iproc.ne.npx .or. xper.eq.1) return
  
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do nflow=1,noutlet
        do j=max(jmino_,outlet(nflow)%jmino),min(jmaxo_,outlet(nflow)%jmaxo)
           do i=imax+1,imaxo
              if (rhoU(i,j,k) .lt. 0.0_WP) then
                 rhoU(i,j,k) = 0.0_WP
                 U(i,j,k) = 0.0_WP
              end if
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine outflow_clip_negative


