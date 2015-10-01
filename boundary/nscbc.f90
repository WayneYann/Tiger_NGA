module nscbc
  use data
  use parallel
  use boundary
  use borders
  implicit none

  ! NSCBC Parameters
  real(WP), parameter :: sigma = 0.28_WP
  real(WP) :: P_inf

end module nscbc


! ================ !
! Initialize NSCBC !
! ================ !
subroutine nscbc_init
  use nscbc
  implicit none

  call parser_read('NSCBC reference pressure',P_inf)

  return
end subroutine nscbc_init


! ======================= !
! NSCBC inlet for density !
! ======================= !
subroutine nscbc_inflow_scalar
  use scalar
  use velocity
  use nscbc
  use metric_velocity_conv
  use time_info
  use memory
  use combustion
  implicit none

  integer :: i,j,k
  real(WP) :: L_1, L_2, L_5
  real(WP) :: a_m,M_m
  real(WP) :: newRHO

  if (ntime.ne.0) then
     call combustion_soundspeed(tmp1)
     call combustion_soundspeed_old(tmp2)
     !$OMP PARALLEL DO PRIVATE(L_1,L_2,L_5,a_m,M_m,newRHO)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (mask(imin,j).ne.1) then
              a_m = 0.5_WP*(tmp1(imin-1,j,k)+tmp2(imin-1,j,k))
              M_m = U(imin,j,k) / a_m

              ! Outgoing acoustic wave
              L_1 = (U(imin,j,k)-a_m) * &
                   ((0.5_WP*(P(imin,j,k)+Pold(imin,j,k))-0.5_WP*(P(imin-1,j,k)+Pold(imin-1,j,k)))*dxmi(imax) - &
                     0.5_WP*(RHO(imin-1,j,k)+RHOold(imin-1,j,k))*a_m* &
                     (U(imin+1,j,k)-U(imin,j,k))*dxi(imin))

              ! Incoming acoustic wave
              L_5 = L_1 - 2.0_WP*0.5_WP*(RHO(imin-1,j,k)+RHOold(imin-1,j,k))*a_m* &
                   (U(imin,j,k)-Uold(imin,j,k))/dt_uvw

              ! Incoming entropic wave (energy or temperature)
              if (isc_E.ne.0) then
                 L_2 = 0.5_WP*Ggas(i,j,k)*(L_5+L_1)+0.5_WP*(RHO(imin-1,j,k)+RHOold(imin-1,j,k))*a_m**2/ &
                      (0.5_WP*(SC(imin-1,j,k,isc_E)+SCold(imin-1,j,k,isc_E)))* &
                      (SC(imin-1,j,k,isc_E)-SCold(imin-1,j,k,isc_E))/dt
              end if
              if (isc_T.ne.0) then
                 L_2 = 0.5_WP*Ggas(i,j,k)*(L_5+L_1)+0.5_WP*(RHO(imin-1,j,k)+RHOold(imin-1,j,k))*a_m**2/ &
                      (0.5_WP*(SC(imin-1,j,k,isc_T)+SCold(imin-1,j,k,isc_T)))* &
                      (SC(imin-1,j,k,isc_T)-SCold(imin-1,j,k,isc_T))/dt
              end if

              ! Compute and update density
              newRHO = RHOold(imin-1,j,k) - dt/a_m**2*(L_2+0.5_WP*(L_1+L_5)) - &
                   dt*sum(divc_v(imin,j,:)*rhoV(imin,j-stc1:j+stc2,k)) - &
                   dt*sum(divc_w(imin,j,:)*rhoW(imin,j,k-stc1:k+stc2))
              do i=imino,imin-1
                 RHO(i,j,k) = newRHO
              end do
           end if
        end do
     end do
     !$OMP END PARALLEL DO
  else
     ! Neumann condition

     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino,imin-1
              RHO(i,j,k) = RHO(imin,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  return
end subroutine nscbc_inflow_scalar


! ========================= !
! NSCBC outlet for velocity !
! ========================= !
subroutine nscbc_outflow_velocity
  use velocity
  use nscbc
  use metric_generic
  use metric_velocity_conv
  use time_info
  use memory
  implicit none

  integer :: i,j,k
  real(WP) :: L_1, L_2, L_3, L_4, L_5
  real(WP) :: a_f,M_f
  real(WP) :: rhox,rhoy,rhoz
  real(WP) :: newU,newV,newW

  ! Only cold for now, need to get the native energy term...
  ! Semi-implicit right now; should probably make fully implicit...

  if (ntime.ne.0) then
     call combustion_soundspeed_old(tmp1)
     !$OMP PARALLEL DO PRIVATE(L_1,L_2,L_3,L_4,L_5,a_f,M_f,rhox,rhoy,rhoz,newU,newV,newW)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (mask(imax,j).ne.1) then
              a_f = sum(interp_sc_x(imax+1,j,:)*tmp1(imax+1-st2:imax+1+st1,j,k))
              M_f = 0.5_WP*(U(imax+1,j,k)+Uold(imax+1,j,k)) / a_f
           
              ! Incoming acoustic wave
              L_1 = sigma*(1.0_WP-M_f**2)*a_f/xL * &
                   (sum(interp_sc_x(imax+1,j,:)*Pold(imax+1-st2:imax+1+st1,j,k)) - P_inf)

              ! Outgoing entropic wave
              L_2 = 0.5_WP*(U(imax+1,j,k)+Uold(imax+1,j,k))*( &
                   a_f**2*(RHOold(imax+1,j,k)-RHOold(imax,j,k))*dxmi(imax)- &
                   (Pold(imax+1,j,k)-Pold(imax,j,k))*dxmi(imax))
              
              ! Outgoing vortical wave (y-component)
              L_3 = sum(interp_Ju_y(imax+1,j,:)*0.5_WP*(U(imax+1-stc2:imax+1+stc1,j,k)+&
                      Uold(imax+1-stc2:imax+1+stc1,j,k)))* & 
                   0.5_WP*((V(imax+1,j,k)+Vold(imax+1,j,k))-(V(imax,j,k)+Vold(imax,j,k)))*dxmi(imax)
           
              ! Outgoing vortical wave (w-component)
              L_4 = sum(interp_Ju_z(imax+1,j,:)*0.5_WP*(U(imax+1-stc2:imax+1+stc1,j,k)+&
                      Uold(imax+1-stc2:imax+1+stc1,j,k)))* &
                   0.5_WP*((W(imax+1,j,k)+Wold(imax+1,j,k))-(W(imax,j,k)+Wold(imax,j,k)))*dxmi(imax)
           
              ! Outgoing acoustic wave
              L_5 = (0.5_WP*(U(imax+1,j,k)+Uold(imax+1,j,k)+a_f) * &
                   ((Pold(imax+1,j,k)-Pold(imax,j,k))*dxmi(imax) + &
                    sum(interp_sc_x(imax+1,j,:)*RHOold(imax+1-st2:imax+1+st1,j,k))* &
                       a_f*(0.5_WP*(U(imax+1,j,k)+Uold(imax+1,j,k))-0.5_WP*(U(imax,j,k)+Uold(imax,j,k)))*dxi(imax)))
              
              ! Compute new velocities (U is momentum)
              newU = rhoUold(imax+1,j,k) - dt_uvw/a_f*(M_f*L_2+0.5_WP*((M_f-1.0_WP)*L_1+(M_f+1.0_WP)*L_5))
              newV = Vold(imax+1,j,k) - dt_uvw*(L_3)
              newW = Wold(imax+1,j,k) - dt_uvw*(L_4)
           
              ! Compute the interpolated densities
              rhox = sum(interp_sc_x(imax+1,j,:)*RHOmid(imax+1-st2:imax+1+st1,j,k))
              rhoy = sum(interp_sc_y(imax+1,j,:)*RHOmid(imax+1,j-st2:j+st1,k))
              rhoz = sum(interp_sc_z(imax+1,j,:)*RHOmid(imax+1,j,k-st2:k+st1))
           
              ! Update the velocities and momenta
              do i=imax+1,imaxo
                 rhoU(i,j,k) = newU
                 U(i,j,k) = newU / rhox
                 V(i,j,k) = newV
                 rhoV(i,j,k) = newV * rhoy
                 W(i,j,k) = newW
                 rhoW(i,j,k) = newW * rhoz
              end do
           end if
        end do
     end do
     !$OMP END PARALLEL DO
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

  return
end subroutine nscbc_outflow_velocity


! ======================== !
! NSCBC outlet for scalars !
! ======================== !
subroutine nscbc_outflow_scalar
  use scalar
  use nscbc
  use metric_generic
  use metric_velocity_conv
  use time_info
  use memory
  use combustion
  implicit none

  integer :: i,j,k,isc
  real(WP) :: L_1, L_2, L_5, L_6
  real(WP) :: a_m,M_m
  real(WP) :: newRHO,newE,newT,newSC

  ! Only cold for now, need to get the native energy term...
  ! Semi-implicit right now; should probably make fully implicit...

  if (ntime.ne.0) then
     call combustion_soundspeed(tmp1)
     call combustion_soundspeed_old(tmp2)
     !$OMP PARALLEL DO PRIVATE(L_1,L_2,L_5,L_6,a_m,M_m,newRHO,newE,newT,newSC)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (mask(imax,j).ne.1) then
              a_m = 0.5_WP*(tmp1(imax+1,j,k)+tmp2(imax+1,j,k))
              M_m = U(imax+1,j,k) / a_m

             ! Incoming acoustic wave
              L_1 = sigma*(1.0_WP-M_m**2)*a_m/xL * &
                   (0.5_WP*(P(imax+1,j,k)+Pold(imax+1,j,k)) - P_inf)

              ! Outgoing entropic wave
              L_2 = U(imax+1,j,k)*( &
                   a_m**2*(0.5_WP*(RHO(imax+1,j,k)+RHOold(imax+1,j,k))-0.5_WP*(RHO(imax,j,k)+RHOold(imax,j,k)))*dxmi(imax)- &
                   (0.5_WP*(P(imax+1,j,k)+Pold(imax+1,j,k))-0.5_WP*(P(imax,j,k)+Pold(imax,j,k)))*dxmi(imax))
              
              ! Outgoing acoustic wave
              L_5 = (U(imax+1,j,k)+a_m) * &
                   ((0.5_WP*(P(imax+1,j,k)+Pold(imax+1,j,k))-0.5_WP*(P(imax,j,k)+Pold(imax,j,k)))*dxmi(imax) + &
                     0.5_WP*(RHO(imax+1,j,k)+RHOold(imax+1,j,k))*a_m* &
                     (U(imax+1,j,k)-U(imax,j,k))*dxi(imax))

              ! Compute and update density
              newRHO = RHOold(imax+1,j,k) - dt/a_m**2*(L_2+0.5_WP*(L_1+L_5))
              do i=imax+1,imaxo
                 RHO(i,j,k) = newRHO
              end do
              
              ! Compute and update energy
              if (isc_E.ne.0) then
                 newE = SCold(imax+1,j,k,isc_E) - dt*0.5_WP*(SC(imax+1,j,k,isc_E)+SCold(imax+1,j,k,isc_E))&
                      /(0.5*(RHO(imax+1,j,k)+RHOold(imax+1,j,k)))/a_m**2* &
                      ((Ggas(i,j,k)-1.0_WP)/2.0_WP*(L_5+L_1)-L_2)
                 do i=imax+1,imaxo
                    SC(i,j,k,isc_E) = newE
                 end do
              end if

              ! Compute and update temperature
              if (isc_T.ne.0) then
                 newT = SCold(imax+1,j,k,isc_T) - dt*0.5_WP*(SC(imax+1,j,k,isc_T)+SCold(imax+1,j,k,isc_T))&
                      /(0.5*(RHO(imax+1,j,k)+RHOold(imax+1,j,k)))/a_m**2* &
                      ((Ggas(i,j,k)-1.0_WP)/2.0_WP*(L_5+L_1)-L_2)
                 do i=imax+1,imaxo
                    SC(i,j,k,isc_T) = newT
                 end do
              end if
           
              ! Additional scalars
              do isc=1,nscalar
                 if (isc.ne.isc_E .and. isc.ne.isc_T) then
                    ! Outgoing scalar "wave"
                    L_6 = U(imax+1,j,k)*& 
                         (0.5_WP*(SC(imax+1,j,k,isc)+SCold(imax+1,j,k,isc))-&
                          0.5_WP*(SC(imax,j,k,isc)+SCold(imax+1,j,k,isc)))*dxmi(imax)
                 
                    ! Compute and update scalar
                    newSC = SCold(imax+1,j,k,isc) - dt*(L_6)
                    do i=imax+1,imaxo
                       SC(i,j,k,isc) = newSC
                    end do
                 end if
              end do
           end if
        end do
     end do
     !$OMP END PARALLEL DO
  else
     ! Neumann condition

     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imax+1,imaxo
              RHO(i,j,k) = RHO(imax,j,k)
              SC(i,j,k,:) = SC(imax,j,k,:)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine nscbc_outflow_scalar
