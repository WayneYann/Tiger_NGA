module monitor_conservation
  use precision
  implicit none
  
  ! Conservation information
  real(WP) :: mass,mom_x,mom_y,mom_z,energy
  real(WP), dimension(:), pointer :: cons_sc
  
end module monitor_conservation


! ========================================== !
! Initialize the monitor/conservation module !
! ========================================== !
subroutine monitor_conservation_init
  use monitor_conservation
  use data
  implicit none
  
  integer :: isc
  
  ! Create a file to monitor at each timestep
  call monitor_create_file_step('conservation',5+nscalar)
  call monitor_set_header(1,'Mass','r')
  call monitor_set_header(2,'Momentum_x','r')
  call monitor_set_header(3,'Momentum_y','r')
  call monitor_set_header(4,'Momentum_z','r')
  call monitor_set_header(5,'Kinetic energy','r')
  if (nscalar.ge.1) then
     ! Allocate
     allocate(cons_sc(nscalar))
     ! Create header
     do isc=1,nscalar
        call monitor_set_header(5+isc,SC_name(isc),'r')
     end do
  end if
  
  return
end subroutine monitor_conservation_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_conservation_compute
  use monitor_conservation
  use data
  use metric_generic
  use parallel
  implicit none
  
  real(WP) :: buf1,buf2,buf3,buf4,buf5
  integer  :: isc,i,j,k
  
  ! Select the right file
  call monitor_select_file('conservation')
  
  ! Mass and energy
  buf1 = 0.0_WP
  buf5 = 0.0_WP
  !$OMP PARALLEL DO REDUCTION(+:buf1,buf5)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + vol(i,j) * RHO(i,j,k)
           buf5 = buf5 + vol(i,j) * 0.5_WP*( &
                sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)*rhoU(i-st1:i+st2,j,k)) + &
                sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)*rhoV(i,j-st1:j+st2,k)) + &
                sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)*rhoW(i,j,k-st1:k+st2)) )
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_sum(buf1,mass  )
  call parallel_sum(buf5,energy)
  
  ! Momentum
  buf2 = 0.0_WP
  buf3 = 0.0_WP
  buf4 = 0.0_WP
  if (icyl.eq.1) then
     !$OMP PARALLEL DO REDUCTION(+:buf2,buf3,buf4)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf2 = buf2 + vol(i,j) * sum(interp_u_xm(i,j,:)*rhoU(i-st1:i+st2,j,k))
              buf3 = buf3 + vol(i,j) * sum(interp_v_ym(i,j,:)*rhoV(i,j-st1:j+st2,k))*cos(zm(k))
              buf4 = buf4 + vol(i,j) * sum(interp_v_ym(i,j,:)*rhoV(i,j-st1:j+st2,k))*sin(zm(k))
              buf3 = buf3 - vol(i,j) * sum(interp_w_zm(i,j,:)*rhoW(i,j,k-st1:k+st2)*sin(z(k-st1:k+st2)))
              buf4 = buf4 + vol(i,j) * sum(interp_w_zm(i,j,:)*rhoW(i,j,k-st1:k+st2)*cos(z(k-st1:k+st2)))
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO REDUCTION(+:buf2,buf3,buf4)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf2 = buf2 + vol(i,j) * sum(interp_u_xm(i,j,:)*rhoU(i-st1:i+st2,j,k))
              buf3 = buf3 + vol(i,j) * sum(interp_v_ym(i,j,:)*rhoV(i,j-st1:j+st2,k))
              buf4 = buf4 + vol(i,j) * sum(interp_w_zm(i,j,:)*rhoW(i,j,k-st1:k+st2))
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  call parallel_sum(buf2,mom_x )
  call parallel_sum(buf3,mom_y )
  call parallel_sum(buf4,mom_z )
  
  ! Transfer to monitor
  call monitor_set_single_value(1,mass)
  call monitor_set_single_value(2,mom_x)
  call monitor_set_single_value(3,mom_y)
  call monitor_set_single_value(4,mom_z)
  call monitor_set_single_value(5,energy)
  
  ! Scalars
  if (nscalar.ge.1) then
     do isc=1,nscalar
        ! Compute int(rhoSC)
        buf1 = 0.0_WP
        !$OMP PARALLEL DO REDUCTION(+:buf1)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 buf1 = buf1 + vol(i,j) * RHO(i,j,k)*SC(i,j,k,isc)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        call parallel_sum(buf1,cons_sc(isc))
        ! Transfer to monitor
        call monitor_set_single_value(5+isc,cons_sc(isc))
     end do
  end if
  
  return
end subroutine monitor_conservation_compute
