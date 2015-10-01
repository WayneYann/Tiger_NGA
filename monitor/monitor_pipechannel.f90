module monitor_pipechannel
  use precision
  use geometry
  use partition
  use borders
  implicit none

  ! Values to monitor
  real(WP), dimension(:), pointer :: Cfx
  real(WP), dimension(:), pointer :: Cfz
  real(WP), dimension(:), pointer :: umean
  real(WP), dimension(:), pointer :: wmean
  real(WP), dimension(:), pointer :: RHOmean

  ! Array of values to transfer
  real(WP), dimension(:), pointer :: mval

end module monitor_pipechannel



! ========================================= !
! Initialize the monitor/pipechannel module !
! ========================================= !
subroutine monitor_pipechannel_init
  use monitor_pipechannel
  use config
  implicit none
  integer :: nflow
  character(7) :: header

  ! Allocate array
  allocate(Cfx(ninlet))
  allocate(Cfz(ninlet))
  allocate(umean(ninlet))
  allocate(wmean(ninlet))
  allocate(RHOmean(ninlet))
  allocate(mval(4*ninlet))
   
  ! Create a file to monitor at each timestep
  if (icyl.eq.1) then
     call monitor_create_file_step('pipe',4*ninlet)
  else
     call monitor_create_file_step('channel',4*ninlet)
  end if
  do nflow=1,ninlet
     write(header,'(a3,i1)') 'Cfx_',nflow
     call monitor_set_header(nflow+0*ninlet,trim(header),'r')
     write(header,'(a3,i1)') 'Cfz_',nflow
     call monitor_set_header(nflow+1*ninlet,trim(header),'r')
     write(header,'(a6,i1)') 'Umean_',nflow
     call monitor_set_header(nflow+2*ninlet,trim(header),'r')
     write(header,'(a6,i1)') 'Wmean_',nflow
     call monitor_set_header(nflow+3*ninlet,trim(header),'r')
  end do
  
  return
end subroutine monitor_pipechannel_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_pipechannel_compute
  use monitor_pipechannel
  use data
  use parser
  use metric_generic
  use metric_velocity_visc
  implicit none
  
  integer :: i,j,k,nflow
  real(WP) :: tmp_x,tmp_z,tmp_s,surf
  
  ! Compute mean densities
  do nflow=1,ninlet
     tmp_x = 0.0_WP
     do j = max(inlet(nflow)%jmin,jmin_),min(inlet(nflow)%jmax,jmax_)
        !$OMP PARALLEL DO REDUCTION(+:tmp_x)
        do k=kmin_,kmax_
           do i=imin_,imax_
              tmp_x = tmp_x + RHO(i,j,k)*dA(j)*dxm(i-1)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
     call parallel_sum(tmp_x,RHOmean(nflow))
     RHOmean(nflow) = RHOmean(nflow)/inlet(nflow)%A
  end do
  RHOmean = RHOmean/(xm(imax)-xm(imin-1))
  
  ! Compute mass flow rate between the walls
  do nflow=1,ninlet
     tmp_x = 0.0_WP
     tmp_z = 0.0_WP

     do j = max(inlet(nflow)%jmin,jmin_),min(inlet(nflow)%jmax,jmax_)
        !$OMP PARALLEL DO REDUCTION(+:tmp_x,tmp_z)
        do k=kmin_,kmax_
           do i=imin_,imax_
              tmp_x = tmp_x + U(i,j,k)*dA(j)*dxm(i-1)
              tmp_z = tmp_z + W(i,j,k)*dA(j)*dxm(i-1)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
     
     call parallel_sum(tmp_x,umean(nflow))
     call parallel_sum(tmp_z,wmean(nflow))
     
     umean(nflow) = umean(nflow)/inlet(nflow)%A
     wmean(nflow) = wmean(nflow)/inlet(nflow)%A
  end do
  
  umean = umean/(xm(imax)-xm(imin-1))
  wmean = wmean/(xm(imax)-xm(imin-1))
  
  ! Compute friction force on the walls
  do nflow=1,ninlet

     ! Lower wall
     tmp_x = 0.0_WP
     tmp_z = 0.0_WP
     tmp_s = 0.0_WP
     if (inlet(nflow)%jmin.ge.jmin_ .and. inlet(nflow)%jmin.le.jmax_) then 
        j = max(inlet(nflow)%jmin,jmin_)
        if (mask(imin_,j-1).eq.1) then
           !$OMP PARALLEL DO REDUCTION(+:tmp_s,tmp_x,tmp_z)
           do k=kmin_,kmax_
              do i=imin_,imax_
                 tmp_s = tmp_s + dz_v(j) * dxm(i-1) 

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
           !$OMP PARALLEL DO REDUCTION(+:tmp_s,tmp_x,tmp_z)
           do k=kmin_,kmax_
              do i=imin_,imax_
                 tmp_s = tmp_s + dz_v(j) * dxm(i-1) 

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

     call parallel_sum(tmp_x,Cfx(nflow))
     call parallel_sum(tmp_z,Cfz(nflow))
     call parallel_sum(tmp_s,surf)

     if (umean(nflow).eq.0.0_WP) then
        Cfx(nflow) = 0.0_WP
     else
        Cfx(nflow) = 2.0_WP*Cfx(nflow)/(RHOmean(nflow)*umean(nflow)**2) / surf
     end if
     if (wmean(nflow).eq.0.0_WP) then
        Cfz(nflow) = 0.0_WP
     else
        Cfz(nflow) = 2.0_WP*Cfz(nflow)/(RHOmean(nflow)*umean(nflow)**2) / surf
     end if
  end do
  
  ! Transfer values to monitor
  if (icyl.eq.1) then
     call monitor_select_file('pipe')
  else
     call monitor_select_file('channel')
  end if
  do nflow=1,ninlet
     mval(nflow+0*ninlet) = Cfx(nflow)
     mval(nflow+1*ninlet) = Cfz(nflow)
     mval(nflow+2*ninlet) = umean(nflow)
     mval(nflow+3*ninlet) = wmean(nflow)
  end do
  call monitor_set_array_values(mval)

  return
end subroutine monitor_pipechannel_compute

