module monitor_hit
  use precision
  implicit none

  ! Homogeneous isotropic turbulence
  real(WP) :: TKE,eps,n_kolmogorov,keta,tau_eddy,L_int,Re_L,eps_favre
  real(WP), dimension(:), pointer :: variance,zmean
  real(WP), dimension(:), pointer :: scalar_stat
  real(WP) :: urms2,nu_favre
  real(WP) :: lambda,Re_lambda
  
end module monitor_hit


! ================================= !
! Initialize the monitor/hit module !
! ================================= !
subroutine monitor_hit_init
  use monitor_hit
  use data
  implicit none
  integer :: isc
  
  ! Create a file to monitor at each timestep
  call monitor_create_file_step('hit',8)
  call monitor_set_header(1,'TKE','r')
  call monitor_set_header(2,'urms','r')
  call monitor_set_header(3,'epsilon','r')
  call monitor_set_header(4,'tau_eddy','r')
  call monitor_set_header(5,'eta','r')
  call monitor_set_header(6,'keta','r')
  call monitor_set_header(7,'Re_turb','r')
  call monitor_set_header(8,'Re_lambda','r')
  
  if (nscalar.ge.1) then
     ! Allocate
     allocate(zmean(nscalar))
     allocate(variance(nscalar))
     allocate(scalar_stat(2*nscalar))
     
     ! Create header
     call monitor_create_file_step('hit_scalar',2*nscalar)
     do isc=1,nscalar
        call monitor_set_header(2*isc-1,'avg_'//SC_name(isc),'r')
        call monitor_set_header(2*isc+0,'var_'//SC_name(isc),'r')
     end do
  end if
  
  return
end subroutine monitor_hit_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_hit_compute
  use monitor_hit
  use data
  use metric_generic
  use memory
  use math
  implicit none
  
  real(WP) :: buf1,buf2
  real(WP) :: RHOmean,div1,div2,div3
  integer  :: isc,i,j,k
  
  ! Favre(urms^2) = <rhoU.U/3> / <RHO>
  buf1 = 0.0_WP
  buf2 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + vol(i,j) * ( &
                sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)*rhoU(i-st1:i+st2,j,k)) + &
                sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)*rhoV(i,j-st1:j+st2,k)) + &
                sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)*rhoW(i,j,k-st1:k+st2)) )
           buf2 = buf2 + RHO(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,urms2)
  urms2 = urms2/vol_total
  call parallel_sum(buf2,RHOmean)
  RHOmean = RHOmean/vol_total
  urms2 = urms2/RHOmean
  urms2 = urms2/3.0_WP
  
  ! TKE = 3/2 * urms^2
  TKE = 1.5_WP*urms2
  
  ! Mean and variance of the scalars
  do isc=1,nscalar
     buf1 = 0.0_WP
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf1 = buf1 + RHO(i,j,k)*SC(i,j,k,isc) * vol(i,j)
           end do
        end do
     end do
     call parallel_sum(buf1,zmean(isc))
     zmean(isc) = zmean(isc)/vol_total
     zmean(isc) = zmean(isc)/RHOmean
     buf1 = 0.0_WP
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf1 = buf1 + RHO(i,j,k)*(SC(i,j,k,isc)-zmean(isc))**2 * vol(i,j)
           end do
        end do
     end do
     call parallel_sum(buf1,variance(isc))
     variance(isc) = variance(isc)/vol_total
     variance(isc) = variance(isc)/RHOmean
     scalar_stat(2*isc-1) = zmean(isc)
     scalar_stat(2*isc)   = variance(isc)
  end do
  
  ! Favre average of nu = <mu>/<RHO>
  buf1 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + VISC(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,nu_favre)
  nu_favre = nu_favre/vol_total
  nu_favre = nu_favre/RHOmean
  
  ! Interpolate velocities
  call interpolate_velocities
  
  ! Compute VISC.S
  !     ( 1 4 6 )
  ! S = ( 4 2 5 )
  !     ( 6 5 3 )
  call strainrate_compute(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  tmp1 = 2.0_WP*VISC*tmp1
  tmp2 = 2.0_WP*VISC*tmp2
  tmp3 = 2.0_WP*VISC*tmp3
  tmp4 = 2.0_WP*VISC*tmp4
  tmp5 = 2.0_WP*VISC*tmp5
  tmp6 = 2.0_WP*VISC*tmp6
  
  ! Compute -U.(div(VISC.S))
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! div(VISC.S)
           div1=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp1(i-1,j,k)+interp_sc_x(i,j,0)*tmp1(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp1(i,j,k)+interp_sc_x(i+1,j,0)*tmp1(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp4(i,j-1,k)+interp_sc_y(i,j,0)*tmp4(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp4(i,j,k)+interp_sc_y(i,j+1,0)*tmp4(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp6(i,j,k-1)+interp_sc_z(i,j,0)*tmp6(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp6(i,j,k)+interp_sc_z(i,j,0)  *tmp6(i,j,k+1))
           div2=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp4(i-1,j,k)+interp_sc_x(i,j,0)*tmp4(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp4(i,j,k)+interp_sc_x(i+1,j,0)*tmp4(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp2(i,j-1,k)+interp_sc_y(i,j,0)*tmp2(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp2(i,j,k)+interp_sc_y(i,j+1,0)*tmp2(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp5(i,j,k-1)+interp_sc_z(i,j,0)*tmp5(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp5(i,j,k)+interp_sc_z(i,j,0)  *tmp5(i,j,k+1))
           div3=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp6(i-1,j,k)+interp_sc_x(i,j,0)*tmp6(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp6(i,j,k)+interp_sc_x(i+1,j,0)*tmp6(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp5(i,j-1,k)+interp_sc_y(i,j,0)*tmp5(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp5(i,j,k)+interp_sc_y(i,j+1,0)*tmp5(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp3(i,j,k-1)+interp_sc_z(i,j,0)*tmp3(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp3(i,j,k)+interp_sc_z(i,j,0)  *tmp3(i,j,k+1))
           ! div(VISC.S)
           div1=+sum(div_u(i,j,:)*tmp1(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp4(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp6(i,j,k-st1:k+st2))
           div2=+sum(div_u(i,j,:)*tmp4(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp2(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp5(i,j,k-st1:k+st2))
           div3=+sum(div_u(i,j,:)*tmp6(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp5(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp3(i,j,k-st1:k+st2))
           ! -U.(div(VISC.S))
           tmp7(i,j,k) = &
                - sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)) * div1 &
                - sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)) * div2 &
                - sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)) * div3
           
        end do
     end do
  end do
  
  ! Epsilon = -<U.(div(VISC.S))> 
  buf1 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + tmp7(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,eps)
  eps = eps/vol_total
  eps_favre = eps/RHOmean
  
  if (eps_favre.eq.0.0_WP .or. nu_favre.eq.0.0_WP) then
     eps_favre = 1.0_WP
     nu_favre  = 1.0_WP
  end if
  
  ! n_kolmogorov
  n_kolmogorov = (nu_favre**3/eps_favre)**(0.25_WP)
  keta = Pi*n_kolmogorov*dzi
  
  ! L_int
  L_int = TKE**(1.5_WP)/eps_favre
  
  ! Re_L
  Re_L = TKE**2/(eps_favre*nu_favre)
  
  ! Lambda
  if (Re_L.le.0.0_WP) then
     lambda = 0.0_WP
  else
     lambda = L_int*sqrt(10.0_WP)*Re_L**(-0.5_WP)
  end if
  
  ! Re_lambda
  Re_lambda = urms2**(0.5_WP)*lambda/nu_favre
  
  ! tau_eddy
  tau_eddy = TKE/eps_favre
  
  ! Transfer values to monitor
  call monitor_select_file('hit')
  call monitor_set_single_value(1,TKE)
  call monitor_set_single_value(2,sqrt(urms2))
  call monitor_set_single_value(3,eps_favre)
  call monitor_set_single_value(4,tau_eddy)
  call monitor_set_single_value(5,n_kolmogorov)
  call monitor_set_single_value(6,keta)
  call monitor_set_single_value(7,Re_L)
  call monitor_set_single_value(8,Re_lambda)
  if (nscalar.ge.1) then
     call monitor_select_file('hit_scalar')
     call monitor_set_array_values(scalar_stat)
  end if
  
  return
end subroutine monitor_hit_compute

