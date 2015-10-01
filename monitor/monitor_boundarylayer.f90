module monitor_boundarylayer
  use precision
  use geometry
  use partition
  use borders
  implicit none
  
  ! boundary layer
  real(WP) :: disp_thick,mom_thick,bl_thick
  real(WP) :: Cf,Vinft
  real(WP), dimension(:), pointer :: Ubase

end module monitor_boundarylayer


! =========================================== !
! Initialize the monitor/boundarylayer module !
! =========================================== !
subroutine monitor_boundarylayer_init
  use monitor_boundarylayer
  implicit none

  ! Allocate array
  allocate(Ubase(jmin:jmax))

  ! Create a file to monitor at each timestep
  call monitor_create_file_step('boundarylayer',5)
  call monitor_set_header(1,'bl_thick','r')
  call monitor_set_header(2,'disp_thick','r')
  call monitor_set_header(3,'mom_thick','r')
  call monitor_set_header(4,'Cf','r')
  call monitor_set_header(5,'Vinft','r')

  return
end subroutine monitor_boundarylayer_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_boundarylayer_compute
  use monitor_boundarylayer
  use data
  use parser
  use metric_generic
  use metric_velocity_visc
  implicit none
  
  real(WP) :: Uinft,tmp
  integer  :: i,j,k
  
  ! Compute boundary layer stats
  Ubase = 0.0_WP
  do j=jmin_,jmax_
     Ubase(j) = sum(U(imin_:imax_,j,kmin_:kmax_))
  end do
  do j=jmin,jmax
     call parallel_sum(Ubase(j),tmp)
     Ubase(j) = tmp
  end do
  Ubase = Ubase/real(nx*nz,WP)
  
  ! Calculate the displacement and momentum thickness
  Uinft = Ubase(jmax) + epsilon(1.0_WP)
  disp_thick = 0.0_WP
  mom_thick = 0.0_WP
  do j=jmin+1,jmax
     disp_thick = disp_thick + (1.0_WP-Ubase(j)/Uinft) * dy(j)
     mom_thick  = mom_thick  + Ubase(j)/Uinft * (1.0_WP-Ubase(j)/Uinft) * dy(j)
  end do
  
  ! Compute boundary layer thickness
  loopj:do j=jmax-1,jmin+1,-1
     if (Ubase(j).le.0.99_WP*Uinft) exit loopj
  end do loopj
  tmp = (0.99_WP*Uinft-Ubase(j))/(Ubase(j+1)-Ubase(j)+epsilon(1.0_WP))
  bl_thick = tmp*ym(j+1)+(1.0_WP-tmp)*ym(j)
  
  ! Wall friction
  Cf = 0.0_WP
  if (jproc.eq.1) then
     j = jmin+1
     do k=kmin_,kmax_
        do i=imin_,imax_
           Cf = Cf + dz_v(j) * dxm(i-1) * &
                sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k))* &
                sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k))
        end do
     end do
  end if
  call parallel_sum(Cf,tmp)
  Cf = tmp / ((x(imax+1)-x(imin))*(z(kmax+1)-z(kmin)))
  Cf = 2.0_WP * Cf / Uinft**2
  
  ! Compute the Vertical exit velocity
  call parallel_sum(sum(V(imin_:imax_,jmax+1,kmin_:kmax_)),Vinft)
  Vinft = Vinft / real(nx*nz,WP)
  
  ! Transfer values to monitor
  call monitor_select_file('boundarylayer')
  call monitor_set_single_value(1,bl_thick)
  call monitor_set_single_value(2,disp_thick)
  call monitor_set_single_value(3,mom_thick)
  call monitor_set_single_value(4,Cf)
  call monitor_set_single_value(5,Vinft)

  return
end subroutine monitor_boundarylayer_compute
