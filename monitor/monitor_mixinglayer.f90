module monitor_mixinglayer
  use precision
  use geometry
  use partition
  implicit none

  ! Mixing layer
  real(WP) :: disp_thick,mom_thick,vort_thick,lambda
  real(WP) :: Re_mom,Re_vort,Re_lambda
  real(WP), dimension(:), pointer :: Ubase

end module monitor_mixinglayer


! ========================================= !
! Initialize the monitor/mixinglayer module !
! ========================================= !
subroutine monitor_mixinglayer_init
  use monitor_mixinglayer
  implicit none
  
  ! Allocate array
  allocate(Ubase(jmin:jmax))

  ! Create a file to monitor at each timestep
  call monitor_create_file_step('mixinglayer',7)
  call monitor_set_header(1,'disp_thick','r')
  call monitor_set_header(2,'mom_thick','r')
  call monitor_set_header(3,'vort_thick','r')
  call monitor_set_header(4,'lambda','r')
  call monitor_set_header(5,'Re_mom','r')
  call monitor_set_header(6,'Re_vort','r')
  call monitor_set_header(7,'Re_lambda','r')

  return
end subroutine monitor_mixinglayer_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_mixinglayer_compute
  use monitor_mixinglayer
  use data
  use parser
  implicit none
    
  real(WP) :: nu,up2,dupdx2,U1,U2,tmp,uprime1,uprime2,max_der
  integer :: i,j,k
  
  ! Compute mixing layer stats
  Ubase = 0.0_WP
  do j=jmin_,jmax_
     Ubase(j) = sum(U(imin_:imax_,j,kmin_:kmax_))
  end do
  do j=jmin,jmax
     call parallel_sum(Ubase(j),tmp)
     Ubase(j) = tmp
  end do
  Ubase = Ubase/real(nx*nz,WP)
  
  ! Calculate the momentum thickness and Re_mom
  U1 = Ubase(jmax)
  U2 = Ubase(jmin)
  
  disp_thick = 0.0_WP
  mom_thick  = 0.0_WP
  do j=jmin,jmax
     disp_thick = disp_thick + (0.5_WP*(U1-U2)-abs(Ubase(j)-0.5_WP*(U1+U2))) * dy(j)
     mom_thick  = mom_thick  + (U1-Ubase(j))*(Ubase(j)-U2)*dy(j)
  end do
  mom_thick  = mom_thick/((U1-U2)**2+epsilon(1.0_WP))
  disp_thick = 2.0_WP*disp_thick/(U1-U2+epsilon(1.0_WP))
  
  tmp = sum(VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_)/RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  call parallel_sum(tmp,nu)
  nu = nu /real(nx*ny*nz,WP)
  re_mom = mom_thick*(U1-U2)/(nu+epsilon(1.0_WP)) 
  
  ! Calculate the vorticity thickness and Re_vort
  max_der = 0.0_WP
  do j=jmin,jmax-1
     tmp = abs(Ubase(j+1)-Ubase(j))*dymi(j)
     if (tmp.gt.max_der) max_der = tmp
  end do
  vort_thick = (U1-U2)/(max_der+epsilon(1.0_WP))
  re_vort = vort_thick*(U1-U2)/(nu+epsilon(1.0_WP)) 
  
  ! Calculate the taylor microscale (lambda) and Re_lambda
  up2    = 0.0_WP
  dupdx2 = 0.0_WP
  j = (nyo-1)/2+1
  if (j.ge.jmin_ .and. j.le.jmax_) then
     do k=kmin_,kmax_
        do i=imin_,imax_
           uprime1 = U(i,j,k)-Ubase(j)
           uprime2 = U(i+1,j,k)-Ubase(j)
           up2 = up2 + ((uprime1+uprime2)/2.0_WP)**2
           dupdx2 = dupdx2 + ((uprime2-uprime1)/dx(i))**2
        end do
     end do
  else
     up2    = 0.0_WP
     dupdx2 = 0.0_WP
  end if
  call parallel_sum(up2,tmp)
  up2 = tmp
  call parallel_sum(dupdx2,tmp)
  dupdx2 = tmp
  lambda = sqrt( up2/(dupdx2+epsilon(1.0_WP)) )
!  re_lambda = lambda*(U1-U2)/nu
  re_lambda = lambda*up2**(0.5)/(nu+epsilon(1.0_WP))
  
  ! Transfer values to monitor
  call monitor_select_file('mixinglayer')
  call monitor_set_single_value(1,disp_thick)
  call monitor_set_single_value(2,mom_thick)
  call monitor_set_single_value(3,vort_thick)
  call monitor_set_single_value(4,lambda)
  call monitor_set_single_value(5,Re_mom)
  call monitor_set_single_value(6,Re_vort)
  call monitor_set_single_value(7,Re_lambda)

  return
end subroutine monitor_mixinglayer_compute

