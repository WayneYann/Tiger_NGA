module monitor_flamespeed
  use precision
  use geometry
  use partition
  use borders
  implicit none
  
  ! List of locations
  integer :: nloc_max = 100
  integer :: nloc
  real(WP), dimension(:), pointer :: tloc
  real(WP), dimension(:), pointer :: xloc
  
  ! Old position
  real(WP) :: xloc_old
  
  ! Values to monitor
  real(WP) :: xloc_
  real(WP) :: Tmin,Tmid,Tmax
  real(WP) :: SL,SL_avg

end module monitor_flamespeed



! ======================================== !
! Initialize the monitor/flamespeed module !
! ======================================== !
subroutine monitor_flamespeed_init
  use monitor_flamespeed
  implicit none

  ! Create a file to monitor at each timestep
  call monitor_create_file_step('flame_speed',6)
  call monitor_set_header(1,'Tmin [K]','r')
  call monitor_set_header(2,'Tmid [K]','r')
  call monitor_set_header(3,'Tmax [K]','r')
  call monitor_set_header(4,'flame_loc','r')
  call monitor_set_header(5,'SL [cm/s]','r')
  call monitor_set_header(6,'SL_avg [cm/s]','r')
  
  ! Set flame position at the end of the domain
  xloc = xL
  
  ! Allocate arrays
  allocate(tloc(nloc_max))
  allocate(xloc(nloc_max))
  nloc = 1
  xloc(1) = xL
  tloc(1) = 0.0_WP
  
  return
end subroutine monitor_flamespeed_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_flamespeed_compute
  use monitor_flamespeed
  use data
  !use combustion
  use masks
  use time_info
  implicit none
  integer  :: i
  real(WP) :: tmp,alpha
  real(WP) :: a11,a12,a22,y1,y2
  
  ! Get min/max temp
  Tmax = -huge(1.0_WP)
  Tmin = +huge(1.0_WP)
  do i=imin_,imax_
     if (mask(i,jmin_).eq.0) then
        !if (T(i,jmin_,kmin_).gt.Tmax) Tmax = T(i,jmin_,kmin_)
        !if (T(i,jmin_,kmin_).lt.Tmin) Tmin = T(i,jmin_,kmin_)
     end if
  end do
  call parallel_max( Tmax, tmp); Tmax =  tmp
  call parallel_max(-Tmin, tmp); Tmin = -tmp
  
  ! Get medium temperature for the flame location
  Tmid = 0.5_WP*(Tmin+Tmax)
  
  ! Find location of flame
  loop: do i=imin_,imax_
     !if (T(i,jmin_,kmin_).gt.Tmid) exit loop
  end do loop
  if (i.gt.imax_) then
     tmp = -1.0_WP
  else
     !alpha = (Tmid-T(i-1,jmin_,kmin_))/(T(i,jmin_,kmin_)-T(i-1,jmin_,kmin_))
     tmp  = alpha*xm(i) + (1.0_WP-alpha)*xm(i-1)
  end if
  call parallel_max(tmp, xloc_)
  
  ! Compute instantenous speed
  SL = 100.0_WP*(xloc_old-xloc_)/dt
  xloc_old = xloc_
  
  ! Store position in array
  nloc = nloc + 1
  if (nloc.gt.nloc_max) then
     nloc = nloc_max/2
     do i=1,nloc
        xloc(i) = xloc(2*i)
        tloc(i) = tloc(2*i)
     end do
  end if
  xloc(nloc) = xloc_
  tloc(nloc) = time
  
  ! Compute average speed
  a11 = sum(tloc(1:nloc)**2)
  a12 = sum(tloc(1:nloc))
  a22 = real(nloc,WP)
  y1  = sum(xloc(1:nloc)*tloc(1:nloc))
  y2  = sum(xloc(1:nloc))
  SL_avg = -100.0_WP*(a22*y1-a12*y2)/(a11*a22-a12**2+epsilon(1.0_WP))
  
  ! Transfer values to monitor
  call monitor_select_file('flame_speed')
  call monitor_set_single_value(1,Tmin)
  call monitor_set_single_value(2,Tmid)
  call monitor_set_single_value(3,Tmax)
  call monitor_set_single_value(4,xloc)
  call monitor_set_single_value(5,SL)
  call monitor_set_single_value(6,SL_avg)

  return
end subroutine monitor_flamespeed_compute

