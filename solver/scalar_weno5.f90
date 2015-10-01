module scalar_weno5
  use scalar
  use data
  use metric_generic
  implicit none

  ! Ultimate weno5
  logical :: ult

  ! Coefficient for the weno 5th
  real(WP), dimension(:,:,:,:), pointer :: weno5_xp,weno5_yp,weno5_zp
  real(WP), dimension(:,:,:,:), pointer :: weno5_xm,weno5_ym,weno5_zm
  
  ! Coefficient for the weno 5th
  ! -> S0 : scheme 0
  ! -> S1 : scheme 1
  ! -> S2 : scheme 2
  real(WP), dimension(:,:,:), pointer :: S0_xp,S0_yp,S0_xm,S0_ym
  real(WP), dimension(:,:,:), pointer :: S1_xp,S1_yp,S1_xm,S1_ym
  real(WP), dimension(:,:,:), pointer :: S2_xp,S2_yp,S2_xm,S2_ym
  real(WP), dimension(:),     pointer :: S0_zp,S1_zp,S2_zp,S0_zm,S1_zm,S2_zm
  
end module scalar_weno5


! ============================ !
! Initialize the weno 5 module !
! ============================ !
subroutine scalar_weno5_init
  use scalar_weno5
  use masks
  use math
  use parser
  use string
  implicit none
  
  integer  :: i,j,k,loc
  real(WP) :: a0,a1,a2
  character(str_medium) :: scheme
  
  ! Ultimate weno5 or not
  ult = .false.
  call parser_read('Scalar scheme',scheme)
  loc = index(scheme,' ')
  if (trim(scheme(loc+1:)).eq.'ultimate') ult = .true.

  ! Allocate the metrics
  allocate(weno5_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-3:+1))
  allocate(weno5_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-3:+1))
  allocate(weno5_zp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-3:+1))
  allocate(weno5_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2:+2))
  allocate(weno5_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2:+2))
  allocate(weno5_zm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,kmin_-st1:kmax_+st2,-2:+2))
  
  allocate(S0_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S1_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S2_xp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S0_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S1_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S2_yp(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-3:+1))
  allocate(S0_zp(-3:+1))
  allocate(S1_zp(-3:+1))
  allocate(S2_zp(-3:+1))
  
  allocate(S0_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S1_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S2_xm(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S0_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S1_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S2_ym(imin_-st1:imax_+st2,jmin_-st1:jmax_+st2,-2:+2))
  allocate(S0_zm(-2:+2))
  allocate(S1_zm(-2:+2))
  allocate(S2_zm(-2:+2))
  
  ! Interpolation at the faces
  S0_xp=0.0_WP;S0_yp=0.0_WP;S0_zp=0.0_WP;S0_xm=0.0_WP;S0_ym=0.0_WP;S0_zm=0.0_WP
  S1_xp=0.0_WP;S1_yp=0.0_WP;S1_zp=0.0_WP;S1_xm=0.0_WP;S1_ym=0.0_WP;S1_zm=0.0_WP
  S2_xp=0.0_WP;S2_yp=0.0_WP;S2_zp=0.0_WP;S2_xm=0.0_WP;S2_ym=0.0_WP;S2_zm=0.0_WP
  
  ! ======================================================================== !
  ! Definition of the stencils                                               !
  !                                                                          !
  !                             SCface                                       !
  !                               |                                          !
  !                               |                                          !
  !   S0_1      S0_2      S0_3    |   xxxx      xxxx                         !
  !   xxxx      S1_2      S1_3    |   S1_4      xxxx                         !
  !   xxxx      xxxx      S2_3    |   S2_4      S2_5                         !
  !                               |                                          !
  !             xxxx      xxxx    |   S0_8      S0_9     S0_10               !
  !             xxxx      S1_7    |   S1_8      S1_9     xxxxx               !
  !             S2_6      S2_7    |   S2_8      xxxx     xxxxx               !
  !                               |                                          !
  !    i-3       i-2       i-1    |    i         i+1      i+2                !
  ! [point 1] [point 2] [point 3] | [point 4] [point 5]             case U>0 !
  !           [point 6] [point 7] | [point 8] [point 9] [point 10]  case U<0 !
  !                               |                                          !
  ! -----------------------------------------------------------------> x,y,z !
  ! ======================================================================== !
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2
        
        ! Along x
        ! S0 +
        call poly3_coeff(x(i-3:i  ),x(i),S0_xp(i,j,-3:-1))
        ! S0 -
        call poly3_coeff(x(i  :i+3),x(i),S0_xm(i,j, 0:+2))
        
        ! S1 +
        call poly3_coeff(x(i-2:i+1),x(i),S1_xp(i,j,-2: 0))
        ! S1 -
        call poly3_coeff(x(i-1:i+2),x(i),S1_xm(i,j,-1:+1))
        
        ! S2 +
        call poly3_coeff(x(i-1:i+2),x(i),S2_xp(i,j,-1:+1))
        ! S2 -
        call poly3_coeff(x(i-2:i+1),x(i),S2_xm(i,j,-2: 0))
        
        ! Along y
        ! S0 +
        call poly3_coeff(y(j-3:j  ),y(j),S0_yp(i,j,-3:-1))
        ! S0 -
        call poly3_coeff(y(j  :j+3),y(j),S0_ym(i,j, 0:+2))
        
        ! S1 +
        call poly3_coeff(y(j-2:j+1),y(j),S1_yp(i,j,-2: 0))
        ! S1 -
        call poly3_coeff(y(j-1:j+2),y(j),S1_ym(i,j,-1:+1))
        
        ! S2 +
        call poly3_coeff(y(j-1:j+2),y(j),S2_yp(i,j,-1:+1))
        ! S2 -
        call poly3_coeff(y(j-2:j+1),y(j),S2_ym(i,j,-2: 0))
        
        ! Along z
        ! S0 +
        S0_zp(-3) =   1.0_WP/3.0_WP
        S0_zp(-2) = - 7.0_WP/6.0_WP
        S0_zp(-1) =  11.0_WP/6.0_WP
        ! S0 -
        S0_zm( 0) =  11.0_WP/6.0_WP
        S0_zm(+1) = - 7.0_WP/6.0_WP
        S0_zm(+2) =   1.0_WP/3.0_WP
        
        ! S1 +
        S1_zp(-2) = - 1.0_WP/6.0_WP
        S1_zp(-1) =   5.0_WP/6.0_WP
        S1_zp( 0) =   1.0_WP/3.0_WP
        ! S1 -
        S1_zm(-1) =   1.0_WP/3.0_WP
        S1_zm( 0) =   5.0_WP/6.0_WP
        S1_zm(+1) = - 1.0_WP/6.0_WP
        
        ! S2 +
        S2_zp(-1) =   1.0_WP/3.0_WP
        S2_zp( 0) =   5.0_WP/6.0_WP
        S2_zp(+1) = - 1.0_WP/6.0_WP
        ! S2 -
        S2_zm(-2) = - 1.0_WP/6.0_WP
        S2_zm(-1) =   5.0_WP/6.0_WP
        S2_zm( 0) =   1.0_WP/3.0_WP
        
        ! Wall BC - x
        if (mask(i-3,j).eq.1) then
           S0_xp(i,j,-2)   = S0_xp(i,j,-2) + S0_xp(i,j,-3)
           S0_xp(i,j,-3)   = 0.0_WP
        end if
        if (mask(i-2,j).eq.1) then
           S0_xp(i,j,:)  = 0.0_WP
           S0_xp(i,j,-1) = 1.0_WP
           S1_xp(i,j,-1) = S1_xp(i,j,-1) + S1_xp(i,j,-2)
           S1_xp(i,j,-2) = 0.0_WP
           S2_xm(i,j,-1) = S2_xm(i,j,-1) + S2_xm(i,j,-2)
           S2_xm(i,j,-2) = 0.0_WP
        end if
        if (mask(i-1,j).eq.1) then
           S0_xp(i,j,:)  = 0.0_WP
           S0_xp(i,j, 0) = 1.0_WP
           S1_xp(i,j,:)  = 0.0_WP
           S1_xp(i,j, 0) = 1.0_WP
           S2_xp(i,j, 0) = S2_xp(i,j,0) + S2_xp(i,j,-1)
           S2_xp(i,j,-1) = 0.0_WP
           S1_xm(i,j, 0) = S1_xm(i,j,0) + S1_xm(i,j,-1)
           S1_xm(i,j,-1) = 0.0_WP
           S2_xm(i,j,:)  = 0.0_WP
           S2_xm(i,j, 0) = 1.0_WP
        end if
        if (mask(i,j).eq.1) then
           S1_xp(i,j,-1) = S1_xp(i,j,-1) + S1_xp(i,j,0)
           S1_xp(i,j, 0) = 0.0_WP
           S2_xp(i,j,:)  = 0.0_WP
           S2_xp(i,j,-1) = 1.0_WP
           S0_xm(i,j,:)  = 0.0_WP
           S0_xm(i,j,-1) = 1.0_WP
           S1_xm(i,j,:)  = 0.0_WP
           S1_xm(i,j,-1) = 1.0_WP
           S2_xm(i,j,-1) = S2_xm(i,j,-1) + S2_xm(i,j,0)
           S2_xm(i,j, 0) = 0.0_WP
        end if
        if (mask(i+1,j).eq.1) then
           S2_xp(i,j, 0) = S2_xp(i,j,0) + S2_xp(i,j,+1)
           S2_xp(i,j,+1) = 0.0_WP
           S0_xm(i,j,:)  = 0.0_WP
           S0_xm(i,j, 0) = 1.0_WP
           S1_xm(i,j, 0) = S1_xm(i,j,0) + S1_xm(i,j,+1)
           S1_xm(i,j,+1) = 0.0_WP
        end if
        if (mask(i+2,j).eq.1) then
           S0_xm(i,j,+1) = S0_xm(i,j,+1) + S0_xm(i,j,+2)
           S0_xm(i,j,+2) = 0.0_WP
        end if
        
        ! Wall BC - y
        if (mask(i,j-3).eq.1) then
           S0_yp(i,j,-2)   = S0_yp(i,j,-2) + S0_yp(i,j,-3)
           S0_yp(i,j,-3)   = 0.0_WP
        end if
        if (mask(i,j-2).eq.1) then
           S0_yp(i,j,:)  = 0.0_WP
           S0_yp(i,j,-1) = 1.0_WP
           S1_yp(i,j,-1) = S1_yp(i,j,-1) + S1_yp(i,j,-2)
           S1_yp(i,j,-2) = 0.0_WP
           S2_ym(i,j,-1) = S2_ym(i,j,-1) + S2_ym(i,j,-2)
           S2_ym(i,j,-2) = 0.0_WP
        end if
        if (mask(i,j-1).eq.1) then
           S0_yp(i,j,:)  = 0.0_WP
           S0_yp(i,j, 0) = 1.0_WP
           S1_yp(i,j,:)  = 0.0_WP
           S1_yp(i,j, 0) = 1.0_WP
           S2_yp(i,j, 0) = S2_yp(i,j,0) + S2_yp(i,j,-1)
           S2_yp(i,j,-1) = 0.0_WP
           S1_ym(i,j, 0) = S1_ym(i,j,0) + S1_ym(i,j,-1)
           S1_ym(i,j,-1) = 0.0_WP
           S2_ym(i,j,:)  = 0.0_WP
           S2_ym(i,j, 0) = 1.0_WP
        end if
        if (mask(i,j).eq.1) then
           S1_yp(i,j,-1) = S1_yp(i,j,-1) + S1_yp(i,j,0)
           S1_yp(i,j, 0) = 0.0_WP
           S2_yp(i,j,:)  = 0.0_WP
           S2_yp(i,j,-1) = 1.0_WP
           S0_ym(i,j,:)  = 0.0_WP
           S0_ym(i,j,-1) = 1.0_WP
           S1_ym(i,j,:)  = 0.0_WP
           S1_ym(i,j,-1) = 1.0_WP
           S2_ym(i,j,-1) = S2_ym(i,j,-1) + S2_ym(i,j,0)
           S2_ym(i,j, 0) = 0.0_WP
        end if
        if (mask(i,j+1).eq.1) then
           S2_yp(i,j, 0) = S2_yp(i,j,0) + S2_yp(i,j,+1)
           S2_yp(i,j,+1) = 0.0_WP
           S0_ym(i,j,:)  = 0.0_WP
           S0_ym(i,j, 0) = 1.0_WP
           S1_ym(i,j, 0) = S1_ym(i,j,0) + S1_ym(i,j,+1)
           S1_ym(i,j,+1) = 0.0_WP
        end if
        if (mask(i,j+2).eq.1) then
           S0_ym(i,j,+1) = S0_ym(i,j,+1) + S0_ym(i,j,+2)
           S0_ym(i,j,+2) = 0.0_WP
        end if
     end do
  end do
  
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
     
     S0_yp(:,jmin,:)  = 0.0_WP
     S0_yp(:,jmin, 0) = 1.0_WP
     S1_yp(:,jmin,:)  = 0.0_WP
     S1_yp(:,jmin, 0) = 1.0_WP
     S2_yp(:,jmin, 0) = S2_yp(:,jmin,0) + S2_yp(:,jmin,-1)
     S2_yp(:,jmin,-1) = 0.0_WP
     S1_ym(:,jmin, 0) = S1_ym(:,jmin,0) + S1_ym(:,jmin,-1)
     S1_ym(:,jmin,-1) = 0.0_WP
     S2_ym(:,jmin,:)  = 0.0_WP
     S2_ym(:,jmin, 0) = 1.0_WP
     
     S0_yp(:,jmin+1,:)  = 0.0_WP
     S0_yp(:,jmin+1,-1) = 1.0_WP
     S1_yp(:,jmin+1,-1) = S1_yp(:,jmin+1,-1) + S1_yp(:,jmin+1,-2)
     S1_yp(:,jmin+1,-2) = 0.0_WP
     S2_ym(:,jmin+1,-1) = S2_ym(:,jmin+1,-1) + S2_ym(:,jmin+1,-2)
     S2_ym(:,jmin+1,-2) = 0.0_WP
     
     S0_yp(:,jmin+2,-2) = S0_yp(:,jmin+2,-2) + S0_yp(:,jmin+2,-3)
     S0_yp(:,jmin+2,-3) = 0.0_WP
     
  end if
  
  if (yper.ne.1 .and. jproc.eq.npy) then
     
     S1_yp(:,jmax+1,-1) = S1_yp(:,jmax+1,-1) + S1_yp(:,jmax+1,0)
     S1_yp(:,jmax+1, 0) = 0.0_WP
     S2_yp(:,jmax+1,:)  = 0.0_WP
     S2_yp(:,jmax+1,-1) = 1.0_WP
     S0_ym(:,jmax+1,:)  = 0.0_WP
     S0_ym(:,jmax+1,-1) = 1.0_WP
     S1_ym(:,jmax+1,:)  = 0.0_WP
     S1_ym(:,jmax+1,-1) = 1.0_WP
     S2_ym(:,jmax+1,-1) = S2_ym(:,jmax+1,-1) + S2_ym(:,jmax+1,0)
     S2_ym(:,jmax+1, 0) = 0.0_WP
     
     S2_yp(:,jmax, 0) = S2_yp(:,jmax,0) + S2_yp(:,jmax,+1)
     S2_yp(:,jmax,+1) = 0.0_WP
     S0_ym(:,jmax,:)  = 0.0_WP
     S0_ym(:,jmax, 0) = 1.0_WP
     S1_ym(:,jmax, 0) = S1_ym(:,jmax,0) + S1_ym(:,jmax,+1)
     S1_ym(:,jmax,+1) = 0.0_WP
     
     S0_ym(:,jmax-1,+1) = S0_ym(:,jmax-1,+1) + S0_ym(:,jmax-1,+2)
     S0_ym(:,jmax-1,+2) = 0.0_WP
     
  end if
  
  ! In z
  ! - Periodic
  !   -> nothing to be done
  
  ! Precompute optimal coefficients
  a0=1.0_WP;a1=6.0_WP;a2=3.0_WP
  do k=kmin_-st1,kmax_+st2
     do j=jmin_-st1,jmax_+st2
        do i=imin_-st1,imax_+st2
           ! Direction x - U>0
           weno5_xp(i,j,k,:) = (a0*S0_xp(i,j,:)+a1*S1_xp(i,j,:)+a2*S2_xp(i,j,:))/(a0+a1+a2)
           ! Direction x - U<0
           weno5_xm(i,j,k,:) = (a0*S0_xm(i,j,:)+a1*S1_xm(i,j,:)+a2*S2_xm(i,j,:))/(a0+a1+a2)
           ! Direction y - V>0
           weno5_yp(i,j,k,:) = (a0*S0_yp(i,j,:)+a1*S1_yp(i,j,:)+a2*S2_yp(i,j,:))/(a0+a1+a2)
           ! Direction y - V<0
           weno5_ym(i,j,k,:) = (a0*S0_ym(i,j,:)+a1*S1_ym(i,j,:)+a2*S2_ym(i,j,:))/(a0+a1+a2)
           ! Direction z - W>0
           weno5_zp(i,j,k,:) = (a0*S0_zp(:)+a1*S1_zp(:)+a2*S2_zp(:))/(a0+a1+a2)
           ! Direction z - W<0
           weno5_zm(i,j,k,:) = (a0*S0_zm(:)+a1*S1_zm(:)+a2*S2_zm(:))/(a0+a1+a2)
        end do
     end do
  end do
  
  return
end subroutine scalar_weno5_init


! ==================================================== !
! Compute the weno 5 coefficients for the given scalar !
! ==================================================== !
subroutine scalar_weno5_coeff(isc)
  use scalar_weno5
  implicit none
  
  integer, intent(in) :: isc
  integer  :: i,j,k
  real(WP) :: r0,r1,r2,a0,a1,a2
  real(WP), parameter :: epsilon = 1.0e-6_WP
  
  ! If ultimate
  if (ult) return

  do k=kmin_-st1,kmax_+st2
     do j=jmin_-st1,jmax_+st2
        do i=imin_-st1,imax_+st2
           ! Direction x - U>0
           r0 = 13.0_WP*(SC(i-3,j,k,isc)-2.0_WP*SC(i-2,j,k,isc)+SC(i-1,j,k,isc))**2+3.0_WP*(SC(i-3,j,k,isc)-4.0_WP*SC(i-2,j,k,isc)+3.0_WP*SC(i-1,j,k,isc))**2
           r1 = 13.0_WP*(SC(i-2,j,k,isc)-2.0_WP*SC(i-1,j,k,isc)+SC(i  ,j,k,isc))**2+3.0_WP*(SC(i-2,j,k,isc)-SC(i  ,j,k,isc))**2
           r2 = 13.0_WP*(SC(i-1,j,k,isc)-2.0_WP*SC(i  ,j,k,isc)+SC(i+1,j,k,isc))**2+3.0_WP*(3.0_WP*SC(i-1,j,k,isc)-4.0_WP*SC(i  ,j,k,isc)+SC(i+1,j,k,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_xp(i,j,k,:) = (a0*S0_xp(i,j,:)+a1*S1_xp(i,j,:)+a2*S2_xp(i,j,:))/(a0+a1+a2)
           
           ! Direction x - U<0
           r0 = 13.0_WP*(SC(i+2,j,k,isc)-2.0_WP*SC(i+1,j,k,isc)+SC(i  ,j,k,isc))**2+3.0_WP*(SC(i+2,j,k,isc)-4.0_WP*SC(i+1,j,k,isc)+3.0_WP*SC(i  ,j,k,isc))**2
           r1 = 13.0_WP*(SC(i+1,j,k,isc)-2.0_WP*SC(i  ,j,k,isc)+SC(i-1,j,k,isc))**2+3.0_WP*(SC(i+1,j,k,isc)-SC(i-1,j,k,isc))**2
           r2 = 13.0_WP*(SC(i  ,j,k,isc)-2.0_WP*SC(i-1,j,k,isc)+SC(i-2,j,k,isc))**2+3.0_WP*(3.0_WP*SC(i  ,j,k,isc)-4.0_WP*SC(i-1,j,k,isc)+SC(i-2,j,k,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_xm(i,j,k,:)= (a0*S0_xm(i,j,:)+a1*S1_xm(i,j,:)+a2*S2_xm(i,j,:))/(a0+a1+a2)
           
           ! Direction y - V>0
           r0 = 13.0_WP*(SC(i,j-3,k,isc)-2.0_WP*SC(i,j-2,k,isc)+SC(i,j-1,k,isc))**2+3.0_WP*(SC(i,j-3,k,isc)-4.0_WP*SC(i,j-2,k,isc)+3.0_WP*SC(i,j-1,k,isc))**2
           r1 = 13.0_WP*(SC(i,j-2,k,isc)-2.0_WP*SC(i,j-1,k,isc)+SC(i,j  ,k,isc))**2+3.0_WP*(SC(i,j-2,k,isc)-SC(i,j  ,k,isc))**2
           r2 = 13.0_WP*(SC(i,j-1,k,isc)-2.0_WP*SC(i,j  ,k,isc)+SC(i,j+1,k,isc))**2+3.0_WP*(3.0_WP*SC(i,j-1,k,isc)-4.0_WP*SC(i,j  ,k,isc)+SC(i,j+1,k,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_yp(i,j,k,:) = (a0*S0_yp(i,j,:)+a1*S1_yp(i,j,:)+a2*S2_yp(i,j,:))/(a0+a1+a2)
           
           ! Direction y - V<0
           r0 = 13.0_WP*(SC(i,j+2,k,isc)-2.0_WP*SC(i,j+1,k,isc)+SC(i,j  ,k,isc))**2+3.0_WP*(SC(i,j+2,k,isc)-4.0_WP*SC(i,j+1,k,isc)+3.0_WP*SC(i,j  ,k,isc))**2
           r1 = 13.0_WP*(SC(i,j+1,k,isc)-2.0_WP*SC(i,j  ,k,isc)+SC(i,j-1,k,isc))**2+3.0_WP*(SC(i,j+1,k,isc)-SC(i,j-1,k,isc))**2
           r2 = 13.0_WP*(SC(i,j  ,k,isc)-2.0_WP*SC(i,j-1,k,isc)+SC(i,j-2,k,isc))**2+3.0_WP*(3.0_WP*SC(i,j  ,k,isc)-4.0_WP*SC(i,j-1,k,isc)+SC(i,j-2,k,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_ym(i,j,k,:)= (a0*S0_ym(i,j,:)+a1*S1_ym(i,j,:)+a2*S2_ym(i,j,:))/(a0+a1+a2)
           
           ! Direction z - W>0
           r0 = 13.0_WP*(SC(i,j,k-3,isc)-2.0_WP*SC(i,j,k-2,isc)+SC(i,j,k-1,isc))**2+3.0_WP*(SC(i,j,k-3,isc)-4.0_WP*SC(i,j,k-2,isc)+3.0_WP*SC(i,j,k-1,isc))**2
           r1 = 13.0_WP*(SC(i,j,k-2,isc)-2.0_WP*SC(i,j,k-1,isc)+SC(i,j,k  ,isc))**2+3.0_WP*(SC(i,j,k-2,isc)-SC(i,j,k  ,isc))**2
           r2 = 13.0_WP*(SC(i,j,k-1,isc)-2.0_WP*SC(i,j,k  ,isc)+SC(i,j,k+1,isc))**2+3.0_WP*(3.0_WP*SC(i,j,k-1,isc)-4.0_WP*SC(i,j,k  ,isc)+SC(i,j,k+1,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_zp(i,j,k,:) = (a0*S0_zp(:)+a1*S1_zp(:)+a2*S2_zp(:))/(a0+a1+a2)
           
           ! Direction z - W<0
           r0 = 13.0_WP*(SC(i,j,k+2,isc)-2.0_WP*SC(i,j,k+1,isc)+SC(i,j,k  ,isc))**2+3.0_WP*(SC(i,j,k+2,isc)-4.0_WP*SC(i,j,k+1,isc)+3.0_WP*SC(i,j,k  ,isc))**2
           r1 = 13.0_WP*(SC(i,j,k+1,isc)-2.0_WP*SC(i,j,k  ,isc)+SC(i,j,k-1,isc))**2+3.0_WP*(SC(i,j,k+1,isc)-SC(i,j,k-1,isc))**2
           r2 = 13.0_WP*(SC(i,j,k  ,isc)-2.0_WP*SC(i,j,k-1,isc)+SC(i,j,k-2,isc))**2+3.0_WP*(3.0_WP*SC(i,j,k  ,isc)-4.0_WP*SC(i,j,k-1,isc)+SC(i,j,k-2,isc))**2
           a0 = 1.0_WP/(epsilon+r0)**2
           a1 = 6.0_WP/(epsilon+r1)**2
           a2 = 3.0_WP/(epsilon+r2)**2
           weno5_zm(i,j,k,:)= (a0*S0_zm(:)+a1*S1_zm(:)+a2*S2_zm(:))/(a0+a1+a2)
        end do
     end do
  end do
  
  return
end subroutine scalar_weno5_coeff



! =========================================================== !
! Compute the residuals of the scalar equations               !
!                                                             !
! - velocity field n+1 stored in U/rhoU                       !
! - scalar field n+1 stored in SCmid                          !
!                                                             !
! 3 working arrays of size (at least) (nx_+1)*(ny_+1)*(nz_+1) !
! =========================================================== !
subroutine scalar_weno5_residual
  use scalar_weno5
  use parallel
  use memory
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: rhs
  
  do isc=1,nscalar
     
     ! Get coefficients first
     call scalar_weno5_coeff(isc)
     
     ! Now compute residuals
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              
              FX(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoUt(i,j,k,isc)+abs(rhoUt(i,j,k,isc)))*sum(weno5_xp(i,j,k,:)*SC(i-3:i+1,j,k,isc)) &
                   - 0.5_WP*(rhoUt(i,j,k,isc)-abs(rhoUt(i,j,k,isc)))*sum(weno5_xm(i,j,k,:)*SC(i-2:i+2,j,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_x(i,j,:)*DIFF(i-st2:i+st1,j,k,isc)) * &
                     sum(grad_x(i,j,:)*SC(i-st2:i+st1,j,k,isc))
              
              FY(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoVt(i,j,k,isc)+abs(rhoVt(i,j,k,isc)))*sum(weno5_yp(i,j,k,:)*SC(i,j-3:j+1,k,isc)) &
                   - 0.5_WP*(rhoVt(i,j,k,isc)-abs(rhoVt(i,j,k,isc)))*sum(weno5_ym(i,j,k,:)*SC(i,j-2:j+2,k,isc)) &
                   ! Viscous term
                   + sum(interp_sc_y(i,j,:)*DIFF(i,j-st2:j+st1,k,isc)) * &
                     sum(grad_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
              
              FZ(i,j,k) = &
                   ! Convective term
                   - 0.5_WP*(rhoWt(i,j,k,isc)+abs(rhoWt(i,j,k,isc)))*sum(weno5_zp(i,j,k,:)*SC(i,j,k-3:k+1,isc)) &
                   - 0.5_WP*(rhoWt(i,j,k,isc)-abs(rhoWt(i,j,k,isc)))*sum(weno5_zm(i,j,k,:)*SC(i,j,k-2:k+2,isc)) &
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
end subroutine scalar_weno5_residual


! =========================================================== !
! Inverse the linear system obtained from the implicit scalar !
! transport equation                                          !
! =========================================================== !
subroutine scalar_weno5_inverse
  use scalar_weno5
  use parallel
  use memory
  use implicit
  use time_info
  implicit none
  
  integer  :: i,j,k,isc
  real(WP) :: conv1,conv2,conv3,conv4
  real(WP) :: visc1,visc2
  real(WP) :: dt2
  
  dt2 = dt/2.0_WP

  ! If purely explicit return
  if ((.not.implicit_any)) return
  
  do isc=1,nscalar
     
     ! Get coefficients first
     call scalar_weno5_coeff(isc)
     
     ! X-direction
     if (implicit_x) then
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 
                 visc1 = sum(interp_sc_x(i  ,j,:)*DIFF(i  -st2:i  +st1,j,k,isc))
                 visc2 = sum(interp_sc_x(i+1,j,:)*DIFF(i+1-st2:i+1+st1,j,k,isc))
                 
                 Ax(j,k,i,:) = 0.0_WP

                 Ax(j,k,i,-1) = - dt2 * div_u(i,j,0)*visc1*grad_x(i,j,-1)
                 
                 Ax(j,k,i, 0) = RHO_s(i,j,k) - dt2*(div_u(i,j,0)*visc1*grad_x(i,j,0) + div_u(i,j,1)*visc2*grad_x(i+1,j,-1))
                 
                 Ax(j,k,i,+1) = - dt2 * div_u(i,j,1)*visc2*grad_x(i+1,j,0)
                 
                 ! JFM 6/13/14
                 ! Use RHO_s, since that's the density at the time to which we're updating
                 Rx(j,k,i) = RHO_s(i,j,k) * ResSC(i,j,k,isc)
                 
              end do
           end do
        end do
     
        if (npx.eq.1) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    
                    conv1 = 0.5_WP*(rhoUt(i  ,j,k,isc)+abs(rhoUt(i  ,j,k,isc)))
                    conv2 = 0.5_WP*(rhoUt(i+1,j,k,isc)+abs(rhoUt(i+1,j,k,isc)))
                    conv3 = 0.5_WP*(rhoUt(i  ,j,k,isc)-abs(rhoUt(i  ,j,k,isc)))
                    conv4 = 0.5_WP*(rhoUt(i+1,j,k,isc)-abs(rhoUt(i+1,j,k,isc)))
                    
                    Ax(j,k,i,-3) = Ax(j,k,i,-3) + dt2*(div_u(i,j,0)*(conv1*weno5_xp(i,j,k,-3)                         )                                                                     )
                    Ax(j,k,i,-2) = Ax(j,k,i,-2) + dt2*(div_u(i,j,0)*(conv1*weno5_xp(i,j,k,-2)+conv3*weno5_xm(i,j,k,-2))+div_u(i,j,1)*(conv2*weno5_xp(i+1,j,k,-3)                           ))
                    Ax(j,k,i,-1) = Ax(j,k,i,-1) + dt2*(div_u(i,j,0)*(conv1*weno5_xp(i,j,k,-1)+conv3*weno5_xm(i,j,k,-1))+div_u(i,j,1)*(conv2*weno5_xp(i+1,j,k,-2)+conv4*weno5_xm(i+1,j,k,-2)))
                    Ax(j,k,i, 0) = Ax(j,k,i, 0) + dt2*(div_u(i,j,0)*(conv1*weno5_xp(i,j,k, 0)+conv3*weno5_xm(i,j,k, 0))+div_u(i,j,1)*(conv2*weno5_xp(i+1,j,k,-1)+conv4*weno5_xm(i+1,j,k,-1)))
                    Ax(j,k,i,+1) = Ax(j,k,i,+1) + dt2*(div_u(i,j,0)*(conv1*weno5_xp(i,j,k,+1)+conv3*weno5_xm(i,j,k,+1))+div_u(i,j,1)*(conv2*weno5_xp(i+1,j,k, 0)+conv4*weno5_xm(i+1,j,k, 0)))
                    Ax(j,k,i,+2) = Ax(j,k,i,+2) + dt2*(div_u(i,j,0)*(                         conv3*weno5_xm(i,j,k,+2))+div_u(i,j,1)*(conv2*weno5_xp(i+1,j,k,+1)+conv4*weno5_xm(i+1,j,k,+1)))
                    Ax(j,k,i,+3) = Ax(j,k,i,+3) + dt2*(                                                                 div_u(i,j,1)*(                           conv4*weno5_xm(i+1,j,k,+2)))
                    
                 end do
              end do
           end do
           call implicit_solve_x(7)
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
                           
                 Ay(i,k,j, 0) = RHO_s(i,j,k) - dt2*(div_v(i,j,0)*visc1*grad_y(i,j,0) + div_v(i,j,1)*visc2*grad_y(i,j+1,-1))
                 
                 Ay(i,k,j,+1) = - dt2 * div_v(i,j,1)*visc2*grad_y(i,j+1,0)
                 
                 Ry(i,k,j) = RHO_s(i,j,k) * Rx(j,k,i)
                 
              end do
           end do
        end do

        if (npy.eq.1) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    
                    conv1 = 0.5_WP*(rhoVt(i,j  ,k,isc)+abs(rhoVt(i,j  ,k,isc)))
                    conv2 = 0.5_WP*(rhoVt(i,j+1,k,isc)+abs(rhoVt(i,j+1,k,isc)))
                    conv3 = 0.5_WP*(rhoVt(i,j  ,k,isc)-abs(rhoVt(i,j  ,k,isc)))
                    conv4 = 0.5_WP*(rhoVt(i,j+1,k,isc)-abs(rhoVt(i,j+1,k,isc)))
                    
                    Ay(i,k,j,-3) = Ay(i,k,j,-3) + dt2*(div_v(i,j,0)*(conv1*weno5_yp(i,j,k,-3)                         )                                                                     )
                    Ay(i,k,j,-2) = Ay(i,k,j,-2) + dt2*(div_v(i,j,0)*(conv1*weno5_yp(i,j,k,-2)+conv3*weno5_ym(i,j,k,-2))+div_v(i,j,1)*(conv2*weno5_yp(i,j+1,k,-3)                           ))
                    Ay(i,k,j,-1) = Ay(i,k,j,-1) + dt2*(div_v(i,j,0)*(conv1*weno5_yp(i,j,k,-1)+conv3*weno5_ym(i,j,k,-1))+div_v(i,j,1)*(conv2*weno5_yp(i,j+1,k,-2)+conv4*weno5_ym(i,j+1,k,-2)))
                    Ay(i,k,j, 0) = Ay(i,k,j, 0) + dt2*(div_v(i,j,0)*(conv1*weno5_yp(i,j,k, 0)+conv3*weno5_ym(i,j,k, 0))+div_v(i,j,1)*(conv2*weno5_yp(i,j+1,k,-1)+conv4*weno5_ym(i,j+1,k,-1)))
                    Ay(i,k,j,+1) = Ay(i,k,j,+1) + dt2*(div_v(i,j,0)*(conv1*weno5_yp(i,j,k,+1)+conv3*weno5_ym(i,j,k,+1))+div_v(i,j,1)*(conv2*weno5_yp(i,j+1,k, 0)+conv4*weno5_ym(i,j+1,k, 0)))
                    Ay(i,k,j,+2) = Ay(i,k,j,+2) + dt2*(div_v(i,j,0)*(                         conv3*weno5_ym(i,j,k,+2))+div_v(i,j,1)*(conv2*weno5_yp(i,j+1,k,+1)+conv4*weno5_ym(i,j+1,k,+1)))
                    Ay(i,k,j,+3) = Ay(i,k,j,+3) + dt2*(                                                                 div_v(i,j,1)*(                           conv4*weno5_ym(i,j+1,k,+2)))
                    
                 end do
              end do
           end do
           if (icyl.eq.1 .and. jproc.eq.1) then
              do k=kmin_,kmax_
                 Ay(:,k,jmin, 0) = Ay(:,k,jmin,0) + Ay(:,k,jmin,-1) + Ay(:,k,jmin,-2) + Ay(:,k,jmin,-3)
                 Ay(:,k,jmin,-1) = 0.0_WP
                 Ay(:,k,jmin,-2) = 0.0_WP
                 Ay(:,k,jmin,-3) = 0.0_WP
                 
                 Ay(:,k,jmin+1,-1) = Ay(:,k,jmin+1,-1) + Ay(:,k,jmin+1,-2) + Ay(:,k,jmin+1,-3)
                 Ay(:,k,jmin+1,-2) = 0.0_WP
                 Ay(:,k,jmin+1,-3) = 0.0_WP
                 
                 Ay(:,k,jmin+2,-2) = Ay(:,k,jmin+2,-2) + Ay(:,k,jmin+2,-3)
                 Ay(:,k,jmin+2,-3) = 0.0_WP
              end do
           end if
           call implicit_solve_y(7)
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
                 
                 Az(i,j,k, 0) = RHO_s(i,j,k) - dt2*(div_w(i,j,0)*visc1*grad_z(i,j,0) + div_w(i,j,1)*visc2*grad_z(i,j,-1))
                 
                 Az(i,j,k,+1) = - dt2 * div_w(i,j,1)*visc2*grad_z(i,j,0)
                 
                 Rz(i,j,k) = RHO_s(i,j,k) * Ry(i,k,j)
                 
              end do
           end do
        end do
        
        if (npz.eq.1) then
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    
                    conv1 = 0.5_WP*(rhoWt(i,j,k  ,isc)+abs(rhoWt(i,j,k  ,isc)))
                    conv2 = 0.5_WP*(rhoWt(i,j,k+1,isc)+abs(rhoWt(i,j,k+1,isc)))
                    conv3 = 0.5_WP*(rhoWt(i,j,k  ,isc)-abs(rhoWt(i,j,k  ,isc)))
                    conv4 = 0.5_WP*(rhoWt(i,j,k+1,isc)-abs(rhoWt(i,j,k+1,isc)))
                    
                    Az(i,j,k,-3) = Az(i,j,k,-3) + dt2*(div_w(i,j,0)*(conv1*weno5_zp(i,j,k,-3)                         )                                                                     )
                    Az(i,j,k,-2) = Az(i,j,k,-2) + dt2*(div_w(i,j,0)*(conv1*weno5_zp(i,j,k,-2)+conv3*weno5_zm(i,j,k,-2))+div_w(i,j,1)*(conv2*weno5_zp(i,j,k+1,-3)                           ))
                    Az(i,j,k,-1) = Az(i,j,k,-1) + dt2*(div_w(i,j,0)*(conv1*weno5_zp(i,j,k,-1)+conv3*weno5_zm(i,j,k,-1))+div_w(i,j,1)*(conv2*weno5_zp(i,j,k+1,-2)+conv4*weno5_zm(i,j,k+1,-2)))
                    Az(i,j,k, 0) = Az(i,j,k, 0) + dt2*(div_w(i,j,0)*(conv1*weno5_zp(i,j,k, 0)+conv3*weno5_zm(i,j,k, 0))+div_w(i,j,1)*(conv2*weno5_zp(i,j,k+1,-1)+conv4*weno5_zm(i,j,k+1,-1)))
                    Az(i,j,k,+1) = Az(i,j,k,+1) + dt2*(div_w(i,j,0)*(conv1*weno5_zp(i,j,k,+1)+conv3*weno5_zm(i,j,k,+1))+div_w(i,j,1)*(conv2*weno5_zp(i,j,k+1, 0)+conv4*weno5_zm(i,j,k+1, 0)))
                    Az(i,j,k,+2) = Az(i,j,k,+2) + dt2*(div_w(i,j,0)*(                         conv3*weno5_zm(i,j,k,+2))+div_w(i,j,1)*(conv2*weno5_zp(i,j,k+1,+1)+conv4*weno5_zm(i,j,k+1,+1)))
                    Az(i,j,k,+3) = Az(i,j,k,+3) + dt2*(                                                                 div_w(i,j,1)*(                           conv4*weno5_zm(i,j,k+1,+2)))
                    
                 end do
              end do
           end do
           call implicit_solve_z(7)
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
end subroutine scalar_weno5_inverse
