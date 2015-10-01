module bicgstab2
  use pressure
  implicit none
  
  real(WP), dimension(:,:,:), pointer :: res,res0
  real(WP), dimension(:,:,:), pointer :: r_,r_hat
  real(WP), dimension(:,:,:), pointer :: s_,s_hat
  real(WP), dimension(:,:,:), pointer :: u_,u_hat
  real(WP), dimension(:,:,:), pointer :: v_,v_hat
  real(WP), dimension(:,:,:), pointer :: y_,y_hat
  real(WP), dimension(:,:,:), pointer :: t_,x_,w_
  
end module bicgstab2


! ================================= !
! Initialize the BiCGStab(2) module !
! ================================= !
subroutine bicgstab2_init
  use bicgstab2
  use partition
  implicit none
  
  ! Allocation of work arrays - internal
  allocate(res (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(res0(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(r_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(s_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(t_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(u_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(v_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(w_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(y_  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  
  ! Allocation of work arrays - external
  allocate(r_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(s_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(u_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(v_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(x_   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(y_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  return
end subroutine bicgstab2_init


! ========================= !
! BiCGStab(2) Linear solver !
! ========================= !
subroutine bicgstab2_solve
  use bicgstab2
  use parallel
  implicit none
  
  logical  :: done
  real(WP) :: rho0,alpha,rho1,gamma,omega1,omega2,beta
  real(WP) :: mu,nu,tau
  real(WP) :: max_res,max_res0
  integer  :: iter
  
  ! Initialization
  rho0   = 1.0_WP
  alpha  = 0.0_WP
  omega2 = 1.0_WP
  u_     = 0.0_WP
  
  ! Compute initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  done = .false.
  iter = 0
  call parallel_max(maxval(abs(res)),max_res0)
  if (max_res0.eq.0.0_WP) return
  
  ! iterate
  do while(.not.done)
     
     iter = iter + 1
     
     rho0 = -omega2*rho0
     
     ! Even BiCG step ------------------------------
     
     ! rho1 = (r0,res)
     call parallel_sum(sum(res*res0),rho1)
     beta = alpha*rho1/rho0
     rho0 = rho1
     
     ! u = res - beta*u
     ! u_hat = D^-1 u
     u_ = res - beta*u_
     call bicgstab2_apply_precond(u_,u_hat)
     
     ! v = A*u_hat
     call pressure_operator(u_hat,v_)
     
     ! gamma = (v,r0)
     call parallel_sum(sum(v_*res0),gamma)
     alpha = rho0/gamma
     
     ! r = res - alpha*v
     ! r_hat = D^-1 r
     r_ = res - alpha*v_
     call bicgstab2_apply_precond(r_,r_hat)
     
     ! s = A*r_hat
     call pressure_operator(r_hat,s_)
     
     ! x = scal + alpha*u_hat
     x_= DP + alpha*u_hat
     
     ! Odd BiCG step -------------------------------
     
     ! rho1 = (r0,s)
     call parallel_sum(sum(s_*res0),rho1)
     beta = alpha*rho1/rho0
     rho0 = rho1
     
     ! v = s - beta*v
     ! v_hat = D^-1 v
     v_ = s_ - beta*v_
     call bicgstab2_apply_precond(v_,v_hat)
     
     ! w = A*v_hat
     call pressure_operator(v_hat,w_)
     
     ! gamma = (w,r0)
     call parallel_sum(sum(w_*res0),gamma)
     alpha = rho0/gamma
     
     ! u = r - beta*u
     ! r = r - alpha*v
     ! s = s - alpha*w
     ! s_hat = D^-1 s
     u_ = r_ -  beta*u_
     r_ = r_ - alpha*v_
     s_ = s_ - alpha*w_
     call bicgstab2_apply_precond(s_,s_hat)
     
     ! t = A*s_hat
     call pressure_operator(s_hat,t_)
     
     ! GCR(2)-part ---------------------------------
     
     ! omega1 = (r,s)
     ! mu = (s,s)
     ! nu = (s,t)
     ! tau = (t,t)
     ! omega2 = (r,t)
     call parallel_sum(sum(r_*s_),omega1)
     call parallel_sum(sum(s_*s_),mu)
     call parallel_sum(sum(s_*t_),nu)
     call parallel_sum(sum(t_*t_),tau)
     call parallel_sum(sum(r_*t_),omega2)
     tau = tau - nu**2/mu
     omega2 = (omega2-nu*omega1/mu)/tau
     omega1 = (omega1-nu*omega2)/mu
     
     ! scal = x + omega1*r + omega2*s + alpha*u
     ! res = r - omega1*s - omega2*t
     y_ = omega1*r_ + omega2*s_ + alpha*u_
     call bicgstab2_apply_precond(y_,y_hat)
     DP  = x_ + y_hat
     res = r_ - omega1*s_ - omega2*t_
     
     ! Check if converged
     call parallel_max(maxval(abs(res)),max_res)
     max_res = max_res/max_res0
     if (max_res.lt.cvg .or. iter.ge.600) done=.true.
     
     if (.not.done) then
        ! u = u - omega1*v - omega2*w
        u_ = u_ - omega1*v_ - omega2*w_
     end if
     
  end do
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab2_solve


! ================================================= !
! Select wich operator to use to get the residuals  !
! Compute the product of vec by the matrix K: K.vec !
! ================================================= !
subroutine bicgstab2_apply_precond(vec,Kvec)
  use bicgstab
  use partition
  implicit none

  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: vec
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: Kvec
  
  select case(trim(pressure_precond))
  case('none')
     Kvec = 0.0_WP
     Kvec(imin_:imax_,jmin_:jmax_,kmin_:kmax_) = vec
  case('diag')
     call pressure_precond_diag(vec,Kvec)
  case('tridiag')
     call pressure_precond_tridiag(vec,Kvec)
  case('icc')
     call pressure_precond_icc(vec,Kvec)
  case default
     call pressure_precond_icc(vec,Kvec)
     !call pressure_precond_diag(vec,Kvec)
     !call die('bicgstab_apply_precond: unknown preconditioner')
  end select
  
  return
end subroutine bicgstab2_apply_precond
