module bicgstab
  use pressure
  implicit none

  ! Maximum residual
  real(WP) :: max_res,max_res0
  ! Iteration number
  integer  :: iter
  
  ! Work variables
  ! --------------
  
  ! RP = B
  ! DP = X
  ! res = B-AX
  real(WP), dimension(:,:,:), pointer :: res,res0,pp,s,vv,t,p_hat,s_hat
  
  ! For 3D bicgstab
  real(WP) :: alpha,rho1,rho2,gamma,omega,beta
  
  ! For bicgstab in X
  real(WP), dimension(:),   pointer :: alpha_x, rho1_x, rho2_x, gamma_x, omega_x, beta_x, max_res_x, max_res0_x, buf1_x, buf2_x
  logical,  dimension(:),   pointer :: done_x
  
  ! For bicgstab in Y
  real(WP), dimension(:),   pointer :: alpha_y, rho1_y, rho2_y, gamma_y, omega_y, beta_y, max_res_y, max_res0_y, buf1_y, buf2_y
  logical,  dimension(:),   pointer :: done_y
  
  ! For bicgstab in Z
  real(WP), dimension(:),   pointer :: alpha_z, rho1_z, rho2_z, gamma_z, omega_z, beta_z, max_res_z, max_res0_z, buf1_z, buf2_z
  logical,  dimension(:),   pointer :: done_z
  
  ! For bicgstab in YZ
  real(WP), dimension(:,:), pointer :: alpha_yz,rho1_yz,rho2_yz,gamma_yz,omega_yz,beta_yz,max_res_yz,max_res0_yz,buf1_yz,buf2_yz
  logical,  dimension(:,:), pointer :: done_yz

  ! For bicgstab in XZ
  real(WP), dimension(:,:), pointer :: alpha_xz,rho1_xz,rho2_xz,gamma_xz,omega_xz,beta_xz,max_res_xz,max_res0_xz,buf1_xz,buf2_xz
  logical,  dimension(:,:), pointer :: done_xz
  
  ! For bicgstab in XY
  real(WP), dimension(:,:), pointer :: alpha_xy,rho1_xy,rho2_xy,gamma_xy,omega_xy,beta_xy,max_res_xy,max_res0_xy,buf1_xy,buf2_xy
  logical,  dimension(:,:), pointer :: done_xy
  
end module bicgstab


! ============================== !
! Initialize the BiCGStab module !
! ============================== !
subroutine bicgstab_init
  use bicgstab
  use partition
  implicit none
  
  ! Allocation of work arrays
  allocate(res (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(res0(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(pp  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(s   (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(vv  (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(t   (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(p_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(s_hat(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Allocation of direction-dependent work variables
  if      (     fft_x .and. .not.fft_y .and. .not.fft_z) then
     allocate(alpha_x   (imin_:imax_))
     allocate(rho1_x    (imin_:imax_))
     allocate(rho2_x    (imin_:imax_))
     allocate(gamma_x   (imin_:imax_))
     allocate(omega_x   (imin_:imax_))
     allocate(beta_x    (imin_:imax_))
     allocate(max_res0_x(imin_:imax_))
     allocate(max_res_x (imin_:imax_))
     allocate(buf1_x    (imin_:imax_))
     allocate(buf2_x    (imin_:imax_))
     allocate(done_x    (imin_:imax_))
  else if (.not.fft_x .and.      fft_y .and. .not.fft_z) then
     allocate(alpha_y   (jmin_:jmax_))
     allocate(rho1_y    (jmin_:jmax_))
     allocate(rho2_y    (jmin_:jmax_))
     allocate(gamma_y   (jmin_:jmax_))
     allocate(omega_y   (jmin_:jmax_))
     allocate(beta_y    (jmin_:jmax_))
     allocate(max_res0_y(jmin_:jmax_))
     allocate(max_res_y (jmin_:jmax_))
     allocate(buf1_y    (jmin_:jmax_))
     allocate(buf2_y    (jmin_:jmax_))
     allocate(done_y    (jmin_:jmax_))
  else if (.not.fft_x .and. .not.fft_y .and.      fft_z) then
     allocate(alpha_z   (kmin_:kmax_))
     allocate(rho1_z    (kmin_:kmax_))
     allocate(rho2_z    (kmin_:kmax_))
     allocate(gamma_z   (kmin_:kmax_))
     allocate(omega_z   (kmin_:kmax_))
     allocate(beta_z    (kmin_:kmax_))
     allocate(max_res0_z(kmin_:kmax_))
     allocate(max_res_z (kmin_:kmax_))
     allocate(buf1_z    (kmin_:kmax_))
     allocate(buf2_z    (kmin_:kmax_))
     allocate(done_z    (kmin_:kmax_))
  else if (.not.fft_x .and.      fft_y .and.      fft_z) then
     allocate(alpha_yz   (jmin_:jmax_,kmin_:kmax_))
     allocate(rho1_yz    (jmin_:jmax_,kmin_:kmax_))
     allocate(rho2_yz    (jmin_:jmax_,kmin_:kmax_))
     allocate(gamma_yz   (jmin_:jmax_,kmin_:kmax_))
     allocate(omega_yz   (jmin_:jmax_,kmin_:kmax_))
     allocate(beta_yz    (jmin_:jmax_,kmin_:kmax_))
     allocate(max_res0_yz(jmin_:jmax_,kmin_:kmax_))
     allocate(max_res_yz (jmin_:jmax_,kmin_:kmax_))
     allocate(buf1_yz    (jmin_:jmax_,kmin_:kmax_))
     allocate(buf2_yz    (jmin_:jmax_,kmin_:kmax_))
     allocate(done_yz    (jmin_:jmax_,kmin_:kmax_))
  else if (     fft_x .and. .not.fft_y .and.      fft_z) then
     allocate(alpha_xz   (imin_:imax_,kmin_:kmax_))
     allocate(rho1_xz    (imin_:imax_,kmin_:kmax_))
     allocate(rho2_xz    (imin_:imax_,kmin_:kmax_))
     allocate(gamma_xz   (imin_:imax_,kmin_:kmax_))
     allocate(omega_xz   (imin_:imax_,kmin_:kmax_))
     allocate(beta_xz    (imin_:imax_,kmin_:kmax_))
     allocate(max_res0_xz(imin_:imax_,kmin_:kmax_))
     allocate(max_res_xz (imin_:imax_,kmin_:kmax_))
     allocate(buf1_xz    (imin_:imax_,kmin_:kmax_))
     allocate(buf2_xz    (imin_:imax_,kmin_:kmax_))
     allocate(done_xz    (imin_:imax_,kmin_:kmax_))
  else if (     fft_x .and.      fft_y .and. .not.fft_z) then
     allocate(alpha_xy   (imin_:imax_,jmin_:jmax_))
     allocate(rho1_xy    (imin_:imax_,jmin_:jmax_))
     allocate(rho2_xy    (imin_:imax_,jmin_:jmax_))
     allocate(gamma_xy   (imin_:imax_,jmin_:jmax_))
     allocate(omega_xy   (imin_:imax_,jmin_:jmax_))
     allocate(beta_xy    (imin_:imax_,jmin_:jmax_))
     allocate(max_res0_xy(imin_:imax_,jmin_:jmax_))
     allocate(max_res_xy (imin_:imax_,jmin_:jmax_))
     allocate(buf1_xy    (imin_:imax_,jmin_:jmax_))
     allocate(buf2_xy    (imin_:imax_,jmin_:jmax_))
     allocate(done_xy    (imin_:imax_,jmin_:jmax_))
  else if (.not.fft_x .and. .not.fft_y .and. .not.fft_z) then
     ! Nothing to do, only scalars
  else if (     fft_x .and.      fft_y .and.      fft_z) then
     ! Nothing to do, only scalars
  else
     call die('bicgstab_init: directional bicgstab solver not yet implemented.')
  end if
  
  return
end subroutine bicgstab_init


! ================================ !
! BiCGStab Linear solver           !
! Bi-Conjugate Gradient Stabilized !
!                                  !
! Solves Ax=b                      !
!                                  !
! RP = B                           !
! DP = X                           !
! res = B-AX                       !
! ================================ !
subroutine bicgstab_solve
  use bicgstab
  use parallel
  implicit none
  
  if      (     fft_x .and. .not.fft_y .and. .not.fft_z) then
     call bicgstab_solve_x
  else if (.not.fft_x .and.      fft_y .and. .not.fft_z) then
     call bicgstab_solve_y
  else if (.not.fft_x .and. .not.fft_y .and.      fft_z) then
     call bicgstab_solve_z
  else if (.not.fft_x .and.      fft_y .and.      fft_z) then
     call bicgstab_solve_yz
  else if (     fft_x .and. .not.fft_y .and.      fft_z) then
     call bicgstab_solve_xz
  else if (     fft_x .and.      fft_y .and. .not.fft_z) then
     call bicgstab_solve_xy
  else if (.not.fft_x .and. .not.fft_y .and. .not.fft_z) then
     call bicgstab_solve_3d
  else if (     fft_x .and.      fft_y .and.      fft_z) then
     call pressure_precond_diag(RP,DP)
  end if
  
  return
end subroutine bicgstab_solve


! ================================================= !
! Select wich operator to use to get the residuals  !
! Compute the product of vec by the matrix K: K.vec !
! ================================================= !
subroutine bicgstab_apply_precond(vec,Kvec)
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
     call die('bicgstab_apply_precond: unknown preconditioner')
  end select
  
  return
end subroutine bicgstab_apply_precond


! ================================== !
! Solver with no periodic directions !
! ================================== !
subroutine bicgstab_solve_3d
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer  :: done
  real(WP) :: buf1,buf2
  
  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1 = 1.0_WP
  rho2 = 1.0_WP
  alpha = 0.0_WP
  omega = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res = 0.0_WP
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max(maxval(abs(res)),max_res0)
  if (max_res0.eq.0.0_WP) return

  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum(sum(res*res0),rho1)
     beta = alpha * rho1 / (rho2 * omega)
     pp = res + beta * (pp - omega * vv)
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum(sum(vv*res0),gamma)
     if (gamma.eq.0.0_WP) exit loop
     alpha = rho1/gamma
     DP = DP  + alpha*p_hat
     s  = res - alpha*vv
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum(sum(s*t),buf1)
     call parallel_sum(sum(t*t),buf2)
     if (buf2.eq.0.0_WP) exit loop
     omega = buf1/buf2
     DP  = DP + omega*s_hat
     res = s  - omega*t
     rho2 = rho1
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max(maxval(abs(res)),max_res)
     max_res = max_res/max_res0
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_3d


! ==================================== !
! Solver with 1 periodic direction : X !
! ==================================== !
subroutine bicgstab_solve_x
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: i

  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1_x = 1.0_WP
  rho2_x = 1.0_WP
  alpha_x = 0.0_WP
  omega_x = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res_x = 0.0_WP
  done_x = .false.
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max_dir(maxval(maxval(abs(res),dim=2),dim=2),max_res0_x,'yz')
  max_res0 = maxval(max_res0_x)
  if (max_res0.eq.0.0_WP) return
  
  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum_dir(sum(sum(res*res0,dim=2),dim=2),rho1_x,'yz')
     do i=imin_,imax_
        if (rho1_x(i).eq.0.0_WP) done_x(i) = .true.
        if (done_x(i)) then
           pp(i,:,:) = 0.0_WP
        else
           beta_x(i) = alpha_x(i) * rho1_x(i) / (rho2_x(i) * omega_x(i))
           pp(i,:,:) = res(i,:,:) + beta_x(i) * (pp(i,:,:) - omega_x(i) * vv(i,:,:))
        end if
     end do
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(sum(vv*res0,dim=2),dim=2),gamma_x,'yz')
     do i=imin_,imax_
        if (gamma_x(i).eq.0.0_WP) done_x(i) = .true.
        if (done_x(i)) then
           s(i,:,:) = 0.0_WP
        else
           alpha_x(i) = rho1_x(i)/gamma_x(i)
           DP(i,:,:) = DP(i,:,:)  + alpha_x(i)*p_hat(i,:,:)
           s(i,:,:) = res(i,:,:) - alpha_x(i)*vv(i,:,:)
        end if
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(sum(s*t,dim=2),dim=2),buf1_x,'yz')
     call parallel_sum_dir(sum(sum(t*t,dim=2),dim=2),buf2_x,'yz')
     do i=imin_,imax_
        if (buf2_x(i).eq.0.0_WP) done_x(i) = .true.
        if (done_x(i)) then
           res(i,:,:) = 0.0_WP
        else
           omega_x(i) = buf1_x(i)/buf2_x(i)
           DP(i,:,:)  = DP(i,:,:) + omega_x(i)*s_hat(i,:,:)
           res(i,:,:) = s(i,:,:)  - omega_x(i)*t(i,:,:)
        end if
        rho2_x(i) = rho1_x(i)
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(maxval(abs(res),dim=2),dim=2),max_res_x,'yz')
     do i=imin_,imax_
        if (done_x(i)) then
           max_res_x(i) = 0.0_WP
        else
           max_res_x(i) = max_res_x(i)/max_res0_x(i)
        end if
     end do
     max_res = maxval(max_res_x)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_x


! ==================================== !
! Solver with 1 periodic direction : Y !
! ==================================== !
subroutine bicgstab_solve_y
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: j

  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1_y = 1.0_WP
  rho2_y = 1.0_WP
  alpha_y = 0.0_WP
  omega_y = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res_y = 0.0_WP
  done_y = .false.
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max_dir(maxval(maxval(abs(res),dim=1),dim=2),max_res0_y,'xz')
  max_res0 = maxval(max_res0_y)
  if (max_res0.eq.0.0_WP) return
  
  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum_dir(sum(sum(res*res0,dim=1),dim=2),rho1_y,'xz')
     do j=jmin_,jmax_
        if (rho1_y(j).eq.0.0_WP) done_y(j) = .true.
        if (done_y(j)) then
           pp(:,j,:) = 0.0_WP
        else
           beta_y(j) = alpha_y(j) * rho1_y(j) / (rho2_y(j) * omega_y(j))
           pp(:,j,:) = res(:,j,:) + beta_y(j) * (pp(:,j,:) - omega_y(j) * vv(:,j,:))
        end if
     end do
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(sum(vv*res0,dim=1),dim=2),gamma_y,'xz')
     do j=jmin_,jmax_
        if (gamma_y(j).eq.0.0_WP) done_y(j) = .true.
        if (done_y(j)) then
           s(:,j,:) = 0.0_WP
        else
           alpha_y(j) = rho1_y(j)/gamma_y(j)
           DP(:,j,:) = DP(:,j,:)  + alpha_y(j)*p_hat(:,j,:)
           s(:,j,:) = res(:,j,:) - alpha_y(j)*vv(:,j,:)
        end if
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(sum(s*t,dim=1),dim=2),buf1_y,'xz')
     call parallel_sum_dir(sum(sum(t*t,dim=1),dim=2),buf2_y,'xz')
     do j=jmin_,jmax_
        if (buf2_y(j).eq.0.0_WP) done_y(j) = .true.
        if (done_y(j)) then
           res(:,j,:) = 0.0_WP
        else
           omega_y(j) = buf1_y(j)/buf2_y(j)
           DP(:,j,:)  = DP(:,j,:) + omega_y(j)*s_hat(:,j,:)
           res(:,j,:) = s(:,j,:)  - omega_y(j)*t(:,j,:)
        end if
        rho2_y(j) = rho1_y(j)
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(maxval(abs(res),dim=1),dim=2),max_res_y,'xz')
     do j=jmin_,jmax_
        if (done_y(j)) then
           max_res_y(j) = 0.0_WP
        else
           max_res_y(j) = max_res_y(j)/max_res0_y(j)
        end if
     end do
     max_res = maxval(max_res_y)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_y


! ==================================== !
! Solver with 1 periodic direction : Z !
! ==================================== !
subroutine bicgstab_solve_z
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: k

  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1_z = 1.0_WP
  rho2_z = 1.0_WP
  alpha_z = 0.0_WP
  omega_z = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res_z = 0.0_WP
  done_z = .false.
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max_dir(maxval(maxval(abs(res),dim=1),dim=1),max_res0_z,'xy')
  max_res0 = maxval(max_res0_z)
  if (max_res0.eq.0.0_WP) return
  
  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum_dir(sum(sum(res*res0,dim=1),dim=1),rho1_z,'xy')
     do k=kmin_,kmax_
        if (rho1_z(k).eq.0.0_WP) done_z(k) = .true.
        if (done_z(k)) then
           pp(:,:,k) = 0.0_WP
        else
           beta_z(k) = alpha_z(k) * rho1_z(k) / (rho2_z(k) * omega_z(k))
           pp(:,:,k) = res(:,:,k) + beta_z(k) * (pp(:,:,k) - omega_z(k) * vv(:,:,k))
        end if
     end do
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(sum(vv*res0,dim=1),dim=1),gamma_z,'xy')
     do k=kmin_,kmax_
        if (gamma_z(k).eq.0.0_WP) done_z(k) = .true.
        if (done_z(k)) then
           s(:,:,k) = 0.0_WP
        else
           alpha_z(k) = rho1_z(k)/gamma_z(k)
           DP(:,:,k) = DP(:,:,k)  + alpha_z(k)*p_hat(:,:,k)
           s(:,:,k)  = res(:,:,k) - alpha_z(k)*vv(:,:,k)
        end if
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(sum(s*t,dim=1),dim=1),buf1_z,'xy')
     call parallel_sum_dir(sum(sum(t*t,dim=1),dim=1),buf2_z,'xy')
     do k=kmin_,kmax_
        if (buf2_z(k).eq.0.0_WP) done_z(k) = .true.
        if (done_z(k)) then
           res(:,:,k) = 0.0_WP
        else
           omega_z(k) = buf1_z(k)/buf2_z(k)
           DP(:,:,k)  = DP(:,:,k) + omega_z(k)*s_hat(:,:,k)
           res(:,:,k) = s(:,:,k)  - omega_z(k)*t(:,:,k)
        end if
        rho2_z(k) = rho1_z(k)
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(maxval(abs(res),dim=1),dim=1),max_res_z,'xy')
     do k=kmin_,kmax_
        if (done_z(k)) then
           max_res_z(k) = 0.0_WP
        else
           max_res_z(k) = max_res_z(k)/max_res0_z(k)
        end if
     end do
     max_res = maxval(max_res_z)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_z


! ====================================== !
! Solver with 2 periodic directions : YZ !
! ====================================== !
subroutine bicgstab_solve_yz
  use bicgstab
  use parallel
  use metric_velocity_conv
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: j,k
    
  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1_yz = 1.0_WP
  rho2_yz = 1.0_WP
  alpha_yz = 0.0_WP
  omega_yz = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res_yz = 0.0_WP
  done_yz = .false.
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max_dir(maxval(abs(res),dim=1),max_res0_yz,'x')
  max_res0 = maxval(max_res0_yz)
  if (max_res0.eq.0.0_WP) return
  
  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum_dir(sum(res*res0,dim=1),rho1_yz,'x')
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (rho1_yz(j,k).eq.0_WP) done_yz(j,k) = .true.
           if (done_yz(j,k)) then
              pp(:,j,k) = 0.0_WP
           else
              beta_yz(j,k) = alpha_yz(j,k) * rho1_yz(j,k) / (rho2_yz(j,k) * omega_yz(j,k))
              pp(:,j,k) = res(:,j,k) + beta_yz(j,k) * (pp(:,j,k) - omega_yz(j,k) * vv(:,j,k))
           end if
        end do
     end do
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(vv*res0,dim=1),gamma_yz,'x')
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (gamma_yz(j,k).eq.0.0_WP) done_yz(j,k) = .true.
           if (done_yz(j,k)) then
              s(:,j,k) = 0.0_WP
           else
              alpha_yz(j,k) = rho1_yz(j,k)/gamma_yz(j,k)
              DP(:,j,k) = DP(:,j,k)  + alpha_yz(j,k)*p_hat(:,j,k)
              s(:,j,k)  = res(:,j,k) - alpha_yz(j,k)*vv(:,j,k)
           end if
        end do
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(s*t,dim=1),buf1_yz,'x')
     call parallel_sum_dir(sum(t*t,dim=1),buf2_yz,'x')
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (buf2_yz(j,k).eq.0.0_WP) done_yz(j,k) = .true.
           if (done_yz(j,k)) then
              res(:,j,k) = 0.0_WP
           else
              omega_yz(j,k) = buf1_yz(j,k)/buf2_yz(j,k)
              DP(:,j,k)  = DP(:,j,k) + omega_yz(j,k)*s_hat(:,j,k)
              res(:,j,k) = s(:,j,k)  - omega_yz(j,k)*t(:,j,k)
           end if
           rho2_yz(j,k) = rho1_yz(j,k)
        end do
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(abs(res),dim=1),max_res_yz,'x')
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (done_yz(j,k)) then
              max_res_yz(j,k) = 0.0_WP
           else
              max_res_yz(j,k) = max_res_yz(j,k)/max_res0_yz(j,k)
           end if
        end do
     end do
     max_res = maxval(max_res_yz)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_yz


! ====================================== !
! Solver with 2 periodic directions : XZ !
! ====================================== !
subroutine bicgstab_solve_xz
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: i,j,k
  
  !$OMP PARALLEL

  ! Initialization
  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           p_hat(i,j,k) = 0.0_WP
           s_hat(i,j,k) = 0.0_WP
          end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           vv(i,j,k) = 0.0_WP
           pp(i,j,k) = 0.0_WP
        end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmin_,kmax_
     do i=imin_,imax_
        rho1_xz(i,k) = 1.0_WP
        rho2_xz(i,k) = 1.0_WP
        alpha_xz(i,k) = 0.0_WP
        omega_xz(i,k) = 1.0_WP
        max_res_xz(i,k) = 0.0_WP
        done_xz(i,k) = .false.
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL

  iter = 0
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           res (i,j,k) = RP(i,j,k) - res(i,j,k)
           res0(i,j,k) = res(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  call parallel_max_dir(maxval(abs(res),dim=2),max_res0_xz,'y')
  max_res0 = maxval(max_res0_xz)
  if (max_res0.eq.0.0_WP) return

  ! Here we go
  loop:do while(done.ne.1)

     ! Increment
     iter = iter + 1

     ! Step 1
     call parallel_sum_dir(sum(res*res0,dim=2),rho1_xz,'y')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           if (rho1_xz(i,k).eq.0.0_WP) done_xz(i,k) = .true.
           if (done_xz(i,k)) then
              pp(i,:,k) = 0.0_WP
           else
              beta_xz(i,k) = alpha_xz(i,k) * rho1_xz(i,k) / (rho2_xz(i,k) * omega_xz(i,k))
              pp(i,:,k) = res(i,:,k) + beta_xz(i,k) * (pp(i,:,k) - omega_xz(i,k) * vv(i,:,k))
           end if
        end do
     end do
     !$OMP END PARALLEL DO
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(vv*res0,dim=2),gamma_xz,'y')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           if (gamma_xz(i,k).eq.0.0_WP) done_xz(i,k) = .true.
           if (done_xz(i,k)) then
              s(i,:,k) = 0.0_WP
           else
              alpha_xz(i,k) = rho1_xz(i,k)/gamma_xz(i,k)
              DP(i,:,k) = DP(i,:,k)  + alpha_xz(i,k)*p_hat(i,:,k)
              s(i,:,k)  = res(i,:,k) - alpha_xz(i,k)*vv(i,:,k)
           end if
        end do
     end do
     !$OMP END PARALLEL DO
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(s*t,dim=2),buf1_xz,'y')
     call parallel_sum_dir(sum(t*t,dim=2),buf2_xz,'y')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           if (buf2_xz(i,k).eq.0.0_WP) done_xz(i,k) = .true.
           if (done_xz(i,k)) then
              res(i,:,k) = 0.0_WP
           else
              omega_xz(i,k) = buf1_xz(i,k)/buf2_xz(i,k)
              DP(i,:,k)  = DP(i,:,k) + omega_xz(i,k)*s_hat(i,:,k)
              res(i,:,k) = s(i,:,k)  - omega_xz(i,k)*t(i,:,k)
           end if
           rho2_xz(i,k) = rho1_xz(i,k)
        end do
     end do
     !$OMP END PARALLEL DO
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(abs(res),dim=2),max_res_xz,'y')
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do i=imin_,imax_
           if (done_xz(i,k)) then
              max_res_xz(i,k) = 0.0_WP
           else
              max_res_xz(i,k) = max_res_xz(i,k)/max_res0_xz(i,k)
           end if
        end do
     end do
     !$OMP END PARALLEL DO

     max_res = maxval(max_res_xz)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter

  return
end subroutine bicgstab_solve_xz


! ====================================== !
! Solver with 2 periodic directions : XY !
! ====================================== !
subroutine bicgstab_solve_xy
  use bicgstab
  use parallel
  implicit none
  
  ! Private stuff
  integer :: done
  integer :: i,j
  
  ! Initialization
  vv = 0.0_WP
  pp = 0.0_WP
  rho1_xy = 1.0_WP
  rho2_xy = 1.0_WP
  alpha_xy = 0.0_WP
  omega_xy = 1.0_WP
  p_hat = 0.0_WP
  s_hat = 0.0_WP
  
  iter = 0
  max_res_xy = 0.0_WP
  done_xy = .false.
  done = 0
  
  ! Compute the initial residual
  call pressure_operator(DP,res)
  res  = RP - res
  res0 = res
  call parallel_max_dir(maxval(abs(res),dim=3),max_res0_xy,'z')
  max_res0 = maxval(max_res0_xy)
  if (max_res0.eq.0.0_WP) return
  
  ! Here we go
  loop:do while(done.ne.1)
     
     ! Increment
     iter = iter + 1
     
     ! Step 1
     call parallel_sum_dir(sum(res*res0,dim=3),rho1_xy,'z')
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (rho1_xy(i,j).eq.0.0_WP) done_xy(i,j) = .true.
           if (done_xy(i,j)) then
              pp(i,j,:) = 0.0_WP
           else
              beta_xy(i,j) = alpha_xy(i,j) * rho1_xy(i,j) / (rho2_xy(i,j) * omega_xy(i,j))
              pp(i,j,:) = res(i,j,:) + beta_xy(i,j) * (pp(i,j,:) - omega_xy(i,j) * vv(i,j,:))
           end if
        end do
     end do
     
     ! Preconditioning for p
     call bicgstab_apply_precond(pp,p_hat)
     
     ! vv = [A] p_hat
     call pressure_operator(p_hat,vv)
     
     ! Step 2
     call parallel_sum_dir(sum(vv*res0,dim=3),gamma_xy,'z')
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (gamma_xy(i,j).eq.0.0_WP) done_xy(i,j) = .true.
           if (done_xy(i,j)) then
              s(i,j,:) = 0.0_WP
           else
              alpha_xy(i,j) = rho1_xy(i,j)/gamma_xy(i,j)
              DP(i,j,:) = DP(i,j,:)  + alpha_xy(i,j)*p_hat(i,j,:)
              s(i,j,:)  = res(i,j,:) - alpha_xy(i,j)*vv(i,j,:)
           end if
        end do
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,s)
     !s = RP - s
     
     ! Preconditioning for s
     call bicgstab_apply_precond(s,s_hat)
     
     ! t = [A] s_hat
     call pressure_operator(s_hat,t)
     
     ! Step 3 & update
     call parallel_sum_dir(sum(s*t,dim=3),buf1_xy,'z')
     call parallel_sum_dir(sum(t*t,dim=3),buf2_xy,'z')
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (buf2_xy(i,j).eq.0.0_WP) done_xy(i,j) = .true.
           if (done_xy(i,j)) then
              res(i,j,:) = 0.0_WP
           else
              omega_xy(i,j) = buf1_xy(i,j)/buf2_xy(i,j)
              DP(i,j,:)  = DP(i,j,:) + omega_xy(i,j)*s_hat(i,j,:)
              res(i,j,:) = s(i,j,:)  - omega_xy(i,j)*t(i,j,:)
           end if
           rho2_xy(i,j) = rho1_xy(i,j)
        end do
     end do
     
     ! TEST - For now - Better convergence
     !call pressure_operator(DP,res)
     !res = RP - res
     
     ! Check if done
     call parallel_max_dir(maxval(abs(res),dim=3),max_res_xy,'z')
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (done_xy(i,j)) then
              max_res_xy(i,j) = 0.0_WP
           else
              max_res_xy(i,j) = max_res_xy(i,j)/max_res0_xy(i,j)
           end if
        end do
     end do
     max_res = maxval(max_res_xy)
     if (max_res.lt.cvg .or. iter.ge.max_iter) done=1
  end do loop
  
  ! Transfer to monitor
  max_resP = max_res
  it_p  = iter
  
  return
end subroutine bicgstab_solve_xy
