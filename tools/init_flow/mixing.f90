module mixing
  use precision
  use param
  implicit none

  ! Lengths of the domain
  real(WP) :: Lx,Ly,Lz

  ! Inflow profile
  character(len=str_medium) :: inflow_type
  ! Convective velocity and velocity difference
  real(WP) :: delta_U
  real(WP) :: U1,U2
  ! Half thickness of the layer
  real(WP) :: thick
  ! Reynolds number of the mean flow
  real(WP) :: Re
    
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX
  
  ! Flow Field stats
  real(WP), dimension(:), pointer :: Ubase   ! [m/s]
  real(WP), dimension(:), pointer :: Up      ! [1/s]
  real(WP), dimension(:), pointer :: Upp     ! [1/m*s]
  
  ! Temporal Orr Sommerfeld
  integer :: nos
  real(WP), dimension(:), pointer :: xos,yos
  real(WP), parameter :: stretching = 3.0_WP
  real(WP), parameter ::  h = 1.0E-6_WP
  real(WP),    dimension(:,:), pointer :: Id,D1,D2,D4,T,G,Tp
  complex(WP), dimension(:,:), pointer :: A,B,C
  complex(WP), dimension(:),   pointer :: eta_tilde,u_tilde,v_tilde,w_tilde
  
  ! Orr-Sommerfeld - lapack solver
  integer :: ierr,lwork
  integer,     dimension(:),   pointer :: iwork
  complex(WP), dimension(:),   pointer :: work
  real(WP),    dimension(:),   pointer :: rwork
  complex(WP), dimension(:),   pointer :: eigval_a,eigval_b  ! a/b = lambda = eigenvalue
  complex(WP), dimension(:,:), pointer :: eigvec_l           ! left eigenvector
  complex(WP), dimension(:,:), pointer :: eigvec_r           ! right eigenvector
  
  ! Complex constants
  complex(WP), parameter :: zero = (0.0_WP,0.0_WP)
  complex(WP), parameter :: one  = (1.0_WP,0.0_WP)
  complex(WP), parameter :: ii   = (0.0_WP,1.0_WP)

  ! Statistics
  real(WP) :: disp_thick,mom_thick,vort_thick,lambda
  real(WP) :: Re_disp,Re_mom,Re_vort,Re_lambda

end module mixing

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine mixing_grid
  use mixing
  use parser
  implicit none
  
  integer :: i,j,k
  
  ! Read in the size of the domain
  call parser_read('Lx',Lx)
  call parser_read('Ly',Ly)
  call parser_read('Lz',Lz)
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  ! Set the periodicity
  xper = 1
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  do i=1,nx+1
     x(i) = (i-1)*Lx/real(nx,WP)
  end do
  do j=1,ny+1
     y(j) = (j-1)*Ly/real(ny,WP) - 0.5_WP*Ly
  end do
  do k=1,nz+1
     z(k) = (k-1)*Lz/real(nz,WP) - 0.5_WP*Lz
  end do

  ! Create the mid points 
  do i=1,nx
     xm(i)= 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)= 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k)= 0.5_WP*(z(k)+z(k+1))
  end do
  
  ! Create the masks
  mask = 0  

  return
end subroutine mixing_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine mixing_data
  use mixing
  use fileio
  use parser
  use math
  implicit none

  integer :: i,j,k,iunit,n
  real(WP) :: alpha_opt,ratio,val,err,grad,rnd,amp_2d,amp_3d
  complex(WP) :: omega
  real(WP), dimension(9) :: beta,alpha,amp
  complex(WP), dimension(:), pointer :: u_fluct,v_fluct,w_fluct

  ! Allocate the array data
  nvar = 5
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Link the pointers
  U => data(:,:,:,1); names(1) = 'U'
  V => data(:,:,:,2); names(2) = 'V'
  W => data(:,:,:,3); names(3) = 'W'
  P => data(:,:,:,4); names(4) = 'P'
  ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  
  ! Read some parameters: mean flow
  call parser_read('Velocity difference',delta_U)
  call parser_read('Half layer thickness',thick)
  U1 = - 0.5_WP*delta_U
  U2 = + 0.5_WP*delta_U
  
  ! *** Mean Flow ***
  V = 0.0_WP
  W = 0.0_WP
  P = 0.0_WP
  ZMIX = 0.0_WP
  do j=1,ny
     U(:,j,:) = 0.5_WP*delta_U * tanh(ym(j)/thick)
     if (ym(j).ge.0.0_WP) ZMIX(:,j,:) = 1.0_WP
  end do

  ! *** Get most unstable mode ***

  ! Initialize the Orr-Sommerfeld solver
  call orr_sommerfeld_temporal_init

  ! Allocate perturbation arrays
  allocate(u_fluct(ny))
  allocate(v_fluct(ny))
  allocate(w_fluct(ny))

  ! Initial guess for omega (Valid for high Reynolds number)
  alpha_opt = 0.40_WP/thick

  ! Newton optimization
  iunit = iopen()
  open(iunit,file='alpha_opt.txt')
  write(*,*)
  write(*,'(a)') "# Finding the most unstable frequency"
  write(iunit,'(a)') "# Finding the most unstable frequency"

  ! Initial point
  call orr_sommerfeld_temporal_get_mode(alpha_opt,0.0_WP,omega,u_fluct,v_fluct,w_fluct)
  val = aimag(omega)
  
  !stop
  
  ! Direction of search
  call orr_sommerfeld_temporal_get_mode(alpha_opt+h,0.0_WP,omega,u_fluct,v_fluct,w_fluct)
  grad = (aimag(omega) - val)/h
  if (grad.gt.0.0_WP) then
     ratio = 2.0_WP
  else
     ratio = 0.5_WP
  end if
  err = 1.0_WP
  write(*,    '(a12,es16.6,a16,es12.4)') "alpha_opt =", alpha_opt, "error =", err
  write(iunit,'(a12,es16.6,a16,es12.4)') "alpha_opt =", alpha_opt, "error =", err
  
  ! Linear search
  do while(err.gt.h)
     call orr_sommerfeld_temporal_get_mode(alpha_opt*ratio+h,0.0_WP,omega,u_fluct,v_fluct,w_fluct)
     grad = aimag(omega)
     call orr_sommerfeld_temporal_get_mode(alpha_opt*ratio,  0.0_WP,omega,u_fluct,v_fluct,w_fluct)
     grad = (grad - aimag(omega))/h

     if (aimag(omega).gt.val) then
        alpha_opt = ratio*alpha_opt
        val = aimag(omega)
        if ((grad.lt.0.0_WP) .and. (ratio.gt.1.0_WP)) ratio = ratio**(-0.1_WP)
        if ((grad.gt.0.0_WP) .and. (ratio.lt.1.0_WP)) ratio = ratio**(-0.1_WP)
     else
        ratio = ratio**(0.1_WP)
     end if

     err = max(ratio-1.0_WP,1.0_WP/ratio-1.0_WP)
     !write(*,'(10es12.4)'),alpha_opt,ratio,grad,aimag(omega),val
     write(*,    '(a12,es16.6,a16,es12.4)') "alpha_opt =", alpha_opt, "error =", err
     write(iunit,'(a12,es16.6,a16,es12.4)') "alpha_opt =", alpha_opt, "error =", err
  end do

  close(iclose(iunit))
  write(*,*)

  write(*,'(a24,es16.6)') "Optimal domain length Lx=", 4.0_WP * twoPi / alpha_opt
  write(*,*)

  !stop

  ! **** Create the flow field ****
  ! 3D - Fluctuations
  beta(1) =  0.0_WP
  beta(2) =  0.0_WP
  beta(3) =  0.0_WP
  beta(4) =  +alpha_opt
  beta(5) =  +alpha_opt/2.0_WP
  beta(6) =  +alpha_opt/4.0_WP
  beta(7) =  -alpha_opt
  beta(8) =  -alpha_opt/2.0_WP
  beta(9) =  -alpha_opt/4.0_WP

  alpha(1) = alpha_opt
  alpha(2) = alpha_opt/2.0_WP
  alpha(3) = alpha_opt/4.0_WP
  alpha(4) = alpha_opt
  alpha(5) = alpha_opt/2.0_WP
  alpha(6) = alpha_opt/4.0_WP
  alpha(7) = alpha_opt
  alpha(8) = alpha_opt/2.0_WP
  alpha(9) = alpha_opt/4.0_WP
  
  call parser_read('Fluct. amplitude 2D',amp_2d)
  call parser_read('Fluct. amplitude 3D',amp_3d)
  amp(1:3) = amp_2d
  amp(4:6) = amp_3d
  write(*,'(a)') "# Setting the modes"
  write(*,'(a6,4a16)') "mode","alpha","beta","omega_r","omega_i"

  ! Get purturbations
  call random_init
  do n=1,9
     call orr_sommerfeld_temporal_get_mode(alpha(n),beta(n),omega,u_fluct,v_fluct,w_fluct)
     write(*,'(i6,4es16.6)') n,alpha(n),beta(n),omega
     if (n.gt.3) then
        call random_number(rnd)
        rnd = rnd*twoPi
     else
        rnd = 0.0_WP
     end if
     do k=1,nz
        do j=1,ny
           do i=1,nx
              U(i,j,k) = U(i,j,k) + amp(n)*real (u_fluct(j)*exp(ii*(rnd+beta(n)*zm(k)-alpha(n)*x (i))))
              V(i,j,k) = V(i,j,k) + amp(n)*real (v_fluct(j)*exp(ii*(rnd+beta(n)*zm(k)-alpha(n)*xm(i))))
              W(i,j,k) = W(i,j,k) + amp(n)*real (w_fluct(j)*exp(ii*(rnd+beta(n)*z(k) -alpha(n)*xm(i))))
           end do
        end do
     end do
  end do
  
  ! Print some statistic on the flow field
  call mixing_stats

  return
end subroutine mixing_data


! ====================================== !
! Initialize the Orr_Sommerfeld solver   !
! Allocate and pre-compute some matrices !
! ====================================== !
subroutine orr_sommerfeld_temporal_init
  use mixing
  use parser
  use fileio
  use math
  implicit none

  integer :: i,j,iunit

  ! Read number of grid points
  call parser_read('Orr-Sommerfeld points',nos)

  ! Get the operating Reynolds number
  call parser_read('Reynolds number',Re)

  ! Create grid
  allocate(xos(nos))
  allocate(yos(nos))
  xos(1)   = +1.0_WP
  yos(1)   = +huge(1.0_WP)
  xos(nos) = -1.0_WP
  yos(nos) = -huge(1.0_WP)
  do j=2,nos-1
     xos(j) = cos(Pi*real(j-1,WP)/real(nos-1,WP))
     yos(j) = stretching * 0.5_WP * thick * log((1.0_WP+xos(j))/(1.0_WP-xos(j)))
  end do
  
  ! Orr-Sommerfeld - lapack solver
  allocate(work(2*nos))
  allocate(rwork(8*nos))
  allocate(eigval_a(nos))
  allocate(eigval_b(nos))
  allocate(eigvec_l(1,nos))
  allocate(eigvec_r(nos,nos))
  allocate(iwork(nos))

  ! Temporary solution from Orr-Sommerfeld
  allocate(u_tilde(nos))
  allocate(v_tilde(nos))
  allocate(w_tilde(nos))
  allocate(eta_tilde(nos))

  ! Setup eigenvalue problem
  allocate(A(nos,nos))
  allocate(B(nos,nos))
  allocate(C(nos,nos))
  
  ! Precompute the stencil for derivations
  allocate(Id(nos,nos))
  allocate(D1(nos,nos))
  allocate(D2(nos,nos))
  allocate(D4(nos,nos))
  
  ! Chebyshev matrices
  allocate(T (nos,nos))
  allocate(G (nos,nos))
  allocate(Tp(nos,nos))
  
  Id = 0.0_WP
  D1 = 0.0_WP
  D2 = 0.0_WP
  D4 = 0.0_WP

  do i=1,nos

     Id(i,i) = 1.0_WP

     do j=1,nos
        ! Chebyshev to real
        T(i,j) = cos(Pi*real((i-1)*(j-1),WP)/real(nos-1,WP))

        ! Real to Chebyshev
        Tp(i,j) = 2.0_WP * cos(Pi*real((i-1)*(j-1),WP)/real(nos-1,WP)) / real(nos-1,WP)
        if (i.eq.1 .or. i.eq.nos) Tp(i,j) = 0.5_WP * Tp(i,j)
        if (j.eq.1 .or. j.eq.nos) Tp(i,j) = 0.5_WP * Tp(i,j)
        
        ! Derivation in Chebyshev space
        if (i.ge.j .or. mod(i+j,2).eq.0) then
           G(i,j) = 0.0_WP
        else
           if (i.eq.1 .or. i.eq.nos) then
              G(i,j) = real(j-1,WP)
           else
              G(i,j) = 2.0_WP*real(j-1,WP)
           end if
        end if
        
     end do
  end do
  
  ! Derivation in real space
  call dgemm('N','N',nos,nos,nos,1.0_WP,T,nos,G,nos,0.0_WP,D1,nos)
  G = D1
  call dgemm('N','N',nos,nos,nos,1.0_WP,G,nos,Tp,nos,0.0_WP,D1,nos)
  
  ! Account for the mapping
  do i=1,nos
     D1(i,:) = D1(i,:) * (1.0_WP-xos(i)**2) / (thick*stretching)
  end do
     
  ! Higher order derivatives
  call dgemm('N','N',nos,nos,nos,1.0_WP,D1,nos,D1,nos,0.0_WP,D2,nos)
  call dgemm('N','N',nos,nos,nos,1.0_WP,D2,nos,D2,nos,0.0_WP,D4,nos)
  
  ! Create base flow
  allocate(Ubase(nos))
  allocate(Up   (nos))
  allocate(Upp  (nos))

  do j=1,nos
     Ubase(j) = 0.5_WP*delta_U * tanh(yos(j)/thick)
     Up(j)    = 0.5_WP*delta_U/thick / cosh(yos(j)/thick)**2
     Upp(j)   = -delta_U/thick**2 * tanh(yos(j)/thick) / cosh(yos(j)/thick)**2
  end do

  ! Save to a file the mean flow
  iunit = iopen()
  open(iunit,file='mean.txt')
  do j=1,nos
     write(iunit,'(i3,4e16.6E3)') j,yos(j),Ubase(j),Up(j),Upp(j)
  end do
  close(iclose(iunit))

  return
end subroutine orr_sommerfeld_temporal_init


! ======================================================= !
! Compute the solution of the Orr-Sommerfeld equation for !
! -> Re    : Reynolds number                              !
! -> alpha : wave number in x                             !
! -> beta  : wave number in z                             !
! ======================================================= !
subroutine orr_sommerfeld_temporal_get_mode(alpha,beta,omega,u_fluct,v_fluct,w_fluct)
  use mixing
  use fileio
  implicit none
  
  real(WP), intent(in) :: alpha,beta
  complex(WP), intent(out) :: omega
  complex(WP), dimension(ny) :: u_fluct,v_fluct,w_fluct

  complex(WP), dimension(:), pointer :: coeff_u,coeff_v,coeff_w
  real(WP) :: speed,speed_,Umin,Umax,norm,xx,num,den,theta,interp
  integer  :: j,jj,mode_index,iunit
  complex(WP) :: omega_

  A = zero
  B = zero

  ! Interior points
  do j=2,nos-1
     A(j,:) = &
          + alpha*Ubase(j)*(D2(j,:)-(alpha**2+beta**2)*Id(j,:)) &
          - alpha*Upp(j)*Id(j,:) &
          +ii*(D4(j,:)-2.0_WP*(alpha**2+beta**2)*D2(j,:)+(beta**2+alpha**2)**2*Id(j,:))/Re
     
     B(j,:) = D2(j,:)-(alpha**2+beta**2)*Id(j,:)
  end do
  
  ! BC - Zero val at inft
  A(1,1)     = one
  A(nos,nos) = one

  ! Get eigenvalues and eigenvectors
  lwork = 2*nos
  call ZGGEV('N','V',nos,A,nos,B,nos,eigval_a,eigval_b,eigvec_l,1,eigvec_r,nos,work,lwork,rwork,ierr)
  if (ierr.ne.0) stop "ZGGEV: unable to compute the eigenvalues/vectors"

  ! Get the most unstable mode
  ! -> phase speed (real part) between min and max of Ubase
  ! -> negative imaginary part
  Umin = minval(Ubase)
  Umax = maxval(Ubase)
  
  mode_index = -1
  omega = -ii*huge(1.0_WP)
  
  ! Detect the modes
  do j=1,nos
     omega_ = eigval_a(j)/(eigval_b(j)+epsilon(abs(eigval_b(j))))
     speed_ = real(omega_)/alpha
     
     !write(*,'(i3,8e16.6e3)') j,omega_,speed_,eigval_a(j),eigval_b(j)

     if (abs(eigval_b(j)).ne.0.0_WP) then
        if (speed_.ge.Umin .and. speed_.le.Umax) then
           if (aimag(omega_).gt.aimag(omega)) then
              mode_index = j
              omega = omega_
              speed = speed_
           end if
        end if
     end if
  end do
  !write(*,'(3i3)') mode_index
  if (mode_index.eq.-1) &
       stop "orr_sommerfeld_get_mode: Could not find a correct mode"

  ! Get v tilde
  v_tilde = eigvec_r(:,mode_index)

  ! Rotate for symmetry on real part
  num = 0.0_WP
  den = 0.0_WP
  do j=1,nos/2
     num = num + (real (v_tilde(j))-real (v_tilde(nos-j+1)))*(aimag(v_tilde(j))-aimag(v_tilde(nos-j+1)))
     den = den + (aimag(v_tilde(j))-aimag(v_tilde(nos-j+1)))*(aimag(v_tilde(j))-aimag(v_tilde(nos-j+1)))
  end do
  theta = atan(num/den)
  v_tilde = exp(ii*theta)*v_tilde

  ! Create the linear system for eta
  C = zero
  do j=2,nos-1
     C(j,:) = (alpha*Ubase(j)-omega)*Id(j,:) + ii*(D2(j,:)-(alpha**2+beta**2)*Id(j,:))/Re
     eta_tilde(j) = -ii*beta*Up(j)*v_tilde(j)
  end do
  
  ! Boundary conditons
  C(1,1)     = one
  C(nos,nos) = one
  eta_tilde(1)   = zero
  eta_tilde(nos) = zero
  
  ! Solve it (linear solver)
  call ZGESV(nos,1,C,nos,iwork,eta_tilde,nos,ierr)
  
  ! Get back the two remaining components
  do j=1,nos
     u_tilde(j) = ii*(alpha*sum(D1(j,:)*v_tilde)-beta*eta_tilde(j)) / (alpha**2+beta**2)
  end do
  w_tilde = (ii*eta_tilde+beta*u_tilde)/alpha
  
  ! High Freq filter
  call high_freq_filter(u_tilde)
  call high_freq_filter(v_tilde)
  call high_freq_filter(w_tilde)

  ! Write the solution to a file
  iunit = iopen()
  open(iunit,file='mode1.txt',form="formatted",iostat=ierr)
  do jj=1,nos
     write(iunit,'(i3,6e14.6)') jj,u_tilde(jj),v_tilde(jj),w_tilde(jj)
  end do
  close(iclose(iunit))     

  ! Compute the Chebyshev coefficients
  coeff_u => work(1:nos)
  coeff_v => work(nos+1:2*nos)
  coeff_w => work(2*nos+1:3*nos)
  do j=1,nos
     coeff_u(j) = sum(Tp(j,:)*u_tilde)
     coeff_v(j) = sum(Tp(j,:)*v_tilde)
     coeff_w(j) = sum(Tp(j,:)*w_tilde)
  end do
  
  ! Interpolate for the grid points
  do j=1,ny
     ! U, W
     xx = tanh(ym(j)/(thick*stretching))
     loop1:do jj=nos-1,1,-1
        if (xx.lt.xos(jj)) exit loop1
     end do loop1
     jj = max(jj,1)
     interp = (xx-xos(jj+1))/(xos(jj)-xos(jj+1))
     u_fluct(j) = interp*u_tilde(jj) + (1.0_WP-interp)*u_tilde(jj+1)
     w_fluct(j) = interp*w_tilde(jj) + (1.0_WP-interp)*w_tilde(jj+1)

     ! V
     xx = tanh(y(j)/(thick*stretching))
     loop2:do jj=nos-1,1,-1
        if (xx.lt.xos(jj)) exit loop2
     end do loop2
     jj = max(jj,1)
     interp = (xx-xos(jj+1))/(xos(jj)-xos(jj+1))
     v_fluct(j) = interp*v_tilde(jj) + (1.0_WP-interp)*v_tilde(jj+1)
  end do

  ! Normalize everything
  norm = maxval(abs(u_fluct)**2+abs(v_fluct)**2+abs(w_fluct)**2)
  u_fluct = u_fluct / norm
  v_fluct = v_fluct / norm
  w_fluct = w_fluct / norm

  ! Write the solution to a file
  iunit = iopen()
  open(iunit,file='mode.txt',form="formatted",iostat=ierr)
  do jj=1,ny
     write(iunit,'(i3,8e14.6)') jj,y(jj),ym(jj),u_fluct(jj),v_fluct(jj),w_fluct(jj)
  end do
  close(iclose(iunit))     

  return
end subroutine orr_sommerfeld_temporal_get_mode


! ======================================= !
! Remove high frequencies from the signal !
! Remove discontinuities at the borders   !
! ======================================= !
subroutine high_freq_filter(uin)
  use mixing
  implicit none
  
  complex(WP), dimension(nos) :: uin
  real(WP), parameter :: eps = 0.9_WP
  complex(WP) :: err
  integer :: j

  eta_tilde(1)   = zero
  eta_tilde(nos) = zero
  do j=2,nos-1
     eta_tilde(j) = 0.5_WP*(uin(j) + 0.5_WP*(uin(j-1)+uin(j+1)))
  end do
  uin = eta_tilde
  
  err = 0.75_WP*uin(2)-0.25_WP*uin(3)
  do j=2,nos/2
     uin(j) = uin(j) - err
     err = eps*err
  end do

  err = 0.75_WP*uin(nos-1)-0.25_WP*uin(nos-2)
  do j=nos-1,nos/2,-1
     uin(j) = uin(j) - err
     err = eps*err
  end do

  return
end subroutine high_freq_filter


! =========================================== !
! Compute some statistics of the mixing layer !
! -> recomputes the Ubase !!                  !
! =========================================== !
subroutine mixing_stats
  use mixing
  implicit none

  real(WP) :: nu,up2,dupdx2,tmp,uprime1,uprime2,max_der
  integer :: i,j,k
  
  ! Compute mixing layer stats
  deallocate(Ubase)
  allocate(Ubase(ny))
  Ubase = 0.0_WP
  do j=1,ny
     Ubase(j) = sum(U(:,j,:))
  end do
  Ubase = Ubase/real(nx*nz,WP)
  
  ! Upper-lower velocities
  U1 = Ubase(ny)
  U2 = Ubase(1)
  
  ! Viscosity taken to be inverse of Reynolds
  nu = 0.5_WP*(U1-U2)*thick/Re
  
  ! Calculate the displacement/momentum thickness
  disp_thick = 0.0_WP
  mom_thick  = 0.0_WP
  do j=1,ny
     disp_thick = disp_thick + (0.5_WP*(U1-U2)-abs(Ubase(j)-0.5_WP*(U1+U2))) * (y(j+1)-y(j))
     mom_thick  = mom_thick  + (U1-Ubase(j))*(Ubase(j)-U2) * (y(j+1)-y(j))
  end do
  mom_thick  = mom_thick/(U1-U2)**2
  disp_thick = 2.0_WP*disp_thick/(U1-U2)
  
  ! Calculate the vorticity thickness
  max_der = 0.0_WP
  do j=1,ny-1
     tmp = abs(Ubase(j+1)-Ubase(j))/(ym(j+1)-ym(j))
     if (tmp.gt.max_der) max_der = tmp
  end do
  vort_thick = (U1-U2)/max_der
  
  ! Calculate the taylor microscale (lambda) and Re_lambda
  up2 = 0.0_WP
  dupdx2 = 0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx-1
           uprime1 = U(i,j,k)-Ubase(j)
           uprime2 = U(i+1,j,k)-Ubase(j)
           up2 = up2 + ((uprime1+uprime2)/2.0_WP)**2
           dupdx2 = dupdx2 + ((uprime2-uprime1)/(x(i+1)-x(i)))**2
        end do
     end do
  end do
  lambda = sqrt( up2/dupdx2 )

  ! Compute the Reynolds numbers
  re_disp   = disp_thick*(U1-U2)/nu 
  re_mom    = mom_thick *(U1-U2)/nu 
  re_vort   = vort_thick*(U1-U2)/nu 
  re_lambda = lambda    *(U1-U2)/nu

  ! Print everything
  write(*,*)
  write(*,'(a14,es16.6,a4,es16.6)') "Disp. thick. ", disp_thick,"  Re ", re_disp
  write(*,'(a14,es16.6,a4,es16.6)') "Mom. thick. ", mom_thick,"  Re ", re_mom
  write(*,'(a14,es16.6,a4,es16.6)') "Vort. thick. ", vort_thick,"  Re ", re_vort
  write(*,'(a14,es16.6,a4,es16.6)') "Lambda ", lambda,"  Re ", re_lambda
  write(*,*)

  return
end subroutine mixing_stats

