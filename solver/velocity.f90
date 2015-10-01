module velocity
  use precision
  use partition
  use geometry
  use time_info
  implicit none
  
  ! Solution at old time step : n
  real(WP), dimension(:,:,:), pointer :: Uold
  real(WP), dimension(:,:,:), pointer :: Vold
  real(WP), dimension(:,:,:), pointer :: Wold
  real(WP), dimension(:,:,:), pointer :: rhoUold
  real(WP), dimension(:,:,:), pointer :: rhoVold
  real(WP), dimension(:,:,:), pointer :: rhoWold
  
  ! Solution at half time step : n+1
  real(WP), dimension(:,:,:), pointer :: RHOmid

  ! Source term between time n and n+1
  real(WP), dimension(:,:,:), pointer :: srcU
  real(WP), dimension(:,:,:), pointer :: srcV
  real(WP), dimension(:,:,:), pointer :: srcW
  real(WP), dimension(:,:,:), pointer :: srcP
  
  ! Residuals between time n and n+1
  real(WP), dimension(:,:,:), pointer :: ResU
  real(WP), dimension(:,:,:), pointer :: ResV
  real(WP), dimension(:,:,:), pointer :: ResW
  
  ! Divergence
  real(WP), dimension(:,:,:), pointer :: divg
  
  ! Values to monitor
  real(WP) :: max_u,max_v,max_w,max_p,max_divg
  real(WP) :: max_resU,max_resV,max_resW

  ! Using velocity relaxation
  logical :: velocity_relax
  real(WP) :: relax_alpha

end module velocity


! ============================================================== !
! Initilize the velocity module                                  !
!                                                                ! 
! -> allocate the arrays                                         !
! -> update ghost cells                                          !
! -> apply boundary conditions                                   !
!                                                                !
! Before: U/V/W OR rhoU/rhoV/rhoW correct only inside the domain !
!         -> imin_:imax_,jmin_:jmax_,kmin_:kmax_                 !
! After : U/V/W/ AND rhoU/rhoV/rhoW correct everywhere           !
!         -> imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_           !
! ============================================================== !
subroutine velocity_init
  use velocity
  use data
  use parser
  implicit none

  integer :: i,j,k
  
  ! Create & Start the timer
  call timing_create('velocity')
  call timing_start ('velocity')

  ! Velocity relaxation
  call parser_read('Velocity relaxation',velocity_relax,.false.)
  call parser_read('Relaxation alpha',relax_alpha,0.5_WP)

  ! Allocate arrays for old solution
  allocate(Uold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Vold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Wold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoUold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoVold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoWold (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Allocate arrays for new time solution
  allocate(RHOmid (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Allocate arrays for source terms
  allocate(srcU (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(srcV (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(srcW (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(srcP (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! Zero the source term
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           srcU(i,j,k) = 0.0_WP
           srcV(i,j,k) = 0.0_WP
           srcW(i,j,k) = 0.0_WP
           srcP(i,j,k) = 0.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Allocate Residuals for momentum equation 
  allocate(ResU (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(ResV (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(ResW (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  
  ! Allocate array for divergence
  allocate(divg (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  
  ! Compute density at n+1 for rho_multiply/divide
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHOmid(i,j,k) = 0.5_WP*(RHO(i,j,k)+RHOold(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! The velocity/momentum fields read from the file were without ghost cells
  ! Update the ghost cells for periodicity and domain decomposition
  ! Take care of the pressure here rather than in pressure.f90
  if (vel_present) then
     call boundary_update_border(U,'+','ym')
     call boundary_update_border(V,'-','y')
     call boundary_update_border(W,'-','ym')
  end if
  if (mom_present) then
     call boundary_update_border(rhoU,'+','ym')
     call boundary_update_border(rhoV,'-','y')
     call boundary_update_border(rhoW,'-','ym')
  end if
  call boundary_update_border(P,'+','ym')
  
  ! Synchronise U/V/W - rhoU/rhoV/rhoW
  if ((mom_present) .and. (.not.vel_present)) call rho_divide
  
  ! Apply boundary conditions on U/V/W
  call boundary_velocity_dirichlet
  call boundary_velocity_outflow
  call boundary_velocity_neumann
  
  ! Synchronise U/V/W - rhoU/rhoV/rhoW
  call rho_multiply
  
  ! Apply boundary conditions on rhoU/rhoV/rhoW
  call boundary_momentum_neumann
  
  ! Predict outlet velocity
  call velocity_predict_outlet
  
  ! Initialize the implicit solver
  call implicit_init
  
  ! Create a new file to monitor at each timestep
  call monitor_create_file_step('velocity',5)
  call monitor_set_header (1,'max_u','r')
  call monitor_set_header (2,'max_v','r')
  call monitor_set_header (3,'max_w','r')
  call monitor_set_header (4,'max_p','r')
  call monitor_set_header (5,'max_divg','r')
  
  ! Create a new file to monitor at each subiteration
  call monitor_create_file_iter('convergence_velocity',3)
  call monitor_set_header (1,'res_u','r')
  call monitor_set_header (2,'res_v','r')
  call monitor_set_header (3,'res_w','r')
  
  ! Stop the timer
  call timing_stop('velocity')
  
  return
end subroutine velocity_init


! ====================================================== !
! PRE-TIMESTEP Routine                                   !
!   -> store old solution                                !
!   -> compute the source term for the momentum equation !
! ====================================================== !
subroutine velocity_prestep
  use velocity
  use data
  implicit none

  integer :: i,j,k
  
  ! Start the timer
  call timing_start('velocity')

  ! Zero the source term
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           srcU(i,j,k) = 0.0_WP
           srcV(i,j,k) = 0.0_WP
           srcW(i,j,k) = 0.0_WP
           srcP(i,j,k) = 0.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Compute all the body forces
  call bodyforce_src
  
  ! Compute SGS source terms
  call sgsmodel_src_vel(srcU,srcV,srcW)
  
  ! Save the old velocity
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           Uold(i,j,k) = U(i,j,k)
           Vold(i,j,k) = V(i,j,k)
           Wold(i,j,k) = W(i,j,k)
           rhoUold(i,j,k) = rhoU(i,j,k)
           rhoVold(i,j,k) = rhoV(i,j,k)
           rhoWold(i,j,k) = rhoW(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Stop the timer
  call timing_stop('velocity')

  return
end subroutine velocity_prestep


! ======================================================= !
! ADVANCE the solution                                    !
!   -> second order in time and space                     !
!   -> explicit prediction                                !
!   -> ADI solver for each directions/components          !
!                                                         !
! U*(n+1,k+1) = U(n) + dt*F(0.5*(U(n+1,k)+U(n))           !
!                    + 0.5*dt*dF/dU*(U(n+1,k+1)-U(n+1,k)) !
! n : time step                                           !
! k : inner loop iteration                                !
! ======================================================= !
subroutine velocity_step
  use velocity
  use data
  use parallel
  implicit none
  integer :: i,j,k
  real(WP) :: alpha
  
  ! Start the timer
  call timing_start('velocity')

  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Compute density at n+1 for rho_multiply/divide
           RHOmid(i,j,k) = 0.5_WP*(RHO(i,j,k)+RHOold(i,j,k))
  
           ! Compute mid point
           ! Store it in the 'n+1' velocity
           rhoU(i,j,k) = 0.5_WP*(rhoU(i,j,k)+rhoUold(i,j,k))
           rhoV(i,j,k) = 0.5_WP*(rhoV(i,j,k)+rhoVold(i,j,k))
           rhoW(i,j,k) = 0.5_WP*(rhoW(i,j,k)+rhoWold(i,j,k))
           U(i,j,k) = 0.5_WP*(U(i,j,k)+Uold(i,j,k))
           V(i,j,k) = 0.5_WP*(V(i,j,k)+Vold(i,j,k))
           W(i,j,k) = 0.5_WP*(W(i,j,k)+Wold(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Solve for U
  call velocity_residuals_u
  call velocity_inverse_u
  
  ! Solve for V
  call velocity_residuals_v
  call velocity_inverse_v
  
  ! Solve for W
  call velocity_residuals_w
  call velocity_inverse_w
  
  ! Update U,V,W
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           U(i,j,k) = 2.0_WP*U(i,j,k)-Uold(i,j,k) + ResU(i,j,k)
           V(i,j,k) = 2.0_WP*V(i,j,k)-Vold(i,j,k) + ResV(i,j,k)
           W(i,j,k) = 2.0_WP*W(i,j,k)-Wold(i,j,k) + ResW(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Set the monitor values
  call parallel_max(maxval(abs(ResU)),max_resU)
  call parallel_max(maxval(abs(ResV)),max_resV)
  call parallel_max(maxval(abs(ResW)),max_resW)

  ! SOR
  if (velocity_relax) then
!  alpha = 0.5e-1_WP !0.5 by default
     alpha = relax_alpha
     U = alpha*U+(1.0_WP-alpha)*Uold
     V = alpha*V+(1.0_WP-alpha)*Vold
     W = alpha*W+(1.0_WP-alpha)*Wold
  end if

  ! Update the physical boundaries
  call boundary_velocity_dirichlet
  call boundary_velocity_neumann
  call boundary_velocity_outflow
  call boundary_update_border(U,'+','ym')
  call boundary_update_border(V,'-','y')
  call boundary_update_border(W,'-','ym')
  
  ! Synchronise U/rhoU
  call rho_multiply
  call boundary_momentum_neumann

  ! Transfer values to monitor
  call monitor_select_file('convergence_velocity')
  call monitor_set_single_value(1,max_resU)
  call monitor_set_single_value(2,max_resV)
  call monitor_set_single_value(3,max_resW)

  ! Stop the timer
  call timing_stop('velocity')
  
  return
end subroutine velocity_step


! ========================================= !
! Compute the density at the next time step !
! ========================================= !
subroutine density_step
  use velocity
  use data
  use metric_velocity_conv
  implicit none

  integer :: i,j,k

  ! Start the timer
  call timing_start('velocity')

  ! Solve for RHO
  call velocity_compute_divergence
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           RHO(i,j,k) = RHOold(i,j,k) - dt*( &
                  sum(divc_u(i,j,:)*rhoU(i-stc1:i+stc2,j,k)) &
                + sum(divc_v(i,j,:)*rhoV(i,j-stc1:j+stc2,k)) &
                + sum(divc_w(i,j,:)*rhoW(i,j,k-stc1:k+stc2)))                
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Update the boundaries
  call boundary_update_border(RHO,'+','ym')

  ! Stop the timer
  call timing_stop('velocity')

  return
end subroutine density_step


! =========================================== !
! Compute and apply the pressure correction   !
! U(n+1,k+1) = U*(n+1,k+1) - dt*grad(DP(n+1)) !
! =========================================== !
subroutine velocity_pressure
  use velocity
  use data
  use pressure
  implicit none
  real(WP) :: dti
  integer :: i,j,k
  
  ! Start the timer
  call timing_start('velocity')

  dti = 1.0_WP/dt_uvw

  ! Adjust mass fluxes to ensure the RHS of the Poisson equation
  ! to be of volumetric integral null
  ! Time step size for scalar to be used here (dRHO/dt)
  call boundary_velocity_massflux
  
  ! Compute divergence
  call velocity_compute_divergence
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           RP(i,j,k) = dti * divg(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Solve the Poisson equation for the pressure
  ! and update the pressure with the correction
  call pressure_step

  ! Apply pressure correction
  call velocity_apply_pressure

  ! Update borders
  call boundary_update_border(rhoU,'+','ym')
  call boundary_update_border(rhoV,'-','y')
  call boundary_update_border(rhoW,'-','ym')

  ! Synchronise U/rhoU
  call rho_divide

  ! Stop the timer
  call timing_stop('velocity')

  return
end subroutine velocity_pressure


! ======================================================== !
! From the predicted new velocity, predict the new density !
! ======================================================== !
subroutine velocity_predict_density
  use velocity
  use data
  implicit none

  integer :: i,j,k
  
  ! No density prediction if running at constant density
  if (constant_density) return
  
  ! Compute the divergence
  call velocity_compute_divergence
  
  ! Predict new density from predicted velocity
  !RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) =         &
  !     + RHOold(imin_:imax_,jmin_:jmax_,kmin_:kmax_) &
  !     + dRHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)   &
  !     - divg * dt
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHO(i,j,k) = RHO(i,j,k) + dRHO(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Update ghost cells
  call boundary_density_neumann
  call boundary_update_border(RHO,'+','ym')
  
  return
end subroutine velocity_predict_density


! ============================================ !
! Compute the error in the continuity equation !
! ============================================ !
subroutine velocity_maxdivergence
  use velocity
  use parallel
  implicit none
  
  ! Compute the divergence
  call velocity_compute_divergence
  
  ! Get the global max
  call parallel_max(abs(divg),max_divg)
  
  return
end subroutine velocity_maxdivergence


! ============================================ !
! Compute the error in the continuity equation !
! ============================================ !
subroutine velocity_compute_divergence
  use velocity
  use data
  use metric_velocity_conv
  implicit none
  integer :: i,j,k
  
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           divg(i,j,k) = &
                + dRHO(i,j,k)/dt &
                + sum(divc_u(i,j,:)*rhoU(i-stc1:i+stc2,j,k)) &
                + sum(divc_v(i,j,:)*rhoV(i,j-stc1:j+stc2,k)) &
                + sum(divc_w(i,j,:)*rhoW(i,j,k-stc1:k+stc2)) &
                - srcP(i,j,k)/dt
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine velocity_compute_divergence


! ========================================= !
! Apply pressure correction on the velocity !
! ========================================= !
subroutine velocity_apply_pressure
  use velocity
  use pressure
  use data
  use metric_velocity_conv
  implicit none
  
  integer :: i,j,k

  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoU(i,j,k) = rhoU(i,j,k) - dt_uvw*sum(grad_Px(i,j,:)*DP(i-stc2:i+stc1,j,k))
           rhoV(i,j,k) = rhoV(i,j,k) - dt_uvw*sum(grad_Py(i,j,:)*DP(i,j-stc2:j+stc1,k))
           rhoW(i,j,k) = rhoW(i,j,k) - dt_uvw*sum(grad_Pz(i,j,:)*DP(i,j,k-stc2:k+stc1))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
 
!!$  do k=kmin_,kmax_
!!$     do i=imin_,imax_
!!$        !rhoV(i,jmin,k) = rhoV(i,jmin,k) - dt_uvw*(DP(i,jmin,k)-DP(i,jmin-1,k))*dymi(jmin-1)
!!$        rhoV(i,jmin,k) = rhoV(i,jmin,k) - 0.25_WP*dt_uvw*(DP(i,jmin+1,k)+DP(i,jmin,k)-DP(i,jmin-1,k)-DP(i,jmin-2,k))*dymi(jmin-1)
!!$     end do
!!$  end do
  
  return
end subroutine velocity_apply_pressure


! ======================================================== !
! Predict the outlet velocity from the continuity equation !
! Assumes dRHO = 0 at outlet - Enforce zero divergence     !
! ======================================================== !
subroutine velocity_predict_outlet
  use velocity
  use parallel
  use data
  use metric_velocity_conv
  use metric_generic
  use masks
  implicit none
  integer  :: i,j,k
  real(WP) :: newU,rhoi
  
  ! Nothing to do if not last proc in x or if periodic
  if (iproc.ne.npx .or. xper.eq.1) return

  !$OMP PARALLEL DO PRIVATE(newU,rhoi)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        if (mask_u(imax+1,j).eq.2) then 
           newU = -1.0_WP / divc_u(imax,j,1) * ( &
                + dRHO(imax,j,k)/dt &
                + sum(divc_u(imax,j,-stc1:0)*rhoU(imax-stc1:imax,j,k)) &
                + sum(divc_v(imax,j,:)*rhoV(imax,j-stc1:j+stc2,k)) &
                + sum(divc_w(imax,j,:)*rhoW(imax,j,k-stc1:k+stc2)) &
                - srcP(imax,j,k)/dt )
           !rhoi = sum(interp_sc_x(imax+1,j,:)*RHOmid(imax-st1:imax+st2,j,k))
           rhoi = RHOmid(imax,j,k)
           do i=imax+1,imaxo
              rhoU(i,j,k) = newU
              U(i,j,k) = newU / rhoi
           end do
        end if
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine velocity_predict_outlet


! ======================================================== !
! Compute the momemtum components from velocity components !
!   -> Assume that U,V,W good even at the BCs              !
!   -> Create rhoU,rhoV,rhoW good also at the boundaries   !
! ======================================================== !
subroutine rho_multiply
  use velocity
  use data
  use metric_generic
  implicit none
  integer :: i,j,k

  !$OMP PARALLEL

  ! rhoU
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoU(i,j,k) = U(i,j,k) * sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
        end do
     end do
  end do
  !$OMP END DO
  
  ! rhoV
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoV(i,j,k) = V(i,j,k) * sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
        end do
     end do
  end do
  !$OMP END DO
  
  ! rhoW
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoW(i,j,k) = W(i,j,k) * sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL

  ! CPU borders and periodicity
  call boundary_update_border(rhoU,'+','ym')
  call boundary_update_border(rhoV,'-','y')
  call boundary_update_border(rhoW,'-','ym')

  return
end subroutine rho_multiply


! ======================================================== !
! Compute the velocity components from momemtum components !
!   -> Assume that rhoU,rhoV,rhoW good even at the BCs     !
!   -> Create U,V,W good also at the boundaries            !
! ======================================================== !
subroutine rho_divide
  use velocity
  use data
  use metric_generic
  implicit none
  integer :: i,j,k

  !$OMP PARALLEL

  ! rhoU
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           U(i,j,k) = rhoU(i,j,k) / sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
        end do
     end do
  end do
  !$OMP END DO
  
  ! rhoV
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           V(i,j,k) = rhoV(i,j,k) / sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
        end do
     end do
  end do
  !$OMP END DO
  
  ! rhoW
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           W(i,j,k) = rhoW(i,j,k) / sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL

  ! CPU borders and periodicity
  call boundary_update_border(U,'+','ym')
  call boundary_update_border(V,'-','y')
  call boundary_update_border(W,'-','ym')

  return
end subroutine rho_divide


! U-momentum component
! --------------------
subroutine velocity_residuals_u
  use velocity
  use data
  use memory
  use metric_velocity_conv
  use metric_velocity_visc
  use metric_generic
  implicit none
  
  ! FX(i,j,k) -> xm(i),ym(j),zm(k)
  ! FY(i,j,k) -> x(i),y(j),zm(k)
  ! FZ(i,j,k) -> x(i),ym(j),z(k)
  real(WP) :: rhs,RHOi
  integer  :: i,j,k,ii,jj,kk,st,n

  !$OMP PARALLEL PRIVATE(rhs,RHOi,i,j,k,n)

  ! Convective part
  !$OMP DO
  do kk=kmin_-stc1,kmax_+stc2
     do jj=jmin_-stc1,jmax_+stc2
        do ii=imin_-stc1,imax_+stc2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           rhoUi(i,j,k) = sum(interp_Ju_xm(i,j,:)*rhoU(i-stc1:i+stc2,j,k))
           
           i = ii; j = jj; k = kk;
           
           rhoVi(i,j,k) = sum(interp_Jv_x(i,j,:)*rhoV(i-stc2:i+stc1,j,k))
           rhoWi(i,j,k) = sum(interp_Jw_x(i,j,:)*rhoW(i-stc2:i+stc1,j,k))
           
        end do
     end do
  end do
  !$OMP END DO
  
  ! Viscous part
  !$OMP DO
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           FX(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_u_x(i,j,:)*U(i-stv1:i+stv2,j,k)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k)) &
                                   + sum(divv_v(i,j,:)*V(i,j-stv1:j+stv2,k)) &
                                   + sum(divv_w(i,j,:)*W(i,j,k-stv1:k+stv2))))
           
           i = ii; j = jj; k = kk;
           
           FY(i,j,k) = &
                + sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k)) * &
                ( sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k)) &
                + sum(grad_v_x(i,j,:)*V(i-stv2:i+stv1,j,k)) )
           
           FZ(i,j,k) = &
                + sum(interp_sc_xz(i,j,:,:)*VISC(i-st2:i+st1,j,k-st2:k+st1)) * &
                ( sum(grad_u_z(i,j,:)*U(i,j,k-stv2:k+stv1)) &
                + sum(grad_w_x(i,j,:)*W(i-stv2:i+stv1,j,k)) )
        end do
     end do
  end do
  !$OMP END DO
  
  ! Residual
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           RHOi = sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
           
           ! Pressure + Viscous terms
           if (.not.compressible) then
              rhs = -sum(grad_Px(i,j,:)*P(i-stc2:i+stc1,j,k))
           else
              rhs = -sum(grad_Px(i,j,:)*0.5_WP*(Pold(i-stc2:i+stc1,j,k)+P(i-stc2:i+stc1,j,k)))
           end if
           rhs = rhs + sum(divv_xx(i,j,:)*FX(i-stv2:i+stv1,j,k)) &
                     + sum(divv_xy(i,j,:)*FY(i,j-stv1:j+stv2,k)) &
                     + sum(divv_xz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) 
           
           ! Convective term
           do st=-stc2,stc1
              n = interp_xx(i,j,st)
              rhs = rhs - divc_xx(i,j,st) * rhoUi(i+st,j,k) * &
                   0.5_WP*(U(i+st+n+1,j,k)+U(i+st-n,j,k))
           end do
           do st=-stc1,stc2              
              n = interp_xy(i,j,st)
              rhs = rhs - divc_xy(i,j,st) * rhoVi(i,j+st,k) * &
                   0.5_WP*(U(i,j+st+n-1,k)+U(i,j+st-n,k))              
              n = interp_xz(i,j,st)
              rhs = rhs - divc_xz(i,j,st) * rhoWi(i,j,k+st) * &
                   0.5_WP*(U(i,j,k+st+n-1)+U(i,j,k+st-n))
           end do
           
           ! Full residual
           ResU(i,j,k) = -2.0_WP*U(i,j,k)+Uold(i,j,k) &
                + ( rhoUold(i,j,k) + dt_uvw*rhs + srcU(i,j,k) ) / RHOi
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL
  
  return
end subroutine velocity_residuals_u


! V-momentum component
! --------------------
subroutine velocity_residuals_v
  use velocity
  use data
  use memory
  use metric_velocity_conv
  use metric_velocity_visc
  use metric_generic
  implicit none
  
  ! FX(i,j,k) -> x(i),y(j),zm(k)
  ! FY(i,j,k) -> xm(i),ym(j),zm(k)
  ! FZ(i,j,k) -> xm(i),y(j),z(k)
  real(WP) :: rhs,RHOi
  integer  :: i,j,k,ii,jj,kk,st,n
  real(WP) :: Acos, Asin
  
  !$OMP PARALLEL PRIVATE(rhs,RHOi,i,j,k,n,Acos,Asin)

  ! Convective part
  !$OMP DO
  do kk=kmin_-stc1,kmax_+stc2
     do jj=jmin_-stc1,jmax_+stc2
        do ii=imin_-stc1,imax_+stc2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           rhoVi(i,j,k) = sum(interp_Jv_ym(i,j,:)*rhoV(i,j-stc1:j+stc2,k))
           
           Fcyl(i,j,k) = &
                - sum(interp_cyl_w_zm(i,j,:)*rhoW(i,j,k-stc1:k+stc2)) * & 
                  sum(interp_cyl_w_zm(i,j,:)*W(i,j,k-stc1:k+stc2))
           
           i = ii; j = jj; k = kk;
           
           rhoUi(i,j,k) = sum(interp_Ju_y(i,j,:)*rhoU(i,j-stc2:j+stc1,k))
           rhoWi(i,j,k) = sum(interp_Jw_y(i,j,:)*rhoW(i,j-stc2:j+stc1,k))
           
        end do
     end do
  end do
  !$OMP END DO
  
  ! Viscous part
  !$OMP DO
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           FY(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_v_y(i,j,:)*V(i,j-stv1:j+stv2,k)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k)) &
                                   + sum(divv_v(i,j,:)*V(i,j-stv1:j+stv2,k)) &
                                   + sum(divv_w(i,j,:)*W(i,j,k-stv1:k+stv2))))
           
           Fcylv(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_w_z(i,j,:)*W(i,j,k-stv1:k+stv2)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k))  &
                                   + sum(divv_v(i,j,:)*V(i,j-stv1:j+stv2,k))  &
                                   + sum(divv_w(i,j,:)*W(i,j,k-stv1:k+stv2))) &
                   + ymi(j)*sum(interpv_cyl_v_ym(i,j,:)*V(i,j-stv1:j+stv2,k)))
           
           i = ii; j = jj; k = kk;
           
           FX(i,j,k) = &
                + sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k)) * &
                ( sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k)) &
                + sum(grad_v_x(i,j,:)*V(i-stv2:i+stv1,j,k)) )
           
           FZ(i,j,k) = &
                + sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k-st2:k+st1)) * &
                ( sum(grad_v_z(i,j,:)*V(i,j,k-stv2:k+stv1)) &
                + sum(grad_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) &
                - yi(j)*sum(interpv_cyl_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) )
        end do
     end do
  end do
  !$OMP END DO
  
  ! Residual
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           RHOi = sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
           
           ! Pressure + Viscous terms
           if (.not.compressible) then
              rhs = -sum(grad_Py(i,j,:)*P   (i,j-stc2:j+stc1,k))
           else
              rhs = -sum(grad_Py(i,j,:)*0.5_WP*(Pold(i,j-stc2:j+stc1,k)+P(i,j-stc2:j+stc1,k)))
           end if
           rhs = rhs + sum(divv_yx(i,j,:)*FX(i-stv1:i+stv2,j,k)) &
                     + sum(divv_yy(i,j,:)*FY(i,j-stv2:j+stv1,k)) &
                     + sum(divv_yz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) 
           
           ! Convective term
           do st=-stc2,stc1
              n = interp_yy(i,j,st)
              rhs = rhs - divc_yy(i,j,st) * rhoVi(i,j+st,k) * &
                   0.5_WP*(V(i,j+st+n+1,k)+V(i,j+st-n,k))
           end do
           do st=-stc1,stc2
              n = interp_yx(i,j,st)
              rhs = rhs - divc_yx(i,j,st) * rhoUi(i+st,j,k) * &
                   0.5_WP*(V(i+st+n-1,j,k)+V(i+st-n,j,k))
              n = interp_yz(i,j,st)
              rhs = rhs - divc_yz(i,j,st) * rhoWi(i,j,k+st) * &
                   0.5_WP*(V(i,j,k+st+n-1)+V(i,j,k+st-n))
           end do
           
           ! Cylindrical term - Convective
           rhs = rhs - ymmi(j)*sum(interp_cyl_F_y(i,j,:)*Fcyl(i,j-stc2:j+stc1,k))
           
           ! Cylindrical term - Viscous
           rhs = rhs - yi(j)*sum(interpv_cyl_F_y(i,j,:)*Fcylv(i,j-stv2:j+stv1,k))
           
           ! Full residual
           ResV(i,j,k) = -2.0_WP*V(i,j,k)+Vold(i,j,k) &
                + ( rhoVold(i,j,k) + dt_uvw*rhs + srcV(i,j,k) ) / RHOi
        end do
     end do
  end do
  !$OMP END DO
  
  ! Morinishi's axis treatment
  if (icyl.eq.1 .and. jproc.eq.1) then
     !$OMP DO
     do i=imin_,imax_
        Acos = 2.0_WP*sum(ResV(i,jmin,:)*cos(zm(kmin:kmax)))/real(nz,WP)
        Asin = 2.0_WP*sum(ResV(i,jmin,:)*sin(zm(kmin:kmax)))/real(nz,WP)
        ResV(i,jmin,:) = Acos*cos(zm(kmin:kmax)) + Asin*sin(zm(kmin:kmax))
     end do
     !$OMP END DO
  end if

  !$OMP END PARALLEL
  
  return
end subroutine velocity_residuals_v


! W-momentum component
! --------------------
subroutine velocity_residuals_w
  use velocity
  use data
  use memory
  use metric_velocity_conv
  use metric_velocity_visc
  use metric_generic
  implicit none
  
  ! FX(i,j,k) -> x(i),ym(j),z(k)
  ! FY(i,j,k) -> xm(i),y(j),z(k)
  ! FZ(i,j,k) -> xm(i),ym(j),zm(k)
  ! Fcyl(i,j,k) -> xm(i),ym(j),zm(k)
  real(WP) :: rhs,RHOi
  integer  :: i,j,k,ii,jj,kk,st,n
  
  !$OMP PARALLEL PRIVATE(rhs,RHOi,i,j,k,n)

  ! Convective part
  !$OMP DO
  do kk=kmin_-stc1,kmax_+stc2
     do jj=jmin_-stc1,jmax_+stc2
        do ii=imin_-stc1,imax_+stc2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           rhoWi(i,j,k) = sum(interp_Jw_zm(i,j,:)*rhoW(i,j,k-stc1:k+stc2))
           
           Fcyl(i,j,k) = &
                - sum(interp_cyl_v_ym(i,j,:)*rhoV(i,j-stc1:j+stc2,k)) * &
                  sum(interp_cyl_w_zm(i,j,:)*W(i,j,k-stc1:k+stc2))
           
           i = ii; j = jj; k = kk;
           
           rhoUi(i,j,k) = sum(interp_Ju_z(i,j,:)*rhoU(i,j,k-stc2:k+stc1))
           rhoVi(i,j,k) = sum(interp_Jv_z(i,j,:)*rhoV(i,j,k-stc2:k+stc1))
           
        end do
     end do
  end do
  !$OMP END DO
  
  ! Viscous part
  !$OMP DO
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           FZ(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_w_z(i,j,:)*W(i,j,k-stv1:k+stv2)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k))  &
                                   + sum(divv_v(i,j,:)*V(i,j-stv1:j+stv2,k))  &
                                   + sum(divv_w(i,j,:)*W(i,j,k-stv1:k+stv2))) &
                   + ymi(j)*sum(interpv_cyl_v_ym(i,j,:)*V(i,j-stv1:j+stv2,k)))
           
           i = ii; j = jj; k = kk;
           
           FX(i,j,k) = &
                + sum(interp_sc_xz(i,j,:,:)*VISC(i-st2:i+st1,j,k-st2:k+st1)) * &
                ( sum(grad_u_z(i,j,:)*U(i,j,k-stv2:k+stv1)) &
                + sum(grad_w_x(i,j,:)*W(i-stv2:i+stv1,j,k)) )
           
           FY(i,j,k) = &
                + sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k-st2:k+st1)) * &
                ( sum(grad_v_z(i,j,:)*V(i,j,k-stv2:k+stv1)) &
                + sum(grad_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) &
                - yi(j)*sum(interpv_cyl_w_y(i,j,:)*W(i,j-stv2:j+stv1,k)) )
        end do
     end do
  end do
  !$OMP END DO
  
  ! Residual
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           RHOi = sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
           
           ! Pressure + Viscous terms
           if (.not.compressible) then
              rhs = -sum(grad_Pz(i,j,:)*P   (i,j,k-stc2:k+stc1))
           else
              rhs = -sum(grad_Pz(i,j,:)*0.5_WP*(Pold(i,j,k-stc2:k+stc1)+Pold(i,j,k-stc2:k+stc1)))
           end if
           rhs = rhs + sum(divv_zx(i,j,:)*FX(i-stv1:i+stv2,j,k)) &
                     + sum(divv_zy(i,j,:)*FY(i,j-stv1:j+stv2,k)) &
                     + sum(divv_zz(i,j,:)*FZ(i,j,k-stv2:k+stv1)) 
           
           ! Convective term
           do st=-stc2,stc1
              n = interp_zz(i,j,st)
              rhs = rhs - divc_zz(i,j,st) * rhoWi(i,j,k+st) * &
                   0.5_WP*(W(i,j,k+st+n+1)+W(i,j,k+st-n))
           end do
           do st=-stc1,stc2
              n = interp_zx(i,j,st)
              rhs = rhs - divc_zx(i,j,st) * rhoUi(i+st,j,k) * &
                   0.5_WP*(W(i+st+n-1,j,k)+W(i+st-n,j,k))
              n = interp_zy(i,j,st)
              rhs = rhs - divc_zy(i,j,st) * rhoVi(i,j+st,k) * &
                   0.5_WP*(W(i,j+st+n-1,k)+W(i,j+st-n,k))
           end do
           
           ! Cylindrical term - Convective
           rhs = rhs + ymi(j)*sum(interp_cyl_F_z(i,j,:)*Fcyl(i,j,k-stc2:k+stc1))
           
           ! Cylindrical term - Viscous
           rhs = rhs + ymi(j)*sum(interpv_cyl_F_ym(i,j,:)*FY(i,j-stv1:j+stv2,k))
           
           ! Full residual
           ResW(i,j,k) = -2.0_WP*W(i,j,k)+Wold(i,j,k) &
                + ( rhoWold(i,j,k) + dt_uvw*rhs + srcW(i,j,k) ) / RHOi
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL
  
  return
end subroutine velocity_residuals_w


! =========================== !
! Monitor the velocity module !
! =========================== !
subroutine velocity_monitor()
  use velocity
  use data
  use parallel
  use masks
  implicit none

  ! Compute min/max
  call parallel_max(abs(U),max_u)
  call parallel_max(abs(V),max_v)
  call parallel_max(abs(W),max_w)
  call parallel_max(abs(P),max_p)
  call velocity_maxdivergence
  
  ! Transfer the values to monitor
  call monitor_select_file('velocity')
  call monitor_set_single_value(1,max_u)
  call monitor_set_single_value(2,max_v)
  call monitor_set_single_value(3,max_w)
  call monitor_set_single_value(4,max_p)
  call monitor_set_single_value(5,max_divg)
  
  return
end subroutine velocity_monitor
