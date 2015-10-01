! ======================================== !
! Pressure solver for low-Mach number flow !
! LLNL's HYPRE libraries :                 !
!  - multigrid :                           !
!         = SMG                            !
!         = PFMG                           !
!  - Krylov-based solvers :                !
!         = PCG                            !
!         = BICGSTAB                       !
!         = GMRES                          !
!         = SMG-PCG                        !
!         = PFMG-PCG                       !
!         = SMG-GMRES                      !
!         = PFMG-GMRES                     !
!         = SMG-BICGSTAB                   !
!         = PFMG-BICGSTAB                  !
!  - hybrid solvers                        !
!                                          !
! Homemade solvers :                       !
!  - bicgstab                              !
!  - bicgstab(2)                           !
! ======================================== !
module pressure
  use precision
  use geometry
  use partition
  use structure
  use string
  use config
  implicit none
  
  ! RHS arrays
  real(WP), dimension(:,:,:), pointer :: RP
  
  ! Solution array
  real(WP), dimension(:,:,:), pointer :: DP
  
  ! Solver parameters
  character(len=str_medium) :: pressure_solver
  character(len=str_medium) :: pressure_precond
  real(WP) :: cvg
  integer  :: max_iter
  logical  :: use_HYPRE
  logical  :: use_HYPRE_AMG

  ! Array of diagonal values for ICC
  real(WP), dimension(:,:,:), pointer :: diag
  
  ! Values to monitor
  real(WP) :: max_resP,int_RP
  integer  :: it_p

  ! Direction for fft transform
  logical :: fft_x
  logical :: fft_y
  logical :: fft_z

end module pressure


! ============================== !
! Initialize the pressure solver !
! ============================== !
subroutine pressure_init
  use pressure
  use parser
  implicit none
  
  integer :: i,j,k

  ! Create & Start the timer
  call timing_create('pressure')
  call timing_start ('pressure')

  ! Allocate the RHS and DP
  allocate(RP(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(DP(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  !$OMP PARALLEL

  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           RP(i,j,k)=0.0_WP
        end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           DP(i,j,k)=0.0_WP
        end do
     end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL
  
  ! Initialize HYPRE structured solvers from input file
  call parser_read("Pressure solver",pressure_solver)
  call parser_read("Pressure precond",pressure_precond,'none')
  call parser_read("Pressure cvg",cvg)
  call parser_read("Pressure iterations",max_iter)
  
  ! Initialize the fft routines
  call fourier_init

  ! Rescale the operator
  call fourier_operator
  call pressure_rescale
  
  ! Detection of the solver
  use_HYPRE    =.false.
  use_HYPRE_AMG=.false.
  select case (trim(pressure_solver))
     ! HYPRE multigrid solvers
  case ('AMG')
     use_HYPRE_AMG=.true.
  case ('SMG')
     use_HYPRE=.true.
  case ('PFMG')
     use_HYPRE=.true.
  case ('PCG')
     use_HYPRE=.true.
  case ('BICGSTAB')
     use_HYPRE=.true.
  case ('GMRES')
     use_HYPRE=.true.
  case ('HYBRID')
     use_HYPRE=.true.
  ! Homemade BiCGStab solver
  case ('bicgstab')
     use_HYPRE=.false.
     call bicgstab_init
  case ('bicgstab(2)')
     use_HYPRE=.false.
     call bicgstab2_init
  ! Unknown solver
  case default
     call die('Pressure_init: unknown pressure solver')
  end select
  
  ! Initialization of Hypre
  if (use_HYPRE) then
     if (vel_conv_order.gt.2) call die('HYPRE requires a 2nd order accurate velocity scheme.')
     call hypre_init_operator
     call hypre_solver_init
  end if
  ! Initialization of Hypre - AMG
  if (use_HYPRE_AMG) then
     call hypre_amg_init_operator
     call hypre_amg_solver_init
  end if
  
  ! Initialization of the preconditioner
  if (.not.use_HYPRE) then
     select case(trim(pressure_precond))
     case('none')
        ! Nothing to do
     case('diag')
        ! Nothing to do
     case('tridiag')
        ! Nothing to do
     case('icc')
        call pressure_icc_init
        call pressure_icc_diag
     case default
        call die('Unknown pressure preconditioner')
     end select
  end if
  
  ! Create new file to monitor
  call monitor_create_file_iter('convergence_pressure',3)
  call monitor_set_header (1,'it_p', 'i')
  call monitor_set_header (2,'res_p','r')
  call monitor_set_header (3,'int_RP','r')

  ! Stop the timer
  call timing_stop('pressure')

  return
end subroutine pressure_init


! ================================================== !
! Compute the scaling coefficients for the laplacian !
! Ensure that the matrix is symmetric                !
! ================================================== !
subroutine pressure_rescale
  use pressure
  use metric_velocity_conv
  implicit none
  integer  :: i,j,k,st

  ! Apply scaling
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           lap(i,j,k,:,:) = -vol(i,j)*lap(i,j,k,:,:)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  do st=-stcp,stcp
     call boundary_update_border(lap(:,:,:,1,st),'+','ym')
     call boundary_update_border(lap(:,:,:,2,st),'+','ym')
     call boundary_update_border(lap(:,:,:,3,st),'+','ym')
  end do
  
  return
end subroutine pressure_rescale


! ============================= !
! Solve the system using HYPRE  !
! P(n+1/2) = P(n-1/2) + DP(n+1) !
! ============================= !
subroutine pressure_step
  use pressure
  use data
  implicit none
  real(WP) :: my_int
  integer :: i,j,k

  ! Initialize the solution
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           DP(i,j,k) = 0.0_WP
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Zero iteration => no pressure solver
  if (max_iter.eq.0) return
  
  ! Start the timer
  call timing_start('pressure')
  
  ! Check the integral of RP
  my_int = 0.0_WP
  !$OMP PARALLEL DO REDUCTION(+:my_int)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           my_int = my_int+RP(i,j,k)*vol(i,j)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_sum(my_int,int_RP)
  int_RP = int_RP/vol_total
  
  ! Prepare the RHS
  call pressure_RHS
  
  ! Solve
  call pressure_SOLVE
  
  ! Update boundaries of DP
  call pressure_DP_BC
  call boundary_update_border(DP,'+','ym')
  
  ! Update the pressure
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           P(i,j,k) = P(i,j,k) + DP(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Pass the values to monitor
  call monitor_select_file('convergence_pressure')
  call monitor_set_single_value(1,real(it_p,WP))
  call monitor_set_single_value(2,max_resP)
  call monitor_set_single_value(3,int_RP)

  ! Stop the timer
  call timing_stop('pressure')

  return
end subroutine pressure_step


! ============================== !
! Transfer the RHS to the solver !
! ============================== !
subroutine pressure_RHS
  use pressure
  implicit none
  integer :: i,j,k

  ! Preprocess the RHS with fft
  call fourier_rhs

  ! Rescale the RHS
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           RP(i,j,k) = -vol(i,j)*RP(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! If necessary, transfer to HYPRE
  if (use_HYPRE) call hypre_rhs_transfer
  if (use_HYPRE_AMG) call hypre_amg_rhs_transfer

  return
end subroutine pressure_RHS


! ========================== !
! Solve the Poisson equation !
! ========================== !
subroutine pressure_SOLVE
  use pressure
  use parallel
  use infnan
  implicit none
  
  ! Call the solver
  if (use_HYPRE) then
     call hypre_solve
     call hypre_sol_transfer
  elseif (use_HYPRE_AMG) then
     call hypre_amg_solve
     call hypre_amg_sol_transfer
  else
     select case (trim(pressure_solver))
     case ('bicgstab')
        call bicgstab_solve
     case ('bicgstab(2)')
        call bicgstab2_solve
     end select
  end if
  
  ! Postprocess the RHS with fft
  call fourier_dp
  
  return
end subroutine pressure_SOLVE


! ========================================= !
! Pressure BC on the pressure correction DP !
! ========================================= !
subroutine pressure_DP_BC
  use pressure
  use parallel
  use masks
  implicit none
  integer :: i,j
  real(WP) :: DPmean,DPtmp
  
  ! Fix DPmean = 0
  DPtmp = 0.0_WP
  !$OMP PARALLEL DO REDUCTION(+:DPtmp)
  do j=jmin_,jmax_
     do i=imin_,imax_
        DPtmp = DPtmp + vol(i,j)*sum(DP(i,j,kmin_:kmax_))
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_sum(DPtmp,DPmean)
  DPmean = DPmean/vol_total
  !$OMP PARALLEL DO
  do j=jmin_,jmax_
     do i=imin_,imax_
        if (mask(i,j).eq.0) DP(i,j,:) = DP(i,j,:) - DPmean
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine pressure_DP_BC


! ============================= !
! Pressure operator: B = lap(A) !
! A is with overlap             !
! B is interior only            !
! COMMUNICATION OF A DONE HERE  !
! ============================= !
subroutine pressure_operator(A,B)
  use pressure
  use metric_velocity_conv
  implicit none
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: A
  real(WP), dimension(imin_ :imax_ ,jmin_ :jmax_ ,kmin_ :kmax_ ) :: B
  integer :: i,j,k

  ! Communicate A
  call boundary_update_border(A,'+','ym')
  
  ! Apply operator
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           B(i,j,k) = &
                sum(A(i-stcp:i+stcp,j,k)*lap(i,j,k,1,:)) + & ! X-direction
                sum(A(i,j-stcp:j+stcp,k)*lap(i,j,k,2,:)) + & ! Y-direction
                sum(A(i,j,k-stcp:k+stcp)*lap(i,j,k,3,:))     ! Z-direction
           
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine pressure_operator


! ================================================ !
! Initialize the Incomplete Cholesky factorization !
! ================================================ !
subroutine pressure_icc_init
  use pressure
  implicit none
  
  ! Allocate the array for the diagonal values
  allocate(diag(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  return
end subroutine pressure_icc_init


! ============================================= !
! Prepare the Incomplete Cholesky factorization !
! ============================================= !
subroutine pressure_icc_diag
  use pressure
  use metric_velocity_conv
  use masks
  implicit none
  integer  :: i,j,k,st,n
  
  ! Preset diagonal
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (mask(i,j).eq.0) then
              diag(i,j,k) = sum(lap(i,j,k,:,0))
           else
              diag(i,j,k) = 1.0_WP
           end if
        end do
     end do
  end do
  call boundary_update_border(diag,'+','ym')  
  
  ! Compute this diagonal
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           do st=-stcp,-1
              !n=-st
              do n=1,stcp
                 diag(i,j,k) = diag(i,j,k) &
                      - lap(i,j,k,1,st)*lap(i+st,j,k,1,n)/diag(i+st,j,k) &
                      - lap(i,j,k,2,st)*lap(i,j+st,k,2,n)/diag(i,j+st,k) &
                      - lap(i,j,k,3,st)*lap(i,j,k+st,3,n)/diag(i,j,k+st)
              end do
           end do
        end do
     end do
  end do

  return
end subroutine pressure_icc_diag


! =========================================================== !
! Pressure preconditioner : Incomplete Cholesky factorization !
! K = (D+L)*D^1*(D+U)                                         !
! with A = L+diag(A)+U,   U=transpose(L)                      !
! =========================================================== !
subroutine pressure_precond_icc(A,B)
  use pressure
  use metric_velocity_conv
  implicit none

  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: B
  integer :: i,j,k,st
  
  ! Communicate everything first
  B(imin_:imax_,jmin_:jmax_,kmin_:kmax_) = A
  call boundary_update_border(B,'+','ym')
  
  ! Lower triangular
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           do st=-stcp,-1
              B(i,j,k) =  B(i,j,k) &
                   - lap(i,j,k,1,st)*B(i+st,j,k) &
                   - lap(i,j,k,2,st)*B(i,j+st,k) &
                   - lap(i,j,k,3,st)*B(i,j,k+st) 
           end do
           B(i,j,k) = B(i,j,k) / diag(i,j,k)
        end do
     end do
  end do
  
  ! Diagonal
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           B(i,j,k) = diag(i,j,k)*B(i,j,k) 
        end do
     end do
  end do
  
  ! Communicate everything
  call boundary_update_border(B,'+','ym')

  ! Upper triangular
  do i=imax_,imin_,-1
     do j=jmax_,jmin_,-1
        do k=kmax_,kmin_,-1
           do st=1,stcp
              B(i,j,k) =  B(i,j,k) &
                   - lap(i,j,k,1,st)*B(i+st,j,k) &
                   - lap(i,j,k,2,st)*B(i,j+st,k) &
                   - lap(i,j,k,3,st)*B(i,j,k+st) 
           end do
           B(i,j,k) = B(i,j,k) / diag(i,j,k)
        end do
     end do
  end do
  
  return
end subroutine pressure_precond_icc


! ================================== !
! Pressure preconditioner : Diagonal !
! ================================== !
subroutine pressure_precond_diag(A,B)
  use pressure
  use metric_velocity_conv
  use masks
  implicit none
  
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: B
  integer :: i,j,k
  
  ! Diagonal
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              B(i,j,k) = A(i,j,k) / sum(lap(i,j,k,:,0))
           else
              B(i,j,k) = 0.0_WP
           end if
        end do
     end do
  end do
  
  return
end subroutine pressure_precond_diag


! ===================================== !
! Pressure preconditioner : tridiagonal !
! ===================================== !
subroutine pressure_precond_tridiag(A,B)
  use pressure
  use metric_velocity_conv
  use metric_generic
  use masks
  use memory
  use parallel
  implicit none
  
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: A
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: B
  integer  :: i,j,k,st
  
  ! Store initial solution
  B(imin_:imax_,jmin_:jmax_,kmin_:kmax_) = A
  call boundary_update_border(B,'+','ym')
  
  ! Select direction
  if (.not.fft_x .and. nx.ne.1) then
     
     ! X-direction
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).eq.0) then
                 do st=-stcp,stcp
                    Ax(j,k,i,st) = lap(i,j,k,1,st)
                 end do
                 Ax(j,k,i,0) = Ax(j,k,i,0) + lap(i,j,k,2,0) + lap(i,j,k,3,0)
                 Rx(j,k,i) = B(i,j,k)
              else
                 Ax(j,k,i,:) = 0.0_WP
                 Ax(j,k,i,0) = 1.0_WP
                 Rx(j,k,i) = 0.0_WP
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call implicit_solve_x(2*stcp+1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              B(i,j,k) = Rx(j,k,i)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call boundary_update_border(B,'+','ym')
     
  else if (.not.fft_y .and. ny.ne.1) then
     
     ! Y-direction
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).eq.0) then
                 do st=-stcp,stcp
                    Ay(i,k,j,st) = lap(i,j,k,2,st)
                 end do
                 Ay(i,k,j,0) = Ay(i,k,j,0) + lap(i,j,k,1,0) + lap(i,j,k,3,0)
                 Ry(i,k,j) = B(i,j,k)
              else
                 Ay(i,k,j,:) = 0.0_WP
                 Ay(i,k,j,0) = 1.0_WP
                 Ry(i,k,j) = 0.0_WP
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call implicit_solve_y(2*stcp+1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              B(i,j,k) = Ry(i,k,j)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call boundary_update_border(B,'+','ym')
     
  else if (.not.fft_z .and. nz.ne.1) then
     
     ! Z-direction
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).eq.0) then
                 do st=-stcp,stcp
                    Az(i,j,k,st) = lap(i,j,k,3,st)
                 end do
                 Az(i,j,k,0) = Az(i,j,k,0) + lap(i,j,k,1,0) + lap(i,j,k,2,0)
                 Rz(i,j,k) = B(i,j,k)
              else
                 Az(i,j,k,:) = 0.0_WP
                 Az(i,j,k,0) = 1.0_WP
                 Rz(i,j,k) = 0.0_WP
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call implicit_solve_z(2*stcp+1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              B(i,j,k) = Rz(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call boundary_update_border(B,'+','ym')
     
  end if
  
  return
end subroutine pressure_precond_tridiag
