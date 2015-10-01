module scalar
  use precision
  use partition
  use geometry
  use time_info
  implicit none
  
  ! Solution at old time step : n
  real(WP), dimension(:,:,:,:), pointer :: SCold
  
  ! Source term between time n and n+1
  real(WP), dimension(:,:,:,:), pointer :: srcSCmid
  real(WP), dimension(:,:,:,:), pointer :: srcSCfull
  
  ! Residual between time n and n+1
  real(WP), dimension(:,:,:,:), pointer :: ResSC
  
  ! Values to monitor
  real(WP), dimension(:), pointer :: max_resSC
  real(WP), dimension(:), pointer :: ext_SC
  
end module scalar


! ==================================================== !
! Initialize the scalar module                         !
!                                                      !
! -> allocate the arrays                               !
! -> update ghost cells                                !
! -> apply boundary conditions                         !
! -> run specific init                                 !
!                                                      !
! Before: SC correct only inside the domain            !
!         -> imin_:imax_,jmin_:jmax_,kmin_:kmax_       !
! After : SC correct everywhere                        !
!         -> imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_ !
! ==================================================== !
subroutine scalar_init
  use scalar
  use parser
  use data
  use implicit
  implicit none
  
  integer :: isc
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Create & Start the timer
  call timing_create('scalar')
  call timing_start ('scalar')
  
  ! Allocate arrays for old solution
  allocate(SCold    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(srcSCmid (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(srcSCfull(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(ResSC    (imin_ :imax_ ,jmin_ :jmax_ ,kmin_ :kmax_ ,nscalar))
  
  ! The scalars fields read from the file were without ghost cells
  ! Update the ghost cells for periodicity and domain decomposition
  do isc=1,nscalar
     call boundary_update_border(SC(:,:,:,isc),'+','ym')
  end do
  
  ! Apply boundary conditions
  call boundary_scalar_dirichlet
  call boundary_scalar_outflow
  call boundary_scalar_neumann
  
  ! Initialize the given scheme
  select case (trim(scalar_scheme))
  case ('quick')
     call scalar_quick_init
  case ('bquick')
     call scalar_bquick_init
  case ('weno3')
     call scalar_weno3_init
  case ('weno5')
     call scalar_weno5_init
  case ('houc')
     call scalar_houc_init
  case default
     call die('Unknown scalar scheme specified')
  end select
  
  ! Create new file to monitor at each iterations
  call monitor_create_file_step('scalar',2*nscalar)
  allocate(ext_SC(2*nscalar))
  do isc=1,nscalar
     call monitor_set_header(2*isc-1,'min_'//SC_name(isc),'r')
     call monitor_set_header(2*isc+0,'max_'//SC_name(isc),'r')
  end do
  
  ! Create new file to monitor at each subiterations
  call monitor_create_file_iter('convergence_scalar',nscalar)
  allocate(max_resSC(nscalar))
  do isc=1,nscalar
     call monitor_set_header(isc,'res_'//SC_name(isc),'r')
  end do
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_init


! ================================================== !
! PRE-TIMESTEP Routine                               !
!                                                    !
! -> Set up the iterative process                    !
! -> Compute the source term for the scalar equation !
! ================================================== !
subroutine scalar_prestep
  use scalar
  use data
  implicit none
  integer :: i,j,k,isc
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')
  
  ! Compute Dirichlet BC
  call boundary_scalar_dirichlet
  
  !$OMP PARALLEL

  ! Save the old velocity
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SCold(i,j,k,isc) = SC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
     
  ! Zero the source term
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              srcSCmid (i,j,k,isc) = 0.0_WP
              srcSCfull(i,j,k,isc) = 0.0_WP
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL
  
  ! Precompute some of the source terms
  call combustion_source_scalar_full(SC,srcSCfull)
  
  ! Source term from SGS
  call sgsmodel_src_sc(srcSCfull)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_prestep


! ========================================================== !
! ADVANCE the solution                                       !
!   -> second order in time                                  !
!   -> variable accuracy in space                            !
!   -> explicit prediction                                   !
!                                                            !
! Z(n+3/2,k+1) = Z(n+1/2) + dt*F(0.5*(Z(n+3/2,k)+Z(n+1/2))   !
!                   + 0.5*dt*dF/dZ*(Z(n+3/2,k+1)-Z(n+3/2,k)) !
! n : time step                                              !
! k : inner loop iteration                                   !
! Velocity field used : best approximation for U(n+1)        !
! ========================================================== !
subroutine scalar_step
  use scalar
  use data
  implicit none
  integer :: i,j,k,isc
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')

  !$OMP PARALLEL
  
  ! Compute mid point
  ! Store it in the 'n+1/2' scalar
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 0.5_WP*(SC(i,j,k,isc)+SCold(i,j,k,isc))
           end do
        end do
     end do
     !$OMP END DO
  end do
  
  ! Compute the combustion source terms
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              srcSCmid(i,j,k,isc) = 0.0_WP
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL

  call combustion_source_scalar_mid(SC,srcSCmid)
  
  select case (trim(scalar_scheme))
  case ('quick')
     call scalar_quick_residual
     call scalar_quick_inverse
  case ('bquick')
     call scalar_bquick_residual
     call scalar_bquick_inverse
  case ('weno3')
     call scalar_weno3_residual
     call scalar_weno3_inverse
  case ('weno5')
     call scalar_weno5_residual
     call scalar_weno5_inverse
  case ('houc')
     call scalar_houc_residual
     call scalar_houc_inverse
  end select
  
  !$OMP PARALLEL

  ! Update the scalars
  do isc=1,nscalar
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              SC(i,j,k,isc) = 2.0_WP*SC(i,j,k,isc)-SCold(i,j,k,isc) + ResSC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
  
  !$OMP END PARALLEL
  
  ! Update the physical boundaries
  call boundary_scalar_neumann
  call boundary_scalar_outflow
  
  ! Update the overlapped cells
  do isc=1,nscalar
     call boundary_update_border(SC(:,:,:,isc),'+','ym')
  end do
  
  ! Compute max of residuals
  do isc=1,nscalar
     call parallel_max(maxval(abs(resSC(:,:,:,isc))),max_resSC(isc))
  end do
  
  ! Transfer values to monitor
  call monitor_select_file('convergence_scalar')
  call monitor_set_array_values(max_resSC)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_step


! =================== !
! Monitor the scalars !
! =================== !
subroutine scalar_monitor
  use scalar
  use data
  implicit none
  integer :: isc
  
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')
  
  ! Get min/max
  do isc=1,nscalar
     call parallel_max(-SC(:,:,:,isc),ext_SC(2*isc-1))
     ext_SC(2*isc-1) = -ext_SC(2*isc-1)
     call parallel_max( SC(:,:,:,isc),ext_SC(2*isc+0))
  end do
  
  ! Transfer values to monitor
  call monitor_select_file('scalar')
  call monitor_set_array_values(ext_SC)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_monitor
