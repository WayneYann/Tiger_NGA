! ======================================== !
! Interface for LLNL's HYPRE libraries :   !
!  - multigrid :                           !
!         = SMG                            !
!         = PFMG                           !
!  - Krylov-based solvers :                !
!         = PCG                            !
!         = BiCGStab                       !
!         = GMRES                          !
!         = SMG-PCG                        !
!         = PFMG-PCG                       !
!         = SMG-GMRES                      !
!         = PFMG-GMRES                     !
!  - hybrid solvers                        !
! ======================================== !
module hypre
  use pressure
  use parallel
  implicit none
  
  ! Dimension of the problem
  integer :: dim
  
  ! HYPRE objects
  integer(kind=8) :: grid
  integer(kind=8) :: stencil
  integer(kind=8) :: matrix
  integer(kind=8) :: rhs
  integer(kind=8) :: sol
  integer(kind=8) :: solver
  integer(kind=8) :: precond
  
  ! Preconditioner info
  integer :: precond_id,relax_type
  
  ! Operator info
  integer,  dimension(:), pointer :: sten_ind
  real(WP), dimension(:), pointer :: val
  integer :: pressure_stencil_size

  !$OMP THREADPRIVATE(val)
  
  ! Solve in walls?
  logical :: solve_wall
  
end module hypre


! ===================================== !
! Initialization for the HYPRE operator !
! ===================================== !
subroutine hypre_init_operator
  use hypre
  use metric_velocity_conv
  use parser
  implicit none
  
  integer :: i,j,k,count,d,st,ierr
  integer, dimension(3) :: lower,upper
  integer, dimension(3) :: offset
  integer, dimension(6) :: matrix_num_ghost
  integer, dimension(3) :: periodicity_hypre
  character(len=str_medium) :: relax
  
  ! Check here if there are sufficient points for communication_border
  dim = 3
  if (nz.eq.1) dim = 2
  if (ny.eq.1) dim = 1
  
  if (zper.eq.1 .and. nz.lt.nover .and. nz.ne.1) &
       call die('hypre_init_operator: nz has to be greater than nover or equal to 1')
  if (yper.eq.1 .and. ny.lt.nover .and. ny.ne.1) &
       call die('hypre_init_operator: ny has to be greater than nover or equal to 1')
  if (xper.eq.1 .and. nx.lt.nover .and. nx.ne.1) &
       call die('hypre_init_operator: nx has to be greater than nover')
 
  ! Do we solve in the walls
  call parser_read('Pressure solve walls',solve_wall,.false.)
  
  ! Set relaxation scheme for PFMG
  if (trim(pressure_precond).eq.'PFMG' .or. trim(pressure_solver).eq.'PFMG') then
     call parser_read('Pressure relaxation',relax,'weighted Jacobi')
     select case (trim(relax))
     case ('Jacobi')
        relax_type=0
     case ('weighted Jacobi')
        relax_type=1
     case ('sym RBGS')
        relax_type=2
     case ('non-sym RBGS')
        relax_type=3
     case default
        call die('hypre: unknow relaxation type')
     end select
  end if

  ! Build the grid object
  call HYPRE_StructGridCreate(comm,dim,grid,ierr)
  if (.not.solve_wall) then
     !$OMP PARALLEL DO PRIVATE(lower,upper,ierr)
     do k=kbmin_,kbmax_
        do j=jbmin_,jbmax_
           do i=ibmin_,ibmax_
              if (bmask(i,j,k).eq.0) then
                 lower(1) = ib1(i)
                 lower(2) = jb1(j)
                 lower(3) = kb1(k)
                 upper(1) = ib2(i)
                 upper(2) = jb2(j)
                 upper(3) = kb2(k)
                 call HYPRE_StructGridSetExtents(grid,lower(1:dim),upper(1:dim),ierr)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     lower(1) = imin_
     lower(2) = jmin_
     lower(3) = kmin_
     upper(1) = imax_
     upper(2) = jmax_
     upper(3) = kmax_
     call HYPRE_StructGridSetExtents(grid,lower(1:dim),upper(1:dim),ierr)
  end if
  periodicity_hypre = 0
  if (periodicity(1).eq.1) periodicity_hypre(1) = nx
  if (periodicity(2).eq.1) periodicity_hypre(2) = ny
  if (periodicity(3).eq.1) periodicity_hypre(3) = nz
  call HYPRE_StructGridSetPeriodic(grid,periodicity_hypre(1:dim),ierr)
  call HYPRE_StructGridAssemble(grid,ierr)
  
  ! Build the stencil
  pressure_stencil_size = 2*dim*stcp+1
  call HYPRE_StructStencilCreate(dim,pressure_stencil_size,stencil,ierr)
  count=0
  offset=0
  call HYPRE_StructStencilSetElement(stencil,count,offset(1:dim),ierr)
  count=count+1
  do d=1,dim
     do st=-stcp,stcp
        if (st.eq.0) cycle
        offset=0
        offset(d) = st
        call HYPRE_StructStencilSetElement(stencil,count,offset(1:dim),ierr)
        count=count+1
     end do
  end do
  
  ! Build the matrix
  call HYPRE_StructMatrixCreate(comm,grid,stencil,matrix,ierr)
  matrix_num_ghost = stcp
  call HYPRE_StructMatrixSetNumGhost(matrix,matrix_num_ghost(1:2*dim),ierr)
  call HYPRE_StructMatrixInitialize(matrix,ierr)
  allocate(sten_ind(pressure_stencil_size))
  !$OMP PARALLEL
  allocate(val(pressure_stencil_size))
  !$OMP END PARALLEL
  do i=1,pressure_stencil_size
     sten_ind(i) = i-1
  end do
  
  ! Prepare RHS
  call HYPRE_StructVectorCreate(comm,grid,rhs,ierr)
  call HYPRE_StructVectorInitialize(rhs,ierr)
  call HYPRE_StructVectorAssemble(rhs,ierr)
  
  ! Create solution vector
  call HYPRE_StructVectorCreate(comm,grid,sol,ierr)
  call HYPRE_StructVectorInitialize(sol,ierr)
  call HYPRE_StructVectorAssemble(sol,ierr)
  
  return
end subroutine hypre_init_operator


! ========================= !
! Update the HYPRE operator !
! ========================= !
subroutine hypre_update_operator
  use hypre
  use metric_velocity_conv
  use masks
  implicit none
  
  integer :: i,j,k,count,d,st,ierr
  integer, dimension(3) :: ind
  
  ! Build the matrix
  !$OMP PARALLEL DO PRIVATE(ind,count,ierr)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! Mesh position
           ind(1) = i
           ind(2) = j
           ind(3) = k
           ! Matrix coefficients
           count=1
           val(count) = sum(lap(i,j,k,1:dim,0))
           count=count+1
           do d=1,dim
              do st=-stcp,stcp
                 if (st.eq.0) cycle
                 val(count) = lap(i,j,k,d,st)
                 count=count+1
              end do
           end do
           if (solve_wall) then 
              ! Test for walls
              if (mask(i,j).ne.0) then
                 ! Identity
                 count=1
                 val(count) = 1.0_WP
                 count=count+1
                 do d=1,dim
                    do st=-stcp,stcp
                       if (st.eq.0) cycle
                       val(count) = 0.0_WP
                       count=count+1
                    end do
                 end do
              end if
           end if
           ! Set the values
           call HYPRE_StructMatrixSetValues(matrix,ind(1:dim),pressure_stencil_size,sten_ind,val,ierr)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call HYPRE_StructMatrixAssemble(matrix,ierr)
  
  return
end subroutine hypre_update_operator


! ========================================== !
! Initialization of the HYPRE preconditioner !
! ========================================== !
subroutine hypre_precond_init
  use hypre
  implicit none
  
  integer  :: ierr
  integer  :: precond_max_iter
  real(WP) :: precond_cvg
  
  select case(trim(pressure_precond))
  case('PFMG')
     call HYPRE_StructPFMGCreate(comm,precond,ierr)
     precond_max_iter=1
     call HYPRE_StructPFMGSetMaxIter(precond,precond_max_iter,ierr)
     precond_cvg=0.0_WP
     call HYPRE_StructPFMGSetTol(precond,precond_cvg,ierr)
     call HYPRE_StructPFMGSetRelaxType(precond,relax_type,ierr)
     precond_id = 1
  case('SMG')
     call HYPRE_StructSMGCreate(comm,precond,ierr)
     precond_max_iter=1
     call HYPRE_StructSMGSetMaxIter(precond,precond_max_iter,ierr)
     precond_cvg=0.0_WP
     call HYPRE_StructSMGSetTol(precond,precond_cvg,ierr)
     precond_id = 0
  case('DIAGONAL')
     precond = 0
     precond_id = 8
  case ('none')
     precond_id = -1
  case default
     call die('Hypre_precond_init: Unknown HYPRE preconditioner')
  end select
  
  return
end subroutine hypre_precond_init


! ================================== !
! Initialization of the HYPRE solver !
! ================================== !
subroutine hypre_solver_init
  use hypre
  implicit none
  
  integer  :: ierr
  integer  :: log_level
  real(WP) :: conv_tol
  integer  :: dscg_max_iter
  integer  :: pcg_max_iter
  integer  :: solver_type

  ! Update the operator
  call hypre_update_operator
  
  ! Setup the preconditioner
  call hypre_precond_init
  
  ! Setup the solver
  log_level = 1 ! 1 for some info (res,it)
  select case(trim(pressure_solver))
  case('PFMG')
     call HYPRE_StructPFMGCreate(comm,solver,ierr)
     call HYPRE_StructPFMGSetMaxIter(solver,max_iter,ierr)
     call HYPRE_StructPFMGSetTol(solver,cvg,ierr)
     call HYPRE_StructPFMGSetRelaxType(solver,relax_type,ierr)
     call HYPRE_StructPFMGSetLogging(solver,log_level,ierr)
     call HYPRE_StructPFMGSetup(solver,matrix,rhs,sol,ierr)
  case('SMG')
     call HYPRE_StructSMGCreate(comm,solver,ierr)
     call HYPRE_StructSMGSetMaxIter(solver,max_iter,ierr)
     call HYPRE_StructSMGSetTol(solver,cvg,ierr)
     call HYPRE_StructSMGSetLogging(solver,log_level,ierr)
     call HYPRE_StructSMGSetup(solver,matrix,rhs,sol,ierr)
  case('PCG')
     call HYPRE_StructPCGCreate(comm,solver,ierr)
     call HYPRE_StructPCGSetMaxIter(solver,max_iter,ierr)
     call HYPRE_StructPCGSetTol(solver,cvg,ierr)
     call HYPRE_StructPCGSetLogging(solver,log_level,ierr)
     if (precond_id.ge.0) call HYPRE_StructPCGSetPrecond(solver,precond_id,precond,ierr)
     call HYPRE_StructPCGSetup(solver,matrix,rhs,sol,ierr)
  case('BICGSTAB')
     call HYPRE_StructBICGSTABCreate(comm,solver,ierr)
     call HYPRE_StructBICGSTABSetMaxIter(solver,max_iter,ierr)
     call HYPRE_StructBICGSTABSetTol(solver,cvg,ierr)
     call HYPRE_StructBICGSTABSetLogging(solver,log_level,ierr)
     if (precond_id.ge.0) call HYPRE_StructBICGSTABSetPrecond(solver,precond_id,precond,ierr)
     call HYPRE_StructBICGSTABSetup(solver,matrix,rhs,sol,ierr)
  case('GMRES')
     call HYPRE_StructGMRESCreate(comm,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,max_iter,ierr)
     call HYPRE_StructGMRESSetTol(solver,cvg,ierr)
     call HYPRE_StructGMRESSetLogging(solver,log_level,ierr)
     if (precond_id.ge.0) call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     call HYPRE_StructGMRESSetup(solver,matrix,rhs,sol,ierr)
  case ('HYBRID')
     call HYPRE_StructHYBRIDCreate(comm,solver,ierr)
     dscg_max_iter=2*max_iter
     call HYPRE_StructHYBRIDSetDSCGMaxIte(solver,dscg_max_iter,ierr)
     pcg_max_iter=max_iter
     call HYPRE_StructHYBRIDSetPCGMaxIter(solver,pcg_max_iter,ierr)
     call HYPRE_StructHYBRIDSetTol(solver,cvg,ierr)
     conv_tol=0.9_WP
     call HYPRE_StructHYBRIDSetConvergenc(solver,conv_tol,ierr)
     call HYPRE_StructHYBRIDSetLogging(solver,log_level,ierr)
     solver_type=1
     call HYPRE_StructHybridSetSolverType(solver,solver_type,ierr)
     if (precond_id.ge.0) call HYPRE_StructHYBRIDSetPrecond(solver,precond_id,precond,ierr)
     call HYPRE_StructHYBRIDSetup(solver,matrix,rhs,sol,ierr)
  end select
  
  return
end subroutine hypre_solver_init


! ========================= !
! Transfer the RHS to HYPRE !
! ========================= !
subroutine hypre_rhs_transfer
  use hypre
  use masks
  implicit none
  
  integer :: i,j,k,ierr
  integer, dimension(3) :: ind
  real(WP) :: value
  
  !$OMP PARALLEL DO PRIVATE(ind,value,ierr)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              ind(1) = i
              ind(2) = j
              ind(3) = k
              value  = RP(i,j,k)
              call HYPRE_StructVectorSetValues(rhs,ind(1:dim),value,ierr)
           else if (solve_wall) then
              ind(1) = i
              ind(2) = j
              ind(3) = k
              value  = 0.0_WP
              call HYPRE_StructVectorSetValues(rhs,ind(1:dim),value,ierr)
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call HYPRE_StructVectorAssemble(rhs,ierr)
  
  return
end subroutine hypre_rhs_transfer


! ============================ !
! Solve the problem with HYPRE !
! ============================ !
subroutine hypre_solve
  use hypre
  implicit none
  
  integer :: ierr
  
  ! Solve
  select case(trim(pressure_solver))
  case('PFMG')
     call HYPRE_StructPFMGSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructPFMGGetNumIteration(solver,it_p,ierr)
     call HYPRE_StructPFMGGetFinalRelativ(solver,max_resP,ierr)
  case('SMG')
     call HYPRE_StructSMGSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructSMGGetNumIterations(solver,it_p,ierr)
     call HYPRE_StructSMGGetFinalRelative(solver,max_resP,ierr)
  case ('PCG')
     call HYPRE_StructPCGSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructPCGGetNumIterations(solver,it_p,ierr)
     call HYPRE_StructPCGGetFinalRelative(solver,max_resP,ierr)
  case ('BICGSTAB')
     call HYPRE_StructBICGSTABSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructBICGSTABGetNumItera(solver,it_p,ierr)
     call HYPRE_StructBICGSTABGetFinalRel(solver,max_resP,ierr)
  case ('GMRES')
     call HYPRE_StructGMRESSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructGMRESGetNumIteratio(solver,it_p,ierr)
     call HYPRE_StructGMRESGetFinalRelati(solver,max_resP,ierr)
  case ('HYBRID')
     call HYPRE_StructHYBRIDSolve(solver,matrix,rhs,sol,ierr)
     call HYPRE_StructHYBRIDGetNumIterati(solver,it_p,ierr)
     call HYPRE_StructHYBRIDGetFinalRelat(solver,max_resP,ierr)
  end select
  
  return
end subroutine hypre_solve


! ================================ !
! Transfer the solution from HYPRE !
! ================================ !
subroutine hypre_sol_transfer
  use hypre
  use masks
  implicit none
  
  integer :: i,j,k,ierr
  integer, dimension(3) :: ind
  real(WP) :: value
  
  !$OMP PARALLEL DO PRIVATE(ind,value,ierr)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              ind(1) = i
              ind(2) = j
              ind(3) = k
              call HYPRE_StructVectorGetValues(sol,ind(1:dim),value,ierr)
              DP(i,j,k) = value
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine hypre_sol_transfer
