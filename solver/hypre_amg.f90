module hypre_amg
  use pressure
  implicit none
  
  ! Number of cells 
  integer :: npcells
  integer :: npcells_
  integer, dimension(:), pointer :: npcells_proc
  
  ! Unique index for all pressure points
  ! -1 : point is not used
  ! >0 : normal point
  integer, dimension(:,:,:), pointer :: p_index
  
  ! Index of first and last local cell 
  integer :: p_min,p_max
  
  ! Arrays NECESSARY to pass arg to HYPRE
  ! Dont change that
  integer :: nrows, ncols
  integer,  dimension(:), pointer :: cols
  integer,  dimension(:), pointer :: rows
  real(WP), dimension(:), pointer :: values

  ! HYPRE objects
  integer(kind=8), parameter :: HYPRE_PARCSR = 5555
  integer(kind=8) :: matrix
  integer(kind=8) :: rhs,sol
  integer(kind=8) :: solver
  integer(kind=8) :: par_matrix, par_rhs, par_sol
  
end module hypre_amg


! =========================================== !
! Prepare the global index of pressure points !
! Following ex5.c from HYPRE examples         !
! =========================================== !
subroutine hypre_amg_prepare
  use hypre_amg
  use metric_velocity_conv
  use parallel
  use masks
  implicit none
  integer  :: i,j,k,ierr
  integer  :: count
  
  ! Count number of local/global pressure cells
  npcells_ = 0
  !$OMP PARALLEL DO REDUCTION(+:npcells_)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if(mask(i,j).eq.0) npcells_ = npcells_ + 1
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_sum(npcells_,npcells)
  
  ! Create an array with nbr cells per cpu
  allocate(npcells_proc(0:nproc))
  call MPI_allgather(npcells_,1,MPI_INTEGER,npcells_proc(1:nproc),1,MPI_INTEGER,comm,ierr)
  npcells_proc(0) = 0
  do i=1,nproc
     npcells_proc(i) = npcells_proc(i) + npcells_proc(i-1)
  end do
  
  ! Create a single index for all pressure points - default -1
  allocate(p_index(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           p_index(i,j,k) = -1
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Start with an offset
  count = npcells_proc(irank-1)-1
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if(mask(i,j).eq.0) then 
              count = count + 1
              p_index(i,j,k) = count
           end if
        end do
     end do
  end do
  
  ! Take care of periodicity and domain decomposition
  call boundary_update_border_int(p_index,'+','ym')
  
  ! Get min/max
  p_max = -1
  p_min = maxval(p_index)
  !$OMP PARALLEL DO REDUCTION(min:p_min) REDUCTION(max:p_max)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (p_index(i,j,k).ne.-1) then
              p_min = min(p_min,p_index(i,j,k))
              p_max = max(p_max,p_index(i,j,k))
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine hypre_amg_prepare


! ================================== !
! Initialize the HYPRE AMG operators !
! ================================== !
subroutine hypre_amg_init_operator
  use hypre_amg
  use parallel
  use metric_velocity_conv
  implicit none
  integer :: i,j,k,st,ierr
  integer  :: ind,tmpi,nst
  real(WP) :: tmpr
  
  par_matrix = 0
  par_rhs = 0
  par_sol = 0
  
  ! Prepare the HYPRE AMG variables
  call hypre_amg_prepare
  
  ! Create the matrix in HYPRE
  call HYPRE_IJMatrixCreate(comm, p_min, p_max, p_min, p_max, matrix, ierr)
  call HYPRE_IJMatrixSetObjectType(matrix, HYPRE_PARCSR, ierr)
  call HYPRE_IJMatrixInitialize(matrix, ierr)
  
  ! Allocate the temporary arrays
  allocate(rows  (1))
  allocate(cols  (3*2*stcp+1))
  allocate(values(3*2*stcp+1))
  
  ! Fill up the matrix, one row at a time
  nrows = 1

  ! Fill up the matrix coefficients
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           if (p_index(i,j,k).ne.-1) then 
              nst = 0
              
              ! Diagonal
              nst = nst + 1
              cols  (nst) = p_index(i,j,k)
              values(nst) = sum(lap(i,j,k,:,0))
              
              do st = -stcp,+stcp
                 if (st.eq.0) cycle
                 
                 ! Left-Right
                 if (nx.ne.1 .and. .not.fft_x) then
                    if(p_index(i+st,j,k).ne.-1) then                     
                       nst = nst + 1
                       cols  (nst) = p_index(i+st,j,k)
                       values(nst) = lap(i,j,k,1,st)
                    end if
                 end if
                 
                 ! Top-Bottom
                 if (ny.ne.1 .and..not. fft_y ) then
                    if (p_index(i,j+st,k).ne.-1) then
                       nst = nst + 1
                       cols  (nst) = p_index(i,j+st,k)
                       values(nst) = lap(i,j,k,2,st)
                    end if
                 end if
                 
                 ! Front-Back
                 if (nz.ne.1 .and. .not. fft_z) then
                    if (p_index(i,j,k+st).ne.-1) then
                       nst = nst + 1
                       cols  (nst) = p_index(i,j,k+st)
                       values(nst) = lap(i,j,k,3,st)
                    end if
                 end if
                 
              end do
              
              ! Sort the points
              do st=1,nst
                 ind = st + minloc(cols(st:nst),1) - 1
                 tmpr = values(st)
                 values(st) = values(ind)
                 values(ind) = tmpr
                 tmpi = cols(st)
                 cols(st) = cols(ind)
                 cols(ind) = tmpi
              end do
              
              ncols = nst
              rows = p_index(i,j,k)
              call HYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, &
                   cols, values, ierr)

           end if
        end do
     end do
  end do
  
  ! Assemble the matrix
  call HYPRE_IJMatrixAssemble(matrix, ierr)
  call HYPRE_IJMatrixGetObject(matrix, par_matrix, ierr)
  
  ! Prepare RHS
  call HYPRE_IJVectorCreate(comm, p_min, p_max, rhs, ierr)
  call HYPRE_IJVectorSetObjectType(rhs, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(rhs, ierr)
  call HYPRE_IJVectorGetObject(rhs, par_rhs, ierr)
  
  ! Create solution vector
  call HYPRE_IJVectorCreate(comm, p_min, p_max, sol, ierr)
  call HYPRE_IJVectorSetObjectType(sol, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(sol, ierr)
  call HYPRE_IJVectorGetObject(sol, par_sol, ierr)
  
  ! Deallocate temporary arrays and reallocate them
  deallocate(rows)
  deallocate(cols)
  deallocate(values)
  allocate(rows  (npcells_))
  allocate(values(npcells_))
  
  !call HYPRE_IJMatrixPrint(matrix, "IJ.out.A", ierr)

  return
end subroutine hypre_amg_init_operator


! ================================ !
! Prepare the solver for HYPRE AMG !
! ================================ !
subroutine hypre_amg_solver_init
  use hypre_amg
  implicit none
  integer :: ierr
  
  ! Create the solver
  solver = 0
  call HYPRE_BoomerAMGCreate(solver, ierr)
  
  ! Set some parameters
  call HYPRE_BoomerAMGSetPrintLevel(solver, 0, ierr)          ! print solve info + parameters 
  call HYPRE_BoomerAMGSetCoarsenType(solver, 6, ierr)         ! Falgout coarsening
  call HYPRE_BoomerAMGSetRelaxType(solver, 3, ierr)           ! G-S/Jacobi hybrid relaxation
                                                              ! 0 WJ,1 SeGS,3 GS,6 SyGS,9 GE
  call HYPRE_BoomerAMGSetNumSweeps(solver, 1, ierr)           ! Sweeeps on each level
  call HYPRE_BoomerAMGSetMaxLevels(solver, 20, ierr)          ! maximum number of levels
  call HYPRE_BoomerAMGSetMaxIter(solver, max_iter, ierr)      ! maximum nbr of iter
  call HYPRE_BoomerAMGSetTol(solver, cvg, ierr)               ! conv. tolerance
  call HYPRE_BoomerAMGSetStrongThrshld(solver, 0.20_WP, ierr) ! strength threshold
  call HYPRE_BoomerAMGSetInterpType(solver, 6, ierr)          ! extended classical modified interpolation

  ! Now setup
  call HYPRE_BoomerAMGSetup(solver, par_matrix, par_rhs, par_sol, ierr)

  return
end subroutine hypre_amg_solver_init


! ========================= !
! Transfer the RHS to HYPRE !
! ========================= !
subroutine hypre_amg_rhs_transfer
  use hypre_amg
  use masks
  implicit none
  integer :: i,j,k,ierr
  integer :: count
  
  ! Set the RHS
  count = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              count = count + 1
              rows  (count) = p_index(i,j,k)
              values(count) = RP(i,j,k)
           end if
        end do
     end do
  end do
  call HYPRE_IJVectorSetValues(rhs, npcells_, rows, values, ierr)
  call HYPRE_IJVectorAssemble (rhs, ierr)
  call HYPRE_IJVectorGetObject(rhs, par_rhs, ierr)
  !call HYPRE_IJVectorPrint(rhs, "IJ.out.b", ierr)
  
  ! Set the initial guess
  count = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              count = count + 1
              rows  (count) = p_index(i,j,k)
              values(count) = 0.0_WP
           end if
        end do
     end do
  end do
  call HYPRE_IJVectorSetValues(sol, npcells_, rows, values, ierr)
  call HYPRE_IJVectorAssemble (sol, ierr)
  call HYPRE_IJVectorGetObject(sol, par_sol, ierr)
  
  return
end subroutine hypre_amg_rhs_transfer


! ================================ !
! Solve the problem with HYPRE AMG !
! ================================ !
subroutine hypre_amg_solve
  use hypre_amg
  implicit none
  integer :: ierr
  
  ! Now solve
  call HYPRE_BoomerAMGSolve(solver, par_matrix, par_rhs, par_sol, ierr)

  ! Run info - needed logging turned on
  call HYPRE_BoomerAMGGetNumIterations(solver, it_p, ierr)
  call HYPRE_BoomerAMGGetFinalReltvRes(solver, max_resP, ierr)
  
  ! Destroy solver
  !call HYPRE_BoomerAMGDestroy(solver, ierr)
  
  return
end subroutine hypre_amg_solve


! ================================ !
! Transfer the solution from HYPRE !
! ================================ !
subroutine hypre_amg_sol_transfer
  use hypre_amg
  use masks
  implicit none
  integer :: i,j,k,ierr
  integer :: count
  
  call HYPRE_IJVectorGetValues(sol, npcells_, rows, values, ierr)
  count = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              count = count + 1
              DP(i,j,k) = values(count)
           end if
        end do
     end do
  end do
  
  return
end subroutine hypre_amg_sol_transfer
