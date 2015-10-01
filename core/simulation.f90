module simulation
  use time_info
  implicit none

  ! Formulation: Compressible or low Mach
  logical :: compressible

  ! Constant density
  logical :: constant_density
  
end module simulation

! ================================ !
! Initialization of the simulation !
! ================================ !
subroutine simulation_init
  use simulation
  use parser
  use partition
  use string
  implicit none

  character(len=str_medium) :: tmp_string

  ! Initialize the monitoring system
  call monitor_pre_init

  ! Initialize the time parameters
  call time_info_init

  ! Problem formulation
  call parser_read('Compressible',compressible,.false.)
  call parser_read('Chemistry model',tmp_string,'none')
  if (trim(tmp_string).eq.'none') constant_density = .true.

  ! Initialize the modules with the previous time step
  call boundary_init
  call scalar_init
  call combustion_init
  call spray_init
  call pollutants_init
  call velocity_init

  ! Pre-interpolate the velocities
  call interpolate_init
  call interpolate_velocities
  call strainrate_init
  
  ! Additional modules
  call bodyforce_init
  
  ! Initialize the pressure
  if (.not.compressible) call pressure_init
  
  ! Initialize SGS models
  call sgsmodel_init
  
  ! Output stuff
  call dump_init
  call stat_init
  call inflow_generation_init
  
  ! Monitor the initial data
  call monitor_post_init
  call monitor_timestep
  
  ! Log
  call monitor_log("SIMULATION INITIALIZED")
  
  return
end subroutine simulation_init


! =========================== !
! Simulation time advancement !
! =========================== !
subroutine simulation_run
  use simulation
  implicit none
  
  do while (.not. done())
     
     ! Increment time information
     call predict_timestepsize
     ntime = ntime + 1
     time  = time  + dt
     
     ! Pre-step routines 1 - physical properties
     call combustion_prestep
     call spray_prestep
     call pollutants_prestep
     call velocity_CFL_centerline
     
     ! Pre-step routines 2 - prepare the subiterations
     call scalar_prestep
     call velocity_prestep
     call velocity_predict_density
     
     !! New subiteration loop for compressible?
     do niter=1,max_niter

        ! Subfilter diffusivity
        call sgsmodel_eddyDIFF
        
        ! INTRODUCE STRANG SPLITTING
        ! Scalar transport equations
        call scalar_step

        ! Equation of state
        call combustion_step

        ! Continuity source terms
        !call spray_step
        !call soot_step
        
        ! Subfilter viscosity
        call sgsmodel_eddyVISC

        ! Momentum equations
        call velocity_step
        
        ! Continuity equation
        if (.not.compressible) call velocity_pressure
        
        ! Monitoring - per sub-iterations
        call monitor_iteration
        
     end do

     ! Advance soot here
!!$     call soot_step
        
     ! Monitoring - per iterations
     call monitor_timestep
     
     ! Write out data and results
     call interpolate_velocities
     call simulation_write(.false.)
     call dump_result
     call dump_statistics
     call inflow_generation_save
     
  end do
  
end subroutine simulation_run


! ==================== !
! Write solution files !
! ==================== !
subroutine simulation_write(flag)
  use simulation
  implicit none
  logical, intent(in) :: flag
  
  call data_write(flag)
  call optdata_write(flag)
  
  return
end subroutine simulation_write


! ================== !
! End the simulation !
! ================== !
subroutine simulation_finalize
  use simulation
  implicit none
  
  ! Write necessary data to files
  call data_finalize
  call optdata_finalize
  
  ! Log
  call monitor_log("END OF SIMULATION")
  
  ! Finalize the monitor module
  call monitor_finalize
  
  return
end subroutine simulation_finalize
