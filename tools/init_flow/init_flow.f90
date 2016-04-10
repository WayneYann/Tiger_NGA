program init_flow
  use string
  use parser
  use param
  implicit none

  character(len=str_medium) :: simulation
  character(len=str_medium) :: input_name
  
  ! Parse the command line
  call commandline_args(input_name)
  
  ! Initialize the parser
  call parser_init()
  
  ! Read the input file
  call parser_parsefile(input_name)
  
  ! Detect the config
  call parser_readchar('Simulation',simulation)
  
  select case (trim(simulation))
  case ('Taylor vortex')
     call taylor_grid()
     call taylor_data()
  case ('von Karman')
     call vonkarman_grid()
     call vonkarman_data()
  case ('vortex convection')
     call conv_vortex_grid
     call conv_vortex_data
  case ('scalar convection')
     call conv_scalar_grid
     call conv_scalar_data
  case('isotropic turbulence')
     call dns_box_grid()
     call dns_box_data()
     call dns_box_chemtable()
  case ('channel')
     call channel_grid()
     call channel_data()
  case ('pipe')
     call pipe_grid()
     call pipe_data()
  case ('cylindrical jet')
     call jetcyl_grid()
     call jetcyl_data()
  case ('cylindrical coflow')
     call jetcyl_coflow_grid
     call jetcyl_coflow_data
  case ('cartesian jet')
     call jetcart_grid
     call jetcart_data
  case ('concentric cartesian jet')
     call jetcart2_grid
     call jetcart2_data
!!$  case ('RATS slot jet')
!!$     call RATSjet_grid
!!$     call RATSjet_data
  case ('mixing layer')
     call mixing_grid
     call mixing_data
  case ('spatial mixing layer')
     call spatial_mixing_grid
     call spatial_mixing_inflow
     call spatial_mixing_data
  case ('laminar flame')
     call laminar_flame_grid
     call laminar_flame_data
  case ('Lamb vortex')
     call lamb_vortex_grid
     call lamb_vortex_data
  case ('density vortex')
     call density_vortex_grid
     call density_vortex_data
     call density_vortex_chemtable
  case ('boundary layer')
     call boundary_layer_grid
     call boundary_layer_data
  case ('gfm test')
     call gfm_test_grid
     call gfm_test_data
  case ('Marmottant')
     call marmottant_grid
     call marmottant_data
  case ('milk crown')
     call milk_crown_grid
     call milk_crown_data
  case ('deformation')
     call deformation_grid
     call deformation_data
  case ('Rayleigh-Taylor')
     call rayleigh_taylor_grid
     call rayleigh_taylor_data
     call rayleigh_taylor_chemtable
  case ('cartesian round jet')
     call round_cart_jet_grid
     call round_cart_jet_data
  case ('diesel DNS')
     call diesel_dns_grid
     call diesel_dns_data
  case ('energy conservation')
     call energy_cons_grid
     call energy_cons_data
  case ('vortex ring')
     call vortex_ring_grid
     call vortex_ring_data
  case ('Rayleigh-Taylor vortex')
     call rt_vortex_grid
     call rt_vortex_data
  case ('Rayleigh instability')
     call rayleigh_grid
     call rayleigh_data
  case ('wave')
     call wave_grid
     call wave_data
  case ('spurious currents')
     call spurious_currents_grid
     call spurious_currents_data
  case ('curvature')
     call curvature_grid
     call curvature_data
  case ('detailed ignition chemistry')
     call ignition_grid
     call ignition_data
  case ('Helium Plume')
     call helium_plume_grid
     call helium_plume_data
     call helium_plume_optdata
  case ('blowing')
     call blowing_grid
     call blowing_data
  case ('air layer')
     call air_layer_grid
     call air_layer_data
  case ('grating')
     call grating_grid
     call grating_data
  case ('bubble')
     call bubble_grid
     call bubble_data
  case ('drop collision')
     call drop_coll_grid
     call drop_coll_data
  case ('glass')
     call glass_grid
     call glass_data
  case ('density ratio')
     call density_ratio_grid
     call density_ratio_data
 ! case ('concentric cartesian jet')
 !    call jetcart2_grid
 !    call jetcart2_data
  case ('uniform')
     call uniform_grid
     call uniform_data
  case ('simple jet')
     call simplejet_grid
     call simplejet_data
!!$  case ('sydney inflow')
!!$     call sydneyinflow_grid
!!$     call sydneyinflow_data
!!$  case ('sydney inflow fine')
!!$     call sydneyinflowfine_grid
!!$     call sydneyinflowfine_data
!!$  case ('sydney fine')
!!$     call sydneyfine_grid
!!$     call sydneyfine_data
!!$  case ('sydney')
!!$     call sydney_grid
!!$     call sydney_data
  case ('bluffbody_SD')
     call bluffbody_SD_grid
     call bluffbody_SD_data
  case default
     print*, "Unknown simulation", simulation
     print*, "Assuming a 2D Gambit mesh file is present"
     call gambit_grid
     call gambit_data
     call gambit_optdata
  end select
  
  ! Write the files
  call param_write_config(simulation)
  call param_write_data
  call param_write_optdata
  call param_write_inflow
  call param_write_chemtable

end program init_flow
     
! ============================== !
! Read the command line argument !
! ============================== !
subroutine commandline_args(input_name)
  implicit none
  integer :: n_input
  !integer, external :: iargc
  integer, intrinsic :: iargc
  character(len=*) :: input_name
  !external :: getarg
  intrinsic :: getarg

  n_input = iargc()
  if (n_input.GT.1) stop "Only supports one command line argument (input file name)"
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
