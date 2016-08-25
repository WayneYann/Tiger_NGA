module finitechem
  use combustion
  use string
  implicit none
  
  ! Clipping values for temperature
  real(WP), parameter :: Tmin =  50.0_WP
  real(WP), parameter :: Tmax = 5000.0_WP
  ! Properties of the species
  real(WP), dimension(:), pointer :: W_sp, Cp_sp, h_sp ! Cp [J/(mol.K)]
  ! Molecular weight [kg/kmol] average over the grid point
  ! Used for temporary calculations
  real(WP)  :: W_mix
  ! Thermo data
  real(WP), parameter :: R_cst =8314.34_WP  ! [J/(kmol.K)]
  ! Number of reactions, species, non_steady species
  integer :: Nreactions, N_tot, N_nons
  ! First mass fraction scalar index
  integer :: isc_sc

  ! Molar diffusion or not
  logical :: molar

  ! Output variables
  ! Chemical and diffusion source term
  !real(WP), dimension(:,:,:,:), pointer ::  diff_src
  ! Diffusion correction term
  real(WP), dimension(:,:,:,:), pointer ::  diff_corr

  ! CVODE tolerances
  real(WP) :: rtol, atol

  ! Variables used to calculate the viscosity, diffusivity, Lewis no
  real(WP),dimension(:),pointer :: kOverEps, mucoeff
  real(WP),dimension(:),pointer :: Lewis
  ! Used for checks throughout code
  logical :: ispresent

  ! External ignition source
  logical :: ignition

  ! Gas radiation
  logical :: gasrad
  integer :: isc_H2O, isc_CO2, isc_CH4, isc_CO

  ! Burner heat loss
  logical :: heatloss
  
  ! Soot radiation
  logical :: finitechem_sootrad

  ! Suppress chemical source term
  logical :: coldflow

  ! Heat release rate and chemical source terms at every grid point
  real(WP), dimension(:,:,:), pointer :: HR_gp
  real(WP), dimension(:,:,:,:), pointer :: chemsrc_gp
  real(WP), dimension(:,:,:,:), pointer :: diff_gp

contains

  ! ======================================================= !
  ! Compute the average Molecular weight in a cell(kg/kmol) !
  ! ======================================================= !
  subroutine finitechem_W(scalar, Wmix)
    implicit none
    real(WP) :: Wmix
    real(WP), dimension(N_nons) :: scalar, scalar_t
    
    scalar_t = min(max(scalar,0.0_WP),1.0_WP)
    Wmix = sum(scalar_t)/sum(scalar_t/W_sp)

    return
  end subroutine finitechem_W

  ! =================================================== !
  ! Compute the average Cp in the mixture (J kg-1 K-1)  !
  ! =================================================== !
  subroutine finitechem_Cp(sol,temp,Cp_mix)
    implicit none
    real(WP) :: temp, Cp_mix
    real(WP), dimension(N_nons) :: sol

    ! Function in mechanism.f file
    call COMPTHERMODATA(h_sp,Cp_sp,temp)
    Cp_mix = sum(sol*Cp_Sp)

  end subroutine finitechem_Cp
  
  ! ================================================================ !
  ! Computes the chemical source term of the system (called by solver)
  ! ================================================================ !
  subroutine finitechem_compute_rhs(sol,rhs)
    use time_info
    implicit none
    
    real(WP), dimension(N_nons+1) :: sol,rhs
    integer :: n
    real(WP), dimension(Nreactions)    :: K_rxn, omega_rxn, M_rxn
    ! Assign workspace for all species conc
    real(WP), dimension(N_tot)      :: conc
    ! Assign cdot for non steady state species
    real(WP), dimension(N_nons) :: conc_dot
    real(WP) :: tmp, P_gp, rho_gp, Cp_mixture

    rhs = 0.0_WP
    conc = 0.0_WP
    conc_dot = 0.0_WP

    ! Compressible or Low-Mach Constant Volume
    ! -- Solving the internal energy form of the temperature equation
    ! Low-Mach Constant Pressure
    ! -- Solving the enthalpy form of the temperature equation
    
!!$    if (.true.) then!compressible .or. is_constant_volume) then
       ! Evaluation thermodynamic pressure at n+1
       P_gp = 0.5_WP*(P_thermo+P_thermo_old)

       ! Density from sum of rho*Y
       ! -- Generalization for soot
       rho_gp = sum(sol(1:N_nons))

       ! Get the W mol and Cp of the mixture
!!$       call finitechem_W(sol(1:N_nons)/rho_gp,W_mix)
       call finitechem_Cp(sol(1:N_nons)/rho_gp,sol(N_nons+1)/rho_gp,Cp_mixture)

       ! Calculate the concentrations of the non steady state species
       do n=1, N_nons
          conc(n)=sol(n)/W_sp(n)
       end do

       ! Obtain the source terms from mechanism.f
       call PRODRATES(conc_dot,omega_rxn,K_rxn,conc,M_rxn,sol(N_nons+1)/rho_gp,P_gp)

       ! Concentration ->  mass fraction for non steady species
       do n=1, N_nons
          rhs(n)=conc_dot(n)*W_sp(n)
       end do
!!$    else
!!$       ! Get the W and Cp of the mixture
!!$       call finitechem_W(sol(1:N_nons),W_mix)
!!$       call finitechem_Cp(sol(1:N_nons),sol(N_nons+1),Cp_mixture)
!!$       rho_gp = P_thermo*W_mix/(R_cst*sol(N_nons+1))
!!$
!!$       ! Calculation the concentrations of the non-steady-state species
!!$       do n=1, N_nons
!!$          conc(n)=rho_gp*sol(n)/(W_sp(n))
!!$       end do
!!$
!!$       ! Obtain the source terms from mechanism.f
!!$       call PRODRATES(conc_dot,omega_rxn,K_rxn,conc,M_rxn,sol(N_nons+1),P_thermo)
!!$
!!$       ! Concentration ->  mass fraction for non steady species
!!$       do n=1, N_nons
!!$          rhs(n)=conc_dot(n)*W_sp(n)/rho_gp
!!$       end do
!!$    end if

    ! The contribution comes from non steady species
    tmp = 0.0_WP
    do n=1, N_nons
       if (.not.compressible .and. .not.is_constant_volume) then
          tmp = tmp - h_sp(n)*rhs(n)
       else
          tmp = tmp - (h_sp(n)-R_cst/W_sp(n)*sol(N_nons+1)/rho_gp)*rhs(n)
       end if
    end do
    ! Chemical Source term for temperature
    if (.not.compressible .and. .not.is_constant_volume) then
       rhs(N_nons+1) = tmp/Cp_mixture
    else
       rhs(N_nons+1) = tmp/(Cp_mixture-R_cst/W_mix)
    end if

  end subroutine finitechem_compute_rhs



end module finitechem

! ================================= !
! Initialize the One Step Chemistry !
! ================================= !
subroutine finitechem_init
  use finitechem
  use parser

  implicit none
  
  integer  :: i, isc
  logical  :: sp_found
  ! Check first scalar index
  character (len=20), dimension(:), pointer :: name_check
  ! List of variables to be outputted
  character(len=str_medium), dimension(:), pointer :: list, list_tot
  character(len=str_medium) :: name, order_check

  ! Get number of species
  call GETNSPECIES(N_tot)

  ! Get number of non steady species
  call GETNSPECS(N_nons)

  allocate (name_check(N_tot))
  allocate(W_sp(N_tot))
  W_sp = 0.0_WP
  
  ! Read the names to find index of first species from mechanism.f file 
  call GETSPECIESNAMES(name_check)
  ! Read the molecular masses from mechanism.f file
  call GETMOLARMASS(W_sp)
  ! Get number of reactions from mechanism.f file
  call GETNREACTIONS(Nreactions)

  ! Gas radiation
  call parser_read('Gas radiation',gasrad,.false.)

  ! Find first species index by comparing scalar names with first one
  sp_found=.false.
  i = 0
  if (N_nons.gt.0) then
     do while (.not.sp_found)
        i = i+1
        if (i.gt.nscalar) then
           call die('Species ' // trim(name_check(1)) // 'not initialized. Run init flow again')
        end if
        order_check=name_check(1)
        if (SC_name(i).eq.trim(order_check(1:str_short))) then
           isc_sc = i
           sp_found=.true.
        end if
     end do
  end if

  ! Check the species order to make sure 
  ! it is the same as the initialization
  if (N_nons.gt.0) then
     do i=isc_sc,isc_sc+N_nons-1
        order_check=name_check(i-isc_sc+1)
        if (trim(SC_name(i)).ne.trim(order_check(1:str_short))) then
           call die('Species ' // trim(name_check(i-isc_sc+1)) &
                // ' not found in initial file')
        end if
        if (.true.) then ! get index for radiation coefficients ! if (gasrad)
           if (trim(SC_name(i)).eq.'H2O') isc_H2O = i
           if (trim(SC_name(i)).eq.'CO2') isc_CO2 = i
           if (trim(SC_name(i)).eq.'CH4') isc_CH4 = i
           if (trim(SC_name(i)).eq.'CO') isc_CO = i
        end if
     end do
  end if
  
  ! Molar diffusion or not
  call parser_read('Molar diffusion',molar,.true.)

  ! External heat source for ignition
  call parser_read('External ignition',ignition,.false.)

  ! Heat loss to the burner
  call parser_read('Burner heat loss',heatloss,.false.)

  ! Soot radiation
  call parser_read('Soot radiation',finitechem_sootrad,.false.)
  
  ! CVODE Tolerances
  call parser_read('Finitechem rel. tol.',rtol,1.0e-6_WP)
  call parser_read('Finitechem abs. tol.',atol,1.0e-12_WP)

  ! Allocate Arrays
  ! First derivative of the scalar field, source terms for diff, chem
  !allocate(diff_src(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,N_nons+1))
  allocate(diff_corr(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,N_nons+1))
  allocate(HR_gp(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(chemsrc_gp(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,3))  ! For H2O, CO, and CO2
  allocate(diff_gp(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,3))  ! For H2O, CO, and CO2
  HR_gp = 3.0_WP
  chemsrc_gp = 1.0e-60_WP

  ! Properties of species
  ! Cp is calculated for the non steady species
  allocate(Cp_sp(N_nons))
  ! h is calculated for only the non steady species 
  allocate(h_sp(N_nons))
  ! kOverOmega is calculated for all the species
  allocate(kOverEps(N_tot))
  ! mucoeff is calculated for all the species
  allocate(mucoeff(N_tot))
  ! Get mucoeff from mechanism.f file
  call GETMUCOEFF(mucoeff)
  ! Get kOverEps from mechanism.f file
  call GETKOVEREPS(kOverEps)
  ! Allocate Lewis numbers
  allocate(Lewis(N_nons))
  ! Get Lewis numbers
  Lewis = 1.0_WP
  do isc=isc_sc,isc_sc+N_nons-1
     call parser_is_defined('Lewis '// SC_name(isc),ispresent)!
     if (ispresent) then
        call parser_read('Lewis '// SC_name(isc),Lewis(isc-isc_sc+1))
     end if
  end do

!!$! SD
!!$tmp1 = 0.0_WP
!!$
!!$do i=imin_,imax_
!!$   do j=jmin_,jmax_
!!$      do k=kmin_,kmax_
!!$         tmp1(i,j,k) = sum(SC(i,j,k,isc_sc:isc_sc+N_nons-1))
!!$         if (mask(i,j).eq.1) tmp1(i,j,k) = 1.0_WP
!!$      end do
!!$   end do
!!$end do
!!$
!!$print*,irank,maxval(tmp1(imin_:imax_,jmin_:jmax_,kmin_:kmax_)),minval(tmp1(imin_:imax_,jmin_:jmax_,kmin_:kmax_))


  return
end subroutine finitechem_init

! ========================================================== !
! Compute the new density from the mass fraction             !
! ========================================================== !
subroutine finitechem_density
  use finitechem
  use masks
  implicit none
  
  integer :: i,j,k
  real(WP):: T_density

  ! Compute the new density fron the equation of state
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (mask(i,j).ne.1) then
              call finitechem_W(SC(i,j,k,isc_sc:isc_sc-1+N_nons), W_mix)
              RHO(i,j,k) = P_thermo*W_mix/(R_cst*SC(i,j,k,isc_T))
           end if
        end do
     end do
  end do

  return
end subroutine finitechem_density


! ========================================================= !
! Compute the viscosity by fitting using FlameMaster method !
! ========================================================= !
subroutine finitechem_viscosity
  use finitechem
  use masks
  implicit none

  integer  :: i,j,k,nsc,j_sc,k_sc
  real(WP) :: temp,temp1,sum_1,sum_2
  real(WP),dimension(N_nons) :: eta
  real(WP),dimension(N_nons,N_nons) :: phi
  real(WP) :: omegamu

  do k = kmino_,kmaxo_
     do j = jmino_,jmaxo_
        do i = imino_,imaxo_
           if (mask(i,j).ne.1) then
              temp=min(max(SC(i,j,k,isc_T ),Tmin),Tmax)
              do nsc=1,N_nons
                 call finitechem_omegamu(temp*kOverEps(nsc),omegamu)
                 eta(nsc)=mucoeff(nsc)*sqrt(temp)/omegamu
              end do
              call finitechem_W(SC(i,j,k,isc_sc:isc_sc+N_nons-1), W_mix)
              ! Use Wilke's method (same as Chemkin Manual 
              ! but obtained from Wilke's paper)
              do k_sc=1,N_nons
                 do j_sc=1,N_nons
                    phi(k_sc,j_sc)=(1.0_WP/sqrt(8.0_WP))&
                         *((1.0_WP+(W_sp(k_sc)/W_sp(j_sc)))**-0.5_WP)
                    temp1=((eta(k_sc)/eta(j_sc))**0.5_WP)&
                         *((W_sp(j_sc)/W_sp(k_sc))**0.25_WP)
                    phi(k_sc,j_sc)=phi(k_sc,j_sc)*((1.0_WP+temp1)**2)
                 end do
              end do
              sum_2=0.0_WP
              do k_sc=1,N_nons
                 sum_1=0.0_WP
                 do j_sc=1,N_nons
                    if (j_sc.ne.k_sc) then
                       sum_1=sum_1+(SC(i,j,k,j_sc+isc_sc-1)*W_mix*&
                            phi(k_sc,j_sc)/W_sp(j_sc))
                    end if
                 end do
                 sum_1=sum_1+SC(i,j,k,k_sc+isc_sc-1)*W_mix/W_sp(k_sc)
                 if (SC(i,j,k,k_sc+isc_sc-1).gt.0.0_WP) then
                    sum_2=sum_2+(SC(i,j,k,k_sc+isc_sc-1)*W_mix*&
                         eta(k_sc)/(W_sp(k_sc)*sum_1))
                 end if
              end do
              VISC(i,j,k)=sum_2
           end if
        end do
     end do
  end do
  call boundary_update_border(VISC,'+','ym')

  return
end subroutine finitechem_viscosity

!=====================================================!
! Compute the viscosity by fitting the transport data !
! obtained using FlameMaster!
!=====================================================!
subroutine finitechem_diffusivity
  use finitechem
  use masks
  implicit none

  integer :: i,j,k,nsc
  real(WP),dimension(N_tot) :: cond, eta
  real(WP) :: temp,sum_1,sum_2,lambda, Cp_mix, omegamu

  do k = kmino_,kmaxo_
     do j = jmino_,jmaxo_
        do i = imino_,imaxo_
           if (mask(i,j).ne.1) then
              call finitechem_W(SC(i,j,k,isc_sc:isc_sc-1+N_nons), W_mix)
              temp=min(max(SC(i,j,k,isc_T ),Tmin),Tmax)
              do nsc=1,N_nons
                 call finitechem_omegamu(temp*kOverEps(nsc),omegamu)
                 eta(nsc)=mucoeff(nsc)*sqrt(temp)/omegamu
                 call COMPTHERMODATA(h_sp,Cp_sp,temp)
                 cond(nsc)=eta(nsc)*(Cp_sp(nsc)+1.25_WP*R_cst/W_sp(nsc))
              end do
              sum_1=0.0_WP
              sum_2=0.0_WP
              do nsc=1,N_nons
                 sum_1=sum_1+(SC(i,j,k,nsc+isc_sc-1)*W_mix/(cond(nsc)*W_sp(nsc)))
                 sum_2=sum_2+(SC(i,j,k,nsc+isc_sc-1)*cond(nsc)*W_mix/W_sp(nsc))
              end do
              ! Average thermal conductivity of the mixture
              lambda=0.5_WP*(sum_2+(1.0_WP/sum_1))
              ! Average Cp based on scalar field
              call finitechem_Cp(SC(i,j,k,isc_sc:isc_sc+N_nons-1),temp,Cp_mix)
              ! Since scalar solver has rho*scalar,coeff is not divided by RHO
              do nsc=isc_sc,isc_sc+N_nons-1
                 DIFF(i,j,k,nsc)=(lambda/Cp_mix)/Lewis(nsc-isc_sc+1)
              end do
              ! Set temperature diffusivity equal to the thermal diffusivity
              DIFF(i,j,k,isc_T)=lambda/Cp_mix
              ! Set mixture fraction diffusivity equal to the thermal diffusivity
              if (isc_ZMIX.ne.0) DIFF(i,j,k,isc_ZMIX)=lambda/Cp_mix
              ! Set enthalpy to the thermal diffusivity
              !! WATCH THIS...
              !DIFF(i,j,k,isc_enthalpy)=lambda/Cp_mix
           end if
        end do
     end do
  end do
  call boundary_update_border(DIFF,'+','ym')
  
 return
end subroutine finitechem_diffusivity

! ==================================================== !
! Compute the omegamu needed for transport calculation !
! Obtained from FlameMaster - Tspecies.c               !
! ==================================================== !
subroutine finitechem_omegamu(temp,omegamu)
  use finitechem
  implicit none
  real(WP) ::  m1 = 3.3530622607_WP
  real(WP) ::  m2 = 2.53272006_WP
  real(WP) ::  m3 = 2.9024238575_WP
  real(WP) ::  m4 = 0.11186138893_WP
  real(WP) ::  m5 = 0.8662326188_WP  ! = -0.1337673812 + 1.0
  real(WP) ::  m6 = 1.3913958626_WP
  real(WP) ::  m7 = 3.158490576_WP
  real(WP) ::  m8 = 0.18973411754_WP
  real(WP) ::  m9 = 0.00018682962894_WP
  real(WP) ::  num, den;
  real(WP) ::  temp, omegamu
  
  num = m1 + temp*(m2 + temp*(m3 + temp*m4))
  den = m5 + temp*(m6 + temp*(m7 + temp*(m8 + temp*m9)))
  omegamu = num/den

  return
end subroutine finitechem_omegamu


! =========================================== !
! Non-chemistry source terms for finitechem   !
! -- Chemical source terms computed elsewhere !
! =========================================== !
subroutine finitechem_source_transport(src)
  use finitechem
  use time_info
  use masks
  implicit none

  integer :: i,j,k, n
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: src  
  real(WP) :: Cp_mix
  
  ! Mass transfer terms
  call finitechem_diffusion

  ! Molar diffusion and correction for species equations
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
              do n = isc_sc, isc_sc+N_nons-1
                 src(i,j,k,n) = src(i,j,k,n) + diff_corr(i,j,k,n)
              end do
           end if
        end do
     end do
  end do
  
  ! Temperature terms: heat capacity gradient and enthalpy flux
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
                 src(i,j,k,isc_T) = src(i,j,k,isc_T) + diff_corr(i,j,k,N_nons+1)
              end if
        end do
     end do
  end do
  
  return
end subroutine finitechem_source_transport


! ======================================= !
! Transport "source" terms for finitechem !
! ======================================= !
subroutine finitechem_diffusion
  use finitechem
  use metric_generic
  use memory
  use time_info
  use masks
  implicit none
  ! Every grid point molecular weight
  real(WP), dimension(:,:,:), pointer :: W_gp, Cp_mix

  ! Correction velocity in 3 dimensions
  real(WP), dimension(:,:,:), pointer :: corr_vel_x
  real(WP), dimension(:,:,:), pointer :: corr_vel_y
  real(WP), dimension(:,:,:), pointer :: corr_vel_z
  integer :: i,j,k, isc

  ! Gas radiation
  real(WP) :: Qrad,alphaH2O,alphaCO2,alphaCH4,alphaCO
  real(WP) :: c0,c1,c2,c3,c4,c5
  real(WP) :: tm1,T2
  real(WP), parameter :: Tb = 298.0_WP
  real(WP), parameter :: radconst = 4.0_WP * 5.67051e-8_WP * R_cst / 101320.0_WP ! [W/(m^2-K^4)]|[atm to Pa]

  ! Burner heat loss
  real(WP) :: Qloss
  
  diff_corr = 0.0_WP
  FX = 0.0_WP
  FY = 0.0_WP
  FZ = 0.0_WP
  Cp_mix => tmp9
  W_gp => tmp10
  corr_vel_x => tmp11
  corr_vel_y => tmp12
  corr_vel_z => tmp13
  corr_vel_x = 0.0_WP
  corr_vel_y = 0.0_WP
  corr_vel_z = 0.0_WP
  diff_gp = 0.0_WP

  ! 1.  Calculate the diffusion flux based on molar diffusion (DIFF * grad(Y*W_gp/W_sp) * W_sp/W_gp)
  ! 2.  Recall that the generic scalar transport equation already has Fickian diffusion
  !     2.1.  Get the net molar diffusion flux(sum over species)
  !     2.2.  Get the molar diffusion correction
  !     2.3.  Ensure zero net molar diffusion flux
  ! 3.  Use the corrected molar diffusion flux to compute enthalpy flux
  ! 4.  cp/cv corrector which combines with the generic diffusion term to get the Fourier term

  ! Compute the mixture molar mass
  W_gp = 1.0_WP
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_ 
           if (mask(i,j).ne.1) then
              call finitechem_W(SC(i,j,k,isc_sc:isc_sc-1+N_nons), W_gp(i,j,k))
           end if
        end do
     end do
  end do

  ! Species diffusion correction velocities
  do isc=isc_sc,isc_sc+N_nons-1
     ! Fickian diffusion
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              if (mask(i,j).ne.1) then
                 corr_vel_x(i,j,k) = corr_vel_x(i,j,k) + &
                      sum(interp_sc_x(i,j,:) * DIFF(i-st2:i+st1,j,k,isc)) * &
                      sum(grad_x(i,j,:) * SC(i-st2:i+st1,j,k,isc))
                 corr_vel_y(i,j,k) = corr_vel_y(i,j,k) + &
                      sum(interp_sc_y(i,j,:) * DIFF(i,j-st2:j+st1,k,isc)) * &
                      sum(grad_y(i,j,:) * SC(i,j-st2:j+st1,k,isc))
                 corr_vel_z(i,j,k) = corr_vel_z(i,j,k) + &
                      sum(interp_sc_z(i,j,:) * DIFF(i,j,k-st2:k+st1,isc)) * &
                      sum(grad_z(i,j,:) * SC(i,j,k-st2:k+st1,isc))
              end if
           end do
        end do
     end do

     ! Molar diffusion
     if (molar) then
        do k=kmin_-st1,kmax_+st2
           do j=jmin_-st1,jmax_+st2
              do i=imin_-st1,imax_+st2
                 if (mask(i,j).ne.1) then
                    corr_vel_x(i,j,k) = corr_vel_x(i,j,k) + &
                         sum(interp_sc_x(i,j,:) * DIFF(i-st2:i+st1,j,k,isc)) * &
                         sum(interp_sc_x(i,j,:) * SC(i-st2:i+st1,j,k,isc) / W_gp(i-st2:i+st1,j,k)) * &
                         sum(grad_x(i,j,:) * W_gp(i-st2:i+st1,j,k))
                    corr_vel_y(i,j,k) = corr_vel_y(i,j,k) + &
                         sum(interp_sc_y(i,j,:) * DIFF(i,j-st2:j+st1,k,isc)) * &
                         sum(interp_sc_y(i,j,:) * SC(i,j-st2:j+st1,k,isc) / W_gp(i,j-st2:j+st1,k)) * &
                         sum(grad_y(i,j,:) * W_gp(i,j-st2:j+st1,k))
                    corr_vel_z(i,j,k) = corr_vel_z(i,j,k) + &
                         sum(interp_sc_z(i,j,:) * DIFF(i,j,k-st2:k-st1,isc)) * &
                         sum(interp_sc_z(i,j,:) * SC(i,j,k-st2:k-st1,isc) / W_gp(i,j,k-st2:k+st1)) * &
                         sum(grad_z(i,j,:) * W_gp(i,j,k-st2:k+st1))
                 end if
              end do
           end do
        end do
     end if
  end do

  ! Diffusion fluxes
  do isc=isc_sc,isc_sc+N_nons-1
     ! Correction velocity
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              if (mask(i,j).ne.1) then             
                 FX(i,j,k) = -corr_vel_x(i,j,k) * &
                      sum(interp_sc_x(i,j,:)*SC(i-st2:j+st1,j,k,isc))
                 FY(i,j,k) = -corr_vel_y(i,j,k) * &
                      sum(interp_sc_y(i,j,:)*SC(i,j-st2:j+st1,k,isc))
                 FZ(i,j,k) = -corr_vel_z(i,j,k) * &
                      sum(interp_sc_z(i,j,:)*SC(i,j,k-st2:k+st1,isc))
              end if
           end do
        end do
     end do

     ! Molar diffusion flux
     if (molar) then
        do k=kmin_-st1,kmax_+st2
           do j=jmin_-st1,jmax_+st2
              do i=imin_-st1,imax_+st2
                 if (mask(i,j).ne.1) then
                    FX(i,j,k) = FX(i,j,k) + sum(interp_sc_x(i,j,:) * DIFF(i-st2:i+st1,j,k,isc)) * &
                         sum(interp_sc_x(i,j,:) * SC(i-st2:i+st1,j,k,isc) / W_gp(i-st2:i+st1,j,k)) * &
                         sum(grad_x(i,j,:) * W_gp(i-st2:i+st1,j,k))
                    FY(i,j,k) = FY(i,j,k) + sum(interp_sc_y(i,j,:) * DIFF(i,j-st2:j+st1,k,isc)) * &
                         sum(interp_sc_y(i,j,:) * SC(i,j-st2:j+st1,k,isc) / W_gp(i,j-st2:j+st1,k)) * &
                         sum(grad_y(i,j,:) * W_gp(i,j-st2:j+st1,k))
                    FZ(i,j,k) = FZ(i,j,k) + sum(interp_sc_z(i,j,:) * DIFF(i,j,k-st2:k+st1,isc)) * &
                         sum(interp_sc_z(i,j,:) * SC(i,j,k-st2:k+st1,isc) / W_gp(i,j,k-st2:k+st1)) * &
                         sum(grad_z(i,j,:) * W_gp(i,j,k-st2:k+st1))
                 end if
              end do
           end do
        end do
     end if

     ! Species source term from molar diffusion flux and correction velocity
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).ne.1) then
                 diff_corr(i,j,k,isc-isc_sc+1) = &
                      sum(div_u(i,j,:) * FX(i-st1:i+st2,j,k)) + &
                      sum(div_v(i,j,:) * FY(i,j-st1:j+st2,k)) + &
                      sum(div_w(i,j,:) * FZ(i,j,k-st1:k+st2))
              end if
           end do
        end do
     end do

     ! Total flux for enthalpy flux
     do k=kmin_-st1,kmax_+st2
        do j=jmin_-st1,jmax_+st2
           do i=imin_-st1,imax_+st2
              if (mask(i,j).ne.1) then
                 FX(i,j,k) = FX(i,j,k) + sum(interp_sc_x(i,j,:) * DIFF(i-st2:i+st1,j,k,isc)) * &
                      sum(grad_x(i,j,:) * SC(i-st2:i+st1,j,k,isc))
                 FY(i,j,k) = FY(i,j,k) + sum(interp_sc_y(i,j,:) * DIFF(i,j-st2:j+st1,k,isc)) * &
                      sum(grad_y(i,j,:) * SC(i,j-st2:j+st1,k,isc))
                 FZ(i,j,k) = FZ(i,j,k) + sum(interp_sc_z(i,j,:) * DIFF(i,j,k-st2:k-st1,isc)) * &
                      sum(grad_z(i,j,:) * SC(i,j,k-st2:k+st1,isc))
              end if
           end do
        end do
     end do

     ! Enthalpy flux for temperature
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).ne.1) then
                 call COMPTHERMODATA(h_sp,Cp_sp,SC(i,j,k,isc_T))
                 ! X
                 diff_corr(i,j,k,N_nons+1) =  diff_corr(i,j,k,N_nons+1) + Cp_sp(isc) * &
                      sum(interp_u_xm(i,j,:) * FX(i-st1:i+st2,j,k)) * &
                      sum(grad_xm(i,j,:) * SC(i-stp:i+stp,j,k,isc_T))

                 ! Y
                 diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + Cp_sp(isc) * &
                      sum(interp_v_ym(i,j,:) * FY(i,j-st1:j+st2,k)) * &
                      sum(grad_ym(i,j,:) * SC(i,j-stp:j+stp,k,isc_T))
                 ! Z
                 diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + Cp_sp(isc) * &
                      sum(interp_w_zm(i,j,:) * FZ(i,j,k-st1:k+st2)) * &
                      sum(grad_zm(i,j,:) * SC(i,j,k-stp:k+stp,isc_T))
              end if
           end do
        end do
     end do
  end do

  ! Include gas radiation
  if (gasrad) then
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              if (mask(i,j).ne.1) then
                 tm1 = 1000.0_WP / SC(i,j,k,isc_T)
                 T2 = SC(i,j,k,isc_T)**2.0_WP
                 ! first H2O
                 c0 = -0.23093_WP
                 c1 = -1.12390_WP
                 c2 =  9.41530_WP
                 c3 = -2.99880_WP
                 c4 =  0.51382_WP
                 c5 = -1.86840E-05_WP
                 alphaH2O = c0 + Tm1 * (c1 + Tm1 * (c2 + Tm1 * (c3 + Tm1 * (c4 + Tm1 * c5)))) ! (atm m)^-1
                 ! then CO2
                 c0 =  18.741_WP
                 c1 = -121.310_WP
                 c2 =  273.500_WP
                 c3 = -194.050_WP
                 c4 =  56.310_WP
                 c5 = -5.8169_WP
                 alphaCO2 = c0 + Tm1 * (c1 + Tm1 * (c2 + Tm1 * (c3 + Tm1 * (c4 + Tm1 * c5)))) ! (atm m)^-1
                 ! then CH4
                 alphaCH4 = 6.6334_WP - 0.0035686_WP * SC(i,j,k,isc_T) + 1.6682e-08_WP * T2 &
                      + 2.5611e-10_WP * T2 * SC(i,j,k,isc_T) - 2.6558e-14_WP * T2 * T2 !(atm m)^-1
                 ! and CO
                 if (SC(i,j,k,isc_T).le.750.0_WP) then
                    c0 = 4.7869_WP
                    c1 = -0.06953_WP
                    c2 = 2.95775e-4_WP
                    c3 = -4.25732e-7_WP
                    c4 = 2.02894e-10_WP
                 else
                    c0 = 10.09_WP
                    c1 = -0.01183_WP
                    c2 = 4.7753e-6_WP
                    c3 = -5.87209e-10_WP
                    c4 = -2.5334e-14_WP
                 end if
                 alphaCO =  c0 + SC(i,j,k,isc_T) * (c1 + SC(i,j,k,isc_T) * &
                      (c2 + SC(i,j,k,isc_T) * (c3 + SC(i,j,k,isc_T) * c4))) ! (atm m)^-1
                 ! Gas radiation term
                 Qrad = -radconst * RHO(i,j,k) * SC(i,j,k,isc_T) * &
                      (SC(i,j,k,isc_T)**4.0_WP - Tb**4.0_WP) * &
                      (alphaCO2 * SC(i,j,k,isc_CO2) / W_sp(isc_CO2) + &
                      alphaH2O * SC(i,j,k,isc_H2O) / W_sp(isc_H2O) + &
                      alphaCH4 * SC(i,j,k,isc_CH4) / W_sp(isc_CH4) + &
                      alphaCO * SC(i,j,k,isc_CO) / W_sp(isc_CO))
                 diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + Qrad
              end if
           end do
        end do
     end do
  end if

  ! Burner heat loss
  ! According to the mesh file, find the first point downstream of the wall
  if (heatloss) then
     do j=48,55
        Qloss = -max(SC(5,j+2,3,isc_T)-298.0_WP, 0.0_WP) * RHO(5,j+2,3) / 5e-6_WP
        diff_corr(5,j+2,3,N_nons+1) = diff_corr(5,j+2,3,N_nons+1) + Qloss
     end do
  end if
  
  ! Compute specific heat
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_ 
           if (mask(i,j).ne.1) then
              call finitechem_Cp(SC(i,j,k,isc_sc:isc_sc-1+N_nons),SC(i,j,k,isc_T),Cp_mix(i,j,k))
              ! cv for compressible or constant volume
              if (compressible .or. is_constant_volume) Cp_mix(i,j,k) = Cp_mix(i,j,k) - R_cst / W_gp(i,j,k)
           end if
        end do
     end do
  end do

  ! cp/cv correction: lambda/(Cp^2) * d(Cp)/d(x) * d(T)/d(x)
  ! cp_mix is actually cv for compressible or constant volume
  ! Note: Dividing enthalpy flux from above by cp (or cv)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
              ! X
              diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + &
                   DIFF(i,j,k,isc_T) * sum(grad_xm(i,j,:) * Cp_mix(i-stp:i+stp,j,k)) * &
                   sum(grad_xm(i,j,:) * SC(i-stp:i+stp,j,k,isc_T))
              ! Y
              diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + &
                   DIFF(i,j,k,isc_T) * sum(grad_ym(i,j,:) * Cp_mix(i,j-stp:j+stp,k)) * &
                   sum(grad_ym(i,j,:) * SC(i,j-stp:j+stp,k,isc_T))
              ! Z
              diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + &
                   DIFF(i,j,k,isc_T) * sum(grad_zm(i,j,:) * Cp_mix(i,j,k-stp:k+stp)) * &
                   sum(grad_zm(i,j,:) * SC(i,j,k-stp:k+stp,isc_T))

              diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) / Cp_mix(i,j,k)
           end if
        end do
     end do
  end do

  ! SD
  ! External ignition source
  if (ignition) then
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ! For coflow flames
              if (mask(i,j).ne.1 .and. x(i).ge.0.00 .and. x(i).le.0.015 .and. y(j).ge.0.0015 .and. y(j).le.0.0055) then
                 diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + max(2000.0_WP * RHO(i,j,k) - SC(i,j,k,isc_T) * RHO(i,j,k), 0.0_WP) / 5e-3_WP
                 diff_corr(i,j,k,5) = diff_corr(i,j,k,5) + max(1.0e-2_WP * RHO(i,j,k) - SC(i,j,k,isc_sc+4) * RHO(i,j,k), 0.0_WP) / 5e-3_WP !OH
              end if
              ! For 1D planar premixed flames
!!$              if (mask(i,j).ne.1 .and. i.eq.101) then
!!$                 diff_corr(i,j,k,N_nons+1) = diff_corr(i,j,k,N_nons+1) + max(1500.0_WP * RHO(i,j,k) - SC(i,j,k,isc_T) * RHO(i,j,k), 0.0_WP) / 2e-6_WP
!!$              end if
           end do
        end do
     end do
  end if
  ! SD end

  ! SD To dump the diffusion terms (including Fickian diffusion)
  corr_vel_x = 0.0_WP
  corr_vel_y = 0.0_WP
  corr_vel_z = 0.0_WP
  
  do isc=isc_sc,isc_sc+N_nons-1
     if (isc.eq.isc_H2O .or. isc.eq.isc_CO .or. isc.eq.isc_CO2) then
        ! Get Fickian diffusion
        do k=kmin_-st1,kmax_+st2
           do j=jmin_-st1,jmax_+st2
              do i=imin_-st1,imax_+st2
                 if (mask(i,j).ne.1) then
                    corr_vel_x(i,j,k) = sum(interp_sc_x(i,j,:) * DIFF(i-st2:i+st1,j,k,isc)) * &
                         sum(grad_x(i,j,:) * SC(i-st2:i+st1,j,k,isc))
                    corr_vel_y(i,j,k) = sum(interp_sc_y(i,j,:) * DIFF(i,j-st2:j+st1,k,isc)) * &
                      sum(grad_y(i,j,:) * SC(i,j-st2:j+st1,k,isc))
                    corr_vel_z(i,j,k) = sum(interp_sc_z(i,j,:) * DIFF(i,j,k-st2:k+st1,isc)) * &
                         sum(grad_z(i,j,:) * SC(i,j,k-st2:k+st1,isc))
                 end if
              end do
           end do
        end do
        
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 if (mask(i,j).ne.1) then
                    if (isc.eq.isc_H2O) then
                       diff_gp(i,j,k,1) = &
                            sum(div_u(i,j,:) * corr_vel_x(i-st1:i+st2,j,k)) + &
                            sum(div_v(i,j,:) * corr_vel_y(i,j-st1:j+st2,k)) + &
                            sum(div_w(i,j,:) * corr_vel_z(i,j,k-st1:k+st2))
                       diff_gp(i,j,k,1) = diff_gp(i,j,k,1) + diff_corr(i,j,k,isc_H2O)
                    else if (isc.eq.isc_CO) then
                       diff_gp(i,j,k,2) = &
                            sum(div_u(i,j,:) * corr_vel_x(i-st1:i+st2,j,k)) + &
                            sum(div_v(i,j,:) * corr_vel_y(i,j-st1:j+st2,k)) + &
                            sum(div_w(i,j,:) * corr_vel_z(i,j,k-st1:k+st2))
                       diff_gp(i,j,k,2) = diff_gp(i,j,k,2) + diff_corr(i,j,k,isc_CO)
                    else
                       diff_gp(i,j,k,3) = &
                            sum(div_u(i,j,:) * corr_vel_x(i-st1:i+st2,j,k)) + &
                            sum(div_v(i,j,:) * corr_vel_y(i,j-st1:j+st2,k)) + &
                            sum(div_w(i,j,:) * corr_vel_z(i,j,k-st1:k+st2))
                       diff_gp(i,j,k,3) = diff_gp(i,j,k,3) + diff_corr(i,j,k,isc_CO2)
                    end if
                 end if
              end do
           end do
        end do
     end if
  end do
  ! SD


  ! Do not divide by RHO since the source term is needed that way
  diff_corr=diff_corr*dt_

  return
end subroutine finitechem_diffusion


! ======================================================= !
! Compute pressure-divergence source term for temperature !
! ======================================================= !
subroutine finitechem_pressuredivergence(src_,SC_)
  use finitechem
  use metric_velocity_conv
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar), intent(inout) :: src_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar), intent(inout) :: SC_
  integer :: i,j,k
  real(WP) :: cp_, W_

  if (compressible) call die('This is not implemented yet...')

  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           call finitechem_W(SC_(i,j,k,isc_sc:isc_sc-1+N_nons),W_)
           call finitechem_Cp(SC_(i,j,k,isc_sc:isc_sc-1+N_nons),SC_(i,j,k,isc_T),cp_)
           src_(i,j,k,isc_T) = src_(i,j,k,isc_T) - dt_/(cp_-R_cst/W_)* &
                0.5_WP*(P_thermo_old+P_thermo)*( &
                sum(divc_u(i,j,:)*U(i-stc1:i+stc2,j,k)) + &
                sum(divc_v(i,j,:)*V(i,j-stc1:j+stc2,k)) + &
                sum(divc_w(i,j,:)*W(i,j,k-stc1:k+stc2)))
        end do
     end do
  end do

  return
end subroutine finitechem_pressuredivergence


! ======================================= !
! Compute the source terms for combustion !
! ======================================= !
subroutine finitechem_source_chemistry(RHO_,SC_)
  use finitechem
  use memory
  use time_info
  use masks 
  use parallel
  
  implicit none
  
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: RHO_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar), intent(inout) :: SC_  
  real(WP), dimension(N_nons+1) :: solold,sol
  integer  :: i,j,k,n
  integer(8), dimension(21) :: iout
  real(WP), dimension(6) :: rout
  integer(8) :: ipar
  real(WP) :: rpar, t_, rho_gp
  integer :: cverr
  logical :: reinit
  real(WP) :: Cp_mix

  call parser_read('Cold flow',coldflow,.false.)
  if(coldflow) return

  ! Initialize NVECTOR
  call FNVINITS(1, int(N_nons+1,8), cverr)

  ! Initialization on first entry
  reinit = .false.

  ! Use CVODE to perform the integration
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
              ! Get initial field
              if (.false.) then!.not.compressible .and. .not.is_constant_volume) then
                 sol(1:N_nons) = max(SC_(i,j,k,isc_sc:isc_sc-1+N_nons),0.0_WP)
                 sol(N_nons+1) = SC_(i,j,k,isc_T)
              else
                 sol(1:N_nons) = max(RHO_(i,j,k)*SC_(i,j,k,isc_sc:isc_sc-1+N_nons),0.0_WP)
                 sol(N_nons+1) = RHO_(i,j,k)*SC_(i,j,k,isc_T)
              end if

              ! Store the old solution
              solold  = sol
              
              if (.not.reinit) then
                 ! Initialize CVODE
                 ! -- Third argument (2): Stiff BDF Integration
                 ! -- Fourth argument (2): Nonlinear Newton solver
                 ! -- Fifth argument (1): Scalar absolute tolerance
                 call FCVMALLOC(0.0_WP, sol, 2, 2, 1, rtol, atol, iout, rout, ipar, rpar, cverr)
              else
                 ! Re-initialize CVODE
                 ! -- Third argument (1): Scalar absolute tolerance
                 call FCVREINIT(0.0_WP, sol, 1, rtol, atol, cverr)
              end if

              ! Set stop time
              call FCVSETRIN("STOP_TIME",dt,cverr)

              ! Set maximum steps
              call FCVSETIIN("MAX_NSTEPS",10000,cverr)

              ! Initialize the dense, direct linear solver
              if (.not.reinit) call FCVDENSE(int(N_nons+1,8), cverr)

              ! Integrate with CVODE
              ! -- Fourth argument (3): Integrate to dt and stop
              ! -- Note: Solving for rho*Y and rho*T

              call FCVODE(dt,t_,sol,1,cverr)

              ! Directly update scalar values with solution from CVODE
              if (.false.) then!.not.compressible .and. .not.is_constant_volume) then
                 SC_(i,j,k,isc_sc:isc_sc-1+N_nons) = sol(1:N_nons)
                 SC_(i,j,k,isc_T)                  = sol(N_nons+1)
              else
                 SC_(i,j,k,isc_sc:isc_sc-1+N_nons) = sol(1:N_nons)/RHO_(i,j,k)
                 SC_(i,j,k,isc_T)                  = sol(N_nons+1)/RHO_(i,j,k)
              end if

              ! SD
              ! Get the heat release
              call finitechem_Cp(SC_(i,j,k,isc_sc:isc_sc-1+N_nons),SC_(i,j,k,isc_T),Cp_mix)
              if (.false.) then!.not.compressible .and. .not.is_constant_volume) then
                 HR_gp(i,j,k) = RHO_(i,j,k) * Cp_mix * (sol(N_nons+1) - solold(N_nons+1)) / dt
              else
                 HR_gp(i,j,k) = Cp_mix * (sol(N_nons+1) - solold(N_nons+1)) / dt
                 ! H2O
                 chemsrc_gp(i,j,k,1) = (sol(isc_H2O-isc_sc+1) - solold(isc_H2O-isc_sc+1)) / dt
                 ! CO
                 chemsrc_gp(i,j,k,2) = (sol(isc_CO-isc_sc+1) - solold(isc_CO-isc_sc+1)) / dt
                 ! CO2
                 chemsrc_gp(i,j,k,3) = (sol(isc_CO2-isc_sc+1) - solold(isc_CO2-isc_sc+1)) / dt
              end if
              
              ! Re-initialize on subsequent entries
              if (.not.reinit) reinit = .true.
           end if
        end do
     end do
  end do

  ! Free CVODE memory
  call FCVFREE

  return
end subroutine finitechem_source_chemistry


! ================================================= !
! Computes the residual of the system for one point !
! -- Fixed subroutine name with CVODE               !
! ================================================= !
subroutine FCVFUN(t_, sol, soldot, ipar, rpar, cverr)
  use finitechem
  implicit none
  
  real(WP), intent(in) :: t_
  real(WP), dimension(N_nons+1) :: sol, soldot
  real(WP) :: rpar
  integer(8) :: ipar
  integer :: cverr
  real(WP) :: P_thermo_

  call finitechem_compute_rhs(sol,soldot)

  cverr = 0

end subroutine FCVFUN
