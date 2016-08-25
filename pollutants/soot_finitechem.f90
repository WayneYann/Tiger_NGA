module soot_finitechem
  use soot
  implicit none

  ! Number of aromatic species
  integer, parameter :: nPAH = 8

  ! Size of naphthalene dimers
  real(WP), parameter :: nbrC = 20.0_WP
  real(WP), parameter :: nbrH = 16.0_WP

  ! Soot relevant constants: Also in soot_hmom
  real(WP), parameter :: SootDensity = 1800.0_WP    ! kg/m^3
  real(WP), parameter :: MolarMassSoot = 12.0e-3_WP ! kg/mol
  real(WP), parameter :: Avogadro = 6.022e23_WP     ! 1/mol
  real(WP), parameter :: Df = 1.8_WP                ! 1
  real(WP), parameter :: chisoot = 1.0e19_WP        ! 1/m^2

  ! PAH indices
  integer, dimension(nPAH) :: i_PAH

  ! PAH carbon number
  integer, dimension(nPAH) :: size_PAH

  ! PAH collision efficiencies
  real(WP), dimension(nPAH) :: gamma_PAH

  ! PAH production rates
  real(WP), dimension(nPAH) :: pahprodrate

  ! Surf/ox indices
  integer :: i_OH,i_H2O,i_H,i_H2,i_C2H2,i_O2,i_CO

  ! Surface chemistry rate parameters
  real(WP), parameter :: A1f_s = 6.72e1_WP / 1e6_WP
  real(WP), parameter :: n1f_s = 3.33_WP
  real(WP), parameter :: E1f_s = 6.09_WP * 1e3_WP / R_uni
  real(WP), parameter :: A1b_s = 6.44e-1_WP / 1e6_WP
  real(WP), parameter :: n1b_s = 3.79_WP
  real(WP), parameter :: E1b_s = 27.96_WP * 1e3_WP / R_uni
  real(WP), parameter :: A2f_s = 1.00e8_WP / 1e6_WP
  real(WP), parameter :: n2f_s = 1.80_WP
  real(WP), parameter :: E2f_s = 68.42_WP * 1e3_WP / R_uni
  real(WP), parameter :: A2b_s = 8.68e4_WP / 1e6_WP
  real(WP), parameter :: n2b_s = 2.36_WP
  real(WP), parameter :: E2b_s = 25.46_WP * 1e3_WP / R_uni
  real(WP), parameter :: A3f_s = 1.13e16_WP
  real(WP), parameter :: n3f_s = -0.06_WP
  real(WP), parameter :: E3f_s = 476.05_WP * 1e3_WP / R_uni
  real(WP), parameter :: A3b_s = 4.17e13_WP / 1e6_WP
  real(WP), parameter :: n3b_s = 0.15_WP
  real(WP), parameter :: E3b_s = 0.00_WP * 1e3_WP / R_uni
  real(WP), parameter :: A4_s = 2.52e9_WP / 1e6_WP
  real(WP), parameter :: n4_s = 1.10_WP
  real(WP), parameter :: E4_s = 17.13_WP * 1e3_WP / R_uni
  real(WP), parameter :: A5_s = 2.20e12_WP / 1e6_WP
  real(WP), parameter :: n5_s = 0.00_WP
  real(WP), parameter :: E5_s = 31.38_WP * 1e3_WP / R_uni
  real(WP), parameter :: eff6 = 0.13_WP

contains

  subroutine soot_finitechem_sootstar(Y_,temp,dens,C_SootStar)
    use finitechem
    implicit none

    real(WP), dimension(N_tot), intent(in) :: Y_
    real(WP), intent(in) :: temp,dens
    real(WP), intent(out) :: C_SootStar
    real(WP) :: C_OH,C_H,C_H2O,C_H2,C_C2H2
    real(WP) :: k1f,k2f,k3f,k1b,k2b,k3b,k4

    C_OH = dens * Y_(i_OH-isc_sc+1) / (W_sp(i_OH-isc_sc+1)/1000.0_WP)
    C_H = dens * Y_(i_H-isc_sc+1) / (W_sp(i_H-isc_sc+1)/1000.0_WP)
    C_H2O = dens * Y_(i_H2O-isc_sc+1) / (W_sp(i_H2O-isc_sc+1)/1000.0_WP)
    C_H2 = dens * Y_(i_H2-isc_sc+1) / (W_sp(i_H2-isc_sc+1)/1000.0_WP)
    C_C2H2 = dens * Y_(i_C2H2-isc_sc+1) / (W_sp(i_C2H2-isc_sc+1)/1000.0_WP)
    
    k1f = A1f_s * temp**n1f_s * exp(-E1f_s/temp)
    k2f = A2f_s * temp**n2f_s * exp(-E2f_s/temp)
    k3f = A3f_s * temp**n3f_s * exp(-E3f_s/temp)
    k1b = A1b_s * temp**n1b_s * exp(-E1b_s/temp)
    k2b = A2b_s * temp**n2b_s * exp(-E2b_s/temp)
    k3b = A3b_s * temp**n3b_s * exp(-E3b_s/temp)
    k4 = A4_s * temp**n4_s * exp(-E4_s/temp)

    C_SootStar = (k1f*C_OH + k2f*C_H + k3f) / (k1b*C_H2O + k2b*C_H2 + k3b*C_H + k4*C_C2H2)
    C_SootStar = C_SootStar / (1.0_WP + C_SootStar)

    return
  end subroutine soot_finitechem_sootstar

  subroutine soot_finitechem_pahprodrates(Y_,temp,dens,srcSC_)
    use finitechem
    implicit none

    real(WP), intent(in), dimension(N_tot) :: Y_
    real(WP), intent(in) :: temp, dens
    real(WP), intent(inout), dimension(nscalar) :: srcSC_
    real(WP) :: aromconc,betadimer
    integer :: i

    do i=1,nPAH
       if (i_PAH(i).ne.-1) then
          aromconc = dens * Y_(i_PAH(i)-isc_sc+1) / (W_sp(i_PAH(i)-isc_sc+1)/1000.0_WP)
          betadimer = 4.0 * sqrt(2.0) * size_PAH(i)**(1.0_WP/6.0_WP) * &
               (6.0_WP/Pi)**(2.0_WP/3.0_WP) * AVOGADRO * &
               sqrt(Pi*R_uni*temp/(2.0_WP*AVOGADRO*SootDensity)) * &
               (MolarMassSoot/(AVOGADRO*SootDensity))**(1.0_WP/6.0_WP)
          srcSC_(i_PAH(i)) = srcSC_(i_PAH(i)) - dt_*(W_sp(i_PAH(i)-isc_sc+1)/1000.0_WP) * & 
               gamma_PAH(i) * betadimer * aromconc * aromconc
       end if
    end do

    return
  end subroutine soot_finitechem_pahprodrates

  subroutine soot_finitechem_sootrad(Y_,temp,srcSC_T,mom2)
    use finitechem
    implicit none

    real(WP), intent(in), dimension(N_tot) :: Y_
    real(WP), intent(in) :: temp,mom2
    real(WP), intent(inout) :: srcSC_T
    real(WP) :: alphas,fvo,Qrad,Cp_mix,W_gp
    real(WP), parameter :: Tb = 293.0_WP
    
    alphas = -3.75e5_WP + 1735.0_WP * temp ! [m^-1]
    fvo = mom2 * MolarMassSoot / SootDensity
    Qrad = -4.0_WP * alphas * fvo * 5.67051e-8_WP * (temp**4.0_WP - Tb**4.0_WP) ! [W/m^3]
      
    !Qrad = Qrad * 0.0_WP  ! This is fake.  SD
    
    call finitechem_Cp(Y_,temp,Cp_mix)
    
    if (compressible .or. is_constant_volume) then
       call finitechem_W(Y_, W_gp)
       Cp_mix = Cp_mix - R_cst / W_gp
    end if
    
    srcSC_T = srcSC_T + dt_ * Qrad / Cp_mix
    
    return
  end subroutine soot_finitechem_sootrad


end module soot_finitechem


subroutine soot_finitechem_init
  use soot_finitechem
  use finitechem
  implicit none

  integer :: i

  i_PAH = -1

  ! Find the PAH indices and surf/ox indices
!!$  do i = 1,N_tot
!!$     select case(trim(SC_name(i)))
!!$        case('A2-C10H8')
!!$           i_PAH(1) = i
!!$        case('A2R5-C12')
!!$           i_PAH(2) = i
!!$        case('P2-C12H1')
!!$           i_PAH(3) = i
!!$        case('A3-C14H1')
!!$           i_PAH(4) = i
!!$        case('A3R5-C16')
!!$           i_PAH(5) = i
!!$        case('A4-C16H1')
!!$           i_PAH(6) = i
!!$        case('FLTN-C16')
!!$           i_PAH(7) = i
!!$        case('A4R5-C18')
!!$           i_PAH(8) = i
!!$        case('OH')
!!$           i_OH = i
!!$        case('H2O')
!!$           i_H2O = i
!!$        case('H')
!!$           i_H = i
!!$        case('H2')
!!$           i_H2 = i
!!$        case('C2H2')
!!$           i_C2H2 = i
!!$        case('O2')
!!$           i_O2 = i
!!$        case('CO')
!!$           i_CO = i
!!$     end select
!!$  end do
  do i = 1,N_tot
     select case(trim(SC_name(i)))
        case('A2XC10H8')
           i_PAH(1) = i
        case('A2R5XC12')
           i_PAH(2) = i
        case('P2XC12H1')
           i_PAH(3) = i
        case('A3XC14H1')
           i_PAH(4) = i
        case('A3R5XC16')
           i_PAH(5) = i
        case('A4XC16H1')
           i_PAH(6) = i
        case('FLTNXC16')
           i_PAH(7) = i
        case('A4R5XC18')
           i_PAH(8) = i
        case('OH')
           i_OH = i
        case('H2O')
           i_H2O = i
        case('H')
           i_H = i
        case('H2')
           i_H2 = i
        case('C2H2')
           i_C2H2 = i
        case('O2')
           i_O2 = i
        case('CO')
           i_CO = i
     end select
  end do

  ! PAH carbon number
  size_PAH(1) = 10.0_WP
  size_PAH(2) = 12.0_WP
  size_PAH(3) = 12.0_WP
  size_PAH(4) = 14.0_WP
  size_PAH(5) = 16.0_WP
  size_PAH(6) = 16.0_WP
  size_PAH(7) = 16.0_WP
  size_PAH(8) = 18.0_WP
  
  ! PAH collision efficiencies
  gamma_PAH(1) = 2.0e-3_WP
  gamma_PAH(2) = 4.0e-3_WP
  gamma_PAH(3) = 8.5e-3_WP
  gamma_PAH(4) = 1.5e-2_WP
  gamma_PAH(5) = 2.5e-2_WP
  gamma_PAH(6) = 2.5e-2_WP
  gamma_PAH(7) = 2.5e-2_WP
  gamma_PAH(8) = 3.9e-2_WP

  ! SD
!!$  gamma_PAH = gamma_PAH * 10.0_WP
!!$  print*,'Gamma',gamma_PAH
  ! SD END

  return
end subroutine soot_finitechem_init


subroutine soot_finitechem_dimerprodrate(Y_,temp,dens,prodRate)
  use soot_finitechem
  use finitechem
  implicit none

  real(WP), dimension(N_tot), intent(in) :: Y_
  real(WP), intent(in) :: temp,dens
  real(WP), intent(out) :: prodRate
  integer :: i
  real(WP) :: aromconc,betadimer

  prodRate = 0.0_WP

  do i=1,nPAH
     if (i_PAH(i).ne.-1) then
        aromconc = dens * Y_(i_PAH(i)-isc_sc+1) / (W_sp(i_PAH(i)-isc_sc+1)/1000.0_WP)
        betadimer = 4.0 * sqrt(2.0) * size_PAH(i)**(1.0_WP/6.0_WP) * &
             (6.0_WP/Pi)**(2.0_WP/3.0_WP) * AVOGADRO * &
             sqrt(Pi*R_uni*temp/(2.0_WP*AVOGADRO*SootDensity)) * &
             (MolarMassSoot/(AVOGADRO*SootDensity))**(1.0_WP/6.0_WP)
        prodRate = prodRate + 0.5_WP * gamma_PAH(i) * & 
             betadimer * aromconc * aromconc
     end if
  end do
     
  return
end subroutine soot_finitechem_dimerprodrate


subroutine soot_finitechem_rhodot(rhodot,Y_,temp,dens)
  use soot_finitechem
  use finitechem
  implicit none

  real(WP), intent(out) :: rhodot
  real(WP), dimension(N_tot), intent(in) :: Y_
  real(WP), intent(in) :: temp,dens
  integer :: i
  real(WP) :: aromconc,betadimer

  rhodot = 0.0_WP

  do i=1,nPAH
     if (i_PAH(i).ne.-1) then
        aromconc = dens * Y_(i_PAH(i)-isc_sc+1) / (W_sp(i_PAH(i)-isc_sc+1)/1000.0_WP)
        betadimer = 4.0 * sqrt(2.0) * size_PAH(i)**(1.0_WP/6.0_WP) * &
             (6.0_WP/Pi)**(2.0_WP/3.0_WP) * AVOGADRO * &
             sqrt(Pi*R_uni*temp/(2.0_WP*AVOGADRO*SootDensity)) * &
             (MolarMassSoot/(AVOGADRO*SootDensity))**(1.0_WP/6.0_WP)
        rhodot = rhodot - (W_sp(i_PAH(i)-isc_sc+1)/1000.0_WP) * & 
             gamma_PAH(i) * betadimer * aromconc * aromconc
     end if
  end do
     
  return
end subroutine soot_finitechem_rhodot


!!$subroutine soot_finitechem_gasprodrates(SC_)
!!$  use soot_finitechem
!!$  use finitechem
!!$  use time_info
!!$  implicit none
!!$
!!$  real(WP), intent(inout), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: SC_
!!$  integer :: i,j,k,n
!!$  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: rhodot
!!$
!!$  do k=kmino_,kmaxo_
!!$     do j=jmino_,jmaxo_
!!$        do i=imino_,imaxo_
!!$           call soot_finitechem_pahprodrates(SC_(i,j,k,isc_sc:isc_sc+N_tot-1),T(i,j,k),RHO(i,j,k),pahprodrate)
!!$
!!$           do n=1,nPAH
!!$              if (i_PAH(n).ne.-1) then
!!$                 SC_(i,j,k,i_PAH(n)) = SC_(i,j,k,i_PAH(n)) + dt*pahprodrate(n)/RHO(i,j,k)
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  call soot_finitechem_rhodot(rhodot,SC_(:,:,:,isc_sc:isc_sc+N_tot-1),SC_(:,:,:,isc_T),RHO,nxo_*nyo_*nzo_)
!!$  SC_(:,:,:,isc_T) = SC_(:,:,:,isc_T) - dt*SC_(:,:,:,isc_T)*rhodot/RHO
!!$
!!$  return
!!$end subroutine soot_finitechem_gasprodrates


!!$subroutine soot_finitechem_nbrC(nbrC_)
!!$  use soot_finitechem
!!$  implicit none
!!$
!!$  real(WP), intent(out) :: nbrC_
!!$
!!$  nbrC_ = nbrC
!!$
!!$  return
!!$end subroutine soot_finitechem_nbrC


!!$subroutine soot_finitechem_nbrH(nbrH_)
!!$  use soot_finitechem
!!$  implicit none
!!$
!!$  real(WP), intent(out) :: nbrH_
!!$
!!$  nbrH_ = nbrH
!!$
!!$  return
!!$end subroutine soot_finitechem_nbrH



subroutine soot_finitechem_ksg(Y_,temp,dens,ksg)
  use soot_finitechem
  use finitechem
  implicit none

  real(WP), intent(in) :: temp,dens
  real(WP), dimension(N_tot), intent(in) :: Y_
  real(WP), intent(out) :: ksg
  real(WP) :: C_C2H2, C_SootStar

  C_C2H2 = dens * Y_(i_C2H2-isc_sc+1) / (W_sp(i_C2H2-isc_sc+1)/1000.0_WP)
  
  call soot_finitechem_sootstar(Y_,temp,dens,C_SootStar)

  ksg = A4_s * temp**n4_s * exp(-E4_s/temp) * C_C2H2 * C_SootStar

  return
end subroutine soot_finitechem_ksg


subroutine soot_finitechem_kox(Y_,temp,dens,kox)
  use soot_finitechem
  use finitechem
  implicit none

  real(WP), intent(in) :: temp,dens
  real(WP), dimension(N_tot), intent(in) :: Y_
  real(WP), intent(out) :: kox
  real(WP) :: k6, C_OH, C_O2, C_SootStar

  C_OH = dens * Y_(i_OH-isc_sc+1) / (W_sp(i_OH-isc_sc+1)/1000.0_WP)
  C_O2 = dens * Y_(i_O2-isc_sc+1) / (W_sp(i_O2-isc_sc+1)/1000.0_WP)

  call soot_finitechem_sootstar(Y_,temp,dens,C_SootStar)

  k6 = 8.94_WP * eff6 * sqrt(temp) * AVOGADRO
  kox = A5_s * temp**n5_s * exp(-E5_s/temp) * C_O2 * C_SootStar + & 
         (0.5_WP/chisoot) * k6 * C_OH

  return
end subroutine soot_finitechem_kox


subroutine soot_finitechem_ko2(Y_,temp,dens,ko2)
  use soot_finitechem
  use finitechem
  implicit none

  real(WP), intent(in) :: temp,dens
  real(WP), dimension(N_tot), intent(in) :: Y_
  real(WP), intent(out) :: ko2
  real(WP) :: C_O2, C_SootStar

  C_O2 = dens * Y_(i_O2-isc_sc+1) / (W_sp(i_O2-isc_sc+1)/1000.0_WP)

  call soot_finitechem_sootstar(Y_,temp,dens,C_SootStar)

  ko2 = A5_s * temp**n5_s * exp(-E5_s/temp) * C_O2 * C_SootStar

  return
end subroutine soot_finitechem_ko2
