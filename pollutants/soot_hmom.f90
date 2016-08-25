module soot_hmom
  use soot
  use pah
  use infnan
  implicit none

  ! ------------ Constants ------------

  ! Soot relevant constants
  real(WP), parameter :: SootDensity = 1800.0_WP    ! kg/m^3
  real(WP), parameter :: MolarMassSoot = 12.0e-3_WP ! kg/mol
  real(WP), parameter :: Avogadro = 6.022e23_WP     ! 1/mol
  real(WP), parameter :: Df = 1.8_WP                ! 1
  real(WP), parameter :: chisoot = 1.7e19_WP        ! 1/m^2

  ! Particle sizes
  real(WP) :: DIMER_NBRC,NUCL_NBRC,NUCL_SURF

  ! Derived constants
  real(WP) :: CarbonToDiam
  real(WP) :: Cfm
  real(WP) :: Ccn
  real(WP) :: lambda
  real(WP) :: av,as

  ! Approximation of small surface area change
  real(WP) :: FitC
  real(WP) :: FitExp

  ! Clipping for HMOM
  real(WP), parameter :: SMALLWEIGHT = 1.0e-20_WP
  real(WP), parameter :: BIG_NP = 1.0e6_WP

  ! ----------  Local environment -----------

  ! Relevant quantities to soot formation
  real(WP) :: DIMER_conc
  real(WP) :: dens,sqrtT,T_MU,MUsqrtW_RHOsqrtT
  real(WP), dimension(:), pointer :: MassFracs
  real(WP) :: totdivg

  !$OMP THREADPRIVATE(DIMER_conc,dens,sqrtT,T_MU,MUsqrtW_RHOsqrtT,MassFracs,totdivg)

  ! Quantities from the gas phase
  real(WP) :: wCoeff,oxCoeff,o2Coeff
  real(WP) :: prodRate

  !$OMP THREADPRIVATE(wCoeff,oxCoeff,o2Coeff,prodRate)

  ! Quantities for gas-phase source terms
  real(WP) :: ypah, rhodot, zsrc, zzsrc, z2rhodot, enthdot, sootrad

  !$OMP THREADPRIVATE(ypah,rhodot,zsrc,zzsrc,z2rhodot,enthdot,sootrad)

  ! ------------- HMOM -------------

  ! HMOM moment indices
  real(WP), dimension(:),   pointer :: frac
  real(WP), dimension(:),   pointer :: kmom
  real(WP), dimension(:,:), pointer :: moments

  ! HMOM moments
  real(WP), dimension(:), pointer :: mom, old_mom

  ! HMOM source terms
  real(WP), dimension(:),   pointer :: src_mom

  !$OMP THREADPRIVATE(kmom,mom,old_mom,src_mom)

contains

  ! Create the moments
  ! ------------------
  subroutine soot_define_moments
    implicit none

    allocate(moments(nMoments,dim))

    if (nMoments.ge.4) then
       ! Total number
       moments(1,1) = 0.0_WP
       moments(1,2) = 0.0_WP
       ! Total volume
       moments(2,1) = 1.0_WP
       moments(2,2) = 0.0_WP
       ! Total surface
       moments(3,1) = 0.0_WP
       moments(3,2) = 1.0_WP
    end if
    if (nMoments.ge.7) then
       ! Total volume*volume
       moments(4,1) = 2.0_WP
       moments(4,2) = 0.0_WP
       ! Total volume*surface
       moments(5,1) = 1.0_WP
       moments(5,2) = 1.0_WP
       ! Total surface*surface
       moments(6,1) = 0.0_WP
       moments(6,2) = 2.0_WP
    end if

    return
  end subroutine soot_define_moments

  ! Compute the fractional moment
  ! -----------------------------
  function soot_fracmom(i,j)
    use soot
    implicit none
    real(WP) :: soot_fracmom
    real(WP) :: m00, m10, m01, m20, m11, m02
    real(WP) :: f1,f2,f3,f4,f5,f6
    real(WP), intent(in) :: i,j 

    m00 = mom(1) - mom(nMoments)
    m10 = mom(2)-NUCL_NBRC*mom(nMoments)
    m01 = mom(3)-NUCL_SURF*mom(nMoments)

    if (nMoments.eq.4) then
       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP) then
          soot_fracmom = mom(1) * (mom(2)/mom(1))**i * (mom(3)/mom(1))**j
       else
          soot_fracmom = mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j + m00*(m10/m00)**i*(m01/m00)**j
       end if
    elseif (nMoments.eq.7) then
       m20 = mom(4) - (NUCL_NBRC*NUCL_NBRC)*mom(nMoments)
       m11 = mom(5) - (NUCL_NBRC*NUCL_SURF)*mom(nMoments)
       m02 = mom(6) - (NUCL_SURF*NUCL_SURF)*mom(nMoments)

       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP .or. &
            m20.lt.1e-25_WP .or. m11.lt.1e-25_WP .or. m02.lt.1e-25) then
          f1 = mom(1)
          f2 = (mom(2)**(2.0_WP)*mom(1)**(-1.5_WP)*mom(4)**(-0.5_WP))**i
          f3 = (mom(3)**(2.0_WP)*mom(1)**(-1.5_WP)*mom(6)**(-0.5_WP))**j
          f4 = (mom(4)**(0.5_WP)*mom(1)**(0.5_WP)*mom(2)**(-1.0_WP))**(i*i)
          f5 = (mom(5)*mom(1)/mom(2)/mom(3))**(i*j)
          f6 = (mom(6)**(0.5_WP)*mom(1)**(0.5_WP)*mom(3)**(-1.0_WP))**(j*j)
          soot_fracmom = f1*f2*f3*f4*f5*f6
       else
          f1 = m00
          f2 = (m10**(2.0_WP)*m00**(-1.5_WP)*m20**(-0.5_WP))**i
          f3 = (m01**(2.0_WP)*m00**(-1.5_WP)*m02**(-0.5_WP))**j
          f4 = (m20**(0.5_WP)*m00**(0.5_WP)*m10**(-1.0_WP))**(i*i)
          f5 = (m11*m00/m10/m01)**(i*j)
          f6 = (m02**(0.5_WP)*m00**(0.5_WP)*m01**(-1.0_WP))**(j*j)
          soot_fracmom = mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j + f1*f2*f3*f4*f5*f6
       end if
    end if

    return
  end function soot_fracmom

  ! Compute the fractional moment of the second mode
  ! ------------------------------------------------
  function soot_fracmomlarge(i,j)
    use soot
    implicit none
    real(WP) :: soot_fracmomlarge
    real(WP) :: m00, m10, m01, m20, m11, m02
    real(WP), intent(in) :: i,j 

    m00 = mom(1) - mom(nMoments)
    m10 = mom(2)-(NUCL_NBRC)*mom(nMoments)
    m01 = mom(3)-(NUCL_SURF)*mom(nMoments)

    if (nMoments.eq.4) then
       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP) then
          soot_fracmomlarge = 1.0e-60_WP * NUCL_NBRC**i * NUCL_SURF**j
       else
          soot_fracmomlarge = soot_fracmom(i,j) - mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j
       end if
    elseif (nMoments.eq.7) then
       m20 = mom(4) - (NUCL_NBRC*NUCL_NBRC)*mom(nMoments)
       m11 = mom(5) - (NUCL_NBRC*NUCL_SURF)*mom(nMoments)
       m02 = mom(6) - (NUCL_SURF*NUCL_SURF)*mom(nMoments)

       if (m00.lt.1e-25_WP .or. m10.lt.1e-25_WP .or. m01.lt.1e-25_WP .or. &
            m20.lt.1e-25_WP .or. m11.lt.1e-25_WP .or. m02.lt.1e-25) then
          soot_fracmomlarge = 1.0e-60_WP * NUCL_NBRC**i * NUCL_SURF**j
       else
          soot_fracmomlarge = soot_fracmom(i,j) - mom(nMoments)*NUCL_NBRC**i*NUCL_SURF**j
       end if
    end if 

    return
  end function soot_fracmomlarge


  ! Compute the collision coefficients
  ! ----------------------------------
  function soot_beta_nucl()
    use soot
    implicit none
    real(WP) :: soot_beta_nucl

    ! 2.2_WP * 4.0_WP * sqrt(2.0_WP) = 12.45
    soot_beta_nucl = 12.45_WP * Cfm * sqrtT * DIMER_NBRC**(1.0_WP/6.0_WP)

    return
  end function soot_beta_nucl

  function soot_beta_cond()
    use soot
    
    implicit none
    real(WP) :: soot_beta_cond
    real(WP) :: dV,S0v,omega

    ! Compute the subfilter intermittency and normalize the moments
    omega = 0.0_WP

    if (use_soot_sgs) omega = min(max(1.0_WP - (mom(1)**2/mom(nMoments+1)),0.0_WP),1.0_WP-1e-8_WP) ! SD 1.0_WP - (mom(1)**2/mom(nMoments+1)
    mom = mom / (1.0_WP-omega)

    dV = DIMER_NBRC
  
    S0v =          soot_fracmom(2.0_WP*av       , 2.0_WP*as       ) * dV**(-3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(       av       ,        as       ) * dV**(-1.0_WP/6.0_WP) + &
                   soot_fracmom(          0.0_WP,           0.0_WP) * dV**( 1.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(2.0_WP*av-1.0_WP, 2.0_WP*as       ) * dV**( 3.0_WP/6.0_WP) + &
                   soot_fracmom(       av-1.0_WP,        as       ) * dV**( 5.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(         -1.0_WP,           0.0_WP) * dV**( 7.0_WP/6.0_WP)

    soot_beta_cond = (1.0_WP-omega) * Cfm * sqrtT * S0v

    ! Unnormalize the moments
    mom = mom * (1.0_WP-omega)

    return
  end function soot_beta_cond


  ! Compute Nucleation source term for moment k
  ! -------------------------------------------
  subroutine soot_nucleation(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src

    src = 0.5_WP * soot_beta_nucl() * DIMER_conc**2 * &
          NUCL_NBRC**k(1) * NUCL_SURF**k(2)

    return
  end subroutine soot_nucleation

  ! Compute Nucleation source term for weight of first mode
  ! -------------------------------------------------------
  subroutine soot_nucleationsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src

    src = 0.5_WP * soot_beta_nucl() * DIMER_conc**2

    return
  end subroutine soot_nucleationsmall

  ! Compute Coagulation source term for moment k
  ! --------------------------------------------
  subroutine soot_coagulation(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: V0,S0,Vi,Si
    real(WP) :: ss_fm,ss_cn,ss
    real(WP) :: sl_fm,sl_cn,sl
    real(WP) :: ll_fm,ll_cn,ll
    real(WP) :: slip

    V0 = NUCL_NBRC
    S0 = NUCL_SURF
    Vi = soot_fracmomlarge(1.0_WP,0.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)
    Si = soot_fracmomlarge(0.0_WP,1.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)

    slip = lambda * MUsqrtW_RHOsqrtT

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       ss = 0.0_WP
    else
       ss_fm = 2.2_WP * Cfm * sqrtT * 2.0_WP**2.5_WP * (2.0_WP**(k(1)+2.0_WP/3.0_WP*k(2)-1.0_WP) - 1.0_WP) * & 
               V0**(k(1)+2.0_WP/3.0_WP*k(2)+1.0_WP/6.0_WP) * mom(nMoments) * mom(nMoments)
       ss_cn = 4.0_WP * Ccn * T_MU * (1.0_WP + 1.257_WP*slip*V0**(-1.0_WP/3.0_WP)) * (2.0_WP**(k(1)+2.0_WP/3.0_WP*k(2)-1.0_WP) - 1.0_WP) * &
               V0**(k(1)+2.0_WP/3.0_WP*k(2)) * mom(nMoments) * mom(nMoments)
       ss = (ss_fm * ss_cn) / (ss_fm + ss_cn)
    end if

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       sl = 0.0_WP
    else
       if ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          sl_fm = -soot_psi_sl(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
          sl_cn = -Ccn*T_MU * (2.000_WP*                                 mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                             NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                             NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                               1.257_WP*slip*NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                               1.257_WP*slip*                            mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                               1.257_WP*slip*NUCL_NBRC**(-2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                               1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))

       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          sl_fm = FitC*soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+1.0_WP,1.0_WP,0.0_WP) - &
                       soot_psi_sl(               0.0_WP,              0.0_WP,0.0_WP,1.0_WP)
          sl_cn = Ccn*T_MU * ( FitC * ( &
                                       2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+1.0_WP) &
                             ) - ( &
                                       2.000_WP*     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                     NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                       1.257_WP*slip*                            mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       elseif ((k(1).lt.2.01_WP .and. k(1).gt.1.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          sl_fm = 2.0_WP*soot_psi_sl(1.0_WP,0.0_WP,1.0_WP,0.0_WP)
          sl_cn = 2.0_WP*Ccn*T_MU * (2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                                                   NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                                   NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                                     1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          +1.0_WP, 0.0_WP   ) + &
                                     1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                     1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                                     1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as))
       elseif ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          sl_fm = FitC*soot_psi_sl(-2.0_WP*FitExp       ,3.0_WP*FitExp+1.0_WP,1.0_WP,0.0_WP) + &
                       soot_psi_sl(               0.0_WP,              1.0_WP,1.0_WP,0.0_WP) + &
                  FitC*soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+1.0_WP,2.0_WP,0.0_WP) - &
                       soot_psi_sl(               0.0_WP,              0.0_WP,1.0_WP,1.0_WP)
          sl_cn = Ccn*T_MU * ( FitC * ( &
                                       2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp,           3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp,-       as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp,           3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp,-       as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp,-2.0_WP*as+3.0_WP*FitExp+1.0_WP) &
                             ) + ( &
                                       2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                                     NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP) &
                             ) + FitC * ( &
                                       2.000_WP*     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                                     NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+1.0_WP) + &
                                       1.257_WP*slip*NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+1.0_WP) &
                             ) - ( &
                                       2.000_WP*     NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                       1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                       1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                       1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.2.01_WP .and. k(2).gt.1.99_WP)) then
          sl_fm = 2.0_WP*FitC*     soot_psi_sl(-2.0_WP*FitExp-1.0_WP,3.0_WP*FitExp+2.0_WP,1.0_WP,0.0_WP) + &
                         FitC*FitC*soot_psi_sl(-4.0_WP*FitExp-2.0_WP,6.0_WP*FitExp+2.0_WP,2.0_WP,0.0_WP) - &
                                   soot_psi_sl(               0.0_WP,              0.0_WP,0.0_WP,2.0_WP)
          sl_cn = 2.0_WP*Ccn*T_MU * ( 2.0_WP*FitC * ( &
                                              2.000_WP*     NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -2.0_WP*FitExp-1.0_WP,           3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(-       av-2.0_WP*FitExp-1.0_WP,-       as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-2.0_WP*FitExp-1.0_WP,        as+3.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-2.0_WP*FitExp-1.0_WP,-2.0_WP*as+3.0_WP*FitExp+2.0_WP) &
                                    ) + FitC*FitC * ( &
                                              2.000_WP*     NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(          -4.0_WP*FitExp-2.0_WP,           6.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av-4.0_WP*FitExp-2.0_WP,-       as+6.0_WP*FitExp+2.0_WP) + &
                                                            NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-4.0_WP*FitExp-2.0_WP,        as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(          -4.0_WP*FitExp-2.0_WP,           6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**2               *mom(nMoments)*soot_fracmomlarge(-       av-4.0_WP*FitExp-2.0_WP,-       as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av-4.0_WP*FitExp-2.0_WP,        as+6.0_WP*FitExp+2.0_WP) + &
                                              1.257_WP*slip*NUCL_NBRC**( 7.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av-4.0_WP*FitExp-2.0_WP,-2.0_WP*as+6.0_WP*FitExp+2.0_WP) &
                                    ) - ( &
                                              2.000_WP*     NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                            NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                                            NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                              1.257_WP*slip*NUCL_NBRC                  *mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                              1.257_WP*slip*NUCL_NBRC**( 4.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                              1.257_WP*slip*NUCL_NBRC**( 2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                                              1.257_WP*slip*NUCL_NBRC**( 5.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)))
       endif
       sl = (sl_fm * sl_cn) / (sl_fm + sl_cn)
    end if

    if ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
       ll = 0.0_WP
    elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
       ll = 0.0_WP
    else
       if ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          ll_fm = -0.5_WP*soot_psi_ll(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
          ll_cn = -0.5*Ccn*T_MU * (2.000_WP*     soot_fracmomlarge( 0.0_WP   , 0.0_WP   )*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                                 soot_fracmomlarge(        av,        as)*soot_fracmomlarge(-       av,-       as) + &
                                                 soot_fracmomlarge(-       av,-       as)*soot_fracmomlarge(        av,        as) + &
                                   1.257_WP*slip*soot_fracmomlarge(-       av,-       as)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                   1.257_WP*slip*soot_fracmomlarge( 0.0_WP   , 0.0_WP   )*soot_fracmomlarge(-       av,-       as) + &
                                   1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as)*soot_fracmomlarge(        av,        as) + &
                                   1.257_WP*slip*soot_fracmomlarge(        av,        as)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))
       elseif ((k(1).lt.2.01_WP .and. k(1).gt.1.99_WP) .and. (k(2).lt.0.01_WP .and. k(2).gt.-0.01_WP)) then
          ll_fm = soot_psi_ll(1.0_WP,0.0_WP,1.0_WP,0.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                                            soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                                            soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(           1.0_WP, 0.0_WP   ) + &
                              1.257_WP*slip*soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(-       av+1.0_WP,-       as) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as)*soot_fracmomlarge(        av+1.0_WP,        as) + &
                              1.257_WP*slip*soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as))
       elseif ((k(1).lt.1.01_WP .and. k(1).gt.0.99_WP) .and. (k(2).lt.1.01_WP .and. k(2).gt.0.99_WP)) then
          ll_fm = soot_psi_ll(1.0_WP,0.0_WP,0.0_WP,1.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                            soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                            soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av+1.0_WP,-       as)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(           1.0_WP, 0.0_WP   )*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av+1.0_WP,-2.0_WP*as)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(        av+1.0_WP,        as)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP))
       elseif ((k(1).lt.0.01_WP .and. k(1).gt.-0.01_WP) .and. (k(2).lt.2.01_WP .and. k(2).gt.1.99_WP)) then
          ll_fm = soot_psi_ll(0.0_WP,1.0_WP,0.0_WP,1.0_WP)
          ll_cn = Ccn*T_MU * (2.000_WP*     soot_fracmomlarge( 0.0_WP   ,           1.0_WP)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                                            soot_fracmomlarge(        av,        as+1.0_WP)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                                            soot_fracmomlarge(-       av,-       as+1.0_WP)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-       av,-       as+1.0_WP)*soot_fracmomlarge( 0.0_WP   ,           1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge( 0.0_WP   ,           1.0_WP)*soot_fracmomlarge(-       av,-       as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP)*soot_fracmomlarge(        av,        as+1.0_WP) + &
                              1.257_WP*slip*soot_fracmomlarge(        av,        as+1.0_WP)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as+1.0_WP))
       endif
       ll = (ll_fm * ll_cn) / (ll_fm + ll_cn)
    end if

    src = ss + sl + ll

    return
  end subroutine soot_coagulation

  ! Compute Coagulation source term for weight of first mode
  ! --------------------------------------------------------
  subroutine soot_coagulationsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: V0,Vi,Si
    real(WP) :: ss_fm,ss_cn,ss
    real(WP) :: sl_fm,sl_cn,sl
    real(WP) :: slip

    V0 = NUCL_NBRC
    Vi = soot_fracmomlarge(1.0_WP,0.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)
    Si = soot_fracmomlarge(0.0_WP,1.0_WP)/soot_fracmomlarge(0.0_WP,0.0_WP)

    slip = lambda * MUsqrtW_RHOsqrtT

    ss_fm = -2.2_WP * Cfm * sqrtT * 2.0_WP**2.5_WP * V0**(1.0_WP/6.0_WP) * mom(nMoments) * mom(nMoments)
    ss_cn = -4.0_WP * Ccn * T_MU * (1.0_WP + 1.257_WP*slip*V0**(-1.0_WP/3.0_WP)) * mom(nMoments) * mom(nMoments)
    ss = (ss_fm * ss_cn) / (ss_fm + ss_cn)

    sl_fm = -soot_psi_sl(0.0_WP,0.0_WP,0.0_WP,0.0_WP)
    sl_cn = -Ccn * T_MU * (2.000_WP*                                 mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                                         NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                                         NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                           1.257_WP*slip*NUCL_NBRC**(-1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge( 0.0_WP   , 0.0_WP   ) + &
                           1.257_WP*slip*                            mom(nMoments)*soot_fracmomlarge(-       av,-       as) + &
                           1.257_WP*slip*NUCL_NBRC**(-2.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(        av,        as) + &
                           1.257_WP*slip*NUCL_NBRC**( 1.0_WP/3.0_WP)*mom(nMoments)*soot_fracmomlarge(-2.0_WP*av,-2.0_WP*as))
    sl = (sl_fm * sl_cn) / (sl_fm + sl_cn)

    src = ss + sl

    return
  end subroutine soot_coagulationsmall

  ! Compute Condensation source term for moment k
  ! ---------------------------------------------
  subroutine soot_condensation(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV
    real(WP) :: S0v,S0s

    dV = DIMER_NBRC

    S0v =          soot_fracmom(k(1)+2.0_WP*av-1.0_WP,k(2)+2.0_WP*as) * dV**( 3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(k(1)+       av-1.0_WP,k(2)+       as) * dV**( 5.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)          -1.0_WP,k(2)          ) * dV**( 7.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)+2.0_WP*av-2.0_WP,k(2)+2.0_WP*as) * dV**( 9.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)+       av-2.0_WP,k(2)+       as) * dV**(11.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)          -2.0_WP,k(2)          ) * dV**(13.0_WP/6.0_WP)

    S0s =          soot_fracmom(k(1)-2.0_WP*FitExp+2.0_WP*av-1.0_WP, k(2)+3.0_WP*FitExp+2.0_WP*as) * dV**( 3.0_WP/6.0_WP) + &
          2.0_WP * soot_fracmom(k(1)-2.0_WP*FitExp+       av-1.0_WP, k(2)+3.0_WP*FitExp+       as) * dV**( 5.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)-2.0_WP*FitExp          -1.0_WP, k(2)+3.0_WP*FitExp          ) * dV**( 7.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)-2.0_WP*FitExp+2.0_WP*av-2.0_WP, k(2)+3.0_WP*FitExp+2.0_WP*as) * dV**( 9.0_WP/6.0_WP) + &
                   soot_fracmom(k(1)-2.0_WP*FitExp+       av-2.0_WP, k(2)+3.0_WP*FitExp+       as) * dV**(11.0_WP/6.0_WP) + &
          0.5_WP * soot_fracmom(k(1)-2.0_WP*FitExp          -2.0_WP, k(2)+3.0_WP*FitExp          ) * dV**(13.0_WP/6.0_WP)

    src = Cfm * sqrtT * (k(1) * S0v + FitC * k(2) * S0s) * DIMER_conc    

    return
  end subroutine soot_condensation

  ! Compute Condensation source term for weight of first mode
  ! ---------------------------------------------------------
  subroutine soot_condensationsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0

    dV = DIMER_NBRC
    V0 = NUCL_NBRC

    src = -Cfm * sqrtT * sqrt(1.0_WP/dV + 1.0_WP/V0) * (dV**(1.0_WP/3.0_WP) + V0**(1.0_WP/3.0_WP))**2.0_WP * DIMER_conc * mom(nMoments)

    return
  end subroutine soot_condensationsmall

  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfacereaction(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV,chicarb
    
    dV = 2.0_WP
    chicarb = chisoot * (36.0_WP*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = wCoeff * chicarb * dV * (k(1) * soot_fracmom(k(1)-1.0_WP, k(2)+1.0_WP) + &
          k(2) * FitC * soot_fracmom(k(1)-1.0_WP-2.0_WP*FitExp,k(2)+1.0_WP+3.0_WP*FitExp))
   
    return
  end subroutine soot_surfacereaction

  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfacereactionsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb

    dV = 2.0_WP
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = -wCoeff * chicarb * V0**(2.0_WP/3.0_WP) * mom(nMoments)

    return
  end subroutine soot_surfacereactionsmall

  ! Compute Surface Oxidation source term for moment k
  ! --------------------------------------------------
  subroutine soot_surfaceoxidation(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb,small,large

    dV = 2.0_WP
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0_WP*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    small = -oxCoeff * chicarb * dV * V0**(k(1)-1.0_WP+(2.0_WP/3.0_WP)*(k(2)+1.0_WP)) * mom(nMoments)
    large = -oxCoeff * chicarb * dV * (k(1)+(2.0_WP/3.0_WP)*k(2)) * soot_fracmomlarge(k(1)-1.0_WP,k(2)+1.0_WP)
   
    src = small + large

    return
  end subroutine soot_surfaceoxidation

  ! Compute Surface Reaction source term for moment k
  ! -------------------------------------------------
  subroutine soot_surfaceoxidationsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: dV,V0,chicarb,small,inter,large

    dV = 2.0_WP
    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    ! Oxidation of the small particles
    small = -oxCoeff * chicarb * dV * V0**(-1.0_WP/3.0_WP) * mom(nMoments)

    ! Oxidation of the larger particles as they become small
    ! Proportionality constant is the ratio between the volumes of the two modes
    inter = V0 / (soot_fracmomlarge(1.0_WP,0.0_WP) / soot_fracmomlarge(0.0_WP,0.0_WP))
    large = inter * oxCoeff * chicarb * dV * soot_fracmomlarge(-1.0_WP,1.0_WP)

    src = small + large

    return
  end subroutine soot_surfaceoxidationsmall

  ! Compute Fragmentation source term for moment k
  ! ----------------------------------------------
  subroutine soot_fragmentation(k,src)
    use soot
    implicit none

    real(WP), dimension(dim), intent(in) :: k
    real(WP), intent(out) :: src
    real(WP) :: chicarb

    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    src = (2.0_WP**(1.0_WP-k(1)-k(2))-1.0_WP) * 2.0_WP*o2Coeff*chicarb*soot_fracmomlarge(k(1)-1.0_WP,k(2)+1.0_WP)

    return
  end subroutine soot_fragmentation

  ! Compute Fragmentation source term for delta function
  ! ----------------------------------------------------
  subroutine soot_fragmentationsmall(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    real(WP) :: V0,chicarb,inter

    V0 = NUCL_NBRC
    chicarb = chisoot * (36.0*Pi)**(1.0_WP/3.0_WP) * (MolarMassSoot/(Avogadro*SootDensity))**(2.0_WP/3.0_WP)

    inter = V0 / (soot_fracmomlarge(1.0_WP,0.0_WP) / soot_fracmomlarge(0.0_WP,0.0_WP))

    src = inter * 2.0_WP*o2Coeff*chicarb*soot_fracmomlarge(-1.0_WP,1.0_WP)

    return
  end subroutine soot_fragmentationsmall

  ! Compute the thermophoretic source term for the number density squared
  subroutine soot_squared_thermophoresis(src)
    use soot
    implicit none

    real(WP), intent(out) :: src
    
    src = mom(nMoments+1) * totdivg

    return
  end subroutine soot_squared_thermophoresis


  ! Compute PsiSL for free molecular coagulation source term
  ! --------------------------------------------------------
  function soot_psi_sl(x,y,a,b)
    use soot
    implicit none

    real(WP), intent(in) :: x,y,a,b
    real(WP) :: soot_psi_sl
    real(WP) :: psi1,psi2

    psi1 =          NUCL_NBRC**( 1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(-1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(-3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y)

    psi2 =          NUCL_NBRC**( 7.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(         -0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**( 5.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av-0.5_WP+x,       as+y) + &
                    NUCL_NBRC**( 3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) + &
                    NUCL_NBRC**( 1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(          0.5_WP+x,          y) + &
           2.0_WP * NUCL_NBRC**(-1.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(       av+0.5_WP+x,       as+y) + &
                    NUCL_NBRC**(-3.0_WP/6.0_WP+a+(2.0_WP/3.0_WP)*b)*mom(nMoments) * soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y)

    soot_psi_sl = 2.2_WP * Cfm * sqrtT * sqrt(psi1*psi2)

    return
  end function soot_psi_sl


  ! Compute PsiLL for free molecular coagulation source term
  ! --------------------------------------------------------
  function soot_psi_ll(x,y,a,b)
    use soot
    implicit none

    real(WP), intent(in) :: x,y,a,b
    real(WP) :: soot_psi_ll
    real(WP) :: psi1,psi2

    psi1 =          soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b)

    psi2 =          soot_fracmomlarge(2.0_WP*av+0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(         -0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av+0.5_WP+x,       as+y) * soot_fracmomlarge(       av-0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         +0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av-0.5_WP+a,2.0_WP*as+b) + &
                    soot_fracmomlarge(2.0_WP*av-0.5_WP+x,2.0_WP*as+y) * soot_fracmomlarge(          0.5_WP+a,          b) + &
           2.0_WP * soot_fracmomlarge(       av-0.5_WP+x,       as+y) * soot_fracmomlarge(       av+0.5_WP+a,       as+b) + &
                    soot_fracmomlarge(         -0.5_WP+x,          y) * soot_fracmomlarge(2.0_WP*av+0.5_WP+a,2.0_WP*as+b)

    soot_psi_ll = 2.2_WP * Cfm * sqrtT * sqrt(psi1*psi2)

    return
  end function soot_psi_ll


  ! Compute PAH concentration
  ! -------------------------
  subroutine soot_pah_molecules
    implicit none

    real(WP) :: betaN,betaC,delta

    ! Nucleation collision coeff
    betaN = soot_beta_nucl()

    ! Condensation collision coeff
    betaC = soot_beta_cond()

    ! Solve quadratic equation
    delta = betaC**2 + 4.0_WP*betaN*prodRate
    DIMER_conc = (sqrt(delta)-betaC)/(2.0_WP*betaN)

    return
  end subroutine soot_pah_molecules

  ! Dispersed-phase source terms
  ! ----------------------------
  subroutine soot_hmom_source_solid
    implicit none

    integer :: i
    real(WP) :: omega, src_tmp

    
    ! Compute dimer concentration
    call soot_pah_molecules 


    ! Compute the subfilter intermittency and normalized moments
    omega = 0.0_WP
    if (use_soot_sgs) omega = min(max(1.0_WP - (mom(1)**2/mom(nMoments+1)),0.0_WP),1.0_WP-1e-8_WP) ! SD 1.0_WP - (mom(1)**2/mom(nMoments+1))
    mom = mom / (1.0_WP-omega)


    ! Source terms for moments
    do i=1,nMoments-1
       kmom = moments(i,:)
       src_mom(i) = 0.0_WP
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = 0.0_WP
       ! Nucleation
       call soot_nucleation(kmom,src_tmp)
       src_mom(i) = src_mom(i) + src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = & 
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
       ! Coagulation
       call soot_coagulation(kmom,src_tmp)
       src_mom(i) = src_mom(i) + (1.0_WP-omega)*src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = & 
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
       ! Condensation
       call soot_condensation(kmom,src_tmp)
       src_mom(i) = src_mom(i) + (1.0_WP-omega)*src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = & 
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
       ! Surface Reaction
       call soot_surfacereaction(kmom,src_tmp)
       src_mom(i) = src_mom(i) + (1.0_WP-omega)*src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = &
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
       ! Surface Oxidation
       call soot_surfaceoxidation(kmom,src_tmp)
       src_mom(i) = src_mom(i) + (1.0_WP-omega)*src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = & 
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
       ! Fragmentation
       call soot_fragmentation(kmom,src_tmp)
       src_mom(i) = src_mom(i) + (1.0_WP-omega)*src_tmp
       if (use_soot_sgs .and. i.eq.1) src_mom(nMoments+1) = & 
            src_mom(nMoments+1) + 2.0_WP*(1.0_WP-omega)*mom(1)*src_tmp
    end do

    ! Source term for weight of first peak
       src_mom(nMoments) = 0.0_WP
       ! Nucleation
       call soot_nucleationsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
       ! Coagulation
       call soot_coagulationsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
       ! Condensation
       call soot_condensationsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
       ! Surface Reaction
       call soot_surfacereactionsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
       ! Surface Oxidation
       call soot_surfaceoxidationsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
       ! Fragmentation
       call soot_fragmentationsmall(src_tmp)
       src_mom(nMoments) = src_mom(nMoments) + (1.0_WP-omega)*src_tmp
    !!!

    ! Get the unnormalized moments back
    mom = mom * (1.0_WP-omega)

    ! Thermophoresis source term for the number density squared
    if (use_soot_sgs) then
       call soot_squared_thermophoresis(src_tmp)
       src_mom(nMoments+1) = src_mom(nMoments+1) + src_tmp
    end if

    return
  end subroutine soot_hmom_source_solid

  ! Gas-phase source terms
  ! ----------------------
  subroutine soot_hmom_source_gas(SC_,srcSC_,srcP_)
    use time_info
    use chemtable
    use finitechem
    use soot_finitechem
    implicit none

    real(WP), dimension(nScalar), intent(in) :: SC_
    real(WP), dimension(nScalar), intent(inout) :: srcSC_
    real(WP), intent(inout) :: srcP_
    
    ! Continuity
    srcP_ = srcP_ + dt_*rhodot

    select case(trim(chemistry))
       case ('chemtable')
          ! Mixture Fraction
          srcSC_(isc_ZMIX) = srcSC_(isc_ZMIX) + dt_*zsrc*dens

          ! Mixture Fraction Squared
          if (isc_ZMIX2.gt.0) then
             ! Mixture fraction source term
             srcSC_(isc_ZMIX2) = srcSC_(isc_ZMIX2) + dt_*2.0_WP*zzsrc*dens
             
             ! Density source term
             srcSC_(isc_ZMIX2) = srcSC_(isc_ZMIX2) - dt_*z2rhodot*dens
          end if
          
          ! Mixture Fraction Variance
          if (isc_ZVAR.gt.0) then
             ! Mixture fraction source term
             srcSC_(isc_ZVAR) = srcSC_(isc_ZVAR) + dt_*2.0_WP*(zzsrc-zsrc*SC_(isc_ZMIX))*dens
             
             ! Density source term
             srcSC_(isc_ZVAR) = srcSC_(isc_ZVAR) - dt_*(z2rhodot-rhodot*SC_(isc_ZMIX)**2)*dens
          end if
          
          ! Enthalpy
          if (isc_ENTH.gt.0) then
             ! Mass transfer
             srcSC_(isc_ENTH) = srcSC_(isc_ENTH) + dt_*enthdot*dens


             ! Soot radiation
             srcSC_(isc_ENTH) = srcSC_(isc_ENTH) - dt_*sootrad*mom(2)*MolarMassSoot/SootDensity*dens

          end if

          ! PAH
          if (use_PAH) then
             ! Rescale with transported value
             srcSC_(isc_PAH) = srcSC_(isc_PAH) + dt_*rhodot*min((SC_(isc_PAH)/ypah)**2,1.0e4_WP) ! M coded it to be 1.0e4_WP
          end if
       case ('finite chem')
          ! PAH equations
          call soot_finitechem_pahprodrates(SC_(isc_sc:isc_sc+N_nons-1),SC_(isc_T),dens,srcSC_)

          ! Soot radiation
          if (finitechem_sootrad) call soot_finitechem_sootrad(SC_(isc_sc:isc_sc+N_nons-1),SC_(isc_T),srcSC_(isc_T),mom(2))
          ! Temperature equation
          srcSC_(isc_T) = srcSC_(isc_T) + dt_*SC_(isc_T)*rhodot !SD
    end select

    return
  end subroutine soot_hmom_source_gas
  

  ! Compute the source terms for moments
  ! ------------------------------------
  subroutine soot_hmom_src_mom
    use soot
    implicit none

    return
  end subroutine soot_hmom_src_mom

end module soot_hmom


! ===================================== !
! INITIALIZE the hmom submodule of soot !
! ===================================== !
subroutine soot_hmom_init
  use soot_hmom
  use finitechem
  use masks
  implicit none

  real(WP) :: tol
  integer :: i,j,k,n
  logical :: soot_reset
  

  ! Initialize the finite rate chemistry portion of soot, if applicable
  select case(trim(chemistry))
     case('finite chem')
        call soot_finitechem_init
        !$OMP PARALLEL
        allocate(MassFracs(N_tot))
        !$OMP END PARALLEL
  end select

  ! Define some parameters
  CarbonToDiam = (6.0_WP*MolarMassSoot/(Pi*SootDensity*Avogadro))**(1.0_WP/3.0_WP)
  Cfm = Avogadro * (8.0_WP*Pi*R_uni/MolarMassSoot)**(0.5_WP) &
       * (0.75_WP*MolarMassSoot/(Avogadro*Pi*SootDensity))**(2.0_WP/3.0_WP)
  ! Need square root temperature for free-molecular regime constant
  Ccn = 8.0_WP * R_uni / 3.0_WP
  ! Need temperature / viscosity for continuum regime constant
  lambda = 3.0_WP * sqrt(Pi/(8.0_WP*R_uni)) / (6.0_WP*MolarMassSoot/(Pi*sootDensity*Avogadro))**(1.0_WP/3.0_WP)
  ! Need viscosity * sqrt molar mass / (density * sqrt tempeature) for slip factor constant
  av = 1.0_WP - (2.0_WP / Df)
  as = (3.0_WP / Df) - 1.0_WP
  FitC = 2.0_WP/3.0_WP
  FitExp = -0.2043_WP

  ! Allocate arrays for HMOM
  allocate(frac(dim))
  !$OMP PARALLEL
  allocate(kmom(dim))
  allocate(    mom(nEquations))
  allocate(old_mom(nEquations))
  allocate(src_mom(nEquations))
  !$OMP END PARALLEL

  ! Define the momenst used for HMOM
  call soot_define_moments

  ! Set dimer size
  DIMER_NBRC = 20.0_WP

  ! Set nucleated particle size
  NUCL_NBRC = 2.0_WP*DIMER_NBRC
  NUCL_SURF = NUCL_NBRC**(2.0_WP/3.0_WP)

  call parser_read('Soot reset', soot_reset,.false.)
  if (soot_reset) then
     SC(:,:,:,isc_mom(1)) = SMALLWEIGHT
     SC(:,:,:,isc_mom(nMoments)) = SMALLWEIGHT
     do n=2,nMoments-1
        SC(:,:,:,isc_mom(n)) = SMALLWEIGHT * NUCL_NBRC**moments(n,1) * NUCL_SURF**moments(n,2)
     end do
  end if

  ! If walls, put some typical values
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (mask(i,j).eq.1) then
              SC(i,j,k,isc_mom(1)) = SMALLWEIGHT
              SC(i,j,k,isc_mom(nMoments)) = SMALLWEIGHT
              do n=2,nMoments-1
                 SC(i,j,k,isc_mom(n)) = SMALLWEIGHT * NUCL_NBRC**moments(n,1) * NUCL_SURF**moments(n,2)
              end do
              if (use_soot_sgs) SC(i,j,k,isc_mom(nMoments+1)) = SMALLWEIGHT**2
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine soot_hmom_init


! ============================================ !
! Compute the source terms for the the scalars !
! - Chemtable model                            !
! ============================================ !
subroutine soot_hmom_source_scalar(SC_,srcSC_,srcP_)
  use soot_hmom
  use metric_generic
  use memory
  use finitechem
  use masks
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(in) :: SC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: srcP_

  integer :: i,j,k,n
  real(WP) :: y_pah

  real(WP) :: t_start, dt_soot, rate
  integer :: step

  real(WP) :: Wmix


  ! If not soot, then return
  if (.not.use_soot) return

  ! Start the timer
  call timing_start('soot')

  select case(trim(chemistry))
      case ('chemtable')
         ! Global chemtable lookups
         ! -- Experience says that is is always faster than local lookups and is more flexible
         call chemtable_lookup('sqrtT',tmp1)
         call chemtable_lookup('T_MU',tmp2)
         call chemtable_lookup('MUsqrtW_RHOsqrtT',tmp3)
         call chemtable_lookup('wCoeff',tmp4)
         call chemtable_lookup('oxCoeff',tmp5)
         call chemtable_lookup('o2Coeff',tmp6)
         call chemtable_lookup('ProdRate',tmp7)
         
         if (use_pah) then
            call chemtable_lookup('Y_PAH',tmp8)

            !$OMP PARALLEL DO
            do k=kmin_,kmax_
               do j=jmin_,jmax_
                  do i=imin_,imax_
                     tmp7(i,j,k) = max(tmp7(i,j,k),0.0_WP) * min((SC_(i,j,k,isc_PAH) / max(tmp8(i,j,k),1.0e-60))**2,1.0e4_WP)
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
         end if

         call chemtable_lookup('RHODOT',tmp9)
         call chemtable_lookup('SRC_ZMIX',tmp10)  ! SD
         call chemtable_lookup('ZSRC_ZMIX',tmp11)
         call chemtable_lookup('Z2RHODOT',tmp12)
         if (isc_ENTH.gt.0) then
            call chemtable_lookup('ENTHDOT',tmp13)
            call chemtable_lookup('SOOT_RAD',tmp14)
         end if
      case ('finite chem')
  end select
  
  !$OMP PARALLEL DO PRIVATE (t_start,rate,dt_soot,step,n)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
              ! Set the local environment
              dens = RHO(i,j,k)
              totdivg = sum(div_u(i,j,:)*Ut(i-st1:i+st2,j,k,isc_mom(1))) + &
                   sum(div_v(i,j,:)*Vt(i,j-st1:j+st2,k,isc_mom(1))) + &
                   sum(div_w(i,j,:)*Wt(i,j,k-st1:k+st2,isc_mom(1)))

              select case(trim(chemistry))
                 case ('chemtable')
                    sqrtT = tmp1(i,j,k)
                    T_MU = tmp2(i,j,k)
                    MUsqrtW_RHOsqrtT = tmp3(i,j,k)
                    wCoeff = tmp4(i,j,k)
                    oxCoeff = tmp5(i,j,k)
                    o2Coeff = tmp6(i,j,k)
                    ProdRate = tmp7(i,j,k)
                 case ('finite chem')
                    sqrtT = sqrt(T(i,j,k))
                    T_MU = T(i,j,k) / VISC(i,j,k)
                    call finitechem_W(SC_(i,j,k,isc_sc:isc_sc+N_nons-1),Wmix)
                    Wmix = Wmix / 1000.0_WP
                    MUsqrtW_RHOsqrtT = VISC(i,j,k)/RHO(i,j,k)*sqrt(Wmix/T(i,j,k))
                    call soot_finitechem_ksg(SC_(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,wCoeff)
                    call soot_finitechem_kox(SC_(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,oxCoeff)
                    call soot_finitechem_ko2(SC_(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,o2Coeff)
                    call soot_finitechem_dimerprodrate(SC_(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,ProdRate)
              end select

              ! Get the values
              do n=1,nEquations
                 mom(n) = SC_(i,j,k,isc_mom(n))
              end do
              call soot_hmom_clip
              do n=1,nEquations
                 mom(n) = mom(n)*dens
                 old_mom(n) = mom(n)
              end do


              t_start = 0.0_WP
              step = 0
              
              do while (t_start.lt.dt_ .and. step.lt.10000)
                 ! Dispersed-phase source terms
                 call soot_hmom_source_solid

                 ! Compute timestep
                 dt_soot = dt_ - t_start
                 !rate = 1.05_WP*maxval(-dt_soot*src_mom/mom) !SD default:1.05
                 rate = 1.05_WP*maxval(abs(-dt_soot*src_mom/mom))
                 if (rate.gt.1.0_WP) then
                    dt_soot = dt_soot / rate
                 end if
!!$              dt_soot = dt_ / 100.0_WP

                 ! Explicit integration
                 mom = mom + dt_soot*src_mom

                 ! Clip
                 mom = mom / dens
                 call soot_hmom_clip
                 mom = mom * dens

                 ! Increment
                 t_start = t_start+dt_soot
                 step = step+1
              end do

!!$              if (mom(1)*Avogadro*1.0e-6_WP.gt.2e+11_WP) then
!!$                 print*,'STEPs: ', step
!!$                 print*,i,j
!!$              end if


              ! If not done, add remaining source term
              if (t_start.lt.dt) then
                 call soot_hmom_source_solid
                 mom = mom + (dt_-t_start)*src_mom

                 mom = mom / dens
                 call soot_hmom_clip
                 mom = mom * dens
              end if

              ! Compute source term
              do n=1,nEquations
                 srcSC_(i,j,k,isc_mom(n)) = srcSC_(i,j,k,isc_mom(n)) + (mom(n)-old_mom(n))
              end do

              select case(trim(chemistry))
                 case ('chemtable')
                    if (use_pah) ypah = tmp8(i,j,k)
                    rhodot = min(tmp9(i,j,k),0.0_WP)
                    zsrc = min(tmp10(i,j,k),0.0_WP)
                    zzsrc = min(tmp11(i,j,k),0.0_WP)
                    z2rhodot = min(tmp12(i,j,k),0.0_WP)
                    if (isc_ENTH.gt.0) then
                       enthdot = min(tmp13(i,j,k),0.0_WP)
                       sootrad = max(tmp14(i,j,k),0.0_WP)
                    end if
                 case ('finite chem')
                    call soot_finitechem_rhodot(rhodot,SC_(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens)
                    ! SD
                    !if (x(i).eq.2e-2_WP .and. j.eq.3) print*,'rhodot',i,j,rhodot,'A2',SC_(i,j,k,37-isc_sc+1),'A2R5',SC_(i,j,k,39-isc_sc+1)
              end select

              ! Evaluate moments at midstep for gas-phase source terms
              ! -- MEM: This is not the midstep.  What should we do here?
              do n=1,nEquations
                 mom(n) = 0.5_WP*(mom(n)+old_mom(n))
              end do

              ! Gas-phase source terms
              call soot_hmom_source_gas(SC_(i,j,k,:),srcSC_(i,j,k,:),srcP_(i,j,k))
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Stop the timer
  call timing_stop('soot')

  return
end subroutine soot_hmom_source_scalar


! ======================== !
! Clip the soot quantities !
! ======================== !
subroutine soot_hmom_clip
  use soot_hmom
  implicit none

  integer :: n

  ! Check for globally small moments
  if ( mom(1).lt.SMALLWEIGHT .or. &
       mom(2).lt.SMALLWEIGHT*NUCL_NBRC .or. &
       mom(3).lt.SMALLWEIGHT*NUCL_SURF ) then
     mom(2) = max(mom(2), SMALLWEIGHT*NUCL_NBRC)
     mom(1) = mom(2) / NUCL_NBRC
     mom(3) = NUCL_SURF * mom(1)
     mom(nMoments) = mom(1)
  end if

  ! Check for size of second mode
  if (mom(2).lt.NUCL_NBRC*mom(1) .or. &
      mom(3).lt.NUCL_SURF*mom(1)) then
     mom(1) = mom(2) / NUCL_NBRC
     mom(3) = mom(1) * NUCL_SURF
  end if

  ! Check for (co)variance of second mode
  if (nMoments.ge.7) then
     mom(4) = max(mom(4),mom(2)*mom(2)/mom(1))
     mom(5) = max(mom(5),mom(2)*mom(3)/mom(1))
     mom(6) = max(mom(6),mom(3)*mom(3)/mom(1))
  end if

  ! Check for small weight of first mode
  if (mom(nMoments).lt.SMALLWEIGHT) then
     do n=1,nMoments-1
        mom(n) = mom(n) + (SMALLWEIGHT-mom(nMoments))*NUCL_NBRC**moments(n,1)*NUCL_SURF**moments(n,2)
     end do
     mom(nMoments) = SMALLWEIGHT
  end if
  if (mom(nMoments).gt.mom(1)) then
     mom(nMoments) = mom(1)
  end if

  ! Clip the number density squared
  if (use_soot_sgs) then
     if (mom(nMoments+1)*dens.lt.(mom(1)*dens)**2) then
        mom(nMoments+1) = mom(1)**2*dens
     end if
  end if

  return
end subroutine soot_hmom_clip


! ============================================== !
! Routine to execute after the transport of soot !
! -> Clean up the variables (clipping,..)        !
! -> Compute some mean statistical quantities    !
! ============================================== !
subroutine soot_hmom_poststep
  use soot_hmom
  use finitechem
  use masks
  use memory
  implicit none

  integer :: i,j,k,n
  real(WP) :: omega,src_tmp,Wmix
  
  select case(trim(chemistry))
  case ('chemtable')
     ! Compute quantities which are only functions of the combustion scalars
     call chemtable_lookup('sqrtT',tmp1)
     call chemtable_lookup('T_MU',tmp2)
     call chemtable_lookup('MUsqrtW_RHOsqrtT',tmp3)
     call chemtable_lookup('wCoeff',tmp4)
     call chemtable_lookup('oxCoeff',tmp5)
     call chemtable_lookup('o2Coeff',tmp6)
     call chemtable_lookup('ProdRate',tmp7)
     if (use_pah) then
        call chemtable_lookup('Y_PAH',tmp8)
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 tmp7(i,j,k) = tmp7(i,j,k) * min((SC(i,j,k,isc_PAH) / tmp8(i,j,k))**2,1.0e4_WP)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case ('finite chem')
  end select
  
  !$OMP PARALLEL DO PRIVATE(omega,src_tmp,Wmix)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_

           ! If walls put some typical values
           if (mask(i,j).eq.1) then
              numdens(i,j,k)  = SMALLWEIGHT*RHO(i,j,k) * Avogadro*1.0E-6_WP
              volfrac(i,j,k)  = SMALLWEIGHT*NUCL_NBRC*RHO(i,j,k) * MolarMassSoot/SootDensity
              partdiam(i,j,k) = NUCL_NBRC/NUCL_SURF * CarbonToDiam*1.0E9_WP
              partaggr(i,j,k) = 1.0_WP
              intermit(i,j,k) = 1.0_WP
              Nsrc_nucl(i,j,k) = 0.0_WP
              Nsrc_coag(i,j,k) = 0.0_WP
              Nsrc_ox  (i,j,k) = 0.0_WP
              Nsrc_frag(i,j,k) = 0.0_WP
              FVsrc_nucl(i,j,k) = 0.0_WP
              FVsrc_cond(i,j,k) = 0.0_WP
              FVsrc_sg  (i,j,k) = 0.0_WP
              FVsrc_ox  (i,j,k) = 0.0_WP
           else

              ! Get the values
              do n=1,nEquations !SD
                 if (SC(i,j,k,isc_ZMIX).gt.soot_z) then
                    mom(n) = SC(i,j,k,isc_mom(n))
                 else
                    mom(n) = 0.0_WP
                 end if
              end do
              call soot_hmom_clip
              do n=1,nEquations
                 SC(i,j,k,isc_mom(n)) = mom(n)
                 mom(n) = mom(n) * RHO(i,j,k)
              end do
              
              ! Compute the subfilter intermittency and normalize the moments
              omega = 0.0_WP
              if (use_soot_sgs) omega = min(max(1.0_WP - (mom(1)**2/mom(nMoments+1)),0.0_WP),1.0_WP-1e-8_WP) ! SD 1.0_WP - (mom(1)**2/mom(nMoments+1))
              
              mom = mom / (1.0_WP - omega)

              ! Number density
              numdens(i,j,k) = (1.0_WP-omega) * soot_fracmom(0.0_WP,0.0_WP) * Avogadro*1.0e-6_WP

              ! Volume Fraction
              volfrac(i,j,k) = (1.0_WP-omega) * soot_fracmom(1.0_WP,0.0_WP) * MolarMassSoot/SootDensity
              
              ! Primary particle diameter
              partdiam(i,j,k) = (1.0_WP-omega) * soot_fracmom(1.0_WP,-1.0_WP) / soot_fracmom(0.0_WP,0.0_WP) * CarbonToDiam*1.0E9_WP
              
              ! Number of primary particles per aggregate
              partaggr(i,j,k) = (1.0_WP-omega) * soot_fracmom(-2.0_WP,3.0_WP) / soot_fracmom(0.0_WP,0.0_WP)

              ! Intermittency
              if (volfrac(i,j,k)/(1.0_WP-omega) .lt. soot_int) then
                 intermit(i,j,k) = 1.0_WP
              else
                 intermit(i,j,k) = omega
              end if

              ! Get the density
              dens = RHO(i,j,k)

              select case(trim(chemistry))
              case ('chemtable')
                 ! Compute quantities which are only functions of the combustion scalars
                 sqrtT = tmp1(i,j,k)
                 T_MU = tmp2(i,j,k)
                 MUsqrtW_RHOsqrtT = tmp3(i,j,k)
                 wCoeff = tmp4(i,j,k)
                 oxCoeff = tmp5(i,j,k)
                 o2Coeff = tmp6(i,j,k)
                 prodRate = tmp7(i,j,k)
              case ('finite chem')
                 sqrtT = sqrt(T(i,j,k))
                 T_MU = T(i,j,k) / VISC(i,j,k)
                 call finitechem_W(SC(i,j,k,isc_sc:isc_sc+N_nons-1),Wmix)
                 Wmix = Wmix / 1000.0_WP
                 MUsqrtW_RHOsqrtT = VISC(i,j,k)/RHO(i,j,k)*sqrt(Wmix/T(i,j,k))
                 call soot_finitechem_ksg(SC(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,wCoeff)
                 call soot_finitechem_kox(SC(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,oxCoeff)
                 call soot_finitechem_ko2(SC(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,o2Coeff)
                 call soot_finitechem_dimerprodrate(SC(i,j,k,isc_sc:isc_sc+N_nons-1),T(i,j,k),dens,ProdRate)
              end select


              ! Rescale dimer production rate for transported 
              !prodRate = local_nbrC*prodRate / DIMER_NBRC
     
              ! Prepare for compute source terms
              call soot_pah_molecules

              ! Number density source terms
              kmom = moments(1,:)
              call soot_nucleation(kmom,src_tmp)
              Nsrc_nucl(i,j,k) = src_tmp * Avogadro*1.0e-6_WP
              call soot_coagulation(kmom,src_tmp)
              Nsrc_coag(i,j,k) = (1.0_WP-omega)*src_tmp * Avogadro*1.0e-6_WP
              call soot_surfaceoxidation(kmom,src_tmp)
              Nsrc_ox  (i,j,k) = (1.0_WP-omega)*src_tmp * Avogadro*1.0e-6_WP
              call soot_fragmentation(kmom,src_tmp)
              Nsrc_frag(i,j,k) = (1.0_WP-omega)*src_tmp * Avogadro*1.0e-6_WP

              ! Volume fraction source terms
              kmom = moments(2,:)
              call soot_nucleation(kmom,src_tmp)
              FVsrc_nucl(i,j,k) = src_tmp * MolarMassSoot/SootDensity
              call soot_condensation(kmom,src_tmp)
              FVsrc_cond(i,j,k) = (1.0_WP-omega)*src_tmp * MolarMassSoot/SootDensity
              call soot_surfacereaction(kmom,src_tmp)
              FVsrc_sg(i,j,k) = (1.0_WP-omega)*src_tmp * MolarMassSoot/SootDensity
              call soot_surfaceoxidation(kmom,src_tmp)
              FVsrc_ox(i,j,k) = (1.0_WP-omega)*src_tmp * MolarMassSoot/SootDensity
              ! SD
              !if (i.eq.261 .and. j.eq.3) print*,'SOOT_HMOM_POSTSTEP','FVsrc_PAH',FVsrc_nucl(i,j,3)+FVsrc_cond(261,3,3),'FVsrc_SG',FVsrc_sg(i,j,3),'FVsrc_OX', FVsrc_ox(i,j,3),'fv',volfrac(i,j,3) ! SD
              ! Unnormalize the moments
              mom = mom * (1.0_WP-omega)
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine soot_hmom_poststep
