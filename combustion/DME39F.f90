!----------------------------------------------------------
! ======= DME39F.f90 =======
!----------------------------------------------------------
subroutine PRODRATES( CDOT, W, K, C, M, TEMP, PRESSURE)
!----------------------------------------------------------
!     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT
!     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),
!     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE
!     JUST WORK SPACE FOR THIS FUNCTION.
!     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN
!     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE 
!     CONCENTRATIONS, WHICH ARE COMPUTED IN THIS FUNCTION.
!     TEMP IS THE TEMPERATURE IN [K] AND
!     PRESSURE IS THE PRESSURE IN [PA].
!     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',
!     'CTCHZERO'
!----------------------------------------------------------
!!$CH2GSGXCH2 => CH2X
!!$CH3OCH2O2  => RO2
!!$CH3OCH2OH  => ROH
!!$CH2OCH2O2H => QOOH
!!$O2CH2OCH2O2H => O2QOOH
!!$HO2CH2OCHO => HO2QHO

  implicit none
      include 'DME39F90.h'
      real(DP) :: CDOT(39), W(350), K(350), &
      C(39), M(20), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM1) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM1) = M(MM1) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM2) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM2) = M(MM2) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM3) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM3) = M(MM3) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM4) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM4) = M(MM4) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM5) = C(SN2) + C(SH) &
	    + 0.78 * C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 11 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM5) = M(MM5) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM6) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM6) = M(MM6) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM7) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM7) = M(MM7) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM8) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM8) = M(MM8) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM9) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM9) = M(MM9) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM10) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.5 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.9 * C(SCO) &
	    + 3.8 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM10) = M(MM10) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM11) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SH2) &
	    + 5 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 2 * C(SCO) &
	    + 3 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM11) = M(MM11) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM12) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM12) = M(MM12) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM0) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SH2) &
	    + C(SH2O) + C(SHO2) &
	    + C(SH2O2) + C(SCO) &
	    + C(SCO2) + C(SHCO) &
	    + C(SCH3) + C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM0) = M(MM0) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM13) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM13) = M(MM13) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM14) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM14) = M(MM14) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM15) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM15) = M(MM15) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM16) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM16) = M(MM16) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM17) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM17) = M(MM17) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM18) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + C(SHO2) &
	    + C(SH2O2) + 1.5 * C(SCO) &
	    + 2 * C(SCO2) + C(SHCO) &
	    + C(SCH3) + 2 * C(SCH4) &
	    + C(SCH2O) + C(SCH3O) &
	    + 3 * C(SC2H6) + C(SCH2OH) &
	    + C(SC2H5) + C(SCH2) &
	    + C(SCH2X) + C(SC2H4) &
	    + C(SCH3HCO) + C(SC2H2) &
	    + C(SC2H3) + C(SCH3OCH3) &
	    + C(SCH3OCH2) + C(SCH3OCH2O) &
	    + C(SCH3OCHO) + C(SCH3OCO)
      M(MM18) = M(MM18) + C(SRO2) + C(SROH) &
	    + C(SQOOH) + C(SO2QOOH) &
	    + C(SHO2QHO) + C(SOCH2OCHO) &
	    + C(SHOCH2OCO) + C(SHOCH2O) &
	    + C(SHCOOH)
      M(MM19) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SHCO) + C(SCH3) &
	    + C(SCH4) + C(SCH2O) &
	    + C(SCH3O) + C(SC2H6) &
	    + C(SCH2OH) + C(SC2H5) &
	    + C(SCH2) + C(SCH2X) &
	    + C(SC2H4) + C(SCH3HCO) &
	    + C(SC2H2) + C(SC2H3) &
	    + C(SCH3OCH3) + C(SCH3OCH2) &
	    + C(SCH3OCH2O) + C(SCH3OCHO) &
	    + C(SCH3OCO) + C(SRO2) &
	    + C(SROH) + C(SQOOH)
      M(MM19) = M(MM19) + C(SO2QOOH) + C(SHO2QHO) &
	    + C(SOCH2OCHO) + C(SHOCH2OCO) &
	    + C(SHOCH2O) + C(SHCOOH)


      K(R1F) = 3.5470000000D+12 &
	   * exp(-0.406 * LT - 69450000 / RT)
      K(R1B) = 7.0530202531D+09 &
	   * exp(0.0289576 * LT + 1111131.16 / RT)
      K(R2F) = 5.0800000000D+01 &
	   * exp(2.67 * LT - 26317000 / RT)
      K(R2B) = 2.9795859804D+01 &
	   * exp(2.63497 * LT - 20550909.81 / RT)
      K(R3F) = 2.1600000000D+05 &
	   * exp(1.51 * LT - 14351000 / RT)
      K(R3B) = 2.2182689717D+06 &
	   * exp(1.41638 * LT - 77059304.43 / RT)
      K(R4F) = 2.9700000000D+03 &
	   * exp(2.02 * LT - 56066000 / RT)
      K(R4B) = 1.6962435366D+02 &
	   * exp(2.07859 * LT + 12408394.61 / RT)
      K(R5F) = 4.5770000000D+16 &
	   * exp(-1.4 * LT - 436726000 / RT)
      K(R5B) = 1.9619085286D+14 &
	   * exp(-1.7488 * LT - 3652636.474 / RT)
      K(R6F) = 6.1650000000D+09 * exp(-0.5 * LT)
      K(R6B) = 4.2424184779D+14 &
	   * exp(-0.621192 * LT - 497868404.5 / RT)
      K(R7F) = 4.7140000000D+12 * exp(-1 * LT)
      K(R7B) = 6.4503598776D+14 &
	   * exp(-0.686234 * LT - 427307273.3 / RT)
      K(R8F) = 3.8000000000D+16 * exp(-2 * LT)
      K(R8B) = 9.1042926083D+19 &
	   * exp(-1.74482 * LT - 495781668 / RT)
      K0TROE = 6.3660000000D+14 &
	   * exp(-1.72 * LT - 2196000 / RT)
      KINFTROE = 1.4750000000D+09 * exp(0.6 * LT)
      FCTROE = 0.2 * EXP( -TEMP / 1e-30 ) &
	   + 0.8 * EXP( -TEMP / 1e+30 )
      K(R9F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K0TROE = 9.0293444696D+17 &
	   * exp(-1.73291 * LT - 206859348.7 / RT)
      KINFTROE = 2.0920959932D+12 &
	   * exp(0.587093 * LT - 204663348.7 / RT)
      K(R9B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K(R10F) = 1.6600000000D+10 * exp(-3443000 / RT)
      K(R10B) = 2.7303646846D+09 &
	   * exp(0.361707 * LT - 231853014.8 / RT)
      K(R11F) = 7.0790000000D+10 * exp(-1234000 / RT)
      K(R11B) = 1.3579714350D+07 &
	   * exp(0.761631 * LT - 153316793.4 / RT)
      K(R12F) = 3.2500000000D+10
      K(R12B) = 3.1353652512D+09 &
	   * exp(0.326673 * LT - 222643924.6 / RT)
      K(R13F) = 2.8900000000D+10 * exp(2079000 / RT)
      K(R13B) = 4.8816975193D+10 &
	   * exp(0.268083 * LT - 289039319.2 / RT)
      K(R14F) = 4.2000000000D+11 * exp(-50133000 / RT)
      K(R14B) = 6.3574616360D+13 &
	   * exp(-0.389273 * LT - 212208672.1 / RT)
      K(R15F) = 1.3000000000D+08 * exp(6817000 / RT)
      K(R15B) = 1.9677857445D+10 &
	   * exp(-0.389273 * LT - 155258672.1 / RT)
      K0TROE = 1.2020000000D+14 * exp(-190372000 / RT)
      KINFTROE = 2.9510000000D+14 * exp(-202631000 / RT)
      FCTROE = 0.5 * EXP( -TEMP / 1e-30 ) &
	   + 0.5 * EXP( -TEMP / 1e+30 )
      K(R16F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM6) )
      K0TROE = 1.0739872977D+05 &
	   * exp(1.16381 * LT + 24284227.4 / RT)
      KINFTROE = 2.6367192308D+05 &
	   * exp(1.16381 * LT + 12025227.4 / RT)
      K(R16B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM6) )
      K(R17F) = 2.4100000000D+10 * exp(-16610000 / RT)
      K(R17B) = 5.1591045676D+04 &
	   * exp(1.41899 * LT - 297735440.6 / RT)
      K(R18F) = 4.8200000000D+10 * exp(-33263000 / RT)
      K(R18B) = 5.2375134409D+07 &
	   * exp(0.75098 * LT - 99597342.69 / RT)
      K(R19F) = 9.5500000000D+03 &
	   * exp(2 * LT - 16610000 / RT)
      K(R19B) = 6.0865850326D+00 &
	   * exp(2.71595 * LT - 77178252.5 / RT)
      K(R20F) = 1.0000000000D+09
      K(R20B) = 1.1159341435D+07 &
	   * exp(0.657356 * LT - 129042647.1 / RT)
      K(R21F) = 5.8000000000D+11 * exp(-39986000 / RT)
      K(R21B) = 6.4724180323D+09 &
	   * exp(0.657356 * LT - 169028647.1 / RT)
      K0TROE = 1.5500000000D+18 &
	   * exp(-2.79 * LT - 17535000 / RT)
      KINFTROE = 1.8000000000D+07 * exp(-9975000 / RT)
      FCTROE = 1 * EXP( -0 / TEMP )
      K(R22F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K0TROE = 1.1858122417D+27 &
	   * exp(-3.75543 * LT - 552312028.5 / RT)
      KINFTROE = 1.3770722807D+16 &
	   * exp(-0.965431 * LT - 544752028.5 / RT)
      K(R22B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K(R23F) = 2.5300000000D+09 * exp(-199577000 / RT)
      K(R23B) = 2.8127059229D+13 &
	   * exp(-0.844239 * LT - 236485624 / RT)
      K(R24F) = 3.0100000000D+10 * exp(-96232000 / RT)
      K(R24B) = 3.2283088892D+13 &
	   * exp(-0.517566 * LT - 355784548.6 / RT)
      K(R25F) = 2.2290000000D+02 &
	   * exp(1.89 * LT + 4848000 / RT)
      K(R25B) = 1.2462349669D+09 &
	   * exp(0.610803 * LT - 102621755.2 / RT)
      K(R26F) = 4.7480000000D+08 &
	   * exp(0.659 * LT - 62233000 / RT)
      K(R26B) = 9.5208413773D+04 &
	   * exp(0.891462 * LT + 1964560.922 / RT)
      K(R27F) = 7.5800000000D+09 * exp(-1715000 / RT)
      K(R27B) = 2.1558741858D+09 &
	   * exp(0.219554 * LT - 142180787.8 / RT)
      K(R28F) = 7.2300000000D+10
      K(R28B) = 3.3822450822D+09 &
	   * exp(0.581262 * LT - 368875802.6 / RT)
      K(R29F) = 3.0200000000D+10
      K(R29B) = 8.2864006839D+08 &
	   * exp(0.546227 * LT - 363109712.4 / RT)
      K(R30F) = 3.0200000000D+10
      K(R30B) = 1.4508889496D+10 &
	   * exp(0.487638 * LT - 431584107 / RT)
      K(R31F) = 3.0000000000D+10
      K(R31B) = 4.6022489291D+15 &
	   * exp(-0.73297 * LT - 470579467.6 / RT)
      K(R32F) = 3.0000000000D+10
      K(R32B) = 6.4520015108D+09 &
	   * exp(-0.285105 * LT - 195354987.7 / RT)
      K(R33F) = 3.0000000000D+09
      K(R33B) = 2.8141850446D+04 &
	   * exp(0.813723 * LT - 304678241.7 / RT)
      K(R34F) = 2.6500000000D+10
      K(R34B) = 4.2149530019D+11 &
	   * exp(0.291083 * LT - 376812980.9 / RT)
      K(R35F) = 3.0000000000D+10
      K(R35B) = 9.8375045521D+11 &
	   * exp(0.111904 * LT - 314381994.1 / RT)
      K(R36F) = 3.3000000000D+36 &
	   * exp(-6.3 * LT - 417982000 / RT)
      K(R36B) = 2.0179707531D+31 &
	   * exp(-6.17944 * LT - 39402444.96 / RT)
      K(R37F) = 3.1000000000D+42 &
	   * exp(-8 * LT - 407982000 / RT)
      K(R37B) = 8.8680758342D+35 &
	   * exp(-7.29818 * LT - 398278247.6 / RT)
      K(R38F) = 5.7400000000D+04 &
	   * exp(1.9 * LT - 11500000 / RT)
      K(R38B) = 8.1887003172D+01 &
	   * exp(2.36936 * LT - 65993808.48 / RT)
      K(R39F) = 1.8100000000D+10 * exp(-12887000 / RT)
      K(R39B) = 1.5145160141D+07 &
	   * exp(0.434323 * LT - 61614718.3 / RT)
      K(R40F) = 3.4300000000D+06 &
	   * exp(1.18 * LT + 1870000 / RT)
      K(R40B) = 5.0252498407D+04 &
	   * exp(1.55573 * LT - 115332112.9 / RT)
      K(R41F) = 1.2300000000D+03 &
	   * exp(3 * LT - 217568000 / RT)
      K(R41B) = 1.0668309991D+01 &
	   * exp(3.10765 * LT - 43651793.69 / RT)
      K(R42F) = 4.1100000000D+01 &
	   * exp(2.5 * LT - 42719000 / RT)
      K(R42B) = 5.3959360590D+01 &
	   * exp(2.21838 * LT - 30878465.8 / RT)
      K(R43F) = 3.6360000000D-09 &
	   * exp(5.42 * LT - 4176000 / RT)
      K(R43B) = 1.7636282446D-09 &
	   * exp(5.59918 * LT - 66606986.74 / RT)
      K(R44F) = 8.4300000000D+10
      K(R44B) = 1.3513165132D+13 &
	   * exp(-0.304279 * LT - 296214081.9 / RT)
      K(R45F) = 1.9900000000D+15 &
	   * exp(-1.57 * LT - 122298000 / RT)
      K(R45B) = 3.5389559985D+17 &
	   * exp(-2.05404 * LT - 7423963.149 / RT)
      K(R46F) = 3.7400000000D+08 * exp(-61254000 / RT)
      K(R46B) = 1.1921064944D+08 &
	   * exp(0.130678 * LT - 286906950.7 / RT)
      K(R47F) = 2.4100000000D+07 &
	   * exp(0.76 * LT + 9728000 / RT)
      K(R47B) = 4.1346990936D+08 &
	   * exp(0.602631 * LT - 98041887.75 / RT)
      K0TROE = 8.0540000000D+25 &
	   * exp(-3.75 * LT - 4107000 / RT)
      KINFTROE = 2.2770000000D+12 &
	   * exp(-0.69 * LT - 732000 / RT)
      FCTROE = 1 * EXP( -TEMP / 570 ) &
	   + 1 * EXP( -1e+30 / TEMP )
      K(R48F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM11) )
      K0TROE = 1.2503013086D+37 &
	   * exp(-5.10728 * LT - 388336266 / RT)
      KINFTROE = 3.5348101313D+23 &
	   * exp(-2.04728 * LT - 384961266 / RT)
      K(R48B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM11) )
      K0TROE = 2.4770000000D+27 &
	   * exp(-4.76 * LT - 10209000 / RT)
      KINFTROE = 1.2700000000D+13 &
	   * exp(-0.63 * LT - 1602000 / RT)
      FCTROE = 0.217 * EXP( -TEMP / 74 ) &
	   + 0.783 * EXP( -TEMP / 2941 ) &
	   + 1 * EXP( -6964 / TEMP )
      K(R49F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM12) )
      K0TROE = 1.9647543084D+32 &
	   * exp(-4.70138 * LT - 451219541.8 / RT)
      KINFTROE = 1.0073629276D+18 &
	   * exp(-0.571379 * LT - 442612541.8 / RT)
      K(R49B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM12) )
      K(R50F) = 5.4700000000D+04 &
	   * exp(1.97 * LT - 46903000 / RT)
      K(R50B) = 1.6088191997D+02 &
	   * exp(2.26018 * LT - 38965821.74 / RT)
      K(R51F) = 3.1500000000D+09 &
	   * exp(0.5 * LT - 43053000 / RT)
      K(R51B) = 5.4340380842D+06 &
	   * exp(0.755145 * LT - 29349731.56 / RT)
      K(R52F) = 5.7200000000D+03 &
	   * exp(1.96 * LT - 11042000 / RT)
      K(R52B) = 1.7277320766D+02 &
	   * exp(2.15656 * LT - 65813126.17 / RT)
      K(R53F) = 3.1600000000D+09
      K(R53B) = 1.7671757301D+11 &
	   * exp(0.0715283 * LT - 236347193 / RT)
      K(R54F) = 1.8100000000D+08 * exp(-77739000 / RT)
      K(R54B) = 4.8991459399D+08 &
	   * exp(-0.460801 * LT - 3467479.058 / RT)
      K(R55F) = 1.0000000000D+11 * exp(-105018000 / RT)
      K(R55B) = 1.6133092170D+07 &
	   * exp(0.3196 * LT + 13210812.05 / RT)
      K(R56F) = 6.0000000000D+09
      K(R56B) = 2.2582448200D+08 &
	   * exp(0.6684 * LT - 314844551.5 / RT)
      K(R57F) = 9.6350000000D+10
      K(R57B) = 1.3268877722D+07 &
	   * exp(0.937644 * LT - 12864379.42 / RT)
      K(R58F) = 4.2000000000D+10
      K(R58B) = 9.2717405987D+08 &
	   * exp(0.633365 * LT - 309078461.3 / RT)
      K(R59F) = 2.4000000000D+10
      K(R59B) = 9.2766563237D+09 &
	   * exp(0.574776 * LT - 377552855.9 / RT)
      K(R60F) = 2.4100000000D+11 * exp(-20991000 / RT)
      K(R60B) = 5.5147298808D+10 &
	   * exp(0.306692 * LT - 107425536.7 / RT)
      K(R61F) = 1.5100000000D+12 * exp(-1 * LT)
      K(R61B) = 3.4552871867D+11 &
	   * exp(-0.693308 * LT - 86434536.68 / RT)
      K(R62F) = 1.2000000000D+10
      K(R62B) = 4.1564533076D+11 &
	   * exp(-0.0825805 * LT - 248510208.8 / RT)
      K(R63F) = 1.5000000000D+10
      K(R63B) = 3.9573817470D+11 &
	   * exp(0.199042 * LT - 260350743 / RT)
      K(R64F) = 8.3000000000D+14 &
	   * exp(-1.2 * LT - 64852000 / RT)
      K(R64B) = 1.0871891022D+10 &
	   * exp(-0.899045 * LT + 21928285.78 / RT)
      K(R65F) = 3.2000000000D+10
      K(R65B) = 3.5780123422D+05 &
	   * exp(0.919 * LT - 44312905.69 / RT)
      K(R66F) = 6.0000000000D+09
      K(R66B) = 1.0754063973D+07 &
	   * exp(0.61472 * LT - 340526987.6 / RT)
      K(R67F) = 1.8000000000D+10
      K(R67B) = 5.6488769410D+08 &
	   * exp(0.556131 * LT - 409001382.2 / RT)
      K(R68F) = 9.0330000000D+10 * exp(-50124000 / RT)
      K(R68B) = 1.6782188533D+09 &
	   * exp(0.288048 * LT - 168007063 / RT)
      K(R69F) = 2.2000000000D+07 * exp(-7314000 / RT)
      K(R69B) = 4.0873258908D+05 &
	   * exp(0.288048 * LT - 125197063 / RT)
      K(R70F) = 3.0000000000D+08
      K(R70B) = 8.4366940079D+08 &
	   * exp(-0.101225 * LT - 279958735.1 / RT)
      K(R71F) = 1.6000000000D+10 * exp(-49371000 / RT)
      K(R71B) = 1.0002342066D+12 &
	   * exp(-0.360197 * LT - 201153660.8 / RT)
      K(R72F) = 4.9900000000D+09 &
	   * exp(0.1 * LT - 44350000 / RT)
      K(R72B) = 7.1770697790D+14 &
	   * exp(-1.06665 * LT - 7072459.094 / RT)
      K(R73F) = 2.4600000000D+03 &
	   * exp(2 * LT - 34602000 / RT)
      K(R73B) = 3.0002823388D+03 &
	   * exp(1.77911 * LT - 57823753.84 / RT)
      K(R74F) = 1.6000000000D+10 * exp(2385000 / RT)
      K(R74B) = 1.4227561628D+10 &
	   * exp(-0.283175 * LT - 58139741.82 / RT)
      K(R75F) = 5.6000000000D+04 &
	   * exp(1.6 * LT - 22677000 / RT)
      K(R75B) = 1.3868879542D+03 &
	   * exp(2.01744 * LT - 54226372.33 / RT)
      K(R76F) = 2.5010000000D+10
      K(R76B) = 8.4953948027D+08 &
	   * exp(0.47973 * LT + 5753615.644 / RT)
      K(R77F) = 4.0000000000D+10
      K(R77B) = 3.6126044532D+17 &
	   * exp(-1.23115 * LT - 275776402 / RT)
      K(R78F) = 1.2000000000D+10 * exp(2385000 / RT)
      K(R78B) = 7.9017837861D+16 &
	   * exp(-1.29344 * LT - 310694390 / RT)
      K(R79F) = 1.6000000000D+10
      K(R79B) = 6.0768947333D+03 &
	   * exp(1.39873 * LT - 38559290.05 / RT)
      K(R80F) = 1.1500000000D+05 &
	   * exp(1.9 * LT - 31506000 / RT)
      K(R80B) = 2.4856705124D+01 &
	   * exp(2.43943 * LT - 43072556.61 / RT)
      K(R81F) = 8.9800000000D+04 &
	   * exp(1.92 * LT - 23807000 / RT)
      K(R81B) = 1.1384507997D+01 &
	   * exp(2.4244 * LT - 29607466.43 / RT)
      K(R82F) = 3.5400000000D+03 &
	   * exp(2.12 * LT - 3640000 / RT)
      K(R82B) = 7.8579531528D+00 &
	   * exp(2.56581 * LT - 77914861.04 / RT)
      K(R83F) = 4.0000000000D+10 * exp(-212966000 / RT)
      K(R83B) = 5.2564573058D+07 &
	   * exp(0.177723 * LT + 3877458.177 / RT)
      K(R84F) = 2.9400000000D+08 * exp(-62509000 / RT)
      K(R84B) = 5.8481019909D+07 &
	   * exp(-0.21155 * LT - 7741213.928 / RT)
      K(R85F) = 6.1400000000D+03 &
	   * exp(1.74 * LT - 43723000 / RT)
      K(R85B) = 4.5122606321D+02 &
	   * exp(1.98925 * LT - 63226734.87 / RT)
      K0TROE = 1.9900000000D+35 &
	   * exp(-7.08 * LT - 27970000 / RT)
      KINFTROE = 5.2100000000D+14 &
	   * exp(-0.99 * LT - 6611000 / RT)
      FCTROE = 0.1578 * EXP( -TEMP / 125 ) &
	   + 0.8422 * EXP( -TEMP / 2219 ) &
	   + 1 * EXP( -6882 / TEMP )
      K(R86F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM13) )
      K0TROE = 2.1478775258D+41 &
	   * exp(-7.27063 * LT - 449476806.9 / RT)
      KINFTROE = 5.6233376430D+20 &
	   * exp(-1.18063 * LT - 428117806.9 / RT)
      K(R86B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM13) )
      K(R87F) = 2.0000000000D+09
      K(R87B) = 3.0285635419D+08 &
	   * exp(0.44657 * LT - 281895010.8 / RT)
      K(R88F) = 1.3200000000D+11
      K(R88B) = 1.4711495303D+08 &
	   * exp(0.862375 * LT - 333491622.8 / RT)
      K(R89F) = 2.0000000000D+07
      K(R89B) = 1.8412981636D+07 &
	   * exp(0.084863 * LT - 53484996.04 / RT)
      K(R90F) = 1.4000000000D+09
      K(R90B) = 9.8081931579D+11 &
	   * exp(-0.0928598 * LT - 270328454.2 / RT)
      K(R91F) = 1.2000000000D+11
      K(R91B) = 2.5971815126D+13 &
	   * exp(0.0418316 * LT - 357309246 / RT)
      K(R92F) = 8.0200000000D+10
      K(R92B) = 1.0134383975D+14 &
	   * exp(-0.536729 * LT - 317349590.2 / RT)
      K0TROE = 7.0000000000D+47 &
	   * exp(-9.31 * LT - 417814000 / RT)
      KINFTROE = 8.0000000000D+12 &
	   * exp(0.44 * LT - 371414000 / RT)
      FCTROE = 0.2655 * EXP( -TEMP / 180 ) &
	   + 0.7345 * EXP( -TEMP / 1035 ) &
	   + 1 * EXP( -5417 / TEMP )
      K(R93F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM14) )
      K0TROE = 1.3620567873D+42 &
	   * exp(-8.99745 * LT - 237600623.6 / RT)
      KINFTROE = 1.5566363283D+07 &
	   * exp(0.75255 * LT - 191200623.6 / RT)
      K(R93B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM14) )
      K0TROE = 1.2000000000D+36 &
	   * exp(-7.62 * LT - 29162000 / RT)
      KINFTROE = 1.0800000000D+09 &
	   * exp(0.454 * LT - 7615000 / RT)
      FCTROE = 0.0247 * EXP( -TEMP / 210 ) &
	   + 0.9753 * EXP( -TEMP / 984 ) &
	   + 1 * EXP( -4374 / TEMP )
      K(R94F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM15) )
      K0TROE = 1.8487437087D+39 &
	   * exp(-7.71777 * LT - 180340352.7 / RT)
      KINFTROE = 1.6638693378D+12 &
	   * exp(0.35623 * LT - 158793352.7 / RT)
      K(R94B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM15) )
      K0TROE = 1.4000000000D+24 &
	   * exp(-3.86 * LT - 13891000 / RT)
      KINFTROE = 6.0800000000D+09 &
	   * exp(0.27 * LT - 1172000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(R95F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM16) )
      K0TROE = 7.6776073818D+29 &
	   * exp(-3.99563 * LT - 480393379.5 / RT)
      KINFTROE = 3.3342752058D+15 &
	   * exp(0.134367 * LT - 467674379.5 / RT)
      K(R95B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM16) )
      K(R96F) = 1.3250000000D+03 &
	   * exp(2.53 * LT - 51212000 / RT)
      K(R96B) = 5.6366383607D-01 &
	   * exp(3.01443 * LT - 17782984.04 / RT)
      K(R97F) = 1.8000000000D+03 &
	   * exp(2 * LT - 10460000 / RT)
      K(R97B) = 7.8638867801D+00 &
	   * exp(2.39081 * LT - 39739288.47 / RT)
      K(R98F) = 2.2700000000D+02 &
	   * exp(2 * LT - 38493000 / RT)
      K(R98B) = 3.2832975054D+01 &
	   * exp(2.19425 * LT - 13001162.3 / RT)
      K(R99F) = 1.9200000000D+04 &
	   * exp(1.83 * LT - 920000 / RT)
      K(R99B) = 2.0159534780D-01 &
	   * exp(2.71516 * LT - 107010420.4 / RT)
      K(R100F) = 5.0000000000D+09
      K(R100B) = 5.4793052566D+10 &
	   * exp(0.0832929 * LT - 348997307.5 / RT)
      K(R101F) = 1.5100000000D+04 &
	   * exp(1.91 * LT - 15648000 / RT)
      K(R101B) = 3.7676751508D+00 &
	   * exp(2.3594 * LT + 23547106.15 / RT)
      K(R102F) = 4.2150000000D+10 * exp(-240998000 / RT)
      K(R102B) = 1.0901576399D+08 &
	   * exp(0.122726 * LT + 20841030.75 / RT)
      K(R103F) = 9.6400000000D+10
      K(R103B) = 1.0286605207D+11 &
	   * exp(0.176917 * LT - 286289003.1 / RT)
      K(R104F) = 1.2100000000D+07 * exp(2494000 / RT)
      K(R104B) = 3.0907166297D+07 &
	   * exp(0.266547 * LT - 97269358.64 / RT)
      K(R105F) = 3.9000000000D+08
      K(R105B) = 1.4149455616D+11 &
	   * exp(-0.113262 * LT - 294226181.4 / RT)
      K(R106F) = 9.6000000000D+08
      K(R106B) = 2.4080304685D+12 &
	   * exp(-0.307517 * LT - 319718019.1 / RT)
      K(R107F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4247000 / RT)
      K(R107B) = 6.1431582523D+11 &
	   * exp(-0.823557 * LT - 375185477.3 / RT)
      K(R108F) = 1.3370000000D+03 &
	   * exp(1.61 * LT + 1607000 / RT)
      K(R108B) = 8.6738893086D+03 &
	   * exp(1.42521 * LT - 56271988.32 / RT)
      K0TROE = 3.8000000000D+34 &
	   * exp(-7.27 * LT - 30208000 / RT)
      KINFTROE = 5.6000000000D+09 * exp(-10042000 / RT)
      FCTROE = 0.2493 * EXP( -TEMP / 98.5 ) &
	   + 0.7507 * EXP( -TEMP / 1302 ) &
	   + 1 * EXP( -4167 / TEMP )
      K(R109F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM17) )
      K0TROE = 8.3078895790D+36 &
	   * exp(-7.09812 * LT - 176992360.4 / RT)
      KINFTROE = 1.2243205695D+12 &
	   * exp(0.171883 * LT - 156826360.4 / RT)
      K(R109B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM17) )
      K(R110F) = 4.0800000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(R110B) = 1.0646326490D-02 &
	   * exp(3.31614 * LT - 198897303.8 / RT)
      K(R111F) = 4.8300000000D-07 &
	   * exp(4 * LT + 8368000 / RT)
      K(R111B) = 8.9104911373D-10 &
	   * exp(4.84011 * LT - 219504326.1 / RT)
      K0TROE = 3.2000000000D+21 &
	   * exp(-3.14 * LT - 5146000 / RT)
      KINFTROE = 2.5000000000D+13 * exp(-0.8 * LT)
      FCTROE = 0.32 * EXP( -TEMP / 78 ) &
	   + 0.68 * EXP( -TEMP / 1995 ) &
	   + 1 * EXP( -5590 / TEMP )
      K(R112F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM18) )
      K0TROE = 3.0957026570D+26 &
	   * exp(-3.30227 * LT - 469378295.6 / RT)
      KINFTROE = 2.4185177008D+18 &
	   * exp(-0.962268 * LT - 464232295.6 / RT)
      K(R112B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM18) )
      K(R113F) = 8.0000000000D+10
      K(R113B) = 7.5862942833D+12 &
	   * exp(-0.345989 * LT - 381866822.5 / RT)
      K(R114F) = 2.0000000000D+10
      K(R114B) = 2.2665974682D+15 &
	   * exp(-0.780313 * LT - 333139104.2 / RT)
      K(R115F) = 5.0000000000D+02 &
	   * exp(2 * LT - 30250000 / RT)
      K(R115B) = 2.0733714238D+05 &
	   * exp(1.48893 * LT - 61408932.1 / RT)
      K(R116F) = 1.3200000000D+10 * exp(-6276000 / RT)
      K(R116B) = 2.4890153347D+09 &
	   * exp(0.0889686 * LT - 317581691.3 / RT)
      K(R117F) = 2.0000000000D+10
      K(R117B) = 4.3480359040D+11 &
	   * exp(-0.018682 * LT - 485221897.6 / RT)
      K(R118F) = 3.2000000000D+10
      K(R118B) = 5.4402246782D+16 &
	   * exp(-1.08087 * LT - 559795321.3 / RT)
      K(R119F) = 9.0000000000D+09 * exp(-2510000 / RT)
      K(R119B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(R120F) = 3.0000000000D+10
      K(R120B) = 2.1872817488D+10 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R121F) = 9.0000000000D+09
      K(R121B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R122F) = 7.0000000000D+09
      K(R122B) = 5.1036574138D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R123F) = 1.5000000000D+10
      K(R123B) = 4.8515577882D+10 &
	   * exp(0.172986 * LT - 788045613 / RT)
      K(R124F) = 1.5000000000D+10
      K(R124B) = 1.0370851892D+12 &
	   * exp(-0.408275 * LT - 419169810.4 / RT)
      K(R125F) = 3.0000000000D+10
      K(R125B) = 2.4788436370D+15 &
	   * exp(-0.842599 * LT - 370442092.1 / RT)
      K(R126F) = 7.0000000000D+10
      K(R126B) = 2.1163554877D+13 &
	   * exp(-0.573354 * LT - 68461920.07 / RT)
      K(R127F) = 2.8000000000D+10
      K(R127B) = 7.7189773425D+05 &
	   * exp(0.259144 * LT - 284411118.3 / RT)
      K(R128F) = 1.2000000000D+10
      K(R128B) = 7.9258453041D+08 &
	   * exp(0.51432 * LT - 780192786.3 / RT)
      K(R129F) = 1.4000000000D+10
      K(R129B) = 2.0690264836D+08 &
	   * exp(0.436598 * LT - 262972337 / RT)
      K(R130F) = 7.0000000000D+15 * exp(-341724000 / RT)
      K(R130B) = 3.7753643217D+04 &
	   * exp(1.51966 * LT + 20713522.45 / RT)
      K(R131F) = 1.6990000000D+42 &
	   * exp(-7.954 * LT - 384119000 / RT)
      K(R131B) = 2.1946663920D+31 &
	   * exp(-6.47478 * LT - 33101635.18 / RT)
      K(R132F) = 6.7100000000D+03 &
	   * exp(2 * LT + 2635000 / RT)
      K(R132B) = 2.5724684790D+02 &
	   * exp(2.10394 * LT - 86734824.77 / RT)
      K(R133F) = 2.9700000000D+04 &
	   * exp(2 * LT - 16877000 / RT)
      K(R133B) = 1.1087242276D+02 &
	   * exp(2.19757 * LT - 43538520.34 / RT)
      K(R134F) = 2.6800000000D-02 &
	   * exp(3.778 * LT - 40297000 / RT)
      K(R134B) = 3.4015900172D-02 &
	   * exp(3.68539 * LT - 74895698.6 / RT)
      K(R135F) = 1.8550000000D-06 &
	   * exp(5.29 * LT + 456000 / RT)
      K(R135B) = 4.0616566920D-09 &
	   * exp(5.45253 * LT - 20439430.15 / RT)
      K(R136F) = 2.0000000000D+10 * exp(-69036000 / RT)
      K(R136B) = 6.8709845478D+10 &
	   * exp(-0.553412 * LT - 29363177.65 / RT)
      K(R137F) = 4.1000000000D+10 * exp(-187903000 / RT)
      K(R137B) = 9.3054713257D+08 &
	   * exp(-0.16414 * LT + 13845494.45 / RT)
      K(R138F) = 1.2000000000D+13 * exp(-107738000 / RT)
      K(R138B) = 1.2688738503D+02 &
	   * exp(1.9314 * LT - 76352192.59 / RT)
      K(R139F) = 2.4100000000D+10
      K(R139B) = 1.9727820584D+10 &
	   * exp(0.452187 * LT - 319631557.4 / RT)
      K(R140F) = 5.4900000000D+00 &
	   * exp(2.8 * LT - 24527000 / RT)
      K(R140B) = 2.0980138029D+00 &
	   * exp(3.07179 * LT - 52359288.14 / RT)
      K(R141F) = 9.0000000000D+09
      K(R141B) = 4.9424652594D+11 &
	   * exp(-0.286081 * LT - 124405355.8 / RT)
      K(R142F) = 1.7450000000D+16 &
	   * exp(-0.66 * LT - 49036000 / RT)
      K(R142B) = 4.2975913579D+12 &
	   * exp(-0.727194 * LT - 40420853 / RT)
      K(R143F) = 1.0000000000D+10 * exp(-207945000 / RT)
      K(R143B) = 2.6565778440D+07 &
	   * exp(-1.12449 * LT - 3300998.939 / RT)
      K(R144F) = 2.3400000000D+04 &
	   * exp(1.61 * LT + 146000 / RT)
      K(R144B) = 1.0500535004D+02 &
	   * exp(0.753594 * LT - 86328318.16 / RT)
      K(R145F) = 1.2200000000D+09 * exp(-71128000 / RT)
      K(R145B) = 4.9058790251D+08 &
	   * exp(-1.51376 * LT - 28559671.04 / RT)
      K(R146F) = 2.3500000000D+02 &
	   * exp(2.5 * LT - 9330000 / RT)
      K(R146B) = 6.0227548826D-02 &
	   * exp(1.70218 * LT - 27329923.54 / RT)
      K(R147F) = 4.5500000000D+03 &
	   * exp(2 * LT - 20920000 / RT)
      K(R147B) = 1.9881379390D+00 &
	   * exp(1.23722 * LT - 44686013.73 / RT)
      K(R148F) = 7.5500000000D-04 &
	   * exp(3.46 * LT - 22933000 / RT)
      K(R148B) = 1.1216623790D-04 &
	   * exp(2.40704 * LT - 54636191.99 / RT)
      K(R149F) = 7.4510000000D+12 &
	   * exp(-1.76 * LT - 71756000 / RT)
      K(R149B) = 6.5430526784D+01 &
	   * exp(1.83191 * LT + 1119894.578 / RT)
      K(R150F) = 1.5140000000D+12 &
	   * exp(-1.78 * LT - 57823000 / RT)
      K(R150B) = 8.3113859279D+02 &
	   * exp(1.45172 * LT - 136729766.3 / RT)
      K(R151F) = 2.0000000000D+09
      K(R151B) = 7.7901495616D+19 &
	   * exp(-1.66077 * LT - 152376875.2 / RT)
      K(R152F) = 1.5970000000D+20 * exp(-4.5 * LT)
      K(R152B) = 2.3471800921D+09 &
	   * exp(-2.52515 * LT + 3362483.601 / RT)
      K(R153F) = 6.8440000000D+19 * exp(-4.5 * LT)
      K(R153B) = 1.2073443805D+09 &
	   * exp(-2.46046 * LT - 422480241.3 / RT)
      K(R154F) = 9.7220000000D+15 &
	   * exp(-1.1 * LT - 86358000 / RT)
      K(R154B) = 3.2115786960D+04 &
	   * exp(0.960117 * LT - 38336724.59 / RT)
      K(R155F) = 5.0000000000D+07 * exp(-2092000 / RT)
      K(R155B) = 1.7465831985D+07 &
	   * exp(-0.0801019 * LT - 198140201.7 / RT)
      K(R156F) = 6.0000000000D+10 * exp(-89956000 / RT)
      K(R156B) = 2.7310450662D+12 &
	   * exp(-0.87687 * LT - 46223662.18 / RT)
      K(R157F) = 1.5000000000D+13 * exp(-85772000 / RT)
      K(R157B) = 2.8515336158D-11 &
	   * exp(4.59973 * LT - 171394605.9 / RT)
      K(R158F) = 7.0000000000D+08
      K(R158B) = 3.0813052680D+19 &
	   * exp(-1.67723 * LT - 152474430.4 / RT)
      K(R159F) = 4.0000000000D+10 * exp(-77404000 / RT)
      K(R159B) = 3.4129896457D+00 &
	   * exp(1.45568 * LT - 248443163.9 / RT)
      K(R160F) = 3.0000000000D+16 * exp(-167360000 / RT)
      K(R160B) = 7.8411183626D+04 &
	   * exp(2.03667 * LT + 19479411.25 / RT)
      K(R161F) = 1.0000000000D+11 * exp(-58576000 / RT)
      K(R161B) = 1.0908995697D+09 &
	   * exp(0.529885 * LT - 85941547.44 / RT)
      K(R162F) = 2.1770000000D+16 &
	   * exp(-2.69 * LT - 71965000 / RT)
      K(R162B) = 1.1879914769D+05 &
	   * exp(-0.563231 * LT - 13891644.72 / RT)
      K(R163F) = 5.3110000000D+15 &
	   * exp(-2.61 * LT - 87069000 / RT)
      K(R163B) = 1.0900751837D+06 &
	   * exp(-0.852332 * LT - 148319468.8 / RT)
      K(R164F) = 1.0000000000D+14 * exp(-62342000 / RT)
      K(R164B) = 1.4245687646D+11 &
	   * exp(-0.102835 * LT - 41360025.96 / RT)
      K(R165F) = 4.5000000000D+12 * exp(-1.11 * LT)
      K(R165B) = 4.1462884067D+21 &
	   * exp(-2.3397 * LT - 106374743.1 / RT)
      K(R166F) = 2.3000000000D+10 * exp(-209200000 / RT)
      K(R166B) = 5.1478464291D-02 &
	   * exp(1.94073 * LT - 176811783 / RT)
      K(R167F) = 1.5000000000D+13 * exp(-238488000 / RT)
      K(R167B) = 1.8277568334D+07 &
	   * exp(0.755153 * LT - 250861233.7 / RT)
      K(R168F) = 4.5930000000D+18 &
	   * exp(-0.46 * LT - 453127000 / RT)
      K(R168B) = 2.1397693475D+07 &
	   * exp(0.993088 * LT + 10845324.08 / RT)
      K(R169F) = 2.6200000000D+03 &
	   * exp(2.06 * LT - 3833000 / RT)
      K(R169B) = 3.2786035279D-02 &
	   * exp(2.72153 * LT - 78914538.1 / RT)
      K(R170F) = 1.8500000000D+04 &
	   * exp(1.51 * LT + 4025000 / RT)
      K(R170B) = 4.1406590841D-08 &
	   * exp(3.45073 * LT + 36413217.05 / RT)
      K(R171F) = 4.2400000000D+03 &
	   * exp(2.1 * LT - 20368000 / RT)
      K(R171B) = 5.1664593155D-03 &
	   * exp(2.85515 * LT - 32741233.67 / RT)
      K(R172F) = 6.0300000000D+10 &
	   * exp(-0.35 * LT - 12502000 / RT)
      K(R172B) = 1.3141792487D-02 &
	   * exp(1.68435 * LT + 82594521.48 / RT)
      K(R173F) = 3.9000000000D-10 &
	   * exp(5.8 * LT - 9205000 / RT)
      K(R173B) = 2.8898944836D-20 &
	   * exp(7.54417 * LT + 77954343.22 / RT)
      K(R174F) = 1.0000000000D+09 * exp(-49873000 / RT)
      K(R174B) = 2.0056686253D-01 &
	   * exp(1.28337 * LT + 111557864.2 / RT)
      K(R175F) = 1.7700000000D+15 &
	   * exp(-1.9 * LT - 12447000 / RT)
      K(R175B) = 2.2625738811D+02 &
	   * exp(0.0993154 * LT + 88415611.66 / RT)

      W(R1F) = K(R1F) * C(SH) * C(SO2)
      W(R1B) = K(R1B) * C(SOH) * C(SO)
      W(R2F) = K(R2F) * C(SO) * C(SH2)
      W(R2B) = K(R2B) * C(SOH) * C(SH)
      W(R3F) = K(R3F) * C(SH2) * C(SOH)
      W(R3B) = K(R3B) * C(SH) * C(SH2O)
      W(R4F) = K(R4F) * C(SO) * C(SH2O)
      W(R4B) = K(R4B) * C(SOH) * C(SOH)
      W(R5F) = K(R5F) * C(SH2) * M(MM1)
      W(R5B) = K(R5B) * C(SH) * C(SH) * M(MM1)
      W(R6F) = K(R6F) * C(SO) * C(SO) * M(MM2)
      W(R6B) = K(R6B) * C(SO2) * M(MM2)
      W(R7F) = K(R7F) * C(SO) * C(SH) * M(MM3)
      W(R7B) = K(R7B) * C(SOH) * M(MM3)
      W(R8F) = K(R8F) * C(SH) * C(SOH) * M(MM4)
      W(R8B) = K(R8B) * C(SH2O) * M(MM4)
      W(R9F) = K(R9F) * C(SH) * C(SO2)
      W(R9B) = K(R9B) * C(SHO2)
      W(R10F) = K(R10F) * C(SHO2) * C(SH)
      W(R10B) = K(R10B) * C(SO2) * C(SH2)
      W(R11F) = K(R11F) * C(SHO2) * C(SH)
      W(R11B) = K(R11B) * C(SOH) * C(SOH)
      W(R12F) = K(R12F) * C(SHO2) * C(SO)
      W(R12B) = K(R12B) * C(SOH) * C(SO2)
      W(R13F) = K(R13F) * C(SHO2) * C(SOH)
      W(R13B) = K(R13B) * C(SO2) * C(SH2O)
      W(R14F) = K(R14F) * C(SHO2) * C(SHO2)
      W(R14B) = K(R14B) * C(SO2) * C(SH2O2)
      W(R15F) = K(R15F) * C(SHO2) * C(SHO2)
      W(R15B) = K(R15B) * C(SO2) * C(SH2O2)
      W(R16F) = K(R16F) * C(SH2O2)
      W(R16B) = K(R16B) * C(SOH) * C(SOH)
      W(R17F) = K(R17F) * C(SH2O2) * C(SH)
      W(R17B) = K(R17B) * C(SOH) * C(SH2O)
      W(R18F) = K(R18F) * C(SH2O2) * C(SH)
      W(R18B) = K(R18B) * C(SH2) * C(SHO2)
      W(R19F) = K(R19F) * C(SH2O2) * C(SO)
      W(R19B) = K(R19B) * C(SHO2) * C(SOH)
      W(R20F) = K(R20F) * C(SH2O2) * C(SOH)
      W(R20B) = K(R20B) * C(SH2O) * C(SHO2)
      W(R21F) = K(R21F) * C(SH2O2) * C(SOH)
      W(R21B) = K(R21B) * C(SH2O) * C(SHO2)
      W(R22F) = K(R22F) * C(SCO) * C(SO)
      W(R22B) = K(R22B) * C(SCO2)
      W(R23F) = K(R23F) * C(SCO) * C(SO2)
      W(R23B) = K(R23B) * C(SO) * C(SCO2)
      W(R24F) = K(R24F) * C(SCO) * C(SHO2)
      W(R24B) = K(R24B) * C(SOH) * C(SCO2)
      W(R25F) = K(R25F) * C(SCO) * C(SOH)
      W(R25B) = K(R25B) * C(SH) * C(SCO2)
      W(R26F) = K(R26F) * C(SHCO) * M(MM8)
      W(R26B) = K(R26B) * C(SCO) * C(SH) * M(MM8)
      W(R27F) = K(R27F) * C(SHCO) * C(SO2)
      W(R27B) = K(R27B) * C(SHO2) * C(SCO)
      W(R28F) = K(R28F) * C(SHCO) * C(SH)
      W(R28B) = K(R28B) * C(SH2) * C(SCO)
      W(R29F) = K(R29F) * C(SHCO) * C(SO)
      W(R29B) = K(R29B) * C(SOH) * C(SCO)
      W(R30F) = K(R30F) * C(SHCO) * C(SOH)
      W(R30B) = K(R30B) * C(SH2O) * C(SCO)
      W(R31F) = K(R31F) * C(SHCO) * C(SO)
      W(R31B) = K(R31B) * C(SH) * C(SCO2)
      W(R32F) = K(R32F) * C(SHCO) * C(SHO2)
      W(R32B) = K(R32B) * C(SH) * C(SOH) * C(SCO2)
      W(R33F) = K(R33F) * C(SHCO) * C(SHCO)
      W(R33B) = K(R33B) * C(SCO) * C(SCO) * C(SH2)
      W(R34F) = K(R34F) * C(SHCO) * C(SCH3)
      W(R34B) = K(R34B) * C(SCH4) * C(SCO)
      W(R35F) = K(R35F) * C(SHCO) * C(SHCO)
      W(R35B) = K(R35B) * C(SCO) * C(SCH2O)
      W(R36F) = K(R36F) * C(SCH2O) * M(MM9)
      W(R36B) = K(R36B) * C(SH) * C(SHCO) * M(MM9)
      W(R37F) = K(R37F) * C(SCH2O) * M(MM10)
      W(R37B) = K(R37B) * C(SH2) * C(SCO) * M(MM10)
      W(R38F) = K(R38F) * C(SCH2O) * C(SH)
      W(R38B) = K(R38B) * C(SH2) * C(SHCO)
      W(R39F) = K(R39F) * C(SCH2O) * C(SO)
      W(R39B) = K(R39B) * C(SOH) * C(SHCO)
      W(R40F) = K(R40F) * C(SCH2O) * C(SOH)
      W(R40B) = K(R40B) * C(SH2O) * C(SHCO)
      W(R41F) = K(R41F) * C(SCH2O) * C(SO2)
      W(R41B) = K(R41B) * C(SHO2) * C(SHCO)
      W(R42F) = K(R42F) * C(SCH2O) * C(SHO2)
      W(R42B) = K(R42B) * C(SH2O2) * C(SHCO)
      W(R43F) = K(R43F) * C(SCH2O) * C(SCH3)
      W(R43B) = K(R43B) * C(SCH4) * C(SHCO)
      W(R44F) = K(R44F) * C(SCH3) * C(SO)
      W(R44B) = K(R44B) * C(SH) * C(SCH2O)
      W(R45F) = K(R45F) * C(SCH3) * C(SO2)
      W(R45B) = K(R45B) * C(SO) * C(SCH3O)
      W(R46F) = K(R46F) * C(SCH3) * C(SO2)
      W(R46B) = K(R46B) * C(SOH) * C(SCH2O)
      W(R47F) = K(R47F) * C(SCH3) * C(SHO2)
      W(R47B) = K(R47B) * C(SOH) * C(SCH3O)
      W(R48F) = K(R48F) * C(SCH3) * C(SCH3)
      W(R48B) = K(R48B) * C(SC2H6)
      W(R49F) = K(R49F) * C(SCH3) * C(SH)
      W(R49B) = K(R49B) * C(SCH4)
      W(R50F) = K(R50F) * C(SCH4) * C(SH)
      W(R50B) = K(R50B) * C(SH2) * C(SCH3)
      W(R51F) = K(R51F) * C(SCH4) * C(SO)
      W(R51B) = K(R51B) * C(SOH) * C(SCH3)
      W(R52F) = K(R52F) * C(SCH4) * C(SOH)
      W(R52B) = K(R52B) * C(SH2O) * C(SCH3)
      W(R53F) = K(R53F) * C(SCH3) * C(SHO2)
      W(R53B) = K(R53B) * C(SO2) * C(SCH4)
      W(R54F) = K(R54F) * C(SCH4) * C(SHO2)
      W(R54B) = K(R54B) * C(SH2O2) * C(SCH3)
      W(R55F) = K(R55F) * C(SCH2OH) * M(MM0)
      W(R55B) = K(R55B) * C(SH) * C(SCH2O) * M(MM0)
      W(R56F) = K(R56F) * C(SCH2OH) * C(SH)
      W(R56B) = K(R56B) * C(SH2) * C(SCH2O)
      W(R57F) = K(R57F) * C(SCH2OH) * C(SH)
      W(R57B) = K(R57B) * C(SOH) * C(SCH3)
      W(R58F) = K(R58F) * C(SCH2OH) * C(SO)
      W(R58B) = K(R58B) * C(SOH) * C(SCH2O)
      W(R59F) = K(R59F) * C(SCH2OH) * C(SOH)
      W(R59B) = K(R59B) * C(SH2O) * C(SCH2O)
      W(R60F) = K(R60F) * C(SCH2OH) * C(SO2)
      W(R60B) = K(R60B) * C(SHO2) * C(SCH2O)
      W(R61F) = K(R61F) * C(SCH2OH) * C(SO2)
      W(R61B) = K(R61B) * C(SHO2) * C(SCH2O)
      W(R62F) = K(R62F) * C(SCH2OH) * C(SHO2)
      W(R62B) = K(R62B) * C(SH2O2) * C(SCH2O)
      W(R63F) = K(R63F) * C(SCH2OH) * C(SHCO)
      W(R63B) = K(R63B) * C(SCH2O) * C(SCH2O)
      W(R64F) = K(R64F) * C(SCH3O) * M(MM0)
      W(R64B) = K(R64B) * C(SH) * C(SCH2O) * M(MM0)
      W(R65F) = K(R65F) * C(SCH3O) * C(SH)
      W(R65B) = K(R65B) * C(SOH) * C(SCH3)
      W(R66F) = K(R66F) * C(SCH3O) * C(SO)
      W(R66B) = K(R66B) * C(SOH) * C(SCH2O)
      W(R67F) = K(R67F) * C(SCH3O) * C(SOH)
      W(R67B) = K(R67B) * C(SH2O) * C(SCH2O)
      W(R68F) = K(R68F) * C(SCH3O) * C(SO2)
      W(R68B) = K(R68B) * C(SHO2) * C(SCH2O)
      W(R69F) = K(R69F) * C(SCH3O) * C(SO2)
      W(R69B) = K(R69B) * C(SHO2) * C(SCH2O)
      W(R70F) = K(R70F) * C(SCH3O) * C(SHO2)
      W(R70B) = K(R70B) * C(SH2O2) * C(SCH2O)
      W(R71F) = K(R71F) * C(SCH3O) * C(SCO)
      W(R71B) = K(R71B) * C(SCO2) * C(SCH3)
      W(R72F) = K(R72F) * C(SCH3) * C(SCH3)
      W(R72B) = K(R72B) * C(SC2H5) * C(SH)
      W(R73F) = K(R73F) * C(SCH4) * C(SCH2)
      W(R73B) = K(R73B) * C(SCH3) * C(SCH3)
      W(R74F) = K(R74F) * C(SCH4) * C(SCH2X)
      W(R74B) = K(R74B) * C(SCH3) * C(SCH3)
      W(R75F) = K(R75F) * C(SCH3) * C(SOH)
      W(R75B) = K(R75B) * C(SH2O) * C(SCH2)
      W(R76F) = K(R76F) * C(SCH3) * C(SOH)
      W(R76B) = K(R76B) * C(SH2O) * C(SCH2X)
      W(R77F) = K(R77F) * C(SCH3) * C(SCH2)
      W(R77B) = K(R77B) * C(SH) * C(SC2H4)
      W(R78F) = K(R78F) * C(SCH3) * C(SCH2X)
      W(R78B) = K(R78B) * C(SH) * C(SC2H4)
      W(R79F) = K(R79F) * C(SCH3O) * C(SH)
      W(R79B) = K(R79B) * C(SH2O) * C(SCH2X)
      W(R80F) = K(R80F) * C(SC2H6) * C(SH)
      W(R80B) = K(R80B) * C(SH2) * C(SC2H5)
      W(R81F) = K(R81F) * C(SC2H6) * C(SO)
      W(R81B) = K(R81B) * C(SOH) * C(SC2H5)
      W(R82F) = K(R82F) * C(SC2H6) * C(SOH)
      W(R82B) = K(R82B) * C(SH2O) * C(SC2H5)
      W(R83F) = K(R83F) * C(SC2H6) * C(SO2)
      W(R83B) = K(R83B) * C(SHO2) * C(SC2H5)
      W(R84F) = K(R84F) * C(SC2H6) * C(SHO2)
      W(R84B) = K(R84B) * C(SH2O2) * C(SC2H5)
      W(R85F) = K(R85F) * C(SC2H6) * C(SCH3)
      W(R85B) = K(R85B) * C(SCH4) * C(SC2H5)
      W(R86F) = K(R86F) * C(SC2H5) * C(SH)
      W(R86B) = K(R86B) * C(SC2H6)
      W(R87F) = K(R87F) * C(SC2H5) * C(SH)
      W(R87B) = K(R87B) * C(SH2) * C(SC2H4)
      W(R88F) = K(R88F) * C(SC2H5) * C(SO)
      W(R88B) = K(R88B) * C(SCH2O) * C(SCH3)
      W(R89F) = K(R89F) * C(SC2H5) * C(SO2)
      W(R89B) = K(R89B) * C(SHO2) * C(SC2H4)
      W(R90F) = K(R90F) * C(SC2H5) * C(SC2H5)
      W(R90B) = K(R90B) * C(SC2H6) * C(SC2H4)
      W(R91F) = K(R91F) * C(SC2H5) * C(SHCO)
      W(R91B) = K(R91B) * C(SCO) * C(SC2H6)
      W(R92F) = K(R92F) * C(SC2H5) * C(SO)
      W(R92B) = K(R92B) * C(SH) * C(SCH3HCO)
      W(R93F) = K(R93F) * C(SC2H4)
      W(R93B) = K(R93B) * C(SC2H2) * C(SH2)
      W(R94F) = K(R94F) * C(SC2H4) * C(SH)
      W(R94B) = K(R94B) * C(SC2H5)
      W(R95F) = K(R95F) * C(SC2H3) * C(SH)
      W(R95B) = K(R95B) * C(SC2H4)
      W(R96F) = K(R96F) * C(SC2H4) * C(SH)
      W(R96B) = K(R96B) * C(SH2) * C(SC2H3)
      W(R97F) = K(R97F) * C(SC2H4) * C(SOH)
      W(R97B) = K(R97B) * C(SH2O) * C(SC2H3)
      W(R98F) = K(R98F) * C(SC2H4) * C(SCH3)
      W(R98B) = K(R98B) * C(SCH4) * C(SC2H3)
      W(R99F) = K(R99F) * C(SC2H4) * C(SO)
      W(R99B) = K(R99B) * C(SHCO) * C(SCH3)
      W(R100F) = K(R100F) * C(SC2H3) * C(SOH)
      W(R100B) = K(R100B) * C(SH2O) * C(SC2H2)
      W(R101F) = K(R101F) * C(SC2H4) * C(SO)
      W(R101B) = K(R101B) * C(SC2H3) * C(SOH)
      W(R102F) = K(R102F) * C(SC2H4) * C(SO2)
      W(R102B) = K(R102B) * C(SHO2) * C(SC2H3)
      W(R103F) = K(R103F) * C(SC2H3) * C(SH)
      W(R103B) = K(R103B) * C(SH2) * C(SC2H2)
      W(R104F) = K(R104F) * C(SC2H3) * C(SH2O2)
      W(R104B) = K(R104B) * C(SHO2) * C(SC2H4)
      W(R105F) = K(R105F) * C(SC2H3) * C(SCH3)
      W(R105B) = K(R105B) * C(SCH4) * C(SC2H2)
      W(R106F) = K(R106F) * C(SC2H3) * C(SC2H3)
      W(R106B) = K(R106B) * C(SC2H2) * C(SC2H4)
      W(R107F) = K(R107F) * C(SC2H3) * C(SO2)
      W(R107B) = K(R107B) * C(SCH2O) * C(SHCO)
      W(R108F) = K(R108F) * C(SC2H3) * C(SO2)
      W(R108B) = K(R108B) * C(SC2H2) * C(SHO2)
      W(R109F) = K(R109F) * C(SC2H2) * C(SH)
      W(R109B) = K(R109B) * C(SC2H3)
      W(R110F) = K(R110F) * C(SC2H2) * C(SO)
      W(R110B) = K(R110B) * C(SCO) * C(SCH2)
      W(R111F) = K(R111F) * C(SC2H2) * C(SOH)
      W(R111B) = K(R111B) * C(SCO) * C(SCH3)
      W(R112F) = K(R112F) * C(SCH2) * C(SH)
      W(R112B) = K(R112B) * C(SCH3)
      W(R113F) = K(R113F) * C(SCH2) * C(SO)
      W(R113B) = K(R113B) * C(SH) * C(SHCO)
      W(R114F) = K(R114F) * C(SCH2) * C(SOH)
      W(R114B) = K(R114B) * C(SH) * C(SCH2O)
      W(R115F) = K(R115F) * C(SCH2) * C(SH2)
      W(R115B) = K(R115B) * C(SCH3) * C(SH)
      W(R116F) = K(R116F) * C(SCH2) * C(SO2)
      W(R116B) = K(R116B) * C(SOH) * C(SHCO)
      W(R117F) = K(R117F) * C(SCH2) * C(SHO2)
      W(R117B) = K(R117B) * C(SOH) * C(SCH2O)
      W(R118F) = K(R118F) * C(SCH2) * C(SCH2)
      W(R118B) = K(R118B) * C(SH2) * C(SC2H2)
      W(R119F) = K(R119F) * C(SCH2X) * M(MM19)
      W(R119B) = K(R119B) * C(SCH2) * M(MM19)
      W(R120F) = K(R120F) * C(SCH2X) * C(SH2O)
      W(R120B) = K(R120B) * C(SH2O) * C(SCH2)
      W(R121F) = K(R121F) * C(SCH2X) * C(SCO)
      W(R121B) = K(R121B) * C(SCO) * C(SCH2)
      W(R122F) = K(R122F) * C(SCH2X) * C(SCO2)
      W(R122B) = K(R122B) * C(SCO2) * C(SCH2)
      W(R123F) = K(R123F) * C(SCH2X) * C(SO)
      W(R123B) = K(R123B) * C(SH2) * C(SCO)
      W(R124F) = K(R124F) * C(SCH2X) * C(SO)
      W(R124B) = K(R124B) * C(SH) * C(SHCO)
      W(R125F) = K(R125F) * C(SCH2X) * C(SOH)
      W(R125B) = K(R125B) * C(SH) * C(SCH2O)
      W(R126F) = K(R126F) * C(SCH2X) * C(SH2)
      W(R126B) = K(R126B) * C(SH) * C(SCH3)
      W(R127F) = K(R127F) * C(SCH2X) * C(SO2)
      W(R127B) = K(R127B) * C(SCO) * C(SOH) * C(SH)
      W(R128F) = K(R128F) * C(SCH2X) * C(SO2)
      W(R128B) = K(R128B) * C(SH2O) * C(SCO)
      W(R129F) = K(R129F) * C(SCH2X) * C(SCO2)
      W(R129B) = K(R129B) * C(SCO) * C(SCH2O)
      W(R130F) = K(R130F) * C(SCH3HCO)
      W(R130B) = K(R130B) * C(SHCO) * C(SCH3)
      W(R131F) = K(R131F) * C(SCH3OCH3)
      W(R131B) = K(R131B) * C(SCH3O) * C(SCH3)
      W(R132F) = K(R132F) * C(SCH3OCH3) * C(SOH)
      W(R132B) = K(R132B) * C(SH2O) * C(SCH3OCH2)
      W(R133F) = K(R133F) * C(SCH3OCH3) * C(SH)
      W(R133B) = K(R133B) * C(SH2) * C(SCH3OCH2)
      W(R134F) = K(R134F) * C(SCH3OCH3) * C(SCH3)
      W(R134B) = K(R134B) * C(SCH4) * C(SCH3OCH2)
      W(R135F) = K(R135F) * C(SCH3OCH3) * C(SO)
      W(R135B) = K(R135B) * C(SOH) * C(SCH3OCH2)
      W(R136F) = K(R136F) * C(SCH3OCH3) * C(SHO2)
      W(R136B) = K(R136B) * C(SH2O2) * C(SCH3OCH2)
      W(R137F) = K(R137F) * C(SCH3OCH3) * C(SO2)
      W(R137B) = K(R137B) * C(SHO2) * C(SCH3OCH2)
      W(R138F) = K(R138F) * C(SCH3OCH2)
      W(R138B) = K(R138B) * C(SCH3) * C(SCH2O)
      W(R139F) = K(R139F) * C(SCH3OCH2) * C(SCH3O)
      W(R139B) = K(R139B) * C(SCH2O) * C(SCH3OCH3)
      W(R140F) = K(R140F) * C(SCH3OCH2) * C(SCH2O)
      W(R140B) = K(R140B) * C(SHCO) * C(SCH3OCH3)
      W(R141F) = K(R141F) * C(SCH3OCH2) * C(SHO2)
      W(R141B) = K(R141B) * C(SOH) * C(SCH3OCH2O)
      W(R142F) = K(R142F) * C(SCH3OCH2O)
      W(R142B) = K(R142B) * C(SH) * C(SCH3OCHO)
      W(R143F) = K(R143F) * C(SCH3OCHO) * C(SO2)
      W(R143B) = K(R143B) * C(SHO2) * C(SCH3OCO)
      W(R144F) = K(R144F) * C(SCH3OCHO) * C(SOH)
      W(R144B) = K(R144B) * C(SH2O) * C(SCH3OCO)
      W(R145F) = K(R145F) * C(SCH3OCHO) * C(SHO2)
      W(R145B) = K(R145B) * C(SH2O2) * C(SCH3OCO)
      W(R146F) = K(R146F) * C(SCH3OCHO) * C(SO)
      W(R146B) = K(R146B) * C(SOH) * C(SCH3OCO)
      W(R147F) = K(R147F) * C(SCH3OCHO) * C(SH)
      W(R147B) = K(R147B) * C(SH2) * C(SCH3OCO)
      W(R148F) = K(R148F) * C(SCH3OCHO) * C(SCH3)
      W(R148B) = K(R148B) * C(SCH4) * C(SCH3OCO)
      W(R149F) = K(R149F) * C(SCH3OCO)
      W(R149B) = K(R149B) * C(SCO) * C(SCH3O)
      W(R150F) = K(R150F) * C(SCH3OCO)
      W(R150B) = K(R150B) * C(SCO2) * C(SCH3)
      W(R151F) = K(R151F) * C(SCH3OCH2) * C(SO2)
      W(R151B) = K(R151B) * C(SRO2)
      W(R152F) = K(R152F) * C(SRO2) * C(SRO2)
      W(R152B) = K(R152B) * C(SCH3OCH2O) * C(SCH3OCH2O) * C(SO2)
      W(R153F) = K(R153F) * C(SRO2) * C(SRO2)
      W(R153B) = K(R153B) * C(SROH) * C(SCH3OCHO) * C(SO2)
      W(R154F) = K(R154F) * C(SCH3OCH2O)
      W(R154B) = K(R154B) * C(SCH2O) * C(SCH3O)
      W(R155F) = K(R155F) * C(SCH3OCH2O) * C(SO2)
      W(R155B) = K(R155B) * C(SHO2) * C(SCH3OCHO)
      W(R156F) = K(R156F) * C(SRO2)
      W(R156B) = K(R156B) * C(SQOOH)
      W(R157F) = K(R157F) * C(SQOOH)
      W(R157B) = K(R157B) * C(SCH2O) * C(SCH2O) * C(SOH)
      W(R158F) = K(R158F) * C(SQOOH) * C(SO2)
      W(R158B) = K(R158B) * C(SO2QOOH)
      W(R159F) = K(R159F) * C(SO2QOOH)
      W(R159B) = K(R159B) * C(SOH) * C(SHO2QHO)
      W(R160F) = K(R160F) * C(SHO2QHO)
      W(R160B) = K(R160B) * C(SOH) * C(SOCH2OCHO)
      W(R161F) = K(R161F) * C(SOCH2OCHO)
      W(R161B) = K(R161B) * C(SHOCH2OCO)
      W(R162F) = K(R162F) * C(SHOCH2OCO)
      W(R162B) = K(R162B) * C(SCO) * C(SHOCH2O)
      W(R163F) = K(R163F) * C(SHOCH2OCO)
      W(R163B) = K(R163B) * C(SCO2) * C(SCH2OH)
      W(R164F) = K(R164F) * C(SHOCH2O)
      W(R164B) = K(R164B) * C(SH) * C(SHCOOH)
      W(R165F) = K(R165F) * C(SCH2O) * C(SOH)
      W(R165B) = K(R165B) * C(SHOCH2O)
      W(R166F) = K(R166F) * C(SHCOOH) * M(MM0)
      W(R166B) = K(R166B) * C(SH2O) * C(SCO) * M(MM0)
      W(R167F) = K(R167F) * C(SHCOOH) * M(MM0)
      W(R167B) = K(R167B) * C(SH2) * C(SCO2) * M(MM0)
      W(R168F) = K(R168F) * C(SHCOOH)
      W(R168B) = K(R168B) * C(SOH) * C(SHCO)
      W(R169F) = K(R169F) * C(SHCOOH) * C(SOH)
      W(R169B) = K(R169B) * C(SH) * C(SCO2) * C(SH2O)
      W(R170F) = K(R170F) * C(SHCOOH) * C(SOH)
      W(R170B) = K(R170B) * C(SOH) * C(SCO) * C(SH2O)
      W(R171F) = K(R171F) * C(SHCOOH) * C(SH)
      W(R171B) = K(R171B) * C(SH) * C(SCO2) * C(SH2)
      W(R172F) = K(R172F) * C(SHCOOH) * C(SH)
      W(R172B) = K(R172B) * C(SOH) * C(SCO) * C(SH2)
      W(R173F) = K(R173F) * C(SHCOOH) * C(SCH3)
      W(R173B) = K(R173B) * C(SOH) * C(SCO) * C(SCH4)
      W(R174F) = K(R174F) * C(SHCOOH) * C(SHO2)
      W(R174B) = K(R174B) * C(SOH) * C(SCO) * C(SH2O2)
      W(R175F) = K(R175F) * C(SHCOOH) * C(SO)
      W(R175B) = K(R175B) * C(SOH) * C(SOH) * C(SCO)


      CDOT(SN2) = 0.0_DP
      CDOT(SH) = - W(R1F) + W(R1B) + W(R2F) & 
	    - W(R2B) + W(R3F) - W(R3B) & 
	    + 2 * W(R5F) - 2 * W(R5B) - W(R7F) & 
	    + W(R7B) - W(R8F) + W(R8B) & 
	    - W(R9F) + W(R9B) - W(R10F) & 
	    + W(R10B) - W(R11F) + W(R11B) & 
	    - W(R17F) + W(R17B) - W(R18F) & 
	    + W(R18B) + W(R25F) - W(R25B) & 
	    + W(R26F) - W(R26B) - W(R28F) & 
	    + W(R28B) + W(R31F) - W(R31B) & 
	    + W(R32F) - W(R32B) + W(R36F) & 
	    - W(R36B) - W(R38F) + W(R38B) & 
	    + W(R44F) - W(R44B) - W(R49F) & 
	    + W(R49B) - W(R50F) + W(R50B) & 
	    + W(R55F) - W(R55B) - W(R56F)
      CDOT(SH) = CDOT(SH) + W(R56B) - W(R57F) + W(R57B) & 
	    + W(R64F) - W(R64B) - W(R65F) & 
	    + W(R65B) + W(R72F) - W(R72B) & 
	    + W(R77F) - W(R77B) + W(R78F) & 
	    - W(R78B) - W(R79F) + W(R79B) & 
	    - W(R80F) + W(R80B) - W(R86F) & 
	    + W(R86B) - W(R87F) + W(R87B) & 
	    + W(R92F) - W(R92B) - W(R94F) & 
	    + W(R94B) - W(R95F) + W(R95B) & 
	    - W(R96F) + W(R96B) - W(R103F) & 
	    + W(R103B) - W(R109F) + W(R109B) & 
	    - W(R112F) + W(R112B) + W(R113F) & 
	    - W(R113B) + W(R114F) - W(R114B) & 
	    + W(R115F) - W(R115B) + W(R124F) & 
	    - W(R124B) + W(R125F) - W(R125B)
      CDOT(SH) = CDOT(SH) + W(R126F) - W(R126B) + W(R127F) & 
	    - W(R127B) - W(R133F) + W(R133B) & 
	    + W(R142F) - W(R142B) - W(R147F) & 
	    + W(R147B) + W(R164F) - W(R164B) & 
	    + W(R169F) - W(R169B) - W(R171F) & 
	    + W(R171F) - W(R171B) + W(R171B) & 
	    - W(R172F) + W(R172B)
      CDOT(SO2) = - W(R1F) + W(R1B) + W(R6F) & 
	    - W(R6B) - W(R9F) + W(R9B) & 
	    + W(R10F) - W(R10B) + W(R12F) & 
	    - W(R12B) + W(R13F) - W(R13B) & 
	    + W(R14F) - W(R14B) + W(R15F) & 
	    - W(R15B) - W(R23F) + W(R23B) & 
	    - W(R27F) + W(R27B) - W(R41F) & 
	    + W(R41B) - W(R45F) + W(R45B) & 
	    - W(R46F) + W(R46B) + W(R53F) & 
	    - W(R53B) - W(R60F) + W(R60B) & 
	    - W(R61F) + W(R61B) - W(R68F) & 
	    + W(R68B) - W(R69F) + W(R69B) & 
	    - W(R83F) + W(R83B) - W(R89F) & 
	    + W(R89B) - W(R102F) + W(R102B) & 
	    - W(R107F) + W(R107B) - W(R108F)
      CDOT(SO2) = CDOT(SO2) + W(R108B) - W(R116F) + W(R116B) & 
	    - W(R127F) + W(R127B) - W(R128F) & 
	    + W(R128B) - W(R137F) + W(R137B) & 
	    - W(R143F) + W(R143B) - W(R151F) & 
	    + W(R151B) + W(R152F) - W(R152B) & 
	    + W(R153F) - W(R153B) - W(R155F) & 
	    + W(R155B) - W(R158F) + W(R158B)
      CDOT(SO) = W(R1F) - W(R1B) - W(R2F) & 
	    + W(R2B) - W(R4F) + W(R4B) & 
	    - 2 * W(R6F) + 2 * W(R6B) - W(R7F) & 
	    + W(R7B) - W(R12F) + W(R12B) & 
	    - W(R19F) + W(R19B) - W(R22F) & 
	    + W(R22B) + W(R23F) - W(R23B) & 
	    - W(R29F) + W(R29B) - W(R31F) & 
	    + W(R31B) - W(R39F) + W(R39B) & 
	    - W(R44F) + W(R44B) + W(R45F) & 
	    - W(R45B) - W(R51F) + W(R51B) & 
	    - W(R58F) + W(R58B) - W(R66F) & 
	    + W(R66B) - W(R81F) + W(R81B) & 
	    - W(R88F) + W(R88B) - W(R92F) & 
	    + W(R92B) - W(R99F) + W(R99B) & 
	    - W(R101F) + W(R101B) - W(R110F)
      CDOT(SO) = CDOT(SO) + W(R110B) - W(R113F) + W(R113B) & 
	    - W(R123F) + W(R123B) - W(R124F) & 
	    + W(R124B) - W(R135F) + W(R135B) & 
	    - W(R146F) + W(R146B) - W(R175F) & 
	    + W(R175B)
      CDOT(SOH) = W(R1F) - W(R1B) + W(R2F) & 
	    - W(R2B) - W(R3F) + W(R3B) & 
	    + 2 * W(R4F) - 2 * W(R4B) + W(R7F) & 
	    - W(R7B) - W(R8F) + W(R8B) & 
	    + 2 * W(R11F) - 2 * W(R11B) + W(R12F) & 
	    - W(R12B) - W(R13F) + W(R13B) & 
	    + 2 * W(R16F) - 2 * W(R16B) + W(R17F) & 
	    - W(R17B) + W(R19F) - W(R19B) & 
	    - W(R20F) + W(R20B) - W(R21F) & 
	    + W(R21B) + W(R24F) - W(R24B) & 
	    - W(R25F) + W(R25B) + W(R29F) & 
	    - W(R29B) - W(R30F) + W(R30B) & 
	    + W(R32F) - W(R32B) + W(R39F) & 
	    - W(R39B) - W(R40F) + W(R40B) & 
	    + W(R46F) - W(R46B) + W(R47F)
      CDOT(SOH) = CDOT(SOH) - W(R47B) + W(R51F) - W(R51B) & 
	    - W(R52F) + W(R52B) + W(R57F) & 
	    - W(R57B) + W(R58F) - W(R58B) & 
	    - W(R59F) + W(R59B) + W(R65F) & 
	    - W(R65B) + W(R66F) - W(R66B) & 
	    - W(R67F) + W(R67B) - W(R75F) & 
	    + W(R75B) - W(R76F) + W(R76B) & 
	    + W(R81F) - W(R81B) - W(R82F) & 
	    + W(R82B) - W(R97F) + W(R97B) & 
	    - W(R100F) + W(R100B) + W(R101F) & 
	    - W(R101B) - W(R111F) + W(R111B) & 
	    - W(R114F) + W(R114B) + W(R116F) & 
	    - W(R116B) + W(R117F) - W(R117B) & 
	    - W(R125F) + W(R125B) + W(R127F) & 
	    - W(R127B) - W(R132F) + W(R132B)
      CDOT(SOH) = CDOT(SOH) + W(R135F) - W(R135B) + W(R141F) & 
	    - W(R141B) - W(R144F) + W(R144B) & 
	    + W(R146F) - W(R146B) + W(R157F) & 
	    - W(R157B) + W(R159F) - W(R159B) & 
	    + W(R160F) - W(R160B) - W(R165F) & 
	    + W(R165B) + W(R168F) - W(R168B) & 
	    - W(R169F) + W(R169B) - W(R170F) & 
	    + W(R170F) - W(R170B) + W(R170B) & 
	    + W(R172F) - W(R172B) + W(R173F) & 
	    - W(R173B) + W(R174F) - W(R174B) & 
	    + 2 * W(R175F) - 2 * W(R175B)
      CDOT(SH2) = - W(R2F) + W(R2B) - W(R3F) & 
	    + W(R3B) - W(R5F) + W(R5B) & 
	    + W(R10F) - W(R10B) + W(R18F) & 
	    - W(R18B) + W(R28F) - W(R28B) & 
	    + W(R33F) - W(R33B) + W(R37F) & 
	    - W(R37B) + W(R38F) - W(R38B) & 
	    + W(R50F) - W(R50B) + W(R56F) & 
	    - W(R56B) + W(R80F) - W(R80B) & 
	    + W(R87F) - W(R87B) + W(R93F) & 
	    - W(R93B) + W(R96F) - W(R96B) & 
	    + W(R103F) - W(R103B) - W(R115F) & 
	    + W(R115B) + W(R118F) - W(R118B) & 
	    + W(R123F) - W(R123B) - W(R126F) & 
	    + W(R126B) + W(R133F) - W(R133B) & 
	    + W(R147F) - W(R147B) + W(R167F)
      CDOT(SH2) = CDOT(SH2) - W(R167B) + W(R171F) - W(R171B) & 
	    + W(R172F) - W(R172B)
      CDOT(SH2O) = W(R3F) - W(R3B) - W(R4F) & 
	    + W(R4B) + W(R8F) - W(R8B) & 
	    + W(R13F) - W(R13B) + W(R17F) & 
	    - W(R17B) + W(R20F) - W(R20B) & 
	    + W(R21F) - W(R21B) + W(R30F) & 
	    - W(R30B) + W(R40F) - W(R40B) & 
	    + W(R52F) - W(R52B) + W(R59F) & 
	    - W(R59B) + W(R67F) - W(R67B) & 
	    + W(R75F) - W(R75B) + W(R76F) & 
	    - W(R76B) + W(R79F) - W(R79B) & 
	    + W(R82F) - W(R82B) + W(R97F) & 
	    - W(R97B) + W(R100F) - W(R100B) & 
	    - W(R120F) + W(R120F) - W(R120B) & 
	    + W(R120B) + W(R128F) - W(R128B) & 
	    + W(R132F) - W(R132B) + W(R144F)
      CDOT(SH2O) = CDOT(SH2O) - W(R144B) + W(R166F) - W(R166B) & 
	    + W(R169F) - W(R169B) + W(R170F) & 
	    - W(R170B)
      CDOT(SHO2) = W(R9F) - W(R9B) - W(R10F) & 
	    + W(R10B) - W(R11F) + W(R11B) & 
	    - W(R12F) + W(R12B) - W(R13F) & 
	    + W(R13B) - 2 * W(R14F) + 2 * W(R14B) & 
	    - 2 * W(R15F) + 2 * W(R15B) + W(R18F) & 
	    - W(R18B) + W(R19F) - W(R19B) & 
	    + W(R20F) - W(R20B) + W(R21F) & 
	    - W(R21B) - W(R24F) + W(R24B) & 
	    + W(R27F) - W(R27B) - W(R32F) & 
	    + W(R32B) + W(R41F) - W(R41B) & 
	    - W(R42F) + W(R42B) - W(R47F) & 
	    + W(R47B) - W(R53F) + W(R53B) & 
	    - W(R54F) + W(R54B) + W(R60F) & 
	    - W(R60B) + W(R61F) - W(R61B) & 
	    - W(R62F) + W(R62B) + W(R68F)
      CDOT(SHO2) = CDOT(SHO2) - W(R68B) + W(R69F) - W(R69B) & 
	    - W(R70F) + W(R70B) + W(R83F) & 
	    - W(R83B) - W(R84F) + W(R84B) & 
	    + W(R89F) - W(R89B) + W(R102F) & 
	    - W(R102B) + W(R104F) - W(R104B) & 
	    + W(R108F) - W(R108B) - W(R117F) & 
	    + W(R117B) - W(R136F) + W(R136B) & 
	    + W(R137F) - W(R137B) - W(R141F) & 
	    + W(R141B) + W(R143F) - W(R143B) & 
	    - W(R145F) + W(R145B) + W(R155F) & 
	    - W(R155B) - W(R174F) + W(R174B)
      CDOT(SH2O2) = W(R14F) - W(R14B) + W(R15F) & 
	    - W(R15B) - W(R16F) + W(R16B) & 
	    - W(R17F) + W(R17B) - W(R18F) & 
	    + W(R18B) - W(R19F) + W(R19B) & 
	    - W(R20F) + W(R20B) - W(R21F) & 
	    + W(R21B) + W(R42F) - W(R42B) & 
	    + W(R54F) - W(R54B) + W(R62F) & 
	    - W(R62B) + W(R70F) - W(R70B) & 
	    + W(R84F) - W(R84B) - W(R104F) & 
	    + W(R104B) + W(R136F) - W(R136B) & 
	    + W(R145F) - W(R145B) + W(R174F) & 
	    - W(R174B)
      CDOT(SCO) = - W(R22F) + W(R22B) - W(R23F) & 
	    + W(R23B) - W(R24F) + W(R24B) & 
	    - W(R25F) + W(R25B) + W(R26F) & 
	    - W(R26B) + W(R27F) - W(R27B) & 
	    + W(R28F) - W(R28B) + W(R29F) & 
	    - W(R29B) + W(R30F) - W(R30B) & 
	    + 2 * W(R33F) - 2 * W(R33B) + W(R34F) & 
	    - W(R34B) + W(R35F) - W(R35B) & 
	    + W(R37F) - W(R37B) - W(R71F) & 
	    + W(R71B) + W(R91F) - W(R91B) & 
	    + W(R110F) - W(R110B) + W(R111F) & 
	    - W(R111B) - W(R121F) + W(R121F) & 
	    - W(R121B) + W(R121B) + W(R123F) & 
	    - W(R123B) + W(R127F) - W(R127B) & 
	    + W(R128F) - W(R128B) + W(R129F)
      CDOT(SCO) = CDOT(SCO) - W(R129B) + W(R149F) - W(R149B) & 
	    + W(R162F) - W(R162B) + W(R166F) & 
	    - W(R166B) + W(R170F) - W(R170B) & 
	    + W(R172F) - W(R172B) + W(R173F) & 
	    - W(R173B) + W(R174F) - W(R174B) & 
	    + W(R175F) - W(R175B)
      CDOT(SCO2) = W(R22F) - W(R22B) + W(R23F) & 
	    - W(R23B) + W(R24F) - W(R24B) & 
	    + W(R25F) - W(R25B) + W(R31F) & 
	    - W(R31B) + W(R32F) - W(R32B) & 
	    + W(R71F) - W(R71B) - W(R122F) & 
	    + W(R122F) - W(R122B) + W(R122B) & 
	    - W(R129F) + W(R129B) + W(R150F) & 
	    - W(R150B) + W(R163F) - W(R163B) & 
	    + W(R167F) - W(R167B) + W(R169F) & 
	    - W(R169B) + W(R171F) - W(R171B)
      CDOT(SHCO) = - W(R26F) + W(R26B) - W(R27F) & 
	    + W(R27B) - W(R28F) + W(R28B) & 
	    - W(R29F) + W(R29B) - W(R30F) & 
	    + W(R30B) - W(R31F) + W(R31B) & 
	    - W(R32F) + W(R32B) - 2 * W(R33F) & 
	    + 2 * W(R33B) - W(R34F) + W(R34B) & 
	    - 2 * W(R35F) + 2 * W(R35B) + W(R36F) & 
	    - W(R36B) + W(R38F) - W(R38B) & 
	    + W(R39F) - W(R39B) + W(R40F) & 
	    - W(R40B) + W(R41F) - W(R41B) & 
	    + W(R42F) - W(R42B) + W(R43F) & 
	    - W(R43B) - W(R63F) + W(R63B) & 
	    - W(R91F) + W(R91B) + W(R99F) & 
	    - W(R99B) + W(R107F) - W(R107B) & 
	    + W(R113F) - W(R113B) + W(R116F)
      CDOT(SHCO) = CDOT(SHCO) - W(R116B) + W(R124F) - W(R124B) & 
	    + W(R130F) - W(R130B) + W(R140F) & 
	    - W(R140B) + W(R168F) - W(R168B)
      CDOT(SCH3) = - W(R34F) + W(R34B) - W(R43F) & 
	    + W(R43B) - W(R44F) + W(R44B) & 
	    - W(R45F) + W(R45B) - W(R46F) & 
	    + W(R46B) - W(R47F) + W(R47B) & 
	    - 2 * W(R48F) + 2 * W(R48B) - W(R49F) & 
	    + W(R49B) + W(R50F) - W(R50B) & 
	    + W(R51F) - W(R51B) + W(R52F) & 
	    - W(R52B) - W(R53F) + W(R53B) & 
	    + W(R54F) - W(R54B) + W(R57F) & 
	    - W(R57B) + W(R65F) - W(R65B) & 
	    + W(R71F) - W(R71B) - 2 * W(R72F) & 
	    + 2 * W(R72B) + 2 * W(R73F) - 2 * W(R73B) & 
	    + 2 * W(R74F) - 2 * W(R74B) - W(R75F) & 
	    + W(R75B) - W(R76F) + W(R76B) & 
	    - W(R77F) + W(R77B) - W(R78F)
      CDOT(SCH3) = CDOT(SCH3) + W(R78B) - W(R85F) + W(R85B) & 
	    + W(R88F) - W(R88B) - W(R98F) & 
	    + W(R98B) + W(R99F) - W(R99B) & 
	    - W(R105F) + W(R105B) + W(R111F) & 
	    - W(R111B) + W(R112F) - W(R112B) & 
	    + W(R115F) - W(R115B) + W(R126F) & 
	    - W(R126B) + W(R130F) - W(R130B) & 
	    + W(R131F) - W(R131B) - W(R134F) & 
	    + W(R134B) + W(R138F) - W(R138B) & 
	    - W(R148F) + W(R148B) + W(R150F) & 
	    - W(R150B) - W(R173F) + W(R173B)
      CDOT(SCH4) = W(R34F) - W(R34B) + W(R43F) & 
	    - W(R43B) + W(R49F) - W(R49B) & 
	    - W(R50F) + W(R50B) - W(R51F) & 
	    + W(R51B) - W(R52F) + W(R52B) & 
	    + W(R53F) - W(R53B) - W(R54F) & 
	    + W(R54B) - W(R73F) + W(R73B) & 
	    - W(R74F) + W(R74B) + W(R85F) & 
	    - W(R85B) + W(R98F) - W(R98B) & 
	    + W(R105F) - W(R105B) + W(R134F) & 
	    - W(R134B) + W(R148F) - W(R148B) & 
	    + W(R173F) - W(R173B)
      CDOT(SCH2O) = W(R35F) - W(R35B) - W(R36F) & 
	    + W(R36B) - W(R37F) + W(R37B) & 
	    - W(R38F) + W(R38B) - W(R39F) & 
	    + W(R39B) - W(R40F) + W(R40B) & 
	    - W(R41F) + W(R41B) - W(R42F) & 
	    + W(R42B) - W(R43F) + W(R43B) & 
	    + W(R44F) - W(R44B) + W(R46F) & 
	    - W(R46B) + W(R55F) - W(R55B) & 
	    + W(R56F) - W(R56B) + W(R58F) & 
	    - W(R58B) + W(R59F) - W(R59B) & 
	    + W(R60F) - W(R60B) + W(R61F) & 
	    - W(R61B) + W(R62F) - W(R62B) & 
	    + 2 * W(R63F) - 2 * W(R63B) + W(R64F) & 
	    - W(R64B) + W(R66F) - W(R66B) & 
	    + W(R67F) - W(R67B) + W(R68F)
      CDOT(SCH2O) = CDOT(SCH2O) - W(R68B) + W(R69F) - W(R69B) & 
	    + W(R70F) - W(R70B) + W(R88F) & 
	    - W(R88B) + W(R107F) - W(R107B) & 
	    + W(R114F) - W(R114B) + W(R117F) & 
	    - W(R117B) + W(R125F) - W(R125B) & 
	    + W(R129F) - W(R129B) + W(R138F) & 
	    - W(R138B) + W(R139F) - W(R139B) & 
	    - W(R140F) + W(R140B) + W(R154F) & 
	    - W(R154B) + 2 * W(R157F) - 2 * W(R157B) & 
	    - W(R165F) + W(R165B)
      CDOT(SCH3O) = W(R45F) - W(R45B) + W(R47F) & 
	    - W(R47B) - W(R64F) + W(R64B) & 
	    - W(R65F) + W(R65B) - W(R66F) & 
	    + W(R66B) - W(R67F) + W(R67B) & 
	    - W(R68F) + W(R68B) - W(R69F) & 
	    + W(R69B) - W(R70F) + W(R70B) & 
	    - W(R71F) + W(R71B) - W(R79F) & 
	    + W(R79B) + W(R131F) - W(R131B) & 
	    - W(R139F) + W(R139B) + W(R149F) & 
	    - W(R149B) + W(R154F) - W(R154B)
      CDOT(SC2H6) = W(R48F) - W(R48B) - W(R80F) & 
	    + W(R80B) - W(R81F) + W(R81B) & 
	    - W(R82F) + W(R82B) - W(R83F) & 
	    + W(R83B) - W(R84F) + W(R84B) & 
	    - W(R85F) + W(R85B) + W(R86F) & 
	    - W(R86B) + W(R90F) - W(R90B) & 
	    + W(R91F) - W(R91B)
      CDOT(SCH2OH) = - W(R55F) + W(R55B) - W(R56F) & 
	    + W(R56B) - W(R57F) + W(R57B) & 
	    - W(R58F) + W(R58B) - W(R59F) & 
	    + W(R59B) - W(R60F) + W(R60B) & 
	    - W(R61F) + W(R61B) - W(R62F) & 
	    + W(R62B) - W(R63F) + W(R63B) & 
	    + W(R163F) - W(R163B)
      CDOT(SC2H5) = W(R72F) - W(R72B) + W(R80F) & 
	    - W(R80B) + W(R81F) - W(R81B) & 
	    + W(R82F) - W(R82B) + W(R83F) & 
	    - W(R83B) + W(R84F) - W(R84B) & 
	    + W(R85F) - W(R85B) - W(R86F) & 
	    + W(R86B) - W(R87F) + W(R87B) & 
	    - W(R88F) + W(R88B) - W(R89F) & 
	    + W(R89B) - 2 * W(R90F) + 2 * W(R90B) & 
	    - W(R91F) + W(R91B) - W(R92F) & 
	    + W(R92B) + W(R94F) - W(R94B)
      CDOT(SCH2) = - W(R73F) + W(R73B) + W(R75F) & 
	    - W(R75B) - W(R77F) + W(R77B) & 
	    + W(R110F) - W(R110B) - W(R112F) & 
	    + W(R112B) - W(R113F) + W(R113B) & 
	    - W(R114F) + W(R114B) - W(R115F) & 
	    + W(R115B) - W(R116F) + W(R116B) & 
	    - W(R117F) + W(R117B) - 2 * W(R118F) & 
	    + 2 * W(R118B) + W(R119F) - W(R119B) & 
	    + W(R120F) - W(R120B) + W(R121F) & 
	    - W(R121B) + W(R122F) - W(R122B)
      CDOT(SCH2X) = - W(R74F) + W(R74B) + W(R76F) & 
	    - W(R76B) - W(R78F) + W(R78B) & 
	    + W(R79F) - W(R79B) - W(R119F) & 
	    + W(R119B) - W(R120F) + W(R120B) & 
	    - W(R121F) + W(R121B) - W(R122F) & 
	    + W(R122B) - W(R123F) + W(R123B) & 
	    - W(R124F) + W(R124B) - W(R125F) & 
	    + W(R125B) - W(R126F) + W(R126B) & 
	    - W(R127F) + W(R127B) - W(R128F) & 
	    + W(R128B) - W(R129F) + W(R129B)
      CDOT(SC2H4) = W(R77F) - W(R77B) + W(R78F) & 
	    - W(R78B) + W(R87F) - W(R87B) & 
	    + W(R89F) - W(R89B) + W(R90F) & 
	    - W(R90B) - W(R93F) + W(R93B) & 
	    - W(R94F) + W(R94B) + W(R95F) & 
	    - W(R95B) - W(R96F) + W(R96B) & 
	    - W(R97F) + W(R97B) - W(R98F) & 
	    + W(R98B) - W(R99F) + W(R99B) & 
	    - W(R101F) + W(R101B) - W(R102F) & 
	    + W(R102B) + W(R104F) - W(R104B) & 
	    + W(R106F) - W(R106B)
      CDOT(SCH3HCO) = W(R92F) - W(R92B) - W(R130F) & 
	    + W(R130B)
      CDOT(SC2H2) = W(R93F) - W(R93B) + W(R100F) & 
	    - W(R100B) + W(R103F) - W(R103B) & 
	    + W(R105F) - W(R105B) + W(R106F) & 
	    - W(R106B) + W(R108F) - W(R108B) & 
	    - W(R109F) + W(R109B) - W(R110F) & 
	    + W(R110B) - W(R111F) + W(R111B) & 
	    + W(R118F) - W(R118B)
      CDOT(SC2H3) = - W(R95F) + W(R95B) + W(R96F) & 
	    - W(R96B) + W(R97F) - W(R97B) & 
	    + W(R98F) - W(R98B) - W(R100F) & 
	    + W(R100B) + W(R101F) - W(R101B) & 
	    + W(R102F) - W(R102B) - W(R103F) & 
	    + W(R103B) - W(R104F) + W(R104B) & 
	    - W(R105F) + W(R105B) - 2 * W(R106F) & 
	    + 2 * W(R106B) - W(R107F) + W(R107B) & 
	    - W(R108F) + W(R108B) + W(R109F) & 
	    - W(R109B)
      CDOT(SCH3OCH3) = - W(R131F) + W(R131B) - W(R132F) & 
	    + W(R132B) - W(R133F) + W(R133B) & 
	    - W(R134F) + W(R134B) - W(R135F) & 
	    + W(R135B) - W(R136F) + W(R136B) & 
	    - W(R137F) + W(R137B) + W(R139F) & 
	    - W(R139B) + W(R140F) - W(R140B)
      CDOT(SCH3OCH2) = W(R132F) - W(R132B) + W(R133F) & 
	    - W(R133B) + W(R134F) - W(R134B) & 
	    + W(R135F) - W(R135B) + W(R136F) & 
	    - W(R136B) + W(R137F) - W(R137B) & 
	    - W(R138F) + W(R138B) - W(R139F) & 
	    + W(R139B) - W(R140F) + W(R140B) & 
	    - W(R141F) + W(R141B) - W(R151F) & 
	    + W(R151B)
      CDOT(SCH3OCH2O) = W(R141F) - W(R141B) - W(R142F) & 
	    + W(R142B) + 2 * W(R152F) - 2 * W(R152B) & 
	    - W(R154F) + W(R154B) - W(R155F) & 
	    + W(R155B)
      CDOT(SCH3OCHO) = W(R142F) - W(R142B) - W(R143F) & 
	    + W(R143B) - W(R144F) + W(R144B) & 
	    - W(R145F) + W(R145B) - W(R146F) & 
	    + W(R146B) - W(R147F) + W(R147B) & 
	    - W(R148F) + W(R148B) + W(R153F) & 
	    - W(R153B) + W(R155F) - W(R155B)
      CDOT(SCH3OCO) = W(R143F) - W(R143B) + W(R144F) & 
	    - W(R144B) + W(R145F) - W(R145B) & 
	    + W(R146F) - W(R146B) + W(R147F) & 
	    - W(R147B) + W(R148F) - W(R148B) & 
	    - W(R149F) + W(R149B) - W(R150F) & 
	    + W(R150B)
      CDOT(SRO2) = W(R151F) - W(R151B) - 2 * W(R152F) & 
	    + 2 * W(R152B) - 2 * W(R153F) + 2 * W(R153B) & 
	    - W(R156F) + W(R156B)
      CDOT(SROH) = W(R153F) - W(R153B)
      CDOT(SQOOH) = W(R156F) - W(R156B) - W(R157F) & 
	    + W(R157B) - W(R158F) + W(R158B)
      CDOT(SO2QOOH) = W(R158F) - W(R158B) - W(R159F) & 
	    + W(R159B)
      CDOT(SHO2QHO) = W(R159F) - W(R159B) - W(R160F) & 
	    + W(R160B)
      CDOT(SOCH2OCHO) = W(R160F) - W(R160B) - W(R161F) & 
	    + W(R161B)
      CDOT(SHOCH2OCO) = W(R161F) - W(R161B) - W(R162F) & 
	    + W(R162B) - W(R163F) + W(R163B)
      CDOT(SHOCH2O) = W(R162F) - W(R162B) - W(R164F) & 
	    + W(R164B) + W(R165F) - W(R165B)
      CDOT(SHCOOH) = W(R164F) - W(R164B) - W(R166F) & 
	    + W(R166B) - W(R167F) + W(R167B) & 
	    - W(R168F) + W(R168B) - W(R169F) & 
	    + W(R169B) - W(R170F) + W(R170B) & 
	    - W(R171F) + W(R171B) - W(R172F) & 
	    + W(R172B) - W(R173F) + W(R173B) & 
	    - W(R174F) + W(R174B) - W(R175F) & 
	    + W(R175B)
      END


      DOUBLE PRECISION FUNCTION GETLINDRATECOEFF( TEMP, PRESSURE, &
	   K0, KINF, FC, CONCIN )

      IMPLICIT NONE
      INTEGER ITROE
      DOUBLE PRECISION TEMP, PRESSURE, K0, KINF, FC
      DOUBLE PRECISION R
      PARAMETER (R = 8314.34)
!     [J / kmole K]
      DOUBLE PRECISION NTMP,CCOEFF,DCOEFF,LGKNULL
      DOUBLE PRECISION KL
      DOUBLE PRECISION F
      DOUBLE PRECISION CONC, CONCIN

      ITROE = 1
      IF (CONCIN.GT.0.0) THEN
        CONC = CONCIN
      ELSE
        CONC = PRESSURE / ( R * TEMP )
      END IF
      NTMP = 0.75 - 1.27 * DLOG10( FC )
      IF (ITROE.EQ.1) THEN
        CCOEFF = - 0.4 - 0.67 * DLOG10( FC )
        DCOEFF = 0.14
        K0 = K0 * CONC / MAX(KINF, 1.0D-60)
        LGKNULL = DLOG10( K0 )
        F=(LGKNULL+CCOEFF)/(NTMP-DCOEFF*(LGKNULL+CCOEFF))
        F = FC**(1.0 / ( F * F + 1.0 ))
        GETLINDRATECOEFF = KINF * F * K0 / ( 1.0 + K0 )
      ELSE
        K0 = K0 * CONC / KINF
        KL = K0 / ( 1.0 + K0 )
        F = DLOG10( K0 ) / NTMP
        F = FC ** ( 1.0 / ( F * F + 1.0 ) )
        GETLINDRATECOEFF = KINF * F * KL
      END IF
      END


      DOUBLE PRECISION FUNCTION CTCHZERO( X )

      IMPLICIT NONE
      DOUBLE PRECISION X
      IF ( X.EQ.0.0D0 ) THEN
         X = 1.0D-20
      END IF
      CTCHZERO = X
      END


      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )
!------------------------------------------------------------------
!	SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0
!------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RAD, A, B, C

      IF ( A.EQ.0.0D0 ) THEN
        IF ( B.EQ.0.0D0 ) THEN
          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '
          SOLVQDRT = 0.0
        ELSE
          SOLVQDRT = C / B
        END IF
      ELSE

        B = B / A
        C = C / A

        RAD = 0.25D0 * B * B + C
        IF ( RAD.GE.0.0D0 ) THEN
          RAD = DSQRT( RAD ) + 0.5 * B
          IF ( RAD.NE.0.0D0 ) THEN
            SOLVQDRT = C / RAD
          ELSE
            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )
          ENDIF
        ELSE
          SOLVQDRT = 0.0D0
        ENDIF
      END IF
      END


      SUBROUTINE GETMOLARMASS( MM )
!------------------------------------------------------------------
!	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE
!------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION MM(39)
      INCLUDE 'DME39F90.h'

      MM(SN2) =  2.80200000e+01
      MM(SH) =  1.00800000e+00
      MM(SO2) =  3.20000000e+01
      MM(SO) =  1.60000000e+01
      MM(SOH) =  1.70080000e+01
      MM(SH2) =  2.01600000e+00
      MM(SH2O) =  1.80160000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SH2O2) =  3.40160000e+01
      MM(SCO) =  2.80100000e+01
      MM(SCO2) =  4.40100000e+01
      MM(SHCO) =  2.90180000e+01
      MM(SCH3) =  1.50340000e+01
      MM(SCH4) =  1.60420000e+01
      MM(SCH2O) =  3.00260000e+01
      MM(SCH3O) =  3.10340000e+01
      MM(SC2H6) =  3.00680000e+01
      MM(SCH2OH) =  3.10340000e+01
      MM(SC2H5) =  2.90600000e+01
      MM(SCH2) =  1.40260000e+01
      MM(SCH2X) =  1.40260000e+01
      MM(SC2H4) =  2.80520000e+01
      MM(SCH3HCO) =  4.40520000e+01
      MM(SC2H2) =  2.60360000e+01
      MM(SC2H3) =  2.70440000e+01
      MM(SCH3OCH3) =  4.60680000e+01
      MM(SCH3OCH2) =  4.50600000e+01
      MM(SCH3OCH2O) =  6.10600000e+01
      MM(SCH3OCHO) =  6.00520000e+01
      MM(SCH3OCO) =  5.90440000e+01
      MM(SRO2) =  7.70600000e+01
      MM(SROH) =  6.20680000e+01
      MM(SQOOH) =  7.70600000e+01
      MM(SO2QOOH) =  1.09060000e+02
      MM(SHO2QHO) =  9.20520000e+01
      MM(SOCH2OCHO) =  7.50440000e+01
      MM(SHOCH2OCO) =  7.50440000e+01
      MM(SHOCH2O) =  4.70340000e+01
      MM(SHCOOH) =  4.60260000e+01

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(39)
      INCLUDE 'DME39F90.h'

      NAMES(SN2)='N2                  '
      NAMES(SH)='H                   '
      NAMES(SO2)='O2                  '
      NAMES(SO)='O                   '
      NAMES(SOH)='OH                  '
      NAMES(SH2)='H2                  '
      NAMES(SH2O)='H2O                 '
      NAMES(SHO2)='HO2                 '
      NAMES(SH2O2)='H2O2                '
      NAMES(SCO)='CO                  '
      NAMES(SCO2)='CO2                 '
      NAMES(SHCO)='HCO                 '
      NAMES(SCH3)='CH3                 '
      NAMES(SCH4)='CH4                 '
      NAMES(SCH2O)='CH2O                '
      NAMES(SCH3O)='CH3O                '
      NAMES(SC2H6)='C2H6                '
      NAMES(SCH2OH)='CH2OH               '
      NAMES(SC2H5)='C2H5                '
      NAMES(SCH2)='CH2                 '
      NAMES(SCH2X)='CH2X          '
      NAMES(SC2H4)='C2H4                '
      NAMES(SCH3HCO)='CH3HCO              '
      NAMES(SC2H2)='C2H2                '
      NAMES(SC2H3)='C2H3                '
      NAMES(SCH3OCH3)='CH3OCH3             '
      NAMES(SCH3OCH2)='CH3OCH2             '
      NAMES(SCH3OCH2O)='CH3OCH2O            '
      NAMES(SCH3OCHO)='CH3OCHO             '
      NAMES(SCH3OCO)='CH3OCO              '
      NAMES(SRO2)='RO2           '
      NAMES(SROH)='ROH           '
      NAMES(SQOOH)='QOOH          '
      NAMES(SO2QOOH)='O2QOOH        '
      NAMES(SHO2QHO)='HO2QHO          '
      NAMES(SOCH2OCHO)='OCH2OCHO            '
      NAMES(SHOCH2OCO)='HOCH2OCO            '
      NAMES(SHOCH2O)='HOCH2O              '
      NAMES(SHCOOH)='HCOOH               '

      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'DME39F90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 350

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'DME39F90.h'

      NSPECIES_NONS = 39
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'DME39F90.h'
      real(DP) :: MUCOEFF(39)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SH2O2) =  1.30193418e-06
      MUCOEFF(SCO) =  1.06039632e-06
      MUCOEFF(SCO2) =  1.25056029e-06
      MUCOEFF(SHCO) =  1.11568663e-06
      MUCOEFF(SCH3) =  7.16749611e-07
      MUCOEFF(SCH4) =  7.61887935e-07
      MUCOEFF(SCH2O) =  1.13489904e-06
      MUCOEFF(SCH3O) =  1.09210283e-06
      MUCOEFF(SC2H6) =  7.73519281e-07
      MUCOEFF(SCH2OH) =  1.09210283e-06
      MUCOEFF(SC2H5) =  7.60443019e-07
      MUCOEFF(SCH2) =  6.92304430e-07
      MUCOEFF(SCH2X) =  6.92304430e-07
      MUCOEFF(SC2H4) =  1.15674186e-06
      MUCOEFF(SCH3HCO) =  1.12408509e-06
      MUCOEFF(SC2H2) =  9.83705677e-07
      MUCOEFF(SC2H3) =  1.00256724e-06
      MUCOEFF(SCH3OCH3) =  8.47347231e-07
      MUCOEFF(SCH3OCH2) =  8.38025684e-07
      MUCOEFF(SCH3OCH2O) =  8.82360133e-07
      MUCOEFF(SCH3OCHO) =  9.32832661e-07
      MUCOEFF(SCH3OCO) =  9.24970520e-07
      MUCOEFF(SRO2) =  1.09591340e-06
      MUCOEFF(SROH) =  9.83548117e-07
      MUCOEFF(SQOOH) =  1.09591340e-06
      MUCOEFF(SO2QOOH) =  1.30375048e-06
      MUCOEFF(SHO2QHO) =  1.19778357e-06
      MUCOEFF(SOCH2OCHO) =  1.08148306e-06
      MUCOEFF(SHOCH2OCO) =  1.08148306e-06
      MUCOEFF(SHOCH2O) =  1.39234784e-06
      MUCOEFF(SHCOOH) =  9.31154667e-07

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'DME39F90.h'
      real(DP) :: KOVEREPS(39)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SH2O2) =  9.31098696e-03
      KOVEREPS(SCO) =  1.01936799e-02
      KOVEREPS(SCO2) =  4.09836066e-03
      KOVEREPS(SHCO) =  2.00803213e-03
      KOVEREPS(SCH3) =  6.94444444e-03
      KOVEREPS(SCH4) =  7.07213579e-03
      KOVEREPS(SCH2O) =  2.00803213e-03
      KOVEREPS(SCH3O) =  2.39808153e-03
      KOVEREPS(SC2H6) =  4.04040404e-03
      KOVEREPS(SCH2OH) =  2.39808153e-03
      KOVEREPS(SC2H5) =  4.04040404e-03
      KOVEREPS(SCH2) =  6.94444444e-03
      KOVEREPS(SCH2X) =  6.94444444e-03
      KOVEREPS(SC2H4) =  4.19463087e-03
      KOVEREPS(SCH3HCO) =  2.29357798e-03
      KOVEREPS(SC2H2) =  3.76931775e-03
      KOVEREPS(SC2H3) =  3.76931775e-03
      KOVEREPS(SCH3OCH3) =  3.03582271e-03
      KOVEREPS(SCH3OCH2) =  3.03582271e-03
      KOVEREPS(SCH3OCH2O) =  2.12359312e-03
      KOVEREPS(SCH3OCHO) =  2.46002460e-03
      KOVEREPS(SCH3OCO) =  2.46002460e-03
      KOVEREPS(SRO2) =  3.03582271e-03
      KOVEREPS(SROH) =  3.03582271e-03
      KOVEREPS(SQOOH) =  3.03582271e-03
      KOVEREPS(SO2QOOH) =  3.03582271e-03
      KOVEREPS(SHO2QHO) =  3.03582271e-03
      KOVEREPS(SOCH2OCHO) =  3.03582271e-03
      KOVEREPS(SHOCH2OCO) =  3.03582271e-03
      KOVEREPS(SHOCH2O) =  2.07555002e-03
      KOVEREPS(SHCOOH) =  2.12494688e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'DME39F90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(350), C(39), M(20)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 39
!------------------------------------------------------------------
      implicit none
      include 'DME39F90.h'
      real(DP) :: H(39), CP(39), T

      IF (T.GT.1000.0_DP) THEN
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  2.92664000e+00_DP + T * (  7.43988500e-04_DP &
	   + T * ( -1.89492033e-07_DP + T * (  2.52426000e-11_DP &
	   + T * ( -1.35067020e-15_DP ) ) ) ) ) -9.22797700e+02_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    2.92664000e+00_DP + T * (  1.48797700e-03_DP &
	   + T * ( -5.68476100e-07_DP + T * (  1.00970400e-10_DP &
	   + T * ( -6.75335100e-15_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  2.54716300e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.69757800e+00_DP + T * (  3.06759850e-04_DP &
	   + T * ( -4.19614000e-08_DP + T * (  4.43820250e-12_DP &
	   + T * ( -2.27287000e-16_DP ) ) ) ) ) -1.23393000e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.69757800e+00_DP + T * (  6.13519700e-04_DP &
	   + T * ( -1.25884200e-07_DP + T * (  1.77528100e-11_DP &
	   + T * ( -1.13643500e-15_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.54206000e+00_DP + T * ( -1.37753100e-05_DP &
	   + T * ( -1.03426767e-09_DP + T * (  1.13776675e-12_DP &
	   + T * ( -8.73610400e-17_DP ) ) ) ) ) +  2.92308000e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.54206000e+00_DP + T * ( -2.75506200e-05_DP &
	   + T * ( -3.10280300e-09_DP + T * (  4.55106700e-12_DP &
	   + T * ( -4.36805200e-16_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  2.86472886e+00_DP + T * (  5.28252240e-04_DP &
	   + T * ( -8.63609193e-08_DP + T * (  7.63046685e-12_DP &
	   + T * ( -2.66391752e-16_DP ) ) ) ) ) +  3.68362875e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    2.86472886e+00_DP + T * (  1.05650448e-03_DP &
	   + T * ( -2.59082758e-07_DP + T * (  3.05218674e-11_DP &
	   + T * ( -1.33195876e-15_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.99142300e+00_DP + T * (  3.50032200e-04_DP &
	   + T * ( -1.87794300e-08_DP + T * ( -2.30789450e-12_DP &
	   + T * (  3.16550400e-16_DP ) ) ) ) ) -8.35034000e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.99142300e+00_DP + T * (  7.00064400e-04_DP &
	   + T * ( -5.63382900e-08_DP + T * ( -9.23157800e-12_DP &
	   + T * (  1.58275200e-15_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  2.67214600e+00_DP + T * (  1.52814650e-03_DP &
	   + T * ( -2.91008667e-07_DP + T * (  3.00249000e-11_DP &
	   + T * ( -1.27832360e-15_DP ) ) ) ) ) -2.98992100e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    2.67214600e+00_DP + T * (  3.05629300e-03_DP &
	   + T * ( -8.73026000e-07_DP + T * (  1.20099600e-10_DP &
	   + T * ( -6.39161800e-15_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.01721090e+00_DP + T * (  1.11991006e-03_DP &
	   + T * ( -2.11219383e-07_DP + T * (  2.85615925e-11_DP &
	   + T * ( -2.15817070e-15_DP ) ) ) ) ) +  1.11856713e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.01721090e+00_DP + T * (  2.23982013e-03_DP &
	   + T * ( -6.33658150e-07_DP + T * (  1.14246370e-10_DP &
	   + T * ( -1.07908535e-14_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.57316700e+00_DP + T * (  2.16806800e-03_DP &
	   + T * ( -4.91563000e-07_DP + T * (  5.87226000e-11_DP &
	   + T * ( -2.86330800e-15_DP ) ) ) ) ) -1.80069600e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.57316700e+00_DP + T * (  4.33613600e-03_DP &
	   + T * ( -1.47468900e-06_DP + T * (  2.34890400e-10_DP &
	   + T * ( -1.43165400e-14_DP ) ) ) ) )
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  3.02507800e+00_DP + T * (  7.21344500e-04_DP &
	   + T * ( -1.87694267e-07_DP + T * (  2.54645250e-11_DP &
	   + T * ( -1.38219040e-15_DP ) ) ) ) ) -1.42683500e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    3.02507800e+00_DP + T * (  1.44268900e-03_DP &
	   + T * ( -5.63082800e-07_DP + T * (  1.01858100e-10_DP &
	   + T * ( -6.91095200e-15_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  4.45362300e+00_DP + T * (  1.57008450e-03_DP &
	   + T * ( -4.26137000e-07_DP + T * (  5.98499250e-11_DP &
	   + T * ( -3.33806600e-15_DP ) ) ) ) ) -4.89669600e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    4.45362300e+00_DP + T * (  3.14016900e-03_DP &
	   + T * ( -1.27841100e-06_DP + T * (  2.39399700e-10_DP &
	   + T * ( -1.66903300e-14_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  3.55727100e+00_DP + T * (  1.67278650e-03_DP &
	   + T * ( -4.45002000e-07_DP + T * (  6.17643250e-11_DP &
	   + T * ( -3.42770200e-15_DP ) ) ) ) ) +  3.91632400e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    3.55727100e+00_DP + T * (  3.34557300e-03_DP &
	   + T * ( -1.33500600e-06_DP + T * (  2.47057300e-10_DP &
	   + T * ( -1.71385100e-14_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  2.97812060e+00_DP + T * (  2.89892600e-03_DP &
	   + T * ( -6.58526667e-07_DP + T * (  7.68244750e-11_DP &
	   + T * ( -3.58348320e-15_DP ) ) ) ) ) +  1.65095130e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    2.97812060e+00_DP + T * (  5.79785200e-03_DP &
	   + T * ( -1.97558000e-06_DP + T * (  3.07297900e-10_DP &
	   + T * ( -1.79174160e-14_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  1.68347900e+00_DP + T * (  5.11862000e-03_DP &
	   + T * ( -1.29170967e-06_DP + T * (  1.69639625e-10_DP &
	   + T * ( -9.00684600e-15_DP ) ) ) ) ) -1.00807900e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    1.68347900e+00_DP + T * (  1.02372400e-02_DP &
	   + T * ( -3.87512900e-06_DP + T * (  6.78558500e-10_DP &
	   + T * ( -4.50342300e-14_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * ( -5.81492155e-01_DP + T * (  7.15608360e-03_DP &
	   + T * ( -2.57147386e-06_DP + T * (  3.71982433e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.41272885e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	   -5.81492155e-01_DP + T * (  1.43121672e-02_DP &
	   + T * ( -7.71442158e-06_DP + T * (  1.48792973e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3O) =  2.67910679e+02_DP * ( &
	   T * (  3.77080000e+00_DP + T * (  3.93574850e-03_DP &
	   + T * ( -8.85461333e-07_DP + T * (  9.86107750e-11_DP &
	   + T * ( -4.22523200e-15_DP ) ) ) ) ) +  1.27832500e+02_DP )
      CP(SCH3O) =  2.67910679e+02_DP * ( &
	    3.77080000e+00_DP + T * (  7.87149700e-03_DP &
	   + T * ( -2.65638400e-06_DP + T * (  3.94443100e-10_DP &
	   + T * ( -2.11261600e-14_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  4.82593800e+00_DP + T * (  6.92021500e-03_DP &
	   + T * ( -1.51908633e-06_DP + T * (  1.68124175e-10_DP &
	   + T * ( -7.19632200e-15_DP ) ) ) ) ) -1.27177900e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    4.82593800e+00_DP + T * (  1.38404300e-02_DP &
	   + T * ( -4.55725900e-06_DP + T * (  6.72496700e-10_DP &
	   + T * ( -3.59816100e-14_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  3.74691030e+00_DP + T * (  4.43230605e-03_DP &
	   + T * ( -1.41935740e-06_DP + T * (  2.52201000e-10_DP &
	   + T * ( -1.89003122e-14_DP ) ) ) ) ) -3.66648240e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    3.74691030e+00_DP + T * (  8.86461210e-03_DP &
	   + T * ( -4.25807220e-06_DP + T * (  1.00880400e-09_DP &
	   + T * ( -9.45015610e-14_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  4.28788140e+00_DP + T * (  6.21694650e-03_DP &
	   + T * ( -1.47130397e-06_DP + T * (  1.76635255e-10_DP &
	   + T * ( -8.40702720e-15_DP ) ) ) ) ) +  1.20564550e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    4.28788140e+00_DP + T * (  1.24338930e-02_DP &
	   + T * ( -4.41391190e-06_DP + T * (  7.06541020e-10_DP &
	   + T * ( -4.20351360e-14_DP ) ) ) ) )
      H(SCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.87410113e+00_DP + T * (  1.82819646e-03_DP &
	   + T * ( -4.69648657e-07_DP + T * (  6.50448872e-11_DP &
	   + T * ( -3.75455134e-15_DP ) ) ) ) ) +  4.62636040e+04_DP )
      CP(SCH2) =  5.92780550e+02_DP * ( &
	    2.87410113e+00_DP + T * (  3.65639292e-03_DP &
	   + T * ( -1.40894597e-06_DP + T * (  2.60179549e-10_DP &
	   + T * ( -1.87727567e-14_DP ) ) ) ) )
      H(SCH2X) =  5.92780550e+02_DP * ( &
	   T * (  2.29203842e+00_DP + T * (  2.32794318e-03_DP &
	   + T * ( -6.70639823e-07_DP + T * (  1.04476500e-10_DP &
	   + T * ( -6.79432730e-15_DP ) ) ) ) ) +  5.09259997e+04_DP )
      CP(SCH2X) =  5.92780550e+02_DP * ( &
	    2.29203842e+00_DP + T * (  4.65588637e-03_DP &
	   + T * ( -2.01191947e-06_DP + T * (  4.17906000e-10_DP &
	   + T * ( -3.39716365e-14_DP ) ) ) ) )
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  2.03611116e+00_DP + T * (  7.32270755e-03_DP &
	   + T * ( -2.23692638e-06_DP + T * (  3.68057308e-10_DP &
	   + T * ( -2.51412122e-14_DP ) ) ) ) ) +  4.93988614e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    2.03611116e+00_DP + T * (  1.46454151e-02_DP &
	   + T * ( -6.71077915e-06_DP + T * (  1.47222923e-09_DP &
	   + T * ( -1.25706061e-13_DP ) ) ) ) )
      H(SCH3HCO) =  1.88739217e+02_DP * ( &
	   T * (  5.40411080e+00_DP + T * (  5.86152950e-03_DP &
	   + T * ( -1.40877123e-06_DP + T * (  1.70931128e-10_DP &
	   + T * ( -8.19697260e-15_DP ) ) ) ) ) -2.25931220e+04_DP )
      CP(SCH3HCO) =  1.88739217e+02_DP * ( &
	    5.40411080e+00_DP + T * (  1.17230590e-02_DP &
	   + T * ( -4.22631370e-06_DP + T * (  6.83724510e-10_DP &
	   + T * ( -4.09848630e-14_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  4.14756964e+00_DP + T * (  2.98083332e-03_DP &
	   + T * ( -7.90982840e-07_DP + T * (  1.16853043e-10_DP &
	   + T * ( -7.22470426e-15_DP ) ) ) ) ) +  2.59359992e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    4.14756964e+00_DP + T * (  5.96166664e-03_DP &
	   + T * ( -2.37294852e-06_DP + T * (  4.67412171e-10_DP &
	   + T * ( -3.61235213e-14_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.01672400e+00_DP + T * (  5.16511460e-03_DP &
	   + T * ( -1.56027450e-06_DP + T * (  2.54408220e-10_DP &
	   + T * ( -1.72521408e-14_DP ) ) ) ) ) +  3.46128739e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.01672400e+00_DP + T * (  1.03302292e-02_DP &
	   + T * ( -4.68082349e-06_DP + T * (  1.01763288e-09_DP &
	   + T * ( -8.62607041e-14_DP ) ) ) ) )
      H(SCH3OCH3) =  1.80479726e+02_DP * ( &
	   T * (  8.30815546e-01_DP + T * (  1.34586631e-02_DP &
	   + T * ( -4.62915923e-06_DP + T * (  8.68787697e-10_DP &
	   + T * ( -6.83413568e-14_DP ) ) ) ) ) -2.34120975e+04_DP )
      CP(SCH3OCH3) =  1.80479726e+02_DP * ( &
	    8.30815546e-01_DP + T * (  2.69173263e-02_DP &
	   + T * ( -1.38874777e-05_DP + T * (  3.47515079e-09_DP &
	   + T * ( -3.41706784e-13_DP ) ) ) ) )
      H(SCH3OCH2) =  1.84517088e+02_DP * ( &
	   T * (  2.88564016e+00_DP + T * (  1.00486064e-02_DP &
	   + T * ( -2.95479066e-06_DP + T * (  3.59394375e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.16039910e+03_DP )
      CP(SCH3OCH2) =  1.84517088e+02_DP * ( &
	    2.88564016e+00_DP + T * (  2.00972128e-02_DP &
	   + T * ( -8.86437198e-06_DP + T * (  1.43757750e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OCH2O) =  1.36166721e+02_DP * ( &
	   T * (  1.95639414e+00_DP + T * (  1.31815389e-02_DP &
	   + T * ( -4.20895463e-06_DP + T * (  5.51447057e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.89171112e+04_DP )
      CP(SCH3OCH2O) =  1.36166721e+02_DP * ( &
	    1.95639414e+00_DP + T * (  2.63630778e-02_DP &
	   + T * ( -1.26268639e-05_DP + T * (  2.20578823e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OCHO) =  1.38452341e+02_DP * ( &
	   T * (  1.62368619e+00_DP + T * (  1.25773482e-02_DP &
	   + T * ( -4.19104960e-06_DP + T * (  5.61335630e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -4.38326475e+04_DP )
      CP(SCH3OCHO) =  1.38452341e+02_DP * ( &
	    1.62368619e+00_DP + T * (  2.51546964e-02_DP &
	   + T * ( -1.25731488e-05_DP + T * (  2.24534252e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OCO) =  1.40816002e+02_DP * ( &
	   T * (  6.56385115e+00_DP + T * (  8.47732580e-03_DP &
	   + T * ( -3.10522179e-06_DP + T * (  4.44166255e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.21467361e+04_DP )
      CP(SCH3OCO) =  1.40816002e+02_DP * ( &
	    6.56385115e+00_DP + T * (  1.69546516e-02_DP &
	   + T * ( -9.31566538e-06_DP + T * (  1.77666502e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SRO2) =  1.07894368e+02_DP * ( &
	   T * (  6.77859461e+00_DP + T * (  1.08441866e-02_DP &
	   + T * ( -3.19700053e-06_DP + T * (  3.91917137e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.05773390e+04_DP )
      CP(SRO2) =  1.07894368e+02_DP * ( &
	    6.77859461e+00_DP + T * (  2.16883731e-02_DP &
	   + T * ( -9.59100160e-06_DP + T * (  1.56766855e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SROH) =  1.33955339e+02_DP * ( &
	   T * (  1.90477229e+00_DP + T * (  1.42141490e-02_DP &
	   + T * ( -4.44548473e-06_DP + T * (  5.75424503e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -4.51404469e+04_DP )
      CP(SROH) =  1.33955339e+02_DP * ( &
	    1.90477229e+00_DP + T * (  2.84282980e-02_DP &
	   + T * ( -1.33364542e-05_DP + T * (  2.30169801e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SQOOH) =  1.07894368e+02_DP * ( &
	   T * (  9.67966551e+00_DP + T * (  9.42218485e-03_DP &
	   + T * ( -2.89713943e-06_DP + T * (  3.67007170e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.61394620e+04_DP )
      CP(SQOOH) =  1.07894368e+02_DP * ( &
	    9.67966551e+00_DP + T * (  1.88443697e-02_DP &
	   + T * ( -8.69141828e-06_DP + T * (  1.46802868e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SO2QOOH) =  7.62363836e+01_DP * ( &
	   T * (  1.34531977e+01_DP + T * (  1.02667791e-02_DP &
	   + T * ( -3.11560479e-06_DP + T * (  3.90688325e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.55012286e+04_DP )
      CP(SO2QOOH) =  7.62363836e+01_DP * ( &
	    1.34531977e+01_DP + T * (  2.05335582e-02_DP &
	   + T * ( -9.34681436e-06_DP + T * (  1.56275330e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SHO2QHO) =  9.03222092e+01_DP * ( &
	   T * (  9.30893518e+00_DP + T * (  1.07302405e-02_DP &
	   + T * ( -3.54804393e-06_DP + T * (  4.72094540e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -5.94810865e+04_DP )
      CP(SHO2QHO) =  9.03222092e+01_DP * ( &
	    9.30893518e+00_DP + T * (  2.14604810e-02_DP &
	   + T * ( -1.06441318e-05_DP + T * (  1.88837816e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SOCH2OCHO) =  1.10792868e+02_DP * ( &
	   T * (  2.64583321e+00_DP + T * (  1.28729906e-02_DP &
	   + T * ( -4.55462327e-06_DP + T * (  6.37131345e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.97056884e+04_DP )
      CP(SOCH2OCHO) =  1.10792868e+02_DP * ( &
	    2.64583321e+00_DP + T * (  2.57459812e-02_DP &
	   + T * ( -1.36638698e-05_DP + T * (  2.54852538e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SHOCH2OCO) =  1.10792868e+02_DP * ( &
	   T * (  3.21081862e+00_DP + T * (  1.16629115e-02_DP &
	   + T * ( -4.01300580e-06_DP + T * (  5.50467563e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -4.33279323e+04_DP )
      CP(SHOCH2OCO) =  1.10792868e+02_DP * ( &
	    3.21081862e+00_DP + T * (  2.33258230e-02_DP &
	   + T * ( -1.20390174e-05_DP + T * (  2.20187025e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SHOCH2O) =  1.76772973e+02_DP * ( &
	   T * (  1.72976110e+00_DP + T * (  8.07434020e-03_DP &
	   + T * ( -2.58564033e-06_DP + T * (  3.42722605e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.29193952e+04_DP )
      CP(SHOCH2O) =  1.76772973e+02_DP * ( &
	    1.72976110e+00_DP + T * (  1.61486804e-02_DP &
	   + T * ( -7.75692098e-06_DP + T * (  1.37089042e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SHCOOH) =  1.80644418e+02_DP * ( &
	   T * (  2.71806922e+00_DP + T * (  6.15775995e-03_DP &
	   + T * ( -2.00992155e-06_DP + T * (  2.65342645e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -4.67812848e+04_DP )
      CP(SHCOOH) =  1.80644418e+02_DP * ( &
	    2.71806922e+00_DP + T * (  1.23155199e-02_DP &
	   + T * ( -6.02976465e-06_DP + T * (  1.06137058e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      ELSE
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  3.29867700e+00_DP + T * (  7.04120000e-04_DP &
	   + T * ( -1.32107400e-06_DP + T * (  1.41037875e-09_DP &
	   + T * ( -4.88971000e-13_DP ) ) ) ) ) -1.02090000e+03_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    3.29867700e+00_DP + T * (  1.40824000e-03_DP &
	   + T * ( -3.96322200e-06_DP + T * (  5.64151500e-09_DP &
	   + T * ( -2.44485500e-12_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  2.54716300e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.21293600e+00_DP + T * (  5.63743000e-04_DP &
	   + T * ( -1.91871667e-07_DP + T * (  3.28469250e-10_DP &
	   + T * ( -1.75371080e-13_DP ) ) ) ) ) -1.00524900e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.21293600e+00_DP + T * (  1.12748600e-03_DP &
	   + T * ( -5.75615000e-07_DP + T * (  1.31387700e-09_DP &
	   + T * ( -8.76855400e-13_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.94642900e+00_DP + T * ( -8.19083000e-04_DP &
	   + T * (  8.07010667e-07_DP + T * ( -4.00710750e-10_DP &
	   + T * (  7.78139200e-14_DP ) ) ) ) ) +  2.91476400e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.94642900e+00_DP + T * ( -1.63816600e-03_DP &
	   + T * (  2.42103200e-06_DP + T * ( -1.60284300e-09_DP &
	   + T * (  3.89069600e-13_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  4.12530561e+00_DP + T * ( -1.61272470e-03_DP &
	   + T * (  2.17588230e-06_DP + T * ( -1.44963411e-09_DP &
	   + T * (  4.12474758e-13_DP ) ) ) ) ) +  3.34630913e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    4.12530561e+00_DP + T * ( -3.22544939e-03_DP &
	   + T * (  6.52764691e-06_DP + T * ( -5.79853643e-09_DP &
	   + T * (  2.06237379e-12_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.29812400e+00_DP + T * (  4.12472100e-04_DP &
	   + T * ( -2.71433833e-07_DP + T * ( -2.36885850e-11_DP &
	   + T * (  8.26974400e-14_DP ) ) ) ) ) -1.01252100e+03_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.29812400e+00_DP + T * (  8.24944200e-04_DP &
	   + T * ( -8.14301500e-07_DP + T * ( -9.47543400e-11_DP &
	   + T * (  4.13487200e-13_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  3.38684200e+00_DP + T * (  1.73749100e-03_DP &
	   + T * ( -2.11823200e-06_DP + T * (  1.74214525e-09_DP &
	   + T * ( -5.01317600e-13_DP ) ) ) ) ) -3.02081100e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    3.38684200e+00_DP + T * (  3.47498200e-03_DP &
	   + T * ( -6.35469600e-06_DP + T * (  6.96858100e-09_DP &
	   + T * ( -2.50658800e-12_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.30179801e+00_DP + T * ( -2.37456025e-03_DP &
	   + T * (  7.05276303e-06_DP + T * ( -6.06909735e-09_DP &
	   + T * (  1.85845025e-12_DP ) ) ) ) ) +  2.94808040e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.30179801e+00_DP + T * ( -4.74912051e-03_DP &
	   + T * (  2.11582891e-05_DP + T * ( -2.42763894e-08_DP &
	   + T * (  9.29225124e-12_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  3.38875400e+00_DP + T * (  3.28461300e-03_DP &
	   + T * ( -4.95004333e-08_DP + T * ( -1.15645150e-09_DP &
	   + T * (  4.94303000e-13_DP ) ) ) ) ) -1.76631500e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    3.38875400e+00_DP + T * (  6.56922600e-03_DP &
	   + T * ( -1.48501300e-07_DP + T * ( -4.62580600e-09_DP &
	   + T * (  2.47151500e-12_DP ) ) ) ) )
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  3.26245200e+00_DP + T * (  7.55970500e-04_DP &
	   + T * ( -1.29391833e-06_DP + T * (  1.39548600e-09_DP &
	   + T * ( -4.94990200e-13_DP ) ) ) ) ) -1.43105400e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    3.26245200e+00_DP + T * (  1.51194100e-03_DP &
	   + T * ( -3.88175500e-06_DP + T * (  5.58194400e-09_DP &
	   + T * ( -2.47495100e-12_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  2.27572500e+00_DP + T * (  4.96103600e-03_DP &
	   + T * ( -3.46970333e-06_DP + T * (  1.71667175e-09_DP &
	   + T * ( -4.23456000e-13_DP ) ) ) ) ) -4.83731400e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    2.27572500e+00_DP + T * (  9.92207200e-03_DP &
	   + T * ( -1.04091100e-05_DP + T * (  6.86668700e-09_DP &
	   + T * ( -2.11728000e-12_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  2.89833000e+00_DP + T * (  3.09957350e-03_DP &
	   + T * ( -3.20769467e-06_DP + T * (  2.72456250e-09_DP &
	   + T * ( -9.14977000e-13_DP ) ) ) ) ) +  4.15992200e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    2.89833000e+00_DP + T * (  6.19914700e-03_DP &
	   + T * ( -9.62308400e-06_DP + T * (  1.08982500e-08_DP &
	   + T * ( -4.57488500e-12_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  3.65717970e+00_DP + T * (  1.06329895e-03_DP &
	   + T * (  1.81946277e-06_DP + T * ( -1.65452507e-09_DP &
	   + T * (  4.93141480e-13_DP ) ) ) ) ) +  1.64227160e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    3.65717970e+00_DP + T * (  2.12659790e-03_DP &
	   + T * (  5.45838830e-06_DP + T * ( -6.61810030e-09_DP &
	   + T * (  2.46570740e-12_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  7.78741500e-01_DP + T * (  8.73834000e-03_DP &
	   + T * ( -9.27803000e-06_DP + T * (  7.62427000e-09_DP &
	   + T * ( -2.44786200e-12_DP ) ) ) ) ) -9.82522900e+03_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    7.78741500e-01_DP + T * (  1.74766800e-02_DP &
	   + T * ( -2.78340900e-05_DP + T * (  3.04970800e-08_DP &
	   + T * ( -1.22393100e-11_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  2.69626120e+00_DP + T * (  2.46307115e-03_DP &
	   + T * (  2.76088313e-07_DP + T * ( -1.37595490e-10_DP &
	   + T * ( -7.92206520e-14_DP ) ) ) ) ) -1.49707930e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    2.69626120e+00_DP + T * (  4.92614230e-03_DP &
	   + T * (  8.28264940e-07_DP + T * ( -5.50381960e-10_DP &
	   + T * ( -3.96103260e-13_DP ) ) ) ) )
      H(SCH3O) =  2.67910679e+02_DP * ( &
	   T * (  2.10620400e+00_DP + T * (  3.60829750e-03_DP &
	   + T * (  1.77949067e-06_DP + T * ( -1.84440900e-09_DP &
	   + T * (  4.15122200e-13_DP ) ) ) ) ) +  9.78601100e+02_DP )
      CP(SCH3O) =  2.67910679e+02_DP * ( &
	    2.10620400e+00_DP + T * (  7.21659500e-03_DP &
	   + T * (  5.33847200e-06_DP + T * ( -7.37763600e-09_DP &
	   + T * (  2.07561100e-12_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  1.46253900e+00_DP + T * (  7.74733500e-03_DP &
	   + T * (  1.92683567e-06_DP + T * ( -3.14458000e-09_DP &
	   + T * (  9.17253400e-13_DP ) ) ) ) ) -1.12391800e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    1.46253900e+00_DP + T * (  1.54946700e-02_DP &
	   + T * (  5.78050700e-06_DP + T * ( -1.25783200e-08_DP &
	   + T * (  4.58626700e-12_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  2.76465589e+00_DP + T * (  5.97253050e-03_DP &
	   + T * ( -2.45923508e-06_DP + T * (  4.83935243e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.43520932e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    2.76465589e+00_DP + T * (  1.19450610e-02_DP &
	   + T * ( -7.37770524e-06_DP + T * (  1.93574097e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  4.30585800e+00_DP + T * ( -2.09168190e-03_DP &
	   + T * (  1.65690900e-05_DP + T * ( -1.49764685e-08_DP &
	   + T * (  4.60969560e-12_DP ) ) ) ) ) +  1.28417140e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    4.30585800e+00_DP + T * ( -4.18336380e-03_DP &
	   + T * (  4.97072700e-05_DP + T * ( -5.99058740e-08_DP &
	   + T * (  2.30484780e-11_DP ) ) ) ) )
      H(SCH2) =  5.92780550e+02_DP * ( &
	   T * (  3.76267867e+00_DP + T * (  4.84436072e-04_DP &
	   + T * (  9.31632803e-07_DP + T * ( -9.62727883e-10_DP &
	   + T * (  3.37483438e-13_DP ) ) ) ) ) +  4.60040401e+04_DP )
      CP(SCH2) =  5.92780550e+02_DP * ( &
	    3.76267867e+00_DP + T * (  9.68872143e-04_DP &
	   + T * (  2.79489841e-06_DP + T * ( -3.85091153e-09_DP &
	   + T * (  1.68741719e-12_DP ) ) ) ) )
      H(SCH2X) =  5.92780550e+02_DP * ( &
	   T * (  4.19860411e+00_DP + T * ( -1.18330710e-03_DP &
	   + T * (  2.74432073e-06_DP + T * ( -1.67203995e-09_DP &
	   + T * (  3.88629474e-13_DP ) ) ) ) ) +  5.04968163e+04_DP )
      CP(SCH2X) =  5.92780550e+02_DP * ( &
	    4.19860411e+00_DP + T * ( -2.36661419e-03_DP &
	   + T * (  8.23296220e-06_DP + T * ( -6.68815981e-09_DP &
	   + T * (  1.94314737e-12_DP ) ) ) ) )
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  3.95920148e+00_DP + T * ( -3.78526124e-03_DP &
	   + T * (  1.90330097e-05_DP + T * ( -1.72897188e-08_DP &
	   + T * (  5.39768746e-12_DP ) ) ) ) ) +  5.08977593e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    3.95920148e+00_DP + T * ( -7.57052247e-03_DP &
	   + T * (  5.70990292e-05_DP + T * ( -6.91588753e-08_DP &
	   + T * (  2.69884373e-11_DP ) ) ) ) )
      H(SCH3HCO) =  1.88739217e+02_DP * ( &
	   T * (  4.72945950e+00_DP + T * ( -1.59664290e-03_DP &
	   + T * (  1.58449737e-05_DP + T * ( -1.43646527e-08_DP &
	   + T * (  4.38622240e-12_DP ) ) ) ) ) -2.15728780e+04_DP )
      CP(SCH3HCO) =  1.88739217e+02_DP * ( &
	    4.72945950e+00_DP + T * ( -3.19328580e-03_DP &
	   + T * (  4.75349210e-05_DP + T * ( -5.74586110e-08_DP &
	   + T * (  2.19311120e-11_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  8.08681094e-01_DP + T * (  1.16807815e-02_DP &
	   + T * ( -1.18390605e-05_DP + T * (  7.00381092e-09_DP &
	   + T * ( -1.70014595e-12_DP ) ) ) ) ) +  2.64289807e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    8.08681094e-01_DP + T * (  2.33615629e-02_DP &
	   + T * ( -3.55171815e-05_DP + T * (  2.80152437e-08_DP &
	   + T * ( -8.50072974e-12_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.21246645e+00_DP + T * (  7.57395810e-04_DP &
	   + T * (  8.64031373e-06_DP + T * ( -8.94144617e-09_DP &
	   + T * (  2.94301746e-12_DP ) ) ) ) ) +  3.48598468e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.21246645e+00_DP + T * (  1.51479162e-03_DP &
	   + T * (  2.59209412e-05_DP + T * ( -3.57657847e-08_DP &
	   + T * (  1.47150873e-11_DP ) ) ) ) )
      H(SCH3OCH3) =  1.80479726e+02_DP * ( &
	   T * (  2.16630705e+00_DP + T * (  1.06687408e-02_DP &
	   + T * ( -2.13085193e-06_DP + T * ( -2.92811453e-11_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.36262464e+04_DP )
      CP(SCH3OCH3) =  1.80479726e+02_DP * ( &
	    2.16630705e+00_DP + T * (  2.13374815e-02_DP &
	   + T * ( -6.39255580e-06_DP + T * ( -1.17124581e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OCH2) =  1.84517088e+02_DP * ( &
	   T * (  2.91327415e+00_DP + T * (  1.01682329e-02_DP &
	   + T * ( -3.19904114e-06_DP + T * (  5.18696313e-10_DP &
	   + T * ( -3.42686724e-14_DP ) ) ) ) ) -1.18844240e+03_DP )
      CP(SCH3OCH2) =  1.84517088e+02_DP * ( &
	    2.91327415e+00_DP + T * (  2.03364659e-02_DP &
	   + T * ( -9.59712342e-06_DP + T * (  2.07478525e-09_DP &
	   + T * ( -1.71343362e-13_DP ) ) ) ) )
      H(SCH3OCH2O) =  1.36166721e+02_DP * ( &
	   T * (  3.25889339e+00_DP + T * (  1.11073180e-02_DP &
	   + T * ( -2.59518780e-06_DP + T * ( -6.03710395e-11_DP &
	   + T * (  9.03828992e-14_DP ) ) ) ) ) -1.92377212e+04_DP )
      CP(SCH3OCH2O) =  1.36166721e+02_DP * ( &
	    3.25889339e+00_DP + T * (  2.22146359e-02_DP &
	   + T * ( -7.78556340e-06_DP + T * ( -2.41484158e-10_DP &
	   + T * (  4.51914496e-13_DP ) ) ) ) )
      H(SCH3OCHO) =  1.38452341e+02_DP * ( &
	   T * (  3.08839783e+00_DP + T * (  1.01880024e-02_DP &
	   + T * ( -2.28259013e-06_DP + T * ( -1.82046551e-10_DP &
	   + T * (  1.12426043e-13_DP ) ) ) ) ) -4.41855167e+04_DP )
      CP(SCH3OCHO) =  1.38452341e+02_DP * ( &
	    3.08839783e+00_DP + T * (  2.03760048e-02_DP &
	   + T * ( -6.84777040e-06_DP + T * ( -7.28186203e-10_DP &
	   + T * (  5.62130216e-13_DP ) ) ) ) )
      H(SCH3OCO) =  1.40816002e+02_DP * ( &
	   T * (  3.94199159e+00_DP + T * (  1.21717442e-02_DP &
	   + T * ( -5.51985200e-06_DP + T * (  1.14634353e-09_DP &
	   + T * ( -6.63591416e-14_DP ) ) ) ) ) -2.14404829e+04_DP )
      CP(SCH3OCO) =  1.40816002e+02_DP * ( &
	    3.94199159e+00_DP + T * (  2.43434884e-02_DP &
	   + T * ( -1.65595560e-05_DP + T * (  4.58537411e-09_DP &
	   + T * ( -3.31795708e-13_DP ) ) ) ) )
      H(SRO2) =  1.07894368e+02_DP * ( &
	   T * (  2.21029612e+00_DP + T * (  1.84438727e-02_DP &
	   + T * ( -9.41871850e-06_DP + T * (  2.89326332e-09_DP &
	   + T * ( -3.94260940e-13_DP ) ) ) ) ) -1.94940940e+04_DP )
      CP(SRO2) =  1.07894368e+02_DP * ( &
	    2.21029612e+00_DP + T * (  3.68877454e-02_DP &
	   + T * ( -2.82561555e-05_DP + T * (  1.15730533e-08_DP &
	   + T * ( -1.97130470e-12_DP ) ) ) ) )
      H(SROH) =  1.33955339e+02_DP * ( &
	   T * (  3.15851876e+00_DP + T * (  1.22162875e-02_DP &
	   + T * ( -2.88994928e-06_DP + T * ( -1.48329832e-11_DP &
	   + T * (  8.72800006e-14_DP ) ) ) ) ) -4.54488899e+04_DP )
      CP(SROH) =  1.33955339e+02_DP * ( &
	    3.15851876e+00_DP + T * (  2.44325751e-02_DP &
	   + T * ( -8.66984784e-06_DP + T * ( -5.93319328e-11_DP &
	   + T * (  4.36400003e-13_DP ) ) ) ) )
      H(SQOOH) =  1.07894368e+02_DP * ( &
	   T * (  2.52895507e+00_DP + T * (  2.12064145e-02_DP &
	   + T * ( -1.24468795e-05_DP + T * (  4.16598332e-09_DP &
	   + T * ( -5.92886624e-13_DP ) ) ) ) ) -1.44293306e+04_DP )
      CP(SQOOH) =  1.07894368e+02_DP * ( &
	    2.52895507e+00_DP + T * (  4.24128290e-02_DP &
	   + T * ( -3.73406386e-05_DP + T * (  1.66639333e-08_DP &
	   + T * ( -2.96443312e-12_DP ) ) ) ) )
      H(SO2QOOH) =  7.62363836e+01_DP * ( &
	   T * (  1.99640551e+00_DP + T * (  2.91613116e-02_DP &
	   + T * ( -1.84419926e-05_DP + T * (  6.49526350e-09_DP &
	   + T * ( -9.54282010e-13_DP ) ) ) ) ) -3.27628742e+04_DP )
      CP(SO2QOOH) =  7.62363836e+01_DP * ( &
	    1.99640551e+00_DP + T * (  5.83226232e-02_DP &
	   + T * ( -5.53259778e-05_DP + T * (  2.59810540e-08_DP &
	   + T * ( -4.77141005e-12_DP ) ) ) ) )
      H(SHO2QHO) =  9.03222092e+01_DP * ( &
	   T * (  3.47935703e+00_DP + T * (  2.01476196e-02_DP &
	   + T * ( -1.10036432e-05_DP + T * (  3.35900293e-09_DP &
	   + T * ( -4.37203160e-13_DP ) ) ) ) ) -5.80629934e+04_DP )
      CP(SHO2QHO) =  9.03222092e+01_DP * ( &
	    3.47935703e+00_DP + T * (  4.02952392e-02_DP &
	   + T * ( -3.30109296e-05_DP + T * (  1.34360117e-08_DP &
	   + T * ( -2.18601580e-12_DP ) ) ) ) )
      H(SOCH2OCHO) =  1.10792868e+02_DP * ( &
	   T * (  5.19690837e+00_DP + T * (  7.94198615e-03_DP &
	   + T * (  1.17846849e-07_DP + T * ( -1.52614231e-09_DP &
	   + T * (  3.89323602e-13_DP ) ) ) ) ) -4.02242792e+04_DP )
      CP(SOCH2OCHO) =  1.10792868e+02_DP * ( &
	    5.19690837e+00_DP + T * (  1.58839723e-02_DP &
	   + T * (  3.53540547e-07_DP + T * ( -6.10456923e-09_DP &
	   + T * (  1.94661801e-12_DP ) ) ) ) )
      H(SHOCH2OCO) =  1.10792868e+02_DP * ( &
	   T * (  6.08180801e+00_DP + T * (  6.43841795e-03_DP &
	   + T * (  6.81398060e-07_DP + T * ( -1.52538730e-09_DP &
	   + T * (  3.59641118e-13_DP ) ) ) ) ) -4.39526183e+04_DP )
      CP(SHOCH2OCO) =  1.10792868e+02_DP * ( &
	    6.08180801e+00_DP + T * (  1.28768359e-02_DP &
	   + T * (  2.04419418e-06_DP + T * ( -6.10154921e-09_DP &
	   + T * (  1.79820559e-12_DP ) ) ) ) )
      H(SHOCH2O) =  1.76772973e+02_DP * ( &
	   T * (  4.11183145e+00_DP + T * (  3.76925348e-03_DP &
	   + T * (  1.25779123e-06_DP + T * ( -1.34686501e-09_DP &
	   + T * (  2.91231774e-13_DP ) ) ) ) ) -2.34414546e+04_DP )
      CP(SHOCH2O) =  1.76772973e+02_DP * ( &
	    4.11183145e+00_DP + T * (  7.53850697e-03_DP &
	   + T * (  3.77337370e-06_DP + T * ( -5.38746005e-09_DP &
	   + T * (  1.45615887e-12_DP ) ) ) ) )
      H(SHCOOH) =  1.80644418e+02_DP * ( &
	   T * (  1.43548185e+00_DP + T * (  8.16815080e-03_DP &
	   + T * ( -3.54191403e-06_DP + T * (  8.30332443e-10_DP &
	   + T * ( -8.04352206e-14_DP ) ) ) ) ) -4.64616504e+04_DP )
      CP(SHCOOH) =  1.80644418e+02_DP * ( &
	    1.43548185e+00_DP + T * (  1.63363016e-02_DP &
	   + T * ( -1.06257421e-05_DP + T * (  3.32132977e-09_DP &
	   + T * ( -4.02176103e-13_DP ) ) ) ) )
      END IF

      END
