!----------------------------------------------------------
! ======= grimech30F.f90 =======
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
      implicit none
      include 'grimech30F90.h'
      real(DP) :: CDOT(36), W(422), K(422), &
      C(36), M(32), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM1) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2.4 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.75 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 3.6 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 15.4 * C(SH2O) + 0.83 * C(SAR)
      M(MM1) = M(MM1) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM2) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM2) = M(MM2) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM3) = C(SN2) + C(SO) &
	    + 6 * C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 3.5 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.5 * C(SAR)
      M(MM3) = M(MM3) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM4) = C(SO) + C(SH) &
	    + C(SOH) + C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 0.75 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SCH4) &
	    + 1.5 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 1.5 * C(SC2H6) &
	    + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO)
      M(MM4) = M(MM4) + C(SC3H8) + C(SC3H7)
      M(MM5) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + C(SHO2) &
	    + C(SH2O2) + C(SCH) &
	    + C(SCO) + C(SCH2) &
	    + C(SHCO) + C(SCH2YXCH2) &
	    + C(SCH3) + C(SCH2O) &
	    + 2 * C(SCH4) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 0.63 * C(SAR) + C(SC) &
	    + C(SHCCOH) + C(SCH2CHO)
      M(MM5) = M(MM5) + C(SCH3CHO) + C(SC3H8) &
	    + C(SC3H7)
      M(MM6) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 0.73 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 3.65 * C(SH2O) + 0.38 * C(SAR)
      M(MM6) = M(MM6) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM7) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM7) = M(MM7) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM8) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 3 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM8) = M(MM8) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM9) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM9) = M(MM9) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM10) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM10) = M(MM10) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM11) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM11) = M(MM11) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM12) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM12) = M(MM12) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM13) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM13) = M(MM13) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM14) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM14) = M(MM14) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM15) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM15) = M(MM15) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM16) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM16) = M(MM16) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM17) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM17) = M(MM17) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM18) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM18) = M(MM18) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM19) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM19) = M(MM19) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM20) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM20) = M(MM20) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM21) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM21) = M(MM21) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM22) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM22) = M(MM22) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM23) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM23) = M(MM23) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM24) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + C(SAR)
      M(MM24) = M(MM24) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM25) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM25) = M(MM25) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM26) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + C(SAR) + C(SC)
      M(MM26) = M(MM26) + C(SHCCOH) + C(SCH2CHO) &
	    + C(SCH3CHO) + C(SC3H8) &
	    + C(SC3H7)
      M(MM27) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM27) = M(MM27) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM37) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM37) = M(MM37) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM38) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM38) = M(MM38) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM39) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM39) = M(MM39) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM40) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM40) = M(MM40) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)
      M(MM41) = C(SN2) + C(SO) &
	    + C(SO2) + C(SH) &
	    + C(SOH) + 2 * C(SH2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SCH2) + C(SHCO) &
	    + C(SCH2YXCH2) + C(SCH3) &
	    + C(SCH2O) + 2 * C(SCH4) &
	    + 2 * C(SCO2) + C(SCH2OH) &
	    + C(SCH3O) + C(SCH3OH) &
	    + C(SC2H) + C(SC2H2) &
	    + C(SHCCO) + C(SC2H3) &
	    + C(SCH2CO) + C(SC2H4) &
	    + C(SC2H5) + 3 * C(SC2H6) &
	    + 6 * C(SH2O) + 0.7 * C(SAR)
      M(MM41) = M(MM41) + C(SC) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SC3H8) + C(SC3H7)


      K(R1F) = 1.2000000000D+11 * exp(-1 * LT)
      K(R1B) = 8.4127616577D+15 &
	   * exp(-1.12163 * LT - 497822466 / RT)
      K(R2F) = 5.0000000000D+11 * exp(-1 * LT)
      K(R2B) = 6.9505130758D+13 &
	   * exp(-0.684828 * LT - 425203810.5 / RT)
      K(R3F) = 3.8700000000D+01 &
	   * exp(2.7 * LT - 26192000 / RT)
      K(R3B) = 2.3016447434D+01 &
	   * exp(2.66494 * LT - 18304322.73 / RT)
      K(R4F) = 2.0000000000D+10
      K(R4B) = 1.9300151207D+09 &
	   * exp(0.32701 * LT - 220535995.2 / RT)
      K(R5F) = 9.6300000000D+03 &
	   * exp(2 * LT - 16736000 / RT)
      K(R5B) = 1.9422606141D+01 &
	   * exp(2.5719 * LT - 75687466.21 / RT)
      K(R6F) = 5.7000000000D+10
      K(R6B) = 4.0039790898D+12 &
	   * exp(-0.0551873 * LT - 739371737.2 / RT)
      K(R7F) = 8.0000000000D+10
      K(R7B) = 7.8223107159D+12 &
	   * exp(-0.345666 * LT - 383372794.3 / RT)
      K(R8F) = 1.5000000000D+10
      K(R8B) = 4.7190290127D+10 &
	   * exp(0.174995 * LT - 787993698.2 / RT)
      K(R9F) = 1.5000000000D+10
      K(R9B) = 1.0693498414D+12 &
	   * exp(-0.407952 * LT - 420675782.2 / RT)
      K(R10F) = 5.0600000000D+10
      K(R10B) = 6.7208157338D+12 &
	   * exp(-0.275761 * LT - 288886347 / RT)
      K(R11F) = 1.0200000000D+06 &
	   * exp(1.5 * LT - 35982000 / RT)
      K(R11B) = 3.4821276777D+02 &
	   * exp(1.9765 * LT - 19030609.55 / RT)
      K0TROE = 6.0200000000D+08 * exp(-12552000 / RT)
      KINFTROE = 1.8000000000D+07 * exp(-9979000 / RT)
      FCTROE = 1 * EXP( -0 / TEMP )
      K(R12F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 4.6182917899D+17 &
	   * exp(-0.963346 * LT - 547255890 / RT)
      KINFTROE = 1.3808845883D+16 &
	   * exp(-0.963346 * LT - 544682890 / RT)
      K(R12B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(R13F) = 3.0000000000D+10
      K(R13B) = 7.8737463365D+08 &
	   * exp(0.547884 * LT - 359430238.7 / RT)
      K(R14F) = 3.0000000000D+10
      K(R14B) = 4.3452967068D+15 &
	   * exp(-0.730634 * LT - 468930318.2 / RT)
      K(R15F) = 3.9000000000D+10 * exp(-14811000 / RT)
      K(R15B) = 3.8627520933D+07 &
	   * exp(0.41439 * LT - 70409403.98 / RT)
      K(R16F) = 1.0000000000D+10
      K(R16B) = 1.3918415136D+08 &
	   * exp(0.696539 * LT - 302440224 / RT)
      K(R17F) = 1.0000000000D+10
      K(R17B) = 1.5736169844D+07 &
	   * exp(0.637464 * LT - 330967678.1 / RT)
      K(R18F) = 3.8800000000D+02 &
	   * exp(2.5 * LT - 12970000 / RT)
      K(R18B) = 1.7260888544D+00 &
	   * exp(2.67595 * LT - 32461855.27 / RT)
      K(R19F) = 1.3000000000D+02 &
	   * exp(2.5 * LT - 20920000 / RT)
      K(R19B) = 5.1152342447D+00 &
	   * exp(2.73503 * LT - 11884401.26 / RT)
      K(R20F) = 5.0000000000D+10
      K(R20B) = 2.3360600807D+07 &
	   * exp(1.01264 * LT - 325621371.1 / RT)
      K(R21F) = 1.3500000000D+04 &
	   * exp(2 * LT - 7950000 / RT)
      K(R21B) = 4.1635672871D+04 &
	   * exp(1.81798 * LT - 89764345.93 / RT)
      K(R22F) = 4.6000000000D+16 &
	   * exp(-1.41 * LT - 121127000 / RT)
      K(R22B) = 9.3604172115D+12 &
	   * exp(-0.848524 * LT + 10152055.94 / RT)
      K(R23F) = 6.9400000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(R23B) = 1.8060199928D-02 &
	   * exp(3.31672 * LT - 198861019.4 / RT)
      K(R24F) = 3.0000000000D+10
      K(R24B) = 1.2732201502D+14 &
	   * exp(-0.600835 * LT - 379074877 / RT)
      K(R25F) = 1.2500000000D+04 &
	   * exp(1.83 * LT - 920000 / RT)
      K(R25B) = 1.4453520437D-01 &
	   * exp(2.70863 * LT - 108375950.5 / RT)
      K(R26F) = 2.2400000000D+10
      K(R26B) = 1.9455952162D+07 &
	   * exp(0.900726 * LT - 325723724.7 / RT)
      K(R27F) = 8.9800000000D+04 &
	   * exp(1.92 * LT - 23807000 / RT)
      K(R27B) = 1.0047450020D+01 &
	   * exp(2.44412 * LT - 27335892.79 / RT)
      K(R28F) = 1.0000000000D+11
      K(R28B) = 1.5577208903D+03 &
	   * exp(1.38578 * LT - 426695895.8 / RT)
      K(R29F) = 1.0000000000D+10 * exp(-33472000 / RT)
      K(R29B) = 4.5588427293D+06 &
	   * exp(0.562107 * LT - 14614041.66 / RT)
      K(R30F) = 1.7500000000D+09 * exp(-5648000 / RT)
      K(R30B) = 3.7150267058D+06 &
	   * exp(0.782321 * LT - 205386794.6 / RT)
      K(R31F) = 2.5000000000D+09 * exp(-199995000 / RT)
      K(R31B) = 2.7356941842D+13 &
	   * exp(-0.84172 * LT - 236876424.1 / RT)
      K(R32F) = 1.0000000000D+11 * exp(-167360000 / RT)
      K(R32B) = 1.0263642436D+09 &
	   * exp(0.0873801 * LT - 2422408.767 / RT)
      K(R33F) = 2.8000000000D+12 * exp(-0.86 * LT)
      K(R33B) = 4.0334267652D+15 &
	   * exp(-0.871839 * LT - 204667815.3 / RT)
      K(R34F) = 2.0800000000D+13 * exp(-1.24 * LT)
      K(R34B) = 2.9962598827D+16 &
	   * exp(-1.25184 * LT - 204667815.3 / RT)
      K(R35F) = 1.1260000000D+13 * exp(-0.76 * LT)
      K(R35B) = 1.6220137634D+16 &
	   * exp(-0.771839 * LT - 204667815.3 / RT)
      K(R36F) = 2.6000000000D+13 * exp(-1.24 * LT)
      K(R36B) = 3.7453248534D+16 &
	   * exp(-1.25184 * LT - 204667815.3 / RT)
      K(R37F) = 7.0000000000D+11 * exp(-0.8 * LT)
      K(R37B) = 1.0083566913D+15 &
	   * exp(-0.811839 * LT - 204667815.3 / RT)
      K(R38F) = 2.6500000000D+13 &
	   * exp(-0.671 * LT - 71300000 / RT)
      K(R38B) = 5.2545483828D+10 &
	   * exp(-0.234202 * LT + 1318655.436 / RT)
      K(R39F) = 1.0000000000D+12 * exp(-1 * LT)
      K(R39B) = 2.3373273118D+14 &
	   * exp(-0.649765 * LT - 433091487.8 / RT)
      K(R40F) = 9.0000000000D+10 * exp(-0.6 * LT)
      K(R40B) = 2.1035945806D+13 &
	   * exp(-0.249765 * LT - 433091487.8 / RT)
      K(R41F) = 6.0000000000D+13 * exp(-1.25 * LT)
      K(R41B) = 1.4023963871D+16 &
	   * exp(-0.899765 * LT - 433091487.8 / RT)
      K(R42F) = 5.5000000000D+14 * exp(-2 * LT)
      K(R42B) = 1.2855300215D+17 &
	   * exp(-1.64977 * LT - 433091487.8 / RT)
      K(R43F) = 2.2000000000D+16 * exp(-2 * LT)
      K(R43B) = 5.7024497877D+19 &
	   * exp(-1.7546 * LT - 497934272.2 / RT)
      K(R44F) = 3.9700000000D+09 * exp(-2807000 / RT)
      K(R44B) = 1.4164546439D+07 &
	   * exp(0.694036 * LT - 223454801.5 / RT)
      K(R45F) = 4.4800000000D+10 * exp(-4469000 / RT)
      K(R45B) = 7.2691127188D+09 &
	   * exp(0.362073 * LT - 232892672.5 / RT)
      K(R46F) = 8.4000000000D+10 * exp(-2657000 / RT)
      K(R46B) = 1.6073095430D+07 &
	   * exp(0.763808 * LT - 150574339.8 / RT)
      K(R47F) = 1.2100000000D+04 &
	   * exp(2 * LT - 21757000 / RT)
      K(R47B) = 4.1033565943D+01 &
	   * exp(2.60696 * LT - 88596143.48 / RT)
      K(R48F) = 1.0000000000D+10 * exp(-15062000 / RT)
      K(R48B) = 7.1960369081D+04 &
	   * exp(1.26593 * LT - 294661267.7 / RT)
      K(R49F) = 1.6500000000D+11
      K(R49B) = 8.4149682406D+10 &
	   * exp(0.23442 * LT - 97396261.59 / RT)
      K0TROE = 1.0400000000D+20 &
	   * exp(-2.76 * LT - 6694000 / RT)
      KINFTROE = 6.0000000000D+11
      FCTROE = 0.438 * EXP( -TEMP / 91 ) &
	   + 0.562 * EXP( -TEMP / 5836 ) &
	   + 1 * EXP( -8552 / TEMP )
      K(R50F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K0TROE = 1.0745374952D+25 &
	   * exp(-2.92912 * LT - 470785853.8 / RT)
      KINFTROE = 6.1992547800D+16 &
	   * exp(-0.169123 * LT - 464091853.8 / RT)
      K(R50B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K(R51F) = 3.0000000000D+10
      K(R51B) = 1.3435867056D+09 &
	   * exp(0.230182 * LT - 48621961.01 / RT)
      K0TROE = 2.6200000000D+27 &
	   * exp(-4.76 * LT - 10209000 / RT)
      KINFTROE = 1.3900000000D+13 &
	   * exp(-0.534 * LT - 2243000 / RT)
      FCTROE = 0.217 * EXP( -TEMP / 74 ) &
	   + 0.783 * EXP( -TEMP / 2941 ) &
	   + 1 * EXP( -6964 / TEMP )
      K(R52F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM8) )
      K0TROE = 1.0668506651D+33 &
	   * exp(-4.92133 * LT - 452364201 / RT)
      KINFTROE = 5.6600092537D+18 &
	   * exp(-0.695325 * LT - 444398201 / RT)
      K(R52B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM8) )
      K(R53F) = 6.6000000000D+05 &
	   * exp(1.62 * LT - 45355000 / RT)
      K(R53B) = 3.7884462370D+02 &
	   * exp(2.13156 * LT - 36291286.82 / RT)
      K0TROE = 2.4700000000D+18 &
	   * exp(-2.57 * LT - 1778000 / RT)
      KINFTROE = 1.0900000000D+09 &
	   * exp(0.48 * LT + 1088000 / RT)
      FCTROE = 0.2176 * EXP( -TEMP / 271 ) &
	   + 0.7824 * EXP( -TEMP / 2755 ) &
	   + 1 * EXP( -6570 / TEMP )
      K(R54F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K0TROE = 3.4666626717D+23 &
	   * exp(-2.66922 * LT - 371383406.5 / RT)
      KINFTROE = 1.5298227985D+14 &
	   * exp(0.380781 * LT - 368517406.5 / RT)
      K(R54B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K(R55F) = 7.3400000000D+10
      K(R55B) = 3.2391338749D+09 &
	   * exp(0.582947 * LT - 367317915.9 / RT)
      K0TROE = 1.2700000000D+26 &
	   * exp(-4.82 * LT - 27322000 / RT)
      KINFTROE = 5.4000000000D+08 &
	   * exp(0.454 * LT - 15062000 / RT)
      FCTROE = 0.2813 * EXP( -TEMP / 103 ) &
	   + 0.7187 * EXP( -TEMP / 1291 ) &
	   + 1 * EXP( -4160 / TEMP )
      K(R56F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM10) )
      K0TROE = 1.2684133243D+30 &
	   * exp(-5.20137 * LT - 150085586.5 / RT)
      KINFTROE = 5.3932535051D+12 &
	   * exp(0.0726333 * LT - 137825586.5 / RT)
      K(R56B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM10) )
      K0TROE = 2.2000000000D+24 &
	   * exp(-4.8 * LT - 23263000 / RT)
      KINFTROE = 5.4000000000D+08 &
	   * exp(0.454 * LT - 10878000 / RT)
      FCTROE = 0.242 * EXP( -TEMP / 94 ) &
	   + 0.758 * EXP( -TEMP / 1555 ) &
	   + 1 * EXP( -4200 / TEMP )
      K(R57F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM11) )
      K0TROE = 1.9434371792D+29 &
	   * exp(-5.12229 * LT - 117499132.5 / RT)
      KINFTROE = 4.7702548944D+13 &
	   * exp(0.131708 * LT - 105114132.5 / RT)
      K(R57B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM11) )
      K(R58F) = 5.7400000000D+04 &
	   * exp(1.9 * LT - 11473000 / RT)
      K(R58B) = 9.5590953891D+01 &
	   * exp(2.34945 * LT - 74959081.25 / RT)
      K0TROE = 4.3600000000D+25 &
	   * exp(-4.65 * LT - 21255000 / RT)
      KINFTROE = 1.0550000000D+09 &
	   * exp(0.5 * LT - 360000 / RT)
      FCTROE = 0.4 * EXP( -TEMP / 100 ) &
	   + 0.6 * EXP( -TEMP / 90000 ) &
	   + 1 * EXP( -10000 / TEMP )
      K(R59F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM12) )
      K0TROE = 1.3623915050D+30 &
	   * exp(-4.51078 * LT - 426966955.3 / RT)
      KINFTROE = 3.2966124719D+13 &
	   * exp(0.639219 * LT - 406071955.3 / RT)
      K(R59B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM12) )
      K(R60F) = 2.0000000000D+10
      K(R60B) = 4.6805022133D+08 &
	   * exp(0.731602 * LT - 310327901.3 / RT)
      K(R61F) = 1.6500000000D+08 &
	   * exp(0.65 * LT + 1188000 / RT)
      K(R61B) = 1.7290289241D+04 &
	   * exp(1.6223 * LT - 12365877.06 / RT)
      K(R62F) = 3.2800000000D+10 &
	   * exp(-0.09 * LT - 2552000 / RT)
      K(R62B) = 1.1826581729D+05 &
	   * exp(1.35911 * LT - 12645307.5 / RT)
      K0TROE = 4.6600000000D+35 &
	   * exp(-7.44 * LT - 58911000 / RT)
      KINFTROE = 2.4300000000D+09 &
	   * exp(0.515 * LT - 209000 / RT)
      FCTROE = 0.3 * EXP( -TEMP / 100 ) &
	   + 0.7 * EXP( -TEMP / 90000 ) &
	   + 1 * EXP( -10000 / TEMP )
      K(R63F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM13) )
      K0TROE = 1.6463061592D+39 &
	   * exp(-7.35986 * LT - 493150409.3 / RT)
      KINFTROE = 8.5848153795D+12 &
	   * exp(0.595144 * LT - 434448409.3 / RT)
      K(R63B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM13) )
      K(R64F) = 4.1500000000D+04 &
	   * exp(1.63 * LT - 8050000 / RT)
      K(R64B) = 4.6919928896D+03 &
	   * exp(1.57093 * LT - 36577454.01 / RT)
      K(R65F) = 2.0000000000D+10
      K(R65B) = 5.2917790613D+07 &
	   * exp(0.672527 * LT - 338855355.3 / RT)
      K(R66F) = 1.5000000000D+09 &
	   * exp(0.5 * LT + 460000 / RT)
      K(R66B) = 1.7771284596D+04 &
	   * exp(1.41323 * LT - 41621331.07 / RT)
      K(R67F) = 2.6200000000D+11 &
	   * exp(-0.23 * LT - 4477000 / RT)
      K(R67B) = 1.0680606960D+05 &
	   * exp(1.16003 * LT - 43097761.51 / RT)
      K(R68F) = 1.7000000000D+04 &
	   * exp(2.1 * LT - 20376000 / RT)
      K(R68B) = 1.2716073149D+02 &
	   * exp(2.31102 * LT - 47755532.54 / RT)
      K(R69F) = 4.2000000000D+03 &
	   * exp(2.1 * LT - 20376000 / RT)
      K(R69B) = 2.7787158476D+02 &
	   * exp(2.37009 * LT - 19228078.53 / RT)
      K0TROE = 3.7500000000D+27 &
	   * exp(-4.8 * LT - 7950000 / RT)
      KINFTROE = 1.0000000000D+14 * exp(-1 * LT)
      FCTROE = 0.3536 * EXP( -TEMP / 132 ) &
	   + 0.6464 * EXP( -TEMP / 1315 ) &
	   + 1 * EXP( -5566 / TEMP )
      K(R70F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM14) )
      K0TROE = 2.5617736443D+33 &
	   * exp(-5.0463 * LT - 564432866.5 / RT)
      KINFTROE = 6.8313963848D+19 &
	   * exp(-1.2463 * LT - 556482866.5 / RT)
      K(R70B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM14) )
      K0TROE = 3.8000000000D+34 &
	   * exp(-7.27 * LT - 30208000 / RT)
      KINFTROE = 5.6000000000D+09 * exp(-10042000 / RT)
      FCTROE = 0.2493 * EXP( -TEMP / 98.5 ) &
	   + 0.7507 * EXP( -TEMP / 1302 ) &
	   + 1 * EXP( -4167 / TEMP )
      K(R71F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM15) )
      K0TROE = 8.4202583646D+36 &
	   * exp(-7.09812 * LT - 177009237.8 / RT)
      KINFTROE = 1.2408801800D+12 &
	   * exp(0.171883 * LT - 156843237.8 / RT)
      K(R71B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM15) )
      K0TROE = 1.4000000000D+24 &
	   * exp(-3.86 * LT - 13891000 / RT)
      KINFTROE = 6.0800000000D+09 &
	   * exp(0.27 * LT - 1172000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(R72F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM16) )
      K0TROE = 7.7814512533D+29 &
	   * exp(-3.99563 * LT - 480410256.8 / RT)
      KINFTROE = 3.3793731157D+15 &
	   * exp(0.134367 * LT - 467691256.8 / RT)
      K(R72B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM16) )
      K(R73F) = 3.0000000000D+10
      K(R73B) = 3.1644553172D+10 &
	   * exp(0.178352 * LT - 286290250 / RT)
      K0TROE = 6.0000000000D+35 &
	   * exp(-7.62 * LT - 29162000 / RT)
      KINFTROE = 5.4000000000D+08 &
	   * exp(0.454 * LT - 7615000 / RT)
      FCTROE = 0.0247 * EXP( -TEMP / 210 ) &
	   + 0.9753 * EXP( -TEMP / 984 ) &
	   + 1 * EXP( -4374 / TEMP )
      K(R74F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM17) )
      K0TROE = 1.1210501411D+39 &
	   * exp(-7.74131 * LT - 180499632.3 / RT)
      KINFTROE = 1.0089451270D+12 &
	   * exp(0.332687 * LT - 158952632.3 / RT)
      K(R74B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM17) )
      K(R75F) = 1.3250000000D+03 &
	   * exp(2.53 * LT - 51212000 / RT)
      K(R75B) = 5.5718940109D-01 &
	   * exp(3.01587 * LT - 17784230.95 / RT)
      K0TROE = 1.9900000000D+35 &
	   * exp(-7.08 * LT - 27970000 / RT)
      KINFTROE = 5.2100000000D+14 &
	   * exp(-0.99 * LT - 6611000 / RT)
      FCTROE = 0.1578 * EXP( -TEMP / 125 ) &
	   + 0.8422 * EXP( -TEMP / 2219 ) &
	   + 1 * EXP( -6882 / TEMP )
      K(R76F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM18) )
      K0TROE = 2.4724095869D+41 &
	   * exp(-7.28894 * LT - 449644917.7 / RT)
      KINFTROE = 6.4729919335D+20 &
	   * exp(-1.19894 * LT - 428285917.7 / RT)
      K(R76B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM18) )
      K(R77F) = 2.0000000000D+09
      K(R77B) = 2.5019333849D+08 &
	   * exp(0.471548 * LT - 281753855.5 / RT)
      K(R78F) = 1.1500000000D+05 &
	   * exp(1.9 * LT - 31506000 / RT)
      K(R78B) = 2.1634657872D+01 &
	   * exp(2.45918 * LT - 42922570.06 / RT)
      K(R79F) = 1.0000000000D+11
      K(R79B) = 1.1573048940D+05 &
	   * exp(1.56102 * LT - 71793685.47 / RT)
      K(R80F) = 5.0000000000D+10 * exp(-33472000 / RT)
      K(R80B) = 3.8326334707D+07 &
	   * exp(0.59717 * LT - 22501718.93 / RT)
      K(R81F) = 1.1300000000D+10 * exp(-14343000 / RT)
      K(R81B) = 3.2307674042D+03 &
	   * exp(1.57654 * LT - 143469758.4 / RT)
      K(R82F) = 1.0000000000D+10
      K(R82B) = 3.7197289740D+08 &
	   * exp(0.411006 * LT - 125789549.2 / RT)
      K0TROE = 5.0700000000D+21 &
	   * exp(-3.42 * LT - 352920000 / RT)
      KINFTROE = 4.3000000000D+04 &
	   * exp(1.5 * LT - 333046000 / RT)
      FCTROE = 0.068 * EXP( -TEMP / 197 ) &
	   + 0.932 * EXP( -TEMP / 1540 ) &
	   + 1 * EXP( -10300 / TEMP )
      K(R83F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM19) )
      K0TROE = 1.6124629783D+28 &
	   * exp(-4.10217 * LT - 355207490.6 / RT)
      KINFTROE = 1.3675721512D+11 &
	   * exp(0.817835 * LT - 335333490.6 / RT)
      K(R83B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM19) )
      K(R84F) = 2.1600000000D+05 &
	   * exp(1.51 * LT - 14351000 / RT)
      K(R84B) = 2.3953722075D+06 &
	   * exp(1.40516 * LT - 79193784.45 / RT)
      K0TROE = 2.3000000000D+12 &
	   * exp(-0.9 * LT + 7113000 / RT)
      KINFTROE = 7.4000000000D+10 * exp(-0.37 * LT)
      FCTROE = 0.2654 * EXP( -TEMP / 94 ) &
	   + 0.7346 * EXP( -TEMP / 1756 ) &
	   + 1 * EXP( -5182 / TEMP )
      K(R85F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM20) )
      K0TROE = 8.2846323983D+20 &
	   * exp(-1.92053 * LT - 211222004.5 / RT)
      KINFTROE = 2.6654904238D+19 &
	   * exp(-1.39053 * LT - 218335004.5 / RT)
      K(R85B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM20) )
      K(R86F) = 3.5700000000D+01 &
	   * exp(2.4 * LT + 8828000 / RT)
      K(R86B) = 6.6567177782D+02 &
	   * exp(2.33023 * LT - 63902461.72 / RT)
      K(R87F) = 1.4500000000D+10 * exp(2092000 / RT)
      K(R87B) = 2.6090995309D+10 &
	   * exp(0.257238 * LT - 291174456.9 / RT)
      K(R88F) = 2.0000000000D+09 * exp(-1787000 / RT)
      K(R88B) = 7.5214771530D+07 &
	   * exp(0.502125 * LT - 133468927.9 / RT)
      K(R89F) = 1.7000000000D+15 * exp(-123051000 / RT)
      K(R89B) = 6.3932555801D+13 &
	   * exp(0.502125 * LT - 254732927.9 / RT)
      K(R90F) = 5.0000000000D+10
      K(R90B) = 1.1579534098D+13 &
	   * exp(-0.254545 * LT - 649863152.8 / RT)
      K(R91F) = 3.0000000000D+10
      K(R91B) = 8.0293064785D+13 &
	   * exp(-0.603071 * LT - 379941498.5 / RT)
      K(R92F) = 2.0000000000D+10
      K(R92B) = 1.9744350049D+15 &
	   * exp(-0.760056 * LT - 327774390.3 / RT)
      K(R93F) = 1.1300000000D+04 &
	   * exp(2 * LT - 12552000 / RT)
      K(R93B) = 7.6976589861D+03 &
	   * exp(2.18763 * LT - 88713757.49 / RT)
      K(R94F) = 3.0000000000D+10
      K(R94B) = 2.1593228251D+15 &
	   * exp(-0.822343 * LT - 365077378.3 / RT)
      K0TROE = 4.0000000000D+30 &
	   * exp(-5.92 * LT - 13138000 / RT)
      KINFTROE = 2.7900000000D+15 &
	   * exp(-1.43 * LT - 5565000 / RT)
      FCTROE = 0.588 * EXP( -TEMP / 195 ) &
	   + 0.412 * EXP( -TEMP / 5900 ) &
	   + 1 * EXP( -6394 / TEMP )
      K(R95F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM21) )
      K0TROE = 1.1927711189D+39 &
	   * exp(-6.75308 * LT - 405296078.2 / RT)
      KINFTROE = 8.3195785542D+23 &
	   * exp(-2.26308 * LT - 397723078.2 / RT)
      K(R95B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM21) )
      K(R96F) = 5.6000000000D+04 &
	   * exp(1.6 * LT - 22677000 / RT)
      K(R96B) = 1.4048779004D+03 &
	   * exp(2.01452 * LT - 56519418.41 / RT)
      K(R97F) = 6.4400000000D+14 &
	   * exp(-1.34 * LT - 5929000 / RT)
      K(R97B) = 2.2159142320D+13 &
	   * exp(-0.863191 * LT - 2468430.436 / RT)
      K(R98F) = 1.0000000000D+05 &
	   * exp(1.6 * LT - 13054000 / RT)
      K(R98B) = 6.3655575377D+02 &
	   * exp(2.00672 * LT - 68833071.27 / RT)
      K(R99F) = 4.7600000000D+04 &
	   * exp(1.228 * LT - 293000 / RT)
      K(R99B) = 2.6269086456D+11 &
	   * exp(-0.0505181 * LT - 109793079.5 / RT)
      K(R100F) = 5.0000000000D+10
      K(R100B) = 2.4469331101D+10 &
	   * exp(0.478112 * LT - 432160700.4 / RT)
      K(R101F) = 3.4300000000D+06 &
	   * exp(1.18 * LT + 1870000 / RT)
      K(R101B) = 6.3345866067D+04 &
	   * exp(1.52462 * LT - 126458865.7 / RT)
      K(R102F) = 5.0000000000D+09
      K(R102B) = 1.2976325138D+09 &
	   * exp(0.626766 * LT - 375170685.8 / RT)
      K(R103F) = 5.0000000000D+09
      K(R103B) = 1.4671042236D+08 &
	   * exp(0.567692 * LT - 403698139.8 / RT)
      K(R104F) = 1.4400000000D+03 &
	   * exp(2 * LT + 3515000 / RT)
      K(R104B) = 1.1944991455D+02 &
	   * exp(2.10618 * LT - 88707316.99 / RT)
      K(R105F) = 6.3000000000D+03 &
	   * exp(2 * LT - 6276000 / RT)
      K(R105B) = 4.6222629959D+03 &
	   * exp(2.16526 * LT - 69970862.98 / RT)
      K(R106F) = 2.0000000000D+10
      K(R106B) = 3.0312687340D+14 &
	   * exp(-0.743493 * LT - 213093401.9 / RT)
      K(R107F) = 2.1800000000D-07 &
	   * exp(4.5 * LT + 4184000 / RT)
      K(R107B) = 1.4748019480D-03 &
	   * exp(3.75588 * LT - 96488304.27 / RT)
      K(R108F) = 5.0400000000D+02 &
	   * exp(2.3 * LT - 56484000 / RT)
      K(R108B) = 9.1663502510D+07 &
	   * exp(1.14487 * LT - 31366755.02 / RT)
      K(R109F) = 3.3700000000D+04 &
	   * exp(2 * LT - 58576000 / RT)
      K(R109B) = 1.2786721446D+02 &
	   * exp(2.4917 * LT - 27405.78234 / RT)
      K(R110F) = 4.8300000000D-07 &
	   * exp(4 * LT + 8368000 / RT)
      K(R110B) = 9.3422522348D-10 &
	   * exp(4.83242 * LT - 221431062.7 / RT)
      K(R111F) = 5.0000000000D+09
      K(R111B) = 5.8488027152D+10 &
	   * exp(0.0735165 * LT - 351133034.5 / RT)
      K(R112F) = 3.6000000000D+03 &
	   * exp(2 * LT - 10460000 / RT)
      K(R112B) = 1.6788377430D+01 &
	   * exp(2.38103 * LT - 41875015.4 / RT)
      K(R113F) = 3.5400000000D+03 &
	   * exp(2.12 * LT - 3640000 / RT)
      K(R113B) = 7.3854116736D+00 &
	   * exp(2.57434 * LT - 79899354.5 / RT)
      K(R114F) = 7.5000000000D+09 * exp(-8368000 / RT)
      K(R114B) = 6.3754053454D+07 &
	   * exp(0.492334 * LT - 62240503.38 / RT)
      K(R115F) = 1.3000000000D+08 * exp(6820000 / RT)
      K(R115B) = 6.2200353343D+09 &
	   * exp(-0.244887 * LT - 154764529 / RT)
      K(R116F) = 4.2000000000D+11 * exp(-50208000 / RT)
      K(R116B) = 2.0095498773D+13 &
	   * exp(-0.244887 * LT - 211792529 / RT)
      K(R117F) = 2.0000000000D+10
      K(R117B) = 3.7780097921D+11 &
	   * exp(0.00375206 * LT - 475691730.1 / RT)
      K(R118F) = 1.0000000000D+09
      K(R118B) = 2.8267421607D+11 &
	   * exp(-0.149487 * LT - 237487385.7 / RT)
      K(R119F) = 3.7800000000D+10
      K(R119B) = 6.1049832144D+11 &
	   * exp(-0.149417 * LT - 105836008.7 / RT)
      K(R120F) = 1.5000000000D+11 * exp(-98742000 / RT)
      K(R120B) = 1.5839793423D+14 &
	   * exp(-0.51471 * LT - 356159419.3 / RT)
      K(R121F) = 5.6000000000D+03 &
	   * exp(2 * LT - 50208000 / RT)
      K(R121B) = 2.7500401862D+03 &
	   * exp(1.84249 * LT - 46854937.77 / RT)
      K(R122F) = 5.8000000000D+10 * exp(-2410000 / RT)
      K(R122B) = 2.6634134985D+10 &
	   * exp(0.182253 * LT - 579654497.4 / RT)
      K(R123F) = 5.0000000000D+10
      K(R123B) = 9.0545188562D+14 &
	   * exp(-1.00978 * LT - 327673077.5 / RT)
      K(R124F) = 5.0000000000D+10
      K(R124B) = 5.9866880368D+15 &
	   * exp(-1.08697 * LT - 420064090.2 / RT)
      K(R125F) = 6.7100000000D+10
      K(R125B) = 3.5609741582D+11 &
	   * exp(-0.166273 * LT - 307322843.1 / RT)
      K(R126F) = 1.0800000000D+11 * exp(-13012000 / RT)
      K(R126B) = 1.7581777781D+12 &
	   * exp(-0.292468 * LT - 1693026.963 / RT)
      K(R127F) = 5.7100000000D+09 * exp(3159000 / RT)
      K(R127B) = 8.2750138742D+14 &
	   * exp(-0.947689 * LT - 248453632.8 / RT)
      K(R128F) = 4.0000000000D+10
      K(R128B) = 1.0797267305D+18 &
	   * exp(-1.3719 * LT - 548460717.8 / RT)
      K(R129F) = 3.0000000000D+10
      K(R129B) = 1.7367155877D+15 &
	   * exp(-1.0309 * LT - 231170101.7 / RT)
      K(R130F) = 6.0000000000D+10
      K(R130B) = 4.7412070602D+15 &
	   * exp(-1.0052 * LT - 255534157.6 / RT)
      K0TROE = 2.6900000000D+22 &
	   * exp(-3.74 * LT - 8100000 / RT)
      KINFTROE = 5.0000000000D+10
      FCTROE = 0.4243 * EXP( -TEMP / 237 ) &
	   + 0.5757 * EXP( -TEMP / 1652 ) &
	   + 1 * EXP( -5069 / TEMP )
      K(R131F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM22) )
      K0TROE = 1.2130524545D+32 &
	   * exp(-5.18097 * LT - 320775841.3 / RT)
      KINFTROE = 2.2547443392D+20 &
	   * exp(-1.44097 * LT - 312675841.3 / RT)
      K(R131B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM22) )
      K(R132F) = 1.9000000000D+11 * exp(-66074000 / RT)
      K(R132B) = 9.2145125176D+07 &
	   * exp(0.675447 * LT - 336515419 / RT)
      K(R133F) = 9.4600000000D+10 * exp(2155000 / RT)
      K(R133B) = 1.7498844055D+17 &
	   * exp(-1.35597 * LT - 319203631.8 / RT)
      K(R134F) = 5.0000000000D+10
      K(R134B) = 1.1388201171D+12 &
	   * exp(0.12683 * LT - 657557391.2 / RT)
      K(R135) = 5.0000000000D+09 * exp(-6276000 / RT)
      K(R136F) = 5.0000000000D+02 &
	   * exp(2 * LT - 30250000 / RT)
      K(R136B) = 2.2102362917D+05 &
	   * exp(1.48064 * LT - 61250366.04 / RT)
      K(R137F) = 1.6000000000D+12 * exp(-49974000 / RT)
      K(R137B) = 2.6529851155D+18 &
	   * exp(-1.07943 * LT - 609753690.8 / RT)
      K(R138F) = 4.0000000000D+10
      K(R138B) = 3.3825283043D+17 &
	   * exp(-1.2243 * LT - 275916843.8 / RT)
      K(R139F) = 2.4600000000D+03 &
	   * exp(2 * LT - 34602000 / RT)
      K(R139B) = 6.2419602881D+02 &
	   * exp(1.9922 * LT - 56538652.86 / RT)
      K0TROE = 2.6900000000D+27 &
	   * exp(-5.11 * LT - 29685000 / RT)
      KINFTROE = 8.1000000000D+08 &
	   * exp(0.5 * LT - 18870000 / RT)
      FCTROE = 0.4093 * EXP( -TEMP / 275 ) &
	   + 0.5907 * EXP( -TEMP / 1226 ) &
	   + 1 * EXP( -5185 / TEMP )
      K(R140F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM23) )
      K0TROE = 9.7210519966D+38 &
	   * exp(-6.85567 * LT - 364650095.4 / RT)
      KINFTROE = 2.9271569209D+20 &
	   * exp(-1.24567 * LT - 353835095.4 / RT)
      K(R140B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM23) )
      K(R141F) = 3.0000000000D+10
      K(R141B) = 3.9791438217D+10 &
	   * exp(0.240947 * LT - 382586114.2 / RT)
      K(R142F) = 1.5000000000D+10 * exp(-2510000 / RT)
      K(R142B) = 1.0936408744D+10 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(R143F) = 9.0000000000D+09 * exp(-2510000 / RT)
      K(R143B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(R144F) = 2.8000000000D+10
      K(R144B) = 7.4729040691D+05 &
	   * exp(0.261558 * LT - 282283554.9 / RT)
      K(R145F) = 1.2000000000D+10
      K(R145B) = 8.3014013421D+08 &
	   * exp(0.506957 * LT - 780217827.2 / RT)
      K(R146F) = 7.0000000000D+10
      K(R146B) = 2.2560577672D+13 &
	   * exp(-0.581645 * LT - 68303354.01 / RT)
      K0TROE = 1.8800000000D+32 &
	   * exp(-6.36 * LT - 21087000 / RT)
      KINFTROE = 4.8200000000D+14 &
	   * exp(-1.16 * LT - 4791000 / RT)
      FCTROE = 0.3973 * EXP( -TEMP / 208 ) &
	   + 0.6027 * EXP( -TEMP / 3922 ) &
	   + 1 * EXP( -10180 / TEMP )
      K(R147F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM24) )
      K0TROE = 1.6292506138D+42 &
	   * exp(-7.66989 * LT - 416705647.8 / RT)
      KINFTROE = 4.1771212546D+24 &
	   * exp(-2.46989 * LT - 400409647.8 / RT)
      K(R147B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM24) )
      K(R148F) = 3.0000000000D+10
      K(R148B) = 2.1872817488D+10 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R149F) = 1.2000000000D+10 * exp(2385000 / RT)
      K(R149B) = 7.3985424247D+16 &
	   * exp(-1.28658 * LT - 310834831.8 / RT)
      K(R150F) = 1.6000000000D+10 * exp(2385000 / RT)
      K(R150B) = 2.9599839164D+09 &
	   * exp(-0.0700846 * LT - 56854640.83 / RT)
      K(R151F) = 9.0000000000D+09
      K(R151B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R152F) = 7.0000000000D+09
      K(R152B) = 5.1036574138D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R153F) = 1.4000000000D+10
      K(R153B) = 1.8259393135D+08 &
	   * exp(0.456176 * LT - 255577298.7 / RT)
      K(R154F) = 4.0000000000D+10 * exp(2301000 / RT)
      K(R154B) = 2.4252938103D+09 &
	   * exp(-0.0224647 * LT - 77418924.07 / RT)
      K(R155F) = 3.5600000000D+10 * exp(-127528000 / RT)
      K(R155B) = 5.9581571870D+12 &
	   * exp(-0.476427 * LT - 12828013.49 / RT)
      K(R156F) = 2.3100000000D+09 * exp(-84998000 / RT)
      K(R156B) = 6.0837726033D+08 &
	   * exp(0.161037 * LT - 301265691.5 / RT)
      K(R157F) = 2.4500000000D+01 &
	   * exp(2.47 * LT - 21673000 / RT)
      K(R157B) = 1.4474473694D+02 &
	   * exp(2.5654 * LT - 97575856.67 / RT)
      K0TROE = 3.4000000000D+35 &
	   * exp(-7.03 * LT - 11556000 / RT)
      KINFTROE = 6.7700000000D+13 &
	   * exp(-1.18 * LT - 2736000 / RT)
      FCTROE = 0.381 * EXP( -TEMP / 73.2 ) &
	   + 0.619 * EXP( -TEMP / 1180 ) &
	   + 1 * EXP( -9999 / TEMP )
      K(R158F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM25) )
      K0TROE = 6.4597135449D+46 &
	   * exp(-8.41543 * LT - 396393540 / RT)
      KINFTROE = 1.2862429617D+25 &
	   * exp(-2.56543 * LT - 387573540 / RT)
      K(R158B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM25) )
      K(R159F) = 6.8400000000D+09 &
	   * exp(0.1 * LT - 44350000 / RT)
      K(R159B) = 1.0459793288D+15 &
	   * exp(-1.07649 * LT - 7512622.297 / RT)
      K(R160F) = 2.6480000000D+10
      K(R160B) = 2.0357930432D+12 &
	   * exp(0.0713867 * LT - 376381629.1 / RT)
      K(R161F) = 3.3200000000D+00 &
	   * exp(2.81 * LT - 24518000 / RT)
      K(R161B) = 9.6322067040D+00 &
	   * exp(2.74789 * LT - 97067794.44 / RT)
      K(R162F) = 3.0000000000D+04 &
	   * exp(1.5 * LT - 41589000 / RT)
      K(R162B) = 3.9093824409D+05 &
	   * exp(1.19946 * LT - 78032245.72 / RT)
      K(R163F) = 1.0000000000D+04 &
	   * exp(1.5 * LT - 41589000 / RT)
      K(R163B) = 1.1525974507D+06 &
	   * exp(1.25853 * LT - 49504791.71 / RT)
      K(R164F) = 2.2700000000D+02 &
	   * exp(2 * LT - 38493000 / RT)
      K(R164B) = 1.6630129560D+02 &
	   * exp(1.97431 * LT - 14128944.14 / RT)
      K(R165F) = 6.1400000000D+03 &
	   * exp(1.74 * LT - 43723000 / RT)
      K(R165B) = 2.0123493121D+03 &
	   * exp(1.78762 * LT - 64203283.24 / RT)
      K(R166F) = 1.5000000000D+15 &
	   * exp(-1 * LT - 71128000 / RT)
      K(R166B) = 2.8320737803D+11 &
	   * exp(-0.767288 * LT - 5354428.13 / RT)
      K(R167F) = 1.8700000000D+14 &
	   * exp(-1 * LT - 71128000 / RT)
      K(R167B) = 3.5306519794D+10 &
	   * exp(-0.767288 * LT - 5354428.13 / RT)
      K(R168F) = 1.3450000000D+10 * exp(-1674000 / RT)
      K(R168B) = 3.6580676525D+09 &
	   * exp(0.220873 * LT - 140568243.4 / RT)
      K(R169F) = 1.8000000000D+10 * exp(-3766000 / RT)
      K(R169B) = 2.5961607220D+09 &
	   * exp(0.369528 * LT - 85670228.82 / RT)
      K(R170F) = 4.2800000000D-16 &
	   * exp(7.6 * LT + 14770000 / RT)
      K(R170B) = 6.9793035518D-18 &
	   * exp(7.91045 * LT - 95661682.83 / RT)
      K(R171F) = 1.0000000000D+10 * exp(3159000 / RT)
      K(R171B) = 2.4794782651D+07 &
	   * exp(0.846372 * LT - 629785214.1 / RT)
      K(R172F) = 5.6800000000D+07 &
	   * exp(0.9 * LT - 8339000 / RT)
      K(R172B) = 1.6601154348D+11 &
	   * exp(0.303461 * LT - 131730378.7 / RT)
      K(R173F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4247000 / RT)
      K(R173B) = 5.5766810183D+11 &
	   * exp(-0.801138 * LT - 369286088.3 / RT)
      K0TROE = 1.5800000000D+48 &
	   * exp(-9.3 * LT - 409195000 / RT)
      KINFTROE = 8.0000000000D+12 &
	   * exp(0.44 * LT - 363046000 / RT)
      FCTROE = 0.2655 * EXP( -TEMP / 180 ) &
	   + 0.7345 * EXP( -TEMP / 1035 ) &
	   + 1 * EXP( -5417 / TEMP )
      K(R174F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM27) )
      K0TROE = 2.9984874433D+42 &
	   * exp(-8.98601 * LT - 228965993.2 / RT)
      KINFTROE = 1.5182214903D+07 &
	   * exp(0.753985 * LT - 182816993.2 / RT)
      K(R174B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM27) )
      K(R175F) = 8.4000000000D+08 * exp(-16213000 / RT)
      K(R175B) = 6.4762207428D+08 &
	   * exp(0.109475 * LT - 69543182.98 / RT)
      K(R176F) = 3.2000000000D+09 * exp(-3573000 / RT)
      K(R176B) = 9.8839182305D-02 &
	   * exp(1.82258 * LT - 357650240.4 / RT)
      K(R177F) = 1.0000000000D+10
      K(R177B) = 5.0507727073D+01 &
	   * exp(1.5678 * LT - 344881549.9 / RT)
      K(R284) = 3.3700000000D+10
      K(R285F) = 6.7000000000D+03 &
	   * exp(1.83 * LT - 920000 / RT)
      K(R285B) = 5.2184108813D+04 &
	   * exp(1.40463 * LT - 58648455.6 / RT)
      K(R286F) = 1.0960000000D+11
      K(R286B) = 1.1538240882D+14 &
	   * exp(-0.512869 * LT - 317161135.9 / RT)
      K(R287F) = 5.0000000000D+12 * exp(-72509000 / RT)
      K(R287B) = 8.9968949343D+12 &
	   * exp(0.257238 * LT - 365775456.9 / RT)
      K(R288) = 8.0000000000D+06 &
	   * exp(0.5 * LT + 7343000 / RT)
      K0TROE = 4.8200000000D+19 &
	   * exp(-2.8 * LT - 2469000 / RT)
      KINFTROE = 1.9700000000D+09 &
	   * exp(0.43 * LT + 1548000 / RT)
      FCTROE = 0.422 * EXP( -TEMP / 122 ) &
	   + 0.578 * EXP( -TEMP / 2535 ) &
	   + 1 * EXP( -9365 / TEMP )
      K(R289F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM37) )
      K0TROE = 8.1072637985D+25 &
	   * exp(-3.26159 * LT - 455241880.8 / RT)
      KINFTROE = 3.3135497268D+15 &
	   * exp(-0.0315917 * LT - 451224880.8 / RT)
      K(R289B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM37) )
      K(R290) = 5.8000000000D+09 * exp(-6276000 / RT)
      K(R291F) = 2.4000000000D+09 * exp(-6276000 / RT)
      K(R291B) = 4.6980064580D+11 &
	   * exp(-0.323258 * LT - 261431734.9 / RT)
      K(R292) = 2.0000000000D+11 * exp(-45978000 / RT)
      K(R293) = 6.8200000000D+07 &
	   * exp(0.25 * LT + 3912000 / RT)
      K(R294F) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 46000 / RT)
      K(R294B) = 1.8710324658D+10 &
	   * exp(-0.149379 * LT - 26471246.48 / RT)
      K(R295F) = 1.3370000000D+03 &
	   * exp(1.61 * LT + 1607000 / RT)
      K(R295B) = 8.6917200738D+03 &
	   * exp(1.42628 * LT - 56259577.53 / RT)
      K(R296F) = 2.9200000000D+09 * exp(-7565000 / RT)
      K(R296B) = 1.6072755775D+06 &
	   * exp(0.523982 * LT - 21998497.9 / RT)
      K(R297) = 2.9200000000D+09 * exp(-7565000 / RT)
      K(R298) = 3.0100000000D+10 * exp(-163804000 / RT)
      K(R299F) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10063000 / RT)
      K(R299B) = 1.8972913589D+03 &
	   * exp(1.71905 * LT - 32384175.17 / RT)
      K(R300) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10063000 / RT)
      K(R301) = 2.3430000000D+07 &
	   * exp(0.73 * LT + 4657000 / RT)
      K(R302) = 3.0100000000D+09 * exp(-49886000 / RT)
      K(R303) = 2.7200000000D+03 &
	   * exp(1.77 * LT - 24769000 / RT)
      K0TROE = 1.0120000000D+36 &
	   * exp(-7.63 * LT - 16125000 / RT)
      KINFTROE = 4.8650000000D+08 &
	   * exp(0.422 * LT + 7343000 / RT)
      FCTROE = 0.535 * EXP( -TEMP / 201 ) &
	   + 0.465 * EXP( -TEMP / 1773 ) &
	   + 1 * EXP( -5333 / TEMP )
      K(R304F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM38) )
      K0TROE = 1.0322723577D+39 &
	   * exp(-7.59017 * LT - 161297835.4 / RT)
      KINFTROE = 4.9624555534D+11 &
	   * exp(0.46183 * LT - 137829835.4 / RT)
      K(R304B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM38) )
      K(R305) = 1.5000000000D+11
      K(R306) = 1.8100000000D+07
      K(R307) = 2.3500000000D+07
      K(R308F) = 2.2000000000D+10
      K(R308B) = 3.2660500844D+04 &
	   * exp(1.304 * LT - 49727494.87 / RT)
      K(R309F) = 1.1000000000D+10
      K(R309B) = 2.5205680886D+09 &
	   * exp(0.310405 * LT - 287918652.4 / RT)
      K(R310F) = 1.2000000000D+10
      K(R310B) = 3.0493428013D+10 &
	   * exp(0.20557 * LT - 352761436.8 / RT)
      K(R311F) = 3.0100000000D+10
      K(R311B) = 4.2643057978D+08 &
	   * exp(0.331703 * LT - 36173617.81 / RT)
      K0TROE = 2.7100000000D+68 &
	   * exp(-16.82 * LT - 54664000 / RT)
      KINFTROE = 9.4300000000D+09
      FCTROE = 0.8473 * EXP( -TEMP / 291 ) &
	   + 0.1527 * EXP( -TEMP / 2742 ) &
	   + 1 * EXP( -7748 / TEMP )
      K(R312F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM39) )
      K0TROE = 1.6125862272D+81 &
	   * exp(-18.5813 * LT - 430275207.3 / RT)
      KINFTROE = 5.6113240305D+22 &
	   * exp(-1.76126 * LT - 375611207.3 / RT)
      K(R312B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM39) )
      K(R313F) = 1.9300000000D+02 &
	   * exp(2.68 * LT - 15548000 / RT)
      K(R313B) = 1.6501412880D-02 &
	   * exp(3.21792 * LT - 17150464.49 / RT)
      K(R314F) = 1.3200000000D+03 &
	   * exp(2.54 * LT - 28267000 / RT)
      K(R314B) = 1.8976251454D-01 &
	   * exp(3.11299 * LT - 37757141.76 / RT)
      K(R315F) = 3.1600000000D+04 &
	   * exp(1.8 * LT - 3908000 / RT)
      K(R315B) = 5.0378221682D+01 &
	   * exp(2.26815 * LT - 78240926.21 / RT)
      K(R316F) = 3.7800000000D-01 &
	   * exp(2.72 * LT - 6276000 / RT)
      K(R316B) = 8.9168034617D+00 &
	   * exp(2.75397 * LT - 63625001.72 / RT)
      K(R317F) = 9.0300000000D-04 &
	   * exp(3.65 * LT - 29932000 / RT)
      K(R317B) = 2.2615544726D-04 &
	   * exp(3.71143 * LT - 48485854.94 / RT)
      K0TROE = 3.0000000000D+57 &
	   * exp(-14.6 * LT - 76023000 / RT)
      KINFTROE = 2.5500000000D+03 &
	   * exp(1.6 * LT - 23849000 / RT)
      FCTROE = 0.8106 * EXP( -TEMP / 277 ) &
	   + 0.1894 * EXP( -TEMP / 8748 ) &
	   + 1 * EXP( -7891 / TEMP )
      K(R318F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM40) )
      K0TROE = 2.0514723075D+67 &
	   * exp(-16.2598 * LT - 179370493.6 / RT)
      KINFTROE = 1.7437514614D+13 &
	   * exp(-0.0598217 * LT - 127196493.6 / RT)
      K(R318B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM40) )
      K(R319F) = 9.6400000000D+10
      K(R319B) = 3.4984678612D+06 &
	   * exp(1.26275 * LT - 336876485.7 / RT)
      K0TROE = 4.4200000000D+55 &
	   * exp(-13.545 * LT - 47518000 / RT)
      KINFTROE = 3.6130000000D+10
      FCTROE = 0.685 * EXP( -TEMP / 369 ) &
	   + 0.315 * EXP( -TEMP / 3285 ) &
	   + 1 * EXP( -6667 / TEMP )
      K(R320F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM41) )
      K0TROE = 7.1862994128D+61 &
	   * exp(-13.7678 * LT - 471119346 / RT)
      KINFTROE = 5.8742307191D+16 &
	   * exp(-0.222751 * LT - 423601346 / RT)
      K(R320B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM41) )
      K(R321F) = 4.0600000000D+03 &
	   * exp(2.19 * LT - 3724000 / RT)
      K(R321B) = 1.1093163569D-03 &
	   * exp(3.72851 * LT - 51714138.7 / RT)
      K(R322F) = 2.4100000000D+10
      K(R322B) = 6.2838833067D+07 &
	   * exp(0.566209 * LT - 34436261.64 / RT)
      K(R323F) = 2.5500000000D+07 &
	   * exp(0.255 * LT + 3946000 / RT)
      K(R323B) = 2.8781079189D+10 &
	   * exp(0.0440876 * LT - 214987530.7 / RT)
      K(R324) = 2.4100000000D+10
      K(R325F) = 1.9270000000D+10 * exp(-0.32 * LT)
      K(R325B) = 8.0515241135D+08 &
	   * exp(0.0420217 * LT - 11152761 / RT)

      W(R1F) = K(R1F) * C(SO) * C(SO) * M(MM1)
      W(R1B) = K(R1B) * C(SO2) * M(MM1)
      W(R2F) = K(R2F) * C(SO) * C(SH) * M(MM2)
      W(R2B) = K(R2B) * C(SOH) * M(MM2)
      W(R3F) = K(R3F) * C(SO) * C(SH2)
      W(R3B) = K(R3B) * C(SOH) * C(SH)
      W(R4F) = K(R4F) * C(SO) * C(SHO2)
      W(R4B) = K(R4B) * C(SO2) * C(SOH)
      W(R5F) = K(R5F) * C(SO) * C(SH2O2)
      W(R5B) = K(R5B) * C(SHO2) * C(SOH)
      W(R6F) = K(R6F) * C(SO) * C(SCH)
      W(R6B) = K(R6B) * C(SCO) * C(SH)
      W(R7F) = K(R7F) * C(SO) * C(SCH2)
      W(R7B) = K(R7B) * C(SHCO) * C(SH)
      W(R8F) = K(R8F) * C(SO) * C(SCH2YXCH2)
      W(R8B) = K(R8B) * C(SCO) * C(SH2)
      W(R9F) = K(R9F) * C(SO) * C(SCH2YXCH2)
      W(R9B) = K(R9B) * C(SHCO) * C(SH)
      W(R10F) = K(R10F) * C(SO) * C(SCH3)
      W(R10B) = K(R10B) * C(SCH2O) * C(SH)
      W(R11F) = K(R11F) * C(SO) * C(SCH4)
      W(R11B) = K(R11B) * C(SCH3) * C(SOH)
      W(R12F) = K(R12F) * C(SO) * C(SCO)
      W(R12B) = K(R12B) * C(SCO2)
      W(R13F) = K(R13F) * C(SO) * C(SHCO)
      W(R13B) = K(R13B) * C(SCO) * C(SOH)
      W(R14F) = K(R14F) * C(SO) * C(SHCO)
      W(R14B) = K(R14B) * C(SCO2) * C(SH)
      W(R15F) = K(R15F) * C(SO) * C(SCH2O)
      W(R15B) = K(R15B) * C(SHCO) * C(SOH)
      W(R16F) = K(R16F) * C(SO) * C(SCH2OH)
      W(R16B) = K(R16B) * C(SCH2O) * C(SOH)
      W(R17F) = K(R17F) * C(SO) * C(SCH3O)
      W(R17B) = K(R17B) * C(SCH2O) * C(SOH)
      W(R18F) = K(R18F) * C(SO) * C(SCH3OH)
      W(R18B) = K(R18B) * C(SCH2OH) * C(SOH)
      W(R19F) = K(R19F) * C(SO) * C(SCH3OH)
      W(R19B) = K(R19B) * C(SCH3O) * C(SOH)
      W(R20F) = K(R20F) * C(SO) * C(SC2H)
      W(R20B) = K(R20B) * C(SCO) * C(SCH)
      W(R21F) = K(R21F) * C(SO) * C(SC2H2)
      W(R21B) = K(R21B) * C(SHCCO) * C(SH)
      W(R22F) = K(R22F) * C(SO) * C(SC2H2)
      W(R22B) = K(R22B) * C(SC2H) * C(SOH)
      W(R23F) = K(R23F) * C(SO) * C(SC2H2)
      W(R23B) = K(R23B) * C(SCH2) * C(SCO)
      W(R24F) = K(R24F) * C(SO) * C(SC2H3)
      W(R24B) = K(R24B) * C(SCH2CO) * C(SH)
      W(R25F) = K(R25F) * C(SO) * C(SC2H4)
      W(R25B) = K(R25B) * C(SHCO) * C(SCH3)
      W(R26F) = K(R26F) * C(SO) * C(SC2H5)
      W(R26B) = K(R26B) * C(SCH2O) * C(SCH3)
      W(R27F) = K(R27F) * C(SO) * C(SC2H6)
      W(R27B) = K(R27B) * C(SC2H5) * C(SOH)
      W(R28F) = K(R28F) * C(SO) * C(SHCCO)
      W(R28B) = K(R28B) * C(SCO) * C(SCO) * C(SH)
      W(R29F) = K(R29F) * C(SO) * C(SCH2CO)
      W(R29B) = K(R29B) * C(SHCCO) * C(SOH)
      W(R30F) = K(R30F) * C(SO) * C(SCH2CO)
      W(R30B) = K(R30B) * C(SCO2) * C(SCH2)
      W(R31F) = K(R31F) * C(SO2) * C(SCO)
      W(R31B) = K(R31B) * C(SCO2) * C(SO)
      W(R32F) = K(R32F) * C(SO2) * C(SCH2O)
      W(R32B) = K(R32B) * C(SHCO) * C(SHO2)
      W(R33F) = K(R33F) * C(SH) * C(SO2) * M(MM4)
      W(R33B) = K(R33B) * C(SHO2) * M(MM4)
      W(R34F) = K(R34F) * C(SH) * C(SO2) * C(SO2)
      W(R34B) = K(R34B) * C(SO2) * C(SHO2)
      W(R35F) = K(R35F) * C(SH) * C(SO2) * C(SH2O)
      W(R35B) = K(R35B) * C(SH2O) * C(SHO2)
      W(R36F) = K(R36F) * C(SH) * C(SO2) * C(SN2)
      W(R36B) = K(R36B) * C(SN2) * C(SHO2)
      W(R37F) = K(R37F) * C(SH) * C(SO2) * C(SAR)
      W(R37B) = K(R37B) * C(SAR) * C(SHO2)
      W(R38F) = K(R38F) * C(SH) * C(SO2)
      W(R38B) = K(R38B) * C(SOH) * C(SO)
      W(R39F) = K(R39F) * C(SH) * C(SH) * M(MM5)
      W(R39B) = K(R39B) * C(SH2) * M(MM5)
      W(R40F) = K(R40F) * C(SH) * C(SH) * C(SH2)
      W(R40B) = K(R40B) * C(SH2) * C(SH2)
      W(R41F) = K(R41F) * C(SH) * C(SH) * C(SH2O)
      W(R41B) = K(R41B) * C(SH2O) * C(SH2)
      W(R42F) = K(R42F) * C(SH) * C(SH) * C(SCO2)
      W(R42B) = K(R42B) * C(SCO2) * C(SH2)
      W(R43F) = K(R43F) * C(SH) * C(SOH) * M(MM6)
      W(R43B) = K(R43B) * C(SH2O) * M(MM6)
      W(R44F) = K(R44F) * C(SH) * C(SHO2)
      W(R44B) = K(R44B) * C(SH2O) * C(SO)
      W(R45F) = K(R45F) * C(SH) * C(SHO2)
      W(R45B) = K(R45B) * C(SH2) * C(SO2)
      W(R46F) = K(R46F) * C(SH) * C(SHO2)
      W(R46B) = K(R46B) * C(SOH) * C(SOH)
      W(R47F) = K(R47F) * C(SH) * C(SH2O2)
      W(R47B) = K(R47B) * C(SH2) * C(SHO2)
      W(R48F) = K(R48F) * C(SH) * C(SH2O2)
      W(R48B) = K(R48B) * C(SH2O) * C(SOH)
      W(R49F) = K(R49F) * C(SH) * C(SCH)
      W(R49B) = K(R49B) * C(SH2) * C(SC)
      W(R50F) = K(R50F) * C(SH) * C(SCH2)
      W(R50B) = K(R50B) * C(SCH3)
      W(R51F) = K(R51F) * C(SH) * C(SCH2YXCH2)
      W(R51B) = K(R51B) * C(SH2) * C(SCH)
      W(R52F) = K(R52F) * C(SH) * C(SCH3)
      W(R52B) = K(R52B) * C(SCH4)
      W(R53F) = K(R53F) * C(SH) * C(SCH4)
      W(R53B) = K(R53B) * C(SH2) * C(SCH3)
      W(R54F) = K(R54F) * C(SH) * C(SHCO)
      W(R54B) = K(R54B) * C(SCH2O)
      W(R55F) = K(R55F) * C(SH) * C(SHCO)
      W(R55B) = K(R55B) * C(SCO) * C(SH2)
      W(R56F) = K(R56F) * C(SH) * C(SCH2O)
      W(R56B) = K(R56B) * C(SCH2OH)
      W(R57F) = K(R57F) * C(SH) * C(SCH2O)
      W(R57B) = K(R57B) * C(SCH3O)
      W(R58F) = K(R58F) * C(SH) * C(SCH2O)
      W(R58B) = K(R58B) * C(SH2) * C(SHCO)
      W(R59F) = K(R59F) * C(SH) * C(SCH2OH)
      W(R59B) = K(R59B) * C(SCH3OH)
      W(R60F) = K(R60F) * C(SH) * C(SCH2OH)
      W(R60B) = K(R60B) * C(SCH2O) * C(SH2)
      W(R61F) = K(R61F) * C(SH) * C(SCH2OH)
      W(R61B) = K(R61B) * C(SCH3) * C(SOH)
      W(R62F) = K(R62F) * C(SH) * C(SCH2OH)
      W(R62B) = K(R62B) * C(SH2O) * C(SCH2YXCH2)
      W(R63F) = K(R63F) * C(SH) * C(SCH3O)
      W(R63B) = K(R63B) * C(SCH3OH)
      W(R64F) = K(R64F) * C(SH) * C(SCH3O)
      W(R64B) = K(R64B) * C(SCH2OH) * C(SH)
      W(R65F) = K(R65F) * C(SH) * C(SCH3O)
      W(R65B) = K(R65B) * C(SCH2O) * C(SH2)
      W(R66F) = K(R66F) * C(SH) * C(SCH3O)
      W(R66B) = K(R66B) * C(SCH3) * C(SOH)
      W(R67F) = K(R67F) * C(SH) * C(SCH3O)
      W(R67B) = K(R67B) * C(SH2O) * C(SCH2YXCH2)
      W(R68F) = K(R68F) * C(SH) * C(SCH3OH)
      W(R68B) = K(R68B) * C(SH2) * C(SCH2OH)
      W(R69F) = K(R69F) * C(SH) * C(SCH3OH)
      W(R69B) = K(R69B) * C(SH2) * C(SCH3O)
      W(R70F) = K(R70F) * C(SH) * C(SC2H)
      W(R70B) = K(R70B) * C(SC2H2)
      W(R71F) = K(R71F) * C(SH) * C(SC2H2)
      W(R71B) = K(R71B) * C(SC2H3)
      W(R72F) = K(R72F) * C(SH) * C(SC2H3)
      W(R72B) = K(R72B) * C(SC2H4)
      W(R73F) = K(R73F) * C(SH) * C(SC2H3)
      W(R73B) = K(R73B) * C(SC2H2) * C(SH2)
      W(R74F) = K(R74F) * C(SH) * C(SC2H4)
      W(R74B) = K(R74B) * C(SC2H5)
      W(R75F) = K(R75F) * C(SH) * C(SC2H4)
      W(R75B) = K(R75B) * C(SH2) * C(SC2H3)
      W(R76F) = K(R76F) * C(SH) * C(SC2H5)
      W(R76B) = K(R76B) * C(SC2H6)
      W(R77F) = K(R77F) * C(SH) * C(SC2H5)
      W(R77B) = K(R77B) * C(SC2H4) * C(SH2)
      W(R78F) = K(R78F) * C(SH) * C(SC2H6)
      W(R78B) = K(R78B) * C(SH2) * C(SC2H5)
      W(R79F) = K(R79F) * C(SH) * C(SHCCO)
      W(R79B) = K(R79B) * C(SCO) * C(SCH2YXCH2)
      W(R80F) = K(R80F) * C(SH) * C(SCH2CO)
      W(R80B) = K(R80B) * C(SH2) * C(SHCCO)
      W(R81F) = K(R81F) * C(SH) * C(SCH2CO)
      W(R81B) = K(R81B) * C(SCO) * C(SCH3)
      W(R82F) = K(R82F) * C(SH) * C(SHCCOH)
      W(R82B) = K(R82B) * C(SCH2CO) * C(SH)
      W(R83F) = K(R83F) * C(SH2) * C(SCO)
      W(R83B) = K(R83B) * C(SCH2O)
      W(R84F) = K(R84F) * C(SOH) * C(SH2)
      W(R84B) = K(R84B) * C(SH2O) * C(SH)
      W(R85F) = K(R85F) * C(SOH) * C(SOH)
      W(R85B) = K(R85B) * C(SH2O2)
      W(R86F) = K(R86F) * C(SOH) * C(SOH)
      W(R86B) = K(R86B) * C(SH2O) * C(SO)
      W(R87F) = K(R87F) * C(SOH) * C(SHO2)
      W(R87B) = K(R87B) * C(SH2O) * C(SO2)
      W(R88F) = K(R88F) * C(SOH) * C(SH2O2)
      W(R88B) = K(R88B) * C(SH2O) * C(SHO2)
      W(R89F) = K(R89F) * C(SOH) * C(SH2O2)
      W(R89B) = K(R89B) * C(SH2O) * C(SHO2)
      W(R90F) = K(R90F) * C(SOH) * C(SC)
      W(R90B) = K(R90B) * C(SCO) * C(SH)
      W(R91F) = K(R91F) * C(SOH) * C(SCH)
      W(R91B) = K(R91B) * C(SHCO) * C(SH)
      W(R92F) = K(R92F) * C(SOH) * C(SCH2)
      W(R92B) = K(R92B) * C(SCH2O) * C(SH)
      W(R93F) = K(R93F) * C(SOH) * C(SCH2)
      W(R93B) = K(R93B) * C(SH2O) * C(SCH)
      W(R94F) = K(R94F) * C(SOH) * C(SCH2YXCH2)
      W(R94B) = K(R94B) * C(SCH2O) * C(SH)
      W(R95F) = K(R95F) * C(SOH) * C(SCH3)
      W(R95B) = K(R95B) * C(SCH3OH)
      W(R96F) = K(R96F) * C(SOH) * C(SCH3)
      W(R96B) = K(R96B) * C(SH2O) * C(SCH2)
      W(R97F) = K(R97F) * C(SOH) * C(SCH3)
      W(R97B) = K(R97B) * C(SH2O) * C(SCH2YXCH2)
      W(R98F) = K(R98F) * C(SOH) * C(SCH4)
      W(R98B) = K(R98B) * C(SH2O) * C(SCH3)
      W(R99F) = K(R99F) * C(SOH) * C(SCO)
      W(R99B) = K(R99B) * C(SCO2) * C(SH)
      W(R100F) = K(R100F) * C(SOH) * C(SHCO)
      W(R100B) = K(R100B) * C(SCO) * C(SH2O)
      W(R101F) = K(R101F) * C(SOH) * C(SCH2O)
      W(R101B) = K(R101B) * C(SH2O) * C(SHCO)
      W(R102F) = K(R102F) * C(SOH) * C(SCH2OH)
      W(R102B) = K(R102B) * C(SCH2O) * C(SH2O)
      W(R103F) = K(R103F) * C(SOH) * C(SCH3O)
      W(R103B) = K(R103B) * C(SCH2O) * C(SH2O)
      W(R104F) = K(R104F) * C(SOH) * C(SCH3OH)
      W(R104B) = K(R104B) * C(SH2O) * C(SCH2OH)
      W(R105F) = K(R105F) * C(SOH) * C(SCH3OH)
      W(R105B) = K(R105B) * C(SH2O) * C(SCH3O)
      W(R106F) = K(R106F) * C(SOH) * C(SC2H)
      W(R106B) = K(R106B) * C(SHCCO) * C(SH)
      W(R107F) = K(R107F) * C(SOH) * C(SC2H2)
      W(R107B) = K(R107B) * C(SCH2CO) * C(SH)
      W(R108F) = K(R108F) * C(SOH) * C(SC2H2)
      W(R108B) = K(R108B) * C(SHCCOH) * C(SH)
      W(R109F) = K(R109F) * C(SOH) * C(SC2H2)
      W(R109B) = K(R109B) * C(SH2O) * C(SC2H)
      W(R110F) = K(R110F) * C(SOH) * C(SC2H2)
      W(R110B) = K(R110B) * C(SCO) * C(SCH3)
      W(R111F) = K(R111F) * C(SOH) * C(SC2H3)
      W(R111B) = K(R111B) * C(SC2H2) * C(SH2O)
      W(R112F) = K(R112F) * C(SOH) * C(SC2H4)
      W(R112B) = K(R112B) * C(SH2O) * C(SC2H3)
      W(R113F) = K(R113F) * C(SOH) * C(SC2H6)
      W(R113B) = K(R113B) * C(SH2O) * C(SC2H5)
      W(R114F) = K(R114F) * C(SOH) * C(SCH2CO)
      W(R114B) = K(R114B) * C(SH2O) * C(SHCCO)
      W(R115F) = K(R115F) * C(SHO2) * C(SHO2)
      W(R115B) = K(R115B) * C(SH2O2) * C(SO2)
      W(R116F) = K(R116F) * C(SHO2) * C(SHO2)
      W(R116B) = K(R116B) * C(SH2O2) * C(SO2)
      W(R117F) = K(R117F) * C(SHO2) * C(SCH2)
      W(R117B) = K(R117B) * C(SCH2O) * C(SOH)
      W(R118F) = K(R118F) * C(SHO2) * C(SCH3)
      W(R118B) = K(R118B) * C(SCH4) * C(SO2)
      W(R119F) = K(R119F) * C(SHO2) * C(SCH3)
      W(R119B) = K(R119B) * C(SCH3O) * C(SOH)
      W(R120F) = K(R120F) * C(SHO2) * C(SCO)
      W(R120B) = K(R120B) * C(SCO2) * C(SOH)
      W(R121F) = K(R121F) * C(SHO2) * C(SCH2O)
      W(R121B) = K(R121B) * C(SH2O2) * C(SHCO)
      W(R122F) = K(R122F) * C(SC) * C(SO2)
      W(R122B) = K(R122B) * C(SCO) * C(SO)
      W(R123F) = K(R123F) * C(SC) * C(SCH2)
      W(R123B) = K(R123B) * C(SC2H) * C(SH)
      W(R124F) = K(R124F) * C(SC) * C(SCH3)
      W(R124B) = K(R124B) * C(SC2H2) * C(SH)
      W(R125F) = K(R125F) * C(SCH) * C(SO2)
      W(R125B) = K(R125B) * C(SHCO) * C(SO)
      W(R126F) = K(R126F) * C(SCH) * C(SH2)
      W(R126B) = K(R126B) * C(SCH2) * C(SH)
      W(R127F) = K(R127F) * C(SCH) * C(SH2O)
      W(R127B) = K(R127B) * C(SCH2O) * C(SH)
      W(R128F) = K(R128F) * C(SCH) * C(SCH2)
      W(R128B) = K(R128B) * C(SC2H2) * C(SH)
      W(R129F) = K(R129F) * C(SCH) * C(SCH3)
      W(R129B) = K(R129B) * C(SC2H3) * C(SH)
      W(R130F) = K(R130F) * C(SCH) * C(SCH4)
      W(R130B) = K(R130B) * C(SC2H4) * C(SH)
      W(R131F) = K(R131F) * C(SCH) * C(SCO)
      W(R131B) = K(R131B) * C(SHCCO)
      W(R132F) = K(R132F) * C(SCH) * C(SCO2)
      W(R132B) = K(R132B) * C(SCO) * C(SHCO)
      W(R133F) = K(R133F) * C(SCH) * C(SCH2O)
      W(R133B) = K(R133B) * C(SCH2CO) * C(SH)
      W(R134F) = K(R134F) * C(SCH) * C(SHCCO)
      W(R134B) = K(R134B) * C(SC2H2) * C(SCO)
      W(R135) = K(R135) * C(SCH2) * C(SO2)
      W(R136F) = K(R136F) * C(SCH2) * C(SH2)
      W(R136B) = K(R136B) * C(SCH3) * C(SH)
      W(R137F) = K(R137F) * C(SCH2) * C(SCH2)
      W(R137B) = K(R137B) * C(SC2H2) * C(SH2)
      W(R138F) = K(R138F) * C(SCH2) * C(SCH3)
      W(R138B) = K(R138B) * C(SC2H4) * C(SH)
      W(R139F) = K(R139F) * C(SCH2) * C(SCH4)
      W(R139B) = K(R139B) * C(SCH3) * C(SCH3)
      W(R140F) = K(R140F) * C(SCH2) * C(SCO)
      W(R140B) = K(R140B) * C(SCH2CO)
      W(R141F) = K(R141F) * C(SCH2) * C(SHCCO)
      W(R141B) = K(R141B) * C(SCO) * C(SC2H3)
      W(R142F) = K(R142F) * C(SCH2YXCH2) * C(SN2)
      W(R142B) = K(R142B) * C(SN2) * C(SCH2)
      W(R143F) = K(R143F) * C(SCH2YXCH2) * C(SAR)
      W(R143B) = K(R143B) * C(SAR) * C(SCH2)
      W(R144F) = K(R144F) * C(SCH2YXCH2) * C(SO2)
      W(R144B) = K(R144B) * C(SCO) * C(SOH) * C(SH)
      W(R145F) = K(R145F) * C(SCH2YXCH2) * C(SO2)
      W(R145B) = K(R145B) * C(SH2O) * C(SCO)
      W(R146F) = K(R146F) * C(SCH2YXCH2) * C(SH2)
      W(R146B) = K(R146B) * C(SH) * C(SCH3)
      W(R147F) = K(R147F) * C(SCH2YXCH2) * C(SH2O)
      W(R147B) = K(R147B) * C(SCH3OH)
      W(R148F) = K(R148F) * C(SCH2YXCH2) * C(SH2O)
      W(R148B) = K(R148B) * C(SH2O) * C(SCH2)
      W(R149F) = K(R149F) * C(SCH2YXCH2) * C(SCH3)
      W(R149B) = K(R149B) * C(SC2H4) * C(SH)
      W(R150F) = K(R150F) * C(SCH2YXCH2) * C(SCH4)
      W(R150B) = K(R150B) * C(SCH3) * C(SCH3)
      W(R151F) = K(R151F) * C(SCH2YXCH2) * C(SCO)
      W(R151B) = K(R151B) * C(SCO) * C(SCH2)
      W(R152F) = K(R152F) * C(SCH2YXCH2) * C(SCO2)
      W(R152B) = K(R152B) * C(SCO2) * C(SCH2)
      W(R153F) = K(R153F) * C(SCH2YXCH2) * C(SCO2)
      W(R153B) = K(R153B) * C(SCH2O) * C(SCO)
      W(R154F) = K(R154F) * C(SCH2YXCH2) * C(SC2H6)
      W(R154B) = K(R154B) * C(SC2H5) * C(SCH3)
      W(R155F) = K(R155F) * C(SCH3) * C(SO2)
      W(R155B) = K(R155B) * C(SCH3O) * C(SO)
      W(R156F) = K(R156F) * C(SCH3) * C(SO2)
      W(R156B) = K(R156B) * C(SCH2O) * C(SOH)
      W(R157F) = K(R157F) * C(SCH3) * C(SH2O2)
      W(R157B) = K(R157B) * C(SCH4) * C(SHO2)
      W(R158F) = K(R158F) * C(SCH3) * C(SCH3)
      W(R158B) = K(R158B) * C(SC2H6)
      W(R159F) = K(R159F) * C(SCH3) * C(SCH3)
      W(R159B) = K(R159B) * C(SC2H5) * C(SH)
      W(R160F) = K(R160F) * C(SCH3) * C(SHCO)
      W(R160B) = K(R160B) * C(SCO) * C(SCH4)
      W(R161F) = K(R161F) * C(SCH3) * C(SCH2O)
      W(R161B) = K(R161B) * C(SCH4) * C(SHCO)
      W(R162F) = K(R162F) * C(SCH3) * C(SCH3OH)
      W(R162B) = K(R162B) * C(SCH4) * C(SCH2OH)
      W(R163F) = K(R163F) * C(SCH3) * C(SCH3OH)
      W(R163B) = K(R163B) * C(SCH4) * C(SCH3O)
      W(R164F) = K(R164F) * C(SCH3) * C(SC2H4)
      W(R164B) = K(R164B) * C(SCH4) * C(SC2H3)
      W(R165F) = K(R165F) * C(SCH3) * C(SC2H6)
      W(R165B) = K(R165B) * C(SCH4) * C(SC2H5)
      W(R166F) = K(R166F) * C(SHCO) * C(SH2O)
      W(R166B) = K(R166B) * C(SH2O) * C(SCO) * C(SH)
      W(R167F) = K(R167F) * C(SHCO) * M(MM26)
      W(R167B) = K(R167B) * C(SCO) * C(SH) * M(MM26)
      W(R168F) = K(R168F) * C(SHCO) * C(SO2)
      W(R168B) = K(R168B) * C(SCO) * C(SHO2)
      W(R169F) = K(R169F) * C(SCH2OH) * C(SO2)
      W(R169B) = K(R169B) * C(SCH2O) * C(SHO2)
      W(R170F) = K(R170F) * C(SCH3O) * C(SO2)
      W(R170B) = K(R170B) * C(SCH2O) * C(SHO2)
      W(R171F) = K(R171F) * C(SC2H) * C(SO2)
      W(R171B) = K(R171B) * C(SCO) * C(SHCO)
      W(R172F) = K(R172F) * C(SC2H) * C(SH2)
      W(R172B) = K(R172B) * C(SC2H2) * C(SH)
      W(R173F) = K(R173F) * C(SC2H3) * C(SO2)
      W(R173B) = K(R173B) * C(SCH2O) * C(SHCO)
      W(R174F) = K(R174F) * C(SC2H4)
      W(R174B) = K(R174B) * C(SC2H2) * C(SH2)
      W(R175F) = K(R175F) * C(SC2H5) * C(SO2)
      W(R175B) = K(R175B) * C(SC2H4) * C(SHO2)
      W(R176F) = K(R176F) * C(SHCCO) * C(SO2)
      W(R176B) = K(R176B) * C(SCO) * C(SCO) * C(SOH)
      W(R177F) = K(R177F) * C(SHCCO) * C(SHCCO)
      W(R177B) = K(R177B) * C(SC2H2) * C(SCO) * C(SCO)
      W(R284) = K(R284) * C(SO) * C(SCH3)
      W(R285F) = K(R285F) * C(SO) * C(SC2H4)
      W(R285B) = K(R285B) * C(SCH2CHO) * C(SH)
      W(R286F) = K(R286F) * C(SO) * C(SC2H5)
      W(R286B) = K(R286B) * C(SCH3CHO) * C(SH)
      W(R287F) = K(R287F) * C(SOH) * C(SHO2)
      W(R287B) = K(R287B) * C(SH2O) * C(SO2)
      W(R288) = K(R288) * C(SOH) * C(SCH3)
      W(R289F) = K(R289F) * C(SCH) * C(SH2)
      W(R289B) = K(R289B) * C(SCH3)
      W(R290) = K(R290) * C(SCH2) * C(SO2)
      W(R291F) = K(R291F) * C(SCH2) * C(SO2)
      W(R291B) = K(R291B) * C(SCH2O) * C(SO)
      W(R292) = K(R292) * C(SCH2) * C(SCH2)
      W(R293) = K(R293) * C(SCH2YXCH2) * C(SH2O)
      W(R294F) = K(R294F) * C(SC2H3) * C(SO2)
      W(R294B) = K(R294B) * C(SCH2CHO) * C(SO)
      W(R295F) = K(R295F) * C(SC2H3) * C(SO2)
      W(R295B) = K(R295B) * C(SC2H2) * C(SHO2)
      W(R296F) = K(R296F) * C(SO) * C(SCH3CHO)
      W(R296B) = K(R296B) * C(SCH2CHO) * C(SOH)
      W(R297) = K(R297) * C(SO) * C(SCH3CHO)
      W(R298) = K(R298) * C(SO2) * C(SCH3CHO)
      W(R299F) = K(R299F) * C(SH) * C(SCH3CHO)
      W(R299B) = K(R299B) * C(SH2) * C(SCH2CHO)
      W(R300) = K(R300) * C(SH) * C(SCH3CHO)
      W(R301) = K(R301) * C(SOH) * C(SCH3CHO)
      W(R302) = K(R302) * C(SHO2) * C(SCH3CHO)
      W(R303) = K(R303) * C(SCH3) * C(SCH3CHO)
      W(R304F) = K(R304F) * C(SH) * C(SCH2CO)
      W(R304B) = K(R304B) * C(SCH2CHO)
      W(R305) = K(R305) * C(SO) * C(SCH2CHO)
      W(R306) = K(R306) * C(SO2) * C(SCH2CHO)
      W(R307) = K(R307) * C(SO2) * C(SCH2CHO)
      W(R308F) = K(R308F) * C(SH) * C(SCH2CHO)
      W(R308B) = K(R308B) * C(SHCO) * C(SCH3)
      W(R309F) = K(R309F) * C(SH) * C(SCH2CHO)
      W(R309B) = K(R309B) * C(SH2) * C(SCH2CO)
      W(R310F) = K(R310F) * C(SOH) * C(SCH2CHO)
      W(R310B) = K(R310B) * C(SCH2CO) * C(SH2O)
      W(R311F) = K(R311F) * C(SOH) * C(SCH2CHO)
      W(R311B) = K(R311B) * C(SCH2OH) * C(SHCO)
      W(R312F) = K(R312F) * C(SCH3) * C(SC2H5)
      W(R312B) = K(R312B) * C(SC3H8)
      W(R313F) = K(R313F) * C(SO) * C(SC3H8)
      W(R313B) = K(R313B) * C(SC3H7) * C(SOH)
      W(R314F) = K(R314F) * C(SH) * C(SC3H8)
      W(R314B) = K(R314B) * C(SH2) * C(SC3H7)
      W(R315F) = K(R315F) * C(SOH) * C(SC3H8)
      W(R315B) = K(R315B) * C(SH2O) * C(SC3H7)
      W(R316F) = K(R316F) * C(SC3H7) * C(SH2O2)
      W(R316B) = K(R316B) * C(SC3H8) * C(SHO2)
      W(R317F) = K(R317F) * C(SCH3) * C(SC3H8)
      W(R317B) = K(R317B) * C(SCH4) * C(SC3H7)
      W(R318F) = K(R318F) * C(SCH3) * C(SC2H4)
      W(R318B) = K(R318B) * C(SC3H7)
      W(R319F) = K(R319F) * C(SO) * C(SC3H7)
      W(R319B) = K(R319B) * C(SCH2O) * C(SC2H5)
      W(R320F) = K(R320F) * C(SH) * C(SC3H7)
      W(R320B) = K(R320B) * C(SC3H8)
      W(R321F) = K(R321F) * C(SH) * C(SC3H7)
      W(R321B) = K(R321B) * C(SC2H5) * C(SCH3)
      W(R322F) = K(R322F) * C(SOH) * C(SC3H7)
      W(R322B) = K(R322B) * C(SCH2OH) * C(SC2H5)
      W(R323F) = K(R323F) * C(SHO2) * C(SC3H7)
      W(R323B) = K(R323B) * C(SC3H8) * C(SO2)
      W(R324) = K(R324) * C(SHO2) * C(SC3H7)
      W(R325F) = K(R325F) * C(SCH3) * C(SC3H7)
      W(R325B) = K(R325B) * C(SC2H5) * C(SC2H5)


      CDOT(SN2) = - W(R36F) + W(R36F) - W(R36B) & 
	    + W(R36B) - W(R142F) + W(R142F) & 
	    - W(R142B) + W(R142B)
      CDOT(SO) = - 2 * W(R1F) + 2 * W(R1B) - W(R2F) & 
	    + W(R2B) - W(R3F) + W(R3B) & 
	    - W(R4F) + W(R4B) - W(R5F) & 
	    + W(R5B) - W(R6F) + W(R6B) & 
	    - W(R7F) + W(R7B) - W(R8F) & 
	    + W(R8B) - W(R9F) + W(R9B) & 
	    - W(R10F) + W(R10B) - W(R11F) & 
	    + W(R11B) - W(R12F) + W(R12B) & 
	    - W(R13F) + W(R13B) - W(R14F) & 
	    + W(R14B) - W(R15F) + W(R15B) & 
	    - W(R16F) + W(R16B) - W(R17F) & 
	    + W(R17B) - W(R18F) + W(R18B) & 
	    - W(R19F) + W(R19B) - W(R20F) & 
	    + W(R20B) - W(R21F) + W(R21B) & 
	    - W(R22F) + W(R22B) - W(R23F)
      CDOT(SO) = CDOT(SO) + W(R23B) - W(R24F) + W(R24B) & 
	    - W(R25F) + W(R25B) - W(R26F) & 
	    + W(R26B) - W(R27F) + W(R27B) & 
	    - W(R28F) + W(R28B) - W(R29F) & 
	    + W(R29B) - W(R30F) + W(R30B) & 
	    + W(R31F) - W(R31B) + W(R38F) & 
	    - W(R38B) + W(R44F) - W(R44B) & 
	    + W(R86F) - W(R86B) + W(R122F) & 
	    - W(R122B) + W(R125F) - W(R125B) & 
	    + W(R155F) - W(R155B) - W(R284) & 
	    - W(R285F) + W(R285B) - W(R286F) & 
	    + W(R286B) + W(R291F) - W(R291B) & 
	    + W(R294F) - W(R294B) - W(R296F) & 
	    + W(R296B) - W(R297) - W(R305) & 
	    - W(R313F) + W(R313B) - W(R319F)
      CDOT(SO) = CDOT(SO) + W(R319B)
      CDOT(SO2) = W(R1F) - W(R1B) + W(R4F) & 
	    - W(R4B) - W(R31F) + W(R31B) & 
	    - W(R32F) + W(R32B) - W(R33F) & 
	    + W(R33B) - 2 * W(R34F) + W(R34F) & 
	    - W(R34B) + 2 * W(R34B) - W(R35F) & 
	    + W(R35B) - W(R36F) + W(R36B) & 
	    - W(R37F) + W(R37B) - W(R38F) & 
	    + W(R38B) + W(R45F) - W(R45B) & 
	    + W(R87F) - W(R87B) + W(R115F) & 
	    - W(R115B) + W(R116F) - W(R116B) & 
	    + W(R118F) - W(R118B) - W(R122F) & 
	    + W(R122B) - W(R125F) + W(R125B) & 
	    - W(R135) - W(R144F) + W(R144B) & 
	    - W(R145F) + W(R145B) - W(R155F) & 
	    + W(R155B) - W(R156F) + W(R156B)
      CDOT(SO2) = CDOT(SO2) - W(R168F) + W(R168B) - W(R169F) & 
	    + W(R169B) - W(R170F) + W(R170B) & 
	    - W(R171F) + W(R171B) - W(R173F) & 
	    + W(R173B) - W(R175F) + W(R175B) & 
	    - W(R176F) + W(R176B) + W(R287F) & 
	    - W(R287B) - W(R290) - W(R291F) & 
	    + W(R291B) - W(R294F) + W(R294B) & 
	    - W(R295F) + W(R295B) - W(R298) & 
	    - W(R306) - W(R307) + W(R323F) & 
	    - W(R323B)
      CDOT(SH) = - W(R2F) + W(R2B) + W(R3F) & 
	    - W(R3B) + W(R6F) - W(R6B) & 
	    + W(R7F) - W(R7B) + W(R9F) & 
	    - W(R9B) + W(R10F) - W(R10B) & 
	    + W(R14F) - W(R14B) + W(R21F) & 
	    - W(R21B) + W(R24F) - W(R24B) & 
	    + W(R28F) - W(R28B) - W(R33F) & 
	    + W(R33B) - W(R34F) + W(R34B) & 
	    - W(R35F) + W(R35B) - W(R36F) & 
	    + W(R36B) - W(R37F) + W(R37B) & 
	    - W(R38F) + W(R38B) - 2 * W(R39F) & 
	    + 2 * W(R39B) - 2 * W(R40F) + 2 * W(R40B) & 
	    - 2 * W(R41F) + 2 * W(R41B) - 2 * W(R42F) & 
	    + 2 * W(R42B) - W(R43F) + W(R43B) & 
	    - W(R44F) + W(R44B) - W(R45F)
      CDOT(SH) = CDOT(SH) + W(R45B) - W(R46F) + W(R46B) & 
	    - W(R47F) + W(R47B) - W(R48F) & 
	    + W(R48B) - W(R49F) + W(R49B) & 
	    - W(R50F) + W(R50B) - W(R51F) & 
	    + W(R51B) - W(R52F) + W(R52B) & 
	    - W(R53F) + W(R53B) - W(R54F) & 
	    + W(R54B) - W(R55F) + W(R55B) & 
	    - W(R56F) + W(R56B) - W(R57F) & 
	    + W(R57B) - W(R58F) + W(R58B) & 
	    - W(R59F) + W(R59B) - W(R60F) & 
	    + W(R60B) - W(R61F) + W(R61B) & 
	    - W(R62F) + W(R62B) - W(R63F) & 
	    + W(R63B) - W(R64F) + W(R64F) & 
	    - W(R64B) + W(R64B) - W(R65F) & 
	    + W(R65B) - W(R66F) + W(R66B)
      CDOT(SH) = CDOT(SH) - W(R67F) + W(R67B) - W(R68F) & 
	    + W(R68B) - W(R69F) + W(R69B) & 
	    - W(R70F) + W(R70B) - W(R71F) & 
	    + W(R71B) - W(R72F) + W(R72B) & 
	    - W(R73F) + W(R73B) - W(R74F) & 
	    + W(R74B) - W(R75F) + W(R75B) & 
	    - W(R76F) + W(R76B) - W(R77F) & 
	    + W(R77B) - W(R78F) + W(R78B) & 
	    - W(R79F) + W(R79B) - W(R80F) & 
	    + W(R80B) - W(R81F) + W(R81B) & 
	    - W(R82F) + W(R82F) - W(R82B) & 
	    + W(R82B) + W(R84F) - W(R84B) & 
	    + W(R90F) - W(R90B) + W(R91F) & 
	    - W(R91B) + W(R92F) - W(R92B) & 
	    + W(R94F) - W(R94B) + W(R99F)
      CDOT(SH) = CDOT(SH) - W(R99B) + W(R106F) - W(R106B) & 
	    + W(R107F) - W(R107B) + W(R108F) & 
	    - W(R108B) + W(R123F) - W(R123B) & 
	    + W(R124F) - W(R124B) + W(R126F) & 
	    - W(R126B) + W(R127F) - W(R127B) & 
	    + W(R128F) - W(R128B) + W(R129F) & 
	    - W(R129B) + W(R130F) - W(R130B) & 
	    + W(R133F) - W(R133B) + W(R135) & 
	    + W(R136F) - W(R136B) + W(R138F) & 
	    - W(R138B) + W(R144F) - W(R144B) & 
	    + W(R146F) - W(R146B) + W(R149F) & 
	    - W(R149B) + W(R159F) - W(R159B) & 
	    + W(R166F) - W(R166B) + W(R167F) & 
	    - W(R167B) + W(R172F) - W(R172B) & 
	    + W(R284) + W(R285F) - W(R285B)
      CDOT(SH) = CDOT(SH) + W(R286F) - W(R286B) + 2 * W(R290) & 
	    + 2 * W(R292) - W(R299F) + W(R299B) & 
	    - W(R300) - W(R304F) + W(R304B) & 
	    + W(R305) - W(R308F) + W(R308B) & 
	    - W(R309F) + W(R309B) - W(R314F) & 
	    + W(R314B) - W(R320F) + W(R320B) & 
	    - W(R321F) + W(R321B)
      CDOT(SOH) = W(R2F) - W(R2B) + W(R3F) & 
	    - W(R3B) + W(R4F) - W(R4B) & 
	    + W(R5F) - W(R5B) + W(R11F) & 
	    - W(R11B) + W(R13F) - W(R13B) & 
	    + W(R15F) - W(R15B) + W(R16F) & 
	    - W(R16B) + W(R17F) - W(R17B) & 
	    + W(R18F) - W(R18B) + W(R19F) & 
	    - W(R19B) + W(R22F) - W(R22B) & 
	    + W(R27F) - W(R27B) + W(R29F) & 
	    - W(R29B) + W(R38F) - W(R38B) & 
	    - W(R43F) + W(R43B) + 2 * W(R46F) & 
	    - 2 * W(R46B) + W(R48F) - W(R48B) & 
	    + W(R61F) - W(R61B) + W(R66F) & 
	    - W(R66B) - W(R84F) + W(R84B) & 
	    - 2 * W(R85F) + 2 * W(R85B) - 2 * W(R86F)
      CDOT(SOH) = CDOT(SOH) + 2 * W(R86B) - W(R87F) + W(R87B) & 
	    - W(R88F) + W(R88B) - W(R89F) & 
	    + W(R89B) - W(R90F) + W(R90B) & 
	    - W(R91F) + W(R91B) - W(R92F) & 
	    + W(R92B) - W(R93F) + W(R93B) & 
	    - W(R94F) + W(R94B) - W(R95F) & 
	    + W(R95B) - W(R96F) + W(R96B) & 
	    - W(R97F) + W(R97B) - W(R98F) & 
	    + W(R98B) - W(R99F) + W(R99B) & 
	    - W(R100F) + W(R100B) - W(R101F) & 
	    + W(R101B) - W(R102F) + W(R102B) & 
	    - W(R103F) + W(R103B) - W(R104F) & 
	    + W(R104B) - W(R105F) + W(R105B) & 
	    - W(R106F) + W(R106B) - W(R107F) & 
	    + W(R107B) - W(R108F) + W(R108B)
      CDOT(SOH) = CDOT(SOH) - W(R109F) + W(R109B) - W(R110F) & 
	    + W(R110B) - W(R111F) + W(R111B) & 
	    - W(R112F) + W(R112B) - W(R113F) & 
	    + W(R113B) - W(R114F) + W(R114B) & 
	    + W(R117F) - W(R117B) + W(R119F) & 
	    - W(R119B) + W(R120F) - W(R120B) & 
	    + W(R135) + W(R144F) - W(R144B) & 
	    + W(R156F) - W(R156B) + W(R176F) & 
	    - W(R176B) - W(R287F) + W(R287B) & 
	    - W(R288) + W(R296F) - W(R296B) & 
	    + W(R297) - W(R301) + W(R306) & 
	    + W(R307) - W(R310F) + W(R310B) & 
	    - W(R311F) + W(R311B) + W(R313F) & 
	    - W(R313B) - W(R315F) + W(R315B) & 
	    - W(R322F) + W(R322B) + W(R324)
      CDOT(SH2) = - W(R3F) + W(R3B) + W(R8F) & 
	    - W(R8B) + W(R39F) - W(R39B) & 
	    - W(R40F) + 2 * W(R40F) - 2 * W(R40B) & 
	    + W(R40B) + W(R41F) - W(R41B) & 
	    + W(R42F) - W(R42B) + W(R45F) & 
	    - W(R45B) + W(R47F) - W(R47B) & 
	    + W(R49F) - W(R49B) + W(R51F) & 
	    - W(R51B) + W(R53F) - W(R53B) & 
	    + W(R55F) - W(R55B) + W(R58F) & 
	    - W(R58B) + W(R60F) - W(R60B) & 
	    + W(R65F) - W(R65B) + W(R68F) & 
	    - W(R68B) + W(R69F) - W(R69B) & 
	    + W(R73F) - W(R73B) + W(R75F) & 
	    - W(R75B) + W(R77F) - W(R77B) & 
	    + W(R78F) - W(R78B) + W(R80F)
      CDOT(SH2) = CDOT(SH2) - W(R80B) - W(R83F) + W(R83B) & 
	    - W(R84F) + W(R84B) - W(R126F) & 
	    + W(R126B) - W(R136F) + W(R136B) & 
	    + W(R137F) - W(R137B) - W(R146F) & 
	    + W(R146B) - W(R172F) + W(R172B) & 
	    + W(R174F) - W(R174B) + W(R284) & 
	    + W(R288) - W(R289F) + W(R289B) & 
	    + W(R293) + W(R299F) - W(R299B) & 
	    + W(R300) + W(R309F) - W(R309B) & 
	    + W(R314F) - W(R314B)
      CDOT(SHO2) = - W(R4F) + W(R4B) + W(R5F) & 
	    - W(R5B) + W(R32F) - W(R32B) & 
	    + W(R33F) - W(R33B) + W(R34F) & 
	    - W(R34B) + W(R35F) - W(R35B) & 
	    + W(R36F) - W(R36B) + W(R37F) & 
	    - W(R37B) - W(R44F) + W(R44B) & 
	    - W(R45F) + W(R45B) - W(R46F) & 
	    + W(R46B) + W(R47F) - W(R47B) & 
	    - W(R87F) + W(R87B) + W(R88F) & 
	    - W(R88B) + W(R89F) - W(R89B) & 
	    - 2 * W(R115F) + 2 * W(R115B) - 2 * W(R116F) & 
	    + 2 * W(R116B) - W(R117F) + W(R117B) & 
	    - W(R118F) + W(R118B) - W(R119F) & 
	    + W(R119B) - W(R120F) + W(R120B) & 
	    - W(R121F) + W(R121B) + W(R157F)
      CDOT(SHO2) = CDOT(SHO2) - W(R157B) + W(R168F) - W(R168B) & 
	    + W(R169F) - W(R169B) + W(R170F) & 
	    - W(R170B) + W(R175F) - W(R175B) & 
	    - W(R287F) + W(R287B) + W(R295F) & 
	    - W(R295B) + W(R298) - W(R302) & 
	    + W(R316F) - W(R316B) - W(R323F) & 
	    + W(R323B) - W(R324)
      CDOT(SH2O2) = - W(R5F) + W(R5B) - W(R47F) & 
	    + W(R47B) - W(R48F) + W(R48B) & 
	    + W(R85F) - W(R85B) - W(R88F) & 
	    + W(R88B) - W(R89F) + W(R89B) & 
	    + W(R115F) - W(R115B) + W(R116F) & 
	    - W(R116B) + W(R121F) - W(R121B) & 
	    - W(R157F) + W(R157B) + W(R302) & 
	    - W(R316F) + W(R316B)
      CDOT(SCH) = - W(R6F) + W(R6B) + W(R20F) & 
	    - W(R20B) - W(R49F) + W(R49B) & 
	    + W(R51F) - W(R51B) - W(R91F) & 
	    + W(R91B) + W(R93F) - W(R93B) & 
	    - W(R125F) + W(R125B) - W(R126F) & 
	    + W(R126B) - W(R127F) + W(R127B) & 
	    - W(R128F) + W(R128B) - W(R129F) & 
	    + W(R129B) - W(R130F) + W(R130B) & 
	    - W(R131F) + W(R131B) - W(R132F) & 
	    + W(R132B) - W(R133F) + W(R133B) & 
	    - W(R134F) + W(R134B) - W(R289F) & 
	    + W(R289B)
      CDOT(SCO) = W(R6F) - W(R6B) + W(R8F) & 
	    - W(R8B) - W(R12F) + W(R12B) & 
	    + W(R13F) - W(R13B) + W(R20F) & 
	    - W(R20B) + W(R23F) - W(R23B) & 
	    + 2 * W(R28F) - 2 * W(R28B) - W(R31F) & 
	    + W(R31B) + W(R55F) - W(R55B) & 
	    + W(R79F) - W(R79B) + W(R81F) & 
	    - W(R81B) - W(R83F) + W(R83B) & 
	    + W(R90F) - W(R90B) - W(R99F) & 
	    + W(R99B) + W(R100F) - W(R100B) & 
	    + W(R110F) - W(R110B) - W(R120F) & 
	    + W(R120B) + W(R122F) - W(R122B) & 
	    - W(R131F) + W(R131B) + W(R132F) & 
	    - W(R132B) + W(R134F) - W(R134B) & 
	    + W(R135) - W(R140F) + W(R140B)
      CDOT(SCO) = CDOT(SCO) + W(R141F) - W(R141B) + W(R144F) & 
	    - W(R144B) + W(R145F) - W(R145B) & 
	    - W(R151F) + W(R151F) - W(R151B) & 
	    + W(R151B) + W(R153F) - W(R153B) & 
	    + W(R160F) - W(R160B) + W(R166F) & 
	    - W(R166B) + W(R167F) - W(R167B) & 
	    + W(R168F) - W(R168B) + W(R171F) & 
	    - W(R171B) + 2 * W(R176F) - 2 * W(R176B) & 
	    + 2 * W(R177F) - 2 * W(R177B) + W(R284) & 
	    + W(R297) + W(R298) + W(R300) & 
	    + W(R301) + W(R302) + W(R303) & 
	    + W(R306)
      CDOT(SCH2) = - W(R7F) + W(R7B) + W(R23F) & 
	    - W(R23B) + W(R30F) - W(R30B) & 
	    - W(R50F) + W(R50B) - W(R92F) & 
	    + W(R92B) - W(R93F) + W(R93B) & 
	    + W(R96F) - W(R96B) - W(R117F) & 
	    + W(R117B) - W(R123F) + W(R123B) & 
	    + W(R126F) - W(R126B) - W(R128F) & 
	    + W(R128B) - W(R135) - W(R136F) & 
	    + W(R136B) - 2 * W(R137F) + 2 * W(R137B) & 
	    - W(R138F) + W(R138B) - W(R139F) & 
	    + W(R139B) - W(R140F) + W(R140B) & 
	    - W(R141F) + W(R141B) + W(R142F) & 
	    - W(R142B) + W(R143F) - W(R143B) & 
	    + W(R148F) - W(R148B) + W(R151F) & 
	    - W(R151B) + W(R152F) - W(R152B)
      CDOT(SCH2) = CDOT(SCH2) - W(R290) - W(R291F) + W(R291B) & 
	    - 2 * W(R292) + W(R305)
      CDOT(SHCO) = W(R7F) - W(R7B) + W(R9F) & 
	    - W(R9B) - W(R13F) + W(R13B) & 
	    - W(R14F) + W(R14B) + W(R15F) & 
	    - W(R15B) + W(R25F) - W(R25B) & 
	    + W(R32F) - W(R32B) - W(R54F) & 
	    + W(R54B) - W(R55F) + W(R55B) & 
	    + W(R58F) - W(R58B) + W(R91F) & 
	    - W(R91B) - W(R100F) + W(R100B) & 
	    + W(R101F) - W(R101B) + W(R121F) & 
	    - W(R121B) + W(R125F) - W(R125B) & 
	    + W(R132F) - W(R132B) - W(R160F) & 
	    + W(R160B) + W(R161F) - W(R161B) & 
	    - W(R166F) + W(R166B) - W(R167F) & 
	    + W(R167B) - W(R168F) + W(R168B) & 
	    + W(R171F) - W(R171B) + W(R173F)
      CDOT(SHCO) = CDOT(SHCO) - W(R173B) + 2 * W(R307) + W(R308F) & 
	    - W(R308B) + W(R311F) - W(R311B)
      CDOT(SCH2YXCH2) = - W(R8F) + W(R8B) - W(R9F) & 
	    + W(R9B) - W(R51F) + W(R51B) & 
	    + W(R62F) - W(R62B) + W(R67F) & 
	    - W(R67B) + W(R79F) - W(R79B) & 
	    - W(R94F) + W(R94B) + W(R97F) & 
	    - W(R97B) - W(R142F) + W(R142B) & 
	    - W(R143F) + W(R143B) - W(R144F) & 
	    + W(R144B) - W(R145F) + W(R145B) & 
	    - W(R146F) + W(R146B) - W(R147F) & 
	    + W(R147B) - W(R148F) + W(R148B) & 
	    - W(R149F) + W(R149B) - W(R150F) & 
	    + W(R150B) - W(R151F) + W(R151B) & 
	    - W(R152F) + W(R152B) - W(R153F) & 
	    + W(R153B) - W(R154F) + W(R154B) & 
	    - W(R293)
      CDOT(SCH3) = - W(R10F) + W(R10B) + W(R11F) & 
	    - W(R11B) + W(R25F) - W(R25B) & 
	    + W(R26F) - W(R26B) + W(R50F) & 
	    - W(R50B) - W(R52F) + W(R52B) & 
	    + W(R53F) - W(R53B) + W(R61F) & 
	    - W(R61B) + W(R66F) - W(R66B) & 
	    + W(R81F) - W(R81B) - W(R95F) & 
	    + W(R95B) - W(R96F) + W(R96B) & 
	    - W(R97F) + W(R97B) + W(R98F) & 
	    - W(R98B) + W(R110F) - W(R110B) & 
	    - W(R118F) + W(R118B) - W(R119F) & 
	    + W(R119B) - W(R124F) + W(R124B) & 
	    - W(R129F) + W(R129B) + W(R136F) & 
	    - W(R136B) - W(R138F) + W(R138B) & 
	    + 2 * W(R139F) - 2 * W(R139B) + W(R146F)
      CDOT(SCH3) = CDOT(SCH3) - W(R146B) - W(R149F) + W(R149B) & 
	    + 2 * W(R150F) - 2 * W(R150B) + W(R154F) & 
	    - W(R154B) - W(R155F) + W(R155B) & 
	    - W(R156F) + W(R156B) - W(R157F) & 
	    + W(R157B) - 2 * W(R158F) + 2 * W(R158B) & 
	    - 2 * W(R159F) + 2 * W(R159B) - W(R160F) & 
	    + W(R160B) - W(R161F) + W(R161B) & 
	    - W(R162F) + W(R162B) - W(R163F) & 
	    + W(R163B) - W(R164F) + W(R164B) & 
	    - W(R165F) + W(R165B) - W(R284) & 
	    - W(R288) + W(R289F) - W(R289B) & 
	    + W(R297) + W(R298) + W(R300) & 
	    + W(R301) + W(R302) - W(R303) & 
	    + W(R303) + W(R308F) - W(R308B) & 
	    - W(R312F) + W(R312B) - W(R317F)
      CDOT(SCH3) = CDOT(SCH3) + W(R317B) - W(R318F) + W(R318B) & 
	    + W(R321F) - W(R321B) - W(R325F) & 
	    + W(R325B)
      CDOT(SCH2O) = W(R10F) - W(R10B) - W(R15F) & 
	    + W(R15B) + W(R16F) - W(R16B) & 
	    + W(R17F) - W(R17B) + W(R26F) & 
	    - W(R26B) - W(R32F) + W(R32B) & 
	    + W(R54F) - W(R54B) - W(R56F) & 
	    + W(R56B) - W(R57F) + W(R57B) & 
	    - W(R58F) + W(R58B) + W(R60F) & 
	    - W(R60B) + W(R65F) - W(R65B) & 
	    + W(R83F) - W(R83B) + W(R92F) & 
	    - W(R92B) + W(R94F) - W(R94B) & 
	    - W(R101F) + W(R101B) + W(R102F) & 
	    - W(R102B) + W(R103F) - W(R103B) & 
	    + W(R117F) - W(R117B) - W(R121F) & 
	    + W(R121B) + W(R127F) - W(R127B) & 
	    - W(R133F) + W(R133B) + W(R153F)
      CDOT(SCH2O) = CDOT(SCH2O) - W(R153B) + W(R156F) - W(R156B) & 
	    - W(R161F) + W(R161B) + W(R169F) & 
	    - W(R169B) + W(R170F) - W(R170B) & 
	    + W(R173F) - W(R173B) + W(R288) & 
	    + W(R291F) - W(R291B) + W(R293) & 
	    + W(R306) + W(R319F) - W(R319B) & 
	    + W(R324)
      CDOT(SCH4) = - W(R11F) + W(R11B) + W(R52F) & 
	    - W(R52B) - W(R53F) + W(R53B) & 
	    - W(R98F) + W(R98B) + W(R118F) & 
	    - W(R118B) - W(R130F) + W(R130B) & 
	    - W(R139F) + W(R139B) - W(R150F) & 
	    + W(R150B) + W(R157F) - W(R157B) & 
	    + W(R160F) - W(R160B) + W(R161F) & 
	    - W(R161B) + W(R162F) - W(R162B) & 
	    + W(R163F) - W(R163B) + W(R164F) & 
	    - W(R164B) + W(R165F) - W(R165B) & 
	    + W(R303) + W(R317F) - W(R317B)
      CDOT(SCO2) = W(R12F) - W(R12B) + W(R14F) & 
	    - W(R14B) + W(R30F) - W(R30B) & 
	    + W(R31F) - W(R31B) - W(R42F) & 
	    + W(R42F) - W(R42B) + W(R42B) & 
	    + W(R99F) - W(R99B) + W(R120F) & 
	    - W(R120B) - W(R132F) + W(R132B) & 
	    - W(R152F) + W(R152F) - W(R152B) & 
	    + W(R152B) - W(R153F) + W(R153B) & 
	    + W(R290) + W(R305)
      CDOT(SCH2OH) = - W(R16F) + W(R16B) + W(R18F) & 
	    - W(R18B) + W(R56F) - W(R56B) & 
	    - W(R59F) + W(R59B) - W(R60F) & 
	    + W(R60B) - W(R61F) + W(R61B) & 
	    - W(R62F) + W(R62B) + W(R64F) & 
	    - W(R64B) + W(R68F) - W(R68B) & 
	    - W(R102F) + W(R102B) + W(R104F) & 
	    - W(R104B) + W(R162F) - W(R162B) & 
	    - W(R169F) + W(R169B) + W(R311F) & 
	    - W(R311B) + W(R322F) - W(R322B)
      CDOT(SCH3O) = - W(R17F) + W(R17B) + W(R19F) & 
	    - W(R19B) + W(R57F) - W(R57B) & 
	    - W(R63F) + W(R63B) - W(R64F) & 
	    + W(R64B) - W(R65F) + W(R65B) & 
	    - W(R66F) + W(R66B) - W(R67F) & 
	    + W(R67B) + W(R69F) - W(R69B) & 
	    - W(R103F) + W(R103B) + W(R105F) & 
	    - W(R105B) + W(R119F) - W(R119B) & 
	    + W(R155F) - W(R155B) + W(R163F) & 
	    - W(R163B) - W(R170F) + W(R170B)
      CDOT(SCH3OH) = - W(R18F) + W(R18B) - W(R19F) & 
	    + W(R19B) + W(R59F) - W(R59B) & 
	    + W(R63F) - W(R63B) - W(R68F) & 
	    + W(R68B) - W(R69F) + W(R69B) & 
	    + W(R95F) - W(R95B) - W(R104F) & 
	    + W(R104B) - W(R105F) + W(R105B) & 
	    + W(R147F) - W(R147B) - W(R162F) & 
	    + W(R162B) - W(R163F) + W(R163B)
      CDOT(SC2H) = - W(R20F) + W(R20B) + W(R22F) & 
	    - W(R22B) - W(R70F) + W(R70B) & 
	    - W(R106F) + W(R106B) + W(R109F) & 
	    - W(R109B) + W(R123F) - W(R123B) & 
	    - W(R171F) + W(R171B) - W(R172F) & 
	    + W(R172B)
      CDOT(SC2H2) = - W(R21F) + W(R21B) - W(R22F) & 
	    + W(R22B) - W(R23F) + W(R23B) & 
	    + W(R70F) - W(R70B) - W(R71F) & 
	    + W(R71B) + W(R73F) - W(R73B) & 
	    - W(R107F) + W(R107B) - W(R108F) & 
	    + W(R108B) - W(R109F) + W(R109B) & 
	    - W(R110F) + W(R110B) + W(R111F) & 
	    - W(R111B) + W(R124F) - W(R124B) & 
	    + W(R128F) - W(R128B) + W(R134F) & 
	    - W(R134B) + W(R137F) - W(R137B) & 
	    + W(R172F) - W(R172B) + W(R174F) & 
	    - W(R174B) + W(R177F) - W(R177B) & 
	    + W(R292) + W(R295F) - W(R295B)
      CDOT(SHCCO) = W(R21F) - W(R21B) - W(R28F) & 
	    + W(R28B) + W(R29F) - W(R29B) & 
	    - W(R79F) + W(R79B) + W(R80F) & 
	    - W(R80B) + W(R106F) - W(R106B) & 
	    + W(R114F) - W(R114B) + W(R131F) & 
	    - W(R131B) - W(R134F) + W(R134B) & 
	    - W(R141F) + W(R141B) - W(R176F) & 
	    + W(R176B) - 2 * W(R177F) + 2 * W(R177B)
      CDOT(SC2H3) = - W(R24F) + W(R24B) + W(R71F) & 
	    - W(R71B) - W(R72F) + W(R72B) & 
	    - W(R73F) + W(R73B) + W(R75F) & 
	    - W(R75B) - W(R111F) + W(R111B) & 
	    + W(R112F) - W(R112B) + W(R129F) & 
	    - W(R129B) + W(R141F) - W(R141B) & 
	    + W(R164F) - W(R164B) - W(R173F) & 
	    + W(R173B) - W(R294F) + W(R294B) & 
	    - W(R295F) + W(R295B)
      CDOT(SCH2CO) = W(R24F) - W(R24B) - W(R29F) & 
	    + W(R29B) - W(R30F) + W(R30B) & 
	    - W(R80F) + W(R80B) - W(R81F) & 
	    + W(R81B) + W(R82F) - W(R82B) & 
	    + W(R107F) - W(R107B) - W(R114F) & 
	    + W(R114B) + W(R133F) - W(R133B) & 
	    + W(R140F) - W(R140B) - W(R304F) & 
	    + W(R304B) + W(R309F) - W(R309B) & 
	    + W(R310F) - W(R310B)
      CDOT(SC2H4) = - W(R25F) + W(R25B) + W(R72F) & 
	    - W(R72B) - W(R74F) + W(R74B) & 
	    - W(R75F) + W(R75B) + W(R77F) & 
	    - W(R77B) - W(R112F) + W(R112B) & 
	    + W(R130F) - W(R130B) + W(R138F) & 
	    - W(R138B) + W(R149F) - W(R149B) & 
	    - W(R164F) + W(R164B) - W(R174F) & 
	    + W(R174B) + W(R175F) - W(R175B) & 
	    - W(R285F) + W(R285B) - W(R318F) & 
	    + W(R318B)
      CDOT(SC2H5) = - W(R26F) + W(R26B) + W(R27F) & 
	    - W(R27B) + W(R74F) - W(R74B) & 
	    - W(R76F) + W(R76B) - W(R77F) & 
	    + W(R77B) + W(R78F) - W(R78B) & 
	    + W(R113F) - W(R113B) + W(R154F) & 
	    - W(R154B) + W(R159F) - W(R159B) & 
	    + W(R165F) - W(R165B) - W(R175F) & 
	    + W(R175B) - W(R286F) + W(R286B) & 
	    - W(R312F) + W(R312B) + W(R319F) & 
	    - W(R319B) + W(R321F) - W(R321B) & 
	    + W(R322F) - W(R322B) + W(R324) & 
	    + 2 * W(R325F) - 2 * W(R325B)
      CDOT(SC2H6) = - W(R27F) + W(R27B) + W(R76F) & 
	    - W(R76B) - W(R78F) + W(R78B) & 
	    - W(R113F) + W(R113B) - W(R154F) & 
	    + W(R154B) + W(R158F) - W(R158B) & 
	    - W(R165F) + W(R165B)
      CDOT(SH2O) = - W(R35F) + W(R35F) - W(R35B) & 
	    + W(R35B) - W(R41F) + W(R41F) & 
	    - W(R41B) + W(R41B) + W(R43F) & 
	    - W(R43B) + W(R44F) - W(R44B) & 
	    + W(R48F) - W(R48B) + W(R62F) & 
	    - W(R62B) + W(R67F) - W(R67B) & 
	    + W(R84F) - W(R84B) + W(R86F) & 
	    - W(R86B) + W(R87F) - W(R87B) & 
	    + W(R88F) - W(R88B) + W(R89F) & 
	    - W(R89B) + W(R93F) - W(R93B) & 
	    + W(R96F) - W(R96B) + W(R97F) & 
	    - W(R97B) + W(R98F) - W(R98B) & 
	    + W(R100F) - W(R100B) + W(R101F) & 
	    - W(R101B) + W(R102F) - W(R102B) & 
	    + W(R103F) - W(R103B) + W(R104F)
      CDOT(SH2O) = CDOT(SH2O) - W(R104B) + W(R105F) - W(R105B) & 
	    + W(R109F) - W(R109B) + W(R111F) & 
	    - W(R111B) + W(R112F) - W(R112B) & 
	    + W(R113F) - W(R113B) + W(R114F) & 
	    - W(R114B) - W(R127F) + W(R127B) & 
	    + W(R145F) - W(R145B) - W(R147F) & 
	    + W(R147B) - W(R148F) + W(R148F) & 
	    - W(R148B) + W(R148B) - W(R166F) & 
	    + W(R166F) - W(R166B) + W(R166B) & 
	    + W(R287F) - W(R287B) - W(R293) & 
	    + W(R301) + W(R310F) - W(R310B) & 
	    + W(R315F) - W(R315B)
      CDOT(SAR) = - W(R37F) + W(R37F) - W(R37B) & 
	    + W(R37B) - W(R143F) + W(R143F) & 
	    - W(R143B) + W(R143B)
      CDOT(SC) = W(R49F) - W(R49B) - W(R90F) & 
	    + W(R90B) - W(R122F) + W(R122B) & 
	    - W(R123F) + W(R123B) - W(R124F) & 
	    + W(R124B)
      CDOT(SHCCOH) = - W(R82F) + W(R82B) + W(R108F) & 
	    - W(R108B)
      CDOT(SCH2CHO) = W(R285F) - W(R285B) + W(R294F) & 
	    - W(R294B) + W(R296F) - W(R296B) & 
	    + W(R299F) - W(R299B) + W(R304F) & 
	    - W(R304B) - W(R305) - W(R306) & 
	    - W(R307) - W(R308F) + W(R308B) & 
	    - W(R309F) + W(R309B) - W(R310F) & 
	    + W(R310B) - W(R311F) + W(R311B)
      CDOT(SCH3CHO) = W(R286F) - W(R286B) - W(R296F) & 
	    + W(R296B) - W(R297) - W(R298) & 
	    - W(R299F) + W(R299B) - W(R300) & 
	    - W(R301) - W(R302) - W(R303)
      CDOT(SC3H8) = W(R312F) - W(R312B) - W(R313F) & 
	    + W(R313B) - W(R314F) + W(R314B) & 
	    - W(R315F) + W(R315B) + W(R316F) & 
	    - W(R316B) - W(R317F) + W(R317B) & 
	    + W(R320F) - W(R320B) + W(R323F) & 
	    - W(R323B)
      CDOT(SC3H7) = W(R313F) - W(R313B) + W(R314F) & 
	    - W(R314B) + W(R315F) - W(R315B) & 
	    - W(R316F) + W(R316B) + W(R317F) & 
	    - W(R317B) + W(R318F) - W(R318B) & 
	    - W(R319F) + W(R319B) - W(R320F) & 
	    + W(R320B) - W(R321F) + W(R321B) & 
	    - W(R322F) + W(R322B) - W(R323F) & 
	    + W(R323B) - W(R324) - W(R325F) & 
	    + W(R325B)
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

      DOUBLE PRECISION MM(36)
      INCLUDE 'grimech30F90.h'

      MM(SN2) =  2.80200000e+01
      MM(SO) =  1.60000000e+01
      MM(SO2) =  3.20000000e+01
      MM(SH) =  1.00800000e+00
      MM(SOH) =  1.70080000e+01
      MM(SH2) =  2.01600000e+00
      MM(SHO2) =  3.30080000e+01
      MM(SH2O2) =  3.40160000e+01
      MM(SCH) =  1.30180000e+01
      MM(SCO) =  2.80100000e+01
      MM(SCH2) =  1.40260000e+01
      MM(SHCO) =  2.90180000e+01
      MM(SCH2YXCH2) =  1.40260000e+01
      MM(SCH3) =  1.50340000e+01
      MM(SCH2O) =  3.00260000e+01
      MM(SCH4) =  1.60420000e+01
      MM(SCO2) =  4.40100000e+01
      MM(SCH2OH) =  3.10340000e+01
      MM(SCH3O) =  3.10340000e+01
      MM(SCH3OH) =  3.20420000e+01
      MM(SC2H) =  2.50280000e+01
      MM(SC2H2) =  2.60360000e+01
      MM(SHCCO) =  4.10280000e+01
      MM(SC2H3) =  2.70440000e+01
      MM(SCH2CO) =  4.20360000e+01
      MM(SC2H4) =  2.80520000e+01
      MM(SC2H5) =  2.90600000e+01
      MM(SC2H6) =  3.00680000e+01
      MM(SH2O) =  1.80160000e+01
      MM(SAR) =  3.99480000e+01
      MM(SC) =  1.20100000e+01
      MM(SHCCOH) =  4.20360000e+01
      MM(SCH2CHO) =  4.30440000e+01
      MM(SCH3CHO) =  4.40520000e+01
      MM(SC3H8) =  4.40940000e+01
      MM(SC3H7) =  4.30860000e+01

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(36)
      INCLUDE 'grimech30F90.h'

      NAMES(SN2)='N2                  '
      NAMES(SO)='O                   '
      NAMES(SO2)='O2                  '
      NAMES(SH)='H                   '
      NAMES(SOH)='OH                  '
      NAMES(SH2)='H2                  '
      NAMES(SHO2)='HO2                 '
      NAMES(SH2O2)='H2O2                '
      NAMES(SCH)='CH                  '
      NAMES(SCO)='CO                  '
      NAMES(SCH2)='CH2                 '
      NAMES(SHCO)='HCO                 '
      NAMES(SCH2YXCH2)='CH2*-CH2            '
      NAMES(SCH3)='CH3                 '
      NAMES(SCH2O)='CH2O                '
      NAMES(SCH4)='CH4                 '
      NAMES(SCO2)='CO2                 '
      NAMES(SCH2OH)='CH2OH               '
      NAMES(SCH3O)='CH3O                '
      NAMES(SCH3OH)='CH3OH               '
      NAMES(SC2H)='C2H                 '
      NAMES(SC2H2)='C2H2                '
      NAMES(SHCCO)='HCCO                '
      NAMES(SC2H3)='C2H3                '
      NAMES(SCH2CO)='CH2CO               '
      NAMES(SC2H4)='C2H4                '
      NAMES(SC2H5)='C2H5                '
      NAMES(SC2H6)='C2H6                '
      NAMES(SH2O)='H2O                 '
      NAMES(SAR)='AR                  '
      NAMES(SC)='C                   '
      NAMES(SHCCOH)='HCCOH               '
      NAMES(SCH2CHO)='CH2CHO              '
      NAMES(SCH3CHO)='CH3CHO              '
      NAMES(SC3H8)='C3H8                '
      NAMES(SC3H7)='C3H7                '

      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'grimech30F90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 422

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'grimech30F90.h'

      NSPECIES_NONS = 36
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'grimech30F90.h'
      real(DP) :: MUCOEFF(36)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SH2O2) =  1.30193418e-06
      MUCOEFF(SCH) =  1.27351520e-06
      MUCOEFF(SCO) =  1.06039632e-06
      MUCOEFF(SCH2) =  6.92304430e-07
      MUCOEFF(SHCO) =  1.11568663e-06
      MUCOEFF(SCH2YXCH2) =  6.92304430e-07
      MUCOEFF(SCH3) =  7.16749611e-07
      MUCOEFF(SCH2O) =  1.13489904e-06
      MUCOEFF(SCH4) =  7.61887935e-07
      MUCOEFF(SCO2) =  1.25056029e-06
      MUCOEFF(SCH2OH) =  1.09210283e-06
      MUCOEFF(SCH3O) =  1.09210283e-06
      MUCOEFF(SCH3OH) =  1.14921582e-06
      MUCOEFF(SC2H) =  7.94406422e-07
      MUCOEFF(SC2H2) =  8.10245830e-07
      MUCOEFF(SHCCO) =  2.73563116e-06
      MUCOEFF(SC2H3) =  8.25781476e-07
      MUCOEFF(SCH2CO) =  1.09806251e-06
      MUCOEFF(SC2H4) =  8.96560348e-07
      MUCOEFF(SC2H5) =  7.77507128e-07
      MUCOEFF(SC2H6) =  7.90876817e-07
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SAR) =  1.52144564e-06
      MUCOEFF(SC) =  8.50486820e-07
      MUCOEFF(SHCCOH) =  1.09806251e-06
      MUCOEFF(SCH2CHO) =  1.11114998e-06
      MUCOEFF(SCH3CHO) =  1.12408509e-06
      MUCOEFF(SC3H8) =  7.14133965e-07
      MUCOEFF(SC3H7) =  7.05924132e-07

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'grimech30F90.h'
      real(DP) :: KOVEREPS(36)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SH2O2) =  9.31098696e-03
      KOVEREPS(SCH) =  1.25000000e-02
      KOVEREPS(SCO) =  1.01936799e-02
      KOVEREPS(SCH2) =  6.94444444e-03
      KOVEREPS(SHCO) =  2.00803213e-03
      KOVEREPS(SCH2YXCH2) =  6.94444444e-03
      KOVEREPS(SCH3) =  6.94444444e-03
      KOVEREPS(SCH2O) =  2.00803213e-03
      KOVEREPS(SCH4) =  7.07213579e-03
      KOVEREPS(SCO2) =  4.09836066e-03
      KOVEREPS(SCH2OH) =  2.39808153e-03
      KOVEREPS(SCH3O) =  2.39808153e-03
      KOVEREPS(SCH3OH) =  2.07555002e-03
      KOVEREPS(SC2H) =  4.78468900e-03
      KOVEREPS(SC2H2) =  4.78468900e-03
      KOVEREPS(SHCCO) =  6.66666667e-03
      KOVEREPS(SC2H3) =  4.78468900e-03
      KOVEREPS(SCH2CO) =  2.29357798e-03
      KOVEREPS(SC2H4) =  3.56125356e-03
      KOVEREPS(SC2H5) =  3.96353547e-03
      KOVEREPS(SC2H6) =  3.96353547e-03
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SAR) =  7.32600733e-03
      KOVEREPS(SC) =  1.40056022e-02
      KOVEREPS(SHCCOH) =  2.29357798e-03
      KOVEREPS(SCH2CHO) =  2.29357798e-03
      KOVEREPS(SCH3CHO) =  2.29357798e-03
      KOVEREPS(SC3H8) =  3.74812594e-03
      KOVEREPS(SC3H7) =  3.74812594e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'grimech30F90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(422), C(36), M(32)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 36
!------------------------------------------------------------------
      implicit none
      include 'grimech30F90.h'
      real(DP) :: H(36), CP(36), T

      IF (T.GT.1000.0_DP) THEN
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  2.92664000e+00_DP + T * (  7.43988400e-04_DP &
	   + T * ( -1.89492000e-07_DP + T * (  2.52425950e-11_DP &
	   + T * ( -1.35067020e-15_DP ) ) ) ) ) -9.22797700e+02_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    2.92664000e+00_DP + T * (  1.48797680e-03_DP &
	   + T * ( -5.68476000e-07_DP + T * (  1.00970380e-10_DP &
	   + T * ( -6.75335100e-15_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.56942078e+00_DP + T * ( -4.29870569e-05_DP &
	   + T * (  1.39828196e-08_DP + T * ( -2.50444497e-12_DP &
	   + T * (  2.45667382e-16_DP ) ) ) ) ) +  2.92175791e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.56942078e+00_DP + T * ( -8.59741137e-05_DP &
	   + T * (  4.19484589e-08_DP + T * ( -1.00177799e-11_DP &
	   + T * (  1.22833691e-15_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.28253784e+00_DP + T * (  7.41543770e-04_DP &
	   + T * ( -2.52655556e-07_DP + T * (  5.23676387e-11_DP &
	   + T * ( -4.33435588e-15_DP ) ) ) ) ) -1.08845772e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.28253784e+00_DP + T * (  1.48308754e-03_DP &
	   + T * ( -7.57966669e-07_DP + T * (  2.09470555e-10_DP &
	   + T * ( -2.16717794e-14_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000001e+00_DP + T * ( -1.15421486e-11_DP &
	   + T * (  5.38539827e-15_DP + T * ( -1.18378809e-18_DP &
	   + T * (  9.96394714e-23_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000001e+00_DP + T * ( -2.30842973e-11_DP &
	   + T * (  1.61561948e-14_DP + T * ( -4.73515235e-18_DP &
	   + T * (  4.98197357e-22_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  3.09288767e+00_DP + T * (  2.74214858e-04_DP &
	   + T * (  4.21684093e-08_DP + T * ( -2.19865389e-11_DP &
	   + T * (  2.34824752e-15_DP ) ) ) ) ) +  3.85865700e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    3.09288767e+00_DP + T * (  5.48429716e-04_DP &
	   + T * (  1.26505228e-07_DP + T * ( -8.79461556e-11_DP &
	   + T * (  1.17412376e-14_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.33727920e+00_DP + T * ( -2.47012365e-05_DP &
	   + T * (  1.66485593e-07_DP + T * ( -4.48915985e-11_DP &
	   + T * (  4.00510752e-15_DP ) ) ) ) ) -9.50158922e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.33727920e+00_DP + T * ( -4.94024731e-05_DP &
	   + T * (  4.99456778e-07_DP + T * ( -1.79566394e-10_DP &
	   + T * (  2.00255376e-14_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.01721090e+00_DP + T * (  1.11991006e-03_DP &
	   + T * ( -2.11219383e-07_DP + T * (  2.85615925e-11_DP &
	   + T * ( -2.15817070e-15_DP ) ) ) ) ) +  1.11856713e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.01721090e+00_DP + T * (  2.23982013e-03_DP &
	   + T * ( -6.33658150e-07_DP + T * (  1.14246370e-10_DP &
	   + T * ( -1.07908535e-14_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.16500285e+00_DP + T * (  2.45415847e-03_DP &
	   + T * ( -6.33797417e-07_DP + T * (  9.27964965e-11_DP &
	   + T * ( -5.75816610e-15_DP ) ) ) ) ) -1.78617877e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.16500285e+00_DP + T * (  4.90831694e-03_DP &
	   + T * ( -1.90139225e-06_DP + T * (  3.71185986e-10_DP &
	   + T * ( -2.87908305e-14_DP ) ) ) ) )
      H(SCH) =  6.38680289e+02_DP * ( &
	   T * (  2.87846473e+00_DP + T * (  4.85456840e-04_DP &
	   + T * (  4.81485517e-08_DP + T * ( -3.26719623e-11_DP &
	   + T * (  3.52158766e-15_DP ) ) ) ) ) +  7.10124364e+04_DP )
      CP(SCH) =  6.38680289e+02_DP * ( &
	    2.87846473e+00_DP + T * (  9.70913681e-04_DP &
	   + T * (  1.44445655e-07_DP + T * ( -1.30687849e-10_DP &
	   + T * (  1.76079383e-14_DP ) ) ) ) )
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  2.71518561e+00_DP + T * (  1.03126372e-03_DP &
	   + T * ( -3.32941924e-07_DP + T * (  5.75132520e-11_DP &
	   + T * ( -4.07295432e-15_DP ) ) ) ) ) -1.41518724e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    2.71518561e+00_DP + T * (  2.06252743e-03_DP &
	   + T * ( -9.98825771e-07_DP + T * (  2.30053008e-10_DP &
	   + T * ( -2.03647716e-14_DP ) ) ) ) )
      H(SCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.87410113e+00_DP + T * (  1.82819646e-03_DP &
	   + T * ( -4.69648657e-07_DP + T * (  6.50448872e-11_DP &
	   + T * ( -3.75455134e-15_DP ) ) ) ) ) +  4.62636040e+04_DP )
      CP(SCH2) =  5.92780550e+02_DP * ( &
	    2.87410113e+00_DP + T * (  3.65639292e-03_DP &
	   + T * ( -1.40894597e-06_DP + T * (  2.60179549e-10_DP &
	   + T * ( -1.87727567e-14_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  2.77217438e+00_DP + T * (  2.47847763e-03_DP &
	   + T * ( -8.28152043e-07_DP + T * (  1.47290445e-10_DP &
	   + T * ( -1.06701742e-14_DP ) ) ) ) ) +  4.01191815e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    2.77217438e+00_DP + T * (  4.95695526e-03_DP &
	   + T * ( -2.48445613e-06_DP + T * (  5.89161778e-10_DP &
	   + T * ( -5.33508711e-14_DP ) ) ) ) )
      H(SCH2YXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.29203842e+00_DP + T * (  2.32794318e-03_DP &
	   + T * ( -6.70639823e-07_DP + T * (  1.04476500e-10_DP &
	   + T * ( -6.79432730e-15_DP ) ) ) ) ) +  5.09259997e+04_DP )
      CP(SCH2YXCH2) =  5.92780550e+02_DP * ( &
	    2.29203842e+00_DP + T * (  4.65588637e-03_DP &
	   + T * ( -2.01191947e-06_DP + T * (  4.17906000e-10_DP &
	   + T * ( -3.39716365e-14_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  2.28571772e+00_DP + T * (  3.61995018e-03_DP &
	   + T * ( -9.95714493e-07_DP + T * (  1.48921161e-10_DP &
	   + T * ( -9.34308788e-15_DP ) ) ) ) ) +  1.67755843e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    2.28571772e+00_DP + T * (  7.23990037e-03_DP &
	   + T * ( -2.98714348e-06_DP + T * (  5.95684644e-10_DP &
	   + T * ( -4.67154394e-14_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  1.76069008e+00_DP + T * (  4.60000041e-03_DP &
	   + T * ( -1.47419604e-06_DP + T * (  2.51603030e-10_DP &
	   + T * ( -1.76771128e-14_DP ) ) ) ) ) -1.39958323e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    1.76069008e+00_DP + T * (  9.20000082e-03_DP &
	   + T * ( -4.42258813e-06_DP + T * (  1.00641212e-09_DP &
	   + T * ( -8.83855640e-14_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  7.48514950e-02_DP + T * (  6.69547335e-03_DP &
	   + T * ( -1.91095270e-06_DP + T * (  3.05731338e-10_DP &
	   + T * ( -2.03630460e-14_DP ) ) ) ) ) -9.46834459e+03_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    7.48514950e-02_DP + T * (  1.33909467e-02_DP &
	   + T * ( -5.73285809e-06_DP + T * (  1.22292535e-09_DP &
	   + T * ( -1.01815230e-13_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  3.85746029e+00_DP + T * (  2.20718513e-03_DP &
	   + T * ( -7.38271347e-07_DP + T * (  1.30872547e-10_DP &
	   + T * ( -9.44168328e-15_DP ) ) ) ) ) -4.87591660e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    3.85746029e+00_DP + T * (  4.41437026e-03_DP &
	   + T * ( -2.21481404e-06_DP + T * (  5.23490188e-10_DP &
	   + T * ( -4.72084164e-14_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  3.69266569e+00_DP + T * (  4.32288399e-03_DP &
	   + T * ( -1.25033707e-06_DP + T * (  1.96808659e-10_DP &
	   + T * ( -1.29710840e-14_DP ) ) ) ) ) -3.24250627e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    3.69266569e+00_DP + T * (  8.64576797e-03_DP &
	   + T * ( -3.75101120e-06_DP + T * (  7.87234636e-10_DP &
	   + T * ( -6.48554201e-14_DP ) ) ) ) )
      H(SCH3O) =  2.67910679e+02_DP * ( &
	   T * (  3.77079900e+00_DP + T * (  3.93574850e-03_DP &
	   + T * ( -8.85461333e-07_DP + T * (  9.86107750e-11_DP &
	   + T * ( -4.22523200e-15_DP ) ) ) ) ) +  1.27832520e+02_DP )
      CP(SCH3O) =  2.67910679e+02_DP * ( &
	    3.77079900e+00_DP + T * (  7.87149700e-03_DP &
	   + T * ( -2.65638400e-06_DP + T * (  3.94443100e-10_DP &
	   + T * ( -2.11261600e-14_DP ) ) ) ) )
      H(SCH3OH) =  2.59482554e+02_DP * ( &
	   T * (  1.78970791e+00_DP + T * (  7.04691460e-03_DP &
	   + T * ( -2.12166945e-06_DP + T * (  3.45427713e-10_DP &
	   + T * ( -2.34120440e-14_DP ) ) ) ) ) -2.53748747e+04_DP )
      CP(SCH3OH) =  2.59482554e+02_DP * ( &
	    1.78970791e+00_DP + T * (  1.40938292e-02_DP &
	   + T * ( -6.36500835e-06_DP + T * (  1.38171085e-09_DP &
	   + T * ( -1.17060220e-13_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  3.16780652e+00_DP + T * (  2.37610951e-03_DP &
	   + T * ( -6.12623590e-07_DP + T * (  7.60475630e-11_DP &
	   + T * ( -3.54465540e-15_DP ) ) ) ) ) +  6.71210650e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    3.16780652e+00_DP + T * (  4.75221902e-03_DP &
	   + T * ( -1.83787077e-06_DP + T * (  3.04190252e-10_DP &
	   + T * ( -1.77232770e-14_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  4.14756964e+00_DP + T * (  2.98083332e-03_DP &
	   + T * ( -7.90982840e-07_DP + T * (  1.16853043e-10_DP &
	   + T * ( -7.22470426e-15_DP ) ) ) ) ) +  2.59359992e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    4.14756964e+00_DP + T * (  5.96166664e-03_DP &
	   + T * ( -2.37294852e-06_DP + T * (  4.67412171e-10_DP &
	   + T * ( -3.61235213e-14_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  5.62820580e+00_DP + T * (  2.04267005e-03_DP &
	   + T * ( -5.31151567e-07_DP + T * (  7.15651300e-11_DP &
	   + T * ( -3.88156640e-15_DP ) ) ) ) ) +  1.93272150e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    5.62820580e+00_DP + T * (  4.08534010e-03_DP &
	   + T * ( -1.59345470e-06_DP + T * (  2.86260520e-10_DP &
	   + T * ( -1.94078320e-14_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.01672400e+00_DP + T * (  5.16511460e-03_DP &
	   + T * ( -1.56027450e-06_DP + T * (  2.54408220e-10_DP &
	   + T * ( -1.72521408e-14_DP ) ) ) ) ) +  3.46128739e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.01672400e+00_DP + T * (  1.03302292e-02_DP &
	   + T * ( -4.68082349e-06_DP + T * (  1.01763288e-09_DP &
	   + T * ( -8.62607041e-14_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  4.51129732e+00_DP + T * (  4.50179872e-03_DP &
	   + T * ( -1.38979878e-06_DP + T * (  2.30836470e-10_DP &
	   + T * ( -1.58967640e-14_DP ) ) ) ) ) -7.55105311e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    4.51129732e+00_DP + T * (  9.00359745e-03_DP &
	   + T * ( -4.16939635e-06_DP + T * (  9.23345882e-10_DP &
	   + T * ( -7.94838201e-14_DP ) ) ) ) )
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  2.03611116e+00_DP + T * (  7.32270755e-03_DP &
	   + T * ( -2.23692638e-06_DP + T * (  3.68057308e-10_DP &
	   + T * ( -2.51412122e-14_DP ) ) ) ) ) +  4.93988614e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    2.03611116e+00_DP + T * (  1.46454151e-02_DP &
	   + T * ( -6.71077915e-06_DP + T * (  1.47222923e-09_DP &
	   + T * ( -1.25706061e-13_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  1.95465642e+00_DP + T * (  8.69863610e-03_DP &
	   + T * ( -2.66068889e-06_DP + T * (  4.38044223e-10_DP &
	   + T * ( -2.99283152e-14_DP ) ) ) ) ) +  1.28575200e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    1.95465642e+00_DP + T * (  1.73972722e-02_DP &
	   + T * ( -7.98206668e-06_DP + T * (  1.75217689e-09_DP &
	   + T * ( -1.49641576e-13_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  1.07188150e+00_DP + T * (  1.08426339e-02_DP &
	   + T * ( -3.34186890e-06_DP + T * (  5.53530003e-10_DP &
	   + T * ( -3.80005780e-14_DP ) ) ) ) ) -1.14263932e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    1.07188150e+00_DP + T * (  2.16852677e-02_DP &
	   + T * ( -1.00256067e-05_DP + T * (  2.21412001e-09_DP &
	   + T * ( -1.90002890e-13_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  3.03399249e+00_DP + T * (  1.08845902e-03_DP &
	   + T * ( -5.46908393e-08_DP + T * ( -2.42604967e-11_DP &
	   + T * (  3.36401984e-15_DP ) ) ) ) ) -3.00042971e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    3.03399249e+00_DP + T * (  2.17691804e-03_DP &
	   + T * ( -1.64072518e-07_DP + T * ( -9.70419870e-11_DP &
	   + T * (  1.68200992e-14_DP ) ) ) ) )
      H(SAR) =  2.08129068e+02_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -7.45375000e+02_DP )
      CP(SAR) =  2.08129068e+02_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC) =  6.92284763e+02_DP * ( &
	   T * (  2.49266888e+00_DP + T * (  2.39944642e-05_DP &
	   + T * ( -2.41445007e-08_DP + T * (  9.35727573e-12_DP &
	   + T * ( -9.74555786e-16_DP ) ) ) ) ) +  8.54512953e+04_DP )
      CP(SC) =  6.92284763e+02_DP * ( &
	    2.49266888e+00_DP + T * (  4.79889284e-05_DP &
	   + T * ( -7.24335020e-08_DP + T * (  3.74291029e-11_DP &
	   + T * ( -4.87277893e-15_DP ) ) ) ) )
      H(SHCCOH) =  1.97790941e+02_DP * ( &
	   T * (  5.92382910e+00_DP + T * (  3.39618000e-03_DP &
	   + T * ( -8.55285467e-07_DP + T * (  1.12469603e-10_DP &
	   + T * ( -5.98802020e-15_DP ) ) ) ) ) +  7.26462600e+03_DP )
      CP(SHCCOH) =  1.97790941e+02_DP * ( &
	    5.92382910e+00_DP + T * (  6.79236000e-03_DP &
	   + T * ( -2.56585640e-06_DP + T * (  4.49878410e-10_DP &
	   + T * ( -2.99401010e-14_DP ) ) ) ) )
      H(SCH2CHO) =  1.93159093e+02_DP * ( &
	   T * (  5.97567000e+00_DP + T * (  4.06529550e-03_DP &
	   + T * ( -9.14541333e-07_DP + T * (  1.01757600e-10_DP &
	   + T * ( -4.35203400e-15_DP ) ) ) ) ) +  4.90321800e+02_DP )
      CP(SCH2CHO) =  1.93159093e+02_DP * ( &
	    5.97567000e+00_DP + T * (  8.13059100e-03_DP &
	   + T * ( -2.74362400e-06_DP + T * (  4.07030400e-10_DP &
	   + T * ( -2.17601700e-14_DP ) ) ) ) )
      H(SCH3CHO) =  1.88739217e+02_DP * ( &
	   T * (  5.40411080e+00_DP + T * (  5.86152950e-03_DP &
	   + T * ( -1.40877123e-06_DP + T * (  1.70931128e-10_DP &
	   + T * ( -8.19697260e-15_DP ) ) ) ) ) -2.25931220e+04_DP )
      CP(SCH3CHO) =  1.88739217e+02_DP * ( &
	    5.40411080e+00_DP + T * (  1.17230590e-02_DP &
	   + T * ( -4.22631370e-06_DP + T * (  6.83724510e-10_DP &
	   + T * ( -4.09848630e-14_DP ) ) ) ) )
      H(SC3H8) =  1.88559441e+02_DP * ( &
	   T * (  7.53413680e+00_DP + T * (  9.43611950e-03_DP &
	   + T * ( -2.09061637e-06_DP + T * (  2.28689123e-10_DP &
	   + T * ( -9.56761380e-15_DP ) ) ) ) ) -1.64675160e+04_DP )
      CP(SC3H8) =  1.88559441e+02_DP * ( &
	    7.53413680e+00_DP + T * (  1.88722390e-02_DP &
	   + T * ( -6.27184910e-06_DP + T * (  9.14756490e-10_DP &
	   + T * ( -4.78380690e-14_DP ) ) ) ) )
      H(SC3H7) =  1.92970803e+02_DP * ( &
	   T * (  7.70269870e+00_DP + T * (  8.02210150e-03_DP &
	   + T * ( -1.76110733e-06_DP + T * (  1.90746475e-10_DP &
	   + T * ( -7.87845680e-15_DP ) ) ) ) ) +  8.29843360e+03_DP )
      CP(SC3H7) =  1.92970803e+02_DP * ( &
	    7.70269870e+00_DP + T * (  1.60442030e-02_DP &
	   + T * ( -5.28332200e-06_DP + T * (  7.62985900e-10_DP &
	   + T * ( -3.93922840e-14_DP ) ) ) ) )
      ELSE
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  3.29867700e+00_DP + T * (  7.04120200e-04_DP &
	   + T * ( -1.32107400e-06_DP + T * (  1.41037875e-09_DP &
	   + T * ( -4.88970800e-13_DP ) ) ) ) ) -1.02089990e+03_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    3.29867700e+00_DP + T * (  1.40824040e-03_DP &
	   + T * ( -3.96322200e-06_DP + T * (  5.64151500e-09_DP &
	   + T * ( -2.44485400e-12_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  3.16826710e+00_DP + T * ( -1.63965942e-03_DP &
	   + T * (  2.21435465e-06_DP + T * ( -1.53201656e-09_DP &
	   + T * (  4.22531942e-13_DP ) ) ) ) ) +  2.91222592e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    3.16826710e+00_DP + T * ( -3.27931884e-03_DP &
	   + T * (  6.64306396e-06_DP + T * ( -6.12806624e-09_DP &
	   + T * (  2.11265971e-12_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.78245636e+00_DP + T * ( -1.49836708e-03_DP &
	   + T * (  3.28243400e-06_DP + T * ( -2.42032377e-09_DP &
	   + T * (  6.48745674e-13_DP ) ) ) ) ) -1.06394356e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.78245636e+00_DP + T * ( -2.99673416e-03_DP &
	   + T * (  9.84730201e-06_DP + T * ( -9.68129509e-09_DP &
	   + T * (  3.24372837e-12_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  3.52666409e-13_DP &
	   + T * ( -6.65306547e-16_DP + T * (  5.75204080e-19_DP &
	   + T * ( -1.85546466e-22_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  7.05332819e-13_DP &
	   + T * ( -1.99591964e-15_DP + T * (  2.30081632e-18_DP &
	   + T * ( -9.27732332e-22_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  3.99201543e+00_DP + T * ( -1.20065876e-03_DP &
	   + T * (  1.53931280e-06_DP + T * ( -9.70283332e-10_DP &
	   + T * (  2.72822940e-13_DP ) ) ) ) ) +  3.61508056e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    3.99201543e+00_DP + T * ( -2.40131752e-03_DP &
	   + T * (  4.61793841e-06_DP + T * ( -3.88113333e-09_DP &
	   + T * (  1.36411470e-12_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.34433112e+00_DP + T * (  3.99026037e-03_DP &
	   + T * ( -6.49271700e-06_DP + T * (  5.03930235e-09_DP &
	   + T * ( -1.47522352e-12_DP ) ) ) ) ) -9.17935173e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.34433112e+00_DP + T * (  7.98052075e-03_DP &
	   + T * ( -1.94781510e-05_DP + T * (  2.01572094e-08_DP &
	   + T * ( -7.37611761e-12_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.30179801e+00_DP + T * ( -2.37456025e-03_DP &
	   + T * (  7.05276303e-06_DP + T * ( -6.06909735e-09_DP &
	   + T * (  1.85845025e-12_DP ) ) ) ) ) +  2.94808040e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.30179801e+00_DP + T * ( -4.74912051e-03_DP &
	   + T * (  2.11582891e-05_DP + T * ( -2.42763894e-08_DP &
	   + T * (  9.29225124e-12_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.27611269e+00_DP + T * ( -2.71411208e-04_DP &
	   + T * (  5.57785670e-06_DP + T * ( -5.39427032e-09_DP &
	   + T * (  1.72490873e-12_DP ) ) ) ) ) -1.77025821e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.27611269e+00_DP + T * ( -5.42822417e-04_DP &
	   + T * (  1.67335701e-05_DP + T * ( -2.15770813e-08_DP &
	   + T * (  8.62454363e-12_DP ) ) ) ) )
      H(SCH) =  6.38680289e+02_DP * ( &
	   T * (  3.48981665e+00_DP + T * (  1.61917771e-04_DP &
	   + T * ( -5.62996883e-07_DP + T * (  7.90543317e-10_DP &
	   + T * ( -2.81218134e-13_DP ) ) ) ) ) +  7.07972934e+04_DP )
      CP(SCH) =  6.38680289e+02_DP * ( &
	    3.48981665e+00_DP + T * (  3.23835541e-04_DP &
	   + T * ( -1.68899065e-06_DP + T * (  3.16217327e-09_DP &
	   + T * ( -1.40609067e-12_DP ) ) ) ) )
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  3.57953347e+00_DP + T * ( -3.05176840e-04_DP &
	   + T * (  3.38938110e-07_DP + T * (  2.26751471e-10_DP &
	   + T * ( -1.80884900e-13_DP ) ) ) ) ) -1.43440860e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    3.57953347e+00_DP + T * ( -6.10353680e-04_DP &
	   + T * (  1.01681433e-06_DP + T * (  9.07005884e-10_DP &
	   + T * ( -9.04424499e-13_DP ) ) ) ) )
      H(SCH2) =  5.92780550e+02_DP * ( &
	   T * (  3.76267867e+00_DP + T * (  4.84436072e-04_DP &
	   + T * (  9.31632803e-07_DP + T * ( -9.62727883e-10_DP &
	   + T * (  3.37483438e-13_DP ) ) ) ) ) +  4.60040401e+04_DP )
      CP(SCH2) =  5.92780550e+02_DP * ( &
	    3.76267867e+00_DP + T * (  9.68872143e-04_DP &
	   + T * (  2.79489841e-06_DP + T * ( -3.85091153e-09_DP &
	   + T * (  1.68741719e-12_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  4.22118584e+00_DP + T * ( -1.62196266e-03_DP &
	   + T * (  4.59331487e-06_DP + T * ( -3.32860233e-09_DP &
	   + T * (  8.67537730e-13_DP ) ) ) ) ) +  3.83956496e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    4.22118584e+00_DP + T * ( -3.24392532e-03_DP &
	   + T * (  1.37799446e-05_DP + T * ( -1.33144093e-08_DP &
	   + T * (  4.33768865e-12_DP ) ) ) ) )
      H(SCH2YXCH2) =  5.92780550e+02_DP * ( &
	   T * (  4.19860411e+00_DP + T * ( -1.18330710e-03_DP &
	   + T * (  2.74432073e-06_DP + T * ( -1.67203995e-09_DP &
	   + T * (  3.88629474e-13_DP ) ) ) ) ) +  5.04968163e+04_DP )
      CP(SCH2YXCH2) =  5.92780550e+02_DP * ( &
	    4.19860411e+00_DP + T * ( -2.36661419e-03_DP &
	   + T * (  8.23296220e-06_DP + T * ( -6.68815981e-09_DP &
	   + T * (  1.94314737e-12_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  3.67359040e+00_DP + T * (  1.00547588e-03_DP &
	   + T * (  1.91007285e-06_DP + T * ( -1.71779356e-09_DP &
	   + T * (  5.08771468e-13_DP ) ) ) ) ) +  1.64449988e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    3.67359040e+00_DP + T * (  2.01095175e-03_DP &
	   + T * (  5.73021856e-06_DP + T * ( -6.87117425e-09_DP &
	   + T * (  2.54385734e-12_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  4.79372315e+00_DP + T * ( -4.95416684e-03_DP &
	   + T * (  1.24406669e-05_DP + T * ( -9.48213152e-09_DP &
	   + T * (  2.63545304e-12_DP ) ) ) ) ) -1.43089567e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    4.79372315e+00_DP + T * ( -9.90833369e-03_DP &
	   + T * (  3.73220008e-05_DP + T * ( -3.79285261e-08_DP &
	   + T * (  1.31772652e-11_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  5.14987613e+00_DP + T * ( -6.83548940e-03_DP &
	   + T * (  1.63933533e-05_DP + T * ( -1.21185757e-08_DP &
	   + T * (  3.33387912e-12_DP ) ) ) ) ) -1.02466476e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    5.14987613e+00_DP + T * ( -1.36709788e-02_DP &
	   + T * (  4.91800599e-05_DP + T * ( -4.84743026e-08_DP &
	   + T * (  1.66693956e-11_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  2.35677352e+00_DP + T * (  4.49229839e-03_DP &
	   + T * ( -2.37452090e-06_DP + T * (  6.14797555e-10_DP &
	   + T * ( -2.87399096e-14_DP ) ) ) ) ) -4.83719697e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    2.35677352e+00_DP + T * (  8.98459677e-03_DP &
	   + T * ( -7.12356269e-06_DP + T * (  2.45919022e-09_DP &
	   + T * ( -1.43699548e-13_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  3.86388918e+00_DP + T * (  2.79836152e-03_DP &
	   + T * (  1.97757264e-06_DP + T * ( -2.61330030e-09_DP &
	   + T * (  8.73934556e-13_DP ) ) ) ) ) -3.19391367e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    3.86388918e+00_DP + T * (  5.59672304e-03_DP &
	   + T * (  5.93271791e-06_DP + T * ( -1.04532012e-08_DP &
	   + T * (  4.36967278e-12_DP ) ) ) ) )
      H(SCH3O) =  2.67910679e+02_DP * ( &
	   T * (  2.10620400e+00_DP + T * (  3.60829750e-03_DP &
	   + T * (  1.77949067e-06_DP + T * ( -1.84440900e-09_DP &
	   + T * (  4.15122000e-13_DP ) ) ) ) ) +  9.78601100e+02_DP )
      CP(SCH3O) =  2.67910679e+02_DP * ( &
	    2.10620400e+00_DP + T * (  7.21659500e-03_DP &
	   + T * (  5.33847200e-06_DP + T * ( -7.37763600e-09_DP &
	   + T * (  2.07561000e-12_DP ) ) ) ) )
      H(SCH3OH) =  2.59482554e+02_DP * ( &
	   T * (  5.71539582e+00_DP + T * ( -7.61545645e-03_DP &
	   + T * (  2.17480385e-05_DP + T * ( -1.77701722e-08_DP &
	   + T * (  5.22705396e-12_DP ) ) ) ) ) -2.56427656e+04_DP )
      CP(SCH3OH) =  2.59482554e+02_DP * ( &
	    5.71539582e+00_DP + T * ( -1.52309129e-02_DP &
	   + T * (  6.52441155e-05_DP + T * ( -7.10806889e-08_DP &
	   + T * (  2.61352698e-11_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  2.88965733e+00_DP + T * (  6.70498055e-03_DP &
	   + T * ( -9.49231670e-06_DP + T * (  7.36977613e-09_DP &
	   + T * ( -2.18663022e-12_DP ) ) ) ) ) +  6.68393932e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    2.88965733e+00_DP + T * (  1.34099611e-02_DP &
	   + T * ( -2.84769501e-05_DP + T * (  2.94791045e-08_DP &
	   + T * ( -1.09331511e-11_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  8.08681094e-01_DP + T * (  1.16807815e-02_DP &
	   + T * ( -1.18390605e-05_DP + T * (  7.00381092e-09_DP &
	   + T * ( -1.70014595e-12_DP ) ) ) ) ) +  2.64289807e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    8.08681094e-01_DP + T * (  2.33615629e-02_DP &
	   + T * ( -3.55171815e-05_DP + T * (  2.80152437e-08_DP &
	   + T * ( -8.50072974e-12_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  2.25172140e+00_DP + T * (  8.82751050e-03_DP &
	   + T * ( -7.90970033e-06_DP + T * (  4.31893975e-09_DP &
	   + T * ( -1.01329622e-12_DP ) ) ) ) ) +  2.00594490e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    2.25172140e+00_DP + T * (  1.76550210e-02_DP &
	   + T * ( -2.37291010e-05_DP + T * (  1.72757590e-08_DP &
	   + T * ( -5.06648110e-12_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.21246645e+00_DP + T * (  7.57395810e-04_DP &
	   + T * (  8.64031373e-06_DP + T * ( -8.94144617e-09_DP &
	   + T * (  2.94301746e-12_DP ) ) ) ) ) +  3.48598468e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.21246645e+00_DP + T * (  1.51479162e-03_DP &
	   + T * (  2.59209412e-05_DP + T * ( -3.57657847e-08_DP &
	   + T * (  1.47150873e-11_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  2.13583630e+00_DP + T * (  9.05943605e-03_DP &
	   + T * ( -5.79824913e-06_DP + T * (  2.33599392e-09_DP &
	   + T * ( -4.02915230e-13_DP ) ) ) ) ) -7.04291804e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    2.13583630e+00_DP + T * (  1.81188721e-02_DP &
	   + T * ( -1.73947474e-05_DP + T * (  9.34397568e-09_DP &
	   + T * ( -2.01457615e-12_DP ) ) ) ) )
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  3.95920148e+00_DP + T * ( -3.78526124e-03_DP &
	   + T * (  1.90330097e-05_DP + T * ( -1.72897188e-08_DP &
	   + T * (  5.39768746e-12_DP ) ) ) ) ) +  5.08977593e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    3.95920148e+00_DP + T * ( -7.57052247e-03_DP &
	   + T * (  5.70990292e-05_DP + T * ( -6.91588753e-08_DP &
	   + T * (  2.69884373e-11_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  4.30646568e+00_DP + T * ( -2.09329446e-03_DP &
	   + T * (  1.65714269e-05_DP + T * ( -1.49781651e-08_DP &
	   + T * (  4.61018008e-12_DP ) ) ) ) ) +  1.28416265e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    4.30646568e+00_DP + T * ( -4.18658892e-03_DP &
	   + T * (  4.97142807e-05_DP + T * ( -5.99126606e-08_DP &
	   + T * (  2.30509004e-11_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  4.29142492e+00_DP + T * ( -2.75077135e-03_DP &
	   + T * (  1.99812763e-05_DP + T * ( -1.77116571e-08_DP &
	   + T * (  5.37371542e-12_DP ) ) ) ) ) -1.15222055e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    4.29142492e+00_DP + T * ( -5.50154270e-03_DP &
	   + T * (  5.99438288e-05_DP + T * ( -7.08466285e-08_DP &
	   + T * (  2.68685771e-11_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  4.19864056e+00_DP + T * ( -1.01821705e-03_DP &
	   + T * (  2.17346737e-06_DP + T * ( -1.37199266e-09_DP &
	   + T * (  3.54395634e-13_DP ) ) ) ) ) -3.02937267e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    4.19864056e+00_DP + T * ( -2.03643410e-03_DP &
	   + T * (  6.52040211e-06_DP + T * ( -5.48797062e-09_DP &
	   + T * (  1.77197817e-12_DP ) ) ) ) )
      H(SAR) =  2.08129068e+02_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -7.45375000e+02_DP )
      CP(SAR) =  2.08129068e+02_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC) =  6.92284763e+02_DP * ( &
	   T * (  2.55423955e+00_DP + T * ( -1.60768862e-04_DP &
	   + T * (  2.44597415e-07_DP + T * ( -1.83058722e-10_DP &
	   + T * (  5.33042892e-14_DP ) ) ) ) ) +  8.54438832e+04_DP )
      CP(SC) =  6.92284763e+02_DP * ( &
	    2.55423955e+00_DP + T * ( -3.21537724e-04_DP &
	   + T * (  7.33792245e-07_DP + T * ( -7.32234889e-10_DP &
	   + T * (  2.66521446e-13_DP ) ) ) ) )
      H(SHCCOH) =  1.97790941e+02_DP * ( &
	   T * (  1.24237330e+00_DP + T * (  1.55361005e-02_DP &
	   + T * ( -1.69556213e-05_DP + T * (  1.07842828e-08_DP &
	   + T * ( -2.80291880e-12_DP ) ) ) ) ) +  8.03161430e+03_DP )
      CP(SHCCOH) =  1.97790941e+02_DP * ( &
	    1.24237330e+00_DP + T * (  3.10722010e-02_DP &
	   + T * ( -5.08668640e-05_DP + T * (  4.31371310e-08_DP &
	   + T * ( -1.40145940e-11_DP ) ) ) ) )
      H(SCH2CHO) =  1.93159093e+02_DP * ( &
	   T * (  3.40906200e+00_DP + T * (  5.36928700e-03_DP &
	   + T * (  6.30497333e-07_DP + T * ( -1.78964575e-09_DP &
	   + T * (  5.73477000e-13_DP ) ) ) ) ) +  1.52147660e+03_DP )
      CP(SCH2CHO) =  1.93159093e+02_DP * ( &
	    3.40906200e+00_DP + T * (  1.07385740e-02_DP &
	   + T * (  1.89149200e-06_DP + T * ( -7.15858300e-09_DP &
	   + T * (  2.86738500e-12_DP ) ) ) ) )
      H(SCH3CHO) =  1.88739217e+02_DP * ( &
	   T * (  4.72945950e+00_DP + T * ( -1.59664290e-03_DP &
	   + T * (  1.58449737e-05_DP + T * ( -1.43646527e-08_DP &
	   + T * (  4.38622240e-12_DP ) ) ) ) ) -2.15728780e+04_DP )
      CP(SCH3CHO) =  1.88739217e+02_DP * ( &
	    4.72945950e+00_DP + T * ( -3.19328580e-03_DP &
	   + T * (  4.75349210e-05_DP + T * ( -5.74586110e-08_DP &
	   + T * (  2.19311120e-11_DP ) ) ) ) )
      H(SC3H8) =  1.88559441e+02_DP * ( &
	   T * (  9.33553810e-01_DP + T * (  1.32122895e-02_DP &
	   + T * (  2.03532423e-06_DP + T * ( -5.49437475e-09_DP &
	   + T * (  1.90298506e-12_DP ) ) ) ) ) -1.39585200e+04_DP )
      CP(SC3H8) =  1.88559441e+02_DP * ( &
	    9.33553810e-01_DP + T * (  2.64245790e-02_DP &
	   + T * (  6.10597270e-06_DP + T * ( -2.19774990e-08_DP &
	   + T * (  9.51492530e-12_DP ) ) ) ) )
      H(SC3H7) =  1.92970803e+02_DP * ( &
	   T * (  1.05155180e+00_DP + T * (  1.29959900e-02_DP &
	   + T * (  7.93351333e-07_DP + T * ( -4.90239225e-09_DP &
	   + T * (  1.87464940e-12_DP ) ) ) ) ) +  1.06318630e+04_DP )
      CP(SC3H7) =  1.92970803e+02_DP * ( &
	    1.05155180e+00_DP + T * (  2.59919800e-02_DP &
	   + T * (  2.38005400e-06_DP + T * ( -1.96095690e-08_DP &
	   + T * (  9.37324700e-12_DP ) ) ) ) )
      END IF

      END
