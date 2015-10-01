!----------------------------------------------------------
! ======= TheSootF.f90 =======
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
      include 'TheSootF90.h'
      real(DP) :: CDOT(158), W(1806), K(1806), &
      C(158), M(10), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM1) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SHO2) &
	    + C(SH2O2) + C(SCO) &
	    + C(SHCO) + C(SC) &
	    + C(SCH) + C(STXCH2) &
	    + C(SCH3) + C(SCH2O) &
	    + C(SHCCO) + C(SC2H) &
	    + C(SCH2CO) + C(SC2H2) &
	    + C(SSXCH2) + 0.63 * C(SAR) &
	    + C(SCH3OH) + C(SCH2OH) &
	    + C(SCH3O) + 2 * C(SCH4) &
	    + C(SCH3O2) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SHCCOH) + C(SCH2CHO)
      M(MM1) = M(MM1) + C(SCH3CHO) + C(SH2C2) &
	    + C(SC2H5O) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SC3H8) &
	    + C(SC3H6) + C(SC3H3) &
	    + C(SPXC3H4) + C(SAXC3H4) &
	    + C(SSXC3H5) + C(SNXC4H3) &
	    + C(SC2H3CHO) + C(SAXC3H5) &
	    + C(SC2O) + C(SC4H4) &
	    + C(SC3H2) + C(SC3H2O) &
	    + C(SC4H2) + C(SIXC4H3) &
	    + C(STXC3H5) + C(SC3H5O) &
	    + C(SC4H) + C(SC8H2) &
	    + C(SC6H2) + C(SC4H6) &
	    + C(SNXC4H5) + C(SIXC4H5) &
	    + C(SA1XC6H6) + C(SNXC7H16)
      M(MM1) = M(MM1) + C(SC5H11) + C(SPXC4H9) &
	    + C(SC7H15) + C(SPXC4H8) &
	    + C(SC5H10) + C(SC7H14) &
	    + C(SC7H15O) + C(SC3H7CHO) &
	    + C(SC4H7) + C(SC7H13) &
	    + C(SC5H9) + C(SC4H7O) &
	    + C(SNXC3H7O) + C(SIXC8H18) &
	    + C(SYXC7H15) + C(SIXC4H8) &
	    + C(SIXC3H7) + C(STXC4H9) &
	    + C(SCXC8H17) + C(SYXC7H14) &
	    + C(SDXC8H17O) + C(SCH3COCH3) &
	    + C(SIXC4H7) + C(SXXC7H13) &
	    + C(SIXC3H5CH) + C(STXC4H9O) &
	    + C(SIXC4H7O) + C(SC5H4CH2) &
	    + C(SA1XXC6H5) + C(SA1C2H2XC)
      M(MM1) = M(MM1) + C(SA1C2H3XC) + C(SA1C2HXC8) &
	    + C(SA1C2HYXC) + C(SA1C2H3YX) &
	    + C(SA2XXC10H) + C(SA2XC10H8) &
	    + C(SA2YXC10H) + C(SA2C2H2AX) &
	    + C(SA2C2H2BX) + C(SA2C2HAXC) &
	    + C(SA2C2HBXC) + C(SA2C2HAYX) &
	    + C(SA2C2HBYX) + C(SA2R5XC12) &
	    + C(SA2R5XXC1) + C(SA2R5C2H2) &
	    + C(SA2R5C2HX) + C(SA2R5C2HY) &
	    + C(SP2XC12H1) + C(SP2XXC12H) &
	    + C(SA3XXC14H) + C(SA3XC14H1) &
	    + C(SA3YXC14H) + C(SA3R5XXC1) &
	    + C(SA3R5XC16) + C(SA4XC16H1) &
	    + C(SA4XXC16H) + C(SA4R5XC18) &
	    + C(SFLTNXC16) + C(SC5H6)
      M(MM1) = M(MM1) + C(SC5H5) + C(STXC5H5O) &
	    + C(SC5H4O) + C(SSXC5H5O) &
	    + C(SC9H8) + C(SC9H7) &
	    + C(SA1CH2XC7) + C(SC9H6O) &
	    + C(SOXC6H4) + C(SA1CH3XC7) &
	    + C(SA1OHXC6H) + C(SHOA1CH3X) &
	    + C(SOA1CH3XC) + C(SA1CH2OXC) &
	    + C(SA1CH2OHX) + C(SA1CHOXC7) &
	    + C(SA1OXC6H5) + C(SA1CH3YXC) &
	    + C(SA1C2H4XC) + C(SA1C2H5XC) &
	    + C(SC8H9O2) + C(SC8H8OOH) &
	    + C(SOC8H7OOH) + C(SA1CH3CH3) &
	    + C(SA1CH3CH2) + C(SA1CH3CHO) &
	    + C(SA2CH3XC1) + C(SA1CHOCH2) &
	    + C(SA1CHOCHO) + C(SA2OHXC10)
      M(MM1) = M(MM1) + C(SA2CH2XC1) + C(SA2CH2OXC) &
	    + C(SA2CHOXC1) + C(SA2OXC10H) &
	    + C(SOC6H4O)
      M(MM2) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6.3 * C(SH2O) + 3.6 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.38 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 2 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM2) = M(MM2) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM2) = M(MM2) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM2) = M(MM2) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM2) = M(MM2) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM2) = M(MM2) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM3) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 12 * C(SH2O) + 3.6 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.7 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 2 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM3) = M(MM3) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM3) = M(MM3) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM3) = M(MM3) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM3) = M(MM3) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM3) = M(MM3) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM4) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.4 * C(SH2) &
	    + 15.4 * C(SH2O) + 3.6 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.83 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 2 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM4) = M(MM4) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM4) = M(MM4) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM4) = M(MM4) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM4) = M(MM4) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM4) = M(MM4) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM5) = C(SN2) + C(SH) &
	    + 0.85 * C(SO2) + C(SO) &
	    + C(SOH) + 0.75 * C(SH2) &
	    + 11.89 * C(SH2O) + 2.18 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.09 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.4 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM5) = M(MM5) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM5) = M(MM5) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM5) = M(MM5) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM5) = M(MM5) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM5) = M(MM5) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM7) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 12 * C(SH2O) + 3.6 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.7 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 2 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM7) = M(MM7) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM7) = M(MM7) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM7) = M(MM7) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM7) = M(MM7) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM7) = M(MM7) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM8) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 3.6 * C(SCO2) + C(SHO2) &
	    + C(SH2O2) + 1.75 * C(SCO) &
	    + C(SHCO) + C(SC) &
	    + C(SCH) + C(STXCH2) &
	    + C(SCH3) + C(SCH2O) &
	    + C(SHCCO) + C(SC2H) &
	    + C(SCH2CO) + C(SC2H2) &
	    + C(SSXCH2) + C(SAR) &
	    + C(SCH3OH) + C(SCH2OH) &
	    + C(SCH3O) + 2 * C(SCH4) &
	    + C(SCH3O2) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5)
      M(MM8) = M(MM8) + C(SHCCOH) + C(SCH2CHO) &
	    + C(SCH3CHO) + C(SH2C2) &
	    + C(SC2H5O) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SC3H8) &
	    + C(SC3H6) + C(SC3H3) &
	    + C(SPXC3H4) + C(SAXC3H4) &
	    + C(SSXC3H5) + C(SNXC4H3) &
	    + C(SC2H3CHO) + C(SAXC3H5) &
	    + C(SC2O) + C(SC4H4) &
	    + C(SC3H2) + C(SC3H2O) &
	    + C(SC4H2) + C(SIXC4H3) &
	    + C(STXC3H5) + C(SC3H5O) &
	    + C(SC4H) + C(SC8H2) &
	    + C(SC6H2) + C(SC4H6) &
	    + C(SNXC4H5) + C(SIXC4H5)
      M(MM8) = M(MM8) + C(SA1XC6H6) + C(SNXC7H16) &
	    + C(SC5H11) + C(SPXC4H9) &
	    + C(SC7H15) + C(SPXC4H8) &
	    + C(SC5H10) + C(SC7H14) &
	    + C(SC7H15O) + C(SC3H7CHO) &
	    + C(SC4H7) + C(SC7H13) &
	    + C(SC5H9) + C(SC4H7O) &
	    + C(SNXC3H7O) + C(SIXC8H18) &
	    + C(SYXC7H15) + C(SIXC4H8) &
	    + C(SIXC3H7) + C(STXC4H9) &
	    + C(SCXC8H17) + C(SYXC7H14) &
	    + C(SDXC8H17O) + C(SCH3COCH3) &
	    + C(SIXC4H7) + C(SXXC7H13) &
	    + C(SIXC3H5CH) + C(STXC4H9O) &
	    + C(SIXC4H7O) + C(SC5H4CH2)
      M(MM8) = M(MM8) + C(SA1XXC6H5) + C(SA1C2H2XC) &
	    + C(SA1C2H3XC) + C(SA1C2HXC8) &
	    + C(SA1C2HYXC) + C(SA1C2H3YX) &
	    + C(SA2XXC10H) + C(SA2XC10H8) &
	    + C(SA2YXC10H) + C(SA2C2H2AX) &
	    + C(SA2C2H2BX) + C(SA2C2HAXC) &
	    + C(SA2C2HBXC) + C(SA2C2HAYX) &
	    + C(SA2C2HBYX) + C(SA2R5XC12) &
	    + C(SA2R5XXC1) + C(SA2R5C2H2) &
	    + C(SA2R5C2HX) + C(SA2R5C2HY) &
	    + C(SP2XC12H1) + C(SP2XXC12H) &
	    + C(SA3XXC14H) + C(SA3XC14H1) &
	    + C(SA3YXC14H) + C(SA3R5XXC1) &
	    + C(SA3R5XC16) + C(SA4XC16H1) &
	    + C(SA4XXC16H) + C(SA4R5XC18)
      M(MM8) = M(MM8) + C(SFLTNXC16) + C(SC5H6) &
	    + C(SC5H5) + C(STXC5H5O) &
	    + C(SC5H4O) + C(SSXC5H5O) &
	    + C(SC9H8) + C(SC9H7) &
	    + C(SA1CH2XC7) + C(SC9H6O) &
	    + C(SOXC6H4) + C(SA1CH3XC7) &
	    + C(SA1OHXC6H) + C(SHOA1CH3X) &
	    + C(SOA1CH3XC) + C(SA1CH2OXC) &
	    + C(SA1CH2OHX) + C(SA1CHOXC7) &
	    + C(SA1OXC6H5) + C(SA1CH3YXC) &
	    + C(SA1C2H4XC) + C(SA1C2H5XC) &
	    + C(SC8H9O2) + C(SC8H8OOH) &
	    + C(SOC8H7OOH) + C(SA1CH3CH3) &
	    + C(SA1CH3CH2) + C(SA1CH3CHO) &
	    + C(SA2CH3XC1) + C(SA1CHOCH2)
      M(MM8) = M(MM8) + C(SA1CHOCHO) + C(SA2OHXC10) &
	    + C(SA2CH2XC1) + C(SA2CH2OXC) &
	    + C(SA2CHOXC1) + C(SA2OXC10H) &
	    + C(SOC6H4O)
      M(MM9) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + 2 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.5 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.7 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 3 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM9) = M(MM9) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM9) = M(MM9) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM9) = M(MM9) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM9) = M(MM9) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM9) = M(MM9) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM0) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SH2) &
	    + C(SH2O) + C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM0) = M(MM0) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM0) = M(MM0) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM0) = M(MM0) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM0) = M(MM0) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM0) = M(MM0) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)
      M(MM6) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6 * C(SH2O) + 3.6 * C(SCO2) &
	    + C(SHO2) + C(SH2O2) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SC) + C(SCH) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SHCCO) &
	    + C(SC2H) + C(SCH2CO) &
	    + C(SC2H2) + C(SSXCH2) &
	    + 0.7 * C(SAR) + C(SCH3OH) &
	    + C(SCH2OH) + C(SCH3O) &
	    + 2 * C(SCH4) + C(SCH3O2) &
	    + C(SC2H3) + C(SC2H4)
      M(MM6) = M(MM6) + C(SC2H5) + C(SHCCOH) &
	    + C(SCH2CHO) + C(SCH3CHO) &
	    + C(SH2C2) + C(SC2H5O) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SC3H8) + C(SC3H6) &
	    + C(SC3H3) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SSXC3H5) &
	    + C(SNXC4H3) + C(SC2H3CHO) &
	    + C(SAXC3H5) + C(SC2O) &
	    + C(SC4H4) + C(SC3H2) &
	    + C(SC3H2O) + C(SC4H2) &
	    + C(SIXC4H3) + C(STXC3H5) &
	    + C(SC3H5O) + C(SC4H) &
	    + C(SC8H2) + C(SC6H2) &
	    + C(SC4H6) + C(SNXC4H5)
      M(MM6) = M(MM6) + C(SIXC4H5) + C(SA1XC6H6) &
	    + C(SNXC7H16) + C(SC5H11) &
	    + C(SPXC4H9) + C(SC7H15) &
	    + C(SPXC4H8) + C(SC5H10) &
	    + C(SC7H14) + C(SC7H15O) &
	    + C(SC3H7CHO) + C(SC4H7) &
	    + C(SC7H13) + C(SC5H9) &
	    + C(SC4H7O) + C(SNXC3H7O) &
	    + C(SIXC8H18) + C(SYXC7H15) &
	    + C(SIXC4H8) + C(SIXC3H7) &
	    + C(STXC4H9) + C(SCXC8H17) &
	    + C(SYXC7H14) + C(SDXC8H17O) &
	    + C(SCH3COCH3) + C(SIXC4H7) &
	    + C(SXXC7H13) + C(SIXC3H5CH) &
	    + C(STXC4H9O) + C(SIXC4H7O)
      M(MM6) = M(MM6) + C(SC5H4CH2) + C(SA1XXC6H5) &
	    + C(SA1C2H2XC) + C(SA1C2H3XC) &
	    + C(SA1C2HXC8) + C(SA1C2HYXC) &
	    + C(SA1C2H3YX) + C(SA2XXC10H) &
	    + C(SA2XC10H8) + C(SA2YXC10H) &
	    + C(SA2C2H2AX) + C(SA2C2H2BX) &
	    + C(SA2C2HAXC) + C(SA2C2HBXC) &
	    + C(SA2C2HAYX) + C(SA2C2HBYX) &
	    + C(SA2R5XC12) + C(SA2R5XXC1) &
	    + C(SA2R5C2H2) + C(SA2R5C2HX) &
	    + C(SA2R5C2HY) + C(SP2XC12H1) &
	    + C(SP2XXC12H) + C(SA3XXC14H) &
	    + C(SA3XC14H1) + C(SA3YXC14H) &
	    + C(SA3R5XXC1) + C(SA3R5XC16) &
	    + C(SA4XC16H1) + C(SA4XXC16H)
      M(MM6) = M(MM6) + C(SA4R5XC18) + C(SFLTNXC16) &
	    + C(SC5H6) + C(SC5H5) &
	    + C(STXC5H5O) + C(SC5H4O) &
	    + C(SSXC5H5O) + C(SC9H8) &
	    + C(SC9H7) + C(SA1CH2XC7) &
	    + C(SC9H6O) + C(SOXC6H4) &
	    + C(SA1CH3XC7) + C(SA1OHXC6H) &
	    + C(SHOA1CH3X) + C(SOA1CH3XC) &
	    + C(SA1CH2OXC) + C(SA1CH2OHX) &
	    + C(SA1CHOXC7) + C(SA1OXC6H5) &
	    + C(SA1CH3YXC) + C(SA1C2H4XC) &
	    + C(SA1C2H5XC) + C(SC8H9O2) &
	    + C(SC8H8OOH) + C(SOC8H7OOH) &
	    + C(SA1CH3CH3) + C(SA1CH3CH2) &
	    + C(SA1CH3CHO) + C(SA2CH3XC1)
      M(MM6) = M(MM6) + C(SA1CHOCH2) + C(SA1CHOCHO) &
	    + C(SA2OHXC10) + C(SA2CH2XC1) &
	    + C(SA2CH2OXC) + C(SA2CHOXC1) &
	    + C(SA2OXC10H) + C(SOC6H4O)


      K(R1F) = 2.6400000000D+13 &
	   * exp(-0.67 * LT - 71300000 / RT)
      K(R1B) = 5.2766713024D+10 &
	   * exp(-0.234291 * LT - 479604.2721 / RT)
      K(R2F) = 4.5900000000D+01 &
	   * exp(2.7 * LT - 26190000 / RT)
      K(R2B) = 2.7517349858D+01 &
	   * exp(2.66385 * LT - 20100582.44 / RT)
      K(R3F) = 1.7300000000D+05 &
	   * exp(1.51 * LT - 14350000 / RT)
      K(R3B) = 1.9032628245D+06 &
	   * exp(1.40625 * LT - 77394524.74 / RT)
      K(R4F) = 3.9700000000D+01 &
	   * exp(2.4 * LT + 8830000 / RT)
      K(R4B) = 7.2853303338D+02 &
	   * exp(2.33241 * LT - 60303942.3 / RT)
      K(R5F) = 1.7800000000D+12 * exp(-1 * LT)
      K(R5B) = 4.1604426150D+14 &
	   * exp(-0.649765 * LT - 433091487.8 / RT)
      K(R6F) = 9.0000000000D+10 * exp(-0.6 * LT)
      K(R6B) = 2.1035945806D+13 &
	   * exp(-0.249765 * LT - 433091487.8 / RT)
      K(R7F) = 5.6200000000D+13 * exp(-1.25 * LT)
      K(R7B) = 1.3135779492D+16 &
	   * exp(-0.899765 * LT - 433091487.8 / RT)
      K(R8F) = 5.5000000000D+14 * exp(-2 * LT)
      K(R8B) = 1.2855300215D+17 &
	   * exp(-1.64977 * LT - 433091487.8 / RT)
      K(R9F) = 4.4000000000D+16 * exp(-2 * LT)
      K(R9B) = 1.1314226588D+20 &
	   * exp(-1.75351 * LT - 496136012.5 / RT)
      K(R10F) = 9.4300000000D+12 * exp(-1 * LT)
      K(R10B) = 1.3213721422D+15 &
	   * exp(-0.685917 * LT - 427002070.2 / RT)
      K(R11F) = 1.2000000000D+11 * exp(-1 * LT)
      K(R11B) = 8.4127616577D+15 &
	   * exp(-1.12163 * LT - 497822466 / RT)
      K0TROE = 6.3300000000D+13 * exp(-1.4 * LT)
      KINFTROE = 5.1200000000D+09 * exp(0.44 * LT)
      FCTROE = 0.5 * EXP( -TEMP / 1e-10 ) &
	   + 0.5 * EXP( -TEMP / 1e+10 )
      K(R12F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K0TROE = 9.1184255084D+16 &
	   * exp(-1.41184 * LT - 204667815.3 / RT)
      KINFTROE = 7.3754089420D+12 &
	   * exp(0.428161 * LT - 204667815.3 / RT)
      K(R12B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K(R13F) = 5.9200000000D+02 &
	   * exp(2.43 * LT - 223850000 / RT)
      K(R13B) = 3.6485333254D+03 &
	   * exp(2.06793 * LT + 4573672.487 / RT)
      K0TROE = 2.0100000000D+11 &
	   * exp(-0.58 * LT + 9590000 / RT)
      KINFTROE = 1.1100000000D+11 * exp(-0.37 * LT)
      FCTROE = 0.2654 * EXP( -TEMP / 94 ) &
	   + 0.7346 * EXP( -TEMP / 1756 ) &
	   + 1 * EXP( -5182 / TEMP )
      K(R14F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.1253840510D+19 &
	   * exp(-1.59836 * LT - 205148485.1 / RT)
      KINFTROE = 3.9349135804D+19 &
	   * exp(-1.38836 * LT - 214738485.1 / RT)
      K(R14B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(R15F) = 3.9700000000D+09 * exp(-2810000 / RT)
      K(R15B) = 1.4164546439D+07 &
	   * exp(0.694036 * LT - 223457801.5 / RT)
      K(R16F) = 7.4900000000D+10 * exp(-2660000 / RT)
      K(R16B) = 1.4562476643D+07 &
	   * exp(0.761631 * LT - 154173859.2 / RT)
      K(R17F) = 4.0000000000D+10
      K(R17B) = 3.8909647883D+09 &
	   * exp(0.325921 * LT - 222334254.9 / RT)
      K(R18F) = 2.3800000000D+10 * exp(2090000 / RT)
      K(R18B) = 4.2484744234D+10 &
	   * exp(0.258327 * LT - 289378197.2 / RT)
      K(R19F) = 1.0000000000D+13 * exp(-72510000 / RT)
      K(R19B) = 1.7850732872D+13 &
	   * exp(0.258327 * LT - 363978197.2 / RT)
      K(R20F) = 1.3000000000D+08 * exp(6820000 / RT)
      K(R20B) = 6.2200353343D+09 &
	   * exp(-0.244887 * LT - 154764529 / RT)
      K(R21F) = 3.6600000000D+11 * exp(-50210000 / RT)
      K(R21B) = 1.7511791788D+13 &
	   * exp(-0.244887 * LT - 211794529 / RT)
      K(R22F) = 6.0500000000D+03 &
	   * exp(2 * LT - 21760000 / RT)
      K(R22B) = 2.0516782971D+01 &
	   * exp(2.60696 * LT - 88599143.48 / RT)
      K(R23F) = 2.4100000000D+10 * exp(-16610000 / RT)
      K(R23B) = 1.7481432523D+05 &
	   * exp(1.26484 * LT - 298007527.4 / RT)
      K(R24F) = 9.6300000000D+03 &
	   * exp(2 * LT - 16610000 / RT)
      K(R24B) = 1.9578260237D+01 &
	   * exp(2.57081 * LT - 77359725.92 / RT)
      K(R25F) = 2.0000000000D+09 * exp(-1790000 / RT)
      K(R25B) = 7.4616787480D+07 &
	   * exp(0.503214 * LT - 131673668.2 / RT)
      K(R26F) = 2.6700000000D+38 &
	   * exp(-7 * LT - 157320000 / RT)
      K(R26B) = 9.9613411288D+36 &
	   * exp(-6.49679 * LT - 287203668.2 / RT)
      K0TROE = 1.1700000000D+18 &
	   * exp(-2.79 * LT - 17540000 / RT)
      KINFTROE = 1.3600000000D+07 * exp(-9980000 / RT)
      FCTROE = 1 * EXP( -0 / TEMP )
      K(R27F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K0TROE = 8.9757498242D+26 &
	   * exp(-3.75335 * LT - 552243890 / RT)
      KINFTROE = 1.0433350223D+16 &
	   * exp(-0.963346 * LT - 544683890 / RT)
      K(R27B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K(R28F) = 8.0000000000D+08 &
	   * exp(0.14 * LT - 30760000 / RT)
      K(R28B) = 4.3798719198D+15 &
	   * exp(-1.13743 * LT - 138461819.8 / RT)
      K(R29F) = 8.7800000000D+07 &
	   * exp(0.03 * LT + 70000 / RT)
      K(R29B) = 4.8069094320D+14 &
	   * exp(-1.24743 * LT - 107631819.8 / RT)
      K(R30F) = 1.1200000000D+09 * exp(-199580000 / RT)
      K(R30B) = 1.2255909945D+13 &
	   * exp(-0.84172 * LT - 236461424.1 / RT)
      K(R31F) = 3.0100000000D+10 * exp(-96230000 / RT)
      K(R31B) = 3.2039914124D+13 &
	   * exp(-0.515799 * LT - 355445679 / RT)
      K(R32F) = 1.2000000000D+11
      K(R32B) = 5.2955867164D+09 &
	   * exp(0.582947 * LT - 367317915.9 / RT)
      K(R33F) = 3.0000000000D+10
      K(R33B) = 7.9368470794D+08 &
	   * exp(0.546795 * LT - 361228498.4 / RT)
      K(R34F) = 3.0000000000D+10
      K(R34B) = 4.3452967068D+15 &
	   * exp(-0.730634 * LT - 468930318.2 / RT)
      K(R35F) = 3.0200000000D+10
      K(R35B) = 1.4661973921D+10 &
	   * exp(0.4792 * LT - 430362440.7 / RT)
      K(R36F) = 1.8700000000D+14 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R36B) = 3.5306519794D+10 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(R37F) = 2.2400000000D+15 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R37B) = 4.2292301785D+11 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(R38F) = 1.2000000000D+07 &
	   * exp(0.81 * LT + 3040000 / RT)
      K(R38B) = 3.2637034817D+06 &
	   * exp(1.03087 * LT - 135854243.4 / RT)
      K(RG01F) = 5.0000000000D+10
      K(RG01B) = 1.1487472704D+13 &
	   * exp(-0.253456 * LT - 648064893.1 / RT)
      K(RG02F) = 5.8000000000D+10 * exp(-2410000 / RT)
      K(RG02B) = 2.6634134985D+10 &
	   * exp(0.182253 * LT - 579654497.4 / RT)
      K(RG03F) = 1.6500000000D+11
      K(RG03B) = 8.4149682406D+10 &
	   * exp(0.23442 * LT - 97396261.59 / RT)
      K(RG04F) = 5.7000000000D+10
      K(RG04B) = 4.0039790897D+12 &
	   * exp(-0.0551873 * LT - 739371737.2 / RT)
      K(RG05F) = 3.0000000000D+10
      K(RG05B) = 7.9654706507D+13 &
	   * exp(-0.601982 * LT - 378143238.8 / RT)
      K(RG06F) = 1.0800000000D+11 * exp(-13010000 / RT)
      K(RG06B) = 1.7581777781D+12 &
	   * exp(-0.292468 * LT - 1691026.963 / RT)
      K0TROE = 4.8200000000D+19 &
	   * exp(-2.8 * LT - 2470000 / RT)
      KINFTROE = 1.9700000000D+09 &
	   * exp(0.43 * LT + 1550000 / RT)
      FCTROE = 0.422 * EXP( -TEMP / 122 ) &
	   + 0.578 * EXP( -TEMP / 2535 ) &
	   + 1 * EXP( -9365 / TEMP )
      K(RG07F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.6936076132D+25 &
	   * exp(-3.25474 * LT - 455400199.9 / RT)
      KINFTROE = 3.1444827797D+15 &
	   * exp(-0.0247361 * LT - 451380199.9 / RT)
      K(RG07B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG08F) = 5.7100000000D+09 * exp(3160000 / RT)
      K(RG08B) = 8.2750138742D+14 &
	   * exp(-0.947689 * LT - 248452632.8 / RT)
      K(RG09F) = 6.7100000000D+10
      K(RG09B) = 3.5609741582D+11 &
	   * exp(-0.166273 * LT - 307322843.1 / RT)
      K0TROE = 2.6900000000D+22 &
	   * exp(-3.74 * LT - 8100000 / RT)
      KINFTROE = 5.0000000000D+10
      FCTROE = 0.4243 * EXP( -TEMP / 237 ) &
	   + 0.5757 * EXP( -TEMP / 1652 ) &
	   + 1 * EXP( -5069 / TEMP )
      K(RG10F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.2130524545D+32 &
	   * exp(-5.18097 * LT - 320775841.3 / RT)
      KINFTROE = 2.2547443392D+20 &
	   * exp(-1.44097 * LT - 312675841.3 / RT)
      K(RG10B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG11F) = 1.9000000000D+11 * exp(-66070000 / RT)
      K(RG11B) = 9.2145125176D+07 &
	   * exp(0.675447 * LT - 336511419 / RT)
      K0TROE = 5.0700000000D+21 &
	   * exp(-3.42 * LT - 352920000 / RT)
      KINFTROE = 4.3000000000D+04 &
	   * exp(1.5 * LT - 333050000 / RT)
      FCTROE = 0.068 * EXP( -TEMP / 197 ) &
	   + 0.932 * EXP( -TEMP / 1540 ) &
	   + 1 * EXP( -10300 / TEMP )
      K(RG12F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.6124629783D+28 &
	   * exp(-4.10217 * LT - 355207490.6 / RT)
      KINFTROE = 1.3675721512D+11 &
	   * exp(0.817835 * LT - 335337490.6 / RT)
      K(RG12B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.4700000000D+18 &
	   * exp(-2.57 * LT - 1780000 / RT)
      KINFTROE = 1.0900000000D+09 &
	   * exp(0.48 * LT + 1090000 / RT)
      FCTROE = 0.2176 * EXP( -TEMP / 271 ) &
	   + 0.7824 * EXP( -TEMP / 2755 ) &
	   + 1 * EXP( -6570 / TEMP )
      K(RG13F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 3.4666626717D+23 &
	   * exp(-2.66922 * LT - 371385406.5 / RT)
      KINFTROE = 1.5298227985D+14 &
	   * exp(0.380781 * LT - 368515406.5 / RT)
      K(RG13B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.0400000000D+20 &
	   * exp(-2.76 * LT - 6690000 / RT)
      KINFTROE = 6.0000000000D+11
      FCTROE = 0.438 * EXP( -TEMP / 91 ) &
	   + 0.562 * EXP( -TEMP / 5836 ) &
	   + 1 * EXP( -8552 / TEMP )
      K(RG14F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.0197114661D+25 &
	   * exp(-2.92227 * LT - 470939173 / RT)
      KINFTROE = 5.8829507657D+16 &
	   * exp(-0.162268 * LT - 464249173 / RT)
      K(RG14B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG15F) = 8.0000000000D+10
      K(RG15B) = 7.8223107159D+12 &
	   * exp(-0.345666 * LT - 383372794.3 / RT)
      K(RG16F) = 2.0000000000D+10
      K(RG16B) = 1.9587375479D+15 &
	   * exp(-0.758967 * LT - 325976130.6 / RT)
      K(RG17F) = 1.1300000000D+04 &
	   * exp(2 * LT - 12550000 / RT)
      K(RG17B) = 7.6364598733D+03 &
	   * exp(2.18872 * LT - 86913497.78 / RT)
      K(RG18F) = 5.0000000000D+02 &
	   * exp(2 * LT - 30250000 / RT)
      K(RG18B) = 2.0974636059D+05 &
	   * exp(1.4875 * LT - 61407685.18 / RT)
      K(RG19) = 5.8000000000D+09 * exp(-6280000 / RT)
      K(RG20F) = 2.4000000000D+09 * exp(-6280000 / RT)
      K(RG20B) = 4.6980064580D+11 &
	   * exp(-0.323258 * LT - 261435734.9 / RT)
      K(RG21) = 5.0000000000D+09 * exp(-6280000 / RT)
      K(RG22F) = 2.0000000000D+10
      K(RG22B) = 3.8082870216D+11 &
	   * exp(0.0026632 * LT - 477489989.8 / RT)
      K(RG23F) = 5.0000000000D+10
      K(RG23B) = 9.0545188563D+14 &
	   * exp(-1.00978 * LT - 327673077.5 / RT)
      K0TROE = 2.6900000000D+27 &
	   * exp(-5.11 * LT - 29690000 / RT)
      KINFTROE = 8.1000000000D+08 &
	   * exp(0.5 * LT - 18870000 / RT)
      FCTROE = 0.4093 * EXP( -TEMP / 275 ) &
	   + 0.5907 * EXP( -TEMP / 1226 ) &
	   + 1 * EXP( -5185 / TEMP )
      K(RG24F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 9.7210519966D+38 &
	   * exp(-6.85567 * LT - 364655095.4 / RT)
      KINFTROE = 2.9271569209D+20 &
	   * exp(-1.24567 * LT - 353835095.4 / RT)
      K(RG24B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG25F) = 4.0000000000D+10
      K(RG25B) = 1.0797267305D+18 &
	   * exp(-1.3719 * LT - 548460717.8 / RT)
      K(RG26F) = 1.6000000000D+12 * exp(-49970000 / RT)
      K(RG26B) = 2.6529851155D+18 &
	   * exp(-1.07943 * LT - 609749690.8 / RT)
      K(RG27) = 2.0000000000D+11 * exp(-45980000 / RT)
      K(RG28F) = 1.5000000000D+10 * exp(-2510000 / RT)
      K(RG28B) = 1.0936408744D+10 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(RG29F) = 9.0000000000D+09 * exp(-2510000 / RT)
      K(RG29B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(RG30F) = 3.0000000000D+10
      K(RG30B) = 1.3435867056D+09 &
	   * exp(0.230182 * LT - 48621961.01 / RT)
      K(RG31F) = 1.5000000000D+10
      K(RG31B) = 4.7190290127D+10 &
	   * exp(0.174995 * LT - 787993698.2 / RT)
      K(RG32F) = 1.5000000000D+10
      K(RG32B) = 1.0693498414D+12 &
	   * exp(-0.407952 * LT - 420675782.2 / RT)
      K(RG33F) = 3.0000000000D+10
      K(RG33B) = 2.1421554446D+15 &
	   * exp(-0.821254 * LT - 363279118.6 / RT)
      K(RG34F) = 7.0000000000D+10
      K(RG34B) = 2.1409471364D+13 &
	   * exp(-0.574789 * LT - 68460673.16 / RT)
      K(RG35F) = 2.8000000000D+10
      K(RG35B) = 7.5327924346D+05 &
	   * exp(0.260469 * LT - 284081814.6 / RT)
      K(RG36F) = 1.2000000000D+10
      K(RG36B) = 8.3014013421D+08 &
	   * exp(0.506957 * LT - 780217827.2 / RT)
      K0TROE = 1.8800000000D+32 &
	   * exp(-6.36 * LT - 21090000 / RT)
      KINFTROE = 4.8200000000D+14 &
	   * exp(-1.16 * LT - 4790000 / RT)
      FCTROE = 0.3973 * EXP( -TEMP / 208 ) &
	   + 0.6027 * EXP( -TEMP / 3922 ) &
	   + 1 * EXP( -10180 / TEMP )
      K(RG37F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.6292506138D+42 &
	   * exp(-7.66989 * LT - 416708647.8 / RT)
      KINFTROE = 4.1771212546D+24 &
	   * exp(-2.46989 * LT - 400408647.8 / RT)
      K(RG37B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG38F) = 3.0000000000D+10
      K(RG38B) = 2.1872817488D+10 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG39) = 6.8200000000D+07 &
	   * exp(0.25 * LT + 3910000 / RT)
      K(RG40F) = 9.0000000000D+09
      K(RG40B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG41F) = 7.0000000000D+09
      K(RG41B) = 5.1036574138D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG42F) = 1.4000000000D+10
      K(RG42B) = 1.8259393135D+08 &
	   * exp(0.456176 * LT - 255577298.7 / RT)
      K0TROE = 1.2700000000D+26 &
	   * exp(-4.82 * LT - 27320000 / RT)
      KINFTROE = 5.4000000000D+08 &
	   * exp(0.45 * LT - 15060000 / RT)
      FCTROE = 0.2813 * EXP( -TEMP / 103 ) &
	   + 0.7187 * EXP( -TEMP / 1291 ) &
	   + 1 * EXP( -4160 / TEMP )
      K(RG43F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.2684133243D+30 &
	   * exp(-5.20137 * LT - 150083586.5 / RT)
      KINFTROE = 5.3932535051D+12 &
	   * exp(0.0686333 * LT - 137823586.5 / RT)
      K(RG43B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.2000000000D+24 &
	   * exp(-4.8 * LT - 23260000 / RT)
      KINFTROE = 5.4000000000D+08 &
	   * exp(0.45 * LT - 10880000 / RT)
      FCTROE = 0.242 * EXP( -TEMP / 94 ) &
	   + 0.758 * EXP( -TEMP / 1555 ) &
	   + 1 * EXP( -4200 / TEMP )
      K(RG44F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.9434371792D+29 &
	   * exp(-5.12229 * LT - 117496132.5 / RT)
      KINFTROE = 4.7702548944D+13 &
	   * exp(0.127708 * LT - 105116132.5 / RT)
      K(RG44B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG45F) = 5.7400000000D+04 &
	   * exp(1.9 * LT - 11470000 / RT)
      K(RG45B) = 9.5590953891D+01 &
	   * exp(2.34945 * LT - 74956081.25 / RT)
      K(RG46F) = 3.9000000000D+10 * exp(-14810000 / RT)
      K(RG46B) = 3.8937084533D+07 &
	   * exp(0.413302 * LT - 72206663.69 / RT)
      K(RG47F) = 3.4300000000D+06 &
	   * exp(1.18 * LT + 1870000 / RT)
      K(RG47B) = 6.2842244016D+04 &
	   * exp(1.52571 * LT - 124660606 / RT)
      K(RG48F) = 1.0000000000D+11 * exp(-167360000 / RT)
      K(RG48B) = 1.0263642436D+09 &
	   * exp(0.0873801 * LT - 2422408.767 / RT)
      K(RG49F) = 5.6000000000D+03 &
	   * exp(2 * LT - 50210000 / RT)
      K(RG49B) = 2.7500401862D+03 &
	   * exp(1.84249 * LT - 46856937.77 / RT)
      K(RG50F) = 9.4600000000D+10 * exp(2160000 / RT)
      K(RG50B) = 1.7498844055D+17 &
	   * exp(-1.35597 * LT - 319198631.8 / RT)
      K0TROE = 3.4700000000D+32 &
	   * exp(-6.3 * LT - 21230000 / RT)
      KINFTROE = 6.9200000000D+10 * exp(0.18 * LT)
      FCTROE = 0.217 * EXP( -TEMP / 74 ) &
	   + 0.783 * EXP( -TEMP / 2941 ) &
	   + 1 * EXP( -6964 / TEMP )
      K(RG51F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K0TROE = 1.5077439481D+38 &
	   * exp(-6.46997 * LT - 463233395.7 / RT)
      KINFTROE = 3.0067977293D+16 &
	   * exp(0.0100263 * LT - 442003395.7 / RT)
      K(RG51B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K(RG52F) = 5.0600000000D+10
      K(RG52B) = 7.0821685788D+12 &
	   * exp(-0.282617 * LT - 288729027.8 / RT)
      K(RG53) = 3.3700000000D+10
      K0TROE = 4.0000000000D+30 &
	   * exp(-5.92 * LT - 13140000 / RT)
      KINFTROE = 2.7900000000D+15 &
	   * exp(-1.43 * LT - 5570000 / RT)
      FCTROE = 0.588 * EXP( -TEMP / 195 ) &
	   + 0.412 * EXP( -TEMP / 5900 ) &
	   + 1 * EXP( -6394 / TEMP )
      K(RG54F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.2469090961D+39 &
	   * exp(-6.75885 * LT - 403342499.3 / RT)
      KINFTROE = 8.6971909454D+23 &
	   * exp(-2.26885 * LT - 395772499.3 / RT)
      K(RG54B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG55F) = 5.6000000000D+04 &
	   * exp(1.6 * LT - 22680000 / RT)
      K(RG55B) = 1.4686430659D+03 &
	   * exp(2.00876 * LT - 54566839.56 / RT)
      K(RG56) = 8.0000000000D+06 * exp(7340000 / RT)
      K(RG57F) = 6.4400000000D+14 &
	   * exp(-1.34 * LT - 5930000 / RT)
      K(RG57B) = 2.3164910420D+13 &
	   * exp(-0.868957 * LT - 513851.5818 / RT)
      K(RG58F) = 1.3800000000D+10 * exp(-127700000 / RT)
      K(RG58B) = 2.4338024093D+12 &
	   * exp(-0.483283 * LT - 12842694.35 / RT)
      K(RG59F) = 5.8700000000D+08 * exp(-57910000 / RT)
      K(RG59B) = 1.6421391660D+08 &
	   * exp(0.153092 * LT - 275818632.1 / RT)
      K0TROE = 3.8200000000D+25 &
	   * exp(-4.89 * LT - 14360000 / RT)
      KINFTROE = 1.0100000000D+05 * exp(1.63 * LT)
      FCTROE = 0.955 * EXP( -TEMP / 880.1 ) &
	   + 0.045 * EXP( -TEMP / 2.5e+09 ) &
	   + 1 * EXP( -1.786e+09 / TEMP )
      K(RG60F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM0) )
      K0TROE = 1.4075759346D+36 &
	   * exp(-6.63941 * LT - 158160807.4 / RT)
      KINFTROE = 3.7216012929D+15 &
	   * exp(-0.119414 * LT - 143800807.4 / RT)
      K(RG60B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM0) )
      K(RG61F) = 1.0000000000D+10 * exp(5020000 / RT)
      K(RG61B) = 5.9178142524D+08 &
	   * exp(0.661222 * LT - 119287047.3 / RT)
      K(RG62) = 1.4000000000D+13 &
	   * exp(-1.61 * LT - 7780000 / RT)
      K(RG63) = 2.4700000000D+08 * exp(6570000 / RT)
      K(RG64) = 1.9900000000D+09 * exp(-48830000 / RT)
      K(RG65F) = 1.0000000000D+10
      K(RG65B) = 1.7155506298D+11 &
	   * exp(-0.157361 * LT - 107476949.3 / RT)
      K(RG66F) = 3.6100000000D+09
      K(RG66B) = 1.0889029554D+12 &
	   * exp(-0.158135 * LT - 237335580.4 / RT)
      K(RG67F) = 2.4500000000D+01 &
	   * exp(2.47 * LT - 21670000 / RT)
      K(RG67B) = 1.5445378643D+02 &
	   * exp(2.55675 * LT - 97421051.39 / RT)
      K(RG68F) = 5.0000000000D+10
      K(RG68B) = 6.3085696118D+15 &
	   * exp(-1.09382 * LT - 419906771 / RT)
      K(RG69F) = 3.0000000000D+10
      K(RG69B) = 1.8300922168D+15 &
	   * exp(-1.03775 * LT - 231012782.6 / RT)
      K(RG70F) = 2.6500000000D+10
      K(RG70B) = 2.1739887746D+12 &
	   * exp(0.0627383 * LT - 376229823.8 / RT)
      K(RG71F) = 3.3200000000D+00 &
	   * exp(2.81 * LT - 24520000 / RT)
      K(RG71B) = 1.0278306684D+01 &
	   * exp(2.73924 * LT - 96917989.16 / RT)
      K(RG72F) = 1.0000000000D+11
      K(RG72B) = 8.9109851476D+17 &
	   * exp(-1.23115 * LT - 275759524.7 / RT)
      K(RG73F) = 1.2000000000D+10 * exp(2390000 / RT)
      K(RG73B) = 7.7963340707D+16 &
	   * exp(-1.29344 * LT - 310672512.6 / RT)
      K(RG74F) = 6.8400000000D+09 &
	   * exp(0.1 * LT - 44350000 / RT)
      K(RG74B) = 1.1614797579D+15 &
	   * exp(-1.0902 * LT - 7197984.006 / RT)
      K0TROE = 4.6600000000D+35 &
	   * exp(-7.44 * LT - 58910000 / RT)
      KINFTROE = 2.4300000000D+09 &
	   * exp(0.52 * LT - 210000 / RT)
      FCTROE = 0.3 * EXP( -TEMP / 100 ) &
	   + 0.7 * EXP( -TEMP / 90000 ) &
	   + 1 * EXP( -10000 / TEMP )
      K(RG75F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.6463061592D+39 &
	   * exp(-7.35986 * LT - 493149409.3 / RT)
      KINFTROE = 8.5848153795D+12 &
	   * exp(0.600144 * LT - 434449409.3 / RT)
      K(RG75B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG76F) = 4.1500000000D+04 &
	   * exp(1.63 * LT - 8050000 / RT)
      K(RG76B) = 4.6919928896D+03 &
	   * exp(1.57093 * LT - 36577454.01 / RT)
      K(RG77F) = 2.0000000000D+10
      K(RG77B) = 5.2917790613D+07 &
	   * exp(0.672527 * LT - 338855355.3 / RT)
      K(RG78F) = 1.5000000000D+09 &
	   * exp(0.5 * LT + 460000 / RT)
      K(RG78B) = 1.6999695549D+04 &
	   * exp(1.41899 * LT - 43576909.93 / RT)
      K(RG79F) = 2.6200000000D+11 &
	   * exp(-0.23 * LT - 4480000 / RT)
      K(RG79B) = 1.0680606960D+05 &
	   * exp(1.16003 * LT - 43100761.51 / RT)
      K(RG80F) = 1.0000000000D+10
      K(RG80B) = 1.5862280588D+07 &
	   * exp(0.636375 * LT - 332765937.8 / RT)
      K(RG81F) = 5.0000000000D+09
      K(RG81B) = 1.4554402259D+08 &
	   * exp(0.568781 * LT - 401899880.1 / RT)
      K(RG82F) = 4.2800000000D-16 &
	   * exp(7.6 * LT + 14770000 / RT)
      K(RG82B) = 6.9793035518D-18 &
	   * exp(7.91045 * LT - 95661682.83 / RT)
      K0TROE = 4.3600000000D+25 &
	   * exp(-4.65 * LT - 21260000 / RT)
      KINFTROE = 1.0600000000D+09 &
	   * exp(0.5 * LT - 360000 / RT)
      FCTROE = 0.4 * EXP( -TEMP / 100 ) &
	   + 0.6 * EXP( -TEMP / 9000 ) &
	   + 1 * EXP( -10000 / TEMP )
      K(RG83F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.3623915050D+30 &
	   * exp(-4.51078 * LT - 426971955.3 / RT)
      KINFTROE = 3.3122362277D+13 &
	   * exp(0.639219 * LT - 406071955.3 / RT)
      K(RG83B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG84F) = 2.0000000000D+10
      K(RG84B) = 4.6805022133D+08 &
	   * exp(0.731602 * LT - 310327901.3 / RT)
      K(RG85F) = 1.6500000000D+08 &
	   * exp(0.65 * LT + 1190000 / RT)
      K(RG85B) = 1.6539583927D+04 &
	   * exp(1.62807 * LT - 14319455.91 / RT)
      K(RG86F) = 3.2800000000D+10 &
	   * exp(-0.09 * LT - 2550000 / RT)
      K(RG86B) = 1.1826581729D+05 &
	   * exp(1.35911 * LT - 12643307.5 / RT)
      K(RG87F) = 1.0000000000D+10
      K(RG87B) = 1.4029958269D+08 &
	   * exp(0.69545 * LT - 304238483.7 / RT)
      K(RG88F) = 5.0000000000D+09
      K(RG88B) = 1.2873158761D+09 &
	   * exp(0.627855 * LT - 373372426.1 / RT)
      K(RG89F) = 1.8000000000D+10 * exp(-3770000 / RT)
      K(RG89B) = 2.5961607220D+09 &
	   * exp(0.369528 * LT - 85674228.82 / RT)
      K(RG90F) = 6.6000000000D+05 &
	   * exp(1.62 * LT - 45360000 / RT)
      K(RG90B) = 3.5503024345D+02 &
	   * exp(2.14021 * LT - 36448092.09 / RT)
      K(RG91F) = 1.0200000000D+06 &
	   * exp(1.5 * LT - 35980000 / RT)
      K(RG91B) = 3.2893910502D+02 &
	   * exp(1.98406 * LT - 20978674.53 / RT)
      K(RG92F) = 1.0000000000D+05 &
	   * exp(1.6 * LT - 13050000 / RT)
      K(RG92B) = 5.9179879483D+02 &
	   * exp(2.01646 * LT - 67182616.83 / RT)
      K(RG93F) = 6.0000000000D+10
      K(RG93B) = 4.6820647566D+15 &
	   * exp(-1.00341 * LT - 255528643.7 / RT)
      K(RG94F) = 2.4600000000D+03 &
	   * exp(2 * LT - 34600000 / RT)
      K(RG94B) = 5.5511242907D+02 &
	   * exp(2.00771 * LT - 56845777.28 / RT)
      K(RG95F) = 1.6000000000D+10 * exp(2390000 / RT)
      K(RG95B) = 2.6323843569D+09 &
	   * exp(-0.0545804 * LT - 57158765.25 / RT)
      K(RG96F) = 1.7000000000D+04 &
	   * exp(2.1 * LT - 20380000 / RT)
      K(RG96B) = 1.2716073149D+02 &
	   * exp(2.31102 * LT - 47759532.54 / RT)
      K(RG97F) = 4.2000000000D+03 &
	   * exp(2.1 * LT - 20380000 / RT)
      K(RG97B) = 2.7787158476D+02 &
	   * exp(2.37009 * LT - 19232078.53 / RT)
      K(RG98F) = 3.8800000000D+02 &
	   * exp(2.5 * LT - 12970000 / RT)
      K(RG98B) = 1.7399218488D+00 &
	   * exp(2.67486 * LT - 34260114.98 / RT)
      K(RG99F) = 1.3000000000D+02 &
	   * exp(2.5 * LT - 20920000 / RT)
      K(RG99B) = 5.1562280825D+00 &
	   * exp(2.73394 * LT - 13682660.97 / RT)
      K(RG100F) = 1.4400000000D+03 &
	   * exp(2 * LT + 3520000 / RT)
      K(RG100B) = 1.1850024546D+02 &
	   * exp(2.10727 * LT - 86904057.28 / RT)
      K(RG101F) = 6.3000000000D+03 &
	   * exp(2 * LT - 6280000 / RT)
      K(RG101B) = 4.5855143694D+03 &
	   * exp(2.16634 * LT - 68176603.27 / RT)
      K(RG102F) = 3.0000000000D+04 &
	   * exp(1.5 * LT - 41590000 / RT)
      K(RG102B) = 4.1716122698D+05 &
	   * exp(1.19081 * LT - 77881440.45 / RT)
      K(RG103F) = 1.0000000000D+04 &
	   * exp(1.5 * LT - 41590000 / RT)
      K(RG103B) = 1.2299102838D+06 &
	   * exp(1.24988 * LT - 49353986.44 / RT)
      K0TROE = 2.6000000000D+27 &
	   * exp(-4.8 * LT - 7950000 / RT)
      KINFTROE = 1.0000000000D+14 * exp(-1 * LT)
      FCTROE = 0.3536 * EXP( -TEMP / 132 ) &
	   + 0.6464 * EXP( -TEMP / 1315 ) &
	   + 1 * EXP( -5566 / TEMP )
      K(RG104F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.7761630600D+33 &
	   * exp(-5.0463 * LT - 564432866.5 / RT)
      KINFTROE = 6.8313963848D+19 &
	   * exp(-1.2463 * LT - 556482866.5 / RT)
      K(RG104B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG105F) = 5.0000000000D+10
      K(RG105B) = 2.3360600807D+07 &
	   * exp(1.01264 * LT - 325621371.1 / RT)
      K(RG106F) = 2.0000000000D+10
      K(RG106B) = 3.0071690750D+14 &
	   * exp(-0.742405 * LT - 211295142.2 / RT)
      K(RG107F) = 1.0000000000D+10 * exp(3160000 / RT)
      K(RG107B) = 2.4794782651D+07 &
	   * exp(0.846372 * LT - 629784214.1 / RT)
      K(RG108F) = 3.3100000000D+03 &
	   * exp(2.26 * LT - 3770000 / RT)
      K(RG108B) = 9.6742642418D+06 &
	   * exp(1.66346 * LT - 127161378.7 / RT)
      K(RG109F) = 1.0000000000D+11
      K(RG109B) = 1.1573048940D+05 &
	   * exp(1.56102 * LT - 71793685.47 / RT)
      K(RG110F) = 1.0000000000D+11
      K(RG110B) = 1.5577208903D+03 &
	   * exp(1.38578 * LT - 426695895.8 / RT)
      K(RG111F) = 4.2000000000D+07 * exp(-3570000 / RT)
      K(RG111B) = 1.3076606325D-03 &
	   * exp(1.82149 * LT - 359445500.1 / RT)
      K(RG112F) = 5.0000000000D+10
      K(RG112B) = 1.1388201171D+12 &
	   * exp(0.12683 * LT - 657557391.2 / RT)
      K(RG113F) = 3.0000000000D+10
      K(RG113B) = 3.9791438217D+10 &
	   * exp(0.240947 * LT - 382586114.2 / RT)
      K(RG114F) = 1.0000000000D+10
      K(RG114B) = 5.0507727073D+01 &
	   * exp(1.5678 * LT - 344881549.9 / RT)
      K0TROE = 6.3400000000D+25 &
	   * exp(-4.66 * LT - 15820000 / RT)
      KINFTROE = 1.7100000000D+07 &
	   * exp(1.27 * LT - 11330000 / RT)
      FCTROE = 0.2122 * EXP( -TEMP / -10212 )
      K(RG115F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.4048536324D+28 &
	   * exp(-4.48812 * LT - 162621237.8 / RT)
      KINFTROE = 3.7891162641D+09 &
	   * exp(1.44188 * LT - 158131237.8 / RT)
      K(RG115B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG116F) = 8.1000000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RG116B) = 2.4981403723D+04 &
	   * exp(1.81798 * LT - 89764345.93 / RT)
      K(RG117F) = 1.2500000000D+04 &
	   * exp(2 * LT - 7950000 / RT)
      K(RG117B) = 3.2529178544D-02 &
	   * exp(3.31672 * LT - 198861019.4 / RT)
      K(RG118F) = 1.8100000000D+10
      K(RG118B) = 8.8241850413D+13 &
	   * exp(-0.560387 * LT - 129480796.2 / RT)
      K(RG119F) = 2.6300000000D+03 &
	   * exp(2.14 * LT - 71380000 / RT)
      K(RG119B) = 9.8996184977D+00 &
	   * exp(2.63279 * LT - 11033146.07 / RT)
      K(RG120F) = 2.4100000000D+03 &
	   * exp(2 * LT - 53190000 / RT)
      K(RG120B) = 4.3482685714D+08 &
	   * exp(0.845959 * LT - 26274495.32 / RT)
      K(RG121F) = 7.5300000000D+03 &
	   * exp(1.55 * LT - 8810000 / RT)
      K(RG121B) = 5.0536550148D+07 &
	   * exp(0.806965 * LT - 107684044.6 / RT)
      K(RG122F) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RG122B) = 2.3307922343D+03 &
	   * exp(1.57037 * LT - 238948122.1 / RT)
      K(RG123F) = 5.0000000000D+10 * exp(-33470000 / RT)
      K(RG123B) = 3.8326334707D+07 &
	   * exp(0.59717 * LT - 22499718.93 / RT)
      K(RG124F) = 1.5000000000D+06 &
	   * exp(1.38 * LT - 2570000 / RT)
      K(RG124B) = 4.0698109539D-01 &
	   * exp(2.9634 * LT - 131854077.6 / RT)
      K(RG125F) = 1.0000000000D+10 * exp(-33470000 / RT)
      K(RG125B) = 4.5953776074D+06 &
	   * exp(0.561018 * LT - 16410301.37 / RT)
      K(RG126F) = 1.7500000000D+09 * exp(-5650000 / RT)
      K(RG126B) = 3.7150267058D+06 &
	   * exp(0.782321 * LT - 205388794.6 / RT)
      K(RG127F) = 7.5000000000D+09 * exp(-8370000 / RT)
      K(RG127B) = 6.3247186168D+07 &
	   * exp(0.493423 * LT - 60444243.67 / RT)
      K(RG128F) = 1.0000000000D+10
      K(RG128B) = 3.7197289740D+08 &
	   * exp(0.411006 * LT - 125789549.2 / RT)
      K0TROE = 1.4000000000D+24 &
	   * exp(-3.86 * LT - 13890000 / RT)
      KINFTROE = 6.0800000000D+09 &
	   * exp(0.27 * LT - 1170000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(RG129F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.7814512532D+29 &
	   * exp(-3.99563 * LT - 480409256.8 / RT)
      KINFTROE = 3.3793731157D+15 &
	   * exp(0.134367 * LT - 467689256.8 / RT)
      K(RG129B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG130F) = 3.0000000000D+10
      K(RG130B) = 3.1644553172D+10 &
	   * exp(0.178352 * LT - 286290250 / RT)
      K(RG131F) = 1.0300000000D+10 &
	   * exp(0.21 * LT + 1790000 / RT)
      K(RG131B) = 6.4367857895D+18 &
	   * exp(-0.90671 * LT - 531013741.4 / RT)
      K(RG132F) = 5.0000000000D+09
      K(RG132B) = 5.8023026639D+10 &
	   * exp(0.0746054 * LT - 349334774.8 / RT)
      K(RG133F) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RG133B) = 8.7112228115D+03 &
	   * exp(1.42628 * LT - 56256577.53 / RT)
      K(RG134F) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RG134B) = 2.7009535860D+12 &
	   * exp(-0.705084 * LT - 35031275.47 / RT)
      K(RG135F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RG135B) = 5.5766810183D+11 &
	   * exp(-0.801138 * LT - 369289088.3 / RT)
      K(RG136F) = 1.3200000000D+34 &
	   * exp(-6.57 * LT - 206930000 / RT)
      K(RG136B) = 8.9644644228D+28 &
	   * exp(-6.05412 * LT - 53201135.59 / RT)
      K(RG137F) = 6.5100000000D+34 &
	   * exp(-6.87 * LT - 197460000 / RT)
      K(RG137B) = 1.1995390279D+23 &
	   * exp(-4.77073 * LT - 173015213.1 / RT)
      K(RG138F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RG138B) = 4.3300717766D+04 &
	   * exp(1.61395 * LT - 328407812.9 / RT)
      K(RG139) = 1.8100000000D+07
      K(RG140) = 2.3500000000D+07
      K(RG141F) = 2.2000000000D+10
      K(RG141B) = 2.1470530275D+02 &
	   * exp(1.86656 * LT - 41328785.02 / RT)
      K(RG142F) = 1.1000000000D+10
      K(RG142B) = 1.7460739608D+07 &
	   * exp(0.86611 * LT - 279362623.4 / RT)
      K(RG143F) = 1.2000000000D+10
      K(RG143B) = 2.0955781346D+08 &
	   * exp(0.762364 * LT - 342407148.1 / RT)
      K(RG144F) = 3.0100000000D+10
      K(RG144B) = 2.9305285012D+06 &
	   * exp(0.888497 * LT - 25819329.1 / RT)
      K(RG145F) = 5.0000000000D+10
      K(RG145B) = 1.3886732031D+22 &
	   * exp(-1.57811 * LT - 365339262.7 / RT)
      K(RG146F) = 2.9200000000D+09 * exp(-7570000 / RT)
      K(RG146B) = 1.5095453499D+08 &
	   * exp(0.0256312 * LT - 27904022.48 / RT)
      K(RG147F) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RG147B) = 1.7677598030D+05 &
	   * exp(1.22178 * LT - 36483440.04 / RT)
      K(RG148) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RG149) = 2.9200000000D+09 * exp(-7570000 / RT)
      K(RG150) = 3.0100000000D+10 * exp(-163800000 / RT)
      K(RG151) = 2.3400000000D+07 &
	   * exp(0.73 * LT + 4660000 / RT)
      K(RG152) = 3.0100000000D+09 * exp(-49890000 / RT)
      K(RG153) = 2.7200000000D+03 &
	   * exp(1.77 * LT - 24770000 / RT)
      K0TROE = 7.0000000000D+47 &
	   * exp(-9.31 * LT - 417980000 / RT)
      KINFTROE = 8.0000000000D+12 &
	   * exp(0.44 * LT - 371540000 / RT)
      FCTROE = 0.265 * EXP( -TEMP / 180 ) &
	   + 0.735 * EXP( -TEMP / 1035 ) &
	   + 1 * EXP( -5417 / TEMP )
      K(RG154F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM6) )
      K0TROE = 2.8487507449D+39 &
	   * exp(-8.45483 * LT - 49499564.99 / RT)
      KINFTROE = 3.2557151370D+04 &
	   * exp(1.29517 * LT - 3059564.99 / RT)
      K(RG154B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM6) )
      K0TROE = 2.0300000000D+33 &
	   * exp(-6.64 * LT - 24140000 / RT)
      KINFTROE = 1.3700000000D+06 &
	   * exp(1.46 * LT - 5670000 / RT)
      FCTROE = 1.569 * EXP( -TEMP / 299 ) &
	   + -0.569 * EXP( -TEMP / -9147 ) &
	   + 1 * EXP( -152.4 / TEMP )
      K(RG155F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 3.7928863106D+36 &
	   * exp(-6.76131 * LT - 175477632.3 / RT)
      KINFTROE = 2.5597311555D+09 &
	   * exp(1.33869 * LT - 157007632.3 / RT)
      K(RG155B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG156F) = 1.2700000000D+02 &
	   * exp(2.75 * LT - 48740000 / RT)
      K(RG156B) = 5.3406078444D-02 &
	   * exp(3.23587 * LT - 15312230.95 / RT)
      K(RG157F) = 7.6600000000D+06 &
	   * exp(0.88 * LT - 4770000 / RT)
      K(RG157B) = 8.6124762219D+09 &
	   * exp(-0.101077 * LT - 71054484.59 / RT)
      K(RG158F) = 7.1500000000D+01 &
	   * exp(2.47 * LT - 3890000 / RT)
      K(RG158B) = 1.1230421728D-03 &
	   * exp(3.41854 * LT - 16859503.17 / RT)
      K(RG159F) = 3.8900000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RG159B) = 4.2684378019D+00 &
	   * exp(2.24549 * LT - 111323269.6 / RT)
      K(RG160F) = 1.3100000000D-04 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RG160B) = 6.0605344705D-07 &
	   * exp(4.58212 * LT - 26016755.69 / RT)
      K(RG161F) = 3.7500000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RG161B) = 1.4704165297D+43 &
	   * exp(-9.27438 * LT - 141884176.2 / RT)
      K(RG162F) = 2.2700000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RG162B) = 1.7745629538D+02 &
	   * exp(1.96566 * LT - 13974138.86 / RT)
      K0TROE = 3.0000000000D+57 &
	   * exp(-14.6 * LT - 76020000 / RT)
      KINFTROE = 2.5500000000D+03 &
	   * exp(1.6 * LT - 23850000 / RT)
      FCTROE = 0.8106 * EXP( -TEMP / 277 ) &
	   + 0.1894 * EXP( -TEMP / 8748 ) &
	   + 1 * EXP( -7891 / TEMP )
      K(RG163F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.1619171325D+67 &
	   * exp(-16.2667 * LT - 179210339.7 / RT)
      KINFTROE = 1.8376295626D+13 &
	   * exp(-0.0666836 * LT - 127040339.7 / RT)
      K(RG163B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.9900000000D+35 &
	   * exp(-7.08 * LT - 27970000 / RT)
      KINFTROE = 5.2100000000D+14 &
	   * exp(-0.99 * LT - 6610000 / RT)
      FCTROE = 0.1578 * EXP( -TEMP / 125 ) &
	   + 0.8422 * EXP( -TEMP / 2219 ) &
	   + 1 * EXP( -6882 / TEMP )
      K(RG164F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.4724095869D+41 &
	   * exp(-7.28894 * LT - 449644917.7 / RT)
      KINFTROE = 6.4729919335D+20 &
	   * exp(-1.19894 * LT - 428284917.7 / RT)
      K(RG164B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG165F) = 2.0000000000D+09
      K(RG165B) = 2.5019333849D+08 &
	   * exp(0.471548 * LT - 281753855.5 / RT)
      K(RG166F) = 1.1800000000D+01 &
	   * exp(2.45 * LT - 12220000 / RT)
      K(RG166B) = 2.7441404727D+03 &
	   * exp(2.40134 * LT - 302885763.4 / RT)
      K(RG167F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RG167B) = 9.3219772183D+18 &
	   * exp(-1.00898 * LT - 386358614.1 / RT)
      K(RG168F) = 1.3200000000D+20 &
	   * exp(-2.02 * LT - 86820000 / RT)
      K(RG168B) = 3.6998605525D+08 &
	   * exp(-0.0734395 * LT - 24692429.7 / RT)
      K(RG169F) = 5.4500000000D+15 &
	   * exp(-0.69 * LT - 93010000 / RT)
      K(RG169B) = 3.0228969809D+10 &
	   * exp(-0.222332 * LT - 26616285.89 / RT)
      K(RG170F) = 2.2900000000D+07 * exp(-3660000 / RT)
      K(RG170B) = 1.8296940316D+05 &
	   * exp(0.455829 * LT - 141934101.2 / RT)
      K(RG171F) = 1.9200000000D+04 &
	   * exp(1.02 * LT + 8510000 / RT)
      K(RG171B) = 1.4802790269D+04 &
	   * exp(1.12947 * LT - 44820182.98 / RT)
      K0TROE = 5.6400000000D+71 &
	   * exp(-15.74 * LT - 413040000 / RT)
      KINFTROE = 1.2900000000D+37 &
	   * exp(-5.84 * LT - 407470000 / RT)
      FCTROE = 0.69 * EXP( -TEMP / 50 ) &
	   + 0.31 * EXP( -TEMP / 3000 ) &
	   + 1 * EXP( -9000 / TEMP )
      K(RG172F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 8.9945859647D+58 &
	   * exp(-13.9719 * LT - 37586111.82 / RT)
      KINFTROE = 2.0572723217D+24 &
	   * exp(-4.07189 * LT - 32016111.82 / RT)
      K(RG172B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 3.7200000000D+62 &
	   * exp(-13.14 * LT - 425010000 / RT)
      KINFTROE = 1.8800000000D+50 &
	   * exp(-9.72 * LT - 449120000 / RT)
      FCTROE = 0.61 * EXP( -TEMP / 100 ) &
	   + 0.39 * EXP( -TEMP / 1900 ) &
	   + 1 * EXP( -6000 / TEMP )
      K(RG173F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.7632750149D+51 &
	   * exp(-11.7409 * LT - 40487098.26 / RT)
      KINFTROE = 8.9111748064D+38 &
	   * exp(-8.32086 * LT - 64597098.26 / RT)
      K(RG173B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG174F) = 1.7000000000D+02 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RG174B) = 3.1981668158D-02 &
	   * exp(3.25918 * LT - 35436570.06 / RT)
      K(RG175F) = 3.1700000000D-02 &
	   * exp(3.8 * LT - 13100000 / RT)
      K(RG175B) = 3.5752414238D-06 &
	   * exp(4.32303 * LT - 18427152.49 / RT)
      K(RG176F) = 1.6100000000D+03 &
	   * exp(2.22 * LT - 3100000 / RT)
      K(RG176B) = 3.3321974605D+00 &
	   * exp(2.67543 * LT - 77561094.8 / RT)
      K(RG177F) = 4.0000000000D+10 * exp(2300000 / RT)
      K(RG177B) = 2.3015482642D+09 &
	   * exp(-0.0156091 * LT - 77577243.21 / RT)
      K(RG178F) = 8.4300000000D+11 * exp(-93120000 / RT)
      K(RG178B) = 2.9482095678D+11 &
	   * exp(0.0389713 * LT - 113448478 / RT)
      K(RG179F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RG179B) = 1.1503526688D+06 &
	   * exp(1.29275 * LT - 335226320.5 / RT)
      K0TROE = 4.4200000000D+55 &
	   * exp(-13.55 * LT - 47520000 / RT)
      KINFTROE = 3.6100000000D+10
      FCTROE = 0.685 * EXP( -TEMP / 369 ) &
	   + 0.315 * EXP( -TEMP / 3285 ) &
	   + 1 * EXP( -6667 / TEMP )
      K(RG180F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.1858177017D+61 &
	   * exp(-13.7727 * LT - 471121180.8 / RT)
      KINFTROE = 5.8689597066D+16 &
	   * exp(-0.222745 * LT - 423601180.8 / RT)
      K(RG180B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG181F) = 2.4100000000D+10
      K(RG181B) = 1.5941539969D+10 &
	   * exp(0.382569 * LT - 359404694.1 / RT)
      K(RG182F) = 3.3100000000D+09 * exp(3220000 / RT)
      K(RG182B) = 3.6997054316D+11 &
	   * exp(-0.0338928 * LT - 302052077.3 / RT)
      K0TROE = 6.2600000000D+32 &
	   * exp(-6.66 * LT - 29290000 / RT)
      KINFTROE = 3.0600000000D+11 &
	   * exp(-0.37 * LT - 16870000 / RT)
      FCTROE = 1 * EXP( -TEMP / 1310 ) &
	   + 1 * EXP( -48097 / TEMP )
      K(RG183F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.4335109175D+36 &
	   * exp(-6.79608 * LT - 166021318.4 / RT)
      KINFTROE = 1.1895436753D+15 &
	   * exp(-0.506081 * LT - 153601318.4 / RT)
      K(RG183B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG184F) = 3.7000000000D+13 &
	   * exp(-1.63 * LT - 14300000 / RT)
      K(RG184B) = 1.3710677099D+13 &
	   * exp(-1.50576 * LT - 82236496.87 / RT)
      K(RG185F) = 5.8000000000D-05 &
	   * exp(4.71 * LT - 25990000 / RT)
      K(RG185B) = 8.3386088338D-09 &
	   * exp(5.28298 * LT - 35480306.99 / RT)
      K(RG186F) = 2.3500000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RG186B) = 2.0254773459D-07 &
	   * exp(4.62683 * LT - 14050889.43 / RT)
      K(RG187F) = 5.3600000000D+03 &
	   * exp(2.01 * LT - 1530000 / RT)
      K(RG187B) = 8.4777978991D+00 &
	   * exp(2.47923 * LT - 74064831.73 / RT)
      K(RG188F) = 9.0300000000D-04 &
	   * exp(3.65 * LT - 29930000 / RT)
      K(RG188B) = 2.4134146435D-04 &
	   * exp(3.70277 * LT - 48332214.89 / RT)
      K(RG189F) = 9.6400000000D+00 &
	   * exp(2.6 * LT - 58200000 / RT)
      K(RG189B) = 4.0868505075D-01 &
	   * exp(2.56602 * LT - 851163.5039 / RT)
      K(RR001F) = 2.4500000000D+12 &
	   * exp(-0.64 * LT - 207940000 / RT)
      K(RR001B) = 5.2538461199D+09 &
	   * exp(-0.0988122 * LT - 19688571.81 / RT)
      K(RR002F) = 3.3000000000D+09
      K(RR002B) = 4.3822111283D+10 &
	   * exp(-0.0661922 * LT - 416043871.6 / RT)
      K(RR003F) = 1.0000000000D+10
      K(RR003B) = 8.9645037186D+07 &
	   * exp(0.318776 * LT - 330486347.7 / RT)
      K(RR004F) = 1.9000000000D+11
      K(RR004B) = 2.5098929533D+16 &
	   * exp(-1.35613 * LT - 105055737.5 / RT)
      K(RR005F) = 3.4600000000D+09 &
	   * exp(0.44 * LT - 22860000 / RT)
      K(RR005B) = 1.8216946136D+04 &
	   * exp(1.58056 * LT - 46357562.82 / RT)
      K(RR006F) = 8.9500000000D+10 &
	   * exp(-0.02 * LT - 47070000 / RT)
      K(RR006B) = 1.2083246981D+05 &
	   * exp(1.2253 * LT - 75469260.17 / RT)
      K(RR007F) = 7.4500000000D+40 &
	   * exp(-10.13 * LT - 77500000 / RT)
      K(RR007B) = 1.5255854924D+50 &
	   * exp(-11.5465 * LT - 186098173.3 / RT)
      K(RR008F) = 7.8000000000D+10
      K(RR008B) = 5.6362388680D+19 &
	   * exp(-1.29954 * LT - 250786354.2 / RT)
      K(RR009F) = 1.0000000000D+08 * exp(-12550000 / RT)
      K(RR009B) = 1.5287954728D+07 &
	   * exp(0.204893 * LT - 189399422.9 / RT)
      K(RR010F) = 1.2100000000D+07 * exp(2490000 / RT)
      K(RR010B) = 9.7578085241D+07 &
	   * exp(0.121092 * LT - 97776912.53 / RT)
      K(RR011F) = 9.0000000000D+10
      K(RR011B) = 9.4447045999D+12 &
	   * exp(0.0970785 * LT - 400745685 / RT)
      K(RR012F) = 1.8000000000D+10
      K(RR012B) = 6.5741679900D+24 &
	   * exp(-2.33737 * LT - 433035302.1 / RT)
      K(RR013F) = 9.0300000000D+09 * exp(3200000 / RT)
      K(RR013B) = 1.7706961729D+13 &
	   * exp(-0.341857 * LT - 292002157.9 / RT)
      K(RR014F) = 4.0400000000D+42 &
	   * exp(-7.67 * LT - 467900000 / RT)
      K(RR014B) = 3.9209380170D+30 &
	   * exp(-6.00376 * LT - 34921721.89 / RT)
      K(RR015F) = 1.9300000000D+15 &
	   * exp(-1.25 * LT - 32090000 / RT)
      K(RR015B) = 1.0584325883D+24 &
	   * exp(-3.18571 * LT - 101188124 / RT)
      K(RR016F) = 5.9300000000D+51 &
	   * exp(-11.76 * LT - 98530000 / RT)
      K(RR016B) = 1.1141410801D+55 &
	   * exp(-11.4905 * LT - 462410154.1 / RT)
      K(RR017F) = 2.4100000000D+10
      K(RR017B) = 3.0422905514D+16 &
	   * exp(-1.37788 * LT - 159986443 / RT)
      K(RR018F) = 5.0000000000D+10
      K(RR018B) = 6.3084963399D+04 &
	   * exp(1.41602 * LT - 18103989.21 / RT)
      K(RR019F) = 5.0000000000D+10
      K(RR019B) = 4.4314188478D+06 &
	   * exp(1.36084 * LT - 757475726.4 / RT)
      K(RR020F) = 2.0000000000D+10
      K(RR020B) = 1.2649965411D+04 &
	   * exp(1.04675 * LT - 330473656.1 / RT)
      K(RR021F) = 2.0000000000D+10
      K(RR021B) = 2.5283980856D+01 &
	   * exp(1.48246 * LT - 259653260.4 / RT)
      K(RR022F) = 5.0000000000D+10
      K(RR022B) = 3.7594731564D+11 &
	   * exp(0.267581 * LT - 384856198.1 / RT)
      K(RR023F) = 3.0000000000D+10
      K(RR023B) = 1.3558457941D+10 &
	   * exp(0.27143 * LT - 165356182 / RT)
      K(RR024F) = 1.0000000000D+10
      K(RR024B) = 3.1185275273D+07 &
	   * exp(0.606272 * LT - 131240969.4 / RT)
      K(RR025F) = 5.0000000000D+09
      K(RR025B) = 1.3533569131D+07 &
	   * exp(0.605334 * LT - 113774621.6 / RT)
      K(RR026F) = 1.0000000000D+09
      K(RR026B) = 2.4177349976D+09 &
	   * exp(0.352248 * LT - 405043602.2 / RT)
      K(RR027F) = 3.6000000000D+10 * exp(-46020000 / RT)
      K(RR027B) = 1.1575885276D+10 &
	   * exp(0.084667 * LT - 66207404.11 / RT)
      K(RR028F) = 9.0000000000D+07
      K(RR028B) = 4.1464938961D+06 &
	   * exp(0.393202 * LT - 92132061.56 / RT)
      K(RR029F) = 7.5000000000D+09 * exp(-54390000 / RT)
      K(RR029B) = 1.0687278636D+10 &
	   * exp(0.0769611 * LT - 52331626.84 / RT)
      K(RR030F) = 4.9000000000D+11 * exp(-0.5 * LT)
      K(RR030B) = 8.1202932548D+08 &
	   * exp(0.176365 * LT - 4176769.025 / RT)
      K(RR031F) = 1.2000000000D+10
      K(RR031B) = 2.6332608966D+16 &
	   * exp(-1.47857 * LT - 108949436.1 / RT)
      K(RR032F) = 4.2200000000D+10 * exp(-259830000 / RT)
      K(RR032B) = 1.0936945048D+08 &
	   * exp(0.123795 * LT + 2021441.533 / RT)
      K(RR033) = 4.9000000000D+09 &
	   * exp(0.42 * LT - 317150000 / RT)
      K(RR034F) = 1.2000000000D+11
      K(RR034B) = 2.8148930110D+13 &
	   * exp(0.0237669 * LT - 355901345.9 / RT)
      K(RR035F) = 3.0000000000D+08
      K(RR035B) = 2.5874528983D+11 &
	   * exp(-0.197106 * LT - 217007102.4 / RT)
      K(RR036F) = 3.0000000000D+08
      K(RR036B) = 1.1066571938D+10 &
	   * exp(-0.135412 * LT - 214914712 / RT)
      K(RR037F) = 3.1000000000D+10
      K(RR037B) = 1.2648832430D+13 &
	   * exp(-0.591432 * LT - 112520403.1 / RT)
      K(RR038F) = 2.6100000000D-01 &
	   * exp(3.37 * LT - 66580000 / RT)
      K(RR038B) = 1.4479008027D-02 &
	   * exp(3.32222 * LT - 11157426.57 / RT)
      K(RR039F) = 1.3600000000D+11
      K(RR039B) = 5.8014628471D+18 &
	   * exp(-0.483536 * LT - 670190892.2 / RT)
      K(RR040F) = 1.0000000000D+10
      K(RR040B) = 5.3108215024D+09 &
	   * exp(0.514343 * LT - 313844820.3 / RT)
      K(RR041F) = 1.2500000000D+08 * exp(-4180000 / RT)
      K(RR041B) = 7.7263247114D+01 &
	   * exp(1.00075 * LT - 263245198.7 / RT)
      K(RR042F) = 5.0000000000D+10
      K(RR042B) = 3.1525378438D+18 &
	   * exp(-1.19135 * LT - 465232193.1 / RT)
      K(RR043F) = 5.0000000000D+10
      K(RR043B) = 2.8554902815D+17 &
	   * exp(-1.30881 * LT - 178465787 / RT)
      K(RR044F) = 5.0000000000D+09
      K(RR044B) = 1.0892535033D+17 &
	   * exp(-1.28933 * LT - 185700190.6 / RT)
      K(RR045F) = 1.0000000000D+10
      K(RR045B) = 4.8188337806D+10 &
	   * exp(0.18992 * LT - 287562460.5 / RT)
      K(RR046F) = 5.0000000000D+10
      K(RR046B) = 1.9579560457D+22 &
	   * exp(-1.55827 * LT - 485826868.1 / RT)
      K(RR200F) = 3.4600000000D+09 &
	   * exp(0.44 * LT - 22860000 / RT)
      K(RR200B) = 6.0360475259D+03 &
	   * exp(1.75196 * LT - 93515998.37 / RT)
      K(RR047) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RR048) = 2.9200000000D+09 * exp(-7570000 / RT)
      K(RR049) = 3.0100000000D+10 * exp(-163800000 / RT)
      K(RR050) = 2.3400000000D+07 &
	   * exp(0.73 * LT + 4660000 / RT)
      K(RR051) = 3.0100000000D+09 * exp(-49890000 / RT)
      K(RR052) = 2.7200000000D+03 &
	   * exp(1.77 * LT - 24770000 / RT)
      K0TROE = 2.8000000000D+24 &
	   * exp(-3.86 * LT - 13890000 / RT)
      KINFTROE = 1.0200000000D+10 &
	   * exp(0.27 * LT - 1170000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(RR053F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 3.8610572228D+29 &
	   * exp(-3.97975 * LT - 439116845.8 / RT)
      KINFTROE = 1.4065279883D+15 &
	   * exp(0.150251 * LT - 426396845.8 / RT)
      K(RR053B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RR054F) = 1.1000000000D+07 &
	   * exp(1.13 * LT - 58280000 / RT)
      K(RR054B) = 1.8645069743D+04 &
	   * exp(1.59998 * LT - 66144642.01 / RT)
      K(RR055F) = 7.9400000000D+26 &
	   * exp(-5.06 * LT - 20340000 / RT)
      K(RR055B) = 8.1610655573D+31 &
	   * exp(-5.06898 * LT - 393338860.7 / RT)
      K(RR056F) = 3.1600000000D+26 &
	   * exp(-5 * LT - 19710000 / RT)
      K(RR056B) = 1.2666373930D+32 &
	   * exp(-5.11373 * LT - 387807163.3 / RT)
      K(RR057F) = 7.5300000000D+03 &
	   * exp(1.55 * LT - 8810000 / RT)
      K(RR057B) = 2.3470297119D+15 &
	   * exp(0.0186034 * LT - 477264514.5 / RT)
      K(RR058F) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR058B) = 1.1463339595D+05 &
	   * exp(1.63305 * LT - 446954897.3 / RT)
      K(RR059F) = 1.1300000000D+02 &
	   * exp(2.28 * LT - 10320000 / RT)
      K(RR059B) = 2.1071838653D+00 &
	   * exp(2.64624 * LT - 81229166.75 / RT)
      K(RR060F) = 1.8800000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RR060B) = 1.0162158826D+33 &
	   * exp(-7.26513 * LT - 287763381.1 / RT)
      K(RR061F) = 1.3800000000D+11
      K(RR061B) = 4.2690327506D+13 &
	   * exp(-0.363787 * LT - 244964046.4 / RT)
      K(RR062F) = 1.7000000000D+02 &
	   * exp(1.7 * LT - 6280000 / RT)
      K(RR062B) = 1.2306728046D+00 &
	   * exp(2.34085 * LT - 349953693.6 / RT)
      K(RR063F) = 8.0000000000D+08
      K(RR063B) = 2.5061786983D-02 &
	   * exp(1.80032 * LT - 121159499.7 / RT)
      K(RR064F) = 3.0000000000D+08
      K(RR064B) = 8.3477672646D+10 &
	   * exp(-0.101893 * LT - 163429348 / RT)
      K(RR065F) = 3.0000000000D+08
      K(RR065B) = 2.1405800537D+10 &
	   * exp(0.00285485 * LT - 168331045.4 / RT)
      K(RR066F) = 4.0000000000D+11 * exp(-175430000 / RT)
      K(RR066B) = 7.6438171011D-03 &
	   * exp(2.23323 * LT - 275388910.5 / RT)
      K(RR067F) = 2.5000000000D+10
      K(RR067B) = 1.8919886866D+12 &
	   * exp(0.11898 * LT - 302323591.4 / RT)
      K(RR068F) = 2.5000000000D+10
      K(RR068B) = 4.8515406766D+11 &
	   * exp(0.223728 * LT - 307225288.8 / RT)
      K(RR069F) = 5.0000000000D+10
      K(RR069B) = 1.2426252739D+16 &
	   * exp(-1.23168 * LT - 222587604.3 / RT)
      K(RR070F) = 5.0000000000D+10
      K(RR070B) = 7.7450588497D+17 &
	   * exp(-1.33185 * LT - 224722517.8 / RT)
      K(RR071F) = 7.7600000000D+39 &
	   * exp(-7.8 * LT - 328220000 / RT)
      K(RR071B) = 1.9898615631D+39 &
	   * exp(-7.69525 * LT - 333121697.4 / RT)
      K(RR073F) = 2.4700000000D+12 &
	   * exp(-0.33 * LT - 26930000 / RT)
      K(RR073B) = 6.3337088411D+11 &
	   * exp(-0.225252 * LT - 31831697.35 / RT)
      K(RR074F) = 2.0100000000D+46 &
	   * exp(-10.77 * LT - 82100000 / RT)
      K(RR074B) = 3.2976487096D+51 &
	   * exp(-11.2885 * LT - 326398622 / RT)
      K(RR075F) = 6.7000000000D+39 &
	   * exp(-12.46 * LT - 68450000 / RT)
      K(RR075B) = 7.2625724752D+42 &
	   * exp(-12.5238 * LT - 222514367.3 / RT)
      K(RR076F) = 8.8300000000D+49 &
	   * exp(-12.36 * LT - 68810000 / RT)
      K(RR076B) = 3.7326325453D+53 &
	   * exp(-12.5285 * LT - 217972669.9 / RT)
      K(RR077F) = 1.5300000000D+46 &
	   * exp(-11.97 * LT - 59180000 / RT)
      K(RR077B) = 1.6495716481D+50 &
	   * exp(-12.246 * LT - 191275736.1 / RT)
      K(RR078F) = 8.5000000000D+01 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RR078B) = 1.9329120587D-01 &
	   * exp(3.05922 * LT - 84112627.12 / RT)
      K(RR079F) = 4.4900000000D+04 &
	   * exp(1.92 * LT - 23810000 / RT)
      K(RR079B) = 6.1211557595D+01 &
	   * exp(2.24307 * LT - 77813209.56 / RT)
      K(RR080F) = 8.0500000000D+02 &
	   * exp(2.22 * LT - 3100000 / RT)
      K(RR080B) = 2.0139176673D+01 &
	   * exp(2.47547 * LT - 126237151.9 / RT)
      K(RR081F) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RR081B) = 1.7839562769D+12 &
	   * exp(-0.16099 * LT - 162124535 / RT)
      K(RR082F) = 1.3000000000D-01 &
	   * exp(3.37 * LT - 66580000 / RT)
      K(RR082B) = 8.7173128474D-02 &
	   * exp(3.12226 * LT - 59833483.64 / RT)
      K(RR083F) = 1.3300000000D+03 &
	   * exp(2.53 * LT - 51210000 / RT)
      K(RR083B) = 7.7554312548D-01 &
	   * exp(2.99397 * LT - 116204324.5 / RT)
      K(RR084F) = 1.3100000000D-04 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RR084B) = 8.4038495518D-07 &
	   * exp(4.56022 * LT - 124438849.2 / RT)
      K(RR085F) = 2.2700000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RR085B) = 2.4607004806D+02 &
	   * exp(1.94376 * LT - 112396232.4 / RT)
      K(RR086F) = 9.7600000000D+07 &
	   * exp(0.12 * LT - 97780000 / RT)
      K(RR086B) = 1.6782252051D+07 &
	   * exp(-0.0229938 * LT - 95935180.99 / RT)
      K(RR087F) = 9.6300000000D+03 &
	   * exp(2.05 * LT - 750000 / RT)
      K(RR087B) = 1.2470015987D-01 &
	   * exp(3.02862 * LT - 90776201.99 / RT)
      K(RR088F) = 4.0500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR088B) = 6.5763711850D-02 &
	   * exp(2.95854 * LT - 113261908.7 / RT)
      K(RR089F) = 6.2500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR089B) = 7.6307688077D-01 &
	   * exp(3.22612 * LT - 498118106.9 / RT)
      K(RR090F) = 1.0000000000D+10
      K(RR090B) = 1.7042930120D+10 &
	   * exp(-0.132572 * LT - 188385703.1 / RT)
      K(RR091F) = 1.0000000000D+10
      K(RR091B) = 6.6463486802D+10 &
	   * exp(-0.237321 * LT - 183484005.8 / RT)
      K(RR092F) = 2.4100000000D+03 &
	   * exp(2 * LT - 53190000 / RT)
      K(RR092B) = 2.2893692008D+03 &
	   * exp(1.98652 * LT - 49772058.13 / RT)
      K(RR093F) = 7.5300000000D+03 &
	   * exp(1.55 * LT - 8810000 / RT)
      K(RR093B) = 2.6607561038D+02 &
	   * exp(1.94752 * LT - 131181607.4 / RT)
      K(RR094F) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR094B) = 2.0838122248D+03 &
	   * exp(1.52072 * LT - 225293668.9 / RT)
      K(RR095) = 4.0900000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RR096) = 5.8400000000D+09 * exp(-7570000 / RT)
      K(RR097) = 2.8900000000D+05 &
	   * exp(1.35 * LT + 6580000 / RT)
      K(RR098) = 4.0900000000D+01 &
	   * exp(2.5 * LT - 42690000 / RT)
      K(RR099) = 3.4900000000D-11 &
	   * exp(6.21 * LT - 6820000 / RT)
      K(RR100F) = 7.0600000000D+56 &
	   * exp(-14.08 * LT - 317430000 / RT)
      K(RR100B) = 4.6645746281D+54 &
	   * exp(-13.6253 * LT - 227195745.3 / RT)
      K(RR101F) = 5.0000000000D+51 &
	   * exp(-13.02 * LT - 306690000 / RT)
      K(RR101B) = 8.4256335535D+49 &
	   * exp(-12.6727 * LT - 199388811.5 / RT)
      K(RR102F) = 1.5000000000D+48 &
	   * exp(-12.71 * LT - 225520000 / RT)
      K(RR102B) = 3.8257490317D+48 &
	   * exp(-12.8174 * LT - 208453066.2 / RT)
      K(RR103F) = 9.5600000000D+00 &
	   * exp(2.8 * LT - 13770000 / RT)
      K(RR103B) = 1.3619748690D-02 &
	   * exp(3.66875 * LT - 202562865.8 / RT)
      K(RR104F) = 6.0300000000D+09
      K(RR104B) = 9.4510742621D+07 &
	   * exp(0.765006 * LT - 251837390.6 / RT)
      K(RR105F) = 4.8600000000D+08 &
	   * exp(-0.32 * LT + 550000 / RT)
      K(RR105B) = 1.2871407986D+09 &
	   * exp(0.0285441 * LT - 197154773.7 / RT)
      K(RR106F) = 2.4100000000D+09
      K(RR106B) = 8.1647196307D+09 &
	   * exp(0.382884 * LT - 222220634.9 / RT)
      K(RR107F) = 9.6400000000D+08 * exp(550000 / RT)
      K(RR107B) = 7.3002216160D+09 &
	   * exp(0.309573 * LT - 176826295.8 / RT)
      K(RR108F) = 1.0000000000D+09
      K(RR108B) = 1.1451889853D+07 &
	   * exp(0.787987 * LT - 119581532.1 / RT)
      K(RR109F) = 2.0600000000D+01 &
	   * exp(2.19 * LT - 73600000 / RT)
      K(RR109B) = 1.8087353764D-01 &
	   * exp(2.69668 * LT - 33969193.32 / RT)
      K(RR110F) = 3.3600000000D+02 &
	   * exp(1.81 * LT - 80290000 / RT)
      K(RR110B) = 4.4602469311D+02 &
	   * exp(1.66065 * LT - 292530403.7 / RT)
      K(RR111) = 9.7100000000D+17 &
	   * exp(-2.7 * LT - 104520000 / RT)
      K(RR112F) = 3.0800000000D+06 &
	   * exp(0.37 * LT - 70750000 / RT)
      K(RR112B) = 7.0070474526D+03 &
	   * exp(1.028 * LT - 325362179.3 / RT)
      K(RR113F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RR113B) = 8.6220550094D+15 &
	   * exp(-0.333429 * LT - 321281479.9 / RT)
      K(RR114F) = 4.2000000000D+29 &
	   * exp(-5.16 * LT - 126050000 / RT)
      K(RR114B) = 4.6528481728D+32 &
	   * exp(-5.70891 * LT - 415200216.9 / RT)
      K(RR115F) = 6.0000000000D+10
      K(RR115B) = 2.1283843105D+10 &
	   * exp(0.502181 * LT - 298106582.2 / RT)
      K(RR116F) = 2.6600000000D+09
      K(RR116B) = 3.4693729364D+09 &
	   * exp(0.281308 * LT - 159212338.8 / RT)
      K(RR117F) = 1.0600000000D+13 &
	   * exp(-0.94 * LT - 10560000 / RT)
      K(RR117B) = 4.0003442177D+12 &
	   * exp(-0.855882 * LT - 58003268.83 / RT)
      K(RR118F) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RR118B) = 9.1127053346D+06 &
	   * exp(0.902702 * LT - 60093366.86 / RT)
      K(RR119F) = 3.3400000000D+09
      K(RR119B) = 1.8467642799D+08 &
	   * exp(0.518782 * LT - 283928817.9 / RT)
      K(RR120F) = 6.0000000000D+10
      K(RR120B) = 7.0278140237D+07 &
	   * exp(0.880152 * LT - 400211007.7 / RT)
      K(RR121) = 5.0000000000D+09
      K(RR122F) = 2.0000000000D+10
      K(RR122B) = 3.2504162150D+01 &
	   * exp(1.3277 * LT - 124722796.7 / RT)
      K(RR123F) = 9.0000000000D+10
      K(RR123B) = 4.8320783021D+12 &
	   * exp(0.0474627 * LT - 388340836.9 / RT)
      K(RR124F) = 1.0000000000D+08
      K(RR124B) = 1.0278827655D+10 &
	   * exp(-0.00142624 * LT - 292840725.8 / RT)
      K(RR125F) = 3.3400000000D+09
      K(RR125B) = 7.2407949316D+07 &
	   * exp(0.62621 * LT - 300995751.7 / RT)
      K(RR126F) = 6.0000000000D+10
      K(RR126B) = 3.5987068607D+06 &
	   * exp(1.26938 * LT - 423845942.6 / RT)
      K(RR127) = 5.0000000000D+09
      K(RR128F) = 2.0000000000D+10
      K(RR128B) = 1.6644286678D+00 &
	   * exp(1.71693 * LT - 148357731.6 / RT)
      K(RR129F) = 9.0000000000D+10
      K(RR129B) = 1.8945616645D+12 &
	   * exp(0.15489 * LT - 405407770.7 / RT)
      K(RR130F) = 1.0000000000D+08
      K(RR130B) = 4.0301236059D+09 &
	   * exp(0.106001 * LT - 309907659.6 / RT)
      K(RR131F) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RR131B) = 4.5663206858D+02 &
	   * exp(1.76671 * LT - 53895145.38 / RT)
      K(RR132F) = 6.7000000000D+02 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RR132B) = 8.9518169802D+01 &
	   * exp(1.87414 * LT - 70962079.17 / RT)
      K(RR133) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RR134) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RR135) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RR136F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RR136B) = 1.1941331843D+11 &
	   * exp(-0.691616 * LT - 403226009 / RT)
      K(RR137F) = 1.9200000000D+04 &
	   * exp(1.02 * LT + 8510000 / RT)
      K(RR137B) = 2.5515359785D+04 &
	   * exp(1.07196 * LT - 42093448.02 / RT)
      K(RR138) = 1.0000000000D+09 * exp(-25100000 / RT)
      K(RR139F) = 1.0000000000D+14 * exp(-121750000 / RT)
      K(RR139B) = 2.4418106508D+11 &
	   * exp(-0.221635 * LT - 81879319.52 / RT)
      K(RR140) = 2.0300000000D+12 &
	   * exp(0.09 * LT - 98580000 / RT)
      K(RR141F) = 8.0000000000D+18 &
	   * exp(-2.39 * LT - 46780000 / RT)
      K(RR141B) = 4.3155004306D+12 &
	   * exp(-0.859397 * LT - 80320978.73 / RT)
      K(RR142F) = 1.2000000000D+05 &
	   * exp(1.6 * LT - 1370000 / RT)
      K(RR142B) = 7.2781615229D+01 &
	   * exp(2.14953 * LT - 101195463.3 / RT)
      K(RR143F) = 3.5000000000D+04 &
	   * exp(1.6 * LT + 4070000 / RT)
      K(RR143B) = 3.5179051144D-02 &
	   * exp(2.82589 * LT - 99932232.35 / RT)
      K(RR144F) = 6.6000000000D+02 &
	   * exp(2.54 * LT - 28270000 / RT)
      K(RR144B) = 8.2106582340D+01 &
	   * exp(2.62077 * LT - 97481333.73 / RT)
      K(RR145F) = 9.6500000000D+01 &
	   * exp(2.68 * LT - 15550000 / RT)
      K(RR145B) = 7.1970624796D+00 &
	   * exp(2.72461 * LT - 78671916.16 / RT)
      K(RR146F) = 2.0000000000D+05 &
	   * exp(1.46 * LT - 2250000 / RT)
      K(RR146B) = 2.7372640710D+05 &
	   * exp(1.43702 * LT - 134505858.5 / RT)
      K(RR147F) = 9.6000000000D+00 &
	   * exp(2.6 * LT - 58200000 / RT)
      K(RR147B) = 3.5216920977D+02 &
	   * exp(2.0738 * LT - 60572190.24 / RT)
      K(RR148F) = 4.5200000000D-04 &
	   * exp(3.65 * LT - 29930000 / RT)
      K(RR148B) = 1.0453243323D-01 &
	   * exp(3.21056 * LT - 108053241.6 / RT)
      K(RR149F) = 4.0000000000D+02 &
	   * exp(2.5 * LT - 40960000 / RT)
      K(RR149B) = 3.2877696006D-01 &
	   * exp(3.03548 * LT - 19937079.03 / RT)
      K(RR150F) = 6.0000000000D+07 &
	   * exp(0.7 * LT - 31920000 / RT)
      K(RR150B) = 2.9565590310D+04 &
	   * exp(1.19933 * LT - 4807661.464 / RT)
      K(RR151F) = 1.1000000000D+03 &
	   * exp(2 * LT - 6070000 / RT)
      K(RR151B) = 9.9468766212D+00 &
	   * exp(2.43174 * LT - 48091603.77 / RT)
      K(RR152F) = 8.4000000000D-04 &
	   * exp(3.5 * LT - 48790000 / RT)
      K(RR152B) = 1.2835099968D-03 &
	   * exp(3.51528 * LT - 36678986.93 / RT)
      K(RR153F) = 6.6500000000D+02 &
	   * exp(2.53 * LT - 51210000 / RT)
      K(RR153B) = 1.3940817680D+00 &
	   * exp(2.95806 * LT - 13120145.24 / RT)
      K(RR154F) = 1.2100000000D+08 &
	   * exp(0.7 * LT - 37490000 / RT)
      K(RR154B) = 1.5207082165D+05 &
	   * exp(1.0919 * LT + 6689272.325 / RT)
      K(RR155F) = 6.5500000000D-05 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RR155B) = 1.5106385521D-06 &
	   * exp(4.52431 * LT - 21354669.98 / RT)
      K(RR156F) = 1.1400000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RR156B) = 4.4427312105D+02 &
	   * exp(1.90785 * LT - 9312053.143 / RT)
      K(RH01F) = 1.0000000000D+10 * exp(3160000 / RT)
      K(RH01B) = 1.5210017608D-01 &
	   * exp(1.77683 * LT - 443157004 / RT)
      K(RH02F) = 6.0000000000D+10
      K(RH02B) = 2.7921110540D+18 &
	   * exp(-0.788396 * LT - 551409473.7 / RT)
      K(RH03F) = 3.2000000000D+06 &
	   * exp(1.8 * LT - 125970000 / RT)
      K(RH03B) = 1.6072671724D+01 &
	   * exp(2.93863 * LT - 7652014.074 / RT)
      K(RH04F) = 4.0000000000D+11 * exp(-225520000 / RT)
      K(RH04B) = 2.2060311607D+17 &
	   * exp(-0.552714 * LT - 398918197.4 / RT)
      K(RH05) = 1.5100000000D+11 * exp(-234300000 / RT)
      K(RH06F) = 1.5100000000D+10 * exp(-178660000 / RT)
      K(RH06B) = 2.5913304529D+13 &
	   * exp(-0.924094 * LT - 205722336 / RT)
      K(RH07F) = 9.5600000000D+09 * exp(-130120000 / RT)
      K(RH07B) = 4.2262757709D+05 &
	   * exp(0.751128 * LT - 236628846.9 / RT)
      K(RH08F) = 2.0600000000D+04 &
	   * exp(2 * LT - 7950000 / RT)
      K(RH08B) = 2.2950590346D-02 &
	   * exp(3.13616 * LT - 282089544 / RT)
      K0TROE = 2.3000000000D+39 &
	   * exp(-8.1 * LT - 10490000 / RT)
      KINFTROE = 4.3100000000D+07 &
	   * exp(1.16 * LT - 7330000 / RT)
      FCTROE = 0.0748 * EXP( -TEMP / -4216 )
      K(RH09F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.2501321284D+42 &
	   * exp(-8.26008 * LT - 193072256.9 / RT)
      KINFTROE = 2.3426389015D+10 &
	   * exp(0.999916 * LT - 189912256.9 / RT)
      K(RH09B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RH10F) = 1.3700000000D+36 &
	   * exp(-7.87 * LT - 64610000 / RT)
      K(RH10B) = 4.7217121748D+38 &
	   * exp(-7.9297 * LT - 199616108.7 / RT)
      K(RH11F) = 9.1500000000D+06 &
	   * exp(1.03 * LT - 90990000 / RT)
      K(RH11B) = 5.0560557267D+02 &
	   * exp(2.06488 * LT - 35716538.81 / RT)
      K(RH12F) = 3.3000000000D+09 &
	   * exp(-0.25 * LT - 9940000 / RT)
      K(RH12B) = 3.6180564088D+06 &
	   * exp(0.452329 * LT - 282304319.6 / RT)
      K(RH13F) = 4.1000000000D+43 &
	   * exp(-9.5 * LT - 221750000 / RT)
      K(RH13B) = 6.4659597175D+43 &
	   * exp(-9.60039 * LT - 269326148.2 / RT)
      K(RH14F) = 2.5000000000D+17 &
	   * exp(-1.67 * LT - 45190000 / RT)
      K(RH14B) = 3.9426583642D+17 &
	   * exp(-1.77039 * LT - 92766148.25 / RT)
      K(RH15F) = 2.0000000000D+44 &
	   * exp(-10.26 * LT - 54680000 / RT)
      K(RH15B) = 7.4803555018D+49 &
	   * exp(-10.4028 * LT - 526163576.5 / RT)
      K(RH16F) = 3.4000000000D+40 &
	   * exp(-9.01 * LT - 50710000 / RT)
      K(RH16B) = 8.0634708732D+45 &
	   * exp(-9.0524 * LT - 474617428.3 / RT)
      K(RH17F) = 6.3000000000D+22 &
	   * exp(-3.34 * LT - 41880000 / RT)
      K(RH17B) = 5.9560076535D+19 &
	   * exp(-2.28677 * LT - 347576512.3 / RT)
      K(RH18F) = 2.8000000000D+20 &
	   * exp(-2.55 * LT - 45100000 / RT)
      K(RH18B) = 1.6785086788D+17 &
	   * exp(-1.39638 * LT - 303220364 / RT)
      K(RH19F) = 1.5000000000D+10
      K(RH19B) = 1.0172597244D+10 &
	   * exp(0.409933 * LT - 298085379.1 / RT)
      K(RH20F) = 3.0000000000D+10
      K(RH20B) = 1.2900683124D+10 &
	   * exp(0.510318 * LT - 250509230.9 / RT)
      K(RH21F) = 2.5000000000D+09
      K(RH21B) = 1.8652337345D+10 &
	   * exp(0.306186 * LT - 361129903.9 / RT)
      K(RH22F) = 5.0000000000D+09
      K(RH22B) = 2.3654518883D+10 &
	   * exp(0.406572 * LT - 313553755.6 / RT)
      K(RH23F) = 6.7000000000D+02 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RH23B) = 2.8003479993D+03 &
	   * exp(1.65786 * LT - 68051706.63 / RT)
      K(RH24F) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RH24B) = 3.5513449818D+03 &
	   * exp(1.75825 * LT - 20475558.39 / RT)
      K(RH25F) = 2.0000000000D+10
      K(RH25B) = 1.6504795180D+07 &
	   * exp(0.970973 * LT - 227513612.4 / RT)
      K(RH26F) = 1.6300000000D+08 * exp(7530000 / RT)
      K(RH26B) = 4.0425156849D+06 &
	   * exp(0.664277 * LT - 360458358.8 / RT)
      K(RH27F) = 1.7000000000D+02 &
	   * exp(1.7 * LT - 6280000 / RT)
      K(RH27B) = 3.6988379195D-10 &
	   * exp(3.9463 * LT - 340861711.7 / RT)
      K(RH28F) = 1.2700000000D+02 &
	   * exp(2.75 * LT - 48740000 / RT)
      K(RH28B) = 7.9365363993D-02 &
	   * exp(3.24302 * LT - 10347911.29 / RT)
      K(RH29F) = 6.3500000000D+01 &
	   * exp(2.75 * LT - 48740000 / RT)
      K(RH29B) = 6.2582103235D-02 &
	   * exp(3.14263 * LT - 57924059.54 / RT)
      K(RH30F) = 6.5500000000D-05 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RH30B) = 4.5032001812D-07 &
	   * exp(4.58927 * LT - 21052436.03 / RT)
      K(RH31F) = 3.2800000000D-05 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RH31B) = 3.5563372177D-07 &
	   * exp(4.48888 * LT - 68628584.28 / RT)
      K(RH32F) = 1.1400000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RH32B) = 1.3243742499D+02 &
	   * exp(1.97281 * LT - 9009819.194 / RT)
      K(RH33F) = 5.6800000000D+01 &
	   * exp(2 * LT - 38490000 / RT)
      K(RH33B) = 1.0406467934D+02 &
	   * exp(1.87242 * LT - 56585967.44 / RT)
      K(RH34F) = 6.2500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RH34B) = 2.9857187952D+00 &
	   * exp(3.10516 * LT - 468923868 / RT)
      K(RH35F) = 3.5800000000D+01 &
	   * exp(2.47 * LT - 3890000 / RT)
      K(RH35B) = 2.3000634086D-04 &
	   * exp(3.47671 * LT - 106918837.3 / RT)
      K(RH36F) = 1.9500000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RH36B) = 1.2309062294D+00 &
	   * exp(2.34618 * LT - 162360276.5 / RT)
      K(RH37F) = 7.8000000000D+10
      K(RH37B) = 2.7708034384D+17 &
	   * exp(-1.42196 * LT - 114055959.4 / RT)
      K(RH38F) = 7.8000000000D+10
      K(RH38B) = 1.8874596234D+19 &
	   * exp(-1.96405 * LT - 108982566.6 / RT)
      K(RH39F) = 7.8000000000D+10
      K(RH39B) = 2.3090595629D+17 &
	   * exp(-1.33851 * LT - 152178000.8 / RT)
      K(RH40F) = 7.8000000000D+10
      K(RH40B) = 2.6650315696D+19 &
	   * exp(-2.06272 * LT - 145380321.9 / RT)
      K(RB00F) = 1.0000000000D+09
      K(RB00B) = 6.8182357335D+24 &
	   * exp(-3.01231 * LT - 364611266.9 / RT)
      K(RB01F) = 1.9000000000D+11
      K(RB01B) = 3.5052571308D+22 &
	   * exp(-1.7372 * LT - 354038492.4 / RT)
      K(RB02F) = 1.3200000000D+09 &
	   * exp(0.16 * LT - 34780000 / RT)
      K(RB02B) = 1.8479981673D+19 &
	   * exp(-1.51022 * LT - 204229236.9 / RT)
      K(RB04F) = 8.4300000000D+10
      K(RB04B) = 3.0917302550D+27 &
	   * exp(-2.77864 * LT - 496077857.7 / RT)
      K(RB05F) = 1.2000000000D+19 &
	   * exp(-2.44 * LT - 57130000 / RT)
      K(RB05B) = 4.1805841138D+28 &
	   * exp(-4.62745 * LT - 125794777.8 / RT)
      K(RB06F) = 2.4000000000D+17 &
	   * exp(-2.04 * LT - 64280000 / RT)
      K(RB06B) = 1.5163415165D+25 &
	   * exp(-3.8821 * LT - 86927999.16 / RT)
      K(RB07F) = 9.6000000000D+08
      K(RB07B) = 2.4080304685D+12 &
	   * exp(-0.307517 * LT - 319718019.1 / RT)
      K0TROE = 2.6000000000D+51 &
	   * exp(-11.94 * LT - 40890000 / RT)
      KINFTROE = 1.5000000000D+09
      FCTROE = 0.825 * EXP( -TEMP / 1340.6 ) &
	   + 0.175 * EXP( -TEMP / 60000 ) &
	   + 1 * EXP( -9769.8 / TEMP )
      K(RB08F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.8696668182D+66 &
	   * exp(-14.3484 * LT - 425256613.8 / RT)
      KINFTROE = 1.0786539336D+24 &
	   * exp(-2.40843 * LT - 384366613.8 / RT)
      K(RB08B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RB09F) = 7.2300000000D+08 * exp(-20920000 / RT)
      K(RB09B) = 2.5734807876D+13 &
	   * exp(-1.1124 * LT - 84019579.62 / RT)
      K(RB10F) = 5.7000000000D+36 &
	   * exp(-6.27 * LT - 470090000 / RT)
      K(RB10B) = 5.4144775118D+29 &
	   * exp(-5.67882 * LT - 42676920.07 / RT)
      K(RB11F) = 5.3000000000D+44 &
	   * exp(-8.62 * LT - 517180000 / RT)
      K(RB11B) = 9.1303543420D+35 &
	   * exp(-7.68347 * LT - 43750141.44 / RT)
      K(RB12F) = 2.5000000000D+15 * exp(-396220000 / RT)
      K(RB12B) = 1.2837469486D+05 &
	   * exp(1.58909 * LT - 205418218.7 / RT)
      K(RB14F) = 8.9400000000D+04 &
	   * exp(1.14 * LT - 51800000 / RT)
      K(RB14B) = 6.2546358130D+14 &
	   * exp(-1.25945 * LT - 63167753.16 / RT)
      K(RB15F) = 2.8300000000D+05 &
	   * exp(1.06 * LT - 46700000 / RT)
      K(RB15B) = 5.0770575440D+14 &
	   * exp(-1.2347 * LT - 62969450.52 / RT)
      K(RB16F) = 1.3300000000D+03 &
	   * exp(2.53 * LT - 51210000 / RT)
      K(RB16B) = 5.3552893087D-04 &
	   * exp(3.81677 * LT - 10871629.23 / RT)
      K(RB17F) = 6.6500000000D+02 &
	   * exp(2.53 * LT - 38660000 / RT)
      K(RB17B) = 1.4764640528D-02 &
	   * exp(3.47142 * LT - 44338407.86 / RT)
      K(RB18F) = 2.2000000000D+08
      K(RB18B) = 9.1137574052D+14 &
	   * exp(-1.25062 * LT - 46427788.33 / RT)
      K(RB19F) = 7.5000000000D+03 &
	   * exp(1.9 * LT - 15650000 / RT)
      K(RB19B) = 9.9828929904D-02 &
	   * exp(2.80526 * LT - 15238990.3 / RT)
      K(RB20F) = 6.2000000000D+03 &
	   * exp(2 * LT - 14350000 / RT)
      K(RB20B) = 2.7464749857D-02 &
	   * exp(3.18302 * LT - 37056153.97 / RT)
      K(RB21F) = 3.1000000000D+03 &
	   * exp(2 * LT - 1800000 / RT)
      K(RB21B) = 7.5720868741D-01 &
	   * exp(2.83767 * LT - 70522932.6 / RT)
      K(RB22F) = 2.0000000000D+11 * exp(-95540000 / RT)
      K(RB22B) = 1.4970623133D+08 &
	   * exp(0.766561 * LT - 64113537.13 / RT)
      K(RB23F) = 1.0000000000D+11 * exp(-82840000 / RT)
      K(RB23B) = 4.1274309621D+09 &
	   * exp(0.421208 * LT - 97430315.77 / RT)
      K(RB24F) = 5.0000000000D+10 * exp(-95540000 / RT)
      K(RB24B) = 4.7875611343D+07 &
	   * exp(0.800901 * LT - 88629398.27 / RT)
      K(RB25F) = 2.5000000000D+10 * exp(-82840000 / RT)
      K(RB25B) = 1.3199402512D+09 &
	   * exp(0.455548 * LT - 121946176.9 / RT)
      K(RB28) = 7.6600000000D+06 &
	   * exp(0.88 * LT - 4770000 / RT)
      K(RB29F) = 7.1500000000D+01 &
	   * exp(2.47 * LT - 3890000 / RT)
      K(RB29B) = 1.4303992598D-06 &
	   * exp(4.58683 * LT - 281251274.7 / RT)
      K(RB30F) = 3.8900000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RB30B) = 3.5475897840D-02 &
	   * exp(2.95278 * LT - 150862792.8 / RT)
      K(RB31F) = 3.7500000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RB31B) = 3.4254390076D+29 &
	   * exp(-6.62052 * LT - 119296129.1 / RT)
      K(RB32F) = 1.3000000000D+48 &
	   * exp(-11.92 * LT - 69040000 / RT)
      K(RB32B) = 1.0193779458D+52 &
	   * exp(-12.2223 * LT - 219503410.5 / RT)
      K(RB33F) = 4.9000000000D+48 &
	   * exp(-11.92 * LT - 74060000 / RT)
      K(RB33B) = 2.1186435579D+54 &
	   * exp(-12.5677 * LT - 270540189.1 / RT)
      K(RB34F) = 1.5000000000D+67 &
	   * exp(-16.89 * LT - 247300000 / RT)
      K(RB34B) = 8.2710604472D+68 &
	   * exp(-17.2354 * LT - 293316778.6 / RT)
      K(RB35F) = 3.1000000000D+23 &
	   * exp(-3.35 * LT - 72900000 / RT)
      K(RB35B) = 1.7093524924D+25 &
	   * exp(-3.69535 * LT - 118916778.6 / RT)
      K(RB36F) = 1.5000000000D+10
      K(RB36B) = 4.4711466211D+08 &
	   * exp(0.652554 * LT - 282628077.3 / RT)
      K(RB37F) = 2.0000000000D+09
      K(RB37B) = 6.5585873964D+08 &
	   * exp(0.548808 * LT - 345672602 / RT)
      K(RB38F) = 5.0000000000D+09
      K(RB38B) = 5.4798862760D+14 &
	   * exp(-0.703823 * LT - 407656286.7 / RT)
      K(RB39F) = 1.2100000000D+07 * exp(2490000 / RT)
      K(RB39B) = 1.0190792609D+11 &
	   * exp(-0.679809 * LT - 104687514.3 / RT)
      K(RB40F) = 6.0000000000D+08
      K(RB40B) = 2.4178165332D+14 &
	   * exp(-0.924696 * LT - 268762043.3 / RT)
      K(RB99F) = 1.0300000000D+10 &
	   * exp(0.21 * LT + 1790000 / RT)
      K(RB99B) = 1.0294910143D+08 &
	   * exp(1.09896 * LT - 553019079.5 / RT)
      K(RB41F) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RB41B) = 2.4616656286D+02 &
	   * exp(1.90048 * LT - 52594404.81 / RT)
      K(RB42) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RB43F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RB43B) = 2.2969084277D+13 &
	   * exp(-1.19719 * LT - 410070984.8 / RT)
      K(RB44F) = 3.0000000000D+10
      K(RB44B) = 1.6217315723D+07 &
	   * exp(0.997907 * LT - 236611298.7 / RT)
      K(RB45F) = 2.0000000000D+10 * exp(-8370000 / RT)
      K(RB45B) = 2.9279079883D+02 &
	   * exp(1.81725 * LT - 51416466.09 / RT)
      K(RB46F) = 4.0000000000D+09
      K(RB46B) = 2.3788681409D+07 &
	   * exp(0.894161 * LT - 299655823.4 / RT)
      K(RB47F) = 5.0000000000D+09
      K(RB47B) = 9.9380598978D+12 &
	   * exp(-0.35847 * LT - 361639508.1 / RT)
      K(RB48F) = 6.0000000000D+08
      K(RB48B) = 4.3848365309D+12 &
	   * exp(-0.579343 * LT - 222745264.6 / RT)
      K(RB49F) = 1.2100000000D+07 * exp(2490000 / RT)
      K(RB49B) = 1.8481534515D+09 &
	   * exp(-0.334456 * LT - 58670735.62 / RT)
      K(RB50F) = 2.1600000000D+07 * exp(-10460000 / RT)
      K(RB50B) = 2.3456044263D+05 &
	   * exp(0.591536 * LT - 355851374.7 / RT)
      K(RB51F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RB51B) = 6.4953481542D+04 &
	   * exp(1.56463 * LT - 330125493.9 / RT)
      K(RB52F) = 1.8400000000D-16 &
	   * exp(7.07 * LT + 15110000 / RT)
      K(RB52B) = 4.1737418083D-06 &
	   * exp(6.131 * LT - 570010743.3 / RT)
      K(RHP00) = 8.1000000000D+77 &
	   * exp(-17.62 * LT - 503750000 / RT)
      K(RHP01) = 1.4200000000D+78 &
	   * exp(-17.71 * LT - 505010000 / RT)
      K(RHP02) = 1.7500000000D+03 &
	   * exp(2.6 * LT - 18250000 / RT)
      K(RHP03) = 1.7200000000D+02 &
	   * exp(2.81 * LT - 9460000 / RT)
      K(RHP04) = 7.4000000000D+05 &
	   * exp(1.5 * LT - 1080000 / RT)
      K(RHP05) = 2.8900000000D+10 &
	   * exp(0.2 * LT - 209660000 / RT)
      K(RHP06) = 7.5700000000D+09 &
	   * exp(0.21 * LT - 73780000 / RT)
      K(RHP08) = 1.2100000000D+78 &
	   * exp(-19.78 * LT - 163270000 / RT)
      K(RHP09) = 1.1600000000D+05 &
	   * exp(0.98 * LT + 6400000 / RT)
      K(RHP10) = 1.8900000000D+12 &
	   * exp(0.02 * LT - 116250000 / RT)
      K(RHP11) = 7.7300000000D+18 &
	   * exp(-1.75 * LT - 133780000 / RT)
      K(RHP12) = 2.5300000000D+18 &
	   * exp(-1.65 * LT - 132560000 / RT)
      K(RHP13) = 2.4900000000D+16 &
	   * exp(-1.18 * LT - 123500000 / RT)
      K(RHP14) = 1.4100000000D+15 &
	   * exp(-0.53 * LT - 156520000 / RT)
      K(RHP15) = 6.6600000000D+08 &
	   * exp(0.45 * LT - 8430000 / RT)
      K(RHP16) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP17) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP18) = 6.1800000000D+16 &
	   * exp(-1.36 * LT - 77480000 / RT)
      K(RHP19) = 2.2300000000D+15 &
	   * exp(-0.7 * LT - 77850000 / RT)
      K(RHP20) = 8.9200000000D+19 &
	   * exp(-2.03 * LT - 87880000 / RT)
      K(RHP21) = 1.1900000000D+22 &
	   * exp(-1.87 * LT - 308010000 / RT)
      K(RHP22) = 1.0800000000D+03 &
	   * exp(3.77 * LT - 278430000 / RT)
      K(RHP23F) = 3.7000000000D+10 * exp(-16320000 / RT)
      K(RHP23B) = 2.8612674099D+09 &
	   * exp(0.278353 * LT - 80220857.27 / RT)
      K(RHP24) = 3.0000000000D+10 * exp(-5150000 / RT)
      K(RHP25) = 9.8700000000D+19 &
	   * exp(-3.62 * LT - 12150000 / RT)
      K(RHP26) = 2.0700000000D+13 &
	   * exp(-1.65 * LT + 7630000 / RT)
      K(RHP27) = 2.5000000000D+13 * exp(-188280000 / RT)
      K(RHP28) = 2.5000000000D+13 * exp(-188280000 / RT)
      K(RHP29) = 7.4600000000D+21 &
	   * exp(-2.61 * LT - 134000000 / RT)
      K(RHP30) = 8.4600000000D+14 &
	   * exp(-0.47 * LT - 157390000 / RT)
      K(RHP31) = 7.1000000000D+09 &
	   * exp(0.12 * LT - 6110000 / RT)
      K(RHP32) = 3.1500000000D-19 &
	   * exp(8.84 * LT - 29730000 / RT)
      K(RHP33) = 9.1700000000D+20 &
	   * exp(-1.63 * LT - 309570000 / RT)
      K(RHP34) = 2.8000000000D+10 * exp(-16740000 / RT)
      K(RHP35) = 2.5400000000D+02 &
	   * exp(2.56 * LT + 4730000 / RT)
      K(RHP36) = 5.1200000000D+03 &
	   * exp(2 * LT + 1250000 / RT)
      K(RHP37) = 2.5000000000D+13 * exp(-188280000 / RT)
      K(RHP38) = 1.3400000000D+15 &
	   * exp(-0.52 * LT - 160330000 / RT)
      K(RHP39) = 8.6100000000D+17 &
	   * exp(-1.4 * LT - 162800000 / RT)
      K(RHP40) = 7.6800000000D+09 &
	   * exp(0.11 * LT - 6190000 / RT)
      K(RHP41) = 8.3900000000D+14 &
	   * exp(-0.64 * LT - 112570000 / RT)
      K(RHP42) = 1.0000000000D+31 &
	   * exp(-5.43 * LT - 177850000 / RT)
      K(RHP43) = 5.0000000000D+15 * exp(-297060000 / RT)
      K(RHP44) = 5.0000000000D+10 * exp(-16320000 / RT)
      K(RHP45) = 2.2500000000D+10 * exp(-9280000 / RT)
      K(RHP46) = 2.7000000000D+10 * exp(-138910000 / RT)
      K(RHP47) = 1.4000000000D+09 * exp(-62340000 / RT)
      K(RHP48) = 1.0000000000D+08 * exp(-30540000 / RT)
      K(RHP49) = 5.0000000000D+10
      K(RHP50) = 3.0000000000D+08
      K(RHP51) = 1.3000000000D+10 * exp(-3560000 / RT)
      K(RHP52) = 1.3000000000D+10 * exp(-3560000 / RT)
      K(RHP53) = 7.2300000000D+02 &
	   * exp(2.34 * LT + 4390000 / RT)
      K(RHP54) = 1.3000000000D+10 * exp(-3560000 / RT)
      K(RHP55) = 1.0000000000D+09
      K(RHP58) = 5.0000000000D+08
      K(RHP56) = 1.0000000000D+09
      K(RHP57) = 1.0000000000D+09
      K(RHP59) = 4.0000000000D+10 * exp(-17570000 / RT)
      K(RHP60) = 2.6900000000D+07 &
	   * exp(0.76 * LT + 1420000 / RT)
      K(RHP61) = 2.0000000000D+10 &
	   * exp(0.5 * LT - 176570000 / RT)
      K(RHP62) = 1.7000000000D+09 * exp(-35310000 / RT)
      K(RHP63) = 2.8000000000D+09 * exp(-56900000 / RT)
      K(RHP64F) = 5.0100000000D+31 &
	   * exp(-5.9 * LT - 162290000 / RT)
      K(RHP64B) = 3.8359636976D+31 &
	   * exp(-6.66828 * LT - 40378380.57 / RT)
      K(RHP65F) = 1.8800000000D+03 &
	   * exp(1.84 * LT - 12800000 / RT)
      K(RHP65B) = 1.6201772977D+14 &
	   * exp(-0.0347199 * LT - 164270220.3 / RT)
      K(RHP67) = 3.1600000000D+10
      K(RHP68) = 1.0000000000D+06
      K(RHP69) = 8.0000000000D+09
      K(RHP70) = 3.9800000000D+09
      K(RHP71) = 6.3100000000D+09
      K(RHP72) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP73) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP74) = 7.9400000000D+14 * exp(-79500000 / RT)
      K(RHP75) = 7.9400000000D+14 * exp(-79500000 / RT)
      K(RHP76) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP77) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RHP78) = 1.3900000000D+16 &
	   * exp(-0.87 * LT - 82720000 / RT)
      K(RHP79) = 2.5100000000D+14 * exp(-97910000 / RT)
      K(RHP86) = 1.0000000000D+09
      K(RHP87) = 1.0000000000D+09
      K(RHP88) = 1.0000000000D+09
      K(RIC00) = 1.0200000000D+49 &
	   * exp(-9.38 * LT - 404100000 / RT)
      K(RIC01) = 5.7500000000D+49 &
	   * exp(-9.66 * LT - 410200000 / RT)
      K(RIC02) = 1.9400000000D+57 &
	   * exp(-11.84 * LT - 414130000 / RT)
      K(RIC03) = 5.1500000000D-02 &
	   * exp(3.92 * LT - 10160000 / RT)
      K(RIC04) = 1.2500000000D+01 &
	   * exp(3.07 * LT - 5820000 / RT)
      K(RIC05) = 1.0300000000D+04 &
	   * exp(1.99 * LT + 1190000 / RT)
      K(RIC06) = 1.0300000000D+08 &
	   * exp(0.84 * LT - 196290000 / RT)
      K(RIC07) = 1.1400000000D-21 &
	   * exp(9.25 * LT + 8890000 / RT)
      K(RIC08) = 9.8500000000D+07 &
	   * exp(0.73 * LT - 70890000 / RT)
      K(RIC09) = 4.9300000000D+07 &
	   * exp(0.73 * LT - 70890000 / RT)
      K(RIC10) = 4.2800000000D+22 &
	   * exp(-2.81 * LT - 127700000 / RT)
      K(RIC11) = 2.5500000000D+39 &
	   * exp(-7.47 * LT - 189480000 / RT)
      K(RIC12) = 4.2200000000D+24 &
	   * exp(-3.34 * LT - 158660000 / RT)
      K(RIC13) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RIC14) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RIC15) = 1.3300000000D+23 &
	   * exp(-2.98 * LT - 64440000 / RT)
      K(RIC16) = 7.9500000000D+33 &
	   * exp(-6 * LT - 97570000 / RT)
      K(RIC17) = 2.6900000000D+20 &
	   * exp(-2.08 * LT - 62990000 / RT)
      K(RIC18) = 7.9200000000D+13 &
	   * exp(-1.01 * LT - 15980000 / RT)
      K(RIC19) = 1.8300000000D+17 &
	   * exp(-0.9 * LT - 165710000 / RT)
      K(RIC20) = 1.9400000000D+18 &
	   * exp(-1.49 * LT - 136770000 / RT)
      K(RIC21) = 1.1200000000D-44 &
	   * exp(15.73 * LT + 65590000 / RT)
      K(RIC22) = 6.5300000000D+59 &
	   * exp(-12.99 * LT - 394690000 / RT)
      K(RIC23) = 2.0900000000D+65 &
	   * exp(-14.94 * LT - 384520000 / RT)
      K(RIC24) = 2.3800000000D-16 &
	   * exp(7.67 * LT + 47650000 / RT)
      K(RIC25) = 3.9300000000D+03 &
	   * exp(2.36 * LT - 92290000 / RT)
      K(RIC26) = 7.7600000000D-12 &
	   * exp(6.18 * LT + 41330000 / RT)
      K(RIC27) = 1.9200000000D+02 &
	   * exp(2.76 * LT - 159330000 / RT)
      K(RIC28) = 1.0000000000D+10
      K(RIC29) = 1.7800000000D+39 &
	   * exp(-7.3 * LT - 156050000 / RT)
      K(RIC30) = 7.2500000000D+39 &
	   * exp(-7.59 * LT - 140030000 / RT)
      K(RIC31) = 3.0900000000D+13 &
	   * exp(0.03 * LT - 58870000 / RT)
      K(RIC32) = 1.3400000000D+71 &
	   * exp(-17.29 * LT - 263240000 / RT)
      K(RIC33) = 7.8400000000D+45 &
	   * exp(-9.64 * LT - 226420000 / RT)
      K(RIC34) = 2.1400000000D+38 &
	   * exp(-8.43 * LT - 60540000 / RT)
      K(RIC35) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RIC36) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RIC37) = 1.9200000000D+66 &
	   * exp(-14.22 * LT - 535970000 / RT)
      K(RIC38F) = 3.0700000000D+55 &
	   * exp(-11.49 * LT - 478230000 / RT)
      K(RIC38B) = 1.1210368967D+52 &
	   * exp(-11.8109 * LT - 110049848.5 / RT)
      K(RIC40) = 5.6800000000D+30 &
	   * exp(-5.72 * LT - 83680000 / RT)
      K(RIC41F) = 3.4000000000D+02 &
	   * exp(2.5 * LT - 10430000 / RT)
      K(RIC41B) = 2.9018835603D+01 &
	   * exp(2.52929 * LT - 75341336.28 / RT)
      K(RIC42) = 1.2100000000D+08 &
	   * exp(0.7 * LT - 31940000 / RT)
      K(RIC43) = 5.2000000000D+03 &
	   * exp(2 * LT + 1250000 / RT)
      K(RIC44) = 4.4200000000D-03 &
	   * exp(3.5 * LT - 23740000 / RT)
      K(RIC45) = 1.9300000000D+01 &
	   * exp(2.6 * LT - 58200000 / RT)
      K(RIC46) = 1.9300000000D+01 &
	   * exp(2.6 * LT - 58200000 / RT)
      K(RIC47) = 1.5800000000D+04 &
	   * exp(1.76 * LT + 5090000 / RT)
      K(RIC48) = 3.3300000000D+04 &
	   * exp(1.76 * LT - 320000 / RT)
      K(RIC49) = 1.0100000000D+18 &
	   * exp(-1.45 * LT - 129040000 / RT)
      K(RIC50) = 5.0000000000D+13 * exp(-121750000 / RT)
      K(RIC51) = 3.0000000000D+07 * exp(-6900000 / RT)
      K(RIC52) = 1.2300000000D+47 &
	   * exp(-9.74 * LT - 310700000 / RT)
      K(RIC57) = 3.1000000000D+10
      K(RIC53) = 6.0300000000D+10
      K(RIC54) = 2.4700000000D+10 &
	   * exp(-0.45 * LT - 96320000 / RT)
      K(RIC55) = 7.1400000000D+12 &
	   * exp(-1.21 * LT - 88070000 / RT)
      K(RIC56) = 7.2900000000D+26 &
	   * exp(-5.71 * LT - 89750000 / RT)
      K(RIC58) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RIC59) = 2.6000000000D+09 * exp(-10880000 / RT)
      K(RIC60) = 7.1800000000D+09 * exp(-5810000 / RT)
      K(RIC61) = 2.6900000000D+07 &
	   * exp(0.76 * LT + 1420000 / RT)
      K(RIC62) = 1.0000000000D+09 * exp(-49870000 / RT)
      K(RIC63) = 3.9800000000D+09 * exp(-36400000 / RT)
      K(RIC64) = 9.7700000000D-09 &
	   * exp(5.36 * LT - 71250000 / RT)
      K(RIC65) = 9.8800000000D+18 &
	   * exp(-1.59 * LT - 168820000 / RT)
      K(RIC66) = 1.7300000000D+10 &
	   * exp(0.03 * LT - 7520000 / RT)
      K(RIC67) = 4.1000000000D+08 * exp(-30140000 / RT)
      K(RIC68) = 7.6500000000D+08 &
	   * exp(-0.06 * LT - 21530000 / RT)
      K(RIC69) = 5.6300000000D+04 &
	   * exp(2 * LT - 32220000 / RT)
      K(RIC70) = 1.1300000000D+11 * exp(-32840000 / RT)
      K(RIC71) = 1.0500000000D+07 &
	   * exp(0.97 * LT - 6640000 / RT)
      K(RIC72) = 1.7000000000D+10 * exp(-85610000 / RT)
      K(RP000F) = 1.5000000000D+09 &
	   * exp(0.5 * LT - 8370000 / RT)
      K(RP000B) = 6.1199043605D+09 &
	   * exp(0.794062 * LT - 150584803.3 / RT)
      K(RP001F) = 4.6200000000D+12 &
	   * exp(-0.89 * LT - 38250000 / RT)
      K(RP001B) = 2.4351103303D+22 &
	   * exp(-2.30141 * LT - 194865690 / RT)
      K(RP002F) = 6.8000000000D+21 &
	   * exp(-3.45 * LT - 85090000 / RT)
      K(RP002B) = 6.5000342420D+29 &
	   * exp(-4.51606 * LT - 195688911.4 / RT)
      K(RP003F) = 1.4500000000D+45 &
	   * exp(-8.9 * LT - 405860000 / RT)
      K(RP003B) = 5.9159075488D+45 &
	   * exp(-8.60594 * LT - 548074803.3 / RT)
      K(RP004F) = 2.2400000000D+68 &
	   * exp(-14.65 * LT - 596540000 / RT)
      K(RP004B) = 5.6381479842D+62 &
	   * exp(-14.1444 * LT - 262706144.6 / RT)
      K(RP005F) = 9.6000000000D+67 &
	   * exp(-17.77 * LT - 130960000 / RT)
      K(RP005B) = 3.7352568967D+81 &
	   * exp(-19.121 * LT - 575688821.6 / RT)
      K(RP006F) = 1.3800000000D+13 &
	   * exp(-1 * LT - 37240000 / RT)
      K(RP006B) = 2.9676257451D+23 &
	   * exp(-2.11735 * LT - 336070493.3 / RT)
      K(RP007F) = 1.6700000000D+20 &
	   * exp(-3.3 * LT - 104430000 / RT)
      K(RP007B) = 6.5129325290D+28 &
	   * exp(-4.07199 * LT - 357243714.7 / RT)
      K(RP008) = 2.1600000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RP009F) = 8.2500000000D+43 &
	   * exp(-10.1 * LT - 70960000 / RT)
      K(RP009B) = 2.9151378138D+59 &
	   * exp(-11.8517 * LT - 535008868.8 / RT)
      K(RP010F) = 1.0700000000D+42 &
	   * exp(-9.57 * LT - 71190000 / RT)
      K(RP010B) = 1.5425608195D+58 &
	   * exp(-11.0277 * LT - 677453672.1 / RT)
      K(RP011F) = 5.7700000000D+34 &
	   * exp(-7 * LT - 131820000 / RT)
      K(RP011B) = 5.1317962822D+44 &
	   * exp(-8.24617 * LT - 262035013.4 / RT)
      K(RK012F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RK012B) = 2.8151364037D+14 &
	   * exp(0.27021 * LT - 191386104.5 / RT)
      K(RP013F) = 6.0000000000D+09
      K(RP013B) = 3.8029710467D+23 &
	   * exp(-2.06437 * LT - 497693591.2 / RT)
      K(RP014F) = 1.8700000000D+04 &
	   * exp(1.47 * LT - 23150000 / RT)
      K(RP014B) = 7.3122065614D+11 &
	   * exp(-0.382874 * LT - 44794932.53 / RT)
      K(RP015F) = 9.4500000000D-06 &
	   * exp(4.47 * LT - 18710000 / RT)
      K(RP015B) = 2.7559045191D-05 &
	   * exp(4.39414 * LT - 28239401.86 / RT)
      K(RP016F) = 4.3000000000D+60 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(RP016B) = 6.4261715627D+54 &
	   * exp(-12.2493 * LT - 132912595.6 / RT)
      K(RK017F) = 1.2900000000D+05 &
	   * exp(1.89 * LT - 73550000 / RT)
      K(RK017B) = 4.5060198908D+01 &
	   * exp(2.47097 * LT - 19964083.4 / RT)
      K(RP018F) = 7.8000000000D+00 &
	   * exp(2.68 * LT - 3070000 / RT)
      K(RP018B) = 2.9974411044D-02 &
	   * exp(3.15722 * LT - 12528608.14 / RT)
      K(RK019F) = 5.9000000000D+10 &
	   * exp(0.54 * LT - 115340000 / RT)
      K(RK019B) = 1.3748876245D+10 &
	   * exp(0.647444 * LT - 97958978.78 / RT)
      K(RK020F) = 7.1800000000D+10 &
	   * exp(1.02 * LT - 161810000 / RT)
      K(RK020B) = 4.1382935254D+06 &
	   * exp(1.32078 * LT - 30384498.09 / RT)
      K(RK021F) = 1.6500000000D+08 &
	   * exp(0.49 * LT - 44480000 / RT)
      K(RK021B) = 2.2227996788D+06 &
	   * exp(1.14102 * LT - 346145985.9 / RT)
      K(RP022F) = 2.5000000000D+09
      K(RP022B) = 3.7051777850D+08 &
	   * exp(0.547273 * LT - 364710510.6 / RT)
      K(RP023F) = 4.3000000000D+60 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(RP023B) = 6.1048303870D+54 &
	   * exp(-12.2599 * LT - 135870254.3 / RT)
      K(RK024F) = 5.2300000000D+03 &
	   * exp(2.36 * LT - 70780000 / RT)
      K(RK024B) = 1.7355069987D+00 &
	   * exp(2.93038 * LT - 20151742.05 / RT)
      K(RP025F) = 1.3400000000D-01 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP025B) = 4.8919514761D-04 &
	   * exp(3.79663 * LT - 18506266.79 / RT)
      K(RP026F) = 3.0100000000D+14 &
	   * exp(0.34 * LT - 465490000 / RT)
      K(RP026B) = 1.8338189281D+09 &
	   * exp(0.452698 * LT + 848724.5271 / RT)
      K(RP027F) = 6.3500000000D+01 &
	   * exp(2.75 * LT - 48740000 / RT)
      K(RP027B) = 9.0423895896D-02 &
	   * exp(3.21293 * LT - 15492763.27 / RT)
      K(RP028F) = 6.5500000000D-05 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RP028B) = 1.0261325190D-06 &
	   * exp(4.55919 * LT - 26197288.01 / RT)
      K(RK100F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RK100B) = 3.3117446183D+14 &
	   * exp(0.979368 * LT - 386582470.9 / RT)
      K(RK102F) = 3.0200000000D+00 &
	   * exp(2.55 * LT - 13310000 / RT)
      K(RK102B) = 3.6178371167D+10 &
	   * exp(1.2202 * LT - 292005478.6 / RT)
      K(RP104F) = 3.8000000000D+04 &
	   * exp(1.62 * LT - 18570000 / RT)
      K(RP104B) = 1.0608179381D+14 &
	   * exp(0.397642 * LT - 279884457.4 / RT)
      K(RP105F) = 3.6000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP105B) = 7.8690280050D+25 &
	   * exp(-3.13503 * LT - 311868721.5 / RT)
      K(RP106F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP106B) = 9.5260034588D+36 &
	   * exp(-6.03013 * LT - 365946869.1 / RT)
      K(RP107F) = 1.2600000000D+01 &
	   * exp(2.61 * LT - 6000000 / RT)
      K(RP107B) = 7.6077074478D+12 &
	   * exp(0.803868 * LT - 279683497.7 / RT)
      K(RP108F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP108B) = 6.5580381913D+54 &
	   * exp(-12.2468 * LT - 121345107.2 / RT)
      K(RK109F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK109B) = 4.7232496168D+01 &
	   * exp(2.45348 * LT - 6416595.001 / RT)
      K(RK110F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK110B) = 1.8883129559D-03 &
	   * exp(3.49973 * LT - 16231119.74 / RT)
      K(RP111F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP111B) = 7.4400104578D+54 &
	   * exp(-12.2625 * LT - 122636363.9 / RT)
      K(RK112F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK112B) = 5.3584662850D+01 &
	   * exp(2.43769 * LT - 7707851.712 / RT)
      K(RK113F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK113B) = 2.1422668990D-03 &
	   * exp(3.48394 * LT - 17522376.45 / RT)
      K(RK114F) = 1.2800000000D+03 &
	   * exp(2.05 * LT - 8080000 / RT)
      K(RK114B) = 2.2606983071D+14 &
	   * exp(0.192174 * LT - 181300070.2 / RT)
      K(RK115F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RK115B) = 2.8326023455D+14 &
	   * exp(0.287919 * LT - 186404982.5 / RT)
      K(RP116F) = 6.0000000000D+09
      K(RP116B) = 4.7823578389D+18 &
	   * exp(-2.02971 * LT - 26418832.42 / RT)
      K(RP117F) = 6.0000000000D+09
      K(RP117B) = 2.3313052517D+18 &
	   * exp(-1.93396 * LT - 26373744.7 / RT)
      K(RP118F) = 7.2000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP118B) = 1.0228650489D+20 &
	   * exp(-2.88623 * LT - 27235427.42 / RT)
      K(RP119F) = 7.2000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP119B) = 5.6568546074D+19 &
	   * exp(-2.80628 * LT - 28481596.41 / RT)
      K(RP120F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP120B) = 1.2133507349D+31 &
	   * exp(-5.78384 * LT - 92841063.38 / RT)
      K(RP121F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP121B) = 5.9148458473D+30 &
	   * exp(-5.6881 * LT - 92795975.65 / RT)
      K(RK122F) = 7.1800000000D+10 &
	   * exp(1.02 * LT - 161810000 / RT)
      K(RK122B) = 2.5916728182D+06 &
	   * exp(1.38384 * LT - 60849308.23 / RT)
      K(RK123F) = 1.6500000000D+08 &
	   * exp(0.49 * LT - 44480000 / RT)
      K(RK123B) = 1.3920640168D+06 &
	   * exp(1.20408 * LT - 376610796 / RT)
      K(RP124F) = 2.5000000000D+09
      K(RP124B) = 2.3204271261D+08 &
	   * exp(0.610331 * LT - 395175320.8 / RT)
      K(RK125F) = 7.1800000000D+10 &
	   * exp(1.02 * LT - 161810000 / RT)
      K(RK125B) = 4.2776375363D+06 &
	   * exp(1.29787 * LT - 58089370.29 / RT)
      K(RK126F) = 1.6500000000D+08 &
	   * exp(0.49 * LT - 44480000 / RT)
      K(RK126B) = 2.2976454626D+06 &
	   * exp(1.11811 * LT - 373850858.1 / RT)
      K(RP127F) = 2.5000000000D+09
      K(RP127B) = 3.8299379864D+08 &
	   * exp(0.524359 * LT - 392415382.8 / RT)
      K(RP128F) = 2.1500000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP128B) = 5.8671354650D+54 &
	   * exp(-12.2367 * LT - 111793778.6 / RT)
      K(RK129F) = 1.3200000000D+05 &
	   * exp(1.88 * LT - 70380000 / RT)
      K(RK129B) = 8.4193995682D+01 &
	   * exp(2.47349 * LT + 4284733.639 / RT)
      K(RP130F) = 6.7200000000D-02 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP130B) = 4.7155149309D-04 &
	   * exp(3.81974 * LT + 5530208.899 / RT)
      K(RP131F) = 2.1500000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP131B) = 5.1217965796D+54 &
	   * exp(-12.225 * LT - 110343473.6 / RT)
      K(RK132F) = 1.3200000000D+05 &
	   * exp(1.88 * LT - 70380000 / RT)
      K(RK132B) = 7.3498306231D+01 &
	   * exp(2.48528 * LT + 5735038.634 / RT)
      K(RP133F) = 6.7200000000D-02 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP133B) = 4.1164735991D-04 &
	   * exp(3.83153 * LT + 6980513.894 / RT)
      K(RK200F) = 2.8800000000D+11 &
	   * exp(0.23 * LT - 71240000 / RT)
      K(RK200B) = 2.1196844222D+08 &
	   * exp(0.811179 * LT - 85000673.5 / RT)
      K(RP201F) = 3.5200000000D+09 * exp(-55900000 / RT)
      K(RP201B) = 7.1773752755D+10 &
	   * exp(0.217336 * LT - 170621365.3 / RT)
      K(RP202F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP202B) = 1.2393355091D+55 &
	   * exp(-12.2259 * LT - 125051622.4 / RT)
      K(RK203F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK203B) = 8.9259787720D+01 &
	   * exp(2.47429 * LT - 10123110.2 / RT)
      K(RK204F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK204B) = 3.5685264864D-03 &
	   * exp(3.52054 * LT - 19937634.94 / RT)
      K(RK205F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RK205B) = 2.4444249624D+16 &
	   * exp(-0.322738 * LT - 188836244.6 / RT)
      K(RP206F) = 6.0000000000D+09
      K(RP206B) = 2.0118251902D+20 &
	   * exp(-2.54462 * LT - 28805006.84 / RT)
      K(RP207F) = 3.6000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP207B) = 4.0658510808D+21 &
	   * exp(-3.38034 * LT - 33328117.05 / RT)
      K(RP208F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP208B) = 5.1042804724D+32 &
	   * exp(-6.29875 * LT - 95227237.8 / RT)
      K(RP209F) = 2.1500000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP209B) = 5.6385261266D+54 &
	   * exp(-12.2371 * LT - 118341439.3 / RT)
      K(RK210F) = 1.3200000000D+05 &
	   * exp(1.88 * LT - 70380000 / RT)
      K(RK210B) = 8.0913428228D+01 &
	   * exp(2.4731 * LT - 2262927.138 / RT)
      K(RP211F) = 6.7200000000D-02 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP211B) = 4.5317777810D-04 &
	   * exp(3.81935 * LT - 1017451.878 / RT)
      K(RK212F) = 7.1800000000D+10 &
	   * exp(1.02 * LT - 161810000 / RT)
      K(RK212B) = 6.0688599637D+04 &
	   * exp(1.88558 * LT - 52706961.72 / RT)
      K(RK213F) = 1.6500000000D+08 &
	   * exp(0.49 * LT - 44480000 / RT)
      K(RK213B) = 3.2597639329D+04 &
	   * exp(1.70581 * LT - 368468449.5 / RT)
      K(RP214F) = 2.5000000000D+09
      K(RP214B) = 5.4336902348D+06 &
	   * exp(1.11207 * LT - 387032974.3 / RT)
      K(RP301F) = 9.5500000000D+08 * exp(-9070000 / RT)
      K(RP301B) = 1.0406159274D+16 &
	   * exp(-1.57654 * LT - 32761800.62 / RT)
      K(RP302F) = 1.3900000000D+10 * exp(-470000 / RT)
      K(RP302B) = 2.4550866312D+23 &
	   * exp(-1.78804 * LT - 500210459.3 / RT)
      K(RP304F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP304B) = 8.0262122096D+55 &
	   * exp(-12.4557 * LT - 126936787.3 / RT)
      K(RP305F) = 4.0100000000D+05 &
	   * exp(1.8 * LT - 68420000 / RT)
      K(RP305B) = 8.7473417297D+02 &
	   * exp(2.17451 * LT - 8898275.049 / RT)
      K(RP306F) = 2.6900000000D-01 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP306B) = 6.4556007352D-03 &
	   * exp(3.60077 * LT - 9612799.79 / RT)
      K(RK401F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RK401B) = 3.2105425876D+14 &
	   * exp(0.97143 * LT - 409641149.6 / RT)
      K(RK403F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RK403B) = 4.0290414867D+14 &
	   * exp(0.96566 * LT - 412605223.6 / RT)
      K(RP405F) = 3.8000000000D+04 &
	   * exp(1.62 * LT - 18570000 / RT)
      K(RP405B) = 4.7450169757D+13 &
	   * exp(0.473091 * LT - 323503486.8 / RT)
      K(RP406F) = 3.8000000000D+04 &
	   * exp(1.62 * LT - 18570000 / RT)
      K(RP406B) = 8.5798776135D+13 &
	   * exp(0.393137 * LT - 322257317.8 / RT)
      K(RP407F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RP407B) = 1.9481403928D+13 &
	   * exp(0.444306 * LT - 286489670.2 / RT)
      K(RP408F) = 6.0000000000D+09
      K(RP408B) = 3.4325764465D+26 &
	   * exp(-1.92589 * LT - 766849162.3 / RT)
      K(RP409F) = 6.0000000000D+09
      K(RP409B) = 4.3076808769D+26 &
	   * exp(-1.93166 * LT - 769813236.2 / RT)
      K(RP410F) = 1.8000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP410B) = 2.8101475820D+25 &
	   * exp(-3.12264 * LT - 325022940.8 / RT)
      K(RP411F) = 1.8000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP411B) = 3.0785672344D+25 &
	   * exp(-3.11662 * LT - 326536709.8 / RT)
      K(RP412F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP412B) = 3.7260182077D+36 &
	   * exp(-6.03025 * LT - 400179905.4 / RT)
      K(RP413F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP413B) = 4.6759329706D+36 &
	   * exp(-6.03603 * LT - 403143979.4 / RT)
      K(RP414F) = 1.2600000000D+01 &
	   * exp(2.61 * LT - 6000000 / RT)
      K(RP414B) = 7.0239241606D+12 &
	   * exp(0.801281 * LT - 318366492.8 / RT)
      K(RP415F) = 1.2600000000D+01 &
	   * exp(2.61 * LT - 6000000 / RT)
      K(RP415B) = 6.1912766332D+12 &
	   * exp(0.817072 * LT - 317075236.1 / RT)
      K(RP416F) = 3.1800000000D+08 * exp(-9070000 / RT)
      K(RP416B) = 6.2938340689D+19 &
	   * exp(-1.89045 * LT - 242726314.2 / RT)
      K(RP417F) = 2.3900000000D+08 * exp(-9070000 / RT)
      K(RP417B) = 1.9527079107D+19 &
	   * exp(-1.90968 * LT - 253355059.9 / RT)
      K(RP418F) = 1.3900000000D+10 * exp(-470000 / RT)
      K(RP418B) = 1.8408529969D+27 &
	   * exp(-2.12118 * LT - 720803718.6 / RT)
      K(RP419F) = 4.3000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP419B) = 8.1270195379D+54 &
	   * exp(-12.2546 * LT - 110170749.6 / RT)
      K(RK420F) = 2.0000000000D+05 &
	   * exp(1.8 * LT - 68420000 / RT)
      K(RK420B) = 8.8351184781D+01 &
	   * exp(2.37567 * LT + 7867762.628 / RT)
      K(RK421F) = 1.3400000000D-01 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RK421B) = 6.5123816232D-04 &
	   * exp(3.80192 * LT + 7153237.887 / RT)
      K(RP422F) = 4.3000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP422B) = 8.1270195379D+54 &
	   * exp(-12.2546 * LT - 110170749.6 / RT)
      K(RK423F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK423B) = 1.1706531984D+02 &
	   * exp(2.44567 * LT + 4757762.628 / RT)
      K(RK424F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK424B) = 4.6801667934D-03 &
	   * exp(3.49192 * LT - 5056762.113 / RT)
      K(RP425F) = 1.3000000000D+11 &
	   * exp(1.08 * LT - 294550000 / RT)
      K(RP425B) = 5.8424565666D-02 &
	   * exp(2.8367 * LT - 18258059.5 / RT)
      K(RP501F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RP501B) = 3.5811034647D+14 &
	   * exp(0.969451 * LT - 413991985.3 / RT)
      K(RP502F) = 3.8000000000D+04 &
	   * exp(1.62 * LT - 18570000 / RT)
      K(RP502B) = 6.8087951892D+11 &
	   * exp(0.969942 * LT - 321562695.3 / RT)
      K(RP503F) = 1.2800000000D+03 &
	   * exp(2.05 * LT - 8080000 / RT)
      K(RP503B) = 7.6582509943D+10 &
	   * exp(0.783907 * LT - 210386999.4 / RT)
      K(RP504F) = 6.0000000000D+09
      K(RP504B) = 3.7865926687D+17 &
	   * exp(-1.08774 * LT - 488597249.4 / RT)
      K(RP505F) = 1.2600000000D+01 &
	   * exp(2.61 * LT - 6000000 / RT)
      K(RP505B) = 4.2399475729D+12 &
	   * exp(0.78322 * LT - 318811875.7 / RT)
      K(RP506F) = 2.1500000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP506B) = 7.1084295844D+54 &
	   * exp(-12.252 * LT - 114867691.1 / RT)
      K(RK507F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK507B) = 2.0478616520D+02 &
	   * exp(2.44819 * LT + 60821.12064 / RT)
      K(RK508F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK508B) = 8.1871677406D-03 &
	   * exp(3.49444 * LT - 9753703.62 / RT)
      K(RK600F) = 1.2800000000D+03 &
	   * exp(2.05 * LT - 8080000 / RT)
      K(RK600B) = 4.3267077058D+11 &
	   * exp(0.882089 * LT - 289585820.2 / RT)
      K(RP601F) = 1.7200000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP601B) = 1.2436066718D+54 &
	   * exp(-12.2341 * LT - 85281739.28 / RT)
      K(RK602F) = 5.3000000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK602B) = 8.9567406666D+01 &
	   * exp(2.46617 * LT + 29646772.93 / RT)
      K(RK603F) = 1.9300000000D+00 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK603B) = 3.5882616284D-03 &
	   * exp(3.51242 * LT + 19832248.19 / RT)
      K(RK700F) = 1.2800000000D+03 &
	   * exp(2.05 * LT - 8080000 / RT)
      K(RK700B) = 8.7943130143D+10 &
	   * exp(0.76659 * LT - 230159310.1 / RT)
      K(RK701F) = 1.2800000000D+03 &
	   * exp(2.05 * LT - 8080000 / RT)
      K(RK701B) = 1.0865488865D+11 &
	   * exp(0.882751 * LT - 279772179 / RT)
      K(RP800) = 6.3700000000D+08 * exp(-9070000 / RT)
      K(RP801) = 9.5500000000D+08 * exp(-9070000 / RT)
      K(RP802) = 1.3900000000D+10 * exp(-470000 / RT)
      K(RCP01F) = 1.7300000000D+68 &
	   * exp(-15.16 * LT - 486900000 / RT)
      K(RCP01B) = 7.5313864907D+61 &
	   * exp(-14.464 * LT - 137417681 / RT)
      K(RCP02F) = 2.8000000000D+10 * exp(-9450000 / RT)
      K(RCP02B) = 2.8490915000D+06 &
	   * exp(1.04627 * LT - 93059168.76 / RT)
      K(RCP03F) = 3.3000000000D+11 * exp(-51650000 / RT)
      K(RCP03B) = 1.1779783274D+04 &
	   * exp(1.13526 * LT - 2998222.535 / RT)
      K(RCP04) = 3.3000000000D+11 * exp(-51650000 / RT)
      K(RCP05F) = 4.7700000000D+01 &
	   * exp(2.71 * LT - 4630000 / RT)
      K(RCP05B) = 2.9097834187D-03 &
	   * exp(3.72012 * LT - 82149751.19 / RT)
      K(RCP06F) = 3.0800000000D+03 * exp(2 * LT)
      K(RCP06B) = 3.4478768376D+00 &
	   * exp(2.94252 * LT - 146653693.5 / RT)
      K(RCP07F) = 1.0000000000D+11 * exp(-155440000 / RT)
      K(RCP07B) = 6.2711180531D+07 &
	   * exp(0.684196 * LT - 10625496.27 / RT)
      K(RCP08F) = 1.1000000000D+01 &
	   * exp(2.6 * LT - 53970000 / RT)
      K(RCP08B) = 3.3005564203D-01 &
	   * exp(3.03931 * LT - 70740025.27 / RT)
      K(RCP09F) = 1.8000000000D-04 * exp(4 * LT)
      K(RCP09B) = 3.4048615419D-05 &
	   * exp(4.52606 * LT - 92521076.66 / RT)
      K(RCP10F) = 6.0000000000D+09
      K(RCP10B) = 1.4518195750D+09 &
	   * exp(0.560401 * LT - 117036937.8 / RT)
      K(RCP11F) = 6.0000000000D+09
      K(RCP11B) = 1.5162412910D+12 &
	   * exp(-0.2405 * LT - 123947539.5 / RT)
      K(RCP12F) = 7.6600000000D+06 &
	   * exp(0.88 * LT - 4770000 / RT)
      K(RCP12B) = 2.1945908043D+10 &
	   * exp(-0.123996 * LT - 111472570.2 / RT)
      K(RCP13) = 3.8900000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RCP14) = 3.7500000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RCP15) = 3.7500000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RCP16F) = 6.8700000000D+52 &
	   * exp(-12.5 * LT - 176000000 / RT)
      K(RCP16B) = 5.5097990876D+64 &
	   * exp(-13.5715 * LT - 487565243.7 / RT)
      K(RCP17) = 6.3900000000D+26 &
	   * exp(-4.03 * LT - 147300000 / RT)
      K(RCP18) = 4.9100000000D+28 &
	   * exp(-4.85 * LT - 103650000 / RT)
      K(RCP19F) = 7.0000000000D+10
      K(RCP19B) = 1.2141817939D+15 &
	   * exp(-1.03629 * LT - 240921557.8 / RT)
      K(RCP20F) = 4.3400000000D+04 &
	   * exp(1.3 * LT - 73920000 / RT)
      K(RCP20B) = 1.5046362505D+06 &
	   * exp(0.699415 * LT - 244021162 / RT)
      K(RCP21F) = 3.1000000000D+10
      K(RCP21B) = 1.5338425964D+13 &
	   * exp(-0.692564 * LT - 18361752.59 / RT)
      K(RCP22) = 1.0200000000D+10
      K(RCP23F) = 2.0000000000D+13 * exp(-125520000 / RT)
      K(RCP23B) = 9.7282907499D+08 &
	   * exp(0.103818 * LT - 72591594.15 / RT)
      K(RCP24) = 1.0000000000D+12 * exp(-150620000 / RT)
      K(RCP25) = 6.6700000000D+09
      K(RCP26) = 3.3300000000D+09
      K(RCP27) = 1.6700000000D+09
      K(RCP28) = 1.4300000000D-16 &
	   * exp(7.6 * LT + 14770000 / RT)
      K(RCP29) = 1.3300000000D+09
      K(RCP30) = 1.2800000000D+04 &
	   * exp(1.02 * LT + 8510000 / RT)
      K(RCP31) = 3.3700000000D+44 &
	   * exp(-8 * LT - 454700000 / RT)
      K(RCP32F) = 2.7400000000D+06 &
	   * exp(1.46 * LT - 5670000 / RT)
      K(RCP32B) = 1.0395869997D+12 &
	   * exp(0.796263 * LT - 220933331.5 / RT)
      K(RCP33F) = 1.0000000000D+10 * exp(-8370000 / RT)
      K(RCP33B) = 2.2229923910D+04 &
	   * exp(1.1293 * LT - 406220625.7 / RT)
      K(RCP34) = 7.8600000000D-04 &
	   * exp(3.07 * LT - 23970000 / RT)
      K(RI00F) = 1.0000000000D+10
      K(RI00B) = 5.2485026204D+25 &
	   * exp(-2.03247 * LT - 521300387.2 / RT)
      K(RI01F) = 1.7300000000D+68 &
	   * exp(-15.16 * LT - 486900000 / RT)
      K(RI01B) = 1.8203611359D+63 &
	   * exp(-14.9678 * LT - 147515954.9 / RT)
      K(RI02F) = 2.8000000000D+10 * exp(-9450000 / RT)
      K(RI02B) = 6.8863488094D+07 &
	   * exp(0.542404 * LT - 103157442.7 / RT)
      K(RI03F) = 3.1600000000D+01 &
	   * exp(2.48 * LT - 46280000 / RT)
      K(RI03B) = 4.5268324245D+08 &
	   * exp(1.36023 * LT - 103588858.9 / RT)
      K(RI05F) = 4.7700000000D+01 &
	   * exp(2.71 * LT - 4630000 / RT)
      K(RI05B) = 7.0330431931D-02 &
	   * exp(3.21625 * LT - 92248025.11 / RT)
      K(RI06F) = 3.0800000000D+03 * exp(2 * LT)
      K(RI06B) = 8.3336328636D+01 &
	   * exp(2.43866 * LT - 156751967.4 / RT)
      K(RI07F) = 1.0000000000D+11 * exp(-155440000 / RT)
      K(RI07B) = 1.5157500677D+09 &
	   * exp(0.180331 * LT - 20723770.18 / RT)
      K(RI08F) = 1.1000000000D+01 &
	   * exp(2.6 * LT - 53970000 / RT)
      K(RI08B) = 7.9775545205D+00 &
	   * exp(2.53544 * LT - 80838299.18 / RT)
      K(RI09F) = 1.8000000000D-04 * exp(4 * LT)
      K(RI09B) = 8.2296634648D-04 &
	   * exp(4.0222 * LT - 102619350.6 / RT)
      K(RI12) = 3.8300000000D+06 &
	   * exp(0.88 * LT - 4770000 / RT)
      K(RI15) = 1.8800000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RI17) = 2.5600000000D+26 &
	   * exp(-4.03 * LT - 147300000 / RT)
      K(RI18) = 1.9600000000D+28 &
	   * exp(-4.85 * LT - 103650000 / RT)
      K(RI19F) = 2.8000000000D+10
      K(RI19B) = 2.1712467711D+13 &
	   * exp(-0.525133 * LT - 246624332.8 / RT)
      K(RI20F) = 1.7400000000D+04 &
	   * exp(1.3 * LT - 73920000 / RT)
      K(RI20B) = 2.6968482979D+04 &
	   * exp(1.21058 * LT - 249723937 / RT)
      K(RI21) = 1.2400000000D+10
      K(RI22) = 4.0800000000D+09
      K(RI23) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RI25) = 2.3800000000D+15 * exp(-377480000 / RT)
      K(RI26) = 1.1900000000D+15 * exp(-377480000 / RT)
      K(RI31) = 3.3700000000D+44 &
	   * exp(-8 * LT - 454700000 / RT)
      K(RI32) = 1.3700000000D+06 &
	   * exp(1.46 * LT - 5670000 / RT)
      K(RT01F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RT01B) = 4.5426319227D-03 &
	   * exp(3.57723 * LT - 50953235.53 / RT)
      K(RT02F) = 1.5600000000D+13 &
	   * exp(0.68 * LT - 373240000 / RT)
      K(RT02B) = 1.2813196505D+10 &
	   * exp(0.254459 * LT + 1780318.268 / RT)
      K(RT03F) = 4.3500000000D+22 &
	   * exp(-1.73 * LT - 436010000 / RT)
      K(RT03B) = 5.2773986505D+10 &
	   * exp(-0.111273 * LT + 6505423.169 / RT)
      K(RT04F) = 5.8300000000D+64 &
	   * exp(-14.15 * LT - 285890000 / RT)
      K(RT04B) = 8.6112522502D+55 &
	   * exp(-12.1057 * LT - 218394895.1 / RT)
      K(RT05F) = 8.2000000000D+14 * exp(-337550000 / RT)
      K(RT05B) = 4.3614434871D+01 &
	   * exp(2.2454 * LT - 46091406.02 / RT)
      K(RT06F) = 2.1800000000D+04 &
	   * exp(2.5 * LT - 192650000 / RT)
      K(RT06B) = 2.5793217860D+04 &
	   * exp(2.06262 * LT - 22297497.04 / RT)
      K(RT07F) = 6.4700000000D-03 &
	   * exp(3.98 * LT - 14160000 / RT)
      K(RT07B) = 1.2421004032D-03 &
	   * exp(3.90469 * LT - 72231169.53 / RT)
      K(RT08F) = 1.7700000000D+02 &
	   * exp(2.39 * LT + 2520000 / RT)
      K(RT08B) = 3.7383361476D+02 &
	   * exp(2.21095 * LT - 118595694.3 / RT)
      K(RT09F) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RT09B) = 6.8369388406D-01 &
	   * exp(2.76966 * LT - 48745764.6 / RT)
      K(RT10F) = 3.1400000000D-02 &
	   * exp(3.37 * LT - 19750000 / RT)
      K(RT10B) = 1.6980894524D+04 &
	   * exp(1.86503 * LT - 21919143.14 / RT)
      K(RT11F) = 1.1800000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RT11B) = 1.3580893683D-04 &
	   * exp(3.97854 * LT - 62631751.96 / RT)
      K(RT12F) = 1.6900000000D+09 &
	   * exp(0.3 * LT - 18420000 / RT)
      K(RT12B) = 1.2806518310D+17 &
	   * exp(-0.890885 * LT - 447591213.4 / RT)
      K(RT13F) = 1.6600000000D+04 &
	   * exp(1.8 * LT - 16630000 / RT)
      K(RT13B) = 1.9339064481D+07 &
	   * exp(0.788084 * LT - 79341310.47 / RT)
      K(RT14F) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RT14B) = 1.5060641068D+14 &
	   * exp(-0.595515 * LT - 160103077.4 / RT)
      K(RT15F) = 9.3300000000D+01 &
	   * exp(2.5 * LT - 61440000 / RT)
      K(RT15B) = 5.2817783980D+03 &
	   * exp(1.81773 * LT - 52672026.04 / RT)
      K(RT16F) = 7.9400000000D+10 * exp(-50000000 / RT)
      K(RT16B) = 1.0571054566D+14 &
	   * exp(-0.637038 * LT - 151028340.4 / RT)
      K(RT17F) = 2.2800000000D+11
      K(RT17B) = 8.0932136961D+17 &
	   * exp(-0.658482 * LT - 335373742.3 / RT)
      K(RT18F) = 2.0000000000D+10
      K(RT18B) = 4.8247335589D+17 &
	   * exp(-0.492206 * LT - 341715698 / RT)
      K(RT19F) = 1.1900000000D+06 &
	   * exp(1.03 * LT + 9410000 / RT)
      K(RT19B) = 5.8610185036D+06 &
	   * exp(0.819065 * LT - 50475531.31 / RT)
      K(RT20) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RT21F) = 3.7600000000D+12 &
	   * exp(-1.55 * LT - 47370000 / RT)
      K(RT21B) = 3.1887853062D+12 &
	   * exp(-1.51538 * LT - 253335838.4 / RT)
      K(RT22F) = 6.2600000000D+34 &
	   * exp(-8.86 * LT - 69380000 / RT)
      K(RT22B) = 2.6802796432D+32 &
	   * exp(-8.19307 * LT - 329695787.4 / RT)
      K0TROE = 4.6600000000D+35 &
	   * exp(-7.44 * LT - 58910000 / RT)
      KINFTROE = 2.4300000000D+09 &
	   * exp(0.52 * LT - 210000 / RT)
      FCTROE = 0.3 * EXP( -TEMP / 100 ) &
	   + 0.7 * EXP( -TEMP / 90000 ) &
	   + 1 * EXP( -10000 / TEMP )
      K(RT23F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 4.4376857912D+38 &
	   * exp(-6.95964 * LT - 492254025.9 / RT)
      KINFTROE = 2.3140722044D+12 &
	   * exp(1.00036 * LT - 433554025.9 / RT)
      K(RT23B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RT24F) = 4.2000000000D+03 &
	   * exp(2.1 * LT - 20380000 / RT)
      K(RT24B) = 1.0308564486D+03 &
	   * exp(1.96988 * LT - 20127461.91 / RT)
      K(RT25F) = 1.3000000000D+02 &
	   * exp(2.5 * LT - 20920000 / RT)
      K(RT25B) = 1.9128731618D+01 &
	   * exp(2.33372 * LT - 14578044.35 / RT)
      K(RT26F) = 6.3000000000D+03 &
	   * exp(2 * LT - 6280000 / RT)
      K(RT26B) = 1.7011480543D+04 &
	   * exp(1.76613 * LT - 69071986.65 / RT)
      K(RT27F) = 1.0000000000D+04 &
	   * exp(1.5 * LT - 41590000 / RT)
      K(RT27B) = 4.5627585431D+06 &
	   * exp(0.849667 * LT - 50249369.82 / RT)
      K(RT28F) = 8.5500000000D+01 &
	   * exp(2.19 * LT - 160000 / RT)
      K(RT28B) = 9.8745294140D+02 &
	   * exp(2.07772 * LT - 106603252.2 / RT)
      K(RT29F) = 5.2600000000D+28 &
	   * exp(-5.08 * LT - 93090000 / RT)
      K(RT29B) = 6.2875420673D+24 &
	   * exp(-4.82261 * LT - 34502491.8 / RT)
      K(RT30F) = 7.2100000000D+33 &
	   * exp(-6.21 * LT - 154180000 / RT)
      K(RT30B) = 4.1991617492D+20 &
	   * exp(-3.78987 * LT - 40040180.59 / RT)
      K(RT31F) = 2.3700000000D+32 &
	   * exp(-6.1 * LT - 120540000 / RT)
      K(RT31B) = 1.5941382341D+20 &
	   * exp(-3.79215 * LT - 112843432.7 / RT)
      K(RT32F) = 1.3300000000D+10
      K(RT32B) = 3.7159198196D+08 &
	   * exp(0.607624 * LT - 374503979.6 / RT)
      K(RT33F) = 6.6700000000D+09
      K(RT33B) = 1.1172090559D+08 &
	   * exp(0.571472 * LT - 368414562 / RT)
      K(RT34F) = 3.3300000000D+09
      K(RT34B) = 1.0235559534D+09 &
	   * exp(0.503877 * LT - 437548504.3 / RT)
      K(RT35F) = 2.8500000000D-16 &
	   * exp(7.6 * LT + 14770000 / RT)
      K(RT35B) = 4.9074531718D-17 &
	   * exp(7.84555 * LT - 131310307.1 / RT)
      K(RT36) = 2.1000000000D+16 * exp(-342000000 / RT)
      K(RT37) = 4.0900000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RT38) = 5.8400000000D+09 * exp(-7570000 / RT)
      K(RT39) = 2.8900000000D+05 &
	   * exp(1.35 * LT + 6580000 / RT)
      K(RT40) = 1.2000000000D+02 &
	   * exp(2.5 * LT - 157130000 / RT)
      K(RT41) = 4.0900000000D+01 &
	   * exp(2.5 * LT - 42690000 / RT)
      K(RT42) = 3.4900000000D-11 &
	   * exp(6.21 * LT - 6820000 / RT)
      K(RT43F) = 1.0100000000D+71 &
	   * exp(-15.92 * LT - 522120000 / RT)
      K(RT43B) = 1.5527601509D+66 &
	   * exp(-15.741 * LT - 155660097.1 / RT)
      K(RT44F) = 1.1500000000D+11 * exp(-51870000 / RT)
      K(RT44B) = 4.1323812035D+08 &
	   * exp(0.529204 * LT - 118501584.9 / RT)
      K(RT45F) = 2.8100000000D+10 * exp(-30760000 / RT)
      K(RT45B) = 6.0534474546D+07 &
	   * exp(0.493052 * LT - 91302167.33 / RT)
      K(RT46F) = 2.9500000000D+03 &
	   * exp(2 * LT + 5490000 / RT)
      K(RT46B) = 1.1662112183D+02 &
	   * exp(2.42546 * LT - 124186109.6 / RT)
      K(RT47) = 2.9000000000D+10 * exp(-152400000 / RT)
      K(RT50F) = 8.3900000000D+60 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(RT50B) = 5.8699915686D+54 &
	   * exp(-12.2592 * LT - 143869255.6 / RT)
      K(RT51F) = 3.9100000000D+05 &
	   * exp(1.8 * LT - 68420000 / RT)
      K(RT51B) = 6.3939878670D+01 &
	   * exp(2.37108 * LT - 25790743.41 / RT)
      K(RT52F) = 2.6000000000D+10 * exp(-61500000 / RT)
      K(RT52B) = 2.5489558805D+06 &
	   * exp(0.534927 * LT - 12781325.85 / RT)
      K(RT53F) = 1.3600000000D+01 &
	   * exp(2.69 * LT - 2590000 / RT)
      K(RT53B) = 2.4467332309D-02 &
	   * exp(3.15733 * LT - 23005268.15 / RT)
      K(RT54F) = 1.7900000000D-05 &
	   * exp(4.46 * LT - 57060000 / RT)
      K(RT54B) = 5.4416006531D-06 &
	   * exp(4.51087 * LT - 23342651.31 / RT)
      K(RT55F) = 1.0000000000D+11
      K(RT55B) = 1.6651442300D+20 &
	   * exp(-1.23276 * LT - 538432054.9 / RT)
      K(RT56F) = 3.0000000000D+10
      K(RT56B) = 3.5650010139D+17 &
	   * exp(-1.54684 * LT - 111429984.6 / RT)
      K(RT57F) = 3.0000000000D+10
      K(RT57B) = 6.9312742321D+13 &
	   * exp(-0.785212 * LT - 262943843.8 / RT)
      K(RT58F) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(RT58B) = 2.0355198398D+22 &
	   * exp(-3.38113 * LT - 70689588.9 / RT)
      K(RT59) = 2.5500000000D+10 &
	   * exp(-0.44 * LT + 6900000 / RT)
      K(RT60) = 2.5500000000D+10 &
	   * exp(-0.44 * LT + 6900000 / RT)
      K(RE01F) = 3.6100000000D+10
      K(RE01B) = 9.9215141237D+17 &
	   * exp(-0.723732 * LT - 430050869.8 / RT)
      K(RE02F) = 2.0000000000D+10
      K(RE02B) = 3.0096212144D+21 &
	   * exp(-1.4418 * LT - 329900840.3 / RT)
      K(RE03F) = 2.0000000000D+10
      K(RE03B) = 1.1999372607D+25 &
	   * exp(-2.29587 * LT - 434547961.2 / RT)
      K(RE04F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RE04B) = 6.2409163294D-06 &
	   * exp(4.25437 * LT - 58920697.52 / RT)
      K(RE05F) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RE05B) = 9.3929607283D-04 &
	   * exp(3.4468 * LT - 56713226.59 / RT)
      K(RE06F) = 4.8300000000D-05 &
	   * exp(4.71 * LT - 25990000 / RT)
      K(RE06B) = 4.1076734557D-10 &
	   * exp(5.78397 * LT - 29030618.02 / RT)
      K(RE07F) = 1.9600000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RE07B) = 9.9930663771D-09 &
	   * exp(5.12782 * LT - 7601200.453 / RT)
      K(RE08F) = 4.4700000000D+03 &
	   * exp(2.01 * LT - 1530000 / RT)
      K(RE08B) = 4.1822401369D-01 &
	   * exp(2.98022 * LT - 67615142.76 / RT)
      K(RE09F) = 8.0300000000D+00 &
	   * exp(2.6 * LT - 58200000 / RT)
      K(RE09B) = 2.0137726410D-02 &
	   * exp(3.06701 * LT + 5598525.469 / RT)
      K(RE10F) = 7.5300000000D-04 &
	   * exp(3.65 * LT - 29930000 / RT)
      K(RE10B) = 1.1904806541D-05 &
	   * exp(4.20376 * LT - 41882525.92 / RT)
      K(RE11F) = 1.7200000000D+11 &
	   * exp(0.78 * LT - 161940000 / RT)
      K(RE11B) = 4.2169346325D+00 &
	   * exp(2.47345 * LT - 6105276.279 / RT)
      K(RE12F) = 3.7900000000D+06 &
	   * exp(2.08 * LT - 134330000 / RT)
      K(RE12B) = 1.0596119199D+04 &
	   * exp(1.84471 * LT - 9669610.664 / RT)
      K(RE13F) = 1.6700000000D+09
      K(RE13B) = 1.0912986807D+09 &
	   * exp(0.114944 * LT - 308431098.5 / RT)
      K(RE14F) = 2.4100000000D+10
      K(RE14B) = 1.7325942283D+11 &
	   * exp(0.0111972 * LT - 371475623.2 / RT)
      K(RE15F) = 3.3100000000D+09 * exp(3220000 / RT)
      K(RE15B) = 4.0209968985D+12 &
	   * exp(-0.405265 * LT - 314123006.4 / RT)
      K(RE16F) = 3.7000000000D+13 &
	   * exp(-1.63 * LT - 14300000 / RT)
      K(RE16B) = 1.4901345826D+14 &
	   * exp(-1.87713 * LT - 94307425.97 / RT)
      K(RE17F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RE17B) = 8.1033372636D+08 &
	   * exp(0.465447 * LT - 387229057.3 / RT)
      K(RE18F) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RE18B) = 2.4565719683D+09 &
	   * exp(0.34697 * LT - 375286263.6 / RT)
      K(RE19) = 7.0000000000D+09 * exp(4180000 / RT)
      K(RE30) = 4.8200000000D+33 &
	   * exp(-8.23 * LT - 21620000 / RT)
      K(RE31) = 3.8800000000D+01 &
	   * exp(1.84 * LT + 2420000 / RT)
      K(RE32) = 1.2100000000D+45 &
	   * exp(-10.15 * LT - 170770000 / RT)
      K(RE33) = 2.4500000000D+25 &
	   * exp(-4.48 * LT - 136430000 / RT)
      K(RE34) = 1.2300000000D+13 &
	   * exp(-1.12 * LT - 113030000 / RT)
      K(RE35) = 5.4500000000D-02 &
	   * exp(3.57 * LT - 67350000 / RT)
      K(RE36) = 1.0000000000D+08
      K(RE37) = 6.3000000000D+14 * exp(-180100000 / RT)
      K(RST01) = 2.4000000000D+14 * exp(-326900000 / RT)
      K(RST02F) = 5.1000000000D+04 &
	   * exp(1.66 * LT - 2750000 / RT)
      K(RST02B) = 3.3224050070D+00 &
	   * exp(2.42996 * LT - 146684040.1 / RT)
      K(RST03F) = 1.1400000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RST03B) = 3.0178171755D+02 &
	   * exp(1.94272 * LT - 14154671.17 / RT)
      K(RST04F) = 1.8800000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RST04B) = 3.7188357721D+32 &
	   * exp(-7.56182 * LT - 104134582.7 / RT)
      K(RST05F) = 1.8800000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RST05B) = 1.2267086362D+32 &
	   * exp(-7.44334 * LT - 116077376.4 / RT)
      K(RST06F) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RST06B) = 8.8911775253D-03 &
	   * exp(3.2153 * LT + 6432403.456 / RT)
      K(RST10F) = 1.3300000000D+10 * exp(-61500000 / RT)
      K(RST10B) = 2.6458805967D+06 &
	   * exp(0.534226 * LT - 4782324.487 / RT)
      K(RST11F) = 1.0300000000D+10 &
	   * exp(0.21 * LT + 1790000 / RT)
      K(RST11B) = 2.0794236065D+07 &
	   * exp(1.09997 * LT - 542709192.8 / RT)
      K(RST12F) = 6.7000000000D+02 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RST12B) = 5.5627276119D+01 &
	   * exp(1.89895 * LT - 71632313.4 / RT)
      K(RST13) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RST14F) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RST14B) = 2.9640252599D+12 &
	   * exp(-1.01221 * LT - 393486532.8 / RT)
      K(RST00F) = 1.0000000000D+10
      K(RST00B) = 7.2340731681D+07 &
	   * exp(0.771407 * LT - 84130331.83 / RT)
      K(RXY00F) = 2.5000000000D+18 &
	   * exp(-0.6 * LT - 396590000 / RT)
      K(RXY00B) = 1.0443593278D+15 &
	   * exp(-1.02266 * LT - 21966408.23 / RT)
      K(RXY01F) = 4.3200000000D+29 &
	   * exp(-3.58 * LT - 460930000 / RT)
      K(RXY01B) = 2.4885817320D+17 &
	   * exp(-1.95607 * LT - 18728819.02 / RT)
      K(RXY02F) = 1.2900000000D-02 &
	   * exp(3.98 * LT - 14160000 / RT)
      K(RXY02B) = 1.2595609433D-03 &
	   * exp(3.90757 * LT - 72627896.03 / RT)
      K(RXY03F) = 2.3600000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RXY03B) = 1.3814507233D-04 &
	   * exp(3.98142 * LT - 63028478.46 / RT)
      K(RXY04F) = 3.5400000000D+02 &
	   * exp(2.39 * LT + 2520000 / RT)
      K(RXY04B) = 3.8026416345D+02 &
	   * exp(2.21383 * LT - 118992420.8 / RT)
      K(RXY05F) = 4.3600000000D+04 &
	   * exp(2.5 * LT - 192650000 / RT)
      K(RXY05B) = 2.6236903331D+04 &
	   * exp(2.0655 * LT - 22694223.54 / RT)
      K(RXY06F) = 1.8700000000D+02 &
	   * exp(2.5 * LT - 61440000 / RT)
      K(RXY06B) = 5.3841505118D+03 &
	   * exp(1.82061 * LT - 53068752.54 / RT)
      K(RXY07F) = 8.4400000000D+11 * exp(-93120000 / RT)
      K(RXY07B) = 1.5319708690D+14 &
	   * exp(-0.592634 * LT - 160499803.9 / RT)
      K(RXY09F) = 4.6200000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RXY09B) = 3.8039484327D-03 &
	   * exp(3.57309 * LT - 50939563.4 / RT)
      K(RXY10F) = 1.5700000000D+00 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RXY10B) = 6.9907410294D-01 &
	   * exp(2.77812 * LT - 49168706.55 / RT)
      K(RXY11) = 1.8200000000D+05 &
	   * exp(1.55 * LT - 12930000 / RT)
      K(RXY12) = 8.2000000000D+14 * exp(-337550000 / RT)
      K(RXY13F) = 5.8300000000D+64 &
	   * exp(-14.15 * LT - 285890000 / RT)
      K(RXY13B) = 8.0394580864D+55 &
	   * exp(-12.1034 * LT - 218312410.8 / RT)
      K(RXY14F) = 4.3700000000D+15 &
	   * exp(-1.34 * LT - 6660000 / RT)
      K(RXY14B) = 2.3683526004D+18 &
	   * exp(-1.74072 * LT - 283813818.9 / RT)
      K(RXY15F) = 5.9900000000D+20 &
	   * exp(-2.47 * LT - 67750000 / RT)
      K(RXY15B) = 1.1561143354D+14 &
	   * exp(-0.706024 * LT - 288901438.6 / RT)
      K(RXY16F) = 1.9700000000D+19 &
	   * exp(-2.36 * LT - 34110000 / RT)
      K(RXY16B) = 3.8721310084D+13 &
	   * exp(-0.71765 * LT - 361376776.5 / RT)
      K(RXY17) = 2.0000000000D+10
      K(RXY18F) = 1.3800000000D-01 &
	   * exp(2.42 * LT - 31130000 / RT)
      K(RXY18B) = 1.4948586364D-01 &
	   * exp(2.45499 * LT - 237463423.1 / RT)
      K(RXY19F) = 6.5700000000D+00 &
	   * exp(1.87 * LT - 20930000 / RT)
      K(RXY19B) = 3.0118535542D-02 &
	   * exp(2.52284 * LT - 282691027.5 / RT)
      K(RXY201) = 2.2800000000D+10 &
	   * exp(-0.31 * LT + 2750000 / RT)
      K(RXY202) = 3.1300000000D+15 &
	   * exp(-1.44 * LT - 58340000 / RT)
      K(RXY203) = 1.0300000000D+14 &
	   * exp(-1.33 * LT - 24700000 / RT)
      K(RXY22) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RXY23F) = 1.2500000000D+18 &
	   * exp(-0.6 * LT - 396590000 / RT)
      K(RXY23B) = 1.0588306890D+15 &
	   * exp(-1.02766 * LT - 19435800.39 / RT)
      K(RXY24) = 2.1600000000D+29 &
	   * exp(-3.58 * LT - 460930000 / RT)
      K(RXY25) = 2.1000000000D+16 * exp(-342000000 / RT)
      K(RXY26F) = 6.4700000000D-03 &
	   * exp(3.98 * LT - 14160000 / RT)
      K(RXY26B) = 1.2809740203D-03 &
	   * exp(3.90257 * LT - 70097288.19 / RT)
      K(RXY27F) = 1.1800000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RXY27B) = 1.4005930548D-04 &
	   * exp(3.97642 * LT - 60497870.63 / RT)
      K(RXY28F) = 1.7700000000D+02 &
	   * exp(2.39 * LT + 2520000 / RT)
      K(RXY28B) = 3.8553336526D+02 &
	   * exp(2.20882 * LT - 116461812.9 / RT)
      K(RXY29F) = 2.1800000000D+04 &
	   * exp(2.5 * LT - 192650000 / RT)
      K(RXY29B) = 2.6600459910D+04 &
	   * exp(2.0605 * LT - 20163615.7 / RT)
      K(RXY30F) = 9.3300000000D+01 &
	   * exp(2.5 * LT - 61440000 / RT)
      K(RXY30B) = 5.4470805191D+03 &
	   * exp(1.81561 * LT - 50538144.7 / RT)
      K(RXY31F) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RXY31B) = 1.5531989112D+14 &
	   * exp(-0.597639 * LT - 157969196.1 / RT)
      K(RXY33) = 4.0900000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RXY34) = 5.8400000000D+09 * exp(-7570000 / RT)
      K(RXY35) = 2.8900000000D+05 &
	   * exp(1.35 * LT + 6580000 / RT)
      K(RXY36) = 1.2000000000D+02 &
	   * exp(2.5 * LT - 157130000 / RT)
      K(RXY37) = 4.0900000000D+01 &
	   * exp(2.5 * LT - 42690000 / RT)
      K(RXY38) = 3.4900000000D-11 &
	   * exp(6.21 * LT - 6820000 / RT)
      K(RXY39F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RXY39B) = 8.3778170592D-06 &
	   * exp(4.21307 * LT - 67532957.61 / RT)
      K(RXY40F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RXY40B) = 2.9278172747D-03 &
	   * exp(3.56984 * LT - 50175252.2 / RT)
      K(RXY41F) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RXY41B) = 1.5357183416D-03 &
	   * exp(3.41811 * LT - 65762100.76 / RT)
      K(RXY42) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RXY43F) = 4.3700000000D+15 &
	   * exp(-1.34 * LT - 6660000 / RT)
      K(RXY43B) = 4.9681868542D+18 &
	   * exp(-1.77556 * LT - 280885485.3 / RT)
      K(RXY44) = 5.9900000000D+20 &
	   * exp(-2.47 * LT - 67750000 / RT)
      K(RXY45F) = 1.9700000000D+19 &
	   * exp(-2.36 * LT - 34110000 / RT)
      K(RXY45B) = 2.9395633250D+13 &
	   * exp(-0.715895 * LT - 363143073.1 / RT)
      K(RXY46) = 2.0000000000D+10
      K(RXY47F) = 1.3800000000D-01 &
	   * exp(2.42 * LT - 31980000 / RT)
      K(RXY47B) = 3.1358240429D-01 &
	   * exp(2.42015 * LT - 235385089.6 / RT)
      K(RXY48) = 1.1900000000D+06 &
	   * exp(1.03 * LT + 9410000 / RT)
      K(RXY50) = 4.2000000000D+16 * exp(-342000000 / RT)
      K(RXY51) = 8.1800000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RXY52) = 1.1700000000D+10 * exp(-7570000 / RT)
      K(RXY53) = 5.7800000000D+05 &
	   * exp(1.35 * LT + 6580000 / RT)
      K(RXY54) = 2.4000000000D+02 &
	   * exp(2.5 * LT - 157130000 / RT)
      K(RXY55) = 8.1800000000D+01 &
	   * exp(2.5 * LT - 42690000 / RT)
      K(RXY56) = 6.9800000000D-11 &
	   * exp(6.21 * LT - 6820000 / RT)
      K(RXY57F) = 4.6200000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RXY57B) = 6.0637611379D-06 &
	   * exp(4.24966 * LT - 72227587.75 / RT)
      K(RXY58) = 1.5700000000D+00 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RN01F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RN01B) = 2.6303314786D-05 &
	   * exp(4.14699 * LT - 50409540.47 / RT)
      K(RN02F) = 7.8300000000D-01 &
	   * exp(2.88 * LT - 13480000 / RT)
      K(RN02B) = 1.0967302857D-03 &
	   * exp(3.74052 * LT - 49731754.41 / RT)
      K(RN04F) = 1.2500000000D+18 &
	   * exp(-0.6 * LT - 396590000 / RT)
      K(RN04B) = 7.6429253495D+12 &
	   * exp(-0.435322 * LT - 9851879.559 / RT)
      K(RN05F) = 3.2000000000D+34 &
	   * exp(-5.02 * LT - 478030000 / RT)
      K(RN05B) = 2.7785868184D+20 &
	   * exp(-2.80977 * LT - 12814647.68 / RT)
      K(RN06F) = 5.8300000000D+64 &
	   * exp(-14.15 * LT - 285890000 / RT)
      K(RN06B) = 8.2792870995D+55 &
	   * exp(-12.1044 * LT - 207412768.1 / RT)
      K(RN07F) = 8.2000000000D+14 * exp(-337550000 / RT)
      K(RN07B) = 6.5241355800D+02 &
	   * exp(1.84259 * LT - 103227017.7 / RT)
      K(RN08F) = 6.4700000000D-03 &
	   * exp(3.98 * LT - 14160000 / RT)
      K(RN08B) = 9.2464157998D-06 &
	   * exp(4.49491 * LT - 60513367.35 / RT)
      K(RN09F) = 1.1800000000D-03 &
	   * exp(4.09 * LT - 10650000 / RT)
      K(RN09B) = 1.0109858237D-06 &
	   * exp(4.56876 * LT - 50913949.79 / RT)
      K(RN10F) = 1.7700000000D+02 &
	   * exp(2.39 * LT + 2520000 / RT)
      K(RN10B) = 2.7828837614D+00 &
	   * exp(2.80117 * LT - 106877892.1 / RT)
      K(RN11F) = 2.1800000000D+04 &
	   * exp(2.5 * LT - 192650000 / RT)
      K(RN11B) = 1.9200929050D+02 &
	   * exp(2.65284 * LT - 10579694.87 / RT)
      K(RN12F) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RN12B) = 1.1211408447D+12 &
	   * exp(-0.00529554 * LT - 148385275.3 / RT)
      K(RN13F) = 9.3300000000D+01 &
	   * exp(2.5 * LT - 61440000 / RT)
      K(RN13B) = 3.9318495595D+01 &
	   * exp(2.40795 * LT - 40954223.87 / RT)
      K(RN14) = 1.1000000000D+10 * exp(-18960000 / RT)
      K(RN15) = 1.4700000000D+10 * exp(-18960000 / RT)
      K(RN16F) = 2.2800000000D+11
      K(RN16B) = 1.0156529597D+18 &
	   * exp(-0.655564 * LT - 318345379.3 / RT)
      K(RN17) = 2.0000000000D+10
      K(RN18F) = 1.1900000000D+06 &
	   * exp(1.03 * LT + 9410000 / RT)
      K(RN18B) = 7.3552497358D+06 &
	   * exp(0.821983 * LT - 33447168.22 / RT)
      K(RN18) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RN20F) = 3.7600000000D+12 &
	   * exp(-1.55 * LT - 47370000 / RT)
      K(RN20B) = 1.7194981983D+14 &
	   * exp(-1.93634 * LT - 261259801.2 / RT)
      K(RN21F) = 3.0800000000D+06 &
	   * exp(0.37 * LT - 70750000 / RT)
      K(RN21B) = 1.9592341791D+03 &
	   * exp(1.34013 * LT - 347359225.5 / RT)
      K(RN22F) = 5.2600000000D+28 &
	   * exp(-5.08 * LT - 93090000 / RT)
      K(RN22B) = 2.7016744652D+26 &
	   * exp(-5.24648 * LT - 59454817.65 / RT)
      K(RN23F) = 7.2100000000D+33 &
	   * exp(-6.21 * LT - 154180000 / RT)
      K(RN23B) = 3.2171025483D+20 &
	   * exp(-3.7915 * LT - 46086416.7 / RT)
      K(RN24F) = 1.3300000000D+10
      K(RN24B) = 1.5966820713D+10 &
	   * exp(0.183753 * LT - 399456305.4 / RT)
      K(RN25F) = 6.6700000000D+09
      K(RN25B) = 4.8005009689D+09 &
	   * exp(0.147601 * LT - 393366887.9 / RT)
      K(RN26F) = 3.3300000000D+09
      K(RN26B) = 4.3980858551D+10 &
	   * exp(0.0800063 * LT - 462500830.2 / RT)
      K(RN27F) = 2.8500000000D-16 &
	   * exp(7.6 * LT + 14770000 / RT)
      K(RN27B) = 2.1086683447D-15 &
	   * exp(7.42168 * LT - 156262633 / RT)
      K(RN28) = 2.1000000000D+16 * exp(-342000000 / RT)
      K(RN29) = 4.0900000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RN30) = 5.8400000000D+09 * exp(-7570000 / RT)
      K(RN31) = 2.8900000000D+05 &
	   * exp(1.35 * LT + 6580000 / RT)
      K(RN32) = 1.2000000000D+02 &
	   * exp(2.5 * LT - 157130000 / RT)
      K(RN33) = 4.0900000000D+01 &
	   * exp(2.5 * LT - 42690000 / RT)
      K(RN34) = 3.4900000000D-11 &
	   * exp(6.21 * LT - 6820000 / RT)
      K(ROX00F) = 1.2900000000D+61 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(ROX00B) = 7.9583820932D+54 &
	   * exp(-12.2685 * LT - 143541341.3 / RT)
      K0TROE = 1.0000000000D+81 &
	   * exp(-18.87 * LT - 376980000 / RT)
      KINFTROE = 4.3000000000D+12 &
	   * exp(0.62 * LT - 323430000 / RT)
      FCTROE = 0.098 * EXP( -TEMP / 696 ) &
	   + 0.902 * EXP( -TEMP / 358 ) &
	   + 1 * EXP( -3856 / TEMP )
      K(ROX01F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 4.1465511326D+75 &
	   * exp(-18.4977 * LT - 50909087.01 / RT)
      KINFTROE = 1.7830169870D+07 &
	   * exp(0.992264 * LT + 2640912.989 / RT)
      K(ROX01B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(ROX02F) = 5.0000000000D+75 &
	   * exp(-19.31 * LT - 284180000 / RT)
      K(ROX02B) = 2.7802587455D+86 &
	   * exp(-20.3484 * LT - 537844017.3 / RT)
      K(ROX03F) = 6.0200000000D+05 &
	   * exp(1.8 * LT - 68420000 / RT)
      K(ROX03B) = 8.6806271172D+01 &
	   * exp(2.36173 * LT - 25462829.1 / RT)
      K(ROX04F) = 2.3400000000D+01 &
	   * exp(2.68 * LT - 3070000 / RT)
      K(ROX04B) = 3.7121295904D-02 &
	   * exp(3.13799 * LT - 23157353.84 / RT)
      K(ROX05F) = 1.3200000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX05B) = 5.8610887211D+04 &
	   * exp(1.73243 * LT - 25122529.07 / RT)
      K(ROX06F) = 4.3400000000D+11 * exp(-269000000 / RT)
      K(ROX06B) = 3.8569228355D+08 &
	   * exp(0.199659 * LT + 2380843.389 / RT)
      K(ROX99F) = 2.7400000000D+05 &
	   * exp(1.55 * LT - 12930000 / RT)
      K(ROX99B) = 2.4543429804D+08 &
	   * exp(0.545147 * LT - 73785671.75 / RT)
      K(ROX07F) = 1.9900000000D+04 &
	   * exp(1.8 * LT - 16630000 / RT)
      K(ROX07B) = 1.7825337704D+07 &
	   * exp(0.795147 * LT - 77485671.75 / RT)
      K(ROX08F) = 2.0300000000D+09 &
	   * exp(0.3 * LT - 18420000 / RT)
      K(ROX08B) = 1.2630305520D+17 &
	   * exp(-0.903488 * LT - 447154599.3 / RT)
      K(ROX09F) = 5.8000000000D+13 &
	   * exp(-0.77 * LT - 63980000 / RT)
      K(ROX09B) = 7.5519642443D+08 &
	   * exp(0.101044 * LT - 368610363.6 / RT)
      K(ROX10F) = 4.0000000000D+10 * exp(-61500000 / RT)
      K(ROX10B) = 3.4578688810D+06 &
	   * exp(0.525581 * LT - 12453411.54 / RT)
      K(ROX11F) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX11B) = 1.7748967608D+22 &
	   * exp(-3.36472 * LT - 69161864.49 / RT)
      K(ROX12F) = 3.0000000000D+10 * exp(-37580000 / RT)
      K(ROX12B) = 5.1597641764D+17 &
	   * exp(-1.69626 * LT - 249923385.2 / RT)
      K(ROX13F) = 1.0000000000D+11
      K(ROX13B) = 1.4519431559D+20 &
	   * exp(-1.21635 * LT - 536904330.4 / RT)
      K(ROX14F) = 3.0000000000D+10
      K(ROX14B) = 3.1085468331D+17 &
	   * exp(-1.53043 * LT - 109902260.2 / RT)
      K(ROX15F) = 3.0000000000D+10
      K(ROX15B) = 6.0438105008D+13 &
	   * exp(-0.768803 * LT - 261416119.4 / RT)
      K(ROX16F) = 3.8900000000D-06 &
	   * exp(4.57 * LT - 22000000 / RT)
      K(ROX16B) = 1.4511636996D-05 &
	   * exp(4.52848 * LT - 56045263 / RT)
      K(ROX17F) = 6.5900000000D+15 &
	   * exp(-0.61 * LT - 310110000 / RT)
      K(ROX17B) = 1.3791119702D+03 &
	   * exp(1.46453 * LT - 186005764.3 / RT)
      K(ROX18F) = 1.0100000000D+71 &
	   * exp(-15.92 * LT - 522120000 / RT)
      K(ROX18B) = 1.4540798147D+66 &
	   * exp(-15.7214 * LT - 154241072.4 / RT)
      K(ROX19F) = 6.8300000000D-02 &
	   * exp(3.4 * LT - 30260000 / RT)
      K(ROX19B) = 2.2983020763D-04 &
	   * exp(3.94887 * LT - 95472560.24 / RT)
      K(ROX20F) = 1.7300000000D-02 &
	   * exp(3.4 * LT + 4780000 / RT)
      K(ROX20B) = 6.4044991233D-04 &
	   * exp(3.84512 * LT - 123477085 / RT)
      K(ROX21F) = 3.7000000000D-07 &
	   * exp(4.7 * LT - 20200000 / RT)
      K(ROX21B) = 2.3145507889D-06 &
	   * exp(4.72866 * LT - 94324468.15 / RT)
      K(ROX22F) = 7.3200000000D+10 * exp(-173220000 / RT)
      K(ROX22B) = 1.5180780383D+09 &
	   * exp(0.186797 * LT - 10008887.75 / RT)
      K(ROX23F) = 2.9000000000D+10 * exp(-152400000 / RT)
      K(ROX23B) = 1.8351605777D-04 &
	   * exp(2.57193 * LT - 46692372.83 / RT)
      K(ROX24F) = 6.3400000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(ROX24B) = 5.2650873562D+13 &
	   * exp(-0.571534 * LT - 171611520.7 / RT)
      K(ROX25F) = 6.5100000000D+04 &
	   * exp(1.3 * LT - 73920000 / RT)
      K(ROX25B) = 1.0805712956D+05 &
	   * exp(1.13417 * LT - 176361125 / RT)
      K(ROX26F) = 7.4000000000D+11 * exp(-246860000 / RT)
      K(ROX26B) = 9.7808623359D-02 &
	   * exp(2.13717 * LT - 208812409.9 / RT)
      K(ROX27F) = 4.3000000000D+06 &
	   * exp(1.45 * LT - 16180000 / RT)
      K(ROX27B) = 2.1563743598D-01 &
	   * exp(2.92344 * LT - 193395741.4 / RT)
      K(ROX28) = 3.0000000000D+10 * exp(-20920000 / RT)
      K(ROX30F) = 1.8300000000D+10 * exp(-18950000 / RT)
      K(ROX30B) = 3.1309737000D+12 &
	   * exp(-0.681196 * LT - 84925002.67 / RT)
      K(ROX31F) = 1.8300000000D+10 * exp(-18950000 / RT)
      K(ROX31B) = 3.1543081311D+17 &
	   * exp(-0.802389 * LT - 449214284.2 / RT)
      K(ROX32F) = 2.2000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX32B) = 2.7062147246D+04 &
	   * exp(2.13353 * LT - 26652213.94 / RT)
      K(ROX33F) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX33B) = 2.7426874706D+21 &
	   * exp(-3.06281 * LT - 96437429.5 / RT)
      K(ROX34F) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX34B) = 2.4175569754D+21 &
	   * exp(-3.04702 * LT - 95146172.79 / RT)
      K(ROX35) = 3.0000000000D+10 * exp(-37580000 / RT)
      K(ROX36) = 3.0000000000D+10 * exp(-37580000 / RT)
      K(ROX37F) = 1.0000000000D+11
      K(ROX37B) = 2.2436382722D+19 &
	   * exp(-0.91444 * LT - 564179895.5 / RT)
      K(ROX38F) = 1.0000000000D+11
      K(ROX38B) = 1.9776673111D+19 &
	   * exp(-0.898649 * LT - 562888638.7 / RT)
      K(ROX39F) = 3.0000000000D+10
      K(ROX39B) = 4.8035314725D+16 &
	   * exp(-1.22852 * LT - 137177825.2 / RT)
      K(ROX40F) = 3.0000000000D+10
      K(ROX40B) = 4.2340992703D+16 &
	   * exp(-1.21273 * LT - 135886568.5 / RT)
      K(ROX41F) = 8.6200000000D+15 &
	   * exp(-0.61 * LT - 310110000 / RT)
      K(ROX41B) = 5.1809578901D+03 &
	   * exp(1.18494 * LT - 220339310.1 / RT)
      K(ROX42F) = 1.0100000000D+71 &
	   * exp(-15.92 * LT - 522120000 / RT)
      K(ROX42B) = 1.0025283853D+66 &
	   * exp(-15.7988 * LT - 157830718.5 / RT)
      K(ROX43F) = 6.8300000000D-02 &
	   * exp(3.4 * LT - 30260000 / RT)
      K(ROX43B) = 1.5845850043D-04 &
	   * exp(3.87143 * LT - 99062206.29 / RT)
      K(ROX44F) = 1.7300000000D-02 &
	   * exp(3.4 * LT + 4780000 / RT)
      K(ROX44B) = 4.4156394303D-04 &
	   * exp(3.76768 * LT - 127066731 / RT)
      K(ROX45F) = 3.7000000000D-07 &
	   * exp(4.7 * LT - 20200000 / RT)
      K(ROX45B) = 1.5957878251D-06 &
	   * exp(4.65122 * LT - 97914114.2 / RT)
      K(ROX46F) = 2.9000000000D+10 * exp(-152400000 / RT)
      K(ROX46B) = 1.8477233288D-02 &
	   * exp(1.86592 * LT - 87534546.43 / RT)
      K(ROX47) = 1.6800000000D+11
      K(ROX48) = 6.5100000000D+04 &
	   * exp(1.3 * LT - 73920000 / RT)
      K(ROX50) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX51) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX52) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX53) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX54) = 8.5700000000D+17 &
	   * exp(-2.27 * LT - 30080000 / RT)
      K(ROX60) = 1.1000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX61) = 1.1000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX62) = 2.2000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX63) = 1.7600000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX64) = 2.2000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX65) = 2.2000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(ROX66) = 2.2000000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)

      CALL COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )

      W(R1F) = K(R1F) * C(SH) * C(SO2)
      W(R1B) = K(R1B) * C(SOH) * C(SO)
      W(R2F) = K(R2F) * C(SO) * C(SH2)
      W(R2B) = K(R2B) * C(SOH) * C(SH)
      W(R3F) = K(R3F) * C(SOH) * C(SH2)
      W(R3B) = K(R3B) * C(SH2O) * C(SH)
      W(R4F) = K(R4F) * C(SOH) * C(SOH)
      W(R4B) = K(R4B) * C(SH2O) * C(SO)
      W(R5F) = K(R5F) * C(SH) * C(SH) * M(MM1)
      W(R5B) = K(R5B) * C(SH2) * M(MM1)
      W(R6F) = K(R6F) * C(SH) * C(SH) * C(SH2)
      W(R6B) = K(R6B) * C(SH2) * C(SH2)
      W(R7F) = K(R7F) * C(SH) * C(SH) * C(SH2O)
      W(R7B) = K(R7B) * C(SH2O) * C(SH2)
      W(R8F) = K(R8F) * C(SH) * C(SH) * C(SCO2)
      W(R8B) = K(R8B) * C(SCO2) * C(SH2)
      W(R9F) = K(R9F) * C(SH) * C(SOH) * M(MM2)
      W(R9B) = K(R9B) * C(SH2O) * M(MM2)
      W(R10F) = K(R10F) * C(SH) * C(SO) * M(MM3)
      W(R10B) = K(R10B) * C(SOH) * M(MM3)
      W(R11F) = K(R11F) * C(SO) * C(SO) * M(MM4)
      W(R11B) = K(R11B) * C(SO2) * M(MM4)
      W(R12F) = K(R12F) * C(SH) * C(SO2)
      W(R12B) = K(R12B) * C(SHO2)
      W(R13F) = K(R13F) * C(SH2) * C(SO2)
      W(R13B) = K(R13B) * C(SH) * C(SHO2)
      W(R14F) = K(R14F) * C(SOH) * C(SOH)
      W(R14B) = K(R14B) * C(SH2O2)
      W(R15F) = K(R15F) * C(SHO2) * C(SH)
      W(R15B) = K(R15B) * C(SH2O) * C(SO)
      W(R16F) = K(R16F) * C(SHO2) * C(SH)
      W(R16B) = K(R16B) * C(SOH) * C(SOH)
      W(R17F) = K(R17F) * C(SHO2) * C(SO)
      W(R17B) = K(R17B) * C(SO2) * C(SOH)
      W(R18F) = K(R18F) * C(SHO2) * C(SOH)
      W(R18B) = K(R18B) * C(SO2) * C(SH2O)
      W(R19F) = K(R19F) * C(SHO2) * C(SOH)
      W(R19B) = K(R19B) * C(SO2) * C(SH2O)
      W(R20F) = K(R20F) * C(SHO2) * C(SHO2)
      W(R20B) = K(R20B) * C(SH2O2) * C(SO2)
      W(R21F) = K(R21F) * C(SHO2) * C(SHO2)
      W(R21B) = K(R21B) * C(SH2O2) * C(SO2)
      W(R22F) = K(R22F) * C(SH2O2) * C(SH)
      W(R22B) = K(R22B) * C(SH2) * C(SHO2)
      W(R23F) = K(R23F) * C(SH2O2) * C(SH)
      W(R23B) = K(R23B) * C(SOH) * C(SH2O)
      W(R24F) = K(R24F) * C(SH2O2) * C(SO)
      W(R24B) = K(R24B) * C(SOH) * C(SHO2)
      W(R25F) = K(R25F) * C(SH2O2) * C(SOH)
      W(R25B) = K(R25B) * C(SH2O) * C(SHO2)
      W(R26F) = K(R26F) * C(SH2O2) * C(SOH)
      W(R26B) = K(R26B) * C(SH2O) * C(SHO2)
      W(R27F) = K(R27F) * C(SCO) * C(SO)
      W(R27B) = K(R27B) * C(SCO2)
      W(R28F) = K(R28F) * C(SCO) * C(SOH)
      W(R28B) = K(R28B) * C(SH) * C(SCO2)
      W(R29F) = K(R29F) * C(SCO) * C(SOH)
      W(R29B) = K(R29B) * C(SH) * C(SCO2)
      W(R30F) = K(R30F) * C(SCO) * C(SO2)
      W(R30B) = K(R30B) * C(SO) * C(SCO2)
      W(R31F) = K(R31F) * C(SCO) * C(SHO2)
      W(R31B) = K(R31B) * C(SOH) * C(SCO2)
      W(R32F) = K(R32F) * C(SHCO) * C(SH)
      W(R32B) = K(R32B) * C(SH2) * C(SCO)
      W(R33F) = K(R33F) * C(SHCO) * C(SO)
      W(R33B) = K(R33B) * C(SOH) * C(SCO)
      W(R34F) = K(R34F) * C(SHCO) * C(SO)
      W(R34B) = K(R34B) * C(SH) * C(SCO2)
      W(R35F) = K(R35F) * C(SHCO) * C(SOH)
      W(R35B) = K(R35B) * C(SH2O) * C(SCO)
      W(R36F) = K(R36F) * C(SHCO) * M(MM8)
      W(R36B) = K(R36B) * C(SH) * C(SCO) * M(MM8)
      W(R37F) = K(R37F) * C(SHCO) * C(SH2O)
      W(R37B) = K(R37B) * C(SH2O) * C(SH) * C(SCO)
      W(R38F) = K(R38F) * C(SHCO) * C(SO2)
      W(R38B) = K(R38B) * C(SHO2) * C(SCO)
      W(RG01F) = K(RG01F) * C(SC) * C(SOH)
      W(RG01B) = K(RG01B) * C(SH) * C(SCO)
      W(RG02F) = K(RG02F) * C(SC) * C(SO2)
      W(RG02B) = K(RG02B) * C(SO) * C(SCO)
      W(RG03F) = K(RG03F) * C(SCH) * C(SH)
      W(RG03B) = K(RG03B) * C(SH2) * C(SC)
      W(RG04F) = K(RG04F) * C(SCH) * C(SO)
      W(RG04B) = K(RG04B) * C(SH) * C(SCO)
      W(RG05F) = K(RG05F) * C(SCH) * C(SOH)
      W(RG05B) = K(RG05B) * C(SH) * C(SHCO)
      W(RG06F) = K(RG06F) * C(SCH) * C(SH2)
      W(RG06B) = K(RG06B) * C(SH) * C(STXCH2)
      W(RG07F) = K(RG07F) * C(SCH) * C(SH2)
      W(RG07B) = K(RG07B) * C(SCH3)
      W(RG08F) = K(RG08F) * C(SCH) * C(SH2O)
      W(RG08B) = K(RG08B) * C(SH) * C(SCH2O)
      W(RG09F) = K(RG09F) * C(SCH) * C(SO2)
      W(RG09B) = K(RG09B) * C(SO) * C(SHCO)
      W(RG10F) = K(RG10F) * C(SCH) * C(SCO)
      W(RG10B) = K(RG10B) * C(SHCCO)
      W(RG11F) = K(RG11F) * C(SCH) * C(SCO2)
      W(RG11B) = K(RG11B) * C(SCO) * C(SHCO)
      W(RG12F) = K(RG12F) * C(SCO) * C(SH2)
      W(RG12B) = K(RG12B) * C(SCH2O)
      W(RG13F) = K(RG13F) * C(SHCO) * C(SH)
      W(RG13B) = K(RG13B) * C(SCH2O)
      W(RG14F) = K(RG14F) * C(STXCH2) * C(SH)
      W(RG14B) = K(RG14B) * C(SCH3)
      W(RG15F) = K(RG15F) * C(STXCH2) * C(SO)
      W(RG15B) = K(RG15B) * C(SH) * C(SHCO)
      W(RG16F) = K(RG16F) * C(STXCH2) * C(SOH)
      W(RG16B) = K(RG16B) * C(SH) * C(SCH2O)
      W(RG17F) = K(RG17F) * C(STXCH2) * C(SOH)
      W(RG17B) = K(RG17B) * C(SH2O) * C(SCH)
      W(RG18F) = K(RG18F) * C(STXCH2) * C(SH2)
      W(RG18B) = K(RG18B) * C(SCH3) * C(SH)
      W(RG19) = K(RG19) * C(STXCH2) * C(SO2)
      W(RG20F) = K(RG20F) * C(STXCH2) * C(SO2)
      W(RG20B) = K(RG20B) * C(SO) * C(SCH2O)
      W(RG21) = K(RG21) * C(STXCH2) * C(SO2)
      W(RG22F) = K(RG22F) * C(STXCH2) * C(SHO2)
      W(RG22B) = K(RG22B) * C(SOH) * C(SCH2O)
      W(RG23F) = K(RG23F) * C(STXCH2) * C(SC)
      W(RG23B) = K(RG23B) * C(SH) * C(SC2H)
      W(RG24F) = K(RG24F) * C(STXCH2) * C(SCO)
      W(RG24B) = K(RG24B) * C(SCH2CO)
      W(RG25F) = K(RG25F) * C(STXCH2) * C(SCH)
      W(RG25B) = K(RG25B) * C(SH) * C(SC2H2)
      W(RG26F) = K(RG26F) * C(STXCH2) * C(STXCH2)
      W(RG26B) = K(RG26B) * C(SH2) * C(SC2H2)
      W(RG27) = K(RG27) * C(STXCH2) * C(STXCH2)
      W(RG28F) = K(RG28F) * C(SSXCH2) * C(SN2)
      W(RG28B) = K(RG28B) * C(SN2) * C(STXCH2)
      W(RG29F) = K(RG29F) * C(SSXCH2) * C(SAR)
      W(RG29B) = K(RG29B) * C(SAR) * C(STXCH2)
      W(RG30F) = K(RG30F) * C(SSXCH2) * C(SH)
      W(RG30B) = K(RG30B) * C(SH2) * C(SCH)
      W(RG31F) = K(RG31F) * C(SSXCH2) * C(SO)
      W(RG31B) = K(RG31B) * C(SH2) * C(SCO)
      W(RG32F) = K(RG32F) * C(SSXCH2) * C(SO)
      W(RG32B) = K(RG32B) * C(SH) * C(SHCO)
      W(RG33F) = K(RG33F) * C(SSXCH2) * C(SOH)
      W(RG33B) = K(RG33B) * C(SH) * C(SCH2O)
      W(RG34F) = K(RG34F) * C(SSXCH2) * C(SH2)
      W(RG34B) = K(RG34B) * C(SH) * C(SCH3)
      W(RG35F) = K(RG35F) * C(SSXCH2) * C(SO2)
      W(RG35B) = K(RG35B) * C(SCO) * C(SOH) * C(SH)
      W(RG36F) = K(RG36F) * C(SSXCH2) * C(SO2)
      W(RG36B) = K(RG36B) * C(SH2O) * C(SCO)
      W(RG37F) = K(RG37F) * C(SSXCH2) * C(SH2O)
      W(RG37B) = K(RG37B) * C(SCH3OH)
      W(RG38F) = K(RG38F) * C(SSXCH2) * C(SH2O)
      W(RG38B) = K(RG38B) * C(SH2O) * C(STXCH2)
      W(RG39) = K(RG39) * C(SSXCH2) * C(SH2O)
      W(RG40F) = K(RG40F) * C(SSXCH2) * C(SCO)
      W(RG40B) = K(RG40B) * C(SCO) * C(STXCH2)
      W(RG41F) = K(RG41F) * C(SSXCH2) * C(SCO2)
      W(RG41B) = K(RG41B) * C(SCO2) * C(STXCH2)
      W(RG42F) = K(RG42F) * C(SSXCH2) * C(SCO2)
      W(RG42B) = K(RG42B) * C(SCO) * C(SCH2O)
      W(RG43F) = K(RG43F) * C(SCH2O) * C(SH)
      W(RG43B) = K(RG43B) * C(SCH2OH)
      W(RG44F) = K(RG44F) * C(SCH2O) * C(SH)
      W(RG44B) = K(RG44B) * C(SCH3O)
      W(RG45F) = K(RG45F) * C(SCH2O) * C(SH)
      W(RG45B) = K(RG45B) * C(SH2) * C(SHCO)
      W(RG46F) = K(RG46F) * C(SCH2O) * C(SO)
      W(RG46B) = K(RG46B) * C(SOH) * C(SHCO)
      W(RG47F) = K(RG47F) * C(SCH2O) * C(SOH)
      W(RG47B) = K(RG47B) * C(SH2O) * C(SHCO)
      W(RG48F) = K(RG48F) * C(SCH2O) * C(SO2)
      W(RG48B) = K(RG48B) * C(SHO2) * C(SHCO)
      W(RG49F) = K(RG49F) * C(SCH2O) * C(SHO2)
      W(RG49B) = K(RG49B) * C(SH2O2) * C(SHCO)
      W(RG50F) = K(RG50F) * C(SCH2O) * C(SCH)
      W(RG50B) = K(RG50B) * C(SH) * C(SCH2CO)
      W(RG51F) = K(RG51F) * C(SCH3) * C(SH)
      W(RG51B) = K(RG51B) * C(SCH4)
      W(RG52F) = K(RG52F) * C(SCH3) * C(SO)
      W(RG52B) = K(RG52B) * C(SH) * C(SCH2O)
      W(RG53) = K(RG53) * C(SCH3) * C(SO)
      W(RG54F) = K(RG54F) * C(SCH3) * C(SOH)
      W(RG54B) = K(RG54B) * C(SCH3OH)
      W(RG55F) = K(RG55F) * C(SCH3) * C(SOH)
      W(RG55B) = K(RG55B) * C(SH2O) * C(STXCH2)
      W(RG56) = K(RG56) * C(SCH3) * C(SOH)
      W(RG57F) = K(RG57F) * C(SCH3) * C(SOH)
      W(RG57B) = K(RG57B) * C(SH2O) * C(SSXCH2)
      W(RG58F) = K(RG58F) * C(SCH3) * C(SO2)
      W(RG58B) = K(RG58B) * C(SO) * C(SCH3O)
      W(RG59F) = K(RG59F) * C(SCH3) * C(SO2)
      W(RG59B) = K(RG59B) * C(SOH) * C(SCH2O)
      W(RG60F) = K(RG60F) * C(SCH3) * C(SO2)
      W(RG60B) = K(RG60B) * C(SCH3O2)
      W(RG61F) = K(RG61F) * C(SCH3O2) * C(SCH3)
      W(RG61B) = K(RG61B) * C(SCH3O) * C(SCH3O)
      W(RG62) = K(RG62) * C(SCH3O2) * C(SCH3O2)
      W(RG63) = K(RG63) * C(SCH3O2) * C(SHO2)
      W(RG64) = K(RG64) * C(SCH3O2) * C(SCH2O)
      W(RG65F) = K(RG65F) * C(SCH3) * C(SHO2)
      W(RG65B) = K(RG65B) * C(SOH) * C(SCH3O)
      W(RG66F) = K(RG66F) * C(SCH3) * C(SHO2)
      W(RG66B) = K(RG66B) * C(SO2) * C(SCH4)
      W(RG67F) = K(RG67F) * C(SCH3) * C(SH2O2)
      W(RG67B) = K(RG67B) * C(SHO2) * C(SCH4)
      W(RG68F) = K(RG68F) * C(SCH3) * C(SC)
      W(RG68B) = K(RG68B) * C(SH) * C(SC2H2)
      W(RG69F) = K(RG69F) * C(SCH3) * C(SCH)
      W(RG69B) = K(RG69B) * C(SH) * C(SC2H3)
      W(RG70F) = K(RG70F) * C(SCH3) * C(SHCO)
      W(RG70B) = K(RG70B) * C(SCO) * C(SCH4)
      W(RG71F) = K(RG71F) * C(SCH3) * C(SCH2O)
      W(RG71B) = K(RG71B) * C(SHCO) * C(SCH4)
      W(RG72F) = K(RG72F) * C(SCH3) * C(STXCH2)
      W(RG72B) = K(RG72B) * C(SH) * C(SC2H4)
      W(RG73F) = K(RG73F) * C(SCH3) * C(SSXCH2)
      W(RG73B) = K(RG73B) * C(SH) * C(SC2H4)
      W(RG74F) = K(RG74F) * C(SCH3) * C(SCH3)
      W(RG74B) = K(RG74B) * C(SH) * C(SC2H5)
      W(RG75F) = K(RG75F) * C(SCH3O) * C(SH)
      W(RG75B) = K(RG75B) * C(SCH3OH)
      W(RG76F) = K(RG76F) * C(SCH3O) * C(SH)
      W(RG76B) = K(RG76B) * C(SH) * C(SCH2OH)
      W(RG77F) = K(RG77F) * C(SCH3O) * C(SH)
      W(RG77B) = K(RG77B) * C(SH2) * C(SCH2O)
      W(RG78F) = K(RG78F) * C(SCH3O) * C(SH)
      W(RG78B) = K(RG78B) * C(SOH) * C(SCH3)
      W(RG79F) = K(RG79F) * C(SCH3O) * C(SH)
      W(RG79B) = K(RG79B) * C(SH2O) * C(SSXCH2)
      W(RG80F) = K(RG80F) * C(SCH3O) * C(SO)
      W(RG80B) = K(RG80B) * C(SOH) * C(SCH2O)
      W(RG81F) = K(RG81F) * C(SCH3O) * C(SOH)
      W(RG81B) = K(RG81B) * C(SH2O) * C(SCH2O)
      W(RG82F) = K(RG82F) * C(SCH3O) * C(SO2)
      W(RG82B) = K(RG82B) * C(SHO2) * C(SCH2O)
      W(RG83F) = K(RG83F) * C(SCH2OH) * C(SH)
      W(RG83B) = K(RG83B) * C(SCH3OH)
      W(RG84F) = K(RG84F) * C(SCH2OH) * C(SH)
      W(RG84B) = K(RG84B) * C(SH2) * C(SCH2O)
      W(RG85F) = K(RG85F) * C(SCH2OH) * C(SH)
      W(RG85B) = K(RG85B) * C(SOH) * C(SCH3)
      W(RG86F) = K(RG86F) * C(SCH2OH) * C(SH)
      W(RG86B) = K(RG86B) * C(SH2O) * C(SSXCH2)
      W(RG87F) = K(RG87F) * C(SCH2OH) * C(SO)
      W(RG87B) = K(RG87B) * C(SOH) * C(SCH2O)
      W(RG88F) = K(RG88F) * C(SCH2OH) * C(SOH)
      W(RG88B) = K(RG88B) * C(SH2O) * C(SCH2O)
      W(RG89F) = K(RG89F) * C(SCH2OH) * C(SO2)
      W(RG89B) = K(RG89B) * C(SHO2) * C(SCH2O)
      W(RG90F) = K(RG90F) * C(SCH4) * C(SH)
      W(RG90B) = K(RG90B) * C(SH2) * C(SCH3)
      W(RG91F) = K(RG91F) * C(SCH4) * C(SO)
      W(RG91B) = K(RG91B) * C(SOH) * C(SCH3)
      W(RG92F) = K(RG92F) * C(SCH4) * C(SOH)
      W(RG92B) = K(RG92B) * C(SH2O) * C(SCH3)
      W(RG93F) = K(RG93F) * C(SCH4) * C(SCH)
      W(RG93B) = K(RG93B) * C(SH) * C(SC2H4)
      W(RG94F) = K(RG94F) * C(SCH4) * C(STXCH2)
      W(RG94B) = K(RG94B) * C(SCH3) * C(SCH3)
      W(RG95F) = K(RG95F) * C(SCH4) * C(SSXCH2)
      W(RG95B) = K(RG95B) * C(SCH3) * C(SCH3)
      W(RG96F) = K(RG96F) * C(SCH3OH) * C(SH)
      W(RG96B) = K(RG96B) * C(SH2) * C(SCH2OH)
      W(RG97F) = K(RG97F) * C(SCH3OH) * C(SH)
      W(RG97B) = K(RG97B) * C(SH2) * C(SCH3O)
      W(RG98F) = K(RG98F) * C(SCH3OH) * C(SO)
      W(RG98B) = K(RG98B) * C(SOH) * C(SCH2OH)
      W(RG99F) = K(RG99F) * C(SCH3OH) * C(SO)
      W(RG99B) = K(RG99B) * C(SOH) * C(SCH3O)
      W(RG100F) = K(RG100F) * C(SCH3OH) * C(SOH)
      W(RG100B) = K(RG100B) * C(SH2O) * C(SCH2OH)
      W(RG101F) = K(RG101F) * C(SCH3OH) * C(SOH)
      W(RG101B) = K(RG101B) * C(SH2O) * C(SCH3O)
      W(RG102F) = K(RG102F) * C(SCH3OH) * C(SCH3)
      W(RG102B) = K(RG102B) * C(SCH4) * C(SCH2OH)
      W(RG103F) = K(RG103F) * C(SCH3OH) * C(SCH3)
      W(RG103B) = K(RG103B) * C(SCH4) * C(SCH3O)
      W(RG104F) = K(RG104F) * C(SC2H) * C(SH)
      W(RG104B) = K(RG104B) * C(SC2H2)
      W(RG105F) = K(RG105F) * C(SC2H) * C(SO)
      W(RG105B) = K(RG105B) * C(SCO) * C(SCH)
      W(RG106F) = K(RG106F) * C(SC2H) * C(SOH)
      W(RG106B) = K(RG106B) * C(SHCCO) * C(SH)
      W(RG107F) = K(RG107F) * C(SC2H) * C(SO2)
      W(RG107B) = K(RG107B) * C(SCO) * C(SHCO)
      W(RG108F) = K(RG108F) * C(SC2H) * C(SH2)
      W(RG108B) = K(RG108B) * C(SH) * C(SC2H2)
      W(RG109F) = K(RG109F) * C(SHCCO) * C(SH)
      W(RG109B) = K(RG109B) * C(SCO) * C(SSXCH2)
      W(RG110F) = K(RG110F) * C(SHCCO) * C(SO)
      W(RG110B) = K(RG110B) * C(SCO) * C(SCO) * C(SH)
      W(RG111F) = K(RG111F) * C(SHCCO) * C(SO2)
      W(RG111B) = K(RG111B) * C(SCO) * C(SCO) * C(SOH)
      W(RG112F) = K(RG112F) * C(SHCCO) * C(SCH)
      W(RG112B) = K(RG112B) * C(SCO) * C(SC2H2)
      W(RG113F) = K(RG113F) * C(SHCCO) * C(STXCH2)
      W(RG113B) = K(RG113B) * C(SCO) * C(SC2H3)
      W(RG114F) = K(RG114F) * C(SHCCO) * C(SHCCO)
      W(RG114B) = K(RG114B) * C(SCO) * C(SCO) * C(SC2H2)
      W(RG115F) = K(RG115F) * C(SC2H2) * C(SH)
      W(RG115B) = K(RG115B) * C(SC2H3)
      W(RG116F) = K(RG116F) * C(SC2H2) * C(SO)
      W(RG116B) = K(RG116B) * C(SH) * C(SHCCO)
      W(RG117F) = K(RG117F) * C(SC2H2) * C(SO)
      W(RG117B) = K(RG117B) * C(SCO) * C(STXCH2)
      W(RG118F) = K(RG118F) * C(SC2H) * C(SOH)
      W(RG118B) = K(RG118B) * C(SO) * C(SC2H2)
      W(RG119F) = K(RG119F) * C(SC2H2) * C(SOH)
      W(RG119B) = K(RG119B) * C(SH2O) * C(SC2H)
      W(RG120F) = K(RG120F) * C(SC2H2) * C(SOH)
      W(RG120B) = K(RG120B) * C(SH) * C(SHCCOH)
      W(RG121F) = K(RG121F) * C(SC2H2) * C(SOH)
      W(RG121B) = K(RG121B) * C(SH) * C(SCH2CO)
      W(RG122F) = K(RG122F) * C(SC2H2) * C(SOH)
      W(RG122B) = K(RG122B) * C(SCO) * C(SCH3)
      W(RG123F) = K(RG123F) * C(SCH2CO) * C(SH)
      W(RG123B) = K(RG123B) * C(SH2) * C(SHCCO)
      W(RG124F) = K(RG124F) * C(SCH2CO) * C(SH)
      W(RG124B) = K(RG124B) * C(SCO) * C(SCH3)
      W(RG125F) = K(RG125F) * C(SCH2CO) * C(SO)
      W(RG125B) = K(RG125B) * C(SOH) * C(SHCCO)
      W(RG126F) = K(RG126F) * C(SCH2CO) * C(SO)
      W(RG126B) = K(RG126B) * C(SCO2) * C(STXCH2)
      W(RG127F) = K(RG127F) * C(SCH2CO) * C(SOH)
      W(RG127B) = K(RG127B) * C(SH2O) * C(SHCCO)
      W(RG128F) = K(RG128F) * C(SHCCOH) * C(SH)
      W(RG128B) = K(RG128B) * C(SH) * C(SCH2CO)
      W(RG129F) = K(RG129F) * C(SC2H3) * C(SH)
      W(RG129B) = K(RG129B) * C(SC2H4)
      W(RG130F) = K(RG130F) * C(SC2H3) * C(SH)
      W(RG130B) = K(RG130B) * C(SH2) * C(SC2H2)
      W(RG131F) = K(RG131F) * C(SC2H3) * C(SO)
      W(RG131B) = K(RG131B) * C(SCH2CHO)
      W(RG132F) = K(RG132F) * C(SC2H3) * C(SOH)
      W(RG132B) = K(RG132B) * C(SH2O) * C(SC2H2)
      W(RG133F) = K(RG133F) * C(SC2H3) * C(SO2)
      W(RG133B) = K(RG133B) * C(SHO2) * C(SC2H2)
      W(RG134F) = K(RG134F) * C(SC2H3) * C(SO2)
      W(RG134B) = K(RG134B) * C(SO) * C(SCH2CHO)
      W(RG135F) = K(RG135F) * C(SC2H3) * C(SO2)
      W(RG135B) = K(RG135B) * C(SCH2O) * C(SHCO)
      W(RG136F) = K(RG136F) * C(SCH2CHO)
      W(RG136B) = K(RG136B) * C(SH) * C(SCH2CO)
      W(RG137F) = K(RG137F) * C(SCH2CHO)
      W(RG137B) = K(RG137B) * C(SCO) * C(SCH3)
      W(RG138F) = K(RG138F) * C(SCH2CHO) * C(SO)
      W(RG138B) = K(RG138B) * C(SHCO) * C(SCH2O)
      W(RG139) = K(RG139) * C(SCH2CHO) * C(SO2)
      W(RG140) = K(RG140) * C(SCH2CHO) * C(SO2)
      W(RG141F) = K(RG141F) * C(SCH2CHO) * C(SH)
      W(RG141B) = K(RG141B) * C(SHCO) * C(SCH3)
      W(RG142F) = K(RG142F) * C(SCH2CHO) * C(SH)
      W(RG142B) = K(RG142B) * C(SH2) * C(SCH2CO)
      W(RG143F) = K(RG143F) * C(SCH2CHO) * C(SOH)
      W(RG143B) = K(RG143B) * C(SCH2CO) * C(SH2O)
      W(RG144F) = K(RG144F) * C(SCH2CHO) * C(SOH)
      W(RG144B) = K(RG144B) * C(SCH2OH) * C(SHCO)
      W(RG145F) = K(RG145F) * C(SCH3) * C(SHCO)
      W(RG145B) = K(RG145B) * C(SCH3CHO)
      W(RG146F) = K(RG146F) * C(SCH3CHO) * C(SO)
      W(RG146B) = K(RG146B) * C(SOH) * C(SCH2CHO)
      W(RG147F) = K(RG147F) * C(SCH3CHO) * C(SH)
      W(RG147B) = K(RG147B) * C(SH2) * C(SCH2CHO)
      W(RG148) = K(RG148) * C(SCH3CHO) * C(SH)
      W(RG149) = K(RG149) * C(SCH3CHO) * C(SO)
      W(RG150) = K(RG150) * C(SCH3CHO) * C(SO2)
      W(RG151) = K(RG151) * C(SCH3CHO) * C(SOH)
      W(RG152) = K(RG152) * C(SCH3CHO) * C(SHO2)
      W(RG153) = K(RG153) * C(SCH3CHO) * C(SCH3)
      W(RG154F) = K(RG154F) * C(SC2H4)
      W(RG154B) = K(RG154B) * C(SH2) * C(SH2C2)
      W(RG155F) = K(RG155F) * C(SC2H4) * C(SH)
      W(RG155B) = K(RG155B) * C(SC2H5)
      W(RG156F) = K(RG156F) * C(SC2H4) * C(SH)
      W(RG156B) = K(RG156B) * C(SH2) * C(SC2H3)
      W(RG157F) = K(RG157F) * C(SC2H4) * C(SO)
      W(RG157B) = K(RG157B) * C(SH) * C(SCH2CHO)
      W(RG158F) = K(RG158F) * C(SC2H4) * C(SO)
      W(RG158B) = K(RG158B) * C(SCH2O) * C(STXCH2)
      W(RG159F) = K(RG159F) * C(SC2H4) * C(SO)
      W(RG159B) = K(RG159B) * C(SHCO) * C(SCH3)
      W(RG160F) = K(RG160F) * C(SC2H4) * C(SOH)
      W(RG160B) = K(RG160B) * C(SH2O) * C(SC2H3)
      W(RG161F) = K(RG161F) * C(SC2H4) * C(SOH)
      W(RG161B) = K(RG161B) * C(SC2H5O)
      W(RG162F) = K(RG162F) * C(SC2H4) * C(SCH3)
      W(RG162B) = K(RG162B) * C(SCH4) * C(SC2H3)
      W(RG163F) = K(RG163F) * C(SC2H4) * C(SCH3)
      W(RG163B) = K(RG163B) * C(SNXC3H7)
      W(RG164F) = K(RG164F) * C(SC2H5) * C(SH)
      W(RG164B) = K(RG164B) * C(SC2H6)
      W(RG165F) = K(RG165F) * C(SC2H5) * C(SH)
      W(RG165B) = K(RG165B) * C(SH2) * C(SC2H4)
      W(RG166F) = K(RG166F) * C(SC2H5) * C(SCH3)
      W(RG166B) = K(RG166B) * C(SCH4) * C(SC2H4)
      W(RG167F) = K(RG167F) * C(SC2H5) * C(SO)
      W(RG167B) = K(RG167B) * C(SC2H5O)
      W(RG168F) = K(RG168F) * C(SC2H5O)
      W(RG168B) = K(RG168B) * C(SCH2O) * C(SCH3)
      W(RG169F) = K(RG169F) * C(SC2H5O)
      W(RG169B) = K(RG169B) * C(SH) * C(SCH3CHO)
      W(RG170F) = K(RG170F) * C(SC2H5O) * C(SO2)
      W(RG170B) = K(RG170B) * C(SHO2) * C(SCH3CHO)
      W(RG171F) = K(RG171F) * C(SC2H5) * C(SO2)
      W(RG171B) = K(RG171B) * C(SHO2) * C(SC2H4)
      W(RG172F) = K(RG172F) * C(SC3H8)
      W(RG172B) = K(RG172B) * C(SCH3) * C(SC2H5)
      W(RG173F) = K(RG173F) * C(SC2H6)
      W(RG173B) = K(RG173B) * C(SCH3) * C(SCH3)
      W(RG174F) = K(RG174F) * C(SC2H6) * C(SH)
      W(RG174B) = K(RG174B) * C(SH2) * C(SC2H5)
      W(RG175F) = K(RG175F) * C(SC2H6) * C(SO)
      W(RG175B) = K(RG175B) * C(SOH) * C(SC2H5)
      W(RG176F) = K(RG176F) * C(SC2H6) * C(SOH)
      W(RG176B) = K(RG176B) * C(SH2O) * C(SC2H5)
      W(RG177F) = K(RG177F) * C(SC2H6) * C(SSXCH2)
      W(RG177B) = K(RG177B) * C(SCH3) * C(SC2H5)
      W(RG178F) = K(RG178F) * C(SC2H6) * C(SCH3)
      W(RG178B) = K(RG178B) * C(SCH4) * C(SC2H5)
      W(RG179F) = K(RG179F) * C(SNXC3H7) * C(SO)
      W(RG179B) = K(RG179B) * C(SCH2O) * C(SC2H5)
      W(RG180F) = K(RG180F) * C(SNXC3H7) * C(SH)
      W(RG180B) = K(RG180B) * C(SC3H8)
      W(RG181F) = K(RG181F) * C(SNXC3H7) * C(SOH)
      W(RG181B) = K(RG181B) * C(SH2O) * C(SC3H6)
      W(RG182F) = K(RG182F) * C(SNXC3H7) * C(SCH3)
      W(RG182B) = K(RG182B) * C(SCH4) * C(SC3H6)
      W(RG183F) = K(RG183F) * C(SC3H6) * C(SH)
      W(RG183B) = K(RG183B) * C(SNXC3H7)
      W(RG184F) = K(RG184F) * C(SNXC3H7) * C(SO2)
      W(RG184B) = K(RG184B) * C(SHO2) * C(SC3H6)
      W(RG185F) = K(RG185F) * C(SC3H8) * C(SH)
      W(RG185B) = K(RG185B) * C(SH2) * C(SNXC3H7)
      W(RG186F) = K(RG186F) * C(SC3H8) * C(SO)
      W(RG186B) = K(RG186B) * C(SOH) * C(SNXC3H7)
      W(RG187F) = K(RG187F) * C(SC3H8) * C(SOH)
      W(RG187B) = K(RG187B) * C(SH2O) * C(SNXC3H7)
      W(RG188F) = K(RG188F) * C(SC3H8) * C(SCH3)
      W(RG188B) = K(RG188B) * C(SCH4) * C(SNXC3H7)
      W(RG189F) = K(RG189F) * C(SC3H8) * C(SHO2)
      W(RG189B) = K(RG189B) * C(SH2O2) * C(SNXC3H7)
      W(RR001F) = K(RR001F) * C(SC2H2) * M(MM6)
      W(RR001B) = K(RR001B) * C(SH2C2) * M(MM6)
      W(RR002F) = K(RR002F) * C(SH2C2) * C(SO2)
      W(RR002B) = K(RR002B) * C(SCO2) * C(STXCH2)
      W(RR003F) = K(RR003F) * C(SH2C2) * C(SO2)
      W(RR003B) = K(RR003B) * C(SHCO) * C(SHCO)
      W(RR004F) = K(RR004F) * C(SC2H2) * C(SSXCH2)
      W(RR004B) = K(RR004B) * C(SH) * C(SC3H3)
      W(RR005F) = K(RR005F) * C(SPXC3H4) * C(SH)
      W(RR005B) = K(RR005B) * C(SCH3) * C(SC2H2)
      W(RR006F) = K(RR006F) * C(SAXC3H4) * C(SH)
      W(RR006B) = K(RR006B) * C(SCH3) * C(SC2H2)
      W(RR007F) = K(RR007F) * C(SC2H2) * C(SCH3)
      W(RR007B) = K(RR007B) * C(SSXC3H5)
      W(RR008F) = K(RR008F) * C(SC2H2) * C(SC2H)
      W(RR008B) = K(RR008B) * C(SNXC4H3)
      W(RR009F) = K(RR009F) * C(SC2H2) * C(SHCCO)
      W(RR009B) = K(RR009B) * C(SCO) * C(SC3H3)
      W(RR010F) = K(RR010F) * C(SC2H3) * C(SH2O2)
      W(RR010B) = K(RR010B) * C(SHO2) * C(SC2H4)
      W(RR011F) = K(RR011F) * C(SC2H3) * C(SHCO)
      W(RR011B) = K(RR011B) * C(SCO) * C(SC2H4)
      W(RR012F) = K(RR012F) * C(SC2H3) * C(SHCO)
      W(RR012B) = K(RR012B) * C(SC2H3CHO)
      W(RR013F) = K(RR013F) * C(SC2H3) * C(SCH3)
      W(RR013B) = K(RR013B) * C(SCH4) * C(SC2H2)
      W(RR014F) = K(RR014F) * C(SC3H6)
      W(RR014B) = K(RR014B) * C(SCH3) * C(SC2H3)
      W(RR015F) = K(RR015F) * C(SC2H3) * C(SCH3)
      W(RR015B) = K(RR015B) * C(SH) * C(SAXC3H5)
      W(RR016F) = K(RR016F) * C(SAXC3H5) * C(SH)
      W(RR016B) = K(RR016B) * C(SC3H6)
      W(RR017F) = K(RR017F) * C(SC2H) * C(SCH3)
      W(RR017B) = K(RR017B) * C(SH) * C(SC3H3)
      W(RR018F) = K(RR018F) * C(SC2O) * C(SH)
      W(RR018B) = K(RR018B) * C(SCO) * C(SCH)
      W(RR019F) = K(RR019F) * C(SC2O) * C(SO)
      W(RR019B) = K(RR019B) * C(SCO) * C(SCO)
      W(RR020F) = K(RR020F) * C(SC2O) * C(SOH)
      W(RR020B) = K(RR020B) * C(SCO) * C(SCO) * C(SH)
      W(RR021F) = K(RR021F) * C(SC2O) * C(SO2)
      W(RR021B) = K(RR021B) * C(SCO) * C(SCO) * C(SO)
      W(RR022F) = K(RR022F) * C(SHCCO) * C(SCH3)
      W(RR022B) = K(RR022B) * C(SCO) * C(SC2H4)
      W(RR023F) = K(RR023F) * C(SHCCO) * C(SOH)
      W(RR023B) = K(RR023B) * C(SH2O) * C(SC2O)
      W(RR024F) = K(RR024F) * C(SHCCO) * C(SOH)
      W(RR024B) = K(RR024B) * C(SHCO) * C(SHCO)
      W(RR025F) = K(RR025F) * C(SCH2CO) * C(SOH)
      W(RR025B) = K(RR025B) * C(SCO) * C(SCH2OH)
      W(RR026F) = K(RR026F) * C(SCH2CO) * C(STXCH2)
      W(RR026B) = K(RR026B) * C(SCO) * C(SC2H4)
      W(RR027F) = K(RR027F) * C(SCH2CO) * C(STXCH2)
      W(RR027B) = K(RR027B) * C(SCH3) * C(SHCCO)
      W(RR028F) = K(RR028F) * C(SCH2CO) * C(SCH3)
      W(RR028B) = K(RR028B) * C(SCO) * C(SC2H5)
      W(RR029F) = K(RR029F) * C(SCH2CO) * C(SCH3)
      W(RR029B) = K(RR029B) * C(SCH4) * C(SHCCO)
      W(RR030F) = K(RR030F) * C(SCH2CHO) * C(SCH3)
      W(RR030B) = K(RR030B) * C(SHCO) * C(SC2H5)
      W(RR031F) = K(RR031F) * C(SC2H4) * C(SC2H)
      W(RR031B) = K(RR031B) * C(SH) * C(SC4H4)
      W(RR032F) = K(RR032F) * C(SC2H4) * C(SO2)
      W(RR032B) = K(RR032B) * C(SHO2) * C(SC2H3)
      W(RR033) = K(RR033) * C(SC2H4) * C(SO2)
      W(RR034F) = K(RR034F) * C(SC2H5) * C(SHCO)
      W(RR034B) = K(RR034B) * C(SCO) * C(SC2H6)
      W(RR035F) = K(RR035F) * C(SC2H5) * C(SHO2)
      W(RR035B) = K(RR035B) * C(SO2) * C(SC2H6)
      W(RR036F) = K(RR036F) * C(SC2H5) * C(SHO2)
      W(RR036B) = K(RR036B) * C(SH2O2) * C(SC2H4)
      W(RR037F) = K(RR037F) * C(SC2H5) * C(SHO2)
      W(RR037B) = K(RR037B) * C(SOH) * C(SC2H5O)
      W(RR038F) = K(RR038F) * C(SC2H6) * C(SHO2)
      W(RR038B) = K(RR038B) * C(SH2O2) * C(SC2H5)
      W(RR039F) = K(RR039F) * C(SC3H2) * C(SO)
      W(RR039B) = K(RR039B) * C(SC3H2O)
      W(RR040F) = K(RR040F) * C(SC3H2) * C(SOH)
      W(RR040B) = K(RR040B) * C(SHCO) * C(SC2H2)
      W(RR041F) = K(RR041F) * C(SC3H2) * C(SO2)
      W(RR041B) = K(RR041B) * C(SH) * C(SCO) * C(SHCCO)
      W(RR042F) = K(RR042F) * C(SC3H2) * C(SCH)
      W(RR042B) = K(RR042B) * C(SH) * C(SC4H2)
      W(RR043F) = K(RR043F) * C(SC3H2) * C(STXCH2)
      W(RR043B) = K(RR043B) * C(SH) * C(SNXC4H3)
      W(RR044F) = K(RR044F) * C(SC3H2) * C(SCH3)
      W(RR044B) = K(RR044B) * C(SH) * C(SC4H4)
      W(RR045F) = K(RR045F) * C(SC3H2) * C(SHCCO)
      W(RR045B) = K(RR045B) * C(SCO) * C(SNXC4H3)
      W(RR046F) = K(RR046F) * C(SC2H) * C(SHCO)
      W(RR046B) = K(RR046B) * C(SC3H2O)
      W(RR200F) = K(RR200F) * C(SC3H2O) * C(SH)
      W(RR200B) = K(RR200B) * C(SHCO) * C(SC2H2)
      W(RR047) = K(RR047) * C(SC3H2O) * C(SH)
      W(RR048) = K(RR048) * C(SC3H2O) * C(SO)
      W(RR049) = K(RR049) * C(SC3H2O) * C(SO2)
      W(RR050) = K(RR050) * C(SC3H2O) * C(SOH)
      W(RR051) = K(RR051) * C(SC3H2O) * C(SHO2)
      W(RR052) = K(RR052) * C(SC3H2O) * C(SCH3)
      W(RR053F) = K(RR053F) * C(SC3H2) * C(SH)
      W(RR053B) = K(RR053B) * C(SC3H3)
      W(RR054F) = K(RR054F) * C(SC3H3) * C(SH)
      W(RR054B) = K(RR054B) * C(SH2) * C(SC3H2)
      W(RR055F) = K(RR055F) * C(SC3H3) * C(SH)
      W(RR055B) = K(RR055B) * C(SPXC3H4)
      W(RR056F) = K(RR056F) * C(SC3H3) * C(SH)
      W(RR056B) = K(RR056B) * C(SAXC3H4)
      W(RR057F) = K(RR057F) * C(SC3H3) * C(SOH)
      W(RR057B) = K(RR057B) * C(SC2H3CHO)
      W(RR058F) = K(RR058F) * C(SC3H3) * C(SOH)
      W(RR058B) = K(RR058B) * C(SCO) * C(SC2H4)
      W(RR059F) = K(RR059F) * C(SC3H3) * C(SOH)
      W(RR059B) = K(RR059B) * C(SH2O) * C(SC3H2)
      W(RR060F) = K(RR060F) * C(SC3H3) * C(SOH)
      W(RR060B) = K(RR060B) * C(SC2H2) * C(SCH2O)
      W(RR061F) = K(RR061F) * C(SC3H3) * C(SO)
      W(RR061B) = K(RR061B) * C(SH) * C(SC3H2O)
      W(RR062F) = K(RR062F) * C(SC3H3) * C(SO2)
      W(RR062B) = K(RR062B) * C(SHCO) * C(SCH2CO)
      W(RR063F) = K(RR063F) * C(SC3H3) * C(SHO2)
      W(RR063B) = K(RR063B) * C(SC2H3) * C(SCO) * C(SOH)
      W(RR064F) = K(RR064F) * C(SC3H3) * C(SHO2)
      W(RR064B) = K(RR064B) * C(SO2) * C(SAXC3H4)
      W(RR065F) = K(RR065F) * C(SC3H3) * C(SHO2)
      W(RR065B) = K(RR065B) * C(SO2) * C(SPXC3H4)
      W(RR066F) = K(RR066F) * C(SPXC3H4) * C(SO2)
      W(RR066B) = K(RR066B) * C(SCO) * C(SHCO) * C(SCH3)
      W(RR067F) = K(RR067F) * C(SC3H3) * C(SHCO)
      W(RR067B) = K(RR067B) * C(SCO) * C(SAXC3H4)
      W(RR068F) = K(RR068F) * C(SC3H3) * C(SHCO)
      W(RR068B) = K(RR068B) * C(SCO) * C(SPXC3H4)
      W(RR069F) = K(RR069F) * C(SC3H3) * C(SCH)
      W(RR069B) = K(RR069B) * C(SH) * C(SIXC4H3)
      W(RR070F) = K(RR070F) * C(SC3H3) * C(STXCH2)
      W(RR070B) = K(RR070B) * C(SH) * C(SC4H4)
      W(RR071F) = K(RR071F) * C(SAXC3H4)
      W(RR071B) = K(RR071B) * C(SPXC3H4)
      W(RR073F) = K(RR073F) * C(SAXC3H4) * C(SH)
      W(RR073B) = K(RR073B) * C(SH) * C(SPXC3H4)
      W(RR074F) = K(RR074F) * C(SAXC3H4) * C(SH)
      W(RR074B) = K(RR074B) * C(SAXC3H5)
      W(RR075F) = K(RR075F) * C(SAXC3H4) * C(SH)
      W(RR075B) = K(RR075B) * C(STXC3H5)
      W(RR076F) = K(RR076F) * C(SPXC3H4) * C(SH)
      W(RR076B) = K(RR076B) * C(STXC3H5)
      W(RR077F) = K(RR077F) * C(SPXC3H4) * C(SH)
      W(RR077B) = K(RR077B) * C(SSXC3H5)
      W(RR078F) = K(RR078F) * C(SPXC3H4) * C(SH)
      W(RR078B) = K(RR078B) * C(SH2) * C(SC3H3)
      W(RR079F) = K(RR079F) * C(SPXC3H4) * C(SO)
      W(RR079B) = K(RR079B) * C(SOH) * C(SC3H3)
      W(RR080F) = K(RR080F) * C(SPXC3H4) * C(SOH)
      W(RR080B) = K(RR080B) * C(SH2O) * C(SC3H3)
      W(RR081F) = K(RR081F) * C(SPXC3H4) * C(SCH3)
      W(RR081B) = K(RR081B) * C(SCH4) * C(SC3H3)
      W(RR082F) = K(RR082F) * C(SPXC3H4) * C(SHO2)
      W(RR082B) = K(RR082B) * C(SH2O2) * C(SC3H3)
      W(RR083F) = K(RR083F) * C(SAXC3H4) * C(SH)
      W(RR083B) = K(RR083B) * C(SH2) * C(SC3H3)
      W(RR084F) = K(RR084F) * C(SAXC3H4) * C(SOH)
      W(RR084B) = K(RR084B) * C(SH2O) * C(SC3H3)
      W(RR085F) = K(RR085F) * C(SAXC3H4) * C(SCH3)
      W(RR085B) = K(RR085B) * C(SCH4) * C(SC3H3)
      W(RR086F) = K(RR086F) * C(SAXC3H4) * C(SHO2)
      W(RR086B) = K(RR086B) * C(SH2O2) * C(SC3H3)
      W(RR087F) = K(RR087F) * C(SAXC3H4) * C(SO)
      W(RR087B) = K(RR087B) * C(STXCH2) * C(SCH2CO)
      W(RR088F) = K(RR088F) * C(SPXC3H4) * C(SO)
      W(RR088B) = K(RR088B) * C(SCH3) * C(SHCCO)
      W(RR089F) = K(RR089F) * C(SPXC3H4) * C(SO)
      W(RR089B) = K(RR089B) * C(SCO) * C(SC2H4)
      W(RR090F) = K(RR090F) * C(SAXC3H4) * C(SC2H)
      W(RR090B) = K(RR090B) * C(SC3H3) * C(SC2H2)
      W(RR091F) = K(RR091F) * C(SPXC3H4) * C(SC2H)
      W(RR091B) = K(RR091B) * C(SC3H3) * C(SC2H2)
      W(RR092F) = K(RR092F) * C(SPXC3H4) * C(SOH)
      W(RR092B) = K(RR092B) * C(SCH3) * C(SHCCOH)
      W(RR093F) = K(RR093F) * C(SPXC3H4) * C(SOH)
      W(RR093B) = K(RR093B) * C(SCH3) * C(SCH2CO)
      W(RR094F) = K(RR094F) * C(SPXC3H4) * C(SOH)
      W(RR094B) = K(RR094B) * C(SCO) * C(SC2H5)
      W(RR095) = K(RR095) * C(SC2H3CHO) * C(SH)
      W(RR096) = K(RR096) * C(SC2H3CHO) * C(SO)
      W(RR097) = K(RR097) * C(SC2H3CHO) * C(SOH)
      W(RR098) = K(RR098) * C(SC2H3CHO) * C(SHO2)
      W(RR099) = K(RR099) * C(SC2H3CHO) * C(SCH3)
      W(RR100F) = K(RR100F) * C(SAXC3H5)
      W(RR100B) = K(RR100B) * C(STXC3H5)
      W(RR101F) = K(RR101F) * C(SAXC3H5)
      W(RR101B) = K(RR101B) * C(SSXC3H5)
      W(RR102F) = K(RR102F) * C(STXC3H5)
      W(RR102B) = K(RR102B) * C(SSXC3H5)
      W(RR103F) = K(RR103F) * C(SAXC3H5) * C(SH)
      W(RR103B) = K(RR103B) * C(SH2) * C(SAXC3H4)
      W(RR104F) = K(RR104F) * C(SAXC3H5) * C(SOH)
      W(RR104B) = K(RR104B) * C(SH2O) * C(SAXC3H4)
      W(RR105F) = K(RR105F) * C(SAXC3H5) * C(SCH3)
      W(RR105B) = K(RR105B) * C(SCH4) * C(SAXC3H4)
      W(RR106F) = K(RR106F) * C(SAXC3H5) * C(SC2H3)
      W(RR106B) = K(RR106B) * C(SC2H4) * C(SAXC3H4)
      W(RR107F) = K(RR107F) * C(SAXC3H5) * C(SC2H5)
      W(RR107B) = K(RR107B) * C(SC2H6) * C(SAXC3H4)
      W(RR108F) = K(RR108F) * C(SAXC3H5) * C(SAXC3H5)
      W(RR108B) = K(RR108B) * C(SC3H6) * C(SAXC3H4)
      W(RR109F) = K(RR109F) * C(SAXC3H5) * C(SO2)
      W(RR109B) = K(RR109B) * C(SHO2) * C(SAXC3H4)
      W(RR110F) = K(RR110F) * C(SAXC3H5) * C(SO2)
      W(RR110B) = K(RR110B) * C(SOH) * C(SC2H3CHO)
      W(RR111) = K(RR111) * C(SAXC3H5) * C(SO2)
      W(RR112F) = K(RR112F) * C(SAXC3H5) * C(SO2)
      W(RR112B) = K(RR112B) * C(SCH2O) * C(SCH2CHO)
      W(RR113F) = K(RR113F) * C(SAXC3H5) * C(SO)
      W(RR113B) = K(RR113B) * C(SC3H5O)
      W(RR114F) = K(RR114F) * C(SAXC3H5) * C(SOH)
      W(RR114B) = K(RR114B) * C(SH2) * C(SC2H3CHO)
      W(RR115F) = K(RR115F) * C(SAXC3H5) * C(SHCO)
      W(RR115B) = K(RR115B) * C(SCO) * C(SC3H6)
      W(RR116F) = K(RR116F) * C(SAXC3H5) * C(SHO2)
      W(RR116B) = K(RR116B) * C(SO2) * C(SC3H6)
      W(RR117F) = K(RR117F) * C(SAXC3H5) * C(SHO2)
      W(RR117B) = K(RR117B) * C(SOH) * C(SC3H5O)
      W(RR118F) = K(RR118F) * C(SAXC3H5) * C(SCH3O2)
      W(RR118B) = K(RR118B) * C(SCH3O) * C(SC3H5O)
      W(RR119F) = K(RR119F) * C(STXC3H5) * C(SH)
      W(RR119B) = K(RR119B) * C(SH2) * C(SPXC3H4)
      W(RR120F) = K(RR120F) * C(STXC3H5) * C(SO)
      W(RR120B) = K(RR120B) * C(SCH2CO) * C(SCH3)
      W(RR121) = K(RR121) * C(STXC3H5) * C(SOH)
      W(RR122F) = K(RR122F) * C(STXC3H5) * C(SHO2)
      W(RR122B) = K(RR122B) * C(SOH) * C(SCH2CO) * C(SCH3)
      W(RR123F) = K(RR123F) * C(STXC3H5) * C(SHCO)
      W(RR123B) = K(RR123B) * C(SCO) * C(SC3H6)
      W(RR124F) = K(RR124F) * C(STXC3H5) * C(SCH3)
      W(RR124B) = K(RR124B) * C(SCH4) * C(SPXC3H4)
      W(RR125F) = K(RR125F) * C(SSXC3H5) * C(SH)
      W(RR125B) = K(RR125B) * C(SH2) * C(SPXC3H4)
      W(RR126F) = K(RR126F) * C(SSXC3H5) * C(SO)
      W(RR126B) = K(RR126B) * C(SHCO) * C(SC2H4)
      W(RR127) = K(RR127) * C(SSXC3H5) * C(SOH)
      W(RR128F) = K(RR128F) * C(SSXC3H5) * C(SHO2)
      W(RR128B) = K(RR128B) * C(SOH) * C(SHCO) * C(SC2H4)
      W(RR129F) = K(RR129F) * C(SSXC3H5) * C(SHCO)
      W(RR129B) = K(RR129B) * C(SCO) * C(SC3H6)
      W(RR130F) = K(RR130F) * C(SSXC3H5) * C(SCH3)
      W(RR130B) = K(RR130B) * C(SCH4) * C(SPXC3H4)
      W(RR131F) = K(RR131F) * C(STXC3H5) * C(SO2)
      W(RR131B) = K(RR131B) * C(SHO2) * C(SPXC3H4)
      W(RR132F) = K(RR132F) * C(SSXC3H5) * C(SO2)
      W(RR132B) = K(RR132B) * C(SHO2) * C(SPXC3H4)
      W(RR133) = K(RR133) * C(STXC3H5) * C(SO2)
      W(RR134) = K(RR134) * C(SSXC3H5) * C(SO2)
      W(RR135) = K(RR135) * C(STXC3H5) * C(SO2)
      W(RR136F) = K(RR136F) * C(SSXC3H5) * C(SO2)
      W(RR136B) = K(RR136B) * C(SHCO) * C(SCH3CHO)
      W(RR137F) = K(RR137F) * C(STXC3H5) * C(SO2)
      W(RR137B) = K(RR137B) * C(SHO2) * C(SAXC3H4)
      W(RR138) = K(RR138) * C(SC3H5O) * C(SO2)
      W(RR139F) = K(RR139F) * C(SC3H5O)
      W(RR139B) = K(RR139B) * C(SH) * C(SC2H3CHO)
      W(RR140) = K(RR140) * C(SC3H5O)
      W(RR141F) = K(RR141F) * C(SC3H6) * C(SH)
      W(RR141B) = K(RR141B) * C(SCH3) * C(SC2H4)
      W(RR142F) = K(RR142F) * C(SC3H6) * C(SO)
      W(RR142B) = K(RR142B) * C(SCH3) * C(SCH2CHO)
      W(RR143F) = K(RR143F) * C(SC3H6) * C(SO)
      W(RR143B) = K(RR143B) * C(SHCO) * C(SC2H5)
      W(RR144F) = K(RR144F) * C(SC3H6) * C(SH)
      W(RR144B) = K(RR144B) * C(SH2) * C(SAXC3H5)
      W(RR145F) = K(RR145F) * C(SC3H6) * C(SO)
      W(RR145B) = K(RR145B) * C(SOH) * C(SAXC3H5)
      W(RR146F) = K(RR146F) * C(SC3H6) * C(SOH)
      W(RR146B) = K(RR146B) * C(SH2O) * C(SAXC3H5)
      W(RR147F) = K(RR147F) * C(SC3H6) * C(SHO2)
      W(RR147B) = K(RR147B) * C(SH2O2) * C(SAXC3H5)
      W(RR148F) = K(RR148F) * C(SC3H6) * C(SCH3)
      W(RR148B) = K(RR148B) * C(SCH4) * C(SAXC3H5)
      W(RR149F) = K(RR149F) * C(SC3H6) * C(SH)
      W(RR149B) = K(RR149B) * C(SH2) * C(STXC3H5)
      W(RR150F) = K(RR150F) * C(SC3H6) * C(SO)
      W(RR150B) = K(RR150B) * C(SOH) * C(STXC3H5)
      W(RR151F) = K(RR151F) * C(SC3H6) * C(SOH)
      W(RR151B) = K(RR151B) * C(SH2O) * C(STXC3H5)
      W(RR152F) = K(RR152F) * C(SC3H6) * C(SCH3)
      W(RR152B) = K(RR152B) * C(SCH4) * C(STXC3H5)
      W(RR153F) = K(RR153F) * C(SC3H6) * C(SH)
      W(RR153B) = K(RR153B) * C(SH2) * C(SSXC3H5)
      W(RR154F) = K(RR154F) * C(SC3H6) * C(SO)
      W(RR154B) = K(RR154B) * C(SOH) * C(SSXC3H5)
      W(RR155F) = K(RR155F) * C(SC3H6) * C(SOH)
      W(RR155B) = K(RR155B) * C(SH2O) * C(SSXC3H5)
      W(RR156F) = K(RR156F) * C(SC3H6) * C(SCH3)
      W(RR156B) = K(RR156B) * C(SCH4) * C(SSXC3H5)
      W(RH01F) = K(RH01F) * C(SC4H) * C(SO2)
      W(RH01B) = K(RH01B) * C(SCO) * C(SCO) * C(SC2H)
      W(RH02F) = K(RH02F) * C(SC4H) * C(SH)
      W(RH02B) = K(RH02B) * C(SC4H2)
      W(RH03F) = K(RH03F) * C(SC4H2) * C(SH)
      W(RH03B) = K(RH03B) * C(SH2) * C(SC4H)
      W(RH04F) = K(RH04F) * C(SC4H2) * C(SH2)
      W(RH04B) = K(RH04B) * C(SC4H4)
      W(RH05) = K(RH05) * C(SC4H2) * C(SC4H2)
      W(RH06F) = K(RH06F) * C(SC4H2) * C(SC4H2)
      W(RH06B) = K(RH06B) * C(SH2) * C(SC8H2)
      W(RH07F) = K(RH07F) * C(SC4H2) * C(SO2)
      W(RH07B) = K(RH07B) * C(SHCCO) * C(SHCCO)
      W(RH08F) = K(RH08F) * C(SC4H2) * C(SO)
      W(RH08B) = K(RH08B) * C(SCO) * C(SC3H2)
      W(RH09F) = K(RH09F) * C(SC4H2) * C(SH)
      W(RH09B) = K(RH09B) * C(SIXC4H3)
      W(RH10F) = K(RH10F) * C(SC4H2) * C(SH)
      W(RH10B) = K(RH10B) * C(SNXC4H3)
      W(RH11F) = K(RH11F) * C(SC4H2) * C(SOH)
      W(RH11B) = K(RH11B) * C(SH2O) * C(SC4H)
      W(RH12F) = K(RH12F) * C(SC4H2) * C(SOH)
      W(RH12B) = K(RH12B) * C(SCO) * C(SC3H3)
      W(RH13F) = K(RH13F) * C(SNXC4H3)
      W(RH13B) = K(RH13B) * C(SIXC4H3)
      W(RH14F) = K(RH14F) * C(SNXC4H3) * C(SH)
      W(RH14B) = K(RH14B) * C(SH) * C(SIXC4H3)
      W(RH15F) = K(RH15F) * C(SNXC4H3) * C(SH)
      W(RH15B) = K(RH15B) * C(SC4H4)
      W(RH16F) = K(RH16F) * C(SIXC4H3) * C(SH)
      W(RH16B) = K(RH16B) * C(SC4H4)
      W(RH17F) = K(RH17F) * C(SNXC4H3) * C(SH)
      W(RH17B) = K(RH17B) * C(SC2H2) * C(SC2H2)
      W(RH18F) = K(RH18F) * C(SIXC4H3) * C(SH)
      W(RH18B) = K(RH18B) * C(SC2H2) * C(SC2H2)
      W(RH19F) = K(RH19F) * C(SNXC4H3) * C(SH)
      W(RH19B) = K(RH19B) * C(SH2) * C(SC4H2)
      W(RH20F) = K(RH20F) * C(SIXC4H3) * C(SH)
      W(RH20B) = K(RH20B) * C(SH2) * C(SC4H2)
      W(RH21F) = K(RH21F) * C(SNXC4H3) * C(SOH)
      W(RH21B) = K(RH21B) * C(SH2O) * C(SC4H2)
      W(RH22F) = K(RH22F) * C(SIXC4H3) * C(SOH)
      W(RH22B) = K(RH22B) * C(SH2O) * C(SC4H2)
      W(RH23F) = K(RH23F) * C(SNXC4H3) * C(SO2)
      W(RH23B) = K(RH23B) * C(SHO2) * C(SC4H2)
      W(RH24F) = K(RH24F) * C(SIXC4H3) * C(SO2)
      W(RH24B) = K(RH24B) * C(SHO2) * C(SC4H2)
      W(RH25F) = K(RH25F) * C(SIXC4H3) * C(SO)
      W(RH25B) = K(RH25B) * C(SC2H) * C(SCH2CO)
      W(RH26F) = K(RH26F) * C(SIXC4H3) * C(SO2)
      W(RH26B) = K(RH26B) * C(SCH2CO) * C(SHCCO)
      W(RH27F) = K(RH27F) * C(SIXC4H3) * C(SO2)
      W(RH27B) = K(RH27B) * C(SCO) * C(SC2H2) * C(SHCO)
      W(RH28F) = K(RH28F) * C(SC4H4) * C(SH)
      W(RH28B) = K(RH28B) * C(SH2) * C(SNXC4H3)
      W(RH29F) = K(RH29F) * C(SC4H4) * C(SH)
      W(RH29B) = K(RH29B) * C(SH2) * C(SIXC4H3)
      W(RH30F) = K(RH30F) * C(SC4H4) * C(SOH)
      W(RH30B) = K(RH30B) * C(SH2O) * C(SNXC4H3)
      W(RH31F) = K(RH31F) * C(SC4H4) * C(SOH)
      W(RH31B) = K(RH31B) * C(SH2O) * C(SIXC4H3)
      W(RH32F) = K(RH32F) * C(SC4H4) * C(SCH3)
      W(RH32B) = K(RH32B) * C(SCH4) * C(SNXC4H3)
      W(RH33F) = K(RH33F) * C(SC4H4) * C(SCH3)
      W(RH33B) = K(RH33B) * C(SCH4) * C(SIXC4H3)
      W(RH34F) = K(RH34F) * C(SC4H4) * C(SO)
      W(RH34B) = K(RH34B) * C(SCO) * C(SAXC3H4)
      W(RH35F) = K(RH35F) * C(SC4H4) * C(SO)
      W(RH35B) = K(RH35B) * C(SCH2O) * C(SC3H2)
      W(RH36F) = K(RH36F) * C(SC4H4) * C(SO)
      W(RH36B) = K(RH36B) * C(SHCO) * C(SC3H3)
      W(RH37F) = K(RH37F) * C(SC4H2) * C(SC2H)
      W(RH37B) = K(RH37B) * C(SH) * C(SC6H2)
      W(RH38F) = K(RH38F) * C(SC2H2) * C(SC4H)
      W(RH38B) = K(RH38B) * C(SH) * C(SC6H2)
      W(RH39F) = K(RH39F) * C(SC6H2) * C(SC2H)
      W(RH39B) = K(RH39B) * C(SH) * C(SC8H2)
      W(RH40F) = K(RH40F) * C(SC4H2) * C(SC4H)
      W(RH40B) = K(RH40B) * C(SH) * C(SC8H2)
      W(RB00F) = K(RB00F) * C(SH2C2) * C(SC2H4)
      W(RB00B) = K(RB00B) * C(SC4H6)
      W(RB01F) = K(RB01F) * C(SH2C2) * C(SC2H2)
      W(RB01B) = K(RB01B) * C(SC4H4)
      W(RB02F) = K(RB02F) * C(SC2H3) * C(SC2H2)
      W(RB02B) = K(RB02B) * C(SNXC4H5)
      W(RB04F) = K(RB04F) * C(SC2H3) * C(SC2H3)
      W(RB04B) = K(RB04B) * C(SC4H6)
      W(RB05F) = K(RB05F) * C(SC2H3) * C(SC2H3)
      W(RB05B) = K(RB05B) * C(SH) * C(SIXC4H5)
      W(RB06F) = K(RB06F) * C(SC2H3) * C(SC2H3)
      W(RB06B) = K(RB06B) * C(SH) * C(SNXC4H5)
      W(RB07F) = K(RB07F) * C(SC2H3) * C(SC2H3)
      W(RB07B) = K(RB07B) * C(SC2H4) * C(SC2H2)
      W(RB08F) = K(RB08F) * C(SC3H3) * C(SCH3)
      W(RB08B) = K(RB08B) * C(SC4H6)
      W(RB09F) = K(RB09F) * C(SC3H6) * C(SC2H3)
      W(RB09B) = K(RB09B) * C(SCH3) * C(SC4H6)
      W(RB10F) = K(RB10F) * C(SC4H6)
      W(RB10B) = K(RB10B) * C(SH) * C(SIXC4H5)
      W(RB11F) = K(RB11F) * C(SC4H6)
      W(RB11B) = K(RB11B) * C(SH) * C(SNXC4H5)
      W(RB12F) = K(RB12F) * C(SC4H6)
      W(RB12B) = K(RB12B) * C(SH2) * C(SC4H4)
      W(RB14F) = K(RB14F) * C(SPXC3H4) * C(SCH3)
      W(RB14B) = K(RB14B) * C(SH) * C(SC4H6)
      W(RB15F) = K(RB15F) * C(SAXC3H4) * C(SCH3)
      W(RB15B) = K(RB15B) * C(SH) * C(SC4H6)
      W(RB16F) = K(RB16F) * C(SC4H6) * C(SH)
      W(RB16B) = K(RB16B) * C(SH2) * C(SNXC4H5)
      W(RB17F) = K(RB17F) * C(SC4H6) * C(SH)
      W(RB17B) = K(RB17B) * C(SH2) * C(SIXC4H5)
      W(RB18F) = K(RB18F) * C(SNXC4H5) * C(SOH)
      W(RB18B) = K(RB18B) * C(SO) * C(SC4H6)
      W(RB19F) = K(RB19F) * C(SC4H6) * C(SO)
      W(RB19B) = K(RB19B) * C(SOH) * C(SIXC4H5)
      W(RB20F) = K(RB20F) * C(SC4H6) * C(SOH)
      W(RB20B) = K(RB20B) * C(SH2O) * C(SNXC4H5)
      W(RB21F) = K(RB21F) * C(SC4H6) * C(SOH)
      W(RB21B) = K(RB21B) * C(SH2O) * C(SIXC4H5)
      W(RB22F) = K(RB22F) * C(SC4H6) * C(SCH3)
      W(RB22B) = K(RB22B) * C(SCH4) * C(SNXC4H5)
      W(RB23F) = K(RB23F) * C(SC4H6) * C(SCH3)
      W(RB23B) = K(RB23B) * C(SCH4) * C(SIXC4H5)
      W(RB24F) = K(RB24F) * C(SC4H6) * C(SC2H3)
      W(RB24B) = K(RB24B) * C(SC2H4) * C(SNXC4H5)
      W(RB25F) = K(RB25F) * C(SC4H6) * C(SC2H3)
      W(RB25B) = K(RB25B) * C(SC2H4) * C(SIXC4H5)
      W(RB28) = K(RB28) * C(SC4H6) * C(SO)
      W(RB29F) = K(RB29F) * C(SC4H6) * C(SO)
      W(RB29B) = K(RB29B) * C(SCH2O) * C(SPXC3H4)
      W(RB30F) = K(RB30F) * C(SC4H6) * C(SO)
      W(RB30B) = K(RB30B) * C(SHCO) * C(SAXC3H5)
      W(RB31F) = K(RB31F) * C(SC4H6) * C(SOH)
      W(RB31B) = K(RB31B) * C(SCH2O) * C(SAXC3H5)
      W(RB32F) = K(RB32F) * C(SC4H4) * C(SH)
      W(RB32B) = K(RB32B) * C(SNXC4H5)
      W(RB33F) = K(RB33F) * C(SC4H4) * C(SH)
      W(RB33B) = K(RB33B) * C(SIXC4H5)
      W(RB34F) = K(RB34F) * C(SNXC4H5)
      W(RB34B) = K(RB34B) * C(SIXC4H5)
      W(RB35F) = K(RB35F) * C(SNXC4H5) * C(SH)
      W(RB35B) = K(RB35B) * C(SH) * C(SIXC4H5)
      W(RB36F) = K(RB36F) * C(SNXC4H5) * C(SH)
      W(RB36B) = K(RB36B) * C(SH2) * C(SC4H4)
      W(RB37F) = K(RB37F) * C(SNXC4H5) * C(SOH)
      W(RB37B) = K(RB37B) * C(SH2O) * C(SC4H4)
      W(RB38F) = K(RB38F) * C(SNXC4H5) * C(SHCO)
      W(RB38B) = K(RB38B) * C(SCO) * C(SC4H6)
      W(RB39F) = K(RB39F) * C(SNXC4H5) * C(SH2O2)
      W(RB39B) = K(RB39B) * C(SHO2) * C(SC4H6)
      W(RB40F) = K(RB40F) * C(SNXC4H5) * C(SHO2)
      W(RB40B) = K(RB40B) * C(SO2) * C(SC4H6)
      W(RB99F) = K(RB99F) * C(SNXC4H5) * C(SO)
      W(RB99B) = K(RB99B) * C(SCO) * C(SAXC3H5)
      W(RB41F) = K(RB41F) * C(SNXC4H5) * C(SO2)
      W(RB41B) = K(RB41B) * C(SHO2) * C(SC4H4)
      W(RB42) = K(RB42) * C(SNXC4H5) * C(SO2)
      W(RB43F) = K(RB43F) * C(SNXC4H5) * C(SO2)
      W(RB43B) = K(RB43B) * C(SC2H3CHO) * C(SHCO)
      W(RB44F) = K(RB44F) * C(SIXC4H5) * C(SH)
      W(RB44B) = K(RB44B) * C(SH2) * C(SC4H4)
      W(RB45F) = K(RB45F) * C(SIXC4H5) * C(SH)
      W(RB45B) = K(RB45B) * C(SCH3) * C(SC3H3)
      W(RB46F) = K(RB46F) * C(SIXC4H5) * C(SOH)
      W(RB46B) = K(RB46B) * C(SH2O) * C(SC4H4)
      W(RB47F) = K(RB47F) * C(SIXC4H5) * C(SHCO)
      W(RB47B) = K(RB47B) * C(SCO) * C(SC4H6)
      W(RB48F) = K(RB48F) * C(SIXC4H5) * C(SHO2)
      W(RB48B) = K(RB48B) * C(SO2) * C(SC4H6)
      W(RB49F) = K(RB49F) * C(SIXC4H5) * C(SH2O2)
      W(RB49B) = K(RB49B) * C(SHO2) * C(SC4H6)
      W(RB50F) = K(RB50F) * C(SIXC4H5) * C(SO2)
      W(RB50B) = K(RB50B) * C(SCH2CHO) * C(SCH2CO)
      W(RB51F) = K(RB51F) * C(SIXC4H5) * C(SO)
      W(RB51B) = K(RB51B) * C(SCH2O) * C(SC3H3)
      W(RB52F) = K(RB52F) * C(SNXC4H5) * C(SC2H3)
      W(RB52B) = K(RB52B) * C(SH2) * C(SA1XC6H6)
      W(RHP00) = K(RHP00) * C(SNXC7H16)
      W(RHP01) = K(RHP01) * C(SNXC7H16)
      W(RHP02) = K(RHP02) * C(SNXC7H16) * C(SH)
      W(RHP03) = K(RHP03) * C(SNXC7H16) * C(SO)
      W(RHP04) = K(RHP04) * C(SNXC7H16) * C(SOH)
      W(RHP05) = K(RHP05) * C(SNXC7H16) * C(SO2)
      W(RHP06) = K(RHP06) * C(SNXC7H16) * C(SHO2)
      W(RHP08) = K(RHP08) * C(SC7H15) * C(SH)
      W(RHP09) = K(RHP09) * C(SC7H15) * C(SHO2)
      W(RHP10) = K(RHP10) * C(SC7H15)
      W(RHP11) = K(RHP11) * C(SC7H15)
      W(RHP12) = K(RHP12) * C(SC7H15)
      W(RHP13) = K(RHP13) * C(SC7H15)
      W(RHP14) = K(RHP14) * C(SC7H15)
      W(RHP15) = K(RHP15) * C(SC7H14) * C(SH)
      W(RHP16) = K(RHP16) * C(SC7H15) * C(SHO2)
      W(RHP17) = K(RHP17) * C(SC7H15) * C(SCH3O2)
      W(RHP18) = K(RHP18) * C(SC7H15O)
      W(RHP19) = K(RHP19) * C(SC7H15O)
      W(RHP20) = K(RHP20) * C(SC7H15O)
      W(RHP21) = K(RHP21) * C(SC7H14)
      W(RHP22) = K(RHP22) * C(SC7H14)
      W(RHP23F) = K(RHP23F) * C(SC7H14) * C(SH)
      W(RHP23B) = K(RHP23B) * C(SH2) * C(SC7H13)
      W(RHP24) = K(RHP24) * C(SC7H14) * C(SOH)
      W(RHP25) = K(RHP25) * C(SC7H14) * C(SOH)
      W(RHP26) = K(RHP26) * C(SC7H14) * C(SOH)
      W(RHP27) = K(RHP27) * C(SC7H13)
      W(RHP28) = K(RHP28) * C(SC7H13)
      W(RHP29) = K(RHP29) * C(SC5H11)
      W(RHP30) = K(RHP30) * C(SC5H11)
      W(RHP31) = K(RHP31) * C(SC5H10) * C(SH)
      W(RHP32) = K(RHP32) * C(SC5H11)
      W(RHP33) = K(RHP33) * C(SC5H10)
      W(RHP34) = K(RHP34) * C(SC5H10) * C(SH)
      W(RHP35) = K(RHP35) * C(SC5H10) * C(SO)
      W(RHP36) = K(RHP36) * C(SC5H10) * C(SOH)
      W(RHP37) = K(RHP37) * C(SC5H9)
      W(RHP38) = K(RHP38) * C(SC5H9)
      W(RHP39) = K(RHP39) * C(SPXC4H9)
      W(RHP40) = K(RHP40) * C(SPXC4H8) * C(SH)
      W(RHP41) = K(RHP41) * C(SPXC4H9)
      W(RHP42) = K(RHP42) * C(SPXC4H9)
      W(RHP43) = K(RHP43) * C(SPXC4H8)
      W(RHP44) = K(RHP44) * C(SPXC4H8) * C(SH)
      W(RHP45) = K(RHP45) * C(SPXC4H8) * C(SOH)
      W(RHP46) = K(RHP46) * C(SPXC4H8) * C(SO2)
      W(RHP47) = K(RHP47) * C(SPXC4H8) * C(SHO2)
      W(RHP48) = K(RHP48) * C(SPXC4H8) * C(SCH3)
      W(RHP49) = K(RHP49) * C(SC4H7) * C(SH)
      W(RHP50) = K(RHP50) * C(SC4H7) * C(SHO2)
      W(RHP51) = K(RHP51) * C(SPXC4H8) * C(SO)
      W(RHP52) = K(RHP52) * C(SPXC4H8) * C(SO)
      W(RHP53) = K(RHP53) * C(SPXC4H8) * C(SO)
      W(RHP54) = K(RHP54) * C(SPXC4H8) * C(SO)
      W(RHP55) = K(RHP55) * C(SPXC4H8) * C(SOH)
      W(RHP58) = K(RHP58) * C(SPXC4H8) * C(SOH)
      W(RHP56) = K(RHP56) * C(SPXC4H8) * C(SOH)
      W(RHP57) = K(RHP57) * C(SPXC4H8) * C(SOH)
      W(RHP59) = K(RHP59) * C(SC3H7CHO) * C(SH)
      W(RHP60) = K(RHP60) * C(SC3H7CHO) * C(SOH)
      W(RHP61) = K(RHP61) * C(SC3H7CHO) * C(SO2)
      W(RHP62) = K(RHP62) * C(SC3H7CHO) * C(SCH3)
      W(RHP63) = K(RHP63) * C(SC3H7CHO) * C(SHO2)
      W(RHP64F) = K(RHP64F) * C(SC4H7)
      W(RHP64B) = K(RHP64B) * C(SH) * C(SC4H6)
      W(RHP65F) = K(RHP65F) * C(SC2H4) * C(SC2H3)
      W(RHP65B) = K(RHP65B) * C(SC4H7)
      W(RHP67) = K(RHP67) * C(SC4H7) * C(SH)
      W(RHP68) = K(RHP68) * C(SC4H7) * C(SO2)
      W(RHP69) = K(RHP69) * C(SC4H7) * C(SCH3)
      W(RHP70) = K(RHP70) * C(SC4H7) * C(SC2H5)
      W(RHP71) = K(RHP71) * C(SC4H7) * C(SAXC3H5)
      W(RHP72) = K(RHP72) * C(SC4H7) * C(SHO2)
      W(RHP73) = K(RHP73) * C(SC4H7) * C(SCH3O2)
      W(RHP74) = K(RHP74) * C(SC4H7O)
      W(RHP75) = K(RHP75) * C(SC4H7O)
      W(RHP76) = K(RHP76) * C(SNXC3H7) * C(SHO2)
      W(RHP77) = K(RHP77) * C(SNXC3H7) * C(SCH3O2)
      W(RHP78) = K(RHP78) * C(SNXC3H7O)
      W(RHP79) = K(RHP79) * C(SNXC3H7O)
      W(RHP86) = K(RHP86) * C(SC4H6) * C(SO)
      W(RHP87) = K(RHP87) * C(SC4H6) * C(SOH)
      W(RHP88) = K(RHP88) * C(SC4H6) * C(SOH)
      W(RIC00) = K(RIC00) * C(SIXC8H18)
      W(RIC01) = K(RIC01) * C(SIXC8H18)
      W(RIC02) = K(RIC02) * C(SIXC8H18)
      W(RIC03) = K(RIC03) * C(SIXC8H18) * C(SH)
      W(RIC04) = K(RIC04) * C(SIXC8H18) * C(SO)
      W(RIC05) = K(RIC05) * C(SIXC8H18) * C(SOH)
      W(RIC06) = K(RIC06) * C(SIXC8H18) * C(SO2)
      W(RIC07) = K(RIC07) * C(SIXC8H18) * C(SCH3)
      W(RIC08) = K(RIC08) * C(SIXC8H18) * C(SHO2)
      W(RIC09) = K(RIC09) * C(SIXC8H18) * C(SCH3O2)
      W(RIC10) = K(RIC10) * C(SCXC8H17)
      W(RIC11) = K(RIC11) * C(SCXC8H17)
      W(RIC12) = K(RIC12) * C(SCXC8H17)
      W(RIC13) = K(RIC13) * C(SCXC8H17) * C(SHO2)
      W(RIC14) = K(RIC14) * C(SCXC8H17) * C(SCH3O2)
      W(RIC15) = K(RIC15) * C(SDXC8H17O)
      W(RIC16) = K(RIC16) * C(SDXC8H17O)
      W(RIC17) = K(RIC17) * C(SDXC8H17O)
      W(RIC18) = K(RIC18) * C(SYXC7H14) * C(SH)
      W(RIC19) = K(RIC19) * C(SYXC7H15)
      W(RIC20) = K(RIC20) * C(SYXC7H15)
      W(RIC21) = K(RIC21) * C(SYXC7H15)
      W(RIC22) = K(RIC22) * C(SYXC7H14)
      W(RIC23) = K(RIC23) * C(SYXC7H14)
      W(RIC24) = K(RIC24) * C(SYXC7H14) * C(SH)
      W(RIC25) = K(RIC25) * C(SXXC7H13) * C(SH2)
      W(RIC26) = K(RIC26) * C(SYXC7H14) * C(SOH)
      W(RIC27) = K(RIC27) * C(SXXC7H13) * C(SH2O)
      W(RIC28) = K(RIC28) * C(SXXC7H13) * C(SHO2)
      W(RIC29) = K(RIC29) * C(STXC4H9O)
      W(RIC30) = K(RIC30) * C(STXC4H9O)
      W(RIC31) = K(RIC31) * C(STXC4H9O)
      W(RIC32) = K(RIC32) * C(STXC4H9)
      W(RIC33) = K(RIC33) * C(STXC4H9)
      W(RIC34) = K(RIC34) * C(SIXC4H8) * C(SH)
      W(RIC35) = K(RIC35) * C(STXC4H9) * C(SHO2)
      W(RIC36) = K(RIC36) * C(STXC4H9) * C(SCH3O2)
      W(RIC37) = K(RIC37) * C(SIXC4H8)
      W(RIC38F) = K(RIC38F) * C(SIXC4H8)
      W(RIC38B) = K(RIC38B) * C(SH) * C(SIXC4H7)
      W(RIC40) = K(RIC40) * C(SIXC4H8) * C(SH)
      W(RIC41F) = K(RIC41F) * C(SIXC4H8) * C(SH)
      W(RIC41B) = K(RIC41B) * C(SH2) * C(SIXC4H7)
      W(RIC42) = K(RIC42) * C(SIXC4H8) * C(SO)
      W(RIC43) = K(RIC43) * C(SIXC4H8) * C(SOH)
      W(RIC44) = K(RIC44) * C(SIXC4H8) * C(SCH3)
      W(RIC45) = K(RIC45) * C(SIXC4H8) * C(SHO2)
      W(RIC46) = K(RIC46) * C(SIXC4H8) * C(SCH3O2)
      W(RIC47) = K(RIC47) * C(SIXC4H8) * C(SO)
      W(RIC48) = K(RIC48) * C(SIXC4H8) * C(SO)
      W(RIC49) = K(RIC49) * C(SIXC4H7O)
      W(RIC50) = K(RIC50) * C(SIXC4H7O)
      W(RIC51) = K(RIC51) * C(SIXC4H7O) * C(SO2)
      W(RIC52) = K(RIC52) * C(SIXC4H7)
      W(RIC57) = K(RIC57) * C(SIXC4H7) * C(SHO2)
      W(RIC53) = K(RIC53) * C(SIXC4H7) * C(SO)
      W(RIC54) = K(RIC54) * C(SIXC4H7) * C(SO2)
      W(RIC55) = K(RIC55) * C(SIXC4H7) * C(SO2)
      W(RIC56) = K(RIC56) * C(SIXC4H7) * C(SO2)
      W(RIC58) = K(RIC58) * C(SIXC4H7) * C(SCH3O2)
      W(RIC59) = K(RIC59) * C(SIXC3H5CH) * C(SH)
      W(RIC60) = K(RIC60) * C(SIXC3H5CH) * C(SO)
      W(RIC61) = K(RIC61) * C(SIXC3H5CH) * C(SOH)
      W(RIC62) = K(RIC62) * C(SIXC3H5CH) * C(SHO2)
      W(RIC63) = K(RIC63) * C(SIXC3H5CH) * C(SCH3)
      W(RIC64) = K(RIC64) * C(SIXC3H7)
      W(RIC65) = K(RIC65) * C(SIXC3H7)
      W(RIC66) = K(RIC66) * C(SC3H6) * C(SH)
      W(RIC67) = K(RIC67) * C(SC2H4) * C(SCH3)
      W(RIC68) = K(RIC68) * C(SIXC3H7) * C(SO2)
      W(RIC69) = K(RIC69) * C(SCH3COCH3) * C(SH)
      W(RIC70) = K(RIC70) * C(SCH3COCH3) * C(SO)
      W(RIC71) = K(RIC71) * C(SCH3COCH3) * C(SOH)
      W(RIC72) = K(RIC72) * C(SCH3COCH3) * C(SHO2)
      W(RP000F) = K(RP000F) * C(SC5H4CH2) * C(SH)
      W(RP000B) = K(RP000B) * C(SH) * C(SA1XC6H6)
      W(RP001F) = K(RP001F) * C(SNXC4H5) * C(SC2H2)
      W(RP001B) = K(RP001B) * C(SH) * C(SC5H4CH2)
      W(RP002F) = K(RP002F) * C(SIXC4H5) * C(SC2H2)
      W(RP002B) = K(RP002B) * C(SH) * C(SC5H4CH2)
      W(RP003F) = K(RP003F) * C(SC5H4CH2)
      W(RP003B) = K(RP003B) * C(SA1XC6H6)
      W(RP004F) = K(RP004F) * C(SC5H4CH2)
      W(RP004B) = K(RP004B) * C(SH) * C(SA1XXC6H5)
      W(RP005F) = K(RP005F) * C(SNXC4H3) * C(SC2H2)
      W(RP005B) = K(RP005B) * C(SA1XXC6H5)
      W(RP006F) = K(RP006F) * C(SNXC4H5) * C(SC2H2)
      W(RP006B) = K(RP006B) * C(SH) * C(SA1XC6H6)
      W(RP007F) = K(RP007F) * C(SIXC4H5) * C(SC2H2)
      W(RP007B) = K(RP007B) * C(SH) * C(SA1XC6H6)
      W(RP008) = K(RP008) * C(SAXC3H5) * C(SC3H3)
      W(RP009F) = K(RP009F) * C(SC3H3) * C(SC3H3)
      W(RP009B) = K(RP009B) * C(SC5H4CH2)
      W(RP010F) = K(RP010F) * C(SC3H3) * C(SC3H3)
      W(RP010B) = K(RP010B) * C(SA1XC6H6)
      W(RP011F) = K(RP011F) * C(SC3H3) * C(SC3H3)
      W(RP011B) = K(RP011B) * C(SH) * C(SA1XXC6H5)
      W(RK012F) = K(RK012F) * C(SA1XXC6H5) * C(SC2H2)
      W(RK012B) = K(RK012B) * C(SA1C2H2XC)
      W(RP013F) = K(RP013F) * C(SA1XXC6H5) * C(SC2H3)
      W(RP013B) = K(RP013B) * C(SA1C2H3XC)
      W(RP014F) = K(RP014F) * C(SA1XC6H6) * C(SC2H3)
      W(RP014B) = K(RP014B) * C(SH) * C(SA1C2H3XC)
      W(RP015F) = K(RP015F) * C(SA1XXC6H5) * C(SC2H4)
      W(RP015B) = K(RP015B) * C(SC2H3) * C(SA1XC6H6)
      W(RP016F) = K(RP016F) * C(SA1C2HXC8)
      W(RP016B) = K(RP016B) * C(SH) * C(SA1C2HYXC)
      W(RK017F) = K(RK017F) * C(SA1C2HXC8) * C(SH)
      W(RK017B) = K(RK017B) * C(SH2) * C(SA1C2HYXC)
      W(RP018F) = K(RP018F) * C(SA1C2HXC8) * C(SOH)
      W(RP018B) = K(RP018B) * C(SH2O) * C(SA1C2HYXC)
      W(RK019F) = K(RK019F) * C(SA1C2H2XC)
      W(RK019B) = K(RK019B) * C(SA1C2H3YX)
      W(RK020F) = K(RK020F) * C(SA1C2H2XC)
      W(RK020B) = K(RK020B) * C(SH) * C(SA1C2HXC8)
      W(RK021F) = K(RK021F) * C(SA1C2H2XC) * C(SH)
      W(RK021B) = K(RK021B) * C(SH2) * C(SA1C2HXC8)
      W(RP022F) = K(RP022F) * C(SA1C2H2XC) * C(SOH)
      W(RP022B) = K(RP022B) * C(SH2O) * C(SA1C2HXC8)
      W(RP023F) = K(RP023F) * C(SA1C2H3XC)
      W(RP023B) = K(RP023B) * C(SH) * C(SA1C2H3YX)
      W(RK024F) = K(RK024F) * C(SA1C2H3XC) * C(SH)
      W(RK024B) = K(RK024B) * C(SH2) * C(SA1C2H3YX)
      W(RP025F) = K(RP025F) * C(SA1C2H3XC) * C(SOH)
      W(RP025B) = K(RP025B) * C(SH2O) * C(SA1C2H3YX)
      W(RP026F) = K(RP026F) * C(SA1C2H3XC)
      W(RP026B) = K(RP026B) * C(SH) * C(SA1C2H2XC)
      W(RP027F) = K(RP027F) * C(SA1C2H3XC) * C(SH)
      W(RP027B) = K(RP027B) * C(SH2) * C(SA1C2H2XC)
      W(RP028F) = K(RP028F) * C(SA1C2H3XC) * C(SOH)
      W(RP028B) = K(RP028B) * C(SH2O) * C(SA1C2H2XC)
      W(RK100F) = K(RK100F) * C(SA1C2HYXC) * C(SC2H2)
      W(RK100B) = K(RK100B) * C(SA2XXC10H)
      W(RK102F) = K(RK102F) * C(SA1C2H3YX) * C(SC2H2)
      W(RK102B) = K(RK102B) * C(SH) * C(SA2XC10H8)
      W(RP104F) = K(RP104F) * C(SA1C2H2XC) * C(SC2H2)
      W(RP104B) = K(RP104B) * C(SH) * C(SA2XC10H8)
      W(RP105F) = K(RP105F) * C(SA1C2HXC8) * C(SC2H3)
      W(RP105B) = K(RP105B) * C(SH) * C(SA2XC10H8)
      W(RP106F) = K(RP106F) * C(SA1C2HYXC) * C(SC2H4)
      W(RP106B) = K(RP106B) * C(SH) * C(SA2XC10H8)
      W(RP107F) = K(RP107F) * C(SA1XXC6H5) * C(SC4H4)
      W(RP107B) = K(RP107B) * C(SH) * C(SA2XC10H8)
      W(RP108F) = K(RP108F) * C(SA2XC10H8)
      W(RP108B) = K(RP108B) * C(SH) * C(SA2XXC10H)
      W(RK109F) = K(RK109F) * C(SA2XC10H8) * C(SH)
      W(RK109B) = K(RK109B) * C(SH2) * C(SA2XXC10H)
      W(RK110F) = K(RK110F) * C(SA2XC10H8) * C(SOH)
      W(RK110B) = K(RK110B) * C(SH2O) * C(SA2XXC10H)
      W(RP111F) = K(RP111F) * C(SA2XC10H8)
      W(RP111B) = K(RP111B) * C(SH) * C(SA2YXC10H)
      W(RK112F) = K(RK112F) * C(SA2XC10H8) * C(SH)
      W(RK112B) = K(RK112B) * C(SH2) * C(SA2YXC10H)
      W(RK113F) = K(RK113F) * C(SA2XC10H8) * C(SOH)
      W(RK113B) = K(RK113B) * C(SH2O) * C(SA2YXC10H)
      W(RK114F) = K(RK114F) * C(SA2XXC10H) * C(SC2H2)
      W(RK114B) = K(RK114B) * C(SA2C2H2AX)
      W(RK115F) = K(RK115F) * C(SA2YXC10H) * C(SC2H2)
      W(RK115B) = K(RK115B) * C(SA2C2H2BX)
      W(RP116F) = K(RP116F) * C(SA2XXC10H) * C(SC2H3)
      W(RP116B) = K(RP116B) * C(SH) * C(SA2C2H2AX)
      W(RP117F) = K(RP117F) * C(SA2YXC10H) * C(SC2H3)
      W(RP117B) = K(RP117B) * C(SH) * C(SA2C2H2BX)
      W(RP118F) = K(RP118F) * C(SA2XC10H8) * C(SC2H3)
      W(RP118B) = K(RP118B) * C(SH2) * C(SA2C2H2AX)
      W(RP119F) = K(RP119F) * C(SA2XC10H8) * C(SC2H3)
      W(RP119B) = K(RP119B) * C(SH2) * C(SA2C2H2BX)
      W(RP120F) = K(RP120F) * C(SA2XXC10H) * C(SC2H4)
      W(RP120B) = K(RP120B) * C(SH2) * C(SA2C2H2AX)
      W(RP121F) = K(RP121F) * C(SA2YXC10H) * C(SC2H4)
      W(RP121B) = K(RP121B) * C(SH2) * C(SA2C2H2BX)
      W(RK122F) = K(RK122F) * C(SA2C2H2AX)
      W(RK122B) = K(RK122B) * C(SH) * C(SA2C2HAXC)
      W(RK123F) = K(RK123F) * C(SA2C2H2AX) * C(SH)
      W(RK123B) = K(RK123B) * C(SH2) * C(SA2C2HAXC)
      W(RP124F) = K(RP124F) * C(SA2C2H2AX) * C(SOH)
      W(RP124B) = K(RP124B) * C(SH2O) * C(SA2C2HAXC)
      W(RK125F) = K(RK125F) * C(SA2C2H2BX)
      W(RK125B) = K(RK125B) * C(SH) * C(SA2C2HBXC)
      W(RK126F) = K(RK126F) * C(SA2C2H2BX) * C(SH)
      W(RK126B) = K(RK126B) * C(SH2) * C(SA2C2HBXC)
      W(RP127F) = K(RP127F) * C(SA2C2H2BX) * C(SOH)
      W(RP127B) = K(RP127B) * C(SH2O) * C(SA2C2HBXC)
      W(RP128F) = K(RP128F) * C(SA2C2HAXC)
      W(RP128B) = K(RP128B) * C(SH) * C(SA2C2HAYX)
      W(RK129F) = K(RK129F) * C(SA2C2HAXC) * C(SH)
      W(RK129B) = K(RK129B) * C(SH2) * C(SA2C2HAYX)
      W(RP130F) = K(RP130F) * C(SA2C2HAXC) * C(SOH)
      W(RP130B) = K(RP130B) * C(SH2O) * C(SA2C2HAYX)
      W(RP131F) = K(RP131F) * C(SA2C2HBXC)
      W(RP131B) = K(RP131B) * C(SH) * C(SA2C2HBYX)
      W(RK132F) = K(RK132F) * C(SA2C2HBXC) * C(SH)
      W(RK132B) = K(RK132B) * C(SH2) * C(SA2C2HBYX)
      W(RP133F) = K(RP133F) * C(SA2C2HBXC) * C(SOH)
      W(RP133B) = K(RP133B) * C(SH2O) * C(SA2C2HBYX)
      W(RK200F) = K(RK200F) * C(SA2C2H2AX)
      W(RK200B) = K(RK200B) * C(SH) * C(SA2R5XC12)
      W(RP201F) = K(RP201F) * C(SA2C2HAXC) * C(SH)
      W(RP201B) = K(RP201B) * C(SH) * C(SA2R5XC12)
      W(RP202F) = K(RP202F) * C(SA2R5XC12)
      W(RP202B) = K(RP202B) * C(SH) * C(SA2R5XXC1)
      W(RK203F) = K(RK203F) * C(SA2R5XC12) * C(SH)
      W(RK203B) = K(RK203B) * C(SH2) * C(SA2R5XXC1)
      W(RK204F) = K(RK204F) * C(SA2R5XC12) * C(SOH)
      W(RK204B) = K(RK204B) * C(SH2O) * C(SA2R5XXC1)
      W(RK205F) = K(RK205F) * C(SA2R5XXC1) * C(SC2H2)
      W(RK205B) = K(RK205B) * C(SA2R5C2H2)
      W(RP206F) = K(RP206F) * C(SA2R5XXC1) * C(SC2H3)
      W(RP206B) = K(RP206B) * C(SH) * C(SA2R5C2H2)
      W(RP207F) = K(RP207F) * C(SA2R5XC12) * C(SC2H3)
      W(RP207B) = K(RP207B) * C(SH2) * C(SA2R5C2H2)
      W(RP208F) = K(RP208F) * C(SA2R5XXC1) * C(SC2H4)
      W(RP208B) = K(RP208B) * C(SH2) * C(SA2R5C2H2)
      W(RP209F) = K(RP209F) * C(SA2R5C2HX)
      W(RP209B) = K(RP209B) * C(SH) * C(SA2R5C2HY)
      W(RK210F) = K(RK210F) * C(SA2R5C2HX) * C(SH)
      W(RK210B) = K(RK210B) * C(SH2) * C(SA2R5C2HY)
      W(RP211F) = K(RP211F) * C(SA2R5C2HX) * C(SOH)
      W(RP211B) = K(RP211B) * C(SH2O) * C(SA2R5C2HY)
      W(RK212F) = K(RK212F) * C(SA2R5C2H2)
      W(RK212B) = K(RK212B) * C(SH) * C(SA2R5C2HX)
      W(RK213F) = K(RK213F) * C(SA2R5C2H2) * C(SH)
      W(RK213B) = K(RK213B) * C(SH2) * C(SA2R5C2HX)
      W(RP214F) = K(RP214F) * C(SA2R5C2H2) * C(SOH)
      W(RP214B) = K(RP214B) * C(SH2O) * C(SA2R5C2HX)
      W(RP301F) = K(RP301F) * C(SA1XC6H6) * C(SA1XXC6H5)
      W(RP301B) = K(RP301B) * C(SH) * C(SP2XC12H1)
      W(RP302F) = K(RP302F) * C(SA1XXC6H5) * C(SA1XXC6H5)
      W(RP302B) = K(RP302B) * C(SP2XC12H1)
      W(RP304F) = K(RP304F) * C(SP2XC12H1)
      W(RP304B) = K(RP304B) * C(SH) * C(SP2XXC12H)
      W(RP305F) = K(RP305F) * C(SP2XC12H1) * C(SH)
      W(RP305B) = K(RP305B) * C(SH2) * C(SP2XXC12H)
      W(RP306F) = K(RP306F) * C(SP2XC12H1) * C(SOH)
      W(RP306B) = K(RP306B) * C(SH2O) * C(SP2XXC12H)
      W(RK401F) = K(RK401F) * C(SA2C2HAYX) * C(SC2H2)
      W(RK401B) = K(RK401B) * C(SA3XXC14H)
      W(RK403F) = K(RK403F) * C(SA2C2HBYX) * C(SC2H2)
      W(RK403B) = K(RK403B) * C(SA3XXC14H)
      W(RP405F) = K(RP405F) * C(SA2C2H2AX) * C(SC2H2)
      W(RP405B) = K(RP405B) * C(SH) * C(SA3XC14H1)
      W(RP406F) = K(RP406F) * C(SA2C2H2BX) * C(SC2H2)
      W(RP406B) = K(RP406B) * C(SH) * C(SA3XC14H1)
      W(RP407F) = K(RP407F) * C(SP2XXC12H) * C(SC2H2)
      W(RP407B) = K(RP407B) * C(SH) * C(SA3XC14H1)
      W(RP408F) = K(RP408F) * C(SA2C2HAYX) * C(SC2H3)
      W(RP408B) = K(RP408B) * C(SA3XC14H1)
      W(RP409F) = K(RP409F) * C(SA2C2HBYX) * C(SC2H3)
      W(RP409B) = K(RP409B) * C(SA3XC14H1)
      W(RP410F) = K(RP410F) * C(SA2C2HAXC) * C(SC2H3)
      W(RP410B) = K(RP410B) * C(SH) * C(SA3XC14H1)
      W(RP411F) = K(RP411F) * C(SA2C2HBXC) * C(SC2H3)
      W(RP411B) = K(RP411B) * C(SH) * C(SA3XC14H1)
      W(RP412F) = K(RP412F) * C(SA2C2HAYX) * C(SC2H4)
      W(RP412B) = K(RP412B) * C(SH) * C(SA3XC14H1)
      W(RP413F) = K(RP413F) * C(SA2C2HBYX) * C(SC2H4)
      W(RP413B) = K(RP413B) * C(SH) * C(SA3XC14H1)
      W(RP414F) = K(RP414F) * C(SA2XXC10H) * C(SC4H4)
      W(RP414B) = K(RP414B) * C(SH) * C(SA3XC14H1)
      W(RP415F) = K(RP415F) * C(SA2YXC10H) * C(SC4H4)
      W(RP415B) = K(RP415B) * C(SH) * C(SA3XC14H1)
      W(RP416F) = K(RP416F) * C(SA1C2HXC8) * C(SA1XXC6H5)
      W(RP416B) = K(RP416B) * C(SH) * C(SA3XC14H1)
      W(RP417F) = K(RP417F) * C(SA1C2HYXC) * C(SA1XC6H6)
      W(RP417B) = K(RP417B) * C(SH) * C(SA3XC14H1)
      W(RP418F) = K(RP418F) * C(SA1C2HYXC) * C(SA1XXC6H5)
      W(RP418B) = K(RP418B) * C(SA3XC14H1)
      W(RP419F) = K(RP419F) * C(SA3XC14H1)
      W(RP419B) = K(RP419B) * C(SH) * C(SA3XXC14H)
      W(RK420F) = K(RK420F) * C(SA3XC14H1) * C(SH)
      W(RK420B) = K(RK420B) * C(SH2) * C(SA3XXC14H)
      W(RK421F) = K(RK421F) * C(SA3XC14H1) * C(SOH)
      W(RK421B) = K(RK421B) * C(SH2O) * C(SA3XXC14H)
      W(RP422F) = K(RP422F) * C(SA3XC14H1)
      W(RP422B) = K(RP422B) * C(SH) * C(SA3YXC14H)
      W(RK423F) = K(RK423F) * C(SA3XC14H1) * C(SH)
      W(RK423B) = K(RK423B) * C(SH2) * C(SA3YXC14H)
      W(RK424F) = K(RK424F) * C(SA3XC14H1) * C(SOH)
      W(RK424B) = K(RK424B) * C(SH2O) * C(SA3YXC14H)
      W(RP425F) = K(RP425F) * C(SA3XXC14H)
      W(RP425B) = K(RP425B) * C(SC2H2) * C(SA2R5XXC1)
      W(RP501F) = K(RP501F) * C(SA2R5C2HY) * C(SC2H2)
      W(RP501B) = K(RP501B) * C(SA3R5XXC1)
      W(RP502F) = K(RP502F) * C(SA2R5C2H2) * C(SC2H2)
      W(RP502B) = K(RP502B) * C(SH) * C(SA3R5XC16)
      W(RP503F) = K(RP503F) * C(SA3YXC14H) * C(SC2H2)
      W(RP503B) = K(RP503B) * C(SH) * C(SA3R5XC16)
      W(RP504F) = K(RP504F) * C(SA3YXC14H) * C(SC2H3)
      W(RP504B) = K(RP504B) * C(SH2) * C(SA3R5XC16)
      W(RP505F) = K(RP505F) * C(SA2R5XXC1) * C(SC4H4)
      W(RP505B) = K(RP505B) * C(SH) * C(SA3R5XC16)
      W(RP506F) = K(RP506F) * C(SA3R5XC16)
      W(RP506B) = K(RP506B) * C(SH) * C(SA3R5XXC1)
      W(RK507F) = K(RK507F) * C(SA3R5XC16) * C(SH)
      W(RK507B) = K(RK507B) * C(SH2) * C(SA3R5XXC1)
      W(RK508F) = K(RK508F) * C(SA3R5XC16) * C(SOH)
      W(RK508B) = K(RK508B) * C(SH2O) * C(SA3R5XXC1)
      W(RK600F) = K(RK600F) * C(SA3XXC14H) * C(SC2H2)
      W(RK600B) = K(RK600B) * C(SH) * C(SA4XC16H1)
      W(RP601F) = K(RP601F) * C(SA4XC16H1)
      W(RP601B) = K(RP601B) * C(SH) * C(SA4XXC16H)
      W(RK602F) = K(RK602F) * C(SA4XC16H1) * C(SH)
      W(RK602B) = K(RK602B) * C(SH2) * C(SA4XXC16H)
      W(RK603F) = K(RK603F) * C(SA4XC16H1) * C(SOH)
      W(RK603B) = K(RK603B) * C(SH2O) * C(SA4XXC16H)
      W(RK700F) = K(RK700F) * C(SA4XXC16H) * C(SC2H2)
      W(RK700B) = K(RK700B) * C(SH) * C(SA4R5XC18)
      W(RK701F) = K(RK701F) * C(SA3R5XXC1) * C(SC2H2)
      W(RK701B) = K(RK701B) * C(SH) * C(SA4R5XC18)
      W(RP800) = K(RP800) * C(SA2XC10H8) * C(SA1XXC6H5)
      W(RP801) = K(RP801) * C(SA2XXC10H) * C(SA1XC6H6)
      W(RP802) = K(RP802) * C(SA2XXC10H) * C(SA1XXC6H5)
      W(RCP01F) = K(RCP01F) * C(SC5H6)
      W(RCP01B) = K(RCP01B) * C(SH) * C(SC5H5)
      W(RCP02F) = K(RCP02F) * C(SC5H6) * C(SH)
      W(RCP02B) = K(RCP02B) * C(SH2) * C(SC5H5)
      W(RCP03F) = K(RCP03F) * C(SC5H6) * C(SH)
      W(RCP03B) = K(RCP03B) * C(SC2H2) * C(SAXC3H5)
      W(RCP04) = K(RCP04) * C(SC5H6) * C(SH)
      W(RCP05F) = K(RCP05F) * C(SC5H6) * C(SO)
      W(RCP05B) = K(RCP05B) * C(SOH) * C(SC5H5)
      W(RCP06F) = K(RCP06F) * C(SC5H6) * C(SOH)
      W(RCP06B) = K(RCP06B) * C(SH2O) * C(SC5H5)
      W(RCP07F) = K(RCP07F) * C(SC5H6) * C(SO2)
      W(RCP07B) = K(RCP07B) * C(SHO2) * C(SC5H5)
      W(RCP08F) = K(RCP08F) * C(SC5H6) * C(SHO2)
      W(RCP08B) = K(RCP08B) * C(SH2O2) * C(SC5H5)
      W(RCP09F) = K(RCP09F) * C(SC5H6) * C(SCH3)
      W(RCP09B) = K(RCP09B) * C(SCH4) * C(SC5H5)
      W(RCP10F) = K(RCP10F) * C(SC5H6) * C(SC2H3)
      W(RCP10B) = K(RCP10B) * C(SC2H4) * C(SC5H5)
      W(RCP11F) = K(RCP11F) * C(SC5H6) * C(SNXC4H5)
      W(RCP11B) = K(RCP11B) * C(SC4H6) * C(SC5H5)
      W(RCP12F) = K(RCP12F) * C(SC5H6) * C(SO)
      W(RCP12B) = K(RCP12B) * C(SH) * C(STXC5H5O)
      W(RCP13) = K(RCP13) * C(SC5H6) * C(SO)
      W(RCP14) = K(RCP14) * C(SC5H6) * C(SOH)
      W(RCP15) = K(RCP15) * C(SC5H6) * C(SOH)
      W(RCP16F) = K(RCP16F) * C(SC3H3) * C(SC2H2)
      W(RCP16B) = K(RCP16B) * C(SC5H5)
      W(RCP17) = K(RCP17) * C(SC5H5) * C(SC5H5)
      W(RCP18) = K(RCP18) * C(SC5H5) * C(SCH3)
      W(RCP19F) = K(RCP19F) * C(SC5H5) * C(SO)
      W(RCP19B) = K(RCP19B) * C(SH) * C(SC5H4O)
      W(RCP20F) = K(RCP20F) * C(SC5H5) * C(SO2)
      W(RCP20B) = K(RCP20B) * C(SOH) * C(SC5H4O)
      W(RCP21F) = K(RCP21F) * C(SC5H5) * C(SHO2)
      W(RCP21B) = K(RCP21B) * C(SOH) * C(SSXC5H5O)
      W(RCP22) = K(RCP22) * C(SC5H5) * C(SOH)
      W(RCP23F) = K(RCP23F) * C(SSXC5H5O)
      W(RCP23B) = K(RCP23B) * C(SH) * C(SC5H4O)
      W(RCP24) = K(RCP24) * C(STXC5H5O)
      W(RCP25) = K(RCP25) * C(SSXC5H5O) * C(SH)
      W(RCP26) = K(RCP26) * C(SSXC5H5O) * C(SO)
      W(RCP27) = K(RCP27) * C(SSXC5H5O) * C(SOH)
      W(RCP28) = K(RCP28) * C(SSXC5H5O) * C(SO2)
      W(RCP29) = K(RCP29) * C(STXC5H5O) * C(SH)
      W(RCP30) = K(RCP30) * C(STXC5H5O) * C(SO2)
      W(RCP31) = K(RCP31) * C(SC5H4O)
      W(RCP32F) = K(RCP32F) * C(SC5H4O) * C(SH)
      W(RCP32B) = K(RCP32B) * C(STXC5H5O)
      W(RCP33F) = K(RCP33F) * C(SC5H4O) * C(SO)
      W(RCP33B) = K(RCP33B) * C(SCO2) * C(SC4H4)
      W(RCP34) = K(RCP34) * C(SC5H5) * C(SC5H6)
      W(RI00F) = K(RI00F) * C(SA1XXC6H5) * C(SC3H3)
      W(RI00B) = K(RI00B) * C(SC9H8)
      W(RI01F) = K(RI01F) * C(SC9H8)
      W(RI01B) = K(RI01B) * C(SH) * C(SC9H7)
      W(RI02F) = K(RI02F) * C(SC9H8) * C(SH)
      W(RI02B) = K(RI02B) * C(SH2) * C(SC9H7)
      W(RI03F) = K(RI03F) * C(SA1CH2XC7) * C(SC2H2)
      W(RI03B) = K(RI03B) * C(SH) * C(SC9H8)
      W(RI05F) = K(RI05F) * C(SC9H8) * C(SO)
      W(RI05B) = K(RI05B) * C(SOH) * C(SC9H7)
      W(RI06F) = K(RI06F) * C(SC9H8) * C(SOH)
      W(RI06B) = K(RI06B) * C(SH2O) * C(SC9H7)
      W(RI07F) = K(RI07F) * C(SC9H8) * C(SO2)
      W(RI07B) = K(RI07B) * C(SHO2) * C(SC9H7)
      W(RI08F) = K(RI08F) * C(SC9H8) * C(SHO2)
      W(RI08B) = K(RI08B) * C(SH2O2) * C(SC9H7)
      W(RI09F) = K(RI09F) * C(SC9H8) * C(SCH3)
      W(RI09B) = K(RI09B) * C(SCH4) * C(SC9H7)
      W(RI12) = K(RI12) * C(SC9H8) * C(SO)
      W(RI15) = K(RI15) * C(SC9H8) * C(SOH)
      W(RI17) = K(RI17) * C(SC9H7) * C(SC5H5)
      W(RI18) = K(RI18) * C(SC9H7) * C(SCH3)
      W(RI19F) = K(RI19F) * C(SC9H7) * C(SO)
      W(RI19B) = K(RI19B) * C(SH) * C(SC9H6O)
      W(RI20F) = K(RI20F) * C(SC9H7) * C(SO2)
      W(RI20B) = K(RI20B) * C(SOH) * C(SC9H6O)
      W(RI21) = K(RI21) * C(SC9H7) * C(SHO2)
      W(RI22) = K(RI22) * C(SC9H7) * C(SOH)
      W(RI23) = K(RI23) * C(SC9H7) * C(SC3H3)
      W(RI25) = K(RI25) * C(SC9H7)
      W(RI26) = K(RI26) * C(SC9H7)
      W(RI31) = K(RI31) * C(SC9H6O)
      W(RI32) = K(RI32) * C(SC9H6O) * C(SH)
      W(RT01F) = K(RT01F) * C(SA1CH3XC7) * C(SH)
      W(RT01B) = K(RT01B) * C(SCH3) * C(SA1XC6H6)
      W(RT02F) = K(RT02F) * C(SA1CH3XC7)
      W(RT02B) = K(RT02B) * C(SH) * C(SA1CH2XC7)
      W(RT03F) = K(RT03F) * C(SA1CH3XC7)
      W(RT03B) = K(RT03B) * C(SCH3) * C(SA1XXC6H5)
      W(RT04F) = K(RT04F) * C(SA1CH2XC7) * C(SH)
      W(RT04B) = K(RT04B) * C(SCH3) * C(SA1XXC6H5)
      W(RT05F) = K(RT05F) * C(SA1CH2XC7)
      W(RT05B) = K(RT05B) * C(SC2H2) * C(SC5H5)
      W(RT06F) = K(RT06F) * C(SA1CH3XC7) * C(SO2)
      W(RT06B) = K(RT06B) * C(SHO2) * C(SA1CH2XC7)
      W(RT07F) = K(RT07F) * C(SA1CH3XC7) * C(SH)
      W(RT07B) = K(RT07B) * C(SH2) * C(SA1CH2XC7)
      W(RT08F) = K(RT08F) * C(SA1CH3XC7) * C(SOH)
      W(RT08B) = K(RT08B) * C(SH2O) * C(SA1CH2XC7)
      W(RT09F) = K(RT09F) * C(SA1CH3XC7) * C(SOH)
      W(RT09B) = K(RT09B) * C(SCH3) * C(SA1OHXC6H)
      W(RT10F) = K(RT10F) * C(SA1CH3XC7) * C(SOH)
      W(RT10B) = K(RT10B) * C(SH) * C(SHOA1CH3X)
      W(RT11F) = K(RT11F) * C(SA1CH3XC7) * C(SO)
      W(RT11B) = K(RT11B) * C(SOH) * C(SA1CH2XC7)
      W(RT12F) = K(RT12F) * C(SA1CH3XC7) * C(SO)
      W(RT12B) = K(RT12B) * C(SHOA1CH3X)
      W(RT13F) = K(RT13F) * C(SA1CH3XC7) * C(SO)
      W(RT13B) = K(RT13B) * C(SH) * C(SOA1CH3XC)
      W(RT14F) = K(RT14F) * C(SA1CH3XC7) * C(SCH3)
      W(RT14B) = K(RT14B) * C(SCH4) * C(SA1CH2XC7)
      W(RT15F) = K(RT15F) * C(SA1CH3XC7) * C(SHO2)
      W(RT15B) = K(RT15B) * C(SH2O2) * C(SA1CH2XC7)
      W(RT16F) = K(RT16F) * C(SA1CH3XC7) * C(SA1XXC6H5)
      W(RT16B) = K(RT16B) * C(SA1XC6H6) * C(SA1CH2XC7)
      W(RT17F) = K(RT17F) * C(SA1CH2XC7) * C(SO)
      W(RT17B) = K(RT17B) * C(SA1CH2OXC)
      W(RT18F) = K(RT18F) * C(SA1CH2XC7) * C(SOH)
      W(RT18B) = K(RT18B) * C(SA1CH2OHX)
      W(RT19F) = K(RT19F) * C(SA1CH2XC7) * C(SHO2)
      W(RT19B) = K(RT19B) * C(SOH) * C(SA1CH2OXC)
      W(RT20) = K(RT20) * C(SA1CH2XC7) * C(SC3H3)
      W(RT21F) = K(RT21F) * C(SA1CH2XC7) * C(SO2)
      W(RT21B) = K(RT21B) * C(SOH) * C(SA1CHOXC7)
      W(RT22F) = K(RT22F) * C(SA1CH2XC7) * C(SO2)
      W(RT22B) = K(RT22B) * C(SCH2O) * C(SA1OXC6H5)
      W(RT23F) = K(RT23F) * C(SA1CH2OXC) * C(SH)
      W(RT23B) = K(RT23B) * C(SA1CH2OHX)
      W(RT24F) = K(RT24F) * C(SA1CH2OHX) * C(SH)
      W(RT24B) = K(RT24B) * C(SH2) * C(SA1CH2OXC)
      W(RT25F) = K(RT25F) * C(SA1CH2OHX) * C(SO)
      W(RT25B) = K(RT25B) * C(SOH) * C(SA1CH2OXC)
      W(RT26F) = K(RT26F) * C(SA1CH2OHX) * C(SOH)
      W(RT26B) = K(RT26B) * C(SH2O) * C(SA1CH2OXC)
      W(RT27F) = K(RT27F) * C(SA1CH2OHX) * C(SCH3)
      W(RT27B) = K(RT27B) * C(SCH4) * C(SA1CH2OXC)
      W(RT28F) = K(RT28F) * C(SA1XXC6H5) * C(SCH2O)
      W(RT28B) = K(RT28B) * C(SHCO) * C(SA1XC6H6)
      W(RT29F) = K(RT29F) * C(SA1CH2OXC)
      W(RT29B) = K(RT29B) * C(SH) * C(SA1CHOXC7)
      W(RT30F) = K(RT30F) * C(SA1CH2OXC)
      W(RT30B) = K(RT30B) * C(SCH2O) * C(SA1XXC6H5)
      W(RT31F) = K(RT31F) * C(SA1CH2OXC)
      W(RT31B) = K(RT31B) * C(SHCO) * C(SA1XC6H6)
      W(RT32F) = K(RT32F) * C(SA1CH2OXC) * C(SH)
      W(RT32B) = K(RT32B) * C(SH2) * C(SA1CHOXC7)
      W(RT33F) = K(RT33F) * C(SA1CH2OXC) * C(SO)
      W(RT33B) = K(RT33B) * C(SOH) * C(SA1CHOXC7)
      W(RT34F) = K(RT34F) * C(SA1CH2OXC) * C(SOH)
      W(RT34B) = K(RT34B) * C(SH2O) * C(SA1CHOXC7)
      W(RT35F) = K(RT35F) * C(SA1CH2OXC) * C(SO2)
      W(RT35B) = K(RT35B) * C(SHO2) * C(SA1CHOXC7)
      W(RT36) = K(RT36) * C(SA1CHOXC7)
      W(RT37) = K(RT37) * C(SA1CHOXC7) * C(SH)
      W(RT38) = K(RT38) * C(SA1CHOXC7) * C(SO)
      W(RT39) = K(RT39) * C(SA1CHOXC7) * C(SOH)
      W(RT40) = K(RT40) * C(SA1CHOXC7) * C(SO2)
      W(RT41) = K(RT41) * C(SA1CHOXC7) * C(SHO2)
      W(RT42) = K(RT42) * C(SA1CHOXC7) * C(SCH3)
      W(RT43F) = K(RT43F) * C(SHOA1CH3X)
      W(RT43B) = K(RT43B) * C(SH) * C(SOA1CH3XC)
      W(RT44F) = K(RT44F) * C(SHOA1CH3X) * C(SH)
      W(RT44B) = K(RT44B) * C(SH2) * C(SOA1CH3XC)
      W(RT45F) = K(RT45F) * C(SHOA1CH3X) * C(SO)
      W(RT45B) = K(RT45B) * C(SOH) * C(SOA1CH3XC)
      W(RT46F) = K(RT46F) * C(SHOA1CH3X) * C(SOH)
      W(RT46B) = K(RT46B) * C(SH2O) * C(SOA1CH3XC)
      W(RT47) = K(RT47) * C(SOA1CH3XC)
      W(RT50F) = K(RT50F) * C(SA1CH3XC7)
      W(RT50B) = K(RT50B) * C(SH) * C(SA1CH3YXC)
      W(RT51F) = K(RT51F) * C(SA1CH3XC7) * C(SH)
      W(RT51B) = K(RT51B) * C(SH2) * C(SA1CH3YXC)
      W(RT52F) = K(RT52F) * C(SA1CH3XC7) * C(SO)
      W(RT52B) = K(RT52B) * C(SOH) * C(SA1CH3YXC)
      W(RT53F) = K(RT53F) * C(SA1CH3XC7) * C(SOH)
      W(RT53B) = K(RT53B) * C(SH2O) * C(SA1CH3YXC)
      W(RT54F) = K(RT54F) * C(SA1CH3XC7) * C(SCH3)
      W(RT54B) = K(RT54B) * C(SCH4) * C(SA1CH3YXC)
      W(RT55F) = K(RT55F) * C(SA1CH3YXC) * C(SO)
      W(RT55B) = K(RT55B) * C(SOA1CH3XC)
      W(RT56F) = K(RT56F) * C(SA1CH3YXC) * C(SOH)
      W(RT56B) = K(RT56B) * C(SH) * C(SOA1CH3XC)
      W(RT57F) = K(RT57F) * C(SA1CH3YXC) * C(SHO2)
      W(RT57B) = K(RT57B) * C(SOH) * C(SOA1CH3XC)
      W(RT58F) = K(RT58F) * C(SA1CH3YXC) * C(SO2)
      W(RT58B) = K(RT58B) * C(SO) * C(SOA1CH3XC)
      W(RT59) = K(RT59) * C(SA1CH3YXC) * C(SO2)
      W(RT60) = K(RT60) * C(SA1CH3YXC) * C(SO2)
      W(RE01F) = K(RE01F) * C(SA1C2H4XC) * C(SH)
      W(RE01B) = K(RE01B) * C(SA1C2H5XC)
      W(RE02F) = K(RE02F) * C(SA1CH2XC7) * C(SCH3)
      W(RE02B) = K(RE02B) * C(SA1C2H5XC)
      W(RE03F) = K(RE03F) * C(SA1XXC6H5) * C(SC2H5)
      W(RE03B) = K(RE03B) * C(SA1C2H5XC)
      W(RE04F) = K(RE04F) * C(SA1C2H5XC) * C(SH)
      W(RE04B) = K(RE04B) * C(SC2H5) * C(SA1XC6H6)
      W(RE05F) = K(RE05F) * C(SA1C2H5XC) * C(SOH)
      W(RE05B) = K(RE05B) * C(SC2H5) * C(SA1OHXC6H)
      W(RE06F) = K(RE06F) * C(SA1C2H5XC) * C(SH)
      W(RE06B) = K(RE06B) * C(SH2) * C(SA1C2H4XC)
      W(RE07F) = K(RE07F) * C(SA1C2H5XC) * C(SO)
      W(RE07B) = K(RE07B) * C(SOH) * C(SA1C2H4XC)
      W(RE08F) = K(RE08F) * C(SA1C2H5XC) * C(SOH)
      W(RE08B) = K(RE08B) * C(SH2O) * C(SA1C2H4XC)
      W(RE09F) = K(RE09F) * C(SA1C2H5XC) * C(SHO2)
      W(RE09B) = K(RE09B) * C(SH2O2) * C(SA1C2H4XC)
      W(RE10F) = K(RE10F) * C(SA1C2H5XC) * C(SCH3)
      W(RE10B) = K(RE10B) * C(SCH4) * C(SA1C2H4XC)
      W(RE11F) = K(RE11F) * C(SA1C2H4XC)
      W(RE11B) = K(RE11B) * C(SC2H4) * C(SA1XXC6H5)
      W(RE12F) = K(RE12F) * C(SA1C2H4XC)
      W(RE12B) = K(RE12B) * C(SH) * C(SA1C2H3XC)
      W(RE13F) = K(RE13F) * C(SA1C2H4XC) * C(SH)
      W(RE13B) = K(RE13B) * C(SH2) * C(SA1C2H3XC)
      W(RE14F) = K(RE14F) * C(SA1C2H4XC) * C(SOH)
      W(RE14B) = K(RE14B) * C(SH2O) * C(SA1C2H3XC)
      W(RE15F) = K(RE15F) * C(SA1C2H4XC) * C(SCH3)
      W(RE15B) = K(RE15B) * C(SCH4) * C(SA1C2H3XC)
      W(RE16F) = K(RE16F) * C(SA1C2H4XC) * C(SO2)
      W(RE16B) = K(RE16B) * C(SHO2) * C(SA1C2H3XC)
      W(RE17F) = K(RE17F) * C(SA1C2H4XC) * C(SO)
      W(RE17B) = K(RE17B) * C(SCH2O) * C(SA1CH2XC7)
      W(RE18F) = K(RE18F) * C(SA1C2H4XC) * C(SO)
      W(RE18B) = K(RE18B) * C(SCH3) * C(SA1CHOXC7)
      W(RE19) = K(RE19) * C(SA1C2H4XC) * C(SHO2)
      W(RE30) = K(RE30) * C(SA1C2H4XC) * C(SO2)
      W(RE31) = K(RE31) * C(SA1C2H4XC) * C(SO2)
      W(RE32) = K(RE32) * C(SC8H9O2)
      W(RE33) = K(RE33) * C(SC8H9O2)
      W(RE34) = K(RE34) * C(SC8H9O2)
      W(RE35) = K(RE35) * C(SC8H9O2)
      W(RE36) = K(RE36) * C(SC8H8OOH) * C(SO2)
      W(RE37) = K(RE37) * C(SOC8H7OOH)
      W(RST01) = K(RST01) * C(SA1C2H3XC)
      W(RST02F) = K(RST02F) * C(SA1C2H3XC) * C(SO)
      W(RST02B) = K(RST02B) * C(SHCO) * C(SA1CH2XC7)
      W(RST03F) = K(RST03F) * C(SA1C2H3XC) * C(SCH3)
      W(RST03B) = K(RST03B) * C(SCH4) * C(SA1C2H2XC)
      W(RST04F) = K(RST04F) * C(SA1C2H3XC) * C(SOH)
      W(RST04B) = K(RST04B) * C(SCH3) * C(SA1CHOXC7)
      W(RST05F) = K(RST05F) * C(SA1C2H3XC) * C(SOH)
      W(RST05B) = K(RST05B) * C(SCH2O) * C(SA1CH2XC7)
      W(RST06F) = K(RST06F) * C(SA1C2H3XC) * C(SOH)
      W(RST06B) = K(RST06B) * C(SA1OHXC6H) * C(SC2H3)
      W(RST10F) = K(RST10F) * C(SA1C2H3XC) * C(SO)
      W(RST10B) = K(RST10B) * C(SOH) * C(SA1C2H3YX)
      W(RST11F) = K(RST11F) * C(SA1C2H2XC) * C(SO)
      W(RST11B) = K(RST11B) * C(SCO) * C(SA1CH2XC7)
      W(RST12F) = K(RST12F) * C(SA1C2H2XC) * C(SO2)
      W(RST12B) = K(RST12B) * C(SHO2) * C(SA1C2HXC8)
      W(RST13) = K(RST13) * C(SA1C2H2XC) * C(SO2)
      W(RST14F) = K(RST14F) * C(SA1C2H2XC) * C(SO2)
      W(RST14B) = K(RST14B) * C(SHCO) * C(SA1CHOXC7)
      W(RST00F) = K(RST00F) * C(SA1C2HXC8) * C(SOH)
      W(RST00B) = K(RST00B) * C(SHCCO) * C(SA1XC6H6)
      W(RXY00F) = K(RXY00F) * C(SA1CH3CH3)
      W(RXY00B) = K(RXY00B) * C(SH) * C(SA1CH3CH2)
      W(RXY01F) = K(RXY01F) * C(SA1CH3CH3)
      W(RXY01B) = K(RXY01B) * C(SCH3) * C(SA1CH3YXC)
      W(RXY02F) = K(RXY02F) * C(SA1CH3CH3) * C(SH)
      W(RXY02B) = K(RXY02B) * C(SH2) * C(SA1CH3CH2)
      W(RXY03F) = K(RXY03F) * C(SA1CH3CH3) * C(SO)
      W(RXY03B) = K(RXY03B) * C(SOH) * C(SA1CH3CH2)
      W(RXY04F) = K(RXY04F) * C(SA1CH3CH3) * C(SOH)
      W(RXY04B) = K(RXY04B) * C(SH2O) * C(SA1CH3CH2)
      W(RXY05F) = K(RXY05F) * C(SA1CH3CH3) * C(SO2)
      W(RXY05B) = K(RXY05B) * C(SHO2) * C(SA1CH3CH2)
      W(RXY06F) = K(RXY06F) * C(SA1CH3CH3) * C(SHO2)
      W(RXY06B) = K(RXY06B) * C(SH2O2) * C(SA1CH3CH2)
      W(RXY07F) = K(RXY07F) * C(SA1CH3CH3) * C(SCH3)
      W(RXY07B) = K(RXY07B) * C(SCH4) * C(SA1CH3CH2)
      W(RXY09F) = K(RXY09F) * C(SA1CH3CH3) * C(SH)
      W(RXY09B) = K(RXY09B) * C(SCH3) * C(SA1CH3XC7)
      W(RXY10F) = K(RXY10F) * C(SA1CH3CH3) * C(SOH)
      W(RXY10B) = K(RXY10B) * C(SCH3) * C(SHOA1CH3X)
      W(RXY11) = K(RXY11) * C(SA1CH3CH3) * C(SO)
      W(RXY12) = K(RXY12) * C(SA1CH3CH2)
      W(RXY13F) = K(RXY13F) * C(SA1CH3CH2) * C(SH)
      W(RXY13B) = K(RXY13B) * C(SCH3) * C(SA1CH3YXC)
      W(RXY14F) = K(RXY14F) * C(SA1CH3CH2) * C(SO)
      W(RXY14B) = K(RXY14B) * C(SH) * C(SA1CH3CHO)
      W(RXY15F) = K(RXY15F) * C(SA1CH3CH2) * C(SO)
      W(RXY15B) = K(RXY15B) * C(SCH2O) * C(SA1CH3YXC)
      W(RXY16F) = K(RXY16F) * C(SA1CH3CH2) * C(SO)
      W(RXY16B) = K(RXY16B) * C(SHCO) * C(SA1CH3XC7)
      W(RXY17) = K(RXY17) * C(SA1CH3CH2) * C(SOH)
      W(RXY18F) = K(RXY18F) * C(SA1CH3CH2) * C(SO2)
      W(RXY18B) = K(RXY18B) * C(SOH) * C(SA1CH3CHO)
      W(RXY19F) = K(RXY19F) * C(SA1CH3CH2) * C(SO2)
      W(RXY19B) = K(RXY19B) * C(SCH2O) * C(SOA1CH3XC)
      W(RXY201) = K(RXY201) * C(SA1CH3CH2) * C(SHO2)
      W(RXY202) = K(RXY202) * C(SA1CH3CH2) * C(SHO2)
      W(RXY203) = K(RXY203) * C(SA1CH3CH2) * C(SHO2)
      W(RXY22) = K(RXY22) * C(SA1CH3CH2) * C(SC3H3)
      W(RXY23F) = K(RXY23F) * C(SA1CH3CHO)
      W(RXY23B) = K(RXY23B) * C(SH) * C(SA1CHOCH2)
      W(RXY24) = K(RXY24) * C(SA1CH3CHO)
      W(RXY25) = K(RXY25) * C(SA1CH3CHO)
      W(RXY26F) = K(RXY26F) * C(SA1CH3CHO) * C(SH)
      W(RXY26B) = K(RXY26B) * C(SH2) * C(SA1CHOCH2)
      W(RXY27F) = K(RXY27F) * C(SA1CH3CHO) * C(SO)
      W(RXY27B) = K(RXY27B) * C(SOH) * C(SA1CHOCH2)
      W(RXY28F) = K(RXY28F) * C(SA1CH3CHO) * C(SOH)
      W(RXY28B) = K(RXY28B) * C(SH2O) * C(SA1CHOCH2)
      W(RXY29F) = K(RXY29F) * C(SA1CH3CHO) * C(SO2)
      W(RXY29B) = K(RXY29B) * C(SHO2) * C(SA1CHOCH2)
      W(RXY30F) = K(RXY30F) * C(SA1CH3CHO) * C(SHO2)
      W(RXY30B) = K(RXY30B) * C(SH2O2) * C(SA1CHOCH2)
      W(RXY31F) = K(RXY31F) * C(SA1CH3CHO) * C(SCH3)
      W(RXY31B) = K(RXY31B) * C(SCH4) * C(SA1CHOCH2)
      W(RXY33) = K(RXY33) * C(SA1CH3CHO) * C(SH)
      W(RXY34) = K(RXY34) * C(SA1CH3CHO) * C(SO)
      W(RXY35) = K(RXY35) * C(SA1CH3CHO) * C(SOH)
      W(RXY36) = K(RXY36) * C(SA1CH3CHO) * C(SO2)
      W(RXY37) = K(RXY37) * C(SA1CH3CHO) * C(SHO2)
      W(RXY38) = K(RXY38) * C(SA1CH3CHO) * C(SCH3)
      W(RXY39F) = K(RXY39F) * C(SA1CH3CHO) * C(SH)
      W(RXY39B) = K(RXY39B) * C(SHCO) * C(SA1CH3XC7)
      W(RXY40F) = K(RXY40F) * C(SA1CH3CHO) * C(SH)
      W(RXY40B) = K(RXY40B) * C(SCH3) * C(SA1CHOXC7)
      W(RXY41F) = K(RXY41F) * C(SA1CH3CHO) * C(SOH)
      W(RXY41B) = K(RXY41B) * C(SHCO) * C(SHOA1CH3X)
      W(RXY42) = K(RXY42) * C(SA1CH3CHO) * C(SOH)
      W(RXY43F) = K(RXY43F) * C(SA1CHOCH2) * C(SO)
      W(RXY43B) = K(RXY43B) * C(SH) * C(SA1CHOCHO)
      W(RXY44) = K(RXY44) * C(SA1CHOCH2) * C(SO)
      W(RXY45F) = K(RXY45F) * C(SA1CHOCH2) * C(SO)
      W(RXY45B) = K(RXY45B) * C(SHCO) * C(SA1CHOXC7)
      W(RXY46) = K(RXY46) * C(SA1CHOCH2) * C(SOH)
      W(RXY47F) = K(RXY47F) * C(SA1CHOCH2) * C(SO2)
      W(RXY47B) = K(RXY47B) * C(SOH) * C(SA1CHOCHO)
      W(RXY48) = K(RXY48) * C(SA1CHOCH2) * C(SHO2)
      W(RXY50) = K(RXY50) * C(SA1CHOCHO)
      W(RXY51) = K(RXY51) * C(SA1CHOCHO) * C(SH)
      W(RXY52) = K(RXY52) * C(SA1CHOCHO) * C(SO)
      W(RXY53) = K(RXY53) * C(SA1CHOCHO) * C(SOH)
      W(RXY54) = K(RXY54) * C(SA1CHOCHO) * C(SO2)
      W(RXY55) = K(RXY55) * C(SA1CHOCHO) * C(SHO2)
      W(RXY56) = K(RXY56) * C(SA1CHOCHO) * C(SCH3)
      W(RXY57F) = K(RXY57F) * C(SA1CHOCHO) * C(SH)
      W(RXY57B) = K(RXY57B) * C(SHCO) * C(SA1CHOXC7)
      W(RXY58) = K(RXY58) * C(SA1CHOCHO) * C(SOH)
      W(RN01F) = K(RN01F) * C(SA2CH3XC1) * C(SH)
      W(RN01B) = K(RN01B) * C(SCH3) * C(SA2XC10H8)
      W(RN02F) = K(RN02F) * C(SA2CH3XC1) * C(SOH)
      W(RN02B) = K(RN02B) * C(SCH3) * C(SA2OHXC10)
      W(RN04F) = K(RN04F) * C(SA2CH3XC1)
      W(RN04B) = K(RN04B) * C(SH) * C(SA2CH2XC1)
      W(RN05F) = K(RN05F) * C(SA2CH3XC1)
      W(RN05B) = K(RN05B) * C(SCH3) * C(SA2XXC10H)
      W(RN06F) = K(RN06F) * C(SA2CH2XC1) * C(SH)
      W(RN06B) = K(RN06B) * C(SCH3) * C(SA2XXC10H)
      W(RN07F) = K(RN07F) * C(SA2CH2XC1)
      W(RN07B) = K(RN07B) * C(SC2H2) * C(SC9H7)
      W(RN08F) = K(RN08F) * C(SA2CH3XC1) * C(SH)
      W(RN08B) = K(RN08B) * C(SH2) * C(SA2CH2XC1)
      W(RN09F) = K(RN09F) * C(SA2CH3XC1) * C(SO)
      W(RN09B) = K(RN09B) * C(SOH) * C(SA2CH2XC1)
      W(RN10F) = K(RN10F) * C(SA2CH3XC1) * C(SOH)
      W(RN10B) = K(RN10B) * C(SH2O) * C(SA2CH2XC1)
      W(RN11F) = K(RN11F) * C(SA2CH3XC1) * C(SO2)
      W(RN11B) = K(RN11B) * C(SHO2) * C(SA2CH2XC1)
      W(RN12F) = K(RN12F) * C(SA2CH3XC1) * C(SCH3)
      W(RN12B) = K(RN12B) * C(SCH4) * C(SA2CH2XC1)
      W(RN13F) = K(RN13F) * C(SA2CH3XC1) * C(SHO2)
      W(RN13B) = K(RN13B) * C(SH2O2) * C(SA2CH2XC1)
      W(RN14) = K(RN14) * C(SA2CH3XC1) * C(SO)
      W(RN15) = K(RN15) * C(SA2CH3XC1) * C(SO)
      W(RN16F) = K(RN16F) * C(SA2CH2XC1) * C(SO)
      W(RN16B) = K(RN16B) * C(SA2CH2OXC)
      W(RN17) = K(RN17) * C(SA2CH2XC1) * C(SOH)
      W(RN18F) = K(RN18F) * C(SA2CH2XC1) * C(SHO2)
      W(RN18B) = K(RN18B) * C(SOH) * C(SA2CH2OXC)
      W(RN18) = K(RN18) * C(SA2CH2XC1) * C(SC3H3)
      W(RN20F) = K(RN20F) * C(SA2CH2XC1) * C(SO2)
      W(RN20B) = K(RN20B) * C(SOH) * C(SA2CHOXC1)
      W(RN21F) = K(RN21F) * C(SA2CH2XC1) * C(SO2)
      W(RN21B) = K(RN21B) * C(SCH2O) * C(SA2OXC10H)
      W(RN22F) = K(RN22F) * C(SA2CH2OXC)
      W(RN22B) = K(RN22B) * C(SH) * C(SA2CHOXC1)
      W(RN23F) = K(RN23F) * C(SA2CH2OXC)
      W(RN23B) = K(RN23B) * C(SCH2O) * C(SA2XXC10H)
      W(RN24F) = K(RN24F) * C(SA2CH2OXC) * C(SH)
      W(RN24B) = K(RN24B) * C(SH2) * C(SA2CHOXC1)
      W(RN25F) = K(RN25F) * C(SA2CH2OXC) * C(SO)
      W(RN25B) = K(RN25B) * C(SOH) * C(SA2CHOXC1)
      W(RN26F) = K(RN26F) * C(SA2CH2OXC) * C(SOH)
      W(RN26B) = K(RN26B) * C(SH2O) * C(SA2CHOXC1)
      W(RN27F) = K(RN27F) * C(SA2CH2OXC) * C(SO2)
      W(RN27B) = K(RN27B) * C(SHO2) * C(SA2CHOXC1)
      W(RN28) = K(RN28) * C(SA2CHOXC1)
      W(RN29) = K(RN29) * C(SA2CHOXC1) * C(SH)
      W(RN30) = K(RN30) * C(SA2CHOXC1) * C(SO)
      W(RN31) = K(RN31) * C(SA2CHOXC1) * C(SOH)
      W(RN32) = K(RN32) * C(SA2CHOXC1) * C(SO2)
      W(RN33) = K(RN33) * C(SA2CHOXC1) * C(SHO2)
      W(RN34) = K(RN34) * C(SA2CHOXC1) * C(SCH3)
      W(ROX00F) = K(ROX00F) * C(SA1XC6H6)
      W(ROX00B) = K(ROX00B) * C(SH) * C(SA1XXC6H5)
      W(ROX01F) = K(ROX01F) * C(SA1XXC6H5)
      W(ROX01B) = K(ROX01B) * C(SH) * C(SOXC6H4)
      W(ROX02F) = K(ROX02F) * C(SC4H2) * C(SC2H2)
      W(ROX02B) = K(ROX02B) * C(SOXC6H4)
      W(ROX03F) = K(ROX03F) * C(SA1XC6H6) * C(SH)
      W(ROX03B) = K(ROX03B) * C(SH2) * C(SA1XXC6H5)
      W(ROX04F) = K(ROX04F) * C(SA1XC6H6) * C(SOH)
      W(ROX04B) = K(ROX04B) * C(SH2O) * C(SA1XXC6H5)
      W(ROX05F) = K(ROX05F) * C(SA1XC6H6) * C(SOH)
      W(ROX05B) = K(ROX05B) * C(SH) * C(SA1OHXC6H)
      W(ROX06F) = K(ROX06F) * C(SA1XC6H6) * C(SO2)
      W(ROX06B) = K(ROX06B) * C(SHO2) * C(SA1XXC6H5)
      W(ROX99F) = K(ROX99F) * C(SA1XC6H6) * C(SO)
      W(ROX99B) = K(ROX99B) * C(SH) * C(SA1OXC6H5)
      W(ROX07F) = K(ROX07F) * C(SA1XC6H6) * C(SO)
      W(ROX07B) = K(ROX07B) * C(SH) * C(SA1OXC6H5)
      W(ROX08F) = K(ROX08F) * C(SA1XC6H6) * C(SO)
      W(ROX08B) = K(ROX08B) * C(SA1OHXC6H)
      W(ROX09F) = K(ROX09F) * C(SA1XC6H6) * C(SO)
      W(ROX09B) = K(ROX09B) * C(SCO) * C(SC5H6)
      W(ROX10F) = K(ROX10F) * C(SA1XC6H6) * C(SO)
      W(ROX10B) = K(ROX10B) * C(SOH) * C(SA1XXC6H5)
      W(ROX11F) = K(ROX11F) * C(SA1XXC6H5) * C(SO2)
      W(ROX11B) = K(ROX11B) * C(SO) * C(SA1OXC6H5)
      W(ROX12F) = K(ROX12F) * C(SA1XXC6H5) * C(SO2)
      W(ROX12B) = K(ROX12B) * C(SH) * C(SOC6H4O)
      W(ROX13F) = K(ROX13F) * C(SA1XXC6H5) * C(SO)
      W(ROX13B) = K(ROX13B) * C(SA1OXC6H5)
      W(ROX14F) = K(ROX14F) * C(SA1XXC6H5) * C(SOH)
      W(ROX14B) = K(ROX14B) * C(SH) * C(SA1OXC6H5)
      W(ROX15F) = K(ROX15F) * C(SA1XXC6H5) * C(SHO2)
      W(ROX15B) = K(ROX15B) * C(SOH) * C(SA1OXC6H5)
      W(ROX16F) = K(ROX16F) * C(SA1XXC6H5) * C(SCH4)
      W(ROX16B) = K(ROX16B) * C(SCH3) * C(SA1XC6H6)
      W(ROX17F) = K(ROX17F) * C(SA1OHXC6H)
      W(ROX17B) = K(ROX17B) * C(SCO) * C(SC5H6)
      W(ROX18F) = K(ROX18F) * C(SA1OHXC6H)
      W(ROX18B) = K(ROX18B) * C(SH) * C(SA1OXC6H5)
      W(ROX19F) = K(ROX19F) * C(SA1OHXC6H) * C(SH)
      W(ROX19B) = K(ROX19B) * C(SH2) * C(SA1OXC6H5)
      W(ROX20F) = K(ROX20F) * C(SA1OHXC6H) * C(SOH)
      W(ROX20B) = K(ROX20B) * C(SH2O) * C(SA1OXC6H5)
      W(ROX21F) = K(ROX21F) * C(SA1OHXC6H) * C(SCH3)
      W(ROX21B) = K(ROX21B) * C(SCH4) * C(SA1OXC6H5)
      W(ROX22F) = K(ROX22F) * C(SA1OHXC6H) * C(SO2)
      W(ROX22B) = K(ROX22B) * C(SHO2) * C(SA1OXC6H5)
      W(ROX23F) = K(ROX23F) * C(SA1OXC6H5)
      W(ROX23B) = K(ROX23B) * C(SC5H5) * C(SCO)
      W(ROX24F) = K(ROX24F) * C(SA1OXC6H5) * C(SO)
      W(ROX24B) = K(ROX24B) * C(SH) * C(SOC6H4O)
      W(ROX25F) = K(ROX25F) * C(SA1OXC6H5) * C(SO2)
      W(ROX25B) = K(ROX25B) * C(SOH) * C(SOC6H4O)
      W(ROX26F) = K(ROX26F) * C(SOC6H4O)
      W(ROX26B) = K(ROX26B) * C(SCO) * C(SC5H4O)
      W(ROX27F) = K(ROX27F) * C(SOC6H4O) * C(SH)
      W(ROX27B) = K(ROX27B) * C(SCO) * C(STXC5H5O)
      W(ROX28) = K(ROX28) * C(SOC6H4O) * C(SO)
      W(ROX30F) = K(ROX30F) * C(SA2XC10H8) * C(SO)
      W(ROX30B) = K(ROX30B) * C(SH) * C(SA2OXC10H)
      W(ROX31F) = K(ROX31F) * C(SA2XC10H8) * C(SO)
      W(ROX31B) = K(ROX31B) * C(SA2OHXC10)
      W(ROX32F) = K(ROX32F) * C(SA2XC10H8) * C(SOH)
      W(ROX32B) = K(ROX32B) * C(SH) * C(SA2OHXC10)
      W(ROX33F) = K(ROX33F) * C(SA2XXC10H) * C(SO2)
      W(ROX33B) = K(ROX33B) * C(SO) * C(SA2OXC10H)
      W(ROX34F) = K(ROX34F) * C(SA2YXC10H) * C(SO2)
      W(ROX34B) = K(ROX34B) * C(SO) * C(SA2OXC10H)
      W(ROX35) = K(ROX35) * C(SA2XXC10H) * C(SO2)
      W(ROX36) = K(ROX36) * C(SA2YXC10H) * C(SO2)
      W(ROX37F) = K(ROX37F) * C(SA2XXC10H) * C(SO)
      W(ROX37B) = K(ROX37B) * C(SA2OXC10H)
      W(ROX38F) = K(ROX38F) * C(SA2YXC10H) * C(SO)
      W(ROX38B) = K(ROX38B) * C(SA2OXC10H)
      W(ROX39F) = K(ROX39F) * C(SA2XXC10H) * C(SOH)
      W(ROX39B) = K(ROX39B) * C(SH) * C(SA2OXC10H)
      W(ROX40F) = K(ROX40F) * C(SA2YXC10H) * C(SOH)
      W(ROX40B) = K(ROX40B) * C(SH) * C(SA2OXC10H)
      W(ROX41F) = K(ROX41F) * C(SA2OHXC10)
      W(ROX41B) = K(ROX41B) * C(SCO) * C(SC9H8)
      W(ROX42F) = K(ROX42F) * C(SA2OHXC10)
      W(ROX42B) = K(ROX42B) * C(SH) * C(SA2OXC10H)
      W(ROX43F) = K(ROX43F) * C(SA2OHXC10) * C(SH)
      W(ROX43B) = K(ROX43B) * C(SH2) * C(SA2OXC10H)
      W(ROX44F) = K(ROX44F) * C(SA2OHXC10) * C(SOH)
      W(ROX44B) = K(ROX44B) * C(SH2O) * C(SA2OXC10H)
      W(ROX45F) = K(ROX45F) * C(SA2OHXC10) * C(SCH3)
      W(ROX45B) = K(ROX45B) * C(SCH4) * C(SA2OXC10H)
      W(ROX46F) = K(ROX46F) * C(SA2OXC10H)
      W(ROX46B) = K(ROX46B) * C(SCO) * C(SC9H7)
      W(ROX47) = K(ROX47) * C(SA2OXC10H) * C(SO)
      W(ROX48) = K(ROX48) * C(SA2OXC10H) * C(SO2)
      W(ROX50) = K(ROX50) * C(SA3XXC14H) * C(SO2)
      W(ROX51) = K(ROX51) * C(SA3YXC14H) * C(SO2)
      W(ROX52) = K(ROX52) * C(SA4XXC16H) * C(SO2)
      W(ROX53) = K(ROX53) * C(SA2R5XXC1) * C(SO2)
      W(ROX54) = K(ROX54) * C(SA3R5XXC1) * C(SO2)
      W(ROX60) = K(ROX60) * C(SA3XC14H1) * C(SOH)
      W(ROX61) = K(ROX61) * C(SA3XC14H1) * C(SOH)
      W(ROX62) = K(ROX62) * C(SA4XC16H1) * C(SOH)
      W(ROX63) = K(ROX63) * C(SA2R5XC12) * C(SOH)
      W(ROX64) = K(ROX64) * C(SA3R5XC16) * C(SOH)
      W(ROX65) = K(ROX65) * C(SA4R5XC18) * C(SOH)
      W(ROX66) = K(ROX66) * C(SFLTNXC16) * C(SOH)


      CDOT(SN2) = - W(RG28F) + W(RG28F) - W(RG28B) & 
	    + W(RG28B)
      CDOT(SH) = - W(R1F) + W(R1B) + W(R2F) & 
	    - W(R2B) + W(R3F) - W(R3B) & 
	    - 2 * W(R5F) + 2 * W(R5B) - 2 * W(R6F) & 
	    + 2 * W(R6B) - 2 * W(R7F) + 2 * W(R7B) & 
	    - 2 * W(R8F) + 2 * W(R8B) - W(R9F) & 
	    + W(R9B) - W(R10F) + W(R10B) & 
	    - W(R12F) + W(R12B) + W(R13F) & 
	    - W(R13B) - W(R15F) + W(R15B) & 
	    - W(R16F) + W(R16B) - W(R22F) & 
	    + W(R22B) - W(R23F) + W(R23B) & 
	    + W(R28F) - W(R28B) + W(R29F) & 
	    - W(R29B) - W(R32F) + W(R32B) & 
	    + W(R34F) - W(R34B) + W(R36F) & 
	    - W(R36B) + W(R37F) - W(R37B) & 
	    + W(RG01F) - W(RG01B) - W(RG03F)
      CDOT(SH) = CDOT(SH) + W(RG03B) + W(RG04F) - W(RG04B) & 
	    + W(RG05F) - W(RG05B) + W(RG06F) & 
	    - W(RG06B) + W(RG08F) - W(RG08B) & 
	    - W(RG13F) + W(RG13B) - W(RG14F) & 
	    + W(RG14B) + W(RG15F) - W(RG15B) & 
	    + W(RG16F) - W(RG16B) + W(RG18F) & 
	    - W(RG18B) + 2 * W(RG19) + W(RG21) & 
	    + W(RG23F) - W(RG23B) + W(RG25F) & 
	    - W(RG25B) + 2 * W(RG27) - W(RG30F) & 
	    + W(RG30B) + W(RG32F) - W(RG32B) & 
	    + W(RG33F) - W(RG33B) + W(RG34F) & 
	    - W(RG34B) + W(RG35F) - W(RG35B) & 
	    - W(RG43F) + W(RG43B) - W(RG44F) & 
	    + W(RG44B) - W(RG45F) + W(RG45B) & 
	    + W(RG50F) - W(RG50B) - W(RG51F)
      CDOT(SH) = CDOT(SH) + W(RG51B) + W(RG52F) - W(RG52B) & 
	    + W(RG53) + W(RG68F) - W(RG68B) & 
	    + W(RG69F) - W(RG69B) + W(RG72F) & 
	    - W(RG72B) + W(RG73F) - W(RG73B) & 
	    + W(RG74F) - W(RG74B) - W(RG75F) & 
	    + W(RG75B) - W(RG76F) + W(RG76F) & 
	    - W(RG76B) + W(RG76B) - W(RG77F) & 
	    + W(RG77B) - W(RG78F) + W(RG78B) & 
	    - W(RG79F) + W(RG79B) - W(RG83F) & 
	    + W(RG83B) - W(RG84F) + W(RG84B) & 
	    - W(RG85F) + W(RG85B) - W(RG86F) & 
	    + W(RG86B) - W(RG90F) + W(RG90B) & 
	    + W(RG93F) - W(RG93B) - W(RG96F) & 
	    + W(RG96B) - W(RG97F) + W(RG97B) & 
	    - W(RG104F) + W(RG104B) + W(RG106F)
      CDOT(SH) = CDOT(SH) - W(RG106B) + W(RG108F) - W(RG108B) & 
	    - W(RG109F) + W(RG109B) + W(RG110F) & 
	    - W(RG110B) - W(RG115F) + W(RG115B) & 
	    + W(RG116F) - W(RG116B) + W(RG120F) & 
	    - W(RG120B) + W(RG121F) - W(RG121B) & 
	    - W(RG123F) + W(RG123B) - W(RG124F) & 
	    + W(RG124B) - W(RG128F) + W(RG128F) & 
	    - W(RG128B) + W(RG128B) - W(RG129F) & 
	    + W(RG129B) - W(RG130F) + W(RG130B) & 
	    + W(RG136F) - W(RG136B) - W(RG141F) & 
	    + W(RG141B) - W(RG142F) + W(RG142B) & 
	    - W(RG147F) + W(RG147B) - W(RG148) & 
	    - W(RG155F) + W(RG155B) - W(RG156F) & 
	    + W(RG156B) + W(RG157F) - W(RG157B) & 
	    - W(RG164F) + W(RG164B) - W(RG165F)
      CDOT(SH) = CDOT(SH) + W(RG165B) + W(RG169F) - W(RG169B) & 
	    - W(RG174F) + W(RG174B) - W(RG180F) & 
	    + W(RG180B) - W(RG183F) + W(RG183B) & 
	    - W(RG185F) + W(RG185B) + W(RR004F) & 
	    - W(RR004B) - W(RR005F) + W(RR005B) & 
	    - W(RR006F) + W(RR006B) + W(RR015F) & 
	    - W(RR015B) - W(RR016F) + W(RR016B) & 
	    + W(RR017F) - W(RR017B) - W(RR018F) & 
	    + W(RR018B) + W(RR020F) - W(RR020B) & 
	    + W(RR031F) - W(RR031B) + W(RR033) & 
	    + W(RR041F) - W(RR041B) + W(RR042F) & 
	    - W(RR042B) + W(RR043F) - W(RR043B) & 
	    + W(RR044F) - W(RR044B) - W(RR200F) & 
	    + W(RR200B) - W(RR047) - W(RR053F) & 
	    + W(RR053B) - W(RR054F) + W(RR054B)
      CDOT(SH) = CDOT(SH) - W(RR055F) + W(RR055B) - W(RR056F) & 
	    + W(RR056B) + W(RR061F) - W(RR061B) & 
	    + W(RR069F) - W(RR069B) + W(RR070F) & 
	    - W(RR070B) - W(RR073F) + W(RR073F) & 
	    - W(RR073B) + W(RR073B) - W(RR074F) & 
	    + W(RR074B) - W(RR075F) + W(RR075B) & 
	    - W(RR076F) + W(RR076B) - W(RR077F) & 
	    + W(RR077B) - W(RR078F) + W(RR078B) & 
	    - W(RR083F) + W(RR083B) - W(RR095) & 
	    - W(RR103F) + W(RR103B) - W(RR119F) & 
	    + W(RR119B) + W(RR121) - W(RR125F) & 
	    + W(RR125B) + W(RR127) + W(RR134) & 
	    + W(RR139F) - W(RR139B) - W(RR141F) & 
	    + W(RR141B) - W(RR144F) + W(RR144B) & 
	    - W(RR149F) + W(RR149B) - W(RR153F)
      CDOT(SH) = CDOT(SH) + W(RR153B) - W(RH02F) + W(RH02B) & 
	    - W(RH03F) + W(RH03B) + 2 * W(RH05) & 
	    - W(RH09F) + W(RH09B) - W(RH10F) & 
	    + W(RH10B) - W(RH14F) + W(RH14F) & 
	    - W(RH14B) + W(RH14B) - W(RH15F) & 
	    + W(RH15B) - W(RH16F) + W(RH16B) & 
	    - W(RH17F) + W(RH17B) - W(RH18F) & 
	    + W(RH18B) - W(RH19F) + W(RH19B) & 
	    - W(RH20F) + W(RH20B) - W(RH28F) & 
	    + W(RH28B) - W(RH29F) + W(RH29B) & 
	    + W(RH37F) - W(RH37B) + W(RH38F) & 
	    - W(RH38B) + W(RH39F) - W(RH39B) & 
	    + W(RH40F) - W(RH40B) + W(RB05F) & 
	    - W(RB05B) + W(RB06F) - W(RB06B) & 
	    + W(RB10F) - W(RB10B) + W(RB11F)
      CDOT(SH) = CDOT(SH) - W(RB11B) + W(RB14F) - W(RB14B) & 
	    + W(RB15F) - W(RB15B) - W(RB16F) & 
	    + W(RB16B) - W(RB17F) + W(RB17B) & 
	    + W(RB28) - W(RB32F) + W(RB32B) & 
	    - W(RB33F) + W(RB33B) - W(RB35F) & 
	    + W(RB35F) - W(RB35B) + W(RB35B) & 
	    - W(RB36F) + W(RB36B) - W(RB44F) & 
	    + W(RB44B) - W(RB45F) + W(RB45B) & 
	    - W(RHP02) - W(RHP08) + W(RHP14) & 
	    - W(RHP15) - W(RHP23F) + W(RHP23B) & 
	    + W(RHP30) - W(RHP31) - W(RHP34) & 
	    + W(RHP39) - W(RHP40) - W(RHP44) & 
	    - W(RHP49) - W(RHP59) + W(RHP64F) & 
	    - W(RHP64B) - W(RHP67) + W(RHP79) & 
	    - W(RIC03) - W(RIC18) + W(RIC19)
      CDOT(SH) = CDOT(SH) - W(RIC24) + W(RIC25) + W(RIC29) & 
	    + W(RIC33) - W(RIC34) + W(RIC38F) & 
	    - W(RIC38B) - W(RIC40) - W(RIC41F) & 
	    + W(RIC41B) + W(RIC50) - W(RIC59) & 
	    + W(RIC65) - W(RIC66) - W(RIC69) & 
	    - W(RP000F) + W(RP000F) - W(RP000B) & 
	    + W(RP000B) + W(RP001F) - W(RP001B) & 
	    + W(RP002F) - W(RP002B) + W(RP004F) & 
	    - W(RP004B) + W(RP006F) - W(RP006B) & 
	    + W(RP007F) - W(RP007B) + 2 * W(RP008) & 
	    + W(RP011F) - W(RP011B) + W(RP014F) & 
	    - W(RP014B) + W(RP016F) - W(RP016B) & 
	    - W(RK017F) + W(RK017B) + W(RK020F) & 
	    - W(RK020B) - W(RK021F) + W(RK021B) & 
	    + W(RP023F) - W(RP023B) - W(RK024F)
      CDOT(SH) = CDOT(SH) + W(RK024B) + W(RP026F) - W(RP026B) & 
	    - W(RP027F) + W(RP027B) + W(RK102F) & 
	    - W(RK102B) + W(RP104F) - W(RP104B) & 
	    + W(RP105F) - W(RP105B) + W(RP106F) & 
	    - W(RP106B) + W(RP107F) - W(RP107B) & 
	    + W(RP108F) - W(RP108B) - W(RK109F) & 
	    + W(RK109B) + W(RP111F) - W(RP111B) & 
	    - W(RK112F) + W(RK112B) + W(RP116F) & 
	    - W(RP116B) + W(RP117F) - W(RP117B) & 
	    + W(RK122F) - W(RK122B) - W(RK123F) & 
	    + W(RK123B) + W(RK125F) - W(RK125B) & 
	    - W(RK126F) + W(RK126B) + W(RP128F) & 
	    - W(RP128B) - W(RK129F) + W(RK129B) & 
	    + W(RP131F) - W(RP131B) - W(RK132F) & 
	    + W(RK132B) + W(RK200F) - W(RK200B)
      CDOT(SH) = CDOT(SH) - W(RP201F) + W(RP201F) - W(RP201B) & 
	    + W(RP201B) + W(RP202F) - W(RP202B) & 
	    - W(RK203F) + W(RK203B) + W(RP206F) & 
	    - W(RP206B) + W(RP209F) - W(RP209B) & 
	    - W(RK210F) + W(RK210B) + W(RK212F) & 
	    - W(RK212B) - W(RK213F) + W(RK213B) & 
	    + W(RP301F) - W(RP301B) + W(RP304F) & 
	    - W(RP304B) - W(RP305F) + W(RP305B) & 
	    + W(RP405F) - W(RP405B) + W(RP406F) & 
	    - W(RP406B) + W(RP407F) - W(RP407B) & 
	    + W(RP410F) - W(RP410B) + W(RP411F) & 
	    - W(RP411B) + W(RP412F) - W(RP412B) & 
	    + W(RP413F) - W(RP413B) + W(RP414F) & 
	    - W(RP414B) + W(RP415F) - W(RP415B) & 
	    + W(RP416F) - W(RP416B) + W(RP417F)
      CDOT(SH) = CDOT(SH) - W(RP417B) + W(RP419F) - W(RP419B) & 
	    - W(RK420F) + W(RK420B) + W(RP422F) & 
	    - W(RP422B) - W(RK423F) + W(RK423B) & 
	    + W(RP502F) - W(RP502B) + W(RP503F) & 
	    - W(RP503B) + W(RP505F) - W(RP505B) & 
	    + W(RP506F) - W(RP506B) - W(RK507F) & 
	    + W(RK507B) + W(RK600F) - W(RK600B) & 
	    + W(RP601F) - W(RP601B) - W(RK602F) & 
	    + W(RK602B) + W(RK700F) - W(RK700B) & 
	    + W(RK701F) - W(RK701B) + W(RP800) & 
	    + W(RP801) + W(RCP01F) - W(RCP01B) & 
	    - W(RCP02F) + W(RCP02B) - W(RCP03F) & 
	    + W(RCP03B) - W(RCP04) + W(RCP12F) & 
	    - W(RCP12B) + 2 * W(RCP17) + 2 * W(RCP18) & 
	    + W(RCP19F) - W(RCP19B) + W(RCP22)
      CDOT(SH) = CDOT(SH) + W(RCP23F) - W(RCP23B) - W(RCP25) & 
	    - W(RCP29) - W(RCP32F) + W(RCP32B) & 
	    + W(RI01F) - W(RI01B) - W(RI02F) & 
	    + W(RI02B) + W(RI03F) - W(RI03B) & 
	    + 2 * W(RI12) + 2 * W(RI17) + 2 * W(RI18) & 
	    + W(RI19F) - W(RI19B) + 2 * W(RI22) & 
	    + 2 * W(RI23) - W(RI32) - W(RT01F) & 
	    + W(RT01B) + W(RT02F) - W(RT02B) & 
	    - W(RT04F) + W(RT04B) - W(RT07F) & 
	    + W(RT07B) + W(RT10F) - W(RT10B) & 
	    + W(RT13F) - W(RT13B) + 2 * W(RT20) & 
	    - W(RT23F) + W(RT23B) - W(RT24F) & 
	    + W(RT24B) + W(RT29F) - W(RT29B) & 
	    - W(RT32F) + W(RT32B) + W(RT36) & 
	    - W(RT37) + W(RT43F) - W(RT43B)
      CDOT(SH) = CDOT(SH) - W(RT44F) + W(RT44B) + W(RT47) & 
	    + W(RT50F) - W(RT50B) - W(RT51F) & 
	    + W(RT51B) + W(RT56F) - W(RT56B) & 
	    + W(RT59) - W(RE01F) + W(RE01B) & 
	    - W(RE04F) + W(RE04B) - W(RE06F) & 
	    + W(RE06B) + W(RE12F) - W(RE12B) & 
	    - W(RE13F) + W(RE13B) + W(RXY00F) & 
	    - W(RXY00B) - W(RXY02F) + W(RXY02B) & 
	    - W(RXY09F) + W(RXY09B) + 2 * W(RXY11) & 
	    + W(RXY12) - W(RXY13F) + W(RXY13B) & 
	    + W(RXY14F) - W(RXY14B) + W(RXY201) & 
	    + 2 * W(RXY22) + W(RXY23F) - W(RXY23B) & 
	    + W(RXY25) - W(RXY26F) + W(RXY26B) & 
	    - W(RXY33) - W(RXY39F) + W(RXY39B) & 
	    - W(RXY40F) + W(RXY40B) + W(RXY42)
      CDOT(SH) = CDOT(SH) + W(RXY43F) - W(RXY43B) + W(RXY48) & 
	    + W(RXY50) - W(RXY51) + W(RXY51) & 
	    - W(RXY57F) + W(RXY57B) + W(RXY58) & 
	    - W(RN01F) + W(RN01B) + W(RN04F) & 
	    - W(RN04B) - W(RN06F) + W(RN06B) & 
	    - W(RN08F) + W(RN08B) + 2 * W(RN14) & 
	    + W(RN17) + 2 * W(RN18) + W(RN22F) & 
	    - W(RN22B) - W(RN24F) + W(RN24B) & 
	    + W(RN28) - W(RN29) + W(ROX00F) & 
	    - W(ROX00B) + W(ROX01F) - W(ROX01B) & 
	    - W(ROX03F) + W(ROX03B) + W(ROX05F) & 
	    - W(ROX05B) + W(ROX99F) - W(ROX99B) & 
	    + W(ROX07F) - W(ROX07B) + W(ROX12F) & 
	    - W(ROX12B) + W(ROX14F) - W(ROX14B) & 
	    + W(ROX18F) - W(ROX18B) - W(ROX19F)
      CDOT(SH) = CDOT(SH) + W(ROX19B) + W(ROX24F) - W(ROX24B) & 
	    - W(ROX27F) + W(ROX27B) + W(ROX30F) & 
	    - W(ROX30B) + W(ROX32F) - W(ROX32B) & 
	    + W(ROX35) + W(ROX36) + W(ROX39F) & 
	    - W(ROX39B) + W(ROX40F) - W(ROX40B) & 
	    + W(ROX42F) - W(ROX42B) - W(ROX43F) & 
	    + W(ROX43B) + W(ROX47)
      CDOT(SO2) = - W(R1F) + W(R1B) + W(R11F) & 
	    - W(R11B) - W(R12F) + W(R12B) & 
	    - W(R13F) + W(R13B) + W(R17F) & 
	    - W(R17B) + W(R18F) - W(R18B) & 
	    + W(R19F) - W(R19B) + W(R20F) & 
	    - W(R20B) + W(R21F) - W(R21B) & 
	    - W(R30F) + W(R30B) - W(R38F) & 
	    + W(R38B) - W(RG02F) + W(RG02B) & 
	    - W(RG09F) + W(RG09B) - W(RG19) & 
	    - W(RG20F) + W(RG20B) - W(RG21) & 
	    - W(RG35F) + W(RG35B) - W(RG36F) & 
	    + W(RG36B) - W(RG48F) + W(RG48B) & 
	    - W(RG58F) + W(RG58B) - W(RG59F) & 
	    + W(RG59B) - W(RG60F) + W(RG60B) & 
	    + W(RG62) + W(RG63) + W(RG66F)
      CDOT(SO2) = CDOT(SO2) - W(RG66B) - W(RG82F) + W(RG82B) & 
	    - W(RG89F) + W(RG89B) - W(RG107F) & 
	    + W(RG107B) - W(RG111F) + W(RG111B) & 
	    - W(RG133F) + W(RG133B) - W(RG134F) & 
	    + W(RG134B) - W(RG135F) + W(RG135B) & 
	    - W(RG139) - W(RG140) - W(RG150) & 
	    - W(RG170F) + W(RG170B) - W(RG171F) & 
	    + W(RG171B) - W(RG184F) + W(RG184B) & 
	    - W(RR002F) + W(RR002B) - W(RR003F) & 
	    + W(RR003B) - W(RR021F) + W(RR021B) & 
	    - W(RR032F) + W(RR032B) - W(RR033) & 
	    + W(RR035F) - W(RR035B) - W(RR041F) & 
	    + W(RR041B) - W(RR049) - W(RR062F) & 
	    + W(RR062B) + W(RR064F) - W(RR064B) & 
	    + W(RR065F) - W(RR065B) - W(RR066F)
      CDOT(SO2) = CDOT(SO2) + W(RR066B) - W(RR109F) + W(RR109B) & 
	    - W(RR110F) + W(RR110B) - W(RR111) & 
	    - W(RR112F) + W(RR112B) + W(RR116F) & 
	    - W(RR116B) - W(RR131F) + W(RR131B) & 
	    - W(RR132F) + W(RR132B) - W(RR133) & 
	    - W(RR134) - W(RR135) - W(RR136F) & 
	    + W(RR136B) - W(RR137F) + W(RR137B) & 
	    - W(RR138) - W(RH01F) + W(RH01B) & 
	    - W(RH07F) + W(RH07B) - W(RH23F) & 
	    + W(RH23B) - W(RH24F) + W(RH24B) & 
	    - W(RH26F) + W(RH26B) - W(RH27F) & 
	    + W(RH27B) + W(RB40F) - W(RB40B) & 
	    - W(RB41F) + W(RB41B) - W(RB42) & 
	    - W(RB43F) + W(RB43B) + W(RB48F) & 
	    - W(RB48B) - W(RB50F) + W(RB50B)
      CDOT(SO2) = CDOT(SO2) - W(RHP05) + W(RHP09) - W(RHP46) & 
	    + W(RHP50) - W(RHP61) - W(RHP68) & 
	    - W(RIC06) - W(RIC51) - W(RIC54) & 
	    - W(RIC55) - W(RIC56) - W(RIC68) & 
	    - W(RCP07F) + W(RCP07B) - W(RCP20F) & 
	    + W(RCP20B) - W(RCP28) - W(RCP30) & 
	    - W(RI07F) + W(RI07B) - W(RI20F) & 
	    + W(RI20B) - W(RT06F) + W(RT06B) & 
	    - W(RT21F) + W(RT21B) - W(RT22F) & 
	    + W(RT22B) - W(RT35F) + W(RT35B) & 
	    - W(RT40) - W(RT58F) + W(RT58B) & 
	    - W(RT59) - W(RT60) - W(RE16F) & 
	    + W(RE16B) - W(RE30) - W(RE31) & 
	    + W(RE32) - W(RE36) - W(RST12F) & 
	    + W(RST12B) - W(RST13) - W(RST14F)
      CDOT(SO2) = CDOT(SO2) + W(RST14B) - W(RXY05F) + W(RXY05B) & 
	    - W(RXY18F) + W(RXY18B) - W(RXY19F) & 
	    + W(RXY19B) - W(RXY29F) + W(RXY29B) & 
	    - W(RXY36) - W(RXY47F) + W(RXY47B) & 
	    - W(RXY54) + W(RXY54) - W(RN11F) & 
	    + W(RN11B) - W(RN20F) + W(RN20B) & 
	    - W(RN21F) + W(RN21B) - W(RN27F) & 
	    + W(RN27B) - W(RN32) - W(ROX06F) & 
	    + W(ROX06B) - W(ROX11F) + W(ROX11B) & 
	    - W(ROX12F) + W(ROX12B) - W(ROX22F) & 
	    + W(ROX22B) - W(ROX25F) + W(ROX25B) & 
	    - W(ROX33F) + W(ROX33B) - W(ROX34F) & 
	    + W(ROX34B) - W(ROX35) - W(ROX36) & 
	    - W(ROX48) - W(ROX50) - W(ROX51) & 
	    - W(ROX52) - W(ROX53) - W(ROX54)
      CDOT(SO) = W(R1F) - W(R1B) - W(R2F) & 
	    + W(R2B) + W(R4F) - W(R4B) & 
	    - W(R10F) + W(R10B) - 2 * W(R11F) & 
	    + 2 * W(R11B) + W(R15F) - W(R15B) & 
	    - W(R17F) + W(R17B) - W(R24F) & 
	    + W(R24B) - W(R27F) + W(R27B) & 
	    + W(R30F) - W(R30B) - W(R33F) & 
	    + W(R33B) - W(R34F) + W(R34B) & 
	    + W(RG02F) - W(RG02B) - W(RG04F) & 
	    + W(RG04B) + W(RG09F) - W(RG09B) & 
	    - W(RG15F) + W(RG15B) + W(RG20F) & 
	    - W(RG20B) - W(RG31F) + W(RG31B) & 
	    - W(RG32F) + W(RG32B) - W(RG46F) & 
	    + W(RG46B) - W(RG52F) + W(RG52B) & 
	    - W(RG53) + W(RG58F) - W(RG58B)
      CDOT(SO) = CDOT(SO) - W(RG80F) + W(RG80B) - W(RG87F) & 
	    + W(RG87B) - W(RG91F) + W(RG91B) & 
	    - W(RG98F) + W(RG98B) - W(RG99F) & 
	    + W(RG99B) - W(RG105F) + W(RG105B) & 
	    - W(RG110F) + W(RG110B) - W(RG116F) & 
	    + W(RG116B) - W(RG117F) + W(RG117B) & 
	    + W(RG118F) - W(RG118B) - W(RG125F) & 
	    + W(RG125B) - W(RG126F) + W(RG126B) & 
	    - W(RG131F) + W(RG131B) + W(RG134F) & 
	    - W(RG134B) - W(RG138F) + W(RG138B) & 
	    - W(RG146F) + W(RG146B) - W(RG149) & 
	    - W(RG157F) + W(RG157B) - W(RG158F) & 
	    + W(RG158B) - W(RG159F) + W(RG159B) & 
	    - W(RG167F) + W(RG167B) - W(RG175F) & 
	    + W(RG175B) - W(RG179F) + W(RG179B)
      CDOT(SO) = CDOT(SO) - W(RG186F) + W(RG186B) - W(RR019F) & 
	    + W(RR019B) + W(RR021F) - W(RR021B) & 
	    - W(RR039F) + W(RR039B) - W(RR048) & 
	    - W(RR061F) + W(RR061B) - W(RR079F) & 
	    + W(RR079B) - W(RR087F) + W(RR087B) & 
	    - W(RR088F) + W(RR088B) - W(RR089F) & 
	    + W(RR089B) - W(RR096) - W(RR113F) & 
	    + W(RR113B) - W(RR120F) + W(RR120B) & 
	    - W(RR126F) + W(RR126B) + W(RR133) & 
	    + W(RR134) - W(RR142F) + W(RR142B) & 
	    - W(RR143F) + W(RR143B) - W(RR145F) & 
	    + W(RR145B) - W(RR150F) + W(RR150B) & 
	    - W(RR154F) + W(RR154B) - W(RH08F) & 
	    + W(RH08B) - W(RH25F) + W(RH25B) & 
	    - W(RH34F) + W(RH34B) - W(RH35F)
      CDOT(SO) = CDOT(SO) + W(RH35B) - W(RH36F) + W(RH36B) & 
	    + W(RB18F) - W(RB18B) - W(RB19F) & 
	    + W(RB19B) - W(RB28) - W(RB29F) & 
	    + W(RB29B) - W(RB30F) + W(RB30B) & 
	    - W(RB99F) + W(RB99B) + W(RB42) & 
	    - W(RB51F) + W(RB51B) - W(RHP03) & 
	    - W(RHP35) - W(RHP51) - W(RHP52) & 
	    - W(RHP53) - W(RHP54) - W(RHP86) & 
	    - W(RIC04) - W(RIC42) - W(RIC47) & 
	    - W(RIC48) - W(RIC53) - W(RIC60) & 
	    - W(RIC70) - W(RCP05F) + W(RCP05B) & 
	    - W(RCP12F) + W(RCP12B) - W(RCP13) & 
	    - W(RCP19F) + W(RCP19B) - W(RCP26) & 
	    - W(RCP33F) + W(RCP33B) - W(RI05F) & 
	    + W(RI05B) - W(RI12) - W(RI19F)
      CDOT(SO) = CDOT(SO) + W(RI19B) - W(RT11F) + W(RT11B) & 
	    - W(RT12F) + W(RT12B) - W(RT13F) & 
	    + W(RT13B) - W(RT17F) + W(RT17B) & 
	    - W(RT25F) + W(RT25B) - W(RT33F) & 
	    + W(RT33B) - W(RT38) - W(RT45F) & 
	    + W(RT45B) - W(RT52F) + W(RT52B) & 
	    - W(RT55F) + W(RT55B) + W(RT58F) & 
	    - W(RT58B) - W(RE07F) + W(RE07B) & 
	    - W(RE17F) + W(RE17B) - W(RE18F) & 
	    + W(RE18B) - W(RST02F) + W(RST02B) & 
	    - W(RST10F) + W(RST10B) - W(RST11F) & 
	    + W(RST11B) + W(RST13) - W(RXY03F) & 
	    + W(RXY03B) - W(RXY11) - W(RXY14F) & 
	    + W(RXY14B) - W(RXY15F) + W(RXY15B) & 
	    - W(RXY16F) + W(RXY16B) - W(RXY27F)
      CDOT(SO) = CDOT(SO) + W(RXY27B) - W(RXY34) - W(RXY43F) & 
	    + W(RXY43B) - W(RXY44) - W(RXY45F) & 
	    + W(RXY45B) - W(RXY52) + W(RXY52) & 
	    - W(RN09F) + W(RN09B) - W(RN14) & 
	    - W(RN15) - W(RN16F) + W(RN16B) & 
	    - W(RN25F) + W(RN25B) - W(RN30) & 
	    - W(ROX99F) + W(ROX99B) - W(ROX07F) & 
	    + W(ROX07B) - W(ROX08F) + W(ROX08B) & 
	    - W(ROX09F) + W(ROX09B) - W(ROX10F) & 
	    + W(ROX10B) + W(ROX11F) - W(ROX11B) & 
	    - W(ROX13F) + W(ROX13B) - W(ROX24F) & 
	    + W(ROX24B) - W(ROX28) - W(ROX30F) & 
	    + W(ROX30B) - W(ROX31F) + W(ROX31B) & 
	    + W(ROX33F) - W(ROX33B) + W(ROX34F) & 
	    - W(ROX34B) - W(ROX37F) + W(ROX37B)
      CDOT(SO) = CDOT(SO) - W(ROX38F) + W(ROX38B) - W(ROX47)
      CDOT(SOH) = W(R1F) - W(R1B) + W(R2F) & 
	    - W(R2B) - W(R3F) + W(R3B) & 
	    - 2 * W(R4F) + 2 * W(R4B) - W(R9F) & 
	    + W(R9B) + W(R10F) - W(R10B) & 
	    - 2 * W(R14F) + 2 * W(R14B) + 2 * W(R16F) & 
	    - 2 * W(R16B) + W(R17F) - W(R17B) & 
	    - W(R18F) + W(R18B) - W(R19F) & 
	    + W(R19B) + W(R23F) - W(R23B) & 
	    + W(R24F) - W(R24B) - W(R25F) & 
	    + W(R25B) - W(R26F) + W(R26B) & 
	    - W(R28F) + W(R28B) - W(R29F) & 
	    + W(R29B) + W(R31F) - W(R31B) & 
	    + W(R33F) - W(R33B) - W(R35F) & 
	    + W(R35B) - W(RG01F) + W(RG01B) & 
	    - W(RG05F) + W(RG05B) - W(RG16F)
      CDOT(SOH) = CDOT(SOH) + W(RG16B) - W(RG17F) + W(RG17B) & 
	    + W(RG21) + W(RG22F) - W(RG22B) & 
	    - W(RG33F) + W(RG33B) + W(RG35F) & 
	    - W(RG35B) + W(RG46F) - W(RG46B) & 
	    - W(RG47F) + W(RG47B) - W(RG54F) & 
	    + W(RG54B) - W(RG55F) + W(RG55B) & 
	    - W(RG56) - W(RG57F) + W(RG57B) & 
	    + W(RG59F) - W(RG59B) + W(RG63) & 
	    + W(RG64) + W(RG65F) - W(RG65B) & 
	    + W(RG78F) - W(RG78B) + W(RG80F) & 
	    - W(RG80B) - W(RG81F) + W(RG81B) & 
	    + W(RG85F) - W(RG85B) + W(RG87F) & 
	    - W(RG87B) - W(RG88F) + W(RG88B) & 
	    + W(RG91F) - W(RG91B) - W(RG92F) & 
	    + W(RG92B) + W(RG98F) - W(RG98B)
      CDOT(SOH) = CDOT(SOH) + W(RG99F) - W(RG99B) - W(RG100F) & 
	    + W(RG100B) - W(RG101F) + W(RG101B) & 
	    - W(RG106F) + W(RG106B) + W(RG111F) & 
	    - W(RG111B) - W(RG118F) + W(RG118B) & 
	    - W(RG119F) + W(RG119B) - W(RG120F) & 
	    + W(RG120B) - W(RG121F) + W(RG121B) & 
	    - W(RG122F) + W(RG122B) + W(RG125F) & 
	    - W(RG125B) - W(RG127F) + W(RG127B) & 
	    - W(RG132F) + W(RG132B) + W(RG139) & 
	    + W(RG140) - W(RG143F) + W(RG143B) & 
	    - W(RG144F) + W(RG144B) + W(RG146F) & 
	    - W(RG146B) + W(RG149) - W(RG151) & 
	    - W(RG160F) + W(RG160B) - W(RG161F) & 
	    + W(RG161B) + W(RG175F) - W(RG175B) & 
	    - W(RG176F) + W(RG176B) - W(RG181F)
      CDOT(SOH) = CDOT(SOH) + W(RG181B) + W(RG186F) - W(RG186B) & 
	    - W(RG187F) + W(RG187B) - W(RR020F) & 
	    + W(RR020B) - W(RR023F) + W(RR023B) & 
	    - W(RR024F) + W(RR024B) - W(RR025F) & 
	    + W(RR025B) + W(RR037F) - W(RR037B) & 
	    - W(RR040F) + W(RR040B) + W(RR048) & 
	    - W(RR050) - W(RR057F) + W(RR057B) & 
	    - W(RR058F) + W(RR058B) - W(RR059F) & 
	    + W(RR059B) - W(RR060F) + W(RR060B) & 
	    + W(RR063F) - W(RR063B) + W(RR079F) & 
	    - W(RR079B) - W(RR080F) + W(RR080B) & 
	    - W(RR084F) + W(RR084B) - W(RR092F) & 
	    + W(RR092B) - W(RR093F) + W(RR093B) & 
	    - W(RR094F) + W(RR094B) + W(RR096) & 
	    - W(RR097) - W(RR104F) + W(RR104B)
      CDOT(SOH) = CDOT(SOH) + W(RR110F) - W(RR110B) + W(RR111) & 
	    - W(RR114F) + W(RR114B) + W(RR117F) & 
	    - W(RR117B) - W(RR121) + W(RR122F) & 
	    - W(RR122B) - W(RR127) + W(RR128F) & 
	    - W(RR128B) + W(RR145F) - W(RR145B) & 
	    - W(RR146F) + W(RR146B) + W(RR150F) & 
	    - W(RR150B) - W(RR151F) + W(RR151B) & 
	    + W(RR154F) - W(RR154B) - W(RR155F) & 
	    + W(RR155B) - W(RH11F) + W(RH11B) & 
	    - W(RH12F) + W(RH12B) - W(RH21F) & 
	    + W(RH21B) - W(RH22F) + W(RH22B) & 
	    - W(RH30F) + W(RH30B) - W(RH31F) & 
	    + W(RH31B) - W(RB18F) + W(RB18B) & 
	    + W(RB19F) - W(RB19B) - W(RB20F) & 
	    + W(RB20B) - W(RB21F) + W(RB21B)
      CDOT(SOH) = CDOT(SOH) - W(RB31F) + W(RB31B) - W(RB37F) & 
	    + W(RB37B) - W(RB46F) + W(RB46B) & 
	    + W(RHP03) - W(RHP04) + W(RHP16) & 
	    - W(RHP24) - W(RHP25) - W(RHP26) & 
	    + W(RHP35) - W(RHP36) - W(RHP45) & 
	    - W(RHP55) - W(RHP58) - W(RHP56) & 
	    - W(RHP57) - W(RHP60) + W(RHP72) & 
	    + W(RHP76) - W(RHP87) - W(RHP88) & 
	    + W(RIC04) - W(RIC05) + W(RIC09) & 
	    + W(RIC13) - W(RIC26) + W(RIC27) & 
	    + W(RIC28) + W(RIC35) + W(RIC42) & 
	    - W(RIC43) + W(RIC46) + W(RIC57) & 
	    + W(RIC54) + W(RIC56) + W(RIC60) & 
	    - W(RIC61) + W(RIC70) - W(RIC71) & 
	    - W(RP018F) + W(RP018B) - W(RP022F)
      CDOT(SOH) = CDOT(SOH) + W(RP022B) - W(RP025F) + W(RP025B) & 
	    - W(RP028F) + W(RP028B) - W(RK110F) & 
	    + W(RK110B) - W(RK113F) + W(RK113B) & 
	    - W(RP124F) + W(RP124B) - W(RP127F) & 
	    + W(RP127B) - W(RP130F) + W(RP130B) & 
	    - W(RP133F) + W(RP133B) - W(RK204F) & 
	    + W(RK204B) - W(RP211F) + W(RP211B) & 
	    - W(RP214F) + W(RP214B) - W(RP306F) & 
	    + W(RP306B) - W(RK421F) + W(RK421B) & 
	    - W(RK424F) + W(RK424B) - W(RK508F) & 
	    + W(RK508B) - W(RK603F) + W(RK603B) & 
	    + W(RCP05F) - W(RCP05B) - W(RCP06F) & 
	    + W(RCP06B) - W(RCP14) - W(RCP15) & 
	    + W(RCP20F) - W(RCP20B) + W(RCP21F) & 
	    - W(RCP21B) - W(RCP22) + W(RCP26)
      CDOT(SOH) = CDOT(SOH) - W(RCP27) + W(RI05F) - W(RI05B) & 
	    - W(RI06F) + W(RI06B) - W(RI15) & 
	    + W(RI20F) - W(RI20B) - W(RI22) & 
	    - W(RT08F) + W(RT08B) - W(RT09F) & 
	    + W(RT09B) - W(RT10F) + W(RT10B) & 
	    + W(RT11F) - W(RT11B) - W(RT18F) & 
	    + W(RT18B) + W(RT19F) - W(RT19B) & 
	    + W(RT21F) - W(RT21B) + W(RT25F) & 
	    - W(RT25B) - W(RT26F) + W(RT26B) & 
	    + W(RT33F) - W(RT33B) - W(RT34F) & 
	    + W(RT34B) + W(RT38) - W(RT39) & 
	    + W(RT45F) - W(RT45B) - W(RT46F) & 
	    + W(RT46B) + W(RT52F) - W(RT52B) & 
	    - W(RT53F) + W(RT53B) - W(RT56F) & 
	    + W(RT56B) + W(RT57F) - W(RT57B)
      CDOT(SOH) = CDOT(SOH) - W(RE05F) + W(RE05B) + W(RE07F) & 
	    - W(RE07B) - W(RE08F) + W(RE08B) & 
	    - W(RE14F) + W(RE14B) + W(RE19) & 
	    + W(RE31) + W(RE34) + W(RE36) & 
	    + W(RE37) - W(RST04F) + W(RST04B) & 
	    - W(RST05F) + W(RST05B) - W(RST06F) & 
	    + W(RST06B) + W(RST10F) - W(RST10B) & 
	    - W(RST00F) + W(RST00B) + W(RXY03F) & 
	    - W(RXY03B) - W(RXY04F) + W(RXY04B) & 
	    - W(RXY10F) + W(RXY10B) - W(RXY17) & 
	    + W(RXY18F) - W(RXY18B) + W(RXY201) & 
	    + W(RXY202) + W(RXY203) + W(RXY27F) & 
	    - W(RXY27B) - W(RXY28F) + W(RXY28B) & 
	    + W(RXY34) - W(RXY35) - W(RXY41F) & 
	    + W(RXY41B) - W(RXY42) - W(RXY46)
      CDOT(SOH) = CDOT(SOH) + W(RXY47F) - W(RXY47B) + W(RXY48) & 
	    - W(RXY53) + W(RXY53) - W(RXY58) & 
	    - W(RN02F) + W(RN02B) + W(RN09F) & 
	    - W(RN09B) - W(RN10F) + W(RN10B) & 
	    - W(RN17) + W(RN18F) - W(RN18B) & 
	    + W(RN20F) - W(RN20B) + W(RN25F) & 
	    - W(RN25B) - W(RN26F) + W(RN26B) & 
	    + W(RN30) - W(RN31) - W(ROX04F) & 
	    + W(ROX04B) - W(ROX05F) + W(ROX05B) & 
	    + W(ROX10F) - W(ROX10B) - W(ROX14F) & 
	    + W(ROX14B) + W(ROX15F) - W(ROX15B) & 
	    - W(ROX20F) + W(ROX20B) + W(ROX25F) & 
	    - W(ROX25B) - W(ROX32F) + W(ROX32B) & 
	    - W(ROX39F) + W(ROX39B) - W(ROX40F) & 
	    + W(ROX40B) - W(ROX44F) + W(ROX44B)
      CDOT(SOH) = CDOT(SOH) + W(ROX48) - W(ROX60) - W(ROX61) & 
	    - W(ROX62) - W(ROX63) - W(ROX64) & 
	    - W(ROX65) - W(ROX66)
      CDOT(SH2) = - W(R2F) + W(R2B) - W(R3F) & 
	    + W(R3B) + W(R5F) - W(R5B) & 
	    - W(R6F) + 2 * W(R6F) - 2 * W(R6B) & 
	    + W(R6B) + W(R7F) - W(R7B) & 
	    + W(R8F) - W(R8B) - W(R13F) & 
	    + W(R13B) + W(R22F) - W(R22B) & 
	    + W(R32F) - W(R32B) + W(RG03F) & 
	    - W(RG03B) - W(RG06F) + W(RG06B) & 
	    - W(RG07F) + W(RG07B) - W(RG12F) & 
	    + W(RG12B) - W(RG18F) + W(RG18B) & 
	    + W(RG26F) - W(RG26B) + W(RG30F) & 
	    - W(RG30B) + W(RG31F) - W(RG31B) & 
	    - W(RG34F) + W(RG34B) + W(RG39) & 
	    + W(RG45F) - W(RG45B) + W(RG53) & 
	    + W(RG56) + W(RG77F) - W(RG77B)
      CDOT(SH2) = CDOT(SH2) + W(RG84F) - W(RG84B) + W(RG90F) & 
	    - W(RG90B) + W(RG96F) - W(RG96B) & 
	    + W(RG97F) - W(RG97B) - W(RG108F) & 
	    + W(RG108B) + W(RG123F) - W(RG123B) & 
	    + W(RG130F) - W(RG130B) + W(RG142F) & 
	    - W(RG142B) + W(RG147F) - W(RG147B) & 
	    + W(RG148) + W(RG154F) - W(RG154B) & 
	    + W(RG156F) - W(RG156B) + W(RG165F) & 
	    - W(RG165B) + W(RG174F) - W(RG174B) & 
	    + W(RG185F) - W(RG185B) + W(RR047) & 
	    + W(RR054F) - W(RR054B) + W(RR078F) & 
	    - W(RR078B) + W(RR083F) - W(RR083B) & 
	    + W(RR095) + W(RR103F) - W(RR103B) & 
	    + W(RR114F) - W(RR114B) + W(RR119F) & 
	    - W(RR119B) + W(RR125F) - W(RR125B)
      CDOT(SH2) = CDOT(SH2) + W(RR144F) - W(RR144B) + W(RR149F) & 
	    - W(RR149B) + W(RR153F) - W(RR153B) & 
	    + W(RH03F) - W(RH03B) - W(RH04F) & 
	    + W(RH04B) + W(RH06F) - W(RH06B) & 
	    + W(RH19F) - W(RH19B) + W(RH20F) & 
	    - W(RH20B) + W(RH28F) - W(RH28B) & 
	    + W(RH29F) - W(RH29B) + W(RB12F) & 
	    - W(RB12B) + W(RB16F) - W(RB16B) & 
	    + W(RB17F) - W(RB17B) + W(RB36F) & 
	    - W(RB36B) + W(RB44F) - W(RB44B) & 
	    + W(RB52F) - W(RB52B) + W(RHP02) & 
	    + W(RHP23F) - W(RHP23B) + W(RHP34) & 
	    + W(RHP44) + W(RHP59) + W(RHP67) & 
	    + W(RIC03) + W(RIC24) - W(RIC25) & 
	    + W(RIC41F) - W(RIC41B) + W(RIC59)
      CDOT(SH2) = CDOT(SH2) + W(RIC69) + W(RK017F) - W(RK017B) & 
	    + W(RK021F) - W(RK021B) + W(RK024F) & 
	    - W(RK024B) + W(RP027F) - W(RP027B) & 
	    + W(RK109F) - W(RK109B) + W(RK112F) & 
	    - W(RK112B) + W(RP118F) - W(RP118B) & 
	    + W(RP119F) - W(RP119B) + W(RP120F) & 
	    - W(RP120B) + W(RP121F) - W(RP121B) & 
	    + W(RK123F) - W(RK123B) + W(RK126F) & 
	    - W(RK126B) + W(RK129F) - W(RK129B) & 
	    + W(RK132F) - W(RK132B) + W(RK203F) & 
	    - W(RK203B) + W(RP207F) - W(RP207B) & 
	    + W(RP208F) - W(RP208B) + W(RK210F) & 
	    - W(RK210B) + W(RK213F) - W(RK213B) & 
	    + W(RP305F) - W(RP305B) + W(RK420F) & 
	    - W(RK420B) + W(RK423F) - W(RK423B)
      CDOT(SH2) = CDOT(SH2) + W(RP504F) - W(RP504B) + W(RK507F) & 
	    - W(RK507B) + W(RK602F) - W(RK602B) & 
	    + W(RP800) + W(RP801) + W(RP802) & 
	    + W(RCP02F) - W(RCP02B) + W(RCP25) & 
	    + W(RCP29) + W(RI02F) - W(RI02B) & 
	    + W(RT07F) - W(RT07B) + W(RT24F) & 
	    - W(RT24B) + W(RT32F) - W(RT32B) & 
	    + W(RT37) + W(RT44F) - W(RT44B) & 
	    + W(RT51F) - W(RT51B) + W(RE06F) & 
	    - W(RE06B) + W(RE13F) - W(RE13B) & 
	    + W(RXY02F) - W(RXY02B) + W(RXY17) & 
	    + W(RXY26F) - W(RXY26B) + W(RXY33) & 
	    + W(RXY46) + W(RN08F) - W(RN08B) & 
	    + W(RN24F) - W(RN24B) + W(RN29) & 
	    + W(ROX03F) - W(ROX03B) + W(ROX19F)
      CDOT(SH2) = CDOT(SH2) - W(ROX19B) + W(ROX43F) - W(ROX43B)
      CDOT(SH2O) = W(R3F) - W(R3B) + W(R4F) & 
	    - W(R4B) - W(R7F) + W(R7F) & 
	    - W(R7B) + W(R7B) + W(R9F) & 
	    - W(R9B) + W(R15F) - W(R15B) & 
	    + W(R18F) - W(R18B) + W(R19F) & 
	    - W(R19B) + W(R23F) - W(R23B) & 
	    + W(R25F) - W(R25B) + W(R26F) & 
	    - W(R26B) + W(R35F) - W(R35B) & 
	    - W(R37F) + W(R37F) - W(R37B) & 
	    + W(R37B) - W(RG08F) + W(RG08B) & 
	    + W(RG17F) - W(RG17B) + W(RG36F) & 
	    - W(RG36B) - W(RG37F) + W(RG37B) & 
	    - W(RG38F) + W(RG38F) - W(RG38B) & 
	    + W(RG38B) - W(RG39) + W(RG47F) & 
	    - W(RG47B) + W(RG55F) - W(RG55B)
      CDOT(SH2O) = CDOT(SH2O) + W(RG57F) - W(RG57B) + W(RG79F) & 
	    - W(RG79B) + W(RG81F) - W(RG81B) & 
	    + W(RG86F) - W(RG86B) + W(RG88F) & 
	    - W(RG88B) + W(RG92F) - W(RG92B) & 
	    + W(RG100F) - W(RG100B) + W(RG101F) & 
	    - W(RG101B) + W(RG119F) - W(RG119B) & 
	    + W(RG127F) - W(RG127B) + W(RG132F) & 
	    - W(RG132B) + W(RG143F) - W(RG143B) & 
	    + W(RG151) + W(RG160F) - W(RG160B) & 
	    + W(RG176F) - W(RG176B) + W(RG181F) & 
	    - W(RG181B) + W(RG187F) - W(RG187B) & 
	    + W(RR023F) - W(RR023B) + W(RR050) & 
	    + W(RR059F) - W(RR059B) + W(RR080F) & 
	    - W(RR080B) + W(RR084F) - W(RR084B) & 
	    + W(RR097) + W(RR104F) - W(RR104B)
      CDOT(SH2O) = CDOT(SH2O) + W(RR146F) - W(RR146B) + W(RR151F) & 
	    - W(RR151B) + W(RR155F) - W(RR155B) & 
	    + W(RH11F) - W(RH11B) + W(RH21F) & 
	    - W(RH21B) + W(RH22F) - W(RH22B) & 
	    + W(RH30F) - W(RH30B) + W(RH31F) & 
	    - W(RH31B) + W(RB20F) - W(RB20B) & 
	    + W(RB21F) - W(RB21B) + W(RB37F) & 
	    - W(RB37B) + W(RB46F) - W(RB46B) & 
	    + W(RHP04) + W(RHP24) + W(RHP36) & 
	    + W(RHP45) + W(RHP60) + W(RIC05) & 
	    + W(RIC26) - W(RIC27) + W(RIC43) & 
	    + W(RIC61) + W(RIC71) + W(RP018F) & 
	    - W(RP018B) + W(RP022F) - W(RP022B) & 
	    + W(RP025F) - W(RP025B) + W(RP028F) & 
	    - W(RP028B) + W(RK110F) - W(RK110B)
      CDOT(SH2O) = CDOT(SH2O) + W(RK113F) - W(RK113B) + W(RP124F) & 
	    - W(RP124B) + W(RP127F) - W(RP127B) & 
	    + W(RP130F) - W(RP130B) + W(RP133F) & 
	    - W(RP133B) + W(RK204F) - W(RK204B) & 
	    + W(RP211F) - W(RP211B) + W(RP214F) & 
	    - W(RP214B) + W(RP306F) - W(RP306B) & 
	    + W(RK421F) - W(RK421B) + W(RK424F) & 
	    - W(RK424B) + W(RK508F) - W(RK508B) & 
	    + W(RK603F) - W(RK603B) + W(RCP06F) & 
	    - W(RCP06B) + W(RCP27) + W(RI06F) & 
	    - W(RI06B) + W(RI21) + W(RT08F) & 
	    - W(RT08B) + W(RT26F) - W(RT26B) & 
	    + W(RT34F) - W(RT34B) + W(RT39) & 
	    + W(RT46F) - W(RT46B) + W(RT53F) & 
	    - W(RT53B) + W(RE08F) - W(RE08B)
      CDOT(SH2O) = CDOT(SH2O) + W(RE14F) - W(RE14B) + W(RXY04F) & 
	    - W(RXY04B) + W(RXY28F) - W(RXY28B) & 
	    + W(RXY35) + W(RN10F) - W(RN10B) & 
	    + W(RN26F) - W(RN26B) + W(RN31) & 
	    + W(ROX04F) - W(ROX04B) + W(ROX20F) & 
	    - W(ROX20B) + W(ROX44F) - W(ROX44B)
      CDOT(SCO2) = - W(R8F) + W(R8F) - W(R8B) & 
	    + W(R8B) + W(R27F) - W(R27B) & 
	    + W(R28F) - W(R28B) + W(R29F) & 
	    - W(R29B) + W(R30F) - W(R30B) & 
	    + W(R31F) - W(R31B) + W(R34F) & 
	    - W(R34B) - W(RG11F) + W(RG11B) & 
	    + W(RG19) - W(RG41F) + W(RG41F) & 
	    - W(RG41B) + W(RG41B) - W(RG42F) & 
	    + W(RG42B) + W(RG126F) - W(RG126B) & 
	    + W(RR002F) - W(RR002B) + W(RR033) & 
	    + W(RCP33F) - W(RCP33B) + W(RT59)
      CDOT(SHO2) = W(R12F) - W(R12B) + W(R13F) & 
	    - W(R13B) - W(R15F) + W(R15B) & 
	    - W(R16F) + W(R16B) - W(R17F) & 
	    + W(R17B) - W(R18F) + W(R18B) & 
	    - W(R19F) + W(R19B) - 2 * W(R20F) & 
	    + 2 * W(R20B) - 2 * W(R21F) + 2 * W(R21B) & 
	    + W(R22F) - W(R22B) + W(R24F) & 
	    - W(R24B) + W(R25F) - W(R25B) & 
	    + W(R26F) - W(R26B) - W(R31F) & 
	    + W(R31B) + W(R38F) - W(R38B) & 
	    - W(RG22F) + W(RG22B) + W(RG48F) & 
	    - W(RG48B) - W(RG49F) + W(RG49B) & 
	    - W(RG63) - W(RG65F) + W(RG65B) & 
	    - W(RG66F) + W(RG66B) + W(RG67F) & 
	    - W(RG67B) + W(RG82F) - W(RG82B)
      CDOT(SHO2) = CDOT(SHO2) + W(RG89F) - W(RG89B) + W(RG133F) & 
	    - W(RG133B) + W(RG150) - W(RG152) & 
	    + W(RG170F) - W(RG170B) + W(RG171F) & 
	    - W(RG171B) + W(RG184F) - W(RG184B) & 
	    - W(RG189F) + W(RG189B) + W(RR010F) & 
	    - W(RR010B) + W(RR032F) - W(RR032B) & 
	    - W(RR035F) + W(RR035B) - W(RR036F) & 
	    + W(RR036B) - W(RR037F) + W(RR037B) & 
	    - W(RR038F) + W(RR038B) + W(RR049) & 
	    - W(RR051) - W(RR063F) + W(RR063B) & 
	    - W(RR064F) + W(RR064B) - W(RR065F) & 
	    + W(RR065B) - W(RR082F) + W(RR082B) & 
	    - W(RR086F) + W(RR086B) - W(RR098) & 
	    + W(RR109F) - W(RR109B) - W(RR116F) & 
	    + W(RR116B) - W(RR117F) + W(RR117B)
      CDOT(SHO2) = CDOT(SHO2) - W(RR122F) + W(RR122B) - W(RR128F) & 
	    + W(RR128B) + W(RR131F) - W(RR131B) & 
	    + W(RR132F) - W(RR132B) + W(RR137F) & 
	    - W(RR137B) + W(RR138) - W(RR147F) & 
	    + W(RR147B) + W(RH23F) - W(RH23B) & 
	    + W(RH24F) - W(RH24B) + W(RB39F) & 
	    - W(RB39B) - W(RB40F) + W(RB40B) & 
	    + W(RB41F) - W(RB41B) - W(RB48F) & 
	    + W(RB48B) + W(RB49F) - W(RB49B) & 
	    + W(RHP05) - W(RHP06) - W(RHP09) & 
	    - W(RHP16) + W(RHP46) - W(RHP47) & 
	    - W(RHP50) + W(RHP61) - W(RHP63) & 
	    + W(RHP68) - W(RHP72) - W(RHP76) & 
	    + W(RIC06) - W(RIC08) - W(RIC13) & 
	    - W(RIC28) - W(RIC35) - W(RIC45)
      CDOT(SHO2) = CDOT(SHO2) + W(RIC51) - W(RIC57) - W(RIC62) & 
	    + W(RIC68) - W(RIC72) + W(RCP07F) & 
	    - W(RCP07B) - W(RCP08F) + W(RCP08B) & 
	    - W(RCP21F) + W(RCP21B) + W(RCP28) & 
	    + W(RCP30) + W(RI07F) - W(RI07B) & 
	    - W(RI08F) + W(RI08B) - W(RI21) & 
	    + W(RT06F) - W(RT06B) - W(RT15F) & 
	    + W(RT15B) - W(RT19F) + W(RT19B) & 
	    + W(RT35F) - W(RT35B) + W(RT40) & 
	    - W(RT41) - W(RT57F) + W(RT57B) & 
	    - W(RE09F) + W(RE09B) + W(RE16F) & 
	    - W(RE16B) - W(RE19) + W(RE33) & 
	    + W(RST12F) - W(RST12B) + W(RXY05F) & 
	    - W(RXY05B) - W(RXY06F) + W(RXY06B) & 
	    - W(RXY201) - W(RXY202) - W(RXY203)
      CDOT(SHO2) = CDOT(SHO2) + W(RXY29F) - W(RXY29B) - W(RXY30F) & 
	    + W(RXY30B) + W(RXY36) - W(RXY37) & 
	    - W(RXY48) - W(RXY55) + W(RXY55) & 
	    + W(RN11F) - W(RN11B) - W(RN13F) & 
	    + W(RN13B) - W(RN18F) + W(RN18B) & 
	    + W(RN27F) - W(RN27B) + W(RN32) & 
	    - W(RN33) + W(ROX06F) - W(ROX06B) & 
	    - W(ROX15F) + W(ROX15B) + W(ROX22F) & 
	    - W(ROX22B)
      CDOT(SH2O2) = W(R14F) - W(R14B) + W(R20F) & 
	    - W(R20B) + W(R21F) - W(R21B) & 
	    - W(R22F) + W(R22B) - W(R23F) & 
	    + W(R23B) - W(R24F) + W(R24B) & 
	    - W(R25F) + W(R25B) - W(R26F) & 
	    + W(R26B) + W(RG49F) - W(RG49B) & 
	    - W(RG67F) + W(RG67B) + W(RG152) & 
	    + W(RG189F) - W(RG189B) - W(RR010F) & 
	    + W(RR010B) + W(RR036F) - W(RR036B) & 
	    + W(RR038F) - W(RR038B) + W(RR051) & 
	    + W(RR082F) - W(RR082B) + W(RR086F) & 
	    - W(RR086B) + W(RR098) + W(RR147F) & 
	    - W(RR147B) - W(RB39F) + W(RB39B) & 
	    - W(RB49F) + W(RB49B) + W(RHP06) & 
	    + W(RHP47) + W(RHP63) + W(RIC08)
      CDOT(SH2O2) = CDOT(SH2O2) + W(RIC45) + W(RIC62) + W(RIC72) & 
	    + W(RCP08F) - W(RCP08B) + W(RI08F) & 
	    - W(RI08B) + W(RT15F) - W(RT15B) & 
	    + W(RT41) + W(RE09F) - W(RE09B) & 
	    + W(RXY06F) - W(RXY06B) + W(RXY30F) & 
	    - W(RXY30B) + W(RXY37) + W(RN13F) & 
	    - W(RN13B) + W(RN33)
      CDOT(SCO) = - W(R27F) + W(R27B) - W(R28F) & 
	    + W(R28B) - W(R29F) + W(R29B) & 
	    - W(R30F) + W(R30B) - W(R31F) & 
	    + W(R31B) + W(R32F) - W(R32B) & 
	    + W(R33F) - W(R33B) + W(R35F) & 
	    - W(R35B) + W(R36F) - W(R36B) & 
	    + W(R37F) - W(R37B) + W(R38F) & 
	    - W(R38B) + W(RG01F) - W(RG01B) & 
	    + W(RG02F) - W(RG02B) + W(RG04F) & 
	    - W(RG04B) - W(RG10F) + W(RG10B) & 
	    + W(RG11F) - W(RG11B) - W(RG12F) & 
	    + W(RG12B) + W(RG21) - W(RG24F) & 
	    + W(RG24B) + W(RG31F) - W(RG31B) & 
	    + W(RG35F) - W(RG35B) + W(RG36F) & 
	    - W(RG36B) - W(RG40F) + W(RG40F)
      CDOT(SCO) = CDOT(SCO) - W(RG40B) + W(RG40B) + W(RG42F) & 
	    - W(RG42B) + W(RG53) + W(RG70F) & 
	    - W(RG70B) + W(RG105F) - W(RG105B) & 
	    + W(RG107F) - W(RG107B) + W(RG109F) & 
	    - W(RG109B) + 2 * W(RG110F) - 2 * W(RG110B) & 
	    + 2 * W(RG111F) - 2 * W(RG111B) + W(RG112F) & 
	    - W(RG112B) + W(RG113F) - W(RG113B) & 
	    + 2 * W(RG114F) - 2 * W(RG114B) + W(RG117F) & 
	    - W(RG117B) + W(RG122F) - W(RG122B) & 
	    + W(RG124F) - W(RG124B) + W(RG137F) & 
	    - W(RG137B) + W(RG139) + W(RG148) & 
	    + W(RG149) + W(RG150) + W(RG151) & 
	    + W(RG152) + W(RG153) + W(RR009F) & 
	    - W(RR009B) + W(RR011F) - W(RR011B) & 
	    + W(RR018F) - W(RR018B) + 2 * W(RR019F)
      CDOT(SCO) = CDOT(SCO) - 2 * W(RR019B) + 2 * W(RR020F) - 2 * W(RR020B) & 
	    + 2 * W(RR021F) - 2 * W(RR021B) + W(RR022F) & 
	    - W(RR022B) + W(RR025F) - W(RR025B) & 
	    + W(RR026F) - W(RR026B) + W(RR028F) & 
	    - W(RR028B) + W(RR034F) - W(RR034B) & 
	    + W(RR041F) - W(RR041B) + W(RR045F) & 
	    - W(RR045B) + W(RR047) + W(RR048) & 
	    + W(RR049) + W(RR050) + W(RR051) & 
	    + W(RR052) + W(RR058F) - W(RR058B) & 
	    + W(RR063F) - W(RR063B) + W(RR066F) & 
	    - W(RR066B) + W(RR067F) - W(RR067B) & 
	    + W(RR068F) - W(RR068B) + W(RR089F) & 
	    - W(RR089B) + W(RR094F) - W(RR094B) & 
	    + W(RR095) + W(RR096) + W(RR097) & 
	    + W(RR098) + W(RR099) + W(RR115F)
      CDOT(SCO) = CDOT(SCO) - W(RR115B) + W(RR123F) - W(RR123B) & 
	    + W(RR129F) - W(RR129B) + W(RR135) & 
	    + 2 * W(RH01F) - 2 * W(RH01B) + W(RH08F) & 
	    - W(RH08B) + W(RH12F) - W(RH12B) & 
	    + W(RH27F) - W(RH27B) + W(RH34F) & 
	    - W(RH34B) + W(RB28) + W(RB38F) & 
	    - W(RB38B) + W(RB99F) - W(RB99B) & 
	    + W(RB42) + W(RB47F) - W(RB47B) & 
	    + W(RHP52) + W(RHP58) + W(RHP59) & 
	    + W(RHP60) + W(RHP61) + W(RHP62) & 
	    + W(RHP63) + W(RIC59) + W(RIC60) & 
	    + W(RIC61) + W(RIC62) + W(RIC63) & 
	    + W(RCP13) + W(RCP24) + W(RCP31) & 
	    + W(RI31) + W(RI32) + W(RT36) & 
	    + W(RT37) + W(RT38) + W(RT39)
      CDOT(SCO) = CDOT(SCO) + W(RT40) + W(RT41) + W(RT42) & 
	    + W(RT47) + 2 * W(RT60) + W(RE37) & 
	    + W(RST11F) - W(RST11B) + W(RST13) & 
	    + W(RXY11) + W(RXY24) + W(RXY25) & 
	    + W(RXY33) + W(RXY34) + W(RXY35) & 
	    + W(RXY36) + W(RXY37) + W(RXY38) & 
	    + W(RXY42) + W(RXY44) + 2 * W(RXY50) & 
	    + W(RXY51) + W(RXY52) + W(RXY53) & 
	    + W(RXY54) + W(RXY55) + W(RXY56) & 
	    + W(RXY58) + W(RN14) + W(RN15) & 
	    + W(RN28) + W(RN29) + W(RN30) & 
	    + W(RN31) + W(RN32) + W(RN33) & 
	    + W(RN34) + W(ROX09F) - W(ROX09B) & 
	    + W(ROX17F) - W(ROX17B) + W(ROX23F) & 
	    - W(ROX23B) + W(ROX26F) - W(ROX26B)
      CDOT(SCO) = CDOT(SCO) + W(ROX27F) - W(ROX27B) + 2 * W(ROX28) & 
	    + W(ROX35) + W(ROX36) + W(ROX41F) & 
	    - W(ROX41B) + W(ROX46F) - W(ROX46B) & 
	    + W(ROX47) + W(ROX48) + 2 * W(ROX50) & 
	    + 2 * W(ROX51) + 2 * W(ROX52) + 2 * W(ROX53) & 
	    + 2 * W(ROX54) + W(ROX60) + W(ROX61)
      CDOT(SHCO) = - W(R32F) + W(R32B) - W(R33F) & 
	    + W(R33B) - W(R34F) + W(R34B) & 
	    - W(R35F) + W(R35B) - W(R36F) & 
	    + W(R36B) - W(R37F) + W(R37B) & 
	    - W(R38F) + W(R38B) + W(RG05F) & 
	    - W(RG05B) + W(RG09F) - W(RG09B) & 
	    + W(RG11F) - W(RG11B) - W(RG13F) & 
	    + W(RG13B) + W(RG15F) - W(RG15B) & 
	    + W(RG32F) - W(RG32B) + W(RG45F) & 
	    - W(RG45B) + W(RG46F) - W(RG46B) & 
	    + W(RG47F) - W(RG47B) + W(RG48F) & 
	    - W(RG48B) + W(RG49F) - W(RG49B) & 
	    + W(RG64) - W(RG70F) + W(RG70B) & 
	    + W(RG71F) - W(RG71B) + W(RG107F) & 
	    - W(RG107B) + W(RG135F) - W(RG135B)
      CDOT(SHCO) = CDOT(SHCO) + W(RG138F) - W(RG138B) + 2 * W(RG140) & 
	    + W(RG141F) - W(RG141B) + W(RG144F) & 
	    - W(RG144B) - W(RG145F) + W(RG145B) & 
	    + W(RG159F) - W(RG159B) + 2 * W(RR003F) & 
	    - 2 * W(RR003B) - W(RR011F) + W(RR011B) & 
	    - W(RR012F) + W(RR012B) + 2 * W(RR024F) & 
	    - 2 * W(RR024B) + W(RR030F) - W(RR030B) & 
	    - W(RR034F) + W(RR034B) + W(RR040F) & 
	    - W(RR040B) - W(RR046F) + W(RR046B) & 
	    + W(RR200F) - W(RR200B) + W(RR062F) & 
	    - W(RR062B) + W(RR066F) - W(RR066B) & 
	    - W(RR067F) + W(RR067B) - W(RR068F) & 
	    + W(RR068B) - W(RR115F) + W(RR115B) & 
	    - W(RR123F) + W(RR123B) + W(RR126F) & 
	    - W(RR126B) + W(RR127) + W(RR128F)
      CDOT(SHCO) = CDOT(SHCO) - W(RR128B) - W(RR129F) + W(RR129B) & 
	    + W(RR136F) - W(RR136B) + W(RR143F) & 
	    - W(RR143B) + W(RH27F) - W(RH27B) & 
	    + W(RH36F) - W(RH36B) + W(RB30F) & 
	    - W(RB30B) - W(RB38F) + W(RB38B) & 
	    + W(RB43F) - W(RB43B) - W(RB47F) & 
	    + W(RB47B) + W(RHP19) + W(RHP26) & 
	    + W(RHP54) + W(RHP56) + W(RHP79) & 
	    + W(RIC16) + W(RIC29) + W(RIC47) & 
	    + W(RCP14) + W(RCP15) + W(RI15) & 
	    + W(RT28F) - W(RT28B) + W(RT31F) & 
	    - W(RT31B) + W(RE31) + W(RE34) & 
	    + W(RST02F) - W(RST02B) + W(RST14F) & 
	    - W(RST14B) + W(RXY16F) - W(RXY16B) & 
	    + W(RXY203) + W(RXY39F) - W(RXY39B)
      CDOT(SHCO) = CDOT(SHCO) + W(RXY41F) - W(RXY41B) + W(RXY45F) & 
	    - W(RXY45B) + W(RXY57F) - W(RXY57B) & 
	    + W(RXY58)
      CDOT(SC) = - W(RG01F) + W(RG01B) - W(RG02F) & 
	    + W(RG02B) + W(RG03F) - W(RG03B) & 
	    - W(RG23F) + W(RG23B) - W(RG68F) & 
	    + W(RG68B)
      CDOT(SCH) = - W(RG03F) + W(RG03B) - W(RG04F) & 
	    + W(RG04B) - W(RG05F) + W(RG05B) & 
	    - W(RG06F) + W(RG06B) - W(RG07F) & 
	    + W(RG07B) - W(RG08F) + W(RG08B) & 
	    - W(RG09F) + W(RG09B) - W(RG10F) & 
	    + W(RG10B) - W(RG11F) + W(RG11B) & 
	    + W(RG17F) - W(RG17B) - W(RG25F) & 
	    + W(RG25B) + W(RG30F) - W(RG30B) & 
	    - W(RG50F) + W(RG50B) - W(RG69F) & 
	    + W(RG69B) - W(RG93F) + W(RG93B) & 
	    + W(RG105F) - W(RG105B) - W(RG112F) & 
	    + W(RG112B) + W(RR018F) - W(RR018B) & 
	    - W(RR042F) + W(RR042B) - W(RR069F) & 
	    + W(RR069B)
      CDOT(STXCH2) = W(RG06F) - W(RG06B) - W(RG14F) & 
	    + W(RG14B) - W(RG15F) + W(RG15B) & 
	    - W(RG16F) + W(RG16B) - W(RG17F) & 
	    + W(RG17B) - W(RG18F) + W(RG18B) & 
	    - W(RG19) - W(RG20F) + W(RG20B) & 
	    - W(RG21) - W(RG22F) + W(RG22B) & 
	    - W(RG23F) + W(RG23B) - W(RG24F) & 
	    + W(RG24B) - W(RG25F) + W(RG25B) & 
	    - 2 * W(RG26F) + 2 * W(RG26B) - 2 * W(RG27) & 
	    + W(RG28F) - W(RG28B) + W(RG29F) & 
	    - W(RG29B) + W(RG38F) - W(RG38B) & 
	    + W(RG40F) - W(RG40B) + W(RG41F) & 
	    - W(RG41B) + W(RG55F) - W(RG55B) & 
	    - W(RG72F) + W(RG72B) - W(RG94F) & 
	    + W(RG94B) - W(RG113F) + W(RG113B)
      CDOT(STXCH2) = CDOT(STXCH2) + W(RG117F) - W(RG117B) + W(RG126F) & 
	    - W(RG126B) + W(RG158F) - W(RG158B) & 
	    + W(RR002F) - W(RR002B) - W(RR026F) & 
	    + W(RR026B) - W(RR027F) + W(RR027B) & 
	    - W(RR043F) + W(RR043B) - W(RR070F) & 
	    + W(RR070B) + W(RR087F) - W(RR087B) & 
	    + W(RHP54)
      CDOT(SCH3) = W(RG07F) - W(RG07B) + W(RG14F) & 
	    - W(RG14B) + W(RG18F) - W(RG18B) & 
	    + W(RG34F) - W(RG34B) - W(RG51F) & 
	    + W(RG51B) - W(RG52F) + W(RG52B) & 
	    - W(RG53) - W(RG54F) + W(RG54B) & 
	    - W(RG55F) + W(RG55B) - W(RG56) & 
	    - W(RG57F) + W(RG57B) - W(RG58F) & 
	    + W(RG58B) - W(RG59F) + W(RG59B) & 
	    - W(RG60F) + W(RG60B) - W(RG61F) & 
	    + W(RG61B) - W(RG65F) + W(RG65B) & 
	    - W(RG66F) + W(RG66B) - W(RG67F) & 
	    + W(RG67B) - W(RG68F) + W(RG68B) & 
	    - W(RG69F) + W(RG69B) - W(RG70F) & 
	    + W(RG70B) - W(RG71F) + W(RG71B) & 
	    - W(RG72F) + W(RG72B) - W(RG73F)
      CDOT(SCH3) = CDOT(SCH3) + W(RG73B) - 2 * W(RG74F) + 2 * W(RG74B) & 
	    + W(RG78F) - W(RG78B) + W(RG85F) & 
	    - W(RG85B) + W(RG90F) - W(RG90B) & 
	    + W(RG91F) - W(RG91B) + W(RG92F) & 
	    - W(RG92B) + 2 * W(RG94F) - 2 * W(RG94B) & 
	    + 2 * W(RG95F) - 2 * W(RG95B) - W(RG102F) & 
	    + W(RG102B) - W(RG103F) + W(RG103B) & 
	    + W(RG122F) - W(RG122B) + W(RG124F) & 
	    - W(RG124B) + W(RG137F) - W(RG137B) & 
	    + W(RG141F) - W(RG141B) - W(RG145F) & 
	    + W(RG145B) + W(RG148) + W(RG149) & 
	    + W(RG150) + W(RG151) + W(RG152) & 
	    - W(RG153) + W(RG153) + W(RG159F) & 
	    - W(RG159B) - W(RG162F) + W(RG162B) & 
	    - W(RG163F) + W(RG163B) - W(RG166F)
      CDOT(SCH3) = CDOT(SCH3) + W(RG166B) + W(RG168F) - W(RG168B) & 
	    + W(RG172F) - W(RG172B) + 2 * W(RG173F) & 
	    - 2 * W(RG173B) + W(RG177F) - W(RG177B) & 
	    - W(RG178F) + W(RG178B) - W(RG182F) & 
	    + W(RG182B) - W(RG188F) + W(RG188B) & 
	    + W(RR005F) - W(RR005B) + W(RR006F) & 
	    - W(RR006B) - W(RR007F) + W(RR007B) & 
	    - W(RR013F) + W(RR013B) + W(RR014F) & 
	    - W(RR014B) - W(RR015F) + W(RR015B) & 
	    - W(RR017F) + W(RR017B) - W(RR022F) & 
	    + W(RR022B) + W(RR027F) - W(RR027B) & 
	    - W(RR028F) + W(RR028B) - W(RR029F) & 
	    + W(RR029B) - W(RR030F) + W(RR030B) & 
	    + W(RR033) - W(RR044F) + W(RR044B) & 
	    - W(RR052) + W(RR066F) - W(RR066B)
      CDOT(SCH3) = CDOT(SCH3) - W(RR081F) + W(RR081B) - W(RR085F) & 
	    + W(RR085B) + W(RR088F) - W(RR088B) & 
	    + W(RR092F) - W(RR092B) + W(RR093F) & 
	    - W(RR093B) - W(RR099) - W(RR105F) & 
	    + W(RR105B) + W(RR120F) - W(RR120B) & 
	    + W(RR121) + W(RR122F) - W(RR122B) & 
	    - W(RR124F) + W(RR124B) - W(RR130F) & 
	    + W(RR130B) + W(RR133) + W(RR135) & 
	    + W(RR141F) - W(RR141B) + W(RR142F) & 
	    - W(RR142B) - W(RR148F) + W(RR148B) & 
	    - W(RR152F) + W(RR152B) - W(RR156F) & 
	    + W(RR156B) - W(RH32F) + W(RH32B) & 
	    - W(RH33F) + W(RH33B) - W(RB08F) & 
	    + W(RB08B) + W(RB09F) - W(RB09B) & 
	    - W(RB14F) + W(RB14B) - W(RB15F)
      CDOT(SCH3) = CDOT(SCH3) + W(RB15B) - W(RB22F) + W(RB22B) & 
	    - W(RB23F) + W(RB23B) + W(RB45F) & 
	    - W(RB45B) + W(RHP38) + W(RHP42) & 
	    + W(RHP43) - W(RHP48) + W(RHP52) & 
	    + W(RHP58) + W(RHP56) - W(RHP62) & 
	    - W(RHP69) + W(RHP75) + W(RIC00) & 
	    + W(RIC01) - W(RIC07) + W(RIC11) & 
	    + W(RIC12) + W(RIC15) + W(RIC31) & 
	    + W(RIC32) + W(RIC37) + W(RIC40) & 
	    - W(RIC44) + 2 * W(RIC48) + W(RIC52) & 
	    + W(RIC55) - W(RIC63) + W(RIC64) & 
	    - W(RIC67) + W(RIC69) + W(RIC70) & 
	    + W(RIC71) + W(RIC72) - W(RCP09F) & 
	    + W(RCP09B) - W(RCP18) + W(RCP34) & 
	    - W(RI09F) + W(RI09B) - W(RI18)
      CDOT(SCH3) = CDOT(SCH3) + W(RT01F) - W(RT01B) + W(RT03F) & 
	    - W(RT03B) + W(RT04F) - W(RT04B) & 
	    + W(RT09F) - W(RT09B) - W(RT14F) & 
	    + W(RT14B) - W(RT27F) + W(RT27B) & 
	    - W(RT42) - W(RT54F) + W(RT54B) & 
	    - W(RE02F) + W(RE02B) - W(RE10F) & 
	    + W(RE10B) - W(RE15F) + W(RE15B) & 
	    + W(RE18F) - W(RE18B) - W(RST03F) & 
	    + W(RST03B) + W(RST04F) - W(RST04B) & 
	    + W(RXY01F) - W(RXY01B) - W(RXY07F) & 
	    + W(RXY07B) + W(RXY09F) - W(RXY09B) & 
	    + W(RXY10F) - W(RXY10B) + W(RXY13F) & 
	    - W(RXY13B) + W(RXY24) - W(RXY31F) & 
	    + W(RXY31B) - W(RXY38) + W(RXY40F) & 
	    - W(RXY40B) + W(RXY42) - W(RXY56)
      CDOT(SCH3) = CDOT(SCH3) + W(RXY56) + W(RN01F) - W(RN01B) & 
	    + W(RN02F) - W(RN02B) + W(RN05F) & 
	    - W(RN05B) + W(RN06F) - W(RN06B) & 
	    - W(RN12F) + W(RN12B) + W(RN15) & 
	    - W(RN34) + W(ROX16F) - W(ROX16B) & 
	    - W(ROX21F) + W(ROX21B) - W(ROX45F) & 
	    + W(ROX45B) + W(ROX60) + W(ROX61)
      CDOT(SCH2O) = W(RG08F) - W(RG08B) + W(RG12F) & 
	    - W(RG12B) + W(RG13F) - W(RG13B) & 
	    + W(RG16F) - W(RG16B) + W(RG20F) & 
	    - W(RG20B) + W(RG22F) - W(RG22B) & 
	    + W(RG33F) - W(RG33B) + W(RG39) & 
	    + W(RG42F) - W(RG42B) - W(RG43F) & 
	    + W(RG43B) - W(RG44F) + W(RG44B) & 
	    - W(RG45F) + W(RG45B) - W(RG46F) & 
	    + W(RG46B) - W(RG47F) + W(RG47B) & 
	    - W(RG48F) + W(RG48B) - W(RG49F) & 
	    + W(RG49B) - W(RG50F) + W(RG50B) & 
	    + W(RG52F) - W(RG52B) + W(RG56) & 
	    + W(RG59F) - W(RG59B) - W(RG64) & 
	    - W(RG71F) + W(RG71B) + W(RG77F) & 
	    - W(RG77B) + W(RG80F) - W(RG80B)
      CDOT(SCH2O) = CDOT(SCH2O) + W(RG81F) - W(RG81B) + W(RG82F) & 
	    - W(RG82B) + W(RG84F) - W(RG84B) & 
	    + W(RG87F) - W(RG87B) + W(RG88F) & 
	    - W(RG88B) + W(RG89F) - W(RG89B) & 
	    + W(RG135F) - W(RG135B) + W(RG138F) & 
	    - W(RG138B) + W(RG139) + W(RG158F) & 
	    - W(RG158B) + W(RG168F) - W(RG168B) & 
	    + W(RG179F) - W(RG179B) + W(RR060F) & 
	    - W(RR060B) + W(RR111) + W(RR112F) & 
	    - W(RR112B) + W(RR135) + W(RR140) & 
	    + W(RH35F) - W(RH35B) + W(RB29F) & 
	    - W(RB29B) + W(RB31F) - W(RB31B) & 
	    + W(RB51F) - W(RB51B) + W(RHP53) & 
	    + W(RHP55) + W(RHP78) + W(RIC17) & 
	    + W(RIC30) + W(RIC49) + W(RIC55)
      CDOT(SCH2O) = CDOT(SCH2O) + W(RIC56) + W(RT22F) - W(RT22B) & 
	    - W(RT28F) + W(RT28B) + W(RT30F) & 
	    - W(RT30B) + W(RE17F) - W(RE17B) & 
	    + W(RE19) + W(RE37) + W(RST05F) & 
	    - W(RST05B) + W(RXY15F) - W(RXY15B) & 
	    + W(RXY19F) - W(RXY19B) + W(RXY202) & 
	    + W(RXY44) + W(RN21F) - W(RN21B) & 
	    + W(RN23F) - W(RN23B)
      CDOT(SHCCO) = W(RG10F) - W(RG10B) + W(RG106F) & 
	    - W(RG106B) - W(RG109F) + W(RG109B) & 
	    - W(RG110F) + W(RG110B) - W(RG111F) & 
	    + W(RG111B) - W(RG112F) + W(RG112B) & 
	    - W(RG113F) + W(RG113B) - 2 * W(RG114F) & 
	    + 2 * W(RG114B) + W(RG116F) - W(RG116B) & 
	    + W(RG123F) - W(RG123B) + W(RG125F) & 
	    - W(RG125B) + W(RG127F) - W(RG127B) & 
	    - W(RR009F) + W(RR009B) - W(RR022F) & 
	    + W(RR022B) - W(RR023F) + W(RR023B) & 
	    - W(RR024F) + W(RR024B) + W(RR027F) & 
	    - W(RR027B) + W(RR029F) - W(RR029B) & 
	    + W(RR041F) - W(RR041B) - W(RR045F) & 
	    + W(RR045B) + W(RR088F) - W(RR088B) & 
	    + 2 * W(RH07F) - 2 * W(RH07B) + W(RH26F)
      CDOT(SHCCO) = CDOT(SHCCO) - W(RH26B) + W(RST00F) - W(RST00B) & 
	    + W(ROX62) + W(ROX63) + W(ROX64) & 
	    + W(ROX65) + W(ROX66)
      CDOT(SC2H) = W(RG23F) - W(RG23B) - W(RG104F) & 
	    + W(RG104B) - W(RG105F) + W(RG105B) & 
	    - W(RG106F) + W(RG106B) - W(RG107F) & 
	    + W(RG107B) - W(RG108F) + W(RG108B) & 
	    - W(RG118F) + W(RG118B) + W(RG119F) & 
	    - W(RG119B) - W(RR008F) + W(RR008B) & 
	    - W(RR017F) + W(RR017B) - W(RR031F) & 
	    + W(RR031B) - W(RR046F) + W(RR046B) & 
	    + W(RR047) + W(RR048) + W(RR049) & 
	    + W(RR050) + W(RR051) + W(RR052) & 
	    - W(RR090F) + W(RR090B) - W(RR091F) & 
	    + W(RR091B) + W(RH01F) - W(RH01B) & 
	    + W(RH25F) - W(RH25B) - W(RH37F) & 
	    + W(RH37B) - W(RH39F) + W(RH39B)
      CDOT(SCH2CO) = W(RG24F) - W(RG24B) + W(RG50F) & 
	    - W(RG50B) + W(RG121F) - W(RG121B) & 
	    - W(RG123F) + W(RG123B) - W(RG124F) & 
	    + W(RG124B) - W(RG125F) + W(RG125B) & 
	    - W(RG126F) + W(RG126B) - W(RG127F) & 
	    + W(RG127B) + W(RG128F) - W(RG128B) & 
	    + W(RG136F) - W(RG136B) + W(RG142F) & 
	    - W(RG142B) + W(RG143F) - W(RG143B) & 
	    - W(RR025F) + W(RR025B) - W(RR026F) & 
	    + W(RR026B) - W(RR027F) + W(RR027B) & 
	    - W(RR028F) + W(RR028B) - W(RR029F) & 
	    + W(RR029B) + W(RR062F) - W(RR062B) & 
	    + W(RR087F) - W(RR087B) + W(RR093F) & 
	    - W(RR093B) + W(RR120F) - W(RR120B) & 
	    + W(RR121) + W(RR122F) - W(RR122B)
      CDOT(SCH2CO) = CDOT(SCH2CO) + W(RR133) + W(RH25F) - W(RH25B) & 
	    + W(RH26F) - W(RH26B) + W(RB50F) & 
	    - W(RB50B) + W(RHP86) + W(RHP87) & 
	    + W(RIC48) + W(RIC55) + W(RIC69) & 
	    + W(RIC70) + W(RIC71) + W(RIC72) & 
	    + W(ROX28)
      CDOT(SC2H2) = W(RG25F) - W(RG25B) + W(RG26F) & 
	    - W(RG26B) + W(RG27) + W(RG68F) & 
	    - W(RG68B) + W(RG104F) - W(RG104B) & 
	    + W(RG108F) - W(RG108B) + W(RG112F) & 
	    - W(RG112B) + W(RG114F) - W(RG114B) & 
	    - W(RG115F) + W(RG115B) - W(RG116F) & 
	    + W(RG116B) - W(RG117F) + W(RG117B) & 
	    + W(RG118F) - W(RG118B) - W(RG119F) & 
	    + W(RG119B) - W(RG120F) + W(RG120B) & 
	    - W(RG121F) + W(RG121B) - W(RG122F) & 
	    + W(RG122B) + W(RG130F) - W(RG130B) & 
	    + W(RG132F) - W(RG132B) + W(RG133F) & 
	    - W(RG133B) - W(RR001F) + W(RR001B) & 
	    - W(RR004F) + W(RR004B) + W(RR005F) & 
	    - W(RR005B) + W(RR006F) - W(RR006B)
      CDOT(SC2H2) = CDOT(SC2H2) - W(RR007F) + W(RR007B) - W(RR008F) & 
	    + W(RR008B) - W(RR009F) + W(RR009B) & 
	    + W(RR013F) - W(RR013B) + W(RR040F) & 
	    - W(RR040B) + W(RR200F) - W(RR200B) & 
	    + W(RR060F) - W(RR060B) + W(RR090F) & 
	    - W(RR090B) + W(RR091F) - W(RR091B) & 
	    + W(RR111) + 2 * W(RH17F) - 2 * W(RH17B) & 
	    + 2 * W(RH18F) - 2 * W(RH18B) + W(RH27F) & 
	    - W(RH27B) - W(RH38F) + W(RH38B) & 
	    - W(RB01F) + W(RB01B) - W(RB02F) & 
	    + W(RB02B) + W(RB07F) - W(RB07B) & 
	    - W(RP001F) + W(RP001B) - W(RP002F) & 
	    + W(RP002B) - W(RP005F) + W(RP005B) & 
	    - W(RP006F) + W(RP006B) - W(RP007F) & 
	    + W(RP007B) - W(RK012F) + W(RK012B)
      CDOT(SC2H2) = CDOT(SC2H2) - W(RK100F) + W(RK100B) - W(RK102F) & 
	    + W(RK102B) - W(RP104F) + W(RP104B) & 
	    - W(RK114F) + W(RK114B) - W(RK115F) & 
	    + W(RK115B) - W(RK205F) + W(RK205B) & 
	    - W(RK401F) + W(RK401B) - W(RK403F) & 
	    + W(RK403B) - W(RP405F) + W(RP405B) & 
	    - W(RP406F) + W(RP406B) - W(RP407F) & 
	    + W(RP407B) + W(RP425F) - W(RP425B) & 
	    - W(RP501F) + W(RP501B) - W(RP502F) & 
	    + W(RP502B) - W(RP503F) + W(RP503B) & 
	    - W(RK600F) + W(RK600B) - W(RK700F) & 
	    + W(RK700B) - W(RK701F) + W(RK701B) & 
	    + W(RCP03F) - W(RCP03B) + W(RCP04) & 
	    + W(RCP13) + W(RCP15) - W(RCP16F) & 
	    + W(RCP16B) + 2 * W(RCP31) - W(RI03F)
      CDOT(SC2H2) = CDOT(SC2H2) + W(RI03B) + W(RI25) + W(RI31) & 
	    + W(RT05F) - W(RT05B) + W(RXY12) & 
	    + W(RN07F) - W(RN07B) - W(ROX02F) & 
	    + W(ROX02B) + W(ROX28)
      CDOT(SSXCH2) = - W(RG28F) + W(RG28B) - W(RG29F) & 
	    + W(RG29B) - W(RG30F) + W(RG30B) & 
	    - W(RG31F) + W(RG31B) - W(RG32F) & 
	    + W(RG32B) - W(RG33F) + W(RG33B) & 
	    - W(RG34F) + W(RG34B) - W(RG35F) & 
	    + W(RG35B) - W(RG36F) + W(RG36B) & 
	    - W(RG37F) + W(RG37B) - W(RG38F) & 
	    + W(RG38B) - W(RG39) - W(RG40F) & 
	    + W(RG40B) - W(RG41F) + W(RG41B) & 
	    - W(RG42F) + W(RG42B) + W(RG57F) & 
	    - W(RG57B) - W(RG73F) + W(RG73B) & 
	    + W(RG79F) - W(RG79B) + W(RG86F) & 
	    - W(RG86B) - W(RG95F) + W(RG95B) & 
	    + W(RG109F) - W(RG109B) - W(RG177F) & 
	    + W(RG177B) - W(RR004F) + W(RR004B)
      CDOT(SAR) = - W(RG29F) + W(RG29F) - W(RG29B) & 
	    + W(RG29B)
      CDOT(SCH3OH) = W(RG37F) - W(RG37B) + W(RG54F) & 
	    - W(RG54B) + W(RG75F) - W(RG75B) & 
	    + W(RG83F) - W(RG83B) - W(RG96F) & 
	    + W(RG96B) - W(RG97F) + W(RG97B) & 
	    - W(RG98F) + W(RG98B) - W(RG99F) & 
	    + W(RG99B) - W(RG100F) + W(RG100B) & 
	    - W(RG101F) + W(RG101B) - W(RG102F) & 
	    + W(RG102B) - W(RG103F) + W(RG103B)
      CDOT(SCH2OH) = W(RG43F) - W(RG43B) + W(RG76F) & 
	    - W(RG76B) - W(RG83F) + W(RG83B) & 
	    - W(RG84F) + W(RG84B) - W(RG85F) & 
	    + W(RG85B) - W(RG86F) + W(RG86B) & 
	    - W(RG87F) + W(RG87B) - W(RG88F) & 
	    + W(RG88B) - W(RG89F) + W(RG89B) & 
	    + W(RG96F) - W(RG96B) + W(RG98F) & 
	    - W(RG98B) + W(RG100F) - W(RG100B) & 
	    + W(RG102F) - W(RG102B) + W(RG144F) & 
	    - W(RG144B) + W(RR025F) - W(RR025B)
      CDOT(SCH3O) = W(RG44F) - W(RG44B) + W(RG58F) & 
	    - W(RG58B) + 2 * W(RG61F) - 2 * W(RG61B) & 
	    + 2 * W(RG62) + W(RG63) + W(RG64) & 
	    + W(RG65F) - W(RG65B) - W(RG75F) & 
	    + W(RG75B) - W(RG76F) + W(RG76B) & 
	    - W(RG77F) + W(RG77B) - W(RG78F) & 
	    + W(RG78B) - W(RG79F) + W(RG79B) & 
	    - W(RG80F) + W(RG80B) - W(RG81F) & 
	    + W(RG81B) - W(RG82F) + W(RG82B) & 
	    + W(RG97F) - W(RG97B) + W(RG99F) & 
	    - W(RG99B) + W(RG101F) - W(RG101B) & 
	    + W(RG103F) - W(RG103B) + W(RR118F) & 
	    - W(RR118B) + W(RHP17) + W(RHP73) & 
	    + W(RHP77) + W(RIC09) + W(RIC14) & 
	    + W(RIC36) + W(RIC46) + W(RIC58)
      CDOT(SCH4) = W(RG51F) - W(RG51B) + W(RG66F) & 
	    - W(RG66B) + W(RG67F) - W(RG67B) & 
	    + W(RG70F) - W(RG70B) + W(RG71F) & 
	    - W(RG71B) - W(RG90F) + W(RG90B) & 
	    - W(RG91F) + W(RG91B) - W(RG92F) & 
	    + W(RG92B) - W(RG93F) + W(RG93B) & 
	    - W(RG94F) + W(RG94B) - W(RG95F) & 
	    + W(RG95B) + W(RG102F) - W(RG102B) & 
	    + W(RG103F) - W(RG103B) + W(RG153) & 
	    + W(RG162F) - W(RG162B) + W(RG166F) & 
	    - W(RG166B) + W(RG178F) - W(RG178B) & 
	    + W(RG182F) - W(RG182B) + W(RG188F) & 
	    - W(RG188B) + W(RR013F) - W(RR013B) & 
	    + W(RR029F) - W(RR029B) + W(RR052) & 
	    + W(RR081F) - W(RR081B) + W(RR085F)
      CDOT(SCH4) = CDOT(SCH4) - W(RR085B) + W(RR099) + W(RR105F) & 
	    - W(RR105B) + W(RR124F) - W(RR124B) & 
	    + W(RR130F) - W(RR130B) + W(RR148F) & 
	    - W(RR148B) + W(RR152F) - W(RR152B) & 
	    + W(RR156F) - W(RR156B) + W(RH32F) & 
	    - W(RH32B) + W(RH33F) - W(RH33B) & 
	    + W(RB22F) - W(RB22B) + W(RB23F) & 
	    - W(RB23B) + W(RHP48) + W(RHP62) & 
	    + W(RHP69) + W(RIC07) + W(RIC44) & 
	    + W(RIC63) + W(RCP09F) - W(RCP09B) & 
	    + W(RI09F) - W(RI09B) + W(RT14F) & 
	    - W(RT14B) + W(RT27F) - W(RT27B) & 
	    + W(RT42) + W(RT54F) - W(RT54B) & 
	    + W(RE10F) - W(RE10B) + W(RE15F) & 
	    - W(RE15B) + W(RST03F) - W(RST03B)
      CDOT(SCH4) = CDOT(SCH4) + W(RXY07F) - W(RXY07B) + W(RXY31F) & 
	    - W(RXY31B) + W(RXY38) + W(RN12F) & 
	    - W(RN12B) + W(RN34) - W(ROX16F) & 
	    + W(ROX16B) + W(ROX21F) - W(ROX21B) & 
	    + W(ROX45F) - W(ROX45B)
      CDOT(SCH3O2) = W(RG60F) - W(RG60B) - W(RG61F) & 
	    + W(RG61B) - 2 * W(RG62) - W(RG63) & 
	    - W(RG64) - W(RR118F) + W(RR118B) & 
	    - W(RHP17) - W(RHP73) - W(RHP77) & 
	    - W(RIC09) - W(RIC14) - W(RIC36) & 
	    - W(RIC46) - W(RIC58)
      CDOT(SC2H3) = W(RG69F) - W(RG69B) + W(RG113F) & 
	    - W(RG113B) + W(RG115F) - W(RG115B) & 
	    - W(RG129F) + W(RG129B) - W(RG130F) & 
	    + W(RG130B) - W(RG131F) + W(RG131B) & 
	    - W(RG132F) + W(RG132B) - W(RG133F) & 
	    + W(RG133B) - W(RG134F) + W(RG134B) & 
	    - W(RG135F) + W(RG135B) + W(RG156F) & 
	    - W(RG156B) + W(RG160F) - W(RG160B) & 
	    + W(RG162F) - W(RG162B) - W(RR010F) & 
	    + W(RR010B) - W(RR011F) + W(RR011B) & 
	    - W(RR012F) + W(RR012B) - W(RR013F) & 
	    + W(RR013B) + W(RR014F) - W(RR014B) & 
	    - W(RR015F) + W(RR015B) + W(RR032F) & 
	    - W(RR032B) + W(RR063F) - W(RR063B) & 
	    + W(RR095) + W(RR096) + W(RR097)
      CDOT(SC2H3) = CDOT(SC2H3) + W(RR098) + W(RR099) - W(RR106F) & 
	    + W(RR106B) + W(RR140) - W(RB02F) & 
	    + W(RB02B) - 2 * W(RB04F) + 2 * W(RB04B) & 
	    - 2 * W(RB05F) + 2 * W(RB05B) - 2 * W(RB06F) & 
	    + 2 * W(RB06B) - 2 * W(RB07F) + 2 * W(RB07B) & 
	    - W(RB09F) + W(RB09B) - W(RB24F) & 
	    + W(RB24B) - W(RB25F) + W(RB25B) & 
	    - W(RB52F) + W(RB52B) - W(RHP65F) & 
	    + W(RHP65B) + W(RHP74) + W(RHP88) & 
	    - W(RP013F) + W(RP013B) - W(RP014F) & 
	    + W(RP014B) + W(RP015F) - W(RP015B) & 
	    - W(RP105F) + W(RP105B) - W(RP116F) & 
	    + W(RP116B) - W(RP117F) + W(RP117B) & 
	    - W(RP118F) + W(RP118B) - W(RP119F) & 
	    + W(RP119B) - W(RP206F) + W(RP206B)
      CDOT(SC2H3) = CDOT(SC2H3) - W(RP207F) + W(RP207B) - W(RP408F) & 
	    + W(RP408B) - W(RP409F) + W(RP409B) & 
	    - W(RP410F) + W(RP410B) - W(RP411F) & 
	    + W(RP411B) - W(RP504F) + W(RP504B) & 
	    - W(RCP10F) + W(RCP10B) + W(RT60) & 
	    + W(RST06F) - W(RST06B)
      CDOT(SC2H4) = W(RG72F) - W(RG72B) + W(RG73F) & 
	    - W(RG73B) + W(RG93F) - W(RG93B) & 
	    + W(RG129F) - W(RG129B) - W(RG154F) & 
	    + W(RG154B) - W(RG155F) + W(RG155B) & 
	    - W(RG156F) + W(RG156B) - W(RG157F) & 
	    + W(RG157B) - W(RG158F) + W(RG158B) & 
	    - W(RG159F) + W(RG159B) - W(RG160F) & 
	    + W(RG160B) - W(RG161F) + W(RG161B) & 
	    - W(RG162F) + W(RG162B) - W(RG163F) & 
	    + W(RG163B) + W(RG165F) - W(RG165B) & 
	    + W(RG166F) - W(RG166B) + W(RG171F) & 
	    - W(RG171B) + W(RR010F) - W(RR010B) & 
	    + W(RR011F) - W(RR011B) + W(RR022F) & 
	    - W(RR022B) + W(RR026F) - W(RR026B) & 
	    - W(RR031F) + W(RR031B) - W(RR032F)
      CDOT(SC2H4) = CDOT(SC2H4) + W(RR032B) - W(RR033) + W(RR036F) & 
	    - W(RR036B) + W(RR058F) - W(RR058B) & 
	    + W(RR089F) - W(RR089B) + W(RR106F) & 
	    - W(RR106B) + W(RR126F) - W(RR126B) & 
	    + W(RR127) + W(RR128F) - W(RR128B) & 
	    + W(RR141F) - W(RR141B) - W(RB00F) & 
	    + W(RB00B) + W(RB07F) - W(RB07B) & 
	    + W(RB24F) - W(RB24B) + W(RB25F) & 
	    - W(RB25B) + W(RHP10) + W(RHP29) & 
	    + W(RHP37) + W(RHP41) + W(RHP51) & 
	    - W(RHP65F) + W(RHP65B) + W(RHP86) & 
	    + W(RIC64) - W(RIC67) - W(RP015F) & 
	    + W(RP015B) - W(RP106F) + W(RP106B) & 
	    - W(RP120F) + W(RP120B) - W(RP121F) & 
	    + W(RP121B) - W(RP208F) + W(RP208B)
      CDOT(SC2H4) = CDOT(SC2H4) - W(RP412F) + W(RP412B) - W(RP413F) & 
	    + W(RP413B) + W(RCP10F) - W(RCP10B) & 
	    + W(RCP13) + W(RCP15) + W(RI15) & 
	    + W(RE11F) - W(RE11B)
      CDOT(SC2H5) = W(RG74F) - W(RG74B) + W(RG155F) & 
	    - W(RG155B) - W(RG164F) + W(RG164B) & 
	    - W(RG165F) + W(RG165B) - W(RG166F) & 
	    + W(RG166B) - W(RG167F) + W(RG167B) & 
	    - W(RG171F) + W(RG171B) + W(RG172F) & 
	    - W(RG172B) + W(RG174F) - W(RG174B) & 
	    + W(RG175F) - W(RG175B) + W(RG176F) & 
	    - W(RG176B) + W(RG177F) - W(RG177B) & 
	    + W(RG178F) - W(RG178B) + W(RG179F) & 
	    - W(RG179B) + W(RR028F) - W(RR028B) & 
	    + W(RR030F) - W(RR030B) - W(RR034F) & 
	    + W(RR034B) - W(RR035F) + W(RR035B) & 
	    - W(RR036F) + W(RR036B) - W(RR037F) & 
	    + W(RR037B) + W(RR038F) - W(RR038B) & 
	    + W(RR094F) - W(RR094B) - W(RR107F)
      CDOT(SC2H5) = CDOT(SC2H5) + W(RR107B) + W(RR143F) - W(RR143B) & 
	    + W(RHP00) + W(RHP13) + W(RHP19) & 
	    + W(RHP26) + W(RHP32) + W(RHP33) & 
	    + W(RHP41) + W(RHP52) + W(RHP54) & 
	    + W(RHP56) + W(RHP57) - W(RHP70) & 
	    + W(RHP78) + W(RHP79) + W(RHP87) & 
	    - W(RE03F) + W(RE03B) + W(RE04F) & 
	    - W(RE04B) + W(RE05F) - W(RE05B)
      CDOT(SHCCOH) = W(RG120F) - W(RG120B) - W(RG128F) & 
	    + W(RG128B) + W(RR092F) - W(RR092B)
      CDOT(SCH2CHO) = W(RG131F) - W(RG131B) + W(RG134F) & 
	    - W(RG134B) - W(RG136F) + W(RG136B) & 
	    - W(RG137F) + W(RG137B) - W(RG138F) & 
	    + W(RG138B) - W(RG139) - W(RG140) & 
	    - W(RG141F) + W(RG141B) - W(RG142F) & 
	    + W(RG142B) - W(RG143F) + W(RG143B) & 
	    - W(RG144F) + W(RG144B) + W(RG146F) & 
	    - W(RG146B) + W(RG147F) - W(RG147B) & 
	    + W(RG157F) - W(RG157B) - W(RR030F) & 
	    + W(RR030B) + W(RR112F) - W(RR112B) & 
	    + W(RR142F) - W(RR142B) + W(RB50F) & 
	    - W(RB50B)
      CDOT(SCH3CHO) = W(RG145F) - W(RG145B) - W(RG146F) & 
	    + W(RG146B) - W(RG147F) + W(RG147B) & 
	    - W(RG148) - W(RG149) - W(RG150) & 
	    - W(RG151) - W(RG152) - W(RG153) & 
	    + W(RG169F) - W(RG169B) + W(RG170F) & 
	    - W(RG170B) + W(RR136F) - W(RR136B) & 
	    + W(RHP20) + W(RHP25) + W(RHP51) & 
	    + W(RHP57) + W(RHP74) + W(RHP88)
      CDOT(SH2C2) = W(RG154F) - W(RG154B) + W(RR001F) & 
	    - W(RR001B) - W(RR002F) + W(RR002B) & 
	    - W(RR003F) + W(RR003B) - W(RB00F) & 
	    + W(RB00B) - W(RB01F) + W(RB01B) & 
	    + W(RST01)
      CDOT(SC2H5O) = W(RG161F) - W(RG161B) + W(RG167F) & 
	    - W(RG167B) - W(RG168F) + W(RG168B) & 
	    - W(RG169F) + W(RG169B) - W(RG170F) & 
	    + W(RG170B) + W(RR037F) - W(RR037B)
      CDOT(SNXC3H7) = W(RG163F) - W(RG163B) - W(RG179F) & 
	    + W(RG179B) - W(RG180F) + W(RG180B) & 
	    - W(RG181F) + W(RG181B) - W(RG182F) & 
	    + W(RG182B) + W(RG183F) - W(RG183B) & 
	    - W(RG184F) + W(RG184B) + W(RG185F) & 
	    - W(RG185B) + W(RG186F) - W(RG186B) & 
	    + W(RG187F) - W(RG187B) + W(RG188F) & 
	    - W(RG188B) + W(RG189F) - W(RG189B) & 
	    + W(RHP01) + W(RHP12) + W(RHP18) & 
	    + W(RHP21) + W(RHP29) + W(RHP55) & 
	    + W(RHP59) + W(RHP60) + W(RHP61) & 
	    + W(RHP62) + W(RHP63) - W(RHP76) & 
	    - W(RHP77)
      CDOT(SC2H6) = W(RG164F) - W(RG164B) - W(RG173F) & 
	    + W(RG173B) - W(RG174F) + W(RG174B) & 
	    - W(RG175F) + W(RG175B) - W(RG176F) & 
	    + W(RG176B) - W(RG177F) + W(RG177B) & 
	    - W(RG178F) + W(RG178B) + W(RR034F) & 
	    - W(RR034B) + W(RR035F) - W(RR035B) & 
	    - W(RR038F) + W(RR038B) + W(RR107F) & 
	    - W(RR107B) + W(RHP58) + W(RHP70)
      CDOT(SC3H8) = - W(RG172F) + W(RG172B) + W(RG180F) & 
	    - W(RG180B) - W(RG185F) + W(RG185B) & 
	    - W(RG186F) + W(RG186B) - W(RG187F) & 
	    + W(RG187B) - W(RG188F) + W(RG188B) & 
	    - W(RG189F) + W(RG189B)
      CDOT(SC3H6) = W(RG181F) - W(RG181B) + W(RG182F) & 
	    - W(RG182B) - W(RG183F) + W(RG183B) & 
	    + W(RG184F) - W(RG184B) - W(RR014F) & 
	    + W(RR014B) + W(RR016F) - W(RR016B) & 
	    + W(RR108F) - W(RR108B) + W(RR115F) & 
	    - W(RR115B) + W(RR116F) - W(RR116B) & 
	    + W(RR123F) - W(RR123B) + W(RR129F) & 
	    - W(RR129B) - W(RR141F) + W(RR141B) & 
	    - W(RR142F) + W(RR142B) - W(RR143F) & 
	    + W(RR143B) - W(RR144F) + W(RR144B) & 
	    - W(RR145F) + W(RR145B) - W(RR146F) & 
	    + W(RR146B) - W(RR147F) + W(RR147B) & 
	    - W(RR148F) + W(RR148B) - W(RR149F) & 
	    + W(RR149B) - W(RR150F) + W(RR150B) & 
	    - W(RR151F) + W(RR151B) - W(RR152F)
      CDOT(SC3H6) = CDOT(SC3H6) + W(RR152B) - W(RR153F) + W(RR153B) & 
	    - W(RR154F) + W(RR154B) - W(RR155F) & 
	    + W(RR155B) - W(RR156F) + W(RR156B) & 
	    - W(RB09F) + W(RB09B) + W(RHP11) & 
	    + W(RHP28) + W(RHP32) + W(RHP42) & 
	    + W(RHP53) + W(RHP71) + W(RIC12) & 
	    + W(RIC21) + W(RIC32) + W(RIC40) & 
	    + W(RIC65) - W(RIC66) + W(RIC68)
      CDOT(SC3H3) = W(RR004F) - W(RR004B) + W(RR009F) & 
	    - W(RR009B) + W(RR017F) - W(RR017B) & 
	    + W(RR053F) - W(RR053B) - W(RR054F) & 
	    + W(RR054B) - W(RR055F) + W(RR055B) & 
	    - W(RR056F) + W(RR056B) - W(RR057F) & 
	    + W(RR057B) - W(RR058F) + W(RR058B) & 
	    - W(RR059F) + W(RR059B) - W(RR060F) & 
	    + W(RR060B) - W(RR061F) + W(RR061B) & 
	    - W(RR062F) + W(RR062B) - W(RR063F) & 
	    + W(RR063B) - W(RR064F) + W(RR064B) & 
	    - W(RR065F) + W(RR065B) - W(RR067F) & 
	    + W(RR067B) - W(RR068F) + W(RR068B) & 
	    - W(RR069F) + W(RR069B) - W(RR070F) & 
	    + W(RR070B) + W(RR078F) - W(RR078B) & 
	    + W(RR079F) - W(RR079B) + W(RR080F)
      CDOT(SC3H3) = CDOT(SC3H3) - W(RR080B) + W(RR081F) - W(RR081B) & 
	    + W(RR082F) - W(RR082B) + W(RR083F) & 
	    - W(RR083B) + W(RR084F) - W(RR084B) & 
	    + W(RR085F) - W(RR085B) + W(RR086F) & 
	    - W(RR086B) + W(RR090F) - W(RR090B) & 
	    + W(RR091F) - W(RR091B) + W(RH12F) & 
	    - W(RH12B) + W(RH36F) - W(RH36B) & 
	    - W(RB08F) + W(RB08B) + W(RB45F) & 
	    - W(RB45B) + W(RB51F) - W(RB51B) & 
	    - W(RP008) - 2 * W(RP009F) + 2 * W(RP009B) & 
	    - 2 * W(RP010F) + 2 * W(RP010B) - 2 * W(RP011F) & 
	    + 2 * W(RP011B) - W(RCP16F) + W(RCP16B) & 
	    - W(RI00F) + W(RI00B) - W(RI23) & 
	    + W(RI25) - W(RT20) - W(RXY22) & 
	    - W(RN18)
      CDOT(SPXC3H4) = - W(RR005F) + W(RR005B) + W(RR055F) & 
	    - W(RR055B) + W(RR065F) - W(RR065B) & 
	    - W(RR066F) + W(RR066B) + W(RR068F) & 
	    - W(RR068B) + W(RR071F) - W(RR071B) & 
	    + W(RR073F) - W(RR073B) - W(RR076F) & 
	    + W(RR076B) - W(RR077F) + W(RR077B) & 
	    - W(RR078F) + W(RR078B) - W(RR079F) & 
	    + W(RR079B) - W(RR080F) + W(RR080B) & 
	    - W(RR081F) + W(RR081B) - W(RR082F) & 
	    + W(RR082B) - W(RR088F) + W(RR088B) & 
	    - W(RR089F) + W(RR089B) - W(RR091F) & 
	    + W(RR091B) - W(RR092F) + W(RR092B) & 
	    - W(RR093F) + W(RR093B) - W(RR094F) & 
	    + W(RR094B) + W(RR119F) - W(RR119B) & 
	    + W(RR124F) - W(RR124B) + W(RR125F)
      CDOT(SPXC3H4) = CDOT(SPXC3H4) - W(RR125B) + W(RR130F) - W(RR130B) & 
	    + W(RR131F) - W(RR131B) + W(RR132F) & 
	    - W(RR132B) - W(RB14F) + W(RB14B) & 
	    + W(RB29F) - W(RB29B) + W(RT60)
      CDOT(SAXC3H4) = - W(RR006F) + W(RR006B) + W(RR056F) & 
	    - W(RR056B) + W(RR064F) - W(RR064B) & 
	    + W(RR067F) - W(RR067B) - W(RR071F) & 
	    + W(RR071B) - W(RR073F) + W(RR073B) & 
	    - W(RR074F) + W(RR074B) - W(RR075F) & 
	    + W(RR075B) - W(RR083F) + W(RR083B) & 
	    - W(RR084F) + W(RR084B) - W(RR085F) & 
	    + W(RR085B) - W(RR086F) + W(RR086B) & 
	    - W(RR087F) + W(RR087B) - W(RR090F) & 
	    + W(RR090B) + W(RR103F) - W(RR103B) & 
	    + W(RR104F) - W(RR104B) + W(RR105F) & 
	    - W(RR105B) + W(RR106F) - W(RR106B) & 
	    + W(RR107F) - W(RR107B) + W(RR108F) & 
	    - W(RR108B) + W(RR109F) - W(RR109B) & 
	    + W(RR137F) - W(RR137B) + W(RH34F)
      CDOT(SAXC3H4) = CDOT(SAXC3H4) - W(RH34B) - W(RB15F) + W(RB15B) & 
	    + W(RIC52) + W(RIC56)
      CDOT(SSXC3H5) = W(RR007F) - W(RR007B) + W(RR077F) & 
	    - W(RR077B) + W(RR101F) - W(RR101B) & 
	    + W(RR102F) - W(RR102B) - W(RR125F) & 
	    + W(RR125B) - W(RR126F) + W(RR126B) & 
	    - W(RR127) - W(RR128F) + W(RR128B) & 
	    - W(RR129F) + W(RR129B) - W(RR130F) & 
	    + W(RR130B) - W(RR132F) + W(RR132B) & 
	    - W(RR134) - W(RR136F) + W(RR136B) & 
	    + W(RR153F) - W(RR153B) + W(RR154F) & 
	    - W(RR154B) + W(RR155F) - W(RR155B) & 
	    + W(RR156F) - W(RR156B) + W(RCP04)
      CDOT(SNXC4H3) = W(RR008F) - W(RR008B) + W(RR043F) & 
	    - W(RR043B) + W(RR045F) - W(RR045B) & 
	    + W(RH10F) - W(RH10B) - W(RH13F) & 
	    + W(RH13B) - W(RH14F) + W(RH14B) & 
	    - W(RH15F) + W(RH15B) - W(RH17F) & 
	    + W(RH17B) - W(RH19F) + W(RH19B) & 
	    - W(RH21F) + W(RH21B) - W(RH23F) & 
	    + W(RH23B) + W(RH28F) - W(RH28B) & 
	    + W(RH30F) - W(RH30B) + W(RH32F) & 
	    - W(RH32B) - W(RP005F) + W(RP005B)
      CDOT(SC2H3CHO) = W(RR012F) - W(RR012B) + W(RR057F) & 
	    - W(RR057B) - W(RR095) - W(RR096) & 
	    - W(RR097) - W(RR098) - W(RR099) & 
	    + W(RR110F) - W(RR110B) + W(RR114F) & 
	    - W(RR114B) + W(RR134) + W(RR138) & 
	    + W(RR139F) - W(RR139B) + W(RB43F) & 
	    - W(RB43B) + W(RHP75)
      CDOT(SAXC3H5) = W(RR015F) - W(RR015B) - W(RR016F) & 
	    + W(RR016B) + W(RR074F) - W(RR074B) & 
	    - W(RR100F) + W(RR100B) - W(RR101F) & 
	    + W(RR101B) - W(RR103F) + W(RR103B) & 
	    - W(RR104F) + W(RR104B) - W(RR105F) & 
	    + W(RR105B) - W(RR106F) + W(RR106B) & 
	    - W(RR107F) + W(RR107B) - 2 * W(RR108F) & 
	    + 2 * W(RR108B) - W(RR109F) + W(RR109B) & 
	    - W(RR110F) + W(RR110B) - W(RR111) & 
	    - W(RR112F) + W(RR112B) - W(RR113F) & 
	    + W(RR113B) - W(RR114F) + W(RR114B) & 
	    - W(RR115F) + W(RR115B) - W(RR116F) & 
	    + W(RR116B) - W(RR117F) + W(RR117B) & 
	    - W(RR118F) + W(RR118B) + W(RR144F) & 
	    - W(RR144B) + W(RR145F) - W(RR145B)
      CDOT(SAXC3H5) = CDOT(SAXC3H5) + W(RR146F) - W(RR146B) + W(RR147F) & 
	    - W(RR147B) + W(RR148F) - W(RR148B) & 
	    + W(RB28) + W(RB30F) - W(RB30B) & 
	    + W(RB31F) - W(RB31B) + W(RB99F) & 
	    - W(RB99B) + W(RB42) + W(RHP22) & 
	    + W(RHP27) + W(RHP33) + W(RHP37) & 
	    + W(RHP43) - W(RHP71) + W(RIC23) & 
	    - W(RP008) + W(RCP03F) - W(RCP03B)
      CDOT(SC2O) = - W(RR018F) + W(RR018B) - W(RR019F) & 
	    + W(RR019B) - W(RR020F) + W(RR020B) & 
	    - W(RR021F) + W(RR021B) + W(RR023F) & 
	    - W(RR023B)
      CDOT(SC4H4) = W(RR031F) - W(RR031B) + W(RR044F) & 
	    - W(RR044B) + W(RR070F) - W(RR070B) & 
	    + W(RH04F) - W(RH04B) + W(RH15F) & 
	    - W(RH15B) + W(RH16F) - W(RH16B) & 
	    - W(RH28F) + W(RH28B) - W(RH29F) & 
	    + W(RH29B) - W(RH30F) + W(RH30B) & 
	    - W(RH31F) + W(RH31B) - W(RH32F) & 
	    + W(RH32B) - W(RH33F) + W(RH33B) & 
	    - W(RH34F) + W(RH34B) - W(RH35F) & 
	    + W(RH35B) - W(RH36F) + W(RH36B) & 
	    + W(RB01F) - W(RB01B) + W(RB12F) & 
	    - W(RB12B) - W(RB32F) + W(RB32B) & 
	    - W(RB33F) + W(RB33B) + W(RB36F) & 
	    - W(RB36B) + W(RB37F) - W(RB37B) & 
	    + W(RB41F) - W(RB41B) + W(RB44F)
      CDOT(SC4H4) = CDOT(SC4H4) - W(RB44B) + W(RB46F) - W(RB46B) & 
	    - W(RP107F) + W(RP107B) - W(RP414F) & 
	    + W(RP414B) - W(RP415F) + W(RP415B) & 
	    - W(RP505F) + W(RP505B) + W(RCP33F) & 
	    - W(RCP33B)
      CDOT(SC3H2) = - W(RR039F) + W(RR039B) - W(RR040F) & 
	    + W(RR040B) - W(RR041F) + W(RR041B) & 
	    - W(RR042F) + W(RR042B) - W(RR043F) & 
	    + W(RR043B) - W(RR044F) + W(RR044B) & 
	    - W(RR045F) + W(RR045B) - W(RR053F) & 
	    + W(RR053B) + W(RR054F) - W(RR054B) & 
	    + W(RR059F) - W(RR059B) + W(RH08F) & 
	    - W(RH08B) + W(RH35F) - W(RH35B)
      CDOT(SC3H2O) = W(RR039F) - W(RR039B) + W(RR046F) & 
	    - W(RR046B) - W(RR200F) + W(RR200B) & 
	    - W(RR047) - W(RR048) - W(RR049) & 
	    - W(RR050) - W(RR051) - W(RR052) & 
	    + W(RR061F) - W(RR061B)
      CDOT(SC4H2) = W(RR042F) - W(RR042B) + W(RH02F) & 
	    - W(RH02B) - W(RH03F) + W(RH03B) & 
	    - W(RH04F) + W(RH04B) - 2 * W(RH05) & 
	    - 2 * W(RH06F) + 2 * W(RH06B) - W(RH07F) & 
	    + W(RH07B) - W(RH08F) + W(RH08B) & 
	    - W(RH09F) + W(RH09B) - W(RH10F) & 
	    + W(RH10B) - W(RH11F) + W(RH11B) & 
	    - W(RH12F) + W(RH12B) + W(RH19F) & 
	    - W(RH19B) + W(RH20F) - W(RH20B) & 
	    + W(RH21F) - W(RH21B) + W(RH22F) & 
	    - W(RH22B) + W(RH23F) - W(RH23B) & 
	    + W(RH24F) - W(RH24B) - W(RH37F) & 
	    + W(RH37B) - W(RH40F) + W(RH40B) & 
	    + W(RI25) + W(RI26) - W(ROX02F) & 
	    + W(ROX02B)
      CDOT(SIXC4H3) = W(RR069F) - W(RR069B) + W(RH09F) & 
	    - W(RH09B) + W(RH13F) - W(RH13B) & 
	    + W(RH14F) - W(RH14B) - W(RH16F) & 
	    + W(RH16B) - W(RH18F) + W(RH18B) & 
	    - W(RH20F) + W(RH20B) - W(RH22F) & 
	    + W(RH22B) - W(RH24F) + W(RH24B) & 
	    - W(RH25F) + W(RH25B) - W(RH26F) & 
	    + W(RH26B) - W(RH27F) + W(RH27B) & 
	    + W(RH29F) - W(RH29B) + W(RH31F) & 
	    - W(RH31B) + W(RH33F) - W(RH33B)
      CDOT(STXC3H5) = W(RR075F) - W(RR075B) + W(RR076F) & 
	    - W(RR076B) + W(RR100F) - W(RR100B) & 
	    - W(RR102F) + W(RR102B) - W(RR119F) & 
	    + W(RR119B) - W(RR120F) + W(RR120B) & 
	    - W(RR121) - W(RR122F) + W(RR122B) & 
	    - W(RR123F) + W(RR123B) - W(RR124F) & 
	    + W(RR124B) - W(RR131F) + W(RR131B) & 
	    - W(RR133) - W(RR135) - W(RR137F) & 
	    + W(RR137B) + W(RR149F) - W(RR149B) & 
	    + W(RR150F) - W(RR150B) + W(RR151F) & 
	    - W(RR151B) + W(RR152F) - W(RR152B) & 
	    + W(RIC37) + W(RIC49) + W(RIC59) & 
	    + W(RIC60) + W(RIC61) + W(RIC62) & 
	    + W(RIC63)
      CDOT(SC3H5O) = W(RR113F) - W(RR113B) + W(RR117F) & 
	    - W(RR117B) + W(RR118F) - W(RR118B) & 
	    - W(RR138) - W(RR139F) + W(RR139B) & 
	    - W(RR140)
      CDOT(SC4H) = - W(RH01F) + W(RH01B) - W(RH02F) & 
	    + W(RH02B) + W(RH03F) - W(RH03B) & 
	    + W(RH11F) - W(RH11B) - W(RH38F) & 
	    + W(RH38B) - W(RH40F) + W(RH40B)
      CDOT(SC8H2) = W(RH05) + W(RH06F) - W(RH06B) & 
	    + W(RH39F) - W(RH39B) + W(RH40F) & 
	    - W(RH40B)
      CDOT(SC6H2) = W(RH37F) - W(RH37B) + W(RH38F) & 
	    - W(RH38B) - W(RH39F) + W(RH39B)
      CDOT(SC4H6) = W(RB00F) - W(RB00B) + W(RB04F) & 
	    - W(RB04B) + W(RB08F) - W(RB08B) & 
	    + W(RB09F) - W(RB09B) - W(RB10F) & 
	    + W(RB10B) - W(RB11F) + W(RB11B) & 
	    - W(RB12F) + W(RB12B) + W(RB14F) & 
	    - W(RB14B) + W(RB15F) - W(RB15B) & 
	    - W(RB16F) + W(RB16B) - W(RB17F) & 
	    + W(RB17B) + W(RB18F) - W(RB18B) & 
	    - W(RB19F) + W(RB19B) - W(RB20F) & 
	    + W(RB20B) - W(RB21F) + W(RB21B) & 
	    - W(RB22F) + W(RB22B) - W(RB23F) & 
	    + W(RB23B) - W(RB24F) + W(RB24B) & 
	    - W(RB25F) + W(RB25B) - W(RB28) & 
	    - W(RB29F) + W(RB29B) - W(RB30F) & 
	    + W(RB30B) - W(RB31F) + W(RB31B)
      CDOT(SC4H6) = CDOT(SC4H6) + W(RB38F) - W(RB38B) + W(RB39F) & 
	    - W(RB39B) + W(RB40F) - W(RB40B) & 
	    + W(RB47F) - W(RB47B) + W(RB48F) & 
	    - W(RB48B) + W(RB49F) - W(RB49B) & 
	    + W(RHP38) + W(RHP64F) - W(RHP64B) & 
	    + W(RHP67) + W(RHP68) + W(RHP69) & 
	    + W(RHP70) + W(RHP71) - W(RHP86) & 
	    - W(RHP87) - W(RHP88) + W(RCP11F) & 
	    - W(RCP11B) + W(RCP14)
      CDOT(SNXC4H5) = W(RB02F) - W(RB02B) + W(RB06F) & 
	    - W(RB06B) + W(RB11F) - W(RB11B) & 
	    + W(RB16F) - W(RB16B) - W(RB18F) & 
	    + W(RB18B) + W(RB20F) - W(RB20B) & 
	    + W(RB22F) - W(RB22B) + W(RB24F) & 
	    - W(RB24B) + W(RB32F) - W(RB32B) & 
	    - W(RB34F) + W(RB34B) - W(RB35F) & 
	    + W(RB35B) - W(RB36F) + W(RB36B) & 
	    - W(RB37F) + W(RB37B) - W(RB38F) & 
	    + W(RB38B) - W(RB39F) + W(RB39B) & 
	    - W(RB40F) + W(RB40B) - W(RB99F) & 
	    + W(RB99B) - W(RB41F) + W(RB41B) & 
	    - W(RB42) - W(RB43F) + W(RB43B) & 
	    - W(RB52F) + W(RB52B) - W(RP001F) & 
	    + W(RP001B) - W(RP006F) + W(RP006B)
      CDOT(SNXC4H5) = CDOT(SNXC4H5) - W(RCP11F) + W(RCP11B) + W(RCP24)
      CDOT(SIXC4H5) = W(RB05F) - W(RB05B) + W(RB10F) & 
	    - W(RB10B) + W(RB17F) - W(RB17B) & 
	    + W(RB19F) - W(RB19B) + W(RB21F) & 
	    - W(RB21B) + W(RB23F) - W(RB23B) & 
	    + W(RB25F) - W(RB25B) + W(RB33F) & 
	    - W(RB33B) + W(RB34F) - W(RB34B) & 
	    + W(RB35F) - W(RB35B) - W(RB44F) & 
	    + W(RB44B) - W(RB45F) + W(RB45B) & 
	    - W(RB46F) + W(RB46B) - W(RB47F) & 
	    + W(RB47B) - W(RB48F) + W(RB48B) & 
	    - W(RB49F) + W(RB49B) - W(RB50F) & 
	    + W(RB50B) - W(RB51F) + W(RB51B) & 
	    - W(RP002F) + W(RP002B) - W(RP007F) & 
	    + W(RP007B)
      CDOT(SA1XC6H6) = W(RB52F) - W(RB52B) + W(RP000F) & 
	    - W(RP000B) + W(RP003F) - W(RP003B) & 
	    + W(RP006F) - W(RP006B) + W(RP007F) & 
	    - W(RP007B) + W(RP010F) - W(RP010B) & 
	    - W(RP014F) + W(RP014B) + W(RP015F) & 
	    - W(RP015B) - W(RP301F) + W(RP301B) & 
	    - W(RP417F) + W(RP417B) - W(RP801) & 
	    + W(RT01F) - W(RT01B) + W(RT16F) & 
	    - W(RT16B) + W(RT28F) - W(RT28B) & 
	    + W(RT31F) - W(RT31B) + W(RE04F) & 
	    - W(RE04B) + W(RST01) + W(RST00F) & 
	    - W(RST00B) + W(RXY12) - W(ROX00F) & 
	    + W(ROX00B) - W(ROX03F) + W(ROX03B) & 
	    - W(ROX04F) + W(ROX04B) - W(ROX05F) & 
	    + W(ROX05B) - W(ROX06F) + W(ROX06B)
      CDOT(SA1XC6H6) = CDOT(SA1XC6H6) - W(ROX99F) + W(ROX99B) - W(ROX07F) & 
	    + W(ROX07B) - W(ROX08F) + W(ROX08B) & 
	    - W(ROX09F) + W(ROX09B) - W(ROX10F) & 
	    + W(ROX10B) + W(ROX16F) - W(ROX16B)
      CDOT(SNXC7H16) = - W(RHP00) - W(RHP01) - W(RHP02) & 
	    - W(RHP03) - W(RHP04) - W(RHP05) & 
	    - W(RHP06) + W(RHP08) + W(RHP09)
      CDOT(SC5H11) = W(RHP00) + W(RHP10) + W(RHP20) & 
	    + W(RHP25) - W(RHP29) - W(RHP30) & 
	    + W(RHP31) - W(RHP32)
      CDOT(SPXC4H9) = W(RHP01) + W(RHP11) + W(RHP19) & 
	    + W(RHP22) + W(RHP26) - W(RHP39) & 
	    + W(RHP40) - W(RHP41) - W(RHP42)
      CDOT(SC7H15) = W(RHP02) + W(RHP03) + W(RHP04) & 
	    + W(RHP05) + W(RHP06) - W(RHP08) & 
	    - W(RHP09) - W(RHP10) - W(RHP11) & 
	    - W(RHP12) - W(RHP13) - W(RHP14) & 
	    + W(RHP15) - W(RHP16) - W(RHP17)
      CDOT(SPXC4H8) = W(RHP12) + W(RHP27) + W(RHP39) & 
	    - W(RHP40) - W(RHP43) - W(RHP44) & 
	    - W(RHP45) - W(RHP46) - W(RHP47) & 
	    - W(RHP48) + W(RHP49) + W(RHP50) & 
	    - W(RHP51) - W(RHP52) - W(RHP53) & 
	    - W(RHP54) - W(RHP55) - W(RHP58) & 
	    - W(RHP56) - W(RHP57)
      CDOT(SC5H10) = W(RHP13) + W(RHP30) - W(RHP31) & 
	    - W(RHP33) - W(RHP34) - W(RHP35) & 
	    - W(RHP36)
      CDOT(SC7H14) = W(RHP14) - W(RHP15) - W(RHP21) & 
	    - W(RHP22) - W(RHP23F) + W(RHP23B) & 
	    - W(RHP24) - W(RHP25) - W(RHP26)
      CDOT(SC7H15O) = W(RHP16) + W(RHP17) - W(RHP18) & 
	    - W(RHP19) - W(RHP20)
      CDOT(SC3H7CHO) = W(RHP18) - W(RHP59) - W(RHP60) & 
	    - W(RHP61) - W(RHP62) - W(RHP63)
      CDOT(SC4H7) = W(RHP21) + W(RHP28) + W(RHP44) & 
	    + W(RHP45) + W(RHP46) + W(RHP47) & 
	    + W(RHP48) - W(RHP49) - W(RHP50) & 
	    - W(RHP64F) + W(RHP64B) + W(RHP65F) & 
	    - W(RHP65B) - W(RHP67) - W(RHP68) & 
	    - W(RHP69) - W(RHP70) - W(RHP71) & 
	    - W(RHP72) - W(RHP73)
      CDOT(SC7H13) = W(RHP23F) - W(RHP23B) + W(RHP24) & 
	    - W(RHP27) - W(RHP28)
      CDOT(SC5H9) = W(RHP34) + W(RHP35) + W(RHP36) & 
	    - W(RHP37) - W(RHP38)
      CDOT(SC4H7O) = W(RHP72) + W(RHP73) - W(RHP74) & 
	    - W(RHP75)
      CDOT(SNXC3H7O) = W(RHP76) + W(RHP77) - W(RHP78) & 
	    - W(RHP79)
      CDOT(SIXC8H18) = - W(RIC00) - W(RIC01) - W(RIC02) & 
	    - W(RIC03) - W(RIC04) - W(RIC05) & 
	    - W(RIC06) - W(RIC07) - W(RIC08) & 
	    - W(RIC09)
      CDOT(SYXC7H15) = W(RIC00) + W(RIC17) + W(RIC18) & 
	    - W(RIC19) - W(RIC20) - W(RIC21)
      CDOT(SIXC4H8) = W(RIC01) + W(RIC10) + W(RIC12) & 
	    + W(RIC15) + W(RIC20) + W(RIC33) & 
	    - W(RIC34) - W(RIC37) - W(RIC38F) & 
	    + W(RIC38B) - W(RIC40) - W(RIC41F) & 
	    + W(RIC41B) - W(RIC42) - W(RIC43) & 
	    - W(RIC44) - W(RIC45) - W(RIC46) & 
	    - W(RIC47) - W(RIC48)
      CDOT(SIXC3H7) = W(RIC01) + W(RIC16) + W(RIC20) & 
	    + W(RIC22) + W(RIC28) + W(RIC29) & 
	    + W(RIC30) + W(RIC47) - W(RIC64) & 
	    - W(RIC65) + W(RIC66) + W(RIC67) & 
	    - W(RIC68)
      CDOT(STXC4H9) = 2 * W(RIC02) + W(RIC10) + W(RIC16) & 
	    + W(RIC21) + W(RIC23) - W(RIC32) & 
	    - W(RIC33) + W(RIC34) - W(RIC35) & 
	    - W(RIC36)
      CDOT(SCXC8H17) = W(RIC03) + W(RIC04) + W(RIC05) & 
	    + W(RIC06) + W(RIC07) + W(RIC08) & 
	    + W(RIC09) - W(RIC10) - W(RIC11) & 
	    - W(RIC12) - W(RIC13) - W(RIC14)
      CDOT(SYXC7H14) = W(RIC11) - W(RIC18) + W(RIC19) & 
	    - W(RIC22) - W(RIC23) - W(RIC24) & 
	    + W(RIC25) - W(RIC26) + W(RIC27)
      CDOT(SDXC8H17O) = W(RIC13) + W(RIC14) - W(RIC15) & 
	    - W(RIC16) - W(RIC17)
      CDOT(SCH3COCH3) = W(RIC15) + W(RIC31) - W(RIC69) & 
	    - W(RIC70) - W(RIC71) - W(RIC72)
      CDOT(SIXC4H7) = W(RIC22) + W(RIC38F) - W(RIC38B) & 
	    + W(RIC41F) - W(RIC41B) + W(RIC42) & 
	    + W(RIC43) + W(RIC44) + W(RIC45) & 
	    + W(RIC46) - W(RIC52) - W(RIC57) & 
	    - W(RIC53) - W(RIC54) - W(RIC55) & 
	    - W(RIC56) - W(RIC58)
      CDOT(SXXC7H13) = W(RIC24) - W(RIC25) + W(RIC26) & 
	    - W(RIC27) - W(RIC28)
      CDOT(SIXC3H5CH) = W(RIC28) + W(RIC50) + W(RIC51) & 
	    + W(RIC54) - W(RIC59) - W(RIC60) & 
	    - W(RIC61) - W(RIC62) - W(RIC63)
      CDOT(STXC4H9O) = - W(RIC29) - W(RIC30) - W(RIC31) & 
	    + W(RIC35) + W(RIC36)
      CDOT(SIXC4H7O) = - W(RIC49) - W(RIC50) - W(RIC51) & 
	    + W(RIC57) + W(RIC53) + W(RIC58)
      CDOT(SC5H4CH2) = - W(RP000F) + W(RP000B) + W(RP001F) & 
	    - W(RP001B) + W(RP002F) - W(RP002B) & 
	    - W(RP003F) + W(RP003B) - W(RP004F) & 
	    + W(RP004B) + W(RP008) + W(RP009F) & 
	    - W(RP009B) + W(RCP18) + W(RT47) & 
	    + W(RT59)
      CDOT(SA1XXC6H5) = W(RP004F) - W(RP004B) + W(RP005F) & 
	    - W(RP005B) + W(RP011F) - W(RP011B) & 
	    - W(RK012F) + W(RK012B) - W(RP013F) & 
	    + W(RP013B) - W(RP015F) + W(RP015B) & 
	    - W(RP107F) + W(RP107B) - W(RP301F) & 
	    + W(RP301B) - 2 * W(RP302F) + 2 * W(RP302B) & 
	    - W(RP416F) + W(RP416B) - W(RP418F) & 
	    + W(RP418B) - W(RP800) - W(RP802) & 
	    - W(RI00F) + W(RI00B) + W(RT03F) & 
	    - W(RT03B) + W(RT04F) - W(RT04B) & 
	    - W(RT16F) + W(RT16B) - W(RT28F) & 
	    + W(RT28B) + W(RT30F) - W(RT30B) & 
	    + W(RT36) + W(RT37) + W(RT38) & 
	    + W(RT39) + W(RT40) + W(RT41) & 
	    + W(RT42) - W(RE03F) + W(RE03B)
      CDOT(SA1XXC6H5) = CDOT(SA1XXC6H5) + W(RE11F) - W(RE11B) + W(RE37) & 
	    + W(RXY24) + W(RXY44) + W(RXY50) & 
	    + W(ROX00F) - W(ROX00B) - W(ROX01F) & 
	    + W(ROX01B) + W(ROX03F) - W(ROX03B) & 
	    + W(ROX04F) - W(ROX04B) + W(ROX06F) & 
	    - W(ROX06B) + W(ROX10F) - W(ROX10B) & 
	    - W(ROX11F) + W(ROX11B) - W(ROX12F) & 
	    + W(ROX12B) - W(ROX13F) + W(ROX13B) & 
	    - W(ROX14F) + W(ROX14B) - W(ROX15F) & 
	    + W(ROX15B) - W(ROX16F) + W(ROX16B)
      CDOT(SA1C2H2XC) = W(RK012F) - W(RK012B) - W(RK019F) & 
	    + W(RK019B) - W(RK020F) + W(RK020B) & 
	    - W(RK021F) + W(RK021B) - W(RP022F) & 
	    + W(RP022B) + W(RP026F) - W(RP026B) & 
	    + W(RP027F) - W(RP027B) + W(RP028F) & 
	    - W(RP028B) - W(RP104F) + W(RP104B) & 
	    + W(RST03F) - W(RST03B) - W(RST11F) & 
	    + W(RST11B) - W(RST12F) + W(RST12B) & 
	    - W(RST13) - W(RST14F) + W(RST14B)
      CDOT(SA1C2H3XC) = W(RP013F) - W(RP013B) + W(RP014F) & 
	    - W(RP014B) - W(RP023F) + W(RP023B) & 
	    - W(RK024F) + W(RK024B) - W(RP025F) & 
	    + W(RP025B) - W(RP026F) + W(RP026B) & 
	    - W(RP027F) + W(RP027B) - W(RP028F) & 
	    + W(RP028B) + W(RE12F) - W(RE12B) & 
	    + W(RE13F) - W(RE13B) + W(RE14F) & 
	    - W(RE14B) + W(RE15F) - W(RE15B) & 
	    + W(RE16F) - W(RE16B) + W(RE33) & 
	    - W(RST01) - W(RST02F) + W(RST02B) & 
	    - W(RST03F) + W(RST03B) - W(RST04F) & 
	    + W(RST04B) - W(RST05F) + W(RST05B) & 
	    - W(RST06F) + W(RST06B) - W(RST10F) & 
	    + W(RST10B)
      CDOT(SA1C2HXC8) = - W(RP016F) + W(RP016B) - W(RK017F) & 
	    + W(RK017B) - W(RP018F) + W(RP018B) & 
	    + W(RK020F) - W(RK020B) + W(RK021F) & 
	    - W(RK021B) + W(RP022F) - W(RP022B) & 
	    - W(RP105F) + W(RP105B) - W(RP416F) & 
	    + W(RP416B) + W(RST12F) - W(RST12B) & 
	    - W(RST00F) + W(RST00B)
      CDOT(SA1C2HYXC) = W(RP016F) - W(RP016B) + W(RK017F) & 
	    - W(RK017B) + W(RP018F) - W(RP018B) & 
	    - W(RK100F) + W(RK100B) - W(RP106F) & 
	    + W(RP106B) - W(RP417F) + W(RP417B) & 
	    - W(RP418F) + W(RP418B)
      CDOT(SA1C2H3YX) = W(RK019F) - W(RK019B) + W(RP023F) & 
	    - W(RP023B) + W(RK024F) - W(RK024B) & 
	    + W(RP025F) - W(RP025B) - W(RK102F) & 
	    + W(RK102B) + W(RI32) + W(RST10F) & 
	    - W(RST10B)
      CDOT(SA2XXC10H) = W(RK100F) - W(RK100B) + W(RP108F) & 
	    - W(RP108B) + W(RK109F) - W(RK109B) & 
	    + W(RK110F) - W(RK110B) - W(RK114F) & 
	    + W(RK114B) - W(RP116F) + W(RP116B) & 
	    - W(RP120F) + W(RP120B) - W(RP414F) & 
	    + W(RP414B) - W(RP801) - W(RP802) & 
	    + W(RN05F) - W(RN05B) + W(RN06F) & 
	    - W(RN06B) + W(RN23F) - W(RN23B) & 
	    + W(RN28) + W(RN29) + W(RN30) & 
	    + W(RN31) + W(RN32) + W(RN33) & 
	    + W(RN34) - W(ROX33F) + W(ROX33B) & 
	    - W(ROX35) - W(ROX37F) + W(ROX37B) & 
	    - W(ROX39F) + W(ROX39B) + W(ROX53)
      CDOT(SA2XC10H8) = W(RK102F) - W(RK102B) + W(RP104F) & 
	    - W(RP104B) + W(RP105F) - W(RP105B) & 
	    + W(RP106F) - W(RP106B) + W(RP107F) & 
	    - W(RP107B) - W(RP108F) + W(RP108B) & 
	    - W(RK109F) + W(RK109B) - W(RK110F) & 
	    + W(RK110B) - W(RP111F) + W(RP111B) & 
	    - W(RK112F) + W(RK112B) - W(RK113F) & 
	    + W(RK113B) - W(RP118F) + W(RP118B) & 
	    - W(RP119F) + W(RP119B) - W(RP800) & 
	    + W(RCP17) + W(RI18) + W(RT20) & 
	    + W(RN01F) - W(RN01B) + W(RN14) & 
	    - W(ROX30F) + W(ROX30B) - W(ROX31F) & 
	    + W(ROX31B) - W(ROX32F) + W(ROX32B) & 
	    + W(ROX63)
      CDOT(SA2YXC10H) = W(RP111F) - W(RP111B) + W(RK112F) & 
	    - W(RK112B) + W(RK113F) - W(RK113B) & 
	    - W(RK115F) + W(RK115B) - W(RP117F) & 
	    + W(RP117B) - W(RP121F) + W(RP121B) & 
	    - W(RP415F) + W(RP415B) - W(ROX34F) & 
	    + W(ROX34B) - W(ROX36) - W(ROX38F) & 
	    + W(ROX38B) - W(ROX40F) + W(ROX40B)
      CDOT(SA2C2H2AX) = W(RK114F) - W(RK114B) + W(RP116F) & 
	    - W(RP116B) + W(RP118F) - W(RP118B) & 
	    + W(RP120F) - W(RP120B) - W(RK122F) & 
	    + W(RK122B) - W(RK123F) + W(RK123B) & 
	    - W(RP124F) + W(RP124B) - W(RK200F) & 
	    + W(RK200B) - W(RP405F) + W(RP405B) & 
	    + W(ROX51)
      CDOT(SA2C2H2BX) = W(RK115F) - W(RK115B) + W(RP117F) & 
	    - W(RP117B) + W(RP119F) - W(RP119B) & 
	    + W(RP121F) - W(RP121B) - W(RK125F) & 
	    + W(RK125B) - W(RK126F) + W(RK126B) & 
	    - W(RP127F) + W(RP127B) - W(RP406F) & 
	    + W(RP406B) + W(ROX50)
      CDOT(SA2C2HAXC) = W(RK122F) - W(RK122B) + W(RK123F) & 
	    - W(RK123B) + W(RP124F) - W(RP124B) & 
	    - W(RP128F) + W(RP128B) - W(RK129F) & 
	    + W(RK129B) - W(RP130F) + W(RP130B) & 
	    - W(RP201F) + W(RP201B) - W(RP410F) & 
	    + W(RP410B) + W(ROX60)
      CDOT(SA2C2HBXC) = W(RK125F) - W(RK125B) + W(RK126F) & 
	    - W(RK126B) + W(RP127F) - W(RP127B) & 
	    - W(RP131F) + W(RP131B) - W(RK132F) & 
	    + W(RK132B) - W(RP133F) + W(RP133B) & 
	    - W(RP411F) + W(RP411B) + W(ROX61)
      CDOT(SA2C2HAYX) = W(RP128F) - W(RP128B) + W(RK129F) & 
	    - W(RK129B) + W(RP130F) - W(RP130B) & 
	    - W(RK401F) + W(RK401B) - W(RP408F) & 
	    + W(RP408B) - W(RP412F) + W(RP412B)
      CDOT(SA2C2HBYX) = W(RP131F) - W(RP131B) + W(RK132F) & 
	    - W(RK132B) + W(RP133F) - W(RP133B) & 
	    - W(RK403F) + W(RK403B) - W(RP409F) & 
	    + W(RP409B) - W(RP413F) + W(RP413B)
      CDOT(SA2R5XC12) = W(RK200F) - W(RK200B) + W(RP201F) & 
	    - W(RP201B) - W(RP202F) + W(RP202B) & 
	    - W(RK203F) + W(RK203B) - W(RK204F) & 
	    + W(RK204B) - W(RP207F) + W(RP207B) & 
	    + W(RI23) - W(ROX63)
      CDOT(SA2R5XXC1) = W(RP202F) - W(RP202B) + W(RK203F) & 
	    - W(RK203B) + W(RK204F) - W(RK204B) & 
	    - W(RK205F) + W(RK205B) - W(RP206F) & 
	    + W(RP206B) - W(RP208F) + W(RP208B) & 
	    + W(RP425F) - W(RP425B) - W(RP505F) & 
	    + W(RP505B) - W(ROX53)
      CDOT(SA2R5C2H2) = W(RK205F) - W(RK205B) + W(RP206F) & 
	    - W(RP206B) + W(RP207F) - W(RP207B) & 
	    + W(RP208F) - W(RP208B) - W(RK212F) & 
	    + W(RK212B) - W(RK213F) + W(RK213B) & 
	    - W(RP214F) + W(RP214B) - W(RP502F) & 
	    + W(RP502B)
      CDOT(SA2R5C2HX) = - W(RP209F) + W(RP209B) - W(RK210F) & 
	    + W(RK210B) - W(RP211F) + W(RP211B) & 
	    + W(RK212F) - W(RK212B) + W(RK213F) & 
	    - W(RK213B) + W(RP214F) - W(RP214B)
      CDOT(SA2R5C2HY) = W(RP209F) - W(RP209B) + W(RK210F) & 
	    - W(RK210B) + W(RP211F) - W(RP211B) & 
	    - W(RP501F) + W(RP501B)
      CDOT(SP2XC12H1) = W(RP301F) - W(RP301B) + W(RP302F) & 
	    - W(RP302B) - W(RP304F) + W(RP304B) & 
	    - W(RP305F) + W(RP305B) - W(RP306F) & 
	    + W(RP306B)
      CDOT(SP2XXC12H) = W(RP304F) - W(RP304B) + W(RP305F) & 
	    - W(RP305B) + W(RP306F) - W(RP306B) & 
	    - W(RP407F) + W(RP407B)
      CDOT(SA3XXC14H) = W(RK401F) - W(RK401B) + W(RK403F) & 
	    - W(RK403B) + W(RP419F) - W(RP419B) & 
	    + W(RK420F) - W(RK420B) + W(RK421F) & 
	    - W(RK421B) - W(RP425F) + W(RP425B) & 
	    - W(RK600F) + W(RK600B) - W(ROX50) & 
	    + W(ROX52) + W(ROX54)
      CDOT(SA3XC14H1) = W(RP405F) - W(RP405B) + W(RP406F) & 
	    - W(RP406B) + W(RP407F) - W(RP407B) & 
	    + W(RP408F) - W(RP408B) + W(RP409F) & 
	    - W(RP409B) + W(RP410F) - W(RP410B) & 
	    + W(RP411F) - W(RP411B) + W(RP412F) & 
	    - W(RP412B) + W(RP413F) - W(RP413B) & 
	    + W(RP414F) - W(RP414B) + W(RP415F) & 
	    - W(RP415B) + W(RP416F) - W(RP416B) & 
	    + W(RP417F) - W(RP417B) + W(RP418F) & 
	    - W(RP418B) - W(RP419F) + W(RP419B) & 
	    - W(RK420F) + W(RK420B) - W(RK421F) & 
	    + W(RK421B) - W(RP422F) + W(RP422B) & 
	    - W(RK423F) + W(RK423B) - W(RK424F) & 
	    + W(RK424B) + W(RI17) + W(RN18) & 
	    - W(ROX60) - W(ROX61) + W(ROX62)
      CDOT(SA3XC14H1) = CDOT(SA3XC14H1) + W(ROX64) + W(ROX66)
      CDOT(SA3YXC14H) = W(RP422F) - W(RP422B) + W(RK423F) & 
	    - W(RK423B) + W(RK424F) - W(RK424B) & 
	    - W(RP503F) + W(RP503B) - W(RP504F) & 
	    + W(RP504B) - W(ROX51)
      CDOT(SA3R5XXC1) = W(RP501F) - W(RP501B) + W(RP506F) & 
	    - W(RP506B) + W(RK507F) - W(RK507B) & 
	    + W(RK508F) - W(RK508B) - W(RK701F) & 
	    + W(RK701B) - W(ROX54)
      CDOT(SA3R5XC16) = W(RP502F) - W(RP502B) + W(RP503F) & 
	    - W(RP503B) + W(RP504F) - W(RP504B) & 
	    + W(RP505F) - W(RP505B) - W(RP506F) & 
	    + W(RP506B) - W(RK507F) + W(RK507B) & 
	    - W(RK508F) + W(RK508B) - W(ROX64)
      CDOT(SA4XC16H1) = W(RK600F) - W(RK600B) - W(RP601F) & 
	    + W(RP601B) - W(RK602F) + W(RK602B) & 
	    - W(RK603F) + W(RK603B) - W(ROX62) & 
	    + W(ROX65)
      CDOT(SA4XXC16H) = W(RP601F) - W(RP601B) + W(RK602F) & 
	    - W(RK602B) + W(RK603F) - W(RK603B) & 
	    - W(RK700F) + W(RK700B) - W(ROX52)
      CDOT(SA4R5XC18) = W(RK700F) - W(RK700B) + W(RK701F) & 
	    - W(RK701B) - W(ROX65)
      CDOT(SFLTNXC16) = W(RP800) + W(RP801) + W(RP802) & 
	    - W(ROX66)
      CDOT(SC5H6) = - W(RCP01F) + W(RCP01B) - W(RCP02F) & 
	    + W(RCP02B) - W(RCP03F) + W(RCP03B) & 
	    - W(RCP04) - W(RCP05F) + W(RCP05B) & 
	    - W(RCP06F) + W(RCP06B) - W(RCP07F) & 
	    + W(RCP07B) - W(RCP08F) + W(RCP08B) & 
	    - W(RCP09F) + W(RCP09B) - W(RCP10F) & 
	    + W(RCP10B) - W(RCP11F) + W(RCP11B) & 
	    - W(RCP12F) + W(RCP12B) - W(RCP13) & 
	    - W(RCP14) - W(RCP15) - W(RCP34) & 
	    + W(ROX09F) - W(ROX09B) + W(ROX17F) & 
	    - W(ROX17B)
      CDOT(SC5H5) = W(RCP01F) - W(RCP01B) + W(RCP02F) & 
	    - W(RCP02B) + W(RCP05F) - W(RCP05B) & 
	    + W(RCP06F) - W(RCP06B) + W(RCP07F) & 
	    - W(RCP07B) + W(RCP08F) - W(RCP08B) & 
	    + W(RCP09F) - W(RCP09B) + W(RCP10F) & 
	    - W(RCP10B) + W(RCP11F) - W(RCP11B) & 
	    + W(RCP16F) - W(RCP16B) - 2 * W(RCP17) & 
	    - W(RCP18) - W(RCP19F) + W(RCP19B) & 
	    - W(RCP20F) + W(RCP20B) - W(RCP21F) & 
	    + W(RCP21B) - W(RCP22) - W(RCP34) & 
	    - W(RI17) + W(RI26) + W(RT05F) & 
	    - W(RT05B) + W(ROX23F) - W(ROX23B)
      CDOT(STXC5H5O) = W(RCP12F) - W(RCP12B) - W(RCP24) & 
	    - W(RCP29) - W(RCP30) + W(RCP32F) & 
	    - W(RCP32B) + W(ROX27F) - W(ROX27B)
      CDOT(SC5H4O) = W(RCP19F) - W(RCP19B) + W(RCP20F) & 
	    - W(RCP20B) + W(RCP23F) - W(RCP23B) & 
	    + W(RCP25) + W(RCP26) + W(RCP27) & 
	    + W(RCP28) + W(RCP29) + W(RCP30) & 
	    - W(RCP31) - W(RCP32F) + W(RCP32B) & 
	    - W(RCP33F) + W(RCP33B) + W(ROX26F) & 
	    - W(ROX26B)
      CDOT(SSXC5H5O) = W(RCP21F) - W(RCP21B) + W(RCP22) & 
	    - W(RCP23F) + W(RCP23B) - W(RCP25) & 
	    - W(RCP26) - W(RCP27) - W(RCP28)
      CDOT(SC9H8) = W(RCP34) + W(RI00F) - W(RI00B) & 
	    - W(RI01F) + W(RI01B) - W(RI02F) & 
	    + W(RI02B) + W(RI03F) - W(RI03B) & 
	    - W(RI05F) + W(RI05B) - W(RI06F) & 
	    + W(RI06B) - W(RI07F) + W(RI07B) & 
	    - W(RI08F) + W(RI08B) - W(RI09F) & 
	    + W(RI09B) - W(RI12) - W(RI15) & 
	    + W(ROX41F) - W(ROX41B)
      CDOT(SC9H7) = W(RI01F) - W(RI01B) + W(RI02F) & 
	    - W(RI02B) + W(RI05F) - W(RI05B) & 
	    + W(RI06F) - W(RI06B) + W(RI07F) & 
	    - W(RI07B) + W(RI08F) - W(RI08B) & 
	    + W(RI09F) - W(RI09B) - W(RI17) & 
	    - W(RI18) - W(RI19F) + W(RI19B) & 
	    - W(RI20F) + W(RI20B) - W(RI21) & 
	    - W(RI22) - W(RI23) - W(RI25) & 
	    - W(RI26) + W(RN07F) - W(RN07B) & 
	    + W(RN15) + W(ROX46F) - W(ROX46B)
      CDOT(SA1CH2XC7) = - W(RI03F) + W(RI03B) + W(RT02F) & 
	    - W(RT02B) - W(RT04F) + W(RT04B) & 
	    - W(RT05F) + W(RT05B) + W(RT06F) & 
	    - W(RT06B) + W(RT07F) - W(RT07B) & 
	    + W(RT08F) - W(RT08B) + W(RT11F) & 
	    - W(RT11B) + W(RT14F) - W(RT14B) & 
	    + W(RT15F) - W(RT15B) + W(RT16F) & 
	    - W(RT16B) - W(RT17F) + W(RT17B) & 
	    - W(RT18F) + W(RT18B) - W(RT19F) & 
	    + W(RT19B) - W(RT20) - W(RT21F) & 
	    + W(RT21B) - W(RT22F) + W(RT22B) & 
	    - W(RE02F) + W(RE02B) + W(RE17F) & 
	    - W(RE17B) + W(RE19) + W(RE31) & 
	    + W(RE34) + W(RST02F) - W(RST02B) & 
	    + W(RST05F) - W(RST05B) + W(RST11F)
      CDOT(SA1CH2XC7) = CDOT(SA1CH2XC7) - W(RST11B) + W(RST13)
      CDOT(SC9H6O) = W(RI12) + W(RI19F) - W(RI19B) & 
	    + W(RI20F) - W(RI20B) + W(RI21) & 
	    + W(RI22) - W(RI31) - W(RI32) & 
	    + W(ROX35) + W(ROX36) + W(ROX47) & 
	    + W(ROX48)
      CDOT(SOXC6H4) = W(RI15) + W(RI31) + W(ROX01F) & 
	    - W(ROX01B) + W(ROX02F) - W(ROX02B)
      CDOT(SA1CH3XC7) = - W(RT01F) + W(RT01B) - W(RT02F) & 
	    + W(RT02B) - W(RT03F) + W(RT03B) & 
	    - W(RT06F) + W(RT06B) - W(RT07F) & 
	    + W(RT07B) - W(RT08F) + W(RT08B) & 
	    - W(RT09F) + W(RT09B) - W(RT10F) & 
	    + W(RT10B) - W(RT11F) + W(RT11B) & 
	    - W(RT12F) + W(RT12B) - W(RT13F) & 
	    + W(RT13B) - W(RT14F) + W(RT14B) & 
	    - W(RT15F) + W(RT15B) - W(RT16F) & 
	    + W(RT16B) - W(RT50F) + W(RT50B) & 
	    - W(RT51F) + W(RT51B) - W(RT52F) & 
	    + W(RT52B) - W(RT53F) + W(RT53B) & 
	    - W(RT54F) + W(RT54B) + W(RXY09F) & 
	    - W(RXY09B) + W(RXY11) + W(RXY16F) & 
	    - W(RXY16B) + W(RXY203) + W(RXY39F)
      CDOT(SA1CH3XC7) = CDOT(SA1CH3XC7) - W(RXY39B)
      CDOT(SA1OHXC6H) = W(RT09F) - W(RT09B) + W(RE05F) & 
	    - W(RE05B) + W(RST06F) - W(RST06B) & 
	    + W(ROX05F) - W(ROX05B) + W(ROX08F) & 
	    - W(ROX08B) - W(ROX17F) + W(ROX17B) & 
	    - W(ROX18F) + W(ROX18B) - W(ROX19F) & 
	    + W(ROX19B) - W(ROX20F) + W(ROX20B) & 
	    - W(ROX21F) + W(ROX21B) - W(ROX22F) & 
	    + W(ROX22B)
      CDOT(SHOA1CH3X) = W(RT10F) - W(RT10B) + W(RT12F) & 
	    - W(RT12B) - W(RT43F) + W(RT43B) & 
	    - W(RT44F) + W(RT44B) - W(RT45F) & 
	    + W(RT45B) - W(RT46F) + W(RT46B) & 
	    + W(RXY10F) - W(RXY10B) + W(RXY41F) & 
	    - W(RXY41B)
      CDOT(SOA1CH3XC) = W(RT13F) - W(RT13B) + W(RT43F) & 
	    - W(RT43B) + W(RT44F) - W(RT44B) & 
	    + W(RT45F) - W(RT45B) + W(RT46F) & 
	    - W(RT46B) - W(RT47) + W(RT55F) & 
	    - W(RT55B) + W(RT56F) - W(RT56B) & 
	    + W(RT57F) - W(RT57B) + W(RT58F) & 
	    - W(RT58B) + W(RXY19F) - W(RXY19B)
      CDOT(SA1CH2OXC) = W(RT17F) - W(RT17B) + W(RT19F) & 
	    - W(RT19B) - W(RT23F) + W(RT23B) & 
	    + W(RT24F) - W(RT24B) + W(RT25F) & 
	    - W(RT25B) + W(RT26F) - W(RT26B) & 
	    + W(RT27F) - W(RT27B) - W(RT29F) & 
	    + W(RT29B) - W(RT30F) + W(RT30B) & 
	    - W(RT31F) + W(RT31B) - W(RT32F) & 
	    + W(RT32B) - W(RT33F) + W(RT33B) & 
	    - W(RT34F) + W(RT34B) - W(RT35F) & 
	    + W(RT35B)
      CDOT(SA1CH2OHX) = W(RT18F) - W(RT18B) + W(RT23F) & 
	    - W(RT23B) - W(RT24F) + W(RT24B) & 
	    - W(RT25F) + W(RT25B) - W(RT26F) & 
	    + W(RT26B) - W(RT27F) + W(RT27B)
      CDOT(SA1CHOXC7) = W(RT21F) - W(RT21B) + W(RT29F) & 
	    - W(RT29B) + W(RT32F) - W(RT32B) & 
	    + W(RT33F) - W(RT33B) + W(RT34F) & 
	    - W(RT34B) + W(RT35F) - W(RT35B) & 
	    - W(RT36) - W(RT37) - W(RT38) & 
	    - W(RT39) - W(RT40) - W(RT41) & 
	    - W(RT42) + W(RE18F) - W(RE18B) & 
	    + W(RST04F) - W(RST04B) + W(RST14F) & 
	    - W(RST14B) + W(RXY40F) - W(RXY40B) & 
	    + W(RXY45F) - W(RXY45B) + W(RXY51) & 
	    + W(RXY52) + W(RXY53) + W(RXY54) & 
	    + W(RXY55) + W(RXY56) + W(RXY57F) & 
	    - W(RXY57B)
      CDOT(SA1OXC6H5) = W(RT22F) - W(RT22B) + W(RXY42) & 
	    + W(RXY58) + W(ROX99F) - W(ROX99B) & 
	    + W(ROX07F) - W(ROX07B) + W(ROX11F) & 
	    - W(ROX11B) + W(ROX13F) - W(ROX13B) & 
	    + W(ROX14F) - W(ROX14B) + W(ROX15F) & 
	    - W(ROX15B) + W(ROX18F) - W(ROX18B) & 
	    + W(ROX19F) - W(ROX19B) + W(ROX20F) & 
	    - W(ROX20B) + W(ROX21F) - W(ROX21B) & 
	    + W(ROX22F) - W(ROX22B) - W(ROX23F) & 
	    + W(ROX23B) - W(ROX24F) + W(ROX24B) & 
	    - W(ROX25F) + W(ROX25B)
      CDOT(SA1CH3YXC) = W(RT50F) - W(RT50B) + W(RT51F) & 
	    - W(RT51B) + W(RT52F) - W(RT52B) & 
	    + W(RT53F) - W(RT53B) + W(RT54F) & 
	    - W(RT54B) - W(RT55F) + W(RT55B) & 
	    - W(RT56F) + W(RT56B) - W(RT57F) & 
	    + W(RT57B) - W(RT58F) + W(RT58B) & 
	    - W(RT59) - W(RT60) + W(RXY01F) & 
	    - W(RXY01B) + W(RXY13F) - W(RXY13B) & 
	    + W(RXY15F) - W(RXY15B) + W(RXY202) & 
	    + W(RXY25) + W(RXY33) + W(RXY34) & 
	    + W(RXY35) + W(RXY36) + W(RXY37) & 
	    + W(RXY38)
      CDOT(SA1C2H4XC) = - W(RE01F) + W(RE01B) + W(RE06F) & 
	    - W(RE06B) + W(RE07F) - W(RE07B) & 
	    + W(RE08F) - W(RE08B) + W(RE09F) & 
	    - W(RE09B) + W(RE10F) - W(RE10B) & 
	    - W(RE11F) + W(RE11B) - W(RE12F) & 
	    + W(RE12B) - W(RE13F) + W(RE13B) & 
	    - W(RE14F) + W(RE14B) - W(RE15F) & 
	    + W(RE15B) - W(RE16F) + W(RE16B) & 
	    - W(RE17F) + W(RE17B) - W(RE18F) & 
	    + W(RE18B) - W(RE19) - W(RE30) & 
	    - W(RE31) + W(RE32)
      CDOT(SA1C2H5XC) = W(RE01F) - W(RE01B) + W(RE02F) & 
	    - W(RE02B) + W(RE03F) - W(RE03B) & 
	    - W(RE04F) + W(RE04B) - W(RE05F) & 
	    + W(RE05B) - W(RE06F) + W(RE06B) & 
	    - W(RE07F) + W(RE07B) - W(RE08F) & 
	    + W(RE08B) - W(RE09F) + W(RE09B) & 
	    - W(RE10F) + W(RE10B)
      CDOT(SC8H9O2) = W(RE30) - W(RE32) - W(RE33) & 
	    - W(RE34) - W(RE35)
      CDOT(SC8H8OOH) = W(RE35) - W(RE36)
      CDOT(SOC8H7OOH) = W(RE36) - W(RE37)
      CDOT(SA1CH3CH3) = - W(RXY00F) + W(RXY00B) - W(RXY01F) & 
	    + W(RXY01B) - W(RXY02F) + W(RXY02B) & 
	    - W(RXY03F) + W(RXY03B) - W(RXY04F) & 
	    + W(RXY04B) - W(RXY05F) + W(RXY05B) & 
	    - W(RXY06F) + W(RXY06B) - W(RXY07F) & 
	    + W(RXY07B) - W(RXY09F) + W(RXY09B) & 
	    - W(RXY10F) + W(RXY10B) - W(RXY11)
      CDOT(SA1CH3CH2) = W(RXY00F) - W(RXY00B) + W(RXY02F) & 
	    - W(RXY02B) + W(RXY03F) - W(RXY03B) & 
	    + W(RXY04F) - W(RXY04B) + W(RXY05F) & 
	    - W(RXY05B) + W(RXY06F) - W(RXY06B) & 
	    + W(RXY07F) - W(RXY07B) - W(RXY12) & 
	    - W(RXY13F) + W(RXY13B) - W(RXY14F) & 
	    + W(RXY14B) - W(RXY15F) + W(RXY15B) & 
	    - W(RXY16F) + W(RXY16B) - W(RXY17) & 
	    - W(RXY18F) + W(RXY18B) - W(RXY19F) & 
	    + W(RXY19B) - W(RXY201) - W(RXY202) & 
	    - W(RXY203) - W(RXY22)
      CDOT(SA1CH3CHO) = W(RXY14F) - W(RXY14B) + W(RXY17) & 
	    + W(RXY18F) - W(RXY18B) + W(RXY201) & 
	    - W(RXY23F) + W(RXY23B) - W(RXY24) & 
	    - W(RXY25) - W(RXY26F) + W(RXY26B) & 
	    - W(RXY27F) + W(RXY27B) - W(RXY28F) & 
	    + W(RXY28B) - W(RXY29F) + W(RXY29B) & 
	    - W(RXY30F) + W(RXY30B) - W(RXY31F) & 
	    + W(RXY31B) - W(RXY33) - W(RXY34) & 
	    - W(RXY35) - W(RXY36) - W(RXY37) & 
	    - W(RXY38) - W(RXY39F) + W(RXY39B) & 
	    - W(RXY40F) + W(RXY40B) - W(RXY41F) & 
	    + W(RXY41B) - W(RXY42)
      CDOT(SA2CH3XC1) = W(RXY22) - W(RN01F) + W(RN01B) & 
	    - W(RN02F) + W(RN02B) - W(RN04F) & 
	    + W(RN04B) - W(RN05F) + W(RN05B) & 
	    - W(RN08F) + W(RN08B) - W(RN09F) & 
	    + W(RN09B) - W(RN10F) + W(RN10B) & 
	    - W(RN11F) + W(RN11B) - W(RN12F) & 
	    + W(RN12B) - W(RN13F) + W(RN13B) & 
	    - W(RN14) - W(RN15)
      CDOT(SA1CHOCH2) = W(RXY23F) - W(RXY23B) + W(RXY26F) & 
	    - W(RXY26B) + W(RXY27F) - W(RXY27B) & 
	    + W(RXY28F) - W(RXY28B) + W(RXY29F) & 
	    - W(RXY29B) + W(RXY30F) - W(RXY30B) & 
	    + W(RXY31F) - W(RXY31B) - W(RXY43F) & 
	    + W(RXY43B) - W(RXY44) - W(RXY45F) & 
	    + W(RXY45B) - W(RXY46) - W(RXY47F) & 
	    + W(RXY47B) - W(RXY48)
      CDOT(SA1CHOCHO) = W(RXY43F) - W(RXY43B) + W(RXY46) & 
	    + W(RXY47F) - W(RXY47B) + W(RXY48) & 
	    - W(RXY50) - W(RXY51) - W(RXY52) & 
	    - W(RXY53) - W(RXY54) - W(RXY55) & 
	    - W(RXY56) - W(RXY57F) + W(RXY57B) & 
	    - W(RXY58)
      CDOT(SA2OHXC10) = W(RN02F) - W(RN02B) + W(ROX31F) & 
	    - W(ROX31B) + W(ROX32F) - W(ROX32B) & 
	    - W(ROX41F) + W(ROX41B) - W(ROX42F) & 
	    + W(ROX42B) - W(ROX43F) + W(ROX43B) & 
	    - W(ROX44F) + W(ROX44B) - W(ROX45F) & 
	    + W(ROX45B)
      CDOT(SA2CH2XC1) = W(RN04F) - W(RN04B) - W(RN06F) & 
	    + W(RN06B) - W(RN07F) + W(RN07B) & 
	    + W(RN08F) - W(RN08B) + W(RN09F) & 
	    - W(RN09B) + W(RN10F) - W(RN10B) & 
	    + W(RN11F) - W(RN11B) + W(RN12F) & 
	    - W(RN12B) + W(RN13F) - W(RN13B) & 
	    - W(RN16F) + W(RN16B) - W(RN17) & 
	    - W(RN18F) + W(RN18B) - W(RN18) & 
	    - W(RN20F) + W(RN20B) - W(RN21F) & 
	    + W(RN21B)
      CDOT(SA2CH2OXC) = W(RN16F) - W(RN16B) + W(RN17) & 
	    + W(RN18F) - W(RN18B) - W(RN22F) & 
	    + W(RN22B) - W(RN23F) + W(RN23B) & 
	    - W(RN24F) + W(RN24B) - W(RN25F) & 
	    + W(RN25B) - W(RN26F) + W(RN26B) & 
	    - W(RN27F) + W(RN27B)
      CDOT(SA2CHOXC1) = W(RN20F) - W(RN20B) + W(RN22F) & 
	    - W(RN22B) + W(RN24F) - W(RN24B) & 
	    + W(RN25F) - W(RN25B) + W(RN26F) & 
	    - W(RN26B) + W(RN27F) - W(RN27B) & 
	    - W(RN28) - W(RN29) - W(RN30) & 
	    - W(RN31) - W(RN32) - W(RN33) & 
	    - W(RN34)
      CDOT(SA2OXC10H) = W(RN21F) - W(RN21B) + W(ROX30F) & 
	    - W(ROX30B) + W(ROX33F) - W(ROX33B) & 
	    + W(ROX34F) - W(ROX34B) + W(ROX37F) & 
	    - W(ROX37B) + W(ROX38F) - W(ROX38B) & 
	    + W(ROX39F) - W(ROX39B) + W(ROX40F) & 
	    - W(ROX40B) + W(ROX42F) - W(ROX42B) & 
	    + W(ROX43F) - W(ROX43B) + W(ROX44F) & 
	    - W(ROX44B) + W(ROX45F) - W(ROX45B) & 
	    - W(ROX46F) + W(ROX46B) - W(ROX47) & 
	    - W(ROX48)
      CDOT(SOC6H4O) = W(ROX12F) - W(ROX12B) + W(ROX24F) & 
	    - W(ROX24B) + W(ROX25F) - W(ROX25B) & 
	    - W(ROX26F) + W(ROX26B) - W(ROX27F) & 
	    + W(ROX27B) - W(ROX28)
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

      DOUBLE PRECISION MM(158)
      INCLUDE 'TheSootF90.h'

      MM(SN2) =  2.80200000e+01
      MM(SH) =  1.00800000e+00
      MM(SO2) =  3.20000000e+01
      MM(SO) =  1.60000000e+01
      MM(SOH) =  1.70080000e+01
      MM(SH2) =  2.01600000e+00
      MM(SH2O) =  1.80160000e+01
      MM(SCO2) =  4.40100000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SH2O2) =  3.40160000e+01
      MM(SCO) =  2.80100000e+01
      MM(SHCO) =  2.90180000e+01
      MM(SC) =  1.20100000e+01
      MM(SCH) =  1.30180000e+01
      MM(STXCH2) =  1.40260000e+01
      MM(SCH3) =  1.50340000e+01
      MM(SCH2O) =  3.00260000e+01
      MM(SHCCO) =  4.10280000e+01
      MM(SC2H) =  2.50280000e+01
      MM(SCH2CO) =  4.20360000e+01
      MM(SC2H2) =  2.60360000e+01
      MM(SSXCH2) =  1.40260000e+01
      MM(SAR) =  3.99480000e+01
      MM(SCH3OH) =  3.20420000e+01
      MM(SCH2OH) =  3.10340000e+01
      MM(SCH3O) =  3.10340000e+01
      MM(SCH4) =  1.60420000e+01
      MM(SCH3O2) =  4.70340000e+01
      MM(SC2H3) =  2.70440000e+01
      MM(SC2H4) =  2.80520000e+01
      MM(SC2H5) =  2.90600000e+01
      MM(SHCCOH) =  4.20360000e+01
      MM(SCH2CHO) =  4.30440000e+01
      MM(SCH3CHO) =  4.40520000e+01
      MM(SH2C2) =  2.60360000e+01
      MM(SC2H5O) =  4.50600000e+01
      MM(SNXC3H7) =  4.30860000e+01
      MM(SC2H6) =  3.00680000e+01
      MM(SC3H8) =  4.40940000e+01
      MM(SC3H6) =  4.20780000e+01
      MM(SC3H3) =  3.90540000e+01
      MM(SPXC3H4) =  4.00620000e+01
      MM(SAXC3H4) =  4.00620000e+01
      MM(SSXC3H5) =  4.10700000e+01
      MM(SNXC4H3) =  5.10640000e+01
      MM(SC2H3CHO) =  5.60620000e+01
      MM(SAXC3H5) =  4.10700000e+01
      MM(SC2O) =  4.00200000e+01
      MM(SC4H4) =  5.20720000e+01
      MM(SC3H2) =  3.80460000e+01
      MM(SC3H2O) =  5.40460000e+01
      MM(SC4H2) =  5.00560000e+01
      MM(SIXC4H3) =  5.10640000e+01
      MM(STXC3H5) =  4.10700000e+01
      MM(SC3H5O) =  5.70700000e+01
      MM(SC4H) =  4.90480000e+01
      MM(SC8H2) =  9.80960000e+01
      MM(SC6H2) =  7.40760000e+01
      MM(SC4H6) =  5.40880000e+01
      MM(SNXC4H5) =  5.30800000e+01
      MM(SIXC4H5) =  5.30800000e+01
      MM(SA1XC6H6) =  7.81080000e+01
      MM(SNXC7H16) =  1.00198000e+02
      MM(SC5H11) =  7.11380000e+01
      MM(SPXC4H9) =  5.71120000e+01
      MM(SC7H15) =  9.91900000e+01
      MM(SPXC4H8) =  5.61040000e+01
      MM(SC5H10) =  7.01300000e+01
      MM(SC7H14) =  9.81820000e+01
      MM(SC7H15O) =  1.15190000e+02
      MM(SC3H7CHO) =  7.21040000e+01
      MM(SC4H7) =  5.50960000e+01
      MM(SC7H13) =  9.71740000e+01
      MM(SC5H9) =  6.91220000e+01
      MM(SC4H7O) =  7.10960000e+01
      MM(SNXC3H7O) =  5.90860000e+01
      MM(SIXC8H18) =  1.14224000e+02
      MM(SYXC7H15) =  9.91900000e+01
      MM(SIXC4H8) =  5.61040000e+01
      MM(SIXC3H7) =  4.30860000e+01
      MM(STXC4H9) =  5.71120000e+01
      MM(SCXC8H17) =  1.13216000e+02
      MM(SYXC7H14) =  9.81820000e+01
      MM(SDXC8H17O) =  1.29216000e+02
      MM(SCH3COCH3) =  5.80780000e+01
      MM(SIXC4H7) =  5.50960000e+01
      MM(SXXC7H13) =  9.71740000e+01
      MM(SIXC3H5CH) =  7.00880000e+01
      MM(STXC4H9O) =  7.31120000e+01
      MM(SIXC4H7O) =  7.10960000e+01
      MM(SC5H4CH2) =  7.81080000e+01
      MM(SA1XXC6H5) =  7.71000000e+01
      MM(SA1C2H2XC) =  1.03136000e+02
      MM(SA1C2H3XC) =  1.04144000e+02
      MM(SA1C2HXC8) =  1.02128000e+02
      MM(SA1C2HYXC) =  1.01120000e+02
      MM(SA1C2H3YX) =  1.03136000e+02
      MM(SA2XXC10H) =  1.27156000e+02
      MM(SA2XC10H8) =  1.28164000e+02
      MM(SA2YXC10H) =  1.27156000e+02
      MM(SA2C2H2AX) =  1.53192000e+02
      MM(SA2C2H2BX) =  1.53192000e+02
      MM(SA2C2HAXC) =  1.52184000e+02
      MM(SA2C2HBXC) =  1.52184000e+02
      MM(SA2C2HAYX) =  1.51176000e+02
      MM(SA2C2HBYX) =  1.51176000e+02
      MM(SA2R5XC12) =  1.52184000e+02
      MM(SA2R5XXC1) =  1.51176000e+02
      MM(SA2R5C2H2) =  1.77212000e+02
      MM(SA2R5C2HX) =  1.76204000e+02
      MM(SA2R5C2HY) =  1.75196000e+02
      MM(SP2XC12H1) =  1.54200000e+02
      MM(SP2XXC12H) =  1.53192000e+02
      MM(SA3XXC14H) =  1.77212000e+02
      MM(SA3XC14H1) =  1.78220000e+02
      MM(SA3YXC14H) =  1.77212000e+02
      MM(SA3R5XXC1) =  2.01232000e+02
      MM(SA3R5XC16) =  2.02240000e+02
      MM(SA4XC16H1) =  2.02240000e+02
      MM(SA4XXC16H) =  2.01232000e+02
      MM(SA4R5XC18) =  2.26260000e+02
      MM(SFLTNXC16) =  2.02240000e+02
      MM(SC5H6) =  6.60980000e+01
      MM(SC5H5) =  6.50900000e+01
      MM(STXC5H5O) =  8.10900000e+01
      MM(SC5H4O) =  8.00820000e+01
      MM(SSXC5H5O) =  8.10900000e+01
      MM(SC9H8) =  1.16154000e+02
      MM(SC9H7) =  1.15146000e+02
      MM(SA1CH2XC7) =  9.11260000e+01
      MM(SC9H6O) =  1.30138000e+02
      MM(SOXC6H4) =  7.60920000e+01
      MM(SA1CH3XC7) =  9.21340000e+01
      MM(SA1OHXC6H) =  9.41080000e+01
      MM(SHOA1CH3X) =  1.08134000e+02
      MM(SOA1CH3XC) =  1.07126000e+02
      MM(SA1CH2OXC) =  1.07126000e+02
      MM(SA1CH2OHX) =  1.08134000e+02
      MM(SA1CHOXC7) =  1.06118000e+02
      MM(SA1OXC6H5) =  9.31000000e+01
      MM(SA1CH3YXC) =  9.11260000e+01
      MM(SA1C2H4XC) =  1.05152000e+02
      MM(SA1C2H5XC) =  1.06160000e+02
      MM(SC8H9O2) =  1.37152000e+02
      MM(SC8H8OOH) =  1.37152000e+02
      MM(SOC8H7OOH) =  1.52144000e+02
      MM(SA1CH3CH3) =  1.06160000e+02
      MM(SA1CH3CH2) =  1.05152000e+02
      MM(SA1CH3CHO) =  1.20144000e+02
      MM(SA2CH3XC1) =  1.42190000e+02
      MM(SA1CHOCH2) =  1.19136000e+02
      MM(SA1CHOCHO) =  1.34128000e+02
      MM(SA2OHXC10) =  1.44164000e+02
      MM(SA2CH2XC1) =  1.41182000e+02
      MM(SA2CH2OXC) =  1.57182000e+02
      MM(SA2CHOXC1) =  1.56174000e+02
      MM(SA2OXC10H) =  1.43156000e+02
      MM(SOC6H4O) =  1.08092000e+02

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(158)
      INCLUDE 'TheSootF90.h'

!!$      NAMES(SN2)='N2                  '
!!$      NAMES(SH)='H                   '
!!$      NAMES(SO2)='O2                  '
!!$      NAMES(SO)='O                   '
!!$      NAMES(SOH)='OH                  '
!!$      NAMES(SH2)='H2                  '
!!$      NAMES(SH2O)='H2O                 '
!!$      NAMES(SCO2)='CO2                 '
!!$      NAMES(SHO2)='HO2                 '
!!$      NAMES(SH2O2)='H2O2                '
!!$      NAMES(SCO)='CO                  '
!!$      NAMES(SHCO)='HCO                 '
!!$      NAMES(SC)='C                   '
!!$      NAMES(SCH)='CH                  '
!!$      NAMES(STXCH2)='T-CH2               '
!!$      NAMES(SCH3)='CH3                 '
!!$      NAMES(SCH2O)='CH2O                '
!!$      NAMES(SHCCO)='HCCO                '
!!$      NAMES(SC2H)='C2H                 '
!!$      NAMES(SCH2CO)='CH2CO               '
!!$      NAMES(SC2H2)='C2H2                '
!!$      NAMES(SSXCH2)='S-CH2               '
!!$      NAMES(SAR)='AR                  '
!!$      NAMES(SCH3OH)='CH3OH               '
!!$      NAMES(SCH2OH)='CH2OH               '
!!$      NAMES(SCH3O)='CH3O                '
!!$      NAMES(SCH4)='CH4                 '
!!$      NAMES(SCH3O2)='CH3O2               '
!!$      NAMES(SC2H3)='C2H3                '
!!$      NAMES(SC2H4)='C2H4                '
!!$      NAMES(SC2H5)='C2H5                '
!!$      NAMES(SHCCOH)='HCCOH               '
!!$      NAMES(SCH2CHO)='CH2CHO              '
!!$      NAMES(SCH3CHO)='CH3CHO              '
!!$      NAMES(SH2C2)='H2C2                '
!!$      NAMES(SC2H5O)='C2H5O               '
!!$      NAMES(SNXC3H7)='N-C3H7              '
!!$      NAMES(SC2H6)='C2H6                '
!!$      NAMES(SC3H8)='C3H8                '
!!$      NAMES(SC3H6)='C3H6                '
!!$      NAMES(SC3H3)='C3H3                '
!!$      NAMES(SPXC3H4)='P-C3H4              '
!!$      NAMES(SAXC3H4)='A-C3H4              '
!!$      NAMES(SSXC3H5)='S-C3H5              '
!!$      NAMES(SNXC4H3)='N-C4H3              '
!!$      NAMES(SC2H3CHO)='C2H3CHO             '
!!$      NAMES(SAXC3H5)='A-C3H5              '
!!$      NAMES(SC2O)='C2O                 '
!!$      NAMES(SC4H4)='C4H4                '
!!$      NAMES(SC3H2)='C3H2                '
!!$      NAMES(SC3H2O)='C3H2O               '
!!$      NAMES(SC4H2)='C4H2                '
!!$      NAMES(SIXC4H3)='I-C4H3              '
!!$      NAMES(STXC3H5)='T-C3H5              '
!!$      NAMES(SC3H5O)='C3H5O               '
!!$      NAMES(SC4H)='C4H                 '
!!$      NAMES(SC8H2)='C8H2                '
!!$      NAMES(SC6H2)='C6H2                '
!!$      NAMES(SC4H6)='C4H6                '
!!$      NAMES(SNXC4H5)='N-C4H5              '
!!$      NAMES(SIXC4H5)='I-C4H5              '
!!$      NAMES(SA1XC6H6)='A1-C6H6             '
!!$      NAMES(SNXC7H16)='N-C7H16             '
!!$      NAMES(SC5H11)='C5H11               '
!!$      NAMES(SPXC4H9)='P-C4H9              '
!!$      NAMES(SC7H15)='C7H15               '
!!$      NAMES(SPXC4H8)='P-C4H8              '
!!$      NAMES(SC5H10)='C5H10               '
!!$      NAMES(SC7H14)='C7H14               '
!!$      NAMES(SC7H15O)='C7H15O              '
!!$      NAMES(SC3H7CHO)='C3H7CHO             '
!!$      NAMES(SC4H7)='C4H7                '
!!$      NAMES(SC7H13)='C7H13               '
!!$      NAMES(SC5H9)='C5H9                '
!!$      NAMES(SC4H7O)='C4H7O               '
!!$      NAMES(SNXC3H7O)='N-C3H7O             '
!!$      NAMES(SIXC8H18)='I-C8H18             '
!!$      NAMES(SYXC7H15)='Y-C7H15             '
!!$      NAMES(SIXC4H8)='I-C4H8              '
!!$      NAMES(SIXC3H7)='I-C3H7              '
!!$      NAMES(STXC4H9)='T-C4H9              '
!!$      NAMES(SCXC8H17)='C-C8H17             '
!!$      NAMES(SYXC7H14)='Y-C7H14             '
!!$      NAMES(SDXC8H17O)='D-C8H17O            '
!!$      NAMES(SCH3COCH3)='CH3COCH3            '
!!$      NAMES(SIXC4H7)='I-C4H7              '
!!$      NAMES(SXXC7H13)='X-C7H13             '
!!$      NAMES(SIXC3H5CH)='I-C3H5CHO           '
!!$      NAMES(STXC4H9O)='T-C4H9O             '
!!$      NAMES(SIXC4H7O)='I-C4H7O             '
!!$      NAMES(SC5H4CH2)='C5H4CH2             '
!!$      NAMES(SA1XXC6H5)='A1--C6H5            '
!!$      NAMES(SA1C2H2XC)='A1C2H2-C8H7         '
!!$      NAMES(SA1C2H3XC)='A1C2H3-C8H8         '
!!$      NAMES(SA1C2HXC8)='A1C2H-C8H6          '
!!$      NAMES(SA1C2HYXC)='A1C2H*-C8H5         '
!!$      NAMES(SA1C2H3YX)='A1C2H3*-C8H7        '
!!$      NAMES(SA2XXC10H)='A2--C10H7           '
!!$      NAMES(SA2XC10H8)='A2-C10H8            '
!!$      NAMES(SA2YXC10H)='A2*-C10H7           '
!!$      NAMES(SA2C2H2AX)='A2C2H2A-C12H9       '
!!$      NAMES(SA2C2H2BX)='A2C2H2B-C12H9       '
!!$      NAMES(SA2C2HAXC)='A2C2HA-C12H8        '
!!$      NAMES(SA2C2HBXC)='A2C2HB-C12H8        '
!!$      NAMES(SA2C2HAYX)='A2C2HA*-C12H7       '
!!$      NAMES(SA2C2HBYX)='A2C2HB*-C12H7       '
!!$      NAMES(SA2R5XC12)='A2R5-C12H8          '
!!$      NAMES(SA2R5XXC1)='A2R5--C12H7         '
!!$      NAMES(SA2R5C2H2)='A2R5C2H2-C14H9      '
!!$      NAMES(SA2R5C2HX)='A2R5C2H-C14H8       '
!!$      NAMES(SA2R5C2HY)='A2R5C2H*-C14H7      '
!!$      NAMES(SP2XC12H1)='P2-C12H10           '
!!$      NAMES(SP2XXC12H)='P2--C12H9           '
!!$      NAMES(SA3XXC14H)='A3--C14H9           '
!!$      NAMES(SA3XC14H1)='A3-C14H10           '
!!$      NAMES(SA3YXC14H)='A3*-C14H9           '
!!$      NAMES(SA3R5XXC1)='A3R5--C16H9         '
!!$      NAMES(SA3R5XC16)='A3R5-C16H10         '
!!$      NAMES(SA4XC16H1)='A4-C16H10           '
!!$      NAMES(SA4XXC16H)='A4--C16H9           '
!!$      NAMES(SA4R5XC18)='A4R5-C18H10         '
!!$      NAMES(SFLTNXC16)='FLTN-C16H10         '
!!$      NAMES(SC5H6)='C5H6                '
!!$      NAMES(SC5H5)='C5H5                '
!!$      NAMES(STXC5H5O)='T-C5H5O             '
!!$      NAMES(SC5H4O)='C5H4O               '
!!$      NAMES(SSXC5H5O)='S-C5H5O             '
!!$      NAMES(SC9H8)='C9H8                '
!!$      NAMES(SC9H7)='C9H7                '
!!$      NAMES(SA1CH2XC7)='A1CH2-C7H7          '
!!$      NAMES(SC9H6O)='C9H6O               '
!!$      NAMES(SOXC6H4)='O-C6H4              '
!!$      NAMES(SA1CH3XC7)='A1CH3-C7H8          '
!!$      NAMES(SA1OHXC6H)='A1OH-C6H6O          '
!!$      NAMES(SHOA1CH3X)='HOA1CH3-C7H8O       '
!!$      NAMES(SOA1CH3XC)='OA1CH3-C7H7O        '
!!$      NAMES(SA1CH2OXC)='A1CH2O-C7H7O        '
!!$      NAMES(SA1CH2OHX)='A1CH2OH-C7H8O       '
!!$      NAMES(SA1CHOXC7)='A1CHO-C7H6O         '
!!$      NAMES(SA1OXC6H5)='A1O-C6H5O           '
!!$      NAMES(SA1CH3YXC)='A1CH3*-C7H7         '
!!$      NAMES(SA1C2H4XC)='A1C2H4-C8H9         '
!!$      NAMES(SA1C2H5XC)='A1C2H5-C8H10        '
!!$      NAMES(SC8H9O2)='C8H9O2              '
!!$      NAMES(SC8H8OOH)='C8H8OOH             '
!!$      NAMES(SOC8H7OOH)='OC8H7OOH            '
!!$      NAMES(SA1CH3CH3)='A1CH3CH3-C8H10      '
!!$      NAMES(SA1CH3CH2)='A1CH3CH2-C8H9       '
!!$      NAMES(SA1CH3CHO)='A1CH3CHO-C8H8O      '
!!$      NAMES(SA2CH3XC1)='A2CH3-C11H10        '
!!$      NAMES(SA1CHOCH2)='A1CHOCH2-C8H7O      '
!!$      NAMES(SA1CHOCHO)='A1CHOCHO-C8H6O2     '
!!$      NAMES(SA2OHXC10)='A2OH-C10H8O         '
!!$      NAMES(SA2CH2XC1)='A2CH2-C11H9         '
!!$      NAMES(SA2CH2OXC)='A2CH2O-C11H9O       '
!!$      NAMES(SA2CHOXC1)='A2CHO-C11H8O        '
!!$      NAMES(SA2OXC10H)='A2O-C10H7O          '
!!$      NAMES(SOC6H4O)='OC6H4O              '


      NAMES(SN2)  = 'N2'
      NAMES(SH)  = 'H'
      NAMES(SO2)  = 'O2'
      NAMES(SO)  = 'O'
      NAMES(SOH)  = 'OH'
      NAMES(SH2)  = 'H2'
      NAMES(SH2O)  = 'H2O'
      NAMES(SCO2)  = 'CO2'
      NAMES(SHO2)  = 'HO2'
      NAMES(SH2O2)  = 'H2O2'
      NAMES(SCO)  = 'CO'
      NAMES(SHCO)  = 'HCO'
      NAMES(SC)  = 'C'
      NAMES(SCH)  = 'CH'
      NAMES(STXCH2)  = 'TXCH2'
      NAMES(SCH3)  = 'CH3'
      NAMES(SCH2O)  = 'CH2O'
      NAMES(SHCCO)  = 'HCCO'
      NAMES(SC2H)  = 'C2H'
      NAMES(SCH2CO)  = 'CH2CO'
      NAMES(SC2H2)  = 'C2H2'
      NAMES(SSXCH2)  = 'SXCH2'
      NAMES(SAR)  = 'AR'
      NAMES(SCH3OH)  = 'CH3OH'
      NAMES(SCH2OH)  = 'CH2OH'
      NAMES(SCH3O)  = 'CH3O'
      NAMES(SCH4)  = 'CH4'
      NAMES(SCH3O2)  = 'CH3O2'
      NAMES(SC2H3)  = 'C2H3'
      NAMES(SC2H4)  = 'C2H4'
      NAMES(SC2H5)  = 'C2H5'
      NAMES(SHCCOH)  = 'HCCOH'
      NAMES(SCH2CHO)  = 'CH2CHO'
      NAMES(SCH3CHO)  = 'CH3CHO'
      NAMES(SH2C2)  = 'H2C2'
      NAMES(SC2H5O)  = 'C2H5O'
      NAMES(SNXC3H7)  = 'NXC3H7'
      NAMES(SC2H6)  = 'C2H6'
      NAMES(SC3H8)  = 'C3H8'
      NAMES(SC3H6)  = 'C3H6'
      NAMES(SC3H3)  = 'C3H3'
      NAMES(SPXC3H4)  = 'PXC3H4'
      NAMES(SAXC3H4)  = 'AXC3H4'
      NAMES(SSXC3H5)  = 'SXC3H5'
      NAMES(SNXC4H3)  = 'NXC4H3'
      NAMES(SC2H3CHO)  = 'C2H3CHO'
      NAMES(SAXC3H5)  = 'AXC3H5'
      NAMES(SC2O)  = 'C2O'
      NAMES(SC4H4)  = 'C4H4'
      NAMES(SC3H2)  = 'C3H2'
      NAMES(SC3H2O)  = 'C3H2O'
      NAMES(SC4H2)  = 'C4H2'
      NAMES(SIXC4H3)  = 'IXC4H3'
      NAMES(STXC3H5)  = 'TXC3H5'
      NAMES(SC3H5O)  = 'C3H5O'
      NAMES(SC4H)  = 'C4H'
      NAMES(SC8H2)  = 'C8H2'
      NAMES(SC6H2)  = 'C6H2'
      NAMES(SC4H6)  = 'C4H6'
      NAMES(SNXC4H5)  = 'NXC4H5'
      NAMES(SIXC4H5)  = 'IXC4H5'
      NAMES(SA1XC6H6)  = 'A1XC6H6'
      NAMES(SNXC7H16)  = 'NXC7H16'
      NAMES(SC5H11)  = 'C5H11'
      NAMES(SPXC4H9)  = 'PXC4H9'
      NAMES(SC7H15)  = 'C7H15'
      NAMES(SPXC4H8)  = 'PXC4H8'
      NAMES(SC5H10)  = 'C5H10'
      NAMES(SC7H14)  = 'C7H14'
      NAMES(SC7H15O)  = 'C7H15O'
      NAMES(SC3H7CHO)  = 'C3H7CHO'
      NAMES(SC4H7)  = 'C4H7'
      NAMES(SC7H13)  = 'C7H13'
      NAMES(SC5H9)  = 'C5H9'
      NAMES(SC4H7O)  = 'C4H7O'
      NAMES(SNXC3H7O)  = 'NXC3H7O'
      NAMES(SIXC8H18)  = 'IXC8H18'
      NAMES(SYXC7H15)  = 'YXC7H15'
      NAMES(SIXC4H8)  = 'IXC4H8'
      NAMES(SIXC3H7)  = 'IXC3H7'
      NAMES(STXC4H9)  = 'TXC4H9'
      NAMES(SCXC8H17)  = 'CXC8H17'
      NAMES(SYXC7H14)  = 'YXC7H14'
      NAMES(SDXC8H17O)  = 'DXC8H17O'
      NAMES(SCH3COCH3)  = 'CH3COCH3'
      NAMES(SIXC4H7)  = 'IXC4H7'
      NAMES(SXXC7H13)  = 'XXC7H13'
      NAMES(SIXC3H5CH)  = 'IXC3H5CH'
      NAMES(STXC4H9O)  = 'TXC4H9O'
      NAMES(SIXC4H7O)  = 'IXC4H7O'
      NAMES(SC5H4CH2)  = 'C5H4CH2'
      NAMES(SA1XXC6H5)  = 'A1XXC6H5'
      NAMES(SA1C2H2XC)  = 'A1C2H2XC'
      NAMES(SA1C2H3XC)  = 'A1C2H3XC'
      NAMES(SA1C2HXC8)  = 'A1C2HXC8'
      NAMES(SA1C2HYXC)  = 'A1C2HYXC'
      NAMES(SA1C2H3YX)  = 'A1C2H3YX'
      NAMES(SA2XXC10H)  = 'A2XXC10H'
      NAMES(SA2XC10H8)  = 'A2XC10H8'
      NAMES(SA2YXC10H)  = 'A2YXC10H'
      NAMES(SA2C2H2AX)  = 'A2C2H2AX'
      NAMES(SA2C2H2BX)  = 'A2C2H2BX'
      NAMES(SA2C2HAXC)  = 'A2C2HAXC'
      NAMES(SA2C2HBXC)  = 'A2C2HBXC'
      NAMES(SA2C2HAYX)  = 'A2C2HAYX'
      NAMES(SA2C2HBYX)  = 'A2C2HBYX'
      NAMES(SA2R5XC12)  = 'A2R5XC12'
      NAMES(SA2R5XXC1)  = 'A2R5XXC1'
      NAMES(SA2R5C2H2)  = 'A2R5C2H2'
      NAMES(SA2R5C2HX)  = 'A2R5C2HX'
      NAMES(SA2R5C2HY)  = 'A2R5C2HY'
      NAMES(SP2XC12H1)  = 'P2XC12H1'
      NAMES(SP2XXC12H)  = 'P2XXC12H'
      NAMES(SA3XXC14H)  = 'A3XXC14H'
      NAMES(SA3XC14H1)  = 'A3XC14H1'
      NAMES(SA3YXC14H)  = 'A3YXC14H'
      NAMES(SA3R5XXC1)  = 'A3R5XXC1'
      NAMES(SA3R5XC16)  = 'A3R5XC16'
      NAMES(SA4XC16H1)  = 'A4XC16H1'
      NAMES(SA4XXC16H)  = 'A4XXC16H'
      NAMES(SA4R5XC18)  = 'A4R5XC18'
      NAMES(SFLTNXC16)  = 'FLTNXC16'
      NAMES(SC5H6)  = 'C5H6'
      NAMES(SC5H5)  = 'C5H5'
      NAMES(STXC5H5O)  = 'TXC5H5O'
      NAMES(SC5H4O)  = 'C5H4O'
      NAMES(SSXC5H5O)  = 'SXC5H5O'
      NAMES(SC9H8)  = 'C9H8'
      NAMES(SC9H7)  = 'C9H7'
      NAMES(SA1CH2XC7)  = 'A1CH2XC7'
      NAMES(SC9H6O)  = 'C9H6O'
      NAMES(SOXC6H4)  = 'OXC6H4'
      NAMES(SA1CH3XC7)  = 'A1CH3XC7'
      NAMES(SA1OHXC6H)  = 'A1OHXC6H'
      NAMES(SHOA1CH3X)  = 'HOA1CH3X'
      NAMES(SOA1CH3XC)  = 'OA1CH3XC'
      NAMES(SA1CH2OXC)  = 'A1CH2OXC'
      NAMES(SA1CH2OHX)  = 'A1CH2OHX'
      NAMES(SA1CHOXC7)  = 'A1CHOXC7'
      NAMES(SA1OXC6H5)  = 'A1OXC6H5'
      NAMES(SA1CH3YXC)  = 'A1CH3YXC'
      NAMES(SA1C2H4XC)  = 'A1C2H4XC'
      NAMES(SA1C2H5XC)  = 'A1C2H5XC'
      NAMES(SC8H9O2)  = 'C8H9O2'
      NAMES(SC8H8OOH)  = 'C8H8OOH'
      NAMES(SOC8H7OOH)  = 'OC8H7OOH'
      NAMES(SA1CH3CH3)  = 'A1CH3CH3'
      NAMES(SA1CH3CH2)  = 'A1CH3CH2'
      NAMES(SA1CH3CHO)  = 'A1CH3CHO'
      NAMES(SA2CH3XC1)  = 'A2CH3XC1'
      NAMES(SA1CHOCH2)  = 'A1CHOCH2'
      NAMES(SA1CHOCHO)  = 'A1CHOCHO'
      NAMES(SA2OHXC10)  = 'A2OHXC10'
      NAMES(SA2CH2XC1)  = 'A2CH2XC1'
      NAMES(SA2CH2OXC)  = 'A2CH2OXC'
      NAMES(SA2CHOXC1)  = 'A2CHOXC1'
      NAMES(SA2OXC10H)  = 'A2OXC10H'
      NAMES(SOC6H4O)  = 'OC6H4O'


      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'TheSootF90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 1806

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'TheSootF90.h'

      NSPECIES_NONS = 158
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'TheSootF90.h'
      real(DP) :: MUCOEFF(158)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SCO2) =  1.25056029e-06
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SH2O2) =  1.30193418e-06
      MUCOEFF(SCO) =  1.06039632e-06
      MUCOEFF(SHCO) =  1.11568663e-06
      MUCOEFF(SC) =  8.50486820e-07
      MUCOEFF(SCH) =  1.27351520e-06
      MUCOEFF(STXCH2) =  6.92304430e-07
      MUCOEFF(SCH3) =  7.16749611e-07
      MUCOEFF(SCH2O) =  1.13489904e-06
      MUCOEFF(SHCCO) =  2.73563116e-06
      MUCOEFF(SC2H) =  7.94406422e-07
      MUCOEFF(SCH2CO) =  1.09806251e-06
      MUCOEFF(SC2H2) =  8.10245830e-07
      MUCOEFF(SSXCH2) =  6.92304430e-07
      MUCOEFF(SAR) =  1.52144564e-06
      MUCOEFF(SCH3OH) =  1.14921582e-06
      MUCOEFF(SCH2OH) =  1.09210283e-06
      MUCOEFF(SCH3O) =  1.09210283e-06
      MUCOEFF(SCH4) =  7.61887935e-07
      MUCOEFF(SCH3O2) =  1.39234784e-06
      MUCOEFF(SC2H3) =  8.25781476e-07
      MUCOEFF(SC2H4) =  8.96560348e-07
      MUCOEFF(SC2H5) =  7.77507128e-07
      MUCOEFF(SHCCOH) =  1.09806251e-06
      MUCOEFF(SCH2CHO) =  1.11114998e-06
      MUCOEFF(SCH3CHO) =  1.12408509e-06
      MUCOEFF(SH2C2) =  8.10245830e-07
      MUCOEFF(SC2H5O) =  9.21331248e-07
      MUCOEFF(SNXC3H7) =  7.05924132e-07
      MUCOEFF(SC2H6) =  7.90876817e-07
      MUCOEFF(SC3H8) =  7.14133965e-07
      MUCOEFF(SC3H6) =  6.97617690e-07
      MUCOEFF(SC3H3) =  7.36234631e-07
      MUCOEFF(SPXC3H4) =  7.45675363e-07
      MUCOEFF(SAXC3H4) =  7.45675363e-07
      MUCOEFF(SSXC3H5) =  6.89211145e-07
      MUCOEFF(SNXC4H3) =  7.10878342e-07
      MUCOEFF(SC2H3CHO) =  8.13052585e-07
      MUCOEFF(SAXC3H5) =  6.89211145e-07
      MUCOEFF(SC2O) =  1.15237034e-06
      MUCOEFF(SC4H4) =  7.17860399e-07
      MUCOEFF(SC3H2) =  9.79454295e-07
      MUCOEFF(SC3H2O) =  8.66094458e-07
      MUCOEFF(SC4H2) =  7.03827024e-07
      MUCOEFF(SIXC4H3) =  7.10878342e-07
      MUCOEFF(STXC3H5) =  6.89211145e-07
      MUCOEFF(SC3H5O) =  8.67975059e-07
      MUCOEFF(SC4H) =  6.96704344e-07
      MUCOEFF(SC8H2) =  8.19457368e-07
      MUCOEFF(SC6H2) =  8.56202770e-07
      MUCOEFF(SC4H6) =  7.31624648e-07
      MUCOEFF(SNXC4H5) =  7.24775199e-07
      MUCOEFF(SIXC4H5) =  7.24775199e-07
      MUCOEFF(SA1XC6H6) =  8.43012086e-07
      MUCOEFF(SNXC7H16) =  6.83360789e-07
      MUCOEFF(SC5H11) =  8.85961433e-07
      MUCOEFF(SPXC4H9) =  7.34680473e-07
      MUCOEFF(SC7H15) =  6.79914768e-07
      MUCOEFF(SPXC4H8) =  7.72325002e-07
      MUCOEFF(SC5H10) =  7.41929850e-07
      MUCOEFF(SC7H14) =  6.94097963e-07
      MUCOEFF(SC7H15O) =  7.17931521e-07
      MUCOEFF(SC3H7CHO) =  9.03389424e-07
      MUCOEFF(SC4H7) =  9.16329178e-07
      MUCOEFF(SC7H13) =  6.90525741e-07
      MUCOEFF(SC5H9) =  7.44969464e-07
      MUCOEFF(SC4H7O) =  8.32363989e-07
      MUCOEFF(SNXC3H7O) =  8.21714561e-07
      MUCOEFF(SIXC8H18) =  6.93454798e-07
      MUCOEFF(SYXC7H15) =  6.98783435e-07
      MUCOEFF(SIXC4H8) =  7.72021504e-07
      MUCOEFF(SIXC3H7) =  7.57312844e-07
      MUCOEFF(STXC4H9) =  7.34680473e-07
      MUCOEFF(SCXC8H17) =  6.90388230e-07
      MUCOEFF(SYXC7H14) =  6.99071939e-07
      MUCOEFF(SDXC8H17O) =  7.16848556e-07
      MUCOEFF(SCH3COCH3) =  8.61252854e-07
      MUCOEFF(SIXC4H7) =  9.16329178e-07
      MUCOEFF(SXXC7H13) =  6.95474118e-07
      MUCOEFF(SIXC3H5CH) =  7.80165984e-07
      MUCOEFF(STXC4H9O) =  8.44082762e-07
      MUCOEFF(SIXC4H7O) =  8.32363989e-07
      MUCOEFF(SC5H4CH2) =  8.43012086e-07
      MUCOEFF(SA1XXC6H5) =  8.37554799e-07
      MUCOEFF(SA1C2H2XC) =  7.53008758e-07
      MUCOEFF(SA1C2H3XC) =  7.56679578e-07
      MUCOEFF(SA1C2HXC8) =  8.24475477e-07
      MUCOEFF(SA1C2HYXC) =  8.20396614e-07
      MUCOEFF(SA1C2H3YX) =  7.53008758e-07
      MUCOEFF(SA2XXC10H) =  7.88113678e-07
      MUCOEFF(SA2XC10H8) =  7.91231307e-07
      MUCOEFF(SA2YXC10H) =  7.88113678e-07
      MUCOEFF(SA2C2H2AX) =  7.89235966e-07
      MUCOEFF(SA2C2H2BX) =  7.89235966e-07
      MUCOEFF(SA2C2HAXC) =  7.86635103e-07
      MUCOEFF(SA2C2HBXC) =  7.86635103e-07
      MUCOEFF(SA2C2HAYX) =  7.84025612e-07
      MUCOEFF(SA2C2HBYX) =  7.84025612e-07
      MUCOEFF(SA2R5XC12) =  7.86635103e-07
      MUCOEFF(SA2R5XXC1) =  7.84025612e-07
      MUCOEFF(SA2R5C2H2) =  7.33542820e-07
      MUCOEFF(SA2R5C2HX) =  7.35675558e-07
      MUCOEFF(SA2R5C2HY) =  7.33568271e-07
      MUCOEFF(SP2XC12H1) =  8.32493507e-07
      MUCOEFF(SP2XXC12H) =  8.29768055e-07
      MUCOEFF(SA3XXC14H) =  7.33542820e-07
      MUCOEFF(SA3XC14H1) =  7.35626095e-07
      MUCOEFF(SA3YXC14H) =  7.33542820e-07
      MUCOEFF(SA3R5XXC1) =  7.14468605e-07
      MUCOEFF(SA3R5XC16) =  7.16255807e-07
      MUCOEFF(SA4XC16H1) =  7.24192099e-07
      MUCOEFF(SA4XXC16H) =  7.22385094e-07
      MUCOEFF(SA4R5XC18) =  7.65991844e-07
      MUCOEFF(SFLTNXC16) =  7.24192099e-07
      MUCOEFF(SC5H6) =  8.02573579e-07
      MUCOEFF(SC5H5) =  7.96430411e-07
      MUCOEFF(STXC5H5O) =  9.24146205e-07
      MUCOEFF(SC5H4O) =  9.18384382e-07
      MUCOEFF(SSXC5H5O) =  7.09635066e-07
      MUCOEFF(SC9H8) =  7.53247192e-07
      MUCOEFF(SC9H7) =  7.49971680e-07
      MUCOEFF(SA1CH2XC7) =  7.89808619e-07
      MUCOEFF(SC9H6O) =  7.97301352e-07
      MUCOEFF(SOXC6H4) =  8.32061720e-07
      MUCOEFF(SA1CH3XC7) =  7.94164881e-07
      MUCOEFF(SA1OHXC6H) =  7.38868658e-07
      MUCOEFF(SHOA1CH3X) =  8.60363246e-07
      MUCOEFF(SOA1CH3XC) =  8.56343804e-07
      MUCOEFF(SA1CH2OXC) =  8.56343804e-07
      MUCOEFF(SA1CH2OHX) =  8.60363246e-07
      MUCOEFF(SA1CHOXC7) =  8.52305406e-07
      MUCOEFF(SA1OXC6H5) =  7.34900958e-07
      MUCOEFF(SA1CH3YXC) =  7.89808619e-07
      MUCOEFF(SA1C2H4XC) =  7.60332675e-07
      MUCOEFF(SA1C2H5XC) =  7.63968304e-07
      MUCOEFF(SC8H9O2) =  8.68352299e-07
      MUCOEFF(SC8H8OOH) =  8.68352299e-07
      MUCOEFF(SOC8H7OOH) =  9.14581265e-07
      MUCOEFF(SA1CH3CH3) =  7.63968304e-07
      MUCOEFF(SA1CH3CH2) =  7.60332675e-07
      MUCOEFF(SA1CH3CHO) =  8.12729323e-07
      MUCOEFF(SA2CH3XC1) =  7.60367161e-07
      MUCOEFF(SA1CHOCH2) =  8.09312770e-07
      MUCOEFF(SA1CHOCHO) =  8.58725995e-07
      MUCOEFF(SA2OHXC10) =  8.39167872e-07
      MUCOEFF(SA2CH2XC1) =  7.57667206e-07
      MUCOEFF(SA2CH2OXC) =  7.99448018e-07
      MUCOEFF(SA2CHOXC1) =  7.96880486e-07
      MUCOEFF(SA2OXC10H) =  8.36228979e-07
      MUCOEFF(SOC6H4O) =  9.91705721e-07

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'TheSootF90.h'
      real(DP) :: KOVEREPS(158)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SCO2) =  4.09836066e-03
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SH2O2) =  9.31098696e-03
      KOVEREPS(SCO) =  1.01936799e-02
      KOVEREPS(SHCO) =  2.00803213e-03
      KOVEREPS(SC) =  1.40056022e-02
      KOVEREPS(SCH) =  1.25000000e-02
      KOVEREPS(STXCH2) =  6.94444444e-03
      KOVEREPS(SCH3) =  6.94444444e-03
      KOVEREPS(SCH2O) =  2.00803213e-03
      KOVEREPS(SHCCO) =  6.66666667e-03
      KOVEREPS(SC2H) =  4.78468900e-03
      KOVEREPS(SCH2CO) =  2.29357798e-03
      KOVEREPS(SC2H2) =  4.78468900e-03
      KOVEREPS(SSXCH2) =  6.94444444e-03
      KOVEREPS(SAR) =  7.32600733e-03
      KOVEREPS(SCH3OH) =  2.07555002e-03
      KOVEREPS(SCH2OH) =  2.39808153e-03
      KOVEREPS(SCH3O) =  2.39808153e-03
      KOVEREPS(SCH4) =  7.07213579e-03
      KOVEREPS(SCH3O2) =  2.07555002e-03
      KOVEREPS(SC2H3) =  4.78468900e-03
      KOVEREPS(SC2H4) =  3.56125356e-03
      KOVEREPS(SC2H5) =  3.96353547e-03
      KOVEREPS(SHCCOH) =  2.29357798e-03
      KOVEREPS(SCH2CHO) =  2.29357798e-03
      KOVEREPS(SCH3CHO) =  2.29357798e-03
      KOVEREPS(SH2C2) =  4.78468900e-03
      KOVEREPS(SC2H5O) =  2.12494688e-03
      KOVEREPS(SNXC3H7) =  3.74812594e-03
      KOVEREPS(SC2H6) =  3.96353547e-03
      KOVEREPS(SC3H8) =  3.74812594e-03
      KOVEREPS(SC3H6) =  3.74812594e-03
      KOVEREPS(SC3H3) =  3.96825397e-03
      KOVEREPS(SPXC3H4) =  3.96825397e-03
      KOVEREPS(SAXC3H4) =  3.96825397e-03
      KOVEREPS(SSXC3H5) =  3.74812594e-03
      KOVEREPS(SNXC4H3) =  2.80112045e-03
      KOVEREPS(SC2H3CHO) =  2.33208955e-03
      KOVEREPS(SAXC3H5) =  3.74812594e-03
      KOVEREPS(SC2O) =  4.30292599e-03
      KOVEREPS(SC4H4) =  2.80112045e-03
      KOVEREPS(SC3H2) =  4.78468900e-03
      KOVEREPS(SC3H2O) =  3.96825397e-03
      KOVEREPS(SC4H2) =  2.80112045e-03
      KOVEREPS(SIXC4H3) =  2.80112045e-03
      KOVEREPS(STXC3H5) =  3.74812594e-03
      KOVEREPS(SC3H5O) =  2.43309002e-03
      KOVEREPS(SC4H) =  2.80112045e-03
      KOVEREPS(SC8H2) =  2.01897840e-03
      KOVEREPS(SC6H2) =  2.80112045e-03
      KOVEREPS(SC4H6) =  2.80112045e-03
      KOVEREPS(SNXC4H5) =  2.80112045e-03
      KOVEREPS(SIXC4H5) =  2.80112045e-03
      KOVEREPS(SA1XC6H6) =  2.15146299e-03
      KOVEREPS(SNXC7H16) =  2.17580505e-03
      KOVEREPS(SC5H11) =  2.26893712e-03
      KOVEREPS(SPXC4H9) =  2.84090909e-03
      KOVEREPS(SC7H15) =  2.17580505e-03
      KOVEREPS(SPXC4H8) =  2.89268152e-03
      KOVEREPS(SC5H10) =  2.58933195e-03
      KOVEREPS(SC7H14) =  2.18435998e-03
      KOVEREPS(SC7H15O) =  1.78253119e-03
      KOVEREPS(SC3H7CHO) =  2.15424386e-03
      KOVEREPS(SC4H7) =  2.81690141e-03
      KOVEREPS(SC7H13) =  2.18435998e-03
      KOVEREPS(SC5H9) =  2.52016129e-03
      KOVEREPS(SC4H7O) =  2.01612903e-03
      KOVEREPS(SNXC3H7O) =  2.07684320e-03
      KOVEREPS(SIXC8H18) =  2.18102508e-03
      KOVEREPS(SYXC7H15) =  2.28675966e-03
      KOVEREPS(SIXC4H8) =  2.90275762e-03
      KOVEREPS(SIXC3H7) =  3.29597891e-03
      KOVEREPS(STXC4H9) =  2.84090909e-03
      KOVEREPS(SCXC8H17) =  2.18102508e-03
      KOVEREPS(SYXC7H14) =  2.27686703e-03
      KOVEREPS(SDXC8H17O) =  1.72028213e-03
      KOVEREPS(SCH3COCH3) =  2.29621125e-03
      KOVEREPS(SIXC4H7) =  2.81690141e-03
      KOVEREPS(SXXC7H13) =  2.27686703e-03
      KOVEREPS(SIXC3H5CH) =  2.29147571e-03
      KOVEREPS(STXC4H9O) =  2.01612903e-03
      KOVEREPS(SIXC4H7O) =  2.01612903e-03
      KOVEREPS(SC5H4CH2) =  2.15146299e-03
      KOVEREPS(SA1XXC6H5) =  2.15146299e-03
      KOVEREPS(SA1C2H2XC) =  1.83083120e-03
      KOVEREPS(SA1C2H3XC) =  1.83083120e-03
      KOVEREPS(SA1C2HXC8) =  1.86706497e-03
      KOVEREPS(SA1C2HYXC) =  1.86706497e-03
      KOVEREPS(SA1C2H3YX) =  1.83083120e-03
      KOVEREPS(SA2XXC10H) =  1.58629442e-03
      KOVEREPS(SA2XC10H8) =  1.58629442e-03
      KOVEREPS(SA2YXC10H) =  1.58629442e-03
      KOVEREPS(SA2C2H2AX) =  1.44279325e-03
      KOVEREPS(SA2C2H2BX) =  1.44279325e-03
      KOVEREPS(SA2C2HAXC) =  1.44279325e-03
      KOVEREPS(SA2C2HBXC) =  1.44279325e-03
      KOVEREPS(SA2C2HAYX) =  1.44279325e-03
      KOVEREPS(SA2C2HBYX) =  1.44279325e-03
      KOVEREPS(SA2R5XC12) =  1.44279325e-03
      KOVEREPS(SA2R5XXC1) =  1.44279325e-03
      KOVEREPS(SA2R5C2H2) =  1.29533679e-03
      KOVEREPS(SA2R5C2HX) =  1.29399586e-03
      KOVEREPS(SA2R5C2HY) =  1.29399586e-03
      KOVEREPS(SP2XC12H1) =  1.47819660e-03
      KOVEREPS(SP2XXC12H) =  1.47819660e-03
      KOVEREPS(SA3XXC14H) =  1.29533679e-03
      KOVEREPS(SA3XC14H1) =  1.29533679e-03
      KOVEREPS(SA3YXC14H) =  1.29533679e-03
      KOVEREPS(SA3R5XXC1) =  1.19402985e-03
      KOVEREPS(SA3R5XC16) =  1.19402985e-03
      KOVEREPS(SA4XC16H1) =  1.19774823e-03
      KOVEREPS(SA4XXC16H) =  1.19774823e-03
      KOVEREPS(SA4R5XC18) =  1.19774823e-03
      KOVEREPS(SFLTNXC16) =  1.19774823e-03
      KOVEREPS(SC5H6) =  2.50000000e-03
      KOVEREPS(SC5H5) =  2.50000000e-03
      KOVEREPS(STXC5H5O) =  2.06611570e-03
      KOVEREPS(SC5H4O) =  2.06611570e-03
      KOVEREPS(SSXC5H5O) =  1.62074554e-03
      KOVEREPS(SC9H8) =  1.58629442e-03
      KOVEREPS(SC9H7) =  1.58629442e-03
      KOVEREPS(SA1CH2XC7) =  2.01897840e-03
      KOVEREPS(SC9H6O) =  1.58629442e-03
      KOVEREPS(SOXC6H4) =  2.15146299e-03
      KOVEREPS(SA1CH3XC7) =  2.01897840e-03
      KOVEREPS(SA1OHXC6H) =  2.43902439e-03
      KOVEREPS(SHOA1CH3X) =  2.01897840e-03
      KOVEREPS(SOA1CH3XC) =  2.01897840e-03
      KOVEREPS(SA1CH2OXC) =  2.01897840e-03
      KOVEREPS(SA1CH2OHX) =  2.01897840e-03
      KOVEREPS(SA1CHOXC7) =  2.01897840e-03
      KOVEREPS(SA1OXC6H5) =  2.43902439e-03
      KOVEREPS(SA1CH3YXC) =  2.01897840e-03
      KOVEREPS(SA1C2H4XC) =  1.83083120e-03
      KOVEREPS(SA1C2H5XC) =  1.83083120e-03
      KOVEREPS(SC8H9O2) =  1.83083120e-03
      KOVEREPS(SC8H8OOH) =  1.83083120e-03
      KOVEREPS(SOC8H7OOH) =  1.83083120e-03
      KOVEREPS(SA1CH3CH3) =  1.83083120e-03
      KOVEREPS(SA1CH3CH2) =  1.83083120e-03
      KOVEREPS(SA1CH3CHO) =  1.83083120e-03
      KOVEREPS(SA2CH3XC1) =  1.44279325e-03
      KOVEREPS(SA1CHOCH2) =  1.83083120e-03
      KOVEREPS(SA1CHOCHO) =  1.83083120e-03
      KOVEREPS(SA2OHXC10) =  1.58629442e-03
      KOVEREPS(SA2CH2XC1) =  1.44279325e-03
      KOVEREPS(SA2CH2OXC) =  1.44279325e-03
      KOVEREPS(SA2CHOXC1) =  1.44279325e-03
      KOVEREPS(SA2OXC10H) =  1.58629442e-03
      KOVEREPS(SOC6H4O) =  2.15146299e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'TheSootF90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(1806), C(158), M(10)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 158
!------------------------------------------------------------------
      implicit none
      include 'TheSootF90.h'
      real(DP) :: H(158), CP(158), T

      IF (T.GT.1000.0_DP) THEN
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  2.92664000e+00_DP + T * (  7.43988400e-04_DP &
	   + T * ( -1.89492000e-07_DP + T * (  2.52425950e-11_DP &
	   + T * ( -1.35067020e-15_DP ) ) ) ) ) -9.22797700e+02_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    2.92664000e+00_DP + T * (  1.48797680e-03_DP &
	   + T * ( -5.68476000e-07_DP + T * (  1.00970380e-10_DP &
	   + T * ( -6.75335100e-15_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000001e+00_DP + T * ( -1.15421486e-11_DP &
	   + T * (  5.38539827e-15_DP + T * ( -1.18378809e-18_DP &
	   + T * (  9.96394714e-23_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000001e+00_DP + T * ( -2.30842973e-11_DP &
	   + T * (  1.61561948e-14_DP + T * ( -4.73515235e-18_DP &
	   + T * (  4.98197357e-22_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.28253784e+00_DP + T * (  7.41543770e-04_DP &
	   + T * ( -2.52655556e-07_DP + T * (  5.23676387e-11_DP &
	   + T * ( -4.33435588e-15_DP ) ) ) ) ) -1.08845772e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.28253784e+00_DP + T * (  1.48308754e-03_DP &
	   + T * ( -7.57966669e-07_DP + T * (  2.09470555e-10_DP &
	   + T * ( -2.16717794e-14_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.56942078e+00_DP + T * ( -4.29870569e-05_DP &
	   + T * (  1.39828196e-08_DP + T * ( -2.50444497e-12_DP &
	   + T * (  2.45667382e-16_DP ) ) ) ) ) +  2.92175791e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.56942078e+00_DP + T * ( -8.59741137e-05_DP &
	   + T * (  4.19484589e-08_DP + T * ( -1.00177799e-11_DP &
	   + T * (  1.22833691e-15_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  2.86472886e+00_DP + T * (  5.28252240e-04_DP &
	   + T * ( -8.63609193e-08_DP + T * (  7.63046685e-12_DP &
	   + T * ( -2.66391752e-16_DP ) ) ) ) ) +  3.71885774e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    2.86472886e+00_DP + T * (  1.05650448e-03_DP &
	   + T * ( -2.59082758e-07_DP + T * (  3.05218674e-11_DP &
	   + T * ( -1.33195876e-15_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.33727920e+00_DP + T * ( -2.47012365e-05_DP &
	   + T * (  1.66485593e-07_DP + T * ( -4.48915985e-11_DP &
	   + T * (  4.00510752e-15_DP ) ) ) ) ) -9.50158922e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.33727920e+00_DP + T * ( -4.94024731e-05_DP &
	   + T * (  4.99456778e-07_DP + T * ( -1.79566394e-10_DP &
	   + T * (  2.00255376e-14_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  3.03399249e+00_DP + T * (  1.08845902e-03_DP &
	   + T * ( -5.46908393e-08_DP + T * ( -2.42604967e-11_DP &
	   + T * (  3.36401984e-15_DP ) ) ) ) ) -3.00042971e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    3.03399249e+00_DP + T * (  2.17691804e-03_DP &
	   + T * ( -1.64072518e-07_DP + T * ( -9.70419870e-11_DP &
	   + T * (  1.68200992e-14_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  3.85746029e+00_DP + T * (  2.20718513e-03_DP &
	   + T * ( -7.38271347e-07_DP + T * (  1.30872547e-10_DP &
	   + T * ( -9.44168328e-15_DP ) ) ) ) ) -4.87591660e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    3.85746029e+00_DP + T * (  4.41437026e-03_DP &
	   + T * ( -2.21481404e-06_DP + T * (  5.23490188e-10_DP &
	   + T * ( -4.72084164e-14_DP ) ) ) ) )
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
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  2.71518561e+00_DP + T * (  1.03126372e-03_DP &
	   + T * ( -3.32941924e-07_DP + T * (  5.75132520e-11_DP &
	   + T * ( -4.07295432e-15_DP ) ) ) ) ) -1.41518724e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    2.71518561e+00_DP + T * (  2.06252743e-03_DP &
	   + T * ( -9.98825771e-07_DP + T * (  2.30053008e-10_DP &
	   + T * ( -2.03647716e-14_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  2.77217438e+00_DP + T * (  2.47847763e-03_DP &
	   + T * ( -8.28152043e-07_DP + T * (  1.47290445e-10_DP &
	   + T * ( -1.06701742e-14_DP ) ) ) ) ) +  4.01191815e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    2.77217438e+00_DP + T * (  4.95695526e-03_DP &
	   + T * ( -2.48445613e-06_DP + T * (  5.89161778e-10_DP &
	   + T * ( -5.33508711e-14_DP ) ) ) ) )
      H(SC) =  6.92284763e+02_DP * ( &
	   T * (  2.49266888e+00_DP + T * (  2.39944642e-05_DP &
	   + T * ( -2.41445007e-08_DP + T * (  9.35727573e-12_DP &
	   + T * ( -9.74555786e-16_DP ) ) ) ) ) +  8.54512953e+04_DP )
      CP(SC) =  6.92284763e+02_DP * ( &
	    2.49266888e+00_DP + T * (  4.79889284e-05_DP &
	   + T * ( -7.24335020e-08_DP + T * (  3.74291029e-11_DP &
	   + T * ( -4.87277893e-15_DP ) ) ) ) )
      H(SCH) =  6.38680289e+02_DP * ( &
	   T * (  2.87846473e+00_DP + T * (  4.85456840e-04_DP &
	   + T * (  4.81485517e-08_DP + T * ( -3.26719623e-11_DP &
	   + T * (  3.52158766e-15_DP ) ) ) ) ) +  7.10124364e+04_DP )
      CP(SCH) =  6.38680289e+02_DP * ( &
	    2.87846473e+00_DP + T * (  9.70913681e-04_DP &
	   + T * (  1.44445655e-07_DP + T * ( -1.30687849e-10_DP &
	   + T * (  1.76079383e-14_DP ) ) ) ) )
      H(STXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.87410113e+00_DP + T * (  1.82819646e-03_DP &
	   + T * ( -4.69648657e-07_DP + T * (  6.50448872e-11_DP &
	   + T * ( -3.75455134e-15_DP ) ) ) ) ) +  4.62636040e+04_DP )
      CP(STXCH2) =  5.92780550e+02_DP * ( &
	    2.87410113e+00_DP + T * (  3.65639292e-03_DP &
	   + T * ( -1.40894597e-06_DP + T * (  2.60179549e-10_DP &
	   + T * ( -1.87727567e-14_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  2.97812060e+00_DP + T * (  2.89892600e-03_DP &
	   + T * ( -6.58526667e-07_DP + T * (  7.68244750e-11_DP &
	   + T * ( -3.58348320e-15_DP ) ) ) ) ) +  1.65095130e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    2.97812060e+00_DP + T * (  5.79785200e-03_DP &
	   + T * ( -1.97558000e-06_DP + T * (  3.07297900e-10_DP &
	   + T * ( -1.79174160e-14_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  1.76069008e+00_DP + T * (  4.60000041e-03_DP &
	   + T * ( -1.47419604e-06_DP + T * (  2.51603030e-10_DP &
	   + T * ( -1.76771128e-14_DP ) ) ) ) ) -1.39958323e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    1.76069008e+00_DP + T * (  9.20000082e-03_DP &
	   + T * ( -4.42258813e-06_DP + T * (  1.00641212e-09_DP &
	   + T * ( -8.83855640e-14_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  5.62820580e+00_DP + T * (  2.04267005e-03_DP &
	   + T * ( -5.31151567e-07_DP + T * (  7.15651300e-11_DP &
	   + T * ( -3.88156640e-15_DP ) ) ) ) ) +  1.93272150e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    5.62820580e+00_DP + T * (  4.08534010e-03_DP &
	   + T * ( -1.59345470e-06_DP + T * (  2.86260520e-10_DP &
	   + T * ( -1.94078320e-14_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  3.16780652e+00_DP + T * (  2.37610951e-03_DP &
	   + T * ( -6.12623590e-07_DP + T * (  7.60475630e-11_DP &
	   + T * ( -3.54465540e-15_DP ) ) ) ) ) +  6.71210650e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    3.16780652e+00_DP + T * (  4.75221902e-03_DP &
	   + T * ( -1.83787077e-06_DP + T * (  3.04190252e-10_DP &
	   + T * ( -1.77232770e-14_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  4.51129732e+00_DP + T * (  4.50179872e-03_DP &
	   + T * ( -1.38979878e-06_DP + T * (  2.30836470e-10_DP &
	   + T * ( -1.58967640e-14_DP ) ) ) ) ) -7.55105311e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    4.51129732e+00_DP + T * (  9.00359745e-03_DP &
	   + T * ( -4.16939635e-06_DP + T * (  9.23345882e-10_DP &
	   + T * ( -7.94838201e-14_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  4.14756964e+00_DP + T * (  2.98083332e-03_DP &
	   + T * ( -7.90982840e-07_DP + T * (  1.16853043e-10_DP &
	   + T * ( -7.22470426e-15_DP ) ) ) ) ) +  2.59359992e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    4.14756964e+00_DP + T * (  5.96166664e-03_DP &
	   + T * ( -2.37294852e-06_DP + T * (  4.67412171e-10_DP &
	   + T * ( -3.61235213e-14_DP ) ) ) ) )
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.29203842e+00_DP + T * (  2.32794318e-03_DP &
	   + T * ( -6.70639823e-07_DP + T * (  1.04476500e-10_DP &
	   + T * ( -6.79432730e-15_DP ) ) ) ) ) +  5.09259997e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    2.29203842e+00_DP + T * (  4.65588637e-03_DP &
	   + T * ( -2.01191947e-06_DP + T * (  4.17906000e-10_DP &
	   + T * ( -3.39716365e-14_DP ) ) ) ) )
      H(SAR) =  2.08129068e+02_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -7.45375000e+02_DP )
      CP(SAR) =  2.08129068e+02_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OH) =  2.59482554e+02_DP * ( &
	   T * (  1.78970791e+00_DP + T * (  7.04691460e-03_DP &
	   + T * ( -2.12166945e-06_DP + T * (  3.45427713e-10_DP &
	   + T * ( -2.34120440e-14_DP ) ) ) ) ) -2.53748747e+04_DP )
      CP(SCH3OH) =  2.59482554e+02_DP * ( &
	    1.78970791e+00_DP + T * (  1.40938292e-02_DP &
	   + T * ( -6.36500835e-06_DP + T * (  1.38171085e-09_DP &
	   + T * ( -1.17060220e-13_DP ) ) ) ) )
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
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  1.65326226e+00_DP + T * (  5.01315495e-03_DP &
	   + T * ( -1.10553746e-06_DP + T * (  1.34120785e-10_DP &
	   + T * ( -6.29393516e-15_DP ) ) ) ) ) -1.00095936e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    1.65326226e+00_DP + T * (  1.00263099e-02_DP &
	   + T * ( -3.31661238e-06_DP + T * (  5.36483138e-10_DP &
	   + T * ( -3.14696758e-14_DP ) ) ) ) )
      H(SCH3O2) =  1.76772973e+02_DP * ( &
	   T * (  5.92505819e+00_DP + T * (  4.50097271e-03_DP &
	   + T * ( -1.08084770e-06_DP + T * (  1.31090679e-10_DP &
	   + T * ( -6.28526006e-15_DP ) ) ) ) ) -1.53258958e+03_DP )
      CP(SCH3O2) =  1.76772973e+02_DP * ( &
	    5.92505819e+00_DP + T * (  9.00194542e-03_DP &
	   + T * ( -3.24254309e-06_DP + T * (  5.24362718e-10_DP &
	   + T * ( -3.14263003e-14_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.01672400e+00_DP + T * (  5.16511460e-03_DP &
	   + T * ( -1.56027450e-06_DP + T * (  2.54408220e-10_DP &
	   + T * ( -1.72521408e-14_DP ) ) ) ) ) +  3.46128739e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.01672400e+00_DP + T * (  1.03302292e-02_DP &
	   + T * ( -4.68082349e-06_DP + T * (  1.01763288e-09_DP &
	   + T * ( -8.62607041e-14_DP ) ) ) ) )
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
      H(SHCCOH) =  1.97790941e+02_DP * ( &
	   T * (  5.92382910e+00_DP + T * (  3.39618000e-03_DP &
	   + T * ( -8.55285467e-07_DP + T * (  1.12469603e-10_DP &
	   + T * ( -5.98802020e-15_DP ) ) ) ) ) +  7.26462600e+03_DP )
      CP(SHCCOH) =  1.97790941e+02_DP * ( &
	    5.92382910e+00_DP + T * (  6.79236000e-03_DP &
	   + T * ( -2.56585640e-06_DP + T * (  4.49878410e-10_DP &
	   + T * ( -2.99401010e-14_DP ) ) ) ) )
      H(SCH2CHO) =  1.93159093e+02_DP * ( &
	   T * (  2.42606357e+00_DP + T * (  8.62000105e-03_DP &
	   + T * ( -3.25710706e-06_DP + T * (  6.66389180e-10_DP &
	   + T * ( -5.64240156e-14_DP ) ) ) ) ) +  8.33106990e+02_DP )
      CP(SCH2CHO) =  1.93159093e+02_DP * ( &
	    2.42606357e+00_DP + T * (  1.72400021e-02_DP &
	   + T * ( -9.77132119e-06_DP + T * (  2.66555672e-09_DP &
	   + T * ( -2.82120078e-13_DP ) ) ) ) )
      H(SCH3CHO) =  1.88739217e+02_DP * ( &
	   T * (  2.68543112e+00_DP + T * (  8.84011865e-03_DP &
	   + T * ( -2.88467580e-06_DP + T * (  5.09201472e-10_DP &
	   + T * ( -3.75261870e-14_DP ) ) ) ) ) -2.21653701e+04_DP )
      CP(SCH3CHO) =  1.88739217e+02_DP * ( &
	    2.68543112e+00_DP + T * (  1.76802373e-02_DP &
	   + T * ( -8.65402739e-06_DP + T * (  2.03680589e-09_DP &
	   + T * ( -1.87630935e-13_DP ) ) ) ) )
      H(SH2C2) =  3.19340144e+02_DP * ( &
	   T * (  4.27803400e+00_DP + T * (  2.37814020e-03_DP &
	   + T * ( -5.43366967e-07_DP + T * (  6.36570150e-11_DP &
	   + T * ( -2.97727580e-15_DP ) ) ) ) ) +  4.83166880e+04_DP )
      CP(SH2C2) =  3.19340144e+02_DP * ( &
	    4.27803400e+00_DP + T * (  4.75628040e-03_DP &
	   + T * ( -1.63010090e-06_DP + T * (  2.54628060e-10_DP &
	   + T * ( -1.48863790e-14_DP ) ) ) ) )
      H(SC2H5O) =  1.84517088e+02_DP * ( &
	   T * (  2.46262349e+00_DP + T * (  1.04751979e-02_DP &
	   + T * ( -3.13097250e-06_DP + T * (  3.91101567e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.83932658e+03_DP )
      CP(SC2H5O) =  1.84517088e+02_DP * ( &
	    2.46262349e+00_DP + T * (  2.09503959e-02_DP &
	   + T * ( -9.39291750e-06_DP + T * (  1.56440627e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SNXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  7.70404050e+00_DP + T * (  8.02077000e-03_DP &
	   + T * ( -1.76053223e-06_DP + T * (  1.90636007e-10_DP &
	   + T * ( -7.87069240e-15_DP ) ) ) ) ) +  8.29795310e+03_DP )
      CP(SNXC3H7) =  1.92970803e+02_DP * ( &
	    7.70404050e+00_DP + T * (  1.60415400e-02_DP &
	   + T * ( -5.28159670e-06_DP + T * (  7.62544030e-10_DP &
	   + T * ( -3.93534620e-14_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  1.07188150e+00_DP + T * (  1.08426339e-02_DP &
	   + T * ( -3.34186890e-06_DP + T * (  5.53530003e-10_DP &
	   + T * ( -3.80005780e-14_DP ) ) ) ) ) -1.14263932e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    1.07188150e+00_DP + T * (  2.16852677e-02_DP &
	   + T * ( -1.00256067e-05_DP + T * (  2.21412001e-09_DP &
	   + T * ( -1.90002890e-13_DP ) ) ) ) )
      H(SC3H8) =  1.88559441e+02_DP * ( &
	   T * (  7.53413680e+00_DP + T * (  9.43611950e-03_DP &
	   + T * ( -2.09061637e-06_DP + T * (  2.28689123e-10_DP &
	   + T * ( -9.56761380e-15_DP ) ) ) ) ) -1.64675160e+04_DP )
      CP(SC3H8) =  1.88559441e+02_DP * ( &
	    7.53413680e+00_DP + T * (  1.88722390e-02_DP &
	   + T * ( -6.27184910e-06_DP + T * (  9.14756490e-10_DP &
	   + T * ( -4.78380690e-14_DP ) ) ) ) )
      H(SC3H6) =  1.97593517e+02_DP * ( &
	   T * (  4.71697982e-01_DP + T * (  1.44756535e-02_DP &
	   + T * ( -5.22006063e-06_DP + T * (  1.02860800e-09_DP &
	   + T * ( -8.46150282e-14_DP ) ) ) ) ) +  1.12603387e+03_DP )
      CP(SC3H6) =  1.97593517e+02_DP * ( &
	    4.71697982e-01_DP + T * (  2.89513070e-02_DP &
	   + T * ( -1.56601819e-05_DP + T * (  4.11443199e-09_DP &
	   + T * ( -4.23075141e-13_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  6.14915291e+00_DP + T * (  4.67031583e-03_DP &
	   + T * ( -1.25018451e-06_DP + T * (  1.72539079e-10_DP &
	   + T * ( -9.21649988e-15_DP ) ) ) ) ) +  3.83854848e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    6.14915291e+00_DP + T * (  9.34063166e-03_DP &
	   + T * ( -3.75055354e-06_DP + T * (  6.90156316e-10_DP &
	   + T * ( -4.60824994e-14_DP ) ) ) ) )
      H(SPXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  2.81460543e+00_DP + T * (  9.27622480e-03_DP &
	   + T * ( -3.18342256e-06_DP + T * (  5.99878425e-10_DP &
	   + T * ( -4.74970514e-14_DP ) ) ) ) ) +  2.07010771e+04_DP )
      CP(SPXC3H4) =  2.07536818e+02_DP * ( &
	    2.81460543e+00_DP + T * (  1.85524496e-02_DP &
	   + T * ( -9.55026768e-06_DP + T * (  2.39951370e-09_DP &
	   + T * ( -2.37485257e-13_DP ) ) ) ) )
      H(SAXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  2.56128757e+00_DP + T * (  9.75400640e-03_DP &
	   + T * ( -3.46871220e-06_DP + T * (  6.75412933e-10_DP &
	   + T * ( -5.50148658e-14_DP ) ) ) ) ) +  2.13894289e+04_DP )
      CP(SAXC3H4) =  2.07536818e+02_DP * ( &
	    2.56128757e+00_DP + T * (  1.95080128e-02_DP &
	   + T * ( -1.04061366e-05_DP + T * (  2.70165173e-09_DP &
	   + T * ( -2.75074329e-13_DP ) ) ) ) )
      H(SSXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  2.02509360e+00_DP + T * (  1.17756625e-02_DP &
	   + T * ( -4.27515187e-06_DP + T * (  8.48948055e-10_DP &
	   + T * ( -7.03589448e-14_DP ) ) ) ) ) +  3.11812042e+04_DP )
      CP(SSXC3H5) =  2.02443146e+02_DP * ( &
	    2.02509360e+00_DP + T * (  2.35513249e-02_DP &
	   + T * ( -1.28254556e-05_DP + T * (  3.39579222e-09_DP &
	   + T * ( -3.51794724e-13_DP ) ) ) ) )
      H(SNXC4H3) =  1.62821949e+02_DP * ( &
	   T * (  7.25330164e+00_DP + T * (  5.97904230e-03_DP &
	   + T * ( -1.75571892e-06_DP + T * (  2.74954688e-10_DP &
	   + T * ( -1.76803350e-14_DP ) ) ) ) ) +  6.28977574e+04_DP )
      CP(SNXC4H3) =  1.62821949e+02_DP * ( &
	    7.25330164e+00_DP + T * (  1.19580846e-02_DP &
	   + T * ( -5.26715675e-06_DP + T * (  1.09981875e-09_DP &
	   + T * ( -8.84016751e-14_DP ) ) ) ) )
      H(SC2H3CHO) =  1.48306161e+02_DP * ( &
	   T * (  5.56154592e+00_DP + T * (  8.96479185e-03_DP &
	   + T * ( -2.67821586e-06_DP + T * (  3.30738438e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.29035886e+04_DP )
      CP(SC2H3CHO) =  1.48306161e+02_DP * ( &
	    5.56154592e+00_DP + T * (  1.79295837e-02_DP &
	   + T * ( -8.03464758e-06_DP + T * (  1.32295375e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SAXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  2.28794927e+00_DP + T * (  1.18200788e-02_DP &
	   + T * ( -4.26304833e-06_DP + T * (  8.42096350e-10_DP &
	   + T * ( -6.94898898e-14_DP ) ) ) ) ) +  1.83033514e+04_DP )
      CP(SAXC3H5) =  2.02443146e+02_DP * ( &
	    2.28794927e+00_DP + T * (  2.36401575e-02_DP &
	   + T * ( -1.27891450e-05_DP + T * (  3.36838540e-09_DP &
	   + T * ( -3.47449449e-13_DP ) ) ) ) )
      H(SC2O) =  2.07754623e+02_DP * ( &
	   T * (  5.42468378e+00_DP + T * (  9.26969725e-04_DP &
	   + T * ( -1.72644319e-07_DP + T * (  1.69411557e-11_DP &
	   + T * ( -7.06630474e-16_DP ) ) ) ) ) +  3.31537194e+04_DP )
      CP(SC2O) =  2.07754623e+02_DP * ( &
	    5.42468378e+00_DP + T * (  1.85393945e-03_DP &
	   + T * ( -5.17932956e-07_DP + T * (  6.77646230e-11_DP &
	   + T * ( -3.53315237e-15_DP ) ) ) ) )
      H(SC4H4) =  1.59670072e+02_DP * ( &
	   T * (  4.97237210e+00_DP + T * (  9.65699520e-03_DP &
	   + T * ( -3.27065503e-06_DP + T * (  6.07512635e-10_DP &
	   + T * ( -4.74199476e-14_DP ) ) ) ) ) +  3.30561454e+04_DP )
      CP(SC4H4) =  1.59670072e+02_DP * ( &
	    4.97237210e+00_DP + T * (  1.93139904e-02_DP &
	   + T * ( -9.81196508e-06_DP + T * (  2.43005054e-09_DP &
	   + T * ( -2.37099738e-13_DP ) ) ) ) )
      H(SC3H2) =  2.18533880e+02_DP * ( &
	   T * (  7.67920588e+00_DP + T * (  1.92780413e-03_DP &
	   + T * ( -2.74476656e-07_DP + T * ( -1.54527148e-11_DP &
	   + T * (  5.86228404e-15_DP ) ) ) ) ) +  6.29136323e+04_DP )
      CP(SC3H2) =  2.18533880e+02_DP * ( &
	    7.67920588e+00_DP + T * (  3.85560826e-03_DP &
	   + T * ( -8.23429967e-07_DP + T * ( -6.18108592e-11_DP &
	   + T * (  2.93114202e-14_DP ) ) ) ) )
      H(SC3H2O) =  1.53838212e+02_DP * ( &
	   T * (  5.51551710e+00_DP + T * (  6.01482820e-03_DP &
	   + T * ( -2.03019663e-06_DP + T * (  3.72165652e-10_DP &
	   + T * ( -2.85176948e-14_DP ) ) ) ) ) +  1.29567538e+04_DP )
      CP(SC3H2O) =  1.53838212e+02_DP * ( &
	    5.51551710e+00_DP + T * (  1.20296564e-02_DP &
	   + T * ( -6.09058988e-06_DP + T * (  1.48866261e-09_DP &
	   + T * ( -1.42588474e-13_DP ) ) ) ) )
      H(SC4H2) =  1.66100767e+02_DP * ( &
	   T * (  9.75839793e+00_DP + T * (  1.89436611e-03_DP &
	   + T * (  1.02047338e-07_DP + T * ( -1.58413756e-10_DP &
	   + T * (  2.25860064e-14_DP ) ) ) ) ) +  5.22698696e+04_DP )
      CP(SC4H2) =  1.66100767e+02_DP * ( &
	    9.75839793e+00_DP + T * (  3.78873223e-03_DP &
	   + T * (  3.06142015e-07_DP + T * ( -6.33655024e-10_DP &
	   + T * (  1.12930032e-13_DP ) ) ) ) )
      H(SIXC4H3) =  1.62821949e+02_DP * ( &
	   T * (  7.29283596e+00_DP + T * (  6.08324745e-03_DP &
	   + T * ( -1.83641769e-06_DP + T * (  2.98228375e-10_DP &
	   + T * ( -2.00986184e-14_DP ) ) ) ) ) +  5.71961011e+04_DP )
      CP(SIXC4H3) =  1.62821949e+02_DP * ( &
	    7.29283596e+00_DP + T * (  1.21664949e-02_DP &
	   + T * ( -5.50925306e-06_DP + T * (  1.19291350e-09_DP &
	   + T * ( -1.00493092e-13_DP ) ) ) ) )
      H(STXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  3.15893724e+00_DP + T * (  1.02324668e-02_DP &
	   + T * ( -3.36492707e-06_DP + T * (  6.02893455e-10_DP &
	   + T * ( -4.53070324e-14_DP ) ) ) ) ) +  2.87351148e+04_DP )
      CP(STXC3H5) =  2.02443146e+02_DP * ( &
	    3.15893724e+00_DP + T * (  2.04649335e-02_DP &
	   + T * ( -1.00947812e-05_DP + T * (  2.41157382e-09_DP &
	   + T * ( -2.26535162e-13_DP ) ) ) ) )
      H(SC3H5O) =  1.45686701e+02_DP * ( &
	   T * (  3.39074577e+00_DP + T * (  1.20650810e-02_DP &
	   + T * ( -3.78836313e-06_DP + T * (  4.94752345e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  9.00757452e+03_DP )
      CP(SC3H5O) =  1.45686701e+02_DP * ( &
	    3.39074577e+00_DP + T * (  2.41301620e-02_DP &
	   + T * ( -1.13650894e-05_DP + T * (  1.97900938e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC4H) =  1.69514353e+02_DP * ( &
	   T * (  7.44964925e+00_DP + T * (  2.61586587e-03_DP &
	   + T * ( -6.54467943e-07_DP + T * (  7.51669528e-11_DP &
	   + T * ( -2.32243092e-15_DP ) ) ) ) ) +  9.30052443e+04_DP )
      CP(SC4H) =  1.69514353e+02_DP * ( &
	    7.44964925e+00_DP + T * (  5.23173174e-03_DP &
	   + T * ( -1.96340383e-06_DP + T * (  3.00667811e-10_DP &
	   + T * ( -1.16121546e-14_DP ) ) ) ) )
      H(SC8H2) =  8.47571766e+01_DP * ( &
	   T * (  1.63586996e+01_DP + T * (  5.42962975e-03_DP &
	   + T * ( -1.30551599e-06_DP + T * (  1.58526758e-10_DP &
	   + T * ( -7.60826312e-15_DP ) ) ) ) ) +  1.02366984e+05_DP )
      CP(SC8H2) =  8.47571766e+01_DP * ( &
	    1.63586996e+01_DP + T * (  1.08592595e-02_DP &
	   + T * ( -3.91654796e-06_DP + T * (  6.34107033e-10_DP &
	   + T * ( -3.80413156e-14_DP ) ) ) ) )
      H(SC6H2) =  1.12240672e+02_DP * ( &
	   T * (  1.25328010e+01_DP + T * (  4.38831605e-03_DP &
	   + T * ( -1.04432053e-06_DP + T * (  1.25929550e-10_DP &
	   + T * ( -6.01438420e-15_DP ) ) ) ) ) +  7.97843380e+04_DP )
      CP(SC6H2) =  1.12240672e+02_DP * ( &
	    1.25328010e+01_DP + T * (  8.77663210e-03_DP &
	   + T * ( -3.13296160e-06_DP + T * (  5.03718200e-10_DP &
	   + T * ( -3.00719210e-14_DP ) ) ) ) )
      H(SC4H6) =  1.53718755e+02_DP * ( &
	   T * ( -8.99531092e+00_DP + T * (  3.00857534e-02_DP &
	   + T * ( -1.40019253e-05_DP + T * (  3.33325140e-09_DP &
	   + T * ( -3.14847380e-13_DP ) ) ) ) ) +  1.49296107e+04_DP )
      CP(SC4H6) =  1.53718755e+02_DP * ( &
	   -8.99531092e+00_DP + T * (  6.01715069e-02_DP &
	   + T * ( -4.20057758e-05_DP + T * (  1.33330056e-08_DP &
	   + T * ( -1.57423690e-12_DP ) ) ) ) )
      H(SNXC4H5) =  1.56637905e+02_DP * ( &
	   T * (  4.87674639e+00_DP + T * (  1.13767150e-02_DP &
	   + T * ( -3.92382327e-06_DP + T * (  7.38128637e-10_DP &
	   + T * ( -5.82913132e-14_DP ) ) ) ) ) +  4.11081097e+04_DP )
      CP(SNXC4H5) =  1.56637905e+02_DP * ( &
	    4.87674639e+00_DP + T * (  2.27534299e-02_DP &
	   + T * ( -1.17714698e-05_DP + T * (  2.95251455e-09_DP &
	   + T * ( -2.91456566e-13_DP ) ) ) ) )
      H(SIXC4H5) =  1.56637905e+02_DP * ( &
	   T * (  4.34643669e+00_DP + T * (  1.22880720e-02_DP &
	   + T * ( -4.36512283e-06_DP + T * (  8.47120312e-10_DP &
	   + T * ( -6.87039266e-14_DP ) ) ) ) ) +  3.58709780e+04_DP )
      CP(SIXC4H5) =  1.56637905e+02_DP * ( &
	    4.34643669e+00_DP + T * (  2.45761440e-02_DP &
	   + T * ( -1.30953685e-05_DP + T * (  3.38848125e-09_DP &
	   + T * ( -3.43519633e-13_DP ) ) ) ) )
      H(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   T * ( -2.06240612e-01_DP + T * (  2.32061220e-02_DP &
	   + T * ( -9.25511787e-06_DP + T * (  1.97227634e-09_DP &
	   + T * ( -1.72073052e-13_DP ) ) ) ) ) +  8.09883905e+03_DP )
      CP(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   -2.06240612e-01_DP + T * (  4.64122440e-02_DP &
	   + T * ( -2.77653536e-05_DP + T * (  7.88910537e-09_DP &
	   + T * ( -8.60365259e-13_DP ) ) ) ) )
      H(SNXC7H16) =  8.29791014e+01_DP * ( &
	   T * (  5.14079241e+00_DP + T * (  3.26539335e-02_DP &
	   + T * ( -9.82758747e-06_DP + T * (  1.23431681e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.72533890e+04_DP )
      CP(SNXC7H16) =  8.29791014e+01_DP * ( &
	    5.14079241e+00_DP + T * (  6.53078671e-02_DP &
	   + T * ( -2.94827624e-05_DP + T * (  4.93726726e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H11) =  1.16876212e+02_DP * ( &
	   T * (  4.88920629e+00_DP + T * (  2.11417269e-02_DP &
	   + T * ( -6.19477000e-06_DP + T * (  7.60311908e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  3.43475468e+03_DP )
      CP(SC5H11) =  1.16876212e+02_DP * ( &
	    4.88920629e+00_DP + T * (  4.22834537e-02_DP &
	   + T * ( -1.85843100e-05_DP + T * (  3.04124763e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SPXC4H9) =  1.45579563e+02_DP * ( &
	   T * (  3.81812330e+00_DP + T * (  1.70395244e-02_DP &
	   + T * ( -4.97118110e-06_DP + T * (  6.08013118e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  6.65901304e+03_DP )
      CP(SPXC4H9) =  1.45579563e+02_DP * ( &
	    3.81812330e+00_DP + T * (  3.40790489e-02_DP &
	   + T * ( -1.49135433e-05_DP + T * (  2.43205247e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC7H15) =  8.38223611e+01_DP * ( &
	   T * (  3.74721159e+00_DP + T * (  3.24672581e-02_DP &
	   + T * ( -1.00447008e-05_DP + T * (  1.29354536e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.37018357e+03_DP )
      CP(SC7H15) =  8.38223611e+01_DP * ( &
	    3.74721159e+00_DP + T * (  6.49345162e-02_DP &
	   + T * ( -3.01341025e-05_DP + T * (  5.17418142e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SPXC4H8) =  1.48195138e+02_DP * ( &
	   T * (  3.04470367e+00_DP + T * (  1.63725883e-02_DP &
	   + T * ( -4.84544123e-06_DP + T * (  5.99360043e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.52177534e+03_DP )
      CP(SPXC4H8) =  1.48195138e+02_DP * ( &
	    3.04470367e+00_DP + T * (  3.27451765e-02_DP &
	   + T * ( -1.45363237e-05_DP + T * (  2.39744017e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H10) =  1.18556110e+02_DP * ( &
	   T * (  3.98580522e+00_DP + T * (  2.06214993e-02_DP &
	   + T * ( -6.14634990e-06_DP + T * (  7.65388102e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -5.70112071e+03_DP )
      CP(SC5H10) =  1.18556110e+02_DP * ( &
	    3.98580522e+00_DP + T * (  4.12429986e-02_DP &
	   + T * ( -1.84390497e-05_DP + T * (  3.06155241e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC7H14) =  8.46829358e+01_DP * ( &
	   T * (  5.45858240e+00_DP + T * (  2.93078813e-02_DP &
	   + T * ( -8.76964720e-06_DP + T * (  1.09504683e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.33463299e+04_DP )
      CP(SC7H14) =  8.46829358e+01_DP * ( &
	    5.45858240e+00_DP + T * (  5.86157625e-02_DP &
	   + T * ( -2.63089416e-05_DP + T * (  4.38018732e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC7H15O) =  7.21793558e+01_DP * ( &
	   T * (  7.08994686e+00_DP + T * (  3.16234273e-02_DP &
	   + T * ( -9.48738963e-06_DP + T * (  1.18579623e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.96708300e+04_DP )
      CP(SC7H15O) =  7.21793558e+01_DP * ( &
	    7.08994686e+00_DP + T * (  6.32468545e-02_DP &
	   + T * ( -2.84621689e-05_DP + T * (  4.74318493e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC3H7CHO) =  1.15310385e+02_DP * ( &
	   T * (  3.99143562e+00_DP + T * (  1.76801297e-02_DP &
	   + T * ( -5.37657257e-06_DP + T * (  6.80257197e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.76351912e+04_DP )
      CP(SC3H7CHO) =  1.15310385e+02_DP * ( &
	    3.99143562e+00_DP + T * (  3.53602593e-02_DP &
	   + T * ( -1.61297177e-05_DP + T * (  2.72102879e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC4H7) =  1.50906418e+02_DP * ( &
	   T * (  8.49073768e+00_DP + T * (  9.55284870e-03_DP &
	   + T * ( -2.24790221e-06_DP + T * (  2.68358167e-10_DP &
	   + T * ( -1.27250367e-14_DP ) ) ) ) ) +  2.04659294e+04_DP )
      CP(SC4H7) =  1.50906418e+02_DP * ( &
	    8.49073768e+00_DP + T * (  1.91056974e-02_DP &
	   + T * ( -6.74370664e-06_DP + T * (  1.07343267e-09_DP &
	   + T * ( -6.36251837e-14_DP ) ) ) ) )
      H(SC7H13) =  8.55613642e+01_DP * ( &
	   T * (  5.78335940e+00_DP + T * (  2.80473858e-02_DP &
	   + T * ( -8.49517907e-06_DP + T * (  1.07237834e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  4.88852256e+03_DP )
      CP(SC7H13) =  8.55613642e+01_DP * ( &
	    5.78335940e+00_DP + T * (  5.60947716e-02_DP &
	   + T * ( -2.54855372e-05_DP + T * (  4.28951338e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H9) =  1.20285003e+02_DP * ( &
	   T * (  3.78447384e+00_DP + T * (  1.96499757e-02_DP &
	   + T * ( -5.93368333e-06_DP + T * (  7.46485712e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  1.12889042e+04_DP )
      CP(SC5H9) =  1.20285003e+02_DP * ( &
	    3.78447384e+00_DP + T * (  3.92999515e-02_DP &
	   + T * ( -1.78010500e-05_DP + T * (  2.98594285e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC4H7O) =  1.16945257e+02_DP * ( &
	   T * (  6.21920403e+00_DP + T * (  1.55186555e-02_DP &
	   + T * ( -4.91383277e-06_DP + T * (  6.44514742e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  2.91790666e+03_DP )
      CP(SC4H7O) =  1.16945257e+02_DP * ( &
	    6.21920403e+00_DP + T * (  3.10373110e-02_DP &
	   + T * ( -1.47414983e-05_DP + T * (  2.57805897e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SNXC3H7O) =  1.40715906e+02_DP * ( &
	   T * (  3.42509077e+00_DP + T * (  1.46569606e-02_DP &
	   + T * ( -4.38194527e-06_DP + T * (  5.47000915e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -7.01563135e+03_DP )
      CP(SNXC3H7O) =  1.40715906e+02_DP * ( &
	    3.42509077e+00_DP + T * (  2.93139212e-02_DP &
	   + T * ( -1.31458358e-05_DP + T * (  2.18800366e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SIXC8H18) =  7.27897815e+01_DP * ( &
	   T * (  8.11399812e+00_DP + T * (  3.59469073e-02_DP &
	   + T * ( -1.08572515e-05_DP + T * (  1.36609200e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.29583815e+04_DP )
      CP(SIXC8H18) =  7.27897815e+01_DP * ( &
	    8.11399812e+00_DP + T * (  7.18938145e-02_DP &
	   + T * ( -3.25717545e-05_DP + T * (  5.46436799e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SYXC7H15) =  8.38223611e+01_DP * ( &
	   T * (  2.64218915e+00_DP + T * (  3.30616507e-02_DP &
	   + T * ( -1.00564725e-05_DP + T * (  1.26896539e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -6.15969086e+03_DP )
      CP(SYXC7H15) =  8.38223611e+01_DP * ( &
	    2.64218915e+00_DP + T * (  6.61233015e-02_DP &
	   + T * ( -3.01694174e-05_DP + T * (  5.07586154e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SIXC4H8) =  1.48195138e+02_DP * ( &
	   T * (  2.86958571e+00_DP + T * (  1.64824603e-02_DP &
	   + T * ( -4.88104033e-06_DP + T * (  6.04074228e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -4.22675236e+03_DP )
      CP(SIXC4H8) =  1.48195138e+02_DP * ( &
	    2.86958571e+00_DP + T * (  3.29649207e-02_DP &
	   + T * ( -1.46431210e-05_DP + T * (  2.41629691e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SIXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  8.06336900e+00_DP + T * (  7.87244000e-03_DP &
	   + T * ( -1.72746400e-06_DP + T * (  1.86931125e-10_DP &
	   + T * ( -7.70884400e-15_DP ) ) ) ) ) +  5.31387100e+03_DP )
      CP(SIXC3H7) =  1.92970803e+02_DP * ( &
	    8.06336900e+00_DP + T * (  1.57448800e-02_DP &
	   + T * ( -5.18239200e-06_DP + T * (  7.47724500e-10_DP &
	   + T * ( -3.85442200e-14_DP ) ) ) ) )
      H(STXC4H9) =  1.45579563e+02_DP * ( &
	   T * ( -1.58631603e+00_DP + T * (  2.11173753e-02_DP &
	   + T * ( -6.44271080e-06_DP + T * (  8.14133765e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  4.56608047e+03_DP )
      CP(STXC4H9) =  1.45579563e+02_DP * ( &
	   -1.58631603e+00_DP + T * (  4.22347506e-02_DP &
	   + T * ( -1.93281324e-05_DP + T * (  3.25653506e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCXC8H17) =  7.34378533e+01_DP * ( &
	   T * (  5.37958758e+00_DP + T * (  3.63998835e-02_DP &
	   + T * ( -1.11696659e-05_DP + T * (  1.42412168e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.02203872e+04_DP )
      CP(SCXC8H17) =  7.34378533e+01_DP * ( &
	    5.37958758e+00_DP + T * (  7.27997671e-02_DP &
	   + T * ( -3.35089976e-05_DP + T * (  5.69648672e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SYXC7H14) =  8.46829358e+01_DP * ( &
	   T * (  5.51949920e+00_DP + T * (  2.91994107e-02_DP &
	   + T * ( -8.66953547e-06_DP + T * (  1.07295509e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.61420030e+04_DP )
      CP(SYXC7H14) =  8.46829358e+01_DP * ( &
	    5.51949920e+00_DP + T * (  5.83988213e-02_DP &
	   + T * ( -2.60086064e-05_DP + T * (  4.29182035e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SDXC8H17O) =  6.43445084e+01_DP * ( &
	   T * (  1.04791401e+01_DP + T * (  3.45965196e-02_DP &
	   + T * ( -1.04618950e-05_DP + T * (  1.31870540e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.54260262e+04_DP )
      CP(SDXC8H17O) =  6.43445084e+01_DP * ( &
	    1.04791401e+01_DP + T * (  6.91930393e-02_DP &
	   + T * ( -3.13856849e-05_DP + T * (  5.27482162e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3COCH3) =  1.43158167e+02_DP * ( &
	   T * (  7.29796974e+00_DP + T * (  8.78284565e-03_DP &
	   + T * ( -2.10559355e-06_DP + T * (  2.55063883e-10_DP &
	   + T * ( -1.22180718e-14_DP ) ) ) ) ) -2.95368927e+04_DP )
      CP(SCH3COCH3) =  1.43158167e+02_DP * ( &
	    7.29796974e+00_DP + T * (  1.75656913e-02_DP &
	   + T * ( -6.31678065e-06_DP + T * (  1.02025553e-09_DP &
	   + T * ( -6.10903592e-14_DP ) ) ) ) )
      H(SIXC4H7) =  1.50906418e+02_DP * ( &
	   T * (  4.86718299e+00_DP + T * (  1.37705580e-02_DP &
	   + T * ( -3.99386733e-06_DP + T * (  4.84570495e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  1.34120718e+04_DP )
      CP(SIXC4H7) =  1.50906418e+02_DP * ( &
	    4.86718299e+00_DP + T * (  2.75411161e-02_DP &
	   + T * ( -1.19816020e-05_DP + T * (  1.93828198e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SXXC7H13) =  8.55613642e+01_DP * ( &
	   T * (  4.05917970e+00_DP + T * (  2.92133917e-02_DP &
	   + T * ( -8.83685540e-06_DP + T * (  1.11062484e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  1.47082184e+02_DP )
      CP(SXXC7H13) =  8.55613642e+01_DP * ( &
	    4.05917970e+00_DP + T * (  5.84267834e-02_DP &
	   + T * ( -2.65105662e-05_DP + T * (  4.44249937e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SIXC3H5CH) =  1.18627154e+02_DP * ( &
	   T * (  7.19597854e+00_DP + T * (  1.24978145e-02_DP &
	   + T * ( -3.68171107e-06_DP + T * (  4.52301075e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.72891601e+04_DP )
      CP(SIXC3H5CH) =  1.18627154e+02_DP * ( &
	    7.19597854e+00_DP + T * (  2.49956291e-02_DP &
	   + T * ( -1.10451332e-05_DP + T * (  1.80920430e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(STXC4H9O) =  1.13720593e+02_DP * ( &
	   T * (  6.29676884e+00_DP + T * (  1.74474283e-02_DP &
	   + T * ( -5.10250397e-06_DP + T * (  6.23539407e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -1.53396915e+04_DP )
      CP(STXC4H9O) =  1.13720593e+02_DP * ( &
	    6.29676884e+00_DP + T * (  3.48948567e-02_DP &
	   + T * ( -1.53075119e-05_DP + T * (  2.49415763e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SIXC4H7O) =  1.16945257e+02_DP * ( &
	   T * (  4.69209202e+00_DP + T * (  1.59142293e-02_DP &
	   + T * ( -4.88631303e-06_DP + T * (  6.25543413e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  4.11802116e+03_DP )
      CP(SIXC4H7O) =  1.16945257e+02_DP * ( &
	    4.69209202e+00_DP + T * (  3.18284586e-02_DP &
	   + T * ( -1.46589391e-05_DP + T * (  2.50217365e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H4CH2) =  1.06446715e+02_DP * ( &
	   T * (  2.78194214e+00_DP + T * (  2.03161008e-02_DP &
	   + T * ( -7.83431090e-06_DP + T * (  1.62764004e-09_DP &
	   + T * ( -1.39361817e-13_DP ) ) ) ) ) +  2.43155607e+04_DP )
      CP(SC5H4CH2) =  1.06446715e+02_DP * ( &
	    2.78194214e+00_DP + T * (  4.06322016e-02_DP &
	   + T * ( -2.35029327e-05_DP + T * (  6.51056017e-09_DP &
	   + T * ( -6.96809087e-13_DP ) ) ) ) )
      H(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   T * (  1.38016336e+00_DP + T * (  2.02016004e-02_DP &
	   + T * ( -8.07502950e-06_DP + T * (  1.72180830e-09_DP &
	   + T * ( -1.50192160e-13_DP ) ) ) ) ) +  3.86973520e+04_DP )
      CP(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	    1.38016336e+00_DP + T * (  4.04032009e-02_DP &
	   + T * ( -2.42250885e-05_DP + T * (  6.88723321e-09_DP &
	   + T * ( -7.50960802e-13_DP ) ) ) ) )
      H(SA1C2H2XC) =  8.06153041e+01_DP * ( &
	   T * (  5.98044803e+00_DP + T * (  2.34715873e-02_DP &
	   + T * ( -8.91261580e-06_DP + T * (  1.82444488e-09_DP &
	   + T * ( -1.54218606e-13_DP ) ) ) ) ) +  4.32864172e+04_DP )
      CP(SA1C2H2XC) =  8.06153041e+01_DP * ( &
	    5.98044803e+00_DP + T * (  4.69431747e-02_DP &
	   + T * ( -2.67378474e-05_DP + T * (  7.29777950e-09_DP &
	   + T * ( -7.71093028e-13_DP ) ) ) ) )
      H(SA1C2H3XC) =  7.98350361e+01_DP * ( &
	   T * (  5.40554217e-01_DP + T * (  3.08651181e-02_DP &
	   + T * ( -1.24649102e-05_DP + T * (  2.67616467e-09_DP &
	   + T * ( -2.34609968e-13_DP ) ) ) ) ) +  1.50413170e+04_DP )
      CP(SA1C2H3XC) =  7.98350361e+01_DP * ( &
	    5.40554217e-01_DP + T * (  6.17302362e-02_DP &
	   + T * ( -3.73947305e-05_DP + T * (  1.07046587e-08_DP &
	   + T * ( -1.17304984e-12_DP ) ) ) ) )
      H(SA1C2HXC8) =  8.14109745e+01_DP * ( &
	   T * (  5.81520488e+00_DP + T * (  2.20436466e-02_DP &
	   + T * ( -8.40179527e-06_DP + T * (  1.72568807e-09_DP &
	   + T * ( -1.46275782e-13_DP ) ) ) ) ) +  3.30271906e+04_DP )
      CP(SA1C2HXC8) =  8.14109745e+01_DP * ( &
	    5.81520488e+00_DP + T * (  4.40872933e-02_DP &
	   + T * ( -2.52053858e-05_DP + T * (  6.90275228e-09_DP &
	   + T * ( -7.31378908e-13_DP ) ) ) ) )
      H(SA1C2HYXC) =  8.22225079e+01_DP * ( &
	   T * (  7.23812069e+00_DP + T * (  1.91906054e-02_DP &
	   + T * ( -7.29502437e-06_DP + T * (  1.49290312e-09_DP &
	   + T * ( -1.26070293e-13_DP ) ) ) ) ) +  6.49528135e+04_DP )
      CP(SA1C2HYXC) =  8.22225079e+01_DP * ( &
	    7.23812069e+00_DP + T * (  3.83812109e-02_DP &
	   + T * ( -2.18850731e-05_DP + T * (  5.97161247e-09_DP &
	   + T * ( -6.30351467e-13_DP ) ) ) ) )
      H(SA1C2H3YX) =  8.06153041e+01_DP * ( &
	   T * (  3.90114779e+00_DP + T * (  2.57947010e-02_DP &
	   + T * ( -1.01693507e-05_DP + T * (  2.13977724e-09_DP &
	   + T * ( -1.84609351e-13_DP ) ) ) ) ) +  4.59935428e+04_DP )
      CP(SA1C2H3YX) =  8.06153041e+01_DP * ( &
	    3.90114779e+00_DP + T * (  5.15894020e-02_DP &
	   + T * ( -3.05080522e-05_DP + T * (  8.55910896e-09_DP &
	   + T * ( -9.23046757e-13_DP ) ) ) ) )
      H(SA2XXC10H) =  6.53869263e+01_DP * ( &
	   T * (  3.22892303e+00_DP + T * (  3.15632243e-02_DP &
	   + T * ( -1.26860794e-05_DP + T * (  2.71135172e-09_DP &
	   + T * ( -2.36685024e-13_DP ) ) ) ) ) +  4.78400840e+04_DP )
      CP(SA2XXC10H) =  6.53869263e+01_DP * ( &
	    3.22892303e+00_DP + T * (  6.31264486e-02_DP &
	   + T * ( -3.80582381e-05_DP + T * (  1.08454069e-08_DP &
	   + T * ( -1.18342512e-12_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * (  1.76826275e+00_DP + T * (  3.44571753e-02_DP &
	   + T * ( -1.38107392e-05_DP + T * (  2.94785772e-09_DP &
	   + T * ( -2.57194122e-13_DP ) ) ) ) ) +  1.45412795e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	    1.76826275e+00_DP + T * (  6.89143506e-02_DP &
	   + T * ( -4.14322176e-05_DP + T * (  1.17914309e-08_DP &
	   + T * ( -1.28597061e-12_DP ) ) ) ) )
      H(SA2YXC10H) =  6.53869263e+01_DP * ( &
	   T * (  3.29950506e+00_DP + T * (  3.15066683e-02_DP &
	   + T * ( -1.26586694e-05_DP + T * (  2.70451890e-09_DP &
	   + T * ( -2.36015394e-13_DP ) ) ) ) ) +  4.76658373e+04_DP )
      CP(SA2YXC10H) =  6.53869263e+01_DP * ( &
	    3.29950506e+00_DP + T * (  6.30133365e-02_DP &
	   + T * ( -3.79760083e-05_DP + T * (  1.08180756e-08_DP &
	   + T * ( -1.18007697e-12_DP ) ) ) ) )
      H(SA2C2H2AX) =  5.42739830e+01_DP * ( &
	   T * (  8.53385239e+00_DP + T * (  3.43771399e-02_DP &
	   + T * ( -1.33916700e-05_DP + T * (  2.78702865e-09_DP &
	   + T * ( -2.38124486e-13_DP ) ) ) ) ) +  5.27583345e+04_DP )
      CP(SA2C2H2AX) =  5.42739830e+01_DP * ( &
	    8.53385239e+00_DP + T * (  6.87542797e-02_DP &
	   + T * ( -4.01750101e-05_DP + T * (  1.11481146e-08_DP &
	   + T * ( -1.19062243e-12_DP ) ) ) ) )
      H(SA2C2H2BX) =  5.42739830e+01_DP * ( &
	   T * (  7.59341132e+00_DP + T * (  3.51388460e-02_DP &
	   + T * ( -1.36912519e-05_DP + T * (  2.85023082e-09_DP &
	   + T * ( -2.43721726e-13_DP ) ) ) ) ) +  5.29378620e+04_DP )
      CP(SA2C2H2BX) =  5.42739830e+01_DP * ( &
	    7.59341132e+00_DP + T * (  7.02776920e-02_DP &
	   + T * ( -4.10737558e-05_DP + T * (  1.14009233e-08_DP &
	   + T * ( -1.21860863e-12_DP ) ) ) ) )
      H(SA2C2HAXC) =  5.46334700e+01_DP * ( &
	   T * (  7.55690939e+00_DP + T * (  3.35536637e-02_DP &
	   + T * ( -1.30965256e-05_DP + T * (  2.73699948e-09_DP &
	   + T * ( -2.34921860e-13_DP ) ) ) ) ) +  3.91372224e+04_DP )
      CP(SA2C2HAXC) =  5.46334700e+01_DP * ( &
	    7.55690939e+00_DP + T * (  6.71073273e-02_DP &
	   + T * ( -3.92895768e-05_DP + T * (  1.09479979e-08_DP &
	   + T * ( -1.17460930e-12_DP ) ) ) ) )
      H(SA2C2HBXC) =  5.46334700e+01_DP * ( &
	   T * (  7.63899557e+00_DP + T * (  3.34667928e-02_DP &
	   + T * ( -1.30500037e-05_DP + T * (  2.72499707e-09_DP &
	   + T * ( -2.33731998e-13_DP ) ) ) ) ) +  3.92947046e+04_DP )
      CP(SA2C2HBXC) =  5.46334700e+01_DP * ( &
	    7.63899557e+00_DP + T * (  6.69335855e-02_DP &
	   + T * ( -3.91500110e-05_DP + T * (  1.08999883e-08_DP &
	   + T * ( -1.16865999e-12_DP ) ) ) ) )
      H(SA2C2HAYX) =  5.49977510e+01_DP * ( &
	   T * (  8.88459555e+00_DP + T * (  3.07883585e-02_DP &
	   + T * ( -1.20329594e-05_DP + T * (  2.51481775e-09_DP &
	   + T * ( -2.15735058e-13_DP ) ) ) ) ) +  7.36259656e+04_DP )
      CP(SA2C2HAYX) =  5.49977510e+01_DP * ( &
	    8.88459555e+00_DP + T * (  6.15767170e-02_DP &
	   + T * ( -3.60988783e-05_DP + T * (  1.00592710e-08_DP &
	   + T * ( -1.07867529e-12_DP ) ) ) ) )
      H(SA2C2HBYX) =  5.49977510e+01_DP * ( &
	   T * (  8.88789581e+00_DP + T * (  3.07752080e-02_DP &
	   + T * ( -1.20249504e-05_DP + T * (  2.51275178e-09_DP &
	   + T * ( -2.15535646e-13_DP ) ) ) ) ) +  7.39803629e+04_DP )
      CP(SA2C2HBYX) =  5.49977510e+01_DP * ( &
	    8.88789581e+00_DP + T * (  6.15504161e-02_DP &
	   + T * ( -3.60748511e-05_DP + T * (  1.00510071e-08_DP &
	   + T * ( -1.07767823e-12_DP ) ) ) ) )
      H(SA2R5XC12) =  5.46334700e+01_DP * ( &
	   T * (  3.65432884e+00_DP + T * (  3.76323618e-02_DP &
	   + T * ( -1.51621650e-05_DP + T * (  3.24488352e-09_DP &
	   + T * ( -2.83461654e-13_DP ) ) ) ) ) +  2.65223472e+04_DP )
      CP(SA2R5XC12) =  5.46334700e+01_DP * ( &
	    3.65432884e+00_DP + T * (  7.52647236e-02_DP &
	   + T * ( -4.54864951e-05_DP + T * (  1.29795341e-08_DP &
	   + T * ( -1.41730827e-12_DP ) ) ) ) )
      H(SA2R5XXC1) =  5.49977510e+01_DP * ( &
	   T * (  4.90108932e+00_DP + T * (  3.49465809e-02_DP &
	   + T * ( -1.41408620e-05_DP + T * (  3.03365740e-09_DP &
	   + T * ( -2.65367420e-13_DP ) ) ) ) ) +  5.94391140e+04_DP )
      CP(SA2R5XXC1) =  5.49977510e+01_DP * ( &
	    4.90108932e+00_DP + T * (  6.98931618e-02_DP &
	   + T * ( -4.24225860e-05_DP + T * (  1.21346296e-08_DP &
	   + T * ( -1.32683710e-12_DP ) ) ) ) )
      H(SA2R5C2H2) =  4.69174774e+01_DP * ( &
	   T * (  7.80126948e+00_DP + T * (  4.07307673e-02_DP &
	   + T * ( -1.62550125e-05_DP + T * (  3.45333605e-09_DP &
	   + T * ( -3.00000482e-13_DP ) ) ) ) ) +  6.50870635e+04_DP )
      CP(SA2R5C2H2) =  4.69174774e+01_DP * ( &
	    7.80126948e+00_DP + T * (  8.14615345e-02_DP &
	   + T * ( -4.87650376e-05_DP + T * (  1.38133442e-08_DP &
	   + T * ( -1.50000241e-12_DP ) ) ) ) )
      H(SA2R5C2HX) =  4.71858755e+01_DP * ( &
	   T * (  9.29417050e+00_DP + T * (  3.68848112e-02_DP &
	   + T * ( -1.45277531e-05_DP + T * (  3.05383987e-09_DP &
	   + T * ( -2.63098636e-13_DP ) ) ) ) ) +  5.14110057e+04_DP )
      CP(SA2R5C2HX) =  4.71858755e+01_DP * ( &
	    9.29417050e+00_DP + T * (  7.37696223e-02_DP &
	   + T * ( -4.35832594e-05_DP + T * (  1.22153595e-08_DP &
	   + T * ( -1.31549318e-12_DP ) ) ) ) )
      H(SA2R5C2HY) =  4.74573620e+01_DP * ( &
	   T * (  1.06117453e+01_DP + T * (  3.41337752e-02_DP &
	   + T * ( -1.34734130e-05_DP + T * (  2.83431000e-09_DP &
	   + T * ( -2.44192732e-13_DP ) ) ) ) ) +  8.51154522e+04_DP )
      CP(SA2R5C2HY) =  4.74573620e+01_DP * ( &
	    1.06117453e+01_DP + T * (  6.82675505e-02_DP &
	   + T * ( -4.04202391e-05_DP + T * (  1.13372400e-08_DP &
	   + T * ( -1.22096366e-12_DP ) ) ) ) )
      H(SP2XC12H1) =  5.39191958e+01_DP * ( &
	   T * (  5.73686527e+00_DP + T * (  3.77329723e-02_DP &
	   + T * ( -1.46228290e-05_DP + T * (  3.04040273e-09_DP &
	   + T * ( -2.60024582e-13_DP ) ) ) ) ) +  1.66022411e+04_DP )
      CP(SP2XC12H1) =  5.39191958e+01_DP * ( &
	    5.73686527e+00_DP + T * (  7.54659445e-02_DP &
	   + T * ( -4.38684869e-05_DP + T * (  1.21616109e-08_DP &
	   + T * ( -1.30012291e-12_DP ) ) ) ) )
      H(SP2XXC12H) =  5.42739830e+01_DP * ( &
	   T * (  3.97670430e+00_DP + T * (  3.89716035e-02_DP &
	   + T * ( -1.58160723e-05_DP + T * (  3.40940205e-09_DP &
	   + T * ( -2.99715664e-13_DP ) ) ) ) ) +  5.02964162e+04_DP )
      CP(SP2XXC12H) =  5.42739830e+01_DP * ( &
	    3.97670430e+00_DP + T * (  7.79432070e-02_DP &
	   + T * ( -4.74482168e-05_DP + T * (  1.36376082e-08_DP &
	   + T * ( -1.49857832e-12_DP ) ) ) ) )
      H(SA3XXC14H) =  4.69174774e+01_DP * ( &
	   T * (  4.71264594e+00_DP + T * (  4.33615001e-02_DP &
	   + T * ( -1.75337211e-05_DP + T * (  3.76148872e-09_DP &
	   + T * ( -3.29127952e-13_DP ) ) ) ) ) +  5.37939644e+04_DP )
      CP(SA3XXC14H) =  4.69174774e+01_DP * ( &
	    4.71264594e+00_DP + T * (  8.67230002e-02_DP &
	   + T * ( -5.26011632e-05_DP + T * (  1.50459549e-08_DP &
	   + T * ( -1.64563976e-12_DP ) ) ) ) )
      H(SA3XC14H1) =  4.66521154e+01_DP * ( &
	   T * (  3.38725839e+00_DP + T * (  4.60942802e-02_DP &
	   + T * ( -1.85762331e-05_DP + T * (  3.97809798e-09_DP &
	   + T * ( -3.47768516e-13_DP ) ) ) ) ) +  1.91061794e+04_DP )
      CP(SA3XC14H1) =  4.66521154e+01_DP * ( &
	    3.38725839e+00_DP + T * (  9.21885604e-02_DP &
	   + T * ( -5.57286994e-05_DP + T * (  1.59123919e-08_DP &
	   + T * ( -1.73884258e-12_DP ) ) ) ) )
      H(SA3YXC14H) =  4.69174774e+01_DP * ( &
	   T * (  4.71264594e+00_DP + T * (  4.33615001e-02_DP &
	   + T * ( -1.75337211e-05_DP + T * (  3.76148872e-09_DP &
	   + T * ( -3.29127952e-13_DP ) ) ) ) ) +  5.37939644e+04_DP )
      CP(SA3YXC14H) =  4.69174774e+01_DP * ( &
	    4.71264594e+00_DP + T * (  8.67230002e-02_DP &
	   + T * ( -5.26011632e-05_DP + T * (  1.50459549e-08_DP &
	   + T * ( -1.64563976e-12_DP ) ) ) ) )
      H(SA3R5XXC1) =  4.13171861e+01_DP * ( &
	   T * (  6.47777347e+00_DP + T * (  4.66668858e-02_DP &
	   + T * ( -1.89542299e-05_DP + T * (  4.07616817e-09_DP &
	   + T * ( -3.57133286e-13_DP ) ) ) ) ) +  6.47489852e+04_DP )
      CP(SA3R5XXC1) =  4.13171861e+01_DP * ( &
	    6.47777347e+00_DP + T * (  9.33337716e-02_DP &
	   + T * ( -5.68626896e-05_DP + T * (  1.63046727e-08_DP &
	   + T * ( -1.78566643e-12_DP ) ) ) ) )
      H(SA3R5XC16) =  4.11112540e+01_DP * ( &
	   T * (  5.07024731e+00_DP + T * (  4.94958652e-02_DP &
	   + T * ( -2.00466753e-05_DP + T * (  4.30511463e-09_DP &
	   + T * ( -3.76952302e-13_DP ) ) ) ) ) +  3.06528296e+04_DP )
      CP(SA3R5XC16) =  4.11112540e+01_DP * ( &
	    5.07024731e+00_DP + T * (  9.89917305e-02_DP &
	   + T * ( -6.01400259e-05_DP + T * (  1.72204585e-08_DP &
	   + T * ( -1.88476151e-12_DP ) ) ) ) )
      H(SA4XC16H1) =  4.11112540e+01_DP * ( &
	   T * (  4.54060055e+00_DP + T * (  4.99057604e-02_DP &
	   + T * ( -2.02125434e-05_DP + T * (  4.33937548e-09_DP &
	   + T * ( -3.79804636e-13_DP ) ) ) ) ) +  2.12755890e+04_DP )
      CP(SA4XC16H1) =  4.11112540e+01_DP * ( &
	    4.54060055e+00_DP + T * (  9.98115207e-02_DP &
	   + T * ( -6.06376301e-05_DP + T * (  1.73575019e-08_DP &
	   + T * ( -1.89902318e-12_DP ) ) ) ) )
      H(SA4XXC16H) =  4.13171861e+01_DP * ( &
	   T * (  5.85098138e+00_DP + T * (  4.71615553e-02_DP &
	   + T * ( -1.91632717e-05_DP + T * (  4.12147133e-09_DP &
	   + T * ( -3.61083018e-13_DP ) ) ) ) ) +  5.89572568e+04_DP )
      CP(SA4XXC16H) =  4.13171861e+01_DP * ( &
	    5.85098138e+00_DP + T * (  9.43231105e-02_DP &
	   + T * ( -5.74898152e-05_DP + T * (  1.64858853e-08_DP &
	   + T * ( -1.80541509e-12_DP ) ) ) ) )
      H(SA4R5XC18) =  3.67468399e+01_DP * ( &
	   T * (  6.20190827e+00_DP + T * (  5.33281375e-02_DP &
	   + T * ( -2.16943100e-05_DP + T * (  4.66940325e-09_DP &
	   + T * ( -4.09294700e-13_DP ) ) ) ) ) +  3.34439422e+04_DP )
      CP(SA4R5XC18) =  3.67468399e+01_DP * ( &
	    6.20190827e+00_DP + T * (  1.06656275e-01_DP &
	   + T * ( -6.50829301e-05_DP + T * (  1.86776130e-08_DP &
	   + T * ( -2.04647350e-12_DP ) ) ) ) )
      H(SFLTNXC16) =  4.11112540e+01_DP * ( &
	   T * (  4.54792547e+00_DP + T * (  4.99994870e-02_DP &
	   + T * ( -2.02900657e-05_DP + T * (  4.36322555e-09_DP &
	   + T * ( -3.82399156e-13_DP ) ) ) ) ) +  2.54780117e+04_DP )
      CP(SFLTNXC16) =  4.11112540e+01_DP * ( &
	    4.54792547e+00_DP + T * (  9.99989740e-02_DP &
	   + T * ( -6.08701972e-05_DP + T * (  1.74529022e-08_DP &
	   + T * ( -1.91199578e-12_DP ) ) ) ) )
      H(SC5H6) =  1.25788072e+02_DP * ( &
	   T * (  2.30537462e-01_DP + T * (  2.04785913e-02_DP &
	   + T * ( -8.05296527e-06_DP + T * (  1.69940870e-09_DP &
	   + T * ( -1.47274884e-13_DP ) ) ) ) ) +  1.43779465e+04_DP )
      CP(SC5H6) =  1.25788072e+02_DP * ( &
	    2.30537462e-01_DP + T * (  4.09571826e-02_DP &
	   + T * ( -2.41588958e-05_DP + T * (  6.79763480e-09_DP &
	   + T * ( -7.36374421e-13_DP ) ) ) ) )
      H(SC5H5) =  1.27736058e+02_DP * ( &
	   T * (  4.21464919e+00_DP + T * (  1.35917364e-02_DP &
	   + T * ( -4.43910697e-06_DP + T * (  7.72450297e-10_DP &
	   + T * ( -5.55759746e-14_DP ) ) ) ) ) +  2.88952416e+04_DP )
      CP(SC5H5) =  1.27736058e+02_DP * ( &
	    4.21464919e+00_DP + T * (  2.71834728e-02_DP &
	   + T * ( -1.33173209e-05_DP + T * (  3.08980119e-09_DP &
	   + T * ( -2.77879873e-13_DP ) ) ) ) )
      H(STXC5H5O) =  1.02532248e+02_DP * ( &
	   T * (  1.26065350e+01_DP + T * (  8.37353350e-03_DP &
	   + T * ( -2.03658623e-06_DP + T * (  2.49186440e-10_DP &
	   + T * ( -1.20223668e-14_DP ) ) ) ) ) +  1.41146570e+03_DP )
      CP(STXC5H5O) =  1.02532248e+02_DP * ( &
	    1.26065350e+01_DP + T * (  1.67470670e-02_DP &
	   + T * ( -6.10975870e-06_DP + T * (  9.96745760e-10_DP &
	   + T * ( -6.01118340e-14_DP ) ) ) ) )
      H(SC5H4O) =  1.03822832e+02_DP * ( &
	   T * (  4.25344911e+00_DP + T * (  1.56819909e-02_DP &
	   + T * ( -6.09546950e-06_DP + T * (  1.27102091e-09_DP &
	   + T * ( -1.08969098e-13_DP ) ) ) ) ) +  3.87579835e+03_DP )
      CP(SC5H4O) =  1.03822832e+02_DP * ( &
	    4.25344911e+00_DP + T * (  3.13639818e-02_DP &
	   + T * ( -1.82864085e-05_DP + T * (  5.08408365e-09_DP &
	   + T * ( -5.44845492e-13_DP ) ) ) ) )
      H(SSXC5H5O) =  1.02532248e+02_DP * ( &
	   T * (  8.54053120e+00_DP + T * (  1.14947550e-02_DP &
	   + T * ( -3.18125210e-06_DP + T * (  4.26540300e-10_DP &
	   + T * ( -1.94918720e-14_DP ) ) ) ) ) +  2.22636990e+04_DP )
      CP(SSXC5H5O) =  1.02532248e+02_DP * ( &
	    8.54053120e+00_DP + T * (  2.29895100e-02_DP &
	   + T * ( -9.54375630e-06_DP + T * (  1.70616120e-09_DP &
	   + T * ( -9.74593600e-14_DP ) ) ) ) )
      H(SC9H8) =  7.15803158e+01_DP * ( &
	   T * (  1.15459802e+00_DP + T * (  3.27112098e-02_DP &
	   + T * ( -1.30835036e-05_DP + T * (  2.78922102e-09_DP &
	   + T * ( -2.43185350e-13_DP ) ) ) ) ) +  1.68166108e+04_DP )
      CP(SC9H8) =  7.15803158e+01_DP * ( &
	    1.15459802e+00_DP + T * (  6.54224196e-02_DP &
	   + T * ( -3.92505107e-05_DP + T * (  1.11568841e-08_DP &
	   + T * ( -1.21592675e-12_DP ) ) ) ) )
      H(SC9H7) =  7.22069373e+01_DP * ( &
	   T * (  3.65597547e+00_DP + T * (  2.87404231e-02_DP &
	   + T * ( -1.14290200e-05_DP + T * (  2.42569698e-09_DP &
	   + T * ( -2.10772824e-13_DP ) ) ) ) ) +  3.06843457e+04_DP )
      CP(SC9H7) =  7.22069373e+01_DP * ( &
	    3.65597547e+00_DP + T * (  5.74808463e-02_DP &
	   + T * ( -3.42870600e-05_DP + T * (  9.70278793e-09_DP &
	   + T * ( -1.05386412e-12_DP ) ) ) ) )
      H(SA1CH2XC7) =  9.12400413e+01_DP * ( &
	   T * (  3.30049696e+00_DP + T * (  2.40027670e-02_DP &
	   + T * ( -9.28143407e-06_DP + T * (  1.93092839e-09_DP &
	   + T * ( -1.65430827e-13_DP ) ) ) ) ) +  2.17498572e+04_DP )
      CP(SA1CH2XC7) =  9.12400413e+01_DP * ( &
	    3.30049696e+00_DP + T * (  4.80055340e-02_DP &
	   + T * ( -2.78443022e-05_DP + T * (  7.72371356e-09_DP &
	   + T * ( -8.27154136e-13_DP ) ) ) ) )
      H(SC9H6O) =  6.38886413e+01_DP * ( &
	   T * (  4.65659248e+00_DP + T * (  2.85027911e-02_DP &
	   + T * ( -1.14391400e-05_DP + T * (  2.44044361e-09_DP &
	   + T * ( -2.12668074e-13_DP ) ) ) ) ) +  4.57857140e+03_DP )
      CP(SC9H6O) =  6.38886413e+01_DP * ( &
	    4.65659248e+00_DP + T * (  5.70055822e-02_DP &
	   + T * ( -3.43174199e-05_DP + T * (  9.76177442e-09_DP &
	   + T * ( -1.06334037e-12_DP ) ) ) ) )
      H(SOXC6H4) =  1.09266940e+02_DP * ( &
	   T * (  2.98618725e+00_DP + T * (  1.68818922e-02_DP &
	   + T * ( -6.67461303e-06_DP + T * (  1.40963421e-09_DP &
	   + T * ( -1.22000829e-13_DP ) ) ) ) ) +  5.12231321e+04_DP )
      CP(SOXC6H4) =  1.09266940e+02_DP * ( &
	    2.98618725e+00_DP + T * (  3.37637843e-02_DP &
	   + T * ( -2.00238391e-05_DP + T * (  5.63853682e-09_DP &
	   + T * ( -6.10004145e-13_DP ) ) ) ) )
      H(SA1CH3XC7) =  9.02418217e+01_DP * ( &
	   T * ( -1.01117220e+00_DP + T * (  2.92650956e-02_DP &
	   + T * ( -1.15865023e-05_DP + T * (  2.45545248e-09_DP &
	   + T * ( -2.13361740e-13_DP ) ) ) ) ) +  3.99363395e+03_DP )
      CP(SA1CH3XC7) =  9.02418217e+01_DP * ( &
	   -1.01117220e+00_DP + T * (  5.85301912e-02_DP &
	   + T * ( -3.47595069e-05_DP + T * (  9.82180993e-09_DP &
	   + T * ( -1.06680870e-12_DP ) ) ) ) )
      H(SA1OHXC6H) =  8.83489183e+01_DP * ( &
	   T * (  9.33151850e-01_DP + T * (  2.53298752e-02_DP &
	   + T * ( -1.05872946e-05_DP + T * (  2.34389723e-09_DP &
	   + T * ( -2.10588162e-13_DP ) ) ) ) ) -1.37575260e+04_DP )
      CP(SA1OHXC6H) =  8.83489183e+01_DP * ( &
	    9.33151850e-01_DP + T * (  5.06597504e-02_DP &
	   + T * ( -3.17618838e-05_DP + T * (  9.37558892e-09_DP &
	   + T * ( -1.05294081e-12_DP ) ) ) ) )
      H(SHOA1CH3X) =  7.68892300e+01_DP * ( &
	   T * ( -6.71722270e-02_DP + T * (  3.15933480e-02_DP &
	   + T * ( -1.30255983e-05_DP + T * (  2.85416352e-09_DP &
	   + T * ( -2.54530110e-13_DP ) ) ) ) ) -1.78554350e+04_DP )
      CP(SHOA1CH3X) =  7.68892300e+01_DP * ( &
	   -6.71722270e-02_DP + T * (  6.31866960e-02_DP &
	   + T * ( -3.90767950e-05_DP + T * (  1.14166541e-08_DP &
	   + T * ( -1.27265055e-12_DP ) ) ) ) )
      H(SOA1CH3XC) =  7.76127177e+01_DP * ( &
	   T * (  4.14521668e+00_DP + T * (  2.51392063e-02_DP &
	   + T * ( -9.69283660e-06_DP + T * (  2.00710898e-09_DP &
	   + T * ( -1.71158850e-13_DP ) ) ) ) ) -1.22524065e+03_DP )
      CP(SOA1CH3XC) =  7.76127177e+01_DP * ( &
	    4.14521668e+00_DP + T * (  5.02784126e-02_DP &
	   + T * ( -2.90785098e-05_DP + T * (  8.02843592e-09_DP &
	   + T * ( -8.55794251e-13_DP ) ) ) ) )
      H(SA1CH2OXC) =  7.76127177e+01_DP * ( &
	   T * (  2.07930551e+00_DP + T * (  2.76359457e-02_DP &
	   + T * ( -1.10434375e-05_DP + T * (  2.35169868e-09_DP &
	   + T * ( -2.04829704e-13_DP ) ) ) ) ) +  1.17599254e+04_DP )
      CP(SA1CH2OXC) =  7.76127177e+01_DP * ( &
	    2.07930551e+00_DP + T * (  5.52718914e-02_DP &
	   + T * ( -3.31303125e-05_DP + T * (  9.40679473e-09_DP &
	   + T * ( -1.02414852e-12_DP ) ) ) ) )
      H(SA1CH2OHX) =  7.68892300e+01_DP * ( &
	   T * (  1.51623145e+01_DP + T * (  1.34685184e-02_DP &
	   + T * ( -3.26763210e-06_DP + T * (  3.98205190e-10_DP &
	   + T * ( -1.91421567e-14_DP ) ) ) ) ) -1.88226234e+04_DP )
      CP(SA1CH2OHX) =  7.68892300e+01_DP * ( &
	    1.51623145e+01_DP + T * (  2.69370369e-02_DP &
	   + T * ( -9.80289631e-06_DP + T * (  1.59282076e-09_DP &
	   + T * ( -9.57107837e-14_DP ) ) ) ) )
      H(SA1CHOXC7) =  7.83499501e+01_DP * ( &
	   T * (  1.87355756e+00_DP + T * (  2.63115775e-02_DP &
	   + T * ( -1.05881654e-05_DP + T * (  2.26600767e-09_DP &
	   + T * ( -1.98061225e-13_DP ) ) ) ) ) -7.23603865e+03_DP )
      CP(SA1CHOXC7) =  7.83499501e+01_DP * ( &
	    1.87355756e+00_DP + T * (  5.26231551e-02_DP &
	   + T * ( -3.17644962e-05_DP + T * (  9.06403069e-09_DP &
	   + T * ( -9.90306123e-13_DP ) ) ) ) )
      H(SA1OXC6H5) =  8.93054780e+01_DP * ( &
	   T * (  3.39256520e+00_DP + T * (  2.08689845e-02_DP &
	   + T * ( -8.32789277e-06_DP + T * (  1.77206751e-09_DP &
	   + T * ( -1.54268601e-13_DP ) ) ) ) ) +  3.60336039e+03_DP )
      CP(SA1OXC6H5) =  8.93054780e+01_DP * ( &
	    3.39256520e+00_DP + T * (  4.17379690e-02_DP &
	   + T * ( -2.49836783e-05_DP + T * (  7.08827005e-09_DP &
	   + T * ( -7.71343006e-13_DP ) ) ) ) )
      H(SA1CH3YXC) =  9.12400413e+01_DP * ( &
	   T * (  5.19780822e-01_DP + T * (  2.63109877e-02_DP &
	   + T * ( -1.04327811e-05_DP + T * (  2.21204344e-09_DP &
	   + T * ( -1.92250440e-13_DP ) ) ) ) ) +  3.45681510e+04_DP )
      CP(SA1CH3YXC) =  9.12400413e+01_DP * ( &
	    5.19780822e-01_DP + T * (  5.26219754e-02_DP &
	   + T * ( -3.12983433e-05_DP + T * (  8.84817377e-09_DP &
	   + T * ( -9.61252202e-13_DP ) ) ) ) )
      H(SA1C2H4XC) =  7.90697276e+01_DP * ( &
	   T * (  1.61326962e+01_DP + T * (  1.41452137e-02_DP &
	   + T * ( -3.39339587e-06_DP + T * (  4.10441592e-10_DP &
	   + T * ( -1.96275066e-14_DP ) ) ) ) ) +  2.08791061e+04_DP )
      CP(SA1C2H4XC) =  7.90697276e+01_DP * ( &
	    1.61326962e+01_DP + T * (  2.82904273e-02_DP &
	   + T * ( -1.01801876e-05_DP + T * (  1.64176637e-09_DP &
	   + T * ( -9.81375329e-14_DP ) ) ) ) )
      H(SA1C2H5XC) =  7.83189525e+01_DP * ( &
	   T * (  1.56901336e+01_DP + T * (  1.61831537e-02_DP &
	   + T * ( -3.89548593e-06_DP + T * (  4.72473905e-10_DP &
	   + T * ( -2.26403582e-14_DP ) ) ) ) ) -4.38669907e+03_DP )
      CP(SA1C2H5XC) =  7.83189525e+01_DP * ( &
	    1.56901336e+01_DP + T * (  3.23663075e-02_DP &
	   + T * ( -1.16864578e-05_DP + T * (  1.88989562e-09_DP &
	   + T * ( -1.13201791e-13_DP ) ) ) ) )
      H(SC8H9O2) =  6.06213544e+01_DP * ( &
	   T * (  1.74946727e+00_DP + T * (  3.66934387e-02_DP &
	   + T * ( -1.47880348e-05_DP + T * (  3.17289810e-09_DP &
	   + T * ( -2.78067684e-13_DP ) ) ) ) ) +  1.71120874e+04_DP )
      CP(SC8H9O2) =  6.06213544e+01_DP * ( &
	    1.74946727e+00_DP + T * (  7.33868774e-02_DP &
	   + T * ( -4.43641045e-05_DP + T * (  1.26915924e-08_DP &
	   + T * ( -1.39033842e-12_DP ) ) ) ) )
      H(SC8H8OOH) =  6.06213544e+01_DP * ( &
	   T * (  1.74946727e+00_DP + T * (  3.66934387e-02_DP &
	   + T * ( -1.47880348e-05_DP + T * (  3.17289810e-09_DP &
	   + T * ( -2.78067684e-13_DP ) ) ) ) ) +  1.71120874e+04_DP )
      CP(SC8H8OOH) =  6.06213544e+01_DP * ( &
	    1.74946727e+00_DP + T * (  7.33868774e-02_DP &
	   + T * ( -4.43641045e-05_DP + T * (  1.26915924e-08_DP &
	   + T * ( -1.39033842e-12_DP ) ) ) ) )
      H(SOC8H7OOH) =  5.46478336e+01_DP * ( &
	   T * (  5.90515216e+00_DP + T * (  3.35364725e-02_DP &
	   + T * ( -1.34427129e-05_DP + T * (  2.87414260e-09_DP &
	   + T * ( -2.51281962e-13_DP ) ) ) ) ) -2.17516302e+04_DP )
      CP(SOC8H7OOH) =  5.46478336e+01_DP * ( &
	    5.90515216e+00_DP + T * (  6.70729450e-02_DP &
	   + T * ( -4.03281386e-05_DP + T * (  1.14965704e-08_DP &
	   + T * ( -1.25640981e-12_DP ) ) ) ) )
      H(SA1CH3CH3) =  7.83189525e+01_DP * ( &
	   T * ( -1.95577967e+00_DP + T * (  3.54776361e-02_DP &
	   + T * ( -1.39991811e-05_DP + T * (  2.95930647e-09_DP &
	   + T * ( -2.56679432e-13_DP ) ) ) ) ) -6.91883225e+01_DP )
      CP(SA1CH3CH3) =  7.83189525e+01_DP * ( &
	   -1.95577967e+00_DP + T * (  7.09552723e-02_DP &
	   + T * ( -4.19975432e-05_DP + T * (  1.18372259e-08_DP &
	   + T * ( -1.28339716e-12_DP ) ) ) ) )
      H(SA1CH3CH2) =  7.90697276e+01_DP * ( &
	   T * (  2.36833258e+00_DP + T * (  3.01967514e-02_DP &
	   + T * ( -1.16835702e-05_DP + T * (  2.43202571e-09_DP &
	   + T * ( -2.08472266e-13_DP ) ) ) ) ) +  1.76344521e+04_DP )
      CP(SA1CH3CH2) =  7.90697276e+01_DP * ( &
	    2.36833258e+00_DP + T * (  6.03935028e-02_DP &
	   + T * ( -3.50507105e-05_DP + T * (  9.72810282e-09_DP &
	   + T * ( -1.04236133e-12_DP ) ) ) ) )
      H(SA1CH3CHO) =  6.92031229e+01_DP * ( &
	   T * (  9.22059379e-01_DP + T * (  3.25311564e-02_DP &
	   + T * ( -1.30067510e-05_DP + T * (  2.77167152e-09_DP &
	   + T * ( -2.41572622e-13_DP ) ) ) ) ) -1.13895901e+04_DP )
      CP(SA1CH3CHO) =  6.92031229e+01_DP * ( &
	    9.22059379e-01_DP + T * (  6.50623128e-02_DP &
	   + T * ( -3.90202529e-05_DP + T * (  1.10866861e-08_DP &
	   + T * ( -1.20786311e-12_DP ) ) ) ) )
      H(SA2CH3XC1) =  5.84734510e+01_DP * ( &
	   T * (  1.43553166e+00_DP + T * (  4.06226654e-02_DP &
	   + T * ( -1.61966967e-05_DP + T * (  3.44489735e-09_DP &
	   + T * ( -2.99849226e-13_DP ) ) ) ) ) +  1.04031044e+04_DP )
      CP(SA2CH3XC1) =  5.84734510e+01_DP * ( &
	    1.43553166e+00_DP + T * (  8.12453307e-02_DP &
	   + T * ( -4.85900900e-05_DP + T * (  1.37795894e-08_DP &
	   + T * ( -1.49924613e-12_DP ) ) ) ) )
      H(SA1CHOCH2) =  6.97886449e+01_DP * ( &
	   T * (  5.25955102e+00_DP + T * (  2.72439788e-02_DP &
	   + T * ( -1.06895200e-05_DP + T * (  2.24419442e-09_DP &
	   + T * ( -1.93360743e-13_DP ) ) ) ) ) +  6.61511721e+03_DP )
      CP(SA1CHOCH2) =  6.97886449e+01_DP * ( &
	    5.25955102e+00_DP + T * (  5.44879576e-02_DP &
	   + T * ( -3.20685601e-05_DP + T * (  8.97677768e-09_DP &
	   + T * ( -9.66803716e-13_DP ) ) ) ) )
      H(SA1CHOCHO) =  6.19881009e+01_DP * ( &
	   T * (  3.92023646e+00_DP + T * (  2.95004798e-02_DP &
	   + T * ( -1.19741621e-05_DP + T * (  2.57441813e-09_DP &
	   + T * ( -2.25560976e-13_DP ) ) ) ) ) -2.20802686e+04_DP )
      CP(SA1CHOCHO) =  6.19881009e+01_DP * ( &
	    3.92023646e+00_DP + T * (  5.90009596e-02_DP &
	   + T * ( -3.59224864e-05_DP + T * (  1.02976725e-08_DP &
	   + T * ( -1.12780488e-12_DP ) ) ) ) )
      H(SA2OHXC10) =  5.76727893e+01_DP * ( &
	   T * (  2.08930252e+01_DP + T * (  1.55280033e-02_DP &
	   + T * ( -3.81358540e-06_DP + T * (  4.69682165e-10_DP &
	   + T * ( -2.27647762e-14_DP ) ) ) ) ) -1.35886443e+04_DP )
      CP(SA2OHXC10) =  5.76727893e+01_DP * ( &
	    2.08930252e+01_DP + T * (  3.10560066e-02_DP &
	   + T * ( -1.14407562e-05_DP + T * (  1.87872866e-09_DP &
	   + T * ( -1.13823881e-13_DP ) ) ) ) )
      H(SA2CH2XC1) =  5.88909351e+01_DP * ( &
	   T * (  4.97463689e+00_DP + T * (  3.55734585e-02_DP &
	   + T * ( -1.40082277e-05_DP + T * (  2.95063847e-09_DP &
	   + T * ( -2.54923676e-13_DP ) ) ) ) ) +  2.96267836e+04_DP )
      CP(SA2CH2XC1) =  5.88909351e+01_DP * ( &
	    4.97463689e+00_DP + T * (  7.11469171e-02_DP &
	   + T * ( -4.20246831e-05_DP + T * (  1.18025539e-08_DP &
	   + T * ( -1.27461838e-12_DP ) ) ) ) )
      H(SA2CH2OXC) =  5.28962604e+01_DP * ( &
	   T * (  3.88859655e+00_DP + T * (  3.90441740e-02_DP &
	   + T * ( -1.56825396e-05_DP + T * (  3.34906190e-09_DP &
	   + T * ( -2.92138322e-13_DP ) ) ) ) ) +  2.16424127e+04_DP )
      CP(SA2CH2OXC) =  5.28962604e+01_DP * ( &
	    3.88859655e+00_DP + T * (  7.80883481e-02_DP &
	   + T * ( -4.70476187e-05_DP + T * (  1.33962476e-08_DP &
	   + T * ( -1.46069161e-12_DP ) ) ) ) )
      H(SA2CHOXC1) =  5.32376708e+01_DP * ( &
	   T * ( -5.55792190e-01_DP + T * (  4.36874222e-02_DP &
	   + T * ( -1.87272784e-05_DP + T * (  4.18083215e-09_DP &
	   + T * ( -3.76185702e-13_DP ) ) ) ) ) +  1.03455228e+03_DP )
      CP(SA2CHOXC1) =  5.32376708e+01_DP * ( &
	   -5.55792190e-01_DP + T * (  8.73748443e-02_DP &
	   + T * ( -5.61818352e-05_DP + T * (  1.67233286e-08_DP &
	   + T * ( -1.88092851e-12_DP ) ) ) ) )
      H(SA2OXC10H) =  5.80788790e+01_DP * ( &
	   T * (  2.10591364e+01_DP + T * (  1.41281535e-02_DP &
	   + T * ( -3.44428953e-06_DP + T * (  4.22167585e-10_DP &
	   + T * ( -2.03949534e-14_DP ) ) ) ) ) +  4.09143507e+03_DP )
      CP(SA2OXC10H) =  5.80788790e+01_DP * ( &
	    2.10591364e+01_DP + T * (  2.82563070e-02_DP &
	   + T * ( -1.03328686e-05_DP + T * (  1.68867034e-09_DP &
	   + T * ( -1.01974767e-13_DP ) ) ) ) )
      H(SOC6H4O) =  7.69191059e+01_DP * ( &
	   T * (  5.70290193e+00_DP + T * (  1.92522500e-02_DP &
	   + T * ( -7.62482393e-06_DP + T * (  1.60822969e-09_DP &
	   + T * ( -1.38851568e-13_DP ) ) ) ) ) -1.40768967e+04_DP )
      CP(SOC6H4O) =  7.69191059e+01_DP * ( &
	    5.70290193e+00_DP + T * (  3.85045001e-02_DP &
	   + T * ( -2.28744718e-05_DP + T * (  6.43291874e-09_DP &
	   + T * ( -6.94257842e-13_DP ) ) ) ) )
      ELSE
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  3.29867700e+00_DP + T * (  7.04120200e-04_DP &
	   + T * ( -1.32107400e-06_DP + T * (  1.41037875e-09_DP &
	   + T * ( -4.88970800e-13_DP ) ) ) ) ) -1.02089990e+03_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    3.29867700e+00_DP + T * (  1.40824040e-03_DP &
	   + T * ( -3.96322200e-06_DP + T * (  5.64151500e-09_DP &
	   + T * ( -2.44485400e-12_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  3.52666409e-13_DP &
	   + T * ( -6.65306547e-16_DP + T * (  5.75204080e-19_DP &
	   + T * ( -1.85546466e-22_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  7.05332819e-13_DP &
	   + T * ( -1.99591964e-15_DP + T * (  2.30081632e-18_DP &
	   + T * ( -9.27732332e-22_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.78245636e+00_DP + T * ( -1.49836708e-03_DP &
	   + T * (  3.28243400e-06_DP + T * ( -2.42032377e-09_DP &
	   + T * (  6.48745674e-13_DP ) ) ) ) ) -1.06394356e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.78245636e+00_DP + T * ( -2.99673416e-03_DP &
	   + T * (  9.84730201e-06_DP + T * ( -9.68129509e-09_DP &
	   + T * (  3.24372837e-12_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  3.16826710e+00_DP + T * ( -1.63965942e-03_DP &
	   + T * (  2.21435465e-06_DP + T * ( -1.53201656e-09_DP &
	   + T * (  4.22531942e-13_DP ) ) ) ) ) +  2.91222592e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    3.16826710e+00_DP + T * ( -3.27931884e-03_DP &
	   + T * (  6.64306396e-06_DP + T * ( -6.12806624e-09_DP &
	   + T * (  2.11265971e-12_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  4.12530561e+00_DP + T * ( -1.61272470e-03_DP &
	   + T * (  2.17588230e-06_DP + T * ( -1.44963411e-09_DP &
	   + T * (  4.12474758e-13_DP ) ) ) ) ) +  3.38153812e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    4.12530561e+00_DP + T * ( -3.22544939e-03_DP &
	   + T * (  6.52764691e-06_DP + T * ( -5.79853643e-09_DP &
	   + T * (  2.06237379e-12_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.34433112e+00_DP + T * (  3.99026037e-03_DP &
	   + T * ( -6.49271700e-06_DP + T * (  5.03930235e-09_DP &
	   + T * ( -1.47522352e-12_DP ) ) ) ) ) -9.17935173e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.34433112e+00_DP + T * (  7.98052075e-03_DP &
	   + T * ( -1.94781510e-05_DP + T * (  2.01572094e-08_DP &
	   + T * ( -7.37611761e-12_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  4.19864056e+00_DP + T * ( -1.01821705e-03_DP &
	   + T * (  2.17346737e-06_DP + T * ( -1.37199266e-09_DP &
	   + T * (  3.54395634e-13_DP ) ) ) ) ) -3.02937267e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    4.19864056e+00_DP + T * ( -2.03643410e-03_DP &
	   + T * (  6.52040211e-06_DP + T * ( -5.48797062e-09_DP &
	   + T * (  1.77197817e-12_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  2.35677352e+00_DP + T * (  4.49229839e-03_DP &
	   + T * ( -2.37452090e-06_DP + T * (  6.14797555e-10_DP &
	   + T * ( -2.87399096e-14_DP ) ) ) ) ) -4.83719697e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    2.35677352e+00_DP + T * (  8.98459677e-03_DP &
	   + T * ( -7.12356269e-06_DP + T * (  2.45919022e-09_DP &
	   + T * ( -1.43699548e-13_DP ) ) ) ) )
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
      H(SCO) =  2.96834702e+02_DP * ( &
	   T * (  3.57953347e+00_DP + T * ( -3.05176840e-04_DP &
	   + T * (  3.38938110e-07_DP + T * (  2.26751471e-10_DP &
	   + T * ( -1.80884900e-13_DP ) ) ) ) ) -1.43440860e+04_DP )
      CP(SCO) =  2.96834702e+02_DP * ( &
	    3.57953347e+00_DP + T * ( -6.10353680e-04_DP &
	   + T * (  1.01681433e-06_DP + T * (  9.07005884e-10_DP &
	   + T * ( -9.04424499e-13_DP ) ) ) ) )
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  4.22118584e+00_DP + T * ( -1.62196266e-03_DP &
	   + T * (  4.59331487e-06_DP + T * ( -3.32860233e-09_DP &
	   + T * (  8.67537730e-13_DP ) ) ) ) ) +  3.83956496e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    4.22118584e+00_DP + T * ( -3.24392532e-03_DP &
	   + T * (  1.37799446e-05_DP + T * ( -1.33144093e-08_DP &
	   + T * (  4.33768865e-12_DP ) ) ) ) )
      H(SC) =  6.92284763e+02_DP * ( &
	   T * (  2.55423955e+00_DP + T * ( -1.60768862e-04_DP &
	   + T * (  2.44597415e-07_DP + T * ( -1.83058722e-10_DP &
	   + T * (  5.33042892e-14_DP ) ) ) ) ) +  8.54438832e+04_DP )
      CP(SC) =  6.92284763e+02_DP * ( &
	    2.55423955e+00_DP + T * ( -3.21537724e-04_DP &
	   + T * (  7.33792245e-07_DP + T * ( -7.32234889e-10_DP &
	   + T * (  2.66521446e-13_DP ) ) ) ) )
      H(SCH) =  6.38680289e+02_DP * ( &
	   T * (  3.48981665e+00_DP + T * (  1.61917771e-04_DP &
	   + T * ( -5.62996883e-07_DP + T * (  7.90543317e-10_DP &
	   + T * ( -2.81218134e-13_DP ) ) ) ) ) +  7.07972934e+04_DP )
      CP(SCH) =  6.38680289e+02_DP * ( &
	    3.48981665e+00_DP + T * (  3.23835541e-04_DP &
	   + T * ( -1.68899065e-06_DP + T * (  3.16217327e-09_DP &
	   + T * ( -1.40609067e-12_DP ) ) ) ) )
      H(STXCH2) =  5.92780550e+02_DP * ( &
	   T * (  3.76267867e+00_DP + T * (  4.84436072e-04_DP &
	   + T * (  9.31632803e-07_DP + T * ( -9.62727883e-10_DP &
	   + T * (  3.37483438e-13_DP ) ) ) ) ) +  4.60040401e+04_DP )
      CP(STXCH2) =  5.92780550e+02_DP * ( &
	    3.76267867e+00_DP + T * (  9.68872143e-04_DP &
	   + T * (  2.79489841e-06_DP + T * ( -3.85091153e-09_DP &
	   + T * (  1.68741719e-12_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  3.65717970e+00_DP + T * (  1.06329895e-03_DP &
	   + T * (  1.81946277e-06_DP + T * ( -1.65452507e-09_DP &
	   + T * (  4.93141480e-13_DP ) ) ) ) ) +  1.64227160e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    3.65717970e+00_DP + T * (  2.12659790e-03_DP &
	   + T * (  5.45838830e-06_DP + T * ( -6.61810030e-09_DP &
	   + T * (  2.46570740e-12_DP ) ) ) ) )
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  4.79372315e+00_DP + T * ( -4.95416684e-03_DP &
	   + T * (  1.24406669e-05_DP + T * ( -9.48213152e-09_DP &
	   + T * (  2.63545304e-12_DP ) ) ) ) ) -1.43089567e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    4.79372315e+00_DP + T * ( -9.90833369e-03_DP &
	   + T * (  3.73220008e-05_DP + T * ( -3.79285261e-08_DP &
	   + T * (  1.31772652e-11_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  2.25172140e+00_DP + T * (  8.82751050e-03_DP &
	   + T * ( -7.90970033e-06_DP + T * (  4.31893975e-09_DP &
	   + T * ( -1.01329622e-12_DP ) ) ) ) ) +  2.00594490e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    2.25172140e+00_DP + T * (  1.76550210e-02_DP &
	   + T * ( -2.37291010e-05_DP + T * (  1.72757590e-08_DP &
	   + T * ( -5.06648110e-12_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  2.88965733e+00_DP + T * (  6.70498055e-03_DP &
	   + T * ( -9.49231670e-06_DP + T * (  7.36977613e-09_DP &
	   + T * ( -2.18663022e-12_DP ) ) ) ) ) +  6.68393932e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    2.88965733e+00_DP + T * (  1.34099611e-02_DP &
	   + T * ( -2.84769501e-05_DP + T * (  2.94791045e-08_DP &
	   + T * ( -1.09331511e-11_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  2.13583630e+00_DP + T * (  9.05943605e-03_DP &
	   + T * ( -5.79824913e-06_DP + T * (  2.33599392e-09_DP &
	   + T * ( -4.02915230e-13_DP ) ) ) ) ) -7.04291804e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    2.13583630e+00_DP + T * (  1.81188721e-02_DP &
	   + T * ( -1.73947474e-05_DP + T * (  9.34397568e-09_DP &
	   + T * ( -2.01457615e-12_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  8.08681094e-01_DP + T * (  1.16807815e-02_DP &
	   + T * ( -1.18390605e-05_DP + T * (  7.00381092e-09_DP &
	   + T * ( -1.70014595e-12_DP ) ) ) ) ) +  2.64289807e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    8.08681094e-01_DP + T * (  2.33615629e-02_DP &
	   + T * ( -3.55171815e-05_DP + T * (  2.80152437e-08_DP &
	   + T * ( -8.50072974e-12_DP ) ) ) ) )
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  4.19860411e+00_DP + T * ( -1.18330710e-03_DP &
	   + T * (  2.74432073e-06_DP + T * ( -1.67203995e-09_DP &
	   + T * (  3.88629474e-13_DP ) ) ) ) ) +  5.04968163e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    4.19860411e+00_DP + T * ( -2.36661419e-03_DP &
	   + T * (  8.23296220e-06_DP + T * ( -6.68815981e-09_DP &
	   + T * (  1.94314737e-12_DP ) ) ) ) )
      H(SAR) =  2.08129068e+02_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -7.45375000e+02_DP )
      CP(SAR) =  2.08129068e+02_DP * ( &
	    2.50000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP + T * (  0.00000000e+00_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SCH3OH) =  2.59482554e+02_DP * ( &
	   T * (  5.71539582e+00_DP + T * ( -7.61545645e-03_DP &
	   + T * (  2.17480385e-05_DP + T * ( -1.77701722e-08_DP &
	   + T * (  5.22705396e-12_DP ) ) ) ) ) -2.56427656e+04_DP )
      CP(SCH3OH) =  2.59482554e+02_DP * ( &
	    5.71539582e+00_DP + T * ( -1.52309129e-02_DP &
	   + T * (  6.52441155e-05_DP + T * ( -7.10806889e-08_DP &
	   + T * (  2.61352698e-11_DP ) ) ) ) )
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
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  5.14911468e+00_DP + T * ( -6.83110045e-03_DP &
	   + T * (  1.63817974e-05_DP + T * ( -1.21061692e-08_DP &
	   + T * (  3.33206882e-12_DP ) ) ) ) ) -1.02465983e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    5.14911468e+00_DP + T * ( -1.36622009e-02_DP &
	   + T * (  4.91453921e-05_DP + T * ( -4.84246767e-08_DP &
	   + T * (  1.66603441e-11_DP ) ) ) ) )
      H(SCH3O2) =  1.76772973e+02_DP * ( &
	   T * (  4.76597792e+00_DP + T * ( -1.75538574e-03_DP &
	   + T * (  1.51464717e-05_DP + T * ( -1.41690932e-08_DP &
	   + T * (  4.43182964e-12_DP ) ) ) ) ) -4.82401289e+02_DP )
      CP(SCH3O2) =  1.76772973e+02_DP * ( &
	    4.76597792e+00_DP + T * ( -3.51077148e-03_DP &
	   + T * (  4.54394152e-05_DP + T * ( -5.66763729e-08_DP &
	   + T * (  2.21591482e-11_DP ) ) ) ) )
      H(SC2H3) =  3.07437509e+02_DP * ( &
	   T * (  3.21246645e+00_DP + T * (  7.57395810e-04_DP &
	   + T * (  8.64031373e-06_DP + T * ( -8.94144617e-09_DP &
	   + T * (  2.94301746e-12_DP ) ) ) ) ) +  3.48598468e+04_DP )
      CP(SC2H3) =  3.07437509e+02_DP * ( &
	    3.21246645e+00_DP + T * (  1.51479162e-03_DP &
	   + T * (  2.59209412e-05_DP + T * ( -3.57657847e-08_DP &
	   + T * (  1.47150873e-11_DP ) ) ) ) )
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
      H(SHCCOH) =  1.97790941e+02_DP * ( &
	   T * (  1.24237330e+00_DP + T * (  1.55361005e-02_DP &
	   + T * ( -1.69556213e-05_DP + T * (  1.07842828e-08_DP &
	   + T * ( -2.80291880e-12_DP ) ) ) ) ) +  8.03161430e+03_DP )
      CP(SHCCOH) =  1.97790941e+02_DP * ( &
	    1.24237330e+00_DP + T * (  3.10722010e-02_DP &
	   + T * ( -5.08668640e-05_DP + T * (  4.31371310e-08_DP &
	   + T * ( -1.40145940e-11_DP ) ) ) ) )
      H(SCH2CHO) =  1.93159093e+02_DP * ( &
	   T * (  1.09685733e+00_DP + T * (  1.10114398e-02_DP &
	   + T * ( -4.81944813e-06_DP + T * (  7.51948945e-10_DP &
	   + T * (  1.21798575e-13_DP ) ) ) ) ) +  1.06943322e+03_DP )
      CP(SCH2CHO) =  1.93159093e+02_DP * ( &
	    1.09685733e+00_DP + T * (  2.20228796e-02_DP &
	   + T * ( -1.44583444e-05_DP + T * (  3.00779578e-09_DP &
	   + T * (  6.08992877e-13_DP ) ) ) ) )
      H(SCH3CHO) =  1.88739217e+02_DP * ( &
	   T * (  1.40653856e+00_DP + T * (  1.08492219e-02_DP &
	   + T * ( -4.91910883e-06_DP + T * (  1.82608869e-09_DP &
	   + T * ( -4.18238934e-13_DP ) ) ) ) ) -2.17973223e+04_DP )
      CP(SCH3CHO) =  1.88739217e+02_DP * ( &
	    1.40653856e+00_DP + T * (  2.16984438e-02_DP &
	   + T * ( -1.47573265e-05_DP + T * (  7.30435478e-09_DP &
	   + T * ( -2.09119467e-12_DP ) ) ) ) )
      H(SH2C2) =  3.19340144e+02_DP * ( &
	   T * (  3.28154830e+00_DP + T * (  3.48823955e-03_DP &
	   + T * ( -7.95174800e-07_DP + T * ( -3.02610800e-10_DP &
	   + T * (  1.96379090e-13_DP ) ) ) ) ) +  4.86217940e+04_DP )
      CP(SH2C2) =  3.19340144e+02_DP * ( &
	    3.28154830e+00_DP + T * (  6.97647910e-03_DP &
	   + T * ( -2.38552440e-06_DP + T * ( -1.21044320e-09_DP &
	   + T * (  9.81895450e-13_DP ) ) ) ) )
      H(SC2H5O) =  1.84517088e+02_DP * ( &
	   T * (  4.94420708e-01_DP + T * (  1.35887217e-02_DP &
	   + T * ( -5.53030033e-06_DP + T * (  1.28801050e-09_DP &
	   + T * ( -1.29699383e-13_DP ) ) ) ) ) -3.35252925e+03_DP )
      CP(SC2H5O) =  1.84517088e+02_DP * ( &
	    4.94420708e-01_DP + T * (  2.71774434e-02_DP &
	   + T * ( -1.65909010e-05_DP + T * (  5.15204200e-09_DP &
	   + T * ( -6.48496915e-13_DP ) ) ) ) )
      H(SNXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  1.04754730e+00_DP + T * (  1.30038970e-02_DP &
	   + T * (  7.85408400e-07_DP + T * ( -4.89807925e-09_DP &
	   + T * (  1.87360232e-12_DP ) ) ) ) ) +  1.06326370e+04_DP )
      CP(SNXC3H7) =  1.92970803e+02_DP * ( &
	    1.04754730e+00_DP + T * (  2.60077940e-02_DP &
	   + T * (  2.35622520e-06_DP + T * ( -1.95923170e-08_DP &
	   + T * (  9.36801160e-12_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  4.29142492e+00_DP + T * ( -2.75077135e-03_DP &
	   + T * (  1.99812763e-05_DP + T * ( -1.77116571e-08_DP &
	   + T * (  5.37371542e-12_DP ) ) ) ) ) -1.15222055e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    4.29142492e+00_DP + T * ( -5.50154270e-03_DP &
	   + T * (  5.99438288e-05_DP + T * ( -7.08466285e-08_DP &
	   + T * (  2.68685771e-11_DP ) ) ) ) )
      H(SC3H8) =  1.88559441e+02_DP * ( &
	   T * (  9.33553810e-01_DP + T * (  1.32122895e-02_DP &
	   + T * (  2.03532423e-06_DP + T * ( -5.49437475e-09_DP &
	   + T * (  1.90298506e-12_DP ) ) ) ) ) -1.39585200e+04_DP )
      CP(SC3H8) =  1.88559441e+02_DP * ( &
	    9.33553810e-01_DP + T * (  2.64245790e-02_DP &
	   + T * (  6.10597270e-06_DP + T * ( -2.19774990e-08_DP &
	   + T * (  9.51492530e-12_DP ) ) ) ) )
      H(SC3H6) =  1.97593517e+02_DP * ( &
	   T * ( -2.29261670e-03_DP + T * (  1.55130533e-02_DP &
	   + T * ( -5.57171827e-06_DP + T * (  4.73985425e-10_DP &
	   + T * (  2.49915830e-13_DP ) ) ) ) ) +  1.13437406e+03_DP )
      CP(SC3H6) =  1.97593517e+02_DP * ( &
	   -2.29261670e-03_DP + T * (  3.10261065e-02_DP &
	   + T * ( -1.67151548e-05_DP + T * (  1.89594170e-09_DP &
	   + T * (  1.24957915e-12_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  1.40299238e+00_DP + T * (  1.50886664e-02_DP &
	   + T * ( -1.32816458e-05_DP + T * (  7.33836572e-09_DP &
	   + T * ( -1.74110916e-12_DP ) ) ) ) ) +  3.93108220e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    1.40299238e+00_DP + T * (  3.01773327e-02_DP &
	   + T * ( -3.98449373e-05_DP + T * (  2.93534629e-08_DP &
	   + T * ( -8.70554579e-12_DP ) ) ) ) )
      H(SPXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  1.46175323e+00_DP + T * (  1.23013301e-02_DP &
	   + T * ( -6.34064650e-06_DP + T * (  2.15090855e-09_DP &
	   + T * ( -3.33458480e-13_DP ) ) ) ) ) +  2.09209793e+04_DP )
      CP(SPXC3H4) =  2.07536818e+02_DP * ( &
	    1.46175323e+00_DP + T * (  2.46026602e-02_DP &
	   + T * ( -1.90219395e-05_DP + T * (  8.60363422e-09_DP &
	   + T * ( -1.66729240e-12_DP ) ) ) ) )
      H(SAXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  3.68928265e-01_DP + T * (  1.44675699e-02_DP &
	   + T * ( -8.14621360e-06_DP + T * (  2.81367915e-09_DP &
	   + T * ( -4.06080524e-13_DP ) ) ) ) ) +  2.17585256e+04_DP )
      CP(SAXC3H4) =  2.07536818e+02_DP * ( &
	    3.68928265e-01_DP + T * (  2.89351397e-02_DP &
	   + T * ( -2.44386408e-05_DP + T * (  1.12547166e-08_DP &
	   + T * ( -2.03040262e-12_DP ) ) ) ) )
      H(SSXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  3.13106581e-01_DP + T * (  1.59384831e-02_DP &
	   + T * ( -8.44733377e-06_DP + T * (  2.57497682e-09_DP &
	   + T * ( -2.70603708e-13_DP ) ) ) ) ) +  3.13767683e+04_DP )
      CP(SSXC3H5) =  2.02443146e+02_DP * ( &
	    3.13106581e-01_DP + T * (  3.18769663e-02_DP &
	   + T * ( -2.53420013e-05_DP + T * (  1.02999073e-08_DP &
	   + T * ( -1.35301854e-12_DP ) ) ) ) )
      H(SNXC4H3) =  1.62821949e+02_DP * ( &
	   T * ( -3.55175031e-02_DP + T * (  2.15254251e-02_DP &
	   + T * ( -1.91909716e-05_DP + T * (  1.03970785e-08_DP &
	   + T * ( -2.41501714e-12_DP ) ) ) ) ) +  6.43506593e+04_DP )
      CP(SNXC4H3) =  1.62821949e+02_DP * ( &
	   -3.55175031e-02_DP + T * (  4.30508503e-02_DP &
	   + T * ( -5.75729147e-05_DP + T * (  4.15883142e-08_DP &
	   + T * ( -1.20750857e-11_DP ) ) ) ) )
      H(SC2H3CHO) =  1.48306161e+02_DP * ( &
	   T * (  2.92355162e-01_DP + T * (  1.77160709e-02_DP &
	   + T * ( -9.83121080e-06_DP + T * (  3.20250310e-09_DP &
	   + T * ( -4.52288216e-13_DP ) ) ) ) ) -1.16521584e+04_DP )
      CP(SC2H3CHO) =  1.48306161e+02_DP * ( &
	    2.92355162e-01_DP + T * (  3.54321417e-02_DP &
	   + T * ( -2.94936324e-05_DP + T * (  1.28100124e-08_DP &
	   + T * ( -2.26144108e-12_DP ) ) ) ) )
      H(SAXC3H5) =  2.02443146e+02_DP * ( &
	   T * ( -1.03516444e+00_DP + T * (  1.87521683e-02_DP &
	   + T * ( -1.08793747e-05_DP + T * (  3.69156533e-09_DP &
	   + T * ( -4.87482308e-13_DP ) ) ) ) ) +  1.88792254e+04_DP )
      CP(SAXC3H5) =  2.02443146e+02_DP * ( &
	   -1.03516444e+00_DP + T * (  3.75043366e-02_DP &
	   + T * ( -3.26381242e-05_DP + T * (  1.47662613e-08_DP &
	   + T * ( -2.43741154e-12_DP ) ) ) ) )
      H(SC2O) =  2.07754623e+02_DP * ( &
	   T * (  2.86278214e+00_DP + T * (  5.98506020e-03_DP &
	   + T * ( -6.02837407e-06_DP + T * (  3.81944325e-09_DP &
	   + T * ( -1.04012633e-12_DP ) ) ) ) ) +  3.37501779e+04_DP )
      CP(SC2O) =  2.07754623e+02_DP * ( &
	    2.86278214e+00_DP + T * (  1.19701204e-02_DP &
	   + T * ( -1.80851222e-05_DP + T * (  1.52777730e-08_DP &
	   + T * ( -5.20063163e-12_DP ) ) ) ) )
      H(SC4H4) =  1.59670072e+02_DP * ( &
	   T * ( -2.31343354e-01_DP + T * (  2.05907248e-02_DP &
	   + T * ( -1.49208019e-05_DP + T * (  6.88585392e-09_DP &
	   + T * ( -1.41275363e-12_DP ) ) ) ) ) +  3.40632704e+04_DP )
      CP(SC4H4) =  1.59670072e+02_DP * ( &
	   -2.31343354e-01_DP + T * (  4.11814497e-02_DP &
	   + T * ( -4.47624056e-05_DP + T * (  2.75434157e-08_DP &
	   + T * ( -7.06376813e-12_DP ) ) ) ) )
      H(SC3H2) =  2.18533880e+02_DP * ( &
	   T * (  4.52861333e+00_DP + T * (  8.87827505e-03_DP &
	   + T * ( -8.49609820e-06_DP + T * (  5.04186573e-09_DP &
	   + T * ( -1.25708941e-12_DP ) ) ) ) ) +  6.35410087e+04_DP )
      CP(SC3H2) =  2.18533880e+02_DP * ( &
	    4.52861333e+00_DP + T * (  1.77565501e-02_DP &
	   + T * ( -2.54882946e-05_DP + T * (  2.01674629e-08_DP &
	   + T * ( -6.28544707e-12_DP ) ) ) ) )
      H(SC3H2O) =  1.53838212e+02_DP * ( &
	   T * (  1.89401982e+00_DP + T * (  1.33150743e-02_DP &
	   + T * ( -9.90617387e-06_DP + T * (  4.85725965e-09_DP &
	   + T * ( -1.08680553e-12_DP ) ) ) ) ) +  1.37271761e+04_DP )
      CP(SC3H2O) =  1.53838212e+02_DP * ( &
	    1.89401982e+00_DP + T * (  2.66301486e-02_DP &
	   + T * ( -2.97185216e-05_DP + T * (  1.94290386e-08_DP &
	   + T * ( -5.43402767e-12_DP ) ) ) ) )
      H(SC4H2) =  1.66100767e+02_DP * ( &
	   T * (  1.73325212e-01_DP + T * (  2.26974515e-02_DP &
	   + T * ( -2.43374610e-05_DP + T * (  1.48812934e-08_DP &
	   + T * ( -3.74969432e-12_DP ) ) ) ) ) +  5.42239385e+04_DP )
      CP(SC4H2) =  1.66100767e+02_DP * ( &
	    1.73325212e-01_DP + T * (  4.53949030e-02_DP &
	   + T * ( -7.30123830e-05_DP + T * (  5.95251736e-08_DP &
	   + T * ( -1.87484716e-11_DP ) ) ) ) )
      H(SIXC4H3) =  1.62821949e+02_DP * ( &
	   T * (  3.02566263e+00_DP + T * (  1.52346812e-02_DP &
	   + T * ( -1.22781728e-05_DP + T * (  6.50088380e-09_DP &
	   + T * ( -1.52430870e-12_DP ) ) ) ) ) +  5.80551505e+04_DP )
      CP(SIXC4H3) =  1.62821949e+02_DP * ( &
	    3.02566263e+00_DP + T * (  3.04693624e-02_DP &
	   + T * ( -3.68345185e-05_DP + T * (  2.60035352e-08_DP &
	   + T * ( -7.62154351e-12_DP ) ) ) ) )
      H(STXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  8.80980628e-01_DP + T * (  1.48180962e-02_DP &
	   + T * ( -8.42418673e-06_DP + T * (  3.59129540e-09_DP &
	   + T * ( -7.79133242e-13_DP ) ) ) ) ) +  2.92321259e+04_DP )
      CP(STXC3H5) =  2.02443146e+02_DP * ( &
	    8.80980628e-01_DP + T * (  2.96361924e-02_DP &
	   + T * ( -2.52725602e-05_DP + T * (  1.43651816e-08_DP &
	   + T * ( -3.89566621e-12_DP ) ) ) ) )
      H(SC3H5O) =  1.45686701e+02_DP * ( &
	   T * (  1.19822582e+00_DP + T * (  1.52789918e-02_DP &
	   + T * ( -6.02100920e-06_DP + T * (  1.21537508e-09_DP &
	   + T * ( -8.39709124e-14_DP ) ) ) ) ) +  9.58217784e+03_DP )
      CP(SC3H5O) =  1.45686701e+02_DP * ( &
	    1.19822582e+00_DP + T * (  3.05579837e-02_DP &
	   + T * ( -1.80630276e-05_DP + T * (  4.86150033e-09_DP &
	   + T * ( -4.19854562e-13_DP ) ) ) ) )
      H(SC4H) =  1.69514353e+02_DP * ( &
	   T * (  3.23559253e+00_DP + T * (  1.13545766e-02_DP &
	   + T * ( -1.06147764e-05_DP + T * (  6.12162010e-09_DP &
	   + T * ( -1.51597288e-12_DP ) ) ) ) ) +  9.39080960e+04_DP )
      CP(SC4H) =  1.69514353e+02_DP * ( &
	    3.23559253e+00_DP + T * (  2.27091533e-02_DP &
	   + T * ( -3.18443291e-05_DP + T * (  2.44864804e-08_DP &
	   + T * ( -7.57986440e-12_DP ) ) ) ) )
      H(SC8H2) =  8.47571766e+01_DP * ( &
	   T * ( -3.26701608e-01_DP + T * (  4.71664338e-02_DP &
	   + T * ( -5.76254613e-05_DP + T * (  3.92041345e-08_DP &
	   + T * ( -1.08097685e-11_DP ) ) ) ) ) +  1.05392079e+05_DP )
      CP(SC8H2) =  8.47571766e+01_DP * ( &
	   -3.26701608e-01_DP + T * (  9.43328676e-02_DP &
	   + T * ( -1.72876384e-04_DP + T * (  1.56816538e-07_DP &
	   + T * ( -5.40488426e-11_DP ) ) ) ) )
      H(SC6H2) =  1.12240672e+02_DP * ( &
	   T * ( -5.41092160e-01_DP + T * (  3.72663140e-02_DP &
	   + T * ( -4.52608400e-05_DP + T * (  3.05665750e-08_DP &
	   + T * ( -8.36504140e-12_DP ) ) ) ) ) +  8.21151320e+04_DP )
      CP(SC6H2) =  1.12240672e+02_DP * ( &
	   -5.41092160e-01_DP + T * (  7.45326280e-02_DP &
	   + T * ( -1.35782520e-04_DP + T * (  1.22266300e-07_DP &
	   + T * ( -4.18252070e-11_DP ) ) ) ) )
      H(SC4H6) =  1.53718755e+02_DP * ( &
	   T * (  4.01336263e+00_DP + T * (  2.22313425e-03_DP &
	   + T * (  2.60227673e-05_DP + T * ( -2.79185323e-08_DP &
	   + T * (  9.21507692e-12_DP ) ) ) ) ) +  1.14807231e+04_DP )
      CP(SC4H6) =  1.53718755e+02_DP * ( &
	    4.01336263e+00_DP + T * (  4.44626850e-03_DP &
	   + T * (  7.80683019e-05_DP + T * ( -1.11674129e-07_DP &
	   + T * (  4.60753846e-11_DP ) ) ) ) )
      H(SNXC4H5) =  1.56637905e+02_DP * ( &
	   T * ( -1.16849950e+00_DP + T * (  2.39503037e-02_DP &
	   + T * ( -1.70792334e-05_DP + T * (  7.65610660e-09_DP &
	   + T * ( -1.51981393e-12_DP ) ) ) ) ) +  4.22787216e+04_DP )
      CP(SNXC4H5) =  1.56637905e+02_DP * ( &
	   -1.16849950e+00_DP + T * (  4.79006074e-02_DP &
	   + T * ( -5.12377002e-05_DP + T * (  3.06244264e-08_DP &
	   + T * ( -7.59906965e-12_DP ) ) ) ) )
      H(SIXC4H5) =  1.56637905e+02_DP * ( &
	   T * ( -3.31905498e-01_DP + T * (  2.20081938e-02_DP &
	   + T * ( -1.42563415e-05_DP + T * (  5.78210790e-09_DP &
	   + T * ( -1.03434304e-12_DP ) ) ) ) ) +  3.67510686e+04_DP )
      CP(SIXC4H5) =  1.56637905e+02_DP * ( &
	   -3.31905498e-01_DP + T * (  4.40163876e-02_DP &
	   + T * ( -4.27690246e-05_DP + T * (  2.31284316e-08_DP &
	   + T * ( -5.17171519e-12_DP ) ) ) ) )
      H(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   T * ( -5.51558393e+00_DP + T * (  3.22726613e-02_DP &
	   + T * ( -1.47134309e-05_DP + T * (  1.86928040e-09_DP &
	   + T * (  6.20564508e-13_DP ) ) ) ) ) +  9.11031457e+03_DP )
      CP(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   -5.51558393e+00_DP + T * (  6.45453225e-02_DP &
	   + T * ( -4.41402928e-05_DP + T * (  7.47712161e-09_DP &
	   + T * (  3.10282254e-12_DP ) ) ) ) )
      H(SNXC7H16) =  8.29791014e+01_DP * ( &
	   T * ( -1.26836187e+00_DP + T * (  4.27177910e-02_DP &
	   + T * ( -1.75115595e-05_DP + T * (  4.07364302e-09_DP &
	   + T * ( -4.04789850e-13_DP ) ) ) ) ) -2.56586565e+04_DP )
      CP(SNXC7H16) =  8.29791014e+01_DP * ( &
	   -1.26836187e+00_DP + T * (  8.54355820e-02_DP &
	   + T * ( -5.25346786e-05_DP + T * (  1.62945721e-08_DP &
	   + T * ( -2.02394925e-12_DP ) ) ) ) )
      H(SC5H11) =  1.16876212e+02_DP * ( &
	   T * ( -9.05255912e-01_DP + T * (  3.05316426e-02_DP &
	   + T * ( -1.36497275e-05_DP + T * (  3.65233675e-09_DP &
	   + T * ( -4.37719230e-13_DP ) ) ) ) ) +  4.83995303e+03_DP )
      CP(SC5H11) =  1.16876212e+02_DP * ( &
	   -9.05255912e-01_DP + T * (  6.10632852e-02_DP &
	   + T * ( -4.09491825e-05_DP + T * (  1.46093470e-08_DP &
	   + T * ( -2.18859615e-12_DP ) ) ) ) )
      H(SPXC4H9) =  1.45579563e+02_DP * ( &
	   T * ( -4.37779725e-01_DP + T * (  2.39486182e-02_DP &
	   + T * ( -1.04674386e-05_DP + T * (  2.74466180e-09_DP &
	   + T * ( -3.24021328e-13_DP ) ) ) ) ) +  7.68945248e+03_DP )
      CP(SPXC4H9) =  1.45579563e+02_DP * ( &
	   -4.37779725e-01_DP + T * (  4.78972364e-02_DP &
	   + T * ( -3.14023159e-05_DP + T * (  1.09786472e-08_DP &
	   + T * ( -1.62010664e-12_DP ) ) ) ) )
      H(SC7H15) =  8.38223611e+01_DP * ( &
	   T * ( -3.79155767e-02_DP + T * (  3.78363285e-02_DP &
	   + T * ( -1.35824545e-05_DP + T * (  2.33169736e-09_DP &
	   + T * ( -9.84721490e-14_DP ) ) ) ) ) -2.35605303e+03_DP )
      CP(SC7H15) =  8.38223611e+01_DP * ( &
	   -3.79155767e-02_DP + T * (  7.56726570e-02_DP &
	   + T * ( -4.07473634e-05_DP + T * (  9.32678943e-09_DP &
	   + T * ( -4.92360745e-13_DP ) ) ) ) )
      H(SPXC4H8) =  1.48195138e+02_DP * ( &
	   T * ( -8.31372089e-01_DP + T * (  2.26290489e-02_DP &
	   + T * ( -9.78861863e-06_DP + T * (  2.50551090e-09_DP &
	   + T * ( -2.86383360e-13_DP ) ) ) ) ) -1.57875035e+03_DP )
      CP(SPXC4H8) =  1.48195138e+02_DP * ( &
	   -8.31372089e-01_DP + T * (  4.52580978e-02_DP &
	   + T * ( -2.93658559e-05_DP + T * (  1.00220436e-08_DP &
	   + T * ( -1.43191680e-12_DP ) ) ) ) )
      H(SC5H10) =  1.18556110e+02_DP * ( &
	   T * ( -1.06223481e+00_DP + T * (  2.87109147e-02_DP &
	   + T * ( -1.24828963e-05_DP + T * (  3.18412472e-09_DP &
	   + T * ( -3.59219578e-13_DP ) ) ) ) ) -4.46546666e+03_DP )
      CP(SC5H10) =  1.18556110e+02_DP * ( &
	   -1.06223481e+00_DP + T * (  5.74218294e-02_DP &
	   + T * ( -3.74486890e-05_DP + T * (  1.27364989e-08_DP &
	   + T * ( -1.79609789e-12_DP ) ) ) ) )
      H(SC7H14) =  8.46829358e+01_DP * ( &
	   T * ( -2.03026994e+00_DP + T * (  4.13162189e-02_DP &
	   + T * ( -1.81838157e-05_DP + T * (  4.69264555e-09_DP &
	   + T * ( -5.35142440e-13_DP ) ) ) ) ) -1.15141029e+04_DP )
      CP(SC7H14) =  8.46829358e+01_DP * ( &
	   -2.03026994e+00_DP + T * (  8.26324377e-02_DP &
	   + T * ( -5.45514471e-05_DP + T * (  1.87705822e-08_DP &
	   + T * ( -2.67571220e-12_DP ) ) ) ) )
      H(SC7H15O) =  7.21793558e+01_DP * ( &
	   T * ( -4.59189934e-01_DP + T * (  4.37232324e-02_DP &
	   + T * ( -1.89671712e-05_DP + T * (  4.80489770e-09_DP &
	   + T * ( -5.37506930e-13_DP ) ) ) ) ) -1.78233113e+04_DP )
      CP(SC7H15O) =  7.21793558e+01_DP * ( &
	   -4.59189934e-01_DP + T * (  8.74464647e-02_DP &
	   + T * ( -5.69015135e-05_DP + T * (  1.92195908e-08_DP &
	   + T * ( -2.68753465e-12_DP ) ) ) ) )
      H(SC3H7CHO) =  1.15310385e+02_DP * ( &
	   T * (  1.87415959e+00_DP + T * (  2.09620158e-02_DP &
	   + T * ( -7.83829263e-06_DP + T * (  1.56728418e-09_DP &
	   + T * ( -1.21888782e-13_DP ) ) ) ) ) -2.71032194e+04_DP )
      CP(SC3H7CHO) =  1.15310385e+02_DP * ( &
	    1.87415959e+00_DP + T * (  4.19240315e-02_DP &
	   + T * ( -2.35148779e-05_DP + T * (  6.26913673e-09_DP &
	   + T * ( -6.09443908e-13_DP ) ) ) ) )
      H(SC4H7) =  1.50906418e+02_DP * ( &
	   T * (  5.07355313e+00_DP + T * (  2.63809664e-03_DP &
	   + T * (  2.07813774e-05_DP + T * ( -2.13550865e-08_DP &
	   + T * (  6.91780062e-12_DP ) ) ) ) ) +  2.24615054e+04_DP )
      CP(SC4H7) =  1.50906418e+02_DP * ( &
	    5.07355313e+00_DP + T * (  5.27619329e-03_DP &
	   + T * (  6.23441322e-05_DP + T * ( -8.54203458e-08_DP &
	   + T * (  3.45890031e-11_DP ) ) ) ) )
      H(SC7H13) =  8.55613642e+01_DP * ( &
	   T * ( -2.01707658e+00_DP + T * (  4.04457779e-02_DP &
	   + T * ( -1.81127968e-05_DP + T * (  4.70165270e-09_DP &
	   + T * ( -5.32118506e-13_DP ) ) ) ) ) +  6.81102828e+03_DP )
      CP(SC7H13) =  8.55613642e+01_DP * ( &
	   -2.01707658e+00_DP + T * (  8.08915559e-02_DP &
	   + T * ( -5.43383904e-05_DP + T * (  1.88066108e-08_DP &
	   + T * ( -2.66059253e-12_DP ) ) ) ) )
      H(SC5H9) =  1.20285003e+02_DP * ( &
	   T * ( -1.38013950e+00_DP + T * (  2.78804243e-02_DP &
	   + T * ( -1.23381309e-05_DP + T * (  3.17209752e-09_DP &
	   + T * ( -3.57077670e-13_DP ) ) ) ) ) +  1.25589824e+04_DP )
      CP(SC5H9) =  1.20285003e+02_DP * ( &
	   -1.38013950e+00_DP + T * (  5.57608487e-02_DP &
	   + T * ( -3.70143928e-05_DP + T * (  1.26883901e-08_DP &
	   + T * ( -1.78538835e-12_DP ) ) ) ) )
      H(SC4H7O) =  1.16945257e+02_DP * ( &
	   T * ( -1.60619192e+00_DP + T * (  2.79281341e-02_DP &
	   + T * ( -1.45198589e-05_DP + T * (  4.26473197e-09_DP &
	   + T * ( -5.31270360e-13_DP ) ) ) ) ) +  4.85090326e+03_DP )
      CP(SC4H7O) =  1.16945257e+02_DP * ( &
	   -1.60619192e+00_DP + T * (  5.58562682e-02_DP &
	   + T * ( -4.35595767e-05_DP + T * (  1.70589279e-08_DP &
	   + T * ( -2.65635180e-12_DP ) ) ) ) )
      H(SNXC3H7O) =  1.40715906e+02_DP * ( &
	   T * (  2.89706514e-01_DP + T * (  1.96537780e-02_DP &
	   + T * ( -8.26897850e-06_DP + T * (  2.01770893e-09_DP &
	   + T * ( -2.15996580e-13_DP ) ) ) ) ) -6.24474269e+03_DP )
      CP(SNXC3H7O) =  1.40715906e+02_DP * ( &
	    2.89706514e-01_DP + T * (  3.93075560e-02_DP &
	   + T * ( -2.48069355e-05_DP + T * (  8.07083573e-09_DP &
	   + T * ( -1.07998290e-12_DP ) ) ) ) )
      H(SIXC8H18) =  7.27897815e+01_DP * ( &
	   T * ( -4.20868893e+00_DP + T * (  5.57202905e-02_DP &
	   + T * ( -2.63782194e-05_DP + T * (  7.31015605e-09_DP &
	   + T * ( -8.87486382e-13_DP ) ) ) ) ) -2.99446875e+04_DP )
      CP(SIXC8H18) =  7.27897815e+01_DP * ( &
	   -4.20868893e+00_DP + T * (  1.11440581e-01_DP &
	   + T * ( -7.91346582e-05_DP + T * (  2.92406242e-08_DP &
	   + T * ( -4.43743191e-12_DP ) ) ) ) )
      H(SYXC7H15) =  8.38223611e+01_DP * ( &
	   T * (  1.30897106e+00_DP + T * (  3.48068221e-02_DP &
	   + T * ( -1.10383352e-05_DP + T * (  1.45722064e-09_DP &
	   + T * (  7.08854628e-15_DP ) ) ) ) ) -5.78512513e+03_DP )
      CP(SYXC7H15) =  8.38223611e+01_DP * ( &
	    1.30897106e+00_DP + T * (  6.96136442e-02_DP &
	   + T * ( -3.31150057e-05_DP + T * (  5.82888256e-09_DP &
	   + T * (  3.54427314e-14_DP ) ) ) ) )
      H(SIXC4H8) =  1.48195138e+02_DP * ( &
	   T * (  9.38433173e-01_DP + T * (  1.95273643e-02_DP &
	   + T * ( -7.21457160e-06_DP + T * (  1.46816769e-09_DP &
	   + T * ( -1.22887096e-13_DP ) ) ) ) ) -3.74817891e+03_DP )
      CP(SIXC4H8) =  1.48195138e+02_DP * ( &
	    9.38433173e-01_DP + T * (  3.90547287e-02_DP &
	   + T * ( -2.16437148e-05_DP + T * (  5.87267077e-09_DP &
	   + T * ( -6.14435479e-13_DP ) ) ) ) )
      H(SIXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  1.71330000e+00_DP + T * (  1.27130800e-02_DP &
	   + T * (  5.26936000e-07_DP + T * ( -4.55321500e-09_DP &
	   + T * (  1.76554200e-12_DP ) ) ) ) ) +  7.53580900e+03_DP )
      CP(SIXC3H7) =  1.92970803e+02_DP * ( &
	    1.71330000e+00_DP + T * (  2.54261600e-02_DP &
	   + T * (  1.58080800e-06_DP + T * ( -1.82128600e-08_DP &
	   + T * (  8.82771000e-12_DP ) ) ) ) )
      H(STXC4H9) =  1.45579563e+02_DP * ( &
	   T * ( -2.73729203e+00_DP + T * (  2.27695172e-02_DP &
	   + T * ( -7.54636463e-06_DP + T * (  1.14237763e-09_DP &
	   + T * ( -3.10643590e-14_DP ) ) ) ) ) +  4.87138887e+03_DP )
      CP(STXC4H9) =  1.45579563e+02_DP * ( &
	   -2.73729203e+00_DP + T * (  4.55390345e-02_DP &
	   + T * ( -2.26390939e-05_DP + T * (  4.56951052e-09_DP &
	   + T * ( -1.55321795e-13_DP ) ) ) ) )
      H(SCXC8H17) =  7.34378533e+01_DP * ( &
	   T * ( -9.73159697e-02_DP + T * (  4.46326862e-02_DP &
	   + T * ( -1.70957938e-05_DP + T * (  3.44101320e-09_DP &
	   + T * ( -2.55576792e-13_DP ) ) ) ) ) -8.81147302e+03_DP )
      CP(SCXC8H17) =  7.34378533e+01_DP * ( &
	   -9.73159697e-02_DP + T * (  8.92653724e-02_DP &
	   + T * ( -5.12873814e-05_DP + T * (  1.37640528e-08_DP &
	   + T * ( -1.27788396e-12_DP ) ) ) ) )
      H(SYXC7H14) =  8.46829358e+01_DP * ( &
	   T * ( -8.42232649e-01_DP + T * (  3.94899149e-02_DP &
	   + T * ( -1.68191492e-05_DP + T * (  4.22335120e-09_DP &
	   + T * ( -4.74403974e-13_DP ) ) ) ) ) -1.45971538e+04_DP )
      CP(SYXC7H14) =  8.46829358e+01_DP * ( &
	   -8.42232649e-01_DP + T * (  7.89798297e-02_DP &
	   + T * ( -5.04574475e-05_DP + T * (  1.68934048e-08_DP &
	   + T * ( -2.37201987e-12_DP ) ) ) ) )
      H(SDXC8H17O) =  6.43445084e+01_DP * ( &
	   T * ( -3.80491312e+00_DP + T * (  5.75923045e-02_DP &
	   + T * ( -2.85871395e-05_DP + T * (  8.29694817e-09_DP &
	   + T * ( -1.04895079e-12_DP ) ) ) ) ) -2.19418053e+04_DP )
      CP(SDXC8H17O) =  6.43445084e+01_DP * ( &
	   -3.80491312e+00_DP + T * (  1.15184609e-01_DP &
	   + T * ( -8.57614185e-05_DP + T * (  3.31877927e-08_DP &
	   + T * ( -5.24475393e-12_DP ) ) ) ) )
      H(SCH3COCH3) =  1.43158167e+02_DP * ( &
	   T * (  5.55638920e+00_DP + T * ( -1.41931773e-03_DP &
	   + T * (  2.35240984e-05_DP + T * ( -2.19532746e-08_DP &
	   + T * (  6.80581902e-12_DP ) ) ) ) ) -2.78325393e+04_DP )
      CP(SCH3COCH3) =  1.43158167e+02_DP * ( &
	    5.55638920e+00_DP + T * ( -2.83863547e-03_DP &
	   + T * (  7.05722951e-05_DP + T * ( -8.78130984e-08_DP &
	   + T * (  3.40290951e-11_DP ) ) ) ) )
      H(SIXC4H7) =  1.50906418e+02_DP * ( &
	   T * ( -7.20881697e-04_DP + T * (  2.18247865e-02_DP &
	   + T * ( -1.05461959e-05_DP + T * (  3.09962457e-09_DP &
	   + T * ( -4.08756720e-13_DP ) ) ) ) ) +  1.45717785e+04_DP )
      CP(SIXC4H7) =  1.50906418e+02_DP * ( &
	   -7.20881697e-04_DP + T * (  4.36495730e-02_DP &
	   + T * ( -3.16385877e-05_DP + T * (  1.23984983e-08_DP &
	   + T * ( -2.04378360e-12_DP ) ) ) ) )
      H(SXXC7H13) =  8.55613642e+01_DP * ( &
	   T * ( -3.06783292e-01_DP + T * (  3.60345573e-02_DP &
	   + T * ( -1.40077972e-05_DP + T * (  3.00125220e-09_DP &
	   + T * ( -2.65210428e-13_DP ) ) ) ) ) +  1.23740449e+03_DP )
      CP(SXXC7H13) =  8.55613642e+01_DP * ( &
	   -3.06783292e-01_DP + T * (  7.20691145e-02_DP &
	   + T * ( -4.20233916e-05_DP + T * (  1.20050088e-08_DP &
	   + T * ( -1.32605214e-12_DP ) ) ) ) )
      H(SIXC3H5CH) =  1.18627154e+02_DP * ( &
	   T * (  6.27183793e-01_DP + T * (  2.33390127e-02_DP &
	   + T * ( -1.24810210e-05_DP + T * (  3.95826355e-09_DP &
	   + T * ( -5.47904310e-13_DP ) ) ) ) ) -1.57203117e+04_DP )
      CP(SIXC3H5CH) =  1.18627154e+02_DP * ( &
	    6.27183793e-01_DP + T * (  4.66780254e-02_DP &
	   + T * ( -3.74430631e-05_DP + T * (  1.58330542e-08_DP &
	   + T * ( -2.73952155e-12_DP ) ) ) ) )
      H(STXC4H9O) =  1.13720593e+02_DP * ( &
	   T * ( -7.70464068e-01_DP + T * (  2.90963330e-02_DP &
	   + T * ( -1.45399307e-05_DP + T * (  4.37355548e-09_DP &
	   + T * ( -5.83672080e-13_DP ) ) ) ) ) -1.36502805e+04_DP )
      CP(STXC4H9O) =  1.13720593e+02_DP * ( &
	   -7.70464068e-01_DP + T * (  5.81926660e-02_DP &
	   + T * ( -4.36197921e-05_DP + T * (  1.74942219e-08_DP &
	   + T * ( -2.91836040e-12_DP ) ) ) ) )
      H(SIXC4H7O) =  1.16945257e+02_DP * ( &
	   T * (  1.74700687e+00_DP + T * (  2.03891718e-02_DP &
	   + T * ( -8.15834143e-06_DP + T * (  1.76625739e-09_DP &
	   + T * ( -1.50314118e-13_DP ) ) ) ) ) +  4.86979233e+03_DP )
      CP(SIXC4H7O) =  1.16945257e+02_DP * ( &
	    1.74700687e+00_DP + T * (  4.07783436e-02_DP &
	   + T * ( -2.44750243e-05_DP + T * (  7.06502958e-09_DP &
	   + T * ( -7.51570589e-13_DP ) ) ) ) )
      H(SC5H4CH2) =  1.06446715e+02_DP * ( &
	   T * ( -5.34007612e+00_DP + T * (  3.58641913e-02_DP &
	   + T * ( -2.15274819e-05_DP + T * (  6.96727892e-09_DP &
	   + T * ( -7.90002910e-13_DP ) ) ) ) ) +  2.58936616e+04_DP )
      CP(SC5H4CH2) =  1.06446715e+02_DP * ( &
	   -5.34007612e+00_DP + T * (  7.17283827e-02_DP &
	   + T * ( -6.45824457e-05_DP + T * (  2.78691157e-08_DP &
	   + T * ( -3.95001455e-12_DP ) ) ) ) )
      H(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   T * ( -4.87654845e+00_DP + T * (  3.13402891e-02_DP &
	   + T * ( -1.62467429e-05_DP + T * (  3.52805717e-09_DP &
	   + T * (  1.03703662e-13_DP ) ) ) ) ) +  3.99269438e+04_DP )
      CP(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   -4.87654845e+00_DP + T * (  6.26805782e-02_DP &
	   + T * ( -4.87402286e-05_DP + T * (  1.41122287e-08_DP &
	   + T * (  5.18518312e-13_DP ) ) ) ) )
      H(SA1C2H2XC) =  8.06153041e+01_DP * ( &
	   T * ( -6.30997035e+00_DP + T * (  4.75453915e-02_DP &
	   + T * ( -3.18566445e-05_DP + T * (  1.24202002e-08_DP &
	   + T * ( -2.03584362e-12_DP ) ) ) ) ) +  4.57329298e+04_DP )
      CP(SA1C2H2XC) =  8.06153041e+01_DP * ( &
	   -6.30997035e+00_DP + T * (  9.50907829e-02_DP &
	   + T * ( -9.55699336e-05_DP + T * (  4.96808010e-08_DP &
	   + T * ( -1.01792181e-11_DP ) ) ) ) )
      H(SA1C2H3XC) =  7.98350361e+01_DP * ( &
	   T * ( -5.38499941e+00_DP + T * (  4.10182578e-02_DP &
	   + T * ( -1.78153959e-05_DP + T * (  1.39773752e-09_DP &
	   + T * (  1.12227810e-12_DP ) ) ) ) ) +  1.60857559e+04_DP )
      CP(SA1C2H3XC) =  7.98350361e+01_DP * ( &
	   -5.38499941e+00_DP + T * (  8.20365155e-02_DP &
	   + T * ( -5.34461878e-05_DP + T * (  5.59095007e-09_DP &
	   + T * (  5.61139050e-12_DP ) ) ) ) )
      H(SA1C2HXC8) =  8.14109745e+01_DP * ( &
	   T * ( -5.21036925e+00_DP + T * (  4.32775972e-02_DP &
	   + T * ( -2.81669161e-05_DP + T * (  1.05480176e-08_DP &
	   + T * ( -1.63353233e-12_DP ) ) ) ) ) +  3.52488620e+04_DP )
      CP(SA1C2HXC8) =  8.14109745e+01_DP * ( &
	   -5.21036925e+00_DP + T * (  8.65551944e-02_DP &
	   + T * ( -8.45007483e-05_DP + T * (  4.21920706e-08_DP &
	   + T * ( -8.16766167e-12_DP ) ) ) ) )
      H(SA1C2HYXC) =  8.22225079e+01_DP * ( &
	   T * ( -4.42757639e+00_DP + T * (  4.18334322e-02_DP &
	   + T * ( -2.90035454e-05_DP + T * (  1.17571415e-08_DP &
	   + T * ( -2.03633970e-12_DP ) ) ) ) ) +  6.73302359e+04_DP )
      CP(SA1C2HYXC) =  8.22225079e+01_DP * ( &
	   -4.42757639e+00_DP + T * (  8.36668645e-02_DP &
	   + T * ( -8.70106362e-05_DP + T * (  4.70285661e-08_DP &
	   + T * ( -1.01816985e-11_DP ) ) ) ) )
      H(SA1C2H3YX) =  8.06153041e+01_DP * ( &
	   T * ( -5.36214520e+00_DP + T * (  4.33516648e-02_DP &
	   + T * ( -2.51432653e-05_DP + T * (  7.52849635e-09_DP &
	   + T * ( -6.81362836e-13_DP ) ) ) ) ) +  4.77818209e+04_DP )
      CP(SA1C2H3YX) =  8.06153041e+01_DP * ( &
	   -5.36214520e+00_DP + T * (  8.67033297e-02_DP &
	   + T * ( -7.54297960e-05_DP + T * (  3.01139854e-08_DP &
	   + T * ( -3.40681418e-12_DP ) ) ) ) )
      H(SA2XXC10H) =  6.53869263e+01_DP * ( &
	   T * ( -8.02718034e+00_DP + T * (  5.14622590e-02_DP &
	   + T * ( -2.78090670e-05_DP + T * (  6.80338458e-09_DP &
	   + T * ( -1.44911911e-13_DP ) ) ) ) ) +  5.01363344e+04_DP )
      CP(SA2XXC10H) =  6.53869263e+01_DP * ( &
	   -8.02718034e+00_DP + T * (  1.02924518e-01_DP &
	   + T * ( -8.34272010e-05_DP + T * (  2.72135383e-08_DP &
	   + T * ( -7.24559554e-13_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * ( -8.72434585e+00_DP + T * (  5.26880040e-02_DP &
	   + T * ( -2.67236897e-05_DP + T * (  5.46364935e-09_DP &
	   + T * (  2.84133212e-13_DP ) ) ) ) ) +  1.66588912e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   -8.72434585e+00_DP + T * (  1.05376008e-01_DP &
	   + T * ( -8.01710690e-05_DP + T * (  2.18545974e-08_DP &
	   + T * (  1.42066606e-12_DP ) ) ) ) )
      H(SA2YXC10H) =  6.53869263e+01_DP * ( &
	   T * ( -8.00768796e+00_DP + T * (  5.15206445e-02_DP &
	   + T * ( -2.79396999e-05_DP + T * (  6.91229315e-09_DP &
	   + T * ( -1.77768442e-13_DP ) ) ) ) ) +  4.99740633e+04_DP )
      CP(SA2YXC10H) =  6.53869263e+01_DP * ( &
	   -8.00768796e+00_DP + T * (  1.03041289e-01_DP &
	   + T * ( -8.38190998e-05_DP + T * (  2.76491726e-08_DP &
	   + T * ( -8.88842208e-13_DP ) ) ) ) )
      H(SA2C2H2AX) =  5.42739830e+01_DP * ( &
	   T * ( -9.26784872e+00_DP + T * (  6.75215650e-02_DP &
	   + T * ( -4.32639473e-05_DP + T * (  1.56805083e-08_DP &
	   + T * ( -2.32697154e-12_DP ) ) ) ) ) +  5.64832554e+04_DP )
      CP(SA2C2H2AX) =  5.42739830e+01_DP * ( &
	   -9.26784872e+00_DP + T * (  1.35043130e-01_DP &
	   + T * ( -1.29791842e-04_DP + T * (  6.27220331e-08_DP &
	   + T * ( -1.16348577e-11_DP ) ) ) ) )
      H(SA2C2H2BX) =  5.42739830e+01_DP * ( &
	   T * ( -9.38518818e+00_DP + T * (  6.74128980e-02_DP &
	   + T * ( -4.29874243e-05_DP + T * (  1.53620107e-08_DP &
	   + T * ( -2.18934312e-12_DP ) ) ) ) ) +  5.63724233e+04_DP )
      CP(SA2C2H2BX) =  5.42739830e+01_DP * ( &
	   -9.38518818e+00_DP + T * (  1.34825796e-01_DP &
	   + T * ( -1.28962273e-04_DP + T * (  6.14480428e-08_DP &
	   + T * ( -1.09467156e-11_DP ) ) ) ) )
      H(SA2C2HAXC) =  5.46334700e+01_DP * ( &
	   T * ( -8.23047877e+00_DP + T * (  6.30265880e-02_DP &
	   + T * ( -3.91663747e-05_DP + T * (  1.34197745e-08_DP &
	   + T * ( -1.77092462e-12_DP ) ) ) ) ) +  4.23747629e+04_DP )
      CP(SA2C2HAXC) =  5.46334700e+01_DP * ( &
	   -8.23047877e+00_DP + T * (  1.26053176e-01_DP &
	   + T * ( -1.17499124e-04_DP + T * (  5.36790980e-08_DP &
	   + T * ( -8.85462308e-12_DP ) ) ) ) )
      H(SA2C2HBXC) =  5.46334700e+01_DP * ( &
	   T * ( -8.22579974e+00_DP + T * (  6.31237755e-02_DP &
	   + T * ( -3.93802473e-05_DP + T * (  1.35996841e-08_DP &
	   + T * ( -1.82516934e-12_DP ) ) ) ) ) +  4.25495114e+04_DP )
      CP(SA2C2HBXC) =  5.46334700e+01_DP * ( &
	   -8.22579974e+00_DP + T * (  1.26247551e-01_DP &
	   + T * ( -1.18140742e-04_DP + T * (  5.43987363e-08_DP &
	   + T * ( -9.12584671e-12_DP ) ) ) ) )
      H(SA2C2HAYX) =  5.49977510e+01_DP * ( &
	   T * ( -7.36585075e+00_DP + T * (  6.12789955e-02_DP &
	   + T * ( -3.95790920e-05_DP + T * (  1.43485070e-08_DP &
	   + T * ( -2.10117262e-12_DP ) ) ) ) ) +  7.69836555e+04_DP )
      CP(SA2C2HAYX) =  5.49977510e+01_DP * ( &
	   -7.36585075e+00_DP + T * (  1.22557991e-01_DP &
	   + T * ( -1.18737276e-04_DP + T * (  5.73940282e-08_DP &
	   + T * ( -1.05058631e-11_DP ) ) ) ) )
      H(SA2C2HBYX) =  5.49977510e+01_DP * ( &
	   T * ( -7.35510706e+00_DP + T * (  6.12626680e-02_DP &
	   + T * ( -3.95804270e-05_DP + T * (  1.43577795e-08_DP &
	   + T * ( -2.10465530e-12_DP ) ) ) ) ) +  7.73354743e+04_DP )
      CP(SA2C2HBYX) =  5.49977510e+01_DP * ( &
	   -7.35510706e+00_DP + T * (  1.22525336e-01_DP &
	   + T * ( -1.18741281e-04_DP + T * (  5.74311178e-08_DP &
	   + T * ( -1.05232765e-11_DP ) ) ) ) )
      H(SA2R5XC12) =  5.46334700e+01_DP * ( &
	   T * ( -1.05497902e+01_DP + T * (  6.27683950e-02_DP &
	   + T * ( -3.45486817e-05_DP + T * (  8.82472825e-09_DP &
	   + T * ( -3.29016768e-13_DP ) ) ) ) ) +  2.94426605e+04_DP )
      CP(SA2R5XC12) =  5.46334700e+01_DP * ( &
	   -1.05497902e+01_DP + T * (  1.25536790e-01_DP &
	   + T * ( -1.03646045e-04_DP + T * (  3.52989130e-08_DP &
	   + T * ( -1.64508384e-12_DP ) ) ) ) )
      H(SA2R5XXC1) =  5.49977510e+01_DP * ( &
	   T * ( -9.79699017e+00_DP + T * (  6.11386065e-02_DP &
	   + T * ( -3.49775030e-05_DP + T * (  9.69866022e-09_DP &
	   + T * ( -6.32579322e-13_DP ) ) ) ) ) +  6.24840181e+04_DP )
      CP(SA2R5XXC1) =  5.49977510e+01_DP * ( &
	   -9.79699017e+00_DP + T * (  1.22277213e-01_DP &
	   + T * ( -1.04932509e-04_DP + T * (  3.87946409e-08_DP &
	   + T * ( -3.16289661e-12_DP ) ) ) ) )
      H(SA2R5C2H2) =  4.69174774e+01_DP * ( &
	   T * ( -9.79888742e+00_DP + T * (  7.27002575e-02_DP &
	   + T * ( -4.27866497e-05_DP + T * (  1.28967134e-08_DP &
	   + T * ( -1.20348472e-12_DP ) ) ) ) ) +  6.87094743e+04_DP )
      CP(SA2R5C2H2) =  4.69174774e+01_DP * ( &
	   -9.79888742e+00_DP + T * (  1.45400515e-01_DP &
	   + T * ( -1.28359949e-04_DP + T * (  5.15868534e-08_DP &
	   + T * ( -6.01742362e-12_DP ) ) ) ) )
      H(SA2R5C2HX) =  4.71858755e+01_DP * ( &
	   T * ( -9.95199604e+00_DP + T * (  7.27067810e-02_DP &
	   + T * ( -4.64002687e-05_DP + T * (  1.63758248e-08_DP &
	   + T * ( -2.27673192e-12_DP ) ) ) ) ) +  5.53993662e+04_DP )
      CP(SA2R5C2HX) =  4.71858755e+01_DP * ( &
	   -9.95199604e+00_DP + T * (  1.45413562e-01_DP &
	   + T * ( -1.39200806e-04_DP + T * (  6.55032994e-08_DP &
	   + T * ( -1.13836596e-11_DP ) ) ) ) )
      H(SA2R5C2HY) =  4.74573620e+01_DP * ( &
	   T * ( -9.09090029e+00_DP + T * (  7.09590820e-02_DP &
	   + T * ( -4.67904557e-05_DP + T * (  1.72745595e-08_DP &
	   + T * ( -2.59576030e-12_DP ) ) ) ) ) +  8.92211518e+04_DP )
      CP(SA2R5C2HY) =  4.74573620e+01_DP * ( &
	   -9.09090029e+00_DP + T * (  1.41918164e-01_DP &
	   + T * ( -1.40371367e-04_DP + T * (  6.90982381e-08_DP &
	   + T * ( -1.29788015e-11_DP ) ) ) ) )
      H(SP2XC12H1) =  5.39191958e+01_DP * ( &
	   T * ( -1.19438051e+01_DP + T * (  7.10815795e-02_DP &
	   + T * ( -4.44991497e-05_DP + T * (  1.55126429e-08_DP &
	   + T * ( -2.11533328e-12_DP ) ) ) ) ) +  2.01936932e+04_DP )
      CP(SP2XC12H1) =  5.39191958e+01_DP * ( &
	   -1.19438051e+01_DP + T * (  1.42163159e-01_DP &
	   + T * ( -1.33497449e-04_DP + T * (  6.20505718e-08_DP &
	   + T * ( -1.05766664e-11_DP ) ) ) ) )
      H(SP2XXC12H) =  5.42739830e+01_DP * ( &
	   T * ( -9.50091731e+00_DP + T * (  6.26052515e-02_DP &
	   + T * ( -3.27906074e-05_DP + T * (  6.89835438e-09_DP &
	   + T * (  3.35908690e-13_DP ) ) ) ) ) +  5.29903482e+04_DP )
      CP(SP2XXC12H) =  5.42739830e+01_DP * ( &
	   -9.50091731e+00_DP + T * (  1.25210503e-01_DP &
	   + T * ( -9.83718223e-05_DP + T * (  2.75934175e-08_DP &
	   + T * (  1.67954345e-12_DP ) ) ) ) )
      H(SA3XXC14H) =  4.69174774e+01_DP * ( &
	   T * ( -1.08881743e+01_DP + T * (  7.05935385e-02_DP &
	   + T * ( -3.78438093e-05_DP + T * (  8.97938392e-09_DP &
	   + T * ( -9.00424064e-14_DP ) ) ) ) ) +  5.70158539e+04_DP )
      CP(SA3XXC14H) =  4.69174774e+01_DP * ( &
	   -1.08881743e+01_DP + T * (  1.41187077e-01_DP &
	   + T * ( -1.13531428e-04_DP + T * (  3.59175357e-08_DP &
	   + T * ( -4.50212032e-13_DP ) ) ) ) )
      H(SA3XC14H1) =  4.66521154e+01_DP * ( &
	   T * ( -1.15461369e+01_DP + T * (  7.18790815e-02_DP &
	   + T * ( -3.69563747e-05_DP + T * (  7.80450353e-09_DP &
	   + T * (  2.91950464e-13_DP ) ) ) ) ) +  2.21687904e+04_DP )
      CP(SA3XC14H1) =  4.66521154e+01_DP * ( &
	   -1.15461369e+01_DP + T * (  1.43758163e-01_DP &
	   + T * ( -1.10869124e-04_DP + T * (  3.12180141e-08_DP &
	   + T * (  1.45975232e-12_DP ) ) ) ) )
      H(SA3YXC14H) =  4.69174774e+01_DP * ( &
	   T * ( -1.08881743e+01_DP + T * (  7.05935385e-02_DP &
	   + T * ( -3.78438093e-05_DP + T * (  8.97938392e-09_DP &
	   + T * ( -9.00424064e-14_DP ) ) ) ) ) +  5.70158539e+04_DP )
      CP(SA3YXC14H) =  4.69174774e+01_DP * ( &
	   -1.08881743e+01_DP + T * (  1.41187077e-01_DP &
	   + T * ( -1.13531428e-04_DP + T * (  3.59175357e-08_DP &
	   + T * ( -4.50212032e-13_DP ) ) ) ) )
      H(SA3R5XXC1) =  4.13171861e+01_DP * ( &
	   T * ( -1.25419370e+01_DP + T * (  8.01855080e-02_DP &
	   + T * ( -4.50096677e-05_DP + T * (  1.19011023e-08_DP &
	   + T * ( -5.87324606e-13_DP ) ) ) ) ) +  6.87107686e+04_DP )
      CP(SA3R5XXC1) =  4.13171861e+01_DP * ( &
	   -1.25419370e+01_DP + T * (  1.60371016e-01_DP &
	   + T * ( -1.35029003e-04_DP + T * (  4.76044093e-08_DP &
	   + T * ( -2.93662303e-12_DP ) ) ) ) )
      H(SA3R5XC16) =  4.11112540e+01_DP * ( &
	   T * ( -1.32241574e+01_DP + T * (  8.14313235e-02_DP &
	   + T * ( -4.39890187e-05_DP + T * (  1.06098466e-08_DP &
	   + T * ( -1.70634345e-13_DP ) ) ) ) ) +  3.44430693e+04_DP )
      CP(SA3R5XC16) =  4.11112540e+01_DP * ( &
	   -1.32241574e+01_DP + T * (  1.62862647e-01_DP &
	   + T * ( -1.31967056e-04_DP + T * (  4.24393865e-08_DP &
	   + T * ( -8.53171724e-13_DP ) ) ) ) )
      H(SA4XC16H1) =  4.11112540e+01_DP * ( &
	   T * ( -1.31524443e+01_DP + T * (  8.04394215e-02_DP &
	   + T * ( -4.25732390e-05_DP + T * (  9.77297245e-09_DP &
	   + T * (  1.48798225e-14_DP ) ) ) ) ) +  2.49673872e+04_DP )
      CP(SA4XC16H1) =  4.11112540e+01_DP * ( &
	   -1.31524443e+01_DP + T * (  1.60878843e-01_DP &
	   + T * ( -1.27719717e-04_DP + T * (  3.90918898e-08_DP &
	   + T * (  7.43991125e-14_DP ) ) ) ) )
      H(SA4XXC16H) =  4.13171861e+01_DP * ( &
	   T * ( -1.23671835e+01_DP + T * (  7.88287575e-02_DP &
	   + T * ( -4.30812110e-05_DP + T * (  1.07157418e-08_DP &
	   + T * ( -3.08983758e-13_DP ) ) ) ) ) +  6.27797890e+04_DP )
      CP(SA4XXC16H) =  4.13171861e+01_DP * ( &
	   -1.23671835e+01_DP + T * (  1.57657515e-01_DP &
	   + T * ( -1.29243633e-04_DP + T * (  4.28629673e-08_DP &
	   + T * ( -1.54491879e-12_DP ) ) ) ) )
      H(SA4R5XC18) =  3.67468399e+01_DP * ( &
	   T * ( -1.47695663e+01_DP + T * (  8.98290110e-02_DP &
	   + T * ( -4.93966333e-05_DP + T * (  1.24430526e-08_DP &
	   + T * ( -4.12874670e-13_DP ) ) ) ) ) +  3.78467972e+04_DP )
      CP(SA4R5XC18) =  3.67468399e+01_DP * ( &
	   -1.47695663e+01_DP + T * (  1.79658022e-01_DP &
	   + T * ( -1.48189900e-04_DP + T * (  4.97722102e-08_DP &
	   + T * ( -2.06437335e-12_DP ) ) ) ) )
      H(SFLTNXC16) =  4.11112540e+01_DP * ( &
	   T * ( -1.29396091e+01_DP + T * (  8.01598610e-02_DP &
	   + T * ( -4.21507513e-05_DP + T * (  9.40122770e-09_DP &
	   + T * (  1.36967841e-13_DP ) ) ) ) ) +  2.91084888e+04_DP )
      CP(SFLTNXC16) =  4.11112540e+01_DP * ( &
	   -1.29396091e+01_DP + T * (  1.60319722e-01_DP &
	   + T * ( -1.26452254e-04_DP + T * (  3.76049108e-08_DP &
	   + T * (  6.84839205e-13_DP ) ) ) ) )
      H(SC5H6) =  1.25788072e+02_DP * ( &
	   T * ( -5.13691194e+00_DP + T * (  3.03476727e-02_DP &
	   + T * ( -1.53517612e-05_DP + T * (  3.21143002e-09_DP &
	   + T * (  1.48242970e-13_DP ) ) ) ) ) +  1.53675713e+04_DP )
      CP(SC5H6) =  1.25788072e+02_DP * ( &
	   -5.13691194e+00_DP + T * (  6.06953453e-02_DP &
	   + T * ( -4.60552837e-05_DP + T * (  1.28457201e-08_DP &
	   + T * (  7.41214852e-13_DP ) ) ) ) )
      H(SC5H5) =  1.27736058e+02_DP * ( &
	   T * ( -7.37844042e+00_DP + T * (  4.86195909e-02_DP &
	   + T * ( -5.65263793e-05_DP + T * (  3.79546668e-08_DP &
	   + T * ( -1.02415096e-11_DP ) ) ) ) ) +  3.05514662e+04_DP )
      CP(SC5H5) =  1.27736058e+02_DP * ( &
	   -7.37844042e+00_DP + T * (  9.72391818e-02_DP &
	   + T * ( -1.69579138e-04_DP + T * (  1.51818667e-07_DP &
	   + T * ( -5.12075479e-11_DP ) ) ) ) )
      H(STXC5H5O) =  1.02532248e+02_DP * ( &
	   T * (  2.30436010e-01_DP + T * (  1.61612860e-02_DP &
	   + T * (  9.63363600e-06_DP + T * ( -1.76701532e-08_DP &
	   + T * (  6.68143480e-12_DP ) ) ) ) ) +  5.55547240e+03_DP )
      CP(STXC5H5O) =  1.02532248e+02_DP * ( &
	    2.30436010e-01_DP + T * (  3.23225720e-02_DP &
	   + T * (  2.89009080e-05_DP + T * ( -7.06806130e-08_DP &
	   + T * (  3.34071740e-11_DP ) ) ) ) )
      H(SC5H4O) =  1.03822832e+02_DP * ( &
	   T * ( -3.64380971e+00_DP + T * (  3.07164598e-02_DP &
	   + T * ( -1.97383079e-05_DP + T * (  7.08083390e-09_DP &
	   + T * ( -1.00545227e-12_DP ) ) ) ) ) +  5.46809680e+03_DP )
      CP(SC5H4O) =  1.03822832e+02_DP * ( &
	   -3.64380971e+00_DP + T * (  6.14329196e-02_DP &
	   + T * ( -5.92149236e-05_DP + T * (  2.83233356e-08_DP &
	   + T * ( -5.02726134e-12_DP ) ) ) ) )
      H(SSXC5H5O) =  1.02532248e+02_DP * ( &
	   T * ( -3.07776000e+00_DP + T * (  2.62908395e-02_DP &
	   + T * ( -9.61883767e-06_DP + T * ( -8.47136975e-10_DP &
	   + T * (  1.26722798e-12_DP ) ) ) ) ) +  2.55104550e+04_DP )
      CP(SSXC5H5O) =  1.02532248e+02_DP * ( &
	   -3.07776000e+00_DP + T * (  5.25816790e-02_DP &
	   + T * ( -2.88565130e-05_DP + T * ( -3.38854790e-09_DP &
	   + T * (  6.33613990e-12_DP ) ) ) ) )
      H(SC9H8) =  7.15803158e+01_DP * ( &
	   T * ( -8.12447817e+00_DP + T * (  4.88828533e-02_DP &
	   + T * ( -2.43478658e-05_DP + T * (  4.70737525e-09_DP &
	   + T * (  3.68066426e-13_DP ) ) ) ) ) +  1.86589996e+04_DP )
      CP(SC9H8) =  7.15803158e+01_DP * ( &
	   -8.12447817e+00_DP + T * (  9.77657067e-02_DP &
	   + T * ( -7.30435974e-05_DP + T * (  1.88295010e-08_DP &
	   + T * (  1.84033213e-12_DP ) ) ) ) )
      H(SC9H7) =  7.22069373e+01_DP * ( &
	   T * ( -8.73685384e+00_DP + T * (  5.17108180e-02_DP &
	   + T * ( -3.07807798e-05_DP + T * (  9.39057395e-09_DP &
	   + T * ( -8.81210540e-13_DP ) ) ) ) ) +  3.31641009e+04_DP )
      CP(SC9H7) =  7.22069373e+01_DP * ( &
	   -8.73685384e+00_DP + T * (  1.03421636e-01_DP &
	   + T * ( -9.23423393e-05_DP + T * (  3.75622958e-08_DP &
	   + T * ( -4.40605270e-12_DP ) ) ) ) )
      H(SA1CH2XC7) =  9.12400413e+01_DP * ( &
	   T * ( -6.07053038e+00_DP + T * (  4.17600754e-02_DP &
	   + T * ( -2.47233361e-05_DP + T * (  7.82884618e-09_DP &
	   + T * ( -8.47341736e-13_DP ) ) ) ) ) +  2.35894712e+04_DP )
      CP(SA1CH2XC7) =  9.12400413e+01_DP * ( &
	   -6.07053038e+00_DP + T * (  8.35201507e-02_DP &
	   + T * ( -7.41700083e-05_DP + T * (  3.13153847e-08_DP &
	   + T * ( -4.23670868e-12_DP ) ) ) ) )
      H(SC9H6O) =  6.38886413e+01_DP * ( &
	   T * ( -6.53928778e+00_DP + T * (  4.84661643e-02_DP &
	   + T * ( -2.72566219e-05_DP + T * (  7.41748685e-09_DP &
	   + T * ( -4.49986784e-13_DP ) ) ) ) ) +  6.88883578e+03_DP )
      CP(SC9H6O) =  6.38886413e+01_DP * ( &
	   -6.53928778e+00_DP + T * (  9.69323286e-02_DP &
	   + T * ( -8.17698656e-05_DP + T * (  2.96699474e-08_DP &
	   + T * ( -2.24993392e-12_DP ) ) ) ) )
      H(SOXC6H4) =  1.09266940e+02_DP * ( &
	   T * ( -3.46229657e+00_DP + T * (  2.87008288e-02_DP &
	   + T * ( -1.64328123e-05_DP + T * (  4.76701207e-09_DP &
	   + T * ( -3.90861418e-13_DP ) ) ) ) ) +  5.25223614e+04_DP )
      CP(SOXC6H4) =  1.09266940e+02_DP * ( &
	   -3.46229657e+00_DP + T * (  5.74016575e-02_DP &
	   + T * ( -4.92984369e-05_DP + T * (  1.90680483e-08_DP &
	   + T * ( -1.95430709e-12_DP ) ) ) ) )
      H(SA1CH3XC7) =  9.02418217e+01_DP * ( &
	   T * ( -4.54072038e+00_DP + T * (  3.42713573e-02_DP &
	   + T * ( -1.19037675e-05_DP + T * ( -1.04849410e-09_DP &
	   + T * (  1.48355959e-12_DP ) ) ) ) ) +  4.64121087e+03_DP )
      CP(SA1CH3XC7) =  9.02418217e+01_DP * ( &
	   -4.54072038e+00_DP + T * (  6.85427145e-02_DP &
	   + T * ( -3.57113024e-05_DP + T * ( -4.19397642e-09_DP &
	   + T * (  7.41779795e-12_DP ) ) ) ) )
      H(SA1OHXC6H) =  8.83489183e+01_DP * ( &
	   T * ( -3.56571190e+00_DP + T * (  3.30067717e-02_DP &
	   + T * ( -1.30985939e-05_DP + T * ( -9.05634825e-10_DP &
	   + T * (  1.72483122e-12_DP ) ) ) ) ) -1.31101467e+04_DP )
      CP(SA1OHXC6H) =  8.83489183e+01_DP * ( &
	   -3.56571190e+00_DP + T * (  6.60135435e-02_DP &
	   + T * ( -3.92957818e-05_DP + T * ( -3.62253930e-09_DP &
	   + T * (  8.62415610e-12_DP ) ) ) ) )
      H(SHOA1CH3X) =  7.68892300e+01_DP * ( &
	   T * ( -2.49882920e+00_DP + T * (  3.45763628e-02_DP &
	   + T * ( -9.60285530e-06_DP + T * ( -4.32281625e-09_DP &
	   + T * (  2.72653336e-12_DP ) ) ) ) ) -1.76336196e+04_DP )
      CP(SHOA1CH3X) =  7.68892300e+01_DP * ( &
	   -2.49882920e+00_DP + T * (  6.91527256e-02_DP &
	   + T * ( -2.88085659e-05_DP + T * ( -1.72912650e-08_DP &
	   + T * (  1.36326668e-11_DP ) ) ) ) )
      H(SOA1CH3XC) =  7.76127177e+01_DP * ( &
	   T * ( -3.88641950e+00_DP + T * (  3.92129808e-02_DP &
	   + T * ( -2.06935634e-05_DP + T * (  5.42917367e-09_DP &
	   + T * ( -3.31557060e-13_DP ) ) ) ) ) +  4.71681386e+02_DP )
      CP(SOA1CH3XC) =  7.76127177e+01_DP * ( &
	   -3.88641950e+00_DP + T * (  7.84259616e-02_DP &
	   + T * ( -6.20806903e-05_DP + T * (  2.17166947e-08_DP &
	   + T * ( -1.65778530e-12_DP ) ) ) ) )
      H(SA1CH2OXC) =  7.76127177e+01_DP * ( &
	   T * ( -4.75332952e+00_DP + T * (  3.92916108e-02_DP &
	   + T * ( -1.85706249e-05_DP + T * (  2.99161712e-09_DP &
	   + T * (  5.03808950e-13_DP ) ) ) ) ) +  1.31155256e+04_DP )
      CP(SA1CH2OXC) =  7.76127177e+01_DP * ( &
	   -4.75332952e+00_DP + T * (  7.85832217e-02_DP &
	   + T * ( -5.57118748e-05_DP + T * (  1.19664685e-08_DP &
	   + T * (  2.51904475e-12_DP ) ) ) ) )
      H(SA1CH2OHX) =  7.68892300e+01_DP * ( &
	   T * (  2.85739935e+00_DP + T * (  1.19385310e-02_DP &
	   + T * (  2.80169339e-05_DP + T * ( -3.34963965e-08_DP &
	   + T * (  1.13987775e-11_DP ) ) ) ) ) -1.37956049e+04_DP )
      CP(SA1CH2OHX) =  7.68892300e+01_DP * ( &
	    2.85739935e+00_DP + T * (  2.38770620e-02_DP &
	   + T * (  8.40508017e-05_DP + T * ( -1.33985586e-07_DP &
	   + T * (  5.69938876e-11_DP ) ) ) ) )
      H(SA1CHOXC7) =  7.83499501e+01_DP * ( &
	   T * ( -3.47171048e+00_DP + T * (  3.46445945e-02_DP &
	   + T * ( -1.44201170e-05_DP + T * (  8.59677740e-10_DP &
	   + T * (  9.62020522e-13_DP ) ) ) ) ) -6.14558774e+03_DP )
      CP(SA1CHOXC7) =  7.83499501e+01_DP * ( &
	   -3.47171048e+00_DP + T * (  6.92891889e-02_DP &
	   + T * ( -4.32603509e-05_DP + T * (  3.43871096e-09_DP &
	   + T * (  4.81010261e-12_DP ) ) ) ) )
      H(SA1OXC6H5) =  8.93054780e+01_DP * ( &
	   T * ( -4.51502441e+00_DP + T * (  3.51975764e-02_DP &
	   + T * ( -1.97518737e-05_DP + T * (  5.32274058e-09_DP &
	   + T * ( -2.90337352e-13_DP ) ) ) ) ) +  5.19173466e+03_DP )
      CP(SA1OXC6H5) =  8.93054780e+01_DP * ( &
	   -4.51502441e+00_DP + T * (  7.03951529e-02_DP &
	   + T * ( -5.92556211e-05_DP + T * (  2.12909623e-08_DP &
	   + T * ( -1.45168676e-12_DP ) ) ) ) )
      H(SA1CH3YXC) =  9.12400413e+01_DP * ( &
	   T * ( -3.91657299e+00_DP + T * (  3.32935047e-02_DP &
	   + T * ( -1.33241679e-05_DP + T * (  5.11701993e-10_DP &
	   + T * (  9.97118780e-13_DP ) ) ) ) ) +  3.54243469e+04_DP )
      CP(SA1CH3YXC) =  9.12400413e+01_DP * ( &
	   -3.91657299e+00_DP + T * (  6.65870094e-02_DP &
	   + T * ( -3.99725038e-05_DP + T * (  2.04680797e-09_DP &
	   + T * (  4.98559390e-12_DP ) ) ) ) )
      H(SA1C2H4XC) =  7.90697276e+01_DP * ( &
	   T * (  7.33299107e-01_DP + T * (  2.29526579e-02_DP &
	   + T * (  1.26085744e-05_DP + T * ( -2.28091853e-08_DP &
	   + T * (  8.51179356e-12_DP ) ) ) ) ) +  2.61572945e+04_DP )
      CP(SA1C2H4XC) =  7.90697276e+01_DP * ( &
	    7.33299107e-01_DP + T * (  4.59053158e-02_DP &
	   + T * (  3.78257231e-05_DP + T * ( -9.12367411e-08_DP &
	   + T * (  4.25589678e-11_DP ) ) ) ) )
      H(SA1C2H5XC) =  7.83189525e+01_DP * ( &
	   T * (  1.24076722e+00_DP + T * (  1.79566415e-02_DP &
	   + T * (  2.51407491e-05_DP + T * ( -3.29760752e-08_DP &
	   + T * (  1.14949361e-11_DP ) ) ) ) ) +  1.18391719e+03_DP )
      CP(SA1C2H5XC) =  7.83189525e+01_DP * ( &
	    1.24076722e+00_DP + T * (  3.59132829e-02_DP &
	   + T * (  7.54222474e-05_DP + T * ( -1.31904301e-07_DP &
	   + T * (  5.74746803e-11_DP ) ) ) ) )
      H(SC8H9O2) =  6.06213544e+01_DP * ( &
	   T * ( -3.87929417e+00_DP + T * (  4.49066033e-02_DP &
	   + T * ( -1.64561836e-05_DP + T * ( -1.01586420e-09_DP &
	   + T * (  1.91431787e-12_DP ) ) ) ) ) +  1.81922098e+04_DP )
      CP(SC8H9O2) =  6.06213544e+01_DP * ( &
	   -3.87929417e+00_DP + T * (  8.98132066e-02_DP &
	   + T * ( -4.93685508e-05_DP + T * ( -4.06345681e-09_DP &
	   + T * (  9.57158935e-12_DP ) ) ) ) )
      H(SC8H8OOH) =  6.06213544e+01_DP * ( &
	   T * ( -3.87929417e+00_DP + T * (  4.49066033e-02_DP &
	   + T * ( -1.64561836e-05_DP + T * ( -1.01586420e-09_DP &
	   + T * (  1.91431787e-12_DP ) ) ) ) ) +  1.81922098e+04_DP )
      CP(SC8H8OOH) =  6.06213544e+01_DP * ( &
	   -3.87929417e+00_DP + T * (  8.98132066e-02_DP &
	   + T * ( -4.93685508e-05_DP + T * ( -4.06345681e-09_DP &
	   + T * (  9.57158935e-12_DP ) ) ) ) )
      H(SOC8H7OOH) =  5.46478336e+01_DP * ( &
	   T * ( -1.93576547e+00_DP + T * (  4.62450358e-02_DP &
	   + T * ( -2.03243198e-05_DP + T * (  1.97462569e-09_DP &
	   + T * (  1.08205393e-12_DP ) ) ) ) ) -2.01714879e+04_DP )
      CP(SOC8H7OOH) =  5.46478336e+01_DP * ( &
	   -1.93576547e+00_DP + T * (  9.24900716e-02_DP &
	   + T * ( -6.09729593e-05_DP + T * (  7.89850275e-09_DP &
	   + T * (  5.41026965e-12_DP ) ) ) ) )
      H(SA1CH3CH3) =  7.83189525e+01_DP * ( &
	   T * ( -3.46066830e+00_DP + T * (  3.58894658e-02_DP &
	   + T * ( -8.52036773e-06_DP + T * ( -4.37176937e-09_DP &
	   + T * (  2.45713912e-12_DP ) ) ) ) ) +  1.62314629e+02_DP )
      CP(SA1CH3CH3) =  7.83189525e+01_DP * ( &
	   -3.46066830e+00_DP + T * (  7.17789316e-02_DP &
	   + T * ( -2.55611032e-05_DP + T * ( -1.74870775e-08_DP &
	   + T * (  1.22856956e-11_DP ) ) ) ) )
      H(SA1CH3CH2) =  7.90697276e+01_DP * ( &
	   T * ( -5.06171538e+00_DP + T * (  4.35354132e-02_DP &
	   + T * ( -2.15249935e-05_DP + T * (  4.61621605e-09_DP &
	   + T * (  9.95743576e-14_DP ) ) ) ) ) +  1.90750247e+04_DP )
      CP(SA1CH3CH2) =  7.90697276e+01_DP * ( &
	   -5.06171538e+00_DP + T * (  8.70708264e-02_DP &
	   + T * ( -6.45749806e-05_DP + T * (  1.84648642e-08_DP &
	   + T * (  4.97871788e-13_DP ) ) ) ) )
      H(SA1CH3CHO) =  6.92031229e+01_DP * ( &
	   T * ( -2.28640538e+00_DP + T * (  3.60390418e-02_DP &
	   + T * ( -1.07859591e-05_DP + T * ( -2.61287950e-09_DP &
	   + T * (  1.97213191e-12_DP ) ) ) ) ) -1.07389561e+04_DP )
      CP(SA1CH3CHO) =  6.92031229e+01_DP * ( &
	   -2.28640538e+00_DP + T * (  7.20780836e-02_DP &
	   + T * ( -3.23578774e-05_DP + T * ( -1.04515180e-08_DP &
	   + T * (  9.86065956e-12_DP ) ) ) ) )
      H(SA2CH3XC1) =  5.84734510e+01_DP * ( &
	   T * ( -7.20788596e+00_DP + T * (  5.46223345e-02_DP &
	   + T * ( -2.37103682e-05_DP + T * (  2.35931837e-09_DP &
	   + T * (  1.20563270e-12_DP ) ) ) ) ) +  1.21406213e+04_DP )
      CP(SA2CH3XC1) =  5.84734510e+01_DP * ( &
	   -7.20788596e+00_DP + T * (  1.09244669e-01_DP &
	   + T * ( -7.11311046e-05_DP + T * (  9.43727348e-09_DP &
	   + T * (  6.02816350e-12_DP ) ) ) ) )
      H(SA1CHOCH2) =  6.97886449e+01_DP * ( &
	   T * ( -3.91654934e+00_DP + T * (  4.37771316e-02_DP &
	   + T * ( -2.39035516e-05_DP + T * (  6.44161300e-09_DP &
	   + T * ( -4.00917698e-13_DP ) ) ) ) ) +  8.48223476e+03_DP )
      CP(SA1CHOCH2) =  6.97886449e+01_DP * ( &
	   -3.91654934e+00_DP + T * (  8.75542631e-02_DP &
	   + T * ( -7.17106548e-05_DP + T * (  2.57664520e-08_DP &
	   + T * ( -2.00458849e-12_DP ) ) ) ) )
      H(SA1CHOCHO) =  6.19881009e+01_DP * ( &
	   T * ( -1.05921174e+00_DP + T * (  3.62775300e-02_DP &
	   + T * ( -1.32481793e-05_DP + T * ( -7.00683277e-10_DP &
	   + T * (  1.44399999e-12_DP ) ) ) ) ) -2.09983131e+04_DP )
      CP(SA1CHOCHO) =  6.19881009e+01_DP * ( &
	   -1.05921174e+00_DP + T * (  7.25550600e-02_DP &
	   + T * ( -3.97445378e-05_DP + T * ( -2.80273311e-09_DP &
	   + T * (  7.21999995e-12_DP ) ) ) ) )
      H(SA2OHXC10) =  5.76727893e+01_DP * ( &
	   T * ( -2.08768263e+00_DP + T * (  3.84049753e-02_DP &
	   + T * ( -5.11976743e-06_DP + T * ( -1.01164408e-08_DP &
	   + T * (  4.67519558e-12_DP ) ) ) ) ) -6.29056385e+03_DP )
      CP(SA2OHXC10) =  5.76727893e+01_DP * ( &
	   -2.08768263e+00_DP + T * (  7.68099506e-02_DP &
	   + T * ( -1.53593023e-05_DP + T * ( -4.04657632e-08_DP &
	   + T * (  2.33759779e-11_DP ) ) ) ) )
      H(SA2CH2XC1) =  5.88909351e+01_DP * ( &
	   T * ( -9.33279584e+00_DP + T * (  6.19204430e-02_DP &
	   + T * ( -3.60778607e-05_DP + T * (  1.08420439e-08_DP &
	   + T * ( -1.00357549e-12_DP ) ) ) ) ) +  3.25141112e+04_DP )
      CP(SA2CH2XC1) =  5.88909351e+01_DP * ( &
	   -9.33279584e+00_DP + T * (  1.23840886e-01_DP &
	   + T * ( -1.08233582e-04_DP + T * (  4.33681755e-08_DP &
	   + T * ( -5.01787746e-12_DP ) ) ) ) )
      H(SA2CH2OXC) =  5.28962604e+01_DP * ( &
	   T * ( -7.78907165e+00_DP + T * (  5.91291930e-02_DP &
	   + T * ( -2.97536620e-05_DP + T * (  6.00588557e-09_DP &
	   + T * (  3.26602332e-13_DP ) ) ) ) ) +  2.40306202e+04_DP )
      CP(SA2CH2OXC) =  5.28962604e+01_DP * ( &
	   -7.78907165e+00_DP + T * (  1.18258386e-01_DP &
	   + T * ( -8.92609860e-05_DP + T * (  2.40235423e-08_DP &
	   + T * (  1.63301166e-12_DP ) ) ) ) )
      H(SA2CHOXC1) =  5.32376708e+01_DP * ( &
	   T * ( -4.87929110e+00_DP + T * (  4.91546815e-02_DP &
	   + T * ( -1.51615113e-05_DP + T * ( -5.04520860e-09_DP &
	   + T * (  3.54298268e-12_DP ) ) ) ) ) +  1.63189711e+03_DP )
      CP(SA2CHOXC1) =  5.32376708e+01_DP * ( &
	   -4.87929110e+00_DP + T * (  9.83093630e-02_DP &
	   + T * ( -4.54845339e-05_DP + T * ( -2.01808344e-08_DP &
	   + T * (  1.77149134e-11_DP ) ) ) ) )
      H(SA2OXC10H) =  5.80788790e+01_DP * ( &
	   T * ( -1.15176448e+00_DP + T * (  3.05677256e-02_DP &
	   + T * (  1.06717028e-05_DP + T * ( -2.48571322e-08_DP &
	   + T * (  9.59980086e-12_DP ) ) ) ) ) +  1.14058756e+04_DP )
      CP(SA2OXC10H) =  5.80788790e+01_DP * ( &
	   -1.15176448e+00_DP + T * (  6.11354512e-02_DP &
	   + T * (  3.20151083e-05_DP + T * ( -9.94285290e-08_DP &
	   + T * (  4.79990043e-11_DP ) ) ) ) )
      H(SOC6H4O) =  7.69191059e+01_DP * ( &
	   T * ( -2.04371804e+00_DP + T * (  3.30482234e-02_DP &
	   + T * ( -1.89325813e-05_DP + T * (  5.62252577e-09_DP &
	   + T * ( -5.34699342e-13_DP ) ) ) ) ) -1.24369410e+04_DP )
      CP(SOC6H4O) =  7.69191059e+01_DP * ( &
	   -2.04371804e+00_DP + T * (  6.60964467e-02_DP &
	   + T * ( -5.67977439e-05_DP + T * (  2.24901031e-08_DP &
	   + T * ( -2.67349671e-12_DP ) ) ) ) )
      END IF

      END
