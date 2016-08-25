!----------------------------------------------------------
! ======= EthyRedF.f90 =======
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
      include 'EthyRedF90.h'
      real(DP) :: CDOT(45), W(279), K(279), &
      C(45), M(7), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM1) = C(SN2) + C(SH) &
	    + C(SO) + C(SOH) &
	    + C(SH2O2) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + C(SC2H2) + C(SC2H3) &
	    + C(SC2H4) + 2 * C(SCH4) &
	    + C(SC2H5) + C(SCH2CO) &
	    + C(SHCCO) + C(SCH2OH) &
	    + C(SC2H) + C(SC3H3) &
	    + C(SNXC4H3) + C(SIXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4) &
	    + C(SA1XC6H6) + C(SA1XXC6H5)
      M(MM1) = M(MM1) + C(SC5H5) + C(SC3H6) &
	    + C(SC4H2) + C(SOXC6H4) &
	    + C(SA2XC10H8) + C(SC9H8) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6) + C(SC9H7)
      M(MM3) = C(SN2) + C(SH) &
	    + 2 * C(SH2) + C(SO) &
	    + C(SOH) + 12 * C(SH2O) &
	    + C(SH2O2) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + 3.6 * C(SCO2) + C(SC2H2) &
	    + C(SC2H3) + C(SC2H4) &
	    + 2 * C(SCH4) + C(SC2H5) &
	    + C(SCH2CO) + C(SHCCO) &
	    + C(SCH2OH) + C(SC2H) &
	    + C(SC3H3) + C(SNXC4H3) &
	    + C(SIXC3H7) + 3 * C(SC2H6)
      M(MM3) = M(MM3) + C(SPXC3H4) + C(SA1XC6H6) &
	    + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H2) &
	    + C(SOXC6H4) + C(SA2XC10H8) &
	    + C(SC9H8) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6) &
	    + C(SC9H7)
      M(MM2) = C(SN2) + C(SH) &
	    + 2 * C(SH2) + C(SO) &
	    + C(SOH) + 6.3 * C(SH2O) &
	    + C(SH2O2) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + 3.6 * C(SCO2) + C(SC2H2) &
	    + C(SC2H3) + C(SC2H4) &
	    + 2 * C(SCH4) + C(SC2H5) &
	    + C(SCH2CO) + C(SHCCO) &
	    + C(SCH2OH) + C(SC2H) &
	    + C(SC3H3) + C(SNXC4H3) &
	    + C(SIXC3H7) + 3 * C(SC2H6)
      M(MM2) = M(MM2) + C(SPXC3H4) + C(SA1XC6H6) &
	    + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H2) &
	    + C(SOXC6H4) + C(SA2XC10H8) &
	    + C(SC9H8) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6) &
	    + C(SC9H7)
      M(MM5) = C(SN2) + C(SH) &
	    + 0.75 * C(SH2) + C(SO) &
	    + C(SOH) + 11.89 * C(SH2O) &
	    + C(SH2O2) + 0.85 * C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.09 * C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + 2.18 * C(SCO2) + C(SC2H2) &
	    + C(SC2H3) + C(SC2H4) &
	    + C(SCH4) + C(SC2H5) &
	    + C(SCH2CO) + C(SHCCO) &
	    + C(SCH2OH) + C(SC2H) &
	    + C(SC3H3) + C(SNXC4H3) &
	    + C(SIXC3H7) + C(SC2H6)
      M(MM5) = M(MM5) + C(SPXC3H4) + C(SA1XC6H6) &
	    + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H2) &
	    + C(SOXC6H4) + C(SA2XC10H8) &
	    + C(SC9H8) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6) &
	    + C(SC9H7)
      M(MM9) = C(SN2) + C(SH) &
	    + 2 * C(SH2) + C(SO) &
	    + C(SOH) + 6 * C(SH2O) &
	    + C(SH2O2) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.5 * C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + 2 * C(SCO2) + C(SC2H2) &
	    + C(SC2H3) + C(SC2H4) &
	    + 3 * C(SCH4) + C(SC2H5) &
	    + C(SCH2CO) + C(SHCCO) &
	    + C(SCH2OH) + C(SC2H) &
	    + C(SC3H3) + C(SNXC4H3) &
	    + C(SIXC3H7) + 3 * C(SC2H6)
      M(MM9) = M(MM9) + C(SPXC3H4) + C(SA1XC6H6) &
	    + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H2) &
	    + C(SOXC6H4) + C(SA2XC10H8) &
	    + C(SC9H8) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6) &
	    + C(SC9H7)
      M(MM7) = C(SN2) + C(SH) &
	    + 2 * C(SH2) + C(SO) &
	    + C(SOH) + 12 * C(SH2O) &
	    + C(SH2O2) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(STXCH2) + C(SCH3) &
	    + C(SCH2O) + C(SSXCH2) &
	    + 3.6 * C(SCO2) + C(SC2H2) &
	    + C(SC2H3) + C(SC2H4) &
	    + 2 * C(SCH4) + C(SC2H5) &
	    + C(SCH2CO) + C(SHCCO) &
	    + C(SCH2OH) + C(SC2H) &
	    + C(SC3H3) + C(SNXC4H3) &
	    + C(SIXC3H7) + 3 * C(SC2H6)
      M(MM7) = M(MM7) + C(SPXC3H4) + C(SA1XC6H6) &
	    + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H2) &
	    + C(SOXC6H4) + C(SA2XC10H8) &
	    + C(SC9H8) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6) &
	    + C(SC9H7)
      M(MM8) = C(SN2) + C(SH) &
	    + 2 * C(SH2) + C(SO) &
	    + C(SOH) + C(SH2O2) &
	    + C(SO2) + C(SHO2) &
	    + C(SCH) + 1.75 * C(SCO) &
	    + C(SHCO) + C(STXCH2) &
	    + C(SCH3) + C(SCH2O) &
	    + C(SSXCH2) + 3.6 * C(SCO2) &
	    + C(SC2H2) + C(SC2H3) &
	    + C(SC2H4) + 2 * C(SCH4) &
	    + C(SC2H5) + C(SCH2CO) &
	    + C(SHCCO) + C(SCH2OH) &
	    + C(SC2H) + C(SC3H3) &
	    + C(SNXC4H3) + C(SIXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4)
      M(MM8) = M(MM8) + C(SA1XC6H6) + C(SA1XXC6H5) &
	    + C(SC5H5) + C(SC3H6) &
	    + C(SC4H2) + C(SOXC6H4) &
	    + C(SA2XC10H8) + C(SC9H8) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6) + C(SC9H7)


      K(R5) = 1.7800000000D+12 * exp(-1 * LT)
      K(R2F) = 4.5900000000D+01 &
	   * exp(2.7 * LT - 26190000 / RT)
      K(R2B) = 2.7517349858D+01 &
	   * exp(2.66385 * LT - 20100582.44 / RT)
      K(R10) = 9.4300000000D+12 * exp(-1 * LT)
      K(R9F) = 4.4000000000D+16 * exp(-2 * LT)
      K(R9B) = 1.1314226588D+20 &
	   * exp(-1.75351 * LT - 496136012.5 / RT)
      K(R3F) = 1.7300000000D+05 &
	   * exp(1.51 * LT - 14350000 / RT)
      K(R3B) = 1.9032628245D+06 &
	   * exp(1.40625 * LT - 77394524.74 / RT)
      K(R4F) = 3.9700000000D+01 &
	   * exp(2.4 * LT + 8830000 / RT)
      K(R4B) = 7.2853303338D+02 &
	   * exp(2.33241 * LT - 60303942.3 / RT)
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
      K0TROE = 6.3300000000D+13 * exp(-1.4 * LT)
      KINFTROE = 5.1200000000D+09 * exp(0.44 * LT)
      FCTROE = 0.5 * EXP( -TEMP / 1e-10 ) &
	   + 0.5 * EXP( -TEMP / 1e+10 )
      K(R12) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K(R1F) = 2.6400000000D+13 &
	   * exp(-0.67 * LT - 71300000 / RT)
      K(R1B) = 5.2766713024D+10 &
	   * exp(-0.234291 * LT - 479604.2721 / RT)
      K(R19F) = 1.0000000000D+13 * exp(-72510000 / RT)
      K(R19B) = 1.7850732872D+13 &
	   * exp(0.258327 * LT - 363978197.2 / RT)
      K(R18F) = 2.3800000000D+10 * exp(2090000 / RT)
      K(R18B) = 4.2484744234D+10 &
	   * exp(0.258327 * LT - 289378197.2 / RT)
      K(R16) = 7.4900000000D+10 * exp(-2660000 / RT)
      K(R17) = 4.0000000000D+10
      K(R15) = 3.9700000000D+09 * exp(-2810000 / RT)
      K(R13) = 3.6490000000D+03 &
	   * exp(2.068 * LT + 4574000 / RT)
      K(R26) = 2.6700000000D+38 &
	   * exp(-7 * LT - 157320000 / RT)
      K(RG04) = 5.7000000000D+10
      K(RG05) = 3.0000000000D+10
      K(RG09) = 6.7100000000D+10
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
      K(RG07) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG08) = 5.7100000000D+09 * exp(3160000 / RT)
      K(RG38F) = 3.0000000000D+10
      K(RG38B) = 2.1872817488D+10 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG39) = 6.8200000000D+07 &
	   * exp(0.25 * LT + 3910000 / RT)
      K(RG35) = 2.8000000000D+10
      K(RG34F) = 7.0000000000D+10
      K(RG34B) = 2.1409471364D+13 &
	   * exp(-0.574789 * LT - 68460673.16 / RT)
      K(RG31) = 1.5000000000D+10
      K(RG36) = 1.2000000000D+10
      K(RG30F) = 3.0000000000D+10
      K(RG30B) = 1.3435867056D+09 &
	   * exp(0.230182 * LT - 48621961.01 / RT)
      K(RG33) = 3.0000000000D+10
      K(RG19) = 5.8000000000D+09 * exp(-6280000 / RT)
      K(RG17F) = 1.1300000000D+04 &
	   * exp(2 * LT - 12550000 / RT)
      K(RG17B) = 7.6364598733D+03 &
	   * exp(2.18872 * LT - 86913497.78 / RT)
      K(RG20) = 2.4000000000D+09 * exp(-6280000 / RT)
      K(RG16) = 2.0000000000D+10
      K(RG25) = 4.0000000000D+10
      K0TROE = 1.0400000000D+20 &
	   * exp(-2.76 * LT - 6690000 / RT)
      KINFTROE = 6.0000000000D+11
      FCTROE = 0.438 * EXP( -TEMP / 91 ) &
	   + 0.562 * EXP( -TEMP / 5836 ) &
	   + 1 * EXP( -8552 / TEMP )
      K(RG14) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG21) = 5.0000000000D+09 * exp(-6280000 / RT)
      K(RG18F) = 5.0000000000D+02 &
	   * exp(2 * LT - 30250000 / RT)
      K(RG18B) = 2.0974636059D+05 &
	   * exp(1.4875 * LT - 61407685.18 / RT)
      K(RG26) = 1.6000000000D+12 * exp(-49970000 / RT)
      K(RG15) = 8.0000000000D+10
      K(RG69) = 3.0000000000D+10
      K(RG52) = 5.0600000000D+10
      K(RG55F) = 5.6000000000D+04 &
	   * exp(1.6 * LT - 22680000 / RT)
      K(RG55B) = 1.4686430659D+03 &
	   * exp(2.00876 * LT - 54566839.56 / RT)
      K(RG72) = 1.0000000000D+11
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
      K(RG74F) = 6.8400000000D+09 &
	   * exp(0.1 * LT - 44350000 / RT)
      K(RG74B) = 1.1614797579D+15 &
	   * exp(-1.0902 * LT - 7197984.006 / RT)
      K(RG53) = 3.3700000000D+10
      K(RG59) = 5.8700000000D+08 * exp(-57910000 / RT)
      K(RG57F) = 6.4400000000D+14 &
	   * exp(-1.34 * LT - 5930000 / RT)
      K(RG57B) = 2.3164910420D+13 &
	   * exp(-0.868957 * LT - 513851.5818 / RT)
      K(RG91F) = 1.0200000000D+06 &
	   * exp(1.5 * LT - 35980000 / RT)
      K(RG91B) = 3.2893910502D+02 &
	   * exp(1.98406 * LT - 20978674.53 / RT)
      K(RG93) = 6.0000000000D+10
      K(RG92F) = 1.0000000000D+05 &
	   * exp(1.6 * LT - 13050000 / RT)
      K(RG92B) = 5.9179879483D+02 &
	   * exp(2.01646 * LT - 67182616.83 / RT)
      K(RG90F) = 6.6000000000D+05 &
	   * exp(1.62 * LT - 45360000 / RT)
      K(RG90B) = 3.5503024345D+02 &
	   * exp(2.14021 * LT - 36448092.09 / RT)
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
      K(R28F) = 8.0000000000D+08 &
	   * exp(0.14 * LT - 30760000 / RT)
      K(R28B) = 4.3798719198D+15 &
	   * exp(-1.13743 * LT - 138461819.8 / RT)
      K(R31) = 3.0100000000D+10 * exp(-96230000 / RT)
      K0TROE = 1.1700000000D+18 &
	   * exp(-2.79 * LT - 17540000 / RT)
      KINFTROE = 1.3600000000D+07 * exp(-9980000 / RT)
      FCTROE = 1 * EXP( -0 / TEMP )
      K(R27) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K(R29F) = 8.7800000000D+07 &
	   * exp(0.03 * LT + 70000 / RT)
      K(R29B) = 4.8069094320D+14 &
	   * exp(-1.24743 * LT - 107631819.8 / RT)
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
      K(RG40F) = 9.0000000000D+09
      K(RG40B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R32) = 1.2000000000D+11
      K(R37F) = 2.2400000000D+15 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R37B) = 4.2292301785D+11 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(R38) = 1.2000000000D+07 &
	   * exp(0.81 * LT + 3040000 / RT)
      K(R33) = 3.0000000000D+10
      K(R35) = 3.0200000000D+10
      K(R34) = 3.0000000000D+10
      K(R36F) = 1.8700000000D+14 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R36B) = 3.5306519794D+10 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(RG70) = 2.6500000000D+10
      K(RG46) = 3.9000000000D+10 * exp(-14810000 / RT)
      K0TROE = 3.4670000000D+23 &
	   * exp(-2.669 * LT - 371385000 / RT)
      KINFTROE = 1.5300000000D+14 &
	   * exp(0.381 * LT - 368515000 / RT)
      FCTROE = 0.2176 * EXP( -TEMP / 271 ) &
	   + 0.7824 * EXP( -TEMP / 2755 ) &
	   + 1 * EXP( -6570 / TEMP )
      K(RG13) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG71) = 3.3200000000D+00 &
	   * exp(2.81 * LT - 24520000 / RT)
      K0TROE = 1.6120000000D+28 &
	   * exp(-4.102 * LT - 355207000 / RT)
      KINFTROE = 1.3680000000D+11 &
	   * exp(0.818 * LT - 335337000 / RT)
      FCTROE = 0.068 * EXP( -TEMP / 197 ) &
	   + 0.932 * EXP( -TEMP / 1540 ) &
	   + 1 * EXP( -10300 / TEMP )
      K(RG12) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG47) = 3.4300000000D+06 &
	   * exp(1.18 * LT + 1870000 / RT)
      K(RG50) = 9.4600000000D+10 * exp(2160000 / RT)
      K(RG45F) = 5.7400000000D+04 &
	   * exp(1.9 * LT - 11470000 / RT)
      K(RG45B) = 9.5590953891D+01 &
	   * exp(2.34945 * LT - 74956081.25 / RT)
      K(RG89) = 1.8000000000D+10 * exp(-3770000 / RT)
      K(RG86F) = 3.2800000000D+10 &
	   * exp(-0.09 * LT - 2550000 / RT)
      K(RG86B) = 1.1826581729D+05 &
	   * exp(1.35911 * LT - 12643307.5 / RT)
      K(RG85F) = 1.6500000000D+08 &
	   * exp(0.65 * LT + 1190000 / RT)
      K(RG85B) = 1.6539583927D+04 &
	   * exp(1.62807 * LT - 14319455.91 / RT)
      K0TROE = 1.2680000000D+30 &
	   * exp(-5.201 * LT - 150084000 / RT)
      KINFTROE = 5.3930000000D+12 &
	   * exp(0.069 * LT - 137824000 / RT)
      FCTROE = 0.2813 * EXP( -TEMP / 103 ) &
	   + 0.7187 * EXP( -TEMP / 1291 ) &
	   + 1 * EXP( -4160 / TEMP )
      K(RG43) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG84) = 2.0000000000D+10
      K(RG41F) = 7.0000000000D+09
      K(RG41B) = 5.1036574138D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG42) = 1.4000000000D+10
      K(RG11) = 1.9000000000D+11 * exp(-66070000 / RT)
      K(RG108F) = 3.3100000000D+03 &
	   * exp(2.26 * LT - 3770000 / RT)
      K(RG108B) = 9.6742642418D+06 &
	   * exp(1.66346 * LT - 127161378.7 / RT)
      K(RG107) = 1.0000000000D+10 * exp(3160000 / RT)
      K(RG105) = 5.0000000000D+10
      K(RG106) = 2.0000000000D+10
      K(RR017F) = 2.4100000000D+10
      K(RR017B) = 3.0422905514D+16 &
	   * exp(-1.37788 * LT - 159986443 / RT)
      K0TROE = 2.6000000000D+27 &
	   * exp(-4.8 * LT - 7950000 / RT)
      KINFTROE = 1.0000000000D+14 * exp(-1 * LT)
      FCTROE = 0.3536 * EXP( -TEMP / 132 ) &
	   + 0.6464 * EXP( -TEMP / 1315 ) &
	   + 1 * EXP( -5566 / TEMP )
      K(RG104) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG118F) = 1.8100000000D+10
      K(RG118B) = 8.8241850413D+13 &
	   * exp(-0.560387 * LT - 129480796.2 / RT)
      K(RG121) = 7.5300000000D+03 &
	   * exp(1.55 * LT - 8810000 / RT)
      K0TROE = 6.3400000000D+25 &
	   * exp(-4.66 * LT - 15820000 / RT)
      KINFTROE = 1.7100000000D+07 &
	   * exp(1.27 * LT - 11330000 / RT)
      FCTROE = 0.2122 * EXP( -TEMP / -10210 )
      K(RG115F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 1.4048536324D+28 &
	   * exp(-4.48812 * LT - 162621237.8 / RT)
      KINFTROE = 3.7891162641D+09 &
	   * exp(1.44188 * LT - 158131237.8 / RT)
      K(RG115B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG117) = 1.2500000000D+04 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR008F) = 7.8000000000D+10
      K(RR008B) = 5.6362388680D+19 &
	   * exp(-1.29954 * LT - 250786354.2 / RT)
      K(RG119F) = 2.6300000000D+03 &
	   * exp(2.14 * LT - 71380000 / RT)
      K(RG119B) = 9.8996184977D+00 &
	   * exp(2.63279 * LT - 11033146.07 / RT)
      K(RG116) = 8.1000000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR004F) = 1.9000000000D+11
      K(RR004B) = 2.5098929533D+16 &
	   * exp(-1.35613 * LT - 105055737.5 / RT)
      K(RG122) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RG130F) = 3.0000000000D+10
      K(RG130B) = 3.1644553172D+10 &
	   * exp(0.178352 * LT - 286290250 / RT)
      K(RB07) = 9.6000000000D+08
      K(RR011) = 9.0000000000D+10
      K(RR013) = 9.0300000000D+09 * exp(3200000 / RT)
      K(RG132) = 5.0000000000D+09
      K(RG133) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RG135) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K0TROE = 1.4000000000D+24 &
	   * exp(-3.86 * LT - 13890000 / RT)
      KINFTROE = 6.0800000000D+09 &
	   * exp(0.27 * LT - 1170000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(RG129F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.7814512533D+29 &
	   * exp(-3.99563 * LT - 480409256.8 / RT)
      KINFTROE = 3.3793731157D+15 &
	   * exp(0.134367 * LT - 467689256.8 / RT)
      K(RG129B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG159) = 3.8900000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RIC67) = 4.1000000000D+08 * exp(-30140000 / RT)
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
      K(RG162F) = 2.2700000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RG162B) = 1.7745629538D+02 &
	   * exp(1.96566 * LT - 13974138.86 / RT)
      K(RG160F) = 1.3100000000D-04 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RG160B) = 6.0605344705D-07 &
	   * exp(4.58212 * LT - 26016755.69 / RT)
      K(RG158) = 7.1500000000D+01 &
	   * exp(2.47 * LT - 3890000 / RT)
      K(RG156F) = 1.2700000000D+02 &
	   * exp(2.75 * LT - 48740000 / RT)
      K(RG156B) = 5.3406078444D-02 &
	   * exp(3.23587 * LT - 15312230.95 / RT)
      K(RG171) = 1.9200000000D+04 &
	   * exp(1.02 * LT + 8510000 / RT)
      K(RR034) = 1.2000000000D+11
      K0TROE = 1.9900000000D+35 &
	   * exp(-7.08 * LT - 27970000 / RT)
      KINFTROE = 5.2100000000D+14 &
	   * exp(-0.99 * LT - 6610000 / RT)
      FCTROE = 0.1578 * EXP( -TEMP / 125 ) &
	   + 0.8422 * EXP( -TEMP / 2219 ) &
	   + 1 * EXP( -6882 / TEMP )
      K(RG164) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG178) = 8.4300000000D+11 * exp(-93120000 / RT)
      K(RG174F) = 1.7000000000D+02 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RG174B) = 3.1981668158D-02 &
	   * exp(3.25918 * LT - 35436570.06 / RT)
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
      K(RG110) = 1.0000000000D+11
      K(RG109F) = 1.0000000000D+11
      K(RG109B) = 1.1573048940D+05 &
	   * exp(1.56102 * LT - 71793685.47 / RT)
      K(RG114) = 1.0000000000D+10
      K(RR009) = 1.0000000000D+08 * exp(-12550000 / RT)
      K(RG112) = 5.0000000000D+10
      K(RG113) = 3.0000000000D+10
      K(RR024) = 1.0000000000D+10
      K(RR022) = 5.0000000000D+10
      K(RG124F) = 1.5000000000D+06 &
	   * exp(1.38 * LT - 2570000 / RT)
      K(RG124B) = 4.0698109539D-01 &
	   * exp(2.9634 * LT - 131854077.6 / RT)
      K(RG123F) = 5.0000000000D+10 * exp(-33470000 / RT)
      K(RG123B) = 3.8326334707D+07 &
	   * exp(0.59717 * LT - 22499718.93 / RT)
      K(RG127F) = 7.5000000000D+09 * exp(-8370000 / RT)
      K(RG127B) = 6.3247186168D+07 &
	   * exp(0.493423 * LT - 60444243.67 / RT)
      K(RR025) = 5.0000000000D+09
      K(RR068) = 2.5000000000D+10
      K(RLPP009F) = 2.2600000000D+43 &
	   * exp(-9.862 * LT - 70549000 / RT)
      K(RLPP009B) = 3.2581191141D+59 &
	   * exp(-11.3197 * LT - 676812672.1 / RT)
      K(RR055F) = 7.9400000000D+26 &
	   * exp(-5.06 * LT - 20340000 / RT)
      K(RR055B) = 8.1610655573D+31 &
	   * exp(-5.06898 * LT - 393338860.7 / RT)
      K(RR058) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR062) = 1.7000000000D+02 &
	   * exp(1.7 * LT - 6280000 / RT)
      K(RP011F) = 5.7700000000D+34 &
	   * exp(-7 * LT - 131820000 / RT)
      K(RP011B) = 5.1317962822D+44 &
	   * exp(-8.24617 * LT - 262035013.4 / RT)
      K(RCP16F) = 6.8700000000D+52 &
	   * exp(-12.5 * LT - 176000000 / RT)
      K(RCP16B) = 5.5097990876D+64 &
	   * exp(-13.5715 * LT - 487565243.7 / RT)
      K(RR089) = 6.2500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR079) = 4.4900000000D+04 &
	   * exp(1.92 * LT - 23810000 / RT)
      K(RR080F) = 8.0500000000D+02 &
	   * exp(2.22 * LT - 3100000 / RT)
      K(RR080B) = 2.0139176673D+01 &
	   * exp(2.47547 * LT - 126237151.9 / RT)
      K(RR005F) = 3.4600000000D+09 &
	   * exp(0.44 * LT - 22860000 / RT)
      K(RR005B) = 1.8216946136D+04 &
	   * exp(1.58056 * LT - 46357562.82 / RT)
      K(RR088) = 4.0500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR081F) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RR081B) = 1.7839562769D+12 &
	   * exp(-0.16099 * LT - 162124535 / RT)
      K(RR078F) = 8.5000000000D+01 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RR078B) = 1.9329120587D-01 &
	   * exp(3.05922 * LT - 84112627.12 / RT)
      K(RR014F) = 4.0400000000D+42 &
	   * exp(-7.67 * LT - 467900000 / RT)
      K(RR014B) = 3.9209380170D+30 &
	   * exp(-6.00376 * LT - 34921721.89 / RT)
      K(RR143) = 3.5000000000D+04 &
	   * exp(1.6 * LT + 4070000 / RT)
      K(RR141F) = 8.0000000000D+18 &
	   * exp(-2.39 * LT - 46780000 / RT)
      K(RR141B) = 4.3155004306D+12 &
	   * exp(-0.859397 * LT - 80320978.73 / RT)
      K(RIC66) = 1.7300000000D+10 &
	   * exp(0.03 * LT - 7520000 / RT)
      K(RIC64) = 9.7700000000D-09 &
	   * exp(5.36 * LT - 71250000 / RT)
      K(RIC65) = 9.8800000000D+18 &
	   * exp(-1.59 * LT - 168820000 / RT)
      K(RH10F) = 1.3700000000D+36 &
	   * exp(-7.87 * LT - 64610000 / RT)
      K(RH10B) = 4.7217121748D+38 &
	   * exp(-7.9297 * LT - 199616108.7 / RT)
      K(RH07) = 9.5600000000D+09 * exp(-130120000 / RT)
      K(RH12) = 3.3000000000D+09 &
	   * exp(-0.25 * LT - 9940000 / RT)
      K(ROX02F) = 5.0000000000D+75 &
	   * exp(-19.31 * LT - 284180000 / RT)
      K(ROX02B) = 2.7802587455D+86 &
	   * exp(-20.3484 * LT - 537844017.3 / RT)
      K(RH17) = 6.3000000000D+22 &
	   * exp(-3.34 * LT - 41880000 / RT)
      K(RH19) = 1.5000000000D+10
      K(RLPCP18) = 4.9100000000D+28 &
	   * exp(-4.85 * LT - 103650000 / RT)
      K(RCP17) = 6.3900000000D+26 &
	   * exp(-4.03 * LT - 147300000 / RT)
      K(RLPP015) = 9.4500000000D-06 &
	   * exp(4.47 * LT - 18710000 / RT)
      K(RI00F) = 1.0000000000D+10
      K(RI00B) = 5.2485026204D+25 &
	   * exp(-2.03247 * LT - 521300387.2 / RT)
      K(RLPOX00) = 2.0330000000D+44 &
	   * exp(-9.396 * LT - 100630000 / RT)
      K(RK012F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RK012B) = 2.8151364037D+14 &
	   * exp(0.27021 * LT - 191386104.5 / RT)
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
      K(RLPOX03F) = 5.5340000000D+05 &
	   * exp(1.811 * LT - 68565000 / RT)
      K(RLPOX03B) = 7.9798323034D+01 &
	   * exp(2.37273 * LT - 25607829.1 / RT)
      K(RLPOX04) = 2.1510000000D+01 &
	   * exp(2.691 * LT - 3215000 / RT)
      K(RI03F) = 3.1600000000D+01 &
	   * exp(2.48 * LT - 46280000 / RT)
      K(RI03B) = 4.5268324245D+08 &
	   * exp(1.36023 * LT - 103588858.9 / RT)
      K(RT04F) = 5.8300000000D+64 &
	   * exp(-14.15 * LT - 285890000 / RT)
      K(RT04B) = 8.6112522502D+55 &
	   * exp(-12.1057 * LT - 218394895.1 / RT)
      K(RT20) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RT05F) = 8.2000000000D+14 * exp(-337550000 / RT)
      K(RT05B) = 4.3614434871D+01 &
	   * exp(2.2454 * LT - 46091406.02 / RT)
      K(RP106F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP106B) = 9.5260034588D+36 &
	   * exp(-6.03013 * LT - 365946869.1 / RT)
      K(RK100F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RK100B) = 3.3117446183D+14 &
	   * exp(0.979368 * LT - 386582470.9 / RT)
      K(RP105F) = 3.6000000000D+14 &
	   * exp(-1.44 * LT - 65930000 / RT)
      K(RP105B) = 7.8690280050D+25 &
	   * exp(-3.13503 * LT - 311868721.5 / RT)
      K(RLPST00) = 1.0000000000D+10
      K(RK017F) = 1.2900000000D+05 &
	   * exp(1.89 * LT - 73550000 / RT)
      K(RK017B) = 4.5060198908D+01 &
	   * exp(2.47097 * LT - 19964083.4 / RT)
      K(RP016F) = 4.3000000000D+60 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(RP016B) = 6.4261715627D+54 &
	   * exp(-12.2493 * LT - 132912595.6 / RT)
      K(RP018F) = 7.8000000000D+00 &
	   * exp(2.68 * LT - 3070000 / RT)
      K(RP018B) = 2.9974411044D-02 &
	   * exp(3.15722 * LT - 12528608.14 / RT)
      K(RK020F) = 7.1800000000D+10 &
	   * exp(1.02 * LT - 161810000 / RT)
      K(RK020B) = 4.1382935254D+06 &
	   * exp(1.32078 * LT - 30384498.09 / RT)
      K(RP104F) = 3.8000000000D+04 &
	   * exp(1.62 * LT - 18570000 / RT)
      K(RP104B) = 1.0608179381D+14 &
	   * exp(0.397642 * LT - 279884457.4 / RT)
      K(RST13) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RI26) = 1.1900000000D+15 * exp(-377480000 / RT)
      K(RI07) = 1.5160000000D+09 &
	   * exp(0.18 * LT - 20724000 / RT)
      K(RI18) = 1.9600000000D+28 &
	   * exp(-4.85 * LT - 103650000 / RT)
      K(RI25) = 2.3800000000D+15 * exp(-377480000 / RT)
      K(RI09F) = 1.8000000000D-04 * exp(4 * LT)
      K(RI09B) = 8.2296634648D-04 &
	   * exp(4.0222 * LT - 102619350.6 / RT)
      K(RI02F) = 2.8000000000D+10 * exp(-9450000 / RT)
      K(RI02B) = 6.8863488094D+07 &
	   * exp(0.542404 * LT - 103157442.7 / RT)
      K(RI05) = 4.7700000000D+01 &
	   * exp(2.71 * LT - 4630000 / RT)
      K(RI01F) = 1.7300000000D+68 &
	   * exp(-15.16 * LT - 486900000 / RT)
      K(RI01B) = 1.8203611359D+63 &
	   * exp(-14.9678 * LT - 147515954.9 / RT)
      K(RI06F) = 3.0800000000D+03 * exp(2 * LT)
      K(RI06B) = 8.3336328636D+01 &
	   * exp(2.43866 * LT - 156751967.4 / RT)
      K(RK110F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK110B) = 1.8883129559D-03 &
	   * exp(3.49973 * LT - 16231119.74 / RT)
      K(RP108F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP108B) = 6.5580381913D+54 &
	   * exp(-12.2468 * LT - 121345107.2 / RT)
      K(RK109F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK109B) = 4.7232496168D+01 &
	   * exp(2.45348 * LT - 6416595.001 / RT)
      K(RG28F) = 1.5000000000D+10 * exp(-2510000 / RT)
      K(RG28B) = 1.0936408744D+10 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)

      W(R5) = K(R5) * C(SH) * C(SH) * M(MM1)
      W(R2F) = K(R2F) * C(SO) * C(SH2)
      W(R2B) = K(R2B) * C(SOH) * C(SH)
      W(R10) = K(R10) * C(SH) * C(SO) * M(MM3)
      W(R9F) = K(R9F) * C(SH) * C(SOH) * M(MM2)
      W(R9B) = K(R9B) * C(SH2O) * M(MM2)
      W(R3F) = K(R3F) * C(SOH) * C(SH2)
      W(R3B) = K(R3B) * C(SH2O) * C(SH)
      W(R4F) = K(R4F) * C(SOH) * C(SOH)
      W(R4B) = K(R4B) * C(SH2O) * C(SO)
      W(R14F) = K(R14F) * C(SOH) * C(SOH)
      W(R14B) = K(R14B) * C(SH2O2)
      W(R12) = K(R12) * C(SH) * C(SO2)
      W(R1F) = K(R1F) * C(SH) * C(SO2)
      W(R1B) = K(R1B) * C(SOH) * C(SO)
      W(R19F) = K(R19F) * C(SHO2) * C(SOH)
      W(R19B) = K(R19B) * C(SO2) * C(SH2O)
      W(R18F) = K(R18F) * C(SHO2) * C(SOH)
      W(R18B) = K(R18B) * C(SO2) * C(SH2O)
      W(R16) = K(R16) * C(SHO2) * C(SH)
      W(R17) = K(R17) * C(SHO2) * C(SO)
      W(R15) = K(R15) * C(SHO2) * C(SH)
      W(R13) = K(R13) * C(SH) * C(SHO2)
      W(R26) = K(R26) * C(SH2O2) * C(SOH)
      W(RG04) = K(RG04) * C(SCH) * C(SO)
      W(RG05) = K(RG05) * C(SCH) * C(SOH)
      W(RG09) = K(RG09) * C(SCH) * C(SO2)
      W(RG06F) = K(RG06F) * C(SCH) * C(SH2)
      W(RG06B) = K(RG06B) * C(SH) * C(STXCH2)
      W(RG07) = K(RG07) * C(SCH) * C(SH2)
      W(RG08) = K(RG08) * C(SCH) * C(SH2O)
      W(RG38F) = K(RG38F) * C(SSXCH2) * C(SH2O)
      W(RG38B) = K(RG38B) * C(SH2O) * C(STXCH2)
      W(RG39) = K(RG39) * C(SSXCH2) * C(SH2O)
      W(RG35) = K(RG35) * C(SSXCH2) * C(SO2)
      W(RG34F) = K(RG34F) * C(SSXCH2) * C(SH2)
      W(RG34B) = K(RG34B) * C(SH) * C(SCH3)
      W(RG31) = K(RG31) * C(SSXCH2) * C(SO)
      W(RG36) = K(RG36) * C(SSXCH2) * C(SO2)
      W(RG30F) = K(RG30F) * C(SSXCH2) * C(SH)
      W(RG30B) = K(RG30B) * C(SH2) * C(SCH)
      W(RG33) = K(RG33) * C(SSXCH2) * C(SOH)
      W(RG19) = K(RG19) * C(STXCH2) * C(SO2)
      W(RG17F) = K(RG17F) * C(STXCH2) * C(SOH)
      W(RG17B) = K(RG17B) * C(SH2O) * C(SCH)
      W(RG20) = K(RG20) * C(STXCH2) * C(SO2)
      W(RG16) = K(RG16) * C(STXCH2) * C(SOH)
      W(RG25) = K(RG25) * C(STXCH2) * C(SCH)
      W(RG14) = K(RG14) * C(STXCH2) * C(SH)
      W(RG21) = K(RG21) * C(STXCH2) * C(SO2)
      W(RG18F) = K(RG18F) * C(STXCH2) * C(SH2)
      W(RG18B) = K(RG18B) * C(SCH3) * C(SH)
      W(RG26) = K(RG26) * C(STXCH2) * C(STXCH2)
      W(RG15) = K(RG15) * C(STXCH2) * C(SO)
      W(RG69) = K(RG69) * C(SCH3) * C(SCH)
      W(RG52) = K(RG52) * C(SCH3) * C(SO)
      W(RG55F) = K(RG55F) * C(SCH3) * C(SOH)
      W(RG55B) = K(RG55B) * C(SH2O) * C(STXCH2)
      W(RG72) = K(RG72) * C(SCH3) * C(STXCH2)
      W(RG51F) = K(RG51F) * C(SCH3) * C(SH)
      W(RG51B) = K(RG51B) * C(SCH4)
      W(RG74F) = K(RG74F) * C(SCH3) * C(SCH3)
      W(RG74B) = K(RG74B) * C(SH) * C(SC2H5)
      W(RG53) = K(RG53) * C(SCH3) * C(SO)
      W(RG59) = K(RG59) * C(SCH3) * C(SO2)
      W(RG57F) = K(RG57F) * C(SCH3) * C(SOH)
      W(RG57B) = K(RG57B) * C(SH2O) * C(SSXCH2)
      W(RG91F) = K(RG91F) * C(SCH4) * C(SO)
      W(RG91B) = K(RG91B) * C(SOH) * C(SCH3)
      W(RG93) = K(RG93) * C(SCH4) * C(SCH)
      W(RG92F) = K(RG92F) * C(SCH4) * C(SOH)
      W(RG92B) = K(RG92B) * C(SH2O) * C(SCH3)
      W(RG90F) = K(RG90F) * C(SCH4) * C(SH)
      W(RG90B) = K(RG90B) * C(SH2) * C(SCH3)
      W(RG24F) = K(RG24F) * C(STXCH2) * C(SCO)
      W(RG24B) = K(RG24B) * C(SCH2CO)
      W(R28F) = K(R28F) * C(SCO) * C(SOH)
      W(R28B) = K(R28B) * C(SH) * C(SCO2)
      W(R31) = K(R31) * C(SCO) * C(SHO2)
      W(R27) = K(R27) * C(SCO) * C(SO)
      W(R29F) = K(R29F) * C(SCO) * C(SOH)
      W(R29B) = K(R29B) * C(SH) * C(SCO2)
      W(RG10F) = K(RG10F) * C(SCH) * C(SCO)
      W(RG10B) = K(RG10B) * C(SHCCO)
      W(RG40F) = K(RG40F) * C(SSXCH2) * C(SCO)
      W(RG40B) = K(RG40B) * C(SCO) * C(STXCH2)
      W(R32) = K(R32) * C(SHCO) * C(SH)
      W(R37F) = K(R37F) * C(SHCO) * C(SH2O)
      W(R37B) = K(R37B) * C(SH2O) * C(SH) * C(SCO)
      W(R38) = K(R38) * C(SHCO) * C(SO2)
      W(R33) = K(R33) * C(SHCO) * C(SO)
      W(R35) = K(R35) * C(SHCO) * C(SOH)
      W(R34) = K(R34) * C(SHCO) * C(SO)
      W(R36F) = K(R36F) * C(SHCO) * M(MM8)
      W(R36B) = K(R36B) * C(SH) * C(SCO) * M(MM8)
      W(RG70) = K(RG70) * C(SCH3) * C(SHCO)
      W(RG46) = K(RG46) * C(SCH2O) * C(SO)
      W(RG13) = K(RG13) * C(SCH2O)
      W(RG71) = K(RG71) * C(SCH3) * C(SCH2O)
      W(RG12) = K(RG12) * C(SCH2O)
      W(RG47) = K(RG47) * C(SCH2O) * C(SOH)
      W(RG50) = K(RG50) * C(SCH2O) * C(SCH)
      W(RG45F) = K(RG45F) * C(SCH2O) * C(SH)
      W(RG45B) = K(RG45B) * C(SH2) * C(SHCO)
      W(RG89) = K(RG89) * C(SCH2OH) * C(SO2)
      W(RG86F) = K(RG86F) * C(SCH2OH) * C(SH)
      W(RG86B) = K(RG86B) * C(SH2O) * C(SSXCH2)
      W(RG85F) = K(RG85F) * C(SCH2OH) * C(SH)
      W(RG85B) = K(RG85B) * C(SOH) * C(SCH3)
      W(RG43) = K(RG43) * C(SCH2OH)
      W(RG84) = K(RG84) * C(SCH2OH) * C(SH)
      W(RG41F) = K(RG41F) * C(SSXCH2) * C(SCO2)
      W(RG41B) = K(RG41B) * C(SCO2) * C(STXCH2)
      W(RG42) = K(RG42) * C(SSXCH2) * C(SCO2)
      W(RG11) = K(RG11) * C(SCH) * C(SCO2)
      W(RG108F) = K(RG108F) * C(SC2H) * C(SH2)
      W(RG108B) = K(RG108B) * C(SH) * C(SC2H2)
      W(RG107) = K(RG107) * C(SC2H) * C(SO2)
      W(RG105) = K(RG105) * C(SC2H) * C(SO)
      W(RG106) = K(RG106) * C(SC2H) * C(SOH)
      W(RR017F) = K(RR017F) * C(SC2H) * C(SCH3)
      W(RR017B) = K(RR017B) * C(SH) * C(SC3H3)
      W(RG104) = K(RG104) * C(SC2H) * C(SH)
      W(RG118F) = K(RG118F) * C(SC2H) * C(SOH)
      W(RG118B) = K(RG118B) * C(SO) * C(SC2H2)
      W(RG121) = K(RG121) * C(SC2H2) * C(SOH)
      W(RG115F) = K(RG115F) * C(SC2H2) * C(SH)
      W(RG115B) = K(RG115B) * C(SC2H3)
      W(RG117) = K(RG117) * C(SC2H2) * C(SO)
      W(RR008F) = K(RR008F) * C(SC2H2) * C(SC2H)
      W(RR008B) = K(RR008B) * C(SNXC4H3)
      W(RG119F) = K(RG119F) * C(SC2H2) * C(SOH)
      W(RG119B) = K(RG119B) * C(SH2O) * C(SC2H)
      W(RG116) = K(RG116) * C(SC2H2) * C(SO)
      W(RR004F) = K(RR004F) * C(SC2H2) * C(SSXCH2)
      W(RR004B) = K(RR004B) * C(SH) * C(SC3H3)
      W(RG122) = K(RG122) * C(SC2H2) * C(SOH)
      W(RG130F) = K(RG130F) * C(SC2H3) * C(SH)
      W(RG130B) = K(RG130B) * C(SH2) * C(SC2H2)
      W(RB07) = K(RB07) * C(SC2H3) * C(SC2H3)
      W(RR011) = K(RR011) * C(SC2H3) * C(SHCO)
      W(RR013) = K(RR013) * C(SC2H3) * C(SCH3)
      W(RG132) = K(RG132) * C(SC2H3) * C(SOH)
      W(RG133) = K(RG133) * C(SC2H3) * C(SO2)
      W(RG135) = K(RG135) * C(SC2H3) * C(SO2)
      W(RG129F) = K(RG129F) * C(SC2H3) * C(SH)
      W(RG129B) = K(RG129B) * C(SC2H4)
      W(RG159) = K(RG159) * C(SC2H4) * C(SO)
      W(RIC67) = K(RIC67) * C(SC2H4) * C(SCH3)
      W(RG155F) = K(RG155F) * C(SC2H4) * C(SH)
      W(RG155B) = K(RG155B) * C(SC2H5)
      W(RG162F) = K(RG162F) * C(SC2H4) * C(SCH3)
      W(RG162B) = K(RG162B) * C(SCH4) * C(SC2H3)
      W(RG160F) = K(RG160F) * C(SC2H4) * C(SOH)
      W(RG160B) = K(RG160B) * C(SH2O) * C(SC2H3)
      W(RG158) = K(RG158) * C(SC2H4) * C(SO)
      W(RG156F) = K(RG156F) * C(SC2H4) * C(SH)
      W(RG156B) = K(RG156B) * C(SH2) * C(SC2H3)
      W(RG171) = K(RG171) * C(SC2H5) * C(SO2)
      W(RR034) = K(RR034) * C(SC2H5) * C(SHCO)
      W(RG164) = K(RG164) * C(SC2H5) * C(SH)
      W(RG178) = K(RG178) * C(SC2H6) * C(SCH3)
      W(RG174F) = K(RG174F) * C(SC2H6) * C(SH)
      W(RG174B) = K(RG174B) * C(SH2) * C(SC2H5)
      W(RG173F) = K(RG173F) * C(SC2H6)
      W(RG173B) = K(RG173B) * C(SCH3) * C(SCH3)
      W(RG110) = K(RG110) * C(SHCCO) * C(SO)
      W(RG109F) = K(RG109F) * C(SHCCO) * C(SH)
      W(RG109B) = K(RG109B) * C(SCO) * C(SSXCH2)
      W(RG114) = K(RG114) * C(SHCCO) * C(SHCCO)
      W(RR009) = K(RR009) * C(SC2H2) * C(SHCCO)
      W(RG112) = K(RG112) * C(SHCCO) * C(SCH)
      W(RG113) = K(RG113) * C(SHCCO) * C(STXCH2)
      W(RR024) = K(RR024) * C(SHCCO) * C(SOH)
      W(RR022) = K(RR022) * C(SHCCO) * C(SCH3)
      W(RG124F) = K(RG124F) * C(SCH2CO) * C(SH)
      W(RG124B) = K(RG124B) * C(SCO) * C(SCH3)
      W(RG123F) = K(RG123F) * C(SCH2CO) * C(SH)
      W(RG123B) = K(RG123B) * C(SH2) * C(SHCCO)
      W(RG127F) = K(RG127F) * C(SCH2CO) * C(SOH)
      W(RG127B) = K(RG127B) * C(SH2O) * C(SHCCO)
      W(RR025) = K(RR025) * C(SCH2CO) * C(SOH)
      W(RR068) = K(RR068) * C(SC3H3) * C(SHCO)
      W(RLPP009F) = K(RLPP009F) * C(SC3H3) * C(SC3H3)
      W(RLPP009B) = K(RLPP009B) * C(SA1XC6H6)
      W(RR055F) = K(RR055F) * C(SC3H3) * C(SH)
      W(RR055B) = K(RR055B) * C(SPXC3H4)
      W(RR058) = K(RR058) * C(SC3H3) * C(SOH)
      W(RR062) = K(RR062) * C(SC3H3) * C(SO2)
      W(RP011F) = K(RP011F) * C(SC3H3) * C(SC3H3)
      W(RP011B) = K(RP011B) * C(SH) * C(SA1XXC6H5)
      W(RCP16F) = K(RCP16F) * C(SC3H3) * C(SC2H2)
      W(RCP16B) = K(RCP16B) * C(SC5H5)
      W(RR089) = K(RR089) * C(SPXC3H4) * C(SO)
      W(RR079) = K(RR079) * C(SPXC3H4) * C(SO)
      W(RR080F) = K(RR080F) * C(SPXC3H4) * C(SOH)
      W(RR080B) = K(RR080B) * C(SH2O) * C(SC3H3)
      W(RR005F) = K(RR005F) * C(SPXC3H4) * C(SH)
      W(RR005B) = K(RR005B) * C(SCH3) * C(SC2H2)
      W(RR088) = K(RR088) * C(SPXC3H4) * C(SO)
      W(RR081F) = K(RR081F) * C(SPXC3H4) * C(SCH3)
      W(RR081B) = K(RR081B) * C(SCH4) * C(SC3H3)
      W(RR078F) = K(RR078F) * C(SPXC3H4) * C(SH)
      W(RR078B) = K(RR078B) * C(SH2) * C(SC3H3)
      W(RR014F) = K(RR014F) * C(SC3H6)
      W(RR014B) = K(RR014B) * C(SCH3) * C(SC2H3)
      W(RR143) = K(RR143) * C(SC3H6) * C(SO)
      W(RR141F) = K(RR141F) * C(SC3H6) * C(SH)
      W(RR141B) = K(RR141B) * C(SCH3) * C(SC2H4)
      W(RIC66) = K(RIC66) * C(SC3H6) * C(SH)
      W(RIC64) = K(RIC64) * C(SIXC3H7)
      W(RIC65) = K(RIC65) * C(SIXC3H7)
      W(RH10F) = K(RH10F) * C(SC4H2) * C(SH)
      W(RH10B) = K(RH10B) * C(SNXC4H3)
      W(RH07) = K(RH07) * C(SC4H2) * C(SO2)
      W(RH12) = K(RH12) * C(SC4H2) * C(SOH)
      W(ROX02F) = K(ROX02F) * C(SC4H2) * C(SC2H2)
      W(ROX02B) = K(ROX02B) * C(SOXC6H4)
      W(RH17) = K(RH17) * C(SNXC4H3) * C(SH)
      W(RH19) = K(RH19) * C(SNXC4H3) * C(SH)
      W(RLPCP18) = K(RLPCP18) * C(SC5H5) * C(SCH3)
      W(RCP17) = K(RCP17) * C(SC5H5) * C(SC5H5)
      W(RLPP015) = K(RLPP015) * C(SA1XXC6H5) * C(SC2H4)
      W(RI00F) = K(RI00F) * C(SA1XXC6H5) * C(SC3H3)
      W(RI00B) = K(RI00B) * C(SC9H8)
      W(RLPOX00) = K(RLPOX00) * C(SH) * C(SA1XXC6H5)
      W(RK012F) = K(RK012F) * C(SA1XXC6H5) * C(SC2H2)
      W(RK012B) = K(RK012B) * C(SA1C2H2XC8H7)
      W(ROX01F) = K(ROX01F) * C(SA1XXC6H5)
      W(ROX01B) = K(ROX01B) * C(SH) * C(SOXC6H4)
      W(RLPOX03F) = K(RLPOX03F) * C(SA1XC6H6) * C(SH)
      W(RLPOX03B) = K(RLPOX03B) * C(SH2) * C(SA1XXC6H5)
      W(RLPOX04) = K(RLPOX04) * C(SA1XC6H6) * C(SOH)
      W(RI03F) = K(RI03F) * C(SA1CH2XC7H7) * C(SC2H2)
      W(RI03B) = K(RI03B) * C(SH) * C(SC9H8)
      W(RT04F) = K(RT04F) * C(SA1CH2XC7H7) * C(SH)
      W(RT04B) = K(RT04B) * C(SCH3) * C(SA1XXC6H5)
      W(RT20) = K(RT20) * C(SA1CH2XC7H7) * C(SC3H3)
      W(RT05F) = K(RT05F) * C(SA1CH2XC7H7)
      W(RT05B) = K(RT05B) * C(SC2H2) * C(SC5H5)
      W(RP106F) = K(RP106F) * C(SA1C2HYXC8H5) * C(SC2H4)
      W(RP106B) = K(RP106B) * C(SH) * C(SA2XC10H8)
      W(RK100F) = K(RK100F) * C(SA1C2HYXC8H5) * C(SC2H2)
      W(RK100B) = K(RK100B) * C(SA2XXC10H7)
      W(RP105F) = K(RP105F) * C(SA1C2HXC8H6) * C(SC2H3)
      W(RP105B) = K(RP105B) * C(SH) * C(SA2XC10H8)
      W(RLPST00) = K(RLPST00) * C(SA1C2HXC8H6) * C(SOH)
      W(RK017F) = K(RK017F) * C(SA1C2HXC8H6) * C(SH)
      W(RK017B) = K(RK017B) * C(SH2) * C(SA1C2HYXC8H5)
      W(RP016F) = K(RP016F) * C(SA1C2HXC8H6)
      W(RP016B) = K(RP016B) * C(SH) * C(SA1C2HYXC8H5)
      W(RP018F) = K(RP018F) * C(SA1C2HXC8H6) * C(SOH)
      W(RP018B) = K(RP018B) * C(SH2O) * C(SA1C2HYXC8H5)
      W(RK020F) = K(RK020F) * C(SA1C2H2XC8H7)
      W(RK020B) = K(RK020B) * C(SH) * C(SA1C2HXC8H6)
      W(RP104F) = K(RP104F) * C(SA1C2H2XC8H7) * C(SC2H2)
      W(RP104B) = K(RP104B) * C(SH) * C(SA2XC10H8)
      W(RST13) = K(RST13) * C(SA1C2H2XC8H7) * C(SO2)
      W(RI26) = K(RI26) * C(SC9H7)
      W(RI07) = K(RI07) * C(SHO2) * C(SC9H7)
      W(RI18) = K(RI18) * C(SC9H7) * C(SCH3)
      W(RI25) = K(RI25) * C(SC9H7)
      W(RI09F) = K(RI09F) * C(SC9H8) * C(SCH3)
      W(RI09B) = K(RI09B) * C(SCH4) * C(SC9H7)
      W(RI02F) = K(RI02F) * C(SC9H8) * C(SH)
      W(RI02B) = K(RI02B) * C(SH2) * C(SC9H7)
      W(RI05) = K(RI05) * C(SC9H8) * C(SO)
      W(RI01F) = K(RI01F) * C(SC9H8)
      W(RI01B) = K(RI01B) * C(SH) * C(SC9H7)
      W(RI06F) = K(RI06F) * C(SC9H8) * C(SOH)
      W(RI06B) = K(RI06B) * C(SH2O) * C(SC9H7)
      W(RK110F) = K(RK110F) * C(SA2XC10H8) * C(SOH)
      W(RK110B) = K(RK110B) * C(SH2O) * C(SA2XXC10H7)
      W(RP108F) = K(RP108F) * C(SA2XC10H8)
      W(RP108B) = K(RP108B) * C(SH) * C(SA2XXC10H7)
      W(RK109F) = K(RK109F) * C(SA2XC10H8) * C(SH)
      W(RK109B) = K(RK109B) * C(SH2) * C(SA2XXC10H7)
      W(RG28F) = K(RG28F) * C(SSXCH2) * C(SN2)
      W(RG28B) = K(RG28B) * C(SN2) * C(STXCH2)


      CDOT(SN2) = - W(RG28F) + W(RG28F) - W(RG28B) & 
	    + W(RG28B)
      CDOT(SH) = - 2 * W(R5) + W(R2F) - W(R2B) & 
	    - W(R10) - W(R9F) + W(R9B) & 
	    + W(R3F) - W(R3B) - W(R12) & 
	    - W(R1F) + W(R1B) - W(R16) & 
	    - W(R15) - W(R13) + W(RG04) & 
	    + W(RG05) + W(RG06F) - W(RG06B) & 
	    + W(RG08) + W(RG35) + W(RG34F) & 
	    - W(RG34B) - W(RG30F) + W(RG30B) & 
	    + W(RG33) + 2 * W(RG19) + W(RG16) & 
	    + W(RG25) - W(RG14) + W(RG21) & 
	    + W(RG18F) - W(RG18B) + W(RG15) & 
	    + W(RG69) + W(RG52) + W(RG72) & 
	    - W(RG51F) + W(RG51B) + W(RG74F) & 
	    - W(RG74B) + W(RG53) + W(RG93) & 
	    - W(RG90F) + W(RG90B) + W(R28F)
      CDOT(SH) = CDOT(SH) - W(R28B) + W(R29F) - W(R29B) & 
	    - W(R32) + W(R37F) - W(R37B) & 
	    + W(R34) + W(R36F) - W(R36B) & 
	    + W(RG13) + W(RG50) - W(RG45F) & 
	    + W(RG45B) - W(RG86F) + W(RG86B) & 
	    - W(RG85F) + W(RG85B) + W(RG43) & 
	    - W(RG84) + W(RG108F) - W(RG108B) & 
	    + W(RG106) + W(RR017F) - W(RR017B) & 
	    - W(RG104) + W(RG121) - W(RG115F) & 
	    + W(RG115B) + W(RG116) + W(RR004F) & 
	    - W(RR004B) - W(RG130F) + W(RG130B) & 
	    - W(RG129F) + W(RG129B) - W(RG155F) & 
	    + W(RG155B) - W(RG156F) + W(RG156B) & 
	    - W(RG164) - W(RG174F) + W(RG174B) & 
	    + W(RG110) - W(RG109F) + W(RG109B)
      CDOT(SH) = CDOT(SH) - W(RG124F) + W(RG124B) - W(RG123F) & 
	    + W(RG123B) - W(RR055F) + W(RR055B) & 
	    + W(RP011F) - W(RP011B) - W(RR005F) & 
	    + W(RR005B) - W(RR078F) + W(RR078B) & 
	    - W(RR141F) + W(RR141B) - W(RIC66) & 
	    + W(RIC65) - W(RH10F) + W(RH10B) & 
	    - W(RH17) - W(RH19) + 2 * W(RLPCP18) & 
	    + 2 * W(RCP17) - W(RLPOX00) + W(ROX01F) & 
	    - W(ROX01B) - W(RLPOX03F) + W(RLPOX03B) & 
	    + W(RI03F) - W(RI03B) - W(RT04F) & 
	    + W(RT04B) + 2 * W(RT20) + W(RP106F) & 
	    - W(RP106B) + W(RP105F) - W(RP105B) & 
	    - W(RK017F) + W(RK017B) + W(RP016F) & 
	    - W(RP016B) + W(RK020F) - W(RK020B) & 
	    + W(RP104F) - W(RP104B) + 2 * W(RI18)
      CDOT(SH) = CDOT(SH) - W(RI02F) + W(RI02B) + W(RI01F) & 
	    - W(RI01B) + W(RP108F) - W(RP108B) & 
	    - W(RK109F) + W(RK109B)
      CDOT(SH2) = W(R5) - W(R2F) + W(R2B) & 
	    - W(R3F) + W(R3B) + W(R13) & 
	    - W(RG06F) + W(RG06B) - W(RG07) & 
	    + W(RG39) - W(RG34F) + W(RG34B) & 
	    + W(RG31) + W(RG30F) - W(RG30B) & 
	    - W(RG18F) + W(RG18B) + W(RG26) & 
	    + W(RG53) + W(RG90F) - W(RG90B) & 
	    + W(R32) + W(RG12) + W(RG45F) & 
	    - W(RG45B) + W(RG84) - W(RG108F) & 
	    + W(RG108B) + W(RG130F) - W(RG130B) & 
	    + W(RG156F) - W(RG156B) + W(RG174F) & 
	    - W(RG174B) + W(RG123F) - W(RG123B) & 
	    + W(RR078F) - W(RR078B) + W(RH19) & 
	    + W(RLPOX03F) - W(RLPOX03B) + W(RK017F) & 
	    - W(RK017B) + W(RI02F) - W(RI02B)
      CDOT(SH2) = CDOT(SH2) + W(RK109F) - W(RK109B)
      CDOT(SO) = - W(R2F) + W(R2B) - W(R10) & 
	    + W(R4F) - W(R4B) + W(R1F) & 
	    - W(R1B) - W(R17) + W(R15) & 
	    - W(RG04) + W(RG09) - W(RG31) & 
	    + W(RG20) - W(RG15) - W(RG52) & 
	    - W(RG53) - W(RG91F) + W(RG91B) & 
	    - W(R27) - W(R33) - W(R34) & 
	    - W(RG46) - W(RG105) + W(RG118F) & 
	    - W(RG118B) - W(RG117) - W(RG116) & 
	    - W(RG159) - W(RG158) - W(RG110) & 
	    - W(RR089) - W(RR079) - W(RR088) & 
	    - W(RR143) + W(RST13) - W(RI05)
      CDOT(SOH) = W(R2F) - W(R2B) + W(R10) & 
	    - W(R9F) + W(R9B) - W(R3F) & 
	    + W(R3B) - 2 * W(R4F) + 2 * W(R4B) & 
	    - 2 * W(R14F) + 2 * W(R14B) + W(R1F) & 
	    - W(R1B) - W(R19F) + W(R19B) & 
	    - W(R18F) + W(R18B) + 2 * W(R16) & 
	    + W(R17) - W(R26) - W(RG05) & 
	    + W(RG35) - W(RG33) - W(RG17F) & 
	    + W(RG17B) - W(RG16) + W(RG21) & 
	    - W(RG55F) + W(RG55B) + W(RG59) & 
	    - W(RG57F) + W(RG57B) + W(RG91F) & 
	    - W(RG91B) - W(RG92F) + W(RG92B) & 
	    - W(R28F) + W(R28B) + W(R31) & 
	    - W(R29F) + W(R29B) + W(R33) & 
	    - W(R35) + W(RG46) - W(RG47)
      CDOT(SOH) = CDOT(SOH) + W(RG85F) - W(RG85B) - W(RG106) & 
	    - W(RG118F) + W(RG118B) - W(RG121) & 
	    - W(RG119F) + W(RG119B) - W(RG122) & 
	    - W(RG132) - W(RG160F) + W(RG160B) & 
	    - W(RR024) - W(RG127F) + W(RG127B) & 
	    - W(RR025) - W(RR058) + W(RR079) & 
	    - W(RR080F) + W(RR080B) - W(RH12) & 
	    - W(RLPOX04) - W(RLPST00) - W(RP018F) & 
	    + W(RP018B) + W(RI05) - W(RI06F) & 
	    + W(RI06B) - W(RK110F) + W(RK110B)
      CDOT(SH2O) = W(R9F) - W(R9B) + W(R3F) & 
	    - W(R3B) + W(R4F) - W(R4B) & 
	    + W(R19F) - W(R19B) + W(R18F) & 
	    - W(R18B) + W(R15) + W(R26) & 
	    - W(RG08) - W(RG38F) + W(RG38F) & 
	    - W(RG38B) + W(RG38B) - W(RG39) & 
	    + W(RG36) + W(RG17F) - W(RG17B) & 
	    + W(RG55F) - W(RG55B) + W(RG57F) & 
	    - W(RG57B) + W(RG92F) - W(RG92B) & 
	    - W(R37F) + W(R37F) - W(R37B) & 
	    + W(R37B) + W(R35) + W(RG47) & 
	    + W(RG86F) - W(RG86B) + W(RG119F) & 
	    - W(RG119B) + W(RG132) + W(RG160F) & 
	    - W(RG160B) + W(RG127F) - W(RG127B) & 
	    + W(RR080F) - W(RR080B) + W(RLPOX04)
      CDOT(SH2O) = CDOT(SH2O) + W(RP018F) - W(RP018B) + W(RI06F) & 
	    - W(RI06B) + W(RK110F) - W(RK110B)
      CDOT(SH2O2) = W(R14F) - W(R14B) - W(R26)
      CDOT(SO2) = - W(R12) - W(R1F) + W(R1B) & 
	    + W(R19F) - W(R19B) + W(R18F) & 
	    - W(R18B) + W(R17) + W(R13) & 
	    - W(RG09) - W(RG35) - W(RG36) & 
	    - W(RG19) - W(RG20) - W(RG21) & 
	    - W(RG59) - W(R38) - W(RG89) & 
	    - W(RG107) - W(RG133) - W(RG135) & 
	    - W(RG171) - W(RR062) - W(RH07) & 
	    - W(RST13) + W(RI07)
      CDOT(SHO2) = W(R12) - W(R19F) + W(R19B) & 
	    - W(R18F) + W(R18B) - W(R16) & 
	    - W(R17) - W(R15) - W(R13) & 
	    + W(R26) - W(R31) + W(R38) & 
	    + W(RG89) + W(RG133) + W(RG171) & 
	    - W(RI07)
      CDOT(SCH) = - W(RG04) - W(RG05) - W(RG09) & 
	    - W(RG06F) + W(RG06B) - W(RG07) & 
	    - W(RG08) + W(RG30F) - W(RG30B) & 
	    + W(RG17F) - W(RG17B) - W(RG25) & 
	    - W(RG69) - W(RG93) - W(RG10F) & 
	    + W(RG10B) - W(RG50) - W(RG11) & 
	    + W(RG105) - W(RG112)
      CDOT(SCO) = W(RG04) + W(RG35) + W(RG31) & 
	    + W(RG36) + W(RG21) + W(RG53) & 
	    - W(RG24F) + W(RG24B) - W(R28F) & 
	    + W(R28B) - W(R31) - W(R27) & 
	    - W(R29F) + W(R29B) - W(RG10F) & 
	    + W(RG10B) - W(RG40F) + W(RG40F) & 
	    - W(RG40B) + W(RG40B) + W(R32) & 
	    + W(R37F) - W(R37B) + W(R38) & 
	    + W(R33) + W(R35) + W(R36F) & 
	    - W(R36B) + W(RG70) + W(RG12) & 
	    + W(RG42) + W(RG11) + W(RG107) & 
	    + W(RG105) + W(RG117) + W(RG122) & 
	    + W(RR011) + W(RR034) + 2 * W(RG110) & 
	    + W(RG109F) - W(RG109B) + 2 * W(RG114) & 
	    + W(RR009) + W(RG112) + W(RG113)
      CDOT(SCO) = CDOT(SCO) + W(RR022) + W(RG124F) - W(RG124B) & 
	    + W(RR025) + W(RR068) + W(RR058) & 
	    + W(RR089) + W(RH12) + W(RST13)
      CDOT(SHCO) = W(RG05) + W(RG09) + W(RG15) & 
	    - W(R32) - W(R37F) + W(R37B) & 
	    - W(R38) - W(R33) - W(R35) & 
	    - W(R34) - W(R36F) + W(R36B) & 
	    - W(RG70) + W(RG46) + W(RG13) & 
	    + W(RG71) + W(RG47) + W(RG45F) & 
	    - W(RG45B) + W(RG11) + W(RG107) & 
	    - W(RR011) + W(RG135) + W(RG159) & 
	    - W(RR034) + 2 * W(RR024) - W(RR068) & 
	    + W(RR062) + W(RR143)
      CDOT(STXCH2) = W(RG06F) - W(RG06B) + W(RG38F) & 
	    - W(RG38B) - W(RG19) - W(RG17F) & 
	    + W(RG17B) - W(RG20) - W(RG16) & 
	    - W(RG25) - W(RG14) - W(RG21) & 
	    - W(RG18F) + W(RG18B) - 2 * W(RG26) & 
	    - W(RG15) + W(RG55F) - W(RG55B) & 
	    - W(RG72) - W(RG24F) + W(RG24B) & 
	    + W(RG40F) - W(RG40B) + W(RG41F) & 
	    - W(RG41B) + W(RG117) + W(RG158) & 
	    - W(RG113) + W(RG28F) - W(RG28B)
      CDOT(SCH3) = W(RG07) + W(RG34F) - W(RG34B) & 
	    + W(RG14) + W(RG18F) - W(RG18B) & 
	    - W(RG69) - W(RG52) - W(RG55F) & 
	    + W(RG55B) - W(RG72) - W(RG51F) & 
	    + W(RG51B) - 2 * W(RG74F) + 2 * W(RG74B) & 
	    - W(RG53) - W(RG59) - W(RG57F) & 
	    + W(RG57B) + W(RG91F) - W(RG91B) & 
	    + W(RG92F) - W(RG92B) + W(RG90F) & 
	    - W(RG90B) - W(RG70) - W(RG71) & 
	    + W(RG85F) - W(RG85B) - W(RR017F) & 
	    + W(RR017B) + W(RG122) - W(RR013) & 
	    + W(RG159) - W(RIC67) - W(RG162F) & 
	    + W(RG162B) - W(RG178) + 2 * W(RG173F) & 
	    - 2 * W(RG173B) - W(RR022) + W(RG124F) & 
	    - W(RG124B) + W(RR005F) - W(RR005B)
      CDOT(SCH3) = CDOT(SCH3) + W(RR088) - W(RR081F) + W(RR081B) & 
	    + W(RR014F) - W(RR014B) + W(RR141F) & 
	    - W(RR141B) + W(RIC64) - W(RLPCP18) & 
	    + W(RT04F) - W(RT04B) - W(RI18) & 
	    - W(RI09F) + W(RI09B)
      CDOT(SCH2O) = W(RG08) + W(RG39) + W(RG33) & 
	    + W(RG20) + W(RG16) + W(RG52) & 
	    + W(RG59) - W(RG46) - W(RG13) & 
	    - W(RG71) - W(RG12) - W(RG47) & 
	    - W(RG50) - W(RG45F) + W(RG45B) & 
	    + W(RG89) + W(RG43) + W(RG84) & 
	    + W(RG42) + W(RG135) + W(RG158)
      CDOT(SSXCH2) = - W(RG38F) + W(RG38B) - W(RG39) & 
	    - W(RG35) - W(RG34F) + W(RG34B) & 
	    - W(RG31) - W(RG36) - W(RG30F) & 
	    + W(RG30B) - W(RG33) + W(RG57F) & 
	    - W(RG57B) - W(RG40F) + W(RG40B) & 
	    + W(RG86F) - W(RG86B) - W(RG41F) & 
	    + W(RG41B) - W(RG42) - W(RR004F) & 
	    + W(RR004B) + W(RG109F) - W(RG109B) & 
	    - W(RG28F) + W(RG28B)
      CDOT(SCO2) = W(RG19) + W(R28F) - W(R28B) & 
	    + W(R31) + W(R27) + W(R29F) & 
	    - W(R29B) + W(R34) - W(RG41F) & 
	    + W(RG41F) - W(RG41B) + W(RG41B) & 
	    - W(RG42) - W(RG11)
      CDOT(SC2H2) = W(RG25) + W(RG26) + W(RG108F) & 
	    - W(RG108B) + W(RG104) + W(RG118F) & 
	    - W(RG118B) - W(RG121) - W(RG115F) & 
	    + W(RG115B) - W(RG117) - W(RR008F) & 
	    + W(RR008B) - W(RG119F) + W(RG119B) & 
	    - W(RG116) - W(RR004F) + W(RR004B) & 
	    - W(RG122) + W(RG130F) - W(RG130B) & 
	    + W(RB07) + W(RR013) + W(RG132) & 
	    + W(RG133) + W(RG114) - W(RR009) & 
	    + W(RG112) - W(RCP16F) + W(RCP16B) & 
	    + W(RR005F) - W(RR005B) - W(ROX02F) & 
	    + W(ROX02B) + 2 * W(RH17) - W(RK012F) & 
	    + W(RK012B) - W(RI03F) + W(RI03B) & 
	    + W(RT05F) - W(RT05B) - W(RK100F) & 
	    + W(RK100B) - W(RP104F) + W(RP104B)
      CDOT(SC2H2) = CDOT(SC2H2) + W(RI25)
      CDOT(SC2H3) = W(RG69) + W(RG115F) - W(RG115B) & 
	    - W(RG130F) + W(RG130B) - 2 * W(RB07) & 
	    - W(RR011) - W(RR013) - W(RG132) & 
	    - W(RG133) - W(RG135) - W(RG129F) & 
	    + W(RG129B) + W(RG162F) - W(RG162B) & 
	    + W(RG160F) - W(RG160B) + W(RG156F) & 
	    - W(RG156B) + W(RG113) + W(RR014F) & 
	    - W(RR014B) + W(RLPP015) - W(RP105F) & 
	    + W(RP105B)
      CDOT(SC2H4) = W(RG72) + W(RG93) + W(RB07) & 
	    + W(RR011) + W(RG129F) - W(RG129B) & 
	    - W(RG159) - W(RIC67) - W(RG155F) & 
	    + W(RG155B) - W(RG162F) + W(RG162B) & 
	    - W(RG160F) + W(RG160B) - W(RG158) & 
	    - W(RG156F) + W(RG156B) + W(RG171) & 
	    + W(RR022) + W(RR058) + W(RR089) & 
	    + W(RR141F) - W(RR141B) + W(RIC64) & 
	    - W(RLPP015) - W(RP106F) + W(RP106B)
      CDOT(SCH4) = W(RG51F) - W(RG51B) - W(RG91F) & 
	    + W(RG91B) - W(RG93) - W(RG92F) & 
	    + W(RG92B) - W(RG90F) + W(RG90B) & 
	    + W(RG70) + W(RG71) + W(RR013) & 
	    + W(RG162F) - W(RG162B) + W(RG178) & 
	    + W(RR081F) - W(RR081B) + W(RI09F) & 
	    - W(RI09B)
      CDOT(SC2H5) = W(RG74F) - W(RG74B) + W(RG155F) & 
	    - W(RG155B) - W(RG171) - W(RR034) & 
	    - W(RG164) + W(RG178) + W(RG174F) & 
	    - W(RG174B) + W(RR143)
      CDOT(SCH2CO) = W(RG24F) - W(RG24B) + W(RG50) & 
	    + W(RG121) - W(RG124F) + W(RG124B) & 
	    - W(RG123F) + W(RG123B) - W(RG127F) & 
	    + W(RG127B) - W(RR025) + W(RR062)
      CDOT(SHCCO) = W(RG10F) - W(RG10B) + W(RG106) & 
	    + W(RG116) - W(RG110) - W(RG109F) & 
	    + W(RG109B) - 2 * W(RG114) - W(RR009) & 
	    - W(RG112) - W(RG113) - W(RR024) & 
	    - W(RR022) + W(RG123F) - W(RG123B) & 
	    + W(RG127F) - W(RG127B) + W(RR088) & 
	    + 2 * W(RH07) + W(RLPST00)
      CDOT(SCH2OH) = - W(RG89) - W(RG86F) + W(RG86B) & 
	    - W(RG85F) + W(RG85B) - W(RG43) & 
	    - W(RG84) + W(RR025)
      CDOT(SC2H) = - W(RG108F) + W(RG108B) - W(RG107) & 
	    - W(RG105) - W(RG106) - W(RR017F) & 
	    + W(RR017B) - W(RG104) - W(RG118F) & 
	    + W(RG118B) - W(RR008F) + W(RR008B) & 
	    + W(RG119F) - W(RG119B)
      CDOT(SC3H3) = W(RR017F) - W(RR017B) + W(RR004F) & 
	    - W(RR004B) + W(RR009) - W(RR068) & 
	    - 2 * W(RLPP009F) + 2 * W(RLPP009B) - W(RR055F) & 
	    + W(RR055B) - W(RR058) - W(RR062) & 
	    - 2 * W(RP011F) + 2 * W(RP011B) - W(RCP16F) & 
	    + W(RCP16B) + W(RR079) + W(RR080F) & 
	    - W(RR080B) + W(RR081F) - W(RR081B) & 
	    + W(RR078F) - W(RR078B) + W(RH12) & 
	    - W(RI00F) + W(RI00B) - W(RT20) & 
	    + W(RI25)
      CDOT(SNXC4H3) = W(RR008F) - W(RR008B) + W(RH10F) & 
	    - W(RH10B) - W(RH17) - W(RH19)
      CDOT(SIXC3H7) = W(RIC67) + W(RIC66) - W(RIC64) & 
	    - W(RIC65)
      CDOT(SC2H6) = W(RR034) + W(RG164) - W(RG178) & 
	    - W(RG174F) + W(RG174B) - W(RG173F) & 
	    + W(RG173B)
      CDOT(SPXC3H4) = W(RR068) + W(RR055F) - W(RR055B) & 
	    - W(RR089) - W(RR079) - W(RR080F) & 
	    + W(RR080B) - W(RR005F) + W(RR005B) & 
	    - W(RR088) - W(RR081F) + W(RR081B) & 
	    - W(RR078F) + W(RR078B)
      CDOT(SA1XC6H6) = W(RLPP009F) - W(RLPP009B) + W(RLPCP18) & 
	    + W(RLPP015) + W(RLPOX00) - W(RLPOX03F) & 
	    + W(RLPOX03B) - W(RLPOX04) + W(RLPST00)
      CDOT(SA1XXC6H5) = W(RP011F) - W(RP011B) - W(RLPP015) & 
	    - W(RI00F) + W(RI00B) - W(RLPOX00) & 
	    - W(RK012F) + W(RK012B) - W(ROX01F) & 
	    + W(ROX01B) + W(RLPOX03F) - W(RLPOX03B) & 
	    + W(RLPOX04) + W(RT04F) - W(RT04B)
      CDOT(SC5H5) = W(RCP16F) - W(RCP16B) - W(RLPCP18) & 
	    - 2 * W(RCP17) + W(RT05F) - W(RT05B) & 
	    + W(RI26)
      CDOT(SC3H6) = - W(RR014F) + W(RR014B) - W(RR143) & 
	    - W(RR141F) + W(RR141B) - W(RIC66) & 
	    + W(RIC65)
      CDOT(SC4H2) = - W(RH10F) + W(RH10B) - W(RH07) & 
	    - W(RH12) - W(ROX02F) + W(ROX02B) & 
	    + W(RH19) + W(RI26) + W(RI25)
      CDOT(SOXC6H4) = W(ROX02F) - W(ROX02B) + W(ROX01F) & 
	    - W(ROX01B)
      CDOT(SA2XC10H8) = W(RCP17) + W(RT20) + W(RP106F) & 
	    - W(RP106B) + W(RP105F) - W(RP105B) & 
	    + W(RP104F) - W(RP104B) + W(RI18) & 
	    - W(RK110F) + W(RK110B) - W(RP108F) & 
	    + W(RP108B) - W(RK109F) + W(RK109B)
      CDOT(SC9H8) = W(RI00F) - W(RI00B) + W(RI03F) & 
	    - W(RI03B) + W(RI07) - W(RI09F) & 
	    + W(RI09B) - W(RI02F) + W(RI02B) & 
	    - W(RI05) - W(RI01F) + W(RI01B) & 
	    - W(RI06F) + W(RI06B)
      CDOT(SA1C2H2XC8H7) = W(RK012F) - W(RK012B) - W(RK020F) & 
	    + W(RK020B) - W(RP104F) + W(RP104B) & 
	    - W(RST13)
      CDOT(SA1CH2XC7H7) = - W(RI03F) + W(RI03B) - W(RT04F) & 
	    + W(RT04B) - W(RT20) - W(RT05F) & 
	    + W(RT05B) + W(RST13)
      CDOT(SA1C2HYXC8H5) = - W(RP106F) + W(RP106B) - W(RK100F) & 
	    + W(RK100B) + W(RK017F) - W(RK017B) & 
	    + W(RP016F) - W(RP016B) + W(RP018F) & 
	    - W(RP018B)
      CDOT(SA2XXC10H7) = W(RK100F) - W(RK100B) + W(RK110F) & 
	    - W(RK110B) + W(RP108F) - W(RP108B) & 
	    + W(RK109F) - W(RK109B)
      CDOT(SA1C2HXC8H6) = - W(RP105F) + W(RP105B) - W(RLPST00) & 
	    - W(RK017F) + W(RK017B) - W(RP016F) & 
	    + W(RP016B) - W(RP018F) + W(RP018B) & 
	    + W(RK020F) - W(RK020B)
      CDOT(SC9H7) = - W(RI26) - W(RI07) - W(RI18) & 
	    - W(RI25) + W(RI09F) - W(RI09B) & 
	    + W(RI02F) - W(RI02B) + W(RI05) & 
	    + W(RI01F) - W(RI01B) + W(RI06F) & 
	    - W(RI06B)
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

      DOUBLE PRECISION MM(45)
      INCLUDE 'EthyRedF90.h'

      MM(SN2) =  2.80200000e+01
      MM(SH) =  1.00800000e+00
      MM(SH2) =  2.01600000e+00
      MM(SO) =  1.60000000e+01
      MM(SOH) =  1.70080000e+01
      MM(SH2O) =  1.80160000e+01
      MM(SH2O2) =  3.40160000e+01
      MM(SO2) =  3.20000000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SCH) =  1.30180000e+01
      MM(SCO) =  2.80100000e+01
      MM(SHCO) =  2.90180000e+01
      MM(STXCH2) =  1.40260000e+01
      MM(SCH3) =  1.50340000e+01
      MM(SCH2O) =  3.00260000e+01
      MM(SSXCH2) =  1.40260000e+01
      MM(SCO2) =  4.40100000e+01
      MM(SC2H2) =  2.60360000e+01
      MM(SC2H3) =  2.70440000e+01
      MM(SC2H4) =  2.80520000e+01
      MM(SCH4) =  1.60420000e+01
      MM(SC2H5) =  2.90600000e+01
      MM(SCH2CO) =  4.20360000e+01
      MM(SHCCO) =  4.10280000e+01
      MM(SCH2OH) =  3.10340000e+01
      MM(SC2H) =  2.50280000e+01
      MM(SC3H3) =  3.90540000e+01
      MM(SNXC4H3) =  5.10640000e+01
      MM(SIXC3H7) =  4.30860000e+01
      MM(SC2H6) =  3.00680000e+01
      MM(SPXC3H4) =  4.00620000e+01
      MM(SA1XC6H6) =  7.81080000e+01
      MM(SA1XXC6H5) =  7.71000000e+01
      MM(SC5H5) =  6.50900000e+01
      MM(SC3H6) =  4.20780000e+01
      MM(SC4H2) =  5.00560000e+01
      MM(SOXC6H4) =  7.60920000e+01
      MM(SA2XC10H8) =  1.28164000e+02
      MM(SC9H8) =  1.16154000e+02
      MM(SA1C2H2XC8H7) =  1.03136000e+02
      MM(SA1CH2XC7H7) =  9.11260000e+01
      MM(SA1C2HYXC8H5) =  1.01120000e+02
      MM(SA2XXC10H7) =  1.27156000e+02
      MM(SA1C2HXC8H6) =  1.02128000e+02
      MM(SC9H7) =  1.15146000e+02

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(45)
      INCLUDE 'EthyRedF90.h'
      
      NAMES(SN2)='N2'
      NAMES(SH)='H'
      NAMES(SH2)='H2'
      NAMES(SO)='O'
      NAMES(SOH)='OH'
      NAMES(SH2O)='H2O'
      NAMES(SH2O2)='H2O2'
      NAMES(SO2)='O2'
      NAMES(SHO2)='HO2'
      NAMES(SCH)='CH'
      NAMES(SCO)='CO'
      NAMES(SHCO)='HCO'
      NAMES(STXCH2)='TXCH2'
      NAMES(SCH3)='CH3'
      NAMES(SCH2O)='CH2O'
      NAMES(SSXCH2)='SXCH2'
      NAMES(SCO2)='CO2'
      NAMES(SC2H2)='C2H2'
      NAMES(SC2H3)='C2H3'
      NAMES(SC2H4)='C2H4'
      NAMES(SCH4)='CH4'
      NAMES(SC2H5)='C2H5'
      NAMES(SCH2CO)='CH2CO'
      NAMES(SHCCO)='HCCO'
      NAMES(SCH2OH)='CH2OH'
      NAMES(SC2H)='C2H'
      NAMES(SC3H3)='C3H3'
      NAMES(SNXC4H3)='NXC4H3'
      NAMES(SIXC3H7)='IXC3H7'
      NAMES(SC2H6)='C2H6'
      NAMES(SPXC3H4)='PXC3H4'
      NAMES(SA1XC6H6)='A1XC6H6'
      NAMES(SA1XXC6H5)='A1XXC6H5'
      NAMES(SC5H5)='C5H5'
      NAMES(SC3H6)='C3H6'
      NAMES(SC4H2)='C4H2'
      NAMES(SOXC6H4)='OXC6H4'
      NAMES(SA2XC10H8)='A2XC10H8'
      NAMES(SC9H8)='C9H8'
      NAMES(SA1C2H2XC8H7)='A1C2H2XC'
      NAMES(SA1CH2XC7H7)='A1CH2XC7'
      NAMES(SA1C2HYXC8H5)='A1C2HYXC'
      NAMES(SA2XXC10H7)='A2XXC10H'
      NAMES(SA1C2HXC8H6)='A1C2HXC8'
      NAMES(SC9H7)='C9H7'
      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'EthyRedF90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 279

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'EthyRedF90.h'

      NSPECIES_NONS = 45
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'EthyRedF90.h'
      real(DP) :: MUCOEFF(45)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SH2O2) =  1.30193418e-06
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SCH) =  1.27351520e-06
      MUCOEFF(SCO) =  1.06039632e-06
      MUCOEFF(SHCO) =  1.11568663e-06
      MUCOEFF(STXCH2) =  6.92304430e-07
      MUCOEFF(SCH3) =  7.16749611e-07
      MUCOEFF(SCH2O) =  1.13489904e-06
      MUCOEFF(SSXCH2) =  6.92304430e-07
      MUCOEFF(SCO2) =  1.25056029e-06
      MUCOEFF(SC2H2) =  8.10245830e-07
      MUCOEFF(SC2H3) =  8.25781476e-07
      MUCOEFF(SC2H4) =  8.96560348e-07
      MUCOEFF(SCH4) =  7.61887935e-07
      MUCOEFF(SC2H5) =  7.77507128e-07
      MUCOEFF(SCH2CO) =  1.09806251e-06
      MUCOEFF(SHCCO) =  2.73563116e-06
      MUCOEFF(SCH2OH) =  1.09210283e-06
      MUCOEFF(SC2H) =  7.94406422e-07
      MUCOEFF(SC3H3) =  7.36234631e-07
      MUCOEFF(SNXC4H3) =  7.10878342e-07
      MUCOEFF(SIXC3H7) =  7.57312844e-07
      MUCOEFF(SC2H6) =  7.90876817e-07
      MUCOEFF(SPXC3H4) =  7.45675363e-07
      MUCOEFF(SA1XC6H6) =  8.43012086e-07
      MUCOEFF(SA1XXC6H5) =  8.37554799e-07
      MUCOEFF(SC5H5) =  7.96430411e-07
      MUCOEFF(SC3H6) =  6.97617690e-07
      MUCOEFF(SC4H2) =  7.03827024e-07
      MUCOEFF(SOXC6H4) =  8.32061720e-07
      MUCOEFF(SA2XC10H8) =  7.91231307e-07
      MUCOEFF(SC9H8) =  7.53247192e-07
      MUCOEFF(SA1C2H2XC8H7) =  7.53008758e-07
      MUCOEFF(SA1CH2XC7H7) =  7.89808619e-07
      MUCOEFF(SA1C2HYXC8H5) =  8.20396614e-07
      MUCOEFF(SA2XXC10H7) =  7.88113678e-07
      MUCOEFF(SA1C2HXC8H6) =  8.24475477e-07
      MUCOEFF(SC9H7) =  7.49971680e-07

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'EthyRedF90.h'
      real(DP) :: KOVEREPS(45)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SH2O2) =  9.31098696e-03
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SCH) =  1.25000000e-02
      KOVEREPS(SCO) =  1.01936799e-02
      KOVEREPS(SHCO) =  2.00803213e-03
      KOVEREPS(STXCH2) =  6.94444444e-03
      KOVEREPS(SCH3) =  6.94444444e-03
      KOVEREPS(SCH2O) =  2.00803213e-03
      KOVEREPS(SSXCH2) =  6.94444444e-03
      KOVEREPS(SCO2) =  4.09836066e-03
      KOVEREPS(SC2H2) =  4.78468900e-03
      KOVEREPS(SC2H3) =  4.78468900e-03
      KOVEREPS(SC2H4) =  3.56125356e-03
      KOVEREPS(SCH4) =  7.07213579e-03
      KOVEREPS(SC2H5) =  3.96353547e-03
      KOVEREPS(SCH2CO) =  2.29357798e-03
      KOVEREPS(SHCCO) =  6.66666667e-03
      KOVEREPS(SCH2OH) =  2.39808153e-03
      KOVEREPS(SC2H) =  4.78468900e-03
      KOVEREPS(SC3H3) =  3.96825397e-03
      KOVEREPS(SNXC4H3) =  2.80112045e-03
      KOVEREPS(SIXC3H7) =  3.29597891e-03
      KOVEREPS(SC2H6) =  3.96353547e-03
      KOVEREPS(SPXC3H4) =  3.96825397e-03
      KOVEREPS(SA1XC6H6) =  2.15146299e-03
      KOVEREPS(SA1XXC6H5) =  2.15146299e-03
      KOVEREPS(SC5H5) =  2.50000000e-03
      KOVEREPS(SC3H6) =  3.74812594e-03
      KOVEREPS(SC4H2) =  2.80112045e-03
      KOVEREPS(SOXC6H4) =  2.15146299e-03
      KOVEREPS(SA2XC10H8) =  1.58629442e-03
      KOVEREPS(SC9H8) =  1.58629442e-03
      KOVEREPS(SA1C2H2XC8H7) =  1.83083120e-03
      KOVEREPS(SA1CH2XC7H7) =  2.01897840e-03
      KOVEREPS(SA1C2HYXC8H5) =  1.86706497e-03
      KOVEREPS(SA2XXC10H7) =  1.58629442e-03
      KOVEREPS(SA1C2HXC8H6) =  1.86706497e-03
      KOVEREPS(SC9H7) =  1.58629442e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'EthyRedF90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(279), C(45), M(7)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 45
!------------------------------------------------------------------
      implicit none
      include 'EthyRedF90.h'
      real(DP) :: H(45), CP(45), T

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
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.33727920e+00_DP + T * ( -2.47012365e-05_DP &
	   + T * (  1.66485593e-07_DP + T * ( -4.48915985e-11_DP &
	   + T * (  4.00510752e-15_DP ) ) ) ) ) -9.50158922e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.33727920e+00_DP + T * ( -4.94024731e-05_DP &
	   + T * (  4.99456778e-07_DP + T * ( -1.79566394e-10_DP &
	   + T * (  2.00255376e-14_DP ) ) ) ) )
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
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  3.03399249e+00_DP + T * (  1.08845902e-03_DP &
	   + T * ( -5.46908393e-08_DP + T * ( -2.42604967e-11_DP &
	   + T * (  3.36401984e-15_DP ) ) ) ) ) -3.00042971e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    3.03399249e+00_DP + T * (  2.17691804e-03_DP &
	   + T * ( -1.64072518e-07_DP + T * ( -9.70419870e-11_DP &
	   + T * (  1.68200992e-14_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.16500285e+00_DP + T * (  2.45415847e-03_DP &
	   + T * ( -6.33797417e-07_DP + T * (  9.27964965e-11_DP &
	   + T * ( -5.75816610e-15_DP ) ) ) ) ) -1.78617877e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.16500285e+00_DP + T * (  4.90831694e-03_DP &
	   + T * ( -1.90139225e-06_DP + T * (  3.71185986e-10_DP &
	   + T * ( -2.87908305e-14_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.28253784e+00_DP + T * (  7.41543770e-04_DP &
	   + T * ( -2.52655556e-07_DP + T * (  5.23676387e-11_DP &
	   + T * ( -4.33435588e-15_DP ) ) ) ) ) -1.08845772e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.28253784e+00_DP + T * (  1.48308754e-03_DP &
	   + T * ( -7.57966669e-07_DP + T * (  2.09470555e-10_DP &
	   + T * ( -2.16717794e-14_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.01721090e+00_DP + T * (  1.11991006e-03_DP &
	   + T * ( -2.11219383e-07_DP + T * (  2.85615925e-11_DP &
	   + T * ( -2.15817070e-15_DP ) ) ) ) ) +  1.11856713e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.01721090e+00_DP + T * (  2.23982013e-03_DP &
	   + T * ( -6.33658150e-07_DP + T * (  1.14246370e-10_DP &
	   + T * ( -1.07908535e-14_DP ) ) ) ) )
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
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  2.77217438e+00_DP + T * (  2.47847763e-03_DP &
	   + T * ( -8.28152043e-07_DP + T * (  1.47290445e-10_DP &
	   + T * ( -1.06701742e-14_DP ) ) ) ) ) +  4.01191815e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    2.77217438e+00_DP + T * (  4.95695526e-03_DP &
	   + T * ( -2.48445613e-06_DP + T * (  5.89161778e-10_DP &
	   + T * ( -5.33508711e-14_DP ) ) ) ) )
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
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.29203842e+00_DP + T * (  2.32794318e-03_DP &
	   + T * ( -6.70639823e-07_DP + T * (  1.04476500e-10_DP &
	   + T * ( -6.79432730e-15_DP ) ) ) ) ) +  5.09259997e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    2.29203842e+00_DP + T * (  4.65588637e-03_DP &
	   + T * ( -2.01191947e-06_DP + T * (  4.17906000e-10_DP &
	   + T * ( -3.39716365e-14_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  3.85746029e+00_DP + T * (  2.20718513e-03_DP &
	   + T * ( -7.38271347e-07_DP + T * (  1.30872547e-10_DP &
	   + T * ( -9.44168328e-15_DP ) ) ) ) ) -4.87591660e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    3.85746029e+00_DP + T * (  4.41437026e-03_DP &
	   + T * ( -2.21481404e-06_DP + T * (  5.23490188e-10_DP &
	   + T * ( -4.72084164e-14_DP ) ) ) ) )
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
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  2.03611116e+00_DP + T * (  7.32270755e-03_DP &
	   + T * ( -2.23692638e-06_DP + T * (  3.68057308e-10_DP &
	   + T * ( -2.51412122e-14_DP ) ) ) ) ) +  4.93988614e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    2.03611116e+00_DP + T * (  1.46454151e-02_DP &
	   + T * ( -6.71077915e-06_DP + T * (  1.47222923e-09_DP &
	   + T * ( -1.25706061e-13_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  1.65326226e+00_DP + T * (  5.01315495e-03_DP &
	   + T * ( -1.10553746e-06_DP + T * (  1.34120785e-10_DP &
	   + T * ( -6.29393516e-15_DP ) ) ) ) ) -1.00095936e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    1.65326226e+00_DP + T * (  1.00263099e-02_DP &
	   + T * ( -3.31661238e-06_DP + T * (  5.36483138e-10_DP &
	   + T * ( -3.14696758e-14_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  1.95465642e+00_DP + T * (  8.69863610e-03_DP &
	   + T * ( -2.66068889e-06_DP + T * (  4.38044223e-10_DP &
	   + T * ( -2.99283152e-14_DP ) ) ) ) ) +  1.28575200e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    1.95465642e+00_DP + T * (  1.73972722e-02_DP &
	   + T * ( -7.98206668e-06_DP + T * (  1.75217689e-09_DP &
	   + T * ( -1.49641576e-13_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  4.51129732e+00_DP + T * (  4.50179872e-03_DP &
	   + T * ( -1.38979878e-06_DP + T * (  2.30836470e-10_DP &
	   + T * ( -1.58967640e-14_DP ) ) ) ) ) -7.55105311e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    4.51129732e+00_DP + T * (  9.00359745e-03_DP &
	   + T * ( -4.16939635e-06_DP + T * (  9.23345882e-10_DP &
	   + T * ( -7.94838201e-14_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  5.62820580e+00_DP + T * (  2.04267005e-03_DP &
	   + T * ( -5.31151567e-07_DP + T * (  7.15651300e-11_DP &
	   + T * ( -3.88156640e-15_DP ) ) ) ) ) +  1.93272150e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    5.62820580e+00_DP + T * (  4.08534010e-03_DP &
	   + T * ( -1.59345470e-06_DP + T * (  2.86260520e-10_DP &
	   + T * ( -1.94078320e-14_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  3.69266569e+00_DP + T * (  4.32288399e-03_DP &
	   + T * ( -1.25033707e-06_DP + T * (  1.96808659e-10_DP &
	   + T * ( -1.29710840e-14_DP ) ) ) ) ) -3.24250627e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    3.69266569e+00_DP + T * (  8.64576797e-03_DP &
	   + T * ( -3.75101120e-06_DP + T * (  7.87234636e-10_DP &
	   + T * ( -6.48554201e-14_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  3.16780652e+00_DP + T * (  2.37610951e-03_DP &
	   + T * ( -6.12623590e-07_DP + T * (  7.60475630e-11_DP &
	   + T * ( -3.54465540e-15_DP ) ) ) ) ) +  6.71210650e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    3.16780652e+00_DP + T * (  4.75221902e-03_DP &
	   + T * ( -1.83787077e-06_DP + T * (  3.04190252e-10_DP &
	   + T * ( -1.77232770e-14_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  6.14915291e+00_DP + T * (  4.67031583e-03_DP &
	   + T * ( -1.25018451e-06_DP + T * (  1.72539079e-10_DP &
	   + T * ( -9.21649988e-15_DP ) ) ) ) ) +  3.83854848e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    6.14915291e+00_DP + T * (  9.34063166e-03_DP &
	   + T * ( -3.75055354e-06_DP + T * (  6.90156316e-10_DP &
	   + T * ( -4.60824994e-14_DP ) ) ) ) )
      H(SNXC4H3) =  1.62821949e+02_DP * ( &
	   T * (  7.25330164e+00_DP + T * (  5.97904230e-03_DP &
	   + T * ( -1.75571892e-06_DP + T * (  2.74954688e-10_DP &
	   + T * ( -1.76803350e-14_DP ) ) ) ) ) +  6.28977574e+04_DP )
      CP(SNXC4H3) =  1.62821949e+02_DP * ( &
	    7.25330164e+00_DP + T * (  1.19580846e-02_DP &
	   + T * ( -5.26715675e-06_DP + T * (  1.09981875e-09_DP &
	   + T * ( -8.84016751e-14_DP ) ) ) ) )
      H(SIXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  8.06336900e+00_DP + T * (  7.87244000e-03_DP &
	   + T * ( -1.72746400e-06_DP + T * (  1.86931125e-10_DP &
	   + T * ( -7.70884400e-15_DP ) ) ) ) ) +  5.31387100e+03_DP )
      CP(SIXC3H7) =  1.92970803e+02_DP * ( &
	    8.06336900e+00_DP + T * (  1.57448800e-02_DP &
	   + T * ( -5.18239200e-06_DP + T * (  7.47724500e-10_DP &
	   + T * ( -3.85442200e-14_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  1.07188150e+00_DP + T * (  1.08426339e-02_DP &
	   + T * ( -3.34186890e-06_DP + T * (  5.53530003e-10_DP &
	   + T * ( -3.80005780e-14_DP ) ) ) ) ) -1.14263932e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    1.07188150e+00_DP + T * (  2.16852677e-02_DP &
	   + T * ( -1.00256067e-05_DP + T * (  2.21412001e-09_DP &
	   + T * ( -1.90002890e-13_DP ) ) ) ) )
      H(SPXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  2.81460543e+00_DP + T * (  9.27622480e-03_DP &
	   + T * ( -3.18342256e-06_DP + T * (  5.99878425e-10_DP &
	   + T * ( -4.74970514e-14_DP ) ) ) ) ) +  2.07010771e+04_DP )
      CP(SPXC3H4) =  2.07536818e+02_DP * ( &
	    2.81460543e+00_DP + T * (  1.85524496e-02_DP &
	   + T * ( -9.55026768e-06_DP + T * (  2.39951370e-09_DP &
	   + T * ( -2.37485257e-13_DP ) ) ) ) )
      H(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   T * ( -2.06240612e-01_DP + T * (  2.32061220e-02_DP &
	   + T * ( -9.25511787e-06_DP + T * (  1.97227634e-09_DP &
	   + T * ( -1.72073052e-13_DP ) ) ) ) ) +  8.09883905e+03_DP )
      CP(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   -2.06240612e-01_DP + T * (  4.64122440e-02_DP &
	   + T * ( -2.77653536e-05_DP + T * (  7.88910537e-09_DP &
	   + T * ( -8.60365259e-13_DP ) ) ) ) )
      H(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   T * (  1.38016336e+00_DP + T * (  2.02016004e-02_DP &
	   + T * ( -8.07502950e-06_DP + T * (  1.72180830e-09_DP &
	   + T * ( -1.50192160e-13_DP ) ) ) ) ) +  3.86973520e+04_DP )
      CP(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	    1.38016336e+00_DP + T * (  4.04032009e-02_DP &
	   + T * ( -2.42250885e-05_DP + T * (  6.88723321e-09_DP &
	   + T * ( -7.50960802e-13_DP ) ) ) ) )
      H(SC5H5) =  1.27736058e+02_DP * ( &
	   T * (  4.21464919e+00_DP + T * (  1.35917364e-02_DP &
	   + T * ( -4.43910697e-06_DP + T * (  7.72450297e-10_DP &
	   + T * ( -5.55759746e-14_DP ) ) ) ) ) +  2.88952416e+04_DP )
      CP(SC5H5) =  1.27736058e+02_DP * ( &
	    4.21464919e+00_DP + T * (  2.71834728e-02_DP &
	   + T * ( -1.33173209e-05_DP + T * (  3.08980119e-09_DP &
	   + T * ( -2.77879873e-13_DP ) ) ) ) )
      H(SC3H6) =  1.97593517e+02_DP * ( &
	   T * (  4.71697982e-01_DP + T * (  1.44756535e-02_DP &
	   + T * ( -5.22006063e-06_DP + T * (  1.02860800e-09_DP &
	   + T * ( -8.46150282e-14_DP ) ) ) ) ) +  1.12603387e+03_DP )
      CP(SC3H6) =  1.97593517e+02_DP * ( &
	    4.71697982e-01_DP + T * (  2.89513070e-02_DP &
	   + T * ( -1.56601819e-05_DP + T * (  4.11443199e-09_DP &
	   + T * ( -4.23075141e-13_DP ) ) ) ) )
      H(SC4H2) =  1.66100767e+02_DP * ( &
	   T * (  9.75839793e+00_DP + T * (  1.89436611e-03_DP &
	   + T * (  1.02047338e-07_DP + T * ( -1.58413756e-10_DP &
	   + T * (  2.25860064e-14_DP ) ) ) ) ) +  5.22698696e+04_DP )
      CP(SC4H2) =  1.66100767e+02_DP * ( &
	    9.75839793e+00_DP + T * (  3.78873223e-03_DP &
	   + T * (  3.06142015e-07_DP + T * ( -6.33655024e-10_DP &
	   + T * (  1.12930032e-13_DP ) ) ) ) )
      H(SOXC6H4) =  1.09266940e+02_DP * ( &
	   T * (  2.98618725e+00_DP + T * (  1.68818922e-02_DP &
	   + T * ( -6.67461303e-06_DP + T * (  1.40963421e-09_DP &
	   + T * ( -1.22000829e-13_DP ) ) ) ) ) +  5.12231321e+04_DP )
      CP(SOXC6H4) =  1.09266940e+02_DP * ( &
	    2.98618725e+00_DP + T * (  3.37637843e-02_DP &
	   + T * ( -2.00238391e-05_DP + T * (  5.63853682e-09_DP &
	   + T * ( -6.10004145e-13_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * (  1.76826275e+00_DP + T * (  3.44571753e-02_DP &
	   + T * ( -1.38107392e-05_DP + T * (  2.94785772e-09_DP &
	   + T * ( -2.57194122e-13_DP ) ) ) ) ) +  1.45412795e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	    1.76826275e+00_DP + T * (  6.89143506e-02_DP &
	   + T * ( -4.14322176e-05_DP + T * (  1.17914309e-08_DP &
	   + T * ( -1.28597061e-12_DP ) ) ) ) )
      H(SC9H8) =  7.15803158e+01_DP * ( &
	   T * (  1.15459802e+00_DP + T * (  3.27112098e-02_DP &
	   + T * ( -1.30835036e-05_DP + T * (  2.78922102e-09_DP &
	   + T * ( -2.43185350e-13_DP ) ) ) ) ) +  1.68166108e+04_DP )
      CP(SC9H8) =  7.15803158e+01_DP * ( &
	    1.15459802e+00_DP + T * (  6.54224196e-02_DP &
	   + T * ( -3.92505107e-05_DP + T * (  1.11568841e-08_DP &
	   + T * ( -1.21592675e-12_DP ) ) ) ) )
      H(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   T * (  5.98044803e+00_DP + T * (  2.34715873e-02_DP &
	   + T * ( -8.91261580e-06_DP + T * (  1.82444488e-09_DP &
	   + T * ( -1.54218606e-13_DP ) ) ) ) ) +  4.32864172e+04_DP )
      CP(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	    5.98044803e+00_DP + T * (  4.69431747e-02_DP &
	   + T * ( -2.67378474e-05_DP + T * (  7.29777950e-09_DP &
	   + T * ( -7.71093028e-13_DP ) ) ) ) )
      H(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   T * (  3.30049696e+00_DP + T * (  2.40027670e-02_DP &
	   + T * ( -9.28143407e-06_DP + T * (  1.93092839e-09_DP &
	   + T * ( -1.65430827e-13_DP ) ) ) ) ) +  2.17498572e+04_DP )
      CP(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	    3.30049696e+00_DP + T * (  4.80055340e-02_DP &
	   + T * ( -2.78443022e-05_DP + T * (  7.72371356e-09_DP &
	   + T * ( -8.27154136e-13_DP ) ) ) ) )
      H(SA1C2HYXC8H5) =  8.22225079e+01_DP * ( &
	   T * (  7.23812069e+00_DP + T * (  1.91906054e-02_DP &
	   + T * ( -7.29502437e-06_DP + T * (  1.49290312e-09_DP &
	   + T * ( -1.26070293e-13_DP ) ) ) ) ) +  6.49528135e+04_DP )
      CP(SA1C2HYXC8H5) =  8.22225079e+01_DP * ( &
	    7.23812069e+00_DP + T * (  3.83812109e-02_DP &
	   + T * ( -2.18850731e-05_DP + T * (  5.97161247e-09_DP &
	   + T * ( -6.30351467e-13_DP ) ) ) ) )
      H(SA2XXC10H7) =  6.53869263e+01_DP * ( &
	   T * (  3.22892303e+00_DP + T * (  3.15632243e-02_DP &
	   + T * ( -1.26860794e-05_DP + T * (  2.71135172e-09_DP &
	   + T * ( -2.36685024e-13_DP ) ) ) ) ) +  4.78400840e+04_DP )
      CP(SA2XXC10H7) =  6.53869263e+01_DP * ( &
	    3.22892303e+00_DP + T * (  6.31264486e-02_DP &
	   + T * ( -3.80582381e-05_DP + T * (  1.08454069e-08_DP &
	   + T * ( -1.18342512e-12_DP ) ) ) ) )
      H(SA1C2HXC8H6) =  8.14109745e+01_DP * ( &
	   T * (  5.81520488e+00_DP + T * (  2.20436466e-02_DP &
	   + T * ( -8.40179527e-06_DP + T * (  1.72568807e-09_DP &
	   + T * ( -1.46275782e-13_DP ) ) ) ) ) +  3.30271906e+04_DP )
      CP(SA1C2HXC8H6) =  8.14109745e+01_DP * ( &
	    5.81520488e+00_DP + T * (  4.40872933e-02_DP &
	   + T * ( -2.52053858e-05_DP + T * (  6.90275228e-09_DP &
	   + T * ( -7.31378908e-13_DP ) ) ) ) )
      H(SC9H7) =  7.22069373e+01_DP * ( &
	   T * (  3.65597547e+00_DP + T * (  2.87404231e-02_DP &
	   + T * ( -1.14290200e-05_DP + T * (  2.42569698e-09_DP &
	   + T * ( -2.10772824e-13_DP ) ) ) ) ) +  3.06843457e+04_DP )
      CP(SC9H7) =  7.22069373e+01_DP * ( &
	    3.65597547e+00_DP + T * (  5.74808463e-02_DP &
	   + T * ( -3.42870600e-05_DP + T * (  9.70278793e-09_DP &
	   + T * ( -1.05386412e-12_DP ) ) ) ) )
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
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.34433112e+00_DP + T * (  3.99026037e-03_DP &
	   + T * ( -6.49271700e-06_DP + T * (  5.03930235e-09_DP &
	   + T * ( -1.47522352e-12_DP ) ) ) ) ) -9.17935173e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.34433112e+00_DP + T * (  7.98052075e-03_DP &
	   + T * ( -1.94781510e-05_DP + T * (  2.01572094e-08_DP &
	   + T * ( -7.37611761e-12_DP ) ) ) ) )
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
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  4.19864056e+00_DP + T * ( -1.01821705e-03_DP &
	   + T * (  2.17346737e-06_DP + T * ( -1.37199266e-09_DP &
	   + T * (  3.54395634e-13_DP ) ) ) ) ) -3.02937267e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    4.19864056e+00_DP + T * ( -2.03643410e-03_DP &
	   + T * (  6.52040211e-06_DP + T * ( -5.48797062e-09_DP &
	   + T * (  1.77197817e-12_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.27611269e+00_DP + T * ( -2.71411208e-04_DP &
	   + T * (  5.57785670e-06_DP + T * ( -5.39427032e-09_DP &
	   + T * (  1.72490873e-12_DP ) ) ) ) ) -1.77025821e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.27611269e+00_DP + T * ( -5.42822417e-04_DP &
	   + T * (  1.67335701e-05_DP + T * ( -2.15770813e-08_DP &
	   + T * (  8.62454363e-12_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.78245636e+00_DP + T * ( -1.49836708e-03_DP &
	   + T * (  3.28243400e-06_DP + T * ( -2.42032377e-09_DP &
	   + T * (  6.48745674e-13_DP ) ) ) ) ) -1.06394356e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.78245636e+00_DP + T * ( -2.99673416e-03_DP &
	   + T * (  9.84730201e-06_DP + T * ( -9.68129509e-09_DP &
	   + T * (  3.24372837e-12_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.30179801e+00_DP + T * ( -2.37456025e-03_DP &
	   + T * (  7.05276303e-06_DP + T * ( -6.06909735e-09_DP &
	   + T * (  1.85845025e-12_DP ) ) ) ) ) +  2.94808040e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.30179801e+00_DP + T * ( -4.74912051e-03_DP &
	   + T * (  2.11582891e-05_DP + T * ( -2.42763894e-08_DP &
	   + T * (  9.29225124e-12_DP ) ) ) ) )
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
      H(SHCO) =  2.86523537e+02_DP * ( &
	   T * (  4.22118584e+00_DP + T * ( -1.62196266e-03_DP &
	   + T * (  4.59331487e-06_DP + T * ( -3.32860233e-09_DP &
	   + T * (  8.67537730e-13_DP ) ) ) ) ) +  3.83956496e+03_DP )
      CP(SHCO) =  2.86523537e+02_DP * ( &
	    4.22118584e+00_DP + T * ( -3.24392532e-03_DP &
	   + T * (  1.37799446e-05_DP + T * ( -1.33144093e-08_DP &
	   + T * (  4.33768865e-12_DP ) ) ) ) )
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
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  4.19860411e+00_DP + T * ( -1.18330710e-03_DP &
	   + T * (  2.74432073e-06_DP + T * ( -1.67203995e-09_DP &
	   + T * (  3.88629474e-13_DP ) ) ) ) ) +  5.04968163e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    4.19860411e+00_DP + T * ( -2.36661419e-03_DP &
	   + T * (  8.23296220e-06_DP + T * ( -6.68815981e-09_DP &
	   + T * (  1.94314737e-12_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  2.35677352e+00_DP + T * (  4.49229839e-03_DP &
	   + T * ( -2.37452090e-06_DP + T * (  6.14797555e-10_DP &
	   + T * ( -2.87399096e-14_DP ) ) ) ) ) -4.83719697e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    2.35677352e+00_DP + T * (  8.98459677e-03_DP &
	   + T * ( -7.12356269e-06_DP + T * (  2.45919022e-09_DP &
	   + T * ( -1.43699548e-13_DP ) ) ) ) )
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
      H(SC2H4) =  2.96390275e+02_DP * ( &
	   T * (  3.95920148e+00_DP + T * ( -3.78526124e-03_DP &
	   + T * (  1.90330097e-05_DP + T * ( -1.72897188e-08_DP &
	   + T * (  5.39768746e-12_DP ) ) ) ) ) +  5.08977593e+03_DP )
      CP(SC2H4) =  2.96390275e+02_DP * ( &
	    3.95920148e+00_DP + T * ( -7.57052247e-03_DP &
	   + T * (  5.70990292e-05_DP + T * ( -6.91588753e-08_DP &
	   + T * (  2.69884373e-11_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  5.14911468e+00_DP + T * ( -6.83110045e-03_DP &
	   + T * (  1.63817974e-05_DP + T * ( -1.21061692e-08_DP &
	   + T * (  3.33206882e-12_DP ) ) ) ) ) -1.02465983e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    5.14911468e+00_DP + T * ( -1.36622009e-02_DP &
	   + T * (  4.91453921e-05_DP + T * ( -4.84246767e-08_DP &
	   + T * (  1.66603441e-11_DP ) ) ) ) )
      H(SC2H5) =  2.86109429e+02_DP * ( &
	   T * (  4.30646568e+00_DP + T * ( -2.09329446e-03_DP &
	   + T * (  1.65714269e-05_DP + T * ( -1.49781651e-08_DP &
	   + T * (  4.61018008e-12_DP ) ) ) ) ) +  1.28416265e+04_DP )
      CP(SC2H5) =  2.86109429e+02_DP * ( &
	    4.30646568e+00_DP + T * ( -4.18658892e-03_DP &
	   + T * (  4.97142807e-05_DP + T * ( -5.99126606e-08_DP &
	   + T * (  2.30509004e-11_DP ) ) ) ) )
      H(SCH2CO) =  1.97790941e+02_DP * ( &
	   T * (  2.13583630e+00_DP + T * (  9.05943605e-03_DP &
	   + T * ( -5.79824913e-06_DP + T * (  2.33599392e-09_DP &
	   + T * ( -4.02915230e-13_DP ) ) ) ) ) -7.04291804e+03_DP )
      CP(SCH2CO) =  1.97790941e+02_DP * ( &
	    2.13583630e+00_DP + T * (  1.81188721e-02_DP &
	   + T * ( -1.73947474e-05_DP + T * (  9.34397568e-09_DP &
	   + T * ( -2.01457615e-12_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  2.25172140e+00_DP + T * (  8.82751050e-03_DP &
	   + T * ( -7.90970033e-06_DP + T * (  4.31893975e-09_DP &
	   + T * ( -1.01329622e-12_DP ) ) ) ) ) +  2.00594490e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    2.25172140e+00_DP + T * (  1.76550210e-02_DP &
	   + T * ( -2.37291010e-05_DP + T * (  1.72757590e-08_DP &
	   + T * ( -5.06648110e-12_DP ) ) ) ) )
      H(SCH2OH) =  2.67910679e+02_DP * ( &
	   T * (  3.86388918e+00_DP + T * (  2.79836152e-03_DP &
	   + T * (  1.97757264e-06_DP + T * ( -2.61330030e-09_DP &
	   + T * (  8.73934556e-13_DP ) ) ) ) ) -3.19391367e+03_DP )
      CP(SCH2OH) =  2.67910679e+02_DP * ( &
	    3.86388918e+00_DP + T * (  5.59672304e-03_DP &
	   + T * (  5.93271791e-06_DP + T * ( -1.04532012e-08_DP &
	   + T * (  4.36967278e-12_DP ) ) ) ) )
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  2.88965733e+00_DP + T * (  6.70498055e-03_DP &
	   + T * ( -9.49231670e-06_DP + T * (  7.36977613e-09_DP &
	   + T * ( -2.18663022e-12_DP ) ) ) ) ) +  6.68393932e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    2.88965733e+00_DP + T * (  1.34099611e-02_DP &
	   + T * ( -2.84769501e-05_DP + T * (  2.94791045e-08_DP &
	   + T * ( -1.09331511e-11_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  1.40299238e+00_DP + T * (  1.50886664e-02_DP &
	   + T * ( -1.32816458e-05_DP + T * (  7.33836572e-09_DP &
	   + T * ( -1.74110916e-12_DP ) ) ) ) ) +  3.93108220e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    1.40299238e+00_DP + T * (  3.01773327e-02_DP &
	   + T * ( -3.98449373e-05_DP + T * (  2.93534629e-08_DP &
	   + T * ( -8.70554579e-12_DP ) ) ) ) )
      H(SNXC4H3) =  1.62821949e+02_DP * ( &
	   T * ( -3.55175031e-02_DP + T * (  2.15254251e-02_DP &
	   + T * ( -1.91909716e-05_DP + T * (  1.03970785e-08_DP &
	   + T * ( -2.41501714e-12_DP ) ) ) ) ) +  6.43506593e+04_DP )
      CP(SNXC4H3) =  1.62821949e+02_DP * ( &
	   -3.55175031e-02_DP + T * (  4.30508503e-02_DP &
	   + T * ( -5.75729147e-05_DP + T * (  4.15883142e-08_DP &
	   + T * ( -1.20750857e-11_DP ) ) ) ) )
      H(SIXC3H7) =  1.92970803e+02_DP * ( &
	   T * (  1.71330000e+00_DP + T * (  1.27130800e-02_DP &
	   + T * (  5.26936000e-07_DP + T * ( -4.55321500e-09_DP &
	   + T * (  1.76554200e-12_DP ) ) ) ) ) +  7.53580900e+03_DP )
      CP(SIXC3H7) =  1.92970803e+02_DP * ( &
	    1.71330000e+00_DP + T * (  2.54261600e-02_DP &
	   + T * (  1.58080800e-06_DP + T * ( -1.82128600e-08_DP &
	   + T * (  8.82771000e-12_DP ) ) ) ) )
      H(SC2H6) =  2.76517893e+02_DP * ( &
	   T * (  4.29142492e+00_DP + T * ( -2.75077135e-03_DP &
	   + T * (  1.99812763e-05_DP + T * ( -1.77116571e-08_DP &
	   + T * (  5.37371542e-12_DP ) ) ) ) ) -1.15222055e+04_DP )
      CP(SC2H6) =  2.76517893e+02_DP * ( &
	    4.29142492e+00_DP + T * ( -5.50154270e-03_DP &
	   + T * (  5.99438288e-05_DP + T * ( -7.08466285e-08_DP &
	   + T * (  2.68685771e-11_DP ) ) ) ) )
      H(SPXC3H4) =  2.07536818e+02_DP * ( &
	   T * (  1.46175323e+00_DP + T * (  1.23013301e-02_DP &
	   + T * ( -6.34064650e-06_DP + T * (  2.15090855e-09_DP &
	   + T * ( -3.33458480e-13_DP ) ) ) ) ) +  2.09209793e+04_DP )
      CP(SPXC3H4) =  2.07536818e+02_DP * ( &
	    1.46175323e+00_DP + T * (  2.46026602e-02_DP &
	   + T * ( -1.90219395e-05_DP + T * (  8.60363422e-09_DP &
	   + T * ( -1.66729240e-12_DP ) ) ) ) )
      H(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   T * ( -5.51558393e+00_DP + T * (  3.22726613e-02_DP &
	   + T * ( -1.47134309e-05_DP + T * (  1.86928040e-09_DP &
	   + T * (  6.20564508e-13_DP ) ) ) ) ) +  9.11031457e+03_DP )
      CP(SA1XC6H6) =  1.06446715e+02_DP * ( &
	   -5.51558393e+00_DP + T * (  6.45453225e-02_DP &
	   + T * ( -4.41402928e-05_DP + T * (  7.47712161e-09_DP &
	   + T * (  3.10282254e-12_DP ) ) ) ) )
      H(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   T * ( -4.87654845e+00_DP + T * (  3.13402891e-02_DP &
	   + T * ( -1.62467429e-05_DP + T * (  3.52805717e-09_DP &
	   + T * (  1.03703662e-13_DP ) ) ) ) ) +  3.99269438e+04_DP )
      CP(SA1XXC6H5) =  1.07838392e+02_DP * ( &
	   -4.87654845e+00_DP + T * (  6.26805782e-02_DP &
	   + T * ( -4.87402286e-05_DP + T * (  1.41122287e-08_DP &
	   + T * (  5.18518312e-13_DP ) ) ) ) )
      H(SC5H5) =  1.27736058e+02_DP * ( &
	   T * ( -7.37844042e+00_DP + T * (  4.86195909e-02_DP &
	   + T * ( -5.65263793e-05_DP + T * (  3.79546668e-08_DP &
	   + T * ( -1.02415096e-11_DP ) ) ) ) ) +  3.05514662e+04_DP )
      CP(SC5H5) =  1.27736058e+02_DP * ( &
	   -7.37844042e+00_DP + T * (  9.72391818e-02_DP &
	   + T * ( -1.69579138e-04_DP + T * (  1.51818667e-07_DP &
	   + T * ( -5.12075479e-11_DP ) ) ) ) )
      H(SC3H6) =  1.97593517e+02_DP * ( &
	   T * ( -2.29261670e-03_DP + T * (  1.55130533e-02_DP &
	   + T * ( -5.57171827e-06_DP + T * (  4.73985425e-10_DP &
	   + T * (  2.49915830e-13_DP ) ) ) ) ) +  1.13437406e+03_DP )
      CP(SC3H6) =  1.97593517e+02_DP * ( &
	   -2.29261670e-03_DP + T * (  3.10261065e-02_DP &
	   + T * ( -1.67151548e-05_DP + T * (  1.89594170e-09_DP &
	   + T * (  1.24957915e-12_DP ) ) ) ) )
      H(SC4H2) =  1.66100767e+02_DP * ( &
	   T * (  1.73325212e-01_DP + T * (  2.26974515e-02_DP &
	   + T * ( -2.43374610e-05_DP + T * (  1.48812934e-08_DP &
	   + T * ( -3.74969432e-12_DP ) ) ) ) ) +  5.42239385e+04_DP )
      CP(SC4H2) =  1.66100767e+02_DP * ( &
	    1.73325212e-01_DP + T * (  4.53949030e-02_DP &
	   + T * ( -7.30123830e-05_DP + T * (  5.95251736e-08_DP &
	   + T * ( -1.87484716e-11_DP ) ) ) ) )
      H(SOXC6H4) =  1.09266940e+02_DP * ( &
	   T * ( -3.46229657e+00_DP + T * (  2.87008288e-02_DP &
	   + T * ( -1.64328123e-05_DP + T * (  4.76701207e-09_DP &
	   + T * ( -3.90861418e-13_DP ) ) ) ) ) +  5.25223614e+04_DP )
      CP(SOXC6H4) =  1.09266940e+02_DP * ( &
	   -3.46229657e+00_DP + T * (  5.74016575e-02_DP &
	   + T * ( -4.92984369e-05_DP + T * (  1.90680483e-08_DP &
	   + T * ( -1.95430709e-12_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * ( -8.72434585e+00_DP + T * (  5.26880040e-02_DP &
	   + T * ( -2.67236897e-05_DP + T * (  5.46364935e-09_DP &
	   + T * (  2.84133212e-13_DP ) ) ) ) ) +  1.66588912e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   -8.72434585e+00_DP + T * (  1.05376008e-01_DP &
	   + T * ( -8.01710690e-05_DP + T * (  2.18545974e-08_DP &
	   + T * (  1.42066606e-12_DP ) ) ) ) )
      H(SC9H8) =  7.15803158e+01_DP * ( &
	   T * ( -8.12447817e+00_DP + T * (  4.88828533e-02_DP &
	   + T * ( -2.43478658e-05_DP + T * (  4.70737525e-09_DP &
	   + T * (  3.68066426e-13_DP ) ) ) ) ) +  1.86589996e+04_DP )
      CP(SC9H8) =  7.15803158e+01_DP * ( &
	   -8.12447817e+00_DP + T * (  9.77657067e-02_DP &
	   + T * ( -7.30435974e-05_DP + T * (  1.88295010e-08_DP &
	   + T * (  1.84033213e-12_DP ) ) ) ) )
      H(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   T * ( -6.30997035e+00_DP + T * (  4.75453915e-02_DP &
	   + T * ( -3.18566445e-05_DP + T * (  1.24202002e-08_DP &
	   + T * ( -2.03584362e-12_DP ) ) ) ) ) +  4.57329298e+04_DP )
      CP(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   -6.30997035e+00_DP + T * (  9.50907829e-02_DP &
	   + T * ( -9.55699336e-05_DP + T * (  4.96808010e-08_DP &
	   + T * ( -1.01792181e-11_DP ) ) ) ) )
      H(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   T * ( -6.07053038e+00_DP + T * (  4.17600754e-02_DP &
	   + T * ( -2.47233361e-05_DP + T * (  7.82884618e-09_DP &
	   + T * ( -8.47341736e-13_DP ) ) ) ) ) +  2.35894712e+04_DP )
      CP(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   -6.07053038e+00_DP + T * (  8.35201507e-02_DP &
	   + T * ( -7.41700083e-05_DP + T * (  3.13153847e-08_DP &
	   + T * ( -4.23670868e-12_DP ) ) ) ) )
      H(SA1C2HYXC8H5) =  8.22225079e+01_DP * ( &
	   T * ( -4.42757639e+00_DP + T * (  4.18334322e-02_DP &
	   + T * ( -2.90035454e-05_DP + T * (  1.17571415e-08_DP &
	   + T * ( -2.03633970e-12_DP ) ) ) ) ) +  6.73302359e+04_DP )
      CP(SA1C2HYXC8H5) =  8.22225079e+01_DP * ( &
	   -4.42757639e+00_DP + T * (  8.36668645e-02_DP &
	   + T * ( -8.70106362e-05_DP + T * (  4.70285661e-08_DP &
	   + T * ( -1.01816985e-11_DP ) ) ) ) )
      H(SA2XXC10H7) =  6.53869263e+01_DP * ( &
	   T * ( -8.02718034e+00_DP + T * (  5.14622590e-02_DP &
	   + T * ( -2.78090670e-05_DP + T * (  6.80338458e-09_DP &
	   + T * ( -1.44911911e-13_DP ) ) ) ) ) +  5.01363344e+04_DP )
      CP(SA2XXC10H7) =  6.53869263e+01_DP * ( &
	   -8.02718034e+00_DP + T * (  1.02924518e-01_DP &
	   + T * ( -8.34272010e-05_DP + T * (  2.72135383e-08_DP &
	   + T * ( -7.24559554e-13_DP ) ) ) ) )
      H(SA1C2HXC8H6) =  8.14109745e+01_DP * ( &
	   T * ( -5.21036925e+00_DP + T * (  4.32775972e-02_DP &
	   + T * ( -2.81669161e-05_DP + T * (  1.05480176e-08_DP &
	   + T * ( -1.63353233e-12_DP ) ) ) ) ) +  3.52488620e+04_DP )
      CP(SA1C2HXC8H6) =  8.14109745e+01_DP * ( &
	   -5.21036925e+00_DP + T * (  8.65551944e-02_DP &
	   + T * ( -8.45007483e-05_DP + T * (  4.21920706e-08_DP &
	   + T * ( -8.16766167e-12_DP ) ) ) ) )
      H(SC9H7) =  7.22069373e+01_DP * ( &
	   T * ( -8.73685384e+00_DP + T * (  5.17108180e-02_DP &
	   + T * ( -3.07807798e-05_DP + T * (  9.39057395e-09_DP &
	   + T * ( -8.81210540e-13_DP ) ) ) ) ) +  3.31641009e+04_DP )
      CP(SC9H7) =  7.22069373e+01_DP * ( &
	   -8.73685384e+00_DP + T * (  1.03421636e-01_DP &
	   + T * ( -9.23423393e-05_DP + T * (  3.75622958e-08_DP &
	   + T * ( -4.40605270e-12_DP ) ) ) ) )
      END IF

      END
