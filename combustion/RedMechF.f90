!----------------------------------------------------------
! ======= RedMechF.f90 =======
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
      include 'RedMechF90.h'
      real(DP) :: CDOT(47), W(290), K(290), &
      C(47), M(6), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM3) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 2 * C(SH2) + C(SH) &
	    + C(SOH) + 12 * C(SH2O) &
	    + C(SO2) + C(SHO2) &
	    + C(SCH) + 1.75 * C(SCO) &
	    + C(SHCO) + C(SCH2O) &
	    + C(SCH3) + 3.6 * C(SCO2) &
	    + 2 * C(SCH4) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SC2H) + C(SHCCO) &
	    + C(SC2H2) + C(SC3H3) &
	    + C(SAXC3H5) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SA1XC6H6)
      M(MM3) = M(MM3) + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H8) &
	    + C(SC5H6) + C(SA2XC10H8) &
	    + C(SC5H10) + C(SC5H11) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1CHOXC7H6O) + C(SA1CH3XC7H8) &
	    + C(SC7H15) + C(SNXC7H16) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6)
      M(MM2) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 2 * C(SH2) + C(SH) &
	    + C(SOH) + 6.3 * C(SH2O) &
	    + C(SO2) + C(SHO2) &
	    + C(SCH) + 1.75 * C(SCO) &
	    + C(SHCO) + C(SCH2O) &
	    + C(SCH3) + 3.6 * C(SCO2) &
	    + 2 * C(SCH4) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SC2H) + C(SHCCO) &
	    + C(SC2H2) + C(SC3H3) &
	    + C(SAXC3H5) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SA1XC6H6)
      M(MM2) = M(MM2) + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H8) &
	    + C(SC5H6) + C(SA2XC10H8) &
	    + C(SC5H10) + C(SC5H11) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1CHOXC7H6O) + C(SA1CH3XC7H8) &
	    + C(SC7H15) + C(SNXC7H16) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6)
      M(MM5) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 0.75 * C(SH2) + C(SH) &
	    + C(SOH) + 11.89 * C(SH2O) &
	    + 0.85 * C(SO2) + C(SHO2) &
	    + C(SCH) + 1.09 * C(SCO) &
	    + C(SHCO) + C(SCH2O) &
	    + C(SCH3) + 2.18 * C(SCO2) &
	    + C(SCH4) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SC2H) + C(SHCCO) &
	    + C(SC2H2) + C(SC3H3) &
	    + C(SAXC3H5) + C(SNXC3H7) &
	    + C(SC2H6) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SA1XC6H6)
      M(MM5) = M(MM5) + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H8) &
	    + C(SC5H6) + C(SA2XC10H8) &
	    + C(SC5H10) + C(SC5H11) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1CHOXC7H6O) + C(SA1CH3XC7H8) &
	    + C(SC7H15) + C(SNXC7H16) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6)
      M(MM9) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 2 * C(SH2) + C(SH) &
	    + C(SOH) + 6 * C(SH2O) &
	    + C(SO2) + C(SHO2) &
	    + C(SCH) + 1.5 * C(SCO) &
	    + C(SHCO) + C(SCH2O) &
	    + C(SCH3) + 2 * C(SCO2) &
	    + 3 * C(SCH4) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SC2H) + C(SHCCO) &
	    + C(SC2H2) + C(SC3H3) &
	    + C(SAXC3H5) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SA1XC6H6)
      M(MM9) = M(MM9) + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H8) &
	    + C(SC5H6) + C(SA2XC10H8) &
	    + C(SC5H10) + C(SC5H11) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1CHOXC7H6O) + C(SA1CH3XC7H8) &
	    + C(SC7H15) + C(SNXC7H16) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6)
      M(MM7) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 2 * C(SH2) + C(SH) &
	    + C(SOH) + 12 * C(SH2O) &
	    + C(SO2) + C(SHO2) &
	    + C(SCH) + 1.75 * C(SCO) &
	    + C(SHCO) + C(SCH2O) &
	    + C(SCH3) + 3.6 * C(SCO2) &
	    + 2 * C(SCH4) + C(SC2H3) &
	    + C(SC2H4) + C(SC2H5) &
	    + C(SC2H) + C(SHCCO) &
	    + C(SC2H2) + C(SC3H3) &
	    + C(SAXC3H5) + C(SNXC3H7) &
	    + 3 * C(SC2H6) + C(SPXC3H4) &
	    + C(SAXC3H4) + C(SA1XC6H6)
      M(MM7) = M(MM7) + C(SA1XXC6H5) + C(SC5H5) &
	    + C(SC3H6) + C(SC4H8) &
	    + C(SC5H6) + C(SA2XC10H8) &
	    + C(SC5H10) + C(SC5H11) &
	    + C(SA1C2H2XC8H7) + C(SA1CH2XC7H7) &
	    + C(SA1CHOXC7H6O) + C(SA1CH3XC7H8) &
	    + C(SC7H15) + C(SNXC7H16) &
	    + C(SA1C2HYXC8H5) + C(SA2XXC10H7) &
	    + C(SA1C2HXC8H6)
      M(MM8) = C(SN2) + C(SSXCH2) &
	    + C(STXCH2) + C(SO) &
	    + 2 * C(SH2) + C(SH) &
	    + C(SOH) + C(SO2) &
	    + C(SHO2) + C(SCH) &
	    + 1.75 * C(SCO) + C(SHCO) &
	    + C(SCH2O) + C(SCH3) &
	    + 3.6 * C(SCO2) + 2 * C(SCH4) &
	    + C(SC2H3) + C(SC2H4) &
	    + C(SC2H5) + C(SC2H) &
	    + C(SHCCO) + C(SC2H2) &
	    + C(SC3H3) + C(SAXC3H5) &
	    + C(SNXC3H7) + 3 * C(SC2H6) &
	    + C(SPXC3H4) + C(SAXC3H4) &
	    + C(SA1XC6H6) + C(SA1XXC6H5)
      M(MM8) = M(MM8) + C(SC5H5) + C(SC3H6) &
	    + C(SC4H8) + C(SC5H6) &
	    + C(SA2XC10H8) + C(SC5H10) &
	    + C(SC5H11) + C(SA1C2H2XC8H7) &
	    + C(SA1CH2XC7H7) + C(SA1CHOXC7H6O) &
	    + C(SA1CH3XC7H8) + C(SC7H15) &
	    + C(SNXC7H16) + C(SA1C2HYXC8H5) &
	    + C(SA2XXC10H7) + C(SA1C2HXC8H6)


      K(RG28F) = 1.5000000000D+10 * exp(-2510000 / RT)
      K(RG28B) = 1.0936408744D+10 &
	   * exp(-0.0622863 * LT - 39812987.98 / RT)
      K(R2F) = 4.5900000000D+01 &
	   * exp(2.7 * LT - 26190000 / RT)
      K(R2B) = 2.7517349858D+01 &
	   * exp(2.66385 * LT - 20100582.44 / RT)
      K(R10) = 9.4300000000D+12 * exp(-1 * LT)
      K(R9) = 4.4000000000D+16 * exp(-2 * LT)
      K(R3F) = 1.7300000000D+05 &
	   * exp(1.51 * LT - 14350000 / RT)
      K(R3B) = 1.9032628245D+06 &
	   * exp(1.40625 * LT - 77394524.74 / RT)
      K(R4F) = 3.9700000000D+01 &
	   * exp(2.4 * LT + 8830000 / RT)
      K(R4B) = 7.2853303338D+02 &
	   * exp(2.33241 * LT - 60303942.3 / RT)
      K0TROE = 6.3300000000D+13 * exp(-1.4 * LT)
      KINFTROE = 5.1200000000D+09 * exp(0.44 * LT)
      FCTROE = 0.5 * EXP( -TEMP / 1e-10 ) &
	   + 0.5 * EXP( -TEMP / 1e+10 ) &
	   + 1 * EXP( -1e+10 / TEMP )
      K(R12) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K(R1F) = 2.6400000000D+13 &
	   * exp(-0.67 * LT - 71300000 / RT)
      K(R1B) = 5.2766713024D+10 &
	   * exp(-0.234291 * LT - 479604.2721 / RT)
      K(R13) = 3.6490000000D+03 &
	   * exp(2.07 * LT + 4570000 / RT)
      K(R15) = 3.9700000000D+09 * exp(-2810000 / RT)
      K(R16) = 7.4900000000D+10 * exp(-2660000 / RT)
      K(R17) = 4.0000000000D+10
      K(R18) = 2.3800000000D+10 * exp(2090000 / RT)
      K(R19) = 1.0000000000D+13 * exp(-72510000 / RT)
      K(RG04) = 5.7000000000D+10
      K(RG05) = 3.0000000000D+10
      K(RG06F) = 1.0800000000D+11 * exp(-13010000 / RT)
      K(RG06B) = 1.7581777781D+12 &
	   * exp(-0.292468 * LT - 1691026.963 / RT)
      K(RG08) = 5.7100000000D+09 * exp(3160000 / RT)
      K(RG09) = 6.7100000000D+10
      K(RG15) = 8.0000000000D+10
      K(RG16) = 2.0000000000D+10
      K(RG17F) = 1.1300000000D+04 &
	   * exp(2 * LT - 12550000 / RT)
      K(RG17B) = 7.6364598733D+03 &
	   * exp(2.18872 * LT - 86913497.78 / RT)
      K(RG18F) = 5.0000000000D+02 &
	   * exp(2 * LT - 30250000 / RT)
      K(RG18B) = 2.0974636059D+05 &
	   * exp(1.4875 * LT - 61407685.18 / RT)
      K(RG19) = 5.8000000000D+09 * exp(-6280000 / RT)
      K(RG20) = 2.4000000000D+09 * exp(-6280000 / RT)
      K(RG21) = 5.0000000000D+09 * exp(-6280000 / RT)
      K(RG34F) = 7.0000000000D+10
      K(RG34B) = 2.1409471364D+13 &
	   * exp(-0.574789 * LT - 68460673.16 / RT)
      K(RG35) = 2.8000000000D+10
      K(RG36) = 1.2000000000D+10
      K0TROE = 1.8800000000D+32 &
	   * exp(-6.36 * LT - 21090000 / RT)
      KINFTROE = 4.8200000000D+14 &
	   * exp(-1.16 * LT - 4790000 / RT)
      FCTROE = 0.3973 * EXP( -TEMP / 208 ) &
	   + 0.6027 * EXP( -TEMP / 3922 ) &
	   + 1 * EXP( -10180 / TEMP )
      K(RG37) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG38F) = 3.0000000000D+10
      K(RG38B) = 2.1872817488D+10 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG39) = 6.8200000000D+07 &
	   * exp(0.25 * LT + 3910000 / RT)
      K0TROE = 3.4700000000D+32 &
	   * exp(-6.3 * LT - 21230000 / RT)
      KINFTROE = 6.9200000000D+10 * exp(0.18 * LT)
      FCTROE = 0.217 * EXP( -TEMP / 74 ) &
	   + 0.873 * EXP( -TEMP / 2941 ) &
	   + 1 * EXP( -6964 / TEMP )
      K(RG51F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K0TROE = 1.5077439481D+38 &
	   * exp(-6.46997 * LT - 463233395.7 / RT)
      KINFTROE = 3.0067977293D+16 &
	   * exp(0.0100263 * LT - 442003395.7 / RT)
      K(RG51B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM9) )
      K(RG52) = 5.0600000000D+10
      K(RG53) = 3.3700000000D+10
      K0TROE = 4.0000000000D+30 &
	   * exp(-5.92 * LT - 13140000 / RT)
      KINFTROE = 2.7900000000D+15 &
	   * exp(-1.43 * LT - 5570000 / RT)
      FCTROE = 0.588 * EXP( -TEMP / 195 ) &
	   + 0.412 * EXP( -TEMP / 5900 ) &
	   + 1 * EXP( -6394 / TEMP )
      K(RG54) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG55F) = 5.6000000000D+04 &
	   * exp(1.6 * LT - 22680000 / RT)
      K(RG55B) = 1.4686430659D+03 &
	   * exp(2.00876 * LT - 54566839.56 / RT)
      K(RG57F) = 6.4400000000D+14 &
	   * exp(-1.34 * LT - 5930000 / RT)
      K(RG57B) = 2.3164910420D+13 &
	   * exp(-0.868957 * LT - 513851.5818 / RT)
      K(RG58) = 1.3800000000D+10 * exp(-127530000 / RT)
      K(RG59) = 5.8700000000D+08 * exp(-57910000 / RT)
      K(RG65) = 1.0000000000D+10
      K(RG66) = 3.6100000000D+09
      K(RG69) = 3.0000000000D+10
      K(RG72) = 1.0000000000D+11
      K(RG74F) = 6.8400000000D+09 &
	   * exp(0.1 * LT - 44350000 / RT)
      K(RG74B) = 1.1614797579D+15 &
	   * exp(-1.0902 * LT - 7197984.006 / RT)
      K(RG95) = 2.6320000000D+09 &
	   * exp(-0.06 * LT - 57160000 / RT)
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
      K(RG93) = 6.0000000000D+10
      K(RG94F) = 2.4600000000D+03 &
	   * exp(2 * LT - 34600000 / RT)
      K(RG94B) = 5.5511242907D+02 &
	   * exp(2.00771 * LT - 56845777.28 / RT)
      K0TROE = 1.1700000000D+18 &
	   * exp(-2.79 * LT - 17540000 / RT)
      KINFTROE = 1.3600000000D+07 * exp(-9980000 / RT)
      FCTROE = 1 * EXP( -0 / TEMP )
      K(R27) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM7) )
      K(R28F) = 8.0000000000D+08 &
	   * exp(0.14 * LT - 30760000 / RT)
      K(R28B) = 4.3798719198D+15 &
	   * exp(-1.13743 * LT - 138461819.8 / RT)
      K(R29F) = 8.7800000000D+07 &
	   * exp(0.03 * LT + 70000 / RT)
      K(R29B) = 4.8069094320D+14 &
	   * exp(-1.24743 * LT - 107631819.8 / RT)
      K(R31) = 3.0100000000D+10 * exp(-96230000 / RT)
      K(RG40F) = 9.0000000000D+09
      K(RG40B) = 6.5618452463D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(R32) = 1.2000000000D+11
      K(R33) = 3.0000000000D+10
      K(R34) = 3.0000000000D+10
      K(R35) = 3.0200000000D+10
      K(R36F) = 1.8700000000D+14 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R36B) = 3.5306519794D+10 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(R37F) = 2.2400000000D+15 &
	   * exp(-1 * LT - 71130000 / RT)
      K(R37B) = 4.2292301785D+11 &
	   * exp(-0.767288 * LT - 5356428.13 / RT)
      K(R38) = 1.2000000000D+07 &
	   * exp(0.81 * LT + 3040000 / RT)
      K(RG70) = 2.6500000000D+10
      K(RG45F) = 5.7400000000D+04 &
	   * exp(1.9 * LT - 11470000 / RT)
      K(RG45B) = 9.5590953891D+01 &
	   * exp(2.34945 * LT - 74956081.25 / RT)
      K(RG46) = 3.9000000000D+10 * exp(-14810000 / RT)
      K(RG47) = 3.4300000000D+06 &
	   * exp(1.18 * LT + 1870000 / RT)
      K(RG71) = 3.3200000000D+00 &
	   * exp(2.81 * LT - 24520000 / RT)
      K(RG11) = 1.9000000000D+11 * exp(-66070000 / RT)
      K(RG41F) = 7.0000000000D+09
      K(RG41B) = 5.1036574138D+09 &
	   * exp(-0.0622863 * LT - 37302987.98 / RT)
      K(RG42) = 1.4000000000D+10
      K(RG105) = 5.0000000000D+10
      K(RG106) = 2.0000000000D+10
      K(RG107) = 1.0000000000D+10 * exp(3160000 / RT)
      K(RG108F) = 3.3100000000D+03 &
	   * exp(2.26 * LT - 3770000 / RT)
      K(RG108B) = 9.6742642417D+06 &
	   * exp(1.66346 * LT - 127161378.7 / RT)
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
      K(RG116) = 8.1000000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RG117) = 1.2500000000D+04 &
	   * exp(2 * LT - 7950000 / RT)
      K(RG118) = 3.3240000000D+13 &
	   * exp(-0.44 * LT - 128440000 / RT)
      K(RG119) = 2.6300000000D+03 &
	   * exp(2.14 * LT - 71380000 / RT)
      K(RG121) = 7.5300000000D+03 &
	   * exp(1.55 * LT - 8810000 / RT)
      K(RG122) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR004F) = 1.9000000000D+11
      K(RR004B) = 2.5098929533D+16 &
	   * exp(-1.35613 * LT - 105055737.5 / RT)
      K0TROE = 1.4000000000D+24 &
	   * exp(-3.86 * LT - 13890000 / RT)
      KINFTROE = 6.0800000000D+09 &
	   * exp(0.27 * LT - 1170000 / RT)
      FCTROE = 0.218 * EXP( -TEMP / 207.5 ) &
	   + 0.782 * EXP( -TEMP / 2663 ) &
	   + 1 * EXP( -6095 / TEMP )
      K(RG129) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG130) = 3.0000000000D+10
      K(RG131) = 1.0300000000D+10 &
	   * exp(0.21 * LT + 1790000 / RT)
      K(RG132) = 5.0000000000D+09
      K(RG133) = 1.3400000000D+03 &
	   * exp(1.61 * LT + 1610000 / RT)
      K(RG134) = 3.0300000000D+08 &
	   * exp(0.29 * LT - 50000 / RT)
      K(RG135) = 4.5800000000D+13 &
	   * exp(-1.39 * LT - 4250000 / RT)
      K(RR011) = 9.0000000000D+10
      K(RR013) = 9.0300000000D+09 * exp(3200000 / RT)
      K(RR015F) = 1.9300000000D+15 &
	   * exp(-1.25 * LT - 32090000 / RT)
      K(RR015B) = 1.0584325883D+24 &
	   * exp(-3.18571 * LT - 101188124 / RT)
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
      K(RG156) = 1.3300000000D+03 &
	   * exp(2.53 * LT - 51210000 / RT)
      K(RG157) = 7.6600000000D+06 &
	   * exp(0.88 * LT - 4770000 / RT)
      K(RG158) = 7.1500000000D+01 &
	   * exp(2.47 * LT - 3890000 / RT)
      K(RG159) = 3.8900000000D+05 &
	   * exp(1.36 * LT - 3710000 / RT)
      K(RG160) = 1.3100000000D-04 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RG161) = 3.7500000000D+33 &
	   * exp(-7.8 * LT - 29540000 / RT)
      K(RG162) = 2.2700000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
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
      K(RG167) = 3.1700000000D+10 &
	   * exp(0.03 * LT + 1650000 / RT)
      K(RG171) = 1.9200000000D+04 &
	   * exp(1.02 * LT + 8510000 / RT)
      K(RR034) = 1.2000000000D+11
      K(RR035) = 3.0000000000D+08
      K(RR037) = 3.1000000000D+10
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
      K(RG174) = 1.7000000000D+02 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RG175) = 8.9800000000D+04 &
	   * exp(1.92 * LT - 23810000 / RT)
      K(RG176) = 1.6100000000D+03 &
	   * exp(2.22 * LT - 3100000 / RT)
      K(RG178) = 8.4300000000D+11 * exp(-93120000 / RT)
      K0TROE = 1.2130000000D+32 &
	   * exp(-5.18 * LT - 320780000 / RT)
      KINFTROE = 2.2550000000D+20 &
	   * exp(-1.44 * LT - 312680000 / RT)
      FCTROE = 0.4243 * EXP( -TEMP / 237 ) &
	   + 0.5757 * EXP( -TEMP / 1652 ) &
	   + 1 * EXP( -5069 / TEMP )
      K(RG10) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RG109) = 1.0000000000D+11
      K(RG110) = 1.0000000000D+11
      K(RG111) = 4.2000000000D+07 * exp(-3570000 / RT)
      K(RG114) = 1.0000000000D+10
      K(RR009) = 1.0000000000D+08 * exp(-12550000 / RT)
      K(RR022) = 5.0000000000D+10
      K(RR024) = 1.0000000000D+10
      K(RR055F) = 7.9400000000D+26 &
	   * exp(-5.06 * LT - 20340000 / RT)
      K(RR055B) = 8.1610655573D+31 &
	   * exp(-5.06898 * LT - 393338860.7 / RT)
      K(RR056F) = 3.1600000000D+26 &
	   * exp(-5 * LT - 19710000 / RT)
      K(RR056B) = 1.2666373930D+32 &
	   * exp(-5.11373 * LT - 387807163.3 / RT)
      K(RR058) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR061) = 6.9500000000D+10
      K(RR062) = 1.7000000000D+02 &
	   * exp(1.7 * LT - 6280000 / RT)
      K(RR063) = 8.0000000000D+08
      K(RLP01F) = 1.8700000000D+43 &
	   * exp(-9.84 * LT - 70310000 / RT)
      K(RLP01B) = 2.6958773200D+59 &
	   * exp(-11.2977 * LT - 676573672.1 / RT)
      K(RP011F) = 5.7700000000D+34 &
	   * exp(-7 * LT - 131820000 / RT)
      K(RP011B) = 5.1317962822D+44 &
	   * exp(-8.24617 * LT - 262035013.4 / RT)
      K(RCP11F) = 2.3500000000D+08 * exp(-41820000 / RT)
      K(RCP11B) = 1.8847202119D+20 &
	   * exp(-1.07148 * LT - 353385243.7 / RT)
      K(RR005F) = 3.4600000000D+09 &
	   * exp(0.44 * LT - 22860000 / RT)
      K(RR005B) = 1.8216946136D+04 &
	   * exp(1.58056 * LT - 46357562.82 / RT)
      K(RR078F) = 8.5000000000D+01 &
	   * exp(2.7 * LT - 24020000 / RT)
      K(RR078B) = 1.9329120587D-01 &
	   * exp(3.05922 * LT - 84112627.12 / RT)
      K(RR080) = 8.0500000000D+02 &
	   * exp(2.22 * LT - 3100000 / RT)
      K(RR081) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RR088) = 4.0500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR089) = 6.2500000000D+03 &
	   * exp(2 * LT - 7950000 / RT)
      K(RR094) = 1.2800000000D+06 &
	   * exp(0.73 * LT - 10790000 / RT)
      K(RR006F) = 8.9500000000D+10 &
	   * exp(-0.02 * LT - 47070000 / RT)
      K(RR006B) = 1.2083246981D+05 &
	   * exp(1.2253 * LT - 75469260.17 / RT)
      K(RR071F) = 7.7600000000D+39 &
	   * exp(-7.8 * LT - 328220000 / RT)
      K(RR071B) = 1.9898615631D+39 &
	   * exp(-7.69525 * LT - 333121697.4 / RT)
      K(RR073F) = 2.4700000000D+12 &
	   * exp(-0.33 * LT - 26930000 / RT)
      K(RR073B) = 6.3337088411D+11 &
	   * exp(-0.225252 * LT - 31831697.35 / RT)
      K(RR083) = 1.3300000000D+03 &
	   * exp(2.53 * LT - 51210000 / RT)
      K(RR084) = 1.3100000000D-04 &
	   * exp(4.2 * LT + 3600000 / RT)
      K(RR085) = 2.2700000000D+02 &
	   * exp(2 * LT - 38490000 / RT)
      K(RR074F) = 2.0100000000D+46 &
	   * exp(-10.77 * LT - 82100000 / RT)
      K(RR074B) = 3.2976487096D+51 &
	   * exp(-11.2885 * LT - 326398622 / RT)
      K(RR103) = 9.5600000000D+00 &
	   * exp(2.8 * LT - 13770000 / RT)
      K(RR104) = 6.0300000000D+09
      K(RR105) = 4.8600000000D+08 &
	   * exp(-0.32 * LT + 550000 / RT)
      K(RR111) = 9.7100000000D+17 &
	   * exp(-2.7 * LT - 104520000 / RT)
      K(RP008) = 2.1600000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K0TROE = 6.2600000000D+32 &
	   * exp(-6.66 * LT - 29290000 / RT)
      KINFTROE = 3.0600000000D+11 &
	   * exp(-0.37 * LT - 16870000 / RT)
      FCTROE = 1 * EXP( -TEMP / 1310 ) &
	   + 1 * EXP( -48100 / TEMP )
      K(RG183F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 2.4335109175D+36 &
	   * exp(-6.79608 * LT - 166021318.4 / RT)
      KINFTROE = 1.1895436753D+15 &
	   * exp(-0.506081 * LT - 153601318.4 / RT)
      K(RG183B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(RR014F) = 4.0400000000D+42 &
	   * exp(-7.67 * LT - 467900000 / RT)
      K(RR014B) = 3.9209380170D+30 &
	   * exp(-6.00376 * LT - 34921721.89 / RT)
      K(RR016F) = 5.9300000000D+51 &
	   * exp(-11.76 * LT - 98530000 / RT)
      K(RR016B) = 1.1141410801D+55 &
	   * exp(-11.4905 * LT - 462410154.1 / RT)
      K(RR141F) = 8.0000000000D+18 &
	   * exp(-2.39 * LT - 46780000 / RT)
      K(RR141B) = 4.3155004306D+12 &
	   * exp(-0.859397 * LT - 80320978.73 / RT)
      K(RR142) = 1.2000000000D+05 &
	   * exp(1.6 * LT - 1370000 / RT)
      K(RR143) = 3.5000000000D+04 &
	   * exp(1.6 * LT + 4070000 / RT)
      K(RR144) = 6.6000000000D+02 &
	   * exp(2.54 * LT - 28270000 / RT)
      K(RR145) = 9.6500000000D+01 &
	   * exp(2.68 * LT - 15550000 / RT)
      K(RR146) = 2.0000000000D+05 &
	   * exp(1.46 * LT - 2250000 / RT)
      K(RR148) = 4.5200000000D-04 &
	   * exp(3.65 * LT - 29930000 / RT)
      K(RHP40) = 7.6800000000D+09 &
	   * exp(0.11 * LT - 6190000 / RT)
      K(RHP43) = 5.0000000000D+15 * exp(-297060000 / RT)
      K(RCP06) = 2.5610000000D+09 &
	   * exp(0.06 * LT - 13040000 / RT)
      K(RCP12) = 6.3900000000D+26 &
	   * exp(-4.03 * LT - 147300000 / RT)
      K(RCP13) = 7.0000000000D+10
      K(RCP16) = 6.7760000000D+26 &
	   * exp(-4.7 * LT - 48780000 / RT)
      K(RCP18) = 1.1900000000D+30 &
	   * exp(-6.52 * LT - 56070000 / RT)
      K(RCP20) = 3.0200000000D+10
      K(RCP01F) = 1.7300000000D+68 &
	   * exp(-15.16 * LT - 486900000 / RT)
      K(RCP01B) = 7.5313864907D+61 &
	   * exp(-14.464 * LT - 137417681 / RT)
      K(RCP02) = 2.8000000000D+10 * exp(-9450000 / RT)
      K(RCP03) = 6.6000000000D+11 * exp(-51650000 / RT)
      K(RCP04) = 4.7700000000D+01 &
	   * exp(2.71 * LT - 4630000 / RT)
      K(RCP05) = 3.0800000000D+03 * exp(2 * LT)
      K(RCP08) = 1.8000000000D-04 * exp(4 * LT)
      K(RHP31) = 7.1000000000D+09 &
	   * exp(0.12 * LT - 6110000 / RT)
      K(RHP33) = 9.1700000000D+20 &
	   * exp(-1.63 * LT - 309570000 / RT)
      K(RHP29) = 7.4600000000D+21 &
	   * exp(-2.61 * LT - 134000000 / RT)
      K(RHP30) = 8.4600000000D+14 &
	   * exp(-0.47 * LT - 157390000 / RT)
      K(RHP32) = 3.1500000000D-19 &
	   * exp(8.84 * LT - 29730000 / RT)
      K(ROX12) = 2.6000000000D+10 * exp(-25610000 / RT)
      K(ROX13) = 3.0000000000D+10 * exp(-37570000 / RT)
      K(ROX14) = 1.0000000000D+11
      K(ROX15) = 3.0000000000D+10
      K(ROX16) = 3.0000000000D+10
      K(RK012F) = 3.2900000000D+03 &
	   * exp(2.05 * LT - 13230000 / RT)
      K(RK012B) = 2.6128801700D+14 &
	   * exp(0.279954 * LT - 191326710.8 / RT)
      K(ROX00F) = 1.2900000000D+61 &
	   * exp(-12.48 * LT - 619590000 / RT)
      K(ROX00B) = 7.9583820932D+54 &
	   * exp(-12.2685 * LT - 143541341.3 / RT)
      K(ROX05) = 6.0200000000D+05 &
	   * exp(1.8 * LT - 68420000 / RT)
      K(ROX06) = 4.0300000000D-01 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(ROX17) = 2.7520000000D-05 &
	   * exp(4.46 * LT - 57060000 / RT)
      K(ROX08) = 2.2200000000D+10 * exp(-18960000 / RT)
      K(ROX11) = 1.3200000000D-01 &
	   * exp(3.25 * LT - 23390000 / RT)
      K(RT05F) = 8.2000000000D+14 * exp(-337550000 / RT)
      K(RT05B) = 4.3614434871D+01 &
	   * exp(2.2454 * LT - 46091406.02 / RT)
      K(RT04F) = 5.8300000000D+64 &
	   * exp(-14.15 * LT - 285890000 / RT)
      K(RT04B) = 8.6112522502D+55 &
	   * exp(-12.1057 * LT - 218394895.1 / RT)
      K(RT16) = 3.3100000000D+11
      K(RT18) = 1.0600000000D+13 &
	   * exp(-0.94 * LT - 10560000 / RT)
      K(RT19) = 4.3200000000D+36 &
	   * exp(-7.74 * LT - 99800000 / RT)
      K(RT01F) = 2.3100000000D+03 &
	   * exp(2.17 * LT - 17420000 / RT)
      K(RT01B) = 4.5426319227D-03 &
	   * exp(3.57723 * LT - 50953235.53 / RT)
      K(RT02F) = 1.2500000000D+18 &
	   * exp(-0.6 * LT - 396590000 / RT)
      K(RT02B) = 1.0266984379D+15 &
	   * exp(-1.02554 * LT - 21569681.73 / RT)
      K(RT03F) = 2.1600000000D+29 &
	   * exp(-3.58 * LT - 460930000 / RT)
      K(RT03B) = 2.6205013989D+17 &
	   * exp(-1.96127 * LT - 18414576.83 / RT)
      K(RT07) = 6.4700000000D-03 &
	   * exp(3.98 * LT - 14160000 / RT)
      K(RT08) = 1.6200000000D+10 * exp(-11590000 / RT)
      K(RT13) = 4.2200000000D+11 * exp(-93120000 / RT)
      K(RHP10) = 1.8900000000D+12 &
	   * exp(0.02 * LT - 116250000 / RT)
      K(RHP11) = 7.7300000000D+18 &
	   * exp(-1.75 * LT - 133780000 / RT)
      K(RHP12) = 2.5300000000D+18 &
	   * exp(-1.65 * LT - 132560000 / RT)
      K(RHP13) = 2.4900000000D+16 &
	   * exp(-1.18 * LT - 123500000 / RT)
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
      K(RHP07) = 1.4600000000D+01 &
	   * exp(2.57 * LT - 29010000 / RT)
      K(RT33) = 2.6100000000D+15 &
	   * exp(0.15 * LT - 337020000 / RT)
      K(RT34) = 2.0500000000D+06 &
	   * exp(1.16 * LT - 10060000 / RT)
      K(RT38) = 2.7200000000D+03 &
	   * exp(1.77 * LT - 24770000 / RT)
      K(RK100F) = 1.3400000000D+01 &
	   * exp(2.5 * LT - 5370000 / RT)
      K(RK100B) = 3.3117446183D+14 &
	   * exp(0.979368 * LT - 386582470.9 / RT)
      K(RP106F) = 3.6200000000D+25 &
	   * exp(-4.24 * LT - 99850000 / RT)
      K(RP106B) = 9.5272754603D+36 &
	   * exp(-6.03013 * LT - 381352625.2 / RT)
      K(RP016F) = 2.1000000000D+60 &
	   * exp(-12.4 * LT - 619550000 / RT)
      K(RP016B) = 3.1383628562D+54 &
	   * exp(-12.1693 * LT - 132872595.6 / RT)
      K(RK017F) = 1.3200000000D+05 &
	   * exp(1.88 * LT - 70380000 / RT)
      K(RK017B) = 4.6108110511D+01 &
	   * exp(2.46097 * LT - 16794083.4 / RT)
      K(RP018F) = 1.3400000000D-01 &
	   * exp(3.33 * LT - 6090000 / RT)
      K(RP018B) = 5.1494501025D-04 &
	   * exp(3.80722 * LT - 15548608.14 / RT)
      K(RLP06F) = 1.3400000000D+17 &
	   * exp(-0.86 * LT - 172540000 / RT)
      K(RLP06B) = 8.3211163704D+12 &
	   * exp(-0.568959 * LT - 41173891.76 / RT)
      K(RLP07F) = 7.0910000000D+10 &
	   * exp(-0.26 * LT - 29300000 / RT)
      K(RLP07B) = 2.1330580802D+20 &
	   * exp(-1.4921 * LT - 306079607.2 / RT)
      K(RP108F) = 8.6000000000D+60 &
	   * exp(-12.48 * LT - 619550000 / RT)
      K(RP108B) = 6.5571626174D+54 &
	   * exp(-12.2468 * LT - 105939351.1 / RT)
      K(RK109F) = 2.6500000000D+05 &
	   * exp(1.87 * LT - 71530000 / RT)
      K(RK109B) = 4.7226190083D+01 &
	   * exp(2.45348 * LT + 8989161.129 / RT)
      K(RK110F) = 9.6300000000D-01 &
	   * exp(3.02 * LT - 18300000 / RT)
      K(RK110B) = 1.8880608443D-03 &
	   * exp(3.49973 * LT - 825363.6114 / RT)

      W(RG28F) = K(RG28F) * C(SSXCH2) * C(SN2)
      W(RG28B) = K(RG28B) * C(SN2) * C(STXCH2)
      W(R2F) = K(R2F) * C(SO) * C(SH2)
      W(R2B) = K(R2B) * C(SOH) * C(SH)
      W(R10) = K(R10) * C(SH) * C(SO) * M(MM3)
      W(R9) = K(R9) * C(SH) * C(SOH) * M(MM2)
      W(R3F) = K(R3F) * C(SOH) * C(SH2)
      W(R3B) = K(R3B) * C(SH2O) * C(SH)
      W(R4F) = K(R4F) * C(SOH) * C(SOH)
      W(R4B) = K(R4B) * C(SH2O) * C(SO)
      W(R12) = K(R12) * C(SH) * C(SO2)
      W(R1F) = K(R1F) * C(SH) * C(SO2)
      W(R1B) = K(R1B) * C(SOH) * C(SO)
      W(R13) = K(R13) * C(SHO2) * C(SH)
      W(R15) = K(R15) * C(SHO2) * C(SH)
      W(R16) = K(R16) * C(SHO2) * C(SH)
      W(R17) = K(R17) * C(SHO2) * C(SO)
      W(R18) = K(R18) * C(SHO2) * C(SOH)
      W(R19) = K(R19) * C(SHO2) * C(SOH)
      W(RG04) = K(RG04) * C(SCH) * C(SO)
      W(RG05) = K(RG05) * C(SCH) * C(SOH)
      W(RG06F) = K(RG06F) * C(SCH) * C(SH2)
      W(RG06B) = K(RG06B) * C(SH) * C(STXCH2)
      W(RG08) = K(RG08) * C(SCH) * C(SH2O)
      W(RG09) = K(RG09) * C(SCH) * C(SO2)
      W(RG15) = K(RG15) * C(STXCH2) * C(SO)
      W(RG16) = K(RG16) * C(STXCH2) * C(SOH)
      W(RG17F) = K(RG17F) * C(STXCH2) * C(SOH)
      W(RG17B) = K(RG17B) * C(SH2O) * C(SCH)
      W(RG18F) = K(RG18F) * C(STXCH2) * C(SH2)
      W(RG18B) = K(RG18B) * C(SCH3) * C(SH)
      W(RG19) = K(RG19) * C(STXCH2) * C(SO2)
      W(RG20) = K(RG20) * C(STXCH2) * C(SO2)
      W(RG21) = K(RG21) * C(STXCH2) * C(SO2)
      W(RG34F) = K(RG34F) * C(SSXCH2) * C(SH2)
      W(RG34B) = K(RG34B) * C(SH) * C(SCH3)
      W(RG35) = K(RG35) * C(SSXCH2) * C(SO2)
      W(RG36) = K(RG36) * C(SSXCH2) * C(SO2)
      W(RG37) = K(RG37) * C(SSXCH2) * C(SH2O)
      W(RG38F) = K(RG38F) * C(SSXCH2) * C(SH2O)
      W(RG38B) = K(RG38B) * C(SH2O) * C(STXCH2)
      W(RG39) = K(RG39) * C(SSXCH2) * C(SH2O)
      W(RG51F) = K(RG51F) * C(SCH3) * C(SH)
      W(RG51B) = K(RG51B) * C(SCH4)
      W(RG52) = K(RG52) * C(SCH3) * C(SO)
      W(RG53) = K(RG53) * C(SCH3) * C(SO)
      W(RG54) = K(RG54) * C(SCH3) * C(SOH)
      W(RG55F) = K(RG55F) * C(SCH3) * C(SOH)
      W(RG55B) = K(RG55B) * C(SH2O) * C(STXCH2)
      W(RG57F) = K(RG57F) * C(SCH3) * C(SOH)
      W(RG57B) = K(RG57B) * C(SH2O) * C(SSXCH2)
      W(RG58) = K(RG58) * C(SCH3) * C(SO2)
      W(RG59) = K(RG59) * C(SCH3) * C(SO2)
      W(RG65) = K(RG65) * C(SCH3) * C(SHO2)
      W(RG66) = K(RG66) * C(SCH3) * C(SHO2)
      W(RG69) = K(RG69) * C(SCH3) * C(SCH)
      W(RG72) = K(RG72) * C(SCH3) * C(STXCH2)
      W(RG74F) = K(RG74F) * C(SCH3) * C(SCH3)
      W(RG74B) = K(RG74B) * C(SH) * C(SC2H5)
      W(RG95) = K(RG95) * C(SCH3) * C(SCH3)
      W(RG90F) = K(RG90F) * C(SCH4) * C(SH)
      W(RG90B) = K(RG90B) * C(SH2) * C(SCH3)
      W(RG91F) = K(RG91F) * C(SCH4) * C(SO)
      W(RG91B) = K(RG91B) * C(SOH) * C(SCH3)
      W(RG92F) = K(RG92F) * C(SCH4) * C(SOH)
      W(RG92B) = K(RG92B) * C(SH2O) * C(SCH3)
      W(RG93) = K(RG93) * C(SCH4) * C(SCH)
      W(RG94F) = K(RG94F) * C(SCH4) * C(STXCH2)
      W(RG94B) = K(RG94B) * C(SCH3) * C(SCH3)
      W(R27) = K(R27) * C(SCO) * C(SO)
      W(R28F) = K(R28F) * C(SCO) * C(SOH)
      W(R28B) = K(R28B) * C(SH) * C(SCO2)
      W(R29F) = K(R29F) * C(SCO) * C(SOH)
      W(R29B) = K(R29B) * C(SH) * C(SCO2)
      W(R31) = K(R31) * C(SCO) * C(SHO2)
      W(RG40F) = K(RG40F) * C(SCO) * C(SSXCH2)
      W(RG40B) = K(RG40B) * C(STXCH2) * C(SCO)
      W(R32) = K(R32) * C(SHCO) * C(SH)
      W(R33) = K(R33) * C(SHCO) * C(SO)
      W(R34) = K(R34) * C(SHCO) * C(SO)
      W(R35) = K(R35) * C(SHCO) * C(SOH)
      W(R36F) = K(R36F) * C(SHCO) * M(MM8)
      W(R36B) = K(R36B) * C(SH) * C(SCO) * M(MM8)
      W(R37F) = K(R37F) * C(SHCO) * C(SH2O)
      W(R37B) = K(R37B) * C(SH2O) * C(SH) * C(SCO)
      W(R38) = K(R38) * C(SHCO) * C(SO2)
      W(RG70) = K(RG70) * C(SCH3) * C(SHCO)
      W(RG45F) = K(RG45F) * C(SCH2O) * C(SH)
      W(RG45B) = K(RG45B) * C(SH2) * C(SHCO)
      W(RG46) = K(RG46) * C(SCH2O) * C(SO)
      W(RG47) = K(RG47) * C(SCH2O) * C(SOH)
      W(RG71) = K(RG71) * C(SCH3) * C(SCH2O)
      W(RG11) = K(RG11) * C(SCO2) * C(SCH)
      W(RG41F) = K(RG41F) * C(SCO2) * C(SSXCH2)
      W(RG41B) = K(RG41B) * C(STXCH2) * C(SCO2)
      W(RG42) = K(RG42) * C(SCO2) * C(SSXCH2)
      W(RG105) = K(RG105) * C(SC2H) * C(SO)
      W(RG106) = K(RG106) * C(SC2H) * C(SOH)
      W(RG107) = K(RG107) * C(SC2H) * C(SO2)
      W(RG108F) = K(RG108F) * C(SC2H) * C(SH2)
      W(RG108B) = K(RG108B) * C(SH) * C(SC2H2)
      W(RG115F) = K(RG115F) * C(SC2H2) * C(SH)
      W(RG115B) = K(RG115B) * C(SC2H3)
      W(RG116) = K(RG116) * C(SC2H2) * C(SO)
      W(RG117) = K(RG117) * C(SC2H2) * C(SO)
      W(RG118) = K(RG118) * C(SC2H2) * C(SO)
      W(RG119) = K(RG119) * C(SC2H2) * C(SOH)
      W(RG121) = K(RG121) * C(SC2H2) * C(SOH)
      W(RG122) = K(RG122) * C(SC2H2) * C(SOH)
      W(RR004F) = K(RR004F) * C(SC2H2) * C(SSXCH2)
      W(RR004B) = K(RR004B) * C(SH) * C(SC3H3)
      W(RG129) = K(RG129) * C(SC2H3) * C(SH)
      W(RG130) = K(RG130) * C(SC2H3) * C(SH)
      W(RG131) = K(RG131) * C(SC2H3) * C(SO)
      W(RG132) = K(RG132) * C(SC2H3) * C(SOH)
      W(RG133) = K(RG133) * C(SC2H3) * C(SO2)
      W(RG134) = K(RG134) * C(SC2H3) * C(SO2)
      W(RG135) = K(RG135) * C(SC2H3) * C(SO2)
      W(RR011) = K(RR011) * C(SC2H3) * C(SHCO)
      W(RR013) = K(RR013) * C(SC2H3) * C(SCH3)
      W(RR015F) = K(RR015F) * C(SC2H3) * C(SCH3)
      W(RR015B) = K(RR015B) * C(SH) * C(SAXC3H5)
      W(RG155F) = K(RG155F) * C(SC2H4) * C(SH)
      W(RG155B) = K(RG155B) * C(SC2H5)
      W(RG156) = K(RG156) * C(SC2H4) * C(SH)
      W(RG157) = K(RG157) * C(SC2H4) * C(SO)
      W(RG158) = K(RG158) * C(SC2H4) * C(SO)
      W(RG159) = K(RG159) * C(SC2H4) * C(SO)
      W(RG160) = K(RG160) * C(SC2H4) * C(SOH)
      W(RG161) = K(RG161) * C(SC2H4) * C(SOH)
      W(RG162) = K(RG162) * C(SC2H4) * C(SCH3)
      W(RG163F) = K(RG163F) * C(SC2H4) * C(SCH3)
      W(RG163B) = K(RG163B) * C(SNXC3H7)
      W(RG164F) = K(RG164F) * C(SC2H5) * C(SH)
      W(RG164B) = K(RG164B) * C(SC2H6)
      W(RG167) = K(RG167) * C(SC2H5) * C(SO)
      W(RG171) = K(RG171) * C(SC2H5) * C(SO2)
      W(RR034) = K(RR034) * C(SC2H5) * C(SHCO)
      W(RR035) = K(RR035) * C(SC2H5) * C(SHO2)
      W(RR037) = K(RR037) * C(SC2H5) * C(SHO2)
      W(RG173F) = K(RG173F) * C(SC2H6)
      W(RG173B) = K(RG173B) * C(SCH3) * C(SCH3)
      W(RG174) = K(RG174) * C(SC2H6) * C(SH)
      W(RG175) = K(RG175) * C(SC2H6) * C(SO)
      W(RG176) = K(RG176) * C(SC2H6) * C(SOH)
      W(RG178) = K(RG178) * C(SC2H6) * C(SCH3)
      W(RG10) = K(RG10) * C(SHCCO)
      W(RG109) = K(RG109) * C(SHCCO) * C(SH)
      W(RG110) = K(RG110) * C(SHCCO) * C(SO)
      W(RG111) = K(RG111) * C(SHCCO) * C(SO2)
      W(RG114) = K(RG114) * C(SHCCO) * C(SHCCO)
      W(RR009) = K(RR009) * C(SHCCO) * C(SC2H2)
      W(RR022) = K(RR022) * C(SHCCO) * C(SCH3)
      W(RR024) = K(RR024) * C(SHCCO) * C(SOH)
      W(RR055F) = K(RR055F) * C(SC3H3) * C(SH)
      W(RR055B) = K(RR055B) * C(SPXC3H4)
      W(RR056F) = K(RR056F) * C(SC3H3) * C(SH)
      W(RR056B) = K(RR056B) * C(SAXC3H4)
      W(RR058) = K(RR058) * C(SC3H3) * C(SOH)
      W(RR061) = K(RR061) * C(SC3H3) * C(SO)
      W(RR062) = K(RR062) * C(SC3H3) * C(SO2)
      W(RR063) = K(RR063) * C(SC3H3) * C(SHO2)
      W(RLP01F) = K(RLP01F) * C(SC3H3) * C(SC3H3)
      W(RLP01B) = K(RLP01B) * C(SA1XC6H6)
      W(RP011F) = K(RP011F) * C(SC3H3) * C(SC3H3)
      W(RP011B) = K(RP011B) * C(SH) * C(SA1XXC6H5)
      W(RCP11F) = K(RCP11F) * C(SC3H3) * C(SC2H2)
      W(RCP11B) = K(RCP11B) * C(SC5H5)
      W(RR005F) = K(RR005F) * C(SPXC3H4) * C(SH)
      W(RR005B) = K(RR005B) * C(SCH3) * C(SC2H2)
      W(RR078F) = K(RR078F) * C(SPXC3H4) * C(SH)
      W(RR078B) = K(RR078B) * C(SH2) * C(SC3H3)
      W(RR080) = K(RR080) * C(SPXC3H4) * C(SOH)
      W(RR081) = K(RR081) * C(SPXC3H4) * C(SCH3)
      W(RR088) = K(RR088) * C(SPXC3H4) * C(SO)
      W(RR089) = K(RR089) * C(SPXC3H4) * C(SO)
      W(RR094) = K(RR094) * C(SPXC3H4) * C(SOH)
      W(RR006F) = K(RR006F) * C(SAXC3H4) * C(SH)
      W(RR006B) = K(RR006B) * C(SCH3) * C(SC2H2)
      W(RR071F) = K(RR071F) * C(SAXC3H4)
      W(RR071B) = K(RR071B) * C(SPXC3H4)
      W(RR073F) = K(RR073F) * C(SAXC3H4) * C(SH)
      W(RR073B) = K(RR073B) * C(SH) * C(SPXC3H4)
      W(RR083) = K(RR083) * C(SAXC3H4) * C(SH)
      W(RR084) = K(RR084) * C(SAXC3H4) * C(SOH)
      W(RR085) = K(RR085) * C(SAXC3H4) * C(SCH3)
      W(RR074F) = K(RR074F) * C(SAXC3H4) * C(SH)
      W(RR074B) = K(RR074B) * C(SAXC3H5)
      W(RR103) = K(RR103) * C(SAXC3H5) * C(SH)
      W(RR104) = K(RR104) * C(SAXC3H5) * C(SOH)
      W(RR105) = K(RR105) * C(SAXC3H5) * C(SCH3)
      W(RR111) = K(RR111) * C(SAXC3H5) * C(SO2)
      W(RP008) = K(RP008) * C(SAXC3H5) * C(SC3H3)
      W(RG183F) = K(RG183F) * C(SC3H6) * C(SH)
      W(RG183B) = K(RG183B) * C(SNXC3H7)
      W(RR014F) = K(RR014F) * C(SC3H6)
      W(RR014B) = K(RR014B) * C(SCH3) * C(SC2H3)
      W(RR016F) = K(RR016F) * C(SAXC3H5) * C(SH)
      W(RR016B) = K(RR016B) * C(SC3H6)
      W(RR141F) = K(RR141F) * C(SC3H6) * C(SH)
      W(RR141B) = K(RR141B) * C(SCH3) * C(SC2H4)
      W(RR142) = K(RR142) * C(SC3H6) * C(SO)
      W(RR143) = K(RR143) * C(SC3H6) * C(SO)
      W(RR144) = K(RR144) * C(SC3H6) * C(SH)
      W(RR145) = K(RR145) * C(SC3H6) * C(SO)
      W(RR146) = K(RR146) * C(SC3H6) * C(SOH)
      W(RR148) = K(RR148) * C(SC3H6) * C(SCH3)
      W(RHP40) = K(RHP40) * C(SC4H8) * C(SH)
      W(RHP43) = K(RHP43) * C(SC4H8)
      W(RCP06) = K(RCP06) * C(SC5H5) * C(SHO2)
      W(RCP12) = K(RCP12) * C(SC5H5) * C(SC5H5)
      W(RCP13) = K(RCP13) * C(SC5H5) * C(SO)
      W(RCP16) = K(RCP16) * C(SC5H5) * C(SHO2)
      W(RCP18) = K(RCP18) * C(SC5H5) * C(SHO2)
      W(RCP20) = K(RCP20) * C(SC5H5) * C(SOH)
      W(RCP01F) = K(RCP01F) * C(SC5H6)
      W(RCP01B) = K(RCP01B) * C(SH) * C(SC5H5)
      W(RCP02) = K(RCP02) * C(SC5H6) * C(SH)
      W(RCP03) = K(RCP03) * C(SC5H6) * C(SH)
      W(RCP04) = K(RCP04) * C(SC5H6) * C(SO)
      W(RCP05) = K(RCP05) * C(SC5H6) * C(SOH)
      W(RCP08) = K(RCP08) * C(SC5H6) * C(SCH3)
      W(RHP31) = K(RHP31) * C(SC5H10) * C(SH)
      W(RHP33) = K(RHP33) * C(SC5H10)
      W(RHP29) = K(RHP29) * C(SC5H11)
      W(RHP30) = K(RHP30) * C(SC5H11)
      W(RHP32) = K(RHP32) * C(SC5H11)
      W(ROX12) = K(ROX12) * C(SA1XXC6H5) * C(SO2)
      W(ROX13) = K(ROX13) * C(SA1XXC6H5) * C(SO2)
      W(ROX14) = K(ROX14) * C(SA1XXC6H5) * C(SO)
      W(ROX15) = K(ROX15) * C(SA1XXC6H5) * C(SOH)
      W(ROX16) = K(ROX16) * C(SA1XXC6H5) * C(SHO2)
      W(RK012F) = K(RK012F) * C(SA1XXC6H5) * C(SC2H2)
      W(RK012B) = K(RK012B) * C(SA1C2H2XC8H7)
      W(ROX00F) = K(ROX00F) * C(SA1XC6H6)
      W(ROX00B) = K(ROX00B) * C(SH) * C(SA1XXC6H5)
      W(ROX05) = K(ROX05) * C(SA1XC6H6) * C(SH)
      W(ROX06) = K(ROX06) * C(SA1XC6H6) * C(SOH)
      W(ROX17) = K(ROX17) * C(SA1XC6H6) * C(SCH3)
      W(ROX08) = K(ROX08) * C(SA1XC6H6) * C(SO)
      W(ROX11) = K(ROX11) * C(SA1XC6H6) * C(SOH)
      W(RT05F) = K(RT05F) * C(SA1CH2XC7H7)
      W(RT05B) = K(RT05B) * C(SC2H2) * C(SC5H5)
      W(RT04F) = K(RT04F) * C(SA1CH2XC7H7) * C(SH)
      W(RT04B) = K(RT04B) * C(SCH3) * C(SA1XXC6H5)
      W(RT16) = K(RT16) * C(SA1CH2XC7H7) * C(SO)
      W(RT18) = K(RT18) * C(SA1CH2XC7H7) * C(SHO2)
      W(RT19) = K(RT19) * C(SA1CH2XC7H7) * C(SC3H3)
      W(RT01F) = K(RT01F) * C(SA1CH3XC7H8) * C(SH)
      W(RT01B) = K(RT01B) * C(SCH3) * C(SA1XC6H6)
      W(RT02F) = K(RT02F) * C(SA1CH3XC7H8)
      W(RT02B) = K(RT02B) * C(SH) * C(SA1CH2XC7H7)
      W(RT03F) = K(RT03F) * C(SA1CH3XC7H8)
      W(RT03B) = K(RT03B) * C(SCH3) * C(SA1XXC6H5)
      W(RT07) = K(RT07) * C(SA1CH3XC7H8) * C(SH)
      W(RT08) = K(RT08) * C(SA1CH3XC7H8) * C(SOH)
      W(RT13) = K(RT13) * C(SA1CH3XC7H8) * C(SCH3)
      W(RHP10) = K(RHP10) * C(SC7H15)
      W(RHP11) = K(RHP11) * C(SC7H15)
      W(RHP12) = K(RHP12) * C(SC7H15)
      W(RHP13) = K(RHP13) * C(SC7H15)
      W(RHP00) = K(RHP00) * C(SNXC7H16)
      W(RHP01) = K(RHP01) * C(SNXC7H16)
      W(RHP02) = K(RHP02) * C(SNXC7H16) * C(SH)
      W(RHP03) = K(RHP03) * C(SNXC7H16) * C(SO)
      W(RHP04) = K(RHP04) * C(SNXC7H16) * C(SOH)
      W(RHP07) = K(RHP07) * C(SNXC7H16) * C(SCH3)
      W(RT33) = K(RT33) * C(SA1CHOXC7H6O)
      W(RT34) = K(RT34) * C(SA1CHOXC7H6O) * C(SH)
      W(RT38) = K(RT38) * C(SA1CHOXC7H6O) * C(SCH3)
      W(RK100F) = K(RK100F) * C(SA1C2HYXC8H5) * C(SC2H2)
      W(RK100B) = K(RK100B) * C(SA2XXC10H7)
      W(RP106F) = K(RP106F) * C(SA1C2HYXC8H5) * C(SC2H4)
      W(RP106B) = K(RP106B) * C(SH) * C(SA2XC10H8)
      W(RP016F) = K(RP016F) * C(SA1C2HXC8H6)
      W(RP016B) = K(RP016B) * C(SH) * C(SA1C2HYXC8H5)
      W(RK017F) = K(RK017F) * C(SA1C2HXC8H6) * C(SH)
      W(RK017B) = K(RK017B) * C(SH2) * C(SA1C2HYXC8H5)
      W(RP018F) = K(RP018F) * C(SA1C2HXC8H6) * C(SOH)
      W(RP018B) = K(RP018B) * C(SH2O) * C(SA1C2HYXC8H5)
      W(RLP06F) = K(RLP06F) * C(SA1C2H2XC8H7)
      W(RLP06B) = K(RLP06B) * C(SH) * C(SA1C2HXC8H6)
      W(RLP07F) = K(RLP07F) * C(SA1C2H2XC8H7) * C(SC2H2)
      W(RLP07B) = K(RLP07B) * C(SH) * C(SA2XC10H8)
      W(RP108F) = K(RP108F) * C(SA2XC10H8)
      W(RP108B) = K(RP108B) * C(SH) * C(SA2XXC10H7)
      W(RK109F) = K(RK109F) * C(SA2XC10H8) * C(SH)
      W(RK109B) = K(RK109B) * C(SH2) * C(SA2XXC10H7)
      W(RK110F) = K(RK110F) * C(SA2XC10H8) * C(SOH)
      W(RK110B) = K(RK110B) * C(SH2O) * C(SA2XXC10H7)


      CDOT(SN2) = - W(RG28F) + W(RG28F) - W(RG28B) & 
	    + W(RG28B)
      CDOT(SSXCH2) = - W(RG28F) + W(RG28B) - W(RG34F) & 
	    + W(RG34B) - W(RG35) - W(RG36) & 
	    - W(RG37) - W(RG38F) + W(RG38B) & 
	    - W(RG39) + W(RG57F) - W(RG57B) & 
	    + W(RG95) - W(RG40F) + W(RG40B) & 
	    - W(RG41F) + W(RG41B) - W(RG42) & 
	    - W(RR004F) + W(RR004B) + W(RG109)
      CDOT(STXCH2) = W(RG28F) - W(RG28B) + W(RG06F) & 
	    - W(RG06B) - W(RG15) - W(RG16) & 
	    - W(RG17F) + W(RG17B) - W(RG18F) & 
	    + W(RG18B) - W(RG19) - W(RG20) & 
	    - W(RG21) + W(RG38F) - W(RG38B) & 
	    + W(RG55F) - W(RG55B) - W(RG72) & 
	    - W(RG94F) + W(RG94B) + W(RG40F) & 
	    - W(RG40B) + W(RG41F) - W(RG41B) & 
	    + W(RG117) + W(RG158)
      CDOT(SO) = - W(R2F) + W(R2B) - W(R10) & 
	    + W(R4F) - W(R4B) + W(R1F) & 
	    - W(R1B) + W(R15) - W(R17) & 
	    - W(RG04) + W(RG09) - W(RG15) & 
	    + W(RG20) - W(RG52) - W(RG53) & 
	    + W(RG58) - W(RG91F) + W(RG91B) & 
	    - W(R27) - W(R33) - W(R34) & 
	    - W(RG46) - W(RG105) - W(RG116) & 
	    - W(RG117) - W(RG118) - W(RG131) & 
	    + W(RG134) - W(RG157) - W(RG158) & 
	    - W(RG159) - W(RG167) - W(RG175) & 
	    - W(RG110) - W(RR061) - W(RR088) & 
	    - W(RR089) - W(RR142) - W(RR143) & 
	    - W(RR145) - W(RCP13) - W(RCP04) & 
	    + W(ROX12) - W(ROX14) - W(ROX08)
      CDOT(SO) = CDOT(SO) - W(RT16) - W(RHP03)
      CDOT(SH2) = - W(R2F) + W(R2B) - W(R3F) & 
	    + W(R3B) + W(R13) - W(RG06F) & 
	    + W(RG06B) - W(RG18F) + W(RG18B) & 
	    - W(RG34F) + W(RG34B) + W(RG37) & 
	    + W(RG39) + W(RG53) + W(RG54) & 
	    + W(RG90F) - W(RG90B) + W(R32) & 
	    + W(RG45F) - W(RG45B) - W(RG108F) & 
	    + W(RG108B) + W(RG130) + W(RG156) & 
	    + W(RG174) + W(RR078F) - W(RR078B) & 
	    + W(RR083) + W(RR103) + W(RR144) & 
	    + W(RCP02) + W(ROX05) + W(RT07) & 
	    + W(RHP02) + W(RT34) + W(RK017F) & 
	    - W(RK017B) + W(RK109F) - W(RK109B)
      CDOT(SH) = W(R2F) - W(R2B) - W(R10) & 
	    - W(R9) + W(R3F) - W(R3B) & 
	    - W(R12) - W(R1F) + W(R1B) & 
	    - W(R13) - W(R15) - W(R16) & 
	    + W(RG04) + W(RG05) + W(RG06F) & 
	    - W(RG06B) + W(RG08) + W(RG15) & 
	    + W(RG16) + W(RG18F) - W(RG18B) & 
	    + 2 * W(RG19) + W(RG21) + W(RG34F) & 
	    - W(RG34B) + W(RG35) - W(RG51F) & 
	    + W(RG51B) + W(RG52) + W(RG53) & 
	    + W(RG58) + W(RG65) + W(RG69) & 
	    + W(RG72) + W(RG74F) - W(RG74B) & 
	    - W(RG90F) + W(RG90B) + W(RG93) & 
	    + W(R28F) - W(R28B) + W(R29F) & 
	    - W(R29B) - W(R32) + W(R34)
      CDOT(SH) = CDOT(SH) + W(R36F) - W(R36B) + W(R37F) & 
	    - W(R37B) - W(RG45F) + W(RG45B) & 
	    + W(RG106) + W(RG108F) - W(RG108B) & 
	    - W(RG115F) + W(RG115B) + W(RG116) & 
	    + W(RR004F) - W(RR004B) - W(RG129) & 
	    - W(RG130) + W(RR015F) - W(RR015B) & 
	    - W(RG155F) + W(RG155B) - W(RG156) & 
	    + W(RG157) - W(RG164F) + W(RG164B) & 
	    - W(RG174) - W(RG109) + W(RG110) & 
	    - W(RR055F) + W(RR055B) - W(RR056F) & 
	    + W(RR056B) + W(RR061) + W(RP011F) & 
	    - W(RP011B) - W(RR005F) + W(RR005B) & 
	    - W(RR078F) + W(RR078B) - W(RR006F) & 
	    + W(RR006B) - W(RR073F) + W(RR073F) & 
	    - W(RR073B) + W(RR073B) - W(RR083)
      CDOT(SH) = CDOT(SH) - W(RR074F) + W(RR074B) - W(RR103) & 
	    + 2 * W(RP008) - W(RG183F) + W(RG183B) & 
	    - W(RR016F) + W(RR016B) - W(RR141F) & 
	    + W(RR141B) - W(RR144) - W(RHP40) & 
	    + 2 * W(RCP12) + W(RCP20) + W(RCP01F) & 
	    - W(RCP01B) - W(RCP02) - W(RCP03) & 
	    - W(RHP31) + W(RHP30) + W(ROX13) & 
	    + W(ROX00F) - W(ROX00B) - W(ROX05) & 
	    + W(ROX08) + W(ROX11) - W(RT04F) & 
	    + W(RT04B) + W(RT16) + W(RT18) & 
	    + 2 * W(RT19) - W(RT01F) + W(RT01B) & 
	    + W(RT02F) - W(RT02B) - W(RT07) & 
	    - W(RHP02) - W(RT34) + W(RP106F) & 
	    - W(RP106B) + W(RP016F) - W(RP016B) & 
	    - W(RK017F) + W(RK017B) + W(RLP06F)
      CDOT(SH) = CDOT(SH) - W(RLP06B) + W(RLP07F) - W(RLP07B) & 
	    + W(RP108F) - W(RP108B) - W(RK109F) & 
	    + W(RK109B)
      CDOT(SOH) = W(R2F) - W(R2B) + W(R10) & 
	    - W(R9) - W(R3F) + W(R3B) & 
	    - 2 * W(R4F) + 2 * W(R4B) + W(R1F) & 
	    - W(R1B) + 2 * W(R16) + W(R17) & 
	    - W(R18) - W(R19) - W(RG05) & 
	    - W(RG16) - W(RG17F) + W(RG17B) & 
	    + W(RG21) + W(RG35) - W(RG54) & 
	    - W(RG55F) + W(RG55B) - W(RG57F) & 
	    + W(RG57B) + W(RG59) + W(RG65) & 
	    + W(RG91F) - W(RG91B) - W(RG92F) & 
	    + W(RG92B) - W(R28F) + W(R28B) & 
	    - W(R29F) + W(R29B) + W(R31) & 
	    + W(R33) - W(R35) + W(RG46) & 
	    - W(RG47) - W(RG106) + W(RG118) & 
	    - W(RG119) - W(RG121) - W(RG122)
      CDOT(SOH) = CDOT(SOH) - W(RG132) - W(RG160) - W(RG161) & 
	    + W(RR037) + W(RG175) - W(RG176) & 
	    + W(RG111) - W(RR024) - W(RR058) & 
	    + W(RR063) - W(RR080) - W(RR094) & 
	    - W(RR084) - W(RR104) + W(RR111) & 
	    + W(RR145) - W(RR146) + W(RCP16) & 
	    - W(RCP20) + W(RCP04) - W(RCP05) & 
	    - W(ROX15) + W(ROX16) - W(ROX06) & 
	    - W(ROX11) + W(RT18) - W(RT08) & 
	    + W(RHP03) - W(RHP04) - W(RP018F) & 
	    + W(RP018B) - W(RK110F) + W(RK110B)
      CDOT(SH2O) = W(R9) + W(R3F) - W(R3B) & 
	    + W(R4F) - W(R4B) + W(R15) & 
	    + W(R18) + W(R19) - W(RG08) & 
	    + W(RG17F) - W(RG17B) + W(RG36) & 
	    - W(RG37) - W(RG38F) + W(RG38F) & 
	    - W(RG38B) + W(RG38B) - W(RG39) & 
	    + W(RG55F) - W(RG55B) + W(RG57F) & 
	    - W(RG57B) + W(RG92F) - W(RG92B) & 
	    + W(R35) - W(R37F) + W(R37F) & 
	    - W(R37B) + W(R37B) + W(RG47) & 
	    + W(RG119) + W(RG132) + W(RG160) & 
	    + W(RG176) + W(RR080) + W(RR084) & 
	    + W(RR104) + W(RR146) + W(RCP18) & 
	    + W(RCP05) + W(ROX06) + W(RT08) & 
	    + W(RHP04) + W(RP018F) - W(RP018B)
      CDOT(SH2O) = CDOT(SH2O) + W(RK110F) - W(RK110B)
      CDOT(SO2) = - W(R12) - W(R1F) + W(R1B) & 
	    + W(R13) + W(R17) + W(R18) & 
	    + W(R19) - W(RG09) - W(RG19) & 
	    - W(RG20) - W(RG21) - W(RG35) & 
	    - W(RG36) - W(RG58) - W(RG59) & 
	    + W(RG66) - W(R38) - W(RG107) & 
	    - W(RG133) - W(RG134) - W(RG135) & 
	    - W(RG171) + W(RR035) - W(RG111) & 
	    - W(RR062) - W(RR111) + W(RCP06) & 
	    - W(ROX12) - W(ROX13)
      CDOT(SHO2) = W(R12) - W(R13) - W(R15) & 
	    - W(R16) - W(R17) - W(R18) & 
	    - W(R19) - W(RG65) - W(RG66) & 
	    - W(R31) + W(R38) + W(RG133) & 
	    + W(RG171) - W(RR035) - W(RR037) & 
	    - W(RR063) - W(RCP06) - W(RCP16) & 
	    - W(RCP18) - W(ROX16) - W(RT18)
      CDOT(SCH) = - W(RG04) - W(RG05) - W(RG06F) & 
	    + W(RG06B) - W(RG08) - W(RG09) & 
	    + W(RG17F) - W(RG17B) - W(RG69) & 
	    - W(RG93) - W(RG11) + W(RG105) & 
	    + W(RG10)
      CDOT(SCO) = W(RG04) + W(RG21) + W(RG35) & 
	    + W(RG36) + W(RG53) - W(R27) & 
	    - W(R28F) + W(R28B) - W(R29F) & 
	    + W(R29B) - W(R31) - W(RG40F) & 
	    + W(RG40F) - W(RG40B) + W(RG40B) & 
	    + W(R32) + W(R33) + W(R35) & 
	    + W(R36F) - W(R36B) + W(R37F) & 
	    - W(R37B) + W(R38) + W(RG70) & 
	    + W(RG11) + W(RG42) + W(RG105) & 
	    + W(RG107) + W(RG117) + W(RG121) & 
	    + W(RG122) + W(RG131) + W(RG134) & 
	    + W(RR011) + W(RG157) + W(RR034) & 
	    + W(RG10) + W(RG109) + 2 * W(RG110) & 
	    + 2 * W(RG111) + 2 * W(RG114) + W(RR009) & 
	    + W(RR022) + W(RR058) + W(RR061)
      CDOT(SCO) = CDOT(SCO) + 2 * W(RR062) + W(RR063) + W(RR089) & 
	    + W(RR094) + W(RR142) + W(RCP13) & 
	    + W(RCP16) + W(RCP18) + W(RCP20) & 
	    + W(ROX12) + 2 * W(ROX13) + W(ROX14) & 
	    + W(ROX15) + W(ROX16) + W(ROX08) & 
	    + W(ROX11) + W(RT34) + W(RT38)
      CDOT(SHCO) = W(RG05) + W(RG09) + W(RG15) & 
	    - W(R32) - W(R33) - W(R34) & 
	    - W(R35) - W(R36F) + W(R36B) & 
	    - W(R37F) + W(R37B) - W(R38) & 
	    - W(RG70) + W(RG45F) - W(RG45B) & 
	    + W(RG46) + W(RG47) + W(RG71) & 
	    + W(RG11) + W(RG107) + W(RG135) & 
	    - W(RR011) + W(RG159) - W(RR034) & 
	    + 2 * W(RR024) + W(RR143) + W(RT33)
      CDOT(SCH2O) = W(RG08) + W(RG16) + W(RG20) & 
	    + W(RG37) + W(RG39) + W(RG52) & 
	    + W(RG54) + W(RG58) + W(RG59) & 
	    + W(RG65) - W(RG45F) + W(RG45B) & 
	    - W(RG46) - W(RG47) - W(RG71) & 
	    + W(RG42) + W(RG135) + W(RG158) & 
	    + W(RG161) + W(RG167) + W(RR037) & 
	    + W(RR111)
      CDOT(SCH3) = W(RG18F) - W(RG18B) + W(RG34F) & 
	    - W(RG34B) - W(RG51F) + W(RG51B) & 
	    - W(RG52) - W(RG53) - W(RG54) & 
	    - W(RG55F) + W(RG55B) - W(RG57F) & 
	    + W(RG57B) - W(RG58) - W(RG59) & 
	    - W(RG65) - W(RG66) - W(RG69) & 
	    - W(RG72) - 2 * W(RG74F) + 2 * W(RG74B) & 
	    - 2 * W(RG95) + W(RG90F) - W(RG90B) & 
	    + W(RG91F) - W(RG91B) + W(RG92F) & 
	    - W(RG92B) + 2 * W(RG94F) - 2 * W(RG94B) & 
	    - W(RG70) - W(RG71) + W(RG121) & 
	    + W(RG122) + W(RG131) + W(RG134) & 
	    - W(RR013) - W(RR015F) + W(RR015B) & 
	    + W(RG157) + W(RG159) + W(RG161) & 
	    - W(RG162) - W(RG163F) + W(RG163B)
      CDOT(SCH3) = CDOT(SCH3) + W(RG167) + W(RR037) + 2 * W(RG173F) & 
	    - 2 * W(RG173B) - W(RG178) - W(RR022) & 
	    + W(RR062) + W(RR005F) - W(RR005B) & 
	    - W(RR081) + W(RR088) + W(RR006F) & 
	    - W(RR006B) - W(RR085) - W(RR105) & 
	    + W(RR014F) - W(RR014B) + W(RR141F) & 
	    - W(RR141B) + 2 * W(RR142) - W(RR148) & 
	    + W(RHP43) - W(RCP08) - W(ROX17) & 
	    + W(RT04F) - W(RT04B) + W(RT01F) & 
	    - W(RT01B) + W(RT03F) - W(RT03B) & 
	    - W(RT13) - W(RHP07) - W(RT38)
      CDOT(SCO2) = W(RG19) + W(R27) + W(R28F) & 
	    - W(R28B) + W(R29F) - W(R29B) & 
	    + W(R31) + W(R34) - W(RG11) & 
	    - W(RG41F) + W(RG41F) - W(RG41B) & 
	    + W(RG41B) - W(RG42)
      CDOT(SCH4) = W(RG51F) - W(RG51B) + W(RG66) & 
	    + W(RG95) - W(RG90F) + W(RG90B) & 
	    - W(RG91F) + W(RG91B) - W(RG92F) & 
	    + W(RG92B) - W(RG93) - W(RG94F) & 
	    + W(RG94B) + W(RG70) + W(RG71) & 
	    + W(RR013) + W(RG162) + W(RG178) & 
	    + W(RR081) + W(RR085) + W(RR105) & 
	    + W(RR148) + W(RCP08) + W(ROX17) & 
	    + W(RT13) + W(RHP07) + W(RT38)
      CDOT(SC2H3) = W(RG69) + W(RG115F) - W(RG115B) & 
	    - W(RG129) - W(RG130) - W(RG131) & 
	    - W(RG132) - W(RG133) - W(RG134) & 
	    - W(RG135) - W(RR011) - W(RR013) & 
	    - W(RR015F) + W(RR015B) + W(RG156) & 
	    + W(RG160) + W(RG162) + W(RR063) & 
	    + W(RR014F) - W(RR014B) + W(RCP13) & 
	    + W(RCP16) + W(RCP20)
      CDOT(SC2H4) = W(RG72) + W(RG93) + W(RG129) & 
	    + W(RR011) - W(RG155F) + W(RG155B) & 
	    - W(RG156) - W(RG157) - W(RG158) & 
	    - W(RG159) - W(RG160) - W(RG161) & 
	    - W(RG162) - W(RG163F) + W(RG163B) & 
	    + W(RG171) + W(RR022) + W(RR058) & 
	    + W(RR089) + W(RR141F) - W(RR141B) & 
	    + W(RHP40) + W(RHP29) + W(RHP10) & 
	    + W(RHP11) + W(RHP01) - W(RP106F) & 
	    + W(RP106B)
      CDOT(SC2H5) = W(RG74F) - W(RG74B) + W(RG155F) & 
	    - W(RG155B) - W(RG164F) + W(RG164B) & 
	    - W(RG167) - W(RG171) - W(RR034) & 
	    - W(RR035) - W(RR037) + W(RG174) & 
	    + W(RG175) + W(RG176) + W(RG178) & 
	    + W(RR094) + W(RR143) + W(RHP40) & 
	    + W(RHP33) + W(RHP32) + W(RHP11) & 
	    + W(RHP13) + W(RHP00) + W(RHP01)
      CDOT(SC2H) = - W(RG105) - W(RG106) - W(RG107) & 
	    - W(RG108F) + W(RG108B) + W(RG118) & 
	    + W(RG119)
      CDOT(SHCCO) = W(RG106) + W(RG116) - W(RG10) & 
	    - W(RG109) - W(RG110) - W(RG111) & 
	    - 2 * W(RG114) - W(RR009) - W(RR022) & 
	    - W(RR024) + W(RR088)
      CDOT(SC2H2) = W(RG108F) - W(RG108B) - W(RG115F) & 
	    + W(RG115B) - W(RG116) - W(RG117) & 
	    - W(RG118) - W(RG119) - W(RG121) & 
	    - W(RG122) - W(RR004F) + W(RR004B) & 
	    + W(RG130) + W(RG132) + W(RG133) & 
	    + W(RR013) + W(RG114) - W(RR009) & 
	    + W(RR061) - W(RCP11F) + W(RCP11B) & 
	    + W(RR005F) - W(RR005B) + W(RR006F) & 
	    - W(RR006B) + W(RR111) + W(RCP13) & 
	    + W(RCP16) + 2 * W(RCP18) + W(RCP20) & 
	    + W(RCP03) + 2 * W(ROX13) - W(RK012F) & 
	    + W(RK012B) + W(RT05F) - W(RT05B) & 
	    - W(RK100F) + W(RK100B) - W(RLP07F) & 
	    + W(RLP07B)
      CDOT(SC3H3) = W(RR004F) - W(RR004B) + W(RR009) & 
	    - W(RR055F) + W(RR055B) - W(RR056F) & 
	    + W(RR056B) - W(RR058) - W(RR061) & 
	    - W(RR062) - W(RR063) - 2 * W(RLP01F) & 
	    + 2 * W(RLP01B) - 2 * W(RP011F) + 2 * W(RP011B) & 
	    - W(RCP11F) + W(RCP11B) + W(RR078F) & 
	    - W(RR078B) + W(RR080) + W(RR081) & 
	    + W(RR083) + W(RR084) + W(RR085) & 
	    - W(RP008) - W(RT19)
      CDOT(SAXC3H5) = W(RR015F) - W(RR015B) + W(RR074F) & 
	    - W(RR074B) - W(RR103) - W(RR104) & 
	    - W(RR105) - W(RR111) - W(RP008) & 
	    - W(RR016F) + W(RR016B) + W(RR144) & 
	    + W(RR145) + W(RR146) + W(RR148) & 
	    + W(RHP43) + W(RCP03) + W(RHP33)
      CDOT(SNXC3H7) = W(RG163F) - W(RG163B) + W(RG183F) & 
	    - W(RG183B) + W(RHP29) + W(RHP12) & 
	    + W(RHP01)
      CDOT(SC2H6) = W(RG164F) - W(RG164B) + W(RR034) & 
	    + W(RR035) - W(RG173F) + W(RG173B) & 
	    - W(RG174) - W(RG175) - W(RG176) & 
	    - W(RG178)
      CDOT(SPXC3H4) = W(RR055F) - W(RR055B) - W(RR005F) & 
	    + W(RR005B) - W(RR078F) + W(RR078B) & 
	    - W(RR080) - W(RR081) - W(RR088) & 
	    - W(RR089) - W(RR094) + W(RR071F) & 
	    - W(RR071B) + W(RR073F) - W(RR073B)
      CDOT(SAXC3H4) = W(RR056F) - W(RR056B) - W(RR006F) & 
	    + W(RR006B) - W(RR071F) + W(RR071B) & 
	    - W(RR073F) + W(RR073B) - W(RR083) & 
	    - W(RR084) - W(RR085) - W(RR074F) & 
	    + W(RR074B) + W(RR103) + W(RR104) & 
	    + W(RR105)
      CDOT(SA1XC6H6) = W(RLP01F) - W(RLP01B) + W(RP008) & 
	    - W(ROX00F) + W(ROX00B) - W(ROX05) & 
	    - W(ROX06) - W(ROX17) - W(ROX08) & 
	    - W(ROX11) + W(RT01F) - W(RT01B)
      CDOT(SA1XXC6H5) = W(RP011F) - W(RP011B) - W(ROX12) & 
	    - W(ROX13) - W(ROX14) - W(ROX15) & 
	    - W(ROX16) - W(RK012F) + W(RK012B) & 
	    + W(ROX00F) - W(ROX00B) + W(ROX05) & 
	    + W(ROX06) + W(ROX17) + W(RT04F) & 
	    - W(RT04B) + W(RT03F) - W(RT03B) & 
	    + W(RT33) + W(RT34) + W(RT38)
      CDOT(SC5H5) = W(RCP11F) - W(RCP11B) - W(RCP06) & 
	    - 2 * W(RCP12) - W(RCP13) - W(RCP16) & 
	    - W(RCP18) - W(RCP20) + W(RCP01F) & 
	    - W(RCP01B) + W(RCP02) + W(RCP04) & 
	    + W(RCP05) + W(RCP08) + W(ROX12) & 
	    + W(ROX14) + W(ROX16) + W(ROX08) & 
	    + W(RT05F) - W(RT05B)
      CDOT(SC3H6) = - W(RG183F) + W(RG183B) - W(RR014F) & 
	    + W(RR014B) + W(RR016F) - W(RR016B) & 
	    - W(RR141F) + W(RR141B) - W(RR142) & 
	    - W(RR143) - W(RR144) - W(RR145) & 
	    - W(RR146) - W(RR148) + W(RHP32) & 
	    + W(RHP11)
      CDOT(SC4H8) = - W(RHP40) - W(RHP43) + W(RHP12)
      CDOT(SC5H6) = W(RCP06) - W(RCP01F) + W(RCP01B) & 
	    - W(RCP02) - W(RCP03) - W(RCP04) & 
	    - W(RCP05) - W(RCP08) + W(ROX15) & 
	    + W(ROX11)
      CDOT(SA2XC10H8) = W(RCP12) + W(RT19) + W(RP106F) & 
	    - W(RP106B) + W(RLP07F) - W(RLP07B) & 
	    - W(RP108F) + W(RP108B) - W(RK109F) & 
	    + W(RK109B) - W(RK110F) + W(RK110B)
      CDOT(SC5H10) = - W(RHP31) - W(RHP33) + W(RHP30) & 
	    + W(RHP13)
      CDOT(SC5H11) = W(RHP31) - W(RHP29) - W(RHP30) & 
	    - W(RHP32) + W(RHP10) + W(RHP00)
      CDOT(SA1C2H2XC8H7) = W(RK012F) - W(RK012B) - W(RLP06F) & 
	    + W(RLP06B) - W(RLP07F) + W(RLP07B)
      CDOT(SA1CH2XC7H7) = - W(RT05F) + W(RT05B) - W(RT04F) & 
	    + W(RT04B) - W(RT16) - W(RT18) & 
	    - W(RT19) + W(RT02F) - W(RT02B) & 
	    + W(RT07) + W(RT08) + W(RT13)
      CDOT(SA1CHOXC7H6O) = W(RT16) + W(RT18) - W(RT33) & 
	    - W(RT34) - W(RT38)
      CDOT(SA1CH3XC7H8) = - W(RT01F) + W(RT01B) - W(RT02F) & 
	    + W(RT02B) - W(RT03F) + W(RT03B) & 
	    - W(RT07) - W(RT08) - W(RT13)
      CDOT(SC7H15) = - W(RHP10) - W(RHP11) - W(RHP12) & 
	    - W(RHP13) + W(RHP02) + W(RHP03) & 
	    + W(RHP04) + W(RHP07)
      CDOT(SNXC7H16) = - W(RHP00) - W(RHP01) - W(RHP02) & 
	    - W(RHP03) - W(RHP04) - W(RHP07)
      CDOT(SA1C2HYXC8H5) = - W(RK100F) + W(RK100B) - W(RP106F) & 
	    + W(RP106B) + W(RP016F) - W(RP016B) & 
	    + W(RK017F) - W(RK017B) + W(RP018F) & 
	    - W(RP018B)
      CDOT(SA2XXC10H7) = W(RK100F) - W(RK100B) + W(RP108F) & 
	    - W(RP108B) + W(RK109F) - W(RK109B) & 
	    + W(RK110F) - W(RK110B)
      CDOT(SA1C2HXC8H6) = - W(RP016F) + W(RP016B) - W(RK017F) & 
	    + W(RK017B) - W(RP018F) + W(RP018B) & 
	    + W(RLP06F) - W(RLP06B)
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

      DOUBLE PRECISION MM(47)
      INCLUDE 'RedMechF90.h'

      MM(SN2) =  2.80200000e+01
      MM(SSXCH2) =  1.40260000e+01
      MM(STXCH2) =  1.40260000e+01
      MM(SO) =  1.60000000e+01
      MM(SH2) =  2.01600000e+00
      MM(SH) =  1.00800000e+00
      MM(SOH) =  1.70080000e+01
      MM(SH2O) =  1.80160000e+01
      MM(SO2) =  3.20000000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SCH) =  1.30180000e+01
      MM(SCO) =  2.80100000e+01
      MM(SHCO) =  2.90180000e+01
      MM(SCH2O) =  3.00260000e+01
      MM(SCH3) =  1.50340000e+01
      MM(SCO2) =  4.40100000e+01
      MM(SCH4) =  1.60420000e+01
      MM(SC2H3) =  2.70440000e+01
      MM(SC2H4) =  2.80520000e+01
      MM(SC2H5) =  2.90600000e+01
      MM(SC2H) =  2.50280000e+01
      MM(SHCCO) =  4.10280000e+01
      MM(SC2H2) =  2.60360000e+01
      MM(SC3H3) =  3.90540000e+01
      MM(SAXC3H5) =  4.10700000e+01
      MM(SNXC3H7) =  4.30860000e+01
      MM(SC2H6) =  3.00680000e+01
      MM(SPXC3H4) =  4.00620000e+01
      MM(SAXC3H4) =  4.00620000e+01
      MM(SA1XC6H6) =  7.81080000e+01
      MM(SA1XXC6H5) =  7.71000000e+01
      MM(SC5H5) =  6.50900000e+01
      MM(SC3H6) =  4.20780000e+01
      MM(SC4H8) =  5.61040000e+01
      MM(SC5H6) =  6.60980000e+01
      MM(SA2XC10H8) =  1.28164000e+02
      MM(SC5H10) =  7.01300000e+01
      MM(SC5H11) =  7.11380000e+01
      MM(SA1C2H2XC8H7) =  1.03136000e+02
      MM(SA1CH2XC7H7) =  9.11260000e+01
      MM(SA1CHOXC7H6O) =  1.06118000e+02
      MM(SA1CH3XC7H8) =  9.21340000e+01
      MM(SC7H15) =  9.91900000e+01
      MM(SNXC7H16) =  1.00198000e+02
      MM(SA1C2HYXC8H5) =  1.01120000e+02
      MM(SA2XXC10H7) =  1.27156000e+02
      MM(SA1C2HXC8H6) =  1.02128000e+02

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(47)
      INCLUDE 'RedMechF90.h'

      NAMES(SN2)  = 'N2'
      NAMES(SSXCH2)  = 'SXCH2'
      NAMES(STXCH2)  = 'TXCH2'
      NAMES(SO)  = 'O'
      NAMES(SH2)  = 'H2'
      NAMES(SH)  = 'H'
      NAMES(SOH)  = 'OH'
      NAMES(SH2O)  = 'H2O'
      NAMES(SO2)  = 'O2'
      NAMES(SHO2)  = 'HO2'
      NAMES(SCH)  = 'CH'
      NAMES(SCO)  = 'CO'
      NAMES(SHCO)  = 'HCO'
      NAMES(SCH2O)  = 'CH2O'
      NAMES(SCH3)  = 'CH3'
      NAMES(SCO2)  = 'CO2'
      NAMES(SCH4)  = 'CH4'
      NAMES(SC2H3)  = 'C2H3'
      NAMES(SC2H4)  = 'C2H4'
      NAMES(SC2H5)  = 'C2H5'
      NAMES(SC2H)  = 'C2H'
      NAMES(SHCCO)  = 'HCCO'
      NAMES(SC2H2)  = 'C2H2'
      NAMES(SC3H3)  = 'C3H3'
      NAMES(SAXC3H5)  = 'AXC3H5'
      NAMES(SNXC3H7)  = 'NXC3H7'
      NAMES(SC2H6)  = 'C2H6'
      NAMES(SPXC3H4)  = 'PXC3H4'
      NAMES(SAXC3H4)  = 'AXC3H4'
      NAMES(SA1XC6H6)  = 'A1XC6H6'
      NAMES(SA1XXC6H5)  = 'A1XXC6H5'
      NAMES(SC5H5)  = 'C5H5'
      NAMES(SC3H6)  = 'C3H6'
      NAMES(SC4H8)  = 'C4H8'
      NAMES(SC5H6)  = 'C5H6'
      NAMES(SA2XC10H8)  = 'A2XC10H8'
      NAMES(SC5H10)  = 'C5H10'
      NAMES(SC5H11)  = 'C5H11'
      NAMES(SA1C2H2XC8H7)  = 'A1C2H2XC'
      NAMES(SA1CH2XC7H7)  = 'A1CH2XC7'
      NAMES(SA1CHOXC7H6O)  = 'A1CHOXC7'
      NAMES(SA1CH3XC7H8)  = 'A1CH3XC7'
      NAMES(SC7H15)  = 'C7H15'
      NAMES(SNXC7H16)  = 'NXC7H16'
      NAMES(SA1C2HYXC8H5)  = 'A1C2HYXC'
      NAMES(SA2XXC10H7)  = 'A2XXC10H'
      NAMES(SA1C2HXC8H6)  = 'A1C2HXC8'


      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'RedMechF90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 290

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'RedMechF90.h'

      NSPECIES_NONS = 47
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'RedMechF90.h'
      real(DP) :: MUCOEFF(47)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SSXCH2) =  6.92304430e-07
      MUCOEFF(STXCH2) =  6.92304430e-07
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SCH) =  1.27351520e-06
      MUCOEFF(SCO) =  1.06039632e-06
      MUCOEFF(SHCO) =  1.11568663e-06
      MUCOEFF(SCH2O) =  1.13489904e-06
      MUCOEFF(SCH3) =  7.16749611e-07
      MUCOEFF(SCO2) =  1.25056029e-06
      MUCOEFF(SCH4) =  7.61887935e-07
      MUCOEFF(SC2H3) =  8.25781476e-07
      MUCOEFF(SC2H4) =  8.96560348e-07
      MUCOEFF(SC2H5) =  7.77507128e-07
      MUCOEFF(SC2H) =  7.94406422e-07
      MUCOEFF(SHCCO) =  2.73563116e-06
      MUCOEFF(SC2H2) =  8.10245830e-07
      MUCOEFF(SC3H3) =  7.36234631e-07
      MUCOEFF(SAXC3H5) =  6.89211145e-07
      MUCOEFF(SNXC3H7) =  7.05924132e-07
      MUCOEFF(SC2H6) =  7.90876817e-07
      MUCOEFF(SPXC3H4) =  7.45675363e-07
      MUCOEFF(SAXC3H4) =  7.45675363e-07
      MUCOEFF(SA1XC6H6) =  8.43012086e-07
      MUCOEFF(SA1XXC6H5) =  8.37554799e-07
      MUCOEFF(SC5H5) =  7.96430411e-07
      MUCOEFF(SC3H6) =  6.97617690e-07
      MUCOEFF(SC4H8) =  7.72325002e-07
      MUCOEFF(SC5H6) =  8.02573579e-07
      MUCOEFF(SA2XC10H8) =  7.91231307e-07
      MUCOEFF(SC5H10) =  7.41929850e-07
      MUCOEFF(SC5H11) =  8.85961433e-07
      MUCOEFF(SA1C2H2XC8H7) =  7.53008758e-07
      MUCOEFF(SA1CH2XC7H7) =  7.89808619e-07
      MUCOEFF(SA1CHOXC7H6O) =  8.52305406e-07
      MUCOEFF(SA1CH3XC7H8) =  7.94164881e-07
      MUCOEFF(SC7H15) =  6.79914768e-07
      MUCOEFF(SNXC7H16) =  6.83360789e-07
      MUCOEFF(SA1C2HYXC8H5) =  8.20396614e-07
      MUCOEFF(SA2XXC10H7) =  7.88113678e-07
      MUCOEFF(SA1C2HXC8H6) =  8.24475477e-07

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'RedMechF90.h'
      real(DP) :: KOVEREPS(47)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SSXCH2) =  6.94444444e-03
      KOVEREPS(STXCH2) =  6.94444444e-03
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SCH) =  1.25000000e-02
      KOVEREPS(SCO) =  1.01936799e-02
      KOVEREPS(SHCO) =  2.00803213e-03
      KOVEREPS(SCH2O) =  2.00803213e-03
      KOVEREPS(SCH3) =  6.94444444e-03
      KOVEREPS(SCO2) =  4.09836066e-03
      KOVEREPS(SCH4) =  7.07213579e-03
      KOVEREPS(SC2H3) =  4.78468900e-03
      KOVEREPS(SC2H4) =  3.56125356e-03
      KOVEREPS(SC2H5) =  3.96353547e-03
      KOVEREPS(SC2H) =  4.78468900e-03
      KOVEREPS(SHCCO) =  6.66666667e-03
      KOVEREPS(SC2H2) =  4.78468900e-03
      KOVEREPS(SC3H3) =  3.96825397e-03
      KOVEREPS(SAXC3H5) =  3.74812594e-03
      KOVEREPS(SNXC3H7) =  3.74812594e-03
      KOVEREPS(SC2H6) =  3.96353547e-03
      KOVEREPS(SPXC3H4) =  3.96825397e-03
      KOVEREPS(SAXC3H4) =  3.96825397e-03
      KOVEREPS(SA1XC6H6) =  2.15146299e-03
      KOVEREPS(SA1XXC6H5) =  2.15146299e-03
      KOVEREPS(SC5H5) =  2.50000000e-03
      KOVEREPS(SC3H6) =  3.74812594e-03
      KOVEREPS(SC4H8) =  2.89268152e-03
      KOVEREPS(SC5H6) =  2.50000000e-03
      KOVEREPS(SA2XC10H8) =  1.58629442e-03
      KOVEREPS(SC5H10) =  2.58933195e-03
      KOVEREPS(SC5H11) =  2.26893712e-03
      KOVEREPS(SA1C2H2XC8H7) =  1.83083120e-03
      KOVEREPS(SA1CH2XC7H7) =  2.01897840e-03
      KOVEREPS(SA1CHOXC7H6O) =  2.01897840e-03
      KOVEREPS(SA1CH3XC7H8) =  2.01897840e-03
      KOVEREPS(SC7H15) =  2.17580505e-03
      KOVEREPS(SNXC7H16) =  2.17580505e-03
      KOVEREPS(SA1C2HYXC8H5) =  1.86706497e-03
      KOVEREPS(SA2XXC10H7) =  1.58629442e-03
      KOVEREPS(SA1C2HXC8H6) =  1.86706497e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'RedMechF90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(290), C(47), M(6)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 47
!------------------------------------------------------------------
      implicit none
      include 'RedMechF90.h'
      real(DP) :: H(47), CP(47), T

      IF (T.GT.1000.0_DP) THEN
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  2.92664000e+00_DP + T * (  7.43988400e-04_DP &
	   + T * ( -1.89492000e-07_DP + T * (  2.52425950e-11_DP &
	   + T * ( -1.35067020e-15_DP ) ) ) ) ) -9.22797700e+02_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    2.92664000e+00_DP + T * (  1.48797680e-03_DP &
	   + T * ( -5.68476000e-07_DP + T * (  1.00970380e-10_DP &
	   + T * ( -6.75335100e-15_DP ) ) ) ) )
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.29203842e+00_DP + T * (  2.32794318e-03_DP &
	   + T * ( -6.70639823e-07_DP + T * (  1.04476500e-10_DP &
	   + T * ( -6.79432730e-15_DP ) ) ) ) ) +  5.09259997e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    2.29203842e+00_DP + T * (  4.65588637e-03_DP &
	   + T * ( -2.01191947e-06_DP + T * (  4.17906000e-10_DP &
	   + T * ( -3.39716365e-14_DP ) ) ) ) )
      H(STXCH2) =  5.92780550e+02_DP * ( &
	   T * (  2.87410113e+00_DP + T * (  1.82819646e-03_DP &
	   + T * ( -4.69648657e-07_DP + T * (  6.50448872e-11_DP &
	   + T * ( -3.75455134e-15_DP ) ) ) ) ) +  4.62636040e+04_DP )
      CP(STXCH2) =  5.92780550e+02_DP * ( &
	    2.87410113e+00_DP + T * (  3.65639292e-03_DP &
	   + T * ( -1.40894597e-06_DP + T * (  2.60179549e-10_DP &
	   + T * ( -1.87727567e-14_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.56942078e+00_DP + T * ( -4.29870569e-05_DP &
	   + T * (  1.39828196e-08_DP + T * ( -2.50444497e-12_DP &
	   + T * (  2.45667382e-16_DP ) ) ) ) ) +  2.92175791e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.56942078e+00_DP + T * ( -8.59741137e-05_DP &
	   + T * (  4.19484589e-08_DP + T * ( -1.00177799e-11_DP &
	   + T * (  1.22833691e-15_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.33727920e+00_DP + T * ( -2.47012365e-05_DP &
	   + T * (  1.66485593e-07_DP + T * ( -4.48915985e-11_DP &
	   + T * (  4.00510752e-15_DP ) ) ) ) ) -9.50158922e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.33727920e+00_DP + T * ( -4.94024731e-05_DP &
	   + T * (  4.99456778e-07_DP + T * ( -1.79566394e-10_DP &
	   + T * (  2.00255376e-14_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000001e+00_DP + T * ( -1.15421486e-11_DP &
	   + T * (  5.38539827e-15_DP + T * ( -1.18378809e-18_DP &
	   + T * (  9.96394714e-23_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000001e+00_DP + T * ( -2.30842973e-11_DP &
	   + T * (  1.61561948e-14_DP + T * ( -4.73515235e-18_DP &
	   + T * (  4.98197357e-22_DP ) ) ) ) )
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
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  1.76069008e+00_DP + T * (  4.60000041e-03_DP &
	   + T * ( -1.47419604e-06_DP + T * (  2.51603030e-10_DP &
	   + T * ( -1.76771128e-14_DP ) ) ) ) ) -1.39958323e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    1.76069008e+00_DP + T * (  9.20000082e-03_DP &
	   + T * ( -4.42258813e-06_DP + T * (  1.00641212e-09_DP &
	   + T * ( -8.83855640e-14_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  2.97812060e+00_DP + T * (  2.89892600e-03_DP &
	   + T * ( -6.58526667e-07_DP + T * (  7.68244750e-11_DP &
	   + T * ( -3.58348320e-15_DP ) ) ) ) ) +  1.65095130e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    2.97812060e+00_DP + T * (  5.79785200e-03_DP &
	   + T * ( -1.97558000e-06_DP + T * (  3.07297900e-10_DP &
	   + T * ( -1.79174160e-14_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  3.85746029e+00_DP + T * (  2.20718513e-03_DP &
	   + T * ( -7.38271347e-07_DP + T * (  1.30872547e-10_DP &
	   + T * ( -9.44168328e-15_DP ) ) ) ) ) -4.87591660e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    3.85746029e+00_DP + T * (  4.41437026e-03_DP &
	   + T * ( -2.21481404e-06_DP + T * (  5.23490188e-10_DP &
	   + T * ( -4.72084164e-14_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  1.65326226e+00_DP + T * (  5.01315495e-03_DP &
	   + T * ( -1.10553746e-06_DP + T * (  1.34120785e-10_DP &
	   + T * ( -6.29393516e-15_DP ) ) ) ) ) -1.00095936e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    1.65326226e+00_DP + T * (  1.00263099e-02_DP &
	   + T * ( -3.31661238e-06_DP + T * (  5.36483138e-10_DP &
	   + T * ( -3.14696758e-14_DP ) ) ) ) )
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
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  3.16780652e+00_DP + T * (  2.37610951e-03_DP &
	   + T * ( -6.12623590e-07_DP + T * (  7.60475630e-11_DP &
	   + T * ( -3.54465540e-15_DP ) ) ) ) ) +  6.71210650e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    3.16780652e+00_DP + T * (  4.75221902e-03_DP &
	   + T * ( -1.83787077e-06_DP + T * (  3.04190252e-10_DP &
	   + T * ( -1.77232770e-14_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  5.62820580e+00_DP + T * (  2.04267005e-03_DP &
	   + T * ( -5.31151567e-07_DP + T * (  7.15651300e-11_DP &
	   + T * ( -3.88156640e-15_DP ) ) ) ) ) +  1.93272150e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    5.62820580e+00_DP + T * (  4.08534010e-03_DP &
	   + T * ( -1.59345470e-06_DP + T * (  2.86260520e-10_DP &
	   + T * ( -1.94078320e-14_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  4.14756964e+00_DP + T * (  2.98083332e-03_DP &
	   + T * ( -7.90982840e-07_DP + T * (  1.16853043e-10_DP &
	   + T * ( -7.22470426e-15_DP ) ) ) ) ) +  2.59359992e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    4.14756964e+00_DP + T * (  5.96166664e-03_DP &
	   + T * ( -2.37294852e-06_DP + T * (  4.67412171e-10_DP &
	   + T * ( -3.61235213e-14_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  6.14915291e+00_DP + T * (  4.67031583e-03_DP &
	   + T * ( -1.25018451e-06_DP + T * (  1.72539079e-10_DP &
	   + T * ( -9.21649988e-15_DP ) ) ) ) ) +  3.83854848e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    6.14915291e+00_DP + T * (  9.34063166e-03_DP &
	   + T * ( -3.75055354e-06_DP + T * (  6.90156316e-10_DP &
	   + T * ( -4.60824994e-14_DP ) ) ) ) )
      H(SAXC3H5) =  2.02443146e+02_DP * ( &
	   T * (  2.28794927e+00_DP + T * (  1.18200788e-02_DP &
	   + T * ( -4.26304833e-06_DP + T * (  8.42096350e-10_DP &
	   + T * ( -6.94898898e-14_DP ) ) ) ) ) +  1.83033514e+04_DP )
      CP(SAXC3H5) =  2.02443146e+02_DP * ( &
	    2.28794927e+00_DP + T * (  2.36401575e-02_DP &
	   + T * ( -1.27891450e-05_DP + T * (  3.36838540e-09_DP &
	   + T * ( -3.47449449e-13_DP ) ) ) ) )
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
      H(SC4H8) =  1.48195138e+02_DP * ( &
	   T * (  3.04470367e+00_DP + T * (  1.63725883e-02_DP &
	   + T * ( -4.84544123e-06_DP + T * (  5.99360043e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.52177534e+03_DP )
      CP(SC4H8) =  1.48195138e+02_DP * ( &
	    3.04470367e+00_DP + T * (  3.27451765e-02_DP &
	   + T * ( -1.45363237e-05_DP + T * (  2.39744017e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H6) =  1.25788072e+02_DP * ( &
	   T * (  2.30537462e-01_DP + T * (  2.04785913e-02_DP &
	   + T * ( -8.05296527e-06_DP + T * (  1.69940870e-09_DP &
	   + T * ( -1.47274884e-13_DP ) ) ) ) ) +  1.43779465e+04_DP )
      CP(SC5H6) =  1.25788072e+02_DP * ( &
	    2.30537462e-01_DP + T * (  4.09571826e-02_DP &
	   + T * ( -2.41588958e-05_DP + T * (  6.79763480e-09_DP &
	   + T * ( -7.36374421e-13_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * (  1.76826275e+00_DP + T * (  3.44571753e-02_DP &
	   + T * ( -1.38107392e-05_DP + T * (  2.94785772e-09_DP &
	   + T * ( -2.57194122e-13_DP ) ) ) ) ) +  1.26883657e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	    1.76826275e+00_DP + T * (  6.89143506e-02_DP &
	   + T * ( -4.14322176e-05_DP + T * (  1.17914309e-08_DP &
	   + T * ( -1.28597061e-12_DP ) ) ) ) )
      H(SC5H10) =  1.18556110e+02_DP * ( &
	   T * (  3.98580522e+00_DP + T * (  2.06214993e-02_DP &
	   + T * ( -6.14634990e-06_DP + T * (  7.65388102e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -5.70112071e+03_DP )
      CP(SC5H10) =  1.18556110e+02_DP * ( &
	    3.98580522e+00_DP + T * (  4.12429986e-02_DP &
	   + T * ( -1.84390497e-05_DP + T * (  3.06155241e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SC5H11) =  1.16876212e+02_DP * ( &
	   T * (  4.88920629e+00_DP + T * (  2.11417269e-02_DP &
	   + T * ( -6.19477000e-06_DP + T * (  7.60311908e-10_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) +  3.43475468e+03_DP )
      CP(SC5H11) =  1.16876212e+02_DP * ( &
	    4.88920629e+00_DP + T * (  4.22834537e-02_DP &
	   + T * ( -1.85843100e-05_DP + T * (  3.04124763e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   T * (  5.85935080e+00_DP + T * (  2.36285729e-02_DP &
	   + T * ( -8.99549110e-06_DP + T * (  1.83827944e-09_DP &
	   + T * ( -1.54980166e-13_DP ) ) ) ) ) +  4.33198974e+04_DP )
      CP(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	    5.85935080e+00_DP + T * (  4.72571459e-02_DP &
	   + T * ( -2.69864733e-05_DP + T * (  7.35311775e-09_DP &
	   + T * ( -7.74900830e-13_DP ) ) ) ) )
      H(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   T * (  3.30049696e+00_DP + T * (  2.40027670e-02_DP &
	   + T * ( -9.28143407e-06_DP + T * (  1.93092839e-09_DP &
	   + T * ( -1.65430827e-13_DP ) ) ) ) ) +  2.17498572e+04_DP )
      CP(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	    3.30049696e+00_DP + T * (  4.80055340e-02_DP &
	   + T * ( -2.78443022e-05_DP + T * (  7.72371356e-09_DP &
	   + T * ( -8.27154136e-13_DP ) ) ) ) )
      H(SA1CHOXC7H6O) =  7.83499501e+01_DP * ( &
	   T * (  1.87355756e+00_DP + T * (  2.63115775e-02_DP &
	   + T * ( -1.05881654e-05_DP + T * (  2.26600767e-09_DP &
	   + T * ( -1.98061225e-13_DP ) ) ) ) ) -7.23603865e+03_DP )
      CP(SA1CHOXC7H6O) =  7.83499501e+01_DP * ( &
	    1.87355756e+00_DP + T * (  5.26231551e-02_DP &
	   + T * ( -3.17644962e-05_DP + T * (  9.06403069e-09_DP &
	   + T * ( -9.90306123e-13_DP ) ) ) ) )
      H(SA1CH3XC7H8) =  9.02418217e+01_DP * ( &
	   T * ( -1.01117220e+00_DP + T * (  2.92650956e-02_DP &
	   + T * ( -1.15865023e-05_DP + T * (  2.45545248e-09_DP &
	   + T * ( -2.13361740e-13_DP ) ) ) ) ) +  3.99363395e+03_DP )
      CP(SA1CH3XC7H8) =  9.02418217e+01_DP * ( &
	   -1.01117220e+00_DP + T * (  5.85301912e-02_DP &
	   + T * ( -3.47595069e-05_DP + T * (  9.82180993e-09_DP &
	   + T * ( -1.06680870e-12_DP ) ) ) ) )
      H(SC7H15) =  8.38223611e+01_DP * ( &
	   T * (  3.74721159e+00_DP + T * (  3.24672581e-02_DP &
	   + T * ( -1.00447008e-05_DP + T * (  1.29354536e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -3.37018357e+03_DP )
      CP(SC7H15) =  8.38223611e+01_DP * ( &
	    3.74721159e+00_DP + T * (  6.49345162e-02_DP &
	   + T * ( -3.01341025e-05_DP + T * (  5.17418142e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
      H(SNXC7H16) =  8.29791014e+01_DP * ( &
	   T * (  5.14079241e+00_DP + T * (  3.26539335e-02_DP &
	   + T * ( -9.82758747e-06_DP + T * (  1.23431681e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) ) -2.72533890e+04_DP )
      CP(SNXC7H16) =  8.29791014e+01_DP * ( &
	    5.14079241e+00_DP + T * (  6.53078671e-02_DP &
	   + T * ( -2.94827624e-05_DP + T * (  4.93726726e-09_DP &
	   + T * (  0.00000000e+00_DP ) ) ) ) )
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
      ELSE
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  3.29867700e+00_DP + T * (  7.04120200e-04_DP &
	   + T * ( -1.32107400e-06_DP + T * (  1.41037875e-09_DP &
	   + T * ( -4.88970800e-13_DP ) ) ) ) ) -1.02089990e+03_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    3.29867700e+00_DP + T * (  1.40824040e-03_DP &
	   + T * ( -3.96322200e-06_DP + T * (  5.64151500e-09_DP &
	   + T * ( -2.44485400e-12_DP ) ) ) ) )
      H(SSXCH2) =  5.92780550e+02_DP * ( &
	   T * (  4.19860411e+00_DP + T * ( -1.18330710e-03_DP &
	   + T * (  2.74432073e-06_DP + T * ( -1.67203995e-09_DP &
	   + T * (  3.88629474e-13_DP ) ) ) ) ) +  5.04968163e+04_DP )
      CP(SSXCH2) =  5.92780550e+02_DP * ( &
	    4.19860411e+00_DP + T * ( -2.36661419e-03_DP &
	   + T * (  8.23296220e-06_DP + T * ( -6.68815981e-09_DP &
	   + T * (  1.94314737e-12_DP ) ) ) ) )
      H(STXCH2) =  5.92780550e+02_DP * ( &
	   T * (  3.76267867e+00_DP + T * (  4.84436072e-04_DP &
	   + T * (  9.31632803e-07_DP + T * ( -9.62727883e-10_DP &
	   + T * (  3.37483438e-13_DP ) ) ) ) ) +  4.60040401e+04_DP )
      CP(STXCH2) =  5.92780550e+02_DP * ( &
	    3.76267867e+00_DP + T * (  9.68872143e-04_DP &
	   + T * (  2.79489841e-06_DP + T * ( -3.85091153e-09_DP &
	   + T * (  1.68741719e-12_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  3.16826710e+00_DP + T * ( -1.63965942e-03_DP &
	   + T * (  2.21435465e-06_DP + T * ( -1.53201656e-09_DP &
	   + T * (  4.22531942e-13_DP ) ) ) ) ) +  2.91222592e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    3.16826710e+00_DP + T * ( -3.27931884e-03_DP &
	   + T * (  6.64306396e-06_DP + T * ( -6.12806624e-09_DP &
	   + T * (  2.11265971e-12_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.34433112e+00_DP + T * (  3.99026037e-03_DP &
	   + T * ( -6.49271700e-06_DP + T * (  5.03930235e-09_DP &
	   + T * ( -1.47522352e-12_DP ) ) ) ) ) -9.17935173e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.34433112e+00_DP + T * (  7.98052075e-03_DP &
	   + T * ( -1.94781510e-05_DP + T * (  2.01572094e-08_DP &
	   + T * ( -7.37611761e-12_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  3.52666409e-13_DP &
	   + T * ( -6.65306547e-16_DP + T * (  5.75204080e-19_DP &
	   + T * ( -1.85546466e-22_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  7.05332819e-13_DP &
	   + T * ( -1.99591964e-15_DP + T * (  2.30081632e-18_DP &
	   + T * ( -9.27732332e-22_DP ) ) ) ) )
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
      H(SCH2O) =  2.76904683e+02_DP * ( &
	   T * (  4.79372315e+00_DP + T * ( -4.95416684e-03_DP &
	   + T * (  1.24406669e-05_DP + T * ( -9.48213152e-09_DP &
	   + T * (  2.63545304e-12_DP ) ) ) ) ) -1.43089567e+04_DP )
      CP(SCH2O) =  2.76904683e+02_DP * ( &
	    4.79372315e+00_DP + T * ( -9.90833369e-03_DP &
	   + T * (  3.73220008e-05_DP + T * ( -3.79285261e-08_DP &
	   + T * (  1.31772652e-11_DP ) ) ) ) )
      H(SCH3) =  5.53035786e+02_DP * ( &
	   T * (  3.65717970e+00_DP + T * (  1.06329895e-03_DP &
	   + T * (  1.81946277e-06_DP + T * ( -1.65452507e-09_DP &
	   + T * (  4.93141480e-13_DP ) ) ) ) ) +  1.64227160e+04_DP )
      CP(SCH3) =  5.53035786e+02_DP * ( &
	    3.65717970e+00_DP + T * (  2.12659790e-03_DP &
	   + T * (  5.45838830e-06_DP + T * ( -6.61810030e-09_DP &
	   + T * (  2.46570740e-12_DP ) ) ) ) )
      H(SCO2) =  1.88919337e+02_DP * ( &
	   T * (  2.35677352e+00_DP + T * (  4.49229839e-03_DP &
	   + T * ( -2.37452090e-06_DP + T * (  6.14797555e-10_DP &
	   + T * ( -2.87399096e-14_DP ) ) ) ) ) -4.83719697e+04_DP )
      CP(SCO2) =  1.88919337e+02_DP * ( &
	    2.35677352e+00_DP + T * (  8.98459677e-03_DP &
	   + T * ( -7.12356269e-06_DP + T * (  2.45919022e-09_DP &
	   + T * ( -1.43699548e-13_DP ) ) ) ) )
      H(SCH4) =  5.18285750e+02_DP * ( &
	   T * (  5.14911468e+00_DP + T * ( -6.83110045e-03_DP &
	   + T * (  1.63817974e-05_DP + T * ( -1.21061692e-08_DP &
	   + T * (  3.33206882e-12_DP ) ) ) ) ) -1.02465983e+04_DP )
      CP(SCH4) =  5.18285750e+02_DP * ( &
	    5.14911468e+00_DP + T * ( -1.36622009e-02_DP &
	   + T * (  4.91453921e-05_DP + T * ( -4.84246767e-08_DP &
	   + T * (  1.66603441e-11_DP ) ) ) ) )
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
      H(SC2H) =  3.32201534e+02_DP * ( &
	   T * (  2.88965733e+00_DP + T * (  6.70498055e-03_DP &
	   + T * ( -9.49231670e-06_DP + T * (  7.36977613e-09_DP &
	   + T * ( -2.18663022e-12_DP ) ) ) ) ) +  6.68393932e+04_DP )
      CP(SC2H) =  3.32201534e+02_DP * ( &
	    2.88965733e+00_DP + T * (  1.34099611e-02_DP &
	   + T * ( -2.84769501e-05_DP + T * (  2.94791045e-08_DP &
	   + T * ( -1.09331511e-11_DP ) ) ) ) )
      H(SHCCO) =  2.02650385e+02_DP * ( &
	   T * (  2.25172140e+00_DP + T * (  8.82751050e-03_DP &
	   + T * ( -7.90970033e-06_DP + T * (  4.31893975e-09_DP &
	   + T * ( -1.01329622e-12_DP ) ) ) ) ) +  2.00594490e+04_DP )
      CP(SHCCO) =  2.02650385e+02_DP * ( &
	    2.25172140e+00_DP + T * (  1.76550210e-02_DP &
	   + T * ( -2.37291010e-05_DP + T * (  1.72757590e-08_DP &
	   + T * ( -5.06648110e-12_DP ) ) ) ) )
      H(SC2H2) =  3.19340144e+02_DP * ( &
	   T * (  8.08681094e-01_DP + T * (  1.16807815e-02_DP &
	   + T * ( -1.18390605e-05_DP + T * (  7.00381092e-09_DP &
	   + T * ( -1.70014595e-12_DP ) ) ) ) ) +  2.64289807e+04_DP )
      CP(SC2H2) =  3.19340144e+02_DP * ( &
	    8.08681094e-01_DP + T * (  2.33615629e-02_DP &
	   + T * ( -3.55171815e-05_DP + T * (  2.80152437e-08_DP &
	   + T * ( -8.50072974e-12_DP ) ) ) ) )
      H(SC3H3) =  2.12893430e+02_DP * ( &
	   T * (  1.40299238e+00_DP + T * (  1.50886664e-02_DP &
	   + T * ( -1.32816458e-05_DP + T * (  7.33836572e-09_DP &
	   + T * ( -1.74110916e-12_DP ) ) ) ) ) +  3.93108220e+04_DP )
      CP(SC3H3) =  2.12893430e+02_DP * ( &
	    1.40299238e+00_DP + T * (  3.01773327e-02_DP &
	   + T * ( -3.98449373e-05_DP + T * (  2.93534629e-08_DP &
	   + T * ( -8.70554579e-12_DP ) ) ) ) )
      H(SAXC3H5) =  2.02443146e+02_DP * ( &
	   T * ( -1.03516444e+00_DP + T * (  1.87521683e-02_DP &
	   + T * ( -1.08793747e-05_DP + T * (  3.69156533e-09_DP &
	   + T * ( -4.87482308e-13_DP ) ) ) ) ) +  1.88792254e+04_DP )
      CP(SAXC3H5) =  2.02443146e+02_DP * ( &
	   -1.03516444e+00_DP + T * (  3.75043366e-02_DP &
	   + T * ( -3.26381242e-05_DP + T * (  1.47662613e-08_DP &
	   + T * ( -2.43741154e-12_DP ) ) ) ) )
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
      H(SC4H8) =  1.48195138e+02_DP * ( &
	   T * ( -8.31372089e-01_DP + T * (  2.26290489e-02_DP &
	   + T * ( -9.78861863e-06_DP + T * (  2.50551090e-09_DP &
	   + T * ( -2.86383360e-13_DP ) ) ) ) ) -1.57875035e+03_DP )
      CP(SC4H8) =  1.48195138e+02_DP * ( &
	   -8.31372089e-01_DP + T * (  4.52580978e-02_DP &
	   + T * ( -2.93658559e-05_DP + T * (  1.00220436e-08_DP &
	   + T * ( -1.43191680e-12_DP ) ) ) ) )
      H(SC5H6) =  1.25788072e+02_DP * ( &
	   T * ( -5.13691194e+00_DP + T * (  3.03476727e-02_DP &
	   + T * ( -1.53517612e-05_DP + T * (  3.21143002e-09_DP &
	   + T * (  1.48242970e-13_DP ) ) ) ) ) +  1.53675713e+04_DP )
      CP(SC5H6) =  1.25788072e+02_DP * ( &
	   -5.13691194e+00_DP + T * (  6.06953453e-02_DP &
	   + T * ( -4.60552837e-05_DP + T * (  1.28457201e-08_DP &
	   + T * (  7.41214852e-13_DP ) ) ) ) )
      H(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   T * ( -8.72434585e+00_DP + T * (  5.26880040e-02_DP &
	   + T * ( -2.67236897e-05_DP + T * (  5.46364935e-09_DP &
	   + T * (  2.84133212e-13_DP ) ) ) ) ) +  1.48059774e+04_DP )
      CP(SA2XC10H8) =  6.48726632e+01_DP * ( &
	   -8.72434585e+00_DP + T * (  1.05376008e-01_DP &
	   + T * ( -8.01710690e-05_DP + T * (  2.18545974e-08_DP &
	   + T * (  1.42066606e-12_DP ) ) ) ) )
      H(SC5H10) =  1.18556110e+02_DP * ( &
	   T * ( -1.06223481e+00_DP + T * (  2.87109147e-02_DP &
	   + T * ( -1.24828963e-05_DP + T * (  3.18412472e-09_DP &
	   + T * ( -3.59219578e-13_DP ) ) ) ) ) -4.46546666e+03_DP )
      CP(SC5H10) =  1.18556110e+02_DP * ( &
	   -1.06223481e+00_DP + T * (  5.74218294e-02_DP &
	   + T * ( -3.74486890e-05_DP + T * (  1.27364989e-08_DP &
	   + T * ( -1.79609789e-12_DP ) ) ) ) )
      H(SC5H11) =  1.16876212e+02_DP * ( &
	   T * ( -9.05255912e-01_DP + T * (  3.05316426e-02_DP &
	   + T * ( -1.36497275e-05_DP + T * (  3.65233675e-09_DP &
	   + T * ( -4.37719230e-13_DP ) ) ) ) ) +  4.83995303e+03_DP )
      CP(SC5H11) =  1.16876212e+02_DP * ( &
	   -9.05255912e-01_DP + T * (  6.10632852e-02_DP &
	   + T * ( -4.09491825e-05_DP + T * (  1.46093470e-08_DP &
	   + T * ( -2.18859615e-12_DP ) ) ) ) )
      H(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   T * ( -6.31199276e+00_DP + T * (  4.75548971e-02_DP &
	   + T * ( -3.18784034e-05_DP + T * (  1.24445052e-08_DP &
	   + T * ( -2.04647434e-12_DP ) ) ) ) ) +  4.57330975e+04_DP )
      CP(SA1C2H2XC8H7) =  8.06153041e+01_DP * ( &
	   -6.31199276e+00_DP + T * (  9.51097942e-02_DP &
	   + T * ( -9.56352102e-05_DP + T * (  4.97780207e-08_DP &
	   + T * ( -1.02323717e-11_DP ) ) ) ) )
      H(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   T * ( -6.07053038e+00_DP + T * (  4.17600754e-02_DP &
	   + T * ( -2.47233361e-05_DP + T * (  7.82884618e-09_DP &
	   + T * ( -8.47341736e-13_DP ) ) ) ) ) +  2.35894712e+04_DP )
      CP(SA1CH2XC7H7) =  9.12400413e+01_DP * ( &
	   -6.07053038e+00_DP + T * (  8.35201507e-02_DP &
	   + T * ( -7.41700083e-05_DP + T * (  3.13153847e-08_DP &
	   + T * ( -4.23670868e-12_DP ) ) ) ) )
      H(SA1CHOXC7H6O) =  7.83499501e+01_DP * ( &
	   T * ( -3.47171048e+00_DP + T * (  3.46445945e-02_DP &
	   + T * ( -1.44201170e-05_DP + T * (  8.59677740e-10_DP &
	   + T * (  9.62020522e-13_DP ) ) ) ) ) -6.14558774e+03_DP )
      CP(SA1CHOXC7H6O) =  7.83499501e+01_DP * ( &
	   -3.47171048e+00_DP + T * (  6.92891889e-02_DP &
	   + T * ( -4.32603509e-05_DP + T * (  3.43871096e-09_DP &
	   + T * (  4.81010261e-12_DP ) ) ) ) )
      H(SA1CH3XC7H8) =  9.02418217e+01_DP * ( &
	   T * ( -4.54072038e+00_DP + T * (  3.42713573e-02_DP &
	   + T * ( -1.19037675e-05_DP + T * ( -1.04849410e-09_DP &
	   + T * (  1.48355959e-12_DP ) ) ) ) ) +  4.64121087e+03_DP )
      CP(SA1CH3XC7H8) =  9.02418217e+01_DP * ( &
	   -4.54072038e+00_DP + T * (  6.85427145e-02_DP &
	   + T * ( -3.57113024e-05_DP + T * ( -4.19397642e-09_DP &
	   + T * (  7.41779795e-12_DP ) ) ) ) )
      H(SC7H15) =  8.38223611e+01_DP * ( &
	   T * ( -3.79155767e-02_DP + T * (  3.78363285e-02_DP &
	   + T * ( -1.35824545e-05_DP + T * (  2.33169736e-09_DP &
	   + T * ( -9.84721490e-14_DP ) ) ) ) ) -2.35605303e+03_DP )
      CP(SC7H15) =  8.38223611e+01_DP * ( &
	   -3.79155767e-02_DP + T * (  7.56726570e-02_DP &
	   + T * ( -4.07473634e-05_DP + T * (  9.32678943e-09_DP &
	   + T * ( -4.92360745e-13_DP ) ) ) ) )
      H(SNXC7H16) =  8.29791014e+01_DP * ( &
	   T * ( -1.26836187e+00_DP + T * (  4.27177910e-02_DP &
	   + T * ( -1.75115595e-05_DP + T * (  4.07364302e-09_DP &
	   + T * ( -4.04789850e-13_DP ) ) ) ) ) -2.56586565e+04_DP )
      CP(SNXC7H16) =  8.29791014e+01_DP * ( &
	   -1.26836187e+00_DP + T * (  8.54355820e-02_DP &
	   + T * ( -5.25346786e-05_DP + T * (  1.62945721e-08_DP &
	   + T * ( -2.02394925e-12_DP ) ) ) ) )
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
      END IF

      END
