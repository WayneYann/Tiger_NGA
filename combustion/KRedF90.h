!-------------------------------------------------------------
! ======= KRedF90.h =======
!-------------------------------------------------------------

     integer, parameter ::   SN2 = 1, SH = 2 & 
     , SO2 = 3, SO = 4 & 
     , SOH = 5, SH2 = 6 & 
     , SH2O = 7, SHO2 = 8 & 
     , SH2O2 = 9, SCO = 10 & 
     , SCO2 = 11, SHCO = 12 & 
     , SCH = 13, STXCH2 = 14 & 
     , SCH2O = 15, SHCCO = 16 & 
     , SCH3 = 17, SCH2CO = 18 & 
     , SC2H2 = 19, SSXCH2 = 20 & 
     , SCH2OH = 21, SCH4 = 22 & 
     , SC2H4 = 23, SC2H = 24 & 
     , SC2H3 = 25, SC3H3 = 26 & 
     , SPXC3H4 = 27, SNXC4H3 = 28 & 
     , SC3H2 = 29, SC4H2 = 30 
     integer, parameter ::   SA1XC6H6 = 31, SA1XXC6H5 = 32 & 
     , SA1C2H2XC8H7 = 33, SA1C2HXC8H6 = 34 & 
     , SA1C2HYXC8H5 = 35, SA2XXC10H7 = 36 & 
     , SA2XC10H8 = 37, SA2C2H2AXC12H9 = 38 & 
     , SA2R5XC12H8 = 39, SC5H5 = 40 & 
     , SC9H8 = 41, SC9H7 = 42 & 
     , SA1CH2XC7H7 = 43, SA1CH3XC7H8 = 44 & 
     , SOXC6H4 = 45, SA2OHXC10H8O = 46 & 
     , SEND = 47 
      integer, parameter ::   R1F = 1, R1B = 2  & 
     , R2F = 3, R2B = 4  & 
     , R3F = 5, R3B = 6  & 
     , R4F = 7, R4B = 8  & 
     , R9F = 9, R9B = 10  & 
     , R10 = 11, R12F = 12  & 
     , R12B = 13, R13 = 14  & 
     , R14F = 15, R14B = 16  & 
     , R16 = 17, R17F = 18  & 
     , R17B = 19, R18F = 20  & 
     , R18B = 21, R19F = 22  & 
     , R19B = 23, R26F = 24  & 
     , R26B = 25, R28F = 26  & 
     , R28B = 27, R29F = 28  & 
     , R29B = 29, R32 = 30 
      integer, parameter ::   R36F = 31, R36B = 32  & 
     , R37F = 33, R37B = 34  & 
     , R38 = 35, RG06F = 36  & 
     , RG06B = 37, RG08 = 38  & 
     , RG09 = 39, RG10 = 40  & 
     , RG11 = 41, RG13 = 42  & 
     , RG15 = 43, RG16 = 44  & 
     , RG17F = 45, RG17B = 46  & 
     , RG18F = 47, RG18B = 48  & 
     , RG19 = 49, RG20 = 50  & 
     , RG21 = 51, RG24F = 52  & 
     , RG24B = 53, RG26 = 54  & 
     , RG28F = 55, RG28B = 56  & 
     , RG34F = 57, RG34B = 58  & 
     , RG35 = 59, RG38F = 60 
      integer, parameter ::   RG38B = 61, RG39 = 62  & 
     , RG40F = 63, RG40B = 64  & 
     , RG41F = 65, RG41B = 66  & 
     , RG42 = 67, RG43 = 68  & 
     , RG45 = 69, RG46 = 70  & 
     , RG47 = 71, RG51F = 72  & 
     , RG51B = 73, RG52 = 74  & 
     , RG53 = 75, RG55F = 76  & 
     , RG55B = 77, RG57F = 78  & 
     , RG57B = 79, RG59 = 80  & 
     , RG71 = 81, RG72 = 82  & 
     , RG85 = 83, RG86 = 84  & 
     , RG90F = 85, RG90B = 86  & 
     , RG91F = 87, RG91B = 88  & 
     , RG92F = 89, RG92B = 90 
      integer, parameter ::   RG93 = 91, RG106 = 92  & 
     , RG107 = 93, RG108F = 94  & 
     , RG108B = 95, RG109 = 96  & 
     , RG110 = 97, RG115F = 98  & 
     , RG115B = 99, RG116 = 100  & 
     , RG117 = 101, RG118F = 102  & 
     , RG118B = 103, RG119F = 104  & 
     , RG119B = 105, RG121 = 106  & 
     , RG122 = 107, RG123F = 108  & 
     , RG123B = 109, RG124F = 110  & 
     , RG124B = 111, RG127F = 112  & 
     , RG127B = 113, RG129F = 114  & 
     , RG129B = 115, RG130F = 116  & 
     , RG130B = 117, RG135 = 118  & 
     , RG156F = 119, RG156B = 120 
      integer, parameter ::   RG158 = 121, RG159 = 122  & 
     , RG160F = 123, RG160B = 124  & 
     , RG162F = 125, RG162B = 126  & 
     , RR004F = 127, RR004B = 128  & 
     , RR005F = 129, RR005B = 130  & 
     , RR008F = 131, RR008B = 132  & 
     , RR009 = 133, RR013 = 134  & 
     , RR022 = 135, RR024 = 136  & 
     , RR025 = 137, RR029F = 138  & 
     , RR029B = 139, RR040 = 140  & 
     , RR041 = 141, RR043 = 142  & 
     , RR053F = 143, RR053B = 144  & 
     , RR054F = 145, RR054B = 146  & 
     , RR055F = 147, RR055B = 148  & 
     , RR058 = 149, RR059F = 150 
      integer, parameter ::   RR059B = 151, RR062 = 152  & 
     , RR078F = 153, RR078B = 154  & 
     , RR080F = 155, RR080B = 156  & 
     , RR081F = 157, RR081B = 158  & 
     , RR089 = 159, RR093 = 160  & 
     , RH07 = 161, RH08 = 162  & 
     , RH10F = 163, RH10B = 164  & 
     , RH12 = 165, RP010F = 166  & 
     , RP010B = 167, RP011F = 168  & 
     , RP011B = 169, RK012F = 170  & 
     , RK012B = 171, RP015F = 172  & 
     , RP015B = 173, RP016F = 174  & 
     , RP016B = 175, RK017F = 176  & 
     , RK017B = 177, RP018F = 178  & 
     , RP018B = 179, RK020F = 180 
      integer, parameter ::   RK020B = 181, RK100F = 182  & 
     , RK100B = 183, RP104F = 184  & 
     , RP104B = 185, RP105 = 186  & 
     , RP106 = 187, RP108F = 188  & 
     , RP108B = 189, RK109F = 190  & 
     , RK109B = 191, RK110F = 192  & 
     , RK110B = 193, RK114F = 194  & 
     , RK114B = 195, RP118 = 196  & 
     , RP120 = 197, RK200F = 198  & 
     , RK200B = 199, RCP16F = 200  & 
     , RCP16B = 201, RCP17 = 202  & 
     , RI00F = 203, RI00B = 204  & 
     , RI01F = 205, RI01B = 206  & 
     , RI02F = 207, RI02B = 208  & 
     , RI03F = 209, RI03B = 210 
      integer, parameter ::   RI06F = 211, RI06B = 212  & 
     , RI09F = 213, RI09B = 214  & 
     , RI18 = 215, RI23 = 216  & 
     , RI25 = 217, RI26 = 218  & 
     , RT01F = 219, RT01B = 220  & 
     , RT02F = 221, RT02B = 222  & 
     , RT03F = 223, RT03B = 224  & 
     , RT04F = 225, RT04B = 226  & 
     , RT05F = 227, RT05B = 228  & 
     , RT07F = 229, RT07B = 230  & 
     , RT08F = 231, RT08B = 232  & 
     , RT14 = 233, RT20 = 234  & 
     , RST00 = 235, ROX00F = 236  & 
     , ROX00B = 237, ROX01F = 238  & 
     , ROX01B = 239, ROX02F = 240 
      integer, parameter ::   ROX02B = 241, ROX03F = 242  & 
     , ROX03B = 243, ROX04F = 244  & 
     , ROX04B = 245, ROX16F = 246  & 
     , ROX16B = 247, ROX31 = 248  & 
     , ROX32F = 249, ROX32B = 250  & 
     , ROX41 = 251, ROX63 = 252 & 
     , REND = 253 
     integer, parameter ::   MM2 = 1, MM3 = 2 & 
     , MM5 = 3, MM8 = 4 & 
     , MM9 = 5 & 
     , MEND = 6 

integer, parameter :: DP=kind(1.0d0)


	real(DP) :: A(REND),N(REND),E(REND)


	integer :: IH

      DATA (A(IH),IH=1,30)  /  2.6400000000D+13, 5.2766713024D+10 & 
     , 4.5900000000D+01, 2.7517349858D+01 & 
     , 1.7300000000D+05, 1.9032628245D+06 & 
     , 3.9700000000D+01, 7.2853303338D+02 & 
     , 4.4000000000D+16, 1.1314226588D+20 & 
     , 9.4300000000D+12, 6.3300000000D+16 & 
     , 9.1184255084D+19, 3.6490000000D+03 & 
     , 2.0100000000D+14, 7.1253840510D+22 & 
     , 7.4900000000D+10, 4.0000000000D+10 & 
     , 3.8909647883D+09, 2.3800000000D+10 & 
     , 4.2484744234D+10, 1.0000000000D+13 & 
     , 1.7850732872D+13, 2.6700000000D+38 & 
     , 9.9613411288D+36, 8.0000000000D+08 & 
     , 4.3798719198D+15, 8.7800000000D+07 & 
     , 4.8069094320D+14, 1.2000000000D+11 /
      DATA (A(IH),IH=31,60)  /  1.8700000000D+14, 3.5306519794D+10 & 
     , 2.2400000000D+15, 4.2292301785D+11 & 
     , 1.2000000000D+07, 1.0800000000D+11 & 
     , 1.7581777781D+12, 5.7100000000D+09 & 
     , 6.7100000000D+10, 1.2130000000D+35 & 
     , 1.9000000000D+11, 3.4670000000D+26 & 
     , 8.0000000000D+10, 2.0000000000D+10 & 
     , 1.1300000000D+04, 7.6364598733D+03 & 
     , 5.0000000000D+02, 2.0974636059D+05 & 
     , 5.8000000000D+09, 2.4000000000D+09 & 
     , 5.0000000000D+09, 2.6900000000D+30 & 
     , 9.7210519966D+41, 1.6000000000D+12 & 
     , 1.5000000000D+10, 1.0936408744D+10 & 
     , 7.0000000000D+10, 2.1409471364D+13 & 
     , 2.8000000000D+10, 3.0000000000D+10 /
      DATA (A(IH),IH=61,90)  /  2.1872817488D+10, 6.8200000000D+07 & 
     , 9.0000000000D+09, 6.5618452463D+09 & 
     , 7.0000000000D+09, 5.1036574138D+09 & 
     , 1.4000000000D+10, 1.2680000000D+33 & 
     , 5.7400000000D+04, 3.9000000000D+10 & 
     , 3.4300000000D+06, 3.4700000000D+35 & 
     , 1.5077439481D+41, 5.0600000000D+10 & 
     , 3.3700000000D+10, 5.6000000000D+04 & 
     , 1.4686430659D+03, 6.4400000000D+14 & 
     , 2.3164910420D+13, 5.8700000000D+08 & 
     , 3.3200000000D+00, 1.0000000000D+11 & 
     , 1.6540000000D+04, 1.1830000000D+05 & 
     , 6.6000000000D+05, 3.5503024345D+02 & 
     , 1.0200000000D+06, 3.2893910502D+02 & 
     , 1.0000000000D+05, 5.9179879483D+02 /
      DATA (A(IH),IH=91,120)  /  6.0000000000D+10, 2.0000000000D+10 & 
     , 1.0000000000D+10, 3.3100000000D+03 & 
     , 9.6742642418D+06, 1.0000000000D+11 & 
     , 1.0000000000D+11, 6.3400000000D+28 & 
     , 1.4048536324D+31, 8.1000000000D+03 & 
     , 1.2500000000D+04, 1.8100000000D+10 & 
     , 8.8241850413D+13, 2.6300000000D+03 & 
     , 9.8996184977D+00, 7.5300000000D+03 & 
     , 1.2800000000D+06, 5.0000000000D+10 & 
     , 3.8326334707D+07, 1.5000000000D+06 & 
     , 4.0698109539D-01, 7.5000000000D+09 & 
     , 6.3247186168D+07, 1.4000000000D+27 & 
     , 7.7814512532D+32, 3.0000000000D+10 & 
     , 3.1644553172D+10, 4.5800000000D+13 & 
     , 1.2700000000D+02, 5.3406078444D-02 /
      DATA (A(IH),IH=121,150)  /  7.1500000000D+01, 3.8900000000D+05 & 
     , 1.3100000000D-04, 6.0605344705D-07 & 
     , 2.2700000000D+02, 1.7745629538D+02 & 
     , 1.9000000000D+11, 2.5098929533D+16 & 
     , 3.4600000000D+09, 1.8216946136D+04 & 
     , 7.8000000000D+10, 5.6362388680D+19 & 
     , 1.0000000000D+08, 9.0300000000D+09 & 
     , 5.0000000000D+10, 1.0000000000D+10 & 
     , 5.0000000000D+09, 7.5000000000D+09 & 
     , 1.0687278636D+10, 1.0000000000D+10 & 
     , 1.2500000000D+08, 5.0000000000D+10 & 
     , 2.8000000000D+27, 3.8610572228D+32 & 
     , 1.1000000000D+07, 1.8645069743D+04 & 
     , 7.9400000000D+26, 8.1610655573D+31 & 
     , 1.2800000000D+06, 1.1300000000D+02 /
      DATA (A(IH),IH=151,180)  /  2.1071838653D+00, 1.7000000000D+02 & 
     , 8.5000000000D+01, 1.9329120587D-01 & 
     , 8.0500000000D+02, 2.0139176673D+01 & 
     , 4.2200000000D+11, 1.7839562769D+12 & 
     , 6.2500000000D+03, 7.5300000000D+03 & 
     , 9.5600000000D+09, 2.0600000000D+04 & 
     , 1.3700000000D+36, 4.7217121748D+38 & 
     , 3.3000000000D+09, 1.0700000000D+42 & 
     , 1.5425608195D+58, 5.7700000000D+34 & 
     , 5.1317962822D+44, 3.2900000000D+03 & 
     , 2.8151364037D+14, 9.4500000000D-06 & 
     , 2.7559045191D-05, 4.3000000000D+60 & 
     , 6.4261715627D+54, 1.2900000000D+05 & 
     , 4.5060198908D+01, 7.8000000000D+00 & 
     , 2.9974411044D-02, 7.1800000000D+10 /
      DATA (A(IH),IH=181,210)  /  4.1382935254D+06, 1.3400000000D+01 & 
     , 3.3117446183D+14, 3.8000000000D+04 & 
     , 1.0608179381D+14, 3.6000000000D+14 & 
     , 3.6200000000D+25, 8.6000000000D+60 & 
     , 6.5580381913D+54, 2.6500000000D+05 & 
     , 4.7232496168D+01, 9.6300000000D-01 & 
     , 1.8883129559D-03, 1.2800000000D+03 & 
     , 2.2606983071D+14, 7.2000000000D+14 & 
     , 3.6200000000D+25, 2.8800000000D+11 & 
     , 2.1196844222D+08, 6.8700000000D+52 & 
     , 5.5097990876D+64, 6.3900000000D+26 & 
     , 1.0000000000D+10, 5.2485026205D+25 & 
     , 1.7300000000D+68, 1.8203611359D+63 & 
     , 2.8000000000D+10, 6.8863488094D+07 & 
     , 3.1600000000D+01, 4.5268324245D+08 /
      DATA (A(IH),IH=211,240)  /  3.0800000000D+03, 8.3336328636D+01 & 
     , 1.8000000000D-04, 8.2296634648D-04 & 
     , 1.9600000000D+28, 4.3200000000D+36 & 
     , 2.3800000000D+15, 1.1900000000D+15 & 
     , 2.3100000000D+03, 4.5426319227D-03 & 
     , 1.5600000000D+13, 1.2813196505D+10 & 
     , 4.3500000000D+22, 5.2773986505D+10 & 
     , 5.8300000000D+64, 8.6112522502D+55 & 
     , 8.2000000000D+14, 4.3614434871D+01 & 
     , 6.4700000000D-03, 1.2421004032D-03 & 
     , 1.7700000000D+02, 3.7383361476D+02 & 
     , 1.5060000000D+14, 4.3200000000D+36 & 
     , 1.0000000000D+10, 1.2900000000D+61 & 
     , 7.9583820932D+54, 1.0000000000D+84 & 
     , 4.1465511326D+78, 5.0000000000D+75 /
      DATA (A(IH),IH=241,252)  /  2.7802587455D+86, 6.0200000000D+05 & 
     , 8.6806271172D+01, 2.3400000000D+01 & 
     , 3.7121295904D-02, 3.8900000000D-06 & 
     , 1.4511636996D-05, 1.8300000000D+10 & 
     , 2.2000000000D-01, 2.7062147246D+04 & 
     , 8.6200000000D+15, 1.7600000000D-01 & 
     /
      DATA (N(IH),IH=1,30)  /  -0.67, -0.234291 & 
     , 2.7, 2.66385 & 
     , 1.51, 1.40625 & 
     , 2.4, 2.33241 & 
     , -2, -1.75351 & 
     , -1, -1.4 & 
     , -1.41184, 2.068 & 
     , -0.58, -1.59836 & 
     , 0, 0 & 
     , 0.325921, 0 & 
     , 0.258327, 0 & 
     , 0.258327, -7 & 
     , -6.49679, 0.14 & 
     , -1.13743, 0.03 & 
     , -1.24743, 0 /
      DATA (N(IH),IH=31,60)  /  -1, -0.767288 & 
     , -1, -0.767288 & 
     , 0.81, 0 & 
     , -0.292468, 0 & 
     , 0, -5.181 & 
     , 0, -2.669 & 
     , 0, 0 & 
     , 2, 2.18872 & 
     , 2, 1.4875 & 
     , 0, 0 & 
     , 0, -5.11 & 
     , -6.85567, 0 & 
     , 0, -0.0622863 & 
     , 0, -0.574789 & 
     , 0, 0 /
      DATA (N(IH),IH=61,90)  /  -0.0622863, 0.25 & 
     , 0, -0.0622863 & 
     , 0, -0.0622863 & 
     , 0, -5.201 & 
     , 1.9, 0 & 
     , 1.18, -6.3 & 
     , -6.46997, 0 & 
     , 0, 1.6 & 
     , 2.00876, -1.34 & 
     , -0.868957, 0 & 
     , 2.81, 0 & 
     , 1.628, 1.359 & 
     , 1.62, 2.14021 & 
     , 1.5, 1.98406 & 
     , 1.6, 2.01646 /
      DATA (N(IH),IH=91,120)  /  0, 0 & 
     , 0, 2.26 & 
     , 1.66346, 0 & 
     , 0, -4.66 & 
     , -4.48812, 2 & 
     , 2, 0 & 
     , -0.560387, 2.14 & 
     , 2.63279, 1.55 & 
     , 0.73, 0 & 
     , 0.59717, 1.38 & 
     , 2.9634, 0 & 
     , 0.493423, -3.86 & 
     , -3.99563, 0 & 
     , 0.178352, -1.39 & 
     , 2.75, 3.23587 /
      DATA (N(IH),IH=121,150)  /  2.47, 1.36 & 
     , 4.2, 4.58212 & 
     , 2, 1.96566 & 
     , 0, -1.35613 & 
     , 0.44, 1.58056 & 
     , 0, -1.29954 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0.0769611, 0 & 
     , 0, 0 & 
     , -3.86, -3.97975 & 
     , 1.13, 1.59998 & 
     , -5.06, -5.06898 & 
     , 0.73, 2.28 /
      DATA (N(IH),IH=151,180)  /  2.64624, 1.7 & 
     , 2.7, 3.05922 & 
     , 2.22, 2.47547 & 
     , 0, -0.16099 & 
     , 2, 1.55 & 
     , 0, 2 & 
     , -7.87, -7.9297 & 
     , -0.25, -9.57 & 
     , -11.0277, -7 & 
     , -8.24617, 2.05 & 
     , 0.27021, 4.47 & 
     , 4.39414, -12.48 & 
     , -12.2493, 1.89 & 
     , 2.47097, 2.68 & 
     , 3.15722, 1.02 /
      DATA (N(IH),IH=181,210)  /  1.32078, 2.5 & 
     , 0.979368, 1.62 & 
     , 0.397642, -1.44 & 
     , -4.24, -12.48 & 
     , -12.2468, 1.87 & 
     , 2.45348, 3.02 & 
     , 3.49973, 2.05 & 
     , 0.192174, -1.44 & 
     , -4.24, 0.23 & 
     , 0.811179, -12.5 & 
     , -13.5715, -4.03 & 
     , 0, -2.03247 & 
     , -15.16, -14.9678 & 
     , 0, 0.542404 & 
     , 2.48, 1.36023 /
      DATA (N(IH),IH=211,240)  /  2, 2.43866 & 
     , 4, 4.0222 & 
     , -4.85, -7.74 & 
     , 0, 0 & 
     , 2.17, 3.57723 & 
     , 0.68, 0.254459 & 
     , -1.73, -0.111273 & 
     , -14.15, -12.1057 & 
     , 0, 2.2454 & 
     , 3.98, 3.90469 & 
     , 2.39, 2.21095 & 
     , -0.596, -7.74 & 
     , 0, -12.48 & 
     , -12.2685, -18.87 & 
     , -18.4977, -19.31 /
      DATA (N(IH),IH=241,252)  /  -20.3484, 1.8 & 
     , 2.36173, 2.68 & 
     , 3.13799, 4.57 & 
     , 4.52848, 0 & 
     , 3.25, 2.13353 & 
     , -0.61, 3.25 & 
     /
      DATA (E(IH),IH=1,30)  /  71300000, 479604.2721 & 
      , 26190000, 20100582.44 & 
      , 14350000, 77394524.74 & 
      , -8830000, 60303942.3 & 
      , 0, 496136012.5 & 
      , 0, 0 & 
      , 204667815.3, -4574000 & 
      , -9590000, 205148485.1 & 
      , 2660000, 0 & 
      , 222334254.9, -2090000 & 
      , 289378197.2, 72510000 & 
      , 363978197.2, 157320000 & 
      , 287203668.2, 30760000 & 
      , 138461819.8, -70000 & 
      , 107631819.8, 0 /
      DATA (E(IH),IH=31,60)  /  71130000, 5356428.13 & 
      , 71130000, 5356428.13 & 
      , -3040000, 13010000 & 
      , 1691026.963, -3160000 & 
      , 0, 320776000 & 
      , 66070000, 371385000 & 
      , 0, 0 & 
      , 12550000, 86913497.78 & 
      , 30250000, 61407685.18 & 
      , 6280000, 6280000 & 
      , 6280000, 29690000 & 
      , 364655095.4, 49970000 & 
      , 2510000, 39812987.98 & 
      , 0, 68460673.16 & 
      , 0, 0 /
      DATA (E(IH),IH=61,90)  /  37302987.98, -3910000 & 
      , 0, 37302987.98 & 
      , 0, 37302987.98 & 
      , 0, 150084000 & 
      , 11470000, 14810000 & 
      , -1870000, 21230000 & 
      , 463233395.7, 0 & 
      , 0, 22680000 & 
      , 54566839.56, 5930000 & 
      , 513851.5818, 57910000 & 
      , 24520000, 0 & 
      , 14319000, 12643000 & 
      , 45360000, 36448092.09 & 
      , 35980000, 20978674.53 & 
      , 13050000, 67182616.83 /
      DATA (E(IH),IH=91,120)  /  0, 0 & 
      , -3160000, 3770000 & 
      , 127161378.7, 0 & 
      , 0, 15820000 & 
      , 162621237.8, 7950000 & 
      , 7950000, 0 & 
      , 129480796.2, 71380000 & 
      , 11033146.07, 8810000 & 
      , 10790000, 33470000 & 
      , 22499718.93, 2570000 & 
      , 131854077.6, 8370000 & 
      , 60444243.67, 13890000 & 
      , 480409256.8, 0 & 
      , 286290250, 4250000 & 
      , 48740000, 15312230.95 /
      DATA (E(IH),IH=121,150)  /  3890000, 3710000 & 
      , -3600000, 26016755.69 & 
      , 38490000, 13974138.86 & 
      , 0, 105055737.5 & 
      , 22860000, 46357562.82 & 
      , 0, 250786354.2 & 
      , 12550000, -3200000 & 
      , 0, 0 & 
      , 0, 54390000 & 
      , 52331626.84, 0 & 
      , 4180000, 0 & 
      , 13890000, 439116845.8 & 
      , 58280000, 66144642.01 & 
      , 20340000, 393338860.7 & 
      , 10790000, 10320000 /
      DATA (E(IH),IH=151,180)  /  81229166.75, 6280000 & 
      , 24020000, 84112627.12 & 
      , 3100000, 126237151.9 & 
      , 93120000, 162124535 & 
      , 7950000, 8810000 & 
      , 130120000, 7950000 & 
      , 64610000, 199616108.7 & 
      , 9940000, 71190000 & 
      , 677453672.1, 131820000 & 
      , 262035013.4, 13230000 & 
      , 191386104.5, 18710000 & 
      , 28239401.86, 619590000 & 
      , 132912595.6, 73550000 & 
      , 19964083.4, 3070000 & 
      , 12528608.14, 161810000 /
      DATA (E(IH),IH=181,210)  /  30384498.09, 5370000 & 
      , 386582470.9, 18570000 & 
      , 279884457.4, 65930000 & 
      , 99850000, 619550000 & 
      , 121345107.2, 71530000 & 
      , 6416595.001, 18300000 & 
      , 16231119.74, 8080000 & 
      , 181300070.2, 65930000 & 
      , 99850000, 71240000 & 
      , 85000673.5, 176000000 & 
      , 487565243.7, 147300000 & 
      , 0, 521300387.2 & 
      , 486900000, 147515954.9 & 
      , 9450000, 103157442.7 & 
      , 46280000, 103588858.9 /
      DATA (E(IH),IH=211,240)  /  0, 156751967.4 & 
      , 0, 102619350.6 & 
      , 103650000, 99800000 & 
      , 377480000, 377480000 & 
      , 17420000, 50953235.53 & 
      , 373240000, -1780318.268 & 
      , 436010000, -6505423.169 & 
      , 285890000, 218394895.1 & 
      , 337550000, 46091406.02 & 
      , 14160000, 72231169.53 & 
      , -2520000, 118595694.3 & 
      , 160103000, 99800000 & 
      , 0, 619590000 & 
      , 143541341.3, 376980000 & 
      , 50909087.01, 284180000 /
      DATA (E(IH),IH=241,252)  /  537844017.3, 68420000 & 
      , 25462829.1, 3070000 & 
      , 23157353.84, 22000000 & 
      , 56045263, 18950000 & 
      , 23390000, 26652213.94 & 
      , 310110000, 23390000 & 
     /