!-------------------------------------------------------------
! ======= DME39F90.h =======
!-------------------------------------------------------------

     integer, parameter ::   SN2 = 1, SH = 2 & 
     , SO2 = 3, SO = 4 & 
     , SOH = 5, SH2 = 6 & 
     , SH2O = 7, SHO2 = 8 & 
     , SH2O2 = 9, SCO = 10 & 
     , SCO2 = 11, SHCO = 12 & 
     , SCH3 = 13, SCH4 = 14 & 
     , SCH2O = 15, SCH3O = 16 & 
     , SC2H6 = 17, SCH2OH = 18 & 
     , SC2H5 = 19, SCH2 = 20 & 
     , SCH2X = 21, SC2H4 = 22 & 
     , SCH3HCO = 23, SC2H2 = 24 & 
     , SC2H3 = 25, SCH3OCH3 = 26 & 
     , SCH3OCH2 = 27, SCH3OCH2O = 28 & 
     , SCH3OCHO = 29, SCH3OCO = 30 
     integer, parameter ::   SRO2 = 31, SROH = 32 & 
     , SQOOH = 33, SO2QOOH = 34 & 
     , SHO2QHO = 35, SOCH2OCHO = 36 & 
     , SHOCH2OCO = 37, SHOCH2O = 38 & 
     , SHCOOH = 39 & 
     , SEND = 40
 
      integer, parameter ::   R1F = 1, R1B = 2  & 
     , R2F = 3, R2B = 4  & 
     , R3F = 5, R3B = 6  & 
     , R4F = 7, R4B = 8  & 
     , R5F = 9, R5B = 10  & 
     , R6F = 11, R6B = 12  & 
     , R7F = 13, R7B = 14  & 
     , R8F = 15, R8B = 16  & 
     , R9F = 17, R9B = 18  & 
     , R10F = 19, R10B = 20  & 
     , R11F = 21, R11B = 22  & 
     , R12F = 23, R12B = 24  & 
     , R13F = 25, R13B = 26  & 
     , R14F = 27, R14B = 28  & 
     , R15F = 29, R15B = 30 
      integer, parameter ::   R16F = 31, R16B = 32  & 
     , R17F = 33, R17B = 34  & 
     , R18F = 35, R18B = 36  & 
     , R19F = 37, R19B = 38  & 
     , R20F = 39, R20B = 40  & 
     , R21F = 41, R21B = 42  & 
     , R22F = 43, R22B = 44  & 
     , R23F = 45, R23B = 46  & 
     , R24F = 47, R24B = 48  & 
     , R25F = 49, R25B = 50  & 
     , R26F = 51, R26B = 52  & 
     , R27F = 53, R27B = 54  & 
     , R28F = 55, R28B = 56  & 
     , R29F = 57, R29B = 58  & 
     , R30F = 59, R30B = 60 
      integer, parameter ::   R31F = 61, R31B = 62  & 
     , R32F = 63, R32B = 64  & 
     , R33F = 65, R33B = 66  & 
     , R34F = 67, R34B = 68  & 
     , R35F = 69, R35B = 70  & 
     , R36F = 71, R36B = 72  & 
     , R37F = 73, R37B = 74  & 
     , R38F = 75, R38B = 76  & 
     , R39F = 77, R39B = 78  & 
     , R40F = 79, R40B = 80  & 
     , R41F = 81, R41B = 82  & 
     , R42F = 83, R42B = 84  & 
     , R43F = 85, R43B = 86  & 
     , R44F = 87, R44B = 88  & 
     , R45F = 89, R45B = 90 
      integer, parameter ::   R46F = 91, R46B = 92  & 
     , R47F = 93, R47B = 94  & 
     , R48F = 95, R48B = 96  & 
     , R49F = 97, R49B = 98  & 
     , R50F = 99, R50B = 100  & 
     , R51F = 101, R51B = 102  & 
     , R52F = 103, R52B = 104  & 
     , R53F = 105, R53B = 106  & 
     , R54F = 107, R54B = 108  & 
     , R55F = 109, R55B = 110  & 
     , R56F = 111, R56B = 112  & 
     , R57F = 113, R57B = 114  & 
     , R58F = 115, R58B = 116  & 
     , R59F = 117, R59B = 118  & 
     , R60F = 119, R60B = 120 
      integer, parameter ::   R61F = 121, R61B = 122  & 
     , R62F = 123, R62B = 124  & 
     , R63F = 125, R63B = 126  & 
     , R64F = 127, R64B = 128  & 
     , R65F = 129, R65B = 130  & 
     , R66F = 131, R66B = 132  & 
     , R67F = 133, R67B = 134  & 
     , R68F = 135, R68B = 136  & 
     , R69F = 137, R69B = 138  & 
     , R70F = 139, R70B = 140  & 
     , R71F = 141, R71B = 142  & 
     , R72F = 143, R72B = 144  & 
     , R73F = 145, R73B = 146  & 
     , R74F = 147, R74B = 148  & 
     , R75F = 149, R75B = 150 
      integer, parameter ::   R76F = 151, R76B = 152  & 
     , R77F = 153, R77B = 154  & 
     , R78F = 155, R78B = 156  & 
     , R79F = 157, R79B = 158  & 
     , R80F = 159, R80B = 160  & 
     , R81F = 161, R81B = 162  & 
     , R82F = 163, R82B = 164  & 
     , R83F = 165, R83B = 166  & 
     , R84F = 167, R84B = 168  & 
     , R85F = 169, R85B = 170  & 
     , R86F = 171, R86B = 172  & 
     , R87F = 173, R87B = 174  & 
     , R88F = 175, R88B = 176  & 
     , R89F = 177, R89B = 178  & 
     , R90F = 179, R90B = 180 
      integer, parameter ::   R91F = 181, R91B = 182  & 
     , R92F = 183, R92B = 184  & 
     , R93F = 185, R93B = 186  & 
     , R94F = 187, R94B = 188  & 
     , R95F = 189, R95B = 190  & 
     , R96F = 191, R96B = 192  & 
     , R97F = 193, R97B = 194  & 
     , R98F = 195, R98B = 196  & 
     , R99F = 197, R99B = 198  & 
     , R100F = 199, R100B = 200  & 
     , R101F = 201, R101B = 202  & 
     , R102F = 203, R102B = 204  & 
     , R103F = 205, R103B = 206  & 
     , R104F = 207, R104B = 208  & 
     , R105F = 209, R105B = 210 
      integer, parameter ::   R106F = 211, R106B = 212  & 
     , R107F = 213, R107B = 214  & 
     , R108F = 215, R108B = 216  & 
     , R109F = 217, R109B = 218  & 
     , R110F = 219, R110B = 220  & 
     , R111F = 221, R111B = 222  & 
     , R112F = 223, R112B = 224  & 
     , R113F = 225, R113B = 226  & 
     , R114F = 227, R114B = 228  & 
     , R115F = 229, R115B = 230  & 
     , R116F = 231, R116B = 232  & 
     , R117F = 233, R117B = 234  & 
     , R118F = 235, R118B = 236  & 
     , R119F = 237, R119B = 238  & 
     , R120F = 239, R120B = 240 
      integer, parameter ::   R121F = 241, R121B = 242  & 
     , R122F = 243, R122B = 244  & 
     , R123F = 245, R123B = 246  & 
     , R124F = 247, R124B = 248  & 
     , R125F = 249, R125B = 250  & 
     , R126F = 251, R126B = 252  & 
     , R127F = 253, R127B = 254  & 
     , R128F = 255, R128B = 256  & 
     , R129F = 257, R129B = 258  & 
     , R130F = 259, R130B = 260  & 
     , R131F = 261, R131B = 262  & 
     , R132F = 263, R132B = 264  & 
     , R133F = 265, R133B = 266  & 
     , R134F = 267, R134B = 268  & 
     , R135F = 269, R135B = 270 
      integer, parameter ::   R136F = 271, R136B = 272  & 
     , R137F = 273, R137B = 274  & 
     , R138F = 275, R138B = 276  & 
     , R139F = 277, R139B = 278  & 
     , R140F = 279, R140B = 280  & 
     , R141F = 281, R141B = 282  & 
     , R142F = 283, R142B = 284  & 
     , R143F = 285, R143B = 286  & 
     , R144F = 287, R144B = 288  & 
     , R145F = 289, R145B = 290  & 
     , R146F = 291, R146B = 292  & 
     , R147F = 293, R147B = 294  & 
     , R148F = 295, R148B = 296  & 
     , R149F = 297, R149B = 298  & 
     , R150F = 299, R150B = 300 
      integer, parameter ::   R151F = 301, R151B = 302  & 
     , R152F = 303, R152B = 304  & 
     , R153F = 305, R153B = 306  & 
     , R154F = 307, R154B = 308  & 
     , R155F = 309, R155B = 310  & 
     , R156F = 311, R156B = 312  & 
     , R157F = 313, R157B = 314  & 
     , R158F = 315, R158B = 316  & 
     , R159F = 317, R159B = 318  & 
     , R160F = 319, R160B = 320  & 
     , R161F = 321, R161B = 322  & 
     , R162F = 323, R162B = 324  & 
     , R163F = 325, R163B = 326  & 
     , R164F = 327, R164B = 328  & 
     , R165F = 329, R165B = 330 
      integer, parameter ::   R166F = 331, R166B = 332  & 
     , R167F = 333, R167B = 334  & 
     , R168F = 335, R168B = 336  & 
     , R169F = 337, R169B = 338  & 
     , R170F = 339, R170B = 340  & 
     , R171F = 341, R171B = 342  & 
     , R172F = 343, R172B = 344  & 
     , R173F = 345, R173B = 346  & 
     , R174F = 347, R174B = 348  & 
     , R175F = 349, R175B = 350 & 
     , REND = 351 
     integer, parameter ::   MM1 = 1, MM2 = 2 & 
     , MM3 = 3, MM4 = 4 & 
     , MM5 = 5, MM6 = 6 & 
     , MM7 = 7, MM8 = 8 & 
     , MM9 = 9, MM10 = 10 & 
     , MM11 = 11, MM12 = 12 & 
     , MM0 = 13, MM13 = 14 & 
     , MM14 = 15, MM15 = 16 & 
     , MM16 = 17, MM17 = 18 & 
     , MM18 = 19, MM19 = 20 & 
     , MEND = 21 

integer, parameter :: DP=kind(1.0d0)


	real(DP) :: A(REND),N(REND),E(REND)


	integer :: IH

      DATA (A(IH),IH=1,30)  /  3.5470000000D+12, 7.0530202531D+09 & 
     , 5.0800000000D+01, 2.9795859804D+01 & 
     , 2.1600000000D+05, 2.2182689717D+06 & 
     , 2.9700000000D+03, 1.6962435366D+02 & 
     , 4.5770000000D+16, 1.9619085286D+14 & 
     , 6.1650000000D+09, 4.2424184779D+14 & 
     , 4.7140000000D+12, 6.4503598776D+14 & 
     , 3.8000000000D+16, 9.1042926083D+19 & 
     , 6.3660000000D+17, 9.0293444696D+20 & 
     , 1.6600000000D+10, 2.7303646846D+09 & 
     , 7.0790000000D+10, 1.3579714350D+07 & 
     , 3.2500000000D+10, 3.1353652512D+09 & 
     , 2.8900000000D+10, 4.8816975193D+10 & 
     , 4.2000000000D+11, 6.3574616360D+13 & 
     , 1.3000000000D+08, 1.9677857445D+10 /
      DATA (A(IH),IH=31,60)  /  1.2020000000D+17, 1.0739872977D+08 & 
     , 2.4100000000D+10, 5.1591045676D+04 & 
     , 4.8200000000D+10, 5.2375134409D+07 & 
     , 9.5500000000D+03, 6.0865850326D+00 & 
     , 1.0000000000D+09, 1.1159341435D+07 & 
     , 5.8000000000D+11, 6.4724180323D+09 & 
     , 1.5500000000D+21, 1.1858122417D+30 & 
     , 2.5300000000D+09, 2.8127059229D+13 & 
     , 3.0100000000D+10, 3.2283088892D+13 & 
     , 2.2290000000D+02, 1.2462349669D+09 & 
     , 4.7480000000D+08, 9.5208413773D+04 & 
     , 7.5800000000D+09, 2.1558741858D+09 & 
     , 7.2300000000D+10, 3.3822450822D+09 & 
     , 3.0200000000D+10, 8.2864006839D+08 & 
     , 3.0200000000D+10, 1.4508889496D+10 /
      DATA (A(IH),IH=61,90)  /  3.0000000000D+10, 4.6022489291D+15 & 
     , 3.0000000000D+10, 6.4520015108D+09 & 
     , 3.0000000000D+09, 2.8141850446D+04 & 
     , 2.6500000000D+10, 4.2149530019D+11 & 
     , 3.0000000000D+10, 9.8375045521D+11 & 
     , 3.3000000000D+36, 2.0179707531D+31 & 
     , 3.1000000000D+42, 8.8680758342D+35 & 
     , 5.7400000000D+04, 8.1887003172D+01 & 
     , 1.8100000000D+10, 1.5145160141D+07 & 
     , 3.4300000000D+06, 5.0252498407D+04 & 
     , 1.2300000000D+03, 1.0668309991D+01 & 
     , 4.1100000000D+01, 5.3959360590D+01 & 
     , 3.6360000000D-09, 1.7636282446D-09 & 
     , 8.4300000000D+10, 1.3513165132D+13 & 
     , 1.9900000000D+15, 3.5389559985D+17 /
      DATA (A(IH),IH=91,120)  /  3.7400000000D+08, 1.1921064944D+08 & 
     , 2.4100000000D+07, 4.1346990936D+08 & 
     , 8.0540000000D+28, 1.2503013086D+40 & 
     , 2.4770000000D+30, 1.9647543084D+35 & 
     , 5.4700000000D+04, 1.6088191997D+02 & 
     , 3.1500000000D+09, 5.4340380842D+06 & 
     , 5.7200000000D+03, 1.7277320766D+02 & 
     , 3.1600000000D+09, 1.7671757301D+11 & 
     , 1.8100000000D+08, 4.8991459399D+08 & 
     , 1.0000000000D+11, 1.6133092170D+07 & 
     , 6.0000000000D+09, 2.2582448200D+08 & 
     , 9.6350000000D+10, 1.3268877722D+07 & 
     , 4.2000000000D+10, 9.2717405987D+08 & 
     , 2.4000000000D+10, 9.2766563237D+09 & 
     , 2.4100000000D+11, 5.5147298808D+10 /
      DATA (A(IH),IH=121,150)  /  1.5100000000D+12, 3.4552871867D+11 & 
     , 1.2000000000D+10, 4.1564533076D+11 & 
     , 1.5000000000D+10, 3.9573817470D+11 & 
     , 8.3000000000D+14, 1.0871891022D+10 & 
     , 3.2000000000D+10, 3.5780123422D+05 & 
     , 6.0000000000D+09, 1.0754063973D+07 & 
     , 1.8000000000D+10, 5.6488769410D+08 & 
     , 9.0330000000D+10, 1.6782188533D+09 & 
     , 2.2000000000D+07, 4.0873258908D+05 & 
     , 3.0000000000D+08, 8.4366940079D+08 & 
     , 1.6000000000D+10, 1.0002342066D+12 & 
     , 4.9900000000D+09, 7.1770697790D+14 & 
     , 2.4600000000D+03, 3.0002823388D+03 & 
     , 1.6000000000D+10, 1.4227561628D+10 & 
     , 5.6000000000D+04, 1.3868879542D+03 /
      DATA (A(IH),IH=151,180)  /  2.5010000000D+10, 8.4953948027D+08 & 
     , 4.0000000000D+10, 3.6126044532D+17 & 
     , 1.2000000000D+10, 7.9017837861D+16 & 
     , 1.6000000000D+10, 6.0768947333D+03 & 
     , 1.1500000000D+05, 2.4856705124D+01 & 
     , 8.9800000000D+04, 1.1384507997D+01 & 
     , 3.5400000000D+03, 7.8579531528D+00 & 
     , 4.0000000000D+10, 5.2564573058D+07 & 
     , 2.9400000000D+08, 5.8481019909D+07 & 
     , 6.1400000000D+03, 4.5122606321D+02 & 
     , 1.9900000000D+38, 2.1478775258D+44 & 
     , 2.0000000000D+09, 3.0285635419D+08 & 
     , 1.3200000000D+11, 1.4711495303D+08 & 
     , 2.0000000000D+07, 1.8412981636D+07 & 
     , 1.4000000000D+09, 9.8081931579D+11 /
      DATA (A(IH),IH=181,210)  /  1.2000000000D+11, 2.5971815126D+13 & 
     , 8.0200000000D+10, 1.0134383975D+14 & 
     , 7.0000000000D+50, 1.3620567873D+45 & 
     , 1.2000000000D+39, 1.8487437087D+42 & 
     , 1.4000000000D+27, 7.6776073818D+32 & 
     , 1.3250000000D+03, 5.6366383607D-01 & 
     , 1.8000000000D+03, 7.8638867801D+00 & 
     , 2.2700000000D+02, 3.2832975054D+01 & 
     , 1.9200000000D+04, 2.0159534780D-01 & 
     , 5.0000000000D+09, 5.4793052566D+10 & 
     , 1.5100000000D+04, 3.7676751508D+00 & 
     , 4.2150000000D+10, 1.0901576399D+08 & 
     , 9.6400000000D+10, 1.0286605207D+11 & 
     , 1.2100000000D+07, 3.0907166297D+07 & 
     , 3.9000000000D+08, 1.4149455616D+11 /
      DATA (A(IH),IH=211,240)  /  9.6000000000D+08, 2.4080304685D+12 & 
     , 4.5800000000D+13, 6.1431582523D+11 & 
     , 1.3370000000D+03, 8.6738893086D+03 & 
     , 3.8000000000D+37, 8.3078895790D+39 & 
     , 4.0800000000D+03, 1.0646326490D-02 & 
     , 4.8300000000D-07, 8.9104911373D-10 & 
     , 3.2000000000D+24, 3.0957026570D+29 & 
     , 8.0000000000D+10, 7.5862942833D+12 & 
     , 2.0000000000D+10, 2.2665974682D+15 & 
     , 5.0000000000D+02, 2.0733714238D+05 & 
     , 1.3200000000D+10, 2.4890153347D+09 & 
     , 2.0000000000D+10, 4.3480359040D+11 & 
     , 3.2000000000D+10, 5.4402246782D+16 & 
     , 9.0000000000D+09, 6.5618452463D+09 & 
     , 3.0000000000D+10, 2.1872817488D+10 /
      DATA (A(IH),IH=241,270)  /  9.0000000000D+09, 6.5618452463D+09 & 
     , 7.0000000000D+09, 5.1036574138D+09 & 
     , 1.5000000000D+10, 4.8515577882D+10 & 
     , 1.5000000000D+10, 1.0370851892D+12 & 
     , 3.0000000000D+10, 2.4788436370D+15 & 
     , 7.0000000000D+10, 2.1163554877D+13 & 
     , 2.8000000000D+10, 7.7189773425D+05 & 
     , 1.2000000000D+10, 7.9258453041D+08 & 
     , 1.4000000000D+10, 2.0690264836D+08 & 
     , 7.0000000000D+15, 3.7753643217D+04 & 
     , 1.6990000000D+42, 2.1946663920D+31 & 
     , 6.7100000000D+03, 2.5724684790D+02 & 
     , 2.9700000000D+04, 1.1087242276D+02 & 
     , 2.6800000000D-02, 3.4015900172D-02 & 
     , 1.8550000000D-06, 4.0616566920D-09 /
      DATA (A(IH),IH=271,300)  /  2.0000000000D+10, 6.8709845478D+10 & 
     , 4.1000000000D+10, 9.3054713257D+08 & 
     , 1.2000000000D+13, 1.2688738503D+02 & 
     , 2.4100000000D+10, 1.9727820584D+10 & 
     , 5.4900000000D+00, 2.0980138029D+00 & 
     , 9.0000000000D+09, 4.9424652594D+11 & 
     , 1.7450000000D+16, 4.2975913579D+12 & 
     , 1.0000000000D+10, 2.6565778440D+07 & 
     , 2.3400000000D+04, 1.0500535004D+02 & 
     , 1.2200000000D+09, 4.9058790251D+08 & 
     , 2.3500000000D+02, 6.0227548826D-02 & 
     , 4.5500000000D+03, 1.9881379390D+00 & 
     , 7.5500000000D-04, 1.1216623790D-04 & 
     , 7.4510000000D+12, 6.5430526784D+01 & 
     , 1.5140000000D+12, 8.3113859279D+02 /
      DATA (A(IH),IH=301,330)  /  2.0000000000D+09, 7.7901495616D+19 & 
     , 1.5970000000D+20, 2.3471800921D+09 & 
     , 6.8440000000D+19, 1.2073443805D+09 & 
     , 9.7220000000D+15, 3.2115786960D+04 & 
     , 5.0000000000D+07, 1.7465831985D+07 & 
     , 6.0000000000D+10, 2.7310450662D+12 & 
     , 1.5000000000D+13, 2.8515336158D-11 & 
     , 7.0000000000D+08, 3.0813052680D+19 & 
     , 4.0000000000D+10, 3.4129896457D+00 & 
     , 3.0000000000D+16, 7.8411183626D+04 & 
     , 1.0000000000D+11, 1.0908995697D+09 & 
     , 2.1770000000D+16, 1.1879914769D+05 & 
     , 5.3110000000D+15, 1.0900751837D+06 & 
     , 1.0000000000D+14, 1.4245687646D+11 & 
     , 4.5000000000D+12, 4.1462884067D+21 /
      DATA (A(IH),IH=331,350)  /  2.3000000000D+10, 5.1478464291D-02 & 
     , 1.5000000000D+13, 1.8277568334D+07 & 
     , 4.5930000000D+18, 2.1397693475D+07 & 
     , 2.6200000000D+03, 3.2786035279D-02 & 
     , 1.8500000000D+04, 4.1406590841D-08 & 
     , 4.2400000000D+03, 5.1664593155D-03 & 
     , 6.0300000000D+10, 1.3141792487D-02 & 
     , 3.9000000000D-10, 2.8898944836D-20 & 
     , 1.0000000000D+09, 2.0056686253D-01 & 
     , 1.7700000000D+15, 2.2625738811D+02 & 
     /
      DATA (N(IH),IH=1,30)  /  -0.406, 0.0289576 & 
     , 2.67, 2.63497 & 
     , 1.51, 1.41638 & 
     , 2.02, 2.07859 & 
     , -1.4, -1.7488 & 
     , -0.5, -0.621192 & 
     , -1, -0.686234 & 
     , -2, -1.74482 & 
     , -1.72, -1.73291 & 
     , 0, 0.361707 & 
     , 0, 0.761631 & 
     , 0, 0.326673 & 
     , 0, 0.268083 & 
     , 0, -0.389273 & 
     , 0, -0.389273 /
      DATA (N(IH),IH=31,60)  /  0, 1.16381 & 
     , 0, 1.41899 & 
     , 0, 0.75098 & 
     , 2, 2.71595 & 
     , 0, 0.657356 & 
     , 0, 0.657356 & 
     , -2.79, -3.75543 & 
     , 0, -0.844239 & 
     , 0, -0.517566 & 
     , 1.89, 0.610803 & 
     , 0.659, 0.891462 & 
     , 0, 0.219554 & 
     , 0, 0.581262 & 
     , 0, 0.546227 & 
     , 0, 0.487638 /
      DATA (N(IH),IH=61,90)  /  0, -0.73297 & 
     , 0, -0.285105 & 
     , 0, 0.813723 & 
     , 0, 0.291083 & 
     , 0, 0.111904 & 
     , -6.3, -6.17944 & 
     , -8, -7.29818 & 
     , 1.9, 2.36936 & 
     , 0, 0.434323 & 
     , 1.18, 1.55573 & 
     , 3, 3.10765 & 
     , 2.5, 2.21838 & 
     , 5.42, 5.59918 & 
     , 0, -0.304279 & 
     , -1.57, -2.05404 /
      DATA (N(IH),IH=91,120)  /  0, 0.130678 & 
     , 0.76, 0.602631 & 
     , -3.75, -5.10728 & 
     , -4.76, -4.70138 & 
     , 1.97, 2.26018 & 
     , 0.5, 0.755145 & 
     , 1.96, 2.15656 & 
     , 0, 0.0715283 & 
     , 0, -0.460801 & 
     , 0, 0.3196 & 
     , 0, 0.6684 & 
     , 0, 0.937644 & 
     , 0, 0.633365 & 
     , 0, 0.574776 & 
     , 0, 0.306692 /
      DATA (N(IH),IH=121,150)  /  -1, -0.693308 & 
     , 0, -0.0825805 & 
     , 0, 0.199042 & 
     , -1.2, -0.899045 & 
     , 0, 0.919 & 
     , 0, 0.61472 & 
     , 0, 0.556131 & 
     , 0, 0.288048 & 
     , 0, 0.288048 & 
     , 0, -0.101225 & 
     , 0, -0.360197 & 
     , 0.1, -1.06665 & 
     , 2, 1.77911 & 
     , 0, -0.283175 & 
     , 1.6, 2.01744 /
      DATA (N(IH),IH=151,180)  /  0, 0.47973 & 
     , 0, -1.23115 & 
     , 0, -1.29344 & 
     , 0, 1.39873 & 
     , 1.9, 2.43943 & 
     , 1.92, 2.4244 & 
     , 2.12, 2.56581 & 
     , 0, 0.177723 & 
     , 0, -0.21155 & 
     , 1.74, 1.98925 & 
     , -7.08, -7.27063 & 
     , 0, 0.44657 & 
     , 0, 0.862375 & 
     , 0, 0.084863 & 
     , 0, -0.0928598 /
      DATA (N(IH),IH=181,210)  /  0, 0.0418316 & 
     , 0, -0.536729 & 
     , -9.31, -8.99745 & 
     , -7.62, -7.71777 & 
     , -3.86, -3.99563 & 
     , 2.53, 3.01443 & 
     , 2, 2.39081 & 
     , 2, 2.19425 & 
     , 1.83, 2.71516 & 
     , 0, 0.0832929 & 
     , 1.91, 2.3594 & 
     , 0, 0.122726 & 
     , 0, 0.176917 & 
     , 0, 0.266547 & 
     , 0, -0.113262 /
      DATA (N(IH),IH=211,240)  /  0, -0.307517 & 
     , -1.39, -0.823557 & 
     , 1.61, 1.42521 & 
     , -7.27, -7.09812 & 
     , 2, 3.31614 & 
     , 4, 4.84011 & 
     , -3.14, -3.30227 & 
     , 0, -0.345989 & 
     , 0, -0.780313 & 
     , 2, 1.48893 & 
     , 0, 0.0889686 & 
     , 0, -0.018682 & 
     , 0, -1.08087 & 
     , 0, -0.0622863 & 
     , 0, -0.0622863 /
      DATA (N(IH),IH=241,270)  /  0, -0.0622863 & 
     , 0, -0.0622863 & 
     , 0, 0.172986 & 
     , 0, -0.408275 & 
     , 0, -0.842599 & 
     , 0, -0.573354 & 
     , 0, 0.259144 & 
     , 0, 0.51432 & 
     , 0, 0.436598 & 
     , 0, 1.51966 & 
     , -7.954, -6.47478 & 
     , 2, 2.10394 & 
     , 2, 2.19757 & 
     , 3.778, 3.68539 & 
     , 5.29, 5.45253 /
      DATA (N(IH),IH=271,300)  /  0, -0.553412 & 
     , 0, -0.16414 & 
     , 0, 1.9314 & 
     , 0, 0.452187 & 
     , 2.8, 3.07179 & 
     , 0, -0.286081 & 
     , -0.66, -0.727194 & 
     , 0, -1.12449 & 
     , 1.61, 0.753594 & 
     , 0, -1.51376 & 
     , 2.5, 1.70218 & 
     , 2, 1.23722 & 
     , 3.46, 2.40704 & 
     , -1.76, 1.83191 & 
     , -1.78, 1.45172 /
      DATA (N(IH),IH=301,330)  /  0, -1.66077 & 
     , -4.5, -2.52515 & 
     , -4.5, -2.46046 & 
     , -1.1, 0.960117 & 
     , 0, -0.0801019 & 
     , 0, -0.87687 & 
     , 0, 4.59973 & 
     , 0, -1.67723 & 
     , 0, 1.45568 & 
     , 0, 2.03667 & 
     , 0, 0.529885 & 
     , -2.69, -0.563231 & 
     , -2.61, -0.852332 & 
     , 0, -0.102835 & 
     , -1.11, -2.3397 /
      DATA (N(IH),IH=331,350)  /  0, 1.94073 & 
     , 0, 0.755153 & 
     , -0.46, 0.993088 & 
     , 2.06, 2.72153 & 
     , 1.51, 3.45073 & 
     , 2.1, 2.85515 & 
     , -0.35, 1.68435 & 
     , 5.8, 7.54417 & 
     , 0, 1.28337 & 
     , -1.9, 0.0993154 & 
     /
      DATA (E(IH),IH=1,30)  /  69450000, -1111131.16 & 
      , 26317000, 20550909.81 & 
      , 14351000, 77059304.43 & 
      , 56066000, -12408394.61 & 
      , 436726000, 3652636.474 & 
      , 0, 497868404.5 & 
      , 0, 427307273.3 & 
      , 0, 495781668 & 
      , 2196000, 206859348.7 & 
      , 3443000, 231853014.8 & 
      , 1234000, 153316793.4 & 
      , 0, 222643924.6 & 
      , -2079000, 289039319.2 & 
      , 50133000, 212208672.1 & 
      , -6817000, 155258672.1 /
      DATA (E(IH),IH=31,60)  /  190372000, -24284227.4 & 
      , 16610000, 297735440.6 & 
      , 33263000, 99597342.69 & 
      , 16610000, 77178252.5 & 
      , 0, 129042647.1 & 
      , 39986000, 169028647.1 & 
      , 17535000, 552312028.5 & 
      , 199577000, 236485624 & 
      , 96232000, 355784548.6 & 
      , -4848000, 102621755.2 & 
      , 62233000, -1964560.922 & 
      , 1715000, 142180787.8 & 
      , 0, 368875802.6 & 
      , 0, 363109712.4 & 
      , 0, 431584107 /
      DATA (E(IH),IH=61,90)  /  0, 470579467.6 & 
      , 0, 195354987.7 & 
      , 0, 304678241.7 & 
      , 0, 376812980.9 & 
      , 0, 314381994.1 & 
      , 417982000, 39402444.96 & 
      , 407982000, 398278247.6 & 
      , 11500000, 65993808.48 & 
      , 12887000, 61614718.3 & 
      , -1870000, 115332112.9 & 
      , 217568000, 43651793.69 & 
      , 42719000, 30878465.8 & 
      , 4176000, 66606986.74 & 
      , 0, 296214081.9 & 
      , 122298000, 7423963.149 /
      DATA (E(IH),IH=91,120)  /  61254000, 286906950.7 & 
      , -9728000, 98041887.75 & 
      , 4107000, 388336266 & 
      , 10209000, 451219541.8 & 
      , 46903000, 38965821.74 & 
      , 43053000, 29349731.56 & 
      , 11042000, 65813126.17 & 
      , 0, 236347193 & 
      , 77739000, 3467479.058 & 
      , 105018000, -13210812.05 & 
      , 0, 314844551.5 & 
      , 0, 12864379.42 & 
      , 0, 309078461.3 & 
      , 0, 377552855.9 & 
      , 20991000, 107425536.7 /
      DATA (E(IH),IH=121,150)  /  0, 86434536.68 & 
      , 0, 248510208.8 & 
      , 0, 260350743 & 
      , 64852000, -21928285.78 & 
      , 0, 44312905.69 & 
      , 0, 340526987.6 & 
      , 0, 409001382.2 & 
      , 50124000, 168007063 & 
      , 7314000, 125197063 & 
      , 0, 279958735.1 & 
      , 49371000, 201153660.8 & 
      , 44350000, 7072459.094 & 
      , 34602000, 57823753.84 & 
      , -2385000, 58139741.82 & 
      , 22677000, 54226372.33 /
      DATA (E(IH),IH=151,180)  /  0, -5753615.644 & 
      , 0, 275776402 & 
      , -2385000, 310694390 & 
      , 0, 38559290.05 & 
      , 31506000, 43072556.61 & 
      , 23807000, 29607466.43 & 
      , 3640000, 77914861.04 & 
      , 212966000, -3877458.177 & 
      , 62509000, 7741213.928 & 
      , 43723000, 63226734.87 & 
      , 27970000, 449476806.9 & 
      , 0, 281895010.8 & 
      , 0, 333491622.8 & 
      , 0, 53484996.04 & 
      , 0, 270328454.2 /
      DATA (E(IH),IH=181,210)  /  0, 357309246 & 
      , 0, 317349590.2 & 
      , 417814000, 237600623.6 & 
      , 29162000, 180340352.7 & 
      , 13891000, 480393379.5 & 
      , 51212000, 17782984.04 & 
      , 10460000, 39739288.47 & 
      , 38493000, 13001162.3 & 
      , 920000, 107010420.4 & 
      , 0, 348997307.5 & 
      , 15648000, -23547106.15 & 
      , 240998000, -20841030.75 & 
      , 0, 286289003.1 & 
      , -2494000, 97269358.64 & 
      , 0, 294226181.4 /
      DATA (E(IH),IH=211,240)  /  0, 319718019.1 & 
      , 4247000, 375185477.3 & 
      , -1607000, 56271988.32 & 
      , 30208000, 176992360.4 & 
      , 7950000, 198897303.8 & 
      , -8368000, 219504326.1 & 
      , 5146000, 469378295.6 & 
      , 0, 381866822.5 & 
      , 0, 333139104.2 & 
      , 30250000, 61408932.1 & 
      , 6276000, 317581691.3 & 
      , 0, 485221897.6 & 
      , 0, 559795321.3 & 
      , 2510000, 39812987.98 & 
      , 0, 37302987.98 /
      DATA (E(IH),IH=241,270)  /  0, 37302987.98 & 
      , 0, 37302987.98 & 
      , 0, 788045613 & 
      , 0, 419169810.4 & 
      , 0, 370442092.1 & 
      , 0, 68461920.07 & 
      , 0, 284411118.3 & 
      , 0, 780192786.3 & 
      , 0, 262972337 & 
      , 341724000, -20713522.45 & 
      , 384119000, 33101635.18 & 
      , -2635000, 86734824.77 & 
      , 16877000, 43538520.34 & 
      , 40297000, 74895698.6 & 
      , -456000, 20439430.15 /
      DATA (E(IH),IH=271,300)  /  69036000, 29363177.65 & 
      , 187903000, -13845494.45 & 
      , 107738000, 76352192.59 & 
      , 0, 319631557.4 & 
      , 24527000, 52359288.14 & 
      , 0, 124405355.8 & 
      , 49036000, 40420853 & 
      , 207945000, 3300998.939 & 
      , -146000, 86328318.16 & 
      , 71128000, 28559671.04 & 
      , 9330000, 27329923.54 & 
      , 20920000, 44686013.73 & 
      , 22933000, 54636191.99 & 
      , 71756000, -1119894.578 & 
      , 57823000, 136729766.3 /
      DATA (E(IH),IH=301,330)  /  0, 152376875.2 & 
      , 0, -3362483.601 & 
      , 0, 422480241.3 & 
      , 86358000, 38336724.59 & 
      , 2092000, 198140201.7 & 
      , 89956000, 46223662.18 & 
      , 85772000, 171394605.9 & 
      , 0, 152474430.4 & 
      , 77404000, 248443163.9 & 
      , 167360000, -19479411.25 & 
      , 58576000, 85941547.44 & 
      , 71965000, 13891644.72 & 
      , 87069000, 148319468.8 & 
      , 62342000, 41360025.96 & 
      , 0, 106374743.1 /
      DATA (E(IH),IH=331,350)  /  209200000, 176811783 & 
      , 238488000, 250861233.7 & 
      , 453127000, -10845324.08 & 
      , 3833000, 78914538.1 & 
      , -4025000, -36413217.05 & 
      , 20368000, 32741233.67 & 
      , 12502000, -82594521.48 & 
      , 9205000, -77954343.22 & 
      , 49873000, -111557864.2 & 
      , 12447000, -88415611.66 & 
     /