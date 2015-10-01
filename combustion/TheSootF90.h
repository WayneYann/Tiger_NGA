!-------------------------------------------------------------
! ======= TheSootF90.h =======
!-------------------------------------------------------------

     integer, parameter ::   SN2 = 1, SH = 2 & 
     , SO2 = 3, SO = 4 & 
     , SOH = 5, SH2 = 6 & 
     , SH2O = 7, SCO2 = 8 & 
     , SHO2 = 9, SH2O2 = 10 & 
     , SCO = 11, SHCO = 12 & 
     , SC = 13, SCH = 14 & 
     , STXCH2 = 15, SCH3 = 16 & 
     , SCH2O = 17, SHCCO = 18 & 
     , SC2H = 19, SCH2CO = 20 & 
     , SC2H2 = 21, SSXCH2 = 22 & 
     , SAR = 23, SCH3OH = 24 & 
     , SCH2OH = 25, SCH3O = 26 & 
     , SCH4 = 27, SCH3O2 = 28 & 
     , SC2H3 = 29, SC2H4 = 30 
     integer, parameter ::   SC2H5 = 31, SHCCOH = 32 & 
     , SCH2CHO = 33, SCH3CHO = 34 & 
     , SH2C2 = 35, SC2H5O = 36 & 
     , SNXC3H7 = 37, SC2H6 = 38 & 
     , SC3H8 = 39, SC3H6 = 40 & 
     , SC3H3 = 41, SPXC3H4 = 42 & 
     , SAXC3H4 = 43, SSXC3H5 = 44 & 
     , SNXC4H3 = 45, SC2H3CHO = 46 & 
     , SAXC3H5 = 47, SC2O = 48 & 
     , SC4H4 = 49, SC3H2 = 50 & 
     , SC3H2O = 51, SC4H2 = 52 & 
     , SIXC4H3 = 53, STXC3H5 = 54 & 
     , SC3H5O = 55, SC4H = 56 & 
     , SC8H2 = 57, SC6H2 = 58 & 
     , SC4H6 = 59, SNXC4H5 = 60 
     integer, parameter ::   SIXC4H5 = 61, SA1XC6H6 = 62 & 
     , SNXC7H16 = 63, SC5H11 = 64 & 
     , SPXC4H9 = 65, SC7H15 = 66 & 
     , SPXC4H8 = 67, SC5H10 = 68 & 
     , SC7H14 = 69, SC7H15O = 70 & 
     , SC3H7CHO = 71, SC4H7 = 72 & 
     , SC7H13 = 73, SC5H9 = 74 & 
     , SC4H7O = 75, SNXC3H7O = 76 & 
     , SIXC8H18 = 77, SYXC7H15 = 78 & 
     , SIXC4H8 = 79, SIXC3H7 = 80 & 
     , STXC4H9 = 81, SCXC8H17 = 82 & 
     , SYXC7H14 = 83, SDXC8H17O = 84 & 
     , SCH3COCH3 = 85, SIXC4H7 = 86 & 
     , SXXC7H13 = 87, SIXC3H5CH = 88 & 
     , STXC4H9O = 89, SIXC4H7O = 90 
     integer, parameter ::   SC5H4CH2 = 91, SA1XXC6H5 = 92 & 
     , SA1C2H2XC = 93, SA1C2H3XC = 94 & 
     , SA1C2HXC8 = 95, SA1C2HYXC = 96 & 
     , SA1C2H3YX = 97, SA2XXC10H = 98 & 
     , SA2XC10H8 = 99, SA2YXC10H = 100 & 
     , SA2C2H2AX = 101, SA2C2H2BX = 102 & 
     , SA2C2HAXC = 103, SA2C2HBXC = 104 & 
     , SA2C2HAYX = 105, SA2C2HBYX = 106 & 
     , SA2R5XC12 = 107, SA2R5XXC1 = 108 & 
     , SA2R5C2H2 = 109, SA2R5C2HX = 110 & 
     , SA2R5C2HY = 111, SP2XC12H1 = 112 & 
     , SP2XXC12H = 113, SA3XXC14H = 114 & 
     , SA3XC14H1 = 115, SA3YXC14H = 116 & 
     , SA3R5XXC1 = 117, SA3R5XC16 = 118 & 
     , SA4XC16H1 = 119, SA4XXC16H = 120 
     integer, parameter ::   SA4R5XC18 = 121, SFLTNXC16 = 122 & 
     , SC5H6 = 123, SC5H5 = 124 & 
     , STXC5H5O = 125, SC5H4O = 126 & 
     , SSXC5H5O = 127, SC9H8 = 128 & 
     , SC9H7 = 129, SA1CH2XC7 = 130 & 
     , SC9H6O = 131, SOXC6H4 = 132 & 
     , SA1CH3XC7 = 133, SA1OHXC6H = 134 & 
     , SHOA1CH3X = 135, SOA1CH3XC = 136 & 
     , SA1CH2OXC = 137, SA1CH2OHX = 138 & 
     , SA1CHOXC7 = 139, SA1OXC6H5 = 140 & 
     , SA1CH3YXC = 141, SA1C2H4XC = 142 & 
     , SA1C2H5XC = 143, SC8H9O2 = 144 & 
     , SC8H8OOH = 145, SOC8H7OOH = 146 & 
     , SA1CH3CH3 = 147, SA1CH3CH2 = 148 & 
     , SA1CH3CHO = 149, SA2CH3XC1 = 150 
     integer, parameter ::   SA1CHOCH2 = 151, SA1CHOCHO = 152 & 
     , SA2OHXC10 = 153, SA2CH2XC1 = 154 & 
     , SA2CH2OXC = 155, SA2CHOXC1 = 156 & 
     , SA2OXC10H = 157, SOC6H4O = 158 & 
     , SEND = 159 
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
     , RG01F = 77, RG01B = 78  & 
     , RG02F = 79, RG02B = 80  & 
     , RG03F = 81, RG03B = 82  & 
     , RG04F = 83, RG04B = 84  & 
     , RG05F = 85, RG05B = 86  & 
     , RG06F = 87, RG06B = 88  & 
     , RG07F = 89, RG07B = 90 
      integer, parameter ::   RG08F = 91, RG08B = 92  & 
     , RG09F = 93, RG09B = 94  & 
     , RG10F = 95, RG10B = 96  & 
     , RG11F = 97, RG11B = 98  & 
     , RG12F = 99, RG12B = 100  & 
     , RG13F = 101, RG13B = 102  & 
     , RG14F = 103, RG14B = 104  & 
     , RG15F = 105, RG15B = 106  & 
     , RG16F = 107, RG16B = 108  & 
     , RG17F = 109, RG17B = 110  & 
     , RG18F = 111, RG18B = 112  & 
     , RG19 = 113, RG20F = 114  & 
     , RG20B = 115, RG21 = 116  & 
     , RG22F = 117, RG22B = 118  & 
     , RG23F = 119, RG23B = 120 
      integer, parameter ::   RG24F = 121, RG24B = 122  & 
     , RG25F = 123, RG25B = 124  & 
     , RG26F = 125, RG26B = 126  & 
     , RG27 = 127, RG28F = 128  & 
     , RG28B = 129, RG29F = 130  & 
     , RG29B = 131, RG30F = 132  & 
     , RG30B = 133, RG31F = 134  & 
     , RG31B = 135, RG32F = 136  & 
     , RG32B = 137, RG33F = 138  & 
     , RG33B = 139, RG34F = 140  & 
     , RG34B = 141, RG35F = 142  & 
     , RG35B = 143, RG36F = 144  & 
     , RG36B = 145, RG37F = 146  & 
     , RG37B = 147, RG38F = 148  & 
     , RG38B = 149, RG39 = 150 
      integer, parameter ::   RG40F = 151, RG40B = 152  & 
     , RG41F = 153, RG41B = 154  & 
     , RG42F = 155, RG42B = 156  & 
     , RG43F = 157, RG43B = 158  & 
     , RG44F = 159, RG44B = 160  & 
     , RG45F = 161, RG45B = 162  & 
     , RG46F = 163, RG46B = 164  & 
     , RG47F = 165, RG47B = 166  & 
     , RG48F = 167, RG48B = 168  & 
     , RG49F = 169, RG49B = 170  & 
     , RG50F = 171, RG50B = 172  & 
     , RG51F = 173, RG51B = 174  & 
     , RG52F = 175, RG52B = 176  & 
     , RG53 = 177, RG54F = 178  & 
     , RG54B = 179, RG55F = 180 
      integer, parameter ::   RG55B = 181, RG56 = 182  & 
     , RG57F = 183, RG57B = 184  & 
     , RG58F = 185, RG58B = 186  & 
     , RG59F = 187, RG59B = 188  & 
     , RG60F = 189, RG60B = 190  & 
     , RG61F = 191, RG61B = 192  & 
     , RG62 = 193, RG63 = 194  & 
     , RG64 = 195, RG65F = 196  & 
     , RG65B = 197, RG66F = 198  & 
     , RG66B = 199, RG67F = 200  & 
     , RG67B = 201, RG68F = 202  & 
     , RG68B = 203, RG69F = 204  & 
     , RG69B = 205, RG70F = 206  & 
     , RG70B = 207, RG71F = 208  & 
     , RG71B = 209, RG72F = 210 
      integer, parameter ::   RG72B = 211, RG73F = 212  & 
     , RG73B = 213, RG74F = 214  & 
     , RG74B = 215, RG75F = 216  & 
     , RG75B = 217, RG76F = 218  & 
     , RG76B = 219, RG77F = 220  & 
     , RG77B = 221, RG78F = 222  & 
     , RG78B = 223, RG79F = 224  & 
     , RG79B = 225, RG80F = 226  & 
     , RG80B = 227, RG81F = 228  & 
     , RG81B = 229, RG82F = 230  & 
     , RG82B = 231, RG83F = 232  & 
     , RG83B = 233, RG84F = 234  & 
     , RG84B = 235, RG85F = 236  & 
     , RG85B = 237, RG86F = 238  & 
     , RG86B = 239, RG87F = 240 
      integer, parameter ::   RG87B = 241, RG88F = 242  & 
     , RG88B = 243, RG89F = 244  & 
     , RG89B = 245, RG90F = 246  & 
     , RG90B = 247, RG91F = 248  & 
     , RG91B = 249, RG92F = 250  & 
     , RG92B = 251, RG93F = 252  & 
     , RG93B = 253, RG94F = 254  & 
     , RG94B = 255, RG95F = 256  & 
     , RG95B = 257, RG96F = 258  & 
     , RG96B = 259, RG97F = 260  & 
     , RG97B = 261, RG98F = 262  & 
     , RG98B = 263, RG99F = 264  & 
     , RG99B = 265, RG100F = 266  & 
     , RG100B = 267, RG101F = 268  & 
     , RG101B = 269, RG102F = 270 
      integer, parameter ::   RG102B = 271, RG103F = 272  & 
     , RG103B = 273, RG104F = 274  & 
     , RG104B = 275, RG105F = 276  & 
     , RG105B = 277, RG106F = 278  & 
     , RG106B = 279, RG107F = 280  & 
     , RG107B = 281, RG108F = 282  & 
     , RG108B = 283, RG109F = 284  & 
     , RG109B = 285, RG110F = 286  & 
     , RG110B = 287, RG111F = 288  & 
     , RG111B = 289, RG112F = 290  & 
     , RG112B = 291, RG113F = 292  & 
     , RG113B = 293, RG114F = 294  & 
     , RG114B = 295, RG115F = 296  & 
     , RG115B = 297, RG116F = 298  & 
     , RG116B = 299, RG117F = 300 
      integer, parameter ::   RG117B = 301, RG118F = 302  & 
     , RG118B = 303, RG119F = 304  & 
     , RG119B = 305, RG120F = 306  & 
     , RG120B = 307, RG121F = 308  & 
     , RG121B = 309, RG122F = 310  & 
     , RG122B = 311, RG123F = 312  & 
     , RG123B = 313, RG124F = 314  & 
     , RG124B = 315, RG125F = 316  & 
     , RG125B = 317, RG126F = 318  & 
     , RG126B = 319, RG127F = 320  & 
     , RG127B = 321, RG128F = 322  & 
     , RG128B = 323, RG129F = 324  & 
     , RG129B = 325, RG130F = 326  & 
     , RG130B = 327, RG131F = 328  & 
     , RG131B = 329, RG132F = 330 
      integer, parameter ::   RG132B = 331, RG133F = 332  & 
     , RG133B = 333, RG134F = 334  & 
     , RG134B = 335, RG135F = 336  & 
     , RG135B = 337, RG136F = 338  & 
     , RG136B = 339, RG137F = 340  & 
     , RG137B = 341, RG138F = 342  & 
     , RG138B = 343, RG139 = 344  & 
     , RG140 = 345, RG141F = 346  & 
     , RG141B = 347, RG142F = 348  & 
     , RG142B = 349, RG143F = 350  & 
     , RG143B = 351, RG144F = 352  & 
     , RG144B = 353, RG145F = 354  & 
     , RG145B = 355, RG146F = 356  & 
     , RG146B = 357, RG147F = 358  & 
     , RG147B = 359, RG148 = 360 
      integer, parameter ::   RG149 = 361, RG150 = 362  & 
     , RG151 = 363, RG152 = 364  & 
     , RG153 = 365, RG154F = 366  & 
     , RG154B = 367, RG155F = 368  & 
     , RG155B = 369, RG156F = 370  & 
     , RG156B = 371, RG157F = 372  & 
     , RG157B = 373, RG158F = 374  & 
     , RG158B = 375, RG159F = 376  & 
     , RG159B = 377, RG160F = 378  & 
     , RG160B = 379, RG161F = 380  & 
     , RG161B = 381, RG162F = 382  & 
     , RG162B = 383, RG163F = 384  & 
     , RG163B = 385, RG164F = 386  & 
     , RG164B = 387, RG165F = 388  & 
     , RG165B = 389, RG166F = 390 
      integer, parameter ::   RG166B = 391, RG167F = 392  & 
     , RG167B = 393, RG168F = 394  & 
     , RG168B = 395, RG169F = 396  & 
     , RG169B = 397, RG170F = 398  & 
     , RG170B = 399, RG171F = 400  & 
     , RG171B = 401, RG172F = 402  & 
     , RG172B = 403, RG173F = 404  & 
     , RG173B = 405, RG174F = 406  & 
     , RG174B = 407, RG175F = 408  & 
     , RG175B = 409, RG176F = 410  & 
     , RG176B = 411, RG177F = 412  & 
     , RG177B = 413, RG178F = 414  & 
     , RG178B = 415, RG179F = 416  & 
     , RG179B = 417, RG180F = 418  & 
     , RG180B = 419, RG181F = 420 
      integer, parameter ::   RG181B = 421, RG182F = 422  & 
     , RG182B = 423, RG183F = 424  & 
     , RG183B = 425, RG184F = 426  & 
     , RG184B = 427, RG185F = 428  & 
     , RG185B = 429, RG186F = 430  & 
     , RG186B = 431, RG187F = 432  & 
     , RG187B = 433, RG188F = 434  & 
     , RG188B = 435, RG189F = 436  & 
     , RG189B = 437, RR001F = 438  & 
     , RR001B = 439, RR002F = 440  & 
     , RR002B = 441, RR003F = 442  & 
     , RR003B = 443, RR004F = 444  & 
     , RR004B = 445, RR005F = 446  & 
     , RR005B = 447, RR006F = 448  & 
     , RR006B = 449, RR007F = 450 
      integer, parameter ::   RR007B = 451, RR008F = 452  & 
     , RR008B = 453, RR009F = 454  & 
     , RR009B = 455, RR010F = 456  & 
     , RR010B = 457, RR011F = 458  & 
     , RR011B = 459, RR012F = 460  & 
     , RR012B = 461, RR013F = 462  & 
     , RR013B = 463, RR014F = 464  & 
     , RR014B = 465, RR015F = 466  & 
     , RR015B = 467, RR016F = 468  & 
     , RR016B = 469, RR017F = 470  & 
     , RR017B = 471, RR018F = 472  & 
     , RR018B = 473, RR019F = 474  & 
     , RR019B = 475, RR020F = 476  & 
     , RR020B = 477, RR021F = 478  & 
     , RR021B = 479, RR022F = 480 
      integer, parameter ::   RR022B = 481, RR023F = 482  & 
     , RR023B = 483, RR024F = 484  & 
     , RR024B = 485, RR025F = 486  & 
     , RR025B = 487, RR026F = 488  & 
     , RR026B = 489, RR027F = 490  & 
     , RR027B = 491, RR028F = 492  & 
     , RR028B = 493, RR029F = 494  & 
     , RR029B = 495, RR030F = 496  & 
     , RR030B = 497, RR031F = 498  & 
     , RR031B = 499, RR032F = 500  & 
     , RR032B = 501, RR033 = 502  & 
     , RR034F = 503, RR034B = 504  & 
     , RR035F = 505, RR035B = 506  & 
     , RR036F = 507, RR036B = 508  & 
     , RR037F = 509, RR037B = 510 
      integer, parameter ::   RR038F = 511, RR038B = 512  & 
     , RR039F = 513, RR039B = 514  & 
     , RR040F = 515, RR040B = 516  & 
     , RR041F = 517, RR041B = 518  & 
     , RR042F = 519, RR042B = 520  & 
     , RR043F = 521, RR043B = 522  & 
     , RR044F = 523, RR044B = 524  & 
     , RR045F = 525, RR045B = 526  & 
     , RR046F = 527, RR046B = 528  & 
     , RR200F = 529, RR200B = 530  & 
     , RR047 = 531, RR048 = 532  & 
     , RR049 = 533, RR050 = 534  & 
     , RR051 = 535, RR052 = 536  & 
     , RR053F = 537, RR053B = 538  & 
     , RR054F = 539, RR054B = 540 
      integer, parameter ::   RR055F = 541, RR055B = 542  & 
     , RR056F = 543, RR056B = 544  & 
     , RR057F = 545, RR057B = 546  & 
     , RR058F = 547, RR058B = 548  & 
     , RR059F = 549, RR059B = 550  & 
     , RR060F = 551, RR060B = 552  & 
     , RR061F = 553, RR061B = 554  & 
     , RR062F = 555, RR062B = 556  & 
     , RR063F = 557, RR063B = 558  & 
     , RR064F = 559, RR064B = 560  & 
     , RR065F = 561, RR065B = 562  & 
     , RR066F = 563, RR066B = 564  & 
     , RR067F = 565, RR067B = 566  & 
     , RR068F = 567, RR068B = 568  & 
     , RR069F = 569, RR069B = 570 
      integer, parameter ::   RR070F = 571, RR070B = 572  & 
     , RR071F = 573, RR071B = 574  & 
     , RR073F = 575, RR073B = 576  & 
     , RR074F = 577, RR074B = 578  & 
     , RR075F = 579, RR075B = 580  & 
     , RR076F = 581, RR076B = 582  & 
     , RR077F = 583, RR077B = 584  & 
     , RR078F = 585, RR078B = 586  & 
     , RR079F = 587, RR079B = 588  & 
     , RR080F = 589, RR080B = 590  & 
     , RR081F = 591, RR081B = 592  & 
     , RR082F = 593, RR082B = 594  & 
     , RR083F = 595, RR083B = 596  & 
     , RR084F = 597, RR084B = 598  & 
     , RR085F = 599, RR085B = 600 
      integer, parameter ::   RR086F = 601, RR086B = 602  & 
     , RR087F = 603, RR087B = 604  & 
     , RR088F = 605, RR088B = 606  & 
     , RR089F = 607, RR089B = 608  & 
     , RR090F = 609, RR090B = 610  & 
     , RR091F = 611, RR091B = 612  & 
     , RR092F = 613, RR092B = 614  & 
     , RR093F = 615, RR093B = 616  & 
     , RR094F = 617, RR094B = 618  & 
     , RR095 = 619, RR096 = 620  & 
     , RR097 = 621, RR098 = 622  & 
     , RR099 = 623, RR100F = 624  & 
     , RR100B = 625, RR101F = 626  & 
     , RR101B = 627, RR102F = 628  & 
     , RR102B = 629, RR103F = 630 
      integer, parameter ::   RR103B = 631, RR104F = 632  & 
     , RR104B = 633, RR105F = 634  & 
     , RR105B = 635, RR106F = 636  & 
     , RR106B = 637, RR107F = 638  & 
     , RR107B = 639, RR108F = 640  & 
     , RR108B = 641, RR109F = 642  & 
     , RR109B = 643, RR110F = 644  & 
     , RR110B = 645, RR111 = 646  & 
     , RR112F = 647, RR112B = 648  & 
     , RR113F = 649, RR113B = 650  & 
     , RR114F = 651, RR114B = 652  & 
     , RR115F = 653, RR115B = 654  & 
     , RR116F = 655, RR116B = 656  & 
     , RR117F = 657, RR117B = 658  & 
     , RR118F = 659, RR118B = 660 
      integer, parameter ::   RR119F = 661, RR119B = 662  & 
     , RR120F = 663, RR120B = 664  & 
     , RR121 = 665, RR122F = 666  & 
     , RR122B = 667, RR123F = 668  & 
     , RR123B = 669, RR124F = 670  & 
     , RR124B = 671, RR125F = 672  & 
     , RR125B = 673, RR126F = 674  & 
     , RR126B = 675, RR127 = 676  & 
     , RR128F = 677, RR128B = 678  & 
     , RR129F = 679, RR129B = 680  & 
     , RR130F = 681, RR130B = 682  & 
     , RR131F = 683, RR131B = 684  & 
     , RR132F = 685, RR132B = 686  & 
     , RR133 = 687, RR134 = 688  & 
     , RR135 = 689, RR136F = 690 
      integer, parameter ::   RR136B = 691, RR137F = 692  & 
     , RR137B = 693, RR138 = 694  & 
     , RR139F = 695, RR139B = 696  & 
     , RR140 = 697, RR141F = 698  & 
     , RR141B = 699, RR142F = 700  & 
     , RR142B = 701, RR143F = 702  & 
     , RR143B = 703, RR144F = 704  & 
     , RR144B = 705, RR145F = 706  & 
     , RR145B = 707, RR146F = 708  & 
     , RR146B = 709, RR147F = 710  & 
     , RR147B = 711, RR148F = 712  & 
     , RR148B = 713, RR149F = 714  & 
     , RR149B = 715, RR150F = 716  & 
     , RR150B = 717, RR151F = 718  & 
     , RR151B = 719, RR152F = 720 
      integer, parameter ::   RR152B = 721, RR153F = 722  & 
     , RR153B = 723, RR154F = 724  & 
     , RR154B = 725, RR155F = 726  & 
     , RR155B = 727, RR156F = 728  & 
     , RR156B = 729, RH01F = 730  & 
     , RH01B = 731, RH02F = 732  & 
     , RH02B = 733, RH03F = 734  & 
     , RH03B = 735, RH04F = 736  & 
     , RH04B = 737, RH05 = 738  & 
     , RH06F = 739, RH06B = 740  & 
     , RH07F = 741, RH07B = 742  & 
     , RH08F = 743, RH08B = 744  & 
     , RH09F = 745, RH09B = 746  & 
     , RH10F = 747, RH10B = 748  & 
     , RH11F = 749, RH11B = 750 
      integer, parameter ::   RH12F = 751, RH12B = 752  & 
     , RH13F = 753, RH13B = 754  & 
     , RH14F = 755, RH14B = 756  & 
     , RH15F = 757, RH15B = 758  & 
     , RH16F = 759, RH16B = 760  & 
     , RH17F = 761, RH17B = 762  & 
     , RH18F = 763, RH18B = 764  & 
     , RH19F = 765, RH19B = 766  & 
     , RH20F = 767, RH20B = 768  & 
     , RH21F = 769, RH21B = 770  & 
     , RH22F = 771, RH22B = 772  & 
     , RH23F = 773, RH23B = 774  & 
     , RH24F = 775, RH24B = 776  & 
     , RH25F = 777, RH25B = 778  & 
     , RH26F = 779, RH26B = 780 
      integer, parameter ::   RH27F = 781, RH27B = 782  & 
     , RH28F = 783, RH28B = 784  & 
     , RH29F = 785, RH29B = 786  & 
     , RH30F = 787, RH30B = 788  & 
     , RH31F = 789, RH31B = 790  & 
     , RH32F = 791, RH32B = 792  & 
     , RH33F = 793, RH33B = 794  & 
     , RH34F = 795, RH34B = 796  & 
     , RH35F = 797, RH35B = 798  & 
     , RH36F = 799, RH36B = 800  & 
     , RH37F = 801, RH37B = 802  & 
     , RH38F = 803, RH38B = 804  & 
     , RH39F = 805, RH39B = 806  & 
     , RH40F = 807, RH40B = 808  & 
     , RB00F = 809, RB00B = 810 
      integer, parameter ::   RB01F = 811, RB01B = 812  & 
     , RB02F = 813, RB02B = 814  & 
     , RB04F = 815, RB04B = 816  & 
     , RB05F = 817, RB05B = 818  & 
     , RB06F = 819, RB06B = 820  & 
     , RB07F = 821, RB07B = 822  & 
     , RB08F = 823, RB08B = 824  & 
     , RB09F = 825, RB09B = 826  & 
     , RB10F = 827, RB10B = 828  & 
     , RB11F = 829, RB11B = 830  & 
     , RB12F = 831, RB12B = 832  & 
     , RB14F = 833, RB14B = 834  & 
     , RB15F = 835, RB15B = 836  & 
     , RB16F = 837, RB16B = 838  & 
     , RB17F = 839, RB17B = 840 
      integer, parameter ::   RB18F = 841, RB18B = 842  & 
     , RB19F = 843, RB19B = 844  & 
     , RB20F = 845, RB20B = 846  & 
     , RB21F = 847, RB21B = 848  & 
     , RB22F = 849, RB22B = 850  & 
     , RB23F = 851, RB23B = 852  & 
     , RB24F = 853, RB24B = 854  & 
     , RB25F = 855, RB25B = 856  & 
     , RB28 = 857, RB29F = 858  & 
     , RB29B = 859, RB30F = 860  & 
     , RB30B = 861, RB31F = 862  & 
     , RB31B = 863, RB32F = 864  & 
     , RB32B = 865, RB33F = 866  & 
     , RB33B = 867, RB34F = 868  & 
     , RB34B = 869, RB35F = 870 
      integer, parameter ::   RB35B = 871, RB36F = 872  & 
     , RB36B = 873, RB37F = 874  & 
     , RB37B = 875, RB38F = 876  & 
     , RB38B = 877, RB39F = 878  & 
     , RB39B = 879, RB40F = 880  & 
     , RB40B = 881, RB99F = 882  & 
     , RB99B = 883, RB41F = 884  & 
     , RB41B = 885, RB42 = 886  & 
     , RB43F = 887, RB43B = 888  & 
     , RB44F = 889, RB44B = 890  & 
     , RB45F = 891, RB45B = 892  & 
     , RB46F = 893, RB46B = 894  & 
     , RB47F = 895, RB47B = 896  & 
     , RB48F = 897, RB48B = 898  & 
     , RB49F = 899, RB49B = 900 
      integer, parameter ::   RB50F = 901, RB50B = 902  & 
     , RB51F = 903, RB51B = 904  & 
     , RB52F = 905, RB52B = 906  & 
     , RHP00 = 907, RHP01 = 908  & 
     , RHP02 = 909, RHP03 = 910  & 
     , RHP04 = 911, RHP05 = 912  & 
     , RHP06 = 913, RHP08 = 914  & 
     , RHP09 = 915, RHP10 = 916  & 
     , RHP11 = 917, RHP12 = 918  & 
     , RHP13 = 919, RHP14 = 920  & 
     , RHP15 = 921, RHP16 = 922  & 
     , RHP17 = 923, RHP18 = 924  & 
     , RHP19 = 925, RHP20 = 926  & 
     , RHP21 = 927, RHP22 = 928  & 
     , RHP23F = 929, RHP23B = 930 
      integer, parameter ::   RHP24 = 931, RHP25 = 932  & 
     , RHP26 = 933, RHP27 = 934  & 
     , RHP28 = 935, RHP29 = 936  & 
     , RHP30 = 937, RHP31 = 938  & 
     , RHP32 = 939, RHP33 = 940  & 
     , RHP34 = 941, RHP35 = 942  & 
     , RHP36 = 943, RHP37 = 944  & 
     , RHP38 = 945, RHP39 = 946  & 
     , RHP40 = 947, RHP41 = 948  & 
     , RHP42 = 949, RHP43 = 950  & 
     , RHP44 = 951, RHP45 = 952  & 
     , RHP46 = 953, RHP47 = 954  & 
     , RHP48 = 955, RHP49 = 956  & 
     , RHP50 = 957, RHP51 = 958  & 
     , RHP52 = 959, RHP53 = 960 
      integer, parameter ::   RHP54 = 961, RHP55 = 962  & 
     , RHP58 = 963, RHP56 = 964  & 
     , RHP57 = 965, RHP59 = 966  & 
     , RHP60 = 967, RHP61 = 968  & 
     , RHP62 = 969, RHP63 = 970  & 
     , RHP64F = 971, RHP64B = 972  & 
     , RHP65F = 973, RHP65B = 974  & 
     , RHP67 = 975, RHP68 = 976  & 
     , RHP69 = 977, RHP70 = 978  & 
     , RHP71 = 979, RHP72 = 980  & 
     , RHP73 = 981, RHP74 = 982  & 
     , RHP75 = 983, RHP76 = 984  & 
     , RHP77 = 985, RHP78 = 986  & 
     , RHP79 = 987, RHP86 = 988  & 
     , RHP87 = 989, RHP88 = 990 
      integer, parameter ::   RIC00 = 991, RIC01 = 992  & 
     , RIC02 = 993, RIC03 = 994  & 
     , RIC04 = 995, RIC05 = 996  & 
     , RIC06 = 997, RIC07 = 998  & 
     , RIC08 = 999, RIC09 = 1000  & 
     , RIC10 = 1001, RIC11 = 1002  & 
     , RIC12 = 1003, RIC13 = 1004  & 
     , RIC14 = 1005, RIC15 = 1006  & 
     , RIC16 = 1007, RIC17 = 1008  & 
     , RIC18 = 1009, RIC19 = 1010  & 
     , RIC20 = 1011, RIC21 = 1012  & 
     , RIC22 = 1013, RIC23 = 1014  & 
     , RIC24 = 1015, RIC25 = 1016  & 
     , RIC26 = 1017, RIC27 = 1018  & 
     , RIC28 = 1019, RIC29 = 1020 
      integer, parameter ::   RIC30 = 1021, RIC31 = 1022  & 
     , RIC32 = 1023, RIC33 = 1024  & 
     , RIC34 = 1025, RIC35 = 1026  & 
     , RIC36 = 1027, RIC37 = 1028  & 
     , RIC38F = 1029, RIC38B = 1030  & 
     , RIC40 = 1031, RIC41F = 1032  & 
     , RIC41B = 1033, RIC42 = 1034  & 
     , RIC43 = 1035, RIC44 = 1036  & 
     , RIC45 = 1037, RIC46 = 1038  & 
     , RIC47 = 1039, RIC48 = 1040  & 
     , RIC49 = 1041, RIC50 = 1042  & 
     , RIC51 = 1043, RIC52 = 1044  & 
     , RIC57 = 1045, RIC53 = 1046  & 
     , RIC54 = 1047, RIC55 = 1048  & 
     , RIC56 = 1049, RIC58 = 1050 
      integer, parameter ::   RIC59 = 1051, RIC60 = 1052  & 
     , RIC61 = 1053, RIC62 = 1054  & 
     , RIC63 = 1055, RIC64 = 1056  & 
     , RIC65 = 1057, RIC66 = 1058  & 
     , RIC67 = 1059, RIC68 = 1060  & 
     , RIC69 = 1061, RIC70 = 1062  & 
     , RIC71 = 1063, RIC72 = 1064  & 
     , RP000F = 1065, RP000B = 1066  & 
     , RP001F = 1067, RP001B = 1068  & 
     , RP002F = 1069, RP002B = 1070  & 
     , RP003F = 1071, RP003B = 1072  & 
     , RP004F = 1073, RP004B = 1074  & 
     , RP005F = 1075, RP005B = 1076  & 
     , RP006F = 1077, RP006B = 1078  & 
     , RP007F = 1079, RP007B = 1080 
      integer, parameter ::   RP008 = 1081, RP009F = 1082  & 
     , RP009B = 1083, RP010F = 1084  & 
     , RP010B = 1085, RP011F = 1086  & 
     , RP011B = 1087, RK012F = 1088  & 
     , RK012B = 1089, RP013F = 1090  & 
     , RP013B = 1091, RP014F = 1092  & 
     , RP014B = 1093, RP015F = 1094  & 
     , RP015B = 1095, RP016F = 1096  & 
     , RP016B = 1097, RK017F = 1098  & 
     , RK017B = 1099, RP018F = 1100  & 
     , RP018B = 1101, RK019F = 1102  & 
     , RK019B = 1103, RK020F = 1104  & 
     , RK020B = 1105, RK021F = 1106  & 
     , RK021B = 1107, RP022F = 1108  & 
     , RP022B = 1109, RP023F = 1110 
      integer, parameter ::   RP023B = 1111, RK024F = 1112  & 
     , RK024B = 1113, RP025F = 1114  & 
     , RP025B = 1115, RP026F = 1116  & 
     , RP026B = 1117, RP027F = 1118  & 
     , RP027B = 1119, RP028F = 1120  & 
     , RP028B = 1121, RK100F = 1122  & 
     , RK100B = 1123, RK102F = 1124  & 
     , RK102B = 1125, RP104F = 1126  & 
     , RP104B = 1127, RP105F = 1128  & 
     , RP105B = 1129, RP106F = 1130  & 
     , RP106B = 1131, RP107F = 1132  & 
     , RP107B = 1133, RP108F = 1134  & 
     , RP108B = 1135, RK109F = 1136  & 
     , RK109B = 1137, RK110F = 1138  & 
     , RK110B = 1139, RP111F = 1140 
      integer, parameter ::   RP111B = 1141, RK112F = 1142  & 
     , RK112B = 1143, RK113F = 1144  & 
     , RK113B = 1145, RK114F = 1146  & 
     , RK114B = 1147, RK115F = 1148  & 
     , RK115B = 1149, RP116F = 1150  & 
     , RP116B = 1151, RP117F = 1152  & 
     , RP117B = 1153, RP118F = 1154  & 
     , RP118B = 1155, RP119F = 1156  & 
     , RP119B = 1157, RP120F = 1158  & 
     , RP120B = 1159, RP121F = 1160  & 
     , RP121B = 1161, RK122F = 1162  & 
     , RK122B = 1163, RK123F = 1164  & 
     , RK123B = 1165, RP124F = 1166  & 
     , RP124B = 1167, RK125F = 1168  & 
     , RK125B = 1169, RK126F = 1170 
      integer, parameter ::   RK126B = 1171, RP127F = 1172  & 
     , RP127B = 1173, RP128F = 1174  & 
     , RP128B = 1175, RK129F = 1176  & 
     , RK129B = 1177, RP130F = 1178  & 
     , RP130B = 1179, RP131F = 1180  & 
     , RP131B = 1181, RK132F = 1182  & 
     , RK132B = 1183, RP133F = 1184  & 
     , RP133B = 1185, RK200F = 1186  & 
     , RK200B = 1187, RP201F = 1188  & 
     , RP201B = 1189, RP202F = 1190  & 
     , RP202B = 1191, RK203F = 1192  & 
     , RK203B = 1193, RK204F = 1194  & 
     , RK204B = 1195, RK205F = 1196  & 
     , RK205B = 1197, RP206F = 1198  & 
     , RP206B = 1199, RP207F = 1200 
      integer, parameter ::   RP207B = 1201, RP208F = 1202  & 
     , RP208B = 1203, RP209F = 1204  & 
     , RP209B = 1205, RK210F = 1206  & 
     , RK210B = 1207, RP211F = 1208  & 
     , RP211B = 1209, RK212F = 1210  & 
     , RK212B = 1211, RK213F = 1212  & 
     , RK213B = 1213, RP214F = 1214  & 
     , RP214B = 1215, RP301F = 1216  & 
     , RP301B = 1217, RP302F = 1218  & 
     , RP302B = 1219, RP304F = 1220  & 
     , RP304B = 1221, RP305F = 1222  & 
     , RP305B = 1223, RP306F = 1224  & 
     , RP306B = 1225, RK401F = 1226  & 
     , RK401B = 1227, RK403F = 1228  & 
     , RK403B = 1229, RP405F = 1230 
      integer, parameter ::   RP405B = 1231, RP406F = 1232  & 
     , RP406B = 1233, RP407F = 1234  & 
     , RP407B = 1235, RP408F = 1236  & 
     , RP408B = 1237, RP409F = 1238  & 
     , RP409B = 1239, RP410F = 1240  & 
     , RP410B = 1241, RP411F = 1242  & 
     , RP411B = 1243, RP412F = 1244  & 
     , RP412B = 1245, RP413F = 1246  & 
     , RP413B = 1247, RP414F = 1248  & 
     , RP414B = 1249, RP415F = 1250  & 
     , RP415B = 1251, RP416F = 1252  & 
     , RP416B = 1253, RP417F = 1254  & 
     , RP417B = 1255, RP418F = 1256  & 
     , RP418B = 1257, RP419F = 1258  & 
     , RP419B = 1259, RK420F = 1260 
      integer, parameter ::   RK420B = 1261, RK421F = 1262  & 
     , RK421B = 1263, RP422F = 1264  & 
     , RP422B = 1265, RK423F = 1266  & 
     , RK423B = 1267, RK424F = 1268  & 
     , RK424B = 1269, RP425F = 1270  & 
     , RP425B = 1271, RP501F = 1272  & 
     , RP501B = 1273, RP502F = 1274  & 
     , RP502B = 1275, RP503F = 1276  & 
     , RP503B = 1277, RP504F = 1278  & 
     , RP504B = 1279, RP505F = 1280  & 
     , RP505B = 1281, RP506F = 1282  & 
     , RP506B = 1283, RK507F = 1284  & 
     , RK507B = 1285, RK508F = 1286  & 
     , RK508B = 1287, RK600F = 1288  & 
     , RK600B = 1289, RP601F = 1290 
      integer, parameter ::   RP601B = 1291, RK602F = 1292  & 
     , RK602B = 1293, RK603F = 1294  & 
     , RK603B = 1295, RK700F = 1296  & 
     , RK700B = 1297, RK701F = 1298  & 
     , RK701B = 1299, RP800 = 1300  & 
     , RP801 = 1301, RP802 = 1302  & 
     , RCP01F = 1303, RCP01B = 1304  & 
     , RCP02F = 1305, RCP02B = 1306  & 
     , RCP03F = 1307, RCP03B = 1308  & 
     , RCP04 = 1309, RCP05F = 1310  & 
     , RCP05B = 1311, RCP06F = 1312  & 
     , RCP06B = 1313, RCP07F = 1314  & 
     , RCP07B = 1315, RCP08F = 1316  & 
     , RCP08B = 1317, RCP09F = 1318  & 
     , RCP09B = 1319, RCP10F = 1320 
      integer, parameter ::   RCP10B = 1321, RCP11F = 1322  & 
     , RCP11B = 1323, RCP12F = 1324  & 
     , RCP12B = 1325, RCP13 = 1326  & 
     , RCP14 = 1327, RCP15 = 1328  & 
     , RCP16F = 1329, RCP16B = 1330  & 
     , RCP17 = 1331, RCP18 = 1332  & 
     , RCP19F = 1333, RCP19B = 1334  & 
     , RCP20F = 1335, RCP20B = 1336  & 
     , RCP21F = 1337, RCP21B = 1338  & 
     , RCP22 = 1339, RCP23F = 1340  & 
     , RCP23B = 1341, RCP24 = 1342  & 
     , RCP25 = 1343, RCP26 = 1344  & 
     , RCP27 = 1345, RCP28 = 1346  & 
     , RCP29 = 1347, RCP30 = 1348  & 
     , RCP31 = 1349, RCP32F = 1350 
      integer, parameter ::   RCP32B = 1351, RCP33F = 1352  & 
     , RCP33B = 1353, RCP34 = 1354  & 
     , RI00F = 1355, RI00B = 1356  & 
     , RI01F = 1357, RI01B = 1358  & 
     , RI02F = 1359, RI02B = 1360  & 
     , RI03F = 1361, RI03B = 1362  & 
     , RI05F = 1363, RI05B = 1364  & 
     , RI06F = 1365, RI06B = 1366  & 
     , RI07F = 1367, RI07B = 1368  & 
     , RI08F = 1369, RI08B = 1370  & 
     , RI09F = 1371, RI09B = 1372  & 
     , RI12 = 1373, RI15 = 1374  & 
     , RI17 = 1375, RI18 = 1376  & 
     , RI19F = 1377, RI19B = 1378  & 
     , RI20F = 1379, RI20B = 1380 
      integer, parameter ::   RI21 = 1381, RI22 = 1382  & 
     , RI23 = 1383, RI25 = 1384  & 
     , RI26 = 1385, RI31 = 1386  & 
     , RI32 = 1387, RT01F = 1388  & 
     , RT01B = 1389, RT02F = 1390  & 
     , RT02B = 1391, RT03F = 1392  & 
     , RT03B = 1393, RT04F = 1394  & 
     , RT04B = 1395, RT05F = 1396  & 
     , RT05B = 1397, RT06F = 1398  & 
     , RT06B = 1399, RT07F = 1400  & 
     , RT07B = 1401, RT08F = 1402  & 
     , RT08B = 1403, RT09F = 1404  & 
     , RT09B = 1405, RT10F = 1406  & 
     , RT10B = 1407, RT11F = 1408  & 
     , RT11B = 1409, RT12F = 1410 
      integer, parameter ::   RT12B = 1411, RT13F = 1412  & 
     , RT13B = 1413, RT14F = 1414  & 
     , RT14B = 1415, RT15F = 1416  & 
     , RT15B = 1417, RT16F = 1418  & 
     , RT16B = 1419, RT17F = 1420  & 
     , RT17B = 1421, RT18F = 1422  & 
     , RT18B = 1423, RT19F = 1424  & 
     , RT19B = 1425, RT20 = 1426  & 
     , RT21F = 1427, RT21B = 1428  & 
     , RT22F = 1429, RT22B = 1430  & 
     , RT23F = 1431, RT23B = 1432  & 
     , RT24F = 1433, RT24B = 1434  & 
     , RT25F = 1435, RT25B = 1436  & 
     , RT26F = 1437, RT26B = 1438  & 
     , RT27F = 1439, RT27B = 1440 
      integer, parameter ::   RT28F = 1441, RT28B = 1442  & 
     , RT29F = 1443, RT29B = 1444  & 
     , RT30F = 1445, RT30B = 1446  & 
     , RT31F = 1447, RT31B = 1448  & 
     , RT32F = 1449, RT32B = 1450  & 
     , RT33F = 1451, RT33B = 1452  & 
     , RT34F = 1453, RT34B = 1454  & 
     , RT35F = 1455, RT35B = 1456  & 
     , RT36 = 1457, RT37 = 1458  & 
     , RT38 = 1459, RT39 = 1460  & 
     , RT40 = 1461, RT41 = 1462  & 
     , RT42 = 1463, RT43F = 1464  & 
     , RT43B = 1465, RT44F = 1466  & 
     , RT44B = 1467, RT45F = 1468  & 
     , RT45B = 1469, RT46F = 1470 
      integer, parameter ::   RT46B = 1471, RT47 = 1472  & 
     , RT50F = 1473, RT50B = 1474  & 
     , RT51F = 1475, RT51B = 1476  & 
     , RT52F = 1477, RT52B = 1478  & 
     , RT53F = 1479, RT53B = 1480  & 
     , RT54F = 1481, RT54B = 1482  & 
     , RT55F = 1483, RT55B = 1484  & 
     , RT56F = 1485, RT56B = 1486  & 
     , RT57F = 1487, RT57B = 1488  & 
     , RT58F = 1489, RT58B = 1490  & 
     , RT59 = 1491, RT60 = 1492  & 
     , RE01F = 1493, RE01B = 1494  & 
     , RE02F = 1495, RE02B = 1496  & 
     , RE03F = 1497, RE03B = 1498  & 
     , RE04F = 1499, RE04B = 1500 
      integer, parameter ::   RE05F = 1501, RE05B = 1502  & 
     , RE06F = 1503, RE06B = 1504  & 
     , RE07F = 1505, RE07B = 1506  & 
     , RE08F = 1507, RE08B = 1508  & 
     , RE09F = 1509, RE09B = 1510  & 
     , RE10F = 1511, RE10B = 1512  & 
     , RE11F = 1513, RE11B = 1514  & 
     , RE12F = 1515, RE12B = 1516  & 
     , RE13F = 1517, RE13B = 1518  & 
     , RE14F = 1519, RE14B = 1520  & 
     , RE15F = 1521, RE15B = 1522  & 
     , RE16F = 1523, RE16B = 1524  & 
     , RE17F = 1525, RE17B = 1526  & 
     , RE18F = 1527, RE18B = 1528  & 
     , RE19 = 1529, RE30 = 1530 
      integer, parameter ::   RE31 = 1531, RE32 = 1532  & 
     , RE33 = 1533, RE34 = 1534  & 
     , RE35 = 1535, RE36 = 1536  & 
     , RE37 = 1537, RST01 = 1538  & 
     , RST02F = 1539, RST02B = 1540  & 
     , RST03F = 1541, RST03B = 1542  & 
     , RST04F = 1543, RST04B = 1544  & 
     , RST05F = 1545, RST05B = 1546  & 
     , RST06F = 1547, RST06B = 1548  & 
     , RST10F = 1549, RST10B = 1550  & 
     , RST11F = 1551, RST11B = 1552  & 
     , RST12F = 1553, RST12B = 1554  & 
     , RST13 = 1555, RST14F = 1556  & 
     , RST14B = 1557, RST00F = 1558  & 
     , RST00B = 1559, RXY00F = 1560 
      integer, parameter ::   RXY00B = 1561, RXY01F = 1562  & 
     , RXY01B = 1563, RXY02F = 1564  & 
     , RXY02B = 1565, RXY03F = 1566  & 
     , RXY03B = 1567, RXY04F = 1568  & 
     , RXY04B = 1569, RXY05F = 1570  & 
     , RXY05B = 1571, RXY06F = 1572  & 
     , RXY06B = 1573, RXY07F = 1574  & 
     , RXY07B = 1575, RXY09F = 1576  & 
     , RXY09B = 1577, RXY10F = 1578  & 
     , RXY10B = 1579, RXY11 = 1580  & 
     , RXY12 = 1581, RXY13F = 1582  & 
     , RXY13B = 1583, RXY14F = 1584  & 
     , RXY14B = 1585, RXY15F = 1586  & 
     , RXY15B = 1587, RXY16F = 1588  & 
     , RXY16B = 1589, RXY17 = 1590 
      integer, parameter ::   RXY18F = 1591, RXY18B = 1592  & 
     , RXY19F = 1593, RXY19B = 1594  & 
     , RXY201 = 1595, RXY202 = 1596  & 
     , RXY203 = 1597, RXY22 = 1598  & 
     , RXY23F = 1599, RXY23B = 1600  & 
     , RXY24 = 1601, RXY25 = 1602  & 
     , RXY26F = 1603, RXY26B = 1604  & 
     , RXY27F = 1605, RXY27B = 1606  & 
     , RXY28F = 1607, RXY28B = 1608  & 
     , RXY29F = 1609, RXY29B = 1610  & 
     , RXY30F = 1611, RXY30B = 1612  & 
     , RXY31F = 1613, RXY31B = 1614  & 
     , RXY33 = 1615, RXY34 = 1616  & 
     , RXY35 = 1617, RXY36 = 1618  & 
     , RXY37 = 1619, RXY38 = 1620 
      integer, parameter ::   RXY39F = 1621, RXY39B = 1622  & 
     , RXY40F = 1623, RXY40B = 1624  & 
     , RXY41F = 1625, RXY41B = 1626  & 
     , RXY42 = 1627, RXY43F = 1628  & 
     , RXY43B = 1629, RXY44 = 1630  & 
     , RXY45F = 1631, RXY45B = 1632  & 
     , RXY46 = 1633, RXY47F = 1634  & 
     , RXY47B = 1635, RXY48 = 1636  & 
     , RXY50 = 1637, RXY51 = 1638  & 
     , RXY52 = 1639, RXY53 = 1640  & 
     , RXY54 = 1641, RXY55 = 1642  & 
     , RXY56 = 1643, RXY57F = 1644  & 
     , RXY57B = 1645, RXY58 = 1646  & 
     , RN01F = 1647, RN01B = 1648  & 
     , RN02F = 1649, RN02B = 1650 
      integer, parameter ::   RN04F = 1651, RN04B = 1652  & 
     , RN05F = 1653, RN05B = 1654  & 
     , RN06F = 1655, RN06B = 1656  & 
     , RN07F = 1657, RN07B = 1658  & 
     , RN08F = 1659, RN08B = 1660  & 
     , RN09F = 1661, RN09B = 1662  & 
     , RN10F = 1663, RN10B = 1664  & 
     , RN11F = 1665, RN11B = 1666  & 
     , RN12F = 1667, RN12B = 1668  & 
     , RN13F = 1669, RN13B = 1670  & 
     , RN14 = 1671, RN15 = 1672  & 
     , RN16F = 1673, RN16B = 1674  & 
     , RN17 = 1675, RN18F = 1676  & 
     , RN18B = 1677, RN18 = 1678  & 
     , RN20F = 1679, RN20B = 1680 
      integer, parameter ::   RN21F = 1681, RN21B = 1682  & 
     , RN22F = 1683, RN22B = 1684  & 
     , RN23F = 1685, RN23B = 1686  & 
     , RN24F = 1687, RN24B = 1688  & 
     , RN25F = 1689, RN25B = 1690  & 
     , RN26F = 1691, RN26B = 1692  & 
     , RN27F = 1693, RN27B = 1694  & 
     , RN28 = 1695, RN29 = 1696  & 
     , RN30 = 1697, RN31 = 1698  & 
     , RN32 = 1699, RN33 = 1700  & 
     , RN34 = 1701, ROX00F = 1702  & 
     , ROX00B = 1703, ROX01F = 1704  & 
     , ROX01B = 1705, ROX02F = 1706  & 
     , ROX02B = 1707, ROX03F = 1708  & 
     , ROX03B = 1709, ROX04F = 1710 
      integer, parameter ::   ROX04B = 1711, ROX05F = 1712  & 
     , ROX05B = 1713, ROX06F = 1714  & 
     , ROX06B = 1715, ROX99F = 1716  & 
     , ROX99B = 1717, ROX07F = 1718  & 
     , ROX07B = 1719, ROX08F = 1720  & 
     , ROX08B = 1721, ROX09F = 1722  & 
     , ROX09B = 1723, ROX10F = 1724  & 
     , ROX10B = 1725, ROX11F = 1726  & 
     , ROX11B = 1727, ROX12F = 1728  & 
     , ROX12B = 1729, ROX13F = 1730  & 
     , ROX13B = 1731, ROX14F = 1732  & 
     , ROX14B = 1733, ROX15F = 1734  & 
     , ROX15B = 1735, ROX16F = 1736  & 
     , ROX16B = 1737, ROX17F = 1738  & 
     , ROX17B = 1739, ROX18F = 1740 
      integer, parameter ::   ROX18B = 1741, ROX19F = 1742  & 
     , ROX19B = 1743, ROX20F = 1744  & 
     , ROX20B = 1745, ROX21F = 1746  & 
     , ROX21B = 1747, ROX22F = 1748  & 
     , ROX22B = 1749, ROX23F = 1750  & 
     , ROX23B = 1751, ROX24F = 1752  & 
     , ROX24B = 1753, ROX25F = 1754  & 
     , ROX25B = 1755, ROX26F = 1756  & 
     , ROX26B = 1757, ROX27F = 1758  & 
     , ROX27B = 1759, ROX28 = 1760  & 
     , ROX30F = 1761, ROX30B = 1762  & 
     , ROX31F = 1763, ROX31B = 1764  & 
     , ROX32F = 1765, ROX32B = 1766  & 
     , ROX33F = 1767, ROX33B = 1768  & 
     , ROX34F = 1769, ROX34B = 1770 
      integer, parameter ::   ROX35 = 1771, ROX36 = 1772  & 
     , ROX37F = 1773, ROX37B = 1774  & 
     , ROX38F = 1775, ROX38B = 1776  & 
     , ROX39F = 1777, ROX39B = 1778  & 
     , ROX40F = 1779, ROX40B = 1780  & 
     , ROX41F = 1781, ROX41B = 1782  & 
     , ROX42F = 1783, ROX42B = 1784  & 
     , ROX43F = 1785, ROX43B = 1786  & 
     , ROX44F = 1787, ROX44B = 1788  & 
     , ROX45F = 1789, ROX45B = 1790  & 
     , ROX46F = 1791, ROX46B = 1792  & 
     , ROX47 = 1793, ROX48 = 1794  & 
     , ROX50 = 1795, ROX51 = 1796  & 
     , ROX52 = 1797, ROX53 = 1798  & 
     , ROX54 = 1799, ROX60 = 1800 
      integer, parameter ::   ROX61 = 1801, ROX62 = 1802  & 
     , ROX63 = 1803, ROX64 = 1804  & 
     , ROX65 = 1805, ROX66 = 1806 & 
     , REND = 1807 
     integer, parameter ::   MM1 = 1, MM2 = 2 & 
     , MM3 = 3, MM4 = 4 & 
     , MM5 = 5, MM7 = 6 & 
     , MM8 = 7, MM9 = 8 & 
     , MM0 = 9, MM6 = 10 & 
     , MEND = 11 

integer, parameter :: DP=kind(1.0d0)


	real(DP) :: A(REND),N(REND),E(REND)


	integer :: IH

      DATA (A(IH),IH=1,30)  /  2.6400000000D+13, 5.2766713024D+10 & 
     , 4.5900000000D+01, 2.7517349858D+01 & 
     , 1.7300000000D+05, 1.9032628245D+06 & 
     , 3.9700000000D+01, 7.2853303338D+02 & 
     , 1.7800000000D+12, 4.1604426150D+14 & 
     , 9.0000000000D+10, 2.1035945806D+13 & 
     , 5.6200000000D+13, 1.3135779492D+16 & 
     , 5.5000000000D+14, 1.2855300215D+17 & 
     , 4.4000000000D+16, 1.1314226588D+20 & 
     , 9.4300000000D+12, 1.3213721422D+15 & 
     , 1.2000000000D+11, 8.4127616577D+15 & 
     , 6.3300000000D+16, 9.1184255084D+19 & 
     , 5.9200000000D+02, 3.6485333254D+03 & 
     , 2.0100000000D+14, 7.1253840510D+22 & 
     , 3.9700000000D+09, 1.4164546439D+07 /
      DATA (A(IH),IH=31,60)  /  7.4900000000D+10, 1.4562476643D+07 & 
     , 4.0000000000D+10, 3.8909647883D+09 & 
     , 2.3800000000D+10, 4.2484744234D+10 & 
     , 1.0000000000D+13, 1.7850732872D+13 & 
     , 1.3000000000D+08, 6.2200353343D+09 & 
     , 3.6600000000D+11, 1.7511791788D+13 & 
     , 6.0500000000D+03, 2.0516782971D+01 & 
     , 2.4100000000D+10, 1.7481432523D+05 & 
     , 9.6300000000D+03, 1.9578260237D+01 & 
     , 2.0000000000D+09, 7.4616787480D+07 & 
     , 2.6700000000D+38, 9.9613411288D+36 & 
     , 1.1700000000D+21, 8.9757498242D+29 & 
     , 8.0000000000D+08, 4.3798719198D+15 & 
     , 8.7800000000D+07, 4.8069094320D+14 & 
     , 1.1200000000D+09, 1.2255909945D+13 /
      DATA (A(IH),IH=61,90)  /  3.0100000000D+10, 3.2039914124D+13 & 
     , 1.2000000000D+11, 5.2955867164D+09 & 
     , 3.0000000000D+10, 7.9368470794D+08 & 
     , 3.0000000000D+10, 4.3452967068D+15 & 
     , 3.0200000000D+10, 1.4661973921D+10 & 
     , 1.8700000000D+14, 3.5306519794D+10 & 
     , 2.2400000000D+15, 4.2292301785D+11 & 
     , 1.2000000000D+07, 3.2637034817D+06 & 
     , 5.0000000000D+10, 1.1487472704D+13 & 
     , 5.8000000000D+10, 2.6634134985D+10 & 
     , 1.6500000000D+11, 8.4149682406D+10 & 
     , 5.7000000000D+10, 4.0039790897D+12 & 
     , 3.0000000000D+10, 7.9654706507D+13 & 
     , 1.0800000000D+11, 1.7581777781D+12 & 
     , 4.8200000000D+22, 7.6936076132D+28 /
      DATA (A(IH),IH=91,120)  /  5.7100000000D+09, 8.2750138742D+14 & 
     , 6.7100000000D+10, 3.5609741582D+11 & 
     , 2.6900000000D+25, 1.2130524545D+35 & 
     , 1.9000000000D+11, 9.2145125176D+07 & 
     , 5.0700000000D+24, 1.6124629783D+31 & 
     , 2.4700000000D+21, 3.4666626717D+26 & 
     , 1.0400000000D+23, 1.0197114661D+28 & 
     , 8.0000000000D+10, 7.8223107159D+12 & 
     , 2.0000000000D+10, 1.9587375479D+15 & 
     , 1.1300000000D+04, 7.6364598733D+03 & 
     , 5.0000000000D+02, 2.0974636059D+05 & 
     , 5.8000000000D+09, 2.4000000000D+09 & 
     , 4.6980064580D+11, 5.0000000000D+09 & 
     , 2.0000000000D+10, 3.8082870216D+11 & 
     , 5.0000000000D+10, 9.0545188563D+14 /
      DATA (A(IH),IH=121,150)  /  2.6900000000D+30, 9.7210519966D+41 & 
     , 4.0000000000D+10, 1.0797267305D+18 & 
     , 1.6000000000D+12, 2.6529851155D+18 & 
     , 2.0000000000D+11, 1.5000000000D+10 & 
     , 1.0936408744D+10, 9.0000000000D+09 & 
     , 6.5618452463D+09, 3.0000000000D+10 & 
     , 1.3435867056D+09, 1.5000000000D+10 & 
     , 4.7190290127D+10, 1.5000000000D+10 & 
     , 1.0693498414D+12, 3.0000000000D+10 & 
     , 2.1421554446D+15, 7.0000000000D+10 & 
     , 2.1409471364D+13, 2.8000000000D+10 & 
     , 7.5327924346D+05, 1.2000000000D+10 & 
     , 8.3014013421D+08, 1.8800000000D+35 & 
     , 1.6292506138D+45, 3.0000000000D+10 & 
     , 2.1872817488D+10, 6.8200000000D+07 /
      DATA (A(IH),IH=151,180)  /  9.0000000000D+09, 6.5618452463D+09 & 
     , 7.0000000000D+09, 5.1036574138D+09 & 
     , 1.4000000000D+10, 1.8259393135D+08 & 
     , 1.2700000000D+29, 1.2684133243D+33 & 
     , 2.2000000000D+27, 1.9434371792D+32 & 
     , 5.7400000000D+04, 9.5590953891D+01 & 
     , 3.9000000000D+10, 3.8937084533D+07 & 
     , 3.4300000000D+06, 6.2842244016D+04 & 
     , 1.0000000000D+11, 1.0263642436D+09 & 
     , 5.6000000000D+03, 2.7500401862D+03 & 
     , 9.4600000000D+10, 1.7498844055D+17 & 
     , 3.4700000000D+35, 1.5077439481D+41 & 
     , 5.0600000000D+10, 7.0821685788D+12 & 
     , 3.3700000000D+10, 4.0000000000D+33 & 
     , 1.2469090961D+42, 5.6000000000D+04 /
      DATA (A(IH),IH=181,210)  /  1.4686430659D+03, 8.0000000000D+06 & 
     , 6.4400000000D+14, 2.3164910420D+13 & 
     , 1.3800000000D+10, 2.4338024093D+12 & 
     , 5.8700000000D+08, 1.6421391660D+08 & 
     , 3.8200000000D+28, 1.4075759346D+39 & 
     , 1.0000000000D+10, 5.9178142524D+08 & 
     , 1.4000000000D+13, 2.4700000000D+08 & 
     , 1.9900000000D+09, 1.0000000000D+10 & 
     , 1.7155506298D+11, 3.6100000000D+09 & 
     , 1.0889029554D+12, 2.4500000000D+01 & 
     , 1.5445378643D+02, 5.0000000000D+10 & 
     , 6.3085696118D+15, 3.0000000000D+10 & 
     , 1.8300922168D+15, 2.6500000000D+10 & 
     , 2.1739887746D+12, 3.3200000000D+00 & 
     , 1.0278306684D+01, 1.0000000000D+11 /
      DATA (A(IH),IH=211,240)  /  8.9109851476D+17, 1.2000000000D+10 & 
     , 7.7963340707D+16, 6.8400000000D+09 & 
     , 1.1614797579D+15, 4.6600000000D+38 & 
     , 1.6463061592D+42, 4.1500000000D+04 & 
     , 4.6919928896D+03, 2.0000000000D+10 & 
     , 5.2917790613D+07, 1.5000000000D+09 & 
     , 1.6999695549D+04, 2.6200000000D+11 & 
     , 1.0680606960D+05, 1.0000000000D+10 & 
     , 1.5862280588D+07, 5.0000000000D+09 & 
     , 1.4554402259D+08, 4.2800000000D-16 & 
     , 6.9793035518D-18, 4.3600000000D+28 & 
     , 1.3623915050D+33, 2.0000000000D+10 & 
     , 4.6805022133D+08, 1.6500000000D+08 & 
     , 1.6539583927D+04, 3.2800000000D+10 & 
     , 1.1826581729D+05, 1.0000000000D+10 /
      DATA (A(IH),IH=241,270)  /  1.4029958269D+08, 5.0000000000D+09 & 
     , 1.2873158761D+09, 1.8000000000D+10 & 
     , 2.5961607220D+09, 6.6000000000D+05 & 
     , 3.5503024345D+02, 1.0200000000D+06 & 
     , 3.2893910502D+02, 1.0000000000D+05 & 
     , 5.9179879483D+02, 6.0000000000D+10 & 
     , 4.6820647566D+15, 2.4600000000D+03 & 
     , 5.5511242907D+02, 1.6000000000D+10 & 
     , 2.6323843569D+09, 1.7000000000D+04 & 
     , 1.2716073149D+02, 4.2000000000D+03 & 
     , 2.7787158476D+02, 3.8800000000D+02 & 
     , 1.7399218488D+00, 1.3000000000D+02 & 
     , 5.1562280825D+00, 1.4400000000D+03 & 
     , 1.1850024546D+02, 6.3000000000D+03 & 
     , 4.5855143694D+03, 3.0000000000D+04 /
      DATA (A(IH),IH=271,300)  /  4.1716122698D+05, 1.0000000000D+04 & 
     , 1.2299102838D+06, 2.6000000000D+30 & 
     , 1.7761630600D+36, 5.0000000000D+10 & 
     , 2.3360600807D+07, 2.0000000000D+10 & 
     , 3.0071690750D+14, 1.0000000000D+10 & 
     , 2.4794782651D+07, 3.3100000000D+03 & 
     , 9.6742642418D+06, 1.0000000000D+11 & 
     , 1.1573048940D+05, 1.0000000000D+11 & 
     , 1.5577208903D+03, 4.2000000000D+07 & 
     , 1.3076606325D-03, 5.0000000000D+10 & 
     , 1.1388201171D+12, 3.0000000000D+10 & 
     , 3.9791438217D+10, 1.0000000000D+10 & 
     , 5.0507727073D+01, 6.3400000000D+28 & 
     , 1.4048536324D+31, 8.1000000000D+03 & 
     , 2.4981403723D+04, 1.2500000000D+04 /
      DATA (A(IH),IH=301,330)  /  3.2529178544D-02, 1.8100000000D+10 & 
     , 8.8241850413D+13, 2.6300000000D+03 & 
     , 9.8996184977D+00, 2.4100000000D+03 & 
     , 4.3482685714D+08, 7.5300000000D+03 & 
     , 5.0536550148D+07, 1.2800000000D+06 & 
     , 2.3307922343D+03, 5.0000000000D+10 & 
     , 3.8326334707D+07, 1.5000000000D+06 & 
     , 4.0698109539D-01, 1.0000000000D+10 & 
     , 4.5953776074D+06, 1.7500000000D+09 & 
     , 3.7150267058D+06, 7.5000000000D+09 & 
     , 6.3247186168D+07, 1.0000000000D+10 & 
     , 3.7197289740D+08, 1.4000000000D+27 & 
     , 7.7814512532D+32, 3.0000000000D+10 & 
     , 3.1644553172D+10, 1.0300000000D+10 & 
     , 6.4367857895D+18, 5.0000000000D+09 /
      DATA (A(IH),IH=331,360)  /  5.8023026639D+10, 1.3400000000D+03 & 
     , 8.7112228115D+03, 3.0300000000D+08 & 
     , 2.7009535860D+12, 4.5800000000D+13 & 
     , 5.5766810183D+11, 1.3200000000D+34 & 
     , 8.9644644228D+28, 6.5100000000D+34 & 
     , 1.1995390279D+23, 3.1700000000D+10 & 
     , 4.3300717766D+04, 1.8100000000D+07 & 
     , 2.3500000000D+07, 2.2000000000D+10 & 
     , 2.1470530275D+02, 1.1000000000D+10 & 
     , 1.7460739608D+07, 1.2000000000D+10 & 
     , 2.0955781346D+08, 3.0100000000D+10 & 
     , 2.9305285012D+06, 5.0000000000D+10 & 
     , 1.3886732031D+22, 2.9200000000D+09 & 
     , 1.5095453499D+08, 2.0500000000D+06 & 
     , 1.7677598030D+05, 2.0500000000D+06 /
      DATA (A(IH),IH=361,390)  /  2.9200000000D+09, 3.0100000000D+10 & 
     , 2.3400000000D+07, 3.0100000000D+09 & 
     , 2.7200000000D+03, 7.0000000000D+50 & 
     , 2.8487507449D+42, 2.0300000000D+36 & 
     , 3.7928863106D+39, 1.2700000000D+02 & 
     , 5.3406078444D-02, 7.6600000000D+06 & 
     , 8.6124762219D+09, 7.1500000000D+01 & 
     , 1.1230421728D-03, 3.8900000000D+05 & 
     , 4.2684378019D+00, 1.3100000000D-04 & 
     , 6.0605344705D-07, 3.7500000000D+33 & 
     , 1.4704165297D+43, 2.2700000000D+02 & 
     , 1.7745629538D+02, 3.0000000000D+60 & 
     , 2.1619171325D+70, 1.9900000000D+38 & 
     , 2.4724095869D+44, 2.0000000000D+09 & 
     , 2.5019333849D+08, 1.1800000000D+01 /
      DATA (A(IH),IH=391,420)  /  2.7441404727D+03, 3.1700000000D+10 & 
     , 9.3219772183D+18, 1.3200000000D+20 & 
     , 3.6998605525D+08, 5.4500000000D+15 & 
     , 3.0228969809D+10, 2.2900000000D+07 & 
     , 1.8296940316D+05, 1.9200000000D+04 & 
     , 1.4802790269D+04, 5.6400000000D+74 & 
     , 8.9945859647D+61, 3.7200000000D+65 & 
     , 1.7632750149D+54, 1.7000000000D+02 & 
     , 3.1981668158D-02, 3.1700000000D-02 & 
     , 3.5752414238D-06, 1.6100000000D+03 & 
     , 3.3321974605D+00, 4.0000000000D+10 & 
     , 2.3015482642D+09, 8.4300000000D+11 & 
     , 2.9482095678D+11, 3.1700000000D+10 & 
     , 1.1503526688D+06, 4.4200000000D+58 & 
     , 7.1858177017D+64, 2.4100000000D+10 /
      DATA (A(IH),IH=421,450)  /  1.5941539969D+10, 3.3100000000D+09 & 
     , 3.6997054316D+11, 6.2600000000D+35 & 
     , 2.4335109175D+39, 3.7000000000D+13 & 
     , 1.3710677099D+13, 5.8000000000D-05 & 
     , 8.3386088338D-09, 2.3500000000D-03 & 
     , 2.0254773459D-07, 5.3600000000D+03 & 
     , 8.4777978991D+00, 9.0300000000D-04 & 
     , 2.4134146435D-04, 9.6400000000D+00 & 
     , 4.0868505075D-01, 2.4500000000D+12 & 
     , 5.2538461199D+09, 3.3000000000D+09 & 
     , 4.3822111283D+10, 1.0000000000D+10 & 
     , 8.9645037186D+07, 1.9000000000D+11 & 
     , 2.5098929533D+16, 3.4600000000D+09 & 
     , 1.8216946136D+04, 8.9500000000D+10 & 
     , 1.2083246981D+05, 7.4500000000D+40 /
      DATA (A(IH),IH=451,480)  /  1.5255854924D+50, 7.8000000000D+10 & 
     , 5.6362388680D+19, 1.0000000000D+08 & 
     , 1.5287954728D+07, 1.2100000000D+07 & 
     , 9.7578085241D+07, 9.0000000000D+10 & 
     , 9.4447045999D+12, 1.8000000000D+10 & 
     , 6.5741679900D+24, 9.0300000000D+09 & 
     , 1.7706961729D+13, 4.0400000000D+42 & 
     , 3.9209380170D+30, 1.9300000000D+15 & 
     , 1.0584325883D+24, 5.9300000000D+51 & 
     , 1.1141410801D+55, 2.4100000000D+10 & 
     , 3.0422905514D+16, 5.0000000000D+10 & 
     , 6.3084963399D+04, 5.0000000000D+10 & 
     , 4.4314188478D+06, 2.0000000000D+10 & 
     , 1.2649965411D+04, 2.0000000000D+10 & 
     , 2.5283980856D+01, 5.0000000000D+10 /
      DATA (A(IH),IH=481,510)  /  3.7594731564D+11, 3.0000000000D+10 & 
     , 1.3558457941D+10, 1.0000000000D+10 & 
     , 3.1185275273D+07, 5.0000000000D+09 & 
     , 1.3533569131D+07, 1.0000000000D+09 & 
     , 2.4177349976D+09, 3.6000000000D+10 & 
     , 1.1575885276D+10, 9.0000000000D+07 & 
     , 4.1464938961D+06, 7.5000000000D+09 & 
     , 1.0687278636D+10, 4.9000000000D+11 & 
     , 8.1202932548D+08, 1.2000000000D+10 & 
     , 2.6332608966D+16, 4.2200000000D+10 & 
     , 1.0936945048D+08, 4.9000000000D+09 & 
     , 1.2000000000D+11, 2.8148930110D+13 & 
     , 3.0000000000D+08, 2.5874528983D+11 & 
     , 3.0000000000D+08, 1.1066571938D+10 & 
     , 3.1000000000D+10, 1.2648832430D+13 /
      DATA (A(IH),IH=511,540)  /  2.6100000000D-01, 1.4479008027D-02 & 
     , 1.3600000000D+11, 5.8014628471D+18 & 
     , 1.0000000000D+10, 5.3108215024D+09 & 
     , 1.2500000000D+08, 7.7263247114D+01 & 
     , 5.0000000000D+10, 3.1525378438D+18 & 
     , 5.0000000000D+10, 2.8554902815D+17 & 
     , 5.0000000000D+09, 1.0892535033D+17 & 
     , 1.0000000000D+10, 4.8188337806D+10 & 
     , 5.0000000000D+10, 1.9579560457D+22 & 
     , 3.4600000000D+09, 6.0360475259D+03 & 
     , 2.0500000000D+06, 2.9200000000D+09 & 
     , 3.0100000000D+10, 2.3400000000D+07 & 
     , 3.0100000000D+09, 2.7200000000D+03 & 
     , 2.8000000000D+27, 3.8610572228D+32 & 
     , 1.1000000000D+07, 1.8645069743D+04 /
      DATA (A(IH),IH=541,570)  /  7.9400000000D+26, 8.1610655573D+31 & 
     , 3.1600000000D+26, 1.2666373930D+32 & 
     , 7.5300000000D+03, 2.3470297119D+15 & 
     , 1.2800000000D+06, 1.1463339595D+05 & 
     , 1.1300000000D+02, 2.1071838653D+00 & 
     , 1.8800000000D+33, 1.0162158826D+33 & 
     , 1.3800000000D+11, 4.2690327506D+13 & 
     , 1.7000000000D+02, 1.2306728046D+00 & 
     , 8.0000000000D+08, 2.5061786983D-02 & 
     , 3.0000000000D+08, 8.3477672646D+10 & 
     , 3.0000000000D+08, 2.1405800537D+10 & 
     , 4.0000000000D+11, 7.6438171011D-03 & 
     , 2.5000000000D+10, 1.8919886866D+12 & 
     , 2.5000000000D+10, 4.8515406766D+11 & 
     , 5.0000000000D+10, 1.2426252739D+16 /
      DATA (A(IH),IH=571,600)  /  5.0000000000D+10, 7.7450588497D+17 & 
     , 7.7600000000D+39, 1.9898615631D+39 & 
     , 2.4700000000D+12, 6.3337088411D+11 & 
     , 2.0100000000D+46, 3.2976487096D+51 & 
     , 6.7000000000D+39, 7.2625724752D+42 & 
     , 8.8300000000D+49, 3.7326325453D+53 & 
     , 1.5300000000D+46, 1.6495716481D+50 & 
     , 8.5000000000D+01, 1.9329120587D-01 & 
     , 4.4900000000D+04, 6.1211557595D+01 & 
     , 8.0500000000D+02, 2.0139176673D+01 & 
     , 4.2200000000D+11, 1.7839562769D+12 & 
     , 1.3000000000D-01, 8.7173128474D-02 & 
     , 1.3300000000D+03, 7.7554312548D-01 & 
     , 1.3100000000D-04, 8.4038495518D-07 & 
     , 2.2700000000D+02, 2.4607004806D+02 /
      DATA (A(IH),IH=601,630)  /  9.7600000000D+07, 1.6782252051D+07 & 
     , 9.6300000000D+03, 1.2470015987D-01 & 
     , 4.0500000000D+03, 6.5763711850D-02 & 
     , 6.2500000000D+03, 7.6307688077D-01 & 
     , 1.0000000000D+10, 1.7042930120D+10 & 
     , 1.0000000000D+10, 6.6463486802D+10 & 
     , 2.4100000000D+03, 2.2893692008D+03 & 
     , 7.5300000000D+03, 2.6607561038D+02 & 
     , 1.2800000000D+06, 2.0838122248D+03 & 
     , 4.0900000000D+06, 5.8400000000D+09 & 
     , 2.8900000000D+05, 4.0900000000D+01 & 
     , 3.4900000000D-11, 7.0600000000D+56 & 
     , 4.6645746281D+54, 5.0000000000D+51 & 
     , 8.4256335535D+49, 1.5000000000D+48 & 
     , 3.8257490317D+48, 9.5600000000D+00 /
      DATA (A(IH),IH=631,660)  /  1.3619748690D-02, 6.0300000000D+09 & 
     , 9.4510742621D+07, 4.8600000000D+08 & 
     , 1.2871407986D+09, 2.4100000000D+09 & 
     , 8.1647196307D+09, 9.6400000000D+08 & 
     , 7.3002216160D+09, 1.0000000000D+09 & 
     , 1.1451889853D+07, 2.0600000000D+01 & 
     , 1.8087353764D-01, 3.3600000000D+02 & 
     , 4.4602469311D+02, 9.7100000000D+17 & 
     , 3.0800000000D+06, 7.0070474526D+03 & 
     , 3.1700000000D+10, 8.6220550094D+15 & 
     , 4.2000000000D+29, 4.6528481728D+32 & 
     , 6.0000000000D+10, 2.1283843105D+10 & 
     , 2.6600000000D+09, 3.4693729364D+09 & 
     , 1.0600000000D+13, 4.0003442177D+12 & 
     , 7.0000000000D+09, 9.1127053346D+06 /
      DATA (A(IH),IH=661,690)  /  3.3400000000D+09, 1.8467642799D+08 & 
     , 6.0000000000D+10, 7.0278140237D+07 & 
     , 5.0000000000D+09, 2.0000000000D+10 & 
     , 3.2504162150D+01, 9.0000000000D+10 & 
     , 4.8320783021D+12, 1.0000000000D+08 & 
     , 1.0278827655D+10, 3.3400000000D+09 & 
     , 7.2407949316D+07, 6.0000000000D+10 & 
     , 3.5987068607D+06, 5.0000000000D+09 & 
     , 2.0000000000D+10, 1.6644286678D+00 & 
     , 9.0000000000D+10, 1.8945616645D+12 & 
     , 1.0000000000D+08, 4.0301236059D+09 & 
     , 1.3400000000D+03, 4.5663206858D+02 & 
     , 6.7000000000D+02, 8.9518169802D+01 & 
     , 3.0300000000D+08, 3.0300000000D+08 & 
     , 4.5800000000D+13, 4.5800000000D+13 /
      DATA (A(IH),IH=691,720)  /  1.1941331843D+11, 1.9200000000D+04 & 
     , 2.5515359785D+04, 1.0000000000D+09 & 
     , 1.0000000000D+14, 2.4418106508D+11 & 
     , 2.0300000000D+12, 8.0000000000D+18 & 
     , 4.3155004306D+12, 1.2000000000D+05 & 
     , 7.2781615229D+01, 3.5000000000D+04 & 
     , 3.5179051144D-02, 6.6000000000D+02 & 
     , 8.2106582340D+01, 9.6500000000D+01 & 
     , 7.1970624796D+00, 2.0000000000D+05 & 
     , 2.7372640710D+05, 9.6000000000D+00 & 
     , 3.5216920977D+02, 4.5200000000D-04 & 
     , 1.0453243323D-01, 4.0000000000D+02 & 
     , 3.2877696006D-01, 6.0000000000D+07 & 
     , 2.9565590310D+04, 1.1000000000D+03 & 
     , 9.9468766212D+00, 8.4000000000D-04 /
      DATA (A(IH),IH=721,750)  /  1.2835099968D-03, 6.6500000000D+02 & 
     , 1.3940817680D+00, 1.2100000000D+08 & 
     , 1.5207082165D+05, 6.5500000000D-05 & 
     , 1.5106385521D-06, 1.1400000000D+02 & 
     , 4.4427312105D+02, 1.0000000000D+10 & 
     , 1.5210017608D-01, 6.0000000000D+10 & 
     , 2.7921110540D+18, 3.2000000000D+06 & 
     , 1.6072671724D+01, 4.0000000000D+11 & 
     , 2.2060311607D+17, 1.5100000000D+11 & 
     , 1.5100000000D+10, 2.5913304529D+13 & 
     , 9.5600000000D+09, 4.2262757709D+05 & 
     , 2.0600000000D+04, 2.2950590346D-02 & 
     , 2.3000000000D+42, 1.2501321284D+45 & 
     , 1.3700000000D+36, 4.7217121748D+38 & 
     , 9.1500000000D+06, 5.0560557267D+02 /
      DATA (A(IH),IH=751,780)  /  3.3000000000D+09, 3.6180564088D+06 & 
     , 4.1000000000D+43, 6.4659597175D+43 & 
     , 2.5000000000D+17, 3.9426583642D+17 & 
     , 2.0000000000D+44, 7.4803555018D+49 & 
     , 3.4000000000D+40, 8.0634708732D+45 & 
     , 6.3000000000D+22, 5.9560076535D+19 & 
     , 2.8000000000D+20, 1.6785086788D+17 & 
     , 1.5000000000D+10, 1.0172597244D+10 & 
     , 3.0000000000D+10, 1.2900683124D+10 & 
     , 2.5000000000D+09, 1.8652337345D+10 & 
     , 5.0000000000D+09, 2.3654518883D+10 & 
     , 6.7000000000D+02, 2.8003479993D+03 & 
     , 1.3400000000D+03, 3.5513449818D+03 & 
     , 2.0000000000D+10, 1.6504795180D+07 & 
     , 1.6300000000D+08, 4.0425156849D+06 /
      DATA (A(IH),IH=781,810)  /  1.7000000000D+02, 3.6988379195D-10 & 
     , 1.2700000000D+02, 7.9365363993D-02 & 
     , 6.3500000000D+01, 6.2582103235D-02 & 
     , 6.5500000000D-05, 4.5032001812D-07 & 
     , 3.2800000000D-05, 3.5563372177D-07 & 
     , 1.1400000000D+02, 1.3243742499D+02 & 
     , 5.6800000000D+01, 1.0406467934D+02 & 
     , 6.2500000000D+03, 2.9857187952D+00 & 
     , 3.5800000000D+01, 2.3000634086D-04 & 
     , 1.9500000000D+05, 1.2309062294D+00 & 
     , 7.8000000000D+10, 2.7708034384D+17 & 
     , 7.8000000000D+10, 1.8874596234D+19 & 
     , 7.8000000000D+10, 2.3090595629D+17 & 
     , 7.8000000000D+10, 2.6650315696D+19 & 
     , 1.0000000000D+09, 6.8182357335D+24 /
      DATA (A(IH),IH=811,840)  /  1.9000000000D+11, 3.5052571308D+22 & 
     , 1.3200000000D+09, 1.8479981673D+19 & 
     , 8.4300000000D+10, 3.0917302550D+27 & 
     , 1.2000000000D+19, 4.1805841138D+28 & 
     , 2.4000000000D+17, 1.5163415165D+25 & 
     , 9.6000000000D+08, 2.4080304685D+12 & 
     , 2.6000000000D+54, 1.8696668182D+69 & 
     , 7.2300000000D+08, 2.5734807876D+13 & 
     , 5.7000000000D+36, 5.4144775118D+29 & 
     , 5.3000000000D+44, 9.1303543420D+35 & 
     , 2.5000000000D+15, 1.2837469486D+05 & 
     , 8.9400000000D+04, 6.2546358130D+14 & 
     , 2.8300000000D+05, 5.0770575440D+14 & 
     , 1.3300000000D+03, 5.3552893087D-04 & 
     , 6.6500000000D+02, 1.4764640528D-02 /
      DATA (A(IH),IH=841,870)  /  2.2000000000D+08, 9.1137574052D+14 & 
     , 7.5000000000D+03, 9.9828929904D-02 & 
     , 6.2000000000D+03, 2.7464749857D-02 & 
     , 3.1000000000D+03, 7.5720868741D-01 & 
     , 2.0000000000D+11, 1.4970623133D+08 & 
     , 1.0000000000D+11, 4.1274309621D+09 & 
     , 5.0000000000D+10, 4.7875611343D+07 & 
     , 2.5000000000D+10, 1.3199402512D+09 & 
     , 7.6600000000D+06, 7.1500000000D+01 & 
     , 1.4303992598D-06, 3.8900000000D+05 & 
     , 3.5475897840D-02, 3.7500000000D+33 & 
     , 3.4254390076D+29, 1.3000000000D+48 & 
     , 1.0193779458D+52, 4.9000000000D+48 & 
     , 2.1186435579D+54, 1.5000000000D+67 & 
     , 8.2710604472D+68, 3.1000000000D+23 /
      DATA (A(IH),IH=871,900)  /  1.7093524924D+25, 1.5000000000D+10 & 
     , 4.4711466211D+08, 2.0000000000D+09 & 
     , 6.5585873964D+08, 5.0000000000D+09 & 
     , 5.4798862760D+14, 1.2100000000D+07 & 
     , 1.0190792609D+11, 6.0000000000D+08 & 
     , 2.4178165332D+14, 1.0300000000D+10 & 
     , 1.0294910143D+08, 1.3400000000D+03 & 
     , 2.4616656286D+02, 3.0300000000D+08 & 
     , 4.5800000000D+13, 2.2969084277D+13 & 
     , 3.0000000000D+10, 1.6217315723D+07 & 
     , 2.0000000000D+10, 2.9279079883D+02 & 
     , 4.0000000000D+09, 2.3788681409D+07 & 
     , 5.0000000000D+09, 9.9380598978D+12 & 
     , 6.0000000000D+08, 4.3848365309D+12 & 
     , 1.2100000000D+07, 1.8481534515D+09 /
      DATA (A(IH),IH=901,930)  /  2.1600000000D+07, 2.3456044263D+05 & 
     , 3.1700000000D+10, 6.4953481542D+04 & 
     , 1.8400000000D-16, 4.1737418083D-06 & 
     , 8.1000000000D+77, 1.4200000000D+78 & 
     , 1.7500000000D+03, 1.7200000000D+02 & 
     , 7.4000000000D+05, 2.8900000000D+10 & 
     , 7.5700000000D+09, 1.2100000000D+78 & 
     , 1.1600000000D+05, 1.8900000000D+12 & 
     , 7.7300000000D+18, 2.5300000000D+18 & 
     , 2.4900000000D+16, 1.4100000000D+15 & 
     , 6.6600000000D+08, 7.0000000000D+09 & 
     , 7.0000000000D+09, 6.1800000000D+16 & 
     , 2.2300000000D+15, 8.9200000000D+19 & 
     , 1.1900000000D+22, 1.0800000000D+03 & 
     , 3.7000000000D+10, 2.8612674099D+09 /
      DATA (A(IH),IH=931,960)  /  3.0000000000D+10, 9.8700000000D+19 & 
     , 2.0700000000D+13, 2.5000000000D+13 & 
     , 2.5000000000D+13, 7.4600000000D+21 & 
     , 8.4600000000D+14, 7.1000000000D+09 & 
     , 3.1500000000D-19, 9.1700000000D+20 & 
     , 2.8000000000D+10, 2.5400000000D+02 & 
     , 5.1200000000D+03, 2.5000000000D+13 & 
     , 1.3400000000D+15, 8.6100000000D+17 & 
     , 7.6800000000D+09, 8.3900000000D+14 & 
     , 1.0000000000D+31, 5.0000000000D+15 & 
     , 5.0000000000D+10, 2.2500000000D+10 & 
     , 2.7000000000D+10, 1.4000000000D+09 & 
     , 1.0000000000D+08, 5.0000000000D+10 & 
     , 3.0000000000D+08, 1.3000000000D+10 & 
     , 1.3000000000D+10, 7.2300000000D+02 /
      DATA (A(IH),IH=961,990)  /  1.3000000000D+10, 1.0000000000D+09 & 
     , 5.0000000000D+08, 1.0000000000D+09 & 
     , 1.0000000000D+09, 4.0000000000D+10 & 
     , 2.6900000000D+07, 2.0000000000D+10 & 
     , 1.7000000000D+09, 2.8000000000D+09 & 
     , 5.0100000000D+31, 3.8359636976D+31 & 
     , 1.8800000000D+03, 1.6201772977D+14 & 
     , 3.1600000000D+10, 1.0000000000D+06 & 
     , 8.0000000000D+09, 3.9800000000D+09 & 
     , 6.3100000000D+09, 7.0000000000D+09 & 
     , 7.0000000000D+09, 7.9400000000D+14 & 
     , 7.9400000000D+14, 7.0000000000D+09 & 
     , 7.0000000000D+09, 1.3900000000D+16 & 
     , 2.5100000000D+14, 1.0000000000D+09 & 
     , 1.0000000000D+09, 1.0000000000D+09 /
      DATA (A(IH),IH=991,1020)  /  1.0200000000D+49, 5.7500000000D+49 & 
     , 1.9400000000D+57, 5.1500000000D-02 & 
     , 1.2500000000D+01, 1.0300000000D+04 & 
     , 1.0300000000D+08, 1.1400000000D-21 & 
     , 9.8500000000D+07, 4.9300000000D+07 & 
     , 4.2800000000D+22, 2.5500000000D+39 & 
     , 4.2200000000D+24, 7.0000000000D+09 & 
     , 7.0000000000D+09, 1.3300000000D+23 & 
     , 7.9500000000D+33, 2.6900000000D+20 & 
     , 7.9200000000D+13, 1.8300000000D+17 & 
     , 1.9400000000D+18, 1.1200000000D-44 & 
     , 6.5300000000D+59, 2.0900000000D+65 & 
     , 2.3800000000D-16, 3.9300000000D+03 & 
     , 7.7600000000D-12, 1.9200000000D+02 & 
     , 1.0000000000D+10, 1.7800000000D+39 /
      DATA (A(IH),IH=1021,1050)  /  7.2500000000D+39, 3.0900000000D+13 & 
     , 1.3400000000D+71, 7.8400000000D+45 & 
     , 2.1400000000D+38, 7.0000000000D+09 & 
     , 7.0000000000D+09, 1.9200000000D+66 & 
     , 3.0700000000D+55, 1.1210368967D+52 & 
     , 5.6800000000D+30, 3.4000000000D+02 & 
     , 2.9018835603D+01, 1.2100000000D+08 & 
     , 5.2000000000D+03, 4.4200000000D-03 & 
     , 1.9300000000D+01, 1.9300000000D+01 & 
     , 1.5800000000D+04, 3.3300000000D+04 & 
     , 1.0100000000D+18, 5.0000000000D+13 & 
     , 3.0000000000D+07, 1.2300000000D+47 & 
     , 3.1000000000D+10, 6.0300000000D+10 & 
     , 2.4700000000D+10, 7.1400000000D+12 & 
     , 7.2900000000D+26, 7.0000000000D+09 /
      DATA (A(IH),IH=1051,1080)  /  2.6000000000D+09, 7.1800000000D+09 & 
     , 2.6900000000D+07, 1.0000000000D+09 & 
     , 3.9800000000D+09, 9.7700000000D-09 & 
     , 9.8800000000D+18, 1.7300000000D+10 & 
     , 4.1000000000D+08, 7.6500000000D+08 & 
     , 5.6300000000D+04, 1.1300000000D+11 & 
     , 1.0500000000D+07, 1.7000000000D+10 & 
     , 1.5000000000D+09, 6.1199043605D+09 & 
     , 4.6200000000D+12, 2.4351103303D+22 & 
     , 6.8000000000D+21, 6.5000342420D+29 & 
     , 1.4500000000D+45, 5.9159075488D+45 & 
     , 2.2400000000D+68, 5.6381479842D+62 & 
     , 9.6000000000D+67, 3.7352568967D+81 & 
     , 1.3800000000D+13, 2.9676257451D+23 & 
     , 1.6700000000D+20, 6.5129325290D+28 /
      DATA (A(IH),IH=1081,1110)  /  2.1600000000D+36, 8.2500000000D+43 & 
     , 2.9151378138D+59, 1.0700000000D+42 & 
     , 1.5425608195D+58, 5.7700000000D+34 & 
     , 5.1317962822D+44, 3.2900000000D+03 & 
     , 2.8151364037D+14, 6.0000000000D+09 & 
     , 3.8029710467D+23, 1.8700000000D+04 & 
     , 7.3122065614D+11, 9.4500000000D-06 & 
     , 2.7559045191D-05, 4.3000000000D+60 & 
     , 6.4261715627D+54, 1.2900000000D+05 & 
     , 4.5060198908D+01, 7.8000000000D+00 & 
     , 2.9974411044D-02, 5.9000000000D+10 & 
     , 1.3748876245D+10, 7.1800000000D+10 & 
     , 4.1382935254D+06, 1.6500000000D+08 & 
     , 2.2227996788D+06, 2.5000000000D+09 & 
     , 3.7051777850D+08, 4.3000000000D+60 /
      DATA (A(IH),IH=1111,1140)  /  6.1048303870D+54, 5.2300000000D+03 & 
     , 1.7355069987D+00, 1.3400000000D-01 & 
     , 4.8919514761D-04, 3.0100000000D+14 & 
     , 1.8338189281D+09, 6.3500000000D+01 & 
     , 9.0423895896D-02, 6.5500000000D-05 & 
     , 1.0261325190D-06, 1.3400000000D+01 & 
     , 3.3117446183D+14, 3.0200000000D+00 & 
     , 3.6178371167D+10, 3.8000000000D+04 & 
     , 1.0608179381D+14, 3.6000000000D+14 & 
     , 7.8690280050D+25, 3.6200000000D+25 & 
     , 9.5260034588D+36, 1.2600000000D+01 & 
     , 7.6077074478D+12, 8.6000000000D+60 & 
     , 6.5580381913D+54, 2.6500000000D+05 & 
     , 4.7232496168D+01, 9.6300000000D-01 & 
     , 1.8883129559D-03, 8.6000000000D+60 /
      DATA (A(IH),IH=1141,1170)  /  7.4400104578D+54, 2.6500000000D+05 & 
     , 5.3584662850D+01, 9.6300000000D-01 & 
     , 2.1422668990D-03, 1.2800000000D+03 & 
     , 2.2606983071D+14, 3.2900000000D+03 & 
     , 2.8326023455D+14, 6.0000000000D+09 & 
     , 4.7823578389D+18, 6.0000000000D+09 & 
     , 2.3313052517D+18, 7.2000000000D+14 & 
     , 1.0228650489D+20, 7.2000000000D+14 & 
     , 5.6568546074D+19, 3.6200000000D+25 & 
     , 1.2133507349D+31, 3.6200000000D+25 & 
     , 5.9148458473D+30, 7.1800000000D+10 & 
     , 2.5916728182D+06, 1.6500000000D+08 & 
     , 1.3920640168D+06, 2.5000000000D+09 & 
     , 2.3204271261D+08, 7.1800000000D+10 & 
     , 4.2776375363D+06, 1.6500000000D+08 /
      DATA (A(IH),IH=1171,1200)  /  2.2976454626D+06, 2.5000000000D+09 & 
     , 3.8299379864D+08, 2.1500000000D+60 & 
     , 5.8671354650D+54, 1.3200000000D+05 & 
     , 8.4193995682D+01, 6.7200000000D-02 & 
     , 4.7155149309D-04, 2.1500000000D+60 & 
     , 5.1217965796D+54, 1.3200000000D+05 & 
     , 7.3498306231D+01, 6.7200000000D-02 & 
     , 4.1164735991D-04, 2.8800000000D+11 & 
     , 2.1196844222D+08, 3.5200000000D+09 & 
     , 7.1773752755D+10, 8.6000000000D+60 & 
     , 1.2393355091D+55, 2.6500000000D+05 & 
     , 8.9259787720D+01, 9.6300000000D-01 & 
     , 3.5685264864D-03, 3.2900000000D+03 & 
     , 2.4444249624D+16, 6.0000000000D+09 & 
     , 2.0118251902D+20, 3.6000000000D+14 /
      DATA (A(IH),IH=1201,1230)  /  4.0658510808D+21, 3.6200000000D+25 & 
     , 5.1042804724D+32, 2.1500000000D+60 & 
     , 5.6385261266D+54, 1.3200000000D+05 & 
     , 8.0913428228D+01, 6.7200000000D-02 & 
     , 4.5317777810D-04, 7.1800000000D+10 & 
     , 6.0688599637D+04, 1.6500000000D+08 & 
     , 3.2597639329D+04, 2.5000000000D+09 & 
     , 5.4336902348D+06, 9.5500000000D+08 & 
     , 1.0406159274D+16, 1.3900000000D+10 & 
     , 2.4550866312D+23, 8.6000000000D+60 & 
     , 8.0262122096D+55, 4.0100000000D+05 & 
     , 8.7473417297D+02, 2.6900000000D-01 & 
     , 6.4556007352D-03, 1.3400000000D+01 & 
     , 3.2105425876D+14, 1.3400000000D+01 & 
     , 4.0290414867D+14, 3.8000000000D+04 /
      DATA (A(IH),IH=1231,1260)  /  4.7450169757D+13, 3.8000000000D+04 & 
     , 8.5798776135D+13, 3.2900000000D+03 & 
     , 1.9481403928D+13, 6.0000000000D+09 & 
     , 3.4325764465D+26, 6.0000000000D+09 & 
     , 4.3076808769D+26, 1.8000000000D+14 & 
     , 2.8101475820D+25, 1.8000000000D+14 & 
     , 3.0785672344D+25, 3.6200000000D+25 & 
     , 3.7260182077D+36, 3.6200000000D+25 & 
     , 4.6759329706D+36, 1.2600000000D+01 & 
     , 7.0239241606D+12, 1.2600000000D+01 & 
     , 6.1912766332D+12, 3.1800000000D+08 & 
     , 6.2938340689D+19, 2.3900000000D+08 & 
     , 1.9527079107D+19, 1.3900000000D+10 & 
     , 1.8408529969D+27, 4.3000000000D+60 & 
     , 8.1270195379D+54, 2.0000000000D+05 /
      DATA (A(IH),IH=1261,1290)  /  8.8351184781D+01, 1.3400000000D-01 & 
     , 6.5123816232D-04, 4.3000000000D+60 & 
     , 8.1270195379D+54, 2.6500000000D+05 & 
     , 1.1706531984D+02, 9.6300000000D-01 & 
     , 4.6801667934D-03, 1.3000000000D+11 & 
     , 5.8424565666D-02, 1.3400000000D+01 & 
     , 3.5811034647D+14, 3.8000000000D+04 & 
     , 6.8087951892D+11, 1.2800000000D+03 & 
     , 7.6582509943D+10, 6.0000000000D+09 & 
     , 3.7865926687D+17, 1.2600000000D+01 & 
     , 4.2399475729D+12, 2.1500000000D+60 & 
     , 7.1084295844D+54, 2.6500000000D+05 & 
     , 2.0478616520D+02, 9.6300000000D-01 & 
     , 8.1871677406D-03, 1.2800000000D+03 & 
     , 4.3267077058D+11, 1.7200000000D+60 /
      DATA (A(IH),IH=1291,1320)  /  1.2436066718D+54, 5.3000000000D+05 & 
     , 8.9567406666D+01, 1.9300000000D+00 & 
     , 3.5882616284D-03, 1.2800000000D+03 & 
     , 8.7943130143D+10, 1.2800000000D+03 & 
     , 1.0865488865D+11, 6.3700000000D+08 & 
     , 9.5500000000D+08, 1.3900000000D+10 & 
     , 1.7300000000D+68, 7.5313864907D+61 & 
     , 2.8000000000D+10, 2.8490915000D+06 & 
     , 3.3000000000D+11, 1.1779783274D+04 & 
     , 3.3000000000D+11, 4.7700000000D+01 & 
     , 2.9097834187D-03, 3.0800000000D+03 & 
     , 3.4478768376D+00, 1.0000000000D+11 & 
     , 6.2711180531D+07, 1.1000000000D+01 & 
     , 3.3005564203D-01, 1.8000000000D-04 & 
     , 3.4048615419D-05, 6.0000000000D+09 /
      DATA (A(IH),IH=1321,1350)  /  1.4518195750D+09, 6.0000000000D+09 & 
     , 1.5162412910D+12, 7.6600000000D+06 & 
     , 2.1945908043D+10, 3.8900000000D+05 & 
     , 3.7500000000D+33, 3.7500000000D+33 & 
     , 6.8700000000D+52, 5.5097990876D+64 & 
     , 6.3900000000D+26, 4.9100000000D+28 & 
     , 7.0000000000D+10, 1.2141817939D+15 & 
     , 4.3400000000D+04, 1.5046362505D+06 & 
     , 3.1000000000D+10, 1.5338425964D+13 & 
     , 1.0200000000D+10, 2.0000000000D+13 & 
     , 9.7282907499D+08, 1.0000000000D+12 & 
     , 6.6700000000D+09, 3.3300000000D+09 & 
     , 1.6700000000D+09, 1.4300000000D-16 & 
     , 1.3300000000D+09, 1.2800000000D+04 & 
     , 3.3700000000D+44, 2.7400000000D+06 /
      DATA (A(IH),IH=1351,1380)  /  1.0395869997D+12, 1.0000000000D+10 & 
     , 2.2229923910D+04, 7.8600000000D-04 & 
     , 1.0000000000D+10, 5.2485026204D+25 & 
     , 1.7300000000D+68, 1.8203611359D+63 & 
     , 2.8000000000D+10, 6.8863488094D+07 & 
     , 3.1600000000D+01, 4.5268324245D+08 & 
     , 4.7700000000D+01, 7.0330431931D-02 & 
     , 3.0800000000D+03, 8.3336328636D+01 & 
     , 1.0000000000D+11, 1.5157500677D+09 & 
     , 1.1000000000D+01, 7.9775545205D+00 & 
     , 1.8000000000D-04, 8.2296634648D-04 & 
     , 3.8300000000D+06, 1.8800000000D+33 & 
     , 2.5600000000D+26, 1.9600000000D+28 & 
     , 2.8000000000D+10, 2.1712467711D+13 & 
     , 1.7400000000D+04, 2.6968482979D+04 /
      DATA (A(IH),IH=1381,1410)  /  1.2400000000D+10, 4.0800000000D+09 & 
     , 4.3200000000D+36, 2.3800000000D+15 & 
     , 1.1900000000D+15, 3.3700000000D+44 & 
     , 1.3700000000D+06, 2.3100000000D+03 & 
     , 4.5426319227D-03, 1.5600000000D+13 & 
     , 1.2813196505D+10, 4.3500000000D+22 & 
     , 5.2773986505D+10, 5.8300000000D+64 & 
     , 8.6112522502D+55, 8.2000000000D+14 & 
     , 4.3614434871D+01, 2.1800000000D+04 & 
     , 2.5793217860D+04, 6.4700000000D-03 & 
     , 1.2421004032D-03, 1.7700000000D+02 & 
     , 3.7383361476D+02, 7.8300000000D-01 & 
     , 6.8369388406D-01, 3.1400000000D-02 & 
     , 1.6980894524D+04, 1.1800000000D-03 & 
     , 1.3580893683D-04, 1.6900000000D+09 /
      DATA (A(IH),IH=1411,1440)  /  1.2806518310D+17, 1.6600000000D+04 & 
     , 1.9339064481D+07, 4.2200000000D+11 & 
     , 1.5060641068D+14, 9.3300000000D+01 & 
     , 5.2817783980D+03, 7.9400000000D+10 & 
     , 1.0571054566D+14, 2.2800000000D+11 & 
     , 8.0932136961D+17, 2.0000000000D+10 & 
     , 4.8247335589D+17, 1.1900000000D+06 & 
     , 5.8610185036D+06, 4.3200000000D+36 & 
     , 3.7600000000D+12, 3.1887853062D+12 & 
     , 6.2600000000D+34, 2.6802796432D+32 & 
     , 4.6600000000D+38, 4.4376857912D+41 & 
     , 4.2000000000D+03, 1.0308564486D+03 & 
     , 1.3000000000D+02, 1.9128731618D+01 & 
     , 6.3000000000D+03, 1.7011480543D+04 & 
     , 1.0000000000D+04, 4.5627585431D+06 /
      DATA (A(IH),IH=1441,1470)  /  8.5500000000D+01, 9.8745294140D+02 & 
     , 5.2600000000D+28, 6.2875420673D+24 & 
     , 7.2100000000D+33, 4.1991617492D+20 & 
     , 2.3700000000D+32, 1.5941382341D+20 & 
     , 1.3300000000D+10, 3.7159198196D+08 & 
     , 6.6700000000D+09, 1.1172090559D+08 & 
     , 3.3300000000D+09, 1.0235559534D+09 & 
     , 2.8500000000D-16, 4.9074531718D-17 & 
     , 2.1000000000D+16, 4.0900000000D+06 & 
     , 5.8400000000D+09, 2.8900000000D+05 & 
     , 1.2000000000D+02, 4.0900000000D+01 & 
     , 3.4900000000D-11, 1.0100000000D+71 & 
     , 1.5527601509D+66, 1.1500000000D+11 & 
     , 4.1323812035D+08, 2.8100000000D+10 & 
     , 6.0534474546D+07, 2.9500000000D+03 /
      DATA (A(IH),IH=1471,1500)  /  1.1662112183D+02, 2.9000000000D+10 & 
     , 8.3900000000D+60, 5.8699915686D+54 & 
     , 3.9100000000D+05, 6.3939878670D+01 & 
     , 2.6000000000D+10, 2.5489558805D+06 & 
     , 1.3600000000D+01, 2.4467332309D-02 & 
     , 1.7900000000D-05, 5.4416006531D-06 & 
     , 1.0000000000D+11, 1.6651442300D+20 & 
     , 3.0000000000D+10, 3.5650010139D+17 & 
     , 3.0000000000D+10, 6.9312742321D+13 & 
     , 8.5700000000D+17, 2.0355198398D+22 & 
     , 2.5500000000D+10, 2.5500000000D+10 & 
     , 3.6100000000D+10, 9.9215141237D+17 & 
     , 2.0000000000D+10, 3.0096212144D+21 & 
     , 2.0000000000D+10, 1.1999372607D+25 & 
     , 2.3100000000D+03, 6.2409163294D-06 /
      DATA (A(IH),IH=1501,1530)  /  7.8300000000D-01, 9.3929607283D-04 & 
     , 4.8300000000D-05, 4.1076734557D-10 & 
     , 1.9600000000D-03, 9.9930663771D-09 & 
     , 4.4700000000D+03, 4.1822401369D-01 & 
     , 8.0300000000D+00, 2.0137726410D-02 & 
     , 7.5300000000D-04, 1.1904806541D-05 & 
     , 1.7200000000D+11, 4.2169346325D+00 & 
     , 3.7900000000D+06, 1.0596119199D+04 & 
     , 1.6700000000D+09, 1.0912986807D+09 & 
     , 2.4100000000D+10, 1.7325942283D+11 & 
     , 3.3100000000D+09, 4.0209968985D+12 & 
     , 3.7000000000D+13, 1.4901345826D+14 & 
     , 3.1700000000D+10, 8.1033372636D+08 & 
     , 3.1700000000D+10, 2.4565719683D+09 & 
     , 7.0000000000D+09, 4.8200000000D+33 /
      DATA (A(IH),IH=1531,1560)  /  3.8800000000D+01, 1.2100000000D+45 & 
     , 2.4500000000D+25, 1.2300000000D+13 & 
     , 5.4500000000D-02, 1.0000000000D+08 & 
     , 6.3000000000D+14, 2.4000000000D+14 & 
     , 5.1000000000D+04, 3.3224050070D+00 & 
     , 1.1400000000D+02, 3.0178171755D+02 & 
     , 1.8800000000D+33, 3.7188357721D+32 & 
     , 1.8800000000D+33, 1.2267086362D+32 & 
     , 7.8300000000D-01, 8.8911775253D-03 & 
     , 1.3300000000D+10, 2.6458805967D+06 & 
     , 1.0300000000D+10, 2.0794236065D+07 & 
     , 6.7000000000D+02, 5.5627276119D+01 & 
     , 3.0300000000D+08, 4.5800000000D+13 & 
     , 2.9640252599D+12, 1.0000000000D+10 & 
     , 7.2340731681D+07, 2.5000000000D+18 /
      DATA (A(IH),IH=1561,1590)  /  1.0443593278D+15, 4.3200000000D+29 & 
     , 2.4885817320D+17, 1.2900000000D-02 & 
     , 1.2595609433D-03, 2.3600000000D-03 & 
     , 1.3814507233D-04, 3.5400000000D+02 & 
     , 3.8026416345D+02, 4.3600000000D+04 & 
     , 2.6236903331D+04, 1.8700000000D+02 & 
     , 5.3841505118D+03, 8.4400000000D+11 & 
     , 1.5319708690D+14, 4.6200000000D+03 & 
     , 3.8039484327D-03, 1.5700000000D+00 & 
     , 6.9907410294D-01, 1.8200000000D+05 & 
     , 8.2000000000D+14, 5.8300000000D+64 & 
     , 8.0394580864D+55, 4.3700000000D+15 & 
     , 2.3683526004D+18, 5.9900000000D+20 & 
     , 1.1561143354D+14, 1.9700000000D+19 & 
     , 3.8721310084D+13, 2.0000000000D+10 /
      DATA (A(IH),IH=1591,1620)  /  1.3800000000D-01, 1.4948586364D-01 & 
     , 6.5700000000D+00, 3.0118535542D-02 & 
     , 2.2800000000D+10, 3.1300000000D+15 & 
     , 1.0300000000D+14, 4.3200000000D+36 & 
     , 1.2500000000D+18, 1.0588306890D+15 & 
     , 2.1600000000D+29, 2.1000000000D+16 & 
     , 6.4700000000D-03, 1.2809740203D-03 & 
     , 1.1800000000D-03, 1.4005930548D-04 & 
     , 1.7700000000D+02, 3.8553336526D+02 & 
     , 2.1800000000D+04, 2.6600459910D+04 & 
     , 9.3300000000D+01, 5.4470805191D+03 & 
     , 4.2200000000D+11, 1.5531989112D+14 & 
     , 4.0900000000D+06, 5.8400000000D+09 & 
     , 2.8900000000D+05, 1.2000000000D+02 & 
     , 4.0900000000D+01, 3.4900000000D-11 /
      DATA (A(IH),IH=1621,1650)  /  2.3100000000D+03, 8.3778170592D-06 & 
     , 2.3100000000D+03, 2.9278172747D-03 & 
     , 7.8300000000D-01, 1.5357183416D-03 & 
     , 7.8300000000D-01, 4.3700000000D+15 & 
     , 4.9681868542D+18, 5.9900000000D+20 & 
     , 1.9700000000D+19, 2.9395633250D+13 & 
     , 2.0000000000D+10, 1.3800000000D-01 & 
     , 3.1358240429D-01, 1.1900000000D+06 & 
     , 4.2000000000D+16, 8.1800000000D+06 & 
     , 1.1700000000D+10, 5.7800000000D+05 & 
     , 2.4000000000D+02, 8.1800000000D+01 & 
     , 6.9800000000D-11, 4.6200000000D+03 & 
     , 6.0637611379D-06, 1.5700000000D+00 & 
     , 2.3100000000D+03, 2.6303314786D-05 & 
     , 7.8300000000D-01, 1.0967302857D-03 /
      DATA (A(IH),IH=1651,1680)  /  1.2500000000D+18, 7.6429253495D+12 & 
     , 3.2000000000D+34, 2.7785868184D+20 & 
     , 5.8300000000D+64, 8.2792870995D+55 & 
     , 8.2000000000D+14, 6.5241355800D+02 & 
     , 6.4700000000D-03, 9.2464157998D-06 & 
     , 1.1800000000D-03, 1.0109858237D-06 & 
     , 1.7700000000D+02, 2.7828837614D+00 & 
     , 2.1800000000D+04, 1.9200929050D+02 & 
     , 4.2200000000D+11, 1.1211408447D+12 & 
     , 9.3300000000D+01, 3.9318495595D+01 & 
     , 1.1000000000D+10, 1.4700000000D+10 & 
     , 2.2800000000D+11, 1.0156529597D+18 & 
     , 2.0000000000D+10, 1.1900000000D+06 & 
     , 7.3552497358D+06, 4.3200000000D+36 & 
     , 3.7600000000D+12, 1.7194981983D+14 /
      DATA (A(IH),IH=1681,1710)  /  3.0800000000D+06, 1.9592341791D+03 & 
     , 5.2600000000D+28, 2.7016744652D+26 & 
     , 7.2100000000D+33, 3.2171025483D+20 & 
     , 1.3300000000D+10, 1.5966820713D+10 & 
     , 6.6700000000D+09, 4.8005009689D+09 & 
     , 3.3300000000D+09, 4.3980858551D+10 & 
     , 2.8500000000D-16, 2.1086683447D-15 & 
     , 2.1000000000D+16, 4.0900000000D+06 & 
     , 5.8400000000D+09, 2.8900000000D+05 & 
     , 1.2000000000D+02, 4.0900000000D+01 & 
     , 3.4900000000D-11, 1.2900000000D+61 & 
     , 7.9583820932D+54, 1.0000000000D+84 & 
     , 4.1465511326D+78, 5.0000000000D+75 & 
     , 2.7802587455D+86, 6.0200000000D+05 & 
     , 8.6806271172D+01, 2.3400000000D+01 /
      DATA (A(IH),IH=1711,1740)  /  3.7121295904D-02, 1.3200000000D-01 & 
     , 5.8610887211D+04, 4.3400000000D+11 & 
     , 3.8569228355D+08, 2.7400000000D+05 & 
     , 2.4543429804D+08, 1.9900000000D+04 & 
     , 1.7825337704D+07, 2.0300000000D+09 & 
     , 1.2630305520D+17, 5.8000000000D+13 & 
     , 7.5519642443D+08, 4.0000000000D+10 & 
     , 3.4578688810D+06, 8.5700000000D+17 & 
     , 1.7748967608D+22, 3.0000000000D+10 & 
     , 5.1597641764D+17, 1.0000000000D+11 & 
     , 1.4519431559D+20, 3.0000000000D+10 & 
     , 3.1085468331D+17, 3.0000000000D+10 & 
     , 6.0438105008D+13, 3.8900000000D-06 & 
     , 1.4511636996D-05, 6.5900000000D+15 & 
     , 1.3791119702D+03, 1.0100000000D+71 /
      DATA (A(IH),IH=1741,1770)  /  1.4540798147D+66, 6.8300000000D-02 & 
     , 2.2983020763D-04, 1.7300000000D-02 & 
     , 6.4044991233D-04, 3.7000000000D-07 & 
     , 2.3145507889D-06, 7.3200000000D+10 & 
     , 1.5180780383D+09, 2.9000000000D+10 & 
     , 1.8351605777D-04, 6.3400000000D+10 & 
     , 5.2650873562D+13, 6.5100000000D+04 & 
     , 1.0805712956D+05, 7.4000000000D+11 & 
     , 9.7808623359D-02, 4.3000000000D+06 & 
     , 2.1563743598D-01, 3.0000000000D+10 & 
     , 1.8300000000D+10, 3.1309737000D+12 & 
     , 1.8300000000D+10, 3.1543081311D+17 & 
     , 2.2000000000D-01, 2.7062147246D+04 & 
     , 8.5700000000D+17, 2.7426874706D+21 & 
     , 8.5700000000D+17, 2.4175569754D+21 /
      DATA (A(IH),IH=1771,1800)  /  3.0000000000D+10, 3.0000000000D+10 & 
     , 1.0000000000D+11, 2.2436382722D+19 & 
     , 1.0000000000D+11, 1.9776673111D+19 & 
     , 3.0000000000D+10, 4.8035314725D+16 & 
     , 3.0000000000D+10, 4.2340992703D+16 & 
     , 8.6200000000D+15, 5.1809578901D+03 & 
     , 1.0100000000D+71, 1.0025283853D+66 & 
     , 6.8300000000D-02, 1.5845850043D-04 & 
     , 1.7300000000D-02, 4.4156394303D-04 & 
     , 3.7000000000D-07, 1.5957878251D-06 & 
     , 2.9000000000D+10, 1.8477233288D-02 & 
     , 1.6800000000D+11, 6.5100000000D+04 & 
     , 8.5700000000D+17, 8.5700000000D+17 & 
     , 8.5700000000D+17, 8.5700000000D+17 & 
     , 8.5700000000D+17, 1.1000000000D-01 /
      DATA (A(IH),IH=1801,1806)  /  1.1000000000D-01, 2.2000000000D-01 & 
     , 1.7600000000D-01, 2.2000000000D-01 & 
     , 2.2000000000D-01, 2.2000000000D-01 & 
     /
      DATA (N(IH),IH=1,30)  /  -0.67, -0.234291 & 
     , 2.7, 2.66385 & 
     , 1.51, 1.40625 & 
     , 2.4, 2.33241 & 
     , -1, -0.649765 & 
     , -0.6, -0.249765 & 
     , -1.25, -0.899765 & 
     , -2, -1.64977 & 
     , -2, -1.75351 & 
     , -1, -0.685917 & 
     , -1, -1.12163 & 
     , -1.4, -1.41184 & 
     , 2.43, 2.06793 & 
     , -0.58, -1.59836 & 
     , 0, 0.694036 /
      DATA (N(IH),IH=31,60)  /  0, 0.761631 & 
     , 0, 0.325921 & 
     , 0, 0.258327 & 
     , 0, 0.258327 & 
     , 0, -0.244887 & 
     , 0, -0.244887 & 
     , 2, 2.60696 & 
     , 0, 1.26484 & 
     , 2, 2.57081 & 
     , 0, 0.503214 & 
     , -7, -6.49679 & 
     , -2.79, -3.75335 & 
     , 0.14, -1.13743 & 
     , 0.03, -1.24743 & 
     , 0, -0.84172 /
      DATA (N(IH),IH=61,90)  /  0, -0.515799 & 
     , 0, 0.582947 & 
     , 0, 0.546795 & 
     , 0, -0.730634 & 
     , 0, 0.4792 & 
     , -1, -0.767288 & 
     , -1, -0.767288 & 
     , 0.81, 1.03087 & 
     , 0, -0.253456 & 
     , 0, 0.182253 & 
     , 0, 0.23442 & 
     , 0, -0.0551873 & 
     , 0, -0.601982 & 
     , 0, -0.292468 & 
     , -2.8, -3.25474 /
      DATA (N(IH),IH=91,120)  /  0, -0.947689 & 
     , 0, -0.166273 & 
     , -3.74, -5.18097 & 
     , 0, 0.675447 & 
     , -3.42, -4.10217 & 
     , -2.57, -2.66922 & 
     , -2.76, -2.92227 & 
     , 0, -0.345666 & 
     , 0, -0.758967 & 
     , 2, 2.18872 & 
     , 2, 1.4875 & 
     , 0, 0 & 
     , -0.323258, 0 & 
     , 0, 0.0026632 & 
     , 0, -1.00978 /
      DATA (N(IH),IH=121,150)  /  -5.11, -6.85567 & 
     , 0, -1.3719 & 
     , 0, -1.07943 & 
     , 0, 0 & 
     , -0.0622863, 0 & 
     , -0.0622863, 0 & 
     , 0.230182, 0 & 
     , 0.174995, 0 & 
     , -0.407952, 0 & 
     , -0.821254, 0 & 
     , -0.574789, 0 & 
     , 0.260469, 0 & 
     , 0.506957, -6.36 & 
     , -7.66989, 0 & 
     , -0.0622863, 0.25 /
      DATA (N(IH),IH=151,180)  /  0, -0.0622863 & 
     , 0, -0.0622863 & 
     , 0, 0.456176 & 
     , -4.82, -5.20137 & 
     , -4.8, -5.12229 & 
     , 1.9, 2.34945 & 
     , 0, 0.413302 & 
     , 1.18, 1.52571 & 
     , 0, 0.0873801 & 
     , 2, 1.84249 & 
     , 0, -1.35597 & 
     , -6.3, -6.46997 & 
     , 0, -0.282617 & 
     , 0, -5.92 & 
     , -6.75885, 1.6 /
      DATA (N(IH),IH=181,210)  /  2.00876, 0 & 
     , -1.34, -0.868957 & 
     , 0, -0.483283 & 
     , 0, 0.153092 & 
     , -4.89, -6.63941 & 
     , 0, 0.661222 & 
     , -1.61, 0 & 
     , 0, 0 & 
     , -0.157361, 0 & 
     , -0.158135, 2.47 & 
     , 2.55675, 0 & 
     , -1.09382, 0 & 
     , -1.03775, 0 & 
     , 0.0627383, 2.81 & 
     , 2.73924, 0 /
      DATA (N(IH),IH=211,240)  /  -1.23115, 0 & 
     , -1.29344, 0.1 & 
     , -1.0902, -7.44 & 
     , -7.35986, 1.63 & 
     , 1.57093, 0 & 
     , 0.672527, 0.5 & 
     , 1.41899, -0.23 & 
     , 1.16003, 0 & 
     , 0.636375, 0 & 
     , 0.568781, 7.6 & 
     , 7.91045, -4.65 & 
     , -4.51078, 0 & 
     , 0.731602, 0.65 & 
     , 1.62807, -0.09 & 
     , 1.35911, 0 /
      DATA (N(IH),IH=241,270)  /  0.69545, 0 & 
     , 0.627855, 0 & 
     , 0.369528, 1.62 & 
     , 2.14021, 1.5 & 
     , 1.98406, 1.6 & 
     , 2.01646, 0 & 
     , -1.00341, 2 & 
     , 2.00771, 0 & 
     , -0.0545804, 2.1 & 
     , 2.31102, 2.1 & 
     , 2.37009, 2.5 & 
     , 2.67486, 2.5 & 
     , 2.73394, 2 & 
     , 2.10727, 2 & 
     , 2.16634, 1.5 /
      DATA (N(IH),IH=271,300)  /  1.19081, 1.5 & 
     , 1.24988, -4.8 & 
     , -5.0463, 0 & 
     , 1.01264, 0 & 
     , -0.742405, 0 & 
     , 0.846372, 2.26 & 
     , 1.66346, 0 & 
     , 1.56102, 0 & 
     , 1.38578, 0 & 
     , 1.82149, 0 & 
     , 0.12683, 0 & 
     , 0.240947, 0 & 
     , 1.5678, -4.66 & 
     , -4.48812, 2 & 
     , 1.81798, 2 /
      DATA (N(IH),IH=301,330)  /  3.31672, 0 & 
     , -0.560387, 2.14 & 
     , 2.63279, 2 & 
     , 0.845959, 1.55 & 
     , 0.806965, 0.73 & 
     , 1.57037, 0 & 
     , 0.59717, 1.38 & 
     , 2.9634, 0 & 
     , 0.561018, 0 & 
     , 0.782321, 0 & 
     , 0.493423, 0 & 
     , 0.411006, -3.86 & 
     , -3.99563, 0 & 
     , 0.178352, 0.21 & 
     , -0.90671, 0 /
      DATA (N(IH),IH=331,360)  /  0.0746054, 1.61 & 
     , 1.42628, 0.29 & 
     , -0.705084, -1.39 & 
     , -0.801138, -6.57 & 
     , -6.05412, -6.87 & 
     , -4.77073, 0.03 & 
     , 1.61395, 0 & 
     , 0, 0 & 
     , 1.86656, 0 & 
     , 0.86611, 0 & 
     , 0.762364, 0 & 
     , 0.888497, 0 & 
     , -1.57811, 0 & 
     , 0.0256312, 1.16 & 
     , 1.22178, 1.16 /
      DATA (N(IH),IH=361,390)  /  0, 0 & 
     , 0.73, 0 & 
     , 1.77, -9.31 & 
     , -8.45483, -6.64 & 
     , -6.76131, 2.75 & 
     , 3.23587, 0.88 & 
     , -0.101077, 2.47 & 
     , 3.41854, 1.36 & 
     , 2.24549, 4.2 & 
     , 4.58212, -7.8 & 
     , -9.27438, 2 & 
     , 1.96566, -14.6 & 
     , -16.2667, -7.08 & 
     , -7.28894, 0 & 
     , 0.471548, 2.45 /
      DATA (N(IH),IH=391,420)  /  2.40134, 0.03 & 
     , -1.00898, -2.02 & 
     , -0.0734395, -0.69 & 
     , -0.222332, 0 & 
     , 0.455829, 1.02 & 
     , 1.12947, -15.74 & 
     , -13.9719, -13.14 & 
     , -11.7409, 2.7 & 
     , 3.25918, 3.8 & 
     , 4.32303, 2.22 & 
     , 2.67543, 0 & 
     , -0.0156091, 0 & 
     , 0.0389713, 0.03 & 
     , 1.29275, -13.55 & 
     , -13.7727, 0 /
      DATA (N(IH),IH=421,450)  /  0.382569, 0 & 
     , -0.0338928, -6.66 & 
     , -6.79608, -1.63 & 
     , -1.50576, 4.71 & 
     , 5.28298, 4.09 & 
     , 4.62683, 2.01 & 
     , 2.47923, 3.65 & 
     , 3.70277, 2.6 & 
     , 2.56602, -0.64 & 
     , -0.0988122, 0 & 
     , -0.0661922, 0 & 
     , 0.318776, 0 & 
     , -1.35613, 0.44 & 
     , 1.58056, -0.02 & 
     , 1.2253, -10.13 /
      DATA (N(IH),IH=451,480)  /  -11.5465, 0 & 
     , -1.29954, 0 & 
     , 0.204893, 0 & 
     , 0.121092, 0 & 
     , 0.0970785, 0 & 
     , -2.33737, 0 & 
     , -0.341857, -7.67 & 
     , -6.00376, -1.25 & 
     , -3.18571, -11.76 & 
     , -11.4905, 0 & 
     , -1.37788, 0 & 
     , 1.41602, 0 & 
     , 1.36084, 0 & 
     , 1.04675, 0 & 
     , 1.48246, 0 /
      DATA (N(IH),IH=481,510)  /  0.267581, 0 & 
     , 0.27143, 0 & 
     , 0.606272, 0 & 
     , 0.605334, 0 & 
     , 0.352248, 0 & 
     , 0.084667, 0 & 
     , 0.393202, 0 & 
     , 0.0769611, -0.5 & 
     , 0.176365, 0 & 
     , -1.47857, 0 & 
     , 0.123795, 0.42 & 
     , 0, 0.0237669 & 
     , 0, -0.197106 & 
     , 0, -0.135412 & 
     , 0, -0.591432 /
      DATA (N(IH),IH=511,540)  /  3.37, 3.32222 & 
     , 0, -0.483536 & 
     , 0, 0.514343 & 
     , 0, 1.00075 & 
     , 0, -1.19135 & 
     , 0, -1.30881 & 
     , 0, -1.28933 & 
     , 0, 0.18992 & 
     , 0, -1.55827 & 
     , 0.44, 1.75196 & 
     , 1.16, 0 & 
     , 0, 0.73 & 
     , 0, 1.77 & 
     , -3.86, -3.97975 & 
     , 1.13, 1.59998 /
      DATA (N(IH),IH=541,570)  /  -5.06, -5.06898 & 
     , -5, -5.11373 & 
     , 1.55, 0.0186034 & 
     , 0.73, 1.63305 & 
     , 2.28, 2.64624 & 
     , -7.8, -7.26513 & 
     , 0, -0.363787 & 
     , 1.7, 2.34085 & 
     , 0, 1.80032 & 
     , 0, -0.101893 & 
     , 0, 0.00285485 & 
     , 0, 2.23323 & 
     , 0, 0.11898 & 
     , 0, 0.223728 & 
     , 0, -1.23168 /
      DATA (N(IH),IH=571,600)  /  0, -1.33185 & 
     , -7.8, -7.69525 & 
     , -0.33, -0.225252 & 
     , -10.77, -11.2885 & 
     , -12.46, -12.5238 & 
     , -12.36, -12.5285 & 
     , -11.97, -12.246 & 
     , 2.7, 3.05922 & 
     , 1.92, 2.24307 & 
     , 2.22, 2.47547 & 
     , 0, -0.16099 & 
     , 3.37, 3.12226 & 
     , 2.53, 2.99397 & 
     , 4.2, 4.56022 & 
     , 2, 1.94376 /
      DATA (N(IH),IH=601,630)  /  0.12, -0.0229938 & 
     , 2.05, 3.02862 & 
     , 2, 2.95854 & 
     , 2, 3.22612 & 
     , 0, -0.132572 & 
     , 0, -0.237321 & 
     , 2, 1.98652 & 
     , 1.55, 1.94752 & 
     , 0.73, 1.52072 & 
     , 1.16, 0 & 
     , 1.35, 2.5 & 
     , 6.21, -14.08 & 
     , -13.6253, -13.02 & 
     , -12.6727, -12.71 & 
     , -12.8174, 2.8 /
      DATA (N(IH),IH=631,660)  /  3.66875, 0 & 
     , 0.765006, -0.32 & 
     , 0.0285441, 0 & 
     , 0.382884, 0 & 
     , 0.309573, 0 & 
     , 0.787987, 2.19 & 
     , 2.69668, 1.81 & 
     , 1.66065, -2.7 & 
     , 0.37, 1.028 & 
     , 0.03, -0.333429 & 
     , -5.16, -5.70891 & 
     , 0, 0.502181 & 
     , 0, 0.281308 & 
     , -0.94, -0.855882 & 
     , 0, 0.902702 /
      DATA (N(IH),IH=661,690)  /  0, 0.518782 & 
     , 0, 0.880152 & 
     , 0, 0 & 
     , 1.3277, 0 & 
     , 0.0474627, 0 & 
     , -0.00142624, 0 & 
     , 0.62621, 0 & 
     , 1.26938, 0 & 
     , 0, 1.71693 & 
     , 0, 0.15489 & 
     , 0, 0.106001 & 
     , 1.61, 1.76671 & 
     , 1.61, 1.87414 & 
     , 0.29, 0.29 & 
     , -1.39, -1.39 /
      DATA (N(IH),IH=691,720)  /  -0.691616, 1.02 & 
     , 1.07196, 0 & 
     , 0, -0.221635 & 
     , 0.09, -2.39 & 
     , -0.859397, 1.6 & 
     , 2.14953, 1.6 & 
     , 2.82589, 2.54 & 
     , 2.62077, 2.68 & 
     , 2.72461, 1.46 & 
     , 1.43702, 2.6 & 
     , 2.0738, 3.65 & 
     , 3.21056, 2.5 & 
     , 3.03548, 0.7 & 
     , 1.19933, 2 & 
     , 2.43174, 3.5 /
      DATA (N(IH),IH=721,750)  /  3.51528, 2.53 & 
     , 2.95806, 0.7 & 
     , 1.0919, 4.2 & 
     , 4.52431, 2 & 
     , 1.90785, 0 & 
     , 1.77683, 0 & 
     , -0.788396, 1.8 & 
     , 2.93863, 0 & 
     , -0.552714, 0 & 
     , 0, -0.924094 & 
     , 0, 0.751128 & 
     , 2, 3.13616 & 
     , -8.1, -8.26008 & 
     , -7.87, -7.9297 & 
     , 1.03, 2.06488 /
      DATA (N(IH),IH=751,780)  /  -0.25, 0.452329 & 
     , -9.5, -9.60039 & 
     , -1.67, -1.77039 & 
     , -10.26, -10.4028 & 
     , -9.01, -9.0524 & 
     , -3.34, -2.28677 & 
     , -2.55, -1.39638 & 
     , 0, 0.409933 & 
     , 0, 0.510318 & 
     , 0, 0.306186 & 
     , 0, 0.406572 & 
     , 1.61, 1.65786 & 
     , 1.61, 1.75825 & 
     , 0, 0.970973 & 
     , 0, 0.664277 /
      DATA (N(IH),IH=781,810)  /  1.7, 3.9463 & 
     , 2.75, 3.24302 & 
     , 2.75, 3.14263 & 
     , 4.2, 4.58927 & 
     , 4.2, 4.48888 & 
     , 2, 1.97281 & 
     , 2, 1.87242 & 
     , 2, 3.10516 & 
     , 2.47, 3.47671 & 
     , 1.36, 2.34618 & 
     , 0, -1.42196 & 
     , 0, -1.96405 & 
     , 0, -1.33851 & 
     , 0, -2.06272 & 
     , 0, -3.01231 /
      DATA (N(IH),IH=811,840)  /  0, -1.7372 & 
     , 0.16, -1.51022 & 
     , 0, -2.77864 & 
     , -2.44, -4.62745 & 
     , -2.04, -3.8821 & 
     , 0, -0.307517 & 
     , -11.94, -14.3484 & 
     , 0, -1.1124 & 
     , -6.27, -5.67882 & 
     , -8.62, -7.68347 & 
     , 0, 1.58909 & 
     , 1.14, -1.25945 & 
     , 1.06, -1.2347 & 
     , 2.53, 3.81677 & 
     , 2.53, 3.47142 /
      DATA (N(IH),IH=841,870)  /  0, -1.25062 & 
     , 1.9, 2.80526 & 
     , 2, 3.18302 & 
     , 2, 2.83767 & 
     , 0, 0.766561 & 
     , 0, 0.421208 & 
     , 0, 0.800901 & 
     , 0, 0.455548 & 
     , 0.88, 2.47 & 
     , 4.58683, 1.36 & 
     , 2.95278, -7.8 & 
     , -6.62052, -11.92 & 
     , -12.2223, -11.92 & 
     , -12.5677, -16.89 & 
     , -17.2354, -3.35 /
      DATA (N(IH),IH=871,900)  /  -3.69535, 0 & 
     , 0.652554, 0 & 
     , 0.548808, 0 & 
     , -0.703823, 0 & 
     , -0.679809, 0 & 
     , -0.924696, 0.21 & 
     , 1.09896, 1.61 & 
     , 1.90048, 0.29 & 
     , -1.39, -1.19719 & 
     , 0, 0.997907 & 
     , 0, 1.81725 & 
     , 0, 0.894161 & 
     , 0, -0.35847 & 
     , 0, -0.579343 & 
     , 0, -0.334456 /
      DATA (N(IH),IH=901,930)  /  0, 0.591536 & 
     , 0.03, 1.56463 & 
     , 7.07, 6.131 & 
     , -17.62, -17.71 & 
     , 2.6, 2.81 & 
     , 1.5, 0.2 & 
     , 0.21, -19.78 & 
     , 0.98, 0.02 & 
     , -1.75, -1.65 & 
     , -1.18, -0.53 & 
     , 0.45, 0 & 
     , 0, -1.36 & 
     , -0.7, -2.03 & 
     , -1.87, 3.77 & 
     , 0, 0.278353 /
      DATA (N(IH),IH=931,960)  /  0, -3.62 & 
     , -1.65, 0 & 
     , 0, -2.61 & 
     , -0.47, 0.12 & 
     , 8.84, -1.63 & 
     , 0, 2.56 & 
     , 2, 0 & 
     , -0.52, -1.4 & 
     , 0.11, -0.64 & 
     , -5.43, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 2.34 /
      DATA (N(IH),IH=961,990)  /  0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0.76, 0.5 & 
     , 0, 0 & 
     , -5.9, -6.66828 & 
     , 1.84, -0.0347199 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, 0 & 
     , 0, -0.87 & 
     , 0, 0 & 
     , 0, 0 /
      DATA (N(IH),IH=991,1020)  /  -9.38, -9.66 & 
     , -11.84, 3.92 & 
     , 3.07, 1.99 & 
     , 0.84, 9.25 & 
     , 0.73, 0.73 & 
     , -2.81, -7.47 & 
     , -3.34, 0 & 
     , 0, -2.98 & 
     , -6, -2.08 & 
     , -1.01, -0.9 & 
     , -1.49, 15.73 & 
     , -12.99, -14.94 & 
     , 7.67, 2.36 & 
     , 6.18, 2.76 & 
     , 0, -7.3 /
      DATA (N(IH),IH=1021,1050)  /  -7.59, 0.03 & 
     , -17.29, -9.64 & 
     , -8.43, 0 & 
     , 0, -14.22 & 
     , -11.49, -11.8109 & 
     , -5.72, 2.5 & 
     , 2.52929, 0.7 & 
     , 2, 3.5 & 
     , 2.6, 2.6 & 
     , 1.76, 1.76 & 
     , -1.45, 0 & 
     , 0, -9.74 & 
     , 0, 0 & 
     , -0.45, -1.21 & 
     , -5.71, 0 /
      DATA (N(IH),IH=1051,1080)  /  0, 0 & 
     , 0.76, 0 & 
     , 0, 5.36 & 
     , -1.59, 0.03 & 
     , 0, -0.06 & 
     , 2, 0 & 
     , 0.97, 0 & 
     , 0.5, 0.794062 & 
     , -0.89, -2.30141 & 
     , -3.45, -4.51606 & 
     , -8.9, -8.60594 & 
     , -14.65, -14.1444 & 
     , -17.77, -19.121 & 
     , -1, -2.11735 & 
     , -3.3, -4.07199 /
      DATA (N(IH),IH=1081,1110)  /  -7.74, -10.1 & 
     , -11.8517, -9.57 & 
     , -11.0277, -7 & 
     , -8.24617, 2.05 & 
     , 0.27021, 0 & 
     , -2.06437, 1.47 & 
     , -0.382874, 4.47 & 
     , 4.39414, -12.48 & 
     , -12.2493, 1.89 & 
     , 2.47097, 2.68 & 
     , 3.15722, 0.54 & 
     , 0.647444, 1.02 & 
     , 1.32078, 0.49 & 
     , 1.14102, 0 & 
     , 0.547273, -12.48 /
      DATA (N(IH),IH=1111,1140)  /  -12.2599, 2.36 & 
     , 2.93038, 3.33 & 
     , 3.79663, 0.34 & 
     , 0.452698, 2.75 & 
     , 3.21293, 4.2 & 
     , 4.55919, 2.5 & 
     , 0.979368, 2.55 & 
     , 1.2202, 1.62 & 
     , 0.397642, -1.44 & 
     , -3.13503, -4.24 & 
     , -6.03013, 2.61 & 
     , 0.803868, -12.48 & 
     , -12.2468, 1.87 & 
     , 2.45348, 3.02 & 
     , 3.49973, -12.48 /
      DATA (N(IH),IH=1141,1170)  /  -12.2625, 1.87 & 
     , 2.43769, 3.02 & 
     , 3.48394, 2.05 & 
     , 0.192174, 2.05 & 
     , 0.287919, 0 & 
     , -2.02971, 0 & 
     , -1.93396, -1.44 & 
     , -2.88623, -1.44 & 
     , -2.80628, -4.24 & 
     , -5.78384, -4.24 & 
     , -5.6881, 1.02 & 
     , 1.38384, 0.49 & 
     , 1.20408, 0 & 
     , 0.610331, 1.02 & 
     , 1.29787, 0.49 /
      DATA (N(IH),IH=1171,1200)  /  1.11811, 0 & 
     , 0.524359, -12.48 & 
     , -12.2367, 1.88 & 
     , 2.47349, 3.33 & 
     , 3.81974, -12.48 & 
     , -12.225, 1.88 & 
     , 2.48528, 3.33 & 
     , 3.83153, 0.23 & 
     , 0.811179, 0 & 
     , 0.217336, -12.48 & 
     , -12.2259, 1.87 & 
     , 2.47429, 3.02 & 
     , 3.52054, 2.05 & 
     , -0.322738, 0 & 
     , -2.54462, -1.44 /
      DATA (N(IH),IH=1201,1230)  /  -3.38034, -4.24 & 
     , -6.29875, -12.48 & 
     , -12.2371, 1.88 & 
     , 2.4731, 3.33 & 
     , 3.81935, 1.02 & 
     , 1.88558, 0.49 & 
     , 1.70581, 0 & 
     , 1.11207, 0 & 
     , -1.57654, 0 & 
     , -1.78804, -12.48 & 
     , -12.4557, 1.8 & 
     , 2.17451, 3.33 & 
     , 3.60077, 2.5 & 
     , 0.97143, 2.5 & 
     , 0.96566, 1.62 /
      DATA (N(IH),IH=1231,1260)  /  0.473091, 1.62 & 
     , 0.393137, 2.05 & 
     , 0.444306, 0 & 
     , -1.92589, 0 & 
     , -1.93166, -1.44 & 
     , -3.12264, -1.44 & 
     , -3.11662, -4.24 & 
     , -6.03025, -4.24 & 
     , -6.03603, 2.61 & 
     , 0.801281, 2.61 & 
     , 0.817072, 0 & 
     , -1.89045, 0 & 
     , -1.90968, 0 & 
     , -2.12118, -12.48 & 
     , -12.2546, 1.8 /
      DATA (N(IH),IH=1261,1290)  /  2.37567, 3.33 & 
     , 3.80192, -12.48 & 
     , -12.2546, 1.87 & 
     , 2.44567, 3.02 & 
     , 3.49192, 1.08 & 
     , 2.8367, 2.5 & 
     , 0.969451, 1.62 & 
     , 0.969942, 2.05 & 
     , 0.783907, 0 & 
     , -1.08774, 2.61 & 
     , 0.78322, -12.48 & 
     , -12.252, 1.87 & 
     , 2.44819, 3.02 & 
     , 3.49444, 2.05 & 
     , 0.882089, -12.48 /
      DATA (N(IH),IH=1291,1320)  /  -12.2341, 1.87 & 
     , 2.46617, 3.02 & 
     , 3.51242, 2.05 & 
     , 0.76659, 2.05 & 
     , 0.882751, 0 & 
     , 0, 0 & 
     , -15.16, -14.464 & 
     , 0, 1.04627 & 
     , 0, 1.13526 & 
     , 0, 2.71 & 
     , 3.72012, 2 & 
     , 2.94252, 0 & 
     , 0.684196, 2.6 & 
     , 3.03931, 4 & 
     , 4.52606, 0 /
      DATA (N(IH),IH=1321,1350)  /  0.560401, 0 & 
     , -0.2405, 0.88 & 
     , -0.123996, 1.36 & 
     , -7.8, -7.8 & 
     , -12.5, -13.5715 & 
     , -4.03, -4.85 & 
     , 0, -1.03629 & 
     , 1.3, 0.699415 & 
     , 0, -0.692564 & 
     , 0, 0 & 
     , 0.103818, 0 & 
     , 0, 0 & 
     , 0, 7.6 & 
     , 0, 1.02 & 
     , -8, 1.46 /
      DATA (N(IH),IH=1351,1380)  /  0.796263, 0 & 
     , 1.1293, 3.07 & 
     , 0, -2.03247 & 
     , -15.16, -14.9678 & 
     , 0, 0.542404 & 
     , 2.48, 1.36023 & 
     , 2.71, 3.21625 & 
     , 2, 2.43866 & 
     , 0, 0.180331 & 
     , 2.6, 2.53544 & 
     , 4, 4.0222 & 
     , 0.88, -7.8 & 
     , -4.03, -4.85 & 
     , 0, -0.525133 & 
     , 1.3, 1.21058 /
      DATA (N(IH),IH=1381,1410)  /  0, 0 & 
     , -7.74, 0 & 
     , 0, -8 & 
     , 1.46, 2.17 & 
     , 3.57723, 0.68 & 
     , 0.254459, -1.73 & 
     , -0.111273, -14.15 & 
     , -12.1057, 0 & 
     , 2.2454, 2.5 & 
     , 2.06262, 3.98 & 
     , 3.90469, 2.39 & 
     , 2.21095, 2.88 & 
     , 2.76966, 3.37 & 
     , 1.86503, 4.09 & 
     , 3.97854, 0.3 /
      DATA (N(IH),IH=1411,1440)  /  -0.890885, 1.8 & 
     , 0.788084, 0 & 
     , -0.595515, 2.5 & 
     , 1.81773, 0 & 
     , -0.637038, 0 & 
     , -0.658482, 0 & 
     , -0.492206, 1.03 & 
     , 0.819065, -7.74 & 
     , -1.55, -1.51538 & 
     , -8.86, -8.19307 & 
     , -7.44, -6.95964 & 
     , 2.1, 1.96988 & 
     , 2.5, 2.33372 & 
     , 2, 1.76613 & 
     , 1.5, 0.849667 /
      DATA (N(IH),IH=1441,1470)  /  2.19, 2.07772 & 
     , -5.08, -4.82261 & 
     , -6.21, -3.78987 & 
     , -6.1, -3.79215 & 
     , 0, 0.607624 & 
     , 0, 0.571472 & 
     , 0, 0.503877 & 
     , 7.6, 7.84555 & 
     , 0, 1.16 & 
     , 0, 1.35 & 
     , 2.5, 2.5 & 
     , 6.21, -15.92 & 
     , -15.741, 0 & 
     , 0.529204, 0 & 
     , 0.493052, 2 /
      DATA (N(IH),IH=1471,1500)  /  2.42546, 0 & 
     , -12.48, -12.2592 & 
     , 1.8, 2.37108 & 
     , 0, 0.534927 & 
     , 2.69, 3.15733 & 
     , 4.46, 4.51087 & 
     , 0, -1.23276 & 
     , 0, -1.54684 & 
     , 0, -0.785212 & 
     , -2.27, -3.38113 & 
     , -0.44, -0.44 & 
     , 0, -0.723732 & 
     , 0, -1.4418 & 
     , 0, -2.29587 & 
     , 2.17, 4.25437 /
      DATA (N(IH),IH=1501,1530)  /  2.88, 3.4468 & 
     , 4.71, 5.78397 & 
     , 4.09, 5.12782 & 
     , 2.01, 2.98022 & 
     , 2.6, 3.06701 & 
     , 3.65, 4.20376 & 
     , 0.78, 2.47345 & 
     , 2.08, 1.84471 & 
     , 0, 0.114944 & 
     , 0, 0.0111972 & 
     , 0, -0.405265 & 
     , -1.63, -1.87713 & 
     , 0.03, 0.465447 & 
     , 0.03, 0.34697 & 
     , 0, -8.23 /
      DATA (N(IH),IH=1531,1560)  /  1.84, -10.15 & 
     , -4.48, -1.12 & 
     , 3.57, 0 & 
     , 0, 0 & 
     , 1.66, 2.42996 & 
     , 2, 1.94272 & 
     , -7.8, -7.56182 & 
     , -7.8, -7.44334 & 
     , 2.88, 3.2153 & 
     , 0, 0.534226 & 
     , 0.21, 1.09997 & 
     , 1.61, 1.89895 & 
     , 0.29, -1.39 & 
     , -1.01221, 0 & 
     , 0.771407, -0.6 /
      DATA (N(IH),IH=1561,1590)  /  -1.02266, -3.58 & 
     , -1.95607, 3.98 & 
     , 3.90757, 4.09 & 
     , 3.98142, 2.39 & 
     , 2.21383, 2.5 & 
     , 2.0655, 2.5 & 
     , 1.82061, 0 & 
     , -0.592634, 2.17 & 
     , 3.57309, 2.88 & 
     , 2.77812, 1.55 & 
     , 0, -14.15 & 
     , -12.1034, -1.34 & 
     , -1.74072, -2.47 & 
     , -0.706024, -2.36 & 
     , -0.71765, 0 /
      DATA (N(IH),IH=1591,1620)  /  2.42, 2.45499 & 
     , 1.87, 2.52284 & 
     , -0.31, -1.44 & 
     , -1.33, -7.74 & 
     , -0.6, -1.02766 & 
     , -3.58, 0 & 
     , 3.98, 3.90257 & 
     , 4.09, 3.97642 & 
     , 2.39, 2.20882 & 
     , 2.5, 2.0605 & 
     , 2.5, 1.81561 & 
     , 0, -0.597639 & 
     , 1.16, 0 & 
     , 1.35, 2.5 & 
     , 2.5, 6.21 /
      DATA (N(IH),IH=1621,1650)  /  2.17, 4.21307 & 
     , 2.17, 3.56984 & 
     , 2.88, 3.41811 & 
     , 2.88, -1.34 & 
     , -1.77556, -2.47 & 
     , -2.36, -0.715895 & 
     , 0, 2.42 & 
     , 2.42015, 1.03 & 
     , 0, 1.16 & 
     , 0, 1.35 & 
     , 2.5, 2.5 & 
     , 6.21, 2.17 & 
     , 4.24966, 2.88 & 
     , 2.17, 4.14699 & 
     , 2.88, 3.74052 /
      DATA (N(IH),IH=1651,1680)  /  -0.6, -0.435322 & 
     , -5.02, -2.80977 & 
     , -14.15, -12.1044 & 
     , 0, 1.84259 & 
     , 3.98, 4.49491 & 
     , 4.09, 4.56876 & 
     , 2.39, 2.80117 & 
     , 2.5, 2.65284 & 
     , 0, -0.00529554 & 
     , 2.5, 2.40795 & 
     , 0, 0 & 
     , 0, -0.655564 & 
     , 0, 1.03 & 
     , 0.821983, -7.74 & 
     , -1.55, -1.93634 /
      DATA (N(IH),IH=1681,1710)  /  0.37, 1.34013 & 
     , -5.08, -5.24648 & 
     , -6.21, -3.7915 & 
     , 0, 0.183753 & 
     , 0, 0.147601 & 
     , 0, 0.0800063 & 
     , 7.6, 7.42168 & 
     , 0, 1.16 & 
     , 0, 1.35 & 
     , 2.5, 2.5 & 
     , 6.21, -12.48 & 
     , -12.2685, -18.87 & 
     , -18.4977, -19.31 & 
     , -20.3484, 1.8 & 
     , 2.36173, 2.68 /
      DATA (N(IH),IH=1711,1740)  /  3.13799, 3.25 & 
     , 1.73243, 0 & 
     , 0.199659, 1.55 & 
     , 0.545147, 1.8 & 
     , 0.795147, 0.3 & 
     , -0.903488, -0.77 & 
     , 0.101044, 0 & 
     , 0.525581, -2.27 & 
     , -3.36472, 0 & 
     , -1.69626, 0 & 
     , -1.21635, 0 & 
     , -1.53043, 0 & 
     , -0.768803, 4.57 & 
     , 4.52848, -0.61 & 
     , 1.46453, -15.92 /
      DATA (N(IH),IH=1741,1770)  /  -15.7214, 3.4 & 
     , 3.94887, 3.4 & 
     , 3.84512, 4.7 & 
     , 4.72866, 0 & 
     , 0.186797, 0 & 
     , 2.57193, 0.03 & 
     , -0.571534, 1.3 & 
     , 1.13417, 0 & 
     , 2.13717, 1.45 & 
     , 2.92344, 0 & 
     , 0, -0.681196 & 
     , 0, -0.802389 & 
     , 3.25, 2.13353 & 
     , -2.27, -3.06281 & 
     , -2.27, -3.04702 /
      DATA (N(IH),IH=1771,1800)  /  0, 0 & 
     , 0, -0.91444 & 
     , 0, -0.898649 & 
     , 0, -1.22852 & 
     , 0, -1.21273 & 
     , -0.61, 1.18494 & 
     , -15.92, -15.7988 & 
     , 3.4, 3.87143 & 
     , 3.4, 3.76768 & 
     , 4.7, 4.65122 & 
     , 0, 1.86592 & 
     , 0, 1.3 & 
     , -2.27, -2.27 & 
     , -2.27, -2.27 & 
     , -2.27, 3.25 /
      DATA (N(IH),IH=1801,1806)  /  3.25, 3.25 & 
     , 3.25, 3.25 & 
     , 3.25, 3.25 & 
     /
      DATA (E(IH),IH=1,30)  /  71300000, 479604.2721 & 
      , 26190000, 20100582.44 & 
      , 14350000, 77394524.74 & 
      , -8830000, 60303942.3 & 
      , 0, 433091487.8 & 
      , 0, 433091487.8 & 
      , 0, 433091487.8 & 
      , 0, 433091487.8 & 
      , 0, 496136012.5 & 
      , 0, 427002070.2 & 
      , 0, 497822466 & 
      , 0, 204667815.3 & 
      , 223850000, -4573672.487 & 
      , -9590000, 205148485.1 & 
      , 2810000, 223457801.5 /
      DATA (E(IH),IH=31,60)  /  2660000, 154173859.2 & 
      , 0, 222334254.9 & 
      , -2090000, 289378197.2 & 
      , 72510000, 363978197.2 & 
      , -6820000, 154764529 & 
      , 50210000, 211794529 & 
      , 21760000, 88599143.48 & 
      , 16610000, 298007527.4 & 
      , 16610000, 77359725.92 & 
      , 1790000, 131673668.2 & 
      , 157320000, 287203668.2 & 
      , 17540000, 552243890 & 
      , 30760000, 138461819.8 & 
      , -70000, 107631819.8 & 
      , 199580000, 236461424.1 /
      DATA (E(IH),IH=61,90)  /  96230000, 355445679 & 
      , 0, 367317915.9 & 
      , 0, 361228498.4 & 
      , 0, 468930318.2 & 
      , 0, 430362440.7 & 
      , 71130000, 5356428.13 & 
      , 71130000, 5356428.13 & 
      , -3040000, 135854243.4 & 
      , 0, 648064893.1 & 
      , 2410000, 579654497.4 & 
      , 0, 97396261.59 & 
      , 0, 739371737.2 & 
      , 0, 378143238.8 & 
      , 13010000, 1691026.963 & 
      , 2470000, 455400199.9 /
      DATA (E(IH),IH=91,120)  /  -3160000, 248452632.8 & 
      , 0, 307322843.1 & 
      , 8100000, 320775841.3 & 
      , 66070000, 336511419 & 
      , 352920000, 355207490.6 & 
      , 1780000, 371385406.5 & 
      , 6690000, 470939173 & 
      , 0, 383372794.3 & 
      , 0, 325976130.6 & 
      , 12550000, 86913497.78 & 
      , 30250000, 61407685.18 & 
      , 6280000, 6280000 & 
      , 261435734.9, 6280000 & 
      , 0, 477489989.8 & 
      , 0, 327673077.5 /
      DATA (E(IH),IH=121,150)  /  29690000, 364655095.4 & 
      , 0, 548460717.8 & 
      , 49970000, 609749690.8 & 
      , 45980000, 2510000 & 
      , 39812987.98, 2510000 & 
      , 39812987.98, 0 & 
      , 48621961.01, 0 & 
      , 787993698.2, 0 & 
      , 420675782.2, 0 & 
      , 363279118.6, 0 & 
      , 68460673.16, 0 & 
      , 284081814.6, 0 & 
      , 780217827.2, 21090000 & 
      , 416708647.8, 0 & 
      , 37302987.98, -3910000 /
      DATA (E(IH),IH=151,180)  /  0, 37302987.98 & 
      , 0, 37302987.98 & 
      , 0, 255577298.7 & 
      , 27320000, 150083586.5 & 
      , 23260000, 117496132.5 & 
      , 11470000, 74956081.25 & 
      , 14810000, 72206663.69 & 
      , -1870000, 124660606 & 
      , 167360000, 2422408.767 & 
      , 50210000, 46856937.77 & 
      , -2160000, 319198631.8 & 
      , 21230000, 463233395.7 & 
      , 0, 288729027.8 & 
      , 0, 13140000 & 
      , 403342499.3, 22680000 /
      DATA (E(IH),IH=181,210)  /  54566839.56, -7340000 & 
      , 5930000, 513851.5818 & 
      , 127700000, 12842694.35 & 
      , 57910000, 275818632.1 & 
      , 14360000, 158160807.4 & 
      , -5020000, 119287047.3 & 
      , 7780000, -6570000 & 
      , 48830000, 0 & 
      , 107476949.3, 0 & 
      , 237335580.4, 21670000 & 
      , 97421051.39, 0 & 
      , 419906771, 0 & 
      , 231012782.6, 0 & 
      , 376229823.8, 24520000 & 
      , 96917989.16, 0 /
      DATA (E(IH),IH=211,240)  /  275759524.7, -2390000 & 
      , 310672512.6, 44350000 & 
      , 7197984.006, 58910000 & 
      , 493149409.3, 8050000 & 
      , 36577454.01, 0 & 
      , 338855355.3, -460000 & 
      , 43576909.93, 4480000 & 
      , 43100761.51, 0 & 
      , 332765937.8, 0 & 
      , 401899880.1, -14770000 & 
      , 95661682.83, 21260000 & 
      , 426971955.3, 0 & 
      , 310327901.3, -1190000 & 
      , 14319455.91, 2550000 & 
      , 12643307.5, 0 /
      DATA (E(IH),IH=241,270)  /  304238483.7, 0 & 
      , 373372426.1, 3770000 & 
      , 85674228.82, 45360000 & 
      , 36448092.09, 35980000 & 
      , 20978674.53, 13050000 & 
      , 67182616.83, 0 & 
      , 255528643.7, 34600000 & 
      , 56845777.28, -2390000 & 
      , 57158765.25, 20380000 & 
      , 47759532.54, 20380000 & 
      , 19232078.53, 12970000 & 
      , 34260114.98, 20920000 & 
      , 13682660.97, -3520000 & 
      , 86904057.28, 6280000 & 
      , 68176603.27, 41590000 /
      DATA (E(IH),IH=271,300)  /  77881440.45, 41590000 & 
      , 49353986.44, 7950000 & 
      , 564432866.5, 0 & 
      , 325621371.1, 0 & 
      , 211295142.2, -3160000 & 
      , 629784214.1, 3770000 & 
      , 127161378.7, 0 & 
      , 71793685.47, 0 & 
      , 426695895.8, 3570000 & 
      , 359445500.1, 0 & 
      , 657557391.2, 0 & 
      , 382586114.2, 0 & 
      , 344881549.9, 15820000 & 
      , 162621237.8, 7950000 & 
      , 89764345.93, 7950000 /
      DATA (E(IH),IH=301,330)  /  198861019.4, 0 & 
      , 129480796.2, 71380000 & 
      , 11033146.07, 53190000 & 
      , 26274495.32, 8810000 & 
      , 107684044.6, 10790000 & 
      , 238948122.1, 33470000 & 
      , 22499718.93, 2570000 & 
      , 131854077.6, 33470000 & 
      , 16410301.37, 5650000 & 
      , 205388794.6, 8370000 & 
      , 60444243.67, 0 & 
      , 125789549.2, 13890000 & 
      , 480409256.8, 0 & 
      , 286290250, -1790000 & 
      , 531013741.4, 0 /
      DATA (E(IH),IH=331,360)  /  349334774.8, -1610000 & 
      , 56256577.53, 50000 & 
      , 35031275.47, 4250000 & 
      , 369289088.3, 206930000 & 
      , 53201135.59, 197460000 & 
      , 173015213.1, -1650000 & 
      , 328407812.9, 0 & 
      , 0, 0 & 
      , 41328785.02, 0 & 
      , 279362623.4, 0 & 
      , 342407148.1, 0 & 
      , 25819329.1, 0 & 
      , 365339262.7, 7570000 & 
      , 27904022.48, 10060000 & 
      , 36483440.04, 10060000 /
      DATA (E(IH),IH=361,390)  /  7570000, 163800000 & 
      , -4660000, 49890000 & 
      , 24770000, 417980000 & 
      , 49499564.99, 24140000 & 
      , 175477632.3, 48740000 & 
      , 15312230.95, 4770000 & 
      , 71054484.59, 3890000 & 
      , 16859503.17, 3710000 & 
      , 111323269.6, -3600000 & 
      , 26016755.69, 29540000 & 
      , 141884176.2, 38490000 & 
      , 13974138.86, 76020000 & 
      , 179210339.7, 27970000 & 
      , 449644917.7, 0 & 
      , 281753855.5, 12220000 /
      DATA (E(IH),IH=391,420)  /  302885763.4, -1650000 & 
      , 386358614.1, 86820000 & 
      , 24692429.7, 93010000 & 
      , 26616285.89, 3660000 & 
      , 141934101.2, -8510000 & 
      , 44820182.98, 413040000 & 
      , 37586111.82, 425010000 & 
      , 40487098.26, 24020000 & 
      , 35436570.06, 13100000 & 
      , 18427152.49, 3100000 & 
      , 77561094.8, -2300000 & 
      , 77577243.21, 93120000 & 
      , 113448478, -1650000 & 
      , 335226320.5, 47520000 & 
      , 471121180.8, 0 /
      DATA (E(IH),IH=421,450)  /  359404694.1, -3220000 & 
      , 302052077.3, 29290000 & 
      , 166021318.4, 14300000 & 
      , 82236496.87, 25990000 & 
      , 35480306.99, 10650000 & 
      , 14050889.43, 1530000 & 
      , 74064831.73, 29930000 & 
      , 48332214.89, 58200000 & 
      , 851163.5039, 207940000 & 
      , 19688571.81, 0 & 
      , 416043871.6, 0 & 
      , 330486347.7, 0 & 
      , 105055737.5, 22860000 & 
      , 46357562.82, 47070000 & 
      , 75469260.17, 77500000 /
      DATA (E(IH),IH=451,480)  /  186098173.3, 0 & 
      , 250786354.2, 12550000 & 
      , 189399422.9, -2490000 & 
      , 97776912.53, 0 & 
      , 400745685, 0 & 
      , 433035302.1, -3200000 & 
      , 292002157.9, 467900000 & 
      , 34921721.89, 32090000 & 
      , 101188124, 98530000 & 
      , 462410154.1, 0 & 
      , 159986443, 0 & 
      , 18103989.21, 0 & 
      , 757475726.4, 0 & 
      , 330473656.1, 0 & 
      , 259653260.4, 0 /
      DATA (E(IH),IH=481,510)  /  384856198.1, 0 & 
      , 165356182, 0 & 
      , 131240969.4, 0 & 
      , 113774621.6, 0 & 
      , 405043602.2, 46020000 & 
      , 66207404.11, 0 & 
      , 92132061.56, 54390000 & 
      , 52331626.84, 0 & 
      , 4176769.025, 0 & 
      , 108949436.1, 259830000 & 
      , -2021441.533, 317150000 & 
      , 0, 355901345.9 & 
      , 0, 217007102.4 & 
      , 0, 214914712 & 
      , 0, 112520403.1 /
      DATA (E(IH),IH=511,540)  /  66580000, 11157426.57 & 
      , 0, 670190892.2 & 
      , 0, 313844820.3 & 
      , 4180000, 263245198.7 & 
      , 0, 465232193.1 & 
      , 0, 178465787 & 
      , 0, 185700190.6 & 
      , 0, 287562460.5 & 
      , 0, 485826868.1 & 
      , 22860000, 93515998.37 & 
      , 10060000, 7570000 & 
      , 163800000, -4660000 & 
      , 49890000, 24770000 & 
      , 13890000, 439116845.8 & 
      , 58280000, 66144642.01 /
      DATA (E(IH),IH=541,570)  /  20340000, 393338860.7 & 
      , 19710000, 387807163.3 & 
      , 8810000, 477264514.5 & 
      , 10790000, 446954897.3 & 
      , 10320000, 81229166.75 & 
      , 29540000, 287763381.1 & 
      , 0, 244964046.4 & 
      , 6280000, 349953693.6 & 
      , 0, 121159499.7 & 
      , 0, 163429348 & 
      , 0, 168331045.4 & 
      , 175430000, 275388910.5 & 
      , 0, 302323591.4 & 
      , 0, 307225288.8 & 
      , 0, 222587604.3 /
      DATA (E(IH),IH=571,600)  /  0, 224722517.8 & 
      , 328220000, 333121697.4 & 
      , 26930000, 31831697.35 & 
      , 82100000, 326398622 & 
      , 68450000, 222514367.3 & 
      , 68810000, 217972669.9 & 
      , 59180000, 191275736.1 & 
      , 24020000, 84112627.12 & 
      , 23810000, 77813209.56 & 
      , 3100000, 126237151.9 & 
      , 93120000, 162124535 & 
      , 66580000, 59833483.64 & 
      , 51210000, 116204324.5 & 
      , -3600000, 124438849.2 & 
      , 38490000, 112396232.4 /
      DATA (E(IH),IH=601,630)  /  97780000, 95935180.99 & 
      , 750000, 90776201.99 & 
      , 7950000, 113261908.7 & 
      , 7950000, 498118106.9 & 
      , 0, 188385703.1 & 
      , 0, 183484005.8 & 
      , 53190000, 49772058.13 & 
      , 8810000, 131181607.4 & 
      , 10790000, 225293668.9 & 
      , 10060000, 7570000 & 
      , -6580000, 42690000 & 
      , 6820000, 317430000 & 
      , 227195745.3, 306690000 & 
      , 199388811.5, 225520000 & 
      , 208453066.2, 13770000 /
      DATA (E(IH),IH=631,660)  /  202562865.8, 0 & 
      , 251837390.6, -550000 & 
      , 197154773.7, 0 & 
      , 222220634.9, -550000 & 
      , 176826295.8, 0 & 
      , 119581532.1, 73600000 & 
      , 33969193.32, 80290000 & 
      , 292530403.7, 104520000 & 
      , 70750000, 325362179.3 & 
      , -1650000, 321281479.9 & 
      , 126050000, 415200216.9 & 
      , 0, 298106582.2 & 
      , 0, 159212338.8 & 
      , 10560000, 58003268.83 & 
      , -4180000, 60093366.86 /
      DATA (E(IH),IH=661,690)  /  0, 283928817.9 & 
      , 0, 400211007.7 & 
      , 0, 0 & 
      , 124722796.7, 0 & 
      , 388340836.9, 0 & 
      , 292840725.8, 0 & 
      , 300995751.7, 0 & 
      , 423845942.6, 0 & 
      , 0, 148357731.6 & 
      , 0, 405407770.7 & 
      , 0, 309907659.6 & 
      , -1610000, 53895145.38 & 
      , -1610000, 70962079.17 & 
      , 50000, 50000 & 
      , 4250000, 4250000 /
      DATA (E(IH),IH=691,720)  /  403226009, -8510000 & 
      , 42093448.02, 25100000 & 
      , 121750000, 81879319.52 & 
      , 98580000, 46780000 & 
      , 80320978.73, 1370000 & 
      , 101195463.3, -4070000 & 
      , 99932232.35, 28270000 & 
      , 97481333.73, 15550000 & 
      , 78671916.16, 2250000 & 
      , 134505858.5, 58200000 & 
      , 60572190.24, 29930000 & 
      , 108053241.6, 40960000 & 
      , 19937079.03, 31920000 & 
      , 4807661.464, 6070000 & 
      , 48091603.77, 48790000 /
      DATA (E(IH),IH=721,750)  /  36678986.93, 51210000 & 
      , 13120145.24, 37490000 & 
      , -6689272.325, -3600000 & 
      , 21354669.98, 38490000 & 
      , 9312053.143, -3160000 & 
      , 443157004, 0 & 
      , 551409473.7, 125970000 & 
      , 7652014.074, 225520000 & 
      , 398918197.4, 234300000 & 
      , 178660000, 205722336 & 
      , 130120000, 236628846.9 & 
      , 7950000, 282089544 & 
      , 10490000, 193072256.9 & 
      , 64610000, 199616108.7 & 
      , 90990000, 35716538.81 /
      DATA (E(IH),IH=751,780)  /  9940000, 282304319.6 & 
      , 221750000, 269326148.2 & 
      , 45190000, 92766148.25 & 
      , 54680000, 526163576.5 & 
      , 50710000, 474617428.3 & 
      , 41880000, 347576512.3 & 
      , 45100000, 303220364 & 
      , 0, 298085379.1 & 
      , 0, 250509230.9 & 
      , 0, 361129903.9 & 
      , 0, 313553755.6 & 
      , -1610000, 68051706.63 & 
      , -1610000, 20475558.39 & 
      , 0, 227513612.4 & 
      , -7530000, 360458358.8 /
      DATA (E(IH),IH=781,810)  /  6280000, 340861711.7 & 
      , 48740000, 10347911.29 & 
      , 48740000, 57924059.54 & 
      , -3600000, 21052436.03 & 
      , -3600000, 68628584.28 & 
      , 38490000, 9009819.194 & 
      , 38490000, 56585967.44 & 
      , 7950000, 468923868 & 
      , 3890000, 106918837.3 & 
      , 3710000, 162360276.5 & 
      , 0, 114055959.4 & 
      , 0, 108982566.6 & 
      , 0, 152178000.8 & 
      , 0, 145380321.9 & 
      , 0, 364611266.9 /
      DATA (E(IH),IH=811,840)  /  0, 354038492.4 & 
      , 34780000, 204229236.9 & 
      , 0, 496077857.7 & 
      , 57130000, 125794777.8 & 
      , 64280000, 86927999.16 & 
      , 0, 319718019.1 & 
      , 40890000, 425256613.8 & 
      , 20920000, 84019579.62 & 
      , 470090000, 42676920.07 & 
      , 517180000, 43750141.44 & 
      , 396220000, 205418218.7 & 
      , 51800000, 63167753.16 & 
      , 46700000, 62969450.52 & 
      , 51210000, 10871629.23 & 
      , 38660000, 44338407.86 /
      DATA (E(IH),IH=841,870)  /  0, 46427788.33 & 
      , 15650000, 15238990.3 & 
      , 14350000, 37056153.97 & 
      , 1800000, 70522932.6 & 
      , 95540000, 64113537.13 & 
      , 82840000, 97430315.77 & 
      , 95540000, 88629398.27 & 
      , 82840000, 121946176.9 & 
      , 4770000, 3890000 & 
      , 281251274.7, 3710000 & 
      , 150862792.8, 29540000 & 
      , 119296129.1, 69040000 & 
      , 219503410.5, 74060000 & 
      , 270540189.1, 247300000 & 
      , 293316778.6, 72900000 /
      DATA (E(IH),IH=871,900)  /  118916778.6, 0 & 
      , 282628077.3, 0 & 
      , 345672602, 0 & 
      , 407656286.7, -2490000 & 
      , 104687514.3, 0 & 
      , 268762043.3, -1790000 & 
      , 553019079.5, -1610000 & 
      , 52594404.81, 50000 & 
      , 4250000, 410070984.8 & 
      , 0, 236611298.7 & 
      , 8370000, 51416466.09 & 
      , 0, 299655823.4 & 
      , 0, 361639508.1 & 
      , 0, 222745264.6 & 
      , -2490000, 58670735.62 /
      DATA (E(IH),IH=901,930)  /  10460000, 355851374.7 & 
      , -1650000, 330125493.9 & 
      , -15110000, 570010743.3 & 
      , 503750000, 505010000 & 
      , 18250000, 9460000 & 
      , 1080000, 209660000 & 
      , 73780000, 163270000 & 
      , -6400000, 116250000 & 
      , 133780000, 132560000 & 
      , 123500000, 156520000 & 
      , 8430000, -4180000 & 
      , -4180000, 77480000 & 
      , 77850000, 87880000 & 
      , 308010000, 278430000 & 
      , 16320000, 80220857.27 /
      DATA (E(IH),IH=931,960)  /  5150000, 12150000 & 
      , -7630000, 188280000 & 
      , 188280000, 134000000 & 
      , 157390000, 6110000 & 
      , 29730000, 309570000 & 
      , 16740000, -4730000 & 
      , -1250000, 188280000 & 
      , 160330000, 162800000 & 
      , 6190000, 112570000 & 
      , 177850000, 297060000 & 
      , 16320000, 9280000 & 
      , 138910000, 62340000 & 
      , 30540000, 0 & 
      , 0, 3560000 & 
      , 3560000, -4390000 /
      DATA (E(IH),IH=961,990)  /  3560000, 0 & 
      , 0, 0 & 
      , 0, 17570000 & 
      , -1420000, 176570000 & 
      , 35310000, 56900000 & 
      , 162290000, 40378380.57 & 
      , 12800000, 164270220.3 & 
      , 0, 0 & 
      , 0, 0 & 
      , 0, -4180000 & 
      , -4180000, 79500000 & 
      , 79500000, -4180000 & 
      , -4180000, 82720000 & 
      , 97910000, 0 & 
      , 0, 0 /
      DATA (E(IH),IH=991,1020)  /  404100000, 410200000 & 
      , 414130000, 10160000 & 
      , 5820000, -1190000 & 
      , 196290000, -8890000 & 
      , 70890000, 70890000 & 
      , 127700000, 189480000 & 
      , 158660000, -4180000 & 
      , -4180000, 64440000 & 
      , 97570000, 62990000 & 
      , 15980000, 165710000 & 
      , 136770000, -65590000 & 
      , 394690000, 384520000 & 
      , -47650000, 92290000 & 
      , -41330000, 159330000 & 
      , 0, 156050000 /
      DATA (E(IH),IH=1021,1050)  /  140030000, 58870000 & 
      , 263240000, 226420000 & 
      , 60540000, -4180000 & 
      , -4180000, 535970000 & 
      , 478230000, 110049848.5 & 
      , 83680000, 10430000 & 
      , 75341336.28, 31940000 & 
      , -1250000, 23740000 & 
      , 58200000, 58200000 & 
      , -5090000, 320000 & 
      , 129040000, 121750000 & 
      , 6900000, 310700000 & 
      , 0, 0 & 
      , 96320000, 88070000 & 
      , 89750000, -4180000 /
      DATA (E(IH),IH=1051,1080)  /  10880000, 5810000 & 
      , -1420000, 49870000 & 
      , 36400000, 71250000 & 
      , 168820000, 7520000 & 
      , 30140000, 21530000 & 
      , 32220000, 32840000 & 
      , 6640000, 85610000 & 
      , 8370000, 150584803.3 & 
      , 38250000, 194865690 & 
      , 85090000, 195688911.4 & 
      , 405860000, 548074803.3 & 
      , 596540000, 262706144.6 & 
      , 130960000, 575688821.6 & 
      , 37240000, 336070493.3 & 
      , 104430000, 357243714.7 /
      DATA (E(IH),IH=1081,1110)  /  99800000, 70960000 & 
      , 535008868.8, 71190000 & 
      , 677453672.1, 131820000 & 
      , 262035013.4, 13230000 & 
      , 191386104.5, 0 & 
      , 497693591.2, 23150000 & 
      , 44794932.53, 18710000 & 
      , 28239401.86, 619590000 & 
      , 132912595.6, 73550000 & 
      , 19964083.4, 3070000 & 
      , 12528608.14, 115340000 & 
      , 97958978.78, 161810000 & 
      , 30384498.09, 44480000 & 
      , 346145985.9, 0 & 
      , 364710510.6, 619590000 /
      DATA (E(IH),IH=1111,1140)  /  135870254.3, 70780000 & 
      , 20151742.05, 6090000 & 
      , 18506266.79, 465490000 & 
      , -848724.5271, 48740000 & 
      , 15492763.27, -3600000 & 
      , 26197288.01, 5370000 & 
      , 386582470.9, 13310000 & 
      , 292005478.6, 18570000 & 
      , 279884457.4, 65930000 & 
      , 311868721.5, 99850000 & 
      , 365946869.1, 6000000 & 
      , 279683497.7, 619550000 & 
      , 121345107.2, 71530000 & 
      , 6416595.001, 18300000 & 
      , 16231119.74, 619550000 /
      DATA (E(IH),IH=1141,1170)  /  122636363.9, 71530000 & 
      , 7707851.712, 18300000 & 
      , 17522376.45, 8080000 & 
      , 181300070.2, 13230000 & 
      , 186404982.5, 0 & 
      , 26418832.42, 0 & 
      , 26373744.7, 65930000 & 
      , 27235427.42, 65930000 & 
      , 28481596.41, 99850000 & 
      , 92841063.38, 99850000 & 
      , 92795975.65, 161810000 & 
      , 60849308.23, 44480000 & 
      , 376610796, 0 & 
      , 395175320.8, 161810000 & 
      , 58089370.29, 44480000 /
      DATA (E(IH),IH=1171,1200)  /  373850858.1, 0 & 
      , 392415382.8, 619550000 & 
      , 111793778.6, 70380000 & 
      , -4284733.639, 6090000 & 
      , -5530208.899, 619550000 & 
      , 110343473.6, 70380000 & 
      , -5735038.634, 6090000 & 
      , -6980513.894, 71240000 & 
      , 85000673.5, 55900000 & 
      , 170621365.3, 619550000 & 
      , 125051622.4, 71530000 & 
      , 10123110.2, 18300000 & 
      , 19937634.94, 13230000 & 
      , 188836244.6, 0 & 
      , 28805006.84, 65930000 /
      DATA (E(IH),IH=1201,1230)  /  33328117.05, 99850000 & 
      , 95227237.8, 619550000 & 
      , 118341439.3, 70380000 & 
      , 2262927.138, 6090000 & 
      , 1017451.878, 161810000 & 
      , 52706961.72, 44480000 & 
      , 368468449.5, 0 & 
      , 387032974.3, 9070000 & 
      , 32761800.62, 470000 & 
      , 500210459.3, 619550000 & 
      , 126936787.3, 68420000 & 
      , 8898275.049, 6090000 & 
      , 9612799.79, 5370000 & 
      , 409641149.6, 5370000 & 
      , 412605223.6, 18570000 /
      DATA (E(IH),IH=1231,1260)  /  323503486.8, 18570000 & 
      , 322257317.8, 13230000 & 
      , 286489670.2, 0 & 
      , 766849162.3, 0 & 
      , 769813236.2, 65930000 & 
      , 325022940.8, 65930000 & 
      , 326536709.8, 99850000 & 
      , 400179905.4, 99850000 & 
      , 403143979.4, 6000000 & 
      , 318366492.8, 6000000 & 
      , 317075236.1, 9070000 & 
      , 242726314.2, 9070000 & 
      , 253355059.9, 470000 & 
      , 720803718.6, 619550000 & 
      , 110170749.6, 68420000 /
      DATA (E(IH),IH=1261,1290)  /  -7867762.628, 6090000 & 
      , -7153237.887, 619550000 & 
      , 110170749.6, 71530000 & 
      , -4757762.628, 18300000 & 
      , 5056762.113, 294550000 & 
      , 18258059.5, 5370000 & 
      , 413991985.3, 18570000 & 
      , 321562695.3, 8080000 & 
      , 210386999.4, 0 & 
      , 488597249.4, 6000000 & 
      , 318811875.7, 619550000 & 
      , 114867691.1, 71530000 & 
      , -60821.12064, 18300000 & 
      , 9753703.62, 8080000 & 
      , 289585820.2, 619550000 /
      DATA (E(IH),IH=1291,1320)  /  85281739.28, 71530000 & 
      , -29646772.93, 18300000 & 
      , -19832248.19, 8080000 & 
      , 230159310.1, 8080000 & 
      , 279772179, 9070000 & 
      , 9070000, 470000 & 
      , 486900000, 137417681 & 
      , 9450000, 93059168.76 & 
      , 51650000, 2998222.535 & 
      , 51650000, 4630000 & 
      , 82149751.19, 0 & 
      , 146653693.5, 155440000 & 
      , 10625496.27, 53970000 & 
      , 70740025.27, 0 & 
      , 92521076.66, 0 /
      DATA (E(IH),IH=1321,1350)  /  117036937.8, 0 & 
      , 123947539.5, 4770000 & 
      , 111472570.2, 3710000 & 
      , 29540000, 29540000 & 
      , 176000000, 487565243.7 & 
      , 147300000, 103650000 & 
      , 0, 240921557.8 & 
      , 73920000, 244021162 & 
      , 0, 18361752.59 & 
      , 0, 125520000 & 
      , 72591594.15, 150620000 & 
      , 0, 0 & 
      , 0, -14770000 & 
      , 0, -8510000 & 
      , 454700000, 5670000 /
      DATA (E(IH),IH=1351,1380)  /  220933331.5, 8370000 & 
      , 406220625.7, 23970000 & 
      , 0, 521300387.2 & 
      , 486900000, 147515954.9 & 
      , 9450000, 103157442.7 & 
      , 46280000, 103588858.9 & 
      , 4630000, 92248025.11 & 
      , 0, 156751967.4 & 
      , 155440000, 20723770.18 & 
      , 53970000, 80838299.18 & 
      , 0, 102619350.6 & 
      , 4770000, 29540000 & 
      , 147300000, 103650000 & 
      , 0, 246624332.8 & 
      , 73920000, 249723937 /
      DATA (E(IH),IH=1381,1410)  /  0, 0 & 
      , 99800000, 377480000 & 
      , 377480000, 454700000 & 
      , 5670000, 17420000 & 
      , 50953235.53, 373240000 & 
      , -1780318.268, 436010000 & 
      , -6505423.169, 285890000 & 
      , 218394895.1, 337550000 & 
      , 46091406.02, 192650000 & 
      , 22297497.04, 14160000 & 
      , 72231169.53, -2520000 & 
      , 118595694.3, 13480000 & 
      , 48745764.6, 19750000 & 
      , 21919143.14, 10650000 & 
      , 62631751.96, 18420000 /
      DATA (E(IH),IH=1411,1440)  /  447591213.4, 16630000 & 
      , 79341310.47, 93120000 & 
      , 160103077.4, 61440000 & 
      , 52672026.04, 50000000 & 
      , 151028340.4, 0 & 
      , 335373742.3, 0 & 
      , 341715698, -9410000 & 
      , 50475531.31, 99800000 & 
      , 47370000, 253335838.4 & 
      , 69380000, 329695787.4 & 
      , 58910000, 492254025.9 & 
      , 20380000, 20127461.91 & 
      , 20920000, 14578044.35 & 
      , 6280000, 69071986.65 & 
      , 41590000, 50249369.82 /
      DATA (E(IH),IH=1441,1470)  /  160000, 106603252.2 & 
      , 93090000, 34502491.8 & 
      , 154180000, 40040180.59 & 
      , 120540000, 112843432.7 & 
      , 0, 374503979.6 & 
      , 0, 368414562 & 
      , 0, 437548504.3 & 
      , -14770000, 131310307.1 & 
      , 342000000, 10060000 & 
      , 7570000, -6580000 & 
      , 157130000, 42690000 & 
      , 6820000, 522120000 & 
      , 155660097.1, 51870000 & 
      , 118501584.9, 30760000 & 
      , 91302167.33, -5490000 /
      DATA (E(IH),IH=1471,1500)  /  124186109.6, 152400000 & 
      , 619590000, 143869255.6 & 
      , 68420000, 25790743.41 & 
      , 61500000, 12781325.85 & 
      , 2590000, 23005268.15 & 
      , 57060000, 23342651.31 & 
      , 0, 538432054.9 & 
      , 0, 111429984.6 & 
      , 0, 262943843.8 & 
      , 30080000, 70689588.9 & 
      , -6900000, -6900000 & 
      , 0, 430050869.8 & 
      , 0, 329900840.3 & 
      , 0, 434547961.2 & 
      , 17420000, 58920697.52 /
      DATA (E(IH),IH=1501,1530)  /  13480000, 56713226.59 & 
      , 25990000, 29030618.02 & 
      , 10650000, 7601200.453 & 
      , 1530000, 67615142.76 & 
      , 58200000, -5598525.469 & 
      , 29930000, 41882525.92 & 
      , 161940000, 6105276.279 & 
      , 134330000, 9669610.664 & 
      , 0, 308431098.5 & 
      , 0, 371475623.2 & 
      , -3220000, 314123006.4 & 
      , 14300000, 94307425.97 & 
      , -1650000, 387229057.3 & 
      , -1650000, 375286263.6 & 
      , -4180000, 21620000 /
      DATA (E(IH),IH=1531,1560)  /  -2420000, 170770000 & 
      , 136430000, 113030000 & 
      , 67350000, 0 & 
      , 180100000, 326900000 & 
      , 2750000, 146684040.1 & 
      , 38490000, 14154671.17 & 
      , 29540000, 104134582.7 & 
      , 29540000, 116077376.4 & 
      , 13480000, -6432403.456 & 
      , 61500000, 4782324.487 & 
      , -1790000, 542709192.8 & 
      , -1610000, 71632313.4 & 
      , 50000, 4250000 & 
      , 393486532.8, 0 & 
      , 84130331.83, 396590000 /
      DATA (E(IH),IH=1561,1590)  /  21966408.23, 460930000 & 
      , 18728819.02, 14160000 & 
      , 72627896.03, 10650000 & 
      , 63028478.46, -2520000 & 
      , 118992420.8, 192650000 & 
      , 22694223.54, 61440000 & 
      , 53068752.54, 93120000 & 
      , 160499803.9, 17420000 & 
      , 50939563.4, 13480000 & 
      , 49168706.55, 12930000 & 
      , 337550000, 285890000 & 
      , 218312410.8, 6660000 & 
      , 283813818.9, 67750000 & 
      , 288901438.6, 34110000 & 
      , 361376776.5, 0 /
      DATA (E(IH),IH=1591,1620)  /  31130000, 237463423.1 & 
      , 20930000, 282691027.5 & 
      , -2750000, 58340000 & 
      , 24700000, 99800000 & 
      , 396590000, 19435800.39 & 
      , 460930000, 342000000 & 
      , 14160000, 70097288.19 & 
      , 10650000, 60497870.63 & 
      , -2520000, 116461812.9 & 
      , 192650000, 20163615.7 & 
      , 61440000, 50538144.7 & 
      , 93120000, 157969196.1 & 
      , 10060000, 7570000 & 
      , -6580000, 157130000 & 
      , 42690000, 6820000 /
      DATA (E(IH),IH=1621,1650)  /  17420000, 67532957.61 & 
      , 17420000, 50175252.2 & 
      , 13480000, 65762100.76 & 
      , 13480000, 6660000 & 
      , 280885485.3, 67750000 & 
      , 34110000, 363143073.1 & 
      , 0, 31980000 & 
      , 235385089.6, -9410000 & 
      , 342000000, 10060000 & 
      , 7570000, -6580000 & 
      , 157130000, 42690000 & 
      , 6820000, 17420000 & 
      , 72227587.75, 13480000 & 
      , 17420000, 50409540.47 & 
      , 13480000, 49731754.41 /
      DATA (E(IH),IH=1651,1680)  /  396590000, 9851879.559 & 
      , 478030000, 12814647.68 & 
      , 285890000, 207412768.1 & 
      , 337550000, 103227017.7 & 
      , 14160000, 60513367.35 & 
      , 10650000, 50913949.79 & 
      , -2520000, 106877892.1 & 
      , 192650000, 10579694.87 & 
      , 93120000, 148385275.3 & 
      , 61440000, 40954223.87 & 
      , 18960000, 18960000 & 
      , 0, 318345379.3 & 
      , 0, -9410000 & 
      , 33447168.22, 99800000 & 
      , 47370000, 261259801.2 /
      DATA (E(IH),IH=1681,1710)  /  70750000, 347359225.5 & 
      , 93090000, 59454817.65 & 
      , 154180000, 46086416.7 & 
      , 0, 399456305.4 & 
      , 0, 393366887.9 & 
      , 0, 462500830.2 & 
      , -14770000, 156262633 & 
      , 342000000, 10060000 & 
      , 7570000, -6580000 & 
      , 157130000, 42690000 & 
      , 6820000, 619590000 & 
      , 143541341.3, 376980000 & 
      , 50909087.01, 284180000 & 
      , 537844017.3, 68420000 & 
      , 25462829.1, 3070000 /
      DATA (E(IH),IH=1711,1740)  /  23157353.84, 23390000 & 
      , 25122529.07, 269000000 & 
      , -2380843.389, 12930000 & 
      , 73785671.75, 16630000 & 
      , 77485671.75, 18420000 & 
      , 447154599.3, 63980000 & 
      , 368610363.6, 61500000 & 
      , 12453411.54, 30080000 & 
      , 69161864.49, 37580000 & 
      , 249923385.2, 0 & 
      , 536904330.4, 0 & 
      , 109902260.2, 0 & 
      , 261416119.4, 22000000 & 
      , 56045263, 310110000 & 
      , 186005764.3, 522120000 /
      DATA (E(IH),IH=1741,1770)  /  154241072.4, 30260000 & 
      , 95472560.24, -4780000 & 
      , 123477085, 20200000 & 
      , 94324468.15, 173220000 & 
      , 10008887.75, 152400000 & 
      , 46692372.83, -1650000 & 
      , 171611520.7, 73920000 & 
      , 176361125, 246860000 & 
      , 208812409.9, 16180000 & 
      , 193395741.4, 20920000 & 
      , 18950000, 84925002.67 & 
      , 18950000, 449214284.2 & 
      , 23390000, 26652213.94 & 
      , 30080000, 96437429.5 & 
      , 30080000, 95146172.79 /
      DATA (E(IH),IH=1771,1800)  /  37580000, 37580000 & 
      , 0, 564179895.5 & 
      , 0, 562888638.7 & 
      , 0, 137177825.2 & 
      , 0, 135886568.5 & 
      , 310110000, 220339310.1 & 
      , 522120000, 157830718.5 & 
      , 30260000, 99062206.29 & 
      , -4780000, 127066731 & 
      , 20200000, 97914114.2 & 
      , 152400000, 87534546.43 & 
      , 0, 73920000 & 
      , 30080000, 30080000 & 
      , 30080000, 30080000 & 
      , 30080000, 23390000 /
      DATA (E(IH),IH=1801,1806)  /  23390000, 23390000 & 
      , 23390000, 23390000 & 
      , 23390000, 23390000 & 
     /
