*IDENT MOSES_CTILE_ESSNOW_SEG
*/ 
*/ MOSES 2.2 + Coastal Tiling for FAMOUS  
*/           + multi-level essery snow  
*/           + D1 space for subgrid hypsometry (25 levels, hardcoded) !!!NOW UNNEC!!?
*/           + JWilliams land tuning param setup (namelist_famous_land_cc.mod)
*/           + Bug fixes. Changing from 25 to 10 elevation classes
*/ 
*/ Based on the MOSES 2.2 mods for vn4.5:
*/   abx0f406, abx1f406, abx2f406, abx3f406, abx4f406, amv1f406,
*/   apa1f406, are1f406, are2f406, are3f406, newdecks, scatter_fix,
*/   update_m22.mdk 
*/
*/ With coastal tiling for MOSES 2.2 derived from:
*/   ang0f503.mf77 (Nic Gedney, Apr 01)
*/  
*/ Also includes:
*/   adb1f406: radiation mod (J.M. Edwards, May 99)
*/   Tibet_snow.mf77
*/   runoff_MII.mf77
*/  
*/ Should be used with: 
*/   ccouple_moses2.2: FAMOUS coupling mod
*/   LANDTOGLOBAL, PTOUV_LAND, PTOUV_SEA: coastal tiling mods
*/   amv2f406: *optional* MOSES 2.2 mod for radiative canopy
*/   arerf406_ctile: reconfiguration mod
*/ 
*/ Decks modified: ADDRES1, ATMPHY1, ATMSTEP1, BLEND_H, BL_CTL1, 
*/   BL_IC7A, BOUYTQ5B, BOUYTQ6A, CLDCTL1, CNTLATM, COMPT2A, CRADINCS, 
*/   C_ROUGH, DARCY5A, DECAY2A, DESCENT, DFPLN3A, DIAG3A, EXFXUV5B, 
*/   FCDCH6A, FCDCH7A, FILL3A, FTSA1A, FXCA3A,G ROWT2A, HYDCON5A, 
*/   HYDROL7A, HYDR_CT1, HYD_IC7A, IMCLTQ7A, INITIAL1, INITVEG1, 
*/   LEAF7A, LOTKA2A, LWRAD3A, NAMSIZE, NVEGPARM, PHIMH6A, PHYSIO7A, 
*/   PPARM1A, PRELIM1, PSLIMS1, RAD_CTL1, READFL1A, ROOTFR7A, RPANCA1A, 
*/   SCREEN7A, SEED, SETMODL1, SFEVAP7A, SFEXCH7A, SFFLUX7A, SFLINT6A, 
*/   SFLINT7A, SFMELT7A, SFREST7A, SFRIB7A, SFSNOW7A, SICEHT5B, 
*/   SMCEXT7A, SOILHT7A, SOILHY5A, SOILMC7A, SPARM1A, STATMPT1, 
*/   STDEV7A, SURFHY7A, SWRAD3A, TRIF, TRIFD2A, TYPPTRA, TYPSIZE, 
*/   U_MODEL1, UPANCIL1, VEG1A, VEG2A, VEG_CTL1, VEG_IC1A, VEG_IC2A
*/ 
*/ Decks added: ALBPFT, ALBSNOW, BDY_EXPL7A, BDY_EXPL8A, BDY_IMPL8A,
*/   CANCAP7A, GAUSS7A, IMBLPT18A, IMBLPT28A, IMSFPT8A, RESTILE7A,
*/   SFEXPL8A, SFIMPL8A, SICEHT7A, SOILEV7A, TILEALB
*/
*/ Annette Osprey 25 January 2008
*/
*/-----
*COMDECK C_LAND_CC
*/-----
      REAL F0(NPFT_1),LAI_MIN(NPFT_1),NL0(NPFT_1),R_GROW(NPFT_1),
     +TLOW(NPFT_1),TUPP(NPFT_1),Q10,V_CRIT_ALPHA,KAPS
C
      COMMON /COMMON_LAND_CC/ F0,LAI_MIN,NL0,R_GROW,
     +TLOW,TUPP,Q10,V_CRIT_ALPHA,KAPS
*/-----
*COMDECK C_EDDY
*/-----
*/Some variables in the EDDY namelist in CNTLOCN are being
*/changed also so read them in too.
      REAL FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH
C
      COMMON /COMMON_EDDY/ FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH
*/-----
*DECLARE ADDRES1
*/-----
*D ABX2F404.79,ABX2F404.80   
! All surface tiles
          ILOUT = NTILES
*/-----
*DECLARE ARGPTRA
*/-----
*I GDR7F405.105
!
     &   JSNOWDEPTH, JRHO_SNOW_GRND, JSNOW_CAN, JSNOW_SOIL_HTF,
     &   JNSNOW, JDS, JSICE, JSLIQ, JTSNOWLAYER, JRHO_SNOW, JRGRAINL,
     &   J_DEEP_ICE_TEMP,
*/-----
*DECLARE ARTPTRA
*/-----
*I GDR7F405.90
!
     &A_SPPTR(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),
     &A_SPPTR(A_IXPTR(61)),A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),
     &A_SPPTR(A_IXPTR(64)),A_SPPTR(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),
     &A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),A_SPPTR(A_IXPTR(69)),
*/-----
*DECLARE ATMPHY1
*/-----
*I APB1F401.7
     &             DOWNWELLING_LW,
*B ABX1F404.278   
      REAL SURF_HTF_TILE(LAND_FIELDDA,NTILESDA)
      REAL EI_TILE(TILE_FIELDDA,NTILESDA)
      REAL DOWNWELLING_LW(P_FIELDDA,BL_LEVELSDA)
*D ARE1F404.16
     &     TILE_FIELDDA,NTILESDA,                                       
*I AJS1F401.177   
     &        NTILESDA,             ! dynamic allocation: NTILES        
*D ARE2F404.3,ARE2F404.8    
*D ARE1F404.20,ARE1F404.26   
     &     ECAN_TILE(TILE_FIELDDA,NTILESDA)! Canopy evaporation from    
!                                          ! land tiles                 
     &    ,MELT_TILE(TILE_FIELDDA,NTILESDA)! Snowmelt on land tiles     
     &    ,LW_DOWN(P_FIELDDA)              ! Downward LW radiation      
     &    ,DOLR(P_FIELDDA)                 ! TOA - surface upward LW rad
     &    ,SW_TILE(TILE_FIELDDA,NTILESDA)  ! SW radiation on land tiles 
     &    ,TILE_FRAC(TILE_FIELDDA,NTILESDA)! Land-surface type fractions
*D ARE1F404.37
*D ARE1F404.39,ARE1F404.52   
*I ATMPHY1.142   
!  Similar for land surface albedos - different values are needed for   
!    direct & diffuse sunlight if the HadCM2 approximate treatment of   
!    sulphate aerosol is being used:                                    
      IF ( L_H2_SULPH ) THEN                                            
         NLALBS = 2                                                     
       ELSE                                                             
         NLALBS = 1                                                     
      ENDIF                                                             
                                                 
*D ARE2F404.10
     &     ,NTILES,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE
     &     ,D1(JLAND_ALB+1),D1(JSICE_ALB+1)                    
*D ARN1F404.92
     &     ,WORKF,BL_LEVELS,L_RADHEAT,RADHEAT_DIM1,NLALBS
*D AWI1F403.120,AWI1F403.127  
*D ARE2F404.12,ARE2F404.17   
      SAL_DIM = 1                                                       
*IF DEF,A03_7A                                                          
! Extra workspace required if MOSES II is used                          
      SAL_DIM = P_FIELD                                                 
*ENDIF                                                                  
*D ARE2F404.18
     &  NTILES,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,
     &  D1(JLAND_ALB+1),D1(JSICE_ALB+1),                       
*I APBBF401.1
     &             DOWNWELLING_LW,
*D ARE1F404.53,ARE1F404.55   
     &            NTILES,TILE_FIELDDA,TILE_PTS,TILE_INDEX,              
     &            DOLR,LW_DOWN,SW_TILE,ECAN_TILE,MELT_TILE,TILE_FRAC,   
*D ARE1F404.56,ARE1F404.57   
     &              NTILES,TILE_FIELDDA,TILE_PTS,TILE_INDEX,            
     &              ECAN_TILE,MELT_TILE,TILE_FRAC,                      
*D AJS1F401.194,AJS1F401.195
      CALL BL_CTL(WORK1,WORK2,EI_TILE,WORK12,WORK3,WORK13,WORK4,
     &            WORK14,WORK6,SURF_HTF_TILE,DOWNWELLING_LW,
     &            WORK9,WORK10,WORK5,
*D AJS1F401.204,AJS1F401.205   
      CALL HYDR_CTL(WORK2,EI_TILE,WORK12,WORK3,WORK13,
     &              WORK4,WORK14,SURF_HTF_TILE,WORK5,WORK6,
*/-----
*DECLARE ATMSTEP1
*/-----
*I GSS2F305.90
     &                    DOWNWELLING_LW,
*B ATMSTEP1.49
      REAL  DOWNWELLING_LW (P_FIELD,BL_LEVELS)
*I APB1F401.1
     &              DOWNWELLING_LW,
*I ACB2F405.41    
      ELSEIF ( (H_SECT(3).EQ."08A") )THEN
        L_RADHEAT = .TRUE.
        RADHEAT_DIM1 = P_FIELDDA
        TILE_FIELD = LAND_FIELD
*D ARE1F404.13
     &              TILE_FIELD,NTILES,
*/-----
*DECLARE BLEND_H
*/-----
*D BLEND_H.5
      PARAMETER (LB = 20.0)                                            
*/-----
*DECLARE BL_CTL1
*/-----
*D AJS1F401.217,AJS1F401.218
      SUBROUTINE BL_CTL(CLOUD_FRACTION,SNOW_SUBLIMATION,EI_TILE,
     &           SNOWMELT,CANOPY_EVAPORATION,EXT,
*D AJS1F401.219
     &           SOIL_EVAPORATION,SURF_HT_FLUX_LAND,SURF_RADFLUX,
     &           SURF_HTF_TILE,DOWNWELLING_LW,
*D ARE1F404.59,ARE1F404.61   
     &           NTILESDA,TILE_FIELDDA,TILE_PTS,TILE_INDEX,             
     &           OLR,LW_DOWN,SW_TILE,ECAN_TILE,MELT_TILE,TILE_FRAC,     
*D AJS1F401.228
     &       SURF_HT_FLUX_LAND(P_FIELDDA),         !                    
     &       SURF_HTF_TILE(LAND_FIELDDA,NTILESDA),
     &       DOWNWELLING_LW(P_FIELDDA,BL_LEVELSDA),
*I AJS1F401.229   
     &       SURF_HT_FLUX(P_FIELDDA),              
     &       SURF_HT_FLUX_SICE(P_FIELDDA),           
*I ARE1F404.65    
     &       NTILESDA,                             ! IN                 
*D ARE1F404.69,ARE1F404.74   
     &       OLR(P_FIELDDA),                       ! IN                 
     &       LW_DOWN(P_FIELDDA),                   ! IN                 
     &       SW_TILE(TILE_FIELDDA,NTILESDA),       ! IN                 
     &       ECAN_TILE(TILE_FIELDDA,NTILESDA),     ! OUT                
     &       MELT_TILE(TILE_FIELDDA,NTILESDA),     ! OUT                
     &       TILE_FRAC(TILE_FIELDDA,NTILESDA)      ! OUT                
*D ARE1F404.77,ARE1F404.78   
     &       EI_TILE(TILE_FIELDDA,NTILESDA),                            
     &       ESOIL_TILE(TILE_FIELDDA,NTILESDA),                         
     &       FTL_TILE(TILE_FIELDDA,NTILESDA),                           
     &       GS_TILE(TILE_FIELDDA,NTILESDA),                            
     &       LE_TILE(TILE_FIELDDA,NTILESDA),                            
*I ARE1F404.80    
     &       Q1P5M_TILE(TILE_FIELDDA,NTILESDA),                         
     &       T1P5M_TILE(TILE_FIELDDA,NTILESDA),                         
     &       RAD_TILE(TILE_FIELDDA,NTILESDA),                           
*D ARE1F404.83
     &       RIB_TILE(TILE_FIELDDA,NTILESDA)                            
*D ARE1F404.85
     &   RHO_ARESIST_TILE(TILE_FIELDDA,NTILESDA),                       
*D ARE1F404.87
     &   ARESIST_TILE(TILE_FIELDDA,NTILESDA),                           
*D ARE1F404.89
     &   RESIST_B_TILE(TILE_FIELDDA,NTILESDA),                          
*D ABX1F405.200
*D ABX1F405.205
     & ,     PLLTILE(NTILESDA)  ! pseudolevel list for surface tiles
     & ,     L_CTILE            ! Coastal tiling switch for  
!                               ! BL_INTCT    
*I APC5F400.3     
     &    .OR.SF(328,3)     ! Needed for T at 1.5 m over land tiles     
*I APC5F400.4     
     &    .OR.SF(329,3)     ! Needed for Q at 1.5 m over land tiles     
*D AJS1F401.256
        SURF_HT_FLUX_LAND(I) = 0.0                                      
*I ACN1F405.19    

! Set coastal tiling flag
      IF ( H_SECT(3).EQ.'07A' .OR. H_SECT(3).EQ.'08A') THEN
        L_CTILE = .TRUE.
      ENDIF                                                 
 
*D AYY1F404.66
     & L_BL_LSPICE,L_MOM,L_MIXLEN,L_CTILE,                              
*D ARN0F405.26
*I ANG1F405.2     
! OUT Coastal tiling diagnostics :
                                   
     & STASHWORK(SI(339,3,im_index)),
     & STASHWORK(SI(391,3,im_index)),STASHWORK(SI(392,3,im_index)),
     & STASHWORK(SI(393,3,im_index)),STASHWORK(SI(394,3,im_index)),
     & STASHWORK(SI(389,3,im_index)),STASHWORK(SI(390,3,im_index)),
     & STASHWORK(SI(347,3,im_index)),STASHWORK(SI(353,3,im_index)),
     & STASHWORK(SI(343,3,im_index)),STASHWORK(SI(381,3,im_index)),
     & STASHWORK(SI(395,3,im_index)),
                                                                       
*I AJS1F401.297
     & D1(JOROG),
*I BL_CTL1.179
     & DOWNWELLING_LW,
*I AJS1F401.307   
     & SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,                      
     & SURF_HTF_TILE,
*D ABX1F405.218,ARE1F404.95   
     & A_INTHD(23),                                                     
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,NTILES,                           
     & D1(JCANHT_PFT),D1(JCAN_WATER_TYP),D1(JCATCH_TYP),                
*D ARE1F404.97
     & LW_DOWN,SW_TILE,D1(JZ0_TYP),                                     
*I ACN1F405.20    
     & D1(JFRAC_LAND),                           
*D ARE1F404.99
     & OLR,D1(JSNODEP_TYP),D1(JTSTAR_TYP),                              
*I ARE1F404.101   
     & D1(JTSTAR_LAND),D1(JTSTAR_SEA),D1(JTSTAR_SICE),                  
*D ARE1F404.103
     & ECAN_TILE,EI_TILE,ESOIL_TILE,FTL_TILE,                           
     & GS_TILE,LE_TILE,MELT_TILE,RAD_TILE,                              
*D ARE1F404.106,ABX1F405.220  
     & RIB_TILE,Q1P5M_TILE,T1P5M_TILE,TILE_INDEX,TILE_PTS,TILE_FRAC,    
     & L_ESSERY_SNOW,
     & D1(JTSNOWLAYER(1)), D1(JNSNOW(1)),D1(JSNOWDEPTH(1)),
     & D1(JSICE(1)), D1(JSLIQ(1)),D1(JDS(1)),  NSMAX,
     & D1(J_DEEP_ICE_TEMP(1)),
*I AJS1F401.403   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.404   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.TRUE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)   
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO2,RESS_SO2,RES_FACTOR)                   
*ENDIF                                                                  
!                                                                       
*I AWO3F405.66    
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.67    
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.TRUE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)   
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_NH3,RESS_NH3,RES_FACTOR)                   
*ENDIF                                                                  
!                                                                       
*I AJS1F401.469   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.470   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_AIT,RESS_SO4_AIT,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AJS1F401.519   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.520   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_ACC,RESS_SO4_ACC,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AJS1F401.569   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.570   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_DIS,RESS_SO4_DIS,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AWO3F405.175   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.176   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_FreshSoot,RESS_Soot,RES_FACTOR)            
*ENDIF                                                                  
!                                                                       
*I AWO3F405.235   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.236   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_AgedSoot,RESS_Soot,RES_FACTOR)             
*ENDIF                                                                  
!                                                                       
*I AWO3F405.282   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.283   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SootInCloud,RESS_Soot,RES_FACTOR)          
*ENDIF                                                                  
!                                                                       
*I APBGF401.39    
          D1(JTSTAR_LAND+TOP_ROW_START+I-2)=0.0                         
          D1(JTSTAR_SEA+TOP_ROW_START+I-2)=0.0                          
          D1(JTSTAR_SICE+TOP_ROW_START+I-2)=0.0                         
*I APBGF401.47    
          D1(JTSTAR_LAND+P_BOT_ROW_START+I-2)=0.0                       
          D1(JTSTAR_SEA+P_BOT_ROW_START+I-2)=0.0                        
          D1(JTSTAR_SICE+P_BOT_ROW_START+I-2)=0.0                       
*I APB2F401.96    
C     
       CALL POLAR(D1(JTSTAR_LAND),D1(JTSTAR_LAND),D1(JTSTAR_LAND),      
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
C     
       CALL POLAR(D1(JTSTAR_SEA),D1(JTSTAR_SEA),D1(JTSTAR_SEA),      
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
C     
      CALL POLAR(D1(JTSTAR_SICE),D1(JTSTAR_SICE),D1(JTSTAR_SICE),       
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
*I APBGF401.57    
          D1(JTSTAR_LAND+TOP_ROW_START+I-2)=                            
     &      D1(JTSTAR_LAND+TOP_ROW_START+ROW_LENGTH+I-2)                
          D1(JTSTAR_SEA+TOP_ROW_START+I-2)=                             
     &      D1(JTSTAR_SEA+TOP_ROW_START+ROW_LENGTH+I-2)                 
          D1(JTSTAR_SICE+TOP_ROW_START+I-2)=                            
     &      D1(JTSTAR_SICE+TOP_ROW_START+ROW_LENGTH+I-2)                
*I APBGF401.66    
          D1(JTSTAR_LAND+P_BOT_ROW_START+I-2)=                          
     &      D1(JTSTAR_LAND+P_BOT_ROW_START-ROW_LENGTH+I-2)              
          D1(JTSTAR_SEA+P_BOT_ROW_START+I-2)=                           
     &      D1(JTSTAR_SEA+P_BOT_ROW_START-ROW_LENGTH+I-2)               
          D1(JTSTAR_SICE+P_BOT_ROW_START+I-2)=                          
     &      D1(JTSTAR_SICE+P_BOT_ROW_START-ROW_LENGTH+I-2)              
*I ABX1F405.444
      DO J=1,NTILES
        DO I=1,TILE_FIELDDA
          IF (tile_frac(I,J).lt.1e-3) SW_TILE(I,J)=0.
        END DO
      END DO
CL ITEM 335 SURFACE NET SW RAD ON TILES                             

      IF (SF(335,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,335,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(335,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),SW_TILE(1,PSLEVEL),                         
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF


CL ITEM 417 _TILE_ FRAC MAP - NOT TYPE like 317!

      IF (SF(417,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,417,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(417,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),TILE_FRAC(1,PSLEVEL),                         
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
        ENDIF

*I BL_CTL1.298   
! Coastal tiling diagnostics

C Item 337:
C Land heat flux from surface to level 1 (land mean) (W/m2)

      IF (SF(337,3)) THEN                                               
                                                                        
        CALL COPYDIAG(STASHWORK(SI(337,3,im_index)),SURF_HT_FLUX_LAND,  
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,337,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
                                                                        
      END IF                                                            

C Item 338:
C Net surface sea-ice heat flux (sea mean) (W/m2) 

      IF (SF(338,3)) THEN                                               
                                                                        
        CALL COPYDIAG(STASHWORK(SI(338,3,im_index)),SURF_HT_FLUX_SICE,  
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,338,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
                                                                        
      END IF                                                            

! End of coastal tiling diagsnostics
                                                                        
*D ARE1F405.18,AJS1F401.704  
*D ABX1F405.221
CL ITEM 287: CANOPY EVAPORATION ON TILES                                
*D ABX1F405.224
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.226
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.232,ABX1F405.233  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.244
CL ITEM 288: TRANSPIRATION + SOIL EVAPORATION ON TILES                  
*D ABX1F405.247
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.249
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.255,ABX1F405.256  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.293
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.295
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.301,ABX1F405.302  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.368
CL ITEM 294: BULK RICHARDSON NUMBER ON TILES                            
*D ABX1F405.371
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.373
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.379,ABX1F405.380  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*I ABX1F405.389   
CL ITEM 314: SURFACE NET RADIATION ON TILES                             
*D ABX1F405.391,ABX1F405.405  
      IF (SF(314,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,314,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.411,ABX1F405.412  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.415,ABX1F405.416  
     &          STASHWORK(SI(314,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),RAD_TILE(1,PSLEVEL),                         
*D ABX1F405.422
*D ABX1F405.426
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.428
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.434,ABX1F405.435  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.524
CL ITEM 321: CANOPY WATER CONTENT ON TILES                              
*D ABX1F405.527
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.529
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.535,ABX1F405.536  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.540
     &           *P_FIELD),D1(JCAN_WATER_TYP+((PSLEVEL-1)*LAND_FIELD)), 
*D ABX1F405.547
CL ITEM 322: CANOPY CAPACITY ON TILES                                   
*D ABX1F405.550
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.552
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.558,ABX1F405.559  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.563
     &           *P_FIELD),D1(JCATCH_TYP+((PSLEVEL-1)*LAND_FIELD)),     
*D ABX1F405.570,ABX1F405.578  
*D ABX1F405.582
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.584
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.590,ABX1F405.591  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*I ABX1F405.623   
                                                                        
CL ITEM 328: TEMPERATURE AT 1.5M OVER TILES                             
                                                                        
      IF (SF(328,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,328,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(328,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),T1P5M_TILE(1,PSLEVEL),                       
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 329: SPECIFIC HUMIDITY AT 1.5M OVER TILES                       
                                                                        
      IF (SF(329,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,329,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(329,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),Q1P5M_TILE(1,PSLEVEL),                       
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 330: SURFACE LATENT HEAT FLUX ON TILES                          
                                                                        
      IF (SF(330,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,330,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(330,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),LE_TILE(1,PSLEVEL),                          
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 331: SUBLIMATION ON TILES                                       
                                                                        
      IF (SF(331,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,331,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(331,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),EI_TILE(1,PSLEVEL),                          
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 332: TOA outward LW radiation after boundary layer              
        IF (SF(332,3)) THEN                                             
        CALL COPYDIAG(STASHWORK(SI(332,3,im_index)),OLR,                
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,332,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
        END IF                                                          
*I TJ181193.16
CL ITEM 402 SURF_HTF_TILE
      IF (SF(402,3)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,402,3,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(402,3,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),surf_htf_tile(1,PSLEVEL),     
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF
*/-----
*DECLARE BL_IC7A
*/-----
*D BL_IC7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I BL_IC7A.65
     & GBM_OROG,
*I BL_IC7A.74
     & DOWNWELLING_LW,
*D BL_IC7A.76,BL_IC7A.77   
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RAD_SICE,                       
     & TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN,L_CTILE,              
*I BL_IC7A.86    
! OUT Coastal tiling diagnostics :              
     & RIB_SSI,TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,VSHR_LAND,VSHR_SSI,
     & E_SSI,EI_SICE,FTL_SSI,RADNET_SICE,FLANDG,                        
                                                                        
*D BL_IC7A.102
     & SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,
     & SURF_HTF_TILE,
     & ZH,T1_SD,Q1_SD,ERROR,                               
*D ABX1F405.754
     & ASTEPS_SINCE_TRIFFID,                                            
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,NTILES,                           
*D BL_IC7A.109
     & FRAC,LW_DOWN,SW_TILE,Z0V_TILE,                                   
*I ACN1F405.114   
     & FLAND,                          
*D BL_IC7A.111
     & OLR,SNOW_TILE,TSTAR_TILE,                                        
*I BL_IC7A.112   
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,                  
*D BL_IC7A.114,ABX1F405.755  
     & ECAN_TILE,EI_TILE,ESOIL_TILE,FTL_TILE,GS_TILE,LE_TILE,MELT_TILE, 
     & RADNET_TILE,G_LEAF,GPP_FT,NPP_FT,RESP_P_FT,RESP_S,RESP_W_FT,     
*D BL_IC7A.117,ABX1F405.756  
     & RIB_TILE,Q1P5M_TILE,T1P5M_TILE,TILE_INDEX,TILE_PTS,TILE_FRAC,
     & L_ESSERY_SNOW,
     & TSNOWLAYER, NSNOW,SNOW_DEPTH, SNOW_SICE_1, SNOW_SLIQ_1,
     & SNOW_DEPTH_1, NSMAX,T_DEEP_ICE,
*I BL_IC7A.186   
     &,ASTEPS_SINCE_TRIFFID      ! IN Number of atmospheric             
!                                !    timesteps since last call         
!                                !    to TRIFFID.                       
*I ACN1F405.116   
     &,NTILES                    ! IN No. of land-surface tiles         
*D BL_IC7A.197,BL_IC7A.198  
     &,CANOPY_TILE(LAND_FIELD,NTILES)                                   
!                                ! IN Surface/canopy water for          
*D BL_IC7A.200
     &,CATCH_TILE(LAND_FIELD,NTILES)                                    
*D BL_IC7A.202
!                                !    land tiles (kg per sq m)          
*I BL_IC7A.203   
     &,FRAC(LAND_FIELD,NTYPE)    ! IN Fractions of surface types.       
     &,FLAND(LAND_FIELD)         ! IN Land fraction on land tiles.    
     &,FLANDG(P_FIELD)           ! IN Land fraction on all tiles.     
*I BL_IC7A.210   
     &,LW_DOWN(P_FIELD)          ! IN Surface downward LW radiation     
!                                !    (W/m2).                           
*D BL_IC7A.227,BL_IC7A.229  
     &,SW_TILE(LAND_FIELD,NTILES)! IN Surface net SW radiation on land  
!                                !    tiles (W/m2).                     
*D BL_IC7A.233
     &,Z0V_TILE(LAND_FIELD,NTILES)!IN Tile roughness lengths (m).       
*D BL_IC7A.266,BL_IC7A.271  
     &,RAD_SICE(P_FIELD)         ! IN Surface net SW and downward LW    
!                                !    radiation on sea-ice fraction     
!                                !    (W/sq m, positive downwards).     
*I ABX1F405.767   
     &,L_ESSERY_SNOW
     &,L_CTILE                  ! IN true if model to be run with 
!                               !    coastal tiling       
*I BL_IC7A.308   
     &,OLR(P_FIELD)              ! IN    TOA - surface upward LW on     
!                                !       last radiation timestep        
!                                ! OUT   Corrected TOA outward LW       
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                                ! INOUT Snow on tiles (kg/m2).         
*I BL_IC7A.312   
     &,TSTAR_LAND(P_FIELD)       ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(P_FIELD)        ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(P_FIELD)       ! INOUT Sea-ice sfc temperature (K). 
*D BL_IC7A.314
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D BL_IC7A.357,BL_IC7A.358  
     &,EI_TILE(LAND_FIELD,NTILES)!OUT EI for land tiles                 
     &,ESOIL_TILE(LAND_FIELD,NTILES)                                    
                                ! OUT ES for land tiles                 
*D BL_IC7A.367
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*I BL_IC7A.368   
     &,E_SSI(P_FIELD)           ! OUT   Surface FQW for mean sea.       
     &,EI_SICE(P_FIELD)         ! OUT   Sea-ice sumblimation            
!                               !       (sea mean).                     
     &,FTL_SSI(P_FIELD)         ! OUT sea mean surface heat flux        
*I BL_IC7A.371   
     &,GS_TILE(LAND_FIELD,NTILES)!OUT Surface conductance for           
!                               !     land tiles.                       
*I BL_IC7A.373   
     &,LE_TILE(LAND_FIELD,NTILES)!OUT Surface latent heat flux for      
!                               !     land tiles (W/m2).                
     &,MELT_TILE(LAND_FIELD,NTILES)                                     
!                               ! OUT Snowmelt on tiles (kg/m2/s).      
*I BL_IC7A.377   
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                               ! OUT Tile surface net radiation.       
*D BL_IC7A.389
     &,RIB_TILE(LAND_FIELD,NTILES)                                      
!                               ! OUT RIB for land tiles.               
     &,SURF_HTF_TILE(LAND_FIELD,NTILES)
     &,RIB_LAND(P_FIELD)        !     Land mean bulk Richardson no.  
!                                        for lowest layer.              
     &,RIB_SSI(P_FIELD)         ! OUT Sea mean bulk Richardson no.   
!                                        for lowest layer.              
*I BL_IC7A.395   
     &,SURF_HT_FLUX_LAND(P_FIELD)                               
!                               ! OUT Net downward heat flux at         
!                               !     surface over land                 
!                               !     fraction of gridbox (W/m2).       
     &,SURF_HT_FLUX_SICE(P_FIELD)                               
!                               ! OUT Net downward heat flux at         
!                               !     surface over sea-ice              
!                               !     fraction of gridbox (W/m2).       
*I BL_IC7A.400   
     &,TAUX_LAND(U_FIELD)       ! OUT W'ly compt of land sfc wind    
!                               !     stress (N/sq m). (On UV-grid    
!                               !     with first and last rows       
!                               !     undefined or, at present,      
!                               !     set to missing data            
     &,TAUX_SSI(U_FIELD)        ! OUT W'ly compt of mean sea sfc wind
!                               !     stress (N/sq m). (On UV-grid    
!                               !     with first and last rows       
!                               !     undefined or, at present,      
!                               !     set to missing data            
*D ABX1F405.784,ABX1F405.787  
     &,TAUY_LAND(U_FIELD)       ! OUT S'ly compt of land sfc wind    
!                               !     stress (N/sq m).  On UV-grid;   
!                               !     comments as per TAUX.          
     &,TAUY_SSI(U_FIELD)        ! OUT S'ly compt of mean sea sfc wind
!                               !     stress (N/sq m).  On UV-grid;   
!                               !     comments as per TAUX.          
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
!                               ! OUT Tile fractions. Equal to surface  
!                               !     type fractions if NTILES=NTYPE,   
!                               !     or 1 if NTILES=1.                 
*I BL_IC7A.405   
     &,VSHR_LAND(P_FIELD)       ! OUT Magnitude of land sfc-to-lowest
!                                     atm level wind shear (m per s).   
     &,VSHR_SSI(P_FIELD)        ! OUT Mag. of mean sea sfc-to-lowest 
!                                     atm level wind shear (m per s).   
*D BL_IC7A.413
     &,RHO_ARESIST_TILE(LAND_FIELD,NTILES)                              
*D BL_IC7A.415
     &,ARESIST_TILE(LAND_FIELD,NTILES)                                  
*D BL_IC7A.417
     &,RESIST_B_TILE(LAND_FIELD,NTILES)                                 
*I BL_IC7A.418   
     &,RADNET_SICE(P_FIELD)     ! OUT Sea-ice surface net radiation.  
*I BL_IC7A.436   
     &,Q1P5M_TILE(LAND_FIELD,NTILES)                                    
!                               ! OUT Q1P5M over land tiles.            
*I BL_IC7A.437   
     &,T1P5M_TILE(LAND_FIELD,NTILES)                                    
!                               ! OUT T1P5M over land tiles.            
*I ARN0F405.180
      INTEGER
     & NSMAX
!
!
      REAL
     & TSNOWLAYER(LAND_FIELD, NTILES, NSMAX)
     & ,SNOW_DEPTH(LAND_FIELD, NTILES)
     & ,NSNOW(LAND_FIELD, NTILES)
     & ,SNOW_SICE_1(LAND_FIELD, NTILES, NSMAX)
     & ,SNOW_SLIQ_1(LAND_FIELD, NTILES, NSMAX)
     & ,SNOW_DEPTH_1(LAND_FIELD, NTILES, NSMAX)
     & ,T_DEEP_ICE(LAND_FIELD, ST_LEVELS)
!
*D BL_IC7A.444,BL_IC7A.445  
     & ECAN_TILE(LAND_FIELD,NTILES)                                     
                      ! OUT ECAN for land tiles                         
*D BL_IC7A.455,BL_IC7A.460  
*D BL_IC7A.480,BL_IC7A.482  
     & CANHC_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Areal heat capacity of canopy   
!                               !       for land tiles (J/K/m2).        
     &,VFRAC_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Fractional canopy coverage for  
!                               !       land tiles.                     
     &,WT_EXT_TILE(LAND_FIELD,SM_LEVELS,NTILES)                         
*D BL_IC7A.485
!                               !       soil layer by each tile.        
*D BL_IC7A.491,BL_IC7A.494  
                                                                        
      REAL                                                              
     & ALPHA1(LAND_FIELD,NTILES)! LOCAL Mean gradient of saturated      
!                               !       specific humidity with respect  
!                               !       to temperature between the      
!                               !       bottom model layer and tile     
!                               !       surfaces                        
     &,ALPHA1_SICE(P_FIELD)     ! LOCAL ALPHA1 for sea-ice.             
     &,ASHTF(P_FIELD)           ! LOCAL Coefficient to calculate        
!                               !       surface heat flux into soil or  
!                               !       sea-ice.                        
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Coefficient to calculate        
!                               !       surface heat flux into land     
!                               !       tiles.                          
     &,DTRDZ(P_FIELD,BL_LEVELS) ! LOCAL -g.dt/dp for model layers.      
     &,FLAKE(LAND_FIELD,NTILES) ! LOCAL Lake fraction.                  
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Surface FQW for land tiles      
     &,FQW_ICE(P_FIELD)         ! LOCAL Surface FQW for sea-ice         
     &,FTL_ICE(P_FIELD)         ! LOCAL Surface FTL for sea-ice         
     &,QW(P_FIELD,BL_LEVELS)    ! LOCAL Total water content             
     &,TL(P_FIELD,BL_LEVELS)    ! LOCAL Ice/liquid water temperature    
     &,TSTAR_TILE_OLD(LAND_FIELD,NTILES)                                
!                               ! LOCAL Tile surface temperatures at    
!                               !       beginning of timestep.          
     &,FRACA(LAND_FIELD,NTILES) ! LOCAL Fraction of surface moisture    
!                               !       flux with only aerodynamic      
!                               !       resistance for snow-free land   
!                               !       tiles.                          
     &,RESFS(LAND_FIELD,NTILES) ! LOCAL Combined soil, stomatal         
!                               !       and aerodynamic resistance      
!                               !       factor for fraction (1-FRACA)   
!                               !       of snow-free land tiles.        
     &,RESFT(LAND_FIELD,NTILES) ! LOCAL Total resistance factor.        
!                               !       FRACA+(1-FRACA)*RESFS for       
!                               !       snow-free land, 1 for snow.     
     &,RHOKH_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Surface exchange coefficients   
!                               !       for land tiles                  
     &,RHOKH_SICE(P_FIELD)      ! LOCAL Surface exchange coefficients   
!                               !       for sea and sea-ice             
     &,RHOKPM(LAND_FIELD,NTILES)! LOCAL Land surface exchange coeff.    
     &,RHOKPM_SICE(P_FIELD)     ! LOCAL Sea-ice surface exchange coeff. 
     &,H_BLEND_OROG(P_FIELD)    ! LOCAL Blending height used as part of 
!                               !       effective roughness scheme      
     &,Z0H(P_FIELD)             ! LOCAL Roughness length for heat and   
!                               !       moisture (m).                   
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Tile roughness lengths for heat 
!                               !       and moisture (m).               
     &,Z0M(P_FIELD)             ! LOCAL Roughness length for            
!                               !       momentum (m).                   
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Tile roughness lengths for      
!                               !       momentum.                       
     &,Z0M_EFF(P_FIELD)         ! LOCAL Effective grid-box roughness    
!                               !       length for momentum             
     &,CDR10M_UV(U_FIELD)       ! LOCAL Ratio of CD's reqd for          
!                               !       calculation of 10 m wind. On    
!                               !       UV-grid; comments as per RHOKM. 
     &,CHR1P5M(LAND_FIELD,NTILES)!LOCAL Ratio of coefffs for            
!                               !       calculation of 1.5m temp for    
!                               !       land tiles.                     
     &,CHR1P5M_SICE(P_FIELD)    ! LOCAL CHR1P5M for sea and sea-ice     
!                               !       (leads ignored).                
     &,CT_CTQ(P_FIELD,BL_LEVELS)! LOCAL Coefficient in T and q          
!                                       tri-diagonal implicit matrix    
     &,CQ_CM(U_FIELD,BL_LEVELS) ! LOCAL Coefficient in U and V          
!                                       tri-diagonal implicit matrix    
     &,DQW(P_FIELD,BL_LEVELS)   ! LOCAL BL increment to q field         
     &,DTL(P_FIELD,BL_LEVELS)   ! LOCAL BL increment to T field         
     &,DU(U_FIELD,BL_LEVELS)    ! LOCAL BL increment to u wind field    
     &,DV(U_FIELD,BL_LEVELS)    ! LOCAL BL increment to v wind field    
     &,RDZ(P_FIELD,BL_LEVELS)   ! LOCAL RDZ(,1) is the reciprocal of    
!                               !       the height of level 1, i.e. of  
!                               !       the middle of layer 1.  For     
!                               !       K > 1, RDZ(,K) is the           
!                               !       reciprocal of the vertical      
!                               !       distance from level K-1 to      
!                               !       level K.                        
     &,RDZUV(U_FIELD,BL_LEVELS) ! LOCAL RDZ (K > 1) on UV-grid.         
!                               !       Comments as per RHOKM (RDZUV).  
     &,FB_SURF(P_FIELD)         ! Surface flux buoyancy over density    
!                               ! (m^2/s^3)                             
!                                                                       
     &,U_S(P_FIELD)             ! Surface friction velocity (m/s)       
     &,TV1_SD(P_FIELD)          ! Standard deviation of turbulent       
!                               ! fluctuations of surface layer         
!                               ! virtual temperature (K).              
     &,E_LAND(P_FIELD)          ! LOCAL FQW over mean land      
     &,EI_LAND(P_FIELD)         ! LOCAL EI over mean land       
     &,FTL_LAND(P_FIELD)        ! LOCAL FTL over mean land      
     &,FLANDG_UV(U_FIELD)       ! Land frac (on UV-grid, with 1st 
!                               ! and last rows undefined or, at 
!                               ! present, set to "missing data")
     &,TSTAR_SSI(P_FIELD)       ! LOCAL Sea mean sfc temperature (K).
                                                                        
*I BL_IC7A.497   
     &,J                        ! LOCAL Tile point index                
*I BL_IC7A.527
      REAL
      !FOR CALL TO Z
     & ZLB(P_FIELD,0:BL_LEVELS)
     &,TV(P_FIELD,BL_LEVELS)
     &,DZL(P_FIELD,BL_LEVELS)
      
      INTEGER K

      REAL
      !FOR DERIVATION OF ELEVATION ADJUSTMENTS
     & GBM_OROG(P_FIELD),ELEVATION(NELEV),delevation,
     & BLDEPTH,DUMMY,
!
     & TGRAD,QGRAD,CFGRAD,QCLGRAD,QCFGRAD,SOILGRAD,ICEGRAD,
     & PGRAD,EXNERGRAD,LWGRAD,
!+seg 3/2019     
     & dLWdT,
!-seg
!
     & T_ELEV(LAND_FIELD,NELEV),Q_ELEV(LAND_FIELD,NELEV),
     & CF_ELEV(LAND_FIELD,NELEV),QCF_ELEV(LAND_FIELD,NELEV),
     & QCL_ELEV(LAND_FIELD,NELEV),
     & tl_elev(LAND_FIELD,NELEV),qw_elev(LAND_FIELD,NELEV),
     & z1_elev(LAND_FIELD,NELEV),
     & tsoil_elev(LAND_FIELD,NELEV),tice_elev(LAND_FIELD,NELEV),
     & p_elev(LAND_FIELD,NELEV),exner_elev(LAND_FIELD,NELEV),
     & downwelling_lw(P_FIELD,BL_LEVELS),lwdown_elev(LAND_FIELD,NELEV)

!+seg These are mid-points. Changed to 10 levels (CESM type) from 25
      data elevation / 100., 300., 550., 850., 1150.,
     &                 1450.,1800.,2250.,2750.,3600./
      data soilgrad /-.001/
      data icegrad /-.001/
*D BL_IC7A.598
*I BL_IC7A.632   
! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:                             
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN                                                   
                                                                        
        DO I=P1,P1+P_POINTS-1
          FLANDG(I)=0.0                                          
          IF(ICE_FRACT(I).LE.0.0)THEN                               
            TSTAR_SSI(I)=TSTAR_SEA(I)                             
          ELSE                                                        
            TSTAR_SSI(I)=ICE_FRACT(I)*TSTAR_SICE(I)             
     &        +(1.0-ICE_FRACT(I))*TSTAR_SEA(I)                    
          ENDIF                                                       
        ENDDO           
        DO L=1,LAND1+LAND_PTS-1                                         
          I = LAND_INDEX(L)                                             
          FLANDG(I)=FLAND(L)                                            
        ENDDO                                                           
                                                                        
      ENDIF                                                             
                                                                        
*D BL_IC7A.634
! Call TILEPTS to calculate TILE_PTS and TILE_INDEX for surface types   
*D BL_IC7A.650
     & N_ROWS,FIRST_ROW,ROW_LENGTH,FLANDG,                              
*D BL_IC7A.653
     & VSHR,VSHR_LAND,VSHR_SSI,Z1                                       
*D BL_IC7A.662
     & P_FIELD,SM_LEVELS,NTILES,TILE_PTS,TILE_INDEX,                    
*I BL_IC7A.665   
     & CANHC_TILE,VFRAC_TILE,FLAKE,                                     
*D ABX1F405.802
     & RESP_P,RESP_P_FT,RESP_S,RESP_W_FT,SMC,WT_EXT_TILE                
*I BL_IC7A.668   
                                                                        
!---------------------------------------------------------------------- 
! If TRIFFID is being used apply any correction to the land-atmosphere  
! fluxes on the first timestep after the last TRIFFID call. Such a      
! correction will typically be associated with a total depletion of     
! carbon or with maintanence of the seed fraction. The corrections      
! are stored in the accumulation variables after the call to TRIFFID.   
! The correction is added to the instantaneous land-atmosphere fluxes   
! (so that the atmospheric carbon budget is corrected) but is not       
! included in the accumulation variables which drive TRIFFID, since     
! this has already been dealt with during the last TRIFFID call.        
!---------------------------------------------------------------------- 
      IF (L_TRIFFID.AND.(ASTEPS_SINCE_TRIFFID.EQ.1)) THEN               
        DO N=1,NPFT                                                     
          DO L=LAND1,LAND1+LAND_PTS-1                                   
            NPP_FT(L,N)=NPP_FT(L,N)+NPP_FT_ACC(L,N)/TIMESTEP            
            RESP_P_FT(L,N)=RESP_P_FT(L,N)-NPP_FT_ACC(L,N)/TIMESTEP      
            NPP_FT_ACC(L,N)=-NPP_FT_ACC(L,N)                            
          ENDDO                                                         
        ENDDO                                                           
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          RESP_S(L)=RESP_S(L)+RESP_S_ACC(L)/TIMESTEP                    
          RESP_S_ACC(L)=-RESP_S_ACC(L)                                  
        ENDDO                                                           
      ENDIF                                                             
*D BL_IC7A.685
*D BL_IC7A.687
! Reset TILE_PTS and TILE_INDEX and set tile fractions to 1 if aggregate
! tiles are used (NTILES=1).                                            
! Otherwise, set tile fractions to surface type fractions.              
*D BL_IC7A.689
      DO N=1,NTILES                                                     
*D BL_IC7A.691
          TILE_FRAC(L,N) = 0.                                           
*D BL_IC7A.695,BL_IC7A.703  
      IF (NTILES.EQ.1) THEN                                             
        TILE_PTS(1) = LAND_PTS                                          
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          TILE_FRAC(L,1) = 1.                                           
          TILE_INDEX(L+1-LAND1,1) = L                                   
        ENDDO                                                           
      ELSEIF (NTILES.EQ.2*NELEV) THEN
        DO L=LAND1,LAND1+LAND_PTS-1
          DO N=1,NTYPE-NELEV
            k=mod(n-1,nelev)+1
            TILE_FRAC(L,K) = TILE_FRAC(L,K)+FRAC(L,N)
          END DO
          DO N=NTYPE-NELEV+1,NTYPE
            k=mod(n-1,nelev)+1
            TILE_FRAC(L,NELEV+K) = TILE_FRAC(L,NELEV+K)+FRAC(L,N)
          END DO
        END DO
        DO N=1,NTILES
          J=0
          DO L=LAND1,LAND1+LAND_PTS-1
            IF (TILE_FRAC(L,N).gt.0) THEN
              J=J+1
              TILE_INDEX(J,N)=L
            END IF
          END DO
          TILE_PTS(N)=J
        END DO

      ELSE                                                              
        DO N=1,NTYPE                                                    
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            TILE_FRAC(L,N) = FRAC(L,N)                                  
          ENDDO                                                         
        ENDDO                                                           
      ENDIF                                                             

! NEW, TO DO interp of T, Q and rad fluxes to right heights!
      DO I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
      ENDDO
      DO K=1,BL_LEVELS
        CALL Z(P_POINTS,EXNER(P1,K),EXNER(P1,K+1),PSTAR(P1),
     &    AKH(K),BKH(K),Q(P1,K),QCF(P1,K),
     &    QCL(P1,K),T(P1,K),ZLB(P1,K-1),TV(P1,K),
     &    ZLB(P1,K),DZL(P1,K),RDZ(P1,K),LTIMER)
      ENDDO

!C      do l=1,land_field
!+seg 3/2019
      tsoil_elev=0.0
      tice_elev=0.0
!-seg
      do l=land1,land1+land_pts-1
        i=land_index(l)
        bldepth=0.
        DO K=1,BL_LEVELS
          bldepth=bldepth+dzl(i,k)
        END DO
!+seg 3/2019
! Calculated version:
!       tgrad=(t(i,BL_LEVELS)-t(i,1))/bldepth
!
! Use fixed temperature lapse rate
! What to use?
!     DALR ~ 0.0098 deg/m  (9.8 deg/km)
! Model:
!     LAPSE=0.0065     !  NEAR SURFACE LAPSE RATE
!     LAPSE_TROP=0.003 !  TROPOPAUSE LAPSE RATE
! Set to latter
!        tgrad=0.003
! Try with 6 deg/km
! Negative gradient
        tgrad=-0.006
!-seg
        qgrad=(q(i,BL_LEVELS)-q(i,1))/bldepth

        !? What else needs to be done?
        cfgrad=(cf(i,BL_LEVELS)-cf(i,1))/bldepth
        qclgrad=(qcl(i,BL_LEVELS)-qcl(i,1))/bldepth
        qcfgrad=(qcf(i,BL_LEVELS)-qcf(i,1))/bldepth

        pgrad=(pstar(I)-(PSTAR(I)*BKH(BL_LEVELS+1) + AKH(BL_LEVELS+1)))
     &                                                     /bldepth
        exnergrad=(exner(i,BL_LEVELS)-exner(i,1))/bldepth
!+seg 3/2019
! Calculate value
!       lwgrad=(downwelling_lw(i,BL_LEVELS)-downwelling_lw(i,1))/bldepth
!
!       dLW/dz = dLW/dT * dT/dz
!       dT/dz  = tgrad
!       dLW/dt = dLWdT
!       
!       lwgrad = dLWdT * tgrad
!
!        dLWdT  = 5.4
!       From later runs of xmvlv.
!       Note: the relationship is sinusoidal. This is average value.
!       17/5/17
        dLWdT  = 3.6
        lwgrad = dLWdT * tgrad
!-seg


        DO K=1,NELEV
          delevation=elevation(k)-gbm_orog(i)
!          delevation=0.

          call test_delevation(t(i,1),delevation,tgrad,
     &                         t_elev(l,k))

          call test_delevation(q(i,1),delevation,qgrad,
     &                         q_elev(l,k))
          call test_delevation(T_DEEP_SOIL(l,1),delevation,soilgrad,
     &                         tsoil_elev(l,k))
          call test_delevation(T_DEEP_ICE(l,1),delevation,icegrad,
     &                         tice_elev(l,k))
          call test_delevation(lw_down(i),delevation,lwgrad,
     &                         lwdown_elev(l,k))

          delevation=0.

          call test_delevation(qcl(i,1),delevation,qclgrad,
     &                         qcl_elev(l,k))
          call test_delevation(qcf(i,1),delevation,qcfgrad,
     &                         qcf_elev(l,k))
          call test_delevation(cf(i,1),delevation,cfgrad,
     &                         cf_elev(l,k))

          call test_delevation(pstar(I),delevation,pgrad,
     &                         p_elev(l,k))
          call test_delevation(exner(i,1),delevation,exnergrad,
     &                         exner_elev(l,k))

        ENDDO

      enddo


*D BL_IC7A.706,BL_IC7A.713  
! Call boundary layer routine                                           
*D BL_IC7A.716
      CALL SF_EXPL (                                                    
*D BL_IC7A.719,BL_IC7A.720  
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
*D BL_IC7A.724
     & AK(1),BK(1),AKH(1),BKH(1),DELTA_AK(1),DELTA_BK(1),               
*D BL_IC7A.732,BL_IC7A.735  
     & NTILES,NELEV,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY_TILE,CATCH_TILE,FLAKE,GS_TILE,HCON,            
     & HO2R2_OROG,FLAND,FLANDG,
     & SNOW_TILE,SIL_OROG_LAND,SMVCST,STHF,STHU,             
     & TILE_FRAC,VFRAC_TILE,Z0V_TILE,                                   
*D BL_IC7A.738
     & ICE_FRACT,U_0,V_0,                                               
*D BL_IC7A.741
     & CF(1,1),QCF(1,1),QCL(1,1),                                       
*D BL_IC7A.744,ABX1F405.835  
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,
     & VSHR,VSHR_LAND,VSHR_SSI,ZH,                 
     & Q(1,1),T(1,1),T_DEEP_SOIL,TI,                                    
     & TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,
     & TSTAR_TILE,U(1,1),V(1,1),                                  
     & L_BL_LSPICE,                                                     

! IN SUBGRID ELEVATION CORRECTIONS
     & t_elev,q_elev,cf_elev,qcl_elev,qcf_elev,
     & p_elev,tsoil_elev,tice_elev,exner_elev,lwdown_elev,
! OUT SUBGRID ELEVATION CORRECTIONS
     & tl_elev,qw_elev,z1_elev,

*D BL_IC7A.748
     & SFME,SQ1P5,ST1P5,SU10,SV10,                                      
*D BL_IC7A.751,BL_IC7A.752  
     & Z0MSEA,                                                          
*D BL_IC7A.755,BL_IC7A.757  
     & CD,CH,E_SEA,QW(1,1),TL(1,1),FQW(1,1),                            
     & FTL(1,1),FTL_TILE,LE_TILE,H_SEA,RADNET_SICE,RADNET_TILE,         
     & RHOKM(1,1),RIB,RIB_TILE,TAUX(1,1),TAUY(1,1),                     
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
*D BL_IC7A.760,BL_IC7A.761  
     & FME,                                                             
*D BL_IC7A.768,BL_IC7A.769  
! OUT data required for 4D-VAR :                                        
     & RHO_CD_MODV1,                                                    
*D BL_IC7A.772,BL_IC7A.774  
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,                                  
                                                                        
! OUT data required elsewhere in boundary layer or surface code         
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ(1,1),FQW_TILE,         
     & FQW_ICE,FTL_ICE,TSTAR_TILE_OLD,FRACA,RESFS,RESFT,                
     & RHOKH(1,1),RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,             
     & Z1,H_BLEND_OROG,Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z0M_EFF,               
     & CDR10M_UV,CHR1P5M,CHR1P5M_SICE,                                  
     & TSNOWLAYER(:,:,1), NSNOW, SNOW_DEPTH, 
     & SNOW_SICE_1(:,:,1),SNOW_SLIQ_1(:,:,1), SNOW_DEPTH_1(:,:,1),
     & L_ESSERY_SNOW,T_DEEP_ICE,
     & FLANDG_UV,                                               
*I BL_IC7A.778   
                                                                        
      CALL BDY_EXPL (                                                   
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,ROW_LENGTH,                                      
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,U_POINTS,U1,                       
                                                                        
! IN values defining vertical grid of model atmosphere :                
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,              
     & EXNER,                                                           
                                                                        
! IN sea/sea-ice data :                                                 
     & U_0,V_0,                                                         
                                                                        
! IN cloud data :                                                       
     & CF,QCF,QCL,CCA,CCB,CCT,                                          
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,RAD_HR,RADHR_DIM1,                                         
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,                                  
     & H_BLEND_OROG,Z0M_EFF,                                            
     & TIMESTEP,L_BL_LSPICE,L_MOM,                                      
                                                                        
! INOUT data :                                                          
     & Q,T,U,V,ZH,                                                      
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & QW,TL,FQW,FTL,                                                   
     & RHOKH,RHOKM,                                                     
     & TAUX,TAUY,ZHT,                                                   
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,     
                                                                        
! OUT data required for tracer mixing :                                 
     & NRML,                                                            
                                                                        
! OUT data required for 4D-VAR :                                        
     & RHO_KM,                                                          
                                                                        
! OUT data required elsewhere in UM system :                            
     & DTRDZ,RDZ,RDZUV,                                                 
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,                                      
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
                                                                        
      CALL SF_IMPL (                                                    
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
     & P_POINTS,P1,LAND1,LAND_PTS,U_POINTS,U1,                          
                                                                        
! IN soil/vegetation/land surface data :                                
     & LAND_INDEX,LAND_MASK,                                            
     & NTILES,NELEV,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY_TILE,FLAKE,SMC,                                
     & TILE_FRAC,WT_EXT_TILE,                                           
     & FLAND,FLANDG,                                                    
                                                                        
! IN sea/sea-ice data :                                                 
     & DI,ICE_FRACT,U_0,V_0,                                            
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,                         
     & T_DEEP_SOIL,QW(1,1),TL(1,1),U(1,1),V(1,1),RHOKM(1,1),            
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,                             
     & DTRDZ(1,1),DU(1,1),DV(1,1),FQW_TILE,FQW_ICE,FTL_ICE,             
     & TSTAR_TILE_OLD,                                                  
     & FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,      
     & Z1,Z0H,Z0H_TILE,Z0M,Z0M_TILE,CDR10M_UV,CHR1P5M,CHR1P5M_SICE,     
     & CT_CTQ(1,1),DQW(1,1),DTL(1,1),CQ_CM(1,1),                        
     & L_NEG_TSTAR,                                                     
     & FLANDG_UV,                                               

! IN SUBGRID ELEVATION CORRECTIONS
     & tl_elev,qw_elev,z1_elev,lwdown_elev,
     & p_elev,tsoil_elev,tice_elev,
                                                                        
! IN STASH flags :-                                                     
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            
                                                                        
! INOUT data :                                                          
     & TI,TSTAR,
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       
     & TSTAR_TILE,SNOW_TILE,                                   
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 
     & E_SEA,FQW(1,1),FTL(1,1),FTL_TILE,H_SEA,OLR,TAUX(1,1),TAUY(1,1),  
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & ECAN,EI_TILE,ESOIL_TILE,                                         
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    
     & SURF_HTF_TILE,
                                                                        
! OUT diagnostic requiring STASH flags :                                
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     
                                                                        
! OUT data required elsewhere in UM system :                            
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,MELT_TILE,                          
     & ERROR,                                                           
! LOGICAL FOR NEW SNOW SCHEME
     & L_ESSERY_SNOW,T_DEEP_ICE,
!
     & TSNOWLAYER(:,:,1), NSNOW,
!  
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
                                                                        
      CALL BDY_IMPL (                                                   
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,P1,U_FIELD,U1,P_POINTS,U_POINTS,ROW_LENGTH,              
                                                                        
! IN values defining vertical grid of model atmosphere :                
     & BL_LEVELS,                                                       
                                                                        
! IN data :                                                             
     & RHOKH,RHOKM,RDZ,RDZUV,                                           
                                                                        
! INOUT data :                                                          
     & Q,T,U,V,QW,TL,FQW,FTL,TAUX,TAUY,                                 
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,                                      
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                

                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        RIB_LAND(I)=0.0                                             
        RIB_SSI(I)=0.0                                              
        FTL_LAND(I)=0.0                                             
        FTL_SSI(I)=0.0                                              
        E_LAND(I)=0.0                                               
        E_SSI(I)=0.0                                                
        EI_LAND(I)=0.0                                              
        EI_SICE(I)=0.0                                              
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                          
          RIB_LAND(I)=RIB_LAND(I) +                                 
     &      RIB_TILE(L,N)*TILE_FRAC(L,N)
          FTL_LAND(I)=FTL_LAND(I) +                                 
     &      FTL_TILE(L,N)*TILE_FRAC(L,N)                                
          E_LAND(I)=E_LAND(I) +                                     
     &      FQW_TILE(L,N)*TILE_FRAC(L,N)                                
          EI_LAND(I)=EI_LAND(I) +                                   
     &      EI_TILE(L,N)*TILE_FRAC(L,N)                                 
        ENDDO
      ENDDO                                                             
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        IF(FLANDG(I).LT.1.0)THEN                                        
            RIB_SSI(I)=(RIB(I)-RIB_LAND(I)*FLANDG(I))        
     &        /(1.0-FLANDG(I))                                          
            FTL_SSI(I)=(FTL(I,1)-FTL_LAND(I)*FLANDG(I))         
     &        /(1.0-FLANDG(I))                                        
            E_SSI(I)=(FQW(I,1)-E_LAND(I)*FLANDG(I))             
     &        /(1.0-FLANDG(I))                                        
            EI_SICE(I)=(EI(I)-EI_LAND(I)*FLANDG(I))             
     &        /(1.0-FLANDG(I))
        ENDIF                                        
      ENDDO
*I BL_IC7A.783

      SUBROUTINE test_delevation(old,dheight,grad,new)

      real old,new
      real dheight,grad

      new=old
      if (abs(old).gt.1e-12) then
        new=old+(dheight*grad)
        if (sign(1.,new).NE.sign(1.,old)) new=old
      endif

      return

      END
                                                                        
*/-----
*DECLARE BOUYTQ5B
*/-----
*D BOUYTQ5B.45,ADM3F404.359  
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
*I ADM3F404.366   
                                                                        
      REAL
     & BT_CLD(P_FIELD,BL_LEVELS) 
!                             ! DUMMY Used in 6A boundary layer scheme
     &,BQ_CLD(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,A_QS(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,A_DQSDT(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,DQSDT(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
*/-----
*DECLARE BOUYTQ6A
*/-----
*D BOUYTQ6A.2
*IF DEF,A03_8A                                                          
*D BOUYTQ6A.44,BOUYTQ6A.46   
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,L_BL_LSPICE,LTIMER
*I BOUYTQ6A.92    
                                                                        
      REAL
     & CF(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme
     &,TL(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme
     &,BF(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme

      LOGICAL
     & L_BL_LSPICE            ! DUMMY Used in 7A boundary layer scheme
*/-----
*DECLARE CLDCTL1
*/-----
*D ARE2F404.22
     &           NTILESDA,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,            
     &           LAND_ALB,SICE_ALB,
*D ARN1F404.102
     &           L_RADHEAT,RADHEAT_DIM1,NLALBS,P_FIELDDA,               
*I @DYALLOC.747   
     &       NTILESDA,     !  COPY OF NTILES                            
     &       TILE_FIELDDA, ! Set to LAND_FIELD for MOSES II, 1 otherwise
*I ARN1F404.105   
     &       NLALBS,     !  Number of fields of land surface albedo
*D ARE2F404.23,ARE2F404.28   
     &       LW_DOWN(P_FIELDDA),            ! Surface downward LW       
C                                           ! (W/sq m)                  
     &       DOLR(P_FIELDDA),               ! TOA - surf upward LW rad  
     &       SW_TILE(TILE_FIELDDA,NTILESDA),! Surface net SW radiation  
C                                           ! on land tiles (W/sq m)    
     &       LAND_ALB(P_FIELDDA,NLALBS),    ! IN: Mean land albedo      
     &       SICE_ALB(P_FIELDDA,NLALBS),    ! IN: Mean sea-ice albedo   
*I @DYALLOC.755   

C
C MOSES 2.2 coastal tiling variables
      REAL
     &       FRACSOLID(P_FIELDDA),     ! Solid surface fraction in gridb
     &       ICE_FRACT(P_FIELDDA),     ! Ice fraction
     &       FLANDG(P_FIELDDA),        ! Land fraction 
     &       SW_NET_LAND(P_FIELDDA),   ! SW net local flux over land
     &       SW_NET_SICE(P_FIELDDA)    ! SW net local flux over sea-ice 

      REAL
     &       ALBSOLID                  ! Mean solid surface albedo

      LOGICAL 
     &       L_CTILE                   ! Coastal tiling switch
                                       ! Set to .TRUE.
                    
*D ARE2F404.29
     &     ,RADINCS((P_FIELDDA*(P_LEVELSDA+2+NTILES)+511)/512*512*2)  
*D ARE2F404.33
     &       N,J,                                                       
*D ARE2F404.37,CLDCTL1.112  
                                                                        
*D ARE2F404.38
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512*2               
C       no. words for LW/SW                                             
*D ARE2F404.41
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C       offset to 2nd RADINCS                                           
*D ARE2F404.43
C zenith angle adjustment and surface net SW on tiles)                  
*I CLDCTL1.387   


        IF ( H_SECT(3) .EQ. '07A' .OR.                                  
     &       H_SECT(3) .EQ. '08A' ) THEN                                
          L_CTILE=.TRUE. 
        ENDIF          

C  Set up GLOBAL fractional land field:                                 
        CALL LAND_TO_GLOBAL                                             
     &    (D1(JLAND),D1(JFRAC_LAND),FLANDG,LAND_PTS,P_FIELDDA)          

C Set up ice fraction field and solid fraction field                    
        DO I=1,P_FIELD
          ICE_FRACT(I) = D1(JICE_FRACTION+I-1)   
          FRACSOLID(I) = FLANDG(I) + (1.0-FLANDG(I))*ICE_FRACT(I)  
        ENDDO     

*I CLDCTL1.394   
        IF (L_CTILE) THEN
          
          DO I=1,P_FIELD
            SW_NET_LAND(I) = 0.0
            SW_NET_SICE(I) = 0.0
            SURF_RADFLUX(I) = 0.0       
          ENDDO

          DO I=FIRST_POINT,LAST_POINT                                   
                                                                        
! land_alb is only set on points where there is incoming sw radiation   
! at the previous timestep, therefore it will be zero over some         
! land points                                                           
            IF (FRACSOLID(I).GT.0.0) THEN
          
              IF (FLANDG(I).GT.0.0.AND.LAND_ALB(I,1).LE.0.0) THEN
                SW_NET_LAND(I) = 0.0
                SW_NET_SICE(I) = 0.0                      
              ELSE                                         
                ALBSOLID = ( FLANDG(I) * LAND_ALB(I,1) +
     &            (1.0-FLANDG(I)) * ICE_FRACT(I) * SICE_ALB(I,1) )  
     &              /FRACSOLID(I)                                
                                                                        
                IF (FLANDG(I).GT.0.0) THEN
                  SW_NET_LAND(I) = RADINCS(I)                
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)   
     &              * (1.0-LAND_ALB(I,1))/(1.0-ALBSOLID)   
                ENDIF
            
                IF (ICE_FRACT(I).GT.0.0) THEN
                  SW_NET_SICE(I) = RADINCS(I)      
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)  
     &              * (1.0-SICE_ALB(I,1))/(1.0-ALBSOLID)  
                
                SURF_RADFLUX(I) = SW_NET_SICE(I) * ICE_FRACT(I)   
                ENDIF                                                   
              ENDIF 
                                                                      
            ENDIF 
          ENDDO                                          
        ELSE ! If not L_CTILE          
          DO I=FIRST_POINT,LAST_POINT                                   
            SURF_RADFLUX(I) =                                           
     &          RADINCS(I) * COS_ZENITH_ANGLE(I) + RADINCS(I+LEN) 
          ENDDO
        ENDIF

*D CLDCTL1.396,CLDCTL1.397  
*I CLDCTL1.398   

C                                                                       
C Output SW land and sea-ice:                                   

        IF(SF(257,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(257,1,im_index)+I-1) = SW_NET_LAND(I)          
          ENDDO                                                         
          CALL EXTDIAG (STASHWORK,si(1,1,im_index),SF(1,1),257,257,INT9,
     &       ROW_LENGTH, STLIST, LEN_STLIST, STINDEX(1,1,1,im_index), 2,
     &       STASH_LEVELS, NUM_STASH_LEVELS+1,                          
     &       STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                     
     &       im_ident,1,                                                
*CALL ARGPPX                                                            
     &     ICODE, CMESSAGE)                                             
        ENDIF                                                           
C                                                                       
C                                                                       
        IF(SF(258,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(258,1,im_index)+I-1) = SW_NET_SICE(I)          
          ENDDO                                                         
          CALL EXTDIAG (STASHWORK,si(1,1,im_index),SF(1,1),258,258,INT9,
     &       ROW_LENGTH, STLIST, LEN_STLIST, STINDEX(1,1,1,im_index), 2,
     &       STASH_LEVELS, NUM_STASH_LEVELS+1,                          
     &       STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                     
     &       im_ident,1,                                                
*CALL ARGPPX                                                            
     &       ICODE, CMESSAGE)                                           
        ENDIF                                                           
*D ARE2F404.46
        IF ( H_SECT(3) .EQ. '07A' .OR.                                  
     &       H_SECT(3) .EQ. '08A' ) THEN                                
*D ABX1F405.97,ABX1F405.103  
C Set the net SW flux on land tiles                                     
          DO N=1,NTILES                                                 
            DO L=LAND1,LAND1+LAND_PTS-1                                 
              I = LAND_LIST(L)                                          
              J = I + (P_LEVELS + 1 + N)*P_FIELD                        
              SW_TILE(L,N) = RADINCS(J)*COS_ZENITH_ANGLE(I)             
            ENDDO                                                       
          ENDDO                                                         
*D ARE2F404.58,ARE1F405.6    
C Set the surface downward and TOA outward LW fluxes                    
           DO I=FIRST_POINT,LAST_POINT                                  
             J = I + (P_LEVELS + 2)*P_FIELD                             
             LW_DOWN(I) = RADINCS(J+LEN)                                
             DOLR(I) = RADINCS(J+LEN+P_FIELD)                           
*I ARE1F405.8     
                                                                        
C Overwrite SURF_RADFLUX for sea-ice with net SW + downward LW          
           DO I=FIRST_POINT,LAST_POINT                                  
             SURF_RADFLUX(I) = SURF_RADFLUX(I) +
     &                           ICE_FRACT(I)*LW_DOWN(I)                
           ENDDO                                                        
                                                                        
*IF DEF,MPP                                                             
          CALL SWAPB_LAND(SW_TILE,LAND_FIELD,P_FIELD,                   
     &                    ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,            
     &                    NTILES,LAND_LIST)                             
          CALL SWAPBOUNDS(LW_DOWN,ROW_LENGTH,P_ROWS,                    
     &                    EW_Halo,NS_Halo,1)                            
          CALL SWAPBOUNDS(DOLR,ROW_LENGTH,P_ROWS,                       
     &                    EW_Halo,NS_Halo,1)                            
*ENDIF                                                                  
*D ABX1F405.110
     &       //'encountered in CLDCTL.'                                 
*/-----
*DECLARE CNTLATM
*/-----
*I CNTLATM.107
     &  L_ESSERY_SNOW,
*I CNTLATM.124
     &  L_ESSERY_SNOW,
*I CNTLATM.154
     &  L_ESSERY_SNOW,
*I CNTLATM.112   
     &   LTLEADS        ,           !  Let leads temp vary if .TRUE.    
*D CNTLATM.136
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR, LTLEADS,           
*D CNTLATM.166
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR, LTLEADS,           
*/-----
*DECLARE COMPT2A
*/-----
*D COMPT2A.62,COMPT2A.66   
     &,FRACN,FRACM                ! WORK Fractions used in the spreading
C                                 !      calculation.
*D COMPT2A.96,COMPT2A.97   

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)                           
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)                           
*D COMPT2A.109
*D COMPT2A.126
*D COMPT2A.148,COMPT2A.149  

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        DENOM = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)           
*D COMPT2A.178,COMPT2A.179  

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)                           
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)                           
*D COMPT2A.191
*D COMPT2A.208
*/-----
*DECLARE CRADINCS
*/-----
*D ARE2F404.69
     &        RADINCS ( (P_FIELD*(P_LEVELS+2+NTILES)+511)/512*512*2 )   
*/-----
*DECLARE C_ROUGH
*/-----
*I ARN1F404.63    
*IF DEF,A03_8A
!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Z0FSEA = roughness length for free convective heat and moisture
!          transport over the sea (m).
!          DUMMY VARIABLE - Only used in 7A boundary layer scheme
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL Z0FSEA,Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0FSEA = 1.3E-3,
     &          Z0HSEA = 4.0E-5,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
!!----------------------------------------------------------------------
*ENDIF
*/-----
*DECLARE CRUNTIMC
*/-----
*I AJX3F405.148
     & ,DZSNOW(10)
*I AJX3F405.149
     & ,DZSNOW

*/-----
*DECLARE DARCY5A
*/-----
*I DARCY5A.24    
     &,                 DWFLUX_DSTHU1,DWFLUX_DSTHU2
*I DARCY5A.43    
!  4.6      2/99     Extended to provide derivatives. Peter Cox
*I DARCY5A.84    
     &,DWFLUX_DSTHU1(NPNTS) ! OUT The rate of change of the explicit
!                           !     flux with STHU1 (kg/m2/s).
     &,DWFLUX_DSTHU2(NPNTS) ! OUT The rate of change of the explicit
!                           !     flux with STHU2 (kg/m2/s).
*I DARCY5A.88    

      REAL
     & DTHK_DTH1,DTHK_DTH2  ! WORK DTHETAK/DTHETA(1:2).
     &,DK_DTH1,DK_DTH2      ! WORK DK/DTHETA(1:2) (kg/m2/s).
     &,PD                   ! WORK Hydraulic potential difference (m).
*I DARCY5A.97    
     &,DK_DTHK(NPNTS)       ! WORK The rate of change of K with THETAK
!                           !      (kg/m2/s).
*I DARCY5A.99    
     &,DPSI_DTH(NPNTS,2)    ! WORK The rate of change of PSI with
!                           !      THETA(1:2) (m).       
*D DARCY5A.121
            DPSI_DTH(I,N)=0.0   
          ELSE
*D DARCY5A.123,DARCY5A.124  
            DPSI_DTH(I,N)=-B(I)*PSI(I,N)/THETA(I,N)
*I DARCY5A.136   
      DTHK_DTH1=DZ2/(DZ1+DZ2)
      DTHK_DTH2=DZ1/(DZ2+DZ1)

*D DARCY5A.140
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K,DK_DTHK
     &,             LTIMER)          
*D DARCY5A.147
        PD=(2.0*(PSI(I,2)-PSI(I,1))/(DZ2+DZ1)+1)                   
        WFLUX(I)=K(I)*PD 

!-----------------------------------------------------------------------
! Calculate the rate of change of WFLUX with respect to the STHU1 and
! STHU2.
!-----------------------------------------------------------------------
        DK_DTH1=DK_DTHK(I)*DTHK_DTH1
        DK_DTH2=DK_DTHK(I)*DTHK_DTH2
        DWFLUX_DSTHU1(I)=DK_DTH1*PD-2*K(I)*DPSI_DTH(I,1)/(DZ1+DZ2)     
        DWFLUX_DSTHU2(I)=DK_DTH2*PD+2*K(I)*DPSI_DTH(I,2)/(DZ1+DZ2)
*/-----
*DECLARE DECAY2A
*/-----
*D DECAY2A.52,DECAY2A.56   
*D DECAY2A.69
*/-----
*DECLARE DESCENT
*/-----
*D DESCENT.4
     + DENOM_MIN                  ! Minimum value for the denominator
C                                 ! of the update equation. Ensures 
C                                 ! that gradient descent does not  
C                                 ! lead to an unstable solution.  
     +,GAMMA_EQ                   ! Inverse timestep for gradient       
*D ABX1F405.1722
      PARAMETER(DENOM_MIN=1.0E-6, GAMMA_EQ = 1.0E-1, ITER_EQ = 10)      
*/-----
*DECLARE DFPLN3A
*/-----
*I DFPLN3A.25    
     &   , T_SOLID, T_SEA, L_CTILE                                      
*D ADB1F401.54,ADB1F401.55   
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
     &   , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                  
*I DFPLN3A.54    
     &   , L_CTILE                                                      
!             Coastal tiling switch                                     
*I DFPLN3A.63    
     &   , T_SOLID(NPD_PROFILE)                                         
!             TEMPERATURES AT SOLID SURFACE                             
     &   , T_SEA(NPD_PROFILE)                                           
!             SURFACE TEMPERATURE OVER OPEN SEA                         
*I DFPLN3A.74    
     &   , FLANDG(NPD_PROFILE)                                          
!             GATHERED LAND FRACTION                                    
*D ADB1F401.61,ADB1F401.64   
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER     
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER    
*D ADB1F401.67
!             FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX            
*I ADB1F401.71    
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!             PLANCK FUNCTION OVER SEA LEADS                            
*D ADB1F401.83,ADB1F401.85   
!             PLANCKIAN OF SOLID SURFACE (GATHERED OVER SOLID POINTS)   
     &   , SEAFRAC(NPD_PROFILE)                                         
!             FRACTION OF OPEN SEA IN GRID-BOX                          
*I DFPLN3A.144   
      IF(L_CTILE)THEN                                                   
!     SOURCE AT THE OPEN SEA SURFACE.                                   
!     CALCULATE OVER ALL POINTS EVEN THOUGH IT IS                       
!     OVERWRITTEN WHERE THERE IS LAND OR SEA-ICE.                       
      DO L=1, N_PROFILE                                                 
         T_RATIO(L)=T_SEA(L)/T_REF_PLANCK                               
         IF(FLANDG(L).GE.1.0.OR.ICE_FRACTION(L).GE.1.0)                 
     &      T_RATIO(L)=T_SOLID(L)/T_REF_PLANCK                          
                                                                        
         PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)                
      ENDDO                                                             
      DO J=N_DEG_FIT-1, 0, -1                                           
         DO L=1, N_PROFILE                                              
            PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)                
     &         +THERMAL_COEFFICIENT(J)                                  
         ENDDO                                                          
      ENDDO                                                             
!                                                                       
C                                                                       
! Initialise to zero:                                                   
      DO LG=1,NPD_PROFILE                                               
        PLANCK_LEADS_SEA(LG)=0.0                                        
      ENDDO                                                             
!     SET THE SOURCE FUNCTION OVER OPEN SEA LEADS.                      
!     DETERMINE THE TEMPERATURE OF THE NON-SEA FRACTION.                
!     CALCULATE THE SOURCE FUNCTION AT POINTS WITH SOLID SURFACE        
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
         PLANCK_LEADS_SEA(LG)=PLANCK_GROUND(LG)                         
         SEAFRAC(LG)=(1.-FLANDG(LG))*(1.0E+00-ICE_FRACTION(LG))         
         T_RATIO(L)=T_SOLID(LG)/T_REF_PLANCK                            
         PLANCK_GROUND_G(L)=THERMAL_COEFFICIENT(N_DEG_FIT)              
      ENDDO                                                             
      DO J=N_DEG_FIT-1, 0, -1                                           
         DO L=1, N_FRAC_SOL_POINT                                       
            PLANCK_GROUND_G(L)=PLANCK_GROUND_G(L)*T_RATIO(L)            
     &         +THERMAL_COEFFICIENT(J)                                  
         ENDDO                                                          
      ENDDO                                                             
!     DETERMINE THE OVERALL PLANCKIAN FUNCTION OF THE SURFACE.          
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
         PLANCK_GROUND(LG)=(1.-SEAFRAC(LG))*PLANCK_GROUND_G(L)          
     &      +PLANCK_LEADS_SEA(LG)*SEAFRAC(LG)                           
      ENDDO                                                             
                                                                        
      ELSE                      !End of L_CTILE                         
                                                                        
*D ADB1F401.102,ADB1F401.103  
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
*D ADB1F401.109
      DO L=1, N_FRAC_SOL_POINT                                          
*D ADB1F401.114
         DO L=1, N_FRAC_SOL_POINT                                       
*D ADB1F401.121,ADB1F401.122  
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
*I ADB1F401.126   
      ENDIF                                                             
                                                                        
*/-----
*DECLARE DIAG3A
*/-----
*D DIAG3A.29
!     Dummy arguments                                                   
*D DIAG3A.31,DIAG3A.32   
     &    N                                                             
!           Length of array                                             
*D DIAG3A.34,DIAG3A.35   
     &    X(N)                                                          
!           Array to be zeroed                                          
*D DIAG3A.37
!     Local variables                                                   
*D DIAG3A.39,DIAG3A.40   
     &    I                                                             
!           loop variable                                               
*D DIAG3A.45
        X(I)=0.0E+00                                                    
*D DIAG3A.72,DIAG3A.78   
     &  , SEA_FLUX                                                      
     &  , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                        
     &  , L_SURF_DOWN_CLR, SURF_DOWN_CLR                                
     &  , L_SURF_UP_CLR, SURF_UP_CLR                                    
     &  , L_FLUX_BELOW_690NM_SURF                                       
     &  , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF                
     &  , L_MOSES_II                                                    
     &  , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        
     &  , NPD_PROFILE                                                   
     & )                                                               
*D DIAG3A.85
!     Dummy arguments                                                   
*D DIAG3A.87
!     Dimensions of arrays                                              
*D DIAG3A.89,DIAG3A.90   
     &    NPD_PROFILE                                                   
!           Maximum number of atmospheric profiles                      
*D DIAG3A.93,DIAG3A.94   
     &    N_PROFILE                                                     
!           Number of atmospheric profiles                              
*D DIAG3A.96
!     Switches for diagnostics:                                         
*D DIAG3A.98,DIAG3A.105  
     &    L_FLUX_BELOW_690NM_SURF                                       
!           Flux below 690nm at surface to be calculated                
     &  , L_MOSES_II                                                   
!           Surface SW fluxes required for MOSES II 
     &  , L_SURFACE_DOWN_FLUX                                           
!           Downward surface flux required                              
     &  , L_SURF_DOWN_CLR                                               
!           Calculate downward clear flux                               
     &  , L_SURF_UP_CLR                                                 
!           Calculate upward clear flux                                 
*D DIAG3A.107
!     Surface fluxes for coupling or diagnostic use                     
*D DIAG3A.109,DIAG3A.118  
     &    SEA_FLUX(NPD_PROFILE)                                         
!           Net downward flux into sea                                  
     &  , SURFACE_DOWN_FLUX(NPD_PROFILE)                                
!           Downward flux at surface                                    
     &  , SURF_DOWN_CLR(NPD_PROFILE)                                    
!           Clear-sky downward flux at surface                          
     &  , SURF_UP_CLR(NPD_PROFILE)                                      
!           Clear-sky upward flux at surface                            
     &  , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                            
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &  , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                          
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &  , SURF_VIS_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam visible flux                   
     &  , SURF_VIS_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse visible flux                       
     &  , SURF_NIR_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam near-infrared flux             
     &  , SURF_NIR_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse near-infrared flux                 
*D DIAG3A.125
        CALL R2_ZERO_1D(N_PROFILE, SURFACE_DOWN_FLUX)                   
*D DIAG3A.129
        CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_CLR)                       
*D DIAG3A.133
        CALL R2_ZERO_1D(N_PROFILE, SURF_UP_CLR)                         
*D DIAG3A.137
        CALL R2_ZERO_1D(N_PROFILE, FLUX_BELOW_690NM_SURF)               
        CALL R2_ZERO_1D(N_PROFILE, FL_SEA_BELOW_690NM_SURF)             
*I DIAG3A.139   
      IF (L_MOSES_II) THEN                                              
        CALL R2_ZERO_1D(N_PROFILE, SURF_VIS_DIR)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_VIS_DIF)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_NIR_DIR)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_NIR_DIF)                        
      ENDIF                                                             
*I ADB1F401.131   
!       4.6             07-05-99                Calculation of fluxes   
!                                               into the sea changed    
!                                               to use albedos for      
!                                               sea. Code for obsolete  
!                                               net solvers removed.    
!                                               (J. M. Edwards)         
!       5.3             22-02-02    Add two new "blue band" fluxes.     
!                                               (N. Gedney)             
*D DIAG3A.163,DIAG3A.176  
      SUBROUTINE R2_COUPLE_DIAG(N_PROFILE, ISOLIR                       
     &  , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                           
     &  , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                               
     &  , FLANDG, ICE_FRACTION                                          
     &  , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                           
     &  , PLANCK_AIR_SURFACE, THERMAL_SOURCE_GROUND                     
     &  , FLUX_DOWN, FLUX_UP, FLUX_DIRECT                               
     &  , FLUX_DOWN_CLEAR, FLUX_UP_CLEAR, FLUX_DIRECT_CLEAR             
     &  , WEIGHT_690NM                                                  
     &  , SEA_FLUX                                                      
     &  , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                        
     &  , L_SURF_DOWN_CLR, SURF_DOWN_CLR                                
     &  , L_SURF_UP_CLR, SURF_UP_CLR                                    
     &  , L_FLUX_BELOW_690NM_SURF                                       
     &  , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF                
     &  , L_MOSES_II, L_CTILE                                           
     &  , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        
     &  , NPD_PROFILE                                                   
     &  )                                                               
*D DIAG3A.183,DIAG3A.184  
!     Comdecks included                                                 
!     Spectral regions                                                  
*D DIAG3A.187
!     Dummy Arguments                                                   
*D DIAG3A.189
!     Dimensions of arrays                                              
*D DIAG3A.191,DIAG3A.192  
     &    NPD_PROFILE                                                   
!           Maximum number of atmospheric profiles                      
*D DIAG3A.195,DIAG3A.198  
     &    N_PROFILE                                                     
!           Number of atmospheric profiles                              
     &  , ISOLIR                                                        
!           Spectral region                                             
*D DIAG3A.200
!     Logical switches for the code                                     
*D DIAG3A.202,DIAG3A.203  
     &    L_NET                                                         
!           Flag for net fluxes                                         
*D DIAG3A.205
!     Switches for diagnostics:                                         
*D DIAG3A.207,DIAG3A.214  
     &    L_FLUX_BELOW_690NM_SURF                                       
!           Flux below 690nm at surface to be calculated                
     &  , L_MOSES_II                                                    
!           Surface SW required for MOSES II                            
     &  , L_CTILE                                                       
!           Switch for coastal tiling                                   
     &  , L_SURFACE_DOWN_FLUX                                           
!           Downward surface flux required                              
     &  , L_SURF_DOWN_CLR                                               
!           Calculate downward clear flux                               
     &  , L_SURF_UP_CLR                                                 
!           Calculate upward clear flux                                 
*D DIAG3A.216
!     Albedos                                                           
*D DIAG3A.218,DIAG3A.225  
     &    ALBEDO_FIELD_DIFF(NPD_PROFILE)                                
!           Diffuse albedo meaned over grid box                         
     &  , ALBEDO_FIELD_DIR(NPD_PROFILE)                                 
!           Direct albedo meaned over grid box                          
     &  , ALBEDO_SEA_DIFF(NPD_PROFILE)                                  
!           Diffuse albedo of open sea                                  
     &  , ALBEDO_SEA_DIR(NPD_PROFILE)                                   
!           Direct albedo meaned of open sea                            
*D DIAG3A.228,ADB1F401.136  
     &    THERMAL_SOURCE_GROUND(NPD_PROFILE)                            
!           Thermal source at ground                                    
     &  , PLANCK_AIR_SURFACE(NPD_PROFILE)                               
!           Planck function at near-surface air temperature in band     
*D ADB1F401.137,ADB1F401.142  
!     Arguments relating to sea ice.                                    
*D ADB1F401.144,ADB1F401.148  
     &    PLANCK_FREEZE_SEA                                             
!           Planck function over freezing sea                           
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!           Planck function over sea leads                              
     &   , FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
     &   , ICE_FRACTION(NPD_PROFILE)                                    
!            FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX!            
*D DIAG3A.232,DIAG3A.233  
     &    WEIGHT_690NM                                                  
!           Weighting applied to band for region below 690 nm           
*D DIAG3A.235
!     Calculated fluxes                                                 
*D DIAG3A.237,DIAG3A.248  
     &    FLUX_DOWN(NPD_PROFILE)                                        
!           Total downward or net flux at surface                       
     &  , FLUX_DIRECT(NPD_PROFILE)                                      
!           Direct solar flux at surface                                
     &  , FLUX_UP(NPD_PROFILE)                                          
!           Upward flux at surface                                      
     &  , FLUX_DOWN_CLEAR(NPD_PROFILE)                                  
!           Total clear-sky downward or net flux at surface             
     &  , FLUX_UP_CLEAR(NPD_PROFILE)                                    
!           Clear-sky upward flux at surface                            
     &  , FLUX_DIRECT_CLEAR(NPD_PROFILE)                                
!           Clear-sky direct solar flux at surface                      
*D DIAG3A.251
!     Surface fluxes for coupling or diagnostic use                     
*D DIAG3A.253,DIAG3A.262  
     &    SEA_FLUX(NPD_PROFILE)                                         
!           Net downward flux into sea                                  
     &  , SURFACE_DOWN_FLUX(NPD_PROFILE)                                
!           Downward flux at surface                                    
     &  , SURF_DOWN_CLR(NPD_PROFILE)                                    
!           Clear-sky downward flux at surface                          
     &  , SURF_UP_CLR(NPD_PROFILE)                                      
!           Clear-sky upward flux at surface                            
     &  , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                            
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &  , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                          
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &  , SURF_VIS_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam visible flux                   
     &  , SURF_VIS_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse visible flux                       
     &  , SURF_NIR_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam near-infrared flux             
     &  , SURF_NIR_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse near-infrared flux                 
*D DIAG3A.265
!     Local variables                                                   
*D DIAG3A.267,DIAG3A.268  
     &    L                                                             
!           Loop variable                                               
*D DIAG3A.272,DIAG3A.313  
!     This is the flux into the sea over the ice-free parts of the      
!     grid-box. The model is that the downward fluxes are uniform       
!     across the grid-box, but the upward fluxes are not. At this       
!     stage no weighting by the actual ice-free area is carried out:    
!     that will be done later.                                          
      IF (ISOLIR.EQ.IP_SOLAR) THEN                                      
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN         
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))                
     &        +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))    
          ENDIF                                                         
        ENDDO                                                           
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                             
        IF(L_CTILE)THEN                                                 
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN         
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +(1.0E+00-ALBEDO_SEA_DIFF(L))                             
     &        *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      
     &        -PLANCK_LEADS_SEA(L))                                     
          ENDIF                                                         
        ENDDO                                                           
        ELSE                                                            
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0) THEN                                    
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +(1.0E+00-ALBEDO_SEA_DIFF(L))                             
     &        *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      
     &        -PLANCK_FREEZE_SEA)                                       
          ENDIF                                                         
        ENDDO                                                           
        ENDIF                                                           
*D DIAG3A.317,DIAG3A.338  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   
     &        +FLUX_DOWN(L)                                             
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   
     &        +FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                       
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.342,DIAG3A.365  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           
     &        +FLUX_DOWN_CLEAR(L)                                       
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           
     &        +FLUX_DOWN_CLEAR(L)+PLANCK_AIR_SURFACE(L)                 
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.369,DIAG3A.393  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               
     &        +FLUX_UP_CLEAR(L)                                         
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               
     &        +FLUX_UP_CLEAR(L)+PLANCK_AIR_SURFACE(L)                   
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.396
!     This diagnostic is available only in the solar region. Over       
!     sea-points it refers only to the flux over the open sea           
!     (see the comments about sea_flux above). Over land, both the      
!     upward and downward fluxes are taken as uniform.                  
*D DIAG3A.398,DIAG3A.406  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
         IF (L_CTILE) THEN                                              
          DO L=1, N_PROFILE                                             
            FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)           
     &          +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))                 
            IF(FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN        
              FL_SEA_BELOW_690NM_SURF(L)                                
     &          =FL_SEA_BELOW_690NM_SURF(L)                             
     &          +WEIGHT_690NM                                           
     &          *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             
     &          +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)                     
     &          -ALBEDO_SEA_DIR(L)))                                    
            ENDIF                                                       
          ENDDO                                                         
         ELSE                                                           
          DO L=1, N_PROFILE                                             
            IF (FLANDG(L).GT.0.0) THEN                                  
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         
     &          +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))                 
*D DIAG3A.408,DIAG3A.411  
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         
     &          +WEIGHT_690NM                                           
     &          *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             
     &          +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))) 
*I DIAG3A.412   
          ENDDO                                                         
*I DIAG3A.413   
        ENDIF                                                           
*I DIAG3A.415   
!     SURFACE SHORTWAVE DIAGNOSTICS REQUIRED FOR MOSES II  
      IF (L_MOSES_II) THEN                                              
        IF (ISOLIR.EQ.IP_SOLAR) THEN                   
          DO L=1, N_PROFILE                                             
            SURF_VIS_DIR(L) = SURF_VIS_DIR(L) +                         
     &                        WEIGHT_690NM*FLUX_DIRECT(L)               
            SURF_NIR_DIR(L) = SURF_NIR_DIR(L) +                         
     &                        (1. - WEIGHT_690NM)*FLUX_DIRECT(L)        
            SURF_VIS_DIF(L) = SURF_VIS_DIF(L) +                         
     &                        WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_DIRECT(L))
            SURF_NIR_DIF(L) = SURF_NIR_DIF(L) +                         
     &                 (1. - WEIGHT_690NM)*(FLUX_DOWN(L)-FLUX_DIRECT(L))
          ENDDO                                                         
        ENDIF                                                           
      ENDIF                                                             
*/-----
*DECLARE EXFXUV5B
*/-----
*D ACB1F405.4
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D EXFXUV5B.138,EXFXUV5B.149  
*/-----
*DECLARE FCDCH6A
*/-----
*D FCDCH6A.2
*IF DEF,A03_8A
*D FCDCH6A.21
!!!   SUBROUTINES FCDCH_SEA AND FCDCH_LAND-----------------------------
*I FCDCH6A.38    

!     SUBROUTINE FCDCH_SEA--------------------------------------------- 
!                                                                       
!     Transfer coefficients for sea, sea-ice and leads                  
!                                                                       
!     ----------------------------------------------------------------- 
*D FCDCH6A.40,FCDCH6A.43   
      SUBROUTINE FCDCH_SEA(
     & P_POINTS,P_FIELD,P1,LAND_MASK,
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,
     & CDV,CHV,V_S,RECIP_L_MO,LTIMER
*D FCDCH6A.54
*D FCDCH6A.70,FCDCH6A.72   
*D FCDCH6A.79,FCDCH6A.81   
*D FCDCH6A.84,FCDCH6A.86   
*I FCDCH6A.89    

      REAL
     & RIB(P_FIELD)  ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(P_FIELD)  ! DUMMY Used in 7A boundary layer scheme
*D FCDCH6A.104
      EXTERNAL TIMER,PHI_M_H_SEA
*D FCDCH6A.124
      INTEGER I     ! Loop counter; horizontal field index.
*D FCDCH6A.130
*D FCDCH6A.142
*D FCDCH6A.144,FCDCH6A.152  
      DO I=P1,P1+P_POINTS-1                                             
        IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.164
        ENDIF  ! LAND_MASK
*D FCDCH6A.166

      CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
*D FCDCH6A.171,FCDCH6A.181  
      DO I=P1,P1+P_POINTS-1                                             
        IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.186
            V_S(I) = BETA *
*D FCDCH6A.188
*D FCDCH6A.194
*D FCDCH6A.196
          CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
*D FCDCH6A.198,FCDCH6A.200  
        ENDIF  ! LAND_MASK
*D FCDCH6A.207,FCDCH6A.217  
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.220
*D FCDCH6A.224
*D FCDCH6A.227
*D FCDCH6A.231
          ENDIF  ! LAND_MASK
*D FCDCH6A.233
        CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
*D FCDCH6A.238
*D FCDCH6A.240,FCDCH6A.249  
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
*D FCDCH6A.251,FCDCH6A.253  
          ENDIF  ! LAND_MASK
*I FCDCH6A.255   

!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        IF ( .NOT. LAND_MASK(I) ) THEN
          CDV(I) = CDV(I) / VSHR(I)
          CHV(I) = CHV(I) / VSHR(I)
        ENDIF  ! LAND_MASK
      ENDDO


      IF (LTIMER) THEN                                                  
        CALL TIMER('FCDCH   ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               


!!  Arguments:--------------------------------------------------------- 
      SUBROUTINE FCDCH_LAND(
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,      
     & CDV,CHV,CDV_STD,V_S,V_S_STD,RECIP_L_MO,LTIMER                    
     &)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      LOGICAL                                                           
     & LTIMER                                                           
                                                                        
                                                                        
      REAL                                                              
     & DB(LAND_FIELD)! IN Buoyancy difference between surface and lowest
!                    !    temperature and humidity level in the         
!                    !    atmosphere (m/s^2).                           
     &,VSHR(P_FIELD) ! IN Wind speed difference between the surface and 
!                    !    the lowest wind level in the atmosphere (m/s).
     +,Z0M(LAND_FIELD)! IN Roughness length for momentum transport (m).
     +,Z0H(LAND_FIELD)! IN Roughness length for heat and moisture (m).
     +,ZH(P_FIELD)   ! IN Depth of boundary layer (m).                  
     +,Z1_UV(P_FIELD)! IN Height of lowest wind level (m).              
     +,Z1_TQ(P_FIELD)! IN Height of lowest temperature and humidity     
!                    !    level (m).                                    
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                    ! IN for adjusting the surface transfer            
!                    !    coefficients to remove form drag effects.     
                                                                        
      REAL                                                              
     & CDV(LAND_FIELD)! OUT Surface transfer coefficient for momentum
!                    !     including orographic form drag (m/s).        
     +,CHV(LAND_FIELD)! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).        
     &,CDV_STD(LAND_FIELD)
!                    ! OUT Surface transfer coefficient for momentum    
!                    !     excluding orographic form drag (m/s).        
     &,V_S(LAND_FIELD)! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).        
     &,V_S_STD(LAND_FIELD)
!                    ! OUT Surface layer scaling velocity               
!                    !     excluding orographic form drag (m/s).        
     &,RECIP_L_MO(LAND_FIELD)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).

      REAL
     & RIB(LAND_FIELD)  ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(LAND_FIELD)  ! DUMMY Used in 7A boundary layer scheme

                                                                        
!*L  Workspace usage----------------------------------------------------
!                                                                       
!     Local work arrays.                                                
!                                                                       
      REAL                                                              
     & PHI_M(LAND_FIELD)! Monin-Obukhov stability function for momentum
!                     ! integrated to the model's lowest wind level.    
     &,PHI_H(LAND_FIELD)! Monin-Obukhov stability function for scalars
!                     ! integrated to the model's lowest temperature    
!                     ! and humidity level.                             
!                                                                       
!*----------------------------------------------------------------------
                                                                        
      EXTERNAL TIMER,PHI_M_H_LAND
                                                                        
!*----------------------------------------------------------------------
!  Common and local constants.                                          
*CALL C_VKMAN                                                           
      REAL BETA,THIRD                                                   
      PARAMETER (                                                       
     & BETA=0.08,   ! Tunable parameter in the surface layer scaling    
!                   ! velocity formula (multiplying the turbulent       
!                   ! convective scaling velocity).                     
     + THIRD=1./3.  ! One third.                                        
     +)                                                                 
      INTEGER N_ITS ! Number of iterations for Monin-Obukhov length     
!                   ! and stability functions.                          
      PARAMETER (                                                       
     & N_ITS=5                                                          
     &)                                                                 
!                                                                       
!  Define local variables                                               
!                                                                       
      INTEGER I,J,L ! Loop counter; horizontal field index.             
      INTEGER IT    ! Iteration loop counter.                           
                                                                        
      REAL                                                              
     & B_FLUX       ! Surface bouyancy flux over air density.           
     &,U_S          ! Surface friction velocity (effective value).      
     &,U_S_STD      ! Surface friction velocity (standard value).       
     &,W_S          ! Surface turbulent convective scaling velocity.    
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('FCDCH   ',3)                                        
      ENDIF                                                             
                                                                        
!                                                                       
!-----------------------------------------------------------------------
!! 1. Set initial values for the iteration.                             
!-----------------------------------------------------------------------
!                                                                       
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        IF (DB(L) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = -VKMAN/(BETA*BETA*BETA*ZH(I))
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = 0.0
        ENDIF
        ENDDO
      CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                   TILE_INDEX,LAND_INDEX,
     &                   RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &                   PHI_M,PHI_H,
     &                   LTIMER)
!                                                                       
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        IF (DB(L) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          V_S_STD(L) = BETA *
     &        SQRT( BETA * ( VKMAN / PHI_H(L) ) * ZH(I) * (-DB(L)) )
          V_S(L) = V_S_STD(L)
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          V_S(L) = ( VKMAN / PHI_M(L) ) * VSHR(I)
          V_S_STD(L) = V_S(L) * WIND_PROFILE_FACTOR(L)
        ENDIF
        CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
        CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
        CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *
     &                        WIND_PROFILE_FACTOR(L)
      ENDDO                                                             
!-----------------------------------------------------------------------
!! 2. Iterate to obtain sucessively better approximations for CD & CH.  
!-----------------------------------------------------------------------
      DO IT = 1,N_ITS                                                   
!                                                                       
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          B_FLUX = -CHV(L) * DB(L)
          U_S = SQRT( CDV(L) * VSHR(I) )
          U_S_STD = SQRT( CDV_STD(L) * VSHR(I) )
          IF (DB(L) .LT. 0.0) THEN
            W_S = (ZH(I) * B_FLUX)**THIRD
            V_S(L) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
            V_S_STD(L) = SQRT(U_S_STD*U_S_STD + BETA*BETA*W_S*W_S)
          ELSE
            V_S(L) = U_S
            V_S_STD(L) = U_S_STD
          ENDIF
          RECIP_L_MO(L) = -VKMAN * B_FLUX /
     &                     (V_S(L)*V_S(L)*V_S(L))
        ENDDO                                                           
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
     &                     RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &                     PHI_M,PHI_H,
     &                     LTIMER)
!                                                                       
                                                                        
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
          CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
          CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *
     &                          WIND_PROFILE_FACTOR(L)
        ENDDO                                                           
      ENDDO ! Iteration loop
                                                                        
!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        CDV(L) = CDV(L) / VSHR(I)
        CDV_STD(L) = CDV_STD(L) / VSHR(I)
        CHV(L) = CHV(L) / VSHR(I)
      ENDDO

*/-----
*DECLARE FCDCH7A
*/-----
*D FCDCH7A.46,FCDCH7A.48   
      SUBROUTINE FCDCH_SEA (P_POINTS,P_FIELD,P1,FLANDG,              
     &                      RIB,DB,VSHR,Z0M,Z0H,Z0F,                    
     &                      ZH,Z1_UV,Z1_TQ,CD,CH,                       
     &                      V_S,RECIP_L_MO,LTIMER)                      
*D FCDCH7A.59
*D FCDCH7A.62
     & FLANDG(P_FIELD)                                          
!                         ! IN Land fraction                         
     &,RIB(P_FIELD)       ! IN Bulk Richardson number.                  
*I FCDCH7A.75    
      REAL                                                              
     & DB(P_FIELD)   ! DUMMY Used in 6A boundary layer scheme           
     &,VSHR(P_FIELD) ! DUMMY Used in 6A boundary layer scheme           
     &,ZH(P_FIELD)   ! DUMMY Used in 6A boundary layer scheme           
     &,V_S(P_FIELD)  ! DUMMY Used in 6A boundary layer scheme           
     &,RECIP_L_MO(P_FIELD)                                              
!                    ! DUMMY Used in 6A boundary layer scheme           
                                                                        
                                                                        
*D FCDCH7A.121
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
*D FCDCH7A.176
        ENDIF ! SEA_MASK                                              
*D FCDCH7A.194,FCDCH7A.195  
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,      
     & CD,CH,CD_STD,V_S,V_S_STD,RECIP_L_MO,LTIMER                       
*D FCDCH7A.220,FCDCH7A.221
     &,Z1_UV(LAND_FIELD)     ! IN Height of lowest uv level (m).
     &,Z1_TQ(LAND_FIELD)     ! IN Height of lowest tq level (m).
*I FCDCH7A.229   
                                                                        
      REAL                                                              
     & DB(LAND_FIELD)       ! DUMMY Used in 6A boundary layer scheme    
     &,VSHR(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme    
     &,ZH(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme    
     &,V_S(LAND_FIELD)      ! DUMMY Used in 6A boundary layer scheme    
     &,V_S_STD(LAND_FIELD)  ! DUMMY Used in 6A boundary layer scheme    
     &,RECIP_L_MO(LAND_FIELD)                                           
!                           ! DUMMY Used in 6A boundary layer scheme    

*D FCDCH7A.288,FCDCH7A.289
          ZETAM = LOG( (Z1_UV(L) + Z0M(L)) / Z0M(L) )
          ZETAH = LOG( (Z1_TQ(L) + Z0M(L)) / Z0H(L) )
*D FCDCH7A.326
            RFZ = CZ * SQRT ( Z1_UV(L) / Z0F(L) )
*/-----
*DECLARE FILL3A
*/-----
*I ADB1F402.149   
!       5.3             25-04-01   Gather land, sea and                 
!                                  sea-ice temperatures and             
!                                  land fraction. Replace TFS           
!                                  with general sea temperature.        
!                                       (N. Gedney)                     
*D ADB1F401.205,ADB1F401.206  
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA
     &   , AB, BB, AC, BC, PEXNER, TAC                    
     &   , P, T, T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS               
*I FILL3A.466   
     &   , TSTAR_SOLID(NPD_FIELD)                                       
!             SOLID SURFACE TEMPERATURES                                
     &   , TSTAR_SEA(NPD_FIELD)                                         
!             OPEN SEA SURFACE TEMPERATURES                             
*I ADB1F401.208   
     &   , T_SOLID(NPD_PROFILE)                                         
!             GATHERED TEMPERATURE OF SOLID SURFACE                     
     &   , T_SEA(NPD_PROFILE)                                           
!             GATHERED OPEN SEA TEMPERATURE                             
*I ADB1F401.213   
            T_SOLID(L)=TSTAR_SOLID(LG)                                  
            T_SEA(L)=TSTAR_SEA(LG)                                      
*/-----
*DECLARE FRUNOF7A
*/-----
*D FRUNOF7A.50
     & SURF_ROFF,NSNOW,L_ESSERY_SNOW
*I FRUNOF7A.68
     &,NSNOW(NPNTS)         ! IN 

      LOGICAL L_ESSERY_SNOW,DEEPSNOW
*I FRUNOF7A.88
        !NEW snow scheme to intercept throughfall on deep snow
        !simple infiltration not apprpriate
        DEEPSNOW=.FALSE.
        IF (L_ESSERY_SNOW) THEN
          IF (NSNOW(I).gt.0) DEEPSNOW=.TRUE.
        ENDIF

!!        IF (.NOT.DEEPSNOW) THEN
*D  FRUNOF7A.89
 ! need to limit very small runoff to prevent blowups
         IF (R(I).GT.1E-99) THEN
*I FRUNOF7A.109
!!        ENDIF

*/-----
*DECLARE FTSA1A
*/-----
*D ARE2F404.294,ARE2F404.304  
*D ARE2F404.307,AJG1F405.36   
     &     LAND, FLANDG, AICE,                                          
     &     TSTAR, TSTAR_SICE, SFA, MDSA, COSZ, S,                       
     &     ALPHAM,ALPHAC,ALPHAB,DTICE,L_MOSES_II,L_SSICE_ALBEDO,        
*D ARE2F404.309
     &     L1, L2, SA_LAND, SA_SICE, SAOS)                              
*D ARE2F404.310,ARE2F404.312  
*D ARE2F404.313,ARE2F404.314  
     &    ,L_MOSES_II                     ! .TRUE. if MOSES II land     
!                                         ! surface is selected.        
*I ADB1F400.16    
     &     FLANDG(L1),                    ! Land fraction               
*D ARE2F404.315
*D ARE2F404.316
     &     TSTAR_SICE(L1),                ! Seaice surface temperature  
*D ARE2F404.317,ARE2F404.318  
*D ARE2F404.319,FTSA1A.50   
     &     SA_LAND(L1),                   ! Surface Albedos for Land.   
!                                         ! (Not output for MOSESII).   
     &     SA_SICE(L1),                   ! Surface Albedos for seaice  
*D FTSA1A.58,FTSA1A.60   
      REAL DSA                           ! Deep-snow albedo (alphasubD) 
                                                                        
*D ARE2F404.321,ARE2F404.328  
     &     SNOW_ALBEDO,                   ! Snow albedo                 
*D ARE2F404.329
      PARAMETER ( MASKD = 0.2 )                                         
*I AWA1F304.1407  
      DO J=1, L2                                                        
        SA_LAND(J)=0.0                                                  
        SA_SICE(J)=0.0                                                  
      ENDDO                                                             
*D ARE2F404.330,ARE2F404.334  
! Land surface albedos are calculated by routine TILEALB if MOSES II    
! is selected                                                           
      IF (.NOT. L_MOSES_II) THEN                                        
*D ARE2F404.336
*D ARE2F404.337,ARE2F404.353  
          IF ( TSTAR(J) .LT. TCLAND ) THEN                              
             DSA = MDSA(J)                                              
           ELSE                                                         
            DSA=MDSA(J) + KLAND * (SFA(J)-MDSA(J)) * (TSTAR(J)-TCLAND)  
*D ARE2F404.355,ARE2F404.384  
          SA_LAND(J) = SFA(J) + (DSA-SFA(J)) * ( 1. - EXP(-MASKD*S(J)) )
          SAOS(J,1) = SA_LAND(J)                                        
          SAOS(J,2) = SA_LAND(J)                                        
*D ARE2F404.387,ARE2F404.407  
      ENDIF                                                             
*D ARE2F404.410
        IF (FLANDG(J).LT.1.0) THEN                                      
*D FTSA1A.103
             SA_SICE(J) = 0.                                            
*D FTSA1A.105,FTSA1A.106  
*D AJG1F405.63,AJG1F405.64   
                 if (tstar_sice(j).gt.tcice) then                       
                   snow_albedo=ice2+ice1*tstar_sice(j)                  
*D AJG1F405.68
                 sa_sice(j)=alphab                                      
*D AJG1F405.71
                 sa_sice(j)=alphab                                      
*D FTSA1A.108,FTSA1A.109  
             IF ( TSTAR_SICE(J) .LT. TCICE ) THEN                       
                SA_SICE(J) = ALPHAC                                     
*D FTSA1A.111
                SA_SICE(J) = ICE1 * TSTAR_SICE(J) + ICE2                
*I FTSA1A.114   
                                                                        
*I FTSA1A.115   
                                                                        
*/-----
*DECLARE FXCA3A
*/-----
*I ADB2F404.557   
!       4.6             10-05-98                Land flag added to      
!                                               the argument list for   
!                                               diagnostics.            
!                                               (J. M. Edwards)         
*I FXCA3A.17    
!                                                                       
*D FXCA3A.38
     &   , P, T, T_GROUND, T_SOLID, T_SEA, T_LEVEL, D_MASS              
*D ADB1F401.433,FXCA3A.80   
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG                      
*D ADB1F401.435
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               
     &   , L_MOSES_II, L_CTILE                                          
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       
*I FXCA3A.232   
     &   , T_SOLID(NPD_PROFILE)                                         
!             TEMPERATURE OF SOLID SURFACE                              
     &   , T_SEA(NPD_PROFILE)                                           
!             SURFACE TEMPERATURE OVER OPEN SEA                         
*I ADB1F401.438   
     &   , L_MOSES_II                                                   
!             SURFACE SW REQUIRED FOR MOSES II                          
     &   , L_CTILE                                                      
!             COASTAL TILING SWITCH                                     
*I FXCA3A.458   
                                                                        
      REAL     !, INTENT(IN)                                            
     &     FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
*D ADB1F401.440,ADB1F401.443  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.450
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &   , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                         
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &   , SURF_VIS_DIR(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX                 
     &   , SURF_VIS_DIF(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIFFUSE VISIBLE FLUX                     
     &   , SURF_NIR_DIR(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX           
     &   , SURF_NIR_DIF(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX               
*I ADB1F401.454   
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!             PLANCK FUNCTION OVER SEA LEADS                            
*D ADB1F401.462
     &   , L_FLUX_BELOW_690NM_SURF                                      
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               
     &   , L_MOSES_II                                                   
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       
*I FXCA3A.1090  
     &         , T_SOLID, T_SEA, L_CTILE                                
*D ADB1F401.464,ADB1F401.465  
     &         , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION       
     &         , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA            
*D FXCA3A.1515
         CALL R2_COUPLE_DIAG(N_PROFILE, ISOLIR                          
*D ADB1F401.467,ADB1F401.468  
     &      , FLANDG, ICE_FRACTION                                      
     &      , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                       
*D FXCA3A.1519,ADB1F401.470  
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+2)                           
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+1)                           
*D FXCA3A.1522,ADB1F401.471  
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+2)                     
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+1)                     
*D ADB1F401.473
     &      , L_FLUX_BELOW_690NM_SURF                                   
     &      , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF            
     &      , L_MOSES_II, L_CTILE                                       
     &      , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF    
*D FXCA3A.1549
*/-----
*DECLARE GROWT2A
*/-----
*D ABX1F405.1629
*D GROWT2A.58,GROWT2A.62   
*D GROWT2A.92
*/-----
*DECLARE HTCOND5A
*/-----
*D HTCOND5A.110,HTCOND5A.116
*D HTCOND5A.123
*D HTCOND5A.141
*/-----
*DECLARE HYDCON5A
*/-----
*I HYDCON5A.23    
     &,                   DK_DTHK
*I HYDCON5A.69    
     &,DK_DTHK(NPNTS)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).
*I HYDCON5A.81    
        DK_DTHK(I)=0.0
*I HYDCON5A.83    
          DK_DTHK(I)=(2*B(I)+3)*KS(I)*(THETAK(I)**(2*B(I)+2))     
*/-----
*DECLARE HYDROL7A
*/-----
*D HYDROL7A.24,HYDROL7A.25   
     &                   GBM_OROG,LICE_PTS,LICE_INDEX,SOIL_PTS,
     &                   SOIL_INDEX,NTILES,TILE_PTS,TILE_INDEX,
*D HYDROL7A.26
     &                   NPNTS,NSHYD,NSMAX,NELEV,
     &                   B,CAN_CPY,CON_RAIN,CON_SNOW,
*D HYDROL7A.27,HYDROL7A.35   
     &                   E_CANOPY,EXT,HCAP,HCON,INFIL_TILE,LS_RAIN,     
     &                   LS_SNOW,MELT_TILE,EI_TILE,SATCON,SATHH,
     &                   SURF_HT_FLUX,surf_htf_tile,TSTAR_TILE,FRAC, 
     &                   TIMESTEP,V_SAT,V_WILT,L_SNOW_ALBEDO,
     &                   STF_HF_SNOW_MELT,STF_SUB_SURF_ROFF,            
     &                   CAN_WCNT,RGRAIN,SMCL,SNOW_TILE,STHF,STHU,TSOIL,
     &                   CAN_WCNT_GB,HF_SNOW_MELT,INFIL,LYING_SNOW,SMC,
     &                   SNOW_MELT,TICE,SNOMLT_SUB_HTF,
*D HYDROL7A.36
     &                   SUB_SURF_ROFF,SURF_ROFF,TOT_TFALL,
! New variables/arrays passed into/out of routine, for Essery snow
! scheme:
     &                   SNOWDEPTH, RHO_SNOW_GRND, SNOW_CAN,
     &              SNOW_SOIL_HTF_P, NSNOW, DS, SICE, SLIQ, TSNOWLAYER,
     &                   RHO_SNOW,RGRAINL,L_ESSERY_SNOW,
     &                   dzsnow, refreezing_layer
     &                  , refreezing_tile, refreezing
     &                  , melt_tiles_before_snow_scheme
     &                  , melt_before_snow_scheme 
     &                  , surf_ht_flux_ice, tfall_tile,
     &                   LTIMER
*D HYDROL7A.48,HYDROL7A.51   
*D HYDROL7A.54,HYDROL7A.56   
*D HYDROL7A.101
     &,NTILES              ! IN Number of tiles.                        
*D HYDROL7A.118,HYDROL7A.119  
     &,TILE_PTS(NTILES)    ! IN Number of tile points.                  
     &,TILE_INDEX(NPNTS,NTILES)                                         
*D HYDROL7A.124,HYDROL7A.126  
     &,CAN_CPY(NPNTS,NTILES)!IN Canopy/surface capacity of              
!                          !    land tiles (kg/m2).           
*D HYDROL7A.129,HYDROL7A.130  
     &,E_CANOPY(NPNTS,NTILES)                                          
!                          ! IN Canopy evaporation from
*I HYDROL7A.135   
     &,INFIL_TILE(NPNTS,NTILES)                                         
!                          ! IN Maximum surface infiltration            
*I HYDROL7A.137   
     &,MELT_TILE(NPNTS,NTILES)
!                          ! IN Snowmelt on tiles (kg/m2/s).
*D HYDROL7A.141,HYDROL7A.148  
     &,SURF_HT_FLUX(NPNTS) ! IN Net downward surface heat flux (W/m2) 
     &,TSTAR_TILE(NPNTS,NTILES)
!                          ! IN Surface temperature (K).                
     &,FRAC(NPNTS,NTILES)  ! IN Tile fractions.
*D HYDROL7A.159,HYDROL7A.160  
     & CAN_WCNT(NPNTS,NTILES)                                          
!                          ! INOUT Canopy water content for
*D HYDROL7A.162,HYDROL7A.163  
     &,RGRAIN(NPNTS,NTILES)! INOUT Snow grain size (microns).           
*D HYDROL7A.166
     &,SNOW_TILE(NPNTS,NTILES)
!                          ! INOUT Snowmass on tiles (kg/m2). 
*D HYDROL7A.174
*I HYDROL7A.179   
     &,HF_SNOW_MELT(NPNTS)  ! OUT Gridbox snowmelt heat flux (W/m2).    
*I HYDROL7A.181   
     &,LYING_SNOW(NPNTS)    ! OUT Gridbox snowmass (kg/m2). 
*I HYDROL7A.173
     &,TICE(NPNTS,NSHYD)
*I HYDROL7A.211
      REAL CATCH_SNOW(NPNTS,NTILES)
      REAL SNOW_GRND(NPNTS,NTILES)
      REAL SURF_HTF_TILE(NPNTS,NTILES)
      REAL refreezing_layer(NPNTS,NTILES,NSMAX)
      REAL refreezing_tile(NPNTS,NTILES)
      REAL refreezing(NPNTS)

      REAL melt_tiles_before_snow_scheme(npnts, ntiles)
      REAL snow_tile_old(npnts, ntiles)
      REAL melt_before_snow_scheme(npnts)
      REAL SURF_HT_FLUX_ICE(NPNTS)
      REAL ice_roff(NPNTS,NTILES)
      REAL melt_roff(NPNTS)

!
      INTEGER NSMAX,NELEV
      LOGICAL L_ESSERY_SNOW

      REAL SMCL1(NPNTS)
      REAL STHF1(NPNTS)
      REAL TSOIL1(NPNTS,nelev)
      REAL TICE1(NPNTS,nelev)
      REAL HCONS(npnts)

      REAL dzsnow(nsmax)

! Additional variables for JULES snow scheme
      REAL
     & SNOWDEPTH(    npnts,ntiles),
! Snow depth on ground on tiles (m)
     & RHO_SNOW_GRND(npnts,ntiles),
! Snowpack bulk density (kg/m3)
     & SNOW_CAN(     npnts,ntiles),
! Snow on the canopy (kg/m2)
     & SNOW_SOIL_HTF_P(npnts,ntiles),
! Surface heat flux beneath snow (W/m2)
     & NSNOW(        npnts,ntiles),
! Number of snow layers on ground on tiles
! NOTE that this is converted to an integer.
     & DS(        npnts,ntiles,nsmax),
! Snow layer thickness (m)
     & SICE(      npnts,ntiles,nsmax),
! Snow layer ice mass on tiles (Kg/m2)
     & SLIQ(      npnts,ntiles,nsmax),
! Snow layer liquid mass on tiles (Kg/m2)
     & TSNOWLAYER(npnts,ntiles,nsmax),
! Snow layer temperature (K)
     & RHO_SNOW(  npnts,ntiles,nsmax),
! Snow layer densities (kg/m3)
     & RGRAINL(   npnts,ntiles,nsmax)
     &,EI_TILE(NPNTS,NTILES)
     &,TFALL_TILE(NPNTS,NTILES)

       REAL DUMMY1(NPNTS),DUMMY2(NPNTS)

       REAL elevation(10)
       REAL delevation,gbm_orog(npnts)
       REAL icegrad,soilgrad
       REAL ICE_FRAC(NPNTS),TSTAR_NONICE(NPNTS)

!+seg These are mid-points. Changed to 10 levels (CESM type) from 25
      data elevation / 100., 300., 550., 850., 1150.,
     &                 1450.,1800.,2250.,2750.,3600./
       data soilgrad /-.001/
       data icegrad /-.001/
! Snow layer grain size on tiles (microns)
*B HYDROL7A.218
     & ,SNOW_INTCTL
*I HYDROL7A.224
!
!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------

        DO I=1,NPNTS
          ice_frac(i)=0
          tstar_nonice(i)=0
          DO N=(ntiles-nelev+1),ntiles
        ! we need these to get soil-only and ice-only grid-box-means
            ice_frac(i)=ice_frac(i)+frac(i,n)
          end do
          DO N=1,(ntiles-nelev)
        ! we need this for the temperature of water infiltrating the soil
            if (ice_frac(i) .lt. 1)
     &      tstar_nonice(i)=tstar_nonice(i)
     &                   +tstar_tile(i,n)*frac(i,n)/(1-ice_frac(i))
          end do
        end do


      IF (L_ESSERY_SNOW) THEN
      CALL HEAT_CON (NPNTS,HCON,STHU,STHF,V_SAT,HCONS,LTIMER)
!
        DO I=1,NPNTS
          SMCL1(I) = SMCL(I,1)
          STHF1(I) = STHF(I,1)
          do n=1,nelev
            delevation=elevation(n)-gbm_orog(i)
            TSOIL1(I,n) = TSOIL(I,1)+(delevation*soilgrad)
            TICE1(I,n) = TICE(I,1)+(delevation*icegrad)
          enddo
        ENDDO
        DO I=1,NPNTS
          melt_before_snow_scheme(i) = 0.0
          DO J=1,NTILES

!            snow_tile_old(I,J) = snow_tile(I,j)
            melt_tiles_before_snow_scheme(i,j) = melt_tile(i,j)
            melt_before_snow_scheme(i) = melt_before_snow_scheme(i)
     &                  +  (frac(i,j) * melt_tile(i,j))
          ENDDO
        END DO
!
        CALL SNOW_INTCTL (NPNTS,TIMESTEP,STF_HF_SNOW_MELT,NTILES,
     &                    TILE_PTS, TILE_INDEX,CATCH_SNOW,CON_SNOW,
     &                    FRAC,LS_SNOW,TFALL_TILE,
     &                     HCAP,HCONS,
     &                     MELT_TILE, EI_TILE,SMCL1,
     &                     STHF1,SURF_HTF_TILE,
     &                     TSOIL1, TICE1, TSTAR_TILE,
     &                     V_SAT,
     &                     RGRAIN,
     &                     SNOW_GRND, SNOW_MELT, SNOW_TILE,
     &                     HF_SNOW_MELT, LYING_SNOW,
     &                     SNOMLT_SUB_HTF, 
     &                     l_snow_albedo, dzsnow, refreezing_layer,
     &                     refreezing_tile, refreezing,NELEV,
! Prognostics for Essery snow scheme:
     &          NSMAX,SNOWDEPTH,RHO_SNOW_GRND,SNOW_CAN,SNOW_SOIL_HTF_P,
     &           NSNOW,DS,SICE,SLIQ,TSNOWLAYER,RHO_SNOW,RGRAINL )
!

      ELSE
*I HYDROL7A.230
      ENDIF
*/-----


*D HYDROL7A.202,HYDROL7A.208  
*D HYDROL7A.212,HYDROL7A.214  
*D HYDROL7A.217
     & SFSNOW,SURF_HYD,SOIL_HYD,SOIL_HTC,ICE_HTC,SOILMC
*D HYDROL7A.226
! Add snowfall to snow mass and update snow grain size            
*D HYDROL7A.228
      CALL SFSNOW(NPNTS,NTILES,TILE_PTS,TILE_INDEX,
     &            CON_SNOW,LS_SNOW,FRAC,TSTAR_TILE,TIMESTEP,
     &            RGRAIN,SNOW_TILE,L_SNOW_ALBEDO,LTIMER)
*D HYDROL7A.231
! Calculate the gridbox-mean snow mass and melt rate 
*D HYDROL7A.233,HYDROL7A.237  
      DO I=1,NPNTS  
        LYING_SNOW(I) = 0.
        SNOW_MELT(I) = 0. 
        SURF_HT_FLUX(I)= 0.
        SURF_HT_FLUX_ICE(I)= 0.
      ENDDO               
! just for soil, hydrology etc.
      DO N=1,(ntiles-nelev)

        DO J=1,TILE_PTS(N)                                              
          I = TILE_INDEX(J,N)

           if (ice_frac(i).lt.1) then !REALLY SHOULD

            if (L_ESSERY_SNOW) THEN
          SURF_HT_FLUX(I)=SURF_HT_FLUX(I)+
     &                  SNOW_SOIL_HTF_P(I,N)*(FRAC(I,N)/(1-ice_frac(i)))
            else
          SURF_HT_FLUX(I)=SURF_HT_FLUX(I)+
     &                  SURF_HTF_TILE(I,N)*(FRAC(I,N)/(1-ice_frac(i)))
            endif
            endif
           ! NOT NORMALISED FOR SOIL_ONLY, FULL GBM REQD
          LYING_SNOW(I) = LYING_SNOW(I) + FRAC(I,N)*SNOW_TILE(I,N)
          if (ice_frac(i).lt.1) !REALLY SHOULD
     &     SNOW_MELT(I) = SNOW_MELT(I)+
     &                    (MELT_TILE(I,N)*FRAC(I,N)/(1-ice_frac(i)))

        END DO

      END DO

      CALL SURF_HYD (NPNTS,NTILES,NELEV,TILE_PTS,TILE_INDEX,
     &               CAN_CPY,E_CANOPY,FRAC,INFIL_TILE,CON_RAIN,LS_RAIN,
     &               MELT_TILE,SNOW_MELT,TIMESTEP,
     &               CAN_WCNT,CAN_WCNT_GB,DSMC_DT,SURF_ROFF,TOT_TFALL
     &               ,tfall_tile,nsnow,l_essery_snow
     &               )

      DO I=1,NPNTS
        !put snow_melt,surf_roff back from soil-only avg to full GBM
        !for diagnostics
        !now that dsmc_dt (soil only) is done
        SNOW_MELT(I) = SNOW_MELT(I)*(1-ice_frac(i))
        SURF_ROFF(I) = SURF_ROFF(I)*(1-ice_frac(i))
      END DO

      !ice separate, different subsurface.
      DO N=(ntiles-nelev+1),NTILES
        DO J=1,TILE_PTS(N)                                              
          I = TILE_INDEX(J,N)

          LYING_SNOW(I) = LYING_SNOW(I) + FRAC(I,N)*SNOW_TILE(I,N)
          SNOW_MELT(I) = SNOW_MELT(I)+(MELT_TILE(I,N)*FRAC(I,N))
          SURF_ROFF(I)=SURF_ROFF(I)+melt_tile(I,n)*FRAC(I,N)

            IF (ice_frac(i).gt.0) THEN !REALLY SHOULD BE IF WE'RE HERE

            if (L_ESSERY_SNOW) then

             SURF_HT_FLUX_ICE(I)=SURF_HT_FLUX_ICE(I)+
     &                        SNOW_SOIL_HTF_P(I,N)*FRAC(I,N)/ice_frac(i)
            else

             SURF_HT_FLUX_ICE(I)=SURF_HT_FLUX_ICE(I)+
     &                        SURF_HTF_TILE(I,N)*FRAC(I,N)/ice_frac(i)
            endif
            ENDIF

        ENDDO                                                           
      ENDDO  

*I HYDROL7A.244   
          SNOMLT_SUB_HTF(I) = 0. 
*D HYDROL7A.247,HYDROL7A.256  
*D HYDROL7A.258,HYDROL7A.266
*I HYDROL7A.271   
      IF (SOIL_PTS.NE.0) THEN 
*D HYDROL7A.274
     &               SUB_SURF_ROFF,SMCL,STHU,SURF_ROFF,W_FLUX,
*I HYDROL7A.275   

      ELSE 

!-----------------------------------------------------------------------
! If required by STASH flag and there are no soil points, 
! set sub-surface runoff to zero.
!-----------------------------------------------------------------------
        IF(STF_SUB_SURF_ROFF) THEN   
          DO I=1,NPNTS  
            SUB_SURF_ROFF(I)=0.0
          ENDDO 
        ENDIF   

      ENDIF    
*D HYDROL7A.281,HYDROL7A.284  
        CALL SOIL_HTC (NPNTS,NSHYD,NTILES,NELEV,SOIL_PTS,SOIL_INDEX,
     &                 TILE_PTS,TILE_INDEX,B,DZSOIL,FRAC,HCAP,HCON,
     &                 SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP,V_SAT,   
     &                 W_FLUX,SMCL,STHU,STHF,TSOIL,NSNOW,L_ESSERY_SNOW,
     &                 TSTAR_NONICE,LTIMER)
*D HYDROL7A.292,HYDROL7A.293
     &                surf_ht_flux_ice,TIMESTEP,                           
     &                TICE,L_ESSERY_SNOW,LTIMER)
*D HYDROL7A.294,HYDROL7A.298  
*/-----
*DECLARE HYDR_CT1
*/-----
*D AJS1F401.757
      SUBROUTINE HYDR_CTL(SNOW_SUBLIMATION,EI_TILE,SNOW_MELT,
*D AJS1F401.759
     &           SURF_HT_FLUX,SURF_HTF_TILE,LS_RAIN,LS_SNOW,
*I HYDR_CT1.62
*CALL CRUNTIMC
      REAL refreezing_layer(land_field, ntiles*nsmax)
      REAL refreezing_tile(land_field, ntiles)
      REAL refreezing(land_field)
      REAL melt_tiles_before_snow_scheme(land_field, ntiles)
      REAL melt_before_snow_scheme(land_field)
      REAL surf_htf_tile(land_field, ntiles)
*I AJS1F401.776
     &      ,SURF_HT_FLUX_ICE(LAND_FIELDDA)
*I ACB1F304.16
     &     ,SUMIFRAC
*I HYDR_CT1.89 
     &     ,N
*I AJS1F401.825
     &  NSMAX,NELEV,
*B ACB1F304.21
      INTEGER PESLEVEL,ESLEVEL
      REAL EI_TILE(TILE_FIELDDA,NTILESDA),GBM_OROG_LAND(LAND_FIELDDA)
      REAL TFALL_TILE(TILE_FIELDDA,NTILESDA)
*I ACB1F304.22
     &,     PLLSNOW(NTILES*NSMAX) ! pseudolevel list for Essery snow tiles

*I GDR6F405.101
        GBM_OROG_LAND(I)=D1(JOROG+LAND_LIST(I)-1)
*D GDR6F405.113
     &   SURF_HT_FLUX_LAND,surf_htf_tile,GBM_OROG_LAND,
     &   D1(JVEG_FRAC),D1(JVOL_SMC_SAT),
*I AJS1F401.835
     &   EI_TILE,
*I AJS1F401.841
     &   D1(J_DEEP_ICE_TEMP(1)),
*I ARE1F404.147
! New variables/arrays passed into/out of routine, for Essery snow
! scheme:
     &  D1(JSNOWDEPTH(1)), D1(JRHO_SNOW_GRND(1)), D1(JSNOW_CAN(1)),
     &  D1(JSNOW_SOIL_HTF(1)), D1(JNSNOW(1)), D1(JDS(1)), D1(JSICE(1)),
     &  D1(JSLIQ(1)), D1(JTSNOWLAYER(1)), D1(JRHO_SNOW(1)),
     &  D1(JRGRAINL(1)), L_ESSERY_SNOW,dzsnow, refreezing_layer, 
     &  refreezing_tile, refreezing,
     &  melt_tiles_before_snow_scheme,
     &  melt_before_snow_scheme,
     &  surf_ht_flux_ice,tfall_tile,
*D ARE1F404.140,AJS1F401.807
      DO I=LAND1,LAND1+LAND_PTS-1
      !sum over ice elevations to determine whether this point has
      !soil and/or ice. Precision tolerances need to be set to avoid
      !v small false positives

        SUMIFRAC=0.
        DO N=(NTYPE_1-1)*NELEV,NTYPE_1*NELEV-1
          SUMIFRAC=SUMIFRAC+D1(JFRAC_TYP+(N*LAND_FIELD)+I-1)
        ENDDO

        IF (SUMIFRAC .GT. 1e-6) THEN
          LICE_PTS=LICE_PTS+1
          LICE_INDEX(LICE_PTS)=I
        ENDIF
        IF (SUMIFRAC .LT. (1.-1e-6)) THEN
          SOIL_PTS=SOIL_PTS+1
          SOIL_INDEX(SOIL_PTS)=I
        END IF

      ENDDO
*I AJS1F400.313
CL ITEM 432 DEEP ICE TEMPERATURES

      IF (SF(432,8)) THEN
        CALL SET_LEVELS_LIST(ST_LEVELS,LEN_STLIST,
     &       STLIST(1,STINDEX(1,432,8,im_index)),
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,ST_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
            CALL FROM_LAND_POINTS(STASHWORK(SI(432,8,im_index)
     &          +(LEVEL_OUT-1) *P_FIELD),D1(J_DEEP_ICE_TEMP(LEVEL)),
     &          LD1(JLAND),P_FIELD,LAND_FIELD)

          END IF
        END DO
      END IF 
      IF(SF(433,8)) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=1,LAND_FIELD
          STASHWORK(si(433,8,im_index)+LAND_LIST(I)-1)=
     &    surf_ht_flux_land(i)
        END DO
      END IF
      IF(SF(434,8)) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=1,LAND_FIELD
          STASHWORK(si(434,8,im_index)+LAND_LIST(I)-1)=
     &    surf_ht_flux_ice(i)
        END DO
      END IF
*B GPB8F405.28

CL ITEM 435 CANOPY THROUGHFALL ON TILES                             
                                                                        
      IF (SF(435,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,435,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(435,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),tfall_tile(1,PSLEVEL),    
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF                                                        
        END DO                                                          
      END IF    

! Copy snow prognostics to diagnostic fields STASH can cope with
!
! Item 376 (Snow depth on ground on tiles (m)):
      IF (SF(376,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,376,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(376,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),D1(JSNOWDEPTH(1)+((PSLEVEL-1)*LAND_FIELD)),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF
C ITEM 377 - Essery scheme SNOWPACK BULK DENSITY
      IF (SF(377,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,377,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(377,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),D1(JRHO_SNOW_GRND(1)
     &           +((PSLEVEL-1)*LAND_FIELD)),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF

C ITEM 378 - Essery scheme SNOW ON THE CANOPY
      IF (SF(378,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,378,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(378,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),D1(JSNOW_CAN(1)
     &           +((PSLEVEL-1)*LAND_FIELD)),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF

C ITEM 379 - Essery scheme SURFACE HEAT FLUX UNDER SNOW
      IF (SF(379,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,379,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(379,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),D1(JSNOW_SOIL_HTF(1)
     &           +((PSLEVEL-1)*LAND_FIELD)),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF

C ITEM 380 - Essery scheme NUMBER OF SNOW LAYERS ON TILES
      IF (SF(380,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,380,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(380,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD),D1(JNSNOW(1)
     &           +((PSLEVEL-1)*LAND_FIELD)),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF

C ITEM 382 - Essery scheme SNOW LYR ICE MASS ON TILES(KG M-2)
C HARDCODE LAYER DIAGNOSTICS TO JUST LAND ICE, PS_LEVEL SPEC
C CAN'T COPE WITH ALL POSS LEVELS, DOESN'T WANT TO DO SUBSET
      PESLEVEL=0
      PSLEVEL_OUT=0
      DO ESLEVEL=1,NSMAX  !ALL LAYERS
      DO PSLEVEL=(NTILES-NELEV+1),NTILES !JUST ICE TILES
        PESLEVEL=PSLEVEL+(ESLEVEL-1)*NTILES !1D COUNTER FOR STEPPING THROUGH D1
        PSLEVEL_OUT=PSLEVEL_OUT+1 
        do I=1,LAND_FIELD
      IF (SF(381,8)) STASHWORK(SI(381,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JDS(1)+(PESLEVEL-1)*LAND_FIELD+I-1)

      IF (SF(382,8)) STASHWORK(SI(382,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JSICE(1)+(PESLEVEL-1)*LAND_FIELD+I-1)

      IF (SF(383,8)) STASHWORK(SI(383,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JSLIQ(1)+(PESLEVEL-1)*LAND_FIELD+I-1)

      IF (SF(384,8)) STASHWORK(SI(384,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JTSNOWLAYER(1)+(PESLEVEL-1)*LAND_FIELD+I-1)

      IF (SF(385,8)) STASHWORK(SI(385,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JRHO_SNOW(1)+(PESLEVEL-1)*LAND_FIELD+I-1)

      IF (SF(386,8)) STASHWORK(SI(386,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = D1(JRGRAINL(1)+(PESLEVEL-1)*LAND_FIELD+I-1)
      IF (SF(390,8)) STASHWORK(SI(390,8,im_index)+
     &          (PSLEVEL_OUT-1)*P_FIELD+LAND_LIST(i)-1
     &                                    )
     &        = refreezing_layer(I,PESLEVEL)
        enddo
      END DO
      END DO


! Item 391 (Refreezing on tiles (kg/m2/s)):
      IF (SF(391,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,391,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(391,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD), refreezing_tile(:,pslevel),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF
!
! Item 392 (Gridbox total refreezing (kg/m2/s)):
      IF (SF(392,8)) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=1,LAND_FIELD
          STASHWORK(SI(392,8,im_index)+LAND_LIST(I)-1) =
     &    refreezing(i)
        ENDDO
      ENDIF
! Item 393 (Snow melt on tiles before call to snow scheme(kg/m2/s)):
      IF (SF(393,8)) THEN
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,
     &       STLIST(1,STINDEX(1,393,8,im_index)),
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,
     &       ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          RETURN
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL FROM_LAND_POINTS (
     &          STASHWORK(SI(393,8,im_index)+(PSLEVEL_OUT-1)
     &           *P_FIELD), melt_tiles_before_snow_scheme(:,pslevel),
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF
        END DO
      END IF
! Item 394 (Gridbox total snow melt before call to snow scheme (kg/m2/s)):
      IF (SF(394,8)) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=1,LAND_FIELD
          STASHWORK(SI(394,8,im_index)+LAND_LIST(I)-1) =
     &    melt_before_snow_scheme(i)
        ENDDO
      ENDIF
*D ARE1F404.113,ARE1F404.114  
     &           NTILESDA,TILE_FIELDDA,TILE_PTS,TILE_INDEX,
     &           CAN_EVAP_TILE,MELT_TILE,TILE_FRAC,  
*I ARE1F404.118   
     &       NTILESDA,                                                  
*D ARE1F404.122,ARE1F404.125  
     &       CAN_EVAP_TILE(TILE_FIELDDA,NTILESDA),                      
     &       MELT_TILE(TILE_FIELDDA,NTILESDA),
     &       TILE_FRAC(TILE_FIELDDA,NTILESDA)
*D ARE2F404.413
*I ACB1F304.20    
     &     ,PSLEVEL     !  loop counter for pseudolevels             
     &     ,PSLEVEL_OUT !  index for pseudolevels sent to STASH
*I ABX1F405.937   
     &,     PLLTILE(NTILESDA) ! pseudolevel list for surface tiles
*D ARE2F404.414,ARE2F404.415  
*D ARE2F404.416
     &   D1(JCANOPY_WATER),D1(JRGRAIN_TYP),L_SNOW_ALBEDO,SNODEP_LAND,
*D ARE1F404.142,ARE1F404.146  
     &   NTILES,TILE_PTS,TILE_INDEX,                                    
     &   D1(JCATCH_TYP),CAN_EVAP_TILE,                                  
     &   TILE_FRAC,D1(JINFIL_TYP),MELT_TILE,           
     &   D1(JTSTAR_TYP),D1(JCAN_WATER_TYP),D1(JSNODEP_TYP),
*D ARE2F404.417,ARE2F404.419  
*I AJS1F401.951   

CL ITEM 236 SNOW AMOUNT ON TILES                             
                                                                        
      IF (SF(236,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,236,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),D1(JSNODEP_TYP+((PSLEVEL-1)*LAND_FIELD)),    
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF                                                        
        END DO                                                          
      END IF    

CL ITEM 237 SNOW MELT RATE ON TILES                             
                                                                        
      IF (SF(237,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,237,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),MELT_TILE(1,PSLEVEL),        
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF

CL ITEM 238 SNOW GRAIN SIZE ON TILES                             
                                                                        
      IF (SF(238,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,238,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),D1(JRGRAIN_TYP+((PSLEVEL-1)*LAND_FIELD)),    
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF
                         
*/-----
*DECLARE HYD_IC7A
*/-----
*I HYD_IC7A.53
     &  NSMAX,NELEV,
*D HYD_IC7A.59,HYD_IC7A.62
     &     ROOTD,SATCON,SATHH,SNOW_SUBLIMATION,EI_TILE,
     &     SOILB,SOIL_EVAPORATION,SURF_HT_FLUX,surf_htf_tile,
     &     GBM_OROG,VFRAC,VSAT,VWILT,TIMESTEP,
     &     CAN_WCNT,RGRAIN,L_SNOW_ALBEDO,
     &     SNODEP,STHF,STHU,
*D HYD_IC7A.64
     &     INFIL,T_DEEP_ICE,STF_HF_SNOW_MELT,
*D HYD_IC7A.71
     &     NTILES,TILE_PTS,TILE_INDEX,                                  
*D HYD_IC7A.73,HYD_IC7A.74   
     &     FRAC,INFIL_TILE,MELT_TILE,TSTAR_TILE,       
     &     CAN_WCNT_TILE,SNOW_TILE,
*D HYD_IC7A.75
! New variables/arrays passed into/out of routine, for Essery snow
! scheme:
     &     SNOWDEPTH,RHO_SNOW_GRND,SNOW_CAN,SNOW_SOIL_HTF,
     &     NSNOW,DS,SICE,SLIQ,TSNOWLAYER,RHO_SNOW,RGRAINL,
     &     L_ESSERY_SNOW,dzsnow, refreezing_layer, 
     &     refreezing_tile, refreezing,
     &     melt_tiles_before_snow_scheme,
     &     melt_before_snow_scheme,
     &     surf_ht_flux_ice,tfall_tile,
*D HYD_IC7A.86,HYD_IC7A.87   
*D HYD_IC7A.102,HYD_IC7A.103  
     &   NTILES,                  ! IN Number of tiles. 
     &   TILE_PTS(NTILES),        ! IN Number of tile points.           
     &   TILE_INDEX(LAND,NTILES)  ! IN Index of tile points.            
*D HYD_IC7A.109,HYD_IC7A.113  
     &  CAN_CPY_TILE(LAND,NTILES),! IN Canopy/surface capacity of       
!                                 !    land tiles (kg/m2)     
     &  CAN_EVAP_TILE(LAND,NTILES),!IN Canopy evaporation from
*D HYD_IC7A.119
     &  FRAC(LAND,NTILES),        ! IN Tile fractions                   
*I HYD_IC7A.121   
     &  INFIL_TILE(LAND,NTILES),  ! IN Maximum surface infiltration     
!                                 !    rate for tiles (kg/m2/s)     
*I HYD_IC7A.123   
     &  MELT_TILE(LAND,NTILES),   ! IN Surface snowmelt on tiles
!                                 !    (kg/m2/s)       
*D HYD_IC7A.127,HYD_IC7A.133  
     &  SURF_HT_FLUX(LAND),       ! INOUT Net downward surface heat flux
     &  SURF_HT_FLUX_ICE(LAND),   ! OUT Net downward surface heat flux
!                                 !    (W/m2)
     &  SURF_HTF_TILE(LAND,NTILES),
     &  TSTAR_TILE(LAND,NTILES),  ! IN Tile surface temperatures (K)    
*D HYD_IC7A.152,HYD_IC7A.156  
     &  CAN_WCNT_TILE(LAND,NTILES),!INOUT Canopy water content for      
!                                 !       land tiles (Kg/m2)  
     &  RGRAIN(LAND,NTILES),      ! INOUT Snow grain size (microns) 
     &  SNOW_TILE(LAND,NTILES),   ! INOUT Snowmass on tiles (kg/m2).
*D HYD_IC7A.163,HYD_IC7A.164  
     &  T_DEEP_SOIL(LAND,ST_LEVELS)!INOUT Deep soil temps. (K) 
     &  ,T_DEEP_ICE(LAND,ST_LEVELS)
*D HYD_IC7A.176
     &  SNODEP(LAND),             ! OUT Snow depth (Kg of water)        
*I HYD_IC7A.207   
     &  SNOW_SUBLIMATION(LAND),
*D HYD_IC7A.210
*I HYD_IC7A.213
      INTEGER NSMAX,NELEV
      LOGICAL L_ESSERY_SNOW
      REAL dzsnow(nsmax)
!
      REAL refreezing_layer(LAND,NTILES,NSMAX)
      REAL refreezing_tile(LAND,NTILES)
      REAL tfall_tile(LAND,NTILES)
      REAL refreezing(LAND)

      REAL melt_tiles_before_snow_scheme(land, ntiles)
      REAL melt_before_snow_scheme(land)
! Additional variables for JULES
      REAL
     & SNOWDEPTH(    land,ntiles),
!                            ! Snow depth on ground on tiles (m)
     & RHO_SNOW_GRND(land,ntiles),
!                            ! Snowpack bulk density (kg/m3)
     & SNOW_CAN(     land,ntiles),
!                            ! Snow on the canopy (kg/m2)
     & SNOW_SOIL_HTF(land,ntiles),
!                            ! Surface heat flux beneath snow (W/m2)
     & NSNOW(        land,ntiles),
!                            ! Number of snow layers on ground on tiles
!                            ! NOTE that this is converted to an integer.
     & DS(        land,ntiles,nsmax),
!                            ! Snow layer thickness (m)
     & SICE(      land,ntiles,nsmax),
!                            ! Snow layer ice mass on tiles (Kg/m2)
     & SLIQ(      land,ntiles,nsmax),
!                            ! Snow layer liquid mass on tiles (Kg/m2)
     & TSNOWLAYER(land,ntiles,nsmax),
!                            ! Snow layer temperature (K)
     & RHO_SNOW(  land,ntiles,nsmax),
!                            ! Snow layer densities (kg/m3)
     & RGRAINL(   land,ntiles,nsmax),
     & EI_TILE(   land,ntiles),
     & gbm_orog(land)
!
*D HYD_IC7A.229
     &     GBM_OROG,LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,NTILES,              
*D HYD_IC7A.230
     &     TILE_PTS,TILE_INDEX,LAND,SM_LEVELS,NSMAX,NELEV,
*D HYD_IC7A.233,HYD_IC7A.241  
     &     EXT,HCAP,HCON,INFIL_TILE,LS_RAIN,LS_SNOW,MELT_TILE,
     &     EI_TILE,SATCON,SATHH,
     &     SURF_HT_FLUX,SURF_HTF_TILE,TSTAR_TILE,                      
     &     FRAC,TIMESTEP,VSAT,VWILT,L_SNOW_ALBEDO,
     &     STF_HF_SNOW_MELT,STF_SUB_SURF_ROFF,                        
     &     CAN_WCNT_TILE,RGRAIN,SMCL,                  
     &     SNOW_TILE,STHF,STHU,T_DEEP_SOIL,  
     &     CAN_WCNT_GB,HF_SNOW_MELT,INFIL,SNODEP,SMC,SNOW_MELT, 
*B HYD_IC7A.242
     &     T_DEEP_ICE,
*D HYD_IC7A.243
                                    ! Snow layer grain size on tiles (microns)
*I HYD_IC7A.245
! New variables/arrays passed into/out of routine, for Essery snow
! scheme:
     &     ,SNOWDEPTH, RHO_SNOW_GRND, SNOW_CAN,
     &     SNOW_SOIL_HTF, NSNOW, DS, SICE, SLIQ, TSNOWLAYER,
     &     RHO_SNOW, RGRAINL,L_ESSERY_SNOW,
     & dzsnow, refreezing_layer, refreezing_tile, refreezing
     &, melt_tiles_before_snow_scheme
     &, melt_before_snow_scheme
     &, surf_ht_flux_ice,tfall_tile
*/-----
*DECLARE ICEHTC5A
*/-----
*I ICEHTC5A.26
     &,L_ESSERY_SNOW
*I ICEHTC5A.55
      REAL RHO_ICE,ICE_HCON,ICE_HCAP
      PARAMETER (
     + RHO_ICE=917.0
     +,ICE_HCON=2.4
     +,ICE_HCAP=1.80E6
     +)
*I ICEHTC5A.76
     &,L_ESSERY_SNOW
*I ICEHTC5A.100
          if (L_ESSERY_SNOW) THEN 
            H_FLUX(I,N)=-ICE_HCON*2.0*(TSOIL(I,N+1)-TSOIL(I,N))
     &                             /(DZ(N+1)+DZ(N))
          !!!!GLIMMER COUPLING - USE ONLY TOP LAYER AS A THERMAL BUFFER, ONLY
          !!!!INFLUENCE IS SURF_HT_FLUX(I)
          !!!!            H_FLUX(I,N)=0.
          !!!!
          else
*I ICEHTC5A.102
          end if
*I ICEHTC5A.120
          if (L_ESSERY_SNOW) THEN 
          TSOIL(I,N)=TSOIL(I,N)
     &     +1.0/(ICE_HCAP*DZ(N))*(H_FLUX(I,N-1)
     &     -H_FLUX(I,N))*TIMESTEP
          else
*I ICEHTC5A.124
          endif

*/-----
*DECLARE INITMIN1
*/-----
*D INITMIN1.109
        DO N=(SOIL_1-1)*NELEV+1,SOIL_1*NELEV
*I INITMIN1.111
        ENDDO
*D INITMIN1.143
        DO N=(SOIL_1-1)*NELEV+1,SOIL_1*NELEV
*I INITMIN1.145
        ENDDO
*/-----
*DECLARE IMCLTQ7A
*/-----
*D IMCLTQ7A.47
     &                      LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS,      
*D IMCLTQ7A.49
     &                      ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,        
*D IMCLTQ7A.53
     &                      SNOW_TILE,TILE_FRAC,                        
*D IMCLTQ7A.70,IMCLTQ7A.72   
     &,NTILES                      ! IN Number of land surface tiles.   
     &,TILE_PTS(NTILES)            ! IN Number of tiles.                
     &,TILE_INDEX(LAND_FIELD,NTILES)!IN Index of tile points.           
*D IMCLTQ7A.75
     & ALPHA1(LAND_FIELD,NTILES)   ! IN Gradient of saturated specific  
*D IMCLTQ7A.81,IMCLTQ7A.84   
!                                  !    heat flux into sea-ice (W/m2/K).
     &,ASHTF_TILE(LAND_FIELD,NTILES)!IN Coefficient to calculate heat
!                                  !    flux into land tiles (W/m2/K).  
*D IMCLTQ7A.91
     &,RESFT(LAND_FIELD,NTILES)    ! IN Total resistance factor         
*D IMCLTQ7A.94
     &,RHOKH_1(LAND_FIELD,NTILES)  ! IN Land surface exchange coeff.    
*D IMCLTQ7A.97
     &,RHOKPM(LAND_FIELD,NTILES)   ! IN Land surface exchange coeff.    
*D IMCLTQ7A.99,IMCLTQ7A.100  
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow on land tiles (kg/m2)
     &,TILE_FRAC(LAND_FIELD,NTILES)! IN Tile fractions.                 
*D IMCLTQ7A.111
     &,FQW_TILE(LAND_FIELD,NTILES) ! INOUT Surface flux of QW for land  
*D IMCLTQ7A.117
     &,FTL_TILE(LAND_FIELD,NTILES) ! INOUT H/Cp for land tiles.         
*I IMCLTQ7A.157   
     &,LAT_HT   ! Latent heat of evaporation for snow-free land 
!               ! or sublimation for snow-covered land and ice.
*D IMCLTQ7A.224,IMCLTQ7A.225  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*D IMCLTQ7A.232
     &                      RESFT(L,N)*RHOKPM(L,N) *
     &                      (CP*RHOKH_1(L,N) + ASHTF_TILE(L,N))
*D IMCLTQ7A.234,IMCLTQ7A.244  
*D IMCLTQ7A.274,IMCLTQ7A.275  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*D IMCLTQ7A.279
          LAT_HT = LC
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS   
          CT1(I) = CT1(I) - LAT_HT*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) * 
*D IMCLTQ7A.282,IMCLTQ7A.283  
     &          RHOKPM(L,N)*( LAT_HT*ALPHA1(L,N)*RESFT(L,N)*RHOKH_1(L,N)
     &                                               + ASHTF_TILE(L,N) )
*D IMCLTQ7A.285,IMCLTQ7A.296  
*D IMCLTQ7A.400,IMCLTQ7A.401  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*I IMCLTQ7A.404   
          LAT_HT = LC
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS
*D IMCLTQ7A.406,IMCLTQ7A.407  
     &                             ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N) *
     &                               (DQW(I,1) - ALPHA1(L,N)*DTL(I,1)) 
     &                                      - ASHTF_TILE(L,N)*DTL(I,1) )
*D IMCLTQ7A.409,IMCLTQ7A.410  
     &                    RHOKPM(L,N)*( CP*RHOKH_1(L,N) * 
     &                               (DQW(I,1) - ALPHA1(L,N)*DTL(I,1)) 
     &                                      + ASHTF_TILE(L,N)*DQW(I,1) )
*D IMCLTQ7A.414,IMCLTQ7A.428  
*/-----
*DECLARE INITIAL1
*/-----
*D ABX1F404.224
      IF (L_VEG_FRACS .AND. STEPim(a_im) == 0 ) THEN
*/-----
*DECLARE INITVEG1
*/-----
*I INITVEG1.47    
!   4.6   16/02/99   Correct initialisation of accumulated leaf turnover
!                    rate when both TRIFFID and phenology in use.
!                    Richard Betts
*D INITVEG1.88,INITVEG1.90   
*D INITVEG1.125,INITVEG1.136  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX, 
     &            D1(JFRAC_TYP),D1(JCANHT_PFT),           
     &            D1(JLAI_PFT),D1(JSAT_SOIL_COND),
     &            D1(JCATCH_TYP),D1(JINFIL_TYP),D1(JZ0_TYP))
*D INITVEG1.187
      ENDIF
*I INITVEG1.188   
      IF (L_PHENOL) THEN
*D ABX1F405.72,ABX1F405.75   
        DO L = 1,LAND_FIELD
          D1(JG_LF_PFT_ACC+L-1) = 0.0
        ENDDO
*/-----
*DECLARE LEAF7A
*/-----
*D LEAF7A.174
        RD(L) = FDC3 * VCM(L)   
*D LEAF7A.189
        RD(L) = FDC3 * VCM(L)  
*D LEAF7A.465
        RD(L) = FDC4 * VCM(L)  
*D LEAF7A.481
        RD(L) = FDC4 * VCM(L)   
*/-----
*DECLARE LOTKA2A
*/-----
*D ABX1F405.1644
     &,                 C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI,PC_S    
*I ABX1F405.1655  
     &,FRAC_AGRIC(LAND_FIELD)     ! IN Fraction of agriculture.    
*D LOTKA2A.53,ABX1F405.1651 
*I LOTKA2A.153   
C----------------------------------------------------------------------
C Exclude non-crop types from agricultural regions
C----------------------------------------------------------------------
*D LOTKA2A.157
          SPACE(L,N)=1.0-NOSOIL(L)-FRAC_AGRIC(L)*(1-CROP(N))
     &                            -FRAC_MIN*(NPFT-K)  
*D LOTKA2A.176,LOTKA2A.177  
          B(L,N) = PC_S(L,N)*SPACE(L,N)/C_VEG(L,N)-G_AREA(N)   
*/-----
*DECLARE LWRAD3A
*/-----
*I ADB2F404.631   
!       4.6             10-05-98                Land flag passed to     
!                                               FLUX_CALC.              
!                                               (J. M. Edwards)         
*I ASK1F405.289   
!                                                                       
*D LWRAD3A.27
     &   , TAC, PEXNER, TSTAR, TSTAR_SOLID, TSTAR_SEA, L_CTILE 
     &   , PSTAR, AB, BB, AC, BC                    
*D LWRAD3A.33
     &   , LAND, FLANDG, ICE_FRACTION                                   
*I LWRAD3A.58
     &   ,bl_levels,downwelling_lw
*I LWRAD3A.219   
     &   , L_CTILE                                                      
!             COASTAL TILING SWITCH                                     

      INTEGER BL_LEVELS
      REAL downwelling_lw(NPD_FIELD,BL_LEVELS)

      REAL      !, INTENT(IN)                                           
     &     FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
*I LWRAD3A.224   
     &   , TSTAR_SOLID(NPD_FIELD)                                       
!             SOLID SURFACE TEMPERATURE                                 
     &   , TSTAR_SEA(NPD_FIELD)                                         
!             OPEN SEA SURFACE TEMPERATURES                             
*D LWRAD3A.226
!             SEA ICE FRACTION OF SEA PORTION OF GRID BOX               
*I LWRAD3A.321   
     &   , T_SOLID(NPD_PROFILE)                                         
!             GATHERED TEMPERATURE OF SOLID SURFACE                     
     &   , T_SEA(NPD_PROFILE)                                           
!             GATHERED OPEN SEA TEMPERATURE                             
*D ADB1F401.537,ADB1F401.540  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.556,ADB1F401.557  
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA 
     &   , AB, BB, AC, BC, PEXNER, TAC                    
     &   , P, T, T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS               
*D ADB1F401.568
     &   , FLANDG                                                     
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION           
*D ADB1F401.574
     &   , P, T, T_SURFACE, T_SOLID, T_SEA, T_BDY, D_MASS               
*D ADB1F401.575,LWRAD3A.629  
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION           
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG, LWSEA               
*D LWRAD3A.631
     &   , L_DUMMY, DUMMY, DUMMY, DUMMY                           
     &   , L_DUMMY, L_CTILE, DUMMY, DUMMY, DUMMY, DUMMY         
*I LWRAD3A.665
C     Store flux for bottom of level i, so i=nlevs+1 is TOA
C     flux_net/up level 0 is the TOA, nlevs the surface
C     for subgridscale orog, LW_DOWN=flux_net(nlev-1)+flux_up(nlev-1)
      do i=1,BL_LEVELS
        do l=1,n_profile
          downwelling_lw(l,i)=flux_net(l,nlevs-i+1)
     &    +flux_up(l,nlevs-i+1)
        enddo
      enddo
*D LWRAD3A.706
         IF (FLANDG(L).EQ.1.0.OR.ICE_FRACTION(L).EQ.1.0) THEN           
*D LWRAD3A.708
         ELSE IF (FLANDG(L).LT.TOL_TEST.AND.
     &     ICE_FRACTION(L).LT.TOL_TEST) THEN                     
*D ADB1F401.576,ADB1F401.577  
!           LWSEA MUST BE SCALED BY THE FRACTION OF OPEN SEA TO         
!           TOTAL SEA FOR CONSISTENCY WITH UPPER LEVELS IN THE MODEL.   
*D LWRAD3A.714
            LWOUT(L, 1)=LWOUT(L, 1)-(1.0E+00-FLANDG(L))*LWSEA(L)        
*D ADB1F401.579
     &   , FLANDG                                                       
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
*I ADB1F401.584   
!              FRACTION OF SEA-ICE IN SEA PART OF GRID BOX.             
     &   , FLANDG(NPD_PROFILE)                                          
!                  LAND FRACTION                                        
                                                                        
*D ADB1F401.587,ADB1F401.590  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.598,ADB1F401.599  
!     SET THE FRACTIONAL OPEN SEA COVERAGE. POINTS ARE REQUIRED WHERE   
!     THIS IS NEITHER 0 NOR 1.                                          
*D ADB1F401.601
       SEARCH_ARRAY(L)=(1.0E+00-FLANDG(L))*(1.0E+00-ICE_FRACTION(L))    
       SEARCH_ARRAY(L)=SEARCH_ARRAY(L)*(1.-SEARCH_ARRAY(L))             
*D GSS2F402.246
      N_FRAC_SOL_POINT=0                                                
*D GSS2F402.249,GSS2F402.250  
          N_FRAC_SOL_POINT                  =N_FRAC_SOL_POINT+1         
          I_FRAC_SOL_POINT(N_FRAC_SOL_POINT)=L                          
*/-----
*DECLARE MICROB7A
*I MICROB7A.24
     &,                   FRAC 
*I MICROB7A.37
     &,FRAC(LAND_FIELD) 
*I MICROB7A.56
*/NSTYPES defines the number of Plant Functional Types (PFTs)
*/as well as various other surface-specific quantities
*CALL NSTYPES
C
*/This calls the previously defined common block
*CALL C_LAND_CC
C
*/The following 6 lines Remove Q10 and KAPS
*/from the microbiology deck. In Ben Booth's original mod
*/it is only Q10 which is deleted in this way. I can't
*/work out why this is but this solution is necessary
*/due to the inclusion of the definition of KAPS
*/within the common deck
*D MICROB7A.63,MICROB7A.68
*D MICROB7A.73
         IF (FRAC(L) .LT. 1.0) THEN
*/-----
*/-----
*DECLARE NAMSIZE
*/-----
*D WRB1F401.628
!
     & W_LEN_EXTCNST, NSMAX, GLOBAL_WAVE,
!
*/-----

*I RB300993.113   
     & ,NTILES
*/-----
*DECLARE NSTYPES
*/-----
*I ABX1F405.1724
     +,nnvg_1
     +,npft_1
     +,ntype_1
     +,soil_1
     +,lake_1
     +,NELEV
*D ABX1F405.1725
!+seg Elevations changed from 25 to 10, Mar 2016
      PARAMETER (NELEV=10)                       ! number of subgrid elevations
      PARAMETER (NNVG_1=4,NPFT_1=5,NTYPE_1=9,SOIL_1=8,lake_1=7)
      PARAMETER (NNVG=nnvg_1*NELEV 
     &          ,NPFT=npft_1*NELEV 
     &          ,NTYPE=ntype_1*NELEV 
     &          )
*/-----
*DECLARE NVEGPARM
*/-----
*D NVEGPARM.4,NVEGPARM.10   
     + ALBSNC_NVG(NNVG_1)           ! Snow-covered albedo.                   NVEGPARM.4      
     +,ALBSNF_NVG(NNVG_1)           ! Snow-free albedo.                      NVEGPARM.5      
     +,CATCH_NVG(NNVG_1)            ! Canopy capacity (kg/m2).               NVEGPARM.6      
     +,GS_NVG(NNVG_1)               ! Surface conductance (m/s).             NVEGPARM.7      
     +,INFIL_NVG(NNVG_1)            ! Infiltration enhancement factor.       NVEGPARM.8      
     +,ROOTD_NVG(NNVG_1)            ! Rootdepth (m).                         NVEGPARM.9      
     +,Z0_NVG(NNVG_1)               ! Roughness length (m).                  NVEGPARM.10     

*D NVEGPARM.16,NVEGPARM.18   
      DATA CATCH_NVG   /  0.50,  0.00,  0.00,  0.00 /    
      DATA GS_NVG      /  0.00,  0.00,  1E-2,   1E6 /
      DATA INFIL_NVG   /  0.10,  0.00,  0.50,  0.00 /
*/-----
*DECLARE PFTPARM
*/-----
*D PFTPARM.7,PFTPARM.16
     + ALBSNC_MAX(NPFT_1)           ! Snow-covered albedo for large LAI.
     +,ALBSNC_MIN(NPFT_1)           ! Snow-covered albedo for zero LAI.   
     +,ALBSNF_MAX(NPFT_1)           ! Snow-free albedo for large LAI.
     +,DZ0V_DH(NPFT_1)              ! Rate of change of vegetation
C                                 ! roughness length with height.
     +,CATCH0(NPFT_1)               ! Minimum canopy capacity (kg/m2).
     +,DCATCH_DLAI(NPFT_1)          ! Rate of change of canopy capacity
C                                 ! with LAI.
     +,INFIL_F(NPFT_1)              ! Infiltration enhancement factor.
     +,KEXT(NPFT_1)                 ! Light extinction coefficient.
     +,ROOTD_FT(NPFT_1)             ! Rootdepth (m).
*/-----
*DECLARE PHIMH6A
*/-----
*D PHIMH6A.2
*IF DEF,A03_8A                                                          
*D PHIMH6A.21
!!!   SUBROUTINES PHI_M_H_SEA AND PHI_M_H_LAND ------------------------
*D PHIMH6A.39,PHIMH6A.40   
      SUBROUTINE PHI_M_H_SEA(
     & P_POINTS,P_FIELD,P1,LAND_MASK,
*D PHIMH6A.52
*D PHIMH6A.84
      INTEGER I       ! Loop counter; horizontal field index.           
*D PHIMH6A.104
      DO I=P1,P1+P_POINTS-1                                             
*D PHIMH6A.106,PHIMH6A.111  
        IF ( .NOT. LAND_MASK(I) ) THEN
*D PHIMH6A.113,PHIMH6A.115  
*D PHIMH6A.162
        ENDIF  ! LAND_MASK
                                                                        
      ENDDO                                                             
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('PHI_M_H ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

!!!                                                                     
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE PHI_M_H_LAND(
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & RECIP_L_MO,Z_UV,Z_TQ,Z0M,Z0H,PHI_M,PHI_H,LTIMER                  
     &)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      LOGICAL                                                           
     & LTIMER                                                           
                                                                        
                                                                        
      REAL                                                              
     & RECIP_L_MO(LAND_FIELD)
!                       ! IN Reciprocal of the Monin-Obukhov length 
!                       !     (m^-1).
     &,Z_UV(P_FIELD)    ! IN Height of wind level above roughness
!                       !    height (m)
     &,Z_TQ(P_FIELD)    ! IN Height of temperature, moisture and scalar
!                       !    lev above the roughness height (m).        
     &,Z0M(LAND_FIELD)  ! IN Roughness length for momentum (m).
     &,Z0H(LAND_FIELD)  ! IN Roughness length for heat/moisture/scalars
!                       !    (m)
!                                                                       
      REAL                                                              
     & PHI_M(LAND_FIELD)! OUT Stability function for momentum.
     &,PHI_H(LAND_FIELD)! OUT Stability function for
!                       !     heat/moisture/scalars.
!                                                                       
!*L  Workspace usage----------------------------------------------------
!    No work areas are required.                                        
!                                                                       
!*----------------------------------------------------------------------
!*L  External subprograms called:                                       
                                                                        
      EXTERNAL TIMER                                                    
                                                                        
!*----------------------------------------------------------------------
!  Common and local physical constants.                                 
!                                                                       
!  None.                                                                
!                                                                       
!  Define local variables.                                              
!                                                                       
      INTEGER I,J,L     ! Loop counter; horizontal field index.
!                                                                       
      REAL                                                              
     & PHI_MN         ! Neutral value of stability function for momentum
     &,PHI_HN         ! Neutral value of stability function for scalars.
     &,ZETA_UV        ! Temporary in calculation of PHI_M.              
     &,ZETA_0M        ! Temporary in calculation of PHI_M.              
     &,ZETA_TQ        ! Temporary in calculation of PHI_H.              
     &,ZETA_0H        ! Temporary in calculation of PHI_H.              
     &,X_UV_SQ        ! Temporary in calculation of PHI_M.              
     &,X_0M_SQ        ! Temporary in calculation of PHI_M.              
     &,X_UV           ! Temporary in calculation of PHI_M.              
     &,X_0M           ! Temporary in calculation of PHI_M.              
     &,Y_TQ           ! Temporary in calculation of PHI_H.              
     &,Y_0H           ! Temporary in calculation of PHI_H.              
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('PHI_M_H ',3)                                        
      ENDIF                                                             
                                                                        
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
!                                                                       
!-----------------------------------------------------------------------
!! 1. Calculate neutral values of PHI_M and PHI_H.                      
!-----------------------------------------------------------------------
!                                                                       
        PHI_MN = LOG( (Z_UV(I) + Z0M(L)) / Z0M(L) )
        PHI_HN = LOG( (Z_TQ(I) + Z0M(L)) / Z0H(L) )
!                                                                       
!-----------------------------------------------------------------------
!! 2. Calculate stability parameters.                                   
!-----------------------------------------------------------------------
!                                                                       
        ZETA_UV = (Z_UV(I) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_TQ = (Z_TQ(I) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_0M = Z0M(L) * RECIP_L_MO(L)
        ZETA_0H = Z0H(L) * RECIP_L_MO(L)
!                                                                       
!-----------------------------------------------------------------------
!! 3. Calculate PHI_M and PHI_H for neutral and stable conditions.      
!-----------------------------------------------------------------------
!                                                                       
        IF (RECIP_L_MO(L) .GE. 0.0) THEN
          PHI_M(L) = PHI_MN + 4.0 * (ZETA_UV - ZETA_0M)
          PHI_H(L) = PHI_HN +
     &               (1.0 + 2.0*ZETA_TQ) * (1.0 + 2.0*ZETA_TQ) -
     &               (1.0 + 2.0*ZETA_0H) * (1.0 + 2.0*ZETA_0H)
!                                                                       
!-----------------------------------------------------------------------
!! 4. Calculate PHI_M and PHI_H for unstable conditions.                
!-----------------------------------------------------------------------
!                                                                       
        ELSE
                                                                        
          X_UV_SQ = SQRT(1.0 - 16.0*ZETA_UV)
          X_0M_SQ = SQRT(1.0 - 16.0*ZETA_0M)
          X_UV = SQRT(X_UV_SQ)
          X_0M = SQRT(X_0M_SQ)
          PHI_M(L) = PHI_MN - 2.0*LOG( (1.0+X_UV) / (1.0+X_0M) )
     &                    - LOG( (1.0+X_UV_SQ) / (1.0+X_0M_SQ) )
     &                    + 2.0*( ATAN(X_UV) - ATAN(X_0M) )
                                                                        
          Y_TQ = SQRT(1.0 - 16.0*ZETA_TQ)
          Y_0H = SQRT(1.0 - 16.0*ZETA_0H)
          PHI_H(L) = PHI_HN - 2.0*LOG( (1.0+Y_TQ) / (1.0+Y_0H) )
                                                                        
        ENDIF
*/-----
*DECLARE PHYSIO7A
*/-----
*D PHYSIO7A.33
     &,                   P_FIELD,NSHYD,NTILES,TILE_PTS,TILE_INDEX      
*D PHYSIO7A.36,ABX1F405.839  
     &,                   V_CRIT_NEW,V_SAT,V_WILT,WIND,Z0_TILE,Z1
     &,                   CANHC_TILE,VFRAC_TILE,FLAKE,G_LEAF,GS,GS_TILE
     &,                   GPP,GPP_FT,NPP,NPP_FT,RESP_P,RESP_P_FT       
     &,                   RESP_S,RESP_W_FT,SMCT,WT_EXT_TILE)          
*I PHYSIO7A.56    
     &,NTILES                     ! IN Number of surface tiles.
*D PHYSIO7A.67
     &,FRAC(LAND_FIELD,NTYPE)     ! IN Surface type fractions.          
*D PHYSIO7A.78
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D PHYSIO7A.80
     &,V_CRIT_NEW(LAND_FIELD)
*D PHYSIO7A.91
     &,Z0_TILE(LAND_FIELD,NTILES) ! IN Tile roughness lengths (m).
*D ABX1F405.840,PHYSIO7A.98   
     & CANHC_TILE(LAND_FIELD,NTILES)
!                                 ! OUT Areal heat capacity of canopy
!                                 !     for land tiles (J/K/m2).
     &,FLAKE(LAND_FIELD,NTILES)   ! OUT Lake fraction.
     &,G_LEAF(LAND_FIELD,NPFT)    ! OUT Leaf turnover rate (/360days).  
     &,GS_TILE(LAND_FIELD,NTILES) ! OUT Surface conductance for
*D PHYSIO7A.113
     &,VFRAC_TILE(LAND_FIELD,NTILES)
!                                 ! OUT Fractional canopy coverage for
!                                 !     land tiles.
     &,WT_EXT_TILE(LAND_FIELD,NSHYD,NTILES)
!                                 ! OUT Fraction of evapotranspiration  
*D PHYSIO7A.115
!                                 !     soil layer by each tile.
*D PHYSIO7A.118
     & CANHC(LAND_FIELD)          ! WORK Canopy heat capacity (J/K/m2).
     &,CH_TYPE(LAND_FIELD,NTYPE)  ! WORK CANHC for surface types. 
     &,F_ROOT(NSHYD)              ! WORK Fraction of roots in each soil 
*I PHYSIO7A.120   
     &,GSOIL(LAND_FIELD)          ! WORK Bare soil conductance.
     &,GS_TYPE(LAND_FIELD,NTYPE)  ! WORK Conductance for surface types.
*I PHYSIO7A.126   
     &,TSTAR(LAND_FIELD)          ! WORK Surface temperature (K).
     &,VFRAC(LAND_FIELD)          ! WORK Fractional canopy coverage. 
     &,VF_TYPE(LAND_FIELD,NTYPE)  ! WORK VFRAC for surface types. 
     &,WT_EXT(LAND_FIELD,NSHYD)   ! WORK Gridbox-mean WT_EXT.  
     &,WT_EXT_TYPE(LAND_FIELD,NSHYD,NTYPE) 
!                                 ! WORK WT_EXT for surface types. 
     &,Z0(LAND_FIELD)             ! WORK Roughness length (m).
     &,FRAC_ICE(LAND_FIELD)       ! WORK GBM ice-surface fraction
     &,elev_frac(LAND_FIELD,NELEV)! WORK ELEV non-ice-surface fraction
*D PHYSIO7A.129
     & I,J,K,L,N,N_1,KK                 ! Loop indices
*I PHYSIO7A.135
*CALL C_LAND_CC
*I PHYSIO7A.142
*/
*/Define the new critical soil moisture concentration
      DO L=1,LAND_FIELD
         V_CRIT_NEW(L) = V_WILT(L) +
     &      V_CRIT_ALPHA * (V_SAT(L) - V_WILT(L))
      ENDDO
*I PHYSIO7A.149   
          DO N=1,NTYPE 
            WT_EXT_TYPE(L,K,N)=0.0
          ENDDO                        
*D PHYSIO7A.174
          GS_TYPE(L,N)=GS(L)
*I PHYSIO7A.185   
        CANHC(L)=0.0
        VFRAC(L)=0.0
      ENDDO

      DO L=LAND1,LAND1+LAND_PTS-1 
        IF (V_CRIT_NEW(L).GT.0.)
     &   GSOIL(L)=GS_NVG(SOIL_1-NPFT_1)*
     *            (STHU(L,1)*V_SAT(L)/V_CRIT_NEW(L))**2
*I PHYSIO7A.193
        n_1=(n-1)/nelev + 1

        IF (NTILES.EQ.1) THEN                                           
          DO L=1,LAND_FIELD
            TSTAR(L) = TSTAR_TILE(L,1)
            Z0(L) = Z0_TILE(L,1)
          ENDDO 
        ELSE IF (NTILES.eq.NELEV*2) THEN 
          k=mod(n-1,nelev)+1
          DO L=1,LAND_FIELD
            TSTAR(L) = TSTAR_TILE(L,K)
            Z0(L) = Z0_TILE(L,K)
          ENDDO
        ELSE
          DO L=1,LAND_FIELD
            TSTAR(L) = TSTAR_TILE(L,N)
            Z0(L) = Z0_TILE(L,N)
          ENDDO
        ENDIF                                                      
*D PHYSIO7A.195
        CALL ROOT_FRAC(NSHYD,DZSOIL,ROOTD_FT(N_1),F_ROOT)

*D PHYSIO7A.198,PHYSIO7A.199  
     &,               F_ROOT,STHU,V_CRIT_NEW,V_SAT,V_WILT 
     &,               WT_EXT_TYPE(1,1,N),FSMC)
*D PHYSIO7A.203
     &,             RIB,WIND,Z0,Z0,Z1,RA)                   
*D PHYSIO7A.206
     &,               TILE_PTS(N),TILE_INDEX(1,N),N_1
*D PHYSIO7A.208
     &,               Q1,RA,TSTAR
*D PHYSIO7A.210
     &,               RESP_W_FT(1,N),GS_TYPE(1,N))

        CALL SOIL_EVAP (LAND_FIELD,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)
     &,                 GSOIL,LAI(1,N),GS_TYPE(1,N),WT_EXT_TYPE(1,1,N))
*D PHYSIO7A.213
     &,                N_1,FSMC,TSTAR,G_LEAF(1,N)) 

        CALL CANCAP (LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N),N_1
     &,              HT(1,N),LAI(1,N),CH_TYPE(1,N),VF_TYPE(1,N))
*D PHYSIO7A.218,PHYSIO7A.219  
! Non-vegetated surface types 
*D PHYSIO7A.221,PHYSIO7A.228  
      DO N=NPFT+1,NTYPE                                             
        n_1=(n-npft-1)/nelev + 1
*D PHYSIO7A.231
          GS_TYPE(L,N) = GS_NVG(N_1)
          DO K=1,NSHYD 
            WT_EXT_TYPE(L,K,N) = 0.                   
          ENDDO                                        
*I PHYSIO7A.232   
      ENDDO 
*I PHYSIO7A.233   
! Copy soil conductance and add bare soil fraction to extraction from
! surface layer
      DO N=(soil_1-1)*nelev+1,(soil_1)*nelev
      DO J=1,TILE_PTS(N)                                              
        L=TILE_INDEX(J,N)                                             
        GS_TYPE(L,N) = GSOIL(L)
        WT_EXT_TYPE(L,1,N) = 1.                          
      ENDDO  
      ENDDO  

!----------------------------------------------------------------------
! Canopy heat capacity and coverage set to 0 for non-vegetated surfaces 
!---------------------------------------------------------------------- 
      DO N=NPFT+1,NTYPE                                               
        DO J=1,TILE_PTS(N)                                              
          L=TILE_INDEX(J,N)
          CH_TYPE(L,N) = 0.
          VF_TYPE(L,N) = 0.
        ENDDO
      ENDDO

      IF (NTILES.EQ.NELEV*2) THEN
      DO L=LAND1,LAND1+LAND_PTS-1
      frac_ice(l)=0.
      DO K=1,NELEV
        elev_frac(l,k)=0.
      ENDDO
      do n=(ntype-nelev)+1,ntype
        frac_ice(l)=frac_ice(l)+frac(l,n)
      enddo
      DO N=1,NTYPE-NELEV
        k=mod(n-1,nelev)+1
        elev_frac(l,k)=elev_frac(l,k)+frac(l,n)
        gs_tile(l,k)=0.
        canhc_tile(l,k)=0.
        vfrac_tile(l,k)=0.
        do kk=1,nshyd
          wt_ext_tile(l,kk,k)=0.
        end do
      END DO
      END DO
      ENDIF
*D PHYSIO7A.234
*I PHYSIO7A.239
     &,             FRAC_ICE
*D PHYSIO7A.245

      IF (NTILES.EQ.1) THEN
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)                                           
            CANHC(L) = CANHC(L) + FRAC(L,N)*CH_TYPE(L,N) 
            GS(L) = GS(L) + FRAC(L,N)*GS_TYPE(L,N)                     
            VFRAC(L) = VFRAC(L) + FRAC(L,N)*VF_TYPE(L,N)
            DO K=1,NSHYD               
              WT_EXT(L,K) = WT_EXT(L,K) + FRAC(L,N)*WT_EXT_TYPE(L,K,N)
            ENDDO                                  
          ENDDO                                                         
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          flake(l,1)=0.
          do n=(lake_1-1)*nelev+1,lake_1*nelev
            flake(l,1)=flake(l,1)+frac(l,n)
          enddo
          GS_TILE(L,1) = 0.
          IF (FLAKE(L,1).LT.1.)
     &      GS_TILE(L,1) = GS(L) / (1. - FLAKE(L,1))
          CANHC_TILE(L,1) = CANHC(L)                                    
          VFRAC_TILE(L,1) = VFRAC(L)
          DO K=1,NSHYD
            WT_EXT_TILE(L,K,1) = WT_EXT(L,K)
          ENDDO                     
        ENDDO
      ELSEIF (NTILES.EQ.NELEV*2) THEN
        DO N=1,NTYPE-NELEV
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)                                           
            if (frac_ice(l) .LT. 1) THEN
            CANHC(L) = CANHC(L) + FRAC(L,N)*CH_TYPE(L,N)/(1-frac_ice(l)) 
            GS(L) = GS(L) + FRAC(L,N)*GS_TYPE(L,N)/(1-frac_ice(l))                     
            VFRAC(L) = VFRAC(L) + FRAC(L,N)*VF_TYPE(L,N)/(1-frac_ice(l))
            DO K=1,NSHYD               
              WT_EXT(L,K) = WT_EXT(L,K) + FRAC(L,N)
     &                     *WT_EXT_TYPE(L,K,N)/(1-frac_ice(l))
            ENDDO                                  
            endif
          ENDDO                                                         
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          do n=1,ntiles
            k=mod(n-1,nelev)+1
            flake(l,n)=0.
            if (n.le.nelev.AND. elev_frac(l,k).gt.0)
     &        flake(l,n)=frac(l,(lake_1-1)*nelev+n)/elev_frac(l,k)
          end do
        ENDDO
        DO N=1,NTYPE-NELEV
          k=mod(n-1,nelev)+1
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)                                           
            if (elev_frac(l,k).gt.0) then
            IF (FLAKE(L,K).LT.1.) 
     &      GS_TILE(L,K) = GS_TILE(L,K)+
     &        FRAC(L,N)*GS_TYPE(L,N)/elev_frac(l,k)/ (1. - FLAKE(L,K))
            CANHC_TILE(L,K) = CANHC_TILE(L,K)+
     &        FRAC(L,N)*CH_TYPE(L,N)/elev_frac(l,k)
            VFRAC_TILE(L,K) = VFRAC_TILE(L,K)+
     &        FRAC(L,N)*VF_TYPE(L,N)/elev_frac(l,k)
            DO KK=1,NSHYD
             WT_EXT_TILE(L,KK,K) = WT_EXT_TILE(L,KK,K)+
     &          FRAC(L,N)*WT_EXT_TYPE(L,KK,K)/elev_frac(l,k)
            END DO                     
          endif
          END DO                     
        END DO                     
        DO N=NTYPE-NELEV+1,NTYPE
           k=mod(n-1,nelev)+1
           DO J=1,TILE_PTS(N)
           L=TILE_INDEX(J,N)
            GS_TILE(L,K+NELEV) = GS_TYPE(L,N)
            CANHC_TILE(L,K+NELEV) = CH_TYPE(L,N)            
            VFRAC_TILE(L,K+NELEV) = VF_TYPE(L,N)
            DO KK=1,NSHYD
              WT_EXT_TILE(L,KK,K+NELEV) = WT_EXT_TYPE(L,KK,N)
            ENDDO
           ENDDO
        ENDDO
      ELSE
        DO N=1,NTYPE                                                    
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)
            FLAKE(L,N) = 0.
            GS_TILE(L,N) = GS_TYPE(L,N) 
            CANHC_TILE(L,N) = CH_TYPE(L,N)                              
            VFRAC_TILE(L,N) = VF_TYPE(L,N)
            DO K=1,NSHYD
              WT_EXT_TILE(L,K,N) = WT_EXT_TYPE(L,K,N)
            ENDDO                        
          ENDDO                                                         
        ENDDO 
        do n=(lake_1-1)*nelev+1,lake_1*nelev
*D PHYSIO7A.248
          FLAKE(L,N) = 1.
*D PHYSIO7A.250
        enddo
      ENDIF
*D PHYSIO7A.268
            !don't do for pure ice points.
          if (FRAC_ICE(L).lt.1)
     &    SMCT(L) = SMCT(L) + MAX( 0. ,
*/-----
*DECLARE PPARM1A
*/-----
*D PPARM1A.26,PPARM1A.28   
     &,                      HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)
*D PPARM1A.40
     &,J,L,N_1                    ! WORK Loop counters.
*D PPARM1A.43,PPARM1A.44   
     & HT(LAND_FIELD)             ! IN Vegetation height (m).           
*D PPARM1A.46,PPARM1A.47   
     &,SATCON(LAND_FIELD)         ! IN Saturated hydraulic conductivity
!                                 !    (kg/m2/s). 
*I PPARM1A.48    
     &,INFIL_T(LAND_FIELD)        ! OUT Maximum surface infiltration 
!                                 !     rate (kg/m2/s).         
*D PPARM1A.50
*B PPARM1A.62 
      N_1=((n-1)/nelev) +1
*D PPARM1A.64,PPARM1A.68   
        Z0_T(L) = DZ0V_DH(N_1) * HT(L)
        CATCH_T(L) = CATCH0(N_1) + DCATCH_DLAI(N_1) * LAI(L)
        INFIL_T(L) = INFIL_F(N_1) * SATCON(L)                
*/-----
*DECLARE PRELIM1
*/-----
*I ABX1F405.2     
!   4.6    19/01/99    Use vegetation sampling frequencies in
!                      atmosphere timesteps, and set start times to
!                      first phenology or TRIFFID timestep.  Requires
!                      call to CNTLGEN to provide information on
!                      timestep length.  Richard Betts
*I GSS3F401.804   
*CALL CNTLGEN
*I PRELIM1.95    
! Local vectors:
      INTEGER SECS_PER_ASTEP ! Number of seconds in atmos timestep
      INTEGER A_PHENOL_STEP  ! Leaf phenology period in atmos timesteps
      INTEGER A_TRIFFID_STEP ! TRIFFID period in atmos timesteps

*I PRELIM1.102   
! 0.1  Convert PHENOL_PERIOD and TRIFFID_PERIOD to atmosphere timesteps.

      SECS_PER_ASTEP = FLOAT(SECS_PER_PERIODim(atmos_sm))/
     &                 FLOAT(STEPS_PER_PERIODim(atmos_sm))

      A_PHENOL_STEP = PHENOL_PERIOD*(86400.0/SECS_PER_ASTEP)
      A_TRIFFID_STEP = TRIFFID_PERIOD*(86400.0/SECS_PER_ASTEP)

*D ABX1F405.5,ABX1F405.6    
          ELSE IF((ITIMA.EQ.14).AND.(A_PHENOL_STEP.NE.1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_PHENOL_STEP)
*D ABX1F405.8
     &      LIST_S(st_start_time_code,NRECS)+A_PHENOL_STEP-IMD
*D ABX1F405.10,ABX1F405.11   
          ELSE IF((ITIMA.EQ.15).AND.(A_TRIFFID_STEP.NE.1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_TRIFFID_STEP)
*D ABX1F405.13
     &      LIST_S(st_start_time_code,NRECS)+A_TRIFFID_STEP-IMD
*D ABX1F405.19
              LIST_S(st_freq_code,NRECS)=A_PHENOL_STEP
*D ABX1F405.21
     &      (MOD(LIST_S(st_freq_code,NRECS),A_PHENOL_STEP).NE.0) THEN
*D ABX1F405.23
     &       'PRELIM: INCORRECT SAMPLING FOR A_PHENOL_STEP . FREQ=',
*D ABX1F405.32
              LIST_S(st_freq_code,NRECS)=A_TRIFFID_STEP
*D ABX1F405.34
     &      (MOD(LIST_S(st_freq_code,NRECS),A_TRIFFID_STEP).NE.0) THEN
*D ABX1F405.36
     &       'PRELIM: INCORRECT SAMPLING FOR A_TRIFFID_STEP . FREQ=',
*D ABX1F405.46
              LIST_S(st_freq_code,NRECS)=A_PHENOL_STEP
*D ABX1F405.48
              LIST_S(st_freq_code,NRECS)=A_TRIFFID_STEP
*/-----
*DECLARE PHENOL1A
*/-----
*D PHENOL1A.57
     & J,L,N_1                    ! Loop counters
*I PHENOL1A.61
       N_1=(n-1)/nelev +1
*D PHENOL1A.67,PHENOL1A.68
         LAI_BAL(L) = (A_WS(N_1)*ETA_SL(N_1)*HT(L)
     &               /A_WL(N_1))**(1.0/(B_WL(N_1)-1))
*D PHENOL1A.79,PHENOL1A.80
        IF (G_LEAF(L).GT.2*G_LEAF_0(N_1)) THEN
          DPHEN = -DTIME_PHEN*G_GROW(N_1)
*D PHENOL1A.84
          DPHEN = DTIME_PHEN*G_GROW(N_1)*(1.0-PHEN(L))
*/-----
*DECLARE PSLIMS1
*/-----
*I PSLIMS1.70
!
      ELSE IF(IPLAST .EQ. 101) THEN
          !New snow scheme: all surface types <times> max no. of snow levels
                  ILAST = NTILES * NSMAX
!
*D ABX1F404.16,ABX1F404.17   
!Direct vegetation parametrization: all surface tiles
        ILAST=NTILES
*/-----
*DECLARE RAD_CTL1
*/-----
*D ARE2F404.72
     &             NTILESDA,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,
     &             LAND_ALB,SICE_ALB,           
*D ARE2F404.95
     &       I,J,L,N,K,
*I  APBBF401.3
     &             downwelling_lw,
*I ARN1F404.123   
     &       NTILESDA,      ! and NTILES                                
*D ARE2F404.74,ARE2F404.75   
     &       SAL_DIM,     ! IN Set to P_FIELD for MOSES II,             
C                         !    1 otherwise                              
     &       TILE_FIELDDA,! IN Set to LAND_FIELD for MOSES II,          
C                         !    1 otherwise                              
*I GSS1F304.768   
     &      DOLR(P_FIELDDA),                                            
     &      LW_DOWN(P_FIELDDA),                                         
     &      SW_TILE(TILE_FIELDDA,NTILESDA),
     &      LAND_ALB(P_FIELDDA,NLALBS),   ! Mean land albedo            
     &      SICE_ALB(P_FIELDDA,NLALBS),   ! Mean sea-ice albedo         
     &      downwelling_lw(P_FIELDDA,BL_LEVELSDA),
*D ARE2F404.76,ARE2F404.81   
*I ARE2F404.83    
*CALL C_0_DG_C                                                          
*I ADB2F404.909   
     &     ,TILE_ALBEDO                                                 
*D ARE2F404.85,ARE2F404.86   
C zenith angle adjustment, net surface SW on tiles and downward LW      
     &      RADINCS((P_FIELDDA*(P_LEVELSDA+2+NTILES)+511)/512*512),   
*D ARE2F404.87,ARE2F404.88   
     &      LAND_ALBEDO(SAL_DIM,4),                                     
*I AWI1F403.143   
     &      FLANDG(P_FIELDDA),            ! Land fraction               
*D ARE2F404.89,ARE2F404.92   
     &     ,ALB_TILE(TILE_FIELDDA,NTILESDA,4)                           
     &     ,ICE_FRACT(P_FIELDDA)                                     
     &     ,TILE_FRAC(TILE_FIELDDA,NTYPE)                               
     &     ,TSTAR_TILE(TILE_FIELDDA,NTILESDA)                           
     &     ,SURF_DOWN_SW(P_FIELDDA,4)                                   
     &     ,T_SOL_RAD(P_FIELDDA)     ! Effective surface radiative temp
     &     ,TSTAR_SICE(P_FIELDDA)    ! Sea-ice sfc temperature (K)
     &     ,FRACSOLID(P_FIELDDA)     ! Solid surface fraction in gridbox
     &     ,SW_NET_LAND(P_FIELDDA)   ! SW net local flux over land
     &     ,SW_NET_SICE(P_FIELDDA)   ! SW net local flux over sea-ice  
     &     ,SW_NET_RTS(P_FIELDDA)    ! net SW on tiles
     &     ,SW_DOWN_RTS(P_FIELDDA)   ! down SW on tiles
*I RAD_CTL1.123   
                                                                        
! Index arrays for MOSES II                                             
      INTEGER                                                           
     & TILE_PTS(NTYPE)               ! Number of land points which      
                                     ! include the nth surface type     
     &,TILE_INDEX(TILE_FIELDDA,NTYPE)! Indices of land points which     
                                     ! include the nth surface type     
*I ADB1F401.772   
     &  ,L_SURF_DOWN_FLUX     !Logical to calculate surface downward    
!                             !fluxes                                   
*I ACN2F405.47    
     &  ,L_CTILE              !coastal tiling switch      
*I ADB1F400.70    
      REAL                                                              
     &       ALBSOLID    !Mean solid surface albedo                     
!                                                                       
*I ADB2F404.914   
      LOGICAL 
     &      LAND0P5(P_FIELDDA) ! LOCAL Land mask 
!                              !   (TRUE if land fraction >0.5)  
      REAL dummy

*D ARE2F404.96,ARE2F404.99   
!  Set index arrays and flags for MOSES II, and copy tile fractions     
!  and surface temperatures from D1                                     
      IF ( H_SECT(3) .EQ. "07A" .OR.                                    
     &     H_SECT(3) .EQ. "08A" ) THEN                                  
        CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,D1(JFRAC_TYP),   
     &               TILE_PTS,TILE_INDEX)                               
        DO N=1,NTYPE                                                    
!          DO POINT=1,TILE_PTS(N)                                        
!            L = TILE_INDEX(POINT,N)                                     
!! USEFUL TO POPULATE ALL POINTS< NOT JUST FRAC>0 ?
          DO L=LAND1,LAND1+LAND_PTS-1
            J = (N-1)*LAND_FIELD + L - 1                                
            TILE_FRAC(L,N) = D1(JFRAC_TYP+J)                            
          ENDDO                                                         
        ENDDO                                                           
        DO N=1,NTILES                                                   
          DO L=LAND1,LAND1+LAND_PTS-1                                   
            J = (N-1)*LAND_FIELD + L - 1                                
            TSTAR_TILE(L,N) = D1(JTSTAR_TYP+J)                          
          ENDDO                                                         
        ENDDO                                                           
        L_MOSES_II=.TRUE.
        L_CTILE=.TRUE.               
        L_SURF_DOWN_FLUX=.TRUE.                                         
*D ARE2F404.101
        L_MOSES_II=.FALSE.                                              
        L_SURF_DOWN_FLUX=SF(235, 1)                                     
*D RAD_CTL1.173
*D ARE2F404.104,ARE2F404.105  
*I AJS1F401.961   
        SW_NET_LAND(I) = 0.0
        SW_NET_SICE(I) = 0.0    
        LAND_ALB(I,1) = 0.0                                         
*I RAD_CTL1.180   
!  Set up GLOBAL fractional land field:                                 
      CALL LAND_TO_GLOBAL                                               
     & (D1(JLAND),D1(JFRAC_LAND),FLANDG,LAND_PTS,P_FIELDDA)   
                                                          
!  Set up FRACSOLID (solid surface fraction in grid-box),
!  ICE_FRACT (ice fraction in grid-box) and
!  LAND0P5 (set to TRUE where land fraction greater than or equal
!  to 0.5) 
      DO I=1,P_FIELD
        ICE_FRACT(I) = D1(JICE_FRACTION+I-1)   
        FRACSOLID(I) = FLANDG(I) + (1.0-FLANDG(I))*ICE_FRACT(I)
        LAND0P5(I) = .FALSE.
      ENDDO     

      DO I=1,P_FIELD                                              
        IF (FLANDG(I).GE.0.5) THEN
          LAND0P5(I) = .TRUE.
        ENDIF                        
      ENDDO  
     
*D ARE2F404.108
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C                                                  !no words for SW incs
*D ARE2F404.109,ARE2F404.110  
*D RAD_CTL1.311,ABX1F405.141  
*D @DYALLOC.3042,ARE2F404.132  
     &      D1(JLAND+JS),FLANDG(FIRST_POINT),
     &      D1(JICE_FRACTION+JS),D1(JTSTAR+JS),D1(JTSTAR_SICE+JS),      
*D ARE2F404.133
*D AJG1F405.31
     &      ALPHAM,ALPHAC,ALPHAB,DTICE,L_MOSES_II,L_SSICE_ALBEDO,       
*D ARE2F404.134
*D ARE2F404.135,GHM5F405.5    
     &      LAND_ALB(FIRST_POINT,1),SICE_ALB(FIRST_POINT,1),            
*I GHM5F405.9     
            IF ( L_MOSES_II ) THEN                                      
!-----------------------------------------------------------------------
! Calculate MOSES II tile albedos, then reset TILE_PTS and TILE_INDEX   
! and set tile fractions to 1 if aggregate tiles are used (NTILES=1).   
!-----------------------------------------------------------------------
              CALL TILE_ALBEDO (                                        
     &          P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_LIST,NTILES,     
     &          TILE_PTS,TILE_INDEX,L_SNOW_ALBEDO,D1(JSOIL_ALB),        
     &          COS_ZENITH_ANGLE,TILE_FRAC,D1(JLAI_PFT),D1(JRGRAIN_TYP),
     &          D1(JSNODEP_TYP),TSTAR_TILE,D1(JZ0_TYP),                 
     &          ALB_TILE,LAND_ALBEDO,                                 
     &          L_ESSERY_SNOW, 
     &          NSMAX, D1(JRHO_SNOW(1)),
     &          D1(JNSNOW(1)),D1(JRHO_SNOW_GRND(1))
     &          )

                                                                        
              IF (NTILES.EQ.1) THEN                                     
                TILE_PTS(1) = LAND_PTS                                  
                DO L=LAND1,LAND1+LAND_PTS-1                             
                  TILE_FRAC(L,1) = 1.                                   
                  TILE_INDEX(L+1-LAND1,1) = L                           
                ENDDO                                                   
              ELSEIF (NTILES.EQ.2*NELEV) THEN
                DO L=LAND1,LAND1+LAND_PTS-1
                  DO N=1,NTYPE
                    TILE_FRAC(L,N) = 0.
                  END DO
                  DO N=1,NTYPE-NELEV
                    J = (N-1)*LAND_FIELD + L - 1                                
                    k=mod(n-1,nelev)+1
                    TILE_FRAC(L,K) = TILE_FRAC(L,K)+D1(JFRAC_TYP+J)
                  END DO
                  DO N=NTYPE-NELEV+1,NTYPE
                    k=mod(n-1,nelev)+1
                    J = (N-1)*LAND_FIELD + L - 1                                
                    TILE_FRAC(L,NELEV+K)=TILE_FRAC(L,NELEV+K)
     &                                  +D1(JFRAC_TYP+J)
                  END DO
                END DO
                DO N=1,NTILES
                  J=0
                  DO L=LAND1,LAND1+LAND_PTS-1
                  IF (TILE_FRAC(L,N).gt.0) THEN
                    J=J+1
                    TILE_INDEX(J,N)=L
                  END IF
                  END DO
                  TILE_PTS(N)=J
                END DO

              ENDIF !NTILES                                                    
                                                                        
            ENDIF !L_MOSES_II

*D ADB1F404.2
*D ADB1F404.4,ARE2F404.137  
     &       (H_SECT(3).EQ."07A").OR.                                   
     &       (H_SECT(3).EQ."08A") ) THEN                                
*D ARE2F404.139
*D ARE2F404.140
         DO LEVEL=0,P_LEVELS+1+NTILES                                   
*D ARE2F404.142
            IF( L_MOSES_II ) THEN                                       
*D GHM5F405.76,ARE2F404.150  
     &        LAND_ALBEDO(FIRST_POINT_SAL,1), L_CTILE,                  
     &        LAND_ALB(FIRST_POINT,1), SICE_ALB(FIRST_POINT,1),
     &        FLANDG(FIRST_POINT), OPEN_SEA_ALBEDO(FIRST_POINT,1),    
     &        D1(JICE_FRACTION+JS),D1(JLAND+JS),LAND0P5(FIRST_POINT),
     &        D1(JSNODEP+JS),         
!                       MOSES II flag and array dimension               
     &        L_MOSES_II, SAL_DIM,                                      
*D ADB1F401.823,ADB1F400.198  
     &        STASHWORK(JS+SI(204,1,im_index)),
     &        STASHWORK(JS+SI(259,1,im_index)), 
     &        STASHWORK(JS+SI(260,1,im_index)), L_FLUX_BELOW_690NM_SURF,
     &        STASHWORK(JS+SI(235,1,im_index)), L_SURF_DOWN_FLUX,       
*I ADB2F404.1008  
     &        SURF_DOWN_SW(FIRST_POINT,1),                              
*I RAD_CTL1.472   
          DO BAND=1,4                                                   
            DO POINT = FIRST_POINT, LAST_POINT                          
              SURF_DOWN_SW(POINT,BAND) = 0.                             
            ENDDO                                                       
          ENDDO                                                         
*D ARE2F404.151
CL Set up downward surface SW components if required for MOSES II       
!   Calculate SW_NET_RTS, SW_DOWN_RTS and LAND_ALBEDO
*D ARE2F404.153,ARE2F404.155  
      IF ( L_MOSES_II ) THEN
*I ARE2F404.156   
        DO I=FIRST_POINT, LAST_POINT
          LAND_ALB(I,1) = 0.0
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_LIST(L)
          SW_DOWN_RTS(I) = 0.0
          SW_NET_RTS(I) = 0.0
        END DO                                        
        DO N=1,NTILES                                                   
          DO I=FIRST_POINT,LAST_POINT                                   
            J = I + (P_LEVELS + 1 + N)*P_FIELD                          
            RADINCS(J) = 0.                                             
          ENDDO                                                         
        ENDDO
               
        DO N=1,NTILES                                                 
          DO POINT=1,TILE_PTS(N)                                      
            L = TILE_INDEX(POINT,N)                                   
            I = LAND_LIST(L)                                          
            J = I + (P_LEVELS + 1 + N)*P_FIELD        
            DO BAND=1,4                                                 
              RADINCS(J) = RADINCS(J) + (1. - ALB_TILE(L,N,BAND)) *     
     &                                              SURF_DOWN_SW(I,BAND)
            ENDDO
            SW_NET_RTS(I)=SW_NET_RTS(I)+RADINCS(J)*TILE_FRAC(L,N)
          ENDDO                                                         
        ENDDO

        IF (L_CTILE) THEN
          DO L=LAND1,LAND1+LAND_PTS-1
            I = LAND_LIST(L)
            DO BAND=1,4
              SW_DOWN_RTS(I) = SW_DOWN_RTS(I) + SURF_DOWN_SW(I,BAND)
            ENDDO
            IF (SW_DOWN_RTS(I).GT.0.0) THEN
              LAND_ALB(I,1) = 1.0-SW_NET_RTS(I)/SW_DOWN_RTS(I)
            ENDIF 
          ENDDO
        ENDIF
                  
      ENDIF                                                        
*D ARE2F404.157,ARE2F404.158  
CL and net surface on tiles                                             
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C                                                                       
*D RAD_CTL1.622,RAD_CTL1.623  
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = NETSW(I)                                    
*D RAD_CTL1.626
!  Mutliply sea field by sea fraction. 
     &                   *(1.-FLANDG(I))       
        END DO 
      ELSE
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = NETSW(I)                                    
     &                   - MEAN_COSZ(I) * RADINCS(I)                    
     &                   - STASHWORK(SI(203,1,im_index)+I-1)  
        ENDDO
      ENDIF                                                   
*D RAD_CTL1.695,RAD_CTL1.697  
! Only set on sea-ice points if MOSES2:                                 
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
*D ARE2F404.159,ARE2F404.160  
! land_alb is only set on points where there is incoming sw radiation   
! at the previous timestep, therefore it will be zero over some         
! land points                                                           
          IF (FRACSOLID(I).GT.0.0) THEN
*D ABX1F405.142
            IF (FLANDG(I).GT.0.0.AND.LAND_ALB(I,1).LE.0.0) THEN
              SW_NET_LAND(I) = 0.0
              SW_NET_SICE(I) = 0.0                                      
            ELSE                                         
              ALBSOLID = ( FLANDG(I) * LAND_ALB(I,1) +
     &          (1.0-FLANDG(I)) * ICE_FRACT(I) * SICE_ALB(I,1) )        
     &            /FRACSOLID(I)                                         
*D ABX1F405.144,ABX1F405.145  
              IF (FLANDG(I).GT.0.0) THEN
                SW_NET_LAND(I) = RADINCS(I)                             
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)                
     &              * (1.0-LAND_ALB(I,1))/(1.0-ALBSOLID)                
              ENDIF
            
              IF (ICE_FRACT(I).GT.0.0) THEN
                SW_NET_SICE(I) = RADINCS(I)                             
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)                
     &              * (1.0-SICE_ALB(I,1))/(1.0-ALBSOLID)                
                                                                        
                SURF_RADFLUX(I) = SW_NET_SICE(I) * ICE_FRACT(I)
              ENDIF                                                     
            ENDIF 
                                                                      
          ENDIF 
        ENDDO                                                           
      ELSE                                                              
*D ABX1F405.147,ABX1F405.160  
          SURF_RADFLUX(I) = RADINCS(I) * COS_ZENITH_ANGLE(I)            
        END DO                                                          
*I RAD_CTL1.698   
      IF ( L_MOSES_II ) THEN                                            
        DO N=1,NTILES                                                   
          DO POINT=1,TILE_PTS(N)                                        
            L = TILE_INDEX(POINT,N)                                     
            I = LAND_LIST(L)                                            
            J = I + (P_LEVELS + 1 + N)*P_FIELD                          
            SW_TILE(L,N) = RADINCS(J) * COS_ZENITH_ANGLE(I)             
          END DO                                                        
        END DO                                                          
      ENDIF                                                             
                                                                        
*D RAD_CTL1.705
       IF(SF(202,1)) THEN                                              
*D RAD_CTL1.712,RAD_CTL1.715  
*I RAD_CTL1.716   
          IF (L_CTILE) THEN                                      
            DO I = FIRST_POINT,LAST_POINT                               
              STASHWORK(SI(201,1,im_index)+I-1)=RADINCS(I)*MEAN_COSZ(I)+
!  Multiply sea field by sea fraction  
     &        STASHWORK(SI(203,1,im_index)+I-1)*(1.-FLANDG(I))          
            END DO                                                      
          ELSE
            DO I = FIRST_POINT,LAST_POINT                               
              STASHWORK(SI(201,1,im_index)+I-1)=RADINCS(I)*MEAN_COSZ(I)+
     &        STASHWORK(SI(203,1,im_index)+I-1)
            END DO
          ENDIF

        END IF 

C Calculate new diagnostics for down sw land and sea-ice portions: 
    
        IF(SF(257,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(257,1,im_index)+I-1) = SW_NET_LAND(I)          
          ENDDO                                                         
        END IF                                                          
                                                                        
        IF(SF(258,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(258,1,im_index)+I-1) = SW_NET_SICE(I)          
          ENDDO                                                         
*D ARE2F404.174
        OFFSET=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512              
*D ARE2F404.175
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512  !no words for L
*I ADB1F400.144
        DO LEVEL=1,BL_LEVELS
          DO I=1,P_FIELD
            downwelling_lw(I,level)=0.0
          ENDDO
        ENDDO
*D ARE2F404.178
          T_SOL_RAD(I) = 0.0 
          TSTAR_SICE(I) = D1(JTSTAR_SICE+I-1)
*D ARE2F404.180
        IF ( L_CTILE ) THEN
          DO I=1,P_FIELD                                            
            T_SOL_RAD(I) = (1.0-FLANDG(I))*                             
     &            ICE_FRACT(I)*TSTAR_SICE(I)**4                         
          ENDDO
          DO N=1,NTILES                                                 
            DO POINT=1,TILE_PTS(N)    
              L = TILE_INDEX(POINT,N)                                   
              I = LAND_LIST(L)                                          
              T_SOL_RAD(I) = T_SOL_RAD(I) + 
     &            FLANDG(I)*TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
            ENDDO
          ENDDO                     
          DO  I=1,P_FIELD
            IF (FRACSOLID(I).GT.0.0) THEN
              T_SOL_RAD(I)=(T_SOL_RAD(I)/FRACSOLID(I))**0.25            
            ENDIF                                                       
          ENDDO
                                                                        
        ELSEIF ( L_MOSES_II ) THEN                                      
*D ARE2F404.183,ARE2F404.184  
            T_SOL_RAD(I) = 0.                                           
*D ARE2F404.186,ABX1F405.163  
          DO N=1,NTILES                                                 
            DO POINT=1,TILE_PTS(N)                                      
              L = TILE_INDEX(POINT,N)                                   
*D ARE2F404.189,ARE2F404.191  
              T_SOL_RAD(I) = T_SOL_RAD(I) + TILE_FRAC(L,N) *            
     &                                      TSTAR_TILE(L,N)**4          
*D ARE2F404.196
            T_SOL_RAD(I) = T_SOL_RAD(I)**0.25                           
*I ARE2F404.197   
        ENDIF                                                           
                                                                        
        IF ( L_MOSES_II ) THEN                                          
           L_SURF_DOWN_FLUX=.TRUE.                                      
        ELSE                                                            
           L_SURF_DOWN_FLUX=SF(207, 2)                                  
*D ARE2F404.200
         DO LEVEL=0,P_LEVELS+1+NTILES                                   
*I APBBF401.69    
           LW_DOWN(I)=0.0                                               
*D ARE2F404.201
     &      D1(JP_EXNER(1)+JS_LOCAL(I)),T_SOL_RAD(FP_LOCAL(I)),         
*D ARE2F404.202
     &        D1(JP_EXNER(1)+JS_LOCAL(I)),D1(JTSTAR+JS_LOCAL(I)),
     &        T_SOL_RAD(FP_LOCAL(I)),D1(JTSTAR_SEA+JS_LOCAL(I)),
     &        L_CTILE,        
*D ADB1F400.322
C Only want the 0.5 threshold LAND mask and fractional land:            
     &        LAND0P5(FP_LOCAL(I)),FLANDG(JS_LOCAL(I)+1),           
*D ADB1F400.334
     &        LW_DOWN(FP_LOCAL(I)),L_SURF_DOWN_FLUX,                    
*I ADB1F400.355
     &        ,BL_LEVELS,downwelling_lw(FP_LOCAL(I),1)
*D ARE2F404.203
C Store downward surface LW radiation flux if required for MOSES II     
C DOLR is TOA outward LW - surface upward LW for land and sea-ice       
*D ARE2F404.205,ARE2F404.207  
      IF ( L_MOSES_II ) THEN                                            
                                                                        
        DO I=FIRST_POINT,LAST_POINT
          DOLR(I) = OLR(I) 
C         If no land in grid box and some sea-ice...                    
          IF (FLANDG(I).EQ.0.0 .AND. ICE_FRACT(I).GT.0.0) THEN 
            DOLR(I) = DOLR(I) - ICE_FRACT(I)*SBCON*TSTAR_SICE(I)**4     
          ENDIF                                                         
        ENDDO                                                           
                                                                        
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          I = LAND_LIST(L)                                              
          DOLR(I) = DOLR(I) - FRACSOLID(I)*SBCON*T_SOL_RAD(I)**4 
        ENDDO                                                           
                                                                        
        DO I=FIRST_POINT,LAST_POINT                                     
          J = I + (P_LEVELS + 2)*P_FIELD + OFFSET                       
          RADINCS(J) = LW_DOWN(I)                                       
          RADINCS(J+P_FIELD) = DOLR(I)                                  
        ENDDO                                                           
                                                                        
      ENDIF                                                             
*D RAD_CTL1.906,RAD_CTL1.907  
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = - OLR(I)                                    
*D RAD_CTL1.910
!  Mutliply sea field by sea fraction
     &                   *(1.-FLANDG(I))     
        END DO
      ELSE
          NET_ATM_FLUX(I) = - OLR(I)                                    
     &                   - RADINCS(I+OFFSET)                            
     &                   - STASHWORK(SI(203,2,im_index)+I-1)
      ENDIF                                                            
*I RAD_CTL1.940   
      IF (L_MOSES_II) THEN                                              
        DO I=FIRST_POINT,LAST_POINT                                     
          SURF_RADFLUX(I) = SURF_RADFLUX(I) +                           
     &                      D1(JICE_FRACTION+I-1)*LW_DOWN(I)            
        ENDDO                                                           
      ELSE                                                              
*I RAD_CTL1.943   
      ENDIF                                                             
*D ARE2F404.211,ARE2F404.235  
*D RAD_CTL1.964,RAD_CTL1.967  
          IF  (L_CTILE) THEN
            DO I = FIRST_POINT,LAST_POINT                               
                STASHWORK(SI(201,2,im_index)+I-1) = RADINCS(I+OFFSET)+
!  Multiply sea field by sea fraction    
     &          STASHWORK(SI(203,2,im_index)+I-1)*(1.-FLANDG(I)) 
            END DO 
          ELSE
            DO I = FIRST_POINT,LAST_POINT                               
                STASHWORK(SI(201,2,im_index)+I-1) = RADINCS(I+OFFSET)+  
     &          STASHWORK(SI(203,2,im_index)+I-1)
            END DO
          ENDIF                                                       
*I RAD_CTL1.1005  
        IF(SF(207,2)) THEN                                              
                                                                        
CL  Downward Surface radiative flux                                     
                                                                        
          CALL COPYDIAG (STASHWORK(SI(207,2,im_index)),LW_DOWN,         
     &        FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                
     &        im_ident,2,207,                                           
*CALL ARGPPX                                                            
     &        ICODE,CMESSAGE)                                           
                                                                        
          IF (ICODE .GT. 0) RETURN                                      
                                                                        
        END IF                                                          
                                                                        
*D ABX1F405.185
        CALL SWAPB_LAND(SW_TILE,LAND_FIELD,P_FIELD,                     
     &                  ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,              
     &                  NTILES,LAND_LIST)                               
        CALL SWAPBOUNDS(downwelling_lw(1,1),ROW_LENGTH,P_ROWS,
     &                  EW_Halo,NS_Halo,bl_levels)
        CALL SWAPBOUNDS(LW_DOWN,ROW_LENGTH,P_ROWS,                      
*D ABX1F405.187
        CALL SWAPBOUNDS(DOLR,ROW_LENGTH,P_ROWS,                         
*I ABX1F405.188   
      ENDIF
      IF (L_CTILE) THEN
        CALL SWAPBOUNDS(LAND_ALB(1,1),ROW_LENGTH,P_ROWS, 
     &                  EW_Halo,NS_Halo,NLALBS)                         
        CALL SWAPBOUNDS(SICE_ALB(1,1),ROW_LENGTH,P_ROWS,
     &                  EW_Halo,NS_Halo,NLALBS)        
*/-----
*DECLARE READFL1A
*/-----
*I GSI1F405.366   
      ELSE IF (LOOKUP(LBPACK,K).EQ.0.AND.
     &         ppxref_grid_type.EQ.ppx_atm_compressed) THEN
! Grid code suggests land-only field but field is unpacked,
! so use grid code for land data on full field
        fake_D1_ADDR(d1_grid_type)=ppx_atm_tland


*/-----
*DECLARE READLSA1
*/-----
*I READLSA1.88
C
*CALL NSTYPES
C
*CALL C_LAND_CC
C
      NAMELIST /LAND_CC/ F0,LAI_MIN,NL0,R_GROW,
     +TLOW,TUPP,Q10,V_CRIT_ALPHA,KAPS
C
*CALL C_EDDY
C
      NAMELIST /EDDY/ FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH

*B AJX3F405.151
     &  ,DZSNOW
*I ADR1F305.198
      DO J = 1,10
         DZSNOW(J) = 1.
      ENDDO
*B READLSA1.132
      WRITE(6,*) 'Print variables from namelists'
      REWIND 5
      READ(5,LAND_CC)
      REWIND 5
      READ(5,EDDY)
      REWIND 5
      WRITE(6,16) 'F0:       ',F0(1),F0(3),F0(3),F0(4),F0(5)
      WRITE(6,16) 'LAI_MIN:  ',LAI_MIN(1),LAI_MIN(2),LAI_MIN(3)
     +                        ,LAI_MIN(4),LAI_MIN(5)
      WRITE(6,16) 'NL0:      ',NL0(1),NL0(2),NL0(3),NL0(4),NL0(5)
      WRITE(6,16) 'R_GROW:   ',R_GROW(1),R_GROW(2),R_GROW(3)
     +                        ,R_GROW(4),R_GROW(5)
      WRITE(6,16) 'TLOW:     ',TLOW(1),TLOW(2),TLOW(3)
     +                        ,TLOW(4),TLOW(5)
      WRITE(6,16) 'TUPP:     ',TUPP(1),TUPP(2),TUPP(3)
     +                        ,TUPP(4),TUPP(5)
      WRITE(6,22) 'Q10:      ',Q10 
      WRITE(6,22) 'V_CRIT_ALPHA:   ',V_CRIT_ALPHA
      WRITE(6,21) 'KAPS:     ',KAPS
      WRITE(6,23) 'RHCRIT:     ',RHCRIT(1),RHCRIT(2)
     +                             ,RHCRIT(3),RHCRIT(4)
     +                             ,RHCRIT(5),RHCRIT(6)
     +                             ,RHCRIT(7),RHCRIT(8)
     +                             ,RHCRIT(9),RHCRIT(10)
     +                             ,RHCRIT(11)
      WRITE(6,22) 'VF1:     ',VF1
      WRITE(6,21) 'CT:     ',CT
      WRITE(6,21) 'CW_SEA:     ',CW_SEA
      WRITE(6,21) 'CW_LAND:     ',CW_LAND
      WRITE(6,21) 'KAY_GWAVE:     ',KAY_GWAVE
      WRITE(6,21) 'KAY_LEE_GWAVE:     ',KAY_LEE_GWAVE
      WRITE(6,21) 'Z0FSEA:     ',Z0FSEA
      WRITE(6,22) 'ALPHAM:     ',ALPHAM
      WRITE(6,24) 'DIFF_COEFF:     ',DIFF_COEFF(1)
     +                              ,DIFF_COEFF(2)
     +                              ,DIFF_COEFF(3)
     +                              ,DIFF_COEFF(4)
     +                              ,DIFF_COEFF(5)
     +                              ,DIFF_COEFF(6)
     +                              ,DIFF_COEFF(7)
     +                              ,DIFF_COEFF(8)
     +                              ,DIFF_COEFF(9)
     +                              ,DIFF_COEFF(10)
     +                              ,DIFF_COEFF(11) 
      WRITE(6,24) 'DIFF_COEFF_Q:   ',DIFF_COEFF_Q(1)
     +                              ,DIFF_COEFF_Q(2)
     +                              ,DIFF_COEFF_Q(3)
     +                              ,DIFF_COEFF_Q(4)
     +                              ,DIFF_COEFF_Q(5)
     +                              ,DIFF_COEFF_Q(6)
     +                              ,DIFF_COEFF_Q(7)
     +                              ,DIFF_COEFF_Q(8)
     +                              ,DIFF_COEFF_Q(9)
     +                              ,DIFF_COEFF_Q(10)
     +                              ,DIFF_COEFF_Q(11)
      WRITE(6,21) 'FNUB_SI:     ',FNUB_SI    
      WRITE(6,21) 'KAPPA0_SI:     ',KAPPA0_SI   
      WRITE(6,21) 'DKAPPA_DZ_SI:     ',DKAPPA_DZ_SI
      WRITE(6,21) 'FNU0_SI:     ',FNU0_SI    
      WRITE(6,21) 'STABLM_SI:     ',STABLM_SI
      WRITE(6,21) 'AHI1_SI:     ',AHI1_SI
      WRITE(6,21) 'AHI2_SI:     ',AHI2_SI
      WRITE(6,21) 'AHI3_SI:     ',AHI3_SI
      WRITE(6,21) 'SLOPE_MAX:     ',SLOPE_MAX
      WRITE(6,21) 'ATHKDF1_SI:     ',ATHKDF1_SI
      WRITE(6,21) 'ATHKDF2_SI:     ',ATHKDF2_SI
      WRITE(6,21) 'ATHKDF3_SI:     ',ATHKDF3_SI
      WRITE(6,21) 'AM0_SI:     ',AM0_SI
      WRITE(6,21) 'AM1_SI:     ',AM1_SI
      WRITE(6,21) 'AH_SI:     ',AH_SI
      WRITE(6,21) 'BM_SI:     ',BM_SI
      WRITE(6,21) 'CRIT_RI:     ',CRIT_RI
      WRITE(6,21) 'MAX_QLARGE_DEPTH:     ',MAX_QLARGE_DEPTH
C
 16    format(A13,5F9.3)
 21    format(A13,e12.4)
 22    format(A20,F5.3)
 23    format(A13,11F9.3)
 24    format(A13,11e12.4)
C
*/-----
*DECLARE ROOTFR7A
*/-----
*D ROOTFR7A.72
      PARAMETER (P=1.0)
*/-----
*DECLARE RPANCA1A
*/-----
*I GDG0F401.1361  
     &                     FLAND_CTILE,                                 
     &                     TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            
     &                     TSTAR_SICE_CTILE,                            

*D RPANCA1A.111,RPANCA1A.113  
     &       ICE_FRACTION(P_FIELD), 
!                               !INOUT  Ice frac of sea part of grid
!                               !       box, updated if requested   
*IF -DEF,RECON                                                     
     &       FLAND_CTILE(LAND_FIELD),                                   
!                                !IN  Fractional land on land pts.   
*ENDIF                                                                  
     &       FLAND_G(P_FIELD),   !WORK Frac land over all points.       
     &       TSTAR(P_FIELD),     !INOUT  TSTAR:updated if requested     
     &       TSTAR_LAND_CTILE(P_FIELD),                                 
!                                !INOUT  as above, but for land.        
     &       TSTAR_SEA_CTILE(P_FIELD),                                  
!                                !INOUT  as above, but for open sea.    
     &       TSTAR_SICE_CTILE(P_FIELD),                                 
!                                !INOUT  as above, but for sea-ice.     
*D RPANCA1A.120
     &       LAND(P_FIELD),      ! WORK LAND mask                       
     &       SEA(P_FIELD),       ! WORK SEA mask                        
     &       LTSTAR_SICE         ! IN TRUE if TSTAR_SICE has been read i
!                                ! from input dump.                     
!                                ! If FALSE set to TSTAR_SEA.           
*I RPANCA1A.132   
*CALL CNTLATM
*I RPANCA1A.165   
     &,      TSTAR_LAND(P_FIELD)!Temporary store for land surface temp.
     &,      TSTAR_SEA(P_FIELD) !as above, but for open sea.           
     &,      TSTAR_SICE(P_FIELD)!as above, but for sea-ice.            
     &,      TSTAR_SSI(P_FIELD) !as above, but for sea mean.           
*I UDG4F402.255   
     &       L,                 ! Land index                        
*I GRS2F404.12    
     &       ,L_CTILE           ! Coastal tiling switch
C                               ! always true for 7A/8A BL
*I RPANCA1A.266   
!     Set coastal tiling flag
      IF ( H_SECT(3).EQ.'07A' .OR. H_SECT(3).EQ.'08A') THEN
        L_CTILE = .TRUE.
      ENDIF                                                 
                                                                        
!     Set up surface temperatures:                                      
                                                                        
       IF(L_CTILE)THEN                                                  
         DO I=1,P_FIELD                                                 
            TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)                           
            TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)                             
            TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)                           
            IF(ICE_FRACTION(I).LE.0.0)THEN                              
              TSTAR_SSI(I)=TSTAR_SEA(I)                                 
            ELSE                                                        
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                
     &          +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)                     
            ENDIF                                                       
         ENDDO                                                          
       ELSE                                                             
         DO I=1,P_FIELD                                                 
            TSTAR_LAND(I)=TSTAR(I)                                      
            TSTAR_SSI(I)=TSTAR(I)                                       
         ENDDO                                                          
       ENDIF                                                            
                                                                        
                                                                        
*I RPANCA1A.332   
! Read in fractional land field                                         
*IF DEF,RECON                                                      
       IF(L_CTILE)THEN                                                  
        FILE=FILEANCIL(48)                                             
        NFTIN=FTNANCIL(FILE)                                            
        WRITE(6,*)'READING IN LAND FRACTION'                            
        CALL READFLDS(NFTIN,1,NLOOKUP(48),LOOKUP(1,LOOKUP_START(FILE)),
     &                LEN1_LOOKUP,ANCIL1,P_FIELD,FIXHD(1,FILE),         
*CALL ARGPPX                                              
     &                      ICODE,CMESSAGE)                             
        WRITE(6,*)'READ IN LAND FRACTION'                               
                                                                        
      DO I=1,P_FIELD                                                    
        FLAND_G(I)=0.0                                                  
        IF(LAND(I))FLAND_G(I)=ANCIL1(I)                                 
! If land or sea fraction is less than machine tolerance print warning  
          IF(LAND(I).AND.FLAND_G(I).LE.1.42E-14)THEN                
           WRITE(6,*)'*****************WARNING********************'     
           WRITE(6,*)'LAND FRACTION IS LESS THAN MACHINE TOLERANCE'     
          ENDIF                                                         
          IF(.NOT.LAND(I).AND.1.0-FLAND_G(I).LE.1.42E-14)THEN       
           WRITE(6,*)'*****************WARNING********************'     
           WRITE(6,*)'SEA FRACTION IS LESS THAN MACHINE TOLERANCE'      
          ENDIF                                                         
!                                                                       
          IF(FLAND_G(I).LE.0.0.AND.LAND(I))THEN                         
           WRITE(6,*)'*ERROR* a) LAND FRAC & LAND MASK ARE INCONSISTENT'
           ICODE = 800                                                  
           CMESSAGE='REPLANCA:ERROR:LAND FRAC & MASK ARE INCONSISTENT'  
          ENDIF                                                         
          IF(FLAND_G(I).GT.0.0.AND..NOT.LAND(I))THEN                    
           WRITE(6,*)'*ERROR* b) LAND FRAC & LAND MASK ARE INCONSISTENT'
           ICODE = 801                                                  
           CMESSAGE='REPLANCA:ERROR:LAND FRAC & MASK ARE INCONSISTENT'  
          ENDIF                                                         
         ENDDO                                                          
                                                                        
        ELSE                     ! Not coastal tiling:                  
         DO I=1,P_FIELD                                                 
          IF(LAND(I))THEN                                               
            FLAND_G(I)=1.0                                              
          ELSE                                                          
            FLAND_G(I)=0.0                                              
          ENDIF                                                         
         ENDDO                                                          
        ENDIF                                                           
*ENDIF                                                                 
*IF -DEF,RECON                                                     
! Set up global fractional land field                                   
         IF(L_CTILE)THEN                                                
           L=0                                                          
           DO I=1,P_FIELD                                               
             FLAND_G(I)=0.0                                             
             IF(LAND(I))THEN                                            
               L=L+1                                                    
               FLAND_G(I)=FLAND_CTILE(L)                                
            ENDIF                                                       
           ENDDO                                                        
         ELSE                                                           
           DO I=1,P_FIELD                                               
             IF(LAND(I))THEN                                            
               FLAND_G(I)=1.0                                           
             ELSE                                                       
               FLAND_G(I)=0.0                                           
            ENDIF                                                       
           ENDDO                                                        
         ENDIF                                                          
!                                                                       
*ENDIF                                                                  
                                                                        
        DO I=1,P_FIELD                                                  
          SEA(I)=.FALSE.                                                
          IF(FLAND_G(I).LT.1.0)SEA(I)=.TRUE.                            
        ENDDO                                                           
*D GRS2F404.52
            IF(SEA(I)) THEN                                             
*D RPANCA1A.874
c      write(6,*)'cmt rpanca1a.873 ',field,lookup_step(field),i2,i1,    
c     +lookup_start(file),nlookup(field),level,nlookups,                
c     +lookup(item_code,i1+LOOKUP_START(FILE)-1),                       
c     +lookup(item_code,i2+LOOKUP_START(FILE)-1)                        
cc                                                                      
cc      do i=1,45                                                       
cc       write(6,*)'cmtloop',i,lookup(i,i1+lookup_start(file)-1),       
cc     +lookup(i,i2+lookup_start(file)-1)                               
cc      enddo                                                           
cmt                                                                     
cmt  Current code checks that data doesn't go beyond end of file        
cmt  but doesn't check that data is the same type. Additional check     
cmt  made which if failed, forces code to go back to start of file.     
cmt  This avoids problems at year end.                                  
cmt                                                                     
        IF (I1.LE.FIXHD(152,FILE)                                       
     +  .AND. (LOOKUP(ITEM_CODE,I1+LOOKUP_START(FILE)-1) .EQ.           
     +         STASHANCIL(FIELD)) ) THEN                                
cmt                                                                     
*I GRS2F404.186   
            IF(.NOT.LTLEADS)THEN                                        
*I RPANCA1A.968   
            ELSE                                                        
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           
     &              TIME,P_FIELD)                                       
            ENDIF                                                       
*D RPANCA1A.1025
            IF(SEA(I).AND.FIELD.EQ.26) THEN                             
*D RPANCA1A.1040
                IF(TSTAR_LAND(I).GT.TM) TSTAR_LAND(I)=TM                
*D RPANCA1A.1053,RPANCA1A.1054 
                IF(TSTAR_LAND(I).GT.TM.AND.ANCIL_DATA(I).GT.0.0) THEN   
                  TSTAR_LAND(I)=TM                                      
*D RPANCA1A.1074
            IF (SEA(I)) THEN                                            
*D RPANCA1A.1081,RPANCA1A.1085 
          IF(.NOT.LTLEADS)THEN                                          
            DO I=1,P_FIELD                                              
              IF(ICE_FRACTION(I).GT.0.0) THEN                           
                TSTAR_SSI(I)=AMIN1(TSTAR_SSI(I),TFS)                    
              ENDIF                                                     
            END DO                                                      
          ENDIF                                                         
*D RPANCA1A.1101,RPANCA1A.1102 
            IF(L_CTILE.OR.ICE_FRACTION(I).EQ.0.0)THEN                   
              IF (SEA(I)) THEN                                          
                IF (L_SSTANOM) THEN                                     
*D RPANCA1A.1104
                TSTAR_SEA(I)=ANCIL_DATA(I)+TSTAR_ANOM(I)                
*D RPANCA1A.1106
                TSTAR_ANOM(I)=TSTAR_SEA(I)-ANCIL_DATA(I)                
*D RPANCA1A.1108,RPANCA1A.1109 
                ELSE                                                    
                  TSTAR_SEA(I)=ANCIL_DATA(I)                            
                END IF                                                  
                IF(ICE_FRACTION(I).EQ.0.0)TSTAR_SSI(I)=TSTAR_SEA(I)     
*D GJT1F304.121
            IF (SEA(I)) THEN                                            
*D TJ240293.40
              D1(D1_ANCILADD(FIELD)+I-1)=TSTAR_LAND(I)                  
                                                                        
*D RPANCA1A.1121
            IF(SEA(I)) THEN                                             
*D RPANCA1A.1141
            IF(SEA(I)) THEN                                             
*I RPANCA1A.1161  
      IF(L_CTILE)THEN                                                   
        DO I=1,P_FIELD                                                  
          IF(SEA(I).AND.ICE_FRACTION(I).GT.0.0)THEN                     
            IF(LTLEADS.OR.LAMIPII)THEN                                  
                                                                        
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                
     &          +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)                      
                                                                        
            ELSE                                                        
                                                                        
              TSTAR_SEA(I)=TFS                                          
              TSTAR_SICE(I)=(TSTAR_SSI(I)                               
     &          -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)     
                                                                        
            ENDIF                                                       
          ENDIF                                                         
C                                                                       
          TSTAR(I)=FLAND_G(I)*TSTAR_LAND(I)                             
     &      +(1.-FLAND_G(I))*TSTAR_SSI(I)                               
        ENDDO                                                           
      ELSE                                                              
        DO I=1,P_FIELD                                                  
          IF(LAND(I))THEN                                               
            TSTAR(I)=TSTAR_LAND(I)                                      
          ELSE                                                          
            TSTAR(I)=TSTAR_SSI(I)                                       
          ENDIF                                                         
        ENDDO                                                           
      ENDIF                                                             
                                                                        
!     Set up surface temperatures:                                      
      IF(L_CTILE)THEN                                                   
        DO I=1,P_FIELD                                                  
          TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)                             
          TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)                               
          TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)                             
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
*/-----
*DECLARE SCREEN7A
*/-----
*D SCREEN7A.24,SCREEN7A.26   
*D SCREEN7A.36,SCREEN7A.37   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                        
*D SCREEN7A.39
     & TILE_FRAC,TL_1,TSTAR_SSI,TSTAR_TILE,                             
*D SCREEN7A.41
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE                                
     & ,nelev,p_elev,QIM_ELEV,TIM_ELEV,z1_elev
*D SCREEN7A.54
     &,NTILES               ! IN Number of tiles per land point.        
*D SCREEN7A.56
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SCREEN7A.58
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
*D SCREEN7A.61,SCREEN7A.62   
     & SQ1P5                ! IN STASH flag for 1.5-metre sp humidity.  
*D SCREEN7A.66
     & FLANDG(P_FIELD)      ! IN Fraction of gridbox which is land.     
     &,CHR1P5M(LAND_FIELD,NTILES)                                       
*D SCREEN7A.74
     &,RESFT(LAND_FIELD,NTILES)                                         
*D SCREEN7A.76
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*D SCREEN7A.80,SCREEN7A.81   
     &,TSTAR_SSI(P_FIELD)   ! IN Sea/sea-ice mean sfc temperature (K).  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SCREEN7A.85
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
*D SCREEN7A.89
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
*I SCREEN7A.95    
     &,Q1P5M_TILE(LAND_FIELD,NTILES)                                    
!                           ! OUT Q1P5M over land tiles.                
*I SCREEN7A.97    
     &,T1P5M_TILE(LAND_FIELD,NTILES)                                    
!                           ! OUT T1P5M over land tiles.                
*D SCREEN7A.105,SCREEN7A.106  
*I SCREEN7A.112
      INTEGER NELEV,K
      REAL p_elev(land_field,nelev),z1_elev(LAND_FIELD,NELEV),
     & tim_elev(land_field,nelev),qim_elev(land_field,nelev)
*D SCREEN7A.129,SCREEN7A.131  
          IF (FLANDG(I).LT.1.0 ) THEN                                 
            T1P5M(I) = (1.-FLANDG(I))*                              
     &        (TSTAR_SSI(I) - GRCP*Z1P5M +                            
     &        CHR1P5M_SICE(I) *  (TL_1(I) - TSTAR_SSI(I) +        
     &          GRCP*(Z1(I)+Z0M(I)-Z0H(I))))                      
*D SCREEN7A.135
        DO N=1,NTILES                                                   
          K=mod(n-1,nelev)+1
          DO L=1,LAND_FIELD                                             
            T1P5M_TILE(L,N) = 0.                                        
          ENDDO                                                         
*D SCREEN7A.139,SCREEN7A.141  
            T1P5M_TILE(L,N) = TSTAR_TILE(L,N) - GRCP*Z1P5M +            
     &            CHR1P5M(L,N)*( TIM_ELEV(L,K) - TSTAR_TILE(L,N) +
     &            GRCP*(Z1_ELEV(L,K)+Z0M_TILE(L,N)-Z0H_TILE(L,N)) )
*D SCREEN7A.142
            T1P5M(I) = T1P5M(I)                                     
     &        + FLANDG(I)*TILE_FRAC(L,N)*T1P5M_TILE(L,N)              
*D SCREEN7A.153
        CALL QSAT(QS(P1),TSTAR_SSI(P1),PSTAR(P1),P_POINTS)              
*D SCREEN7A.156
          IF (FLANDG(I).LT.1.0 ) THEN                                 
*D SCREEN7A.158
            Q1P5M(I) = (1.-FLANDG(I))*                              
     &        (QW_1(I) + CER1P5M*( QW_1(I) - QS(I) ))             
*D SCREEN7A.167
        DO N=1,NTILES                                                   
          k=mod(n-1,nelev)+1
          DO L=1,LAND_FIELD                                             
            Q1P5M_TILE(L,N) = 0.                                        
          ENDDO                                                         
*D SCREEN7A.169
     &              P_ELEV(LAND1,K),LAND_PTS)
*D SCREEN7A.174,SCREEN7A.175  
            Q1P5M_TILE(L,N) = QIM_ELEV(L,K) + 
     &                        CER1P5M*( QIM_ELEV(L,K) - QS_TILE(L) )
            Q1P5M(I) = Q1P5M(I)                                     
     &        + FLANDG(I)*TILE_FRAC(L,N)*Q1P5M_TILE(L,N)
*/-----
*DECLARE SEED
*/-----
*D SEED.4,SEED.5    
     + FRAC_MIN                   ! Minimum areal fraction for PFTs.    
     +,FRAC_SEED                  ! "Seed" fraction for PFTs.
      PARAMETER(FRAC_MIN = 1.0E-6, FRAC_SEED = 0.01)                 
*/-----
*DECLARE SETMODL1
*/-----
*I GAV0F405.6     
!   4.6     22/01/99   Remove lines which overwrite H_VERS to 0 for
!                      section 19 for all internal models.
!                      Richard Betts
*D GSS3F401.1081,GSS3F401.1083 
*/-----
*DECLARE SFEVAP7A
*/-----
*D SFEVAP7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFEVAP7A.50,SFEVAP7A.55   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,NSHYD,LTIMER,FLAND,               
     & ASHTF_TILE,CANOPY,DTRDZ_1,FLAKE,FRACA,SNOW_TILE,RESFS,           
     & RESFT,RHOKH_1,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,                
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                        
     & ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT           
*D SFEVAP7A.68
     &,NTILES                ! IN Number of tiles per land point.       
*D SFEVAP7A.70
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SFEVAP7A.72
     &,TILE_PTS(NTILES)      ! IN Number of tile points.                
*D SFEVAP7A.79,SFEVAP7A.84   
     & FLAND(LAND_FIELD)     ! IN Fraction of gridbox which is land.    
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Coefficient to calculate surface      
!                            !    heat flux into soil.                  
     &,CANOPY(LAND_FIELD,NTILES)                                        
!                            ! IN Surface/canopy water on land          
!                            !    tiles (kg/m2).                        
*D SFEVAP7A.86
     &,FLAKE(LAND_FIELD,NTILES)                                         
!                            ! IN Lake fraction.                        
     &,FRACA(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.89,SFEVAP7A.91   
!                            !    for land tiles.                       
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                            ! IN Lying snow amount on tiles (kg/m2).   
     &,RESFS(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.94,SFEVAP7A.95   
!                            !    of land tiles.                        
     &,RESFT(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.98
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
*D SFEVAP7A.100
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*D SFEVAP7A.103
     &,WT_EXT_TILE(LAND_FIELD,NSHYD,NTILES)                             
*D SFEVAP7A.105
!                            !    extracted from each soil layer        
!                            !    by each tile.                         
*D SFEVAP7A.110
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
*D SFEVAP7A.113
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*D SFEVAP7A.115,SFEVAP7A.119  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SFEVAP7A.126,SFEVAP7A.127  
     &,ECAN_TILE(LAND_FIELD,NTILES)                                     
!                            ! OUT ECAN for land tiles.                 
     &,ELAKE_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT Lake evaporation.                    
*D SFEVAP7A.131,SFEVAP7A.132  
     &,ESOIL_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT ESOIL for land tiles.                
     &,EI_TILE(LAND_FIELD,NTILES)                                       
!                            ! OUT Sublimation from snow or land-ice    
!                            !     (kg per sq m per s).                 
*I SFEVAP7A.136   
*CALL C_0_DG_C                                                          
*D SFEVAP7A.143,SFEVAP7A.145  
     &,E_TILE_OLD(LAND_FIELD,NTILES)                                    
!                            ! Surface moisture flux before adjustment. 
     &,LE_TILE_OLD(LAND_FIELD,NTILES)                                   
!                            ! Surf latent heat flux before adjustment. 
*D SFEVAP7A.152,SFEVAP7A.153  
                                                                        
*D SFEVAP7A.165
      DO N=1,NTILES                                                     
*D SFEVAP7A.168
          E_TILE_OLD(L,N) = FQW_TILE(L,N)                               
          IF (SNOW_TILE(L,N) .GT. 0.) THEN                              
            LE_TILE_OLD(L,N) = (LC + LF)*FQW_TILE(L,N)                  
          ELSE                                                          
            LE_TILE_OLD(L,N) = LC*FQW_TILE(L,N)                         
          ENDIF                                                         
*D SFEVAP7A.172
      DO N=1,NTILES                                                     
*I SFEVAP7A.175   
          ELAKE_TILE(L,N) = 0.                                          
          EI_TILE(L,N) = 0.                                             
*D SFEVAP7A.179,SFEVAP7A.182  
*D SFEVAP7A.184
! Sublimation from snow-covered land tiles                              
*D SFEVAP7A.186,SFEVAP7A.197  
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          IF (SNOW_TILE(L,N) .GT. 0.) THEN                              
            EI_TILE(L,N) =  FQW_TILE(L,N)                               
            EDT = EI_TILE(L,N)*TIMESTEP                                 
            IF ( EDT .GT. SNOW_TILE(L,N) )                              
     &        EI_TILE(L,N) = SNOW_TILE(L,N) / TIMESTEP                  
            FQW_TILE(L,N) = FQW_TILE(L,N) -  EI_TILE(L,N)               
          ENDIF                                                         
        ENDDO                                                           
*D SFEVAP7A.202
*D SFEVAP7A.209
      DO N=1,NTILES                                                     
*D SFEVAP7A.214,SFEVAP7A.215  
            ECAN_TILE(L,N) = (1. - FLAKE(L,N)) *                        
     &                       FRACA(L,N) * FQW_TILE(L,N) / RESFT(L,N)    
            ESOIL_TILE(L,N) = (1. - FLAKE(L,N)) *                       
     &                        (1. - FRACA(L,N))*RESFS(L,N)*FQW_TILE(L,N)
*I SFEVAP7A.216   
            ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N) / RESFT(L,N)     
*D SFEVAP7A.219
              ESOIL_TILE(L,N) =  (1. - FLAKE(L,N)) *                    
     &                           (1. - FRACA(L,N)*CANOPY(L,N)/EDT) *    
*D SFEVAP7A.223,SFEVAP7A.225  
          ELSEIF (SNOW_TILE(L,N).LE.0.) THEN                            
            IF (TSTAR_TILE(L,N).GE.TM) THEN                             
              ECAN_TILE(L,N) = (1. - FLAKE(L,N))*FQW_TILE(L,N)          
              ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N)                
            ELSE                                                        
              EI_TILE(L,N) =  FQW_TILE(L,N)                             
            ENDIF                                                       
*D SFEVAP7A.239
          DO N=1,NTILES                                                 
*D SFEVAP7A.249
          EXT(L,K) = 0.                                                 
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      DO K=1,NSHYD                                                      
        DO N=1,NTILES                                                   
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            EXT(L,K) = EXT(L,K) + TILE_FRAC(L,N)*WT_EXT_TILE(L,K,N)     
     &                                          *ESOIL_TILE(L,N)        
          ENDDO                                                         
*D SFEVAP7A.262,SFEVAP7A.263  
      DO N=1,NTILES                                                     
*D ABX1F405.912,SFEVAP7A.266  
          DIFF_LAT_HTF = (LC + LF)*EI_TILE(L,N) + LC*ECAN_TILE(L,N)     
     &                    + LC*ESOIL_TILE(L,N) + LC*ELAKE_TILE(L,N)     
     &                    - LE_TILE_OLD(L,N)                            
*D SFEVAP7A.268
     &                        ( 1. + ASHTF_TILE(L,N)/(CP*RHOKH_1(L,N)) )
*D SFEVAP7A.270
          DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF_TILE(L,N)   
*D SFEVAP7A.273
          DFQW(L) = DFQW(L) + TILE_FRAC(L,N)*( ECAN_TILE(L,N) +         
     &                  ESOIL_TILE(L,N) + EI_TILE(L,N) + ELAKE_TILE(L,N)
     &                  - E_TILE_OLD(L,N) )                             
*D SFEVAP7A.275,SFEVAP7A.289  
*D SFEVAP7A.298,SFEVAP7A.301  
        FTL_1(I) = FTL_1(I) + FLAND(L)*DFTL(L)                      
        FQW_1(I) = FQW_1(I) + FLAND(L)*DFQW(L)                      
*/-----
*DECLARE SFEXCH7A
*/-----
*D SFEXCH7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFEXCH7A.49,SFEXCH7A.54   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,NELEV, 
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLAND,FLANDG,
     & BQ_1,BT_1,CANHC_TILE,CANOPY,CATCH,DZSOIL,FLAKE,GC,HCONS,         
     & HO2R2_OROG,ICE_FRACT,SNOW_TILE,PSTAR,QW_1,RADNET,RADNET_TILE,    
     & SIL_OROG,SMVCST,TILE_FRAC,TIMESTEP,T_1,Q_1,QCF_1,QCL_1,          
     & TL_1,TI,
     & TS1,snow_depth,nsnow,hcons_snow,snow_depth_1,l_essery_snow,
     & TSTAR_TILE,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,            
     & VFRAC_TILE,VSHR_LAND,VSHR_SSI,ZH,Z0_TILE,Z1_UV,Z1_TQ,LAND_MASK,  
*D SFEXCH7A.56
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,CD,CH,CDR10M,                
*D SFEXCH7A.59
     & RESFS,RESFT,RIB,RIB_TILE,                                        
     & FB_SURF,U_S,Q1_SD,T1_SD,TV1_SD,Z0M_EFF,                          
*D SFEXCH7A.62
     & RHO_CD_MODV1,RHOKH_1,RHOKH_1_SICE,RHOKM_1,RHOKM_LAND,RHOKM_SSI,
     & RHOKPM,RHOKPM_SICE,    
*I SFEXCH7A.63
     & ,T_ELEV,Q_ELEV,QCF_ELEV,QCL_ELEV,P_ELEV
     & ,TL_ELEV,QW_ELEV,BQ_ELEV,BT_ELEV,Z1_ELEV
*D SFEXCH7A.76
     &,NTILES                ! IN Number of land tiles per land point.  
     &,NELEV                 ! IN Number of elevation classes for tiles
*D SFEXCH7A.78
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SFEXCH7A.80
     &,TILE_PTS(NTILES)      ! IN Number of tile points.                
*D SFEXCH7A.87
     &,CANHC_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Areal heat capacity of canopy for     
!                            !    land tiles (J/K/m2).                  
     &,CANOPY(LAND_FIELD,NTILES)                                        
*D SFEXCH7A.90
     &,CATCH(LAND_FIELD,NTILES)                                         
*D SFEXCH7A.92
!                            !    of land tiles (kg/m2).                
*D SFEXCH7A.95
     &,FLAKE(LAND_FIELD,NTILES)                                         
!                            ! IN Lake fraction.                        
     &,GC(LAND_FIELD,NTILES) ! IN "Stomatal" conductance to evaporation 
*D SFEXCH7A.102
     &,FLAND(LAND_FIELD)     ! IN Land fraction on land tiles.          
     &,FLANDG(P_FIELD)       ! IN Land fraction on all tiles.           
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                            ! IN Lying snow on land tiles (kg/m2).     
     &,SNOW_DEPTH(LAND_FIELD,NTILES)                                     
     &,SNOW_DEPTH_1(LAND_FIELD,NTILES)                                     
     &,NSNOW(LAND_FIELD,NTILES)
     &,HCONS_SNOW(LAND_FIELD,NTILES)
     &,HCONS_SURF(LAND_FIELD)
     &,DZ_SURF(LAND_FIELD)
*D SFEXCH7A.104
*D SFEXCH7A.107,SFEXCH7A.110  
     &,RADNET(P_FIELD)       ! IN Sea-ice net surface radiation (W/m2)  
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                            ! IN Land tile net surface radiation (W/m2)
*D SFEXCH7A.115
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*I SFEXCH7A.117   
     &,T_1(P_FIELD)          ! IN Atmospheric temperature (K).          
     &,Q_1(P_FIELD)          ! IN Specific humidity ( kg/kg air).       
     &,QCF_1(P_FIELD)        ! IN Cloud ice (kg per kg air)             
     &,QCL_1(P_FIELD)        ! IN Cloud liquid water (kg                
!                            !    per kg air).                          
*D SFEXCH7A.122
     &,TS1(LAND_FIELD,NTILES)       ! IN Temperature of top soil or land-ice
!
*D SFEXCH7A.124,SFEXCH7A.126  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SFEXCH7A.128,SFEXCH7A.129  
     &,TSTAR_LAND(P_FIELD)   ! IN Land mean surface temperature (K).    
     &,TSTAR_SEA(P_FIELD)    ! IN Open sea surface temperature (K).     
     &,TSTAR_SICE(P_FIELD)   ! IN Sea-ice surface temperature (K).      
     &,TSTAR_SSI(P_FIELD)    ! IN Mean sea surface temperature (K).     
     &,VFRAC_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Fractional canopy coverage for        
!                            !    land tiles.                           
     &,VSHR_LAND(P_FIELD)    ! IN Magnitude of land sfc-to-lowest-level 
*D SFEXCH7A.131
     &,VSHR_SSI(P_FIELD)     ! IN Mag. of mean sea sfc-to-lowest-level  
!                            !    wind shear                            
     &,ZH(P_FIELD)           ! IN Height above surface of top of        
!                            !    boundary layer (metres).              
     &,Z0_TILE(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.133
*I SFEXCH7A.148
     &,L_ESSERY_SNOW
*D SFEXCH7A.159
     & ALPHA1(LAND_FIELD,NTILES)                                        
*D SFEXCH7A.166,SFEXCH7A.168  
!                            !     flux into sea-ice (W/m2/K)           
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT Coefficient to calculate surface heat
!                            !     flux into land tiles (W/m2/K)        
*I SFEXCH7A.170   
     &,CD_SSI(P_FIELD)       ! OUT Bulk transfer coefficient for        
!                            !      momentum over sea mean.             
*D SFEXCH7A.172
     &,CH_SSI(P_FIELD)       ! OUT Bulk transfer coefficient for heat   
!                            !    and/or moisture over sea mean.        
*D SFEXCH7A.179
     &,CHR1P5M(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.190
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.195
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.198
     &,FRACA(LAND_FIELD,NTILES)                                         
*D SFEXCH7A.201
!                            !     for land tiles.                      
*I SFEXCH7A.206   
     &,RESFS(LAND_FIELD,NTILES)                                         
!                            ! OUT Combined soil, stomatal and          
!                            !     aerodynamic resistance factor for    
!                            !     fraction 1-FRACA of land tiles       
     &,RESFT(LAND_FIELD,NTILES)                                         
!                            ! OUT Total resistance factor              
!                            !     FRACA+(1-FRACA)*RESFS for snow-free  
!                            !     tiles, 1 for snow and land-ice.      
     &,RIB(P_FIELD)          ! OUT Mean bulk Richardson number for      
!                            !     lowest layer                         
     &,RIB_TILE(LAND_FIELD,NTILES)                                      
!                            ! OUT RIB for land tiles.                  
     &,FB_SURF(P_FIELD)      ! OUT Surface flux buoyancy over           
!                            !     density (m^2/s^3)                    
     &,U_S(P_FIELD)          ! OUT Surface friction velocity (m/s)      
*D SFEXCH7A.210,SFEXCH7A.221  
*I SFEXCH7A.224   
     &,TV1_SD(P_FIELD)       ! OUT Standard deviation of turbulent      
!                            !     fluctuations of surface layer        
!                            !     virtual temperature (K).             
*D SFEXCH7A.229
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.233
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.238
     &,RHO_ARESIST_TILE(LAND_FIELD,NTILES)                              
*D SFEXCH7A.240
     &,ARESIST_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.242
     &,RESIST_B_TILE(LAND_FIELD,NTILES)                                 
*D SFEXCH7A.249
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.257
     &,RHOKM_LAND(P_FIELD)   ! OUT For land momentum. NB: This is output
!                            !     on UV-grid, but with the first and   
!                            !      last rows set to "missing data".    
     &,RHOKM_SSI(P_FIELD)    ! OUT For mean sea mom. NB: This is output 
!                            !     on UV-grid, but with the first and   
!                            !     last rows set to "missing data".     
     &,RHOKPM(LAND_FIELD,NTILES)                                        
*I SFEXCH7A.271   
*CALL C_EPSLON                                                          
*D SFEXCH7A.284
*D SFEXCH7A.287
*D SFEXCH7A.296
      EXTERNAL SF_OROG,SF_OROG_GB,QSAT,SF_RESIST,TIMER,                 
     &         SFL_INT_LAND,SFL_INT_SEA,                                
*D SFEXCH7A.316
*I SFEXCH7A.320   
     &,Z0M_SEA(P_FIELD)            ! Open sea roughness length for      
!                                  ! momentum transport.                
     &,DB_SEA(P_FIELD)             ! Buoyancy difference for sea points 
     &,V_S_SEA(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for sea points (m/s).              
     &,RECIP_L_MO_SEA(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for sea points (m^-1).      
*I SFEXCH7A.324   
     &,CD_LAND(P_FIELD)            ! Bulk transfer coefficient for      
!                                  !      momentum over land.           
*D SFEXCH7A.333
*I SFEXCH7A.335   
     &,DB_ICE(P_FIELD)             ! Buoyancy difference for sea ice    
     &,V_S_ICE(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for sea ice (m/s).                 
     &,V_S_MIZ(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for marginal sea ice (m/s).        
     &,RECIP_L_MO_ICE(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for sea ice (m^-1).         
     &,RECIP_L_MO_MIZ(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for marginal sea ice (m^-1).
     &,RHO_ARESIST_LAND(P_FIELD)   ! Land mean of rho_aresist_tile      
!
*D SFEXCH7A.342
     & CD_STD(LAND_FIELD,NTILES)   ! Local drag coefficient for calc    
*D SFEXCH7A.344,SFEXCH7A.345  
     &,CD_TILE(LAND_FIELD,NTILES)  ! Drag coefficient                   
     &,CH_TILE(LAND_FIELD,NTILES)  ! Transfer coefficient for heat and  
*I SFEXCH7A.350   
     &,FZ0(LAND_FIELD)             ! Aggregation function for Z0.       
*D SFEXCH7A.352,SFEXCH7A.353  
     &,QSTAR_TILE(LAND_FIELD,NTILES)!Surface saturated sp humidity.     
     &,RHOKM_1_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.355
     &,WIND_PROFILE_FACTOR(LAND_FIELD,NTILES)                           
*D SFEXCH7A.359,SFEXCH7A.360  
     &,Z0_GB(LAND_FIELD)           ! GBM roughness length               
     &,Z0M_EFF_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.362
     &,Z0F_TILE(LAND_FIELD,NTILES) !Roughness length for free convective
*I SFEXCH7A.363   
     &,DB_TILE(LAND_FIELD,NTILES)  ! Buoyancy difference for surface    
!                                  ! tile                               
     &,V_S_TILE(LAND_FIELD,NTILES) ! Surface layer scaling velocity     
!                                  ! for tiles (m/s).                   
     &,V_S_STD(LAND_FIELD,NTILES)  ! Surface layer scaling velocity     
!                                  ! for tiles excluding orographic     
!                                  ! form drag (m/s).                   
     &,RECIP_L_MO_TILE(LAND_FIELD,NTILES)                               
!                                  ! Reciprocal of the Monin-Obukhov    
!                                  ! length for tiles (m^-1).           
*I SFEXCH7A.371
     &,N_1,K
*I SFEXCH7A.377
      REAL
     & T_ELEV(LAND_FIELD,NELEV)
     &,Q_ELEV(LAND_FIELD,NELEV)
     &,QCF_ELEV(LAND_FIELD,NELEV)
     &,QCL_ELEV(LAND_FIELD,NELEV)
     &,P_ELEV(LAND_FIELD,NELEV)
     &,TL_ELEV(LAND_FIELD,NELEV)
     &,QW_ELEV(LAND_FIELD,NELEV)
     &,BQ_ELEV(LAND_FIELD,NELEV)
     &,BT_ELEV(LAND_FIELD,NELEV)
     &,Z1_ELEV(LAND_FIELD,NELEV)

      REAL
     & QS_ELEV(LAND_FIELD,NELEV)
     &,FB_SURF_T,TV1_SD_T
*D ABX1F405.902
      DO N=1,NTILES                                                     
*D SFEXCH7A.389
        IF ( ICE_FRACT(I).GT.0.0 .AND. FLANDG(I).LT.1.0 ) THEN         
*D SFEXCH7A.400,SFEXCH7A.402  
        RHOSTAR(I) = PSTAR(I) / 
     &    ( R*(FLANDG(I)*TSTAR_LAND(I) +                            
     &    (1.-FLANDG(I))*TSTAR_SSI(I)) )                            
*D SFEXCH7A.404,SFEXCH7A.408  
*D SFEXCH7A.412
      CALL QSAT(QSTAR_ICE(P1),TSTAR_SICE(P1),PSTAR(P1),P_POINTS)        
*D SFEXCH7A.417
      DO K=1,NELEV                                                     
        CALL QSAT(QS_ELEV(LAND1,K),TL_ELEV(LAND1,K),
     &            P_ELEV(LAND1,K),LAND_PTS)
      ENDDO
      DO N=1,NTILES                                                     
        k=mod(n-1,nelev)+1
*D SFEXCH7A.419
     &            P_ELEV(LAND1,K),LAND_PTS)
*I SFEXCH7A.437   
        DB_SEA(I) = 0.                                                  
        DB_ICE(I) = 0.                                                  
*D SFEXCH7A.443
      DO N=1,NTILES                                                     
        N_1=(n-1)/nelev + 1
*D SFEXCH7A.446,SFEXCH7A.450  
          IF ( SNOW_TILE(L,N).GT.0. .AND. 
     &         N.le.(ntiles-nelev)) THEN
            Z0 = Z0_TILE(L,N) - 4.0E-4*SNOW_TILE(L,N)                   
            if (L_ESSERY_SNOW) Z0 = Z0_TILE(L,N) - 0.1*SNOW_DEPTH(L,N)
            ZETA1 = MIN( 5.0E-4 , Z0_TILE(L,N)  )                       
*I SFEXCH7A.451   
          ELSE                                                          
            Z0M_TILE(L,N) = Z0_TILE(L,N)                                
*D ABX1F405.909,ABX1F405.910
           Z0H_TILE(L,N) = Z0H_Z0M(N_1)*Z0M_TILE(L,N)
           Z0F_TILE(L,N) = Z0H_Z0M(N_1)*Z0M_TILE(L,N)
*I SFEXCH7A.455   
          DB_TILE(L,N) = 0.                                             
*D SFEXCH7A.459
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.463
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG,Z0M_TILE(1,N),Z1_ELEV(1,K),
*D SFEXCH7A.470
! of Richardson number. RESFT=1 for snow and land-ice.                  
*D SFEXCH7A.473,SFEXCH7A.474  
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.478,SFEXCH7A.479
          ZETAM = LOG ( (Z1_ELEV(L,K) + Z0M_TILE(L,N))/Z0M_TILE(L,N) )
          ZETAH = LOG ( (Z1_ELEV(L,K) + Z0M_TILE(L,N))/Z0H_TILE(L,N) )
*D SFEXCH7A.481
          DQ(L) = QW_ELEV(L,K) - QSTAR_TILE(L,N)
*D SFEXCH7A.486,SFEXCH7A.487  
     &   CANOPY(1,N),CATCH(1,N),CHN,DQ,EPDT,FLAKE(1,N),GC(1,N),         
     &   SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),RESFS(1,N),RESFT(1,N),     
     &   LTIMER                                                         
*D SFEXCH7A.489,SFEXCH7A.494  
*D SFEXCH7A.503,SFEXCH7A.506  
     & P_POINTS,P_FIELD,P1,FLANDG,NSICE,SICE_INDEX,                  
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    
     & TSTAR_SEA,VSHR_SSI,Z0_ICE,Z0H_SEA,Z0_ICE,Z0MSEA,Z1_TQ,Z1_UV,     
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             
*D SFEXCH7A.510
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.513,SFEXCH7A.515  
     &   BQ_ELEV(1,K),BT_ELEV(1,K),QSTAR_TILE(1,N),
     &   QW_ELEV(1,K),RESFT(1,N),TL_ELEV(1,K),
     &   TSTAR_TILE(1,N),VSHR_LAND,Z0H_TILE(1,N),Z0M_TILE(1,N),
     &   Z1_ELEV(1,K),Z1_ELEV(1,K),
     &   RIB_TILE(1,N),DB_TILE(1,N),LTIMER
*D SFEXCH7A.524
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.528
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG,Z0M_TILE(1,N),Z1_ELEV(1,K),
*D SFEXCH7A.538,SFEXCH7A.540  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_ICE,DB_ICE,VSHR_SSI,Z0_ICE,Z0_ICE,Z0_ICE,      
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_ICE,CH_ICE,V_S_ICE,                             
     &               RECIP_L_MO_ICE,LTIMER)                             
*D SFEXCH7A.543,SFEXCH7A.545  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_ICE,DB_ICE,VSHR_SSI,Z0_MIZ,Z0_MIZ,Z0_MIZ,      
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_MIZ,CH_MIZ,V_S_MIZ,                             
     &               RECIP_L_MO_MIZ,LTIMER)                             
*D SFEXCH7A.548,SFEXCH7A.550  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_SEA,DB_SEA,VSHR_SSI,Z0MSEA,Z0H_SEA,Z0F_SEA,    
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_SEA,CH_SEA,V_S_SEA,                             
     &               RECIP_L_MO_SEA,LTIMER)                             
*D SFEXCH7A.553
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.556,SFEXCH7A.558  
     &   RIB_TILE(1,N),DB_TILE(1,N),VSHR_LAND,                          
     &   Z0M_EFF_TILE(1,N),Z0H_TILE(1,N),Z0F_TILE(1,N),ZH,              
     &   Z1_ELEV(1,K),Z1_ELEV(1,K),WIND_PROFILE_FACTOR(1,N),
     &   CD_TILE(1,N),CH_TILE(1,N),CD_STD(1,N),                         
     &   V_S_TILE(1,N),V_S_STD(1,N),RECIP_L_MO_TILE(1,N),LTIMER         
*D SFEXCH7A.563,SFEXCH7A.564  
!!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles       
*D SFEXCH7A.567
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.571,SFEXCH7A.572
          DQ(L) = QW_ELEV(L,K) - QSTAR_TILE(L,N)
          EPDT(L) = - RHOSTAR(I)*CH_TILE(L,N)*VSHR_LAND(I)
     &      *DQ(L)*TIMESTEP    
*D SFEXCH7A.576,SFEXCH7A.578  
     &   CANOPY(1,N),CATCH(1,N),CH_TILE(1,N),DQ,EPDT,FLAKE(1,N),        
     &   GC(1,N),SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),                   
     &   RESFS(1,N),RESFT(1,N),                                         
     &   LTIMER)                                                        
*I SFEXCH7A.585   
        CD_LAND(I) = 0.                                               
        CD_SSI(I) = 0.                                                
        CH_SSI(I) = 0.                                                
*D SFEXCH7A.590
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D SFEXCH7A.592,SFEXCH7A.595  
            CD_SSI(I) = ( ICE_FRACT(I)*CD_MIZ(I) +                
     &              (0.7-ICE_FRACT(I))*CD_SEA(I) ) / 0.7  ! P2430.5 
            CH_SSI(I) = ( ICE_FRACT(I)*CH_MIZ(I) +                
     &              (0.7-ICE_FRACT(I))*CH_SEA(I) ) / 0.7  ! P2430.4 
*D SFEXCH7A.597,SFEXCH7A.600  
            CD_SSI(I) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +          
     &              (ICE_FRACT(I)-0.7)*CD_ICE(I) ) / 0.3  ! P2430.7 
            CH_SSI(I) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +          
     &              (ICE_FRACT(I)-0.7)*CH_ICE(I) ) / 0.3  ! P2430.7 
*I SFEXCH7A.601   
          CD(I)=(1.-FLANDG(I))*CD_SSI(I)                          
          CH(I)=(1.-FLANDG(I))*CH_SSI(I)                          
*D SFEXCH7A.607,SFEXCH7A.608  
      DO N=1,NTILES                                                     
       DO J=1,TILE_PTS(N)                                              
*D SFEXCH7A.611,SFEXCH7A.612  
          CD_LAND(I) = CD_LAND(I) + TILE_FRAC(L,N)*CD_TILE(L,N)     
          CD(I) = CD(I) + FLANDG(I)* TILE_FRAC(L,N)*CD_TILE(L,N)  
          CH(I) = CH(I) + FLANDG(I)*TILE_FRAC(L,N)*CH_TILE(L,N) 
*I SFEXCH7A.625   
        RHO_ARESIST_LAND(I)=0.0                                      
*D SFEXCH7A.628
        RHOKM_LAND(I) = 0.                                            
        RHOKM_SSI(I) = 0.                                             
*D SFEXCH7A.631,SFEXCH7A.636  
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
          RHOKM_SSI(I) = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)   ! P243.124  
          RHOKH_1_SICE(I) = RHOSTAR(I)*CH_SSI(I)*VSHR_SSI(I)    
!                                                           ! P243.125  
          RHO_ARESIST(I) = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)     
          ARESIST(I) =  1. / (CD_SSI(I) * VSHR_SSI(I))            
          RESIST_B(I)= (CD_SSI(I)/CH_SSI(I) - 1.0) * ARESIST(I) 
*D SFEXCH7A.642
      DO N=1,NTILES                                                     
*D SFEXCH7A.651,SFEXCH7A.655  
          RHOKM_1_TILE(L,N) = RHOSTAR(I)*CD_TILE(L,N)*VSHR_LAND(I)  
!                                                         ! P243.124    
          RHOKM_LAND(I) = RHOKM_LAND(I) +                           
     &           TILE_FRAC(L,N)*RHOKM_1_TILE(L,N)                       
          RHOKH_1(L,N) = RHOSTAR(I)*CH_TILE(L,N)*VSHR_LAND(I)       
!                                                         ! P243.125    
          RHO_ARESIST_TILE(L,N) = RHOSTAR(I) * CD_STD(L,N)            
     &                * VSHR_LAND(I)                                  
          ARESIST_TILE(L,N) = 1. / ( CD_STD(L,N) * VSHR_LAND(I) )     
*I SFEXCH7A.657   
          IF (RESIST_B_TILE(L,N) .LT. 0.) RESIST_B_TILE(L,N) = 0.       
          RHO_ARESIST_LAND(I) = RHO_ARESIST_LAND(I) +               
     &                     TILE_FRAC(L,N)*RHO_ARESIST_TILE(L,N)         
*I SFEXCH7A.660   
      DO L=LAND1,LAND1+LAND_PTS-1                                       
        I = LAND_INDEX(L)                                               
        RHO_ARESIST(I) = FLANDG(I)*RHO_ARESIST_LAND(I) +          
     &                      (1.0-FLANDG(I))*RHO_ARESIST(I)          
        ARESIST(I) = RHOSTAR(I) / RHO_ARESIST(I)                        
      ENDDO                                                             
                                                                        
*I SFEXCH7A.661   
        RHOKM_1(I)= FLANDG(I) * RHOKM_LAND(I) +                   
     &                     (1.0-FLANDG(I)) * RHOKM_SSI(I)           
*D SFEXCH7A.667,SFEXCH7A.668  
!!  moisture.                                                           
*D SFEXCH7A.674
*D SFEXCH7A.677
! Sea and sea-ice                                                       
      CALL SF_FLUX_SEA (                                                
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,FLANDG,                  
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1_SICE,      
     & TI,TL_1,TSTAR_SICE,TSTAR_SEA,Z0_ICE,Z0_ICE,Z0H_SEA,Z0MSEA,Z1_TQ, 
     & ALPHA1_SICE,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,       
     & RHOKPM_SICE,LTIMER                                               
     & )                                                                
                                                                        
! Land tiles                                                            
      DO N=1,NTILES                                                     
*D SFEXCH7A.684,SFEXCH7A.694  
      DO N=1,NTILES                                                     
        k=mod(n-1,nelev)+1
        DO L = 1,LAND_FIELD
            hcons_surf(l)=hcons(l)
            if (n.gt.(ntiles-nelev)) hcons_surf(l)=0.24 !!ice_hcon
            dz_surf(l)=dzsoil
           if (l_essery_snow) then
             if (nsnow(l,n).gt.0) then
               hcons_surf(l)=hcons_snow(l,n)
               dz_surf(l)=snow_depth_1(l,n)
             endif
          endif
        ENDDO
        CALL SF_FLUX_LAND (N,                                             
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),FLAND,
     &   LAND_INDEX,TILE_INDEX(1,N),     
     &   CANHC_TILE(1,N),DZ_SURF,HCONS_SURF,
     &   QS_ELEV(1,K),QSTAR_TILE(1,N),QW_ELEV(1,K),
     &   RADNET_TILE(1,N),RESFT(1,N),RHOKH_1(1,N),SMVCST,SNOW_TILE(1,N),
     &   TILE_FRAC(1,N),TIMESTEP,TL_ELEV(1,K),TS1(1,N),TSTAR_TILE(1,N),
     &   VFRAC_TILE(1,N),Z0H_TILE(1,N),Z0M_EFF_TILE(1,N),Z1_ELEV(1,K), 
     &   FQW_1,FTL_1,                                                   
     &   ALPHA1(1,N),ASHTF_TILE(1,N),FQW_TILE(1,N),FTL_TILE(1,N),       
     &   RHOKPM(1,N),snow_depth(1,n),nsnow(1,n),l_essery_snow,
     &   LTIMER
     & )                                                                
*D SFEXCH7A.696,SFEXCH7A.727  
*D SFEXCH7A.742,SFEXCH7A.743  
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_SSI,RHOSTAR,VSHR_SSI,      
*D SFEXCH7A.749
      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
*D SFEXCH7A.751,SFEXCH7A.753
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),
     &   LAND_INDEX,TILE_INDEX(1,N),FLAND,     
     &   BQ_ELEV(1,k),BT_ELEV(1,k),
     &   FQW_TILE(1,N),FTL_TILE(1,N),RHOKM_1_TILE(1,N),
     &   RHOSTAR,VSHR_LAND,Z0M_TILE(1,N),Z1_ELEV(1,k),TILE_FRAC(1,N),
*I SFEXCH7A.756   
                                                                        
!-----------------------------------------------------------------------
!! Calculate scaling parameters required for new boundary layer scheme  
!-----------------------------------------------------------------------
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        U_S(I) = SQRT(FLANDG(I)*CD_LAND(I)*VSHR_LAND(I)         
     &    +(1.-FLANDG(I))*CD_SSI(I)*VSHR_SSI(I))                  
        FB_SURF(I) = G * ( BT_1(I)*FTL_1(I) +                           
     &                     BQ_1(I)*FQW_1(I) ) / RHOSTAR(I)              
        IF (FB_SURF(I) .GT. 0.0) THEN                                   
          TV1_SD(I) = T_1(I) *                                          
     &                ( 1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I) ) *
     &                ( BT_1(I)*T1_SD(I) + BQ_1(I)*Q1_SD(I) )           
          IF (TV1_SD(I) .LT. 0.0) THEN                                  
            TV1_SD(I)=0.0                                               
          ENDIF                                                         
        ELSE                                                            
            TV1_SD(I)=0.0                                               
        ENDIF                                                           
      ENDDO                                                             

!! This is now a bit inconsistent, given the use of elevations - 
!! eg BT_1 is not nec. eq mean(BT_ELEV), but FTL_1 has been correctly
!! calc'd from FTL_TILE (incl. BT_ELEV etc.) so adjust
      DO N=1,NTILES                                                     
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        FB_SURF(I)=0.
        TV1_SD(I)=0.0
      ENDDO
      ENDDO

      DO N=1,NTILES                                                     
        K=mod(n-1,nelev)+1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          FB_SURF_T = G * ( BT_ELEV(L,K)*FTL_TILE(L,N) +
     &                     BQ_ELEV(L,K)*FQW_TILE(L,N) ) / RHOSTAR(I)
          IF (FB_SURF_T .GT. 0.0) THEN
            TV1_SD_T = T_ELEV(L,K) *
     &  ( 1.0 + C_VIRTUAL*Q_ELEV(L,K) - QCL_ELEV(L,K) - QCF_ELEV(L,K) )*
     &                ( BT_ELEV(L,K)*T1_SD(I) + BQ_ELEV(L,K)*Q1_SD(I) )
            IF (TV1_SD_T .LT. 0.0) THEN
              TV1_SD_T=0.0
            ENDIF 
          ELSE
            TV1_SD_T=0.0
          ENDIF

          TV1_SD(I)=TV1_SD(I)+(FLAND(L)*TILE_FRAC(L,N)*TV1_SD_T)
          FB_SURF(I)=FB_SURF(I)+(FLAND(L)*TILE_FRAC(L,N)*FB_SURF_T)

        ENDDO
      ENDDO
                                                                        
*D SFEXCH7A.776,SFEXCH7A.777  
        IF (FLANDG(I).LT.1.0) THEN                                    
          TAU = RHOKM_SSI(I) * VSHR_SSI(I)             ! P243.130   
*D SFEXCH7A.779
     &      TAU = RHOSTAR(I) * CD_SEA(I) * VSHR_SSI(I) * VSHR_SSI(I)  
*D SFEXCH7A.801,SFEXCH7A.802  
        IF ( .NOT.LAND_MASK(I) )H_BLEND_OROG(I) = H_BLEND_MIN       
        IF (FLANDG(I).LT.1.0) THEN                                    
*D SFEXCH7A.817
      DO L = LAND1,LAND1+LAND_PTS-1                                     
        Z0_GB(L) = 0.                                                   
        FZ0(L) = 0.
      ENDDO
                                                                        
! Copy sea/sea-ice roughness lengths and Richardson numbers onto the    
! the land point tiles.                                                 
                                                                        
! Weight so that gridbox mean can be calculated:                        
      DO L=LAND1,LAND1+LAND_PTS-1          
        I = LAND_INDEX(L)                         
        IF(FLAND(L).LT.1.0)THEN                                         
          RIB(I) = (1.0-FLAND(L))*RIB(I)                            
          FZ0(L) = (1.0-FLAND(L)) / (LOG(LB/Z0M(I))**2)               
        ENDIF                                                           
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
*D SFEXCH7A.821
          RIB(I) = RIB(I) + FLAND(L)*TILE_FRAC(L,N)*RIB_TILE(L,N)  
          FZ0(L) = FZ0(L) 
     &      + FLAND(L)*TILE_FRAC(L,N) / (LOG(LB/Z0M_TILE(L,N))**2)
*D SFEXCH7A.824,SFEXCH7A.825  
*D SFEXCH7A.827,SFEXCH7A.833  
        Z0_GB(L) = LB * EXP( - SQRT(1./FZ0(L)) )                        
*D SFEXCH7A.852
          IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).LE.0. )           
     &      Z0M_SEA(I) = Z0MSEA(I)  
          IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0. ) THEN
*I SFEXCH7A.856   
            Z0M_SEA(I) = Z0_ICE(I)                                  
            V_S_SEA(I) = V_S_ICE(I)                                     
            RECIP_L_MO_SEA(I) = RECIP_L_MO_ICE(I)                       
*D SFEXCH7A.861,SFEXCH7A.864  
     &   P_POINTS,P_FIELD,P1,FLANDG,                                 
     &   VSHR_SSI,CD_SEA,CH_SEA,RIB,Z0M_SEA,Z0H_SEA,Z0F_SEA,Z1_UV,      
     &   RECIP_L_MO_SEA,V_S_SEA,                                        
     &   SU10,SV10,ST1P5,SQ1P5,                                         
     &   CDR10M,CHR1P5M_SICE,LTIMER                                     
*D SFEXCH7A.868
        DO N=1,NTILES                                                   
          k=mod(n-1,nelev)+1
*D SFEXCH7A.870,SFEXCH7A.875  
     &     P_FIELD,LAND_FIELD,TILE_PTS(N),
     &     TILE_INDEX(1,N),LAND_INDEX,FLANDG,   
     &     VSHR_LAND,CD_STD(1,N),CD_TILE(1,N),CH_TILE(1,N),
     &     RIB_TILE(1,N),    
     &     TILE_FRAC(1,N),WIND_PROFILE_FACTOR(1,N),                     
     &     Z0M_EFF_TILE(1,N),Z0M_TILE(1,N),Z0H_TILE(1,N),               
     &     Z0F_TILE(1,N),Z1_ELEV(1,k),RECIP_L_MO_TILE(1,N),
     &     V_S_TILE(1,N),V_S_STD(1,N),                                  
     &     SU10,SV10,ST1P5,SQ1P5,                                       
     &     CDR10M,CHR1P5M(1,N),LTIMER                                   
*/-----
*DECLARE SFFLUX7A
*/-----
*D SFFLUX7A.34
      SUBROUTINE SF_FLUX_LAND (TYP,
*D SFFLUX7A.35,SFFLUX7A.37   
     & P_FIELD,LAND_FIELD,TILE_PTS,FLAND,LAND_INDEX,TILE_INDEX,         
     & CANHC,DZ_SURF,HCONS,QS1,QSTAR,QW_1,RADNET,RESFT,RHOKH_1,SMVCST,   
     & SNOW,TILE_FRAC,TIMESTEP,TL_1,TS1,TSTAR,VFRAC,Z0H,Z0M_EFF,Z1_TQ,  
*D SFFLUX7A.39
     & ALPHA1,ASHTF,FQW_1,FTL_1,RHOKPM,
     & snow_depth,nsnow,l_essery_snow,
     & LTIMER
*I SFFLUX7A.49
     &,TYP
*I SFFLUX7A.52
     & ,L_essery_snow              ! IN Logical for new snow scheme
     & ,OLDDEEPSNOW
*D SFFLUX7A.55,SFFLUX7A.58   
     & FLAND(LAND_FIELD)                                                
     &,CANHC(LAND_FIELD)   ! IN Areal heat capacity of canopy (J/K/m2). 
     &,DZ_SURF(LAND_FIELD) ! IN Soil or land-ice surface layer          
!                          !    thickness (m).                          
     &,HCONS(LAND_FIELD)   ! IN Soil thermal conductivity (W/m/K).      
     &,QS1(LAND_FIELD)        ! IN Sat. specific humidity
*D SFFLUX7A.60
     &,QW_1(LAND_FIELD)       ! IN Total water content of lowest
*D SFFLUX7A.62
     &,RADNET(LAND_FIELD)  ! IN Net surface radiation (W/m2) positive   
*I SFFLUX7A.65    
     &,SMVCST(LAND_FIELD)  ! IN Volumetric saturation point             
!                          !    - zero at land-ice points.              
     &,SNOW(LAND_FIELD)    ! IN Lying snow amount (kg/m2).              
     &,SNOW_DEPTH(LAND_FIELD)    ! IN Lying snow amount (m)
     &,NSNOW(LAND_FIELD)
*I SFFLUX7A.67    
     &,TIMESTEP            ! IN Timestep (s).                           
*D SFFLUX7A.68    
     &,TL_1(LAND_FIELD)       ! IN Liquid/frozen water temperature for
*I SFFLUX7A.71    
     &,VFRAC(LAND_FIELD)   ! IN Fractional canopy coverage.           
*D SFFLUX7A.74
     &,Z1_TQ(LAND_FIELD)      ! IN Height of lowest atmospheric level (m).
*D SFFLUX7A.81
     & ASHTF(LAND_FIELD)   ! OUT Coefficient to calculate surface       
!                          !     heat flux into soil (W/m2/K).          
     &,ALPHA1(LAND_FIELD)  ! OUT Gradient of saturated specific humidity
*I SFFLUX7A.92    
*CALL CSIGMA                                                            
*CALL C_SOILH                                                           
*D SFFLUX7A.94
*I SFFLUX7A.109   
     &,DS_RATIO            ! 2 * snowdepth / depth of top soil layer.   
*I SFFLUX7A.110   
     &,LH                  ! Latent heat (J/K/kg).                      
*D SFFLUX7A.126,SFFLUX7A.134
        D_T = TSTAR(L) - TL_1(L)
        IF (D_T .GT. 0.05 .OR. D_T .LT. -0.05) THEN
          ALPHA1(L) = (QSTAR(L) - QS1(L)) / D_T
        ELSEIF (TL_1(L) .GT. TM) THEN
          ALPHA1(L) = EPSILON*LC*QS1(L)*(1. + C_VIRTUAL*QS1(L))
     &                                            / ( R*TL_1(L)*TL_1(L))
        ELSE
          ALPHA1(L) = EPSILON*LS*QS1(L)*(1. + C_VIRTUAL*QS1(L))
     &                                            / ( R*TL_1(L)*TL_1(L))
*I SFFLUX7A.139   
        ASHTF(L) = 2.0 * HCONS(L) / DZ_SURF(L)                              


        !IF (SNOW(L).GT.0.0 .AND. SMVCST(L).NE.0.) THEN                  
        !valid over ice as well now ice_htc is ice, not snow
        !assumes DZSNOW(1)=DZSOIL(1)=appropriate interaction depth
        OLDDEEPSNOW=.FALSE.
        IF (SNOW(L).GT.0.0) OLDDEEPSNOW=.TRUE.
        IF (L_ESSERY_SNOW) THEN
          IF (NSNOW(l).GT.0) OLDDEEPSNOW=.FALSE.
        ENDIF

        IF (OLDDEEPSNOW) THEN

          if (l_essery_snow) then
            DS_RATIO = 2.0 *snow_depth(L)/DZ_SURF(L)
          else
            DS_RATIO = 2.0 * SNOW(L) / (RHO_SNOW * DZ_SURF(L))
         endif
          IF (DS_RATIO.LE.1.0) THEN                                     
            ASHTF(L) =  ASHTF(L) /                                      
     &                         (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
          ELSE                                                          
            ASHTF(L) =  ASHTF(L)*SNOW_HCON / HCONS(L)                   
          ENDIF                                                         
        ENDIF                                                           
        ASHTF(L) = (1. - VFRAC(L))*ASHTF(L) + CANHC(L)/TIMESTEP +       
     &             4*(1. + VFRAC(L))*SBCON*TS1(L)**3                    
      ENDDO                                                             
                                                                        
      DO J=1,TILE_PTS                                                   
        L = TILE_INDEX(J)                                               
*D SFFLUX7A.142,SFFLUX7A.153
        LH = LC                                                         
        IF (SNOW(L) .GT. 0.) LH = LS                                    
        RHOKPM(L) = RHOKH_1(L) / ( ASHTF(L)  +                          
     &                         RHOKH_1(L)*(LH*ALPHA1(L)*RESFT(L) + CP) )
        RAD_REDUC = RADNET(L) - ASHTF(L) * ( TL_1(L) - TS1(L)           
     &                         + GRCP*(Z1_TQ(L) + Z0M_EFF(L) - Z0H(L)) )
     &                        + CANHC(L)*(TSTAR(L) - TS1(L)) / TIMESTEP 
        DQ1 = QS1(L) - QW_1(L) +
     &                   GRCP*ALPHA1(L)*(Z1_TQ(L) + Z0M_EFF(L) - Z0H(L))
        FQW_1(L) = RESFT(L)*RHOKPM(L)*( ALPHA1(L)*RAD_REDUC
     &                                + (CP*RHOKH_1(L) + ASHTF(L))*DQ1 )
        FTL_1(L) = RHOKPM(L)*(RAD_REDUC - LH*RESFT(L)*RHOKH_1(L)*DQ1)

        FTL_1_GB(I) = FTL_1_GB(I) + FLAND(L)*TILE_FRAC(L)*FTL_1(L)      
        FQW_1_GB(I) = FQW_1_GB(I) + FLAND(L)*TILE_FRAC(L)*FQW_1(L)
*D SFFLUX7A.171,SFFLUX7A.175  
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,FLANDG,                  
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1,TI,        
     & TL_1,TSTAR_SICE,TSTAR_SEA,Z0H_ICE,Z0M_ICE,Z0H_SEA,Z0M_SEA,Z1_TQ, 
     & ALPHA1,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,RHOKPM,     
     & LTIMER)                                                          
*D SFFLUX7A.188
*D SFFLUX7A.191,SFFLUX7A.192  
     & FLANDG(P_FIELD)                                          
*D SFFLUX7A.205
     &,TSTAR_SICE(P_FIELD)  ! IN Sea-ice surface temperature (K).       
*I SFFLUX7A.218   
     &,ASHTF(P_FIELD)      ! OUT Coefficient to calculate surface       
!                          !     heat flux into sea-ice (W/m2/K).       
*I SFFLUX7A.231   
*CALL C_KAPPAI                                                          
*I SFFLUX7A.233   
*CALL CSIGMA                                                            
*D SFFLUX7A.273
        D_T = TSTAR_SICE(I) - TL_1(I)                                   
*I SFFLUX7A.282   
        ASHTF(I) = 2 * KAPPAI / DE + 4*SBCON*TI(I)**3                   
*D SFFLUX7A.286
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
*D SFFLUX7A.310,SFFLUX7A.311  
          FTL_1(I) = (1.-FLANDG(I))*(FTL_ICE(I) + H_SEA(I) / CP)        
          FQW_1(I) = (1.-FLANDG(I))*(FQW_ICE(I) + E_SEA(I))             
*/-----
*DECLARE SFLINT6A
*/-----
*D SFLINT6A.2
*IF DEF,A03_8A                                                          
*D SFLINT6A.21
!!!  SUBROUTINES SFL_INT_SEA AND SFL_INT_LAND--------------------------
*D SFLINT6A.45,SFLINT6A.50   
      SUBROUTINE SFL_INT_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK
     &,VSHR,CD,CH,RIB,Z0M,Z0H,Z0F,Z1
     &,RECIP_L_MO,V_S
     +,SU10,SV10,ST1P5,SQ1P5
     &,CDR10M,CHR1P5M,LTIMER
*D ARN0F405.1829,ARN0F405.1830 
     &,VSHR(P_FIELD)     ! IN Wind speed difference between the
!                        !    surface and the lowest wind level in
!                        !    the atmosphere (m/s).
*D ARN0F405.1831,SFLINT6A.67   
*D SFLINT6A.72
                                                                        
      LOGICAL                                                           
     & LAND_MASK(P_FIELD)        ! IN T for land points, F otherwise.
     &,SU10                      ! IN 10m U-wind diagnostic flag        
     &,SV10                      ! IN 10m V-wind diagnostic flag        
     &,ST1P5                     ! IN screen temp diagnostic flag       
     &,SQ1P5                     ! IN screen specific humidity          
!                                !    diagnostic flag                   
     +,LTIMER                    ! IN TIMER diagnostics flag            
! Output variables                                                      
!                                                                       
      REAL                                                              
     + CDR10M(P_FIELD)   ! OUT interpolation coefficicent for 10m wind  
     +,CHR1P5M(P_FIELD)  ! OUT Interpolation coefficient for 1.5m       
!                        !     temperature                              

      REAL
     & RIB(P_FIELD)    ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(P_FIELD)    ! DUMMY Used in 7A boundary layer scheme
     &,Z1(P_FIELD)     ! DUMMY Used in 7A boundary layer scheme

!*                                                                      
!*L---------------------------------------------------------------------
      EXTERNAL TIMER , PHI_M_H                                          
!*                                                                      
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-                              
*CALL C_VKMAN                                                           
      REAL Z_OBS_TQ,Z_OBS_WIND                                          
      PARAMETER (                                                       
     + Z_OBS_TQ = 1.5    ! Height of screen observations of temperature 
!                        ! and humidity.                                
     +,Z_OBS_WIND = 10.0 ! Height of surface wind observations.         
     +)                                                                 
!                                                                       
!  Define local storage.                                                
!                                                                       
!  (a) Local work arrays.                                               
!                                                                       
      REAL                                                              
     & Z_WIND(P_FIELD)     ! Height of wind observations.               
     &,Z_TEMP(P_FIELD)     ! Height of temperature and humidity         
!                          ! observations.                              
     &,PHI_M_OBS(P_FIELD)  ! Monin-Obukhov stability function for       
!                          ! momentum integrated to the wind observation
!                          ! height.                                    
     &,PHI_H_OBS(P_FIELD)  ! Monin-Obukhov stability function for       
!                          ! scalars integrated to their observation    
!                          ! height.                                    
!                                                                       
!  (b) Scalars.                                                         
!                                                                       
      INTEGER                                                           
     + I       ! Loop counter (horizontal field index).                 
!*                                                                      
      IF (LTIMER) THEN                                                  
        CALL TIMER('SFL_INT   ',3)                                      
      ENDIF                                                             
!                                                                       
!-----------------------------------------------------------------------
!! 1. If diagnostics required calculate M-O stability functions at      
!!    observation heights.                                              
!-----------------------------------------------------------------------
                                                                        
      IF (SU10 .OR. SV10 .OR. ST1P5 .OR. SQ1P5) THEN                    
        DO I=P1,P1+P_POINTS-1                                           
          Z_WIND(I) = Z_OBS_WIND                                        
          Z_TEMP(I) = Z_OBS_TQ + Z0H(I) - Z0M(I)
        ENDDO                                                           
        CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
     &                    RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,
     &                    PHI_M_OBS,PHI_H_OBS,LTIMER)
      ENDIF                                                             
                                                                        
!-----------------------------------------------------------------------
!! 2. If diagnostics required calculate interpolation coefficient       
!!    for 1.5m screen temperature and specific humidity.                
!-----------------------------------------------------------------------
!                                                                       
      IF (ST1P5 .OR. SQ1P5) THEN                                        
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CHR1P5M(I) = CH(I) * VSHR(I) * PHI_H_OBS(I)/(VKMAN*V_S(I))
          ENDIF
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
!-----------------------------------------------------------------------
!! 3. If diagnostics required calculate interpolation coefficient       
!!    for 10m winds.                                                    
!-----------------------------------------------------------------------
!                                                                       
      IF ( SU10 .OR. SV10 ) THEN
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CDR10M(I) = CD(I) * VSHR(I) * PHI_M_OBS(I)/(VKMAN*V_S(I))
          ENDIF
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
      IF (LTIMER) THEN                                                  
        CALL TIMER('SFL_INT ',4)                                        
      ENDIF                                                             
      RETURN                                                            
      END                                                               

!!!                                                                     
!!!---------------------------------------------------------------------
!*L  Arguments :-                                                       
      SUBROUTINE SFL_INT_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX
     &,VSHR,CD_STD,CD,CH,RIB,TILE_FRAC,WIND_PROFILE_FACTOR
     &,Z0M,Z0M_STD,Z0H,Z0F,Z1
     &,RECIP_L_MO,V_S,V_S_STD
     &,SU10,SV10,ST1P5,SQ1P5
     &,CDR10M,CHR1P5M,LTIMER
     +)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      REAL                                                              
     + Z0M(LAND_FIELD)      ! IN Roughness length for momentum (m).     
     +,Z0H(LAND_FIELD)      ! IN Roughness length for heat and          
!                        !    moisture (m).                             
     &,Z0M_STD(LAND_FIELD)  ! IN Roughness length for momentum without  
!                        !    orographic component (m).                 
     &,VSHR(P_FIELD)        ! IN Wind speed difference between the
!                           !    surface and the lowest wind level in
!                           !    the atmosphere (m/s).
     &,CD(LAND_FIELD)       ! IN Surface drag coefficient.              
     &,CH(LAND_FIELD)       ! IN Surface transfer coefficient for heat a
!                        !    moisture.                                 
     &,CD_STD(LAND_FIELD)   ! IN Surface drag coefficient excluding     
!                        !    orographic from drag.                     
     &,TILE_FRAC(LAND_FIELD)                                            
!                        ! IN Tile fraction.                            
     &,RECIP_L_MO(LAND_FIELD)                                           
!                        ! IN Reciprocal of the Monin-Obukhov length (m)
     &,V_S(LAND_FIELD)      ! IN Surface layer scaling velocity includin
!                        !    orographic form drag (m/s).
     &,V_S_STD(LAND_FIELD)  ! IN Surface layer scaling velocity excludin
*D SFLINT6A.86
     +,CHR1P5M(LAND_FIELD)  ! OUT Interpolation coefficient for 1.5m    
*D SFLINT6A.88,SFLINT6A.89   

      REAL
     & RIB(LAND_FIELD)   ! DUMMY Used in 7A boundary layer scheme
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                        ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(LAND_FIELD)   ! DUMMY Used in 7A boundary layer scheme
     &,Z1(P_FIELD)       ! DUMMY Used in 7A boundary layer scheme

*D SFLINT6A.112
     &,PHI_M_OBS(LAND_FIELD)  ! Monin-Obukhov stability function for    
*D SFLINT6A.115
     &,PHI_H_OBS(LAND_FIELD)  ! Monin-Obukhov stability function for    
*D SFLINT6A.118,SFLINT6A.119  
*D SFLINT6A.124,SFLINT6A.126  
     + I,J,L       ! Loop counter (horizontal field index).
*D SFLINT6A.138,SFLINT6A.139  
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
*D SFLINT6A.141,SFLINT6A.142  
          Z_TEMP(I) = Z_OBS_TQ + Z0H(L) - Z0M(L)
*D SFLINT6A.144
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
*D SFLINT6A.155,ARN0F405.1838 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CHR1P5M(L) = CH(L) * VSHR(I) * PHI_H_OBS(L)/(VKMAN*V_S_STD(L))
*D SFLINT6A.166,ARN0F405.1841 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CDR10M(I) = CDR10M(I) + TILE_FRAC(L) *
     &                CD(L) * VSHR(I) * PHI_M_OBS(L)/(VKMAN*V_S(L))
*D ARN0F405.1844
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          Z_TEMP(I) = Z_OBS_TQ + Z0H(L) - Z0M_STD(L)
        ENDDO
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
*D ARN0F405.1847,ARN0F405.1848 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CDR10M(I) = CDR10M(I) + TILE_FRAC(L) *
     &                CD_STD(L) * VSHR(I) * PHI_M_OBS(L)/
     &                     (VKMAN*V_S_STD(L))
*/-----
*DECLARE SFLINT7A
*/-----
*D SFLINT7A.60,SFLINT7A.64   
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,FLANDG,        
     & VSHR,CD_STD,CD,CH,RIB,TILE_FRAC,WIND_PROFILE_FACTOR,             
     & Z0M_EFF,Z0M,Z0H,Z0F,Z1,                                          
     & RECIP_L_MO,V_S,V_S_STD,                                          
     & SU10,SV10,ST1P5,SQ1P5,                                           
     & CDR10M,CHR1P5M,LTIMER                                            
*D SFLINT7A.79
     & FLANDG(P_FIELD)   ! IN Land fraction                             
     &,CD_STD(LAND_FIELD)! IN Surface drag coefficient for shear stress 
*D SFLINT7A.100
     &,Z1(LAND_FIELD)       ! IN Height of lowest model level (m).
*I SFLINT7A.115   
      REAL                                                              
     & VSHR(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme  
     &,RECIP_L_MO(LAND_FIELD) ! DUMMY Used in 6A boundary layer scheme  
     &,V_S(LAND_FIELD)        ! DUMMY Used in 6A boundary layer scheme  
     &,V_S_STD(LAND_FIELD)    ! DUMMY Used in 6A boundary layer scheme  
                                                                        
*D SFLINT7A.171
          Z1E(L) = Z1(L) + Z0M_EFF(L)
*D SFLINT7A.177
          Z_OVER_Z1 = Z10M  / Z1(L)
*D SFLINT7A.198
     &                FLANDG(I)*TILE_FRAC(L)*CD10*WIND_PROFILE_FACTOR(L)
*D SFLINT7A.213
          Z1E(L) = Z1(L) + Z0M_EFF(L)
*D SFLINT7A.254,SFLINT7A.256  
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
     & VSHR,CD,CH,RIB,Z0M,Z0H,Z0F,Z1,                                   
     & RECIP_L_MO,V_S,                                                  
     & SU10,SV10,ST1P5,SQ1P5,                                           
     & CDR10M,CHR1P5M,LTIMER                                            
*D SFLINT7A.267
     & FLANDG(P_FIELD)    ! IN Land fraction.                         
     &,CD(P_FIELD)        ! IN Effective surface drag coefficient,      
*D SFLINT7A.281,SFLINT7A.282  
     & SU10               ! IN 10m U-wind diagnostic flag               
*I SFLINT7A.293   
      REAL                                                              
     & VSHR(P_FIELD)       ! DUMMY Used in 6A boundary layer scheme     
     &,RECIP_L_MO(P_FIELD) ! DUMMY Used in 6A boundary layer scheme     
     &,V_S(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme     
                                                                        
*D SFLINT7A.347
          IF ( FLANDG(I).LT.1.0 ) THEN                                
*I SFLINT7A.371   
            CDR10M(I) = CDR10M(I)*(1.0-FLANDG(I))
                                                                       
*D SFLINT7A.383
          IF ( FLANDG(I).LT.1.0 ) THEN                                
*/-----
*DECLARE SFMELT7A
*/-----
*D SFMELT7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I SFMELT7A.25    
C                                                                       
C           Modified for MOSES II. RE 4/4/00                            
C                                                                       
*D SFMELT7A.28,SFMELT7A.34   
     & POINTS,P_FIELD,P1,LAND_FIELD,NTILES,LAND_INDEX                   
     &,TILE_INDEX,TILE_PTS,LAND_MASK,LTIMER,SIMLT,SMLT,FLANDG           
     &,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_1,ICE_FRACT            
     &,RHOKH_1,RHOKH_1_SICE,TILE_FRAC,TIMESTEP                          
     &,EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE                     
     &,TSTAR_SEA,TSTAR_SSI,TSTAR_TILE,SNOW_TILE                         
     &,EI_LAND,EI_SICE                                                  
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE                  
     &,L_ESSERY_SNOW 
*D SFMELT7A.45,SFMELT7A.47   
     &,LAND_INDEX(P_FIELD)  ! IN Index of land points.                  
     &,NTILES               ! IN Number of tiles per land point.        
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
!                           ! IN Index of tile points.                  
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
!
*I SFMELT7A.53
     &, L_ESSERY_SNOW
!
*D SFMELT7A.56
     & FLANDG(P_FIELD)      ! IN Fraction of gridbox which is land.     
     &,ALPHA1(LAND_FIELD,NTILES)                                        
!                           ! IN Gradients of saturated specific        
*D SFMELT7A.59
!                           !    and land tile surfaces.                
*D SFMELT7A.63,SFMELT7A.64   
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                           ! IN Coefficient to calculate surface       
!                           !    heat flux into soil.                   
*D SFMELT7A.68,SFMELT7A.69   
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
!                           ! IN Surface exchange coeffs for land tiles.
*D SFMELT7A.72,SFMELT7A.73   
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
!                           ! IN Tile fractions.                        
*D SFMELT7A.77
     & EI_TILE(LAND_FIELD,NTILES)                                       
!                           ! INOUT Sublimation for land tiles (kg/m2/s)
     &,FQW_1(P_FIELD)       ! INOUT GBM surface moisture flux (kg/m2/s).
*D SFMELT7A.79
*D SFMELT7A.81,SFMELT7A.88   
     &,FTL_ICE(P_FIELD)     ! INOUT FTL for sea-ice.                    
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
!                           ! INOUT FTL for land tiles.                 
     &,TSTAR_SEA(P_FIELD)   ! IN Open sea surface temperature (K).      
     &,TSTAR_SSI(P_FIELD)   ! INOUT Sea mean surface temperature (K).   
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
!                           ! INOUT Land tile surface temperatures (K). 
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                           ! INOUT Lying snow on land tiles (kg/m2).   
*D SFMELT7A.91,SFMELT7A.92   
     & EI_LAND(P_FIELD)     ! OUT Sublimation from lying snow (kg/m2/s).
     &,EI_SICE(P_FIELD)     ! OUT Sublimation from sea-ice (kg/m2/s).   
     &,MELT_TILE(LAND_FIELD,NTILES)                                     
!                           ! OUT Surface snowmelt on tiles (kg/m2/s).  
*D SFMELT7A.98
     &,SNOWMELT(P_FIELD)    ! OUT GBM surface snowmelt (kg/m2/s).       
*D SFMELT7A.101
*I SFMELT7A.120   
     &,N                    ! Loop counter - tile index.                
cccccc Tibet Snow mod cccccc
      REAL MASKD
      PARAMETER( MASKD = 0.1 )
cccccccccccccccccccccccccccc
*D SFMELT7A.130
        EI_LAND(I) = 0.0                                              
        EI_SICE(I) = 0.0                                              
*D SFMELT7A.133,SFMELT7A.135  
      DO N=1,NTILES                                                     
        DO L=1,LAND_FIELD                                               
          MELT_TILE(L,N) = 0.                                           
        ENDDO                                                           
      ENDDO                                                             
                                                                        
*D SFMELT7A.137
!  Melt snow on land tiles if TSTAR_TILE is greater than TM.            
*D SFMELT7A.139,SFMELT7A.153  
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          SNOW_MAX = MAX( 0.0, SNOW_TILE(L,N) - EI_TILE(L,N)*TIMESTEP ) 
          IF ( SNOW_MAX.GT.0.0 .AND. TSTAR_TILE(L,N).GT.TM ) THEN       
            LCMELT = (CP + LC*ALPHA1(L,N))*RHOKH_1(L,N)                 
     &               + ASHTF_TILE(L,N)                                  
            LSMELT = LCMELT + LF*ALPHA1(L,N)*RHOKH_1(L,N)               
cccccc Tibet Snow mod cccccc
c           DTSTAR = - MIN( TSTAR_TILE(L,N) - TM ,                      
c    &                      LF*SNOW_MAX / (LCMELT*TIMESTEP) )           
            DTSTAR = - MIN( (TSTAR_TILE(L,N) - TM) *
     &                      (1.0 - EXP(-MASKD*SNOW_MAX)) ,
     &                      LF*SNOW_MAX / (LCMELT*TIMESTEP) )           
cccccccccccccccccccccccccccc
            MELT_TILE(L,N) = - LSMELT*DTSTAR / LF                       
            DFTL = CP*RHOKH_1(L,N)*DTSTAR                               
            DFQW = ALPHA1(L,N)*RHOKH_1(L,N)*DTSTAR                      
            FTL_TILE(L,N) = FTL_TILE(L,N) + DFTL                        
            EI_TILE(L,N) = EI_TILE(L,N) + DFQW                          
            TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR                  
*D SFMELT7A.157,SFMELT7A.166  
            DFTL = TILE_FRAC(L,N)*DFTL                                  
            DFQW = TILE_FRAC(L,N)*DFQW                                  
            FTL_1(I) = FTL_1(I) + FLANDG(I)*DFTL                  
            FQW_1(I) = FQW_1(I) + FLANDG(I)*DFQW                  
          ENDIF                                                         
          EI_LAND(I) = EI_LAND(I) + TILE_FRAC(L,N)*EI_TILE(L,N)     
        ENDDO                                                           
*I SFMELT7A.168   
!-----------------------------------------------------------------------
!  Increment snow by sublimation and melt                               
!-----------------------------------------------------------------------
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             

!         IF (L_ESSERY_SNOW) THEN
! ALL ADJUSTMENTS TO SNOW MASS DONE IN SEPARATE SNOW ROUTINES NOW
!           SNOW_TILE(L,N) = SNOW_TILE(L,N) - 
!     &                     (EI_TILE(L,N) * TIMESTEP)
!         ELSE
          IF (.NOT.L_ESSERY_SNOW) THEN
          SNOW_TILE(L,N) = SNOW_TILE(L,N) -                             
     &                     (EI_TILE(L,N) + MELT_TILE(L,N))*TIMESTEP     
          END IF
          SNOWMELT(I) = SNOWMELT(I) + TILE_FRAC(L,N)*MELT_TILE(L,N)     
        ENDDO                                                           
      ENDDO                                                             
      IF (SMLT) THEN                                                    
        DO I=P1,P1+POINTS-1                                             
          SNOMLT_SURF_HTF(I) = LF*SNOWMELT(I)                           
        ENDDO                                                           
      ENDIF                                                             
                                                                        
*D SFMELT7A.170
        IF ( FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.0 ) THEN      
*D SFMELT7A.172
!   Melt sea-ice if TSTAR > TSTARMAX                                    
*D SFMELT7A.174,SFMELT7A.176  
          EI_SICE(I) = FQW_ICE(I)                                   
          TSTARMAX = ICE_FRACT(I)*TM                                  
     &         + (1.0 - ICE_FRACT(I))*TSTAR_SEA(I)                  
          IF ( TSTAR_SSI(I) .GT. TSTARMAX ) THEN                      
*D SFMELT7A.178,SFMELT7A.179  
     &                       + ICE_FRACT(I)*GAMMA(1)*DTRDZ_1(I) )   
            DTSTAR = TSTARMAX - TSTAR_SSI(I)                          
*D SFMELT7A.181
     &                                                 + ASHTF(I)     
*D SFMELT7A.184,SFMELT7A.185  
            TSTAR_SSI(I) = TSTARMAX                                   
*D SFMELT7A.187,SFMELT7A.196  
            FTL_1(I) = FTL_1(I) + (1.0-FLANDG(I))*DFTL            
            FQW_1(I) = FQW_1(I) + (1.0-FLANDG(I))*DFQW            
            EI_SICE(I) = EI_SICE(I) + DFQW                          
            FTL_ICE(I) = FTL_ICE(I) + DFTL                          
            FQW_ICE(I) = FQW_ICE(I) + DFQW                          
*/-----
*DECLARE SFOROG7A
*/-----
*D  SFOROG7A.56
     &,Z1(LAND_FIELD)          ! IN Height of lowest atmospheric level (m).
*D SFOROG7A.114
          ZETA2 = LOG ( 1.0 + Z1(L)/Z0M(L) )
*D SFOROG7A.117
          H_BLEND_OROG = MAX ( Z1(L) / (1.0 - ZETA2) ,
*/-----
*DECLARE SFREST7A
*/-----
*D SFREST7A.30
     & CANOPY,CATCH,CH,DQ,EPDT,FLAKE,GC,SNOW,VSHR,                      
*I SFREST7A.56    
     &,FLAKE(LAND_FIELD)   ! IN Lake fraction.                          
*I SFREST7A.58    
     &,SNOW(LAND_FIELD)    ! IN Lying snow amount (kg per sq metre).    
*I SFREST7A.75    
cccccc Tibet Snow mod cccccccc
      REAL MASKD
      PARAMETER( MASKD = 0.1 )
cccccccccccccccccccccccccccccc
                                                                        
*D SFREST7A.93,SFREST7A.94   
! Set to 1 for negative moisture flux or snow-covered land              
! (no surface/stomatal resistance to condensation).                     
*D SFREST7A.97,SFREST7A.98   
cccccc Tibet Snow mod cccccccc
        IF(SNOW(L) .GT. 0.) FRACA(L) = 1.0 - EXP(-MASKD*SNOW(L))
cccccccccccccccccccccccccccccc
        IF (DQ(L).LT.0. .AND. SNOW(L).LE.0.) FRACA(L) = 0.0             
        IF (DQ(L).LT.0. .AND. SNOW(L).LE.0. .AND. CATCH(L).GT.0.)       
*D SFREST7A.106
        RESFT(L) = FLAKE(L) + (1. - FLAKE(L)) *                         
     &                        ( FRACA(L) + (1. - FRACA(L))*RESFS(L) )   
*/-----
*DECLARE SFRIB7A
*/-----
*D SFRIB7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFRIB7A.37
     & RIB,DB,LTIMER                                                    
*D SFRIB7A.53,SFRIB7A.55
     & BQ_1(LAND_FIELD)       ! IN A buoyancy parameter for lowest atm
!                          !    level. ("beta-q twiddle").
     &,BT_1(LAND_FIELD)       ! IN A buoyancy parameter for lowest atm
*D SFRIB7A.58
     &,QW_1(LAND_FIELD)       ! IN Total water content of lowest
*D SFRIB7A.61
     &,TL_1(LAND_FIELD)       ! IN Liquid/frozen water temperature for
*D SFRIB7A.68,SFRIB7A.69
     &,Z1_TQ(LAND_FIELD)      ! IN Height of lowest TQ level (m).
     &,Z1_UV(LAND_FIELD)      ! IN Height of lowest UV level (m).
*I SFRIB7A.72    
     &,DB(LAND_FIELD)      ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level.          
*D SFRIB7A.102,SFRIB7A.104
        DTEMP(L) = TL_1(L) - TSTAR(L) + (G/CP)*(Z1_TQ(L)+Z0M(L)-Z0H(L))
!                                                             ! P243.118
        DQ(L) = QW_1(L) - QSTAR(L)                            ! P243.119
*D SFRIB7A.113,SFRIB7A.114  
        DB(L) = G*(BT_1(L)*DTEMP(L) + BQ_1(L)*RESFT(L)*DQ(L))           
        RIB(L) = Z1_UV(L)*DB(L) / ( VSHR(I)*VSHR(I) )                   
*D SFRIB7A.130,SFRIB7A.131  
     & P_POINTS,P_FIELD,P1,FLANDG,NSICE,SICE_INDEX,                  
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    
*D SFRIB7A.133
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             
*D SFRIB7A.147
*D SFRIB7A.150
     & FLANDG(P_FIELD)     ! IN Land fraction on all pts.               
     &,BQ_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm     
*D SFRIB7A.163
     &,TSTAR_SICE(P_FIELD)  ! IN Surface temperature of sea-ice (K).    
*I SFRIB7A.183   
     &,DB_SEA(P_FIELD)     ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level over      
!                          !     sea or sea-ice leads.                  
     &,DB_ICE(P_FIELD)     ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level over      
!                          !     sea-ice.                               
*D SFRIB7A.207
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D SFRIB7A.212,SFRIB7A.213  
          DB_SEA(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )                  
          RIB_SEA(I) = Z1_UV(I)*DB_SEA(I) / ( VSHR(I)*VSHR(I) )         
*D SFRIB7A.220
        DTEMP = TL_1(I) - TSTAR_SICE(I)                                 
*D SFRIB7A.223,SFRIB7A.224  
        DB_ICE(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )                    
        RIB_ICE(I) = Z1_UV(I)*DB_ICE(I) / ( VSHR(I) * VSHR(I) )         
*/-----
*DECLARE SFSNOW7A
*/-----
*D SFSNOW7A.23,SFSNOW7A.28   
CLL  Purpose:  adds the large-scale and convective snowfall to the
CLL            snowdepth;  
*D SFSNOW7A.45,SFSNOW7A.49   
     & NPNTS,NTILES,TILE_PTS,TILE_INDEX,
     & CONV_SNOW,LS_SNOW,TILE_FRAC,TSTAR_TILE,TIMESTEP,      
     & RGRAIN,SNOW_TILE,L_SNOW_ALBEDO,
     & LTIMER)            
*D SFSNOW7A.55,SFSNOW7A.56   
     &,NTILES               ! IN Number of tiles.
     &,TILE_PTS(NTILES)     ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                           ! IN Index of tile points.  
*D SFSNOW7A.61,SFSNOW7A.67   
     &,TILE_FRAC(NPNTS,NTILES)
                            ! IN Tile fractions.                        
     &,TSTAR_TILE(NPNTS,NTILES)
!                           ! IN Tile surface temperatures (K).
*D SFSNOW7A.69
*D SFSNOW7A.72,SFSNOW7A.80   
     & RGRAIN(NPNTS,NTILES) ! INOUT Snow grain size (microns).
     &,SNOW_TILE(NPNTS,NTILES)
!                           ! INOUT Snow on the ground (kg/m2). 
*D SFSNOW7A.83,SFSNOW7A.84   
     & L_SNOW_ALBEDO        ! IN Flag for prognostic snow albedo.       
*D SFSNOW7A.87,SFSNOW7A.89   
*CALL C_0_DG_C
*D SFSNOW7A.93,SFSNOW7A.95   
     & R0                  ! Grain size for fresh snow (microns).       
*D SFSNOW7A.98,SFSNOW7A.100  
     &,SNOWFALL(NPNTS)     ! Snowfall in timestep (kg/m2).
*D SFSNOW7A.102
      INTEGER I,J,N        ! Loop counters.                             
*D SFSNOW7A.109,SFSNOW7A.154  
*D SFSNOW7A.161,SFSNOW7A.162  
        SNOWFALL(I) = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I))
      ENDDO

      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N) 
          SNOW_TILE(I,N) = SNOW_TILE(I,N) + SNOWFALL(I)
        ENDDO
*D SFSNOW7A.169,SFSNOW7A.171  
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N) 
          IF ( SNOW_TILE(I,N) .GT. 0.) THEN                             
*D SFSNOW7A.173,SFSNOW7A.174  
            IF (TSTAR_TILE(I,N) .LT. TM) THEN                           
              IF (RGRAIN(I,N) .LT. 150.) THEN                           
*D SFSNOW7A.177
                RATE = 0.23E6*EXP(-3.7E4/(8.13451*TSTAR_TILE(I,N)))     
*D SFSNOW7A.180,SFSNOW7A.183  
            RGRAIN(I,N) = SQRT( RGRAIN(I,N)**2 
     &                          + (RATE/3.14159)*TIMESTEP )  
     &                              - (RGRAIN(I,N) - R0)*SNOWFALL(I)/2.5
            RGRAIN(I,N) = MIN( RMAX, RGRAIN(I,N) )                      
            RGRAIN(I,N) = MAX( R0, RGRAIN(I,N) )                        
*D SFSNOW7A.185
            RGRAIN(I,N) = R0                                            
*I SFSNOW7A.187   
      ENDDO                                                           
*/-----
*DECLARE SFSTOM7A
*/-----
*D SFSTOM7A.123
      PARAMETER (ITER = 3)

*/-----
*DECLARE SICEHT5B
*/-----
*I SICEHT5B.44    
!!!
!!! Note: 
!!!   Routine renamed to avoid clash with SICEHT7A which replaces 
!!!   this deck in MOSES2.2 7A boundary layer. 
!!!
*D SICEHT5B.48
      SUBROUTINE SICE_HTF_5B (
*/-----
*DECLARE SIEVE7A
*/-----
*I SIEVE7A.50
     & ,TFALL_TILE,SOIL
*I SIEVE7A.58
     &,N
*I SIEVE7A.70
     &,TFALL_TILE(NPNTS)     ! INOUT Cummulative canopy throughfall

      LOGICAL SOIL
*D SIEVE7A.95
        TFALL_TILE(I) = TFALL_TILE(I)+TFALL(I)
        TOT_TFALL(I) = TOT_TFALL(I) + FRAC(I)*TFALL(I)
*/-----
*DECLARE SMCEXT7A
*/-----
*D SMCEXT7A.24
     &,                   F_ROOT,STHU,V_CRIT_NEW,V_SAT,V_WILT
*D SMCEXT7A.67
*D SMCEXT7A.71
     &,V_CRIT_NEW(NPNTS)
*D SMCEXT7A.120
     &               /(V_CRIT_NEW(I)-V_WILT(I))
*D SMCEXT7A.137,SMCEXT7A.138  
     &      WT_EXT(I,N) = F_ROOT(N)*FSMC_L(I,N)/FSMC(I)
  
*/-----
*DECLARE SOILHT7A
*/-----
*D SOILHT7A.24,SOILHT7A.26   
     & NPNTS,NSHYD,NTILES,NELEV,SOIL_PTS,SOIL_INDEX,TILE_PTS,TILE_INDEX
     &,BEXP,DZ,FRAC,HCAP,HCON,SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP
     &,V_SAT,W_FLUX,SMCL,STHU,STHF,TSOIL                        
     &,NSNOW,L_ESSERY_SNOW,TSTAR_NONICE
*I SOILHT7A.37    
!     GAUSS - to solve matrix equation for temperature increments       
*I GPB8F405.40    
!  4.6     10/99     Implicit updating. P. Cox, R. Essery  
*I SOILHT7A.63    
     &,NTILES               ! IN Number of tiles.       
     &,NELEV
*I SOILHT7A.72    
     &,TILE_PTS(NTILES)     ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                           ! IN Index of tile points.                  
*D SOILHT7A.75
     & BEXP(NPNTS)          ! IN Clapp-Hornberger exponent.             
*I SOILHT7A.76    
     &,FRAC(NPNTS,NTILES)   ! IN Tile fractions.    
*D SOILHT7A.83,SOILHT7A.84   
     &,SNOW_TILE(NPNTS,NTILES)
!                           ! IN Lying snow on tiles (kg/m2). 
     &,NSNOW(NPNTS,NTILES)  ! IN number of snow layers
     &,SURF_HT_FLUX(NPNTS)  ! IN Net downward surface heat flux (W/m2). 
     &,TSTAR_NONICE(NPNTS)  ! IN Used for temp of infiltrating water. Use
                            ! of TSOIL(1) previously led to +ve feedbacks 
*D SOILHT7A.108
     &,ITER_PTS             ! WORK Number of soil points which require  
*I SOILHT7A.116   
     &,SI_TILE              ! WORK Tile snow insulation factor.
     &,SNOW_DEPTH           ! WORK Depth of lying snow (m). 
     &,ICE_FRAC(NPNTS)
*D SOILHT7A.123
     & ITER_INDEX(NPNTS)    ! WORK Array of soil points which require   
*D SOILHT7A.133,SOILHT7A.134  
*D SOILHT7A.142
     &,DTSLMAX(NPNTS,NSHYD) ! WORK Maximum value of DTSL (K/timestep).  
     &,DTSLMIN(NPNTS,NSHYD) ! WORK Minimum value of DTSL (K/timestep).
     &,HCAPT(NPNTS,NSHYD)  ! WORK The total volumetric heat capacity   
*I SOILHT7A.149   
     &,HADV(NPNTS,NSHYD)    ! WORK Heat flux due to moisture advection
!                           !      (W/m2).
     &,SIFACT(NPNTS)        ! WORK Snow insulation factor.
*I SOILHT7A.168   
      LOGICAL
     & ITER(NPNTS),L_ESSERY_SNOW ! WORK .T. on points requiring iterations.
     & ,DEEPSNOW

!-----------------------------------------------------------------------
! Variables required for the implicit calculation.
!-----------------------------------------------------------------------
      REAL
     & DHFLUX_DTSL1(NPNTS,0:NSHYD),DHFLUX_DTSL2(NPNTS,0:NSHYD)
     &,DHADV_DTSL0(NPNTS,NSHYD),DHADV_DTSL1(NPNTS,NSHYD)
     &,DHADV_DTSL2(NPNTS,NSHYD) ! WORK Rate of change of the explicit
!                           ! fluxes with the layer temperatures
!                           ! (W/m2/K).  
     &,A(NPNTS,NSHYD),B(NPNTS,NSHYD),C(NPNTS,NSHYD),D(NPNTS,NSHYD)
!                           ! WORK Matrix elements.
     &,GAMCON               ! WORK Forward timestep weighting constant.
! Local parameters:
      REAL
     & GAMMA                ! Forward timestep weighting.
      PARAMETER (GAMMA=1.0)    
*D SOILHT7A.172
     & HEAT_CON,GAUSS                                                   
*D SOILHT7A.190
        TSL(I,0)=MAX(0.,TSTAR_NONICE(I)-ZERODEGC)
*I SOILHT7A.216   
!-----------------------------------------------------------------------
! Initialise the array of points for which calculations are required.   
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        ITER(I)=.FALSE.
      ENDDO

      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)
        ITER(I)=.TRUE.                                                  
      ENDDO                  
                                                                        
*I SOILHT7A.238   
! Calculate the snow insulation factor                               
!--------------------------------------------------------------------
      DO I=1,NPNTS                                               
        SIFACT(I) = 0.                                              
        ICE_FRAC(I) = 0.
        DO N=NTILES-NELEV+1,NTILES
          ICE_FRAC(I) = ICE_FRAC(I)+FRAC(I,N)
        ENDDO
      ENDDO         
!Essery snow NSNOW>0 has already absorbed the heat flux
      DO N=1,NTILES-NELEV
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)

          DEEPSNOW=.FALSE.
          IF (L_ESSERY_SNOW) THEN
            IF (NSNOW(I,N).gt.0) DEEPSNOW=.TRUE.
          ENDIF
          IF (.NOT.DEEPSNOW) THEN

          SI_TILE = 1.
          IF (SNOW_TILE(I,N).GT.0.) THEN
            SNOW_DEPTH = SNOW_TILE(I,N) / RHO_SNOW
            IF (SNOW_DEPTH .LE. (0.5*DZ(1))) THEN
              SI_TILE = 1. / ( 1. + 2*SNOW_DEPTH/(DZ(1) + DZ(2)) )
            ELSE 
              SI_TILE =(DZ(1) + DZ(2)) /
     &                 ( HC(I,1)*(2*SNOW_DEPTH - DZ(1))/SNOW_HCON
     &                   + 2*DZ(1) + DZ(2) )
            ENDIF
          ENDIF

          if (ICE_FRAC(I).lt.1) 
          !Normalise by non-ice fracftion. ice-frac(i) really should 
          !be <1 if we're here!
     &     SIFACT(I) = SIFACT(I) + FRAC(I,N)*SI_TILE/(1-ICE_FRAC(I))
         ENDIF
        ENDDO
      ENDDO
                                                                        
!--------------------------------------------------------------------   
*I SOILHT7A.245   
          DHFLUX_DTSL1(I,N)=HC(I,N)*2.0/(DZ(N+1)+DZ(N))
          DHFLUX_DTSL2(I,N)=-HC(I,N)*2.0/(DZ(N+1)+DZ(N))                
*D ARE1F405.44
        H_FLUX(I,0) = SURF_HT_FLUX(I)
        H_FLUX(I,1) = SIFACT(I)*H_FLUX(I,1)
        DHFLUX_DTSL1(I,NSHYD)=0.0
        DHFLUX_DTSL2(I,NSHYD)=0.0
        DHFLUX_DTSL1(I,0)=0.0
        DHFLUX_DTSL2(I,0)=0.0
        DHFLUX_DTSL1(I,1)=SIFACT(I)*DHFLUX_DTSL1(I,1)
        DHFLUX_DTSL2(I,1)=SIFACT(I)*DHFLUX_DTSL2(I,1)               
*D SOILHT7A.262
          HADV(I,N)=HCAPW*DZ(N)* 
*I SOILHT7A.264   
          DHADV_DTSL0(I,N)=HCAPW*DZ(N)*W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))    
          DHADV_DTSL1(I,N)=HCAPW*DZ(N)*                                 
     &                     (-W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))
     &                      +W_FLUX(I,N)/(DZ(N)+DZ(N+1)))            
          DHADV_DTSL2(I,N)=-HCAPW*DZ(N)*W_FLUX(I,N)/(DZ(N)+DZ(N+1))     
*D SOILHT7A.271
        HADV(I,1)=HCAPW*DZ(1)*
*D SOILHT7A.274
        DHADV_DTSL0(I,1)=0.0                
        DHADV_DTSL1(I,1)=HCAPW*DZ(1)*                                   
     &    (-W_FLUX(I,0)/DZ(1)+W_FLUX(I,1)/(DZ(1)+DZ(2)))            
        DHADV_DTSL2(I,1)=-HCAPW*DZ(1)*W_FLUX(I,1)/(DZ(1)+DZ(2))
        HADV(I,NSHYD)=HCAPW*DZ(NSHYD)*
*I SOILHT7A.276   
        DHADV_DTSL0(I,NSHYD)=HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))                
        DHADV_DTSL1(I,NSHYD)=-HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))
        DHADV_DTSL2(I,NSHYD)=0.0                           
*D ACB2F405.15,ACB2F405.17   
! unfrozen. Check that (SMCLSAT/SMCL)**BEXP will not overflow when SMCL 
! is very small. The function EPSILON  gives the number of type (real)  
! of SMCL that is negligeable compared to 1.                            
*D SOILHT7A.289
     &             *(SMCLSAT(I,N)/SMCL(I,N))**(BEXP(I))                 
*D SOILHT7A.295,SOILHT7A.297  
          DHSL0(I,N)=TIMESTEP*(H_FLUX(I,N-1)-H_FLUX(I,N)+HADV(I,N))     
*I SOILHT7A.300   
        ENDDO
      ENDDO                                          
                                                                        
*D SOILHT7A.302
! Iteration loop                                                     
*I SOILHT7A.303   
      DO M=1,MMAX                                                       
      

!-----------------------------------------------------------------------
! Define the array of points which fail to meet the flux criterion.     
!-----------------------------------------------------------------------
      ITER_PTS=0                                                        
      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)                                                 
                                                                        
        IF (ITER(I)) THEN                                   
          ITER_PTS=ITER_PTS+1                                 
          ITER_INDEX(ITER_PTS)=I                                        
        ENDIF                                                           
        ITER(I)=.FALSE.
                                                                        
      ENDDO                                                             

!-----------------------------------------------------------------------
! Update calculations at these points.  
!-----------------------------------------------------------------------
      DO N=1,NSHYD

        DO J=1,ITER_PTS
          I=ITER_INDEX(J)

*I SOILHT7A.307   
          DTSLMAX(I,N)=1.0E4-TSL(I,N)
          DTSLMIN(I,N)=-ZERODEGC-TSL(I,N)

*I SOILHT7A.309   
            DTSLMIN(I,N)=TMAX(I,N)-TSL(I,N)  
*D SOILHT7A.315,SOILHT7A.316  
     &               /(BEXP(I)*SATHH(I)*RHO_WATER*DZ(N))              
     &            *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I)-1.0)     
            DTSLMAX(I,N)=TMAX(I,N)-TSL(I,N)
*D SOILHT7A.319
          HCAPT(I,N)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I,N)/DZ(N)         
*D SOILHT7A.323,SOILHT7A.324  
*D SOILHT7A.327,SOILHT7A.328  
! Calculate the matrix elements required for the implicit update.
*D SOILHT7A.330,SOILHT7A.372  
          GAMCON=GAMMA*TIMESTEP/(HCAPT(I,N)*DZ(N))
          A(I,N)=-GAMCON*(DHFLUX_DTSL1(I,N-1)+DHADV_DTSL0(I,N))
          B(I,N)=1.0-GAMCON*(DHFLUX_DTSL2(I,N-1)-DHFLUX_DTSL1(I,N)
     &                                          +DHADV_DTSL1(I,N))
          C(I,N)=GAMCON*(DHFLUX_DTSL2(I,N)+DHADV_DTSL2(I,N))
          D(I,N)=1.0/(HCAPT(I,N)*DZ(N))*DHSL(I,N)                       
*D SOILHT7A.375,SOILHT7A.389  
      ENDDO
*D SOILHT7A.393
! Solve the triadiagonal matrix equation.
*D SOILHT7A.395,SOILHT7A.398  
      CALL GAUSS(NSHYD,NPNTS,ITER_PTS,ITER_INDEX,A,B,C,D
     &,          DTSLMIN,DTSLMAX,DTSL)
*D SOILHT7A.400,SOILHT7A.402  
!-----------------------------------------------------------------------
! Diagnose the implicit DHSL
!-----------------------------------------------------------------------
      DO N=2,NSHYD-1
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)
          DHSL(I,N)=DHSL(I,N)-DZ(N)*HCAPT(I,N)*(A(I,N)*DTSL(I,N-1)
     &                    +(B(I,N)-1)*DTSL(I,N)+C(I,N)*DTSL(I,N+1)) 
        ENDDO
      ENDDO
*D SOILHT7A.404,SOILHT7A.413  
      DO J=1,ITER_PTS
        I=ITER_INDEX(J)
        DHSL(I,1)=DHSL(I,1)-DZ(1)*HCAPT(I,1)*(
     &                  +(B(I,1)-1)*DTSL(I,1)+C(I,1)*DTSL(I,2)) 
        DHSL(I,NSHYD)=DHSL(I,NSHYD)-DZ(NSHYD)*HCAPT(I,NSHYD)*
     &  (A(I,NSHYD)*DTSL(I,NSHYD-1)+(B(I,NSHYD)-1)*DTSL(I,NSHYD)) 
      ENDDO
*D SOILHT7A.415,SOILHT7A.417  
!-----------------------------------------------------------------------
! Update the layer temperatures
!-----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)
*D SOILHT7A.419
*D SOILHT7A.436,SOILHT7A.445  
*D SOILHT7A.453
     &                *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I))     
*D SOILHT7A.465
!-----------------------------------------------------------------------
! Calculate the error in flux conservation                              
!-----------------------------------------------------------------------
          CEACUR(I)=ABS(DHSL(I,N))/TIMESTEP                          

          IF (CEACUR(I) .GT. FACUR) THEN                                
            ITER(I)=.TRUE.
          ENDIF

*I SOILHT7A.466   
      ENDDO
                                                                        
!-----------------------------------------------------------------------
! End of iteration loop  
!-----------------------------------------------------------------------
*/-----
*DECLARE SOILHY5A
*/-----
*D SOILHY5A.23,SOILHY5A.25   
      SUBROUTINE SOIL_HYD (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX              
     &,                    BEXP,DZ,EXT,FW,KS,SATHH,TIMESTEP,V_SAT       
     &,                    SLOW_RUNOFF,SMCL,STHU,SURF_ROFF,W_FLUX       
     &,                    STF_SLOW_RUNOFF,LTIMER                       
*I SOILHY5A.48    
!  4.6      2/99     Modified for implicit updating.  Peter Cox         
*D SOILHY5A.82
     & BEXP(NPNTS)          ! IN Clapp-Hornberger exponent.             
*I SOILHY5A.111   
     &,SURF_ROFF(NPNTS)     ! INOUT Surface runoff (kg/m2/s).           
*I SOILHY5A.117   
      REAL                                                              
     & GAMCON               ! WORK Constant (s/mm).                     
                                                                        
*D SOILHY5A.120,SOILHY5A.124  
     & A(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n-1).   
     &,B(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n).     
     &,C(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n+1).   
     &,D(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the RHS of the equation.          
     &,DSMCL(NPNTS,NSHYD)   ! WORK Soil moisture increment              
!                           !      (kg/m2/timestep).                    
     &,DSTHU(NPNTS,NSHYD)   ! WORK Increment to STHU (/timestep).       
     &,DSTHUMIN(NPNTS,NSHYD)! WORK Minimum value of DSTHU.              
     &,DSTHUMAX(NPNTS,NSHYD)! WORK Maximum value of DSTHU.              
     &,DWFLUX_DSTHU1(NPNTS,NSHYD) ! WORK The rate of change of the expli
!                           !      flux with STHU1 (kg/m2/s).           
     &,DWFLUX_DSTHU2(NPNTS,NSHYD) ! WORK The rate of change of the expli
!                           !      flux with STHU2 (kg/m2/s).           
*D SOILHY5A.127,SOILHY5A.128  
                                                                        
! Local parameters:                                                     
      REAL                                                              
     & GAMMA                ! Forward timestep weighting.               
      PARAMETER (GAMMA=1.0)                                             
*D SOILHY5A.141,GRB0F405.538  
*I SOILHY5A.146   
          DSTHUMIN(I,N)=-STHU(I,N)                                      
          DSTHUMAX(I,N)=1.0-SMCL(I,N)/SMCLSAT(I,N)                      
          DWFLUX_DSTHU1(I,N)=0.0                                        
          DWFLUX_DSTHU2(I,N)=0.0                                        
*D SOILHY5A.159
! Calculate the Darcian fluxes and their dependencies on the soil       
! moisture contents.                                                    
*I SOILHY5A.160   
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KS,STHU(1,NSHYD)     
     &,             W_FLUX(1,NSHYD),DWFLUX_DSTHU1(1,NSHYD)              
     &,             LTIMER)                                             
                                                                        
*I SOILHY5A.161   
        CALL DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KS,SATHH             
     &,             STHU(1,N-1),DZ(N-1),STHU(1,N),DZ(N),W_FLUX(1,N-1)   
     &,             DWFLUX_DSTHU1(1,N-1),DWFLUX_DSTHU2(1,N-1)           
     &,             LTIMER)                                             
       ENDDO                                                            
                                                                        
!-----------------------------------------------------------------------
! Calculate the explicit increments.                                    
!-----------------------------------------------------------------------
ccccc runoff_MII mod cccc
c     DO N=NSHYD,1,-1
      DO N=1,NSHYD,1
ccccccccccccccccccccccccc
*D SOILHY5A.164,SOILHY5A.166  
          DSMCL(I,N)=(W_FLUX(I,N-1)-W_FLUX(I,N)-EXT(I,N))*TIMESTEP      
*D SOILHY5A.169,SOILHY5A.170  
! Limit the explicit fluxes to prevent supersaturation.                 
*D SOILHY5A.172,SOILHY5A.275  
          IF (DSMCL(I,N).GT.(SMCLSAT(I,N)-SMCL(I,N))) THEN              
            DSMCL(I,N)=SMCLSAT(I,N)-SMCL(I,N)                           
ccccc runoff_MII mod cccc
c           W_FLUX(I,N-1)=DSMCL(I,N)/TIMESTEP+W_FLUX(I,N)+EXT(I,N)
            W_FLUX(I,N)=W_FLUX(I,N-1)-DSMCL(I,N)/TIMESTEP-EXT(I,N)
ccccccccccccccccccccccccc
*D SOILHY5A.282,SOILHY5A.283  
! Calculate the matrix elements required for the implicit update.       
*D SOILHY5A.285,GRB0F405.548  
*D SOILHY5A.288,SOILHY5A.289  
        GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,1)                              
        A(I,1)=0.0                                                      
        B(I,1)=1.0+GAMCON*DWFLUX_DSTHU1(I,1)                            
        C(I,1)=GAMCON*DWFLUX_DSTHU2(I,1)                                
        D(I,1)=DSMCL(I,1)/SMCLSAT(I,1)                                  
      ENDDO                                                             
*D SOILHY5A.291
      DO N=2,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,N)                            
          A(I,N)=-GAMCON*DWFLUX_DSTHU1(I,N-1)                           
          B(I,N)=1.0-GAMCON*(DWFLUX_DSTHU2(I,N-1)-DWFLUX_DSTHU1(I,N))   
          C(I,N)=GAMCON*DWFLUX_DSTHU2(I,N)                              
          D(I,N)=DSMCL(I,N)/SMCLSAT(I,N)                                
        ENDDO                                                           
      ENDDO                                                             
*D SOILHY5A.293,SOILHY5A.295  
!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.                               
!-----------------------------------------------------------------------
      CALL GAUSS(NSHYD,NPNTS,SOIL_PTS,SOIL_INDEX,A,B,C,D                
     &,          DSTHUMIN,DSTHUMAX,DSTHU)                               
*D SOILHY5A.297
!-----------------------------------------------------------------------
! Diagnose the implicit fluxes.                                         
!-----------------------------------------------------------------------
      DO N=1,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          DSMCL(I,N)=DSTHU(I,N)*SMCLSAT(I,N)                            
          W_FLUX(I,N)=W_FLUX(I,N-1)-EXT(I,N)-DSMCL(I,N)/TIMESTEP        
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
! Update the prognostic variables.                                      
!-----------------------------------------------------------------------
      DO N=1,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          SMCLU(I,N)=SMCLU(I,N)+DSMCL(I,N)                              
          SMCL(I,N)=SMCL(I,N)+DSMCL(I,N)                                
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)                             
        ENDDO                                                           
*I SOILHY5A.311   
                                                                        
!-----------------------------------------------------------------------
! Update surface runoff diagnostic.                                     
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)                                                 
        SURF_ROFF(I)=SURF_ROFF(I)+(FW(I)-W_FLUX(I,0))                   
      ENDDO                                                             
*/-----
*DECLARE SOILMC7A
*/-----
*D SOILMC7A.77
            SMC(I) = SMC(I) + RHO_WATER * ( ZSMC - Z1 ) * 
*/-----
*DECLARE SPARM1A
*/-----
*D SPARM1A.27,SPARM1A.29   
      SUBROUTINE SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES
     &,                 TILE_PTS,TILE_INDEX,FRAC,HT,LAI,SATCON
     &,                 CATCH_TILE,INFIL_TILE,Z0_TILE)           
*I SPARM1A.38    
     &,NTILES                ! IN Number of surface tiles.
*D SPARM1A.45,SPARM1A.46   
     & FRAC(LAND_FIELD,NTYPE)     ! IN Fractional cover of each       
*I SPARM1A.49    
     &,SATCON(LAND_FIELD)         ! IN Saturated hydraulic
!                                 !    conductivity (kg/m2/s).        
*D SPARM1A.52,SPARM1A.58   
     & CATCH_TILE(LAND_FIELD,NTILES)! OUT Canopy capacity for each tile 
!                                   !     (kg/m2).
     &,INFIL_TILE(LAND_FIELD,NTILES)! OUT Maximum surface infiltration 
!                                   !     rate for each tile (kg/m2/s).
     &,Z0_TILE(LAND_FIELD,NTILES)   ! OUT Roughness length for each     
!                                   !     tile (m).                     
     &,FZ0_TILE(LAND_FIELD,NTILES)
*D SPARM1A.60,SPARM1A.63   
     & CATCH(LAND_FIELD)          ! WORK GBM canopy capacity (kg/m2). 
     &,CATCH_T(LAND_FIELD,NTYPE)  ! WORK Capacities for types.
*I SPARM1A.64    
     &,INFIL(LAND_FIELD)          ! WORK GBM infiltration rate(kg/m2/s).
     &,INFIL_T(LAND_FIELD,NTYPE)  ! WORK Infiltration rates for types.
     &,Z0(LAND_FIELD)             ! WORK GBM roughness length (m).
     &,Z0_T(LAND_FIELD,NTYPE)     ! WORK Roughness lengths for types.   
     &,FRAC_LAKE
     &,FRAC_ICE(LAND_FIELD)
     &,elev_frac(LAND_FIELD,NELEV)
     &,frac_lake_elev(LAND_FIELD,NELEV)
*D SPARM1A.67,SPARM1A.74   
     & J,L,N,N_1,K                  ! WORK Loop counters
*D SPARM1A.85,SPARM1A.87   
     &,                 HT(1,N),LAI(1,N),SATCON
     &,                 CATCH_T(1,N),INFIL_T(1,N),Z0_T(1,N))
*I SPARM1A.93
        n_1=(n-npft-1)/nelev + 1
*D SPARM1A.96,SPARM1A.99   
          CATCH_T(L,N) = CATCH_NVG(N_1)
          INFIL_T(L,N) = INFIL_NVG(N_1)*SATCON(L)
          Z0_T(L,N) = Z0_NVG(N_1)
*D SPARM1A.103,SPARM1A.106  
      IF (NTILES .EQ. 1) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
        DO L=1,LAND_FIELD
          CATCH(L) = 0.0
          FZ0(L) = 0.0
          INFIL(L) = 0.0
          Z0(L) = 0.0
*D SPARM1A.108
*I SPARM1A.109   
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            FZ0(L) = FZ0(L) + FRAC(L,N) / (LOG(LB / Z0_T(L,N)))**2
          ENDDO
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
        ENDDO

        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH(L) = CATCH(L) + FRAC(L,N) * CATCH_T(L,N)
            INFIL(L) = INFIL(L) + FRAC(L,N) * INFIL_T(L,N)
          ENDDO
        ENDDO

        DO L=1,LAND_FIELD
!         Canopy capacity is average over non-lake surface types
          CATCH_TILE(L,1) = 0.
          frac_lake=0.
          do n=(lake_1-1)*nelev+1,lake_1*nelev
            frac_lake=frac_lake+frac(l,n)
          enddo
          IF (FRAC_LAKE .LT. 1.)
     &      CATCH_TILE(L,1) = CATCH(L) / (1. - FRAC_LAKE)
          INFIL_TILE(L,1) = INFIL(L)
          Z0_TILE(L,1) = Z0(L)
        ENDDO  
      ELSEIF (NTILES .EQ. NELEV*2) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
        DO L=1,LAND_FIELD
          CATCH(L) = 0.0
          FZ0(L) = 0.0
          INFIL(L) = 0.0
          Z0(L) = 0.0
          DO N=1,NTILES
            FZ0_TILE(L,N)=0.
            catch_tile(l,n)=0.
            infil_tile(l,n)=0.
            z0_tile(l,n)=0.
            k=mod(n-1,nelev)+1
            elev_frac(l,k)=0.
          ENDDO
        ENDDO

        DO L=1,LAND_FIELD
        FRAC_ICE(L)=0.
        DO N=NTYPE-NELEV+1,NTYPE
          FRAC_ICE(L)=FRAC_ICE(L)+FRAC(L,N)
        END DO
          DO N=1,NTYPE-NELEV
            k=mod(n-1,nelev)+1
            elev_frac(l,k)=elev_frac(l,k)+frac(l,n)
          END DO
        END DO
        DO N=1,NTYPE-NELEV
          k=mod(n-1,nelev)+1
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            if (FRAC_ICE(L).lt.1)
     &      FZ0(L) = FZ0(L) + (FRAC(L,N)/(1-FRAC_ICE(L))) / 
     &                        (LOG(LB / Z0_T(L,N)))**2
            if (elev_frac(L,k).gt.0)
     &      FZ0_TILE(L,K)= FZ0_TILE(L,K) + (FRAC(L,N)/elev_frac(l,k)) / 
     &                        (LOG(LB / Z0_T(L,N)))**2
          ENDDO
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
          DO N=1,NELEV
            Z0_TILE(L,N)=LB * EXP(-SQRT(1. / FZ0_TILE(L,N)))
          ENDDO
        ENDDO

        DO N=1,NTYPE-NELEV
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            if (FRAC_ICE(L).lt.1) THEN
            CATCH(L) = CATCH(L) + (FRAC(L,N)/(1-FRAC_ICE(L))) * 
     &                            CATCH_T(L,N)
            INFIL(L) = INFIL(L) + (FRAC(L,N)/(1-FRAC_ICE(L))) * 
     &                            INFIL_T(L,N)
            ENDIF
          ENDDO
        ENDDO

        DO L=LAND1,LAND1+LAND_PTS-1
!         Canopy capacity is average over non-lake surface types
          do n=1,nelev
            frac_lake_elev(l,n)=0.
            if (elev_frac(l,n).gt.0)
     &     frac_lake_elev(l,n)=frac(l,(lake_1-1)*nelev+n)/elev_frac(l,n)
          end do
        ENDDO

        do n=1,ntype-nelev
          k=mod(n-1,nelev)+1
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            IF (ELEV_FRAC(l,K) .gt. 0) THEN
              IF (FRAC_LAKE_ELEV(l,k) .LT. 1.)
     &           CATCH_TILE(L,K) = CATCH_TILE(L,K)+
     &    FRAC(L,N)*CATCH_T(L,N)/elev_frac(l,k)/(1.-frac_lake_elev(l,K))

            INFIL_TILE(L,k) =  INFIL_TILE(L,K)+
     &          FRAC(L,N)*INFIL_T(L,N)/elev_frac(l,k)
            ENDIF
          end do
        end do

        DO N=NTYPE-NELEV+1,NTYPE
          k=mod(n-1,nelev)+1
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH_TILE(L,nelev+k) = CATCH_T(L,N)
            INFIL_TILE(L,nelev+k) = INFIL_T(L,N)
            Z0_TILE(L,nelev+k) = Z0_T(L,N)
          ENDDO
        ENDDO

      ELSE
*D SPARM1A.111
! Copy surface-type arrays to tiles if separate tiles used
*D SPARM1A.113,SPARM1A.118  
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH_TILE(L,N) = CATCH_T(L,N)
            INFIL_TILE(L,N) = INFIL_T(L,N)
            Z0_TILE(L,N) = Z0_T(L,N)
          ENDDO
        ENDDO
*D SPARM1A.120,SPARM1A.141  
      ENDIF
*/-----
*DECLARE SPINDEX
*/-----
*D GDR7F405.77 
      PARAMETER(A_IXPTR_LEN = 69)
!
*/-----

*DECLARE STATMPT1
*/-----
*I STATMPT1.105
!
! JULES version 2 prognostics
      J_DEEP_ICE_TEMP(1) = SI(243,Sect_No,im_index)
      JSNOWDEPTH(1)     = SI(376,Sect_No,im_index) ! Snow depth on ground on    tiles (m)
      JRHO_SNOW_GRND(1) = SI(377,Sect_No,im_index) ! Snowpack bulk density (kg/ m3)
      JSNOW_CAN(1)      = SI(378,Sect_No,im_index) ! Snow on the canopy (kg/m2)
      JSNOW_SOIL_HTF(1) = SI(379,Sect_No,im_index) ! Surface heat flux beneath  snow (W/m2)
      JNSNOW(1)         = SI(380,Sect_No,im_index) ! Number of snow layers on   ground on tiles
      JDS(1)            = SI(381,Sect_No,im_index) ! Snow layer thickness (m)
      JSICE(1)          = SI(382,Sect_No,im_index) ! Snow layer ice mass on     tiles (Kg/m2)
      JSLIQ(1)          = SI(383,Sect_No,im_index) ! Snow layer liquid mass on  tiles (Kg/m2)
      JTSNOWLAYER(1)    = SI(384,Sect_No,im_index) ! Snow layer temperature (K)
      JRHO_SNOW(1)      = SI(385,Sect_No,im_index) ! Snow layer densities (kg/  m3)
      JRGRAINL(1)       = SI(386,Sect_No,im_index) ! Snow layer grain size on   tiles (microns)
!
      DO LEV=2,NTILES
        JSNOWDEPTH(LEV) = JSNOWDEPTH(LEV-1) + LAND_FIELD
        JRHO_SNOW_GRND(LEV) = JRHO_SNOW_GRND(LEV-1) + LAND_FIELD
        JSNOW_CAN(LEV) = JSNOW_CAN(LEV-1) + LAND_FIELD
        JSNOW_SOIL_HTF(LEV) = JSNOW_SOIL_HTF(LEV-1) + LAND_FIELD
        JNSNOW(LEV) = JNSNOW(LEV-1) + LAND_FIELD
      END DO
      DO LEV=2,NTILES*NSMAX
        JDS(LEV) = JDS(LEV-1) + LAND_FIELD
        JSICE(LEV) = JSICE(LEV-1) + LAND_FIELD
        JSLIQ(LEV) = JSLIQ(LEV-1) + LAND_FIELD
        JTSNOWLAYER(LEV) = JTSNOWLAYER(LEV-1) + LAND_FIELD
        JRHO_SNOW(LEV) = JRHO_SNOW(LEV-1) + LAND_FIELD
        JRGRAINL(LEV) = JRGRAINL(LEV-1) + LAND_FIELD
      END DO
*I STATMPT1.112
      DO LEV=2,ST_LEVELS
        J_DEEP_ICE_TEMP(LEV)=J_DEEP_ICE_TEMP(LEV-1)+LAND_FIELD
      ENDDO
!
*/-----

*I GSM3F404.7     

! Coastal tiling fields
      JFRAC_LAND     = SI(505,Sect_No,im_index)                         
      JTSTAR_LAND    = SI(506,Sect_No,im_index)                         
      JTSTAR_SEA     = SI(507,Sect_No,im_index)                         
      JTSTAR_SICE    = SI(508,Sect_No,im_index)                         
      JSICE_ALB      = SI(509,Sect_No,im_index)                         
      JLAND_ALB      = SI(510,Sect_No,im_index)                         
                             
*D ABX1F404.74,ABX1F404.78   
      JCAN_WATER_TYP= SI(229,Sect_No,im_index) ! Canopy water content   
C                                              ! on tiles               
      JCATCH_TYP    = SI(230,Sect_No,im_index) ! Canopy capacity on     
C                                              ! tiles                  
      JRGRAIN_TYP   = SI(231,Sect_No,im_index) ! Snow grain size on     
C                                              ! tiles                  
*I ABX1F404.81    
      JSNODEP_TYP   = SI(235,Sect_No,im_index) ! Tiled snow depth 
      JSNODEP_TYP_LC= SI(435,Sect_No,im_index) ! Tiled snow depth      
      JINFIL_TYP    = SI(236,Sect_No,im_index) ! Max tile infilt rate   
*/-----
*DECLARE STDEV7A
*/-----
*D STDEV7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D STDEV7A.34
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
*D STDEV7A.49
*D STDEV7A.52
     & FLANDG(P_FIELD)       ! IN Land fraction.                   
     &,BQ_1(P_FIELD)         ! IN Buoyancy parameter.                   
*D STDEV7A.91
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D STDEV7A.100,STDEV7A.101  
            T1_SD(I) = MAX ( 0.0 , 
     &          (1.-FLANDG(I))*1.93*FTL_1(I) / (RHOSTAR(I)*WS1) )   
            Q1_SD(I) = MAX ( 0.0 , 
     &          (1.-FLANDG(I))*1.93*FQW_1(I) / (RHOSTAR(I)*WS1) )   
*D STDEV7A.118,STDEV7A.119  
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,FLAND,         
     & BQ_1,BT_1,FQW_1,FTL_1,RHOKM_1,RHOSTAR,VSHR,Z0M,Z1_TQ,TILE_FRAC,  
*D STDEV7A.136,STDEV7A.137
     & FLAND(LAND_FIELD)                                                
     &,BQ_1(LAND_FIELD)         ! IN Buoyancy parameter.                   
     &,BT_1(LAND_FIELD)         ! IN Buoyancy parameter.                   
*D STDEV7A.145   
     &,Z1_TQ(LAND_FIELD)        ! IN Height of lowest tq level.
     &,TILE_FRAC(LAND_FIELD) ! IN Tile fraction.                        
*D STDEV7A.178,STDEV7A.179
        VSF1_CUBED = 1.25*G*(Z1_TQ(L) + Z0M(L)) *
     &             ( BT_1(L)*FTL_1(L) + BQ_1(L)*FQW_1(L) )/RHOSTAR(I)
*D STDEV7A.182,STDEV7A.188  
          T1_SD(I) = T1_SD(I) + MAX ( 0.0 ,                             
     &      FLAND(L)*TILE_FRAC(L)*1.93*FTL_1(L) / (RHOSTAR(I)*WS1) )    
          Q1_SD(I) = Q1_SD(I) + MAX ( 0.0 ,                             
     &      FLAND(L)*TILE_FRAC(L)*1.93*FQW_1(L) / (RHOSTAR(I)*WS1) )    
*/-----
*DECLARE SURFHY7A
*/-----
*D SURFHY7A.55,SURFHY7A.56   
      SUBROUTINE SURF_HYD (NPNTS,NTILES,NELEV,TILE_PTS,TILE_INDEX,            
*D SURFHY7A.58
     &                     LS_RAIN,MELT_TILE,SNOW_MELT,TIMESTEP,        
*D SURFHY7A.60
     &                     CAN_WCNT_GB,DSMC_DT,SURF_ROFF,TOT_TFALL
     &                     ,tfall_tile,nsnow,l_essery_snow
     &                     )
*D SURFHY7A.66,SURFHY7A.68   
     &,NTILES               ! IN Number of tiles.                       
     &,NELEV
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
     &,TILE_INDEX(NPNTS,NTILES)                                         
*D SURFHY7A.70,SURFHY7A.71   
*D SURFHY7A.74,SURFHY7A.79   
     & CAN_CPY(NPNTS,NTILES)! IN Canopy capacity for land tiles (kg/m2).
     &,E_CANOPY(NPNTS,NTILES)!IN Canopy evaporation (kg/m2/s).
     &,FRAC(NPNTS,NTILES)   ! IN Tile fractions.                        
     &,INFIL(NPNTS,NTILES)  ! IN Infiltration rate (kg/m2/s).           
*D SURFHY7A.82,SURFHY7A.83   
     &,MELT_TILE(NPNTS,NTILES)
!                           ! IN Snow melt on tiles (kg/m2/s).        
     &,SNOW_MELT(NPNTS)     ! IN GBM snow melt (kg/m2/s).               
*I SURFHY7A.84
     &,NSNOW(NPNTS,NTILES)

      LOGICAL L_ESSERY_SNOW,SOIL
*D SURFHY7A.87,SURFHY7A.89   
     & CAN_WCNT(NPNTS,NTILES)!INOUT Tile canopy water contents (kg/m2).
     &,TFALL_TILE(NPNTS,NTILES)     ! OUT Cumulative canopy throughfall
*I SURFHY7A.101
     &,ICE_FRAC(NPNTS)
*I  SURFHY7A.115
        ice_frac(i) = 0.0
        DO N=1,NTILES                                                   
          TFALL_TILE(I,N) =0.
        END DO
        DO N=NTILES-nelev+1,ntiles                                                   
          ice_frac(i) = ice_frac(i) + frac(i,n)
        END DO
*D SURFHY7A.118,SURFHY7A.135  
*D SURFHY7A.137
      DO N=1,NTILES                                                   
*D SURFHY7A.146
! THIS IS ALL JUST FOR THE SOIL TILES NOW?!
      DO N=1,(NTILES-NELEV)

        SOIL=.TRUE.

! Surface runoff of snowmelt, assumed to cover 100% of tile
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,            
     &                CAN_CPY(1,N),CAN_CPY(1,N),INFIL(1,N),   
     &                MELT_TILE(1,N),FRAC(1,N),TIMESTEP,               
     &                SURF_ROFF,NSNOW(1,N),L_ESSERY_SNOW)

*D SURFHY7A.162
     &              CAN_WCNT(1,N),TOT_TFALL,TFALL_TILE(1,N),SOIL)
        IF (SOIL) THEN
*D SURFHY7A.166
     &                SURF_ROFF,NSNOW(1,N),L_ESSERY_SNOW)
        ENDIF
*D  SURFHY7A.172
     &              CAN_WCNT(1,N),TOT_TFALL,TFALL_TILE(1,N),SOIL)
        IF (SOIL) THEN
*D SURFHY7A.176
     &                SURF_ROFF,NSNOW(1,N),L_ESSERY_SNOW)
        ENDIF
*D SURFHY7A.182
     &              CAN_WCNT(1,N),TOT_TFALL,TFALL_TILE(1,N),SOIL)
        IF (SOIL) THEN
*D SURFHY7A.186
     &                SURF_ROFF,NSNOW(1,N),L_ESSERY_SNOW)
        ENDIF
*D SURFHY7A.188
        DO J=1,TILE_PTS(N)                                              
          I = TILE_INDEX(J,N)                                       
*D SURFHY7A.194,SURFHY7A.196
      DO I=1,NPNTS
!        DSMC_DT(I) = TOT_TFALL(I) + SNOW_MELT(I) - SURF_ROFF(I)
! these GBMs have only been done on SOIL points, so normalise by GBM
! soil fraction to get the correct soil gbm
        if (ice_frac(i) .lt. 1) then !really should be on these points!
          TOT_TFALL(I)=TOT_TFALL(I)/(1-ICE_FRAC(i))
          SURF_ROFF(I)=SURF_ROFF(I)/(1-ICE_FRAC(i))
        endif
        DSMC_DT(I) = TOT_TFALL(I) + SNOW_MELT(I) - SURF_ROFF(I)

      ENDDO
*/-----
*DECLARE SWRAD3A
*/-----
*I ADB2F404.1504  
!       4.6             10-05-98                Land flag passed to     
!                                               FLUX_CALC.              
!                                               (J. M. Edwards)         
*D ARE2F404.239,ARE2F404.241  
     &   , LAND_ALBEDO, L_CTILE                                         
     &   , LAND_ALB, SICE_ALB, FLANDG                                   
     &   , OPEN_SEA_ALBEDO, ICE_FRACTION, LAND, LAND0P5, LYING_SNOW
!                       MOSES II flag and array dimension               
     &   , L_MOSES_II, SAL_DIM                                          
*D SWRAD3A.46
     &   , FLUX_BELOW_690NM_SURF,  FL_SOLID_BELOW_690NM_SURF            
     &   , FL_SEA_BELOW_690NM_SURF, L_FLUX_BELOW_690NM_SURF             
*I SWRAD3A.74    
     &   , SURF_DOWN_SW                                                 
*D ARE2F404.242,ARE2F404.243  
     &   , LAND0P5(NPD_FIELD)                                           
!             LAND MASK (TRUE if land fraction >0.5)                    
c     &   , SEA(NPD_FIELD)                                              
c!             SEA MASK                                           
*D ARE2F404.246
!             DIMENSION FOR LAND_ALBEDO (MOSES II)                      
*D SWRAD3A.255,ADB1F401.1045 
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX            
     &   , LAND_ALBEDO(SAL_DIM,4)                                       
!             MOSES II LAND SURFACE ALBEDO FIELDS                       
     &   , LAND_ALB(NPD_FIELD)                                          
!             SURFACE ALBEDO OF LAND                                    
     &   , SICE_ALB(NPD_FIELD)                                          
!             SURFACE ALBEDO OF SEA-ICE                                 
     &   , FLANDG(NPD_FIELD)                                            
!             Fractional land                                           
     &   , FLANDG_G(NPD_PROFILE)                                        
!             Gathered Fractional land                                  
     &   , ICE_FRACTION_G(NPD_PROFILE)                                  
!             Gathered Fractional sea-ice                               
*I SWRAD3A.276   
!             WEIGHTED BY (OPEN SEA)/(TOTAL SEA) FRACTION               
*I SWRAD3A.278   
     &   , SURF_DOWN_SW(NPD_FIELD, 4)                                   
!             SURFACE DOWNWARD SHORTWAVE RADIATION COMPONENTS           
!             (*,1) - DIRECT BEAM VISIBLE                               
!             (*,2) - DIFFUSE VISIBLE                                   
!             (*,3) - DIRECT BEAM NEAR-IR                               
!             (*,4) - DIFFUSE NEAR-IR                                   
*I ADB1F401.1048  
     &   , L_MOSES_II                                                   
!             SURFACE SW FLUXES REQUIRED FOR MOSES II                   
     &   , L_CTILE                                                      
!             SWITCH FOR COASTAL TILING                                 
*I ADB1F401.1050  
!             NB: ONLY USED FOR NON MOSESII RUNS.                       
     &   , FL_SOLID_BELOW_690NM_SURF(NPD_FIELD)                         
!             SOLID SFC NET SURFACE FLUX BELOW 690NM                    
     &   , FL_SEA_BELOW_690NM_SURF(NPD_FIELD)                           
!             OPEN SEA NET SURFACE FLUX BELOW 690NM                     
!             (AT POINTS WHERE THERE THIS IS SEA-ICE THIS IS            
!             WEIGHTED BY THE FRACTION OF OPEN SEA.)                    
*I SWRAD3A.432   
     &   , LAND0P5_G(NPD_PROFILE)                                       
!             GATHERED LAND MASK (TRUE if land fraction >0.5)           
*D SWRAD3A.485
!             GATHERED GRID BOX MEAN DOWNWARD SURFACE FLUX              
!             BELOW 690NM                                               
     &   , FL_SEA_BELOW_690NM_SURF_G(NPD_PROFILE)                       
!             GATHERED OPEN SEA NET SURFACE FLUX BELOW 690NM            
*I SWRAD3A.486   
!     TEMPORARY FIELDS ASSOCIATED WITH 690NM FLUX OVER                  
!     MEAN SOLID SURF.                                                  
      REAL                                                              
     &     ALBEDOSOLID                                                  
     &    ,FRACSOLID                                                    
     &    ,FLUXSOLID                                                    
!                                                                       
     &   , SURF_VIS_DIR_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX        
     &   , SURF_VIS_DIF_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIFFUSE VISIBLE FLUX            
     &   , SURF_NIR_DIR_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX  
     &   , SURF_NIR_DIF_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX      
!                                                                       
*D ADB1F401.1069,ADB1F401.1070 
     &     N_FRAC_SOL_POINT                                             
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
*D SWRAD3A.530,SWRAD3A.531  
*D ARE2F404.252
      IF ( L_FLUX_BELOW_690NM_SURF .OR. L_MOSES_II ) THEN               
*D ARE2F404.264,ARE2F404.265  
     &   , L_MICROPHYSICS, L_MOSES_II, SAL_DIM, L_CTILE                 
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               
     &   , LAND_ALB, SICE_ALB                                           
     &   , FLANDG, ICE_FRACTION                                         
     &   , LAND_ALBEDO, WEIGHT_690NM                                    
*D SWRAD3A.545
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G                          
*D ADB2F404.1571
     &   , PSTAR, DUMMY, DUMMY, DUMMY
     &   , AB, BB, AC, BC                                 
*D ADB2F404.1573
     &   , P, T, DUMMY, DUMMY, DUMMY, DUMMY, D_MASS                     
*D ADB1F402.725
     &      , LAND0P5, LYING_SNOW, PSTAR, AB, BB, TRINDX                
*D AYY1F404.372
     &   , L_CLOUD_WATER_PARTITION,  LAND0P5_G                          
*D ADB2F404.1604
     &   , P, T, DUMMY, DUMMY, DUMMY, DUMMY, D_MASS                     
*D ADB1F401.1105,SWRAD3A.760  
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION_G           
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G, FLANDG_G                
*D SWRAD3A.764
     &   , FLUX_BELOW_690NM_SURF_G, FL_SEA_BELOW_690NM_SURF_G           
     &   , L_MOSES_II, L_CTILE                                          
     &   , SURF_VIS_DIR_G,SURF_VIS_DIF_G,SURF_NIR_DIR_G,SURF_NIR_DIF_G  
*I SWRAD3A.795   
      IF (L_MOSES_II) THEN                                              
         DO I=1, 4                                                      
            CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_SW(1, I))              
         ENDDO                                                          
      ENDIF                                                             
*I SWRAD3A.891   
        IF(L_CTILE)THEN                                                 
           CALL R2_ZERO_1D(N_PROFILE, FL_SOLID_BELOW_690NM_SURF)        
           CALL R2_ZERO_1D(N_PROFILE, FL_SEA_BELOW_690NM_SURF)          
          DO L=1, NLIT                                                  
            IF (FLANDG(LIST(L)).LT.1.0) THEN                            
              FL_SEA_BELOW_690NM_SURF(LIST(L))                          
     &          =FL_SEA_BELOW_690NM_SURF_G(L)                           
     &         *(1.0E+00-ICE_FRACTION(LIST(L)))                         
            ENDIF                                                       
            FRACSOLID=FLANDG(LIST(L))                                   
     &        +(1.-FLANDG(LIST(L)))*ICE_FRACTION(LIST(L))               
            IF(FRACSOLID.GT.0.0)THEN                                    
              FL_SOLID_BELOW_690NM_SURF(LIST(L))=                       
     &         (FLUX_BELOW_690NM_SURF_G(L)-                             
     &         (1.-FLANDG(LIST(L)))*FL_SEA_BELOW_690NM_SURF(LIST(L)))   
     &         /FRACSOLID                                               
            ENDIF                                                       
          ENDDO                                                         
        ELSE                                                            
*I ADB1F401.1111  
                                                                        
*I SWRAD3A.895   
        ENDIF                                                           
                                                                        
      ENDIF                                                             
!                                                                       
!                                                                       
!     COMPONENTS OF DOWNWARD FLUX AT THE SURFACE FOR MOSES II           
!                                                                       
      IF (L_MOSES_II) THEN                                              
        DO L=1, NLIT                                                    
           SURF_DOWN_SW(LIST(L),1) = SURF_VIS_DIR_G(L)                  
           SURF_DOWN_SW(LIST(L),2) = SURF_VIS_DIF_G(L)                  
           SURF_DOWN_SW(LIST(L),3) = SURF_NIR_DIR_G(L)                  
           SURF_DOWN_SW(LIST(L),4) = SURF_NIR_DIF_G(L)                  
        ENDDO                                                           
*I SWRAD3A.957   
! Weight open sea flux with open sea fraction over total sea            
      IF(L_CTILE)THEN                                                   
        DO L=1, NLIT                                                    
          IF (FLANDG(LIST(L)).LT.1.0) THEN                              
            SWSEA(LIST(L))=(1.0E+00-ICE_FRACTION(LIST(L)))              
     &         *SEA_FLUX_G(L)                                           
            SWOUT(LIST(L), 1)=SWOUT(LIST(L), 1)                         
     &        -(1.0E+00-FLANDG(LIST(L)))*SWSEA(LIST(L))                 
           ELSE                                                         
            SWSEA(LIST(L))=0.0                                          
           ENDIF                                                        
         ENDDO                                                          
       ENDIF                                                            
       IF(.NOT.L_CTILE)THEN                                             
*I SWRAD3A.965   
      ENDIF                                                             
*D AJS1F401.1424
!     TOTAL DOWNWARD FLUX OF PHOTOSYTHETICALLY ACTIVE RADIATION         
*I AJS1F401.1428  
        IF (L_MOSES_II) THEN                                            
          DO L=1, NLIT                                                  
            SWOUT(LIST(L),NLEVS+2) = SURF_VIS_DIR_G(L) +                
     &                               SURF_VIS_DIF_G(L)                  
          ENDDO                                                         
        ELSE       
*D AJS1F401.1430,AJS1F401.1432 
             IF(LAND(L))THEN                                            
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             
     &           (1 - LAND_ALB(L))                                      
              ELSE                                                      
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             
     &           (1 - SICE_ALB(L))                                      
              ENDIF                                                     
                                                                        
          ENDDO                                                         
        ENDIF                                                           
*I SWRAD3A.976   
!                                                                       
!     DIVIDE SURFACE DOWNWARD SW COMPONENTS BY COSINE OF SOLAR ZENITH   
!     ANGLE FOR MOSES II.                                               
      IF (L_MOSES_II) THEN                                              
        DO I=1, 4                                                       
          DO L=1, N_PROFILE                                             
            SURF_DOWN_SW(L, I) = SURF_DOWN_SW(L, I) /                   
     &                          (COSZIN(L)*LIT(L)+TOL_MACHINE)          
          ENDDO                                                         
        ENDDO                                                           
      ENDIF                                                             
*D ARE2F404.273,ARE2F404.274  
     &   , L_MICROPHYSICS, L_MOSES_II, SAL_DIM, L_CTILE                 
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               
     &   , LAND_ALB, SICE_ALB                                           
     &   , FLANDG, ICE_FRACTION                                         
     &   , LAND_ALBEDO, WEIGHT_690NM                                    
*D SWRAD3A.1008
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                              
*D ARE2F404.276
!             DIMENSION OF LAND_ALBEDO                                  
*I SWRAD3A.1057  
     &   , LAND0P5(NPD_FIELD)                                           
!             LAND MASK (TRUE if land fraction >0.5)                    
*D SWRAD3A.1061,ARE2F404.280  
     &   , FLANDG(NPD_FIELD)                                            
     &   , FLANDG_G(NPD_PROFILE)                                        
!             GATHERED LAND FRACTION                                    
     &   , ICE_FRACTION_G(NPD_PROFILE)                                  
!             GATHERED SEA-ICE FRACTION IN SEA PORTION OF GRID BOX      
     &   , LAND_ALB(NPD_FIELD)                                          
     &   , SICE_ALB(NPD_FIELD)    
     &   , LAND_ALBEDO(SAL_DIM,4)          
!             MOSES II LAND SURFACE ALBEDO FIELDS                       
*D SWRAD3A.1064
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX            
*D ARE2F404.283,ARE2F404.284  
     &   , L_MOSES_II                                                   
     &   , L_CTILE                                                      
!             FLAG FOR COASTAL TILING                                   
*I SWRAD3A.1084  
     &   , LAND0P5_G(NPD_PROFILE)                                       
!             GATHERED LAND MASK (TRUE if land fraction >0.5)           
*I SWRAD3A.1113  
            LAND0P5_G(L)=LAND0P5(LIST(L))                               
            FLANDG_G(L)=FLANDG(LIST(L))                                 
            ICE_FRACTION_G(L)=ICE_FRACTION(LIST(L))                     
*D SWRAD3A.1119,SWRAD3A.1121 
!     THERE IS A COMBINATION OF OPEN SEA, SEA-ICE AND LAND. SEPARATE    
!     ALBEDOS ARE PROVIDED FOR FOR OPEN SEA. BAND-DEPENDENT COPIES      
!     OF THE ALBEDOS MUST BE MADE FOR CALCULATING COUPLING FLUXES.      
*D SWRAD3A.1128
            ALBEDO_SEA_DIR(L, I)=0.0E+00                                
            ALBEDO_SEA_DIFF(L, I)=0.0E+00                               
            ALBEDO_FIELD_DIFF(L, I)=0.0E+00                             
            ALBEDO_FIELD_DIR(L, I)=0.0E+00                              
                                                                        
            IF (FLANDG(LIST(L)).LT.1.0) THEN                            
*D SWRAD3A.1130
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              
*D SWRAD3A.1134
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              
*D SWRAD3A.1139,ARE2F404.285  
            ENDIF                                                       
            IF (FLANDG(LIST(L)).GT.0.0) THEN                            
               IF ( L_MOSES_II ) THEN                                   
              IF ( L_CTILE ) THEN                                       
                ALBEDO_FIELD_DIFF(L,I) =                                
     &            (1.-FLANDG_G(L))*ALBEDO_FIELD_DIFF(L, I) +            
     &            FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)   
     &               + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4))   
                ALBEDO_FIELD_DIR(L,I) =                                 
     &            (1.-FLANDG_G(L))*ALBEDO_FIELD_DIR(L, I) +             
     &            FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)   
     &               + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3))   
               ELSE                                                     
*D ARE2F404.287,ARE2F404.288  
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4)
*D ARE2F404.290,ARE2F404.291  
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3)
               ENDIF                                                    
*D SWRAD3A.1140,SWRAD3A.1141 
! For non MOSES_II cannot have coastal tiling, therefore                
! must be completely land:                                              
               ALBEDO_FIELD_DIFF(L, I)=LAND_ALB(LIST(L))                
               ALBEDO_FIELD_DIR(L, I)=LAND_ALB(LIST(L))                 
*D SWRAD3A.1142,SWRAD3A.1143 
*/-----
*DECLARE TRIF
*/-----
*D TRIF.6,TRIF.62
*CALL C_LAND_CC

      INTEGER
     + C3(NPFT_1)                   ! 1 for C3 Plants, 0 for C4 Plants.
     +,CROP(NPFT_1)                 ! 1 for crop type, 0 for non-crop.
     +,ORIENT(NPFT_1)               ! 1 for horizontal, 0 for spherical.

      REAL
     + ALPHA(NPFT_1)                ! Quantum efficiency
C                                 ! (mol CO2/mol PAR photons).
     +,ALNIR(NPFT_1)                ! Leaf reflection coefficient for
C                                 ! near infra-red.                  
     +,ALPAR(NPFT_1)                ! Leaf reflection coefficient for
C                                 ! PAR.                  
     +,A_WL(NPFT_1)                 ! Allometric coefficient relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index (kg C/m2).
     +,A_WS(NPFT_1)                 ! Woody biomass as a multiple of
C                                 ! live stem biomass.
     +,B_WL(NPFT_1)                 ! Allometric exponent relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index.
     +,DGL_DM(NPFT_1)               ! Rate of change of leaf turnover
C                                 ! rate with moisture availability.
     +,DGL_DT(NPFT_1)               ! Rate of change of leaf turnover
C                                 ! rate with temperature (/K)
     +,DQCRIT(NPFT_1)               ! Critical humidity deficit
C                                 ! (kg H2O/kg air).
     +,ETA_SL(NPFT_1)               ! Live stemwood coefficient
C                                 ! (kg C/m/LAI).
     +,FSMC_OF(NPFT_1)              ! Moisture availability below
C                                 ! which leaves are dropped.
     +,GLMIN(NPFT_1)                ! Minimum leaf conductance for H2O
     +,G_AREA(NPFT_1)               ! Disturbance rate (/360days).
     +,G_GROW(NPFT_1)               ! Rate of leaf growth (/360days).
     +,G_LEAF_0(NPFT_1)             ! Minimum turnover rate for leaves
!                                 ! (/360days).
     +,G_ROOT(NPFT_1)               ! Turnover rate for root biomass
!                                 ! (/360days).
     +,G_WOOD(NPFT_1)               ! Turnover rate for woody biomass
!                                 ! (/360days).
     +,KPAR(NPFT_1)                 ! PAR Extinction coefficient
C                                 ! (m2 leaf/m2 ground).
     +,LAI_MAX(NPFT_1)              ! Maximum projected LAI.
C                                 ! (kg N/kg C).
     +,NR_NL(NPFT_1)                ! Ratio of root nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,NS_NL(NPFT_1)                ! Ratio of stem nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,OMEGA(NPFT_1)                ! Leaf scattering coefficient
C                                 ! for PAR.
     +,OMNIR(NPFT_1)                ! Leaf scattering coefficient for
C                                 ! near infra-red.
     +,SIGL(NPFT_1)                 ! Specific density of leaf carbon
C                                 ! (kg C/m2 leaf).
     +,TLEAF_OF(NPFT_1)             ! Temperature below which leaves are
C                                 ! dropped.
C                                 ! dropped.


*I TRIF.67    
      DATA CROP    /      0,     0,     1,     1,     0 /
      DATA ORIENT  /      0,     0,     0,     0,     0 /  
*I TRIF.68    
      DATA ALNIR   /   0.45,  0.35,  0.58,  0.58,  0.58 /               
      DATA ALPAR   /   0.10,  0.07,  0.10,  0.10,  0.10 /
*D ABX1F405.1732,TRIF.73   
      DATA DGL_DM  /    0.0,   0.0,   0.0,   0.0,   0.0 /               
      DATA DGL_DT  /    9.0,   9.0,   0.0,   0.0,   9.0 /               
*D ABX1F405.1733
      DATA FSMC_OF /   0.00,  0.00,  0.00,  0.00,  0.00 /               
*D TRIF.76
*D TRIF.79
      DATA G_AREA  /  0.005, 0.004,  0.25,  0.25,  0.05 /               
*D ABX1F405.1734,TRIF.82   
      DATA G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25 /               
      DATA G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25 /               
*D TRIF.85,TRIF.87   
      DATA LAI_MAX /   9.00,  9.00,  4.00,  4.00,  4.00 /               
*I TRIF.90    
      DATA OMNIR   /   0.70,  0.45,  0.83,  0.83,  0.83 /             
*D TRIF.91    
*D TRIF.92,ABX1F405.1737 
      DATA SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500 /               
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15 /               
*D ABX1F405.1738
*D ABX1F405.1739

*DECLARE TRIFD2A
*/-----
*D ABX1F405.1591
     &,                   FRAC_VS,FRAC_AGRIC,G_LEAF,NPP,RESP_S,RESP_W   
*D TRIFD2A.61,ABX1F405.1600 
     &,FRAC_AGRIC(LAND_FIELD)     ! IN Fraction of agriculture.    
*I TRIFD2A.89    
     &,FRAC_FLUX                  ! WORK PFT fraction to be used
C                                 !      in the calculation of
C                                 !      the gridbox mean fluxes.
*D ABX1F405.1618
     &,           C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI_BAL,PC_S     
*D TRIFD2A.165
C type                                              
*D TRIFD2A.172,TRIFD2A.174  
          FRAC_FLUX=FRAC(L,N)-(1.0-FORW)*DFRAC(L,N)
          LIT_C(L,N) = NPP(L,N)-GAMMA/FRAC_FLUX*(C_VEG(L,N)*FRAC(L,N)
     &               -(C_VEG(L,N)-DCVEG(L,N))*(FRAC(L,N)-DFRAC(L,N)))
          LIT_C_T(L) = LIT_C_T(L)+FRAC_FLUX*LIT_C(L,N)                  
*/-----
*DECLARE TYPPTRA
*/-----
*I AJS1F401.5
     &             J_DEEP_ICE_TEMP(ST_LEVELS),
*I TYPPTRA.105
      INTEGER
     &   JSNOWDEPTH(NTILES),
     &   JRHO_SNOW_GRND(NTILES),
     &   JSNOW_CAN(NTILES),
     &   JSNOW_SOIL_HTF(NTILES),
     &   JNSNOW(NTILES),
     &   JDS(NTILES*NSMAX),
     &   JSICE(NTILES*NSMAX),
     &   JSLIQ(NTILES*NSMAX),
     &   JTSNOWLAYER(NTILES*NSMAX),
     &   JRHO_SNOW(NTILES*NSMAX),
     &   JRGRAINL(NTILES*NSMAX)
!
*/-----

*D ABX1F404.24
*I TYPPTRA.36    
     &       JFRAC_LAND,             ! Land fraction in grid box
     &       JTSTAR_LAND,            ! Land surface temperature
     &       JTSTAR_SEA,             ! Sea surface temperature
     &       JTSTAR_SICE,            ! Sea-ice surface temperature
     &       JSICE_ALB,              ! Sea-ice albedo
     &       JLAND_ALB,              ! Mean land albedo
*D ABX1F404.41,ABX1F404.43   
     &       JCAN_WATER_TYP,         ! Canopy water content on tiles    
     &       JCATCH_TYP,             ! Canopy capacity on tiles         
     &       JINFIL_TYP,             ! Max infiltration rate on tiles   
     &       JRGRAIN_TYP,            ! Snow grain size on tiles         
     &       JSNODEP_TYP,            ! Snow depth on tiles
     &       JSNODEP_TYP_LC,         ! Snow depth on tiles (dump acc)
*D ABX1F404.46,AJS1F401.24   
     &  JSNSOOT, JTSTAR_ANOM, 
     &  JFRAC_LAND, JTSTAR_LAND, JTSTAR_SEA, JTSTAR_SICE,
     &  JSICE_ALB, JLAND_ALB,
     &  JZH, JZ0, JLAND, JICE_FRACTION,                    
*D ABX1F404.50,ABX1F404.51   
     &  JRSP_S_ACC, JTSNOW, JCAN_WATER_TYP, JCATCH_TYP, JINFIL_TYP,     
     &  JRGRAIN_TYP, JSNODEP_TYP, JSNODEP_TYP_LC, JTSTAR_TYP, JZ0_TYP                   
*/-----
*DECLARE TYPSIZE
*/-----
*I TYPSIZE.25    
     &      ,NTILES               ! IN: No of land surface tiles
*I TYPSIZE.127   
     & NTILES,   
*I TYPSIZE.83
!
      INTEGER NSMAX    ! IN: Max number of snow layers
!
*I AJX1F404.13
     & ,NSMAX
!
*/-----
*DECLARE U_MODEL1
*/-----
*D ARE2F404.530
     &    RADINCS ( (P_FIELDDA*(P_LEVELSDA+2+NTILES)+511)/512*512*2 )
     &   ,DOWNWELLING_LW(P_FIELDDA,BL_LEVELS)
*I GGH3F401.44
     &                  DOWNWELLING_LW,
*/-----
*DECLARE UPANCIL1
*/-----
*I GDG0F401.1492  
     &                D1(JFRAC_LAND),                                   
     &                D1(JTSTAR_LAND),D1(JTSTAR_SEA),                   
     &                D1(JTSTAR_SICE),                                  
*/-----
*DECLARE UMINDEX1
*/-----
*/-----
*I UMINDEX1.235
*CALL NSTYPES
*D GDR7F405.74
!
      A_IXPTR(58)=A_IXPTR(57) + P_LEVELs       ! JSNOWDEPTH
      A_IXPTR(59)=A_IXPTR(58) + NTYPE          ! JRHO_SNOW_GRND
      A_IXPTR(60)=A_IXPTR(59) + NTYPE          ! JSNOW_CAN
      A_IXPTR(61)=A_IXPTR(60) + NTYPE          ! JSNOW_SOIL_HTF
      A_IXPTR(62)=A_IXPTR(61) + NTYPE          ! JNSNOW
      A_IXPTR(63)=A_IXPTR(62) + NTYPE          ! JDS
      A_IXPTR(64)=A_IXPTR(63) + (NSMAX*NTYPE)  ! JSICE
      A_IXPTR(65)=A_IXPTR(64) + (NSMAX*NTYPE)  ! JSLIQ
      A_IXPTR(66)=A_IXPTR(65) + (NSMAX*NTYPE)  ! JTSNOWLAYER
      A_IXPTR(67)=A_IXPTR(66) + (NSMAX*NTYPE)  ! JRHO_SNOW
      A_IXPTR(68)=A_IXPTR(67) + (NSMAX*NTYPE)  ! JRGRAINL
      A_IXPTR(69)=A_IXPTR(68) + (NSMAX*NTYPE)  ! J_DEEP_ICE_TEMP
      A_SPPTR_LEN = A_IXPTR(69) + ST_LEVELS
!
*/-----
*DECLARE VEG1A
*/-----
*D ABX3F405.29
     &,              LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH 
*D VEG1A.26
     &,              ATIMESTEP,SATCON                                
*D ABX1F405.1337
     &,              CATCH_T,INFIL_T,Z0_T                    
*I VEG1A.63    
     &,NTILES                ! IN Number of land-surface tiles.
*D VEG1A.81,VEG1A.82   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).        
*D VEG1A.90,VEG1A.93   
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
!                                   !     (kg/m2).  
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG1A.94,VEG1A.97   
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D VEG1A.114,VEG1A.115  
*D ABX1F405.1348,VEG1A.132  
      DO N=1,NTILES
*I VEG1A.134   
          INFIL_T(L,N)=0.0  
          Z0_T(L,N)=0.0 
*D VEG1A.193,VEG1A.203  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX  
     &,           FRAC,HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)
*D ABX1F405.1352
*D ABX1F405.1354
*/-----
*DECLARE VEG2A
*/-----
*D ABX3F405.74
     &,              LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH
*D VEG2A.29
     &,              ATIMESTEP,FRACA,SATCON                    
*D VEG2A.33
     &,              CATCH_T,INFIL_T,Z0_T                    
*I VEG2A.69    
     &,NTILES                ! IN Number of land-surface tiles.         
*D VEG2A.94,VEG2A.97   
     & ATIMESTEP                    ! IN Atmospheric timestep (s). 
     &,FRACA(LAND_FIELD)            ! IN Fraction of agriculture. 
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
C                                   !    conductivity (kg/m2/s).        
*D VEG2A.114,VEG2A.121  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
C                                   !     (kg/m2). 
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D VEG2A.150,ABX1F405.1445 
     &,RATIO                        ! WORK Ratio of fractional 
!                                   !      coverage before to that 
!                                   !      after TRIFFID. 
     &,DCS(LAND_FIELD)              ! WORK Change in soil carbon
!                                   !      (kg C/m2).
     &,FRAC_AGRIC(LAND_FIELD)       ! WORK Fraction of agriculture as
!                                   !      seen by TRIFFID.          
     &,FRAC_OLD(LAND_FIELD,NTYPE)   ! WORK Fractions of surface types 
!                                   !      before the call to TRIFFID.
     &,FRAC_VS(LAND_FIELD)          ! WORK Total fraction of gridbox    
!                                   !      covered by veg or soil.      
*D ABX1F405.1454,VEG2A.165  
*D VEG2A.170,VEG2A.173  

      LOGICAL
     & AGRIC                        ! .T. for TRIFFID to see agriculture
!                                   ! .F. for natural vegetation.
      PARAMETER (AGRIC=.TRUE.)
*D VEG2A.190
      DO N=1,NTILES                                                     
*I ABX1F405.1459  
          CATCH_T(L,N)=0.0
          INFIL_T(L,N)=0.0
*D VEG2A.196
*D VEG2A.198,VEG2A.207  
        RESP_S_DR(L)=0.0
*D ABX1F405.1481
*D ABX1F405.1483
*I ABX1F405.1462  
        FRAC_AGRIC(L) = 0.0
*D ABX3F405.96,ABX3F405.99   
*D ABX3F405.113
      CALL SWAPB_LAND(FRACA,LAND_FIELD,P_FIELD,                        
*D VEG2A.313
C Define the agricultural regions.
*I VEG2A.314   
        IF (AGRIC) THEN
*D VEG2A.316
            FRAC_AGRIC(L)=FRACA(L)
*I VEG2A.317   
        ENDIF
*D ABX1F405.1490,ABX1F405.1492 
C-----------------------------------------------------------------------
C Take copies of TRIFFID input variables for output as diagnostics.     
C-----------------------------------------------------------------------
*I ABX1F405.1497  
            FRAC_OLD(L,N)=FRAC(L,N)
*I ABX1F405.1501  
          DCS(L)=CS(L)
*D VEG2A.327
          FORW=0.0        
*D ABX1F405.1511
     &,                 FRAC_VS,FRAC_AGRIC,G_LEAF_DR,NPP_DR,RESP_S_DR   
*I ABX1F405.1522  

*D VEG2A.337
C Reset the accumulation fluxes to zero in equilibrium mode.
*I VEG2A.338   
        IF (L_TRIF_EQ) THEN

*D VEG2A.346
              RESP_W_AC(L,N)=0.0      
*I VEG2A.349   
C-----------------------------------------------------------------------
C Reset the accumulation fluxes to the TRIFFID-diagnosed corrections
C in dynamic mode. Such corrections will be typically associated with 
C the total depletion of a carbon reservoir or with the maintenance
C of the plant seed fraction.
C-----------------------------------------------------------------------
        ELSE

          DO L=LAND1,LAND1+LAND_PTS-1                                   
            DCS(L)=CS(L)-DCS(L)
            RESP_S_DR(L)=LIT_C_MN(L)-GAMMA*DCS(L)
            RESP_S_AC(L)=(RESP_S_DR(L)-RESP_S_DR_OUT(L))/GAM_TRIF   
          ENDDO                                                         
                                                                        
          DO N=1,NPFT                                                   
            DO J=1,TILE_PTS(N)                                          
              L=TILE_INDEX(J,N)                                         
              RATIO=FRAC_OLD(L,N)/FRAC(L,N)
              NPP_AC(L,N)=RATIO*(NPP_DR(L,N)-NPP_DR_OUT(L,N))/GAM_TRIF  
              RESP_W_AC(L,N)=RATIO*(RESP_W_DR(L,N)-RESP_W_DR_OUT(L,N))
     &                             /GAM_TRIF   
            ENDDO                                                       
          ENDDO                                                         

        ENDIF
                                                                        
C-----------------------------------------------------------------------
C Reset the accumulated leaf turnover rates to zero.
C-----------------------------------------------------------------------
*D VEG2A.372,VEG2A.382  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX  
     &,           FRAC,HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)           
*/-----
*DECLARE VEG_CTL1
*/-----
*I VEG_CTL1.93    
     &,I                   ! Loop counter for STASHWORK
*I VEG_CTL1.124   
! Initialise STASH workspace to zero to prevent STASH failure
      DO I=1,INT3
         STASHWORK(I)=0.0
       ENDDO

*D ABX3F405.8
     &            LAND_PTS,LAND_LIST,NTILES,P_ROWS,ROW_LENGTH,          
*D VEG_CTL1.130
     &            SECS_PER_STEPim(atmos_im),D1(JDISTURB),
     &            D1(JSAT_SOIL_COND), 
*D VEG_CTL1.134,VEG_CTL1.135  
     &            D1(JCANHT_PFT),D1(JCATCH_TYP),D1(JINFIL_TYP),     
     &            D1(JZ0_TYP),                                  
*D ABX1F405.1302
       CALL STASH(a_sm,a_im,19,STASHWORK,
*/-----
*DECLARE VEG_IC1A
*/-----
*D ABX3F405.13
     &,                 LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS
     &,                 ROW_LENGTH     
*D VEG_IC1A.28
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                 
*D VEG_IC1A.32
     &,                 CATCH_T,INFIL_T,Z0_T                 
*I VEG_IC1A.67    
     &,NTILES                ! IN Number of land-surface tiles.
*D VEG_IC1A.88,VEG_IC1A.89   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).      
*I VEG_IC1A.91    
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).        
*D VEG_IC1A.108,VEG_IC1A.111  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
!                                   !     (kg/m2).
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG_IC1A.112,VEG_IC1A.115  
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D ABX3F405.25
     &,        LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH       
*D VEG_IC1A.127
     &,        ATIMESTEP,SATCON                                       
*D ABX1F405.1334
     &,        CATCH_T,INFIL_T,Z0_T                          
*/-----
*DECLARE VEG_IC2A
*/-----
*D ABX3F405.58
     &,                 LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS
     &,                 ROW_LENGTH          
*D VEG_IC2A.29
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                 
*D VEG_IC2A.33
     &,                 CATCH_T,INFIL_T,Z0_T                 
*I VEG_IC2A.68    
     &,NTILES                       ! IN Number of land-surface tiles.  
*D VEG_IC2A.89,VEG_IC2A.90   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).      
*I VEG_IC2A.92    
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).       
*D VEG_IC2A.111,VEG_IC2A.112  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
C                                   !     (kg/m2). 
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG_IC2A.113,VEG_IC2A.116  
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D ABX3F405.70
     &,        LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH       
*D VEG_IC2A.131
     &,        ATIMESTEP,FRAC_DISTURB,SATCON                          
*D VEG_IC2A.135
     &,        CATCH_T,INFIL_T,Z0_T                          
*/-----
*DECLARE VERSION
*/-----
*D GSS1F400.1151
      PARAMETER(NPSLEVP=500)      ! pseudo levels list

*/-----
*DECLARE VSHRZ7A
*/-----
*D VSHRZ7A.33
     & N_ROWS,FIRST_ROW,ROW_LENGTH,FLANDG, 
*D VSHRZ7A.36
     & VSHR,VSHR_LAND,VSHR_SSI,Z1         
*D ABX1F405.836
     & FLANDG(P_FIELD)             ! IN Land fraction on all tiles
     &,AKH(2)                      ! IN Hybrid 'A' for layer 1.         
*I VSHRZ7A.72    
     &,VSHR_LAND(P_FIELD)          ! OUT Magnitude of surface-to-lowest 
!                                  !     atm level wind shear (m per s).
     &,VSHR_SSI(P_FIELD)           ! OUT Magnitude of surface-to-lowest 
!                                  !     atm level wind shear (m per s).
*I VSHRZ7A.149   
      IF(FLANDG(I).LT.1.0)THEN                                      
*D VSHRZ7A.153
        VSHR_SSI(I) = SQRT(VSHR2)                                   
      ELSE                                                            
        VSHR_SSI(I) = 0.0                                           
      ENDIF                                                           
C                                                                       
      IF(FLANDG(I).GT.0.0)THEN                                      
        VSHR2 = MAX (1.0E-6 , U_1_P(I)*U_1_P(I)                     
     &    + V_1_P(I)*V_1_P(I))                                      
        VSHR_LAND(I) = SQRT(VSHR2)                                    
      ELSE                                                            
        VSHR_LAND(I) = 0.0                                          
      ENDIF                                                           
C                                                                       
      VSHR(I)= FLANDG(I)*VSHR_LAND(I)                           
     &  + (1.0 - FLANDG(I))*VSHR_SSI(I)                           
                                                                 
