* */ Work by Jeff Ridley, Taro Hosoe, Ian Rutt, Oliver Browne, 
*/ Jonathan Gregory, Rupert Gladstone, Robin S Smith, William Roberts,
*/ and Steve George
*/
*/ Includes some code which is Met Office Crown Copyright, so please don't
*/ use it without permission from jonathan.gregory@metoffice.gov.uk.
*/
*/ glsize & lasize are UM variables, for global/local data sizes, they *
*/ relate to the UM grid, not the ISM grid.                            *
*/ JMG 7.12.08
*/ Delete *DEFs and comment out updating of runoff.
*/ RMG 27.02.2009 onwards
*/ see svn change log for details.
*/ http://source.ggy.bris.ac.uk/websvn/listing.php?repname=glimmer-UM&path=%2F&sc=0
*/ 
*/ RSS 08.01.2013 - substantial rewrite to pass Essery-snowscheme SMB
*/ on subgrid levels rather than PDD fields
*/
*/-----
*DECLARE INITIAL1
*/-----  
*B GDG0F401.825
C icesheet coupling needs to be able to force rereading of iceberg ancil
C but we don't want that here
     &                  .FALSE.,0.,
*B GDG0F401.830
C icesheet coupling needs to be able to force rereading of iceberg ancil
C but we don't want that here
     &                  .FALSE.,0.,

*/-----
*DECLARE UPANCIL1
*/-----  
*B GDG0F401.1485
     &                    iceberg_force,iceberg_scale,
*I GDR3F305.318
      LOGICAL ICEBERG_FORCE
      REAL ICEBERG_SCALE
*B GDG0F401.1521
     &                ICEBERG_FORCE,ICEBERG_SCALE,
*/-----
*DECLARE RPANCO1A
*/-----
*B GSS1F304.607
     & ICEBERG_FORCE,ICEBERG_SCALE,
*I RPANCO1A.216
      LOGICAL ICEBERG_FORCE
      REAL ICEBERG_SCALE
*B RPANCO1A.254     
     &                     .OR.(ICEBERG_FORCE .AND. FIELD.eq.24)
*B GMB1F304.166
        IF (FIELD.ne.24) THEN
*I RPANCO1A.790
        ELSE
          write(6,*)"REPLANCO ICBG",ICEBERG_FORCE, ICEBERG_SCALE
          IF (ICEBERG_FORCE) THEN
            DO I=1,LEN_FLD
               D1(POS_STRT+I-1)=ANCIL_DATA(I)
               IF (abs(ANCIL_DATA(I)-RMDI).gt.1)
     &           D1(POS_STRT+I-1)=ANCIL_DATA(I)*ICEBERG_SCALE
            END DO
          ELSE
            write(6,*)"skipping nh icecalv read until forced"
          ENDIF
        ENDIF
*/-----
*DECLARE SFCADD
*/-----
*D ORH1F305.402
     &,                fluxcorh, fluxcorw,fluxicbg
*D OJL1F405.1
     +  ,salref,qfusion)
*I ORH1F305.410
     &,   fluxicbg (IMT_FLX) ! IN  Salt flux correction (SI)
*D OJL1F405.2
        real salref,qfusion
*D  OJL1F405.16
          fluxtoth(I)=fluxcorh(I)+((fluxicbg(i)+fluxcorw(i))*qfusion)
*D OJL1F405.9
          TA(i,1,2)=TA(i,1,2)+con_salt*(fluxcorw(I)+fluxicbg(i))*salref
*D OJL1F405.11
          TA(i,1,2)=TA(i,1,2)+con_salt*(fluxcorw(I)+fluxicbg(i))
     &                                *(0.035+SREF(I))
*B OJP0F404.957
      IF (L_FLUXCORR) VTCO2_FLUX(i)=VTCO2_FLUX(i)
     &  + con_salt*(fluxicbg(i))*TB(i,1,TCO2_TRACER)/C2DTTS
*B OJP0F404.963 
      IF (L_FLUXCORR) VALK_FLUX(i)=VALK_FLUX(i)
     &  + con_salt*(fluxicbg(i))*TB(i,1,ALK_TRACER)/C2DTTS
*/---------
*DECLARE TRACER
*/----------
*D OJT0F304.54
     &,fluxcorh,fluxcorw,fluxicbg
*I ORH1F305.586
     +,  fluxicbg  (IMT_FLX) ! IN  Water flux correction (kg/m2/s)
*D ORH1F305.1318
     +,              fluxcorh,fluxcorw,fluxicbg
*D OJL1F405.19
     +   ,salref,qfusion)
*D OJT0F304.64
     +,              fluxcorh,fluxcorw,fluxicbg
*D OJL1F405.20
     +   ,salref,qfusion)
*/----------
*DECLARE CAOPTR
*/----------
*I CCN1F405.22
     &,   JA_ICERUNOFF          ! greenland icesheet runoff
*I CCN1F405.23
     &,   JA_ICERUNOFF
*/
*DECLARE SWAPA2O2
*/
*I SWAPA2O2.149
     & ,ICEBERGS(G_P_FIELD)         ! greenland iceberg melt rate
*I SWAPA2O2.163
     & ,ICERUNOFF(G_P_FIELD)        ! runoff from greenland icesheet
*I SWAPA2O2.305
      
*I SWAPA2O2.364
      CALL GATHER_FIELD(D1(JA_ICERUNOFF),ICERUNOFF,
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)   
      IF(info.NE.0) THEN      ! Check return code  
         CMESSAGE='SWAPA2O : ERROR in gather of ICERUNOFF' 
         ICODE=31
         GO TO 999
      ENDIF

*I SWAPA2O2.406
              !... for ocean points 
!JMG          IF(.not.atmos_landmask(I)) LSRAIN(I) = LSRAIN(I)+ICERUNOFF(I)
*/
*/**********************************************************************
*/
*DECLARE U_MODEL1
*B U_MODEL1.70
      logical iceberg_force
      real iceberg_scale

      interface
        subroutine glim_initialise(mype,
     +    A_LEN_REALHD,
     +    A_REALHD,
     +    nx,ny,
     +    iiday,iimonth,iiyear 
     +     )


        integer,intent(in)                  :: A_LEN_REALHD,mype
        integer(kind=8),optional,intent(in) :: nx,ny
                        ! nx ny  Size of normal global grid
        REAL(kind=8),optional, INTENT(IN),dimension(A_LEN_REALHD):: 
     +      A_REALHD						
        ! A_REALHD contains grid resolution info
        integer(kind=8),optional            :: iiday,iimonth,iiyear

        end subroutine glim_initialise
      end interface

      interface
        subroutine glim_finalise(mype,
     +    A_LEN_REALHD,
     +    iiday,iimonth,iiyear 
     +     )


        integer,intent(in)                  :: A_LEN_REALHD,mype
        integer(kind=8),optional            :: iiday,iimonth,iiyear

      end subroutine glim_finalise
      end interface
*I U_MODEL1.77
      iceberg_force=.FALSE.
      iceberg_scale=1e-3
*I NF171193.10
     &   iceberg_force,iceberg_scale,
*B GDG0F401.1478
     &                   iceberg_force,iceberg_scale,
*I U_MODEL1.86
      ! -----calls to glimmer initialisation---

      call glim_initialise(mype,
     +                    A_LEN_REALHD,
     +                    A_REALHD=A_SPDUM(A_IXDUM(6)),
     +                    nx=glsize(1),
     +                    ny=glsize(2),
     +                    iiday=i_day,
     +                    iimonth=i_month,
     +                    iiyear=i_year
     +                   )
      write(6,*)'GLIM_INITIALISE: u_model1 glsize(1)=',glsize(1),
     +            'glsize(2)=',glsize(2)
      write(6,*)'GLIM_INITIALISE: u_model1 A_SPDUM(A_IXDUM(6))=',
     +            A_SPDUM(A_IXDUM(6))

! we need another call to close and de-allocate
*I  U_MODEL1.276 
      write(6,*)'Calling GLIM_FINALISE'
      call glim_finalise( mype,
     +                    A_LEN_REALHD,
     +                    iiday=i_day,
     +                    iimonth=i_month,
     +                    iiyear=i_year
     +                    )

*/--------------
*DECLARE ATMSTEP1
*/--------------
*I NF171193.13
     &   iceberg_force,iceberg_scale,
*B ATMSTEP1.75
      LOGICAL ICEBERG_FORCE
      REAL ICEBERG_SCALE

*B ARB0F400.72
     &              iceberg_force,iceberg_scale,

*/--------------
*DECLARE ATMPHY1
*/--------------
*B GSS1F304.1184
     &             iceberg_force,iceberg_scale,
*B ATMPHY1.54

      logical iceberg_force
      real iceberg_scale


*I ATMPHY1.395
      IF(LTIMER) THEN     
        CALL TIMER('GLIM_CTL ',3)
      END IF 

      CALL GLIM_CTL(MY_PROC_ID, 
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS   
*CALL ARGCONA
*CALL ARGPTRA  
     & dzsnow,iceberg_force,iceberg_scale,
     & ICODE)

      IF (L_WRIT_PHY .AND.
     &    (A_STEP.LE.T_WRITD1_END .OR. T_WRITD1_END .EQ. 0)) THEN
     
      IF (A_STEP.EQ.T_WRITD1_START .OR.
     &    WRITD1_TEST.GT.WRITD1_TEST_PREV) THEN

      CALL DUMPCTL(
*CALL ARGSIZE 
*CALL ARGD1  
*CALL ARGDUMA
*CALL ARGDUMO
*CALL ARGDUMW
*CALL ARGCONA
*CALL ARGPTRA
*CALL ARGSTS 
*CALL ARGPPX 
     &          atmos_sm,0,.TRUE.,'af_glimctl',a_step,
     &          ICODE,CMESSAGE)

      END IF

      END IF

      IF(LTIMER) THEN     
        CALL TIMER('GLIM_CTL ',4)
      END IF 

      IF(ICODE.GT.0) RETURN 
*DECLARE OCNFRST1
*D OJT0F304.7
     &,D1(joc_anom_heat),D1(joc_anom_salt),D1(jousr_anc1)
*DECLARE OCN_CTL
*D OJT0F304.10
     &,fluxcorh,fluxcorw,fluxicbg
*I ORH1F305.5544
     &,fluxicbg(IMT_FLX,JMT_flx)
*D ORH1F305.5595
     &,fluxcorh,fluxcorw,fluxicbg
*DECLARE ROWCALC
*D ORH5F401.33
     &,CARYSALT, anomiceh, fluxcorh, fluxcorw,fluxicbg
*I ORH1F305.2009
     &,fluxicbg(IMT_FLX,JMT_flx)
*D ORH1F305.2324
     &,fluxcorh(1,J_FLUX),fluxcorw(1,J_FLUX),fluxicbg(1,J_FLUX)
*DECLARE ROW_CTL
*D OJT0F304.30
     &,fluxcorh,fluxcorw,fluxicbg
*I ORH1F305.510
     & ,fluxicbg(IMT_FLX,JMT_FLX)
*D ORH3F403.253
     &,ICY,FLXTOICE,CARYSALT,anomiceh,fluxcorh,fluxcorw,fluxicbg
*DECLARE BLOKCALC
*D ORH3F403.244
     &,ICY,FLXTOICE,CARYSALT,anomiceh,fluxcorh,fluxcorw,fluxicbg
*I ORH1F305.3577
     &,fluxicbg(*)
*D ORH5F401.7
     &,CARYSALT, anomiceh, fluxcorh, fluxcorw,fluxicbg
*DECLARE BLOKCNTL
*D ORH5F401.16
     &,CARYSALT, anomiceh, fluxcorh, fluxcorw,fluxicbg
*I ORH1F305.3427
     &,fluxicbg(*)
*D ORH5F401.26
     &,CARYSALT, anomiceh, fluxcorh, fluxcorw,fluxicbg

*/
*DECK GLIM_CTL
*/
      SUBROUTINE GLIM_CTL(MY_PROC_ID,
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGCONA
*CALL ARGPTRA
     & dzsnow,iceberg_force,iceberg_scale,
     & ICODE)

      IMPLICIT NONE

!----------------------------------------------------------------------!
! Controls the call to the icesheet model. The icesheet model, Glimmer !
! is run on processor one at the end of each year. UM fields needed to !
! drive it must be gathered & those updated scattered back to the MMP  !
! afterwards.                                                          !
!----------------------------------------------------------------------!

!RMG need to work out which of these/ updates can be removed post restructuring 
!into new gather and scatter subroutines.

*CALL CMAXSIZE     ! contains array size limits for other blocks
*CALL CSUBMODL     ! contains size information for CTIME and a_im
*CALL TYPSIZE      ! contains info for dynamical sizing of field arrays
*CALL TYPD1        ! contains information on D1 addressing array
*CALL TYPSTS       ! STASH related variables
*CALL TYPDUMA      ! DUMP headers for atmospheric model
*CALL TYPPTRA      ! D1 pointers for atmospheric model
*CALL PARVARS      ! contains size information about local & global 
                   !  fields & nproc
*CALL CTIME        ! contains time information
*CALL CAOPTR       ! common block contains pointer JA_FASTRUNOFF
*CALL C_MDI        ! common block containing missing data indicators
*CALL GCCOM        ! allows cross processor talking

*CALL NSTYPES      ! surface/elevation type sizes
*CALL C_A
*CALL TYPCONA

!-- parameters -------!
      INTEGER, PARAMETER    :: gather_pe     = 0
      INTEGER, PARAMETER    :: ICE_COUPLING_MONTH     = 12
!+seg Define ICE_RESET_SNOWMASS
      REAL, PARAMETER       :: GLIMMER_ICE_DENSITY    = 910.
      REAL, PARAMETER       :: GLIMMSET_RS_SDEPTH     = 200.  ! Increased from 100 4/19
      REAL, PARAMETER       :: ICE_RESET_SNOWMASS  = 
     &      GLIMMER_ICE_DENSITY * GLIMMSET_RS_SDEPTH 
!-seg
!+ From Robin 31-5-17 REFAC1.51
      REAL, PARAMETER       :: glimmer_ice_grain=2500.
      REAL, PARAMETER       :: glimmer_ice_temp=273.15
!- REFAC1.51
!-- local variables --!
      CHARACTER (LEN=80) :: CMESSAGE        !  Error return message  

      INTEGER    :: info,
     &              htime,                     !
     &              my_proc_id,                !
     &              ICODE,                     ! Error return code
     &              G_P_FIELD,                 ! global size of p_
     &              G_LAND_FIELD               ! global size of la

      logical    :: ocn_interface   ! is ocean coupling used?
      logical    :: l_ice_couple_step
!+seg From Robin 31-5-17 REFAC1.52
! But, with 10 levels rather than RS 25
      REAL       :: forcing_time_step
      REAL       :: elevation(nelev)
      data elevation / 100., 300., 550., 850., 1150.,
     &                 1450.,1800.,2250.,2750.,3600./
!-seg REFAC1.52

      real       :: calving     !RMG always zero
      real       :: dummy

!SINGLE PE VERSIONS OF THINGS TO GO INTO GLIMMER
      real       :: ice_snowmass_global(glsize(1),glsize(2),nelev)
      real       :: ice_smb_global(glsize(1),glsize(2),nelev)
      real       :: ice_areas_global(glsize(1),glsize(2),nelev)
      real       :: nonice_snowmass_global(glsize(1),glsize(2),nelev)
      real       :: nonice_snowdepth_global(glsize(1),glsize(2),nelev)
      real       :: icestemp_global(glsize(1),glsize(2),nelev)

!SINGLE PE VERSIONS OF THINGS THAT COME OUT OF GLIMMER
      real       :: landfrac_global(glsize(1),glsize(2),nelev)
      real       :: icefrac_global(glsize(1),glsize(2),nelev)
      real       :: icehflux_global(glsize(1),glsize(2),nelev)
      real       :: icehflux_local(lasize(1),lasize(2),nelev)
      real       :: icerunoff_global(glsize(1),glsize(2),nelev)
      real       :: liqrunoff_global(glsize(1),glsize(2),nelev)
!+seg 04/19
      real       :: orogMN_global(glsize(1),glsize(2))
!-seg
      real       :: orogSD_global(glsize(1),glsize(2))
      real       :: orogHO_global(glsize(1),glsize(2))
      real       :: orogSIL_global(glsize(1),glsize(2))
      real       :: orogXX_global(glsize(1),glsize(2))
      real       :: orogXY_global(glsize(1),glsize(2))
      real       :: orogYY_global(glsize(1),glsize(2))
      real       :: waterout_global(glsize(1),glsize(2))
      real       :: land_area_global(glsize(1),glsize(2))

      real       :: dsice,dzsnow(10),delta_vol,ice_smb_total
      real       :: ice_areas_total
!+seg
      real       :: dfrac
      integer    :: dindex
!-seg
      INTEGER    :: n,nl,l,i,j,ns

      real       :: iceberg_scale
      logical    :: iceberg_force

      !----------------------------------------------------------------!
      ! Set up timing characteristics                                  !
      !----------------------------------------------------------------!

      htime            = I_YEAR*360*24 + I_DAY_NUMBER*24 + I_HOUR
                                                       !abstime in hours

      ! initialise vars relating to error reporting 
      cmessage = " "
      icode    = 0
      iceberg_force=.false.

!+seg
      !!HARD-CODED TO ANNUAL!! MATCH WITH GLIM_HST
      forcing_time_step=60.*60.*24.*360.
!-seg
      !---------------------------------------------------------------!
      ! Is this the right time to do something
      L_ICE_COUPLE_STEP=.FALSE.
      if (     I_MONTH.eq.ICE_COUPLING_MONTH 
     &    .AND.PREVIOUS_TIME(2) /= I_MONTH
     &   ) then

        L_ICE_COUPLE_STEP=.TRUE.
        iceberg_force=.TRUE.

      endif

      IF (.NOT.L_ICE_COUPLE_STEP .AND. PREVIOUS_TIME(3) /= I_DAY) 
     & write(6,*)"not coupling - want 1st tstep of month"
     & ,ICE_COUPLING_MONTH

!+seg
! Mass reset (above) to 910 * 100 kg
!        ice_reset_snowmass=0.
!        do ns=1,10
!          ice_reset_snowmass=ice_reset_snowmass
!     &                     +(2*dzsnow(ns))*GLIMMER_ICE_DENSITY
!        end do
!
!-seg
      IF (L_ICE_COUPLE_STEP) THEN

      write(6,*)"Coupling to ice sheet"
      !---------------------------------------------------------------!
        Call gather_ice_fields(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     & icefrac_global,landfrac_global,
     & ice_snowmass_global,nonice_snowmass_global,
     & icestemp_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
     & land_area_global,
!+seg 04/19 New field
     & elevation,
!-seg
     & cos_p_latitude,ice_reset_snowmass
     & )

!      !MAKE ICEMASS INTO AN SMB TERM, RESET UM ICEMASS
!      !---------------------------------------------------------------!
! ASSUME NO SIGNIFICANT CHANGES TO SNOWGRAIN, SNOWTEMP FROM RESETTING OF MASS IN
! BOTTOM LAYER. GBM SNOWMASS WANTED, THRESHOLD CHECK IN LWRAD BUT
! NOT GOING TO BE AN ISSUE IN OUR CASE
! NEED TO CHANGE SNOWMASS(NSMAX),SNOWDEPTH(TILE),SNOWMASS(TILE),SNOW_DENSITY(NSMAX)

        write(6,*)"back from gather ice"
        if (mype == gather_pe) then

        ice_smb_total=0.
        ice_areas_total=0.

            do j=1,glsize(2)
            do i=1,glsize(1)
          do nl=1,nelev
               !SMB WANTED IN kg/m2/s, go from mass/yr->mass/s
        ice_smb_global(i,j,nl)=(ice_snowmass_global(i,j,nl)-
     &                                 ICE_RESET_SNOWMASS)/
     &                                 forcing_time_step  ! seg Defined above

        ice_areas_global(i,j,nl)= land_area_global(i,j)
     &                * (icefrac_global(i,j,nl)*landfrac_global(i,j,nl))


        nonice_snowdepth_global(i,j,nl)=
     &                nonice_snowmass_global(i,j,nl)/GLIMMER_ICE_DENSITY
     !!
!           nonice_snowdepth_global(i,j,nl)=1.
!           ice_smb_global(i,j,nl)=GLIMMER_ICE_DENSITY/(60.*60.*24.*360.)
     !!
        if (j.lt.17) then
!ie just NH, harwired for FAMOUS, ignore Antarctica
        ice_areas_total=ice_areas_total+ice_areas_global(i,j,nl)
        ice_smb_total=ice_smb_total
     &               +(ice_snowmass_global(i,j,nl)-ICE_RESET_SNOWMASS)
!     &               +GLIMMER_ICE_DENSITY
     &               *ice_areas_global(i,j,nl) !kg.yr

        endif

            end do
            end do
          end do

        endif

        DO N=NTILES-NELEV+1,NTILES
         nl=mod(n-1,nelev)+1
         do l=1,land_field

           ns=int(D1(JNSNOW(n)+l-1))
           dsice=D1(JSNODEP_TYP+((n-1)*land_field)+l-1)-
     &                     ICE_RESET_SNOWMASS

           !ON SNOW LAYERS
!+seg
           ! Prog 216, size 225
           dfrac=D1(JFRAC_TYP+((ntype-nelev+nl-1)*land_field)+l-1)
           if (dfrac .gt. 0) then !+dfrac
            if (dfrac .gt. 0 .and. ns .lt. nsmax) then !+dfrac_err
             WRITE(6,*) 'GLIM_CTL: Bad snow layers with ice fraction'
             WRITE(6,*) 'DFRAC = ',dfrac,', ns = ',ns
            endif                                      !-dfrac_err
            dindex = (ntiles*ns)-(2*nelev)+n
            ! Only do next bit if valid index calculated
            if (dindex > 0) then  !+dindex
!-seg
             if (dsice .gt. 
     &        D1(JSICE((ntiles*ns)-(2*nelev)+n)+l-1)
     &        ) then
             write(6,*)"glim_ctl: WARNING, not enough ice mass in ",
     & "lowest snow layer to subtract the SMB we're passing to Glimmer",
     & l,n,ns,dsice,D1(JSICE((ntiles*ns)-(2*nelev)+n)+l-1)
     &,D1(JSNODEP_TYP+((n-1)*land_field)+l-1)
             endif

!           D1(JSICE(n*nsmax)+l-1)=D1(JSICE(n*nsmax)+l-1)-dsice
           !!or does layer addressing run layer major (think it does)
           D1(JSICE((ntiles*ns)-(2*nelev)+n)+l-1)=
     &         D1(JSICE((ntiles*ns)-(2*nelev)+n)+l-1)-dsice
!!           D1(JDS(n*nsmax)+l-1)=D1(JDS(n*nsmax)+l-1)+(dsice*GLIMMER_ICE_DENSITY)
!!           D1(JRHO_SNOW(n*nsmax)+l-1)=(D1(JSICE(n*nsmax)+l-1)+D1(JSLIQ(n*nsmax)+l-1))/D1(JDS(n*nsmax)+l-1)

           !ON SNOWPACK AVG (TILES)
             D1(JSNODEP_TYP+((n-1)*land_field)+l-1)=ICE_RESET_SNOWMASS
!+seg Sign change for from + to -(dsice...             
             D1(JSNOWDEPTH(n)+l-1)=D1(JSNOWDEPTH(n)+l-1)-
     &                           (dsice/GLIMMER_ICE_DENSITY)
!           D1(JSNOWDEPTH(n)+l-1)=ICE_RESET_SNOWMASS/GLIMMER_ICE_DENSITY
             D1(JRHO_SNOW_GRND(n)+l-1)=
     &                      D1(JSNODEP_TYP+((n-1)*land_field)+l-1)/
     &                      D1(JSNOWDEPTH(n)+l-1)
!+seg
            else 
              write(6,*)"glim_ctl: WARNING, invalid index",dindex
            end if !-dindex        
           end if  !-dfrac
!-seg
        end  do
        END  DO

        if (mype == gather_pe) then
!
!          !! The making of new ice points from nonice snowdepth
!          !! on the Glimmer grid has to happen AFTER 3D interp
!          !! of the snowdepth field, so will be in GLINT

          Call glim_intctl(
     &               glsize(1),glsize(2),nelev,
! INTO GLIMMER
     &               htime,
     &               ice_smb_global,
     &               ice_areas_global,
     &               nonice_snowdepth_global,
     &               icestemp_global,
! INOUT OF GLIMMER
     &               landfrac_global,
     &               icefrac_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
! OUT OF GLIMMER
     &               waterout_global,
!     &               icerunoff_global,
!     &               liqrunoff_global,
     &               icehflux_global
     &               ,delta_vol
     &               )

         write(6,*)"back from gl_ic", delta_vol,ice_areas_total

         endif ! mype==gather_pe

        call gc_gsync(nproc,info)

!
!
!      !---------------------------------------------------------------!
          Call scatter_ice_fields(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     & landfrac_global,icefrac_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
     & waterout_global,
     & nonice_snowdepth_global, nonice_snowmass_global,
     & icehflux_global,icehflux_local
!      &, icerunoff_global, liqrunoff_global
     &, delta_vol,ice_smb_total,land_area_global
     &, iceberg_scale,dzsnow,
!+seg 04/19 New fields
     &   glimmer_ice_density,glimmer_ice_grain,glimmer_ice_temp,
     &   forcing_time_step
!-seg
     & )

           Call glim_setsurf(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     &   icehflux_local,
     &   ICODE)

      ENDIF !L_ICE_COUPLING_STEP

!!----------------------------------------------------------------------!
! Error trap.                                                          !
!----------------------------------------------------------------------!
!RMG remove this?

 999  CONTINUE                                                          

      if (ICODE.NE.0) then
         write(6,*)'GLIM_CTL ',CMESSAGE,ICODE
      endif

      RETURN
                                                             
      END SUBROUTINE glim_ctl

*/**********************************************************************

*DECK GLIM_GATHSCAT

      SUBROUTINE gather_ice_fields(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     & icefrac_global,landfrac_global,
     & ice_snowmass_global,nonice_snowmass_global,
     & icestemp_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
     & land_area_global,
!+seg 04/19 New field
     & elevation, 
!-seg
     & cos_p_latitude,ice_reset_snowmass
     & )

        IMPLICIT NONE

*CALL CMAXSIZE     ! contains array size limits for other blocks
*CALL CSUBMODL     ! contains size information for CTIME and a_im
*CALL TYPSIZE      ! contains info for dynamical sizing of field arrays
*CALL TYPD1        ! contains information on D1 addressing array
*CALL TYPSTS       ! STASH related variables
*CALL TYPDUMA      ! DUMP headers for atmospheric model
*CALL TYPPTRA      ! D1 pointers for atmospheric model
*CALL PARVARS      ! contains size information about local & global 
                   !  fields & nproc
*CALL CTIME        ! contains time information
*CALL CAOPTR       ! common block contains pointer JA_FASTRUNOFF
*CALL C_MDI        ! common block containing missing data indicators
*CALL GCCOM        ! allows cross processor talking
*CALL NSTYPES      ! surface/elevation type sizes

*CALL AMAXSIZE
*CALL ATM_LSM
*CALL C_A

!SINGLE PE VERSIONS OF THINGS TO GO INTO GLIMMER
      real       :: ice_snowmass_global(glsize(1),glsize(2),nelev)
      real       :: nonice_snowmass_global(glsize(1),glsize(2),nelev)
      real       :: icestemp_global(glsize(1),glsize(2),nelev)
      real       :: land_area_global(glsize(1),glsize(2))
      real       :: fland_global(glsize(1),glsize(2))
      real       :: cos_p_latitude(p_field)
      real       :: cos_p_latitude_global(glsize(1)*glsize(2))

!SINGLE PE VERSIONS OF THINGS THAT COME OUT OF GLIMMER
      real       :: landfrac_global(glsize(1),glsize(2),nelev)
      real       :: icefrac_global(glsize(1),glsize(2),nelev)
!+seg 04/19 New field
      real       :: orogMN_global(glsize(1),glsize(2))
!-seg
      real       :: orogSD_global(glsize(1),glsize(2))
      real       :: orogHO_global(glsize(1),glsize(2))
      real       :: orogSIL_global(glsize(1),glsize(2))
      real       :: orogXX_global(glsize(1),glsize(2))
      real       :: orogXY_global(glsize(1),glsize(2))
      real       :: orogYY_global(glsize(1),glsize(2))

      real       :: icestemp(glsize(1),glsize(2))


      INTEGER, PARAMETER :: gather_pe     = 0
      CHARACTER (LEN=80) :: CMESSAGE        !  Error return message

      INTEGER    :: snow_btemp,STASHMACRO_TAG
      INTEGER    :: idummy,info
      INTEGER    :: nl,i,j,n

      real :: snowmass_all(glsize(1),glsize(2),ntiles)
      real :: landfrac_all(glsize(1),glsize(2),ntype)

      real :: land_loc(lasize(1)*lasize(2))
      real :: fracsum

      REAL delevation,icegrad,soilgrad,pi,ice_reset_snowmass
!+seg 04/19 Elevations now defined elsewhere, now argument in call to subroutine
      REAL elevation(nelev)
!-seg
      data soilgrad /-.001/
      data icegrad /-.001/

      write(6,*) "starting gather ice fields"

!+seg 04/19 Straight D1 array as it is not on land-grid
        CALL GATHER_FIELD(D1(JOROG),
     &  orogMN_global,
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
!-seg

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_SD),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogSD_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_SIL),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogSIL_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_HO2),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogHO_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_XX),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogXX_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_XY),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogXY_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_yy),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  orogYY_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

c    !      write(6,*) "do land_area_global"
c! gather and do calc on 1 pe rather than local calc then gather to 
c! save having to explicitly make a land_index to relate fland to cos_p 

        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JFRAC_LAND),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &  fland_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL GATHER_FIELD(cos_p_latitude, ! LAND AREA FRACTIONS
     &  cos_p_latitude_global,
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)

        if (mype.eq.gather_pe) then
        do j=1,glsize(2)
        do i=1,glsize(1)
          if (fland_global(i,j).lt.0) fland_global(i,j)=0.

          pi=3.1415927
c!area=Re^2 * dx * 2 * sin(dy/2) *cos(y) dx,dy in lat
c!http://ferret.pmel.noaa.gov/Ferret/faq/averages-integrals-on-the-sphere
          land_area_global(i,j)=
     &                (A**2)*(2*pi)/glsize(1)*2*sin(pi/2./(glsize(2)-1))
     &              * COS_P_LATITUDE_global(((j-1)*row_length)+i)
     &              * fland_global(i,j)
        end do
        end do
        end if

      STASHMACRO_TAG=88

        ! Call up snowpack bottom buffer temp
      CALL FINDPTR(ATMOS_IM, 8,432,IMDI,IMDI,IMDI,IMDI
     &     ,IMDI,IMDI,IMDI,IMDI,IMDI,IMDI
     &     ,IMDI,IMDI,IMDI,STASHMACRO_TAG,IMDI,snow_btemp,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,
     &     INFO,CMESSAGE)

!       write(6,*)"found 8432 with pointer",snow_btemp

!      ???CALL FROM_LAND_POINTS NEEDED, IS IN STASH??
!
      CALL GATHER_FIELD(D1(snow_btemp), ! mean ice temp
     &     icestemp,
     &     lasize(1),lasize(2),glsize(1),glsize(2),
     &     gather_pe,GC_ALL_PROC_GROUP,info)
      IF (MYPE.eq.GATHER_PE) THEN
        do j=1,glsize(2)
        do i=1,glsize(1)
          if (abs(icestemp(i,j)-RMDI).lt.1) then
            icestemp(i,j)=0
          else
            icestemp(i,j)=icestemp(i,j)-273.15
          endif
!          write(300,*)icestemp(i,j)
        end do
        end do
!        write(300,*)"-"
      ENDIF
      do nl=1,nelev
        do j=1,glsize(2)




        do i=1,glsize(1)
!+seg
!         delevation=elevation(n)-D1(JOROG+(j*glsize(1))+i-1)
          delevation=elevation(nl)-D1(JOROG+(j*glsize(1))+i-1)
!-seg
          icestemp_global(i,j,nl)=icestemp(i,j)+(delevation*icegrad)
        end do
        end do

      end do
!      IF (MYPE.eq.GATHER_PE) THEN
!        do j=1,glsize(2)
!        do i=1,glsize(1)
!          write(400,*)icestemp(i,j)
!        end do
!        end do
!        write(400,*)"-"
!      ENDIF

      do nl=1,ntype
        CALL FROM_LAND_POINTS(LAND_LOC,
     &   D1(JFRAC_TYP+(nl-1)*LAND_FIELD),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC, ! LAND AREA FRACTIONS
     &  landfrac_all(1,1,nl),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
      end do

!+seg New scheme 12Sep17
      !!!!!!!!!!!!!!!!FOR SNOW ACCEL, CALL ROUTINE HERE. WANT TO ADJUST 
      !!!!!!!!!!!!!!!!NONICE SNOW BEFORE IT'S GATHERED, REALLY
      !!!!!!!!!!!!!!!!NEEDS TO PASS BITS OF D1 INTO A SEP ROUTINE WHICH
      !!!!!!!!!!!!!!!!CAN THEN USE THE Glimmer MODULES TO INTERROGATE
      !!!!!!!!!!!!!!!!THE GLIMMER NAMELIST

      CALL nonice_snow_accel(mype
     &                      ,gather_pe
     &                      ,nproc
     &                      ,land_field
     &                      ,ntiles
     &                      ,nelev
     &                      ,nsmax
     &                      ,D1(JRGRAIN_TYP)
     &                      ,D1(JSNODEP_TYP)
     &                      ,D1(JSNODEP_TYP_LC)
     &                      ,D1(JSNOWDEPTH(1))
     &                      ,D1(JRHO_SNOW_GRND(1))
     &                      ,D1(JNSNOW(1))
     &                      ,D1(JDS(1))
     &                      ,D1(JSICE(1))
     &                      ,D1(JSLIQ(1))
     &                      ,D1(JRHO_SNOW(1))
     &                      ,D1(JRGRAINL(1))
     &                      )
!-seg 12Sep17

      do nl=1,ntiles
        CALL FROM_LAND_POINTS(LAND_LOC,
     &    D1(JSNODEP_TYP+(nl-1)*land_field),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC, ! snowmass on tiles
     &    snowmass_all(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
!        CALL FROM_LAND_POINTS(LAND_LOC,
!     &    D1(JRHO_SNOW_GRND(nl)),
!     &    atmos_landmask_local,
!     &    lasize(1)*lasize(2),idummy)
!        CALL GATHER_FIELD(LAND_LOC, ! snowdensity on tiles
!     &    snowdens_all(1,1,nl),
!     &    lasize(1),lasize(2),glsize(1),glsize(2),
!     &    gather_pe,GC_ALL_PROC_GROUP,info)
      enddo

      IF (MYPE.eq.GATHER_PE) THEN
!      do j=1,glsize(2)
!      do i=1,glsize(1)
!      fracsum=0.
!      do nl=1,ntype
!       if (abs(landfrac_all(i,j,nl)).le.1) then
!         fracsum=fracsum+landfrac_all(i,j,nl)
!       endif
!      end do
!      write(6,*)i,j,fracsum
!      end do
!      end do

      do j=1,glsize(2)
      do i=1,glsize(1)

      do nl=1,nelev
        landfrac_global(i,j,nl)=0.
        nonice_snowmass_global(i,j,nl)=0.








      end do

      do n=1,ntype
        nl=mod(n-1,nelev)+1
        landfrac_global(i,j,nl)=landfrac_global(i,j,nl)
     &                         +landfrac_all(i,j,n)
      end do

! GLIMMER expects fraction of /elevation/ that's ice, not fraction of column.
! I think.
      do nl=1,nelev
        if (landfrac_global(i,j,nl).gt.0) then
         icefrac_global(i,j,nl)=landfrac_all(i,j,ntype-nelev+nl)/
     &                         landfrac_global(i,j,nl)
         ice_snowmass_global(i,j,nl)=snowmass_all(i,j,ntiles-nelev+nl)
        else
          icefrac_global(i,j,nl)=0.
          ice_snowmass_global(i,j,nl)=ice_reset_snowmass
        endif
      end do

! Aggregate non-ice tile snow depths on elevations, normalising the raw column
! fractions by the fraction of non-ice-land at that elevation
!+seg 04/19 Correction to indexing
      do nl=1,(ntiles-nelev)
!-seg
        if (icefrac_global(i,j,nl).lt.1 .AND.
     &      landfrac_global(i,j,nl).gt.0
     &     ) then
!+seg 04/19 Correction to indexing
        nonice_snowmass_global(i,j,nl)=snowmass_all(i,j,nl)
!-seg
        else
          nonice_snowmass_global(i,j,nl)=0.
        endif
      end do

      end do !j
      end do !i
      ENDIF !MYPE

      write(6,*) "finishing gather ice fields"

      return

      END SUBROUTINE !gather_ice_fields

      SUBROUTINE scatter_ice_fields(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     & landfrac_global,icefrac_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
     & waterout_global,
     & nonice_snowdepth_global, nonice_snowmass_global,
     & icehflux_global, icehflux_local 
!     & ,icerunoff_global, liqrunoff_global
     & ,delta_vol,ice_smb_total,land_area_global
     & ,iceberg_scale,dzsnow,
!+seg 04/19 New fields from RS
     &   glimmer_ice_density,glimmer_ice_grain,glimmer_ice_temp,
     &   forcing_time_step
!-seg
     & )

        IMPLICIT NONE

*CALL CMAXSIZE     ! contains array size limits for other blocks
*CALL CSUBMODL     ! contains size information for CTIME and a_im
*CALL TYPSIZE      ! contains info for dynamical sizing of field arrays
*CALL TYPD1        ! contains information on D1 addressing array
*CALL TYPSTS       ! STASH related variables
*CALL TYPDUMA      ! DUMP headers for atmospheric model
*CALL TYPPTRA      ! D1 pointers for atmospheric model
*CALL PARVARS      ! contains size information about local & global 
                   !  fields & nproc
*CALL CTIME        ! contains time information
*CALL CAOPTR       ! common block contains pointer JA_FASTRUNOFF
*CALL C_MDI        ! common block containing missing data indicators
*CALL GCCOM        ! allows cross processor talking
*CALL NSTYPES      ! surface/elevation type sizes

*CALL AMAXSIZE
*CALL ATM_LSM

      real       :: landfrac_global(glsize(1),glsize(2),nelev)
      real       :: icefrac_global(glsize(1),glsize(2),nelev)
      real       :: nonice_snowdepth_global(glsize(1),glsize(2),nelev)
      real       :: nonice_snowmass_global(glsize(1),glsize(2),nelev)
!+seg 15/05/19
      real       :: orogMN_global(glsize(1),glsize(2))
!-seg
      real       :: orogSD_global(glsize(1),glsize(2))
      real       :: orogHO_global(glsize(1),glsize(2))
      real       :: orogSIL_global(glsize(1),glsize(2))
      real       :: orogXX_global(glsize(1),glsize(2))
      real       :: orogXY_global(glsize(1),glsize(2))
      real       :: orogYY_global(glsize(1),glsize(2))

      real       :: waterout_global(glsize(1),glsize(2))
      real       :: icehflux_global(glsize(1),glsize(2),nelev)
      real       :: icehflux_local(lasize(1),lasize(2),nelev)
      real       :: icerunoff_global(glsize(1),glsize(2),nelev)
      real       :: liqrunoff_global(glsize(1),glsize(2),nelev)

      real       :: land_area_global(glsize(1),glsize(2))


      INTEGER, PARAMETER :: gather_pe     = 0
      CHARACTER (LEN=80) :: CMESSAGE        !  Error return message

      INTEGER    :: idummy,info
      INTEGER    :: nl,i,j,n
      real    :: dzsnow(nsmax)

!+seg From Robin 31-5-17 REFAC1.43
      REAL     :: forcing_time_step
      REAL     :: glimmer_ice_density,glimmer_ice_grain,glimmer_ice_temp
!-seG

      real :: land_loc(lasize(1),lasize(2))
      real :: glimfrac, fracsum,
     &        icefrac_local(land_field,nelev),
     &        landfrac_local(land_field,nelev),
     &        nisnowdepth_local(land_field,nelev),
     &        nisnowmass_local(land_field,nelev),
!
     &        icearea_local(nelev),
     &        onilandfrac_local(nelev),
     &        nnilandfrac_local(land_field,nelev),
     &        nnilandfrac_global(glsize(1),glsize(2),nelev)

      real :: dsice,delta_vol,ice_smb_total,delta_nis
     &       ,ice_runoff_total,conserv_global,ice_runoff_A
!+seg 04/19 GLIMMER_ICE_DENSITY,GLIMMER_ICE_GRAIN,GLIMMER_ICE_TEMP 
!     parameters removed
      real            :: iceberg_scale
      real            :: oldmass,newmass

      write(6,*)"in scatter fields",mype

!      write(6,*)"skipping some of scatter fields",mype
!      goto 6969

!+seg 04/19 Scattering new field
        CALL SCATTER_FIELD(D1(JOROG),
     &  orogMN_global,
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(D1(JOROG),lasize(1),lasize(2),offx,offy,1)
!-seg

        CALL SCATTER_FIELD(LAND_LOC,
     &  orogSD_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_SD),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)


        CALL SCATTER_FIELD(LAND_LOC,
     &  orogHO_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_HO2),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(LAND_LOC,
     &  orogSIL_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_SIL),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)
        CALL SCATTER_FIELD(LAND_LOC,
     &  orogXX_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_XX),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(LAND_LOC,
     &  orogXY_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_XY),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(LAND_LOC,
     &  orogYY_global(1,1),
     &  lasize(1),lasize(2),glsize(1),glsize(2),
     &  gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &   D1(JOROG_GRAD_YY),
     &   atmos_landmask_local,
     &   lasize(1)*lasize(2),idummy)


      do nl=1,nelev
        CALL SCATTER_FIELD(LAND_LOC,
     &    icefrac_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &    icefrac_local(1,nl),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(LAND_LOC,
     &    landfrac_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &    landfrac_local(1,nl),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(LAND_LOC,
     &    nonice_snowmass_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &    nisnowmass_local(1,nl),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)



        CALL SCATTER_FIELD(LAND_LOC,
     &    nonice_snowdepth_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(LAND_LOC,lasize(1),lasize(2),offx,offy,1)
        CALL TO_LAND_POINTS(LAND_LOC,
     &    nisnowdepth_local(1,nl),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)

        CALL SCATTER_FIELD(icehflux_local(1,1,nl),
     &    icehflux_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
        CALL SWAPBOUNDS(icehflux_local(1,1,nl),
     &                  lasize(1),lasize(2),offx,offy,1)



      enddo

!+seg From Robin 31-5-17 REFAC1.46
! Major changes/new section
      CALL update_non_ice_snow(land_field, ntiles, nelev,
     &                               nsmax,dzsnow,
     &                               D1(JNSNOW(1)),
     &                               D1(JSNODEP_TYP),
     &                               D1(JSNODEP_TYP_LC), !29Oct17
     &                               D1(JSNOWDEPTH(1)),
     &                               D1(JRHO_SNOW_GRND(1)),
     &                               D1(JRGRAIN_TYP),
     &                               D1(JDS(1)),
     &                               D1(JRHO_SNOW(1)),
     &                               D1(JSICE(1)),
     &                               D1(JSLIQ(1)),
     &                               D1(JRGRAINL(1)),
     &                               D1(JTSNOWLAYER(1)),
     &                               D1(JTSTAR_TYP),
     &                               nisnowdepth_local,
     &                               glimmer_ice_density,
     &                               glimmer_ice_grain,
     &                               glimmer_ice_temp)



      CALL update_tile_fractions(land_field,nelev,ntype,
     &                       D1(JFRAC_TYP),
     &                       landfrac_local,
     &                       icefrac_local,
     &                       nnilandfrac_local)

! Old method for doing this inline deleted
!-seg



      write(6,*)"done new areas, doing conservation"
      !DO CONSERVATION
      !NEED A GLOBAL nnilandfrac_local - at this point easiest to just
      !GATHER WHAT WE'VE ALREADY GOT!

      do nl=1,nelev
      nnilandfrac_global(:,:,nl)=0.
        CALL FROM_LAND_POINTS(LAND_LOC,
     &    nnilandfrac_local(1,nl),
     &    atmos_landmask_local,
     &    lasize(1)*lasize(2),idummy)
        CALL GATHER_FIELD(LAND_LOC,
     &    nnilandfrac_global(1,1,nl),
     &    lasize(1),lasize(2),glsize(1),glsize(2),
     &    gather_pe,GC_ALL_PROC_GROUP,info)
      end do

      if (mype == gather_pe) then

        delta_vol=delta_vol*GLIMMER_ICE_DENSITY

        delta_nis=0.
        do nl=1,nelev
          do j=1,glsize(2)
          do i=1,glsize(1)
          if (nnilandfrac_global(i,j,nl).lt.0)
     &             nnilandfrac_global(i,j,nl)=0.
          if (nnilandfrac_global(i,j,nl).gt.1) then
            write(6,*)"nnifl oops:",nnilandfrac_global(i,j,nl)
          end if
!change in NIS
          delta_nis=delta_nis+(
     &                  GLIMMER_ICE_DENSITY
     &                * nonice_snowdepth_global(i,j,nl)
     &                * land_area_global(i,j)
     &                * nnilandfrac_global(i,j,nl)
     &                )
!          if (abs(nonice_snowdepth_global(i,j,nl)).gt.0) then
!            write(6,*)"nis:",nonice_snowdepth_global(i,j,nl)
!     &,nnilandfrac_global(i,j,nl)
!          endif
          end do
          end do
        end do

        ice_runoff_total=0.
        do j=1,glsize(2)
        do i=1,glsize(1)
!I've plumbed GBM fw_out into my calve_out, which is simply m/yr (water)
!not mm/s (water) as documented for the usual Glimmer runoff
          ice_runoff_total=ice_runoff_total+(
     &                  waterout_global(i,j)
     &                * land_area_global(i,j)
     &                * 1000.
     &                )
        end do
        end do

       write(6,*)"TOTAL WATER CONSERVATION, kg/yr: "
       conserv_global=delta_vol+delta_nis+ice_runoff_total-ice_smb_total
       write(6,'(a,e11.5)')"Error (-ve means can increase runoff): "
     &  ,conserv_global
       write(6,'(a,e11.5,f10.3)')"GCM SMB: ",ice_smb_total
     &                      ,conserv_global/ice_smb_total*100
       write(6,'(a,e11.5,f10.3)')"DELTA ICESHEET ",delta_vol
     &    ,conserv_global/delta_vol*100
       write(6,'(a,e11.5,f10.3)')"DELTA NI SNOW ",delta_nis
     &    ,conserv_global/delta_nis*100
       write(6,'(a,e11.5,f10.3)')"ICE CALVING ",ice_runoff_total
     &                         ,conserv_global/ice_runoff_total*100.

        ice_runoff_A=0.
        conserv_global=conserv_global*-1
        do j=1,glsize(2)
        do i=1,glsize(1)
          waterout_global(i,j)=waterout_global(i,j)+waterout_global(i,j)
     &                       *(conserv_global/ice_runoff_total)

          ice_runoff_A=ice_runoff_A+(
     &                  waterout_global(i,j)
     &                * land_area_global(i,j)
     &                * 1000.
     &                )
        end do
        end do
       conserv_global=delta_vol+delta_nis+ice_runoff_A-ice_smb_total
       write(6,'(a,2e11.3)')"ADJUST ICE, CONSERV ",ice_runoff_A
     &                         ,conserv_global

c! NH icecalv pattern we're scaling integrates up to 1Sv. Need to 
c! change units again
!+seg 04/19 Use pre-calculated forcing_time_step
       iceberg_scale=ice_runoff_A/1000./forcing_time_step/1e6
!-seg

       end if

       call gc_rbcast(1,1,0,nproc,info,iceberg_scale)
       write(6,*)"NH ICEBERG SCALING FACTOR is",iceberg_scale

       write(6,*) "finishing scatter ice fields"

      return



      END SUBROUTINE scatter_ice_fields

!+seg From Robin 31-5-17 REFAC1.48
! Completely new section

      SUBROUTINE update_non_ice_snow(land_field, ntiles, nelev,
     &                               nsmax,dzsnow,
     &                               nsnow,
     &                               snow_tile,
     &                               snow_tile_last_couple,  !29Oct17
     &                               snowdepth,
     &                               rho_snow_grnd,
     &                               rgrain,
     &                               ds,
     &                               rho_snow,
     &                               sice,
     &                               sliq,
     &                               rgrainl,
     &                               tsnow,
     &                               tstar_tile,
     &                               nisnowdepth_local_elev,
     &                               glimmer_ice_density,
     &                               glimmer_ice_grain,
     &                               glimmer_ice_temp)

      INTEGER :: land_field,ntiles,nelev,nsmax

      REAL :: dzsnow(nsmax)

      REAL :: nsnow(land_field,ntiles),
     &        snow_tile(land_field,ntiles),
     &        snow_tile_last_couple(land_field,ntiles),
     &        snowdepth(land_field,ntiles),
     &        rho_snow_grnd(land_field,ntiles),
     &        rgrain(land_field,ntiles),
     &        tstar_tile(land_field,ntiles)

      REAL :: ds(land_field,ntiles,nsmax),
     &        rho_snow(land_field,ntiles,nsmax),
     &        sice(land_field,ntiles,nsmax),
     &        sliq(land_field,ntiles,nsmax),
     &        rgrainl(land_field,ntiles,nsmax),
     &        tsnow(land_field,ntiles,nsmax)

      REAL :: nisnowdepth_local_elev(land_field,nelev)

      REAL :: glimmer_ice_density, glimmer_ice_grain, glimmer_ice_temp

!LOCAL
      INTEGER :: i,l,nl

      REAL :: dsice, newmass, oldmass
!+ seg 15/6/17
!!      write (6,*) 'GVALS ',glimmer_ice_density, glimmer_ice_grain,
!!     &  glimmer_ice_temp
! DEAL WITH CHANGES TO NON-ICE SNOW THAT GLIMMER MAY HAVE MADE
! BY CREATING/DELETING ICESHEET POINTS
      DO l=1,land_field
        DO n=1,ntiles-nelev
          nl=mod(n-1,nelev)+1
          !nisnowdepth is now just the /changes/ that Glimmer has made
          dsice=nisnowdepth_local_elev(l,nl)
          IF (abs(dsice) .gt. 1e-9) THEN
             IF (snowdepth(l,n)+dsice .LE. dzsnow(1)) THEN
               !We're not multilayer snow - at least, not now
               nsnow(l,n)=0.
               !NSNOW=0 has no depth
               snowdepth(l,n)=0.
               oldmass=snow_tile(l,n)
               newmass=dsice*glimmer_ice_density
               IF (dsice.gt.0) THEN
                 !adding ice to nsnow=0 snow
                 snow_tile(l,n)=max(0.,oldmass+newmass)
                 rgrain(l,n)=rgrain(l,n)
     &                               *(oldmass/(oldmass+newmass)) +
     &                       glimmer_ice_grain
     &                               *(newmass/(oldmass+newmass))
                 rho_snow_grnd(l,n)=rho_snow_grnd(l,n)
     &                              *(oldmass/(oldmass+newmass)) +
     &                       glimmer_ice_density
     &                              *(newmass/(oldmass+newmass))
               ELSE
               !subtracting mass from nsnow=0 snow
                 !newmass=dsice*rho_snow_grnd(l,n)
                 !still want to subtract at Glimmer density, no?
                 snow_tile(l,n)=max(0.,oldmass+newmass)
               END IF

             ELSE
               !We're dealing with multilayer snow
               !Will need to call a version of relayersnow to deal with
               !the changes to the snowpack
               IF (nsnow(l,n).EQ.0.) THEN
                 !initialise some stuff from the nsnow=0 GBMs first
                 nsnow(l,n)=1
                 ds(l,n,1)=snow_tile(l,n)/rho_snow_grnd(l,n)
                 sice(l,n,1)=snow_tile(l,n)
                 sliq(l,n,1)=0.
                 rgrainl(l,n,1)=rgrain(l,n)
                 rho_snow(l,n,1)=rho_snow_grnd(l,n)
                 tsnow(l,n,1)=tstar_tile(l,n)
                 !+seg
                 ! 
                 if (isnan(tstar_tile(l,n)) .or. tstar_tile(l,n) < 150
     &               .or. tstar_tile(l,n) > 320) then
                   write(6,*) 'TSTAR is bad'
                 endif
               END IF

             END IF !single or multilayer changes
          END IF !|dsice| > 0
          !!IGNORE FULL GBM CHANGES HERE TOO - IF THERE'S SIGNIFICANT GL
          !INVOLVEMENT IN THE BOX THERE'S A LARGE MASS OF SNOW HERE SOME
          !SO WE'RE UNLIKELY TO CROSS THE small LWRAD TRESHOLD FOR A
          !SNOW-COVERED SURFACE. ONLY AFFECTS 1 TIMESTEP ANYWAY.
        end do
      END DO

      !a modified version of relayer snow for changing multilayer
      !snow on non-ice tiles/elevs
      CALL glim_altersnowlayer(
     &      land_field,ntiles,nelev,nsmax,dzsnow
     &     ,nsnow,nisnowdepth_local_elev,snow_tile
     &     ,snowdepth,rho_snow_grnd,rgrain
     &     ,ds,rho_snow,sice
     &     ,sliq,rgrainl,tsnow
     &     ,glimmer_ice_density,glimmer_ice_grain,glimmer_ice_temp
     &     ,.TRUE.
     &)

      !+SEG Changes from RS 29Oct17
      !We want to update the record of the nonice_snow at coupling here
      !so that the acceleration only takes account of how the UM climate
      !has evolved the snowpack
      DO n = 1,ntiles - nelev
        DO l = 1,land_field
            snow_tile_last_couple(l,n) = snow_tile(l,n)
        END DO !loop on points
      END DO !loop on types
      !-SEG

      RETURN

      END SUBROUTINE update_non_ice_snow

      SUBROUTINE update_tile_fractions(land_field,nelev,ntype,
     &                       tilefrac_type,
     &                       landfrac_local_elev,
     &                       icefrac_local_elev,
     &                       nnilandfrac_local_elev)

!INTENT IN
      INTEGER :: land_field, nelev, ntype

      REAL :: landfrac_local_elev(land_field,nelev)
      REAL :: icefrac_local_elev(land_field,nelev)

!INTENT IN OUT
      REAL tilefrac_type(land_field,ntype)

!INTENT OUT
      REAL :: nnilandfrac_local_elev(land_field,nelev)

! LOCAL
      INTEGER :: l,nl, max_type

      REAL    :: glimfrac, sumfrac, max_frac
      REAL    :: fraction_threshold, gridbox_frac_disc
      REAL    :: icearea_local_elev(nelev)
      REAL    :: nicearea_local_elev(nelev)
      REAL    :: onilandfrac_local_elev(nelev)

! UPDATE UM TILE FRACTIONS WITH LAND/ICEFRAC. In splicing into the UM-sized
! field, Glimmer has essentially already multiplied icefrac /and/ landfrac
! by the total fraction of the UM box that falls onto the local grid (glimfrac)
! hence no glimfrac* factor on the arrays from Glimmer, just a (1-glimfrac)
! the original UM arrays, and the extra /glimfrac in the icefrac calc
! to avoid double counting it

      DO l=1,land_field

        glimfrac=0.
        DO nl=1,nelev
          nnilandfrac_local_elev(l,nl)=0.
        ! TOTAL FRAC OF THIS GRIDBOX THAT GLIMMER MAY HAVE CHANGED
        ! glimfrac=1 implies UM box is entirely in Glimmer domain
        ! glimfrac=0. is entirely outside
          glimfrac=glimfrac+landfrac_local_elev(l,nl)
        END DO

        IF (glimfrac.GT.0) THEN
        ! ONLY CONSIDER AREAS THAT MIGHT HAVE CHANGED
          !ICE FRACTIONS
          DO nl=1,nelev
            !glimmer ice fraction expressed in UM terms`
            icearea_local_elev(nl)=(icefrac_local_elev(l,nl)/glimfrac)*
     &                             landfrac_local_elev(l,nl)

            !UM ice fraction is Glimmer frac, possibly spliced into 
            !the pre-existing UM ice frac if Glimmer doesn't 
            !cover the whole gridbox
            tilefrac_type(l,ntype-nelev+nl)=icearea_local_elev(nl)+
     &                (1-glimfrac)*tilefrac_type(l,ntype-nelev+nl)
          END DO

          !NON_ICE-FRACTIONS. THESE HAVE TO BE UNAGGREGATED
          onilandfrac_local_elev(:)=0.
          DO n=1,ntype-nelev
            nl=mod(n-1,nelev)+1
            !UM's pre-existing aggregated non-ice fraction in this elevation 
            onilandfrac_local_elev(nl)=onilandfrac_local_elev(nl) +
     &                               tilefrac_type(l,n)
          END DO

          DO nl=1,nelev
            !glimmer non-ice fraction expressed in UM terms`
            nicearea_local_elev(nl)=landfrac_local_elev(l,nl) -
     &                              icearea_local_elev(nl)

            !UM new aggregated non-ice fraction is glimmer frac,
            !possibly spliced into the pre-existing frac if Glimmer 
            !doesn't cover the whole gridbox
            nnilandfrac_local_elev(l,nl)=nicearea_local_elev(nl) +
     &                (1-glimfrac)*onilandfrac_local_elev(nl)

            IF (nnilandfrac_local_elev(l,nl).GT.1)
            !SANITY CHECK
     &        write(6,*)"nnifl oops:",nnilandfrac_local_elev(l,nl)
          END DO

          !DIS-AGGREGATE, SOMEHOW
          DO n=1,ntype-nelev
            nl=MOD(n-1,nelev)+1

            IF (onilandfrac_local_elev(nl).gt.0) THEN
            !some non-ice fraction existed before. We can just
            !scale the non-ice fractions in each elevation up 
            !by the relative change
              tilefrac_type(l,n)=tilefrac_type(l,n) *
     &          nnilandfrac_local_elev(l,nl)/onilandfrac_local_elev(nl)
            ELSE
            !there was no non-ice fraction before. Proper initialisation is a
            !bit of an issue. For the fractions, assume it's all bare soil
              IF (nnilandfrac_local_elev(l,nl).gt.0)
     &                tilefrac_type(l,((soil_1-1)*nelev)+nl)=
     &                nnilandfrac_local_elev(l,nl)
            ENDIF
          END DO

        END IF !glimfrac

        !Final sanity check and rounding of fractions - we may have 
        !float(0.) issues from the various scalings above
        sumfrac=0.
        fraction_threshold=1e-4
        DO n=1,ntype
          !round off
          IF (tilefrac_type(l,n).LT.fraction_threshold)
     &        tilefrac_type(l,n)=0.
          sumfrac=sumfrac+tilefrac_type(l,n)
        END DO

        !do we need to make up for the things we've rounded off?
        gridbox_frac_disc=sumfrac -1.
        IF (abs(gridbox_frac_disc) .LT. 1e-9) THEN
          gridbox_frac_disc=0.
        ELSE IF ( gridbox_frac_disc .GT. 0.) THEN
          WRITE(6,*) "update_frac fail",l,sumfrac
          !+ seg
          ! Remove extra fraction
          max_frac=0.
          max_type=1.
          DO n=1,ntype
            IF (tilefrac_type(l,n).GT. max_frac) THEN
              max_frac=tilefrac_type(l,n)
              max_type=n
            END IF
          END DO
          tilefrac_type(l,max_type)=tilefrac_type(l,max_type)
     &                 -gridbox_frac_disc
          !- seg
        ELSE
        !simple: add the lost fractions to the tile with the most
          max_frac=0.
          max_type=1.
          DO n=1,ntype
            IF (tilefrac_type(l,n).GT. max_frac) THEN
              max_frac=tilefrac_type(l,n)
              max_type=n
            END IF
          END DO
          tilefrac_type(l,max_type)=tilefrac_type(l,max_type)
     &                 -gridbox_frac_disc
        END IF

      END DO !land_field

      RETURN

      END SUBROUTINE update_tile_fractions

!-seg  REFAC1.48


 
*DECK GLIM_SETSURF

      SUBROUTINE GLIM_SETSURF(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGDUMA
*CALL ARGSTS
*CALL ARGPTRA
     & icehflux_local,
     & ICODE )

        IMPLICIT NONE

*CALL CMAXSIZE     ! contains array size limits for other blocks
*CALL CSUBMODL     ! contains size information for CTIME and a_im
*CALL TYPSIZE      ! contains info for dynamical sizing of field arrays
*CALL TYPD1        ! contains information on D1 addressing array
*CALL TYPSTS       ! STASH related variables
*CALL TYPDUMA      ! DUMP headers for atmospheric model
*CALL TYPPTRA      ! D1 pointers for atmospheric model
*CALL PARVARS      ! contains size information about local & global 
                   !  fields & nproc
*CALL CTIME        ! contains time information
*CALL CAOPTR       ! common block contains pointer JA_FASTRUNOFF
*CALL C_MDI        ! common block containing missing data indicators
*CALL GCCOM        ! allows cross processor talking
*CALL NSTYPES      ! surface/elevation type sizes

*CALL AMAXSIZE
*CALL ATM_LSM

*CALL SOIL_THICK

      integer :: i,j,n,nl,idummy,l,icode

      integer :: lice_pts,lice_index(land_field),land_index(p_field)

      real :: orog_land(land_field)
      real :: surf_ht_flux(land_field)
      real :: orog_local(lasize(1)*lasize(2))
      real :: icehflux_local(lasize(1)*lasize(2),nelev)
      real :: elevation(nelev)
      real :: sumifrac(land_field)

!+seg These are mid-points. Changed to 10 levels (CESM type) from 25
      data elevation / 100., 300., 550., 850., 1150.,
     &                 1450.,1800.,2250.,2750.,3600./


      write(6,*)"in setsurf fields",land_field,lasize(1)*lasize(2)
      lice_pts=0
      do i=1,land_field
        orog_land(i)=0
        do n=1,ntype
          nl=mod(n-1,nelev)+1
          orog_land(i)=orog_land(i)
     &                 +(D1(JFRAC_TYP+(n-1)*LAND_FIELD+i-1)
     &                 *elevation(nl))
        end do
        sumifrac=0.
        do n=ntype-nelev,ntype-1
          SUMIFRAC(i)=SUMIFRAC(i)+D1(JFRAC_TYP+(N*LAND_FIELD)+I-1)
        end do
        if (sumifrac(i).gt.1e-6) then
          lice_pts=lice_pts+1
          LICE_INDEX(LICE_PTS)=I
        endif
      end do

      l=0
!      orog_local(:)=-99
      do i=1,p_field
        if (D1(JLAND+I-1).ne.0) THEN
          l=l+1
          land_index(l)=i
!+seg
! This is already reset directly from GLIMMER.
!          D1(JOROG+i-1)=orog_land(l)
!-seg
          

          surf_ht_flux(l)=0.
          if (sumifrac(l).gt.0) then
          do nl=1,nelev
            surf_ht_flux(l)=surf_ht_flux(l)+
     &                      icehflux_local(i,nl)*
     &                D1(JFRAC_TYP+((NTYPE-NELEV+NL-1)*LAND_FIELD)+l-1)/
     &                      sumifrac(l)
          end do
          endif
        endif
      end do

      !!!AGGREGATE ICESHEET SURFACE HEAT FLUXES INTO A GBM USING FRACS
      !!!DO A CALL TO THE MODIFIED ICE_HTC WITH THESE AS SURF_HT_FLUX

      !!! Remove, as per discussion with Robin (22/6/2015)
      !!! Investigate later
      !!call ice_htc(land_field,sm_levels,lice_pts,LICE_INDEX,DZSOIL,
      !!&             surf_ht_flux,60.*60.*24.*360.,
      !!&             D1(J_DEEP_ICE_TEMP),.TRUE.,.FALSE.)
        
      write(6,*) "Call to ice_htc commented out"
      write(6,*) "finishing setsurf after glint"

      return

      END SUBROUTINE

*DECK GLIM_ALTERSNOWLAYER
! *****************************COPYRIGHT*******************************
! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
! *****************************COPYRIGHT*******************************

! Description:
!     Version of relayersnow designed to cope with 
!     changes at the _bottom_ of the snowpack as mass is handed to/from
!     an icesheet model.
!     Relayersnow says:
!     Redivide snowpack after changes in depth, conserving mass and ener

! Subroutine Interface:
SUBROUTINE glim_altersnowlayer(land_field,ntiles,nelev,nsmax,dzsnow   & 
          ,nsnow,ni_snowdepth,snowmass,snowdepth,rho_snow_grnd,rgrain &
          ,ds,rho_snow,sice,sliq,rgrainl,tsnow                        &
          ,glimmer_ice_density,glimmer_ice_grain,glimmer_ice_temp     &
          ,l_snow_albedo  ) 


USE c_0_dg_c, ONLY :                                              &
!  imported scalar parameters
 tm   !   temperature at which fresh water freezes and ice melts (K)

USE c_perma, ONLY :                                               &
!  imported scalar parameters
 hcapi                                                            &
            !  Specific heat capacity of ice (J/kg/K)
,hcapw      !  Specific heat capacity of water (J/kg/K)

USE rad_param, ONLY :                                             &
!  imported scalars with intent(in)
 r0         !  Grain size for fresh snow (microns)

IMPLICIT NONE


INTEGER :: &
 land_field& 
,ntiles    &
,nelev    &
,nsmax    

REAL :: &
 dzsnow(nsmax)                  &
,ni_snowdepth(land_field,nelev)        &
!
,nsnow(land_field,ntiles)       &
,snowmass(land_field,ntiles)    &
,snowdepth(land_field,ntiles)   &
,rho_snow_grnd(land_field,ntiles) &
,rgrain(land_field,ntiles) &
,ds(land_field,ntiles,nsmax)    &
,rho_snow(land_field,ntiles,nsmax)      &
,sice(land_field,ntiles,nsmax)  &
,sliq(land_field,ntiles,nsmax)  &
,rgrainl(land_field,ntiles,nsmax)       &
,tsnow(land_field,ntiles,nsmax) &
,glimmer_ice_density    &
,glimmer_ice_grain      &
!,glimmer_ice_temp(land_field,nelev)
,glimmer_ice_temp

LOGICAL :: &
 l_snow_albedo

! LOCAL
REAL, PARAMETER :: thin_snow_limit = 1.0e-12
                  ! Maximum snow thickness (m) that is neglected
                  ! during relayering. All contributions
                  ! (mass, energy etc) from that snow are
                  ! neglected.


INTEGER :: &
 l,n,k,nl  &
 ,nold,old,new  &
 ,iznew,izz

REAL :: &
 d0(nsmax)      &
,dsice     &
,oldremains     &
,newremains(nsmax)     &
,csnow  &
,e(nsmax)       &
,u(nsmax)       &
,r(nsmax)       &
,s(nsmax)       &
,w(nsmax)       &
,mass_n &
,wt
 

!We'll be passing in fields that are (landfield,ntiles[,nsmax]) size. 
!We want to narrow it down to: real tiles (tile index may be out of date!)                             
!                              |dsice| > 0
!                              nsnow > 0


 DO l=1,land_field
DO n=1,ntiles-nelev
 nl=mod(n-1,nelev)+1

 IF (nsnow(l,n).gt.0                                                   &
   .AND. abs(ni_snowdepth(l,nl)).gt.1e-6                               &
    ) THEN

   nold = nsnow(l,n)
   d0(1:nold) = ds(l,n,1:nold)

  ! use ice density for incoming changes, pack density
  ! for outgoing. Used to be lowest layer density, but
  ! this can result in more mass being taken (seg, Sep/17)
  if (ni_snowdepth(l,nl) .gt. 0) then
    d0(nold) = d0(nold) + ni_snowdepth(l,nl)
    dsice=ni_snowdepth(l,nl)*GLIMMER_ICE_DENSITY
  else
    d0(nold) = d0(nold) + ni_snowdepth(l,nl)
    dsice=ni_snowdepth(l,nl)*rho_snow_grnd(l,n) ! SEG: RS 4/Sep/17
  endif

!-----------------------------------------------------------------------
! Calculate tile avg. snowdepth
!-----------------------------------------------------------------------
  snowmass(l,n)=max(snowmass(l,n)+dsice,0.) ! SEG: RS 4/Sep/17

! 1D version of layer snow
  nsnow(l,n) = 0  ! SEG: RS 4/Sep/17
  ds(l,n,:) = 0.0 ! SEG: RS 4/Sep/17
  snowdepth(l,n) = 0.

  if (snowmass(l,n) .gt. 1e-6) then ! SEG: RS 4/Sep/17
  do k=1,nold
    snowdepth(l,n) = snowdepth(l,n) + d0(k)
  end do

!-----------------------------------------------------------------------
! Divide snowpack into new layers
!-----------------------------------------------------------------------
! 1D version of layer snow
  nsnow(l,n) = 0
  ds(l,n,:) = 0.0

!   Only divide into layers if depth is >= a threshold.
! layersnow 1D
  oldremains = snowdepth(l,n)

  do k=1,nsmax
    ds(l,n,k) = dzsnow(k)
    oldremains = oldremains - dzsnow(k)
    if ( oldremains <= dzsnow(k) .OR. k == nsmax) then
      ds(l,n,k) = ds(l,n,k) + oldremains
      oldremains=0
      exit
    end if
  end do
  nsnow(l,n) = k


!-----------------------------------------------------------------------
! Store previous snow layer energy and mass contents
!-----------------------------------------------------------------------
  do k=1,nold
    csnow = sice(l,n,k)*hcapi + sliq(l,n,k)*hcapw
    e(k) = csnow * ( tsnow(l,n,k) - tm )
    r(k) = rgrainl(l,n,k)
    s(k) = sice(l,n,k)
    w(k) = sliq(l,n,k)
  end do

  if (dsice .lt. 0 ) then
  ! subtract mass from layers until you have something left
    oldremains=dsice
    do k=nold,1,-1
      mass_n=s(k)+w(k)
      if (mass_n+oldremains .gt. 0 .and. ds(l,n,k).gt.0) then
          s(k)=s(k)+s(k)*oldremains/(mass_n)
          w(k)=w(k)+w(k)*oldremains/(mass_n)
          csnow = s(k)*hcapi + w(k)*hcapw
          e(k) = csnow * ( tsnow(l,n,k) - tm )
          EXIT
      else
        s(k)=0.
        w(k)=0.
        e(k)=0.
        oldremains=oldremains+mass_n
      endif
    end do
    u(:) = e(:)
    sice(l,n,:) = s(:)
    sliq(l,n,:) = w(:)
    rgrainl(l,n,:) = r(:)

  else

  ! add the new details to the old bottom layer for rearrangement
  k=nold
  csnow = dsice * hcapi
  r(k)=r(k)*s(k)/(s(k)+ dsice) + &
       glimmer_ice_grain*dsice/(s(k)+dsice)
  s(k)=s(k)+dsice
  e(k)=e(k)+csnow * ( glimmer_ice_temp - tm )


!-----------------------------------------------------------------------
! Initialise accumulations for new layer values.
!-----------------------------------------------------------------------
  u(:) = 0.
  sice(l,n,:) = 0.
  sliq(l,n,:) = 0.
  rgrainl(l,n,:) = 0.

!-----------------------------------------------------------------------
! Set the state of the new layers.
!-----------------------------------------------------------------------
!    Initialise with all new layers empty.
  newremains(1:nsnow(l,n)) = ds(l,n,1:nsnow(l,n))
!     Start by filling top new layer.
  iznew = 1

!     Loop over the old layers.
    do old=1,nold

!       All of this old layer remains to be reassigned to new layer(s).
      oldremains = d0(old)

!       Point to first new layer with remaining space.
      izz = iznew

!       Loop over new layers with remaining space.
      do new=izz,nsnow(l,n)

        if ( oldremains > newremains(new) ) then
!-----------------------------------------------------------------------
! The remaining depth in the new layer will be exhausted by some or
! all of the remaining depth from the old layer.
!-----------------------------------------------------------------------

!           Decrement old layer by the remaining space in new layer.
          oldremains = oldremains - newremains(new)

!           Add properties from old layer to accumulation for new layer.
!           Note that wt is <= 1 since here we have oldRemains>newRemain
!           and oldRemains <= d0.
          if ( d0(old) > thin_snow_limit ) then
            wt =  newremains(new) / d0(old)
            u(new) = u(new) + e(old) * wt
            sice(l,n,new) = sice(l,n,new) + s(old) * wt
            sliq(l,n,new) = sliq(l,n,new) + w(old) * wt
            if (l_snow_albedo) rgrainl(l,n,new) = rgrainl(l,n,new) +  &
                             r(old) * newremains(new)
          end if

!           Update the pointer to the next new layer with space.
          izz = new + 1

        else

!-----------------------------------------------------------------------
! The old layer will be exhausted by this increment.
!-----------------------------------------------------------------------
!           Decrement available space in the new layer.
          newremains(new) = newremains(new) - oldremains
!           Add properties from old layer to accumulation for new layer.
          if ( d0(old) > thin_snow_limit ) then
            wt = oldremains /  d0(old)
            u(new) = u(new) + e(old) * wt
            sice(l,n,new) = sice(l,n,new) + s(old) * wt
            sliq(l,n,new) = sliq(l,n,new) + w(old) * wt
            if (l_snow_albedo) rgrainl(l,n,new) = rgrainl(l,n,new) + &
                                                r(old) * oldremains
          end if
!           Proceed to the next old layer by exiting from the new layer 
          exit
        end if
      end do  !  new layers
!       Update pointer to the next new layer with space.
      iznew = izz
    end do  !  old layers
    do k=1,nsnow(l,n)
      if (l_snow_albedo) rgrainl(l,n,k) = rgrainl(l,n,k) / ds(l,n,k)
    end do

    end if !(sice gt 0)

!-----------------------------------------------------------------------
! Diagnose layer temperatures and densities.
!-----------------------------------------------------------------------
    do k=1,nsnow(l,n)
      csnow = sice(l,n,k)*hcapi + sliq(l,n,k)*hcapw
      tsnow(l,n,k) = tm + u(k) / csnow
      rho_snow(l,n,k) = ( sice(l,n,k) + sliq(l,n,k) ) / ds(l,n,k)
    end do
    end if !(mass >0 safety check) ! SEG: RS 4/Sep/17

!-----------------------------------------------------------------------
! Snow surface grain size for radiative calculations
!-----------------------------------------------------------------------
    if (l_snow_albedo) rgrain(l,n) = rgrainl(l,n,1)

!-----------------------------------------------------------------------
! Diagnose bulk density of pack.
!-----------------------------------------------------------------------
!+SEG: RS 4/Sep/17
    if (snowdepth(l,n).gt.0) then
      rho_snow_grnd(l,n) = snowmass(l,n) / snowdepth(l,n)
    else
      rho_snow_grnd(l,n) = 250.
    endif
!-SEG: RS 4/Sep/17


!-----------------------------------------------------------------------
! Set values for unused snow layers.
! Note: not needed for algorithm, but clearer to follow.
!-----------------------------------------------------------------------
  if (  nsnow(l,n) < nsmax ) then
    k = nsnow(l,n) + 1
    if (l_snow_albedo) rgrainl(l,n,k:) = r0
    !rho_snow(l,n,k:) = 0.0
    rho_snow(l,n,k:) = 250.0 ! SEG: RS 4/Sep/17
    sice(l,n,k:) = 0.
    sliq(l,n,k:) = 0.
    tsnow(l,n,k:) = tm
  end if

 END IF !VALID POINT

 END DO  !  L (points)
END DO  !  N (TILES/ELEVS)

END SUBROUTINE glim_altersnowlayer

*DECK GLIM_IC

      SUBROUTINE GLIM_INTCTL(
     &               nx,ny,nelev,
!
     &               htime,
     &               ice_smb_global,
     &               ice_areas_global,
     &               nonice_snowdepth_global,
     &               icestemp_global,
!
     &               landfrac_global,
     &               icefrac_global,
!+seg 04/19 New field
     & orogMN_global,
!-seg
     & orogSD_global,orogHO_global,orogSIL_global,
     & orogXX_global,orogXY_global,orogYY_global,
     &               waterout_global,
!     &               icerunoff_global,
!     &               liqrunoff_global,
     &               icehflux_global
     &               ,delta_vol
     &                      )

      use glim_famtype
      use glimmer_paramets

      implicit none

      integer nx,ny,nelev
      integer xpts, ypts
      integer htime

!+seg 04/19 Number of masks
      integer,parameter :: nmask = 3
!-seg
!+seg 31/aug/17 
      integer :: nylim
!-seg


!SINGLE PE VERSIONS OF THINGS TO GO INTO GLIMMER
      real(kind=8),intent(inout) :: ice_smb_global(nx,ny,nelev)
      real(kind=8),intent(inout) :: ice_areas_global(nx,ny,nelev)
      real(kind=8),intent(inout) :: nonice_snowdepth_global(nx,ny,nelev)
      real(kind=8),intent(in) :: icestemp_global(nx,ny,nelev)

!SINGLE PE VERSIONS OF THINGS THAT COME OUT OF GLIMMER
      real(kind=8),intent(inout) :: landfrac_global(nx,ny,nelev)
      real(kind=8),intent(inout) :: icefrac_global(nx,ny,nelev)
!+seg 04/19 New field
      real(kind=8),intent(inout) :: orogMN_global(nx,ny)
!-seg
      real(kind=8),intent(inout) :: orogSD_global(nx,ny)
      real(kind=8),intent(inout) :: orogHO_global(nx,ny)
      real(kind=8),intent(inout) :: orogSIL_global(nx,ny)
      real(kind=8),intent(inout) :: orogXX_global(nx,ny)
      real(kind=8),intent(inout) :: orogXY_global(nx,ny)
      real(kind=8),intent(inout) :: orogYY_global(nx,ny)
      real(kind=8),intent(out) :: waterout_global(nx,ny)
!      real(kind=8),intent(out) :: icerunoff_global(nx,ny,nelev)
!      real(kind=8),intent(out) :: liqrunoff_global(nx,ny,nelev)
      real(kind=8),intent(out) :: icehflux_global(nx,ny,nelev)
      real(kind=8) :: delta_vol


      integer    :: i,j,l
!+seg Changed from 25 to nelev, mar 2016
      real       :: elevation(10),glimfrac

!+seg 04/19 tmp_orog_mn now used (rather than dum_orog_mn)
      real(kind=8) :: tmp_orog_mn(nx,ny),
!-seg
     &                tmp_orog_SD(nx,ny),
     &                tmp_orog_SIL(nx,ny),
     &                tmp_orog_H(nx,ny),
     &                tmp_sigma_XX(nx,ny),
     &                tmp_sigma_XY(nx,ny),
     &                tmp_sigma_YY(nx,ny),
     &                dum_icemask(nx,ny)

!+seg 31/Aug/17
      real(kind=8) :: smb_leva(nelev)

      integer(kind=4)    :: Fmask(nx,ny)
      integer(kind=4)    :: Fmask3D(nmask,nx,ny,nelev)
!-seg

      logical    :: tmp_landmask(nx,ny)

!+seg These are mid-points. Changed to 10 levels (CESM type) from 25
      data elevation / 100., 300., 550., 850., 1150.,
     &                 1450.,1800.,2250.,2750.,3600./

      write(6,*)"in glim_intctl"
!
      ice_smb(:,:,:)=ice_smb_global(:,:,:)
      ice_areas(:,:,:)=ice_areas_global(:,:,:)
      icestemp(:,:,:)=icestemp_global(:,:,:)
      nonice_snowdepth(:,:,:)=nonice_snowdepth_global(:,:,:)

c     !ice_smb(:,:,:)=icefrac_global(:,:,:)
c!      ice_smb(:,:,:)=0.
c!     nonice_snowdepth(:,:,:)=0.

! Glimmer doesn't work with area-fractions-per-elevation like
! we do - it has fixed elevation classes, but still expects to 
! be given the mean orographic height of the land in each class. 
! For us, that's just the midpoint of the class right now
      do l=1,nelev
        faketopo(:,:,l)=elevation(l)
      end do

!+seg 23/6/17
! Create masks. This could be done completely on GLIMMER side
! perhaps.
      call create_masks(nx,ny,nelev,icefrac_global,landfrac_global,
     &                          Fmask,Fmask3D,htime)
!-seg

!+seg  31/Aug/17
! Create average SMB for each level (for infill in interpolation)
! If this was a global analysis, would be tricky (maybe)
! Create a NH value
      nylim=ny/2
      call create_smb_leva(nx,ny,nylim,nelev,Fmask3D(0,:,:,:),
     &  ice_smb_global,smb_leva)
!-seg

      delta_vol=ice_vol
      UM_time=htime

      write(6,*)"GOING INTO GLINT: non-ice snow turned off (gsdep)"
!          !! The making of new ice points from nonice snowdepth
!          !! on the Glimmer grid has to happen AFTER 3D interp
!          !! of the snowdepth field, so will be in GLINT
          CALL GLINT_GCM(
!! REQ'D FIELDS
     &               ice_sheet
     &               ,UM_time
!! OPTIONAL FIELDS, IN
     &               ,qsmb=ice_SMB
     &               ,gareas=ice_areas
     &               ,tsfc=icestemp
     &               ,topo=faketopo
!! OPTIONAL FIELDS, INOUT
     &               ,gfrac=icefrac
     &               ,glfrac=landfrac
!+seg In this mod, nonice snow turned off 
!     &               ,gsdep=nonice_snowdepth
!-seg
!! OPTIONAL FIELDS, OUT
     &               ,gcalv=fw_out
!     &               ,grofi=icerunoff
!     &               ,grofl=liqrunoff
     &               ,ghflx=icehflux
     &               ,ice_volume=ice_vol
!+SEG 2/3/17
!+SEG 22/6/17 - Need to define these masks using ice fraction
! Fmask - mask from column sum FMask3D IN or INOUT?
     &               ,Fmask=Fmask
     &               ,Fmask3D=Fmask3D
!-SEG
!+seg 31/Aug/17
     &               ,qsmb_lev=smb_leva
!-seg
     &               )
      write(6,*)"BACK FROM GLINT!"

      icefrac_global(:,:,:)=icefrac(:,:,:)
      landfrac_global(:,:,:)=landfrac(:,:,:)
      nonice_snowdepth_global(:,:,:)=nonice_snowdepth(:,:,:)
!      icerunoff_global(:,:,:)=icerunoff(:,:,:)
!      liqrunoff_global(:,:,:)=liqrunoff_global(:,:,:)
      icehflux_global(:,:,:)=icehflux(:,:,:)
      waterout_global(:,:)=fw_out(:,:)
      delta_vol=ice_vol-delta_vol

      write(6,*)"glim_ic delta icevol:",delta_vol,sum(waterout_global)

!!      icerunoff_global(:,:,:)=0.
!!      liqrunoff_global(:,:,:)=0.
!       icehflux_global(:,:,:)=0.
!       waterout_global(:,:)=0.

      do i = 1 , ice_sheet%ninstances
        xpts = ice_sheet%instances(i)%model%general%ewn
        ypts = ice_sheet%instances(i)%model%general%nsn

        write(6,*)"CALLING OROGRAPHY",i,xpts,ypts
c!        write(6,*)ice_sheet%instances(i)%model%projection%stere
        call orography(xpts,             ! IN   nos of ISM cols
     &                 ypts,             ! IN   nos of ISM rows
     &                 lat_dy,           ! IN   latitude increment
     &                 lon_dx,           ! IN   longitude increment
     &                 lat_0,            ! IN   latitude reference pt
     &                 lon_0,            ! IN   longitude reference pt
     &                 ny,               ! IN   number of UM rows
     &                 nx,               ! IN   number of UM columns
     &                 ice_sheet%instances(i)%model%geometry%usrf*thk0,
     &                 ice_sheet%instances(i)%model%geometry%thck*thk0,
                                         ! IN   ISM surface elevation
     &                 ice_sheet%instances(i)%model%projection,
                                         ! IN   ISM projection
     &                 ice_sheet%instances(i)%lgrid,
                                         ! IN   ISM grid
!+seg 04/19 active orog_mn
     &                 tmp_orog_mn,      ! OUT                      |
!-seg
     &                 tmp_orog_SD,      ! OUT  Std.Dev. of orog    |
     &                 tmp_orog_SIL,     ! OUT  orog silhouette     |  on
     &                 tmp_orog_H,       ! OUT  orog peak-to-peak   |  UM
     &                 tmp_sigma_XX,     ! OUT  sigma_xx            | grid
     &                 tmp_sigma_XY,     ! OUT  sigma_xy            |
     &                 tmp_sigma_YY,     ! OUT  sigma_yy            |
     &                 dum_icemask)      !                          |
        write(6,*)"BACK FROM OROGRAPHY",i

!! NHEM stereographic projection produces artifacts at edges, for some 
!! reason. Easiest to cap for sanity, rather than delve into the morass above
!! JG fixing properly...?
        tmp_landmask(:,:)=.TRUE.
        where(orogSD_global .lt.-1e9) tmp_landmask=.FALSE.

        where(tmp_orog_SD   .ne. 0.0
     &  .AND. tmp_orog_SD  .lt. maxval(orogSD_global)*1e3 
     &  .AND. tmp_orog_SD  .gt. minval(orogSD_global,tmp_landmask)*1e3 
     &    )  orogSD_global  = tmp_orog_SD

        where(tmp_orog_sil  .ne. 0.0
     & .AND. tmp_orog_sil  .lt. maxval(orogSIL_global)*1e3 
     & .AND. tmp_orog_sil  .gt. minval(orogSIL_global,tmp_landmask)*1e3 
     &    )  orogSIL_global = tmp_orog_SIL

        where(tmp_orog_h    .ne. 0.0
     &  .AND. tmp_orog_h  .lt. maxval(orogHO_global)*1e3 
     &  .AND. tmp_orog_h  .gt. minval(orogHO_global,tmp_landmask)*1e3 
     &    )  orogHO_global   = tmp_orog_H

        where(tmp_sigma_XX  .ne. 0.0
     &  .AND. tmp_sigma_XX  .lt. maxval(orogXX_global)*1e3 
     &  .AND. tmp_sigma_XX  .gt. minval(orogXX_global,tmp_landmask)*1e3 
     &    )  orogXX_global = tmp_sigma_XX

        where(tmp_sigma_XY  .ne. 0.0
     &  .AND. tmp_sigma_XY  .lt. maxval(orogXY_global)*1e3 
     &  .AND. tmp_sigma_XY  .gt. minval(orogXY_global,tmp_landmask)*1e3 
     &    )  orogXY_global = tmp_sigma_XY

        where(tmp_sigma_YY  .ne. 0.0 
     &  .AND. tmp_sigma_YY  .lt. maxval(orogYY_global)*1e3 
     &  .AND. tmp_sigma_YY  .gt. minval(orogYY_global,tmp_landmask)*1e3 
     &       )  orogYY_global  = tmp_sigma_YY

!+seg 04/19 New field
        where(tmp_orog_mn   .ne. 0.0
     &  .AND. tmp_orog_mn  .lt. maxval(orogMN_global)*1e3
     &    )  orogMN_global  = tmp_orog_mn
!-seg

      end do

      RETURN

      END


      subroutine glim_write_restart(yr)
      use glim_famtype
      use glint_main
      integer yr ! this is potentially not an integer, rather the time stamp that the UM uses
      call glim_force_write_restart(ice_sheet,yr)
      return 
      end

!+seg 04/19 New routines for mask creation
      subroutine create_masks(nx,ny,nelev,icefrac_global,
     &           landfrac_global,Fmask,Fmask3D,htime)

!+
! Name:   create_masks
! Author: S.E.George
! Date:   23-Jun-17
!-
      implicit none

      integer,parameter   :: nmask=3
      integer,parameter   :: mask_ice=1  ! Ice fractional point
      integer,parameter   :: mask_snw=2  ! None-ice fractional point
      integer,parameter   :: mask_lnd=3  ! 3D land mask

      integer(kind=4),parameter   :: rtoi=1000000
      integer,intent(in)  :: nx     ! E-W points
      integer,intent(in)  :: ny     ! N-S points
      integer,intent(in)  :: nelev  ! Nof of elevations

      real(kind=8),intent(in) :: landfrac_global(nx,ny,nelev)
      real(kind=8),intent(in) :: icefrac_global(nx,ny,nelev)

      integer(kind=4),intent(out)     :: Fmask(nx,ny)
      integer(kind=4),intent(out)     :: Fmask3D(nmask,nx,ny,nelev)

      integer,intent(in)  :: htime
!     Local variable

      real(kind=8)        :: dlandfrac_global(nx,ny,nelev)
      real(kind=8)        :: dicefrac_global(nx,ny,nelev)
      integer(kind=4)     :: ilandfrac_global(nx,ny,nelev)
      integer(kind=4)     :: iicefrac_global(nx,ny,nelev)
!      integer(kind=4)     :: ilandfrac_global_sum(nx,ny)
!      integer(kind=4)     :: iicefrac_global_sum(nx,ny)

      integer             :: i, j, k
      integer             :: ice_frac
      integer             :: lnd_frac

! Convert real fractions to integer
! landfrac_global & icefrac_global have missing values (large 
! negative). Need to remove before conversion to integer (these
! are intent(in), so local

      dlandfrac_global=landfrac_global
      dicefrac_global =icefrac_global
      where(dlandfrac_global<0.0) dlandfrac_global = 0.0
      where(dicefrac_global<0.0) dicefrac_global   = 0.0
! Just to make sure, do the other end
      where(dlandfrac_global>1.0) dlandfrac_global = 1.0
      where(dicefrac_global>1.0) dicefrac_global   = 1.0

      ilandfrac_global = dlandfrac_global * rtoi
      iicefrac_global  = dicefrac_global  * rtoi

      write(6,*) 'MaximinfracL',maxval(dlandfrac_global),
     &  maxval(ilandfrac_global),
     &  minval(dlandfrac_global),minval(ilandfrac_global)
      write(6,*) 'MaxminfracI',maxval(dicefrac_global),
     &  maxval(iicefrac_global),
     &  minval(dicefrac_global),minval(iicefrac_global)

! Create masks

      Fmask=0
      Fmask3D=0

      loop2: do i=1,nx
        loop3: do j=1,ny

!          if (iicefrac_global_sum(i,j) > 0) then
!             Fmask(i,j) = 1
!          endif
          loop4 : do k=1,nelev
            ice_frac=iicefrac_global(i,j,k)
            lnd_frac=ilandfrac_global(i,j,k)
            if (lnd_frac == 0) then
              if (ice_frac > 0) then
! Warning message if there is an ice fraction, but no land                
               write(6,*)
     &          'WARNING - create_masks: ice fraction with no land',
     &          i,j,k
              endif
! Set masks to zero
              Fmask3D(mask_ice,i,j,k)=0
              Fmask3D(mask_snw,i,j,k)=0
              Fmask3D(mask_lnd,i,j,k)=0
            else
              Fmask3D(mask_lnd,i,j,k)=1
              if (ice_frac >= rtoi) then ! i.e. ice fraction >=1.0
                Fmask3D(mask_ice,i,j,k)=1 ! All ice
                Fmask3D(mask_snw,i,j,k)=0 ! No non-ice
              else if (ice_frac == 0) then
                Fmask3D(mask_ice,i,j,k)=0 ! No ice
                Fmask3D(mask_snw,i,j,k)=1 ! All non-ice
              else ! Mixed
                Fmask3D(mask_ice,i,j,k)=1
                Fmask3D(mask_snw,i,j,k)=1
              endif
            endif

          end do loop4
          if (sum(Fmask3D(mask_ice,i,j,:)) > 0) then
            Fmask(i,j)=1
          else
            Fmask(i,j)=0
          endif

        end do loop3
      end do loop2

      end subroutine create_masks

      subroutine create_smb_leva(nx,ny,nylim,nelev,mask3D,smb3D,
     & smb_leva)

      implicit none

      integer,intent(in)  :: nx     ! E-W points
      integer,intent(in)  :: ny     ! N-S points
      integer,intent(in)  :: nelev  ! No of elevations
      integer,intent(in)  :: nylim  ! range of latitudes (from NP=1)

      integer(kind=4), intent(in) ::  mask3D(nx,ny,nelev)
      real(kind=8),intent(in)     ::  smb3D(nx,ny,nelev)

      real(kind=8),intent(out)    :: smb_leva(nelev)

      integer  :: i,j,k
      integer  :: icount
      real(kind=8) :: smb_sum


      do k=1,nelev
        icount=0
        smb_sum=0.
        do i=1,nx
          do j=1,nylim
            if (mask3D(i,j,k) == 1) then
              smb_sum=smb_sum+smb3D(i,j,k)
              icount=icount+1
            endif
          end do
        end do
        if (icount == 0) then
          smb_leva(k) =0.0
        else
          smb_leva(k)=smb_sum/icount
        endif
      end do

      write (6,*) 'create_smb_leva vals: smb_leva ',smb_leva

      end subroutine create_smb_leva
!-seg

!+seg 04/19 Accel code from Robin
      SUBROUTINE nonice_snow_accel(mype
     &                            ,gather_pe
     &                            ,nproc
     &                            ,land_field
     &                            ,ntiles
     &                            ,nelev
     &                            ,nsmax
     &                            ,snow_grain_tile
     &                            ,snow_mass_tile
     &                            ,snow_mass_tile_last_couple
     &                            ,snow_depth_tile
     &                            ,snow_rho_tile
     &                            ,snow_n_tile
     &                            ,snow_depth_layer
     &                            ,snow_ice_layer
     &                            ,snow_liq_layer
     &                            ,snow_rho_layer
     &                            ,snow_grain_layer
     &                            )

      USE glim_famtype
      USE glimmer_paramets

      IMPLICIT NONE

      INTEGER :: mype
     &          ,gather_pe
     &          ,nproc,info

      INTEGER :: ntiles
     &          ,nelev
     &          ,nsmax

      INTEGER :: land_field

      REAL(KIND=8) :: snow_grain_tile(land_field,ntiles)
     &               ,snow_mass_tile(land_field,ntiles)
     &               ,snow_mass_tile_last_couple(land_field,ntiles)
     &               ,snow_depth_tile(land_field,ntiles)
     &               ,snow_rho_tile(land_field,ntiles)
     &               ,snow_n_tile(land_field,ntiles)

      REAL(KIND=8) :: snow_depth_layer(land_field,ntiles,nsmax)
     &               ,snow_ice_layer(land_field,ntiles,nsmax)
     &               ,snow_liq_layer(land_field,ntiles,nsmax)
     &               ,snow_rho_layer(land_field,ntiles,nsmax)
     &               ,snow_grain_layer(land_field,ntiles,nsmax)

      INTEGER :: accel_factor

      REAL    :: mass_difference
     &         , accel_mass_diff

      INTEGER :: l,n

      IF (mype .EQ. gather_pe) THEN
        accel_factor = ice_sheet%instances(1)%ice_tstep_multiply
      END IF
      CALL gc_ibcast(1,1,0,nproc,info,accel_factor)

      !!!!!!!!!!!!!!SKIP ARTIFICALLY
      !accel_factor=0
      !!!!!!!!!!!!!!
      write(6,*)"rss",accel_factor
      IF (accel_factor .GT. 1) THEN
        DO n = 1,ntiles - nelev
          DO l = 1,land_field
          !This would be better over tile_ts, but we'd have to construct the
          !indexing anew ourselves from the frac_typs...
            IF (INT(snow_n_tile(l,n)) .EQ. nsmax
     &          .AND.(snow_mass_tile(l,n).GT.9000)
     &          ) THEN
            !We've got a lower layer to deal with, and ~10m of ice-density snow
              mass_difference = snow_mass_tile(l,n)
     &                        - snow_mass_tile_last_couple(l,n)


              accel_mass_diff = mass_difference * accel_factor

              IF ( accel_mass_diff .GT. 0 .OR.
     &            (accel_mass_diff .LT. 0 .AND.
     &            (abs(accel_mass_diff) .LT. snow_ice_layer(l,n,nsmax)))
     &           ) THEN
              !change the bottom layer. None of the other properties should
              !change(?)
               snow_ice_layer(l,n,nsmax)   = snow_ice_layer(l,n,nsmax)
     &                                      + accel_mass_diff

               snow_depth_layer(l,n,nsmax) = snow_depth_layer(l,n,nsmax) 
     &                     + accel_mass_diff / snow_rho_layer(l,n,nsmax)

               snow_mass_tile(l,n)  = snow_mass_tile(l,n)
     &                              + accel_mass_diff

               snow_depth_tile(l,n) = snow_depth_tile(l,n)
     &                     + accel_mass_diff / snow_rho_layer(l,n,nsmax)

              ELSE
              !we're removing mass and there's more than we have solid in the bottom
              !layer. Just nix the whole snowpack. rho and grain need non-zero
              !dummy values
                  snow_mass_tile(l,n)   = 0.
                  snow_depth_tile(l,n)  = 0.
                  snow_n_tile(l,n)      = 0
                  snow_rho_tile(l,n)    = 250.
                  snow_grain_tile(l,n)  = 50.

                  snow_ice_layer(l,n,:)   = 0.
                  snow_liq_layer(l,n,:)   = 0.
                  snow_depth_layer(l,n,:) = 0.
                  snow_rho_layer(l,n,:)   = 250.
                  snow_grain_layer(l,n,:) = 50.
              ENDIF !can we change bottom layer?

            ENDIF !nsmax layers

          END DO !loop on points
        END DO !loop on types
      ENDIF !are we accelerating?
      !makes sense to keep this field up to date in D1 even if we're not using it
      !in *this* run
      DO n = 1,ntiles - nelev
        DO l = 1,land_field
            snow_mass_tile_last_couple(l,n) = snow_mass_tile(l,n)
        END DO !loop on points
      END DO !loop on types

      END SUBROUTINE nonice_snow_accel
!-seg

*DECK GLIM_INIT

!----------------------------------------------------------------------!
!                                                                      !
! Modified from GLIM_AB_INITIO by RMG, 17/03/2009                      !
!----------------------------------------------------------------------!
subroutine glim_initialise(mype,            &
     A_LEN_REALHD, & ! IN
     A_REALHD,     & ! IN,OPT
     nx,           & ! IN
     ny,           & ! IN
     iiday,        & ! IN,OPT
     iimonth,      & ! IN,OPT
     iiyear        ) ! IN,OPT

  use glint_main
  use glimmer_paramets
  use glim_famtype
  use glimmer_config
  use glimmer_log

  use parallel
  use mpi_mod

  implicit none

  !-- arguments ---------------------------------!
  integer(kind=8),         intent(IN) :: mype
  integer(kind=8),         intent(IN) :: A_LEN_REALHD
  integer(kind=8),optional,intent(IN) :: nx, ny ! Global grid size 
  integer(kind=8),optional,intent(IN) :: iiday,iimonth,iiyear
  real(kind=8),optional,intent(in),dimension(A_LEN_REALHD):: A_REALHD

  !-- local variables --------------!  
  character(fname_length)             :: rst_filename 
  ! name of restart file
  logical                             :: rst_exist    
  ! restart file existence 

  integer ICECOMM,IERROR

  if (mype .eq. 0) then
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, 42, 0, ICECOMM, IERROR)
  else
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,0,ICECOMM,IERROR)
  endif

  if (mype .eq. 0) then

    call parallel_set_info(ICECOMM,0)

    start_mode = HOTSTART
    write(6,*)'temporary override: start_mode reset to hotstarts'

    select case(start_mode)

    case(HOTSTART)
      write(6,*)'GLIM_INITIALISE: calling glim_hotstart'
      call glim_hotstart(                                &
          A_LEN_REALHD,                                  &
          A_REALHD,                                      &
          nx,                                            &
          ny,                                            &
          iiday,                                         &
          iimonth,                                       &
          iiyear)
      write(6,*)'GLIM_INITIALISE: glim_hotstart finished'

    case default
      write(6,*)'GLIM_INITIALISE: unknown start_mode:',start_mode
      stop
 
    end select

  end if

end subroutine glim_initialise

*/
*/**********************************************************************
*/

!----------------------------------------------------------------------!
!                                                                      !
! Modified from GLIM_AB_INITIO by RMG, 17/03/2009                      !
!----------------------------------------------------------------------!
subroutine glim_finalise(mype,             &
     A_LEN_REALHD, & ! IN
     iiday,        & ! IN,OPT
     iimonth,      & ! IN,OPT
     iiyear        ) ! IN,OPT

  use glint_main
  use glimmer_paramets
  use glim_famtype
  use glimmer_config
  use glimmer_log

  implicit none

  !-- arguments ---------------------------------!
  integer(kind=8),         intent(IN) :: mype
  integer(kind=8),         intent(IN) :: A_LEN_REALHD
  integer(kind=8),optional,intent(IN) :: iiday,iimonth,iiyear

  !-- local variables --------------!  
  character(20)                       :: gnmlfile='glimmer.nml'


  if (mype .eq. 0) then
    write(6,*)'GLIM_FINALISE: Ending Glimmer.......'
  
    select case(start_mode)
    case(HOTSTART)
     !write out namelist file
     !tranfer dates vars into strings using internal file writing
      write(6,*)'GLIM_FINALISE: Ending hotstarted glimmer'

    case default
     write(6,*)'GLIM_FINALISE: start_mode not recognised:' &
          ,start_mode
     stop

    end select

    write(6,*)'GLIM_FINALISE: Terminating Glimmer'
    call end_glint(ice_sheet)

    write(6,*)'GLIM_FINALISE: Deallocating variables'
    deallocate(temp,precip)
    deallocate(lats,lons,orog)
    deallocate(coverage,orog_out,albedo, &
       ice_frac,fw_out)
    deallocate(cov_orog,fw_in)

    deallocate(ice_smb,ice_areas                  &
              ,icestemp,nonice_snowdepth,faketopo &
              ,landfrac,icefrac,icehflux          &
              )

    write(6,*)'GLIM_FINALISE: finished.'
    write(6,*)''

  end if

end subroutine glim_finalise

*/
*/**********************************************************************
*/

*DECK GLIM_HST

!----------------------------------------------------------------------!
!                                                                      !
! Modified from GLIM_AB_INITIO by RMG, 17/03/2009                      !
!----------------------------------------------------------------------!
subroutine glim_hotstart (  &
     A_LEN_REALHD,          & ! IN
     A_REALHD,              & ! IN
     nx,                    & ! IN
     ny,                    & ! IN
     iiday,                 & ! IN
     iimonth,               & ! IN
     iiyear)                  ! IN

  use glint_main
  use glimmer_paramets
  use glim_famtype
  use glimmer_config
  use glimmer_log

  implicit none

  !-- arguments ---------------------------------!
  integer(kind=8),intent(IN) :: A_LEN_REALHD
  integer(kind=8),intent(IN) :: nx, ny ! Global grid size 
  integer(kind=8),intent(IN) :: iiday,iimonth,iiyear

  real(kind=8),intent(in),dimension(A_LEN_REALHD):: A_REALHD

  !-- local variables --------------!  
  character(fname_length)             :: rst_filename 
  ! name of restart file
  logical                             :: rst_exist    
  ! restart file existence 

  character(fname_length)             :: paramfile='top.config'
  character(20)                       :: gnmlfile='glimmer.nml'
  integer                             :: i,j,forcing_time_step
  integer                             :: now_time,couple_this_year
  integer                             :: couple_next_year
  integer                             :: ICE_COUPLING_MONTH
  ! Array index counters

  !--------------------------------------------------------------------!
  !ojhb: Sort out loading of the variously named hotstart files.       !
  !ojhb: It requires consistent use of a given *TAG*:                  !
  !ojhb:                                                               !
  !ojhb:   1. a script renames hotstart filenames to:                  !
  !ojhb:                RUNID.*TAG*.hot.DATE.nc                        !
  !ojhb:   2. in the top.config we have..                              !
  !ojhb:        [GLINT instance]                                       !
  !ojhb:                       tag : *TAG*                             !
  !ojhb:   3. in instance specific configs we have..                   !
  !ojhb:              [CF input]      #hotstart input..                !
  !ojhb:                     tag : *TAG*_hotstart_in                   !
  !ojhb:              [CF output]     #hotstart output..               !
  !ojhb:                    name : *TAG*.hot.nc                        !
  !ojhb:                                                               !
  type(ConfigData),dimension(:),allocatable :: xtr_cnfgs               !
  type(ConfigSection),pointer               :: tmp_config=>null()      !
  type(ConfigSection),pointer               :: tmp_section=>null()     !
  integer                                   :: nos_instances           !
  character(22)                             :: instance_tag            !
  character(5)                              :: UM_runid
  character(len=100)                        :: message   ! log-writing !
  !--------------------------------------------------------------------!

!+SEG 1/3/17 New mask arrays for smoothing
  integer :: mask(nx,ny)
  integer :: mask3D(3,nx,ny,10)

  mask(:,:)=1
  mask3D(:,:,:,:)=1
!-SEG 1/3/17

  ! Allocate arrays appropriately
  lat_0=A_REALHD(3); lon_dx=A_REALHD(1)
  lon_0=A_REALHD(4); lat_dy=A_REALHD(2)

  call open_log(unit=220)
  write(6,*)'GLIM_HOTSTART: --------realhd when passed in----------'
  write(6,*)'GLIM_HOTSTART: lat_0=',A_REALHD(3), 'lon_dx=', &
       A_REALHD(1), 'lon_0=',A_REALHD(4), 'lat_dy=',A_REALHD(2)

! Insist on signs of increments.
  if (lon_dx<0 .or. lat_dy<0) then
    write(6,*) 'UM longitude and latitude steps must both be positive'
    call abort
  endif

  write(6,*)'GLIM_HOTSTART: A_LEN_REALHD=',A_LEN_REALHD
  write(6,*)'GLIM_HOTSTART: glsize(1)[nx]=',nx,'glsize(2)[ny]',ny
  write(6,*)'GLIM_HOTSTART: ---------------------------------------'
  write(6,*)'GLIM_HOTSTART: Allocating variables' 

  allocate(temp(nx,ny),     precip(nx,ny),     &
       lats(ny),        lons(nx),          &
       orog(nx,ny),     coverage(nx,ny),   &
       orog_out(nx,ny), albedo(nx,ny),     &
       ice_frac(nx,ny), fw_out(nx,ny),     &
       cov_orog(nx,ny), fw_in(nx,ny)       )

!+seg 04/19 Changes to memorry allocation.
!     Previously 25 elevation classes, now 10
  allocate( ice_smb(nx,ny,10), &
            ice_areas(nx,ny,10), &
            icestemp(nx,ny,10), &
            nonice_snowdepth(nx,ny,10), &
            faketopo(nx,ny,10),  &
            !
            landfrac(nx,ny,10),&
            icefrac(nx,ny,10),&
!            icerunoff(nx,ny,25), &
!            liqrunoff(nx,ny,25),  &
            icehflux(nx,ny,10) &
           )
!-seg

  ! Set array contents to zero
  temp     = 0.0
  precip   = 0.0
  lats     = 0.0
  lons     = 0.0
  orog     = 0.0
  albedo   = 0.0
  orog_out = 0.0

  ! Generate the UM grid.    
  write(6,*)'GLIM_HOTSTART: Loading UM grid into Glimmer..'

  do j=1,ny
     lats(j)=lat_0-(j-1)*lat_dy
  end do
  write(6,*)'GLIM_HOTSTART: Climate lats:'
  write(6,'(es12.3)') lats

  do i=1,nx
     lons(i)=lon_0+(i-1)*lon_dx
  end do
  write(6,*)'GLIM_HOTSTART: Climate lons:'
  write(6,'(es12.3)') lons

  write(6,*)'GLIM_HOTSTART: ..End loading UM grid'     
  write(6,*)'GLIM_HOTSTART: Initialising Glimmer before entering'//&
       ' the UM main time loop'

  ! Initialise the ice model
  write(6,*)'GLIM_HOTSTART: The paramfile is ',paramfile
  write(6,*)'GLIM_HOTSTART: paramfile var length:',fname_length

  !******************************************************************
  ! The forcing_time_step is how often we supply temp/precip to     *
  ! Glint, which in this coupling scheme is hard-coded to monthly:  *
  ! RMG: Monthly?  Are we sure?
  !        ice_control -> glim_control -> glint                     *
  ! The first call to glint will be at the end of the first day, ie *
  ! start_time =                                                    *
  !     the current UM date (in hours) plus one forcing-time-step.  *
  !******************************************************************

  !!forcing_time_step = 720 !ojhb - monthly submissions, not 24!
  forcing_time_step = 8640 !rssmith

  ICE_COUPLING_MONTH=12
  now_time=( 360*iiyear + 30*(iimonth-1) + iiday ) *24
  couple_this_year=(360*iiyear + 30*(ICE_COUPLING_MONTH-1) + 1) *24
  couple_next_year=(360*(iiyear+1) + 30*(ICE_COUPLING_MONTH-1) +1) *24
  if (now_time .ge. couple_this_year) then
    start_time=couple_next_year
  else
    start_time = couple_this_year
  endif

  write(6,*)'GLIM_HOTSTART: start_time supplied to initialise_'//&
       'glint is:',start_time
  write(6,*)'GLIM_HOTSTART: iiyear:',iiyear
  write(6,*)'GLIM_HOTSTART: iimonth:',iimonth
  write(6,*)'GLIM_HOTSTART: iiday:',iiday

  call initialise_glint_gcm(ice_sheet,             &
       lats,                  &
       lons,                  &
       forcing_time_step,     &
       (/paramfile/),         &
       daysinyear=360,        &
       start_time=start_time, &
       gcm_debug=.TRUE.,      &
!+seg 04/19 25->10 elevations and mask definition
       glc_nec=10,            &
       glc_nmsk=3,            &
!-seg
       ice_volume=ice_vol,    &
       um_time=int(iiyear),   &
!+seg 04/19 mask arrays
       gmask=mask,            &
       gmask3D=mask3D         &
!-seg
       )

       write(6,*)"glim_hst: ice_vol", ice_vol


  write(6,*)'GLIM_HOTSTART: End of the Glimmer initialisation'

end subroutine glim_hotstart
*DECK GLIM_FAMTYP
!this is to share info with the routine in the time loop
MODULE glim_famtype

  use glimmer_global
  use glint_main

  implicit none 

  integer, parameter :: um_int_kind = selected_int_kind(16)
  ! integer kind for consistency with UM (needed for calling 
  ! fort_get_env, an external c routine).

  type(glint_params),save :: ice_sheet    ! the derived type variable

  ! Arrays which hold the global fields used as input to Glimmer ------!

  real(rk),dimension(:,:),allocatable :: temp     ! Temperature     degC
  real(rk),dimension(:,:),allocatable :: precip   ! Precipitation   mm/s
  real(rk),dimension(:,:),allocatable :: zonwind  ! Zonal wind       m/s
  real(rk),dimension(:,:),allocatable :: merwind  ! Meridional wind  m/s
  real(rk),dimension(:,:),allocatable :: orog     ! Orography          m


  ! Arrays which hold information about the ice model instances -------!

  real(rk),dimension(:,:),allocatable :: coverage ! Coverage map for
                                                  !  normal global grid
  real(rk),dimension(:,:),allocatable :: cov_orog ! Coverage map for
                                                  !  orography grid

  
  ! Arrays which hold output from the model ---------------------------!
  ! These are all on the normal global grid, except for the orography
  
  real(rk),dimension(:,:),allocatable :: albedo   ! Fractional albedo
  real(rk),dimension(:,:),allocatable :: orog_out ! Output orography (m)
  real(rk),dimension(:,:),allocatable :: ice_frac ! Ice coverage fractio
  real(rk),dimension(:,:),allocatable :: fw_out   ! Freshwater output fl
  real(rk),dimension(:,:),allocatable :: fw_in    ! Freshwater input flu

  real(rk),dimension(:,:,:),allocatable :: ice_smb
  real(rk),dimension(:,:,:),allocatable :: ice_areas
  real(rk),dimension(:,:,:),allocatable :: icestemp
  real(rk),dimension(:,:,:),allocatable :: nonice_snowdepth
  real(rk),dimension(:,:,:),allocatable :: faketopo

  real(rk),dimension(:,:,:),allocatable :: landfrac
  real(rk),dimension(:,:,:),allocatable :: icefrac
  real(rk),dimension(:,:,:),allocatable :: icerunoff
  real(rk),dimension(:,:,:),allocatable :: liqrunoff
  real(rk),dimension(:,:,:),allocatable :: icehflux



  ! Arrays which hold information about the global grid ---------------!

  real(rk),dimension(:),  allocatable :: lats      ! Latitudes of normal
                                                   !  global gridpoints
  real(rk),dimension(:),  allocatable :: lons      ! Longitudes of norma
                                                   !  global gridpoints

  
  ! Scalars which hold information about the global grid --------------!

  real(rk) :: lat_dy,lon_dx,lat_0,lon_0

  
  ! Scalar model outputs ----------------------------------------------!

  real(rk) :: twin     ! Timestep-integrated input water flux (kg)
  real(rk) :: twout    ! Timestep-integrated output water flux (kg)
  real(rk) :: ice_vol  ! Total ice volume (m^3)


  ! Hotstarting/restarting Glint --------------------------------------!

  character(80) :: ghotfname
  character(2)  :: cday,cmonth
  character(4)  :: cyear

  character(8)  :: filedate
  character(5)  :: filerunid
  namelist/glimmerhs/filedate,filerunid

  integer, parameter :: HOTSTART = 1
  integer, parameter :: RESTART  = 2

  ! start_mode = HOTSTART means start GLIMMER from config file(s) and 
  ! hotstart file(s).
  ! start_mode = RESTART means start from restart file(s), akin to UM 
  ! start dumps.
  !
  integer,save  :: start_mode

  ! Other variables ---------------------------------------------------!

  integer       :: UM_time     ! Current time (hrs)//used to be real(rk)
  integer       :: start_time  ! time model glimmer is (re)started.
  logical       :: out         ! Outputs set flag  

END MODULE glim_famtype

*DECK GLIM_OROGRAPH
SUBROUTINE orography(  ISM_grid_x,      & ! IN
                       ISM_grid_y,      & ! IN
                       UM_lat_inc,      & ! IN UM latitude increment
                       UM_lon_inc,      & ! IN UM longitude increment
                       UM_lat_start,    & ! IN UM reference latitude
                       UM_lon_start,    & ! IN UM reference longitude
                       UM_grid_rows,    & ! IN nos of rows in UM grid
                       UM_row_length,   & ! IN nos of cols in UM grid
                       ISM_surface,     & ! IN surface topography from ISM
                       ISM_thickness,   & ! IN ice thickness from ISM
                       proj_details,    & ! IN
                       ISM_grid_details,& ! IN
                       UM_orography,    & ! I/O orographic mean on UM gr
                       UM_orog_SD,      & ! I/O SD of source on UM grid
                       UM_orog_AS,      & ! I/O orographic siluette
                       UM_orog_H,       & ! I/O orog peak-to-peak
                       UM_sigma_XX,     & ! I/O sigma_xx on UM grid
                       UM_sigma_XY,     & ! I/O sigma_xy on UM grid
                       UM_sigma_YY,     & ! I/O sigma_yy on UM grid
                       UM_icemask       ) ! I/O ice mask on climate grid

use glimmer_coordinates ! for type(coordsystem_type)
use glimmer_map_types   ! for type(glinmap_proj)
use glimmer_map_trans   ! for the loncorrect routine

use glim_famtype

use glimmer_physcon, only : rhoi

implicit none


!--- arguments ---------!
integer, intent(IN) :: ISM_grid_x,                               &
                       ISM_grid_y,                               &
                       UM_grid_rows,                             &
                       UM_row_length

real, intent(IN) ::    UM_lat_inc,                               &
                       UM_lon_inc,                               &
                       UM_lat_start,                             &
                       UM_lon_start,                             &
                       ISM_surface(ISM_grid_x,ISM_grid_y),       &
                       ISM_thickness(ISM_grid_x,ISM_grid_y)

type(glimmap_proj),    intent(IN) :: proj_details  
type(coordsystem_type),intent(IN) :: ISM_grid_details

real, dimension(um_row_length,um_grid_rows),intent(INOUT) ::     &
                       UM_orography,                             &
                       UM_orog_SD,                               &
                       UM_sigma_XX,                              &
                       UM_sigma_XY,                              &
                       UM_sigma_YY,                              &
                       UM_orog_AS,                               &
                       UM_orog_H,                                &
                       UM_icemask

!--- parameters ---------!
real, parameter ::   hires_lon_inc = 0.2d0,   &
                     hires_lat_inc = 0.2,     &
                     land_threshold = 0.5d0  ! not used
    ! fraction of hires-land cover for update to be applied to GCM cell
real, parameter :: maskd=0.2 ! snow albedo masking depth, as in FTSA

!--- local variables ----!
integer           :: error, ipoint,    & !
                     direction,        & !
                     points_targ,      & !
                     hires_nx,         & !
                     hires_ny,         & !
                     climate_nx,       & !
                     climate_ny,       & !
                     i, j, k, ii, jj,  & !
                     xx, yy,           & !
                     points_srce,      & !
                     st_row,st_col,    & !
                     extra_cells         !

logical           :: LGRAD,            & ! T if gradient fields required
                     REMOVE_MEAN_GRAD, & ! T to remove mean gradient  
                     LFILT_ISO,        & ! T for isotropic resolution on
                                         !  the sphere - implementation
                                         !   would reduce impact of ISM 
                                         !    topo.
                     GLOBAL              !

real  ::             hires_north_limit,   hires_east_limit,      &
                     hires_south_limit,   hires_west_limit,      &
                     ISM_north_limit,     ISM_east_limit,        &
                     ISM_south_limit,     ISM_west_limit,        &
                     climate_north_limit, climate_east_limit,    &
                     climate_south_limit, climate_west_limit,    &
                     maxcells,         & ! nos of highres cells in a
                                         ! climate grid box
                     cell_area           ! area of a climate grid cell

!real ::              ISM_surface(ISM_grid_x,ISM_grid_y)

real, allocatable :: model_lats_y(:,:), model_lons_x(:,:)

real, allocatable :: hires_lons_x(:),     & ! \ coordinates of cell-
                     hires_lats_y(:),     & ! / centres on hires grid.
                     hires_surface(:,:),  & ! mapped surface height
                     hires_ice_mask(:,:), & ! ice mask on hires grid
                     ice_counts(:,:),     & ! ice points count
                     ice_mask(:,:),       & ! ice mask on climate grid
                     OROG_COUNTS(:,:),    & ! land points count
                     OROG_MN(:,:),        & ! Mean surface orography
                     OROG_SD(:,:),        & ! standard deviation 
                     OROG_AS(:,:),        & ! A/S            
                     OROG_H (:,:),        & ! h              
                     SIGMA_XX(:,:),       & ! sigma_xx         
                     SIGMA_XY(:,:),       & ! sigma_xy      
                     SIGMA_YY(:,:)          ! sigma_yy  
integer, allocatable ::   ISM_ice_mask(:,:)

! JMG 14.3.14, instead of hard-coding where used
real,parameter :: RMDI = -1d0*(32768.)**2    ! same as in box_sum

!real,dimension(:,:),allocatable            :: main_usrf_temp,junk

integer inDim(2),inDim3(3)
integer ierror


!-----------------------End of delarations ----------------------------!
!
! This routine:
!  1. finds the lat/lon limits of the ISM domain.
!  2. finds what UM cells cover these limits, and uses this information
!     to create a 'target' grid, which we aim to map the ISM data to.
!  3. establishes a high-resolution lat/lon grid (the 'source')
!     that is slightly larger than the 'target' grid
!  4. use nearest-neighbour to map from the ISM cartesian to the source
!     grid.
!  5. use get_orog_mean_sd & friends to map from source to target grids.
!     get_orog_mean_sd & friends will ignore any GCM-cells that contain
!     mdi, so we apply mdi to all hires-pts outside of the ISM domain.
!     This avoids handing back spurious data to the atmosphere.
!     It is important your ISM domain is designed such that there are full
!     GCM cells covering the region of interest.
!  6. Build an ice-mask based on looking where the proportion of high-
!     resolution points in each GCM-cell that are ice-covered exceeds
!     a threshold. (present zero, ie Any ice implies ice-covered).
!  7. Put the information in the correct place in the temporary,
!     global-sized arrays.
!
!----------------------------------------------------------------------!

!if (main_task) then

write(6,*)'in glim_ororgparhy'


allocate (ISM_ice_mask(ISM_grid_x,ISM_grid_y))
where(ISM_thickness>maskd/rhoi)
  ISM_ice_mask=1
elsewhere
  ISM_ice_mask=0
endwhere


LGRAD             = .TRUE.
REMOVE_MEAN_GRAD  = .FALSE.
LFILT_ISO         = .FALSE. !OJHB - Filtiso requires a Global input.
GLOBAL            = .FALSE.

allocate( model_lats_y(ISM_grid_x,ISM_grid_y), &
          model_lons_x(ISM_grid_x,ISM_grid_y)  )

UM_icemask        = 0.0   !these
UM_orography      = 0.0   !are
UM_orog_SD        = 0.0   !updated externally
UM_sigma_XX       = 0.0   !one
UM_sigma_XY       = 0.0   !instance
UM_sigma_YY       = 0.0   !at a 
UM_orog_AS        = 0.0   !time
UM_orog_H         = 0.0   !so don't really get wiped!!
st_row            = 0 
st_col            = 0
model_lats_y(:,:) = 0.0
model_lons_x(:,:) = 0.0


!----------------------------------------------------------------------!
! 1. Translate x-y coordinates to lat-long                             !
!----------------------------------------------------------------------!
write(6,*)'OROGRAPHY: calling translate_coordinates'
ipoint = 0
direction = -1
do j = 1, ISM_grid_y
   do i = 1, ISM_grid_x
      ipoint = ipoint+1
      call translate_coordinates( i,                    &  !IN
                                  j,                    &  !IN
                                  model_lats_y(i,j),    &  !OUT
                                  model_lons_x(i,j),    &  !OUT
                                  direction,            &  !IN
                                  proj_details,         &  !IN
                                  ISM_grid_details      )  !IN
   enddo
enddo
write(6,*)'OROGRAPHY: coordinates translated'

!----------------------------------------------------------------------!
! Find e/w/n/s limits for Lat/Lon on the ISM grid in terms of cells of !
!  size um_[lat/lon]_inc (0.2degrees at present).                      !
! We need to be careful if the ISM domain lies across the longitude    !
!  wrap of the climate model (e.g. GMT or the intnl. dateline, etc.    !
!  dependent upon climate model config)                                !
! We assume that the ISM domain does not go over the poles; then the   !
!  corners of the ISM domain are (almost) always the furthest pts e/w. !
!  (Exception: If the domain straddles the equator, which is         ) !
!  (           unlikely & would anyway be taken care of by the       ) !
!  (           extension which we apply to climate grid.             ) !
!----------------------------------------------------------------------!
write(6,*)'OROGRAPHY: Find the limits of lat/lon on the ISM grid '// &
          'in terms of cells of size: '
write(6,*)'OROGRAPHY:  ',um_lat_inc,' x ',um_lon_inc

! first guess..
ISM_west_limit  = minval( model_lons_x(1,1:ISM_grid_y) )
! check if dateline/GMT (0/360 divide) goes through that edge..
if ( abs( model_lons_x(1,ISM_grid_y) - model_lons_x(1,1) ) > 180 ) then
   ! the left boundary does cross, correct our guess..
   ISM_west_limit = max(model_lons_x(1,ISM_grid_y),model_lons_x(1,1))
end if

! and similarly for the right boundary..
ISM_east_limit  = maxval( model_lons_x(ISM_grid_x,1:ISM_grid_y) )
if ( abs( model_lons_x(1,ISM_grid_y) - model_lons_x(1,1) ) > 180 ) then
   ISM_east_limit = max(model_lons_x(ISM_grid_x,ISM_grid_y),&
                                         model_lons_x(ISM_grid_x,1) )
end if

ISM_east_limit = loncorrect(ISM_east_limit,0.0)
ISM_west_limit = loncorrect(ISM_west_limit,0.0)
ISM_north_limit = MAXVAL( model_lats_y )
ISM_south_limit = MINVAL( model_lats_y )

write(6,*)'OROGRAPHY: ISM w/e limits',ISM_west_limit,ISM_east_limit
write(6,*)'OROGRAPHY: ISM s/n limits',ISM_south_limit,ISM_north_limit


!----------------------------------------------------------------------!
! 2. Find the UM grid row/column that contains the ISM limits.  Simply !
!    check to see if the distance between the two cell centres is less !
!    than half a UM cell.                                              !
! Note: The UM starts at the north pole, ie +90, so we count downwards.!
!----------------------------------------------------------------------!
! Assume that UM_lon_start is the west limit and UM_lat_start the
! the north limit.                        !
! JMG 13.3.14: Replace i with i-1 in these loops because UM_lon/lat_start
! are the first gridpoint, not the zeroth (unlike PP headers).
climate_east_limit = 0.
climate_west_limit = 0.
! loop over an entire UM row..
do i=1,UM_row_length
  ! look for the box that is centred less than half a box from the ISM..
  ! east..
  if ( abs(UM_lon_start+(UM_lon_inc*(i-1)) - ISM_east_limit ) < &
                                             0.5*UM_lon_inc ) then
    climate_east_limit =  loncorrect( UM_lon_start+ &
      UM_lon_inc*(i-1+ 0.5) , 0.0 )
  end if

  ! west..
  if ( abs(UM_lon_start+(UM_lon_inc*(i-1)) - ISM_west_limit ) < &
                                             0.5*UM_lon_inc ) then
    climate_west_limit = loncorrect( UM_lon_start + &
      UM_lon_inc*(i-1- 0.5), 0.0 )
  end if
end do


! loop over an entire UM column..
do i=1,um_grid_rows
  ! look for the box that is centred less than half a box from the ISM..
  ! north..
  ! look for the box that is centred less than half a box from the ISM..
  ! north..
  if ( abs(UM_lat_start-(UM_lat_inc*(i-1)) - ISM_north_limit ) < &
                                             0.5*UM_lat_inc ) then
    climate_north_limit = UM_lat_start-UM_lat_inc*(i-1-0.5)
  end if
  ! south..
  if ( abs(um_lat_start-(UM_lat_inc*(i-1)) - ISM_south_limit ) < &
                                             0.5*abs(um_lat_inc) ) then
    climate_south_limit = UM_lat_start-UM_lat_inc*(i-1+0.5)
  end if
end do

write(6,*)'OROGRAPHY: The GCM boxes that cover the ISM domain have'//&
          ' these limits..'
write(6,*)'OROGRAPHY: w/e:',climate_west_limit,climate_east_limit
write(6,*)'OROGRAPHY: s/n:',climate_south_limit,climate_north_limit

!----------------------------------------------------------------------!
! How many GCM-boxes n/s and e/w does our climate domain need to be?   !
!----------------------------------------------------------------------!
if (climate_east_limit < climate_west_limit) then
  ! Take care if we cross GMT..
  extra_cells = int(360d0/UM_lon_inc)
else
  extra_cells = 0
end if
climate_nx = int(  (climate_east_limit-climate_west_limit) / & 
                                              UM_lon_inc + extra_cells )
if ( climate_north_limit/climate_south_limit < 0 ) then
  ! take care if domain crosses equator..
  extra_cells = 1
else
  extra_cells = 0
end if
climate_ny = int(  (climate_north_limit-climate_south_limit) / & 
                                              UM_lat_inc + extra_cells )

write(6,*)'OROGRAPHY: The number of climate-boxes that fill this domain'
write(6,*)'OROGRAPHY: are (ew,ns):',climate_nx,climate_ny

points_targ = climate_ny * climate_nx

!----------------------------------------------------------------------!
! Dimension the output grids the same size as the 'climate' grid       !
!----------------------------------------------------------------------!
allocate( ice_mask(climate_nx,climate_ny),    &
          ice_counts(climate_nx,climate_ny),  &
          OROG_COUNTS(climate_nx,climate_ny), &
          OROG_MN(climate_nx,climate_ny),     &
          OROG_SD(climate_nx,climate_ny),     &
          OROG_AS(climate_nx,climate_ny),     &
          OROG_H(climate_nx,climate_ny),      &
          SIGMA_XX(climate_nx,climate_ny),    &
          SIGMA_XY(climate_nx,climate_ny),    &
          SIGMA_YY(climate_nx,climate_ny)     )

ice_mask(:,:)    = 0.0
ice_counts(:,:)  = 0.0
OROG_COUNTS(:,:) = 0.0
OROG_MN(:,:)     = 0.0
OROG_SD(:,:)     = 0.0
OROG_AS(:,:)     = 0.0
OROG_H(:,:)      = 0.0
SIGMA_XX(:,:)    = 0.0
SIGMA_XY(:,:)    = 0.0
SIGMA_YY (:,:)   = 0.0

!----------------------------------------------------------------------!
! 3. To create the source grid (aka high-res), we look for the western !
!    and northern edges and number of points needed to at-least span   !
!    the target (aka climate) grid.                                    !
! The hi-res grid has points at multiples of the grid spacing.
!----------------------------------------------------------------------!
! JMG 13.3.14. Change ceiling to floor and hires_lat_inc to hires_lon_inc
! in the following statement.
hires_west_limit = loncorrect(hires_lon_inc * &
                          floor(climate_west_limit / hires_lon_inc)  &
                                 - 0.5*abs(hires_lon_inc) , 0.0 )
hires_north_limit = abs(hires_lat_inc) * &
                ceiling(abs(climate_north_limit) / abs(hires_lat_inc)) &
                                 + 0.5*abs(hires_lat_inc)

! Take care in case we cross GMT..
if (climate_east_limit < hires_west_limit) then
  ! Crossing GMT...
  hires_nx = ceiling( ( climate_east_limit +360. - hires_west_limit ) &
                                     / hires_lon_inc ) 
else
  ! No extra cells
  hires_nx = ceiling( ( climate_east_limit - hires_west_limit ) &
                                     / hires_lon_inc )
end if
hires_east_limit = loncorrect( hires_west_limit  + &
                                 hires_nx * hires_lon_inc , 0.0 )

if ( climate_north_limit/climate_south_limit < 0 ) then
  hires_ny = ceiling( (hires_north_limit - climate_south_limit ) &
                              / hires_lat_inc ) + 1
else
  hires_ny = ceiling( (hires_north_limit - climate_south_limit ) &
                              / hires_lat_inc )
end if
hires_south_limit = hires_north_limit - hires_ny*hires_lat_inc

write(6,*)'OROGRAPHY: The high-res boxes that cover those GCM-cells'//&
                 ' have these limits:'
write(6,*)'OROGRAPHY: w/e: ',hires_west_limit,hires_east_limit
write(6,*)'OROGRAPHY: s/n: ',hires_south_limit,hires_north_limit
write(6,*)'OROGRAPHY: The number of hires boxes that fill this space'
write(6,*)'OROGRAPHY: are:',hires_nx,hires_ny

points_srce    = hires_nx * hires_ny

! Now allocate the arrays for the 'source' grid..
allocate( hires_surface(hires_nx,hires_ny), &
          hires_ice_mask(hires_nx,hires_ny),     &
          hires_lons_x(hires_nx),           &
          hires_lats_y(hires_ny)            )
hires_surface(:,:)  = RMDI
hires_ice_mask(:,:) = 0.0
hires_lons_x(:)     = 0.0
hires_lats_y(:)     = 0.0

! Fill the hires lons/lats arrays with cell-centre (!) values:
do i = 1, hires_nx
   hires_lons_x(i) = loncorrect( hires_west_limit +  0.5*hires_lon_inc &
                                       + (i-1)*hires_lon_inc, 0.0 )
end do
do i = 1, hires_ny
   hires_lats_y(i) = hires_south_limit + 0.5*hires_lat_inc &
                                            + (i-1)*hires_lat_inc
end do

!----------------------------------------------------------------------!
! 4. For each pt on the UM (high-res) grid: find nearest pt in ISM     !
!    grid coordinates, check that that pt exists (ie lies within the   !
!    ISM grid) and if so set the hi-res topography and ice mask to it  !
!----------------------------------------------------------------------!
do i=1,hires_nx
   do j=1,hires_ny
      direction = 1  ! => from lat-lon to grid point

      call translate_coordinates(  xx,              &  !OUT
                                   yy,              &  !OUT
                                   hires_lats_y(j), &  !IN
                                   hires_lons_x(i), &  !IN
                                   direction,       &  !IN
                                   proj_details,    &  !IN
                                   ISM_grid_details )  !IN

      ! check the index returned lies within our ISM grid
      if ( (xx >= 1) .and. (xx <= ISM_grid_x) .and. &
            (yy >= 1) .and. (yy <= ISM_grid_y) ) then
         hires_surface(i,j) = ISM_surface(xx,yy)
         hires_ice_mask(i,j) = ISM_ice_mask(xx,yy)
      end if
   end do
end do

! get_orog_mean_sd wants grid-cell centres..
climate_west_limit = loncorrect(climate_west_limit+0.5*UM_lon_inc, 0.0)
! JMG 13.3.14
climate_north_limit=climate_north_limit-0.5*UM_lat_inc

!----------------------------------------------------------------------!
! 5.0  Use high res UM grid model (source) to generate means and       !
!      standard devs for climate resolution grid (target).             !
!----------------------------------------------------------------------!
write(6,*)'OJHB: Inputs to get_orog_mean_sd:'
write(6,*)hires_lats_y(hires_ny),climate_north_limit
write(6,*)hires_lons_x(1),climate_west_limit

call get_orog_mean_sd (    &
          hires_surface,   & ! IN   surface  (source data)
          hires_ice_mask,  & ! IN   ice mask (source data)
          OROG_MN,         & ! OUT  orographic mean on climate grid
          OROG_COUNTS,     & ! OUT  number of land points in boxes
          ice_counts,      & ! OUT  number of ice points in grid cells
          OROG_SD,         & ! OUT  SD of source on climate grid
          SIGMA_XX,        & ! OUT  sigma_xx on climate grid
          SIGMA_XY,        & ! OUT  sigma_xy on climate grid
          SIGMA_YY,        & ! OUT  sigma_yy on climate grid
          OROG_AS,         & ! OUT  orographic X-section on climate grid
          OROG_H,          & ! OUT  orographic peak-peak on climate grid
          points_srce,     & ! IN   number of points on source grid
          hires_ny,        & ! IN   number of rows on source grid
          hires_nx,        & ! IN   number of columns on source grid
          points_targ,     & ! IN   number of points on target grid
          climate_ny,      & ! IN   number of rows on target grid
          climate_nx,      & ! IN   number of columns on target grid
          hires_lat_inc,   & ! IN   lat increment on source grid
          hires_lon_inc,   & ! IN   lon increment on source grid
          UM_lat_inc,      & ! IN   lat increment on target grid
          UM_lon_inc,      & ! IN   lon increment on target grid
! JMG 13.3.14: exchange climate and hires pairs of lines below; they were
! the wrong way round
      climate_north_limit, & ! IN   lat origin on target grid (N)
       climate_west_limit, & ! IN   lon origin on target grid (W)
   hires_lats_y(hires_ny), & ! IN lat origin on source grid (N)
          hires_lons_x(1), & ! IN   lon origin on source grid (W)
          LFILT_ISO,       & ! IN   want isotropic filter?
          LGRAD,           & ! IN   want gradient fields?
         REMOVE_MEAN_GRAD, & ! IN   want gradient fields?
          GLOBAL           ) ! IN   doing global calcs?

! Put climate_west_limit back the way it should be..
climate_west_limit = loncorrect(climate_west_limit-0.5*UM_lon_inc,0.0)


!----------------------------------------------------------------------!
! 6. Calculate the ice mask from the counts of ice covered grid cells. !
!    Take ANY ice (in a land cell) to indicate an ice point.           !
!----------------------------------------------------------------------!
!where( (OROG_MN .gt. 0.) .AND. (ice_counts .gt. 0.))
!   ice_mask = 1.0
!elsewhere
!   ice_mask = 0.0
!endwhere
ice_mask=ice_counts
! JMG. ice_counts<0 are missing data
where (orog_mn.eq.0.or.ice_counts.lt.0) ice_mask=0.0
write(6,*)'OROGRAPHY: ice mask set.'


!----------------------------------------------------------------------!
! 7. Find the UM grid cell location in the global grid for the         !
!    retrieved orography parameters and put them into the empty dummy  !
!    arrays to be passed back to the control code.                     !
!----------------------------------------------------------------------!
write(6,*)'OROGRAPHY: writing to Dummy Global-sized arrays.'

st_row = NINT( (UM_lat_start - climate_north_limit) / UM_lat_inc ) + 1
st_col = NINT( (climate_west_limit - UM_lon_start ) / UM_lon_inc ) + 1
jj=climate_ny
do j=st_row,st_row+climate_ny-1
   ii=1
   do i=st_col,st_col+climate_nx-1
      if (i > UM_row_length) then
         k=i-UM_row_length
      else
         k=i
      end if
      UM_icemask(k,j)   = ice_mask(ii,jj)
      UM_orography(k,j) = OROG_MN(ii,jj)
      UM_orog_SD(k,j)   = OROG_SD(ii,jj)
      UM_sigma_XX(k,j)  = SIGMA_XX(ii,jj)
      UM_sigma_XY(k,j)  = SIGMA_XY(ii,jj)
      UM_sigma_YY(k,j)  = SIGMA_YY(ii,jj)
      UM_orog_AS(k,j)   = OROG_AS(ii,jj)
      UM_orog_H(k,j)    = OROG_H(ii,jj)
      ii=ii+1
   end do
   jj=jj-1
end do

!----------------------------------------!
! Deallocate reassigned arrays
write(6,*)'OROGRAPHY: deallocating arrays'
DEALLOCATE (  ice_mask,       &
              ice_counts,     &
              OROG_COUNTS,    &
              OROG_MN,        &
              OROG_SD,        &
              OROG_AS,        &
              OROG_H,         &
              SIGMA_XX,       &
              SIGMA_XY,       &
              SIGMA_YY,       &
              hires_surface,  &
              hires_ice_mask, &
              hires_lons_x,   &
              hires_lats_y    )

write(6,*)'OROGRAPHY: finished.'
write(6,*)''

END SUBROUTINE orography
*/
*/**********************************************************************
*/
SUBROUTINE get_orog_mean_sd( DATA_SRCE,            & ! In
                             DATA_ICE,             & ! In
                             OROG_MN,              & ! Out
                             COUNTS,               & ! Out
                             COUNTS_ICE,           & ! Out
                             OROG_SD,              & ! Out
                             SIGMA_XX,             & ! Out
                             SIGMA_XY,             & ! Out
                             SIGMA_YY,             & ! Out
                             OROG_AS,              & ! Out
                             OROG_H,               & ! Out
                             POINTS_SRCE,          & ! In
                             POINTS_PHI_SRCE,      & ! In
                             POINTS_LAMBDA_SRCE,   & ! In 
                             POINTS_TARG,          & ! In
                             POINTS_PHI_TARG,      & ! In
                             POINTS_LAMBDA_TARG,   & ! In 
                             DELTA_PHI_SRCE,       & ! In
                             DELTA_LAMBDA_SRCE,    & ! In
                             DELTA_PHI_TARG,       & ! In
                             DELTA_LAMBDA_TARG,    & ! In
                             START_PHI_TARG,       & ! In
                             START_LAMBDA_TARG,    & ! In
                             START_PHI_SRCE,       & ! In
                             START_LAMBDA_SRCE,    & ! In
                             LFILT_ISO,            & ! In
                             LGRAD,                & ! In
                             REMOVE_MEAN_GRAD,     & ! In
                             GLOBAL                ) ! In
!----------------------------------------------------------------------!
! Routine to calculate the orographic mean and standard deviation from !
! a high resolution topographic data set to a lower resolution.        !
! Standard deviation gradients are also calculated.                    !
!                                                                      !
! This routine is based on that used in the generation of orography    !
! ancillary files, written by Clive Jones, which are normally held     !
! fixed within a model run. By coding this inline we can allow         !
! topographic change during a model run.                               !
!                                                                      !
! Options, controlled by logicals passed from the call list, include   !
! the facility to smooth the orography in various ways, presumably to  !
! prevent gravity waves.                                               !
!                                                                      !
! The subroutines called by get_orog_mean_sd are also only used in the !
! ancillary file generation code, and may need to be altered to run    !
! inline....                                                           !
!----------------------------------------------------------------------!
implicit none

!-- arguments -------------------!
integer,intent(IN):: POINTS_LAMBDA_SRCE,& ! No. of cols on source grid
                     POINTS_PHI_SRCE,   & ! No. of rows on source grid
                     POINTS_SRCE,       & ! No. of pts on source grid
                     POINTS_LAMBDA_TARG,& ! No. of cols on target grid
                     POINTS_PHI_TARG,   & ! No. of rows on target grid
                     POINTS_TARG          ! No. of points on target grid

real, intent(IN)  :: DELTA_LAMBDA_TARG,                             &
                              ! Longitude increment of target grid (deg)
                     DELTA_PHI_TARG,                                & 
                              ! Latitude increment of target grid (deg)
                     START_PHI_TARG,    & ! start latitude of centre of
                                          !  1st grid-box in target area
                     START_LAMBDA_TARG, & ! start longitude of centre of
                                          !  1st grid-box in target area
                     DELTA_PHI_SRCE,                                & 
                              ! Latitude increment of source grid
                     DELTA_LAMBDA_SRCE,                             & 
                              ! Longitude increment of source grid 
                     START_LAMBDA_SRCE,& ! start longitude of centre of
                                         !  1st grid-box in source area
                     START_PHI_SRCE,   & ! start latitude of centre of
                                         !  1st grid-box in source area
                     DATA_SRCE(POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE), &  
                                        ! IN surface source data
                     DATA_ICE(POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE)
                                        ! IN ice cover data

logical,intent(IN) :: LFILT_ISO,     & ! T for isotropic resolution on
                                       ! the sphere - application would
                                       ! reduce impact of ISM topo.
                      LGRAD,         & ! T if gradient fields reqd
                      REMOVE_MEAN_GRAD, & ! T to remove mean gradient  
                      GLOBAL           ! T if global area required


real,intent(OUT) :: OROG_MN(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),      &
                              ! mean orography    
                    OROG_SD(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),      &
                              ! standard deviation 
                    SIGMA_XX(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),     &
                              ! sigma_xx  
                    SIGMA_XY(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),     &
                              ! sigma_xy     
                    SIGMA_YY(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),     &
                              ! sigma_yy 
                    OROG_AS(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),      &
                              ! orography cross section 
                    OROG_H(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),       &
                              ! orography peak to trough heights 
                    COUNTS(POINTS_LAMBDA_TARG,POINTS_PHI_TARG),       &
                              ! count of land cover weights in grid cell
                    COUNTS_ICE(POINTS_LAMBDA_TARG,POINTS_PHI_TARG)
                              ! fractional area of ice coverage

!--- local variables ------------!
integer :: IGRID,                       & ! Grid indicators:
           IGRID_SRCE,                  & !   1=>p-grid, 2=>u-grid
           IGRID_GRAD,                  & !
           I_L(POINTS_LAMBDA_TARG+1),   & ! Index of source point to
                                          !  left of grid box
           J_T(POINTS_PHI_TARG+1),      & ! Index of 1st source pt
                                          !  below topof grid box
           I_L_G(POINTS_LAMBDA_TARG+1), & ! Index of source pt to
                                          !  left of grid box 
           J_T_G(POINTS_PHI_TARG+1),    & ! Index of 1st source pt
                                          !  below top of grid box
           POINTS_LAMBDA_GRAD             ! No. of cols on gradient grid


real :: as_coefficient, & ! Resolution-dependent ratio of silhouette
                          ! orography to orographic SD
        LONG_L(POINTS_LAMBDA_TARG +1),  & ! Left longitude of grid box
                                          !  in units of delta_long
                                          !   source
        COLAT_T(POINTS_PHI_TARG+1),     & ! Colatitude of top of grid
                                          !  box in units of delta_lat
                                          !   source
        LONG_L_G(POINTS_LAMBDA_TARG+1), & ! Left longitude of grid box
        COLAT_T_G(POINTS_PHI_TARG+1),   & ! Colatitude of top of grid
                                          !  box  
        DATA_SRCE_SD(POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE), &
                                          ! source data squared       
        AREA_BOX,                       & ! area of grid box in sq
                                          !  units of source data.
        START_LAT_GRAD,                 & !
        START_LONG_GRAD,                & !
        ICE_WEIGHTS(POINTS_LAMBDA_TARG,POINTS_PHI_TARG), & !
                                          ! weights for ice area average
        counts_sd(POINTS_LAMBDA_TARG,POINTS_PHI_TARG)   !
                              ! count of orog-sd weights in grid cell

logical :: target_mask(POINTS_LAMBDA_TARG,POINTS_PHI_TARG)
                                          ! remove target-cells with
                                          ! source-mdi in.
real,parameter :: RMDI = -1d0*(32768.)**2 ! same as in box_sum, JMG 14.3.14

                                                                        
!-- External routines called --
external :: BOX_SUM,    & ! to calculate coefficients for area weighting
            BOX_BND,    &
            INTERPBACK, & ! to remove mean gradient
            OROG_GRAD     ! to calculate orographic components.


 !---------------------------------------------------------------------!
 ! 0.0 Initialisation                                                  !
 !---------------------------------------------------------------------!

if (delta_lambda_targ.eq.7.5.and.delta_phi_targ.eq.5.0) then
  as_coefficient = 5.86778e-05 ! FAMOUS
elseif (delta_lambda_targ.eq.3.75.and.delta_phi_targ.eq.2.5) then
  as_coefficient = 8.4454e-05 ! HadCM3
else
  write(6,*) 'UM-Glimmer coupling can''t process this atmosphere grid'
  stop
endif
print *,'as_coefficient = ',as_coefficient

 COUNTS(:,:)       = 0.0
 COUNTS_ICE(:,:)   = 0.0
 counts_sd(:,:)    = 0.0
 OROG_MN(:,:)      = 0.0
 target_mask(:,:)  = .false.
 OROG_AS(:,:)      = 0.0
 OROG_H(:,:)       = 0.0
 ICE_WEIGHTS(:,:)  = 0.0
 OROG_SD(:,:)      = 0.0
 SIGMA_XX(:,:)     = 0.0
 SIGMA_XY(:,:)     = 0.0
 SIGMA_YY(:,:)     = 0.0
 AREA_BOX          = 0.0
 IGRID             = 1
 IGRID_SRCE        = 1
 IGRID_GRAD        = 1

 if (global) then
   points_lambda_grad=points_lambda_srce
 else
   points_lambda_grad=points_lambda_srce+1
 endif
! JMG 13.3.14. Corrected to derive from SRCE rather than TARG grid!
 START_LAT_GRAD  = START_PHI_SRCE   + 0.5*abs(DELTA_PHI_SRCE)
 START_LONG_GRAD = START_LAMBDA_SRCE - 0.5*DELTA_LAMBDA_SRCE

 !---------------------------------------------------------------------!
 ! 1.0 Set up arrays of indexes and boundaries of target grid boxes    !
 !     relative to a source grid for use by BOX_SUM so that area       !
 !     weighted means can be calculated.                               !
 !---------------------------------------------------------------------!
 call BOX_BND(                &
          I_L,                & ! OUT
          LONG_L,             & ! OUT
          J_T,                & ! OUT
          COLAT_T,            & ! OUT
          AREA_BOX,           & ! OUT
          POINTS_LAMBDA_TARG, & ! IN  row length target
          POINTS_PHI_TARG,    & ! IN  rows target
	  POINTS_LAMBDA_SRCE, & ! IN  row length source
          POINTS_PHI_SRCE,    & ! IN  rows source
          DELTA_LAMBDA_TARG,  & ! IN  longitude increment of target grid
          DELTA_PHI_TARG,     & ! IN  latitude increment of target grid
          START_LAMBDA_TARG,  & ! IN  start longitude of centre of 1st
                                !      grid-box in target area
          START_PHI_TARG,     & ! IN  start latitude of centre of first
                                !      grid-box in target area
	  DELTA_LAMBDA_SRCE,  & ! IN  longitude increment of source grid
	  DELTA_PHI_SRCE,     & ! IN  latitude increment of source grid
          START_LAMBDA_SRCE,  & ! IN
	  START_PHI_SRCE,     & ! IN
          IGRID,              & ! IN Grid indicator 1=p-grid,2=u-grid
	  IGRID_SRCE,         & ! IN Grid indicator 1=p-grid,2=u-grid
	  GLOBAL              )


 !---------------------------------------------------------------------!
 ! 2.0 Sum source boxes (whole & partial) contributions to grid boxes. !
 !---------------------------------------------------------------------!
 call BOX_SUM2( POINTS_LAMBDA_SRCE,   & ! IN
                POINTS_PHI_SRCE,      & ! IN
                POINTS_LAMBDA_TARG,   & ! IN
                POINTS_PHI_TARG,      & ! IN
                LONG_L,               & ! IN
                COLAT_T,              & ! IN
                I_L,                  & ! IN
                J_T,                  & ! IN
                GLOBAL,               & ! IN
                OROG_MN,              & ! OUT
                COUNTS,               & ! OUT
                DATA_SRCE             ) ! IN
 write(6,*)'GET_OROG_MEAN_SD: completed orog-mean interpolation'


 !---------------------------------------------------------------------!
 ! 2.1 Create a mask for target cells containing mdi from source cells !
 !     Replace mdi with zero in orog_mn before performing filtering.   !
 !---------------------------------------------------------------------!
 where(OROG_MN .eq. RMDI)
   target_mask = .true.
   orog_mn     = 0.0
 end where


 !---------------------------------------------------------------------!
 ! 2.2 Repeat sumation of boxes but this time for the ice mask         !
 !---------------------------------------------------------------------!
 call BOX_SUM2( POINTS_LAMBDA_SRCE,   & ! IN
                POINTS_PHI_SRCE,      & ! IN
	        POINTS_LAMBDA_TARG,   & ! IN
	        POINTS_PHI_TARG,      & ! IN
	        LONG_L,               & ! IN
	        COLAT_T,              & ! IN
	        I_L,                  & ! IN
	        J_T,                  & ! IN
	        GLOBAL,               & ! IN
                COUNTS_ICE,           & ! OUT
	        ICE_WEIGHTS,          & ! OUT
	        DATA_ICE              ) ! IN
 write(6,*)'GET_OROG_MEAN_SD: completed ice mask interpolation'


 !---------------------------------------------------------------------!
 ! 3.0 Most of the original code was not needed.  Only the subset      !
 !     that filters isotropically has been retained.                   !
 !---------------------------------------------------------------------!
 if ( LFILT_ISO ) then
   ! Averages EW polewards of 48 in order to retain
   !  approx. isotropic resolution as meridians converge
   call FILTISO( OROG_MN,                 & ! INOUT
                 POINTS_TARG,             & ! IN
		 POINTS_LAMBDA_TARG,      & ! IN
		 POINTS_PHI_TARG          ) ! IN
 end if


 !---------------------------------------------------------------------!
 ! 4.0 Construct gradient fields if dictated by logical flag LGRAD     !
 !---------------------------------------------------------------------!
 if (LGRAD) then
                                                                        
    if(REMOVE_MEAN_GRAD) then                                           
      call INTERPBACK( DATA_SRCE,           & ! INOUT
                       OROG_MN,             & !
		       START_LAMBDA_SRCE,   & !              
                       START_PHI_SRCE,      & !
		       DELTA_LAMBDA_TARG,   & !
		       DELTA_PHI_TARG,      & !
                       DELTA_LAMBDA_SRCE,   & !
		       DELTA_PHI_SRCE,      & !
		       POINTS_LAMBDA_TARG,  & !
                       POINTS_PHI_TARG,     & !
		       POINTS_LAMBDA_SRCE,  & !
		       POINTS_PHI_SRCE,     & !
                       POINTS_TARG,         & !
		       POINTS_SRCE,         & !
                       START_LAMBDA_TARG,   & !
		       START_PHI_TARG,      & !
		       GLOBAL               )        
    end if                                              
                                                     
! A new set of indices are derived for lumping the GRAD grid onto the TARG
! UM grid. The GRAD grid is for the gradients of topography on the hi-res
! grid, and is half a hires gridbox larger in each direction.
    call BOX_BND(      I_L_G,               &
                       LONG_L_G,            &
		       J_T_G,               &
		       COLAT_T_G,           &
		       AREA_BOX,            &           
                       POINTS_LAMBDA_TARG,  &
		       POINTS_PHI_TARG,     &
		       POINTS_LAMBDA_GRAD,  &        
                       POINTS_PHI_SRCE+1,   &                           
                       DELTA_LAMBDA_TARG,   &
		       DELTA_PHI_TARG,      &
		       START_LAMBDA_TARG,   &         
                       START_PHI_TARG,      &                           
                       DELTA_LAMBDA_SRCE,   &
		       DELTA_PHI_SRCE,      &
		       START_LONG_GRAD,     &
                       START_LAT_GRAD,      &
		       IGRID,               &
		       IGRID_GRAD,          &
		       GLOBAL               )  

    write(6,*)'GET_OROG_MEAN_SD: Box_Bnd output'
    write(6,*)'GET_OROG_MEAN_SD: ',DELTA_PHI_SRCE,  &
                               DELTA_LAMBDA_SRCE,AREA_BOX            
                                                                        
    call OROG_GRAD(    DATA_SRCE,           &
                       SIGMA_XX,            &
		       SIGMA_YY,            &
		       SIGMA_XY,            &
                       POINTS_LAMBDA_SRCE,  &
		       POINTS_PHI_SRCE,     &
		       POINTS_LAMBDA_GRAD,  &
                       DELTA_PHI_SRCE,      &
		       DELTA_LAMBDA_SRCE,   &
		       START_LAT_GRAD,      &
                       POINTS_LAMBDA_TARG,  &
		       POINTS_PHI_TARG,     &
		       LONG_L_G,            &
                       COLAT_T_G,           &
		       I_L_G,               &
		       J_T_G,               &
		       GLOBAL               )                           
                                                                        
 end if  ! LGRAD                                                  


 !---------------------------------------------------------------------!
 ! 5.0 Construct standard deviations from squared topographic heights  !
 !     to obtain sums of squares                                       !
 !---------------------------------------------------------------------!
 write(6,*)'GET_OROG_MEAN_SD: completed gradients calculations'
                          
 DATA_SRCE_SD(:,:)=DATA_SRCE(:,:)*DATA_SRCE(:,:)     
                                            
                                                                     
 call BOX_SUM2( POINTS_LAMBDA_SRCE,  &
                POINTS_PHI_SRCE,     &
	        POINTS_LAMBDA_TARG,  & 
                POINTS_PHI_TARG,     &
	        LONG_L,              &
	        COLAT_T,             &
	        I_L,                 &
	        J_T,                 &
	        GLOBAL,              &
                OROG_SD,             &
                counts_sd,           &
	        DATA_SRCE_SD)

 if ( sum(counts_sd) /= sum(counts) ) &
   write(6,*)'GET_OROG_MEAN_SD: Mdi differ in orog_mn & orog_sd'
 write(6,*)'GET_OROG_MEAN_SD: completed stdev calculations'

 !---------------------------------------------------------------------!
 ! 5.1 Now form grid box means                                         !
 !---------------------------------------------------------------------!
 if (LGRAD.AND.REMOVE_MEAN_GRAD) then
    OROG_SD(:,:) =   OROG_SD(:,:) * AREA_BOX/(AREA_BOX -1.0)
 else
    OROG_SD(:,:) = ( OROG_SD(:,:) - OROG_MN(:,:)*OROG_MN(:,:) )   &
                      * AREA_BOX/(AREA_BOX -1.0)
 end if                         
 where( OROG_SD(:,:) .lt. 0.0 ) OROG_SD(:,:) = 0.0
 OROG_SD(:,:) = SQRT( OROG_SD(:,:) )
 

 !---------------------------------------------------------------------!
 ! 5.2 Translate SD into other parameters. Climate model uses peak to  !
 !     trough high = SD heights and aerodynamic cross-section as a     !
 !     a linear factor of SD.                                          !
 !---------------------------------------------------------------------!
 OROG_H(:,:)  = OROG_SD(:,:)
 OROG_AS(:,:) = OROG_SD(:,:)*as_coefficient


 !---------------------------------------------------------------------!
 ! 6.0 Set all target arrays to zero where mdi in source               !
 !---------------------------------------------------------------------!
 where(target_mask)
   orog_sd  = 0.0
   sigma_xx = 0.0
   sigma_xy = 0.0
   sigma_yy = 0.0
   orog_h   = 0.0
   orog_as  = 0.0
 end where


 write(6,*)'GET_OROG_MEAN_SD: finished.'
 write(6,*)''                  

END SUBROUTINE get_orog_mean_sd
*/
*/**********************************************************************
*/
SUBROUTINE translate_coordinates(i,          & ! grid location in x/lon 
                                 j,          & ! grid location in y/lat 
				 phi,        & ! latitude  (degrees)
				 lambda,     & ! longitude (degrees)
				 direction,  & ! transform direction
                                 proj,       & ! projection pointer   
                                 grid        ) ! grid detail pointer

use glimmer_map_types   ! for type(glinmap_proj)
use glimmer_coordinates ! for type(coordsystem_type)
use glimmer_map_trans   ! for the xy<-->ll routines

implicit none

!-- arguments --------!
integer,intent(INOUT) :: i, j      ! grid point
real,   intent(INOUT) :: phi       ! latitude = +ve in north hemisphere
real,   intent(INOUT) :: lambda    ! longitude - negative to the west
integer,intent(IN)    :: direction !  1 = phi-lambda to x-y
                                   ! -1 = x-y to phi-lambda
type(glimmap_proj),    intent(IN) :: proj ! projection details
type(coordsystem_type),intent(IN) :: grid ! grid definition

!-- local variables --!
real                  :: x, y      ! ice grid space (m)


if (direction == -1) then
  call glimmap_xy_to_ll(lambda,phi,real(i),real(j),proj,grid)
else if (direction == 1) then
  call glimmap_ll_to_xy(lambda,phi,x,y,proj,grid)
  i = nint(x)
  j = nint(y)
else
   write(6,*)'TRANSLATE_COORDINATES: bad option for direction: ',&
             direction
   STOP
end if

END SUBROUTINE translate_coordinates

*/
*/******************************************************************
*/

*DECK OROGGRAD
! ******************************COPYRIGHT****************************** 
! (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved. 
!                                                                       
! Use, duplication or disclosure of this code is subject to the         
! restrictions as set forth in the contract.                            
!                                                                       
!                Meteorological Office                                  
!                London Road                                            
!                BRACKNELL                                              
!                Berkshire UK                                           
!                RG12 2SZ                                               
!                                                                       
! If no contract has been raised with this copy of the code, the use,   
! duplication or disclosure of it is strictly prohibited.  Permission   
! to do so must first be obtained in writing from the Head of Numerical 
! Modelling at the above address.                                       
! ******************************COPYRIGHT****************************** 
!    SUBROUTINE OROG_GRAD                                               
!                                                                       
!    Routine calculates SIGMA_XX, SIGMA_YY, SIGMA_XY on target          
!    grid. This involves using TSIGMA as a temporary store for          
!    XX,YY,XY on a grid identical to input data grid. Using this        
!    square grid loses some resolution and may cause some               
!    with BOX_SUM.                                                      
!                                                                       
!    NOT SUITABLE FOR SINGLE COLUMN USE                                 
!                                                                       
!    SUITABLE FOR ROTATED GRIDS                                         
!                                                                       
!    ORIGINAL VERSION FOR CRAY Y-MP/IBM                                 
!    WRITTEN 3/93 J.R.MITCHELL                                          
!                                                                       
!    CODE REVIEWED BY C.WILSON ??/??/??                                 
!                                                                       
!    *   VERSION NO. 1 DATED 3/93                                       
!    *          COSMOS DSN MEH1.MJUM.GWAVE(GRADFDIF)                    
!    *   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 3  
!    *   VERSION 1, DATED 3/93                                          
!                                                                       
!    *   SYSTEM TASK:  P65                                              
!                                                                       
!    PURPOSE:  To calculate directional gradient orography approximated 
!              onto an indentical grid to input interpolated orography  
!              using finite differencing between points on either side  
!              of the calculation point (GRAD=T)                        
!                                                                       
!    *   DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION 'P65'              
!    *                    BY C.WILSON VERSION 1 DATED 27/02/90          
!    *                     J.MITCHELL VERSION 2 DATED 3/93              
!  JMG 14.3.14: Modifications to take notice of missing data
!                                                                       
!                                                                       
!  END-------------------------------------------------------------     

SUBROUTINE OROG_GRAD( TOPG,            & ! IN
                      SIGMA_XX,        & ! OUT
                      SIGMA_YY,        & ! OUT
                      SIGMA_XY,        & ! OUT
                      TOPG_ROW_LENGTH, & ! IN
                      TOPG_ROWS,       & ! IN
                      GRAD_ROW_LENGTH, & ! IN          
                      TOPG_DELTA_LAT,  & ! IN
                      TOPG_DELTA_LONG, & ! IN
                      START_LAT_GRAD,  & ! IN    
                      ROW_LENGTH,      & ! IN
                      ROWS,            & ! IN
                      LONG_L,          & ! IN
                      COLAT_T,         & ! IN
                      I_L,             & ! IN
                      J_T,             & ! IN
                      GLOBAL           ) ! IN
                                                                        
implicit none                                                     
                                                                        
integer,intent(IN) ::      &
        TOPG_ROWS,         & ! No of rows in TOPG and GRAD grid          
        TOPG_ROW_LENGTH,   & ! No of points per row in TOPG grid         
        GRAD_ROW_LENGTH,   & ! No of points per row in GRAD grid         
        ROWS,              & ! No of rows in Target grid                 
        ROW_LENGTH,        & ! No of points per row in Target grid       
        I_L(ROW_LENGTH+1), & ! Index of point to left of grid box        
        J_T(ROWS+1)          ! First point below top of grid box         
                                                                        
 real,intent(IN) :: &
         TOPG(TOPG_ROW_LENGTH,TOPG_ROWS),& ! Bi-linear mean
                                  ! interpolated orographic heights  
         LONG_L(ROW_LENGTH +1), & ! Left longitude of grid box in    
                                  !  units of TOPG_DELTA_LONG.        
         COLAT_T(ROWS+1),       & ! Colatitude of top of grid box    
                                  !  in units of TOPG_DELTA_LAT.      
         TOPG_DELTA_LAT,        & ! Lat. resolution of TOPG grid     
         TOPG_DELTA_LONG,       & ! Long resolution of TOPG grid     
         START_LAT_GRAD           ! Start lat of gradient grid       

real,intent(OUT) :: &
        SIGMA_XX(ROW_LENGTH,ROWS), & ! (dh/dx)**2 on target grid        
        SIGMA_YY(ROW_LENGTH,ROWS), & ! (dh/dy)**2 on target grid        
        SIGMA_XY(ROW_LENGTH,ROWS)    ! (dh/dx)*(dh/dy) on target grid   
                                                                        
logical,intent(IN) :: GLOBAL 
                 ! .true. => Cyclic wrap around poles & zero meridian
                                                                        
! Local Arrays (Dynamic allocated)                                   
real :: TWO_DELTA_X(TOPG_ROWS+1),           & ! the distance between
                                              !  points of row J (m)
        DH_DX(GRAD_ROW_LENGTH,TOPG_ROWS+1), & !
        DH_DY(GRAD_ROW_LENGTH,TOPG_ROWS+1), & !
        TSIGMA(GRAD_ROW_LENGTH,TOPG_ROWS+1)   !
                                                                        
!-- Local Variables --
real :: DELTA_X,      & ! Distance between points (m)                
        DELTA_Y,      & ! Distance between rows (m)                  
        LAT_ROW,      & ! Latitude of row                            
        TWO_DELTA_Y     ! Twice distance between rows (m)            

integer,parameter :: IMDI =-32768 ! Integer missing data indicator (mdi)

real,parameter :: A = 6371229.0,                & ! Earth radius.
                  PI=3.14159265358979323846,    & ! Pi
                  PI_OVER_180 =PI/180.0,        & ! Conversion factor
                                                  !  degrees to radians  
                  RECIP_PI_OVER_180 = 180.0/PI, & ! Conversion factor
                                                  !  radians to degrees

                  RMDI = -1d0*(32768.)**2 ! same as in box_sum, JMG 14.3.14
       
integer :: I,J,ioffset    ! Loop indices
                                                                        
  ! *************************************************************
  ! 1.0 Calculate DX(J) and DY assuming A constant and then             
  !     calculate DH_DX(I,J) and DH_DY(I,J) between points              
  ! *************************************************************

  if (global) then
    ioffset=0
  else
    ioffset=1
  endif
                                                    
  DELTA_Y = A * TOPG_DELTA_LAT * PI_OVER_180                        
  TWO_DELTA_Y = 2.0 * DELTA_Y                                    

  do J=1,TOPG_ROWS-1     

    LAT_ROW  = (START_LAT_GRAD - J*TOPG_DELTA_LAT) * PI_OVER_180    
    DELTA_X  = TOPG_DELTA_LONG * PI_OVER_180 * A * COS(LAT_ROW)     
    TWO_DELTA_X(J+1) = 2.0 * DELTA_X

    do I=1,TOPG_ROW_LENGTH-1
      if (TOPG(I,J)/=RMDI .and. TOPG(I,J+1)/=RMDI &
      .and. TOPG(I+1,J)/=RMDI .and. TOPG(I+1,J+1)/=RMDI) then
        DH_DX(I+ioffset,J+1) = ( (TOPG(I+1,J)  -TOPG(I,J)  ) +    &
                          (TOPG(I+1,J+1)-TOPG(I,J+1))  ) /        &
                          TWO_DELTA_X(J+1)       
        DH_DY(I+ioffset,J+1) = ( (TOPG(I,J)   - TOPG(I,J+1)   ) + &
                            (TOPG(I+1,J) - TOPG(I+1,J+1) ) ) /    &
                          TWO_DELTA_Y
      else
        dh_dx(I+ioffset,J+1)=RMDI
        dh_dy(I+ioffset,J+1)=RMDI
      endif
    end do                                                           
  end do 
                                                                        
  ! *******************************************************************
  ! 1.1 Wrap around zero meridian when I=TOPG_ROW_LENGTH, GLOBAL=.TRUE.  
  !  ******************************************************************
                                                                       
  if (GLOBAL) then                                                  
                                                                        
    do J=1,TOPG_ROWS-1                                              


      if (TOPG(1,J)/=RMDI .and. TOPG(TOPG_ROW_LENGTH,J)/=RMDI &
      .and. TOPG(1,J+1)/=RMDI .and. TOPG(TOPG_ROW_LENGTH,J+1)/=RMDI)  &
      then
        DH_DX(TOPG_ROW_LENGTH,J+1)=                                   &
              ( (TOPG(1,J)  -TOPG(TOPG_ROW_LENGTH,J)   ) +            &     
                (TOPG(1,J+1)-TOPG(TOPG_ROW_LENGTH,J+1) ) ) /          &
              TWO_DELTA_X(J+1)                                         
                                                                        
        DH_DY(TOPG_ROW_LENGTH,J+1) =                                  &
           ( (TOPG(TOPG_ROW_LENGTH,J) - TOPG(TOPG_ROW_LENGTH,J+1) ) + &
                (TOPG(1,J) - TOPG(1,J+1) ) ) / TWO_DELTA_Y
      else
        DH_DX(TOPG_ROW_LENGTH,J+1)=RMDI
        DH_DY(TOPG_ROW_LENGTH,J+1)=RMDI
      endif
                                                                         
    end do                                                           
                                                                        
  else                                                              
                                                                        
    do J=2,TOPG_ROWS                                                
                                                                        
      ! Western Boundary                                              
      DH_DX(1,J) = DH_DX(2,J)                                       
      DH_DY(1,J) = DH_DY(2,J)                                       
                                                                        
      ! Eastern Boundary                                              
      DH_DX(GRAD_ROW_LENGTH,J) = DH_DX(GRAD_ROW_LENGTH-1,J)         
      DH_DY(GRAD_ROW_LENGTH,J) = DH_DY(GRAD_ROW_LENGTH-1,J)         
                                                                        
    end do                                                           
                                                                        
  end if                                                             
                                                                       
  ! ******************************************************************
  ! 1.2 Find gradients over poles. Since the polar rows of TARGET grid   
  !     are ignored in the main routine a limited area is on B-grid      
  !     these are here for completion                                    
  ! ******************************************************************
                                                                        
  do I=1,GRAD_ROW_LENGTH                                            
    DH_DX(I,1) = DH_DX(I,2)                       ! Top row         
    DH_DY(I,1) = DH_DY(I,2)                       ! Top row         
    DH_DX(I,TOPG_ROWS+1) = DH_DX(I,TOPG_ROWS)     ! Bottom row      
    DH_DY(I,TOPG_ROWS+1) = DH_DY(I,TOPG_ROWS)     ! Bottom row      
  end do                                                             
                                                                        
  ! ********************************************************************
  ! 2.0 Calculate SIGMA_XX, SIGMA_YY, SIGMA_XY on meaned target grid     
  !     using TSIGMA held temporarily on meaned full grid in XX,YY or XY 
  ! ********************************************************************
                                                                        
  ! Calculate SIGMA_XX on source grid                                  
                                                                        
  do J=1,TOPG_ROWS+1                                                
    do I=1,GRAD_ROW_LENGTH                                          
      if (DH_DX(I,J)/=RMDI) then
        TSIGMA(I,J)=DH_DX(I,J)*DH_DX(I,J)                             
      else
        TSIGMA(I,J)=RMDI
      endif
    end do                                                           
  end do                                                             
                                                                        
  ! Get SIGMA_XX on target grid                                        
  
  call BOX_SUM( GRAD_ROW_LENGTH, & !
                TOPG_ROWS+1,     & !
                ROW_LENGTH,      & !
                ROWS,            & !                
                LONG_L,          & !
                COLAT_T,         & !
                I_L,             & !
                J_T,             & !
                GLOBAL,          & !
                SIGMA_XX,        & !
                TSIGMA           ) !                
                                                                        
  ! Calculate SIGMA_YY on source grid                                  

  do J=1,TOPG_ROWS+1                                                
    do I=1,GRAD_ROW_LENGTH   
      if (DH_DY(I,J)/=RMDI) then                                       
        TSIGMA(I,J)=DH_DY(I,J)*DH_DY(I,J)                             
      else
        TSIGMA(I,J)=RMDI
      endif
    end do     ! I=1,TOPG_ROW_LENGTH                                 
  end do       ! J=1,TOPG_ROWS+1   
                                                                        
                                                                        
  ! Get SIGMA_YY on target grid                                        
   
  call BOX_SUM( GRAD_ROW_LENGTH, & !
                TOPG_ROWS+1,     & !
                ROW_LENGTH,      & !
                ROWS,            & !                
                LONG_L,          & !
                COLAT_T,         & !
                I_L,             & !
                J_T,             & !
                GLOBAL,          & !
                SIGMA_YY,        & !
                TSIGMA           ) !                
                                                                        
  ! Calculate SIGMA_XY on source grid                                  
                                                                        
  do J=1,TOPG_ROWS+1                                                
    do I=1,GRAD_ROW_LENGTH                                          
      if (DH_DX(I,J)/=RMDI .and. DH_DY(I,J)/=RMDI) then
        TSIGMA(I,J)=DH_DX(I,J)*DH_DY(I,J)                             
      else
        TSIGMA(I,J)=RMDI
      endif
    end do   ! I=1,TOPG_ROW_LENGTH                                 
  end do    ! J=1,TOPG_ROWS+1 
                                                                        
  ! Get SIGMA_XY on target grid

  call BOX_SUM( GRAD_ROW_LENGTH, & !
                TOPG_ROWS+1,     & !
                ROW_LENGTH,      & !
                ROWS,            & !                
                LONG_L,          & !
                COLAT_T,         & !
                I_L,             & !
                J_T,             & !
                GLOBAL,          & !
                SIGMA_XY,        & !
                TSIGMA           ) !                
                                                                        
END SUBROUTINE OROG_GRAD

*/
*/*********************************************************************
*/

! ******************************COPYRIGHT****************************** 
! (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved. 
!                                                                       
! Use, duplication or disclosure of this code is subject to the         
! restrictions as set forth in the contract.                            
!                                                                       
!                Meteorological Office                                  
!                London Road                                            
!                BRACKNELL                                              
!                Berkshire UK                                           
!                RG12 2SZ                                               
!                                                                       
! If no contract has been raised with this copy of the code, the use,   
! duplication or disclosure of it is strictly prohibited.  Permission   
! to do so must first be obtained in writing from the Head of Numerical 
! Modelling at the above address.                                       
! ******************************COPYRIGHT****************************** 
!    SUBROUTINE INTERPBACK                                              
!                                                                       
!    Routine calculates bi-linear mean interpolated orography data      
!    using means of target grid                                         
!    NOT SUITABLE FOR SINGLE COLUMN USE                                 
!                                                                     
!                                                                     
!    SUITABLE FOR ROTATED GRIDS                                       
!                                                                     
!    ORIGINAL VERSION FOR CRAY Y-MP/IBM                               
!    WRITTEN 3/93 J.R.MITCHELL                                        
!                                                                     
!    CODE REVIEWED BY C.WILSON ??/??/??                               
!                                                                     
!    *   VERSION NO. 1 DATED 27/02/90                                 
!    *          COSMOS DSN MEH1.MJUM.GWAVE(OROGFDIF)                  
!    *   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!    *   VERSION 1, DATED 15/01/90                                   
!                                                                    
!    *   SYSTEM TASK:  P65                                           
!                                                                    
!    PURPOSE:  To calculate mean interpolated orography on identical 
!              grid                                                  
!                                                                    
!                                                                    
!  END-------------------------------------------------------------  
                                                                     
SUBROUTINE INTERPBACK( DATA_SRCE,          & ! INOUT
                       OROG_MN,            & ! 
                       START_LONG_SRCE,    & !
                       START_LAT_SRCE,     & !
                       DELTA_LAMBDA_TARG,  & !
                       DELTA_PHI_TARG,     & !
                       DELTA_LAMBDA_SRCE,  & !
                       DELTA_PHI_SRCE,     & !
                       POINTS_LAMBDA_TARG, & !
                       POINTS_PHI_TARG,    & !
                       POINTS_LAMBDA_SRCE, & !
                       POINTS_PHI_SRCE,    & !
                       POINTS_TARG,        & !
                       POINTS_SRCE,        & !
                       LAMBDA_ORIGIN_TARG, & ! 
                       PHI_ORIGIN_TARG,    & ! 
                       GLOBAL              ) ! IN
                                                                     
implicit none                                                  

integer :: POINTS_LAMBDA_SRCE,& ! No. of cols on source grid             
           POINTS_PHI_SRCE,   & ! No. of rows on source grid             
           POINTS_SRCE,       & ! No. of points on source grid           
           POINTS_LAMBDA_TARG,& ! No. of cols on target grid             
           POINTS_PHI_TARG,   & ! No. of rows on target grid             
           POINTS_TARG          ! No. of points on target grid           

real :: START_LAT_SRCE,     & ! temporary stores for target grid details   
        START_LONG_SRCE,    &                                               
        DELTA_LAMBDA_TARG,  & ! longitude resolution of target grid       
        DELTA_PHI_TARG,     & ! latitude resolution of target grid        
        DELTA_LAMBDA_SRCE,  & ! longitude resolution of source grid       
        DELTA_PHI_SRCE,     & ! latitude resolution of source grid        
        LAMBDA_ORIGIN_TARG, & ! longitude origin of traget grid          
        PHI_ORIGIN_TARG       ! phi origin of traget grid                   
                                                                     
real :: OROG_MN(POINTS_LAMBDA_TARG,POINTS_PHI_TARG), & ! mean orography  
        DATA_SRCE(POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE)  ! IN source data
                                                                       
logical, intent(IN) :: GLOBAL  ! T if global
                                                                       
!-- Local variables --
integer :: INDEX_B_L(POINTS_LAMBDA_SRCE), & ! Bottom left gather index        
           INDEX_B_R(POINTS_LAMBDA_SRCE), & ! Bottom right gather index       
           IPHI, ILAMBDA                    ! Loop counters                   
                                                                      
real :: PHI_MN(POINTS_PHI_TARG),      & ! True Lats of target mean grid
       LAMBDA_MN(POINTS_LAMBDA_TARG), & ! True Lons of target mean grid
      WEIGHT_T_R(POINTS_LAMBDA_SRCE), & ! Top right interpolation weight 
     WEIGHT_B_R(POINTS_LAMBDA_SRCE),  & ! Bot right interpolation weight 
     WEIGHT_T_L(POINTS_LAMBDA_SRCE),  & ! Top left interpolation weight  
     WEIGHT_B_L(POINTS_LAMBDA_SRCE),  & ! Bot left interpolation weight  
     PHI_TOPG(POINTS_LAMBDA_SRCE),    & ! lats (true) orography data       
     LAMBDA_TOPG(POINTS_LAMBDA_SRCE), & ! Lons(true) orography data 
     DATA_TEMP(POINTS_LAMBDA_SRCE)      ! data interpolated to high res  
                                                                       
! External routines called                                            
EXTERNAL ::  H_INT_CO, & ! to calculate coefficients for interpolation
             H_INT_BL    ! to perform bi-linear interpolation
!----------------------------------------------------------------------

     ! **************************************************************
     ! 1.0 Calculate arrays for interpolation subroutines, indexed by      
     !     box centres.                                                    
     ! **************************************************************
                                                  
     do ILAMBDA=1,POINTS_LAMBDA_SRCE            
       LAMBDA_TOPG(ILAMBDA) = START_LONG_SRCE + &
                                  (ILAMBDA-1)*DELTA_LAMBDA_SRCE  
     end do                                      
     do ILAMBDA=1,POINTS_LAMBDA_TARG            
       LAMBDA_MN(ILAMBDA) = LAMBDA_ORIGIN_TARG + &
                                 (ILAMBDA-1)*DELTA_LAMBDA_TARG
     end do                    
     do IPHI=1,POINTS_PHI_TARG
       PHI_MN(IPHI) = PHI_ORIGIN_TARG - (IPHI-1)*abs(DELTA_PHI_TARG)

     end do                                                          
                                                                       
     ! *************************************************************
     ! 1.1 Find weights to interpolate mean grid onto data grid then
     !     interpolate                                                     
     ! *************************************************************
                                                                       
     do IPHI=1,POINTS_PHI_SRCE                                      
                                                                       
       PHI_TOPG(1) = START_LAT_SRCE - (IPHI-1)*DELTA_PHI_SRCE       
       do ILAMBDA=2,POINTS_LAMBDA_SRCE                              
         PHI_TOPG(ILAMBDA) = PHI_TOPG(1)                            
       end do                                                        
                                                                       
       call H_INT_CO( INDEX_B_L,          &
                      INDEX_B_R,          &
                      WEIGHT_T_R,         &
                      WEIGHT_B_R,         &
                      WEIGHT_T_L,         &
                      WEIGHT_B_L,         &
                      LAMBDA_MN,          &
                      PHI_MN,             &
                      LAMBDA_TOPG,        &
                      PHI_TOPG,           &
                      POINTS_LAMBDA_TARG, &
                      POINTS_PHI_TARG,    &
                      POINTS_LAMBDA_SRCE, &
                      GLOBAL              )     
                                                                       
       call H_INT_BL( POINTS_PHI_TARG,    &
                      POINTS_LAMBDA_TARG, &
                      POINTS_LAMBDA_SRCE, &
                      INDEX_B_L,          &
                      INDEX_B_R,          &
                      OROG_MN,            &
                      WEIGHT_B_L,         &
                      WEIGHT_B_R,         &
                      WEIGHT_T_L,         &
                      WEIGHT_T_R,         &
                      DATA_TEMP           )
                                                                      
       ! ***********************************************************
       ! 2.0 Subract means interpolated onto data grid from the data          
       ! ***********************************************************
                                                                       
       do ILAMBDA=1,POINTS_LAMBDA_SRCE                              
         DATA_SRCE(ILAMBDA,IPHI)=    &
                  DATA_SRCE(ILAMBDA,IPHI)-DATA_TEMP(ILAMBDA)        
       end do                                                        
                                                                        
     end do                                                            

END SUBROUTINE INTERPBACK

*/
*/**********************************************************************
*/

! ******************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved
!                                                                     
! Use, duplication or disclosure of this code is subject to the       
! restrictions as set forth in the contract.                          
!                                                                     
!                Meteorological Office                                
!                London Road                                          
!                BRACKNELL                                            
!                Berkshire UK                                         
!                RG12 2SZ                                             
!                                                                     
! If no contract has been raised with this copy of the code, the use, 
! duplication or disclosure of it is strictly prohibited.  Permission 
! to do so must first be obtained in writing from the Head of Numerical
! Modelling at the above address.                                    
! ******************************COPYRIGHT******************************
! SUBROUTINE FILTISO--------------------------------------------- 
!                                                                  
!  
!  Purpose: Averages EW polewards of 48 in order to retain approx. 
!           isotropic resolution as meridians converge.
!
!  Documentation: None                                             
!                                                                  
!
! ---------------------------------------------------------------

SUBROUTINE FILTISO( P_DATA, P_FIELD, ROW_LENGTH, ROWS )
                                                                      
implicit none                                                   
                                                                      
integer,intent(IN) :: ROWS,       & ! Number of rows to be updated.       
                      ROW_LENGTH, & ! Number of points per row            
                      P_FIELD       ! Number of points in input field     
                                                                      
real,intent(INOUT) :: P_DATA(P_FIELD) ! Data on p points                    

! -- Define local variables --
real :: WORK(P_FIELD), COSFACTOR(ROWS)

integer :: I,       & !     Horizontal loop index
           ROW_NO,  & !     Number of current loop
           ROW_NOL, & ! }   ROW NOS TO LEFT FOR POINT TO LEFT/RIGHT
           ROW_NOR, & ! }   OF CURRENT POINT
           OFFSET,  & !\  
           OFFSETL, & ! } OFFSET OF PTS TO LEFT/RIGHT OF CURRENT ONE
           OFFSETR    !/  

real :: COSLAT, &   ! INVERSE OF COS LATITUDE ->0 AS NUMBER OF POINTS
                    ! ADDED IS INCREASED
        PI  
! ---------------------------------------------------------------------

      PI = 2.0*ASIN(1.0)

      ! SET POLAR ROWS
      do I=1,ROW_LENGTH
        WORK(I)=P_DATA(I)
        WORK(P_FIELD-ROW_LENGTH+I)= P_DATA(P_FIELD-ROW_LENGTH+I)
      end do

      ! work out for each latitude over how many grid-points
      ! average should be made

      do I = ROW_LENGTH+1, P_FIELD-2*ROW_LENGTH+1, ROW_LENGTH
        ROW_NO      = (I-1)/ROW_LENGTH
        COSFACTOR(ROW_NO)= 2./(3.*ABS(SIN( ROW_NO*PI / (ROWS-1) )))
        if (COSFACTOR(ROW_NO).LT.1.) COSFACTOR(ROW_NO)=1.
      end do

      ! and now perform averaging

      do I = ROW_LENGTH+1, P_FIELD-ROW_LENGTH
        WORK(I)=P_DATA(I)
        ROW_NO   = (I-1)/ROW_LENGTH
        COSLAT   = COSFACTOR(ROW_NO)
        OFFSET   = 0
        do while (COSLAT.GT.1.)
         OFFSET  = OFFSET + 1
         OFFSETL = -OFFSET
         OFFSETR = OFFSET
         ROW_NOL = (I-1+OFFSETL)/ROW_LENGTH
         ROW_NOR = (I-1+OFFSETR)/ROW_LENGTH
         if (ROW_NOL .ne. ROW_NO) OFFSETL = OFFSETL+ROW_LENGTH
         if (ROW_NOR .ne. ROW_NO) OFFSETR = OFFSETR-ROW_LENGTH
         if (COSLAT .gt. 3.) then
           WORK(I)=WORK(I)+P_DATA(I+OFFSETL)+P_DATA(I+OFFSETR)
         else
           WORK(I)=WORK(I)+((COSLAT-1.)/2.*   &
                      (P_DATA(I+OFFSETL)+P_DATA(I+OFFSETR)))
         end if
         COSLAT = COSLAT-2.
        end do
        WORK(I) = WORK(I)/COSFACTOR(ROW_NO)
      end do

      ! ---------------------------------------------------
      ! 2.0  Copy field for output
      ! ---------------------------------------------------

      do I=1,P_FIELD
        P_DATA(I)=WORK(I)
      end do

END SUBROUTINE FILTISO

*/
*/**********************************************************************
*/

*DECK BXBNDSUM
! (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved. 
!                                                                       
! Use, duplication or disclosure of this code is subject to the         
! restrictions as set forth in the contract.                            
!                                                                       
!                Meteorological Office                                  
!                London Road                                            
!                BRACKNELL                                              
!                Berkshire UK                                           
!                RG12 2SZ                                               
!                                                                       
! If no contract has been raised with this copy of the code, the use,   
! duplication or disclosure of it is strictly prohibited.  Permission   
! to do so must first be obtained in writing from the Head of Numerical 
! Modelling at the above address.                                       
! ******************************COPYRIGHT*******************************
!                                                                       
!    SUBROUTINE BOX_BND                                                 
!                                                                       
!    Routine sets up arrays of indexes and boundaries of target grid    
!    boxes relative to a source grid for use  by BOX_SUM so that
!    area weighted means can be calculated.                             
!                                                                       
!    NOT SUITABLE FOR SINGLE COLUMN USE                                 
!                                                                       
!    SUITABLE FOR ROTATED GRIDS                                         
!                                                                       
!    ORIGINAL VERSION FOR CRAY Y-MP/IBM                                 
!    WRITTEN 06/09/91 BY C. WILSON                                      
!                                                                       
!    CODE REVIEWED BY R.SMITH ??/??/??                                  
!                                                                       
!    VERSION NO. 1 DATED 06/09/91                                       
!           COSMOS DSN MS15.CWUM.JOBS(BOXBND1)                        
!    PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,    
!    VERSION 1, DATED 12/09/89                                         
!   History:
!   Version   Date     Comment                                           
!   -------   ----     -------                                           
!   4.0      12/04/95  Imported into Unified model. D.M. Goddard        
!   4.1      12/06/96  Corrections to longitude indexes at boundary.     
!                      D.M. Goddard                                      
!   4.5      12/10/98  Stops LONG_L and COLAT_T becoming negative in     
!                      top row, preventing out of bound in index array.  
!                      Author D.M Goddard                                
!                                                                      
!    SYSTEM TASK:  S1 (part,extension for area mean interpolation)      
!                                                                       
!    PURPOSE: To set up for grid-boxes on a target grid the longitude,  
!             colatitude,and indexes of the overlapping source grid-    
!             boxes at left hand side and top of target grid-boxes.     
!             Both grids are regular lat-long with the same pole        
!             and orientation. Either may be a 'p-grid' (ie with        
!             half-size boxes at the poles) or a 'u-grid' (ie with      
!             regular size boxes everywhere.)                           
!                                                                       
!             NB The units used are "source grid-box lengths"           
!             The area of the target grid_boxes in squared "source      
!             grid-box units" is also returned.                         
!                                                                       
!    DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1                     
!                    BY A.DICKINSON/C WILSON VERSION ??DATED ??/??/91   
!                                                                       
!                                                                     
!  ---------------------------------------------------------------

SUBROUTINE BOX_BND( I_L,             & ! IN
                    LONG_L,          & ! OUT
                    J_T,             & ! IN
                    COLAT_T,         & ! OUT
                    AREA_BOX,        & ! OUT
                    ROW_LENGTH,      & ! IN
                    ROWS,            & ! IN
                    ROW_LENGTH_SRCE, & ! IN
                    ROWS_SRCE,       & ! IN                
                    DELTA_LONG,      & ! IN
                    DELTA_LAT,       & ! IN
                    START_LONG,      & ! IN
                    START_LAT,       & ! IN                 
                    DELTA_LONG_SRCE, & ! IN
                    DELTA_LAT_SRCE,  & ! IN
                    START_LONG_SRCE, & ! IN
                    START_LAT_SRCE,  & ! IN
                    IGRID,           & ! IN
                    IGRID_SRCE,      & ! IN
                    GLOBAL           ) ! IN                                    
                                                                       
implicit none                                                     

integer,intent(IN)::ROW_LENGTH,      & ! Nos of points per target-row
                    ROWS,            & ! Nos of rows of target
                    ROW_LENGTH_SRCE, & ! Nos of pts per row source area 
                    ROWS_SRCE,       & ! Nos of rows of source area         
                    IGRID,           & ! Grid indicators:
                    IGRID_SRCE         !      1=>p-grid, 2=>u-grid

real,intent(IN) ::          &
           DELTA_LONG,      & ! Longitude increment of target grid (deg)  
           DELTA_LAT,       & ! Latitude increment of target grid (deg)    
           START_LONG,      & ! start longitude of centre of first grid-   
                              !  box in target area                         
           START_LAT,       & ! start latitude of centre of first grid-   
                              !   box in target area                         
           DELTA_LAT_SRCE,  & ! Latitude increment of source grid         
           DELTA_LONG_SRCE, & ! Longitude increment of source grid         
           START_LONG_SRCE, & ! Start longitude of centre of first grid-   
                              !  box in source area                         
           START_LAT_SRCE     ! Start latitude of centre of first grid     
                              !  box in source area

logical,intent(IN) :: GLOBAL  ! true if global area required

integer,intent(OUT) :: I_L(ROW_LENGTH+1), &
               ! Index of source box overlapping LHS of target grid-box.
                       J_T(ROWS+1)
               ! Index of source box overlapping top of target grid-box.

real,intent(OUT) ::    LONG_L(ROW_LENGTH +1), &
      ! Left longitude of target grid-box (units of DELTA_LONG_SRCE)
                       COLAT_T(ROWS+1), &
      ! Colatitude of top of target grid-box (units of DELTA_LAT_SRCE)
                       AREA_BOX
      ! area of grid box in sq units of source grid      

!-- DEFINE LOCAL VARIABLES --
real :: EW_BOX, & ! length of grid box in units of source               
                             ! (DELTA_LONG_SRCE)                        
        NS_BOX, & ! height of grid box in units of source               
                             ! (DELTA_LAT_SRCE)                         
        START_LONG_BOX,    & ! start long of first grid box left edge     
        START_COLAT_BOX,   & ! start colat of first grid box top edge     
        START_LONG_BOX_SRCE, & ! start longitude of first grid box left     
                               ! edge ( source)                             
        START_COLAT_BOX_SRCE, & ! start colatitude of first grid box
                                ! top edge (source)
        LONG_OFFSET,       & ! start longitude difference                 
        COLAT_OFFSET         ! start colatitude difference                
                                                                        
integer :: I,J  ! loop counters                                       
real    :: P1,P2                                                        
logical :: LNER                                                      
      LNER(P1,P2) = ((ABS(P1-P2)) .GT. (1.E-5*ABS(P1+P2)))              
                                                                        
 ! ****************************************************************
 ! 1.0 Set target gridbox length,height and area in units of source
 !     grid. Set start 'longitude' and 'colatitude' of first grid         
 !     box also in source units.                                            
 ! ****************************************************************
                                                                        
 EW_BOX = DELTA_LONG/DELTA_LONG_SRCE                                
 NS_BOX = abs(DELTA_LAT)/DELTA_LAT_SRCE!th abs                      
 AREA_BOX = EW_BOX * NS_BOX                                           
                                                                        
 ! **********************************************************
 ! 1.1 Set start colatitude of top and start longitude of LHS
 !     of first boxes on both target and source grids.
 ! **********************************************************
 write(6,*)'BOX_BND: DELTA_LAT_SRCE, DELTA_LAT :', &
            DELTA_LAT_SRCE,DELTA_LAT   
      write(6,*)'BOX_BND: DELTA_LONG_SRCE,DELTA_LONG:', &
            DELTA_LONG_SRCE,DELTA_LONG                               
      write(6,*)'BOX_BND: START_LAT_SRCE, START_LAT :', &
            START_LAT_SRCE,START_LAT  
      write(6,*)'BOX_BND: START_LONG_SRCE,START_LONG:', &
            START_LONG_SRCE,START_LONG                              
      START_LONG_BOX = START_LONG - 0.5*DELTA_LONG                      
      START_COLAT_BOX = (90. - START_LAT) - 0.5*abs(DELTA_LAT)
      START_LONG_BOX_SRCE = START_LONG_SRCE - 0.5*DELTA_LONG_SRCE       
      START_COLAT_BOX_SRCE = (90. - START_LAT_SRCE) - 0.5*DELTA_LAT_SRCE

     if (GLOBAL) then
       if (IGRID_SRCE .eq. 1 .and. LNER(START_LAT_SRCE,90.)) then            
         write(6,*)'BOX_BND: source grid not global'                   
         write(6,*)'BOX_BND: (lat starts too far from nth pole)'
         STOP                                                           
       endif                                                            
       if(IGRID_SRCE.EQ.2 .and. &
         LNER(START_LAT_SRCE,(90.-DELTA_LAT_SRCE*0.5))) then            
         write(6,*)'BOX_BND: source grid not global'                   
         write(6,*)'BOX_BND: (lat starts too far from pole)'
       endif                                                            
     endif                                                             
! JMG 17.3.14: Formerly there was a check here on the separation of the
! N limits of the two grids. It triggered a warning but the condition
! usually occurs with our NH grid and appears harmless.
                                                                        
      LONG_OFFSET  =(START_LONG_BOX-START_LONG_BOX_SRCE)/      &        
                                  DELTA_LONG_SRCE                       
      COLAT_OFFSET = (START_COLAT_BOX - START_COLAT_BOX_SRCE)/  &       
                     DELTA_LAT_SRCE                                     
                                                                        
      if (.not.GLOBAL) then                                             
        if (LONG_OFFSET.LT.0.0) then                                      
          write(6,*)'BOX_BND: long_offset=',LONG_OFFSET,&
                    ': Reset to 0.0'
          LONG_OFFSET = 0.0                                               
        endif                                                             
        if (COLAT_OFFSET.LT.0.0) then                                     
          write(6,*)'BOX_BND: colat_offset=',COLAT_OFFSET,&
                    ': Reset to 0.0'
          COLAT_OFFSET = 0.0
        endif
      endif
      !  ************************************************************
      ! 2.0 Set grid box left longitudes, top colatitudes and indices        
      ! *************************************************************
                                                                        
      do I=1,ROW_LENGTH + 1
        LONG_L(I) = LONG_OFFSET + (I-1)*EW_BOX                          
        if (GLOBAL .and. LONG_L(I) .lt. 0.0) then
          LONG_L(I) = LONG_L(I) + real(ROW_LENGTH_SRCE)                     
        else if (GLOBAL .and. LONG_L(I) .ge. ROW_LENGTH_SRCE) then
          LONG_L(I) = LONG_L(I) - real(ROW_LENGTH_SRCE)                     
        end if
        I_L(I) = LONG_L(I) +1                                           
      end do

      if (LONG_L(1) .lt. 0.0) LONG_L(1) = 0.0                                 
                                                                        
      write(6,*)'BOX_BND: this is I_L..'
      write(6,*) I_L                                                    
      write(6,*)'BOX_BND: this is LONG_L..'
      write(6,*) LONG_L                                                 
                                                                        
      COLAT_T(1) = COLAT_OFFSET
      if (GLOBAL .and. IGRID .eq. 1) then
        COLAT_T(1) = (0.0  - START_COLAT_BOX_SRCE)/DELTA_LAT_SRCE       
      end if
      if (COLAT_T(1) .lt. 0.0) COLAT_T(1)=0.0                            
      J_T(1) = COLAT_T(1) + 1                                           
      do J=2,ROWS+1                                                 
        COLAT_T(J) = COLAT_OFFSET  + (J-1)*NS_BOX                       
        J_T(J) = COLAT_T(J) + 1                                         
      end do

      ! ROWS+1 ie bottom boundary                                          
      if (GLOBAL) then
        if (IGRID .eq. 1) &
                   COLAT_T(ROWS+1) = COLAT_T(ROWS+1)-0.5*NS_BOX      
        J_T(ROWS+1) = ROWS_SRCE                                         
      else                                                              
        if (J_T(ROWS+1) .gt. ROWS_SRCE) then
          if (COLAT_T(ROWS+1) .gt. real(ROWS_SRCE)) then
            write(6,*)'BOX_BND: target area larger than source area'   
            STOP                                                        
          else
            J_T(ROWS+1) = ROWS_SRCE                                     
          end if
        end if
      end if

      write(6,*)'BOX_BND: finished.'
      write(6,*)''

 END SUBROUTINE BOX_BND

*/
*/**********************************************************************
*/

! (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
!
! Use, duplication or disclosure of this code is subject to the
! restrictions as set forth in the contract.
!
!                Meteorological Office
!                London Road                                            
!                BRACKNELL                                              
!                Berkshire UK                                          
!                RG12 2SZ                                               
!                                                                       
! If no contract has been raised with this copy of the code, the use,   
! duplication or disclosure of it is strictly prohibited.  Permission   
! to do so must first be obtained in writing from the Head of Numerical 
! Modelling at the above address.                                       
! ******************************COPYRIGHT****************************** 
!
!    SUBROUTINE BOX_SUM2                                                
!                                                                       
!    NOT SUITABLE FOR SINGLE COLUMN USE                                 
!                                                                       
!                                                                       
!    SUITABLE FOR ROTATED GRIDS                                         
!                                                                       
!    ORIGINAL VERSION FOR CRAY Y-MP/IBM                                 
!    WRITTEN 12/07/91 BY C. WILSON                                      
!                                                                       
!    CODE REVIEWED BY R.SMITH ??/??/??                                  
!                                                                       
!    VERSION NO. 2 DATED 16/09/91                                       
!           COSMOS DSN MS15.CWUM.JOBS(BOXSUM2)                          
!    PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,     
!    VERSION 1, DATED 12/09/89                                          
! History:                                                              
! Version   Date     Comment                                            
! -------   ----     -------                                            
! 4.0      12/04/95  Imported into Unified model. D.M. Goddard          
! 4.1      12/06/96  Corrections for zonal means. D.M. Goddard          
! 4.4      30/09/97  Corrections for portable model. D.M. Goddard       
!                                                                       
!    SYSTEM TASK:  S1 (part,extension for area mean interpolation)      
!                                                                       
!    PURPOSE:                                                           
!    Routine sums contributions from gridboxes for source data on a     
!    regular lat-long grid to form means for gridboxes of a regular     
!    lat-long grid specified as target.                                 
!    Both grids are defined with the same pole and orientation;         
!    the original data must be interpolated onto a rotated              
!    grid ,if the target grid is a rotated grid, BEFORE calling this    
!    routine.                                                           
!    The algorithms are general and will cope with either a finer       
!    or coarser resolution source grid.                                 
!   !> Modified to only sum box if contents greater than zero and then
!   !> return the contents of 'WEIGHT' which holds the number of counts
!   !> per box which had values greater than zero. This allows a 
!   !> determination of the number of land points in a box for purposes
!   !> of generating a land mask. On return, all coastal boxes must be
!   !> weighted for the fraction of ocean surface (height=0) in the box.
!  
!    DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1
!                    BY A.DICKINSON/C WILSON VERSION ??DATED ??/??/91
!  
!  
! -------------------------------------------------------------

SUBROUTINE BOX_SUM2( SOURCE_ROW_LENGTH, & ! IN
                     SOURCE_ROWS,       & ! IN
                     ROW_LENGTH,        & ! IN
                     ROWS,              & ! IN
                     LONG_L,            & ! IN
                     COLAT_T,           & ! IN
                     I_L,               & ! IN
                     J_T,               & ! IN
                     GLOBAL,            & ! IN
                     BOXSUM,            & ! OUT
                     WEIGHT,            & ! OUT
                     SOURCE             ) ! IN
                                                                        
implicit none

integer,intent(IN) :: &
           SOURCE_ROW_LENGTH, & ! Nos of pts per row (source data)
                                !  on rotated grid if necessary          
           SOURCE_ROWS,       & ! Nos of rows of source data         
                                !  on rotated grid if necessary          
           ROW_LENGTH,        & ! Nos of pts per row target area  
           ROWS,              & ! Nos of rows of target area         
           I_L(ROW_LENGTH+1), & ! Index of 1st source box to
                                !  overlap with lhs of target box
           J_T(ROWS+1)          ! Index of first source gridbox to
                                !  overlap top of target gridbox
                                                                       
   ! N.B.I_L(I)   is the first source gridbox to overlap LH side
   !                of target box I of a row.
   !     I_L(I+1) is the last source gridbox to overlap RH side
   !                of target box I of a row.
   !     J_T(J)   is the first source gridbox to overlap top of
   !                target box on row J.
   !     J_T(J+1) is the last source gridbox to overlap bottom of
   !                target box on row J.
   !                                                                       
   ! REAL value of:-                                                      
   !     I_L(I)  is also used to measure the 'longitude' of the
   !               RHS of the source gridbox.
   !     J_T(J)  is also used to measure the 'colatitude' of
   !               the bottom of the source gridbox.
                                                                        
real,intent(IN) :: &
         SOURCE(SOURCE_ROW_LENGTH,SOURCE_ROWS), & ! source data
         LONG_L(ROW_LENGTH +1),                 & 
            !  Left longitude of cell (units of source cell EW length) 
         COLAT_T(ROWS +1)
            ! Colatitude of top of cell (units of source cell NS length)

real, intent(OUT) :: &
        BOXSUM(ROW_LENGTH,ROWS), & ! Sum of data on target grid     
        WEIGHT(ROW_LENGTH,ROWS)    ! total counts of land per box
                                                                        
logical,intent(IN) :: GLOBAL ! true if global area required          


!-- DEFINE LOCAL VARIABLES --
real :: EW_SUM(ROW_LENGTH),     & ! summed WE source data               
        EW_WEIGHT(ROW_LENGTH),  & ! summed WE weights for source data   
        BOX_WEIGHT(ROW_LENGTH), & ! summed weights for target boxes    
        RH_BOX

real,parameter :: RMDI = -1d0*(32768.)**2    ! same mdi as in box_sum

integer :: I,J,I1,I2,IT,J1,J2,JT,K ! loop counters                   

logical :: MDI(ROW_LENGTH)
    ! track which target cells contain (any!) source cells with mdi.
!-------------------------------------------------------------------

! ******************************************************************
! 1.0 Sum source boxes (whole & partial) contributions to target box 
! ******************************************************************
                                                                        
do J=1,ROWS ! Loop over target Rows

  ! The top/bottom source-rows that cover the target row. 
  J1 = J_T(J)                                                     
  J2 = J_T(J+1)                                                   
                   
  ! Loop over target columns: give target-cell zero initial weighting
  do I=1,ROW_LENGTH
    BOX_WEIGHT(I)=0.0
    MDI(I) = .FALSE.
  end do                                                          
                                                                        
  do JT=J1,J2 ! Loop over the source Rows within target-cell row
                                                                        
    ! ***********************************************************
    ! 1.1 Sum  EW (whole & partial) contributions to target cells
    ! ***********************************************************

    ! Loop over target columns
    do I=1,ROW_LENGTH 
      ! Initial ew-slice value/weight of zero.
      EW_SUM(I)    = 0.0                                               
      EW_WEIGHT(I) = 0.0
    end do

    ! Loop over target columns
    do I=1,ROW_LENGTH

      ! The left/right source columns that cover the target column.
      I1 = I_L(I)
      I2 = I_L(I+1)

      ! If the source box spans zero longitude need to split summation
      if (I1.GT.I2 .and. GLOBAL) then

        ! Loop over source-cells in ew-slice.
        do IT=I1,SOURCE_ROW_LENGTH

          ! Check source isn't missing data and target-cell doesn't
          ! already contain missing data..
          if ( (SOURCE(IT,JT) .eq. RMDI) .or. MDI(I) ) then 

            MDI(I) = .TRUE.

          else  ! haven't found mdi in target-cell yet, so continue..

            if (IT .eq. I1) then ! Left side partial contribution

              RH_BOX = LONG_L(I+1)                                  
              if ( RH_BOX .lt. LONG_L(I) ) &
                                     RH_BOX = RH_BOX + SOURCE_ROW_LENGTH
              EW_WEIGHT(I) = EW_WEIGHT(I) +                &
                                    ( min(real(I1),RH_BOX) - LONG_L(I) )
              EW_SUM(I)    = EW_SUM(I) + SOURCE(I1,JT) *   &
                                    ( min(real(I1),RH_BOX) - LONG_L(I) )

            else ! Whole contributions

              EW_WEIGHT(I) = EW_WEIGHT(I) + 1.0                    
              EW_SUM(I)    = EW_SUM(I)    + SOURCE(IT,JT)               

            end if  ! close check on left/right/whole contributions.

          end if  ! close check on missing-data in source

        end do  ! close loop over source-cells in ew-slice

        ! Loop over source-cells in ew-slice.
        do IT=1,I2

          ! Check source isn't missing data and target-cell doesn't
          ! already contain missing data..
          if ( (SOURCE(IT,JT) .eq. RMDI) .or. MDI(I) ) then 

            MDI(I) = .TRUE.

          else  ! haven't found mdi in target-cell yet, so continue..

            if (IT .eq. I2) then ! Right side partial contribution

              EW_WEIGHT(I) = EW_WEIGHT(I) + (LONG_L(I+1)-(I2-1))    
              EW_SUM(I)    = EW_SUM(I)    + SOURCE(I2,JT) *      &   
                                                  ( LONG_L(I+1)-(I2-1) )

            else ! Whole contributions

              EW_WEIGHT(I) = EW_WEIGHT(I) + 1.0                    
              EW_SUM(I)    = EW_SUM(I)    + SOURCE(IT,JT)               
 
            end if  ! close check on left/right/whole contributions.
          end if  ! close check on missing-data in source
        end do  ! close loop over source-cells in ew-slice

      else if (I1.LT.I2) then ! no zero meridian crossing

        do IT=I1,I2  ! Loop over source-cells in ew-slice

          ! Check source isn't missing data and target-cell doesn't
          ! already contain missing data..
          if ( (SOURCE(IT,JT) .eq. RMDI) .or. MDI(I) ) then 

            MDI(i) = .TRUE.

          else ! haven't found mdi in target-cell yet, so continue..

            if (IT .eq. I1) then ! Left side partial contribution

              EW_WEIGHT(I) = EW_WEIGHT(I) +                   &
                               ( min(real(I1),LONG_L(I+1)) - LONG_L(I) )
              EW_SUM(I)    = EW_SUM(I)    + SOURCE(I1,JT) *   & 
                               ( min(real(I1),LONG_L(I+1)) - LONG_L(I) )
                                                                        
            else if (IT .eq. I2) then ! Right side partial contrib

              EW_WEIGHT(I) = EW_WEIGHT(I)+ (LONG_L(I+1)-(I2-1))     
              EW_SUM(I)    = EW_SUM(I)   + SOURCE(I2,JT) * &
                                                  ( LONG_L(I+1)-(I2-1) )
                                                                        
            else ! Whole contributions

               EW_WEIGHT(I) = EW_WEIGHT(I) + 1.0                      
               EW_SUM(I)    = EW_SUM(I)    + SOURCE(IT,JT)
                                                                       
            end if  ! close check on left/right/whole contributions.
          end if  ! close check on missing-data in source
        end do  ! close loop over source-cells in ew-slice
                                                                        
      else ! Zonal mean no need to average in EW direction

        do K=1,ROW_LENGTH
          EW_WEIGHT(K)=1.0
          EW_SUM(K)=SOURCE(1,JT)
        end do

      end if ! Close check on crossing meridian.

    end do ! End loop over target columns

    ! *****************************************************************
    ! 1.3 Add summed EW box contributions into rows J target grid boxes  
    ! *****************************************************************
                                                                        
    if (JT .eq. J1) then         ! Top row

      do I=1,ROW_LENGTH   ! Loop over target columns
        BOXSUM(I,J)   = EW_SUM(I) *                      &
                             ( min(real(J1),COLAT_T(J+1)) - COLAT_T(J) )
        BOX_WEIGHT(I) = BOX_WEIGHT(I) + EW_WEIGHT(I) *   &
                             ( min(real(J1),COLAT_T(J+1)) - COLAT_T(J) )
      end do

    else if (JT .eq. J2) then  ! Bottom of row J

      do I=1,ROW_LENGTH   ! Loop over target columns
        BOXSUM(I,J)   = BOXSUM(I,J) + (1-(J2 - COLAT_T(J+1))) *EW_SUM(I)
        BOX_WEIGHT(I) = BOX_WEIGHT(I) +   &                       
                                  (1-(J2 - COLAT_T(J+1))) * EW_WEIGHT(I)
      end do

    else                       ! Whole contributions to row J

      do I=1,ROW_LENGTH   ! Loop over target columns
        BOXSUM(I,J)   = BOXSUM(I,J)   + EW_SUM(I)                     
        BOX_WEIGHT(I) = BOX_WEIGHT(I) + EW_WEIGHT(I) 
      end do

    end if
                                                                        
  end do ! End loop over Source rows
                                                                        
  do I=1,ROW_LENGTH ! Loop over target columns.

    ! Check we have some weighting, and that the target 
    ! cell does not contain any missing data.
    if ( (BOX_WEIGHT(I) .ne. 0.0) .and. .not.MDI(I) ) then
      BOXSUM(I,J) = BOXSUM(I,J) / BOX_WEIGHT(I) 
      WEIGHT(I,J) = BOX_WEIGHT(I)
    else                                                         
      BOXSUM(I,J) = RMDI
    end if                                                      

  end do ! End loop over target columns.
                                                                        
end do ! End loop over Target rows

write(6,*)'BOX_SUM2: finished.'
write(6,*)''                                    
                                                                        
END SUBROUTINE BOX_SUM2

*/
*/********************************************************************
*/

