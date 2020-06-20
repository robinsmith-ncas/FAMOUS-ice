*ID CO2_COUPLE
*/ Enable coupling of ocean-atmosphere CO2 for FAMOUS. 4.5.1 code as it was
*/ only worked for HadCM3LC, where ocean and atmosphere grids are identical
*/  - coupling needs to follow other FAMOUS fields in terms of interpolation,
*/ area averaging and coastal tiling. In addition, the GATHER_FIELD in 
*/ SWAPO2A just hung on HPCx/Ruby for nprocs>1
*/
*DECLARE INITA2O1
*D CCN1F405.77,CCN1F405.80
*DECLARE SWAPO2A2
*D CCN1F405.251,253
c
c GATHER_FIELD doesn't appear to work for input lengths that aren't
c glsize(1),glsize(2),lasize(1),lasize(2) - this call has been changed 
c to a GENERAL_GATHER. We need the D1 index number to do this
c 
        do i=1,N_OBJ_D1_MAX
          if (d1_addr(d1_address,i,ocean_sm).eq.JO_co2flux) exit
        end do

        CALL GENERAL_GATHER_FIELD(D1(JO_co2flux),CO2FLUX,
     &  (lasize(1)-2)*lasize(2),(glsize(1)-2)*glsize(2),
     &   D1_ADDR(d1_object_type,i,ocean_sm),gather_pe,
     &   ICODE,CMESSAGE)
*D CCN1F405.262
c the resultant GATHERed field appears to be shifted one j line
c up relative to other imt*jmt fields - so shift things down
          do j=1,g_jmt-1
*D CCN1F405.265 
              O_CO2FLUX(i+(j-1)*g_imt) = CO2FLUX(i+(j)*(g_imt-2))
*I CCN1F405.267
          do i=1,g_imt-2
            O_CO2FLUX(i+(g_jmt-1)*g_imt) = RMDI
          enddo  ! i
*I CCN1F405.271
c nasty - at the beginning of a model run, the pure-diagnostic CO2 flux
c field may not have any data. It's initialised as RMDI, but these aren't 
c detected over the sea and they destroy the atmospheric pCO2.
          do i=1,g_imt*g_jmt
            if (stepim(o_im).eq.0
     &     .and.o_flddepc(i).gt.0
     &     .and.o_co2flux(i).eq.RMDI
     &         ) then
              o_co2flux(i)=0.
              write(6,*)"RESETTING OCEAN CO2 FLUX TO 0"
            end if
          end do
*I CCN1F405.302 
c need to multiply by the sea-fraction to end up conserving area total
c co2 passed between O and A
        if (mype.eq.gather_pe) then
        DO I=1,G_P_FIELD
          A_CO2FLUX(I)=A_CO2FLUX(I)*(1-FLAND_GLO(I))
        END DO
        end if
*DECLARE TRANA2O1
*I CCN1F405.185
     +,ATMCO2_INV(IMT,JMT)   ! atm. CO2
*D CCN1F405.187
*D CCN1F405.191,CCN1F405.193
*B CCN1F405.199
*IF -DEF,TRANGRID 
*I CCN1F405.200
*ELSE
c follow method used for other fields to get a conserved atmos CO2 
c on the finer ocean grid.
c
        CALL H_INT_BL(JROWS,ICOLS,SEAPOINTS,INDEXBL,INDEXBR,atmco2
     &  ,             WEIGHTBL,WEIGHTBR,WEIGHTTL,WEIGHTTR,WORK)
        CALL POST_H_INT(NCOASTAL,INDEXA,INDEXO,ATPOINTS,atmco2,AMINT
     &  ,SEAPOINTS,WORK,OCPOINT,RMDI,OCPOINTS,atmco2_OUT)
        call do_areaver(irt,jmt,imt,invert,atmco2_OUT,icols,jrows
     &  ,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &  ,atmco2,adjustment,icode,cmessage)
        if (invert.eqv..true.) then
          do j=1,jmt
          do i=1,imt
            atmco2_inv(i,j)=atmco2_out(i,jmt+1-j)
          enddo
          enddo
          call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &      count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &      atmco2_inv,dummy,icode,cmessage)
          do j=1,jmt
          do i=1,imt
            atmco2_out(i,j)=atmco2_inv(i,jmt+1-j)
          enddo
          enddo
        else
          call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &      count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &      atmco2_out,dummy,icode,cmessage)
        endif
*ENDIF
*I CCN1F405.207
            else
              ATMCO2_OUT(i,j) = RMDI
*D CCN1F405.212 

*DECLARE TRANO2A1
*D CCN1F405.326
*B CCN1F405.332
*IF -DEF,TRANGRID
*I CCN1F405.333
*ELSE
c do simple area-averaging onto atmos grid 
c
        CALL COPYO2A(IMT,JMT,co2fluxin,.FALSE.,OMASK
     &,.FALSE.,INVERT,IMT,WORK_O)

        CALL DO_AREAVER(IRT,JMT,
     &                  IMT,.FALSE.,work_o,
     &                  ICOLS,JROWS,
     &                  COUNT_A,BASE_A,
     &                  ICOLS,.FALSE.,AMASKTP,
     &                  POINT_O,WEIGHT,0,
     &                  co2fluxout,DUMMY,
     &                  ICODE,CMESSAGE)

*ENDIF
*D CCN1F405.343
