#ifdef CPRIBM
@PROCESS ALIAS_SIZE(107374182)
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_timestep.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glint_timestep
  !*FD timestep of a GLINT instance

  use glint_type
  use glint_constants

  implicit none

  private
  public glint_i_tstep, glint_i_tstep_gcm, &
         get_i_upscaled_fields, get_i_upscaled_fields_gcm

  interface glint_lapserate
     module procedure glint_lapserate_dp, glint_lapserate_sp
  end interface

contains

  subroutine glint_i_tstep(time,            instance,       &
                           g_temp,          g_temp_range,   &
                           g_precip,        g_zonwind,      &
                           g_merwind,       g_humid,        &
                           g_lwdown,        g_swdown,       &
                           g_airpress,      g_orog,         &
                           g_orog_out,      g_albedo,       &
                           g_ice_frac,      g_veg_frac,     &
                           g_snowice_frac,  g_snowveg_frac, &
                           g_snow_depth,                    &
                           g_water_in,      g_water_out,    &
                           t_win,           t_wout,         &
                           ice_vol,         out_f,          &
                           orogflag,        ice_tstep)

    !*FD Performs time-step of an ice model instance. 
    !*FD Note that input quantities here are accumulated/average totals since the
    !*FD last call.
    !*FD Global output arrays are only valid on the main task.
    !
    
    use glimmer_paramets
    use glimmer_physcon, only: rhow,rhoi
    use glimmer_log
    use glide
    use glissade
    use glide_io
    use glint_mbal_coupling
    use glint_io
    use glint_mbal_io
    use glint_routing
    use glide_diagnostics
    use parallel, only: tasks, main_task
#ifdef BISICLES_CDRIVER
    use bisicles_cdriver
#endif
    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         !*FD Current time in hours
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)

    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)

    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_veg_frac   !*FD Output veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowice_frac !*FD Output snow-ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowveg_frac !*FD Output snow-veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snow_depth !*FD Output snow depth (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (m)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output   
    logical,                intent(in)   :: orogflag     !*FD Set if we have new global orog
    logical,                intent(out)  :: ice_tstep    !*FD Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(rk),dimension(:,:),pointer :: routing_temp => null() ! temporary array for flow routing
    real(rk),dimension(:,:),pointer :: accum_temp   => null() ! temporary array for accumulation
    real(rk),dimension(:,:),pointer :: ablat_temp   => null() ! temporary array for ablation
    integer, dimension(:,:),pointer :: fudge_mask   => null() ! temporary array for fudging
    real(sp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(sp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux
    
!+seg
    real(sp),dimension(:,:),pointer :: thck_seg => null()
    real(sp)   :: tempth_seg
!-seg
    real(rk) :: start_volume,end_volume,flux_fudge            ! note: only valid for single-task runs
    integer :: i, j, k, ii, jj, nx, ny, il, jl, ig, jg

    logical :: gcm_smb   ! true if getting sfc mass balance from a GCM

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep, time Shush=', time
       write(stdout,*) 'next_time =', instance%next_time
    end if

    ! Check whether we're doing anything this time.

    if (time /= instance%next_time) then
       return
    else
       instance%next_time = instance%next_time + instance%mbal_tstep
    end if

    ! Assume we always need this, as it's too complicated to work out when we do and don't

    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)
!+seg
    call coordsystem_allocate(instance%lgrid, thck_seg)
!-seg

    ice_tstep = .false.

    ! Downscale input fields from global to local grid
    ! This subroutine computes instance%acab and instance%artm, the key inputs to GLIDE.

    call glint_downscaling(instance,                  &
                           g_temp,     g_temp_range,  &
                           g_precip,   g_orog,        &
                           g_zonwind,  g_merwind,     &
                           g_humid,    g_lwdown,      &
                           g_swdown,   g_airpress,    &
                           orogflag)

    ! ------------------------------------------------------------------------  
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

    call glide_get_usurf(instance%model, instance%local_orog)
#ifdef BISICLES_CDRIVER
    if (instance%model%options%whichdycore /= DYCORE_BISICLES_CDRIVER) then
       !no need to do this with BISICLES, and it doesn't work in parallel
#endif
       call glint_remove_bath(instance%local_orog,1,1)
#ifdef BISICLES_CDRIVER
    end if
#endif

    ! ------------------------------------------------------------------------  
    ! Adjust the surface temperatures using the lapse-rate, by reducing to
    ! sea-level and then back up to high-res orography
    ! ------------------------------------------------------------------------  

    call glint_lapserate(instance%artm, real(instance%global_orog,rk), real(-instance%data_lapse_rate,rk))
    call glint_lapserate(instance%artm, real(instance%local_orog,rk),  real(instance%lapse_rate,rk))

    ! Process the precipitation field if necessary ---------------------------
    ! and convert from mm/s to m/s

    call glint_calc_precip(instance)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Do accumulation --------------------------------------------------------

    call glint_accumulate(instance%mbal_accum, time, instance%artm, instance%arng, instance%prcp, &
                          instance%snowd, instance%siced, instance%xwind, instance%ywind, &
                          instance%local_orog, real(thck_temp,rk), instance%humid,    &
                          instance%swdown, instance%lwdown, instance%airpress)

    ! Initialise water budget quantities to zero. These will be over-ridden if
    ! there's an ice-model time-step

    t_win=0.0       ; t_wout=0.0
    g_water_out=0.0 ; g_water_in=0.0

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'time =', time
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
    end if

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep
       end if

       ice_tstep = .true.

       ! Prepare arrays for water budgeting

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
          call coordsystem_allocate(instance%lgrid, accum_temp)
          call coordsystem_allocate(instance%lgrid, ablat_temp)
          accum_temp = 0.0
          ablat_temp = 0.0
       end if

       ! Calculate the initial ice volume (scaled and converted to water equivalent)
       ! start_volume is only valid for single-task runs (this is checked in the place
       ! where it is used)

       call glide_get_thk(instance%model, thck_temp)
       thck_temp = thck_temp*real(rhoi/rhow)
       start_volume = sum(thck_temp)

       ! ---------------------------------------------------------------------
       ! Timestepping for the dynamic ice sheet model
       ! ---------------------------------------------------------------------

       do i = 1, instance%n_icetstep

          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'Ice sheet timestep, iteration new =', i
          end if

          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp*real(rhoi/rhow)

          ! Get latest upper-surface elevation (needed for masking)
          call glide_get_usurf(instance%model, instance%local_orog)

#ifdef BISICLES_CDRIVER
          if (instance%model%options%whichdycore /= DYCORE_BISICLES_CDRIVER) then
             !no need to do this with BISICLES, and it doesn't work in parallel
#endif
             call glint_remove_bath(instance%local_orog,1,1)
#ifdef BISICLES_CDRIVER
          end if
#endif

          ! Get the mass-balance, as m water/year 
          call glint_get_mbal(instance%mbal_accum, instance%artm, instance%prcp, instance%ablt, &
                              instance%acab, instance%snowd, instance%siced, instance%mbal_accum_time)

          ! Mask out non-accumulation in ice-free areas
          ! SEG 16Sep17
          ! Getting rid of SMB that might be added to non-ice box
          !where(thck_temp <= 0.0 .and. instance%acab < 0.0)
          where(thck_temp <= 0.0) 
             instance%acab = 0.0
             instance%ablt = instance%prcp
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is divided 
          !  by scale_acab = scyr*thk0/tim0 and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          !TODO - Change to dp
          call glide_set_acab(instance%model, instance%acab*real(rhow/rhoi))
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag_global
             jl = instance%model%numerics%jdiag_global
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
             write (stdout,*) 'evolve_ice  ',instance%evolve_ice
             write (stdout,*) 'dycore  ',instance%model%options%whichdycore
          end if

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.0)
             instance%acab = -thck_temp
             instance%ablt =  thck_temp
          end where

          instance%glide_time = instance%glide_time + instance%model%numerics%tinc

          ! call the dynamic ice sheet model (provided the ice is allowed to evolve)

          if (instance%evolve_ice == EVOLVE_ICE_TRUE) then

             if (instance%model%options%whichdycore == DYCORE_GLIDE) then

!+seg
                if (GLC_DEBUG .and. tasks==1) then
                  call glide_get_thk(instance%model, thck_seg)
                  write(stdout,*) 'Thick before p1 ',sum(thck_seg)
                endif
!-seg
                call glide_tstep_p1(instance%model,instance%glide_time)
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                call glide_get_thk(instance%model, thck_seg)
                write(stdout,*) 'Thick before p2 ',sum(thck_seg)
                endif
!-seg
                call glide_tstep_p2(instance%model)
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                call glide_get_thk(instance%model, thck_seg)
                write(stdout,*) 'Thick before p3 ',sum(thck_seg)
                endif
!-seg
                call glide_tstep_p3(instance%model)
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                call glide_get_thk(instance%model, thck_seg)
                write(stdout,*) 'Thick after p3 ',sum(thck_seg)
                endif
!-seg
#ifdef BISICLES_CDRIVER                
             else if (instance%model%options%whichdycore == DYCORE_BISICLES_CDRIVER) then

                call bisicles_tstep(instance%model, instance%glide_time)
                instance%model%numerics%time = instance%glide_time
#endif
             else   ! glam/glissade dycore

                call glissade_tstep(instance%model,instance%glide_time)

             endif

          endif   ! evolve_ice

          ! Add the calved ice to the ablation field

          call glide_get_calving(instance%model, calve_temp)
          calve_temp = calve_temp * real(rhoi/rhow)

          instance%ablt = instance%ablt + calve_temp/instance%model%numerics%tinc
          instance%acab = instance%acab - calve_temp/instance%model%numerics%tinc

          ! Accumulate for water-budgeting
          if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
             accum_temp = accum_temp + instance%prcp*instance%model%numerics%tinc
             ablat_temp = ablat_temp + instance%ablt*instance%model%numerics%tinc
          endif

          ! write ice sheet diagnostics at specified interval

          call glide_write_diagnostics(instance%model,                  &
                                       instance%model%numerics%time,    &
                                       tstep_count = instance%model%numerics%timecounter)

          ! write netCDf output

          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

       end do   ! n_icestep

       ! Calculate flux fudge factor --------------------------------------------

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then

          ! WJS (1-15-13): I am pretty sure (but not positive) that the stuff in this
          ! conditional will only work right with a single task
          if (tasks > 1) then
             call write_log('The sums in the computation of a flux fudge factor only work with a single task', &
                            GM_FATAL, __FILE__, __LINE__)
          end if
          
          call coordsystem_allocate(instance%lgrid,fudge_mask)

          call glide_get_thk(instance%model,thck_temp)
          end_volume = sum(thck_temp)

          where (thck_temp > 0.0)
             fudge_mask = 1
          elsewhere
             fudge_mask = 0
          endwhere

          flux_fudge = (start_volume + sum(accum_temp) - sum(ablat_temp) - end_volume) / sum(fudge_mask)

          ! Apply fudge_factor

          where(thck_temp > 0.0)
             ablat_temp = ablat_temp + flux_fudge
          endwhere
          
          deallocate(fudge_mask)
          fudge_mask => null()

       endif

       ! Upscale water flux fields ----------------------------------------------
       ! First water input (i.e. mass balance + ablation)

       if (out_f%water_in) then
          call coordsystem_allocate(instance%lgrid, upscale_temp)

          where (thck_temp > 0.0)
             upscale_temp = accum_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call mean_to_global(instance%ups,   &
                              upscale_temp,   &
                              g_water_in,     &
                              instance%out_mask)
          deallocate(upscale_temp)
          upscale_temp => null()
       endif

       ! Now water output (i.e. ablation) - and do routing

       if (out_f%water_out) then
          ! WJS (1-15-13): The flow_router routine (called bolew) currently seems to
          ! assume that it's working on the full (non-decomposed) domain. I'm not sure
          ! what the best way is to fix this, so for now we only allow this code to be
          ! executed if tasks==1.
          if (tasks > 1) then
             call write_log('water_out computation assumes a single task', &
                            GM_FATAL, __FILE__, __LINE__)
          end if

          call coordsystem_allocate(instance%lgrid, upscale_temp)
          call coordsystem_allocate(instance%lgrid, routing_temp)

          where (thck_temp > 0.0)
             upscale_temp = ablat_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call glide_get_usurf(instance%model, instance%local_orog)
          call flow_router(instance%local_orog, &
                           upscale_temp, &
                           routing_temp, &
                           instance%out_mask, &
                           real(instance%lgrid%delta%pt(1),rk), &
                           real(instance%lgrid%delta%pt(2),rk))

          call mean_to_global(instance%ups,   &
                              routing_temp,   &
                              g_water_out,    &
                              instance%out_mask)

          deallocate(upscale_temp,routing_temp)
          upscale_temp => null()
          routing_temp => null()

       endif

       ! Sum water fluxes and convert if necessary ------------------------------

       if (out_f%total_win) then
          if (tasks > 1) call write_log('t_win sum assumes a single task', &
                                        GM_FATAL, __FILE__, __LINE__)

          t_win  = sum(accum_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

       if (out_f%total_wout) then
          if (tasks > 1) call write_log('t_wout sum assumes a single task', &
                                        GM_FATAL, __FILE__, __LINE__)

          t_wout = sum(ablat_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

    end if  ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! We now upscale all fields at once...
    !TODO - This subroutine is called here and also from subroutine glint.
    !       Are both calls needed?

    call get_i_upscaled_fields(instance, g_orog_out, g_albedo, g_ice_frac, g_veg_frac, &
                               g_snowice_frac, g_snowveg_frac, g_snow_depth)

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
       if (tasks > 1) call write_log('ice_vol sum assumes a single task', &
                                     GM_FATAL, __FILE__, __LINE__)

       call glide_get_thk(instance%model, thck_temp)
       ice_vol = sum(thck_temp) * instance%lgrid%delta%pt(1)* &
                                  instance%lgrid%delta%pt(2)
    endif

    ! Tidy up ----------------------------------------------------------------

    if (associated(accum_temp)) then 
       deallocate(accum_temp)
       accum_temp => null()
    end if

    if (associated(ablat_temp)) then
       deallocate(ablat_temp)
       ablat_temp => null()
    end if

    if (associated(calve_temp)) then
       deallocate(calve_temp)
       calve_temp => null()
    end if

    if (associated(thck_temp)) then
       deallocate(thck_temp)
       thck_temp => null()
    endif

  end subroutine glint_i_tstep

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_tstep_gcm(time,            instance,       &
                               ice_tstep,                       &
                               qsmb_g,          tsfc_g,         &
                               topo_g,          gmask,          &
                               gmask3D,                         &
                               gfrac,           gtopo,          &
                               grofi,           grofl,          &
                               ghflx,                           &
                               ice_vol,                         &
                               areas_g,                         &
                               sdep_g,                          &
                               sdep_l,                          &
                               calv_g,                          &
                               qsmb_lev,       & !seg 31/aug/17
                               seg_diag) ! 1/sep/17
    ! SEG gmask3D argument added 27/2/17

    ! Performs time-step of an ice model instance. 
    ! Input quantities here are accumulated/average totals since the last call.
    ! Global output arrays are only valid on the main task.
    !
    use glimmer_paramets
    use glimmer_physcon, only: rhow, rhoi
    use glimmer_log
    use glide
    use glissade
    use glide_io
    use glint_mbal_coupling
    use glint_io
    use glint_mbal_io
    use glide_diagnostics
    use parallel, only: tasks, main_task
#ifdef BISICLES_CDRIVER
    use bisicles_cdriver
#endif

    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         ! Current time in hours
    type(glint_instance), intent(inout)  :: instance     ! Model instance
    logical,                intent(out)  :: ice_tstep    ! Set if we have done an ice time step

    real(rk),dimension(:,:,:),intent(in)  :: qsmb_g    ! Depth of new ice (m)
    real(rk),dimension(:,:,:),intent(in)  :: tsfc_g    ! Surface temperature (C)
    real(rk),dimension(:,:,:),intent(in)  :: topo_g    ! Surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(in)  :: areas_g   ! GCM grid areas on NECs for conservation
    real(rk), dimension(:,:,:),optional,intent(in)  :: sdep_g ! Nonice-snow on GCM NEC grid for downscaling

    integer, dimension(:,:),  optional,intent(in)  :: gmask     ! = 1 where global data are valid, else = 0
    integer, dimension(:,:,:,:),optional,intent(in):: gmask3D   ! = 1 where global data are valid, else = 0  SEG 27/2/17


    real(rk),dimension(:,:,:),optional,intent(out) :: gfrac     ! ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(out) :: gtopo     ! surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofi     ! ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofl     ! liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: ghflx     ! heat flux (W/m^2, positive down)

    real(rk),dimension(:,:),optional,intent(out) :: calv_g       ! icesheet calving flux on GCM grid
    real(rk),dimension(:,:),optional,intent(out) :: sdep_l      ! Nonice-snow depth on icesheet grid
    real(rk),               optional,intent(out) :: ice_vol     ! icesheet volume for this instance

    real(rk),dimension(:),optional,intent(in)    :: qsmb_lev    ! average qsmb fill value seg 31/aug/17

    logical,optional,intent(in) :: seg_diag ! 1/sep/17 seg
    
    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    !TODO - Are these needed?
    real(rk),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(sp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(sp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux
    real(rk),dimension(:,:),pointer :: calv_l       => null()
!+seg
    real(sp),dimension(:,:),pointer :: ice_wall => null() 
    real(sp),dimension(:,:),pointer :: thck_seg => null()
    real(sp),dimension(:,:),pointer :: acab_seg => null()
    real(sp),dimension(:,:),pointer :: mask_seg => null()

    real(sp)   :: tempth_seg

!-seg


    integer :: i, j, k, nx, ny, il, jl, ig, jg
    real(rk), parameter             :: depth_snow2ice =50.    !depth at which nonice-snow accumulation becomes a new icesheet point - and vice-versa
    !+ SEG
    character (len=12) :: cseg1
    !-SEG

    !if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep, time =', time
       write(stdout,*) 'next_time =', instance%next_time
    !end if

    ! Zero outputs
    !TODO - Change to dp
    if (present(gfrac)) gfrac(:,:,:) = 0._rk
    if (present(gtopo)) gtopo(:,:,:) = 0._rk
    if (present(grofi)) grofi(:,:,:) = 0._rk
    if (present(grofl)) grofl(:,:,:) = 0._rk
    if (present(ghflx)) ghflx(:,:,:) = 0._rk
    if (present(ice_vol)) ice_vol = 0._rk
    if (present(calv_g)) calv_g(:,:) = 0._rk
    if (present(sdep_l)) sdep_l(:,:) = 0._rk

    ! Check whether we're doing anything this time.

    write(stdout,*) 'Doing the check'
    if (time /= instance%next_time) then
       return
    else
       instance%next_time = instance%next_time + instance%mbal_tstep
    end if
    write(stdout,*) 'This is a time to do things'
    !TODO - Are these needed?
    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)
!+seg
    call coordsystem_allocate(instance%lgrid, ice_wall)
    call coordsystem_allocate(instance%lgrid, thck_seg)
    call coordsystem_allocate(instance%lgrid, acab_seg)
    call coordsystem_allocate(instance%lgrid, mask_seg)

    if (present(calv_g)) call coordsystem_allocate(instance%lgrid, calv_l)
    


    ice_tstep = .false.

    ! Downscale input fields from global to local grid
    ! This subroutine computes instance%acab and instance%artm, the key inputs to GLIDE.

    write(stdout,*) 'Calling downscaling'
    write(cseg1,'(i12.12)') time
       call glint_downscaling_gcm (instance,              &
                                   qsmb_g,      tsfc_g,   &
                                   topo_g,      cseg1,    &
                                   gmask=gmask,           &
                                   gmask3D=gmask3D,       &
                                   sdep_g=sdep_g,         & 
                                   sdep_l=sdep_l,         &
                                   areas_g=areas_g ,      &
                                   qsmb_lev=qsmb_lev ) !seg 31/aug/17)
      !                             gmask,    &
      !                             gmask3D,               &
      !                             sdep_g,      sdep_l, areas_g , &
      !                             qsmb_lev ) !seg 31/aug/17)
       ! SEG gmask3D argument added 27/2/17




       ! If snowpack on non-ice fraction has got too big, make some new ice
       ! (GCM coverage fractions will get updated during coupling)
       !!  ADD TO ACAB INSTEAD? CAN'T DO THE REVERSE OF THAT AT END OF STEP FOR ICE->SNOW TRANSFORM
       write(stdout,*) 'sdep_g if block'
       if (present(sdep_g)) then
         call glide_get_thk(instance%model,thck_temp)
         !+SEG
         write(stdout,*) 'sdep_l shape ',shape(sdep_l)
         write(stdout,*) 'sdep_g shape ',shape(sdep_g)
         write(stdout,*) 'thck_temp shape ',shape(thck_temp)
         write(stdout,*) 'max g and l',maxval(sdep_g),maxval(sdep_l)
         write(stdout,*) 'min g and l',minval(sdep_g),minval(sdep_l)
         write(stdout,*) 'depth_snow2ice = ',depth_snow2ice
         !STOP
         !-SEG

         !+SEG 18/Oct/17
         ! sdep_l is interpolated to have data everywhere
         ! Should be zero at ice and ocean points
         ! Above output is unadjusted field
         where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask) .OR. &
                thck_temp > 0)
           sdep_l=0.0
         elsewhere
           sdep_l=sdep_l
         endwhere
         !-SEG

! Turn deep snow into ice
! Convert snow depth into snow-depth anom
! In ice depth? no snow -> ice density convert for sdep_l?
         where (sdep_l > depth_snow2ice .AND. thck_temp==0.0)
           thck_temp=sdep_l
           sdep_l=-thck_temp  
         elsewhere
           thck_temp=thck_temp
           sdep_l=0.
         endwhere
         call glide_set_thk(instance%model,thck_temp)
      !! DO I NEED TO CHANGE THE MASK AS WELL?
       end if


    ! ------------------------------------------------------------------------  
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

!TODO: Determine if glint_remove_bath is needed in a CESM run. If so, fix it to work with
!      multiple tasks. 

!!    call glide_get_usurf(instance%model, instance%local_orog)
!!    call glint_remove_bath(instance%local_orog,1,1)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Accumulate acab and artm

    call glint_accumulate_gcm(instance%mbal_accum,   time,        &
                              instance%acab,         instance%artm)

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'time =', time
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
    end if

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep
       end if

       ice_tstep = .true.

       ! ---------------------------------------------------------------------
       ! Timestepping for ice sheet model
       ! ---------------------------------------------------------------------
!+seg
       call glint_get_mbal_gcm(instance%mbal_accum, instance%mbal_accum_time,  &
                                 acab_seg,       instance%artm)

       mask_seg = 1.0
       where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
         acab_seg = 0.0
         mask_seg = 0.0
       endwhere
!+seg
     
       do i = 1, instance%n_icetstep

          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'Ice sheet timestep gcm, iteration =', i
          end if

          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp * real(rhoi/rhow)

          !TODO: Determine if glint_remove_bath is needed in a CESM run. If so, fix it to work with
          !      multiple tasks.  (And decide whether the call is needed both here and above) 
          ! Get latest upper-surface elevation (needed for masking)
!!          call glide_get_usurf(instance%model, instance%local_orog)
!!          call glint_remove_bath(instance%local_orog,1,1)

          call glint_get_mbal_gcm(instance%mbal_accum, instance%mbal_accum_time,  &
                                  instance%acab,       instance%artm)

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.0
          endwhere
          if (GLC_DEBUG .and. main_task) then
             ! Should be the same?
             write (stdout,*) 'ACAB compare 0',i,sum(acab_seg),sum(instance%acab)
          endif

          ! Mask out non-accumulation in ice-free areas

          !TODO - Change to dp
          ! Mask out non-accumulation in ice-free areas
          ! SEG 13Oct17
          !  - changed this previously in the _non_ gcm routine. bugger.
          ! Getting rid of SMB that might be added
          ! to non-ice box
          where(thck_temp <= 0.0) 
             instance%acab = 0.0
          end where
          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'ACAB compare 1',i,sum(acab_seg),sum(instance%acab)
          endif

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is multiplied 
          !  by tim0/(scyr*thk0) and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          ! Adjust glint acab and ablt for output
          ! SEG - moved it 16Oct17
          ! BUT: 
          where (instance%acab < -thck_temp .and. thck_temp > 0.0)
             instance%acab = -thck_temp
          end where
 
          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'ACAB compare 2',i,sum(acab_seg),sum(instance%acab)
          endif

          !TODO - Just rhow/rhoi without 'real'?
          ! Converting to ice density (from water)
          call glide_set_acab(instance%model, instance%acab*real(rhow/rhoi))
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag_global
             jl = instance%model%numerics%jdiag_global
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
             write (stdout,*) 'acab (m/y), artm (C) =', instance%model%climate%acab(il,jl)*rhow/rhoi, instance%model%climate%artm(il,jl)
             write (stdout,*) 'Doing ice walls. Do I need densiity change?'
          end if

          ! Adjust glint acab and ablt for output
          !where (instance%acab < -thck_temp .and. thck_temp > 0.0)
          !   instance%acab = -thck_temp
          !end where

          instance%glide_time = instance%glide_time + instance%model%numerics%tinc

          ! call the dynamic ice sheet model (provided the ice is allowed to evolve)

          if (instance%evolve_ice == EVOLVE_ICE_TRUE) then

             if (instance%model%options%whichdycore == DYCORE_GLIDE) then
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                  call glide_get_thk(instance%model, thck_seg)
                  ! Why is this not changing much in the loop?
                  write(stdout,*) 'Thick before p1 ',sum(thck_seg),i
                  write(stdout,*) 'ACAB1 before p1 ',sum(instance%acab),i
                  call glide_get_acab(instance%model, thck_seg)
                  write(stdout,*) 'ACAB2 before p1 ',sum(thck_seg),i
                endif
!-seg
 
                call glide_tstep_p1(instance%model, instance%glide_time)

!+seg
                if (GLC_DEBUG .and. tasks==1) then
                  call glide_get_thk(instance%model, thck_seg)
                  write(stdout,*) 'Thick before p2 ',sum(thck_seg)
                endif
!-seg
                call glide_tstep_p2(instance%model)
 
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                  call glide_get_thk(instance%model, thck_seg)
                  write(stdout,*) 'Thick before p3 ',sum(thck_seg)
                endif
!-seg
                call glide_tstep_p3(instance%model)
!+seg
                if (GLC_DEBUG .and. tasks==1) then
                  call glide_get_thk(instance%model, thck_seg)
                  write(stdout,*) 'Thick after p3 ',sum(thck_seg)
                endif
!-seg
#ifdef BISICLES_CDRIVER                    
             else if (instance%model%options%whichdycore == DYCORE_BISICLES_CDRIVER) then
                call bisicles_set_header_int(instance%model,"timestamp",time)
                call bisicles_tstep(instance%model, instance%glide_time)

                ! increment time counter - this is from glide_tstep_p3
                instance%model%numerics%timecounter = instance%model%numerics%timecounter + 1
                ! this time counter was advanced within glide_tstep_p1 
                instance%model%numerics%time = instance%glide_time

#endif
             else   ! glam/glissade dycore

                call glissade_tstep(instance%model, instance%glide_time)

             endif

          endif  ! evolve_ice

          ! Add the calved ice to the ablation field

          !TODO - Use this to compute the solid ice runoff,grofi?
          !       Also add basal melting (bmlt) to the liquid runoff, grofl.

          !+seg 30Oct17
          ! Old style Fortran
          ! Turned off 6Nov17
!>!          write(stdout,*) 'Doing ice walls'
!>!          call glide_get_thk(instance%model, thck_temp)
!>!          call do_ice_wall(instance%model%geometry%thkmask,instance%model%geometry%usrf, &
!>!             thck_temp,ice_wall)
!>!          call glide_set_thk(instance%model,thck_temp)
!>!          if (present(calv_g)) then
!>!            call glide_get_calving(instance%model, calve_temp) ! Get present calving
!>!            calve_temp=calve_temp+(ice_wall*thk0)              ! Add ice walls
!>!            call glide_set_calving(instance%model, calve_temp) ! Update
!>!            calv_l=calv_l+(ice_wall*thk0 * real(rhoi/rhow))
!>!            write(111) calv_l
!>!          end if   
          !-seg

          !+seg 17Oct17
          ! Give too-thin ice back to GCM as snowpack on non-ice fraction -----------
          ! This would be ablated if sent back to FAMOUS (1yr timestep)
          ! Use SMB ablation to reduce it?
          if (present(sdep_g)) then
            call glide_get_thk(instance%model, thck_temp)
            call glide_get_thk(instance%model, thck_seg) ! Copy to compute mass diff
            ! This is an anomaly
            ! Adding thck_temp reduces negative anom
            where (thck_temp < depth_snow2ice)
              sdep_l=sdep_l+thck_temp
              thck_temp=0.0
            elsewhere
              thck_temp=thck_temp
            endwhere

            ! Adjust sdep_l
            if (GLC_DEBUG .and. tasks==1) &
              write(stdout,*) 'sdep 1 ',i,sum(sdep_l),cseg1,minval(sdep_l)
! Turn it back into ice again if accum over accell
! The ACAB should then be recalculate (above) to remove ice
! if tahts the case
!+ SEG 20Oct17
! Turned off. 6Nov17
!>!            if (present(seg_diag) .and. main_task) &  !++ seg
!>!              write(216) thck_temp*real(rhoi/rhow)   !++ seg f5
!>!
!>!            if (GLC_DEBUG .and. tasks==1) &
!>!               write(stdout,*) 'sdep max 1',i,maxval(sdep_l)
!>!            where (sdep_l > depth_snow2ice .AND. thck_temp==0.0)
!>!              thck_temp=sdep_l
!>!              sdep_l=0
!>!            elsewhere
!>!              thck_temp=thck_temp
!>!              sdep_l=sdep_l
!>!            endwhere


            if (GLC_DEBUG .and. tasks==1) &
               write(stdout,*) 'sdep max 2',i,maxval(sdep_l)
!            if (i > 1) then
!              ! Negative acab so add (--), reduce anomaly
!              where(acab_seg < 0)
!                sdep_l=sdep_l-acab_seg
!              elsewhere
!                sdep_l=sdep_l
!              endwhere
!              ! Anomaly now, so mostly negative
!              !where(sdep_l < 0)
!              !  sdep_l=0
!              !elsewhere
!              !  sdep_l=sdep_l
!              !endwhere
!            end if  
            if (GLC_DEBUG .and. tasks==1) &
              write(stdout,*) 'sdep 2 ',i,sum(sdep_l)
     
            if (GLC_DEBUG .and. tasks==1) then
              write(stdout,*) 'Icestep ',i,', snow2ice mdiff ',&
                    sum(thck_seg)-sum(thck_temp)
              write(stdout,*) 'Icestep ',i,', sum snow anom ',&
                    sum(sdep_l)
            call glide_set_thk(instance%model,thck_temp)

            endif
          endif
          !-seg
          ! write ice sheet diagnostics at specified interval (model%numerics%dt_diag)

          call glide_write_diagnostics(instance%model,                  &
                                       instance%model%numerics%time,    &
                                       tstep_count = instance%model%numerics%timecounter)

          ! write netCDF output

          print*,'at write io'
          print*,' 1 version of time', instance%model%numerics%timecounter
          print*,' another version of time',instance%model%numerics%time
          
          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

          print*,' IO done'

       end do   ! instance%n_icetstep
 
       ! Give too-thin ice back to GCM as snowpack on non-ice fraction -----------
!       if (present(sdep_g)) then
!         call glide_get_thk(instance%model,thck_temp)
!!+seg
!         if (GLC_DEBUG .and. tasks==1) then
!            write(stdout,*) 'Thick before snow2ice ',sum(thck_temp)
!         endif
!!-seg
!         tempth_seg=sum(sdep_l)
!!         where (thck_temp < depth_snow2ice)
!!           sdep_l=sdep_l+thck_temp
!!           thck_temp=0.0
!!         elsewhere
!!           thck_temp=thck_temp
!!         endwhere
!!         call glide_set_thk(instance%model,thck_temp)
!!+seg
!         if (GLC_DEBUG .and. tasks==1) then
!            call glide_get_thk(instance%model, thck_seg)
!            write(stdout,*) 'Thick after snow2ice ',sum(thck_seg)
!            write(stdout,*) 'Depth sum diff ',sum(sdep_l)-tempth_seg,maxval(sdep_l)
!         endif
!!-seg
!!       !! DO I NEED TO CHANGE THE MASK AS WELL?
!       endif

    end if   ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)
    if (present(calv_g)) then
       ! Upscale calving flux to GCM mean grid as in old GLIMMER
       call mean_to_global(instance%ups,   &
            calv_l, calv_g)
    end if

    !Simple calc. of instance volume - wil only work if all on main task
    call glide_get_thk(instance%model, thck_temp)
    ice_vol = sum(thck_temp) * instance%lgrid%delta%pt(1)*instance%lgrid%delta%pt(2)

    ! Deallocate

    if (associated(calve_temp)) then
       deallocate(calve_temp)
       calve_temp => null()
    end if

    if (associated(thck_temp)) then
       deallocate(thck_temp)
       thck_temp => null()
    endif

    if (associated(calv_l)) then
       deallocate(calv_l)
       calv_l => null()
    endif

  end subroutine glint_i_tstep_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Rewrite to support multiple tasks?

  subroutine glint_remove_bath(orog,x,y)

    ! Sets ocean areas to zero height, working recursively from
    ! a known ocean point.

    use glimmer_log
    use parallel, only : tasks

    real(sp),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    ! Currently, this routine is called assuming point 1,1 is ocean... this won't be true
    ! when running on multiple processors, with a distributed grid
    ! This can't be made a fatal error, because this is currently called even if we have
    ! more than one task... the hope is just that the returned data aren't needed in CESM.
    if (tasks > 1) then
       call write_log('Use of glint_remove_bath currently assumes the use of only one task', &
                      GM_WARNING, __FILE__, __LINE__)
    end if

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y) < 0.0) orog(x,y)=0.0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(sp),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi=(/ -1,1,0,0 /)
    integer,dimension(4) :: yi=(/ 0,0,-1,1 /)
    integer :: ns=4,i

    do i=1,ns
       if (x+xi(i) <= nx.and.x+xi(i) > 0.and. &
            y+yi(i) <= ny.and.y+yi(i) > 0) then
          if (orog(x+xi(i),y+yi(i)) < 0.0) then
             orog(x+xi(i),y+yi(i))=0.0
             call glint_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
          endif
       endif
    enddo

  end subroutine glint_find_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_lapserate_dp(temp,topo,lr)

    !*FD Corrects the temperature field
    !*FD for height, using a constant lapse rate.
    !*FD
    !*FD This the double-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(dp),dimension(:,:), intent(inout) :: temp !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_dp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Remove when we switch to dp

  subroutine glint_lapserate_sp(temp,topo,lr)

    !*FD Corrects the temperature field for height, using a constant lapse rate.
    !*FD
    !*FD This is the single-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(sp),dimension(:,:),intent(inout) :: temp  !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_sp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_calc_precip(instance)

    use glint_precip_param
    use glimmer_log

    !*FD Process precip if necessary

    type(glint_instance) :: instance

    select case (instance%whichprecip)

    case(1)
       ! Do nothing to the precip field

    case(2)
       ! Use the Roe/Lindzen parameterisation
       call glint_precip(instance%prcp, &
                         instance%xwind, &
                         instance%ywind, &
                         instance%artm, &
                         instance%local_orog, &
                         real(instance%lgrid%delta%pt(1),rk), &
                         real(instance%lgrid%delta%pt(2),rk), &
                         fixed_a=.true.)

    case default

       call write_log('Invalid value of whichprecip',GM_FATAL,__FILE__,__LINE__)

    end select

    ! Convert from mm/s to m/s - very important!

    instance%prcp = instance%prcp*0.001

  end subroutine glint_calc_precip

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling(instance,                  &
                               g_temp,     g_temp_range,  &
                               g_precip,   g_orog,        &
                               g_zonwind,  g_merwind,     &
                               g_humid,    g_lwdown,      &
                               g_swdown,   g_airpress,    &
                               orogflag)

    use glint_interp

    !*FD Downscale global fields to the local ice sheet grid

    type(glint_instance) :: instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)
    logical,                intent(in)   :: orogflag

    call interp_to_local(instance%lgrid_fulldomain,g_temp,      instance%downs,localsp=instance%artm)
    call interp_to_local(instance%lgrid_fulldomain,g_temp_range,instance%downs,localsp=instance%arng,z_constrain=.true.)
    call interp_to_local(instance%lgrid_fulldomain,g_precip,    instance%downs,localsp=instance%prcp,z_constrain=.true.)

    if (instance%whichacab==3) then
       call interp_to_local(instance%lgrid_fulldomain,g_humid,   instance%downs,localrk=instance%humid,z_constrain=.true.)
       call interp_to_local(instance%lgrid_fulldomain,g_lwdown,  instance%downs,localrk=instance%lwdown)
       call interp_to_local(instance%lgrid_fulldomain,g_swdown,  instance%downs,localrk=instance%swdown)
       call interp_to_local(instance%lgrid_fulldomain,g_airpress,instance%downs,localrk=instance%airpress,z_constrain=.true.)
    end if

    if (orogflag) call interp_to_local(instance%lgrid_fulldomain,g_orog,instance%downs,localdp=instance%global_orog,z_constrain=.true.)

    if (instance%whichprecip==2 .or. instance%whichacab==3) &
         call interp_wind_to_local(instance%lgrid_fulldomain,g_zonwind,g_merwind,instance%downs,instance%xwind,instance%ywind)

  end subroutine glint_downscaling

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling_gcm (instance,            &
                                    qsmb_g,     tsfc_g,  &
                                    topo_g,    cseg1,    &
                                    gmask,   &
                                    gmask3D,             &
                                    sdep_g,     sdep_l, areas_g , &
                                    qsmb_lev) !seg 31/aug/17
    ! SEG gmask3D argument added 27/2/17
 
    use glimmer_paramets, only: thk0, GLC_DEBUG

    use glint_type
    use glint_interp, only: interp_to_local
    use parallel, only: tasks
    use glide_io, only: glide_get_thk

!    interface
!     subroutine create_halo_seg2(data_in,mask_in,data_out,mask_out,niter_override)
!     real(dp), dimension(:,:), intent(in)  :: data_in
!     integer,  dimension(:,:), intent(in)  :: mask_in
!     real(dp), dimension(:,:), intent(out) :: data_out
!     integer,  dimension(:,:), intent(out) :: mask_out
!     integer, intent(in), optional         :: niter_override
!     end subroutine create_halo_seg2
!   
!     subroutine hval_index(inarr,outval)
!   
!     integer,intent(in)  :: inarr(:)
!     integer,intent(out) :: outval
!     end subroutine hval_index
!
!     subroutine hval_type(inarr,ipos,htype)
!     integer,intent(in)   :: inarr(:)
!     integer,intent(in)   :: ipos
!     integer,intent(out)  :: htype
!     end subroutine hval_type
!
!    end interface

    ! Downscale fields from the global grid (with multiple elevation classes)
    ! to the local ice sheet grid.
    ! 
    ! This routine is used for downscaling when the surface mass balance is
    ! computed in the GCM land surface model.

    type(glint_instance), intent(inout) :: instance
    real(dp),dimension(:,:,:),intent(in) :: qsmb_g       ! Surface mass balance (m)
    real(dp),dimension(:,:,:),intent(in) :: tsfc_g       ! Surface temperature (C)
    real(dp),dimension(:,:,:),intent(in) :: topo_g       ! Surface elevation (m)
    character (len=12), intent(in)  :: cseg1 
    integer ,dimension(:,:),  intent(in),optional :: gmask ! = 1 where global data are valid
                                                           ! = 0 elsewhere
    integer ,dimension(:,:,:,:),intent(in),optional :: gmask3D
    real(dp),dimension(:,:,:),  intent(in),optional  :: areas_g !GCM box areas, for conservation
    real(dp),dimension(:,:,:),  intent(in),optional  :: sdep_g !GCM nonice-tile snowdepth
    real(dp),dimension(:,:),    intent(out),optional :: sdep_l !icesheet nonice-tile snowdepth

    real(dp),dimension(:), intent(in), optional :: qsmb_lev !seg 31/aug/17
    real(dp), parameter :: maskval = 0.0_dp    ! value written to masked out gridcells
    real(dp)            :: temp_maskval ! seg 31/aug/17
    integer ::       &
       nec,          &      ! number of elevation classes
       nxl, nyl             ! local grid dimensions

    integer :: i, j, n, ig, jg
 
    real(dp), dimension(:,:,:), allocatable ::   &
       qsmb_l,    &! interpolation of global mass balance to local grid
       tsfc_l,    &! interpolation of global sfc temperature to local grid
       topo_l,    &! interpolation of global topography in each elev class to local grid
       sdep_lE     ! interpolation of global nonice-tile snowdepth         
!+seg
    real(dp), dimension (:,:,:), allocatable ::   &
       r_imask_l, &
       r_nmask_l, &
       r_lmask_l
     
    integer(sp), dimension (:,:,:), allocatable ::   &
       i_imask_l, &
       i_nmask_l, &
       i_lmask_l

    real(dp), dimension (:,:), allocatable :: dummy_mask

    real, dimension(:,:), allocatable ::  thck_temp

    real(dp) :: fact, usrf

    real(dp) :: smb_total,smb_col_total,smb_interp_total,smb_adjust,gnumloc
    real(dp) :: dxl,dyl
    integer :: nxg, nyg, il, jl
!+seg
    integer :: iice,inice,ilnd
    integer :: htype
    integer,parameter :: imsk=1
    integer,parameter :: nimsk=2
    integer,parameter :: lmsk=3  ! This was set at 2, bugfix 27Sep17

    real(dp),dimension(:,:),allocatable :: ddata_seg  !30/aug/17
    integer ,dimension(:,:),allocatable :: dmask_seg  !30/aug/17
    integer ,dimension(:,:),allocatable :: dmaskni_seg  !18/Oct/17
    integer ,dimension(:,:),allocatable :: dmasktm_seg  !18/Oct/17
    logical                             :: do_halo = .true. !30/aug/17

    real(dp),dimension(:,:,:),allocatable :: qsmb_g_temp
    real(dp),dimension(:,:,:),allocatable :: sdep_g_temp
 !   integer ,dimension(:,:,:,:),allocatable :: gmask3D_temp
 ! changed back to 3D 29Sep17
    integer ,dimension(:,:,:),allocatable :: gmask3D_temp
!-seg
    real(dp), parameter :: lapse = 0.0065_dp   ! atm lapse rate, deg/m
                                               ! used only for extrapolating temperature outside
                                               !  the range provided by the climate model
    nxg = size(qsmb_g,1)
    nyg = size(qsmb_g,2)
    nec = size(qsmb_g,3)
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    allocate(qsmb_l(nxl,nyl,nec))
    allocate(tsfc_l(nxl,nyl,nec))
    allocate(topo_l(nxl,nyl,nec))
    allocate(thck_temp(nxl,nyl))
    if (present(sdep_g)) allocate(sdep_lE(nxl,nyl,nec))
    allocate(r_imask_l(nxl,nyl,nec))
    allocate(r_nmask_l(nxl,nyl,nec))
    allocate(r_lmask_l(nxl,nyl,nec))
    allocate(i_imask_l(nxl,nyl,nec))
    allocate(i_nmask_l(nxl,nyl,nec))
    allocate(i_lmask_l(nxl,nyl,nec))
    allocate(dummy_mask(nxl,nyl))

    do_halo = .true. ! To make sure
    if (do_halo) then  !30/aug/17
      allocate(ddata_seg(nxg,nyg))
      allocate(dmask_seg(nxg,nyg))
      allocate(dmaskni_seg(nxg,nyg))
      allocate(dmasktm_seg(nxg,nyg))
      allocate(qsmb_g_temp(nxg,nyg,nec))
!+seg 21Oct      
      allocate(sdep_g_temp(nxg,nyg,nec))
!      allocate(gmask3D_temp(3,nxg,nyg,nec))
! Changed back to 3D 29Sep17
      allocate(gmask3D_temp(nxg,nyg,nec))

      gmask3D_temp=0
      write(stdout,*) 'glint_timestep: doing SEG HALO v2'
    else
      write(stdout,*) 'glint_timestep: SEG HALO turned off'
    endif

    if (present(sdep_g)) &
       write(stdout,*) 'SDEP_G in the house'
    !   Downscale global fields for each elevation class to local grid.

    if (present(gmask3D)) then   ! set local field = maskval where the global field is masked out
                               ! (i.e., where instance%downs%lmask = 0).
    ! SEG New option 27/2/17
    ! SEG increased dimensions of gmask3D 28/2/17
    ! Indexing the mask number:should they be params? msk_ice = 1, msk_nice = 2 ?
       do n = 1, nec

          if(do_halo) then
            write(stdout,*) 'glint_timestep: do_halo',imsk
            call create_halo_seg2(qsmb_g(:,:,n),gmask3D(imsk,:,:,n),ddata_seg,dmask_seg,4)
            qsmb_g_temp(:,:,n)=ddata_seg
            !gmask3D_temp(:,:,n)=dmask_seg
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=qsmb_l(:,:,n), &
                               gmask = dmask_seg, &
                               maskval=maskval,seg_verb=1)

          else 
            if (present(qsmb_lev)) then
              temp_maskval=qsmb_lev(n)
              write(stdout,*) 'glint_timestep: level, smb_infill ',n,temp_maskval
            else
              temp_maskval=maskval
            endif

            call interp_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n), &
                               gmask = gmask3D(1,:,:,n), &
                               maskval=temp_maskval,seg_verb=1)
          endif
          write(stdout,*) 'SMB Interpolation gm with 3D mask',n
          write(stdout,*) 'Input max',maxval(qsmb_g(:,:,n))
          write(stdout,*) 'Output max',maxval(qsmb_l(:,:,n))

          if(do_halo) then
            write(stdout,*) 'glint_timestep: do_halo temp', &
              maxval(tsfc_g(:,:,n)),minval(tsfc_g(:,:,n))
            call create_halo_seg2(tsfc_g(:,:,n),gmask3D(imsk,:,:,n),ddata_seg,dmasktm_seg,4) ! Doing temp on imsk
                                                                                             ! otherwise interpolate
                                                                                             ! to non-ice (colder
                                                                                             ! regions)
            ! dmasktm_seg should be the same as dmask_seg, as we are working on
            ! same input mask
            write(stdout,*) 'glint_timestep: do_halo temp2', &
              maxval(ddata_seg),minval(ddata_seg)
            !gmask3D_temp(lmsk,:,:,n)=dmask_seg  !  This should be same at calc in smb step
            
            
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = dmasktm_seg, &
            !                   gmask = gmask3D_temp(:,:,n), &
                               maskval=maskval,seg_verb=10)
            !                   gmask = gmask3D_temp(lmsk,:,:,n), &
            !                   maskval=maskval,seg_verb=10)
            write(stdout,*) 'glint_timestep: do_halo temp3', &
              maxval(tsfc_l(:,:,n)),minval(tsfc_l(:,:,n)),lmsk

          else  
            call interp_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = gmask3D(lmsk,:,:,n), &
                               maskval=maskval,seg_verb=2)
          endif
          ! Not masking this, as global field
          call interp_to_local(instance%lgrid_fulldomain, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))

          if (present(sdep_g)) then
            if(do_halo) then
              call create_halo_seg2(sdep_g(:,:,n),gmask3D(nimsk,:,:,n),ddata_seg,dmaskni_seg,5)
              sdep_g_temp(:,:,n)=ddata_seg     ! This is no longer qsmb
              gmask3D_temp(:,:,n)=dmaskni_seg

              call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, &
                               localdp=sdep_lE(:,:,n), &
                               gmask = dmaskni_seg, &
                               maskval=maskval,seg_verb=15)

            else  

              call interp_to_local(instance%lgrid_fulldomain, sdep_g(:,:,n), instance%downs, &
                               localdp=sdep_lE(:,:,n), &
                               gmask = gmask3D(nimsk,:,:,n), &
                               maskval=maskval,seg_verb=3)
            endif
          endif
          ! Need to also interpolate the mask onto a local mask

          ! Ice mask
          !dummy_mask=0.0
          ! Convert input mask from integer to real
          !ddata_seg=gmask3D_temp(:,:,n)
          ddata_seg=0.0
          ddata_seg=dmask_seg

          if(do_halo) then
            write(stdout,*) 'Smearing a mask',maxval(gmask3D_temp(:,:,n)),&
               minval(gmask3D_temp(:,:,n)),imsk,n
            write(stdout,*) 'sizes 1',size(ddata_seg),size(ddata_seg,1),& 
               size(ddata_seg,2)
            write(stdout,*) 'sizes 2',size(gmask3D_temp(:,:,n)), &
               size(gmask3D_temp(:,:,n),1),&
               size(gmask3D_temp(:,:,n),2),&
               imsk,n
            ! Failing here?! 
            !dummy_mask=gmask3D_temp(:,:,n)
            write(stdout,*) 'Done this'
            ! Using the 'temporary' smeared mask   15Sep17 SEG
            !call interp_to_local(instance%lgrid_fulldomain, dummy_mask, instance%downs, localdp=r_imask_l(:,:,n), &
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=r_imask_l(:,:,n), &
                 gmask = dmask_seg, maskval=maskval,  seg_verb=7)

!                 gmask = gmask3D_temp(:,:,n), maskval=maskval, &
!                 seg_verb=7)
          else  
            !dummy_mask=gmask3D(imsk,:,:,n)
            !call interp_to_local(instance%lgrid_fulldomain, dummy_mask, instance%downs, localdp=r_imask_l(:,:,n), &
            ddata_seg=gmask3D(imsk,:,:,n)
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=r_imask_l(:,:,n), &
                               gmask = gmask3D(imsk,:,:,n), maskval=maskval, &
                               seg_verb=4)
          endif


          ! Non-Ice mask
          !dummy_mask=0.0
          !dummy_mask=gmask3D(nimsk,:,:,n)
          !call interp_to_local(instance%lgrid_fulldomain, dummy_mask, instance%downs, localdp=r_nmask_l(:,:,n), &
          ddata_seg=0.0
          if(do_halo) then
            ddata_seg=dmaskni_seg
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=r_nmask_l(:,:,n), &
                               gmask = dmaskni_seg, maskval=maskval, &
                               seg_verb=5 )


          else  
          ! Convert input mask from integer to real
            ddata_seg=gmask3D(nimsk,:,:,n)
            call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=r_nmask_l(:,:,n), &
                               gmask = gmask3D(nimsk,:,:,n), maskval=maskval, &
                               seg_verb=5 )

          endif
          ! All mask
!          dummy_mask=0.0
!          if(do_halo) then
!          write(stdout,*) 'Smearing an other  mask'
!            dummy_mask=gmask3D_temp(lmsk,:,:,n)
!            call interp_to_local(instance%lgrid_fulldomain, dummy_mask, instance%downs, localdp=r_lmask_l(:,:,n), &
!                               gmask = gmask3D_temp(lmsk,:,:,n), &
!                               maskval=maskval,seg_verb=12)
!          else
!          dummy_mask=gmask3D(lmsk,:,:,n)
!          call interp_to_local(instance%lgrid_fulldomain, dummy_mask, instance%downs, localdp=r_lmask_l(:,:,n), &
          ddata_seg=0.0
          ! Convert input mask from integer to real
          ddata_seg=gmask3D(lmsk,:,:,n)
          call interp_to_local(instance%lgrid_fulldomain, ddata_seg, instance%downs, localdp=r_lmask_l(:,:,n), &
                               gmask = gmask3D(lmsk,:,:,n), &
                               maskval=maskval,seg_verb=6)
!          endif

       enddo
       ! Convert to a 0/1 integer mask
       i_imask_l=0
       i_nmask_l=0
       i_lmask_l=0
       where(r_imask_l > 1e-4) i_imask_l=1
       where(r_nmask_l > 1e-4) i_nmask_l=1
       where(r_lmask_l > 1e-4) i_lmask_l=1

    else if (present(gmask)) then   ! set local field = maskval where the global field is masked out
                               ! (i.e., where instance%downs%lmask = 0).
       do n = 1, nec
          call interp_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          write(stdout,*) 'SMB Interpolation gm ',n
          write(stdout,*) 'Input max',maxval(qsmb_g(:,:,n))
          write(stdout,*) 'Output max',maxval(qsmb_l(:,:,n))

          call interp_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          ! Not masking this, as global field
          call interp_to_local(instance%lgrid_fulldomain, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))
          if (present(sdep_g)) &
          call interp_to_local(instance%lgrid_fulldomain, sdep_g(:,:,n), instance%downs, localdp=sdep_lE(:,:,n), &
                               gmask = gmask, maskval=maskval)
       enddo

    else    ! global field values are assumed to be valid everywhere

       do n = 1, nec
          call interp_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n))
          write(stdout,*) 'SMB Interpolation ',n
          write(stdout,*) 'Input max',maxval(qsmb_g(:,:,n))
          write(stdout,*) 'Output max',maxval(qsmb_l(:,:,n))
          call interp_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n))
          call interp_to_local(instance%lgrid_fulldomain, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))
          if (present(sdep_g)) &
          call interp_to_local(instance%lgrid_fulldomain, sdep_g(:,:,n), instance%downs, localdp=sdep_lE(:,:,n))
       enddo

    endif

    ! The following output only works correctly if running with a single task
    if (GLC_DEBUG .and. tasks==1) then
       ig = iglint_global   ! in glint_type; make sure values are appropriate
       jg = jglint_global
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolate fields to local grid'
       write (stdout,*) 'Global cell =', ig, jg
       do n = 1, nec
          write(stdout,*) n, topo_g(ig,jg, n), qsmb_g(ig,jg, n)
       enddo

       do j = 1, nyl
       do i = 1, nxl
           if ( (instance%downs%xloc(i,j,1) == ig .and. instance%downs%yloc(i,j,1) == jg) .or.  &
                (instance%downs%xloc(i,j,2) == ig .and. instance%downs%yloc(i,j,2) == jg) .or.  &
                (instance%downs%xloc(i,j,3) == ig .and. instance%downs%yloc(i,j,3) == jg) .or.  &
                (instance%downs%xloc(i,j,4) == ig .and. instance%downs%yloc(i,j,4) == jg) ) then
               write(stdout,*) i, j, thk0 * instance%model%geometry%usrf(i,j)
           endif
       enddo
       enddo
    
       i = instance%model%numerics%idiag_global
       j = instance%model%numerics%jdiag_global
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolated to local cells: i, j =', i, j
       do n = 1, nec
          write (stdout,*) ' '
          write (stdout,*) 'n =', n
          write (stdout,*) 'qsmb_l =', qsmb_l(i,j,n)
          write (stdout,*) 'tsfc_l =', tsfc_l(i,j,n)
          write (stdout,*) 'topo_l =', topo_l(i,j,n)
       enddo

    end if ! GLC_DEBUG

!   Interpolate tsfc and qsmb to local topography using values in the neighboring 
!    elevation classes.
!   If the local topography is outside the bounds of the global elevations classes,
!    extrapolate the temperature using the prescribed lapse rate.

    !+ SEG
    ! Writie out instance%model%geometry%usrf here
    !- SEG

    !+SEG redo vertical interp

    !write(stdout,*) 'Starting v interp'

    do j = 1, nyl
    do i = 1, nxl

       usrf = instance%model%geometry%usrf(i,j) * thk0   ! actual sfc elevation (m)

       if (present(gmask3D)) then
         ! Need to find the highest point with valid fraction in vertical
         ! SEG 6/17
         ! i_imask_l & i_lmask_l make have been smeared 27Sep17
         ! NOT i_lmask_l 29Sep17
         ! Interp temp using ice mask
         call hval_index(i_imask_l(i,j,:),iice)
         call hval_index(i_nmask_l(i,j,:),inice)
         call hval_index(i_lmask_l(i,j,:),ilnd)

       endif
       !write(stdout,*) ' After hval_index',i,j
       if (usrf <= topo_l(i,j,1)) then
          write(stdout,*) ' Before p1',i,j
          if (present(sdep_g)) &
           sdep_l(i,j)=sdep_lE(i,j,1)
          !write(stdout,*) ' Before p1.1',i,j
          instance%acab(i,j) = qsmb_l(i,j,1)
          instance%artm(i,j) = tsfc_l(i,j,1) + lapse*(topo_l(i,j,1)-usrf)
! The following works because topo_l interpolated from input gcm 'faketopo'
! i.e. mid-point values _everywhere_ 
! SEG 6/17
         !write(stdout,*) ' After p1',i,j
       elseif (usrf > topo_l(i,j,nec)) then
          !write(stdout,*) ' Before p2',i,j
          if (present(gmask3D)) then 
            ! assign nearest valid data
            instance%acab(i,j) = qsmb_l(i,j,iice)
            if (present(sdep_g)) &
             sdep_l(i,j)=sdep_lE(i,j,inice)
            instance%artm(i,j) = tsfc_l(i,j,ilnd) -  &
                lapse*(usrf-topo_l(i,j,ilnd))
            !write(stdout,*) ' After p2',i,j
          else
            instance%acab(i,j) = qsmb_l(i,j,nec)
            if (present(sdep_g)) &
             sdep_l(i,j)=sdep_lE(i,j,nec)
            instance%artm(i,j) = tsfc_l(i,j,nec) - lapse*(usrf-topo_l(i,j,nec))
            !write(stdout,*) ' After p3',i,j
          endif
       else
          !write(stdout,*) ' In the chufter',i,j
          do n = 2, nec
             if (usrf > topo_l(i,j,n-1) .and. usrf <= topo_l(i,j,n)) then
                fact = (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
                ! DIAG
                if (fact > 1) then
                  write(stdout,*) 'We have big fact'
                endif
! SEG 6/17
! Only interpolate if to valid points, otherwise use nearest neighbour
! 
! hval_type return
!    n   n-1   htype
!    1    1    1
!    1    0    2
!    0    1    3
!    0    0    0
! SEG 1/9/17
! Use level averages if no value for interpolation
! Otherwise we get jumps (if using qsmb_lev)

                if (present(gmask3D)) then
                  ! Changed from n < iice 1Oct17
                  if (n <= iice) then
                     ! ICE & TEMP
                     !write(stdout,*) 'hval + ',i,j
                     call hval_type(i_imask_l(i,j,:),n,htype)
                     !write(stdout,*) 'hval - ',i,j
                             
                     if ((htype == 1) .or. (htype == 0)) then 
                       ! DIAG
                       !if (htype == 0) then
                       !  write(stdout,*) 'We have an htype of 0',i,j,n
                       !endif
                       ! if htype is 0, then there should be nothing to interpolate
                       instance%acab(i,j) = fact*qsmb_l(i,j,n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                       instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                     else if (htype == 2) then
                       if (present(qsmb_lev)) then !1/Sep/17
                         instance%acab(i,j) = fact*qsmb_lev(n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                       else  
                         instance%acab(i,j) = qsmb_l(i,j,n) ! Do not interpolate if point below has
                                                            ! no fraction
                       endif                                     
                       instance%artm(i,j) = tsfc_l(i,j,n)
                       ! DIAG
                       !write(stdout,*) 'We have an htype of 2',i,j,n

                     else ! htype == 3
                       if (present(qsmb_lev)) then !1/Sep/17
                         instance%acab(i,j) = fact*qsmb_lev(n-1) + (1._dp-fact)*qsmb_lev(n)
                       else
                         instance%acab(i,j) = qsmb_l(i,j,n-1)
                       endif
                       instance%artm(i,j) = tsfc_l(i,j,n-1)
                     endif
                  else 
                     !write(stdout,*) 'Setting def at ',i,j,iice
                     instance%acab(i,j) = qsmb_l(i,j,iice)
                     instance%artm(i,j) = tsfc_l(i,j,iice)
                  endif

!                  if (n < ilnd) then
!                     ! TEMP (LAND)
!                     call hval_type(i_lmask_l(i,j,:),n,htype)
!                     if ((htype == 1) .or. (htype == 0)) then
!                       instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
!                     else if (htype == 2) then
!                       instance%artm(i,j) = tsfc_l(i,j,n)
!                     else
!                       instance%artm(i,j) = tsfc_l(i,j,n-1)
!                     endif
!                  else
!                     instance%artm(i,j) = tsfc_l(i,j,ilnd)
!                  endif   

                  !write(stdout,*) ' After p4',i,j

                  ! Changed from n < inice 1Oct17
                  if (present(sdep_g)) then 
                   if (n <= inice) then
                     ! snowcover (non-ice)
                     call hval_type(i_nmask_l(i,j,:),n,htype)
                     if ((htype == 1) .or. (htype == 0)) then
                       sdep_l(i,j)  = fact*sdep_lE(i,j,n-1)+ (1._dp-fact)*sdep_lE(i,j,n)
                     else if  (htype == 2) then
                       sdep_l(i,j)  = sdep_lE(i,j,n)
                     else
                       sdep_l(i,j)  = sdep_lE(i,j,n-1)
                     endif
                   else
                     sdep_l(i,j)  = sdep_lE(i,j,inice)
                   endif   
                  endif 
                  !write(stdout,*) ' After p5',i,j
                else
                  ! Non gmask3D case
                  instance%acab(i,j) = fact*qsmb_l(i,j,n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                  instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                  if (present(sdep_g)) &
                   sdep_l(i,j)        = fact*sdep_lE(i,j,n-1)+ (1._dp-fact)*sdep_lE(i,j,n)
                endif
                exit
             endif
          enddo
       endif   ! usrf

       ! SEG 18/Oct/17
       ! Get rid of sdep_l in areas of ice or ocean

       ! The following output only works correctly if running with a single task
       if (GLC_DEBUG .and. tasks==1) then
          if (i==instance%model%numerics%idiag_global .and. j==instance%model%numerics%jdiag_global) then
             n = 4  
             write (stdout,*) ' '
             write (stdout,*) 'Interpolated values, i, j, n =', i, j, n
             write (stdout,*) 'usrf =', usrf
             write (stdout,*) 'acab =', instance%acab(i,j)
             write (stdout,*) 'artm =', instance%artm(i,j)
             write (stdout,*) 'topo(n-1) =', topo_l(i,j,n-1)
             write (stdout,*) 'topo(n) =', topo_l(i,j,n)
             write (stdout,*) 'qsmb(n-1) =', qsmb_l(i,j,n-1)
             write (stdout,*) 'qsmb(n) =', qsmb_l(i,j,n)
             write (stdout,*) 'tsfc(n-1) =', tsfc_l(i,j,n-1)
             write (stdout,*) 'tsfc(n) =', tsfc_l(i,j,n)
             write (stdout,*) 'fact = ', (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
          endif
         !write(stdout,*) ' After p6',i,j
       end if

    enddo  ! i
    enddo  ! j

    !write(stdout,*) 'Ending v interp'
    write(stdout,*) 'ACAB Check 1 (max/min)', & 
      maxval(instance%acab),minval(instance%acab)

! Surface fields
!    write(stdout,*) 'Surface fields sizes '
!    write(stdout,*) 'acab ',size(instance%acab,1),size(instance%acab,2)
!    write(stdout,*) 'artm ',size(instance%artm,1),size(instance%artm,2)
!    write(stdout,*) 'sdep_l ',size(instance%sdep_l,1),size(instance%sdep_l,2)

    if (present(areas_g)) then
    write(6,*)"adjust SMB to conserve wrt each GCM gridbox. NASTY! But ..."
    write(6,*)"We are not doing it here. So there!"

!    !!apply mask as glint_timestep will later do
!    call glide_get_thk(instance%model,thck_temp)
!    where(thck_temp<=0.0 .and. instance%acab<0.0)
!      instance%acab = 0.0
!    end where
!    where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
!       instance%acab = 0.0
!    endwhere
!
!    smb_total=0
!    !NH only
!    do jg = 1, nyg/2
!    do ig = 1, nxg
!      smb_col_total=0.
!      smb_interp_total=0.
!      gnumloc=0.
!      dxl=instance%lgrid%delta%pt(1)
!      dyl=instance%lgrid%delta%pt(2)
!      !fractional areas are all in areas_g, simple integral for GCM NEC column total
!      do n = 1, nec
!         smb_col_total=smb_col_total+areas_g(ig,jg,n)*qsmb_g(ig,jg,n)
!         if (isnan(smb_col_total)) then
!           write(stdout,*) '+SMBC ISNAN',ig,jg,n
!           write(stdout,*) areas_g(ig,jg,n),qsmb_g(ig,jg,n)
!           write(stdout,*) '-SMBC ISNAN'
!         endif
!      end do
!      smb_total=smb_total+smb_col_total
!      if (isnan(smb_total) .or. (abs(smb_total)  > 100) ) then
!        write(stdout,*) '+SMB ISNAN',ig,jg
!        write(stdout,*) smb_col_total, smb_total
!        write(stdout,*) smb_col_total,areas_g(ig,jg,n),qsmb_g(ig,jg,n)
!        write(stdout,*) '-SMB ISNAN'
!      end if
!      !find all icesheet points that contribute to this GCM column and integrate.
!      do jl=1,nyl
!      do il=1,nxl
!         if (instance%ups%gboxx(il,jl).eq.ig .AND. instance%ups%gboxy(il,jl).eq.jg &
!             .AND.abs(instance%acab(il,jl)).gt.0                                   &
!            ) then
!           smb_interp_total=smb_interp_total+instance%acab(il,jl)*dxl*dyl
!           gnumloc=gnumloc+1.
!         endif
!      end do
!      end do
!      !adjust each non-zero icesheet point so that they match the GCM column integral
!      smb_adjust=smb_col_total-smb_interp_total
!      if (abs(smb_adjust).gt.0) then
!       do jl=1,nyl
!       do il=1,nxl
!          if (instance%ups%gboxx(il,jl).eq.ig .AND. instance%ups%gboxy(il,jl).eq.jg &
!             .AND.abs(instance%acab(il,jl)).gt.0                                   &
!             ) then
!            instance%acab(il,jl)=instance%acab(il,jl)+(smb_adjust/(dxl*dyl)/gnumloc)
!          endif
!       end do
!       end do
!      endif
!
!    end do !GCM gridbox loop
!    end do
!
!    write(stdout,*) 'ACAB Check 2 (max/min)', & 
!      maxval(instance%acab),minval(instance%acab),smb_total
!
!    !The above will miss points that are ice in the GCM but don't exist in the icesheet - 
!    !coastline/boundary mismatches etc can cause this. Do a final, full domain correction
!    smb_interp_total=0
!    do jl = 1, nyl
!    do il = 1, nxl
!       smb_interp_total=smb_interp_total+instance%acab(il,jl)*dxl*dyl
!    end do
!    end do
!    if (abs(smb_interp_total) .gt. 1e-12) then
!      write(stdout,*) 'ACAB mismatch components: '
!      write(stdout,*) 'ACAB-a ',&
!        maxval(instance%acab),smb_total,smb_interp_total
!      instance%acab=instance%acab*(smb_total/smb_interp_total)
!      write(stdout,*) 'ACAB-b ',&
!        maxval(instance%acab),(smb_total/smb_interp_total)
!    endif

    end if! areas_g

    write(stdout,*) 'ACAB Check 3 (max/min)', & 
      maxval(instance%acab),minval(instance%acab)

    !deallocate(qsmb_l, tsfc_l, topo_l, sdep_lE)
    deallocate(qsmb_l, tsfc_l, topo_l)
    if (present(sdep_g)) &
      deallocate(sdep_lE)
    !+seg
    deallocate(r_imask_l,r_nmask_l,r_lmask_l)
    deallocate(i_imask_l,i_nmask_l,i_lmask_l)
    deallocate(dummy_mask)
    if (do_halo) then  !30/aug/17
      deallocate(ddata_seg)
      deallocate(dmask_seg)
      deallocate(dmaskni_seg)
      deallocate(dmasktm_seg)
      deallocate(qsmb_g_temp)
      deallocate(gmask3D_temp)

      deallocate(sdep_g_temp)
    endif

  end subroutine glint_downscaling_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Is there a better place to put this subroutine?

  subroutine get_i_upscaled_fields(instance,                    &
                                   orog,         albedo,        &
                                   ice_frac,     veg_frac,      &
                                   snowice_frac, snowveg_frac,  &
                                   snow_depth)

    !*FD Upscales and returns certain fields
    !*FD Output fields are only valid on the main task
    !*FD 
    !*FD \begin{itemize}
    !*FD \item \texttt{orog} --- the orographic elevation (m)
    !*FD \item \texttt{albedo} --- the albedo of ice/snow (this is only a notional value --- need to do
    !*FD some work here)
    !*FD \item \texttt{ice\_frac} --- The fraction covered by ice
    !*FD \item \texttt{veg\_frac} --- The fraction of exposed vegetation
    !*FD \item \texttt{snowice\_frac} --- The fraction of snow-covered ice
    !*FD \item \texttt{snowveg\_frac} --- The fraction of snow-covered vegetation
    !*FD \item \texttt{snow_depth} --- The mean snow-depth over those parts covered in snow (m w.e.)
    !*FD \end{itemize}

    use glimmer_paramets

    ! Arguments ----------------------------------------------------------------------------------------

    type(glint_instance),   intent(in)  :: instance      !*FD the model instance

    real(rk),dimension(:,:),intent(out) :: orog          !*FD the orographic elevation (m)
    real(rk),dimension(:,:),intent(out) :: albedo        !*FD the albedo of ice/snow
    real(rk),dimension(:,:),intent(out) :: ice_frac      !*FD The fraction covered by ice
    real(rk),dimension(:,:),intent(out) :: veg_frac      !*FD The fraction of exposed vegetation
    real(rk),dimension(:,:),intent(out) :: snowice_frac  !*FD The fraction of snow-covered ice
    real(rk),dimension(:,:),intent(out) :: snowveg_frac  !*FD The fraction of snow-covered vegetation
    real(rk),dimension(:,:),intent(out) :: snow_depth    !*FD The mean snow-depth over those 
    !*FD parts covered in snow (m w.e.)

    ! Internal variables -------------------------------------------------------------------------------

    real(rk),dimension(:,:),pointer :: temp => null()

    ! --------------------------------------------------------------------------------------------------
    ! Orography

    call mean_to_global(instance%ups_orog, &
                        instance%model%geometry%usrf, &
                        orog,    &
                        instance%out_mask)
    orog=thk0*orog

    call coordsystem_allocate(instance%lgrid,temp)

    !TODO - Change to dp
    ! Ice-no-snow fraction
    where (instance%mbal_accum%snowd==0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        ice_frac,    &
                        instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowice_frac,    &
                        instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck<=(10.0/thk0))
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowveg_frac,    &
                        instance%out_mask)

    ! Remainder is veg only
    veg_frac=1.0-ice_frac-snowice_frac-snowveg_frac

    ! Snow depth

    call mean_to_global(instance%ups, &
                        instance%mbal_accum%snowd, &
                        snow_depth,    &
                        instance%out_mask)

    ! Albedo

    where ((ice_frac+snowice_frac)>0.0)
       albedo=instance%ice_albedo
    elsewhere
       albedo=0.0
    endwhere

    deallocate(temp)
    temp => null()

  end subroutine get_i_upscaled_fields

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Is there a better place to put this subroutine?

  subroutine get_i_upscaled_fields_gcm(instance,    nec,      &
                                       nxl,         nyl,      &
                                       nxg,         nyg,      &
                                       gfrac,       gtopo,    &
                                       grofi,       grofl,    &
                                       ghflx,       glfrac,   &
                                       gsdep,       lsdep)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! Output fields are only valid on the main task.
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.

    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glimmer_log
    use parallel, only: tasks, main_task

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance),     intent(in)  :: instance      ! the model instance
    integer,                  intent(in)  :: nec           ! number of elevation classes
    integer,                  intent(in)  :: nxl,nyl       ! local grid dimensions 
    integer,                  intent(in)  :: nxg,nyg       ! global grid dimensions 

    !TODO - Should these be inout?
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
    real(dp),dimension(nxg,nyg,nec),intent(out),optional :: glfrac   ! fraction of gridbox area in each mec
    real(dp),dimension(nxg,nyg,nec),intent(out),optional :: gsdep   ! anomalies to nonice-snow field on GCM grid
    real(dp),dimension(nxl,nyl),intent(in),optional :: lsdep   ! anomalies to nonice-snow field on icesheet

 
    ! Internal variables ----------------------------------------------------------------------
 
    real(dp),dimension(nxl,nyl) :: local_field
    real(dp),dimension(nxl,nyl) :: local_topo   ! local surface elevation (m)
    real(dp),dimension(nxl,nyl) :: local_thck   ! local ice thickness (m)

    !TODO - Put this parameter elsewhere?  Make it equal to thkmin?
    real(dp), parameter :: min_thck = 1.0_dp    ! min thickness (m) for setting gfrac = 1

    integer :: i, j            ! indices
 
    integer :: il, jl, ig, jg
    character(len=100) :: message

    !TODO - Pass in topomax as an argument instead of hardwiring it here
    real(dp), dimension(0:nec) :: topomax   ! upper elevation limit of each class

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! Note: These must be consistent with the values in the GCM.  Better to pass as an argument.
    if (nec == 1) then
       topomax = (/ 0._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 3) then
       topomax = (/ 0._dp,  1000._dp,  2000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 5) then
       topomax = (/ 0._dp,   500._dp,  1000._dp,  1500._dp,  2000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 10) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   700._dp,  1000._dp,  1300._dp,  &
                            1600._dp,  2000._dp,  2500._dp,  3000._dp, 10000._dp /)
    elseif (nec == 25) then
       topomax = (/ 0._dp,   200._dp,   400._dp, 600._dp,  800._dp,  1000._dp,  &
                            1200._dp,  1400._dp, 1600._dp, 1800._dp,  2000._dp, &
                            2200._dp,  2400._dp, 2600._dp, 2800._dp,  3000._dp, &
                            3200._dp,  3400._dp, 3600._dp, 3800._dp,  4000._dp, &
                            4200._dp,  4400._dp, 4600._dp, 4800._dp,  5000._dp /)
    elseif (nec == 36) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   600._dp,   800._dp,  &
                 1000._dp,  1200._dp,  1400._dp,  1600._dp,  1800._dp,  &
                 2000._dp,  2200._dp,  2400._dp,  2600._dp,  2800._dp,  &
                 3000._dp,  3200._dp,  3400._dp,  3600._dp,  3800._dp,  &
                 4000._dp,  4200._dp,  4400._dp,  4600._dp,  4800._dp,  &
                 5000._dp,  5200._dp,  5400._dp,  5600._dp,  5800._dp,  &
                 6000._dp,  6200._dp,  6400._dp,  6600._dp,  6800._dp,  &
                 7000._dp, 10000._dp /)
    else
       if (GLC_DEBUG .and. main_task) then
          write(message,'(a6,i3)') 'nec =', nec
          call write_log(trim(message), GM_DIAGNOSTIC)
       end if
       call write_log('ERROR: Current supported values of nec (no. of elevation classes) are 1, 3, 5, 10, 25 or 36', &
                       GM_FATAL,__FILE__,__LINE__)
    endif

    local_topo(:,:) = thk0 * instance%model%geometry%usrf(:,:)
    local_thck(:,:) = thk0 * instance%model%geometry%thck(:,:)

    ! The following output only works correctly if running with a single task
    if (GLC_DEBUG .and. tasks==1) then
       ig = iglint_global    ! defined in glint_type
       jg = jglint_global
       il = instance%model%numerics%idiag_global
       jl = instance%model%numerics%jdiag_global
       write(stdout,*) 'In get_i_upscaled_fields_gcm'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'topo =', local_topo(il,jl) 
       write(stdout,*) 'thck =', local_thck(il,jl) 
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)
    end if

    write(stdout,*) ' RSS  range(local_topo) = ', minval(local_topo) , maxval(local_topo)
    write(stdout,*) ' RSS  range(model%geometry%usrf) = ', minval(instance%model%geometry%usrf) , maxval(instance%model%geometry%usrf)
    write(stdout,*) ' RSS  range(model%geometry%thck) = ', minval(instance%model%geometry%thck) , maxval(instance%model%geometry%thck)

    ! temporary field: = 1 where ice thickness exceeds threshold, else = 0
    !TODO - Use > 0 as threshold?

    do j = 1, nyl
    do i = 1, nxl
       if (local_thck(i,j) > min_thck) then
          local_field(i,j) = 1.d0
       else
          local_field(i,j) = 0.d0
       endif
    enddo
    enddo
    write(stdout,*) 'After temp 1' ! seg
    ! ice fraction
    !TODO - gfrac should be fraction of total grid cell with ice in each elevation class
    !       Currently is ice-covered fraction of cells in a given elevation class

    call mean_to_global_mec(instance%ups,                       &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            local_field,        gfrac,          &
                            local_topo,         instance%out_mask)

    write(stdout,*) 'After temp 2' ! seg
    if (present(glfrac)) then
      ! use full grid. Land-sea mask taken care of by out_mask, here
      ! and in splice
      local_field(:,:) = 1.d0

      call mean_to_global_mec(instance%ups,                     &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            local_field,        glfrac,          &
                            local_topo,         instance%out_mask, &
                            colfrac=.TRUE.)
    end if

    !write(stdout,*) 'After temp 3' ! seg
    if (present(gsdep)) then
      !write(stdout,*) 'Shouldnt be in here'
      call mean_to_global_mec(instance%ups,                     &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            lsdep,              gsdep,          &
                            local_topo,         instance%out_mask)

      write(6,*)sum(lsdep),sum(gsdep)
    end if


    ! surface elevation

    !write(stdout,*) 'Before temp 4' ! seg
    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_topo,          gtopo,     &
                            local_topo,          instance%out_mask)

    !write(stdout,*) 'After temp 4' ! seg
    !TODO - For upscaling, need to copy the appropriate Glide fields into the local_field array

    ! ice runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofi,     &
                            local_topo,          instance%out_mask)

    ! liquid runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofl,     &
                            local_topo,          instance%out_mask)

    ! heat flux

!    local_field(:,:) = 0._dp
!JUST GET INSTANTANEOUS FOR NOW? rssmith
     local_field(:,:) = instance%model%temper%ucondflx(:,:)


    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         ghflx,     &
                            local_topo,          instance%out_mask)
    
    if (GLC_DEBUG .and. main_task) then

!       write(stdout,*) ' '
!       write(stdout,*) 'global ifrac:'
!       do n = 1, nec
!          write(stdout,*) n, gfrac(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gtopo:'
!       do n = 1, nec
!          write(stdout,*) n, gtopo(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofi:'
!       do n = 1, nec
!          write(stdout,*) n, grofi(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofl:'
!       do n = 1, nec
!          write(stdout,*) n, grofl(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global ghflx:'
!       do n = 1, nec
!          write(stdout,*) n, ghflx(ig, jg, n)
!       enddo

    end if

  end subroutine get_i_upscaled_fields_gcm

  subroutine hval_index(inarr,outval)
  !
  ! SEG 28/6/17

  implicit none

  integer,intent(in)  :: inarr(:)
  integer,intent(out) :: outval

  ! local

  integer :: nvals,j

  nvals=size(inarr,1)

  hval1: do j=nvals,1,-1   ! Get peak value
    if (inarr(j) == 1) then
      outval=j
      exit
    endif

  end do hval1

  end subroutine hval_index

  subroutine hval_type(inarr,ipos,htype)

  !
  ! SEG 28/6/17
  implicit none

  integer,intent(in)   :: inarr(:)
  integer,intent(in)   :: ipos
  integer,intent(out)  :: htype

  ! local
  integer      :: iposd
  integer      :: nsize

  nsize=size(inarr,1)
  iposd = ipos - 1 
  !if ((iposd < 2)  .or. (ipos > nsize)) then
  !    call write_log('ERROR: In/out array indexe in hval_type',&
  !                                GM_FATAL,__FILE__,__LINE__)
  !endif

  htype = 0
  if ((inarr(ipos) == 1) .and. (inarr(iposd) == 1)) then
    htype = 1
  else if ((inarr(ipos) == 1) .and. (inarr(iposd) == 0)) then
    htype = 2
  else if ((inarr(ipos) == 0) .and. (inarr(iposd) == 1)) then
    htype = 3
  else
    htype = 0
  endif

  end subroutine hval_type

  subroutine create_halo_seg2(data_in,mask_in,data_out,mask_out,niter_override)
  !SEG 20/sep/17
  ! This version should create an average value from nearby points
  ! 
  implicit none

  real(dp), dimension(:,:), intent(in)  :: data_in
  integer ,dimension(:,:), intent(in)   :: mask_in
  real(dp), dimension(:,:), intent(out) :: data_out
  integer ,dimension(:,:), intent(out)  :: mask_out
  integer, intent(in), optional         :: niter_override


  ! local
  real(dp),dimension(:,:),allocatable   :: data_work
  integer ,dimension(:,:),allocatable  :: mask_work
  integer     :: niter
  integer     :: nx,ny
  integer     :: i,j,ii,jj,tt
  integer     :: xneg,xpos,yneg,ypos
  integer     :: mcount
  real(dp)    :: accum

  ! Decide how many time to iterate
  ! Maybe should check for valid number?
  if (present(niter_override)) then
    niter = niter_override
  else 
    niter = 2
  endif

  ! Data sizes
  nx=size(data_in,1)
  ny=size(data_in,2)
  write(stdout,*) 'Halo v2 check ',nx,ny,niter,maxval(data_in)

  allocate(data_work(nx,ny))
  allocate(mask_work(nx,ny))

  data_out=data_in
  mask_out=mask_in

  do tt=1,niter
 
    data_work=data_out
    mask_work=mask_out

    do i=1,nx
      do j=1,ny
        if (mask_out(i,j) < 1) then ! No value here, try to make one
          ! Create box around point
          xneg = max(i-1,1)
          xpos = min(i+1,nx)
          yneg = max(j-1,1)
          ypos = min(j+1,ny)
          accum=0.0
          mcount=0
          do ii = xneg,xpos
            do jj = yneg, ypos
              if (mask_out(ii,jj) > 0) then
                 mcount=mcount+1
                 accum=accum+data_out(ii,jj)
              endif
            end do
          end do 
          if (mcount > 0) then
            ! Create data value, and update mask
            data_work(i,j) = accum/mcount ! average value
            mask_work(i,j) = 1
          end if
        end if
      end do  
    end do    

    data_out=data_work
    mask_out=mask_work

  end do

  deallocate(data_work)
  deallocate(mask_work)

  end subroutine create_halo_seg2


  subroutine create_halo_seg(data_in,mask_in,data_out,mask_out)
  !SEG 30/aug/17
  implicit none

  real(dp), dimension(:,:), intent(in)  :: data_in
  integer ,dimension(:,:), intent(in)   :: mask_in
  real(dp), dimension(:,:), intent(out) :: data_out
  integer ,dimension(:,:), intent(out)  :: mask_out

  ! local
  integer nx,ny
  integer i,j,ii,jj
  integer xneg,xpos,yneg,ypos
  integer hcount

  data_out=data_in
  mask_out=mask_in
  !
  ! Create a halo of replicated data around point

  nx=size(data_in,1)
  ny=size(data_in,2)
  hcount = 0

  do i=1,nx
    do j=1,ny
      ! Could do the following on 'mask_out' but this would
      ! iterate the change. Unstable?
      if (mask_in(i,j) > 0) then ! We have a valid point
        ! Have to consider edge or central points
        xneg = i-1
        xpos = i+1
        yneg = j-1
        ypos = j+1
        ! Go around the halo
        do ii = xneg,xpos
          do jj = yneg, ypos
            if ( (ii .gt. 0) .and. (ii .le. nx) .and. &
                 (jj .gt. 0) .and. (jj .le. ny) ) then
                 if (mask_in(ii,jj) == 0) then
                    mask_out(ii,jj) = 1
                    data_out(ii,jj) = data_in(i,j)
                    hcount = hcount+1
                 endif
            endif
          end do
        end do
      endif  
    end do    
  end do
  write(stdout,*) 'Halo count seg ',hcount

  end subroutine create_halo_seg

  subroutine do_ice_wall(thkmask,usrf,thck,ice_wall)

  ! seg - serial only
  integer,dimension(:,:),intent(in)     :: thkmask
  real(dp),dimension(:,:),intent(in)    :: usrf
  real(sp),dimension(:,:),intent(inout) :: thck
  real(sp),dimension(:,:),intent(out)   :: ice_wall
  
  ! local
  real(sp), parameter  :: hcliff = 100.
  real(sp)             :: ice_thick_asl
  integer nx,ny 
  integer i,j,ii,jj
  integer xneg,xpos,yneg,ypos
  integer hcount
  real    usflocal

  nx=size(thck,1)
  ny=size(thck,2)

  ice_wall = 0.
  do i=1,nx
    do j=1,ny
      if (usrf(i,j) > 0) then ! We are above sealevel
        if ((usrf(i,j)-thck(i,j)) < 0) then ! grounded under sea
          ice_thick_asl = usrf(i,j)
        else
          ice_thick_asl = thck(i,j)
        endif
        if (ice_thick_asl > hcliff) then
          xneg = i-1
          xpos = i+1
          yneg = j-1
          ypos = j+1
          ! Go around the halo
          hcount=0
          do ii = xneg,xpos
            do jj = yneg, ypos
              if ( (ii .gt. 0) .and. (ii .le. nx) .and. &
                   (jj .gt. 0) .and. (jj .le. ny) ) then
                   ! GLIDE_IS_OCEAN - ocean and no ice
                   if (GLIDE_IS_OCEAN(thkmask(ii,jj)) .and. (hcount == 0) .and. &
                      (ii /= i) .and. (jj /= j)) then
                      ! Dump ice into first ocean point
                      ice_wall(ii,jj)=ice_wall(ii,jj)+ ice_thick_asl
                      ! Remove ice cliff
                      thck(i,j) = thck(i,j)-ice_thick_asl
                      hcount = 1
                   end if  
              endif
            end do
          end do
        endif  
      endif
    end do    
  end do
  write(stdout,*) 'Ice wall sum ',sum(ice_wall)
  end subroutine do_ice_wall

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_timestep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


