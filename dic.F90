
!-----------------------------------------------------
!
! University of Victoria DIC/alkalinity model
!
! Contributors: Hakase Hayashida, Nadja Steiner, Eric Mortenson, Adam Mohanan
!
! Model description: The model will be described in a paper.
!

! ericnotes: this module should take arguments meant to affect the dic level
!   at a point in space, this would mean biological take up and loss, as well
!   as air-sea exchange, and sat. constants for the different variables, and
!   alkalinity...
!------------------------------------------------------

! 

#include "fabm_driver.h"

module uvic_dic

   use fabm_types
   use fabm_expressions
   use fabm_standard_variables

   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_dic
      ! Add variable identifiers and parameters here.
! Declare prognostic variables      
      type (type_state_variable_id) :: id_dic, id_eco_dic, id_phys_dic, id_alk
! Declare horizontal prognostic variables      
! Declare diagnostic variables
      type (type_diagnostic_variable_id) :: id_fdic1, id_falk1
      type (type_diagnostic_variable_id) :: id_alk_check, id_dic_check
      type (type_diagnostic_variable_id) :: id_hplus, id_hco3, id_co3, id_co2sw
      type (type_diagnostic_variable_id) :: id_newdic, id_newalk, id_newalk2
      type (type_horizontal_diagnostic_variable_id) :: id_pco2sw, id_pH
      type (type_horizontal_diagnostic_variable_id) :: id_co2flux, id_fIA_co2
      type (type_horizontal_diagnostic_variable_id) :: id_fdic_ice, id_falk_ice
! Declare environmental variables
      type (type_global_dependency_id) :: id_dt
      type (type_dependency_id) :: id_temp, id_sal, id_press
      type (type_dependency_id) :: id_depth
      type (type_horizontal_dependency_id) :: id_ice_hi
      type (type_dependency_id) :: id_stand_dic, id_stand_alk
      type (type_dependency_id) :: id_no3SW, id_fno3phy, id_fde1nh4, id_fde2nh4, id_fde3nh4, id_fnh4phy, id_fnh4no3
      type (type_dependency_id) :: id_Tice

! Declare horizontal environmental variables
      type (type_horizontal_dependency_id) :: id_wind, id_Tatm, id_Patm
      type (type_horizontal_dependency_id) :: id_precip, id_evap
      type (type_horizontal_dependency_id) :: id_botmelt, id_botgrowth, id_topmelt, id_topgrowth, id_termelt
      type (type_horizontal_dependency_id) :: id_fgrowIA, id_fno3upIA, id_fnitIA, id_fCO2
!      type (type_horizontal_dependency_id) :: id_fgrowIA, id_fCO2, id_fno3upIA
      type (type_horizontal_dependency_id) :: id_fno3up
      type (type_horizontal_dependency_id) :: id_stemp, id_surftemp
! Declare namelist parameters
      real(rk) :: dic_0, alk_0, dic_sw, alk_sw, dic_ice, alk_ice, ik_diff, ik_on, ice_on, IA_on, tplv, btlv, prop2sw, prop2sw_melt
! Declare anything else used in all procedures
      real(rk) :: spd = 86400.0_rk ! sec/day
      integer :: icepump, IApump

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface

      ! Reference model procedures here.
   end type

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!     initialize subroutine          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine initialize(self,configunit)
      class (type_uvic_dic), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
! Declare namelist parameters
   
   real(rk) :: dic_0, alk_0, dic_sw, alk_sw, dic_ice, alk_ice, ik_diff, ik_on, ice_on, IA_on, tplv, btlv, prop2sw, prop2sw_melt
   integer  :: icepump, IApump

!jpnote yaml variables 
 
   call self%get_parameter(self%dic_0, 'dic_0', 'initial DIC in water column','mmol/m3', default=2190.0_rk)
   call self%get_parameter(self%alk_0 , 'alk_0 ', 'initial TA in water column','mmol/m3', default=2100.0_rk)
   call self%get_parameter(self%dic_sw, 'dic_sw', 'dic_sw','', default=2100.0_rk)
   call self%get_parameter(self%alk_sw, 'alk_sw', 'alk_sw','', default=2200.0_rk)
   call self%get_parameter(self%dic_ice, 'dic_ice', '[DIC] for ice','mmol/m3', default=400.0_rk)  !400.0_rk 
   call self%get_parameter(self%alk_ice, 'alk_ice', '[TA] for growing ice','mmol/m3', default=500.0_rk)   !500.0_rk
   call self%get_parameter(self%ik_diff, 'ik_diff', 'difference (in melting ice) in [DIC] and 2*[TA] for ice with ikaite precipitaion','', default=50.0_rk)
   call self%get_parameter(self%ik_on, 'ik_on', '(0 or 1), turns off ikaite pump','0 or 1', default=1.0_rk)
   call self%get_parameter(self%ice_on, 'ice_on', '(0 or 1), turns off ice carbon pump','0 or 1', default=1.0_rk)
   call self%get_parameter(self%IA_on, 'IA_on', '(0 or 1), turns off ice algae carbon pump (not used now, bc can do same by ia_0=ia_b=0 in uvic_icealgae)','0 or 1', default=1.0_rk)
   call self%get_parameter(self%tplv, 'tplv', 'top of brine-associated DIC depth','m', default=40.0_rk)
   call self%get_parameter(self%btlv, 'btlv', 'bottom of brine-associated DIC depth','m', default=50.0_rk)
   call self%get_parameter(self%prop2sw, 'prop2sw', 'proportion of DIC rejected that is released into the ocean (the remainder presumably into the atmosphere)','', default=0.99_rk) 
   
   prop2sw_melt=0.975_rk !SA: 1.0, def., 0.95, 0.9, 0.5  !right now, default = 0.975
   self%prop2sw_melt = prop2sw_melt
   !#endif 

! Register prognostic variables 
   call self%register_state_variable(self%id_dic,'dic','mmol DIC/m^3','dissolved inorganic carbon',initial_value=self%dic_0)
   call self%register_state_variable(self%id_alk,'alk','mmol(eq)/m^3','alkalinity',initial_value=self%alk_0,minimum=0.0_rk)
   call self%register_state_variable(self%id_eco_dic,'eco_dic','mmol DIC/m^3','eco DIC')

! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_fdic1,'fdic1','mmol DIC/m**3/day','1st DIC flux')
      call self%register_diagnostic_variable(self%id_falk1,'falk1','mmol(eq)/m**3/day','1st alk flux')
      call self%register_diagnostic_variable(self%id_dic_check,'dic_check','dic_check-units','dic_check_lngnm')
      call self%register_diagnostic_variable(self%id_alk_check,'alk_check','alk_check-units','alk_check_lngnm')
      call self%register_diagnostic_variable(self%id_hplus,'hplus','mmol H+/m**3','H+ concentration in seawater')
      call self%register_diagnostic_variable(self%id_hco3,'hco3','mmol HCO_3/m**3','carbonate concentration in seawater')
      call self%register_diagnostic_variable(self%id_co3,'co3','mmol CO_3/m**3','bicarbonate concentration in seawater')
      call self%register_diagnostic_variable(self%id_co2sw,'co2sw','mmol CO2/m**3','CO_2 concentration in seawater')
      call self%register_diagnostic_variable(self%id_newdic,'newdic','mmol C/m**3','dic concentration in seawater')
      call self%register_diagnostic_variable(self%id_newalk,'newalk','mmol (eq)/m**3','approx. alkalinity in seawater')
      call self%register_diagnostic_variable(self%id_newalk2,'newalk2','mmol (eq)/m**3','tot. alkalinity')

      call self%register_horizontal_diagnostic_variable(self%id_pH,'pH','-log(H)','pH/acidity',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_pco2sw,'pco2sw','micro atm ','partial P of CO2 in seawater',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_co2flux,'co2flux','mmol CO_2/m**2/s','air-sea CO2 flux',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fIA_co2,'fIA_co2','mmol CO_2/m**2/s','ice-algal CO2 flux',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdic_ice,'fdic_ice','mmol CO_2/m**2/s','ice-growth/melt CO2 flux',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_falk_ice,'falk_ice','[mmol-eq]/m**2/s','ice-growth/melt alk flux',source=source_do_horizontal)

! Register environmental variables

      call self%register_global_dependency(self%id_dt,standard_variables%timestep) 
      call self%register_horizontal_dependency(self%id_Tatm,standard_variables%surface_temperature)
      call self%register_horizontal_dependency(self%id_Patm,standard_variables%surface_air_pressure)
      call self%register_horizontal_dependency(self%id_wind,standard_variables%wind_speed)
      call self%register_horizontal_dependency(self%id_ice_hi,standard_variables%sea_ice_thickness)
      call self%register_horizontal_dependency(self%id_surftemp,standard_variables%surface_ice_temperature)
!ericmod 30jan2018 taking out precip as a shared variable. it's not used and causes an error
!      call self%register_horizontal_dependency(self%id_precip,standard_variables%precip)
      call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt)
      call self%register_horizontal_dependency(self%id_topmelt,standard_variables%topmelt)
      call self%register_horizontal_dependency(self%id_topgrowth,standard_variables%topgrowth)
      call self%register_horizontal_dependency(self%id_termelt,standard_variables%termelt)
      call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)

      call self%register_dependency(self%id_depth,standard_variables%depth)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_press,standard_variables%pressure)
      call self%register_dependency(self%id_sal,standard_variables%practical_salinity)

! Register dic from fabm standard variables for comparison of physically driven dic
! ericmod 5dec14, taking out the stand_dic from standard_variables call
! IF YOU WANT TO BRING IT BACK, you need to implement pml in gotm-cases/resolute/fabm.nml
!     call self%register_dependency(self%id_stand_dic, standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
!     call self%register_dependency(self%id_stand_alk, standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_dependency(self%id_stemp,'uvic_eco_stemp','','')
      call self%request_coupling(self%id_stemp,'uvic_eco_stemp')
      call self%register_horizontal_dependency(self%id_fgrowIA,'uvic_icealgae_fgrow','','')
      call self%request_coupling(self%id_fgrowIA,'uvic_icealgae_fgrow')

!ericmod 5sep2017, adding no3 uptake and ice algal nitrification rate for TA calculation
      call self%register_horizontal_dependency(self%id_fno3upIA,'uvic_icealgae_fno3up','','')
      call self%request_coupling(self%id_fno3upIA,'uvic_icealgae_fno3up')

      call self%register_horizontal_dependency(self%id_fnitIA,'uvic_icealgae_fnit','','')
      call self%request_coupling(self%id_fnitIA,'uvic_icealgae_fnit')

      call self%register_dependency(self%id_fno3phy,'uvic_eco_fno3phy','','')
      call self%request_coupling(self%id_fno3phy,'uvic_eco_fno3phy')
      call self%register_dependency(self%id_fde1nh4,'uvic_eco_fde1nh4','','')
      call self%request_coupling(self%id_fde1nh4,'uvic_eco_fde1nh4')
      call self%register_dependency(self%id_fde2nh4,'uvic_eco_fde2nh4','','')
      call self%request_coupling(self%id_fde2nh4,'uvic_eco_fde2nh4')
      call self%register_dependency(self%id_fnh4phy,'uvic_eco_fnh4phy','','')
      call self%request_coupling(self%id_fnh4phy,'uvic_eco_fnh4phy')
      call self%register_dependency(self%id_fnh4no3,'uvic_eco_fnh4no3','','')
      call self%request_coupling(self%id_fnh4no3,'uvic_eco_fnh4no3')

!ericmod 5dec14 , taking out pml calc
!      call self%register_dependency(self%id_stand_alk,'pml_carbonate_alk')
!      call self%request_coupling(self%id_stand_alk,'pml_carbonate_alk')

   end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!          do subroutine                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uvic_dic),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
   INTEGER  :: it
   INTEGER  :: itmax = 200
   real(rk) :: dic, alk, eco_dic   !state variables
   real(rk) :: fdic1, falk1, dic_check, alk_check  !fluxes
   real(rk) :: hplus, hco3, co3, co2sw, pco2sw
   real(rk) :: pH, boh4, oh, newdic, newalk, newalk2
   real(rk) :: tk
   real(rk) :: stand_dic, sal, temp, Tice, ice_hi                    !phys environmental variables
   real(rk) :: fno3phy,fde1nh4,fde2nh4,fde3nh4,fnh4phy,fnh4no3    !bio environmental variables   
   real(rk) :: spd = 86400.0_rk ! sec/day
   real(rk) :: dt
   real(rk) :: k1_carb, k2_carb, pK1, pK2 !carbonate chemistry constants
   real(rk) :: kb, kw, kw1, kw2, a1, a2, a3,a4, a5, a6
   real(rk) :: btot
   real(rk) :: Z
   real(rk) :: T0C=273.15_rk
   real(rk) :: CA, A, CO2diss, Kr, error, emax, epsilon
   real(rk) :: alpha, TTT, lnk0
   data epsilon /1.e-9/

   real(rk) :: botgrowth, botmelt, topgrowth, topmelt, termelt, fdic_ice, falk_ice
   real(rk) :: depth
!ericmod 2017jun08 adding coupling of ice algal growth to do routine
!   real(rk) :: fgrowIA
   real(rk) :: fno3upIA, fnh4upIA, fgrowIA                   ! env. variables for IA growth

   real(rk) :: dflv
! Declaring nml var's
   real(rk) :: dic_0, alk_0, dic_sw, alk_sw, dic_ice, alk_ice, ik_diff, ice_on, ik_on, IA_on, tplv, btlv, prop2sw, prop2sw_melt

   _LOOP_BEGIN_
! retrieve prognostic/state variables
   _GET_(self%id_dic,dic)
   _GET_(self%id_alk,alk)
   _GET_(self%id_eco_dic,eco_dic)


! ericmod 25aug16 adding new variables for ice brine rejection in do loop 
!!!               (need to take it OUT of surface loop)
   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_HORIZONTAL_(self%id_botmelt,botmelt)
   _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
   _GET_HORIZONTAL_(self%id_topmelt,topmelt)
   _GET_HORIZONTAL_(self%id_topgrowth,topgrowth)
   _GET_HORIZONTAL_(self%id_termelt,termelt)
!end ericmod 25aug16

! ericmod 2015jun08, adding...
!   _GET_HORIZONTAL_(self%id_fgrowIA, fgrowIA)

! retrieve environmental variables
   _GET_(self%id_depth,depth)
!   _GET_(self%id_stand_alk,stand_alk)
   _GET_(self%id_fno3phy,fno3phy)
   _GET_(self%id_fde1nh4,fde1nh4)
   _GET_(self%id_fde2nh4,fde2nh4)
   _GET_(self%id_fnh4phy,fnh4phy)
   _GET_(self%id_fnh4no3,fnh4no3)
   _GET_(self%id_temp,temp)
   _GET_(self%id_Tice, Tice)
   _GET_(self%id_sal,sal)
   _GET_GLOBAL_(self%id_dt,dt)

   prop2sw=self%prop2sw
   prop2sw_melt=self%prop2sw_melt


!!! carbonate chemistry equations...!!!(iterative part c''d out in this (non-surface) subroutine
   tk = temp+T0C
! 1st, calculating dissociation constants
   pK1 = 3633.86/tk - 61.2172 + 9.6777*log(tk) -       &
         0.011555*sal + 0.0001152*sal*sal
   pK2 = 471.78/tk + 25.9290 - 3.16967*log(tk) -       &
         0.01781*sal + 0.0001122*sal*sal
   k1_carb = 10**(-pK1) ! On the pHTOT scale in mol/kg-SW
   k2_carb = 10**(-pK2) ! On the pHTOT scale in mol/kg-SW

   kb=exp((-8966.9 - 2890.53*sal**0.5 - 77.942*sal + 1.728*sal**1.5 - 0.0996*sal**2)/tk  &
         + 148.0248 + 137.1942*sal**0.5 +  1.62142*sal                                   &
        - (24.4344 + 25.085*sal**0.5 + 0.2474*sal) * log(tk) + 0.053105*sal**0.5*tk);

   kw1 = 148.96502 - 13847.26/tk - 23.65218*log(tk)
   kw2 = (448.67/tk -5.977 + 1.0495*log(tk))*sal**0.5 - 0.01615*sal
   kw  = exp(kw1 + kw2)
   btot = (0.000232/10.811)*(sal/1.80655)

   TTT=tk/100.
   lnK0 = -60.2409 + 93.4517/TTT + 23.3585*log(TTT)
   lnK0 = lnK0 + sal*(0.023517-0.023656*TTT+0.0047036*TTT*TTT)
   alpha = exp(lnK0)   !alpha = k0


! iterative method to solve for hplus and alpha (k_0) in order to ultimately get pco2sw
! code came from 2013fall/res_13/gotm_carbon/src/extras/bio/bio_fasham.F90
   Kr = k1_carb/k2_carb
   A = alk

! ericmod 28jan15 commenting out the technique and seeing if everything still runs/compiles
!                taking it out bc not needed and huge drain of computational time...
!   do it = 1, itmax
!      Z = ((dic*Kr)**2 + A*Kr          &
!          *(2.*dic - A)                &
!          *(4.-Kr))**0.5
!          CO2diss = dic - A +          &
!          (A*Kr - dic*Kr - 4*A + Z)/   &
!          (2.*(Kr - 4.))
!      hplus = CO2diss*k1_carb/(2.*A) +       &
!             ((CO2diss*k1_carb)**2 +        &
!             8.*A*CO2diss*k1_carb*k2_carb)**0.5/ &
!             (2.*A)
!      CA = alk - btot*kb/(kb+hplus) - Kw/hplus + hplus
!
!      emax = 1.e12
!      error = abs(CA - A)
!      error = min(error,emax)
!      emax = error
!      A = CA
!      if(error .lt. epsilon) exit
!  enddo
!!!!!!!!end of carbonate chemistry stuff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   co2sw=dic/( 1 + k1_carb/hplus + k1_carb*k2_carb/(hplus*hplus) )
!   co3=dic/( 1 + hplus*hplus/(k1_carb*k2_carb) + hplus/k2_carb )
!   hco3=dic/( hplus/k1_carb + 1 + k2_carb/hplus )

!   pH=-log(hplus)
!   boh4=kb*btot/(hplus+kb)
!   oh=kw/hplus
!   newdic=hco3+co3+co2sw
!   newalk=hco3+2*co3
!   newalk2=newalk+boh4+oh-hplus

!   write(*,*) 'dic, newdic', dic,newdic
!   write(*,*) 'alk, newalk, newalk2', alk, newalk, newalk2



!!!PELAGIC BIOLOGY sink for DIC and then TA
! UNITS NOTE: fno3phy =[mmol N/m3/day], so fdic = (fno3phy+...)*(mol C/mol N)
!    for TA, no need for 106/16 conversion 

!ericmod 30aug2017 adding nh4 uptake!!!
   fdic1=(-fnh4phy-fno3phy+fde1nh4+fde2nh4)*106.0_rk/16.0_rk/self%spd
!   fdic1=(-fno3phy+fde1nh4+fde2nh4)*106.0_rk/16.0_rk/self%spd

! NEWEST ericmod, ADDING THE CHANGE OF BELOW BIO CHANGES TO ALKALINITY
! based on Wolf-Gladrow et al., 2007 (from Goldman and Brewer, 1980)

! Process:                             Alkalinity change:  
! 1 mol uptake of nitrate or nitrite   +1 mol
! 1 mol uptake of ammonium             -1 mol
! 1 mol nitrification                  -2 mol
! 1 mol remineralization to ammonium   +1 mol
! 1 mol uptake of phosphorus           +1 mol
! 1 mol uptake of sulfur               +2 mol
! remineralization does the inverse

! ericmod 20jun2017: the last 2 uptake to Alk relations occur on a 1/16 and 2*2.4/16 ratios coincident with 
!                         nitrate (and ammonium? No, right???) uptake, and lead to a change in these values:
!                         +fno3phy-fnh4phy -> (1+5.8/16)*(fno3phy)-fnh4phy
!                         CORRESPONDINGLY in Wolfe-Gladrow2007 is a drop during remineralization which should change:
!                         +fde1nh4 +fde2nh4 -> (1+5.8/16)*(fde1nh4+fde2nh4)
! BUT BIG QUESTION!! Should we include this alkalinity effect due to P and S uptake and remineralization
! AND BIG Q #2!! Does this lead to a conservation of alkalinity?!?
! A to BIG Q#2:[N]: 1 mol fnh4phy       , 1 mol fno3phy  -> 2 mol detnh4    -> 1 mol nh4no3
!   change in [TA]:  -1 mol             , + 1 mol        -> +2 mol          -> -2 mol       = net 0 change in TA
!     case II [TA]:  -1 + 5.8/16 mol    , + 1+5.8/16 mol -> +2(1+5.8/16) mol-> -2 mol       = net change in TA of 3*5.8/16 
!ericmod 2017jul13: adding the +5.8/16 to nh4 uptake:
!old:
!   falk1=( (1+5.8/16)*fno3phy-fnh4phy-fnh4no3*2.0_rk+(1+5.8/16)*(fde1nh4+fde2nh4) )/self%spd
!new:
   falk1=( (1+5.8/16)*fno3phy+(-1+5.8/16)*fnh4phy-fnh4no3*2.0_rk+(1+5.8/16)*(fde1nh4+fde2nh4) )/self%spd

! Below is old comments from above 2017jun20 correction
! we take that into account by assuming a N:P ratio of 16:1 
! and a N:S ratio of 16:2.4
! In total this gives a N:Alk-ratio of 16:(1+2*2.4) = 16:5.8 = 1:0.3625 = 2.76

!   falk1=(fno3phy-fnh4phy-fnh4no3*2.0_rk+fde1nh4+fde2nh4)/self%spd
! end ericmod 2017jun20


   if(ice_hi.gt.0.0_rk) then  
! as in Roullet20000 and Moreau2016, dic_sw, dic_ice, ta_sw, ta_ice are ref. values
      fdic_ice=(prop2sw_melt*(-botmelt-topmelt-termelt)+prop2sw*(topgrowth+botgrowth))*(self%dic_sw-self%dic_ice)
      falk_ice=(-botmelt-topmelt-termelt+topgrowth+botgrowth)*(self%alk_sw-self%alk_ice)
      if (botmelt+topmelt+termelt.gt.topgrowth+botgrowth) then           !when ICE is MELTING
         if (self%ik_on.eq.1) then  !adding the on/off switch for ikaite precipitation 
             fdic_ice=(prop2sw_melt*(-botmelt-topmelt-termelt)+prop2sw*(topgrowth+botgrowth)) * (self%dic_sw - (self%dic_ice-self%ik_diff) )
             falk_ice=(-botmelt-topmelt-termelt+topgrowth+botgrowth) * (self%alk_sw - (self%alk_ice - self%ik_diff*2) )
         endif
! ericmod 11sep16 (nested inside 25aug16 ericmod)
! adding ikaite dissolution in the water column (at below 50m)
      elseif (botmelt+topmelt+termelt.lt.topgrowth+botgrowth) then       !when ICE is GROWING
         if (self%ik_on.eq.1) then  !adding the on/off switch for ikaite precipitation/dissolution
            if (depth.lt.1.4) then     !depth of top of ikaite uniform dissolution layer
               fdic_ice = (prop2sw_melt*(-botmelt-topmelt-termelt)+prop2sw*(topgrowth+botgrowth))*(self%ik_diff)
               falk_ice = (-botmelt-topmelt-termelt+(topgrowth+botgrowth))*(self%ik_diff*2)
           endif
         endif
      endif
   endif

! ericmod 9apr18
! DIC release at depth of brine-injected dic layer during ice growth, and at surf during ice melt
 tplv=self%tplv
 btlv=self%btlv
 dflv=btlv-tplv
  if(ice_hi.gt.0.0_rk) then
      if (depth.gt.tplv .AND. depth.lt.btlv .AND. topgrowth+botgrowth.gt.botmelt+topmelt+termelt) then      !when ICE is GROWING
         fdic_ice = (prop2sw*(-botmelt-topmelt-termelt)+prop2sw*(topgrowth+botgrowth)) * (self%dic_sw-self%dic_ice-self%ik_diff)
         falk_ice = (-botmelt-topmelt-termelt+(topgrowth+botgrowth)) * (self%alk_sw-self%alk_ice-self%ik_diff*2)
         !_SET_ODE_(self%id_dic, fdic_ice/dflv + fdic1) !jpnote changing set_ode to add source
         _ADD_SOURCE_(self%id_dic, fdic_ice/dflv + fdic1)
        ! _SET_ODE_(self%id_alk, falk_ice/dflv + falk1)
         _ADD_SOURCE_(self%id_alk, falk_ice/dflv + falk1)
      elseif (depth.lt.1.4_rk .AND. topgrowth+botgrowth.lt.botmelt+topmelt+termelt) then    !ICE MELT (dilution at surf., only ikaite here)
       !!!!I************ the dilution difference btw ice and sw dic/ta is done in surf subrout below...
       !!!!I************ also note, that the present formulation for ikaite is dissolved and precipitated at the same place, so net 0 effect!!!!
       !!!! BELOW HAVE BEEN COMMENTED OUT AND REPLACED W =0, BC ALREADY DONE IN SURF SUBROUTINE (BELOW)
!         fdic_ice=(prop2sw_melt*(-botmelt-topmelt-termelt)+prop2swi_melt*(topgrowth+botgrowth)) * (-self%dic_sw-self%dic_ice-self%ik_diff )
!         falk_ice=(-botmelt-topmelt-termelt+topgrowth+botgrowth) * (-self%alk_sw-self%alk_ice-self%ik_diff*2 )
         fdic_ice=0
         falk_ice=0
        ! _SET_ODE_(self%id_dic, fdic_ice + fdic1)
         _ADD_SOURCE_(self%id_dic, fdic_ice + fdic1)
         !_SET_ODE_(self%id_alk, falk_ice + falk1)
         _ADD_SOURCE_(self%id_alk, falk_ice + falk1)
      else
        ! _SET_ODE_(self%id_dic, fdic1)
         _ADD_SOURCE_(self%id_dic, fdic1)
        ! _SET_ODE_(self%id_alk, falk1)
         _ADD_SOURCE_(self%id_alk, falk1)
      endif
  else if(ice_hi.eq.0.0_rk) then
     ! _SET_ODE_(self%id_dic,fdic1)
      _ADD_SOURCE_(self%id_dic,fdic1)
      !_SET_ODE_(self%id_alk,falk1)
      _ADD_SOURCE_(self%id_alk,falk1)
  endif
! end ericmod 9apr18

! 2017ericmod ...
!  want to make: fIA_CO2 = fgrowIA*k_IA2N*106/16*0.03_rk/self%spd

   !_SET_SURFACE_EXCHANGE_(self%id_dic,-1.0) !fgrowIA*106/16*0.03/86400)



   !_SET_ODE_(self%id_eco_dic,fdic1)
   _ADD_SOURCE_(self%id_eco_dic,fdic1)

! save diagnostic variables
   _SET_DIAGNOSTIC_(self%id_fdic1,fdic1)
   _SET_DIAGNOSTIC_(self%id_falk1,falk1)
   _SET_DIAGNOSTIC_(self%id_hplus,hplus)
   _SET_DIAGNOSTIC_(self%id_co2sw,co2sw)
   _SET_DIAGNOSTIC_(self%id_co3,co3)
   _SET_DIAGNOSTIC_(self%id_hco3,hco3)
! ericmmod 5dec14 changed stand_dic below to newdic
   _SET_DIAGNOSTIC_(self%id_dic_check,newdic)
! ericmmod 5dec14 changed stand_alk below to newalk2
   _SET_DIAGNOSTIC_(self%id_alk_check,newalk2)
!   _SET_DIAGNOSTIC_(self%id_pH,pH)
   _SET_DIAGNOSTIC_(self%id_newalk,newalk)
   _SET_DIAGNOSTIC_(self%id_newdic,newdic)
   _LOOP_END_
   end subroutine do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            do_surface subroutine                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALSO, where ice algae should be, delta_DIC = -k*growth_IA
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_dic),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_
   INTEGER  :: it
   INTEGER  :: itmax = 200
   real(rk) :: dic, alk, eric2, dt  !state variables
   real(rk) :: btot
   real(rk) :: fdic1  !fluxes
   real(rk) :: Tatm, Patm, pco2atm, wind, precip        ! env. variables in atm
   real(rk) :: temp, sal, press, Tice                 ! env. variables at sea surf.
   real(rk) :: fgrowIA, fnh4upIA, fno3upIA, fnitIA, ice_hi                  ! env. variables for IA growth
   real(rk) :: botgrowth, botmelt, topgrowth, topmelt, termelt
   real(rk) :: pco2sw, pH                          ! diag. variable
   real(rk) :: stemp, rho,T,T2,T3,T4,T5,S15,S,S2,S3,x,p2,p,K
   real(rk) :: sdic, Tkelv, lnh2o, ph2o
   real(rk) :: p_atm, tak
   real(rk) :: pres_atm, pco2_a, cosat_o, VmCO2, sc, pv
   real(rk) :: prop2sw, prop2sw_melt
   real(rk) :: surftemp 
   real(rk) :: k_gasex, delta_e, co2flux, fIA_co2, fIA_ta, fIA_alk, fdic_ice, falk_ice
   real(rk) :: ik_ice, k_ik
   real(rk) :: aco2(4)
   real(rk) :: T0C=273.15_rk
   real(rk) :: tk, TTT, pK1, pK2, k1_carb, k2_carb, lnk0
   real(rk) :: kb, kw1, kw2, kw, Z, A, CA, Kr, CO2diss, hplus, error, emax 
   real(rk) :: alpha          !=0.065_rk  !???should be converted to...???
   real(rk) :: co2_thinice, thinice_lim
   data aco2/2073.1,125.62,3.6276,0.043219/
   real,parameter            :: U_CO2=49.0
   real(rk) :: k_IA2N=0.355_rk  ! 0.355 mmol N/1 mg chl-a, conversion of mg chl-a to mmol N
   real(rk) :: k_N2C=106.0_rk/16.0_rk !C/N , conversion of N to C
   real(rk) :: spd = 86400.0_rk ! sec/day
   real(rk) :: epsilon
   data epsilon /1.e-9/
   integer :: icepump=0
 
!sdic -> mol/kg
!tak -> mol/kg
!alpha2 -> mol/kg/uatm
!co2sw_out -> uatm
!  Register namelist parameters


   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_dic,dic)
   _GET_(self%id_alk,alk)
   _GET_HORIZONTAL_(self%id_Tatm,Tatm)
   _GET_HORIZONTAL_(self%id_Patm,Patm)
   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_HORIZONTAL_(self%id_wind,wind)
! ericmod 6feb15, adding ice bottom melt/growth
   _GET_HORIZONTAL_(self%id_botmelt,botmelt)
   _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
   _GET_HORIZONTAL_(self%id_topmelt,topmelt)
   _GET_HORIZONTAL_(self%id_topgrowth,topgrowth)
   _GET_HORIZONTAL_(self%id_termelt,termelt)
! ericmod 30jan2018 removing precip variable everywhere in dic.F90
!   _GET_HORIZONTAL_(self%id_precip,precip)
   _GET_HORIZONTAL_(self%id_stemp,stemp)
   _GET_(self%id_temp,temp)
   _GET_HORIZONTAL_(self%id_surftemp, surftemp) ! ericmod 1sep15 adding ice surf temp
   _GET_(self%id_Tice, Tice)
   _GET_(self%id_sal,sal)
   _GET_(self%id_press,press)

   _GET_HORIZONTAL_(self%id_fgrowIA,fgrowIA)
!ericmod 2017jul12 adding the no3 uptake of icealgae as a dependency, in order to differentiate
! between new and recycled PP production for opposite effect on TA
   _GET_HORIZONTAL_(self%id_fno3upIA,fno3upIA)
   _GET_HORIZONTAL_(self%id_fnitIA, fnitIA)

   _GET_GLOBAL_(self%id_dt,dt)
   
   prop2sw=self%prop2sw
   prop2sw_melt=self%prop2sw_melt


      pres_atm=Patm*9.86923267e-06 !unit conversion: from Pa to atm
      pco2_a = pres_atm*400    !partial press of co2 = 0.040% of atm
!ericmod14mar  print2scrn to see what the size of pco2_a (as compared to pco2_sw = 400)
!write(*,*) '(1)Patm(Pa), pres_atm(atm), pco2_a(uatm) =', Patm, pres_atm, pco2_a
      T=temp
      S=sal  
!      p=0_rk  !what is this???
   T2 = T*T
   T3 = T*T2
   T4 = T2*T2
   T5 = T*T4
   S15= S**1.5
   S2 = S*S
   S3 = S*S2

!calculating density at surf
   x=999.842594+6.793952e-02*T-9.09529e-03*T2+1.001685e-04*T3
   x=x-1.120083e-06*T4+6.536332e-09*T5
   x=x+S*(0.824493-4.0899e-03*T+7.6438e-05*T2-8.2467e-07*T3)
   x=x+S*5.3875e-09*T4
   x=x+sqrt(S3)*(-5.72466e-03+1.0227e-04*T-1.6546e-06*T2)
   x=x+4.8314e-04*S2
   rho=x
! end of density calculation

! convert surface alk to mol/g from umol/L
!   tak=alk*1.e-9/rho
! convert surface dic to mol/kg from umol/L
   sdic=dic    !*1.e-9/rho

!!!!!!!!! carbonate chemistry equations...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   tk = temp+T0C
   pK1 = 3633.86/tk - 61.2172 + 9.6777*log(tk) -       &
         0.011555*sal + 0.0001152*sal*sal
   pK2 = 471.78/tk + 25.9290 - 3.16967*log(tk) -       &
         0.01781*sal + 0.0001122*sal*sal
   k1_carb = 10**(-pK1) ! On the pHTOT scale in mol/kg-SW
   k2_carb = 10**(-pK2) ! On the pHTOT scale in mol/kg-SW

   kb=exp((-8966.9 - 2890.53*sal**0.5 - 77.942*sal + 1.728*sal**1.5 - 0.0996*sal**2)/tk  &
         + 148.0248 + 137.1942*sal**0.5 +  1.62142*sal                                   &
        - (24.4344 + 25.085*sal**0.5 + 0.2474*sal) * log(tk) + 0.053105*sal**0.5*tk);

   kw1 = 148.96502 - 13847.26/tk - 23.65218*log(tk)
   kw2 = (448.67/tk -5.977 + 1.0495*log(tk))*sal**0.5 - 0.01615*sal
   kw  = exp(kw1 + kw2)
   btot = (0.000232/10.811)*(sal/1.80655)
! iterative method to solve for hplus and alpha (k0) in order to ultimately get pco2sw
! code came from 2013fall/res_13/gotm_carbon/src/extras/bio/bio_fasham.F90
   TTT=tk/100.
   lnK0 = -60.2409 + 93.4517/TTT + 23.3585*log(TTT)
   lnK0 = lnK0 + sal*(0.023517-0.023656*TTT+0.0047036*TTT*TTT)
   alpha = exp(lnK0)

   Kr = k1_carb/k2_carb
   A = alk
   do it = 1, itmax
      Z = ((sdic*Kr)**2 + A*Kr          &
          *(2.*sdic - A)                &
          *(4.-Kr))**0.5
          CO2diss = sdic - A +          &
          (A*Kr - sdic*Kr - 4*A + Z)/   &
          (2.*(Kr - 4.))
      hplus = CO2diss*k1_carb/(2.*A) +       &
             ((CO2diss*k1_carb)**2 +        &
             8.*A*CO2diss*k1_carb*k2_carb)**0.5/ &
             (2.*A)
      CA = alk - btot*kb/(kb+hplus) - Kw/hplus + hplus

      emax = 1.e12
      error = abs(CA - A)
      error = min(error,emax)
      emax = error
      A = CA
      if(error .lt. epsilon) exit
  enddo
! below will give pco2sw in units of micro-atm 
!  pco2sw = sdic/(alpha*1.e-6*(1.+k1_carb/hplus+(k1_carb*k2_carb)/hplus**2))
! dropped 1.e-6 in pco2sw eq here and in sdic def. above
  pco2sw = sdic/(alpha*(1.+k1_carb/hplus+(k1_carb*k2_carb)/hplus**2))
  pH=-log10(hplus)
!!!!!!!!end of carbonate chemistry stuff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!   calculating air-sea co2flux INTO seawater   !!!!!
! units are microatm*m/s, need to convert to mmol/m2/s for surface exchange
! this script was clipped from GOTM carbon stuff given to me by Nadja, ref's for bubbles etc. are at ....?
   tk=temp+T0C  !unit conversion celsius to kelvin
   lnh2o=24.4543-67.4509*100./tk-4.8489*log(tk*0.01)-0.000544*sal
   ph2o=exp(lnh2o)
   !end ericmod 27nov13
   pco2_a=pco2_a*(pres_atm-ph2o)/(1.-ph2o) !3
   cosat_o=1.e-6*alpha*VmCO2*rho        !4     
   cosat_o=cosat_o*tk/T0C*1./pres_atm      !5
!   k_gasex=pv(wspd,tsurf,cosat_o,aco2) !6
   sc= aco2(1) - aco2(2)*temp + aco2(3)*temp*temp - aco2(4)*temp*temp*temp
   pv=0.31*wind*wind/sqrt(sc/660.)
   pv=0.24*pv     ! convert to m/d
   k_gasex=pv
   delta_e=0.01*(wind/U_CO2)*(wind/U_CO2)  !7
   co2flux=alpha*(pco2_a*(1+delta_e)-pco2sw)*k_gasex/spd !units = mmol/m2/s

!ericmod 6sep2017 adding cubic pv as in Wanninkhof and McGillis 1999
   pv = 0.0283*wind*wind*wind/sqrt(sc/660.)*0.24
   co2flux=alpha*(pco2_a-pco2sw)*pv/spd

! and changing to Wanninkhof 1992
   pv = 0.39*wind*wind/sqrt(sc/660.)*0.24
   co2flux=alpha*(pco2_a-pco2sw)*pv/spd
   
!write(*,*) 'a,delp,pv,flux,wind',alpha, (pco2_a-pco2sw), pv, spd*co2flux, wind
!end ericmod 6sep2017




!ericmod14mar  print2scrn to see what the size of pco2_a (as compared to pco2_sw = 400)
!write(*,*) 'pco2_a, delta_e =', pco2_a, delta_e

! adding  fdic_ice during ICEMELT, which should be in units of dic concentration (mmol/m3) * melt rate (m/d, right???), fdic_ice=mmol/m2/d or (is this /s, not /d ???)
   if(ice_hi.gt.0.0) then
      if (botmelt+topmelt+termelt.lt.topgrowth+botgrowth) then   !ICE GROWTH, NO surface flux (DIC goes to brine depth)
         fdic_ice=0.0_rk
         falk_ice=0.0_rk
      elseif (botmelt+topmelt+termelt.gt.topgrowth+botgrowth) then   !ICE MELT
         fdic_ice=(prop2sw_melt*(-botmelt-topmelt-termelt)+prop2sw_melt*(topgrowth+botgrowth)) * (self%dic_sw-self%dic_ice-self%ik_diff)
         falk_ice=(-botmelt-topmelt-termelt+topgrowth+botgrowth) * (self%alk_sw-self%alk_ice-self%ik_diff*2)
      endif
   endif

! CO2 uptake by ICE ALGAL GROWTH
! UNITS NOTE: fgrowIA = [mmol N/m3/dy], so fIA_CO2 = fgrowIA*k_N2C*d_IA/d_wc_surf/spd
!                       k_N2C=106/12
!           kIA2N= 0.355 is NOT NEEDED!!! (i checked, and this wasn't very important bc ice algae are still relatively unimportant
! pre-22feb2017:   fIA_co2=fgrowIA*k_IA2N*k_N2C*0.03_rk/self%spd
! pre-22feb2017:   fIA_alk= N/A!!!!

!2017jun08 changing DIC uptake by the above notes to put it into correct C units
!   fIA_co2=fgrowIA*k_IA2N*0.03_rk/self%spd
   fIA_co2=fgrowIA*k_N2C*0.03_rk/self%spd
! and the below is wrong, up to 11jul2017
!   fIA_alk=fgrowIA*k_IA2N*0.03_rk/self%spd
! correctd alkalinity changes include: 
! 1. sulfate and phosphate uptake (at Redfield and molar ratios of (1*1 + 2*2.4):16 mole nitrate = 5.8/16)
! 2. ratio for nitrogen to carbon (k_N2C=6.625) rather than k_IA2N = 0.355
!  that all being said,  TA uptake due to ice algal growth is still only a minor contributor to carbonate system

!ericmod 2017jul12 adding the no3 uptake of icealgae as a dependency, in order to differentiate
! between new and recycled PP production for opposite effect on TA
! old:  
!   fIA_alk=fgrowIA*k_N2C*(1.0_rk + 5.8_rk/16_rk)*0.03_rk/self%spd
! new: (note: (fgrowIA-fno3upIA) should = fnh4upIA ), AND (note: dropped k_N2C conversion of Nitrogen to Carbon)
!   fIA_alk= ( fno3upIA*(1.0_rk + 5.8_rk/16.0_rk) + (fgrowIA-fno3upIA)*(-1.0_rk + 5.8_rk/16.0_rk) )*0.03_rk/self%spd
! end ericmod 2017jul12

!ericmod 4sep2017 adding ice algal contribution to TA (here and with coupling to icealgae.F90 code above)
   fnh4upIA = (fgrowIA-fno3upIA)
!   fIA_ta = (fnh4upIA*(-1+5.8/16) + fno3upIA*(+1+5.8/16) )*0.03_rk/self%spd
   fIA_ta = (fnh4upIA*(-1+5.8/16) + fno3upIA*(+1+5.8/16) - fnitIA*2 )*0.03_rk/self%spd
!end ericmod 4sep2017
! pelagic for comparison: falk1=( (1+5.8/16)*fno3phy+(-1+5.8/16)*fnh4phy-fnh4no3*2.0_rk+(1+5.8/16)*(fde1nh4+fde2nh4) )/self%spd

! This is a way to allow co2 through thin ice, not used now (because thinice_lim=0)
   co2_thinice=0.5_rk !this is proportion of openwater co2 flux, 
                      !i.e. 1(0.5) = 100% (50%) icefree co2 flux
   thinice_lim=0.0_rk !max depth of thinice that is permeable to co2 flux btw air-sea


!switch to turn on/off ice algal carbon sink pump
   if (ice_hi.lt.0.1_rk) then
     fIA_co2 = 0.0_rk
     fIA_ta  = 0.0_rk
     fdic_ice = 0.0_rk
     falk_ice = 0.0_rk
   endif

!ericmod 2017jun08 fixes of ice algal flux
!write(*,*) fIA_co2
!NO ICE, airsea flux only
   if(ice_hi.eq.0.0_rk) then
     ! _SET_SURFACE_EXCHANGE_(self%id_dic, co2flux)
      _ADD_SURFACE_FLUX_(self%id_dic, co2flux)
   else if((ice_hi.gt.0) .AND. (ice_hi.lt.thinice_lim)) then
     ! _SET_SURFACE_EXCHANGE_(self%id_dic, co2_thinice*co2flux - fIA_co2 + fdic_ice)
      _ADD_SURFACE_FLUX_(self%id_dic, co2_thinice*co2flux - fIA_co2 + fdic_ice)
      !_SET_SURFACE_EXCHANGE_(self%id_alk, +fIA_ta + falk_ice)
      _ADD_SURFACE_FLUX_(self%id_alk, +fIA_ta + falk_ice)
   else if(ice_hi.gt.thinice_lim) then                                  !dic flux due to ice melt/growth
     ! _SET_SURFACE_EXCHANGE_(self%id_dic, -fIA_co2 + fdic_ice)
      _ADD_SURFACE_FLUX_(self%id_dic, -fIA_co2 + fdic_ice)
      !_SET_SURFACE_EXCHANGE_(self%id_alk, +fIA_ta + falk_ice)
      _ADD_SURFACE_FLUX_(self%id_alk, +fIA_ta + falk_ice)
      co2flux=0                                                         !killing airsea exchange of co2 when ice is present
   endif


   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_co2flux, co2flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fIA_co2, fIA_co2)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdic_ice, fdic_ice)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pH,pH)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pco2sw,pco2sw)
   _HORIZONTAL_LOOP_END_


   end subroutine do_surface




end module uvic_dic

