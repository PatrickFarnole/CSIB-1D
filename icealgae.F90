!jp copy from mortenson
!-----------------------------------------------------
!
! University of Victoria icealgae model
!
! Contributors: Hakase Hayashida, Nadja Steiner, Eric Mortenson, Adam Mohanan
!
! Model description: The model will be described in a paper.
!
!------------------------------------------------------

#include "fabm_driver.h"

module uvic_icealgae

   use fabm_types
   use fabm_expressions  
   use fabm_standard_variables

   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_icealgae
! Declare horizontal prognostic variables      
      type (type_surface_state_variable_id) :: id_no3,id_sil,id_ia,id_nh4
! Declare horizontal diagnostic variables
      type (type_horizontal_diagnostic_variable_id) :: id_chl,id_chlia,id_fgrow,id_fgraze,id_fmort,id_fmort2,id_fmelt,id_fpond,id_fpondno3,id_fpondnh4,id_fpondsil,id_fno3up,id_fsilup,id_fskelno3,id_fskelnh4,id_fskelsil,id_ier,id_lice,id_llig,id_lno3,id_lsil,id_hnu,id_fnit
! Declare environmental variables
      type (type_dependency_id) :: id_ph2,id_no3SW,id_nh4SW,id_silSW,id_u,id_v
      !mortenson
      !type (type_horizontal_dependency_id) :: id_airt,id_ice_hs,id_par,id_temp,id_ice_hi,id_botmelt,id_botgrowth,id_topmelt,id_termelt,id_Amelt
      !hayashida
      type (type_horizontal_dependency_id) :: id_airt,id_ice_hs,id_par,id_temp,id_ice_hi,id_botmelt,id_botgrowth,id_topmelt,id_termelt,id_Amelt,id_ts !,id_bvf
      type (type_global_dependency_id) :: id_dt 

! Declare model parameters
      real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,rnit,skno3_0,sknh4_0,sksil_0,ia_0,ia_b,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,t_sens,nu,md_no3,md_sil,chl2n,sil2n
! Declare anything else used in all procedures
      real(rk) :: spd = 86400.0_rk ! Seconds Per Day (spd)

      contains

      procedure :: initialize
      procedure :: do_surface
      
   end type

contains

   subroutine initialize(self,configunit)
   class (type_uvic_icealgae), intent(inout), target :: self
   integer, intent(in)                          :: configunit
! Declare yaml parameters
  real(rk) :: ia_0,skno3_0,sknh4_0,sksil_0
   
   !read in vals from fabm.yaml 

   call self%get_parameter(self%r_pond, 'r_pond', '','melt pond drainage rate',default=0.0175_rk,scale_factor=1.0_rk/self%spd)
  ! self%r_pond = self%r_pond/self%spd
   call self%get_parameter(self%fmethod, 'fmethod','', 'method for ice-ocean flux',default=0.0_rk)
   call self%get_parameter(self%fflush , 'fflush','', 'method for flushing', default=0.0_rk)
   call self%get_parameter(self%drag , 'drag','-', 'drag coefficient at the ice-water interface ', default=0.005_rk)
   call self%get_parameter(self%f_graze, 'f_graze','-', 'fraction of ice algal growth lost due to grazing ', default=0.1_rk)
   call self%get_parameter(self%zia, 'zia','m', 'ice algal layer thickness ', default=0.03_rk)
   call self%get_parameter(self%ac_ia, 'ac_ia','', 'specific light attenuation coefficient for ice algae', default=0.03_rk)
   call self%get_parameter(self%rnit , 'rnit','per day', 'nitrification rate ', default=0.1_rk,scale_factor=1.0_rk/self%spd)
  ! self%rnit  = self%rnit/self%spd
 !  call self%get_parameter(ia_0 , 'ia_0','mmol-N/m3', 'ia initial value ', default=0.16_rk)
   call self%get_parameter(self%ia_b , 'ia_b','mmol-N/m3',  'ia background value ',default=0.01_rk)
!   call self%get_parameter(skno3_0, 'skno3_0','mmol/m3', 'no3 initial value ', default=2.0_rk)
 !  call self%get_parameter(sknh4_0, 'sknh4_0','mmol/m3', 'nh4 initial value ', default=0.01_rk)
 !  call self%get_parameter(sksil_0, 'sksil_0','mmol/m3', 'sil initial value ', default=5.0_rk)
   call self%get_parameter(self%ks_no3, 'ks_no3','mmol/m3', 'no3 half-saturation value ',default=1.0_rk)
   call self%get_parameter(self%ks_sil, 'ks_sil','mmol/m3', 'sil half-saturation value ', default=4.0_rk)
   call self%get_parameter(self%maxg, 'maxg','d-1', 'maximum specific growth rate ', default=0.8511_rk,scale_factor=1.0_rk/self%spd)
  !! self%maxg = self%maxg/ self%spd
   call self%get_parameter(self%mort , 'mort','d-1', 'linear mortality rate', default=0.05_rk,scale_factor=1.0_rk/self%spd)
  ! self%mort = self%mort/ self%spd
   call self%get_parameter(self%mort2, 'mort2','d-1',  'quadratic mortality rate ',default=0.05_rk,scale_factor=1.0_rk/self%spd)
  ! self%mort2 = self%mort2/ self%spd
   call self%get_parameter(self%crit_melt, 'crit_melt','m d-1', 'critical melt rate [m d-1]', default=0.015_rk,scale_factor=1.0_rk/self%spd)
  ! self%crit_melt = self%crit_melt / self%spd
   call self%get_parameter(self%lcompp, 'lcompp','umol m-2 s-1', '# compensation intensity ', default=0.4_rk)
   call self%get_parameter(self%rpp , 'rpp','[W m-2]-1', 'ratio of photosynthetic parameters (alpha and pbm) [W m-2]-1', default=0.1_rk)
   !call self%get_parameter(self%rpi , 'rpi ', 'ratio of photoinhibition parameters (beta and pbm)', default=0)
   call self%get_parameter(self%t_sens , 't_sens','deg.C-1', 'temperature sensitivity ', default=0.0633_rk)
   call self%get_parameter(self%nu , 'nu','', 'kinematic viscosity?', default=1.86e-6_rk)
   call self%get_parameter(self%md_no3, 'md_no3','', 'molecular diffusion coefficient for nitrate', default=0.47e-9_rk)
   call self%get_parameter(self%md_sil , 'md_sil','', 'molecular diffusion coefficient for dissolved silica', default=0.47e-9_rk)
   call self%get_parameter(self%chl2n , 'chl2n','', 'chl to nitrogen ratio', default=2.8_rk)
   call self%get_parameter(self%sil2n , 'sil2n','', 'silicon to nitrogen ratio', default=1.7_rk)
      
   ! Register prognostic variables
#if 0
      call self%register_state_variable(self%id_no3,'no3','mmol m-3','skel. NO_3',initial_value=skno3_0,minimum=0.0_rk)
      call self%register_state_variable(self%id_sil,'sil','mmol m-3','skel. Si',initial_value=sksil_0,minimum=0.0_rk) 
      call self%register_state_variable(self%id_ia,'ia','mmol m-3','Ice algae',initial_value=ia_0,minimum=0.0_rk) 
      call self%register_state_variable(self%id_nh4,'nh4','mmol m-3','NH4',initial_value=sknh4_0,minimum=0.0_rk) 
#endif
      call self%register_state_variable(self%id_no3,'no3','mmol m-3','skel. NO_3',minimum=0.0_rk)
      call self%register_state_variable(self%id_sil,'sil','mmol m-3','skel. Si',minimum=0.0_rk) 
      call self%register_state_variable(self%id_ia,'ia','mmol m-3','Ice algae',minimum=0.0_rk) 
      call self%register_state_variable(self%id_nh4,'nh4','mmol m-3','NH4',minimum=0.0_rk) 

! Register diagnostic variables
      call self%register_horizontal_diagnostic_variable(self%id_chl,'chl','mg m-3','Ice algae in per cubic meter',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_chlia,'chlia','mg m-2','Ice algae in per square meter',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fgrow,'fgrow','mmol-N m-3 d-1','ice algal growth rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fgraze,'fgraze','mmol-N m-3 d-1','ice algal grazing rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fmort,'fmort','mmol-N m-3 d-1','ice algal mortality rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fmort2,'fmort2','mmol-N m-3 d-1','ice algal quadratic mortality rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fmelt,'fmelt','mmol-N m-3 d-1','ice algal loss rate due to ice melting',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fpond,'fpond','mmol-N m-3 d-1','ice algal loss rate due to meltpond drainage',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fpondno3,'fpondno3','mmol-N m-3 d-1','no3 loss rate due to meltpond drainage',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fpondnh4,'fpondnh4','mmol-N m-3 d-1','nh4 loss rate due to meltpond drainage',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fpondsil,'fpondsil','mmol-N m-3 d-1','sil loss rate due to meltpond drainage',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fnit,'fnit','mmol-N m-3 d-1','nitrification in skeletal layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fno3up,'fno3up','mmol m-3 d-1','ice algal no3 uptake rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fsilup,'fsilup','mmol m-3 d-1','ice algal sil uptake rate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fskelno3,'fskelno3','mmol m-3 d-1','no3 flux to skeletal layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fskelnh4,'fskelnh4','mmol m-3 d-1','nh4 flux to skeletal layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fskelsil,'fskelsil','mmol m-3 d-1','sil flux to skeletal layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_ier,'ier','m d-1','Ice evolution rate (+ means growth, - means melting)',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_lice,'lice','-','Limitation due to ice growth/melt',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_llig,'llig','-','Limitation due to light',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_lno3,'lno3','-','Limitation due to nitrate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_lsil,'lsil','-','Limitation due to silicate',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_hnu,'hnu','m','Molecular sublayer thickness',source=source_do_horizontal)
      ! Register environmental variables
     
      call self%register_horizontal_dependency(self%id_temp,standard_variables%sea_ice_temperature) 
     
      !hayashida
      call self%register_horizontal_dependency(self%id_ts,standard_variables%surface_ice_temperature)
      !call self%register_horizontal_dependency(self%id_bvf,standard_variables%brine_volume_fraction) 

      call self%register_horizontal_dependency(self%id_ice_hi,standard_variables%sea_ice_thickness)
      call self%register_horizontal_dependency(self%id_ice_hs,standard_variables%snow_thickness)
      call self%register_horizontal_dependency(self%id_par,standard_variables%lowest_ice_layer_PAR)
      call self%register_horizontal_dependency(self%id_airt,standard_variables%surface_temperature)
      call self%register_horizontal_dependency(self%id_topmelt,standard_variables%topmelt)
      call self%register_horizontal_dependency(self%id_termelt,standard_variables%termelt)
      call self%register_horizontal_dependency(self%id_Amelt,standard_variables%f_melt)
      call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt)
      call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)

      call self%register_dependency(self%id_u,  standard_variables%zonal_current)
      call self%register_dependency(self%id_v,standard_variables%meridional_current) 
      call self%register_global_dependency(self%id_dt,standard_variables%timestep) 

     ! call self%register_dependency 

      call self%register_dependency(self%id_ph2,'uvic_eco_ph2','','')
      call self%register_dependency(self%id_no3SW,'uvic_eco_no3','','') 
      call self%register_dependency(self%id_nh4SW,'uvic_eco_nh4','','')
      call self%register_dependency(self%id_silSW,'uvic_eco_sil','','')
      call self%request_coupling(self%id_ph2,'uvic_eco_ph2')
      call self%request_coupling(self%id_no3SW,'uvic_eco_no3')
      call self%request_coupling(self%id_nh4SW,'uvic_eco_nh4')
      call self%request_coupling(self%id_silSW,'uvic_eco_sil')

   end subroutine initialize

! Calling and saving the surface variables.
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_icealgae),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk) :: ia,no3,nh4,sil,temp,par,ice_hi,ice_hs,k_snow,albedo,airt,ph2,no3SW,nh4SW,silSW,topmelt,termelt,botmelt,botgrowth,u,v,utaui,lice,llig,lno3,lsil,mum,hnu,grow
   !hayashida 
   real(rk) :: ts  !bvf,

   real(rk) :: fgrow,fgraze,fmort,fmort2,fmelt,fpond,fpondno3,fpondnh4,fpondsil,fnit,fno3up,fnh4up,fsilup,fskelno3,fskelnh4,fskelsil,ier,dt,Amelt
   _HORIZONTAL_LOOP_BEGIN_
   _GET_HORIZONTAL_(self%id_ia,ia)
   _GET_HORIZONTAL_(self%id_no3,no3)
   _GET_HORIZONTAL_(self%id_nh4,nh4)
   _GET_HORIZONTAL_(self%id_sil,sil)
   _GET_HORIZONTAL_(self%id_temp,temp)

!hayashida
   _GET_HORIZONTAL_(self%id_ts,ts)
   !_GET_HORIZONTAL_(self%id_bvf,bvf)

   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_HORIZONTAL_(self%id_ice_hs,ice_hs)
   _GET_HORIZONTAL_(self%id_par,par)
   _GET_HORIZONTAL_(self%id_airt,airt)
   _GET_(self%id_ph2,ph2)
   _GET_(self%id_no3SW,no3SW)
   _GET_(self%id_nh4SW,nh4SW)
   _GET_(self%id_silSW,silSW) 
   _GET_HORIZONTAL_(self%id_topmelt,topmelt)
   _GET_HORIZONTAL_(self%id_termelt,termelt)
   _GET_HORIZONTAL_(self%id_Amelt,Amelt)
   _GET_HORIZONTAL_(self%id_botmelt,botmelt)
   _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
   _GET_(self%id_u,u)
   _GET_(self%id_v,v)
   _GET_GLOBAL_(self%id_dt,dt)
   if (self%fflush.eq.1) then ! brine fflushing; permeability (5% criterion) has not been impelemented...
    botmelt = botmelt + topmelt
   else if (self%fflush.eq.2) then
    botmelt = botmelt + topmelt + termelt
   end if
   ier = max(botmelt,botgrowth)
   if (botmelt.gt.botgrowth) ier = -ier
!--------------Ice Growth Limitation----------------
   lice = max(0.0,1.0_rk-botmelt/self%crit_melt)
!----------------Light Limitation------------------
   if (par.lt.self%lcompp) then
    llig = 0
   else
    llig = tanh(self%rpp*par)
   end if
!---------------Nutrient Limitation----------------
!  lno3 = no3/(self%ks_no3+no3)
   lno3 = (no3+nh4)/(self%ks_no3+no3+nh4)
   lsil = sil/(self%ks_sil+sil)
!---------------total growth rate----------------
   mum = self%maxg*log(2.0_rk)*exp(self%t_sens*(temp-273.15_rk))
   grow = mum * min(lno3,lsil,llig,lice)   
!--------------heat flux & assoc'd ice loss---------------
   fgrow = grow*ia
   fgraze = self%f_graze*fgrow
   fmort = self%mort*ia*exp(self%t_sens*(temp-273.15_rk))
   fmort2 = self%mort2*ia**2
   !mortenson
   fmelt = min(0.0,ier*ia/self%zia)
   !hayashida
   !fmelt = min(0.0,913./1000.*ier*ia/self%zia) !use hakase  !double check (make sure its consistent -- ratio) 




   fpond = Amelt*self%r_pond*ia/self%zia 
   fpondno3 = Amelt*self%r_pond*no3/self%zia 
   fpondnh4 = Amelt*self%r_pond*nh4/self%zia 
   fpondsil = Amelt*self%r_pond*sil/self%zia 
!---------bio-uptake/depletion of nutrients in skel layer-------!
   fnit     = self%rnit*nh4/(1+par)
!  fno3up   = grow*ia
   fno3up   = (0.2/(0.2+nh4))*(no3/(no3+nh4))*grow*ia
   fnh4up   = grow*ia - fno3up
   fsilup   = grow*self%sil2n*ia
   utaui    = sqrt(self%drag)*sqrt(u**2+v**2)
   hnu      = self%nu/utaui

!mortenson
   fskelno3 = self%md_no3/hnu*(no3SW-no3)/self%zia+min(0.0,ier*no3/self%zia)  !probably this bc more conservative 
   fskelnh4 = self%md_no3/hnu*(nh4SW-nh4)/self%zia+min(0.0,ier*nh4/self%zia) 
   fskelsil = self%md_sil/hnu*(silSW-sil)/self%zia+min(0.0,ier*sil/self%zia)
!hayashida
   !fskelno3 = self%md_no3/hnu*(no3SW-no3)/self%zia
   !fskelnh4 = self%md_no3/hnu*(nh4SW-nh4)/self%zia
   !fskelsil = self%md_sil/hnu*(silSW-sil)/self%zia

   if (self%fmethod.eq.1) then ! ice-ocean flux based on melting/growth.
    fskelno3 = max(0.0,ier*(no3SW-no3)/self%zia)+min(0.0,ier*no3/self%zia)
    fskelnh4 = max(0.0,ier*(nh4SW-nh4)/self%zia)+min(0.0,ier*nh4/self%zia)
    fskelsil = max(0.0,ier*(silSW-sil)/self%zia)+min(0.0,ier*sil/self%zia)
    fmelt    = max(0.0,ier*ph2/self%zia)+min(0.0,ier*ia/self%zia)
   elseif (self%fmethod.eq.2) then ! combination of diffusion and melting/growth
    fskelno3 = fskelno3+max(0.0,ier*(no3SW-no3)/self%zia)
    fskelnh4 = fskelnh4+max(0.0,ier*(nh4SW-nh4)/self%zia)
    fskelsil = fskelsil+max(0.0,ier*(silSW-sil)/self%zia)
   !mortenson
    fmelt    = max(0.0,ier*ph2/self%zia)+min(0.0,ier*ia/self%zia)  !use eric
    !hayashida
    !fmelt    = fmelt+max(0.0,ier*ph2/self%zia)
   endif
   if (ice_hi .lt. 0.1_rk) then ! Ice algal layer is absent.
    fskelno3=0.0_rk
    fskelsil=0.0_rk
    fgrow   =0.0_rk
    fgraze  =0.0_rk
    fmort   =0.0_rk
    fmort2  =0.0_rk
    fmelt   =0.0_rk
    fpond   =0.0_rk
    fpondno3=0.0_rk
    fpondnh4=0.0_rk
    fpondsil=0.0_rk
    fnit    =0.0_rk
   ! _SET_SURFACE_ODE_(self%id_no3,(no3SW-no3)/dt)
    _ADD_SURFACE_SOURCE_(self%id_no3,(no3SW-no3)/dt)
   ! _SET_SURFACE_ODE_(self%id_nh4,(nh4SW-nh4)/dt)
    _ADD_SURFACE_SOURCE_(self%id_nh4,(nh4SW-nh4)/dt)
   ! _SET_SURFACE_ODE_(self%id_sil,(silSW-sil)/dt)
    _ADD_SURFACE_SOURCE_(self%id_sil,(silSW-sil)/dt)
   ! _SET_SURFACE_ODE_(self%id_ia,(ph2-ia)/dt)
    _ADD_SURFACE_SOURCE_(self%id_ia,(ph2-ia)/dt)
   else
!mortenson
 
    !_SET_SURFACE_ODE_(self%id_no3,-fno3up+fskelno3+fnit-fpondno3)
    _ADD_SURFACE_SOURCE_(self%id_no3,-fno3up+fskelno3+fnit-fpondno3)
!   _SET_SURFACE_ODE_(self%id_nh4,fmort+fskelnh4-fnit-fpondnh4)
    !_SET_SURFACE_ODE_(self%id_nh4,-fnh4up+0.3*fmort+fskelnh4-fnit-fpondnh4)
    _ADD_SURFACE_SOURCE_(self%id_nh4,-fnh4up+0.3*fmort+fskelnh4-fnit-fpondnh4)
   ! _SET_SURFACE_ODE_(self%id_sil,-fsilup+fskelsil-fpondsil)
    _ADD_SURFACE_SOURCE_(self%id_sil,-fsilup+fskelsil-fpondsil)

!hayashida
#if 0
   _SET_SURFACE_ODE_(self%id_no3,-fno3up+fskelno3+fnit-fpondno3+min(0.0,913./1000.*ier*no3/self%zia))  !probably this one 
   _SET_SURFACE_ODE_(self%id_nh4,-fnh4up+0.3*fmort+fskelnh4-fnit-fpondnh4+min(0.0,913./1000.*ier*nh4/self%zia))
   _SET_SURFACE_ODE_(self%id_sil,-fsilup+fskelsil-fpondsil+min(0.0,913./1000.*ier*sil/self%zia))
#endif
   
    if (ia.lt.self%ia_b) then
    ! _SET_SURFACE_ODE_(self%id_ia,fgrow)
     _ADD_SURFACE_SOURCE_(self%id_ia,fgrow)
     fgrow=0.0_rk
     fgraze=0.0_rk
     fmort=0.0_rk
     fmort2=0.0_rk
     fpond=0.0_rk
    else
     !_SET_SURFACE_ODE_(self%id_ia,fgrow-fgraze-fmort-fmort2+fmelt-fpond)
     _ADD_SURFACE_SOURCE_(self%id_ia,fgrow-fgraze-fmort-fmort2+fmelt-fpond)
    endif
   endif
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_chl,ia*self%chl2n)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_chlia,ia*self%chl2n*self%zia)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fgrow,fgrow*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fgraze,fgraze*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fmort,fmort*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fmort2,fmort2*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fmelt,fmelt*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fpond,fpond*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fpondno3,fpondno3*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fpondnh4,fpondnh4*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fpondsil,fpondsil*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fnit,fnit*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fskelno3,fskelno3*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fskelnh4,fskelnh4*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fskelsil,fskelsil*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ier,ier*self%spd)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fno3up,fno3up*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fsilup,fsilup*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lice,lice)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_llig,llig)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lno3,lno3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lsil,lsil)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_hnu,hnu) 
   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
end module uvic_icealgae