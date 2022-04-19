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
      type (type_surface_state_variable_id) :: id_no3, id_sil,id_ia,id_nh4
! Declare horizontal diagnostic variables
      type (type_horizontal_diagnostic_variable_id) :: id_chl,id_chlia,id_fgrow,id_fgraze,id_fmort,id_fmort2,id_fmelt,id_fpond,id_fpondno3,id_fpondnh4,id_fpondsil,id_fno3up,id_fsilup,id_fskelno3,id_fskelnh4,id_fskelsil,id_ier,id_lice,id_llig,id_lno3,id_lsil,id_hnu,id_fnit
! Declare environmental variables
      type (type_dependency_id) :: id_ph2,id_no3SW,id_nh4SW,id_silSW,id_u,id_v
      type (type_horizontal_dependency_id) :: id_airt,id_ice_hs,id_par,id_temp,id_ice_hi,id_botmelt,id_botgrowth,id_topmelt,id_termelt,id_Amelt,id_ts,id_bvf
      type (type_global_dependency_id) :: id_dt
! Declare namelist parameters
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
! Declare namelist parameters
   real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,rnit,skno3_0,sknh4_0,sksil_0,ia_0,ia_b,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,t_sens,nu,md_no3,md_sil,chl2n,sil2n
! Define the namelist
  ! namelist /uvic_icealgae/ r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,t_sens,nu,md_no3,md_sil,chl2n,sil2n
!  Initialize parameters to default values.
   r_pond    = 0.0175_rk
   fmethod   = 0.0_rk
   fflush     = 0
   drag      = 0.005_rk
   f_graze   = 0.1_rk
   zia       = 0.03_rk
   ac_ia     = 0.03_rk
  ! ia_0      = 0.1_rk
   ia_b      = 0.16_rk
   rnit      = 0.01_rk
  ! skno3_0   = 2.0_rk 
 !  sknh4_0   = 0.01_rk 
  ! sksil_0   = 5.0_rk
   ks_no3    = 1.0_rk
   ks_sil    = 4.0_rk
   maxg      = 0.8511_rk
   mort      = 0.05_rk
   mort2     = 0.05_rk
   crit_melt = 0.015_rk
   lcompp    = 0.4_rk
   rpp       = 0.1_rk
   t_sens    = 0.0633_rk
   nu        = 1.86e-6_rk
   md_no3    = 0.47e-9_rk
   md_sil    = 0.47e-9_rk
   chl2n     = 2.8_rk
   sil2n     = 1.7_rk

   r_pond    = 0.0175
   fmethod   = 0
   fflush    = 0
   drag      = 0.0054
   f_graze   = 0.0 !L2005
   zia       = 0.03
   ac_ia     = 0.017 !mcdonald2015
   rnit      = 0.01
 !  ia_0      = 1.0  !321.0, ! ca. 900 mg-chl m-3
   ia_b      = 1.0  !0.03-0.16 based on Michel's reply.
  ! skno3_0   = 7.2
  ! sknh4_0   = 0.01
  ! sksil_0   = 14.7
   ks_no3    = 1.0
   ks_sil    = 4.0
   maxg      = 0.85
   mort      = 0.03
   mort2     = 0.00015
   crit_melt = 0.015
   lcompp    = 0.0 ! 1.0 !mock&gradinger99
   rpp       = 2.0
   t_sens    = 0.0633
   nu        = 1.85e-6
   md_no3    = 0.47e-9
   md_sil    = 0.47e-9
   chl2n     = 3.533 !Lavoie 2.80134,
   sil2n     = 1.7

#if 0
!  Read namelist parameters
   read(configunit,nml=uvic_icealgae)
!  Register namelist parameters
#endif 
   self%ac_ia = ac_ia
   self%rnit  = rnit/self%spd
   self%ia_b = ia_b
   self%f_graze = f_graze
   self%r_pond = r_pond/self%spd
   self%fmethod = fmethod
   self%fflush   = fflush
   self%drag = drag
   self%zia = zia
   self%ks_no3 =ks_no3
   self%ks_sil =ks_sil
   self%maxg = maxg/ self%spd
   self%mort = mort/ self%spd
   self%mort2 = mort2/ self%spd
   self%crit_melt = crit_melt / self%spd
   self%lcompp = lcompp
   self%rpp = rpp
   self%t_sens = t_sens
   self%nu     = nu
   self%md_no3 = md_no3
   self%md_sil = md_sil
   self%chl2n  = chl2n
   self%sil2n  = sil2n
! Register prognostic variables
      call self%register_state_variable(self%id_no3,'no3','mmol m-3','skel. NO_3',minimum=0.0_rk)  !initial_value=skno3_0
      call self%register_state_variable(self%id_sil,'sil','mmol m-3','skel. Si',minimum=0.0_rk) !initial_value=sksil_0
      call self%register_state_variable(self%id_ia,'ia','mmol m-3','Ice algae',minimum=0.0_rk)  !initial_value=ia_0
      call self%register_state_variable(self%id_nh4,'nh4','mmol m-3','NH4',minimum=0.0_rk)  !initial_value=sknh4_0
! Register diagnostic variables
      call self%register_horizontal_diagnostic_variable(self%id_chl,'chl','mg m-3','Ice algae in per cubic meter',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_chlia,'chlia','mg m-2','Ice algae in per square meter',source=source_do_horizontal) !jpnote 
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
      call self%register_horizontal_dependency(self%id_ts,standard_variables%surface_ice_temperature)
     ! call self%register_horizontal_dependency(self%id_bvf,standard_variables%brine_volume_fraction)
      call self%register_horizontal_dependency(self%id_ice_hi,standard_variables%sea_ice_thickness)
      call self%register_horizontal_dependency(self%id_ice_hs,standard_variables%snow_thickness)
       call self%register_horizontal_dependency(self%id_par,standard_variables%lowest_ice_layer_PAR)
      call self%register_horizontal_dependency(self%id_airt,standard_variables%surface_temperature)
      call self%register_horizontal_dependency(self%id_topmelt,standard_variables%topmelt)
      call self%register_horizontal_dependency(self%id_termelt,standard_variables%termelt)
      call self%register_horizontal_dependency(self%id_Amelt,standard_variables%f_melt)
      call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt)
      call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)
      call self%register_dependency(self%id_u,standard_variables%zonal_current)
      call self%register_dependency(self%id_v,standard_variables%meridional_current)
      call self%register_global_dependency(self%id_dt,standard_variables%timestep)
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
   real(rk) :: bvf,ts,fgrow,fgraze,fmort,fmort2,fmelt,fpond,fpondno3,fpondnh4,fpondsil,fnit,fno3up,fnh4up,fsilup,fskelno3,fskelnh4,fskelsil,ier,dt,Amelt
   _HORIZONTAL_LOOP_BEGIN_
   _GET_HORIZONTAL_(self%id_ia,ia)
   _GET_HORIZONTAL_(self%id_no3,no3)
   _GET_HORIZONTAL_(self%id_nh4,nh4)
   _GET_HORIZONTAL_(self%id_sil,sil)
   _GET_HORIZONTAL_(self%id_temp,temp)
   _GET_HORIZONTAL_(self%id_ts,ts)
  ! _GET_HORIZONTAL_(self%id_bvf,bvf)
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
   fmelt = min(0.0,913./1000.*ier*ia/self%zia)
   fpond = Amelt*self%r_pond*ia/self%zia 
   fpondno3 = Amelt*self%r_pond*no3/self%zia 
   fpondnh4 = Amelt*self%r_pond*nh4/self%zia 
   fpondsil = Amelt*self%r_pond*sil/self%zia 
!   if (bvf .gt. 5) then !Brine drainage
!    fpond = fpond + ia/20./86400. !brine drainage rate constant (20 days)
!    fpondno3 = fpondno3 + no3/20./86400.
!    fpondnh4 = fpondnh4 + nh4/20./86400.
!    fpondsil = fpondsil + sil/20./86400.
!   endif

!---------bio-uptake/depletion of nutrients in skel layer-------!
   fnit     = self%rnit*nh4/(1+par)
!  fno3up   = grow*ia
   fno3up   = (0.2/(0.2+nh4))*(no3/(no3+nh4))*grow*ia
   fnh4up   = grow*ia - fno3up
   fsilup   = grow*self%sil2n*ia
   utaui    = sqrt(self%drag)*sqrt(u**2+v**2)
   hnu      = self%nu/utaui
   fskelno3 = self%md_no3/hnu*(no3SW-no3)/self%zia
   fskelnh4 = self%md_no3/hnu*(nh4SW-nh4)/self%zia
   fskelsil = self%md_sil/hnu*(silSW-sil)/self%zia
   if (self%fmethod.eq.1) then ! ice-ocean flux based on melting/growth.
    fskelno3 = max(0.0,ier*(no3SW-no3)/self%zia)+min(0.0,ier*no3/self%zia)
    fskelnh4 = max(0.0,ier*(nh4SW-nh4)/self%zia)+min(0.0,ier*nh4/self%zia)
    fskelsil = max(0.0,ier*(silSW-sil)/self%zia)+min(0.0,ier*sil/self%zia)
    fmelt    = max(0.0,ier*ph2/self%zia)+min(0.0,ier*ia/self%zia)
   elseif (self%fmethod.eq.2) then ! combination of diffusion and melting/growth
    fskelno3 = fskelno3+max(0.0,ier*(no3SW-no3)/self%zia)
    fskelnh4 = fskelnh4+max(0.0,ier*(nh4SW-nh4)/self%zia)
    fskelsil = fskelsil+max(0.0,ier*(silSW-sil)/self%zia)
    fmelt    = fmelt+max(0.0,ier*ph2/self%zia)
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
    _SET_SURFACE_ODE_(self%id_no3,(no3SW-no3)/dt)
    _SET_SURFACE_ODE_(self%id_nh4,(nh4SW-nh4)/dt)
    _SET_SURFACE_ODE_(self%id_sil,(silSW-sil)/dt)
    _SET_SURFACE_ODE_(self%id_ia,(ph2-ia)/dt)
   else
    _SET_SURFACE_ODE_(self%id_no3,-fno3up+fskelno3+fnit-fpondno3+min(0.0,913./1000.*ier*no3/self%zia))
    _SET_SURFACE_ODE_(self%id_nh4,-fnh4up+0.3*fmort+fskelnh4-fnit-fpondnh4+min(0.0,913./1000.*ier*nh4/self%zia))
    _SET_SURFACE_ODE_(self%id_sil,-fsilup+fskelsil-fpondsil+min(0.0,913./1000.*ier*sil/self%zia))
    if (ia.lt.self%ia_b) then
     _SET_SURFACE_ODE_(self%id_ia,fgrow)
     fgrow=0.0_rk
     fgraze=0.0_rk
     fmort=0.0_rk
     fmort2=0.0_rk
     fpond=0.0_rk
    else
     _SET_SURFACE_ODE_(self%id_ia,fgrow-fgraze-fmort-fmort2+fmelt-fpond)
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