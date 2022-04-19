!-----------------------------------------------------
!
! University of Victoria DMS(P) model
!
! Originally written for GOTM by Nadja Steiner
! Adapted for FABM by Hakase Hayashida
!
!------------------------------------------------------
#include "fabm_driver.h"

module uvic_dms

   use fabm_types
   use fabm_expressions

   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_dms
! Declare prognostic variables      
      type (type_state_variable_id) :: id_dms,id_dmspd
! Declare diagnostic variables
      type (type_diagnostic_variable_id) :: id_dmspp,id_dmspp1,id_dmspp2,id_fph1lys,id_fph2lys,id_fzo1spd,id_fph1spd,id_fzo2spd,id_fph2spd,id_fspdbac,id_fspdph1,id_fspdph2,id_fph1dms,id_fph2dms,id_fspddms,id_fbacdms,id_fdmsbac,id_fdmspho
      type (type_horizontal_diagnostic_variable_id) :: id_fdmsair
! Declare environmental variables
      type (type_dependency_id) :: id_temp,id_par,id_fph1zo1,id_fph2zo2,id_fph1nh4,id_fph2nh4,id_ph1,id_ph2,id_chl1,id_chl2,id_nf,id_sf,id_fnutph1,id_fnutph2,id_density
! Declare horizontal environmental variables
      type (type_horizontal_dependency_id) :: id_botgrowth,id_botmelt,id_stemp,id_spar,id_patm,id_wind,id_ice_hi,id_fdmspd3ia,id_fdms3ia,id_fspdmel,id_fspdpon,id_fdmsmel,id_fdmspon,id_icespd,id_icedms
! Declare namelist parameters
      real(rk) :: f_ow,dms_0,dmspd_0,q1,q2,yield,k_ly1,k_ly2,f_sl1,f_sl2,f_ex1,f_ex2,k_enz,k_in1,k_in2,h_dmspd,h_dms,k_dmspd,k_dms,k_pho
      integer  :: flux
      real(rk) :: zia
! Declare anything else used in every procedure
      real(rk) :: spd = 86400.0_rk ! Seconds Per Day (spd)
      logical :: use_icedms != .true.

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      
   end type

contains

   subroutine initialize(self,configunit)
   class (type_uvic_dms), intent(inout), target :: self
   integer, intent(in)                          :: configunit
! Declare namelist parameters
   real(rk) :: f_ow,dms_0,dmspd_0,q1,q2,yield,k_ly1,k_ly2,f_sl1,f_sl2,f_ex1,f_ex2,k_enz,k_in1,k_in2,h_dmspd,h_dms,k_dmspd,k_dms,k_pho
   integer  :: flux
   real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,rpi,t_sens,nu,md_no3,md_sil,chl2n,sil2n

!  Initialize parameters to default values
#if 0
!for testing purposes
   f_ow    = 0.0_rk
   dms_0   = 2.0_rk
   dmspd_0 = 5.0_rk
   q1      = 0.45_rk
   q2      = 0.01_rk
   yield   = 0.1_rk
   k_ly1  = 0.05_rk
   k_ly2  = 0.05_rk
   f_sl1   = 1.0_rk
   f_sl2   = 1.0_rk
   f_ex1   = 1.0_rk
   f_ex2   = 1.0_rk
   k_enz   = 0.01_rk
   k_in1   = 0_rk
   k_in2   = 0_rk 
   h_dmspd = 3.0_rk
   h_dms   = 3.0_rk
   k_dmspd = 0_rk 
   k_dms   = 0_rk
   k_pho   = 0_rk
   flux    = 0

!#if 0
   f_ow        = 0.0
   dms_0       = 0.1
   dmspd_0     = 0.1
   q1          = 100 !0.15, !Stefels 2007
   q2          = 9.5 !6 !0.015 Galindo 2014
   yield       = 0.20 !Galindo 2015
   k_ly1       = 0.03
   k_ly2       = 0.03
   f_sl1       = 0.3
   f_sl2       = 0.3
   f_ex1       = 0.05
   f_ex2       = 0.05
   k_enz       = 0.02 
   k_in1       = 0
   k_in2       = 0
   h_dmspd     = 0
   h_dms       = 0
   k_dmspd     = 5
   k_dms       = 0.5
   k_pho       = 0.1
   flux        = 0

   zia = 0.03

!  Register namelist parameters
   self%f_ow    = f_ow
   self%q1      = q1
   self%q2      = q2
   self%yield   = yield
   self%k_ly1  = k_ly1 /self%spd
   self%k_ly2  = k_ly2 /self%spd
   self%f_sl1   = f_sl1
   self%f_sl2   = f_sl2
   self%f_ex1   = f_ex1
   self%f_ex2   = f_ex2
   self%k_enz   = k_enz /self%spd
   self%k_in1   = k_in1 /self%spd
   self%k_in2   = k_in2 /self%spd
   self%h_dmspd = h_dmspd
   self%h_dms   = h_dms
   self%k_dmspd = k_dmspd /self%spd
   self%k_dms   = k_dms /self%spd
   self%k_pho   = k_pho /self%spd
   self%flux    = flux
!if(self%use_icedms) then 
   self%zia = zia
#endif
!----

   call self%get_parameter(self%use_icedms, 'use_icedms', '', 'use use_icedms', default=.false.)
   call self%get_parameter(self%f_ow,'f_ow','','fraction of open water during ice melting (e.g. leads)', default=0.0_rk)
   !call self%get_parameter(dms_0,'dms_0','nmol/L','dms initial value ', default=2.0_rk)
   !call self%get_parameter(dmspd_0,'dmspd_0','nmol/L','dmspd initial value', default=5.0_rk)
   call self%get_parameter(self%q1,'q1','mol-S:mol-N','S:N ratio for ph1', default=0.45_rk)
   call self%get_parameter(self%q2,'q2','mol-S:mol-N','S:N ratio for ph2', default=0.01_rk)
   call self%get_parameter(self%yield,'yield','per day','bacterial lyase rate constant', default=0.1_rk)
   call self%get_parameter(self%k_ly1,'k_ly1','','k_ly1', default=0.05_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%k_ly2,'k_ly2','','k_ly2', default=0.05_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%f_sl1,'f_sl1','','f_sl1', default=1.0_rk)
   call self%get_parameter(self%f_sl2,'f_sl2','','f_sl2', default=1.0_rk)
   call self%get_parameter(self%f_ex1,'f_ex1','','f_ex1', default=1.0_rk)
   call self%get_parameter(self%f_ex2,'f_ex2','','f_ex2', default=1.0_rk)
   call self%get_parameter(self%k_enz,'k_enz','d-1','rate constant for enzymatic cleavage (free lyase)', default=0.01_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%k_in1,'k_in1','','k_in1', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%k_in2,'k_in2','','k_in2', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%h_dmspd,'h_dmspd','','half-saturation constant for bacterial dmspd uptake', default=3.0_rk)
   call self%get_parameter(self%h_dms,'h_dms','','half-saturation constant for bacterial dms uptake', default=3.0_rk)
   call self%get_parameter(self%k_dmspd,'k_dmspd','','k_dmspd', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%k_dms,'k_dms','','k_dms', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%k_pho,'k_pho','','k_pho', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%flux,'flux','','air-sea gas transfer velocity parameterization', default=0)
   ! read in icealgae model vars 
   if (self%use_icedms) then
      call self%get_parameter(self%zia, 'zia','m', 'ice algal layer thickness', default=0.03_rk) ! zia = 0.03_rk
   endif

! Register prognostic variables
      call self%register_state_variable(self%id_dms,'dms','nmol/L','Dimethylsulfide',minimum=0.0_rk) !initial_value=dms_0
      call self%register_state_variable(self%id_dmspd,'dmspd','nmol/L','Dissolved dimethylsulfoniopropionate',minimum=0.0_rk) !initial_value=dmspd_0
! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_dmspp,'dmspp','nmol/L','Particulate dimethylsulfoniopropionate in ph1 and ph2')
      call self%register_diagnostic_variable(self%id_dmspp1,'dmspp1','nmol/L','Particulate dimethylsulfoniopropionate in ph1')
      call self%register_diagnostic_variable(self%id_dmspp2,'dmspp2','nmol/L','Particulate dimethylsulfoniopropionate in ph2')
      call self%register_diagnostic_variable(self%id_fph1lys,'fph1lys','nM d-1','ph1 cell lysis')
      call self%register_diagnostic_variable(self%id_fph2lys,'fph2lys','nM d-1','ph2 cell lysis')
      call self%register_diagnostic_variable(self%id_fzo1spd,'fzo1spd','nM d-1','zo1 sloppy feeding')
      call self%register_diagnostic_variable(self%id_fph1spd,'fph1spd','nM d-1','ph1 exudation/cell lysis')
      call self%register_diagnostic_variable(self%id_fzo2spd,'fzo2spd','nM d-1','zo2 sloppy feeding')
      call self%register_diagnostic_variable(self%id_fph2spd,'fph2spd','nM d-1','ph2 exudation/cell lysis')
      call self%register_diagnostic_variable(self%id_fspdbac,'fspdbac','nM d-1','DMSPd bacterial consumption')
      call self%register_diagnostic_variable(self%id_fspdph1,'fspdph1','nM d-1','ph1 extracellular lyase')
      call self%register_diagnostic_variable(self%id_fspdph2,'fspdph2','nM d-1','ph2 extracellular lyase')
      call self%register_diagnostic_variable(self%id_fph1dms,'fph1dms','nM d-1','ph1 intracellular lyase')
      call self%register_diagnostic_variable(self%id_fph2dms,'fph2dms','nM d-1','ph2 intracellular lyase')
      call self%register_diagnostic_variable(self%id_fspddms,'fspddms','nM d-1','enzymatic cleavage (free lyase)')
      call self%register_diagnostic_variable(self%id_fbacdms,'fbacdms','nM d-1','bacterial lyase (yield)')
      call self%register_diagnostic_variable(self%id_fdmsbac,'fdmsbac','nM d-1','DMS bacterial consumption')
      call self%register_diagnostic_variable(self%id_fdmspho,'fdmspho','nM d-1','photolysis')
      call self%register_horizontal_diagnostic_variable(self%id_fdmsair,'fdmsair','umol m-2 d-1','Sea-to-air DMS flux',source=source_do_horizontal)
! Register environmental variables
      call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt)
      call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_patm,standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_spar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_wind,standard_variables%wind_speed)
      call self%register_dependency(self%id_ice_hi,standard_variables%sea_ice_thickness)
      call self%register_dependency(self%id_stemp,'uvic_eco_stemp','','')
      call self%register_dependency(self%id_ph1,'uvic_eco_ph1','','')
      call self%register_dependency(self%id_ph2,'uvic_eco_ph2','','')
      call self%register_dependency(self%id_chl1,'uvic_eco_chl1','','')
      call self%register_dependency(self%id_chl2,'uvic_eco_chl2','','')
      call self%register_dependency(self%id_nf,'uvic_eco_nf','','')
      call self%register_dependency(self%id_sf,'uvic_eco_sf','','')
      call self%register_dependency(self%id_fph1zo1,'uvic_eco_fph1zo1','','')
      call self%register_dependency(self%id_fph1nh4,'uvic_eco_fph1nh4','','')
      call self%register_dependency(self%id_fph2zo2,'uvic_eco_fph2zo2','','')
      call self%register_dependency(self%id_fph2nh4,'uvic_eco_fph2nh4','','')
      call self%register_dependency(self%id_fnutph1,'uvic_eco_fnutph1','','')
      call self%register_dependency(self%id_fnutph2,'uvic_eco_fnutph2','','')
      call self%request_coupling(self%id_stemp,'uvic_eco_stemp')
      call self%request_coupling(self%id_ph1,'uvic_eco_ph1')
      call self%request_coupling(self%id_ph2,'uvic_eco_ph2')
      call self%request_coupling(self%id_chl1,'uvic_eco_chl1')
      call self%request_coupling(self%id_chl2,'uvic_eco_chl2')
      call self%request_coupling(self%id_nf,'uvic_eco_nf')
      call self%request_coupling(self%id_sf,'uvic_eco_sf')
      call self%request_coupling(self%id_fph1zo1,'uvic_eco_fph1zo1')
      call self%request_coupling(self%id_fph1nh4,'uvic_eco_fph1nh4')
      call self%request_coupling(self%id_fph2zo2,'uvic_eco_fph2zo2')
      call self%request_coupling(self%id_fph2nh4,'uvic_eco_fph2nh4')
      call self%request_coupling(self%id_fnutph1,'uvic_eco_fnutph1')
      call self%request_coupling(self%id_fnutph2,'uvic_eco_fnutph2')
 
      if(self%use_icedms) then 
       call self%register_dependency(self%id_fdmspd3ia,'uvic_icedms_fspdaiw','','')  !jpnote: runs if are commented ... ? 
       call self%register_dependency(self%id_fdms3ia,'uvic_icedms_fdmsaiw','','')
      ! call self%register_dependency(self%id_fspdmel,'uvic_icedms_fspdmel','','')   !this
      ! call self%register_dependency(self%id_fspdpon,'uvic_icedms_fspdpon','','')   !this
      ! call self%register_dependency(self%id_fdmsmel,'uvic_icedms_fdmsmel','','') !this
      ! call self%register_dependency(self%id_fdmspon,'uvic_icedms_fdmspon','','')  !this
       call self%request_coupling(self%id_fdmspd3ia,'uvic_icedms_fspdaiw')
       call self%request_coupling(self%id_fdms3ia,'uvic_icedms_fdmsaiw')
       !call self%request_coupling(self%id_fspdmel,'uvic_icedms_fspdmel')  !this
       !call self%request_coupling(self%id_fspdpon,'uvic_icedms_fspdpon')  !this
       !call self%request_coupling(self%id_fdmsmel,'uvic_icedms_fdmsmel') !this
       !call self%request_coupling(self%id_fdmspon,'uvic_icedms_fdmspon') !this
       call self%register_dependency(self%id_icespd,'uvic_icedms_dmspd','','')
       call self%request_coupling(self%id_icespd,'uvic_icedms_dmspd')
       call self%register_dependency(self%id_icedms,'uvic_icedms_dms','','')
       call self%request_coupling(self%id_icedms,'uvic_icedms_dms')  
       call self%register_dependency(self%id_density,standard_variables%density)
      endif
   end subroutine initialize
   
   subroutine do(self,_ARGUMENTS_DO_)

   class (type_uvic_dms),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
! Declare prognostic variables
   real(rk) :: dms,dmspd,dmspp1,dmspp2
! Declare diagnostic variables   
   real(rk) :: fph1lys,fph2lys,fzo1spd,fph1spd,fzo2spd,fph2spd,fspdbac,fspdph1,fspdph2,fph1dms,fph2dms,fspddms,fbacdms,fdmsbac,fdmspho
! Declare environmental parameters
   real(rk) :: ice_hi,temp,par,ph1,ph2,chl1,chl2,nf,sf,fph1zo1,fph2zo2,fph1nh4,fph2nh4,stemp,spar,fnutph1,fnutph2
! Declare anything else used in this subroutine
! molecular weight of N ( for mmol-N to mg-N conversion)
!   real(rk), parameter       :: wN = 14.0_rk
! conversion factor 1mg-S/m3=31.25nmol/L DMS or DMSP
!   real(rk), parameter       :: S2mol = 31.25_rk
   
   _LOOP_BEGIN_

! Retrieve prognostic variables
   _GET_(self%id_dms,dms)
   _GET_(self%id_dmspd,dmspd)
! Retrieve environmental variables
   _GET_(self%id_temp,temp)
   _GET_HORIZONTAL_(self%id_stemp,stemp)
   _GET_(self%id_par,par)
   _GET_HORIZONTAL_(self%id_spar,spar)
   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_(self%id_chl1,chl1)
   _GET_(self%id_chl2,chl2)
   _GET_(self%id_ph1,ph1)
   _GET_(self%id_ph2,ph2) 
   _GET_(self%id_nf,nf)
   _GET_(self%id_sf,sf)
   _GET_(self%id_fph1zo1,fph1zo1)
   _GET_(self%id_fph2zo2,fph2zo2)
   _GET_(self%id_fph1nh4,fph1nh4)
   _GET_(self%id_fph2nh4,fph2nh4)
   _GET_(self%id_fnutph1,fnutph1)
   _GET_(self%id_fnutph2,fnutph2)
   fph1zo1 = fph1zo1 / self%spd
   fph2zo2 = fph2zo2 / self%spd
   fph1nh4 = fph1nh4 / self%spd
   fph2nh4 = fph2nh4 / self%spd
   fnutph1 = fnutph1 / self%spd
   fnutph2 = fnutph2 / self%spd
   dmspp1 = self%q1*chl1
   dmspp2 = self%q2*chl2
! DMS and DMSPd budget components
   fph1lys = self%k_ly1*self%q1*chl1*(1/(nf+0.1))
   fph2lys = self%k_ly2*self%q2*chl2*(1/(min(nf,sf)+0.1))
   fzo1spd = self%f_sl1*self%q1*fph1zo1*(chl1/ph1)
   fzo2spd = self%f_sl2*self%q2*fph2zo2*(chl2/ph2)
   fph1spd = self%q1*fnutph1/ph1*chl1*(self%f_ex1+(1-self%f_ex1)*(1-nf))
   fph2spd = self%q2*fnutph2/ph2*chl2*(self%f_ex2+(1-self%f_ex2)*(1-min(nf,sf)))
   fspdbac = self%k_dmspd*dmspd
!  fspdbac = self%k_dmspd*(dmspd/(dmspd+self%h_dmspd))*dmspd
!  fspdbac = self%k_dmspd*10./(par+10.)*(dmspd/(dmspd+self%h_dmspd))*dmspd
   fph1dms = self%k_in1*self%q1*dmspp1
   fph2dms = self%k_in2*self%q2*dmspp2
   fspddms = self%k_enz*dmspd
   fbacdms = self%yield*fspdbac
   fdmsbac = self%k_dms*dms
!  fdmsbac = self%k_dms*(dms/(dms+self%h_dms))*dms
!  fdmsbac = self%k_dms*10./(par+10.)*(dms/(dms+self%h_dms))*dms
   fdmspho = self%k_pho*par/(par+1)*dms
! Compute and save prognostic variables
   _SET_ODE_(self%id_dms,fph1dms+fph2dms+fspdph1+fspdph2+fspddms+fbacdms-fdmsbac-fdmspho)
   _SET_ODE_(self%id_dmspd,fph1lys+fph2lys+fzo1spd+fph1spd+fzo2spd+fph2spd-fspddms-fspdbac-fspdph1-fspdph2)
! Save diagnostic variables
   _SET_DIAGNOSTIC_(self%id_fph1dms,fph1dms*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2dms,fph2dms*self%spd)
   _SET_DIAGNOSTIC_(self%id_fspdph1,fspdph1*self%spd)
   _SET_DIAGNOSTIC_(self%id_fspdph2,fspdph2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fbacdms,fbacdms*self%spd)
   _SET_DIAGNOSTIC_(self%id_fspddms,fspddms*self%spd)
   _SET_DIAGNOSTIC_(self%id_fdmsbac,fdmsbac*self%spd)
   _SET_DIAGNOSTIC_(self%id_fdmspho,fdmspho*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph1lys,fph1lys*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2lys,fph2lys*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo1spd,fzo1spd*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph1spd,fph1spd*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo2spd,fzo2spd*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2spd,fph2spd*self%spd)
   _SET_DIAGNOSTIC_(self%id_fspdbac,fspdbac*self%spd)
   _SET_DIAGNOSTIC_(self%id_dmspp,dmspp1+dmspp2)
   _SET_DIAGNOSTIC_(self%id_dmspp1,dmspp1)
   _SET_DIAGNOSTIC_(self%id_dmspp2,dmspp2)

   _LOOP_END_

   end subroutine do

! Surface Flux
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_dms),intent(in) :: self
   real (rk) :: botgrowth,botmelt,dms,stemp,patm,dmsair,wind,ice_hi,tracd_fascm
   real (rk) :: adms(4),ddms(2)
   real (rk) :: k_600,v_dms,henry,pres_atm
   real (rk) :: osc,sc,rkb,rk0,W
   real (rk) :: fdmspd3ia,fdms3ia,fspdmel,fdmsmel,fspdpon,fdmspon,icespd,icedms,density,dmspd
   data adms/2674.0,147.12,3.726,0.038/
! Saltzman provides Ostwald solubility, need to be transferred into Bunsen coefficients for consistency
   data ddms/-10.1794,3761.33/

   _DECLARE_ARGUMENTS_DO_SURFACE_

   _HORIZONTAL_LOOP_BEGIN_
   
   _GET_(self%id_dms,dms)
   _GET_HORIZONTAL_(self%id_stemp,stemp)
   _GET_HORIZONTAL_(self%id_wind,wind)
   _GET_HORIZONTAL_(self%id_patm,patm)
   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_HORIZONTAL_(self%id_botmelt,botmelt)
   _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
   pres_atm = patm*9.87e-6 ! converting from [Pa] to [atm]
   pres_atm = 10000*9.87e-6 ! converting from [Pa] to [atm]
! Question: treatment of DMS_atm (0 for now)
   tracd_fascm = 0.0_rk
   sc = adms(1)-adms(2)*stemp+adms(3)*stemp**2-adms(4)*stemp**3
   select case (self%flux)
    case (0) !if (self%flux.eq.0) then
     k_600 = 0.222_rk*wind**2+0.333_rk*wind
    case (1) !else if (self%flux.eq.1) then
     if (wind.lt.3.6_rk) then
      rk0 = 0.17_rk*wind
     else if (wind.ge.3.6_rk) then
      rk0 = 2.85_rk*wind-9.65_rk
     end if
    !rk0 = 1.57e-4*wind/sqrt(sc/(cd*600.)) !Jaehne, 1987 ? woolf 1990
     W = 3.84e-6*wind**3.41_rk  !Monahan & O'Muircheartaigh 1980
    !W = 2.98e-7*wind**4.04  !Zhao & Toba, 2001
     osc = exp(ddms(1) + ddms(2)/(stemp+273.15_rk))/pres_atm ! Ostwald solublity coefficient based on Saltzman 1993
     rkb = 2450_rk*W/(osc*(1.0_rk+(14_rk*osc*sqrt(sc))**(-1.0_rk/1.2_rk))**1.2_rk)
     k_600 = rk0 + rkb
   case default
     k_600 = 0.0_rk
   end select !end if
   v_dms = 0.24_rk*k_600/sqrt(sc/600._rk) ! convert from [cm/h] to [m/d]
   henry = exp(12.64_rk-3547_rk/(stemp+273.15_rk))/(0.08205_rk*(stemp+273.15_rk))
   dmsair = (dms-tracd_fascm/henry)*v_dms/self%spd !sea-to-air flux
   if(ice_hi.gt.0.0_rk) then 
!    if(botmelt.gt.0.0_rk) then
     dmsair=self%f_ow*dmsair
!    else
!     dmsair=0.0_rk
!    end if
   end if
   _SET_SURFACE_EXCHANGE_(self%id_dms,-dmsair)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmsair,dmsair*self%spd)
 
   if(self%use_icedms) then  !jpnote
    _GET_HORIZONTAL_(self%id_fdmspd3ia,fdmspd3ia)
    _GET_HORIZONTAL_(self%id_fdms3ia,fdms3ia)
    _GET_HORIZONTAL_(self%id_fspdmel,fspdmel)
    _GET_HORIZONTAL_(self%id_fspdpon,fspdpon)
    _GET_HORIZONTAL_(self%id_fdmsmel,fdmsmel)
    _GET_HORIZONTAL_(self%id_fdmspon,fdmspon)
    _GET_HORIZONTAL_(self%id_icespd,icespd)
    _GET_HORIZONTAL_(self%id_icedms,icedms)
    _GET_(self%id_dmspd,dmspd)
    _GET_(self%id_density,density)
!HH0: concentration/dilution effect
!    _SET_SURFACE_EXCHANGE_(self%id_dmspd,(-fspdmel+fspdpon)/self%spd*self%zia-fdmspd3ia/self%spd*self%zia)
!    _SET_SURFACE_EXCHANGE_(self%id_dms,(-fdmsmel+fdmspon)/self%spd*self%zia-fdms3ia/self%spd*self%zia)
    _SET_SURFACE_EXCHANGE_(self%id_dmspd,913./density*botgrowth*(dmspd)+(1000./density*icespd-dmspd)*(-913./density*fspdmel+1000./density*fspdpon)/icespd/self%spd*self%zia-fdmspd3ia/self%spd*self%zia)
    _SET_SURFACE_EXCHANGE_(self%id_dms,913./density*botgrowth*(dms)+(1000./density*icedms-dms)*(-913./density*fdmsmel+1000./density*fdmspon)/icedms/self%spd*self%zia-fdms3ia/self%spd*self%zia)
!HH1
   end if
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module uvic_dms
