!-----------------------------------------------------
!
! University of Victoria DMS model in sea ice
!
! Contributors: Hakase Hayashida, Nadja Steiner, Eric Mortenson, Adam Mohanan
!
! Model description: The model will be described in a paper.
!
!------------------------------------------------------

#include "fabm_driver.h"

module uvic_icedms

   use fabm_types
   use fabm_expressions
   use fabm_standard_variables

   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_icedms
! Declare horizontal prognostic variables      
      type (type_surface_state_variable_id) :: id_dmspd,id_dms
! Declare horizontal diagnostic variables
      type (type_horizontal_diagnostic_variable_id) :: id_dmspp,id_flysspd,id_fexuspd,id_fspdbac,id_fspddms,id_fspdaiw,id_fsppdms,id_fbacdms,id_fdmsbac,id_fdmspho,id_fdmsaiw,id_fspdpon,id_fdmspon,id_fspdmel,id_fdmsmel
! Declare environmental variables
      type (type_dependency_id) :: id_dmspdSW,id_dmsSW,id_temp
      type (type_horizontal_dependency_id) :: id_botgrowth,id_bvf,id_ts,id_tb,id_par,id_ice_hi,id_fpond,id_fmelt,id_fgrow,id_fmort,id_fgraze,id_hnu,id_ia,id_chl,id_lno3,id_lsil,id_lice,id_llig
      type (type_global_dependency_id) :: id_dt
! Declare namelist parameters
      real(rk) :: dmsi_0,dmspdi_0,qi,h_dmspdi,h_dmsi,k_dmspdi,k_dmsi,yieldi,k_lyi,f_exi,k_exi,k_ini,k_phoi,zia
! Declare anything else used in all procedures
      real(rk) :: spd = 86400.0_rk ! Seconds Per Day (spd) 

      contains

      procedure :: initialize
      procedure :: do_surface
      
   end type

contains

   subroutine initialize(self,configunit)
   class (type_uvic_icedms), intent(inout), target :: self
   integer, intent(in)                          :: configunit
! Declare yaml parameters
   real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,rpi,t_sens,nu,md_no3,md_sil,chl2n,sil2n
   real(rk) :: dmsi_0,dmspdi_0,qi,h_dmspdi,h_dmsi,k_dmspdi,k_dmsi,yieldi,k_lyi,f_exi,k_exi,k_ini,k_phoi
! Define the yaml parameters


#if 0
!!  Initialize parameters to default values.
   dmsi_0   = 1.0_rk 
   dmspdi_0 = 1.0_rk
   qi       = 9.43_rk
   h_dmspdi = 3.0_rk
   h_dmsi   = 3.0_rk
   k_dmspdi = 5.7_rk
   k_dmsi   = 5.7_rk
   yieldi   = 0.03_rk
   k_lyi    = 1_rk
   f_exi    = 1_rk
   k_exi    = 0_rk
   k_ini    = 0_rk
   k_phoi   = 0.01_rk
   zia      = 0.03_rk
  
  
   !dmsi_0   = 0.1
   !dmspdi_0 = 0.1
   qi  = 9.5 !Galindo 2014
   h_dmspdi = 0
   h_dmsi   = 0
   k_dmspdi = 1 !Galindo 2015
   k_dmsi   = 0.2
   yieldi   = 0.2 !Galindo 2015
   f_exi    = 0.05
   k_lyi    = 0.03
   k_exi    = 0.02
   k_ini    = 0
   k_phoi   = 0.1 ! Taalba 2012

!  Register namelist parameters
   self%qi       = qi
   self%h_dmspdi = h_dmspdi
   self%h_dmsi   = h_dmsi
   self%k_dmspdi = k_dmspdi / self%spd
   self%k_dmsi   = k_dmsi /self%spd
   self%yieldi   = yieldi
   self%k_lyi    = k_lyi /self%spd
   self%f_exi    = f_exi
   self%k_exi    = k_exi /self%spd
   self%k_ini    = k_ini /self%spd
   self%k_phoi   = k_phoi /self%spd
   self%zia      = zia
#endif

   call self%get_parameter(self%qi, 'qi','nmol-S:ug-chla', 'DMSPp:Chl-a ratio for ice algae', default=9.43_rk)
   call self%get_parameter(self%h_dmspdi, 'h_dmspdi','', 'half-saturation constant for bacterial dmspd uptake', default=3.0_rk)
   call self%get_parameter(self%h_dmsi, 'h_dmsi','', 'half-saturation constant for bacterial dms uptake', default=3.0_rk)
   call self%get_parameter(self%k_dmspdi, 'k_dmspdi','per day', 'dmspd loss rate constant', default=5.7_rk,scale_factor=1.0_rk/self%spd)
   !self%k_dmspdi = k_dmspdi / self%spd
   call self%get_parameter(self%k_dmsi, 'k_dmsi','per day', 'dms loss rate constant', default=5.7_rk,scale_factor=1.0_rk/self%spd)
   !self%k_dmsi   = k_dmsi /self%spd
   call self%get_parameter(self%yieldi, 'yieldi','', 'dms yield', default=0.03_rk)
   call self%get_parameter(self%k_lyi, 'k_lyi','', 'fraction of sloppy feeding', default=1.0_rk,scale_factor=1.0_rk/self%spd)
   !self%k_lyi    = k_lyi /self%spd
   call self%get_parameter(self%f_exi, 'f_exi','', 'fraction of exludation/cell lysis', default=1.0_rk)
   call self%get_parameter(self%k_exi, 'k_exi','', 'extracellular lyase rate constant', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   !self%k_exi    = k_exi /self%spd
   call self%get_parameter(self%k_ini, 'k_ini','', 'intracellular lyase rate constant', default=0.0_rk,scale_factor=1.0_rk/self%spd)
   !self%k_ini    = k_ini /self%spd
   call self%get_parameter(self%k_phoi, 'k_phoi','', 'photolysis rate constant', default=0.01_rk,scale_factor=1.0_rk/self%spd)
   !self%k_phoi   = k_phoi /self%spd
   call self%get_parameter(self%zia, 'zia','m', 'ice algal layer thickness ', default=0.03_rk)

#if 0
      ! icedms vars from yaml 
      call self%get_parameter(self%qi, 'qi','nmol-S:ug-chla', 'DMSPp:Chl-a ratio for ice algae', default=9.43_rk)
      call self%get_parameter(self%h_dmspdi, 'h_dmspdi','', 'half-saturation constant for bacterial dmspd uptake', default=3.0_rk)
      call self%get_parameter(self%h_dmsi, 'h_dmsi','', 'half-saturation constant for bacterial dms uptake', default=3.0_rk)
      call self%get_parameter(self%k_dmspdi, 'k_dmspdi','per day', 'dmspd loss rate constant', default=5.7_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%k_dmsi, 'k_dmsi','per day', 'dms loss rate constant', default=5.7_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%yieldi, 'yieldi','', 'dms yield', default=0.03_rk)
      call self%get_parameter(self%k_lyi, 'k_lyi','', 'fraction of sloppy feeding', default=1.0_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%f_exi, 'f_exi','', 'fraction of exludation/cell lysis', default=1.0_rk)
      call self%get_parameter(self%k_exi, 'k_exi','', 'extracellular lyase rate constant', default=0.0_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%k_ini, 'k_ini','', 'intracellular lyase rate constant', default=0.0_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%k_phoi, 'k_phoi','', 'photolysis rate constant', default=0.01_rk,scale_factor=1.0_rk/self%spd)
      call self%get_parameter(self%zia, 'zia','m', 'ice algal layer thickness ', default=0.03_rk)
#endif
! Register prognostic variables
      call self%register_state_variable(self%id_dms,'dms','nmol/L','ice DMS',minimum=0.0_rk)  !initial_value=dmsi_0, #initial value read from yaml 
      call self%register_state_variable(self%id_dmspd,'dmspd','nmol/L','ice DMSPd',minimum=0.0_rk)  !initial_value=dmspdi_0
! Register diagnostic variables
      call self%register_horizontal_diagnostic_variable(self%id_dmspp,'dmspp','nmol/L','ice DMSPp',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_flysspd,'flysspd','nM per day','exudation/cell lysis',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fexuspd,'fexuspd','nM per day','sloppy feeding',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fspdbac,'fspdbac','nM per day','DMSPd bacterial consumption',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fspddms,'fspddms','nM per day','extracellular lyase',source=source_do_horizontal)      
      call self%register_horizontal_diagnostic_variable(self%id_fspdaiw,'fspdaiw','nM per day','DMSPd exchange at ice-water interface',source=source_do_horizontal)      
      call self%register_horizontal_diagnostic_variable(self%id_fsppdms,'fsppdms','nM per day','intracellular lyase',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fbacdms,'fbacdms','nM per day','bacterial lyase (yield)',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdmsbac,'fdmsbac','nM per day','DMS bacterial consumption',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdmspho,'fdmspho','nM per day','photolysis',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdmsaiw,'fdmsaiw','nM per day','DMS exchange at ice-water interface',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fspdmel,'fspdmel','nM per day','DMSPd release due to bottom ice melting',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fspdpon,'fspdpon','nM per day','DMSPd release due to meltpond drainage',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdmsmel,'fdmsmel','nM per day','DMS release due to bottom ice melting',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fdmspon,'fdmspon','nM per day','DMS release due to meltpond drainage',source=source_do_horizontal)
! Register environmental variables
      call self%register_horizontal_dependency(self%id_ts,standard_variables%surface_ice_temperature)
      call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)
  !call self%register_horizontal_dependency(self%id_bvf,standard_variables%brine_volume_fraction)
      call self%register_horizontal_dependency(self%id_tb,standard_variables%sea_ice_temperature)
      call self%register_horizontal_dependency(self%id_par,standard_variables%lowest_ice_layer_PAR)
      call self%register_horizontal_dependency(self%id_ice_hi,standard_variables%sea_ice_thickness)
      call self%register_global_dependency(self%id_dt,standard_variables%timestep)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_horizontal_dependency(self%id_ia,'uvic_icealgae_ia','','')      
      call self%register_horizontal_dependency(self%id_chl,'uvic_icealgae_chl','','')      
      call self%register_horizontal_dependency(self%id_fgrow,'uvic_icealgae_fgrow','','')
      call self%register_horizontal_dependency(self%id_fgraze,'uvic_icealgae_fgraze','','')
      call self%register_horizontal_dependency(self%id_fmort,'uvic_icealgae_fmort','','')
      call self%register_horizontal_dependency(self%id_fmelt,'uvic_icealgae_fmelt','','')
      call self%register_horizontal_dependency(self%id_fpond,'uvic_icealgae_fpond','','')
      call self%register_horizontal_dependency(self%id_hnu,'uvic_icealgae_hnu','','')
      call self%register_horizontal_dependency(self%id_lno3,'uvic_icealgae_lno3','','')
      call self%register_horizontal_dependency(self%id_lsil,'uvic_icealgae_lsil','','')
      call self%register_horizontal_dependency(self%id_lice,'uvic_icealgae_lice','','')
      call self%register_horizontal_dependency(self%id_llig,'uvic_icealgae_llig','','')
      call self%register_dependency(self%id_dmsSW,'uvic_dms_dms','','')
      call self%register_dependency(self%id_dmspdSW,'uvic_dms_dmspd','','')
      call self%request_coupling(self%id_ia,'uvic_icealgae_ia')      
      call self%request_coupling(self%id_chl,'uvic_icealgae_chl')      
      call self%request_coupling(self%id_fgrow,'uvic_icealgae_fgrow')
      call self%request_coupling(self%id_fgraze,'uvic_icealgae_fgraze')
      call self%request_coupling(self%id_fmort,'uvic_icealgae_fmort')
      call self%request_coupling(self%id_fmelt,'uvic_icealgae_fmelt')
      call self%request_coupling(self%id_fpond,'uvic_icealgae_fpond')
      call self%request_coupling(self%id_hnu,'uvic_icealgae_hnu')
      call self%request_coupling(self%id_lno3,'uvic_icealgae_lno3')
      call self%request_coupling(self%id_lsil,'uvic_icealgae_lsil')
      call self%request_coupling(self%id_lice,'uvic_icealgae_lice')
      call self%request_coupling(self%id_llig,'uvic_icealgae_llig')
      call self%request_coupling(self%id_dmsSW,'uvic_dms_dms')
      call self%request_coupling(self%id_dmspdSW,'uvic_dms_dmspd')

   end subroutine initialize

! Calling and saving the surface variables.
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_icedms),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk) :: dmspp,dmspd,dmspdSW,flysspd,fexuspd,fspdbac,fspddms,fspdaiw,fspdpon,fspdmel
   real(rk) :: dms,dmsSW,fsppdms,fbacdms,fdmsbac,fdmspho,fdmsaiw,fdmspon,fdmsmel
   real(rk) :: botgrowth,bvf,ts,tb,par,ia,chl,fgrow,fmort,fgraze,ice_hi,dt,fpond,fmelt,temp,hnu,md_dms,lno3,lsil,lice,llig
   _HORIZONTAL_LOOP_BEGIN_
   _GET_HORIZONTAL_(self%id_dms,dms)
   _GET_HORIZONTAL_(self%id_dmspd,dmspd)
   _GET_HORIZONTAL_(self%id_par,par)
   _GET_HORIZONTAL_(self%id_ice_hi,ice_hi)
   _GET_HORIZONTAL_(self%id_ia,ia)
   _GET_HORIZONTAL_(self%id_chl,chl)
   _GET_HORIZONTAL_(self%id_fmelt,fmelt)
   _GET_HORIZONTAL_(self%id_fpond,fpond)
   _GET_HORIZONTAL_(self%id_hnu,hnu)
   _GET_HORIZONTAL_(self%id_lno3,lno3)
   _GET_HORIZONTAL_(self%id_lsil,lsil)
   _GET_HORIZONTAL_(self%id_lice,lice)
   _GET_HORIZONTAL_(self%id_llig,llig)
   _GET_HORIZONTAL_(self%id_fmort,fmort)
   _GET_HORIZONTAL_(self%id_fgraze,fgraze)
   _GET_HORIZONTAL_(self%id_fgrow,fgrow)
   _GET_HORIZONTAL_(self%id_ts,ts)
  ! _GET_HORIZONTAL_(self%id_bvf,bvf)
   _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
   _GET_HORIZONTAL_(self%id_tb,tb)
   fmort  = fmort/self%spd
   fmelt  = fmelt/self%spd
   fpond  = fpond/self%spd
   fgraze = fgraze/self%spd
   fgrow = fgrow/self%spd
   _GET_(self%id_dmsSW,dmsSW)
   _GET_(self%id_dmspdSW,dmspdSW)
   _GET_(self%id_temp,temp)
   _GET_GLOBAL_(self%id_dt,dt)
   dmspp   = self%qi*chl
   md_dms  = 0.0_rk !0.02/10000*exp(-18.1/(0.008314*(temp+273.15)))
!  flysspd = self%k_lyi*(1/(min(lno3,lsil)+0.1))*(1/(lice+0.1))*dmspp
!  fexuspd = (1/(lice+0.1))*(self%f_exi+(1-min(lno3,lsil)*(1-self%f_exi)))*fgrow/ia*dmspp
   flysspd = self%k_lyi*(1/(min(lno3,lsil)+0.1))*dmspp
   fexuspd = (self%f_exi+(1-min(lno3,lsil)*(1-self%f_exi)))*fgrow/ia*dmspp  !jpnote: not used for now -- test ? 
!  fspdbac = self%k_dmspdi*(dmspd/(dmspd+self%h_dmspdi))*dmspd
!  fspdbac = self%k_dmspdi*10./(10.+par)*(dmspd/(dmspd+self%h_dmspdi))*dmspd
!  fspdbac = self%k_dmspdi*10./(10.+par)*dmspd
   fspdbac = self%k_dmspdi*dmspd!*(0.55-0.45*(tanh((chl-650.)/150.)))
   fspddms = self%k_exi*dmspd
   fspdaiw = md_dms/hnu*(dmspdSW-dmspd)/self%zia
   fsppdms = self%k_ini*dmspp
   fbacdms = self%yieldi*fspdbac
   fdmsbac = self%k_dmsi*dms!*(0.55-0.45*(tanh((chl-650.)/150.)))
!  fdmsbac = self%k_dmsi*10./(10.+par)*(dms/(dms+self%h_dmsi))*dms
!  fdmsbac = self%k_dmsi*10./(10.+par)*dms
   fdmspho = self%k_phoi*(par/(1+par))*dms
   fdmsaiw = md_dms/hnu*(dmsSW-dms)/self%zia
   fspdmel = 913./1000.*fmelt/ia*dmspd  !jpnote in em: fspdmel = fmelt/ia*dmspd
   fspdpon = fpond/ia*dmspd
   fdmsmel = 913./1000.*fmelt/ia*dms !jpnote in em: fdmsmel = fmelt/ia*dms
   fdmspon = fpond/ia*dms
!   if (bvf .gt. 5) then ! brine drainage
!    fspdpon = fspdpon + dmspd/20./86400.
!    fdmspon = fdmspon + dms/20./86400.
!   endif
   if (ice_hi .lt. 0.1_rk) then ! Ice algal layer is absent.
    fspdaiw = 0.0_rk
    fdmsaiw = 0.0_rk
    fspdmel = 0.0_rk
    fspdpon = 0.0_rk
    fdmsmel = 0.0_rk
    fdmspon = 0.0_rk
    _SET_SURFACE_ODE_(self%id_dmspd,(dmspdSW-dmspd)/dt)
    _SET_SURFACE_ODE_(self%id_dms,(dmsSW-dms)/dt)
   else
     _SET_SURFACE_ODE_(self%id_dmspd,flysspd+fexuspd-fspdbac-fspddms+fspdaiw+fspdmel-fspdpon)
     _SET_SURFACE_ODE_(self%id_dms,fsppdms+fbacdms+fspddms-fdmsbac-fdmspho+fdmsaiw+fdmsmel-fdmspon)
   endif
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dmspp,dmspp)   
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_flysspd,flysspd*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fexuspd,fexuspd*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fspdbac,fspdbac*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fspddms,fspddms*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fspdaiw,fspdaiw*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fsppdms,fsppdms*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fbacdms,fbacdms*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmsbac,fdmsbac*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmspho,fdmspho*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmsaiw,fdmsaiw*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fspdmel,fspdmel*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fspdpon,fspdpon*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmsmel,fdmsmel*self%spd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fdmspon,fdmspon*self%spd)
   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
end module uvic_icedms