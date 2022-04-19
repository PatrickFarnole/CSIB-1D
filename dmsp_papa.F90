#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: uvic_dmps_papa --- DMS(P) model based on Steiner & Denman (2008)
! taken from GOTM written by N. Steiner and adapted for FABM by H. Hayashida
!
! !INTERFACE:
   module uvic_dmsp_papa
!
! !DESCRIPTION:
! Modifications from Steiner & Denman 2008:
! 1) no parametrization with ultraviolet radiation (we can't, since don't have atmospheric model).
! 3) DMS in air is set to zero (hence, the air-sea DMS flux is always upwards).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public             :: type_uvic_dmsp_papa
!     Variable identifiers
      type (type_state_variable_id)                 :: id_dms,id_dmspd
      type (type_dependency_id)                     :: id_par,id_temp,id_ph1_zo1,id_ph1_nh4,id_ph2_zo2,id_ph2_nh4,id_ph1,id_ph2
      type (type_horizontal_dependency_id)          :: id_patm,id_wind
      type (type_diagnostic_variable_id)            :: id_dmspp1,id_dmspp2,id_zo1_dms,id_ph1_dms,id_zo2_dms,id_ph2_dms,id_yie_dms,id_dms_bac,id_dms_pho,id_zo1_dmspd,id_ph1_dmspd,id_zo2_dmspd,id_ph2_dmspd,id_dmspd_dms,id_dmspd_bac
      type (type_horizontal_diagnostic_variable_id) :: id_dms_air
!     Model parameters
      real(rk) :: q1,q2,d1,d2,d3,d4,d5,d6,d7,dg,pr,ad
      integer  :: mf,la
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the dmsp namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uvic_dmsp_papa), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk) :: q1,q2,d1,d2,d3,d4,d5,d6,d7,dg,pr,ad
   integer  :: mf,la
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%q1,'q1','mg-S per mg-N','S:N ratio for ph1',default=0.45_rk)
   call self%get_parameter(self%q2,'q2','mg-S per mg-N','S:N ratio for ph2',default=0.0_rk)
   call self%get_parameter(self%d1,'d1','-','proportion of dmspp released due to grazing as dmspd',default=0.8_rk)
   call self%get_parameter(self%d2,'d2','-','proportion of dmspp released due to lysis as dmspd',default=0.9_rk)
   call self%get_parameter(self%d3,'d3','-','proportion of dmspp released due to grazing as dms',default=0.08_rk)
   call self%get_parameter(self%d4,'d4','-','proportion of dmspp released due to lysis as dms',default=0.05_rk)
   call self%get_parameter(self%d5,'d5','d-1','free dmspd-lyase rate',default=0.04_rk,scale_factor=d_per_s)
   call self%get_parameter(self%d6,'d6','d-1','total microbial degradation rate on dmspd',default=5.5_rk,scale_factor=d_per_s)
   call self%get_parameter(self%d7,'d7','d-1','bacterial degradation rate on dms',default=0.35_rk,scale_factor=d_per_s)
   call self%get_parameter(self%dg,'dg','-','dms yield fraction',default=0.2_rk)
   call self%get_parameter(self%pr,'pr','d-1','photolysis rate',default=0.25_rk,scale_factor=d_per_s)
   call self%get_parameter(self%mf,'mf','-','method for air-sea DMS flux (0=Steiner&Denman 2008, 1=Nightingale 2000, 2=zero flux',default=0)

   ! Register state variables
   call self%register_state_variable(self%id_dms,'dms','nmol/L','Dimethylsulfide',minimum=0.0_rk)
   call self%register_state_variable(self%id_dmspd,'dmspd','nmol/L','Dissolved dimethylsulphoniopropionate',minimum=0.0_rk)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dmspp1,'dmspp1','nmol/L','Particulate DMSP in ph1')
   call self%register_diagnostic_variable(self%id_dmspp2,'dmspp2','nmol/L','particulate DMSP in ph2')
   call self%register_diagnostic_variable(self%id_zo1_dms,'zo1_dms','nmol/L/d','dms release via zo1 grazing')
   call self%register_diagnostic_variable(self%id_ph1_dms,'ph1_dms','nmol/L/d','dms release via ph1 lysis')
   call self%register_diagnostic_variable(self%id_zo2_dms,'zo2_dms','nmol/L/d','dms release via zo2 grazing')
   call self%register_diagnostic_variable(self%id_ph2_dms,'ph2_dms','nmol/L/d','dms release via ph2 lysis')
   call self%register_diagnostic_variable(self%id_yie_dms,'yie_dms','nmol/L/d','dms yield')
   call self%register_diagnostic_variable(self%id_dms_bac,'dms_bac','nmol/L/d','bacterial dms uptake')
   call self%register_diagnostic_variable(self%id_dms_pho,'dms_pho','nmol/L/d','photolysis')
   call self%register_diagnostic_variable(self%id_zo1_dmspd,'zo1_dmspd','nmol/L/d','dmspd release via zo1 grazing')
   call self%register_diagnostic_variable(self%id_ph1_dmspd,'ph1_dmspd','nmol/L/d','dmspd release via ph1 lysis')
   call self%register_diagnostic_variable(self%id_zo2_dmspd,'zo2_dmspd','nmol/L/d','dmspd release via zo2 grazing')
   call self%register_diagnostic_variable(self%id_ph2_dmspd,'ph2_dmspd','nmol/L/d','dmspd release via ph2 lysis')
   call self%register_diagnostic_variable(self%id_dmspd_dms,'dmspd_dms','nmol/L/d','dmspd lyase')
   call self%register_diagnostic_variable(self%id_dmspd_bac,'dmspd_bac','nmol/L/d','bacterial dmspd uptake')

   ! Register horizontal diagnostic variables
   call self%register_diagnostic_variable(self%id_dms_air,'dms_air','umol/m2/d','air-sea dms flux')

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_ph1_zo1,'uvic_npzd_papa_ph1_zo1','','')
   call self%register_dependency(self%id_ph1_nh4,'uvic_npzd_papa_ph1_nh4','','')
   call self%register_dependency(self%id_ph2_zo2,'uvic_npzd_papa_ph2_zo2','','')
   call self%register_dependency(self%id_ph2_nh4,'uvic_npzd_papa_ph2_nh4','','')
   call self%register_dependency(self%id_ph1,'uvic_npzd_papa_ph1','','')
   call self%register_dependency(self%id_ph2,'uvic_npzd_papa_ph2','','')
   call self%request_coupling(self%id_ph1_zo1,'uvic_npzd_papa_ph1_zo1')
   call self%request_coupling(self%id_ph1_nh4,'uvic_npzd_papa_ph1_nh4')
   call self%request_coupling(self%id_ph2_zo2,'uvic_npzd_papa_ph2_zo2')
   call self%request_coupling(self%id_ph2_nh4,'uvic_npzd_papa_ph2_nh4')
   call self%request_coupling(self%id_ph1,'uvic_npzd_papa_ph1')
   call self%request_coupling(self%id_ph2,'uvic_npzd_papa_ph2')

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
! ???
!
! !INPUT PARAMETERS:
   class (type_uvic_dmsp_papa),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: dms,dmspd
   real(rk)                   :: par,temp,ph1,ph2,ph1_zo1,ph1_nh4,ph2_zo2,ph2_nh4
   real(rk)                   :: dmspp1,dmspp2,zo1_dms,ph1_dms,zo2_dms,ph2_dms,yie_dms,dms_bac,dms_pho,zo1_dmspd,ph1_dmspd,zo2_dmspd,ph2_dmspd,dmspd_dms,dmspd_bac
   real(rk), parameter        :: wN = 14._rk,S2mol = 31.25_rk,spd = 86400._rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_dms,dms)
   _GET_(self%id_dmspd,dmspd)

   ! Retrieve environmental variables
   _GET_(self%id_temp,temp)
   _GET_(self%id_par,par)
   _GET_(self%id_ph1,ph1)
   _GET_(self%id_ph2,ph2)
   _GET_(self%id_ph1_zo1,ph1_zo1)
   _GET_(self%id_ph2_zo2,ph2_zo2)
   _GET_(self%id_ph1_nh4,ph1_nh4)
   _GET_(self%id_ph2_nh4,ph2_nh4)
   ph1_zo1 = ph1_zo1 / spd
   ph2_zo2 = ph2_zo2 / spd
   ph1_nh4 = ph1_nh4 / spd
   ph2_nh4 = ph2_nh4 / spd
   
   ! DMS and DMSPd budget components
   zo1_dms = wN*S2mol*ph1_zo1*self%d3*self%q1
   ph1_dms = wN*S2mol*ph1_nh4*self%d4*self%q1
   zo2_dms = wN*S2mol*ph2_zo2*self%d3*self%q2
   ph2_dms = wN*S2mol*ph2_nh4*self%d4*self%q2
   yie_dms = dmspd*self%dg*self%d6
   dms_bac = dms*self%d7
   dms_pho = dms*self%pr*par/4.0_rk
   zo1_dmspd = wN*S2mol*ph1_zo1*self%d1*self%q1
   ph1_dmspd = wN*S2mol*ph1_nh4*self%d2*self%q1
   zo2_dmspd = wN*S2mol*ph2_zo2*self%d1*self%q2
   ph2_dmspd = wN*S2mol*ph2_nh4*self%d2*self%q2
   dmspd_dms = dmspd*self%d5
   dmspd_bac = dmspd*self%d6
   dmspp1 = wN*S2mol*ph1*self%q1
   dmspp2 = wN*S2mol*ph2*self%q2

   ! Compute and save prognostic variables
   _SET_ODE_(self%id_dms,zo1_dms+ph1_dms+zo2_dms+ph2_dms+yie_dms+dmspd_dms-dms_bac-dms_pho)
   _SET_ODE_(self%id_dmspd,zo1_dmspd+ph1_dmspd+zo2_dmspd+ph2_dmspd-dmspd_dms-dmspd_bac)
   
   ! Save diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dmspp1,dmspp1)
   _SET_DIAGNOSTIC_(self%id_dmspp2,dmspp2)
   _SET_DIAGNOSTIC_(self%id_zo1_dms,zo1_dms*spd)
   _SET_DIAGNOSTIC_(self%id_ph1_dms,ph1_dms*spd)
   _SET_DIAGNOSTIC_(self%id_zo2_dms,zo2_dms*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_dms,ph2_dms*spd)
   _SET_DIAGNOSTIC_(self%id_yie_dms,yie_dms*spd)
   _SET_DIAGNOSTIC_(self%id_dms_bac,dms_bac*spd)
   _SET_DIAGNOSTIC_(self%id_dms_pho,dms_pho*spd)
   _SET_DIAGNOSTIC_(self%id_zo1_dmspd,zo1_dmspd*spd)
   _SET_DIAGNOSTIC_(self%id_ph1_dmspd,ph1_dmspd*spd)
   _SET_DIAGNOSTIC_(self%id_zo2_dmspd,zo2_dmspd*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_dmspd,ph2_dmspd*spd)
   _SET_DIAGNOSTIC_(self%id_dmspd_dms,dmspd_dms*spd)
   _SET_DIAGNOSTIC_(self%id_dmspd_bac,dmspd_bac*spd)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

! Calling and saving the surface variables.
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_dmsp_papa),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk) :: dmsa,dms,temp,wind,patm,pres_atm,adms(4),ddms(2),rk0,sc,osc,W,rkb,k_600,k_dms,henry,dms_air
   real(rk) :: spd = 86400.
   data adms/2674.0,147.12,3.726,0.038/
   data ddms/-10.1794,3761.33/

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_dms,dms)
   _GET_(self%id_temp,temp)
   _GET_HORIZONTAL_(self%id_wind,wind)
   _GET_HORIZONTAL_(self%id_patm,patm)
   pres_atm = patm*9.87e-6 ! converting from [Pa] to [atm]
   sc = adms(1)-adms(2)*temp+adms(3)*temp**2-adms(4)*temp**3
   select case (self%mf)
    case (0) ! Liss & Merlivat 1986
     if (wind.lt.3.6_rk) then
      rk0 = 0.17_rk*wind
     else if (wind.ge.3.6_rk) then
      rk0 = 2.85_rk*wind-9.65_rk
     end if
     W = 3.84e-6*wind**3.41_rk  ! Monahan & O'Muircheartaigh 1980
     osc = exp(ddms(1) + ddms(2)/(temp+273.15_rk))/pres_atm ! Ostwald solublity coefficient based on Saltzman 1993
     rkb = 2450._rk*W/(osc*(1._rk+(14._rk*osc*sqrt(sc))**(-1.0_rk/1.2_rk))**1.2_rk)
     k_600 = rk0 + rkb
    case (1) ! Nightingale 2000
     k_600 = 0.222_rk*wind**2+0.333_rk*wind
    case (2) ! No air-sea flux
     k_600 = 0.0_rk
   end select !end if
   k_dms = 0.24_rk*k_600/sqrt(sc/600._rk) ! convert from [cm/h] to [m/d]
   henry = exp(12.64_rk-3547_rk/(temp+273.15_rk))/(0.08205_rk*(temp+273.15_rk))
   dmsa = 0.0_rk
   dms_air = (dms-dmsa/henry)*k_dms/spd
   _SET_SURFACE_EXCHANGE_(self%id_dms,-dms_air)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dms_air,dms_air*spd)  

   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
!-----------------------------------------------------------------------

   end module uvic_dmsp_papa
