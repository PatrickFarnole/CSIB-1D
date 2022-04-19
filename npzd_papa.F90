#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: uvic_npzd_papa --- NPZD biogeochemical model based upon
! Monahan and Denman (2004),
! with minor modifications by Steiner et al. (2006)
! taken from GOTM written by N. Steiner and adapted for FABM by H. Hayashida
!
! !INTERFACE:
   module uvic_npzd_papa
!
! !DESCRIPTION:
! Further modifications in this version:
! 1) Nitrification parameterization is based on Kawamiya et al. 1992.
! 2) Nitrate and Silicate are not replenished with deep water.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uvic_npzd_papa
!     Variable identifiers
      type (type_state_variable_id)        :: id_ph1,id_ph2,id_zo1,id_no3,id_nh4,id_det,id_sil,id_psi
      type (type_dependency_id)            :: id_par,id_temp
      type (type_global_dependency_id)     :: id_day
      type (type_diagnostic_variable_id)   :: id_nut_ph1,id_ph1_zo1,id_ph1_nh4,id_nut_ph2,id_ph2_zo2,id_ph2_nh4,id_ph2_det,id_ph2_agg,id_fo1_zo1,id_zo1_zo2,id_zo1_nh4,id_no3_phy,id_nh4_no3,id_nh4_phy,id_det_nh4,id_zo2_nh4,id_fo1_det,id_det_zo1,id_ph2_psi,id_psi_sil,id_sil_ph2,id_lfe1,id_lfe2,id_llig,id_lnut,id_lsil
!     Model parameters
      real(rk) :: kc,alpha,num,kn,a,ks,lfe1,lfe2,mpa,mpd,beta,rm,rc,kp,kz,ga1,ga2,pd,mza,mca,wp2,wd,re,resi,q10p,q10z,q10b,rsin,nit0,knit
   contains
      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uvic_npzd_papa), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk) :: kc,alpha,num,kn,a,ks,lfe1,lfe2,mpa,mpd,beta,rm,rc,kp,kz,ga1,ga2,pd,mza,mca,wp2,wd,re,resi,q10p,q10z,q10b,rsin,nit0,knit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%kc,'kc','m-1 (mmolN m-3)-1','PAR attenuation coefficient for ph1+ph2+det')
   call self%get_parameter(self%alpha,'alpha','d-1 (W m-2)-1','Initial slope of P-I curve',scale_factor=d_per_s)
   call self%get_parameter(self%num,'num','d-1','maximum phytoplankton specific growth rate',scale_factor=d_per_s)
   call self%get_parameter(self%kn,'kn','mmolN m-3','nitrogen uptake half saturation constant')
   call self%get_parameter(self%a,'a','mmolN m-3','scale factor for nitrogen preference function')
   call self%get_parameter(self%ks,'ks','mmolSi m-3','silicic acid half saturation constant')
   call self%get_parameter(self%lfe1,'lfe1','-','iron limitation value for ph1')
   call self%get_parameter(self%lfe2,'lfe2','-','iron limitation value for ph2')
   call self%get_parameter(self%mpa,'mpa','d-1','phytoplankton excretion rate (to nh4)',scale_factor=d_per_s)
   call self%get_parameter(self%mpd,'mpd','d-1','ph1 mortality rate (to det)',scale_factor=d_per_s)
   call self%get_parameter(self%beta,'beta','d-1 (mmolN m-3)-1','ph2 aggregation rate',scale_factor=d_per_s)
   call self%get_parameter(self%rm,'rm','d-1','zo1 maximum specific grazing rate',scale_factor=d_per_s)
   call self%get_parameter(self%rc,'rc','d-1','zo2 maximum specific grazing rate',scale_factor=d_per_s)
   call self%get_parameter(self%kp,'kp','mmolN m-3','zo1 grazing half saturation constant')
   call self%get_parameter(self%kz,'kz','mmolN m-3','zo2 grazing half saturation constant')
   call self%get_parameter(self%ga1,'ga1','-','zo1 grazing assimilation fraction')
   call self%get_parameter(self%ga2,'ga2','-','zo2 grazing assimilation fraction')
   call self%get_parameter(self%pd,'pd','-','zo1 grazing preference for det')
   call self%get_parameter(self%mza,'mza','d-1','zo1 excretion rate',scale_factor=d_per_s)
   call self%get_parameter(self%mca,'mca','-','zo2 grazing fraction transferred to nh4')
   call self%get_parameter(wp2,'wp2','m d-1','sinking speed for ph2',scale_factor=d_per_s)
   call self%get_parameter(wd,'wd','m d-1','sinking speed for det',scale_factor=d_per_s)
   call self%get_parameter(self%re,'re','d-1','remineralization rate for nitrogen',scale_factor=d_per_s)
   call self%get_parameter(self%resi,'resi','d-1','remineralization rate for silicon',scale_factor=d_per_s)
   call self%get_parameter(self%q10p,'q10p','-','temperature sensivity for phytoplankton')
   call self%get_parameter(self%q10z,'q10z','-','temperature sensivity for zooplankton')
   call self%get_parameter(self%q10b,'q10b','-','temperature sensivity for bacteria')
   call self%get_parameter(self%rsin,'rsin','molSi:molN','ph2 uptake ratio for silicon to nitrogen')
   call self%get_parameter(self%nit0,'nit0','d-1','nitrification rate at 0 degree celcius',scale_factor=d_per_s)
   call self%get_parameter(self%knit,'knit','degC-1','temperature sensivity for nitrification')
   
   ! Register state variables
   call self%register_state_variable(self%id_ph1,'ph1','umol/L','Small phytoplankton',minimum=0.0_rk)
   call self%register_state_variable(self%id_ph2,'ph2','umol/L','Large phytoplankton (Diatoms)',minimum=0.0_rk,vertical_movement=wp2)
   call self%register_state_variable(self%id_zo1,'zo1','umol/L','Microzooplankton',minimum=0.0_rk)
   call self%register_state_variable(self%id_no3,'no3','umol/L','Nitrate',minimum=0.0_rk)                         
   call self%register_state_variable(self%id_nh4,'nh4','umol/L','Ammonium',minimum=0.0_rk)
   call self%register_state_variable(self%id_det,'det','umol/L','Detritus',minimum=0.0_rk,vertical_movement=wd)                         
   call self%register_state_variable(self%id_sil,'sil','umol/L','Silicate',minimum=0.0_rk)                         
   call self%register_state_variable(self%id_psi,'psi','umol/L','Particulate silica',minimum=0.0_rk,vertical_movement=wd)                         

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_nut_ph1,'nut_ph1','umol/L/d','primary production by ph1')
   call self%register_diagnostic_variable(self%id_ph1_zo1,'ph1_zo1','umol/L/d','zo1 grazing on ph1')
   call self%register_diagnostic_variable(self%id_ph1_nh4,'ph1_nh4','umol/L/d','excretion by ph1')
   call self%register_diagnostic_variable(self%id_nut_ph2,'nut_ph2','umol/L/d','primary production by ph2')
   call self%register_diagnostic_variable(self%id_ph2_zo2,'ph2_zo2','umol/L/d','zo2 grazing on ph2')
   call self%register_diagnostic_variable(self%id_ph2_nh4,'ph2_nh4','umol/L/d','excretion by ph2 to nh4')
   call self%register_diagnostic_variable(self%id_ph2_det,'ph2_det','umol/L/d','excretion by ph2 to det')
   call self%register_diagnostic_variable(self%id_ph2_agg,'ph2_agg','umol/L/d','aggregation loss of ph2')
   call self%register_diagnostic_variable(self%id_fo1_zo1,'fo1_zo1','umol/L/d','assimilated grazing by zo1')
   call self%register_diagnostic_variable(self%id_zo1_zo2,'zo1_zo2','umol/L/d','zo2 grazing on zo1')
   call self%register_diagnostic_variable(self%id_zo1_nh4,'zo1_nh4','umol/L/d','excretion by zo1')
   call self%register_diagnostic_variable(self%id_no3_phy,'no3_phy','umol/L/d','no3 uptake by ph1 & ph2')
   call self%register_diagnostic_variable(self%id_nh4_no3,'nh4_no3','umol/L/d','nitrification')
   call self%register_diagnostic_variable(self%id_nh4_phy,'nh4_phy','umol/L/d','nh4 uptake by ph1 & ph2')
   call self%register_diagnostic_variable(self%id_det_nh4,'det_nh4','umol/L/d','remineralization to nh4')
   call self%register_diagnostic_variable(self%id_zo2_nh4,'zo2_nh4','umol/L/d','excretion by zo2')
   call self%register_diagnostic_variable(self%id_fo1_det,'fo1_det','umol/L/d','sloppy feeding/egestion by zo1')
   call self%register_diagnostic_variable(self%id_det_zo1,'det_zo1','umol/L/d','zo1 grazing on det')
   call self%register_diagnostic_variable(self%id_ph2_psi,'ph2_psi','umol/L/d','silicon content of ph2')
   call self%register_diagnostic_variable(self%id_psi_sil,'psi_sil','umol/L/d','remineralization to sil')
   call self%register_diagnostic_variable(self%id_sil_ph2,'sil_ph2','umol/L/d','sil uptake by ph2')
   call self%register_diagnostic_variable(self%id_lfe1,'lfe1','-','iron limitation on ph1')
   call self%register_diagnostic_variable(self%id_lfe2,'lfe2','-','iron limitaiton on ph2')
   call self%register_diagnostic_variable(self%id_llig,'llig','-','light limitation')
   call self%register_diagnostic_variable(self%id_lnut,'lnut','-','nitrogen limitation')
   call self%register_diagnostic_variable(self%id_lsil,'lsil','-','silicate limitation')

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_global_dependency(self%id_day,standard_variables%number_of_days_since_start_of_the_year)

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_uvic_npzd_papa),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: ph1,ph2,zo1,zo2,no3,nh4,det,sil,psi
   real(rk)                   :: par,temp,day
   real(rk)                   :: texp,qp,qz,qb,lnut,llig,lsil,tnu,trm,trc,tmpa,tmpd,tmza,tre,trsi
   real(rk)                   :: nut_ph1,ph1_zo1,ph1_nh4,nut_ph2,ph2_zo2,ph2_nh4,ph2_det,ph2_agg,fo1_zo1,zo1_zo2,zo1_nh4,no3_phy,nh4_no3,nh4_phy,det_nh4,zo2_nh4,fo1_det,det_zo1,ph2_psi,psi_sil,sil_ph2
   real(rk), parameter        :: spd = 86400.
   real(rk), dimension(365)   :: zo2_365
   !Daily mesozooplankton spline-interpolated from Fulton 1956-80 series
   !7.4e-4 converts annual series to nitrogen biomass units
   zo2_365 = 7.4e-4*(/ 14.4667,14.3579,14.2537,14.1541,14.0592,13.9689,13.8834,13.8027,13.7269,13.6561,13.5904,13.5298,13.4745,13.4245,13.3800,13.3411,13.3082,13.2813,13.2609,13.2473,13.2406,13.2413,13.2497,13.2661,13.2908,13.3243,13.3668,13.4188,13.4807,13.5528,13.6357,13.7297,13.8353,13.9530,14.0832,14.2266,14.3836,14.5547,14.7406,14.9418,15.1589,15.3926,15.6436,15.9125,16.2000,16.5067,16.8327,17.1778,17.5419,17.9250,18.3271,18.7481,19.1880,19.6470,20.1251,20.6224,21.1390,21.6751,22.2309,22.8066,23.4024,24.0187,24.6557,25.3138,25.9934,26.6947,27.4184,28.1647,28.9343,29.7276,30.5452,31.3877,32.2558,33.1500,34.0712,35.0200,35.9977,37.0072,38.0519,39.1354,40.2612,41.4328,42.6540,43.9285,45.2601,46.6526,48.1102,49.6367,51.2365,52.9139,54.6731,56.5187,58.4555,60.4880,62.6213,64.8603,67.2104,69.6768,72.2651,74.9810,77.8304,80.8194,83.9542,87.2414,90.6877,94.3000,98.0818,102.0221,106.1063,110.3200,114.6491,119.0798,123.5980,128.1901,132.8426,137.5418,142.2745,147.0270,151.7862,156.5386,161.2708,165.9696,170.6213,175.2126,179.7298,184.1591,188.4869,192.6991,196.7814,200.7197,204.4993,208.1055,211.5234,214.7375,217.7323,220.4919,223.0000,225.2449,227.2336,228.9776,230.4881,231.7760,232.8515,233.7251,234.4066,234.9055,235.2314,235.3932,235.3997,235.2599,234.9820,234.5741,234.0445,233.4010,232.6514,231.8029,230.8633,229.8398,228.7395,227.5694,226.3366,225.0478,223.7099,222.3295,220.9133,219.4680,218.0000,216.5148,215.0138,213.4970,211.9646,210.4168,208.8537,207.2754,205.6820,204.0735,202.4500,200.8114,199.1577,197.4890,195.8051,194.1061,192.3917,190.6618,188.9164,187.1552,185.3781,183.5848,181.7752,179.9489,178.1057,176.2452,174.3671,172.4712,170.5569,168.6239,166.6718,164.7000,162.7085,160.6986,158.6721,156.6305,154.5757,152.5093,150.4328,148.3479,146.2562,144.1594,142.0589,139.9564,137.8535,135.7516,133.6526,131.5578,129.4689,127.3875,125.3152,123.2537,121.2045,119.1693,117.1499,115.1479,113.1650,111.2030,109.2636,107.3488,105.4603,103.6000,101.7694,99.9683,98.1961,96.4520,94.7356,93.0464,91.3838,89.7474,88.1368,86.5516,84.9915,83.4562,81.9454,80.4590,78.9966,77.5583,76.1437,74.7529,73.3858,72.0423,70.7226,69.4265,68.1543,66.9061,65.6819,64.4820,63.3065,62.1559,61.0303,59.9300,58.8553,57.8055,56.7798,55.7774,54.7977,53.8400,52.9035,51.9877,51.0919,50.2155,49.3582,48.5192,47.6982,46.8946,46.1080,45.3381,44.5844,43.8465,43.1241,42.4170,41.7247,41.0471,40.3838,39.7346,39.0993,38.4777,37.8697,37.2750,36.6936,36.1253,35.5700,35.0276,34.4977,33.9800,33.4739,32.9793,32.4957,32.0227,31.5601,31.1076,30.6649,30.2316,29.8076,29.3926,28.9864,28.5886,28.1992,27.8179,27.4446,27.0790,26.7209,26.3704,26.0271,25.6910,25.3619,25.0398,24.7245,24.4160,24.1141,23.8188,23.5300,23.2477,22.9715,22.7014,22.4370,22.1782,21.9246,21.6762,21.4327,21.1939,20.9596,20.7297,20.5040,20.2823,20.0645,19.8503,19.6398,19.4326,19.2287,19.0279,18.8302,18.6353,18.4432,18.2537,18.0668,17.8822,17.7000,17.5200,17.3421,17.1662,16.9922,16.8200,16.6496,16.4812,16.3150,16.1512,15.9900,15.8318,15.6767,15.5250,15.3769,15.2327,15.0928,14.9572,14.8264,14.7005,14.5800 /)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_ph1,ph1)
   _GET_(self%id_ph2,ph2)
   _GET_(self%id_zo1,zo1)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_det,det)
   _GET_(self%id_sil,sil)
   _GET_(self%id_psi,psi)
   ! Retrieve environmental variables
   _GET_(self%id_temp,temp)
   _GET_(self%id_par,par)
   _GET_GLOBAL_(self%id_day,day)

   zo2 = zo2_365(nint(day))

   texp = (temp-10._rk)/10._rk
   qp = self%q10p**texp
   qz = self%q10z**texp
   qb = self%q10b**texp
   tnu = self%num*qp
   trm = self%rm*qz*(ph1+self%pd*det)**2/(self%kp**2+(ph1+self%pd*det)**2)
   trc = self%rc*qz*(ph2+zo1)**2/(self%kz**2+(ph2+zo1)**2)
   tmpa = self%mpa*qp
   tmpd = self%mpd*qp
   tmza = self%mza*qz
   tre = self%re*qb
   trsi = self%resi*qb
   llig = 1.0_rk-exp(-self%alpha*par/tnu)
   lnut = (no3+nh4)/(self%kn+no3+nh4)
   lsil = sil/(self%ks+sil)

   nut_ph1 = tnu*ph1*min(lnut,llig,self%lfe1)
   ph1_zo1 = trm*ph1/(ph1+self%pd*det)*zo1
   ph1_nh4 = tmpa*ph1
   nut_ph2 = tnu*ph2*min(lnut,llig,self%lfe2,lsil)
   ph2_zo2 = trc*ph2/(ph2+zo1)*zo2
   ph2_nh4 = tmpa*ph2
   ph2_det = tmpd*ph2
   ph2_agg = self%beta*ph2**2
   fo1_zo1 = self%ga1*trm*zo1
   zo1_zo2 = trc*zo1/(ph2+zo1)*zo2
   zo1_nh4 = tmza*zo1
   no3_phy = tnu*(ph1+ph2)*self%a/(self%a+nh4)*no3/(no3+nh4)
   nh4_no3 = self%nit0*exp(self%knit*temp)*nh4
   nh4_phy = tnu*(ph1+ph2)-no3_phy
   det_nh4 = tre*det
   zo2_nh4 = self%mca*trc*zo2
   ph2_det = tmpd*ph2
   fo1_det = (1-self%ga1)*trm*zo1
   det_zo1 = trm*self%pd*det/(ph1+self%pd*det)*zo1
   ph2_psi = self%rsin*tmpd*ph2
   psi_sil = trsi*psi
   sil_ph2 = self%rsin*nut_ph2
   if (nut_ph1+nut_ph2.gt.0.0) then ! pp reduced if no3_phy or nh4_phy is zero
    if ((nh4_phy.eq.0.0).or.(no3_phy.eq.0.0)) then
     nut_ph1=nut_ph1*(no3_phy+nh4_phy)/(nut_ph1+nut_ph2)
     nut_ph2=nut_ph2*(no3_phy+nh4_phy)/(nut_ph1+nut_ph2)
    endif
   endif

   ! Compute and save prognostic variables

   _SET_ODE_(self%id_ph1,nut_ph1-ph1_zo1-ph1_nh4)
   _SET_ODE_(self%id_ph2,nut_ph2-ph2_zo2-ph2_nh4-ph2_det-ph2_agg)
   _SET_ODE_(self%id_zo1,fo1_zo1-zo1_zo2-zo1_nh4)
   _SET_ODE_(self%id_no3,-no3_phy+nh4_no3)
   _SET_ODE_(self%id_nh4,-nh4_phy+det_nh4+ph1_nh4+ph2_nh4+zo1_nh4+zo2_nh4-nh4_no3)
   _SET_ODE_(self%id_det,ph2_det+fo1_det-det_zo1-det_nh4)
   _SET_ODE_(self%id_psi,ph2_psi-psi_sil)
   _SET_ODE_(self%id_sil,-sil_ph2+psi_sil)
   
   ! Save diagnostic variables
   _SET_DIAGNOSTIC_(self%id_nut_ph1,nut_ph1*spd)
   _SET_DIAGNOSTIC_(self%id_ph1_zo1,ph1_zo1*spd)
   _SET_DIAGNOSTIC_(self%id_ph1_nh4,ph1_nh4*spd)
   _SET_DIAGNOSTIC_(self%id_nut_ph2,nut_ph2*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_zo2,ph2_zo2*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_nh4,ph2_nh4*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_det,ph2_det*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_agg,ph2_agg*spd)
   _SET_DIAGNOSTIC_(self%id_fo1_zo1,fo1_zo1*spd)
   _SET_DIAGNOSTIC_(self%id_zo1_zo2,zo1_zo2*spd)
   _SET_DIAGNOSTIC_(self%id_zo1_nh4,zo1_nh4*spd)
   _SET_DIAGNOSTIC_(self%id_no3_phy,no3_phy*spd)
   _SET_DIAGNOSTIC_(self%id_nh4_no3,nh4_no3*spd)
   _SET_DIAGNOSTIC_(self%id_nh4_phy,nh4_phy*spd)
   _SET_DIAGNOSTIC_(self%id_det_nh4,det_nh4*spd)
   _SET_DIAGNOSTIC_(self%id_zo2_nh4,zo2_nh4*spd)
   _SET_DIAGNOSTIC_(self%id_fo1_det,fo1_det*spd)
   _SET_DIAGNOSTIC_(self%id_det_zo1,det_zo1*spd)
   _SET_DIAGNOSTIC_(self%id_ph2_psi,ph2_psi*spd)
   _SET_DIAGNOSTIC_(self%id_psi_sil,psi_sil*spd)
   _SET_DIAGNOSTIC_(self%id_sil_ph2,sil_ph2*spd)
   _SET_DIAGNOSTIC_(self%id_lfe1,self%lfe1)
   _SET_DIAGNOSTIC_(self%id_lfe2,self%lfe2)
   _SET_DIAGNOSTIC_(self%id_llig,llig)
   _SET_DIAGNOSTIC_(self%id_lnut,lnut)
   _SET_DIAGNOSTIC_(self%id_lsil,lsil)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_uvic_npzd_papa), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: ph1,ph2,det
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_ph1,ph1) ! small phytoplankton
   _GET_(self%id_ph2,ph2) ! large phytoplankton
   _GET_(self%id_det,det) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(ph1+ph2+det))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------

   end module uvic_npzd_papa

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
