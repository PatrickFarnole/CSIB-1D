!-----------------------------------------------------
!
! University of Victoria Ecosystem model
!
! Contributors: Hakase Hayashida, Nadja Steiner, Eric Mortenson, Adam Mohanan
!
! Model description: The model will be described in a paper.
!
!------------------------------------------------------

#include "fabm_driver.h"

module uvic_eco

   use fabm_types
   use fabm_expressions

   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_eco
! Declare prognostic variables      
      type (type_state_variable_id) :: id_ph1,id_ph2,id_zo1,id_zo2,id_nh4,id_no3,id_de1,id_de2,id_bsi,id_sil
! Declare diagnostic variables
      type (type_diagnostic_variable_id) :: id_chl,id_chl1,id_chl2,id_fnutph1,id_fnutph2,id_fph1nh4,id_fph2nh4,id_fph2sil,id_fph1zo1,id_fph2zo2,id_fph2de2,id_fzo1nh4,id_fzo1zo2,id_fzo2htl,id_fnh4phy,id_fzo2nh4,id_fde1nh4,id_fde2nh4,id_fnh4no3,id_fno3phy,id_foo1de1,id_foo2de2,id_fde1zo1,id_fde2zo2,id_foo2bsi,id_fph2bsi,id_fbsisil,id_fsilph2,id_foo1zo1,id_foo2zo2,id_lf1,id_lf2,id_nf,id_sf,id_f1f,id_f2f,id_par_diag
! Declare horizontal diagnostic variables
      type (type_horizontal_diagnostic_variable_id) :: id_stemp,id_fialde2,id_fialbsi,id_fialph2
! Declare environmental variables
      type (type_dependency_id) :: id_temp,id_par,id_density
! Declare horizontal environmental variables
      type (type_horizontal_dependency_id) :: id_botmelt,id_botgrowth,id_ia,id_iceno3,id_icenh4,id_icesil,id_fmelt,id_fpond,id_fpondno3,id_fpondnh4,id_fpondsil,id_fmort,id_fmort2,id_fskelno3,id_fskelnh4,id_fskelsil
! Declare namelist parameters
      real(rk) :: ac,f_seed,ph1_0,ph2_0,zo1_0,zo2_0,no3_0,nh4_0,de1_0,de2_0,bsi_0,sil_0,w1,w2,mu1,mu2,kn,rpp1,rpp2,mp1,mp2,gz1,kz1,az1,az2,mz1,rc,pp1,pp2,pd1,pd2,pz1,gz2,kz2,mz2,rd1,rd2,rd3,rpf,rn0,knt,qp,qz,qb,agg,rsin,ks,pmin
! Declare anything else used in all procedures
      real(rk) :: spd = 86400.0_rk ! Seconds Per Day (spd)
      real(rk) :: zia,ac_ia
      !character(64),dimension(12) :: models
      logical :: use_icealgae = .true.

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      
   end type

contains

   subroutine initialize(self,configunit)
   class (type_uvic_eco), intent(inout), target :: self
   integer, intent(in)                          :: configunit
! Declare namelist parameters
   real(rk) :: ac,f_seed,ph1_0,ph2_0,zo1_0,zo2_0,no3_0,nh4_0,de1_0,de2_0,bsi_0,sil_0,w1,w2,mu1,mu2,kn,rpp1,rpp2,mp1,mp2,gz1,kz1,az1,az2,mz1,rc,pp1,pp2,pd1,pd2,pz1,gz2,kz2,mz2,rd1,rd2,rd3,rpf,rn0,knt,qp,qz,qb,agg,rsin,ks,pmin
   real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,rpi,t_sens,nu,md_no3,md_sil,chl2n,sil2n
 !  character(64),dimension(12) :: models
! Define the namelist
   !namelist /fabm_nml/ models
  ! namelist /uvic_eco/ ac,f_seed,ph1_0,ph2_0,zo1_0,zo2_0,no3_0,nh4_0,de1_0,de2_0,bsi_0,sil_0,w1,w2,mu1,mu2,kn,rpp1,rpp2,mp1,mp2,gz1,kz1,az1,az2,mz1,rc,pp1,pp2,pd1,pd2,pz1,gz2,kz2,mz2,rd1,rd2,rd3,rpf,rn0,knt,qp,qz,qb,agg,rsin,ks,pmin
  ! namelist /uvic_icealgae/ r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,rpi,t_sens,nu,md_no3,md_sil,chl2n,sil2n
!  Initialize parameters to default values.
   ac     = 0.03_rk
   f_seed = 0.0_rk
  ! ph1_0 = 1.0_rk 
  ! ph2_0 = 0.5_rk
  ! zo1_0 = 0.2_rk
  ! zo2_0 = 0.1_rk
  ! no3_0 = 10.0_rk
  ! nh4_0 = 10.0_rk
  ! de1_0 = 1.0_rk
  ! de2_0 = 1.0_rk
  ! bsi_0 = 1.0_rk
  ! sil_0 = 5.0_rk
   w1    = 6.0_rk
   w2    = 6.0_rk
   mu1   = 2.0_rk
   mu2   = 2.0_rk
   kn    = 0.1_rk
   rpp1  = 0.05_rk
   rpp2  = 0.05_rk
   mp1   = 0.05_rk
   mp2   = 0.05_rk
   gz1   = 1.3_rk
   kz1   = 0.6_rk
   az1   = 0.7_rk
   az2   = 0.1_rk
   mz1   = 0.1_rk
   rc    = 0.003_rk
   pp1   = 1.0_rk
   pp2   = 1.0_rk
   pd1   = 0.5_rk
   pd2   = 0.5_rk
   pz1   = 1.0_rk
   gz2   = 0.8_rk  
   kz2   = 0.75_rk 
   mz2   = 0.3_rk    
   rd1   = 0.1_rk
   rd2   = 0.1_rk
   rd3   = 0.1_rk
   rpf   = 0.2_rk     
   qp    = 2.0_rk        
   qz    = 3.0_rk   
   qb    = 3.0_rk
   agg   = 0.07_rk
   rn0   = 0.03_rk
   knt   = 0.0693_rk
   rsin  = 2.0_rk
   ks    = 2.0_rk
   pmin  = 0.01_rk
   zia = 0.03_rk
   ac_ia = 0.007_rk


   ac     = 0.03
   f_seed  = 0.1
  ! ph1_0   = 0.01
  ! ph2_0   = 0.01
  ! zo1_0   = 0.01
  ! zo2_0   = 0.01
  ! nh4_0   = 0.01
  ! no3_0   = 7.2
  ! de1_0   = 0.01
  ! de2_0   = 0.01
  ! bsi_0   = 0.01
  ! sil_0   = 14.7
   w1      = -1 !L09
   w2      = -50 !L09 6,Steiner 2006
   mu1     = 0.5 !0.341 Paranjape 1987
   mu2     = 2.0 !steiner06,1.5,dupont (eppley at 10deg.~1.6) !0.567, !Fig.3: Suzuki and Takahashi 1995
   kn      = 1 !L09
   rpp1    = 0.3 !0.07125 Cota & Smith 1991
   rpp2    = 0.3 !0.045 Lavoie 2009
   mp1     = 0.03 ! L09
   mp2     = 0.03 ! L09
   mz1     = 0.03 !0.377, !Forest 2011
   mz2     = 0.03 !L09
   gz1     = 1  !1.3,!steiner06, !paranjape 1987
   gz2     = 1  !0.8,!steiner06, !0.5 L09
   kz1     = 1  !0.6,!1, !same as below
   kz2     = 1  !0.75,!1, !L09
   az1     = 0.7 !L09
   az2     = 0.7 !same as above
   rc      = 0.1
   pp1     = 1
   pp2     = 1
   pd1     = 0.5
   pd2     = 0.5 !L09
   pz1     = 0.5
   rd1     = 0.03
   rd2     = 0.3 !L09
   rd3     = 0.01
   rpf     = 0.2 !steiner06
   rn0     = 0.03 !0.01 steiner06, 0.03 Kawamiya
   knt     = 0.0693 !Kawamiya 
   qp      = 1.55  ! Table 2 Suzuki and Takahashi 1995  
   qz      = 2.14 !PISCES  
   qb      = 1.9 !PISCES
   agg     = 0.05  !steiner06
   rsin    = 1.7 ! L09
   ks      = 4 !L09
   pmin    = 0.01

#if 0 
   !  Read namelist parameters
   read(configunit,uvic_eco)
   close(configunit)
   open(configunit,file='fabm.nml')
   read(configunit,fabm_nml)
   self%models = models
   if(any(self%models.eq.'uvic_icealgae'))then
    zia = 0.03_rk
    ac_ia = 0.007_rk
    read(configunit,uvic_icealgae)
    self%zia = zia
    self%ac_ia = ac_ia
   endif
   close(configunit)
   open(configunit,file='fabm.nml')
#endif 

   !  Register namelist parameters which will be used in other routines
   self%ac     = ac
   self%f_seed = f_seed
   self%w1  = w1 / self%spd
   self%w2  = w2 / self%spd
   self%mu1 = mu1 / self%spd
   self%mu2 = mu2 / self%spd
   self%kn  = kn
   self%rpp1 = rpp1
   self%rpp2 = rpp2
   self%mp1 = mp1 / self%spd
   self%mp2 = mp2 / self%spd
   self%gz1 = gz1 / self%spd
   self%kz1 = kz1
   self%az1 = az1
   self%az2 = az2
   self%mz1 = mz1 / self%spd
   self%rc  = rc / self%spd
   self%pp1 = pp1
   self%pp2 = pp2
   self%pd1 = pd1
   self%pd2 = pd2
   self%pz1 = pz1
   self%gz2 = gz2 / self%spd
   self%kz2 = kz2
   self%mz2 = mz2 / self%spd
   self%rd1 = rd1 / self%spd
   self%rd2 = rd2 / self%spd
   self%rd3 = rd3 / self%spd
   self%rpf = rpf
   self%rn0 = rn0 / self%spd
   self%knt = knt
   self%qp  = qp
   self%qz  = qz
   self%qb  = qb
   self%agg = agg / self%spd
   self%rsin= rsin
   self%ks  = ks
   self%pmin= pmin

   self%zia = zia
   self%ac_ia = ac_ia
! Register prognostic variables
      call self%register_state_variable(self%id_ph1,'ph1','umol/L','Small phytoplankton (Flagellates)',minimum=0.0_rk) !,initial_value=ph1_0,
      call self%register_state_variable(self%id_ph2,'ph2','umol/L','Large phytoplankton (Diatoms)',minimum=0.0_rk) !initial_value=ph2_0
      call self%register_state_variable(self%id_zo1,'zo1','umol/L','Microzooplankton',minimum=0.0_rk) !initial_value=zo1_0
      call self%register_state_variable(self%id_zo2,'zo2','umol/L','Mesozooplankton',minimum=0.0_rk)    !initial_value=zo2_0
      call self%register_state_variable(self%id_no3,'no3','umol/L','Nitrate',minimum=0.0_rk)      !initial_value=no3_0,                    
      call self%register_state_variable(self%id_nh4,'nh4','umol/L','Ammonium',minimum=0.0_rk) !initial_value=nh4_0
      call self%register_state_variable(self%id_de1,'de1','umol/L','Small Detritus',minimum=0.0_rk,vertical_movement=self%w1)      !initial_value=de1_0                    
      call self%register_state_variable(self%id_de2,'de2','umol/L','Large Detritus',minimum=0.0_rk,vertical_movement=self%w2)        !initial_value=de2_0                  
      call self%register_state_variable(self%id_bsi,'bsi','umol/L','Biogenic Silica',minimum=0.0_rk,vertical_movement=self%w2)    !initial_value=bsi_0                      
      call self%register_state_variable(self%id_sil,'sil','umol/L','Silicate',minimum=0.0_rk)         !initial_value=sil_0                 
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,self%id_ph1,scale_factor=self%ac)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,self%id_ph2,scale_factor=self%ac)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,self%id_de1,scale_factor=self%ac)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,self%id_de2,scale_factor=self%ac)
! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_chl,'chl','mg-chla m-3','Chlorophyll')
      call self%register_diagnostic_variable(self%id_chl1,'chl1','mg-chla m-3','Chlorophyll in ph1')
      call self%register_diagnostic_variable(self%id_chl2,'chl2','mg-chla m-3','Chlorophyll in ph2')
      call self%register_diagnostic_variable(self%id_fnutph1,'fnutph1','umol/L/day','no3,nh4-->ph1 (Primary Production)')
      call self%register_diagnostic_variable(self%id_fnutph2,'fnutph2','umol/L/day','no3,nh4-->ph2 (Primary production)')
      call self%register_diagnostic_variable(self%id_fph1nh4,'fph1nh4','umol/L/day','ph1-->nh4 (ph1 mortality)')
      call self%register_diagnostic_variable(self%id_fph2nh4,'fph2nh4','umol/L/day','ph2-->nh4 (ph2 mortality)')
      call self%register_diagnostic_variable(self%id_fph2sil,'fph2sil','umol/L/day','ph2-->sil (ph2 mortality)')
      call self%register_diagnostic_variable(self%id_fph1zo1,'fph1zo1','umol/L/day','ph1-->zo1 (Ingestion)')
      call self%register_diagnostic_variable(self%id_fph2zo2,'fph2zo2','umol/L/day','ph2-->zo2 (Ingestion)')
      call self%register_diagnostic_variable(self%id_fph2de2,'fph2de2','umol/L/day','ph2-->de2 (Aggregation)')
      call self%register_diagnostic_variable(self%id_fph2bsi,'fph2bsi','umol/L/day','ph2-->bsi (Aggregation)')
      call self%register_diagnostic_variable(self%id_fzo1nh4,'fzo1nh4','umol/L/day','zo1-->nh4 (Excretion)')
      call self%register_diagnostic_variable(self%id_fzo1zo2,'fzo1zo2','umol/L/day','zo1-->zo2 (Ingestion)')
      call self%register_diagnostic_variable(self%id_fzo2htl,'fzo2htl','umol/L/day','zo2-->htl (Closure)')
      call self%register_diagnostic_variable(self%id_fnh4phy,'fnh4phy','umol/L/day','nh4-->ph1,ph2 (nh4 uptake)')
      call self%register_diagnostic_variable(self%id_fzo2nh4,'fzo2nh4','umol/L/day','zo2-->nh4 (Excretion)')
      call self%register_diagnostic_variable(self%id_fde1nh4,'fde1nh4','umol/L/day','de1-->nh4 (Remineralization)')
      call self%register_diagnostic_variable(self%id_fde2nh4,'fde2nh4','umol/L/day','de2-->nh4 (Remineralization)')
      call self%register_diagnostic_variable(self%id_fnh4no3,'fnh4no3','umol/L/day','nh4-->no3 (Nitrification)')
      call self%register_diagnostic_variable(self%id_fno3phy,'fno3phy','umol/L/day','no3-->ph1,ph2 (no3 uptake)')
      call self%register_diagnostic_variable(self%id_foo1de1,'foo1de1','umol/L/day','food1-->de1 (Sloppy feeding/Egestion)')
      call self%register_diagnostic_variable(self%id_foo2de2,'foo2de2','umol/L/day','food2-->de2 (Sloppy feeding/Egestion)')
      call self%register_diagnostic_variable(self%id_fde1zo1,'fde1zo1','umol/L/day','de1-->zo1 (Grazing)')
      call self%register_diagnostic_variable(self%id_fde2zo2,'fde2zo2','umol/L/day','de2-->zo2 (Grazing)')
      call self%register_diagnostic_variable(self%id_fbsisil,'fbsisil','umol/L/day','bsi-->sil (Remineralization)')
      call self%register_diagnostic_variable(self%id_fsilph2,'fsilph2','umol/L/day','sil-->ph2 (Uptake)')
      call self%register_diagnostic_variable(self%id_foo1zo1,'foo1zo1','umol/L/day','ph1,det-->zo1 (zo1 grazing)')
      call self%register_diagnostic_variable(self%id_foo2zo2,'foo2zo2','umol/L/day','ph2,zo1-->zo2 (zo2 grazing)')
      call self%register_diagnostic_variable(self%id_lf1,'lf1','-','Light limitation term for ph1 growth')
      call self%register_diagnostic_variable(self%id_lf2,'lf2','-','Light limitation term for ph2 growth')
      call self%register_diagnostic_variable(self%id_nf,'nf','-','no3+nh4 limitation term for ph1&ph2 growth')
      call self%register_diagnostic_variable(self%id_sf,'sf','-','sil limitation term for ph2 growth')
      call self%register_diagnostic_variable(self%id_f1f,'f1f','-','Iron limitation term for ph1 growth')
      call self%register_diagnostic_variable(self%id_f2f,'f2f','-','Iron limitation term for ph2 growth')
      call self%register_diagnostic_variable(self%id_par_diag,'PAR','W/m2','Photosynthetically Available Radiation')
      call self%register_horizontal_diagnostic_variable(self%id_stemp,'stemp','deg.C','Temperature in the top layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialde2,'fialde2','mmol m-2 d-1','ial --> de2 (ice algal flux',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialbsi,'fialbsi','mmol m-2 d-1','ial --> bsi (ice algal flux',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialph2,'fialph2','mmol m-2 d-1','ial --> ph2 (ice algal flux',source=source_do_horizontal)
! Register environmental variables
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
! Register horizontal environmental variables
     ! if(any(models.eq.'uvic_icealgae'))then
      if (self%use_icealgae) then
       call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)
       call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt)
       call self%register_dependency(self%id_ia,'uvic_icealgae_ia','','')
       call self%register_dependency(self%id_iceno3,'uvic_icealgae_no3','','')
       call self%register_dependency(self%id_icenh4,'uvic_icealgae_nh4','','')
       call self%register_dependency(self%id_icesil,'uvic_icealgae_sil','','')
       call self%register_dependency(self%id_fmelt,'uvic_icealgae_fmelt','','')
       call self%register_dependency(self%id_fpond,'uvic_icealgae_fpond','','')
       call self%register_dependency(self%id_fpondno3,'uvic_icealgae_fpondno3','','')
       call self%register_dependency(self%id_fpondnh4,'uvic_icealgae_fpondnh4','','')
       call self%register_dependency(self%id_fpondsil,'uvic_icealgae_fpondsil','','')
       call self%register_dependency(self%id_fmort,'uvic_icealgae_fmort','','')
       call self%register_dependency(self%id_fmort2,'uvic_icealgae_fmort2','','')
       call self%register_dependency(self%id_fskelno3,'uvic_icealgae_fskelno3','','')
       call self%register_dependency(self%id_fskelnh4,'uvic_icealgae_fskelnh4','','')
       call self%register_dependency(self%id_fskelsil,'uvic_icealgae_fskelsil','','')
       call self%request_coupling(self%id_ia,'uvic_icealgae_ia')
       call self%request_coupling(self%id_iceno3,'uvic_icealgae_no3')
       call self%request_coupling(self%id_icenh4,'uvic_icealgae_nh4')
       call self%request_coupling(self%id_icesil,'uvic_icealgae_sil')
       call self%request_coupling(self%id_fmelt,'uvic_icealgae_fmelt')
       call self%request_coupling(self%id_fmort,'uvic_icealgae_fmort')
       call self%request_coupling(self%id_fmort2,'uvic_icealgae_fmort2')
       call self%request_coupling(self%id_fpond,'uvic_icealgae_fpond')
       call self%request_coupling(self%id_fpondno3,'uvic_icealgae_fpondno3')
       call self%request_coupling(self%id_fpondnh4,'uvic_icealgae_fpondnh4')
       call self%request_coupling(self%id_fpondsil,'uvic_icealgae_fpondsil')
       call self%request_coupling(self%id_fskelno3,'uvic_icealgae_fskelno3')
       call self%request_coupling(self%id_fskelnh4,'uvic_icealgae_fskelnh4')
       call self%request_coupling(self%id_fskelsil,'uvic_icealgae_fskelsil')
       call self%register_dependency(self%id_density,standard_variables%density)
      endif
   end subroutine initialize
   
   subroutine do(self,_ARGUMENTS_DO_)

   class (type_uvic_eco),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
! Declare prognostic variables
   real(rk) :: ph1,ph2,zo1,zo2,no3,nh4,de1,de2,bsi,sil
! Declare diagnostic variables   
   real(rk) :: foo2bsi,fnutph1,fnutph2,fph1nh4,fph2nh4,fph2sil,fph1zo1,fph2zo2,fph2de2,fph2bsi,fzo1nh4,fzo1zo2,fzo2htl,fnh4phy,fzo2nh4,fde1nh4,fde2nh4,fnh4no3,fno3phy,foo1de1,foo2de2,fde1zo1,fde2zo2,fbsisil,fsilph2,foo1zo1,foo2zo2
! Declare environmental variables
   real(rk) :: temp,par,ia
! Declare anything else used in this subroutine
   real(rk) :: pht,nut,texp,mu1q,mu2q,mp1q,mp2q,gz1q,mz1q,mz2q,gz2q,rd1q,rd2q,rd3q,food1,food2,graz1,graz2,lf1,lf2,nf,sf,lim1,lim2,ppt,rpi

   _LOOP_BEGIN_

! Retrieve prognostic variables
   _GET_(self%id_ph1,ph1)
   _GET_(self%id_ph2,ph2)
   _GET_(self%id_zo1,zo1)
   _GET_(self%id_zo2,zo2)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_de1,de1)
   _GET_(self%id_de2,de2)
   _GET_(self%id_bsi,bsi)
   _GET_(self%id_sil,sil)
! Retrieve environmental variables
   _GET_(self%id_temp,temp)
   _GET_(self%id_par,par)
   !if(any(self%models.eq.'uvic_icealgae'))then
   if (self%use_icealgae) then
    _GET_HORIZONTAL_(self%id_ia,ia)
    par = par * exp(-self%ac_ia*ia*self%zia)
   endif
   pht = ph1 + ph2 ! Total phytoplankton biomass
   nut = nh4 + no3 ! Total nutrient concentration
   texp = (temp-10._rk)/10._rk
   mu1q  = self%mu1*self%qp**texp ! ph1 growth rate
   mu2q  = self%mu2*self%qp**texp ! ph2 growth rate
   mp1q = self%mp1*self%qp**texp ! ph1 exretion rate
   mp2q = self%mp2*self%qp**texp ! ph2 exretion rate
   gz1q  = self%gz1*self%qz**texp ! zo1 grazing rate
   mz1q = self%mz1*self%qz**texp ! zo1 excretion rate
   mz2q = self%mz2*self%qz**texp ! zo2 excretion rate
   gz2q  = self%gz2*self%qz**texp ! zo2 grazing rate
   rd1q = self%rd1*self%qb**texp ! de1 remineralization rate
   rd2q = self%rd2*self%qb**texp ! de2 remineralization rate
   rd3q = self%rd3*self%qb**texp ! bsi remineralization rate
   food1 = self%pp1*ph1 + self%pd1*de1 ! Food source for zo1
   food2 = self%pp2*ph2 + self%pz1*zo1 + self%pd2*de2 ! Food source for zo2
   graz1 = gz1q*zo1*(food1**2)/((self%kz1**2) + (food1**2))
   graz2 = gz2q*zo2*(food2**2)/((self%kz2**2) + (food2**2))
   fph1zo1 = graz1*self%pp1*ph1/food1                    
   fde1zo1 = graz1*self%pd1*de1/food1
   fde2zo2 = graz2*self%pd2*de2/food2
   fph2zo2 = graz2*self%pp2*ph2/food2
   fzo1zo2 = graz2*self%pz1*zo1/food2
   lf1 = 1.0_rk - exp(-self%rpp1*par)
   lf2 = 1.0_rk - exp(-self%rpp2*par)
   nf = nut/(self%kn + nut)
   sf = sil/(self%ks + sil)
   lim1 = min(nf,lf1)
   lim2 = min(nf,sf,lf2)
   fnutph1 = lim1*mu1q*ph1
   fnutph2 = lim2*mu2q*ph2
   ppt = fnutph1 + fnutph2
   fph1nh4 = mp1q*ph1
   fph2nh4 = mp2q*ph2
   fph2sil = self%rsin*mp2q*ph2
   fph2de2 = self%agg*(ph2**2)
   fph2bsi = self%rsin*fph2de2
   fzo1nh4 = mz1q*zo1
   fzo2nh4 = mz2q*zo2
   foo1de1 = (1.0_rk-self%az1)*graz1
   foo2de2 = (1.0_rk-self%az2)*graz2
   foo2bsi = (1.0_rk-self%az2)*graz2*self%rsin   
   fde1nh4 = rd1q*de1
   fde2nh4 = rd2q*de2
   fbsisil = rd3q*bsi
   fsilph2 = self%rsin*fnutph2
   foo1zo1 = self%az1*graz1
   foo2zo2 = self%az2*graz2
   rpi = self%rpf/(self%rpf+nh4) ! relative preference index of phy for nh4 uptake relative to no3
   fno3phy = rpi*ppt*no3/nut
   fnh4phy = ppt - fno3phy
   fnh4no3 = self%rn0*exp(self%knt*temp)*nh4/(1+par)
   fzo2htl = self%rc*zo2**2
!   fhtlnh4 = self%htlnh4*fzo2htl
!   fhtldet = self%htldet*fzo2htl
!  if (ppt.gt.0.0) then ! pp reduced if fno3phy or fnh4phy set to zero
!   if ((fnh4phy.eq.0.0).or.(fno3phy.eq.0.0)) then
!    fnutph1=fnutph1*(fno3phy+fnh4phy)/ppt
!    fnutph2=fnutph2*(fno3phy+fnh4phy)/ppt
!    ppt=fno3phy+fnh4phy
!   endif
!  endif
      
! Compute and save prognostic variables
   if (ph1.lt.self%pmin) then
    _SET_ODE_(self%id_ph1,fnutph1)
   else
    _SET_ODE_(self%id_ph1,fnutph1-fph1zo1-fph1nh4)
   endif
   if (ph2.lt.self%pmin) then
    _SET_ODE_(self%id_ph2,fnutph2)
   else
    _SET_ODE_(self%id_ph2,fnutph2-fph2zo2-fph2nh4-fph2de2)
   endif
   if (zo1.lt.self%pmin) then
    _SET_ODE_(self%id_zo1,foo1zo1)
   else
    _SET_ODE_(self%id_zo1,foo1zo1-fzo1zo2-fzo1nh4)
   endif
   if (zo2.lt.self%pmin) then
    _SET_ODE_(self%id_zo2,foo2zo2)
   else
    _SET_ODE_(self%id_zo2,foo2zo2-fzo2nh4-fzo2htl)
   endif
   _SET_ODE_(self%id_no3,fnh4no3-fno3phy)
!  _SET_ODE_(self%id_nh4,fde1nh4+fde2nh4+fph1nh4+fph2nh4+fzo1nh4+fzo2nh4-fnh4no3-fnh4phy)
   _SET_ODE_(self%id_nh4,fde1nh4+fde2nh4+fzo1nh4+fzo2nh4-fnh4no3-fnh4phy)
!  _SET_ODE_(self%id_de1,foo1de1-fde1zo1-fde1nh4)
   _SET_ODE_(self%id_de1,foo1de1-fde1zo1-fde1nh4+fph1nh4)
!  _SET_ODE_(self%id_de2,foo2de2+fph2de2-fde2zo2-fde2nh4)
   _SET_ODE_(self%id_de2,foo2de2+fph2de2-fde2zo2-fde2nh4+fph2nh4)
   _SET_ODE_(self%id_bsi,foo2bsi+fph2bsi+self%rsin*fph2nh4-self%rsin*fde2zo2-fbsisil)
   _SET_ODE_(self%id_sil,-fsilph2+fph2sil+fbsisil)
! Save diagnostic variables  
   _SET_DIAGNOSTIC_(self%id_chl,ph1*0.795+ph2*3.533)
   _SET_DIAGNOSTIC_(self%id_chl1,ph1*0.795)
   _SET_DIAGNOSTIC_(self%id_chl2,ph2*3.533)
   _SET_DIAGNOSTIC_(self%id_fnutph1,fnutph1*self%spd)
   _SET_DIAGNOSTIC_(self%id_fnutph2,fnutph2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph1nh4,fph1nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2nh4,fph2nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2sil,fph2sil*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph1zo1,fph1zo1*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2zo2,fph2zo2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2de2,fph2de2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo1nh4,fzo1nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo1zo2,fzo1zo2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo2htl,fzo2htl*self%spd)
   _SET_DIAGNOSTIC_(self%id_fnh4phy,fnh4phy*self%spd)
   _SET_DIAGNOSTIC_(self%id_fzo2nh4,fzo2nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fde1nh4,fde1nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fde2nh4,fde2nh4*self%spd)
   _SET_DIAGNOSTIC_(self%id_fnh4no3,fnh4no3*self%spd)
   _SET_DIAGNOSTIC_(self%id_fno3phy,fno3phy*self%spd)
   _SET_DIAGNOSTIC_(self%id_foo1de1,foo1de1*self%spd)
   _SET_DIAGNOSTIC_(self%id_foo2de2,foo2de2*self%spd)
   _SET_DIAGNOSTIC_(self%id_fde1zo1,fde1zo1*self%spd)
   _SET_DIAGNOSTIC_(self%id_fde2zo2,fde2zo2*self%spd)
   _SET_DIAGNOSTIC_(self%id_foo2bsi,foo2bsi*self%spd)
   _SET_DIAGNOSTIC_(self%id_fph2bsi,fph2bsi*self%spd)   
   _SET_DIAGNOSTIC_(self%id_fbsisil,fbsisil*self%spd)
   _SET_DIAGNOSTIC_(self%id_fsilph2,fsilph2*self%spd)
   _SET_DIAGNOSTIC_(self%id_foo1zo1,foo1zo1*self%spd)
   _SET_DIAGNOSTIC_(self%id_foo2zo2,foo2zo2*self%spd)
   _SET_DIAGNOSTIC_(self%id_lf1,lf1)
   _SET_DIAGNOSTIC_(self%id_lf2,lf2)
   _SET_DIAGNOSTIC_(self%id_nf,nf)
   _SET_DIAGNOSTIC_(self%id_sf,sf)
   _SET_DIAGNOSTIC_(self%id_par_diag,par)

   _LOOP_END_

   end subroutine do

! Calling and saving the surface variables.
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uvic_eco),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk) :: botmelt,botgrowth,ia,iceno3,icenh4,icesil,fialde2,fialbsi,fialph2,stemp,fmelt,fpond,fpondno3,fpondnh4,fpondsil,fmort,fmort2,fskelno3,fskelnh4,fskelsil
   real(rk) :: ph1,ph2,zo1,zo2,no3,nh4,de1,de2,bsi,sil,density

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_temp,stemp)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_stemp,stemp)
  ! if(any(self%models.eq.'uvic_icealgae'))then
   if (self%use_icealgae) then
    _GET_HORIZONTAL_(self%id_botgrowth,botgrowth)
    _GET_HORIZONTAL_(self%id_botmelt,botmelt)
    _GET_HORIZONTAL_(self%id_ia,ia)
    _GET_HORIZONTAL_(self%id_iceno3,iceno3)
    _GET_HORIZONTAL_(self%id_icenh4,icenh4)
    _GET_HORIZONTAL_(self%id_icesil,icesil)
    _GET_(self%id_ph2,ph2)
    _GET_(self%id_no3,no3)
    _GET_(self%id_nh4,nh4)
    _GET_(self%id_sil,sil)
    _GET_(self%id_ph1,ph1)
    _GET_(self%id_zo1,zo1)
    _GET_(self%id_zo2,zo2)
    _GET_(self%id_de1,de1)
    _GET_(self%id_bsi,bsi)
    _GET_(self%id_density,density)
    _GET_HORIZONTAL_(self%id_fmelt,fmelt)
    _GET_HORIZONTAL_(self%id_fmort,fmort)
    _GET_HORIZONTAL_(self%id_fmort2,fmort2)
    _GET_HORIZONTAL_(self%id_fpond,fpond)
    _GET_HORIZONTAL_(self%id_fpondno3,fpondno3)
    _GET_HORIZONTAL_(self%id_fpondnh4,fpondnh4)
    _GET_HORIZONTAL_(self%id_fpondsil,fpondsil)
    _GET_HORIZONTAL_(self%id_fskelno3,fskelno3)
    _GET_HORIZONTAL_(self%id_fskelnh4,fskelnh4)
    _GET_HORIZONTAL_(self%id_fskelsil,fskelsil)
    ! [fmelt,fmort,fmort2,fskelno3,fskelnh4,fskelsil] = mmol m-3 d-1 in icealgae model, which needs to be converted to umol m-2 s-1 here for surface exchange equation.
    fmelt=fmelt/self%spd*self%zia
    fpond=fpond/self%spd*self%zia
    fpondno3=fpondno3/self%spd*self%zia
    fpondnh4=fpondnh4/self%spd*self%zia
    fpondsil=fpondsil/self%spd*self%zia
    fmort=fmort/self%spd*self%zia
    fmort2=fmort2/self%spd*self%zia
    fskelno3=fskelno3/self%spd*self%zia
    fskelnh4=fskelnh4/self%spd*self%zia
    fskelsil=fskelsil/self%spd*self%zia
    fialde2=min(0.0,(1-self%f_seed)*fmelt)-(1-self%f_seed)*fpond
    fialbsi=(min(0.0,(1-self%f_seed)*fmelt)-(1-self%f_seed)*fpond)*self%rsin
    fialph2=max(0.0,fmelt)+min(0.0,self%f_seed*fmelt)-self%f_seed*fpond
    _SET_SURFACE_EXCHANGE_(self%id_de2,-fialde2+0.7*fmort+fmort2-de2*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*de2)
    _SET_SURFACE_EXCHANGE_(self%id_bsi,-fialbsi+self%rsin*(0.7*fmort+fmort2)-bsi*(913./density*botmelt+abs(fpondsil/icesil))+913./density*botgrowth*bsi)
!HH0: conc/dilution effect
!    _SET_SURFACE_EXCHANGE_(self%id_ph2,-fialph2)
!    _SET_SURFACE_EXCHANGE_(self%id_no3,-fskelno3+fpondno3)
!    _SET_SURFACE_EXCHANGE_(self%id_nh4,-fskelnh4+fpondnh4)
!    _SET_SURFACE_EXCHANGE_(self%id_sil,-fskelsil+fpondsil)
    _SET_SURFACE_EXCHANGE_(self%id_ph2,(1000./density*ia-ph2)*913./1000.*abs(fialph2)/ia+913./density*botgrowth*ph2)
    _SET_SURFACE_EXCHANGE_(self%id_no3,-fskelno3-(no3-1000./density*iceno3)*(913./density*botmelt+1000./density*abs(fpondno3/iceno3))+913./density*botgrowth*(no3-iceno3))
    _SET_SURFACE_EXCHANGE_(self%id_nh4,-fskelnh4-(nh4-1000./density*icenh4)*(913./density*botmelt+1000./density*abs(fpondnh4/icenh4))+913./density*botgrowth*(nh4-icenh4))
    _SET_SURFACE_EXCHANGE_(self%id_sil,-fskelsil-(sil-1000./density*icesil)*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*(sil-icesil))
    _SET_SURFACE_EXCHANGE_(self%id_ph1,-ph1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*ph1)
    _SET_SURFACE_EXCHANGE_(self%id_de1,-de1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*de1)
    _SET_SURFACE_EXCHANGE_(self%id_zo1,-zo1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*zo1)
    _SET_SURFACE_EXCHANGE_(self%id_zo2,-zo2*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*zo2)
!HH1
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialde2,fialde2*self%spd)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialbsi,fialbsi*self%spd)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialph2,fialph2*self%spd)

   endif
   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
end module uvic_eco
