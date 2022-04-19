!jpnote: my copy and paste of this file from mortenson
!-----------------------------------------------------
!
! University of Victoria Ecosystem model
!
! Contributors: Hakase Hayashida, Nadja Steiner, Eric Mortenson, Adam Mohanan
!
! Model description: The model will be described in a paper.

! JPnote: eventually add model description
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
      !mortenson 
      !type (type_dependency_id) :: id_temp,id_par
     !hayashida ! 
      type (type_dependency_id) :: id_temp,id_par,id_density
! Declare horizontal environmental variables
      !mortenson
      !type (type_horizontal_dependency_id) :: id_ia,id_fmelt,id_fpond,id_fpondno3,id_fpondnh4,id_fpondsil,id_fmort,id_fmort2,id_fskelno3,id_fskelnh4,id_fskelsil
     ! hayashida  !
      type (type_horizontal_dependency_id) :: id_botmelt,id_botgrowth,id_ia,id_iceno3,id_icenh4,id_icesil,id_fmelt,id_fpond,id_fpondno3,id_fpondnh4,id_fpondsil,id_fmort,id_fmort2,id_fskelno3,id_fskelnh4,id_fskelsil
! Declare namelist parameters jpnote: yaml paramteters 
      real(rk) :: ac,f_seed,ph1_0,ph2_0,zo1_0,zo2_0,no3_0,nh4_0,de1_0,de2_0,bsi_0,sil_0,w1,w2,mu1,mu2,kn,rpp1,rpp2,mp1,mp2,gz1,kz1,az1,az2,mz1,rc,pp1,pp2,pd1,pd2,pz1,gz2,kz2,mz2,rd1,rd2,rd3,rpf,rn0,knt,qp,qz,qb,agg,rsin,ks,pmin
! Declare anything else used in all procedures
      real(rk) :: spd = 86400.0_rk ! Seconds Per Day (spd) 
      real(rk) :: zia,ac_ia 
      logical :: use_icealgae !jpnote added 

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

   real(rk) :: r_pond,fmethod,fflush,drag,f_graze,zia,ac_ia,ia_0,ia_b,rnit,skno3_0,sknh4_0,sksil_0,ks_no3,ks_sil,maxg,mort,mort2,crit_melt,lcompp,rpp,rpi,t_sens,nu,md_no3,md_sil,chl2n,sil2n
   real(rk) :: ac,f_seed,ph1_0,ph2_0,zo1_0,zo2_0,no3_0,nh4_0,de1_0,de2_0,bsi_0,sil_0,w1,w2,mu1,mu2,kn,rpp1,rpp2,mp1,mp2,gz1,kz1,az1,az2,mz1,rc,pp1,pp2,pd1,pd2,pz1,gz2,kz2,mz2,rd1,rd2,rd3,rpf,rn0,knt,qp,qz,qb,agg,rsin,ks,pmin


!read in from yaml - based on values from fabm.nml 
   !uvic_eco
   call self%get_parameter(self%use_icealgae, 'use_icealgae', '', 'use icealgae', default=.false.)
   call self%get_parameter(self%ac,'ac','m-1','light attenuation coefficient', default=0.03_rk)
   call self%get_parameter(self%f_seed, 'f_seed','-', 'fraction of ice algal fux as ph2 seeding', default=0.0_rk)
!jpnote: done under initliazation and then register_state_variable)
  ! call self%get_parameter(ph1_0, 'ph1_0','umol/L','ph1 initial value', default=1.0_rk )
  ! call self%get_parameter(ph2_0 , 'ph2_0 ','umol/L', 'ph2 initial value', default=0.5_rk)
  ! call self%get_parameter(zo1_0, 'zo1_0','umol/L', 'zo1 initial value', default=0.2_rk)
  !! call self%get_parameter(zo2_0, 'zo2_0','umol/L', 'zo2 initial value', default=0.1_rk)
   !call self%get_parameter(nh4_0, 'nh4_0','umol/L', 'nh4 initial value', default=10.0_rk)
   !call self%get_parameter(no3_0, 'no3_0','umol/L', 'no3 initial value', default=10.0_rk)
   !call self%get_parameter(de1_0, 'de1_0','umol/L','de1 initial value ', default=1.0_rk)
   !call self%get_parameter(de2_0, 'de2_0','umol/L', 'de2 initial value', default=1.0_rk)
   !call self%get_parameter(bsi_0, 'bsi_0','umol/L', 'bsi initial value', default=1.0_rk)
   !call self%get_parameter(sil_0 , 'sil_0 ','umol/L', 'sil initial value', default=5.0_rk)
  
   call self%get_parameter(self%w1, 'w1 ','m/d', 'de1 sinking rate',default=6.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%w2, 'w2 ','m/d', 'de2 sinking rate', default=6.0_rk, scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%mu1, 'mu1','1/d', 'ph1 maximum growth rate', default=2.0_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%mu2, 'mu2 ','1/d', 'ph2 maximum growth rate', default=2.0_rk, scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%kn, 'kn','umol/L', ' no3 & nh4 half saturation constant', default=0.1_rk)
   !call self%get_parameter(self%al1, 'al1','m2/W/time', 'initial slope of P-I curve ([time] is same unit as pmax)', default=0)
   !call self%get_parameter(self%al2 , 'al2 ','m2/W/time', 'initial slope of P-I curve ([time] is same unit as pmax)', default=0)
   call self%get_parameter(self%rpp1, 'rpp1','1/time', 'maximum photosynthetic rate ([time] is same unit as alpha)', default=0.05_rk)
   call self%get_parameter(self%rpp2, 'rpp2 ','1/time', 'maximum photosynthetic rate ([time] is same unit as alpha)', default=0.05_rk)
   call self%get_parameter(self%mp1, 'mp1','1/d', 'ph1 excretion rate ', default=0.05_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%mp2, 'mp2 ','1/d', 'ph2 excretion rate', default=0.05_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%mz1, 'mz1','1/d', 'zo1 excretion rate ', default=0.1_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%mz2, 'mz2','1/d', 'zo2 excretion rate ', default=0.3_rk,scale_factor=1.0_rk/self%spd)
 
   call self%get_parameter(self%gz1, 'gz1','1/d', 'zo1 maximum grazing rate ', default=1.3_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%gz2, 'gz2','1/d', 'zo2 maximum grazing rate ', default=0.8_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%kz1, 'kz1','umol/L', 'zo1 grazing half saturation constant ',default=0.6_rk)
   call self%get_parameter(self%kz2, 'kz2','umol/L', 'zo2 grazing half saturation constant ', default=0.75_rk)
   call self%get_parameter(self%az1, 'az1','-', 'zo1 assimilation fraction ',default=0.7_rk)
   call self%get_parameter(self%az2, 'az2','-', 'zo2 assimilation fraction ',default=0.1_rk)
   call self%get_parameter(self%rc, 'rc','1/d', 'closure mortality rate', default=0.003_rk,scale_factor=1.0_rk/self%spd)
   ! call self%get_parameter(self%htlnh4, 'htlnh4','-', 'closure loss fraction to nh4', default=0)
  ! call self%get_parameter(self%htldet, 'htldet','-', 'closure loss fraction to det', default=0)
   call self%get_parameter(self%pp1, 'pp1','-', 'ph1 fraction as food for zo1', default=1.0_rk)
   call self%get_parameter(self%pp2, 'pp2','-', 'ph2 fraction as food for zo2', default=1.0_rk)
   call self%get_parameter(self%pd1, 'pd1','-', 'de1 fraction as food for zo1', default=0.5_rk)
   call self%get_parameter(self%pd2, 'pd2','-', 'de2 fraction as food for zo2', default=0.5_rk)
   call self%get_parameter(self%pz1, 'pz1','-', 'zo1 fraction as food for zo2', default=1.0_rk)
   call self%get_parameter(self%rd1, 'rd1','1/d', ' de1 remineralization rate', default=0.1_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%rd2, 'rd2','1/d', 'de2 remineralization rate', default=0.1_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%rd3, 'rd3','1/d', 'bsi remineralization rate', default=0.1_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%rpf, 'rpf','-', 'scale factor for nitrogen preference function', default=0.2_rk )
   call self%get_parameter(self%rn0, 'rn0','1/d', 'nitrification rate at 0 deg.C', default=0.03_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%knt, 'knt','1/deg.C', 'nitrification temperature coefficient', default=0.0693_rk)
   call self%get_parameter(self%qp, 'qp','-', 'ph1 & ph2 Q10 factor', default=2.0_rk)
   call self%get_parameter(self%qz, 'qz','-', 'zo1 & zo2 Q10 factor', default=3.0_rk )
   call self%get_parameter(self%qb, 'qb','-', 'bacteria Q10 factor', default=3.0_rk)
   call self%get_parameter(self%agg, 'agg', '1/d', 'ph2 aggregation rate',default=0.07_rk,scale_factor=1.0_rk/self%spd)
   call self%get_parameter(self%rsin, 'rsin','mol-Si/mol-N', 'ph2 Si:N ratio', default=2.0_rk)
   call self%get_parameter(self%ks, 'ks','umol/L', 'si half saturation constant', default=2.0_rk)
   call self%get_parameter(self%pmin, 'pmin','umol-N/L', 'background plankton concentration', default=0.01_rk)

   if (self%use_icealgae) then!read in icealgae model vars 
      call self%get_parameter(self%zia, 'zia','m', 'ice algal layer thickness', default=0.03_rk) ! zia = 0.03_rk
      call self%get_parameter(self%ac_ia, 'ac_ia','', 'specific light attenuation coefficient for ice algae', default=0.007_rk) !ac_ia = 0.007_rk
   endif


! Register prognostic variables
!jpnote these register_state_vars replace the yaml read? put them under intiliazation: 
#if 0
      call self%register_state_variable(self%id_ph1,'ph1','umol/L','Small phytoplankton (Flagellates)',initial_value=ph1_0,minimum=0.0_rk) !jpnote change initial values to self%
      call self%register_state_variable(self%id_ph2,'ph2','umol/L','Large phytoplankton (Diatoms)',initial_value=ph2_0,minimum=0.0_rk)
      call self%register_state_variable(self%id_zo1,'zo1','umol/L','Microzooplankton',initial_value=zo1_0,minimum=0.0_rk)
      call self%register_state_variable(self%id_zo2,'zo2','umol/L','Mesozooplankton',initial_value=zo2_0,minimum=0.0_rk)   
      call self%register_state_variable(self%id_no3,'no3','umol/L','Nitrate',initial_value=no3_0,minimum=0.0_rk)                         
      call self%register_state_variable(self%id_nh4,'nh4','umol/L','Ammonium',initial_value=nh4_0,minimum=0.0_rk)
      call self%register_state_variable(self%id_de1,'de1','umol/L','Small Detritus',initial_value=de1_0,minimum=0.0_rk,vertical_movement=self%w1)                         
      call self%register_state_variable(self%id_de2,'de2','umol/L','Large Detritus',initial_value=de2_0,minimum=0.0_rk,vertical_movement=self%w2)                         
      call self%register_state_variable(self%id_bsi,'bsi','umol/L','Biogenic Silica',initial_value=bsi_0,minimum=0.0_rk,vertical_movement=self%w2)                         
      call self%register_state_variable(self%id_sil,'sil','umol/L','Silicate',initial_value=sil_0,minimum=0.0_rk)  

#endif
!initial values in yaml --> 
      call self%register_state_variable(self%id_ph1,'ph1','umol/L','Small phytoplankton (Flagellates)',minimum=0.0_rk) !jpnote change initial values to self%
      call self%register_state_variable(self%id_ph2,'ph2','umol/L','Large phytoplankton (Diatoms)',minimum=0.0_rk)
      call self%register_state_variable(self%id_zo1,'zo1','umol/L','Microzooplankton',minimum=0.0_rk)
      call self%register_state_variable(self%id_zo2,'zo2','umol/L','Mesozooplankton',minimum=0.0_rk)   
      call self%register_state_variable(self%id_no3,'no3','umol/L','Nitrate',minimum=0.0_rk)                         
      call self%register_state_variable(self%id_nh4,'nh4','umol/L','Ammonium',minimum=0.0_rk)
      call self%register_state_variable(self%id_de1,'de1','umol/L','Small Detritus',minimum=0.0_rk,vertical_movement=self%w1)                         
      call self%register_state_variable(self%id_de2,'de2','umol/L','Large Detritus',minimum=0.0_rk,vertical_movement=self%w2)                         
      call self%register_state_variable(self%id_bsi,'bsi','umol/L','Biogenic Silica',minimum=0.0_rk,vertical_movement=self%w2)                         
      call self%register_state_variable(self%id_sil,'sil','umol/L','Silicate',minimum=0.0_rk)  

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
      call self%register_diagnostic_variable(self%id_fbsisil,'fbsisil','umol/L/day','bsi-->sil (Remineralization)') !**
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
      call self%register_horizontal_diagnostic_variable(self%id_stemp,'stemp','degC','Temperature in the top layer',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialde2,'fialde2','mmol m-2 d-1','ial --> de2 (ice algal flux)',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialbsi,'fialbsi','mmol m-2 d-1','ial --> bsi (ice algal flux)',source=source_do_horizontal)
      call self%register_horizontal_diagnostic_variable(self%id_fialph2,'fialph2','mmol m-2 d-1','ial --> ph2 (ice algal flux)',source=source_do_horizontal)
! Register environmental variables
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

! Register horizontal environmental variables
      !if(any(models.eq.'uvic_icealgae'))then jpnote changed to 
      if (self%use_icealgae) then 

      !hayashida  !
       call self%register_horizontal_dependency(self%id_botgrowth,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_grow)
       call self%register_horizontal_dependency(self%id_botmelt,standard_variables%tendency_of_sea_ice_thickness_due_to_thermodynamics_melt) 
      ! call self%register_dependency(self%id_ia,'uvic_icealgae_ia','','') !repeat do not need
       call self%register_dependency(self%id_iceno3,'uvic_icealgae_no3','','')
       call self%register_dependency(self%id_icenh4,'uvic_icealgae_nh4','','')
       call self%register_dependency(self%id_icesil,'uvic_icealgae_sil','','')

       call self%register_dependency(self%id_ia,'uvic_icealgae_ia','','')
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

       !hayashida !
       call self%request_coupling(self%id_iceno3,'uvic_icealgae_no3')
       !hayashida !
       call self%request_coupling(self%id_icenh4,'uvic_icealgae_nh4')
       !hayashida !
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
       !hayashida !
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
  ! if(any(self%models.eq.'uvic_icealgae'))then jpnote changed to 
   if(self%use_icealgae) then 
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

   ! jpnote _SET_ODE_ needs to be changed 
   if (ph1.lt.self%pmin) then
   ! _SET_ODE_(self%id_ph1,fnutph1)
    _ADD_SOURCE_(self%id_ph1,fnutph1)
   else
   ! _SET_ODE_(self%id_ph1,fnutph1-fph1zo1-fph1nh4)
    _ADD_SOURCE_(self%id_ph1,fnutph1-fph1zo1-fph1nh4)
   endif
   if (ph2.lt.self%pmin) then
   ! _SET_ODE_(self%id_ph2,fnutph2)
    _ADD_SOURCE_(self%id_ph2,fnutph2)
   else
   ! _SET_ODE_(self%id_ph2,fnutph2-fph2zo2-fph2nh4-fph2de2)
    _ADD_SOURCE_(self%id_ph2,fnutph2-fph2zo2-fph2nh4-fph2de2)
   endif
   if (zo1.lt.self%pmin) then
   ! _SET_ODE_(self%id_zo1,foo1zo1)
    _ADD_SOURCE_(self%id_zo1,foo1zo1)
   else
  !  _SET_ODE_(self%id_zo1,foo1zo1-fzo1zo2-fzo1nh4)
    _ADD_SOURCE_(self%id_zo1,foo1zo1-fzo1zo2-fzo1nh4)
   endif
   if (zo2.lt.self%pmin) then
   ! _SET_ODE_(self%id_zo2,foo2zo2)
    _ADD_SOURCE_(self%id_zo2,foo2zo2)
   else
   ! _SET_ODE_(self%id_zo2,foo2zo2-fzo2nh4-fzo2htl)
    _ADD_SOURCE_(self%id_zo2,foo2zo2-fzo2nh4-fzo2htl)
   endif


  !_SET_ODE_(self%id_no3,fnh4no3-fno3phy) 
   _ADD_SOURCE_(self%id_no3,fnh4no3-fno3phy)
!  _SET_ODE_(self%id_nh4,fde1nh4+fde2nh4+fph1nh4+fph2nh4+fzo1nh4+fzo2nh4-fnh4no3-fnh4phy)
    ! _SET_ODE_(self%id_nh4,fde1nh4+fde2nh4+fzo1nh4+fzo2nh4-fnh4no3-fnh4phy)
   _ADD_SOURCE_(self%id_nh4,fde1nh4+fde2nh4+fzo1nh4+fzo2nh4-fnh4no3-fnh4phy)
!  _SET_ODE_(self%id_de1,foo1de1-fde1zo1-fde1nh4)
     !_SET_ODE_(self%id_de1,foo1de1-fde1zo1-fde1nh4+fph1nh4)
   _ADD_SOURCE_(self%id_de1,foo1de1-fde1zo1-fde1nh4+fph1nh4)
!  _SET_ODE_(self%id_de2,foo2de2+fph2de2-fde2zo2-fde2nh4)
   ! _SET_ODE_(self%id_de2,foo2de2+fph2de2-fde2zo2-fde2nh4+fph2nh4)
   _ADD_SOURCE_(self%id_de2,foo2de2+fph2de2-fde2zo2-fde2nh4+fph2nh4)
   !_SET_ODE_(self%id_bsi,foo2bsi+fph2bsi+self%rsin*fph2nh4-self%rsin*fde2zo2-fbsisil)
   _ADD_SOURCE_(self%id_bsi,foo2bsi+fph2bsi+self%rsin*fph2nh4-self%rsin*fde2zo2-fbsisil)
    ! _SET_ODE_(self%id_sil,-fsilph2+fph2sil+fbsisil)
   _ADD_SOURCE_(self%id_sil,-fsilph2+fph2sil+fbsisil)


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
   !mortenson
   !real(rk) :: fialde2,fialbsi,fialph2,stemp,fmelt,fpond,fpondno3,fpondnh4,fpondsil,fmort,fmort2,fskelno3,fskelnh4,fskelsil
  
   !hayashida !
   real(rk) :: botmelt,botgrowth,ia,iceno3,icenh4,icesil,fialde2,fialbsi,fialph2,stemp,fmelt,fpond,fpondno3,fpondnh4,fpondsil,fmort,fmort2,fskelno3,fskelnh4,fskelsil
   !hayashida !
   real(rk) :: ph1,ph2,zo1,zo2,no3,nh4,de1,de2,bsi,sil,density

   logical :: use_icealgae  !jpnote 

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_temp,stemp)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_stemp,stemp)
  ! if(any(self%models.eq.'uvic_icealgae'))then   jpnote : changed to 
   if(self%use_icealgae) then 

   !hayashida !

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

!mortenson --> 

!   _SET_SURFACE_EXCHANGE_(self%id_de2,-fialde2+(1-self%f_seed)*fmort2)
    !_SET_SURFACE_EXCHANGE_(self%id_de2,-fialde2+0.7*fmort+fmort2) !jpnote change set bottom/surface exchange to addbottom/surface flux
    _ADD_SURFACE_FLUX_(self%id_de2,-fialde2+0.7*fmort+fmort2)
!   _SET_SURFACE_EXCHANGE_(self%id_bsi,-fialbsi)
    !_SET_SURFACE_EXCHANGE_(self%id_bsi,-fialbsi+self%rsin*(0.7*fmort+fmort2))
    _ADD_SURFACE_FLUX_(self%id_bsi,-fialbsi+self%rsin*(0.7*fmort+fmort2))
!   _SET_SURFACE_EXCHANGE_(self%id_ph2,-fialph2+self%f_seed*fmort2)
    !_SET_SURFACE_EXCHANGE_(self%id_ph2,-fialph2)
    _ADD_SURFACE_FLUX_(self%id_ph2,-fialph2)
    !_SET_SURFACE_EXCHANGE_(self%id_no3,-fskelno3+fpondno3)
    _ADD_SURFACE_FLUX_(self%id_no3,-fskelno3+fpondno3)
    !_SET_SURFACE_EXCHANGE_(self%id_nh4,-fskelnh4+fpondnh4)
    _ADD_SURFACE_FLUX_(self%id_nh4,-fskelnh4+fpondnh4)
    !_SET_SURFACE_EXCHANGE_(self%id_sil,-fskelsil+fpondsil)
    _ADD_SURFACE_FLUX_(self%id_sil,-fskelsil+fpondsil)

!hayashida --> 
#if 0
    _SET_SURFACE_EXCHANGE_(self%id_de2,-fialde2+0.7*fmort+fmort2-de2*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*de2)
    _SET_SURFACE_EXCHANGE_(self%id_bsi,-fialbsi+self%rsin*(0.7*fmort+fmort2)-bsi*(913./density*botmelt+abs(fpondsil/icesil))+913./density*botgrowth*bsi)
!!HH0: conc/dilution effect
!!    _SET_SURFACE_EXCHANGE_(self%id_ph2,-fialph2)
!!    _SET_SURFACE_EXCHANGE_(self%id_no3,-fskelno3+fpondno3)
!!   _SET_SURFACE_EXCHANGE_(self%id_nh4,-fskelnh4+fpondnh4)
!!    _SET_SURFACE_EXCHANGE_(self%id_sil,-fskelsil+fpondsil)
    _SET_SURFACE_EXCHANGE_(self%id_ph2,(1000./density*ia-ph2)*913./1000.*abs(fialph2)/ia+913./density*botgrowth*ph2)       !probably this 
    _SET_SURFACE_EXCHANGE_(self%id_no3,-fskelno3-(no3-1000./density*iceno3)*(913./density*botmelt+1000./density*abs(fpondno3/iceno3))+913./density*botgrowth*(no3-iceno3))
    _SET_SURFACE_EXCHANGE_(self%id_nh4,-fskelnh4-(nh4-1000./density*icenh4)*(913./density*botmelt+1000./density*abs(fpondnh4/icenh4))+913./density*botgrowth*(nh4-icenh4))
    _SET_SURFACE_EXCHANGE_(self%id_sil,-fskelsil-(sil-1000./density*icesil)*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*(sil-icesil))
    _SET_SURFACE_EXCHANGE_(self%id_ph1,-ph1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*ph1)
    _SET_SURFACE_EXCHANGE_(self%id_de1,-de1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*de1)
    _SET_SURFACE_EXCHANGE_(self%id_zo1,-zo1*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*zo1)
    _SET_SURFACE_EXCHANGE_(self%id_zo2,-zo2*(913./density*botmelt+1000./density*abs(fpondsil/icesil))+913./density*botgrowth*zo2)
#endif

    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialde2,fialde2*self%spd)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialbsi,fialbsi*self%spd)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fialph2,fialph2*self%spd)

   endif
   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
end module uvic_eco
