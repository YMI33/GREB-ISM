!
!--------------------------------------------------------------------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model - Ice Sheet Model (ISM)
!--------------------------------------------------------------------------------------------------------------------
!
!   Authors: Dietmar Dommenget, Janine Flöter, Tobias Bayr, Christian Stassen and Zhiang Xie
!
!   References:	- Conceptual Understanding of Climate Change with a
!              	Globally Resolved Energy Balance Model
!               by Dietmar Dommenget and Janine Flöter, J. Clim Dyn (2011) 37: 2143.
!               doi:10.1007/s00382-011-1026-0

!               - A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0
!               by Christian Stassen, Dietmar Dommenget & Nicholas Loveday.
!               Geosci. Model Dev., 12, 425-440, https://doi.org/10.5194/gmd-12-425-2019, 2019.
!
!               - The Monash Simple Climate Model experiments (MSCM-DB v1.0): 
!               An interactive database of mean climate, climate change, and scenario simulations,
!               by Dommenget, Di., Nice, K., Bayr, T., Kasang, Di., Stassen, C. and Rezny, M., 
!               Geosci. Model Dev., 12(6), 2155–2179, doi:10.5194/gmd-12-2155-2019, 2019.
!
!	           - A coupled ice sheet model for the Global Resolved Energy Balance (GREB-ISM) model 
!               for global simulations on time-scales of 100 kyrs
!	            of the mean climate, climate change and scenarios simulations
!	            by Zhiang Xie, Dietmar Dommenget, Felicity S. McCormack, Andrew N. Mackintosh 
!	            submitted to Geoscientific Model Development
!
!
!  input fields: The GREB model needs the following fields to be specified before
!                the main subroutine greb_model is called:
!
!    Tclim(xdim,ydim,nstep_yr):   mean Tsurf                       [K]
!    uclim(xdim,ydim,nstep_yr):   mean zonal wind speed            [m/s]
!    vclim(xdim,ydim,nstep_yr):   mean meridional wind speed       [m/s]
!    qclim(xdim,ydim,nstep_yr):   mean atmospheric humidity        [kg/kg]
!  cldclim(xdim,ydim,nstep_yr):   total cloud cover                [0-1]
! swetclim(xdim,ydim,nstep_yr):   soil wetness, fraction of total  [0-1]
!   Toclim(xdim,ydim,nstep_yr):   mean deep ocean temperature      [K]
!  mldclim(xdim,ydim,nstep_yr):   mean ocean mixed layer depth     [m]
!precipclim(xdim,ydim,nstep_yr):  mean precipitation               [mm/dy] 
!   b_rock(xdim,ydim):            bed rock data                    [m]
!iceH_clim(xdim,ydim,nstep_yr):   reference ice thickness          [m] 
! sw_solar(ydim,nstep_yr):        24hrs mean solar radiation       [W/m^2]
!
!
! possible experiments:
!
!  log_exp 200s/300s: paleo ice sheet experiments
!
!  log_exp = 200 pi-control CO2= CO2_SCN
!  log_exp = 205 SW-osc-20k 340ppm
!  log_exp = 206 SW-osc-50k 340ppm
!  log_exp = 207 SW-osc-100k 340ppm
!
!  log_exp = 301 paper exp:   10m icesheet
!  log_exp = 302 paper exp: 1000m topo

!  log_exp = 310 experiment: equalibrium experiment CO2= CO2_SCN (standalone)
!  log_exp = 311 experiment: transition experiment forced by CO2 and orbital forcing [250kyr] (standalone)

!  log_exp = 330 experiment: transition experiment forced by CO2 and orbital forcing [250kyr]
!  log_exp = 331 experiment: Eemian scenario (CO2 267.1646, YEAR -127 000) 
!  log_exp = 332 experiment: LGM scenario (CO2 189.6219, YEAR -24 000) 
!
!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: dt        = 12*3600           ! time step [s]
  integer, parameter :: dt_crcl   = 0.5*3600          ! time step circulation [s]

  integer, parameter :: ndt_days  = 24*3600/dt        ! number of timesteps per day
  integer, parameter :: nstep_yr  = ndays_yr*ndt_days ! number of timesteps per year
  integer            :: time_flux = 0                 ! length of integration for flux correction [yrs]
  integer            :: time_ctrl = 0                 ! length of integration for control run  [yrs]
  integer            :: time_scnr = 0                 ! length of integration for scenario run [yrs]
  integer            :: time_out  = 1                 ! length of output time step for scenario run [yrs]
  integer            :: ipx       = 1                 ! points for diagonstic print outs
  integer            :: ipy       = 1                 ! points for diagonstic print outs
  integer, parameter, dimension(12) :: jday_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat
  real               :: dt_ice    = ndays_yr*86400    ! time step for ice sheet model [s]
  integer            :: iprec     = 8                 ! real precision 

  integer            :: ireal     = 4                 ! record length for IO (machine dependent)
! 						        ireal = 4 for Mac Book Pro and Ubuntu Linux, 1 is another option
  integer            :: log_restart = 0               ! control restarts
  integer            :: dt_rstrt   = 10               ! restart frequency [yrs calculated]
  integer            :: dt_out     = 1                ! output frequency  [yrs calculated]

  namelist / numerics / time_flux, time_ctrl, time_scnr, log_restart, dt_rstrt, dt_out

end module mo_numerics


!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  integer  :: log_exp   = 0                    ! process control logics for experiments (see header)
! switches that are the same for both deconstructions
  integer  :: log_hdif        = 1              ! process control Diffusion of heat
  integer  :: log_hadv        = 1              ! process control Advection of heat
  integer  :: log_vdif        = 1              ! process control Diffusion of vapor
  integer  :: log_vadv        = 1              ! process control Advection of vapor
! switches for the hydrological cycle
  integer  :: log_rain        = 0              ! process control precipitation parameterisation
  integer  :: log_eva         = 0              ! process control evaporation parameterisation
  integer  :: log_conv        = 0              ! process control advection parameterisation
  integer  :: log_clim        = 0              ! process control for reference climatology
! switches for ice sheet model exp.
  integer  :: log_ice_sheet   = 1              ! process control for ice sheet
  integer  :: log_ice_dyn     = 1              ! process control for ice sheet dynamic process
  integer  :: log_ice_ithk    = 1              ! ice thickness initial condition (0=ice free; 1=start from current ice thickness)
  integer  :: log_ice_ant     = 1              ! fix ice thickness on Antarctica (0=fix; 1=move)
  integer  :: log_decouple    = 0              ! thermodynamic coupling switch (1 for decouple, 0 for coupling)
! switches for GREB-ice sheet coupling
  integer  :: log_ice_cpl     = 1              ! process control for GREB coupling to ice sheet
  integer  :: log_ice_topo    = 1              ! process control for ice topograph 
  integer  :: log_ice_slv     = 1              ! process control for sea level effect 

! physical parameter (natural constants)
  parameter( pi        = 3.1416 )
  parameter( sig       = 5.6704e-8 )           ! stefan-boltzmann constant [W/m^2/K^4]
  parameter( R_gas     = 8.314 )               ! universal gas constant [J/mol/K]
  parameter( rho_ocean = 999.1 )               ! density of water at T=15C [kg/m^2]
  parameter( rho_land  = 2600. )               ! density of solid rock [kg/m^2]
  parameter( rho_air   = 1.2 )                 ! density of air at 20C at NN [kg/m^2]
  parameter( rho_ice   = 910.)                 ! density of incompressible polycrystalline ice [kg/m^3]
  parameter( cp_ocean  = 4186. )               ! specific heat capacity of water at T=15C [J/kg/K]
  parameter( cp_land   = cp_ocean/4.5 )        ! specific heat capacity of dry land [J/kg/K]
  parameter( cp_air    = 1005. )               ! specific heat capacity of air      [J/kg/K]
  parameter( cp_ice    = 2009. )               ! specific heat capacity of ice [J/kg/K]
  parameter( eps       = 1. )                  ! emissivity for IR
  parameter( grav      = 9.81  )               ! gravitational acceleration [m/s^2]
  real :: S0_var       = 100.                  ! variation of solar constant   [%]

! physical parameter (model values)
  parameter( d_ocean   = 50. )                     ! depth of ocean column [m]
  parameter( d_land    = 2. )                      ! depth of land column  [m]
  parameter( d_air     = 5000. )                   ! depth of air column   [m]
  parameter( d_ice_max = 1. )                      ! maximum shortwave radiation penetarting depth of ice column   [m]
  parameter( d_ice_mov = 10. )                     ! maximum shortwave radiation penetarting depth of ice column   [m]
  parameter( cap_ocean = cp_ocean*rho_ocean )      ! heat capacity 1m ocean  [J/K/m^2]
  parameter( cap_land  = cp_land*rho_land*d_land ) ! heat capacity land   [J/K/m^2]
  parameter( cap_air   = cp_air*rho_air*d_air )    ! heat capacity air    [J/K/m^2]
  parameter( ce        = 2e-3  )                   ! laten heat transfer coefficient for ocean
  parameter( cq_latent = 2.257e6 )                 ! latent heat of condensation/evapoartion f water [J/kg]
  parameter( ci_latent = 3.335e5 )                 ! latent heat of condensation/fusion f ice [J/kg]
  parameter( cq_rain   = -0.1/24./3600. )          ! decrease in air water vapor due to rain [1/s]
  parameter( z_air     = 8400. )                   ! scaling height atmos. heat, CO2
  parameter( z_vapor   = 5000. )                   ! scaling height atmos. water vapor diffusion
  real :: ct_sens   = 22.5                         ! coupling for sensible heat
  real :: da_ice    = 0.25                         ! albedo diff for ice covered points
  real :: a_no_ice  = 0.1                          ! albedo for non-ice covered points
  real :: a_cloud   = 0.35                         ! albedo for clouds
  real :: Tl_ice1   = 273.15-10.                   ! temperature range of land snow-albedo feedback [K]
  real :: Tl_ice2   = 273.15                       ! temperature range of land snow-albedo feedback [K]
  real :: To_ice1   = 273.15-7.                    ! temperature range of ocean ice-albedo feedback [K]
  real :: To_ice2   = 273.15-1.7                   ! temperature range of ocean ice-albedo feedback [K]
  real :: ice_Tcons = 273.15-10.                   ! critical temperature for constitutive equation [K]
  real :: co_turb   = 5.0                          ! turbolent mixing to deep ocean [W/K/m^2]
  real :: kappa     = 5*8e5                        ! atmos. diffusion coefficient [m^2/s]
  real :: r_qviwv = 2.6736e3                       ! regres. factor between viwv and q_air  [kg/m^3]
  real :: CO2_scn = 0.                             ! CO2 in scenario run
  real :: Todeep0   = 273.15+2.                    ! deep abyssal ocean initial temperature
  ! ice sheet model parameters
  real   :: c_lapse   = -0.006                     ! air temperature lapse rate in atmosphere [K/m]
  real*8 :: ice_svlm  = 1./rho_ice                 ! specific volume of ice/snow [m^3/kg]
  real*8 :: ice_kappa = 2.1                        ! ice sheet diffusion coefficient [W/K/m] 
  real*8 :: beta_melt = 8.7d-4                     ! change of melting point with ice depth [K/m]
  real*8 :: geoh_flux = 4.2d-2                     ! geoheat flux on bottom [W/m^2] 
  real*8 :: enh_fact  = 3.                         ! enhance factor for ice flow [1]
  real*8 :: slid_law  = 6d4/365/86400              ! no basal slide
  real*8 :: ice_shlf  = 7000./365/86400            ! ice shelf velocity maximum [m/s]
  real*8 :: shlf_eta  = 2.e14                      ! ice shelf viscosity [Pa s] 
  real*8 :: ref_lat   = 76.875                     ! reference latitude for polar filter [o]
 
 ! paleo runs with orbital forcing
  integer  :: kry_start = -1000                    ! start year of orbital forcing file 
  integer  :: kry_end   = 0                        ! end year in [kyr] has to be <=0  
  integer  :: dtyr_acc  = 1                        ! acceleration in sw_solar forcing 

  
 ! index for paleoclimate simulation
  real     ::  ind_lgm    = -13.99                  ! LGM temperature index (-21 ka)
  real     ::  ind_pre    = -0.89                   ! today temperature index (0 ka)
  real     ::  dT_g                                 ! Greenland temperature difference in paleoclimate
  real     ::  dT_a                                 ! Antarctica temperature difference in paleoclimate
  ! physical paramter (rainfall)
  real :: c_q, c_rq, c_omega, c_omegastd

! parameter emisivity
  real, dimension(10)   :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
&                                    0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)
  real*8, dimension(2)  :: actene = (/1.39d5,6d4/)       ! constitutive equation coefficients, activate energy (>-10oC and <-10oC, Glen's flow law, unit: J/mol)
  real*8, dimension(2)  :: A_fact = (/1.96d3,3.985d-13/) ! constitutive equation coefficients, softness parameter (>-10oC and <-10oC, Glen's flow law, unit: 1/Pa^3/s)

! declare climate fields
  real*4, dimension(xdim,ydim)          ::  z_topo0      ! surface elevation (input initial data) [m]
  real*4, dimension(xdim,ydim)          ::  b_rock0      ! bed rock height (input initial data) [m]
  real*4, dimension(xdim,ydim)          ::  z_ocean      ! ocean depth [m]
  real*4, dimension(xdim,ydim,nstep_yr) ::  iceH_clim0   ! reference ice thickness climatology (input initial data) [m]
  real*4, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, omegaclim, omegastdclim, wsclim
  real*4, dimension(xdim,ydim,nstep_yr) ::  qclim, mldclim, Toclim, cldclim, precipclim,iceH_clim    
  real*4, dimension(xdim,ydim,nstep_yr) ::  Tclim_lgm, pclim, pclim_lgm 
  real*4, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, precip_correct
  real*4, dimension(xdim,ydim,nstep_yr) ::  swetclim, dTrad
  real*4, dimension(ydim,nstep_yr)      ::  sw_solar, sw_solar_ctrl, sw_solar_scnr
  real, dimension(xdim,ydim)            ::  co2_part      = 1.0
  real, dimension(xdim,ydim)            ::  co2_part_scn  = 1.0

! declare field change with time
  real*4, dimension(xdim,ydim)          :: z_topo          ! surface elevation [m] 
  real*4, dimension(xdim,ydim)          :: mask            ! land-sea mask (1 for land, -1 for ocean)
  real*8, dimension(xdim,ydim)          :: glacier_type    ! grid type from ice sheet model perspective (1 for grounded, 0 for floating, -1 for ocean)
  real*8, dimension(xdim, ydim)         :: b_rock          ! bed rock height [m]
  real, dimension(xdim,ydim)            :: cap_surf        ! surface heat capacity [J/K/m^2] 
  real, dimension(xdim, ydim)           :: wz_air          ! topography integration parameter for air mass 
  real, dimension(xdim, ydim)           :: wz_vapor        ! topography integration parameter for water vapor 
  real                                  :: ssh             ! sea level change compared with reference ice thickness [m]
! declare constant fields
  integer jday, ityr, iyear, it_start

! Mike: declare some program constants
  real, dimension(xdim,ydim,nstep_yr)   :: uclim_m, uclim_p
  real, dimension(xdim,ydim,nstep_yr)   :: vclim_m, vclim_p

  real :: t0, t1, t2

  namelist / physics / log_exp, ct_sens, da_ice, a_no_ice, a_cloud, co_turb, kappa, 	&
&                      p_emi, Tl_ice1, Tl_ice2, To_ice1, To_ice2, r_qviwv,          	&
&                      kry_start, kry_end, dtyr_acc, Todeep0, co2_scn, S0_var, &
&                      log_hdif, log_hadv, log_vdif, log_vadv, &
& 		               log_rain, log_eva, log_conv, log_clim,        &
&                      log_ice_sheet, log_ice_ithk, log_ice_cpl,  & 
&                      log_ice_dyn, log_ice_topo, log_ice_slv, log_ice_ant

end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim

 ! declare diagnostic fields
  real*4, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
&                                          ftmn, fqmn, Tomn, ice_Hmn, ice_Tsmn

! declare output fields
  real*4, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, prmm, evamm, qcrclmm, &
&                                          ice_Hm0, ice_Tsmm, term_massmm, term_hadvmm,    & 
&								  	       term_calvmm, albdmn, xdum
  real*4, dimension(xdim,ydim,4)        :: ice_Tmm, ice_vxm0, ice_vym0
  real*4, dimension(xdim,ydim,12)       :: Tmn_ctrl, Tamn_ctrl, Tomn_ctrl, ice_Hmn_ctrl, ice_Tsmn_ctrl, ice_mask_ctrl
  real*4, dimension(xdim,ydim,12)       :: qmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl 
  real*4, dimension(xdim,ydim,12)       :: term_mass_ctrl, term_hadv_ctrl, term_calv_ctrl 
  real*8, dimension(xdim,ydim)          :: ice_Tsm 

end module mo_diagnostics

!+++++++++++++++++++++++++++++++++++++++
program  greb_main
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics
  use mo_diagnostics

! declare temporary fields
  real, dimension(xdim,ydim)     :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1,       &
&                                 ts_ini, ta_ini, q_ini, to_ini
  real*8, dimension(xdim,ydim)   :: ice_Ts0, ice_Ts1, ice_H1, ice_H0, ice_Hini, & 
&                                 ogrid, ice_Tsini
  real*8, dimension(xdim,ydim,4) :: ice_T1, ice_T0, ice_Tini
  integer year, yrs_calc, year0, dt_out0
 
  real*4, dimension(xdim,ydim)   :: ice_out
  integer                        :: t
  real                           :: CO2 

  call greb_initialization
  
  ! open output files
  open(101,file='./control.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(102,file='./scenario.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(103,file='./scenario.gmean.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal)

  ! no ice sheet simulation
  if (log_ice_sheet ==  0) log_ice_cpl = 0
  if (log_ice_cpl   ==  0) then
     log_ice_topo = 0;
     log_ice_slv  = 0;
  end if

  if(log_ice_ithk==0) iceH_clim = 0.                      ! initial value ice sheet thickness

! initialize fields
  Ts_ini   = Tclim(:,:,nstep_yr)                          ! initial value temp. surf
  Ta_ini   = Ts_ini - c_lapse*z_topo                      ! initial value atm. temp.
  To_ini   = Toclim(:,:,nstep_yr)                         ! initial value temp. surf
  q_ini    = qclim(:,:,nstep_yr)                          ! initial value atmos water vapor
  ice_Tsini= Ts_ini                                       ! initial value ice surface temperaturer
  ice_Hini = iceH_clim(:,:,nstep_yr)                      ! initial value ice sheet thickness
  do it=1,4 
     ice_Tini(:,:,it) = ice_Tsini                         ! initial value ice strata temperaturer
  end do    									          ! initialize fields / ice sheet 
 
  ix =93; iy=45
  print*,'ini: ',ix,iy,b_rock0(ix,iy), z_topo(ix,iy), ice_Hini(ix,iy), mask(ix,iy)
  print*,'ice clim: ',log_ice_ithk, iceH_clim(17,2,nstep_yr),iceH_clim(17,41,nstep_yr),iceH_clim(86,46,nstep_yr)
  print*,'ice height: ',ice_Hini(17,2),ice_Hini(17,41),ice_Hini(86,46)

  CO2_ctrl = 340.0
! decon mean state switch
  if (log_exp .gt. 310 .and. log_exp .lt. 330 )  CO2_ctrl = 280.  ! paleo transition experiments 

  sw_solar = sw_solar_ctrl

! define some program constants
  wz_air   = exp(-z_topo/z_air)
  wz_vapor = exp(-z_topo/z_vapor)
  where (uclim(:,:,:) >= 0.0)
     uclim_m = uclim
     uclim_p = 0.0
  elsewhere
     uclim_m = 0.0
     uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0)
     vclim_m = vclim
     vclim_p = 0.0
  elsewhere
     vclim_m = 0.0
     vclim_p = vclim
  end where
          
  ! initialize the output variables
  Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; prmm=0.
  ice_Tsmm=0.; term_hadvmm=0.; term_massmm=0.; term_calvmm=0.
  ice_Tmm=0.; ice_Hm0=0.; albdmn=0.

  ! initialize the rainfall parameterisation
  select case( log_rain )
    case(-1) ! Original GREB
      c_q=1.; c_rq= 0.; c_omega=0.; c_omegastd=0.
    case(1) ! Adding relative humidity (rq)
      c_q=-1.391649; c_rq=3.018774; c_omega= 0.; c_omegastd=0.
    case(2) ! Adding omega
      c_q=0.862162; c_rq=0.; c_omega=-29.02096; c_omegastd=0.
    case(3) ! Adding rq and omega
      c_q=-0.2685845; c_rq=1.4591853; c_omega=-26.9858807; c_omegastd=0.
    case(0) ! Best GREB
        c_q=-1.88; c_rq=2.25; c_omega=-17.69; c_omegastd=59.07 ! Rainfall parameters (ERA-Interim)
      if (log_clim == 1) then
        c_q=-1.27; c_rq=1.99; c_omega=-16.54; c_omegastd=21.15 ! Rainfall parameters (NCEP)
      end if
  end select

  print*,'log_restart = ',log_restart
   
  dt_out0 = dt_out
   
  if(log_restart == 0) then
    !------------------------------------------------------------
    ! newstart
    !------------------------------------------------------------
     
      dt_out     = 1    ! output frequency  [yrs calculated] - 1 year for ctrl runs.
      year0      = 1
      ! year/1000 + kry_start = how many years before
      call orbital_forcing(-(1000*kry_start+1))    ! sw_solar forcing
   
   !------------------------------------------------------------
   ! compute Q-flux corrections
   !------------------------------------------------------------
      print*,'% flux correction ', CO2_ctrl
      1001 format (A4,     T8, A10,   T20, A10,    T32, A15,         T50, A12,      T65, A12,     T80, A15) !TB
      print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
      call qflux_correction(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, ice_Tsini, ice_Hini, ice_Tini)

    !  write flux corrections
     open(31,file='Tsurf_flux_correction.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
     do irec=1, nstep_yr
        write(31,rec=irec)  TF_correct(:,:,irec)
     end do
     
    !------------------------------------------------------------
    ! control run
    !------------------------------------------------------------
      print*,'% CONTROL RUN CO2=',CO2_ctrl,'  time=', time_ctrl,'yr'
      print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
      Ts1 = Ts_ini; Ta1 = Ta_ini; To1 = To_ini; q1 = q_ini; ice_Ts1 = Ts1; ice_H1 = ice_Hini  
      do it=1,4 
          ice_T1(:,:,it)=Ts1 
      end do    ! initialize fields / ice sheet 
      year=1970; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
      do it=1, time_ctrl*nstep_yr                                             ! main time loop
      
         ityr = mod((it-1),nstep_yr)+1           ! time step in year

          call time_loop(it, isrec, year, CO2_ctrl, irec, mon, 101, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, &
    &                   ice_Ts0, ice_H0, ice_T0, ice_Ts1, ice_H1, ice_T1)
          Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0; ice_H1=ice_H0
          ice_Ts1=ice_Ts0;ice_H1=ice_H0;ice_T1=ice_T0

      end do
      dt_out     = 1000 ! output frequency  [yrs calculated] - 1000 years for ctrl runs.
  end if
  
    !------------------------------------------------------------
    ! scenario run
    !------------------------------------------------------------
  print*,'% PALEO RUN: ',log_exp,'  time=', time_scnr,'kyrs'
  print*,'log_restart = ',log_restart
  dt_out = dt_out0

  print*,'% SCENARIO EXP: ',log_exp,'  time=', time_scnr,'yr'
  print 1001, "YEAR", "CO2[ppm]", "SW[W/m^2]", "global mean[C]", "Trop Pac[C]", "Hamburg[C]", "North Pole[C]" !TB
  Ts1 = Ts_ini; Ta1 = Ta_ini; q1 = q_ini; To1 = To_ini; ice_H1 = ice_Hini; ice_Ts1 = Ts1                     ! initialize field / ice sheet
  do it=1,4 
    ice_T1(:,:,it)=Ts1 
  end do    ! initialize fields / ice sheet 

  year=1.; iyear=1; CO2=CO2_SCN; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  it_start=1;

  yrs_calc = 0

!  check for restart
  if(log_restart == 1) then
     call restart_read(it, year, ts0, ta0, to0, q0, ice_Ts0, ice_H0, ice_T0)
     ityr = mod((it-1),nstep_yr)+1           ! time step in year
     Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
     ice_Ts1=ice_Ts0; ice_H1=ice_H0; ice_T1=ice_T0
  if( ityr == 1 .and. (log_exp .ne. 310) ) &
  &   call sealevel(ice_H0, Ts0, To0)
     print*, 'Do restart', year 
  end if  
  
  ! sw_solar forcing
  call orbital_forcing(year)  
  if (log_exp .eq. 200) call orbital_forcing(-(1000*kry_start+1)) ! EXP: PI-COTNTORL RUN
  if (log_exp .eq. 301) call orbital_forcing(-(1000*kry_start+1)) ! EXP: PI-COTNTORL RUN
  if (log_exp .eq. 302) call orbital_forcing(-(1000*kry_start+1)) ! EXP: PI-COTNTORL RUN
  if( log_exp .eq. 310  ) then 
      sw_solar = sw_solar_ctrl 
      Ts1 = Tclim(:,:,ityr)
  end if
  if( log_exp .eq. 311  ) then 
      sw_solar = sw_solar_ctrl 
      if(mod(year,100) .eq. 1 .and. ityr .eq. 1) read (53,*) t, dT_g
      if(mod(year,100) .eq. 1 .and. ityr .eq. 1) read (54,*) t, dT_a
      Ts1 = Tclim(:,:,ityr) + (Tclim_lgm(:,:,ityr) - Tclim(:,:,ityr))*(dT_g-ind_pre)/(ind_lgm-ind_pre)
  end if
  if (log_exp .eq. 331) call orbital_forcing(1000*(1000-127)+1)
  if (log_exp .eq. 332) call orbital_forcing(1000*(1000-24)+1)

! paper experiment: 10m ice sheet
! X is varying   Lon = 40 to 80   X = 11.6667 to 22.3333
! Y is varying   Lat = 40 to 70   Y = 35.1667 to 43.1667
  if (log_exp .eq. 301) ice_H1(1:12,35:43) = 10.0 ! europe

  ! paper experiment: 1000m topo 
  if (log_exp .eq. 302) b_rock(17:33,35:42) = b_rock(17:33,35:42) +1000.0 ! central asia

  do it=it_start, time_scnr*nstep_yr                                              ! main time loop

    ityr = mod((it-1),nstep_yr)+1           ! time step in year

!    write restart
     call restart_write(it, year, yrs_calc, ts1, ta1, to1, q1, ice_Ts1, ice_H1, ice_T1)

     CO2=CO2_SCN
     if(log_exp == 330 .and. mod(year,100) .eq. 1 .and. ityr .eq. 1) read(53,*) t,CO2_SCN
     if(log_exp == 331) CO2_SCN=267.1646
     if(log_exp == 332) CO2_SCN=189.6219

     call time_loop(it,isrec, year, CO2, irec, mon, 102, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, &
&                   ice_Ts0, ice_H0, ice_T0, ice_Ts1, ice_H1, ice_T1)
     Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0
     ice_Ts1=ice_Ts0; ice_H1=ice_H0; ice_T1=ice_T0
     if (mod(it,nstep_yr) == 0) year     = year     +dtyr_acc
     if (mod(it,nstep_yr) == 0) yrs_calc = yrs_calc +1
     if (mod(it,nstep_yr) == 0) call orbital_forcing(year)
     if (mod(it,nstep_yr) == 0 .and. log_exp == 200) call orbital_forcing(-(1000*kry_start+1))
     if (mod(it,nstep_yr) == 0 .and. log_exp == 301) call orbital_forcing(-(1000*kry_start+1))
     if (mod(it,nstep_yr) == 0 .and. log_exp == 302) call orbital_forcing(-(1000*kry_start+1))
     if( log_exp .eq. 310  ) then 
         sw_solar = sw_solar_ctrl 
         Ts1 = Tclim(:,:,ityr)
     end if
     if( log_exp .eq. 311  ) then 
         sw_solar = sw_solar_ctrl 
         if(mod(year,100) .eq. 1 .and. ityr .eq. 1) read (53,*) t, dT_g
         if(mod(year,100) .eq. 1 .and. ityr .eq. 1) read (54,*) t, dT_a
         Ts1 = Tclim(:,:,ityr) + (Tclim_lgm(:,:,ityr) - Tclim(:,:,ityr))*(dT_g-ind_pre)/(ind_lgm-ind_pre)
     end if
     if (mod(it,nstep_yr) == 0 .and. log_exp == 331) call orbital_forcing((1000-127)*1000+1)
     if (mod(it,nstep_yr) == 0 .and. log_exp == 332) call orbital_forcing((1000-24)*1000+1)

! TEST sensitivity experiment: SW oscillation
    if (log_exp == 205) period =  20000
    if (log_exp == 206) period =  50000
    if (log_exp == 207) period = 100000
    if (log_exp == 208) period = 10000
    if (log_exp >= 205 .and. log_exp <= 208 ) then
       sw_solar = sw_solar_ctrl 
       swmag0 =0.1
       swmag = swmag0; SW_solar(33:48,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(33:48,:)
       swmag = 0.75*swmag0; SW_solar(32,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(32,:)
       swmag = 0.60*swmag0; SW_solar(31,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(31,:)
       swmag = 0.50*swmag0; SW_solar(30,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(30,:)
       swmag = 0.25*swmag0; SW_solar(29,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(29,:)
       swmag = 0.10*swmag0; SW_solar(28,:) = (1+swmag*(sin(2*3.14*year/period)))*sw_solar_ctrl(28,:)
    end if
    
  end do

  call restart_write(it, year, yrs_calc, ts1, ta1, to1, q1, ice_Ts1, ice_H1, ice_T1)

end program greb_main

!+++++++++++++++++++++++++++++++++++++++
subroutine greb_initialization
!+++++++++++++++++++++++++++++++++++++++
  ! initialisation of the greb model

  USE mo_numerics
  USE mo_physics
  use mo_diagnostics, only: xdum

  ! declare output fields
  real, dimension(xdim,ydim,ndays_yr) :: Tc1, Ta1, q1, ap1
  real, dimension(xdim,ydim,ndays_yr) :: Tc2, Ta2, q2, ap2

  integer, dimension(ndays_yr)::  t = (/(i,i=1,ndays_yr)/) ! jday index

  100 FORMAT('climate: ',F9.2, 5E12.4)

print*,'% start climate shell'

!+++++++++++++++++++++++++++++ read input file  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! open input files
  open(10,file='namelist')

  ! read namelist
  read(10,numerics)
  read(10,physics)

  ! read climatologies
  open(11,file='tclim',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(12,file='uwind',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(13,file='vwind',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(14,file='qclim',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(15,file='cloud',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(16,file='moist',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(17,file='toclim',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(18,file='mldclim',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(19,file='ztopo',  	  ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(21,file='solar',       ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*ydim*nstep_yr)
  open(22,file='abswind',     ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(23,file='omclim',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(24,file='omstdv',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(25,file='precip',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

  open(51,file='bedrock',     ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(52,file='iceclm',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

 ! read climatologies
 do n=1,nstep_yr
    read(11,rec=n) tclim(:,:,n)
    read(12,rec=n) uclim(:,:,n)
    read(13,rec=n) vclim(:,:,n)
    read(14,rec=n) qclim(:,:,n)
    read(15,rec=n) cldclim(:,:,n)
    read(16,rec=n) swetclim(:,:,n)
    read(17,rec=n) Toclim(:,:,n)
    read(18,rec=n) mldclim(:,:,n)
    read(22,rec=n) wsclim(:,:,n)
    read(23,rec=n) omegaclim(:,:,n)
    read(24,rec=n) omegastdclim(:,:,n)
    read(25,rec=n) precipclim(:,:,n)
    read(52,rec=n) iceH_clim(:,:,n)
  end do
  
   ! read fix data
   read(19,rec=1)  z_topo0
   read(21,rec=1)  sw_solar_ctrl
   read(51,rec=1)  b_rock0

!+++++++++++++++++++++++++++++ parameter initialization  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! define land mask
   z_topo = z_topo0
   where ( z_topo < 0 ) z_topo = -0.1  ! ocean points
   mask = 1
   where ( z_topo < 0 ) mask = -1  ! ocean points
   glacier_type = mask

   ! fix mismatch between Tsurf and iceH_clim (glaciers)
   do n=1,nstep_yr
      where ( iceH_clim (:,:,n) > 1.0 .and.  tclim(:,:,n)-273.15 >= -1.0 )               &
&           tclim(:,:,n) = 273.15 -1.0  ! ice = freezing
   end do

   ! fix Tsurf at n-pole < -1.7 to have sea ice
   where( tclim(:,47:48,:)-273.15 > -1.7 ) tclim(:,47:48,:) = 273.15-1.7

   ! fix bed rock over ocean points
   where( mask < 0 ) b_rock0 = z_topo0
  
   ! forcing data
   if(log_exp == 330 ) open(53,file='co2forcing') 
   if(  log_exp .eq. 311) then
      open(53,file='delta_ts_Greenland')
      open(54,file='delta_ts_Antarctica')
      open(55,file='lgm_tsurf' , ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
      open(56,file='lgm_precip', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

      do n=1,nstep_yr
          read(55,rec=n) Tclim_lgm(:,:,n)
          read(56,rec=n) pclim_lgm(:,:,n)
          pclim(:,:,n) = precipclim(:,:,n)/86400.
      end do
   end if

  ! start greb_model run
  print*,'ice ini: ',iceH_clim(86,46,1),iceH_clim(86,46,360),b_rock(86,46), z_topo(86,46), mask(86,46)
  print*,'% time flux/control/scenario: ', time_flux, time_ctrl, time_scnr

  ! set ocean depth
  z_ocean=0
  do i=1,nstep_yr
     where(mldclim(:,:,i).gt.z_ocean) z_ocean = mldclim(:,:,i)
  end do
  z_ocean = 3.0*z_ocean

  dTrad = -0.16*Tclim -5. ! offset Tatmos-rad
  
  ! ice sheet set topography (marine ice scheme)
  b_rock = b_rock0
  z_topo = - 0.1
  ! grounded ice 
  where(b_rock >= -iceH_clim(:,:,nstep_yr)*rho_ice/rho_ocean) z_topo = b_rock + iceH_clim(:,:,nstep_yr)
  ! floating ice 
  where(iceH_clim(:,:,nstep_yr)>0. .and. b_rock <  -iceH_clim(:,:,nstep_yr)*rho_ice/rho_ocean) &
  z_topo = iceH_clim(:,:,nstep_yr)*(1.-rho_ice/rho_ocean)
 
! heat capacity global [J/K/m^2]
  where (mask  > 0.) cap_surf = cap_land
  where (mask  < 0.) cap_surf = cap_ocean*mldclim(:,:,1)
end subroutine greb_initialization

!+++++++++++++++++++++++++++++++++++++++
subroutine time_loop(it, isrec, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, &
&                   ice_Ts0, ice_H0, ice_T0, ice_Ts1, ice_H1, ice_T1)
!+++++++++++++++++++++++++++++++++++++++
! main time loop

  use mo_numerics
  use mo_physics

  integer year
  real, dimension(xdim,ydim)   :: Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, sw,       &
&                                 ice_cover, Q_sens, Q_lat, Q_lat_air, dq_eva,   &
&                                 dq_rain, dTa_crcl, dq_crcl, dq, dT_ocean, dTo, &
&                                 LW_surf, LWair_down, LWair_up, em,             &
&                                 a_surf, Q_ice, Q_sice, pmon            
  real*8, dimension(xdim,ydim) :: ice_H0,  ice_Ts0, Fn_surf, ice_H1, ice_Ts1,    & 
&                                 dice_hadv, dice_mass, term_calv, dice_h
  real*8, dimension(xdim,ydim,4) :: ice_T1, ice_T0
  real*8, dimension(xdim,ydim,4) :: ice_vx3    ! ice strata velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim,4) :: ice_vy3    ! ice strata velocity at zonal (x) and meridian (y) direction [m/s]

  937 format ('tendencies: ',I5, 10F10.2, E15.2) !TB

  jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
  ityr = mod((it-1),nstep_yr)+1           ! time step in year

  if(log_exp == 310 .or. log_exp == 311) then
      Ta1  = Ts1 - c_lapse*z_topo
      To1  = Toclim(:,:,ityr)
      q1   = qclim(:,:,ityr)
  endif
  call tendencies(CO2, Ts1, Ta1, To1, q1, ice_H1, SW, LW_surf, Q_lat,         &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl,             &
&                    dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em, a_surf, Q_sice,dice_h)

   if(log_exp == 310) then
       dq_rain = -precipclim(:,:,ityr)/(3600*24*r_qviwv*wz_vapor) ! unit: kg/s
   elseif(log_exp == 311) then
       pmon    = pclim(:,:,ityr) + (pclim_lgm(:,:,ityr)-pclim(:,:,ityr))*(dT_g-ind_pre)/(ind_lgm-ind_pre)
       where(pmon < 0) pmon = 0.
       dq_rain = -pmon/(wz_vapor*r_qviwv) ! unit: kg/s
   endif

  ! ice sheet : spread and ablation
  if(log_ice_sheet == 1) then
      ice_Ts1  = Ts1
      Fn_surf  = sw + LW_surf - LWair_down + Q_sens + Q_lat + TF_correct(:,:,ityr)/2.
      call ice_sheet(it, ionum, irec, mon, ice_H1, ice_T1, ice_Ts1, Ta1, dT_ocean, Fn_surf,&
                     dq_rain, ice_H0, ice_T0, ice_Ts0, &
                     dice_mass, dice_hadv, term_calv, Q_ice, year)

   ! no dice_mass, dice_hadv for ocean & Tsurf > To_ice2 & ice_H1 == 0
      where(mask < 0 .and.  Ts1 > To_ice2) dice_hadv = 0.
      where(mask < 0 .and.  ice_h1 == 0)   dice_hadv = 0.
      if(log_ice_ant == 0 ) then
          Q_ice(:,1:20)     = 0.
      end if
  else
     dice_mass = 0.
     dice_hadv = 0.
     Q_ice     = 0.
     ice_H1    = iceH_clim(:,:,ityr)
  end if

  if(log_ice_cpl == 0) Q_ice = 0.0

  ! equation (28) in XIE2021, surface temperature, equation (1) for old version
  Ts0  = Ts1  +dT_ocean + dt*( SW +LW_surf -LWair_down +Q_lat +Q_ice +Q_sice +Q_sens   &
&                              +TF_correct(:,:,ityr)) / cap_surf

  ! sea ice conditions: Tsurf <= To_ice2 if sea ice present
  ! where(mask < 0 .and. ice_H1 > 0.01 .and. Ts0 > To_ice2)  Ts0 = To_ice2

  Tmin_limit = 40 ! no very low Tsurf/Tatmoss;  numerical stability
  where(Ts0 .le. Tmin_limit )     Ts0 = Tmin_limit ! no very low Tsurf;  numerical stability

  ! air temperature
  ! equation (2) in XIE2021, air temperature
  Ta0  = Ta1 +dTa_crcl +dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens )/cap_air
  where(Ta0 .le. Tmin_limit )     Ta0 = Tmin_limit ! no very low Tatmos;  numerical stability

  ! deep ocean temperature
  ! equation (3) in XIE2021, ocean temperature
  where(mask < 0) To0  = To1 +dTo +ToF_correct(:,:,ityr)

  ! air water vapor
  ! equation (4) in XIE2021, surface humidity
  dq = dt*(dq_eva+dq_rain) +dq_crcl + qF_correct(:,:,ityr)
  where(dq .le. -q1 ) dq = -0.9*q1     ! no negative q;  numerical stability
  where(dq .gt.   0.020 ) dq = 0.020   ! no hugh q;  numerical stability

  q0 = q1 + dq

!  ice height
!  equation (7)/(31) in XIE2021
  ice_H0 = ice_H1 +dice_h +dice_mass +dice_hadv
  where (ice_H0 < 0.) ice_H0 = 0.0
  do k = 1,4
      where(ice_H0 <= d_ice_mov .or. (ice_H1 <= d_ice_mov .and. ice_H0 > d_ice_mov)) ice_T0(:,:,k) = Ts0 
  end do
  if(log_ice_ant == 0 ) then
      ice_H0(:,1:20)    = iceH_clim(:,1:20,ityr)
  end if

  if( ityr == 1 .and. (log_exp .ne. 310) ) &
  &   call sealevel(ice_H0, Ts0, To0)
    
  ! write output
  call output(it, year, ionum, irec, mon, ts0, ta0, to0, q0, dq_rain, &
                  ice_Ts0, ice_H0, dice_mass, dice_hadv, term_calv, ice_T0, &
                  ice_vx3, ice_vy3, a_surf)
  ! diagnostics: annual means plots
  call diagonstics(it, year, CO2, ts0, ta0, to0, q0, ice_H0, SW, lw_surf, q_lat, q_sens)

end subroutine time_loop

!+++++++++++++++++++++++++++++++++++++++
subroutine tendencies(CO2, Ts1, Ta1, To1, q1, ice_h1, SW, LW_surf, Q_lat, Q_sens, Q_lat_air,  &
&                     dq_eva, dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, & 
&                     LWair_down, LWair_up, em, a_surf, Q_sice, dice_h)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts1, Ta1, To1, q1, sw, LWair_up,    &
&                                LWair_down, em, Q_sens, Q_lat, Q_lat_air,      &
&                                dq_eva, dq_rain, dTa_crcl, dq_crcl, LW_surf,   &
&                                dT_ocean, dTo, a_surf, Ta_scl, Q_sice
  real*8,dimension(xdim,ydim) :: ice_H1,dice_h

    ! topograph parameters
    wz_air   = exp(-z_topo/z_air)
    wz_vapor = exp(-z_topo/z_vapor)
!$omp parallel sections
!$omp section
    ! SW radiation model
    call SWradiation(Ts1, sw, ice_h1, a_surf)
    ! LW radiation model
    call LWradiation(Ts1, Ta1, q1, CO2, LW_surf, LWair_up, LWair_down, em)
    ! sensible heat flux
    call senseheat(Ta1,Ts1,Q_sens) 

    ! hydro. model
    call hydro(Ts1, q1, Q_lat, Q_lat_air, dq_eva, dq_rain, ice_H1)
    
    ! deep ocean interaction
    call deep_ocean( Ts1, To1, dT_ocean, dTo)

    ! sea ice
    call seaice(Ts1, ice_H1, SW, LW_surf, LWair_down, Q_lat, &
                  Q_sens, dT_ocean, Q_sice, dice_h)
    
!$omp section
    ! atmos. circulation
    call circulation( Ta1, dTa_crcl, z_air, wz_air)       ! air temp
!$omp section
    call circulation( q1,  dq_crcl, z_vapor, wz_vapor)   ! atmos water vapor
!$omp end parallel sections

end subroutine tendencies

!+++++++++++++++++++++++++++++++++++++++
 subroutine  qflux_correction(CO2_ctrl, Ts1, Ta1, q1, To1, ice_Ts1, ice_H1, ice_T1)
!+++++++++++++++++++++++++++++++++++++++
!  compute flux correction values

  USE mo_numerics
  USE mo_physics

! declare temporary fields
  real*4, external             :: gmean
  real, dimension(xdim,ydim)   :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1, sw, ice_cover, &
&                               Q_sens, Q_lat, Q_lat_air, dq_eva, dq_rain, LW_surf,  &
&                               LWair_down, LWair_up, em, dTa_crcl, dq_crcl, dTs,    &
&                               dTa, dq, T_error, dT_ocean, dTo, &
&                               Q_ice, a_surf,q_zonal, Q_sice 
  real*8,dimension(xdim,ydim)  :: Fn_surf, ice_H1, ice_Ts1, dice_mass, dice_hadv,      &
&								term_calv, ice_vxm, ice_vym, ice_H0, ice_Ts0,   &
&								dice_h, xx

  real*8,dimension(xdim,ydim,4) :: ice_T1, ice_T0
  integer year
  
  1019 format ('qflux: ',I5, 6F8.2, F9.3, F8.3) !TB

  year = 0

! time loop
  do it=1, time_flux*ndt_days*ndays_yr
     
     jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
     ityr = mod((it-1),nstep_yr)+1           ! time step in year

     if(log_exp == 310 .or. log_exp == 311) then
         Ta1  = Ts1 - c_lapse*z_topo
         To1  = Toclim(:,:,ityr)
         q1   = qclim(:,:,ityr)
     endif

     call tendencies(CO2_ctrl, Ts1, Ta1, To1, q1, ice_H1, SW, LW_surf, Q_lat,        &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl, dTa_crcl,          &
&                    dT_ocean, dTo, LWair_down, LWair_up, em, a_surf, Q_sice,dice_h)

     ! ice sheet : spread and ablation
     if(log_ice_sheet == 1) then
        ice_Ts1  = Ts1
        Fn_surf  = sw + LW_surf - LWair_down + Q_sens + Q_lat + TF_correct(:,:,ityr)/2.
        call ice_sheet(it, ionum, irec, mon, ice_H1, ice_T1, ice_Ts1, Ta1, dT_ocean    &
&                    , Fn_surf, dq_rain, ice_H0, ice_T0, ice_Ts0             &
&                    , dice_mass, dice_hadv, term_calv, Q_ice, year)


       !! DD FIX: no dice_mass, dice_hadv for ocean & Tsurf > To_ice2 & ice_H1 == 0
        where(mask < 0 .and.  Ts1 > To_ice2) dice_hadv = 0.
        where(mask < 0 .and.  ice_h1 < 0.01) dice_hadv = 0.
        if(log_ice_ant == 0 ) then
            Q_ice(:,1:20)     = 0.
        end if
     else
        dice_mass = 0.
        dice_hadv = 0.
        Q_ice     = 0.
        ice_H1    = iceH_clim(:,:,ityr)
     end if

    if(log_ice_cpl == 0) Q_ice = 0.0

    ! equation (28) in XIE2021, surface temperature
    ! surface temperature without heat flux correction
    dTs = dt*( sw +LW_surf -LWair_down +Q_lat +Q_sens +Q_ice +Q_sice) / cap_surf
    Ts0  = Ts1 +dTs +dT_ocean
    ! air temperature
    dTa = dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens)/cap_air
    Ta0  = Ta1 + dTa +dTa_crcl
    ! deep ocean temperature without heat flux correction
    To0  = To1 +dTo
    ! air water vapor without flux correction
    dq = dt*(dq_eva+dq_rain)
    q0 = q1 +dq +dq_crcl

    ! heat flux correction Tsurf
    T_error              = Tclim(:,:,ityr) -Ts0 ! error relative to Tclim
    TF_correct(:,:,ityr) = T_error*cap_surf/dt ! heat flux in [W/m^2]

    Ts0  = Ts1 +dTs +dT_ocean +TF_correct(:,:,ityr)*dt/ cap_surf

    ! sea ice conditions: Tsurf <= To_ice2 if sea ice present
     !where(mask < 0 .and. ice_H1 > 0.01 .and. Ts0 > To_ice2)  Ts0 = To_ice2
        
    ! heat flux correction deep ocean
    where(mask < 0) ToF_correct(:,:,ityr) = Toclim(:,:,ityr) -To0  ! heat flux in [K/dt]
    To0  = To1 +dTo +ToF_correct(:,:,ityr)

    ! precip flux corrections (land only)
    do i=1,ydim
      q_zonal(:,i) = sum(q1(:,i))/xdim
    end do
    do i=1,7; q_zonal(:,i) = sum(q_zonal(1,i:8))/(8-i+1); end do ! smooth antarctica to s-ocean
    do i=ydim-6,ydim; q_zonal(:,i) = sum(q_zonal(1,ydim-7:i))/(8-(ydim-i)); end do ! smooth n-pole
    dq_rain = dq_rain-q_zonal*precip_correct(:,:,ityr)
    precip_correct(:,:,ityr) = -precipclim(:,:,ityr)/(3600*24*r_qviwv*wz_vapor) - dq_rain
    precip_correct(:,:,ityr) = precip_correct(:,:,ityr)/q_zonal
    where(mask < 0) precip_correct(:,:,ityr) = 0.0

    ! water vapor flux correction
    qF_correct(:,:,ityr) = qclim(:,:,ityr) -q0 
    
    ! air water vapor with flux correction
    q0 = q1 + dq +dq_crcl + qF_correct(:,:,ityr)
 
    ! ice height
    !  equation (7)/(30) in XIE2021
    ice_h0 = ice_h1 +dice_h +dice_mass +dice_hadv
    where ( ice_h0 < 0. ) ice_h0 = 0.0
    do k = 1,4
        where(ice_H0 <= d_ice_mov .or. (ice_H1 <= d_ice_mov .and. ice_H0 > d_ice_mov)) ice_T0(:,:,k) = Ts0 
    end do
    if(log_ice_ant == 0 ) then
        ice_H0(:,1:20)    = iceH_clim(:,1:20,ityr)
    end if

    if( ityr == 1 .and. (log_exp .ne. 310) ) &
    &   call sealevel(ice_H0, Ts0, To0)

    ! diagnostics: annual means plots
    call diagonstics(it, 0, CO2_ctrl, ts0, ta0, to0, q0, ice_H0, sw, lw_surf, q_lat, q_sens)

    ! memory
    Ts1=Ts0; Ta1=Ta0; q1=q0;  To1=To0;
    ice_H1=ice_H0; ice_T1=ice_T0; ice_Ts1=ice_Ts0

  end do
  1003 format ("On global average a heat flux correction of ", F8.2," W/m^2") !TB
  1004 format ("and a water vapour correction of ", F8.4, " g/kg is applied each time step") !TB
  print 1003, gmean(sum(TF_correct,3)/nstep_yr) !TB
  print 1004, gmean(sum(qF_correct,3)/nstep_yr)*100 !TB
end subroutine qflux_correction

!+++++++++++++++++++++++++++++++++++++++
subroutine diagonstics(it, year, CO2, ts0, ta0, to0, q0, ice_cover, sw, lw_surf, q_lat, q_sens)
!+++++++++++++++++++++++++++++++++++++++
!    diagonstics plots

  USE mo_numerics,    ONLY: ndays_yr, xdim, ydim, ipx ,ipy, ndt_days, nstep_yr
  USE mo_physics,     ONLY: ityr, log_exp, TF_correct, qF_correct, cap_surf, Tclim
  use mo_diagnostics

! declare temporary fields 
  real*4, external            :: gmean
  real*8, dimension(xdim,ydim):: ice_cover
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, sw, Q_sens, Q_lat,  LW_surf
  integer year

  ! diagnostics: annual means
  tsmn=tsmn+Ts0; tamn=tamn+ta0; tomn=tomn+to0; qmn=qmn+q0; 
  swmn=swmn+sw;  lwmn=lwmn+LW_surf; qlatmn=qlatmn+q_lat; qsensmn=qsensmn+Q_sens;
  ftmn=ftmn+TF_correct(:,:,ityr); fqmn=fqmn+qF_correct(:,:,ityr);
  if ( ityr == nstep_yr ) then
     tsmn    = tsmn/nstep_yr;      tamn = tamn/nstep_yr;    tomn = tomn/nstep_yr;
     qmn     = qmn/nstep_yr;
     swmn = swmn/nstep_yr;    lwmn = lwmn/nstep_yr;
     qlatmn  = qlatmn/nstep_yr; qsensmn = qsensmn/nstep_yr; ftmn = ftmn/nstep_yr;
     fqmn    = fqmn/nstep_yr;
     1000 format (I7.0, T8, F10.1, T20, F10.2, T37, F10.6, T52, F10.6, T67, F10.6, T85, F10.6) !TB
     1001 format (I7.0, T8, F10.1, T20, F10.3, T37, F10.3, T52, F10.3, T67, F10.3, T85, F10.3) !TB
     if(log_exp .ne. 310 .and. log_exp .ne. 311) then
         print 1000, year, CO2, gmean(swmn), gmean(tsmn)-273.15, tsmn(48,24)-273.15,tsmn(4,39)-273.15, tsmn(1,48)-273.15 !TB
     else
         print 1001, year, Ts0(65,43),ice_cover(65,43),ice_cover(33,43),ice_cover(25,33),ice_cover(86,43),ice_cover(33,1) !TB
     endif

     tsmn=0.; tamn=0.; qmn=0.; swmn=0.;        ! reset annual mean values
     lwmn=0.; qlatmn=0.; qsensmn=0.; ftmn=0.; fqmn=0.;  ! reset annual mean values
  end if

end subroutine diagonstics

!+++++++++++++++++++++++++++++++++++++++
subroutine restart_write(it, year, yrs_calc, ts, ta, to, q, ice_Ts, ice_H, ice_T)
!+++++++++++++++++++++++++++++++++++++++
!    write restart file

  USE mo_numerics,     ONLY: xdim, ydim, nstep_yr, dt_rstrt, ireal, iprec
  USE mo_physics,      ONLY: sw_solar, kry_start, TF_correct, qF_correct, ToF_correct &
&						   , precip_correct, iceH_clim0

  integer year, yrs_calc

  real, dimension(xdim,ydim)::   Ts, Ta, q, To, ice_Ts4, ice_H4
  real, dimension(xdim,ydim,4):: ice_T4
  real*8, dimension(xdim,ydim)::   ice_Ts, ice_H
  real*8, dimension(xdim,ydim,4):: ice_T

  if (mod(it-1,nstep_yr) == 0 .and. mod(yrs_calc,dt_rstrt) == 0) then
     xyr  = kry_start + (year-1)/1000.
     print *,'write restart',it,xyr,yrs_calc,dt_rstrt
     open(31,file='restart.txt')
     write (31,*) it, year
     close(31)

     ice_Ts4 = ice_Ts; ice_H4 = ice_H; ice_T4 = ice_T
     open(31,file='restart.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
     write(31,rec=1) Ts
     write(31,rec=2) Ta
     write(31,rec=3) q
     write(31,rec=4) To
     write(31,rec=5) ice_Ts4
     write(31,rec=6) ice_H4
     write(31,rec=7) ice_T4(:,:,1)
     write(31,rec=8) ice_T4(:,:,2)
     write(31,rec=9) ice_T4(:,:,3)
     write(31,rec=10) ice_T4(:,:,4)     
   	 do irec=1, nstep_yr
       	write(31,rec=4*irec+8)   TF_correct(:,:,irec)
       	write(31,rec=4*irec+9)   qF_correct(:,:,irec)
       	write(31,rec=4*irec+10)  ToF_correct(:,:,irec)
       	write(31,rec=4*irec+11)  precip_correct(:,:,irec)
   	 end do
     close(31)

  end if
     
end subroutine restart_write

!+++++++++++++++++++++++++++++++++++++++
subroutine restart_read(it, year, ts, ta, to, q, ice_Ts, ice_H, ice_T)
!+++++++++++++++++++++++++++++++++++++++
!    read restart file

  USE mo_numerics,     ONLY: xdim, ydim, nstep_yr, dt_rstrt, ireal, iprec
  USE mo_physics,      ONLY: log_exp, sw_solar, kry_start, TF_correct, qF_correct, ToF_correct &
&                          , precip_correct, CO2_SCN  

  integer year, t

  real,   dimension(xdim,ydim)   :: Ts, Ta, q, To, ice_Ts4, ice_H4
  real,   dimension(xdim,ydim,4) :: ice_T4
  real*8, dimension(xdim,ydim)   :: ice_Ts, ice_H
  real*8, dimension(xdim,ydim,4) :: ice_T
  real                           :: temp_r
  integer                        :: temp_i
  
  print *,'read restart'
  
  open(31,file='restart_in.txt')
  read(31,*) it, year
  close(31)

  open(31,file='restart_in.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  read(31,rec=1) Ts
  read(31,rec=2) Ta
  read(31,rec=3) q
  read(31,rec=4) To
  read(31,rec=5) ice_Ts4
  read(31,rec=6) ice_H4
  read(31,rec=7) ice_T4(:,:,1)
  read(31,rec=8) ice_T4(:,:,2)
  read(31,rec=9) ice_T4(:,:,3)
  read(31,rec=10) ice_T4(:,:,4)     
  do irec=1, nstep_yr
     read(31,rec=4*irec+8)  TF_correct(:,:,irec)
     read(31,rec=4*irec+9)  qF_correct(:,:,irec)
     read(31,rec=4*irec+10)  ToF_correct(:,:,irec)
     read(31,rec=4*irec+11)  precip_correct(:,:,irec)
  end do
  close(31)
  ice_Ts = ice_Ts4; ice_H = ice_H4; ice_T = ice_T4
  print *,'restart, Tsurf:', Ts(48,24), Ts(30,40)

  call orbital_forcing(year)
  
  print *,'restart year: ',year
    
  if(log_exp==330) then
    do irec=1,year/100
        read(53,*) t, CO2_SCN 
    end do
  end if 
  if(log_exp==311) then
    do irec=1,year/100
        read(53,*) temp_i,temp_r
    end do
    do irec=1,year/100
        read(54,*) temp_i,temp_r
    end do
  endif

end subroutine restart_read

!+++++++++++++++++++++++++++++++++++++++
subroutine output(it, year, iunit, irec, mon, ts0, ta0, to0, q0, dq_rain, &
                      ice_Ts0, ice_H0, term_mass, term_hadv, term_calv, ice_T0, &
                      ice_vx3, ice_vy3, a_surf)
!+++++++++++++++++++++++++++++++++++++++
! output variables in model
  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, time_scnr &
                           , time_ctrl, ireal, dt, dt_out 
  USE mo_physics,      ONLY: jday, r_qviwv, wz_vapor, iyear, glacier_type, mask, z_topo, ssh
  use mo_diagnostics,  ONLY: ice_Tsmm, ice_Tsmn_ctrl, ice_Hmn_ctrl, ice_mask_ctrl &
                           , Tmn_ctrl, Tamn_ctrl, Tomn_ctrl, qmn_ctrl, prmn_ctrl &
                           , Tmm, Tamm, Tomm, qmm, prmm, albdmn & 
                           , term_massmm, term_mass_ctrl, term_hadvmm, term_hadv_ctrl & 
                           , term_calvmm, term_calv_ctrl, ice_vxm0, ice_vym0, ice_Tmm, ice_Hm0
  implicit none
  real*4, external              :: gmean
  real,   dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, dq_rain, dq_eva, dq_crcl
  real,   dimension(xdim,ydim)  :: a_surf 
  real*8, dimension(xdim,ydim)  :: ice_H0, ice_Ts0, ice_mask 
  real*8, dimension(xdim,ydim)  :: term_hadv, term_mass, term_calv
  real*8, dimension(xdim,ydim,4):: ice_T0                ! ice strata temperature (current) [k] 
  real*8, dimension(xdim,ydim,4):: ice_vx3, ice_vy3      ! ice strata velocity at zonal (x) and meridian (y) direction [m/s]
  integer, parameter          :: nvar  = 19              ! number of output variable
  integer, parameter          :: ngvar = 7               ! number of output variable
  integer                     :: i,it,irec,mon,iyrec     ! work variable for count, controled by subroutine output
  integer                     :: ndm                     ! total time for mean calculation
  integer                     :: iunit                   ! written file uint
  integer                     :: year                    ! date variable year

  ! diagnostics: monthly means
  Tmm=Tmm+Ts0; Tamm=Tamm+ta0; Tomm=Tomm+to0; qmm=qmm+q0; 
  albdmn=albdmn+a_surf
  prmm=prmm+dq_rain*(r_qviwv*wz_vapor);          ! kg/m2/s
  ice_Tsmm = ice_Tsmm+ice_Ts0
  ice_Hm0=ice_H0;ice_Tmm = ice_Tmm+ice_T0 
  term_massmm = term_massmm + term_mass
  term_hadvmm = term_hadvmm + term_hadv
  term_calvmm = term_calvmm + term_calv
  
! control output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. iunit == 101 ) then
     ndm=jday_mon(mon)*ndt_days
     if (it/float(ndt_days)  > 365*(time_ctrl-1)) then
     ! for 301 experiment
     !if( mod(year, dt_out) == 0) then
         write(iunit,rec=nvar*irec+1)  Tmm/ndm
         write(iunit,rec=nvar*irec+2)  Tamm/ndm
         write(iunit,rec=nvar*irec+3)  Tomm/ndm
         write(iunit,rec=nvar*irec+4)  qmm/ndm
         write(iunit,rec=nvar*irec+5)  real(mask, kind=4)
         write(iunit,rec=nvar*irec+6)  prmm/ndm
         write(iunit,rec=nvar*irec+7)  albdmn/ndm
         write(iunit,rec=nvar*irec+8)  real(glacier_type, kind=4)
         write(iunit,rec=nvar*irec+9)  ice_Hm0
         write(iunit,rec=nvar*irec+10) z_topo
         write(iunit,rec=nvar*irec+11) term_massmm
         write(iunit,rec=nvar*irec+12) term_hadvmm
         write(iunit,rec=nvar*irec+13) term_calvmm
         write(iunit,rec=nvar*irec+14) ice_vxm0(:,:,4)
         write(iunit,rec=nvar*irec+15) ice_vym0(:,:,4)
         do i=1,4
              write(iunit,rec=nvar*irec+i+15) ice_Tmm(:,:,5-i)/ndm
         end do
         irec=irec+1; ! set in outside loop
     end if
     Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; prmm=0.
     ice_Tsmm=0.; term_massmm=0.; term_hadvmm=0.; term_calvmm=0.
     ice_Tmm=0.; ice_Hm0=0.; albdmn=0.
     mon=mon+1; if (mon==13) mon=1 ! set in outside loop
  end if

! scenario output
  if ( jday == sum(jday_mon(1:mon)) &
      .and. it/float(ndt_days) == nint(it/float(ndt_days)) ) then
      if(iunit == 102) then
          if( mod(year, dt_out) == 0) then
              ndm=jday_mon(mon)*ndt_days
              write(iunit,rec=nvar*irec+1)  Tmm/ndm
              write(iunit,rec=nvar*irec+2)  Tamm/ndm
              write(iunit,rec=nvar*irec+3)  Tomm/ndm
              write(iunit,rec=nvar*irec+4)  qmm/ndm
              write(iunit,rec=nvar*irec+5)  real(mask, kind=4)
              write(iunit,rec=nvar*irec+6)  prmm/ndm
              write(iunit,rec=nvar*irec+7)  albdmn/ndm
              write(iunit,rec=nvar*irec+8)  real(glacier_type, kind=4)
              write(iunit,rec=nvar*irec+9)  ice_Hm0
              write(iunit,rec=nvar*irec+10) z_topo
              write(iunit,rec=nvar*irec+11) term_massmm
              write(iunit,rec=nvar*irec+12) term_hadvmm
              write(iunit,rec=nvar*irec+13) term_calvmm
              write(iunit,rec=nvar*irec+14) ice_vxm0(:,:,4)
              write(iunit,rec=nvar*irec+15) ice_vym0(:,:,4)
              do i=1,4
                  write(iunit,rec=nvar*irec+i+15) ice_Tmm(:,:,5-i)/ndm
              end do
              
              ! time series output
              iyrec=floor(float(irec/12))
              write(103,rec=ngvar*irec+1) real(gmean(Tmm/ndm), kind=4)
              write(103,rec=ngvar*irec+2) real(gmean(Tamm/ndm), kind=4)
              write(103,rec=ngvar*irec+3) real(gmean(Tomm/ndm) , kind=4)
              write(103,rec=ngvar*irec+4) real(gmean(qmm/ndm), kind=4)
              write(103,rec=ngvar*irec+5) real(gmean(albdmn/ndm), kind=4)
              write(103,rec=ngvar*irec+6) real(gmean(prmm/ndm), kind=4)
              write(103,rec=ngvar*irec+7) real(ssh, kind=4) 
              irec=irec+1; ! set in outside loop
          end if
          Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; prmm=0.
          ice_Tsmm=0.; term_hadvmm=0.; term_massmm=0.; term_calvmm=0.
          ice_Tmm=0.; ice_Hm0=0.; albdmn=0.
          mon=mon+1; if (mon==13) mon=1 ! set in outside loop
       end if
  end if

end subroutine output
