!--------------------------------------------------------------------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model - Ice Sheet Model (ISM) Idealized Benchmark Experiment
!--------------------------------------------------------------------------------------------------------------------
!
!   Authors: Zhiang Xie
!
!   EISMINT I/II experiment

!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 144         ! field dimensions
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
  real               :: time_stp   = 1
  real               :: T          = 40000

  namelist / numerics / time_flux, time_ctrl, time_scnr, log_restart, time_out, dt_out, T, time_stp

end module mo_numerics


!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  integer  :: log_exp   = 0              ! process control logics for expiments (see header)
! deconstruct mean state (dmc) switches
  integer  :: log_cloud_dmc   = 1              ! process control clouds
  integer  :: log_ocean_dmc   = 1              ! process control ocean
  integer  :: log_atmos_dmc   = 1              ! process control Atmosphere
  integer  :: log_co2_dmc     = 1              ! process control CO2
  integer  :: log_hydro_dmc   = 1              ! process control hydro
  integer  :: log_qflux_dmc   = 1              ! process control qflux corrections
! deconstruct 2xco2 (drsp) switches
  integer  :: log_topo_drsp   = 1              ! process control for topo
  integer  :: log_cloud_drsp  = 1              ! process control for clouds
  integer  :: log_humid_drsp  = 1              ! process control for humidity clim
  integer  :: log_ocean_drsp  = 1              ! process control for ocean
  integer  :: log_hydro_drsp  = 1              ! process control for hydro
! switches that are the same for both deconstructions
  integer  :: log_ice         = 1              ! process control ice-albedo
  integer  :: log_snow        = 1              ! process control snow-albedo (ice sheet)
  integer  :: log_hdif        = 1              ! process control Diffusion of heat
  integer  :: log_hadv        = 1              ! process control Advection of heat
  integer  :: log_vdif        = 1              ! process control Diffusion of vapor
  integer  :: log_vadv        = 1              ! process control Advection of vapor
! switches for the hydrological cycle
  integer  :: log_rain        = 0              ! process control precipitation parameterisation
  integer  :: log_eva         = 0              ! process control evaporation parameterisation
  integer  :: log_conv        = 0              ! process control advection parameterisation
  integer  :: log_clim        = 0              ! process control for reference climatology
! switches for external forcing files
  integer  :: log_tsurf_ext   = 0              ! process control evaporation parameterisation
  integer  :: log_hwind_ext   = 0              ! process control advection parameterisation
  integer  :: log_omega_ext   = 0              ! process control for reference climatology
! switches for ice sheet model
  integer  :: log_ice_sheet   = 1              ! process control for ice sheet
  integer  :: log_ice_dyn     = 1              ! process control for ice sheet dynamic process
  integer  :: log_ice_ithk    = 1              ! ice thickness initial condition (0=ice free; 1=start from current ice thickness)
  integer  :: log_decouple    = 0              ! thermodynamic coupling switch (1 for decouple, 0 for coupling)
! switches for GREB-ice sheet coupling
  integer  :: log_ice_cpl     = 1              ! process control for GREB coupling to ice sheet
  integer  :: log_ice_topo    = 1              ! process control for ice topograph 
  integer  :: log_ice_slv     = 1              ! process control for sea level effect 
! switches for 
  integer  :: log_bound_cond  = 0              ! boundary condition control 
  integer  :: log_time_depen  = 0              ! time dependency control 

! parameters for scenarios
  real     :: dradius   = 0.		 ! deviations from actual earth radius in %

! physical parameter (natural constants)
  parameter( pi        = 3.1416 )
  parameter( sig       = 5.6704e-8 )     ! stefan-boltzmann constant [W/m^2/K^4]
  parameter( R_gas     = 8.314 )         ! universal gas constant [J/mol/K]
  parameter( rho_ocean = 999.1 )         ! density of water at T=15C [kg/m^2]
  parameter( rho_land  = 2600. )         ! density of solid rock [kg/m^2]
  parameter( rho_air   = 1.2 )           ! density of air at 20C at NN
  parameter( rho_ice   = 910.)           ! density of incompressible polycrystalline ice [kg/m^3]
  parameter( cp_ocean  = 4186. )         ! specific heat capacity of water at T=15C [J/kg/K]
  parameter( cp_land   = cp_ocean/4.5 )  ! specific heat capacity of dry land [J/kg/K]
  parameter( cp_air    = 1005. )         ! specific heat capacity of air      [J/kg/K]
  parameter( cp_ice    = 2009. )         ! specific heat capacity of ice [J/kg/K]
  parameter( eps       = 1. )            ! emissivity for IR
  real :: S0_var       = 100.            ! variation of solar constant   [%]

! physical parameter (model values)
  parameter( d_ocean   = 50. )                     ! depth of ocean column [m]
  parameter( d_land    = 2. )                      ! depth of land column  [m]
  parameter( d_air     = 5000. )                   ! depth of air column   [m]
  parameter( d_ice_max = 1. )                      ! maximum shortwave radiation penetarting depth of ice column   [m]
  parameter( d_ice_mov = 10. )                      ! maximum shortwave radiation penetarting depth of ice column   [m]
  parameter( cap_ocean = cp_ocean*rho_ocean )      ! heat capacity 1m ocean  [J/K/m^2]
  parameter( cap_land  = cp_land*rho_land*d_land ) ! heat capacity land   [J/K/m^2]
  parameter( cap_air   = cp_air*rho_air*d_air )    ! heat capacity air    [J/K/m^2]
  real :: ice_svlm  = 1./rho_ice                   ! specific volume of ice/snow [m3/kg]
  real :: ct_sens   = 22.5                         ! coupling for sensible heat
  real :: da_ice    = 0.25                         ! albedo diff for ice covered points
  real :: a_no_ice  = 0.1                          ! albedo for non-ice covered points
  real :: a_cloud   = 0.35                         ! albedo for clouds
  real :: a_snow    = 0.8                          ! process control snow-albedo (ice sheet)
  real :: Tl_ice1   = 273.15-10.                   ! temperature range of land snow-albedo feedback
  real :: Tl_ice2   = 273.15                       ! temperature range of land snow-albedo feedback
  real :: To_ice1   = 273.15-7.                    ! temperature range of ocean ice-albedo feedback
  real :: To_ice2   = 273.15-1.7                   ! temperature range of ocean ice-albedo feedback
  real :: co_turb   = 5.0                          ! turbolent mixing to deep ocean [W/K/m^2]
! original real :: kappa     = 8e5                          ! atmos. diffusion coefficient [m^2/s]
  real :: kappa     = 5*8e5                        ! atmos. diffusion coefficient [m^2/s]
  real :: ice_kappa = 2.1                          ! ice sheet diffusion coefficient [W/K/m] 
  real :: rock_kappa= 3.                           ! heat conductivity for lithosphere [W/K/m] 
  real :: beta_melt = 8.7e-4                       ! change of melting point with ice depth [K/m]
  real :: geoh_flux = 4.2e-2                       ! geoheat flux on bottom [W/m2] 
  real :: enh_fact  = 1.                           ! enhance factor for ice flow [1]
  !real :: slid_law  = 6e4/365/86400                ! no basal slide
  real :: slid_law  = 0.                ! no basal slide
  real :: ice_Tcons = 273.15-10.                   ! critical temperature for constitutive equation [K]
  ! real :: ref_lat   = 84.375                       ! reference latitude for polar filter [o]
  real :: ref_lat   = 76.875                       ! reference latitude for polar filter [o]
  parameter( ce        = 2e-3  )                   ! laten heat transfer coefficient for ocean
  parameter( cq_latent = 2.257e6 )                 ! latent heat of condensation/evapoartion f water [J/kg]
  parameter( ci_latent = 3.335e5 )                 ! latent heat of condensation/fusion f ice [J/kg]
  parameter( cq_rain   = -0.1/24./3600. )          ! decrease in air water vapor due to rain [1/s]
  parameter( z_air     = 8400. )                   ! scaling height atmos. heat, CO2
  parameter( z_vapor   = 5000. )                   ! scaling height atmos. water vapor diffusion
  parameter( c_lapse   = -0.0098 )                 ! diabatic temperature decrease rate in atmosphere [K/m]
  parameter( grav      = 9.81  )                   ! gravitational acceleration [m/s^2]
  real :: r_qviwv   = 2.6736e3                     ! regres. factor between viwv and q_air  [kg/m^3]
  real :: co2_scn = 0.                             ! CO2 in scenario run
 ! DiDi 
  real :: Todeep0   = 273.15+2.                    ! deep abyssal ocean initial temperature
 
 ! paleo runs with orbital forcing
  integer  :: kry_start = -1000                    ! start year in [kyr] has to be <=0   
  integer  :: kry_end   = 0                        ! end year in [kyr] has to be <=0  
  integer  :: dtyr_acc  = 1                        ! acceleration in sw_solar forcing 

  ! physical paramter (rainfall)
  real :: c_q, c_rq, c_omega, c_omegastd

! parameter emisivity
  real, dimension(10) :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
&                                             0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)
  real, dimension(2)  :: actene = (/1.39e5,6e4/)       ! constitutive equation coefficients, activate energy (>-10oC and <-10oC, Glen's flow law, unit: J mol-1)
  real, dimension(2)  :: A_fact = (/1.96e3,3.985e-13/) ! constitutive equation coefficients, softness parameter (>-10oC and <-10oC, Glen's flow law, unit: Pa-3 s-1)

! declare climate fields
  real*4, dimension(xdim,ydim)          ::  z_topo0, z_topo, glacier, z_ocean, mask
  real*4, dimension(xdim,ydim,nstep_yr) ::  iceH_clim
  real*4, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, omegaclim, omegastdclim, wsclim
  real*4, dimension(xdim,ydim,nstep_yr) ::  qclim, mldclim, Toclim, cldclim, precipclim
  real*4, dimension(xdim,ydim,nstep_yr) ::  Fn_surfclim, Taclim, dT_oceanclim
  real*4, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, precip_correct
  real*4, dimension(xdim,ydim,nstep_yr) ::  swetclim, dTrad
  real*4, dimension(ydim,nstep_yr)      ::  sw_solar, sw_solar_ctrl, sw_solar_scnr
  real, dimension(xdim,ydim)          ::  co2_part      = 1.0
  real, dimension(xdim,ydim)          ::  co2_part_scn  = 1.0

! declare anomaly fields for enso and climate change
  real*4, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_enso     = 0.
  real*4, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_enso     = 0.
  real*4, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_enso     = 0.
  real*4, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_enso = 0.
  real*4, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_enso    = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_cc   = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_cc      = 0.

! declare constant fields
  real, dimension(xdim,ydim)          ::  cap_surf
! ice sheet : capacity of snow  
  real, dimension(xdim,ydim)          ::  cap_ice
  real, dimension(xdim,ydim)          ::  glacier_type
  real                                ::  shlf_eta, ice_shlf 
  integer jday, ityr, iyear, it_start

! Mike: declare some program constants
  real, dimension(xdim,ydim,nstep_yr) :: uclim_m, uclim_p
  real, dimension(xdim,ydim,nstep_yr) :: vclim_m, vclim_p
! ice sheet related topograph variable
  real, dimension(xdim, ydim)         :: wz_air, wz_vapor
  real*4, dimension(xdim, ydim)       :: b_rock, b_rock0
  real                                :: ssh    ! sea level change [m]


  real :: t0, t1, t2, xtest

  namelist / physics / log_exp, ct_sens, da_ice, a_no_ice, a_cloud, co_turb, kappa, 	&
&                      p_emi, Tl_ice1, Tl_ice2, To_ice1, To_ice2, r_qviwv,          	&
&		                   log_cloud_dmc, log_ocean_dmc, log_atmos_dmc, log_co2_dmc,    &
&                      log_hydro_dmc, log_qflux_dmc, 					                &
&                      log_topo_drsp, log_cloud_drsp, log_humid_drsp, log_hydro_drsp,   &
&                      log_ocean_drsp, log_ice, log_hdif, log_hadv, log_vdif, log_vadv, &
& 		                 S0_var, dradius, log_rain, log_eva, log_conv, log_clim,        &
&                      log_tsurf_ext, log_hwind_ext, log_omega_ext,                     & 
&                      log_ice_sheet, log_ice_ithk, log_restart_r, log_restart_w,       & 
&                      kry_start, kry_end, dtyr_acc, Todeep0, log_ice_cpl, co2_scn,     &
&                      log_ice_dyn, log_ice_topo, log_ice_slv,  &
&                      log_bound_cond, log_time_depen, log_decouple

end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim

 ! declare diagnostic fields
  real*4, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
&                                        ftmn, fqmn, icmn, Tomn, ice_Hmn, ice_Tsmn

! declare output fields
  real*4, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, icmm, prmm, evamm, qcrclmm, &
&                                          ice_Hm0, ice_Tsmm, term_massmm, term_hadvmm,    & 
&								  	       term_calvmm, albdmn, xdum
  real*4, dimension(xdim,ydim,4)        :: ice_Tmm, ice_vxm0, ice_vym0, ice_vzm0
  real*4, dimension(xdim,ydim,12)       :: Tmn_ctrl, Tamn_ctrl, Tomn_ctrl, ice_Hmn_ctrl, ice_Tsmn_ctrl, ice_mask_ctrl
  real*4, dimension(xdim,ydim,12)       :: qmn_ctrl, icmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl 
  real*4, dimension(xdim,ydim,12)       :: term_mass_ctrl, term_hadv_ctrl, term_calv_ctrl 
  real, dimension(xdim,ydim)            :: ice_Tsm 

end module mo_diagnostics

include "ice-sheet.f90"

!+++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++
! EISMINT experiments  
program EISMINT_exp
  use mo_numerics
  use mo_physics
  use mo_diagnostics
  implicit none

  integer                        :: i,j,k
  integer                        :: irec, it, year
  integer                        :: output_count
  integer                        :: it_max
  integer                        :: nvar = 13
  integer, parameter         :: kdim = 4                ! vertical layer
  real, dimension(xdim,ydim)     :: ice_H1, ice_H0
  real, dimension(xdim,ydim)     :: thk_layer           ! ice thickness for specific layer to bottom [m]
  real, dimension(xdim,ydim)     :: thk_adv             ! ice thickness advection at layer [m] 
  real, dimension(xdim,ydim)     :: term_mass           ! mass balance
  real, dimension(xdim,ydim)     :: term_calv           ! mass balance
  real, dimension(xdim,ydim)     :: ice_Ts1             ! ice surface temperature
  real, dimension(xdim,ydim)     :: ice_zs              ! ice surface elevation
  real, dimension(xdim,ydim)     :: term_hadv           ! ice thickness tendency due to advection [m]
  real, dimension(xdim,ydim)     :: termx, termy        ! ice advection term x/y direction (work variable)
  real, dimension(xdim,ydim)     :: term_advx           ! ice advection term x direction (work variable)
  real, dimension(xdim,ydim)     :: term_advy           ! ice advection term y direction (work variable)
  real, dimension(xdim,ydim)     :: dice_Ts             ! ice surface temperature tendency
  real, dimension(xdim,ydim,4)   :: divg                ! horizontal divergence for half year time step [1] 
  real, dimension(xdim,ydim,4)   :: adv_cord            ! advection term due to coordinate transformation (dvx/dz dh/dx + dvy/dz dh/dy) [s-1] 
  real, dimension(xdim,ydim,4)   :: term_sig            ! ice strata temperature tendency due to strain rate [K]
  real, dimension(xdim,ydim,4)   :: term_dif            ! ice strata temperature tendency due to diffusion [K]
  real, dimension(xdim,ydim,4)   :: term_tadv           ! ice strata temperature tendency due to advection [K]
  real, dimension(xdim,ydim,4)   :: term_vadv           ! ice strata temperature tendency due to vertical advection [K]
  real, dimension(xdim,ydim,4)   :: dice_T              ! ice strata temperature tendency [K]
  real, dimension(xdim,ydim,4)   :: ice_T1              ! ice temperature
  real, dimension(xdim,ydim,4)   :: ice_T0              ! ice temperature
  real, dimension(xdim,ydim,4)   :: ice_Td              ! ice temperature for diffusion calculation
  real, dimension(xdim,ydim,4)   :: ice_Tmt             ! ice strata temperature for melting point [k]
  real, dimension(xdim,ydim,4)   :: ice_Tcoef           ! ice temperature
  real, dimension(xdim,ydim,4)   :: ice_vz 
  real, dimension(xdim,ydim,4)   :: ice_vzg                ! vertical grid average of ice velocity at vertiacl direction [m/s]
  real, dimension(xdim,ydim+1)   :: ice_vx, ice_vy         ! ice velocity
  real, dimension(xdim,ydim+1)   :: ice_vxvm, ice_vyvm     ! ice mean velocity
  real, dimension(xdim,ydim,kdim)   :: ice_vx3             ! ice strata velocity at zonal (x) direction [m/s]
  real, dimension(xdim,ydim,kdim)   :: ice_vy3             ! ice strata velocity at meridian (y) direction [m/s]
  real, dimension(xdim,ydim,kdim)   :: sigma_x                 ! x-z component of stress tensor [Pa]
  real, dimension(xdim,ydim,kdim)   :: sigma_y                 ! y-z component of stress tensor [Pa]
  real, dimension(xdim,ydim,kdim)   :: dvxdz                   ! vertical velocity gradient at zonal direction 
  real, dimension(xdim,ydim,kdim)   :: dvydz                   ! vertical velocity gradient at meridian direction 
  real, dimension(ydim)          :: lat, dxlat
  real                           :: dx,dy,dyy,deg
  real,parameter                 :: gs_layer = 1/sqrt(3.)
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)  ! index for latitude grid
  real,dimension(4)              :: zeta = (/-1.,-gs_layer,gs_layer,1./) ! model layer for output 
  real,dimension(4)              :: zeta_c                                 ! model layer for numeric calculation

  open(10, file='namelist')
  read(10,numerics)
  read(10,physics)
  close(10)
  open(102,file='./scenario.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

  ! deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx  = dlon; dy=dlat; dyy=dy*deg
  dyy = 50000
  deg = dyy/dy
  dt_ice = dt_ice*time_stp
  if(log_decouple==0) dyy=25000
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  !zeta_c = zeta_o
  !zeta_c(1) = (zeta_o(1)+zeta_o(2))/2.
  
  ice_Ts1 = 0.; term_mass = 0.; b_rock = 0.; glacier_type = -1.; ice_H1 = 0.
  irec = 0.; year = 0; output_count = 0
  it_max = time_scnr/time_stp
  
  if(log_decouple==0) then
      call boundary_condition_I(term_mass, ice_Ts1, ice_H1, year, dyy)
  else
      call boundary_condition_II(term_mass, ice_Ts1, ice_H1, year, dyy)
  end if

  do k=1,4
      ice_T1(:,:,k)=ice_Ts1
  end do
  
  ! read data for EISMINT I
  if(log_time_depen==1.and.log_decouple==1) then
    if(log_bound_cond==0) then
      open(103,file='../EISMINT_I_exp1.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    else 
      open(103,file='../EISMINT_I_exp2.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    end if
    read(103,rec=(time_scnr/time_out-1)*nvar+1) ice_Hm0
    do k=1,4
      read(103,rec=(time_scnr/time_out-1)*nvar+1+k) ice_Tmm(:,:,5-k)
    end do
    close(103)
    ice_H1 = ice_Hm0
    ice_T1 = ice_Tmm
  end if
  
  ! read data for EISMINT II
  if(log_restart_r==0 .and. log_bound_cond > 1 .and. log_decouple==0) then
    open(103,file='../EISMINT_II_exp1.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
    read(103,rec=(time_scnr/time_out-1)*nvar+1) ice_Hm0
    do k=1,4
      read(103,rec=(time_scnr/time_out-1)*nvar+1+k) ice_Tmm(:,:,5-k)
    end do
    ice_H1 = ice_Hm0
    ice_T1 = ice_Tmm
  end if

  if(log_restart_r==1) then
      open(103,file='./restart.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
      read(103,rec=1) ice_Hm0
      !print*,ice_H1(1,1)
      do k=1,4
          read(103,rec=1+k) ice_Tmm(:,:,5-k)
      end do
      close(103)
    ice_H1 = ice_Hm0
    ice_T1 = ice_Tmm
      open(104,file='./iteration_info',form='formatted')
      read(104,*) it_start,year
      close(104)
  end if
  
  do it=1,it_max

      if(log_decouple==1) call boundary_condition_I(term_mass, ice_Ts1, ice_H1, year, dyy)
      if(log_decouple==0) call boundary_condition_II(term_mass, ice_Ts1, ice_H1, year, dyy)
      ice_zs    = b_rock + ice_H1
    
      term_hadv = 0.; term_calv = 0.; dice_Ts = 0.
      term_tadv = 0.; term_dif  = 0.; term_sig  = 0.; dice_T  = 0.; term_vadv = 0.
    !+++++++++++++++++++++++++++++ dynamic core  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! ice temperature profile estimate and vertical diffusion
      ice_T1(:,:,4) = ice_Ts1
      do j = 1,xdim
         do i = 1,ydim
            call ice_regression_coef(ice_T1(j,i,:), ice_Tcoef(j,i,:))
         end do
      end do
      
      do k = 1,4
          ! varible initialization 
          thk_layer = 0.; 
          ice_vx = 0.; ice_vy = 0.; ice_vxvm = 0.; ice_vyvm = 0.
          ! diagnostic variable in ice sheet model, including strata/vertical mean ice horizontal velocity, strain rate heating
          call ice_sheet_diagnostics(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, &
          &                          ice_vxvm, ice_vyvm, zeta(k), sigma_x(:,:,k), sigma_y(:,:,k))
          ice_vx3(2:xdim,:,k) = (ice_vx(1:xdim-1,1:ydim)+ice_vx(2:xdim,1:ydim))/2.;
          ice_vx3(1,:,k) = (ice_vx(xdim,1:ydim)+ice_vx(1,1:ydim))/2.
          ice_vy3(:,:,k) = (ice_vy(:,1:ydim)+ice_vy(:,2:ydim+1))/2.
          ice_vym0(:,:,k) = ice_vy(:,2:ydim+1)
          ! ice thickness advection
          thk_layer = ice_H1*(zeta(k)+1)/2.
          call ice_advection(ice_vxvm, ice_vyvm, thk_layer, divg(:,:,k), dyy, dxlat, lat)
          ! advection term for ice thickness equation
          if(k==4) term_hadv = divg(:,:,k)
          call ice_advection_upwind(ice_vx, ice_vy, thk_layer, thk_layer, term_advx, term_advy, &
          &                         dyy, dxlat, lat, ice_H1)
          thk_adv = - (term_advx + term_advy)
          ! vertical velocity
          ice_vz(:,:,k) = (divg(:,:,k) - thk_adv)/dt_ice
          if(k > 1) ice_vzg(:,:,k) = (ice_vz(:,:,k-1)+ice_vz(:,:,k))/2.
          
          ! strain rate heating (stress form)
          ! term_sig(:,:,k) = dt_ice*term_sig(:,:,k)/(cp_ice*rho_ice)
          ! ice temperature advection
          call ice_advection_upwind(ice_vx, ice_vy, ice_T1(:,:,k), ice_T1(:,:,k), term_advx, term_advy, &
          &                         dyy, dxlat, lat, ice_H1)
          term_tadv(:,:,k) = - (term_advx + term_advy)
      end do
      ! strain rate heating (velocity form)
      where(ice_H1> d_ice_mov) dvxdz(:,:,1) = (ice_vx3(:,:,2) - ice_vx3(:,:,1))/(zeta(2)-zeta(1))*2./ice_H1
      where(ice_H1> d_ice_mov) dvydz(:,:,1) = (ice_vy3(:,:,2) - ice_vy3(:,:,1))/(zeta(2)-zeta(1))*2./ice_H1
      do k = 2,kdim-1
          where(ice_H1> d_ice_mov) dvxdz(:,:,k) = (ice_vx3(:,:,k+1) - ice_vx3(:,:,k-1))/(zeta(k+1)-zeta(k-1))*2./ice_H1
          where(ice_H1> d_ice_mov) dvydz(:,:,k) = (ice_vy3(:,:,k+1) - ice_vy3(:,:,k-1))/(zeta(k+1)-zeta(k-1))*2./ice_H1
      end do
      term_sig = dt_ice*(sigma_x*dvxdz+sigma_y*dvydz)/(cp_ice*rho_ice)
      ! vertical diffusion and convection
      ice_Td = ice_T1 + term_sig + term_tadv 
      do j = 1,xdim
         do i = 1,ydim
            call ice_temperature_convection_diffusion(ice_Td(j,i,:), ice_H1(j,i),  ice_vz(j,i,:),ice_vzg(j,i,:), &
            &                                         zeta, term_dif(j,i,:))
         end do
      end do

      ice_H0 = ice_H1 + term_mass*time_stp + term_hadv
      where(glacier_type < 0.) ice_H0=0.
      where(ice_H0 < 0.) ice_H0=0.
    
      ! ice temperature equation
      dice_T = term_tadv + term_sig + term_dif
      do k = 1,4
          ice_Tmt(:,:,k) = Tl_ice2 - beta_melt*(1-zeta(k))*ice_H0/2.
      end do
      where(dice_T > ice_Tmt-ice_T1) dice_T = 0.9*(ice_Tmt-ice_T1) ! numeric stability
      ice_T0 = ice_T1 + dice_T
      do k = 1,4
          where(ice_H0 <= d_ice_mov .or. (ice_H1 <= d_ice_mov .and. ice_H0 > d_ice_mov)) ice_T0(:,:,k) = ice_Ts1 
      end do

      ice_H1 = ice_H0
      ice_T1 = ice_T0

      year = year + time_stp
      output_count = output_count + time_stp
      ice_Hm0=ice_H0; ice_Tmm=ice_T0; ice_vzm0=ice_vz;
      if(output_count >= time_out) then
          write(102,rec=irec*nvar+1) ice_Hm0
          do k=1,4
              write(102,rec=irec*nvar+1+k) ice_Tmm(:,:,5-k)
          end do
          do k=1,4
              write(102,rec=irec*nvar+5+k) ice_vym0(:,:,5-k)
          end do
          do k=1,4
              write(102,rec=irec*nvar+9+k) ice_vzm0(:,:,5-k)
          end do
          irec = irec + 1
          output_count = 0
      end if
      print*,'year',year,'ice thickness 1', ice_H0(1,1),'ice thickness 2', ice_H0(1,16)
      print*,'16-20 vel ',ice_vym0(1,16:20,4)*86400*365
      print*,'21-25 vel ',ice_vym0(1,21:25,4)*86400*365
      print*,'ts1', ice_Tmm(1,1,:)
      print*,'ts8', ice_Tmm(1,16,:)
      print*,'vz1',ice_vzm0(1,1,:)*86400*365 
      print*,'vz2',ice_vzm0(1,16,:)*86400*365
      print*,'vy1',ice_vym0(1,1,:)*86400*365
      print*,'vy8',ice_vym0(1,16,:)*86400*365
      print*,'mass   ',term_mass(1,1),term_mass(1,16) 
      print*,'thk adv',term_hadv(1,1),term_hadv(1,16)
      print*,'hadv1', term_tadv(1,1,:)
      print*,'hadv8', term_tadv(1,16,:)
      print*,'diff1', term_dif(1,1,:)
      print*,'diff8', term_dif(1,16,:)
      print*,'sigma1',term_sig(1,1,:)
      print*,'sigma8',term_sig(1,16,:)
  end do
  close(102)
  if(log_restart_w==1) then
      ice_Hm0=ice_H0; ice_Tmm=ice_T0; ice_vzm0=ice_vz;
      open(105,file='./restart.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
      write(105,rec=1) ice_Hm0
      do k=1,4
          write(105,rec=1+k) ice_Tmm(:,:,5-k)
      end do
      close(105)
      open(106,file='./iteration_info',form='formatted')
      write(106,*) it,year
      close(106)
  end if
end program

!+++++++++++++++++++++++++++++++++++++++
subroutine boundary_condition_I(term_mass, ice_Ts1, ice_H1, year, dyy)
!+++++++++++++++++++++++++++++++++++++++
use mo_numerics
use mo_physics
implicit none
integer i,j,k
  integer                        :: year
  real, dimension(xdim,ydim)     :: ice_Ts1                                     ! ice surface temperature
  real, dimension(xdim,ydim)     :: term_mass                                   ! mass balance
  real, dimension(xdim,ydim)     :: ice_H1                                   ! ice thickness
  real, dimension(ydim)          :: lat, dxlat
  real                           :: dyy
  real                           :: s,Rel,d
  real,parameter                 :: gs_layer = 1/sqrt(3.)
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)  ! index for latitude grid
  s = 1d-2; Rel=450;
  if(log_time_depen==0) then
      if(log_bound_cond==0) then
          do j=1,xdim
             do i=1,16
                term_mass(j,i) = 0.3
                ice_Ts1(j,i) = 239 + 8d-8*(i*dyy/1000)**3
                glacier_type(j,i) = 1.
             end do
          end do
      else
          do j=1,xdim
             do i=1,16
                term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
                ice_Ts1(j,i) = 270 - 0.01*ice_H1(j,i)
                glacier_type(j,i) = 1.
             end do
          end do
      end if
  else
      if(log_bound_cond==0) then
          do j=1,xdim
             do i=1,16
                term_mass(j,i) = (0.3+0.2*sin(2*pi/T*year))
                ice_Ts1(j,i) = 239 + 8d-8*(i*dyy/1000)**3 + 10*sin(2*pi*year/T)
                glacier_type(j,i) = 1.
             end do
          end do
      else
          do j=1,xdim
             do i=1,16
                Rel = 450 + 100*sin(2*pi/T*year)
                term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
                ice_Ts1(j,i) = 270 - 0.01*ice_H1(j,i)
                glacier_type(j,i) = 1.
             end do
          end do
      end if
  end if
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine boundary_condition_II(term_mass, ice_Ts1, ice_H1, year, dyy)
!+++++++++++++++++++++++++++++++++++++++
use mo_numerics
use mo_physics
implicit none
integer i,j,k
  integer                        :: year
  real, dimension(xdim,ydim)     :: ice_Ts1                                     ! ice surface temperature
  real, dimension(xdim,ydim)     :: term_mass                                   ! mass balance
  real, dimension(xdim,ydim)     :: ice_H1                                   ! ice thickness
real, dimension(ydim)          :: lat, dxlat
real                           :: dyy
real                           :: s,Rel,d
real,parameter                 :: gs_layer = 1/sqrt(3.)
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)  ! index for latitude grid
  s = 1d-2; Rel=450;
  if(log_bound_cond==1) then
      do j=1,xdim
         do i=1,31
            term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
            ice_Ts1(j,i) = 238.15 + 1.67e-2*(i*dyy/1000)
            glacier_type(j,i) = 1.
         end do
      end do
  else if(log_bound_cond==2) then
      do j=1,xdim
         do i=1,31
            term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
            ice_Ts1(j,i) = 243.15 + 1.67e-2*(i*dyy/1000)
            glacier_type(j,i) = 1.
         end do
      end do
  else if(log_bound_cond==3) then
      Rel=425;
      do j=1,xdim
         do i=1,31
            term_mass(j,i) = min(0.25, s*(Rel-i*dyy/1000))
            ice_Ts1(j,i) = 238.15 + 1.67e-2*(i*dyy/1000)
            glacier_type(j,i) = 1.
         end do
      end do
  else if(log_bound_cond==4) then
      Rel=425;
      do j=1,xdim
         do i=1,31
            term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
            ice_Ts1(j,i) = 238.15 + 1.67e-2*(i*dyy/1000)
            glacier_type(j,i) = 1.
         end do
      end do
  else if(log_bound_cond==5) then
      do j=1,xdim
         do i=1,31
            term_mass(j,i) = min(0.5, s*(Rel-i*dyy/1000))
            ice_Ts1(j,i) = 223.15 + 1.67e-2*(i*dyy/1000)
            glacier_type(j,i) = 1.
         end do
      end do
  end if
end subroutine

!TB
!+++++++++++++++++++++++++++++++++++++++
function gmean(data)
!+++++++++++++++++++++++++++++++++++++++

use mo_numerics,		ONLY: xdim, ydim, dlat

! declare variables
real*4, dimension(xdim,ydim) 	:: data, w
real*4                          :: gmean
real, dimension(ydim)		:: lat

do i=1,ydim
lat(i) = -90-(dlat*0.5)+i*dlat
end do
do i=1,xdim
w(i,:) = cos(2.*3.14159*lat/360.)
end do

gmean = sum(data*w)/sum(w)

end function

