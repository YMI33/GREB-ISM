! ice sheet model for GREB
! author: Zhiang Xie
! start 26 Sep 2019

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_sheet(it, ionum, irec, mon, ice_H1, ice_T1, ice_Ts1, Ta1, dT_ocean,     &
&                    Fn_surf, dq_rain, ice_H0, ice_T0, ice_Ts0, term_mass,   & 
&                    term_hadv, term_calv, ice_fus, year)
!+++++++++++++++++++++++++++++++++++++++
! ice sheet model main routine 
! output: ice thickness budget (term_mass, term_hadv); calving diagnostic (term_calv); ice latent heat (ice_fus)
! implicit output: surface elevation (z_topo); glacier type (glacier_type); ice surface velocity (ice_vxm0, ice_vym0) 
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
                         dlon, dlat, jday_mon, nstep_yr, dt_ice
  USE mo_physics, ONLY: pi, cp_ice, rho_ice, rho_ocean,d_ice_max, d_ice_mov, glacier_type, &
                        log_exp, Tl_ice2, jday, ityr, iceH_clim, b_rock, z_topo, mask, &
                        beta_melt, log_ice_dyn, log_ice_topo, To_ice2
  USE mo_diagnostics
  implicit none
  real,parameter             :: gs_layer = 1/sqrt(3.)   
  integer, parameter         :: kdim = 4                ! vertical layer
  real*8,dimension(kdim)     :: zeta = (/-1e0,-gs_layer,gs_layer,1e0/) ! model layer for output 
  integer                    :: ice_iunit               ! file unit of ice sheet data
  integer                    :: ionum                   ! file unit of GREB
  integer                    :: i, j, k, it, irec       ! work variable for count
  integer                    :: mon, year               ! date variable
  
  real*8                       :: deg                      ! degree-length conversion [m/deg]
  real*8                       :: dx                       ! degree difference for zonal direction
  real*8                       :: dy                       ! degree difference for meridian direction
  real*8                       :: dyy                      ! meridional grid width [m]
  integer, dimension(ydim)     :: ilat = (/(i,i=1,ydim)/)  ! index for latitude grid
  real*8, dimension(ydim)      :: lat                      ! latitude [deg]
  real*8, dimension(ydim)      :: dxlat                    ! zonal grid width [m]

  ! input variable
  real*8, dimension(xdim,ydim) :: Fn_surf                 ! surface net heat flux [J/m2]
  real,   dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real,   dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real,   dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]

  real*8, dimension(xdim,ydim) :: ice_H0                  ! ice thickness (forward) [m]
  real*8, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (current) [m]
  real*8, dimension(xdim,ydim) :: ice_zs                  ! ice surface height [m]
  real*8, dimension(xdim,ydim) :: thk_layer               ! ice thickness for specific layer to bottom [m]
  real*8, dimension(xdim,ydim) :: term_advx               ! ice advection term x direction (work variable)
  real*8, dimension(xdim,ydim) :: term_advy               ! ice advection term y direction (work variable)
  real*8, dimension(xdim,ydim) :: thk_adv                 ! ice thickness advection at layer [m] 

  real*8, dimension(xdim,ydim) :: ice_Ts0                 ! ice surface temperature (forward) [K]
  real*8, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature (current) [K]
  real*8, dimension(xdim,ydim) :: ice_Tmtb                ! basal ice temperature for melting point [k] 

  real*8, dimension(xdim,ydim) :: term_mass               ! ice thickness tendency due to mass balance [m] 
  real*8, dimension(xdim,ydim) :: term_calv               ! ice thickness tendency due to calving [m] 
  real*8, dimension(xdim,ydim) :: term_hadv               ! ice thickness tendency due to advection [m]
  real*8, dimension(xdim,ydim) :: dice_Ts                 ! ice surface temperature tendency due to energy balance [K] 
  real,   dimension(xdim,ydim) :: ice_fus                 ! ice melting latent heat [W/m2]

  real*8, dimension(xdim,ydim+1)   :: ice_vx, ice_vy      ! ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim+1)   :: ice_vxvm, ice_vyvm  ! vertical mean of ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim,kdim)   :: ice_T1              ! ice strata temperature (current) [k] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_Td              ! ice strata temperature (used for convection-diffusion calculation) [k] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_T0              ! ice strata temperature (forward) [k] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_Tmt             ! ice strata temperature for melting point [k] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_Tcoef           ! regression coefficient of ice layer temperature [k] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_vx3             ! ice strata velocity at zonal (x) direction [m/s]
  real*8, dimension(xdim,ydim,kdim)   :: ice_vy3             ! ice strata velocity at meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim,kdim)   :: sigma_x                 ! x-z component of stress tensor [Pa]
  real*8, dimension(xdim,ydim,kdim)   :: sigma_y                 ! y-z component of stress tensor [Pa]
  real*8, dimension(xdim,ydim,kdim)   :: dvxdz                   ! vertical velocity gradient at zonal direction 
  real*8, dimension(xdim,ydim,kdim)   :: dvydz                   ! vertical velocity gradient at meridian direction 

  real*8, dimension(xdim,ydim,kdim)   :: divg                ! total percentage of mass change due to horizontal divergence for half year time step [1] 
  real*8, dimension(xdim,ydim,kdim)   :: adv_cord            ! advection term due to coordinate transformation (dvx/dz dh/dx + dvy/dz dh/dy) [s-1] 
  real*8, dimension(xdim,ydim,kdim)   :: ice_vz              ! ice velocity at vertiacl direction [m/s]
  real*8, dimension(xdim,ydim,kdim)   :: ice_vzg             ! vertical grid average of ice velocity at vertiacl direction [m/s]
  real*8, dimension(xdim,ydim,kdim)   :: term_sig            ! ice strata temperature tendency due to strain rate [K]
  real*8, dimension(xdim,ydim,kdim)   :: term_dif            ! ice strata temperature tendency due to diffusion [K]
  real*8, dimension(xdim,ydim,kdim)   :: term_tadv           ! ice strata temperature tendency due to advection [K]
  real*8, dimension(xdim,ydim,kdim)   :: term_vadv           ! ice strata temperature tendency due to vertical advection [K]
  real*8, dimension(xdim,ydim,kdim)   :: dice_T              ! ice strata temperature tendency [K]
  
  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)

  ice_zs = 0.
  if(log_ice_topo .eq. 1) then
      where(glacier_type >  0) ice_zs = b_rock + ice_H1
      where(glacier_type == 0) ice_zs = ice_H1*(1.-rho_ice/rho_ocean)
  else
      where(glacier_type >  0) ice_zs = b_rock + iceH_clim(:,:,ityr)
      where(glacier_type == 0) ice_zs = iceH_clim(:,:,ityr)*(1.-rho_ice/rho_ocean)
  end if

  term_hadv = 0.; term_mass = 0.; term_calv = 0.; dice_Ts = 0.
  term_tadv = 0.; term_dif  = 0.; term_sig  = 0.; dice_T  = 0.; term_vadv = 0.
  ice_Tsm = ice_Tsm + ice_Ts1
!+++++++++++++++++++++++++++++ dynamic core  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ( jday == sum(jday_mon) .and. log_ice_dyn == 1 &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) ) then
      dvxdz = 0.; dvydz = 0.; sigma_x = 0.; sigma_y = 0.
      ! ice temperature profile estimate, vertical diffusion and pressure-based melting point
      ice_T1(:,:,4) = ice_Tsm/nstep_yr
      ice_Tmtb = Tl_ice2 - beta_melt*ice_H1
      do j = 1,xdim
         do i = 1,ydim
            call ice_regression_coef(ice_T1(j,i,:), ice_Tcoef(j,i,:))
         end do
      end do
     
      ! ice temperature tendency
      do k = 1,4
          ! varible initialization 
          thk_layer = 0.; 
          ice_vx = 0.; ice_vy = 0.; ice_vxvm = 0.; ice_vyvm = 0.
          ! diagnostic variable in ice sheet model, including strata/vertical mean ice horizontal velocity, strain rate heating
          call ice_sheet_diagnostics(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, &
          &                          ice_vxvm, ice_vyvm, zeta(k), sigma_x(:,:,k), sigma_y(:,:,k))
          ! ice shelf velocity
          if(k==4) call ice_velocity_ssa(ice_zs,ice_H1, dyy, dxlat, ice_vx, ice_vy, ice_vxvm, ice_vyvm)
          ! store ice velocity variables
          ice_vx3(2:xdim,:,k) = (ice_vx(1:xdim-1,1:ydim)+ice_vx(2:xdim,1:ydim))/2.;
          ice_vx3(1,:,k) = (ice_vx(xdim,1:ydim)+ice_vx(1,1:ydim))/2.
          ice_vy3(:,:,k) = (ice_vy(:,1:ydim)+ice_vy(:,2:ydim+1))/2.
          ice_vxm0 = ice_vx3; ice_vym0 = ice_vy3
          ! ice thickness advection
          thk_layer = ice_H1*(zeta(k)+1)/2.
          call ice_advection(ice_vxvm, ice_vyvm, thk_layer, divg(:,:,k), dyy, dxlat, lat)
          ! advection term for ice thickness equation
          if(k==4) term_hadv = divg(:,:,k)
          ! vertical velocity
          call ice_advection_upwind(ice_vx, ice_vy, thk_layer, thk_layer, term_advx, term_advy, &
          &                         dyy, dxlat, lat, ice_H1)
          thk_adv = - (term_advx + term_advy)
          ice_vz(:,:,k) = (divg(:,:,k) - thk_adv)/dt_ice
          if(k > 1) ice_vzg(:,:,k) = (ice_vz(:,:,k-1)+ice_vz(:,:,k))/2.
          
          ! strain rate heating (stress from)
          ! term_sig(:,:,k) = dt_ice*term_sig(:,:,k)/(cp_ice*rho_ice)
          ! ice temperature advection
          call ice_advection_upwind(ice_vx, ice_vy, ice_T1(:,:,k), ice_T1(:,:,k), term_advx, term_advy, &
          &                         dyy, dxlat, lat, ice_H1)
          term_tadv(:,:,k) = - (term_advx + term_advy)
      end do
      ! vertical velocity shear for strain rate calculation
      where(ice_H1> d_ice_mov) dvxdz(:,:,1) = (ice_vx3(:,:,2) - ice_vx3(:,:,1))/(zeta(2)-zeta(1))*2./ice_H1
      where(ice_H1> d_ice_mov) dvydz(:,:,1) = (ice_vy3(:,:,2) - ice_vy3(:,:,1))/(zeta(2)-zeta(1))*2./ice_H1
      do k = 2,kdim-1
          where(ice_H1> d_ice_mov) dvxdz(:,:,k) = (ice_vx3(:,:,k+1) - ice_vx3(:,:,k-1))/(zeta(k+1)-zeta(k-1))*2./ice_H1
          where(ice_H1> d_ice_mov) dvydz(:,:,k) = (ice_vy3(:,:,k+1) - ice_vy3(:,:,k-1))/(zeta(k+1)-zeta(k-1))*2./ice_H1
      end do
      ! strain rate heating (velocity form)
      term_sig = dt_ice*(sigma_x*dvxdz+sigma_y*dvydz)/(cp_ice*rho_ice)
      ! vertical diffusion and convection
      ice_Td = ice_T1 + term_sig + term_tadv
      do j = 1,xdim
         do i = 1,ydim
            call ice_temperature_convection_diffusion(ice_Td(j,i,:), ice_H1(j,i), ice_vz(j,i,:), ice_vzg(j,i,:), &
            &                                          zeta, term_dif(j,i,:))
         end do
      end do
      ice_Tsm  = 0.
      ! calving, ice advected to ocean grid 
      where(mask < 0 ) term_calv = - term_hadv

  end if
!+++++++++++++++++++++++++++++ dynamic core  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! mass balance
  call ice_mass_balance(ice_H1, ice_Ts1, Ta1, dq_rain, Fn_surf, dT_ocean, dice_Ts, term_mass, ice_fus)

  ice_Ts0 = ice_Ts1
  
  ! ice temperature equation
  dice_T = term_tadv + term_dif + term_sig 
  do k = 1,4
      ice_Tmt(:,:,k) = Tl_ice2 - beta_melt*(1-zeta(k))*ice_H1/2.
  end do
  where(dice_T > ice_Tmt-ice_T1) dice_T = 0.9*(ice_Tmt-ice_T1) ! numeric stability
  ice_T0 = ice_T1 + dice_T
  ! issue
  ice_T1  = ice_T0

  call ice_topography(ice_zs, ice_H1) ! topography parameter change 

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_mass_balance(ice_H1, ice_Ts1, Ta1, dq_rain, Fn_surf, dT_ocean, dice_Ts, term_mass, ice_fus)
!+++++++++++++++++++++++++++++++++++++++
! mass balance: ice accumulation and ice fusion
! output: ice mass balance (term_mass); ice fusion latent heat (ice_fus)
  USE mo_numerics, ONLY: xdim, ydim, dt
  USE mo_physics,  ONLY: mask, z_topo, rho_ice, cp_ice, d_ice_max, cap_surf,c_lapse
  implicit none
  real*8, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature (current) [K]
  real*8, dimension(xdim,ydim) :: Fn_surf                 ! surface net heat flux [J/m2]
  real,   dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]

  real*8, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (current) [m]
  real*8, dimension(xdim,ydim) :: ice_snf                 ! snowfall accumulation rate [m/s]
  real*8, dimension(xdim,ydim) :: ice_melt                ! ice melting rate (negative -> melting) [m/s]
  real,   dimension(xdim,ydim) :: ice_fus                 ! ice melting latent heat [W/m2]

  real,   dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real*8, dimension(xdim,ydim) :: Ta_scl                  ! in situ air temperature based on the topography [K]
  real,   dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]

  real*8                       :: Tmin_limit              ! no very low Tsurf/Tatmoss;  numerical stability
  real*8, dimension(xdim,ydim) :: dice_Ts                 ! ice temperature total tendency [K]
  real*8, dimension(xdim,ydim) :: term_mass               ! ice thickness tendency due to mass balance [m] 

  dice_Ts = 0.; term_mass = 0.
  Ta_scl = Ta1 + c_lapse*z_topo
  ! ice accumulation, thickness increase by snowfall
  call ice_accumulation(ice_snf, ice_Ts1, Ta_scl, dq_rain) 
  ! ice fusion
  call ice_fusion( ice_fus, ice_melt, ice_Ts1, ice_H1, Fn_surf, dT_ocean) 
  
  ! mass balance 
  where(mask > 0.) term_mass = dt* (ice_snf + ice_melt)
  
  ! ice surface temperature
  dice_Ts  = dt*(Fn_surf + ice_fus ) / cap_surf 

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_accumulation(ice_snf, ice_Ts1, Ta1, dq_rain) 
!+++++++++++++++++++++++++++++++++++++++
! ice sheet : accumulation process
! output: snowfall rate (ice_snf)
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics,  ONLY: r_qviwv, Tl_ice1, Tl_ice2, ice_svlm, wz_vapor
  implicit none
  real*8                       :: T_range                 ! temperature range for partial snowfall [K]
  real*8, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature [K] 
  real*8, dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real*8, dimension(xdim,ydim) :: ice_snf                 ! snowfall accumulation rate [kg/m2/s]
  real*8, dimension(xdim,ydim) :: snf_rate                ! snowfall percentage of precipitation 
  real,   dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]
  ice_snf = 0.; snf_rate = 0.
  ! snowfall rate
  T_range = 2.
  where(ice_Ts1 - Tl_ice2 <= -T_range) snf_rate = 1.
  ! partial snowfall drop when temperature between -2 to 2 ˚C
  where(abs(ice_Ts1- Tl_ice2) < T_range) snf_rate = (-(ice_Ts1-Tl_ice2)/T_range+1)/2.
  ! ice sheet accumulation from snowfall, kg/m2/s   
  where((snf_rate > 0.) .and. (Ta1 <= Tl_ice2)) ice_snf = - snf_rate*ice_svlm*dq_rain*r_qviwv*wz_vapor; 
end subroutine ice_accumulation

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_fusion( ice_fus, ice_melt, ice_Ts1, ice_H1, Fn_surf, dT_ocean) 
!+++++++++++++++++++++++++++++++++++++++
! ice fusion: ice melting latent heat and ice thickness change due to melting
! output: ice thickness budget melting (ice_melt), ice fusion latent heat (ice_fus)
  USE mo_numerics, ONLY: xdim, ydim, dt
  USE mo_physics,  ONLY: ci_latent, cp_ice, Tl_ice2, cap_surf, cap_land,    &
&						 rho_ice, d_ice_max, mask, glacier_type, beta_melt, rho_ocean
  implicit none
  real*8, dimension(xdim,ydim) :: ice_Tse                 ! ice surface temperature estimate [K]
  real*8, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature [K]
  real*8, dimension(xdim,ydim) :: ice_H1                  ! ice thickness [m]
  real,   dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real*8, dimension(xdim,ydim) :: ice_melt                ! ice melting rate [m/s]
  real,   dimension(xdim,ydim) :: ice_fus                 ! ice melting latent heat [W/m2]
  real*8, dimension(xdim,ydim) :: heat_tsurf              ! heat flux power for heating surface temperature above frozen point [W/m2]
  real*8, dimension(xdim,ydim) :: hmax_melt               ! heat flux power maximum for total snow melts [W/m2] 
  real*8, dimension(xdim,ydim) :: Fn_surf                 ! surface heat flux without fusion [W/m2] 
  
  ice_melt  = 0.; ice_fus   = 0.
  heat_tsurf= 0.; hmax_melt = 0.
 
  ! snow heat capacity
!  where((ice_H1 <  d_ice_max).and.(ice_H1 > 0.)) cap_ice = ice_H1*cp_ice*rho_ice ! heat capacity of snow [J/K/m^2]
!  where(ice_H1 >= d_ice_max) cap_ice = d_ice_max*cp_ice*rho_ice ! heat capacity limitation of snow [J/K/m^2]
  
   ice_Tse  = ice_Ts1 + dT_ocean + dt*Fn_surf / cap_surf 
   where( mask > 0 .and. ice_H1 > 0.)
       ! the extra heat flux after the ice reaches melting point
       heat_tsurf = (ice_Tse - Tl_ice2) * cap_surf / dt
       ! the heat flux needed to melt whole ice column
       hmax_melt  = ci_latent * rho_ice * ice_H1 / dt
   end where
   ! surface snow totally melts away
   where( (mask > 0) .and. (heat_tsurf > 0 ) .and. (heat_tsurf >  hmax_melt) .and. (ice_H1 >0.))
         ice_melt    = - ice_H1 / dt
         ice_fus     = - hmax_melt 
   end where
   ! surface snow partially melts
   where( (mask > 0) .and. (heat_tsurf > 0 ) .and. (heat_tsurf <= hmax_melt) .and. (ice_H1 >0.))
         ice_melt    = - heat_tsurf / (rho_ice*ci_latent)
         ice_fus     = - heat_tsurf 
   end where

  where(ice_fus > 0.) ice_fus = 0.

end subroutine ice_fusion

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_sheet_diagnostics(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, &
                               & ice_vxm, ice_vym, zeta, sigma_x, sigma_y)
!+++++++++++++++++++++++++++++++++++++++
  ! diagnostic variable in ice sheet model (SIA), including ice velocity, strain rate heating
  ! output: ice strata velocity (ice_vx, ice_vy), ice vertical mean velocity (ice_vxm, ice_vym)
  !         strain rate tensor component (sigmax, sigma_y)

  USE mo_numerics, ONLY: xdim, ydim, dlat
  USE mo_physics, ONLY: pi, rho_ice, grav, beta_melt, Tl_ice2, b_rock, glacier_type, log_decouple, &
  &                     d_ice_mov, slid_law, ice_shlf
  ! NOTE: variable at point (j-1,i),(j,i-1), (j,i) is represented by _jm1, _im1 and _c respectively
  implicit none

  integer, external      :: cycle_ind                      ! cyclic zonal index 
  integer                :: i, j, ind_zonal_real           ! work variable for count
  integer                :: k1, k2, k3                     ! work variable for count
  real*8                 :: zeta                           ! z-direction coordinate variable (for base layer correction)
  real*8                 :: zh                             ! percentage of ice above the selected layer
  real*8                 :: zr                             ! percentage of ice below the selected layer

  real*8                 :: T_layer                        ! layer temperature  [K]
  real*8                 :: iceH_c,iceH_jm1, iceH_im1      ! ice thickness      [m]  
  real*8                 :: iceH_jp1, iceH_ip1             ! ice thickness      [m]  
  real*8                 :: T_gs3, T_gsjm3, T_gsim3        ! temperature at Gauss point for Gaussian intergration [K]
  real*8                 :: T_gs4, T_gsjm4, T_gsim4        ! temperature at Gauss point for Gaussian intergration [K]
  real*8                 :: gs_zeta3, gs_zeta4             ! Gaussian-Jacob quadrature node for intergration (zeta coordinate) [1] 
  real*8                 :: Tm_layer, Tmt, Tmt_jm, Tmt_im  ! pressure melting temperature [K]
  real*8                 :: iceZ_jm1, iceZ_im1, iceZ_c     ! ice surface height [m]
  real*8                 :: iceZ_jp1, iceZ_ip1             ! ice surface height [m]
  real*8                 :: glc_jm1, glc_im1, glc_c        ! glacier type (-1 for ocean, 0 for floating ice, 1 for grounded ice) 
  real*8                 :: sighrz, sighrzjm, sighrzim     ! rho*g*|| del h ||  [kg/(m2 s2)]
  real*8                 :: delh                           ! || del h ||  [1]
  real*8                 :: coef, coef_jm, coef_im         ! regression coefficient for strata temperature [K] 
  real*8                 :: term_poly3, term_polyjm3, term_polyim3 ! polynomial coefficient and integration parameter [m4]
  real*8                 :: term_poly4, term_polyjm4, term_polyim4 ! polynomial coefficient and integration parameter [m4]
  real*8                 :: cons_int3, cons_intjm3, cons_intim3    ! integration result for strata velocity calculation [(m3 s)/kg]
  real*8                 :: cons_int4, cons_intjm4, cons_intim4    ! integration result for mean velocity calculation [(m3 s)/kg]
  real*8                 :: dZdx, dZdy, dZ                 ! ice surface height gradient
  real*8                 :: coef_vgrd                      ! coefficient for velocity gradient calculation [m/s]
  real*8                 :: dvxdz, dvydz                   ! horizontal velocity vertical gradient [s-1]
  real*8                 :: vxb, vyb                       ! horizontal velocity at basal layer [m/s]
  real*8                 :: shlf_vel                       ! shelf velocity [m/s] 

  real*8, dimension(ydim)      :: lat                      ! latitude [deg]
  real*8, dimension(ydim)      :: dxlat                    ! zonal grid width [m]
  real*8                       :: dyy                      ! meridional grid width [m]
  integer,dimension(ydim)      :: ilat = (/(i,i=1,ydim)/)  ! index for latitude grid

  ! coefficient for Gaussian-Jacob integration: gs_node for node, gs_weight for weight, *m for whole column integration, int(f(x)(1+x)**alpha*(1-x)**beta)
  real*8,dimension(4)   :: coef_poly  = (/1,3,3,1/)       
  ! for column mean velocity, alpha = 0; beta = 1~4
  real*8,dimension(5,4) :: gs_node4   = reshape((/-0.920380285897062,-0.603973164252784,-0.124050379505228, &
                                                   0.390928546707272,0.802929828402347,&
                                              -0.930842120163570,-0.653039358456609,-0.220227225868961, &
                                              0.268666945261774,0.702108425894033,&
                                              -0.938871406986534,-0.691210299947676,-0.297140008474123, &
                                              0.165799576355599,0.607575985206580,&
                                              -0.945228837068531,-0.721772687352885,-0.360118465482019, &
                                              0.0781326509369691,0.520415910395038/),(/5,4/))
  real*8,dimension(5,4) :: gs_weight4 = reshape((/0.387126360906607,0.668698552377478,0.585547948338679, &
                                                  0.295635480290467,0.0629916580867690,&
                                               0.654118274286167,1.00959169519929,0.713601289772721, &
                                               0.256444805783696,0.0329106016247922,&
                                               1.13216844954808,1.60506407143933,0.967119251259600, &
                                               0.271317197812349,0.0243310299406412, &
                                               1.99536955485696,2.64680681926207,1.40878302417225, &
                                               0.326737335215111,0.0223032664936029/),(/5,4/))
  ! for strata velocity, alpha = 0, beta = 0~3
  real*8,dimension(5,4) :: gs_node3   = reshape((/-0.906179845938664,-0.538469310105683,0., &
                                                   0.538469310105683,0.906179845938664,&
                                                  -0.920380285897062,-0.603973164252784,-0.124050379505228, &
                                                   0.390928546707272,0.802929828402347,&
                                                  -0.930842120163570,-0.653039358456609,-0.220227225868961, &
                                                   0.268666945261774,0.702108425894033,&
                                                  -0.938871406986534,-0.691210299947676,-0.297140008474123, &
                                                   0.165799576355599,0.607575985206580/),(/5,4/))
  real*8,dimension(5,4) :: gs_weight3 = reshape((/0.236926885056189,0.478628670499367,0.568888888888889, &
                                                  0.478628670499367,0.236926885056189,&
                                                  0.387126360906607,0.668698552377478,0.585547948338679, &
                                                  0.295635480290467,0.0629916580867690,&
                                                  0.654118274286167,1.00959169519929,0.713601289772721, &
                                                  0.256444805783696,0.0329106016247922,&
                                                  1.13216844954808,1.60506407143933,0.967119251259600, &
                                                  0.271317197812349,0.0243310299406412/),(/5,4/))
  real*8, dimension(xdim,ydim) :: ice_H0                  ! ice thickness (forward) [m]
  real*8, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (current) [m]
  real*8, dimension(xdim,ydim) :: ice_zs                  ! ice surface height [m]
  real*8, dimension(xdim,ydim) :: sigma_x                 ! x-z component of stress tensor [Pa]
  real*8, dimension(xdim,ydim) :: sigma_y                 ! y-z component of stress tensor [Pa]
  real*8, dimension(xdim,ydim) :: adv_cord                ! advection term due to coordinate transformation (dvx/dz dh/dx + dvy/dz dh/dy) [s-1] 

  real*8, dimension(xdim,ydim+1)   :: ice_vx, ice_vy      ! ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim+1)   :: ice_vxm, ice_vym    ! vertical mean of ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim,4)   :: ice_Tcoef           ! regression coefficient of ice layer temperature [k] 
  
  sigma_x = 0.
  sigma_y = 0.  

  lat = dlat*ilat-dlat/2.-90.;  
  zh = (zeta+1)/2
  zr = 1-zh
  
  do i = 1,ydim+1
      do j = 1,xdim
          call meridion_shift(ice_H1, j, i, iceH_c)
          call meridion_shift(ice_H1, j, i+1, iceH_ip1)
          call meridion_shift(ice_H1, j, i-1, iceH_im1)
          ind_zonal_real = cycle_ind(j+1,xdim) 
          call meridion_shift(ice_H1, ind_zonal_real, i, iceH_jp1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_H1, ind_zonal_real, i, iceH_jm1)
          if(iceH_c <= d_ice_mov .and. iceH_im1 <= d_ice_mov .and. iceH_jm1 <= d_ice_mov) cycle

          ! effective stress horizontal part, vertical part is left for intergal
          sighrz      = sigma_horizon(ice_zs, j  , i, dxlat, dyy)
          sighrzjm    = sigma_horizon(ice_zs, j-1, i, dxlat, dyy)
          sighrzim    = sigma_horizon(ice_zs, j, i-1, dxlat, dyy)
          ! vertical intergal part
          cons_int3 = 0.; cons_intjm3 = 0.; cons_intim3 = 0.
          cons_int4 = 0.; cons_intjm4 = 0.; cons_intim4 = 0.
          vxb = 0.; vyb = 0.
          do k1 = 1,4
              do k2 = 1,5
                  ! temperature at Gauss node 3 for strata velocity calculation, 4 for mean velocity calculation
                  T_gs3 = 0.; T_gsjm3 = 0.; T_gsim3 = 0.; term_poly3 = 0.
                  T_gs4 = 0.; T_gsjm4 = 0.; T_gsim4 = 0.; term_poly4 = 0.
                  gs_zeta3 = (gs_node3(k2,k1)+1)*zh - 1. 
                  gs_zeta4 = (gs_node4(k2,k1)+1)*zh - 1. 
                  do k3 = 1,4
                      call meridion_shift(ice_Tcoef(:,:,k3), j, i, coef  )
                      T_gs3   = T_gs3   + coef*gs_zeta3**(k3-1)
                      T_gs4   = T_gs4   + coef*gs_zeta4**(k3-1)

                      ind_zonal_real = cycle_ind(j-1,xdim) 
                      call meridion_shift(ice_Tcoef(:,:,k3), ind_zonal_real, i, coef_jm)
                      T_gsjm3 = T_gsjm3 + coef_jm*gs_zeta3**(k3-1)
                      T_gsjm4 = T_gsjm4 + coef_jm*gs_zeta4**(k3-1)

                      call meridion_shift(ice_Tcoef(:,:,k3), j, i-1, coef_im)
                      T_gsim3 = T_gsim3 + coef_im*gs_zeta3**(k3-1)
                      T_gsim4 = T_gsim4 + coef_im*gs_zeta4**(k3-1)
                  end do

                  ! pressure melting temperature
                  Tmt    = Tl_ice2 - beta_melt*(1-gs_zeta3)*iceH_c  /2.
                  Tmt_im = Tl_ice2 - beta_melt*(1-gs_zeta3)*iceH_im1/2.
                  Tmt_jm = Tl_ice2 - beta_melt*(1-gs_zeta3)*iceH_jm1/2.
                  if(T_gs3   > Tmt     ) T_gs3   = Tmt
                  if(T_gsim3 > Tmt_im  ) T_gsim3 = Tmt_im
                  if(T_gsjm3 > Tmt_jm  ) T_gsjm3 = Tmt_jm

                  ! polynomial coefficient
                  term_poly3 = 0.; term_polyjm3 = 0.; term_polyim3 = 0.
                  if(iceH_c   /= 0) term_poly3    = (iceH_c  /2.)**4*gs_weight3(k2,k1) &
                  &                               * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  if(iceH_jm1 /= 0) term_polyjm3  = (iceH_jm1/2.)**4*gs_weight3(k2,k1) &
                  &                               * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  if(iceH_im1 /= 0) term_polyim3  = (iceH_im1/2.)**4*gs_weight3(k2,k1) &
                  &                               * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  
                  ! Gauss-Jacobi quadrature
                  cons_int3    = cons_int3    + constitutive_equation(T_gs3  , Tmt   )*term_poly3  *sighrz  **2
                  cons_intjm3  = cons_intjm3  + constitutive_equation(T_gsjm3, Tmt_jm)*term_polyjm3*sighrzjm**2
                  cons_intim3  = cons_intim3  + constitutive_equation(T_gsim3, Tmt_im)*term_polyim3*sighrzim**2
                  
                  ! pressure melting temperature
                  Tmt    = Tl_ice2 - beta_melt*(1-gs_zeta4)*iceH_c  /2.
                  Tmt_im = Tl_ice2 - beta_melt*(1-gs_zeta4)*iceH_im1/2.
                  Tmt_jm = Tl_ice2 - beta_melt*(1-gs_zeta4)*iceH_jm1/2.
                  if(T_gs4   > Tmt     ) T_gs4   = Tmt
                  if(T_gsim4 > Tmt_im  ) T_gsim4 = Tmt_im
                  if(T_gsjm4 > Tmt_jm  ) T_gsjm4 = Tmt_jm
                  
                  ! polynomial coefficient
                  term_poly4 = 0.; term_polyjm4 = 0.; term_polyim4 = 0.
                  if(iceH_c   /= 0) term_poly4   = (iceH_c  /2.)**4*gs_weight4(k2,k1)/2. &
                  &                              * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  if(iceH_jm1 /= 0) term_polyjm4 = (iceH_jm1/2.)**4*gs_weight4(k2,k1)/2. &
                  &                              * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  if(iceH_im1 /= 0) term_polyim4 = (iceH_im1/2.)**4*gs_weight4(k2,k1)/2. &
                  &                              * coef_poly(k1)*(2.*zr)**(4-k1)*zh**k1
                  
                  ! Gauss-Jacobi quadrature
                  cons_int4   = cons_int4   + constitutive_equation(T_gs4  , Tmt   )*term_poly4  *sighrz  **2
                  cons_intjm4 = cons_intjm4 + constitutive_equation(T_gsjm4, Tmt_jm)*term_polyjm4*sighrzjm**2
                  cons_intim4 = cons_intim4 + constitutive_equation(T_gsim4, Tmt_im)*term_polyim4*sighrzim**2
              end do
          end do

          ! layer temperature
          T_layer = 0.
          do k1 =1,4 
              call meridion_shift(ice_Tcoef(:,:,k1), j, i, coef  )
              T_layer   = T_layer   + coef*zeta**(k1-1)
          end do
          Tm_layer  = Tl_ice2 - beta_melt*(1-zeta)*iceH_c/2.
          if(T_layer   > Tm_layer) then 
              T_layer   = Tm_layer
          end if

          ! ice surface gradient
          call meridion_shift(ice_zs, j, i  , iceZ_c)
          call meridion_shift(ice_zs, j, i+1, iceZ_ip1)
          call meridion_shift(ice_zs, j, i-1, iceZ_im1)
          ind_zonal_real = cycle_ind(j+1,xdim) 
          call meridion_shift(ice_zs, ind_zonal_real, i, iceZ_jp1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_zs, ind_zonal_real, i, iceZ_jm1)
          if (i>ydim) then
              dZdx = 0.
          else
              dZ = (iceZ_c - iceZ_jm1)
              dZdx = dZ/dxlat(i)
          end if
          dZ = (iceZ_c - iceZ_im1)
          dZdy = dZ/dyy
         
          ! basal velcocity according to ice
          delh = (sighrz+sighrzim)/(2*rho_ice*grav)
          vyb = - slid_law*dZdy*max(iceH_c, iceH_im1)*delh**2
          delh = (sighrz+sighrzjm)/(2*rho_ice*grav)
          vxb = - slid_law*dZdx*max(iceH_c, iceH_jm1)*delh**2

          ice_vx(j,i)   = -2.*rho_ice*grav*dZdx*(cons_int3+cons_intjm3)/2. + vxb
          ice_vy(j,i)   = -2.*rho_ice*grav*dZdy*(cons_int3+cons_intim3)/2. + vyb
          ice_vxm(j,i)  = -2.*rho_ice*grav*dZdx*(cons_int4+cons_intjm4)/2. + vxb
          ice_vym(j,i)  = -2.*rho_ice*grav*dZdy*(cons_int4+cons_intim4)/2. + vyb

          if(i<ydim+1) then
              ! layer height gradient 
              dZ    = (iceZ_jp1 - iceZ_jm1)
              dZdx  = dZ/2./dxlat(i)
              dZ    = (iceZ_ip1 - iceZ_im1)
              dZdy  = dZ/2./dyy
              ! stress tensor component
              sigma_x(j,i)  = -rho_ice*grav*zr*iceH_c*dZdx
              sigma_y(j,i)  = -rho_ice*grav*zr*iceH_c*dZdy
          end if
      end do
  end do

  contains 
       real function constitutive_equation(T_Gauss, T_melt)
       ! constitutive equation based on Glen's flow law
       ! output: consititutive function (SIA)
       USE mo_physics, ONLY: A_fact, actene, R_gas, ice_Tcons, enh_fact, TL_ice2
       implicit none
       real*8               :: T_Gauss
       real*8               :: T_prime
       real*8               :: T_melt
       real*8               :: c
       T_prime = T_Gauss - T_melt + Tl_ice2
       if(T_prime > ice_Tcons) then
           constitutive_equation =  A_fact(1)*exp(-actene(1)/R_gas/T_prime)*enh_fact
       else
           constitutive_equation =  A_fact(2)*exp(-actene(2)/R_gas/T_prime)*enh_fact
       end if
       
       if(log_decouple==1) then
           constitutive_equation = 1d-16/86400/365 ! from Huybrechts et al. (1996)
       end if
       end function 

       real function sigma_horizon(ice_zs, ind_zonal, ind_merid, dxlat, dy)
       ! output: effective strain rate horizontal part: rho g |grad h|
       USE mo_numerics, ONLY: xdim, ydim
       USE mo_physics, ONLY: rho_ice, grav
       implicit none
       integer, external  :: cycle_ind
       integer            :: ind_zonal, ind_merid, ind_zonal_real
       real*8             :: iceZ_jp1, iceZ_jm1, iceZ_ip1, iceZ_im1, iceZ_c   ! surface elevation at different grid [m]
       real*8             :: dZdx_c, dZdy_c, dZ                               ! elevation gradient
       real*8             :: dx, dy                                           ! spatial interval [m]
       real*8, dimension(ydim)      :: dxlat
       real*8, dimension(xdim,ydim) :: ice_zs
       ! zonal
       ind_zonal_real = cycle_ind(ind_zonal-1,xdim)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid, iceZ_jm1)
       ind_zonal_real = cycle_ind(ind_zonal+1,xdim)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid, iceZ_jp1)
       ind_zonal_real = cycle_ind(ind_zonal,xdim)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid  , iceZ_c  )
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid-1, iceZ_im1)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid+1, iceZ_ip1)
       if   (ind_merid.gt.ydim) then
            dx = dxlat(2*ydim+1-ind_merid)
       else if(ind_merid.lt.1 ) then
            dx = dxlat(1-ind_merid)
       else 
            dx = dxlat(ind_merid)
       end if
       ! central difference version
       dZ = abs(iceZ_jp1-iceZ_c) + abs(iceZ_c-iceZ_jm1)
       dZdx_c = dZ/(2.*dx)
       dZ = abs(iceZ_ip1-iceZ_c) + abs(iceZ_c-iceZ_im1)
       dZdy_c = dZ/(2.*dy)
       ! average version
       sigma_horizon  = rho_ice*grav*sqrt(dZdx_c**2+dZdy_c**2)
       end function 

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_velocity_ssa(ice_zs,ice_H1, dyy, dxlat, ice_vx, ice_vy, ice_vxm, ice_vym)
!+++++++++++++++++++++++++++++++++++++++
  ! ice velocity calculation for SSA
  ! output: ice strata velocity (ice_vx, ice_vy), ice vertical mean velocity (ice_vxm, ice_vym)

  USE mo_numerics, ONLY: xdim, ydim, dlat, dlon
  USE mo_physics, ONLY: pi, rho_ice, rho_ocean, grav, glacier_type, ice_shlf, b_rock, shlf_eta
  ! NOTE: variable at point (j-1,i),(j,i-1), (j,i) is represented by _jm1, _im1 and _c respectively
  implicit none
  integer, external               :: cycle_ind
  integer                         :: i, j, ind_zonal_real           ! work variable for count
  integer                         :: it                             ! work variable for count
  integer                         :: itmax                          ! iteration maximum
  real*8                          :: omega                          ! accelerate factor for iteration solver core
  real*8                          :: glc_jm1, glc_im1, glc_c        ! glacier type (-1 for ocean, 0 for floating ice, 1 for grounded ice) 
  real*8                          :: glc_jp1, glc_ip1               ! glacier type (-1 for ocean, 0 for floating ice, 1 for grounded ice) 
  real*8                          :: iceZ_jm1, iceZ_im1, iceZ_c     ! surface elevation in different grid [m] 
  real*8                          :: vx_jm1, vx_jp1, vx_c, vx_ip1   ! zonal velocity value in different grid [m/s] 
  real*8                          :: vx_jp1im1, vx_jp1ip1, vx_im1   ! zonal velocity value in different grid [m/s]
  real*8                          :: vy_im1, vy_ip1, vy_c, vy_jp1   ! meridian velocity value in different grid [m/s] 
  real*8                          :: vy_jm1ip1, vy_jp1ip1, vy_jm1   ! meridian velocity value in different grid [m/s]
  real*8                          :: wlat_im1, wlat_ip1, wlat_c     ! area weight due to latitude change
  real*8                          :: res                            ! velocity residuels [m/s]
  real*8                          :: termZx, termZy                 ! surface gradient term  
  real*8                          :: bfact                          ! boundary factor for grounded ice interface. Refence: Winkelmann et al. 2011 (PISM-PIK)
  real*8                          :: deg                            ! degree to radial
  real*8                          :: dyy                            ! meridional grid width [m]
  real*8                          :: dxx                            ! zonal grid width [m]
  real*8                          :: ddy2                           ! differential operator [m -2]
  real*8                          :: vis_c                          ! viscosity [Pa s]
  real*8                          :: iceH_c,iceH_jm1im1,iceH_jmhimh ! ice thickness in different grid [m] 
  real*8                          :: iceH_jphimh, iceH_jmhiph       ! ice thickness in different grid [m]  
  real*8                          :: iceH_jm1ip1, iceH_ip1, iceH_jm1! ice thickness in different grid [m] 
  real*8                          :: iceH_jp1im1, iceH_jp1, iceH_im1! ice thickness in different grid [m] 
  integer, dimension(ydim)        :: ilat = (/(i,i=1,ydim)/)        ! index for latitude grid
  real*8, dimension(ydim)         :: lat                            ! latitude [deg]
  real*8, dimension(ydim)         :: dxlat                          ! zonal grid width [m]
  real*8, dimension(ydim)         :: wlat                           ! latitude weight
  real*8, dimension(ydim+1)       :: ddx2, ddxdy                    ! differential operator [m -2]
  real*8, dimension(xdim,ydim)    :: ice_zs                         ! ice surface height [m]
  real*8, dimension(xdim,ydim)    :: ice_H1                         ! ice thickness (current) [m]
  real*8, dimension(xdim,ydim+1)  :: ice_vx, ice_vy                 ! ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim+1)  :: ice_vxm, ice_vym               ! vertical mean of ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim+1)  :: vx_ssa, vy_ssa                 ! Shallow Shelf Approximated ice velocity at zonal (x) and meridian (y) direction [m/s]
  real*8, dimension(xdim,ydim+1)  :: vx0_ssa, vy0_ssa               ! Shallow Shelf Approximated ice velocity at zonal (x) and meridian (y) direction [m/s]
  integer,dimension(xdim,ydim+1)  :: shlf_mark                      ! ice shelf mark
  
  deg  = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  lat  = dlat*ilat-dlat/2.-90.; 
  wlat(1:ydim) = cos(2.*pi/360.*lat);
  wlat         = 1.

  vx_ssa   = ice_vxm; vy_ssa   = ice_vym; shlf_mark = -1
  vx0_ssa  = vx_ssa;  vy0_ssa  = vy_ssa 
  itmax  = 200
  ! matrix coefficient
  ddx2(1:ydim) = 1/dxlat**2; ddx2(ydim+1) = 1/dxlat(ydim)**2; ddy2 = 1/dyy**2; 
  ddxdy(1:ydim)= 1/(dxlat*dyy); ddxdy(ydim+1)= 1/(dxlat(ydim)*dyy)

  do it = 1,itmax 
      do i = 1,ydim+1
          do j = 1,xdim
              call meridion_shift(glacier_type, j, i  , glc_c)
              call meridion_shift(glacier_type, j, i-1, glc_im1)
              ind_zonal_real = cycle_ind(j-1,xdim) 
              call meridion_shift(glacier_type, ind_zonal_real, i, glc_jm1)
              if(glc_c /= 0 .and. glc_jm1 /= 0 .and. glc_im1 /= 0) cycle
              shlf_mark(j,i) = 1
              ! area weight for deravative
              if(i >= 2 .and. i <= ydim-1) then
                  wlat_c = wlat(i); wlat_im1 = wlat(i-1); wlat_ip1 = wlat(i+1)
              end if
              if( i > ydim-1) then
                  wlat_c = wlat(2*ydim-i); wlat_im1 = wlat(i-1); wlat_ip1 = wlat(2*ydim-i+1)
              end if
              if( i < 2) then
                  wlat_c = wlat(i); wlat_im1 = wlat(2-i); wlat_ip1 = wlat(i+1)
              end if
              ! velocity at neibor grids (zonal)
              call meridion_shift(vx_ssa, j             , i  , vx_c)
              call meridion_shift(vx_ssa, j             , i-1, vx_im1)
              call meridion_shift(vx_ssa, j             , i+1, vx_ip1)
              ind_zonal_real = cycle_ind(j-1,xdim) 
              call meridion_shift(vx_ssa, ind_zonal_real, i  , vx_jm1)
              ind_zonal_real = cycle_ind(j+1,xdim) 
              call meridion_shift(vx_ssa, ind_zonal_real, i  , vx_jp1)
              call meridion_shift(vx_ssa, ind_zonal_real, i-1, vx_jp1im1)
              call meridion_shift(vx_ssa, ind_zonal_real, i+1, vx_jp1ip1)
    
              ! velocity at neibor grids (meridian)
              call meridion_shift(vy_ssa, j             , i  , vy_c)
              call meridion_shift(vy_ssa, j             , i-1, vy_im1)
              call meridion_shift(vy_ssa, j             , i+1, vy_ip1)
              ind_zonal_real = cycle_ind(j-1,xdim) 
              call meridion_shift(vy_ssa, ind_zonal_real, i  , vy_jm1)
              call meridion_shift(vy_ssa, ind_zonal_real, i+1, vy_jm1ip1)
              ind_zonal_real = cycle_ind(j+1,xdim) 
              call meridion_shift(vy_ssa, ind_zonal_real, i  , vy_jp1)
              call meridion_shift(vy_ssa, ind_zonal_real, i+1, vy_jp1ip1)

              ! surface elevation 
              call meridion_shift(ice_zs, j             , i  , iceZ_c)
              call meridion_shift(ice_zs, j             , i-1, iceZ_im1)
              ind_zonal_real = cycle_ind(j-1,xdim) 
              call meridion_shift(ice_zs, ind_zonal_real, i  , iceZ_jm1)
             
              ! ice thickness 
              call meridion_shift(ice_H1, j             , i  , iceH_c)
              call meridion_shift(ice_H1, j             , i-1, iceH_im1)
              call meridion_shift(ice_H1, j             , i+1, iceH_ip1)
              ind_zonal_real = cycle_ind(j-1,xdim) 
              call meridion_shift(ice_H1, ind_zonal_real, i  , iceH_jm1)
              call meridion_shift(ice_H1, ind_zonal_real, i-1, iceH_jm1im1)
              call meridion_shift(ice_H1, ind_zonal_real, i+1, iceH_jm1ip1)
              ind_zonal_real = cycle_ind(j+1,xdim) 
              call meridion_shift(ice_H1, ind_zonal_real, i  , iceH_jp1)
              call meridion_shift(ice_H1, ind_zonal_real, i-1, iceH_jp1im1)

              iceH_jmhiph = (iceH_jm1+iceH_jm1ip1+iceH_c+iceH_ip1)/4.
              iceH_jmhimh = (iceH_jm1+iceH_jm1im1+iceH_c+iceH_im1)/4.
              iceH_jphimh = (iceH_jp1+iceH_jp1im1+iceH_c+iceH_im1)/4.
              
              vis_c = shlf_eta 

              termZx         = rho_ice*grav*(iceH_c+iceH_jm1)*(iceZ_c-iceZ_jm1)/(2*dxlat(i)*vis_c)
              termZy         = rho_ice*grav*(iceH_c+iceH_im1)*(iceZ_c-iceZ_im1)/(2*dyy*vis_c)
    
              ! Gauss–Seidel method for linear system solution
              if(glc_c == 0 .and. glc_jm1 == 0 ) then
                  vx_ssa(j,i) = (4*ddx2(i)*(vx_jp1*iceH_c + vx_jm1*iceH_jm1)         &
                  &           + ddy2      *(vx_ip1*iceH_jmhiph + vx_im1*iceH_jmhimh) &
                  &           + 2*ddxdy(i)*(vy_ip1*iceH_c + vy_jm1*iceH_jm1 - vy_c*iceH_c - vy_jm1ip1*iceH_jm1)&
                  &           + ddxdy(i)  *(vy_ip1*iceH_jmhiph + vy_jm1*iceH_jmhimh  &
                  &           - vy_c*iceH_jmhimh - vy_jm1ip1*iceH_jmhiph) - termZx)  &
                  &           / (4*ddx2(i)*(iceH_c+iceH_jm1)+ddy2*(iceH_jmhiph+iceH_jmhimh))
              end if
              if(glc_c == 0 .and. glc_im1 == 0 ) then
                  vy_ssa(j,i) = (4*ddy2   *(vy_ip1*iceH_c + vy_im1*iceH_im1)         &
                  &           + ddx2(i)   *(vy_jp1*iceH_jphimh + vy_jm1*iceH_jmhimh) &
                  &           + 2*ddxdy(i)*(vx_jp1*iceH_c + vx_im1*iceH_im1 - vx_c*iceH_c - vx_jp1im1*iceH_im1)&
                  &           + ddxdy(i)  *(vx_jp1*iceH_jphimh + vx_im1*iceH_jmhimh  &
                  &           - vx_c*iceH_jmhimh - vx_jp1im1*iceH_jphimh) - termZy)  &
                  &           / (4*ddy2   *(iceH_c+iceH_im1)+ddx2(i)*(iceH_jphimh+iceH_jmhimh))
              end if
              
!              ! boundary condition  
              termZx         = rho_ice*grav*max(iceH_c,iceH_jm1)/vis_c/2
              termZx         = sign(termZx, iceZ_c-iceZ_jm1)
              termZy         = rho_ice*grav*max(iceH_c,iceH_im1)/vis_c/2
              termZy         = sign(termZy, iceZ_c-iceZ_im1)
              
              if( (glc_c == 0 .or. glc_jm1 == 0) .and. (glc_c < 0 .or. glc_jm1 < 0) ) then
                  bfact       = (rho_ocean - rho_ice)/rho_ocean 
                  vx_ssa(j,i) = (4/dxlat(i)*vx_jp1 + 2/dyy*(vy_ip1-vy_c) - termZx*bfact)/(4/dxlat(i))
              end if
              if( (glc_c == 0 .or. glc_im1 == 0) .and. (glc_c < 0 .or. glc_im1 < 0) ) then
                  bfact       = (rho_ocean - rho_ice)/rho_ocean 
                  vy_ssa(j,i) = (4/dyy*vy_ip1 + 2/dxlat(i)*(vx_jp1-vx_c) - termZy*bfact)/(4/dyy)
              end if
              
          end do
      end do
      res = max(maxval(abs(vx_ssa-vx0_ssa)),maxval(abs(vy_ssa-vy0_ssa)))*86400*365
      if(res < 0.1) then
      print*,'iteration time for SSA veloticy: ',it
      exit
      end if
      omega = 0
      vx0_ssa = vx_ssa*(1-omega)+vx0_ssa*omega; vy0_ssa = vy_ssa*(1-omega)+vy0_ssa*omega
  end do

  where(abs(vx_ssa) > ice_shlf) vx_ssa = vx_ssa/abs(vx_ssa)*ice_shlf 
  where(abs(vy_ssa) > ice_shlf) vy_ssa = vy_ssa/abs(vy_ssa)*ice_shlf 

  where(shlf_mark == 1) 
      ice_vxm = vx_ssa 
      ice_vym = vy_ssa
      ice_vx  = vx_ssa
      ice_vy  = vy_ssa
  end where

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_regression_coef(T1, Tcoef)
!+++++++++++++++++++++++++++++++++++++++
! regression coefficient calcualtion for vertical profile;
! output: regression coefficient for vertical temperature profile
   implicit none
   integer             :: k1, k2
   real*8                :: sum_var
   real*8,dimension(4,4) :: rg_coef   = reshape((/-0.250, 0.750, 0.750, -0.250,&
                                             0.250,	-1.299038105676658,	1.299038105676658, -0.250,&
                                             0.750,	-0.750,	-0.750,	0.750,&
                                             -0.750, 1.299038105676658,	-1.299038105676658,	0.750/),(/4,4/))
   real*8, dimension(4)  :: T1, Tcoef
   do k1 = 1,4
       sum_var = 0
       do k2 = 1,4
            sum_var = sum_var + rg_coef(k2,k1)*T1(k2)
       end do
       Tcoef(k1) = sum_var
   end do
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_temperature_convection_diffusion(T_ini, H1_ini, w_layer,  w, zeta, dif)
!+++++++++++++++++++++++++++++++++++++++
! solution for implicit form of convection-diffusion problem
! output: temperature budget vertical diffusion term (dif)
  USE mo_numerics, ONLY: xdim, ydim, dt_ice
  USE mo_physics, ONLY: d_ice_mov, ice_kappa, cp_ice, rho_ice, geoh_flux, Tl_ice2,beta_melt
  implicit none
  integer, parameter   :: kdim     = 4                   ! vertical layer
 
  integer              :: num                            ! iteration number
  integer              :: k,i                            ! work variable for count
  integer              :: indi,indj                      ! work variable for count
  real*8                 :: dt                             ! iteration time step
  real*8                 :: H1_ini                         ! input ice thickness [m]
  real*8                 :: Tmt_base                       ! melting point at base [K]

  real*8                 :: wp1                            ! vertical velocity used for plus 1 grid [m/s]
  real*8                 :: wm1                            ! vertical velocity used for minus 1 grid [m/s]

  real*8,dimension(kdim) :: zeta                           ! z-direction coordinate variable 
  real*8,dimension(kdim) :: dzeta                          ! z grid width 
  real*8,dimension(kdim) :: w_layer                        ! vertical velocity at layer [m/s] 
  real*8,dimension(kdim) :: w                              ! vertical velocity at half grid [m/s] 

  real*8,dimension(kdim) :: coefp1                         ! coefficient for implicit scheme 
  real*8,dimension(kdim) :: coefct                         ! coefficient for implicit scheme  
  real*8,dimension(kdim) :: coefm1                         ! coefficient for implicit scheme  
  real*8,dimension(kdim) :: coef_f                         ! coefficient for implicit scheme  
  real*8,dimension(kdim) :: kdtdx                          ! kappa*dt/dx 

  real*8,dimension(kdim) :: T_ini                          ! input ice strata temperature [K]
  real*8,dimension(kdim) :: ice_T0                         ! ice strata temperature for iteration (forward) [k]
  real*8,dimension(kdim) :: ice_T1                         ! ice strata temperature for iteration (backward) [k]
  real*8,dimension(kdim) :: dif                            ! ice strata temperature tendency due to diffusion [K]
  
  do k=2,kdim
     dzeta(k)=zeta(k)-zeta(k-1)
  end do
  dzeta(1)=0
  !dzeta(1)=dzeta(2)

  dt=dt_ice; num=1
  
  do k=1,kdim-1
     kdtdx(k) = 2*ice_kappa*dt/(dzeta(k+1)+dzeta(k))*(2./H1_ini)/(cp_ice*rho_ice)
  end do

  
  ice_T1 = T_ini
  do i=1,num
      ! coefficient for matrix
      do k=2,kdim-1
          wp1=0.; wm1=0.
          if(w_layer(k)>0.) then
              wm1 = w(k)
          elseif(w_layer(k)<0.) then
              wp1 = w(k+1)
          end if
          coefp1(k) = -(kdtdx(k)/dzeta(k+1)-dt/dzeta(k+1)*wp1)*(2./H1_ini)
          coefct(k) = 1 + (kdtdx(k)*(1./dzeta(k+1)+1./dzeta(k)) - dt*(wp1/dzeta(k+1)-wm1/dzeta(k)))*(2./H1_ini)
          coefm1(k) = -(kdtdx(k)/dzeta(k)+dt/dzeta(k)*wm1)*(2./H1_ini)
          coef_f(k) = ice_T1(k) 
      end do
    
      k=kdim
      coefp1(k)=0.; coefct(k)=1.; coefm1(k)=0.; coef_f(k)=ice_T1(k)
      k=1
      wp1=0.; wm1=0.
      if(w_layer(k+1)>0.) then
          wm1 = w(k)
      elseif(w_layer(k+1)<0.) then
          wp1 = w(k+1)
      end if

      ! boundary condition I
      coefp1(k) = -(kdtdx(k)/dzeta(k+1)-dt/dzeta(k+1)*wp1)*(2./H1_ini)
      coefct(k) = 1 + (kdtdx(k)/dzeta(k+1) - dt*wp1/dzeta(k+1))*(2./H1_ini) 
      coefm1(k) = 0. 
      coef_f(k) = ice_T1(k) + (kdtdx(k)+dt*wm1)*geoh_flux/ice_kappa 

      ! boundary condition II
      !coefp1(k) = -(kdtdx(k)*(1./dzeta(k+1)+1./dzeta(k)) - dt*(wp1/dzeta(k+1)-wm1/dzeta(k)))*(2./H1_ini)
      !coefct(k) = 1 + (kdtdx(k)*(1./dzeta(k+1)+1./dzeta(k)) - dt*(wp1/dzeta(k+1)-wm1/dzeta(k)))*(2./H1_ini)
      !coefm1(k) = 0. 
      !coef_f(k) = ice_T1(k) + (kdtdx(k)+dt*wm1)*(dzeta(k+1)+dzeta(k))/dzeta(k)*geoh_flux/ice_kappa 
    
      ! solution core
      call Tridiagonal_matrix_algorithm(coefm1,coefct,coefp1,ice_T0,coef_f)
      ice_T1 = ice_T0
  end do
  dif = ice_T0-T_ini
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine Tridiagonal_matrix_algorithm(coff_0,coff_1,coff_2,x,f)
!+++++++++++++++++++++++++++++++++++++++
! According to "The Theory of Splines and Their Applications" (1976) P15 by Ahlberg et al.
! solving tridiagonal matrix equation set with periodic condition
! output: matrix solution (x)
  implicit none
  integer, parameter    :: kdim     = 4                   ! vertical layer
  integer               :: k, n
  real*8,dimension(kdim)  :: coff_0, coff_1, coff_2, x, f  ! input coefficient
  real*8,dimension(kdim)  :: p, q, u, s                    ! variables for solving matrix
  real*8,dimension(kdim)  :: t, v                          ! variables for substitute
  n = kdim
  p(1)=coff_1(1);q(1)=-coff_2(1)/p(1);u(1)=f(1)/p(1);s(1)=-coff_0(1)/p(1);   ! initialization
  ! x(k) = q(k)*x(k+1)+s(k)*x(n)+u(k)
  do k= 2,n-1
      p(k) = coff_0(k)*q(k-1)+coff_1(k)    ! coefficient for diagonal element
      q(k) = -coff_2(k)/p(k)               ! coefficient for diagonal element next step
      u(k) = (f(k)-coff_0(k)*u(k-1))/p(k)  ! coefficient for offset
      s(k) = -coff_0(k)*s(k-1)/p(k)        ! coefficient for xn 
  end do
  ! x(k) = t(k)*x(n)+v(k)
  t(n)=1;v(n)=0 ! initialization
  do k= n-1,1,-1
      t(k) = q(k)*t(k+1)+s(k)              ! coefficient for xn
      v(k) = q(k)*v(k+1)+u(k)              ! coefficient for offset
  end do
  ! xn
  x(n) = (f(n)-v(n-1)*coff_0(n)-v(1)*coff_2(n))/(coff_1(n)+coff_0(n)*t(n-1)+coff_2(n)*t(1))
  ! xk
  do k= 1,n-1
      x(k) = t(k)*x(n) + v(k)
  end do
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_advection(vx, vy, var1, term_adv, dyy, dxlat, lat)
!+++++++++++++++++++++++++++++++++++++++
! FFSL advection scheme. 
! Lin, S. J., and R. B. Rood, 1996: Multidimensional flux-form semi-Lagrangian transport schemes.
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_2d.cpp
! output: advection term for ice thickness bedget (term_adv)
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
                         dlon, dlat, dt_ice
  USE mo_physics, ONLY: pi, ref_lat
  implicit none
  integer, external  :: cycle_ind
  integer, parameter :: num_spec = 128  ! number of spectrum points for fast Fourier transform
  integer, parameter :: num_ord  = 7    ! order of spectrum in pow of 2
  integer            :: i, j

  real*8, dimension(xdim,ydim)     :: var1, varx, vary, term_adv
  real*8, dimension(xdim,ydim)     :: termx, termy ! ice advection term x/y direction (work variable)
  real*8, dimension(xdim,ydim+1)   :: crantx, cranty, fu, fv, vx, vy
  real*8, dimension(ydim)          :: lat, dxlat
  real*8, dimension(num_spec)      :: fr_r, fr_i   ! Fourier coefficient
  real*8, dimension(num_spec)      :: signal_r, signal_i ! signal data
  real*8, dimension(num_spec)      :: filter       ! weight for wave filter
  integer, dimension(num_spec)   :: ispec = (/(i,i=1,num_spec)/)  ! index for spectrum
  real*8                           :: dyy

  fu = 0.; fv = 0.; crantx = 0.; cranty = 0.
  
  ! inner operator
  call ice_advection_upwind(vx, vy, var1, var1, termx, termy, dyy, dxlat, lat, var1)
  varx = var1 - termy/2.
  vary = var1 - termx/2.

  ! outer operator
  do j = 1,xdim
        do i = 1,ydim
            crantx(j,i) = vx(j,i)*dt_ice/dxlat(i)
            cranty(j,i) = vy(j,i)*dt_ice/dyy
            vary(j,i) = vary(j,i)*cos(lat(i)/180*pi)
        end do
        i = ydim+1
        cranty(j,i) = vy(j,i)*dt_ice/dyy
  end do
  
  ! advection term at x/y direction
  termx = 0.; termy = 0.
  call flux_operator(varx, vary, crantx, cranty, fu, fv)
  do j = 1,xdim
      do i = 1,ydim
          termx(j,i) = - (fu(cycle_ind(j+1,xdim),i)-fu(j,i)) 
          termy(j,i) = - (fv(j,i+1)-fv(j,i))/cos(lat(i)/180*pi)
      end do   
  end do

  ! wave filter for south pole
  do i = 1,2
      fr_r = 0.; fr_i = 0.
      signal_r = 0.; signal_i = 0.
      signal_r(1:xdim) = termx(:,i)
      ! Fourier transform
      call FFT(signal_r, signal_i, num_spec, num_ord, fr_r, fr_i, 0)
      do j = 1,num_spec
          fr_r(j) = min(1.,(cos(lat(i)/180*pi)/cos(ref_lat/180*pi)/sin((j-1)*pi/num_spec/2))**2)*fr_r(j)
          fr_i(j) = min(1.,(cos(lat(i)/180*pi)/cos(ref_lat/180*pi)/sin((j-1)*pi/num_spec/2))**2)*fr_i(j)
      end do
      ! Fourier inverse transform
      call FFT(fr_r, fr_i, num_spec, num_ord, signal_r, signal_i, 1)
      termx(:,i) = signal_r(1:xdim)
  end do
  
  term_adv = termx + termy

end subroutine

subroutine ice_advection_upwind(vx, vy, varx, vary, termx, termy, dyy, dxlat, lat, ice_H1)
! Finite difference upwind scheme 
! output: advection budget in x, y direction (termx, termy)
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
                         dlon, dlat, dt_ice
  USE mo_physics, ONLY: pi
  implicit none
  integer, external  :: cycle_ind
  integer            :: i, j, indp1, indm1
  
  real*8               :: vxm, vym
  real*8               :: crantx, cranty 
  real*8               :: crantxl, crantxr, crantyu, crantyd 
  real*8               :: varyp1, varym1
  real*8               :: termxl, termxr, termyu, termyd
  real*8, dimension(xdim,ydim)     :: varx, vary
  real*8, dimension(xdim,ydim)     :: termx, termy
  real*8, dimension(xdim,ydim)     :: ice_H1 
  real*8, dimension(xdim,ydim+1)   :: vx, vy
  real*8, dimension(ydim)          :: lat, dxlat
  real*8                           :: dyy

  termx  = 0.; termy  = 0.
  do j = 1,xdim
        do i = 1,ydim
            crantx = 0.; cranty = 0.
            if(ice_H1(j,i) <= 0.) cycle ! skip no ice grid
            indp1 = cycle_ind(j+1,xdim); indm1 = cycle_ind(j-1,xdim)
            ! no temperature change for divergence grid
            if(vx(j,i)<0. .and. vx(indp1,i  )>0.) termx(j,i) = 0.
            if(vy(j,i)<0. .and. vy(j    ,i+1)>0.) termy(j,i) = 0.
            ! temperature change for convergence grid
            if(vx(j,i)>0. .and. vx(indp1,i  )<0.) then
                crantxl = vx(j    ,i)*dt_ice/dxlat(i)
                termxl = crantxl*(varx(j    ,i)-varx(indm1,i))
                crantxr = vx(indp1,i)*dt_ice/dxlat(i)
                termxr = crantxr*(varx(indp1,i)-varx(j    ,i))
                termx(j,i)  = (termxl*abs(crantxl)+termxr*abs(crantxr))/(abs(crantxl)+abs(crantxr))
            end if

            if(vy(j,i)>0. .and. vy(j,i+1)<0) then
                crantyd = vy(j, i  )*dt_ice/dyy
                call meridion_shift(vary, j, i-1, varym1)
                termyd  = crantyd*(vary(j,i) - varym1   ) 
                crantyu = vy(j, i+1)*dt_ice/dyy
                call meridion_shift(vary, j, i+1, varyp1)
                termyu  = crantyu*(varyp1    - vary(j,i))
                termy(j,i)   = (termyd*abs(crantyd)+termyu*abs(crantyu))/(abs(crantyd)+abs(crantyu))
            end if
            
            ! advection for consistent flow grid
            if(vx(j,i)*vx(indp1,i) > 0.) then
                if(vx(j,i) > 0.) then
                    crantx = vx(j    ,i)*dt_ice/dxlat(i)
                    termx(j,i)  = crantx*(varx(j    ,i)-varx(indm1,i))
                else
                    crantx = vx(indp1,i)*dt_ice/dxlat(i)
                    termx(j,i)  = crantx*(varx(indp1,i)-varx(j    ,i))
                end if
            end if
           
            if(vy(j,i)*vy(j,i+1) > 0.) then
                if(vy(j,i) > 0.) then
                    cranty = vy(j, i  )*dt_ice/dyy
                    call meridion_shift(vary, j, i-1, varym1)
                    termy(j,i)  = cranty*(vary(j,i) - varym1   ) 
                else
                    cranty = vy(j, i+1)*dt_ice/dyy
                    call meridion_shift(vary, j, i+1, varyp1)
                    termy(j,i)  = cranty*(varyp1    - vary(j,i))
                end if
            end if

        end do
  end do
     
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine flux_operator(ice_Hx, ice_Hy, crantx, cranty, fu, fv)
!+++++++++++++++++++++++++++++++++++++++
! flux operator for advection scheme
! output: boundary flux budget (fu, fv)
  USE mo_numerics, ONLY: xdim, ydim
  implicit none
  integer, external             :: cycle_ind
  integer                       :: i, j, k 
  real*8                          :: iceH_indp2, iceH_indp1, iceH_ind, iceH_indm1, iceH_indm2 
  real*8                          :: xl, xr, dx1, dx2, dx3, fl, df, f6 
  real*8                          :: yl, yr, dy1, dy2, dy3
  integer                       :: crant_intx, crant_inty, adv_ind 
  real*8                          :: crant_frax, crant_fray 
  real*8, dimension(xdim,ydim)    :: ice_Hx, ice_Hy
  real*8, dimension(xdim,ydim+1)  :: ice_vx, ice_vy, crantx, cranty, fu, fv
  
  do i = 1,ydim
      do j = 1,xdim
      if(crantx(j,i)==0. .and. cranty(j,i)==0.) cycle
          ! x direction
          crant_intx = int(crantx(j,i))
          crant_frax = crantx(j,i) - float(crant_intx)
          
          ! integer flux
          if(crant_intx .ge. 1) then
              do k = 1,crant_intx
                 fu(j,i) = fu(j,i) + ice_Hx(cycle_ind(j-k,xdim),i) 
              end do
          else if(crant_intx .le. -1) then
              do k = 1,-crant_intx
                 fu(j,i) = fu(j,i) - ice_Hx(cycle_ind(j-1+k,xdim),i)
              end do
          endif
          ! fraction flux
          if(crantx(j,i) > 0) then
              adv_ind = cycle_ind(j-1-crant_intx, xdim)
          else
              adv_ind = cycle_ind(j-crant_intx, xdim)
          end if
          if(crant_frax > 0) then
              xl = 1-crant_frax; xr = 1
          else
              xl = 0; xr = -crant_frax
          end if
          iceH_indm2 = ice_Hx(cycle_ind(adv_ind-2,xdim), i); iceH_indm1 = ice_Hx(cycle_ind(adv_ind-1,xdim), i); 
          iceH_ind   = ice_Hx(adv_ind, i);
          iceH_indp2 = ice_Hx(cycle_ind(adv_ind+2,xdim), i); iceH_indp1 = ice_Hx(cycle_ind(adv_ind+1,xdim), i);
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dx1 = xr-xl; dx2 = xr*xr-xl*xl; dx3 = xr*xr*xr-xl*xl*xl;
          fu(j,i) = fu(j,i) + sign( fl*dx1+0.5*df*dx2+f6*(0.5*dx2-dx3/3.0), crant_frax);
          ! y direction
          crant_inty = int(cranty(j,i))
          crant_fray = cranty(j,i) - float(crant_inty)

          ! integer flux
          if(crant_inty .ge. 1) then
              do k = 1,crant_inty
                 call meridion_shift(ice_Hy, j, i-k, iceH_ind)
                 fv(j,i) = fv(j,i) + iceH_ind
              end do
          else if(crant_inty .le. -1) then
              do k = 1,-crant_inty
                 call meridion_shift(ice_Hy, j, i-1+k, iceH_ind)
                 fv(j,i) = fv(j,i) - iceH_ind
              end do
          endif
          ! fraction flux
          if(cranty(j,i) > 0) then
              adv_ind = i-1-crant_inty
          else
              adv_ind = i-crant_inty
          end if
          if(crant_fray > 0) then
              yl = 1-crant_fray; yr = 1
          else
              yl = 0; yr = -crant_fray
          end if
         
          if ((adv_ind>ydim-2).or.(adv_ind<2)) then
              call meridion_shift(ice_Hy, j, adv_ind-2, iceH_indm2)
              call meridion_shift(ice_Hy, j, adv_ind-1, iceH_indm1)
              call meridion_shift(ice_Hy, j, adv_ind  , iceH_ind  )
              call meridion_shift(ice_Hy, j, adv_ind+1, iceH_indp1)
              call meridion_shift(ice_Hy, j, adv_ind+2, iceH_indp2)
          else
              iceH_indm2 = ice_Hy(j,adv_ind-2); iceH_indm1 = ice_Hy(j,adv_ind-1); iceH_ind = ice_Hy(j,adv_ind);
              iceH_indp2 = ice_Hy(j,adv_ind+2); iceH_indp1 = ice_Hy(j,adv_ind+1);
          end if
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dy1 = yr-yl; dy2 = yr*yr-yl*yl; dy3 = yr*yr*yr-yl*yl*yl;
          fv(j,i) = fv(j,i) + sign(fl*dy1+0.5*df*dy2+f6*(0.5*dy2-dy3/3.0),crant_fray);         
      end do
  end do
  i = ydim + 1
  do j = 1,xdim
          ! y direction
          crant_inty = int(cranty(j,i))
          crant_fray = cranty(j,i) - float(crant_inty)

          ! integer flux
          if(crant_inty .ge. 1) then
              do k = 1,crant_inty
                 call meridion_shift(ice_Hy, j, i-k, iceH_ind)
                 fv(j,i) = fv(j,i) + iceH_ind
              end do
          else if(crant_inty .le. -1) then
              do k = 1,-crant_inty
                 call meridion_shift(ice_Hy, j, i-1+k, iceH_ind)
                 fv(j,i) = fv(j,i) - iceH_ind
              end do
          endif
          ! fraction flux
          if(cranty(j,i) > 0) then
              adv_ind = i-1-crant_inty
          else
              adv_ind = i-crant_inty
          end if
          if(crant_fray > 0) then
              yl = 1-crant_fray; yr = 1
          else
              yl = 0; yr = -crant_fray
          end if
         
          if ((adv_ind>ydim-2).or.(adv_ind<2)) then
              call meridion_shift(ice_Hy, j, adv_ind-2, iceH_indm2)
              call meridion_shift(ice_Hy, j, adv_ind-1, iceH_indm1)
              call meridion_shift(ice_Hy, j, adv_ind  , iceH_ind  )
              call meridion_shift(ice_Hy, j, adv_ind+1, iceH_indp1)
              call meridion_shift(ice_Hy, j, adv_ind+2, iceH_indp2)
          else
              iceH_indm2 = ice_Hy(j,adv_ind-2); iceH_indm1 = ice_Hy(j,adv_ind-1); iceH_ind = ice_Hy(j,adv_ind);
              iceH_indp2 = ice_Hy(j,adv_ind+2); iceH_indp1 = ice_Hy(j,adv_ind+1);
          end if
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dy1 = yr-yl; dy2 = yr*yr-yl*yl; dy3 = yr*yr*yr-yl*yl*yl;
          fv(j,i) = fv(j,i) + sign(fl*dy1+0.5*df*dy2+f6*(0.5*dy2-dy3/3.0), crant_fray);         
  end do
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)
!+++++++++++++++++++++++++++++++++++++++
! boundary flux for FFSL scheme
! coefficient for PPM (fl, df, f6)
  implicit none
  real*8 :: fm2, fm1, f, fp1, fp2 
  real*8 :: dfl, df, dfr, fl, fr, f6
  dfl     = mismatch(fm2, fm1, f)
  df      = mismatch(fm1, f  , fp1)
  dfr     = mismatch(f  , fp1, fp2)
  fl      = 0.5*(fm1+f)+(dfl-df)/6.
  fr      = 0.5*(fp1+f)+(df-dfr)/6.

  fl = f-sign(min(abs(df), abs(fl-f)), df);
  fr = f+sign(min(abs(df), abs(fr-f)), df);
  f6 = 6*f-3*(fl+fr);
  df = fr-fl;

  contains
      real*8 function mismatch(fm1, f, fp1)  
      ! the slope calculation in (Colella & Woodward, 1984)
      real*8 fm1, f, fp1, df, dfMin, dfMax
      if ((fp1-f)*(f-fm1) .gt. 0) then
          df       = (fp1 - fm1)*0.5
          dfMin    = f - min(fm1, f , fp1)
          dfMax    = max(fm1, f, fp1) - f
          mismatch = sign(min(abs(df), 2*dfMin, 2*dfMax), df) ! monontonic
          !mismatch = sign(min(abs(df), 2*f), df) ! positive defined
      else
          mismatch = 0
      end if
      end function mismatch
end

integer function cycle_ind(x,xdim)
! zonal grid shift, cyclic boundary condition
! output: zonal index 
  implicit none
  integer :: x, xdim
  if(x < 1) then
    cycle_ind = xdim + x ! westmost
   else if(x > xdim) then
    cycle_ind = x - xdim ! eastmost
   else
    cycle_ind = x ! otherwise
   endif
end function

subroutine meridion_shift(ice_H1, ind_zonal, ind_merid, ice_nounb)
! meridian grid shift (Allen et al., 1991)
! output: meridional index
 USE mo_numerics, ONLY: xdim, ydim
 implicit none
 integer, external          :: cycle_ind
 integer                    :: ind_zonal, ind_merid
 real*8                       :: ice_nounb
 real*8, dimension(xdim,ydim) :: ice_H1
 if   (ind_merid.gt.ydim) then
      ice_nounb = ice_H1(cycle_ind(ind_zonal+xdim/2,xdim),2*ydim+1-ind_merid) ! North Pole, move to opposite longitude 
 else if(ind_merid.lt.1 ) then
      ice_nounb = ice_H1(cycle_ind(ind_zonal+xdim/2,xdim), 1-ind_merid) ! South Pole, move to opposite longitude 
 else 
      ice_nounb = ice_H1(ind_zonal,ind_merid) ! otherwise
 end if
end subroutine

subroutine FFT(PR,PI,N,K,FR,FI,L)
! Fast Fourier Transform 
! PR, PI is input real and image part of data/Fourier coefficient
! FR, FI is output real and image part of Fourier coefficient/data
! N is the sample size, it is equal to 2^K
! L is flag for tansform direction, 0 for Fourier transform, 1 for inverse
implicit none
real*8,dimension(N) :: PR, PI ,FR ,FI
integer           :: N, K, L
real*8              :: P, Q, S, VR, VI, PODDR, PODDI
integer           :: IT, M, IS, I, J, NV
do IT=0,N-1
    M=IT
    IS=0
    do I=0,K-1
            J=M/2
            IS=2*IS+(M-2*J)
            M=J
    end do
    FR(IT+1)=PR(IS+1)
    FI(IT+1)=PI(IS+1)
end do
PR(1)=1.0
PI(1)=0.0
PR(2)=COS(6.283185306/N)
PI(2)=-SIN(6.283185306/N)
IF (L.NE.0) PI(2)=-PI(2)
do I=3,N
    P=PR(I-1)*PR(2)
    Q=PI(I-1)*PI(2)
    S=(PR(I-1)+PI(I-1))*(PR(2)+PI(2))
    PR(I)=P-Q
    PI(I)=S-P-Q
end do
do IT=0,N-2,2
    VR=FR(IT+1)
    VI=FI(IT+1)
    FR(IT+1)=VR+FR(IT+2)
    FI(IT+1)=VI+FI(IT+2)
    FR(IT+2)=VR-FR(IT+2)
    FI(IT+2)=VI-FI(IT+2)
end do
M=N/2
NV=2
do I=K-2,0,-1
    M=M/2
    NV=2*NV
    do IT=0,(M-1)*NV,NV
        do J=0,(NV/2)-1
            P=PR(M*J+1)*FR(IT+J+1+NV/2)
            Q=PI(M*J+1)*FI(IT+J+1+NV/2)
            S=PR(M*J+1)+PI(M*J+1)
            S=S*(FR(IT+J+1+NV/2)+FI(IT+J+1+NV/2))
            PODDR=P-Q
            PODDI=S-P-Q
            FR(IT+J+1+NV/2)=FR(IT+J+1)-PODDR
            FI(IT+J+1+NV/2)=FI(IT+J+1)-PODDI
            FR(IT+J+1)=FR(IT+J+1)+PODDR
            FI(IT+J+1)=FI(IT+J+1)+PODDI
        end do
    end do
end do
if (L.NE.0) THEN
    do I=1,N
        FR(I)=FR(I)/N
        FI(I)=FI(I)/N
    end do
end if
RETURN
end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_topography(ice_zs, ice_H1)
!+++++++++++++++++++++++++++++++++++++++
! topograhy, glacier type for grid 
! implicit output: surface elevation (z_topo), glacier type (glacier_type) 
  USE mo_numerics,     ONLY: xdim, ydim
  USE mo_physics,      ONLY: b_rock, z_topo, mask, rho_ice, rho_ocean, glacier_type, log_exp, d_ice_mov
  implicit none
  integer, external          :: cycle_ind
  real*8, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (forward) [m]
  real*8, dimension(xdim,ydim) :: ice_zs                  ! ice surface height [m]
  real*8, dimension(xdim,ydim) :: bed0                    ! bed rock data (real*8) [m]
  integer                      :: i, j                    ! work variable
  ! surface elevation
  z_topo = -0.1
  where(ice_zs  > 0.) z_topo = ice_zs

  ! glacier type
  bed0 = b_rock
  do j = 1,xdim
      do i = 1,ydim
          ! ice sheet reaches the bed rock - grounded ice
          if(bed0(j,i) >= -ice_H1(j,i)*rho_ice/rho_ocean) then
              glacier_type(j,i) = 1
          else
          ! ice sheet cannot reach the bed rock
              if(ice_H1(j,i) > d_ice_mov) then
          ! near to grounded ice - floating ice
                  glacier_type(j,i) = 0
              else
          ! open ocean
                  glacier_type(j,i) = -1
              end if
          end if
      end do
  end do
  if(log_exp == 299) glacier_type = 1

end subroutine ice_topography

