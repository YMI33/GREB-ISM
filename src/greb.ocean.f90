!+++++++++++++++++++++++++++++++++++++++
subroutine seaice(Tsurf, ice_H1, SW, LW_surf, LWair_down, Q_lat,  &
&                     Q_sens, dT_ocean, Q_sice, dice_h)
!+++++++++++++++++++++++++++++++++++++++
! prognostic sea ice model

  USE mo_numerics,    ONLY: xdim, ydim,nstep_yr, dt
  USE mo_physics,     ONLY: ityr, mask, cap_surf, cap_land, cap_ocean,             &
&                           log_exp, To_ice1, To_ice2, mldclim,                    &
&                           TF_correct, ci_latent, rho_ice 

! declare temporary fields
  integer it
  real, dimension(xdim,ydim)   :: Tsurf, SW, LW_surf, LWair_down,    &
&  								  Q_lat, Q_sens, dT_ocean, Q_sice, ice_Tse
  real*8,dimension(xdim,ydim)  :: ice_H1, netheat, dtx, dice_h, &
&								 hmax_melt, hmax_ice, cool_tsurf, ice_x,    &
&                                wice, heat_tsurf

  1389 format ('sea ice ', F6.2, F8.4, 4F10.2, F9.3, F8.3) !TB
  1390 format ('sea ice ', I4, F6.2, F6.4, 5F6.2, F9.6) !TB

  !----------------------------------------
  ! heat capacity over oceans (change)
  where(mask < 0. .and. ice_h1 >  0.2 ) cap_surf = cap_land         ! sea ice
  where(mask < 0. .and. ice_h1 <= 0.2 .and. ice_h1 > 0 )  &
  & cap_surf = cap_land*ice_H1/0.2 + cap_ocean*mldclim(:,:,ityr)*(1-ice_H1/0.2) ! sea ice
  where(mask < 0. .and. ice_h1 == 0.) cap_surf = cap_ocean*mldclim(:,:,ityr) ! open ocean

  netheat = SW +LW_surf -LWair_down +Q_lat +Q_sens +TF_correct(:,:,ityr)/2.
  dtx     = dT_ocean +dt*netheat/cap_surf      ! potential change in Tsurf

  dice_h = 0.0; Q_sice = 0.0; heat_tsurf = 0.0; hmax_melt = 0.0;
  ice_Tse  = Tsurf + dtx
  where( mask < 0)
      ! the extra heat flux after the ice reaches melting point
      heat_tsurf = (ice_Tse - To_ice2) * cap_surf / dt
      ! the heat flux needed to melt whole ice column
      where(ice_H1 > 0.) hmax_melt  = ci_latent * rho_ice * ice_H1 / dt
  end where
  ! equation (32) in XIE2021
  ! sea ice totally melts away
  where( (mask < 0) .and. (heat_tsurf >= 0 ) .and. (heat_tsurf >  hmax_melt) .and. (ice_H1 >0.))
        Q_sice      = - hmax_melt 
  end where
  ! sea ice partially melts
  where( (mask < 0) .and. (heat_tsurf >= 0 ) .and. (heat_tsurf <= hmax_melt) .and. (ice_H1 >0.))
        Q_sice      = - heat_tsurf 
  end where
  ! sea ice forms 
  where( (mask < 0) .and. (heat_tsurf < 0 ) .and. (netheat < 0) .and. (ice_H1 < 0.5 ) )
        Q_sice      = - netheat
  end where
  ! equation (31) in XIE2021
  dice_h      = dt*Q_sice / (rho_ice*ci_latent)

  !----------------------------------------
  ! sea ice tendencies due to transport
  cicediff = 1/4.*1./(30.*24.*3600.) ! time scale of ice diffusion [1/sec]

  do i=1,xdim
     ip = i+1; if (ip > xdim) ip = 1
     in = i-1; if (in <    1) in = xdim
     do j=2,ydim-1
        wz =  1./cos(2.*3.1416/360.*(3.75*j-90.-3.75/2.)) ! zonal weight
        jp = j+1; jn = j-1;
        if (mask(i,j) < 0. .and. ice_H1(i,j) > 0.) then ! sea ice
            dice_dif = 0.0
            if (mask(ip,j) < 0.) dice_dif = wz*cicediff*(ice_H1(ip,j)-ice_H1(i,j))
            if (mask(in,j) < 0.) dice_dif = dice_dif +wz*cicediff*(ice_H1(in,j)-ice_H1(i,j))
            if (mask(i,jp) < 0. .and. j < ydim)  dice_dif = dice_dif +cicediff*(ice_H1(i,jp)-ice_H1(i,j))
            if (mask(i,jn) < 0.) dice_dif = dice_dif +cicediff*(ice_H1(i,jn)-ice_H1(i,j))
             dice_dif = dt*dice_dif
            dice_h(i,j) = dice_h(i,j) +dice_dif
        end if
     end do     
  end do

end subroutine seaice

!+++++++++++++++++++++++++++++++++++++++
subroutine deep_ocean(Ts, To, dT_ocean, dTo)
!+++++++++++++++++++++++++++++++++++++++
!              deep ocean model

  USE mo_numerics,    ONLY: xdim, ydim, nstep_yr, dt
  USE mo_physics,     ONLY: ityr, mask, mldclim, log_exp, To_ice2,       &
&                           cap_ocean, co_turb, z_ocean 

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts, To, dT_ocean, dTo, dmld, Tx
  dT_ocean = 0.0;  dTo     = 0.0

  10 format ('deep ocean: ',10E10.2) !TB

  if (ityr >  1) dmld = mldclim(:,:,ityr)-mldclim(:,:,ityr-1)
  if (ityr == 1) dmld = mldclim(:,:,ityr)-mldclim(:,:,nstep_yr)

! entrainment & detrainment
  where ( mask < 0 .and. Ts >= To_ice2 .and. dmld < 0)       &
&       dTo      = -dmld/(z_ocean-mldclim(:,:,ityr))*(Ts-To)
  where ( mask < 0 .and. Ts >= To_ice2 .and. dmld > 0)       &
&       dT_ocean =  dmld/mldclim(:,:,ityr)*(To-Ts)

  c_effmix = 0.5
  dTo      = c_effmix*dTo
  dT_ocean = c_effmix*dT_ocean

! turbulent mixing
  Tx = max(To_ice2,Ts)
  where ( mask < 0 ) dTo      = dTo      + dt*co_turb*(Tx-To)/(cap_ocean*(z_ocean-mldclim(:,:,ityr)))
  where ( mask < 0 ) dT_ocean = dT_ocean + dt*co_turb*(To-Tx)/(cap_ocean*mldclim(:,:,ityr))

end subroutine deep_ocean

!+++++++++++++++++++++++++++++++++++++++
subroutine sealevel(ice_h, Ts, To)
!+++++++++++++++++++++++++++++++++++++++
! compute changes in topography and land-sea mask due to ice sheet heights/mass

  USE mo_numerics,    ONLY: xdim, ydim, nstep_yr
  USE mo_physics,     ONLY: ityr, mask, z_topo, cap_surf, cap_land, cap_ocean, mldclim, &
&                           glacier_type, TF_correct, b_rock, b_rock0, rho_ice, rho_ocean, &
&                           swetclim, iceH_clim, ssh, log_exp, log_ice_slv

! declare temporary fields
  implicit none
  real*8, external              :: gsum 
  real,   dimension(xdim,ydim)  :: Ts
  real,   dimension(xdim,ydim)  :: To
  real*8, dimension(xdim,ydim)  :: ice_h
  real*8, dimension(xdim,ydim)  :: ice_hl ! ice sheet thickness at land [m]
  real*8, dimension(xdim,ydim)  :: ice_h0 ! ice sheet thickness reference [m]
  real*8, dimension(xdim,ydim)  :: ogrid  ! ocean grid 
  integer i, ix, iy
   
10 format ('SSH: ',15F10.2) !TB
   ix =85; iy=9

  ! equation (36) in XIE2021
  ogrid = 0.; ice_hl = 0.
  where(b_rock  <  -ice_h *rho_ice/rho_ocean) ogrid  = 1 
  where(b_rock  >= -ice_h *rho_ice/rho_ocean) ice_hl = ice_h 
  ice_h0 = iceH_clim(:,:,ityr)
  where(b_rock0 <  -ice_h0*rho_ice/rho_ocean) ice_h0 = 0.
  ssh    = - (gsum(ice_hl-ice_h0)*rho_ice)/(rho_ocean*gsum(ogrid))
  if(log_ice_slv == 0) ssh = 0.
  b_rock = b_rock0 - ssh 
  
  ! change land sea mask
  do i=1, nstep_yr
     where(mask < 0 .and. b_rock > 0.0 ) swetclim(:,:,i) = 0.3                    ! change ocean -> land
     where(mask > 0 .and. b_rock < -1.0 .and. ice_h < 0.2 ) swetclim(:,:,i) = 1.0 ! change land -> ocean
     where(mask > 0 .and. b_rock < -1.0 .and. ice_h < 0.2 ) To = Ts               ! change land -> ocean
  end do
  where(mask < 0 .and. b_rock > 0.0 )                    mask =  1 ! change ocean -> land
  where(mask > 0 .and. b_rock < -1.0 .and. ice_h < 0.2 ) mask = -1 ! change land -> ocean
  
  where(mask > 0.) cap_surf = cap_land
  where(mask < 0.) cap_surf = cap_ocean*mldclim(:,:,ityr)
  
end subroutine sealevel

!+++++++++++++++++++++++++++++++++++++++
function gsum(data)
!+++++++++++++++++++++++++++++++++++++++
 
use mo_numerics,		ONLY: xdim, ydim, dlat, dlon
use mo_physics,         ONLY: pi
 
! declare variables
real*8, dimension(xdim,ydim)  :: data, w
real*8, dimension(ydim)	      :: lat
real*8                        :: gsum
 
do i=1,ydim
lat(i) = -90-(dlat*0.5)+i*dlat
end do
do i=1,xdim
w(i,:) = cos(2.*3.14159*lat/360.)
end do
 
gsum = sum(data*w)
 
end function
