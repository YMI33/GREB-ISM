!+++++++++++++++++++++++++++++++++++++++
subroutine SWradiation(Tsurf, sw, ice_cover, a_surf)
!+++++++++++++++++++++++++++++++++++++++
!    SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, sw_solar,da_ice, a_no_ice, a_cloud, mask    &
&                         , Tl_ice1, Tl_ice2, To_ice1, To_ice2 & !,glacier       
&                         , cldclim, log_exp, S0_var  &
&						  , log_ice_cpl, z_topo

! declare temporary fields
  real,   dimension(xdim,ydim)  :: Tsurf, sw, albedo, a_surf, a_atmos
  real*8, dimension(xdim,ydim)  :: ice_cover

! atmos albedo
  a_atmos=cldclim(:,:,ityr)*a_cloud

! surface albedo
   a_surf = a_no_ice          											 ! no ice

  if(log_ice_cpl == 0) then

     ! Land:  ice -> albedo linear function of T_surf
   	 where(mask > 0. .and. Tsurf <= Tl_ice1) a_surf = a_no_ice+da_ice   ! ice
   	 where(mask > 0. .and. Tsurf > Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&         a_surf = a_no_ice +da_ice*(1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))

    ! Ocean: ice -> albedo/heat capacity linear function of T_surf
     where(mask < 0. .and. Tsurf <= To_ice1) a_surf = a_no_ice+da_ice      ! ice
     where(mask < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       a_surf = a_no_ice+da_ice*(1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

     ! glacier -> no albedo changes
     ! where(glacier > 0.5) a_surf = a_no_ice+da_ice
  end if 

  if(log_ice_cpl == 1) then
     ! with ice sheet model
  
     ! Land: snow cover 
     where(mask > 0. .and. ice_cover >= 0.02) a_surf = a_no_ice+da_ice   ! ice
     where(mask > 0. .and. ice_cover > 0. .and. ice_cover < 0.02 ) &
&       a_surf =  a_no_ice +da_ice*ice_cover/0.02

     ! Ocean: ice -> albedo/heat capacity linear function of T_surf
     where(mask < 0. .and. ice_cover >= 0.5) a_surf = a_no_ice+da_ice      ! ice
     where(mask < 0. .and. ice_cover  > 0.0 .and. ice_cover < 0.5 ) &
&          a_surf = a_no_ice + da_ice*ice_cover/0.5
     
  end if 
 
! SW flux
  albedo=a_surf+a_atmos-a_surf*a_atmos
  forall (i=1:xdim)
     sw(i,:)=0.01*S0_var*SW_solar(:,ityr)*(1-albedo(i,:))
  end forall

end subroutine SWradiation

!+++++++++++++++++++++++++++++++++++++++
subroutine senseheat(Ta1,Ts1,Q_sens)
!+++++++++++++++++++++++++++++++++++++++
! compute the sensible heat exchange between Tsurf and Tatmos

  use mo_numerics, only: xdim, ydim
  use mo_physics,  only: ct_sens, z_topo, c_lapse

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts1, Ta1, Q_sens, Ta_scl

  Ta_scl = Ta1 + c_lapse*z_topo

  Q_sens = ct_sens*(Ta_scl-Ts1)

end subroutine senseheat

!+++++++++++++++++++++++++++++++++++++++
subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
!+++++++++++++++++++++++++++++++++++++++
! new approach with LW atmos

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: sig, eps, qclim, cldclim, z_topo, jday, ityr,         &
&                           r_qviwv, z_air, z_vapor, dTrad, p_emi,                &
&                           co2_part

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tair, q, LWsurf, LWair, e_co2, e_cloud,   &
&                                LWair_up, LWair_down, e_vapor, em


  e_co2   = exp(-z_topo/z_air)*co2_part*CO2   ! CO2
  e_vapor = exp(-z_topo/z_air)*r_qviwv*q      ! water vapor
  e_cloud = cldclim(:,:,ityr)                 ! clouds

! total
  em      = p_emi(4)*log( p_emi(1)*e_co2 +p_emi(2)*e_vapor +p_emi(3) ) +p_emi(7)   &
&          +p_emi(5)*log( p_emi(1)*e_co2   +p_emi(3) )                             &
&          +p_emi(6)*log( p_emi(2)*e_vapor +p_emi(3) )
  em      = (p_emi(8)-e_cloud)/p_emi(9)*(em-p_emi(10))+p_emi(10)

  LWsurf      = -sig*Tsurf**4
  LWair_down  = -em*sig*(Tair+dTrad(:,:,ityr))**4
  LWair_up    = LWair_down

end subroutine LWradiation

!+++++++++++++++++++++++++++++++++++++++
subroutine hydro(Tsurf, q, Qlat, Qlat_air, dq_eva, dq_rain, ice_H1)
!+++++++++++++++++++++++++++++++++++++++
!    hydrological model for latent heat and water vapor

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: rho_air, uclim, vclim, z_topo, swetclim, ityr,   &
&                           ce, cq_latent, cq_rain, z_air, r_qviwv, &
&                           omegaclim, omegastdclim, wsclim,  wz_vapor,      &
&                           c_q, c_rq, c_omega, c_omegastd, log_exp,  &    ! Rainfall parameters
&						    mask, precip_correct

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tskin, q, Qlat, Qlat_air, qs, dq_eva, &
&                                dq_rain, abswind, rq, q_zonal
  real*8,dimension(xdim,ydim) :: ice_H1 

  Qlat=0.; Qlat_air=0.; dq_eva=0.; dq_rain=0.

  abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2)
  where(mask > 0. ) abswind = sqrt(abswind**2 +2.0**2) ! land
  where(mask < 0. ) abswind = sqrt(abswind**2 +3.0**2) ! ocean

! saturated humiditiy (max. air water vapor)
  qs = 3.75e-3*exp(17.08085*(Tsurf-273.15)/(Tsurf-273.15+234.175));
  qs = qs*exp(-z_topo/z_air) ! scale qs by topography

! relative humidity
  rq = q/qs

! latent heat flux surface
  if ( log_eva == -1 ) then
    abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2)
    where(mask > 0. ) abswind = sqrt(abswind**2 + 2.0**2) !< land turbulent wind
    where(mask < 0. ) abswind = sqrt(abswind**2 + 3.0**2) !< ocean turbulent wind
    Qlat   = (q-qs)*abswind*cq_latent*rho_air*ce*swetclim(:,:,ityr) ! latend heat flux
  else if ( log_eva == 1 ) then
    where(mask > 0. ) abswind = sqrt(abswind**2 + 144.**2) ! land turbulent wind
    where(mask < 0. ) abswind = sqrt(abswind**2 + 7.1**2) ! ocean turbulent wind
    where(mask > 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.04*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(mask < 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.73*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  else if ( log_eva == 2 ) then
    abswind = wsclim(:,:,ityr) ! use the wind speed climatology
    where(mask > 0. ) abswind = sqrt(abswind**2 + 9.0**2) ! land turbulent wind
    where(mask < 0. ) abswind = sqrt(abswind**2 + 4.0**2) ! ocean turbulent wind
    where(mask > 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.56*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(mask < 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.79*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  else if ( log_eva == 0 ) then
    where(mask > 0. ) Tskin = Tsurf + 5. ! skin temperature land
    where(mask < 0. ) Tskin = Tsurf + 1. ! skin temperature ocean
    qs = 3.75e-3*exp(17.08085*(Tskin-273.15)/(Tskin-273.15+234.175)) ! re-calculate saturation pressure
    qs = qs*exp(-z_topo/z_air) ! scale qs by topography
    where(mask > 0. ) abswind = sqrt(wsclim(:,:,ityr)**2 + 11.5**2) ! land turbulent wind
    where(mask < 0. ) abswind = sqrt(wsclim(:,:,ityr)**2 + 5.4**2) ! ocean turbulent wind
    where(mask > 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.25*ce*swetclim(:,:,ityr) ! latend heat flux land
    where(mask < 0. ) Qlat    = (q-qs)*abswind*cq_latent*rho_air*0.58*ce*swetclim(:,:,ityr) ! latend heat flux ocean
  end if
! change in water vapor
  dq_eva  = -Qlat/cq_latent/r_qviwv  ! evaporation

! precipitation -> Eq. 11 in Stassen et al 2019
! Parameters in unused terms are set to zero
  dq_rain = (c_q + c_rq*rq + c_omega*omegaclim(:,:,ityr) + c_omegastd*omegastdclim(:,:,ityr))*cq_rain*q    ! unit: kg/s
  where(dq_rain >= -0.0015 / (wz_vapor * r_qviwv * 86400.)) dq_rain = -0.0015 / (wz_vapor * r_qviwv * 86400.) !Avoid negative rainfall (dq_rain is negative means positive rainfall!)

  ! precip flux corrections (land only)
  do i=1,ydim
      q_zonal(:,i) = sum(q(:,i))/xdim
  end do
  do i=1,7; q_zonal(:,i) = sum(q_zonal(1,i:8))/(8-i+1); end do ! smooth antarctica to s-ocean
  do i=ydim-6,ydim; q_zonal(:,i) = sum(q_zonal(1,ydim-7:i))/(8-(ydim-i)); end do ! smooth n-pole

  dq_rain = dq_rain + q_zonal*precip_correct(:,:,ityr)

! latent heat flux atmos
  Qlat_air = -dq_rain*cq_latent*r_qviwv
  
end subroutine hydro

!+++++++++++++++++++++++++++++++++++++++
subroutine circulation(X_in, dX_crcl, h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
! circulation with shorter time step

  USE mo_numerics,  ONLY: xdim, ydim, dt, dt_crcl
  USE mo_physics,   ONLY: z_vapor, z_air, &
&                         log_vdif, log_vadv, log_hdif, log_hadv, log_conv
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: X_in, wz
  real,                       intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_crcl

  real, dimension(xdim,ydim) :: X, dx_diffuse, dx_advec, dx_conv
  integer time, tt

  dX_crcl    = 0.0

  dx_diffuse = 0.0
  dx_advec   = 0.0
  dx_conv    = 0.0
  time=max(1,nint(float(dt)/dt_crcl))

  X = X_in;
  do tt=1, time   ! time loop circulation
! dmc & decon2xco2 switch
     if (log_vdif == 1 .and. h_scl .eq. z_vapor) call diffusion(X, dx_diffuse, h_scl, wz)
     if (log_vadv == 1 .and. h_scl .eq. z_vapor) call advection(X, dx_advec, h_scl, wz)
     if (log_conv == 0 .and. h_scl .eq. z_vapor) call convergence(X, dx_conv)
     if (log_hdif == 1 .and. h_scl .eq. z_air)   call diffusion(X, dx_diffuse, h_scl, wz)
     if (log_hadv == 1 .and. h_scl .eq. z_air)   call advection(X, dx_advec, h_scl, wz)
     X = X + dx_diffuse + dx_advec + dx_conv
  end do           ! time loop
  dX_crcl = X - X_in

end subroutine circulation

!+++++++++++++++++++++++++++++++++++++++
subroutine diffusion(T1, dX_diffuse,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    diffusion

  USE mo_numerics,   ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,    ONLY: pi, z_topo, kappa, z_vapor
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_diffuse

  integer :: i
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/)
  real, dimension(ydim)      :: lat, dxlat, ccx
  real, dimension(xdim)      :: T1h, dTxh
  real, dimension(xdim,ydim) :: dTx, dTy

  real    :: deg, dd, dx, dy, dyy, ccy, ccx2
  integer :: j, k, km1, kp1, jm1, jm2, jm3, jp1, jp2, jp3
  integer :: time2, dtdff2, tt2

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=kappa*dt_crcl/dyy**2
  ccx=kappa*dt_crcl/dxlat**2

     ! latitudinal
     do k=1, ydim
        km1=k-1;  kp1=k+1
        if ( k>=2 .and. k<=ydim-1)   dTy(:,k)=ccy*(                                      &
&                         wz(:,km1)*(T1(:,km1)-T1(:,k)) +wz(:,kp1)*(T1(:,kp1)-T1(:,k)) )
        if ( k==1 )                  dTy(:,k)=ccy*wz(:,kp1)*(-T1(:,k)+T1(:,kp1))
        if ( k==ydim )               dTy(:,k)=ccy*wz(:,km1)*(T1(:,km1)-T1(:,k))
        ! longitudinal
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 2
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 3
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           do j=4, xdim-3              ! longitudinal
              jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
              ! 3.order solution: stable unitl 84degree (dx=2.5degree, a=5e5)
              dTx(j,k)=ccx(k)*(                                                           &
&               10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&               +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&               +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&               +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&               +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           end do
           j = xdim-2
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1;
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(1.*dxlat(k)**2/kappa))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=kappa*dtdff2/dxlat(k)**2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
              j = 1
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 2
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 3
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim;
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)=ccx2*(                                                           &
&                      10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                      +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                      +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                      +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                      +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.

                end do           ! longitudinal
              j = xdim-2
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim-1
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability
                T1h=T1h+dTxh
            end do               ! additional time loop
            dTx(:,k)=T1h-T1(:,k)
        end if
    end do          ! y-loop
    dX_diffuse = wz * (dTx + dTy);

end subroutine diffusion

!+++++++++++++++++++++++++++++++++++++++
subroutine advection(T1, dX_advec,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    advection after DD

  USE mo_numerics, ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,  ONLY: pi, z_topo, uclim, vclim, ityr, z_vapor
  USE mo_physics,  ONLY: uclim_m, uclim_p, vclim_m, vclim_p
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_advec

  integer :: i
  integer, dimension(ydim):: ilat = (/(i,i=1,ydim)/)
  real, dimension(ydim) :: lat, dxlat, ccx
  real, dimension(xdim) :: T1h, dTxh
  real, dimension(xdim,ydim) :: ddx, T, dTx, dTy
  integer time2, dtdff2, tt2

  real    :: deg, dx, dy, dd, dyy, ccy, ccx2
  integer :: j, k, km1, km2, kp1, kp2, jm1, jm2, jm3, jp1, jp2, jp3

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=dt_crcl/dyy/2.
  ccx=dt_crcl/dxlat/2.

     ! latitudinal
     k=1
     kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                     vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
     end do
     k=2
     km1=k-1; kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1)))          &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) )/3. )
     end do
     do k=3, ydim-2
        km1=k-1; kp1=k+1; km2=k-2; kp2=k+2
        do j = 1, xdim
           dTy(j,k) = ccy * (                                                     &
&                       -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))        &
&                                           +wz(j,km2)*(T1(j,k)-T1(j,km2)) )      &
&                      + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))        &
&                                           +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
        end do
     end do
     k=ydim-1
     km1=k-1; kp1=k+1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) )/3.      &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1)) ) )
     end do
     k=ydim
     km1=k-1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) ) )/3.
     end do

     ! longitudinal
     do k=1, ydim
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jm1 = xdim; jm2 = xdim-1; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = 2
           jm1 = j-1; jm2 = xdim; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           do j=3, xdim-2              ! longitudinal
                jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2
                dTx(j,k)= ccx(k) * (                                                  &
&                           -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))        &
&                                               +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )      &
&                          + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))        &
&                                               +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           end do
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jp1 = j+1; jp2 = 1
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = xdim
           jm1 = j-1; jm2 = j-2; jp1 = 1; jp2 = 2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.

        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(dxlat(k)/10.0/1.))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=dtdff2/dxlat(k)/2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
                j = 1
                jm1=xdim; jm2=xdim-1; jm3=xdim-2; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 2
                jm1=j-1; jm2=xdim; jm3=xdim-1; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 3
                jm1=j-1; jm2=j-2; jm3=xdim; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)= ccx2 * (                                                          &
&                            -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )          &
&                                                 +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )          &
&                                                 +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )        &
&                           + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )          &
&                                                 +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )          &
&                                                 +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                end do           ! longitudinal
                j = xdim-2
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim-1; jp2=xdim-1; jp3=1
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim-1
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim; jp2=1; jp3=2
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=1; jp2=2; jp3=3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            &
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            &
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability
                T1h = T1h + dTxh
            end do               ! additional time loop
            dTx(:,k) = T1h - T1(:,k)
        end if
    end do          ! y-loop
    dX_advec = dTx + dTy;

end subroutine advection

!+++++++++++++++++++++++++++++++++++++++
subroutine convergence(T1, div)
!+++++++++++++++++++++++++++++++++++++++
! Calculates divergence (convergence) of a given field (i.e. spec. hum.) when omega is known
! Eq. 18 in Stassen et al 2019
  use mo_numerics, only: xdim, ydim, nstep_yr, dlon, dlat, dt_crcl
  use mo_physics,  only: ityr, rho_air, grav, pi, z_vapor, omegaclim
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1
  real, dimension(xdim,ydim), intent(out) :: div

  real    :: w
  integer :: i, j

  do j=1,ydim
    do i=1,xdim
      !< Vertical velocity omega (Pa/s) to m/s
      w = -omegaclim(i,j,ityr) / (rho_air*grav)
      !< Convergence
      div(i,j) = T1(i,j) * w * dt_crcl / z_vapor * 2.5
    end do
  end do

end subroutine convergence

!+++++++++++++++++++++++++++++++++++++++
subroutine orbital_forcing(year)
!+++++++++++++++++++++++++++++++++++++++
! meridianal changed solar radiation at top due to orbital forcing 

  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, nstep_yr, dlat
  USE mo_physics,      ONLY: jday, SW_solar, kry_start

! declare temporary fields
  integer year, kyear, day, day0, ilat
  real    So, lat, lambda
  real, dimension(ydim,nstep_yr)      ::  sw_solar0, sw_solar1

!  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, ice_cover
  real, parameter ::  pi=3.14159265359

  open(23,file='orbit')

! find current kyrs closest to current year. 
  xyr  = kry_start + (year-1)/1000.
  kyr0 = floor(xyr)
  wg1  = xyr-floor(xyr)

!  print *,'orbital forc.:', year, kry_start, xyr, kyr0, wg1
  kyear=year

  inyear = 9999
  do while (inyear .ne. kyr0)  ! Go through orbital.paras until correct year is found.
  	read(23,*) inyear, ecc, omega, epsilon
!    print *,'read orbital:',inyear, ecc, omega, epsilon
  end do
!  print *,'orbital para.:',xyr, ecc, omega, epsilon
   
  call orbital_solar(ecc, omega, epsilon)
  SW_solar0 = sw_solar
  read(23,*) inyear, ecc, omega, epsilon
  call orbital_solar(ecc, omega, epsilon)
  SW_solar1 = sw_solar

  sw_solar = (1-wg1)*SW_solar0 + wg1*sw_solar1

  close(23)
  
end subroutine orbital_forcing

!+++++++++++++++++++++++++++++++++++++++
subroutine orbital_solar(ecc, omega, epsilon)
!+++++++++++++++++++++++++++++++++++++++
! orbital forcing paramters  

  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, nstep_yr, dlat
  USE mo_physics,      ONLY: jday, SW_solar

! declare temporary fields
  integer day, day0, ilat
  real    So, lat, lambda
!  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, ice_cover
  real, parameter ::  pi=3.14159265359

!  print *,'orbital para.:',kyear, ecc, omega, epsilon
  
! add 180 degrees to omega (see lambda definition, Berger 1978 Appendix)
  omega   = omega+180;                   ! longitude of perihelion (precession angle)
!  omega   = unwrap(omega*pi/180)*180/pi;  % remove discontinuities (360 degree jumps)
  if (omega .gt. 360) omega = omega-360
  omega   = omega* pi/180;                
  epsilon = epsilon* pi/180;             

  So       = 1365.                          ! solar constant (W/m^2)
  day0     = nint((31+28+21)*360/365.24)    ! reference day for hybers code start date of year
  day0     = nint(float(nstep_yr*day0/360)) ! convert to greb time steps per year

! === Calculate insolation ===
  do iday=1,nstep_yr
     day = iday-day0
     do ilat=1, ydim
!        print*,iday, ilat
        lat    = dlat*ilat -90 -dlat/2
        rlat   = lat*pi/180                      ! latitude
!       lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
        lambda = day*2*pi/nstep_yr               ! lambda=0 for spring equinox
        delta  = asin(sin(epsilon)*sin(lambda))  !  declination of the sun
        Ho     = acos(-tan(rlat)*tan(delta))     ! hour angle at sunrise/sunset
!       no sunrise or no sunset: Berger 1978 eqn (8),(9)
        if (abs(rlat) .ge. pi/2 - abs(delta)  .and. rlat*delta .gt. 0  ) Ho = pi
        if (abs(rlat) .ge. pi/2 - abs(delta)  .and. rlat*delta .le. 0  ) Ho = 0
!       Insolation: Berger 1978 eq (10)
        SW_solar(ilat,iday) =  So/pi*(1+ecc*cos(lambda-omega))**2/(1-ecc**2)**2         &
&                            *( Ho*sin(rlat)*sin(delta) + cos(rlat)*cos(delta)*sin(Ho))

     end do
  end do
  
end subroutine orbital_solar

!TB
!+++++++++++++++++++++++++++++++++++++++
function gmean(data)
!+++++++++++++++++++++++++++++++++++++++
! global mean for data

use mo_numerics,		ONLY: xdim, ydim, dlat

! declare variables
real*4, dimension(xdim,ydim) 	:: data, w
real*4, dimension(ydim)		    :: lat
real*4                          :: gmean

do i=1,ydim
lat(i) = -90-(dlat*0.5)+i*dlat
end do
do i=1,xdim
w(i,:) = cos(2.*3.14159*lat/360.)
end do

gmean = sum(data*w)/sum(w)

end function

