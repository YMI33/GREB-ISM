#!/bin/csh
# GREB with ice-sheet newstart job-script for standalone experiment
# author: Zhiang Xie, Dietmar Dommenget
#

set COMPILE  = 1

# ----------------------------------------------------------
#      INPUT
# ----------------------------------------------------------
# GREB experiment set up exp = 310 -> dynamic equalibrium run for standalone experiment
# GREB experiment set up exp = 311 -> dynamic equalibrium run for transition experiment
set EXP        = 310 # ice sheet standalone run 

# name of experiment
set NAME       = benchmark_pictrl
# set NAME = benchmark_tran

# exp start time of orbital forcing [kyrs before today]
set KYRSTART   = -200
if ($EXP == 311) set KYRSTART = -250

# exp end time of orbital forcing [kyrs before today]
set KYREND     = 0

# acceleration of orbital forcing 1 -> normal 10 -> 10x faster ...
set DTACC      = 1

# length of qflux correction run [yrs]
set TIME_QFLX  = 15

# length of control run [yrs]
set TIME_CTRL  = 3

# length of single GREB model run [yrs]
set TDIM_RUN   = 1000

# number of GREB model run in one job (per restart/newstart)
set NRUN       = 1

# frequency of grep output in [yrs of simulation]
set DT_OUT     = $TDIM_RUN

# frequency of restart output
set DT_RSTRT   = $TDIM_RUN

# ICE SHEET 
set log_ice_ithk  = 1  # ice thickness initial condition (0=ice free; 1=start from current ice thickness)
set log_ice_topo  = 1  # topogragh coupling control
set log_ice_slv   = 1  # sea level control
set log_ice_ant   = 1

# directories
set WDIR=/Users/zxie0012/Documents/GREB-ISM
set JDIR=$WDIR/job-script
set MDIR=$WDIR/src
set IDIR=$WDIR/input/

mkdir  $WDIR/experiments/$NAME
mkdir  $WDIR/experiments/$NAME/work
mkdir  $WDIR/experiments/$NAME/restarts

#-------------------------------------------------------
# compile GREB model
#-------------------------------------------------------

if ($COMPILE == 1) then
echo ''
echo 'compile GREB model'
cd $WDIR/experiments/$NAME/work/
set ICE   = $MDIR/ice-sheet.f90
set MAIN  = $MDIR/greb.main.f90
set OCEN  = $MDIR/greb.ocean.f90
set ATMO  = $MDIR/greb.atmosphere.f90
\rm -f greb.x *.mod

gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops $MAIN $ATMO $OCEN $ICE -o greb.x

endif 


#-------------------------------------------------------
# run GREB model
#-------------------------------------------------------
echo ''
echo 'run GREB model'
echo ''
cd $WDIR/experiments/$NAME/work/

#--------------------------------------------
# DEFAULT INPUT DATA
set TCLIM=${IDIR}'/erainterim.tsurf.1979-2015.clim.bin'
set MASK=${IDIR}'/global.topography.t31.gad'
set QCLIM=${IDIR}'/erainterim.atmospheric_humidity.clim.bin'
set MOIST=${IDIR}'/ncep.soil_moisture.clim.bin'
set SOLAR=${IDIR}'/solar_radiation.clim.bin'
set UWIND=${IDIR}'/erainterim.zonal_wind.850hpa.clim.bin'
set VWIND=${IDIR}'/erainterim.meridional_wind.850hpa.clim.bin'
set MLD=${IDIR}'/woce.ocean_mixed_layer_depth.clim.bin'
set CLOUD=${IDIR}'/isccp.cloud_cover.clim.bin'
set TOCLM=${IDIR}'/Tocean.clim.bin'
set AWIND=${IDIR}'/erainterim.windspeed.850hpa.clim.bin'
set OMCLM=${IDIR}'/erainterim.omega.vertmean.clim.bin'
set OMSTD=${IDIR}'/erainterim.omega_std.vertmean.clim.bin'
set BROCK=${IDIR}'/bedmachine.bed.rock.bin'
set ICCLM=${IDIR}'/ice.height.first-guess.clim730.bin'
set ORBIT=${IDIR}'/orbital.parameters.last5mill.yrs.nocomments.txt'
set PRECI=${IDIR}'/precip.NCEP-DOE.730clim.gad'

set TS_G='nodts_greenland'
set TS_A='nodts_antarctica'
touch nodts_greenland
touch nodts_antarctica
if ( $EXP == 311 ) set TS_G=${IDIR}'/icesheet_input/grip_temp_nosmooth_250kyr.dat'
if ( $EXP == 311 ) set TS_A=${IDIR}'/icesheet_input/vostok_temp_nosmooth_250kyr.dat'

rm -f tclim ztopo qclim moist solar uwind vwind mldclim cloud orbit toclim abswind
rm -f omclim omstdv iceclm bedrock test.gad precip delta_ts_Greenland delta_ts_Antarctica

ln -s $TCLIM tclim
ln -s $MASK  ztopo
ln -s $QCLIM qclim
ln -s $MOIST moist
ln -s $SOLAR solar
ln -s $UWIND uwind
ln -s $VWIND vwind
ln -s $MLD   mldclim
ln -s $CLOUD cloud
ln -s $TOCLM toclim
ln -s $AWIND abswind
ln -s $OMCLM omclim
ln -s $OMSTD omstdv
ln -s $BROCK bedrock
ln -s $ICCLM iceclm
ln -s $ORBIT orbit
ln -s $PRECI precip
ln -s $TS_G delta_ts_Greenland
ln -s $TS_A delta_ts_Antarctica


rm -f ToF_correct qF_correct TF_correct

if ( $TIME_QFLX == 0) then
ln -s  ../ToF_correct.$QFILES  ToF_correct
ln -s  ../qF_correct.$QFILES   qF_correct
ln -s  ../TF_correct.$QFILES   TF_correct
ln -s  ../Ta_ini.$QFILES       Ta_ini
endif


# ------------   newstart   ----------------
echo ''
echo ' new start'
echo ''

set RESTART  = 0
set xyr = $KYRSTART
set TTYR = 0
@ TTYR = -1 * $KYRSTART  
@ TTYR = $TTYR + $xyr  
set BASE = greb.exp-$EXP.$NAME.$TTYR
if ($TTYR < 10000) set BASE = greb.exp-$EXP.$NAME.0$TTYR
if ($TTYR <  1000) set BASE = greb.exp-$EXP.$NAME.00$TTYR
if ($TTYR <   100) set BASE = greb.exp-$EXP.$NAME.000$TTYR
if ($TTYR <    10) set BASE = greb.exp-$EXP.$NAME.0000$TTYR

#  namelist
cat >namelist <<EOF
&NUMERICS
time_flux   = $TIME_QFLX
time_ctrl   = $TIME_CTRL
time_scnr   = $TDIM_RUN
log_restart = $RESTART
dt_rstrt    = $DT_RSTRT
dt_out      = $DT_OUT
/
&PHYSICS
 log_exp    = $EXP      ! complete GREB model; 2xCO2 forcing
 dtyr_acc   = $DTACC    ! acceleration in sw_solar forcing
 kry_start  = $KYRSTART ! historical start date [kyrs]
 kry_end    = $KYREND   ! historical end date [kyrs]
 Todeep0    =275.15     ! deep ocen temp
 log_ice_ithk=$log_ice_ithk   ! ice thickness initial condition (0=ice free; 1=start from current ice thickness)
 log_ice_topo=$log_ice_topo ! topograph control
 log_ice_slv =$log_ice_slv  ! sea level control
 log_ice_ant =$log_ice_ant  ! antarctica switch
/
EOF
rm -f control scenario

setenv OMP_NUM_THREADS 3
setenv KMP_AFFINITY verbose,none
unlimit stacksize

./greb.x 

cat >job.restart.txt <<EOF
 xyr       $xyr
 NRUN      $NRUN
 TDIM_RUN  $TDIM_RUN
 DTACC     $DTACC
EOF
cp namelist     ../greb.exp-${EXP}.${NAME}.namelist.txt
mv scenario.bin     ../${BASE}.scenario.bin
mv control.bin      ../${BASE}.control.bin
mv scenario.gmean.bin ../${BASE}.scenario.gmean.bin
cp restart.txt restart_in.txt
cp restart.bin restart_in.bin
mv restart.txt  ../restarts/${BASE}.restart.txt
mv restart.bin  ../restarts/${BASE}.restart.bin
set RESTART    = 1
#
cd ../
@ TDIM = 12 * $TIME_CTRL
cat > ${BASE}.control.ctl <<EOF
dset ^${BASE}.control.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   4 levels 1 0.57 -0.57 -1
tdef $TDIM linear 15jan0000 1mo
vars 16
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
mask   1 0 land-sea mask
precip 1 0 precipitation
albd   1 0 surface albedo
glcier 1 0 ice surface temperature
iceh 1 0 ice thickness
zs 1 0 ice surface height
mass 1 0 mass balance
adv 1 0 advection term
calv 1 0 calving
vx 1 0 ice flow zonal velocity 
vy 1 0 ice flow meridianal velocity
tice 4 0 ice temperature in different layers
endvars
EOF
#
@ TDIM = 12 * $TDIM_RUN
cat > ${BASE}.scenario.ctl <<EOF
dset ^${BASE}.scenario.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   4 levels 1 0.57 -0.57 -1
tdef $TDIM linear 15jan0000 1mo
vars 16
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
mask   1 0 land-sea mask
precip 1 0 precipitation
albd   1 0 surface albedo
glcier 1 0 ice surface temperature
iceh   1 0 ice thickness
zs 1 0 ice surface height
mass 1 0 mass balance
adv 1 0 advection term
calv 1 0 calving
vx 1 0 ice flow zonal velocity 
vy 1 0 ice flow meridianal velocity
tice 4 0 ice temperature in different layers
endvars
EOF
#
cat > ${BASE}.scenario.gmean.ctl <<EOF
dset ^${BASE}.scenario.gmean.bin
undef 9.e27
xdef  1 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  1 linear 1 1
tdef  $TDIM linear 15jan0001 1mo
vars 7 
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
albd   1 0 albdo
precip 1 0 precipitation
slv    1 0 sea level 
endvars 
EOF

@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 

# overall ctl file

@ TDIM = 12 * $TDIM_RUN
cat > greb.exp-$EXP.$NAME.scenario.ctl <<EOF
dset ^greb.exp-$EXP.$NAME.0%y4.scenario.bin
undef 9.e27
options template
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef  4  levels 1 0.57 -0.57 -1
tdef $TDIM linear 15jan0001 1mo
vars 16
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
mask   1 0 land-sea mask
precip 1 0 precipitation
albd   1 0 surface albedo
glcier  1 0 ice surface temperature
iceh   1 0 ice thickness
zs 1 0 ice surface height
mass 1 0 mass balance
adv 1 0 advection term
calv 1 0 calving
vx 1 0 ice flow zonal velocity 
vy 1 0 ice flow meridianal velocity
tice 4 0 ice temperature in different layers
endvars
EOF
#
cat > greb.exp-$EXP.$NAME.scenario.gmean.ctl <<EOF
dset ^greb.exp-$EXP.$NAME.0%y4.scenario.gmean.bin
undef 9.e27
options template
xdef  1 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  1 linear 1 1
tdef  $TDIM linear 15jan0001 1mo
vars 7 
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
albd   1 0 albdo
precip 1 0 precipitation
slv    1 0 sea level 
endvars 
EOF
# ------------ restart loop ----------------
echo ''
echo ' restart loop'
echo ''


set RESTART = 1
set loop    = 1
while ( $loop <= $NRUN )
echo 'run:' $loop ' ' $xyr'kyr'
cd $WDIR/experiments/$NAME/work/
set TTYR = 0
@ TTYR = -1 * $KYRSTART  
@ TTYR = $TTYR + $xyr  
set BASE = greb.exp-$EXP.$NAME.$TTYR
if ($TTYR < 10000) set BASE = greb.exp-$EXP.$NAME.0$TTYR
if ($TTYR <  1000) set BASE = greb.exp-$EXP.$NAME.00$TTYR
if ($TTYR <   100) set BASE = greb.exp-$EXP.$NAME.000$TTYR
if ($TTYR <    10) set BASE = greb.exp-$EXP.$NAME.0000$TTYR
#  namelist
cat >namelist <<EOF
&NUMERICS
time_flux   = $TIME_QFLX
time_ctrl   = $TIME_CTRL
time_scnr   = $TDIM_RUN
log_restart = $RESTART
dt_rstrt    = $DT_RSTRT
dt_out      = $DT_OUT
/
&PHYSICS
 log_exp    = $EXP      ! complete GREB model; 2xCO2 forcing
 dtyr_acc   = $DTACC    ! acceleration in sw_solar forcing
 kry_start  = $KYRSTART ! historical start date [kyrs]
 kry_end    = $KYREND   ! historical end date [kyrs]
 Todeep0    =275.15     ! deep ocen temp
 log_ice_ithk=$log_ice_ithk   ! ice thickness initial condition (0=ice free; 1=start from current ice thickness)
 log_ice_topo=$log_ice_topo ! topograph control
 log_ice_slv =$log_ice_slv  ! sea level control
 log_ice_ant =$log_ice_ant  ! antarctica switch
/
EOF
rm -f control scenario
#date
./greb.x > greb.out.txt
#date
cat >job.restart.txt <<EOF
 xyr       $xyr
 NRUN      $NRUN
 TDIM_RUN  $TDIM_RUN
 DTACC     $DTACC
EOF
mv scenario.bin     ../${BASE}.scenario.bin
mv scenario.gmean.bin ../${BASE}.scenario.gmean.bin
cp restart.txt restart_in.txt
cp restart.bin restart_in.bin
mv restart.txt  ../restarts/${BASE}.restart.txt
mv restart.bin  ../restarts/${BASE}.restart.bin
#
cd ../
#
@ TDIM = 12 * $TDIM_RUN
cat > ${BASE}.scenario.ctl <<EOF
dset ^${BASE}.scenario.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   4 levels 1 0.57 -0.57 -1
tdef $TDIM linear 15jan0001 1mo
vars 16
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
mask   1 0 land-sea mask
precip 1 0 precipitation
albd   1 0 surface albedo
glcier  1 0 ice surface temperature
iceh   1 0 ice thickness
zs 1 0 ice surface height
mass 1 0 mass balance
adv 1 0 advection term
calv 1 0 calving
vx 1 0 ice flow zonal velocity 
vy 1 0 ice flow meridianal velocity
tice 4 0 ice temperature in different layers
endvars
EOF
#
cat > ${BASE}.scenario.gmean.ctl <<EOF
dset ^${BASE}.scenario.gmean.bin
undef 9.e27
xdef  1 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef   4 levels 1 0.57 -0.57 -1
tdef  $TDIM linear 15jan0001 1mo
vars 7 
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
albd   1 0 albdo
precip 1 0 precipitation
slv    1 0 sea level 
endvars 
EOF
#
@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 

@ loop ++
end
# ------------ end run loop ----------------
date

# cd $WDIR
# qsub ./run.greb.icealone.pictrl.restart.csh

exit
