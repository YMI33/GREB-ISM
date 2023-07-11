#!/bin/csh
#
# author: Zhiang Xie 
#
# GREB with ice-sheet newstart job-script for NCI GADI

# load module for compiler
module load intel/2018.4
module load mpi/intel/2018.4

set COMPILE  = 1

# ----------------------------------------------------------
#      INPUT
# ----------------------------------------------------------


# GREB experiment set up exp = 200 -> dynamic equalibrium run 
# GREB experiment set up exp = 201 -> reduce solar radiation run 
# GREB experiment set up exp = 202 -> 95% radiation control to 100% radiation
# GREB experiment set up exp = 203 -> 280 ppm control to 340 ppm
# GREB experiment set up exp = 205 -> 20kyrs SW
# GREB experiment set up exp = 206 -> 20kyrs SW
# GREB experiment set up exp = 207 -> 20kyrs SW
set EXP        = 311
set CO2        = 280 

# different scenario run with end of transition run
set NAME = test 

# LGM and PRE index setting
# GRIP
# set LGM = -13.99
# set PRE = -0.89
# 50%
set LGM = -6.26
set PRE = 0.0


# scenario setting
set HEAT = 1
set TOPO = 1
set SLV  = 1
set ALBD = 1
set PREP = 1

# physical options
set LOG_VDIF  = 1	# vapour diffusion
set LOG_VADV  = 1	# vapour advection
set LOG_HDIF  = 1	# heat diffusion
set LOG_HADV  = 1	# heat advection

# initial condition setting
set INI_CTRL   = 0 # set scenario start from control results (1) or initial condition (0)
set INI_FILE   = 0 # set scenario start from initial file (1) or input file (0)

# exp start time of orbital forcing [kyrs before today]
set KYRSTART   = -250

# exp end time of orbital forcing [kyrs before today]
set KYREND     = 0

# acceleration of orbital forcing 1 -> normal 10 -> 10x faster ...
set DTACC      = 1

# length of qflux correction run [yrs]
set TIME_QFLX  = 20

# length of control run [yrs]
# set TIME_CTRL  = 30000
set TIME_CTRL  = 10

# length of single GREB model run [yrs]
set TDIM_RUN   = 1000

# number of GREB model run in one job (per restart/newstart)
set NRUN       = 1

# frequency of grep output in [yrs of simulation]
set DT_OUT     = $TDIM_RUN

# frequency of restart output
set DT_RSTRT   = $TDIM_RUN

# directories
# set WDIR=/work/ess-xieza/GREB-ISM/
set WDIR=/Users/zxie0012/Documents/GREB-ISM
set MDIR=$WDIR/
set JDIR=$WDIR/job-script/
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
set ICE   = $MDIR/src/ice-sheet.f90
set MAIN  = $MDIR/src/greb.main.f90
set OCEN  = $MDIR/src/greb.ocean.f90
set ATMO  = $MDIR/src/greb.atmosphere.f90
\rm -f greb.x *.mod
# NCI RAIJIN
# scott wales RAIJIN:
#ifort -O3 -xHost -g -traceback -qopt-report -fp-model fast  $MODEL $SHELL -o ./greb.x

# GADI (Zhiang):
# gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops $MODEL $SHELL -o ./greb.x
# -fopenmp NOT WORKING FOR ME ON GADI: mpirun greb.x?

# gadi 15.04.20:
# gfortran -march=native -ffast-math -O3 -funroll-loops $MODEL $SHELL $ICE -o ./greb.x

# Ang's code 29.04.20:
# gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops -fdefault-real-8  $MODEL $SHELL $ICE  -o greb.x
gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops $MAIN $ATMO $OCEN $ICE -o greb.x
#
endif 

unlimit stacksize

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
#set ICCLM=${IDIR}'/ice.height.steady.state.clim730.bin'
set ICCLM=${IDIR}'/ice.height.first-guess.clim730.bin'
set ORBIT=${IDIR}'/orbital.parameters.last5mill.yrs.nocomments.txt'
set PRECI=${IDIR}'/precip.NCEP-DOE.730clim.gad'
set LGM_PRECI=${IDIR}'/precip.AWI-ESM_LGM.730clim.gad'
set LGM_TCLIM=${IDIR}'/ts.AWI-ESM_LGM.730clim.gad'
set INITF=/scratch/w40/zx2351/ice-sheet-development/experiments/sny_ctrl-clim/greb.exp-200.sny_ctrl-clim.00030.scenario.bin

set TS_G='nodts_greenland'
set TS_A='nodts_antarctica'
touch nodts_greenland
touch nodts_antarctica
if ( $EXP == 311 ) set TS_G=${IDIR}'/icesheet_input/ts_2myrs_50_250k.dat'
if ( $EXP == 311 ) set TS_A=${IDIR}'/icesheet_input/ts_2myrs_50_250k.dat'
if ( $EXP == 311 ) set FCO2=${IDIR}'/CO2/co2_2myrs_250k_Berends.dat'

rm -f tclim ztopo qclim moist solar uwind vwind mldclim cloud orbit toclim abswind lgm_precip lgm_tsurf init_file
rm -f omclim omstdv iceclm bedrock precip delta_ts_Greenland delta_ts_Antarctica co2forcing

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
ln -s $INITF init_file
ln -s $TS_G delta_ts_Greenland
ln -s $TS_A delta_ts_Antarctica
ln -s $LGM_PRECI lgm_precip
ln -s $LGM_TCLIM lgm_tsurf
ln -s $FCO2 co2forcing 

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
@ xyr = $KYRSTART + 1
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
log_ctrlini = $INI_CTRL
log_fileini = $INI_FILE
dt_rstrt    = $DT_RSTRT
dt_out      = $DT_OUT
/
&PHYSICS
 log_exp     = $EXP         ! complete GREB model; 2xCO2 forcing
 dtyr_acc    = $DTACC       ! acceleration in sw_solar forcing
 kry_start   = $KYRSTART    ! historical start date [kyrs]
 kry_end     = $KYREND      ! historical end date [kyrs]
 Todeep0     = 275.15       ! deep ocen temp
 co2_scn     = $CO2         ! co2 concentration [ppm]
 log_ice_slv = $SLV         ! sea level option
 log_ice_topo= $TOPO        ! lapse rate topograph option
 log_ice_heat= $HEAT        ! fix climatology ice thickness
 log_ice_albd= $ALBD        ! fix albd
 log_ice_prep= $PREP        ! fix precipitation
 log_hdif    = $LOG_HDIF
 log_hadv  	 = $LOG_HADV
 log_vdif  	 = $LOG_VDIF
 log_vadv    = $LOG_VADV
 ind_lgm     = $LGM         ! LGM temperature index (-21 ka)
 ind_pre     = $PRE         ! today temperature index (0 ka)
/
EOF
rm -f control scenario

setenv OMP_NUM_THREADS 3
setenv KMP_AFFINITY verbose,none

date
./greb.x > greb_control.txt
#date
cat >job.restart.txt <<EOF
 xyr       $xyr
 NRUN      $NRUN
 TDIM_RUN  $TDIM_RUN
 DTACC     $DTACC
EOF
cp namelist     ../greb.exp-${EXP}.${NAME}.namelist.txt
mv control.bin      ../${BASE}.control.bin
mv control.gmean.bin      ../${BASE}.control.gmean.bin
mv scenario.bin     ../${BASE}.scenario.bin
mv scenario.gmean.bin ../${BASE}.scenario.gmean.bin
cp restart.txt restart_in.txt
cp restart.bin restart_in.bin
mv restart.txt  ../restarts/${BASE}.restart.txt
mv restart.bin  ../restarts/${BASE}.restart.bin
set RESTART    = 1
#
cd ../
@ TDIM = 12 * $TIME_CTRL / 1000
cat > ${BASE}.control.ctl <<EOF
dset ^${BASE}.control.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   1 levels 1 
tdef 12000 linear 15jan0001 1mo
vars 19
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
sw  1 0 shortwave radiation [w m-2]
lw  1 0 longwave radiation [w m-2]
sh  1 0 sensible heat [w m-2]
lh  1 0 latent heat [w m-2]
ih  1 0 ice latent heat [w m-2]
pdd 1 0 positive degree days 
endvars
EOF

cat > ${BASE}.control.gmean.ctl <<EOF
dset ^${BASE}.control.gmean.bin
undef 9.e27
xdef  1 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  1 linear 1 1
tdef 12000 linear 15jan0001 1mo
vars 7 
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
sw     1 0 shortwave raditaion 
precip 1 0 precipitation
slv    1 0 sea level 
endvars 
EOF
# overall ctl files
cat > greb.exp-$EXP.$NAME.scenario.ctl <<EOF
dset ^greb.exp-$EXP.$NAME.0%y4.scenario.bin
options template
undef 9.e27
xdef  96 linear 0 3.75
ydef  48 linear -88.125 3.75
zdef   1 levels 1 
tdef 12000 linear 15jan0001 1mo
vars 19
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
mask   1 0 land-sea mask
precip 1 0 precipitation
albd   1 0 surface albedo
glcier 1 0 glacier type 
iceh 1 0 ice thickness
zs 1 0 ice surface height
mass 1 0 mass balance
adv 1 0 advection term
calv 1 0 calving
sw  1 0 shortwave radiation [w m-2]
lw  1 0 longwave radiation [w m-2]
sh  1 0 sensible heat [w m-2]
lh  1 0 latent heat [w m-2]
ih  1 0 ice latent heat [w m-2]
pdd 1 0 positive degree days 
endvars
EOF
#
cat > greb.exp-$EXP.$NAME.scenario.gmean.ctl <<EOF
dset ^greb.exp-$EXP.$NAME.0%y4.scenario.gmean.bin
options template
undef 9.e27
xdef  1 linear 0 3.75
ydef  1 linear -88.125 3.75
zdef  1 linear 1 1
tdef 12000 linear 15jan0001 1mo
vars 7 
tsurf  1 0 surface temperature
tatmos 1 0 atmosphere temperature
tocean 1 0 ocean temperature
vapor  1 0 water vapor
sw     1 0 shortwave radiation
precip 1 0 precipitation
slv    1 0 sea level 
endvars 
EOF

@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 
# ------------ restart loop ----------------
echo ''
echo ' restart loop'
echo ''
cd $WDIR/experiments/$NAME/work/
set RESTART = 1
#  namelist
cat >namelist <<EOF
&NUMERICS
time_flux   = $TIME_QFLX
time_ctrl   = $TIME_CTRL
time_scnr   = $TDIM_RUN
log_restart = $RESTART
log_fileini = $INI_FILE
log_ctrlini = $INI_CTRL
dt_rstrt    = $DT_RSTRT
dt_out      = $DT_OUT
/
&PHYSICS
 log_exp     = $EXP         ! complete GREB model; 2xCO2 forcing
 dtyr_acc    = $DTACC       ! acceleration in sw_solar forcing
 kry_start   = $KYRSTART    ! historical start date [kyrs]
 kry_end     = $KYREND      ! historical end date [kyrs]
 Todeep0     = 275.15       ! deep ocen temp
 co2_scn     = $CO2         ! co2 concentration [ppm]
 log_ice_slv = $SLV         ! sea level option
 log_ice_topo= $TOPO        ! lapse rate topograph option
 log_ice_heat= $HEAT        ! fix climatology ice thickness
 log_ice_albd= $ALBD        ! fix albd
 log_ice_prep= $PREP        ! fix precipitation
 log_hdif    = $LOG_HDIF
 log_hadv  	 = $LOG_HADV
 log_vdif  	 = $LOG_VDIF
 log_vadv    = $LOG_VADV
 ind_lgm     = $LGM         ! LGM temperature index (-21 ka)
 ind_pre     = $PRE         ! today temperature index (0 ka)
/
EOF
rm -f control scenario

cd ${JDIR}

# create restart file
# custom setting part for restart file
cat > greb-ice-sheet.${NAME}.restart.com <<EOF
#!/bin/csh
#
# author: Zhiang Xie 
#
# GREB with ice-sheet restart job-script for NCI GADI

# ----------------------------------------------------------
#      INPUT
# ----------------------------------------------------------

# name of experiment
set CO2  = $CO2
set NAME = $NAME
 
set EXP        = $EXP

set NRUN       = $NRUN

# exp start time of orbital forcing [kyrs before today]
set KYRSTART   = $KYRSTART

# directories
set WDIR=$WDIR
set MDIR=$MDIR
set IDIR=$IDIR
set JDIR=$JDIR

EOF

# running seting for restart file
cat >> greb-ice-sheet.${NAME}.restart.com <<"EOF"

unlimit stacksize

#-------------------------------------------------------
# run GREB model
#-------------------------------------------------------
echo ''
echo 'run GREB model'
echo ''
cd $WDIR/experiments/$NAME/work/

# ------------ restart loop ----------------
echo ''
echo ' restart loop'
echo ''
cd $WDIR/experiments/$NAME/work/
set XDUM=`grep -i 'xyr' job.restart.txt`
set xyr = $XDUM[2]
set XDUM=`grep -i 'NRUN' job.restart.txt`
set NRUN = $XDUM[2]
set XDUM=`grep -i 'TDIM_RUN' job.restart.txt`
set TDIM_RUN = $XDUM[2]
set XDUM=`grep -i 'DTACC' job.restart.txt`
set DTACC = $XDUM[2]
#
@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 
#
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
#
rm -f control scenario
setenv OMP_NUM_THREADS 3
setenv KMP_AFFINITY verbose,none
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
#
@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 

@ loop ++
end
# ------------ end run loop ----------------
date

# resubmit job
#
cd ${JDIR}
#if ($xyr<0) bsub < run_greb.${NAME}.bash 

"EOF"

cat > run_greb.${NAME}.bash <<EOF
#!/bin/bash
#BSUB -J ${NAME} 
#BSUB -q ser 
#BSUB -n 3
#BSUB -R "span[ptile=3]" 
#BSUB -e ./elog/error.log 
#BSUB -o ./elog/run_log 
#
csh ./greb-ice-sheet.${NAME}.restart.com

EOF

chmod +x greb-ice-sheet.${NAME}.restart.com 
chmod +x run_greb.${NAME}.bash 
./run_greb.${NAME}.bash 

