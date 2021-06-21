#!/bin/csh
# GREB with ice-sheet restart job-script for standalone experiment
# author: Zhiang Xie 
#

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

# exp end time of orbital forcing [kyrs before today]
set KYREND     = 0

# acceleration of orbital forcing 1 -> normal 10 -> 10x faster ...
set DTACC      = 1

# length of single GREB model run [yrs]
set TDIM_RUN   = 1000

# number of GREB model run in one job (per restart/newstart)
set NRUN       = 50
if ($EXP == 311) set NRUN = 50

# frequency of grep output in [yrs of simulation]
set DT_OUT     = $TDIM_RUN

# frequency of restart output
set DT_RSTRT   = $TDIM_RUN

# ICE SHEET 
set log_ice_sheet = 1  # force vertical velocity omega with external file (0=no forcing; 1=forcing)
set log_ice_ithk  = 1  # ice thickness initial condition (0=ice free; 1=start from current ice thickness)
set log_ice_topo  = 1  # topogragh coupling control
set log_ice_slv   = 1  # sea level control
set log_ice_ant   = 1
# directories
set WDIR=/scratch/w48/zx2351/ice-sheet-development
set MDIR=$WDIR
set IDIR=$WDIR/input/
set ODIR=./


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

@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 

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
rm -f control scenario
#date
setenv OMP_NUM_THREADS 3
setenv KMP_AFFINITY verbose,none
unlimit stacksize

./greb.x > greb.out.txt

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
zdef   1 linear 1 1
tdef $TDIM linear 15jan0  1mo
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
#
@ xyr = $xyr +  $TDIM_RUN * $DTACC / 1000 

@ loop ++
end
# ------------ end run loop ----------------
# cd $WDIR
# if(${xyr}<0) csh ./run.greb.icealone.pictrl.restart.csh

exit
