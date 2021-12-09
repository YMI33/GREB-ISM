#!/bin/csh
# GREB-ISM job-script for EISMINT I/II experiments
# to use this script, the user needs to set $WDIR to their own directory
#
# author: Zhiang Xie 
#  
# create work directory if does not already exist
set NAME = EISMINT

set WDIR=/Users/zxie0012/Documents/GREB-ISM
set MDIR=$WDIR/src
set ODIR=$WDIR/experiments/$NAME

mkdir  $WDIR/experiments/$NAME
mkdir  $WDIR/experiments/$NAME/work
mkdir  $WDIR/experiments/$NAME/restarts

#####################
# BEGIN USER INPUT! #
#####################
set i=1
# length of sensitivity experiment in years
set YEARS=200000
# output time step
set TIMEOUT=1000
while($i<7)

# settings for scenario
# log_decouple = 1 -> EISMINT I
# log_decouple = 0 -> EISMINT II
set log_decouple   = 1 # 0 for dynamic and thermodynamic coupled experiment, 1 for decoupled

# scenario number from list below
# Exp. Value    EISMINT I         EISMINT II
#     1         Fix. Steady       Exp. A
#     2         Mov. Steady       Exp. B
#     3         Fix. Tran. (20k)  Exp. C
#     4         Mov. Tran. (20k)  Exp. D
#     5         Fix. Tran. (40k)    / 
#     6         Fix. Tran. (40k)    /

set EXP=$i

# if scenario is forced climate change (EXP 230) or forced ENSO (EXP 240 or 241)
set log_restart_w=0 # control for restart file writing
set log_restart_r=0 # control for restart file reading
set log_bound_cond = 0 # 0 for fixed boundary; 1 for moving boundary
set log_time_depen = 0 # 0 for steady run; 1 for time dependent run 
set period   = 40000 # period for sin function
set time_stp = 10 # time step for time advance [a]

if( ${log_decouple} == 1 ) then
    if( $EXP == 1 ) then
    set log_bound_cond = 0 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 0 # 0 for steady run; 1 for time dependent run 
    endif
    
    if( $EXP == 2 ) then
    set log_bound_cond = 1 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 0 # 0 for steady run; 1 for time dependent run 
    endif
    
    if( $EXP == 3 ) then
    set log_bound_cond = 0 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 1 # 0 for steady run; 1 for time dependent run 
    set period = 20000 # period for sin function
    endif
    
    if( $EXP == 4 ) then
    set log_bound_cond = 1 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 1 # 0 for steady run; 1 for time dependent run 
    set period = 20000 # period for sin function
    endif
    
    if( $EXP == 5 ) then
    set log_bound_cond = 0 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 1 # 0 for steady run; 1 for time dependent run 
    set period = 40000 # period for sin function
    endif
    
    if( $EXP == 6 ) then
    set log_bound_cond = 1 # 0 for fixed boundary; 1 for moving boundary
    set log_time_depen = 1 # 0 for steady run; 1 for time dependent run 
    set period = 40000 # period for sin function
    endif
endif

if( ${log_decouple} == 0 ) then
    set log_bound_cond = $EXP
endif

# for EXP='100', give here the name of input CO2 forcing data set without '.txt'
set CO2input=none

### compile GREB model (uncomment one of these three options)
### gfortran compiler (Linux (e.g. Ubuntu), Unix or MacBook Air)
gfortran -fopenmp -march=native -O3 -ffast-math -funroll-loops -fdefault-real-8 ${MDIR}/model_verify.f90 -o greb.x

###################
# END USER INPUT! #
###################

setenv OMP_NUM_THREADS 6
setenv KMP_AFFINITY verbose,none
unlimit stacksize

# move complied files to work directory
mv greb.x $ODIR/work/
mv *.mod $ODIR/work/

if( $log_restart_r == 1 ) then
cp $ODIR/restarts/restart.bin ${ODIR}/work/
cp $ODIR/restarts/iteration_info ${ODIR}/work/ 
endif

# change to work directory
cd ${ODIR}/work/


#  generate namelist
cat >namelist <<EOF
&NUMERICS
time_flux = 3  		! length of flux corrections run [yrs]
time_ctrl = 3 		! length of control run [yrs]
time_stp  = $time_stp   ! time step for time advance [yrs]
time_scnr = $YEARS  	! length of scenario run [yrs]
time_out  = $TIMEOUT    ! length of time output [yrs]
T         = $period     ! period of sin function
/
&PHYSICS
 log_exp = $EXP 	! sensitivity run as set above
 log_bound_cond = $log_bound_cond
 log_time_depen = $log_time_depen
 log_decouple  = $log_decouple
/
EOF

echo 'experiment: '$EXP
# run model
time ./greb.x > EISMINT_run_log

# postprocessing
# rename scenario run output and move it to output folder

# calculate months of scenario run for header file
@ TIMES = $YEARS / $TIMEOUT 
#

if( ${log_decouple} == 1 ) set FILENAME=EISMINT_I_exp$EXP
if( ${log_decouple} == 0 ) set FILENAME=EISMINT_II_exp${EXP}

# rename scenario run output and move it to output folder
mv scenario.bin $WDIR/experiments/$NAME/${FILENAME}.bin
if( $log_restart_w == 1 ) then
mv restart.bin $ODIR/restarts/
mv iteration_info $ODIR/restarts/
endif
# scenario run
cat >${ODIR}/${FILENAME}.ctl <<EOF
dset ^${FILENAME}.bin
undef 9.e27
xdef  96 linear 0 3.75
ydef  144 linear -89.375 1.25
zdef   4 levels 1 0.5774 -0.5774 -1
tdef $TIMES linear 15jan0001  ${time_stp}yr
vars 4
h 1 0 ice thickness
temp 4 0 temperature
vy 4 0 ice flow meridianal velocity
vz 4 0 vertical velocity
endvars
EOF
@ i++
cd .. 
end


