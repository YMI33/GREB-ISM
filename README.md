# GREB-ISM
Code for GREB-ISM

## File systems
* input: binary files for model input
* job-script: csh script for running model
* src: source code for model
* experiment: directory for storing experiment results

## Model running instruction
1. edit *GREB-ISM/job-script/greb-ice-sheet.exp-scenario.newstart.com*
2. change environment variable **WDIR** to the absolute path of your GREB-ISM 
3. change experiment setting (EXP:experiment number, KYRSTART:start date of experiment, etc.)
4. run *GREB-ISM/job-script/greb-ice-sheet.exp-scenario.newstart.com*
5. after the job finished, a restart file *greb-ice-sheet.${NAME}.restart.com* will be generated and you can continue the experiment by rerunning the restart file.
6. you can access the experiment output in *GREB-ISM/experiment*, which are binary files in [GrADS format](http://cola.gmu.edu/grads/gadoc/gadocindex.htmli). You can easily read the data by CTL files in the same directory by GrADS. 


