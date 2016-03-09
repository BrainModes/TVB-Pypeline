#!/bin/bash

######## README ##########
# Use this file to setup your local cluster environment in case your're not allowed to install the binaries system-wide
# (i.e. you dont have admin privileges on the machine)

#Set Paths FREESUFER
FREESURFER_HOME=/absolute/path/to/freesurfer

#Set Paths for FSL
FSLDIR=/absolute/path/to/fsl

#Set Paths for MRTrix
MRTrixDIR=/absolute/path/to/mrtrix/

#DCM2nii & More
PATH=${PATH}:/absolute/path/to/dcm2niiBinaries

# Add Python Packages
PYTHONPATH=${PYTHONPATH}:absolute/path/to/TVB-Pypeline-Folder

###################################################################
# Additional stuff required to be executed for using the toolboxes.
# Edit on own risk from here
###################################################################
SUBJECTS_DIR=${FREESURFER_HOME}/subjects
FUNCTIONALS_DIR=${FREESURFER_HOME}/sessions
PATH=${PATH}:${FREESURFER_HOME}/bin
export FREESURFER_HOME SUBJECTS_DIR FUNCTIONALS_DIR PATH
source ${FREESURFER_HOME}/FreeSurferEnv.sh
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh 

. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

LD_LIBRARY_PATH=${MRTrixDIR}/lib/
export LD_LIBRARY_PATH
PATH=${MRTrixDIR}/bin:${PATH}
export PATH

export PYTHONPATH