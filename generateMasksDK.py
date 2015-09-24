#!/usr/bin/env python
#
# =============================================================================
# Authors: Michael Schirner, Simon Rothmeier, Petra Ritter
# BrainModes Research Group (head: P. Ritter)
# Charite University Medicine Berlin & Max Planck Institute Leipzig, Germany
# Correspondence: petra.ritter@charite.de
#
# When using this code please cite as follows:
# Schirner M, Rothmeier S, Jirsa V, McIntosh AR, Ritter P (in prep)
# Constructing subject-specific Virtual Brains from multimodal neuroimaging
#
# This software is distributed under the terms of the GNU General Public License
# as published by the Free Software Foundation. Further details on the GPL
# license can be found at http://www.gnu.org/copyleft/gpl.html.
# =============================================================================
import sys
import os
import numpy as np
#import scipy as sp
import nibabel as nib

#Get Input Parameters
#subPath = sys.argv[1]
#pathOnCluster = sys.argv[2]

#Debug
subPath = '/Users/srothmei/Desktop/charite/toronto/QL_20120814/'
pathOnCluster = '/home/petra/Simon/Pipeline/'

mask_output_folder = subPath + 'mrtrix_68/masks_68/'

#Set the desired number of seedpoints per voxel
seedsPerVoxel = 200
#seedsPerVoxel = 1000

#Mask chunk size
chunkSize = 100000/seedsPerVoxel

#Load the wmborder image
wmborder = nib.load(subPath + 'calc_images/wmoutline2diff_1mm.nii.gz')

#Extract and Save the Affine Matrix for later use
affine_matrix = np.linalg.inv(wmborder.affine)
np.save(mask_output_folder+'affine_matrix.npy',affine_matrix)

#High-Res GM-WM-border
wmparc = nib.load(subPath + 'calc_images/wmparc2diff_1mm.nii.gz')
wmparc_data = wmparc.get_data()
#Remove unwanted areas
wmparc_data[wmparc_data > 3000] = wmparc_data[wmparc_data > 3000] - 2000
wmparc_data[wmparc_data == 1004] = 0
wmparc_data[wmparc_data == 2000] = 0
wmparc_data[wmparc_data == 2004] = 0
wmparc_data[wmparc_data < 1001] = 0
wmparc_data[wmparc_data > 2035] = 0
#Save the image
nib.save(wmparc,mask_output_folder+'wmparcMask_1mm.nii.gz')

### Create the parcellated GWI
wmborder_data = wmborder.get_data()
wmborder_data[wmborder_data > 0] = wmparc_data[wmborder_data > 0]
nib.save(wmborder,mask_output_folder+'gmwmborder_1mm.nii.gz')
#Save it as "normal" File
np.save(mask_output_folder+'wmborder.npy',wmborder_data)

### Create Seed- and Target-Masks
counter = 0
for i in range(1001,1004)+range(1005,1036)+range(2001,2004)+range(2005,2036):
    print('Processing RegionID ' + str(i))

    tmpimg = wmborder_data.copy()
    tmpimg[tmpimg != i] = 0
    tmpimg[tmpimg > 0] = 1
    maskvoxel = np.ravel(tmpimg,'F').nonzero()
    #Calculate the number of masks to genrate for the current region
    nummasks = int(np.floor(np.shape(maskvoxel)[1] / chunkSize))
    for j in range(nummasks):
        wmparc_data = zeros_like()
