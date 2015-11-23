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


def generate_masks(subPath, mask_output_folder, wmoutline2diff_1mm, wmparc2diff_1mm, seedsPerVoxel=None):

    # import sys
    # import os
    import numpy as np
    import nibabel as nib
    # import scipy as sp

    # Debug
    # subPath = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/'
    # pathOnCluster = '/home/petra/Simon/Pipeline/'
    # mask_output_folder = subPath + 'mrtrix_68/masks_68/'

    # Set the desired number of seedpoints per voxel
    if seedsPerVoxel is None:
        seedsPerVoxel = 200
        # seedsPerVoxel = 1000

    # Mask chunk size
    chunkSize = 100000 / seedsPerVoxel

    # Load the wmborder image
    # wmborder = nib.load(calc_images_folder + '/wmoutline2diff_1mm.nii.gz')
    wmborder = nib.load(wmoutline2diff_1mm)

    # Extract and Save the Affine Matrix for later use
    affine_matrix = np.linalg.inv(wmborder.affine)
    np.save(mask_output_folder + 'affine_matrix.npy', affine_matrix)

    # High-Res GM-WM-border
    # wmparc = nib.load(calc_images_folder + '/wmparc2diff_1mm.nii.gz')
    wmparc = nib.load(wmparc2diff_1mm)
    wmparc_data = wmparc.get_data()
    # Remove unwanted areas
    wmparc_data[wmparc_data > 3000] -= 2000
    wmparc_data[wmparc_data == 1004] = 0
    wmparc_data[wmparc_data == 2000] = 0
    wmparc_data[wmparc_data == 2004] = 0
    wmparc_data[wmparc_data < 1001] = 0
    wmparc_data[wmparc_data > 2035] = 0
    # Save the image
    nib.save(wmparc, mask_output_folder + 'wmparcMask_1mm.nii.gz')

    ### Create the parcellated GWI
    wmborder_data = wmborder.get_data()
    wmborder_data[wmborder_data > 0] = wmparc_data[wmborder_data > 0]
    nib.save(wmborder, mask_output_folder + 'gmwmborder_1mm.nii.gz')
    # Save it as "normal" File
    np.save(mask_output_folder + 'wmborder.npy', wmborder_data)

    ### Create Seed- and Target-Masks
    numseeds = np.zeros((1, 3))

    for i in range(1001, 1004) + range(1005, 1036) + range(2001, 2004) + range(2005, 2036):
        print('Processing RegionID ' + str(i))

        tmpimg = wmborder_data.copy()
        tmpimg[tmpimg != i] = 0
        tmpimg[tmpimg > 0] = 1
        # maskvoxel = np.ravel(tmpimg, 'F').nonzero()
        maskvoxel = np.nonzero(tmpimg)

        # Calculate the number of masks to generate for the current region
        nummasks = int(np.floor(float(np.shape(maskvoxel)[1]) / chunkSize))
        for j in range(nummasks):
            wmparc_data = np.zeros_like(tmpimg)
            indexRange = range(chunkSize * j, chunkSize * (j + 1))
            wmparc_data[maskvoxel[0][indexRange], maskvoxel[1][indexRange], maskvoxel[2][indexRange]] = 1
            # Save the generated Mask
            wmparc = nib.Nifti1Image(wmparc_data, wmparc.affine, wmparc.header)
            nib.save(wmparc, mask_output_folder + 'seedmask' + str(i) + str(j + 1) + '_1mm.nii.gz')

            # Store the number of voxels for the current masks and also the ID
            tmpfind = np.array([[int(str(i) + str(j + 1)), np.shape(indexRange)[0], i]])
            numseeds = np.concatenate((numseeds, tmpfind))

        # Now the remaining mask for the current region....
        wmparc_data = np.zeros_like(tmpimg)
        indexRange = range(chunkSize * nummasks, np.shape(maskvoxel)[1])
        wmparc_data[maskvoxel[0][indexRange], maskvoxel[1][indexRange], maskvoxel[2][indexRange]] = 1
        # Save the mask
        wmparc = nib.Nifti1Image(wmparc_data, wmparc.affine, wmparc.header)
        nib.save(wmparc, mask_output_folder + 'seedmask' + str(i) + str(nummasks + 1) + '_1mm.nii.gz')
        # Store the number of voxels for the current masks and also the ID
        tmpfind = np.array([[int(str(i) + str(nummasks + 1)), np.shape(indexRange)[0], i]])
        numseeds = np.concatenate((numseeds, tmpfind))

        # Create the corresponding target mask
        wmparc_data = wmborder_data.copy()
        wmparc_data[wmparc_data == i] = 0
        wmparc_data[wmparc_data > 0] = 1
        wmparc = nib.Nifti1Image(wmparc_data, wmparc.affine, wmparc.header)
        nib.save(wmparc, mask_output_folder + 'targetmask' + str(i) + '_1mm.nii.gz')

    # Finalize the seedcount array
    numseeds = numseeds[numseeds[:, 2] != 0]
    numseeds[:, 1] = numseeds[:, 1] * seedsPerVoxel
    # Store the array as an ASCII file
    np.savetxt(mask_output_folder + 'seedcount.txt', numseeds.astype(int), delimiter=' ', fmt='%1i')

    # Generate Batch File
    with open(mask_output_folder + 'batch_track.sh', 'w') as f:
        for roiid in range(np.shape(numseeds)[0]):
            f.write('{0} {1} {2} {3}\n'.format(subPath, str(int(numseeds[roiid, 0])), str(int(numseeds[roiid, 1])),
                                               str(int(numseeds[roiid, 2]))))

    f.close()
