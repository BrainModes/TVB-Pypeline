
# coding: utf-8

# ## Tracking Thresholding using Morris et al (2008)
# MultiVoxel!

# In[2]:

import numpy as np
import os
from scipy.stats import poisson
from scipy.ndimage.measurements import label
from scipy.stats import mode
import nibabel as nib
from nipype.interfaces import mrtrix as mrt


# In[45]:

theDir = "/Users/srothmei/Desktop/charite/toronto/Adalberto/debug/"
csd_file = theDir + "csd8.nii.gz"
csd_file_random = theDir + "csd_random.nii.gz"

multiVoxelTracks = theDir + 'QL_10011_tracks.tck'
multiVoxelTracksRandom = theDir + 'QL_10011_tracks_random.tck'

singleVoxelTDI = theDir + 'QL_10011_tracks_multiVox_TDI.nii.gz'
singleVoxelRandomTDI = theDir + 'QL_10011_tracks_multiVox_random_TDI.nii.gz'

seedmask = theDir + 'seedmask10011_1mm.nii.gz'
tmpMaskFileName = theDir + 'seedmask1001_tmp.nii.gz'

singleVoxelTracks = theDir + 'QL_10011_tracks_singleVox.tck'
singleVoxelTracksRandom = theDir + 'QL_10011_tracks_singleVox_random.tck'

tracksPerVoxel = 200
trackPerVoxelRandom = tracksPerVoxel * 20
confidenceThreshold = 0.05 # P-Value


# ## Generate non-informative fODF data i.e. fODF data which has equi-distributed directionality information along all possible directions, leaving a sphere

# In[3]:

csd8 = nib.load(csd_file)
csd8_data = csd8.get_data()


# In[4]:

for i in range(csd8_data.shape[0]):
    for j in range(csd8_data.shape[1]):
        for k in range(csd8_data.shape[2]):
            if np.sum(csd8_data[i,j,k,:] == np.zeros_like(csd8_data[30,40,32,:])) < 45:
                csd8_data[i,j,k,:] = np.zeros_like(csd8_data[30,40,32,:])
                csd8_data[i,j,k,0] = 1.



csd8 = nib.Nifti1Image(csd8_data, csd8.affine, csd8.header)
nib.save(csd8, csd_file_random)


# ## Loop over the voxels in the seedmask


# Load the seedmask image
seedmask_image = nib.load(seedmask)
seedmask_data = seedmask_image.get_data().astype("int8")

# Get the indices inside the cube of the seedvoxels
seed_coords = np.transpose(np.nonzero(seedmask_data))

# Init the combined P-Map for later usage
raw_P_map_union = np.zeros_like(seedmask_data, dtype='float64')


# Loop over those indices
counter = 1
for seed_vox in seed_coords:

    print counter
    counter += 1
    
    # First generate a new temporary seedmask for the current voxel
    seedmask_data = np.zeros_like(seedmask_data)
    seedmask_data[seed_coords[0, 0], seed_coords[0, 1], seed_coords[0, 2]] = 1
    # Save the mask
    tmpMask = nib.Nifti1Image(seedmask_data, seedmask_image.affine, seedmask_image.header)
    nib.save(tmpMask, tmpMaskFileName)

    # ## Now perform fiber tracking for a single-voxel ROI using the random-CSD-map.
    # While using the random map, the number of tracks-per-voxel is 20x higher to achive statistical power
    print "Tracking...."

    tracker = mrt.StreamlineTrack()
    tracker.inputs.inputmodel = 'SD_PROB'
    tracker.inputs.stop = True
    tracker.inputs.minimum_tract_length = 30
    tracker.inputs.no_mask_interpolation = True
    tracker.inputs.step_size = 0.2
    tracker.inputs.unidirectional = True #Maybe off?
    tracker.inputs.seed_file = tmpMaskFileName
    tracker.inputs.include_file = theDir + 'targetmask1001_1mm.nii.gz'
    tracker.inputs.mask_file = theDir + 'wmmask_68_1mm.nii.gz'

    tracker.inputs.in_file = csd_file_random
    tracker.inputs.out_file = singleVoxelTracksRandom
    tracker.inputs.desired_number_of_tracks = trackPerVoxelRandom
    tracker.run()

    print "Extracting...."
    # Extract the tracks from the single Voxel (informed) case
    tracksFilter = mrt.FilterTracks()
    tracksFilter.inputs.in_file = multiVoxelTracks
    tracksFilter.inputs.no_mask_interpolation = True
    tracksFilter.inputs.include_file = tmpMaskFileName
    tracksFilter.inputs.invert = False
    tracksFilter.inputs.out_file = singleVoxelTracks
    tracksFilter.run()
    #
    #tracksFilter.inputs.in_file = multiVoxelTracksRandom
    #tracksFilter.inputs.out_file = singleVoxelTracksRandom
    #tracksFilter.run()


    # ## Generate maps of the connection probability

    tdi = mrt.Tracks2Prob()
    tdi.inputs.fraction = False
    tdi.inputs.template_file = seedmask

    tdi.inputs.in_file = singleVoxelTracks
    tdi.inputs.out_filename = singleVoxelTDI
    tdi.run()

    tdi.inputs.in_file = singleVoxelTracksRandom
    tdi.inputs.out_filename = singleVoxelRandomTDI
    tdi.run()


    # ## Now calculate the Z-Map

    tdi_informed = nib.load(singleVoxelTDI)
    tdi_informed_data = tdi_informed.get_data()

    tdi_random = nib.load(singleVoxelRandomTDI)
    tdi_random_data = tdi_random.get_data()

    # Remove the TDI files since they algorithm throws an error when trying to overwrite them with new iterations
    os.remove(singleVoxelTDI)
    os.remove(singleVoxelRandomTDI)

    # Equation (4)
    meanV = (tracksPerVoxel * tdi_random_data) / float(trackPerVoxelRandom)

    # Equation (5)
    stdV = np.sqrt(meanV)


    # ## Compute Stat Measures

    #Quick and dirty loop
    Z_Map = np.zeros_like(meanV, dtype='float64')
    #raw_P_Map = np.zeros_like(Z_Map)
    for x in range(np.shape(meanV)[0]):
        for y in range(np.shape(meanV)[1]):
            for z in range(np.shape(meanV)[2]):
                if stdV[x,y,z] > 0.0:
                    # Equation (7)
                    tmp = (tdi_informed_data[x,y,z] - meanV[x,y,z]) / stdV[x,y,z]
                    #if tmp > 0.0:
                    #    Z_Map[x,y,z] = np.log(tmp)
                    #else:
                    Z_Map[x,y,z] = tmp
                    #Omit the negative Z-values
                    if tmp >= 0.0:
                        k = tdi_informed_data[x,y,z]
                        mu = meanV[x,y,z]
                        # Only update the value if its larger
                        raw_P_map_union[x,y,z] = max(raw_P_map_union[x,y,z], 1 - poisson.pmf(k, mu))

# ### End the Voxel-Loop here! ###

# Thresholded P-map
P_Map_thresholded = np.zeros_like(raw_P_map_union, dtype="int16")
P_Map_thresholded[raw_P_map_union >= 1 - (confidenceThreshold / np.shape(seed_coords)[0])] = 1

                    
# Some debug/visual stuff
#zmapImage = nib.Nifti1Image(Z_Map, tdi_informed.affine, tdi_informed.header)
#nib.save(zmapImage, theDir + 'Zmap.nii.gz')

#pmapImage = nib.Nifti1Image(raw_P_Map, tdi_informed.affine, tdi_informed.header)
#nib.save(pmapImage, theDir + 'Pmap_raw.nii.gz')

#pmapThres = nib.Nifti1Image(P_Map_thresholded, tdi_informed.affine, tdi_informed.header)
#nib.save(pmapThres, theDir + 'Pmap_thresholded.nii.gz')


# ## False positive correction
# Check if the thresholded voxels have a connection to the seed voxel. If not remove them

# In[ ]:

structure = np.ones((3,3,3))
tmp = np.zeros_like(P_Map_thresholded, dtype="int16")
tmp, bar = label(P_Map_thresholded, structure)
modal_val, modal_count = mode(tmp[tmp>0], axis=None)
P_Map_thresholded_fpCorr = np.zeros_like(P_Map_thresholded)
P_Map_thresholded_fpCorr[tmp == modal_val] = 1

# Debugging / Testing
pmapThresFP = nib.Nifti1Image(P_Map_thresholded_fpCorr, tdi_informed.affine, tdi_informed.header)
nib.save(pmapThresFP, theDir + 'Pmap_thresholded_FPcorr.nii.gz')


# ## Apply the P-Map onto the tracks

# In[ ]:

# First invert the mask to apply it with MRTrix's tracks_filter
P_Map_thresholded_fpCorr_inv = np.invert(P_Map_thresholded_fpCorr.astype(bool)).astype(int)
# Now save it to use it!
pmapThresFP_inv = nib.Nifti1Image(P_Map_thresholded_fpCorr_inv, tdi_informed.affine, tdi_informed.header)
nib.save(pmapThresFP_inv, theDir + 'Pmap_thresholded_FPcorr_inv.nii.gz')

tracksFilter = mrt.FilterTracks()
tracksFilter.inputs.in_file = multiVoxelTracks
tracksFilter.inputs.no_mask_interpolation = True
tracksFilter.inputs.exclude_file = theDir + 'Pmap_thresholded_FPcorr_inv.nii.gz'
tracksFilter.inputs.invert = False
tracksFilter.inputs.out_file = theDir + 'QL_10011_tracks_filt.tck'

tracksFilter.run()

