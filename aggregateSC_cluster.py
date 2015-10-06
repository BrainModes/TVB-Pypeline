#!/bin/python
# -*- coding: utf-8 -*-
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

# import sys
# import os
from scipy import io
import numpy as np
import json
from collections import defaultdict


def agg_sc(path, sub_id, steplength=0.2):
    # Debug
    # path = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/mrtrix_68/'
    #  sub_id = 'AJ_20140516_1600'

    # Build some parameters...
    wmborder_file = path + 'masks_68/wmborder.npy'
    tracksPath = path + 'tracks_68/'

    # Load the wmborder file
    wmborder = np.load(wmborder_file)

    # Define the ROI-Range
    region_table = range(1001, 1004) + range(1005, 1036) + range(2001, 2004) + range(2005, 2036)

    counter = 0
    # Generate ROI-ID to voxel hashtable
    print('Generate ROI-ID to voxel hashtable...')
    inverse_region_table = np.zeros((1, np.max(region_table)))
    region_id_table = np.array((0, 0))  # Init Variable
    for regid in region_table:
        counter += 1
        # Transfer table between DK-Numbering and Matrix Numbering...
        inverse_region_table[:, regid - 1] = counter

        tmpids = np.transpose(np.asarray(np.nonzero(wmborder == regid)))
        tmpids = np.ravel_multi_index((tmpids[:, 0], tmpids[:, 1], tmpids[:, 2]), np.shape(wmborder), order='F')
        tmpids = np.vstack((np.ones_like(tmpids) * regid, tmpids))
        region_id_table = np.vstack((region_id_table, np.transpose(tmpids)))
    region_id_table = region_id_table[1:, :]

    # Init storage
    SC_cap_agg_tmp = defaultdict(list)
    SC_dist_agg_tmp = defaultdict(list)

    # Init output
    SC_cap_agg_bwflav1 = np.zeros((len(region_table), len(region_table)))
    SC_cap_agg_bwflav2 = np.zeros_like(SC_cap_agg_bwflav1)
    SC_cap_agg_counts = np.zeros_like(SC_cap_agg_bwflav1)
    # SC_dist_agg = defaultdict(list)
    SC_dist_agg = defaultdict(lambda: defaultdict(list))  # Create a matrix containing lists!
    SC_dist_agg_mean = np.zeros_like(SC_cap_agg_bwflav1)
    SC_dist_agg_median = np.zeros_like(SC_cap_agg_bwflav1)

    # Now loop over the regions....
    for roi in range(len(region_table)):
        print('Processing ROI: ' + str(roi + 1))

        print('Load the SC_row_files...')
        with open(tracksPath + 'SC_cap_row_' + str(roi + 1) + sub_id + '.json', 'r') as inFile:
            SC_cap_row = json.load(inFile)
            inFile.close()

        with open(tracksPath + 'SC_dist_row_' + str(roi + 1) + sub_id + '.json', 'r') as inFile:
            SC_dist_row = json.load(inFile)
            inFile.close()

        # Loop over all border voxels of the current roi
        for seedVoxelID in SC_cap_row.keys():
            # Insert findings from the current roi into the storage matrix
            SC_cap_agg_tmp[seedVoxelID].extend(SC_cap_row[seedVoxelID])
            SC_dist_agg_tmp[seedVoxelID].extend(SC_dist_row[seedVoxelID])

    # Now build the capacity matrices by looking at the previously build storage matrix
    # SC_cap_counts_tmp = SC_cap_agg_tmp.copy()
    for seedVoxelID in SC_cap_agg_tmp.keys():
        # First check for uniqueness of the tracks...
        uniqueTMP, uniqueIndices, uniqueInverse = np.unique(SC_cap_agg_tmp[seedVoxelID], return_index=True,
                                                            return_inverse=True)
        SC_cap_agg_tmp[seedVoxelID] = list(uniqueTMP)
        # Only take the length of unique tracks into account + apply the steplength!
        SC_dist_agg_tmp[seedVoxelID] = list(SC_dist_agg_tmp[seedVoxelID][j] * steplength for j in uniqueIndices)

        # Now build the acutal capacity metrics...
        # First find out in which region he current seed-voxel lies
        seedRegionID = np.flatnonzero(region_table == region_id_table[seedVoxelID, 0])[0]
        # Find out which regions are reached by the current seedVoxel
        targetRegionIDs = inverse_region_table[0, region_id_table[SC_cap_agg_tmp[seedVoxelID], 0]]

        # Insert a connection into the matrices at each target-region
        # Note that this is by definition symetric without explicitly defining it here since the connections
        # have been created symetrically by the computeSCRows-script!
        for i in range(len(targetRegionIDs)):
            SC_cap_agg_counts[seedRegionID, targetRegionIDs[i]] += sum(uniqueInverse == i)
            SC_cap_agg_bwflav1[seedRegionID, targetRegionIDs[i]] += 1
            SC_cap_agg_bwflav2[seedRegionID, targetRegionIDs[i]] += 1. / len(targetRegionIDs)
            SC_dist_agg[seedRegionID][targetRegionIDs[i]].extend(SC_dist_agg_tmp[seedVoxelID][targetRegionIDs[i]])

    # Finalize the distance metrics and also the counts
    for i in SC_dist_agg.keys():
        for j in SC_dist_agg[i].keys():
            SC_dist_agg_mean[i, j] = np.mean(SC_dist_agg[i][j])
            SC_dist_agg_median[i, j] = np.median(SC_dist_agg[i][j])

    # Normalize the capacity matrices
    numberOfTracks = np.sum(SC_cap_agg_counts)
    avgSeedingVoxels = 1  # TODO: Average the seeding voxels over 50 subjects!
    SC_cap_agg_counts_norm = SC_cap_agg_counts / float(numberOfTracks) * avgSeedingVoxels
    SC_cap_agg_bwflav1_norm = SC_cap_agg_bwflav1 / float(numberOfTracks) * avgSeedingVoxels
    SC_cap_agg_bwflav2_norm = SC_cap_agg_bwflav2 / float(numberOfTracks) * avgSeedingVoxels

    # Log
    SC_cap_agg_counts_norm_log = np.log(SC_cap_agg_counts_norm + 1)
    SC_cap_agg_bwflav1_norm_log = np.log(SC_cap_agg_bwflav1_norm + 1)
    SC_cap_agg_bwflav2_norm_log = np.log(SC_cap_agg_bwflav2_norm + 1)

    # Save the output into json and MATLAB format
    print('Storing....')
    # For matlab only the "simple" matrices, not the structs....
    # TODO: Look for a way to store this also in MATLAB format!
    theMatlabDict = {'SC_cap_agg_counts': SC_cap_agg_counts,
                     'SC_cap_agg_bwflav1': SC_cap_agg_bwflav1,
                     'SC_cap_agg_bwflav2': SC_cap_agg_bwflav2,
                     'SC_cap_agg_counts_norm': SC_cap_agg_counts_norm,
                     'SC_cap_agg_bwflav1_norm': SC_cap_agg_bwflav1_norm,
                     'SC_cap_agg_bwflav2_norm': SC_cap_agg_bwflav2_norm,
                     'SC_cap_agg_counts_norm_log': SC_cap_agg_counts_norm_log,
                     'SC_cap_agg_bwflav1_norm_log': SC_cap_agg_bwflav1_norm_log,
                     'SC_cap_agg_bwflav2_norm_log': SC_cap_agg_bwflav2_norm_log,
                     'SC_dist_agg_mean': SC_dist_agg_mean,
                     'SC_dist_agg_median': SC_dist_agg_median}
    io.savemat(tracksPath + sub_id + '_SC.mat', theMatlabDict)

    # Save the more complex output plus the matlab-one into json
    # First convert all items to lists...
    for i in theMatlabDict.keys():
        theMatlabDict[i] = theMatlabDict[i].tolist()

    with open(tracksPath + sub_id + '_SC.json', 'w') as outfile:
        json.dump(theMatlabDict, outfile, sort_keys=True, indent=2)
        outfile.close()

    # Now store the big distance matrix into json format
    with open(tracksPath + sub_id + '_SC_voxelwiseDistanceMatrix_counts.json', 'w') as outfile:
        json.dump(SC_dist_agg, outfile, sort_keys=True, indent=2)
        outfile.close()

    with open(tracksPath + sub_id + '_SC_voxelwiseDistanceMatrix_distinctConnections.json', 'w') as outfile:
        json.dump(SC_dist_agg_tmp, outfile, sort_keys=True, indent=2)
        outfile.close()

    with open(tracksPath + sub_id + '_SC_voxelwiseCapacityMatrix_distinctConnections.json', 'w') as outfile:
        json.dump(SC_cap_agg_tmp, outfile, sort_keys=True, indent=2)
        outfile.close()

    print('Done!')
