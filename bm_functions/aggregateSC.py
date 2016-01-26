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


def aggregate_connectivity(sub_id, wmborder, tracksPath, cap_row_files, dist_row_files, steplength=None):
    # import sys
    # import os
    from scipy import io
    import numpy as np
    import json
    import logging
    from collections import defaultdict

    if steplength is None:
        steplength = 0.2

    # Debug
    # path = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/mrtrix_68/'
    #  sub_id = 'AJ_20140516_1600'

    # Build some parameters...
    # wmborder_file = path + 'masks_68/wmborder.npy'
    # tracksPath = path + 'tracks_68/'

    # Init the logger
    logger = logging.getLogger('interface')
    logger.setLevel('INFO')

    #load the wmborder
    wmborder = np.load(wmborder)

    # Define the ROI-Range
    # region_table = range(1001, 1004) + range(1005, 1036) + range(2001, 2004) + range(2005, 2036)
    region_table = np.unique(wmborder[wmborder > 0])

    counter = 0
    # Generate ROI-ID to voxel hashtable
    logger.info('Generate ROI-ID to voxel hashtable...')
    inverse_region_table = np.zeros((1, np.max(region_table) + 1))
    region_id_table = np.array((0, 0))  # Init Variable
    for regid in region_table:
        counter += 1
        # Transfer table between DK-Numbering and Matrix Numbering...
        inverse_region_table[:, regid] = counter

        # tmpids = np.transpose(np.asarray(np.nonzero(wmborder == regid)))
        tmpids = np.ravel_multi_index(np.nonzero(wmborder == regid), wmborder.shape, order='F')
        tmpids.sort()
        tmpids = np.vstack((np.ones_like(tmpids) * regid, tmpids))
        region_id_table = np.vstack((region_id_table, np.transpose(tmpids)))
    region_id_table = region_id_table[1:, :].astype(int)

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
    logger.info('Processing SC_files... ')
    #for roi in range(len(region_table)):
    cap_row_files.sort()
    dist_row_files.sort()
    for cap_file, dist_file in zip(cap_row_files, dist_row_files):
        #logger.info('Processing ROI: ' + str(roi + 1))

        #logger.info('Load the SC_row_files...')
        # with open(tracksPath + 'SC_cap_row_' + str(roi + 1) + sub_id + '.json', 'r') as inFile:
        with open(cap_file, 'r') as inFile:
            SC_cap_row = json.load(inFile)
            inFile.close()

        # with open(tracksPath + 'SC_dist_row_' + str(roi + 1) + sub_id + '.json', 'r') as inFile:
        with open(dist_file, 'r') as inFile:
            SC_dist_row = json.load(inFile)
            inFile.close()

        # Loop over all border voxels of the current roi
        for seedVoxelID in SC_cap_row.keys():
            # Insert findings from the current roi into the storage matrix
            SC_cap_agg_tmp[seedVoxelID].extend(SC_cap_row[seedVoxelID])
            SC_dist_agg_tmp[seedVoxelID].extend(SC_dist_row[seedVoxelID])

    # Now build the capacity matrices by looking at the previously build storage matrix
    # SC_cap_counts_tmp = SC_cap_agg_tmp.copy()
    logger.info('Now build the capacity matrices...')
    # for seedVoxelID in SC_cap_agg_tmp.keys():
    theKeys = SC_cap_agg_tmp.keys()
    for key in range(len(theKeys)):
        seedVoxelID = theKeys[key]
        logger.debug(str(key) + ' / ' + str(len(theKeys)))
        # First check for uniqueness of the tracks...
        uniqueTMP, uniqueIndices, uniqueInverse = np.unique(SC_cap_agg_tmp[seedVoxelID], return_index=True,
                                                            return_inverse=True)
        # Find out how often each unique values is occuring
        binCount = np.bincount(uniqueInverse)

        SC_cap_agg_tmp[seedVoxelID] = list(uniqueTMP)
        # Only take the length of unique tracks into account + apply the steplength!
        SC_dist_agg_tmp[seedVoxelID] = list(SC_dist_agg_tmp[seedVoxelID][j] * steplength for j in uniqueIndices)

        # Now build the acutal capacity metrics...
        # First find out in which region he current seed-voxel lies
        seedRegionID = np.flatnonzero(region_table == region_id_table[seedVoxelID, 0])[0]
        # Find out which regions are reached by the current seedVoxel
        targetRegionIDs = inverse_region_table[0, region_id_table[SC_cap_agg_tmp[seedVoxelID], 0] - 1]

        # Insert a connection into the matrices at each target-region
        # Note that this is by definition symetric without explicitly defining it here since the connections
        # have been created symetrically by the computeSCRows-script!
        logger.debug('Inserting stuff for the following No. of target-regions: ' + str(len(targetRegionIDs)))
        taRegLen = len(targetRegionIDs)
        for i in range(taRegLen):
            target = int(targetRegionIDs[i] - 1)
            SC_cap_agg_counts[seedRegionID, target] += binCount[i]
            SC_cap_agg_bwflav1[seedRegionID, target] += 1
            SC_cap_agg_bwflav2[seedRegionID, target] += 1. / taRegLen
            SC_dist_agg[seedRegionID][target].extend([SC_dist_agg_tmp[seedVoxelID][i]])

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
    logger.info('Storing....')
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
    SC_matrix_filename_matlab = tracksPath + sub_id + '_SC.mat'
    io.savemat(SC_matrix_filename_matlab, theMatlabDict)

    # Save the more complex output plus the matlab-one into json
    # First convert all items to lists...
    for i in theMatlabDict.keys():
        theMatlabDict[i] = theMatlabDict[i].tolist()

    SC_matrix_filename_json = tracksPath + sub_id + '_SC.json'
    with open(SC_matrix_filename_json, 'w') as outfile:
        json.dump(theMatlabDict, outfile, sort_keys=True, indent=2)
        outfile.close()

    # Now store the big distance matrix into json format
    SC_matrix_voxelwise_filename = tracksPath + sub_id + '_SC_voxelwiseDistanceMatrix_counts.json'
    with open(SC_matrix_voxelwise_filename, 'w') as outfile:
        json.dump(SC_dist_agg, outfile, sort_keys=True, indent=2)
        outfile.close()

    SC_matrix_voxelwise_distinctConnections_distance_filename = tracksPath + sub_id + '_SC_voxelwiseDistanceMatrix_distinctConnections.json'
    with open(SC_matrix_voxelwise_distinctConnections_distance_filename, 'w') as outfile:
        json.dump(SC_dist_agg_tmp, outfile, sort_keys=True, indent=2)
        outfile.close()

    SC_matrix_voxelwise_distinctConnections_capacity_filename = tracksPath + sub_id + '_SC_voxelwiseCapacityMatrix_distinctConnections.json'
    with open(SC_matrix_voxelwise_distinctConnections_capacity_filename, 'w') as outfile:
        json.dump(SC_cap_agg_tmp, outfile, sort_keys=True, indent=2)
        outfile.close()

    logger.info('Done!')

    return SC_matrix_filename_matlab, SC_matrix_filename_json, SC_matrix_voxelwise_filename, \
           SC_matrix_voxelwise_distinctConnections_capacity_filename,\
           SC_matrix_voxelwise_distinctConnections_distance_filename
