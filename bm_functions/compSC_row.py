#!/bin/python
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


def tck2voxel_cluster(tck, affineMatrix):
    from numpy import hstack, round, dot, ones, shape, transpose

    # Loop over all tracks in the data structure
    for i in range(shape(tck)[0]):
        # Transform the coordinates by multiplying the affine matrix with the coords
        # => First add a column of 1 to the right side of the matrix
        zw = hstack((tck[i], ones((shape(tck[i])[0], 1))))
        # => Second multiply the actual matrices
        zw = round(dot(zw, transpose(affineMatrix)))
        # => Third store the result excluding the righthand-side column
        tck[i] = zw[:, :3].astype('int16')

    return tck


def compute_connectivity_row(roi, subid, affine_matrix, wmborder, tracksPath, track_files):
    # import sys
    import os
    # import fnmatch
    import numpy as np
    import logging
    # import scipy as
    import re
    import json
    #from nipype.interfaces.mrtrix import convert as trk
    #import nibabel.trackvis as trk
    from collections import defaultdict
    from os.path import basename

    # Debug
    # path = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/mrtrix_68/'
    # roi = 2
    # subID = 'AJ_20140516_1600'

    # Calling Kenny Loggins
    dangerZone = logging.getLogger('interface')
    dangerZone.setLevel('WARNING')

    # Extract parameters
    #wmborder_file = path + 'masks_68/wmborder.npy'
    #tracksPath = path + 'tracks_68/'

    # User feedback
    print('Computing SC for ROI ' + str(roi))

    # Load the affine matrix
    #affine_matrix = np.load(path + 'masks_68/affine_matrix.npy')
    # Load the wmborder file
    #wmborder = np.load(wmborder_file)

    # Define the ROI-Range
    # region_table = range(1001, 1004) + range(1005, 1036) + range(2001, 2004) + range(2005, 2036)
    region_table = np.unique(wmborder[wmborder > 0]).astype(int)
    region_table = region_table.tolist()

    # Generate ROI-ID to voxel hashtable
    print('Generate ROI-ID to voxel hashtable...')
    region_id_table = np.array((0, 0))  # Init Variable
    for regid in region_table:
        tmpids = np.ravel_multi_index(np.nonzero(wmborder == regid), wmborder.shape, order='F')
        tmpids.sort()
        tmpids = np.vstack((np.ones_like(tmpids) * regid, tmpids))
        region_id_table = np.vstack((region_id_table, np.transpose(tmpids)))
    region_id_table = region_id_table[1:, :].astype(int)

    # Init storage
    SC_cap = defaultdict(list)
    SC_dist = defaultdict(list)

    # Count the number of failure tracks
    qualityMetrics = {'off_seed': 0, 'too_short': 0, 'good_tracks': 0, 'wrong_seed': 0,
                      'expected_tracks': np.count_nonzero(wmborder == region_table[roi - 1]) * 200, 'wrong_target': 0,
                      'generated_tracks': 0}

    ####### Examine the region ##################################################################

    # Check for track-files in the dircetory
    tileFiles = list()
    # for trackFile in os.listdir(tracksPath):

    for trackFile in track_files:
        # Select only files that have max. of 2 trailing number depicting the ordering...
        # i.e. don't select files like '10012_subID.tck' when processing for region 10...
        if re.search("^" + str(region_table[roi - 1]) + "\d{1,2}_.*\.npy$", basename(trackFile)):
            tileFiles.append(trackFile)

    # Loop over all track-files
    for tile in tileFiles:

        # First check if the file is sane considering the filesize...
        if os.stat(tile).st_size > 2000:
            # Load the file
            #tracks_header, tracks = trk.read_mrtrix_tracks(tile, as_generator=False)
            tracks_header = np.load(tile).item()
            print(tile + ': Tracks loaded .....')
            # Transform the coordinates from mm to voxel
            tracks_header['tracks'] = tck2voxel_cluster(tracks_header['tracks'], affine_matrix)

            # Step Size
            if 'step_size' in tracks_header.keys():
                step_size = tracks_header['step_size']
            elif 'step_length' in tracks_header.keys():
                step_size = tracks_header['step_length']
            else:
                step_size = 0.2

            # QA
            qualityMetrics['generated_tracks'] += np.shape(tracks_header['tracks'])[0]

            # Loop over the tracks themselves
            for trackind in range(np.shape(tracks_header['tracks'])[0]):
                # Find the "actual" seed-voxel: sometimes a track starts in a seed
                # voxel then heads into the wm and crosses another voxel of the
                # seeding-region. In this case we consider the last voxel on the
                # track path belonging to the seed region as the actual seed voxel.
                # Then we check whether the remaining path length is at least 10 mm
                # long.

                # Look which regions the track actually touches
                inRegIDs = wmborder[tracks_header['tracks'][trackind][:, 0], tracks_header['tracks'][trackind][:, 1], tracks_header['tracks'][trackind][:, 2]]
                # Create linear indices for all regions that are non-zero valued, EXCLUDING the end-point of the track
                inRegLin = np.flatnonzero(inRegIDs[:-1])
                # Create linear indices relative to the whole wmborder-img for all the steps in the track-path
                pathIndices = np.ravel_multi_index(
                    (tracks_header['tracks'][trackind][:, 0], tracks_header['tracks'][trackind][:, 1], tracks_header['tracks'][trackind][:, 2]), np.shape(wmborder),
                    order='F')

                if inRegLin.size > 0:  # Check if the path start on the border....
                    # Compute the track-length from the endpoint to the last point in the starting region
                    trackLen = np.shape(tracks_header['tracks'][trackind])[0] - inRegLin[-1] + 1
                    # Check if the track has a minimum length (8mm)
                    if trackLen > 8 / float(step_size):
                        # Check if the path has a valid endpoint
                        if inRegIDs[-1] > 0:
                            # Check if the region of the seedpoint matches the examined region
                            if region_table[roi - 1] == inRegIDs[inRegLin[-1]]:
                                # "[...] when you have eliminated the impossible, whatever remains, however improbable,
                                # must be the truth"
                                qualityMetrics['good_tracks'] += 1

                                # Find the entry in the storage-dict in which we'll store the current findings
                                seedID = np.flatnonzero(region_id_table[:, 1] == pathIndices[inRegLin[-1]])[0]
                                # Find the correspdonging entry on which the track refers to in the storage dict
                                targetID = np.flatnonzero(region_id_table[:, 1] == pathIndices[-1])[0]

                                # Write into the capacity storage dict...
                                SC_cap[seedID].append(targetID)  # Add a Connection from Seedvoxel to Targetvoxel
                                SC_cap[targetID].append(seedID)  # Add a Connection from Targetvoxel to Seedvoxel

                                # Distance computation
                                SC_dist[seedID].append(trackLen)
                                SC_dist[targetID].append(trackLen)

                            else:
                                qualityMetrics['wrong_seed'] += 1
                        else:
                            qualityMetrics['wrong_target'] += 1
                    else:
                        qualityMetrics['too_short'] += 1
                else:
                    qualityMetrics['off_seed'] += 1

    for i in SC_cap.keys():
        # Filter out the redundant connections i.e. just count distinct connections
        uniqueTMP, uniqueIndices = np.unique(SC_cap[i], return_index=True)
        SC_cap[i] = list(uniqueTMP)
        # Only take the length of the distinct connections into account
        SC_dist[i] = list(SC_dist[i][j] for j in uniqueIndices)

    # Store everything there is to store
    print('Storing....')
    # np.savez_compressed(outfile, SC_cap=SC_cap, SC_dist=SC_dist, qualityMetrics=qualityMetrics)
    SC_cap_row_filename = tracksPath + 'SC_cap_row_' + str(roi) + subid + '.json'
    with open(SC_cap_row_filename, 'w') as outfile:
        json.dump(SC_cap, outfile, sort_keys=True, indent=2)
        outfile.close()

    SC_dist_row_filename = tracksPath + 'SC_dist_row_' + str(roi) + subid + '.json'
    with open(SC_dist_row_filename, 'w') as outfile:
        json.dump(SC_dist, outfile, sort_keys=True, indent=2)
        outfile.close()

    with open(tracksPath + 'qualityMetrics_' + str(roi) + subid + '.json', 'w') as outfile:
        json.dump(qualityMetrics, outfile, sort_keys=True, indent=2)
        outfile.close()

    # Releasing Mr Loggins...
    dangerZone.setLevel('NOTSET')
    print('Done!')

    return SC_cap_row_filename, SC_dist_row_filename
