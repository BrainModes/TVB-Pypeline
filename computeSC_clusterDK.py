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

import sys
import os
import fnmatch
import csv
import numpy as np
# import scipy as sp

# Debug
path = '/Users/srothmei/Desktop/charite/toronto/QL_20120814/mrtrix_68/'
wmborder_file = path + 'masks_68/wmborder.npy'
roi = 2
subID = 'QL_20120814'
tracksPath = path + 'tracks_68/'
outfile = tracksPath + 'SC_row_' + str(roi) + subID + '.mat'

# User feedback
print('Computing SC for ROI ' + str(roi))

# Load the affine matrix
affineMatrix = np.load(path + 'masks_68/affine_matrix.npy')
# Load the wmborder file
wmborder = np.load(wmborder_file)

# Define the ROI-Range
region_table = range(1001, 1004) + range(1005, 1036) + range(2001, 2004) + range(2005, 2036)
# Generate ROI-ID to voxel hashtable
region_id_table = np.array((0, 0, 0, 0))  # Init Variable
for regid in region_table:
    tmpids = np.asarray(np.nonzero(wmborder == regid))
    tmpids = np.vstack((np.ones(np.shape(tmpids)[1]) * regid, tmpids))
    region_id_table = np.vstack((region_id_table, np.transpose(tmpids)))
region_id_table = region_id_table[1:, :]
SC_cap

# Count the number of failure tracks
qualityMetrics = {'off_seed': 0, 'too_short': 0, 'good_tracks': 0, 'wrong_seed': 0, 'expected_tracks': 0,
                  'wrong_target': 0, 'generated_tracks': 0}

# Examine the region
