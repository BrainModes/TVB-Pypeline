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
import fnmatch
import csv
import numpy as np
# import scipy as sp
from scipy import io

# Debug
path = '/Users/srothmei/Desktop/charite/toronto/QL_20120814/bold'
subName = 'QL_20120814'

# Read the subID_ROIts.dat file
for datFile in os.listdir(path):
    if datFile.endswith('_ROIts.dat'):
        break
# Process the data-rows into a single matrix (timepoints X regions)
with open(path + '/' + datFile, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    for row in reader:
        tmp = np.asarray(map(float, filter(bool, row)))  # Remove empty entries and convert to float and array
        if reader.line_num == 1:  # Init the array
            fMRI = tmp
            continue
        fMRI = np.vstack((fMRI, tmp))

# Read the ROI-table
for statFile in os.listdir(path):
    if fnmatch.fnmatch(statFile, 'aparc*stats_cleared*'):
        break
with open(path + '/' + statFile, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    for row in reader:
        tmp = np.asarray(map(float, filter(bool, row)))  # Remove empty entries and convert to float and array
        if reader.line_num == 1:  # Init the array
            ROI_ID_table = tmp
            continue
        ROI_ID_table = np.vstack((ROI_ID_table, tmp))
# Clear the ROI-table and leave only the Desikan entries
start1 = np.shape(ROI_ID_table)[0] - 69
stop1 = start1 + 34
start2 = stop1 + 1
stop2 = np.shape(ROI_ID_table)[0]
fMRI_DK68 = fMRI[:, range(start1, stop1)+range(start2, stop2)]

# Compute the FC
FC_cc = np.corrcoef(np.transpose(fMRI_DK68))
# Correct for possible NaN values
FC_cc[np.isnan(FC_cc)] = 0

# Store the stuff (FC-matrix, timeseries, ROI_ID_table)
theDict = {subName + '_ROIts': fMRI, subName + '_ROIts_DK68':  fMRI_DK68, 'FC_cc':FC_cc, 'ROI_ID_table': ROI_ID_table}
io.savemat(path + '/' + subName + '_fMRI_new.mat', theDict)
