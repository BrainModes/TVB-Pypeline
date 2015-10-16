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
import os
from scipy import io
import numpy as np
import json
import logging
from collections import defaultdict
from nibabel.freesurfer import io as fs


# Debug
subID = 'QL_20120814'
subFolder = '/Users/srothmei/Desktop/charite/toronto/QL_20120814/'
SC_matrix = 'QL_20120814_SC.mat'
reconallFolder = 'recon_all'

# Create the results folder
os.mkdir(subFolder + 'results/')

# Load the SC matrix
SC = io.loadmat(subFolder + '/mrtrix_68/tracks_68/' + SC_matrix)
weights = SC['SC_cap_agg_bwflav2']

# Load the required things compiuted previously by FREESURFER
lh_vert, lh_faces = fs.read_geometry(subFolder + '/' + reconallFolder + '/surf/lh.pial')
rh_vert, rh_faces = fs.read_geometry(subFolder + '/' + reconallFolder + '/surf/rh.pial')
