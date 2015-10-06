#!/usr/bin/env python
#This function converts the Scanner coordinates from MRTRix's tck-Files to
#Voxel Coordinates using the affine Transformation matrix inside the header
#of a reference image (e.g. the Brainmask used for tracking!)
#INPUT:
#tck - The Struct obtained via the MRTrix MATLAB function read_mrtrix_tracks
#refimage - The NIFTI File with the affine Matrix in it's header
#OUTPUT:
#The tck-Struct with transformed coordinates
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

import numpy as np

def tck2voxel_cluster(tck, affineMatrix):

    # Loop over all tracks in the data structure
    for i in range(np.shape(tck)[0]):
        # Transform the coordinates by multiplying the affine matrix with the coords
        # => First add a column of 1 to the right side of the matrix
        zw = np.hstack((tck[i], np.ones((np.shape(tck[i])[0], 1))))
        # => Second multiply the actual matrices
        zw = np.round(np.dot(zw, np.transpose(affineMatrix)))
        # => Third store the result excluding the righthand-side column
        tck[i] = zw[:, :3]

    return tck
