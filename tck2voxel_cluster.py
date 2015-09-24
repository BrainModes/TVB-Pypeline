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
# Charité University Medicine Berlin & Max Planck Institute Leipzig, Germany
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
def tck2voxel_cluster(tck, affine_matrix):
    for ii in tck:
        print(ii['data'])
    return tck
