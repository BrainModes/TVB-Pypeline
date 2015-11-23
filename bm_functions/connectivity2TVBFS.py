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
import zipfile
from nibabel.freesurfer import io as fs


# Define function to remove free boundary vertices from the mesh
def removeFB(vertices, faces, normals, labels):
    # Look for all vertex-indices that occur only once in the faces-def.-matrix
    uniqueTMP, uniqueIndices, uniqueInverse = np.unique(faces, return_index=True, return_inverse=True)
    binCount = np.bincount(uniqueInverse)
    # Remove the found vertices
    tmp = binCount == 1
    vertices = vertices[np.invert(tmp), :]
    normals = normals[np.invert(tmp), :]

    # Now remove the  corresponding faces and labels
    FBtri = np.flatnonzero(tmp)
    for i in range(FBtri.shape[0]):
        faces == np.delete(faces, FBtri[i], axis=0)
        # Move the indexing of the other faces
        faces[faces > FBtri[i]] -= 1

        # Remove the vertices from the labeling table
        labels = np.delete(labels, FBtri[i], axis=0)
        labels[labels > FBtri[i]] -= 1

        FBtri -= 1

    return vertices, faces, normals, labels


# Define function to calculate vertex normals
def calcVertNormals(vertices, faces):
    vertNormals = np.zeros(vertices.shape, dtype=vertices.dtype)
    triangles = vertices[faces]
    # Calculate face-normals using cross product for each triangle
    faceNorm = np.cross(triangles[::, 1] - triangles[::, 0], triangles[::, 2] - triangles[::, 0])
    # Normalize the vectors to length 1
    vecLen = np.sqrt(faceNorm[:, 0] ** 2 + faceNorm[:, 1] ** 2 + faceNorm[:, 2] ** 2)
    faceNorm /= np.transpose(np.vstack((vecLen, vecLen, vecLen)))

    # Insert the values
    vertNormals[faces[:, 0]] += faceNorm
    vertNormals[faces[:, 1]] += faceNorm
    vertNormals[faces[:, 2]] += faceNorm
    # Final normalization
    vecLen = np.sqrt(vertNormals[:, 0] ** 2 + vertNormals[:, 1] ** 2 + vertNormals[:, 2] ** 2)
    vertNormals /= np.transpose(np.vstack((vecLen, vecLen, vecLen)))

    return vertNormals


# Debug
# subID = 'AJ_20140516_1600'
# subFolder = '/Users/srothmei/Desktop/charite/toronto/AJ_20140516_1600/'
# SC_matrix = 'AJ_20140516_1600_SC.mat'
# reconallFolder = 'recon_all'

def connectivity_2_tvb_fs(subID, subFolder, SC_matrix, reconallFolder='recon_all'):

    # Create the results folder
    if not os.path.exists(subFolder + 'results/'):
        os.mkdir(subFolder + 'results/')

    # Load the SC matrix
    SC = io.loadmat(subFolder + '/mrtrix_68/tracks_68/' + SC_matrix)
    weights = SC['SC_cap_agg_bwflav2']
    delay = SC['SC_dist_agg_mean']

    # Load the required things compiuted previously by FREESURFER
    lh_vert, lh_faces = fs.read_geometry(subFolder + '/' + reconallFolder + '/surf/lh.pial')
    rh_vert, rh_faces = fs.read_geometry(subFolder + '/' + reconallFolder + '/surf/rh.pial')
    cortexMesh = {'vertices': np.vstack((lh_vert, rh_vert)),
                  'faces': np.vstack((lh_faces, rh_faces + np.shape(lh_vert)[0]))}

    # Calculate vertex-normals
    cortexMesh['vertexNormals'] = calcVertNormals(cortexMesh['vertices'], cortexMesh['faces'])

    # Load annotation tables
    lh_labels, lh_ctab, lh_names = fs.read_annot(subFolder + '/' + reconallFolder + '/label/lh.aparc.annot')
    rh_labels, rh_ctab, rh_names = fs.read_annot(subFolder + '/' + reconallFolder + '/label/rh.aparc.annot')
    # Remove the CC i.e. correct labeling
    lh_labels[lh_labels > 3] -= 1
    rh_labels[rh_labels > 3] -= 1
    # Combine into single vectors
    rh_labels += np.max(lh_labels)
    rh_labels[rh_labels == np.min(rh_labels)] = -1
    cortexMesh['labels'] = np.hstack((lh_labels, rh_labels))
    # Store label name-strings
    tmp = lh_names[1:4] + lh_names[5:]
    tmp_lh = ['lh_' + s for s in tmp]
    tmp_rh = ['rh_' + s for s in tmp]
    cortexMesh['labelNames'] = tmp_lh + tmp_rh

    # Do the TVB mesh clean
    cortexMesh['vertices'], cortexMesh['faces'], cortexMesh['vertexNormals'], cortexMesh['labels'] = removeFB(
        cortexMesh['vertices'], cortexMesh['faces'], cortexMesh['vertexNormals'], cortexMesh['labels'])

    # Now finally start storing things....
    # ############

    # Define the filenames
    filenames = ['weights.txt', 'centres.txt', 'tract.txt', 'orientation.txt', 'area.txt', 'cortical.txt', 'hemisphere.txt']

    # 1.) Weights
    np.savetxt(subFolder + 'results/' + filenames[0], weights, delimiter=' ', fmt='%1i')

    # 2.) Position
    # Calc region centers
    # centers = np.zeros((weights.shape[0], 3))
    with open(subFolder + 'results/' + filenames[1], 'w') as f:
        for i in range(weights.shape[0]):
            # First get all vertices corresponding to a certain region
            regionVertices = cortexMesh['vertices'][cortexMesh['labels'] == i + 1]
            # Compute the mean of each region
            tmp = np.mean(regionVertices, axis=0)
            # Now look for the nearest neighbors
            idx = np.sum(np.abs(cortexMesh['vertices'] - tmp), axis=1).argmin()
            # Define the nearest vertex as center
            # centers[i, :] = cortexMesh['vertices'][idx, :]
            center = cortexMesh['vertices'][idx, :]
            # Write file
            f.write('{0} {1} {2} {3}\n'.format(cortexMesh['labelNames'][i], str(center[0]), str(center[1]), str(center[2])))
        f.close()

    # 3.) Tract
    np.savetxt(subFolder + 'results/' + filenames[2], delay, delimiter=' ', fmt='%1i')

    # 4.) Orientation
    with open(subFolder + 'results/' + filenames[3], 'w') as f:
        for i in range(weights.shape[0]):
            # Get all vertex-normals correspodning to the vertices of the current region
            regionVertexNormals = cortexMesh['vertexNormals'][cortexMesh['labels'] == i + 1]
            # Compute mean vector
            orientation = np.mean(regionVertexNormals, axis=0)
            # Normalize it
            orientation /= np.sqrt(orientation[0]**2 + orientation[1]**2 + orientation[2]**2)
            # Write to file
            f.write('{0} {1} {2}\n'.format(str(orientation[0]), str(orientation[1]), str(orientation[2])))
        f.close()

    # 5.) Area
    # I'm not quite sure how to get to the exact value for the surface in mm^2
    # so for now i just count the surface vertices corresponding to each region
    # EDIT: According to the TVB Dokumentation, this attribute is not mandatory
    # for the Input!
    with open(subFolder + 'results/' + filenames[4], 'w') as f:
        for i in range(weights.shape[0]):
            area = np.count_nonzero(cortexMesh['labels'] == i)
            f.write('{0}\n'.format(str(area)))
        f.close()

    # 6.) Cortical
    # Since in the default atlas all areas are cortical
    cortical = np.ones((68, 1))
    np.savetxt(subFolder + 'results/' + filenames[5], cortical, delimiter=' ', fmt='%1i')

    # 7.) Hemisphere
    # Again hard coded for Desikan-Killany Mask!
    # TODO: Make this flexible!
    hemisphere = np.vstack((np.zeros((34, 1)), np.ones((34, 1))))
    np.savetxt(subFolder + 'results/' + filenames[6], hemisphere, delimiter=' ', fmt='%1i')

    # Assemble the Zip-File
    zf = zipfile.ZipFile(subFolder + 'results/' + subID + '_Connectivity.zip', mode='w')
    for fname in filenames:
        zf.write(subFolder + 'results/' + fname)
        os.remove(subFolder + 'results/' + fname)
    zf.close()
