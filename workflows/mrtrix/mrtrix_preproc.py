
# coding: utf-8

# # Pipeline MRTrix 0.2 Preprocessing Workflow

# This workflow performs the specific preprocessing of the dwMRI data in order to perform tracking using <b>MRTrix 0.2</b> as direct legacy conversion from the old pipeline

# ### Handle imports

# In[1]:

from nipype.interfaces import mrtrix
from nipype import Node, Workflow, Function
from nipype.interfaces.utility import IdentityInterface

import logging


# ### Start the logging

logger = logging.getLogger('interface')
logger.setLevel(logging.INFO)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)


# ### Define parameters
# TODO: Probably want to pull this parameter setting to the top-level script mrtrix_main.py
# Threshold for FA during single fiber voxel mask estimation
absolute_threshold_value = 0.7
# Number of erosion passes during single voxel fiber mask estimation
number_of_passes = 1

# Define default filenames
fileNames = {
    'gradFile': 'btable.b',
    'dtFile': 'dt.mif',
    'faFile': 'fa.mif',
    'evFile': 'ev.mif',
    'singleFiberFile': 'singleFibers.mif',
    'responseFile': 'response.txt',
    'csdFile': 'CSD.mif'
            }


# ### Define input and outpute nodes

inputNode = Node(IdentityInterface(fields = ['dwi_file',
                                            'bval_file',
                                            'bvec_file',
                                            'wmmask',
                                            'tracking_dir']), 
                 name = 'input_node')


outputNode = Node(IdentityInterface(fields = ['spherical_harmonics_image']), 
                  name = 'output_node')


# ### Utility functions

def fileNameBuilder(path, fname):
    return path + '/' + fname

# Get the maximum possbile spherical-harmonics order 
# out of the data because mrtrix itself is not capable of doing this safely


def estimateMaxHarmOrder(bval_file):
    import numpy as np

    with open(bval_file, 'r') as f:
        tmp = f.read()
        tmp = np.asarray(tmp.split())
    
    f.close()

    theRowCount = np.count_nonzero(tmp > '0')

    theOrder = 0.5 * (np.sqrt(8 * theRowCount + 1) - 3)

    # Return the near-lowest multiple of 2
    return theOrder - (theOrder % 2)


def multiplyMRTrix(in1, in2, out_file):
    from nipype.interfaces import mrtrix

    mrmult = mrtrix.MRMultiply()
    mrmult.inputs.in_files = [in1, in2]
    mrmult.inputs.out_filename = out_file

    return out_file


# ### MRTrix specific preprocessing

# First convert the FSL-like input of bval and bvec into mrtrix format
fsl2mrtrixNode = Node(mrtrix.FSL2MRTrix(), name = 'fsl_2_mrtrix')

# Diffusion tensor images
dwi2tensorNode = Node(mrtrix.DWI2Tensor(), name = 'dwi_2_tensor')


# Fractional anisotropy (FA) map
tensor2faNode = Node(mrtrix.Tensor2FractionalAnisotropy(), name = 'tensor_2_FA')

# Remove noisy background by multiplying the FA Image with the binary brainmask
mrmultNode = Node(Function(input_names = ['in1', 'in2', 'out_file'],
                           output_names = ['out_file'],
                           function = multiplyMRTrix),
                  name = 'mrmult')

# Eigenvector (EV) map
tensor2vectorNode = Node(mrtrix.Tensor2Vector(), name = 'tensor_2_vector')

# Scale the EV map by the FA Image
scaleEvNode = mrmultNode.clone('scale_ev')

# Mask of single-fibre voxels
erodeNode = Node(mrtrix.Erode(), name = 'erode_wmmask')
erodeNode.inputs.number_of_passes = number_of_passes

cleanFaNode = mrmultNode.clone('multiplyFA_Mask')

thresholdFANode = Node(mrtrix.Threshold(), name = 'threshold_FA')
thresholdFANode.inputs.absolute_threshold_value = absolute_threshold_value

# Response function coefficient
estResponseNode = Node(mrtrix.EstimateResponseForSH(), name = 'estimate_deconv_response')

# CSD computation
csdNode = Node(mrtrix.ConstrainedSphericalDeconvolution(), name = 'compute_CSD')


# ### Connect the Nodes in the workflow

wf = Workflow(name = 'MRTrix_preproc')

wf.connect([
        (inputNode, fsl2mrtrixNode, [('bval_file', 'bval_file'),
                                        ('bvec_file', 'bvec_file'),
                                        (('tracking_dir', fileNameBuilder, fileNames['gradFile']), 'out_encoding_file')]),
        (inputNode, dwi2tensorNode, [('dwi_file', 'in_file'),
                                    (('tracking_dir', fileNameBuilder, fileNames['dtFile']), 'out_filename')]),
        (fsl2mrtrixNode, dwi2tensorNode, [('encoding_file', 'encoding_file')]),
        (dwi2tensorNode, tensor2faNode, [('tensor', 'in_file')]),
        (inputNode, tensor2faNode, [(('tracking_dir', fileNameBuilder, fileNames['faFile']), 'out_filename')]),
        (tensor2faNode, mrmultNode, [('FA', 'in1')]),
        (inputNode, mrmultNode, [('wmmask', 'in2')]),
        (inputNode, mrmultNode, [(('tracking_dir', fileNameBuilder, fileNames['faFile']), 'out_file')]),
        (dwi2tensorNode, tensor2vectorNode, [('tensor', 'in_file')]),
        (inputNode, tensor2vectorNode, [(('tracking_dir', fileNameBuilder, fileNames['evFile']), 'out_filename')]),
        (tensor2vectorNode, scaleEvNode, [('vector', 'in1')]),
        (mrmultNode, scaleEvNode, [('out_file', 'in2')]),
        (inputNode, scaleEvNode, [(('tracking_dir', fileNameBuilder, fileNames['evFile']), 'out_file')]),
        (inputNode, erodeNode, [('wmmask', 'in_file'),
                               (('tracking_dir', fileNameBuilder, 'tmp_singleFiberFile.mif'), 'out_filename')]),
        (erodeNode, cleanFaNode, [('out_file', 'in1')]),
        (inputNode, cleanFaNode, [(('tracking_dir', fileNameBuilder, 'tmp_singleFiberFile.mif'), 'out_file')]),
        (mrmultNode, cleanFaNode, [('out_file', 'in2')]),
        (cleanFaNode, thresholdFANode, [('out_file', 'in_file')]),
        (inputNode, thresholdFANode, [(('tracking_dir', fileNameBuilder, fileNames['singleFiberFile']),
                                       'out_filename')]),
        (inputNode, estResponseNode, [(('bval_file', estimateMaxHarmOrder), 'maximum_harmonic_order'),
                                     (('tracking_dir', fileNameBuilder, fileNames['responseFile']), 'out_filename'),
                                     ('dwi_file', 'in_file')]),
        (fsl2mrtrixNode, estResponseNode, [('encoding_file', 'encoding_file')]),
        (thresholdFANode, estResponseNode, [('out_file', 'mask_image')]),
        (inputNode, csdNode, [(('tracking_dir', fileNameBuilder, fileNames['csdFile']), 'out_filename'),
                             ('dwi_file', 'in_file'),
                             (('bval_file', estimateMaxHarmOrder), 'maximum_harmonic_order'),
                             ('wmmask', 'mask_image')]),
        (fsl2mrtrixNode, csdNode, [('encoding_file', 'encoding_file')]),
        (estResponseNode, csdNode, [('response', 'response_file')]),
        (csdNode, outputNode, [('spherical_harmonics_image', 'spherical_harmonics_image')])
        ])

