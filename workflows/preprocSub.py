
# coding: utf-8

# # Pipeline Preprocessing Workflow

# This script handles all the preprocessing of the rawdata for structural MRI and also diffusion weighted MRI
# The tools used in this workflow are <b>FSL</b>, <b>FREESURFER</b> & <b>Dicom2Nii</b>

# #### Handle the imports
from nipype import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces import freesurfer, fsl
from nipype.interfaces.dcm2nii import Dcm2nii
from bm_functions import mri_convert_bm

import logging
from multiprocessing import cpu_count


# #### Start the logging
logger = logging.getLogger('workflow')
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


# ### Utility Nodes
def extractB0(dwMriFile):
    from nipype import utils as nputils
    from nipype.interfaces import fsl
    import os, shutil
    # This function is used to extract the b0 image out of the 4D series of diffusion images by splitting the
    # series, copying the first image and merging the images together into a single file afterwards
    pth, fname, ext = nputils.filemanip.split_filename(dwMriFile)
    currPth = os.path.realpath('.')

    # Change dir because FSL split takes no input for storing things somewhere else....
    os.chdir(pth)    

    splitter = fsl.utils.Split(dimension = 't', out_base_name = 'tmp_', in_file = dwMriFile)
    res = splitter.run()
    b0 = res.outputs.out_files[0]

    # Copy the b0 image and rename it accordingly
    shutil.copy(b0, 'no_diffusion.nii.gz')
    b0 = os.path.abspath('no_diffusion.nii.gz')

    # Delete the sliced files...
    for f in res.outputs.out_files:
        os.remove(f)

    # Change back into previous directory
    os.chdir(currPth)
    
    return b0


def selectFromList(inList, index):
    try:
        return inList[index]
    except TypeError:
        return inList


# ## Set parameters and build variables
reconallFolderName = 'recon_all'  # Define what the output folder of recon-all should be named
# Predefine some filenames
fileNames = {'wmSurf_lh': 'lh_white.nii.gz',
             'wmSurf_rh': 'rh_white.nii.gz',
             'wmSurf': 'wm_outline.nii.gz',
             'wmoutline2diff_1mm': 'wmoutline2diff_1mm.nii.gz',
             'wmoutline2diff': 'wmoutline2diff.nii.gz',
             'wmparc2diff_1mm': 'wmparc2diff_1mm.nii.gz',
             'aparc+aseg': 'aparc+aseg.nii.gz',
             'aparc+aseg2diff_1mm': 'aparc_aseg2diff_1mm.nii.gz',
             'aparc+aseg2diff': 'aparc_aseg2diff.nii.gz',
             'lowresWmMask': 'wmmask.nii.gz',
             'highresWmMask': 'wmmask_1mm.nii.gz',
             'bval_file': 'bvals.dat',
             'bvec_file': 'bvecs.dat',
             'dwi_file': 'dwi.nii.gz',
             'no_diffusion_image': 'lowb.nii.gz',
             'brainmask': 'brainmask.nii.gz'}


def fileNameBuilder(path, fname):
    return path + '/' + fname


# def rawdataChecker(input_file):
#     # If the input is a single DCM-file instead of a multi-dim-NifTI, we have to fetch all the other files in the series
#     if input_file.endswith('.dcm'):
#         from nipype.interfaces.io import DataFinder
#         from os import path
#         from nipype import Node
#
#         # Setup a datafinder to find the paths to the specific DICOM files
#         t1FinderNode = Node(DataFinder(), name = 't1Finder')
#         t1FinderNode.inputs.match_regex = '.*\.dcm'
#         t1FinderNode.inputs.root_paths = path.split(input_file)[0]
#
#         return t1FinderNode.run().outputs.out_paths
#     else:
#         return input_file  # If other datatype just return the same path


def pathBuilder(subject_folder, subject_id):
    subPath = subject_folder + '/' # Build full path to subject folder
    
    def makeMyDir(dirName):
        import os
        if not os.path.exists(dirName):
            os.makedirs(dirName)
     
    # Path definitions
    dwiPreprocFolder = subPath + '/diff_processed/'
    trackingFolder = subPath + '/tractography/'
    calc_images = subPath + '/calc_images/'
    mask_folder = trackingFolder + '/masks/'
    tracks_folder = trackingFolder + '/tracks/'
    
    # Make folders
    makeMyDir(dwiPreprocFolder)
    makeMyDir(trackingFolder)
    makeMyDir(calc_images)
    
    makeMyDir(mask_folder)
    makeMyDir(tracks_folder)
    
    # RawData Structure
    rawdataFolder = subPath + '/RAWDATA' # Define the path to the folder holding the rawdata dicom-files
    T1RawFolder = rawdataFolder + '/MPRAGE/' # The T1 rawdata folder
    dwiRawFolder = rawdataFolder + '/DTI' # The dwMRI rawdata folder
    fmriRawFolder = rawdataFolder + '/BOLD-EPI/' # The fMRI rawdata folder
    
    return subPath, rawdataFolder, T1RawFolder, dwiRawFolder, fmriRawFolder, dwiPreprocFolder, trackingFolder, calc_images, mask_folder, tracks_folder

pathBuildingNode = Node(Function(input_names = ['subject_folder', 'subject_id'],
                                output_names = ['subPath', 'rawdataFolder', 'T1RawFolder',
                                                'dwiRawFolder', 'fmriRawFolder', 
                                                'dwiPreprocFolder', 'trackingFolder', 'calc_images',
                                               'mask_folder', 'tracks_folder'],
                                function = pathBuilder),
                        name='pathBuilder')


# ### Define Inputnode and Outputnode
inputNode = Node(IdentityInterface(fields=['subject_id', 'subject_folder', 'structural_rawdata', 'diffusion_rawdata',
                                           'b_val', 'b_vec']),
                 mandatory_inputs = True,
                 name='input_node')


mergedOutputs = pathBuildingNode.outputs.copyable_trait_names()
mergedOutputs.extend(fileNames.keys())

outputNode = Node(IdentityInterface(fields = mergedOutputs), mandatory_inputs = False, name = 'output_node')


# ### Structural Data (T1) preprocessing

# Set recon-all parameters
reconallNode = Node(freesurfer.preprocess.ReconAll(), name = 'reconall')
# reconallNode.inputs.T1_files = firstFile
# reconallNode.inputs.subjects_dir = subPath
reconallNode.inputs.subject_id = reconallFolderName
reconallNode.inputs.directive = 'all'
reconallNode.inputs.openmp = cpu_count()
# reconallNode.inputs.args = '-notal-check'

# OAR Workaround
# reconallNode.plugin_args = {'overwrite': True, 'oarsub_args': '-l nodes=1,walltime=16:00:00'}

# Convert the T1 mgz image to nifti format for later usage
# mriConverter = Node(freesurfer.preprocess.MRIConvert(), name = 'convertAparcAseg')
# mriConverter.inputs.out_type = 'niigz'
# mriConverter.inputs.out_orientation = 'RAS'
mriConverter = Node(Function(input_names = ['in_file', 'out_file'],
                            output_names = ['out_file'],
                            function = mri_convert_bm),
                   name = 'convertAparcAseg')

# Convert the Brainmask file
# brainmaskConv = Node(freesurfer.preprocess.MRIConvert(), name = 'convertBrainmask')
# brainmaskConv.inputs.out_type = 'niigz'
# brainmaskConv.inputs.out_orientation = 'RAS'
brainmaskConv = mriConverter.clone('convertBrainmask')


# ### Diffusion Data (dwMRI) preprocessing
# First extract the diffusion vectors and the pulse intensity (bvec and bval)
# Use dcm2nii for this task
dcm2niiNode = Node(Dcm2nii(), name = 'dcm2niiAndBvecs')
dcm2niiNode.inputs.gzip_output = True
dcm2niiNode.inputs.date_in_filename = False
dcm2niiNode.inputs.events_in_filename = False


# Extract the first image of the DTI series i.e. the b0 image
extrctB0Node = Node(Function(input_names = ['dwMriFile'], output_names = ['b0'],
                             function = extractB0), name = 'Extract_b0')


# Perform the registration between subject T1 space and dwMRI space
bbregNode = Node(freesurfer.preprocess.BBRegister(), name = 'BBRegister')
bbregNode.inputs.init = "fsl"
bbregNode.inputs.contrast_type = "t2"
bbregNode.inputs.epi_mask = True
bbregNode.inputs.out_fsl_file = True
bbregNode.inputs.args = "--tol1d 1e-3"
#bbregNode.inputs.subject_id = reconallFolderName


# ### Surface2Vol

# Transform Left Hemisphere
surf2volNode_lh = Node(freesurfer.utils.Surface2VolTransform(), name = 'surf2vol_lh')
surf2volNode_lh.inputs.hemi = 'lh'
surf2volNode_lh.inputs.mkmask = True
#surf2volNode_lh.inputs.subject_id = reconallFolderName
surf2volNode_lh.inputs.vertexvol_file = 'test'

# Transform right hemisphere
surf2volNode_rh = surf2volNode_lh.clone('surf2vol_rh')
surf2volNode_rh.inputs.hemi = 'rh'

# Merge the hemispheres
mergeHemisNode = Node(fsl.BinaryMaths(), name = 'mergeHemis')
mergeHemisNode.inputs.operation = 'add'
mergeHemisNode.inputs.output_type = 'NIFTI_GZ'


# ### Registration

# Rotate high-res (1mm) WM-border to match dwi data w/o resampling
applyReg_anat2diff_1mm = Node(freesurfer.ApplyVolTransform(), name = 'wmoutline2diff_1mm')
applyReg_anat2diff_1mm.inputs.inverse = True
applyReg_anat2diff_1mm.inputs.interp = 'nearest'
applyReg_anat2diff_1mm.inputs.no_resample = True
applyReg_anat2diff_1mm.inputs.args = '--no-save-reg'

# Rotate high-res (1mm) WM-border to match dwi data with resampling
applyReg_anat2diff = applyReg_anat2diff_1mm.clone('wmoutline2diff')
applyReg_anat2diff.inputs.no_resample = False
applyReg_anat2diff.inputs.interp = 'trilin'

# Filter out low voxels produced by trilin. interp.
thresholdNode = Node(fsl.maths.Threshold(), name = 'remove_interp_residuals')
thresholdNode.inputs.thresh = 0.1
thresholdNode.inputs.output_type = 'NIFTI_GZ'

# Binarize
binarizeNode = Node(fsl.maths.UnaryMaths(), name = 'binarize')
binarizeNode.inputs.operation = 'bin'
binarizeNode.inputs.output_type = 'NIFTI_GZ'

# Rotate high-res (1mm) wmparc to match dwi data w/o resampling
applyReg_wmparc2diff_1mm = applyReg_anat2diff_1mm.clone('wmparc2diff_1mm')

# Rotate high-res (1mm) aparc+aseg to match dwi data w/o resampling
applyReg_aparc2diff_1mm = applyReg_anat2diff_1mm.clone('aparc2diff_1mm')

# Rotate high-res (1mm) aparc+aseg to match dwi data with resampling
applyReg_aparc2diff = applyReg_anat2diff.clone('aparc2diff')
applyReg_aparc2diff.inputs.interp = 'nearest'


# ### Create Whitematter Brainmasks

def extractWhitematter(input_image, wmoutline_image, output_image):
    from nipype.interfaces import fsl

    maths = fsl.MultiImageMaths()
    maths.inputs.in_file = input_image
    maths.inputs.out_file = output_image
    maths.inputs.op_string = ' -add %s -uthr 41 -thr 41'
    maths.inputs.operand_files = wmoutline_image
    maths.run()
    
    #maths.inputs.in_file = output_image
    maths.inputs.op_string = '-uthr 2 -thr 2 -add %s'
    maths.inputs.operand_files = output_image
    maths.run()
    
    maths.inputs.op_string = '-uthr 255 -thr 251 -add %s -add %s -bin'
    maths.inputs.operand_files = [output_image, wmoutline_image]
    maths.run()
    
    return output_image

# Create the Nodes
wmmask_lowres = Node(Function(input_names = ['input_image', 'wmoutline_image', 'output_image'],
                             output_names = ['output_image'],
                             function = extractWhitematter), name = 'extract_WM_lowres')

wmmask_highres = wmmask_lowres.clone('extract_WM_highres')


# TODO: REMOVE ME!
### Debug the crap out of this!
# def myAmazingDebugFunction(T1_files, subjects_dir):
#     reconAllPath = '/home/petra/Simon/TVB-Pypeline/subjects/QL_20120814/recon_all/mri/'
#     T1 = reconAllPath + 'T1.mgz'
#     aparc_aseg = [reconAllPath + 'aparc+aseg.mgz']
#     wmparc = reconAllPath + 'wmparc.mgz'
#     subject_id = 'recon_all'
#     brainmask = reconAllPath + 'brainmask.mgz'
#
#     return subject_id, T1, aparc_aseg, wmparc, brainmask
#
# reconallNode = Node(Function(input_names = ['T1_files', 'subjects_dir'],
#                             output_names = ['subject_id', 'T1', 'aparc_aseg', 'wmparc', 'brainmask'],
#                             function = myAmazingDebugFunction),
#                    name = 'recon_debug_all')


# ### Connect the Nodes
wf = Workflow(name = 'preprocSub')

# Input strings to pathbuilder
wf.connect([(inputNode, pathBuildingNode, [('subject_id', 'subject_id'),
                                          ('subject_folder', 'subject_folder')])])


# T1-Rawdata-path into dataFinder to find T1 DICOMs
#wf.connect(pathBuildingNode, 'T1RawFolder', t1FinderNode, 'root_paths')

# T1 DICOM-paths into recon_all
#wf.connect(t1FinderNode, 'out_paths', reconallNode, 'T1_files')
#wf.connect([(inputNode, reconallNode, [(('structural_rawdata', rawdataChecker), 'T1_files')])])
wf.connect([(inputNode, reconallNode, [('structural_rawdata', 'T1_files')])])

# Subject path into recon-all
wf.connect(pathBuildingNode, 'subPath', reconallNode, 'subjects_dir')

# aparc+aseg into mriConverter
wf.connect([(reconallNode, mriConverter, [(('aparc_aseg', selectFromList, 0), 'in_file')])])
wf.connect([(pathBuildingNode, mriConverter, [(('calc_images', fileNameBuilder, fileNames['aparc+aseg']), 'out_file')])])

# Brainmask converter
wf.connect([(reconallNode, brainmaskConv, [('brainmask', 'in_file')]),
            (pathBuildingNode, brainmaskConv, [(('calc_images', fileNameBuilder, fileNames['brainmask']), 'out_file')])])

# dcm2nii
# wf.connect([(pathBuildingNode, dcm2niiNode, [('dwiRawFolder', 'source_dir'),
# wf.connect([(inputNode, dcm2niiNode, [(('diffusion_rawdata', rawdataChecker), 'source_names')]),
wf.connect([(inputNode, dcm2niiNode, [('diffusion_rawdata', 'source_names')]),
             (pathBuildingNode, dcm2niiNode, [('dwiPreprocFolder', 'output_dir')])])

# B0 extraction
wf.connect(dcm2niiNode, 'converted_files', extrctB0Node, 'dwMriFile')

# Register b0 into T1 space
wf.connect([(extrctB0Node, bbregNode, [('b0', 'source_file')]),
           (pathBuildingNode, bbregNode, [('subPath', 'subjects_dir')]),
            (reconallNode, bbregNode, [('subject_id', 'subject_id')])])

# Tranform lh whitematter surface to voxel-space
wf.connect([(pathBuildingNode, surf2volNode_lh, [('subPath', 'subjects_dir'),
                                                (('calc_images', fileNameBuilder, fileNames['wmSurf_lh']),
                                                 'transformed_file')]),
           (reconallNode, surf2volNode_lh, [('T1', 'template_file'),
                                            ('subject_id', 'subject_id')])])
# Tranform rh whitematter...
wf.connect([(pathBuildingNode, surf2volNode_rh, [('subPath', 'subjects_dir'),
                                                (('calc_images', fileNameBuilder, fileNames['wmSurf_rh']),
                                                 'transformed_file')]),
           (reconallNode, surf2volNode_rh, [('T1', 'template_file'),
                                            ('subject_id', 'subject_id')])])

# Merge the hemispheres
wf.connect([(surf2volNode_lh, mergeHemisNode, [('transformed_file', 'in_file')]),
           (surf2volNode_rh, mergeHemisNode, [('transformed_file', 'operand_file')]),
           (pathBuildingNode, mergeHemisNode, [(('calc_images', fileNameBuilder, fileNames['wmSurf']), 'out_file')])])

# Registrations --__--__--__--
wf.connect([(extrctB0Node, applyReg_anat2diff_1mm, [('b0', 'source_file')]),
           (mergeHemisNode, applyReg_anat2diff_1mm, [('out_file', 'target_file')]),
           (bbregNode, applyReg_anat2diff_1mm, [('out_reg_file', 'reg_file')]),
           (pathBuildingNode, applyReg_anat2diff_1mm,
            [(('calc_images',fileNameBuilder ,fileNames['wmoutline2diff_1mm']), 'transformed_file')])])

wf.connect([(extrctB0Node, applyReg_anat2diff, [('b0', 'source_file')]),
           (mergeHemisNode, applyReg_anat2diff, [('out_file', 'target_file')]),
           (bbregNode, applyReg_anat2diff, [('out_reg_file', 'reg_file')]),
           (pathBuildingNode, applyReg_anat2diff,
            [(('calc_images',fileNameBuilder ,fileNames['wmoutline2diff']), 'transformed_file')])])

wf.connect([(applyReg_anat2diff, thresholdNode, [('transformed_file', 'in_file'),
                                                ('transformed_file', 'out_file')])])
wf.connect([(thresholdNode, binarizeNode, [('out_file', 'in_file'),
                                          ('out_file', 'out_file')])])

wf.connect([(extrctB0Node, applyReg_wmparc2diff_1mm, [('b0', 'source_file')]),
           (reconallNode, applyReg_wmparc2diff_1mm, [('wmparc', 'target_file')]),
           (bbregNode, applyReg_wmparc2diff_1mm, [('out_reg_file', 'reg_file')]),
           (pathBuildingNode, applyReg_wmparc2diff_1mm,
            [(('calc_images',fileNameBuilder ,fileNames['wmparc2diff_1mm']), 'transformed_file')])])

wf.connect([(extrctB0Node, applyReg_aparc2diff_1mm, [('b0', 'source_file')]),
           (reconallNode, applyReg_aparc2diff_1mm, [(('aparc_aseg', selectFromList, 0), 'target_file')]),
           (bbregNode, applyReg_aparc2diff_1mm, [('out_reg_file', 'reg_file')]),
           (pathBuildingNode, applyReg_aparc2diff_1mm,
            [(('calc_images', fileNameBuilder, fileNames['aparc+aseg2diff_1mm']), 'transformed_file')])])

wf.connect([(extrctB0Node, applyReg_aparc2diff, [('b0', 'source_file')]),
           (reconallNode, applyReg_aparc2diff, [(('aparc_aseg', selectFromList, 0), 'target_file')]),
           (bbregNode, applyReg_aparc2diff, [('out_reg_file', 'reg_file')]),
           (pathBuildingNode, applyReg_aparc2diff,
            [(('calc_images',fileNameBuilder ,fileNames['aparc+aseg2diff']), 'transformed_file')])])

# Mask creation
wf.connect([(applyReg_aparc2diff, wmmask_lowres, [('transformed_file', 'input_image')]),
           (binarizeNode, wmmask_lowres, [('out_file', 'wmoutline_image')]),
           (pathBuildingNode, wmmask_lowres, [(('calc_images', fileNameBuilder, fileNames['lowresWmMask']),
                                               'output_image')])])

wf.connect([(applyReg_aparc2diff_1mm, wmmask_highres, [('transformed_file', 'input_image')]),
           (applyReg_anat2diff_1mm, wmmask_highres, [('transformed_file', 'wmoutline_image')]),
           (pathBuildingNode, wmmask_highres, [(('calc_images', fileNameBuilder, fileNames['highresWmMask']),
                                               'output_image')])])

# Now forward all relevant stuff to output node
for i in pathBuildingNode.outputs.copyable_trait_names():
    wf.connect([(pathBuildingNode, outputNode, [(i, i)])])
    
wf.connect([(wmmask_highres, outputNode, [('output_image', 'highresWmMask')]),
           (wmmask_lowres, outputNode, [('output_image', 'lowresWmMask')]),
           (applyReg_aparc2diff, outputNode, [('transformed_file', 'aparc+aseg2diff')]),
           (applyReg_aparc2diff_1mm, outputNode, [('transformed_file', 'aparc+aseg2diff_1mm')]),
           (mergeHemisNode, outputNode, [('out_file', 'wmSurf')]),
           (applyReg_anat2diff_1mm, outputNode, [('transformed_file', 'wmoutline2diff_1mm')]),
           (applyReg_anat2diff, outputNode, [('transformed_file', 'wmoutline2diff')]),
           (applyReg_wmparc2diff_1mm, outputNode, [('transformed_file', 'wmparc2diff_1mm')]),
            (extrctB0Node, outputNode, [('b0', 'no_diffusion_image')]),
            (dcm2niiNode, outputNode, [('converted_files', 'dwi_file')]),
            (brainmaskConv, outputNode, [('out_file', 'brainmask')])
           ])

# If the bvecs/bvals have been set by the user, dont overwrite them by the values of dcm2nii
# Note that both files have to be set thus we only test for one
if inputNode.inputs.b_vec is None:
    wf.connect([(inputNode, outputNode, [('b_val', 'bval_file'),
                                       ('b_vec', 'bvec_file')])])
else:
    wf.connect([(dcm2niiNode, outputNode, [('bvals', 'bval_file'),
                                        ('bvecs', 'bvec_file')])])


#wf.write_graph("preproc_workflow_graph.dot", graph2use = 'exec')
#from IPython.display import Image
#Image(filename="workflow_graph.dot.png")



