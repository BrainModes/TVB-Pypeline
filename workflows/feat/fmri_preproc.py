
# coding: utf-8

# ### Handle the imports
from nipype import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces import freesurfer, fsl
# from nipype.interfaces.io import DataFinder

from bm_functions import compute_functional_connectivity, mri_convert_bm
from default_feat_config import gen_default_feat_config

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


# ### Define Input- and Output-Node

# Inputs:
# -----------
# raw_files : Path-string to the folder containing the raw DICOM files
# subject_folder : Path-String to the folder containg the already processed subject data i.e. the path where 
#                the results of this step are also stored
inputNode = Node(IdentityInterface(fields = ['subID',
                                            'raw_files',
                                             'subject_folder',
                                            'parcellation_mask',
                                            'brainmask']), 
                 name = 'inputNode')


outputNode = Node(IdentityInterface(fields = ['mat_file']),
                  name = 'output_node')


# ### Define Filenames

fileNames = {
    'bold_file': 'bold.nii.gz',
    'parcellation_2_func': 'parc_2_func.nii.gz',
    'func_2_anat': 'exfunc2anat_6DOF.nii.gz',
    'func_2_anat_mat': 'exfunc2anat_6DOF.mat',
    'anat_2_func_mat': 'anat2exfunc_6DOF.mat',
    'segstat_sum_file': 'segstat_summary.txt',
    'avgwf_file': 'segstat_average_ROI_timeseries.dat'
}


# ### Utility Functions

def folder_maker(path_name, folder_name=None):
    
    if folder_name is None:
        folder_name = 'bold'
    
    import os
    if not os.path.exists(path_name + '/' + folder_name):
        os.makedirs(path_name + '/' + folder_name)
        
    return path_name + '/' + folder_name + '/'

folderMaker = Node(Function(input_names = ['path_name', 'folder_name'],
                                output_names = ['folder_path'],
                                function = folder_maker),
                                name = 'folder_maker')


def fileNameBuilder(path, fname):
    return path + fname       


def selectFromList(inList, index):
    try:
        return inList[index]
    except TypeError:
        return inList

# ### Convert Images from various formats to NifTi
# rawFinderNode = Node(DataFinder(match_regex = '.*\.dcm'), name = 'DICOM_Finder')

#convertNode = Node(freesurfer.preprocess.MRIConvert(), name = 'DICOM2Nii')
#convertNode.inputs.out_type = 'niigz'
#convertNode.inputs.out_orientation = 'RAS'

convertNode = Node(Function(input_names = ['in_file', 'out_file'],
                            output_names = ['out_file'],
                            function = mri_convert_bm),
                   name = 'DICOM_to_NifTi')


# ### Run FSLs feat
def run_feat(bold_file, bold_folder, brainmask_file, gen_feat_config):
    from nipype.interfaces.fsl import ImageStats, FEAT, Info
    # from bm_functions import gen_default_feat_config
    from numpy import shape
    from textwrap import dedent

    fslFilename = bold_folder + 'feat.fsf'

    # Get the number of voxels in the 4D file
    statComp = ImageStats()
    statComp.inputs.in_file = bold_file
    statComp.inputs.op_string = '-v'

    numVox = int(statComp.run().outputs.out_stat[0])

    # Get the number of raw volumes
    statComp.inputs.split_4d = True

    numVol = shape(statComp.run().outputs.out_stat)[0]

    # Generate the file
    standard_T1_brain = Info.standard_image('MNI152_T1_2mm_brain')
    theString = gen_feat_config(bold_folder, bold_file, brainmask_file, standard_T1_brain, numVox, numVol)
    with open(fslFilename,'w') as out_file:
        out_file.write(dedent(theString))
    out_file.close()   

    # Run feat using the previously manipulated config
    runFeat = FEAT(fsf_file = fslFilename)
    # Run and pass back the foldername
    return runFeat.run().outputs.feat_dir


featNode = Node(Function(input_names=['bold_file', 'bold_folder', 'brainmask_file', 'gen_feat_config'],
                        output_names=['feat_dir'],
                        function=run_feat),
                name = 'FSL_feat')

featNode.inputs.gen_feat_config = gen_default_feat_config


# ## Generate parcellated ROI-Timeseries

# Register example-func to freesurfer brainmask
exfunc2anat = Node(fsl.FLIRT(bins=256, searchr_x=[90,90], searchr_y=[90,90], searchr_z=[90,90],
                       cost='corratio', interp='trilinear', dof=6),
                  name = 'Func_2_Anat')

# invert transformation
invt = Node(fsl.ConvertXFM(invert_xfm=True),
            name = 'invert_transf')

# transform roimask to functional space using FLIRT (using Nearest Neighbor Interpolation for roimask)
roimask2func = Node(fsl.FLIRT(padding_size=0, interp='nearestneighbour', apply_xfm=True),
                    name = 'roimask_2_func')

# Export average region time-series
ss = Node(freesurfer.SegStats(), name = 'SegStats')


# ### Preprocess the data before computing the FC

def segstat_shaping(aparc_stats):
    import os, re
    
    # Remove all comment lines from the files (important for later MATLAB/OCTAVE import!)
    clearedFileName = os.path.dirname(aparc_stats) + '/aparc_stats_tmp.txt'

    with open(clearedFileName,'w') as out_file:
        with open(aparc_stats,'r') as in_file:
            for line in in_file:
                if line.startswith('#'):
                    line = ''
                else:
                    line = re.sub('Seg', '', line.strip()) + '\n'
                out_file.write(line)
    out_file.close()
    in_file.close()
    
    return clearedFileName

segstatPost = Node(Function(input_names = ['aparc_stats'],
                           output_names = ['clearedFileName'],
                           function = segstat_shaping),
                  name = 'segstat_Postprocesser')


# ### Compute the FC
compFCNode = Node(Function(input_names = ['path', 'subName', 'avgwf_txt_file', 'summary_file_cleared'],
                          output_names = ['matfile_name'],
                          function = compute_functional_connectivity), 
                 name = 'compute_FC')


# ## Define the Workflow
wf = Workflow('fMRI_Processing')

# wf.connect([(inputNode, rawFinderNode, [('raw_files', 'root_paths')])])

wf.connect([(inputNode, folderMaker, [('subject_folder', 'path_name')])])


# wf.connect([(rawFinderNode, convertNode, [(('out_paths', selectFromList, 0), 'in_file')]),
wf.connect([(inputNode, convertNode, [(('raw_files', selectFromList, 0), 'in_file')]),
           (folderMaker, convertNode, [(('folder_path', fileNameBuilder, fileNames['bold_file']), 'out_file')])])

wf.connect([(convertNode, featNode, [('out_file', 'bold_file')]),
           (folderMaker, featNode, [('folder_path', 'bold_folder')]),
           (inputNode, featNode, [('brainmask', 'brainmask_file')])])

wf.connect([(featNode, exfunc2anat, [(('feat_dir', fileNameBuilder, 'mean_func.nii.gz'), 'in_file')]),
           (inputNode, exfunc2anat, [('brainmask', 'reference')]),
           (folderMaker, exfunc2anat, [(('folder_path', fileNameBuilder, fileNames['func_2_anat']), 'out_file'),
                                      (('folder_path', fileNameBuilder, fileNames['func_2_anat_mat']), 'out_matrix_file')])])

wf.connect([(exfunc2anat, invt, [('out_matrix_file', 'in_file')]),
          (folderMaker, invt, [(('folder_path', fileNameBuilder, fileNames['anat_2_func_mat']), 'out_file')])])

wf.connect([(inputNode, roimask2func, [('parcellation_mask', 'in_file')]),
           (invt, roimask2func, [('out_file', 'in_matrix_file')]),
           (featNode, roimask2func, [(('feat_dir', fileNameBuilder, 'mean_func.nii.gz'), 'reference')]),
           (folderMaker, roimask2func, [(('folder_path', fileNameBuilder, fileNames['parcellation_2_func']), 'out_file')])])

wf.connect([(roimask2func, ss, [('out_file', 'segmentation_file')]),
           (folderMaker, ss, [(('folder_path', fileNameBuilder, fileNames['avgwf_file']), 'avgwf_txt_file')]),
           (featNode, ss, [(('feat_dir', fileNameBuilder, 'filtered_func_data.nii.gz'), 'in_file')]),
           (folderMaker, ss, [(('folder_path', fileNameBuilder, fileNames['segstat_sum_file']), 'summary_file')])])

wf.connect([(ss, segstatPost, [('summary_file', 'aparc_stats')])])

wf.connect([(segstatPost, compFCNode, [('clearedFileName' , 'summary_file_cleared')]),
           (ss, compFCNode, [('avgwf_txt_file', 'avgwf_txt_file')]),
           (inputNode, compFCNode, [('subID', 'subName')]),
           (folderMaker, compFCNode, [('folder_path', 'path')])])

wf.connect([(compFCNode, outputNode, [('matfile_name', 'mat_file')])])

## Draw the Graph
#wf.write_graph("/Users/srothmei/Desktop/charite/temp/TVB_workflow_graph.dot", graph2use = 'exec', simple_form=False)
#from IPython.display import Image
#Image(filename="./TVB_workflow_graph.dot.png")




