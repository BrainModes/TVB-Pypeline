
# coding: utf-8

# # Pipeline Main Workflow

# This file defines the main (scaffold) workflow fo the pipeline.
# The tractorgraphy building block is intended to be freely exchangeable

# Add the path of our functions to the syspath
import sys, os, getopt
funcPath = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
sys.path.append(funcPath)

# Start of Debug-Stuff
# from nipype import config
# config.enable_debug_mode()
# End of Debug-Stuff

from nipype import Node, Workflow, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
import bm_functions as brainmodes
import logging
from multiprocessing import cpu_count


# ### Inputs parameters
subject_id = None
subject_folder = None
structural_rawdata = None
diffusion_rawdata = None
functional_rawdata = None

usageString = '''+++ TVB Automated Processing Pipeline +++

    This script invokes the TVB automated processing Pipeline.
    Usage:  -s <SUBJECT-ID> -r <SUBJECT-DIR> -a <T1-DATA-DIR> -d <DIFFUSION-DATA-DIR>

    Obligatory inputs are:
        -s --sub-id <SUBJECT-ID>            :   The Identifier of the Subject
        -r --sub-dir <SUBJECT-DIR>          :   The absolute path to the folder where you want the results to be stored
        -a --structural-rawdata <IMG-PATH>  :   Absolute path to your structural T1 anatomical data
        -d --diffusion-rawdata <IMG-PATH>   :   Absolute path to your diffusion weighted MRI data

    Optional inputs are:
        -f --functional-rawdata <IMG-PATH>  :   Absolute path to your functional MRI data (BOLD)

    >> Important Note on Image-Data:
    You can input your data in various file formats, e.g. DICOM or NifTi.
    If you have a series of Images for the same modality, e.g. 192 DICOM images for your T1 data,
    just take the path to the first image in the series and make sure that there are no images
    from different series inside the very same folder.
                '''
try:
    opts, args = getopt.getopt(sys.argv[1:], 'hs:d:r:a:f:', ['sub-id=', 'sub-dir=', 'structural-rawdata=',
                                                             'diffusion-rawdata=', 'functional-rawdata=', 'help'])
except getopt.GetoptError:
    sys.exit( "Invalid argument used! See --help for usage of this script!") # Exit script with error
for opt, arg in opts:
    if opt in ('--help', '-h'):
        sys.exit(usageString)
    elif opt in ('-s', '--sub-id'):
        subject_id = arg
    elif opt in ('-r', '--sub-dir'):
        subject_folder = arg
    elif opt in ('-a', '--structural-rawdata'):
        structural_rawdata = arg
    elif opt in ('-d', '--diffusion-rawdata'):
        diffusion_rawdata = arg
    elif opt in ('-f', '--functional-rawdata'):
        functional_rawdata = arg
# Check if all the obligatory inputs are set
if subject_id is None:
    sys.exit('ERROR: Subject ID must not be empty!')
elif subject_folder is None or not os.access(subject_folder, os.W_OK):
    sys.exit('ERROR: No Results-Path set or not writeable!')
elif structural_rawdata is None or not os.path.isfile(structural_rawdata):
    sys.exit('ERROR: No T1 filepath set or file doesnt exist!')
elif diffusion_rawdata is None or not os.path.isfile(diffusion_rawdata):
    sys.exit('ERROR: No diffusion MRI filepath set or file doesnt exist!')
elif functional_rawdata is not None and not os.path.isfile(structural_rawdata):
    sys.exit('ERROR: Functional rawdata-file doesnt exist!')


# ### Setup
inputNode = Node(IdentityInterface(fields = ['subject_folder', 'subject_id', 'structural_rawdata',
                                             'diffusion_rawdata', 'functional_rawdata']),
                name = 'input_node')

inputNode.inputs.subject_folder = subject_folder
inputNode.inputs.subject_id = subject_id
inputNode.inputs.structural_rawdata = structural_rawdata
inputNode.inputs.diffusion_rawdata = diffusion_rawdata
inputNode.inputs.functional_rawdata = functional_rawdata

# ### Logging
logging.basicConfig(filename = subject_folder + subject_id + '/pipeline.log', level=logging.DEBUG)


# ### Utiliy functions
def roiRange(number_of_rois):
    return range(1,number_of_rois + 1)


# ## Preprocessing
import preprocSub as preprocessing


# ## Functional processing
import feat as funcProc


# ## Tractography-Mask generation
maskGenNode = Node(Function(input_names = ['subPath',
                                           'mask_output_folder',
                                           'wmoutline2diff_1mm',
                                           'wmparc2diff_1mm',
                                           'seedsPerVoxel'],
                           output_names = ['seed_target_masks', 'seed_count',
                                           'number_of_rois', 'affine_matrix',
                                          'wmborder_data'],
                           function = brainmodes.generate_masks),
                   name = 'generate_masks')

maskGenNode.inputs.seedsPerVoxel = 200


# ## Tracking
import mrtrix as mrtrix


# ## Connectivity
connectivityRowNode = MapNode(Function(input_names = ['roi', 
                                                      'subid', 
                                                      'affine_matrix',
                                                      'wmborder', 
                                                      'tracksPath', 'track_files'],
                                   output_names = ['SC_cap_row_filename', 'SC_dist_row_filename'],
                                   function = brainmodes.compute_connectivity_row),
                            name = 'comp_SC_row',
                            iterfield = ['roi'])

aggregateConnectivityNode = Node(Function(input_names = ['sub_id',
                                                         'wmborder',
                                                         'tracksPath',
                                                         'cap_row_files',
                                                         'dist_row_files',
                                                         'steplength'],
                                         output_names = ['SC_matrix_filename_matlab',
                                                         'SC_matrix_filename_json', 'SC_matrix_voxelwise_filename',
                                                         'SC_matrix_voxelwise_distinctConnections_capacity_filename',
                                                         'SC_matrix_voxelwise_distinctConnections_distance_filename'],
                                         function = brainmodes.aggregate_connectivity),
                                name = 'aggregate_SC')

# Incorporate the step size used during tracking (default 0.2 mm)
aggregateConnectivityNode.inputs.steplength = mrtrix.mrtrix_tracking.trackingNode.inputs.step_size


# ## TVB formatting
# TODO: Implement this


# ## Build the Workflow
wf = Workflow(name = 'TVB_pipeline', base_dir = subject_folder + subject_id + '/')

# This is the part where we plug all the different toolboxes together (for sytax see Nipype Docs)

# Connect the Input to the Preprocessing step
wf.connect([(inputNode, preprocessing.wf, [('subject_folder', 'input_node.subject_folder'),
                                           ('subject_id', 'input_node.subject_id'),
                                           ('structural_rawdata', 'structural_rawdata'),
                                           ('diffusion_rawdata', 'diffusion_rawdata')])])

# Mask Generation
wf.connect([(preprocessing.wf, maskGenNode, [('output_node.subPath', 'subPath'),
                                            ('output_node.mask_folder', 'mask_output_folder'),
                                            ('output_node.wmoutline2diff_1mm', 'wmoutline2diff_1mm'),
                                            ('output_node.wmparc2diff_1mm', 'wmparc2diff_1mm')])])

# fMRI Processing
if functional_rawdata is not None:  # Connect only if fMRI input data was set
    wf.connect([(preprocessing.wf, funcProc.wf, [('output_node.brainmask', 'inputNode.brainmask'),
                                                ('output_node.aparc+aseg', 'inputNode.parcellation_mask')]),
                (inputNode, funcProc.wf, [('subject_folder', 'inputNode.subject_folder'),
                                          ('functional_rawdata', 'inputNode.raw_files')
                                        ('subject_id', 'inputNode.subID')])])

# MRTrix TRacking
wf.connect([(maskGenNode, mrtrix.mrtrix_main.wf, [('seed_target_masks', 'input_node.seed_target_masks'),
                                                    ('seed_count', 'input_node.seed_count')]),
             (preprocessing.wf, mrtrix.mrtrix_main.wf, [('output_node.bval_file', 'input_node.bval_file'),
                                                        ('output_node.bvec_file', 'input_node.bvec_file'),
                                                        ('output_node.dwi_file', 'input_node.dwi_file'),
                                                        ('output_node.trackingFolder', 'input_node.tracking_dir'),
                                                        ('output_node.tracks_folder', 'input_node.tracks_dir'),
                                                        ('output_node.highresWmMask', 'input_node.wmmask_1mm'),
                                                        ('output_node.lowresWmMask', 'input_node.wmmask')])])


# BrainModes Connecticitivy Row
wf.connect([(maskGenNode, connectivityRowNode, [(('number_of_rois', roiRange), 'roi'),
                                                ('affine_matrix', 'affine_matrix'),
                                                ('wmborder_data', 'wmborder')]),
            (inputNode, connectivityRowNode, [('subject_id', 'subid')]),
            (preprocessing.wf, connectivityRowNode, [('output_node.tracks_folder', 'tracksPath')]),
            (mrtrix.mrtrix_main.wf, connectivityRowNode, [('output_node.trk_files', 'track_files')])])

# BrainModes Connecticity Aggregation
wf.connect([(inputNode, aggregateConnectivityNode, [('subject_id', 'sub_id')]),
            (maskGenNode, aggregateConnectivityNode, [('wmborder_data', 'wmborder')]),
            (preprocessing.wf, aggregateConnectivityNode, [('output_node.tracks_folder', 'tracksPath')]),
            (connectivityRowNode, aggregateConnectivityNode, [('SC_cap_row_filename', 'cap_row_files'),
                                                                ('SC_dist_row_filename', 'dist_row_files')])])


# ## Draw the Graph
wf.write_graph(subject_folder + subject_id + "/TVB_workflow_graph.dot", graph2use = 'colored')
# from IPython.display import Image
# Image(filename="./TVB_workflow_graph.dot.png")

# ## Run the Workflow
wf.run(plugin='MultiProc', plugin_args={'n_procs': cpu_count()})
#wf.run(plugin='OAR', plugin_args={'oarsub_args': '-l walltime=04:00:00'})
wf.run()




