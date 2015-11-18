
# coding: utf-8

# # MRTRix 0.2 main workflow

# This workflow connects the sub-workflows related to dwMRI processing and tractography using MRTrix 0.2,
# thus serving as a so-called "scaffold workflow".

# In[12]:

from nipype import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function

import mrtrix_preproc as preproc
import mrtrix_tracking as tracking

import logging


# ### Start the logging

# In[3]:

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


# ### Define Input- and Outputnode

# In[5]:

inputNode = Node(IdentityInterface(fields = ['dwi_file',
                                            'b0_image',
                                            'bval_file',
                                            'bvec_file',
                                            'wmmask',
                                            'wmmask_1mm',
                                            'tracking_dir']), 
                 name = 'input_node')


# ### Define the Workflow

# In[6]:

wf = Workflow(name = 'scaffoldWorkflow')

wf.connect([
    (inputNode, preproc.inputNode, [('b0_image', '')])
    ])


# In[ ]:



