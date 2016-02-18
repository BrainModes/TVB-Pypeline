# TVB-Pypeline
Work in Progress!

At the moment this pipeline simply tries so resemble our other, shell-script based, pipeline at http://github.com/BrainModes/TVB-empirical-data-pipeline

The sripts are all included in the notebooks folder, currently as ipython notebooks and py-files.

TVB_pipeline is the "master-control-script" which is later intended to be run via the console
(i.e. python TVB_pipeline.py )

To run it on a specific cluster architecture, simply edit the plugin-type in the master control script

----------

## TODO-List:
+ Write a new Documentation!
+ Implement fMRI processing based on the code here: https://github.com/BrainModes/TVB-empirical-data-pipeline/blob/NSG/fmriFC.sh
+ Implement Tractography Thresholding into the MRTrix module. Possibly trying to include the method described in [Morris et al. (2008)](http://www.sciencedirect.com/science/article/pii/S1053811908007301)
++ Alternatively one could also dig up the old hard-threshold code since the short-range tracking-flares are rendered meaningles anyway by our aggregation method!
+ Make the file-sorting of the user-data more sophisticated. This means that the pipeline should be able to somehow recognize which kinds of data-sets (e.g. fMRI, T1, dwMRI) is included in the user data and then route the particular folder-paths onto the corresponding processing-nodes inside the pipeline. This might be achieved through using nipype's [SelectFiles interface](http://nipy.org/nipype/users/select_files.html)
+ Include some example workflows for different cluster scenarios, realized through e.g. controll-scripts written in BASH
+ Re-Implement Multishell-Tracking using FSLs bedpostx as in https://github.com/BrainModes/TVB-empirical-data-pipeline/tree/multiShell
