# TVB-Pypeline - - Work in Progress!

This project maps our current automatized MRI processing pipeline (http://github.com/BrainModes/TVB-empirical-data-pipeline) 
to Python using Nipype, making the used toolboxes inside easily exchangeable.

For a general overview about the pipeline see [Schirner, Rothmeier et al. (2015)](http://www.sciencedirect.com/science/article/pii/S1053811915002505)

Please note that this pipeline does extensive analysis and is thus computationally heavy. TEsting was carried out on a High-Performance-Clustercomputer using >100 CPU Cores.

----------

## Installation:
The Pipeline uses Nipype which depends mainly on **Python 2.7**. The following list gives an overview about the Python toolboxes which are used in the current state of the Pipeline. See the corresponding Doc-Pages for installation and dependency resolving.
+ [Nipype](http://nipy.org/nipype/users/install.html)
+ [Nibabel](http://nipy.org/nibabel/installation.html#installation)
+ [Dipy](http://nipy.org/dipy/installation.html)

Since Nipype/Python also perform as a wrapper for Toolboxes invoked through the Shell-Interface, you also have to make sure the toolboxes you want to use are installed on your system and their binaries/libs are included in the Shell's searchpath.

For **preprocessing**, the following toolboxes are used:
+ [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
+ [FREESURFER](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)

When it comes to **fiber tractography**, there is a vast number of available tools for that. Their usage also highly depens on how your dwMRI-Data was recorded. One of the main parting points is the number of different diffusion-gradient strengths applied during the measurement (i.e. the number of different **b-values**). If the dataset has only a single value greater than zero, one talks about **single-shell data**. As soon as more than one value (>0) is involved, the data is called **multi-shell data**

Currently, we tested two toolboxes for tractography, one for each of the aforementioned scenarios:
 + [MRTrix 0.2.12](http://jdtournier.github.io/mrtrix-0.2/index.html): Single-Shell Tracking
 + [FSL](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT): Multi-Shell Tracking (Not yet implemented in the Python-Pipeline!)


##### Install the Pipeline
Download the files from the GitHub Repository and unpack the files on your workstation/cluster. 
Next you might want to adapt the plugin-settings to your local architecture (see below).

------
## Running the Pipeline
To run it on a specific cluster architecture, simply edit the plugin-type in the master control script [TVB_pipeline.py](https://github.com/srothmei/TVB-Pypeline/blob/master/workflows/TVB_pipeline.py).
Locate the following code block at the end of the file
```python
# ## Run the Workflow
#wf.run(plugin='MultiProc', plugin_args={'n_procs': cpu_count()})
wf.run(plugin='OAR', plugin_args={'oarsub_args': '-l walltime=04:00:00'})
wf.run()
```
As you can see, plugins are used to handle different situations considering the environment in which the pipeline is intended to be run, e.g. different job schedulers on High-Performance-Clustercomputer or local installations on a multicore workstation.
For an overview about the available plugins see the [Doc-Page about Plugins](http://nipy.org/nipype/users/plugins.html). Since this page is sometimes a bit outdated (e.g. the OAR plugin is not yet listed), see also https://github.com/nipy/nipype/tree/master/nipype/pipeline/plugins


----------

## TODO-List:
+ ~~Write a new Documentation!~~
+ ~~Implement fMRI processing based on the code here: https://github.com/BrainModes/TVB-empirical-data-pipeline/blob/NSG/fmriFC.sh~~
+ Implement Tractography Thresholding into the MRTrix module. Possibly trying to include the method described in [Morris et al. (2008)](http://www.sciencedirect.com/science/article/pii/S1053811908007301). Alternatively one could also dig up the old hard-threshold code since the short-range tracking-flares are rendered meaningles anyway by our aggregation method!
+ Make the file-sorting of the user-data more sophisticated. This means that the pipeline should be able to somehow recognize which kinds of data-sets (e.g. fMRI, T1, dwMRI) is included in the user data and then route the particular folder-paths onto the corresponding processing-nodes inside the pipeline. This might be achieved through using nipype's [SelectFiles interface](http://nipy.org/nipype/users/select_files.html)
+ Include some example workflows for different cluster scenarios, realized through e.g. controll-scripts written in BASH
+ Re-Implement Multishell-Tracking using FSLs bedpostx as in https://github.com/BrainModes/TVB-empirical-data-pipeline/tree/multiShell
+ Implement the formatting of the results into a TVB-ZIP-File as in https://github.com/BrainModes/TVB-empirical-data-pipeline/blob/NSG/matlab_scripts/connectivity2TVBFS.m
