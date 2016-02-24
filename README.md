# TVB-Pypeline - - Work in Progress!

This project maps our current automatized MRI processing pipeline (http://github.com/BrainModes/TVB-empirical-data-pipeline) 
to Python using Nipype, making the used toolboxes inside easily exchangeable.

For a general overview about the pipeline see [Schirner, Rothmeier et al. (2015)](http://www.sciencedirect.com/science/article/pii/S1053811915002505)

![Pipeline Overview](https://github.com/srothmei/TVB-Pypeline/blob/master/doc/overview.png "Graphical Pipeline-Overview")

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

## Preparing your rawdata
Looking at the TODO-List in the bottom-section of this manual, you can see that the organization of the users raw-data is still a bit inflexible considering the fact that the pipeline requires a certain folder-schema. Currently, you need to **precisely stick to the following naming conventions**:
```bash
/home/myUserName/pipeline/subjects/
|-- Sub1/
|   |-- RAWDATA/
|   |   |-- MPRAGE/
|   |   |   |-- Maybe/Some/SubFolders
|   |   |   |   |-- Arbitrary-Image-Names-001.dcm
|   |   |   |   |-- Arbitrary-Image-Names-002.dcm
|   |   |   |   |-- ...
|   |   |-- DTI/ 
|   |   |-- BOLD-EPI/ 
```
Inside the several folders for the different imaging modalities, the number of subfolder doesnt matter.
Note that the pipeline currently only support DICOM data as input

##### Using fMRI data is optional, i.e. if you dont include that data into your RAWDATA-folder, you still get the structural and dwMRI data processed!

## Running the Pipeline

To finally run the pipeline, locate the **TVB_pipeline.py** script using your systems Shell and pass the subjects ID and the absolute path to the folder holding your subjects RAWDATA-folder (see above).
```bash
python /home/myUser/pipeline/TVB_pipeline.tyb --sub-id <SUBJECT-ID> --sub-dir <SUBJECT-DIR>
```

The log-files are stored into a subfolder of your **SUBJECT-DIR** called TVB_pipeline.

## The Results of the Pipeline

Among several intermediate results, like a full **FREESURFER recon_all** dataset, there are also several datasets which are in-house developed. THe generation is described in the aformentioned research article. The following tables can bee seen as a reference linking the explanations in the paper to the file- and variable-names which are generated by the pipeline-code.

##### Diffusion-MRI:
The results are by default stored into **\<SUBJECT-DIR>/tractography/tracks/\<SUBJECT-ID>_SC.mat** (MATLAB/Octave file) and also in JSON format **\<SUBJECT-DIR>/tractography/tracks/\<SUBJECT-ID>_SC.json**
Those files include several matrices representing different metrics:

| Variable-Name | Type of Data | Refered to in the paper as |
| ------------- |:-------------:|:-------------:|
| SC_cap_agg_counts | Region-wise Capacity Matrix using the number of tracts found between different regions | Raw Counts |
|SC_cap_agg_bwflav1| Region-wise Capacity Matrix using the number of distinct connections found on single-voxel level | Distinct Connections |
|SC_cap_agg_bwflav1_norm|Same data as above but normalized to the range between 0 and 1|||
|SC_cap_agg_bwflav2| Region-wise Capacity Matrix using the number of distinct connections found on single-voxel level. Each strength entry is weighted by the total number of connections leaving the corresponding brain area | Weighted Distinct Connections |
|SC_cap_agg_bwflav2_norm|Same data as above but normalized to the range between 0 and 1|||
|SC_dist_\<mean/mode/median>_agg|The mean/mode/median distance between all distinct tracks connecting the individual brain regions|SC Distances|
|SC_dist_var_agg|The variance of the distance between all distinct tracks connecting the individual brain regions||

##### Functional-MRI:
By default, resulting data will be stored into ** \<SUBEJCT-DIR>/bold/**. Results feature a run of FSLs feat pipeline and also regionswise timeseries stored into the file ** \<SUBJECT-ID>_fMRI.mat **. As for the SC-file descbried above, this MATLAB/Octave file stores various things:

| Variable-Name | Type of Data | 
| ------------- |:-------------:|
|ROI_ID_table|Various numbers from FREESURFERs **mri_segstat**. The headers have been removed. They can be found in the following file: **\<SUBEJCT-DIR>/bold/segstat_summary.txt**|

----------

## TODO-List:
+ ~~Write a new Documentation!~~
+ ~~Implement fMRI processing based on the code here: https://github.com/BrainModes/TVB-empirical-data-pipeline/blob/NSG/fmriFC.sh~~
+ Implement Tractography Thresholding into the MRTrix module. Possibly trying to include the method described in [Morris et al. (2008)](http://www.sciencedirect.com/science/article/pii/S1053811908007301). Alternatively one could also dig up the old hard-threshold code since the short-range tracking-flares are rendered meaningles anyway by our aggregation method!
+ Make the file-sorting of the user-data more sophisticated. This means that the pipeline should be able to somehow recognize which kinds of data-sets (e.g. fMRI, T1, dwMRI) is included in the user data and then route the particular folder-paths onto the corresponding processing-nodes inside the pipeline. This might be achieved through using nipype's [SelectFiles interface](http://nipy.org/nipype/users/select_files.html)
+ Include some example workflows for different cluster scenarios, realized through e.g. controll-scripts written in BASH
+ Re-Implement Multishell-Tracking using FSLs bedpostx as in https://github.com/BrainModes/TVB-empirical-data-pipeline/tree/multiShell
+ Implement the formatting of the results into a TVB-ZIP-File as in https://github.com/BrainModes/TVB-empirical-data-pipeline/blob/NSG/matlab_scripts/connectivity2TVBFS.m
+ Check if Non-DICOM data works as input
+ Add a Doc-Section about the resulting data
+ Support multiple runs of fMRI (e.g. bold1; bold2; bold3; ...)
