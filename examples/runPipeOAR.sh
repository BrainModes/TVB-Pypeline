#!/bin/bash

# Use this script to run the pipeline as a Python job in the backhground which will manage the execution of the
# pipeline-graph. For usage, please refer to the main Python-file TVB_pipeline.py!

### Check Input ###
export usage=$(python -B ../workflows/TVB_pipeline.py --help 2>&1)
export subID=none
export subPath=none
export struct=none
export diff=none
export func=none
while [ $# -gt 0 ]
do
    case "$1" in
	-s) subID="$2"; shift;;
	-r) subPath="$2"; shift;;
	-a) struct="$2"; shift;;
	-d) diff="$2"; shift;;
	-f) func="$2"; shift;;
	-*) echo >&2 \
	    ${usage}
	    exit 1;;
	*)  break;;	# terminate while loop
    esac
    shift
done

#Check if the parameters have been set
if [ "$subID" == "none" ]
	then
		echo >&2 \
		"-s SubjectID is missing! +++"
	    exit 1;
elif [ "$subPath" == "none" ]
	then
		echo >&2 \
		"-r /absolute/path/to/subject is missing! +++"
 	    exit 1;
elif [ "$struct" == "none" ]
	then
		echo >&2 \
		"-a /absolute/path/to/structural/rawData is missing! +++"
 	    exit 1;
elif [ "$diff" == "none" ]
	then
		echo >&2 \
		"-d /absolute/path/to/diffusion/rawData is missing! +++"
 	    exit 1;
fi
### Check Input ###


# source ./pipeSetup.sh
if [ "$func" == "none" ]
then
    nohup python -B ../workflows/TVB_pipeline.py -s ${subID} -r ${subPath} -a ${struct} -d ${diff} \
            > ${subPath}/pipeLog.txt &
else
    nohup python -B ../workflows/TVB_pipeline.py -s ${subID} -r ${subPath} -a ${struct} -d ${diff} -f ${func} \
            > ${subPath}/pipeLog.txt &
fi