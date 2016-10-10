# LSP RNAseq pipeline

## Converting FASTQ files to count matrices

1. Login to Orchestra using the instructions available at https://wiki.med.harvard.edu/Orchestra/NewUserGuide

2. Use `vim ~/.bashrc` or `emacs ~/.bashrc` to open your `.bashrc` file and add the following lines to have the path to bcbio in your orchestra setup:
    ```
    # Environment variables for running bcbio
    export PATH=/opt/bcbio/centos/bin:$PATH
    unset PYTHONHOME
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    module load dev/java/jdk1.7
    ```
    Save and exit.
3. Move your fastq files to the following location: `/groups/<groupname>/<yourname>/<projectname>/`
    (For example, this could be `/groups/springer/sarah/Isog/`)

4. Move your yaml file to this directory as well. An example .yaml file is available in this GitHub repository as `yaml_example.yaml`.

5. Make a csv file describing your samples
	example sample_description.csv

	
#Run this line of code
##NOTE - (*fq or *fastq depending on how files named)

bcbio_nextgen.py -w template yaml_singleEnd.yaml sample_description.csv *.fq 

##NOTE2 - FASTQ FILES CAN NOT be named _1 or _anynumber.fastq or will cause FATAL ERROR


#Running this sets up the data for bcbio creates project directory and some sub directories.
work, final, config directories are created


#Move into the "work" directory
#Run the lsf from here 


bsub < submit_bcbio.lsf


vim submit_bcbio.lsf
##This is what submission looks like - if need to edit see notes below

#!/bin/bash

#BSUB -q priority
#BSUB -J bcbio_mov10
#BSUB -n 1
#BSUB -W 100:0
#BSUB -R "rusage[mem=10000]"
#BSUB -e mov10_project.err
#BSUB -o mov10_project.out

bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380


## sometimes times out or memory issue - can just re-submit  with different settings - will look at where left off and start from there.
#alternate submission settings:


bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380

##note hg19 - genome build 37


http://bcbio-nextgen.readthedocs.io/en/latest/index.html
https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates
