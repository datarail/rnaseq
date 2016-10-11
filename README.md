# Laboratory of Systems Pharmacology - RNAseq pipeline

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

5. Make a .csv file describing your samples. Use `sample_description.csv` in this GitHub repository as a template.
	
6. Run the following line of code:
    ```
    bcbio_nextgen.py -w template yaml_example.yaml sample_description.csv *.fq 
    ```
    - NOTE: Use `*.fq` or `*.fastq` depending on how the files are named
    - NOTE: fastq filenames CANNOT end in _1 or _anynumber.fastq or it will cause FATAL ERROR

7. Running this setup of the data for bcbio creates project directory and the following subdirectories: `work`, `final`, `config`

8. Move into the `work` directory. Create the following submission file using `emacs submit_bcbio.lsf` or `vim submit_bcbio.lsf`:
    ```
    #!/bin/sh
    
    #BSUB -q priority
    #BSUB -J bcbio_mov10
    #BSUB -n 1
    #BSUB -W 100:0
    #BSUB -R "rusage[mem=10000]"
    #BSUB -e mov10_project.err
    #BSUB -o mov10_project.out
    
    bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380
    ```

9. Submit the job to Orchestra using the following command on the newly-created submission file: `bsub < submit_bcbio.lsf`
    - The job will sometimes time out or have memory issue. In this case, one can just re-submit the job with different settings. The job will look at where the previous run left off and start from there.
    - Use the following alternate submission settings: `bcbio_nextgen.py ../config/mov10_project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380`

10. Note: The instructions assume hg19 - genome build 37. Additional reading:
    - http://bcbio-nextgen.readthedocs.io/en/latest/index.html
    - https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates

## Running differential expression on the count matrices

Instructions on how to run our edgeR / DESeq2 pipeline here
