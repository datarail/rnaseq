# Laboratory of Systems Pharmacology - RNAseq pipeline

## Converting FASTQ files to count matrices

0. If this is your first time using a linux command shell, make sure you type the commands in exactly as is. Things like spaces and quotes matter. The easiest way to do this is to copy and paste the commands from this README into your command shell.

1. Login to Orchestra using the instructions available at https://wiki.med.harvard.edu/Orchestra/NewUserGuide

2. Familiarize yourself with `nano`, a text editor in linux. Specifically, take a look at keyboard shortcuts like `CTRL+O` to save, `CTRL+X` to exit

3. Use `nano ~/.bashrc` to open your `.bashrc` file and add the following lines to have the path to bcbio in your orchestra setup:
    ```
    # Environment variables for running bcbio
    export PATH=/opt/bcbio/centos/bin:$PATH
    unset PYTHONHOME
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    module load dev/java/jdk1.7
    ```
    Save and exit (`CTRL+O` and `CTRL+X` in `nano`).    
    This only needs to be done once.

4. Change into the project directory: `cd /groups/<groupname>/<yourname>/<projectname>/`
    (For example, this could be `/groups/springer/sarah/Isog/`)
 
5. Move your fastq files to this directory. (Using `mv` command on Orchestra.)

6. Move your yaml file to this directory as well. An example .yaml file is available in this GitHub repository as `yaml_example.yaml`:
    - OPTION 1: Open a new file with a text editor. Copy and paste the content of `yaml_example.yaml` to that file. Save and exit.
    - OPTION 2: Download the file directly from GitHub using the following command: `wget https://raw.githubusercontent.com/sorgerlab/rnaseq/master/yaml_example.yaml`

7. Make a .csv file describing your samples. Use `sample_description.csv` in this GitHub repository as a template. As before, you can either create a new file and copy-and-paste OR using `wget https://raw.githubusercontent.com/sorgerlab/rnaseq/master/sample_description.csv`
	
8. Run the following line of code:
    ```
    bcbio_nextgen.py -w template yaml_example.yaml sample_description.csv *.fq 
    ```
    - NOTE: Use `*.fq` or `*.fastq` depending on how the files are named
    - NOTE: fastq filenames CANNOT end in _1 or _anynumber.fastq or it will cause FATAL ERROR

9. Running this setup creates a project directory and the following subdirectories inside of that: `work`, `final`, `config`
    - NOTE: The project directory will be named `sample_description`. You can control that by substituting a different name in all filenames and commands above.

9. Descend into the `work` subdirectory of your project (by default, this is done using `cd sample_description/work`). Using a text editor (e.g., `nano submit_bcbio.lsf`) create the following submission file:
    ```
    #!/bin/sh
    
    #BSUB -q priority
    #BSUB -J bcbio_project
    #BSUB -n 1
    #BSUB -W 100:0
    #BSUB -R "rusage[mem=10000]"
    #BSUB -e project.err
    #BSUB -o project.out
    
    bcbio_nextgen.py ../config/sample_description.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380
    ```
    - NOTE: If you changed the name of your project from `sample_description` to something else, use that new name when specifying `../config/sample_description.yaml` above. The path has to point to an actual file in the `config` subdirectory.

10. Submit the job to Orchestra using the following command on the newly-created submission file: 
    ```
    bsub < submit_bcbio.lsf
    ```
    - The job will sometimes time out or have memory issue. In this case, one can just re-submit the job with different settings. The job will look at where the previous run left off and start from there.
    - Use the following alternate submission settings: `bcbio_nextgen.py ../config/yaml_example.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380`
    - Bcbio will give you a .count file as count table generated from featureCounts. This could be used as direct input of run_de.R script. However, you might also interested in using  alignment-independent quantification results from salmon for differential analysis. To get count table from salmon you need to first combine the .sf files for each sample and merge the counts of same transcript to gene level by running following code in `R`:
		```
		source('de_rnaseq.R')
		setwd(/path/to/bcbio/final)
		get_sf()
		sf <- read.table('combined.sf',header = T,as.is = T)
		tx2gene <- read.csv('tx2gene.csv',header = F,as.is = T)
		count_table <- sf2tpm(sf,tx2gene,count=T)
		write.table(count_table,'combined.counts',quote=F,row.names=F,sep='\t')

		```

11. Note: The instructions assume hg19 - genome build 37. Additional reading:
    - http://bcbio-nextgen.readthedocs.io/en/latest/index.html
    - https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates

## Running differential expression on the count matrices

To run the next part of the analysis, you will need a few R packages. To install it on Orchestra, start up `R`, by typing the following commands in your shell:

```
module load stats/R/3.3.1
R
```

This gets you into an `R`. You can now install the packages by running the following command in `R`:
```
install.package( c("ggplot2", "reshape", "gplots", "edgeR", "CHBUtils", "pheatmap",
              "DESeq2", "tximport", "DT", "DEGreport", "dplyr", "rmarkdown") )
```
    
Go get a cup of coffee!

Once you obtained the count matrices and finish installing `Rmarkdown`, the basic differential expression pipeline can be run via the following command:

    bcbio-rnaseq summarize path-to-project-summary-yaml -f "~celltype"
    
where `path-to-project-summary-yaml` refers to a file that you should be able to find in your `final` subdirectory (which is on the same level as the `config` and `work` subdirectories covered above). If you didn't change the name of `sample_description.csv`, the summary yaml file will live at `/groups/<groupname>/<yourname>/<projectname>/sample_description/final/<date>_sample_description/project-summary.yaml`.

The `bcbio-rnaseq` command above will run the DESeq2 differential expression pipeline and produce a report in the R markdown format (http://rmarkdown.rstudio.com/). The formula that you provide after the `-f` flag must be in the appropriate format (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) and reference the column names in your `sample_description.csv` file. In the example above, we provide `"~celltype"`, which implies that we are interested in performing a differential expression comparison along the third column of our `sample_description.csv`.

Note that if you omit the formula, the script will perform basic quality control instead:

    bcbio-rnaseq summarize path-to-project-summary-yaml

## run_de.R - Run edgeR differential expression analysis from bcbio count results and sample annotation
Usage: run_de.R [options]

Example:
```
Rscript run_de.R -c path/to/rnaseq.count -a path/to/annotation.txt -o path/to/output

Rscript run_de.R -c example/input/test.count -a example/input/group.tsv -o example/output
```

Options:

	-c CHARACTER, --count=CHARACTER
		Path to .count file from bcbio output, which is a ensumbl ID by sample ID matrix

	-a CHARACTER, --annotation=CHARACTER
 		Path to annotation information file, which is a dataframe with 3 columns: group, condition and control
               		group: contains information which treatment samples will be compared against control cases in each group
               		condition: indicates type of treatment, replicates have same condition
               		control: TRUE for controls and FALSE for treatments
               		order of samples in annotation must be the same as samples in count table

	-o CHARACTER, --output=CHARACTER
  		Path to save differential analysis results

	-p TRUE/FALSE, --pairwise=TRUE/FALSE
  		If the P-values and FDR are given pairwise or as ANOVA-like test for any differences

	-s TRUE/FALSE, --symbol=TRUE/FALSE
  		If gene symbols will be added to the output

	-h, --help
  		Show this help message and exit

The output consists 3 files, each of them is a n by p matrix with n genes and p type of treatments:

logFC.csv: log2 fold change

pval.csv:  P-values

fdr.csv:   FDR
