# Laboratory of Systems Pharmacology - RNAseq pipeline

This is a step-by-step guide on running the `bcbio` pipeline on O2. Old orchestra instructions are still available [here](https://github.com/datarail/rnaseq/tree/master/Orchestra).

## Prerequisites

* If this is your first time using a linux command shell, make sure you type the commands in exactly as is. Things like spaces and quotes matter. The easiest way to do this is to copy and paste the commands from this README into your command shell.

* Login to O2 using the instructions available at https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic

* Familiarize yourself with [nano](https://www.nano-editor.org/dist/v2.2/nano.html), a text editor in linux. Specifically, take a look at keyboard shortcuts like `CTRL+O` to save, `CTRL+X` to exit

* Use `nano ~/.bashrc` to open your `.bashrc` file and add the following lines needed for running the `bcbio` pipeline:
    ```
    # Environment variables for running bcbio
    export PATH=/n/app/bcbio/dev/anaconda/bin/:/n/app/bcbio/tools/bin:$PATH
    alias 'soar=srun -p interactive -t 0-12:00 -n 2 --pty /bin/bash'
    ```
    Save and exit (`CTRL+O` and `CTRL+X` in `nano`). **Note that this only needs to be done once.**
    The `export` line tells your environment where `bcbio` is located, while the second line defines a new command `soar`, which you will use to request an interactive shell after logging into O2.
    
* Log out of O2 and log back in. This will ensure that the new lines in `.bashrc` take effect. (Alternatively, you can type `source ~/.bashrc` without logging out.)

* When you first log into O2, you are given a login shell. While you will be able to submit jobs (like `bcbio` below) from the login shell, you will not be able to run differential expression scripts. To do the latter, you need to request an interactive shell using the new `soar` command you just defined. Simply type `soar` and you should get an output that looks something like this:

    ```
    [as773@login05 ~]$ soar
    srun: job 6094734 queued and waiting for resources
    srun: job 6094734 has been allocated resources
    [as773@compute-a-16-70 ~]$
    ```
    where `as773` is replaced by your eCommons ID.
    
* Finally, if you intend on running the R scripts in this repository (e.g., `merge-sf.R` and `edge.R`), please install the required packages by executing the following commands within R:
    ```
    source("https://bioconductor.org/biocLite.R")
    biocLite( c("tidyverse","stringr","magrittr","edgeR") )
    ```

## Converting FASTQ files to count matrices

* Change into your directory on the `groups` partition: `cd /n/groups/<groupname>/<yourname>/`
    (For example, this could be `/n/groups/lsp/Sokolov/`)
    
* Clone the contents of this GitHub repository by executing `git clone https://github.com/datarail/rnaseq.git`. This will download all the necessary configuration files and `R` scripts into a newly-created `rnaseq` directory. Feel free to rename this to your project name using `mv rnaseq <... project name ...>`.

* Change into the newly created directory: `cd rnaseq` (or whatever you renamed it to in the previous step.)

* Move your fastq files to this directory. (Using `mv` command.)

* Examine `yaml_O2.yaml` and ensure that it has the correct specification. Most of the time, the default content should suffice. Change it only if instructed by the sequencing core.

* Edit `sample_description.csv` (using `nano sample_description.csv`) to match your fastq files. Columns `samplename` and `description` are required. The first should match the filenames of your fastq files. The second column is often set to be the plate position, but can be any abbreviated identifier of your sample. All other columns are optional. The `sample_description.csv` file in this GitHub repository includes an `Index` column as an example. Work with the sequencing core that produced your fastq files to determine if you need additional columns.

    - **NOTE: If you intend to use `hbc/bcbioRNASeq` pipeline to perform differential expression analysis (see below), it is HIGHLY recommended to include all relevant metadata for your samples into `sample_description.csv`.**

* Run the following line of code:
    ```
    bcbio_nextgen.py -w template yaml_O2.yaml sample_description.csv *.fastq
    ```
    - NOTE: Use `*.fq` or `*.fastq` depending on how the files are named.
    - NOTE: bcbio can accept gzipped and bzipped files. Replace `*.fastq` with `*.fastq.gz` or `*.fastq.bz2` accordingly.
    - NOTE: fastq filenames CANNOT end in _1 or _anynumber.fastq or it will cause FATAL ERROR

* Running this setup creates a project directory and the following subdirectories inside of that: `work`, `config`.
    - NOTE: The project directory will be named `sample_description`. You can control that by substituting a different name in all filenames and commands above.
    
* Descend into the `work` subdirectory of your project (by default, this is done using `cd sample_description/work`). 
    - NOTE: If you changed the name of your project from `sample_description` to something else, edit `submit_bcbio.sh` and specify the correct name in place of `../config/sample_description.yaml`.The path has to point to an actual file in the `config` subdirectory.

* From the `work` directory, submit the job to O2 using the following command:
    ```
    sbatch ../../submit_bcbio.sh
    ```
* Bcbio will give you a .count file as count table generated from featureCounts. This could be used as direct input to the `edge.R` script (see Running Differential Expression section below). The counts file will be located in `../final/2017-10-19_sample_description/combined.counts`, but depending on the current date and the contents of your `sample_description.csv`, the actual path may vary.

* You might also interested in using alignment-independent quantification results from salmon for differential analysis. To get a count table from salmon you need `quant.sf` files that are scattered across individual-sample directories. Running the following on the command line will pull all the `.sf` files together and output a counts table and a TPM (transcripts-per-million) table:
    ```
    Rscript merge-sf.R
    ```
    - NOTE: `merge-sf.R` is located in your main project directory (where your FASTQ files are located)
    - NOTE: The resulting counts and TPM tables will be in the `final` subdirectory

## Running Differential Expression analysis

This GitHub repository includes a simple command-line interface for EdgeR. Simply type `Rscript edge.R` to get information on how to run it. The instructions include an example usage, which uses a counts table and metadata file available in the `example/` subdirectory of this repository. The third argument to the script specifies which column in the meta file should be used for defining differential expression contrasts. Finally, the fourth argument denotes which value in that column should be used as control. For example, to compare samples across timepoints, using 0h as control, one can do
```
Rscript edge.R example/test.count example/meta.tsv Timepoint 0h
```
Likewise, to perform a comparison along the Treatment column, one would do
```
Rscript edge.R example/test.count example/meta.tsv Treatment Control
```
If the specified column contains more than two categories (such as Timepoint above), edgeR will compare each category to the control. These results will reside in `logFC` columns in the output file. Additionally, edgeR will perform an F-test-like comparison to determine if any genes differ across ALL categories. These results are then captured by the `F`, `PValue` and `FDR` columns.

The metadata file provides some flexibility in defining the comparison. For example, to perform a time-series differential expression analysis for one cell line at a time, one can simply break up their metadata file into multiple files, each limited to the cell line of interest.

### Advanced differential expression analysis

A more sophisticated differential expression analysis is available as part of the bcbio pipeline. If you remembered to include metadata in your `sample_description.csv` file prior to running the aligner, the usage of the follow-up RNAseq analysis is straightforward. For more information, see https://github.com/hbc/bcbioRNASeq

(NOTE: If your `sample_description.csv` file did not include metadata, running bcbioRNAseq is still possible, but requires you to specify the metadata information explicitly when calling `loadRNASeq()`.)
