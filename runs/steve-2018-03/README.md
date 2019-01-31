1. Move FASTQ files to their own directory, just to organize things better. They are in `fastq`.
2. Make a metadata directory and download the newest transcriptome and gene annotation:
```bash
mkdir -p metadata
cd metadata
# download latest cDNA build
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
# download latest gene annotation
wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz
gunzip Homo_sapiens.GRCh38.91.gtf.gz
```
3. Point to the transcriptome in the `transcriptome_fasta` field in the YAML file.
4. Point to the GTF in the `transcriptome_gtf` field in the YANL file.
5. Run the bcbio-nextgen templating system:
```bash
bcbio_nextgen.py -w template SCRB-seq.yaml NeuroDrugs.csv fastq
```

6. Copy the submission script to NeuroDrugs.
```bash
cp submit-bcbio.sh NeuroDrugs/work
```
7. Kick off run:
```bash
cd NeuroDrugs/work
sbatch submit-bcbio.sh
```
