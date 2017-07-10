source('../common/variant_calling.R')

#After Bcbio variant calling, run this script under /final
dir.create('./vcf_201604/',showWarnings = F)
dir.create('./vcf_201604_filter/',showWarnings = F)
system('find ./ -name *.vcf.gz | xargs cp -t ./vcf_201604/')
system('gzip -d ./vcf_201604/*')
vcf_filter('./vcf_201604/','./vcf_201604_filter/')

#Get vep files, require ensembl-tools-release-85 with variant effect predictor installed under ~/ensembl-tools-release-85/, change if using under different setting 
vcf2vep('./vcf_201604_filter/','./vep_201604/')

#Run this line after vep conversion: get final results
mut_table <- vep_parser('./vep_201604/')
mut_table$SAMPLE <- gsub('_R-freebayes','',mut_table$SAMPLE)
write.table(mut_table,'201604_variantion.tsv',sep = '\t')

#for reference data(mcf7 fastq files from Encode), move the .vcf file into ~/mcf7_encode and run:
dir.create('./vcf_mcf7_f/',showWarnings = F)
dir.create('./vcf_mcf7_f/',showWarnings = F)
vcf_filter('./mcf7_encode/','~/vcf_mcf7_f/')
vcf2vep('./vcf_mcf7_f/','./vep_mcf7/')
mut_table <- vep_parser('./vep_mcf7/')
write.table(mut_table,'mcf7_ref_variantion.tsv',sep = '\t')
