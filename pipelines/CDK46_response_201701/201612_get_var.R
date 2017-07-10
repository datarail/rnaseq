source('../common/variant_calling.R')

#After Bcbio variant calling, run this script under /final
dir.create('./vcf_201612/',showWarnings = F)
dir.create('./vcf_201612_filter/',showWarnings = F)
system('find ./ -name *.vcf.gz | xargs cp -t ./vcf_201612/')
system('gzip -d ./vcf_201612/*')
vcf_filter('./vcf_201612/','./vcf_201612_filter/')

#Get vep files, require ensembl-tools-release-85 with variant effect predictor installed under ~/ensembl-tools-release-85/, change if using under different setting 
vcf2vep('./vcf_201612_filter/','./vep_201612/')

#Run this line after vep conversion: get final results
mut_table <- vep_parser('./vep_201612/')
mut_table[,1] <- gsub('.+S(\\d+)_.+','\\1',mut_table[,1])
write.table(mut_table,'201612_variantion.tsv',sep = '\t')
