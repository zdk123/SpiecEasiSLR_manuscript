library(phyloseq)
source('QMP_fixed.R')


qmp <- import_biom('working_dir/pick_otus/otu_table.biom')
sample_data(qmp) <- import_qiime_sample_data('mapping_file.txt')

qmp_ggcopy <- read.delim('picrust_copyno.txt', row.names=1)

qmp_copyadj <- qmp
qmp_copyadj@otu_table@.Data <- qmp@otu_table@.Data / qmp_ggcopy$copy


abscount <- qmp_copyadj@sam_data$Average_cell_count_per_gram_frozen
qmp_copyadj_rare <- qmp_copyadj
qmp_copyadj_rare@otu_table@.Data <- t(rarefy_even_sampling_depth(t(qmp_copyadj@otu_table@.Data),
            t(abscount)))

ranks <- colnames(tax_table(qmp))[2:6]
qmp_taxsum <- lapply(ranks, tax_glom, physeq=qmp, NArm=FALSE)
names(qmp_taxsum) <- ranks

qmp_copyadj_taxsum <- lapply(ranks, tax_glom, physeq=qmp_copyadj, NArm=FALSE)
names(qmp_copyadj_taxsum) <- ranks

qmp_copyadj_rare_taxsum <- lapply(ranks, tax_glom, physeq=qmp_copyadj_rare, NArm=FALSE)
names(qmp_copyadj_rare_taxsum) <- ranks

save("qmp", "qmp_copyadj", "qmp_copyadj_rare", "qmp_copyadj_rare_taxsum",
     "qmp_copyadj_taxsum", "qmp_ggcopy", "qmp_taxsum", file='QMPphy.RData')

## subset healthy samples
phy <- subset_samples(qmp_copyadj_taxsum$Rank6, Health_status=='Healthy')

save(phy, file='simulator/QMPphyfilt.RData')
