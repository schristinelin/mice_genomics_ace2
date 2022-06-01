library(GEOquery)
library(gplots)
library(data.table)
library(flowCore)
library(dplyr)

fpath <- dirname(rstudioapi::getSourceEditorContext()$path)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205014
# Cytokines and lipid mediators of inflammation in lungs of SARS-CoV-2 infected mice

geo_id <- 'GSE205014'

fps <- getGEOSuppFiles("GSE205014", fetch_files = T)
file_name <- rownames(fps)

if (length(file_name) > 1) {
  print('more than 1 files present, please specify file name')
} else {
  ace2_df = fread(file_name)
}

ace2_num <- ace2_df %>% select(-c('id', 'ensembl_gene', 'entrez_id', 'transcript_type'))
ace2_num[,2:ncol(ace2_num)] <- lapply(ace2_num[,2:ncol(ace2_num)],as.numeric)
subset_df <- aggregate(.~symbol, ace2_num, mean)
rownames(subset_df) <- subset_df$symbol
subset_df <- subset_df %>% select(-('symbol'))


# data is normalized counts
# GR8 1-4 -> Control groups (non-infected) Rep 1-4
# PROT2GR1-1 TO PROT2GR1-5 -> Mice lung expression, day 5 of infection, Wuhan variant Rep 1-5
# PROT2GR3-1 to PROT2GR3-5 -> Mice lung expression, day 7 of infection, Wuhan variant Rep 1-5
# PROT6GR2-6 to PROT6GR2-10 -> Mice lung expression, day 5 of infection, Delta variant Rep 1-5
# PROT6GR2-11 to PROT6GR2-15 -> Mice lung expression, day 7 of infection, Delta variant Rep 1-5

# initiate pdf to save plots
pdf(file = paste0(fpath, "/mice_genomics_heatmap.pdf"))

heatmap.2(as.matrix(subset_df), main = 'Heat map of expressions')

# close pdf
dev.off()