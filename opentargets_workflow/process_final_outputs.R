library(tidyverse)
library(biomaRt)
args = commandArgs(trailingOnly=TRUE)
resfile=args[2]
qtlthres=args[3]
snplist=args[4]
indexfile=args[5]

data <- read_tsv(resfile)
data$SNPid <- paste("chr",data$chr_id,":",data$position,"_",data$ref_allele,"_", data$alt_allele, sep="")

varindex <- read_tsv(indexfile)
varindex$SNPid <- paste("chr",varindex$chr_id,":",varindex$position,"_",varindex$ref_allele,"_", varindex$alt_allele, sep="")

snplist <- read.table(snplist, head=F)
data <- data[data$SNPid %in% snplist$V1,]
varindex <- varindex[varindex$SNPid %in% snplist$V1,]

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','band'), mart = ensembl)

varindex$alt_af_NFE <- sapply(strsplit(varindex$af,","), `[`, 6)
varindex$alt_af_NFE <- sapply(strsplit(varindex$alt_af_NFE," "), `[`, 3)
varindex$alt_af_NFE <- as.numeric(as.character(varindex$alt_af_NFE), digits=4)
varindex <- subset(varindex, select=c("SNPid","rs_id","alt_af_NFE","most_severe_consequence"))

data <- data %>%
    left_join(genes, by=c("gene_id" = "ensembl_gene_id")) %>%
    left_join(varindex, by="SNPid")

data <- data[order(data$SNPid, -(data$overall_score), -(data$qtl_score)),]

write_tsv(data, file="all_results.tsv")

##

data2 <- data %>%
    dplyr::select(band, chr_id, position, SNPid, rs_id, alt_af_NFE, most_severe_consequence, source_id, hgnc_symbol, overall_score)

snplist <- unique(data$SNPid)

for(i in 1) {

    holding <- data2[data2$SNPid == snplist[i],]
    holding <- holding[holding$overall_score == max(holding$overall_score),]

    data_top <- holding %>%
    group_by(hgnc_symbol) %>%
    slice_head()

}

for(i in 2:length(snplist)) {

    holding <- data2[data2$SNPid == snplist[i],]
    holding <- holding[holding$overall_score == max(holding$overall_score),]

    holding_top <- holding %>%
    group_by(hgnc_symbol) %>%
    slice_head()

    data_top <- rbind(data_top, holding_top)
}

data_top <- data_top %>%
    dplyr::select(band, chr_id, position, SNPid, rs_id, alt_af_NFE, hgnc_symbol, overall_score, most_severe_consequence) %>%
    rename(V2G_gene = hgnc_symbol) %>%
    rename(V2G_score = overall_score) %>%
    mutate(band = paste(chr_id, band, sep="")) %>%
    rename(Locus = band)

##

data_vep <- data[data$source_id == 'vep',]

data_vep <- data_vep[order(-(data_vep$fpred_max_score), -(data_vep$overall_score)),]

data_vep <- data_vep %>%
    group_by(SNPid) %>%
    slice_head() %>%
    dplyr::select(SNPid, fpred_max_label, hgnc_symbol) %>%
    rename(VEP_gene = hgnc_symbol) %>%
    rename(VEP_most_severe = fpred_max_label)

##

full_result <- data_top %>%
    left_join(data_vep, by="SNPid") %>%
    rename(GEL_variant_ID = SNPid) %>%
    arrange(chr_id, position) %>%
    dplyr::select(-chr_id, -position)

write_tsv(full_result, file="test_out.tsv", quote=F)

v2g_list <- full_result %>%
    dplyr::select(V2G_gene)

write_tsv(v2g_list, file="v2g_genelist.tsv", col_names=F, quote=F)
