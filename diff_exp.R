if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('pheatmap')) install.packages('pheatmap'); library('pheatmap')
if (!require('stringi')) install.packages('stringi'); library('stringi')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('apeglm')) install.packages('apeglm'); library('apeglm')

# reading and parsing count data
init_sample <- read.delim("~/RNAseq_counts/template.counts", header=FALSE)
calico_data_challenge_samplesheet <- read.delim("~/calico_data_challenge_samplesheet.tsv")

fnames <- calico_data_challenge_samplesheet$name
calico_data_challenge_samplesheet$name <- fnames

fnames <- gsub('.{9}$', '', fnames)
sample_cts <- init_sample

for (f in fnames) 
  {
  name <- paste("~/Challenge/RNAseq_counts/", f, ".counts", sep="")
  data_col <- read.delim(name, header=FALSE)
  names(data_col)[names(data_col) == 'V2'] <- gsub('.{13}$', "", f)
  sample_cts <- merge(sample_cts, data_col, by='V1', all=TRUE)
}
sample_cts$V2 <- NULL
rownames(sample_cts) <- sample_cts$V1
sample_cts$V1 <- NULL

# remove first five rows (unaligned, ambiguous etc)
sample_cts <- tail(sample_cts, -5)
metadata <- calico_data_challenge_samplesheet[, c("name", "strain", "age", "batch")]
metadata$name <- gsub('.{13}$', "", metadata$name)
metadata$strain <- str_remove_all(metadata$strain, "âˆ†")
metadata$age <- factor(metadata$age)

metadata <- mutate_if(metadata, is.character, as.factor)
data_deseq <- DESeqDataSetFromMatrix(countData = sample_cts, colData = metadata, design=~batch+strain+age)

# EDA to find out major factors affecting variance

# PCA of samples
rld_mat <- rlog(data_deseq)
pcaData <- plotPCA(rld_mat, intgroup=c("strain", "age"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age, shape=strain)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave("pca.png")

# check effect of time using LRT
dds_lrt_time <- DESeq(data_deseq, test="LRT", reduced = ~ batch+strain)
res20v0 <- lfcShrink(dds_lrt_time, coef="age_20_vs_0", type="apeglm")
res40v0 <- lfcShrink(dds_lrt_time, coef="age_40_vs_0", type="apeglm")

# get list of significant genes (p-value < 0.001)
padj.cutoff <- 0.001

sig_genes_40v0 <- res40v0 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

sig_genes_20v0 <- res20v0 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# create heatmap of top 20 differentially expressed genes with age (40vs0)
sorted_genes <- sig_genes_40v0[order(sig_genes_40v0$padj),]
topgenes <- head(sorted_genes$gene,20 )
mat <- assay(rld_mat)[topgenes,]
mat <- mat - rowMeans(mat)
colnames(mat) <- metadata$age
pheatmap(mat)
ggsave("age_heatmap.png")

# get upregulated (20v0 positive, 40v0 > 20v0)
up_age <- subset(sig_genes_40v0, sig_genes_20v0$log2FoldChange > 0 & sig_genes_40v0$log2FoldChange > sig_genes_20v0$log2FoldChange)
# get downregulated (20v0 negative, 40v0 < 20v0)
down_age <- subset(sig_genes_40v0, sig_genes_20v0$log2FoldChange < 0 & sig_genes_40v0$log2FoldChange < sig_genes_20v0$log2FoldChange)
up_age[, ] <- lapply(up_age[, ], as.character)
write.table(up_age, file = "Upregulated_age.tsv",
            sep = "\t", row.names = F)
write.table(down_age, file="Downregulated_age.tsv",
            sep= "\t", row.names=F)

## Strain analysis with LRT
dds_lrt_strain <- DESeq(data_deseq, test="LRT", reduced = ~ batch+age)
resultsNames(dds_lrt_strain)
strain_data <- c("strain_FOB1_vs_DBY1200", "strain_SIR2_vs_DBY1200", "strain_UBR2_vs_DBY1200")

for (strain in strain_data) {
    res <- lfcShrink(dds_lrt_strain, coef=strain, type="apeglm")
    sig_genes <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < padj.cutoff)
    
    up <- subset(sig_genes, sig_genes$log2FoldChange > 0)
    down <- subset(sig_genes, sig_genes$log2FoldChange < 0)
    
    up_f <- paste("Upregulated_strain_", strsplit(strain, split = "_")[[1]][2], ".tsv", sep="")
    down_f <- paste("Downregulated_strain_", strsplit(strain, split = "_")[[1]][2], ".tsv", sep="")
    
    write.table(up, file = up_f, sep = "\t", row.names = F)
    write.table(down, file = down_f, sep= "\t", row.names=F)
  
}

# do this for all three strains
# sig_genes_FOB1 <- resFOB1 %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble() %>% 
#   filter(padj < padj.cutoff)


# modeling interaction effects
data_new <- DESeqDataSetFromMatrix(countData = sample_cts, colData = metadata, design=~batch+strain+age+strain:age)
dds_lrt_interaction <- DESeq(data_new, test="LRT", reduced = ~ batch+strain+age)
res_inter <- results(dds_lrt_interaction)

# heatmap to find clusters of differential expression profiles
betas <- coef(dds_lrt_interaction)
topGenes <- head(order(res_inter$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

ggsave("strain_age_heatmap.png")

interaction_terms <- c("strainFOB1.age40", "strainSIR2.age40", "strainUBR2.age40")
for (term in interaction_terms) {
  res <- results(dds_lrt_interaction, name=term)
  
  sig_genes <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < padj.cutoff)
  
  up <- subset(sig_genes, sig_genes$log2FoldChange > 0)
  down <- subset(sig_genes, sig_genes$log2FoldChange < 0)
  name <- strsplit(term, split="[.]")
  
  name <- stri_join_list(name, sep="+")
  
  up_f <- paste("Upregulated_strain+age_", name,".tsv", sep="")
  down_f <- paste("Downregulated_strain+age_", name, ".tsv", sep="")
  
  write.table(up, file = up_f, sep = "\t", row.names = F)
  write.table(down, file = down_f, sep= "\t", row.names=F)
  
}
