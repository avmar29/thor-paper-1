## Load packages
library(tidyverse)
library(patchwork)
library(tximport)
library(DESeq2)
library(edgeR)
library(limma)
library(extrafont)
library(RColorBrewer)
library(khroma)
library(viridis)
library(vsn)
library(pcaExplorer)
library(ggrepel)
library(pheatmap)
library(ggrepel)
library(ragg)
library(xlsx)

## Set colors
# universal palettes
kpal_bright <- colour("bright")(7)
# set palettes specifically for our data
mypal_prot <- colour("bright")(7)[c(1,4,2)]
mypal_treat <- colour("discreterainbow")(23)
mypal_br_short <- brewer.pal(9, "Set1")
mypal_br_long <- c(rbind(brewer.pal(9, "Set3")[-2], brewer.pal(8, "Dark2")), brewer.pal(12, "Set3")[10:12])
mypal_rainbow_fun <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
## Set thresholds for DEA
p_th <- 0.05
lfc_th <- 0.5
## Load gene sets
thor_list <- read.delim("THOR_Prioritised_Consensus_Summary_Table.txt")
thor_kd_targets <- read.delim("THOR_Knockdown_Genes_R1-2.tsv")


#### SEQUENCING BATCH 1 ####

## Set working dir
setwd("thor_kd_r1")
## Load samplesheet
# "samplesheet_kd_r1.csv", a samplesheet file prepared for 'nf-core/rnaseq' pipeline
sample_info <- read_csv("samplesheet_kd_r1.csv") %>% 
  mutate(sample = as.character(sample))
coldata <- data.frame(files = list.files(path = "star_salmon/", pattern = "quant.sf", recursive = TRUE)) %>% 
  mutate(names = str_replace_all(files, "/.*", ""), 
         files = paste0("star_salmon/", files)) %>% 
  left_join(sample_info, by = c("names" = "sample")) %>% 
  select(!starts_with("fastq")) %>% 
  mutate(sample_num = as.numeric(names), 
         sample_id = paste0("X", names),
         group = paste(Protocol, Treatment, sep = "_"), .after = names) %>% 
  arrange(sample_num)
rownames(coldata) <- coldata$names
file.exists(coldata$files)
table(coldata$Treatment, coldata$Protocol)
## Remove outlier - X32
coldata <- coldata %>% 
  filter(names != "32")
## Import Salmon output using tximport
## Load Tx to Gene table
tx2gene <- read.delim("star_salmon/salmon_tx2gene.tsv", sep = "\t", header = FALSE)
colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
gene_info <- unique(tx2gene[, -1])
rownames(gene_info) <- gene_info$gene_id
## Load TPM table and count stats
ds_tpm <- read.delim("star_salmon/salmon.merged.gene_tpm.tsv", 
                     sep = "\t", header = TRUE, row.names = 1)
ds_tpm <- ds_tpm[, coldata$sample_id]
ds_tpm_mean <- t(apply(ds_tpm, 1, function(x) tapply(x, coldata$Treatment, mean)))
## Collect files info
files <- coldata$files
names(files) <- coldata$names
## Load salmon data and make DESeq object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# Combinatorial LM
dds <- DESeqDataSetFromTximport(txi, coldata, ~ 0 + Treatment + Protocol)
## Pre-filtering with DESeq2
keep <- rowSums(counts(dds)) >= 10; table(keep)
## Add gene names
rowData(dds)$geneSymbol <- gene_info[rownames(dds), "gene_name"]
table(duplicated(rowData(dds)$geneSymbol))
## Perform DESeq analysis 
dds <- DESeq(dds)
resultsNames(dds)
## Data transformations and plots
library(vsn)
dds_vsd <- vst(dds)
meanSdPlot(assay(dds_vsd))
plotPCA(dds_vsd, intgroup = "group")
## Manual PCA
rv <- rowVars(assay(dds_vsd))
sel <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
dds_pca <- prcomp(t(assay(dds_vsd)[sel, ]))
percentVar <- dds_pca$sdev^2 / sum(dds_pca$sdev^2)
scree_plot <- data.frame(PC = 1:length(percentVar),
                         Variance = percentVar * 100)
scree_plot[1:15, ] %>%  
  ggplot(aes(x = PC, y = Variance)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(x = "Principle component", y = "% of variance explained")
library(ggrepel)
pca_data <- DESeq2::plotPCA(dds_vsd, ntop = 500,
                            intgroup = c("Protocol", "Treatment", "Experiment_ID"), 
                            returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# color by protocol
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = Treatment)) +
  #geom_point(aes(color = Type, fill = Protocol), shape = 21, size = 4, alpha = 0.85) + 
  geom_point(aes(shape = Type, fill = Protocol), size = 3, alpha = 0.8) + 
  scale_shape_manual(values = c(21, 24)) + 
  #scale_color_manual(values = mypal_br_short[2:1]) + 
  scale_fill_manual(values = as.character(mypal_prot)) + 
  geom_text_repel(size = 2.5, color = "grey50") + 
  labs(title = paste("Top 500 genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
# color by treatment
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = Treatment)) +
  geom_point(aes(shape = Type, fill = Treatment), size = 3, alpha = 0.85) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_manual(values = as.character(mypal_br_long)) + 
  geom_text_repel(size = 2.5, color = "grey50") + 
  labs(title = paste("Top 500 genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
# color by THOR experiment (batch)
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = Treatment)) +
  geom_point(aes(shape = Type, fill = THOR_Exp), size = 3, alpha = 0.85) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_manual(values = as.character(mypal_br_long)) + 
  geom_text_repel(size = 2.5, color = "grey50") + 
  labs(title = paste("Top 500 genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)"),
       fill = "THOR Exp. ID") + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
## Exclude sample #25 (THOR 6 (08.06.21) - Cholesterol CTR)
r1_dds <- dds[, colnames(dds) != "25"]
r1_vsd <- vst(r1_dds)
r1_counts <- counts(r1_dds)
r1_coldata <- coldata[colnames(r1_dds), ] %>% 
  select(-group) %>% 
  mutate(Treatment = str_replace_all(Treatment, "CTR", "Control")) %>% 
  mutate(Experiment_Short = gsub("HOR ", "", Experiment_ID), .after = Experiment_Date) %>% 
  mutate(Group = paste(Experiment_Short, Treatment, sep = "_"), .after = Treatment)
## Additional gene filtering
# exclude pseudoautosome genes
grep("PAR", rownames(r1_counts), value = T)
r1_counts <- r1_counts[- grep("PAR", rownames(r1_counts)), ]
r1_tpm <- as.matrix(ds_tpm[rownames(r1_counts), paste0("X", colnames(r1_dds))])
r1_tpm_mean <- t(apply(r1_tpm, 1, function(x) tapply(x, r1_coldata$Treatment, mean)))
# investigate duplicated genes
tmp_anno <- gene_info[rownames(r1_counts), ]
table(duplicated(tmp_anno$gene_name))
x <- tmp_anno[duplicated(tmp_anno$gene_name), ]
tmp_dup <- tmp_anno[tmp_anno$gene_name %in% x$gene_name, ] %>% arrange(gene_name)
tmp_dup <- cbind(tmp_dup, r1_tpm_mean[tmp_dup$gene_id, ])
## Save data
save(r1_counts, r1_coldata, r1_dds, r1_tpm, file = "r1_data.rda")
### Simple comparisons within batch
source("deg_functions.R")
r1_results = r1_object <- list()
grcol <- "Treatment"
ds_names <- sort(unique(r1_coldata$group[r1_coldata$Treatment != "CTR"]))
for (i in 1:length(ds_names)) {
  sel_exp <- unique(r1_coldata$Experiment_ID[r1_coldata$group == ds_names[i]])
  sel_gene <- gsub(".+_", "", ds_names[i])
  pheno_table <- r1_coldata %>% 
    filter(Experiment_ID == sel_exp & Treatment %in% c("CTR", sel_gene))
  pheno_table$Treatment <- factor(pheno_table$Treatment, levels = c("CTR", sel_gene))
  count_table <- na.omit(r1_counts[, pheno_table$names])
  res <- list()
  res$deseq <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "deseq2", keep_object = TRUE)
  res$edger <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "edger", keep_object = TRUE)
  res$limma <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "limma", keep_object = TRUE)
  res <- purrr::transpose(res)
  r1_object[[i]] <- res$object
  res_out <- map(res[names(res) != "object"], function(x) data.frame(
    GeneID = rownames(x[[1]]), 
    GeneName = gene_info[rownames(x[[1]]), "gene_name"], 
    do.call("cbind", x)))[[1]]
  # add TPM
  tpm_table <- ds_tpm[res_out$GeneID, pheno_table$sample_id]
  tpm_mean <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, mean)))
  tpm_sd <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, sd)))
  r1_results[[i]] <- data.frame(res_out[, 1:2],
                                TPM.Mean = tpm_mean, TPM.SD = tpm_sd,
                                res_out[, 7:ncol(res_out)])
}; names(r1_results) = names(r1_object) <- ds_names
rm(sel_exp, sel_gene, pheno_table, count_table, res, res_out, tpm_table, tpm_mean, tpm_sd)
save(r1_results, r1_object, file = "r1_3methods_degs.RData")
## Save results to Excel files
options(java.parameters = "-Xmx4096m")
library(xlsx)
kd_names <- sort(unique(r1_coldata$Treatment[r1_coldata$Treatment != "CTR"]))
dir.create("dea_results")
for(i in 1:length(kd_names)) {
  # save txt
  write.table(r1_results[[paste0("Baseline_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Baseline.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(r1_results[[paste0("Cholesterol_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Cholesterol.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(r1_results[[paste0("Stretch_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Stretch.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  # save xlsx
  write.xlsx(r1_results[[paste0("Baseline_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Baseline", append = TRUE)
  write.xlsx(r1_results[[paste0("Cholesterol_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Cholesterol", append = TRUE)
  write.xlsx(r1_results[[paste0("Stretch_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Stretch", append = TRUE)
}
## Get significant results
fdr_th <- 0.05
lfc_th <- 0.5
r1_signif <- list()
r1_signif$deseq <- r1_results %>% 
  map(~ filter(.x, deseq.padj < fdr_th, abs(deseq.log2FoldChange) >= lfc_th))
#lapply(r1_signif$deseq, nrow)
r1_signif$edger <- r1_results %>% 
  map(~ filter(.x, edger.FDR < fdr_th, abs(edger.logFC) >= lfc_th))
r1_signif$limma <- r1_results %>% 
  map(~ filter(.x, limma.adj.P.Val < fdr_th, abs(limma.logFC) >= lfc_th))
## Count number of DEGs
r1_deg_stat <- data.frame(group = names(r1_results),
                          deseq_all = map_dbl(r1_signif$deseq, nrow),
                          deseq_up = map_dbl(r1_signif$deseq, ~ sum(.x$deseq.log2FoldChange > 0)),
                          deseq_down = map_dbl(r1_signif$deseq, ~ sum(.x$deseq.log2FoldChange < 0)),
                          edger_all = map_dbl(r1_signif$edger, nrow),
                          edger_up = map_dbl(r1_signif$edger, ~ sum(.x$edger.logFC > 0)),
                          edger_down = map_dbl(r1_signif$edger, ~ sum(.x$edger.logFC < 0)),
                          limma_all = map_dbl(r1_signif$limma, nrow),
                          limma_up = map_dbl(r1_signif$limma, ~ sum(.x$limma.logFC > 0)),
                          limma_down = map_dbl(r1_signif$limma, ~ sum(.x$limma.logFC < 0)))
r1_deg_stat %>% 
  separate(group, into = c("Protocol", "Treatment"), sep = "_") %>% 
  write_delim(file = "res_deg_number_fdr0.05_lfc0.5.tsv", delim = "\t")
## Automatic GSEA with WebGestalt
library(WebGestaltR)
View(listGeneSet("hsapiens"))
mymethod <- "deseq"
mycolumns <- c("GeneName", "deseq.log2FoldChange")
mypath <- paste0(mymethod, "_gsea_webgestalt")
dir.create(mypath)
x <- "Baseline_MFGE8"
for (x in names(r1_signif[[mymethod]])[-c(1:41)]) {
  print(x)
  xpath <- file.path(mypath, x)
  dir.create(xpath)
  xdat <- r1_signif[[mymethod]][[x]][, mycolumns]
  xdat <- xdat %>% arrange(GeneName)
  xfile <- paste0(xpath, "/", x, ".rnk")
  write_delim(xdat, file = xfile, delim = "\t", col_names = FALSE)
  # Gene Ontology - BP
  WebGestaltR(outputDirectory = xpath, projectName = "GO-BP",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Biological_Process_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - MF
  WebGestaltR(outputDirectory = xpath, projectName = "GO-MF",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Molecular_Function_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - CC
  WebGestaltR(outputDirectory = xpath, projectName = "GO-CC",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Cellular_Component_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # KEGG Pathway
  WebGestaltR(outputDirectory = xpath, projectName = "KEGG",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_KEGG",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Reactome
  WebGestaltR(outputDirectory = xpath, projectName = "Reactome",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Reactome",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # WikiPathway
  WebGestaltR(outputDirectory = xpath, projectName = "WikiPathway",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Wikipathway",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
}; rm(x, xpath, xdat, xfile)


#### SEQUENCING BATCH 2 ####

## Set working dir
setwd("thor_kd_r2")
## Load samplesheet
# "samplesheet_kd_r2.csv", a samplesheet file prepared for 'nf-core/rnaseq' pipeline
## Load samplesheet
sample_info <- read_csv("samplesheet_kd_r2.csv") %>% 
  mutate(sample = as.character(sample))
coldata <- data.frame(files = list.files(path = "star_salmon/", pattern = "quant.sf", recursive = TRUE)) %>% 
  mutate(names = str_replace_all(files, "/.*", ""), 
         files = paste0("star_salmon/", files)) %>% 
  left_join(sample_info, by = c("names" = "sample")) %>% 
  select(!starts_with("fastq")) %>% 
  mutate(sample_num = as.numeric(names), 
         sample_id = paste0("X", names),
         group = paste(Protocol, Treatment, sep = "_"), .after = names) %>% 
  arrange(sample_num)
rownames(coldata) <- coldata$names
file.exists(coldata$files)
table(coldata$Treatment, coldata$Protocol)
## Remove outlier - X340
coldata <- coldata %>% 
  filter(names != "340")
## Load Tx to Gene table
tx2gene <- read.delim("star_salmon/salmon_tx2gene.tsv", sep = "\t", header = FALSE)
colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
gene_info <- unique(tx2gene[, -1])
rownames(gene_info) <- gene_info$gene_id
## Load TPM table and count stats
ds_tpm <- read.delim("star_salmon/salmon.merged.gene_tpm.tsv", 
                     sep = "\t", header = TRUE, row.names = 1)
ds_tpm <- ds_tpm[, coldata$sample_id]
ds_tpm_mean <- t(apply(ds_tpm, 1, function(x) tapply(x, coldata$Treatment, mean)))
## Collect files info
files <- coldata$files
names(files) <- coldata$names
## Load Salmon data and make DESeq object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# Combinatorial LM
dds <- DESeqDataSetFromTximport(txi, coldata, ~ 0 + Treatment + Protocol)
## Pre-filtering with DESeq2
keep <- rowSums(counts(dds)) >= 10; table(keep)
## Pre-filtering with edgeR (remove low expressed genes)
keep <- filterByExpr(y = counts(dds), 
                     group = paste(coldata$Treatment, coldata$Protocol, sep = "_"),
                     min.count = 10); table(keep)
dds <- dds[keep, ]
## Add gene names
rowData(dds)$geneSymbol <- gene_info[rownames(dds), "gene_name"]
table(duplicated(rowData(dds)$geneSymbol))
## Perform DESeq analysis 
dds <- DESeq(dds)
resultsNames(dds)
## Data transformations and plots
dds_vsd <- vst(dds)
meanSdPlot(assay(dds_vsd))
## PCA
plotPCA(dds_vsd, intgroup = "group") + theme_bw()
plotPCA(dds_vsd, intgroup = "Protocol") + theme_bw() + labs(title = "Protocol") + scale_color_manual(values = mypal_br_short[c(2,5,1)])
plotPCA(dds_vsd, intgroup = "Treatment") + theme_bw() + labs(title = "siRNA Treatment") + scale_color_brewer(palette = "Set2")
plotPCA(dds_vsd, intgroup = "Experiment_ID") + theme_bw() + labs(title = "Experimental batch") + scale_color_brewer(palette = "Set2")
## Heatmap
top_var_500 <- order(rowVars(assay(dds_vsd)), decreasing = TRUE)[1:500]
vsd_cor <- cor(assay(dds_vsd))
vsd_cor <- cor(assay(dds_vsd)[top_var_500, ])
ann_colors <- list(Protocol = c(Baseline = mypal_br_short[2],
                                Cholesterol = mypal_br_short[5],
                                Stretch = mypal_br_short[1]))
                   #Experiment_ID = as.character(mypal_br_long[1:7]),
                   #Treatment = as.character(mypal_treat[1:7]))
p <- pheatmap(vsd_cor, 
              #clustering_method = "ward.D2", 
              labels_row = paste(coldata$names, coldata$Treatment),
              labels_col = paste(coldata$names, coldata$Treatment),
              fontsize_row = 6,
              fontsize_col = 6, #angle_col = 315, 
              annotation_row = coldata[, c("Experiment_ID", "Protocol", "Treatment")],
              annotation_colors = ann_colors)
## Manual PCA
# Screeplot
rv <- rowVars(assay(dds_vsd))
sel <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
dds_pca <- prcomp(t(assay(dds_vsd)[sel, ]))
percentVar <- dds_pca$sdev^2 / sum(dds_pca$sdev^2)
scree_plot <- data.frame(PC = 1:length(percentVar),
                         Variance = percentVar * 100)
scree_plot[1:15, ] %>%  
  ggplot(aes(x = PC, y = Variance)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(x = "Principle component", y = "% of variance explained")
# PCA plots
ntop_genes <- nrow(dds_vsd)
ntop_genes <- 500
pca_data <- DESeq2::plotPCA(dds_vsd, 
                            ntop = ntop_genes,
                            #ntop = 500,
                            intgroup = c("Protocol", "Treatment", "Experiment_ID"), 
                            returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# color by protocol
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = Treatment)) +
  #geom_point(aes(color = Type, fill = Protocol), shape = 21, size = 4, alpha = 0.85) + 
  geom_point(aes(shape = Type, fill = Protocol), size = 2, alpha = 0.8) + 
  scale_shape_manual(values = c(21, 24)) + 
  #scale_color_manual(values = mypal_br_short[2:1]) + 
  scale_fill_manual(values = as.character(mypal_prot)) + 
  geom_text_repel(size = 2, color = "grey50", max.overlaps = 30) + 
  labs(title = paste("PCA based on", ntop_genes, "genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
file_base <- "pca_all-16834_color-protocol"
ggsave(filename = paste0(file_base, ".svg"), height = 7, width = 12)
ggsave(filename = paste0(file_base, ".tiff"), height = 7, width = 12, 
       dpi = 300, compression = "lzw")
# color by treatment
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = name)) +
  geom_point(aes(shape = Type, fill = Treatment), size = 2, alpha = 0.85) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_bright() + 
  #scale_fill_manual(values = as.character(mypal_treat)) + 
  geom_text_repel(size = 2, color = "grey50", max.overlaps = 30) + 
  labs(title = paste("PCA based on", ntop_genes, "genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
file_base <- "pca_all-16834_color-treatment"
ggsave(filename = paste0(file_base, ".svg"), height = 7, width = 12)
ggsave(filename = paste0(file_base, ".tiff"), height = 7, width = 12, 
       dpi = 300, compression = "lzw")
# color by THOR experiment (batch)
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = paste(name, "/", Treatment))) +
  geom_point(aes(shape = Type, fill = THOR_Exp), size = 2, alpha = 0.85) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_bright() + 
  geom_text_repel(size = 2, color = "grey50", max.overlaps = 30) + 
  labs(title = paste("PCA based on", ntop_genes, "genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)"),
       fill = "THOR Exp. ID") + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
file_base <- "pca_top500_color-batch"
ggsave(filename = paste0(file_base, ".svg"), height = 7, width = 12)
ggsave(filename = paste0(file_base, ".tiff"), height = 7, width = 12, 
       dpi = 300, compression = "lzw")
## Exclude sample #25 (THOR 6 (08.06.21) - Cholesterol CTR)
r2_dds <- dds
r2_vsd <- vst(r2_dds)
r2_counts <- counts(r2_dds)
r2_coldata <- coldata[colnames(r2_dds), ] %>% 
  #select(-group) %>% 
  mutate(Treatment = str_replace_all(Treatment, "CTR", "Control")) %>% 
  mutate(Experiment_Short = gsub("HOR ", "", Experiment_ID), .after = Experiment_Date) %>% 
  mutate(Group = paste(Experiment_Short, Treatment, sep = "_"), .after = Treatment)
## Additional gene filtering
# exclude pseudoautosome genes
grep("PAR", rownames(r2_counts), value = T)
r2_counts <- r2_counts[- grep("PAR", rownames(r2_counts)), ]
r2_tpm <- as.matrix(ds_tpm[rownames(r2_counts), paste0("X", colnames(r2_dds))])
r2_tpm_mean <- t(apply(r2_tpm, 1, function(x) tapply(x, r2_coldata$Treatment, mean)))
# investigate duplicated genes
tmp_anno <- gene_info[rownames(r2_counts), ]
table(duplicated(tmp_anno$gene_name))
x <- tmp_anno[duplicated(tmp_anno$gene_name), ]
tmp_dup <- tmp_anno[tmp_anno$gene_name %in% x$gene_name, ] %>% arrange(gene_name)
tmp_dup <- cbind(tmp_dup, r2_tpm_mean[tmp_dup$gene_id, ])
## Save data
save(r2_counts, r2_coldata, r2_dds, r2_tpm, file = "r2_data.rda")
### Simple comparisons within batch
source("deg_functions.R")
r2_results = r2_object <- list()
grcol <- "Treatment"
ds_names <- sort(unique(r2_coldata$group[r2_coldata$Treatment != "Control"]))
for (i in 1:length(ds_names)) {
  sel_exp <- unique(r2_coldata$Experiment_ID[r2_coldata$group == ds_names[i]])
  sel_gene <- gsub(".+_", "", ds_names[i])
  pheno_table <- r2_coldata %>% 
    filter(Experiment_ID == sel_exp & Treatment %in% c("Control", sel_gene))
  pheno_table$Treatment <- factor(pheno_table$Treatment, levels = c("Control", sel_gene))
  count_table <- na.omit(r2_counts[, pheno_table$names])
  res <- list()
  res$deseq <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "deseq2", keep_object = TRUE)
  res$edger <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "edger", keep_object = TRUE)
  res$limma <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "limma", keep_object = TRUE)
  res <- purrr::transpose(res)
  r2_object[[i]] <- res$object
  res_out <- map(res[names(res) != "object"], function(x) data.frame(
    GeneID = rownames(x[[1]]), 
    GeneName = gene_info[rownames(x[[1]]), "gene_name"], 
    do.call("cbind", x)))[[1]]
  # add TPM
  tpm_table <- ds_tpm[res_out$GeneID, pheno_table$sample_id]
  tpm_mean <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, mean)))
  tpm_sd <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, sd)))
  r2_results[[i]] <- data.frame(res_out[, 1:2],
                                TPM.Mean = tpm_mean, TPM.SD = tpm_sd,
                                res_out[, 7:ncol(res_out)])
}; names(r2_results) = names(r2_object) <- ds_names
rm(sel_exp, sel_gene, pheno_table, count_table, res, res_out, tpm_table, tpm_mean, tpm_sd)
save(r2_results, r2_object, file = "r2_3methods_degs.rda")
## Save results to Excel files
options(java.parameters = "-Xmx4096m")
library(xlsx)
kd_names <- sort(unique(r2_coldata$Treatment[r2_coldata$Treatment != "Control"]))
dir.create("dea_results")
for(i in 1:length(kd_names)) {
  # save txt
  write.table(r2_results[[paste0("Baseline_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Baseline.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(r2_results[[paste0("Cholesterol_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Cholesterol.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(r2_results[[paste0("Stretch_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Stretch.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  # save xlsx
  write.xlsx(r2_results[[paste0("Baseline_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Baseline", append = TRUE)
  write.xlsx(r2_results[[paste0("Cholesterol_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Cholesterol", append = TRUE)
  write.xlsx(r2_results[[paste0("Stretch_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Stretch", append = TRUE)
}
## Get significant results
fdr_th <- 0.05
lfc_th <- 0.5
r2_signif <- list()
r2_signif$deseq <- r2_results %>% 
  map(~ filter(.x, deseq.padj < fdr_th, abs(deseq.log2FoldChange) >= lfc_th))
#lapply(r2_signif$deseq, nrow)
r2_signif$edger <- r2_results %>% 
  map(~ filter(.x, edger.FDR < fdr_th, abs(edger.logFC) >= lfc_th))
r2_signif$limma <- r2_results %>% 
  map(~ filter(.x, limma.adj.P.Val < fdr_th, abs(limma.logFC) >= lfc_th))

## Count number of DEGs
r2_deg_stat <- data.frame(group = names(r2_results),
                          deseq_all = map_dbl(r2_signif$deseq, nrow),
                          deseq_up = map_dbl(r2_signif$deseq, ~ sum(.x$deseq.log2FoldChange > 0)),
                          deseq_down = map_dbl(r2_signif$deseq, ~ sum(.x$deseq.log2FoldChange < 0)),
                          edger_all = map_dbl(r2_signif$edger, nrow),
                          edger_up = map_dbl(r2_signif$edger, ~ sum(.x$edger.logFC > 0)),
                          edger_down = map_dbl(r2_signif$edger, ~ sum(.x$edger.logFC < 0)),
                          limma_all = map_dbl(r2_signif$limma, nrow),
                          limma_up = map_dbl(r2_signif$limma, ~ sum(.x$limma.logFC > 0)),
                          limma_down = map_dbl(r2_signif$limma, ~ sum(.x$limma.logFC < 0)))
r2_deg_stat %>% 
  separate(group, into = c("Protocol", "Treatment"), sep = "_") %>% 
  write_delim(file = "res_deg_number_fdr0.05_lfc0.5.tsv", delim = "\t")
## Automatic GSEA with WebGestalt 
library(WebGestaltR)
View(listGeneSet("hsapiens"))
mymethod <- "deseq"
mycolumns <- c("GeneName", "deseq.log2FoldChange")
mypath <- paste0(mymethod, "_gsea_webgestalt")
dir.create(mypath)
myds <- names(r2_signif[[mymethod]])
x <- myds[18]
for (x in myds) {
  print(x)
  xpath <- file.path(mypath, x)
  dir.create(xpath)
  xdat <- r2_signif[[mymethod]][[x]][, mycolumns]
  xdat <- xdat %>% arrange(GeneName)
  xfile <- paste0(xpath, "/", x, ".rnk")
  write_delim(xdat, file = xfile, delim = "\t", col_names = FALSE)
  # Gene Ontology - BP
  WebGestaltR(outputDirectory = xpath, projectName = "GO-BP",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Biological_Process_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - MF
  WebGestaltR(outputDirectory = xpath, projectName = "GO-MF",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Molecular_Function_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - CC
  WebGestaltR(outputDirectory = xpath, projectName = "GO-CC",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Cellular_Component_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # KEGG Pathway
  WebGestaltR(outputDirectory = xpath, projectName = "KEGG",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_KEGG",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Reactome
  WebGestaltR(outputDirectory = xpath, projectName = "Reactome",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Reactome",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # WikiPathway
  WebGestaltR(outputDirectory = xpath, projectName = "WikiPathway",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Wikipathway",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
}; rm(x, xpath, xdat, xfile)


#### SEQUENCING BATCH 3 ####

## Set working dir
setwd("thor_kd_r3")
## Load samplesheet
# "samplesheet_kd_r3.csv", a samplesheet file prepared for 'nf-core/rnaseq' pipeline
sample_info <- read_csv("samplesheet_kd_r3.csv") %>% 
  mutate(sample = as.character(sample))
coldata <- data.frame(files = list.files(path = "star_salmon/", pattern = "quant.sf", recursive = TRUE)) %>% 
  mutate(names = str_replace_all(files, "/.*", ""), 
         files = paste0("star_salmon/", files)) %>% 
  left_join(sample_info, by = c("names" = "sample")) %>% 
  select(!starts_with("fastq")) %>% 
  mutate(sample_num = as.numeric(names), 
         sample_id = paste0("X", names),
         group = paste(Protocol, Treatment, sep = "_"), .after = names) %>% 
  arrange(sample_num)
rownames(coldata) <- coldata$names
file.exists(coldata$files)
table(coldata$Treatment, coldata$Protocol)
## Remove outlier - X356
coldata <- coldata %>% 
  filter(names != "356")
## Load Tx to Gene table
tx2gene <- read.delim("star_salmon/salmon_tx2gene.tsv", sep = "\t", header = FALSE)
colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
gene_info <- unique(tx2gene[, -1])
rownames(gene_info) <- gene_info$gene_id
## Load TPM table and count stats
ds_tpm <- read.delim("star_salmon/salmon.merged.gene_tpm.tsv", 
                     sep = "\t", header = TRUE, row.names = 1)
ds_tpm <- ds_tpm[, coldata$sample_id]
ds_tpm_mean <- t(apply(ds_tpm, 1, function(x) tapply(x, coldata$Treatment, mean)))
## Collect files info
files <- coldata$files
names(files) <- coldata$names
## Load salmon data and make DESeq object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
## Combinatorial LM
dds <- DESeqDataSetFromTximport(txi, coldata, ~ 0 + Treatment)
## Pre-filtering with DESeq2
keep <- rowSums(counts(dds)) >= 10; table(keep)
## Pre-filtering with edgeR (remove low expressed genes)
keep <- filterByExpr(y = counts(dds), group = paste(coldata$Treatment, coldata$Protocol, sep = "_"), min.count = 10); table(keep)
dds <- dds[keep, ]
## Add gene names
rowData(dds)$geneSymbol <- gene_info[rownames(dds), "gene_name"]
table(duplicated(rowData(dds)$geneSymbol))
## Perform DESeq analysis 
dds <- DESeq(dds)
resultsNames(dds)
## Data transformations and plots
dds_vsd <- vst(dds)
meanSdPlot(assay(dds_vsd))
plotPCA(dds_vsd, intgroup = "group") + theme_bw()
plotPCA(dds_vsd, intgroup = "Treatment", returnData = TRUE) %>% 
  ggplot(aes(PC1, PC2, color = group, label = name)) + 
  geom_point(size = 3) + 
  geom_text_repel(color = "dimgrey") + 
  theme_bw() + 
  labs(title = "siRNA Treatment") + 
  scale_color_brewer(palette = "Set2")
## Explore using pcaExplorer
pcaExplorer(dds = dds, dst = dds_vsd, annotation = gene_info)
## Heatmap
top_var_500 <- order(rowVars(assay(dds_vsd)), decreasing = TRUE)[1:500]
vsd_cor <- cor(assay(dds_vsd))
vsd_cor <- cor(assay(dds_vsd)[top_var_500, ])
ann_colors <- list(Protocol = c(Baseline = mypal_br_short[2],
                                Cholesterol = mypal_br_short[5],
                                Stretch = mypal_br_short[1]))
                   #Experiment_ID = as.character(mypal_br_long[1:7]),
                   #Treatment = as.character(mypal_treat[1:7]))
p <- pheatmap(vsd_cor, 
              #clustering_method = "ward.D2", 
              labels_row = paste(coldata$names, coldata$Treatment),
              labels_col = paste(coldata$names, coldata$Treatment),
              fontsize_row = 8,
              fontsize_col = 8, #angle_col = 315, 
              annotation_row = coldata[, c("Experiment_ID", "Protocol", "Treatment")],
              annotation_colors = ann_colors)
## Manual PCA
# Screeplot
rv <- rowVars(assay(dds_vsd))
sel <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
dds_pca <- prcomp(t(assay(dds_vsd)[sel, ]))
percentVar <- dds_pca$sdev^2 / sum(dds_pca$sdev^2)
scree_plot <- data.frame(PC = 1:length(percentVar),
                         Variance = percentVar * 100)
scree_plot[1:15, ] %>%  
  ggplot(aes(x = PC, y = Variance)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(x = "Principle component", y = "% of variance explained")
# PCA plots
ntop_genes <- nrow(dds_vsd)
ntop_genes <- 500
pca_data <- DESeq2::plotPCA(dds_vsd,ntop = ntop_genes, intgroup = c("Protocol", "Treatment", "Experiment_ID"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
# color by treatment
pca_data %>% 
  mutate(Type = ifelse(Treatment == "CTR", "Control", "KD"),
         THOR_Exp = factor(as.numeric(gsub("THOR ", "", Experiment_ID)))) %>% 
  ggplot(aes(PC1, PC2, label = name)) + 
  geom_point(aes(shape = Type, color = Treatment), size = 5, alpha = 0.8) + 
  geom_text_repel(size = 4, color = "grey50", max.overlaps = 30) + 
  labs(title = paste("PCA based on", ntop_genes, "genes"),
       x = paste0("PC1 (", percent_var[1], "% of variance)"),
       y = paste0("PC2 (", percent_var[2], "% of variance)")) + 
  theme_bw()
## DEA using 3 methods
r3_dds <- dds
r3_vsd <- vst(r3_dds)
r3_counts <- counts(r3_dds)
r3_coldata <- coldata[colnames(r3_dds), ] %>% 
  mutate(Treatment = recode(Treatment, CTR = "Control"))
## Additional gene filtering
# exclude pseudoautosome genes
grep("PAR", rownames(r3_counts), value = T)
r3_counts <- r3_counts[- grep("PAR", rownames(r3_counts)), ]
r3_tpm <- as.matrix(ds_tpm[rownames(r3_counts), paste0("X", colnames(r3_dds))])
r3_tpm_mean <- t(apply(r3_tpm, 1, function(x) tapply(x, r3_coldata$Treatment, mean)))
# investigate duplicated genes
tmp_anno <- gene_info[rownames(r3_counts), ]
table(duplicated(tmp_anno$gene_name))
x <- tmp_anno[duplicated(tmp_anno$gene_name), ]
tmp_dup <- tmp_anno[tmp_anno$gene_name %in% x$gene_name, ] %>% arrange(gene_name)
tmp_dup <- cbind(tmp_dup, r3_tpm_mean[tmp_dup$gene_id, ])
## Save data
save(r3_counts, r3_coldata, r3_dds, r3_tpm, file = "r3_data.rda")
## Export count matrix
# keep only the most variable genes
#library(matrixStats)
r3_tpm_tsd <- rowVars(as.matrix(r3_tpm))
tmp_tpm <- r3_tpm[order(r3_tpm_tsd, decreasing = TRUE), ]
tmp_fd <- gene_info[rownames(tmp_tpm), ]
tmp_fd <- tmp_fd[!(duplicated(tmp_fd$gene_name)), ]
tmp_pd <- r3_coldata %>% 
  select(sample_id, Experiment_ID:Initial_Conc) %>% 
  mutate(Sample_ID = sample_id, .before = 1) %>% 
  select(-sample_id)
colnames(tmp_pd)[1] <- "Sample_ID"
tmp_tpm <- data.frame(Gene = tmp_fd$gene_name,
                      tmp_tpm[tmp_fd$gene_id, tmp_pd$Sample_ID])
write.table(tmp_tpm, file = "thor_rnaseq_r3_tpm.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tmp_pd, file = "thor_rnaseq_r3_metadata.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
## Simple comparisons within batch
source("deg_functions.R")
r3_results = r3_object <- list()
grcol <- "Treatment"
ds_names <- sort(unique(r3_coldata$group[r3_coldata$Treatment != "Control"]))
for (i in 1:length(ds_names)) {
  sel_exp <- unique(r3_coldata$Experiment_ID[r3_coldata$group == ds_names[i]])
  sel_gene <- gsub(".+_", "", ds_names[i])
  pheno_table <- r3_coldata %>% 
    filter(Experiment_ID == sel_exp & Treatment %in% c("Control", sel_gene))
  pheno_table$Treatment <- factor(pheno_table$Treatment, levels = c("Control", sel_gene))
  count_table <- na.omit(r3_counts[, pheno_table$names])
  res <- list()
  res$deseq <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "deseq2", keep_object = TRUE)
  res$edger <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "edger", keep_object = TRUE)
  res$limma <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                                method = "limma", keep_object = TRUE)
  res <- purrr::transpose(res)
  r3_object[[i]] <- res$object
  res_out <- map(res[names(res) != "object"], function(x) data.frame(
    GeneID = rownames(x[[1]]), 
    GeneName = gene_info[rownames(x[[1]]), "gene_name"], 
    do.call("cbind", x)))[[1]]
  # add TPM
  tpm_table <- ds_tpm[res_out$GeneID, pheno_table$sample_id]
  tpm_mean <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, mean)))
  tpm_sd <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Treatment, sd)))
  r3_results[[i]] <- data.frame(res_out[, 1:2],
                                TPM.Mean = tpm_mean, TPM.SD = tpm_sd,
                                res_out[, 7:ncol(res_out)])
}; names(r3_results) = names(r3_object) <- ds_names
rm(sel_exp, sel_gene, pheno_table, count_table, res, res_out, tpm_table, tpm_mean, tpm_sd)
save(r3_results, r3_object, file = "r3_3methods_degs.rda")
## Save results to Excel tables
options(java.parameters = "-Xmx4096m")
library(xlsx)
kd_names <- sort(unique(r3_coldata$Treatment[r3_coldata$Treatment != "Control"]))
dir.create("dea_results")
for(i in 1:length(kd_names)) {
  # save txt
  write.table(r3_results[[paste0("Baseline_", kd_names[i])]], 
              file = paste0("dea_results/", kd_names[i], "_Baseline.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  # save xlsx
  write.xlsx(r3_results[[paste0("Baseline_", kd_names[i])]], 
             file = paste0("dea_results/", kd_names[i], "_Results.xlsx"), 
             row.names = FALSE, sheetName = "Baseline", append = TRUE)
}
## Get significant results
fdr_th <- 0.05
lfc_th <- 0.5
r3_signif <- list()
r3_signif$deseq <- r3_results %>% 
  map(~ filter(.x, deseq.padj < fdr_th, abs(deseq.log2FoldChange) >= lfc_th))
#lapply(r3_signif$deseq, nrow)
r3_signif$edger <- r3_results %>% 
  map(~ filter(.x, edger.FDR < fdr_th, abs(edger.logFC) >= lfc_th))
r3_signif$limma <- r3_results %>% 
  map(~ filter(.x, limma.adj.P.Val < fdr_th, abs(limma.logFC) >= lfc_th))
## Count number of DEGs
r3_deg_stat <- data.frame(group = names(r3_results),
                          deseq_all = map_dbl(r3_signif$deseq, nrow),
                          deseq_up = map_dbl(r3_signif$deseq, ~ sum(.x$deseq.log2FoldChange > 0)),
                          deseq_down = map_dbl(r3_signif$deseq, ~ sum(.x$deseq.log2FoldChange < 0)),
                          edger_all = map_dbl(r3_signif$edger, nrow),
                          edger_up = map_dbl(r3_signif$edger, ~ sum(.x$edger.logFC > 0)),
                          edger_down = map_dbl(r3_signif$edger, ~ sum(.x$edger.logFC < 0)),
                          limma_all = map_dbl(r3_signif$limma, nrow),
                          limma_up = map_dbl(r3_signif$limma, ~ sum(.x$limma.logFC > 0)),
                          limma_down = map_dbl(r3_signif$limma, ~ sum(.x$limma.logFC < 0)))
r3_deg_stat %>% 
  separate(group, into = c("Protocol", "Treatment"), sep = "_") %>% 
  write_delim(file = "res_deg_number_fdr0.05_lfc0.5.tsv", delim = "\t")
## Automatic GSEA with WebGestalt 
library(WebGestaltR)
View(listGeneSet("hsapiens"))
mymethod <- "deseq"
mycolumns <- c("GeneName", "deseq.log2FoldChange")
mypath <- paste0(mymethod, "_gsea_webgestalt")
dir.create(mypath)
myds <- names(r3_signif[[mymethod]])
#x <- myds[18]
for (x in myds) {
  print(x)
  xpath <- file.path(mypath, x)
  dir.create(xpath)
  xdat <- r3_signif[[mymethod]][[x]][, mycolumns]
  xdat <- xdat %>% arrange(GeneName)
  xfile <- paste0(xpath, "/", x, ".rnk")
  write_delim(xdat, file = xfile, delim = "\t", col_names = FALSE)
  # Gene Ontology - BP
  WebGestaltR(outputDirectory = xpath, projectName = "GO-BP",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Biological_Process_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - MF
  WebGestaltR(outputDirectory = xpath, projectName = "GO-MF",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Molecular_Function_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Gene Ontology - CC
  WebGestaltR(outputDirectory = xpath, projectName = "GO-CC",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "geneontology_Cellular_Component_noRedundant",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # KEGG Pathway
  WebGestaltR(outputDirectory = xpath, projectName = "KEGG",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_KEGG",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # Reactome
  WebGestaltR(outputDirectory = xpath, projectName = "Reactome",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Reactome",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
  # WikiPathway
  WebGestaltR(outputDirectory = xpath, projectName = "WikiPathway",
              enrichMethod = "GSEA", organism = "hsapiens", 
              enrichDatabase = "pathway_Wikipathway",
              interestGeneFile = xfile, interestGeneType = "genesymbol",
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05)
}; rm(x, xpath, xdat, xfile)


#### BASELINE EXPERIMENTS ####

## Load  info about control (R0) samples
tmp <- read_csv(file = "smc_kd_zero/samplesheet.csv", col_types = "c?")
r0_coldata <- data.frame(files = list.files(path = "smc_kd_zero/star_salmon/", pattern = "quant.sf", recursive = TRUE)) %>% 
  mutate(names = str_replace_all(files, "/.*", ""), 
         files = paste0("smc_kd_zero/star_salmon/", files)) %>% 
  left_join(tmp, by = c("names" = "sample")) %>% 
  select(!starts_with("fastq")) %>% 
  mutate(sample_num = as.numeric(names), 
         sample_id = paste0("X", names),
         group = paste(Protocol, Treatment, sep = "_"), .after = names) %>% 
  arrange(sample_num) %>% 
  mutate(Experiment_Short = str_replace_all(Experiment_ID, "HOR ", ""), .after = Experiment_Date) %>% 
  mutate(Group = paste(Experiment_Short, Treatment, sep = "/"), .after = Treatment) %>% 
  mutate(Seq_Batch = "R0")
file.exists(r0_coldata$files)
## Load Tx to Gene table
tx2gene <- read.delim("star_salmon/salmon_tx2gene.tsv", sep = "\t", header = FALSE)
colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
gene_info <- unique(tx2gene[, -1])
rownames(gene_info) <- gene_info$gene_id
## Load TPM table and count stats
ds_tpm <- read.delim("star_salmon/salmon.merged.gene_tpm.tsv", 
                     sep = "\t", header = TRUE, row.names = 1)
ds_tpm <- ds_tpm[, r0_coldata$sample_id]
## Gather files info
files <- r0_coldata$files
names(files) <- r0_coldata$sample_id
## Load salmon data and make DESeq object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# Combinatorial LM
dds <- DESeqDataSetFromTximport(txi, r0_coldata, ~ 0 + Experiment_Short + Protocol)
colData(dds)$Exp_Protocol <- paste(colData(dds)$Experiment_Short, colData(dds)$Protocol)
## Add gene names
rowData(dds)$geneSymbol <- gene_info[rownames(dds), "gene_name"]
## VST normalisation
dds_vsd <- vst(dds)
## PCA
plotPCA(dds_vsd, intgroup = "Protocol")
plotPCA(dds_vsd, intgroup = "Exp_Protocol") + theme_bw()
library(vsn)
library(ggrepel)
dds_vsd <- vst(dds[, which(colData(dds)$Experiment_Short == "T02")])
rownames(dds_vsd) <- make.unique(rowData(dds_vsd)$geneSymbol)
meanSdPlot(assay(dds_vsd))
## Fast PCA
dds_pca <- plotPCA(dds_vsd, intgroup = "Protocol")
dds_pca$data$Experiment = colData(dds_vsd)$Experiment_ID
p <- dds_pca$data %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Protocol, shape = Experiment), 
             color = "black", #shape = 21, 
             alpha = 1, size = 4) + 
  scale_fill_manual(values = as.character(mypal_prot[c(1,2,3)])) + 
  scale_shape_manual(values = c(21,24)) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(x = dds_pca$labels$x, y = dds_pca$labels$y, 
       title = "PCA based on top 500 most variable genes") + 
  theme_bw()
p
## Get DEGs 
r0_results = r0_object <- list()
# THOR 01: Cholesterol vs Baseline 
pheno_table <- r0_coldata %>% 
  filter(Experiment_Short == "T01")
rownames(pheno_table) <- pheno_table$sample_id
count_table <- txi$counts[, pheno_table$sample_id]
grcol <- "Protocol"
res <- list()
res$deseq <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "deseq2", keep_object = TRUE)
res$edger <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "edger", keep_object = TRUE)
res$limma <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "limma", keep_object = TRUE)
res <- purrr::transpose(res)
r0_object[["Cholesterol_vs_Baseline"]] <- res$object
res_out <- map(res[names(res) != "object"], function(x) data.frame(
  GeneID = rownames(x[[1]]), 
  GeneName = gene_info[rownames(x[[1]]), "gene_name"], 
  do.call("cbind", x)))[[1]]
# add TPM
tpm_table <- ds_tpm[res_out$GeneID, pheno_table$sample_id]
tpm_mean <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Protocol, mean)))
tpm_sd <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Protocol, sd)))
r0_results[["Cholesterol_vs_Baseline"]] <- data.frame(res_out[, 1:2],
                              TPM.Mean = tpm_mean, TPM.SD = tpm_sd,
                              res_out[, 7:ncol(res_out)])
rm(pheno_table, count_table, res, res_out, tpm_table, tpm_mean, tpm_sd)
# THOR 02: Stretch vs Baseline 
pheno_table <- r0_coldata %>% 
  filter(Experiment_Short == "T02")
rownames(pheno_table) <- pheno_table$sample_id
count_table <- txi$counts[, pheno_table$sample_id]
grcol <- "Protocol"
res <- list()
res$deseq <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "deseq2", keep_object = TRUE)
res$edger <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "edger", keep_object = TRUE)
res$limma <- getDEGFromCounts(count_table, pheno_table, group_column = grcol, 
                              method = "limma", keep_object = TRUE)
res <- purrr::transpose(res)
r0_object[["Stretch_vs_Baseline"]] <- res$object
res_out <- map(res[names(res) != "object"], function(x) data.frame(
  GeneID = rownames(x[[1]]), 
  GeneName = gene_info[rownames(x[[1]]), "gene_name"], 
  do.call("cbind", x)))[[1]]
# add TPM
tpm_table <- ds_tpm[res_out$GeneID, pheno_table$sample_id]
tpm_mean <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Protocol, mean)))
tpm_sd <- t(apply(tpm_table, 1, function(x) tapply(x, pheno_table$Protocol, sd)))
r0_results[["Stretch_vs_Baseline"]] <- data.frame(res_out[, 1:2],
                                                      TPM.Mean = tpm_mean, TPM.SD = tpm_sd,
                                                      res_out[, 7:ncol(res_out)])
rm(pheno_table, count_table, res, res_out, tpm_table, tpm_mean, tpm_sd)
## Save results
r0_dds <- dds
r0_tpm <- ds_tpm
r0_counts <- counts(r0_dds)
save(r0_counts, r0_coldata, r0_dds, r0_tpm, file = "r0_data.rda")
save(r0_results, r0_object, file = "r0_3methods_degs.rda")
## Save DEGs to Excel
options(java.parameters = "-Xmx4096m")
library(xlsx)
write.xlsx(r0_results$Cholesterol_vs_Baseline, 
           file = paste0("Untreated_DEA_Results.xlsx"), 
           row.names = FALSE, sheetName = "Cholesterol_vs_Baseline", append = TRUE)
write.xlsx(r0_results$Stretch_vs_Baseline, 
           file = paste0("Untreated_DEA_Results.xlsx"), 
           row.names = FALSE, sheetName = "Stretch_vs_Baseline", append = TRUE)


#### COMBINE SEQUENCING BATCHES ####

## Load sequencing batches 1 and 2
load("smc_kd_round_1/r1_data.rda")
r1_coldata <- r1_coldata %>% 
  mutate(Experiment_Short = str_replace_all(Experiment_ID, "HOR ", ""), .after = Experiment_Date) %>% 
  mutate(Group = paste(Experiment_Short, Treatment, sep = "/"), .after = Treatment)
load("smc_kd_round_2/r2_data.rda")
r2_coldata <- r2_coldata %>%
  mutate(Treatment = str_replace_all(Treatment, "Control", "CTR"),
         Group = paste(Experiment_Short, Treatment, sep = "/"))
## Load remade results for THOR 30 experiment (sequencing batch 3)
load("smc_kd_round_3/r3_data.rda")
r3_coldata <- r3_coldata %>% 
  mutate(Treatment = str_replace_all(Treatment, "Control", "CTR")) %>% 
  mutate(Experiment_Short = str_replace_all(Experiment_ID, "HOR ", ""), .after = Experiment_Date) %>% 
  mutate(Group = paste(Experiment_Short, Treatment, sep = "/"), .after = Treatment)
## Combine sample info
r1_coldata <- r1_coldata %>% 
  select(1:15) %>% 
  mutate(files = paste0("smc_kd_round_1/", files), Seq_Batch = "R1")
r2_coldata <- r2_coldata %>% 
  select(1:15) %>% 
  mutate(files = paste0("smc_kd_round_2/", files), Seq_Batch = "R2")
r3_coldata <- r3_coldata %>% 
  select(1:15) %>% 
  mutate(files = paste0("smc_kd_round_3/", files), Seq_Batch = "R3")
rc_coldata <- rbind(r0_coldata, r1_coldata, r2_coldata, r3_coldata) %>% 
  mutate(sample_id = paste(Experiment_Short, sample_num, sep = "_"))
## Create combined matrix and DESeq object
## Gather files info
files <- rc_coldata$files
names(files) <- rc_coldata$sample_id
rc_coldata <- rc_coldata %>% 
  select(c(3,4,8:16)) %>% 
  mutate(Exp_Group = Group,
         Group = paste(Experiment_Short, Protocol, Treatment, sep = "_"),
         num = 1:n())
rownames(rc_coldata) <- rc_coldata$sample_id
## Load Tx to Gene table
tx2gene <- read.delim("smc_kd_round_1/star_salmon/salmon_tx2gene.tsv", sep = "\t", header = FALSE)
colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
gene_info <- unique(tx2gene[, -1])
rownames(gene_info) <- gene_info$gene_id
## Load salmon data and make DESeq object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
## Make DESeq2 object
rc_dds <- DESeqDataSetFromTximport(txi, rc_coldata, ~ 0 + Group)
## save some objects
saveRDS(rc_dds, file = "thor_rnaseq_r0r1r2r3_deseq2_object.rds")
rc_dds <- readRDS("thor_rnaseq_r0r1r2r3_deseq2_object.rds")
rc_dds_vsd <- vst(rc_dds)
## Exclude THOR30, TCF21 and PLPP3 
sample_info <- sample_info %>% 
  filter(!(Treatment %in% c("PLPP3", "TCF21") | Experiment_ID == "THOR 30"))
count_matrix <- count_matrix[, sample_info$sample_id]
tpm_matrix <- tpm_matrix[, sample_info$sample_id]
rc_dds <- rc_dds[, sample_info$sample_id]
rc_dds_vsd <- vst(rc_dds)
rc_coldata <- sample_info
## Save data matrices for publication
keep <- !(rowSums(count_matrix == 0) == ncol(count_matrix))
table(keep)
count_matrix <- count_matrix[keep, ]
tpm_matrix <- tpm_matrix[keep, ]
## create dir for Zenodo
dir.create("zenodo_upload")
## Save count matrix
data.frame(count_matrix) %>% 
  rownames_to_column("Gene") %>% 
  write_tsv(file = "zenodo_upload/count_matrix.tsv")
## Save TPM matrix
data.frame(tpm_matrix) %>% 
  rownames_to_column("Gene") %>% 
  write_tsv(file = "zenodo_upload/tpm_matrix.tsv")
## Correct and save sample information
meta_info <- sample_info %>% 
  select(sample_id, Seq_Batch, Experiment_ID, Protocol, Treatment)
colnames(meta_info)[1:2] <- c("Sample_ID", "Sequencing_Batch")
write_tsv(meta_info, file = "zenodo_upload/sample_metadata.tsv")
## Manual PCA
ntop_genes <- 5000
rv <- rowVars(assay(rc_dds_vsd))
sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop_genes, length(rv)))]
dds_pca <- prcomp(t(assay(rc_dds_vsd)[sel, ]))
pca_data <- cbind(dds_pca$x, rc_coldata[, -1])
percent_var <- dds_pca$sdev^2 / sum(dds_pca$sdev^2)
scree_plot <- data.frame(PC = 1:length(percent_var),
                         Variance = percent_var * 100)
## Screeplot
scree_plot[1:15, ] %>%  
  ggplot(aes(x = PC, y = Variance)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(x = "Principle component", y = "% of variance explained")
## 2D plot
library(ggfx)
plot(pca_data[, 1:5], pch = 19, cex = 0.5, col = mypal_prot[factor(pca_data$Protocol)])
pca_data %>% 
  mutate(Type = ifelse(Treatment %in% c("CTR", "UNTRTD"), "Control", "KD")) %>% 
  ggplot(aes(PC1, PC2, label = Exp_Group)) +
  geom_point(aes(shape = Type, fill = Protocol), size = 2, alpha = 0.8) + 
  scale_shape_manual(values = c(21, 24)) + 
  #scale_color_manual(values = mypal_br_short[2:1]) + 
  scale_fill_manual(values = as.character(mypal_prot)) + 
  geom_text_repel(size = 1.5, color = "grey50", max.overlaps = 10, segment.size = 0.2) + 
  #with_outer_glow(geom_text(aes(label = Treatment), size = 2, colour = "grey50", hjust = 0.5, vjust = 0.5, check_overlap = TRUE), colour = 'white', sigma = 0.5, expand = 0.5) + 
  labs(title = paste("PCA based on top", ntop_genes, "genes"),
       x = paste0("PC1 (", round(percent_var[1] * 100), "% of variance)"),
       y = paste0("PC2 (", round(percent_var[2] * 100), "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))

file_base <- "r0r1r2r3_pca_top5000-color-protocol"
ggsave(filename = paste0(file_base, "_v1.tiff"), height = 8, width = 12, 
       dpi = 300, compression = "lzw")
## tSNE based on PCA (1:30 PCs)
library(Rtsne)
pca_tsne <- Rtsne(pca_data[, 1:15])
plot(pca_tsne$Y, pch = 19, cex = 1, col = mypal_rainbow_fun(25)[factor(pca_data$Experiment_ID)])
plot(pca_tsne$Y, pch = 19, cex = 1, col = mypal_rainbow_fun(5)[factor(pca_data$Protocol)])
text(pca_tsne$Y+0.5, labels = pca_data$Group, cex = 0.5)
## Extract count matrix
tmp_ctm <- txi$counts
## Keep only non-zero genes
table(rowMax(tmp_ctm) > 0)
tmp_ctm <- tmp_ctm[rowMax(tmp_ctm) > 0, ]
## Rename from EnsID to HGNC symbol 
tmp_var <- rowVars(tmp_ctm)
tmp_ctm <- tmp_ctm[order(tmp_var, decreasing = TRUE), ]
tmp_anno <- gene_info[rownames(tmp_ctm), ]
table(duplicated(tmp_anno$gene_name))
# remove dups keeping only the most variable genes 
count_matrix <- tmp_ctm[!(duplicated(tmp_anno$gene_name)), ]
gene_annotation <- gene_info[rownames(count_matrix), ]
rownames(count_matrix) <- gene_annotation$gene_name
count_matrix <- count_matrix[order(rownames(count_matrix)), ]
gene_annotation <- gene_annotation[order(gene_annotation$gene_name), ]
rownames(gene_annotation) <- gene_annotation$gene_name
colnames(gene_annotation)[1] <- "ensembl_gene_id_version"
gene_annotation$ensembl_gene_id <- gsub("\\..*", "", gene_annotation$ensembl_gene_id_version)
gene_annotation <- gene_annotation[, c(2,3,1)]
# add gene info using biomart
library(biomaRt)
listMarts()
listDatasets(useMart("ENSEMBL_MART_FUNCGEN"))
mart <- useMart("ENSEMBL_MART_ENSEMBL")
#View(listDatasets(mart))
martds <- useDataset("hsapiens_gene_ensembl", mart)
#View(listFilters(martds)); View(listAttributes(martds))
tmpbm <- getBM(filters = "ensembl_gene_id", 
               values = gene_annotation$ensembl_gene_id,
               attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", 
                              "chromosome_name", "start_position", "end_position",
                              "description"), 
               mart = martds)
gene_annotation2 <- left_join(gene_annotation, tmpbm, by = c("ensembl_gene_id" = "ensembl_gene_id"))
gene_annotation2 <- gene_annotation2 %>% 
  arrange(entrezgene_id) %>% 
  filter(!duplicated(ensembl_gene_id))
rownames(gene_annotation2) <- gene_annotation2$ensembl_gene_id
gene_annotation2 <- gene_annotation2[gene_annotation$ensembl_gene_id, ]
gene_annotation <- gene_annotation2; rm(gene_annotation2)
# extract TPM matrix and sample info
tpm_matrix <- txi$abundance[gene_annotation$ensembl_gene_id_version, ]
rownames(tpm_matrix) <- gene_annotation$gene_name
sample_info <- rc_coldata
# save data
save(count_matrix, tpm_matrix, sample_info, gene_annotation,
     file = "thor_rnaseq_r0r1r2r3_counts_dups_filtered.rda")
load("thor_rnaseq_r0r1r2r3_counts_dups_filtered.rda")
## Fix some columns in sample info for ShinyCell
sample_info2 <- sample_info %>% 
  mutate(sample_id = str_replace_all(sample_id, "_", "S"),
         Group = str_replace_all(Group, "_", " "),
         Exp_Group = str_replace_all(Exp_Group, "/", " "))
## Make a Seurat object
library(Seurat)
rc_seuobj <- CreateSeuratObject(counts = count_matrix, project = "thor_rnaseq",
                                meta.data = sample_info2)
rc_seuobj$sample_type <- ifelse(rc_seuobj$Treatment == "UNTRTD", "Control", 
                                ifelse(rc_seuobj$Treatment == "CTR", "siControl", "siTarget"))

Idents(rc_seuobj) <- "Protocol"
rc_seuobj$mito_genes_pct <- PercentageFeatureSet(rc_seuobj, pattern = "^MT-")
rc_seuobj$ribo_genes_pct <- PercentageFeatureSet(rc_seuobj, "^RP[SL]")
VlnPlot(rc_seuobj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "Protocol")
VlnPlot(rc_seuobj, features = c("mito_genes_pct", "ribo_genes_pct"), 
        group.by = "Protocol")
FeatureScatter(rc_seuobj, feature1 = "nCount_RNA", feature2 = "mito_genes_pct",
               group.by = "Protocol")
FeatureScatter(rc_seuobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "Protocol")
rc_seuobj <- NormalizeData(rc_seuobj)
rc_seuobj <- FindVariableFeatures(rc_seuobj, nfeatures = 5000)
rc_seuobj <- ScaleData(rc_seuobj, features = rownames(rc_seuobj))
rc_seuobj <- RunPCA(rc_seuobj)
# select the number of components
rc_seuobj <- JackStraw(rc_seuobj, num.replicate = 100, dims = 50)
rc_seuobj <- ScoreJackStraw(rc_seuobj, dims = 1:50)
JackStrawPlot(rc_seuobj, dims = 1:50) ## at maximum, we can use up to 36 PCs
ElbowPlot(rc_seuobj, ndims = 50)
rc_seuobj <- FindNeighbors(rc_seuobj, dims = 1:40)
rc_seuobj <- FindClusters(rc_seuobj, resolution = 1)
rc_seuobj <- RunUMAP(rc_seuobj, dims = 1:40)
# save seurat object
saveRDS(rc_seuobj, file = "thor_rnaseq_seuobj_2.rds")
rc_seuobj <- readRDS("thor_rnaseq_seuobj_2.rds")
## Interactive visualization using ShinyCell
library(ShinyCell)
# Fix metadata names
colnames(rc_seuobj@meta.data)[c(5,11,15)] <- c("Sample_Name", "Protocol_Treatmnet", "Experiment_Treatment")
head(rc_seuobj@meta.data)
# Create and save app
scConf <- createConfig(rc_seuobj, maxLevels = 100)
makeShinyApp(rc_seuobj, scConf, gene.mapping = TRUE,
             shiny.title = "THOR RNA-seq Results (2021-22)")
# Save metadata
data.frame(rc_seuobj@meta.data) %>% 
  rownames_to_column("Sample_ID") %>% 
  write_tsv(file = "thor_rnaseq_seuobj_metadata.tsv")
## Visualize PCA
print(rc_seuobj[["pca"]], dims = 1:2, nfeatures = 500)
VizDimLoadings(rc_seuobj, reduction = "pca", ncol = 2,
               dims = 1:2, nfeatures = 80,  balanced = TRUE)
DimPlot(rc_seuobj, reduction = "pca", group.by = "Protocol")
DimHeatmap(rc_seuobj, dims = 2, reduction = "pca", 
           nfeatures = 80, balanced = TRUE) + scale_fill_viridis()
DoHeatmap(rc_seuobj, group.by = "Exp_Group", assay = "RNA", 
          features = thor_kd_targets$gene_name) + 
  scale_fill_gradientn(colours = mypal_rainbow_fun(100)) + NoLegend()
## Combine DEA statistics data 
## Define limits
p_th <- 0.05
lfc_th <- 0.5
## Load DEA results
load("smc_kd_round_1/r1_3methods_deg.rda")
load("smc_kd_round_2/r2_3methods_deg.rda")
load("smc_kd_round_3/r3_3methods_deg.rda")
load("smc_kd_zero/r0_3methods_deg.rda")
## Replace remade datasets (PLPP3, CRISPLD2, LOXL1)
r2_results$Baseline_CRISPLD2 <- r3_results$Baseline_CRISPLD2
r2_results$Baseline_LOXL1 <- r3_results$Baseline_LOXL1
r2_results$Baseline_PLPP3 <- r3_results$Baseline_PLPP3
## Rename R0 results
names(r0_results) <- c("Cholesterol_UNTRTD", "Stretch_UNTRTD")
## Make a unified comparisons info
rc_dea_info <- data.frame(Group = c(names(r0_results), 
                                           names(r1_results), 
                                           names(r2_results)),
                                 Seq_Batch = c(rep("R0", length(r0_results)), 
                                               rep("R1", length(r1_results)), 
                                               rep("R2", length(r2_results)))
                                 )
rc_dea_info <- rc_dea_info %>% 
  mutate(X = Group) %>% 
  separate(X, sep = "_", into = c("Protocol", "Treatment")) %>% 
  mutate(Comparison = paste0(Treatment, "-KD vs siControl"))
rc_dea_info$Comparison[1:2] <- c("Cholesterol vs Baseline", "Stretch vs Baseline")
tmp <- rc_coldata %>% 
  filter(Treatment != "CTR" | group != "Baseline_UNTRTD") %>% 
  dplyr::select(group, Experiment_ID, Experiment_Date, Experiment_Short) %>% 
  distinct_all()
rc_dea_info <- left_join(rc_dea_info, tmp, by = c("Group" = "group"))
rm(tmp)
## Save combined DEA results
rc_dea_results <- c(r0_results, r1_results, r2_results)
save(rc_dea_results, rc_dea_info, file = "thor_rnaseq_r0r1r2r3_dea_results_list.rda")
## Get DESeq2 Wald statistics
#rc_stat <- tx2gene[!duplicated(tx2gene$gene_id), -1]
#colnames(rc_stat) <- c("GeneID", "GeneName")
rc_stat <- data.frame(GeneID = gene_annotation$ensembl_gene_id_version, 
                      GeneName = gene_annotation$gene_name)
tmp0 <- r0_results %>% 
  map_dfc(~ left_join(rc_stat, .x[, c("GeneID", "deseq.stat")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp1 <- r1_results %>% 
  map_dfc(~ left_join(rc_stat, .x[, c("GeneID", "deseq.stat")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp2 <- r2_results %>% 
  map_dfc(~ left_join(rc_stat, .x[, c("GeneID", "deseq.stat")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
rc_stat <- cbind(rc_stat, tmp0, tmp1, tmp2)
rownames(rc_stat) <- rc_stat$GeneID
keep <- apply(rc_stat[, -c(1,2)], 1, function(x) !(all(is.na(x))))
table(keep)
#View(rc_stat[!keep, ])
rc_stat <- rc_stat[keep, ]
table(duplicated(rc_stat$GeneID))
## Get comparison info
rc_stat_info <- data.frame(Group = colnames(rc_stat[, -c(1,2)]))
rc_stat_info <- rc_stat_info %>% 
  separate(Group, c("Protocol", "Treatment"), sep = "_", remove = FALSE)
# remove Baseline_TCF21
rc_stat_info <- rc_stat_info[rc_stat_info$Group != "Baseline_TCF21", ]
rc_stat <- rc_stat[, colnames(rc_stat) != "Baseline_TCF21"]
## Make fold change matrix
tmp0 <- r0_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.log2FoldChange")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp1 <- r1_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.log2FoldChange")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp2 <- r2_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.log2FoldChange")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
rc_fcdat <- cbind(rc_stat[, 1:2], tmp0, tmp1, tmp2) 
## Make p.adj values matrix
tmp0 <- r0_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.padj")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp1 <- r1_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.padj")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
tmp2 <- r2_results %>% 
  map_dfc(~ left_join(rc_stat[, 1:2], .x[, c("GeneID", "deseq.padj")], by = c("GeneID" = "GeneID"))[, -c(1,2)])
rc_padjdat <- cbind(rc_stat[, 1:2], tmp0, tmp1, tmp2) 
rc_psigndat <- data.frame(rc_padjdat[, 1:2], rc_padjdat[, -c(1,2)] < p_th)
## Save objects
save(rc_stat, rc_stat_info, rc_fcdat, rc_padjdat, file = "thor_rnaseq_r0r1r2r3_dea_results.rda")
## PCA
rc_stat_pca <- prcomp(t(tmp[, -c(1,2)]), center = FALSE, scale = FALSE)
summary(rc_stat_pca)
rc_stat_pca_var <- rc_stat_pca$sdev^2 / sum(rc_stat_pca$sdev^2)
scree_plot <- data.frame(PC = 1:length(rc_stat_pca_var), 
                         Variance = rc_stat_pca_var * 100)
# Screeplot
scree_plot[1:15, ] %>%  
  ggplot(aes(x = PC, y = Variance)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(x = "Principle component", y = "% of variance explained")
# 2D plot
rc_stat_info <- cbind(rc_stat_info, rc_stat_pca$x)
plot(rc_stat_info[, 4:8], pch = 19, cex = 0.5, col = mypal_prot[factor(rc_stat_info$Protocol)])
p1 <- rc_stat_info %>% 
  ggplot(aes(PC1, PC2, label = Treatment)) +
  geom_point(aes(fill = Protocol), shape = 21, size = 3, alpha = 0.7) + 
  scale_fill_manual(values = as.character(mypal_prot)) + 
  geom_text_repel(size = 3, color = "grey50", max.overlaps = 30) + 
  labs(title = paste("PCA based on DESeq2 Wald statistic"),
       x = paste0("PC1 (", round(rc_stat_pca_var[1] * 100), "% of variance)"),
       y = paste0("PC2 (", round(rc_stat_pca_var[2] * 100), "% of variance)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        text = element_text(family = "PT Sans")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
p1
ggsave(p1, filename = "deseq2_stat_pca2.tiff", 
       height = 8, width = 12, 
       dpi = 300, compression = "lzw")
       