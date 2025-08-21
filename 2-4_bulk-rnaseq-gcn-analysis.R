## Load packages
library(tidyverse)
library(ggrepel)
library(DESeq2)
library(RColorBrewer)
library(MetBrewer)
library(patchwork)
library(gridExtra)

## Set color palettes
mypal <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
mypal_br_paired <- RColorBrewer::brewer.pal(12, "Paired")
mypal_prot <- khroma::colour("bright")(4)[c(1,4,2)]
mypal_mbren <- met.brewer("Signac")
brpal_s1 <- brewer.pal(9, "Set1")
brpal_s2 <- brewer.pal(8, "Set2")
brpal_s3 <- brewer.pal(12, "Set3")
brpal_sp <- rev(brewer.pal(11, "Spectral"))
brpal_rb <- rev(brewer.pal(n = 9, name = 'RdBu'))
brpal_byr <- rev(brewer.pal(n = 9, name = 'RdYlBu'))
brpal_reds <- brewer.pal(9, "Reds")
mypal_rb100 <- colorRampPalette(brpal_rb)(100)

### Load data
## Load RNAseq counts
load(file = "thor_rnaseq_r0r1r2r3_counts_dups_filtered.rda")
## Load RNAseq DEA results
load(file = "thor_rnaseq_r0r1r2r3_dea_results_list.rda")
load(file = "thor_rnaseq_r0r1r2r3_dea_results.rda")
rc_dds <- readRDS("thor_rnaseq_r0r1r2r3_deseq2_object.rds")
## Load THOR R4 (IFN, TNF) results
r4_results <- readRDS("thor_r4_results_list.rds")
## Load gene annotation
gencode39_anno <- read_tsv("gencode_v39_genes_annotation.tsv")
gencode39_anno$tag[is.na(gencode39_anno$tag)] <- "none"
gencode39_anno <- as.data.frame(gencode39_anno)
rownames(gencode39_anno) <- gencode39_anno$gene_id
rowData(rc_dds) <- DataFrame(gencode39_anno[rownames(rc_dds), 1:3])
## Exclude THOR30, TCF21 and PLPP3 
sample_info_filt <- sample_info %>% 
  filter(!(Treatment %in% c("PLPP3", "TCF21") | Experiment_ID == "THOR 30"))
rc_dds <- rc_dds[, !(colData(rc_dds)$Treatment %in% c("PLPP3", "TCF21") | colData(rc_dds)$Experiment_ID == "THOR 30")]
rc_stat_info <- rc_stat_info %>% 
  filter(!(Treatment %in% c("PLPP3", "TCF21")))
rc_stat <- cbind(rc_stat[, 1:2], rc_stat[, rc_stat_info$Group])
rc_fcdat <- cbind(rc_fcdat[, 1:2], rc_fcdat[, rc_stat_info$Group])
rc_padjdat <- cbind(rc_padjdat[, 1:2], rc_padjdat[, rc_stat_info$Group])
rc_dea_info <- rc_dea_info %>% 
  filter(!(Treatment %in% c("PLPP3", "TCF21") | Experiment_ID == "THOR 30"))
rc_dea_results <- rc_dea_results[rc_dea_info$Group]
### Define thresholds
p_th <- 0.05
lfc_th <- log2(1.5)


#### All control samples WGCNA (using BioNERO) ####

## Load packages
library(BioNERO)
library(ggrepel)
library(limma)
## Filter the data
all_ctr_samples <- sample_info_filt %>% 
  dplyr::filter(Treatment == "CTR") %>% 
  mutate(Protocol = fct_relevel(Protocol, "Cholesterol"))
all_ctr_dds <- rc_dds[, all_ctr_samples$sample_id]
colData(all_ctr_dds)$Protocol <- fct_relevel(colData(all_ctr_dds)$Protocol, "Cholesterol")
## remove lowly expressed genes: less than 10 counts in 90% samples or more 
keep <- rowSums(counts(all_ctr_dds) >= 10) >= floor(0.1 * ncol(all_ctr_dds))
table(keep)
all_ctr_dds <- all_ctr_dds[keep, ]
rowData(all_ctr_dds) <- DataFrame(gencode39_anno[rownames(all_ctr_dds), 1:3])
rownames(all_ctr_dds) <- make.unique(rowData(all_ctr_dds)$gene_name)
## remove ribosomal genes
table(grepl("^RP[LS]", rownames(all_ctr_dds)))
all_ctr_dds <- all_ctr_dds[!(grepl("^RP[LS]", rownames(all_ctr_dds))), ]
## remove mitochondrial genes
table(grepl("^MT-", rownames(all_ctr_dds)))
all_ctr_dds <- all_ctr_dds[!(grepl("^MT-", rownames(all_ctr_dds))), ]
all_ctr_dds
## VST normalization
design(all_ctr_dds) = ~ 0 + factor(Group)
all_ctr_vst <- vst(all_ctr_dds, blind = FALSE)
p11 <- plotPCA(all_ctr_vst, intgroup = "Protocol") + 
  ggtitle("Coaxed cellular state") + theme_bw() #+ Seurat::NoLegend()
p12 <- plotPCA(all_ctr_vst, intgroup = "Experiment_Short") + 
  ggtitle("Experimental batch") + theme_bw()
p11 / p12
## Advanced PCA
library(PCAtools)
library(ggthemes)
mypca <- pca(assay(all_ctr_vst), metadata = colData(all_ctr_vst))
screeplot(mypca, components = 1:30)
pairsplot(mypca, components = 1:6, colby = "Protocol")
pairsplot(mypca, components = 1:5, colby = "Seq_Batch")
pairsplot(mypca, components = 1:5, colby = "Experiment_Short")
biplot(mypca, colby = "Protocol", legendPosition = "right")
mysize <- 4
p12 <- biplot(mypca, x = "PC5", y = "PC6", 
            colby = "Protocol", shape = "Seq_Batch", 
            colkey = as.character(mypal_prot)[c(2,1,3)], 
            showLoadings = TRUE, ntopLoadings = 10, 
            sizeLoadingsNames = 4, colLoadingsNames = "darkred", 
            colLoadingsArrows = "darkred", 
            drawConnectorsLoadings = TRUE, colConnectorsLoadings = NA, 
            encircle = TRUE, encircleFill = TRUE, #shapekey = c(16,15),
            #lab = mypca$metadata$Treatment, labSize = 2, 
            lab = NULL, 
            pointSize = mysize, legendPosition = "right") + 
  theme_few() + 
  labs(title = paste("PCA based on", nrow(all_ctr_vst), "genes"))
p12
plotloadings(mypca, components = 1:6, labSize = 3, shapeSizeRange = c(5,5),
             col = brpal_rb[c(1,5,9)])
eigencorplot(mypca, metavars = c('Protocol','Seq_Batch','Experiment_Short'), 
             col = brpal_rb)
## Save/load objects
#dir.create("Gene_modules")
save(all_ctr_dds, all_ctr_vst, all_ctr_vstm, file = "Gene_modules/all_ctr_deseq.rda")
#load("Gene_modules/all_ctr_deseq.rda")
## Get DEGs for comparisons between conditions before batch correction
## (it is not possible to remove this batch effect)
dds <- all_ctr_dds
colData(dds)[, "Protocol"] <- factor(colData(dds)[, "Protocol"])
colData(dds)[, "Experiment_Short"] <- factor(colData(dds)[, "Experiment_Short"])
full_model <- model.matrix(~ 0 + Protocol, colData(dds))
#View(full_model)
dds <- DESeq(dds, full = full_model)
resultsNames(dds)
all_ctr_res <- list()
all_ctr_res$Chol_vs_Base <- lfcShrink(
  dds, contrast = list("ProtocolCholesterol", "ProtocolBaseline"), 
  type = "ashr", saveCols = 1:3) %>% as.data.frame() %>% 
  mutate(mlog10_pvalue = - log10(pvalue), mlog10_padj = -log10(padj)) %>% 
  select(gene_id, gene_name, log2FoldChange, lfcSE, baseMean, pvalue, padj, 
         mlog10_pvalue, mlog10_padj, gene_type)
all_ctr_res$Base_vs_Stretch <- lfcShrink(
  dds, contrast = list("ProtocolBaseline", "ProtocolStretch"), 
  type = "ashr", saveCols = 1:3) %>% as.data.frame() %>% 
  mutate(mlog10_pvalue = - log10(pvalue), mlog10_padj = -log10(padj)) %>% 
  select(gene_id, gene_name, log2FoldChange, lfcSE, baseMean, pvalue, padj, 
         mlog10_pvalue, mlog10_padj, gene_type)
all_ctr_res$Chol_vs_Stretch <- lfcShrink(
  dds, contrast = list("ProtocolCholesterol", "ProtocolStretch"), 
  type = "ashr", saveCols = 1:3) %>% as.data.frame() %>% 
  mutate(mlog10_pvalue = - log10(pvalue), mlog10_padj = -log10(padj)) %>% 
  select(gene_id, gene_name, log2FoldChange, lfcSE, baseMean, pvalue, padj, 
         mlog10_pvalue, mlog10_padj, gene_type)
writexl::write_xlsx(all_ctr_res, path = "Gene_modules/siControl_DEA/Batch-corrected_Condition_Comparison.xlsx")

#### GCN analysis ####
library(BioNERO)
## Plot the variance to define the top number of genes to keep by the inflexion point
rv <- matrixStats::rowVars(as.matrix(assay(all_ctr_vst)))
rv2 <- data.frame(Seq = seq(1:nrow(all_ctr_vst)), rowVars = rv[order(rv, decreasing = TRUE)])
my_genes <- c("TAGLN", "LCN2", "ABCG1", "ACTA2", "RGS5", "SOX9")
my_coord <- data.frame(
  gene = my_genes,
  coord = sapply(my_genes, \(x) {
    grep(paste0("^", x, "$"), rownames(rv2))[1]
  })
)
theme_set(theme_bw(base_size = 14) + theme(text = element_text(family = "Arial")))
ggplot(rv2, aes(x = Seq, y = rowVars)) + geom_line() + scale_y_log10() + 
  geom_vline(xintercept = my_coord$coord, color = "red", linetype = "dashed") + 
  geom_text_repel(data = my_coord, aes(x = coord, y = 10, label = gene), 
                  color = "red", fontface = "italic") + 
  labs(title = "VST-normalized gene counts ordered by variance", 
       x = NULL, y = "Log10-scaled variance")
#dev.off()
## Preprocess data (and keep top 5000 genes with highest variances)
all_ctr_filt <- exp_preprocess(all_ctr_vst, 
                               Zk_filtering = TRUE,
                               remove_confounders = FALSE,
                               variance_filter = TRUE, n = 5000)
#colData(all_ctr_filt) <- colData(all_ctr_filt)[, c("Protocol", "Experiment_Short", "Initial_Conc")]
dim(all_ctr_filt)
## Samples considered as outliers
colData(all_ctr_vst)[!(colnames(all_ctr_vst) %in% colnames(all_ctr_filt)), ]
# Heatmap of sample correlations
p <- plot_heatmap(all_ctr_filt, type = "samplecor", 
                  coldata_cols = c("Experiment_Short", "Protocol"),
                  show_rownames = FALSE)
p
# Heatmap of gene expression (here, only the first 50 genes)
p <- plot_heatmap(all_ctr_filt[1:80, ], type = "expr", log_trans = TRUE,
                  coldata_cols = c("Experiment_Short", "Protocol"),
                  show_rownames = TRUE, show_colnames = FALSE)
p
## PCA
library(patchwork)
(plot_PCA(all_ctr_vst, metadata_cols = c("Experiment_Short", "Protocol"), size = 4) + 
    ggtitle("Before BioNERO processing")) / 
  (plot_PCA(all_ctr_filt, metadata_cols = c("Experiment_Short", "Protocol"), size = 4) + 
      ggtitle("After BioNERO processing"))
## Pick the soft threshold beta
library(doParallel)
registerDoParallel(cores = 1)
sft <- SFT_fit(all_ctr_filt, net_type = "signed hybrid", 
               cor_method = "pearson", rsquared = 0.8)
sft$power
sft$plot #& geom_hline(yintercept = 0.85, color = "red")
power <- 7
## Make up a GCN
net <- exp2gcn(all_ctr_filt, net_type = "signed hybrid", verbose = TRUE,
               SFTpower = power, cor_method = "pearson")
## Get table of modules
modules_size <- net$genes_and_modules %>% 
  dplyr::count(Modules) %>% 
  arrange(-n)
modules_size$Alt_name <- paste0("M", sprintf("%02d", 1:nrow(modules_size)))
rownames(modules_size) <- modules_size$Modules
## Get hub genes
hubs <- get_hubs_gcn(all_ctr_filt, net)
head(hubs)
# Explore dendro and modules
plot_dendro_and_colors(net)
# Eigengene networks
plot_eigengene_network(net)
# Number of genes per module.
plot_ngenes_per_module(net)
## Assessing module stability
p_ms <- module_stability(all_ctr_filt, net, nRuns = 30)
## Module-trait associations
all_ctr_metadata <- colData(all_ctr_filt)
colData(all_ctr_filt) <- as.data.frame(all_ctr_metadata) %>% 
  dplyr::select(Protocol, Experiment_Short, Seq_Batch, Initial_Conc) %>% 
  mutate(#Protocol = fct_relevel(factor(Protocol), "Stretch"), 
         ProtocolNum = as.numeric(as.character(fct_recode(
           Protocol, "1" = "Stretch", "2" = "Baseline", "3" = "Cholesterol"))),
         Experiment_Short = as.numeric(factor(Experiment_Short)),
         Seq_Batch = as.numeric(factor(Seq_Batch))) %>% 
  DataFrame(.)
MEtrait <- module_trait_cor(exp = all_ctr_filt, MEs = net$MEs, 
                            metadata_cols = c("Protocol", "Experiment_Short",
                                              "Seq_Batch"))
head(MEtrait)
plot_module_trait_cor(MEtrait)
View(net$genes_and_modules)
## Visualizing module expression profile
plot_expression_profile(
  exp = all_ctr_filt, 
  net = net, 
  plot_module = TRUE, 
  modulename = "cyan"
)
plot_expression_profile(
  exp = all_ctr_filt, 
  net = net, 
  plot_module = TRUE, 
  modulename = "darkslateblue"
)
plot_expression_profile(
  exp = all_ctr_filt, 
  net = net, 
  plot_module = TRUE, 
  modulename = "white"
)
plot_expression_profile(
  exp = all_ctr_filt, 
  net = net, 
  plot_module = TRUE, 
  modulename = "brown"
)
plot_expression_profile(
  exp = all_ctr_filt, 
  net = net, 
  plot_module = TRUE, 
  modulename = "darkgrey"
)
p_exprof <- map(modules_size$Modules[modules_size$Modules != "grey"], ~ {
  plot_expression_profile(
    exp = all_ctr_filt, 
    net = net, 
    plot_module = TRUE, 
    modulename = .x
  ) + ggtitle(toupper(.x))
})
p_exprof_out <- wrap_plots(p_exprof, ncol = 6)
ggsave(p_exprof_out, file = "Gene_modules/all_ctr_sft-7_wgcna/all_modules_expr_profile.png",
       width = 24, height = 16)

## Save WGCNA results
write_rds(list(all_ctr_filt, net), file = "Gene_modules/all_ctr_sft-7_wgcna.rds")


