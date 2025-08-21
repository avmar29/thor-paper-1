library(tidyverse)
library(ggrepel)
library(DESeq2)
library(Seurat)
library(scCustomize)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)

#### Define initial parameters ####

## Set color palettes
mypal <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
mypal_br_paired <- RColorBrewer::brewer.pal(12, "Paired")
mypal_prot <- khroma::colour("bright")(4)[c(1,4,2)]
mypal_dpsc <- DiscretePalette_scCustomize(num_colors = 24, palette = "varibow")
mypal_mbren <- MetBrewer::met.brewer("Signac")
mypal_mviz <- microViz::distinct_palette(pal = "kelly")
## Define thresholds
p_th <- 0.05
lfc_th <- log2(1.5)
## Define seed
set.seed(1)

#### Load data ####

## Load RNAseq counts
load("thor_rnaseq_r0r1r2r3_counts_dups_filtered.rda")
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

## Process DESeq object
dds <- rc_dds
dds <- rc_dds[, !(rc_dds$Experiment_ID %in% c("THOR 01", "THOR 02"))]
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep, ]
#dds <- rc_dds
## Add sample info
colData(dds)$Sample_Type <- ifelse(colData(dds)$Treatment == "UNTRTD", NA, 
  ifelse(colData(dds)$Treatment == "CTR", "Control", "Target KD"))
## Annotate object
gene_annotation <- as.data.frame(gencode39_anno)
gene_annotation <- gene_annotation[rownames(dds), ]
rowRanges(dds) <- GRanges(gene_annotation$hg38_chr,
                          IRanges(start = gene_annotation$hg38_start,
                                  end = gene_annotation$hg38_end), 
                          strand = gene_annotation$strand, 
                          gene_annotation[, c("gene_id", "gene_name", 
                                              "gene_type", "hgnc_id")])
table(duplicated(rowData(dds)$gene_name))
## Order the genes by decreasing variance and remove duplicates
#rowData(dds)$row_num <- 1:nrow(dds)
dds <- dds[order(rowVars(assay(dds)), decreasing = TRUE), ]
dds <- dds[!duplicated(rowData(dds)$gene_name), ]
#dds <- dds[order(rowData(dds)$row_num), ]
rownames(dds) <- rowData(dds)$gene_name
## Get gene counts 
rc_counts <- counts(dds, normalized = FALSE, replaced = FALSE)

## Load scRNA-seq data
seu_comb <- readRDS("seu_comb.rds") # integrated scRNA-seq dataset
sample_info <- read_tsv("scrnaseq_sample_info.tsv") # sample information for scRNA-seq dataset
## Adjust meta data
tmp <- left_join(
  seu_comb@meta.data, sample_info, 
  by = join_by(geo_id == sample_id)
)
rownames(tmp) <- rownames(seu_comb@meta.data)
seu_comb@meta.data <- tmp; rm(tmp)
seu_comb$new_annotation <- factor(seu_comb$cell_subtype)
levels(seu_comb$new_annotation) <- c(
  "T cell", "T cell", "Contractile SMC", "Endothelial cell",
  "Macrophage", "Fibromyocyte", "T cell", "Macrophage",
  "B cell", "T cell", "Monocyte", "Pericyte",
  "Transitional SMC", "Fibroblast", "Mast cell", "Osteochondrogenic SMC",
  "Macrophage", "Dendritic cell", "Undefined cell", "Neuron cell", 
  "Proliferating IC", "Undefined SMC", "Endothelial cell", "Plasma cell"
)
sort(table(seu_comb$new_annotation), decreasing = TRUE)
#seu_comb$cell_type[seu_comb$new_annotation == "Fibroblast 1"] <- "Fibroblast"
seu_comb$new_annotation2 <- seu_comb$new_annotation
levels(seu_comb$new_annotation2)[c(4,7,13,18)] <- "Macro_Mono"
seu_comb$new_annotation2 <- gsub("(.*) SMC", "\\1", seu_comb$new_annotation2)
sce_comb <- as.SingleCellExperiment(seu_comb)
## Get mesenchymal clusters only
seu_smc <- subset(seu_comb, cell_type %in% c("Smooth muscle cell", "Fibroblast"))
#DimPlot_scCustom(seu_smc, group.by = "new_annotation2", label = TRUE) + NoLegend()
## Additionally remove some noisy cells on UMAP
umap1 <- Embeddings(seu_smc, reduction = "umap")
seu_smc <- seu_smc[, Cells(seu_smc)[which(umap1[, "UMAP_1"] > 4)]]
sce_smc <- as.SingleCellExperiment(seu_smc)
p1 <- DimPlot(seu_smc, group.by = "new_annotation2", label = TRUE) + NoLegend()
p2 <- DimPlot(seu_smc, group.by = "sample_id", label = FALSE)
p1 + p2

#### Make pseudo-bulk samples ####
## All dataset
seu_bulk <- AggregateExpression(seu_comb, group.by = "sample_id")$RNA
seu_bulk_cc <- as.data.frame.matrix(table(seu_comb$new_annotation, seu_comb$sample_id))
seu_bulk_cf <- apply(seu_bulk_cc, 2, \(x) x/sum(x))
#rownames(seu_bulk_cf)[16] <- "Proliferating IC"
seu_bulk_clust <- hclust(dist(t(seu_bulk_cf)))
# correlation b/w cell cluster fractions
pheatmap(cor(t(seu_bulk_cf), method = "kendall"), 
         display_numbers = TRUE, name = "Cor coeff", 
         main = "Correlation between cell clusters by cell fraction in human plaque samples")
## Mesenchymal cluster
seu_smc_bulk <- AggregateExpression(seu_smc, group.by = "sample_id")$RNA
seu_smc_cc <- as.data.frame.matrix(table(factor(seu_smc$new_annotation2), 
                                         factor(seu_smc$sample_id)))
seu_smc_cf <- apply(seu_smc_cc, 2, \(x) x/sum(x))
seu_smc_clust <- hclust(dist(t(seu_smc_cf)))
plot(seu_smc_clust)
# correlation b/w cell cluster fractions
pheatmap(cor(t(seu_smc_cf), method = "kendall"), 
         display_numbers = TRUE, name = "Cor coeff", 
         main = "Correlation b/w mesenchymal cell clusters by cell fraction in human plaque samples")
library(ggpubr)
ggscatter(data.frame(t(seu_bulk_cf)), 
          x = "Macrophage", y = "Osteochondrogenic.SMC", 
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          title = "Correlation between cell cluster fractions"
) + stat_cor(method = "pearson", label.x = 0, label.y = 0.07)

#### Init helper functions ####

plotCellTypes <- function(mtx, id.order = NULL, title = NULL) {
  if(is.data.frame(mtx)) { mtx <- mtx }
  else { mtx <- data.frame(mtx) }
  if(is.null(id.order)) { id.order = sort(rownames(mtx), decreasing = TRUE) }
  cfdat <- rownames_to_column(mtx, "sample_id") %>% 
    pivot_longer(-sample_id, names_to = "cell_cluster", values_to = "fraction") %>% 
    mutate(type = "Pseudo-bulk", cell_cluster = factor(cell_cluster),
           sample_id = factor(sample_id, levels = id.order))
  p <- cfdat %>% 
    ggplot(aes(fraction, sample_id, fill = cell_cluster)) + 
    geom_bar(stat="identity", width = 0.8) + 
    scale_fill_manual(values = mypal_mviz) + 
    theme_bw() + 
    labs(x = NULL, y = NULL, title = title)
  return(p)
}


#### CTD with MuSiC ####
library(MuSiC)

### Estimate cell content in pseudo-bulk samples
seu_music <- music_prop(
  bulk.mtx = as.matrix(seu_bulk), sc.sce = sce_comb, 
  clusters = 'new_annotation', samples = 'sample_id'
)
## Compare with pseudo-bulk fractions
seu_cf_real <- data.frame(t(seu_bulk_cf)) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "cell_cluster", values_to = "fraction") %>% 
  mutate(
    type = "Pseudo-bulk", 
    cell_cluster = factor(cell_cluster),
    sample_id = factor(sample_id, 
                       levels = seu_bulk_clust$labels[seu_bulk_clust$order])
    )
p0_real <- seu_cf_real %>% 
  ggplot(aes(fraction, sample_id, fill = cell_cluster)) + 
  geom_bar(stat="identity", width = 0.8) + 
  scale_fill_manual(values = mypal_mviz) + 
  theme_bw() + 
  labs(x = NULL, y = NULL, title = "Real pseudo-bulk proportions")

seu_cf_music <- data.frame(seu_music$Est.prop.weighted) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "cell_cluster", values_to = "fraction") %>% 
  mutate(
    type = "MuSiC predicted proportions", 
    cell_cluster = factor(cell_cluster, levels = levels(seu_cf_real$cell_cluster)),
    sample_id = factor(sample_id, 
                       levels = seu_bulk_clust$labels[seu_bulk_clust$order])
  )
p1_music <- seu_cf_music %>% 
  ggplot(aes(fraction, sample_id, fill = cell_cluster)) + 
  geom_bar(stat="identity", width = 0.8) + 
  scale_fill_manual(values = mypal_mviz) + 
  theme_bw() + theme(axis.text.y.left = element_blank(),
                     axis.ticks.y = element_blank()) + 
  labs(x = NULL, y = NULL, title = "Predicted by MuSiC")

(p0_real + NoLegend()) | p1_music

### Estimate cell content in mesenchymal cell cluster
smc_music <- music_prop(
  bulk.mtx = as.matrix(seu_smc_bulk), sc.sce = sce_smc, 
  clusters = 'new_annotation2', samples = 'sample_id'
)
## Compare with pseudo-bulk fractions
smc_cf_real <- data.frame(t(seu_smc_cf)) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "cell_cluster", values_to = "fraction") %>% 
  mutate(
    type = "Pseudo-bulk", 
    cell_cluster = factor(cell_cluster),
    sample_id = factor(sample_id, 
                       levels = seu_smc_clust$labels[seu_smc_clust$order])
  )
p0_smc_real <- smc_cf_real %>% 
  ggplot(aes(fraction, sample_id, fill = cell_cluster)) + 
  geom_bar(stat="identity", width = 0.8) + 
  scale_fill_manual(values = mypal_mviz) + 
  theme_bw() + 
  labs(x = NULL, y = NULL, title = "Real pseudo-bulk proportions")

smc_cf_music <- data.frame(smc_music$Est.prop.weighted) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "cell_cluster", values_to = "fraction") %>% 
  mutate(
    type = "MuSiC predicted proportions", 
    cell_cluster = factor(cell_cluster, levels = levels(smc_cf_real$cell_cluster)),
    sample_id = factor(sample_id, 
                       levels = seu_smc_clust$labels[seu_smc_clust$order])
  )
p1_smc_music <- smc_cf_music %>% 
  ggplot(aes(fraction, sample_id, fill = cell_cluster)) + 
  geom_bar(stat="identity", width = 0.8) + 
  scale_fill_manual(values = mypal_mviz) + 
  theme_bw() + theme(axis.text.y.left = element_blank(),
                     axis.ticks.y = element_blank()) + 
  labs(x = NULL, y = NULL, title = "MuSiC predicted proportions")

(p0_smc_real + NoLegend()) | p1_smc_music

## Estimate cell content in THOR samples
rc_music <- music_prop(
  bulk.mtx = rc_counts, sc.sce = sce_smc, 
  clusters = 'new_annotation2', samples = 'sample_id', 
  select.ct = unique(sce_smc$new_annotation2)[c(4,2,5,1,7,3)]
  )
## Save the results
write_rds(
  list(scpb_athero = seu_music, scpb_smc = smc_music, thor_rc = rc_music),
  file = "thor_music_deconv.rds"
  )

jitter.fig = Jitter_Est(
  list(data.matrix(rc_music$Est.prop.weighted),
       data.matrix(rc_music$Est.prop.allgene)), 
  method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jitter.fig

rc_music_clust <- hclust(dist(rc_music$Est.prop.weighted))
plotCellTypes(rc_music$Est.prop.weighted, 
              id.order = rc_music_clust$labels[rc_music_clust$order])

#### CTD with BisqueRNA (chosen for publication as the most precise) ####

library(BisqueRNA)
### Make ExpressionSet objects
Idents(seu_comb) <- "new_annotation"
eset_comb <- SeuratToExpressionSet(
  seu_comb, delimiter = "_", position = 1, version = "v3")
pData(eset_comb) <- cbind(pData(eset_comb), 
                          seu_comb@meta.data[colnames(eset_comb), ])
Idents(seu_smc) <- "new_annotation2"
eset_smc <- SeuratToExpressionSet(
  seu_smc, delimiter = "_", position = 1, version = "v3")
pData(eset_smc) <- cbind(pData(eset_smc), 
                          seu_smc@meta.data[colnames(eset_smc), ])
eset_smc$cell_type <- as.character(eset_smc$cell_type)

### Entire dataset
seu_bisque <- ReferenceBasedDecomposition(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_bulk)), 
  sc.eset = eset_comb, use.overlap = FALSE)
(p0_real + NoLegend()) + 
  (plotCellTypes(t(seu_bisque$bulk.props), 
                 id.order = seu_bulk_clust$labels[seu_bulk_clust$order]) + 
     labs(title = "BisqueRNA predicted proportions") + 
     theme(axis.text.y = element_blank()))
### Mesenchymal clusters
smc_bisque <- ReferenceBasedDecomposition(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_smc_bulk)), 
  sc.eset = eset_smc, use.overlap = FALSE
  )
(p0_smc_real + NoLegend()) + 
  (plotCellTypes(t(smc_bisque$bulk.props), 
                 id.order = seu_smc_clust$labels[seu_smc_clust$order]) + 
     labs(title = "BisqueRNA predicted proportions") + 
     theme(axis.text.y = element_blank()))

## Remove undefined SMC cluster
eset_smc2 <- eset_smc[, pData(eset_smc)$new_annotation2 != "Undefined"]

## Estimate cell content in THOR samples
rc_bisque <- ReferenceBasedDecomposition(
  bulk.eset = ExpressionSet(assayData = as.matrix(rc_counts)), 
  sc.eset = eset_smc2, use.overlap = FALSE
  )
  
## Save the results
write_rds(
  list(scpb_athero = seu_bisque, scpb_smc = smc_bisque, thor_rc = rc_bisque),
  file = "thor_bisquerna_deconv.rds"
)

#### CTD with SCDC ####

library(SCDC)

### Entire dataset
seu_scdc1 <- SCDC_prop_subcl_marker(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_bulk)), 
  sc.eset = eset_comb, ct.varname = "new_annotation", 
  fl.varname = "cell_type", sample = "sample_id", 
  ct.fl.sub = unique(eset_comb$cell_type), select.marker = T
)
seu_scdc2 <- SCDC_prop(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_bulk)), 
  sc.eset = eset_comb, ct.varname = "new_annotation", sample = "sample_id", 
  ct.sub = unique(eset_comb$new_annotation)
)

seu_scdc = seu_scdc2
(p0_real + NoLegend()) + 
  (plotCellTypes(seu_scdc$prop.est, 
                 id.order = seu_bulk_clust$labels[seu_bulk_clust$order]) + 
     labs(title = "SCDC predicted proportions") + 
     theme(axis.text.y = element_blank()))

### Mesenchymal clusters
smc_scdc1 <- SCDC_prop_subcl_marker(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_smc_bulk)), 
  sc.eset = eset_smc, ct.varname = "new_annotation2", 
  fl.varname = "cell_type", sample = "sample_id", 
  ct.fl.sub = unique(eset_smc$cell_type), select.marker = T
)
smc_scdc2 <- SCDC_prop(
  bulk.eset = ExpressionSet(assayData = as.matrix(seu_smc_bulk)), 
  sc.eset = eset_smc, ct.varname = "new_annotation2", sample = "sample_id", 
  ct.sub = unique(eset_smc$new_annotation2)
)

smc_scdc = smc_scdc2
(p0_smc_real + NoLegend()) + 
  (plotCellTypes(smc_scdc$prop.est, 
                 id.order = seu_bulk_clust$labels[seu_bulk_clust$order]) + 
     labs(title = "SCDC predicted proportions") + 
     theme(axis.text.y = element_blank()))

write_rds(
  list(scpb_athero_scdc_direct = seu_scdc2, scpb_smc_scdc_direct = smc_scdc2,
       scpb_athero_scdc_tree = seu_scdc1, scpb_smc_scdc_tree = smc_scdc1),
  file = "thor_scdc_deconv.rds"
)


#### CTD with DWLS (-- takes several hours!) ####

library(DWLS)


#### CTD with SQUID (-- based on DWLS, takes several hours!) ####

library(DWLS)
library(SQUID)
source("cell_type_deconvolution/SQUID_helper_functions.R")

### Entire dataset
eset_comb_meta <- data.frame(
  cell.id = colnames(eset_comb),
  sample.id = pData(eset_comb)$sample_id,
  cellType = pData(eset_comb)$new_annotation
  )
seu_squid1 <- SQUID(B = as.matrix(seu_bulk), 
                    scC = exprs(eset_comb), scMeta = eset_comb_meta, 
                    pB = NULL, P = NULL, LeaveOneOut = FALSE)
seu_squid = seu_squid1
(p0_real + NoLegend()) + 
  (plotCellTypes(
    as.matrix(column_to_rownames(
      as.data.frame(pivot_wider(
      seu_squid[, c("sample.id", "cellType", "observed_fraction")], 
      names_from = "cellType", values_from = "observed_fraction")), 
      "sample.id")), 
    id.order = seu_bulk_clust$labels[seu_bulk_clust$order]) + 
     labs(title = "SQUID predicted proportions")) + 
     theme(axis.text.y = element_blank())

### Mesenchymal clusters
eset_smc_meta <- data.frame(
  cell.id = colnames(eset_smc),
  sample.id = pData(eset_smc)$sample_id,
  cellType = pData(eset_smc)$new_annotation2
)
smc_squid1 <- SQUID(B = as.matrix(seu_smc_bulk), 
                    scC = exprs(eset_smc), scMeta = eset_smc_meta, 
                    pB = NULL, P = NULL, LeaveOneOut = FALSE)
smc_squid = smc_squid1
(p0_real + NoLegend()) + 
  (plotCellTypes(smc_squid[], 
                 id.order = seu_bulk_clust$labels[seu_bulk_clust$order]) + 
     labs(title = "SQUID predicted proportions") + 
     theme(axis.text.y = element_blank()))




