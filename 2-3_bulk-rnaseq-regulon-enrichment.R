## Load packages
library(tidyverse)
library(RcisTarget)
library(ggrepel)
library(DESeq2)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(yaml) 
library(SummarizedExperiment)
library(DT)
## Set color palettes
mypal_brs1 <- RColorBrewer::brewer.pal(n = 12, name = "Set1")
mypal_brs3 <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
mypal_br_paired <- RColorBrewer::brewer.pal(12, "Paired")
mypal_prot <- khroma::colour("bright")(4)[c(1,4,2)]
### Load data
## Load RNAseq counts
load("/hor_rnaseq_r0r1r2r3_counts_dups_filtered.rda")
## Load RNAseq DEA results
load(file = "smc_kd_combined/thor_rnaseq_r0r1r2r3_dea_results_list.rda")
load(file = "smc_kd_combined/thor_rnaseq_r0r1r2r3_dea_results.rda")
## Load DESeq2 object
rc_dds <- readRDS("d_combined/thor_rnaseq_r0r1r2r3_deseq2_object.rds")
## Load THOR R4 (IFN, TNF) results
r4_results <- readRDS("thor_r4_results_list.rds")
names(r4_results) <- paste0("Baseline_", names(r4_results))
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

### Define thresholds
p_th <- 0.05
lfc_th <- log2(1.5)

##### Perform motif enrichment using RcisTarget #####

## Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_hgnc)
## Load motif rankings
# Downloaded from cistarget DB: https://resources.aertslab.org/cistarget/
motifRankings <- importRankings("RcisTarget/homo_sapiens/hg38/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motifGenes <- colnames(motifRankings)

## Up-regulated genes
rc_deg_up <- c(
  map(rc_dea_results, ~ {
    filter(.x, deseq.padj < p_th, deseq.log2FoldChange >= lfc_th) %>% 
      arrange(GeneName) %>% 
      pull(GeneName) %>% 
      unique(.)
    }),
  map(r4_results[1:2], ~ {
    filter(.x, padj < p_th, log2FoldChange >= lfc_th) %>% 
      arrange(gene_name) %>% 
      pull(gene_name) %>% 
      unique(.)
  })
)
names(rc_deg_up) <- paste0(names(rc_deg_up), "_UP")
## Down-regulated genes
rc_deg_down <- c(
  map(rc_dea_results, ~ {
    filter(.x, deseq.padj < p_th, deseq.log2FoldChange <= -lfc_th) %>% 
      arrange(GeneName) %>% 
      pull(GeneName) %>% 
      unique(.)
  }),
  map(r4_results[1:2], ~ {
    filter(.x, padj < p_th, log2FoldChange <= -lfc_th) %>% 
      arrange(gene_name) %>% 
      pull(gene_name) %>% 
      unique(.)
  })
)
names(rc_deg_down) <- paste0(names(rc_deg_down), "_DOWN")
## Combine UP + DOWN
rc_deg_all <- c(rc_deg_up, rc_deg_down)
rc_deg_all <- rc_deg_all[sort(names(rc_deg_all))]
sapply(rc_deg_all, length)
rc_deg_tab <- imap_dfr(rc_deg_all, ~ data.frame(
  GeneSet = .y, 
  Genes = paste(.x, collapse = "; ")
  )) %>% 
  mutate(tmpcol = GeneSet, .after = GeneSet) %>% 
  separate(tmpcol, sep = "_", into = c("Protocol", "Condition", "Direction"))

# motif enrichment analysis
motif_enrich_genes <- cisTarget(rc_deg_all, motifRankings,
                                motifAnnot = motifAnnotations)
## Get interactive table with logos
motif_enrich_tab <- addLogo(motif_enrich_genes, ) %>% 
  mutate(tmpcol = geneSet, logo = logo, 
         Source = gsub("_.*", "", motif),
         Motif_ID = gsub("^.+_", "", motif),
         .after = "geneSet") %>% 
  separate(tmpcol, sep = "_", into = c("Protocol", "Condition", "Direction")) %>% 
  mutate(Directed_NES = ifelse(Direction == "DOWN", -NES, NES), .after = "NES") 
## Show motif table
DT::datatable(motif_enrich_tab[1:10, -c("enrichedGenes", "TF_lowConf")], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))
## Replace logo at the end of the table
motif_enrich_tab <- motif_enrich_tab %>% 
  dplyr::select(2:ncol(motif_enrich_tab), 1)
## Get annotated highly confident TFs
annotated_tfs <- lapply(split(motif_enrich_tab$TF_highConf,
                              motif_enrich_tab$geneSet),
                        function(x) {
                          genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                          genesSplit <- unique(unlist(strsplit(genes, "; ")))
                          return(genesSplit)
                        })
annotated_tfs_collapsed <- imap_dfr(
  annotated_tfs, 
  ~ data.frame(GeneSet = .y, TF_highConf = paste(.x, collapse = ";"))
  ) %>% 
  mutate(tmpcol = GeneSet, .after = GeneSet) %>% 
  separate(tmpcol, sep = "_", into = c("Protocol", "Condition", "Direction"))

## Get and split highly confident TFs and NES
annotated_tfs_tab <- motif_enrich_tab %>% 
  separate_rows(TF_highConf, sep = ";") %>% 
  distinct_all() %>% 
  mutate(directAnnotation = grepl("directAnnotation", TF_highConf), 
         TF_highConf_tidy = str_trim(gsub("\\(.*\\).*", "", TF_highConf, fixed = FALSE), "both"),
         .after = "TF_highConf") %>% 
  arrange(-directAnnotation) %>% 
  distinct(geneSet, motif, TF_highConf_tidy, .keep_all = TRUE) %>% 
  arrange(geneSet)

## Save results
save(motif_enrich_genes, motif_enrich_tab, annotated_tfs, 
     file = "motif_enrichment.rda")
## Export to Excel
writexl::write_xlsx(list(DEG_list = rc_deg_tab,
                         Motif_enrichment = motif_enrich_tab,
                         Annotated_TF_list = annotated_tfs_collapsed,
                         Annotated_TF_table = annotated_tfs_tab),
                    path = "RcisTarget/motif_enrichment_results.xlsx")


load("motif_enrichment.rda")
## Get and split highly confident TFs and NES
annotated_tfs_tab <- motif_enrich_tab %>% 
  separate_rows(TF_highConf, sep = ";") %>% 
  distinct_all() %>% 
  mutate(directAnnotation = grepl("directAnnotation", TF_highConf), 
         TF_highConf_tidy = str_trim(gsub("\\(.*\\).*", "", TF_highConf, fixed = FALSE), "both"),
         .after = "TF_highConf") %>% 
  arrange(-directAnnotation) %>% 
  distinct(geneSet, motif, TF_highConf_tidy, .keep_all = TRUE) %>% 
  arrange(geneSet)

#top_n = 50
top_n = 50
## Top 50 most frequently regulated regulons in every Protocol
tf_top50 <- annotated_tfs_tab %>% 
  filter(Source == "transfac", 
         !(Condition %in% c("UNTRTD", "TNFvsCTR", "IFNvsCTR"))) %>% 
  group_by(Protocol) %>% 
  dplyr::count(TF_highConf_tidy) %>% 
  arrange(-n) %>% 
  slice_head(n = top_n) %>% 
  ungroup()
## Which of them are dysregulated in all Protocols?
tf_ubiq <- tf_top50 %>% 
  arrange(TF_highConf_tidy) %>% 
  group_by(TF_highConf_tidy) %>% 
  pivot_wider(names_from = Protocol, values_from = n) %>% 
  mutate(Ubiq = !(any(is.na(Baseline), is.na(Cholesterol), is.na(Stretch)))) %>% 
  filter(Ubiq) %>% 
  arrange(TF_highConf_tidy) %>% 
  pull(TF_highConf_tidy)
## Which genes are in the content of these most frequent regulons?
tf_ubiq_df <- annotated_tfs_tab %>% 
  filter(Source %in% c("transfac"), TF_highConf_tidy %in% tf_ubiq) %>% 
  mutate(enrichedGenes = str_split(enrichedGenes, ";"))
write_tsv(tf_ubiq_genes, file = "RcisTarget/regulons_ubiquitously_identified.tsv")
tf_ubiq_genes <- as.data.frame(table(unlist(tf_ubiq_df$enrichedGenes))) %>% 
  arrange(-Freq)
colnames(tf_ubiq_genes) <- c("Gene", "Frequency")
tf_ubiq_genes$Gene <- as.character(tf_ubiq_genes$Gene)
write_tsv(tf_ubiq_genes, file = "RcisTarget/genes_most_frequent_regulon_frequency.tsv")
tf_ubiq_genes_top500 <- map_dfr(tf_ubiq_genes$Gene[1:500], ~ {
  xi <- grep(.x, tf_ubiq_df$enrichedGenes)
  tmp <- tf_ubiq_df[xi, ] %>% 
    mutate(Gene = .x, kdds = paste(Protocol, Condition, sep = "_")) %>% 
    group_by(Gene = Gene, Direction, Motif_ID) %>% 
    reframe(Putative_TFs = paste(unique(TF_highConf_tidy), collapse = ","),,
            Knockdowns = paste(unique(kdds), collapse = ","))
  return(tmp)
})
write_tsv(tf_ubiq_genes_top500, file = "RcisTarget/top500_genes_in_most_frequently_enriched_regulons.tsv")
