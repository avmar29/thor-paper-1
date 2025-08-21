## Define function for searching for DEGs
getDEGFromCounts <- function(count_table, pheno_table, group_column = NULL, 
                             filter_method = "edger", min_count = 10, min_prop = 0.7,
                             method = "deseq2", alpha_th = 0.1, p_adjust_method = "BH", 
                             collect_removed = FALSE, keep_object = FALSE, verbose = FALSE) {
  ## set a group variable
  if (is.null(group_column)) { 
    group_column <- colnames(pheno_table)[ncol(pheno_table)]
  }
  group_var <- factor(pheno_table[, group_column])
  pheno_table[, group_column] <- group_var
  ## set group levels and comparison combinations
  group_lvls <- levels(group_var)
  comp_combs <- combn(length(group_lvls), 2)
  comp_combs <- comp_combs[2:1, ]
  ## hotfix !!!
  if(!any(class(comp_combs) == "matrix")) { comp_combs <- matrix(comp_combs) }
  comp_contr <- paste0(group_lvls[comp_combs[1,]], "-", group_lvls[comp_combs[2,]])
  ## set and check count table
  count_table <- round(count_table)
  ct_class <- apply(count_table, 2, class)
  if (!all(ct_class == "numeric" | ct_class == "integer")) {
    stop("Count table contain non-integer values!")
  }
  
  ## remove lowly expressed genes
  cat("1. Filtering features with low expression...\n")
  require(edgeR)
  require(data.table)
  cpm_table <- cpm(count_table)
  cpm_th <- round(mean(min_count*1e6/colSums(count_table)), 1)
  cat("CPM threshold to keep genes with", min_count, "counts at least for libraries with average size of", 
      round(mean(colSums(count_table)/1e6)), "million reads is", cpm_th, "\n")
  if (filter_method == "experimental") {
    min_n <- min(table(group_var))
    keep_table <- data.table(t(cpm_table))[, lapply(.SD, function(x) sum(x >= cpm_th) >= floor(min_prop * min_n)), by = group_var][, -1]
    keep <- apply(keep_table, 2, any)
    cat("Used experimental function with the following assumption:\n",
        "* CPM must not be less than", cpm_th, 
        "in as minimum as", floor(min_prop * min_n), "samples in any compared group. \n")
  } else if (filter_method == "edger") {
    cat("Used function 'filterByExpr' from the 'edgeR' package (with parameter min.count =", 
        min_count, "and others by default).  \n")
    keep <- filterByExpr(y = count_table, group = group_var, min.count = min_count)
  } else {
    keep <- rep(TRUE, nrow(count_table))
  }
  cat("Out of", nrow(count_table), "features:", 
      table(keep)[1], "were filtered out and", 
      table(keep)[2], "remained.\n")
  count_table <- count_table[keep, ]
  cpm_table <- cpm_table[keep, ]
  cpm_mean <- t(apply(cpm_table, 1, function(x) tapply(x, group_var, mean)))
  cpm_sd <- t(apply(cpm_table, 1, function(x) tapply(x, group_var, sd)))
  ## Find DEGs
  cat("2. Searching for differentially expressed genes...\n")
  output <- list()
  
  ## DESeq2
  if (method == "deseq2") {
    cat("Applied method: Fitting to Negative Binomial distribution \nand Wald test for significance of estimated coefficients \n(default workflow using 'DESeq2' package) \n")
    require(DESeq2)
    deseq_ds <- DESeqDataSetFromMatrix(countData = count_table,
                                       colData = DataFrame(pheno_table),
                                       design = as.formula(paste("~", group_column)))
    deseq_ds <- estimateSizeFactors(deseq_ds)
    deseq_ds <- DESeq(deseq_ds)
    for (i in 1:ncol(comp_combs)) {
      res <- results(deseq_ds, alpha = alpha_th, pAdjustMethod = p_adjust_method,
                     contrast = c(group_column, group_lvls[comp_combs[, i]]),
                     independentFiltering = FALSE)
      res <- as.data.frame(res[order(res$pvalue),])
      res <- data.frame(CPM.Mean = cpm_mean[rownames(res), group_lvls[comp_combs[, i]]],
                        CPM.SD = cpm_sd[rownames(res), group_lvls[comp_combs[, i]]],
                        res)
      #output[[i]] <- as.data.frame(res[order(res$pvalue),])
      output[[i]] <- as.data.frame(res[order(rownames(res)),])
    }; names(output) <- comp_contr
    ## get some stats
    if (verbose) {
      cat("Number of DEGs with unadjusted p-value less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$pvalue < alpha_th))[2,])
      cat("Number of DEGs with p-value adjusted by", p_adjust_method, "method less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$padj < alpha_th))[2,])
    }
    if (keep_object) {
      output$object <- deseq_ds
    }
  } 
  
  ## edgeR
  else if (method == "edger") {
    cat("Applied method: Exact Tests for differences between two groups of Negative-Binomial counts \nproposed by Robinson and Smyth (2008) \n(default workflow using 'edgeR' package) \n")
    require(edgeR)
    edger_ds <- DGEList(count_table, group = group_var)
    edger_ds <- calcNormFactors(edger_ds)
    edger_ds <- estimateDisp(edger_ds)
    for (i in 1:ncol(comp_combs)) {
      tmp <- topTags(exactTest(edger_ds, pair = rev(comp_combs[, i])), n = Inf)$table
      output[[i]] <- tmp[order(rownames(tmp)), ]
    }; names(output) <- comp_contr
    ## get some stats
    if (verbose) {
      cat("Number of DEGs with unadjusted p-value less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$PValue < alpha_th))[2,])
      cat("Number of DEGs with p-value adjusted by", p_adjust_method, "method less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$FDR < alpha_th))[2,])
    }
    if (keep_object) {
      output$object <- edger_ds
    }
  } 
  
  ## voom + limma
  else if (method == "limma") {
    cat("Applied method: Voom transformation, calculation of variance weights \nfollowed by GLM fitting with Empirical Bayes smoothing of standard errors \n('voom + limma' approach) \n")
    require(limma)
    edger_ds <- DGEList(count_table, group = group_var)
    edger_ds <- calcNormFactors(edger_ds)
    edger_ds <- estimateDisp(edger_ds)
    mm <- model.matrix(~0 + group_var)
    colnames(mm) <- group_lvls
    y <- voom(edger_ds, mm, plot = F)
    fit <- lmFit(y, mm)
    contr <- makeContrasts(contrasts = comp_contr, levels = group_lvls)
    eb_fit <- eBayes(contrasts.fit(fit, contr))
    #results <- decideTests(eb_fit, p.value = alpha_th)
    #summary(results)
    for (i in 1:ncol(comp_combs)) {
      #output[[i]] <- topTable(eb_fit, coef = i, sort.by = "P", n = Inf)
      tmp <- topTable(eb_fit, coef = i, sort.by = "none", n = Inf)
      output[[i]] <- tmp[order(rownames(tmp)), ]
    }; names(output) <- comp_contr
    ## get some stats
    if (verbose) {
      cat("Number of DEGs with unadjusted p-value less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$P.Value < alpha_th))[2,])
      cat("Number of DEGs with p-value adjusted by", p_adjust_method, "method less than", alpha_th, "\n")
      print(sapply(output, function(x) table(x$adj.P.Val < alpha_th))[2,])
    }
    if (keep_object) {
      output$object <- edger_ds
    }
  }
  
  ## include filtered out features
  if (collect_removed) {
    output$removed_features <- count_table[!keep, ]
  }
  return(output)
}






