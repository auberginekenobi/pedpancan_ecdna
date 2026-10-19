##################################################
# Data processing functions for DE analysis.     #
##################################################

Sys.setenv(LANGUAGE = "en")

library(tidyverse)
library(readxl)
library(edgeR)
library(variancePartition)
library(BiocParallel)
# library(org.Hs.eg.db) library conflict with dplyr::select


gex_metadata <- function(path='../../data/source/gene-expression/rna-metadata.tsv'){
    ## Import RNA-seq library metadata written by notebooks/gex/1-collate-gex-metadata.ipynb, one row per library.
    tbl <- read_tsv(path, col_types = cols(
        rna_id = col_character(),
        cohort = col_character(),
        patient_id = col_character(),
        event_id = col_character(),
        rna_selection = col_character(),
        strandedness = col_character(),
        read_length = col_character(),
        seq_center = col_character(),
        cancer_type = col_character(),
        cancer_subclass = col_character(),
        tumor_history = col_character(),
        sex = col_character(),
        age_at_collection = col_double(),
        age_at_diagnosis = col_double(),
        ecDNA_any = col_logical()
    ))
    return(tbl)
}

harmonize_gene_symbols <- function(symbols){
    ## Map gene symbols to current official symbols (org.Hs.eg.db). Official symbols are kept; previous symbols and aliases that
    ## map to exactly one official symbol are renamed. An alias that maps to several official symbols is renamed to its only
    ## candidate not already present among the other symbols (eg. SEPT2 -> SEPTIN2 rather than SEPTIN6, which SEPT6 maps to).
    ## Unmappable (eg. clone names) or still-ambiguous symbols are kept as they are.
    db <- org.Hs.eg.db::org.Hs.eg.db
    official <- AnnotationDbi::keys(db, keytype = 'SYMBOL')
    aliases <- intersect(setdiff(symbols, official), AnnotationDbi::keys(db, keytype = 'ALIAS'))
    map <- suppressMessages(AnnotationDbi::select(db, keys = aliases, columns = 'SYMBOL', keytype = 'ALIAS')) %>%
        distinct(ALIAS, SYMBOL)
    unique_map <- map %>% group_by(ALIAS) %>% filter(n() == 1) %>% ungroup()
    present <- union(intersect(symbols, official), unique_map$SYMBOL)
    resolved <- map %>% group_by(ALIAS) %>% filter(n() > 1, !SYMBOL %in% present) %>% filter(n() == 1) %>% ungroup()
    map <- bind_rows(unique_map, resolved)
    return(coalesce(map$SYMBOL[match(symbols, map$ALIAS)], symbols))
}

gex_counts <- function(rna_tbl,
                              sj_path='../../data/source/gene-expression/sjcloud/sjc-gex-counts.tsv.gz',
                              cbtn_path='../../data/source/gene-expression/openpedcan/openpedcan-gex-expected-counts.tsv.gz'){
    ## Import raw counts for the RNA-seq libraries in rna_tbl from the SJ and CBTN matrices, joined on current gene symbol.
    ## Gene symbols are first harmonized to current official symbols (harmonize_gene_symbols). Where renaming creates duplicate
    ## symbols within a matrix, keep the row already named with the current symbol; drop any remaining duplicates.
    ## Returns a list:
    ##   counts:    matrix of genes measured in both matrices x rna_tbl$rna_id, in rna_tbl order. CBTN RSEM expected counts are not integers.
    ##   sj_only:   genes in the SJ matrix with no match in the CBTN matrix, dropped
    ##   cbtn_only: genes in the CBTN matrix with no match in the SJ matrix, dropped
    ## rna_tbl: from gex_metadata()
    read_matrix <- function(path, ids, label){
        m <- read_tsv(path, col_select = c(gene, all_of(ids)), col_types = cols(gene = col_character(), .default = col_double())) %>%
            column_to_rownames('gene') %>%
            as.matrix()
        genes <- tibble(original = rownames(m), symbol = harmonize_gene_symbols(rownames(m)), row = seq_len(nrow(m))) %>%
            group_by(symbol) %>%
            filter(n() == 1 | original == symbol) %>%
            filter(n() == 1) %>%
            ungroup()
        message(label, ': renamed ', sum(genes$original != genes$symbol), ' genes to current symbols; dropped ',
                nrow(m) - nrow(genes), ' genes with duplicate symbols after renaming.')
        m <- m[genes$row, ]
        rownames(m) <- genes$symbol
        return(m)
    }
    sj <- read_matrix(sj_path, rna_tbl$rna_id[rna_tbl$cohort == 'SJC'], 'SJ')
    cbtn <- read_matrix(cbtn_path, rna_tbl$rna_id[rna_tbl$cohort == 'CBTN'], 'CBTN')
    shared <- intersect(rownames(sj), rownames(cbtn))
    message(length(shared), ' genes in both matrices; dropped ', nrow(sj) - length(shared), ' SJ-only and ',
            nrow(cbtn) - length(shared), ' CBTN-only genes.')
    counts <- cbind(sj[shared, ], cbtn[shared, ])[, rna_tbl$rna_id]
    return(list(counts = counts,
                sj_only = setdiff(rownames(sj), shared),
                cbtn_only = setdiff(rownames(cbtn), shared)))
}

prepare_gex_covariates <- function(rna_tbl){
    ## Code the model covariates for the RNA-seq libraries (see notebooks/gex/model-design.md).
    ## Drops libraries without tumor_history. Stops if any covariate has a value outside the expected levels.
    tbl <- rna_tbl %>%
        filter(!is.na(tumor_history)) %>%
        mutate(
            molecular_subtype = ifelse(is.na(cancer_subclass), as.character(cancer_type), paste(cancer_type, cancer_subclass, sep = ".")),
            sex = factor(replace_na(sex, 'unknown'), levels = c('Female','Male','unknown')),
            tumor_history = factor(if_else(tumor_history %in% c('Progressive','Recurrence','Relapse','Second Malignancy'), 'Relapse', tumor_history),
                                   levels = c('Diagnosis','Relapse','Autopsy','Metastasis')),
            cohort = factor(cohort, levels = c('CBTN','SJC')),
            rna_selection = factor(rna_selection, levels = c('Total','PolyA')),
            read_length_sj = factor(if_else(cohort == 'SJC', read_length, '126'), levels = c('126','101','NEC')),
            seq_center_cbtn = factor(if_else(cohort == 'CBTN', replace_na(seq_center, 'NotAvailable'), 'NantOmics'),
                                     levels = c('NantOmics','BGI@CHOP','Ashion','NotAvailable')))
    covariates <- c('molecular_subtype','sex','tumor_history','cohort','rna_selection','read_length_sj','seq_center_cbtn','patient_id','event_id')
    missing <- colSums(is.na(tbl[covariates]))
    if (any(missing > 0)) stop('Covariates with missing or unexpected values: ', paste(names(missing)[missing > 0], collapse = ', '))
    return(tbl)
}

generate_dream_cache <- function(rna,cache_path){
    ## Compute voom weights for all expressed genes with voomWithDreamWeights, and save to cache_path. Takes over an hour.
    ## The random effects must match those in test_gene.
    ## rna: prepare_gex_covariates output
    counts <- gex_counts(rna)$counts
    dge <- edgeR::DGEList(counts)
    keep <- rowSums(cpm(dge) >= 1) >= 10 # expressed (>= 1 CPM) in at least 10 libraries
    dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
    formula <- ~ molecular_subtype + sex + tumor_history + cohort + rna_selection + read_length_sj + seq_center_cbtn + (1 | patient_id) + (1 | event_id)
    vobj <- voomWithDreamWeights(dge, formula, column_to_rownames(rna, 'rna_id'), BPPARAM = MulticoreParam(8))
    saveRDS(vobj, cache_path)
    return(vobj)
}

annotate_exposure <- function(gene_name, rna_tbl, biosamples_tbl, genes_tbl){
    ## Annotate each RNA-seq library with the amplification status of one gene across the WGS biosamples of its collection event.
    ## Status of each WGS biosample:
    ##   ecDNA:            a full (untruncated) copy of the gene is on ecDNA
    ##   intrachromosomal: otherwise, a full copy is on a BFB, linear or complex non-cyclic amplicon
    ##   excluded:         the gene is amplified only in unknown features, or only truncated
    ##   none:             the gene is not amplified
    ## Exposure of each event (and so of its RNA-seq libraries) is the highest-priority status among its WGS biosamples:
    ##   ecDNA > intrachromosomal > excluded > none
    ## exposure_concordant: all WGS biosamples from the event have the same status.
    ## Excluded libraries are dropped from tests; see notebooks/gex/model-design.md for the amplified / ecDNA indicators.
    ## rna_tbl: rna-metadata.tsv; biosamples_tbl: S.T.2; genes_tbl: S.T.5
    if (!gene_name %in% genes_tbl$gene) {
        stop(paste0(gene_name, " is not amplified in any biosample in S.T.5; check the gene name."))
    }
    amp <- genes_tbl %>%
        filter(gene == gene_name) %>%
        mutate(full = is.na(truncated)) %>%
        group_by(sample_name) %>%
        summarise(status = case_when(
            any(full & startsWith(feature, 'ecDNA')) ~ 'ecDNA',
            any(full & !startsWith(feature, 'unknown')) ~ 'intrachromosomal',
            TRUE ~ 'excluded'))
    priority <- c('ecDNA','intrachromosomal','excluded','none')
    events <- biosamples_tbl %>%
        mutate(event_id = coalesce(external_sample_id, biosample_id)) %>%
        left_join(amp, by = c('biosample_id' = 'sample_name')) %>%
        mutate(status = replace_na(status, 'none')) %>%
        group_by(event_id) %>%
        summarise(exposure = priority[min(match(status, priority))],
                  exposure_concordant = n_distinct(status) == 1)
    rna_tbl %>%
        left_join(events, by = 'event_id') %>%
        mutate(exposure = factor(exposure, levels = c('none','ecDNA','intrachromosomal','excluded')))
}

fit_gene <- function(gene_name, vobj, rna_tbl, biosamples_tbl, genes_tbl){
    ## Fit the dream model for one gene by its amplification status (see notebooks/gex/model-design.md). Libraries with excluded
    ## amplifications are dropped. Summarize the result with test_gene and gene_case_excess.
    ## Returns a list:
    ##   gene:    gene_name
    ##   fit:     dream fit; NULL if no library has the gene amplified
    ##   data:    model data, one row per library, with rna_id as row names
    ##   log_cpm: voom log-CPM of the gene in each library of data
    ##   n:       number of libraries by exposure (none, ecDNA, intrachromosomal, excluded), before dropping excluded libraries
    ##   coefs:   fit coefficient for each estimable contrast (ecDNA_vs_none, intra_vs_none, ecDNA_vs_intra)
    ## vobj: voomWithDreamWeights output; rna_tbl: prepare_gex_covariates output; biosamples_tbl: S.T.2; genes_tbl: S.T.5
    if (!gene_name %in% rownames(vobj)) stop(gene_name, ' is not in the voom object (not measured in both cohorts, or filtered as lowly expressed).')
    x <- annotate_exposure(gene_name, rna_tbl, biosamples_tbl, genes_tbl)
    n <- table(x$exposure)
    x <- x %>%
        filter(exposure != 'excluded') %>%
        mutate(amplified = exposure %in% c('ecDNA','intrachromosomal'),
               ecDNA = exposure == 'ecDNA') %>%
        column_to_rownames('rna_id')
    # NB the cached voom weights may have been computed with (1 | patient_id) only; recompute the cache with generate_dream_cache before final results
    covariates <- 'molecular_subtype + sex + tumor_history + cohort + rna_selection + read_length_sj + seq_center_cbtn + (1 | patient_id) + (1 | event_id)'
    fit <- NULL
    coefs <- c()
    if (n[['ecDNA']] > 0 && n[['intrachromosomal']] > 0) {
        form <- as.formula(paste('~ amplified + ecDNA +', covariates))
        # ecDNA vs none is the sum of the two coefficients; makeContrastsDream warns that it does not sum to 0, which is expected
        L <- suppressWarnings(makeContrastsDream(form, x, contrasts = c(ecDNA_vs_none = 'amplifiedTRUE + ecDNATRUE')))
        fit <- dream(vobj[gene_name, rownames(x)], form, x, L = L, ddf = 'Satterthwaite')
        coefs <- c(ecDNA_vs_none = 'ecDNA_vs_none', intra_vs_none = 'amplifiedTRUE', ecDNA_vs_intra = 'ecDNATRUE')
    } else if (n[['ecDNA']] + n[['intrachromosomal']] > 0) {
        form <- as.formula(paste('~ amplified +', covariates))
        fit <- dream(vobj[gene_name, rownames(x)], form, x, ddf = 'Satterthwaite')
        coefs <- if (n[['ecDNA']] > 0) c(ecDNA_vs_none = 'amplifiedTRUE') else c(intra_vs_none = 'amplifiedTRUE')
    }
    return(list(gene = gene_name, fit = fit, data = x, log_cpm = vobj$E[gene_name, rownames(x)], n = n, coefs = coefs))
}

fit_genes <- function(gene_names, vobj, rna_tbl, biosamples_tbl, genes_tbl){
    ## Run fit_gene on each gene in gene_names, with a progress bar; return a list of fits named by gene.
    ## Genes that cannot be fitted (eg. not in the voom object, or not amplified in S.T.5) get list(gene, error) with the reason.
    gene_names %>%
        set_names() %>%
        map(function(gene_name) {
            tryCatch(fit_gene(gene_name, vobj, rna_tbl, biosamples_tbl, genes_tbl),
                     error = function(e) list(gene = gene_name, error = conditionMessage(e)))
        }, .progress = TRUE)
}

test_gene <- function(m){
    ## Contrasts for one gene from fit_gene, as a one-row tibble, each with logFC, SE, df, t and p (Satterthwaite):
    ##   ecDNA_vs_none, intra_vs_none, ecDNA_vs_intra
    ## Contrasts that cannot be estimated (no ecDNA or no intrachromosomal cases) are NA; genes that could not be fitted give
    ## NA results and the reason in `error`.
    count <- function(exposure) if (is.null(m$n)) NA_integer_ else m$n[[exposure]]
    result <- tibble(gene = m$gene,
                     n_none = count('none'), n_ecDNA = count('ecDNA'), n_intra = count('intrachromosomal'), n_excluded = count('excluded'))
    for (contrast in c('ecDNA_vs_none','intra_vs_none','ecDNA_vs_intra')) {
        coef <- m$coefs[contrast]
        estimable <- !is.null(coef) && !is.na(coef)
        result[[paste0(contrast, '_logFC')]] <- if (estimable) m$fit$coefficients[1, coef] else NA_real_
        result[[paste0(contrast, '_SE')]] <- if (estimable) m$fit$stdev.unscaled[1, coef] * m$fit$sigma[1] else NA_real_
        result[[paste0(contrast, '_df')]] <- if (estimable) m$fit$df.residual[1, coef] else NA_real_
        result[[paste0(contrast, '_t')]] <- if (estimable) m$fit$t[1, coef] else NA_real_
        result[[paste0(contrast, '_p')]] <- if (estimable) m$fit$p.value[1, coef] else NA_real_
    }
    result$error <- if (is.null(m$error)) NA_character_ else m$error
    return(result)
}

gene_case_excess <- function(m){
    ## Expression of one gene in each library relative to the model's prediction for an unamplified library with the same covariates:
    ## excess = observed log-CPM - fixed-effect prediction with amplified==FALSE, ecDNA==FALSE (fit_gene; excluded libraries dropped).
    ## For amplified libraries this is the per-case excess over the unamplified baseline; unamplified libraries scatter around 0.
    ## Random effects are not subtracted, so the excess includes patient- and event-level deviations.
    ## Returns no rows for genes that could not be fitted or are not amplified in any library.
    if (is.null(m$fit)) return(tibble())
    X0 <- m$fit$design
    X0[, intersect(c('amplifiedTRUE','ecDNATRUE'), colnames(X0))] <- 0
    m$data %>%
        rownames_to_column('rna_id') %>%
        transmute(gene = m$gene, rna_id, patient_id, cohort, molecular_subtype, exposure, exposure_concordant,
                  log_cpm = m$log_cpm,
                  predicted_none = drop(X0 %*% m$fit$coefficients[1, colnames(X0)]),
                  excess = log_cpm - predicted_none)
}
