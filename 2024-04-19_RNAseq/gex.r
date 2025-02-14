# required for filter_protein_coding
suppressMessages(library(rtracklayer))
suppressMessages(library(plyranges))

##################################################
# Gene expression data loading and preprocessing #
##################################################

filter_genes <- function(gene_matrix){
    # Helper function coarsely removes sets of genes we probably aren't interested in:
    # Not in HUGO, noncoding RNA, antisense, Y paralogues, etc.
    patterns <- "^ENSG|^LINC|-AS\\d+$|^PAR_Y_|-DT$|-OT\\d+$|Y_RNA$|^Metazoa_SRP$|^Vault$|^RNU\\d+-\\d+P$|$RN7SKP\\d+$|^MIR\\d+"
    rows_to_remove <- grep(patterns, rownames(gene_matrix))

    # Subset the matrix to exclude these rows
    filtered_data <- gene_matrix[-rows_to_remove, ]
    return (filtered_data)
}

load_gencode_genes <- function(path_to_gencode) {
    # Load gencode gene annotations from a .gtf file.
    # I'm using the basic primary gene annotation file from GENCODE 47:
    # https://www.gencodegenes.org/human/release_47.html
    g <- rtracklayer::import(path_to_gencode) %>%
        filter(type=='gene') %>%
        filter(gene_type=='protein_coding')
    g$stable_id = sub("\\..*", "", g$gene_id)
    return(g)
}

filter_protein_coding <- function(gene_matrix,path_to_gencode){
    # given a matrix with rownames of the form ENSG00000000003.15_TSPAN6,
    # return a subset consisting of only rows where the stable gene ID ENSG00000000003 
    # is annotated as a protein-coding gene.
    message('filtering for protein-coding genes in gencode 47...')
    g <- load_gencode_genes(path_to_gencode)
    # need to match on stable ENSG IDs
    versioned_ids <- sub("_.*", "", rownames(gene_matrix))
    stable_ids <- sub("\\..*", "", versioned_ids)
    keep <- stable_ids %in% g$stable_id
    filtered_data <- gene_matrix[keep, ]
    return(filtered_data)
}

load_inputs <- function(cts_path,annot_path,path_to_gencode){
    # Returns a list with named elements data$annot and data$cts
    message("loading gene expression data from ", cts_path, " ...")
    cts = as.matrix(read.csv(cts_path,sep='\t',row.names="Gene_or_Transcript_ID",check.names=FALSE))
    annot = read.csv(annot_path,row.names=1) %>%
        filter(!cohort=='PNOC') # drop 1 sample annotated as PNOC

    cts = round(cts) # counts must be integer
    cts <- cts[, rownames(annot)] # only include samples with metadata
    cts = filter_protein_coding(cts,path_to_gencode) # remove noncoding genes
    new_rownames <- sub("^ENSG\\d*\\.\\d*_", "", rownames(cts)) # Remove ENSG prefixes
    rownames(cts) <- new_rownames
    annot$age_at_diagnosis <- (annot$age_at_diagnosis - mean(annot$age_at_diagnosis)) / sd(annot$age_at_diagnosis)
    annot$amplified <- annot$amplicon_class %in% c('ecDNA','chromosomal')
    annot$ecDNA <- annot$amplicon_class == 'ecDNA'
    #annot <- annot[,c("sex","tumor_history","age_at_diagnosis","extent_of_tumor_resection","cancer_type","amplicon_class",'amplified','ecDNA')]

    stopifnot(all(rownames(annot) == colnames(cts)))
    return(list(cts=cts,annot=annot))
}

filter_tumor_types <- function(data,tumor_types){
    # filter inputs to include only the selected tumor types.
    cts = data$cts
    annot = data$annot
    annot <- annot %>% filter(cancer_type %in% tumor_types)
    cts <- cts[, rownames(annot)] # only include samples with metadata
    stopifnot(all(rownames(annot) == colnames(cts)))
    return(list(cts=cts,annot=annot))
}

setup_preprocess_dge <- function(cts,formula,filterbyexp=TRUE){
    dge <- DGEList(cts)
    # remove low-expressed genes
    if (filterbyexp) {
        design <- model.matrix(formula)
        keep <- filterByExpr(dge,design=design)
        paste("Removing",length(keep)-sum(keep),"of",length(keep),"genes with low or invariant gene expression") %>% message
        dge <- dge[keep,]
    }
    # Normalize by library size
    dge <- calcNormFactors(dge,method="TMM")
    return(dge)
}

#########################
# GSEA file conversions #
#########################

write_gene_sets <- function(filepath,set_list){
    # Open a connection to the output file
    file_conn <- file(filepath, open = "wt")
    # Write each set to the file
    for (i in names(set_list)) {
      # Combine the name, description, and vector into a single row
      row <- c(i, "custom_gene_set", set_list[[i]])
      # Write the row to the file, separated by tabs
      writeLines(paste(row, collapse = "\t"), file_conn)
    }
    
    # Close the file connection
    close(file_conn)
}

# Read the tab-separated file
read_gmt = function(file_path){
    data <- readLines(file_path)
    
    # Initialize an empty list to store named vectors
    named_vectors <- list()
    
    # Process each line
    for (line in data) {
      # Split the line by tab
      parts <- strsplit(line, "\t")[[1]]
      
      # The first part is the vector name
      vector_name <- parts[1]
      
      # The third part and onwards are the vector values
      vector_values <- (parts[3:length(parts)])
      
      # Create the named vector and add it to the list
      named_vectors[[vector_name]] <- vector_values
    }
    return(named_vectors)
}


write_gct <- function(es, address) {
  con=file(address)
  open(con, open="w")
  writeLines("#1.2", con)
  ann.col <- ncol(es)
  ann.row <- nrow(es)
  writeLines(sprintf("%s\t%s", ann.row, ann.col), con)
  writeLines(paste0(c("NAME", "Description", colnames(es)), collapse="\t"), con)
  ann.col.table <- cbind(rownames(es),rep(NA,ann.row),es)
  write.table(ann.col.table, file=con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  close(con)
}

write_cls <- function(es, address) {
  con=file(address)
  open(con, open="w")
  writeLines(paste0(c(length(es),2,1), collapse="\t"), con)
  writeLines(paste0("# ",paste0(unique(es), collapse="\t"), collapse=""), con)
  writeLines(paste0(es, collapse="\t"), con)
  close(con)
}

############
# Plotting #
############

write_plot <- function(plt,outfile,width,height){
    pdf.options(encoding='ISOLatin2.enc')
    #pdfName = paste(outfile, ".pdf", sep="")
    pngName = paste(outfile, ".png", sep="")
    svgName = paste(outfile, ".svg", sep = "")
    #ggsave(path="figures", filename=pdfName, device="pdf", width=width, height=height, units='in')
    ggsave(path="figures", device="png", filename=pngName, width=width, height=height, units='in')
    ggsave(path="figures", device="svg", filename=svgName, width=width, height=height, units='in')

}