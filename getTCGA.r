library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

getTCGA <- function(disease = "GBM", clinical = TRUE, cvars = NULL,
                    workflows = c("HTSeq - FPKM", "STAR - Counts")) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Please install the TCGAbiolinks package.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install the SummarizedExperiment package.")
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Please install the biomaRt package.")
  }
  
  project <- paste0("TCGA-", toupper(disease))
  message("Searching for RNASeq2-like data for ", project, "...")
  
  for (workflow in workflows) {
    message("Trying workflow: ", workflow)
    query <- tryCatch({
      TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = workflow
      )
    }, error = function(e) NULL)
    
    if (!is.null(query) && length(query$results[[1]]) > 0) {
      message("Found data using workflow: ", workflow)
      TCGAbiolinks::GDCdownload(query)
      exp_data <- TCGAbiolinks::GDCprepare(query)
      
      # Map Ensembl IDs to gene symbols
      ensembl_ids <- gsub("\\..*", "", rownames(exp_data))  # remove version numbers
      mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      gene_map <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart
      )
      gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
      gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
      
      rownames(exp_data) <- ensembl_ids
      matched <- match(rownames(exp_data), gene_map$ensembl_gene_id)
      gene_symbols <- gene_map$hgnc_symbol[matched]
      valid <- !is.na(gene_symbols) & gene_symbols != ""
      
      exp_data <- exp_data[valid, ]
      gene_symbols <- gene_symbols[valid]
      
      # Select assay: FPKM preferred, fallback to TPM
      available_assays <- SummarizedExperiment::assayNames(exp_data)
      if ("fpkm_unstrand" %in% available_assays) {
        assay.type <- "fpkm_unstrand"
        message("Using assay: fpkm_unstrand")
      } else if ("tpm_unstrand" %in% available_assays) {
        assay.type <- "tpm_unstrand"
        message("FPKM not available. Falling back to assay: tpm_unstrand")
      } else {
        stop("Neither FPKM nor TPM assays are available. Available assays: ",
             paste(available_assays, collapse = ", "))
      }
      
      exp_matrix <- SummarizedExperiment::assay(exp_data, assay.type)
      rownames(exp_matrix) <- gene_symbols
      
      if (!clinical) {
        return(list(dat = exp_data, workflow = workflow, assay.type = assay.type))
      }
      
      message("Retrieving clinical data...")
      cli <- TCGAbiolinks::GDCquery_clinic(project = project, type = "clinical")
      
      exp_barcodes <- substr(colnames(exp_matrix), 1, 12)
      cli_barcodes <- cli$submitter_id
      common <- intersect(exp_barcodes, cli_barcodes)
      
      exp_matrix <- exp_matrix[, exp_barcodes %in% common]
      cli <- cli[cli$submitter_id %in% common, ]
      
      # Ensure matching order
      cli <- cli[match(substr(colnames(exp_matrix), 1, 12), cli$submitter_id), ]
      stopifnot(all(substr(colnames(exp_matrix), 1, 12) == cli$submitter_id))
      
      # Aggregate by gene symbol
      agg_expr <- rowsum(exp_matrix, group = rownames(exp_matrix))
      
      # Add OS_time and OS_event
      cli$OS_time <- as.numeric(ifelse(
        is.na(cli$days_to_death),
        cli$days_to_last_follow_up,
        cli$days_to_death
      ))
      cli$OS_event <- ifelse(cli$vital_status == "Dead", 1, 0)
      
      # Handle cvars
      if (!is.null(cvars)) {
        available_vars <- colnames(cli)
        valid_cvars <- cvars[cvars %in% available_vars]
        missing_cvars <- setdiff(cvars, valid_cvars)
        
        if (length(valid_cvars) > 0) {
          cli <- cli[, unique(c(valid_cvars, "OS_time", "OS_event")), drop = FALSE]
          if (length(missing_cvars) > 0) {
            message("Note: The following clinical variables were not found and will be skipped: ",
                    paste(missing_cvars, collapse = ", "))
          }
        } else {
          cli <- cli[, c("OS_time", "OS_event"), drop = FALSE]
        }
      } else {
        cli <- cli[, c("OS_time", "OS_event"), drop = FALSE]
      }
      
      # Create merged data frame
      merged.dat <- as.data.frame(t(agg_expr))
      merged.dat <- cbind(merged.dat, cli)
      
      sample_ids <- substr(colnames(exp_matrix), 1, 12)
      if (any(duplicated(sample_ids))) {
        message("Duplicate sample IDs found. Making row names unique.")
        sample_ids <- make.unique(sample_ids)
      }
      rownames(merged.dat) <- sample_ids
      
      return(list(dat = exp_data, clinical = cli, merged.dat = merged.dat,
                  workflow = workflow, assay.type = assay.type))
    }
  }
  
  stop("No expression data found for the specified project using the provided workflows.")
}
