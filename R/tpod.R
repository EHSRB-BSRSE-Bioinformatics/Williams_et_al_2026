# Mock-dose benchmark concentration analysis.

tpod_designs <- function() {
  doses <- list(
    `4` = c(0, 0.01, 1, 100),
    `6` = c(0, 0.01, 0.1, 1, 10, 100),
    `8` = c(0, 0.01, 0.05, 0.1, 0.5, 1, 10, 100),
    `12` = c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  )
  group_counts <- c(6L, 4L, 8L, 12L)
  designs <- list()
  # The design order determines the seeded draws in the collaborator's script.
  for (total in c(24L, 48L, 72L)) {
    for (groups in group_counts) {
      replicates <- total / groups
      design <- paste0("d", groups, "n", replicates)
      designs[[design]] <- rep(doses[[as.character(groups)]], each = replicates)
    }
  }
  # The supplied d6n8 code reassigns positions 31 and 32 from dose 1 to 10.
  designs$d6n8[31:32] <- 10
  designs
}

prepare_tpod_expression <- function(counts, metadata) {
  required <- c("Sample_ID", "chemical", "cell_line")
  if (!all(required %in% names(metadata)) ||
      anyNA(metadata$Sample_ID) || anyDuplicated(metadata$Sample_ID)) {
    stop("Metadata must have unique Sample_ID values and chemical and cell_line columns.")
  }
  library_sizes <- colSums(counts)
  metadata_index <- match(colnames(counts), metadata$Sample_ID)
  controls <- which(metadata$chemical[metadata_index] == "DMSO" &
    library_sizes >= 0.1 * median(library_sizes))
  adjusted_counts <- 1e6 * (counts[, controls, drop = FALSE] + 0.5)
  expression <- log2(sweep(adjusted_counts, 2L, library_sizes[controls] + 1, "/"))
  expression[expression < log2(5)] <- log2(5)
  list(
    expression = expression,
    cell_line = metadata$cell_line[metadata_index[controls]]
  )
}

make_tpod_assignments <- function(sample_ids, designs, n_trials, seed = 90214L) {
  schedules <- tpod_designs()
  if (!length(designs) || !all(designs %in% names(schedules)) ||
      anyDuplicated(designs)) {
    stop("Select unique designs from ", paste(names(schedules), collapse = ", "), ".")
  }
  if (n_trials < 1L) stop("The number of trials must be positive.")
  if (length(sample_ids) < max(lengths(schedules[designs]))) {
    stop("Too few eligible DMSO samples for the selected designs.")
  }

  set.seed(seed)
  assignments <- lapply(names(schedules), function(design) {
    dose <- schedules[[design]]
    # Draw for all designs so a subset uses the same draws as a full run.
    trials <- lapply(seq_len(n_trials), function(trial) {
      indices <- head(order(runif(length(sample_ids))), length(dose))
      data.frame(trial = trial, sample_id = sample_ids[indices], dose = dose)
    })
    do.call(rbind, trials)
  })
  names(assignments) <- names(schedules)
  assignments[designs]
}

write_tpod_inputs <- function(expression, assignments, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (trial in unique(assignments$trial)) {
    samples <- assignments[assignments$trial == trial, ]
    values <- rbind(samples$dose, expression[, samples$sample_id, drop = FALSE])
    input <- data.frame(ProbeID = c("Dose", rownames(expression)),
                        values, check.names = FALSE)
    write.table(input, file.path(output_dir, paste0("Exp ", trial, ".txt")),
                sep = "\t", row.names = FALSE)
  }
}

tpod_export_columns <- c(
  analysis = "^Analysis$", probe = "^Probe ID$",
  bmd = "^Best BMD$", bmdl = "^Best BMDL$", bmdu = "^Best BMDU$",
  rsq = "^Best RSquared$", fit_p = "^Best fitPValue$",
  trend_p = "Prefilter.*P-Value", fold_change = "Max Fold Change Absolute Value"
)

tpod_pathway_collections <- c("KEGG", "Reactome", "Hallmark", "GO")

summarize_tpod_fpr <- function(export, annotation, pathways = NULL,
                               expected_trials = 100L) {
  if (!nrow(export)) stop("The BMDExpress export has no data rows.")
  column_names <- trimws(names(export))
  columns <- lapply(tpod_export_columns, function(pattern) {
    matches <- grep(pattern, column_names, ignore.case = TRUE)
    if (length(matches) > 1L) stop("Ambiguous export column matching ", pattern, ".")
    if (length(matches)) export[[matches]]
  })
  missing <- names(columns)[lengths(columns) == 0L]
  optional <- c("trend_p", "fold_change")
  required_missing <- setdiff(missing, optional)
  if (length(required_missing)) {
    stop("Missing required BMDExpress fields: ", paste(required_missing, collapse = ", "))
  }
  if (length(missing)) {
    warning("Skipping prefilters for absent export fields: ", paste(missing, collapse = ", "))
  }

  analysis <- trimws(columns$analysis)
  valid_analysis <- !is.na(analysis) & nzchar(analysis) & analysis != "NA"
  trial_rows <- split(which(valid_analysis), analysis[valid_analysis])
  if (!length(trial_rows) || length(trial_rows) != expected_trials) {
    stop("Expected ", expected_trials, " analysis sets, found ", length(trial_rows), ".")
  }
  names(annotation) <- trimws(names(annotation))
  if (!all(c("Probe_Name", "Enzrez_ID") %in% names(annotation))) {
    stop("Annotation must contain Probe_Name and Enzrez_ID columns.")
  }
  entrez <- trimws(as.character(annotation$Enzrez_ID))
  annotated <- !is.na(entrez) & nzchar(entrez) & entrez != "NA"
  probe_ids <- split(entrez[annotated], annotation$Probe_Name[annotated], drop = TRUE)
  columns$probe <- as.character(columns$probe)
  reported_probes <- columns$probe[valid_analysis]
  reported_probes <- reported_probes[!is.na(reported_probes) &
    nzchar(reported_probes) & reported_probes != "NA"]
  if (length(reported_probes) && !any(reported_probes %in% names(probe_ids))) {
    stop("No export probes matched the Entrez annotation.")
  }

  collections <- list()
  if (!is.null(pathways)) {
    fields <- c("collection", "pathway", "entrez_id")
    if (!all(fields %in% names(pathways)) ||
        anyNA(pathways[fields]) ||
        any(!nzchar(as.matrix(pathways[fields]))) ||
        !setequal(unique(pathways$collection), tpod_pathway_collections)) {
      stop("Pathway table must have collection, pathway and entrez_id for KEGG, Reactome, Hallmark and GO.")
    }
    collections <- lapply(tpod_pathway_collections, function(collection) {
      members <- pathways[pathways$collection == collection, ]
      sets <- lapply(split(as.character(members$entrez_id), members$pathway), unique)
      sets[lengths(sets) >= 30L]
    })
    names(collections) <- paste0("Pathway-", tpod_pathway_collections)
  }

  numeric_fields <- setdiff(names(columns), c("analysis", "probe", missing))
  columns[numeric_fields] <- lapply(columns[numeric_fields], function(values) {
    suppressWarnings(as.numeric(as.character(values)))
  })
  base_filter <- with(columns,
    bmd > 0 & bmd < 100 & fit_p > 0 & bmdl > 0 &
      bmd / bmdl < 20 & bmdu / bmdl < 40
  )
  # This grid order follows the original nested loops.
  grid <- expand.grid(rsq = c(0.6, 0), p = c(0.05, 0.01), fc = c(1.5, 2))
  endpoints <- c("Gene (rank 25)", names(collections))
  summaries <- lapply(seq_len(nrow(grid)), function(i) {
    setting <- grid[i, ]
    keep <- base_filter & columns$rsq >= setting$rsq
    if (!is.null(columns$trend_p)) keep <- keep & columns$trend_p < setting$p
    if (!is.null(columns$fold_change)) keep <- keep & columns$fold_change >= setting$fc

    calls <- vapply(trial_rows, function(rows) {
      probes <- unique(columns$probe[rows[which(keep[rows])]])
      ids <- unlist(probe_ids[probes], use.names = FALSE)
      ids <- unique(trimws(unlist(strsplit(as.character(ids), "[|;,]"))))
      ids <- ids[nzchar(ids) & ids != "NA"]
      pathway_calls <- vapply(collections, function(sets) {
        if (length(ids) < 3L) return(FALSE)
        any(vapply(sets, function(genes) {
          overlap <- sum(ids %in% genes)
          overlap >= 3L && overlap / length(genes) >= 0.05
        }, logical(1L)))
      }, logical(1L))
      c(length(ids) >= 25L, pathway_calls)
    }, logical(length(endpoints)))
    n_total <- length(trial_rows)
    n_positive <- rowSums(matrix(calls, nrow = length(endpoints)))
    rsq_on <- setting$rsq > 0
    data.frame(
      R_Squared_Filter = if (rsq_on) "R-Sq > 0.6" else "No R-Sq",
      Fold_Change_Cutoff = setting$fc,
      P_Value_Cutoff = setting$p,
      tPoD = endpoints,
      n_false_positives = n_positive,
      n_total = n_total,
      FPR = n_positive / n_total,
      Param = paste(if (rsq_on) "R2>=0.6" else "No R2",
                    paste0("FC", setting$fc), paste0("p", setting$p), sep = " | ")
    )
  })
  do.call(rbind, summaries)
}
