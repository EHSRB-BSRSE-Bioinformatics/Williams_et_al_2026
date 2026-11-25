source("R/revision/reactome_ora.R")
args <- commandArgs(TRUE)
mode <- if ("--full" %in% args) "full" else if ("--pilot-parallel" %in% args) "pilot_parallel_org3.17" else "pilot_org3.17"
cores <- as.integer(Sys.getenv("REACTOME_ORA_CORES", "1"))
stopifnot(!is.na(cores), cores >= 1L)
if (mode == "full") suppressPackageStartupMessages(library(DESeq2))
root <- "data/large/reactome_ora_v1"
out <- file.path(root, mode)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
annotation_file <- file.path(root,"annotation/frozen_annotation_org3.17.rds")
a <- readRDS(annotation_file)
compiled <- compile_reactome(a$pathways)
metadata <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
identity <- trial_identity(names(metadata))
reference <- readRDS("data/reference/deg_reference.rds")
annotation_md5 <- unname(tools::md5sum(annotation_file))
source_md5 <- unname(tools::md5sum(c("R/revision/reactome_ora.R", "scripts/revision/run_reactome_ora.R")))
validate_checkpoint <- function(path, expected) {
  x <- readRDS(path)
  stopifnot(identical(x$schema,"reactome_ora_v1"), identical(names(x$rows),expected),
            identical(x$annotation_md5,annotation_md5), identical(x$source_md5,source_md5),
            !any(vapply(x$rows,function(z) z$result$status == "failed",logical(1))))
}
stems <- c(TK6 = "tk6", `MCF-7` = "mcf7", HepaRG = "heparg", `Liver spheroids` = "spheroids")
archive <- "../2023_Meier_FDR/bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10"
for (cell in names(stems)) {
  stem <- stems[[cell]]
  ids <- identity$run[identity$cell_line == cell]
  if (mode != "full") {
    ids <- intersect(ids, names(reference)); trials <- reference[ids]
  } else {
    paths <- file.path(out, sprintf("%s_%02d.rds", stem, 3:10))
    if (all(file.exists(paths))) {
      for (n in 3:10) validate_checkpoint(paths[n-2L],intersect(ids,identity$run[identity$n==n]))
      message("Completed cell, verified and skipping: ", cell); next
    }
    message(Sys.time(), " reading DEG archive: ", cell)
    trials <- readRDS(file.path(archive, paste0(stem, ".RDS")))
    stopifnot(length(trials) == length(ids))
    if (!is.null(names(trials))) stopifnot(identical(names(trials), ids))
    names(trials) <- ids
  }
  # No model refitting and no use of either post-test intermediate.
  for (n in 3:10) {
    path <- file.path(out, sprintf("%s_%02d.rds", stem, n))
    current <- intersect(ids, identity$run[identity$n == n])
    if (!length(current)) next
    if (file.exists(path)) { validate_checkpoint(path,current); next }
    start <- proc.time()
    rows <- parallel::mclapply(current, function(id) {
      tryCatch({
        input <- historical_ora_input(trials[[id]], a$alias_map)
        result <- reactome_hypergeom(input, a$pathways, compiled = compiled)
        list(run = id, input = input, result = result)
      }, error = function(e) list(run = id, result = list(status = "failed", error = conditionMessage(e))))
    }, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
    stopifnot(all(vapply(rows,is.list,logical(1))))
    names(rows) <- current
    payload <- list(schema = "reactome_ora_v1", rows = rows,
                    elapsed = (proc.time() - start)[["elapsed"]], session = sessionInfo(),
                    annotation_md5 = annotation_md5, source_md5 = source_md5, cores = cores)
    temp <- paste0(path, ".tmp")
    saveRDS(payload, temp, compress = "gzip")
    stopifnot(file.rename(temp, path))
    message(Sys.time(), " checkpoint ", path, " trials=", length(rows), " elapsed=", payload$elapsed)
  }
  rm(trials); gc()
}
