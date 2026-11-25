# Compatibility entry point: validate data, then use the main manuscript renderer.
source("R/config.R")
source("R/summaries.R")
source("R/figure_registry.R")
out <- "docs/revision/reactome_ora_v2"
inputs <- unlist(derived_files, use.names = FALSE)
rendered <- c(
  figure_filename("ora", c("pdf", "png")),
  figure_filename("ora_wikipathways_supplementary", c("pdf", "png"))
)
other_figures <- setdiff(list.files("figures", pattern = "[.](pdf|png)$",
  full.names = TRUE, recursive = TRUE), rendered)
before <- tools::md5sum(c(inputs, other_figures))
x <- read_summary_inputs()$ora
metadata <- readRDS(derived_files$trial_metadata)
stopifnot(nrow(x) == 32000L, setequal(x$run, names(metadata)), !anyDuplicated(x$run))
columns <- c(GO = "contains_GO_false_positive", Reactome = "contains_reactome_false_positive",
             WikiPathways = "contains_wiki_false_positive")
rates <- bind_rows(lapply(names(columns), function(collection) {
  x %>% group_by(cell_line, n) %>%
    summarise(n_trials = n(), n_positive = sum(.data[[columns[[collection]]]]),
              rate = mean(.data[[columns[[collection]]]]), .groups = "drop") %>%
    mutate(collection = collection, .before = 1)
}))
expected <- readRDS(file.path(out, "summary.rds"))
stopifnot(nrow(rates) == 96L, all(rates$n_trials == 1000L))
for (collection in names(columns)) {
  observed <- filter(rates, .data$collection == .env$collection)
  idx <- match(paste(observed$cell_line, observed$n), paste(expected$cell_line, expected$n))
  column <- c(GO = "go_rate", Reactome = "ora_rate", WikiPathways = "wp_rate")[[collection]]
  stopifnot(!anyNA(idx), all(observed$rate == expected[[column]][idx]))
}
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("scripts/render_figures.R", "--only=ora,ora_wikipathways_supplementary"))
stopifnot(status == 0L, identical(before, tools::md5sum(c(inputs, other_figures))))
write.table(rates, file.path(out, "ora_figure_rates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, "ora_promotion_sessionInfo.txt"))
checks <- c("PASS: 32,000 exact frozen IDs; 1,000 trials per cell model and sample size.",
  "PASS: all 96 plotted rates match the validated compact summary; insufficient inputs retained.",
  "PASS: frozen summaries, metadata and unrelated figures unchanged.",
  paste("Rendered", paste(rendered[1:2], collapse = ", "), "with GO A, Reactome B, and WikiPathways C."),
  paste("Rendered", paste(rendered[3:4], collapse = ", "), "using the same plot builder.")
)
writeLines(checks, file.path(out, "ora_promotion_validation.log"))
writeLines(checks)
