library(dplyr)

make_trial_metadata <- function(metadata,
                                cell_lines = cell_line_levels,
                                sample_sizes = sample_size_levels,
                                n_trials = n_trials_default,
                                seed = 1L) {
  set.seed(seed)

  metadata <- metadata %>%
    filter(chemical == "DMSO") %>%
    select(-any_of("group"))

  trial_metadata <- list()
  trial_names <- character()
  run_number <- 1L

  for (cell_line_name in cell_lines) {
    available_samples <- metadata %>%
      filter(cell_line == cell_line_name) %>%
      pull(Sample_ID)

    for (n_per_group in sample_sizes) {
      required_n <- n_per_group * 2L
      if (length(available_samples) < required_n) {
        stop(
          "Not enough DMSO samples for ",
          cell_line_name,
          " at n = ",
          n_per_group,
          "."
        )
      }

      for (trial_id in seq_len(n_trials)) {
        sampled <- sample(available_samples, size = required_n, replace = FALSE)
        trial <- data.frame(
          sample_id = sampled,
          group = factor(rep(c("A", "B"), each = n_per_group), levels = c("A", "B"))
        )

        trial_name <- paste(cell_line_name, n_per_group, n_per_group, trial_id, sep = "_")
        trial_metadata[[run_number]] <- trial
        trial_names[[run_number]] <- trial_name
        run_number <- run_number + 1L
      }
    }
  }

  names(trial_metadata) <- trial_names
  trial_metadata
}

