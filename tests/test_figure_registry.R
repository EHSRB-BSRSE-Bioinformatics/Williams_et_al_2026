#!/usr/bin/env Rscript
script_file <- sub(
    "--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
    setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}
source("R/config.R")
source("R/figure_registry.R")

expect_error <- function(expr, pattern) {
    error <- tryCatch(
        {
            force(expr)
            NULL
        },
        error = identity
    )
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

registry <- figure_registry()
stopifnot(
    nrow(registry) == 30L,
    !anyDuplicated(registry$key),
    !anyDuplicated(file.path(registry$directory, registry$stem)),
    identical(
        registry$number[registry$category == "main"],
        paste0("Fig_", 2:6)
    ),
    identical(
        registry$number[registry$category == "supplementary"],
        paste0("Fig_S", 1:13)
    ),
    !"Fig_1" %in% registry$number
)
main_keys <- c("deg", "ora", "gsea", "tgx_ddi", "tpod_fdr_main")
for (i in seq_along(main_keys)) {
    stopifnot(
        figure_filename(main_keys[i], "png") ==
            paste0("figures/Fig_", i + 1L, ".png")
    )
}
stopifnot(
    figure_filename("s8_mean_variance", "pdf") ==
        "figures/supplementary/Fig_S8.pdf",
    figure_filename("s9", "png") == "figures/supplementary/Fig_S9.png",
    grepl("S9", figure_spec("s8_mean_variance")$note),
    grepl("false positives", figure_spec("s9")$note)
)
for (i in seq_len(nrow(registry))) {
    key <- registry$key[i]
    stopifnot(identical(
        figure_filename(key, c("pdf", "png")),
        figure_filename(registry$legacy_stem[i], c("pdf", "png"))
    ))
    if (nzchar(registry$number[i])) {
        stopifnot(identical(
            figure_filename(key),
            figure_filename(registry$number[i])
        ))
    } else {
        stopifnot(
            figure_filename(key, "png") ==
                paste0("figures/exploratory/", registry$legacy_stem[i], ".png")
        )
    }
}
expect_error(figure_filename("Fig_1", "png"), "Unknown")
expect_error(figure_filename("not_registered", "pdf"), "Unknown")
expect_error(figure_filename(NA_character_, "png"), "non-missing")
expect_error(figure_filename("deg", "../png"), "extension")

test_organizer <- function() {
    root <- tempfile("figure-registry-")
    dir.create(file.path(root, "figures"), recursive = TRUE)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    writeLines(
        "legacy DEG",
        file.path(root, "figures", "deg_monte_carlo_false_positive_rates.png")
    )
    writeLines(
        "legacy S8",
        file.path(root, "figures", "supplementary_figure_s8_mean_variance.pdf")
    )
    writeLines("unregistered", file.path(root, "figures", "diagnostic.png"))
    plan <- organize_figure_outputs(root, dry_run = TRUE)
    stopifnot(
        sum(plan$action == "move") == 3L,
        !dir.exists(file.path(root, "figures", "supplementary")),
        !file.exists(figure_path("deg", "png", root, create_dir = FALSE))
    )
    organize_figure_outputs(root)
    stopifnot(
        file.exists(figure_path("deg", "png", root)),
        file.exists(figure_path("s8_mean_variance", "pdf", root)),
        file.exists(file.path(root, "figures", "exploratory", "diagnostic.png"))
    )
    before <- tools::md5sum(list.files(
        file.path(root, "figures"),
        recursive = TRUE,
        full.names = TRUE
    ))
    organize_figure_outputs(root)
    after <- tools::md5sum(list.files(
        file.path(root, "figures"),
        recursive = TRUE,
        full.names = TRUE
    ))
    stopifnot(identical(before, after))
    manifest <- write_figure_manifest(root)
    stopifnot(
        nrow(manifest) == 61L,
        manifest$status[manifest$number == "Fig_1"] == "not_provided",
        all(manifest$status[manifest$number == "Fig_S9"] == "pending"),
        manifest$status[manifest$filename == "figures/Fig_2.png"] == "present"
    )
    manifest_hash <- tools::md5sum(file.path(
        root,
        "docs",
        "figure_manifest.tsv"
    ))
    write_figure_manifest(root)
    stopifnot(identical(
        manifest_hash,
        tools::md5sum(file.path(root, "docs", "figure_manifest.tsv"))
    ))

    # Identical duplicates are safe; conflicting destinations abort the whole plan.
    source <- file.path(
        root,
        "figures",
        "deg_monte_carlo_false_positive_rates.png"
    )
    file.copy(figure_path("deg", "png", root), source)
    plan <- organize_figure_outputs(root)
    stopifnot(
        "remove_identical_duplicate" %in% plan$action,
        !file.exists(source)
    )
    writeLines("different content", source)
    untouched <- file.path(root, "figures", "aaa_new.pdf")
    writeLines("must not move before conflict", untouched)
    before <- tools::md5sum(c(
        source,
        figure_path("deg", "png", root),
        untouched
    ))
    expect_error(organize_figure_outputs(root), "Conflicting")
    stopifnot(identical(before, tools::md5sum(names(before))))
}
test_organizer()

# Parse, but do not execute expensive plotting or load analysis packages.
producers <- unique(registry$script)
stopifnot(
    all(file.exists(producers)),
    identical(figure_spec("s1")$script, figure_spec("s3")$script),
    identical(figure_spec("s7")$script, figure_spec("s9")$script)
)
manifest <- read.delim("docs/figure_manifest.tsv", na.strings = NULL)
registered <- manifest[manifest$filename != "", ]
stopifnot(
    identical(
        registered$filename,
        unlist(Map(figure_filename, registry$key, registry$extensions), use.names = FALSE)
    ),
    identical(registered$script, rep(registry$script, lengths(registry$extensions))),
    all(registered$script_status == "available"),
    all(registered$status[registered$category != "exploratory"] == "present"),
    all(file.exists(registered$filename[registered$status == "present"]))
)
for (path in c(
    producers,
    "R/figure_registry.R",
    "scripts/organize_figures.R",
    list.files("R/revision", pattern = "\\.R$", full.names = TRUE)
)) {
    parse(path)
}
for (path in producers) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    stopifnot(
        grepl('source("R/figure_registry.R")', text, fixed = TRUE),
        grepl("figure_path|figure_filename", text)
    )
    # Old manuscript names must only appear in the registry, not producers.
    stopifnot(
        !any(vapply(
            registry$legacy_stem[1:18],
            function(stem) {
                grepl(paste0('"', stem, '"'), text, fixed = TRUE)
            },
            logical(1)
        ))
    )
}
ast <- parse("scripts/render_figures.R")
assignment <- Filter(
    function(expr) {
        is.call(expr) &&
            identical(expr[[1]], as.name("<-")) &&
            identical(expr[[2]], as.name("figure_specs"))
    },
    as.list(ast)
)[[1]]
keys <- names(as.list(assignment[[3]]))[-1]
stopifnot(setequal(keys, registry$key[registry$script == "scripts/render_figures.R"]))
for (path in producers) {
    walk <- function(expr) {
        if (missing(expr)) {
            return(invisible(NULL))
        }
        if (!is.call(expr)) {
            return(invisible(NULL))
        }
        if (
            as.character(expr[[1]])[1] %in%
                c("figure_path", "figure_filename", "figure_spec") &&
                length(expr) > 1L &&
                is.character(expr[[2]])
        ) {
            figure_spec(expr[[2]])
        }
        for (child in as.list(expr)[-1]) {
            walk(child)
        }
        invisible(NULL)
    }
    for (expr in parse(path)) {
        walk(expr)
    }
}
cat(
    "Figure registry, safe organizer, manifest, and producer parsing checks passed.\n"
)
