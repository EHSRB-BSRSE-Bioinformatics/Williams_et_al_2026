archive_relative_path <- function(path, root) {
    prefix <- paste0(root, .Platform$file.sep)
    in_archive <- startsWith(path, prefix)
    path[in_archive] <- substring(path[in_archive], nchar(prefix) + 1L)
    path
}

public_session_info <- function(lines = capture.output(sessionInfo())) {
    libraries <- grepl("^(BLAS|LAPACK)", lines)
    lines[libraries] <- gsub("/[^;[:space:]]*/", "", lines[libraries])
    lines
}
