# Shared configuration; all data and outputs stay outside the code repository.
canonical_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/|^[A-Za-z]:[/\\\\]", path))
    stop("Use an absolute path: ", path, call. = FALSE)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  parent <- dirname(path)
  if (identical(path, parent)) stop("Cannot resolve path: ", path)
  file.path(canonical_path(parent), basename(path))
}
inside_repo <- function(path, repo) identical(path, repo) || startsWith(path, paste0(repo, "/"))
load_config <- function(path, repo) {
  cfg <- dget(path)
  if (!is.list(cfg)) stop("Configuration must evaluate to a list.")
  if (is.null(cfg$output_dir)) stop("Set output_dir in configuration.")
  for (key in c("output_dir", "tcga_dir", "scrna_raw_dir", "scrna_rds", "copykat_dir",
                "tpm_acute", "tpm_sustained", "hallmark_rds")) {
    if (is.null(cfg[[key]])) next
    cfg[[key]] <- canonical_path(cfg[[key]])
    if (inside_repo(cfg[[key]], repo)) stop(key, " must be outside the code repository.")
  }
  cfg$workers <- if (is.null(cfg$workers)) 1L else as.integer(cfg$workers)
  if (length(cfg$workers) != 1L || is.na(cfg$workers) || cfg$workers < 1L)
    stop("workers must be a positive integer.")
  if (is.null(cfg$font_family)) cfg$font_family <- "Arial"
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg
}
required_input <- function(cfg, key, fallback = NULL) {
  path <- cfg[[key]]
  if (is.null(path)) path <- fallback
  if (is.null(path) || !file.exists(path))
    stop("Missing input for ", key, ". See docs/INPUTS.md.", call. = FALSE)
  path
}
output_path <- function(cfg, subdir) {
  path <- file.path(cfg$output_dir, subdir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
