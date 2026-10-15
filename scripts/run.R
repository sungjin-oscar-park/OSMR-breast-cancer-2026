#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
arg_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
REPO <- normalizePath(file.path(dirname(arg_file), ".."), mustWork = TRUE)
steps <- c("figure1-bulk", "figure1-preprocess", "figure1-copykat", "figure1-single-cell",
           "figure4", "figure4-analyze", "figure4-plot")
if (length(args) != 2L || !args[1] %in% steps) {
  cat("Usage: Rscript scripts/run.R <step> /absolute/path/to/config.R\nSteps:\n",
      paste(steps, collapse = "\n"), "\n")
  quit(status = 1L)
}
source(file.path(REPO, "R/config.R"))
config_file <- normalizePath(args[2], mustWork = TRUE)
cfg <- load_config(config_file, REPO)
source(file.path(REPO, "R/gene_sets.R"))
# Legacy libraries may create Rplots.pdf in the working directory.
# Run in the external output root to prevent such files appearing in the repo.
setwd(cfg$output_dir)
run_file <- function(path) sys.source(file.path(REPO, path), envir = .GlobalEnv)
step <- args[1]
if (step == "figure1-bulk") run_file("figure1/bulk.R")
if (step == "figure1-preprocess") run_file("figure1/sc_preprocess.R")
if (step == "figure1-copykat") run_file("figure1/copykat.R")
if (step == "figure1-single-cell") run_file("figure1/single_cell.R")
if (step %in% c("figure4", "figure4-analyze")) {
  run_file("figure4/preprocess_and_fit.R")
  run_file("figure4/analyze.R")
}
if (step %in% c("figure4", "figure4-plot")) run_file("figure4/plot.R")
run_dir <- output_path(cfg, file.path("provenance", step))
capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
invisible(file.copy(config_file, file.path(run_dir, "config.R"), overwrite = TRUE))
cat("Completed ", step, ": ", cfg$output_dir, "\n", sep = "")
