# =============================================================================
# cohorts.R  --  the GEO cohorts the pipeline knows, and how a script picks one
#
# Sourced by download_geo.R and preprocess_geo.R. A cohort must be named
# explicitly on the command line, so a run can never silently use the wrong
# dataset:
#   --gse40279            or   --gse GSE40279
#   --gse87571            or   --gse GSE87571
#
# Roles: GSE40279 is the DEVELOPMENT cohort (tune feature subsets, aggregation,
# sex handling here). GSE87571 is the held-out VALIDATION cohort: run it once,
# with the settings frozen on GSE40279. Never tune on it.
# =============================================================================

COHORTS <- list(
  GSE40279 = list(role = "development",
                  desc = "Hannum 2013 (PMID 23177740): whole blood, 450K, n=656, ages 19-101"),
  GSE87571 = list(role = "validation",
                  desc = "Johansson 2013: whole blood, 450K, ~730 samples, ages ~14-94 (use --min-age 19 to match GSE40279)")
)

resolve_cohort <- function(args) {
  hits <- character(0)
  i <- which(args == "--gse")
  if (length(i)) hits <- c(hits, toupper(args[i[1] + 1]))
  flag <- grep("^--gse[0-9]+$", args, value = TRUE)
  if (length(flag)) hits <- c(hits, toupper(sub("^--", "", flag)))
  hits <- unique(hits)
  known <- paste(sprintf("  --%s   %s cohort: %s", tolower(names(COHORTS)),
                         sapply(COHORTS, `[[`, "role"), sapply(COHORTS, `[[`, "desc")), collapse = "\n")
  if (!length(hits)) stop("Name the cohort explicitly:\n", known)
  if (length(hits) > 1) stop("More than one cohort given (", paste(hits, collapse = ", "), "); run one at a time.")
  if (!hits %in% names(COHORTS)) stop("Unknown cohort ", hits, ". Known:\n", known)
  c(list(id = hits, slug = tolower(hits)), COHORTS[[hits]])
}
