#' Phenotypes of selected traits in one QtlDataset context
#'
#' Wraps \code{pecotmr::molecularTraits()} with the behaviour of the pecotmr 0.6.x
#' \code{getPhenotypes(qd, contexts = context, traitId = trait_ids)} call: an unknown
#' context is an error, while traits missing from the context only raise a warning and
#' the traits present are returned in \code{trait_ids} order.
#'
#' @param qd A \code{pecotmr::QtlDataset}.
#' @param context Character, one context name.
#' @param trait_ids Character vector of trait (row) IDs; NULL returns every trait.
#' @return The context's \code{SummarizedExperiment} restricted to \code{trait_ids}.
#' @export
get_context_phenotypes <- function(qd, context, trait_ids = NULL) {
  available <- pecotmr::contexts(qd)
  if (!context %in% available) {
    stop("Unknown context(s): ", context, ". Available: ", paste(available, collapse = ", "))
  }
  se <- pecotmr::molecularTraits(qd, context = context)[[context]]
  if (is.null(trait_ids)) return(se)
  missing <- setdiff(trait_ids, rownames(se))
  if (length(missing) > 0L) {
    warning(sprintf("context '%s' is missing trait(s): %s", context, paste(missing, collapse = ", ")))
  }
  se[intersect(trait_ids, rownames(se)), , drop = FALSE]
}

#' Phenotype covariates of one QtlDataset context
#'
#' Samples x covariates matrix, as \code{getPhenotypeCovariates(qd, contexts = context)[[context]]}
#' returned in pecotmr 0.6.x.
#'
#' @param qd A \code{pecotmr::QtlDataset}.
#' @param context Character, one context name.
#' @return A samples x covariates matrix.
#' @export
get_context_covariates <- function(qd, context) {
  cov <- pecotmr::phenotypeCovariates(qd, context = context)
  if (is.list(cov) && !is.data.frame(cov)) cov <- cov[[context]]
  cov
}

#' Parse a genomic region string
#'
#' The pecotmr 0.6.11 \code{parseRegion()} (removed from later pecotmr), e.g.
#' \code{"chr1:100-200"} to a one-row data frame with \code{chrom}, \code{start}, \code{end}.
#'
#' @param region Region string or data frame accepted by pecotmr 0.6.11 \code{parseRegion()}.
#' @return A data frame with columns \code{chrom}, \code{start}, \code{end}.
#' @export
parse_region <- function(region) {
  get("parseRegion", envir = .pecotmr_legacy_ld)(region)
}
