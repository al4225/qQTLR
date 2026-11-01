#' Extract cis genotypes with the pecotmr 0.6.x QtlDataset genotype QC
#'
#' Reads a PLINK1 genotype panel with \code{pecotmr::readGenotypes()} and applies
#' the variant/sample QC that pecotmr 0.6.x applied when a QtlDataset was built
#' (\code{mafCutoff}, \code{macCutoff}, \code{imissCutoff}, \code{xvarCutoff},
#' \code{keepIndel}, \code{keepSamples}). Later pecotmr releases dropped these
#' arguments from \code{loadQtlDatasetFromManifest()} and compute allele
#' frequencies only over the dataset's phenotyped samples; this keeps the 0.6.x
#' behaviour, where MAF and the MAC threshold use every genotyped sample.
#'
#' @param genotype_prefix Stem of a PLINK1 fileset (\code{.bed}/\code{.bim}/\code{.fam}),
#'   or a \code{RangedSummarizedExperiment} already returned by \code{readGenotypes()}.
#' @param region A \code{GRanges} of length one; variants with start inside it are kept.
#' @param maf_cutoff,mac_cutoff,imiss_cutoff,xvar_cutoff Numeric QC thresholds
#'   (0 disables a filter).
#' @param keep_indel Logical; if FALSE, variants whose A1 or A2 is not a single base are dropped.
#' @param keep_samples Optional character vector of samples to restrict to.
#' @return A list with \code{geno} (samples x variants dosage of A1, missing calls
#'   mean-imputed), \code{maf} and \code{af} (named by variant ID, computed before
#'   imputation), as returned by \code{getGenotypes()}/\code{getMaf()} in pecotmr 0.6.x.
#' @export
extract_qtl_genotypes <- function(genotype_prefix, region, maf_cutoff = 0, mac_cutoff = 0,
                                  imiss_cutoff = 0, xvar_cutoff = 0, keep_indel = TRUE,
                                  keep_samples = NULL) {
  panel <- if (methods::is(genotype_prefix, "SummarizedExperiment")) {
    genotype_prefix
  } else {
    pecotmr::readGenotypes(plink1Prefix = genotype_prefix)
  }
  rr <- SummarizedExperiment::rowRanges(panel)
  seq_match <- sub("^chr", "", as.character(GenomicRanges::seqnames(rr))) ==
    sub("^chr", "", as.character(GenomicRanges::seqnames(region))[1])
  snp_idx <- which(seq_match & GenomicRanges::start(rr) >= GenomicRanges::start(region)[1] &
                   GenomicRanges::start(rr) <= GenomicRanges::end(region)[1])
  empty <- list(geno = matrix(numeric(0), nrow = 0L, ncol = 0L,
                              dimnames = list(character(0), character(0))),
                maf = setNames(numeric(0), character(0)),
                af = setNames(numeric(0), character(0)))
  if (length(snp_idx) == 0L) return(empty)

  # Drop indels before materializing the dosage.
  if (!isTRUE(keep_indel)) {
    a1 <- as.character(S4Vectors::mcols(rr)$A1[snp_idx])
    a2 <- as.character(S4Vectors::mcols(rr)$A2[snp_idx])
    snp_idx <- snp_idx[which(nchar(a1) == 1L & nchar(a2) == 1L)]
    if (length(snp_idx) == 0L) return(empty)
  }

  # samples x variants dosage, missing calls left as NA
  dosage <- t(as.matrix(SummarizedExperiment::assay(panel[snp_idx, ], 1L)))
  keep <- rownames(dosage)
  if (length(keep_samples) > 0L) keep <- intersect(keep, keep_samples)
  dosage <- dosage[keep, , drop = FALSE]

  # Per-sample missingness filter
  if (imiss_cutoff > 0 && nrow(dosage) > 0L && ncol(dosage) > 0L) {
    dosage <- dosage[rowMeans(is.na(dosage)) <= imiss_cutoff, , drop = FALSE]
  }

  # Per-variant MAF / MAC / X-variance filters against the remaining sample count;
  # frequencies use the non-missing calls only.
  n_samp <- nrow(dosage)
  n_obs <- colSums(!is.na(dosage))
  sum_d <- colSums(dosage, na.rm = TRUE)
  p <- ifelse(n_obs > 0L, sum_d / (2 * n_obs), NA_real_)
  af <- p
  maf <- pmin(p, 1 - p)
  effective_maf <- max(maf_cutoff, if (n_samp > 0L) mac_cutoff / (2 * n_samp) else 0)
  keep_var <- !is.na(maf) & maf >= effective_maf
  if (xvar_cutoff > 0 && n_samp > 1L) {
    mu <- ifelse(n_obs > 0L, sum_d / n_obs, 0)
    centered <- sweep(dosage, 2L, mu, FUN = "-")
    centered[is.na(centered)] <- 0
    keep_var <- keep_var & colSums(centered * centered) / (n_samp - 1L) >= xvar_cutoff
  }
  dosage <- dosage[, keep_var, drop = FALSE]
  maf <- setNames(maf[keep_var], colnames(dosage))
  af <- setNames(af[keep_var], colnames(dosage))

  # Mean-impute remaining missing calls (frequencies were computed before imputation).
  if (anyNA(dosage)) {
    for (j in seq_len(ncol(dosage))) {
      col <- dosage[, j]
      na <- is.na(col)
      if (any(na)) {
        col[na] <- mean(col[!na])
        dosage[, j] <- col
      }
    }
  }
  list(geno = dosage, maf = maf, af = af)
}
