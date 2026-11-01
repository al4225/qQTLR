# Legacy TWAS harmonization code vendored from pecotmr.
#
# Source: StatFunGen/pecotmr commit 4843e8bd (2026-04-18, PR #472), the pecotmr
# version qQTLR was split from. quantile_twas_pipeline() and
# load_quantile_twas_weights() were written against these functions. Later
# pecotmr releases renamed them (#514), then removed or rewrote them for the S4 API
# (#501, #521; e.g. harmonizeTwas() is now a deprecation stub returning NULL).
#
# The functions are defined inside .pecotmr_legacy_twas_defs() and evaluated
# into a private environment by .onLoad (see zzz.R), so their names never mask
# qQTLR or pecotmr bindings. Access them with .legacy_twas("name").
# This is the transitive closure of harmonize_twas, twas_analysis, twas_z,
# find_data, clean_context_names and build_twas_score_row; the code is copied as-is
# apart from dropped roxygen comments and fixes marked "# qQTLR fix:".

.pecotmr_legacy_twas_defs <- function() {
  # ---- adjust_susie_weights (pecotmr R/susie_wrapper.R)
  adjust_susie_weights <- function(twas_weights_results, keep_variants, run_allele_qc = TRUE,
                                   variable_name_obj = c("susie_results", context, "variant_names"),
                                   susie_obj = c("susie_results", context, "susie_result_trimmed"),
                                   twas_weights_table = c("weights", context), LD_variants, match_min_prop = 0.2) {
    # Intersect the rownames of weights with keep_variants
    twas_weights_variants <- get_nested_element(twas_weights_results, variable_name_obj)
    # Normalize to canonical format (with chr prefix)
    twas_weights_variants <- normalize_variant_id(twas_weights_variants)
    # allele flip twas weights matrix variants name
    if (run_allele_qc) {
      weights_matrix <- get_nested_element(twas_weights_results, twas_weights_table)
      if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weights_matrix))) {
        weights_matrix <- cbind(parse_variant_id(twas_weights_variants), weights_matrix)
      }
      weights_matrix_qced <- allele_qc(weights_matrix, LD_variants, colnames(weights_matrix)[!colnames(weights_matrix) %in% c(
        "chrom",
        "pos", "A2", "A1"
      )], match_min_prop = match_min_prop)
      # allele_qc now outputs canonical variant_ids (with chr prefix)
      original_idx <- match(weights_matrix_qced$qc_summary$variants_id_original, twas_weights_variants)
      intersected_indices <- original_idx[weights_matrix_qced$qc_summary$keep == TRUE]
    } else {
      # Normalize keep_variants to canonical format for matching
      keep_variants_normalized <- normalize_variant_id(keep_variants)
      intersected_variants <- intersect(twas_weights_variants, keep_variants_normalized)
      intersected_indices <- match(intersected_variants, twas_weights_variants)
    }
    if (length(intersected_indices) == 0) {
      stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
    }
    # Subset lbf_matrix, mu, and x_column_scale_factors
    lbf_matrix <- get_nested_element(twas_weights_results, c(susie_obj, "lbf_variable"))
    mu <- get_nested_element(twas_weights_results, c(susie_obj, "mu"))
    x_column_scal_factors <- get_nested_element(twas_weights_results, c(susie_obj, "X_column_scale_factors"))
  
    lbf_matrix_subset <- lbf_matrix[, intersected_indices, drop = FALSE]
    mu_subset <- mu[, intersected_indices, drop = FALSE]
    x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]
  
    # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
    adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
    adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset
    # allele_qc now outputs canonical variant_ids (with chr prefix) -- no need to add chr
    remained_variants_ids <- if (run_allele_qc) {
      weights_matrix_qced$target_data_qced$variant_id
    } else {
      intersected_variants
    }
    return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = remained_variants_ids))
  }

  # ---- allele_qc (pecotmr R/allele_qc.R)
  allele_qc <- function(target_data, ref_variants, col_to_flip = NULL,
                       match_min_prop = 0.2, remove_dups = TRUE,
                       remove_indels = FALSE, remove_strand_ambiguous = TRUE,
                       flip_strand = FALSE, remove_unmatched = TRUE, ...) {
  	strand_flip <- function(ref) {
  	  # qQTLR fix: DNAStringSet() errors on non-IUPAC alleles such as the "*" spanning
  	  # deletion (present in ROSMAP genotypes); leave those unchanged so they simply fail
  	  # to match, and reverse-complement the rest exactly as before.
  	  out <- ref
  	  ok <- !is.na(ref) & !grepl("[^ACGTMRWSYKVHDBN.+-]", ref)
  	  out[ok] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ref[ok])))
  	  out
  	}
  
    # helper to sanitize column names to avoid NA/empty names that break dplyr verbs
    sanitize_names <- function(df) {
  	 nm <- colnames(df)
  	 if (is.null(nm)) {
  	   nm <- rep("unnamed", ncol(df))
  	 }
  	 empty_idx <- is.na(nm) | nm == ""
  	 if (any(empty_idx)) {
  	   # assign stable placeholder names for empties
  	   nm[empty_idx] <- paste0("unnamed_", seq_len(sum(empty_idx)))
  	 }
  	 # ensure names are unique and syntactic
  	 nm <- base::make.unique(nm, sep = "_")
  	 colnames(df) <- nm
  	 df
    }
  
    # check if the pattern is ATCG/DI
    check_ATCG <- function(vec) {
  	grepl("^[ATCGDI]+$", vec)
    }
  
    # transform all inputs to dataframe
    if (is.data.frame(target_data)) {
  	 if (ncol(target_data) > 4 && all(c("chrom", "pos", "A2", "A1") %in% names(target_data))) {
  		# Extract variant columns and standardize
  		variant_cols <- c("chrom", "pos", "A2", "A1")
  		variant_df <- target_data %>% select(all_of(variant_cols))
  		other_cols <- target_data %>% select(-all_of(variant_cols))
  		target_data <- cbind(variant_id_to_df(variant_df), other_cols)
  	 } else {
  		target_data <- variant_id_to_df(target_data)
  	 }
    } else {
  		target_data <- variant_id_to_df(target_data)
    }
    ref_variants <- variant_id_to_df(ref_variants)
  
    # Remove redundant columns that would conflict with the merge.
    # Keep A2/A1 in target_data — they become A2.target/A1.target after merge.
    columns_to_remove <- c("chromosome", "position", "ref", "alt", "variant_id")
    if (any(columns_to_remove %in% colnames(target_data))) {
  	 target_data <- select(target_data, -any_of(columns_to_remove))
    }
    match_result <- inner_join(target_data, ref_variants, by = c("chrom", "pos"), suffix = c(".target", ".ref")) %>% as.data.frame()
  
    # sanitize names after merge as well (merge can introduce empty names in edge cases)
    match_result <- sanitize_names(match_result)
  
    if (nrow(match_result) == 0) {
  	warning("No matching variants found between target data and reference variants.") 
  	return(list(target_data_qced = match_result, qc_summary = match_result))
    }
      # match target & ref by chrom and position
    match_result = match_result %>%
  	mutate(variants_id_original = format_variant_id(chrom, pos, A2.target, A1.target)) %>%
  	mutate(variants_id_qced = format_variant_id(chrom, pos, A2.ref, A1.ref)) %>%
  	# upper case target/reference A1 A2
  	mutate(across(c(A1.target, A2.target, A1.ref, A2.ref), toupper)) %>%
  	mutate(flip1.ref = strand_flip(A1.ref), flip2.ref = strand_flip(A2.ref)) %>%
  	# these pairings are ambiguous: because we cannot tell it's an sign flip / strand flip
  	mutate(strand_unambiguous = if_else((A1.target == "A" & A2.target == "T") | (A1.target == "T" & A2.target == "A") |
  	  (A1.target == "C" & A2.target == "G") | (A1.target == "G" & A2.target == "C"), FALSE, TRUE)) %>%
  	# filter out non-ATCG coded alleles
  	mutate(non_ATCG = !(check_ATCG(A1.target) & check_ATCG(A2.target))) %>%
  	# exact match should be kept all the time
  	mutate(exact_match = A1.target == A1.ref & A2.target == A2.ref) %>%
  	mutate(sign_flip = ((A1.target == A2.ref & A2.target == A1.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
  	mutate(strand_flip = ((A1.target == flip1.ref & A2.target == flip2.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
  	mutate(INDEL = (A2.target == "I" | A2.target == "D" | nchar(A2.target) > 1 | nchar(A1.target) > 1)) %>%
  	mutate(ID_match = ((A2.target == "D" | A2.target == "I") & (nchar(A1.ref) > 1 | nchar(A2.ref) > 1)))
  
    # if not remove, then this should'nt be a condition to filter out any variants
    if (!remove_strand_ambiguous) {
  	match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
    }
    # if all strand flip is un-ambigous, then we know ambigous cases are indeed a strand flip
    # not a sign flip, then we infer there is no ambigous in the whole dataset, and keep those ambigous ones
    if (!any(match_result$strand_flip & match_result$strand_unambiguous)) {
  	match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
    }
  
    # To keep variants: if it's a strand flip, we will keep those unambiguous (because if ambiguous, cannot know if it's strand / sign flip, so discard all)
    # or exact match or indel match (ID_match)
    # If not a strand flip, then we will keep those that are exact match / those are sign flip / INDEL matched
    match_result <- match_result %>% mutate(keep = if_else(strand_flip, true = (strand_unambiguous | exact_match | ID_match), false =
  	(exact_match | sign_flip | ID_match)
    ))
  
    if (remove_indels) {
  	match_result <- match_result %>% mutate(keep = if_else(INDEL, FALSE, keep))
    }
  
    # Flip the signs of col_to_flip for sign-flipped variants
    if (!is.null(col_to_flip)) {
  	missing <- setdiff(col_to_flip, colnames(match_result))
  	if (length(missing) > 0) {
  	  stop("Column(s) '", paste(missing, collapse = "', '"), "' not found in target_data.")
  	}
  	match_result[match_result$sign_flip, col_to_flip] <- -1 * match_result[match_result$sign_flip, col_to_flip]
    }
    # flip the strands if there is a strand flip
    if (flip_strand) {
  	strand_flipped_indices <- which(match_result$strand_flip)
  	match_result[strand_flipped_indices, "A1.target"] <- strand_flip(match_result[strand_flipped_indices, "A1.target"])
  	match_result[strand_flipped_indices, "A2.target"] <- strand_flip(match_result[strand_flipped_indices, "A2.target"])
    }
  
    # Remove all unnecessary columns used to determine qc status
    # Finally keep those variants with FLAG keep = TRUE
    result <- match_result[match_result$keep, , drop = FALSE]
  	
    if (remove_dups) {
  	dups <- duplicated(result[, c("chrom", "pos", "variants_id_qced")])
  	if (any(dups)) {
  	  n_removed <- sum(dups)
  	  warning(sprintf("Removed %d duplicate variant(s), keeping first occurrence.", n_removed))
  	  result <- result[!dups, , drop = FALSE]
  	}
    }
  
    qc_cols <- c("flip1.ref", "flip2.ref", "strand_unambiguous", "non_ATCG",
                 "exact_match", "sign_flip", "strand_flip", "INDEL", "ID_match", "keep")
    result <- result %>%
  	select(-any_of(qc_cols), -A1.target, -A2.target) %>%
  	rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
  
    if (!remove_unmatched) {
  	match_variant <- result %>% pull(variants_id_original)
  	match_result <- match_result %>%
  	  select(-any_of(qc_cols), -variants_id_original, -A1.target, -A2.target) %>%
  	  rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
  	target_data <- target_data %>% mutate(variant_id = format_variant_id(chrom, pos, A2, A1))
  	if (length(setdiff(target_data %>% pull(variant_id), match_variant)) > 0) {
  	  unmatch_data <- target_data %>% filter(!variant_id %in% match_variant)
  	  result <- rbind(result, unmatch_data %>% mutate(variants_id_original = variant_id))
  	  result <- result[match(target_data$variant_id, result$variants_id_original), ] %>% select(-variants_id_original)
  	}
    }
  
    min_match <- match_min_prop * nrow(ref_variants)
    if (nrow(result) < min_match) {
  	stop("Not enough variants have been matched.")
    }
  
    # Error if duplicated variant IDs remain (same variant with different values)
    if (any(duplicated(result$variant_id))) {
  	stop("Duplicated variants with different values found. Please check the input data and determine which to keep.")
    }
  
    return(list(target_data_qced = result, qc_summary = match_result))
  }

  # ---- build_twas_score_row (pecotmr R/twas.R)
  build_twas_score_row <- function(twas_rs, weight_db, context, study) {
    if (is.null(twas_rs)) return(data.frame())
    data.frame(
      gwas_study   = study,
      method       = sub("_[^_]+$", "", names(twas_rs)),
      twas_z       = find_data(twas_rs, c(2, "z")),
      twas_pval    = find_data(twas_rs, c(2, "pval")),
      context      = context,
      molecular_id = weight_db
    )
  }

  # ---- clean_context_names (pecotmr R/file_utils.R)
  clean_context_names <- function(context, gene) {
    # Remove gene name if it matches the last part of the context
    gene <- gene[order(-nchar(unique(gene)))]
    for (gene_id in gene) {
      context <- gsub(paste0("_", gene_id), "", context)
    }
    return(context)
  }

  # ---- compute_LD (pecotmr R/misc.R)
  compute_LD <- function(X, method = c("sample", "population", "gcta"),
                         trim_samples = FALSE, shrinkage = 0) {
    if (is.null(X)) {
      stop("X must be provided.")
    }
    method <- match.arg(method)
    nms <- colnames(X)
  
    if (method == "sample") {
      # ---- Standard sample correlation (N-1 denominator) ----
      # Mean impute only if NAs exist (PLINK2 data typically has none)
      X_imp <- X
      if (anyNA(X_imp)) {
        col_means <- colMeans(X_imp, na.rm = TRUE)
        na_pos <- which(is.na(X_imp), arr.ind = TRUE)
        X_imp[na_pos] <- col_means[na_pos[, 2]]
      }
      if (requireNamespace("Rfast", quietly = TRUE)) {
        # large=FALSE uses tcrossprod internally, ~40x faster than large=TRUE
        R <- Rfast::cora(X_imp, large = FALSE)
      } else {
        R <- cor(X_imp)
      }
    } else if (method == "population") {
      # ---- Population variance (N denominator, GCTA-style) ----
      # Optionally trim trailing samples to a multiple of 4 (matches .bed processing)
      if (trim_samples) {
        N_kept <- (nrow(X) %/% 4L) * 4L
        if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
      }
      N <- nrow(X)
      # Per-SNP means from non-missing values
      col_means <- colMeans(X, na.rm = TRUE)
      # Population variance: E[X^2] - E[X]^2
      col_vars <- colMeans(X^2, na.rm = TRUE) - col_means^2
      # Center; set NA -> 0 so missing pairs don't contribute to cross-products.
      # NOTE: the covariance divides by total N (not pairwise non-missing counts),
      # which is an approximation that assumes uniform missingness across SNPs.
      # With heterogeneous missingness, correlations between high-missing and
      # low-missing columns will be slightly deflated. This matches the GCTA
      # convention and is standard for PLINK-style LD computation.
      if (anyNA(X)) {
        na_rates <- colMeans(is.na(X))
        if (max(na_rates) - min(na_rates) > 0.1) {
          warning("Population LD method with heterogeneous missingness ",
                  "(max NA rate ", round(max(na_rates), 3),
                  ", min ", round(min(na_rates), 3),
                  "): correlations may be biased. Consider using method='sample' ",
                  "which handles missingness via mean imputation.")
        }
      }
      X_c <- sweep(X, 2, col_means)
      X_c[is.na(X_c)] <- 0
      # Covariance with N denominator
      cov_mat <- crossprod(X_c) / N
      # Correlation
      sd_vec <- sqrt(col_vars)
      R <- cov_mat / outer(sd_vec, sd_vec)
    } else {
      # ---- GCTA per-pair missing data correction ----
      # Matches the DENTIST binary's calcLDFromBfile_gcta formula exactly.
      # Unlike "population" which divides by total N, this method tracks
      # per-pair missing counts and applies a correction term.
      if (trim_samples) {
        N_kept <- (nrow(X) %/% 4L) * 4L
        if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
      }
      N <- nrow(X)
      p <- ncol(X)
  
      # Marginal statistics from non-missing values
      col_means <- colMeans(X, na.rm = TRUE)
      col_mean_sq <- colMeans(X^2, na.rm = TRUE)
      col_vars <- col_mean_sq - col_means^2
  
      # Build indicator matrix for non-missing values
      not_na <- !is.na(X)
      # Replace NA with 0 for cross-product computation
      X_zero <- X
      X_zero[is.na(X_zero)] <- 0
  
      # Per-pair non-missing counts: not_na'not_na gives count of jointly observed
      pair_counts <- crossprod(not_na * 1.0)
      n_missing <- N - pair_counts
  
      # Per-pair sums: sum of X_i over samples where both i and j are observed
      # For the correction term we need E_i2 = sum_i_pair / N (pair-specific mean)
      # X_zero' %*% not_na gives, for each (i,j), sum of X_i where j is not missing
      pair_sums <- crossprod(X_zero, not_na * 1.0)
  
      # Cross-product sum: sum(X_i * X_j) over jointly non-missing samples
      sum_XY <- crossprod(X_zero)
  
      # GCTA correction formula:
      # E_i2[i,j] = pair_sums[i,j] / N  (mean of SNP i restricted to non-missing-j samples, divided by N)
      # cov = sum_XY/N + E[i]*E[j]*(N-m)/N - E[i]*E_j2 - E_i2*E[j]
      E_i2 <- pair_sums / N  # p x p: row i, col j = sum of X_i where j non-missing, / N
      E_j2 <- t(E_i2)        # transposed version
  
      cov_mat <- sum_XY / N +
        outer(col_means, col_means) * (pair_counts / N) -
        col_means * E_j2 -
        E_i2 * rep(col_means, each = p)
  
      # Correlation
      sd_vec <- sqrt(col_vars)
      sd_outer <- outer(sd_vec, sd_vec)
      R <- matrix(0.001, p, p)
      valid <- sd_outer > 0
      R[valid] <- cov_mat[valid] / sd_outer[valid]
    }
  
    # Ensure clean output
    diag(R) <- 1.0
    R[is.na(R) | is.nan(R)] <- 0
  
    # Optional shrinkage toward identity: R_s = (1 - shrinkage) * R + shrinkage * I
    # Used e.g. by lassosum (Mak et al 2017) to regularize LD for RSS methods.
    if (shrinkage > 0 && shrinkage <= 1) {
      R <- (1 - shrinkage) * R + shrinkage * diag(nrow(R))
    }
  
    colnames(R) <- rownames(R) <- nms
    R
  }

  # ---- create_LD_matrix (pecotmr R/LD.R)
  create_LD_matrix <- function(LD_matrices, variants) {
    # Merge variant lists, deduplicating boundary overlaps
    merge_variants <- function(variant_list) {
      merged <- character(0)
      for (v in variant_list) {
        ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
        if (length(ids) == 0) next
        if (length(merged) > 0 && tail(merged, 1) == ids[1]) ids <- ids[-1]
        merged <- c(merged, ids)
      }
      merged
    }
  
    all_variants <- merge_variants(variants)
    combined <- matrix(0, nrow = length(all_variants), ncol = length(all_variants),
                       dimnames = list(all_variants, all_variants))
  
    # Place each block into the combined matrix
    for (i in seq_along(LD_matrices)) {
      v <- rownames(LD_matrices[[i]])
      idx <- match(v, all_variants)
      combined[idx, idx] <- LD_matrices[[i]]
    }
    combined
  }

  # ---- detect_variant_convention (pecotmr R/misc.R)
  detect_variant_convention <- function(ids) {
    # Find first non-NA element
    first_id <- ids[!is.na(ids)][1]
    if (is.na(first_id) || length(first_id) == 0) {
      return(list(has_chr = FALSE, allele_sep = ":", has_build = FALSE, example = NA_character_))
    }
    has_chr <- grepl("^chr", first_id)
    # Detect build suffix like :b38 or _b38 at end
    has_build <- grepl("(:|_)b[0-9]+$", first_id)
    id_clean <- strip_build_suffix(first_id)
    # Detect allele separator: check if variant uses underscores between allele fields
    # This catches both full underscore ("1_100_A_G") and mixed ("chr1:100_A_G") formats
    allele_sep <- if (grepl("_[ATCGID*]+_[ATCGID*]+$", id_clean)) "_" else ":"
    list(has_chr = has_chr, allele_sep = allele_sep, has_build = has_build, example = first_id)
  }

  # ---- extract_file_paths (pecotmr R/LD.R)
  extract_file_paths <- function(genomic_data, intersection_rows, column_to_extract) {
    if (!column_to_extract %in% names(genomic_data)) {
      stop("Column '", column_to_extract, "' not found in genomic data.")
    }
    idx <- which(genomic_data$chrom == intersection_rows$start_row$chrom &
                 genomic_data$start >= intersection_rows$start_row$start &
                 genomic_data$start <= intersection_rows$end_row$start)
    genomic_data[[column_to_extract]][idx]
  }

  # ---- extract_LD_for_region (pecotmr R/LD.R)
  extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
    extracted <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)
  
    if (!is.null(extract_coordinates)) {
      extract_coordinates <- extract_coordinates %>%
        mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
        select(chrom, pos)
      extracted <- extracted %>%
        mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
        merge(extract_coordinates, by = c("chrom", "pos"))
      keep_cols <- intersect(c("chrom", "variants", "pos", "GD", "A1", "A2",
                               "variance", "allele_freq", "n_nomiss"), names(extracted))
      extracted <- select(extracted, all_of(keep_cols))
    }
  
    mat <- LD_matrix[extracted$variants, extracted$variants, drop = FALSE]
    list(extracted_LD_matrix = mat, extracted_LD_variants = extracted)
  }

  # ---- find_data (pecotmr R/misc.R)
  find_data <- function(x, depth_obj, show_path = FALSE, rm_null = TRUE, rm_dup = FALSE, docall = c, last_obj = NULL) {
    depth <- as.integer(depth_obj[1])
    list_name <- if (length(depth_obj) > 1) depth_obj[2:length(depth_obj)] else NULL
    if (depth == 1 || depth == 0) {
      if (!is.null(list_name)) {
        if (list_name[1] %in% names(x)) {
          if (any(grepl("^[0-9]+$", list_name))) { # list names, indx name, list names
            second_depth <- which(grepl("^[0-9]+$", list_name))[1]
            data <- get_nested_element(x, list_name[1:second_depth[1] - 1])
            remaining_path <- list_name[second_depth:length(list_name)]
            return(find_data(data, remaining_path,
              show_path = show_path,
              rm_null = rm_null, rm_dup = rm_dup, last_obj = names(data)
            ))
          }
          return(get_nested_element(x, list_name))
        }
      } else {
        return(x)
      }
    } else if (is.list(x)) {
      result <- lapply(x, find_data,
        depth_obj = c(depth - 1, list_name), show_path = show_path,
        rm_null = rm_null, rm_dup = rm_dup, last_obj = names(x)
      )
      shared_list_names <- list()
      if (isTRUE(rm_null)) {
        result <- result[!sapply(result, is.null)]
        result <- result[!sapply(result, function(x) length(x) == 0)]
      }
      if (isTRUE(rm_dup)) {
        unique_result <- list()
        unique_counter <- 1
        for (i in seq_along(result)) {
          duplicate_found <- FALSE
          for (j in seq_along(unique_result)) {
            if (identical(result[[i]], unique_result[[j]])) {
              duplicate_found <- TRUE
              shared_list_names[[paste0("unique_list_", j)]] <- c(shared_list_names[[paste0("unique_list_", j)]], names(result)[i])
              break
            }
          }
          if (!duplicate_found) {
            unique_name <- paste0("unique_list_", unique_counter)
            unique_result[[names(result)[i]]] <- result[[i]]
            shared_list_names[[unique_name]] <- names(result)[i]
            unique_counter <- unique_counter + 1
          }
        }
        result <- unique_result
      }
  
      if (isTRUE(show_path)) {
        if (length(shared_list_names) > 0 & depth == 2) result$shared_list_names <- shared_list_names
        return(result) # Carry original list structure
      } else {
        flat_result <- do.call(docall, unname(result))
        if (length(shared_list_names) > 0 & depth == 2) {
          names(result) <- paste0("unique_list_", seq_along(result))
          result$shared_list_names <- shared_list_names
          return(result)
        } else {
          return(flat_result) # Only return values
        }
      }
    } else {
      message(paste0("list ", depth_obj[length(depth_obj)], " is not found in ", last_obj, ".  \n"))
    }
  }

  # ---- find_intersection_rows (pecotmr R/LD.R)
  find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
    chrom_data <- genomic_data %>% filter(chrom == region_chrom)
    if (nrow(chrom_data) == 0) stop("No data for chromosome ", region_chrom)
  
    # Clamp query to available range
    region_start <- max(region_start, min(chrom_data$start))
    region_end   <- min(region_end,   max(chrom_data$end))
  
    start_row <- genomic_data %>%
      filter(chrom == region_chrom, start <= region_start, end > region_start) %>%
      slice(1)
    end_row <- genomic_data %>%
      filter(chrom == region_chrom, start < region_end, end >= region_end) %>%
      arrange(desc(end)) %>%
      slice(1)
  
    if (nrow(start_row) == 0 || nrow(end_row) == 0) {
      stop("Region ", region_chrom, ":", region_start, "-", region_end,
           " is not covered by any rows in the LD metadata.")
    }
    list(start_row = start_row, end_row = end_row)
  }

  # ---- find_valid_file_path (pecotmr R/misc.R)
  find_valid_file_path <- function(reference_file_path, target_file_path) {
    # Check if the reference file path exits
    try_reference <- function() {
      if (file.exists(reference_file_path)) {
        return(reference_file_path)
      } else {
        return(NULL)
      }
    }
    # Check if the target file path exists
    try_target <- function() {
      if (file.exists(target_file_path)) {
        return(target_file_path)
      } else {
        # If not, construct a new target path by combining the directory of the reference file path with the target file path
        target_full_path <- file.path(dirname(reference_file_path), target_file_path)
        if (file.exists(target_full_path)) {
          return(target_full_path)
        } else {
          return(NULL)
        }
      }
    }
  
    target_result <- try_target()
    if (!is.null(target_result)) {
      return(target_result)
    }
  
    reference_result <- try_reference()
    if (!is.null(reference_result)) {
      return(reference_result)
    }
  
    stop(sprintf(
      "Both reference and target file paths do not work. Tried paths: '%s' and '%s'",
      reference_file_path, file.path(dirname(reference_file_path), target_file_path)
    ))
  }

  # ---- find_valid_file_paths (pecotmr R/misc.R)
  find_valid_file_paths <- function(reference_file_path, target_file_paths) sapply(target_file_paths, function(x) find_valid_file_path(reference_file_path, x))

  # ---- format_variant_id (pecotmr R/misc.R)
  format_variant_id <- function(chrom, pos, A2, A1, chr_prefix = TRUE, allele_sep = ":", convention = NULL) {
    # If convention is provided, use it to determine format automatically
    if (!is.null(convention)) {
      chr_prefix <- convention$has_chr
      allele_sep <- if (!is.null(convention$allele_sep)) convention$allele_sep else ":"
    }
    # Strip any existing chr prefix to normalize, then re-add if requested
    chrom_clean <- strip_chr_prefix(as.character(chrom))
    if (chr_prefix) {
      paste0("chr", chrom_clean, ":", pos, allele_sep, A2, allele_sep, A1)
    } else {
      paste0(chrom_clean, ":", pos, allele_sep, A2, allele_sep, A1)
    }
  }

  # ---- get_nested_element (pecotmr R/misc.R)
  get_nested_element <- function(nested_list, name_vector) {
    if (is.null(name_vector)) {
      return(NULL)
    }
    name_vector <- name_vector[name_vector!='']
    current_element <- nested_list
    for (name in name_vector) {
      if (is.null(current_element[[name]])) {
        stop("Element not found in the list")
      }
      current_element <- current_element[[name]]
    }
    return(current_element)
  }

  # ---- get_regional_ld_meta (pecotmr R/LD.R)
  get_regional_ld_meta <- function(ld_reference_meta_file, region, complete_coverage_required = FALSE) {
    genomic_data <- vroom(ld_reference_meta_file)
    region <- parse_region(region)
    # Set column names
    names(genomic_data) <- c("chrom", "start", "end", "path")
    names(region) <- c("chrom", "start", "end")
  
    # Treat start=0, end=0 as "covers all regions" (used for whole-chromosome PLINK files)
    whole_chrom <- genomic_data$start == 0 & genomic_data$end == 0
    if (any(whole_chrom)) genomic_data$end[whole_chrom] <- Inf
  
    # Order and deduplicate regions
    genomic_data <- order_dedup_regions(genomic_data)
    region <- order_dedup_regions(region)
  
    # Process file paths
    file_path <- genomic_data$path %>%
      str_split(",", simplify = TRUE) %>%
      data.frame() %>%
      `colnames<-`(if (ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))
  
    genomic_data <- cbind(genomic_data, file_path) %>% select(-path)
  
    # Find intersection rows
    intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)
  
    # Validate region
    if (complete_coverage_required) {
      validate_selected_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)
    }
  
    # Extract file paths
    LD_paths <- find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "LD_file_path"))
    bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
      find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "bim_file_path"))
    } else {
      NULL
    }
  
    return(list(
      intersections = list(
        start_index = intersection_rows$start_row,
        end_index = intersection_rows$end_row,
        LD_file_paths = LD_paths,
        bim_file_paths = bim_paths
      ),
      ld_meta_data = genomic_data,
      region = region
    ))
  }

  # ---- harmonize_gwas (pecotmr R/twas.R)
  harmonize_gwas <- function(gwas_file, query_region, ld_variants, col_to_flip=NULL, match_min_prop=0, column_file_path=NULL, comment_string="#"){
      if(is.null(gwas_file)| is.na(gwas_file)) stop("No GWAS file path provided. ")
      if (!is.null(column_file_path)) {  
        rss_result <- load_rss_data(  
          sumstat_path = gwas_file,  
          column_file_path = column_file_path,  
          region = query_region,
          comment_string = comment_string
        )  
        gwas_data_sumstats <- rss_result$sumstats  
      } else {  
        gwas_data_sumstats <- as.data.frame(tabix_region(gwas_file, query_region))  
      }
      if (nrow(gwas_data_sumstats) == 0) {
          if (length(names(gwas_file))==0) names(gwas_file) <- gwas_file
          warning(paste0("No GWAS summary statistics found for the region of ", query_region, " in ", names(gwas_file), ". "))
          return(NULL)
      }
      if (colnames(gwas_data_sumstats)[1] == "#chrom") colnames(gwas_data_sumstats)[1] <- "chrom" # colname update for tabix
      
      # Check if sumstats has z-scores or (beta and se)
      if (!is.null(gwas_data_sumstats$z)) {
        # z-scores already present, nothing to do
      } else if (!is.null(gwas_data_sumstats$beta) && !is.null(gwas_data_sumstats$se)) {  
        gwas_data_sumstats$z <- gwas_data_sumstats$beta / gwas_data_sumstats$se  
      } else {  
        stop("gwas_data_sumstats should have 'z' or ('beta' and 'se') columns")  
      }
      # check for overlapping variants
      if (!any(gwas_data_sumstats$pos %in% gsub("\\:.*$", "", sub("^.*?\\:", "", ld_variants)))) return(NULL)
      gwas_allele_flip <- allele_qc(gwas_data_sumstats, ld_variants, col_to_flip=col_to_flip, match_min_prop = match_min_prop)
      gwas_data_sumstats <- gwas_allele_flip$target_data_qced # post-qc gwas data that is flipped and corrected - gwas study level
      gwas_data_sumstats <- gwas_data_sumstats[!is.na(gwas_data_sumstats$z) & !is.infinite(gwas_data_sumstats$z), ]
      return(gwas_data_sumstats)
  }

  # ---- harmonize_twas (pecotmr R/twas.R)
  harmonize_twas <- function(twas_weights_data, ld_meta_file_path, gwas_meta_file,
                             ld_reference_sample_size, column_file_path = NULL, comment_string = "#") {
    # Function to group contexts based on start and end positions
    group_contexts_by_region <- function(twas_weights_data, molecular_id, chrom, tolerance = 5000) {
      region_info_df <- do.call(rbind, lapply(names(twas_weights_data$weights), function(context) {
        wgt_range <- parse_variant_id(rownames(twas_weights_data[["weights"]][[context]]))$pos
        data.frame(context = context, start = min(wgt_range), end = max(wgt_range))
      }))
      if (nrow(region_info_df) == 1) {
        # Handle case with only one context
        single_context_group <- list(
          context_group_1 = list(
            contexts = region_info_df$context,
            query_region = paste0(chrom, ":", region_info_df$start, "-", region_info_df$end),
            all_variants = unique(rownames(twas_weights_data[["weights"]][[region_info_df$context]]))
          )
        )
        return(single_context_group)
      }
      # Calculate distance matrix and perform hierarchical clustering
      clusters <- cutree(hclust(dist(region_info_df[, c("start", "end")])), h = tolerance)
      # Group contexts and determine query regions
      region_groups <- split(region_info_df, clusters) %>%
        lapply(function(group) {
          list(contexts = group$context, query_region = paste0(
            chrom, ":", min(group$start),
            "-", max(group$end)
          ))
        })
      # Create IRanges objects and merge overlapping intervals
      intervals <- IRanges(start = unlist(lapply(region_groups, function(context_group) {
        as.numeric(gsub(
          "^.*:\\s*|\\s*-.*$", "",
          context_group$query_region
        ))
      })), end = unlist(lapply(
        region_groups,
        function(context_group) as.numeric(sub("^.*?\\-", "", context_group$query_region))
      )))
      reduced_intervals <- reduce(intervals)
  
      # Find which original groups are merged, and update region_groups lists
      overlaps <- findOverlaps(intervals, reduced_intervals)
      # Create merged groups based on overlap mapping
      merged_groups <- lapply(seq_along(reduced_intervals), function(i) {
        context_indices <- queryHits(overlaps)[subjectHits(overlaps) == i]
        merged_contexts <- unlist(lapply(context_indices, function(idx) region_groups[[idx]]$contexts))
        list(contexts = merged_contexts, query_region = paste0(chrom, ":", start(reduced_intervals[i]), "-", end(reduced_intervals[i])))
      })
      names(merged_groups) <- paste0("context_group_", seq_along(merged_groups))
      # add variant names for coordinate extraction
      for (group in names(merged_groups)) {
        contexts <- merged_groups[[group]]$contexts
        merged_groups[[group]]$all_variants <- unique(do.call(c, lapply(
          contexts,
          function(context) {
            rownames(twas_weights_data[["weights"]][[context]])
          }
        )))
      }
      return(merged_groups)
    }
  
    # Step 1: load TWAS weights data
    molecular_ids <- names(twas_weights_data)
    chrom <- as.integer(parse_number(gsub(":.*$", "", rownames(twas_weights_data[[1]]$weights[[1]])[1])))
    gwas_meta_df <- as.data.frame(vroom(gwas_meta_file))
    gwas_files <- unique(gwas_meta_df$file_path[gwas_meta_df$chrom == chrom])
    names(gwas_files) <- unique(gwas_meta_df$study_id[gwas_meta_df$chrom == chrom])
    results <- list()
  
    # Step 2: Load LD for all events/genes by clustered context region
    for (molecular_id in molecular_ids) {
      twas_weights_data[[molecular_id]][["variant_names"]] <- lapply(twas_weights_data[[molecular_id]]$weights, function(x) rownames(x))
    }
    region_variants <- variant_id_to_df(unique(do.call(c, find_data(twas_weights_data, c(2, "variant_names")))))
    region_of_interest <- data.frame(chrom = chrom, start = min(region_variants$pos), end = max(region_variants$pos))
    LD_list <- load_LD_matrix(ld_meta_file_path, region_of_interest, region_variants,
                              n_sample = ld_reference_sample_size)
    # remove duplicate variants
    dup_idx <- which(duplicated(LD_list$LD_variants))
    if (length(dup_idx) >= 1) {
      LD_list$LD_variants <- LD_list$LD_variants[-dup_idx]
      LD_list$LD_matrix <- LD_list$LD_matrix[-dup_idx, -dup_idx]
      LD_list$ref_panel <- LD_list$ref_panel[-dup_idx, ]
    }
  
    # loop through genes/events:
    for (molecular_id in molecular_ids) {
      mol_data <- twas_weights_data[[molecular_id]]
      mol_res <- list(chrom = chrom, variant_names = list())
      mol_res[["data_type"]] <- if ("data_type" %in% names(mol_data)) mol_data$data_type
  
      # group contexts based on the variant position
      context_clusters <- group_contexts_by_region(mol_data, molecular_id, chrom, tolerance = 5000)
  
      # loop through contexts: grouping contexts can be useful during TWAS data harmonization to stratify variants for LD loading
      for (context_group in names(context_clusters)) {
        cluster <- context_clusters[[context_group]]
        contexts <- cluster$contexts
        query_region <- cluster$query_region
        region_of_interest <- region_to_df(query_region)
        all_variants <- variant_id_to_df(cluster$all_variants)
  
        # Step 3: load GWAS data for clustered context groups
        for (study in names(gwas_files)) {
          gwas_file <- gwas_files[study]
          gwas_data_sumstats <- harmonize_gwas(gwas_file, query_region=query_region,
                                                LD_list$LD_variants, c("beta", "z"),
                                                match_min_prop = 0, column_file_path = column_file_path, comment_string = comment_string)
          if(is.null(gwas_data_sumstats)) next
          # loop through context within the context group:
          for (context in contexts) {
            weights_matrix <- mol_data[["weights"]][[context]]
  
            # Step 4: harmonize weights, flip allele
            weights_matrix <- cbind(variant_id_to_df(rownames(weights_matrix)), weights_matrix)
            weights_matrix_qced <- allele_qc(weights_matrix, LD_list$LD_variants,
              colnames(weights_matrix)[!colnames(weights_matrix) %in% c("chrom", "pos", "A2", "A1")],
              match_min_prop = 0
            )
            qced_data <- weights_matrix_qced$target_data_qced
            weights_matrix_subset <- as.matrix(qced_data[, !colnames(qced_data) %in% c(
              "chrom", "pos", "A2", "A1", "variant_id", "variants_id_original"
            ), drop = FALSE])
            rownames(weights_matrix_subset) <- qced_data$variant_id
  
            # intersect post-qc gwas and post-qc weight variants (all now in canonical chr-prefix format)
            gwas_LD_variants <- intersect(gwas_data_sumstats$variant_id, LD_list$LD_variants)
            weights_matrix_subset <- weights_matrix_subset[rownames(weights_matrix_subset) %in% gwas_data_sumstats$variant_id, , drop = FALSE]
            if (nrow(weights_matrix_subset) == 0) next
            postqc_weight_variants <- rownames(weights_matrix_subset)
  
            # Step 5: adjust SuSiE weights based on available variants
            if ("susie_weights" %in% colnames(mol_data[["weights"]][[context]])) {
              adjusted_susie_weights <- adjust_susie_weights(mol_data,
                keep_variants = postqc_weight_variants, run_allele_qc = TRUE,
                variable_name_obj = c("variant_names", context),
                susie_obj = c("susie_results", context),
                twas_weights_table = c("weights", context), postqc_weight_variants, match_min_prop = 0
              )
              weights_matrix_subset <- cbind(
                susie_weights = setNames(adjusted_susie_weights$adjusted_susie_weights, adjusted_susie_weights$remained_variants_ids),
                weights_matrix_subset[adjusted_susie_weights$remained_variants_ids, !colnames(weights_matrix_subset) %in% "susie_weights", drop = FALSE]
              )
              susie_intermediate <- mol_data$susie_results[[context]][c("pip", "cs_variants", "cs_purity")]
              names(susie_intermediate[["pip"]]) <- rownames(weights_matrix) # original variants that is not qced yet
              pip <- susie_intermediate[["pip"]]
              pip_qced <- allele_qc(cbind(parse_variant_id(names(pip)), pip), LD_list$LD_variants, "pip", match_min_prop = 0)
              susie_intermediate[["pip"]] <- abs(pip_qced$target_data_qced$pip)
              names(susie_intermediate[["pip"]]) <- pip_qced$target_data_qced$variant_id
              susie_intermediate[["cs_variants"]] <- lapply(susie_intermediate[["cs_variants"]], function(x) {
                variant_qc <- allele_qc(x, LD_list$LD_variants, match_min_prop = 0)
                variant_qc$target_data_qced$variant_id[variant_qc$target_data_qced$variant_id %in% postqc_weight_variants]
              })
              mol_res[["susie_weights_intermediate_qced"]][[context]] <- susie_intermediate
            }
            rm(weights_matrix) # context specific original weight matrix
            gc()
  
            if (nrow(weights_matrix_subset) == 0) {
                warning("weights_matrix_subset is empty. Skipping this context.")
                next
            }
            mol_res[["variant_names"]][[context]][[study]] <- rownames(weights_matrix_subset)
  
            # Step 6: scale weights by variance (from ref_panel, populated by load_LD_matrix)
            variance <- LD_list$ref_panel$variance[match(rownames(weights_matrix_subset), LD_list$ref_panel$variant_id)]
            mol_res[["weights_qced"]][[context]][[study]] <- list(scaled_weights = weights_matrix_subset * sqrt(variance), weights = weights_matrix_subset)
          }
          # Combine gwas sumstat across different context for a single context group (all variant_ids now in canonical format)
          gwas_data_sumstats <- gwas_data_sumstats[gwas_data_sumstats$variant_id %in% unique(find_data(mol_res[["variant_names"]], c(2, study))), , drop = FALSE]
          mol_res[["gwas_qced"]][[study]] <- rbind(mol_res[["gwas_qced"]][[study]], gwas_data_sumstats)
          gwas_qced <- mol_res[["gwas_qced"]][[study]]
          mol_res[["gwas_qced"]][[study]] <- gwas_qced[!duplicated(gwas_qced[, c("variant_id", "z")]), ]
        }
      }
      twas_weights_data[[molecular_id]] <- NULL
      # extract LD matrix for variants intersect with gwas and twas weights at molecular_id level
      all_molecular_variants <- unique(find_data(mol_res[["gwas_qced"]], c(2, "variant_id")))
      if (is.null(all_molecular_variants)) {
        results[[molecular_id]] <- NULL
      } else {
        # All variant IDs are now in canonical chr-prefix format
        var_indx <- match(all_molecular_variants, LD_list$LD_variants)
        # qQTLR fix: drop = FALSE keeps dimnames when a gene has a single variant;
        # otherwise twas_analysis() finds no variants in the 1x1 LD matrix and skips the gene.
        mol_res[["LD"]] <- as.matrix(LD_list$LD_matrix[var_indx, var_indx, drop = FALSE])
        results[[molecular_id]] <- mol_res
      }
    }
    # return results
    return(list(twas_data_qced = results, ref_panel = LD_list$ref_panel))
  }

  # ---- has_plink1_files (pecotmr R/LD.R)
  has_plink1_files <- function(prefix) {
    file.exists(paste0(prefix, ".bed")) &&
      file.exists(paste0(prefix, ".bim")) &&
      file.exists(paste0(prefix, ".fam"))
  }

  # ---- has_plink2_files (pecotmr R/LD.R)
  has_plink2_files <- function(prefix) {
    file.exists(paste0(prefix, ".pgen")) &&
      (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst"))) &&
      file.exists(paste0(prefix, ".psam"))
  }

  # ---- invert_minmax_scaling (pecotmr R/file_utils.R)
  invert_minmax_scaling <- function(X, u_min, u_max) {
    if (length(u_min) != ncol(X) || length(u_max) != ncol(X)) {
      stop("Length of u_min/u_max (", length(u_min), ") must equal ncol(X) (", ncol(X), ")")
    }
    denom <- u_max - u_min
    denom[denom == 0] <- 1  # monomorphic: scaling was identity
    # Invert: U_original = U_scaled * (u_max - u_min) / 2 + u_min
    sweep(sweep(X, 2, denom / 2, "*"), 2, u_min, "+")
  }

  # ---- is_snp_alleles (pecotmr R/misc.R)
  is_snp_alleles <- function(a1, a2) {
    nchar(a1) == 1L & nchar(a2) == 1L &
      grepl("^[ATCG]$", a1) & grepl("^[ATCG]$", a2)
  }

  # ---- lbf_to_alpha (pecotmr R/susie_wrapper.R)
  lbf_to_alpha <- function(lbf) {
    alpha_matrix <- t(apply(as.matrix(lbf), 1, lbf_to_alpha_vector))
    if (ncol(lbf) == 1) alpha_matrix <- matrix(alpha_matrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
    return(alpha_matrix)
  }

  # ---- load_LD_from_blocks (pecotmr R/LD.R)
  load_LD_from_blocks <- function(LD_meta_file_path, region, extract_coordinates = NULL, n_sample = NULL) {
    # Intersect LD metadata with specified regions
    intersected_LD_files <- get_regional_ld_meta(LD_meta_file_path, region)
  
    LD_file_paths <- intersected_LD_files$intersections$LD_file_paths
    bim_file_paths <- intersected_LD_files$intersections$bim_file_paths
  
    extracted_LD_matrices_list <- list()
    extracted_LD_variants_list <- list()
    block_chroms <- character(length(LD_file_paths))
  
    for (j in seq_along(LD_file_paths)) {
      LD_matrix_processed <- process_LD_matrix(LD_file_paths[j], bim_file_paths[j])
      extracted_LD_list <- extract_LD_for_region(
        LD_matrix = LD_matrix_processed$LD_matrix,
        variants = LD_matrix_processed$LD_variants,
        region = intersected_LD_files$region,
        extract_coordinates = extract_coordinates
      )
      extracted_LD_matrices_list[[j]] <- extracted_LD_list$extracted_LD_matrix
      extracted_LD_variants_list[[j]] <- extracted_LD_list$extracted_LD_variants
      if (nrow(extracted_LD_variants_list[[j]]) > 0) {
        block_chroms[j] <- as.character(extracted_LD_variants_list[[j]]$chrom[1])
      } else {
        block_chroms[j] <- as.character(intersected_LD_files$region$chrom)
      }
    }
  
    # Filter out empty blocks before combining
    non_empty <- sapply(extracted_LD_variants_list, function(v) nrow(v) > 0)
    if (!any(non_empty)) {
      stop("No variants found in any LD block for the specified region.")
    }
    if (any(!non_empty)) {
      message(paste(
        "Removing", sum(!non_empty), "empty LD block(s) with no variants in the region."
      ))
      extracted_LD_matrices_list <- extracted_LD_matrices_list[non_empty]
      extracted_LD_variants_list <- extracted_LD_variants_list[non_empty]
      block_chroms <- block_chroms[non_empty]
      LD_file_paths <- LD_file_paths[non_empty]
    }
  
    LD_matrix <- create_LD_matrix(
      LD_matrices = extracted_LD_matrices_list,
      variants = extracted_LD_variants_list
    )
    LD_variants <- rownames(LD_matrix)
  
    block_variants <- lapply(extracted_LD_variants_list, function(v) v$variants)
    block_positions <- lapply(extracted_LD_variants_list, function(v) v$pos)
    block_metadata <- data.frame(
      block_id = seq_along(LD_file_paths),
      chrom = block_chroms,
      block_start = sapply(block_positions, min),
      block_end = sapply(block_positions, max),
      size = sapply(block_variants, length),
      start_idx = sapply(block_variants, function(v) min(match(v, LD_variants))),
      end_idx = sapply(block_variants, function(v) max(match(v, LD_variants))),
      stringsAsFactors = FALSE
    )
  
    rm(extracted_LD_matrices_list)
  
    ref_panel <- parse_variant_id(rownames(LD_matrix))
    merged_variant_list <- do.call(rbind, extracted_LD_variants_list)
    ref_panel$variant_id <- rownames(LD_matrix)
  
    if ("allele_freq" %in% colnames(merged_variant_list)) {
      ref_panel$allele_freq <- merged_variant_list$allele_freq[match(rownames(LD_matrix), merged_variant_list$variants)]
    }
    if ("variance" %in% colnames(merged_variant_list)) {
      ref_panel$variance <- merged_variant_list$variance[match(rownames(LD_matrix), merged_variant_list$variants)]
    }
    if ("n_nomiss" %in% colnames(merged_variant_list)) {
      ref_panel$n_nomiss <- merged_variant_list$n_nomiss[match(rownames(LD_matrix), merged_variant_list$variants)]
    }
  
    # Compute variance from n_sample + allele_freq if not already present
    if (!is.null(n_sample) && (!"variance" %in% colnames(ref_panel) || all(is.na(ref_panel$variance)))) {
      if ("allele_freq" %in% colnames(ref_panel)) {
        p <- ref_panel$allele_freq
        ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
        ref_panel$n_nomiss <- n_sample
      }
    }
  
    list(
      LD_variants = LD_variants,
      LD_matrix = LD_matrix,
      ref_panel = ref_panel,
      block_metadata = block_metadata,
      is_genotype = FALSE
    )
  }

  # ---- load_LD_from_genotype (pecotmr R/LD.R)
  load_LD_from_genotype <- function(prefix, region, source_type,
                                    return_genotype = FALSE, n_sample = NULL) {
    # Load genotype matrix and variant info
    result <- if (source_type == "plink2") {
      load_plink2_data(prefix, region = region)
    } else {
      load_plink1_data(prefix, region = region)
    }
    X <- result$X
    variant_info <- result$variant_info
  
    # Normalize variant IDs to canonical format (chr:pos:A2:A1)
    variant_ids <- normalize_variant_id(
      format_variant_id(variant_info$chrom, variant_info$pos, variant_info$A2, variant_info$A1)
    )
    colnames(X) <- variant_ids
  
    # Build ref_panel
    ref_panel <- parse_variant_id(variant_ids)
    ref_panel$variant_id <- variant_ids
  
    # Load allele frequency from .afreq file (required for PLINK sources)
    afreq <- read_afreq(prefix)
    if (is.null(afreq)) {
      stop("Allele frequency file (.afreq or .afreq.zst) not found at prefix: ", prefix,
           "\n  The .afreq file is required for PLINK genotype LD sources.")
    }
    freq_match <- match(variant_info$id, afreq$id)
    n_unmatched <- sum(is.na(freq_match))
    if (n_unmatched > 0) {
      warning(n_unmatched, " out of ", length(freq_match),
              " variants have no allele frequency in .afreq file.")
    }
    ref_panel$allele_freq <- afreq$alt_freq[freq_match]
  
    # Compute variance if sample size provided
    if (!is.null(n_sample)) {
      p <- ref_panel$allele_freq
      ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
      ref_panel$n_nomiss <- n_sample
    }
  
    # Block metadata (single block spanning the loaded region)
    positions <- variant_info$pos
    block_metadata <- data.frame(
      block_id = 1L,
      chrom = as.character(variant_info$chrom[1]),
      block_start = min(positions),
      block_end = max(positions),
      size = length(variant_ids),
      start_idx = 1L,
      end_idx = length(variant_ids),
      stringsAsFactors = FALSE
    )
  
    if (return_genotype) {
      return(list(
        LD_variants = variant_ids,
        LD_matrix = X,
        ref_panel = ref_panel,
        block_metadata = block_metadata,
        is_genotype = TRUE
      ))
    }
  
    R <- compute_LD(X, method = "sample")
  
    list(
      LD_variants = variant_ids,
      LD_matrix = R,
      ref_panel = ref_panel,
      block_metadata = block_metadata,
      is_genotype = FALSE
    )
  }

  # ---- load_LD_matrix (pecotmr R/LD.R)
  load_LD_matrix <- function(LD_meta_file_path, region, extract_coordinates = NULL,
                             return_genotype = FALSE, n_sample = NULL) {
    source <- resolve_ld_source(LD_meta_file_path)
    is_plink <- source$type %in% c("plink2", "plink1")
  
    # "auto": return X for PLINK, R for pre-computed
    if (identical(return_genotype, "auto")) return_genotype <- is_plink
  
    if (is_plink) {
      prefix <- resolve_plink_prefix_for_region(source$meta_path, region, source$type)
      return(load_LD_from_genotype(prefix, region, source$type,
                                   return_genotype = return_genotype,
                                   n_sample = n_sample))
    }
  
    # Pre-computed LD blocks (.cor.xz)
    if (return_genotype) {
      stop("return_genotype=TRUE requires PLINK genotype files, not pre-computed LD matrices.")
    }
    load_LD_from_blocks(source$meta_path, region, extract_coordinates, n_sample = n_sample)
  }

  # ---- load_plink1_data (pecotmr R/file_utils.R)
  load_plink1_data <- function(prefix, region = NULL, keep_indel = TRUE, keep_variants_path = NULL) {
    bed_file <- paste0(prefix, ".bed")
    bim_file <- paste0(prefix, ".bim")
    fam_file <- paste0(prefix, ".fam")
    if (!all(file.exists(bed_file, bim_file, fam_file))) {
      stop("PLINK1 fileset (.bed/.bim/.fam) not found at prefix: ", prefix)
    }
    if (!requireNamespace("snpStats", quietly = TRUE)) {
      stop("snpStats is required. Install from https://bioconductor.org/packages/release/bioc/html/snpStats.html")
    }
  
    # --- Region filter via bim ---
    if (!is.null(region)) {
      parsed <- parse_region(region)
      col_types <- list(col_character(), col_character(), col_guess(), col_integer(), col_guess(), col_guess())
      bim_data <- vroom(bim_file, col_names = FALSE, col_types = col_types)
      bim_data$X1 <- strip_chr_prefix(bim_data$X1)
      snp_ids <- bim_data$X2[bim_data$X1 == parsed$chrom &
                              bim_data$X4 >= parsed$start &
                              bim_data$X4 <= parsed$end]
      if (length(snp_ids) == 0) {
        stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
      }
    } else {
      snp_ids <- NULL
    }
  
    geno <- snpStats::read.plink(prefix, select.snps = snp_ids)
    geno_map <- geno$map
  
    # --- Build variant_info from snpStats $map ---
    # snpStats: allele.1 = effect allele (A1), allele.2 = reference allele (A2)
    variant_info <- data.frame(
      chrom = as.character(geno_map$chromosome),
      id    = rownames(geno_map),
      pos   = geno_map$position,
      A2    = as.character(geno_map$allele.2),
      A1    = as.character(geno_map$allele.1),
      stringsAsFactors = FALSE
    )
  
    # --- Dosage matrix: 2 - snpStats encoding gives A1 dosage ---
    X <- 2 - as(geno$genotypes, "numeric")
    rownames(X) <- rownames(geno$genotypes)
    colnames(X) <- variant_info$id
  
    # --- Post-filters: indels and variant whitelist ---
    if (!keep_indel) {
      snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
      X <- X[, snp_mask, drop = FALSE]
      variant_info <- variant_info[snp_mask, , drop = FALSE]
    }
    if (!is.null(keep_variants_path)) {
      keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
      X <- X[, keep_idx, drop = FALSE]
      variant_info <- variant_info[keep_idx, , drop = FALSE]
    }
  
    return(list(X = X, variant_info = variant_info))
  }

  # ---- load_plink2_data (pecotmr R/file_utils.R)
  load_plink2_data <- function(prefix, region = NULL, keep_indel = TRUE, keep_variants_path = NULL) {
    if (!requireNamespace("pgenlibr", quietly = TRUE)) {
      stop("pgenlibr is required. Install from https://cran.r-project.org/web/packages/pgenlibr/index.html")
    }
  
    paths <- resolve_plink2_paths(prefix)
  
    # --- Read variant info from .pvar as text (pgenlibr::NewPvar is unreliable) ---
    all_variant_info <- read_pvar_text(paths$pvar)
  
    variant_idx <- seq_len(nrow(all_variant_info))
    if (!is.null(region)) {
      parsed <- parse_region(region)
      in_region <- strip_chr_prefix(all_variant_info$chrom) == parsed$chrom &
                   all_variant_info$pos >= parsed$start &
                   all_variant_info$pos <= parsed$end
      variant_idx <- which(in_region)
      if (length(variant_idx) == 0) {
        stop(NoSNPsError(paste("No variants found in region", region)))
      }
    }
    variant_info <- all_variant_info[variant_idx, , drop = FALSE]
  
    # --- Read samples from .psam ---
    psam <- as.data.frame(vroom(paths$psam, delim = "\t", show_col_types = FALSE))
    colnames(psam)[1:2] <- c("FID", "IID")
  
    # --- Read genotype dosage via pgenlibr ---
    pgen <- pgenlibr::NewPgen(paths$pgen)
    on.exit(pgenlibr::ClosePgen(pgen), add = TRUE)
    X <- pgenlibr::ReadList(pgen, variant_subset = variant_idx, meanimpute = FALSE)
    rownames(X) <- psam$IID
    colnames(X) <- variant_info$id
  
    # --- Attach allele frequency and rescale stochastic genotypes ---
    afreq <- read_afreq(prefix)
    if (!is.null(afreq)) {
      afreq_cols <- intersect(c("id", "alt_freq", "obs_ct", "u_min", "u_max"), colnames(afreq))
      variant_info <- merge(variant_info, afreq[, afreq_cols, drop = FALSE],
                            by = "id", all.x = TRUE, sort = FALSE)
    }
    # Detect stochastic genotype (non-integer dosage, e.g., from rss_ld_sketch).
    # If U_MIN/U_MAX are in .afreq, invert the [0,2] min-max scaling exactly.
    is_stochastic <- !all(X == round(X), na.rm = TRUE)
    if (is_stochastic) {
      if ("u_min" %in% colnames(variant_info) && "u_max" %in% colnames(variant_info)) {
        idx <- match(colnames(X), variant_info$id)
        X <- invert_minmax_scaling(X, variant_info$u_min[idx], variant_info$u_max[idx])
        message("Stochastic genotype detected: restored original scale via U_MIN/U_MAX from .afreq")
      } else {
        warning("Non-integer genotype values detected (possible stochastic genotype from rss_ld_sketch). ",
                "Provide a .afreq file with U_MIN/U_MAX columns to restore original scale. ",
                "Without inversion, correlation structure is preserved but U'U/B may not approximate R.")
      }
    }
  
    # --- Post-filters: indels and variant whitelist ---
    if (!keep_indel) {
      snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
      X <- X[, snp_mask, drop = FALSE]
      variant_info <- variant_info[snp_mask, , drop = FALSE]
    }
    if (!is.null(keep_variants_path)) {
      keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
      X <- X[, keep_idx, drop = FALSE]
      variant_info <- variant_info[keep_idx, , drop = FALSE]
    }
  
    list(X = X, variant_info = variant_info)
  }

  # ---- load_rss_data (pecotmr R/file_utils.R)
  load_rss_data <- function(sumstat_path, column_file_path, n_sample = 0, n_case = 0, n_control = 0, region = NULL,
                            extract_region_name = NULL, region_name_col = NULL, comment_string = "#") {
    # Validate input files exist
    if (!file.exists(sumstat_path)) {
      stop("Summary statistics file not found: ", sumstat_path)
    }
    if (!file.exists(column_file_path)) {
      stop("Column mapping file not found: ", column_file_path)
    }
    # Read and preprocess column mapping
    if (is.null(comment_string)) {
      column_data <- read.table(column_file_path,
        header = FALSE, sep = ":",
        comment.char = "", # This tells R not to treat any character as comment
        stringsAsFactors = FALSE
      ) %>%
        rename(standard = V1, original = V2)
    } else {
      column_data <- read.table(column_file_path,
        header = FALSE, sep = ":",
        comment.char = comment_string,
        stringsAsFactors = FALSE
      ) %>%
        rename(standard = V1, original = V2)
    }
    # Initialize sumstats variable
    sumstats <- NULL
    var_y <- NULL
    sumstats <- load_tsv_region(file_path = sumstat_path, region = region, extract_region_name = extract_region_name, region_name_col = region_name_col)
    
    # To keep a log message
    n_variants <- nrow(sumstats)
    if (n_variants == 0){
        message(paste0("No variants in region ", region, "."))
        return(list(sumstats = sumstats, n = NULL, var_y = NULL))
    } else {
        message(paste0("Region ", region, " include ", n_variants, " in input sumstats."))
    }
    
    # Standardize column names based on mapping
    for (name in colnames(sumstats)) {
      if (name %in% column_data$original) {
        index <- which(column_data$original == name)
        colnames(sumstats)[colnames(sumstats) == name] <- column_data$standard[index]
      }
    }
    if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in%
      colnames(sumstats))) {
      sumstats$z <- sumstats$beta / sumstats$se
    }
    if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
      sumstats$beta <- sumstats$z
      sumstats$se <- 1
    }
    for (col in c("n_sample", "n_case", "n_control")) {
      if (col %in% colnames(sumstats)) {
        sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]],
          na.rm = TRUE
        )
      }
    }
    if (n_sample != 0 && (n_case + n_control) != 0) {
      stop("Please provide sample size, or case number with control number, but not both")
    } else if (n_sample != 0) {
      n <- n_sample
    } else if ((n_case + n_control) != 0) {
      n <- n_case + n_control
      phi <- n_case / n
      var_y <- 1 / (phi * (1 - phi))
    } else {
      if ("n_sample" %in% colnames(sumstats)) {
        n <- median(sumstats$n_sample)
      } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
        n <- median(sumstats$n_case + sumstats$n_control)
        phi <- median(sumstats$n_case / n)
        var_y <- 1 / (phi * (1 - phi))
      } else {
        warning("Sample size and variance of Y could not be determined from the summary statistics.")
        n <- NULL
      }
    }
    return(list(sumstats = sumstats, n = n, var_y = var_y))
  }

  # ---- load_tsv_region (pecotmr R/file_utils.R)
  load_tsv_region <- function(file_path, region = NULL, extract_region_name = NULL, region_name_col = NULL) {
    sumstats <- NULL
    cmd <- NULL
  
    if (!is.null(region)) {
      if (grepl("^chr", region)) {
        region <- strip_chr_prefix(region)
      }
    }
  
    if (grepl("\\.gz$", file_path)) {
      if (is.null(sumstats) || nrow(sumstats) == 0) {
        # Determine the appropriate command based on provided parameters
        if (!is.null(extract_region_name) && !is.null(region) && !is.null(region_name_col)) {
          # Both region and filter specified
          cmd <- paste0(
            "zcat ", file_path, " | head -1 && tabix ", file_path, " ", region,
            " | awk '$", region_name_col, " ~ /", extract_region_name, "/'"
          )
        } else if (!is.null(extract_region_name) && is.null(region) && !is.null(region_name_col)) {
          # Only filter specified, no region
          cmd <- paste0("zcat ", file_path, " | awk '$", region_name_col, " ~ /", extract_region_name, "/'")
        } else if (!is.null(region) && (is.null(region_name_col) || is.null(extract_region_name))) {
          # Only region specified, no filter
          cmd <- paste0("zcat ", file_path, " | head -1 && tabix ", file_path, " ", region)
        } else {
          # Neither region nor filter specified - read gz directly
          cmd <- NULL
        }
  
        sumstats <- tryCatch(
          {
            if (is.null(cmd)) {
              as.data.frame(vroom(file_path, show_col_types = FALSE))
            } else {
              as.data.frame(vroom(pipe(cmd), delim = "\t", show_col_types = FALSE))
            }
          },
          error = function(e) {
            stop("Data read error. Please make sure this gz file is tabix-indexed and the specified filter column exists.")
          }
        )
      }
    } else {
      warning("Not a tabix-indexed gz file, loading the entire dataset.")
      sumstats <- as.data.frame(vroom(file_path, show_col_types = FALSE))
      # Apply filter if specified
      if (!is.null(extract_region_name) && !is.null(region_name_col)) {
        keep_index <- which(str_detect(sumstats[[region_name_col]], extract_region_name))
        sumstats <- sumstats[keep_index, ]
      }
    }
  
    return(sumstats)
  }

  # ---- match_variants_to_keep (pecotmr R/file_utils.R)
  match_variants_to_keep <- function(variant_info, keep_variants_path) {
    keep_raw <- as.data.frame(vroom(keep_variants_path, show_col_types = FALSE))
    # parse_variant_id handles chr prefix, underscore/colon formats, build suffixes,
    # and returns data.frame with integer chrom/pos and character A2/A1
    keep_variants <- parse_variant_id(
      if ("chrom" %in% names(keep_raw) & "pos" %in% names(keep_raw)) keep_raw else keep_raw[[1]]
    )
    vi_chrom <- as.integer(strip_chr_prefix(variant_info$chrom))
    paste0(vi_chrom, ":", variant_info$pos) %in% paste0(keep_variants$chrom, ":", keep_variants$pos)
  }

  # ---- normalize_variant_id (pecotmr R/misc.R)
  normalize_variant_id <- function(ids, chr_prefix = TRUE, convention = NULL) {
    parsed <- parse_variant_id(ids)
    if (!is.null(convention)) {
      format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, convention = convention)
    } else {
      format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, chr_prefix = chr_prefix)
    }
  }

  # ---- NoSNPsError (pecotmr R/file_utils.R)
  NoSNPsError <- function(message) {
    structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
  }

  # ---- order_dedup_regions (pecotmr R/LD.R)
  order_dedup_regions <- function(df) {
    df$chrom <- as.integer(strip_chr_prefix(df$chrom))
    df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
      arrange(chrom, start)
    df
  }

  # ---- parse_region (pecotmr R/misc.R)
  parse_region <- function(region) {
    if (!is.character(region) || length(region) != 1) {
      return(region)
    }
  
    if (!grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", region)) {
      stop("Input string format must be 'chr:start-end'.")
    }
    parts <- str_split(region, "[:-]")[[1]]
    df <- data.frame(
      chrom = strip_chr_prefix(parts[1]),
      start = as.integer(parts[2]),
      end = as.integer(parts[3])
    )
  
    return(df)
  }

  # ---- parse_variant_id (pecotmr R/misc.R)
  parse_variant_id <- function(ids) {
    # Handle data.frame input
    if (is.data.frame(ids)) {
      if (all(c("chrom", "pos", "A2", "A1") %in% names(ids))) {
        # Already has correct column names
      } else if (all(c("chrom", "pos", "A1", "A2") %in% names(ids))) {
        # Has A1/A2 but need to check they're in the right semantic order
        # (A2 = ref, A1 = alt/effect) -- keep as-is since column names are explicit
      } else if (ncol(ids) >= 4) {
        # Assume positional: chrom, pos, A2, A1
        names(ids)[1:4] <- c("chrom", "pos", "A2", "A1")
      }
      # Detect convention from chrom column before converting
      conv <- list(
        has_chr = any(grepl("^chr", as.character(ids$chrom))),
        allele_sep = ":", has_build = FALSE, example = NA_character_
      )
      ids$chrom <- as.integer(strip_chr_prefix(as.character(ids$chrom)))
      ids$pos <- as.integer(ids$pos)
      attr(ids, "convention") <- conv
      return(ids)
    }
  
    # Detect convention before parsing
    convention <- detect_variant_convention(ids)
  
    # Normalize: convert underscores to colons, strip build suffix
    normalized <- gsub("_", ":", ids)
    normalized <- strip_build_suffix(normalized)
  
    # Split into exactly 4 fields using strcapture (vectorized, no list overhead)
    data <- strcapture(
      "^([^:]+):([^:]+):([^:]+):([^:]+)",
      normalized,
      proto = data.frame(chrom = character(), pos = character(),
                         A2 = character(), A1 = character(),
                         stringsAsFactors = FALSE)
    )
  
    data$chrom <- as.integer(strip_chr_prefix(data$chrom))
    data$pos <- as.integer(data$pos)
  
    attr(data, "convention") <- convention
    return(data)
  }

  # ---- process_LD_matrix (pecotmr R/LD.R)
  process_LD_matrix <- function(LD_file_path, snp_file_path = NULL) {
    # Read .cor.xz matrix
    LD_file_con <- xzfile(LD_file_path)
    LD_matrix <- scan(LD_file_con, quiet = TRUE)
    close(LD_file_con)
    LD_matrix <- matrix(LD_matrix, ncol = sqrt(length(LD_matrix)), byrow = TRUE)
  
    # Auto-detect variant metadata file: .bim (PLINK1) or .pvar/.pvar.zst (PLINK2)
    if (is.null(snp_file_path)) {
      candidates <- paste0(LD_file_path, c(".bim", ".pvar", ".pvar.zst"))
      found <- candidates[file.exists(candidates)]
      if (length(found) == 0) stop("No variant file found for: ", LD_file_path,
                                    " (tried .bim, .pvar, .pvar.zst)")
      snp_file_path <- found[1]
    }
  
    # Detect format by extension or header
    is_pvar <- grepl("\\.(pvar|pvar\\.zst)$", snp_file_path)
    if (!is_pvar) {
      # Check if first line starts with #CHROM (pvar header)
      first_line <- readLines(snp_file_path, n = 1)
      is_pvar <- grepl("^#CHROM", first_line)
    }
  
    if (is_pvar) {
      # PLINK2 .pvar format: read via existing read_pvar_text()
      LD_variants <- read_pvar_text(snp_file_path)
      # read_pvar_text returns: chrom, id, pos, A2 (REF), A1 (ALT)
      LD_variants <- LD_variants %>%
        mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom))),
               variants = normalize_variant_id(id)) %>%
        rename(GD = pos)  # reuse pos column for ordering
      LD_variants$GD <- LD_variants$pos <- as.integer(
        sapply(LD_variants$variants, function(v) strsplit(v, ":")[[1]][2]))
    } else {
      # PLINK1 .bim format
      LD_variants <- read.table(snp_file_path)
      col_names <- if (ncol(LD_variants) == 9) {
        c("chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")
      } else if (ncol(LD_variants) == 6) {
        c("chrom", "variants", "GD", "pos", "A1", "A2")
      } else {
        stop("Unexpected number of columns (", ncol(LD_variants), ") in bim file: ", snp_file_path)
      }
      LD_variants <- LD_variants %>%
        setNames(col_names) %>%
        mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom))),
               variants = normalize_variant_id(variants))
    }
  
    # Label and symmetrize the matrix
    colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
    if (all(LD_matrix[lower.tri(LD_matrix)] == 0)) {
      LD_matrix[lower.tri(LD_matrix)] <- t(LD_matrix)[lower.tri(LD_matrix)]
    } else {
      LD_matrix[upper.tri(LD_matrix)] <- t(LD_matrix)[upper.tri(LD_matrix)]
    }
  
    # Order variants by genomic position
    pos_order <- order(sapply(LD_variants$variants, function(v) as.integer(strsplit(v, ":")[[1]][2])))
    LD_variants <- LD_variants[pos_order, ]
    LD_matrix <- LD_matrix[LD_variants$variants, LD_variants$variants]
  
    list(LD_matrix = LD_matrix, LD_variants = LD_variants)
  }

  # ---- read_afreq (pecotmr R/file_utils.R)
  read_afreq <- function(prefix) {
    afreq_zst <- paste0(prefix, ".afreq.zst")
    afreq_plain <- paste0(prefix, ".afreq")
    if (file.exists(afreq_zst)) {
      if (Sys.which("zstd") == "") stop("zstd CLI is required to read .afreq.zst files")
      af <- as.data.frame(vroom(pipe(paste0("zstd -dcq ", shQuote(afreq_zst))),
                                delim = "\t", show_col_types = FALSE))
    } else if (file.exists(afreq_plain)) {
      af <- as.data.frame(vroom(afreq_plain, delim = "\t", show_col_types = FALSE))
    } else {
      return(NULL)
    }
    # PLINK2 .afreq: REF = A2, ALT = A1, ALT_FREQS = A1 (effect allele) frequency
    af <- rename(af,
      "chrom" = "#CHROM", "id" = "ID",
      "A2" = "REF", "A1" = "ALT",
      "alt_freq" = "ALT_FREQS", "obs_ct" = "OBS_CT"
    )
    cols <- c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct")
    # Stochastic genotype .afreq includes U_MIN/U_MAX for exact min-max inversion
    if ("U_MIN" %in% colnames(af)) {
      af <- rename(af, "u_min" = "U_MIN", "u_max" = "U_MAX")
      cols <- c(cols, "u_min", "u_max")
    }
    af <- select(af, all_of(cols))
    return(af)
  }

  # ---- read_pvar_text (pecotmr R/file_utils.R)
  read_pvar_text <- function(pvar_path) {
    # For .zst files, decompress to temp if archive package unavailable
    if (grepl("\\.zst$", pvar_path) && !requireNamespace("archive", quietly = TRUE)) {
      tmp <- tempfile(fileext = ".pvar")
      on.exit(unlink(tmp), add = TRUE)
      ret <- system2("zstd", c("-dq", shQuote(pvar_path), "-o", shQuote(tmp)))
      if (ret != 0) stop("Failed to decompress ", pvar_path, ". Install the 'archive' R package or 'zstd' CLI.")
      pvar_path <- tmp
    }
    df <- as.data.frame(vroom(pvar_path, comment = "##", show_col_types = FALSE))
    colnames(df)[1] <- "chrom"
    data.frame(
      chrom = as.character(df$chrom),
      id    = as.character(df$ID),
      pos   = as.integer(df$POS),
      A2    = as.character(df$REF),
      A1    = as.character(df$ALT),
      stringsAsFactors = FALSE
    )
  }

  # ---- region_to_df (pecotmr R/misc.R)
  region_to_df <- function(ld_region_id, colnames = c("chrom", "start", "end")) {
    region_of_interest <- as.data.frame(do.call(rbind, lapply(strsplit(ld_region_id, "[_:-]"), function(x) as.integer(strip_chr_prefix(x)))))
    colnames(region_of_interest) <- colnames
    return(region_of_interest)
  }

  # ---- resolve_ld_source (pecotmr R/LD.R)
  resolve_ld_source <- function(path) {
    if (!file.exists(path)) {
      stop("LD metadata file not found: ", path,
           "\n  Expected: a TSV file with columns chrom, start, end, path.")
    }
  
    # Peek at first row to determine underlying data type
    meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
    if (ncol(meta) < 4) stop("LD metadata file must have at least 4 columns (chrom, start, end, path): ", path)
    colnames(meta)[1:4] <- c("chrom", "start", "end", "path")
    raw_path <- gsub(",.*$", "", meta$path[1])  # strip comma-separated bim path
    resolved <- file.path(dirname(path), raw_path)
  
    if (has_plink2_files(resolved)) return(list(type = "plink2", data_path = resolved, meta_path = path))
    if (has_plink1_files(resolved)) return(list(type = "plink1", data_path = resolved, meta_path = path))
  
    # Pre-computed .cor.xz blocks — verify not using 0:0 sentinel
    if (!is.na(meta$start) && !is.na(meta$end) && meta$start == 0 && meta$end == 0) {
      stop("Metadata has start=0, end=0 but path does not resolve to PLINK files: ", resolved,
           "\n  The 0:0 sentinel is only valid for whole-chromosome PLINK genotype files.")
    }
  
    list(type = "precomputed", meta_path = path)
  }

  # ---- resolve_plink_prefix_for_region (pecotmr R/LD.R)
  resolve_plink_prefix_for_region <- function(meta_path, region, source_type) {
    parsed <- parse_region(region)
    meta <- as.data.frame(vroom(meta_path, show_col_types = FALSE))
    colnames(meta) <- c("chrom", "start", "end", "path")
    meta$chrom <- as.integer(strip_chr_prefix(meta$chrom))
    query_chrom <- as.integer(strip_chr_prefix(parsed$chrom))
  
    matching <- meta[meta$chrom == query_chrom, , drop = FALSE]
    if (nrow(matching) == 0) {
      stop("No entry for chromosome ", query_chrom, " in metadata file: ", meta_path)
    }
    raw_path <- gsub(",.*$", "", matching$path[1])
    file.path(dirname(meta_path), raw_path)
  }

  # ---- resolve_plink2_paths (pecotmr R/file_utils.R)
  resolve_plink2_paths <- function(prefix) {
    pgen <- paste0(prefix, ".pgen")
    if (!file.exists(pgen)) {
      stop("PLINK2 .pgen file not found at: ", pgen,
           "\n  Note: .pgen must be uncompressed (plink2 does not compress .pgen).")
    }
    # Prefer plain .pvar (fast, no extra deps); fall back to .pvar.zst
    pvar <- if (file.exists(paste0(prefix, ".pvar"))) {
      paste0(prefix, ".pvar")
    } else if (file.exists(paste0(prefix, ".pvar.zst"))) {
      paste0(prefix, ".pvar.zst")
    } else {
      stop("PLINK2 .pvar[.zst] file not found at prefix: ", prefix)
    }
    psam <- paste0(prefix, ".psam")
    if (!file.exists(psam)) {
      stop("PLINK2 .psam file not found at: ", psam,
           "\n  Note: .psam must be uncompressed (plink2 does not compress .psam).")
    }
    list(pgen = pgen, pvar = pvar, psam = psam)
  }

  # ---- strip_build_suffix (pecotmr R/misc.R)
  strip_build_suffix <- function(x) sub("(:|_)b[0-9]+$", "", x)

  # ---- strip_chr_prefix (pecotmr R/misc.R)
  strip_chr_prefix <- function(x) sub("^chr", "", x)

  # ---- tabix_region (pecotmr R/file_utils.R)
  tabix_region <- function(file, region, tabix_header = "auto", target = "", target_column_index = "") {
    if (!file.exists(file)) {
      stop("Input file does not exist: ", file)
    }
    cmd_output <- tryCatch(
      {
        use_col_names <- if (identical(tabix_header, FALSE)) FALSE else TRUE
        as.data.frame(vroom(pipe(paste0("tabix -h ", file, " ", region)),
                            delim = "\t", col_names = use_col_names, show_col_types = FALSE))
      },
      error = function(e) NULL
    )
  
    if (!is.null(cmd_output) && target != "" && target_column_index != "") {
      cmd_output <- cmd_output %>%
        filter(str_detect(.[[target_column_index]], target))
    } else if (!is.null(cmd_output) && target != "") {
      cmd_output <- cmd_output %>%
        mutate(text = apply(., 1, function(row) paste(row, collapse = "_"))) %>%
        filter(str_detect(text, target)) %>%
        select(-text)
    }
  
    if (is.null(cmd_output) || nrow(cmd_output) == 0) {
      return(tibble())
    }

    # qQTLR fix: with few rows vroom guesses an all-"T"/"F" allele column as
    # logical, so single-variant queries lost T/F alleles in allele_qc().
    for (allele_col in intersect(c("A1", "A2"), names(cmd_output))) {
      if (is.logical(cmd_output[[allele_col]])) {
        cmd_output[[allele_col]] <- ifelse(cmd_output[[allele_col]], "T", "F")
      }
    }

    cmd_output %>%
      as_tibble() %>%
      mutate(
        !!names(.)[1] := as.character(.[[1]]),
        !!names(.)[2] := as.numeric(.[[2]])
      )
  }

  # ---- twas_analysis (pecotmr R/twas.R)
  twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
    # Extract gwas_sumstats
    gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs, gwas_sumstats_db$variant_id), ]
    # Validate that the GWAS subset is not empty
    if (nrow(gwas_sumstats_subset) == 0 | all(is.na(gwas_sumstats_subset))) {
      warning("No GWAS summary statistics found for the specified variants.")
      return(NULL)
    }      
    # Check if extract_variants_objs are in the rownames of LD_matrix
    valid_indices <- extract_variants_objs %in% rownames(LD_matrix)
    if (!any(valid_indices)) {
      warning("None of the specified variants are present in the LD matrix. Skipping this context.")
      return(NULL)
    }    
    # Extract only the valid indices from extract_variants_objs
    valid_variants_objs <- extract_variants_objs[valid_indices]
    # Extract LD_matrix subset using valid indices
    LD_matrix_subset <- LD_matrix[valid_variants_objs, valid_variants_objs]
    # Extract weight matrix subset using valid indices
    weights_matrix <- weights_matrix[valid_variants_objs, , drop = FALSE]
    # Caculate the z score and pvalue of each gene
    twas_z_pval <- apply(
      as.matrix(weights_matrix), 2,
      function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset)
    )
    return(twas_z_pval)
  }

  # ---- twas_z (pecotmr R/twas.R)
  twas_z <- function(weights, z, R = NULL, X = NULL) {
    # Check that weights and z-scores have the same length
    if (length(weights) != length(z)) {
      stop("Weights and z-scores must have the same length.")
    }
  
    if (is.null(R)) R <- compute_LD(X)
  
    stat <- t(weights) %*% z
    denom <- t(weights) %*% R %*% weights
    zscore <- stat / sqrt(denom)
    pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)
  
    return(list(z = zscore, pval = pval))
  }

  # ---- validate_selected_region (pecotmr R/LD.R)
  validate_selected_region <- function(start_row, end_row, region_start, region_end) {
    if (start_row$start > region_start || end_row$end < region_end) {
      stop("Region ", region_start, "-", region_end, " is not fully covered by the LD metadata ",
           "(available: ", start_row$start, "-", end_row$end, ").")
    }
  }

  # ---- variant_id_to_df (pecotmr R/misc.R)
  variant_id_to_df <- function(variant_id) {
    parse_variant_id(variant_id)
  }

  invisible(NULL)
}

# Fetch a function from the legacy pecotmr TWAS environment.
.legacy_twas <- function(name) get(name, envir = .pecotmr_legacy, inherits = FALSE)

# Imports of pecotmr 4843e8bd (its NAMESPACE): whole packages, then importFrom()
# entries. .onLoad rebuilds them as the parent of the legacy environment so
# the vendored code resolves names (e.g. dplyr::filter, not stats::filter) as it did in pecotmr.
.pecotmr_legacy_imports <- list(
  whole = c("dplyr", "tibble", "tidyr"),
  from = list(
    coloc = "coloc.bf_bf",
    doFuture = "registerDoFuture",
    dplyr = c("across", "all_of", "any_of", "arrange", "as_tibble", "bind_rows", "case_when", "distinct", "everything", "filter", "full_join", "group_by", "if_else", "inner_join", "intersect", "left_join", "mutate", "n", "pull", "rename", "row_number", "select", "slice", "summarise", "ungroup"),
    furrr = c("furrr_options", "future_map"),
    future = c("availableCores", "multicore", "plan"),
    IRanges = c("end", "findOverlaps", "IRanges", "reduce", "start"),
    magrittr = "%>%",
    matrixStats = c("colSds", "colVars"),
    purrr = c("compact", "exec", "map", "map_int", "map_lgl", "map2", "pmap"),
    readr = c("col_character", "col_guess", "col_integer", "cols", "parse_number", "read_delim"),
    rlang = "!!!",
    S4Vectors = c("queryHits", "subjectHits"),
    stats = c("as.dist", "coef", "cor", "cutree", "hclust", "lm.fit", "pchisq", "pnorm", "predict", "sd", "setNames", "var"),
    stringr = c("str_detect", "str_remove", "str_replace", "str_split"),
    susieR = c("coef.susie", "get_cs_correlation", "mr.ash", "mr.ash.rss", "susie", "susie_get_cs", "susie_rss", "univariate_regression"),
    tibble = c("as_tibble", "tibble"),
    tidyr = c("replace_na", "separate"),
    tools = "file_path_sans_ext",
    utils = c("read.table", "tail"),
    vctrs = "vec_duplicate_detect",
    vroom = "vroom"
  )
)
