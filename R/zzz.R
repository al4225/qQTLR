# Lazy-bind unexported or renamed pecotmr helpers into qQTLR's namespace.

compute_qvalues              <- NULL
pval_cauchy                  <- NULL
drop_collinear_columns       <- NULL
parse_variant_id             <- NULL
filter_variants_by_ld_reference <- NULL
ld_prune_by_correlation      <- NULL
ld_clump_by_score            <- NULL
.pecotmr_legacy              <- NULL
.pecotmr_legacy_ld           <- NULL

# Evaluate vendored pecotmr code (defs_fn) in a private environment whose parent
# rebuilds that pecotmr version's NAMESPACE imports (whole packages first, then
# importFrom), so names resolve as they did there regardless of what is attached.
# Packages the vendored functions never call are only Suggests: skip if absent.
.build_legacy_env <- function(defs_fn, imports, parent) {
  legacy_imports <- new.env(parent = parent)
  for (pkg in imports$whole) {
    if (!requireNamespace(pkg, quietly = TRUE)) next
    for (fn in getNamespaceExports(pkg)) {
      assign(fn, getExportedValue(pkg, fn), envir = legacy_imports)
    }
  }
  for (pkg in names(imports$from)) {
    if (!requireNamespace(pkg, quietly = TRUE)) next
    for (fn in imports$from[[pkg]]) {
      assign(fn, getExportedValue(pkg, fn), envir = legacy_imports)
    }
  }
  legacy <- new.env(parent = legacy_imports)
  eval(body(defs_fn), envir = legacy)
  legacy
}

.onLoad <- function(libname, pkgname) {
  ns_self <- asNamespace(pkgname)

  # Helper: bind a renamed pecotmr function under a legacy snake_case name.
  bind_from_pecotmr <- function(new_name, old_name, fallback_msg) {
    fn <- tryCatch(
      utils::getFromNamespace(new_name, "pecotmr"),
      error = function(e) NULL
    )
    if (!is.null(fn)) {
      assign(old_name, fn, envir = ns_self)
    } else {
      msg <- fallback_msg
      assign(old_name, function(...) stop(msg), envir = ns_self)
    }
  }

  # pval_cauchy: Cauchy combination test, renamed to pvalAcat in current pecotmr.
  bind_from_pecotmr("pvalAcat", "pval_cauchy",
    "'pvalAcat' is not available in the installed version of pecotmr")

  # drop_collinear_columns: renamed to dropCollinearColumns in current pecotmr.
  bind_from_pecotmr("dropCollinearColumns", "drop_collinear_columns",
    "'dropCollinearColumns' is not available in the installed version of pecotmr")

  # ld_prune_by_correlation: renamed to ldPruneByCorrelation.
  bind_from_pecotmr("ldPruneByCorrelation", "ld_prune_by_correlation",
    "'ldPruneByCorrelation' is not available in the installed version of pecotmr")

  # ld_clump_by_score: renamed to ldClumpByScore.
  bind_from_pecotmr("ldClumpByScore", "ld_clump_by_score",
    "'ldClumpByScore' is not available in the installed version of pecotmr")

  # Legacy pecotmr TWAS helpers (legacy_pecotmr_twas.R, pecotmr 4843e8bd):
  # pecotmr's harmonizeTwas() is now a deprecation stub, so quantile TWAS uses
  # the vendored code instead.
  assign(".pecotmr_legacy",
         .build_legacy_env(.pecotmr_legacy_twas_defs, .pecotmr_legacy_imports, ns_self),
         envir = ns_self)

  # Legacy pecotmr 0.6.11 LD-reference filter and variant-ID parser
  # (legacy_pecotmr_ldfilter.R). Current filterVariantsByLdReference() only reads
  # LD from genotype files and errors on pre-computed .cor.xz blocks, and current
  # parseVariantId() returns chrom as "chr2" instead of "2".
  legacy_ld <- .build_legacy_env(.pecotmr_legacy_ld_defs, .pecotmr_legacy_ld_imports, ns_self)
  assign(".pecotmr_legacy_ld", legacy_ld, envir = ns_self)
  assign("filter_variants_by_ld_reference",
         get("filterVariantsByLdReference", envir = legacy_ld), envir = ns_self)
  assign("parse_variant_id", get("parseVariantId", envir = legacy_ld), envir = ns_self)

  # compute_qvalues: q-value estimation via qvalue package, with FDR fallback.
  assign("compute_qvalues", function(pvals) {
    tryCatch(
      qvalue::qvalue(pvals)$qvalues,
      error = function(e) p.adjust(pvals, method = "fdr")
    )
  }, envir = ns_self)
}
