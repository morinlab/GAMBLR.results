#' @title Registry of `collate_*` functions backed by the SQLite database.
#'
#' @description Deliberately explicit, not name-pattern-discovered: several
#' `collate_*`-named files in this package are dormant or disabled (see
#' `collate_results()`'s own commented-out calls to `collate_csr_results()`
#' and `collate_pga()`, and the fact that `collate_derived_results()`,
#' `collate_battenberg_purity()`, `collate_extra_metadata()`, and
#' `collate_hnrph1_mutations()` aren't called from it at all). Matching
#' `^collate_` by name would silently resurrect those. Add a function here
#' only once it's been deliberately migrated to the core/wrapper split (see
#' `collate_ssm_results.R` for the pattern) and confirmed working.
#'
#' Each entry:
#' \describe{
#'   \item{core_fn}{Name of the "core" function -- takes a sample scope in,
#'     returns only its own new columns, keyed by `sample_id` (and
#'     `seq_type` where applicable). Does not read or write the database
#'     itself; `collate_results_db()` handles that generically for every
#'     registered function.}
#'   \item{metadata_arg}{The core function's parameter name for the sample
#'     scope it's given. Not uniform across functions -- most use
#'     `sample_table`, `collate_lymphgen`'s core uses
#'     `these_samples_metadata` -- so this is recorded explicitly per
#'     function rather than assumed.}
#'   \item{extra_args}{Named list of additional fixed arguments to pass to
#'     `core_fn` on every call (e.g. `collate_sbs_results`'s
#'     `sbs_manipulation`). Empty for functions that need nothing beyond
#'     the sample scope.}
#' }
#'
#' @keywords internal
#' @noRd
collate_registry <- list(
  ssm_results = list(
    core_fn = "compute_ssm_results_core",
    metadata_arg = "sample_table",
    extra_args = list()
  ),
  lymphgen = list(
    core_fn = "compute_lymphgen_core",
    metadata_arg = "these_samples_metadata",
    extra_args = list()
  )
  # Remaining registered-but-not-yet-migrated functions (see
  # docker/RELEASING.md-style planning notes / the CONTRIBUTING.md "Future
  # direction" section): sv_results, curated_sv_results, ashm_results,
  # sbs_results, qc_results, dlbclass, battenberg_purity, csr_results.
  # Each needs the same core/wrapper extraction as ssm_results and
  # lymphgen before being added here.
)
