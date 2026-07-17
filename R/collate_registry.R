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
#'     `core_fn` on every call, overriding its defaults (e.g. forcing a
#'     specific `sbs_manipulation` for `compute_sbs_results_core`). Empty
#'     for every function currently registered -- each one's own defaults
#'     already match what `collate_results()` used.}
#'   \item{batchable}{Whether `collate_results_db()` should split a large
#'     "missing" set into `batch_size`-sized chunks for this function.
#'     Defaults to `TRUE` (via \code{isTRUE()} when absent) for functions
#'     like `compute_ssm_results_core()`, whose per-sample file reads make
#'     memory genuinely scale with how many samples are requested at once.
#'     Set explicitly to `FALSE` for functions whose cost is dominated by
#'     a single, whole-cohort shared read independent of scope size (e.g.
#'     `compute_curated_sv_results_core()`'s curated `.tsv` files,
#'     `compute_sv_results_core()`'s `get_combined_sv()`/`annotate_sv()`) --
#'     chunking those just repeats the same expensive read once per batch
#'     for no benefit, since the read isn't scoped to the batch at all.}
#'   \item{batch_size}{Optional override of `collate_results_db()`'s own
#'     `batch_size` argument, for this function only. Ignored when
#'     `batchable` is `FALSE`. Use this for a function whose per-sample
#'     cost is high enough that even the global default is too much
#'     memory/work per call (e.g. `ssm_results`, via
#'     `get_ssm_by_samples()`).}
#' }
#'
#' @keywords internal
#' @noRd
collate_registry <- list(
  ssm_results = list(
    core_fn = "compute_ssm_results_core",
    metadata_arg = "sample_table",
    extra_args = list(),
    batch_size = 20
  ),
  lymphgen = list(
    core_fn = "compute_lymphgen_core",
    metadata_arg = "these_samples_metadata",
    extra_args = list()
  ),
  dlbclass = list(
    core_fn = "compute_dlbclass_core",
    metadata_arg = "sample_table",
    extra_args = list()
  ),
  qc_results = list(
    core_fn = "compute_qc_results_core",
    metadata_arg = "sample_table",
    extra_args = list()
  ),
  ashm_results = list(
    core_fn = "compute_ashm_results_core",
    metadata_arg = "sample_table",
    extra_args = list()
  ),
  curated_sv_results = list(
    core_fn = "compute_curated_sv_results_core",
    metadata_arg = "sample_table",
    extra_args = list(),
    batchable = FALSE
  ),
  sv_results = list(
    core_fn = "compute_sv_results_core",
    metadata_arg = "sample_table",
    extra_args = list(),
    batchable = FALSE
  ),
  sbs_results = list(
    core_fn = "compute_sbs_results_core",
    metadata_arg = "sample_table",
    extra_args = list()
  )
  # battenberg_purity and csr_results are deliberately not yet registered --
  # pending the user's own independent verification that these
  # revived-from-dormant functions still work correctly against current
  # data (see collate_battenberg_purity.R / collate_csr_results.R).
)
