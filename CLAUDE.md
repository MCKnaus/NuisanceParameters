# CLAUDE.md

Guidance for Claude Code working in the **NuisanceParameters** R package. Read this at the start of every session — it captures everything you need to be productive without re-reading the codebase.

---

## What this package is and who it's for

**NuisanceParameters** estimates conditional expectations (nuisance parameters) used in causal inference under the double/debiased machine learning (DML) framework of Chernozhukov et al. (2018). Targets it can produce:

- `Y.hat`  ≈ `E[Y | X]`
- `Y.hat.d` ≈ `E[Y | D = d, X]`
- `Y.hat.z` ≈ `E[Y | Z = z, X]`
- `D.hat`  ≈ `P[D | X]`     (binary or multivalued treatment)
- `D.hat.z` ≈ `P[D | Z = z, X]`
- `Z.hat`  ≈ `P[Z | X]`

Estimates are obtained by cross-fitted ensemble (stacking) over a wide library of base learners. A defining feature: the package is designed so trained models expose smoother matrices `S` (such that `Y_hat = S Y`), enabling the **OutcomeWeights** companion package to recover the linear outcome weights `omega_i` behind each estimator (see [Knaus, 2024](https://arxiv.org/abs/2411.11559)).

**Audience.** Applied econometricians and causal-inference researchers who want flexible ML-based nuisance estimates that plug directly into DML target-parameter estimators (e.g., `causalDML::DML_aipw`). The package is one of a trilogy: `NuisanceParameters` → `MLeffects` → `OutcomeWeights`.

The current version supports binary and multivalued treatments. Continuous treatments are *not* yet supported.

---

## Architecture at a glance

Main call chain:

```
nuisance_parameters()                 # entry point (R/nuisance_parameters.R)
  └── nuisance_cf()                   # cross-fitting loop per nuisance parameter
        ├── ensemble_short()          # short-stacking branch  (R/ensemble.R)
        └── ensemble()                # standard-stacking branch (R/ensemble.R)
              └── ensemble_core()     # fits the M base learners on a training split
                    └── <method>_fit() + predict.<method>_fit()    # R/ml_wrapper.R
                          (or ovo_fit / ovr_fit for multivalued treatments)
```

Key source files:

- `R/nuisance_parameters.R` — entry point; orchestrates cross-fitting, stratification, clustering, and shape of the returned object.
- `R/ensemble.R` — short-stacking (cross-fitting-based weights) and standard-stacking (extra CV layer for weights).
- `R/ml_wrapper.R` — ~1300 lines; exactly one `*_fit()` constructor + one `predict.*_fit()` method per learner. **Adding a new learner means adding both** and registering it in `create_method()`'s `method = c(...)` choices.
- `R/multiclass_OvO_OvR.R` — One-vs-One and One-vs-Rest wrappers for multivalued treatments.
- `R/utils.R` — `prep_cf_matrix()`, `prep_indicator_mat()`, `ens_weights_maker()`, `add_nupa()`, `create_method()`, `process_methods()`, plot methods for ensemble weights, progress bar setup.
- `R/utils_cleandata.R` — `design_matrix()` (interactions/polynomials/logs) and `data_screen()` (correlation-based screening).
- `R/tune_learners.r` — full-sample / on-the-fold hyperparameter tuning hooks (currently for `forest_grf`, `xgboost`, `ranger`).
- `R/data.R` — roxygen docs for bundled datasets (e.g., `pa_reemployment`).
- `R/zzz.R` — `utils::globalVariables()` declarations to silence `R CMD check` NOTEs.

Exported user API (see `NAMESPACE`): `nuisance_parameters`, `create_method`, `add_nupa`, `design_matrix`, `data_screen`, plus `plot.ens_weights_short` / `plot.ens_weights_stand` and the `predict.*` S3 methods.

### Stacking modes

- **Short stacking** (`stacking = "short"`): one set of weights per nuisance parameter, fit by NNLS (or other) on cross-fitted base-learner predictions. Cheap: cost ≈ F × M.
- **Standard stacking** (`stacking = <integer ≥ 2>`): an inner CV loop inside each cross-fitting fold determines fold-specific weights. Cost ≈ F × V × M. Returns one weight vector per fold.

### Weight-combination methods (`ensemble_type`)

`"nnls"` (default), `"bfgs"` (multivalued only — softmax-constrained Brier-loss optimisation), `"singlebest"`, `"ols"` (continuous only — falls back to `nnls` for multivalued), `"average"`.

### Method-to-task compatibility

`process_methods()` automatically drops classification learners from continuous-outcome NuPas and regression learners from multivalued-treatment NuPas. When extending, keep the `classif_method` / `regr_method` lists in that function in sync with new learners.

---

## Documentation style

All public and internal functions are documented with **roxygen2** (`#'` blocks above each function). Don't edit `NAMESPACE` or files in `man/` by hand — they're regenerated.

**Voice.** Concise, declarative, third-person ("Estimates …", "Generates …", "Creates …"). One-sentence summary on the first line, then a blank `#'`, then a paragraph of detail.

**Conventions to follow when writing new roxygen blocks:**

- First line is a short title; second paragraph elaborates.
- `@param` descriptions start with the type (e.g., "Logical.", "Integer.", "Numeric matrix of …", "Character vector specifying …"), then describe purpose, then defaults if relevant.
- Use `\code{}` for code/identifiers, `\link{}` to cross-reference functions, `\pkg{}` for package names, `\eqn{}` / `\deqn{}` for math.
- For enumerated argument options, use `\describe{ \item{\code{"opt"}}{description} }` (see `ensemble_type` docs in `nuisance_parameters.R` for the canonical example — copy that block when adding a new option-style argument).
- Document the return shape under `@return` with `\describe{}/\item{}` or `\itemize{}` when it's a list.
- Examples: wrap heavy ones in `\dontrun{}` or `\donttest{}` and gate optional dependencies with `if (requireNamespace("pkg", quietly = TRUE))`.
- Mark non-exported helpers with `@keywords internal`. Mark exported S3 methods with both `@method <gen> <class>` and `@exportS3Method` (or `@export` for plot methods).
- Run `devtools::document()` after touching any roxygen comment. Never edit `NAMESPACE` or `man/*.Rd` directly.

**Vignettes** (`vignettes/*.Rmd`). YAML header always lists three authors (Knaus, Glaisner, Rakov) and uses `rmarkdown::html_vignette`. Style is conversational and instructional: a brief introduction with a bulleted list of what the vignette covers, then alternating prose and code chunks. Cross-reference other vignettes by their italicised titles (*Stacking*, *Hyperparameter Tuning*, *Multivalued Treatment*, *Saving models*). Inline math uses `$…$`, display math uses `$$…$$`. Reference papers with full markdown links (e.g., `[Knaus (2024)](https://arxiv.org/abs/2411.11559)`). Examples consistently build on the **hdm** `pension` dataset or the bundled `pa_reemployment` dataset.

---

## R coding conventions in this codebase

This is a base-R codebase. Match the existing style:

- **Indentation.** 2 spaces. No tabs.
- **Naming.** `snake_case` for functions, variables, and list elements. Class names follow `<method>_fit` / `ensemble`, `ensemble_short`, `ensemble_core`, `ens.learner`, `NuisanceParameters`. Internal helpers use the same casing as exported ones.
- **No tidyverse in package code.** Use base R: `lapply`, `vapply`, `apply`, `do.call`, `Reduce`, `stats::lm.fit`, `stats::model.matrix`, etc. `ggplot2` is the only pipeline package imported (for the two `plot.*` methods) and it's referenced explicitly with `ggplot2::`.
- **Namespace qualification.** Always namespace-qualify calls into base packages other than `base` and `methods` (`stats::`, `utils::`, `grDevices::`, `Matrix::`). For optional Suggests-only learners, always guard with `requireNamespace("pkg", quietly = TRUE)` and `stop("The 'pkg' package is not installed. … Install it with: install.packages('pkg')", call. = FALSE)`.
- **Argument validation.** Use `match.arg()` for enumerated string arguments. Throw early with informative `stop()` messages. For nuisance-parameter / option subsets, see `nuisance_parameters.R` lines 95–127 as the template.
- **Comments.** Use short `## Section` banners to separate logical blocks inside long functions (e.g., `## Checks`, `## Preps`, `## Initialize objects`, `## Outcome NuPa`). Inline `#` comments explain *why* a non-obvious step is taken — don't restate what the code does.
- **Verbosity control.** Functions take `quiet = TRUE` as the default; user-visible progress comes from the `progress` package via `setup_pb()` / `update_progress()` when `quiet = FALSE`. Use `message()` (not `cat()` or `print()`) for informational output.
- **Adding a learner.**
  1. Write `<method>_fit(X, Y, arguments = list(), ...)` returning a `structure(list(...), class = "<method>_fit")`.
  2. Write `predict.<method>_fit(object, Xnew = NULL, ...)` returning a numeric vector (continuous/binary) or matrix (multinomial).
  3. Register the method name in the `method = c(...)` choices of `create_method()` and document it under `@details` with the same bullet style as the existing learners.
  4. Add it to `classif_method` or `regr_method` in `process_methods()` so the auto-filter knows where it applies.
  5. If it requires a Suggests-only package, gate with `requireNamespace()` as described above.
  6. If it supports tuning, add a sibling `<method>_tune()` function and the corresponding branch in `tune_learners()`.
  7. Add the `S3method(predict, <method>_fit)` line by re-running `devtools::document()`.

---

## Testing patterns

- `testthat` edition 3 (configured in `DESCRIPTION`).
- One file per topic in `tests/testthat/` (e.g., `test-ensemble.R`, `test-nuisance_parameters.R`, `test-ml_wrapper.R`).
- Synthetic data is generated with `mvtnorm::rmvnorm()` and a fixed `set.seed()`.
- Integration tests in `test-nuisance_parameters.R` exercise the full pipeline: cross-fitting correctness, RMSE-vs-naive thresholds, and timing comparisons (`t_short < t_standard`). Compatibility with `OutcomeWeights::get_smoother_weights` is verified via `S %*% Y ≈ Y_hat` to a tolerance of `1e-5` to `1e-7`.
- Use `skip_if_not_installed("pkg")` for tests that depend on Suggests packages (especially `OutcomeWeights`, `glmnet`, `xgboost`, `grf`).
- Don't treat the test suite as canonical — coverage is uneven and we plan to expand it. When fixing or extending behaviour, add a focused unit test alongside the integration check.

Run a single file with `testthat::test_file("tests/testthat/test-<name>.R")`.

---

## How to validate changes

**Always use the no-vignettes check** when verifying your work — vignettes are slow and pull in heavy Suggests dependencies that aren't always installed:

```r
devtools::check(vignettes = FALSE)
```

Do **not** run a full `devtools::check()` unless explicitly asked. After any roxygen change, also run:

```r
devtools::document()
devtools::test()
```

For interactive iteration use `devtools::load_all()`.

---

## Git workflow — non-negotiable

- **Never push directly to `main`.** Every change goes through a pull request.
- Work on a feature/fix branch; open a PR targeting `main` and let it be reviewed/merged there.
- Don't force-push to `main` and don't run `git push` to `main` from any branch under any circumstance.
- Confirm with the user before any destructive git operation (force pushes on feature branches, `git reset --hard`, branch deletion, etc.).
