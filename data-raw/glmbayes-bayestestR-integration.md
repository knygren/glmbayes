# bayestestR / insight Integration for glmbayes

**Status:** Planning document — not yet implemented
**Scope:** `glmb`, `lmb` object classes
**Owner:** Kjell Nygren
**Related roadmap item:** "bayestestR integration" (glmbayes GitHub README, open TODO)

---

## 1. Goal

Enable `glmb`/`lmb` fitted objects to work directly with `bayestestR`
reporting functions — `hdi()`, `eti()`, `describe_posterior()`,
`point_estimate()`, `p_direction()`, `rope()`, etc. — without requiring
users to manually extract draws first.

Target end state:

```r
fit <- glmb(y ~ x1 + x2, family = poisson(), pfamily = ..., data = df)

hdi(fit)
describe_posterior(fit)
p_direction(fit)
rope(fit)
```

---

## 2. Key mechanism: this is an `insight` problem, not a `bayestestR` problem

Most of bayestestR's *model-object* support (as opposed to its
`numeric`/`data.frame` methods, which already work today on any raw
draws matrix) is not implemented in bayestestR itself. bayestestR calls
generic accessor functions from the **`insight`** package:

- `insight::get_parameters()`
- `insight::find_parameters()`
- `insight::model_info()`
- `insight::get_priors()`
- `insight::n_parameters()`
- `insight::clean_parameters()`

If `insight` can extract draws and metadata from a `glmb`/`lmb` object,
**bayestestR support follows almost for free.**

Because S3 dispatch is based on class, not package ownership, these
methods can live inside `glmbayes` itself — there is no need to submit
a PR to `insight` or `bayestestR`. Just implement e.g.
`get_parameters.glmb <- function(x, ...) {...}` in glmbayes, export it,
and register it via `S3method(insight::get_parameters, glmb)` in
`NAMESPACE`.

---

## 3. Current object structure (from `lmb.R` / `glmb.R`)

The fitted object (`outlist`) already contains most of what's needed.
Relevant fields:

| Field | Description | Relevant to |
|---|---|---|
| `coefficients` | draws matrix, **iterations × parameters** | `get_parameters()` — this is the big one, already in the right shape |
| `coef.means` | posterior mean per parameter | sanity check / `point_estimate()` |
| `coef.mode` | posterior mode | MAP comparison |
| `dispersion` | dispersion draws | needs to be treated as its own "parameter class" |
| `fitted.values`, `residuals`, `linear.predictors` | — | not needed for bayestestR, useful for `pp_check`-style diagnostics later |
| `deviance`, `pD`, `Dbar`, `Dthetabar`, `DIC` | model fit diagnostics | could map to `bayestestR`/`performance` diagnostic reporting eventually, not priority |
| `Prior`, `pfamily` | prior specification | `get_priors()` |
| `family`, `formula`, `terms`, `model`, `call` | standard model metadata | `model_info()`, `find_parameters()` |
| `class(outlist)` | `c("lmb","glmb","glm","lm")` | dispatch — see gotcha in §6 |

**Good news:** `coefficients` is already an iterations × parameters
matrix, which is exactly the shape `insight::get_parameters()` is
expected to return (as a data frame). This is the single most
important piece of plumbing and it already exists.

---

## 4. Suggested order of attack

### Phase 1 — Minimum viable integration (unlocks most bayestestR functions)

Implement, in a new file e.g. `R/insight-methods.R`:

1. **`model_info.glmb <- function(x, ...)`**
   - Must return a list with (at minimum) `is_bayesian = TRUE`,
     plus family/link info pulled from `x$family` or `x$famfunc`.
   - This is the flag that tells `insight`/`bayestestR` to treat the
     object as a posterior, not a frequentist fit. Without this,
     dispatch could silently fall through to `glm`/`lm` semantics.

2. **`get_parameters.glmb <- function(x, ...)`**
   - Return `as.data.frame(x$coefficients)`.
   - Column names should match what `find_parameters()` will report
     (Phase 2) — worth deciding parameter naming convention now so
     you don't rename later.
   - Decide: does `dispersion` get included as a "parameter" here, or
     handled as a separate component? (See Phase 2, "component"
     concept.)

3. Write `lmb.glmb <- function(x, ...) UseMethod(...)` is not needed —
   just confirm `get_parameters.lmb` is *not* separately required if
   `lmb` objects inherit `glmb` in the class vector (they do:
   `c("lmb","glmb","glm","lm")`), so `get_parameters.glmb` will be
   found for both. **Verify this dispatch explicitly with a test**
   rather than assuming — see §7.

**Once Phase 1 is done, these should already work:**
`hdi()`, `eti()`, `spi()`, `bci()`, `point_estimate()`, `map_estimate()`,
`p_direction()`, `rope()`, `p_significance()`, `p_map()`, `p_rope()`.

---

### Phase 2 — Clean parameter grouping

4. **`find_parameters.glmb <- function(x, effects = "all", ...)`**
   - Return a named list separating parameter groups — at minimum
     `conditional` (the regression coefficients) vs. `dispersion` (or
     whatever `insight`'s convention calls a nuisance/auxiliary
     parameter — check `find_parameters.brmsfit` or
     `find_parameters.stanreg` source for the expected list-naming
     convention before inventing your own).
   - This is what makes `describe_posterior(fit)` produce a report
     that clearly separates regression coefficients from dispersion,
     instead of dumping them into one flat table.

5. **`n_parameters.glmb`** — trivial wrapper, likely just
   `ncol(x$coefficients)`. Include only if `insight` doesn't already
   derive it correctly by default (test before writing).

6. **`clean_parameters.glmb`** — optional; mainly affects pretty-
   printing/grouping in `describe_posterior()` output tables. Lower
   priority — implement only if default output looks messy.

---

### Phase 3 — Priors (optional, needed for `check_prior()`/`describe_prior()`)

7. **`get_priors.glmb <- function(x, ...)`**
   - Needs to translate `x$Prior` / `x$pfamily` into whatever
     structure `insight::get_priors()` expects (typically a data
     frame with distribution family, location, scale per parameter).
   - This is more work than Phases 1–2 because your prior
     specification system (`pfamily`, conjugate GLM priors, etc.) is
     bespoke and doesn't map 1:1 onto `rstanarm`/`brms` prior syntax.
   - **Recommendation: defer until Phase 1–2 are working and tested.**
     Not needed for the core reporting functions (`hdi`, `rope`,
     `p_direction`, etc.) — only for prior-sensitivity-specific
     functions (`check_prior()`, `describe_prior()`,
     `sensitivity_to_prior()`).

---

### Explicitly out of scope / low priority

- **`bayesfactor_models()`, `bayesfactor_inclusion()`,
  `bayesfactor_restricted()`** — these need either bridge sampling or
  Savage-Dickey density ratios, which assume access to a normalized
  marginal likelihood or a specific prior/posterior density
  comparison. Your iid accept-reject sampler doesn't naturally
  produce this. Possible future direction, not a Phase 1–3 concern.
- **`si()` (Support Intervals)** — likelihood-ratio based; would need
  the likelihood function, not just draws. Doable later given your
  subgradient/envelope machinery already touches the likelihood, but
  not part of this pass.

---

## 5. Naming/shape decisions to make before writing code

These are worth deciding *before* implementation so Phase 1 and
Phase 2 don't require renaming:

- [ ] Should `dispersion` draws be appended as an extra column in the
      same `get_parameters()` data frame, or excluded from Phase 1
      and only introduced properly in Phase 2 via `find_parameters()`
      grouping? (Recommendation: include in Phase 1 as a plain column
      for simplicity, then formalize the split in Phase 2.)
- [ ] Parameter naming convention — should match `coef.means`/
      `coefficients` column names already used elsewhere in the
      package, to avoid confusing users cross-referencing `summary()`
      output against `describe_posterior()` output.
- [ ] Decide how `rlmb`/two-block Gibbs sampler output (once merged)
      will map to this same interface — ideally the same
      `get_parameters.glmb`/`.lmb` methods should work unchanged if
      the returned object also stores draws in a `coefficients`-style
      matrix. Worth checking this now so the two-block sampler and
      the iid sampler don't end up needing separate insight methods.

---

## 6. Gotcha: class vector and dispatch order

Class vector is currently:

```r
class(outlist) <- c(outlist$class, "lmb", "glmb", "glm", "lm")
```

Because `insight` and other packages *already* define
`get_parameters.glm`, `model_info.glm`, etc. (treating them as
frequentist fits), there's a real risk that generic dispatch finds
the wrong method if your new `.glmb`/`.lmb` methods aren't registered
correctly, or if some intermediate code calls `NextMethod()`
unexpectedly.

**Action item:** after writing `model_info.glmb`, explicitly test:

```r
insight::model_info(fit)$is_bayesian  # must be TRUE, not fall through to glm
```

Do this before writing anything else — it's the cheapest possible
smoke test that dispatch is wired correctly.

---

## 7. Testing plan

1. Fit a small `lmb()` model and a small `glmb()` model (reuse
   existing example datasets — Boston, Dobson RCT example, or
   menarche — already used in your vignettes).
2. Confirm dispatch:
   ```r
   insight::model_info(fit)$is_bayesian
   insight::get_parameters(fit)  # should be draws x params data.frame
   ```
3. Run each bayestestR function against the fitted object directly
   and compare output against running the same function on
   `as.data.frame(fit$coefficients)` manually — results should match
   exactly, since Phase 1 is just a thin wrapper around the same
   matrix.
4. Once Phase 2 is in, confirm `describe_posterior(fit)` shows
   dispersion separated from regression coefficients in the printed
   table.
5. Add these as formal package tests (`testthat`) under
   `tests/testthat/test-insight-methods.R` so CRAN checks catch
   regressions if `coefficients` structure changes later (e.g. when
   the two-block Gibbs sampler output is merged).

---

## 8. Suggested file/NAMESPACE additions

- New file: `R/insight-methods.R` — houses all `.glmb`/`.lmb` methods
  from Phases 1–3.
- `NAMESPACE` additions (roxygen: `@method`, `@export`, plus explicit
  `S3method()` registration for the `insight` generics since they're
  not in `Depends`):
  ```r
  #' @method get_parameters glmb
  #' @export
  get_parameters.glmb <- function(x, ...) { ... }
  ```
- `DESCRIPTION`: add `insight` to `Suggests` (already likely there
  given `bayestestR`'s own dependency chain) — confirm whether it
  needs to be `Imports` instead if methods are exported
  unconditionally, or keep as `Suggests` with methods registered
  conditionally via `.onLoad()` + `requireNamespace("insight")` if you
  want to avoid a hard dependency.

---

## 9. Summary of priority order

1. `model_info.glmb` + `get_parameters.glmb` → unlocks HDI/ETI/SPI/BCI,
   point estimates, pd, ROPE, p_significance immediately.
2. `find_parameters.glmb` (+ `n_parameters`, `clean_parameters` as
   needed) → clean grouped output in `describe_posterior()`.
3. `get_priors.glmb` → only if prior-sensitivity reporting
   (`check_prior`, `describe_prior`) is wanted.
4. Bayes factors / support intervals → deferred, likely a separate
   future project given the methodological gap (no marginal
   likelihood from iid accept-reject sampling).
