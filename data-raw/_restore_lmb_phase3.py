"""Restore glmbayes R/lmb.R from legacy snapshot (Phase 3)."""
from pathlib import Path

legacy = Path(r"C:\Rpackages\glmbayes\legacy_code\core_migration_phase2\R\lmb.R")
out = Path(r"C:\Rpackages\glmbayes\R\lmb.R")
lines = legacy.read_text(encoding="utf-8").splitlines(keepends=True)

# 1-based inclusive slices to keep
keep = (
    list(range(1, 280))  # roxygen + lmb()
    + list(range(546, 564))  # print.lmb
    + list(range(566, 824))  # .uni_lmb
    + list(range(857, 981))  # .mlmb_* helpers (not multi_prior_setup)
    + list(range(1087, 1134))  # normalize + validate helpers
)

body = "".join(lines[i - 1] for i in keep)
body = body.replace(
    "  sim <- rlmb(",
    "  sim <- glmbayesCore::rlmb(",
    1,
)
body = body.replace(
    '      cat("No simbounds returned in sim.\\n")',
    '      warning("No simbounds returned in sim.", call. = FALSE)',
    1,
)
header = (
    "## Phase 3: formula layer in glmbayes; sampling via glmbayesCore::rlmb().\n"
    "## Portions follow stats::lm(); see inst/COPYRIGHTS.\n\n"
)
out.write_text(header + body, encoding="utf-8")
print("Wrote", out, "bytes", out.stat().st_size)
