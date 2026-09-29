"""Restore glmbayes R/glmb.R from legacy snapshot (Phase 3)."""
from pathlib import Path

legacy = Path(r"C:\Rpackages\glmbayes\legacy_code\core_migration_phase2\R\glmb.R")
out = Path(r"C:\Rpackages\glmbayes\R\glmb.R")
lines = legacy.read_text(encoding="utf-8").splitlines(keepends=True)

keep = list(range(1, 481)) + list(range(484, 504))  # glmb + print.glmb

body = "".join(lines[i - 1] for i in keep)
body = body.replace(
    "    sim<-rglmb(n=n,y=y,x=x,family=family,pfamily=pfamily,offset=offset,",
    "    sim <- glmbayesCore::rglmb(n=n,y=y,x=x,family=family,pfamily=pfamily,offset=offset,",
    1,
)
header = (
    "## Phase 3: formula layer in glmbayes; sampling via glmbayesCore::rglmb().\n"
    "## Portions follow stats::glm(); see inst/COPYRIGHTS.\n\n"
)
out.write_text(header + body, encoding="utf-8")
print("Wrote", out, "bytes", out.stat().st_size)
