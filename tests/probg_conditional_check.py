#!/usr/bin/env python
"""Is the sampled M consistent with the conditional implied by its own genealogy?

Background: with uninformative data the marginal of M must be its prior, but a
2-population 4+4 all-'?' run converges (from both directions -- see
overdispersed_start_check.py) to E[M] ~ 0.46 x the prior mean. Both conditionals
have separately tested exact:

  * p(G | theta, M) matches msprime at the chain's own parameter values
  * p(M | G) matches the analytic truncated Gamma when the tree is frozen

Two exact conditionals cannot produce a wrong marginal, so one of them is not
the conditional of the same joint density. This test decides which, by making
probg report the sufficient statistics it actually uses.

A -DDEBUGMIG build prints, every MIGRATE_PROBG_EVERY calls, one line per
migration parameter:

    PROBG  call  locus  heat  param  n_p  L_p  M_p  sumprob

n_p is the number of migration events attributed to parameter p, L_p the pure
migration opportunity (the coefficient of M_p in the waiting term). So the
conditional probg itself implies is

    p(M_p | G)  =  Gamma(shape n_p + 1, rate L_p + 1/scale) truncated to [0,hi]

Two readouts:

  * PIT: u = F(M_p) should be Uniform(0,1) if the M the chain holds really is a
    draw from the conditional probg scores. Skew locates the disagreement.
  * n_p / (L_p * M_p): the tree kernel puts down migration events at rate M_p
    over opportunity L_p, so this ratio must average 1. A stable departure from
    1 is a rate/opportunity convention mismatch between the tree kernel and
    probg, and its size is the factor to look for in the code.

Usage:
    cd src && make CPPFLAGS="-DPRETTY -DNEWVERSION -D_REENTRANT -DDEBUGMIG" migrate-n
    mv migrate-n migrate-n-debugmig && git checkout migrate-n 2>/dev/null
    python tests/probg_conditional_check.py src/migrate-n-debugmig [scale]
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

MIGRATE = str(Path(sys.argv[1]).resolve()) if len(sys.argv) > 1 else str(
    Path.home() / "src/migrate-5.0.7/src/migrate-n-debugmig")
SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 20000.0
THETA = 0.05
SITES = 10
NPOP = 2
NPERPOP = 4
SAMPLE = 20000
EVERY = 200
SEEDS = (131, 231, 331)

PARM = """\
################################################################################
# Parmfile for Migrate
menu=NO
datatype=SequenceData
datamodel=JC69
freqs-from-data=NO:0.25,0.25,0.25,0.25
infile=infile
progress=NO
print-data=NO
pdf-outfile=NO
heating=NO
replicate=NO
long-chains=1
long-inc=20
long-sample={sample}
burn-in=5000
auto-tune=NO
use-M=YES
custom-migration={{**}}
updatefreq= tree:0.70 parameter:0.30 scaler:0.00
theta=OWN:{theta}
bayes-priors= THETA * * UNIFORMPRIOR: {tlo:.7f} {thi:.7f} 0.0000005
bayes-priors= MIG * * EXPPRIOR: 0.000000 {scale:.6f} {scale:.6f}
random-seed=OWN:{seed}
outfile=outfile
bayes-allfile=NO
"""


def run(seed):
    """Return array of (param, n_p, L_p, M_p) from one instrumented run."""
    with tempfile.TemporaryDirectory() as wd:
        wd = Path(wd)
        with open(wd / "infile", "w") as f:
            f.write("   %i 1 . allmiss\n%i\n" % (NPOP, SITES))
            for p in range(1, NPOP + 1):
                f.write("   %i Pop%i\n" % (NPERPOP, p))
                for i in range(NPERPOP):
                    f.write("%-10s%s\n" % ("p%si%i" % (p, i), "?" * SITES))
        (wd / "parmfile").write_text(PARM.format(
            theta=THETA, tlo=THETA * 0.99998, thi=THETA * 1.00002,
            scale=SCALE, seed=seed, sample=SAMPLE))
        proc = subprocess.run(
            [MIGRATE, "parmfile", "-nomenu"], cwd=wd,
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
            env={**os.environ, "MIGRATE_PROBG_EVERY": str(EVERY)}, check=True)
    rows = []
    for line in proc.stderr.splitlines():
        if not line.startswith("PROBG\t"):
            continue
        f = line.split("\t")
        # PROBG call locus heat param n_p L_p M_p sumprob
        rows.append((float(f[4]), float(f[5]), float(f[6]), float(f[7])))
    return rows


def pit(n, L, M, scale, hi):
    """Truncated-Gamma probability integral transform of M given (n, L)."""
    rate = L + 1.0 / scale
    cdf = stats.gamma.cdf(M, a=n + 1.0, scale=1.0 / rate)
    norm = stats.gamma.cdf(hi, a=n + 1.0, scale=1.0 / rate)
    return np.where(norm > 0, cdf / np.maximum(norm, 1e-300), np.nan)


if __name__ == "__main__":
    rows = []
    for s in SEEDS:
        rows += run(s)
    if not rows:
        sys.exit("no PROBG lines -- is this a -DDEBUGMIG build?")
    A = np.array(rows)
    param, n, L, M = A[:, 0].astype(int), A[:, 1], A[:, 2], A[:, 3]

    print("probg conditional check   prior truncExp(scale %g) on [0,%g]"
          % (SCALE, SCALE))
    print("%d records from %d seeds, every %d probg calls\n"
          % (len(A), len(SEEDS), EVERY))
    hdr = ("param    n_p mean    L_p mean     M_p mean   (n+1)/(L+1/s)"
           "   n/(L*M)   PIT mean    KS p")
    print(hdr)
    print("-" * len(hdr))
    for p in sorted(set(param)):
        k = param == p
        nk, Lk, Mk = n[k], L[k], M[k]
        cond = (nk + 1.0) / (Lk + 1.0 / SCALE)
        ok = (Lk > 0) & (Mk > 0)
        ratio = np.median(nk[ok] / (Lk[ok] * Mk[ok]))
        u = pit(nk, Lk, Mk, SCALE, SCALE)
        u = u[np.isfinite(u)]
        ks = stats.kstest(u, "uniform").pvalue if len(u) > 2 else float("nan")
        print("%5d %11.1f %11.5g %12.1f %14.1f %9.3f %10.3f %8.2g"
              % (p, nk.mean(), Lk.mean(), Mk.mean(), cond.mean(), ratio,
                 u.mean(), ks))

    print()
    print("expected if the tree kernel and probg agree:")
    print("  n/(L*M) ~ 1.0     PIT mean ~ 0.5     KS p not tiny")
    print("a stable n/(L*M) away from 1 is the convention mismatch, and its")
    print("value is the factor to hunt for in the migration rate/opportunity.")
