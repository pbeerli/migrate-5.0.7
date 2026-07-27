#!/usr/bin/env python
"""Do exchangeable migration parameters actually behave exchangeably?

In a fully symmetric two-population model -- equal thetas pinned, identical
priors, equal sample sizes -- M_2->1 and M_1->2 are exchangeable. Two things
must therefore hold, and they fail in different ways:

  1. Across chains, neither parameter is favoured. mean log(M21/M12) = 0.
     A violation is an indexing or direction bug.
  2. Within each chain, the fraction of samples with M21 > M12 concentrates
     at 0.5. A violation means the chain is stuck on one side -- it is a
     TRAPPING diagnostic, and it is the one that matters here.

Test 2 is what settled the long-standing "M comes back below its prior"
question. Between-start diagnostics could not: every start converges to the same
wrong answer and longer runs do not move it, which looks like a wrong stationary
distribution but is equally consistent with all starts falling into a trap basin
quickly (see overdispersed_start_check.py). The within-chain statistic
distinguishes them, because trapping shows up directly as chains that never
cross M21 = M12.

Measured, 40 independent seeds, 20000 samples each:

    prior scale      100      1000     20000
    pooled/target  1.0023    0.9990    0.4767
    chains locked    0/16      0/16     11/40
    frac range   .482-.523 .389-.565 .006-1.000
    mean log(M21/M12)                  +0.0101  (t=0.04, p=0.97)

So there is no direction bias at any scale -- only trapping, and only once the
genealogies carry many migration events. The default prior scale here is small
enough that the test passes; run it at 20000 to reproduce the failure.

**Use at least ~40 seeds before believing test 1 in the trapping regime.** Locked
chains are each stuck on one side, so per-seed values are strongly correlated and
a handful of seeds can agree by chance: 8 seeds at scale 20000 gave
mean log(M21/M12) = +0.95 with p = 0.014, against +0.010 and p = 0.97 at 40
seeds. A spurious "reproducible M21/M12 asymmetry" was claimed and retracted
during this investigation for exactly this reason (that time the three runs also
shared a seed). Test 2 is reliable at small seed counts; test 1 is not.

Usage:
    python tests/exchangeability_check.py <path-to-migrate-n> [prior_scale] [n_seeds]

    scale 100 (default) passes in seconds; scale 20000 reproduces the failure
    and takes a couple of minutes at 8 seeds.
"""
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

THETA = 0.05
SITES = 10
NPOP = 2
NPERPOP = 4
SAMPLE = 20000
TOL_MEAN = 0.05      # fractional deviation allowed on the pooled mean
TOL_LOCKED = 0.0     # chains allowed to be locked below 0.1 / above 0.9

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
burn-in=50000
auto-tune=NO
use-M=YES
custom-migration={{**}}
updatefreq= tree:0.70 parameter:0.30 scaler:0.00
theta=OWN:{theta}
bayes-priors= THETA * * UNIFORMPRIOR: {tlo:.7f} {thi:.7f} 0.0000005
bayes-priors= MIG * * EXPPRIOR: 0.000000 {scale:.6f} {scale:.6f}
random-seed=OWN:{seed}
outfile=outfile
bayes-allfile=YES:1:bayesall
"""

COLS = ("Steps Locus Replicate lnPost lnDataL lnPrbGParam lnPrior "
        "treeintervals treelength Theta_1 Theta_2 M_2_1 M_1_2").split()


def run(job):
    """One independent chain. Everything is passed in, so this is spawn-safe."""
    migrate, scale, sample, seed = job
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
            scale=scale, seed=seed, sample=sample))
        subprocess.run([migrate, "parmfile", "-nomenu"], cwd=wd,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=True)
        rows = []
        for line in (wd / "bayesall").read_text().splitlines():
            if line.startswith("#") or line.startswith("Steps"):
                continue
            p = line.split("\t")
            if len(p) >= len(COLS):
                rows.append([float(x) for x in p[:len(COLS)]])
    A = np.array(rows)
    d = {c: A[:, i] for i, c in enumerate(COLS)}
    m21, m12 = d["M_2_1"], d["M_1_2"]
    return (float(m21.mean()), float(m12.mean()),
            float(np.mean(np.log(m21 + 1) - np.log(m12 + 1))),
            float(np.mean(m21 > m12)))


def truncated_exponential_mean(scale, lo, hi):
    z = np.exp(-lo / scale) - np.exp(-hi / scale)
    return ((lo + scale) * np.exp(-lo / scale)
            - (hi + scale) * np.exp(-hi / scale)) / z


if __name__ == "__main__":
    migrate = str(Path(sys.argv[1]).resolve()) if len(sys.argv) > 1 else str(
        Path.home() / "bin/migrate-n")
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    nseed = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    seeds = [101 + 7919 * i for i in range(nseed)]      # distinct, well spread

    jobs = [(migrate, scale, SAMPLE, s) for s in seeds]
    with ProcessPoolExecutor(max_workers=min(10, os.cpu_count() or 1)) as ex:
        res = list(ex.map(run, jobs))

    m21 = np.array([r[0] for r in res])
    m12 = np.array([r[1] for r in res])
    dl = np.array([r[2] for r in res])
    fr = np.array([r[3] for r in res])
    target = truncated_exponential_mean(scale, 0.0, scale)
    ratio = 0.5 * (m21.mean() + m12.mean()) / target
    t = stats.ttest_1samp(dl, 0.0)
    locked = int(((fr < 0.1) | (fr > 0.9)).sum())

    print("%d independent seeds, symmetric 2-pop model, prior mean %.1f\n"
          % (len(res), target))
    print("  mean M_2->1 %10.1f   mean M_1->2 %10.1f" % (m21.mean(), m12.mean()))
    print("  pooled / target %8.4f" % ratio)
    print()
    print("1. systematic asymmetry (indexing bug?)")
    print("   mean log(M21/M12) %+.4f   t %+.2f   p %.3f   seeds M21>M12 %d/%d"
          % (dl.mean(), t.statistic, t.pvalue, int((dl > 0).sum()), len(dl)))
    print()
    print("2. within-chain locking (trapping?)")
    print("   frac(M21>M12): min %.3f  median %.3f  max %.3f"
          % (fr.min(), np.median(fr), fr.max()))
    print("   chains locked (<0.1 or >0.9): %d / %d" % (locked, len(fr)))

    ok_mean = abs(ratio - 1.0) <= TOL_MEAN
    ok_lock = locked <= TOL_LOCKED
    ok_sym = t.pvalue > 0.01
    print()
    if ok_mean and ok_lock and ok_sym:
        print("PASS")
        sys.exit(0)
    if not ok_mean:
        print("FAIL: pooled/target = %.4f, outside 1 +/- %.2f" % (ratio, TOL_MEAN))
    if not ok_lock:
        print("FAIL: %d chains locked on one side of M21 = M12 (trapping)" % locked)
    if not ok_sym:
        print("FAIL: systematic direction asymmetry, p = %.4g" % t.pvalue)
    sys.exit(1)
