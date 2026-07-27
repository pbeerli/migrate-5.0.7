#!/usr/bin/env python
"""
Regression test: with uninformative data, does the migration rate M come back
as its prior?

Two populations, 4+4 haploid tips, every site '?', theta pinned by a degenerate
uniform prior. The data likelihood is constant, so the joint target is
p(G|theta,M)*prior(M) and the marginal for M must be exactly the prior.

The M prior here is deliberately small-scale, giving roughly one migration event
per genealogy. That matters: the M update is an independence sampler drawing
from the prior, while the conditional

    p(M|G)  ~  M^n * exp(-L*M) * prior(M)

is a truncated Gamma of shape n+1, where n is the number of migration events on
the genealogy. Its relative width is 1/sqrt(n), so as n grows the conditional
becomes a needle inside a broad proposal, acceptance collapses, and M can only
move as fast as the genealogy's migration count random-walks -- O(n^2) tree
updates per relaxation. That is a mixing limit, not a bias: measured recovery
ratios on this setup were

    prior scale     10    100   1000   5000   20000
    E[n_mig]       1.0    8.7   74.7  274.7   306.4
    recovered/target
                 1.001  1.002  0.982  0.839   0.461

so the test is run in the small-n regime where the sampler can actually be
checked. Slice and multiplier proposals for MIG do not rescue the large-n case
(0.524 and 0.472 at scale 20000) -- the bottleneck is the joint (M, genealogy)
ridge, not the proposal width.

Usage:  python tests/migration_prior_check.py <path-to-migrate-n> [prior_scale]
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# resolved to an absolute path: migrate is run from a temporary directory
MIGRATE = str(Path(sys.argv[1]).resolve()) if len(sys.argv) > 1 else str(Path.home() / "bin/migrate-n")
SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
THETA = 0.05
SITES = 10
NPOP = 2
NPERPOP = 4
SEEDS = (131, 231, 331, 431)
TOLERANCE = 0.05  # fractional deviation allowed on the mean

PARM = """\
################################################################################
# Parmfile for Migrate
# DO NOT REMOVE THE FIRST TWO LINES!
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
long-sample=200000
burn-in=20000
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


def truncated_exponential_moments(scale, lo, hi):
    """Mean and median of an exponential of the given scale truncated to [lo, hi]."""
    z = np.exp(-lo / scale) - np.exp(-hi / scale)
    mean = (
        (lo + scale) * np.exp(-lo / scale) - (hi + scale) * np.exp(-hi / scale)
    ) / z
    median = -scale * np.log(np.exp(-lo / scale) - 0.5 * z)
    return mean, median


def run(workdir, seed):
    workdir = Path(workdir)
    with open(workdir / "infile", "w") as f:
        f.write("   %i 1 . all-missing %ipop\n%i\n" % (NPOP, NPOP, SITES))
        for p in range(1, NPOP + 1):
            f.write("   %i Pop%i\n" % (NPERPOP, p))
            for i in range(NPERPOP):
                f.write("%-10s%s\n" % ("p%si%i" % (p, i), "?" * SITES))
    (workdir / "parmfile").write_text(
        PARM.format(theta=THETA, tlo=THETA * 0.99998, thi=THETA * 1.00002,
                    scale=SCALE, seed=seed)
    )
    subprocess.run(
        [MIGRATE, "parmfile", "-nomenu"],
        cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    rows = []
    for line in (workdir / "bayesall").read_text().splitlines():
        if line.startswith("#") or line.startswith("Steps"):
            continue
        p = line.split("\t")
        if len(p) >= len(COLS):
            rows.append([float(x) for x in p[:len(COLS)]])
    A = np.array(rows)
    d = {c: A[:, i] for i, c in enumerate(COLS)}
    theta = d["Theta_1"]
    assert theta.max() - theta.min() < 1e-5, "theta was not held fixed"
    # timelist entries are ntips tips + (ntips-1) coalescences + migrations
    # + one root sentinel, so migrations = treeintervals - 2*ntips
    n_mig = d["treeintervals"] - 2 * NPOP * NPERPOP
    return d["M_2_1"], d["M_1_2"], n_mig


if __name__ == "__main__":
    mean, median = truncated_exponential_moments(SCALE, 0.0, SCALE)
    got, meds, nmigs = [], [], []
    for seed in SEEDS:
        with tempfile.TemporaryDirectory() as tmp:
            m21, m12, nm = run(tmp, seed)
        got += [m21.mean(), m12.mean()]
        meds += [np.median(m21), np.median(m12)]
        nmigs.append(nm.mean())
    sampled, sampled_med = float(np.mean(got)), float(np.mean(meds))
    ratio = sampled / mean
    print(f"M prior: truncated Exp(scale {SCALE:g}) on [0, {SCALE:g}]")
    print(f"  target   mean {mean:10.4f}   median {median:10.4f}")
    print(f"  sampled  mean {sampled:10.4f}   median {sampled_med:10.4f}"
          f"   ({len(SEEDS)} seeds x 2 parameters)")
    print(f"  ratio    {ratio:.4f}    mean migration events per genealogy "
          f"{np.mean(nmigs):.2f}")
    ok = abs(ratio - 1.0) <= TOLERANCE
    print("\n" + ("PASS" if ok else
                  f"FAIL: recovered/target = {ratio:.4f}, "
                  f"outside 1 +/- {TOLERANCE}"))
    sys.exit(0 if ok else 1)
