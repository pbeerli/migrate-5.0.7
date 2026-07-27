#!/usr/bin/env python
"""
Regression test: does migrate's genealogy update actually sample the Kingman
coalescent?

Setup: one population, n haploid tips, every site '?', theta pinned by a
degenerate uniform prior, tree updates only (no parameter / scaler moves).
The data likelihood is then constant, so acceptlike() always returns TRUE and
the recorded genealogies are exactly the stationary distribution of the tree
proposal kernel.

Under migrate's convention (eventtime(), mcmc1.c: per-pair coalescence rate
2/theta, so rate_k = k(k-1)/theta) the total tree length is

    T = theta * sum_{j=1}^{n-1} E_j / j,   E_j iid Exp(1)

giving  E[T]/theta = H_{n-1}.  CV and skew are convention-free: they do not
depend on migrate's rate constant, so they test the sampler without requiring
any assumption about the parameterisation.

Usage:  python coalescent_check.py <path-to-migrate-n> [n_tips]
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# resolved to an absolute path: migrate is run from a temporary directory
MIGRATE = str(Path(sys.argv[1]).resolve()) if len(sys.argv) > 1 else str(Path.home() / "bin/migrate-n")
NTIPS = int(sys.argv[2]) if len(sys.argv) > 2 else 6
THETA = 0.05
SITES = 50

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
long-sample=100000
burn-in=20000
auto-tune=NO
updatefreq= tree:1.00 parameter:0.00 scaler:0.00
theta=OWN:{theta}
bayes-priors= THETA * * UNIFORMPRIOR: {lo:.7f} {hi:.7f} 0.0000005
random-seed=OWN:777
outfile=outfile
bayes-allfile=YES:1:bayesall
"""


def run(workdir):
    workdir = Path(workdir)
    with open(workdir / "infile", "w") as f:
        f.write("   1 1 . all-missing 1pop\n%i\n   %i Pop1\n" % (SITES, NTIPS))
        for i in range(NTIPS):
            f.write("%-10s%s\n" % ("ind%i" % i, "?" * SITES))
    (workdir / "parmfile").write_text(
        PARM.format(theta=THETA, lo=THETA * 0.99998, hi=THETA * 1.00002)
    )
    subprocess.run(
        [MIGRATE, "parmfile", "-nomenu"],
        cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    rows, hdr = [], None
    for line in (workdir / "bayesall").read_text().splitlines():
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if p[0] == "Steps":
            hdr = p
            continue
        rows.append([float(x) for x in p[: len(hdr)]])
    k = min(len(r) for r in rows)
    A = np.array([r[:k] for r in rows])
    d = {h: A[:, i] for i, h in enumerate(hdr[:k])}
    theta = d["Theta_1"]
    assert theta.max() - theta.min() < 1e-5, "theta was not held fixed"
    return d["treelength"] / theta.mean()


def describe(x, name):
    hi = np.quantile(x, 0.8)
    skew = ((x - x.mean()) ** 3).mean() / x.std() ** 3
    print(
        f"{name:24s} {x.mean():10.4f} {x.std()/x.mean():8.4f} {skew:8.4f} "
        f"{np.quantile(x,0.99):8.4f} {(x[x>hi]-hi).mean():12.4f}"
    )


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    exact = (rng.exponential(size=(2_000_000, NTIPS - 1)) / np.arange(1, NTIPS)).sum(1)
    with tempfile.TemporaryDirectory() as tmp:
        obs = run(tmp)
    print(f"n = {NTIPS} tips, {len(obs)} recorded genealogies\n")
    print(f"{'':24s} {'mean T/th':>10s} {'CV':>8s} {'skew':>8s} {'q99':>8s} {'tail excess':>12s}")
    describe(exact, "EXACT Kingman")
    describe(obs, "migrate")
    print(
        "\nCV, skew and tail excess are scale-free: they must match regardless of\n"
        "migrate's rate convention. Persistent disagreement means the tree update's\n"
        "stationary distribution is not the coalescent."
    )
