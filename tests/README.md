# Sampler regression tests

Both tests run migrate on data where every site is `?`. The likelihood is then
constant, `acceptlike()` always returns TRUE, and what the chain records is
exactly the stationary distribution of the proposal machinery — no likelihood in
the way. That exposes errors which are invisible on real data, because there the
likelihood masks a mis-specified prior kernel.

## coalescent_check.py

Does the genealogy update sample the Kingman coalescent?

One population, `n` haploid tips, theta pinned by a degenerate uniform prior,
tree updates only.

    python tests/coalescent_check.py src/migrate-n [n_tips]

Compares total tree length against the exact Kingman distribution
`T = theta * sum_{j=1}^{n-1} E_j / j`, `E_j` iid Exp(1). CV, skew and tail
excess are scale-free: they must match regardless of migrate's rate convention,
so the test needs no assumption about the parameterisation.

Expected (n=6): mean T/theta 2.283, CV 0.530, skew 1.345, tail excess 1.04.

Written after `beyond_last_node()` was found to run two coalescence clocks for a
single pair of lineages, coalescing at 4/theta instead of 2/theta. Two lineages
form one pair, so only one clock may run; migration and speciation clocks stay
per-lineage. Before the fix: mean T/theta 1.933, CV 0.441, skew 0.903.

## migration_prior_check.py

Does the migration rate M come back as its prior?

Two populations, 4+4 tips, theta pinned, M free. The marginal for M must be
exactly the prior.

    python tests/migration_prior_check.py src/migrate-n [prior_scale]

Exits non-zero if the recovered mean is more than 5% off. The default prior
scale is deliberately small (~1 migration event per genealogy) — see below.

### Why the small prior scale

The M update is an independence sampler drawing from the prior, while the
conditional is

    p(M|G)  ~  M^n * exp(-L*M) * prior(M)

a truncated Gamma of shape n+1, where n is the number of migration events on the
genealogy and L the branch length in the emigrating population. Its relative
width is 1/sqrt(n). As n grows the conditional becomes a needle inside a broad
proposal, acceptance collapses, and M can then only move as fast as the
genealogy's migration count random-walks — O(n^2) tree updates per relaxation.

This is a mixing limit, not a bias. Measured on this setup:

    prior scale        10     100    1000    5000   20000
    E[n_mig]          1.0     8.7    74.7   274.7   306.4
    M acceptance    0.661   0.432   0.203   0.110   0.069
    recovered/target
                    1.001   1.002   0.982   0.839   0.461

The components were each verified exact: the tree kernel matches msprime at
M = 100/1000/4000 and asymmetric 25/400; with the tree frozen the sampled M
matches the analytic truncated Gamma (mean 12272.1 vs 12272.5, sd 751.4 vs
755.4); `probg_treetimes` has the exact form `A - L*M + n*log M` with `sum(L)`
equal to the tree length and `sum(n)` to the true migration count, correctly
paired per deme.

Switching the MIG proposal does not rescue the large-n case — SLICE gives 0.524
and MULTIPLIER 0.472 at scale 20000, against 0.461 for the default. Nor would a
Gibbs draw from the Gamma conditional: the frozen-tree test shows the M update
already reaches its conditional essentially exactly. The bottleneck is the joint
(M, genealogy) ridge, so fixing it needs a move that changes M and the migration
history together, not a better one-dimensional proposal.

Practical consequence for real runs: whenever a locus carries many migration
events, M needs long chains. The M acceptance ratio in the outfile is the
diagnostic — 0.07 means the chain is not exploring the prior range.

### Counting migration events

Timelist entries are `ntips` tips + `ntips-1` coalescences + migrations + **one
root sentinel**, so `n_mig = treeintervals - 2*ntips`. Getting this off by one
is easy, and it was the source of a spurious "probg drops one migration event"
conclusion during the investigation.
