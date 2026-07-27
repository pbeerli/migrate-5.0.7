# Sampler regression tests

These tests run migrate on data where every site is `?`. The likelihood is then
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

This is a mixing limit, not a bias — specifically **metastable trapping**, see
`overdispersed_start_check.py` below for the evidence. Measured on this setup:

    prior scale        10     100    1000    5000   20000
    E[n_mig]          1.0     8.7    74.7   274.7   306.4
    M acceptance    0.661   0.432   0.203   0.110   0.069
    recovered/target
                    1.001   1.002   0.982   0.839   0.461

Independently reproduced with 16 fresh seeds per scale: 1.0023 at scale 100 and
0.9990 at scale 1000, against 0.4767 at scale 20000 over 40 seeds.

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
(M, genealogy) ridge, so no one-dimensional proposal can help.

A move that changes M and the migration history *together* is the right shape,
but the obvious construction does not work. Proposing M' and redrawing the whole
history by stochastic mapping under Q(M') makes the CTMC path densities cancel,
leaving `alpha = [C(S,H')/C(S,H)] * [P_CTMC(M')/P_CTMC(M)]`. Prototyped and
measured: acceptance 0.98 / 0.91 / 0.053 / 0.0000 at 2 / 4 / 6 / 8 tips, and
worse with more demes. The proposal is coalescent-blind, so it destroys the
even-spread-across-demes structure that C demands; the mismatch is
theta-independent and grows with lineage count. A partial redraw trades
acceptance back against M's step size and hits the same wall. Making this work
needs a coalescent-aware proposal (BASTA/MASCOT-style), which is a research
project, not a patch.

Practical consequence for real runs: whenever a locus carries many migration
events, M needs long chains. The M acceptance ratio in the outfile is the
diagnostic — 0.07 means the chain is not exploring the prior range.

## overdispersed_start_check.py

Is the large-n deficit a slow transient, or the wrong stationary distribution?

    python tests/overdispersed_start_check.py [recorded_samples]

Starts M at 50, 5000 and 19500 (three seeds each) with **burn-in disabled**, and
prints the mean of M in each fifth of the chain. All nine runs meet inside the
first fifth and then show no trend; the grand mean is 3836 against a target of
8360 (ratio 0.459), and ten times more samples does not move it.

**Read this correctly.** Start-independence rules out a slow transient. It does
*not* rule out trapping: if every start falls into a trap basin quickly, all
starts agree, longer runs do not help, and the answer is still wrong. That
distinction was got wrong once during the investigation. The decisive evidence is
a *within-chain* diagnostic, not a between-start one:

In the fully symmetric 4+4 model M_2->1 and M_1->2 are exchangeable, so the
fraction of samples with M21 > M12 must concentrate at 0.5. Over 40 independent
seeds at prior scale 20000 that fraction ranges from **0.006 to 1.000**, with
11/40 chains locked below 0.1 or above 0.9 and only 10/40 in 0.4-0.6 — while the
mean of log(M21/M12) is +0.010 (t = 0.04, p = 0.97), so there is no systematic
direction bias, just chains stuck on one side. At prior scale 100 and 1000 the
same statistic gives 0/16 locked and a range of 0.48-0.52. Trapping, measured.

## probg_conditional_check.py

Does the M the chain holds match the conditional `probg_treetimes` scores?

Needs an instrumented build (`#ifdef DEBUGMIG` in `src/bayes.c`):

    make CPPFLAGS="-DPRETTY -DNEWVERSION -D_REENTRANT -DDEBUGMIG" migrate-n
    mv migrate-n migrate-n-debugmig && make migrate-n
    python tests/probg_conditional_check.py src/migrate-n-debugmig [scale]

The build emits, every `MIGRATE_PROBG_EVERY` calls (default 5000, 0 disables),
one line per migration parameter:

    PROBG  call  locus  heat  param  n_p  L_p  M_p  sumprob

`n_p` is the migration count attributed to parameter p and `L_p` the **pure**
migration opportunity — the coefficient of `M_p` in the waiting term, with no M
folded into it, so it is directly usable as a Gamma rate. (An earlier version of
this scaffolding folded M in and could only report `n/(L*M)`, which cannot be
compared against a sampled M.)

Two readouts:

- `n_p / (L_p * M_p)` must average 1, since the tree kernel lays down migration
  events at rate `M_p` over opportunity `L_p`. Measured **0.993** and 0.989 — the
  tree kernel and probg agree on the rate convention, ruling out the
  factor-of-two mismatch that was suspected.
- The PIT `F(M_p; n_p+1, L_p + 1/scale)` should be Uniform(0,1). Measured mean
  0.525. The reported KS p-value is not meaningful here — the samples are
  autocorrelated.

Conditional means track the sampled M (5990 vs 6041, 2281 vs 2634), and at fixed
M the genealogies match msprime **per deme** (519.8/519.8 migrations against
515.9/515.9; tree length 0.2618 against 0.2593). So `p(G|theta,M)` and `p(M|G)`
are both exact, which is why the wrong marginal has to be non-convergence.

### Counting migration events

Timelist entries are `ntips` tips + `ntips-1` coalescences + migrations + **one
root sentinel**, so `n_mig = treeintervals - 2*ntips`. Getting this off by one
is easy, and it was the source of a spurious "probg drops one migration event"
conclusion during the investigation.
