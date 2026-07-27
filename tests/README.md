# Sampler regression tests

## coalescent_check.py

Checks that the genealogy update actually samples the Kingman coalescent.

One population, `n` haploid tips, every site `?`, theta pinned by a degenerate
uniform prior, tree updates only. The data likelihood is then constant, so
`acceptlike()` always returns TRUE and the recorded genealogies are exactly the
stationary distribution of the tree proposal kernel — no likelihood, no
parameter moves, nothing else in the way.

    python tests/coalescent_check.py src/migrate-n [n_tips]

Compares total tree length against the exact Kingman distribution
`T = theta * sum_{j=1}^{n-1} E_j / j`, `E_j` iid Exp(1). CV, skew and tail
excess are scale-free: they must match regardless of migrate's rate convention,
so the test needs no assumption about the parameterisation.

Expected (n=6): mean T/theta 2.283, CV 0.530, skew 1.345, tail excess 1.04.

This test catches errors that are invisible with real data, where the
likelihood masks a mis-specified prior kernel. It was written after
`beyond_last_node()` was found to run two coalescence clocks for a single pair
of lineages, coalescing at 4/theta instead of 2/theta — which made genealogies
too short and too shallow (mean T/theta 1.933, CV 0.441, skew 0.903).
