/*------------------------------------------------------------------------
  Joint multi-locus combination of Bayesian posteriors (2026-10-01; ported
  from migrate-codex-7 7.0.67 to migrate 5.0.7 6.1.42 (skyline=PARAM
  multipliers added 6.1.71, crossing-time integration and sharp splits 6.1.109): without log-scale histogram bins, the markdown report and the
  numbers-of-migrants (S/M) models, which keep the product of the
  per-locus marginals here).

  bayes_combine_loci() multiplies, parameter by parameter, the per-locus
  MARGINAL posterior histograms. Every locus has already integrated out its
  own nuisance parameters, so that product is not the posterior of
  parameters shared by all loci; on simulated two-deme data it was biased,
  worse with more loci, and its intervals undercovered in a calibration test
  with truths drawn from the prior (tests/sim).

  For every sampled genealogy G of a locus, p(G|phi) is the structured
  coalescent density of probg_treetimes_local():

    sum_i [ c_i (log 2 - log(mu Theta_i)) - A_i(g_i) / (mu Theta_i) + g_i T_i ]
    + sum_ji [ m_ji log M_ji - M_ji S_i / mu ]

  with, per population i, c_i coalescences at total age T_i,
  A_i(g) = sum k_i(k_i-1) (e^{g t1} - e^{g t0}) / g (= sum k_i(k_i-1) dt
  without growth), S_i = sum k_i dt, and m_ji migration events of type ji.
  Each locus posterior is approximated by the mixture over its sampled
  genealogies, p_l(phi|D_l) ~ mean_s p(G_s|phi) prior(phi) / Z(G_s), and the
  combined posterior

    p(phi|D) ~ prior(phi) * prod_l mean_s [ p(G_ls|phi) / Z(G_ls) ]

  is sampled with a single-site Metropolis sampler. Z(G) splits into
  one-dimensional integrals (Theta without growth, each M) and, for every
  growth group, a one-dimensional integral over g of the product of the
  group's Theta integrals; all are numerical under the run's own priors.
  The joint samples replace the combined ("All") histogram before
  calc_hpd_credibility(), so every report uses them. Handled: divergence
  (split models 'd'/'D': per genealogy and split, the waiting and
  split-event terms D(mean, std) of probg_treetimes_local() on a
  JC_NM x JC_NS grid over the prior ranges, computed only for kept
  genealogies; Z gets a two-dimensional block per split), Theta and M
  with every grouping of the connection matrix that bayes->map expresses
  ('*', '0', constant 'c', symmetric 's', mean 'm', letter groups),
  exponential growth (each growing population with its own Theta), and
  use-M=NO (xNm, the rate being xNm/Theta of the receiving population; with
  free entries and no growth: Z(G) then is, per population, an integral
  over Theta of the product of the inner xNm integrals), tied splits t/T,
  estimated locus rates, and skyline=PARAM (free '*' multipliers of Theta
  and M without growth or xNm: per segment statistics, the multipliers are
  extra coordinates of the joint sampler, and Z(G) integrates each free
  entry's multipliers by a recursion on the log value, jc_sky_logR()).
  Mittag-Leffler alpha (joint_ml.c, docs/joint_mlalpha_design.md): a
  genealogy keeps its spells instead of sums, p(G|phi) is the per-population
  spell density of mlh_logp_perpop(), and Z(G) integrates each population's
  Theta and immigration rates (reduced to their sum) at Gauss-Legendre
  nodes of its alpha; up to three populations with '*'/'0' matrices and
  uniform M priors. A growing population whose Theta has free skyline
  multipliers restarts its growth at every segment start (mlf_hastings.c):
  per segment the row keeps log A_is(g) = log sum k(k-1) int e^{g(t-t_s)} dt
  on the growth grid, T_i becomes sum (t - t_s) over its coalescences, and
  Z(G) runs the multiplier recursion at every JC_SKYG-th growth node (2026-11-17).
  Other models (ML with other forces) keep the marginal-histogram product.

  Storage follows the posterior samples (bayes-allfile): by default the
  statistics stay in memory; in the low-memory file mode (has_bayesmdimfile)
  every process instead writes "<outfile>.joint.<locus>.<replicate>", and
  at the end the process that evaluates a locus reads and thins them (needs
  a shared file system, like the posterior-sample files). In memory, each
  process keeps a thinned buffer per locus (world->jointstats, at most
  2*jc_maxsamples() rows: when it fills, every other row is dropped and the
  recording stride doubles). Under MPI, replicate workers send their
  buffers to the locus worker with the other replicate results
  (mpi_send_replicate(), jc_pack_buffer()/jc_unpack_buffer()). The rows
  never go to the master: the locus workers evaluate their own loci for
  every proposal of the master's sampler (MIGMPI_JC, jc_worker_service(),
  see jc_eval_start()); a serial run uses its own buffers directly.
------------------------------------------------------------------------*/
#include "migration.h"
#include "world.h"
#include "bayes.h"
#include "tools.h"
#include "random.h"
#include "sighandler.h"
#include "joint_combine.h"
#include "speciate.h"
#include "joint_ml.h"
#include "mlf_hastings.h"
#ifdef MPI
#include "migrate_mpi.h"
extern const MPI_Datatype mpisizeof;
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define JC_MAXSAMPLES 4000   /* per locus used in the combination for two populations
                                 (at 200 loci, 1000 gave 95% coverage of 0.4 for an M,
                                 4000 0.8); scaled by (numpop-1)^2, see jc_maxsamples() */
#define JC_GRID 600          /* points for a one-dimensional normalizer */
#define JC_TGRID 200         /* Theta points inside a growth integral */
#define JC_NG 129            /* growth grid */
#define JC_SKYG 4            /* growth with skyline multipliers: Z(G) at every 4th node */
#define JC_NM 33             /* divergence-time mean grid, local to each genealogy (quadratic) */
#define JC_DZ 8.0            /* local window: centre +- JC_DZ std */
#define JC_DOUT (-1e30)      /* divergence term outside a row's window */
#define JC_NS 17             /* divergence-time std grid (log-spaced cells) */
#define JC_NGATE 257         /* gated migration: lineage time before the split mean, on a
                                grid over the prior range of that mean */
#define JC_BURNIN 2000       /* sweeps */
#define JC_BURNWARM 500      /* sweeps of a chain started from the main run's last state
                                (error blocks, stepping stones) */
#define JC_SWEEPS 50000      /* main run: at most this many sweeps ... */
#define JC_PROPBUDGET 450000 /* ... within this many proposals (50000 sweeps of 9 parameters) */
#define JC_MINSWEEPS 5000    /* ... and at least this many sweeps */
#define JC_MAXACTIVE 50      /* more free parameters: no joint combination (the reweighting
                                over fixed genealogies collapses; restrict the model) */
#define JC_INCR_UNITS 32     /* more terms per genealogy: update its sum by differences */
#define JC_REFRESH 1000      /* ... and recompute the sums after this many accepted steps */
#define JC_SSRUNGS 16        /* stepping stones for the joint marginal likelihood (24 x 1500
                                halved the error of log C, 0.06-0.14 instead of 0.11-0.26,
                                at about 1.4x the time; 2026-10-06) */
#define JC_SSBURN 300        /* sweeps per rung: burn-in, samples */
#define JC_SSN 600
#define JC_NBLOCK 4          /* blocks of every locus' genealogies for the Monte Carlo error */
#define JC_MCSWEEPS 3000     /* sweeps per block (12500 gave the same errors: the block
                                differences come from the genealogies, not the chain) */

/* the run options the joint combination needs from option_fmt (5.0.7 has
   no canonical options): set by jc_set_run_options() from main() */
static char jc_outfile[LINESIZE] = "outfile";
static boolean jc_recover = FALSE;

/// remembers the outfile name and recover=YES for the joint combination
void
jc_set_run_options (const char *outfilename, boolean recover)
{
  if (outfilename != NULL && outfilename[0] != '\0')
    snprintf (jc_outfile, LINESIZE, "%s", outfilename);
  jc_recover = recover;
}

/* M or xNm: migration_is_M() in world.c */

/* ------------------------------------------------------------- layout */

typedef struct
{
  long numpop, numpop2, np;
  long nrow;                /* doubles per genealogy */
  long off_c, off_K, off_S, off_m;
  long off_mu;              /* the locus rate of the genealogy (estimated rate: per sample) */
  boolean rate;             /* mutation=ESTIMATE: the shared rate is a joint coordinate */
  long rate_p;              /* ... its parameter index */
  long off_lr, nrk;         /* ... log p(D|G,r) on nrk log-spaced rates, then log Lambda(G) */
  double rlo, rhi;          /* ... log of the first and last grid rate */
  long off_rep;             /* the replicate that sampled it (last entry of the row) */
  long nslot;               /* growth groups */
  long *pop_slot;           /* per population: growth group or -1 */
  long *slot_param;         /* per group: its parameter index */
  long *off_t;              /* per population: offset of T_i (growing only) */
  long *off_A;              /* per population: offset of A_i on the grid */
  long *off_sA;             /* per population: offset of log A_is per segment on the grid
                               (growth with free skyline multipliers), or -1 */
  double *ggrid;            /* per group: JC_NG growth values */
  long *rep;                /* per matrix entry: its parameter, or -1 (fixed) */
  long *mslot;              /* per matrix entry: its migration-count slot in a row, or -1
                               (only entries with a parameter: a 50-population stepping
                               stone keeps 100 of 2450) */
  long nmig;                /* migration-count slots per row (and per segment) */
  double *fixed;            /* per matrix entry: value when fixed */
  boolean usem;             /* migration parameters are M (else xNm) */
  boolean *xnm;             /* per matrix entry: parameter is xNm (use-M=NO, S/M) */
  long nsplit;              /* divergence models */
  long *split_model;        /* per split: index into world->species_model */
  long *split_pmu, *split_psig; /* per split: parameter indices (psig -1: fixed) */
  long *split_ns;           /* per split: std grid points (1 when fixed) */
  long *split_rep;          /* per split: first split of its tied group (t/T), else itself */
  long *split_par;          /* per group (its first split): the group of the split that
                               created its ancestor population, -1 if none or unordered
                               (divergence-order TREE/RANK: the mean split time is at most
                               that of the parent) */
  double *split_sfix;       /* per split: the fixed std */
  long *off_D;              /* per split: offset of the D(mu, sigma) grid */
  double *mgrid, *sgrid;    /* per split: JC_NM means, JC_NS stds */
  double *swidth;           /* per split: width of every std cell */
  double *mlo, *mhi;        /* per split: prior range of the mean */
  /* skyline=PARAM: segment s = [segt[s], segt[s+1]); the free ('*') matrix
     entries have multipliers tau_s for s >= 1 (tau_0 = 1), sampled as extra
     coordinates np .. nphi-1 of phi */
  long nseg;                /* 1 without segments */
  double *segt;             /* nseg + 1 boundaries */
  long off_sc, off_sK, off_sS, off_sm; /* per segment: c, K, S per population, m per entry */
  long nfree, ntau, nphi;
  long *sky_j;              /* per matrix entry: free-entry index or -1 */
  long *sky_k;              /* per free entry: its matrix entry */
  int *sky_type;            /* per free entry: SKYPRIOR_RANDOMWALK or SKYPRIOR_LOGUNIFORM */
  double *sky_a, *sky_b;    /* per free entry: sigma, or min max of tau */
  double *sky_lo, *sky_hi;  /* per free entry: the prior range of its value (parameter x tau) */
  /* Mittag-Leffler alpha: a genealogy's spells travel beside its row */
  boolean ml;
  long off_nsp;             /* per population: its number of spells */
  long off_inh;             /* the inheritance scalar of the locus */
  long off_lc;              /* ML: log sites of the locus (mittag-leffler-scale=LOCUS) */
  int mlscale;              /* ML: mittag-leffler-scale */
  long *apar;               /* per population: its alpha parameter, or -1 (fixed) */
  double *afix;             /* per population: its fixed alpha (1 without ML) */
  /* gated migration (forward in time out of a daughter population: backward,
     a lineage may move into the daughter only up to the mean split time D,
     species_mig_gate()): the waiting term of the entry is M S(D) / mu with
     S(D) the lineage time of the receiving population before D, and D must
     be older than the entry's last migration event */
  long ngate;
  long *gate_k;            /* per matrix entry: the split whose mean gates it, or -1 */
  long *off_G;             /* per matrix entry: offset of [last event age, S on the grid] */
  double *gate_b, *gate_c; /* per matrix entry: upper bound and log density of its
                              uniform M prior (Z(G) integrates M in closed form) */
  /* integrated crossing time (MIGRATE_JC_RB, docs/joint_split_rb_design.md):
     one crossing lineage per genealogy is integrated over its segment */
  boolean rb;
  long off_R;               /* the block: split (-1 none), a, b, pieces, then per
                               piece its start, 2 k_d / inh, 2 k_a / inh */
  long *split_to, *split_from; /* per split: daughter and ancestor population */
  double *split_b1;         /* per split: upper bound of the split time (s->max) */
  boolean *rb_ok;           /* per split: its crossings can be integrated */
  boolean *sharp_pu;        /* sharp splits: a daughter or ancestor with its own
                               population unit, evaluated on the moved row
                               (set by jc_active) */
  boolean sharp;            /* sharp splits (S = 0): all crossings of a split move
                               together, integrated over their common window */
  long *sharp_par;          /* sharp: per split the split that created its ancestor, -1 */
  long rbw;                 /* width of the off_R block */
  boolean *mzero;           /* per matrix entry: no migration term ('0', 'd', 't') */
} jc_layout;

#define JC_RBMAXP 24         /* pieces of a crossing segment at most */
#define JC_RBMAXSKY 32       /* skyline=PARAM: free entries at most */
#define JC_RBW (5 + 3 * JC_RBMAXP)
#define JC_RBN(R) ((R)[4 + 3 * JC_RBMAXP])   /* crossings that move (1 without sharp) */
#define JC_SHARPMAXN 128     /* sharp split: crossing lineages at most */
#define JC_SHARP_REL 1e-5    /* sharp split: the S prior maximum below this x the D maximum */
#define JC_SHARPMAXI 96      /* sharp splits: intervals of a row's count timeline at most */
/* sharp row block: status (0 ok, -1 no crossing, -2 unusable), intervals,
   start, inheritance scalar; per split A, B, n, c0, chain; then per interval
   its end and the lineage count of every population */
#define JC_SH_SPLIT(k) (4 + 5 * (k))
#define JC_SH_TL(ly) (4 + 5 * (ly)->nsplit)

/* phi index of tau_s (s >= 1) of free entry j */
#define JC_TAU(ly, j, s) ((ly)->np + (j) * ((ly)->nseg - 1) + (s) - 1)

static void
jc_layout_free (jc_layout *ly)
{
  myfree (ly->pop_slot);
  myfree (ly->slot_param);
  myfree (ly->off_t);
  myfree (ly->off_A);
  myfree (ly->off_sA);
  myfree (ly->ggrid);
  myfree (ly->rep);
  myfree (ly->mslot);
  myfree (ly->xnm);
  myfree (ly->fixed);
  myfree (ly->split_model);
  myfree (ly->split_pmu);
  myfree (ly->split_psig);
  myfree (ly->split_ns);
  myfree (ly->split_rep);
  myfree (ly->split_par);
  myfree (ly->split_sfix);
  myfree (ly->off_D);
  myfree (ly->mgrid);
  myfree (ly->sgrid);
  myfree (ly->swidth);
  myfree (ly->mlo);
  myfree (ly->mhi);
  myfree (ly->segt);
  myfree (ly->sky_j);
  myfree (ly->sky_k);
  myfree (ly->sky_type);
  myfree (ly->sky_a);
  myfree (ly->sky_b);
  myfree (ly->sky_lo);
  myfree (ly->sky_hi);
  myfree (ly->apar);
  myfree (ly->afix);
  myfree (ly->gate_k);
  myfree (ly->off_G);
  myfree (ly->gate_b);
  myfree (ly->gate_c);
  myfree (ly->split_to);
  myfree (ly->split_from);
  myfree (ly->split_b1);
  myfree (ly->rb_ok);
  myfree (ly->sharp_pu);
  myfree (ly->sharp_par);
  myfree (ly->mzero);
}

static boolean jc_nodiv = FALSE;   /* jc_log_normalizer() without the divergence blocks (sharp splits) */
static boolean jc_rb_normal = TRUE;   /* integrated crossing: normal split times (else exponential) */

/* experimental: integrate the crossing time (MIGRATE_JC_RB=1) */
static boolean
jc_rb_mode (void)
{
  const char *e = getenv ("MIGRATE_JC_RB");
  return e != NULL && atol (e) != 0;
}

/* mutation=ESTIMATE estimates one rate for all loci. The joint treats it as
   a coordinate when the per-locus prior p(G|params) does not depend on it
   in reference time: Theta and M only (no divergence, growth or skyline
   segments, whose times the sampler compares with the rate-scaled node
   ages). Otherwise the rate stays a per-locus nuisance parameter. */
#define JC_RATEK 24          /* rate grid points (MIGRATE_JC_RATEK) */
#define JC_RATESUB 8         /* Lambda: integration points per grid interval */
static boolean
jc_rate_shared (world_fmt *world)
{
  const long p = world->numpop2;
  if (world->bayes == NULL || !world->bayes->mu || world->has_mlalpha || world->has_speciation
      || world->has_growth || world->timeelements > 2 || getenv ("MIGRATE_JC_NORATE") != NULL)
    return FALSE;
  return world->bayes->maxparam[p] > world->bayes->minparam[p] && world->bayes->maxparam[p] > 0.0;
}

static void
jc_layout_make (world_fmt *world, jc_layout *ly)
{
  const long numpop = world->numpop;
  const long npx = world->numparamcumvec[SPLITSTDPRIOR];
  const long npg = world->numparamcumvec[GROWTHPRIOR];
  long pop, s, g, p0, p;
  memset (ly, 0, sizeof (*ly));
  ly->numpop = numpop;
  ly->numpop2 = world->numpop2;
  ly->np = world->numparam;
  ly->off_c = 0;
  ly->off_K = numpop;
  ly->off_S = 2 * numpop;
  ly->off_m = 3 * numpop;
  ly->usem = world->options->usem;
  ly->rep = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->xnm = (boolean *) mycalloc ((size_t) world->numpop2, sizeof (boolean));
  for (p = numpop; p < world->numpop2; p++)
    ly->xnm[p] = !migration_is_M (world, p);
  ly->fixed = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  for (p = 0; p < world->numpop2; p++)
    {
      const long r = world->bayes->map[p][1];
      ly->rep[p] = (r >= 0 && r < world->numpop2 && strchr ("0c", world->options->custm2[p]) == NULL) ? r : -1;
      ly->fixed[p] = world->param0[p];
    }
  ly->mslot = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->nmig = 0;
  for (p = 0; p < world->numpop2; p++)
    ly->mslot[p] = (p >= numpop && ly->rep[p] >= 0) ? ly->nmig++ : -1;
  ly->nrow = 3 * numpop + ly->nmig;
  ly->pop_slot = (long *) mycalloc ((size_t) numpop, sizeof (long));
  ly->off_t = (long *) mycalloc ((size_t) numpop, sizeof (long));
  ly->off_A = (long *) mycalloc ((size_t) numpop, sizeof (long));
  for (pop = 0; pop < numpop; pop++)
    {
      ly->pop_slot[pop] = -1;
      if (world->has_growth && world->options->growpops[pop] != 0)
        {
          ly->pop_slot[pop] = world->options->growpops[pop] - 1;
          if (ly->pop_slot[pop] + 1 > ly->nslot)
            ly->nslot = ly->pop_slot[pop] + 1;
        }
    }
  ly->slot_param = (long *) mycalloc ((size_t) (ly->nslot + 1), sizeof (long));
  ly->ggrid = (double *) mycalloc ((size_t) ((ly->nslot + 1) * JC_NG), sizeof (double));
  for (s = 0; s < ly->nslot; s++)
    ly->slot_param[s] = -1;
  for (p0 = npx; world->has_growth && p0 < npg; p0++)
    {
      if (shortcut (p0, world, &p) || p < npx || p >= npg)
        continue;
      s = world->options->growpops[p - npx] - 1;
      if (s >= 0 && s < ly->nslot && ly->slot_param[s] < 0)
        ly->slot_param[s] = p;
    }
  for (s = 0; s < ly->nslot; s++)
    {
      const long sp = ly->slot_param[s];
      const double lo = sp >= 0 ? world->bayes->minparam[sp] : -1.0;
      const double hi = sp >= 0 ? world->bayes->maxparam[sp] : 1.0;
      for (g = 0; g < JC_NG; g++)
        ly->ggrid[s * JC_NG + g] = g == JC_NG - 1 ? hi : lo + (hi - lo) * (double) g / (JC_NG - 1);
    }
  for (pop = 0; pop < numpop; pop++)
    if (ly->pop_slot[pop] >= 0)
      {
        ly->off_t[pop] = ly->nrow;
        ly->off_A[pop] = ly->nrow + 1;
        ly->nrow += 1 + JC_NG;
      }
  /* divergence: one D(mu, sigma) grid per split model */
  const long nsm = world->has_speciation ? world->species_model_size : 0;
  ly->split_model = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_pmu = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_psig = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_ns = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_rep = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_par = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->split_sfix = (double *) mycalloc ((size_t) (nsm + 1), sizeof (double));
  ly->off_D = (long *) mycalloc ((size_t) (nsm + 1), sizeof (long));
  ly->mgrid = (double *) mycalloc ((size_t) ((nsm + 1) * JC_NM), sizeof (double));
  ly->sgrid = (double *) mycalloc ((size_t) ((nsm + 1) * JC_NS), sizeof (double));
  ly->swidth = (double *) mycalloc ((size_t) ((nsm + 1) * JC_NS), sizeof (double));
  ly->mlo = (double *) mycalloc ((size_t) (nsm + 1), sizeof (double));
  ly->mhi = (double *) mycalloc ((size_t) (nsm + 1), sizeof (double));
  for (s = 0; s < nsm; s++)
    {
      species_fmt *sp = &world->species_model[s];
      const long pm = sp->paramindex_mu, ps = sp->paramindex_sigma;
      const long k = ly->nsplit;
      ly->split_model[k] = s;
      ly->split_pmu[k] = pm;
      ly->split_psig[k] = (ps >= 0 && ps < ly->np && world->bayes->map[ps][1] == ps
                           && world->bayes->maxparam[ps] > world->bayes->minparam[ps]) ? ps : -1;
      ly->split_ns[k] = ly->split_psig[k] >= 0 ? JC_NS : 1;
      /* tied splits (t/T): every split of the same ancestor shares the D and
         S of the first one (the sampler copies them, speciate.c) */
      ly->split_rep[k] = k;
      if (sp->type == 't')
        {
          long j;
          for (j = 0; j < k; j++)
            if (world->species_model[ly->split_model[j]].type == 't'
                && world->species_model[ly->split_model[j]].from == sp->from)
              {
                ly->split_rep[k] = ly->split_rep[j];
                ly->split_pmu[k] = ly->split_pmu[j];
                ly->split_psig[k] = ly->split_psig[j];
                ly->split_ns[k] = ly->split_ns[j];
                break;
              }
        }
      ly->split_sfix[k] = (ps >= 0 && ps < ly->np) ? world->param0[ps] : sp->sigma;
      ly->mlo[k] = world->bayes->minparam[ly->split_pmu[k]];
      ly->mhi[k] = world->bayes->maxparam[ly->split_pmu[k]];
      /* std: log-spaced rows from 1e-3 of the upper bound (or the prior's
         own lower bound) up to the upper bound itself, so that no std of the
         prior range lies above the last row (a constant beyond it tilted
         every genealogy's weight towards large std); for a small std the
         divergence term is a narrow peak in the mean and needs this
         resolution */
      if (ly->split_psig[k] < 0)
        {
          ly->sgrid[k * JC_NS] = ly->split_sfix[k];
          ly->swidth[k * JC_NS] = 1.0;
        }
      else
        {
          const long psr = ly->split_psig[k];
          const double shi = world->bayes->maxparam[psr];
          const double slo = world->bayes->minparam[psr] > 1e-3 * shi ? world->bayes->minparam[psr] : 1e-3 * shi;
          for (g = 0; g < JC_NS; g++)
            {
              ly->sgrid[k * JC_NS + g] = slo * pow (shi / slo, (double) g / (JC_NS - 1));
              ly->swidth[k * JC_NS + g] = ly->sgrid[k * JC_NS + g] * log (shi / slo) / (JC_NS - 1);
            }
        }
      ly->off_D[k] = ly->nrow;           /* centre, split-event count, grid */
      ly->nrow += 2 + JC_NM * ly->split_ns[k];
      ly->nsplit++;
    }
  /* divergence-order TREE or RANK: the parent group of a group is the group
     of the split that must be older (species_order_partner: TREE the split
     that created the ancestor population, RANK the split listed before);
     each group has at most one parent, so the ordered groups form trees */
  for (s = 0; s < ly->nsplit; s++)
    {
      ly->split_par[s] = -1;
      if (ly->split_rep[s] != s)
        continue;
      long m, q;
      for (m = 0; m < ly->nsplit && ly->split_par[s] < 0; m++)
        {
          if (ly->split_rep[m] != s)
            continue;
          const long partner = species_order_partner (world, ly->split_model[m]);
          if (partner < 0)
            continue;
          for (q = 0; q < ly->nsplit; q++)
            if (ly->split_model[q] == partner && ly->split_rep[q] != s)
              {
                ly->split_par[s] = ly->split_rep[q];
                break;
              }
        }
    }
  /* gated migration entries: the gate is the first split whose daughter is
     the entry's source population (as species_mig_gate()) */
  ly->gate_k = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->off_G = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->gate_b = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  ly->gate_c = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  for (p = 0; p < world->numpop2; p++)
    {
      long from, to, k;
      ly->gate_k[p] = -1;
      if (p < numpop || ly->rep[p] < 0 || strchr ("0dDtT", world->options->custm2[p]) != NULL)
        continue;
      m2mm (p, numpop, &from, &to);
      for (k = 0; k < ly->nsplit; k++)
        if (world->species_model[ly->split_model[k]].to == from)
          {
            const long r = ly->rep[p];
            ly->gate_k[p] = k;
            ly->off_G[p] = ly->nrow;
            ly->nrow += 1 + JC_NGATE;
            ly->gate_b[p] = world->bayes->maxparam[r];
            ly->gate_c[p] = scaling_prior_joint (world, r, 0.5 * world->bayes->maxparam[r]);
            ly->ngate++;
            break;
          }
    }
  /* skyline=PARAM segments */
  ly->nseg = world->timeelements > 2 ? world->timeelements - 1 : 1;
  ly->sky_j = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->sky_k = (long *) mycalloc ((size_t) world->numpop2, sizeof (long));
  ly->sky_type = (int *) mycalloc ((size_t) world->numpop2, sizeof (int));
  ly->sky_a = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  ly->sky_b = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  ly->sky_lo = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  ly->sky_hi = (double *) mycalloc ((size_t) world->numpop2, sizeof (double));
  ly->segt = (double *) mycalloc ((size_t) (ly->nseg + 1), sizeof (double));
  for (p = 0; p < world->numpop2; p++)
    ly->sky_j[p] = -1;
  ly->segt[ly->nseg] = HUGE_VAL;
  if (ly->nseg > 1)
    {
      for (s = 0; s < ly->nseg; s++)
        ly->segt[s] = world->times[s];
      for (p = 0; p < world->numpop2; p++)
        if (world->skymodel2 != NULL && world->skymodel2[p] == '*' && ly->rep[p] >= 0)
          {
            const long j = ly->nfree++;
            ly->sky_j[p] = j;
            ly->sky_k[j] = p;
            ly->sky_type[j] = world->skyptype[p];
            ly->sky_a[j] = world->skypa[p];
            ly->sky_b[j] = world->skypb[p];
            ly->sky_lo[j] = world->bayes->minparam[ly->rep[p]];
            ly->sky_hi[j] = world->bayes->maxparam[ly->rep[p]];
          }
      ly->off_sc = ly->nrow;
      ly->nrow += ly->nseg * numpop;
      ly->off_sK = ly->nrow;
      ly->nrow += ly->nseg * numpop;
      ly->off_sS = ly->nrow;
      ly->nrow += ly->nseg * numpop;
      ly->off_sm = ly->nrow;
      ly->nrow += ly->nseg * ly->nmig;
    }
  ly->off_sA = (long *) mycalloc ((size_t) numpop, sizeof (long));
  for (pop = 0; pop < numpop; pop++)
    {
      ly->off_sA[pop] = -1;
      if (ly->pop_slot[pop] >= 0 && ly->sky_j[pop] >= 0)
        {
          ly->off_sA[pop] = ly->nrow;
          ly->nrow += ly->nseg * JC_NG;
        }
    }
  ly->ntau = ly->nfree * (ly->nseg - 1);
  ly->nphi = ly->np + ly->ntau;
  /* Mittag-Leffler alpha: spell counts and the inheritance scalar */
  ly->apar = (long *) mycalloc ((size_t) numpop, sizeof (long));
  ly->afix = (double *) mycalloc ((size_t) numpop, sizeof (double));
  for (pop = 0; pop < numpop; pop++)
    {
      ly->apar[pop] = -1;
      ly->afix[pop] = 1.0;
    }
  if (world->has_mlalpha)
    {
      ly->ml = TRUE;
      for (pop = 0; pop < numpop && pop < world->options->mlalphapops_numalloc; pop++)
        {
          const long gw = world->options->mlalphapops[pop];
          if (gw <= 0)
            continue;
          ly->afix[pop] = world->mlalpha[gw - 1];
          if (world->tri_mlalpha != FIXED)
            {
              /* populations sharing an alpha ({1 1}) use the parameter of
                 the first of them (pick_parameter()); bayes->map is the
                 identity here and pointed the second one at an unused
                 slot without a prior */
              long r = -1;
              ly->apar[pop] = (!pick_parameter (npg + pop, &r, world, GROWTHPRIOR)
                               && r >= npg && r < ly->np) ? r : -1;
            }
        }
      ly->off_nsp = ly->nrow;
      ly->nrow += numpop;
      ly->off_inh = ly->nrow;
      ly->nrow += 1;
      ly->off_lc = ly->nrow;
      ly->nrow += 1;
      ly->mlscale = world->options->mlscale;
    }
  /* integrated crossing time: normal or exponential split times, no growth,
     skyline segments or Mittag-Leffler; an ancestor that is itself a
     daughter is excluded (its own split terms would depend on the crossing) */
  ly->split_to = (long *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (long));
  ly->split_from = (long *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (long));
  ly->split_b1 = (double *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (double));
  ly->rb_ok = (boolean *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (boolean));
  ly->sharp_pu = (boolean *) mycalloc ((size_t) (ly->numpop + 1), sizeof (boolean));
  ly->mzero = (boolean *) mycalloc ((size_t) world->numpop2, sizeof (boolean));
  for (p = 0; p < world->numpop2; p++)
    ly->mzero[p] = p < numpop || strchr ("0dt", world->options->custm2[p]) != NULL;
  for (s = 0; s < ly->nsplit; s++)
    {
      const species_fmt *sp = &world->species_model[ly->split_model[s]];
      ly->split_to[s] = sp->to;
      ly->split_from[s] = sp->from;
      ly->split_b1[s] = sp->max;
    }
  if (jc_rb_mode () && ly->nsplit > 0 && !ly->ml && ly->nslot == 0 && ly->nfree <= JC_RBMAXSKY
      && (world->species_model_dist == NORMAL_DIST || world->species_model_dist == EXP_DIST))
    {
      long q;
      for (s = 0; s < ly->nsplit; s++)
        {
          ly->rb_ok[s] = ly->mhi[s] > ly->mlo[s];   /* not for a fixed D */
          for (q = 0; q < ly->nsplit; q++)
            if (ly->split_to[q] == ly->split_from[s])
              ly->rb_ok[s] = FALSE;
          if (ly->rb_ok[s])
            ly->rb = TRUE;
        }
    }
  /* sharp splits (docs/joint_split_rb_design.md, "Sharp splits"): S fixed at
     or bounded near 0 (splitstd=0), every lineage crosses at D. A
     genealogy's split terms are then point masses; all crossings of a split
     slide together within their common window without changing the
     topology, so the joint integrates one time per split (exact). Every
     split must be sharp (normal split model); tied splits share one time,
     growth moves the growth-weighted integrals, and gated migration into a
     daughter of the same split time is an ordinary entry (below). Not with
     skyline segments, Mittag-Leffler or other gated entries.
     MIGRATE_JC_SHARP=0 turns it off; MIGRATE_JC_SHARPCHECK=1 checks the
     terms against the sampler's p(G|params). */
  ly->sharp_par = (long *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (long));
  for (s = 0; s < ly->nsplit; s++)
    {
      long q;
      ly->sharp_par[s] = -1;
      for (q = 0; q < ly->nsplit; q++)
        if (q != s && ly->split_to[q] == ly->split_from[s])
          ly->sharp_par[s] = q;
    }
  /* normal split times only: the exponential split model has no spread (S is
     fixed at 0 there) and its crossings are never all at D */
  /* gated migration (out of a daughter, forward in time, only after its
     split): with every crossing at D a daughter holds no lineage above D, so
     an entry whose receiving population is a daughter of the same split time
     (the gating split's tied group) is exposed over all of that
     population's lineage time: an ordinary migration entry. Other gated
     entries keep the sharp terms off */
  boolean gates_ok = TRUE;
  for (p = 0; p < world->numpop2 && ly->ngate > 0; p++)
    if (ly->gate_k[p] >= 0)
      {
        long from, to, q;
        boolean same = FALSE;
        m2mm (p, numpop, &from, &to);
        for (q = 0; q < ly->nsplit; q++)
          if (ly->split_to[q] == to && ly->split_rep[q] == ly->split_rep[ly->gate_k[p]])
            same = TRUE;
        if (!same)
          gates_ok = FALSE;
      }
  if (!ly->rb && ly->nsplit >= 1 && !ly->ml && ly->nseg == 1 && gates_ok
      && world->species_model_dist == NORMAL_DIST
      && !(getenv ("MIGRATE_JC_SHARP") != NULL && atol (getenv ("MIGRATE_JC_SHARP")) == 0))
    {
      boolean all = TRUE;
      for (s = 0; s < ly->nsplit; s++)
        {
          const double smax = ly->split_psig[s] >= 0 ? world->bayes->maxparam[ly->split_psig[s]] : ly->split_sfix[s];
          if (!(smax <= JC_SHARP_REL * ly->mhi[s]) || !(ly->mhi[s] > ly->mlo[s]))
            all = FALSE;
        }
      if (all)
        {
          ly->sharp = TRUE;
          ly->rb = TRUE;
          for (s = 0; s < ly->nsplit; s++)
            ly->rb_ok[s] = TRUE;
          for (p = 0; p < world->numpop2; p++)   /* the gates become ordinary entries */
            if (ly->gate_k[p] >= 0)
              {
                ly->gate_k[p] = -1;
                ly->ngate--;
              }
        }
    }
  if (ly->rb)
    {
      jc_rb_normal = world->species_model_dist == NORMAL_DIST;
      ly->off_R = ly->nrow;
      ly->rbw = ly->sharp ? JC_SH_TL (ly) + JC_SHARPMAXI * (1 + numpop) : JC_RBW;
      ly->nrow += ly->rbw;
    }
  ly->off_mu = ly->nrow;
  ly->nrow += 1;
  if (jc_rate_shared (world))
    {   /* the shared rate (docs/joint_shared_rate_design.md): the data
           likelihood of the genealogy on a log grid over the rate prior */
      const long p = world->numpop2;
      const double hi = world->bayes->maxparam[p];
      const double lo = world->bayes->minparam[p] > 0.0 ? world->bayes->minparam[p] : 1e-3 * hi;
      ly->rate = TRUE;
      ly->rate_p = p;
      ly->nrk = getenv ("MIGRATE_JC_RATEK") ? atol (getenv ("MIGRATE_JC_RATEK")) : JC_RATEK;
      if (ly->nrk < 4)
        ly->nrk = 4;
      ly->rlo = log (lo);
      ly->rhi = log (hi);
      ly->off_lr = ly->nrow;
      ly->nrow += ly->nrk + 1;
    }
  ly->off_rep = ly->nrow;   /* the replicate, last: the Monte Carlo error blocks */
  ly->nrow += 1;
}

/* the segment of age t: the last s with segt[s] <= t (mlh_seg()) */
static long
jc_seg (const jc_layout *ly, double t)
{
  long s = ly->nseg - 1;
  while (s > 0 && t < ly->segt[s])
    s--;
  return s;
}

static double jc_ml_term (const jc_layout *ly, const double *st, const double *sp, long pop, const double *phi);

/* ------------------------------------------------------------- support */

static long jc_nactive (world_fmt *world);
static double jc_div_term (const jc_layout *ly, const double *st, long k, double mu, double sigma);
static double jc_full_term (const jc_layout *ly, const boolean *active, const double *st, const double *sp,
                            const double *phi, double mu);

/* gated migration entry e: the age of grid point g of its split's mean */
static double
jc_gate_x (const jc_layout *ly, long e, long g)
{
  const long k = ly->gate_k[e];
  return ly->mlo[k] + (ly->mhi[k] - ly->mlo[k]) * (double) g / (JC_NGATE - 1);
}

/* S_e(D): lineage time of the receiving population before D (linear on the
   grid; S is piecewise linear in D) */
static double
jc_gate_S (const jc_layout *ly, const double *st, long e, double D)
{
  const long k = ly->gate_k[e];
  const double *S = st + ly->off_G[e] + 1;
  const double lo = ly->mlo[k], hi = ly->mhi[k];
  double x = (hi > lo) ? (D - lo) / (hi - lo) * (JC_NGATE - 1) : 0.0;
  long j;
  if (x <= 0.0)
    return S[0];
  if (x >= JC_NGATE - 1)
    return S[JC_NGATE - 1];
  j = (long) floor (x);
  x -= (double) j;
  return S[j] + x * (S[j + 1] - S[j]);
}

/* migration term of a gated entry e at rate x and split mean D: zero density
   when D is younger than its last migration event (bayes.c: a move into a
   daughter before it exists, forward in time) */
static double
jc_gate_term (const jc_layout *ly, const double *st, long e, double x, double D, double mu)
{
  if (x <= 0.0 || D < st[ly->off_G[e]])
    return -HUGE_VAL;
  const double m = ly->mslot[e] >= 0 ? st[ly->off_m + ly->mslot[e]] : 0.0;
  return (m > 0.0 ? m * log (x) : 0.0) - x * jc_gate_S (ly, st, e, D) / mu;
}

/* log int prior(M) M^m e^{-M s} dM for the uniform prior on [0, b]:
   e^c gamma(m + 1, b s) / s^(m + 1) */
static double
jc_gate_logI (const jc_layout *ly, long e, double m, double s)
{
  const double b = ly->gate_b[e], c = ly->gate_c[e];
  if (!(s * b > 1e-12))
    return c + (m + 1.0) * log (b) - log (m + 1.0);
  return c + lgamma (m + 1.0) + logincompletegamma (b * s, m + 1.0) - (m + 1.0) * log (s);
}

/* the gated M integrals of one genealogy on the grid of their split's mean
   (set by jc_div_normalizer() for the genealogy it integrates) */
static double *jc_gate_cur = NULL;

static void
jc_gate_table (const jc_layout *ly, const double *st, double mu)
{
  long e, g, n = 0;
  jc_gate_cur = (double *) myrealloc (jc_gate_cur, sizeof (double) * (size_t) ((ly->ngate > 0 ? ly->ngate : 1) * JC_NGATE));
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->gate_k[e] >= 0)
      {
        const double m = ly->mslot[e] >= 0 ? st[ly->off_m + ly->mslot[e]] : 0.0;
        for (g = 0; g < JC_NGATE; g++)
          jc_gate_cur[n * JC_NGATE + g] = jc_gate_logI (ly, e, m, st[ly->off_G[e] + 1 + g] / mu);
        n++;
      }
}

/* sum of the gated M integrals of the entries gated by split group k at
   mean D (-HUGE_VAL below an entry's last event) */
static double
jc_gate_group (const jc_layout *ly, const double *st, long k, double D)
{
  long e, n = 0;
  double v = 0.0;
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->gate_k[e] >= 0)
      {
        if (ly->split_rep[ly->gate_k[e]] == k)
          {
            const double *I = jc_gate_cur + n * JC_NGATE;
            const long kk = ly->gate_k[e];
            if (D < st[ly->off_G[e]])
              return -HUGE_VAL;
            double x = (ly->mhi[kk] > ly->mlo[kk]) ? (D - ly->mlo[kk]) / (ly->mhi[kk] - ly->mlo[kk]) * (JC_NGATE - 1) : 0.0;
            long j;
            if (x <= 0.0)
              v += I[0];
            else if (x >= JC_NGATE - 1)
              v += I[JC_NGATE - 1];
            else
              {
                j = (long) floor (x);
                x -= (double) j;
                v += I[j] + x * (I[j + 1] - I[j]);
              }
          }
        n++;
      }
  return v;
}

/* TRUE when matrix entry p is the parameter of a gated entry */
static boolean
jc_gate_param (const jc_layout *ly, long p)
{
  long e;
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->gate_k[e] >= 0 && ly->rep[e] == p)
      return TRUE;
  return FALSE;
}

/* experiments: MIGRATE_JC_MAXACTIVE raises the parameter limit */
static long
jc_maxactive (void)
{
  return getenv ("MIGRATE_JC_MAXACTIVE") ? atol (getenv ("MIGRATE_JC_MAXACTIVE")) : JC_MAXACTIVE;
}

/* Mittag-Leffler alpha (joint_ml.c): up to three populations, each with its
   own Theta and free ('*') or absent ('0') immigration routes with M
   parameters, uniform M priors U(0, max) with one max per receiving
   population; no growth, skyline segments or divergence */
static boolean
jc_ml_supported (world_fmt *world)
{
  const long numpop = world->numpop;
  long i, e, from, to;
  if (numpop > 3 || world->has_growth || world->timeelements > 2 || world->has_speciation
      || world->bayes->mu)
    return FALSE;
  /* migration out of a daughter (forward in time): not with Mittag-Leffler */
  if (species_mig_gated_model (world))
    return FALSE;
  for (i = 0; i < world->numpop2; i++)
    {
      const char c = world->options->custm2[i];
      if (i < numpop ? c != '*' : strchr ("*0", c) == NULL)
        return FALSE;
      if (c == '*' && world->bayes->map[i][1] != i)
        return FALSE;
      if (i >= numpop && c == '*')
        {
          if (!migration_is_M (world, i) || log_prior_1 == NULL
              || log_prior_1[i] != (logprior1ptr) log_prior_uni1 || world->bayes->minparam[i] != 0.0)
            return FALSE;
          m2mm (i, numpop, &from, &to);
          for (e = mstart (to, numpop); e < mend (to, numpop); e++)
            if (world->options->custm2[e] == '*' && world->bayes->maxparam[e] != world->bayes->maxparam[i])
              return FALSE;
        }
    }
  if (world->tri_mlalpha != FIXED)
    {
      const long npg = world->numparamcumvec[GROWTHPRIOR];
      for (i = 0; i < numpop && i < world->options->mlalphapops_numalloc; i++)
        if (world->options->mlalphapops[i] > 0)
          {
            long r = -1;
            if (pick_parameter (npg + i, &r, world, GROWTHPRIOR) || r < npg || r >= world->numparam || world->bayes->minparam[r] <= 0.0
                || world->bayes->maxparam[r] > 1.0)
              return FALSE;
          }
    }
  if (jc_nactive (world) > jc_maxactive ())
    return FALSE;
  return TRUE;
}

/// TRUE when the model is one the joint combination handles
boolean
jc_supported (world_fmt *world)
{
  long i;
  if (!world->options->bayes_infer)
    return FALSE;
  if (world->has_mlalpha)
    return jc_ml_supported (world);
  /* skyline=PARAM: free ('*') Theta and M multipliers and growth restarted
     per segment, without xNm parameters. Splits (normal or exponential,
     tied or not) use the same split terms as without segments; the
     per-population density of mlf_hastings.c covers them (checked with
     MIGRATE_JC_GATECHECK) */
  if (world->timeelements > 2)
    {
      if (!world->options->skyline_param || world->skymodel2 == NULL)
        return FALSE;
      for (i = 0; i < world->numpop2; i++)
        {
          if (!strchr ("*cg", world->skymodel2[i]))
            return FALSE;
          if (i >= world->numpop && !migration_is_M (world, i))
            return FALSE;
          if (world->skymodel2[i] == '*'
              && (world->bayes->map[i][1] < 0 || world->bayes->map[i][1] >= world->numpop2
                  || strchr ("0c", world->options->custm2[i]) != NULL))
            return FALSE;
        }
    }
  if (world->has_speciation)
    for (i = 0; i < world->species_model_size; i++)
      {
        const species_fmt *sp = &world->species_model[i];
        if (!strchr ("dt", sp->type) || sp->paramindex_mu < 0 || sp->paramindex_mu >= world->numparam)
          return FALSE;   /* (type is lower case: d/D -> 'd', tied t/T -> 't') */
        /* D free, or fixed (a prior of zero width: bayes_fill() maps it out) */
        if (world->bayes->map[sp->paramindex_mu][1] != sp->paramindex_mu
            && (world->bayes->map[sp->paramindex_mu][1] != INVALID
                || world->bayes->maxparam[sp->paramindex_mu] > world->bayes->minparam[sp->paramindex_mu]))
          return FALSE;
      }
  /* migration out of a daughter (forward in time) waits only up to the
     mean split time, which couples its M with D: handled for an entry with
     its own M ('*', not shared) under a uniform prior from 0, without
     skyline segments of that M (Z(G) integrates M at each D in closed form) */
  if (species_mig_gated_model (world))
    for (i = world->numpop; i < world->numpop2; i++)
      {
        long from, to, j;
        if (strchr ("0dDtT", world->options->custm2[i]) != NULL)
          continue;
        m2mm (i, world->numpop, &from, &to);
        if (species_mig_gate (world, from) >= (double) HUGE)
          continue;
        /* skyline=PARAM is fine when the gated M itself has no segments:
           its term uses the lineage time before D over all segments */
        if (world->options->custm2[i] != '*' || world->bayes->map[i][1] != i
            || (world->timeelements > 2 && (world->skymodel2 == NULL || world->skymodel2[i] == '*'))
            || !migration_is_M (world, i) || log_prior_1 == NULL
            || log_prior_1[i] != (logprior1ptr) log_prior_uni1 || world->bayes->minparam[i] != 0.0)
          return FALSE;
        for (j = 0; j < world->numpop2; j++)
          if (j != i && world->bayes->map[j][1] == i)
            return FALSE;
      }
  for (i = 0; i < world->numpop2; i++)
    {
      const char c = world->options->custm2[i];
      /* '*', '0', constant, symmetric, mean (of M or of xNm), letter groups,
         divergence d/D; t/T are recorded as d with a tied species model */
      if (!(strchr ("*0csmSMdDtT", c) || (c >= 'a' && c <= 'z')))
        return FALSE;
    }
  if (world->has_growth)
    for (i = 0; i < world->numpop; i++)
      if (world->options->growpops[i] != 0 && world->bayes->map[i][1] != i)
        return FALSE;   /* a growing population needs its own Theta */
  /* xNm parameters (use-M=NO, S/M) tie a migration parameter to the Theta of
     its target: such a population needs its own Theta and no growth, and
     may be reached by at most one xNm that several entries share (else the
     normalizer is no longer a product of small integrals) */
  {
    const long numpop = world->numpop, numpop2 = world->numpop2;
    long *shared = (long *) mycalloc ((size_t) numpop, sizeof (long));
    long *cnt = (long *) mycalloc ((size_t) numpop2, sizeof (long));
    long e, from, to;
    boolean ok = TRUE;
    for (e = numpop; e < numpop2; e++)
      if (!migration_is_M (world, e) && world->bayes->map[e][1] >= numpop)
        cnt[world->bayes->map[e][1]]++;
    for (e = numpop; e < numpop2 && ok; e++)
      {
        const long r = world->bayes->map[e][1];
        if (migration_is_M (world, e) || r < numpop)
          continue;
        m2mm (e, numpop, &from, &to);
        if (world->bayes->map[to][1] != to
            || (world->has_growth && world->options->growpops[to] != 0))
          ok = FALSE;
        if (cnt[r] > 1)
          {
            if (shared[to] != 0 && shared[to] != r + 1)
              ok = FALSE;
            shared[to] = r + 1;
          }
      }
    myfree (shared);
    myfree (cnt);
    if (!ok)
      return FALSE;
  }
  if (jc_nactive (world) > jc_maxactive ())
    return FALSE;
  return TRUE;
}

static boolean
jc_file_mode (world_fmt *world)
{
  return world->options->has_bayesmdimfile;
}

/// TRUE when the statistics travel in memory (not the bayes-allfile file mode)
boolean
jc_in_memory (world_fmt *world)
{
  return jc_supported (world) && !jc_file_mode (world);
}

static const char *
jc_prefix (world_fmt *world)
{
  (void) world;
  return jc_outfile;
}

static void
jc_filename (world_fmt *world, long locus, long rep, char *name)
{
  snprintf (name, LINESIZE, "%s.joint.%ld.%ld", jc_prefix (world), locus, rep);
}

/* wall clock (debug timing) */
static double
jc_now (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

/* ------------------------------------------------------------- storage */

/* free parameters of the joint combination (with the skyline multipliers);
   fixed for a run */
static boolean *jc_active (world_fmt *world, const jc_layout *ly);
static long jc_nactive_cache = -1;
static long jc_impossible = 0;   /* loci at zero density when the main joint run gave up */

static long
jc_nactive (world_fmt *world)
{
  if (jc_nactive_cache < 0)
    {
      jc_layout ly;
      long p, n = 0;
      jc_layout_make (world, &ly);
      boolean *active = jc_active (world, &ly);
      for (p = 0; p < ly.nphi; p++)
        if (active[p])
          n++;
      myfree (active);
      jc_layout_free (&ly);
      jc_nactive_cache = n;
    }
  return jc_nactive_cache;
}

/* genealogies kept per locus: JC_MAXSAMPLES (sqrt(k) - 1)^2 for k free
   parameters, i.e. 4000 (numpop-1)^2 for a full migration matrix (more
   parameters need more genealogies: three populations at 4000 left most M
   rows with a reweighting ESS below 100), fewer for wide rows (the
   divergence grids) so that a worker's memory stays bounded; the option
   joint-samples replaces both with a fraction of the recorded genealogies
   or a number */
static long
jc_maxsamples (world_fmt *world, long nrow)
{
  long m = (long) (1000000 / (nrow > 0 ? nrow : 1));
  const double r = sqrt ((double) jc_nactive (world)) - 1.0;
  long cap = (long) (JC_MAXSAMPLES * (r > 1.0 ? r * r : 1.0));
  const double js = world->options->joint_samples;
  if (js > 0.0)   /* joint-samples: the user's number, no memory bound */
    {
      const long reps = (world->options->replicate && world->options->replicatenum > 1)
                        ? world->options->replicatenum : 1;
      cap = (long) (js > 1.0 ? js : js * (double) world->options->lsteps * (double) reps);
      m = cap;
    }
  if (getenv ("MIGRATE_JC_MAXSAMPLES") != NULL)   /* experiments */
    {
      cap = atol (getenv ("MIGRATE_JC_MAXSAMPLES"));
      m = cap;
    }
  if (m > cap)
    m = cap;
  if (m < 200)
    m = 200;
  return m;
}

typedef struct
{
  long n;        /* rows stored */
  long seen;     /* samples offered since the last kept row */
  long stride;   /* keep every stride-th sample */
  long cap;
  double *rows;  /* n x nrow */
  long *nspr;    /* per row: its number of spells (Mittag-Leffler) */
  double *sp;    /* the spells of all rows, two doubles each */
  long nsp, capsp;
  double *zc;    /* log Z per stored row (NAN: not yet computed), kept across
                    the chains of a combination (main, error blocks, stepping
                    stones): it depends on neither phi nor the block */
  long nzc;
} jc_buffer;

typedef struct
{
  long loci;
  long nrow;
  jc_buffer *b;
  /* results of the combination, filled by jc_combine() */
  boolean done;
  long np;        /* coordinates: the parameters, then the skyline multipliers */
  long numparam, nseg, nfree;
  long *sky_k;    /* per free skyline entry: its matrix entry */
  boolean *active;
  double *trace;  /* np x nsweeps */
  long nsweeps;   /* length of the main run (trace) */
  boolean has_mcerr;
  double *mcerr;  /* np: Monte Carlo error of the median (split halves) */
  double *mchalf; /* nblock x np: the medians of the blocks */
  long nblock;
  boolean byrep;  /* blocks are replicates */
  boolean has_logc;
  double logc, logc_err; /* joint scaling factor of the marginal likelihood */
  boolean has_boot;
  double *boot_err;      /* np: bootstrap error of the median (genealogy sampling) */
  boolean *boot_flag;    /* np: large genealogy-sampling error (a warning: Joint is kept) */
  double boot_ess;       /* median reweighting ESS of the JC_BOOTK trace points */
  double logc_boot;      /* bootstrap error of the scaling factor */
  long nloci_used, tmin, tmax;
} jc_store;

static FILE *jc_file = NULL;
static long jc_file_locus = -1;
static long jc_file_rep = -1;

static jc_store *
jc_get_store (world_fmt *world, long nrow)
{
  jc_store *js = (jc_store *) world->jointstats;
  long l;
  if (js == NULL)
    {
      js = (jc_store *) mycalloc (1, sizeof (jc_store));
      js->loci = world->loci;
      js->nrow = nrow;
      js->b = (jc_buffer *) mycalloc ((size_t) world->loci, sizeof (jc_buffer));
      for (l = 0; l < world->loci; l++)
        js->b[l].stride = 1;
      world->jointstats = js;
    }
  return js;
}

static void
jc_append (jc_store *js, long locus, const double *row, const double *sp, long nsp)
{
  jc_buffer *b = &js->b[locus];
  if (b->n >= b->cap)
    {
      b->cap = b->cap ? 2 * b->cap : 64;
      b->rows = (double *) myrealloc (b->rows, sizeof (double) * (size_t) (b->cap * js->nrow));
      b->nspr = (long *) myrealloc (b->nspr, sizeof (long) * (size_t) b->cap);
    }
  memcpy (b->rows + b->n * js->nrow, row, sizeof (double) * (size_t) js->nrow);
  b->nspr[b->n] = nsp;
  if (nsp > 0)
    {
      if (b->nsp + nsp > b->capsp)
        {
          b->capsp = 2 * (b->nsp + nsp) + 1024;
          b->sp = (double *) myrealloc (b->sp, sizeof (double) * (size_t) (2 * b->capsp));
        }
      memcpy (b->sp + 2 * b->nsp, sp, sizeof (double) * (size_t) (2 * nsp));
      b->nsp += nsp;
    }
  b->n++;
}

/* log of (e^{g t1} - e^{g t0}) / g (= t1 - t0 for g = 0), without
   overflow: g t reaches thousands for microsatellite time scales */
static double
jc_log_growth_dt (double g, double t0, double t1)
{
  const double dt = t1 - t0;
  if (dt <= 0.0)
    return -HUGE_VAL;
  if (fabs (g) < 1e-12)
    return log (dt);
  if (g > 0.0)
    return g * t0 + log (expm1 (g * dt)) - log (g);
  return g * t0 + log (-expm1 (g * dt)) - log (-g);
}

/* log(e^a + e^b) */
static double
jc_logadd (double a, double b)
{
  if (a == -HUGE_VAL)
    return b;
  if (b == -HUGE_VAL)
    return a;
  return a > b ? a + log1p (exp (b - a)) : b + log1p (exp (a - b));
}

/* integrated crossing time, while a genealogy is recorded: the crossing
   lineage of split k sits in the ancestor on all of its segment (a, b) in
   the row ("base row"); jc_rb_count() gives the lineage counts of interval
   i of the time list accordingly, jc_div_exact() leaves out the lineage's
   split terms (the unit JC_U_SPLITRB carries them, integrated over c) */
typedef struct
{
  long k;          /* split, -1: none */
  long ie;         /* time-list index of its 'd' event */
  double a, c0, b; /* segment and the sampled crossing time */
  long d, an;      /* daughter, ancestor */
  long n;          /* crossings that move (sharp split: all of the split) */
  double cross[JC_SHARPMAXN]; /* their sampled crossing times */
} jc_rbcur_t;
static jc_rbcur_t jc_rbc = { -1, -1, 0.0, 0.0, 0.0, -1, -1, 0, { 0.0 } };

static long
jc_rb_count (vtlist *tl, long i, long pop)
{
  const long k = tl[i].lineages[pop];
  long j, m = 0;
  if (jc_rbc.k < 0 || tl[i - 1].age < jc_rbc.a || (pop != jc_rbc.d && pop != jc_rbc.an))
    return k;
  for (j = 0; j < jc_rbc.n; j++)   /* moving lineages still in the daughter here */
    if (tl[i].age <= jc_rbc.cross[j])
      m++;
  return pop == jc_rbc.d ? k - m : k + m;
}

/* the time of a neighbour node: the half records of a migration or split
   node can disagree (one half kept an old time), the top half is right */
static double
jc_ntime (node *n)
{
  return showtop (n)->tyme;
}

/* finds the crossing to integrate (the first split with exactly one
   crossing whose segment has at most JC_RBMAXP pieces) and fills the row's
   block: per piece of (a, b) its start and the slopes 2 k_d / inh and
   2 k_a / inh of the other lineages' counts (sets jc_rbc) */
static void
jc_rb_find (world_fmt *world, const jc_layout *ly, double *st, double inh)
{
  const long T = world->treetimes->T;
  vtlist *tl = world->treetimes->tl;
  double *R = st + ly->off_R;
  long k, i;
  jc_rbc.k = -1;
  R[0] = -1.0;
  for (k = 0; k < ly->nsplit; k++)
    {
      if (!ly->rb_ok[k])
        continue;
      species_fmt *sp = &world->species_model[ly->split_model[k]];
      long n = 0, ie = -1;
      for (i = 1; i < T - 1; i++)
        if (tl[i].eventnode->type == 'd'
            && get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                        world->species_model, world->species_model_size) == sp)
          {
            n++;
            ie = i;
          }
      if (n != 1)
        continue;
      node *E = tl[ie].eventnode, *n1 = E->back, *n2 = E->next->back;


      if (n1 == NULL || n2 == NULL)
        continue;
      const double t1 = jc_ntime (n1), t2 = jc_ntime (n2);
      node *up = t1 > t2 ? n1 : n2;
      const double c0 = tl[ie].age, a = t1 > t2 ? t2 : t1;
      double b = t1 > t2 ? t1 : t2;
      if (up->type == 'r' || !(a <= c0 && c0 <= b))
        continue;
      if (world->species_model_dist == NORMAL_DIST && b > ly->split_b1[k])
        b = ly->split_b1[k];
      if (!(b > a))
        continue;
      const long d = ly->split_to[k], an = ly->split_from[k];
      long np = 0;
      boolean ok = TRUE;
      for (i = 1; i < T && ok; i++)
        {
          const double t0 = tl[i - 1].age, t1 = tl[i].age;
          const double u0 = t0 > a ? t0 : a, v = t1 < b ? t1 : b;
          if (!(v > u0))
            continue;
          const long kd = tl[i].lineages[d] - (t1 <= c0 ? 1 : 0);
          const long ka = tl[i].lineages[an] - (t0 >= c0 ? 1 : 0);
          if (kd < 0 || ka < 0)
            {
              ok = FALSE;
              break;
            }
          const double cd = 2.0 * (double) kd / inh, ca = 2.0 * (double) ka / inh;
          /* skyline=PARAM: a piece lies inside one segment */
          long sg = jc_seg (ly, u0);
          double u = u0;
          while (u < v)
            {
              const double ve = (sg + 1 < ly->nseg && ly->segt[sg + 1] < v) ? ly->segt[sg + 1] : v;
              if (!(np > 0 && R[4 + 3 * (np - 1) + 1] == cd && R[4 + 3 * (np - 1) + 2] == ca
                    && jc_seg (ly, R[4 + 3 * (np - 1)]) == sg))
                {   /* else: same counts and segment, one piece */
                  if (np == JC_RBMAXP)
                    {
                      ok = FALSE;
                      break;
                    }
                  R[4 + 3 * np] = u;
                  R[4 + 3 * np + 1] = cd;
                  R[4 + 3 * np + 2] = ca;
                  np++;
                }
              u = ve;
              sg++;
            }
        }
      if (!ok || np == 0)
        continue;
      R[0] = (double) k;
      R[1] = a;
      R[2] = b;
      R[3] = (double) np;
      jc_rbc.k = k;
      jc_rbc.ie = ie;
      jc_rbc.a = a;
      jc_rbc.c0 = c0;
      jc_rbc.n = 1;
      jc_rbc.cross[0] = c0;
      JC_RBN (R) = 1.0;
      jc_rbc.b = b;
      jc_rbc.d = d;
      jc_rbc.an = an;
      return;
    }
}

/* sharp splits: per split its crossing lineages (n), their mean crossing
   time c0 (the sampled genealogy is the base row) and the window (A, B):
   above every crossing lineage's lower event, below every upper one; a
   neighbour that is itself a crossing of another split moves too and gives
   instead the order c(child) <= c(parent) ("chain"). Then the lineage
   counts of every population on the time range that the crossings can
   cover */
static void
jc_sharp_find (world_fmt *world, const jc_layout *ly, double *st)
{
  const long T = world->treetimes->T, P = ly->numpop, K = ly->nsplit;
  vtlist *tl = world->treetimes->tl;
  double *R = st + ly->off_R;
  double tlo = HUGE_VAL, thi = -HUGE_VAL;
  long i, k, p, ni = 0;
  boolean any = FALSE;
  R[0] = -1.0;
  R[1] = 0.0;
  R[3] = world->options->inheritance_scalars[world->locus];
  for (k = 0; k < K; k++)
    {
      double *Sk = R + JC_SH_SPLIT (k);
      Sk[0] = -1e300;
      Sk[1] = 1e300;
      Sk[2] = 0.0;
      Sk[3] = 0.0;
      Sk[4] = 0.0;
    }
  for (i = 1; i < T - 1; i++)
    if (tl[i].eventnode->type == 'd')
      {
        species_fmt *sp = get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                                   world->species_model, world->species_model_size);
        node *E = tl[i].eventnode, *n1 = E->back, *n2 = E->next->back;
        for (k = 0; k < K && &world->species_model[ly->split_model[k]] != sp; k++)
          ;
        if (k == K || n1 == NULL || n2 == NULL)
          continue;
        double *Sk = R + JC_SH_SPLIT (k);
        const double t1 = jc_ntime (n1), t2 = jc_ntime (n2);
        node *up = t1 > t2 ? n1 : n2, *dn = t1 > t2 ? n2 : n1;
        const double tup = t1 > t2 ? t1 : t2, tdn = t1 > t2 ? t2 : t1;
        if (dn->type != 'd' && tdn > Sk[0])
          Sk[0] = tdn;
        if (up->type == 'd')
          Sk[4] = 1.0;
        else if (tup < Sk[1])
          Sk[1] = tup;
        Sk[2] += 1.0;
        Sk[3] += tl[i].age;
      }
  /* every other event in a daughter population (a migration out of it, a
     coalescence, with migration after the split) keeps the split above it:
     the daughter holds no lineage past its split */
  for (i = 1; i < T - 1; i++)
    {
      const node *E = tl[i].eventnode;
      if (E->type == 'd')
        continue;
      for (k = 0; k < K; k++)
        {
          double *Sk = R + JC_SH_SPLIT (k);
          if ((E->pop == ly->split_to[k] || E->actualpop == ly->split_to[k]) && tl[i].age > Sk[0])
            Sk[0] = tl[i].age;
        }
    }
  for (k = 0; k < K; k++)
    {
      double *Sk = R + JC_SH_SPLIT (k);
      if (Sk[2] == 0.0)
        continue;
      any = TRUE;
      Sk[3] /= Sk[2];
      const double lo = Sk[0] > ly->mlo[k] ? Sk[0] : ly->mlo[k], hi = Sk[1] < ly->mhi[k] ? Sk[1] : ly->mhi[k];
      if (lo < tlo) tlo = lo;
      if (Sk[3] < tlo) tlo = Sk[3];
      if (hi > thi) thi = hi;
      if (Sk[3] > thi) thi = Sk[3];
    }
  if (!any)
    {   /* no crossing (with migration after the split every lineage can
           leave the daughter before it): D only has to lie above the
           daughter's events */
      for (k = 0; k < K; k++)
        if (R[JC_SH_SPLIT (k)] > ly->mlo[k])
          {
            R[0] = 0.0;
            R[1] = 0.0;
            R[2] = 0.0;
            return;
          }
      return;
    }
  R[2] = tlo;
  double *tlr = R + JC_SH_TL (ly);
  for (i = 1; i < T; i++)
    {
      const double t0 = tl[i - 1].age, t1 = tl[i].age;
      if (t1 <= tlo || t0 >= thi)
        continue;
      const double v = t1 < thi ? t1 : thi;
      boolean same = ni > 0;
      for (p = 0; p < P && same; p++)
        same = tlr[(ni - 1) * (1 + P) + 1 + p] == (double) tl[i].lineages[p];
      if (same)
        {
          tlr[(ni - 1) * (1 + P)] = v;
          continue;
        }
      if (ni == JC_SHARPMAXI)
        {
          R[0] = -2.0;
          return;
        }
      tlr[ni * (1 + P)] = v;
      for (p = 0; p < P; p++)
        tlr[ni * (1 + P) + 1 + p] = (double) tl[i].lineages[p];
      ni++;
    }
  R[0] = 0.0;
  R[1] = (double) ni;
}

/* the D blocks of sharp splits: empty (no split event, flat); the split
   terms are in the unit JC_U_SPLITRB */
static void
jc_sharp_clear (const jc_layout *ly, double *st)
{
  long j, k;
  if (!ly->sharp)
    return;
  for (k = 0; k < ly->nsplit; k++)
    {
      st[ly->off_D[k]] = 0.0;
      st[ly->off_D[k] + 1] = 0.0;
      for (j = 0; j < ly->split_ns[k] * JC_NM; j++)
        st[ly->off_D[k] + 2 + j] = 0.0;
    }
}

/* sharp splits: are the crossing times c inside this row's windows and in
   the order its chains require? */
static boolean
jc_sharp_ok (const jc_layout *ly, const double *R, const double *c)
{
  long k;
  for (k = 0; k < ly->nsplit; k++)
    {
      const double *Sk = R + JC_SH_SPLIT (k);
      if (c[k] < Sk[0])   /* also without crossings: above the daughter's events */
        return FALSE;
      if (Sk[2] == 0.0)
        continue;
      if (c[k] > Sk[1])
        return FALSE;
      if (Sk[4] > 0.0 && ly->sharp_par[k] >= 0 && c[k] > c[ly->sharp_par[k]])
        return FALSE;
    }
  return TRUE;
}

/* sharp splits: the changes of the coalescence integrals K (divided by the
   inheritance scalar) and lineage times S of every population when the
   crossings of each split move from c0 to c. With growth (gp: the growth rate
   of every population, NULL without) a piece (a, b) of a growing population
   adds its count change times (e^{gb} - e^{ga}) / g, the growth-weighted
   integral of int k(k-1) e^{gt} dt */
static void
jc_sharp_delta (const jc_layout *ly, const double *R, const double *c, const double *gp, double *dK, double *dS)
{
  const long P = ly->numpop, K = ly->nsplit, ni = (long) R[1];
  const double *tlr = R + JC_SH_TL (ly);
  long i, k, p, q, nc;
  double u = R[2], cut[64];
  for (p = 0; p < P; p++)
    dK[p] = dS[p] = 0.0;
  for (i = 0; i < ni; i++)
    {
      const double v = tlr[i * (1 + P)];
      const double *cnt = tlr + i * (1 + P) + 1;
      nc = 0;
      for (k = 0; k < K && nc < 62; k++)
        {
          const double *Sk = R + JC_SH_SPLIT (k);
          if (Sk[2] == 0.0)
            continue;
          if (Sk[3] > u && Sk[3] < v) cut[nc++] = Sk[3];
          if (c[k] > u && c[k] < v) cut[nc++] = c[k];
        }
      for (q = 1; q < nc; q++)
        {
          long j;
          for (j = q; j > 0 && cut[j - 1] > cut[j]; j--)
            {
              const double t = cut[j];
              cut[j] = cut[j - 1];
              cut[j - 1] = t;
            }
        }
      cut[nc++] = v;
      double a = u;
      for (q = 0; q < nc; q++)
        {
          const double b = cut[q], len = b - a, m = 0.5 * (a + b);
          if (len > 0.0)
            {
              double add[64];
              for (p = 0; p < P; p++)
                add[p] = 0.0;
              for (k = 0; k < K; k++)
                {
                  const double *Sk = R + JC_SH_SPLIT (k);
                  double sgn = 0.0;
                  if (Sk[2] == 0.0)
                    continue;
                  if (Sk[3] <= m && m < c[k])
                    sgn = 1.0;
                  else if (c[k] <= m && m < Sk[3])
                    sgn = -1.0;
                  add[ly->split_to[k]] += Sk[2] * sgn;
                  add[ly->split_from[k]] -= Sk[2] * sgn;
                }
              for (p = 0; p < P; p++)
                if (add[p] != 0.0)
                  {
                    const double n0 = cnt[p], n1 = cnt[p] + add[p];
                    const double g = gp != NULL ? gp[p] : 0.0;
                    const double w = g != 0.0 ? exp (g * a) * expm1 (g * len) / g : len;
                    dK[p] += w * (n1 * (n1 - 1.0) - n0 * (n0 - 1.0));
                    dS[p] += len * add[p];
                  }
            }
          a = b;
        }
      u = v;
    }
  for (p = 0; p < P; p++)
    dK[p] /= R[3];
}

static double jc_value (const jc_layout *ly, const double *phi, long e);
static double jc_rate (const jc_layout *ly, const double *phi, long e);

/* the log-term change of dK, dS at phi: every population and entry (mode
   0), or only the fixed or inactive ones (mode 2, the terms the normalizer
   does not integrate) */
static double
jc_sharp_lin (const jc_layout *ly, const boolean *active, const double *phi, const double *dK,
              const double *dS, double mu, int mode)
{
  long p, e, from, to;
  double v = 0.0;
  for (p = 0; p < ly->numpop; p++)
    if (dK[p] != 0.0)
      {
        if (mode == 2 && ly->rep[p] >= 0 && active[ly->rep[p]])
          continue;
        const double th = jc_value (ly, phi, p);
        if (!(th > 0.0))
          return -HUGE_VAL;
        v -= dK[p] / (mu * th);
      }
  for (e = ly->numpop; e < ly->numpop2; e++)
    {
      if (ly->mzero[e])
        continue;
      m2mm (e, ly->numpop, &from, &to);
      if (dS[to] == 0.0 || (mode == 2 && ly->rep[e] >= 0 && active[ly->rep[e]]))
        continue;
      const double x = jc_rate (ly, phi, e);
      if (x > 0.0)
        v -= x * dS[to] / mu;
    }
  return v;
}

/* sharp splits: the row st with every split's crossings moved to c, into
   mv (the normalizer's input): coalescence integrals K, lineage times S and,
   with growth, the growth integrals on every point of the growth grid,
   log(e^A + delta) (the moved total stays positive). dK, dS: the plain
   changes */
static void
jc_sharp_move_row (const jc_layout *ly, const double *st, const double *c, double *mv, double *dK, double *dS)
{
  const double *R = st + ly->off_R;
  long p;
  jc_sharp_delta (ly, R, c, NULL, dK, dS);
  memcpy (mv, st, sizeof (double) * (size_t) ly->nrow);
  for (p = 0; p < ly->numpop; p++)
    {
      mv[ly->off_K + p] += dK[p];
      mv[ly->off_S + p] += dS[p];
    }
  if (ly->nslot > 0)
    {
      double gp[64], dKg[64], dSg[64];
      long j;
      for (j = 0; j < JC_NG; j++)
        {
          for (p = 0; p < ly->numpop; p++)
            gp[p] = ly->pop_slot[p] >= 0 ? ly->ggrid[ly->pop_slot[p] * JC_NG + j] : 0.0;
          jc_sharp_delta (ly, R, c, gp, dKg, dSg);
          for (p = 0; p < ly->numpop; p++)
            if (ly->pop_slot[p] >= 0 && dKg[p] != 0.0)
              {
                double *A = mv + ly->off_A[p] + j;
                const double tot = exp (*A) + dKg[p];
                *A = tot > 0.0 ? log (tot) : -HUGE_VAL;
              }
        }
    }
}

/* sharp splits (unit JC_U_SPLITRB): the row with every split's crossings at
   its D, relative to the sampled row; zero outside the windows */
static double
jc_sharp_term (const jc_layout *ly, const double *st, const double *phi, double mu)
{
  const double *R = st + ly->off_R;
  double c[64], dK[64], dS[64];
  long k;
  if (R[0] == -1.0)
    return 0.0;
  if (R[0] < 0.0)
    return JC_DOUT;
  for (k = 0; k < ly->nsplit; k++)
    c[k] = phi[ly->split_pmu[k]];
  if (!jc_sharp_ok (ly, R, c))
    return JC_DOUT;
  double gp[64];
  long p;
  for (p = 0; p < ly->numpop; p++)
    gp[p] = ly->pop_slot[p] >= 0 ? phi[ly->slot_param[ly->pop_slot[p]]] : 0.0;
  jc_sharp_delta (ly, R, c, gp, dK, dS);
  for (p = 0; p < ly->numpop; p++)
    if (ly->sharp_pu[p])
      dK[p] = 0.0;   /* in the population's own unit (jc_pop_term_phi) */
  const double v = jc_sharp_lin (ly, NULL, phi, dK, dS, mu, 0);
  return v > JC_DOUT ? v : JC_DOUT;
}

/* window of the local mean grid of split k, std row b: the genealogy's
   split-event centre +- JC_DZ std, clipped to the prior range (the whole
   range when the genealogy has no split event for it) */
static void
jc_row_window (const jc_layout *ly, const double *st, long k, long b, double *L, double *U)
{
  const double lo = ly->mlo[k], hi = ly->mhi[k];
  const double c = st[ly->off_D[k]], nd = st[ly->off_D[k] + 1];
  const double w = JC_DZ * ly->sgrid[k * JC_NS + b];
  *L = lo;
  *U = hi;
  (void) nd;
  if (c > 0.0)
    {
      if (c - w > lo)
        *L = c - w;
      if (c + w < hi)
        *U = c + w;
      if (*U - *L < 1e-12 * (hi - lo))
        {
          *L = lo;
          *U = hi;
        }
    }
}

/* divergence log-term of every split for the current genealogy: the waiting
   terms of the descendant population's lineages and the split-event
   densities, exactly as probg_treetimes_local() adds them */
static double
jc_div_exact (world_fmt *world, long sm, double mu, double sigma)
{
  const long T = world->treetimes->T;
  vtlist *tl = world->treetimes->tl;
  species_fmt *sp = &world->species_model[sm];
  double D = 0.0;
  long i;
  for (i = 1; i < T - 1; i++)
    {
      const double t0 = tl[i - 1].age, t1 = tl[i].age;
      /* an interval ending at a (dated) tip has waiting terms only */
      const long kto = jc_rb_count (tl, i, sp->to);
      if (kto > 0)
        D += (double) kto * (*log_prob_wait_speciate) (t0, t1, mu, sigma, sp);
      if (i == jc_rbc.ie && jc_rbc.k >= 0)
        continue;   /* the integrated crossing */
      if (tl[i].eventnode->type == 'd'
          && get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                      world->species_model, world->species_model_size) == sp)
        D += (*log_point_prob_speciate) (world->species_model_dist == NORMALSHORTCUT_DIST ? 0.5 * (t0 + t1) : t1,
                                         mu, sigma, sp);
    }
  return D;
}

/* per split: the centre and count of this genealogy's split events, then
   D on each std row's local mean grid */
static void
jc_record_divergence (world_fmt *world, const jc_layout *ly, double *st)
{
  const long T = world->treetimes->T;
  vtlist *tl = world->treetimes->tl;
  long i, k, a, b;
  for (k = 0; k < ly->nsplit; k++)
    {
      species_fmt *sp = &world->species_model[ly->split_model[k]];
      double sum = 0.0, n = 0.0;
      for (i = 1; i < T - 1; i++)
        if (tl[i].eventnode->type == 'd' && !(i == jc_rbc.ie && jc_rbc.k >= 0)
            && get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                        world->species_model, world->species_model_size) == sp)
          {
            sum += tl[i].age;
            n += 1.0;
          }
      if (n == 0.0)
        {   /* no split event: centre on the last time a lineage is in the
               derived population (the term falls off below it, is flat above) */
          for (i = 1; i < T - 1; i++)
            if (jc_rb_count (tl, i, sp->to) > 0 && tl[i].age > sum)
              sum = tl[i].age;
          st[ly->off_D[k]] = sum;
        }
      else
        st[ly->off_D[k]] = sum / n;
      st[ly->off_D[k] + 1] = n;
      double *D = st + ly->off_D[k] + 2;
      const long ns = ly->split_ns[k];
      for (b = 0; b < ns; b++)
        {
          double L, U;
          jc_row_window (ly, st, k, b, &L, &U);
          for (a = 0; a < JC_NM; a++)
            D[b * JC_NM + a] = jc_div_exact (world, ly->split_model[k],
                                             L + (U - L) * (double) a / (JC_NM - 1),
                                             ly->sgrid[k * JC_NS + b]);
        }
    }
}

static void jc_rb_check (world_fmt *world, const jc_layout *ly, const double *st);

extern void construct_tymelist (world_fmt * world, timelist_fmt * timevector);

/* MIGRATE_JC_SHARPCHECK: the sharp terms against the sampler. The crossings
   of every split of the cold-chain genealogy are moved to the sampler's D0
   and to a c' inside the row's window (D set to c'); the difference of
   probg_treetimes() between the two must equal the joint's full term at
   D = c' minus D = D0 (the split-event spikes are the same in both, every
   crossing exactly at its D). docs/joint_split_rb_design.md, "Sharp splits
   with tied splits, growth and gated migration" */
static void
jc_sharp_check (world_fmt *world, const jc_layout *ly, const double *st)
{
  static double worst = 0.0;
  static long nchk = 0;
  const double *R = st + ly->off_R;
  const long K = ly->nsplit, T = world->treetimes->T;
  vtlist *tl = world->treetimes->tl;
  node **nd;
  double *old, cnew[64], d0[64];
  long nn = 0, i, k, j;
  if (!ly->sharp || R[0] != 0.0 || K > 64)
    return;
  nd = (node **) mycalloc ((size_t) T, sizeof (node *));
  old = (double *) mycalloc ((size_t) T, sizeof (double));
  long *ks = (long *) mycalloc ((size_t) T, sizeof (long));
  for (i = 1; i < T - 1; i++)
    if (tl[i].eventnode->type == 'd')
      {
        species_fmt *sp = get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                                   world->species_model, world->species_model_size);
        for (k = 0; k < K && &world->species_model[ly->split_model[k]] != sp; k++)
          ;
        if (k == K)
          continue;
        nd[nn] = showtop (tl[i].eventnode);
        old[nn] = nd[nn]->tyme;
        ks[nn] = k;
        nn++;
      }
  /* c' per split group: uniform in the window (parents first via min) */
  for (k = 0; k < K; k++)
    {
      const double *Sk = R + JC_SH_SPLIT (k);
      d0[k] = world->param0[ly->split_pmu[k]];
      const double lo = Sk[0] > ly->mlo[k] ? Sk[0] : ly->mlo[k];
      const double hi = Sk[2] > 0.0 && Sk[1] < ly->mhi[k] ? Sk[1] : ly->mhi[k];
      /* its own generator: the check must not change the chain */
      static unsigned long long xs = 88172645463325252ULL;
      xs ^= xs << 13; xs ^= xs >> 7; xs ^= xs << 17;
      cnew[k] = ly->split_rep[k] == k ? lo + (hi - lo) * ((double) (xs >> 11) / 9007199254740992.0) : cnew[ly->split_rep[k]];
    }
  for (k = 0; k < K; k++)
    if (ly->split_rep[k] != k)
      cnew[k] = cnew[ly->split_rep[k]];
  if (!jc_sharp_ok (ly, R, cnew) || !jc_sharp_ok (ly, R, d0))
    goto done;
  double p[2], jt[2];
  boolean *act = (boolean *) mycalloc ((size_t) ly->nphi, sizeof (boolean));
  double *phi = (double *) mycalloc ((size_t) ly->nphi, sizeof (double));
  for (j = 0; j < ly->nphi; j++)
    act[j] = TRUE;
  if (ly->rate)
    act[ly->rate_p] = FALSE;
  long side;
  for (side = 0; side < 2; side++)
    {
      const double *c = side == 0 ? d0 : cnew;
      for (j = 0; j < nn; j++)
        adjust_time_all (nd[j], c[ks[j]]);
      for (k = 0; k < K; k++)
        {   /* a tied partner keeps its own copy of D (speciate.c) */
          world->param0[ly->split_pmu[k]] = c[k];
          world->param0[world->species_model[ly->split_model[k]].paramindex_mu] = c[k];
          world->species_model[ly->split_model[k]].mu = c[k];
        }
      construct_tymelist (world, &world->treetimes[0]);
      p[side] = probg_treetimes (world);
      memcpy (phi, world->param0, sizeof (double) * (size_t) ly->np);
      jt[side] = jc_full_term (ly, act, st, NULL, phi, st[ly->off_mu]);
      if (side == 1)
        {   /* the normalizer's moved row against the moved genealogy's own
               statistics (K, S, growth integrals) */
          static double wK = 0.0, wA = 0.0;
          double *mv = (double *) mycalloc ((size_t) ly->nrow, sizeof (double));
          double dK[64], dS[64], Kd[64], Sd[64];
          const double inh = R[3];
          vtlist *tl2 = world->treetimes->tl;
          const long T2 = world->treetimes->T;
          long pp, ii, jj;
          jc_sharp_move_row (ly, st, cnew, mv, dK, dS);
          for (pp = 0; pp < ly->numpop; pp++)
            Kd[pp] = Sd[pp] = 0.0;
          for (ii = 1; ii < T2 - 1; ii++)
            for (pp = 0; pp < ly->numpop; pp++)
              {
                const double k = (double) jc_rb_count (tl2, ii, pp), dt = tl2[ii].age - tl2[ii - 1].age;
                Kd[pp] += k * (k - 1.0) * dt / inh;
                Sd[pp] += k * dt;
              }
          for (pp = 0; pp < ly->numpop; pp++)
            {
              const double eK = fabs (mv[ly->off_K + pp] - Kd[pp]) + fabs (mv[ly->off_S + pp] - Sd[pp]);
              if (eK > wK) wK = eK;
              if (ly->pop_slot[pp] >= 0)
                for (jj = 0; jj < JC_NG; jj += 16)
                  {
                    const double gg = ly->ggrid[ly->pop_slot[pp] * JC_NG + jj];
                    double a = -HUGE_VAL;
                    for (ii = 1; ii < T2 - 1; ii++)
                      {
                        const double k = (double) jc_rb_count (tl2, ii, pp);
                        if (k > 1.0 && tl2[ii].age > tl2[ii - 1].age)
                          a = jc_logadd (a, log (k * (k - 1.0)) + jc_log_growth_dt (gg, tl2[ii - 1].age, tl2[ii].age));
                      }
                    a -= log (inh);
                    const double eA = fabs (mv[ly->off_A[pp] + jj] - a);
                    if (eA > wA) wA = eA;
                  }
            }
          fprintf (stderr, "JC_SHARPROW worst |K|+|S| %.3g, worst log growth integral %.3g\n", wK, wA);
          myfree (mv);
        }
    }
  for (j = 0; j < nn; j++)
    adjust_time_all (nd[j], old[j]);
  for (k = 0; k < K; k++)
    {
      world->param0[ly->split_pmu[k]] = d0[k];
      world->param0[world->species_model[ly->split_model[k]].paramindex_mu] = d0[k];
      world->species_model[ly->split_model[k]].mu = d0[k];
    }
  construct_tymelist (world, &world->treetimes[0]);
  {
    const double diff = (p[1] - p[0]) - (jt[1] - jt[0]);
    nchk++;
    if (fabs (diff) > worst)
      worst = fabs (diff);
    fprintf (stderr, "JC_SHARPCHECK locus %ld: sampler %.6f joint %.6f diff %.3g (worst %.3g, %ld checks)\n",
             world->locus + 1, p[1] - p[0], jt[1] - jt[0], diff, worst, nchk);
  }
  myfree (phi);
  myfree (act);
 done:
  myfree (nd);
  myfree (old);
  myfree (ks);
}

/// Per-genealogy statistics of the current cold-chain genealogy (called
/// right after bayes_save()).
/* grid rate k of the shared-rate grid */
static double
jc_rate_x (const jc_layout *ly, long k)
{
  return exp (ly->rlo + (ly->rhi - ly->rlo) * (double) k / (double) (ly->nrk - 1));
}

/* log p(D|G,r) at grid position u (0..K-1, in log r): cubic through the
   four nearest grid points (Lagrange; at the ends the first/last four) */
static double
jc_rate_interp (const double *x, long K, double u)
{
  long j;
  double a, b, c, d;
  if (u <= 0.0)
    return x[0];
  if (u >= (double) (K - 1))
    return x[K - 1];
  j = (long) floor (u) - 1;
  if (j < 0)
    j = 0;
  if (j > K - 4)
    j = K - 4;
  u -= (double) j;   /* the points are at 0, 1, 2, 3 */
  a = -(u - 1.0) * (u - 2.0) * (u - 3.0) / 6.0;
  b = u * (u - 2.0) * (u - 3.0) / 2.0;
  c = -u * (u - 1.0) * (u - 3.0) / 2.0;
  d = u * (u - 1.0) * (u - 2.0) / 6.0;
  return a * x[j] + b * x[j + 1] + c * x[j + 2] + d * x[j + 3];
}

/* the shared rate: log p(D|G,r_k) of this genealogy on the rate grid (the
   genealogy in reference time, its node ages rescaled to each r_k as the
   sampler's rate move does), then log Lambda(G) = log of the integral of
   p(D|G,r) p_s(r) dr over the grid (p_s: the sampling prior of the rate;
   trapezoid in log r). The tree is restored afterwards. */
static void
jc_record_rate (world_fmt *world, const jc_layout *ly, double *st)
{
  const long K = ly->nrk, locus = world->locus;
  const double r0 = world->options->mu_rates[locus];
  const double like0 = world->likelihood[world->G];
  double *x = st + ly->off_lr;
  double *v = (double *) mycalloc ((size_t) K, sizeof (double));
  double cur = r0, h = (ly->rhi - ly->rlo) / (double) (K - 1);
  long k;
  if (!ly->rate)
    {
      myfree (v);
      return;
    }
  for (k = 0; k < K; k++)
    {
      const double r = jc_rate_x (ly, k);
      world->options->mu_rates[locus] = r;
      world->options->lmu_rates[locus] = log (r);
      recalc_timelist (world, r, cur);
      cur = r;
      /* NODATA: the sampler ignores the data likelihood (prioralone) */
      x[k] = world->options->prioralone ? 0.0 : world->likelihood[world->G];
    }
  world->options->mu_rates[locus] = r0;
  world->options->lmu_rates[locus] = log (r0);
  recalc_timelist (world, r0, cur);
  {   /* Lambda: trapezoid in log r on JC_RATESUB points per grid interval,
         p(D|G,r) from the interpolation */
    const long n = (K - 1) * JC_RATESUB + 1;
    double m = -HUGE_VAL, sum = 0.0;
    v = (double *) myrealloc (v, sizeof (double) * (size_t) n);
    for (k = 0; k < n; k++)
      {
        const double u = (double) k / JC_RATESUB, lr = ly->rlo + u * h, r = exp (lr);
        v[k] = jc_rate_interp (x, K, u) + scaling_prior_joint (world, ly->rate_p, r) + lr;
        if (v[k] > m)
          m = v[k];
      }
    for (k = 0; k < n; k++)
      sum += (k == 0 || k == n - 1 ? 0.5 : 1.0) * exp (v[k] - m);
    x[K] = m + log (sum * h / JC_RATESUB);
  }
  if (getenv ("MIGRATE_JC_RATECHECK") != NULL)
    {   /* restored likelihood; interpolated log L at the sampled rate */
      static double worst = 0.0, worsti = 0.0;
      const double li = jc_rate_interp (x, K, (log (r0) - ly->rlo) / h);
      if (fabs (world->likelihood[world->G] - like0) > worst)
        worst = fabs (world->likelihood[world->G] - like0);
      if (fabs (li - like0) > worsti)
        worsti = fabs (li - like0);
      fprintf (stderr, "JC_RATECHECK locus %ld r %g: restore |dL| worst %g, interpolated - true %g (worst %g)\n",
               locus + 1, r0, worst, li - like0, worsti);
    }
  myfree (v);
}

/* the shared-rate term of a genealogy: log p(D|G,r) - log Lambda(G),
   p(D|G,r) cubic in log r between the grid points (clamped at the ends;
   the joint keeps r inside the grid) */
static double
jc_rate_term (const jc_layout *ly, const double *st, double r)
{
  const double *x = st + ly->off_lr;
  const long K = ly->nrk;
  return jc_rate_interp (x, K, (log (r) - ly->rlo) / (ly->rhi - ly->rlo) * (double) (K - 1)) - x[K];
}

void
jc_record_sample (world_fmt *world)
{
  static double *spells = NULL;
  static long spalloc = 0;
  long nsp = 0;
  if (!jc_supported (world))
    return;
  jc_layout ly;
  jc_layout_make (world, &ly);
  const long numpop = world->numpop;
  const long T = world->treetimes->T;
  vtlist *tl = world->treetimes->tl;
  double *st = (double *) mycalloc ((size_t) ly.nrow, sizeof (double));
  long i, pop, g;
  for (pop = 0; pop < numpop; pop++)   /* log A_i(g), accumulated below */
    if (ly.pop_slot[pop] >= 0)
      {
        for (g = 0; g < JC_NG; g++)
          st[ly.off_A[pop] + g] = -HUGE_VAL;
        if (ly.off_sA[pop] >= 0)
          for (g = 0; g < ly.nseg * JC_NG; g++)
            st[ly.off_sA[pop] + g] = -HUGE_VAL;
      }
  long *kadj = (long *) mycalloc ((size_t) numpop, sizeof (long));
  jc_rbc.k = -1;
  if (ly.rb)
    {
      if (ly.sharp)
        jc_sharp_find (world, &ly, st);
      else
        jc_rb_find (world, &ly, st, world->options->inheritance_scalars[world->locus]);
    }
  if (getenv ("MIGRATE_JC_WINDOW") != NULL)
    {   /* diagnostics (S = 0 design): per split, the window in which all its
           crossings could slide together without changing the topology:
           above every crossing lineage's lower event, below every upper one */
      long k;
      for (k = 0; k < ly.nsplit; k++)
        {
          species_fmt *sp = &world->species_model[ly.split_model[k]];
          long n = 0;
          double amax = -1.0, bmin = HUGE_VAL, cmin = HUGE_VAL, cmax = -1.0;
          for (i = 1; i < T - 1; i++)
            if (tl[i].eventnode->type == 'd'
                && get_fixed_species_model (tl[i].eventnode->pop, tl[i].eventnode->actualpop,
                                            world->species_model, world->species_model_size) == sp)
              {
                node *E = tl[i].eventnode, *n1 = E->back, *n2 = E->next->back;
                if (n1 == NULL || n2 == NULL)
                  continue;
                const double t1 = jc_ntime (n1), t2 = jc_ntime (n2);
                n++;
                if ((t1 < t2 ? t1 : t2) > amax) amax = t1 < t2 ? t1 : t2;
                if ((t1 > t2 ? t1 : t2) < bmin) bmin = t1 > t2 ? t1 : t2;
                if (tl[i].age < cmin) cmin = tl[i].age;
                if (tl[i].age > cmax) cmax = tl[i].age;
              }
          fprintf (stderr, "JCWIN %ld %ld %ld %.6g %.6g %.6g %.6g\n", world->locus + 1, k, n, amax, bmin, cmin, cmax);
        }
    }
  for (i = 1; i < T - 1; i++)
    {
      const double t0 = tl[i - 1].age, t1 = tl[i].age;
      const double dt = t1 - t0;
      const long *k = kadj;
      const char type = tl[i].eventnode->type;
      for (pop = 0; pop < numpop; pop++)
        kadj[pop] = jc_rb_count (tl, i, pop);
      /* intervals ending at a tip count (dated tips): waiting terms only */
      if (ly.nseg > 1)
        {   /* skyline=PARAM: the waiting terms per segment */
          long sg = jc_seg (&ly, t0);
          double a = t0;
          while (a < t1)
            {
              const double b = (sg + 1 < ly.nseg && ly.segt[sg + 1] < t1) ? ly.segt[sg + 1] : t1;
              for (pop = 0; pop < numpop; pop++)
                {
                  st[ly.off_sK + sg * numpop + pop] += (double) k[pop] * (double) (k[pop] - 1) * (b - a);
                  st[ly.off_sS + sg * numpop + pop] += (double) k[pop] * (b - a);
                  if (ly.off_sA[pop] >= 0 && k[pop] > 1 && b > a)
                    {   /* growth restarted at the segment start */
                      const double *gg = ly.ggrid + ly.pop_slot[pop] * JC_NG;
                      const double lk = log ((double) k[pop] * (double) (k[pop] - 1));
                      double *As = st + ly.off_sA[pop] + sg * JC_NG;
                      for (g = 0; g < JC_NG; g++)
                        As[g] = jc_logadd (As[g], lk + jc_log_growth_dt (gg[g], a - ly.segt[sg], b - ly.segt[sg]));
                    }
                }
              a = b;
              sg++;
            }
        }
      for (pop = 0; pop < numpop; pop++)
        {
          const double kk = (double) k[pop] * (double) (k[pop] - 1);
          st[ly.off_K + pop] += kk * dt;
          st[ly.off_S + pop] += (double) k[pop] * dt;
          if (ly.pop_slot[pop] >= 0 && kk > 0.0)
            {
              const double *gg = ly.ggrid + ly.pop_slot[pop] * JC_NG;
              for (g = 0; g < JC_NG; g++)
                st[ly.off_A[pop] + g] = jc_logadd (st[ly.off_A[pop] + g],
                                                   log (kk) + jc_log_growth_dt (gg[g], t0, t1));
            }
        }
      if (type == 'i')
        {
          const long xp = tl[i].eventnode->actualpop;
          st[ly.off_c + xp] += 1.0;
          if (ly.nseg > 1)
            st[ly.off_sc + jc_seg (&ly, t1) * numpop + xp] += 1.0;
          if (ly.pop_slot[xp] >= 0)   /* restarted growth: the age within its segment */
            st[ly.off_t[xp]] += ly.off_sA[xp] >= 0 ? t1 - ly.segt[jc_seg (&ly, t1)] : t1;
        }
      else if (type == 'm')
        {
          const long j = m2mmm (tl[i].eventnode->pop, tl[i].eventnode->actualpop, numpop);
          if (j >= numpop && j < world->numpop2 && ly.mslot[j] >= 0)
            {
              st[ly.off_m + ly.mslot[j]] += 1.0;
              if (ly.nseg > 1)
                st[ly.off_sm + jc_seg (&ly, t1) * ly.nmig + ly.mslot[j]] += 1.0;
            }
        }
    }
  /* gated migration entries: last event age and S(D) on the grid */
  for (pop = ly.numpop; pop < ly.numpop2; pop++)
    if (ly.gate_k[pop] >= 0)
      {
        long from, to;
        double last = 0.0;
        double *S = st + ly.off_G[pop] + 1;
        m2mm (pop, numpop, &from, &to);
        for (i = 1; i < T - 1; i++)
          if (tl[i].eventnode->type == 'm'
              && m2mmm (tl[i].eventnode->pop, tl[i].eventnode->actualpop, numpop) == pop
              && tl[i].age > last)
            last = tl[i].age;
        st[ly.off_G[pop]] = last;
        for (g = 0; g < JC_NGATE; g++)
          {
            const double x = jc_gate_x (&ly, pop, g);
            double v = 0.0;
            for (i = 1; i < T - 1; i++)
              {
                const double t0 = tl[i - 1].age, t1 = tl[i].age;
                if (t0 >= x)
                  break;
                v += (double) jc_rb_count (tl, i, to) * ((t1 < x ? t1 : x) - t0);
              }
            S[g] = v;
          }
      }
  /* the locus coalesces at inheritance scalar x the reference Theta: scale
     the coalescence integrals by 1/scalar; the remaining factor scalar^-c
     is the same for every Theta and cancels in p(G|phi)/Z(G) */
  /* the locus rate with which this genealogy was sampled: an estimated rate
     modifier is a nuisance parameter of the locus, integrated out by the
     mixture over its genealogies */
  st[ly.off_mu] = world->options->mu_rates[world->locus];
  st[ly.off_rep] = (double) world->rep;
  const double inh = world->options->inheritance_scalars[world->locus];
  if (ly.ml)
    {   /* the spells, sorted by population */
      nsp = mlh_spells (world, &spells, &spalloc);
      for (pop = 0; pop < numpop; pop++)
        st[ly.off_nsp + pop] = 0.0;
      for (i = 0; i < nsp; i++)
        st[ly.off_nsp + JML_POP (spells[2 * i], numpop)] += 1.0;
      st[ly.off_inh] = inh;
      st[ly.off_lc] = log (mlh_locus_sites (world, world->locus));
      if (getenv ("MIGRATE_JC_MLCHECK") != NULL)
        {   /* the spells reproduce the sampler's p(G|params) */
          static double worst = 0.0;
          const double a = mlh_spells_logp (world, spells, nsp), b = mlh_probg_treetimes (world);
          if (fabs (a - b) > worst)
            {
              worst = fabs (a - b);
              fprintf (stderr, "JC_MLCHECK spells vs probg: worst |diff| %g (logp %g, %ld spells)\n", worst, b, nsp);
            }
          {   /* the joint's own term (its parameter mapping and function
                 table) at the sampler's parameters */
            static double worst2 = 0.0;
            double *ps = (double *) mycalloc ((size_t) (2 * (nsp > 0 ? nsp : 1)), sizeof (double));
            double *phi = (double *) mycalloc ((size_t) ly.nphi, sizeof (double));
            double amin = 1.0, amax = 0.0, v = 0.0;
            long q, off = 0;
            memcpy (phi, world->param0, sizeof (double) * (size_t) ly.np);
            for (pop = 0; pop < numpop; pop++)
              if (ly.apar[pop] >= 0)
                {
                  phi[ly.apar[pop]] = mlh_alpha (world, pop);
                  if (world->bayes->minparam[ly.apar[pop]] < amin) amin = world->bayes->minparam[ly.apar[pop]];
                  if (world->bayes->maxparam[ly.apar[pop]] > amax) amax = world->bayes->maxparam[ly.apar[pop]];
                }
              else if (ly.afix[pop] < 1.0)
                {
                  if (ly.afix[pop] < amin) amin = ly.afix[pop];
                  if (ly.afix[pop] > amax) amax = ly.afix[pop];
                }
            if (amin < 1.0)
              jml_table_init (amin, amax);
            for (pop = 0; pop < numpop; pop++)
              {   /* prepare each population's spells (sorted, log dt) */
                const long n = (long) st[ly.off_nsp + pop];
                memcpy (ps + 2 * off, spells + 2 * off, sizeof (double) * (size_t) (2 * n));
                jml_spells_prepare (ps + 2 * off, n, numpop);
                off += n;
              }
            for (q = 0; q < numpop; q++)
              v += jc_ml_term (&ly, st, ps, q, phi);
            if (fabs (v - b) > worst2)
              {
                worst2 = fabs (v - b);
                fprintf (stderr, "JC_MLCHECK joint term vs probg: worst |diff| %g (logp %g)\n", worst2, b);
              }
            myfree (ps);
            myfree (phi);
          }
        }
    }
  if (inh != 1.0)
    for (pop = 0; pop < numpop; pop++)
      {
        st[ly.off_K + pop] /= inh;
        for (g = 0; g < ly.nseg && ly.nseg > 1; g++)
          st[ly.off_sK + g * numpop + pop] /= inh;
        if (ly.pop_slot[pop] >= 0)
          for (g = 0; g < JC_NG; g++)
            st[ly.off_A[pop] + g] -= log (inh);
        if (ly.off_sA[pop] >= 0)
          for (g = 0; g < ly.nseg * JC_NG; g++)
            st[ly.off_sA[pop] + g] -= log (inh);
      }
  if (jc_file_mode (world))
    {
      /* thin while writing: about 2 * jc_maxsamples() rows per file */
      static long jc_file_count = 0;
      long fstride = world->options->lsteps / (2 * jc_maxsamples (world, ly.nrow));
      if (fstride < 1)
        fstride = 1;
      if (jc_file == NULL || jc_file_locus != world->locus || jc_file_rep != world->rep)
        jc_file_count = 0;
      if ((jc_file_count++ % fstride) != 0)
        {
          myfree (st);
          myfree (kadj);
          jc_rbc.k = -1;
          jc_layout_free (&ly);
          return;
        }
      jc_record_divergence (world, &ly, st);
      jc_sharp_clear (&ly, st);
      jc_record_rate (world, &ly, st);   /* last: it rescales the tree (tl) */
      if (jc_file == NULL || jc_file_locus != world->locus || jc_file_rep != world->rep)
        {
          char name[LINESIZE];
          if (jc_file != NULL)
            fclose (jc_file);
          jc_filename (world, world->locus, world->rep, name);
          /* recover=YES continues a replicate that a crash interrupted:
             keep the rows written before (a fresh run starts the file) */
          jc_file = fopen (name, jc_recover ? "a" : "w");
          jc_file_locus = world->locus;
          jc_file_rep = world->rep;
        }
      if (jc_file != NULL)
        {
          for (i = 0; i < ly.nrow; i++)
            fprintf (jc_file, i ? " %.17g" : "%.17g", st[i]);
          for (i = 0; i < 2 * nsp; i++)   /* Mittag-Leffler: the spells */
            fprintf (jc_file, " %.17g", spells[i]);
          fprintf (jc_file, "\n");
          fflush (jc_file);
        }
    }
  else
    {
      jc_store *js = jc_get_store (world, ly.nrow);
      jc_buffer *b = &js->b[world->locus];
      if (++b->seen >= b->stride)
        {
          b->seen = 0;
          jc_record_divergence (world, &ly, st);   /* only for kept rows */
          jc_sharp_clear (&ly, st);
          if (getenv ("MIGRATE_JC_RBCHECK") != NULL && ly.rb && st[ly.off_R] >= 0.0)
            jc_rb_check (world, &ly, st);
          if (getenv ("MIGRATE_JC_SHARPCHECK") != NULL)
            jc_sharp_check (world, &ly, st);
          if (getenv ("MIGRATE_JC_GATECHECK") != NULL && !ly.ml)
            {   /* the joint's term at the sampler's parameters against
                   probg_treetimes(): the difference must not depend on G */
              static double dmin = HUGE_VAL, dmax = -HUGE_VAL;
              boolean *act = (boolean *) mycalloc ((size_t) ly.nphi, sizeof (boolean));
              long q;
              for (q = 0; q < ly.nphi; q++)
                act[q] = TRUE;
              if (ly.rate)
                act[ly.rate_p] = FALSE;   /* not part of p(G|params) */
              /* phi: the parameters, then the skyline multipliers (tau) */
              double *phi = (double *) mycalloc ((size_t) ly.nphi, sizeof (double));
              long jj, sg2;
              memcpy (phi, world->param0, sizeof (double) * (size_t) ly.np);
              for (jj = 0; jj < ly.nfree; jj++)
                for (sg2 = 1; sg2 < ly.nseg; sg2++)
                  phi[JC_TAU (&ly, jj, sg2)] = world->timek[sg2 * world->numpop2 + ly.sky_k[jj]];
              const double d = jc_full_term (&ly, act, st, NULL, phi, st[ly.off_mu])
                - probg_treetimes (world);
              myfree (phi);
              if (d < dmin) dmin = d;
              if (d > dmax) dmax = d;
              fprintf (stderr, "JC_GATECHECK locus %ld: joint - probg = %.6f (range %.6f)\n",
                       world->locus + 1, d, dmax - dmin);
              for (q = 0; q < ly.nsplit; q++)
                {
                  const double mu0 = world->param0[ly.split_pmu[q]];
                  const double sg0 = ly.split_psig[q] >= 0 ? world->param0[ly.split_psig[q]] : ly.split_sfix[q];
                  fprintf (stderr, "JC_DIVCHECK locus %ld split %ld grid %g exact %g\n", world->locus + 1, q,
                           jc_div_term (&ly, st, q, mu0, sg0), jc_div_exact (world, ly.split_model[q], mu0, sg0));
                }
              for (q = ly.numpop; q < ly.numpop2; q++)
                if (ly.gate_k[q] >= 0)
                  {   /* the interpolated S(D) against the exact one */
                    long fr, to2, ii;
                    const double D = world->param0[ly.split_pmu[ly.gate_k[q]]];
                    double ex = 0.0;
                    m2mm (q, numpop, &fr, &to2);
                    for (ii = 1; ii < T - 1; ii++)
                      {
                        const double a0 = tl[ii - 1].age, a1 = tl[ii].age;
                        if (a0 >= D)
                          break;
                        ex += (double) tl[ii].lineages[to2] * ((a1 < D ? a1 : D) - a0);
                      }
                    fprintf (stderr, "JC_GATESCHECK locus %ld entry %ld D %g last %g S grid %g exact %g M %g\n",
                             world->locus + 1, q, D, st[ly.off_G[q]], jc_gate_S (&ly, st, q, D), ex, world->param0[q]);
                  }
              myfree (act);
            }
          jc_record_rate (world, &ly, st);   /* last: it rescales the tree (tl) */
          jc_append (js, world->locus, st, spells, nsp);
          if (b->n >= 2 * jc_maxsamples (world, js->nrow))
            {
              /* thin: keep every other row, record half as often */
              long src = 0, dst = 0;
              for (i = 0; i < b->n; i++)
                {
                  const long k = b->nspr[i];
                  if (i % 2 == 0)
                    {
                      if (k > 0 && dst != src)
                        memmove (b->sp + 2 * dst, b->sp + 2 * src, sizeof (double) * (size_t) (2 * k));
                      b->nspr[i / 2] = k;
                      dst += k;
                    }
                  src += k;
                }
              b->nsp = dst;
              for (i = 0; i < b->n / 2; i++)
                memmove (b->rows + i * js->nrow, b->rows + 2 * i * js->nrow,
                         sizeof (double) * (size_t) js->nrow);
              b->n /= 2;
              b->stride *= 2;
            }
        }
    }
  myfree (st);
  myfree (kadj);
  jc_rbc.k = -1;
  jc_layout_free (&ly);
}

/// drop this rank's statistics for one locus (a replicate worker, after
/// sending them to the locus worker)
void
jc_clear_locus (world_fmt *world, long locus)
{
  jc_store *js = (jc_store *) world->jointstats;
  if (js == NULL || locus < 0 || locus >= js->loci)
    return;
  myfree (js->b[locus].rows);
  js->b[locus].rows = NULL;
  myfree (js->b[locus].nspr);
  js->b[locus].nspr = NULL;
  myfree (js->b[locus].sp);
  js->b[locus].sp = NULL;
  js->b[locus].nsp = js->b[locus].capsp = 0;
  myfree (js->b[locus].zc);
  js->b[locus].zc = NULL;
  js->b[locus].nzc = 0;
  js->b[locus].n = js->b[locus].cap = js->b[locus].seen = 0;
  js->b[locus].stride = 1;
}

/// MPI: pack this rank's statistics for one locus: {nrow, n, rows...,
/// total spells, spells per row..., spells...}
long
jc_pack_buffer (MYREAL **buffer, world_fmt *world, long locus, long maxrep, long numpop)
{
  (void) maxrep;
  (void) numpop;
  jc_store *js = (jc_store *) world->jointstats;
  const long n = (js != NULL) ? js->b[locus].n : 0;
  const long nr = (js != NULL) ? js->nrow : 0;
  const long tot = (js != NULL) ? js->b[locus].nsp : 0;
  const long size = 2 + n * nr + 1 + n + 2 * tot;
  long i;
  *buffer = (MYREAL *) myrealloc (*buffer, sizeof (MYREAL) * (size_t) (size + 1));
  (*buffer)[0] = (MYREAL) nr;
  (*buffer)[1] = (MYREAL) n;
  if (n > 0)
    memcpy (*buffer + 2, js->b[locus].rows, sizeof (MYREAL) * (size_t) (n * nr));
  (*buffer)[2 + n * nr] = (MYREAL) tot;
  for (i = 0; i < n; i++)
    (*buffer)[3 + n * nr + i] = (MYREAL) js->b[locus].nspr[i];
  if (tot > 0)
    memcpy (*buffer + 3 + n * nr + n, js->b[locus].sp, sizeof (MYREAL) * (size_t) (2 * tot));
  return size;
}

/// MPI: add one locus' statistics received from another rank
void
jc_unpack_buffer (MYREAL *buffer, world_fmt *world, long locus, long maxrep, long numpop)
{
  (void) maxrep;
  (void) numpop;
  const long nr = (long) buffer[0];
  const long n = (long) buffer[1];
  long i;
  if (nr <= 0 || n <= 0 || locus < 0 || locus >= world->loci)
    return;
  jc_store *js = jc_get_store (world, nr);
  const MYREAL *cnt = buffer + 3 + n * nr;
  const MYREAL *sp = cnt + n;
  long off = 0;
  for (i = 0; i < n; i++)
    {
      const long k = (long) cnt[i];
      jc_append (js, locus, buffer + 2 + i * nr, sp + 2 * off, k);
      off += k;
    }
}

/* ------------------------------------------------------------- combine */

typedef struct
{
  long n;          /* samples kept */
  long total;      /* samples stored before thinning */
  double mu;       /* locus rate */
  double *st;      /* n x nrow */
  double *cur;     /* n: current log p(G|phi) - log Z per sample */
  double *logz;    /* n: log Z per sample */
  double *unit;    /* n x nunit: the terms of jc_full_term() per sample */
  long locus;      /* the locus number (seeds its bootstrap) */
  double *sp;      /* Mittag-Leffler: the spells of the kept rows */
  long *spo;       /* n: offset of each row's spells in sp */
  long *src;       /* n: index of each row among all stored rows (log Z cache) */
  long *ax, na;    /* the active rows (skipping negligible genealogies) */
  double sin, m0;  /* the inactive rows: sum exp(cur - m0) at the last refresh */
  double *h;       /* per split: bandwidth of the split-time smoothing (0: none) */
} jc_locus;

/* Smoothing of the split-time terms (MIGRATE_JC_SMOOTH=FLOOR or KERNEL; off
   by default, 2026-11-16). With a split-time std S near zero every stored
   genealogy pins the split time: its term in D is a spike of width S around
   its split events, a locus' mixture over its genealogies is a comb of
   spikes, and the joint sampler finds a D where one spike per locus lines up
   (one genealogy per locus carries the weight; d80 of cookbook 11a). The
   spikes are an artifact of a finite sample of genealogies. KERNEL convolves
   each genealogy's split-event density with a Gaussian kernel of width h in
   the mean: for the normal split-event terms that is the std
   sqrt(S^2 + h^2); FLOOR uses max(S, h). The same std enters p(G|phi) and
   its normalizer Z(G), so every genealogy still contributes a normalised
   density (a kernel density estimate over the genealogies). h per locus and
   split: Silverman's rule on the split-event centres of the locus' stored
   genealogies (times MIGRATE_JC_HSCALE), from all of its rows, so that the
   blocks and the cached Z(G) use the same h. */
static int jc_smooth = -1;          /* 0 off, 1 FLOOR, 2 KERNEL */
static double jc_hscale = 1.0;
static const double *jc_hcur = NULL; /* h of the locus being evaluated */

static int
jc_smooth_mode (void)
{
  if (jc_smooth < 0)
    {
      const char *e = getenv ("MIGRATE_JC_SMOOTH");
      jc_smooth = 0;
      if (e != NULL && (e[0] == 'F' || e[0] == 'f' || e[0] == '1'))
        jc_smooth = 1;
      if (e != NULL && (e[0] == 'K' || e[0] == 'k' || e[0] == '2'))
        jc_smooth = 2;
      if (getenv ("MIGRATE_JC_HSCALE") != NULL && atof (getenv ("MIGRATE_JC_HSCALE")) > 0.0)
        jc_hscale = atof (getenv ("MIGRATE_JC_HSCALE"));
    }
  return jc_smooth;
}

static int jc_cmp_double (const void *a, const void *b);

/* the bandwidths h[k] of a locus from m rows (rows + idx[i] * nrow, or
   consecutive rows when idx is NULL) */
static void
jc_locus_bandwidth (const jc_layout *ly, const double *rows, const long *idx, long m, double *h)
{
  long k, i;
  for (k = 0; k < ly->nsplit; k++)
    h[k] = 0.0;
  if (jc_smooth_mode () == 0 || m < 2)
    return;
  double *c = (double *) mycalloc ((size_t) m, sizeof (double));
  for (k = 0; k < ly->nsplit; k++)
    {
      long nc = 0;
      double mean = 0.0, var = 0.0;
      for (i = 0; i < m; i++)
        {
          const double *st = rows + (idx != NULL ? idx[i] : i) * ly->nrow;
          if (st[ly->off_D[k] + 1] > 0.0)
            c[nc++] = st[ly->off_D[k]];
        }
      if (nc < 2)
        continue;
      for (i = 0; i < nc; i++)
        mean += c[i] / nc;
      for (i = 0; i < nc; i++)
        var += (c[i] - mean) * (c[i] - mean) / (nc - 1);
      qsort (c, (size_t) nc, sizeof (double), jc_cmp_double);
      const double iqr = (c[(3 * nc) / 4] - c[nc / 4]) / 1.34;
      double s = sqrt (var);
      if (iqr > 0.0 && iqr < s)
        s = iqr;
      h[k] = jc_hscale * 0.9 * s * pow ((double) nc, -0.2);
      if (getenv ("MIGRATE_JC_DIAG") != NULL)
        fprintf (stderr, "JCBANDWIDTH split %ld: %ld genealogies with split events of %ld, h %.3g\n",
                 k, nc, m, h[k]);
    }
  myfree (c);
}

static double
jc_logsumexp (const double *v, long n)
{
  double mx = -HUGE_VAL, s = 0.0;
  long i;
  for (i = 0; i < n; i++)
    if (v[i] > mx)
      mx = v[i];
  if (mx == -HUGE_VAL)
    return mx;
  for (i = 0; i < n; i++)
    s += exp (v[i] - mx);
  return mx + log (s);
}

/* log A_i(g) for a growing population: quadratic interpolation of the
   grid of log A (A itself changes by orders of magnitude between nodes) */
static double
jc_logA_grid (const jc_layout *ly, const double *A, long pop, double g)
{
  const double *gg = ly->ggrid + ly->pop_slot[pop] * JC_NG;
  if (A[0] == -HUGE_VAL)   /* no interval with two lineages */
    return -HUGE_VAL;
  const double h = gg[1] - gg[0];
  double x = (g - gg[0]) / h;
  long j = (long) floor (x + 0.5);
  if (j < 1)
    j = 1;
  if (j > JC_NG - 2)
    j = JC_NG - 2;
  x -= (double) j;
  return A[j] + 0.5 * x * (A[j + 1] - A[j - 1]) + 0.5 * x * x * (A[j + 1] - 2.0 * A[j] + A[j - 1]);
}

static double
jc_logA (const jc_layout *ly, const double *st, long pop, double g)
{
  return jc_logA_grid (ly, st + ly->off_A[pop], pop, g);
}

/* Z(G) of a growth group with skyline multipliers: the growth node the
   segment terms use (-1: interpolate at the given g) */
static long jc_zgidx = -1;

/* log of the last argument: within a proposal all genealogies of a locus
   share the argument (the same number as calling log() each time) */
static double
jc_log_memo (double x, double *lastx, double *lastlog)
{
  if (x != *lastx)
    {
      *lastx = x;
      *lastlog = log (x);
    }
  return *lastlog;
}

/* log term of population pop with size theta and growth g (ignored when the
   population does not grow); logA is log A_i(g), NAN means: look it up; dA
   is added to the coalescence integral A (sharp splits: the crossings moved
   to D; the moved integral is never negative) */
static double
jc_pop_term_dA (const jc_layout *ly, const double *st, long pop, double theta, double g,
                double mu, double logA, double dA)
{
  static double lx = -1.0, ll = 0.0;
  if (theta <= 0.0)
    return -HUGE_VAL;
  const double c = st[ly->off_c + pop];
  double A = st[ly->off_K + pop], tg = 0.0;
  if (ly->pop_slot[pop] >= 0)
    {
      A = exp (isnan (logA) ? jc_logA (ly, st, pop, g) : logA);   /* inf: zero density */
      tg = g * st[ly->off_t[pop]];
    }
  if (dA != 0.0)
    {
      A += dA;
      if (!(A > 0.0))
        A = 0.0;
    }
  return c * (LOG2 - jc_log_memo (mu * theta, &lx, &ll)) - A / (mu * theta) + tg;
}

static double
jc_pop_term (const jc_layout *ly, const double *st, long pop, double theta, double g,
             double mu, double logA)
{
  return jc_pop_term_dA (ly, st, pop, theta, g, mu, logA, 0.0);
}

/* the population term of the joint (full term and units): with sharp
   splits a daughter or ancestor (sharp_pu) is evaluated on the row with
   the crossings at D, in one piece; as base term plus a separate change
   (unit JC_U_SPLITRB) the two cancel to round-off garbage where the
   growth-weighted integrals are huge (stepping-stone rungs near the prior,
   2026-11-27) */
static double
jc_pop_term_phi (const jc_layout *ly, const double *st, long pop, const double *phi, double mu)
{
  const long s = ly->pop_slot[pop];
  const double g = s >= 0 ? phi[ly->slot_param[s]] : 0.0;
  double dA = 0.0;
  if (ly->sharp && ly->sharp_pu[pop] && st[ly->off_R] == 0.0)
    {
      double c[64], gp[64], dK[64], dS[64];
      long k, p;
      for (k = 0; k < ly->nsplit; k++)
        c[k] = phi[ly->split_pmu[k]];
      for (p = 0; p < ly->numpop; p++)
        gp[p] = ly->pop_slot[p] >= 0 ? phi[ly->slot_param[ly->pop_slot[p]]] : 0.0;
      jc_sharp_delta (ly, st + ly->off_R, c, gp, dK, dS);
      dA = dK[pop];
    }
  return jc_pop_term_dA (ly, st, pop, phi[ly->rep[pop]], g, mu, NAN, dA);
}

/* migration term of matrix entry e at rate x (the rate, M) */
static double
jc_mig_term (const jc_layout *ly, const double *st, long e, double x, double mu)
{
  static double lx = -1.0, ll = 0.0;
  long from, to;
  if (x <= 0.0)
    return -HUGE_VAL;
  m2mm (e, ly->numpop, &from, &to);
  const double m = ly->mslot[e] >= 0 ? st[ly->off_m + ly->mslot[e]] : 0.0;
  return (m > 0.0 ? m * jc_log_memo (x, &lx, &ll) : 0.0) - x * st[ly->off_S + to] / mu;
}

/* skyline=PARAM: population pop in segment sg with size v (= Theta tau) */
static double
jc_pop_seg_term (const jc_layout *ly, const double *st, long pop, long sg, double v, double g, double mu)
{
  static double lx = -1.0, ll = 0.0;
  if (v <= 0.0)
    return -HUGE_VAL;
  const double c = st[ly->off_sc + sg * ly->numpop + pop];
  if (ly->off_sA[pop] >= 0)
    {   /* growth restarted at the segment start; segment 0 carries g sum (t - t_s) */
      const double *As = st + ly->off_sA[pop] + sg * JC_NG;
      double lA;
      if (jc_zgidx >= 0)
        {
          lA = As[jc_zgidx];
          g = ly->ggrid[ly->pop_slot[pop] * JC_NG + jc_zgidx];
        }
      else
        lA = jc_logA_grid (ly, As, pop, g);
      return c * (LOG2 - jc_log_memo (mu * v, &lx, &ll)) - exp (lA) / (mu * v)
        + (sg == 0 ? g * st[ly->off_t[pop]] : 0.0);
    }
  return c * (LOG2 - jc_log_memo (mu * v, &lx, &ll)) - st[ly->off_sK + sg * ly->numpop + pop] / (mu * v);
}

/* skyline=PARAM: migration entry e in segment sg at rate x (= M tau) */
static double
jc_mig_seg_term (const jc_layout *ly, const double *st, long e, long sg, double x, double mu)
{
  static double lx = -1.0, ll = 0.0;
  long from, to;
  if (x <= 0.0)
    return -HUGE_VAL;
  m2mm (e, ly->numpop, &from, &to);
  const double m = ly->mslot[e] >= 0 ? st[ly->off_sm + sg * ly->nmig + ly->mslot[e]] : 0.0;
  return (m > 0.0 ? m * jc_log_memo (x, &lx, &ll) : 0.0) - x * st[ly->off_sS + sg * ly->numpop + to] / mu;
}

/* term of a free skyline entry e (Theta or M) in segment sg at value x */
static double
jc_sky_seg_term (const jc_layout *ly, const double *st, long e, long sg, double x, double mu)
{
  return e < ly->numpop ? jc_pop_seg_term (ly, st, e, sg, x, 0.0, mu) : jc_mig_seg_term (ly, st, e, sg, x, mu);
}

/* log prior of the multipliers (on log tau, as bayes_update_timeparam()),
   -HUGE_VAL outside: |log tau| <= SKYLINE_HISTMAX and the log-uniform
   bounds, as in the sampler */
static double
jc_sky_logprior (const jc_layout *ly, const double *phi)
{
  double lp = 0.0;
  long j, sg;
  for (j = 0; j < ly->nfree; j++)
    {
      double prev = 0.0;
      for (sg = 1; sg < ly->nseg; sg++)
        {
          const double tau = phi[JC_TAU (ly, j, sg)];
          const double u = log (tau);
          if (!(tau > 0.0) || fabs (u) > SKYLINE_HISTMAX)
            return -HUGE_VAL;
          if (ly->sky_type[j] == SKYPRIOR_LOGUNIFORM)
            {
              if (u < log (ly->sky_a[j]) || u > log (ly->sky_b[j]))
                return -HUGE_VAL;
            }
          else
            lp += -0.5 * (u - prev) * (u - prev) / (ly->sky_a[j] * ly->sky_a[j]);
          prev = u;
        }
    }
  return lp;
}

/* D of split k on std row b at mean mu: linear on the row's local grid; zero
   density (JC_DOUT) outside the window, where the true term is at least ~30
   log units below the peak (any extrapolation upward made the target
   unbounded) */
static double
jc_row_value (const jc_layout *ly, const double *st, long k, long b, double mu)
{
  const double *D = st + ly->off_D[k] + 2 + b * JC_NM;
  double L, U;
  jc_row_window (ly, st, k, b, &L, &U);
  if (mu < L)
    return JC_DOUT;
  if (mu > U)   /* without a split event the term is flat above its window */
    return (st[ly->off_D[k] + 1] == 0.0) ? D[JC_NM - 1] : JC_DOUT;
  if (!(U > L))   /* a fixed D: every node is at D */
    return D[0];
  /* quadratic through the nearest node and its neighbours (the term is
     close to quadratic in mu, exactly so for the split-event densities),
     clamped to the nodes' range so that it cannot overshoot */
  const double h = (U - L) / (JC_NM - 1);
  double x = (mu - L) / h;
  long j = (long) floor (x + 0.5);
  if (j < 1)
    j = 1;
  if (j > JC_NM - 2)
    j = JC_NM - 2;
  x -= (double) j;
  double v = D[j] + 0.5 * x * (D[j + 1] - D[j - 1]) + 0.5 * x * x * (D[j + 1] - 2.0 * D[j] + D[j - 1]);
  double lo3 = D[j - 1], hi3 = D[j - 1];
  if (D[j] < lo3) lo3 = D[j];
  if (D[j] > hi3) hi3 = D[j];
  if (D[j + 1] < lo3) lo3 = D[j + 1];
  if (D[j + 1] > hi3) hi3 = D[j + 1];
  return v < lo3 ? lo3 : (v > hi3 ? hi3 : v);
}

/* D(mu, sigma) of split k between std rows: the -n_d log(sigma) of the n_d
   normal split-event densities is taken out exactly, the rest is
   interpolated linearly in 1/sigma^2 between the two neighbouring rows
   (convex weights; exact for the normal split-event terms). Outside the std
   grid the nearest row is used. */
static double
jc_div_term (const jc_layout *ly, const double *st, long k, double mu, double sigma)
{
  const long ns = ly->split_ns[k];
  if (ns == 1)
    return jc_row_value (ly, st, k, 0, mu);
  if (jc_hcur != NULL && jc_hcur[k] > 0.0)   /* smoothing (see jc_locus) */
    sigma = jc_smooth == 1 ? (sigma > jc_hcur[k] ? sigma : jc_hcur[k])
      : sqrt (sigma * sigma + jc_hcur[k] * jc_hcur[k]);
  const double *sg = ly->sgrid + k * JC_NS;
  if (sigma <= sg[0])
    return jc_row_value (ly, st, k, 0, mu);
  if (sigma >= sg[ns - 1])
    return jc_row_value (ly, st, k, ns - 1, mu);
  long b = (long) floor ((log (sigma) - log (sg[0])) / (log (sg[1]) - log (sg[0])));
  if (b < 0)
    b = 0;
  if (b > ns - 2)
    b = ns - 2;
  const double nd = st[ly->off_D[k] + 1];
  const double u = 1.0 / (sigma * sigma), u0 = 1.0 / (sg[b] * sg[b]), u1 = 1.0 / (sg[b + 1] * sg[b + 1]);
  const double w = (u - u1) / (u0 - u1);
  const double d0 = jc_row_value (ly, st, k, b, mu), d1 = jc_row_value (ly, st, k, b + 1, mu);
  if (d0 <= JC_DOUT || d1 <= JC_DOUT)
    return JC_DOUT;
  return w * (d0 + nd * log (sg[b])) + (1.0 - w) * (d1 + nd * log (sg[b + 1])) - nd * log (sigma);
}

static double
jc_value (const jc_layout *ly, const double *phi, long e)
{
  return ly->rep[e] >= 0 ? phi[ly->rep[e]] : ly->fixed[e];
}

/* migration rate of entry e: M, or xNm / Theta_to */
static double
jc_rate (const jc_layout *ly, const double *phi, long e)
{
  long from, to;
  if (!ly->xnm[e])
    return jc_value (ly, phi, e);
  m2mm (e, ly->numpop, &from, &to);
  return jc_value (ly, phi, e) / jc_value (ly, phi, to);
}

/* ------------------------------------- integrated crossing time (RB) */

/* log of the upper normal tail Q(z) = 1 - Phi(z), also far in the tail */
static double
jc_logQ (double z)
{
  if (z < 30.0)
    return log (0.5 * erfc (z / M_SQRT2));
  const double z2 = 1.0 / (z * z);
  return -0.5 * z * z - log (z) - 0.5 * log (2.0 * M_PI) + log1p (-z2 + 3.0 * z2 * z2 - 15.0 * z2 * z2 * z2);
}

/* log(Phi(z2) - Phi(z1)) for z1 < z2 */
static double
jc_logphidiff (double z1, double z2)
{
  if (!(z2 > z1))
    return -HUGE_VAL;
  if (z1 > 0.0)
    {
      const double q1 = jc_logQ (z1), q2 = jc_logQ (z2);
      return q1 + log1p (-exp (q2 - q1));
    }
  if (z2 < 0.0)
    {
      const double q1 = jc_logQ (-z2), q2 = jc_logQ (-z1);
      return q1 + log1p (-exp (q2 - q1));
    }
  return log1p (-(exp (jc_logQ (-z1)) + exp (jc_logQ (z2))));
}

/* slope of log p(G_c | phi) in c on a piece with slopes cd = 2 k_d / inh,
   ca = 2 k_a / inh of the other lineages, at c (gated migration counts
   below its split mean only). mode 0: every term; 2: only the terms that
   jc_full_term() leaves out (fixed or inactive parameters: they depend on
   c and enter numerator and normalizer alike) */
/* skyline=PARAM: the multiplier of matrix entry e in segment sg (1 without) */
static double
jc_rb_tau (const jc_layout *ly, const double *phi, long e, long sg)
{
  return (sg > 0 && ly->sky_j[e] >= 0) ? phi[JC_TAU (ly, ly->sky_j[e], sg)] : 1.0;
}

static double
jc_rb_slope (const jc_layout *ly, const boolean *active, const double *phi, long k,
             double cd, double ca, double c, double mu, int mode, double nmov)
{
  const long d = ly->split_to[k], an = ly->split_from[k];
  /* the multipliers are coordinates of phi, not of the parameter vector
     that mode 2 gets (free skyline entries have active parameters) */
  const long sg = mode == 2 ? 0 : jc_seg (ly, c);
  double lam = 0.0;
  long e, from, to, w;
  for (w = 0; w < 2; w++)
    {
      const long pop = w ? an : d;
      const boolean act = ly->rep[pop] >= 0 && active != NULL && active[ly->rep[pop]];
      if (mode == 2 && act)
        continue;
      const double th = jc_value (ly, phi, pop) * jc_rb_tau (ly, phi, pop, sg);
      if (!(th > 0.0))
        return -HUGE_VAL;
      lam += (w ? ca : -cd) / (mu * th);
    }
  for (e = ly->numpop; e < ly->numpop2; e++)
    {
      if (ly->mzero[e])
        continue;
      m2mm (e, ly->numpop, &from, &to);
      if (to != d && to != an)
        continue;
      const boolean act = ly->rep[e] >= 0 && active != NULL && active[ly->rep[e]];
      if (mode == 2 && act)
        continue;
      if (ly->gate_k[e] >= 0 && mode != 2 && !(c < phi[ly->split_pmu[ly->gate_k[e]]]))
        continue;
      const double x = jc_rate (ly, phi, e) * jc_rb_tau (ly, phi, e, sg);
      if (x > 0.0)
        lam += nmov * (to == d ? -x : x) / mu;
    }
  return lam;
}

/* log int_u^v f(c) e^{lam (c - u)} dc of the split density f of split k
   (normal: N(c; D, S), without the 1/S'(a) of the survival; exponential:
   e^{-(c - a)/D} / D) */
static double
jc_rb_piece (const jc_layout *ly, const double *phi, long k, double a, double u, double v, double lam,
             boolean normal)
{
  const double D = phi[ly->split_pmu[k]];
  if (!(D > 0.0) || !(v > u))
    return -HUGE_VAL;
  if (normal)
    {
      const double S = ly->split_psig[k] >= 0 ? phi[ly->split_psig[k]] : ly->split_sfix[k];
      if (!(S > 0.0))
        return -HUGE_VAL;
      const double m = D + lam * S * S;
      return lam * (D - u) + 0.5 * lam * lam * S * S + jc_logphidiff ((u - m) / S, (v - m) / S);
    }
  return -log (D) - (u - a) / D + jc_log_growth_dt (lam - 1.0 / D, 0.0, v - u);
}

/* the cut points of piece (u, v) where a gated entry into d or the
   ancestor stops (its split mean); returns their number, sorted */
static long
jc_rb_cuts (const jc_layout *ly, const double *phi, long k, double u, double v, double *cut)
{
  const long d = ly->split_to[k], an = ly->split_from[k];
  long e, from, to, n = 0, i, j;
  for (e = ly->numpop; e < ly->numpop2 && n < 8; e++)
    if (ly->gate_k[e] >= 0 && !ly->mzero[e])
      {
        m2mm (e, ly->numpop, &from, &to);
        const double g = phi[ly->split_pmu[ly->gate_k[e]]];
        if ((to == d || to == an) && g > u && g < v)
          cut[n++] = g;
      }
  for (i = 1; i < n; i++)
    for (j = i; j > 0 && cut[j - 1] > cut[j]; j--)
      {
        const double t = cut[j];
        cut[j] = cut[j - 1];
        cut[j - 1] = t;
      }
  return n;
}

/* the integrated crossing (unit JC_U_SPLITRB): log of
   int_a^b p(G_c | phi) / p(G_a | phi) f(c) dc over the crossing lineage's
   split terms (its waiting terms from a and its split event at c) */
static double
jc_rb_term (const jc_layout *ly, const double *st, const double *phi, double mu)
{
  if (!ly->rb)
    return 0.0;
  if (ly->sharp)
    return jc_sharp_term (ly, st, phi, mu);
  const double *R = st + ly->off_R;
  const long k = (long) R[0];
  if (k < 0)
    return 0.0;
  const double a = R[1], b = R[2];
  const long np = (long) R[3];
  const long d = ly->split_to[k], an = ly->split_from[k];
  double acc = -HUGE_VAL, L = 0.0, cut[10];
  long j, q;
  /* the slope's parts, once per call (jc_rb_slope() with mode 0) */
  const double thd = jc_value (ly, phi, d), tha = jc_value (ly, phi, an);
  if (!(thd > 0.0) || !(tha > 0.0))
    return JC_DOUT;
  const double wd = 1.0 / (mu * thd), wa = 1.0 / (mu * tha);
  double r0 = 0.0, gD[8], gv[8], sv[JC_RBMAXSKY];
  long ng = 0, ns = 0, se[JC_RBMAXSKY], e, from, to;
  for (e = ly->numpop; e < ly->numpop2; e++)
    {
      if (ly->mzero[e])
        continue;
      m2mm (e, ly->numpop, &from, &to);
      if (to != d && to != an)
        continue;
      const double x = jc_rate (ly, phi, e);
      if (!(x > 0.0))
        continue;
      const double r = (to == d ? -x : x) / mu;
      if (ly->gate_k[e] >= 0 && ng < 8)
        {
          gD[ng] = phi[ly->split_pmu[ly->gate_k[e]]];
          gv[ng++] = r;
        }
      else if (ly->nseg > 1 && ly->sky_j[e] >= 0)
        {   /* skyline=PARAM: times its multiplier in the piece's segment */
          se[ns] = e;
          sv[ns++] = r;
        }
      else
        r0 += r;
    }
  for (j = 0; j < np; j++)
    {
      const double u = R[4 + 3 * j], v = j + 1 < np ? R[4 + 3 * (j + 1)] : b;
      const long nc = jc_rb_cuts (ly, phi, k, u, v, cut);
      const long sg = ly->nseg > 1 ? jc_seg (ly, u) : 0;
      double lam0 = -R[4 + 3 * j + 1] * wd / jc_rb_tau (ly, phi, d, sg)
        + R[4 + 3 * j + 2] * wa / jc_rb_tau (ly, phi, an, sg) + r0;
      for (e = 0; e < ns; e++)
        lam0 += sv[e] * jc_rb_tau (ly, phi, se[e], sg);
      double s0 = u;
      for (q = 0; q <= nc; q++)
        {
          const double s1 = q < nc ? cut[q] : v;
          const double mid = 0.5 * (s0 + s1);
          double lam = lam0;
          for (e = 0; e < ng; e++)
            if (mid < gD[e])
              lam += gv[e];
          acc = jc_logadd (acc, L + jc_rb_piece (ly, phi, k, a, s0, s1, lam, jc_rb_normal));
          L += lam * (s1 - s0);
          s0 = s1;
        }
    }
  if (jc_rb_normal)
    {   /* the survival to a: S'(a) = Phi((b1 - D)/S) - Phi((a - D)/S) */
      const double D = phi[ly->split_pmu[k]];
      const double S = ly->split_psig[k] >= 0 ? phi[ly->split_psig[k]] : ly->split_sfix[k];
      acc -= jc_logphidiff ((a - D) / S, (ly->split_b1[k] - D) / S);
    }
  return acc > JC_DOUT ? acc : JC_DOUT;
}

/* log of the trapezoid integral of exp(v) over x */
static double
jc_logtrapz (const double *v, const double *x, long n)
{
  double *acc = (double *) mycalloc ((size_t) n, sizeof (double));
  long g;
  for (g = 0; g + 1 < n; g++)
    {
      const double a = v[g], b = v[g + 1];
      const double mx = a > b ? a : b;
      acc[g] = (mx == -HUGE_VAL) ? -HUGE_VAL
        : mx + log (0.5 * (exp (a - mx) + exp (b - mx)) * (x[g + 1] - x[g]));
    }
  double r = jc_logsumexp (acc, n - 1);
  myfree (acc);
  return r;
}

/* divergence-order TREE/RANK: TRUE when no group is older than its parent group
   (phi[p] read as x, when p >= 0) */
static boolean
jc_order_ok (const jc_layout *ly, const double *phi, long p, double x)
{
  long k;
  for (k = 0; k < ly->nsplit; k++)
    if (ly->split_par[k] >= 0)
      {
        const long pc = ly->split_pmu[k], pp = ly->split_pmu[ly->split_par[k]];
        const double vc = pc == p ? x : phi[pc], vp = pp == p ? x : phi[pp];
        if (vc > vp)
          return FALSE;
      }
  return TRUE;
}

/* divergence normalizer: integrates exactly the interpolant jc_div_term()
   that the sampler uses, over the whole prior range. In std: JC_DREF points
   per log cell (and the strip below the grid); in the mean, at each of
   them: a fine grid over the union of the two neighbouring rows' windows
   plus the rest of the prior range (where the term is tiny but not zero). */
#define JC_DREF 12
#define JC_DMU 257
typedef struct
{
  world_fmt *world;
  long ns;                      /* std points per split */
  double *sv, *sw, *ps_s;       /* stds, widths, log prior: nsplit x ns */
  double *cv;                   /* splitstd-cv: per split JC_CVN fields of p(S|D), NULL: off */
} jc_divquad;

/* splitstd-cv: log p(S|D) = log p_c(S/D) - log D. A closed-form kernel per
   c prior plus the constant of the exact prior (splitcv_logprior()); the
   kernel is checked against the exact prior at setup and replaced by it
   (kind -1) when it does not match. Fields: kind, min, max, a, scale b,
   mean, sd, constant K, S parameter */
#define JC_CVN 9
static double
jc_cv_kernel (const double *c, double x)
{
  switch ((int) c[0])
    {
    case UNIFORMPRIOR:
      return 0.0;
    case LOGUNIFORMPRIOR:
      return x > 0.0 ? -log (x) : -HUGE_VAL;
    case EXPPRIOR:
    case WEXPPRIOR:
    case GAMMAPRIOR:
    case WGAMMAPRIOR:
      if (x <= 0.0)
        return c[3] == 1.0 ? 0.0 : -HUGE_VAL;
      return (c[3] - 1.0) * log (x) - x / c[4];
    case NORMALPRIOR:
    case WNORMALPRIOR:
      {
        const double z = (x - c[5]) / c[6];
        return -0.5 * z * z;
      }
    default:
      return -HUGE_VAL;
    }
}

static double
jc_cv_lp (world_fmt *world, const double *c, double S, double D)
{
  double x, v;
  if (!(D > 0.0))
    return -HUGE_VAL;
  if (c[0] < 0.0)
    {
      v = splitcv_logprior (world, (long) c[8], S, D);
      return v > (double) -HUGE ? v : -HUGE_VAL;
    }
  x = S / D;
  if (x <= c[1] || x > c[2])
    return -HUGE_VAL;
  v = jc_cv_kernel (c, x);
  return v > -HUGE_VAL ? v + c[7] - log (D) : -HUGE_VAL;
}

static void
jc_cv_make (world_fmt *world, long ps, double *c)
{
  MYREAL f[7];
  long t;
  const int kind = splitcv_cprior (world, ps, f);
  c[0] = (double) kind;
  c[1] = f[0];
  c[2] = f[1];
  c[3] = f[4];
  c[4] = f[6];
  c[5] = f[3];
  c[6] = f[4];
  c[8] = (double) ps;
  {
    const double xr = c[1] + 0.37 * (c[2] - c[1]);
    const double k = jc_cv_kernel (c, xr);
    const MYREAL e = splitcv_logprior (world, ps, xr, 1.0);
    c[7] = (k > -HUGE_VAL && e > (double) -HUGE) ? e - k : 0.0;
    if (!(k > -HUGE_VAL && e > (double) -HUGE))
      c[0] = -1.0;
  }
  for (t = 1; t < 8 && c[0] >= 0.0; t++)
    {   /* the kernel must reproduce the exact prior at other points and D */
      const double x = c[1] + (c[2] - c[1]) * t / 8.0, D = 0.003 * t;
      const MYREAL e = splitcv_logprior (world, ps, x * D, D);
      const double k = jc_cv_lp (world, c, x * D, D);
      if (!(e > (double) -HUGE) || fabs (e - k) > 1e-6 * (1.0 + fabs (e)))
        c[0] = -1.0;
    }
}

static void
jc_divquad_make (world_fmt *world, const jc_layout *ly, jc_divquad *q)
{
  long k, g;
  memset (q, 0, sizeof (*q));
  q->world = world;
  q->ns = JC_DREF * JC_NS + 1;
  q->sv = (double *) mycalloc ((size_t) ((ly->nsplit + 1) * q->ns), sizeof (double));
  q->sw = (double *) mycalloc ((size_t) ((ly->nsplit + 1) * q->ns), sizeof (double));
  q->ps_s = (double *) mycalloc ((size_t) ((ly->nsplit + 1) * q->ns), sizeof (double));
  q->cv = NULL;
  for (k = 0; k < ly->nsplit; k++)
    {
      const long ps = ly->split_psig[k];
      if (ps < 0)
        continue;
      if (splitcv_is_s (world, ps, NULL))
        {   /* splitstd-cv: p(S|D) enters at every D (jc_div_mu_integral(),
               jc_div_order_grid()); the S weights carry no prior */
          if (q->cv == NULL)
            {
              long j;
              q->cv = (double *) mycalloc ((size_t) ((ly->nsplit + 1) * JC_CVN), sizeof (double));
              for (j = 0; j <= ly->nsplit; j++)
                q->cv[j * JC_CVN] = -2.0;   /* no c prior */
            }
          jc_cv_make (world, ps, q->cv + k * JC_CVN);
        }
      const double shi = world->bayes->maxparam[ps], smin = world->bayes->minparam[ps];
      const double slo = smin > 1e-3 * shi ? smin : 1e-3 * shi;
      q->sv[k * q->ns] = 0.5 * (smin + slo);
      q->sw[k * q->ns] = slo - smin;
      for (g = 1; g < q->ns; g++)
        {
          const double e0 = slo * pow (shi / slo, (double) (g - 1) / (q->ns - 1));
          const double e1 = slo * pow (shi / slo, (double) g / (q->ns - 1));
          q->sv[k * q->ns + g] = sqrt (e0 * e1);
          q->sw[k * q->ns + g] = e1 - e0;
        }
      for (g = 0; g < q->ns; g++)
        q->ps_s[k * q->ns + g] = q->sw[k * q->ns + g] > 0.0
          ? (q->cv != NULL && q->cv[k * JC_CVN] > -2.0 ? 0.0 : scaling_prior_joint (world, ps, q->sv[k * q->ns + g]))
          : -HUGE_VAL;
    }
}

/* splitstd-cv: log p(S|D) of split group k, 0 without it */
static double
jc_div_cv (const jc_divquad *q, long k, double sigma, double D)
{
  if (q->cv == NULL || q->cv[k * JC_CVN] <= -2.0)
    return 0.0;
  return jc_cv_lp (q->world, q->cv + k * JC_CVN, sigma, D);
}

static void
jc_divquad_free (jc_divquad *q)
{
  myfree (q->sv); myfree (q->sw); myfree (q->ps_s);
  if (q->cv != NULL)
    myfree (q->cv);
}

/* log int prior(mu) e^{D(mu, sigma)} dmu at one std, for split k */
static double
jc_div_mu_integral (const jc_divquad *q, const jc_layout *ly, const double *st, long k, double sigma)
{
  /* k is the first split of its group: tied splits (t/T) share the mean and
     std, so their terms are integrated together */
  const long pm = ly->split_pmu[k];
  const double lo = ly->mlo[k], hi = ly->mhi[k];
  double L = hi, U = lo;
  long m;
  if (!(hi > lo))
    {   /* a fixed D: the per-locus chains sampled at D, no integral and no prior */
      double v = (ly->ngate > 0 ? jc_gate_group (ly, st, k, lo) : 0.0) + jc_div_cv (q, k, sigma, lo);
      for (m = 0; m < ly->nsplit; m++)
        if (ly->split_rep[m] == k)
          v += jc_div_term (ly, st, m, lo, sigma);
      return v;
    }
  /* the union of the windows of all members' rows around sigma */
  for (m = 0; m < ly->nsplit; m++)
    {
      if (ly->split_rep[m] != k)
        continue;
      if (ly->split_ns[m] > 1)
        {
          const double *sg = ly->sgrid + m * JC_NS;
          double y = (log (sigma) - log (sg[0])) / (log (sg[1]) - log (sg[0]));
          long b = (long) floor (y < 0.0 ? 0.0 : y), j;
          if (b > ly->split_ns[m] - 2)
            b = ly->split_ns[m] - 2;
          for (j = b; j <= b + 1; j++)
            {
              double Lj, Uj;
              if (j < 0 || j >= ly->split_ns[m])
                continue;
              jc_row_window (ly, st, m, j, &Lj, &Uj);
              if (Lj < L)
                L = Lj;
              if (Uj > U)
                U = Uj;
            }
        }
      else
        {
          double Lj, Uj;
          jc_row_window (ly, st, m, 0, &Lj, &Uj);
          if (Lj < L)
            L = Lj;
          if (Uj > U)
            U = Uj;
        }
    }
  if (L > U)
    {
      L = lo;
      U = hi;
    }
  /* fine grid inside [L, U], coarse outside it */
  double x[JC_DMU + 64], v[JC_DMU + 64];
  long n = 0, g;
  const long nout = 31;
  if (L > lo)
    for (g = 0; g < nout; g++)
      x[n++] = lo + (L - lo) * (double) g / nout;
  for (g = 0; g < JC_DMU; g++)
    x[n++] = L + (U - L) * (double) g / (JC_DMU - 1);
  if (U < hi)
    for (g = 1; g <= nout; g++)
      x[n++] = U + (hi - U) * (double) g / nout;
  for (g = 0; g < n; g++)
    {
      /* inside the prior range despite rounding (see the Theta/M grids) */
      if (x[g] < lo)
        x[g] = lo;
      if (x[g] > hi)
        x[g] = hi;
      v[g] = scaling_prior_joint (q->world, pm, x[g]) + jc_div_cv (q, k, sigma, x[g]);
      if (ly->ngate > 0)
        v[g] += jc_gate_group (ly, st, k, x[g]);
      for (m = 0; m < ly->nsplit; m++)
        if (ly->split_rep[m] == k)
          v[g] += jc_div_term (ly, st, m, x[g], sigma);
    }
  return jc_logtrapz (v, x, n);
}

/* ordered splits (divergence-order TREE/RANK): the mean of group k on one grid
   for all stds: JC_DMU points over the union of the windows of every member
   and std row, JC_ORDP points over each row's own window (the rows span a
   factor 1000 in std, so the narrow peaks of the small stds need their own
   points), and the rest of the prior range (coarse); sorted. H[g] = log
   prior(mu) + log int prior(sigma) e^{D(mu, sigma)} dsigma */
#define JC_ORDP 65
static int jc_cmp_double (const void *a, const void *b);
#define JC_ORDOUT 31
static long
jc_div_order_cap (const jc_layout *ly)
{
  long k, m, most = 1;
  for (k = 0; k < ly->nsplit; k++)
    {
      long rows = 0;
      for (m = 0; m < ly->nsplit; m++)
        if (ly->split_rep[m] == k)
          rows += ly->split_ns[m];
      if (rows > most)
        most = rows;
    }
  return 2 * JC_ORDOUT + JC_DMU + most * JC_ORDP;
}

static long
jc_div_order_grid (const jc_divquad *q, const jc_layout *ly, const double *st, long k,
                   double *x, double *H)
{
  const long pm = ly->split_pmu[k];
  const double lo = ly->mlo[k], hi = ly->mhi[k];
  double L = hi, U = lo;
  long m, j, n = 0, g, b;
  const long nout = JC_ORDOUT;
  for (m = 0; m < ly->nsplit; m++)
    if (ly->split_rep[m] == k)
      for (j = 0; j < ly->split_ns[m]; j++)
        {
          double Lj, Uj;
          jc_row_window (ly, st, m, j, &Lj, &Uj);
          if (Lj < L)
            L = Lj;
          if (Uj > U)
            U = Uj;
          for (g = 0; g < JC_ORDP; g++)
            x[n++] = Lj + (Uj - Lj) * (double) g / (JC_ORDP - 1);
        }
  if (L > U)
    {
      L = lo;
      U = hi;
    }
  if (L > lo)
    for (g = 0; g < nout; g++)
      x[n++] = lo + (L - lo) * (double) g / nout;
  for (g = 0; g < JC_DMU; g++)
    x[n++] = L + (U - L) * (double) g / (JC_DMU - 1);
  if (U < hi)
    for (g = 1; g <= nout; g++)
      x[n++] = U + (hi - U) * (double) g / nout;
  for (g = 0; g < n; g++)
    {
      if (x[g] < lo)
        x[g] = lo;
      if (x[g] > hi)
        x[g] = hi;
    }
  if (!(hi > lo))   /* a fixed D: the single point D */
    {
      x[0] = lo;
      n = 1;
    }
  qsort (x, (size_t) n, sizeof (double), jc_cmp_double);
  {   /* drop repeated points */
    long u = 0;
    for (g = 0; g < n; g++)
      if (u == 0 || x[g] > x[u - 1])
        x[u++] = x[g];
    n = u;
  }
  double *row = (double *) mycalloc ((size_t) q->ns, sizeof (double));
  for (g = 0; g < n; g++)
    {
      if (x[g] < lo)
        x[g] = lo;
      if (x[g] > hi)
        x[g] = hi;
      /* a fixed D (one point, n = 1 below) has no prior density */
      const double pr = (hi > lo ? scaling_prior_joint (q->world, pm, x[g]) : 0.0)
        + (ly->ngate > 0 ? jc_gate_group (ly, st, k, x[g]) : 0.0);
      if (ly->split_psig[k] < 0)
        {
          H[g] = pr;
          for (m = 0; m < ly->nsplit; m++)
            if (ly->split_rep[m] == k)
              H[g] += jc_div_term (ly, st, m, x[g], ly->split_sfix[k]);
          continue;
        }
      for (b = 0; b < q->ns; b++)
        {
          if (!(q->sw[k * q->ns + b] > 0.0))
            {
              row[b] = -HUGE_VAL;
              continue;
            }
          double v = q->ps_s[k * q->ns + b] + log (q->sw[k * q->ns + b])
            + jc_div_cv (q, k, q->sv[k * q->ns + b], x[g]);
          for (m = 0; m < ly->nsplit; m++)
            if (ly->split_rep[m] == k)
              v += jc_div_term (ly, st, m, x[g], q->sv[k * q->ns + b]);
          row[b] = v;
        }
      H[g] = pr + jc_logsumexp (row, q->ns);
    }
  myfree (row);
  return n;
}

/* log int_{x[0]}^{a} e^{F} on the grid x (trapezoid, the last cell up to a
   interpolated linearly in e^F) */
static double
jc_log_cum_at (const double *x, const double *F, const double *cum, long n, double a)
{
  long i;
  if (a <= x[0])
    return -HUGE_VAL;
  if (a >= x[n - 1])
    return cum[n - 1];
  for (i = 0; i + 1 < n && x[i + 1] <= a; i++)
    ;
  const double t = (x[i + 1] > x[i]) ? (a - x[i]) / (x[i + 1] - x[i]) : 0.0;
  const double mx = F[i] > F[i + 1] ? F[i] : F[i + 1];
  if (mx == -HUGE_VAL)
    return cum[i];
  const double fi = exp (F[i] - mx), fa = fi + (exp (F[i + 1] - mx) - fi) * t;
  const double part = 0.5 * (fi + fa) * (a - x[i]);
  if (!(part > 0.0))
    return cum[i];
  const double lp = mx + log (part);
  const double m2 = cum[i] > lp ? cum[i] : lp;
  return m2 == -HUGE_VAL ? m2 : m2 + log (exp (cum[i] - m2) + exp (lp - m2));
}

/* the ordered groups: from the tips to the root, F_k(mu) = H_k(mu) + sum over
   the child groups c of log int_{lo_c}^{mu} e^{F_c}; a root adds log int e^{F}.
   Returns the sum over the ordered trees; done[k] marks every group used. */
static double
jc_div_order_normalizer (const jc_layout *ly, const double *st, const jc_divquad *q, boolean *done)
{
  const long ns = ly->nsplit, nmax = jc_div_order_cap (ly);
  double *x = (double *) mycalloc ((size_t) (ns * nmax), sizeof (double));
  double *F = (double *) mycalloc ((size_t) (ns * nmax), sizeof (double));
  double *cum = (double *) mycalloc ((size_t) (ns * nmax), sizeof (double));
  long *n = (long *) mycalloc ((size_t) ns, sizeof (long));
  boolean *inorder = (boolean *) mycalloc ((size_t) ns, sizeof (boolean));
  boolean *ready = (boolean *) mycalloc ((size_t) ns, sizeof (boolean));
  double lz = 0.0;
  long k, c, g, left;
  for (k = 0; k < ns; k++)
    if (ly->split_rep[k] == k && ly->split_par[k] >= 0)
      inorder[k] = inorder[ly->split_par[k]] = TRUE;
  for (k = 0; k < ns; k++)
    if (inorder[k])
      n[k] = jc_div_order_grid (q, ly, st, k, x + k * nmax, F + k * nmax);
  do
    {
      left = 0;
      for (k = 0; k < ns; k++)
        {
          if (!inorder[k] || ready[k])
            continue;
          boolean wait = FALSE;
          for (c = 0; c < ns; c++)
            if (inorder[c] && ly->split_par[c] == k && !ready[c])
              wait = TRUE;
          if (wait)
            {
              left++;
              continue;
            }
          double *xk = x + k * nmax, *Fk = F + k * nmax, *ck = cum + k * nmax;
          for (c = 0; c < ns; c++)
            if (inorder[c] && ly->split_par[c] == k)
              for (g = 0; g < n[k]; g++)
                Fk[g] += jc_log_cum_at (x + c * nmax, F + c * nmax, cum + c * nmax, n[c], xk[g]);
          /* a fixed D is one point: its mass is e^F there (the cumulative
             below it is empty, at and above it the mass itself) */
          ck[0] = n[k] == 1 ? Fk[0] : -HUGE_VAL;
          for (g = 1; g < n[k]; g++)
            {   /* cumulative trapezoid: log of the integral up to xk[g] */
              const double a = Fk[g - 1], bb = Fk[g], w = xk[g] - xk[g - 1];
              const double mx = a > bb ? a : bb;
              const double cell = (mx == -HUGE_VAL || !(w > 0.0)) ? -HUGE_VAL
                : mx + log (0.5 * (exp (a - mx) + exp (bb - mx)) * w);
              const double m2 = ck[g - 1] > cell ? ck[g - 1] : cell;
              ck[g] = m2 == -HUGE_VAL ? m2 : m2 + log (exp (ck[g - 1] - m2) + exp (cell - m2));
            }
          ready[k] = TRUE;
          done[k] = TRUE;
          if (ly->split_par[k] < 0)
            lz += n[k] == 1 ? Fk[0] : jc_logtrapz (Fk, xk, n[k]);
        }
    }
  while (left > 0);
  myfree (x); myfree (F); myfree (cum); myfree (n); myfree (inorder); myfree (ready);
  return lz;
}

/* log of int int prior(mu) prior(sigma) e^{D(mu, sigma)} for every split */
static double
jc_div_normalizer (const jc_layout *ly, const double *st, const jc_divquad *q)
{
  double lz = 0.0;
  long k, b;
  double *row = (double *) mycalloc ((size_t) q->ns, sizeof (double));
  boolean *done = (boolean *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (boolean));
  lz += jc_div_order_normalizer (ly, st, q, done);
  for (k = 0; k < ly->nsplit; k++)
    {
      if (ly->split_rep[k] != k || done[k])   /* integrated with its group's first split, or ordered */
        continue;
      if (ly->split_psig[k] < 0)
        {
          lz += jc_div_mu_integral (q, ly, st, k, ly->split_sfix[k]);
          continue;
        }
      for (b = 0; b < q->ns; b++)
        row[b] = q->sw[k * q->ns + b] > 0.0
          ? q->ps_s[k * q->ns + b] + log (q->sw[k * q->ns + b])
            + jc_div_mu_integral (q, ly, st, k, q->sv[k * q->ns + b])
          : -HUGE_VAL;
      lz += jc_logsumexp (row, q->ns);
    }
  myfree (row);
  myfree (done);
  return lz;
}

/* skyline=PARAM: log R_e(x0) of free entry e on the log-spaced grid x of
   its parameter: the integral over the multipliers tau_1.. of their prior
   times the entry's terms in every segment, the base value x0 = x[i]:

     R_e(x0) = e^{f_0(x0)} int prod_{s>=1} K(u_{s-1}, u_s) e^{f_s(x0 e^{u_s})} du

   in w_s = log(x0 tau_s). The base value x0 is the parameter and stays in
   its prior range (the grid x); a segment value x0 tau_s may leave it (no
   range rule since 7.0.102, as in the sampler): the segments s >= 1 use the
   grid extended by SKYLINE_HISTMAX in log on both sides (|log tau| <=
   SKYLINE_HISTMAX). Before, the segment values were confined to the prior
   range here while the samplers allowed more, so genealogies that call for
   a segment value above the range got a Z(G) too small and dominated the
   joint (two populations with an M skyline: M x tau ran to about 8000 with
   the prior of M up to 1000). The random walk on log tau is a random walk
   on w that starts at w_0 = log x0 (backward recursion with a banded
   Gaussian kernel, rows normalised: a std below the grid spacing pins the
   segments together); the log-uniform prior is independent per segment,
   w_s within w_0 + [log a, log b]. The hard limit |log tau| <=
   SKYLINE_HISTMAX is part of the log-uniform windows; for the random walk
   it is not part of R (its prior mass is negligible unless sigma is large)
   but the joint sampler enforces it. Constant factors of the prior cancel
   in p(G)/Z(G). */
static void
jc_sky_logR (const jc_layout *ly, const double *st, long e, double mu, const double *x, double *R)
{
  const long n = JC_TGRID, j = ly->sky_j[e];
  const double lx0 = log (x[0]);
  const double h = (log (x[n - 1]) - lx0) / (double) (n - 1);
  const long E = (long) ceil (SKYLINE_HISTMAX / h);   /* extension on each side */
  const long ne = n + 2 * E;
  double *xe = (double *) mycalloc ((size_t) ne, sizeof (double));
  double *f = (double *) mycalloc ((size_t) ne, sizeof (double));
  double *B = (double *) mycalloc ((size_t) ne, sizeof (double));
  double *C = (double *) mycalloc ((size_t) ne, sizeof (double));
  long i, q, sg;
  for (i = 0; i < ne; i++)
    xe[i] = (i >= E && i < E + n) ? x[i - E] : exp (lx0 + (double) (i - E) * h);
  if (ly->sky_type[j] == SKYPRIOR_LOGUNIFORM)
    {
      const double la = log (ly->sky_a[j]) > -SKYLINE_HISTMAX ? log (ly->sky_a[j]) : -SKYLINE_HISTMAX;
      const double lb = log (ly->sky_b[j]) < SKYLINE_HISTMAX ? log (ly->sky_b[j]) : SKYLINE_HISTMAX;
      for (i = 0; i < n; i++)
        R[i] = jc_sky_seg_term (ly, st, e, 0, x[i], mu);
      for (sg = 1; sg < ly->nseg; sg++)
        {
          double fmax = -HUGE_VAL;
          for (i = 0; i < ne; i++)
            {
              f[i] = jc_sky_seg_term (ly, st, e, sg, xe[i], mu);
              if (f[i] > fmax)
                fmax = f[i];
            }
          C[0] = 0.0;   /* cumulative trapezoid of e^{f - fmax} in w */
          for (i = 1; i < ne; i++)
            C[i] = C[i - 1] + 0.5 * h * (exp (f[i - 1] - fmax) + exp (f[i] - fmax));
          for (i = 0; i < n; i++)
            {
              /* window [w_0 + log a, w_0 + log b] in extended-grid units */
              double a = (double) (i + E) + la / h, b = (double) (i + E) + lb / h, ca, cb;
              if (a < 0.0) a = 0.0;
              if (b > (double) (ne - 1)) b = (double) (ne - 1);
              if (b <= a || fmax == -HUGE_VAL)
                {
                  R[i] = -HUGE_VAL;
                  continue;
                }
              q = (long) a; if (q > ne - 2) q = ne - 2;
              ca = C[q] + (a - q) * (C[q + 1] - C[q]);
              q = (long) b; if (q > ne - 2) q = ne - 2;
              cb = C[q] + (b - q) * (C[q + 1] - C[q]);
              R[i] += (cb > ca) ? fmax + log (cb - ca) : -HUGE_VAL;
            }
        }
    }
  else
    {
      /* random walk with std sigma on w */
      const double sigma = ly->sky_a[j];
      long band = (long) ceil (6.0 * sigma / h);
      if (band > ne - 1)
        band = ne - 1;
      double *kw = (double *) mycalloc ((size_t) (2 * band + 1), sizeof (double));
      double ksum = 0.0;
      for (q = -band; q <= band; q++)
        ksum += (kw[q + band] = exp (-0.5 * (q * h / sigma) * (q * h / sigma)));
      for (q = 0; q <= 2 * band; q++)
        kw[q] /= ksum;
      /* B_s(w) = f_s(w) + log int K(w, w') e^{B_{s+1}(w')} dw', s = nseg-1 .. 1,
         on the extended grid; then R(x0) = f_0(x0) + log int K e^{B_1} */
      for (i = 0; i < ne; i++)
        B[i] = jc_sky_seg_term (ly, st, e, ly->nseg - 1, xe[i], mu);
      for (sg = ly->nseg - 2; sg >= 0; sg--)
        {
          double bmax = -HUGE_VAL;
          const long from = sg == 0 ? E : 0, to = sg == 0 ? E + n : ne;
          for (i = 0; i < ne; i++)
            if (B[i] > bmax)
              bmax = B[i];
          for (i = 0; i < ne; i++)
            C[i] = bmax == -HUGE_VAL ? 0.0 : exp (B[i] - bmax);
          for (i = from; i < to; i++)
            {
              double sum = 0.0;
              const long lo = i - band < 0 ? 0 : i - band, hi = i + band > ne - 1 ? ne - 1 : i + band;
              for (q = lo; q <= hi; q++)
                sum += kw[q - i + band] * C[q];
              f[i] = jc_sky_seg_term (ly, st, e, sg, xe[i], mu) + (sum > 0.0 ? bmax + log (sum) : -HUGE_VAL);
            }
          for (i = 0; i < ne; i++)
            B[i] = (i >= from && i < to) ? f[i] : -HUGE_VAL;
        }
      for (i = 0; i < n; i++)
        R[i] = B[i + E];
      myfree (kw);
    }
  myfree (xe);
  myfree (f);
  myfree (B);
  myfree (C);
}

/* TRUE when parameter p is the parameter of a free skyline entry */
static boolean
jc_sky_param (const jc_layout *ly, long p)
{
  long j;
  for (j = 0; j < ly->nfree; j++)
    if (ly->rep[ly->sky_k[j]] == p)
      return TRUE;
  return FALSE;
}

/* log normalizer of one genealogy: product of independent blocks */
static double
jc_log_normalizer (const jc_layout *ly, const boolean *active,
                   const double *st, double mu, const double *tgrid, const double *tpri,
                   const double *tgrid_s, const double *tpri_s, const double *gpri,
                   const jc_divquad *dq)
{
  const long numpop = ly->numpop, numpop2 = ly->numpop2;
  double v[JC_GRID], w[JC_NG], u[JC_TGRID], z[JC_TGRID];
  double lz = 0.0;
  long p, g, s, pop, t, e, q;
  /* xNm parameters (use-M=NO, custom-migration S/M): M_ij = xNm/Theta_j ties
     the parameter to the Theta of its target population. A population with
     xNm entries gets R_pop(Theta) = prior + its coalescence term + every
     unshared xNm entry already integrated over its xNm; an xNm shared by
     several entries (S pair, M group, use-M=NO s/m) is integrated outside
     the Theta integrals of all populations it reaches (jc_supported()
     admits only populations reached by at most one shared xNm) */
  boolean *xpop = (boolean *) mycalloc ((size_t) numpop, sizeof (boolean));
  boolean *inmulti = (boolean *) mycalloc ((size_t) numpop, sizeof (boolean));
  long *nent = (long *) mycalloc ((size_t) numpop2, sizeof (long));
  double *R = NULL;
  for (e = numpop; e < numpop2; e++)
    if (ly->xnm[e] && ly->rep[e] >= 0 && active[ly->rep[e]])
      {
        long from, to;
        m2mm (e, numpop, &from, &to);
        xpop[to] = TRUE;
        nent[ly->rep[e]]++;
      }
  for (e = numpop; e < numpop2; e++)
    if (ly->xnm[e] && ly->rep[e] >= 0 && active[ly->rep[e]] && nent[ly->rep[e]] > 1)
      {
        long from, to;
        m2mm (e, numpop, &from, &to);
        inmulti[to] = TRUE;
      }
  for (pop = 0; pop < numpop; pop++)
    if (xpop[pop])
      {
        if (R == NULL)
          R = (double *) mycalloc ((size_t) (numpop * JC_TGRID), sizeof (double));
        for (t = 0; t < JC_TGRID; t++)
          {
            const double th = tgrid_s[pop * JC_TGRID + t];
            double *Rt = R + pop * JC_TGRID + t;
            *Rt = (active[pop] ? tpri_s[pop * JC_TGRID + t] : 0.0)
              + jc_pop_term (ly, st, pop, th, 0.0, mu, NAN);
            for (e = mstart (pop, numpop); e < mend (pop, numpop); e++)
              if (ly->xnm[e] && ly->rep[e] >= 0 && active[ly->rep[e]] && nent[ly->rep[e]] == 1)
                {
                  const long r = ly->rep[e];
                  for (q = 0; q < JC_TGRID; q++)
                    z[q] = tpri_s[r * JC_TGRID + q] + jc_mig_term (ly, st, e, tgrid_s[r * JC_TGRID + q] / th, mu);
                  *Rt += jc_logtrapz (z, tgrid_s + r * JC_TGRID, JC_TGRID);
                }
          }
        if (!inmulti[pop])
          lz += jc_logtrapz (R + pop * JC_TGRID, tgrid_s + pop * JC_TGRID, JC_TGRID);
      }
  /* shared xNm parameters: int prior(x) prod_pop int R_pop(Theta) e^{sum h_e(x/Theta)} */
  for (p = numpop; p < numpop2; p++)
    if (nent[p] > 1)
      {
        for (q = 0; q < JC_TGRID; q++)
          {
            const double x = tgrid_s[p * JC_TGRID + q];
            z[q] = tpri_s[p * JC_TGRID + q];
            for (pop = 0; pop < numpop; pop++)
              {
                boolean hit = FALSE;
                for (t = 0; t < JC_TGRID; t++)
                  u[t] = R[pop * JC_TGRID + t];
                for (e = mstart (pop, numpop); e < mend (pop, numpop); e++)
                  if (ly->xnm[e] && ly->rep[e] == p)
                    {
                      hit = TRUE;
                      for (t = 0; t < JC_TGRID; t++)
                        u[t] += jc_mig_term (ly, st, e, x / tgrid_s[pop * JC_TGRID + t], mu);
                    }
                if (hit)
                  z[q] += jc_logtrapz (u, tgrid_s + pop * JC_TGRID, JC_TGRID);
              }
          }
        lz += jc_logtrapz (z, tgrid_s + p * JC_TGRID, JC_TGRID);
      }
  myfree (R);
  myfree (nent);
  myfree (inmulti);
  /* skyline=PARAM: a parameter of free entries is integrated on its
     JC_TGRID grid, each free entry contributing R_e (jc_sky_logR()) */
  for (p = 0; p < numpop2 && ly->nfree > 0; p++)
    if (active[p] && jc_sky_param (ly, p) && !(p < numpop && ly->off_sA[p] >= 0))
      {
        const double *x = tgrid_s + p * JC_TGRID;
        double Re[JC_TGRID];
        for (t = 0; t < JC_TGRID; t++)
          z[t] = tpri_s[p * JC_TGRID + t];
        for (e = 0; e < numpop2; e++)
          if (ly->rep[e] == p)
            {
              if (ly->sky_j[e] >= 0)
                {
                  jc_sky_logR (ly, st, e, mu, x, Re);
                  for (t = 0; t < JC_TGRID; t++)
                    z[t] += Re[t];
                }
              else
                for (t = 0; t < JC_TGRID; t++)
                  z[t] += e < numpop ? jc_pop_term (ly, st, e, x[t], 0.0, mu, NAN)
                    : jc_mig_term (ly, st, e, x[t], mu);
            }
        lz += jc_logtrapz (z, x, JC_TGRID);
      }
  /* Theta parameters without growth: all populations that share it */
  for (p = 0; p < numpop; p++)
    if (active[p] && !xpop[p] && !(ly->pop_slot[p] >= 0 && ly->rep[p] == p) && !jc_sky_param (ly, p))
      {
        for (g = 0; g < JC_GRID; g++)
          {
            v[g] = tpri[p * JC_GRID + g];
            for (pop = 0; pop < numpop; pop++)
              if (ly->rep[pop] == p)
                v[g] += jc_pop_term (ly, st, pop, tgrid[p * JC_GRID + g], 0.0, mu, NAN);
          }
        lz += jc_logtrapz (v, tgrid + p * JC_GRID, JC_GRID);
      }
  myfree (xpop);
  /* M parameters: all entries that share it (xNm ones are done above) */
  for (p = numpop; p < numpop2; p++)
    if (active[p] && !ly->xnm[p] && !jc_sky_param (ly, p) && !jc_gate_param (ly, p))
      {
        for (g = 0; g < JC_GRID; g++)
          {
            v[g] = tpri[p * JC_GRID + g];
            for (e = numpop; e < numpop2; e++)
              if (ly->rep[e] == p)
                v[g] += jc_mig_term (ly, st, e, tgrid[p * JC_GRID + g], mu);
          }
        lz += jc_logtrapz (v, tgrid + p * JC_GRID, JC_GRID);
      }
  /* growth groups: int prior(g) prod_pop int prior(Theta) e^{f(Theta, g)};
     a population with skyline multipliers contributes R(Theta) of its
     multiplier recursion at each growth node (every JC_SKYG-th node) */
  for (s = 0; s < ly->nslot; s++)
    {
      const double *gg = ly->ggrid + s * JC_NG;
      long stride = 1, ng = 0;
      double wg[JC_NG], xg[JC_NG];
      for (pop = 0; pop < numpop; pop++)
        if (ly->pop_slot[pop] == s && active[pop] && ly->off_sA[pop] >= 0)
          stride = JC_SKYG;
      for (g = 0; g < JC_NG; g += stride)
        {
          w[g] = gpri[s * JC_NG + g];
          for (pop = 0; pop < numpop; pop++)
            if (ly->pop_slot[pop] == s && active[pop])
              {
                if (ly->off_sA[pop] >= 0)
                  {
                    double Re[JC_TGRID];
                    jc_zgidx = g;
                    jc_sky_logR (ly, st, pop, mu, tgrid_s + pop * JC_TGRID, Re);
                    jc_zgidx = -1;
                    for (t = 0; t < JC_TGRID; t++)
                      u[t] = tpri_s[pop * JC_TGRID + t] + Re[t];
                  }
                else
                  {
                    const double Ag = st[ly->off_A[pop] + g];
                    for (t = 0; t < JC_TGRID; t++)
                      u[t] = tpri_s[pop * JC_TGRID + t]
                        + jc_pop_term (ly, st, pop, tgrid_s[pop * JC_TGRID + t], gg[g], mu, Ag);
                  }
                w[g] += jc_logtrapz (u, tgrid_s + pop * JC_TGRID, JC_TGRID);
              }
          wg[ng] = w[g];
          xg[ng++] = gg[g];
        }
      lz += jc_logtrapz (wg, xg, ng);
    }
  if (jc_nodiv)   /* sharp splits: the divergence part is integrated outside */
    return lz;
  if (ly->ngate > 0)   /* the gated M integrals, used inside the divergence integrals */
    jc_gate_table (ly, st, mu);
  lz += jc_div_normalizer (ly, st, dq);
  return lz;
}

/* integrated crossing: the row of the genealogy with the crossing at c
   (from the base row, where it is at a): K and S of daughter and ancestor,
   the gated lineage times, and the divergence grid of the split with the
   lineage's waiting terms from a and its split event at c */
static void
jc_rb_move (world_fmt *world, const jc_layout *ly, const double *base, double c, double *out)
{
  const double *R = base + ly->off_R;
  const long k = (long) R[0], np = (long) R[3];
  const double a = R[1], b = R[2];
  const long d = ly->split_to[k], an = ly->split_from[k];
  long j, e, g, from, to;
  memcpy (out, base, sizeof (double) * (size_t) ly->nrow);
  for (j = 0; j < np; j++)
    {
      const double u = R[4 + 3 * j], v = j + 1 < np ? R[4 + 3 * (j + 1)] : b;
      const double len = (v < c ? v : c) - u;
      if (len <= 0.0)
        break;
      out[ly->off_K + d] += R[4 + 3 * j + 1] * len;
      out[ly->off_K + an] -= R[4 + 3 * j + 2] * len;
      if (ly->nseg > 1)
        {   /* skyline=PARAM: the piece's segment */
          const long sg = jc_seg (ly, u);
          out[ly->off_sK + sg * ly->numpop + d] += R[4 + 3 * j + 1] * len;
          out[ly->off_sK + sg * ly->numpop + an] -= R[4 + 3 * j + 2] * len;
          out[ly->off_sS + sg * ly->numpop + d] += len;
          out[ly->off_sS + sg * ly->numpop + an] -= len;
        }
    }
  out[ly->off_S + d] += c - a;
  out[ly->off_S + an] -= c - a;
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->gate_k[e] >= 0)
      {
        m2mm (e, ly->numpop, &from, &to);
        if (to != d && to != an)
          continue;
        for (g = 0; g < JC_NGATE; g++)
          {
            const double x = jc_gate_x (ly, e, g);
            const double ov = (c < x ? c : x) - a;
            if (ov > 0.0)
              out[ly->off_G[e] + 1 + g] += to == d ? ov : -ov;
          }
      }
  species_fmt *sp = &world->species_model[ly->split_model[k]];
  out[ly->off_D[k]] = c;
  out[ly->off_D[k] + 1] = base[ly->off_D[k] + 1] + 1.0;
  for (j = 0; j < ly->split_ns[k]; j++)
    {
      double L, U;
      const double sig = ly->sgrid[k * JC_NS + j];
      double *Dv = out + ly->off_D[k] + 2 + j * JC_NM;
      jc_row_window (ly, out, k, j, &L, &U);
      for (g = 0; g < JC_NM; g++)
        {
          const double m = L + (U - L) * (double) g / (JC_NM - 1);
          const double rest = jc_row_value (ly, base, k, j, m);
          Dv[g] = rest <= JC_DOUT ? JC_DOUT
            : rest + (*log_prob_wait_speciate) (a, c, m, sig, sp) + (*log_point_prob_speciate) (c, m, sig, sp);
          if (!(Dv[g] > JC_DOUT))
            Dv[g] = JC_DOUT;
        }
    }
  out[ly->off_R] = -1.0;
}

static void jc_gl_nodes (long n, double *x, double *w);

/* sharp splits: log of the S prior masses (the empty D blocks integrate
   them in the other rows' normalizers) */
static double
jc_sharp_logmS (world_fmt *world, const jc_layout *ly)
{
  double lm = 0.0;
  long k, i;
  for (k = 0; k < ly->nsplit; k++)
    {
      const long ps = ly->split_psig[k];
      if (ps < 0 || ly->split_rep[k] != k)   /* a tied group shares one S */
        continue;
      const double lo = world->bayes->minparam[ps], hi = world->bayes->maxparam[ps];
      if (!(hi > lo))
        continue;
      const long N = 201;
      double v[201], x[201];
      for (i = 0; i < N; i++)
        {
          x[i] = lo + (hi - lo) * (double) i / (double) (N - 1);
          v[i] = scaling_prior_joint (world, ps, x[i]);
        }
      lm += jc_logtrapz (v, x, N);
    }
  return lm;
}

typedef struct
{
  world_fmt *world;
  const jc_layout *ly;
  const boolean *active;
  const double *st;
  double mu;
  const double *tgrid, *tpri, *tgrid_s, *tpri_s, *gpri;
  const jc_divquad *dq;
  long order[64];
  long nlev;          /* integration levels: one per tied group */
  double c[64];
  double *mv;
  long Q;
  double x[64], w[64];
} jc_sharpq;

/* nested Gauss-Legendre over the splits in order (parents first): a
   child's range ends at its parent's time when they are ordered (chain or
   divergence-order) */
static double
jc_sharp_nest (jc_sharpq *z, long level)
{
  const jc_layout *ly = z->ly;
  const double *R = z->st + ly->off_R;
  if (level == z->nlev)
    {
      double dK[64], dS[64], v = 0.0;
      long k, p;
      jc_sharp_move_row (ly, z->st, z->c, z->mv, dK, dS);
      (void) p;
      for (k = 0; k < ly->nsplit; k++)
        if (ly->split_rep[k] == k)   /* a tied group shares one D */
          v += scaling_prior_joint (z->world, ly->split_pmu[k], z->c[k]);
      v += jc_sharp_lin (ly, z->active, z->world->param0, dK, dS, z->mu, 2);
      jc_nodiv = TRUE;
      v += jc_log_normalizer (ly, z->active, z->mv, z->mu, z->tgrid, z->tpri, z->tgrid_s, z->tpri_s, z->gpri, z->dq);
      jc_nodiv = FALSE;
      return v;
    }
  const long k = z->order[level];
  const double *Sk = R + JC_SH_SPLIT (k);
  double lo = -HUGE_VAL, hi = HUGE_VAL;
  long m;
  boolean anyx = FALSE;
  /* a tied group (t): one time for all its splits, inside every member's
     window and below every member's parent */
  for (m = 0; m < ly->nsplit; m++)
    {
      if (ly->split_rep[m] != k)
        continue;
      const double *Sm = R + JC_SH_SPLIT (m);
      const double l = Sm[0] > ly->mlo[m] ? Sm[0] : ly->mlo[m];
      const double h = Sm[2] > 0.0 && Sm[1] < ly->mhi[m] ? Sm[1] : ly->mhi[m];
      const long par = ly->sharp_par[m];
      if (l > lo) lo = l;
      if (h < hi) hi = h;
      if (par >= 0 && (ly->split_par[m] >= 0 || Sm[4] > 0.0) && z->c[par] < hi)
        hi = z->c[par];
      if (Sm[2] > 0.0)
        anyx = TRUE;
    }
  if (!(hi > lo))
    return -HUGE_VAL;
  /* one split: Gauss-Legendre per piece between the row's breakpoints */
  double cuts[JC_SHARPMAXI + 10];
  long nc = 0, i, q;
  cuts[nc++] = lo;
  if (ly->nsplit == 1)
    {
      const double *tlr = R + JC_SH_TL (ly);
      for (i = 0; i < (long) R[1]; i++)
        {
          const double t = tlr[i * (1 + ly->numpop)];
          if (t > lo && t < hi)
            cuts[nc++] = t;
        }
    }
  (void) Sk;
  if (!anyx && lo > 0.0)
    for (i = 1; i < 8; i++)   /* no crossing: the prior alone, over a wide range */
      cuts[nc++] = lo * pow (hi / lo, (double) i / 8.0);
  cuts[nc++] = hi;
  double acc = -HUGE_VAL;
  for (i = 0; i + 1 < nc; i++)
    {
      const double u = cuts[i], v = cuts[i + 1];
      if (!(v > u))
        continue;
      for (q = 0; q < z->Q; q++)
        {
          z->c[k] = 0.5 * (u + v) + 0.5 * (v - u) * z->x[q];
          for (m = 0; m < ly->nsplit; m++)
            if (ly->split_rep[m] == k)
              z->c[m] = z->c[k];
          acc = jc_logadd (acc, log (0.5 * (v - u) * z->w[q]) + jc_sharp_nest (z, level + 1));
        }
    }
  return acc;
}

/* sharp splits: log int p(D = c) Z(G_c) e^{F(c)} dc over the splits'
   windows (F: the fixed parameters' terms), times the S prior masses; a row
   without crossings has the usual normalizer */
static double
jc_sharp_lognorm (world_fmt *world, const jc_layout *ly, const boolean *active,
                  const double *st, double mu, const double *tgrid, const double *tpri,
                  const double *tgrid_s, const double *tpri_s, const double *gpri,
                  const jc_divquad *dq)
{
  const double *R = st + ly->off_R;
  jc_sharpq z;
  long k, n = 0;
  if (R[0] == -1.0)
    return jc_log_normalizer (ly, active, st, mu, tgrid, tpri, tgrid_s, tpri_s, gpri, dq);
  if (R[0] < 0.0)
    return 0.0;
  memset (&z, 0, sizeof (z));
  z.world = world;
  z.ly = ly;
  z.active = active;
  z.st = st;
  z.mu = mu;
  z.tgrid = tgrid; z.tpri = tpri; z.tgrid_s = tgrid_s; z.tpri_s = tpri_s; z.gpri = gpri;
  z.dq = dq;
  z.Q = getenv ("MIGRATE_JC_SHARPQ") ? atol (getenv ("MIGRATE_JC_SHARPQ")) : 8;
  if (z.Q < 2) z.Q = 2;
  if (z.Q > 64) z.Q = 64;
  jc_gl_nodes (z.Q, z.x, z.w);
  /* parents first; one level per tied group (its first split), entered
     when the groups of all its members' parents are done */
  boolean *done = (boolean *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (boolean));
  long nrep = 0;
  for (k = 0; k < ly->nsplit; k++)
    if (ly->split_rep[k] == k)
      nrep++;
  while (n < nrep)
    {
      const long n0 = n;
      for (k = 0; k < ly->nsplit; k++)
        {
          long m;
          boolean ready = ly->split_rep[k] == k && !done[k];
          for (m = 0; m < ly->nsplit && ready; m++)
            if (ly->split_rep[m] == k && ly->sharp_par[m] >= 0 && !done[ly->split_rep[ly->sharp_par[m]]])
              ready = FALSE;
          if (ready)
            {
              z.order[n++] = k;
              done[k] = TRUE;
            }
        }
      if (n == n0)
        break;
    }
  myfree (done);
  z.nlev = n;
  z.mv = (double *) mycalloc ((size_t) ly->nrow, sizeof (double));
  const double acc = n == nrep ? jc_sharp_nest (&z, 0) : -HUGE_VAL;
  myfree (z.mv);
  return acc > -HUGE_VAL ? acc + jc_sharp_logmS (world, ly) : 0.0;
}

/* Gauss-Legendre nodes and weights on [-1, 1] */
static void
jc_gl_nodes (long n, double *x, double *w)
{
  long i, j;
  for (i = 0; i < (n + 1) / 2; i++)
    {
      double z = cos (M_PI * ((double) i + 0.75) / ((double) n + 0.5)), pp = 1.0, z1;
      do
        {
          double p1 = 1.0, p2 = 0.0;
          for (j = 0; j < n; j++)
            {
              const double p3 = p2;
              p2 = p1;
              p1 = ((2.0 * (double) j + 1.0) * z * p2 - (double) j * p3) / ((double) j + 1.0);
            }
          pp = (double) n * (z * p1 - p2) / (z * z - 1.0);
          z1 = z;
          z = z1 - p1 / pp;
        }
      while (fabs (z - z1) > 1e-15);
      x[i] = -z;
      x[n - 1 - i] = z;
      w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
}

/* Gauss-Legendre nodes of the crossing's normalizer over the whole segment
   (MIGRATE_JC_RBQ; 16 nodes: log Z within 0.007 of 64 nodes, median 2e-5,
   twopop-div-sbc 2026-11-18; MIGRATE_JC_RBPIECES=1: per piece) */
static long
jc_rb_q (void)
{
  const char *e = getenv ("MIGRATE_JC_RBQ");
  long q = e != NULL ? atol (e) : 16;
  return q < 2 ? 2 : (q > 64 ? 64 : q);
}

/* integrated crossing: log int_a^b Z(G_c) e^{F(c)} dc, Gauss-Legendre over
   (a, b) cut at the prior bounds of the split mean (Z(G_c) is smooth, with
   kinks at the piece ends); F(c)
   are the terms of the fixed parameters, relative to c = a */
static double
jc_rb_lognorm (world_fmt *world, const jc_layout *ly, const boolean *active,
               const double *st, double mu, const double *tgrid, const double *tpri,
               const double *tgrid_s, const double *tpri_s, const double *gpri,
               const jc_divquad *dq)
{
  const double *R = st + ly->off_R;
  const long k = (long) R[0], np = (long) R[3];
  const double a = R[1], b = R[2];
  const long Q = jc_rb_q ();
  double x[64], w[64], cuts[2 * JC_RBMAXP + 4];
  double *mv = (double *) mycalloc ((size_t) ly->nrow, sizeof (double));
  const double *phi0 = world->param0;
  long j, q, n = 0, i;
  double acc = -HUGE_VAL;
  jc_gl_nodes (Q, x, w);
  for (j = 0; j < np; j++)
    if (j == 0 || getenv ("MIGRATE_JC_RBPIECES") != NULL)
      cuts[n++] = R[4 + 3 * j];
  if (ly->mlo[k] > a && ly->mlo[k] < b)
    cuts[n++] = ly->mlo[k];
  if (ly->mhi[k] > a && ly->mhi[k] < b)
    cuts[n++] = ly->mhi[k];
  cuts[n++] = b;
  for (i = 1; i < n; i++)
    for (j = i; j > 0 && cuts[j - 1] > cuts[j]; j--)
      {
        const double t = cuts[j];
        cuts[j] = cuts[j - 1];
        cuts[j - 1] = t;
      }
  for (i = 0; i + 1 < n; i++)
    {
      const double u = cuts[i], v = cuts[i + 1];
      if (!(v > u))
        continue;
      for (q = 0; q < Q; q++)
        {
          const double c = 0.5 * (u + v) + 0.5 * (v - u) * x[q];
          double F = 0.0;
          for (j = 0; j < np; j++)
            {   /* fixed-parameter terms from a to c */
              const double pu = R[4 + 3 * j], pv = j + 1 < np ? R[4 + 3 * (j + 1)] : b;
              const double len = (pv < c ? pv : c) - pu;
              if (len <= 0.0)
                break;
              F += len * jc_rb_slope (ly, active, phi0, k, R[4 + 3 * j + 1], R[4 + 3 * j + 2],
                                      pu + 0.5 * len, mu, 2, 1.0);
            }
          jc_rb_move (world, ly, st, c, mv);
          acc = jc_logadd (acc, log (0.5 * (v - u) * w[q]) + F
                           + jc_log_normalizer (ly, active, mv, mu, tgrid, tpri, tgrid_s, tpri_s, gpri, dq));
        }
    }
  myfree (mv);
  return acc;
}

/* MIGRATE_JC_RBCHECK: the closed form of the integrated crossing against a
   numerical integral over c of the rows moved to c (the coalescence and
   migration terms of every population and entry, at the sampler's
   parameters) plus the sampler's own split terms */
static double
jc_rb_check_terms (const jc_layout *ly, const double *st, const double *phi, double mu)
{
  double v = 0.0;
  long pop, e, sg;
  for (pop = 0; pop < ly->numpop; pop++)
    if (ly->nseg > 1 && ly->sky_j[pop] >= 0)
      for (sg = 0; sg < ly->nseg; sg++)
        v += jc_pop_seg_term (ly, st, pop, sg, jc_value (ly, phi, pop) * jc_rb_tau (ly, phi, pop, sg), 0.0, mu);
    else
      v += jc_pop_term (ly, st, pop, jc_value (ly, phi, pop), 0.0, mu, NAN);
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->nseg > 1 && ly->sky_j[e] >= 0)
      for (sg = 0; sg < ly->nseg; sg++)
        v += jc_mig_seg_term (ly, st, e, sg, jc_rate (ly, phi, e) * jc_rb_tau (ly, phi, e, sg), mu);
    else if (!ly->mzero[e] && jc_rate (ly, phi, e) > 0.0)
      v += ly->gate_k[e] >= 0
        ? jc_gate_term (ly, st, e, jc_rate (ly, phi, e), phi[ly->split_pmu[ly->gate_k[e]]], mu)
        : jc_mig_term (ly, st, e, jc_rate (ly, phi, e), mu);
  return v;
}

static void
jc_rb_check (world_fmt *world, const jc_layout *ly, const double *st)
{
  const double *R = st + ly->off_R;
  const long k = (long) R[0];
  const double a = R[1], b = R[2], mu = st[ly->off_mu];
  double *phi = (double *) mycalloc ((size_t) ly->nphi, sizeof (double));
  long jj, sg2;
  memcpy (phi, world->param0, sizeof (double) * (size_t) ly->np);
  for (jj = 0; jj < ly->nfree; jj++)   /* the sampler's skyline multipliers */
    for (sg2 = 1; sg2 < ly->nseg; sg2++)
      phi[JC_TAU (ly, jj, sg2)] = world->timek[sg2 * world->numpop2 + ly->sky_k[jj]];
  const double D = phi[ly->split_pmu[k]];
  const double S = ly->split_psig[k] >= 0 ? phi[ly->split_psig[k]] : ly->split_sfix[k];
  species_fmt *sp = &world->species_model[ly->split_model[k]];
  const long n1 = 4000, n2 = 4000;
  double *x = (double *) mycalloc ((size_t) (n1 + n2 + 2), sizeof (double));
  double *f = (double *) mycalloc ((size_t) (n1 + n2 + 2), sizeof (double));
  double *mv = (double *) mycalloc ((size_t) ly->nrow, sizeof (double));
  long n = 0, i, j;
  for (i = 0; i <= n1; i++)
    x[n++] = a + (b - a) * (double) i / n1;
  for (i = 0; i <= n2 && jc_rb_normal; i++)
    {
      const double y = D - 12.0 * S + 24.0 * S * (double) i / n2;
      if (y > a && y < b)
        x[n++] = y;
    }
  for (i = 1; i < n; i++)
    for (j = i; j > 0 && x[j - 1] > x[j]; j--)
      {
        const double t = x[j];
        x[j] = x[j - 1];
        x[j - 1] = t;
      }
  const double t0 = jc_rb_check_terms (ly, st, phi, mu);
  for (i = 0; i < n; i++)
    {
      jc_rb_move (world, ly, st, x[i], mv);
      f[i] = jc_rb_check_terms (ly, mv, phi, mu) - t0
        + (x[i] > a ? (*log_prob_wait_speciate) (a, x[i], D, S, sp) : 0.0)
        + (*log_point_prob_speciate) (x[i], D, S, sp);
    }
  const double num = jc_logtrapz (f, x, n);
  const double cf = jc_rb_term (ly, st, phi, mu);
  fprintf (stderr, "JC_RBCHECK locus %ld split %ld a %g b %g c0 %g pieces %ld D %g S %g: closed %.6f numeric %.6f diff %.2e\n",
           world->locus + 1, k, a, b, jc_rbc.c0, (long) R[3], D, S, cf, num, cf - num);
  myfree (x);
  myfree (f);
  myfree (mv);
  myfree (phi);
}

/* the genealogies a combination uses: 0 all, k = 1..JC_NBLOCK the k-th
   block of every locus (set on every rank before jc_local_setup()) */
static long jc_part = 0;
static long jc_nb = JC_NBLOCK;      /* number of blocks */
static boolean jc_byrep = FALSE;    /* blocks are replicates (r mod jc_nb) */

/* blocks: with replicate=YES:N (N >= 2) the replicates -- independent
   chains, so the error includes differences between runs -- grouped into
   at most JC_MAXREPBLOCK blocks; otherwise JC_NBLOCK consecutive stretches */
#define JC_MAXREPBLOCK 10
static long
jc_nblocks (world_fmt *world)
{
  const long reps = world->options->replicate ? world->options->replicatenum : 1;
  if (reps >= 2)
    return reps < JC_MAXREPBLOCK ? reps : JC_MAXREPBLOCK;
  return JC_NBLOCK;
}

static void
jc_set_part (world_fmt *world, long part)
{
  const long reps = world->options->replicate ? world->options->replicatenum : 1;
  jc_part = part;
  jc_byrep = reps >= 2;
  jc_nb = jc_nblocks (world);
}

/* contiguous blocks: the range of row indices of the current part */
static void
jc_part_range (long n, long *a, long *len)
{
  if (jc_part <= 0 || jc_byrep)
    {
      *a = 0;
      *len = n;
      return;
    }
  *a = (n * (jc_part - 1)) / jc_nb;
  *len = (n * jc_part) / jc_nb - *a;
}

/* does row number idx (of n) belong to the current part */
static boolean
jc_row_in_part (const double *row, long nrow, long idx, long a, long len)
{
  if (jc_part <= 0)
    return TRUE;
  if (jc_byrep)
    return ((long) row[nrow - 1]) % jc_nb == jc_part - 1;
  return idx >= a && idx < a + len;
}

/* the spells of row st (Mittag-Leffler) */
static long
jc_row_nsp (const jc_layout *ly, const double *st)
{
  long pop, n = 0;
  if (!ly->ml)
    return 0;
  for (pop = 0; pop < ly->numpop; pop++)
    n += (long) st[ly->off_nsp + pop];
  return n;
}

/* appends n spells to L (capacity in spells: *cap) */
static void
jc_locus_add_spells (jc_locus *L, long *nsp, long *cap, const double *sp, long n)
{
  if (*nsp + n > *cap)
    {
      *cap = 2 * (*nsp + n) + 1024;
      L->sp = (double *) myrealloc (L->sp, sizeof (double) * (size_t) (2 * *cap));
    }
  if (n > 0)
    memcpy (L->sp + 2 * *nsp, sp, sizeof (double) * (size_t) (2 * n));
  *nsp += n;
}

static long
jc_read_locus_files (world_fmt *world, const jc_layout *ly, long locus, jc_locus *L)
{
  const long nrow = ly->nrow;
  long spn = 0, spcap = 0, spbuf = 0;
  double *spl = NULL;
  const long reps = (world->options->replicate && world->options->replicatenum > 0)
    ? world->options->replicatenum : 1;
  long n = 0, r, i, row = 0, kept = 0, ra = 0, rlen = 0, mi = 0;
  double *buf = (double *) mycalloc ((size_t) nrow, sizeof (double));
  for (int pass = 0; pass < 2; pass++)
    {
      long keep = 0;
      if (pass == 1)
        {
          if (n == 0)
            break;
          jc_part_range (n, &ra, &rlen);
          keep = n < jc_maxsamples (world, nrow) ? n : jc_maxsamples (world, nrow);   /* rows of the full combination */
          L->st = (double *) mycalloc ((size_t) (keep * nrow), sizeof (double));
          L->src = (long *) mycalloc ((size_t) (keep + 1), sizeof (long));
          if (ly->ml)
            L->spo = (long *) mycalloc ((size_t) (keep + 1), sizeof (long));
        }
      row = 0;
      for (r = 0; r < reps; r++)
        {
          char name[LINESIZE];
          FILE *f;
          jc_filename (world, locus, r, name);
          if (getenv ("MIGRATE_JC_READPREFIX") != NULL)   /* experiments: rows from elsewhere */
            snprintf (name, LINESIZE, "%s.joint.%ld.%ld", getenv ("MIGRATE_JC_READPREFIX"), locus, r);
          if ((f = fopen (name, "r")) == NULL)
            continue;
          for (;;)
            {
              long k;
              for (i = 0; i < nrow; i++)
                if (fscanf (f, "%lf", &buf[i]) != 1)
                  break;
              if (i < nrow)
                break;
              k = jc_row_nsp (ly, buf);   /* Mittag-Leffler: the spells follow */
              if (k > spbuf)
                {
                  spbuf = 2 * k + 256;
                  spl = (double *) myrealloc (spl, sizeof (double) * (size_t) (2 * spbuf));
                }
              for (i = 0; i < 2 * k; i++)
                if (fscanf (f, "%lf", &spl[i]) != 1)
                  break;
              if (i < 2 * k)
                break;
              if (pass == 0)
                n++;
              else if (mi < keep && row == (mi * n) / keep)
                {   /* a row of the full combination; a block keeps its own */
                  if (jc_row_in_part (buf, nrow, row, ra, rlen))
                    {
                      if (ly->ml)
                        {
                          L->spo[kept] = spn;
                          jc_locus_add_spells (L, &spn, &spcap, spl, k);
                        }
                      L->src[kept] = row;
                      memcpy (L->st + (kept++) * nrow, buf, sizeof (double) * (size_t) nrow);
                    }
                  mi++;
                }
              row++;
            }
          fclose (f);
        }
    }
  myfree (buf);
  myfree (spl);
  if (kept == 0 && L->st != NULL)
    {
      myfree (L->st);
      L->st = NULL;
    }
  L->total = n;
  L->n = kept;
  L->mu = world->options->mu_rates[locus];
  if (kept > 0)
    {
      L->h = (double *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (double));
      jc_locus_bandwidth (ly, L->st, NULL, kept, L->h);
    }
  if (ly->ml)
    for (i = 0; i < kept; i++)
      jml_spells_prepare (L->sp + 2 * L->spo[i], jc_row_nsp (ly, L->st + i * nrow), ly->numpop);
  return kept;
}

static long
jc_read_locus (world_fmt *world, const jc_layout *ly, long locus, jc_locus *L)
{
  jc_store *js = (jc_store *) world->jointstats;
  const long nrow = ly->nrow;
  long i, spn = 0, spcap = 0;
  long *cum = NULL;
  if (jc_file_mode (world))
    return jc_read_locus_files (world, ly, locus, L);
  if (js == NULL || js->nrow != nrow || js->b[locus].n == 0)
    return 0;
  const long n = js->b[locus].n;
  /* the rows of the full combination (evenly spaced); a block uses those
     of them that belong to it, so that the blocks partition them */
  const long kfull = n < jc_maxsamples (world, nrow) ? n : jc_maxsamples (world, nrow);
  long ra, rlen, keep = 0;
  jc_part_range (n, &ra, &rlen);
  L->st = (double *) mycalloc ((size_t) (kfull * nrow), sizeof (double));
  L->src = (long *) mycalloc ((size_t) (kfull + 1), sizeof (long));
  if (ly->ml)
    {   /* offsets of the buffer rows' spells */
      cum = (long *) mycalloc ((size_t) (n + 1), sizeof (long));
      for (i = 0; i < n; i++)
        cum[i + 1] = cum[i] + js->b[locus].nspr[i];
      L->spo = (long *) mycalloc ((size_t) (kfull + 1), sizeof (long));
    }
  for (i = 0; i < kfull; i++)
    {
      const long row = (i * n) / kfull;
      if (jc_row_in_part (js->b[locus].rows + row * nrow, nrow, row, ra, rlen))
        {
          if (ly->ml)
            {
              L->spo[keep] = spn;
              jc_locus_add_spells (L, &spn, &spcap, js->b[locus].sp + 2 * cum[row], js->b[locus].nspr[row]);
            }
          L->src[keep] = row;
          memcpy (L->st + (keep++) * nrow, js->b[locus].rows + row * nrow, sizeof (double) * (size_t) nrow);
        }
    }
  myfree (cum);
  if (keep >= 1)
    {   /* smoothing bandwidths from all rows of the full combination */
      long *idx = (long *) mycalloc ((size_t) (kfull + 1), sizeof (long));
      for (i = 0; i < kfull; i++)
        idx[i] = (i * n) / kfull;
      L->h = (double *) mycalloc ((size_t) (ly->nsplit + 1), sizeof (double));
      jc_locus_bandwidth (ly, js->b[locus].rows, idx, kfull, L->h);
      myfree (idx);
    }
  if (keep < 1)
    {
      myfree (L->st);
      L->st = NULL;
      myfree (L->spo);
      L->spo = NULL;
      myfree (L->sp);
      L->sp = NULL;
      return 0;
    }
  L->n = keep;
  L->total = n;
  L->mu = world->options->mu_rates[locus];
  if (ly->ml)
    for (i = 0; i < keep; i++)
      jml_spells_prepare (L->sp + 2 * L->spo[i], jc_row_nsp (ly, L->st + i * nrow), ly->numpop);
  return keep;
}

/* Mittag-Leffler: the table at one alpha, cached by value (a proposal
   alternates between the current and the proposed alpha) */
#define JC_MLCACHE 16
static const jml_alpha *
jc_ml_table (double alpha)
{
  static jml_alpha cache[JC_MLCACHE];
  static boolean valid[JC_MLCACHE];
  static long next = 0;
  long i;
  for (i = 0; i < JC_MLCACHE; i++)
    if (valid[i] && cache[i].alpha == alpha)
      return &cache[i];
  i = next;
  next = (next + 1) % JC_MLCACHE;
  jml_alpha_set (&cache[i], alpha);
  valid[i] = TRUE;
  return &cache[i];
}

/* the spells of a kept row */
#define JC_SP(Ll, s) ((Ll)->spo != NULL ? (Ll)->sp + 2 * (Ll)->spo[s] : NULL)

/* Mittag-Leffler: log density of population pop's spells at phi */
static double
jc_ml_term (const jc_layout *ly, const double *st, const double *sp, long pop, const double *phi)
{
  double M[3] = {0.0, 0.0, 0.0};
  long q, e, from, to, off = 0;
  for (q = 0; q < pop; q++)
    off += (long) st[ly->off_nsp + q];
  for (e = mstart (pop, ly->numpop); e < mend (pop, ly->numpop); e++)
    if (ly->rep[e] >= 0)
      {
        m2mm (e, ly->numpop, &from, &to);
        M[from] = jc_value (ly, phi, e);
      }
  const double alpha = ly->apar[pop] >= 0 ? phi[ly->apar[pop]] : ly->afix[pop];
  return jml_pop_logp (sp + 2 * off, (long) st[ly->off_nsp + pop], ly->numpop, st[ly->off_mu],
                       st[ly->off_inh], jc_value (ly, phi, pop), M, jc_ml_table (alpha),
                       ly->mlscale, st[ly->off_lc]);
}

/* Mittag-Leffler: log Z(G) of a row. Populations with their own alpha
   parameter, or sharing one, are integrated over it at Gauss-Legendre
   nodes (the product of their inner integrals at each node); populations
   with a fixed alpha need the inner integral only. I: per population the
   inner integral's settings and warm start (kept from row to row). */
typedef struct
{
  world_fmt *world;
  long p;
} jc_mlprior;

static double
jc_ml_logprior (void *ctx, double x)
{
  const jc_mlprior *c = (const jc_mlprior *) ctx;
  return scaling_prior_joint (c->world, c->p, x);
}

static double
jc_ml_lognorm (world_fmt *world, const jc_layout *ly, const double *st, const double *sp, jml_inner *I)
{
  static jml_alpha *nodetab = NULL;   /* JML_NALPHA tables per alpha parameter */
  static double *nodex = NULL, *nodew = NULL;
  static long *nodep = NULL, nnodep = 0;
  const long numpop = ly->numpop;
  long pop, q, j, off[4];
  boolean done[3] = {FALSE, FALSE, FALSE};
  double lz = 0.0;
  off[0] = 0;
  for (pop = 0; pop < numpop; pop++)
    off[pop + 1] = off[pop] + (long) st[ly->off_nsp + pop];
  for (pop = 0; pop < numpop; pop++)
    {
      if (done[pop])
        continue;
      const long p = ly->apar[pop];
      if (p < 0)
        {
          lz += jml_inner_log (sp + 2 * off[pop], off[pop + 1] - off[pop], numpop, pop, st[ly->off_mu],
                               st[ly->off_inh], jc_ml_table (ly->afix[pop]), ly->mlscale, st[ly->off_lc], &I[pop]);
          done[pop] = TRUE;
          continue;
        }
      /* the node tables of parameter p (built once) */
      long slot = -1;
      for (q = 0; q < nnodep; q++)
        if (nodep[q] == p)
          slot = q;
      if (slot < 0)
        {
          slot = nnodep++;
          nodep = (long *) myrealloc (nodep, sizeof (long) * (size_t) nnodep);
          nodex = (double *) myrealloc (nodex, sizeof (double) * (size_t) (nnodep * JML_NALPHA));
          nodew = (double *) myrealloc (nodew, sizeof (double) * (size_t) (nnodep * JML_NALPHA));
          nodetab = (jml_alpha *) myrealloc (nodetab, sizeof (jml_alpha) * (size_t) (nnodep * JML_NALPHA));
          nodep[slot] = p;
          jml_alpha_nodes (world->bayes->minparam[p], world->bayes->maxparam[p],
                           nodex + slot * JML_NALPHA, nodew + slot * JML_NALPHA);
          for (j = 0; j < JML_NALPHA; j++)
            {
              memset (&nodetab[slot * JML_NALPHA + j], 0, sizeof (jml_alpha));
              jml_alpha_set (&nodetab[slot * JML_NALPHA + j], nodex[slot * JML_NALPHA + j]);
            }
        }
      double vals[JML_NALPHA];
      long order[JML_NALPHA];
      for (j = 0; j < JML_NALPHA; j++)
        order[j] = j;
      for (j = 1; j < JML_NALPHA; j++)   /* from alpha near 1 down: warm starts */
        for (q = j; q > 0 && nodex[slot * JML_NALPHA + order[q]] > nodex[slot * JML_NALPHA + order[q - 1]]; q--)
          {
            const long tmp = order[q];
            order[q] = order[q - 1];
            order[q - 1] = tmp;
          }
      for (j = 0; j < JML_NALPHA; j++)
        {
          const long k = slot * JML_NALPHA + order[j];
          double v = log (nodew[k]) + scaling_prior_joint (world, p, nodex[k]);
          for (q = pop; q < numpop; q++)
            if (ly->apar[q] == p)
              v += jml_inner_log (sp + 2 * off[q], off[q + 1] - off[q], numpop, q, st[ly->off_mu],
                                  st[ly->off_inh], &nodetab[k], ly->mlscale, st[ly->off_lc], &I[q]);
          vals[j] = v;
        }
      for (q = pop; q < numpop; q++)
        if (ly->apar[q] == p)
          done[q] = TRUE;
      lz += jc_logsumexp (vals, JML_NALPHA);
    }
  return lz;
}

/* log p(G|phi) up to terms that do not depend on any estimated parameter */
static double
jc_full_term (const jc_layout *ly, const boolean *active, const double *st, const double *sp,
              const double *phi, double mu)
{
  double v = 0.0;
  long pop, e;
  long sg;
  if (ly->ml)
    {
      for (pop = 0; pop < ly->numpop; pop++)
        v += jc_ml_term (ly, st, sp, pop, phi);
      return v;
    }
  for (pop = 0; pop < ly->numpop; pop++)
    if (ly->rep[pop] >= 0 && active[ly->rep[pop]])
      {
        const long s = ly->pop_slot[pop];
        const long j = ly->sky_j[pop];
        if (j >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            v += jc_pop_seg_term (ly, st, pop, sg, phi[ly->rep[pop]] * (sg ? phi[JC_TAU (ly, j, sg)] : 1.0),
                                  s >= 0 ? phi[ly->slot_param[s]] : 0.0, mu);
        else
          v += jc_pop_term_phi (ly, st, pop, phi, mu);
      }
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->rep[e] >= 0 && active[ly->rep[e]])
      {
        const long j = ly->sky_j[e];
        if (j >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            v += jc_mig_seg_term (ly, st, e, sg, phi[ly->rep[e]] * (sg ? phi[JC_TAU (ly, j, sg)] : 1.0), mu);
        else
          v += ly->gate_k[e] >= 0
            ? jc_gate_term (ly, st, e, jc_rate (ly, phi, e), phi[ly->split_pmu[ly->gate_k[e]]], mu)
            : jc_mig_term (ly, st, e, jc_rate (ly, phi, e), mu);
      }
  for (e = 0; e < ly->nsplit; e++)
    v += jc_div_term (ly, st, e, phi[ly->split_pmu[e]],
                      ly->split_psig[e] >= 0 ? phi[ly->split_psig[e]] : ly->split_sfix[e]);
  if (ly->rb)
    v += jc_rb_term (ly, st, phi, mu);
  if (ly->rate && active[ly->rate_p])
    v += jc_rate_term (ly, st, phi[ly->rate_p]);
  return v;
}


/* The terms of jc_full_term() one by one ("units", in its order): a
   proposal of one parameter recomputes only the units that depend on it
   and sums all units again in the same order, so the sum is the same
   number jc_full_term() gives (no old + new - old updates: a divergence
   term at JC_DOUT would swallow the others) */
enum { JC_U_POP, JC_U_POPSEG, JC_U_MIG, JC_U_MIGSEG, JC_U_SPLIT, JC_U_ML, JC_U_SPLITRB, JC_U_RATE };
typedef struct
{
  long n;
  int *kind;
  long *idx, *sg;
  long *ndep;     /* per parameter: number of dependent units */
  long **dep;     /* per parameter: the dependent units */
} jc_units;

static void
jc_units_add (jc_units *U, int kind, long idx, long sg)
{
  U->kind = (int *) myrealloc (U->kind, sizeof (int) * (size_t) (U->n + 1));
  U->idx = (long *) myrealloc (U->idx, sizeof (long) * (size_t) (U->n + 1));
  U->sg = (long *) myrealloc (U->sg, sizeof (long) * (size_t) (U->n + 1));
  U->kind[U->n] = kind;
  U->idx[U->n] = idx;
  U->sg[U->n] = sg;
  U->n++;
}

static void
jc_units_dep (jc_units *U, long nphi, long p, long u)
{
  if (p < 0 || p >= nphi)
    return;
  if (U->ndep[p] > 0 && U->dep[p][U->ndep[p] - 1] == u)
    return;
  U->dep[p] = (long *) myrealloc (U->dep[p], sizeof (long) * (size_t) (U->ndep[p] + 1));
  U->dep[p][U->ndep[p]++] = u;
}

static void
jc_units_make (const jc_layout *ly, const boolean *active, jc_units *U)
{
  long pop, e, sg, u;
  memset (U, 0, sizeof (*U));
  if (ly->ml)   /* Mittag-Leffler: one term per population */
    {
      for (pop = 0; pop < ly->numpop; pop++)
        jc_units_add (U, JC_U_ML, pop, 0);
      U->ndep = (long *) mycalloc ((size_t) ly->nphi, sizeof (long));
      U->dep = (long **) mycalloc ((size_t) ly->nphi, sizeof (long *));
      for (u = 0; u < U->n; u++)
        {
          pop = U->idx[u];
          jc_units_dep (U, ly->nphi, ly->rep[pop], u);
          for (e = mstart (pop, ly->numpop); e < mend (pop, ly->numpop); e++)
            jc_units_dep (U, ly->nphi, ly->rep[e], u);
          jc_units_dep (U, ly->nphi, ly->apar[pop], u);
        }
      return;
    }
  for (pop = 0; pop < ly->numpop; pop++)
    if (ly->rep[pop] >= 0 && active[ly->rep[pop]])
      {
        if (ly->sky_j[pop] >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            jc_units_add (U, JC_U_POPSEG, pop, sg);
        else
          jc_units_add (U, JC_U_POP, pop, 0);
      }
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->rep[e] >= 0 && active[ly->rep[e]])
      {
        if (ly->sky_j[e] >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            jc_units_add (U, JC_U_MIGSEG, e, sg);
        else
          jc_units_add (U, JC_U_MIG, e, 0);
      }
  for (e = 0; e < ly->nsplit; e++)
    jc_units_add (U, JC_U_SPLIT, e, 0);
  if (ly->rb)
    jc_units_add (U, JC_U_SPLITRB, 0, 0);
  if (ly->rate && active[ly->rate_p])
    jc_units_add (U, JC_U_RATE, 0, 0);
  U->ndep = (long *) mycalloc ((size_t) ly->nphi, sizeof (long));
  U->dep = (long **) mycalloc ((size_t) ly->nphi, sizeof (long *));
  for (u = 0; u < U->n; u++)
    {
      const long i = U->idx[u];
      switch (U->kind[u])
        {
        case JC_U_POP:
          jc_units_dep (U, ly->nphi, ly->rep[i], u);
          if (ly->pop_slot[i] >= 0)
            jc_units_dep (U, ly->nphi, ly->slot_param[ly->pop_slot[i]], u);
          if (ly->sharp_pu[i])
            {   /* evaluated on the moved row: also the split times */
              long k;
              for (k = 0; k < ly->nsplit; k++)
                if (ly->split_to[k] == i || ly->split_from[k] == i)
                  jc_units_dep (U, ly->nphi, ly->split_pmu[k], u);
            }
          break;
        case JC_U_POPSEG:
        case JC_U_MIGSEG:
          jc_units_dep (U, ly->nphi, ly->rep[i], u);
          if (U->kind[u] == JC_U_POPSEG && ly->off_sA[i] >= 0)
            jc_units_dep (U, ly->nphi, ly->slot_param[ly->pop_slot[i]], u);
          if (U->sg[u] > 0)
            jc_units_dep (U, ly->nphi, JC_TAU (ly, ly->sky_j[i], U->sg[u]), u);
          break;
        case JC_U_MIG:
          jc_units_dep (U, ly->nphi, ly->rep[i], u);
          if (ly->gate_k[i] >= 0)
            jc_units_dep (U, ly->nphi, ly->split_pmu[ly->gate_k[i]], u);
          if (ly->xnm[i])
            {
              long from, to;
              m2mm (i, ly->numpop, &from, &to);
              jc_units_dep (U, ly->nphi, ly->rep[to], u);
            }
          break;
        case JC_U_SPLIT:
          jc_units_dep (U, ly->nphi, ly->split_pmu[i], u);
          jc_units_dep (U, ly->nphi, ly->split_psig[i], u);
          break;
        case JC_U_RATE:
          jc_units_dep (U, ly->nphi, ly->rate_p, u);
          break;
        case JC_U_SPLITRB:
          {   /* any split's crossing: its D, S, the Theta of daughter and
                 ancestor, the rates into them (and their gates) */
            long k, f, from, to;
            boolean *dep = (boolean *) mycalloc ((size_t) ly->nphi, sizeof (boolean));
            for (k = 0; k < ly->nsplit; k++)
              if (ly->rb_ok[k])
                {
                  const long d = ly->split_to[k], an = ly->split_from[k];
                  if (ly->split_pmu[k] >= 0) dep[ly->split_pmu[k]] = TRUE;
                  if (ly->split_psig[k] >= 0) dep[ly->split_psig[k]] = TRUE;
                  if (ly->rep[d] >= 0) dep[ly->rep[d]] = TRUE;
                  if (ly->rep[an] >= 0) dep[ly->rep[an]] = TRUE;
                  /* growth of daughter and ancestor (sharp terms move
                     growth-weighted integrals) */
                  if (ly->pop_slot[d] >= 0) dep[ly->slot_param[ly->pop_slot[d]]] = TRUE;
                  if (ly->pop_slot[an] >= 0) dep[ly->slot_param[ly->pop_slot[an]]] = TRUE;
                  for (f = 0; f < ly->numpop2; f++)
                    {   /* skyline multipliers of d, the ancestor, the entries into them */
                      long sg;
                      if (ly->sky_j[f] < 0)
                        continue;
                      m2mm (f, ly->numpop, &from, &to);
                      if (f < ly->numpop ? (f != d && f != an) : (to != d && to != an))
                        continue;
                      for (sg = 1; sg < ly->nseg; sg++)
                        dep[JC_TAU (ly, ly->sky_j[f], sg)] = TRUE;
                    }
                  for (f = ly->numpop; f < ly->numpop2; f++)
                    {
                      m2mm (f, ly->numpop, &from, &to);
                      if (ly->mzero[f] || (to != d && to != an))
                        continue;
                      if (ly->rep[f] >= 0) dep[ly->rep[f]] = TRUE;
                      if (ly->gate_k[f] >= 0) dep[ly->split_pmu[ly->gate_k[f]]] = TRUE;
                    }
                }
            for (k = 0; k < ly->nphi; k++)
              if (dep[k])
                jc_units_dep (U, ly->nphi, k, u);
            myfree (dep);
          }
          break;
        }
    }
}

static void
jc_units_free (jc_units *U, long nphi)
{
  long p;
  for (p = 0; p < nphi && U->dep != NULL; p++)
    myfree (U->dep[p]);
  myfree (U->dep);
  myfree (U->ndep);
  myfree (U->kind);
  myfree (U->idx);
  myfree (U->sg);
  memset (U, 0, sizeof (*U));
}

/* one unit: the same expression as in jc_full_term() */
static double
jc_unit_term (const jc_layout *ly, const jc_units *U, long u, const double *st, const double *sp,
              const double *phi, double mu)
{
  const long i = U->idx[u], sg = U->sg[u];
  switch (U->kind[u])
    {
    case JC_U_ML:
      return jc_ml_term (ly, st, sp, i, phi);
    case JC_U_POP:
      return jc_pop_term_phi (ly, st, i, phi, mu);
    case JC_U_POPSEG:
      return jc_pop_seg_term (ly, st, i, sg, phi[ly->rep[i]] * (sg ? phi[JC_TAU (ly, ly->sky_j[i], sg)] : 1.0),
                              ly->pop_slot[i] >= 0 ? phi[ly->slot_param[ly->pop_slot[i]]] : 0.0, mu);
    case JC_U_MIG:
      if (ly->gate_k[i] >= 0)
        return jc_gate_term (ly, st, i, jc_rate (ly, phi, i), phi[ly->split_pmu[ly->gate_k[i]]], mu);
      return jc_mig_term (ly, st, i, jc_rate (ly, phi, i), mu);
    case JC_U_MIGSEG:
      return jc_mig_seg_term (ly, st, i, sg, phi[ly->rep[i]] * (sg ? phi[JC_TAU (ly, ly->sky_j[i], sg)] : 1.0), mu);
    case JC_U_SPLITRB:
      return jc_rb_term (ly, st, phi, mu);
    case JC_U_RATE:
      return jc_rate_term (ly, st, phi[ly->rate_p]);
    default:
      return jc_div_term (ly, st, i, phi[ly->split_pmu[i]],
                          ly->split_psig[i] >= 0 ? phi[ly->split_psig[i]] : ly->split_sfix[i]);
    }
}

/* ---------------------------------------------------- local evaluation */

/* The loci one process evaluates: all loci in a serial run, the loci this
   locus worker collected under MPI. Everything here is per process; only
   sums of log f_l travel. */
typedef struct
{
  jc_layout ly;
  boolean *active;
  long nl, maxn;
  jc_locus *L;
  double *lterm, *newlterm, *newcur;
  jc_units U;
  double *newunit;   /* nl x maxn x (units of the pending parameter) */
  long pending;      /* the parameter of the pending proposal (-2: a block proposal) */
  boolean skip;      /* skip negligible genealogies between refreshes (allowed) */
  boolean skipon;    /* ... now: only in the sampling phase of a chain at beta = 1 */
  long rint;         /* proposals between refreshes (adapted) */
  long nprop;        /* proposals since the last refresh */
  double *newall;    /* block proposal: nl x maxn x nunit, every unit */
  long naccept;      /* accepted steps since the sums were last recomputed */
  long incr;         /* more terms than this: sums updated by differences */
  long tmin, tmax;
  double lognz;   /* log of the prior mass the normalizers integrate (no-event row) */
} jc_local;

static void
jc_local_free (jc_local *J)
{
  long l;
  for (l = 0; l < J->nl; l++)
    {
      myfree (J->L[l].st);
      myfree (J->L[l].sp);
      myfree (J->L[l].spo);
      myfree (J->L[l].src);
      myfree (J->L[l].ax);
      myfree (J->L[l].cur);
      myfree (J->L[l].logz);
      myfree (J->L[l].unit);
      myfree (J->L[l].h);
    }
  jc_hcur = NULL;
  myfree (J->L);
  myfree (J->newunit);
  myfree (J->newall);
  jc_units_free (&J->U, J->ly.nphi);
  myfree (J->lterm);
  myfree (J->newlterm);
  myfree (J->newcur);
  myfree (J->active);
  jc_layout_free (&J->ly);
  memset (J, 0, sizeof (*J));
}

/* active parameters: identical on every rank (priors and map are shared) */
static boolean *
jc_active (world_fmt *world, const jc_layout *ly)
{
  boolean *active = (boolean *) mycalloc ((size_t) ly->nphi, sizeof (boolean));
  long p, s;
  for (p = ly->np; p < ly->nphi; p++)
    active[p] = TRUE;   /* skyline multipliers */
  for (p = 0; p < ly->numpop2; p++)
    if (ly->rep[p] >= 0)
      active[ly->rep[p]] = world->bayes->maxparam[ly->rep[p]] > world->bayes->minparam[ly->rep[p]];
  for (s = 0; s < ly->nslot; s++)
    active[ly->slot_param[s]] = world->bayes->maxparam[ly->slot_param[s]] > world->bayes->minparam[ly->slot_param[s]];
  for (s = 0; s < ly->numpop && ly->ml; s++)   /* Mittag-Leffler alpha */
    if (ly->apar[s] >= 0)
      active[ly->apar[s]] = world->bayes->maxparam[ly->apar[s]] > world->bayes->minparam[ly->apar[s]];
  for (s = 0; s < ly->nsplit; s++)
    {
      active[ly->split_pmu[s]] = world->bayes->maxparam[ly->split_pmu[s]] > world->bayes->minparam[ly->split_pmu[s]];
      if (ly->split_psig[s] >= 0)
        active[ly->split_psig[s]] = TRUE;
    }
  if (ly->rate)
    active[ly->rate_p] = TRUE;
  for (p = 0; p < ly->numpop; p++)
    {
      boolean aff = FALSE;
      for (s = 0; s < ly->nsplit && ly->sharp; s++)
        if (ly->split_to[s] == p || ly->split_from[s] == p)
          aff = TRUE;
      ly->sharp_pu[p] = aff && ly->rep[p] >= 0 && active[ly->rep[p]] && ly->sky_j[p] < 0;
    }
  return active;
}

/* length of phi: the parameters and the skyline multipliers */
static long
jc_phi_size (world_fmt *world)
{
  jc_layout ly;
  jc_layout_make (world, &ly);
  const long n = ly.nphi;
  jc_layout_free (&ly);
  return n;
}

/* Skipping negligible genealogies (Mittag-Leffler joint, 2026-11-01): a
   proposal evaluates only the rows whose log weight was within JC_SKIPT of
   the locus' largest at the last refresh; the others enter with their stale
   values as one sum. Every JC_SKIPR proposals all rows are recomputed at
   the current state (the refresh). The skipped weights are below e^-40 of
   the largest; MIGRATE_JC_DEBUG reports the largest change of a locus term
   that a refresh finds. Experiments: MIGRATE_JC_SKIP=0/1. */
#define JC_SKIPT 40.0
#define JC_SKIPR 20          /* proposals between refreshes at first; adapted between */
#define JC_SKIPRMIN 2        /* ... JC_SKIPRMIN and JC_SKIPRMAX: halved when a refresh */
#define JC_SKIPRMAX 40       /* changes a locus term by more than 1e-3, longer after */
                             /* changes below 1e-5 */

static boolean
jc_skip_mode (const jc_layout *ly)
{
  const char *e = getenv ("MIGRATE_JC_SKIP");
  (void) ly;
  /* off: with the ML M posteriors' heavy tails the chain sometimes moved far
     within a few proposals, and refreshes found locus terms changed by up to
     2.3 log units (2026-11-01); only as an experiment */
  return e != NULL && atol (e) != 0;
}

/* floor of a locus term log f_l: a locus whose every genealogy is impossible
   at phi (divergence terms at JC_DOUT, or minus infinity) has log f_l near
   -1e30; the chains add differences of these terms to a running sum (logF),
   and a visit there and back left round-off of 1e15 and more in it (the
   prior rungs of jc_ss, 2026-11-08: joint scaling factor 3e10 for d-only
   models). exp(JC_LTERM_FLOOR) is still zero. */
#define JC_LTERM_FLOOR (-1e8)
static double
jc_lterm_floor (double x)
{
  return (x > JC_LTERM_FLOOR) ? x : JC_LTERM_FLOOR;
}

/* log of sum exp(v[ax]) over the active rows plus the inactive sum */
static double
jc_lterm_active (const jc_locus *Ll, const double *v)
{
  double mx = Ll->sin > 0.0 ? Ll->m0 : -HUGE_VAL, s = 0.0;
  long q;
  for (q = 0; q < Ll->na; q++)
    if (v[Ll->ax[q]] > mx)
      mx = v[Ll->ax[q]];
  if (mx == -HUGE_VAL)
    return JC_LTERM_FLOOR;
  for (q = 0; q < Ll->na; q++)
    s += exp (v[Ll->ax[q]] - mx);
  if (Ll->sin > 0.0)
    s += Ll->sin * exp (Ll->m0 - mx);
  return jc_lterm_floor (mx + log (s));
}

/* the log Z cache entry of row sidx of locus L (NULL: none) */
static double *
jc_zcache (world_fmt *world, const jc_layout *ly, const jc_locus *L, long sidx)
{
  jc_store *js = jc_get_store (world, ly->nrow);
  jc_buffer *b;
  long i;
  if (L->src == NULL || L->locus < 0 || L->locus >= js->loci)
    return NULL;
  b = &js->b[L->locus];
  if (b->nzc < L->total)
    {
      b->zc = (double *) myrealloc (b->zc, sizeof (double) * (size_t) L->total);
      for (i = b->nzc; i < L->total; i++)
        b->zc[i] = NAN;
      b->nzc = L->total;
    }
  return (L->src[sidx] >= 0 && L->src[sidx] < b->nzc) ? &b->zc[L->src[sidx]] : NULL;
}

/* builds the local object for the given loci at phi; returns the local sum
   of log f_l */
static double
jc_local_setup (world_fmt *world, jc_local *J, const long *loci, long nloci, const double *phi)
{
  long l, p, s, g, i;
  memset (J, 0, sizeof (*J));
  jc_layout_make (world, &J->ly);
  const jc_layout *ly = &J->ly;
  const long numpop = ly->numpop;
  J->active = jc_active (world, ly);
  J->L = (jc_locus *) mycalloc ((size_t) (nloci > 0 ? nloci : 1), sizeof (jc_locus));
  for (i = 0; i < nloci; i++)
    {
      const long locus = loci[i];
      if (locus < 0 || locus >= world->loci || world->data->skiploci[locus])
        continue;
      J->L[J->nl].locus = locus;
      if (jc_read_locus (world, ly, locus, &J->L[J->nl]) > 0)
        J->nl++;
    }
  if (J->nl == 0)
    return 0.0;
  /* grids and log priors */
  double *tgrid = (double *) mycalloc ((size_t) (ly->numpop2 * JC_GRID), sizeof (double));
  double *tpri = (double *) mycalloc ((size_t) (ly->numpop2 * JC_GRID), sizeof (double));
  double *tgrid_s = (double *) mycalloc ((size_t) (ly->numpop2 * JC_TGRID), sizeof (double));
  double *tpri_s = (double *) mycalloc ((size_t) (ly->numpop2 * JC_TGRID), sizeof (double));
  double *gpri = (double *) mycalloc ((size_t) ((ly->nslot + 1) * JC_NG), sizeof (double));
  for (p = 0; p < ly->numpop2; p++)
    {
      if (!J->active[p])
        continue;
      const double hi = world->bayes->maxparam[p];
      const double lo = world->bayes->minparam[p] > hi * 1e-8 ? world->bayes->minparam[p] : hi * 1e-8;
      for (g = 0; g < JC_GRID; g++)
        {
          /* the ends exactly at the bounds: exp(log(x)) can fall just
             outside the prior range, which then gave the end point prior
             density zero and dropped half of the end interval (a log Z
             error up to 0.04 that depended on the genealogy) */
          tgrid[p * JC_GRID + g] = g == 0 ? lo : (g == JC_GRID - 1 ? hi
            : exp (log (lo) + (log (hi) - log (lo)) * (double) g / (JC_GRID - 1)));
          tpri[p * JC_GRID + g] = scaling_prior_joint (world, p, tgrid[p * JC_GRID + g]);
        }
      for (g = 0; g < JC_TGRID; g++)
          {
            tgrid_s[p * JC_TGRID + g] = g == 0 ? lo : (g == JC_TGRID - 1 ? hi
              : exp (log (lo) + (log (hi) - log (lo)) * (double) g / (JC_TGRID - 1)));
            tpri_s[p * JC_TGRID + g] = scaling_prior_joint (world, p, tgrid_s[p * JC_TGRID + g]);
          }
    }
  for (s = 0; s < ly->nslot; s++)
    for (g = 0; g < JC_NG; g++)
      gpri[s * JC_NG + g] = scaling_prior_joint (world, ly->slot_param[s], ly->ggrid[s * JC_NG + g]);
  jc_divquad dq;
  jc_divquad_make (world, ly, &dq);
  /* Mittag-Leffler: the function table and the inner integrals' settings */
  jml_inner mlI[3];
  jc_mlprior mlctx[3];
  memset (mlI, 0, sizeof (mlI));
  if (ly->ml)
    {
      double amin = 1.0, amax = 0.0;
      long pop, e, from, to;
      for (pop = 0; pop < numpop; pop++)
        {
          const long ap = ly->apar[pop];
          const double lo = ap >= 0 ? world->bayes->minparam[ap] : ly->afix[pop];
          const double hi = ap >= 0 ? world->bayes->maxparam[ap] : ly->afix[pop];
          if (lo < amin) amin = lo;
          if (hi > amax) amax = hi;
          mlctx[pop].world = world;
          mlctx[pop].p = ly->rep[pop];
          mlI[pop].tlo = world->bayes->minparam[ly->rep[pop]];
          mlI[pop].thi = world->bayes->maxparam[ly->rep[pop]];
          mlI[pop].logprior = jc_ml_logprior;
          mlI[pop].ctx = &mlctx[pop];
          for (e = mstart (pop, numpop); e < mend (pop, numpop); e++)
            if (ly->rep[e] >= 0)
              {
                m2mm (e, numpop, &from, &to);
                mlI[pop].mig_to[mlI[pop].r++] = from;
                mlI[pop].mmax = world->bayes->maxparam[ly->rep[e]];
              }
        }
      if (amin < 1.0)
        jml_table_init (amin, amax);
    }
  /* per-sample normalizers and current terms */
  double total = 0.0;
  J->lterm = (double *) mycalloc ((size_t) J->nl, sizeof (double));
  J->newlterm = (double *) mycalloc ((size_t) J->nl, sizeof (double));
  J->tmin = J->tmax = J->L[0].total;
  jc_units_make (ly, J->active, &J->U);
  J->pending = -1;
  /* experiments: MIGRATE_JC_INCR sets the term count above which sums are
     updated by differences */
  J->incr = getenv ("MIGRATE_JC_INCR") ? atol (getenv ("MIGRATE_JC_INCR")) : JC_INCR_UNITS;
  const long nu = J->U.n;
  for (l = 0; l < J->nl; l++)
    {
      long sidx, u, nrb = 0, nrbp = 0;
      const double t0rb = jc_now ();
      jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      Ll->cur = (double *) mycalloc ((size_t) Ll->n, sizeof (double));
      Ll->logz = (double *) mycalloc ((size_t) Ll->n, sizeof (double));
      Ll->unit = (double *) mycalloc ((size_t) (Ll->n * (nu > 0 ? nu : 1)), sizeof (double));
      for (sidx = 0; sidx < Ll->n; sidx++)
        {
          const double *st = Ll->st + sidx * ly->nrow;
          double *un = Ll->unit + sidx * nu;
          double v = 0.0;
          double *zc = jc_zcache (world, ly, Ll, sidx);
          if (zc != NULL && !isnan (*zc))
            Ll->logz[sidx] = *zc;
          else
            {
              Ll->logz[sidx] = ly->ml ? jc_ml_lognorm (world, ly, st, JC_SP (Ll, sidx), mlI)
                : (ly->rb && st[ly->off_R] >= 0.0)
                ? (ly->sharp ? jc_sharp_lognorm (world, ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq)
                   : jc_rb_lognorm (world, ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq))
                : jc_log_normalizer (ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
              if (zc != NULL)
                *zc = Ll->logz[sidx];
            }
          for (u = 0; u < nu; u++)
            {
              un[u] = jc_unit_term (ly, &J->U, u, st, JC_SP (Ll, sidx), phi, st[ly->off_mu]);
              v += un[u];
            }
          Ll->cur[sidx] = v - Ll->logz[sidx];
          if (ly->rb && st[ly->off_R] >= 0.0)
            {
              nrb++;
              nrbp += (long) st[ly->off_R + 3];
              if (getenv ("MIGRATE_JC_RBQCHECK") != NULL && nrb <= 20)
                {   /* the normalizer's quadrature: these nodes against 2x and 4x */
                  char buf[32], *old = getenv ("MIGRATE_JC_RBQ") ? strdup (getenv ("MIGRATE_JC_RBQ")) : NULL;
                  const long q0 = jc_rb_q ();
                  double z2, z4;
                  snprintf (buf, sizeof buf, "%ld", 2 * q0);
                  setenv ("MIGRATE_JC_RBQ", buf, 1);
                  z2 = jc_rb_lognorm (world, ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
                  snprintf (buf, sizeof buf, "%ld", 4 * q0);
                  setenv ("MIGRATE_JC_RBQ", buf, 1);
                  z4 = jc_rb_lognorm (world, ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
                  if (old != NULL)
                    {
                      setenv ("MIGRATE_JC_RBQ", old, 1);
                      free (old);
                    }
                  else
                    unsetenv ("MIGRATE_JC_RBQ");
                  fprintf (stderr, "JC_RBQ locus %ld row %ld pieces %ld: Q %ld %.6f, 2Q %+.2e, 4Q %+.2e\n", Ll->locus + 1, sidx,
                           (long) st[ly->off_R + 3], q0, Ll->logz[sidx], z2 - Ll->logz[sidx], z4 - Ll->logz[sidx]);
                }
            }
        }
      if (ly->rb && getenv ("MIGRATE_JC_RBDEBUG") != NULL)
        fprintf (stderr, "JC_RB locus %ld: %ld rows, %ld with an integrated crossing (%.1f pieces), setup %.1f s\n",
                 Ll->locus + 1, Ll->n, nrb, nrb ? (double) nrbp / (double) nrb : 0.0, jc_now () - t0rb);
      J->lterm[l] = jc_lterm_floor (jc_logsumexp (Ll->cur, Ll->n));
      total += J->lterm[l];
      if (Ll->n > J->maxn)
        J->maxn = Ll->n;
      if (Ll->total < J->tmin)
        J->tmin = Ll->total;
      if (Ll->total > J->tmax)
        J->tmax = Ll->total;
    }
  J->newcur = (double *) mycalloc ((size_t) (J->nl * J->maxn), sizeof (double));
  /* skipping negligible genealogies: at first every row is active */
  J->skip = jc_skip_mode (ly);
  J->rint = JC_SKIPR;
  for (l = 0; l < J->nl; l++)
    {
      long sidx;
      J->L[l].ax = (long *) mycalloc ((size_t) (J->L[l].n + 1), sizeof (long));
      for (sidx = 0; sidx < J->L[l].n; sidx++)
        J->L[l].ax[sidx] = sidx;
      J->L[l].na = J->L[l].n;
      J->L[l].sin = 0.0;
      J->L[l].m0 = 0.0;
    }
  {
    long maxdep = 1;
    for (p = 0; p < ly->nphi; p++)
      if (J->U.ndep[p] > maxdep)
        maxdep = J->U.ndep[p];
    J->newunit = (double *) mycalloc ((size_t) (J->nl * J->maxn * maxdep), sizeof (double));
  }
  {   /* the prior mass the normalizers integrate (a row without events:
         Z is then the integral of the prior over the grids and, with
         skyline=PARAM, the range rule); every f_l is relative to it */
    double *z0 = (double *) mycalloc ((size_t) ly->nrow, sizeof (double));
    long pop;
    for (pop = 0; pop < numpop; pop++)
      if (ly->pop_slot[pop] >= 0)
        {
          for (g = 0; g < JC_NG; g++)
            z0[ly->off_A[pop] + g] = -HUGE_VAL;
          if (ly->off_sA[pop] >= 0)
            for (g = 0; g < ly->nseg * JC_NG; g++)
              z0[ly->off_sA[pop] + g] = -HUGE_VAL;
        }
    z0[ly->off_mu] = 1.0;
    J->lognz = jc_log_normalizer (ly, J->active, z0, 1.0, tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
    myfree (z0);
    if (ly->ml)
      {   /* the alpha priors at the nodes of jc_ml_lognorm() */
        long q, j;
        for (p = 0; p < ly->np; p++)
          {
            boolean used = FALSE;
            for (q = 0; q < numpop; q++)
              if (ly->apar[q] == p)
                used = TRUE;
            if (!used)
              continue;
            double x[JML_NALPHA], w[JML_NALPHA], v[JML_NALPHA];
            jml_alpha_nodes (world->bayes->minparam[p], world->bayes->maxparam[p], x, w);
            for (j = 0; j < JML_NALPHA; j++)
              v[j] = log (w[j]) + scaling_prior_joint (world, p, x[j]);
            J->lognz += jc_logsumexp (v, JML_NALPHA);
          }
      }
  }
  myfree (tgrid); myfree (tpri); myfree (tgrid_s); myfree (tpri_s); myfree (gpri);
  jc_divquad_free (&dq);
  return total;
}

/* refresh: every row's terms at the current state phi, the active sets;
   adapts the refresh interval to the largest change of a locus term */
static void
jc_local_refresh (jc_local *J, const double *phi)
{
  const long nu = J->U.n;
  long l, sidx, u;
  double change = 0.0;
  for (l = 0; l < J->nl; l++)
    {
      jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      double mx = -HUGE_VAL;
      for (sidx = 0; sidx < Ll->n; sidx++)
        {
          const double *st = Ll->st + sidx * J->ly.nrow;
          double *un = Ll->unit + sidx * nu, v = 0.0;
          for (u = 0; u < nu; u++)
            {
              un[u] = jc_unit_term (&J->ly, &J->U, u, st, JC_SP (Ll, sidx), phi, st[J->ly.off_mu]);
              v += un[u];
            }
          Ll->cur[sidx] = v - Ll->logz[sidx];
          if (Ll->cur[sidx] > mx)
            mx = Ll->cur[sidx];
        }
      const double exact = jc_lterm_floor (jc_logsumexp (Ll->cur, Ll->n));
      if (isfinite (exact) && isfinite (J->lterm[l]) && fabs (exact - J->lterm[l]) > change)
        change = fabs (exact - J->lterm[l]);
      if (getenv ("MIGRATE_JC_DEBUG") != NULL && fabs (exact - J->lterm[l]) > 1e-3)
        fprintf (stderr, "JCSKIP refresh: locus term changed by %g (active %ld of %ld, interval %ld)\n",
                 fabs (exact - J->lterm[l]), Ll->na, Ll->n, J->rint);
      J->lterm[l] = exact;
      Ll->na = 0;
      Ll->sin = 0.0;
      Ll->m0 = mx;
      for (sidx = 0; sidx < Ll->n; sidx++)
        if (Ll->cur[sidx] >= mx - JC_SKIPT)
          Ll->ax[Ll->na++] = sidx;
        else if (Ll->cur[sidx] > -HUGE_VAL)
          Ll->sin += exp (Ll->cur[sidx] - mx);
    }
  if (J->skipon)
    {
      if (change > 1e-3)
        J->rint = J->rint / 2 > JC_SKIPRMIN ? J->rint / 2 : JC_SKIPRMIN;
      else if (change < 1e-5)
        J->rint = J->rint + 2 < JC_SKIPRMAX ? J->rint + 2 : JC_SKIPRMAX;
    }
}

/* skipping on or off (the chain's sampling phase at beta = 1); either way
   the rows are refreshed at phi first, and off makes every row active */
static void
jc_local_skip (jc_local *J, const double *phi, boolean on)
{
  long l, s;
  if (!J->skip || J->nl == 0 || on == J->skipon)
    return;
  jc_local_refresh (J, phi);
  J->nprop = 0;
  J->skipon = on;
  if (!on)
    for (l = 0; l < J->nl; l++)
      {
        for (s = 0; s < J->L[l].n; s++)
          J->L[l].ax[s] = s;
        J->L[l].na = J->L[l].n;
        J->L[l].sin = 0.0;
      }
}

/* proposal phi[p] -> x: local change of sum log f_l (state kept pending) */
static double
jc_local_propose (jc_local *J, double *phi, long p, double x)
{
  double delta = 0.0;
  if (J->nl == 0)   /* a worker without genealogies (its loci have none) */
    {
      J->pending = p;
      return 0.0;
    }
  if (J->skipon && ++J->nprop >= J->rint)
    {
      jc_local_refresh (J, phi);
      J->nprop = 0;
    }
  const long nu = J->U.n;
  const long nd = (p >= 0 && p < J->ly.nphi) ? J->U.ndep[p] : 0;
  const long *dep = nd > 0 ? J->U.dep[p] : NULL;
  long l, sidx, k, u;
  long *map = (long *) mycalloc ((size_t) (nu > 0 ? nu : 1), sizeof (long));
  for (u = 0; u < nu; u++)
    map[u] = -1;
  for (k = 0; k < nd; k++)
    map[dep[k]] = k;
  J->pending = p;
  for (l = 0; l < J->nl; l++)
    {
      double *nc = J->newcur + l * J->maxn;
      double *nw = J->newunit + l * J->maxn * (nd > 0 ? nd : 1);
      jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      const double keep = phi[p];
      long q;
      phi[p] = x;
      for (q = 0; q < Ll->na; q++)
        {
          sidx = Ll->ax[q];
          const double *st = Ll->st + sidx * J->ly.nrow;
          const double *un = Ll->unit + sidx * nu;
          double *nws = nw + sidx * nd;
          double v = 0.0;
          for (k = 0; k < nd; k++)
            nws[k] = jc_unit_term (&J->ly, &J->U, dep[k], st, JC_SP (Ll, sidx), phi, st[J->ly.off_mu]);
          if (nu > J->incr)
            {   /* many terms: the difference, unless a term is (near) a
                   sentinel or minus infinity, where it would swallow the rest */
              boolean ok = Ll->cur[sidx] > 0.1 * JC_DOUT;
              for (k = 0; k < nd && ok; k++)
                ok = nws[k] > 0.1 * JC_DOUT && un[dep[k]] > 0.1 * JC_DOUT;
              if (ok)
                {
                  v = Ll->cur[sidx];
                  for (k = 0; k < nd; k++)
                    v += nws[k] - un[dep[k]];
                  nc[sidx] = v;
                  continue;
                }
            }
          for (u = 0; u < nu; u++)
            v += map[u] >= 0 ? nws[map[u]] : un[u];
          nc[sidx] = v - Ll->logz[sidx];
        }
      phi[p] = keep;
      J->newlterm[l] = jc_lterm_active (Ll, nc);
      delta += J->newlterm[l] - J->lterm[l];
    }
  myfree (map);
  return delta;
}

/* block proposal: phi -> phin in every coordinate at once (all units are
   recomputed); returns the local change of sum log f_l */
static double
jc_local_propose_block (jc_local *J, const double *phi, const double *phin)
{
  double delta = 0.0;
  const long nu = J->U.n;
  long l, sidx, u, q;
  J->pending = -2;
  if (J->nl == 0)
    return 0.0;
  if (J->skipon && ++J->nprop >= J->rint)
    {
      jc_local_refresh (J, phi);
      J->nprop = 0;
    }
  if (J->newall == NULL)
    J->newall = (double *) mycalloc ((size_t) (J->nl * J->maxn * (nu > 0 ? nu : 1)), sizeof (double));
  for (l = 0; l < J->nl; l++)
    {
      jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      double *nc = J->newcur + l * J->maxn;
      double *nw = J->newall + l * J->maxn * nu;
      for (q = 0; q < Ll->na; q++)
        {
          sidx = Ll->ax[q];
          const double *st = Ll->st + sidx * J->ly.nrow;
          double v = 0.0;
          for (u = 0; u < nu; u++)
            {
              nw[sidx * nu + u] = jc_unit_term (&J->ly, &J->U, u, st, JC_SP (Ll, sidx), phin, st[J->ly.off_mu]);
              v += nw[sidx * nu + u];
            }
          nc[sidx] = v - Ll->logz[sidx];
        }
      J->newlterm[l] = jc_lterm_active (Ll, nc);
      delta += J->newlterm[l] - J->lterm[l];
    }
  return delta;
}

static void
jc_local_accept (jc_local *J)
{
  const long p = J->pending;
  if (p == -2 && J->nl > 0)
    {   /* a block proposal: every unit */
      const long nu = J->U.n;
      long l;
      for (l = 0; l < J->nl; l++)
        {
          jc_locus *Ll = &J->L[l];
          jc_hcur = Ll->h;
          long q;
          for (q = 0; q < Ll->na; q++)   /* the active rows (all without skipping) */
            {
              const long s = Ll->ax[q];
              Ll->cur[s] = J->newcur[l * J->maxn + s];
              memcpy (Ll->unit + s * nu, J->newall + (l * J->maxn + s) * nu, sizeof (double) * (size_t) nu);
            }
          J->lterm[l] = J->newlterm[l];
        }
      J->pending = -1;
      return;
    }
  if (J->nl == 0)
    {
      J->pending = -1;
      return;
    }
  const long nu = J->U.n;
  const long nd = (p >= 0 && p < J->ly.nphi) ? J->U.ndep[p] : 0;
  long l, sidx, k;
  for (l = 0; l < J->nl; l++)
    {
      jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      const double *nw = J->newunit + l * J->maxn * (nd > 0 ? nd : 1);
      long q;
      J->lterm[l] = J->newlterm[l];
      for (q = 0; q < Ll->na; q++)   /* the active rows (all without skipping) */
        {
          sidx = Ll->ax[q];
          Ll->cur[sidx] = J->newcur[l * J->maxn + sidx];
          for (k = 0; k < nd; k++)
            Ll->unit[sidx * nu + J->U.dep[p][k]] = nw[sidx * nd + k];
        }
    }
  J->pending = -1;
  if (nu > J->incr && ++J->naccept >= JC_REFRESH)
    {   /* the sums from the terms again, so that rounding cannot pile up */
      long u;
      J->naccept = 0;
      for (l = 0; l < J->nl; l++)
        {
          jc_locus *Ll = &J->L[l];
          jc_hcur = Ll->h;
          for (sidx = 0; sidx < Ll->n; sidx++)
            {
              double v = 0.0;
              for (u = 0; u < nu; u++)
                v += Ll->unit[sidx * nu + u];
              Ll->cur[sidx] = v - Ll->logz[sidx];
            }
          J->lterm[l] = jc_lterm_floor (jc_logsumexp (Ll->cur, Ll->n));
        }
    }
}

/* ------------------------------------------------- bootstrap guard */

/* The genealogy-sampling noise of the joint estimates: every locus' term
   f_l is a mean over its stored genealogies, and with few genealogies that
   fit a region of the parameter space (small importance-sampling ESS) the
   product over loci follows noise. A block bootstrap over each locus'
   genealogies, evaluated at points of the joint trace, gives for every
   replicate b the change D_bk of sum_l log f_l at point k; reweighting the
   trace points by exp(D_bk) gives each parameter's bootstrap error
   (jc_boot_errors()). The bootstrap of a locus is seeded by its locus
   number, so the result does not depend on how loci are distributed. */
#define JC_BOOTK 200         /* trace points */
#define JC_BOOTB 100         /* bootstrap replicates */
#define JC_BOOTBLOCK 10      /* consecutive genealogies per bootstrap block */

/* splitmix64: a small deterministic generator for the bootstrap */
static double
jc_boot_unif (unsigned long long *x)
{
  unsigned long long z = (*x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= z >> 31;
  return (double) (z >> 11) / 9007199254740992.0;
}

/* adds to D (B x K) the bootstrap changes of this process' loci at the K
   points pts (K x nphi) */
static void
jc_local_boot (jc_local *J, const double *pts, long K, long B, double *D)
{
  const jc_layout *ly = &J->ly;
  long l, s, k, b;
  for (l = 0; l < J->nl; l++)
    {
      const jc_locus *Ll = &J->L[l];
      jc_hcur = Ll->h;
      const long n = Ll->n;
      const long bl = n < JC_BOOTBLOCK ? 1 : JC_BOOTBLOCK;
      const long nb = (n + bl - 1) / bl;
      double *t = (double *) mycalloc ((size_t) (n * K), sizeof (double));
      double *orig = (double *) mycalloc ((size_t) K, sizeof (double));
      double *col = (double *) mycalloc ((size_t) n, sizeof (double));
      long *cnt = (long *) mycalloc ((size_t) n, sizeof (long));
      for (s = 0; s < n; s++)
        {
          const double *st = Ll->st + s * ly->nrow;
          for (k = 0; k < K; k++)
            t[s * K + k] = jc_full_term (ly, J->active, st, JC_SP (Ll, s), pts + k * ly->nphi, st[ly->off_mu]) - Ll->logz[s];
        }
      for (k = 0; k < K; k++)
        {
          for (s = 0; s < n; s++)
            col[s] = t[s * K + k];
          orig[k] = jc_lterm_floor (jc_logsumexp (col, n) - log ((double) n));
        }
      for (b = 0; b < B; b++)
        {
          unsigned long long seed = (unsigned long long) (Ll->locus + 1) * 1000003ULL + (unsigned long long) b * 7919ULL;
          long total = 0, i;
          memset (cnt, 0, sizeof (long) * (size_t) n);
          for (i = 0; i < nb; i++)
            {
              long a = (long) (jc_boot_unif (&seed) * (double) (n - bl + 1));
              long j;
              if (a > n - bl)
                a = n - bl;
              for (j = 0; j < bl; j++)
                cnt[a + j]++;
              total += bl;
            }
          for (k = 0; k < K; k++)
            {
              long m = 0;
              for (s = 0; s < n; s++)
                if (cnt[s] > 0)
                  col[m++] = log ((double) cnt[s]) + t[s * K + k];
              /* the floor as in the chains: a replicate without a genealogy
                 compatible with pts (sharp splits: no window holds D) */
              D[b * K + k] += jc_lterm_floor (jc_logsumexp (col, m) - log ((double) total)) - orig[k];
            }
        }
      myfree (t);
      myfree (orig);
      myfree (col);
      myfree (cnt);
    }
}

/* ------------------------------------------------------- MPI transport */

#ifdef MPI
#define JC_OP_STEP 1.0
#define JC_OP_END 2.0
#define JC_OP_BOOT 3.0
#define JC_OP_BLOCK 4.0
#define JC_OP_SKIP 5.0

/* master and the locus workers that answer result requests */
static long
jc_nworkers (world_fmt *world)
{
  long n = world->loci + 1;
  if (n > numcpu)
    n = numcpu;
  return n - 1;
}

static MPI_Comm
jc_make_comm (long nworkers)
{
  MPI_Group all, grp;
  MPI_Comm comm;
  int *ranks = (int *) mycalloc ((size_t) (nworkers + 1), sizeof (int));
  long i;
  for (i = 0; i <= nworkers; i++)
    ranks[i] = (int) i;
  MPI_Comm_group (comm_world, &all);
  MPI_Group_incl (all, (int) (nworkers + 1), ranks, &grp);
  MPI_Comm_create_group (comm_world, grp, 4711, &comm);
  MPI_Group_free (&grp);
  MPI_Group_free (&all);
  myfree (ranks);
  return comm;
}

/// Worker side (called from mpi_worker_event_loop() on MIGMPI_JC): evaluates
/// this rank's loci for the master's joint sampler until it ends.
void
jc_worker_service (world_fmt *world, long part)
{
  jc_set_part (world, part);
  const long np = jc_phi_size (world);
  MPI_Comm comm = jc_make_comm (jc_nworkers (world));
  double *phi = (double *) mycalloc ((size_t) np, sizeof (double));
  double pkt[4], vals[5], dummy[5];
  long *loci = (long *) mycalloc ((size_t) (locidone > 0 ? locidone : 1), sizeof (long));
  long i, pending_p = -1;
  double pending_x = 0.0;
  double *pending_phi = (double *) mycalloc ((size_t) np, sizeof (double));
  jc_local J;
  for (i = 0; i < locidone; i++)
    loci[i] = world->who[i];
  MPI_Bcast (phi, (int) np, MPI_DOUBLE, 0, comm);
  vals[0] = jc_local_setup (world, &J, loci, locidone, phi);
  vals[1] = (double) J.nl;
  vals[2] = vals[3] = 0.0;
  for (i = 0; i < J.nl; i++)
    {
      vals[2] += (double) J.L[i].n;
      vals[3] += log ((double) J.L[i].n);
    }
  vals[4] = J.nl > 0 ? (double) J.nl * J.lognz : 0.0;
  MPI_Reduce (vals, dummy, 5, MPI_DOUBLE, MPI_SUM, 0, comm);
  vals[0] = J.nl ? (double) J.tmin : HUGE_VAL;
  vals[1] = J.nl ? -(double) J.tmax : HUGE_VAL;
  MPI_Reduce (vals, dummy, 2, MPI_DOUBLE, MPI_MIN, 0, comm);
  for (;;)
    {
      /* pkt: op, accept the previous proposal, parameter, value */
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, comm);
      if (pkt[1] > 0.5 && pending_p >= 0)
        {
          jc_local_accept (&J);
          phi[pending_p] = pending_x;
        }
      else if (pkt[1] > 0.5 && pending_p == -2)
        {
          jc_local_accept (&J);
          memcpy (phi, pending_phi, sizeof (double) * (size_t) np);
        }
      pending_p = -1;
      if (pkt[0] == JC_OP_END)
        break;
      if (pkt[0] == JC_OP_BOOT)
        {   /* bootstrap guard: K points, B replicates */
          const long K = (long) pkt[2], B = (long) pkt[3];
          double *pts = (double *) mycalloc ((size_t) (K * np), sizeof (double));
          double *D = (double *) mycalloc ((size_t) (B * K), sizeof (double));
          MPI_Bcast (pts, (int) (K * np), MPI_DOUBLE, 0, comm);
          jc_local_boot (&J, pts, K, B, D);
          MPI_Reduce (D, NULL, (int) (B * K), MPI_DOUBLE, MPI_SUM, 0, comm);
          myfree (pts);
          myfree (D);
          continue;
        }
      if (pkt[0] == JC_OP_SKIP)
        {   /* skipping negligible genealogies on (pkt[2] = 1) or off */
          jc_local_skip (&J, phi, pkt[2] > 0.5);
          continue;
        }
      if (pkt[0] == JC_OP_BLOCK)
        {   /* every coordinate at once */
          MPI_Bcast (pending_phi, (int) np, MPI_DOUBLE, 0, comm);
          pending_p = -2;
          double d = jc_local_propose_block (&J, phi, pending_phi), dd = 0.0;
          MPI_Reduce (&d, &dd, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
          continue;
        }
      pending_p = (long) pkt[2];
      pending_x = pkt[3];
      double d = jc_local_propose (&J, phi, pending_p, pending_x), dd = 0.0;
      MPI_Reduce (&d, &dd, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
  jc_local_free (&J);
  myfree (pending_phi);
  myfree (phi);
  myfree (loci);
  MPI_Comm_free (&comm);
}
#endif

/* ------------------------------------------------------- master sampler */

/* the evaluator the master talks to: a local object (serial) or the locus
   workers (MPI) */
typedef struct
{
  boolean mpi;
  jc_local local;
#ifdef MPI
  MPI_Comm comm;
#endif
  long nl, tmin, tmax;
  double rows;   /* genealogies used, summed over loci */
  double logF;   /* sum_l log f_l at the current phi (f_l: the mean over genealogies) */
  double lognz;  /* sum over loci of the log prior mass of the normalizers */
  boolean accept_prev;
} jc_eval;

static void
jc_eval_start (world_fmt *world, jc_eval *E, const double *phi, long part)
{
  memset (E, 0, sizeof (*E));
  jc_set_part (world, part);
#ifdef MPI
  if (numcpu > 1)
    {
      const long nw = jc_nworkers (world);
      const long numelem2 = 2 * (world->numpop2 + (world->options->gamma ? 1 : 0));
      MYREAL *temp = (MYREAL *) mycalloc ((size_t) (numelem2 + 2), sizeof (MYREAL));
      double vals[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 }, sums[5], mins[2];
      long w;
      temp[0] = (MYREAL) MIGMPI_JC;
      temp[1] = (MYREAL) part;   /* which genealogies (jc_part) */
      for (w = 1; w <= nw; w++)
        MYMPISEND (temp, numelem2 + 2, mpisizeof, (MYINT) w, (MYINT) w, comm_world);
      myfree (temp);
      E->mpi = TRUE;
      E->comm = jc_make_comm (nw);
      MPI_Bcast ((void *) phi, (int) jc_phi_size (world), MPI_DOUBLE, 0, E->comm);
      MPI_Reduce (vals, sums, 5, MPI_DOUBLE, MPI_SUM, 0, E->comm);
      E->rows = sums[2];
      E->logF = sums[0] - sums[3];
      E->lognz = sums[4];
      vals[0] = vals[1] = HUGE_VAL;
      MPI_Reduce (vals, mins, 2, MPI_DOUBLE, MPI_MIN, 0, E->comm);
      E->nl = (long) (sums[1] + 0.5);
      E->tmin = (long) mins[0];
      E->tmax = (long) (-mins[1]);
      return;
    }
#endif
  {
    long *loci = (long *) mycalloc ((size_t) world->loci, sizeof (long));
    long l;
    for (l = 0; l < world->loci; l++)
      loci[l] = l;
    E->logF = jc_local_setup (world, &E->local, loci, world->loci, phi);
    myfree (loci);
    E->nl = E->local.nl;
    for (l = 0; l < E->local.nl; l++)
      {
        E->rows += (double) E->local.L[l].n;
        E->logF -= log ((double) E->local.L[l].n);
      }
    E->lognz = E->local.nl > 0 ? (double) E->local.nl * E->local.lognz : 0.0;
    E->tmin = E->local.tmin;
    E->tmax = E->local.tmax;
    if (getenv ("MIGRATE_JC_PROFILE") != NULL)
      {   /* diagnostics: sweep parameter MIGRATE_JC_PROFILE over its range (others
             at phi) and print per locus log f_l + log Z0 (0 under NODATA) and the
             importance-sampling ESS of its genealogies; then stop */
        jc_local *J = &E->local;
        const long p = atol (getenv ("MIGRATE_JC_PROFILE"));
        const long K = 41;
        double *ph = (double *) mycalloc ((size_t) J->ly.nphi, sizeof (double));
        long k, q;
        memcpy (ph, phi, sizeof (double) * (size_t) J->ly.nphi);
        fprintf (stderr, "JCPROFILE phi0:");
        for (q = 0; q < J->ly.nphi; q++)
          fprintf (stderr, " %g", phi[q]);
        fprintf (stderr, "\nJCPROFILE param %ld lognz %g\n", p, J->lognz);
        const double hi = world->bayes->maxparam[p];
        const double lo = world->bayes->minparam[p] > hi * 1e-4 ? world->bayes->minparam[p] : hi * 1e-4;
        for (k = 0; k < K; k++)
          {
            const double x = exp (log (lo) + (log (hi) - log (lo)) * (double) k / (K - 1));
            double sum = 0.0;
            jc_local_propose (J, ph, p, x);
            fprintf (stderr, "JCPROFILE x %.6g", x);
            for (l = 0; l < J->nl; l++)
              {
                const double *nc = J->newcur + l * J->maxn;
                const long n = J->L[l].n;
                double mx = -HUGE_VAL, s1 = 0.0, s2 = 0.0;
                long s;
                for (s = 0; s < n; s++)
                  if (nc[s] > mx)
                    mx = nc[s];
                for (s = 0; s < n; s++)
                  {
                    const double w = exp (nc[s] - mx);
                    s1 += w;
                    s2 += w * w;
                  }
                const double lf = J->newlterm[l] - log ((double) n) + J->lognz;
                sum += lf;
                fprintf (stderr, " %.3f:%.0f", lf, s2 > 0.0 ? s1 * s1 / s2 : 0.0);
              }
            fprintf (stderr, " SUM %.3f\n", sum);
            J->pending = -1;
          }
        myfree (ph);
        exit (0);
      }
  }
}

/* sends the previous decision and proposes phi[p] -> x; returns the change of
   sum_l log f_l */
static double
jc_eval_step (jc_eval *E, double *phi, boolean accept_prev, long p, double x)
{
#ifdef MPI
  if (E->mpi)
    {
      double pkt[4] = { JC_OP_STEP, accept_prev ? 1.0 : 0.0, (double) p, x }, d = 0.0, dd = 0.0;
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, E->comm);
      MPI_Reduce (&d, &dd, 1, MPI_DOUBLE, MPI_SUM, 0, E->comm);
      return dd;
    }
#endif
  if (accept_prev)
    jc_local_accept (&E->local);
  return jc_local_propose (&E->local, phi, p, x);
}

/* skipping on or off on the evaluator (a pending decision is settled first) */
static void
jc_eval_skip (jc_eval *E, const double *phi, boolean on)
{
#ifdef MPI
  if (E->mpi)
    {
      double pkt[4] = { JC_OP_SKIP, E->accept_prev ? 1.0 : 0.0, on ? 1.0 : 0.0, 0.0 };
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, E->comm);
      E->accept_prev = FALSE;
      return;
    }
#endif
  if (E->accept_prev)
    jc_local_accept (&E->local);
  E->accept_prev = FALSE;
  jc_local_skip (&E->local, phi, on);
}

/* sends the previous decision and proposes phi -> phin in every coordinate;
   returns the change of sum_l log f_l */
static double
jc_eval_block (jc_eval *E, boolean accept_prev, const double *phi, const double *phin, long nphi)
{
#ifdef MPI
  if (E->mpi)
    {
      double pkt[4] = { JC_OP_BLOCK, accept_prev ? 1.0 : 0.0, 0.0, 0.0 }, d = 0.0, dd = 0.0;
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, E->comm);
      MPI_Bcast ((void *) phin, (int) nphi, MPI_DOUBLE, 0, E->comm);
      MPI_Reduce (&d, &dd, 1, MPI_DOUBLE, MPI_SUM, 0, E->comm);
      return dd;
    }
#endif
  (void) nphi;
  (void) phi;
  if (accept_prev)
    jc_local_accept (&E->local);
  return jc_local_propose_block (&E->local, phi, phin);
}

/* bootstrap guard on the evaluator: D (B x K, zeroed) gets the summed
   bootstrap changes at the K points pts (K x nphi); a pending proposal is
   settled first */
static void
jc_eval_boot (jc_eval *E, long nphi, const double *pts, long K, long B, double *D)
{
#ifdef MPI
  if (E->mpi)
    {
      double pkt[4] = { JC_OP_BOOT, E->accept_prev ? 1.0 : 0.0, (double) K, (double) B };
      double *zero = (double *) mycalloc ((size_t) (B * K), sizeof (double));
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, E->comm);
      MPI_Bcast ((void *) pts, (int) (K * nphi), MPI_DOUBLE, 0, E->comm);
      MPI_Reduce (zero, D, (int) (B * K), MPI_DOUBLE, MPI_SUM, 0, E->comm);
      myfree (zero);
      E->accept_prev = FALSE;
      return;
    }
#endif
  (void) nphi;
  if (E->accept_prev)
    jc_local_accept (&E->local);
  E->accept_prev = FALSE;
  jc_local_boot (&E->local, pts, K, B, D);
}

static void
jc_eval_end (jc_eval *E, boolean accept_prev)
{
#ifdef MPI
  if (E->mpi)
    {
      double pkt[4] = { JC_OP_END, accept_prev ? 1.0 : 0.0, 0.0, 0.0 };
      MPI_Bcast (pkt, 4, MPI_DOUBLE, 0, E->comm);
      MPI_Comm_free (&E->comm);
      return;
    }
#endif
  if (accept_prev)
    jc_local_accept (&E->local);
  jc_local_free (&E->local);
}

static int
jc_cmp_double (const void *a, const void *b)
{
  const double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

/* sweeps of the single-site Metropolis sampler of prior(phi) prod_l
   f_l(phi)^beta on an evaluator that is set up; burn sweeps first (they
   tune step). Records phi (trace, nphi x sweeps) and sum_l log f_l (logftr)
   when given. Theta and M move multiplicatively, growth additively. */
/* Block proposals (Mittag-Leffler joint, 2026-11-01): after a single-site
   burn-in, a multivariate normal step in z (log of the positive parameters,
   the others as they are) with the covariance of the burn-in's second half,
   scaled towards 23% acceptance; the single-site sampler mixed slowly
   (autocorrelation 25-70 sweeps on three ML loci). The covariance is kept
   for the later chains of the run (stepping stones, error blocks). */
static double *jc_bL = NULL;     /* Cholesky factor, d x d */
static long jc_bd = 0;
static double jc_bscale = 0.0;

/* option joint-block (AUTO: Mittag-Leffler joints, where Theta and alpha
   are correlated; on two-population Theta/M data blocks gained only about
   10% ESS per second); experiments: MIGRATE_JC_BLOCK=0/1. Not with skyline
   multipliers (their range rule is not in the block step). */
static boolean
jc_block_mode (world_fmt *world, const jc_layout *ly)
{
  const char *e = getenv ("MIGRATE_JC_BLOCK");
  if (ly->nfree > 0)
    return FALSE;
  if (e != NULL)
    return atol (e) != 0;
  if (world->options->joint_block == JOINT_BLOCK_YES)
    return TRUE;
  if (world->options->joint_block == JOINT_BLOCK_NO)
    return FALSE;
  return ly->ml;
}

static double
jc_gauss (void)
{
  double u1 = UNIF_RANDUM ();
  while (u1 <= 1e-300)
    u1 = UNIF_RANDUM ();
  return sqrt (-2.0 * log (u1)) * cos (2.0 * PI * UNIF_RANDUM ());
}

/* covariance of the recorded z (n x d) -> jc_bL */
static void
jc_block_covariance (const double *z, long n, long d)
{
  double *m = (double *) mycalloc ((size_t) d, sizeof (double));
  double *C = (double *) mycalloc ((size_t) (d * d), sizeof (double));
  long i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++)
      m[j] += z[i * d + j] / (double) n;
  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++)
      for (k = 0; k <= j; k++)
        C[j * d + k] += (z[i * d + j] - m[j]) * (z[i * d + k] - m[k]) / (double) (n - 1);
  for (j = 0; j < d; j++)
    C[j * d + j] += 1e-10 + 1e-6 * C[j * d + j];
  myfree (jc_bL);
  jc_bL = (double *) mycalloc ((size_t) (d * d), sizeof (double));
  for (j = 0; j < d; j++)   /* Cholesky */
    for (k = 0; k <= j; k++)
      {
        double s = C[j * d + k];
        for (i = 0; i < k; i++)
          s -= jc_bL[j * d + i] * jc_bL[k * d + i];
        if (j == k)
          jc_bL[j * d + j] = sqrt (s > 1e-300 ? s : 1e-300);
        else
          jc_bL[j * d + k] = s / jc_bL[k * d + k];
      }
  jc_bd = d;
  jc_bscale = 2.38 / sqrt ((double) d);
  myfree (m);
  myfree (C);
}

/* effective sample size of x[0..n) (initial positive sequence of
   autocorrelation pairs). Not inlined: inlined into jc_chain() (Apple clang,
   -O3 -flto=thin, 2026-11-01) the autocorrelation loop was dropped and every
   ESS came out as n; the cause was not found. */
static __attribute__ ((noinline)) double
jc_ess (const double *x, long n)
{
  double m = 0.0, v = 0.0, s = 0.0;
  long i, k;
  for (i = 0; i < n; i++)
    m += x[i] / (double) n;
  for (i = 0; i < n; i++)
    v += (x[i] - m) * (x[i] - m);
  if (!(v > 0.0))
    return (double) n;
  for (k = 1; k + 1 < n / 2; k += 2)
    {
      double c1 = 0.0, c2 = 0.0;
      for (i = 0; i + k + 1 < n; i++)
        {
          c1 += (x[i] - m) * (x[i + k] - m);
          c2 += (x[i] - m) * (x[i + k + 1] - m);
        }
      if (c1 + c2 < 0.0)
        break;
      s += (c1 + c2) / v;
    }
  return (double) n / (1.0 + 2.0 * s);
}

/* the main run stops when every active coordinate has this ESS (option
   joint-ess, default JOINT_ESS_DEFAULT; experiments: MIGRATE_JC_ESSSTOP;
   checked every JC_ESSCHECK sweeps, at most the full run) */
#define JC_ESSCHECK 5000
static long
jc_essstop (world_fmt *world)
{
  if (getenv ("MIGRATE_JC_ESSSTOP") != NULL)
    return atol (getenv ("MIGRATE_JC_ESSSTOP"));
  return world->options->joint_ess;
}


/* the shared rate: the joint keeps r inside the union of the per-locus 99%
   ranges (set by jc_combine(); 0: not set) */
static double jc_rate_lo = 0.0, jc_rate_hi = 0.0;

/* the range of coordinate p in the joint chain (a parameter) */
static void
jc_bounds (world_fmt *world, const jc_layout *ly, long p, double *lo, double *hi)
{
  *lo = world->bayes->minparam[p];
  *hi = world->bayes->maxparam[p];
  if (ly->rate && p == ly->rate_p && jc_rate_hi > jc_rate_lo)
    {
      *lo = jc_rate_lo;
      *hi = jc_rate_hi;
    }
}

/* returns the sweeps of the main run (fewer than sweeps when the ESS stop
   ended it; the trace is then compacted to that length) */
/* splitstd-cv: sum over the split groups of log p(S|D) at phi (-HUGE_VAL
   outside the c range); the S coordinates then have no prior of their own
   in jc_target_prior() */
static double
jc_cv_logprior (world_fmt *world, const jc_layout *ly, const double *phi)
{
  double v = 0.0;
  long k;
  if (!world->options->splitstd_cv)
    return 0.0;
  for (k = 0; k < ly->nsplit; k++)
    {
      const long ps = ly->split_psig[k];
      if (ly->split_rep[k] != k || ps < 0 || !splitcv_is_s (world, ps, NULL))
        continue;
      const double D = ly->mhi[k] > ly->mlo[k] ? phi[ly->split_pmu[k]] : ly->mlo[k];
      const MYREAL lp = splitcv_logprior (world, ps, phi[ps], D);
      if (!(lp > (double) -HUGE))
        return -HUGE_VAL;
      v += lp;
    }
  return v;
}

/* splitstd-cv: a start with S outside [cmin D, cmax D] moves S to cmean D */
static void
jc_cv_repair (world_fmt *world, const jc_layout *ly, double *phi)
{
  long k;
  if (!world->options->splitstd_cv)
    return;
  for (k = 0; k < ly->nsplit; k++)
    {
      const long ps = ly->split_psig[k];
      MYREAL cmin, cmax, cmean;
      if (ps < 0 || !splitcv_is_s (world, ps, NULL))
        continue;
      const double D = ly->mhi[k] > ly->mlo[k] ? phi[ly->split_pmu[k]] : ly->mlo[k];
      splitcv_cbounds (ps, &cmin, &cmax, &cmean);
      if (D > 0.0 && !(phi[ps] > cmin * D && phi[ps] <= cmax * D))
        phi[ps] = cmean * D;
    }
}

/* the prior of coordinate p in the joint sampler (splitstd-cv: S in jc_cv_logprior()) */
static double
jc_target_prior (world_fmt *world, long p, double x)
{
  if (splitcv_is_s (world, p, NULL))
    return 0.0;
  return scaling_prior_joint (world, p, x);
}

static long
jc_chain (world_fmt *world, const jc_layout *ly, const boolean *active, jc_eval *E, double *phi,
          double *step, double beta, long burn, long sweeps, double *trace, double *logftr)
{
  const long essstop = (trace != NULL) ? jc_essstop (world) : 0;
  long done = sweeps;
  const long np = ly->np, nphi = ly->nphi;
  const long npos = world->numparamcumvec[SPLITSTDPRIOR];
  long *acc = (long *) mycalloc ((size_t) nphi, sizeof (long));
  long *tri = (long *) mycalloc ((size_t) nphi, sizeof (long));
  long p, sweep;
  double skyprior = jc_sky_logprior (ly, phi);
  const boolean cvon = world->options->splitstd_cv;
  double cvprior = jc_cv_logprior (world, ly, phi);
  /* block proposals: the active coordinates */
  const boolean block = jc_block_mode (world, ly);
  long d = 0, q, k, nrec = 0, bacc = 0, btri = 0, nb = 1;
  long *ai = (long *) mycalloc ((size_t) nphi, sizeof (long));
  double *zrec = NULL, *z = NULL, *zn = NULL, *phin = NULL, *eps = NULL;
  for (p = 0; p < nphi; p++)
    if (active[p])
      ai[d++] = p;
  if (block)
    {
      z = (double *) mycalloc ((size_t) (d + 1), sizeof (double));
      zn = (double *) mycalloc ((size_t) (d + 1), sizeof (double));
      eps = (double *) mycalloc ((size_t) (d + 1), sizeof (double));
      phin = (double *) mycalloc ((size_t) nphi, sizeof (double));
      if (jc_bL == NULL && burn >= 1000)
        zrec = (double *) mycalloc ((size_t) ((burn - burn / 2) * d + 1), sizeof (double));
      /* block steps per sweep: about the cost of a single-site sweep */
      nb = ly->ml ? (d + ly->numpop - 1) / ly->numpop : 1;
      if (nb < 1)
        nb = 1;
    }
  boolean skipping = FALSE;
  for (sweep = 0; sweep < burn + sweeps; sweep++)
    {
      if (sweep == burn && beta == 1.0 && sweeps > 0)
        {   /* sampling at beta = 1: skip negligible genealogies (jc_local_skip) */
          jc_eval_skip (E, phi, TRUE);
          skipping = TRUE;
        }
      if (block && jc_bL != NULL && jc_bd == d)
        {
          long b;
          for (b = 0; b < nb; b++)
            {
              boolean ok = TRUE;
              double delta = 0.0;
              for (q = 0; q < d; q++)
                {
                  const long pp = ai[q];
                  z[q] = pp < npos ? log (phi[pp]) : phi[pp];
                  eps[q] = jc_gauss ();
                }
              memcpy (phin, phi, sizeof (double) * (size_t) nphi);
              for (q = 0; q < d && ok; q++)
                {
                  const long pp = ai[q];
                  double s = 0.0;
                  for (k = 0; k <= q; k++)
                    s += jc_bL[q * d + k] * eps[k];
                  zn[q] = z[q] + jc_bscale * s;
                  phin[pp] = pp < npos ? exp (zn[q]) : zn[q];
                  double blo = world->bayes->minparam[pp], bhi = world->bayes->maxparam[pp];
                  if (pp < np)
                    jc_bounds (world, ly, pp, &blo, &bhi);
                  if (!(phin[pp] >= blo && phin[pp] <= bhi) || (pp < npos && !(phin[pp] > 0.0)))
                    ok = FALSE;
                  else
                    delta += jc_target_prior (world, pp, phin[pp]) - jc_target_prior (world, pp, phi[pp])
                      + (pp < npos ? zn[q] - z[q] : 0.0);
                }
              btri++;
              if (ok && !jc_order_ok (ly, phin, -1, 0.0))   /* divergence-order TREE/RANK */
                ok = FALSE;
              double newcvb = cvprior;
              if (ok && cvon)
                {
                  newcvb = jc_cv_logprior (world, ly, phin);
                  if (newcvb > -HUGE_VAL)
                    delta += newcvb - cvprior;
                  else
                    ok = FALSE;
                }
              if (ok)
                {
                  const double dd = jc_eval_block (E, E->accept_prev, phi, phin, nphi);
                  delta += beta * dd;
                  E->accept_prev = log (UNIF_RANDUM ()) < delta;
                  if (E->accept_prev)
                    {
                      bacc++;
                      cvprior = newcvb;
                      memcpy (phi, phin, sizeof (double) * (size_t) nphi);
                      E->logF += dd;
                    }
                }
              if (sweep < burn && btri % 50 == 0)
                {
                  const double rate = (double) bacc / (double) btri;
                  jc_bscale *= rate > 0.234 ? 1.2 : 0.8;
                  bacc = btri = 0;
                }
            }
        }
      else
      for (p = 0; p < nphi; p++)
        {
          if (!active[p])
            continue;
          /* a multiplier: symmetric in log tau, its prior is on log tau */
          const boolean tau = p >= np;
          double lo = exp (-SKYLINE_HISTMAX), hi = exp (SKYLINE_HISTMAX);
          if (!tau)
            jc_bounds (world, ly, p, &lo, &hi);
          const boolean mult = tau || p < npos;
          const double x = mult ? phi[p] * exp (step[p] * (UNIF_RANDUM () - 0.5) * 2.0)
                                : phi[p] + step[p] * (UNIF_RANDUM () - 0.5) * 2.0;
          double newsky = skyprior, newcv = cvprior;
          tri[p]++;
          if (x >= lo && x <= hi && (!mult || x > 0.0) && ly->nfree > 0)
            {
              const double keep = phi[p];
              phi[p] = x;
              newsky = jc_sky_logprior (ly, phi);
              phi[p] = keep;
            }
          if (cvon && x >= lo && x <= hi && (!mult || x > 0.0) && !tau)
            {   /* splitstd-cv: a D or S step changes p(S|D) */
              const double keep = phi[p];
              phi[p] = x;
              newcv = jc_cv_logprior (world, ly, phi);
              phi[p] = keep;
            }
          if (x >= lo && x <= hi && (!mult || x > 0.0) && newsky > -HUGE_VAL && newcv > -HUGE_VAL
              && jc_order_ok (ly, phi, p, x))   /* divergence-order TREE/RANK */
            {
              double delta = (tau ? 0.0 : jc_target_prior (world, p, x) - jc_target_prior (world, p, phi[p])
                              + (mult ? log (x / phi[p]) : 0.0))
                + newsky - skyprior + (newcv - cvprior);
              const double d = jc_eval_step (E, phi, E->accept_prev, p, x);
              delta += beta * d;
              E->accept_prev = log (UNIF_RANDUM ()) < delta;
              if (E->accept_prev)
                {
                  acc[p]++;
                  phi[p] = x;
                  skyprior = newsky;
                  cvprior = newcv;
                  E->logF += d;
                }
            }
          if (sweep < burn && tri[p] % 50 == 0)
            {
              const double rate = (double) acc[p] / (double) tri[p];
              step[p] *= rate > 0.44 ? 1.2 : 0.8;
              acc[p] = tri[p] = 0;
            }
        }
      if (zrec != NULL && sweep >= burn / 2 && sweep < burn)
        {   /* the second half of the burn-in: covariance of the block proposal */
          for (q = 0; q < d; q++)
            zrec[nrec * d + q] = ai[q] < npos ? log (phi[ai[q]]) : phi[ai[q]];
          nrec++;
          if (sweep == burn - 1)
            {
              jc_block_covariance (zrec, nrec, d);
              myfree (zrec);
              zrec = NULL;
            }
        }
      if (sweep >= burn)
        {
          if (trace != NULL)
            for (p = 0; p < nphi; p++)
              trace[p * sweeps + (sweep - burn)] = phi[p];
          if (logftr != NULL)
            logftr[sweep - burn] = E->logF;
          if (essstop > 0 && (sweep - burn + 1) % JC_ESSCHECK == 0 && sweep + 1 < burn + sweeps)
            {   /* ESS stop: every active coordinate */
              const long n = sweep - burn + 1;
              double emin = HUGE_VAL;
              for (p = 0; p < nphi; p++)
                if (active[p])
                  {
                    const double e = jc_ess (trace + p * sweeps, n);
                    if (e < emin)
                      emin = e;
                  }
              if (getenv ("MIGRATE_JC_DEBUG") != NULL)
                fprintf (stderr, "JCESS after %ld sweeps: min ESS %.0f (target %ld)\n", n, emin, essstop);
              if (emin >= (double) essstop)
                {
                  done = n;
                  break;
                }
            }
        }
    }
  if (skipping)
    jc_eval_skip (E, phi, FALSE);
  if (block && getenv ("MIGRATE_JC_DEBUG") != NULL)
    fprintf (stderr, "JCBLOCK chain beta %g burn %ld sweeps %ld: d %ld, %ld block steps per sweep, scale %g\n",
             beta, burn, sweeps, d, nb, jc_bscale);
  if (trace != NULL && done < sweeps)
    for (p = 0; p < nphi; p++)   /* compact to stride done */
      memmove (trace + p * done, trace + p * sweeps, sizeof (double) * (size_t) done);
  myfree (acc);
  myfree (tri);
  myfree (ai);
  myfree (zrec);
  myfree (z);
  myfree (zn);
  myfree (eps);
  myfree (phin);
  return done;
}

/* initial step sizes */
static double *
jc_steps (world_fmt *world, const jc_layout *ly)
{
  const long npos = world->numparamcumvec[SPLITSTDPRIOR];
  double *step = (double *) mycalloc ((size_t) ly->nphi, sizeof (double));
  long p;
  for (p = 0; p < ly->nphi; p++)
    step[p] = (p < npos || p >= ly->np) ? 0.3 : 0.1 * (world->bayes->maxparam[p] - world->bayes->minparam[p]);
  return step;
}

/* bootstrap errors from D (JC_BOOTB x JC_BOOTK): every replicate reweights
   the JC_BOOTK trace points by exp(D_bk); the error of a parameter is the
   sd of its reweighted medians; a parameter is flagged when that error
   exceeds JC_BOOTFLAG posterior sds or the reweighting collapses (median
   ESS below JC_BOOTK / 4); tests (2026-10-03): data sets whose Joint
   estimates were noise had ratios 0.43-1.02 and ESS 42-77, good ones
   0.13-0.25 and ESS 185-196 -- a heuristic, the bootstrap underestimates
   the differences between runs */
#define JC_BOOTFLAG 0.35
static void
jc_boot_errors (jc_store *js, const double *D)
{
  const long np = js->np, K = JC_BOOTK, B = JC_BOOTB;
  double *w = (double *) mycalloc ((size_t) (B * K), sizeof (double));
  double *ess = (double *) mycalloc ((size_t) B, sizeof (double));
  double *med = (double *) mycalloc ((size_t) B, sizeof (double));
  long *ord = (long *) mycalloc ((size_t) K, sizeof (long));
  double *val = (double *) mycalloc ((size_t) K, sizeof (double));
  long b, k, p, i;
  myfree (js->boot_err);
  myfree (js->boot_flag);
  js->boot_err = (double *) mycalloc ((size_t) np, sizeof (double));
  js->boot_flag = (boolean *) mycalloc ((size_t) np, sizeof (boolean));
  for (b = 0; b < B; b++)
    {
      double mx = -HUGE_VAL, sum = 0.0, s2 = 0.0;
      for (k = 0; k < K; k++)
        if (D[b * K + k] > mx)
          mx = D[b * K + k];
      for (k = 0; k < K; k++)
        sum += (w[b * K + k] = exp (D[b * K + k] - mx));
      for (k = 0; k < K; k++)
        {
          w[b * K + k] /= sum;
          s2 += w[b * K + k] * w[b * K + k];
        }
      ess[b] = 1.0 / s2;
    }
  memcpy (med, ess, sizeof (double) * (size_t) B);
  qsort (med, (size_t) B, sizeof (double), jc_cmp_double);
  js->boot_ess = med[B / 2];
  {   /* the scaling factor of the marginal likelihood: C_b / C is the
         posterior mean of exp(D_b), the trace points being posterior draws */
    double m = 0.0, v = 0.0;
    for (b = 0; b < B; b++)
      {
        double mx = -HUGE_VAL, sum = 0.0;
        for (k = 0; k < K; k++)
          if (D[b * K + k] > mx)
            mx = D[b * K + k];
        for (k = 0; k < K; k++)
          sum += exp (D[b * K + k] - mx);
        med[b] = mx + log (sum / K);
        m += med[b] / B;
      }
    for (b = 0; b < B; b++)
      v += (med[b] - m) * (med[b] - m) / (B - 1);
    js->logc_boot = sqrt (v);
  }
  for (p = 0; p < np; p++)
    {
      double m = 0.0, v = 0.0, mean = 0.0, var = 0.0;
      if (!js->active[p])
        continue;
      for (k = 0; k < K; k++)
        val[k] = js->trace[p * js->nsweeps + (long) ((k + 0.5) * js->nsweeps / K)];
      /* points in increasing order of the parameter */
      for (k = 0; k < K; k++)
        ord[k] = k;
      for (k = 1; k < K; k++)
        for (i = k; i > 0 && val[ord[i - 1]] > val[ord[i]]; i--)
          {
            long t = ord[i];
            ord[i] = ord[i - 1];
            ord[i - 1] = t;
          }
      for (b = 0; b < B; b++)
        {
          double c = 0.0;
          for (k = 0; k < K; k++)
            {
              c += w[b * K + ord[k]];
              if (c >= 0.5)
                break;
            }
          med[b] = val[ord[k < K ? k : K - 1]];
          m += med[b] / B;
        }
      for (b = 0; b < B; b++)
        v += (med[b] - m) * (med[b] - m) / (B - 1);
      js->boot_err[p] = sqrt (v);
      for (i = 0; i < js->nsweeps; i++)
        mean += js->trace[p * js->nsweeps + i] / js->nsweeps;
      for (i = 0; i < js->nsweeps; i++)
        var += (js->trace[p * js->nsweeps + i] - mean) * (js->trace[p * js->nsweeps + i] - mean) / js->nsweeps;
      js->boot_flag[p] = (var > 0.0 && js->boot_err[p] > JC_BOOTFLAG * sqrt (var))
        || js->boot_ess < 0.25 * K;
    }
  js->has_boot = TRUE;
  myfree (w);
  myfree (ess);
  myfree (med);
  myfree (ord);
  myfree (val);
}

/* one joint chain from phi0 over the genealogies of part (0: all, k: the
   k-th block of every locus); trace is nphi x sweeps */
static boolean
jc_run (world_fmt *world, const jc_layout *ly, const boolean *active, const double *phi0,
        long part, long burn, long sweeps, double *trace, long *nl, long *tmin, long *tmax, double *rows,
        double *bootD, long *used, double *phiend)
{
  const long nphi = ly->nphi;
  double *phi = (double *) mycalloc ((size_t) nphi, sizeof (double));
  long p, l;
  jc_eval E;
  memcpy (phi, phi0, sizeof (double) * (size_t) nphi);
  jc_cv_repair (world, ly, phi);   /* splitstd-cv: S inside [cmin D, cmax D] */
  const double rt0 = jc_now ();
  jc_eval_start (world, &E, phi, part);
  const double rt1 = jc_now ();
  *nl = E.nl;
  *rows = E.rows;
  *tmin = E.tmin;
  *tmax = E.tmax;
  if (E.nl < 2)
    {
      jc_eval_end (&E, FALSE);
      myfree (phi);
      return FALSE;
    }
  double *step = jc_steps (world, ly);
  /* burn-in, then sampling (the same chain as in one call). A state where a
     locus has no stored genealogy compatible with phi (every row at JC_DOUT
     or minus infinity: its term at the floor) has zero density; the chain
     never accepts a step into one, but it can start in one and, with small
     divergence stds, find no single-coordinate step out of it. Then it
     wandered among impossible states, compared only on the other loci, and
     reported them as a narrow posterior (2026-11-16). Up to two more burn-ins,
     then the combination is given up (All is reported). */
  jc_chain (world, ly, active, &E, phi, step, 1.0, burn, 0, NULL, NULL);
  {
    long extra;
    for (extra = 0; extra < 2 && E.logF < 0.5 * JC_LTERM_FLOOR; extra++)
      jc_chain (world, ly, active, &E, phi, step, 1.0, burn, 0, NULL, NULL);
  }
  if (E.logF < 0.5 * JC_LTERM_FLOOR)
    {
      const long nimp = (long) floor (E.logF / JC_LTERM_FLOOR + 0.5);
      if (part == 0)
        jc_impossible = nimp;
      if (part == 0)
        record_warnings (world, "Joint combination dropped: after %ld burn-in sweeps %ld loci had no stored "
                         "genealogy compatible with the joint parameters (zero density); All is reported. "
                         "Typical cause: a split-time std near zero, where every genealogy fixes the split time.",
                         3 * burn, nimp);
      if (getenv ("MIGRATE_JC_DEBUG") != NULL || getenv ("MIGRATE_JC_DIAG") != NULL)
        fprintf (stderr, "JCIMPOSSIBLE part %ld: %ld loci at zero density after %ld burn-in sweeps\n",
                 part, nimp, 3 * burn);
      jc_eval_end (&E, E.accept_prev);
      myfree (step);
      myfree (phi);
      return FALSE;
    }
  sweeps = jc_chain (world, ly, active, &E, phi, step, 1.0, 0, sweeps, trace, NULL);
  if (getenv ("MIGRATE_JC_DEBUG") != NULL)
    fprintf (stderr, "JCTIME run part %ld: setup %.1f s, chain %.1f s (%ld + %ld sweeps)\n",
             part, rt1 - rt0, jc_now () - rt1, burn, sweeps);
  if (used != NULL)
    *used = sweeps;
  if (phiend != NULL)
    memcpy (phiend, phi, sizeof (double) * (size_t) nphi);
  if (part == 0 && getenv ("MIGRATE_JC_TRACE") != NULL)
    {   /* experiments: the main trace, one sweep per line */
      FILE *tf = fopen (getenv ("MIGRATE_JC_TRACE"), "w");
      long sw;
      for (sw = 0; tf != NULL && sw < sweeps; sw++)
        {
          for (p = 0; p < nphi; p++)
            fprintf (tf, p ? " %.10g" : "%.10g", trace[p * sweeps + sw]);
          fprintf (tf, "\n");
        }
      if (tf != NULL)
        fclose (tf);
    }
  if (bootD != NULL)
    {   /* bootstrap guard at JC_BOOTK evenly spaced points of the trace */
      double *pts = (double *) mycalloc ((size_t) (JC_BOOTK * nphi), sizeof (double));
      long k;
      for (k = 0; k < JC_BOOTK; k++)
        for (p = 0; p < nphi; p++)
          pts[k * nphi + p] = trace[p * sweeps + (long) ((k + 0.5) * sweeps / JC_BOOTK)];
      memset (bootD, 0, sizeof (double) * (size_t) (JC_BOOTB * JC_BOOTK));
      jc_eval_boot (&E, nphi, pts, JC_BOOTK, JC_BOOTB, bootD);
      myfree (pts);
    }
  if (part == 0 && getenv ("MIGRATE_JC_DIAG") != NULL && !E.mpi)
    {   /* diagnostic: importance-sampling ESS of each locus' mixture over
           its genealogies, at 20 points of the joint trace */
      jc_local *J = &E.local;
      double *x = (double *) mycalloc ((size_t) nphi, sizeof (double));
      double *w = (double *) mycalloc ((size_t) (J->maxn + 1), sizeof (double));
      long t, sidx;
      for (l = 0; l < J->nl; l++)
        {
          double esum = 0.0, emin = HUGE_VAL;
          jc_hcur = J->L[l].h;
          for (t = 0; t < 20; t++)
            {
              const long it = (sweeps / 20) * t + sweeps / 40;
              double mx = -HUGE_VAL, s1 = 0.0, s2 = 0.0;
              for (p = 0; p < nphi; p++)
                x[p] = trace[p * sweeps + it];
              for (sidx = 0; sidx < J->L[l].n; sidx++)
                {
                  const double *st = J->L[l].st + sidx * ly->nrow;
                  w[sidx] = jc_full_term (&J->ly, J->active, st, JC_SP (&J->L[l], sidx), x, st[J->ly.off_mu]) - J->L[l].logz[sidx];
                  if (w[sidx] > mx)
                    mx = w[sidx];
                }
              for (sidx = 0; sidx < J->L[l].n; sidx++)
                {
                  const double e = exp (w[sidx] - mx);
                  s1 += e;
                  s2 += e * e;
                }
              esum += s1 * s1 / s2;
              if (s1 * s1 / s2 < emin)
                emin = s1 * s1 / s2;
            }
          fprintf (stderr, "JCDIAG locus %ld n %ld IS-ESS mean %.1f min %.1f\n", l, J->L[l].n, esum / 20.0, emin);
          {   /* the dominant genealogy at the last trace point: its share,
                 split-event centre/count, log Z against the locus median */
            long top = 0, k2;
            double mx = -HUGE_VAL, s1 = 0.0;
            for (p = 0; p < nphi; p++)
              x[p] = trace[p * sweeps + sweeps - 1];
            for (sidx = 0; sidx < J->L[l].n; sidx++)
              {
                const double *st = J->L[l].st + sidx * ly->nrow;
                w[sidx] = jc_full_term (&J->ly, J->active, st, JC_SP (&J->L[l], sidx), x, st[J->ly.off_mu]) - J->L[l].logz[sidx];
                if (w[sidx] > mx)
                  {
                    mx = w[sidx];
                    top = sidx;
                  }
              }
            for (sidx = 0; sidx < J->L[l].n; sidx++)
              s1 += exp (w[sidx] - mx);
            double *lz = (double *) mycalloc ((size_t) (J->L[l].n + 1), sizeof (double));
            memcpy (lz, J->L[l].logz, sizeof (double) * (size_t) J->L[l].n);
            qsort (lz, (size_t) J->L[l].n, sizeof (double), jc_cmp_double);
            const double *st = J->L[l].st + top * ly->nrow;
            fprintf (stderr, "JCTOP locus %ld share %.3f logz %.2f (median %.2f, min %.2f) term %.2f",
                     l, 1.0 / s1, J->L[l].logz[top], lz[J->L[l].n / 2], lz[0],
                     jc_full_term (&J->ly, J->active, st, JC_SP (&J->L[l], top), x, st[J->ly.off_mu]));
            for (k2 = 0; k2 < ly->nsplit; k2++)
              fprintf (stderr, " split%ld c %.6f n %.0f", k2, st[ly->off_D[k2]], st[ly->off_D[k2] + 1]);
            fprintf (stderr, "\n");
            myfree (lz);
          }
        }
      myfree (x);
      myfree (w);
    }
  jc_eval_end (&E, E.accept_prev);
  myfree (phi);
  myfree (step);
  return TRUE;
}

/* the joint scaling factor of the multi-locus marginal likelihood,
   log int prior(phi) prod_l f_l(phi) dphi (f_l = p_l(phi|D_l)/prior(phi)),
   by stepping-stone sampling over a power beta on prod_l f_l: rungs
   beta_k = (k/K)^4, log C = sum_k log E_{beta_k}[F^(beta_{k+1} - beta_k)],
   each expectation from a chain at beta_k (from the posterior down to the
   prior, warm-started); the error from 5 batches per rung */
static boolean
jc_ss (world_fmt *world, const jc_layout *ly, const boolean *active, const double *phi0,
       double *logc, double *err)
{
  const long nphi = ly->nphi;
  /* experiments: MIGRATE_JC_SSRUNGS, MIGRATE_JC_SSN override the ladder */
  const long nr = getenv ("MIGRATE_JC_SSRUNGS") ? atol (getenv ("MIGRATE_JC_SSRUNGS")) : JC_SSRUNGS;
  const long ns = getenv ("MIGRATE_JC_SSN") ? atol (getenv ("MIGRATE_JC_SSN")) : JC_SSN;
  double *phi = (double *) mycalloc ((size_t) nphi, sizeof (double));
  double *lf = (double *) mycalloc ((size_t) ns, sizeof (double));
  double var = 0.0, sum = 0.0;
  long k, i, b;
  jc_eval E;
  memcpy (phi, phi0, sizeof (double) * (size_t) nphi);
  jc_eval_start (world, &E, phi, 0);
  if (E.nl < 2)
    {
      jc_eval_end (&E, FALSE);
      myfree (phi);
      myfree (lf);
      return FALSE;
    }
  double *step = jc_steps (world, ly);
  jc_chain (world, ly, active, &E, phi, step, 1.0, JC_BURNWARM, 0, NULL, NULL);   /* phi0: the main run's last state */
  for (k = nr - 1; k >= 0; k--)
    {
      const double b0 = pow ((double) k / nr, 4.0);
      const double db = pow ((double) (k + 1) / nr, 4.0) - b0;
      double rb[5], mx = -HUGE_VAL, acc = 0.0, m = 0.0, v = 0.0;
      jc_chain (world, ly, active, &E, phi, step, b0, JC_SSBURN, ns, NULL, lf);
      for (i = 0; i < ns; i++)
        if (db * lf[i] > mx)
          mx = db * lf[i];
      for (i = 0; i < ns; i++)
        acc += exp (db * lf[i] - mx);
      sum += mx + log (acc / ns);
      for (b = 0; b < 5; b++)
        {
          double a2 = 0.0;
          for (i = b * (ns / 5); i < (b + 1) * (ns / 5); i++)
            a2 += exp (db * lf[i] - mx);
          rb[b] = mx + log (a2 / (ns / 5));
          m += rb[b] / 5.0;
        }
      for (b = 0; b < 5; b++)
        v += (rb[b] - m) * (rb[b] - m) / 4.0;
      var += v / 5.0;
    }
  jc_eval_end (&E, E.accept_prev);
  /* f_l is relative to the prior mass N_Z the normalizers integrate (with
     skyline=PARAM the range rule truncates it), the per-locus marginal
     likelihoods to the normalised prior: + L log N_Z */
  *logc = sum + E.lognz;
  *err = sqrt (var);
  myfree (phi);
  myfree (lf);
  myfree (step);
  return isfinite (sum);
}

/// Runs the joint combination once (master); TRUE when results are available.
/// The genealogy statistics never leave the process that holds them.
boolean
jc_combine (world_fmt *world)
{
  jc_store *js = (jc_store *) world->jointstats;
  jc_layout ly;
  long l, p, s;
  if (world->loci < 2 || !jc_supported (world))
    return FALSE;
  if (js != NULL && js->done)
    return TRUE;
  myfree (jc_bL);   /* block proposals: a new covariance for this run */
  jc_bL = NULL;
  jc_bd = 0;
  jc_layout_make (world, &ly);
  const long np = ly.np, nphi = ly.nphi;
  /* positive parameters (Theta, M, divergence) get multiplicative moves,
     growth rates additive ones */
  const long npos = world->numparamcumvec[SPLITSTDPRIOR];
  for (s = 0; s < ly.nslot; s++)
    if (ly.slot_param[s] < 0)
      {
        jc_layout_free (&ly);
        return FALSE;
      }
  boolean *active = jc_active (world, &ly);
  /* start: mean over loci of the per-locus posterior means */
  double *phi = (double *) mycalloc ((size_t) nphi, sizeof (double));
  for (p = 0; p < np; p++)
    {
      phi[p] = world->param0[p];
      if (!active[p])
        continue;
      double sum = 0.0;
      long cnt = 0;
      for (l = 0; l < world->loci; l++)
        if (!world->data->skiploci[l])
          {
            sum += world->bayes->histogram[l].means[p];
            cnt++;
          }
      const double lo = world->bayes->minparam[p], hi = world->bayes->maxparam[p];
      phi[p] = cnt ? sum / (double) cnt : 0.5 * (lo + hi);
      if (phi[p] <= lo || phi[p] >= hi)
        phi[p] = 0.5 * (lo + hi);
      if (p < npos && phi[p] <= 0.0)
        phi[p] = 0.01 * hi;
    }
  /* the shared rate: the union over loci of their 99% ranges, from each
     locus' 95% interval widened on the log scale by 2.576 / 1.96, inside the
     rate grid */
  jc_rate_lo = jc_rate_hi = 0.0;
  if (ly.rate)
    {
      const long p = ly.rate_p;
      const double glo = exp (ly.rlo), ghi = exp (ly.rhi);
      double ulo = HUGE_VAL, uhi = 0.0;
      for (l = 0; l < world->loci; l++)
        if (!world->data->skiploci[l])
          {
            const bayeshistogram_fmt *h = &world->bayes->histogram[l];
            double a = h->cred95l[p], b = h->cred95u[p];
            if (!(a > 0.0) || a < glo)
              a = glo;
            if (!(b > a))
              b = ghi;
            const double c = 0.5 * (log (a) + log (b)), w = 0.5 * (log (b) - log (a)) * 2.576 / 1.96;
            if (exp (c - w) < ulo)
              ulo = exp (c - w);
            if (exp (c + w) > uhi)
              uhi = exp (c + w);
          }
      jc_rate_lo = ulo > glo ? ulo : glo;
      jc_rate_hi = uhi < ghi ? uhi : ghi;
      if (!(jc_rate_hi > jc_rate_lo))
        {
          jc_rate_lo = glo;
          jc_rate_hi = ghi;
        }
      if (phi[p] <= jc_rate_lo || phi[p] >= jc_rate_hi)
        phi[p] = sqrt (jc_rate_lo * jc_rate_hi);
      if (getenv ("MIGRATE_JC_DEBUG") != NULL)
        fprintf (stderr, "JCRATE grid [%g, %g] (%ld points), joint range [%g, %g], start %g\n",
                 glo, ghi, ly.nrk, jc_rate_lo, jc_rate_hi, phi[p]);
    }
  /* divergence-order TREE/RANK: a start inside the order (the per-locus means
     already are when the per-locus chains were ordered) */
  for (long pass = 0; pass <= ly.nsplit && !jc_order_ok (&ly, phi, -1, 0.0); pass++)
    for (long k = 0; k < ly.nsplit; k++)
      if (ly.split_par[k] >= 0 && phi[ly.split_pmu[k]] > phi[ly.split_pmu[ly.split_par[k]]])
        {   /* a fixed D stays: then the parent moves up to it */
          if (active[ly.split_pmu[k]])
            phi[ly.split_pmu[k]] = phi[ly.split_pmu[ly.split_par[k]]];
          else if (active[ly.split_pmu[ly.split_par[k]]])
            phi[ly.split_pmu[ly.split_par[k]]] = phi[ly.split_pmu[k]];
        }
  /* skyline multipliers: 1 (or the geometric mean of log-uniform bounds
     that exclude 1), moved inside the range rule */
  for (long j = 0; j < ly.nfree; j++)
    for (s = 1; s < ly.nseg; s++)
      {
        double tau = 1.0;
        const double base = phi[ly.rep[ly.sky_k[j]]];
        if (ly.sky_type[j] == SKYPRIOR_LOGUNIFORM && (ly.sky_a[j] > 1.0 || ly.sky_b[j] < 1.0))
          tau = sqrt (ly.sky_a[j] * ly.sky_b[j]);
        if (base * tau > ly.sky_hi[j])
          tau = 0.99 * ly.sky_hi[j] / base;
        if (base * tau < ly.sky_lo[j])
          tau = 1.01 * ly.sky_lo[j] / base;
        phi[JC_TAU (&ly, j, s)] = tau;
      }
  js = jc_get_store (world, ly.nrow);
  myfree (js->trace);
  myfree (js->active);
  myfree (js->sky_k);
  myfree (js->mcerr);
  myfree (js->mchalf);
  {   /* main run: a budget of proposals, so that a large model does not
         take sweeps x parameters x genealogies x loci evaluations */
    long nact = 0;
    for (p = 0; p < nphi; p++)
      if (active[p])
        nact++;
    js->nsweeps = JC_PROPBUDGET / (nact > 0 ? nact : 1);
    if (js->nsweeps > JC_SWEEPS)
      js->nsweeps = JC_SWEEPS;
    if (js->nsweeps < JC_MINSWEEPS)
      js->nsweeps = JC_MINSWEEPS;
    if (getenv ("MIGRATE_JC_SWEEPS") != NULL)   /* experiments */
      js->nsweeps = atol (getenv ("MIGRATE_JC_SWEEPS"));
  }
  js->trace = (double *) mycalloc ((size_t) (nphi * js->nsweeps), sizeof (double));
  js->active = active;
  js->np = nphi;
  js->numparam = np;
  js->nseg = ly.nseg;
  js->nfree = ly.nfree;
  js->sky_k = (long *) mycalloc ((size_t) (ly.nfree + 1), sizeof (long));
  memcpy (js->sky_k, ly.sky_k, sizeof (long) * (size_t) ly.nfree);
  long nl, tmin, tmax;
  double rows_all, rows_blk = 0.0, rb;
  double *bootD = (double *) mycalloc ((size_t) (JC_BOOTB * JC_BOOTK), sizeof (double));
  long used = js->nsweeps;
  const double jt0 = jc_now ();
  double *phiw = (double *) mycalloc ((size_t) nphi, sizeof (double));   /* the main run's last state */
  if (!jc_run (world, &ly, active, phi, 0, JC_BURNIN, js->nsweeps, js->trace, &nl, &tmin, &tmax, &rows_all, bootD, &used, phiw))
    {
      myfree (phi);
      myfree (phiw);
      myfree (bootD);
      jc_layout_free (&ly);
      return FALSE;
    }
  js->nsweeps = used;   /* the ESS stop may end the main run early (trace compacted) */
  const double jt1 = jc_now ();
  jc_boot_errors (js, bootD);
  myfree (bootD);
  if (js->has_boot && js->boot_ess < 0.25 * JC_BOOTK)
    record_warnings (world, "Joint combination collapsed: the median reweighting ESS is %.0f (below %d), "
                     "a few stored genealogies carry almost all the weight. Every Joint value is flagged "
                     "(Joint*, the tables and plots use All) and the joint marginal likelihood is not reliable: "
                     "compare models with the row for independent parameters",
                     js->boot_ess, JC_BOOTK / 4);
  else if (js->has_boot)
    {   /* split models: the split-time terms of a stored genealogy are narrow
           spikes at its split events; the product over loci can then lock on
           a few of them, and more genealogies do not help (2026-11-18: three
           seeds of a 20-locus run gave non-overlapping joint intervals of D,
           at 2000 and at 20000 genealogies per locus) */
      long k, nf = 0;
      for (k = 0; k < ly.nsplit; k++)
        if (ly.split_pmu[k] >= 0 && ly.split_pmu[k] < js->np && js->boot_flag[ly.split_pmu[k]])
          nf++;
      if (nf > 0)
        record_warnings (world, "Joint split time flagged (Joint*) for %ld split%s: in split models the joint "
                         "estimates of the split time D and its spread S can be far off even with a reweighting "
                         "ESS above %d (independent runs disagreed beyond their intervals); use the All values "
                         "for D and S", nf, nf > 1 ? "s" : "", JC_BOOTK / 4);
    }
  /* Monte Carlo error: the same combination on each block of every locus'
     genealogies (its replicates, or consecutive stretches of its chain);
     sd(block medians) scaled by the share of genealogies a block holds is
     about the standard error of the median from all */
  {
    const long nb = jc_nblocks (world);
    /* experiments: MIGRATE_JC_MCSWEEPS sets the sweeps per block */
    const long hs = getenv ("MIGRATE_JC_MCSWEEPS") ? atol (getenv ("MIGRATE_JC_MCSWEEPS")) : JC_MCSWEEPS;
    double *th = (double *) mycalloc ((size_t) (nphi * hs), sizeof (double));
    double *srt = (double *) mycalloc ((size_t) (hs > js->nsweeps ? hs : js->nsweeps), sizeof (double));
    long part, n1, t1, t2;
    js->mcerr = (double *) mycalloc ((size_t) nphi, sizeof (double));
    js->mchalf = (double *) mycalloc ((size_t) (nb * nphi), sizeof (double));
    js->has_mcerr = TRUE;
    js->nblock = nb;
    js->byrep = (world->options->replicate ? world->options->replicatenum : 1) >= 2;
    for (part = 1; part <= nb && js->has_mcerr; part++)
      {
        if (!jc_run (world, &ly, active, phiw, part, JC_BURNWARM, hs, th, &n1, &t1, &t2, &rb, NULL, NULL, NULL) || n1 < 2)
          {
            js->has_mcerr = FALSE;
            break;
          }
        rows_blk += rb / nb;
        if (getenv ("MIGRATE_JC_DIAG") != NULL)
          fprintf (stderr, "JCDIAG block %ld rows %.0f (all %.0f)\n", part, rb, rows_all);
        for (p = 0; p < nphi; p++)
          if (active[p])
            {
              memcpy (srt, th + p * hs, sizeof (double) * (size_t) hs);
              qsort (srt, (size_t) hs, sizeof (double), jc_cmp_double);
              js->mchalf[(part - 1) * nphi + p] = srt[hs / 2];
            }
      }
    for (p = 0; p < nphi && js->has_mcerr; p++)
      {
        double m = 0.0, v = 0.0;
        long k;
        for (k = 0; k < nb; k++)
          m += js->mchalf[k * nphi + p] / nb;
        for (k = 0; k < nb; k++)
          v += (js->mchalf[k * nphi + p] - m) * (js->mchalf[k * nphi + p] - m) / (nb - 1);
        /* the blocks partition the genealogies of the full combination: a
           block carries rows_blk / rows_all (about 1/nb) of its information */
        js->mcerr[p] = sqrt (v * rows_blk / rows_all);
      }
    myfree (th);
    myfree (srt);
  }
  /* the joint marginal likelihood (its scaling factor) */
  const double jt2 = jc_now ();
  js->has_logc = jc_ss (world, &ly, active, phiw, &js->logc, &js->logc_err);
  if (getenv ("MIGRATE_JC_DEBUG") != NULL)
    fprintf (stderr, "JCTIME main+boot %.1f s, error blocks %.1f s, stepping stones %.1f s\n",
             jt1 - jt0, jt2 - jt1, jc_now () - jt2);
  if (getenv ("MIGRATE_JC_DIAG") != NULL && js->has_logc)
    fprintf (stderr, "JCDIAG logC %.4f err %.4f\n", js->logc, js->logc_err);
  if (js->has_logc && js->has_boot)   /* plus the genealogy-sampling error */
    js->logc_err = sqrt (js->logc_err * js->logc_err + js->logc_boot * js->logc_boot);
  /* tied splits: their own mean/std indices report the shared values */
  for (long k = 0; k < ly.nsplit; k++)
    if (ly.split_rep[k] != k)
      {
        const species_fmt *sp = &world->species_model[ly.split_model[k]];
        const long pairs[2][2] = {{sp->paramindex_mu, ly.split_pmu[k]},
                                  {sp->paramindex_sigma, ly.split_psig[k]}};
        long t;
        for (t = 0; t < 2; t++)
          if (pairs[t][0] >= 0 && pairs[t][0] < np && pairs[t][1] >= 0 && pairs[t][0] != pairs[t][1])
            {
              memcpy (js->trace + pairs[t][0] * js->nsweeps, js->trace + pairs[t][1] * js->nsweeps,
                      sizeof (double) * (size_t) js->nsweeps);
              if (js->has_boot)
                {
                  js->boot_err[pairs[t][0]] = js->boot_err[pairs[t][1]];
                  js->boot_flag[pairs[t][0]] = js->boot_flag[pairs[t][1]];
                }
              if (js->has_mcerr)
                {
                  long kb;
                  js->mcerr[pairs[t][0]] = js->mcerr[pairs[t][1]];
                  for (kb = 0; kb < js->nblock; kb++)
                    js->mchalf[kb * nphi + pairs[t][0]] = js->mchalf[kb * nphi + pairs[t][1]];
                }
              active[pairs[t][0]] = TRUE;
            }
      }
  js->nloci_used = nl;
  js->tmin = tmin;
  js->tmax = tmax;
  js->done = TRUE;
  myfree (phi);
  myfree (phiw);
  jc_layout_free (&ly);
  return TRUE;
}

/* The joint posterior conditions on the stored genealogies; their sampling
   error is not part of it. With e the larger of the bootstrap error of the
   median (jc_boot_errors()) and the Monte Carlo error from the blocks of
   genealogies (mcerr; the bootstrap blocks of 10 genealogies miss longer
   autocorrelation and gave 2-3 times smaller errors), the trace of
   parameter p is widened around its median by c = sqrt(1 + (e / sd)^2), on
   the log scale (the values are positive) and kept inside the prior range,
   so that its spread is about sqrt(sd^2 + e^2). In a calibration with truths
   drawn from the prior (50 loci, 40 replicates) this raised the 95%
   coverage of M from 0.90 / 0.80 to about 0.93 / 0.85: the measured
   genealogy error explains only part of the shortfall. */
static void
jc_widen (world_fmt *world, const jc_store *js, long p, double *x)
{
  double mean = 0.0, var = 0.0, med, c;
  long i;
  double e = 0.0;
  if (js->has_boot && js->boot_err != NULL)
    e = js->boot_err[p];
  if (js->has_mcerr && js->mcerr != NULL && js->mcerr[p] > e)
    e = js->mcerr[p];
  if (e <= 0.0)
    return;
  for (i = 0; i < js->nsweeps; i++)
    {
      if (x[i] <= 0.0)
        return;
      mean += x[i] / js->nsweeps;
    }
  for (i = 0; i < js->nsweeps; i++)
    var += (x[i] - mean) * (x[i] - mean) / js->nsweeps;
  if (var <= 0.0)
    return;
  c = sqrt (1.0 + e * e / var);
  {
    double *srt = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
    memcpy (srt, x, sizeof (double) * (size_t) js->nsweeps);
    qsort (srt, (size_t) js->nsweeps, sizeof (double), jc_cmp_double);
    med = srt[js->nsweeps / 2];
    myfree (srt);
  }
  const double lmed = log (med);
  const double plo = world->bayes->minparam[p], phi = world->bayes->maxparam[p];
  for (i = 0; i < js->nsweeps; i++)
    {
      double y = exp (lmed + c * (log (x[i]) - lmed));
      if (y < plo)
        y = plo;
      if (y > phi)
        y = phi;
      x[i] = y;
    }
}

/// Replaces the combined ("All") histogram of every handled parameter with
/// the joint samples, in calc_hpd_credibility()'s bin layout; the caller
/// then computes modes, medians and HPD intervals as usual.
boolean
jc_fill_histogram (world_fmt *world, bayeshistogram_fmt *hist)
{
  jc_store *js = (jc_store *) world->jointstats;
  long pa0, pa, off = 0, i;
  if (js == NULL || !js->done)
    return FALSE;
  for (pa0 = 0; pa0 < world->numparam; pa0++)
    {
      if (shortcut (pa0, world, &pa) || pa < pa0)
        continue;
      const long nb = hist->bins[pa];
      /* a flagged parameter (bootstrap guard) keeps its joint estimate: the
         product of the marginals is worse, and more so with more loci (in a
         coverage study at 50 and 200 loci its 95% intervals missed the truth
         for most parameters, while Joint covered 0.8-1.0) */
      if (pa < js->np && js->active[pa] && nb > 0)
        {
          const double lo = hist->minima[pa];
          const double *x0 = js->trace + pa * js->nsweeps;
          double *x = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
          double mean = 0.0;
          memcpy (x, x0, sizeof (double) * (size_t) js->nsweeps);
          jc_widen (world, js, pa, x);
          memset (hist->results + off, 0, sizeof (double) * (size_t) nb);
          for (i = 0; i < js->nsweeps; i++)
            {
              long b = bayes_hist_bin (world->bayes, pa, lo, x[i]);
              if (b < 0)
                b = 0;
              if (b >= nb)
                b = nb - 1;
              hist->results[off + b] += 1.0 / (double) js->nsweeps;
              mean += x[i];
            }
          hist->means[pa] = mean / (double) js->nsweeps;
          myfree (x);
        }
      off += nb;
    }
  return TRUE;
}

/* name of joint coordinate p: a parameter, or a skyline multiplier */
static void
jc_param_name (world_fmt *world, const jc_store *js, long p, char *name)
{
  const long npx = world->numparamcumvec[SPLITSTDPRIOR];
  const long npg = world->numparamcumvec[GROWTHPRIOR];
  char buf[LINESIZE];
  long i, n = 0;
  if (p >= js->numparam)
    {
      const long ns = js->nseg - 1, j = (p - js->numparam) / ns, sg = (p - js->numparam) % ns + 1;
      jc_param_name (world, js, js->sky_k[j], buf);
      snprintf (name, LINESIZE, "tau(%s) seg %ld", buf, sg + 1);
      return;
    }
  if (p >= npx && p < npg)
    snprintf (buf, LINESIZE, "Growth_%ld", world->options->growpops[p - npx]);
  else if (p < npx)
    set_paramstr (buf, p, world);
  else
    snprintf (buf, LINESIZE, "Param_%ld", p + 1);
  for (i = 0; buf[i] != '\0'; i++)   /* without padding */
    if (buf[i] != ' ')
      name[n++] = buf[i];
  name[n] = '\0';
}

/// rows of the Monte Carlo error table of the joint estimates (one per
/// sampled joint coordinate); returns their number, 0 without the table
long
jc_mcerr_rows (world_fmt *world, jc_mcerr_row **rows)
{
  jc_store *js = (jc_store *) world->jointstats;
  double *srt;
  long p, i, n = 0;
  *rows = NULL;
  if (js == NULL || !js->done || !js->has_mcerr)
    return 0;
  srt = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
  *rows = (jc_mcerr_row *) mycalloc ((size_t) (js->np + 1), sizeof (jc_mcerr_row));
  for (p = 0; p < js->np; p++)
    {
      long j, kb;
      jc_mcerr_row *r;
      if (!js->active[p])
        continue;
      if (p < js->numparam && (shortcut (p, world, &j) || j != p))
        continue;
      r = &(*rows)[n++];
      const double *tr = js->trace + p * js->nsweeps;
      double mean = 0.0, var = 0.0;
      for (i = 0; i < js->nsweeps; i++)
        mean += tr[i] / js->nsweeps;
      for (i = 0; i < js->nsweeps; i++)
        var += (tr[i] - mean) * (tr[i] - mean) / js->nsweeps;
      memcpy (srt, tr, sizeof (double) * (size_t) js->nsweeps);
      qsort (srt, (size_t) js->nsweeps, sizeof (double), jc_cmp_double);
      r->median = srt[js->nsweeps / 2];
      r->sd = sqrt (var);
      r->mcerr = js->mcerr[p];
      r->ratio = r->sd > 0.0 ? r->mcerr / r->sd : 0.0;
      r->booterr = js->has_boot ? js->boot_err[p] : -1.0;
      r->flagged = js->has_boot && js->boot_flag[p];
      r->blo = HUGE_VAL;
      r->bhi = -HUGE_VAL;
      for (kb = 0; kb < js->nblock; kb++)
        {
          if (js->mchalf[kb * js->np + p] < r->blo)
            r->blo = js->mchalf[kb * js->np + p];
          if (js->mchalf[kb * js->np + p] > r->bhi)
            r->bhi = js->mchalf[kb * js->np + p];
        }
      jc_param_name (world, js, p, r->name);
    }
  myfree (srt);
  return n;
}

/// how the Monte Carlo errors were obtained: blocks, whether they are
/// replicates, the bootstrap's median reweighting ESS (-1 without it)
void
jc_mcerr_info (world_fmt *world, long *nblock, boolean *byrep, double *bootess)
{
  jc_store *js = (jc_store *) world->jointstats;
  *nblock = js ? js->nblock : 0;
  *byrep = js ? js->byrep : FALSE;
  *bootess = (js && js->has_boot) ? js->boot_ess : -1.0;
}

/* the Monte Carlo error table of the joint estimates */
static void
jc_print_mcerr (world_fmt *world, const jc_store *js, FILE *out)
{
  jc_mcerr_row *rows;
  long n = jc_mcerr_rows (world, &rows), i, flagged = 0, nf = 0;
  FPRINTF (out, "\nMonte Carlo error of the joint (All) estimates: the genealogies used by the\n"
                "combination are split into %ld %s, and the\n"
                "combination is repeated on each block. MC error = sd(block medians) x\n"
                "sqrt(genealogies per block / genealogies used). The blocks of one run are\n"
                "correlated, so this is a lower bound: in tests, independent runs differed\n"
                "%s.\n",
           js->nblock, js->byrep ? "groups of replicates" : "consecutive stretches of each chain",
           js->byrep ? "1.3-1.4 times more" : "1.7-1.8 times more (1.3-1.4 with replicate=YES:n)");
  FPRINTF (out, "A ratio MC error / posterior sd above 0.25 (*) asks for longer runs or more\n"
                "replicates.\n\n");
  if (js->has_boot)
    FPRINTF (out, "Boot err: the error of the median from the genealogy sampling alone (block\n"
                  "bootstrap over each locus' genealogies, %d replicates reweighting %d trace points;\n"
                  "median reweighting ESS %.0f). Above %.2f posterior sd, or with a reweighting ESS\n"
                  "below %d, the Joint estimate has a large genealogy-sampling error (*): it is\n"
                  "shown in the Joint* row but should not be trusted, and the plots and the other\n"
                  "tables use the All values for it; more genealogies per locus (longer chains or\n"
                  "replicates) are needed to confirm it.\n\n",
             JC_BOOTB, JC_BOOTK, js->boot_ess, JC_BOOTFLAG, JC_BOOTK / 4);
  FPRINTF (out, "Parameter                    Median     MC error   Blocks: lowest  highest  Post. sd   Ratio   Boot err\n");
  FPRINTF (out, "------------------------------------------------------------------------------------------------------\n");
  for (i = 0; i < n; i++)
    {
      const jc_mcerr_row *r = &rows[i];
      if (r->ratio > 0.25)
        flagged++;
      if (r->flagged)
        nf++;
      FPRINTF (out, "%-26.26s %10.5g %10.5g %10.5g %10.5g %10.5g %7.3f%s %10.5g%s\n", r->name, r->median,
               r->mcerr, r->blo, r->bhi, r->sd, r->ratio, r->ratio > 0.25 ? "*" : " ",
               r->booterr > 0.0 ? r->booterr : 0.0, r->flagged ? "  * noisy" : "");
    }
  if (flagged)
    FPRINTF (out, "(*) %ld estimate%s with a large Monte Carlo error\n", flagged, flagged > 1 ? "s" : "");
  if (nf)
    FPRINTF (out, "(*) %ld Joint estimate%s with a large genealogy-sampling error (Joint* rows;\n"
                  "    the plots and the other tables use All for %s)\n",
             nf, nf > 1 ? "s" : "", nf > 1 ? "them" : "it");
  myfree (rows);
}

/* splitstd-cv: posterior of c = S/D per split from the joint samples (the
   per-locus and All rows hold S and D only as marginals, so c needs the
   joint trace) */
static void
jc_print_splitcv (world_fmt *world, FILE *out)
{
  jc_store *js = (jc_store *) world->jointstats;
  const double pq[5] = { 0.025, 0.25, 0.5, 0.75, 0.975 };
  long k, i, a, nflag = 0;
  if (!world->options->splitstd_cv || !world->has_speciation || world->species_model_size < 1)
    return;
  if (js == NULL || !js->done || js->trace == NULL || js->nsweeps < 2)
    {
      FPRINTF (out, "\nSplit-time spread relative to the split time, c = S/D: needs the joint\n"
                    "multi-locus combination, which did not run (the S and D rows are marginals)\n");
      return;
    }
  FPRINTF (out, "\nSplit-time spread relative to the split time, c = S/D (splitstd-cv; joint\n"
                "samples; S = std. dev. of the lineage split times around their mean D,\n"
                "backward in time)\n");
  FPRINTF (out, "Split        2.5%%       25%%     median      75%%     97.5%%       mean\n");
  FPRINTF (out, "----------------------------------------------------------------------\n");
  double *c = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
  for (k = 0; k < world->species_model_size; k++)
    {
      const species_fmt *s = &world->species_model[k];
      const long pm = s->paramindex_mu, ps = s->paramindex_sigma;
      double mean = 0.0, q[5];
      char name[32];
      boolean flag;
      if (pm < 0 || ps < 0 || pm >= js->np || ps >= js->np)
        continue;
      for (i = 0; i < js->nsweeps; i++)
        {
          const double D = js->trace[pm * js->nsweeps + i];
          c[i] = D > 0.0 ? js->trace[ps * js->nsweeps + i] / D : 0.0;
          mean += c[i] / (double) js->nsweeps;
        }
      qsort (c, (size_t) js->nsweeps, sizeof (double), jc_cmp_double);
      for (a = 0; a < 5; a++)
        q[a] = c[(long) (pq[a] * (js->nsweeps - 1) + 0.5)];
      flag = jc_param_flagged (world, pm) || jc_param_flagged (world, ps);
      if (flag)
        nflag++;
      snprintf (name, sizeof (name), "c_%li->%li%s", s->from + 1, s->to + 1, flag ? "*" : "");
      FPRINTF (out, "%-10.10s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", name,
               q[0], q[1], q[2], q[3], q[4], mean);
    }
  myfree (c);
  if (nflag)
    FPRINTF (out, "(*) S or D is a Joint* estimate (large genealogy-sampling error)\n");
}

/// One line under the Bayesian estimates table saying how "All" was combined
void
jc_print_note (world_fmt *world, FILE *out)
{
  jc_store *js = (jc_store *) world->jointstats;
  if (out == NULL || world->loci < 2)
    return;
  if (js != NULL && js->done)
    FPRINTF (out, "(All)   the product of the per-locus marginal posteriors (each parameter on its own)\n"
                  "(Joint) the joint multi-locus combination over %li loci (%li-%li stored\n"
                  "        genealogies per locus, up to %li used per locus); the plots and the\n"
                  "        other tables use it, except for Joint* parameters (they use All)\n",
             js->nloci_used, js->tmin, js->tmax, jc_maxsamples (world, js->nrow));
  if (js != NULL && js->done)
    {
      /* the split-time spread S: a stored genealogy's split term is a product
         of spikes of width S, one per crossing lineage, so few genealogies fit
         a small S and the product over loci drifts to larger S (2026-11-25:
         NODATA, 20 loci, Joint S 2-4x the prior). All is no alternative: one
         locus barely informs S, and the product of near-flat marginals was
         far off (02-divergence two.txt: truth 0.001, Joint 0.0037, All 0.11) */
      long k;
      boolean hass = FALSE;
      for (k = 0; k < world->species_model_size && world->has_speciation && !world->options->splitstd_zero; k++)
        if (world->species_model[k].paramindex_sigma >= 0
            && world->bayes->maxparam[world->species_model[k].paramindex_sigma]
               > world->bayes->minparam[world->species_model[k].paramindex_sigma])
          hass = TRUE;
      if (hass)
        FPRINTF (out, "        S (spread of the split times): the Joint value tends to be too large\n"
                      "        (few stored genealogies fit a small S); the per-locus and All rows of\n"
                      "        S are mostly the prior, read S as an upper estimate\n");
      if (js->has_mcerr)
        jc_print_mcerr (world, js, out);
    }
  else if (jc_impossible > 0)
    FPRINTF (out, "(All) combined as the product of the per-locus marginal posteriors\n"
                  "      (the joint multi-locus combination was dropped: %li loci had no stored\n"
                  "      genealogy compatible with its parameters, see the warnings)\n",
             jc_impossible);
  else if (jc_nactive_cache > jc_maxactive ())
    FPRINTF (out, "(All) combined as the product of the per-locus marginal posteriors\n"
                  "      (the joint multi-locus combination is limited to %li free parameters,\n"
                  "      this model has %li; a restricted migration model can use it)\n",
             jc_maxactive (), jc_nactive_cache);
  else
    FPRINTF (out, "(All) combined as the product of the per-locus marginal posteriors\n"
                  "      (the joint multi-locus combination does not handle this model yet)\n");
  jc_print_splitcv (world, out);
}


/// the joint scaling factor of the multi-locus marginal likelihood and its
/// Monte Carlo error; FALSE when the joint combination did not provide it
boolean
jc_joint_scaling (world_fmt *world, double *logc, double *err)
{
  jc_store *js = (jc_store *) world->jointstats;
  if (js == NULL || !js->done || !js->has_logc)
    return FALSE;
  *logc = js->logc;
  *err = js->logc_err;
  return TRUE;
}

/// TRUE when the reweighting of the joint combination collapsed (median
/// bootstrap reweighting ESS below a quarter of the trace points): every
/// parameter is then flagged and the joint marginal likelihood is not
/// reliable; *ess gets that median ESS
boolean
jc_collapsed (world_fmt *world, double *ess)
{
  jc_store *js = (jc_store *) world->jointstats;
  if (js == NULL || !js->done || !js->has_boot)
    return FALSE;
  if (ess != NULL)
    *ess = js->boot_ess;
  return js->boot_ess < 0.25 * JC_BOOTK;
}

/// TRUE when the joint estimate of parameter p has a large genealogy-sampling error (bootstrap
/// guard) and the product of the marginals ("All") is used instead
boolean
jc_param_flagged (world_fmt *world, long p)
{
  jc_store *js = (jc_store *) world->jointstats;
  if (js == NULL || !js->done || !js->has_boot || p < 0 || p >= js->np)
    return FALSE;
  return js->boot_flag[p];
}

/// skyline=PARAM: quantiles (2.5%, median, 97.5%) of the multiplier of
/// matrix entry k in segment sg and of the segment value (parameter x
/// multiplier) from the joint multi-locus samples; FALSE when the joint
/// combination did not run for this model
boolean
jc_sky_quantiles (world_fmt *world, long k, long sg, double *qt, double *qv)
{
  jc_store *js = (jc_store *) world->jointstats;
  const double pq[3] = { 0.025, 0.5, 0.975 };
  long j, i, a;
  if (js == NULL || !js->done || js->nseg < 2 || sg < 1 || sg >= js->nseg)
    return FALSE;
  for (j = 0; j < js->nfree; j++)
    if (js->sky_k[j] == k)
      break;
  if (j == js->nfree)
    return FALSE;
  const long pt = js->numparam + j * (js->nseg - 1) + sg - 1;
  const long pb = world->bayes->map[k][1];
  /* a parameter or multiplier with a large genealogy-sampling error
     (Joint*): the caller uses the product of the per-locus posteriors */
  if (jc_param_flagged (world, pb) || jc_param_flagged (world, pt))
    return FALSE;
  double *t = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
  double *v = (double *) mycalloc ((size_t) js->nsweeps, sizeof (double));
  for (i = 0; i < js->nsweeps; i++)
    {
      t[i] = js->trace[pt * js->nsweeps + i];
      v[i] = t[i] * js->trace[pb * js->nsweeps + i];
    }
  qsort (t, (size_t) js->nsweeps, sizeof (double), jc_cmp_double);
  qsort (v, (size_t) js->nsweeps, sizeof (double), jc_cmp_double);
  for (a = 0; a < 3; a++)
    {
      i = (long) (pq[a] * (js->nsweeps - 1) + 0.5);
      qt[a] = t[i];
      qv[a] = v[i];
    }
  myfree (t);
  myfree (v);
  return TRUE;
}
