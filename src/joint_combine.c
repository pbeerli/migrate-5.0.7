/*------------------------------------------------------------------------
  Joint multi-locus combination of Bayesian posteriors (2026-10-01; ported
  from migrate-codex-7 7.0.67 to migrate 5.0.7 6.1.42: without skyline=PARAM
  multipliers, log-scale histogram bins, the markdown report and the
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
  Other models (Mittag-Leffler, growth with skyline segments) keep the
  marginal-histogram product.

  Storage follows the posterior samples (bayes-allfile): by default the
  statistics stay in memory; in the low-memory file mode (has_bayesmdimfile)
  every process instead writes "<outfile>.joint.<locus>.<replicate>", and
  at the end the process that evaluates a locus reads and thins them (needs
  a shared file system, like the posterior-sample files). In memory, each
  process keeps a thinned buffer per locus (world->jointstats, at most
  2*JC_MAXSAMPLES rows: when it fills, every other row is dropped and the
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
#ifdef MPI
#include "migrate_mpi.h"
extern const MPI_Datatype mpisizeof;
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SKYLINE_HISTMAX
#define SKYLINE_HISTMAX 8.0      /* |log tau| limit (no skyline=PARAM in 5.0.7) */
#define SKYPRIOR_RANDOMWALK 1
#define SKYPRIOR_LOGUNIFORM 2
#endif
#define JC_MAXSAMPLES 4000   /* per locus used in the combination (at 200 loci,
                                 1000 gave 95% coverage of 0.4 for an M, 4000 0.8) */
#define JC_GRID 600          /* points for a one-dimensional normalizer */
#define JC_TGRID 200         /* Theta points inside a growth integral */
#define JC_NG 129            /* growth grid */
#define JC_NM 33             /* divergence-time mean grid, local to each genealogy (quadratic) */
#define JC_DZ 8.0            /* local window: centre +- JC_DZ std */
#define JC_DOUT (-1e30)      /* divergence term outside a row's window */
#define JC_NS 17             /* divergence-time std grid (log-spaced cells) */
#define JC_BURNIN 2000       /* sweeps */
#define JC_SWEEPS 50000
#define JC_SSRUNGS 24        /* stepping stones for the joint marginal likelihood */
#define JC_SSBURN 300        /* sweeps per rung: burn-in, samples */
#define JC_SSN 1500
#define JC_NBLOCK 4          /* blocks of every locus' genealogies for the Monte Carlo error */

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

/* M or xNm: 5.0.7 has use-M only (S/M models are not handled here) */
static boolean
migration_is_M (world_fmt *world, long i)
{
  (void) i;
  return world->options->usem;
}

/* ------------------------------------------------------------- layout */

typedef struct
{
  long numpop, numpop2, np;
  long nrow;                /* doubles per genealogy */
  long off_c, off_K, off_S, off_m;
  long off_mu;              /* the locus rate of the genealogy (estimated rate: per sample) */
  long off_rep;             /* the replicate that sampled it (last entry of the row) */
  long nslot;               /* growth groups */
  long *pop_slot;           /* per population: growth group or -1 */
  long *slot_param;         /* per group: its parameter index */
  long *off_t;              /* per population: offset of T_i (growing only) */
  long *off_A;              /* per population: offset of A_i on the grid */
  double *ggrid;            /* per group: JC_NG growth values */
  long *rep;                /* per matrix entry: its parameter, or -1 (fixed) */
  double *fixed;            /* per matrix entry: value when fixed */
  boolean usem;             /* migration parameters are M (else xNm) */
  boolean *xnm;             /* per matrix entry: parameter is xNm (use-M=NO, S/M) */
  long nsplit;              /* divergence models */
  long *split_model;        /* per split: index into world->species_model */
  long *split_pmu, *split_psig; /* per split: parameter indices (psig -1: fixed) */
  long *split_ns;           /* per split: std grid points (1 when fixed) */
  long *split_rep;          /* per split: first split of its tied group (t/T), else itself */
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
} jc_layout;

/* phi index of tau_s (s >= 1) of free entry j */
#define JC_TAU(ly, j, s) ((ly)->np + (j) * ((ly)->nseg - 1) + (s) - 1)

static void
jc_layout_free (jc_layout *ly)
{
  myfree (ly->pop_slot);
  myfree (ly->slot_param);
  myfree (ly->off_t);
  myfree (ly->off_A);
  myfree (ly->ggrid);
  myfree (ly->rep);
  myfree (ly->xnm);
  myfree (ly->fixed);
  myfree (ly->split_model);
  myfree (ly->split_pmu);
  myfree (ly->split_psig);
  myfree (ly->split_ns);
  myfree (ly->split_rep);
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
  ly->nrow = 3 * numpop + (world->numpop2 - numpop);
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
  /* skyline=PARAM segments: none in 5.0.7 (the code stays general) */
  ly->nseg = 1;
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
  ly->ntau = ly->nfree * (ly->nseg - 1);
  ly->nphi = ly->np + ly->ntau;
  ly->off_mu = ly->nrow;
  ly->nrow += 1;
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

/* ------------------------------------------------------------- support */

/// TRUE when the model is one the joint combination handles
boolean
jc_supported (world_fmt *world)
{
  long i;
  if (!world->options->bayes_infer || world->has_mlalpha)
    return FALSE;
  /* a skyline with time-varying parameters is not handled here */
  if (world->timeelements > 2)
    return FALSE;
  if (world->has_speciation)
    for (i = 0; i < world->species_model_size; i++)
      {
        const species_fmt *sp = &world->species_model[i];
        if (!strchr ("dt", sp->type) || sp->paramindex_mu < 0 || sp->paramindex_mu >= world->numparam
            || world->bayes->map[sp->paramindex_mu][1] != sp->paramindex_mu)
          return FALSE;   /* (type is lower case: d/D -> 'd', tied t/T -> 't') */
      }

  for (i = 0; i < world->numpop2; i++)
    {
      const char c = world->options->custm2[i];
      /* '*', '0', constant, symmetric, mean (of M or of xNm), letter groups,
         divergence d/D; t/T are recorded as d with a tied species model */
      /* no numbers-of-migrants S/M: their parameter differs in 5.0.7 */
      if (!(strchr ("*0csmdDtT", c) || (c >= 'a' && c <= 'z')))
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

/* ------------------------------------------------------------- storage */

/* genealogies kept per locus: JC_MAXSAMPLES, fewer for wide rows (the
   divergence grids) so that a worker's memory stays bounded */
static long
jc_maxsamples (long nrow)
{
  long m = (long) (1000000 / (nrow > 0 ? nrow : 1));
  long cap = JC_MAXSAMPLES;
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
  double *trace;  /* np x JC_SWEEPS */
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
jc_append (jc_store *js, long locus, const double *row)
{
  jc_buffer *b = &js->b[locus];
  if (b->n >= b->cap)
    {
      b->cap = b->cap ? 2 * b->cap : 64;
      b->rows = (double *) myrealloc (b->rows, sizeof (double) * (size_t) (b->cap * js->nrow));
    }
  memcpy (b->rows + b->n * js->nrow, row, sizeof (double) * (size_t) js->nrow);
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
      const long kto = tl[i].lineages[sp->to];
      if (kto > 0)
        D += (double) kto * (*log_prob_wait_speciate) (t0, t1, mu, sigma, sp);
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
        if (tl[i].eventnode->type == 'd'
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
            if (tl[i].lineages[sp->to] > 0 && tl[i].age > sum)
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

/// Per-genealogy statistics of the current cold-chain genealogy (called
/// right after bayes_save()).
void
jc_record_sample (world_fmt *world)
{
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
      for (g = 0; g < JC_NG; g++)
        st[ly.off_A[pop] + g] = -HUGE_VAL;
  for (i = 1; i < T - 1; i++)
    {
      const double t0 = tl[i - 1].age, t1 = tl[i].age;
      const double dt = t1 - t0;
      const long *k = tl[i].lineages;
      const char type = tl[i].eventnode->type;
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
          if (ly.pop_slot[xp] >= 0)
            st[ly.off_t[xp]] += t1;
        }
      else if (type == 'm')
        {
          const long j = m2mmm (tl[i].eventnode->pop, tl[i].eventnode->actualpop, numpop);
          if (j >= numpop && j < world->numpop2)
            {
              st[ly.off_m + j - numpop] += 1.0;
              if (ly.nseg > 1)
                st[ly.off_sm + jc_seg (&ly, t1) * (world->numpop2 - numpop) + j - numpop] += 1.0;
            }
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
  if (inh != 1.0)
    for (pop = 0; pop < numpop; pop++)
      {
        st[ly.off_K + pop] /= inh;
        for (g = 0; g < ly.nseg && ly.nseg > 1; g++)
          st[ly.off_sK + g * numpop + pop] /= inh;
        if (ly.pop_slot[pop] >= 0)
          for (g = 0; g < JC_NG; g++)
            st[ly.off_A[pop] + g] -= log (inh);
      }
  if (jc_file_mode (world))
    {
      /* thin while writing: about 2 * jc_maxsamples() rows per file */
      static long jc_file_count = 0;
      long fstride = world->options->lsteps / (2 * jc_maxsamples (ly.nrow));
      if (fstride < 1)
        fstride = 1;
      if (jc_file == NULL || jc_file_locus != world->locus || jc_file_rep != world->rep)
        jc_file_count = 0;
      if ((jc_file_count++ % fstride) != 0)
        {
          myfree (st);
          jc_layout_free (&ly);
          return;
        }
      jc_record_divergence (world, &ly, st);
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
          jc_append (js, world->locus, st);
          if (b->n >= 2 * jc_maxsamples (js->nrow))
            {
              /* thin: keep every other row, record half as often */
              for (i = 0; i < b->n / 2; i++)
                memmove (b->rows + i * js->nrow, b->rows + 2 * i * js->nrow,
                         sizeof (double) * (size_t) js->nrow);
              b->n /= 2;
              b->stride *= 2;
            }
        }
    }
  myfree (st);
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
  js->b[locus].n = js->b[locus].cap = js->b[locus].seen = 0;
  js->b[locus].stride = 1;
}

/// MPI: pack this rank's statistics for one locus: {nrow, n, rows...}
long
jc_pack_buffer (MYREAL **buffer, world_fmt *world, long locus, long maxrep, long numpop)
{
  (void) maxrep;
  (void) numpop;
  jc_store *js = (jc_store *) world->jointstats;
  const long n = (js != NULL) ? js->b[locus].n : 0;
  const long nr = (js != NULL) ? js->nrow : 0;
  const long size = 2 + n * nr;
  *buffer = (MYREAL *) myrealloc (*buffer, sizeof (MYREAL) * (size_t) (size + 1));
  (*buffer)[0] = (MYREAL) nr;
  (*buffer)[1] = (MYREAL) n;
  if (n > 0)
    memcpy (*buffer + 2, js->b[locus].rows, sizeof (MYREAL) * (size_t) (n * nr));
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
  for (i = 0; i < n; i++)
    jc_append (js, locus, buffer + 2 + i * nr);
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
} jc_locus;

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
jc_logA (const jc_layout *ly, const double *st, long pop, double g)
{
  const double *gg = ly->ggrid + ly->pop_slot[pop] * JC_NG;
  const double *A = st + ly->off_A[pop];
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
   population does not grow); logA is log A_i(g), NAN means: look it up */
static double
jc_pop_term (const jc_layout *ly, const double *st, long pop, double theta, double g,
             double mu, double logA)
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
  return c * (LOG2 - jc_log_memo (mu * theta, &lx, &ll)) - A / (mu * theta) + tg;
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
  const double m = st[ly->off_m + e - ly->numpop];
  return (m > 0.0 ? m * jc_log_memo (x, &lx, &ll) : 0.0) - x * st[ly->off_S + to] / mu;
}

/* skyline=PARAM: population pop in segment sg with size v (= Theta tau) */
static double
jc_pop_seg_term (const jc_layout *ly, const double *st, long pop, long sg, double v, double mu)
{
  static double lx = -1.0, ll = 0.0;
  if (v <= 0.0)
    return -HUGE_VAL;
  const double c = st[ly->off_sc + sg * ly->numpop + pop];
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
  const double m = st[ly->off_sm + sg * (ly->numpop2 - ly->numpop) + e - ly->numpop];
  return (m > 0.0 ? m * jc_log_memo (x, &lx, &ll) : 0.0) - x * st[ly->off_sS + sg * ly->numpop + to] / mu;
}

/* term of a free skyline entry e (Theta or M) in segment sg at value x */
static double
jc_sky_seg_term (const jc_layout *ly, const double *st, long e, long sg, double x, double mu)
{
  return e < ly->numpop ? jc_pop_seg_term (ly, st, e, sg, x, mu) : jc_mig_seg_term (ly, st, e, sg, x, mu);
}

/* log prior of the multipliers (on log tau, as bayes_update_timeparam()),
   -HUGE_VAL outside: |log tau| <= SKYLINE_HISTMAX, the log-uniform bounds,
   and every segment value (parameter x tau) within the parameter's prior
   range (mlh_skyline_in_range()) */
static double
jc_sky_logprior (const jc_layout *ly, const double *phi)
{
  double lp = 0.0;
  long j, sg;
  for (j = 0; j < ly->nfree; j++)
    {
      const double base = phi[ly->rep[ly->sky_k[j]]];
      double prev = 0.0;
      for (sg = 1; sg < ly->nseg; sg++)
        {
          const double tau = phi[JC_TAU (ly, j, sg)];
          const double u = log (tau), v = base * tau;
          if (!(tau > 0.0) || fabs (u) > SKYLINE_HISTMAX || v < ly->sky_lo[j] || v > ly->sky_hi[j])
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
} jc_divquad;

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
  for (k = 0; k < ly->nsplit; k++)
    {
      const long ps = ly->split_psig[k];
      if (ps < 0)
        continue;
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
          ? scaling_prior (world, ps, q->sv[k * q->ns + g]) : -HUGE_VAL;
    }
}

static void
jc_divquad_free (jc_divquad *q)
{
  myfree (q->sv); myfree (q->sw); myfree (q->ps_s);
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
      v[g] = scaling_prior (q->world, pm, x[g]);
      for (m = 0; m < ly->nsplit; m++)
        if (ly->split_rep[m] == k)
          v[g] += jc_div_term (ly, st, m, x[g], sigma);
    }
  return jc_logtrapz (v, x, n);
}

/* log of int int prior(mu) prior(sigma) e^{D(mu, sigma)} for every split */
static double
jc_div_normalizer (const jc_layout *ly, const double *st, const jc_divquad *q)
{
  double lz = 0.0;
  long k, b;
  double *row = (double *) mycalloc ((size_t) q->ns, sizeof (double));
  for (k = 0; k < ly->nsplit; k++)
    {
      if (ly->split_rep[k] != k)   /* integrated with its group's first split */
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
  return lz;
}

/* skyline=PARAM: log R_e(x0) of free entry e on the log-spaced grid x of
   its parameter: the integral over the multipliers tau_1.. of their prior
   times the entry's terms in every segment, the base value x0 = x[i]:

     R_e(x0) = e^{f_0(x0)} int prod_{s>=1} K(u_{s-1}, u_s) e^{f_s(x0 e^{u_s})} du

   In w_s = log(x0 tau_s) the range rule (every segment value in the prior
   range) is the grid itself; the random walk on log tau is a random walk on
   w that starts at w_0 = log x0 (backward recursion with a banded Gaussian
   kernel, rows normalised: a std below the grid spacing pins the segments
   together); the log-uniform prior is independent per segment, w_s within
   w_0 + [log a, log b]. The hard limit |log tau| <= SKYLINE_HISTMAX is part
   of the log-uniform windows; for the random walk it is not part of R (its
   prior mass is negligible unless sigma is large) but the joint sampler
   enforces it. Constant factors of the prior cancel in p(G)/Z(G). */
static void
jc_sky_logR (const jc_layout *ly, const double *st, long e, double mu, const double *x, double *R)
{
  const long n = JC_TGRID, j = ly->sky_j[e];
  const double h = (log (x[n - 1]) - log (x[0])) / (double) (n - 1);
  double f[JC_TGRID], B[JC_TGRID], C[JC_TGRID];
  long i, q, sg;
  if (ly->sky_type[j] == SKYPRIOR_LOGUNIFORM)
    {
      const double la = log (ly->sky_a[j]) > -SKYLINE_HISTMAX ? log (ly->sky_a[j]) : -SKYLINE_HISTMAX;
      const double lb = log (ly->sky_b[j]) < SKYLINE_HISTMAX ? log (ly->sky_b[j]) : SKYLINE_HISTMAX;
      for (i = 0; i < n; i++)
        R[i] = jc_sky_seg_term (ly, st, e, 0, x[i], mu);
      for (sg = 1; sg < ly->nseg; sg++)
        {
          double fmax = -HUGE_VAL;
          for (i = 0; i < n; i++)
            {
              f[i] = jc_sky_seg_term (ly, st, e, sg, x[i], mu);
              if (f[i] > fmax)
                fmax = f[i];
            }
          C[0] = 0.0;   /* cumulative trapezoid of e^{f - fmax} in w */
          for (i = 1; i < n; i++)
            C[i] = C[i - 1] + 0.5 * h * (exp (f[i - 1] - fmax) + exp (f[i] - fmax));
          for (i = 0; i < n; i++)
            {
              /* window [i + la/h, i + lb/h] in grid units, clipped to the grid */
              double a = (double) i + la / h, b = (double) i + lb / h, ca, cb;
              if (a < 0.0) a = 0.0;
              if (b > (double) (n - 1)) b = (double) (n - 1);
              if (b <= a || fmax == -HUGE_VAL)
                {
                  R[i] = -HUGE_VAL;
                  continue;
                }
              q = (long) a; if (q > n - 2) q = n - 2;
              ca = C[q] + (a - q) * (C[q + 1] - C[q]);
              q = (long) b; if (q > n - 2) q = n - 2;
              cb = C[q] + (b - q) * (C[q + 1] - C[q]);
              R[i] += (cb > ca) ? fmax + log (cb - ca) : -HUGE_VAL;
            }
        }
      return;
    }
  /* random walk with std sigma on w */
  const double sigma = ly->sky_a[j];
  long band = (long) ceil (6.0 * sigma / h);
  if (band > n - 1)
    band = n - 1;
  double *kw = (double *) mycalloc ((size_t) (2 * band + 1), sizeof (double));
  double ksum = 0.0;
  for (q = -band; q <= band; q++)
    ksum += (kw[q + band] = exp (-0.5 * (q * h / sigma) * (q * h / sigma)));
  for (q = 0; q <= 2 * band; q++)
    kw[q] /= ksum;
  for (i = 0; i < n; i++)
    B[i] = jc_sky_seg_term (ly, st, e, ly->nseg - 1, x[i], mu);
  for (sg = ly->nseg - 2; sg >= 0; sg--)
    {
      double bmax = -HUGE_VAL;
      for (i = 0; i < n; i++)
        if (B[i] > bmax)
          bmax = B[i];
      for (i = 0; i < n; i++)
        C[i] = bmax == -HUGE_VAL ? 0.0 : exp (B[i] - bmax);
      for (i = 0; i < n; i++)
        {
          double sum = 0.0;
          const long lo = i - band < 0 ? 0 : i - band, hi = i + band > n - 1 ? n - 1 : i + band;
          for (q = lo; q <= hi; q++)
            sum += kw[q - i + band] * C[q];
          B[i] = jc_sky_seg_term (ly, st, e, sg, x[i], mu) + (sum > 0.0 ? bmax + log (sum) : -HUGE_VAL);
        }
    }
  memcpy (R, B, sizeof (double) * (size_t) n);
  myfree (kw);
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
    if (active[p] && jc_sky_param (ly, p))
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
    if (active[p] && !ly->xnm[p] && !jc_sky_param (ly, p))
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
  /* growth groups: int prior(g) prod_pop int prior(Theta) e^{f(Theta, g)} */
  for (s = 0; s < ly->nslot; s++)
    {
      const double *gg = ly->ggrid + s * JC_NG;
      for (g = 0; g < JC_NG; g++)
        {
          w[g] = gpri[s * JC_NG + g];
          for (pop = 0; pop < numpop; pop++)
            if (ly->pop_slot[pop] == s && active[pop])
              {
                const double Ag = st[ly->off_A[pop] + g];
                for (t = 0; t < JC_TGRID; t++)
                  u[t] = tpri_s[pop * JC_TGRID + t]
                    + jc_pop_term (ly, st, pop, tgrid_s[pop * JC_TGRID + t], gg[g], mu, Ag);
                w[g] += jc_logtrapz (u, tgrid_s + pop * JC_TGRID, JC_TGRID);
              }
        }
      lz += jc_logtrapz (w, gg, JC_NG);
    }
  lz += jc_div_normalizer (ly, st, dq);
  return lz;
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

static long
jc_read_locus_files (world_fmt *world, long locus, long nrow, jc_locus *L)
{
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
          keep = n < jc_maxsamples (nrow) ? n : jc_maxsamples (nrow);   /* rows of the full combination */
          L->st = (double *) mycalloc ((size_t) (keep * nrow), sizeof (double));
        }
      row = 0;
      for (r = 0; r < reps; r++)
        {
          char name[LINESIZE];
          FILE *f;
          jc_filename (world, locus, r, name);
          if ((f = fopen (name, "r")) == NULL)
            continue;
          for (;;)
            {
              for (i = 0; i < nrow; i++)
                if (fscanf (f, "%lf", &buf[i]) != 1)
                  break;
              if (i < nrow)
                break;
              if (pass == 0)
                n++;
              else if (mi < keep && row == (mi * n) / keep)
                {   /* a row of the full combination; a block keeps its own */
                  if (jc_row_in_part (buf, nrow, row, ra, rlen))
                    memcpy (L->st + (kept++) * nrow, buf, sizeof (double) * (size_t) nrow);
                  mi++;
                }
              row++;
            }
          fclose (f);
        }
    }
  myfree (buf);
  if (kept == 0 && L->st != NULL)
    {
      myfree (L->st);
      L->st = NULL;
    }
  L->total = n;
  L->n = kept;
  L->mu = world->options->mu_rates[locus];
  return kept;
}

static long
jc_read_locus (world_fmt *world, long locus, long nrow, jc_locus *L)
{
  jc_store *js = (jc_store *) world->jointstats;
  long i;
  if (jc_file_mode (world))
    return jc_read_locus_files (world, locus, nrow, L);
  if (js == NULL || js->nrow != nrow || js->b[locus].n == 0)
    return 0;
  const long n = js->b[locus].n;
  /* the rows of the full combination (evenly spaced); a block uses those
     of them that belong to it, so that the blocks partition them */
  const long kfull = n < jc_maxsamples (nrow) ? n : jc_maxsamples (nrow);
  long ra, rlen, keep = 0;
  jc_part_range (n, &ra, &rlen);
  L->st = (double *) mycalloc ((size_t) (kfull * nrow), sizeof (double));
  for (i = 0; i < kfull; i++)
    {
      const long row = (i * n) / kfull;
      if (jc_row_in_part (js->b[locus].rows + row * nrow, nrow, row, ra, rlen))
        memcpy (L->st + (keep++) * nrow, js->b[locus].rows + row * nrow, sizeof (double) * (size_t) nrow);
    }
  if (keep < 1)
    {
      myfree (L->st);
      L->st = NULL;
      return 0;
    }
  L->n = keep;
  L->total = n;
  L->mu = world->options->mu_rates[locus];
  return keep;
}

/* log p(G|phi) up to terms that do not depend on any estimated parameter */
static double
jc_full_term (const jc_layout *ly, const boolean *active, const double *st, const double *phi, double mu)
{
  double v = 0.0;
  long pop, e;
  long sg;
  for (pop = 0; pop < ly->numpop; pop++)
    if (ly->rep[pop] >= 0 && active[ly->rep[pop]])
      {
        const long s = ly->pop_slot[pop];
        const long j = ly->sky_j[pop];
        if (j >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            v += jc_pop_seg_term (ly, st, pop, sg, phi[ly->rep[pop]] * (sg ? phi[JC_TAU (ly, j, sg)] : 1.0), mu);
        else
          v += jc_pop_term (ly, st, pop, phi[ly->rep[pop]], s >= 0 ? phi[ly->slot_param[s]] : 0.0, mu, NAN);
      }
  for (e = ly->numpop; e < ly->numpop2; e++)
    if (ly->rep[e] >= 0 && active[ly->rep[e]])
      {
        const long j = ly->sky_j[e];
        if (j >= 0)
          for (sg = 0; sg < ly->nseg; sg++)
            v += jc_mig_seg_term (ly, st, e, sg, phi[ly->rep[e]] * (sg ? phi[JC_TAU (ly, j, sg)] : 1.0), mu);
        else
          v += jc_mig_term (ly, st, e, jc_rate (ly, phi, e), mu);
      }
  for (e = 0; e < ly->nsplit; e++)
    v += jc_div_term (ly, st, e, phi[ly->split_pmu[e]],
                      ly->split_psig[e] >= 0 ? phi[ly->split_psig[e]] : ly->split_sfix[e]);
  return v;
}


/* The terms of jc_full_term() one by one ("units", in its order): a
   proposal of one parameter recomputes only the units that depend on it
   and sums all units again in the same order, so the sum is the same
   number jc_full_term() gives (no old + new - old updates: a divergence
   term at JC_DOUT would swallow the others) */
enum { JC_U_POP, JC_U_POPSEG, JC_U_MIG, JC_U_MIGSEG, JC_U_SPLIT };
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
          break;
        case JC_U_POPSEG:
        case JC_U_MIGSEG:
          jc_units_dep (U, ly->nphi, ly->rep[i], u);
          if (U->sg[u] > 0)
            jc_units_dep (U, ly->nphi, JC_TAU (ly, ly->sky_j[i], U->sg[u]), u);
          break;
        case JC_U_MIG:
          jc_units_dep (U, ly->nphi, ly->rep[i], u);
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
jc_unit_term (const jc_layout *ly, const jc_units *U, long u, const double *st, const double *phi, double mu)
{
  const long i = U->idx[u], sg = U->sg[u];
  switch (U->kind[u])
    {
    case JC_U_POP:
      {
        const long s = ly->pop_slot[i];
        return jc_pop_term (ly, st, i, phi[ly->rep[i]], s >= 0 ? phi[ly->slot_param[s]] : 0.0, mu, NAN);
      }
    case JC_U_POPSEG:
      return jc_pop_seg_term (ly, st, i, sg, phi[ly->rep[i]] * (sg ? phi[JC_TAU (ly, ly->sky_j[i], sg)] : 1.0), mu);
    case JC_U_MIG:
      return jc_mig_term (ly, st, i, jc_rate (ly, phi, i), mu);
    case JC_U_MIGSEG:
      return jc_mig_seg_term (ly, st, i, sg, phi[ly->rep[i]] * (sg ? phi[JC_TAU (ly, ly->sky_j[i], sg)] : 1.0), mu);
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
  long pending;      /* the parameter of the pending proposal */
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
      myfree (J->L[l].cur);
      myfree (J->L[l].logz);
      myfree (J->L[l].unit);
    }
  myfree (J->L);
  myfree (J->newunit);
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
  for (s = 0; s < ly->nsplit; s++)
    {
      active[ly->split_pmu[s]] = world->bayes->maxparam[ly->split_pmu[s]] > world->bayes->minparam[ly->split_pmu[s]];
      if (ly->split_psig[s] >= 0)
        active[ly->split_psig[s]] = TRUE;
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
      if (jc_read_locus (world, locus, ly->nrow, &J->L[J->nl]) > 0)
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
          tpri[p * JC_GRID + g] = scaling_prior (world, p, tgrid[p * JC_GRID + g]);
        }
      for (g = 0; g < JC_TGRID; g++)
          {
            tgrid_s[p * JC_TGRID + g] = g == 0 ? lo : (g == JC_TGRID - 1 ? hi
              : exp (log (lo) + (log (hi) - log (lo)) * (double) g / (JC_TGRID - 1)));
            tpri_s[p * JC_TGRID + g] = scaling_prior (world, p, tgrid_s[p * JC_TGRID + g]);
          }
    }
  for (s = 0; s < ly->nslot; s++)
    for (g = 0; g < JC_NG; g++)
      gpri[s * JC_NG + g] = scaling_prior (world, ly->slot_param[s], ly->ggrid[s * JC_NG + g]);
  jc_divquad dq;
  jc_divquad_make (world, ly, &dq);
  /* per-sample normalizers and current terms */
  double total = 0.0;
  J->lterm = (double *) mycalloc ((size_t) J->nl, sizeof (double));
  J->newlterm = (double *) mycalloc ((size_t) J->nl, sizeof (double));
  J->tmin = J->tmax = J->L[0].total;
  jc_units_make (ly, J->active, &J->U);
  J->pending = -1;
  const long nu = J->U.n;
  for (l = 0; l < J->nl; l++)
    {
      long sidx, u;
      jc_locus *Ll = &J->L[l];
      Ll->cur = (double *) mycalloc ((size_t) Ll->n, sizeof (double));
      Ll->logz = (double *) mycalloc ((size_t) Ll->n, sizeof (double));
      Ll->unit = (double *) mycalloc ((size_t) (Ll->n * (nu > 0 ? nu : 1)), sizeof (double));
      for (sidx = 0; sidx < Ll->n; sidx++)
        {
          const double *st = Ll->st + sidx * ly->nrow;
          double *un = Ll->unit + sidx * nu;
          double v = 0.0;
          Ll->logz[sidx] = jc_log_normalizer (ly, J->active, st, st[ly->off_mu], tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
          for (u = 0; u < nu; u++)
            {
              un[u] = jc_unit_term (ly, &J->U, u, st, phi, st[ly->off_mu]);
              v += un[u];
            }
          Ll->cur[sidx] = v - Ll->logz[sidx];
        }
      J->lterm[l] = jc_logsumexp (Ll->cur, Ll->n);
      total += J->lterm[l];
      if (Ll->n > J->maxn)
        J->maxn = Ll->n;
      if (Ll->total < J->tmin)
        J->tmin = Ll->total;
      if (Ll->total > J->tmax)
        J->tmax = Ll->total;
    }
  J->newcur = (double *) mycalloc ((size_t) (J->nl * J->maxn), sizeof (double));
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
        for (g = 0; g < JC_NG; g++)
          z0[ly->off_A[pop] + g] = -HUGE_VAL;
    z0[ly->off_mu] = 1.0;
    J->lognz = jc_log_normalizer (ly, J->active, z0, 1.0, tgrid, tpri, tgrid_s, tpri_s, gpri, &dq);
    myfree (z0);
  }
  myfree (tgrid); myfree (tpri); myfree (tgrid_s); myfree (tpri_s); myfree (gpri);
  jc_divquad_free (&dq);
  return total;
}

/* proposal phi[p] -> x: local change of sum log f_l (state kept pending) */
static double
jc_local_propose (jc_local *J, double *phi, long p, double x)
{
  double delta = 0.0;
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
      const double keep = phi[p];
      phi[p] = x;
      for (sidx = 0; sidx < Ll->n; sidx++)
        {
          const double *st = Ll->st + sidx * J->ly.nrow;
          const double *un = Ll->unit + sidx * nu;
          double *nws = nw + sidx * nd;
          double v = 0.0;
          for (k = 0; k < nd; k++)
            nws[k] = jc_unit_term (&J->ly, &J->U, dep[k], st, phi, st[J->ly.off_mu]);
          for (u = 0; u < nu; u++)
            v += map[u] >= 0 ? nws[map[u]] : un[u];
          nc[sidx] = v - Ll->logz[sidx];
        }
      phi[p] = keep;
      J->newlterm[l] = jc_logsumexp (nc, Ll->n);
      delta += J->newlterm[l] - J->lterm[l];
    }
  myfree (map);
  return delta;
}

static void
jc_local_accept (jc_local *J)
{
  const long p = J->pending;
  const long nu = J->U.n;
  const long nd = (p >= 0 && p < J->ly.nphi) ? J->U.ndep[p] : 0;
  long l, sidx, k;
  for (l = 0; l < J->nl; l++)
    {
      jc_locus *Ll = &J->L[l];
      const double *nw = J->newunit + l * J->maxn * (nd > 0 ? nd : 1);
      memcpy (Ll->cur, J->newcur + l * J->maxn, sizeof (double) * (size_t) Ll->n);
      J->lterm[l] = J->newlterm[l];
      for (sidx = 0; sidx < Ll->n; sidx++)
        for (k = 0; k < nd; k++)
          Ll->unit[sidx * nu + J->U.dep[p][k]] = nw[sidx * nd + k];
    }
  J->pending = -1;
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
            t[s * K + k] = jc_full_term (ly, J->active, st, pts + k * ly->nphi, st[ly->off_mu]) - Ll->logz[s];
        }
      for (k = 0; k < K; k++)
        {
          for (s = 0; s < n; s++)
            col[s] = t[s * K + k];
          orig[k] = jc_logsumexp (col, n) - log ((double) n);
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
              D[b * K + k] += jc_logsumexp (col, m) - log ((double) total) - orig[k];
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
      pending_p = (long) pkt[2];
      pending_x = pkt[3];
      double d = jc_local_propose (&J, phi, pending_p, pending_x), dd = 0.0;
      MPI_Reduce (&d, &dd, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
  jc_local_free (&J);
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
static void
jc_chain (world_fmt *world, const jc_layout *ly, const boolean *active, jc_eval *E, double *phi,
          double *step, double beta, long burn, long sweeps, double *trace, double *logftr)
{
  const long np = ly->np, nphi = ly->nphi;
  const long npos = world->numparamcumvec[SPLITSTDPRIOR];
  long *acc = (long *) mycalloc ((size_t) nphi, sizeof (long));
  long *tri = (long *) mycalloc ((size_t) nphi, sizeof (long));
  long p, sweep;
  double skyprior = jc_sky_logprior (ly, phi);
  for (sweep = 0; sweep < burn + sweeps; sweep++)
    {
      for (p = 0; p < nphi; p++)
        {
          if (!active[p])
            continue;
          /* a multiplier: symmetric in log tau, its prior is on log tau */
          const boolean tau = p >= np;
          const double lo = tau ? exp (-SKYLINE_HISTMAX) : world->bayes->minparam[p];
          const double hi = tau ? exp (SKYLINE_HISTMAX) : world->bayes->maxparam[p];
          const boolean mult = tau || p < npos;
          const double x = mult ? phi[p] * exp (step[p] * (UNIF_RANDUM () - 0.5) * 2.0)
                                : phi[p] + step[p] * (UNIF_RANDUM () - 0.5) * 2.0;
          double newsky = skyprior;
          tri[p]++;
          if (x >= lo && x <= hi && (!mult || x > 0.0) && ly->nfree > 0)
            {
              const double keep = phi[p];
              phi[p] = x;
              newsky = jc_sky_logprior (ly, phi);
              phi[p] = keep;
            }
          if (x >= lo && x <= hi && (!mult || x > 0.0) && newsky > -HUGE_VAL)
            {
              double delta = (tau ? 0.0 : scaling_prior (world, p, x) - scaling_prior (world, p, phi[p])
                              + (mult ? log (x / phi[p]) : 0.0))
                + newsky - skyprior;
              const double d = jc_eval_step (E, phi, E->accept_prev, p, x);
              delta += beta * d;
              E->accept_prev = log (UNIF_RANDUM ()) < delta;
              if (E->accept_prev)
                {
                  acc[p]++;
                  phi[p] = x;
                  skyprior = newsky;
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
      if (sweep >= burn)
        {
          if (trace != NULL)
            for (p = 0; p < nphi; p++)
              trace[p * sweeps + (sweep - burn)] = phi[p];
          if (logftr != NULL)
            logftr[sweep - burn] = E->logF;
        }
    }
  myfree (acc);
  myfree (tri);
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
        val[k] = js->trace[p * JC_SWEEPS + (long) ((k + 0.5) * JC_SWEEPS / K)];
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
      for (i = 0; i < JC_SWEEPS; i++)
        mean += js->trace[p * JC_SWEEPS + i] / JC_SWEEPS;
      for (i = 0; i < JC_SWEEPS; i++)
        var += (js->trace[p * JC_SWEEPS + i] - mean) * (js->trace[p * JC_SWEEPS + i] - mean) / JC_SWEEPS;
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
        long part, long sweeps, double *trace, long *nl, long *tmin, long *tmax, double *rows,
        double *bootD)
{
  const long nphi = ly->nphi;
  double *phi = (double *) mycalloc ((size_t) nphi, sizeof (double));
  long p, l;
  jc_eval E;
  memcpy (phi, phi0, sizeof (double) * (size_t) nphi);
  jc_eval_start (world, &E, phi, part);
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
  jc_chain (world, ly, active, &E, phi, step, 1.0, JC_BURNIN, sweeps, trace, NULL);
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
          for (t = 0; t < 20; t++)
            {
              const long it = (sweeps / 20) * t + sweeps / 40;
              double mx = -HUGE_VAL, s1 = 0.0, s2 = 0.0;
              for (p = 0; p < nphi; p++)
                x[p] = trace[p * sweeps + it];
              for (sidx = 0; sidx < J->L[l].n; sidx++)
                {
                  const double *st = J->L[l].st + sidx * ly->nrow;
                  w[sidx] = jc_full_term (&J->ly, J->active, st, x, st[J->ly.off_mu]) - J->L[l].logz[sidx];
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
  jc_chain (world, ly, active, &E, phi, step, 1.0, JC_BURNIN, 0, NULL, NULL);
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
  js->trace = (double *) mycalloc ((size_t) (nphi * JC_SWEEPS), sizeof (double));
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
  if (!jc_run (world, &ly, active, phi, 0, JC_SWEEPS, js->trace, &nl, &tmin, &tmax, &rows_all, bootD))
    {
      myfree (phi);
      myfree (bootD);
      jc_layout_free (&ly);
      return FALSE;
    }
  jc_boot_errors (js, bootD);
  myfree (bootD);
  /* Monte Carlo error: the same combination on each block of every locus'
     genealogies (its replicates, or consecutive stretches of its chain);
     sd(block medians) scaled by the share of genealogies a block holds is
     about the standard error of the median from all */
  {
    const long nb = jc_nblocks (world);
    const long hs = JC_SWEEPS / nb > 5000 ? JC_SWEEPS / nb : 5000;
    double *th = (double *) mycalloc ((size_t) (nphi * hs), sizeof (double));
    double *srt = (double *) mycalloc ((size_t) JC_SWEEPS, sizeof (double));
    long part, n1, t1, t2;
    js->mcerr = (double *) mycalloc ((size_t) nphi, sizeof (double));
    js->mchalf = (double *) mycalloc ((size_t) (nb * nphi), sizeof (double));
    js->has_mcerr = TRUE;
    js->nblock = nb;
    js->byrep = (world->options->replicate ? world->options->replicatenum : 1) >= 2;
    for (part = 1; part <= nb && js->has_mcerr; part++)
      {
        if (!jc_run (world, &ly, active, phi, part, hs, th, &n1, &t1, &t2, &rb, NULL) || n1 < 2)
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
  js->has_logc = jc_ss (world, &ly, active, phi, &js->logc, &js->logc_err);
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
              memcpy (js->trace + pairs[t][0] * JC_SWEEPS, js->trace + pairs[t][1] * JC_SWEEPS,
                      sizeof (double) * (size_t) JC_SWEEPS);
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
  for (i = 0; i < JC_SWEEPS; i++)
    {
      if (x[i] <= 0.0)
        return;
      mean += x[i] / JC_SWEEPS;
    }
  for (i = 0; i < JC_SWEEPS; i++)
    var += (x[i] - mean) * (x[i] - mean) / JC_SWEEPS;
  if (var <= 0.0)
    return;
  c = sqrt (1.0 + e * e / var);
  {
    double *srt = (double *) mycalloc ((size_t) JC_SWEEPS, sizeof (double));
    memcpy (srt, x, sizeof (double) * (size_t) JC_SWEEPS);
    qsort (srt, (size_t) JC_SWEEPS, sizeof (double), jc_cmp_double);
    med = srt[JC_SWEEPS / 2];
    myfree (srt);
  }
  const double lmed = log (med);
  const double plo = world->bayes->minparam[p], phi = world->bayes->maxparam[p];
  for (i = 0; i < JC_SWEEPS; i++)
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
          const double *x0 = js->trace + pa * JC_SWEEPS;
          double *x = (double *) mycalloc ((size_t) JC_SWEEPS, sizeof (double));
          double mean = 0.0;
          memcpy (x, x0, sizeof (double) * (size_t) JC_SWEEPS);
          jc_widen (world, js, pa, x);
          memset (hist->results + off, 0, sizeof (double) * (size_t) nb);
          for (i = 0; i < JC_SWEEPS; i++)
            {
              long b = (long) floor ((x[i] - lo) / world->bayes->deltahist[pa]);
              if (b < 0)
                b = 0;
              if (b >= nb)
                b = nb - 1;
              hist->results[off + b] += 1.0 / (double) JC_SWEEPS;
              mean += x[i];
            }
          hist->means[pa] = mean / (double) JC_SWEEPS;
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
  srt = (double *) mycalloc ((size_t) JC_SWEEPS, sizeof (double));
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
      const double *tr = js->trace + p * JC_SWEEPS;
      double mean = 0.0, var = 0.0;
      for (i = 0; i < JC_SWEEPS; i++)
        mean += tr[i] / JC_SWEEPS;
      for (i = 0; i < JC_SWEEPS; i++)
        var += (tr[i] - mean) * (tr[i] - mean) / JC_SWEEPS;
      memcpy (srt, tr, sizeof (double) * (size_t) JC_SWEEPS);
      qsort (srt, (size_t) JC_SWEEPS, sizeof (double), jc_cmp_double);
      r->median = srt[JC_SWEEPS / 2];
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
                  "below %d, the Joint estimate has a large genealogy-sampling error (*): more\n"
                  "genealogies per locus (longer chains or replicates) are needed to confirm it.\n\n",
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
    FPRINTF (out, "(*) %ld Joint estimate%s with a large genealogy-sampling error (Joint* rows)\n",
             nf, nf > 1 ? "s" : "");
  myfree (rows);
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
                  "        other tables use it\n",
             js->nloci_used, js->tmin, js->tmax, jc_maxsamples (js->nrow));
  if (js != NULL && js->done && js->has_mcerr)
    jc_print_mcerr (world, js, out);
  else
    FPRINTF (out, "(All) combined as the product of the per-locus marginal posteriors\n"
                  "      (the joint multi-locus combination does not handle this model yet)\n");
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
