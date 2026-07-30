/*------------------------------------------------------
 Joint local-M + windowed migration-history update

 Bayesian   R O U T I N E S

 send questions concerning this software to:
 Peter Beerli
 beerli@fsu.edu

 Copyright 2026 Peter Beerli

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies
 or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
/*! \file windowupdate.c

 Joint local-M + windowed migration-history MCMC move.

 bayes_update() (bayes.c) moves M alone at fixed genealogy G; tree_update()
 (mcmc1.c/speciate.c) moves G alone at fixed (Theta,M) -- a Gibbs alternation
 (see the comment above scaler_update() in bayes.c). scaler_update() already
 addresses one mixing bottleneck this produces (the Theta*M ridge, via an
 exact joint rescaling). This file addresses a DIFFERENT one: probg_treetimes()
 shows p(M|G) is Gamma(n+1, L) with n the discrete migration-event count on G,
 so M's conditional width is ~1/sqrt(n) -- at high true migration rates the
 tree accumulates enough events that M becomes pinned to whatever the CURRENT
 count implies, and can only drift as fast as tree_update()'s single-branch
 random walk moves that count.

 The fix prototyped in Python before this port (full derivation and
 acceptance-rate validation in this repo's project memory under
 "TPSC-as-proposal", not re-derived here): propose M' by a small LOCAL step
 (not an independence sampler -- see below for why), then jointly redraw the
 migration history on a small, randomly chosen WINDOW of branches, leaving
 every other branch's history untouched. Both pieces of that design are load
 -bearing, not incidental:
   - windowing (vs. redrawing the whole tree) keeps the move's acceptance
     rate viable -- a naive whole-tree joint redraw collapses to ~0
     acceptance by 6-8 tips (already documented as a dead end in this repo's
     tests/README.md, independently discovered in the Python prototype too);
   - a LOCAL M step (vs. an independence sampler drawing M' from anywhere in
     the prior) is required once the window is narrow: branches OUTSIDE the
     window keep their existing migration-event counts, tuned to the OLD M,
     and a large M jump makes those untouched counts a poor fit to the NEW M
     -- this is invisible when (as in scaler_update, or an earlier
     whole-tree-redraw prototype) every branch gets refreshed, but becomes
     the dominant rejection source the moment some branches are left alone.

 This file currently implements window SELECTION and TOUCHED-BRANCH
 identification only (project task: "window selection and touched-branch
 identification"). The statistical core -- the CTMC bridge sampler, the
 mean-field forward pass used to choose each window node's new population,
 and the exact Hastings ratio -- is not yet implemented; see
 windowed_joint_update()'s TODO. Kept in its own file rather than folded
 into bayes.c because that core is substantial new numerical code with no
 existing counterpart elsewhere in migrate (confirmed by survey: migrate has
 no endpoint-conditioned CTMC bridge sampler and no mean-field message
 passing anywhere else in the codebase), not an extension of anything
 already in bayes.c.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bayes.h"
#include "definitions.h"
#include "migration.h"
#include "random.h"
#include "sighandler.h"
#include "tools.h"
#include "tree.h"
#include "world.h"

/* erase_migr_nodes() (mcmc2.c) strips all migration/divergence nodes off a
   branch, re-splicing its two logical endpoints directly. Not exposed in
   any header (only forward-declared inside mcmc2.c itself), but it has
   external linkage there, so it can be declared and used here too. */
extern long erase_migr_nodes (world_fmt * world, node * up);
#include "windowupdate.h"

/* Logical child of a coalescent node's non-top nodelet, skipping over any
   migration/divergence nodes spliced onto that branch -- the downward
   counterpart of crawlback() (tools.c), which does the same walking toward
   the ancestor. There is no existing "walk toward the child, skipping
   migration nodes" utility elsewhere in migrate (checked): every existing
   user of the migration-node chain either walks upward (crawlback(),
   erase_migr_nodes()) or is the chain-building code itself
   (add_migration()/insert_migr_node()), which is exactly what this mirrors
   in reverse. add_migration() sets `tmp->next->back = p` for the
   CHILD-facing nodelet of a newly spliced migration node `tmp`, which is
   the invariant this relies on: `theNode->next->back` is always the way
   further down, `theNode->back` the way further up. */
static node *
crawlforward (const node * theNode)
{
  node *tmp = theNode->back;
  while (tmp->type == 'm' || tmp->type == 'd')
    {
      tmp = tmp->next->back;
    }
  return tmp;
}

/* Pick up to window_size DISTINCT internal coalescent nodes ('i' type) as
   the move's window. world->nodep[0 .. sumtips-1] are the tips,
   world->nodep[sumtips .. 2*sumtips-2] are the sumtips-1 internal
   coalescent nodes (allocateinterior(), tree.c) -- the root itself is a
   separate, distinctly-allocated 'r' sentinel (world->root), not part of
   this range, so it is never picked (there is no migration history "above"
   the root to redraw anyway). Simple rejection sampling: window sizes used
   in practice (a handful of branches) make collisions rare enough that this
   is not worth a reservoir sampler.

   Returns the number of DISTINCT nodes actually placed in window[] --
   window_size itself if the tree has at least that many internal nodes,
   fewer only for a tree too small to have that many (sumtips-1 <
   window_size), never more. */
long
window_pick (world_fmt * world, node ** window, long window_size)
{
  const long ninternal = world->sumtips - 1;
  const long lo = world->sumtips;
  const long hi = 2 * world->sumtips - 2;      /* inclusive */
  long want = window_size;
  long got = 0;
  long tries;
  long maxtries;

  if (want > ninternal)
    want = ninternal;
  if (want <= 0)
    return 0;

  /* generous but bounded retry budget: with `want` << ninternal (the
     intended regime) a handful of tries per slot suffices; falling back to
     a full linear scan below guarantees termination even if `want` is close
     to ninternal, where rejection sampling alone could stall. */
  maxtries = 20 * want + 20;
  for (tries = 0; got < want && tries < maxtries; tries++)
    {
      node *cand = world->nodep[RANDINT (lo, hi)];
      boolean dup = FALSE;
      long j;
      for (j = 0; j < got; j++)
        {
          if (window[j] == cand)
            {
              dup = TRUE;
              break;
            }
        }
      if (!dup)
        {
          window[got++] = cand;
        }
    }
  if (got < want)
    {
      /* rejection sampling stalled (want close to ninternal) -- finish
         deterministically by scanning the remaining candidates in order. */
      long idx;
      for (idx = lo; idx <= hi && got < want; idx++)
        {
          node *cand = world->nodep[idx];
          boolean dup = FALSE;
          long j;
          for (j = 0; j < got; j++)
            {
              if (window[j] == cand)
                {
                  dup = TRUE;
                  break;
                }
            }
          if (!dup)
            {
              window[got++] = cand;
            }
        }
    }
  return got;
}

/* Collect the distinct branches touched by a window: for each window node
   v, the branch ABOVE v (v itself is that branch's young end) and the two
   branches BELOW v leading to its logical children (found via
   crawlforward() on v's two non-top nodelets, skipping any migration nodes
   already on those branches). Branches are identified by their young-end
   node pointer, matching insert_migr_node()'s/erase_migr_nodes()'s `up`
   argument convention.

   Deduplicates: if a window node's parent (or child) is ALSO a window node,
   the branch between them would otherwise be collected twice (once as the
   first node's "above" branch, once as the second node's "below" branch).
   touched_alloc must be at least 3*window_n (the maximum possible before
   dedup); returns the actual distinct count. */
long
window_touched_branches (node ** window, long window_n,
                          node ** touched, long touched_alloc)
{
  long n = 0;
  long i;

  for (i = 0; i < window_n; i++)
    {
      node *v = window[i];
      node *cands[3];
      long k;

      cands[0] = v;                        /* branch above v */
      cands[1] = crawlforward (v->next);           /* branch to child 1 */
      cands[2] = crawlforward (v->next->next);     /* branch to child 2 */

      for (k = 0; k < 3; k++)
        {
          boolean dup = FALSE;
          long j;
          for (j = 0; j < n; j++)
            {
              if (touched[j] == cands[k])
                {
                  dup = TRUE;
                  break;
                }
            }
          if (!dup)
            {
              if (n >= touched_alloc)
                error ("window_touched_branches: touched[] too small");
              touched[n++] = cands[k];
            }
        }
    }
  return n;
}

/*========================================================================
  Migration-rate generator matrix (Q) and the local-M proposal step.
  ========================================================================*/

/* Builds the full numpop x numpop migration generator matrix (row-major,
   d = world->numpop) from world->param0, matching EXACTLY how
   probg_treetimes() scores a migration event (bayes.c ~line 944-954): the
   rate for a from->to event is param0[m2mmm(from,to,d)] if world->options
   ->usem, else param0[m2mmm(from,to,d)]/param0[to] (the 4Nm-style
   convention, dividing by the destination population's Theta). This
   distinction matters: using the wrong convention here would make the
   bridge sampler/Hastings ratio disagree with what probg_treetimes()
   actually scores, silently biasing the chain.

   `which`/`which_value` let the caller substitute ONE migration parameter's
   value (the one being locally perturbed by this move) without mutating
   world->param0 itself -- pass which=-1 to build Q from param0 as-is. */
static void
build_Q (world_fmt * world, long which, MYREAL which_value, MYREAL * Q)
{
  const long d = world->numpop;
  const boolean usem = world->options->usem;
  const MYREAL *param0 = world->param0;
  long from, to;

  for (from = 0; from < d; from++)
    {
      MYREAL rowsum = 0.0;
      for (to = 0; to < d; to++)
        {
          MYREAL rate;
          if (from == to)
            continue;
          long idx = m2mmm (from, to, d);
          MYREAL mval = (idx == which) ? which_value : param0[idx];
          rate = usem ? mval : mval / param0[to];
          Q[from * d + to] = rate;
          rowsum += rate;
        }
      Q[from * d + from] = -rowsum;
    }
}

/* Local (small-step) proposal for ONE migration parameter, mirroring
   bayes_update()'s parameter-selection convention (uniform pick over the
   migration block [numpop,numpop2), filtered through shortcut() to skip
   constrained/tied entries exactly as bayes_update does) but with a
   DIFFERENT step shape: a step uniform in log-space (log(Mp/M) ~
   Uniform(-window_delta,window_delta)), matching scaler_update's own
   log-uniform convention (bayes.c ~line 2167-2170) rather than a Gaussian
   random walk, so no new RNG primitive is needed beyond UNIF_RANDUM().

   Why local at all, and not scaler_update's/bayes_update's independence
   sampler drawing from anywhere in the prior: this move leaves branches
   OUTSIDE the window untouched, and a large M jump would make their
   existing (old-M-tuned) migration-event counts a poor fit to the new M --
   invisible when every branch gets refreshed (as in scaler_update, or an
   independence-sampler variant that failed for exactly this reason in the
   Python prototype), but the dominant rejection source once some branches
   are deliberately left alone. See file header.

   Returns FALSE if the drawn Mp falls outside [bayes->minparam[which],
   bayes->maxparam[which]] (an immediate reject, matching the rest of
   migrate's Metropolis kernels -- out-of-support proposals have zero prior
   density, equivalent to always-reject); TRUE otherwise, with *which_out,
   *Mold_out, *Mnew_out, *logpriorjac_out filled in. logpriorjac_out is the
   COMPLETE prior+Jacobian contribution to the log Hastings ratio: migrate's
   own log_prior_ratio[which] (bayes.h/bayes.c, correct for whatever prior
   kind -- uniform/exp/gamma/... -- this parameter actually has, reused
   rather than re-derived) plus the log(Mp)-log(M) Jacobian for the
   log-uniform step (symmetric in log M, so its own q-ratio is 1, but the
   prior is stated in M-space and this move works in log-M-space -- same
   derivation as the earlier local-M Python prototype's
   move_windowed_bridge_localM). */
static boolean
local_M_propose (world_fmt * world, long *which_out, MYREAL * Mold_out,
                  MYREAL * Mnew_out, MYREAL * logpriorjac_out)
{
  const long numpop = world->numpop;
  const long numpop2 = world->numpop2;
  const MYREAL delta = world->options->window_delta;
  bayes_fmt *bayes = world->bayes;
  long which;
  long w;
  MYREAL Mold, Mnew, logratio;

  if (numpop2 <= numpop)
    return FALSE;               /* no migration parameters at all (numpop==1) */

  which = numpop + RANDINT (0L, numpop2 - numpop - 1);
  while (shortcut (which, world, &w))
    {
      which = numpop + RANDINT (0L, numpop2 - numpop - 1);
    }

  Mold = world->param0[which];
  logratio = delta * (2.0 * UNIF_RANDUM () - 1.0);   /* Uniform(-delta,delta) */
  Mnew = Mold * EXP (logratio);

  if (Mnew <= bayes->minparam[which] || Mnew >= bayes->maxparam[which])
    return FALSE;

  *which_out = which;
  *Mold_out = Mold;
  *Mnew_out = Mnew;
  *logpriorjac_out = (*log_prior_ratio[which]) (Mnew, Mold, bayes, which)
                      + LOG (Mnew) - LOG (Mold);
  return TRUE;
}

/*========================================================================
  Dense-matrix helpers and the CTMC uniformization bridge sampler.

  Nothing in migrate builds/exponentiates a dense generator matrix or
  samples an endpoint-conditioned CTMC path today (checked in the C
  primitives survey before starting this file) -- migration histories are
  instead grown one event at a time, competing against coalescent hazard,
  by tree_update()'s eventtime_single()/time_to_coalmig(). This move needs
  something different: a branch's population at BOTH ends is already fixed
  (by neighbouring window/non-window nodes), and the task is to redraw a
  plausible migration path between two fixed endpoints under a candidate M
  -- an endpoint-conditioned CTMC bridge, not an open-ended competing-hazard
  draw. Ported from the validated Python prototype's `_bridge()`
  (tpsc_windowed_bridge_move.py), generalized from 2 states to arbitrary
  numpop.
  ========================================================================*/

static void
dmat_eye (MYREAL * A, long d)
{
  long i, j;
  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      A[i * d + j] = (i == j) ? 1.0 : 0.0;
}

static void
dmat_mul (const MYREAL * A, const MYREAL * B, MYREAL * C, long d)
{
  long i, j, k;
  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      {
        MYREAL s = 0.0;
        for (k = 0; k < d; k++)
          s += A[i * d + k] * B[k * d + j];
        C[i * d + j] = s;
      }
}

/* Uniformization context for one branch's bridge draw: R = I + Q/mu (mu >=
   max_i(-Q[i][i]), so R's entries are all >= 0 and each row sums to 1 -- a
   genuine DTMC transition matrix), and its powers R^0..R^nmax, nmax chosen
   generously past the mean of the implied Poisson(mu*T) event count (same
   truncation as the Python prototype: mean + 10*sd + 20). Computed once per
   branch and shared between the sampler and the bridge-density evaluator
   (both need Rpow), rather than recomputed per query. */
typedef struct
{
  long d;
  long nmax;
  MYREAL *Rpow;                 /* (nmax+1) x d x d, row-major */
} unif_fmt;

#define RPOW(u,n,i,j) ((u)->Rpow[(n)*(u)->d*(u)->d + (i)*(u)->d + (j)])

static void
unif_build (const MYREAL * Q, long d, MYREAL T, MYREAL mu, unif_fmt * u)
{
  MYREAL lam = mu * T;
  MYREAL *R;
  long n, i, j;

  u->d = d;
  u->nmax = (long) (lam + 10.0 * sqrt ((double) lam) + 20.0);
  u->Rpow = (MYREAL *) mycalloc ((u->nmax + 1) * d * d, sizeof (MYREAL));
  R = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      R[i * d + j] = ((i == j) ? 1.0 : 0.0) + Q[i * d + j] / mu;
  dmat_eye (u->Rpow, d);        /* Rpow[0] = I */
  for (n = 1; n <= u->nmax; n++)
    dmat_mul (&u->Rpow[(n - 1) * d * d], R, &u->Rpow[n * d * d], d);
  myfree (R);
}

static void
unif_free (unif_fmt * u)
{
  myfree (u->Rpow);
}

/* mu must satisfy mu >= max_i(-Q[i][i]); the caller-visible entry points
   below derive it the same way the Python prototype did (max diagonal *
   1.05 + a small epsilon, to stay strictly clear of the boundary R>=0
   needs). Factored out so the sampler and the transition-matrix validator
   below always agree on mu for the same Q. */
static MYREAL
uniformization_rate (const MYREAL * Q, long d)
{
  long i;
  MYREAL mu = 0.0;
  for (i = 0; i < d; i++)
    if (-Q[i * d + i] > mu)
      mu = -Q[i * d + i];
  return mu * 1.05 + 1e-9;
}

/* P(T) = exp(-lam) * sum_n lam^n/n! * R^n, lam = mu*T -- the same
   uniformization series used by the sampler, evaluated in full rather than
   conditioned on an endpoint. This is used ONLY for validating the sampler
   below against an independent closed form (2-state case); the move itself
   never needs a full transition matrix, only bridge draws. */
static void
ctmc_transition_matrix (const unif_fmt * u, MYREAL lam, MYREAL * P)
{
  long n, i, j;
  MYREAL logfact = 0.0;
  long d = u->d;

  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      P[i * d + j] = 0.0;
  for (n = 0; n <= u->nmax; n++)
    {
      MYREAL logpois, w;
      if (n == 0)
        logpois = -lam;
      else
        {
          logfact += LOG ((MYREAL) n);
          logpois = -lam + n * LOG (lam) - logfact;
        }
      w = EXP (logpois);
      for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
          P[i * d + j] += w * RPOW (u, n, i, j);
    }
}

/* Exact CTMC sample-path log-density for one branch: diagonal survival
   between events plus the off-diagonal jump-rate factor at each event.
   Unconditional on the endpoint, exactly like the Python prototype's
   log_path_density -- already fully specifies the realized trajectory, no
   separate endpoint-marginal division/multiplication needed. `path` times
   are ABSOLUTE; t0 is the branch's own start time.

   migr_table_fmt field convention (confirmed against migrate()/
   insert_migr_node(), mcmc1.c/tools.c -- opposite of the naive
   forward-time reading): for ANY event, `.to` is the population on the
   YOUNGER (tip-ward) side and `.from` is the population on the OLDER
   (root-ward) side, regardless of array position -- migrate()'s own
   construction confirms this exactly (`array[0].to == up->pop`, i.e. the
   FIRST event's `.to` equals the branch's young end). So walking forward
   in time, the running state `s` must equal path[i].to just before
   processing event i, and transitions TO path[i].from afterward. */
static MYREAL
log_path_density (const MYREAL * Q, long d, MYREAL T, long start,
                   const migr_table_fmt * path, long path_n, MYREAL t0)
{
  long s = start;
  MYREAL t_prev = 0.0;
  MYREAL logp = 0.0;
  long i;

  for (i = 0; i < path_n; i++)
    {
      MYREAL t = path[i].time - t0;
      MYREAL dt = t - t_prev;
      logp += Q[s * d + s] * dt;
      logp += LOG (Q[s * d + path[i].from]);
      s = path[i].from;
      t_prev = t;
    }
  logp += Q[s * d + s] * (T - t_prev);
  return logp;
}

/* Endpoint-conditional bridge density of a specific realized path:
   log_path_density minus the log endpoint-transition probability. Needed
   for the move's Hastings ratio (task not yet reached -- see
   windowed_joint_update()'s TODO) and used below to cross-validate the
   sampler. */
static MYREAL
log_bridge_density (const MYREAL * Q, long d, MYREAL T, long start, long end,
                     const migr_table_fmt * path, long path_n, MYREAL t0)
{
  MYREAL mu = uniformization_rate (Q, d);
  unif_fmt u;
  MYREAL *P;
  MYREAL logZ;
  MYREAL logpd;

  unif_build (Q, d, T, mu, &u);
  P = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  ctmc_transition_matrix (&u, mu * T, P);
  logZ = LOG (P[start * d + end] > 1e-300 ? P[start * d + end] : 1e-300);
  myfree (P);
  unif_free (&u);

  logpd = log_path_density (Q, d, T, start, path, path_n, t0);
  return logpd - logZ;
}

/* Samples a CTMC bridge on (t0, t0+T] from state `a` to state `b` under
   generator Q (d x d, row-major): draws the number of "virtual"
   uniformization events, then the intermediate virtual states via the
   standard forward/backward chain-rule bridge, then assigns each virtual
   event an i.i.d. Uniform(t0,t0+T) time (a property of the underlying
   Poisson process, independent of the state chain), and finally collapses
   consecutive equal states (virtual "self-jumps", an artifact of
   uniformization with no CTMC-event meaning) into the real migration
   events. Direct port of the Python prototype's `_bridge()`, generalized
   from 2 states to arbitrary d.

   Writes the result into *out_n and returns a mycalloc'd migr_table_fmt
   array of that length (event='m' throughout -- divergence ('d') events
   are not part of this move); *out_n==0 means a real zero-event bridge (a
   valid, common outcome, e.g. a>b's uniformization simply never left a on
   this draw) and the returned pointer is NULL in that case, matching how
   callers already treat "no events on this branch" elsewhere in this file.
   Returns NULL with *out_n==-1 only on the pathological failure of the
   forward chain-rule step finding zero total probability at every
   candidate state (should not happen for an irreducible Q with mu chosen
   as above; guarded rather than assumed). */
static migr_table_fmt *
ctmc_bridge (const MYREAL * Q, long d, MYREAL T, long a, long b, MYREAL t0,
             long *out_n)
{
  MYREAL mu = uniformization_rate (Q, d);
  unif_fmt u;
  MYREAL lam;
  MYREAL *w;
  MYREAL logfact = 0.0;
  MYREAL maxlw = -HUGE;
  MYREAL *logw;
  MYREAL tot;
  MYREAL r, cum;
  long n, i, j, chosen_n = -1;
  long *states;
  MYREAL *ts;
  migr_table_fmt *out;
  long nout;
  long prev;

  unif_build (Q, d, T, mu, &u);
  lam = mu * T;

  logw = (MYREAL *) mycalloc (u.nmax + 1, sizeof (MYREAL));
  for (n = 0; n <= u.nmax; n++)
    {
      MYREAL rab = RPOW (&u, n, a, b);
      MYREAL logpois;
      if (n == 0)
        logpois = -lam;
      else
        {
          logfact += LOG ((MYREAL) n);
          logpois = -lam + n * LOG (lam) - logfact;
        }
      logw[n] = logpois + LOG (rab > 1e-300 ? rab : 1e-300);
      if (logw[n] > maxlw)
        maxlw = logw[n];
    }
  w = (MYREAL *) mycalloc (u.nmax + 1, sizeof (MYREAL));
  tot = 0.0;
  for (n = 0; n <= u.nmax; n++)
    {
      w[n] = EXP (logw[n] - maxlw);
      tot += w[n];
    }
  myfree (logw);

  r = UNIF_RANDUM () * tot;
  cum = 0.0;
  for (n = 0; n <= u.nmax; n++)
    {
      cum += w[n];
      if (cum >= r)
        {
          chosen_n = n;
          break;
        }
    }
  if (chosen_n < 0)
    chosen_n = u.nmax;          /* rounding fallback, matches the Python np.random.choice edge */
  myfree (w);

  if (chosen_n == 0)
    {
      unif_free (&u);
      *out_n = 0;
      return NULL;
    }

  states = (long *) mycalloc (chosen_n + 1, sizeof (long));
  states[0] = a;
  for (i = 0; i < chosen_n - 1; i++)
    {
      MYREAL *pr = (MYREAL *) mycalloc (d, sizeof (MYREAL));
      MYREAL prtot = 0.0;
      MYREAL rr, c;
      long s = states[i];
      long chosen_s = -1;

      for (j = 0; j < d; j++)
        {
          pr[j] = RPOW (&u, 1, s, j) * RPOW (&u, chosen_n - i - 1, j, b);
          prtot += pr[j];
        }
      if (prtot <= 0.0)
        {
          myfree (pr);
          myfree (states);
          unif_free (&u);
          *out_n = -1;
          return NULL;
        }
      rr = UNIF_RANDUM () * prtot;
      c = 0.0;
      for (j = 0; j < d; j++)
        {
          c += pr[j];
          if (c >= rr)
            {
              chosen_s = j;
              break;
            }
        }
      if (chosen_s < 0)
        chosen_s = d - 1;
      myfree (pr);
      states[i + 1] = chosen_s;
    }
  states[chosen_n] = b;
  unif_free (&u);

  ts = (MYREAL *) mycalloc (chosen_n, sizeof (MYREAL));
  for (i = 0; i < chosen_n; i++)
    ts[i] = (MYREAL) (UNIF_RANDUM () * T);
  for (i = 1; i < chosen_n; i++)      /* insertion sort: chosen_n is small */
    {
      MYREAL key = ts[i];
      long kk = i - 1;
      while (kk >= 0 && ts[kk] > key)
        {
          ts[kk + 1] = ts[kk];
          kk--;
        }
      ts[kk + 1] = key;
    }

  out = (migr_table_fmt *) mycalloc (chosen_n, sizeof (migr_table_fmt));
  nout = 0;
  prev = a;
  for (i = 1; i <= chosen_n; i++)
    {
      if (states[i] != prev)
        {
          /* migr_table_fmt convention (see log_path_density's comment):
             .to = younger/tip-ward side = prev (state before this local
             transition), .from = older/root-ward side = states[i] (state
             after) -- opposite of the naive forward-time reading. */
          out[nout].from = states[i];
          out[nout].to = prev;
          out[nout].time = t0 + ts[i - 1];
          out[nout].event = 'm';
          nout++;
          prev = states[i];
        }
    }
  myfree (ts);
  myfree (states);
  *out_n = nout;
  if (nout == 0)
    {
      myfree (out);
      return NULL;
    }
  return out;
}

/* Self-test for the matrix/bridge machinery above, run once (guarded by a
   static flag) when MIGRATE_WINDOW_DEBUG is set, mirroring the validation
   discipline used throughout the Python prototyping (verify before
   trusting -- project memory "feedback-verify-before-proposing-fixes").
   Three checks, in increasing order of what they cover:
     1. ctmc_transition_matrix (the uniformization series) against an
        INDEPENDENT closed form for a symmetric 2-state CTMC:
        P(T)[a,a] = 0.5*(1+exp(-2*M*T)). Validates the Poisson-weighted
        Rpow machinery on its own.
     2. log_bridge_density for the trivial empty-path case against the same
        closed form algebraically simplified -- log P_path(no event) -
        log P(T)[a,a] = -M*T - log(0.5*(1+exp(-2*M*T))). Validates
        log_path_density/log_bridge_density independently of the sampler.
     3. Monte Carlo: draw many bridge samples for the SAME (a=0,b=0,T) and
        check the empirical fraction that come back with zero REAL events
        (after collapsing virtual self-jumps) against exp(log_bridge_density
        for the empty path) from check 2 -- ties the sampler itself to the
        already-validated density formula, the strongest check available
        without a full brute-force enumeration (impractical here: even a
        2-state bridge has unboundedly many possible event counts). */
/*========================================================================
  Mean-field ("bystander-extended") forward pass and windowed FFBS.

  Ported from the Python prototype's forward_bystander()/
  sample_or_score_labels() (tpsc_windowed_bridge_move.py), generalized from
  2 populations to arbitrary numpop. Operates purely on the genealogy's
  coalescent skeleton (topology + node times, via crawlforward() to skip
  any existing migration nodes) -- never on the migration history currently
  sitting on any branch, since the point of this pass is to inform a
  REPLACEMENT for some of that history, not condition on it.
  ========================================================================*/

/* row-vector (length d) times d x d matrix -> row-vector; the operation
   dmat_mul() (square matrix x square matrix) doesn't cover, needed for
   propagating a population-probability message through one interval. */
static void
vec_mat_mul (const MYREAL * v, const MYREAL * M, MYREAL * out, long d)
{
  long i, j;
  for (j = 0; j < d; j++)
    {
      MYREAL s = 0.0;
      for (i = 0; i < d; i++)
        s += v[i] * M[i * d + j];
      out[j] = s;
    }
}

/* Pointer -> nodep[] index, by linear scan. world->nodep[] holds every tip
   (indices [0,sumtips)) and internal coalescent node (indices
   [sumtips,2*sumtips-2]) as allocated -- every node pointer this file ever
   sees (tips, crawlforward() results, window/touched-branch entries) is one
   of exactly these pointers, so this always finds a match. Not indexed by
   node->id: confirmed by the earlier C-primitives survey that ->id is
   unset/vestigial for these node types outside #ifdef MESS debug builds,
   so it cannot be trusted as a pre-existing index. Tree sizes here are
   small (tens to a few hundred tips), so an O(sumtips) scan, called O(1)
   times per processed node, is not worth optimizing away. */
static long
node_index (world_fmt * world, node * p)
{
  long i;
  long hi = 2 * world->sumtips - 1;
  for (i = 0; i < hi; i++)
    if (world->nodep[i] == p)
      return i;
  error ("node_index: node not found in world->nodep[]");
  return -1;
}

static boolean
in_window (node * v, node ** window, long window_n)
{
  long i;
  for (i = 0; i < window_n; i++)
    if (window[i] == v)
      return TRUE;
  return FALSE;
}

static long
categorical_sample (const MYREAL * p, long d)
{
  MYREAL r = UNIF_RANDUM ();
  MYREAL cum = 0.0;
  long k;
  for (k = 0; k < d; k++)
    {
      cum += p[k];
      if (cum >= r)
        return k;
    }
  return d - 1;
}

/* internal coalescent nodes (world->nodep[sumtips..2*sumtips-2]), sorted by
   increasing tyme, returned as an array of NODEP INDICES (not pointers) so
   a node's identity survives the sort with no further pointer->index
   lookup needed -- order_idx[ninternal-1] is always the MRCA (the internal
   node with the largest tyme), which is this pass's "root" for FFBS
   purposes (distinct from world->root, migrate's separate 'r' sentinel
   ABOVE the MRCA, which is never part of this range -- see the C
   primitives survey). Caller must myfree() the result. */
static long *
sorted_internal_indices (world_fmt * world)
{
  const long sumtips = world->sumtips;
  const long ninternal = sumtips - 1;
  long *order_idx = (long *) mycalloc (ninternal, sizeof (long));
  long i;
  for (i = 0; i < ninternal; i++)
    order_idx[i] = sumtips + i;
  for (i = 1; i < ninternal; i++)
    {
      long key = order_idx[i];
      MYREAL keytyme = world->nodep[key]->tyme;
      long kk = i - 1;
      while (kk >= 0 && world->nodep[order_idx[kk]]->tyme > keytyme)
        {
          order_idx[kk + 1] = order_idx[kk];
          kk--;
        }
      order_idx[kk + 1] = key;
    }
  return order_idx;
}

typedef struct
{
  long d;
  long nnodes;                  /* 2*sumtips-1: valid range of world->nodep[] */
  MYREAL *chk0;                 /* [nnodes*d]: birth message (tip: one-hot at
                                    tip->pop; internal: the merged, renormalized
                                    distribution at coalescence -- doubles as
                                    the Python prototype's separate
                                    merge_part[v], which is identical to
                                    chk0[v] for every internal v there too) */
  MYREAL *Tacc;                 /* [nnodes*d*d]: accumulated effective
                                    (migration + bystander-reweight) operator
                                    from birth to the lineage's own death;
                                    identity if it never survived an interval */
} forward_ctx;

static void
forward_ctx_free (forward_ctx * ctx)
{
  myfree (ctx->chk0);
  myfree (ctx->Tacc);
}

/* Mean-field forward pass: propagates each lineage's population-membership
   probability continuously via Q between coalescent events, reweighting for
   competition with EVERY other currently-active lineage's coalescence
   hazard (not just the pair that eventually merges -- the "bystander"
   extension over a naive TPSC-style pruning pass, validated in the Python
   prototype to matter a great deal for acceptance at realistic lineage
   counts). Because the reweighting a lineage experiences between two
   coalescent events depends only on OTHER lineages' already-computed
   messages, not on its own state, its whole per-interval operator
   (migrate, then diagonally reweight) is a linear map independent of the
   lineage's own distribution -- so it composes via ordinary matrix
   multiplication into a single accumulated Tacc, with no chained
   multi-checkpoint forward-filtering needed (this simplification, and why
   an earlier more complicated design was abandoned, is recorded in the
   Python prototype's own file header/project memory). */
static void
forward_bystander (world_fmt * world, const MYREAL * Q, forward_ctx * ctx)
{
  const long d = world->numpop;
  const long sumtips = world->sumtips;
  const long ninternal = sumtips - 1;
  const long nnodes = 2 * sumtips - 1;
  MYREAL *invtheta;
  long *order_idx;
  node **active_v;
  long *active_idx;
  MYREAL *active_msg;
  long nactive;
  MYREAL t_prev;
  MYREAL *Pdt, *step_op, *migrated, *c_vec, *pv, *tmp_mat, *raw_vec;
  long i, j, k, t;

  ctx->d = d;
  ctx->nnodes = nnodes;
  ctx->chk0 = (MYREAL *) mycalloc (nnodes * d, sizeof (MYREAL));
  ctx->Tacc = (MYREAL *) mycalloc (nnodes * d * d, sizeof (MYREAL));

  invtheta = (MYREAL *) mycalloc (d, sizeof (MYREAL));
  for (i = 0; i < d; i++)
    invtheta[i] = 1.0 / world->param0[i];

  for (i = 0; i < sumtips; i++)
    {
      node *tip = world->nodep[i];
      for (j = 0; j < d; j++)
        ctx->chk0[i * d + j] = (j == tip->pop) ? 1.0 : 0.0;
      dmat_eye (&ctx->Tacc[i * d * d], d);
    }

  order_idx = sorted_internal_indices (world);

  active_v = (node **) mycalloc (sumtips, sizeof (node *));
  active_idx = (long *) mycalloc (sumtips, sizeof (long));
  active_msg = (MYREAL *) mycalloc (sumtips * d, sizeof (MYREAL));
  for (i = 0; i < sumtips; i++)
    {
      active_v[i] = world->nodep[i];
      active_idx[i] = i;
      for (j = 0; j < d; j++)
        active_msg[i * d + j] = ctx->chk0[i * d + j];
    }
  nactive = sumtips;
  t_prev = 0.0;

  Pdt = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  step_op = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  migrated = (MYREAL *) mycalloc (sumtips * d, sizeof (MYREAL));
  c_vec = (MYREAL *) mycalloc (d, sizeof (MYREAL));
  pv = (MYREAL *) mycalloc (d, sizeof (MYREAL));
  tmp_mat = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  raw_vec = (MYREAL *) mycalloc (d, sizeof (MYREAL));

  for (t = 0; t < ninternal; t++)
    {
      long vidx = order_idx[t];
      node *v = world->nodep[vidx];
      node *ca = crawlforward (v->next);
      node *cb = crawlforward (v->next->next);
      MYREAL dt = v->tyme - t_prev;
      long slot_a = -1, slot_b = -1;
      MYREAL s;

      if (dt > 1e-14)
        {
          MYREAL mu = uniformization_rate (Q, d);
          unif_fmt u;
          unif_build (Q, d, dt, mu, &u);
          ctmc_transition_matrix (&u, mu * dt, Pdt);
          unif_free (&u);
        }
      else
        {
          dmat_eye (Pdt, d);
        }

      for (i = 0; i < nactive; i++)
        {
          if (active_v[i] == ca)
            slot_a = i;
          if (active_v[i] == cb)
            slot_b = i;
        }
      if (slot_a < 0 || slot_b < 0)
        error ("forward_bystander: merging child not found in active set");

      for (i = 0; i < nactive; i++)
        vec_mat_mul (&active_msg[i * d], Pdt, &migrated[i * d], d);

      for (i = 0; i < nactive; i++)
        {
          boolean i_is_pair = (i == slot_a || i == slot_b);
          for (k = 0; k < d; k++)
            c_vec[k] = 0.0;
          for (j = 0; j < nactive; j++)
            {
              if (j == i)
                continue;
              if (i_is_pair && (j == slot_a || j == slot_b))
                continue;
              for (k = 0; k < d; k++)
                c_vec[k] += migrated[j * d + k];
            }
          for (k = 0; k < d; k++)
            c_vec[k] *= invtheta[k];

          for (j = 0; j < d; j++)
            for (k = 0; k < d; k++)
              step_op[j * d + k] = Pdt[j * d + k] * EXP (-dt * c_vec[k]);

          {
            MYREAL *cur = &ctx->Tacc[active_idx[i] * d * d];
            dmat_mul (cur, step_op, tmp_mat, d);
            memcpy (cur, tmp_mat, (size_t) (d * d) * sizeof (MYREAL));
          }

          vec_mat_mul (&active_msg[i * d], step_op, raw_vec, d);
          s = 0.0;
          for (k = 0; k < d; k++)
            s += raw_vec[k];
          if (s > 0.0)
            for (k = 0; k < d; k++)
              active_msg[i * d + k] = raw_vec[k] / s;
          else
            for (k = 0; k < d; k++)
              active_msg[i * d + k] = migrated[i * d + k];
        }

      for (k = 0; k < d; k++)
        pv[k] = active_msg[slot_a * d + k] * active_msg[slot_b * d + k]
                * invtheta[k];
      s = 0.0;
      for (k = 0; k < d; k++)
        s += pv[k];
      for (k = 0; k < d; k++)
        ctx->chk0[vidx * d + k] = pv[k] / s;
      dmat_eye (&ctx->Tacc[vidx * d * d], d);

      active_v[slot_a] = v;
      active_idx[slot_a] = vidx;
      for (k = 0; k < d; k++)
        active_msg[slot_a * d + k] = ctx->chk0[vidx * d + k];
      if (slot_b != nactive - 1)
        {
          active_v[slot_b] = active_v[nactive - 1];
          active_idx[slot_b] = active_idx[nactive - 1];
          for (k = 0; k < d; k++)
            active_msg[slot_b * d + k] = active_msg[(nactive - 1) * d + k];
        }
      nactive -= 1;

      t_prev = v->tyme;
    }

  myfree (invtheta);
  myfree (order_idx);
  myfree (active_v);
  myfree (active_idx);
  myfree (active_msg);
  myfree (Pdt);
  myfree (step_op);
  myfree (migrated);
  myfree (c_vec);
  myfree (pv);
  myfree (tmp_mat);
  myfree (raw_vec);
}

/* Windowed forward-filtering/backward-sampling over internal-node
   populations, given forward_bystander()'s chk0/Tacc. Only nodes in
   `window` get a fresh draw (sample=TRUE) or are scored against their OWN
   current population (sample=FALSE, for the reverse-direction Hastings
   term); every other internal node keeps its CURRENT population unchanged
   -- read directly from the tree, not resampled, and contributing no logq
   term at all, since it is not part of what this move proposes. Tips are
   never touched. Traverses internal nodes in DECREASING time order (root/
   MRCA first), matching the Python prototype's sample_or_score_labels():
   a node's own chosen population must already be known (by the time its
   children are processed) from having been set either at the root's own
   handling or by ITS parent's processing earlier in this same loop, since
   a parent always has a larger tyme than its children.

   Writes every internal node's chosen population into new_pop[] (indexed
   by nodep index; entries for tips and for indices outside
   [sumtips,2*sumtips-2] are left untouched). Returns the summed log
   density of whatever got FRESHLY evaluated at window nodes (0 if the
   window is empty), or -HUGE if a zero-probability dead end is hit -- the
   caller must check for -HUGE before trusting new_pop, exactly like the
   Python prototype's `return None, -np.inf`. */
static MYREAL
sample_or_score_labels (world_fmt * world, const forward_ctx * ctx,
                         node ** window, long window_n, boolean sample,
                         long *new_pop)
{
  const long d = ctx->d;
  const long sumtips = world->sumtips;
  const long ninternal = sumtips - 1;
  long *order_idx = sorted_internal_indices (world);
  MYREAL logq = 0.0;
  MYREAL *w = (MYREAL *) mycalloc (d, sizeof (MYREAL));
  long t;

  for (t = ninternal - 1; t >= 0; t--)
    {
      long vidx = order_idx[t];
      node *v = world->nodep[vidx];
      node *kids[2];
      long kk;

      if (t == ninternal - 1)         /* MRCA: this pass's "root", handled directly */
        {
          if (in_window (v, window, window_n))
            {
              MYREAL s = 0.0;
              long x, state;
              for (x = 0; x < d; x++)
                {
                  w[x] = ctx->chk0[vidx * d + x];
                  s += w[x];
                }
              for (x = 0; x < d; x++)
                w[x] /= s;
              state = sample ? categorical_sample (w, d) : v->pop;
              logq += LOG (w[state] > 1e-300 ? w[state] : 1e-300);
              new_pop[vidx] = state;
            }
          else
            {
              new_pop[vidx] = v->pop;
            }
        }
      /* else: new_pop[vidx] was already set when v's own parent was
         processed earlier in this loop (parents precede children here) */

      kids[0] = crawlforward (v->next);
      kids[1] = crawlforward (v->next->next);
      for (kk = 0; kk < 2; kk++)
        {
          node *c = kids[kk];
          long cidx;
          if (c->type != 'i')
            continue;            /* tip: observed, never resampled */
          cidx = node_index (world, c);
          if (in_window (c, window, window_n))
            {
              MYREAL tot = 0.0;
              long x, state;
              for (x = 0; x < d; x++)
                {
                  w[x] = ctx->chk0[cidx * d + x]
                         * ctx->Tacc[cidx * d * d + x * d + new_pop[vidx]];
                  tot += w[x];
                }
              if (tot <= 0.0)
                {
                  myfree (w);
                  myfree (order_idx);
                  return -HUGE;
                }
              for (x = 0; x < d; x++)
                w[x] /= tot;
              state = sample ? categorical_sample (w, d) : c->pop;
              logq += LOG (w[state] > 1e-300 ? w[state] : 1e-300);
              new_pop[cidx] = state;
            }
          else
            {
              new_pop[cidx] = c->pop;
            }
        }
    }

  myfree (w);
  myfree (order_idx);
  return logq;
}

/* Brute-force normalization check on a REAL, currently-running tree's
   window: enumerate every one of d^window_n possible population
   assignments for the window nodes, score each via sample_or_score_labels
   (sample=FALSE, temporarily forcing the tree to that assignment), and sum
   exp(logq) -- must equal 1.0, since sample_or_score_labels's sample=TRUE
   mode is supposed to draw exactly this conditional distribution over the
   window given everything else fixed. This is the live-tree counterpart of
   the Python prototype's check_windowed_normalization() (which used a
   hand-built 3-tip synthetic tree); using the actual running genealogy
   instead avoids needing to hand-construct a synthetic world_fmt/node tree
   in C purely for testing. Restores the tree's original populations before
   returning either way. Silently skipped (no output) if d^window_n would
   be impractically large for a live per-move check. */
static void
window_ffbs_selftest (world_fmt * world, const forward_ctx * ctx,
                       node ** window, long window_n)
{
  const long d = ctx->d;
  long ncombo = 1;
  long *saved_pop;
  long *saved_actualpop;
  long *new_pop;
  long i, combo;
  MYREAL total = 0.0;
  boolean ok;

  for (i = 0; i < window_n; i++)
    {
      ncombo *= d;
      if (ncombo > 5000)
        return;                  /* too large for a live check; skip silently */
    }
  if (window_n == 0)
    return;

  saved_pop = (long *) mycalloc (window_n, sizeof (long));
  saved_actualpop = (long *) mycalloc (window_n, sizeof (long));
  for (i = 0; i < window_n; i++)
    {
      saved_pop[i] = window[i]->pop;
      saved_actualpop[i] = window[i]->actualpop;
    }

  new_pop = (long *) mycalloc (ctx->nnodes, sizeof (long));
  for (combo = 0; combo < ncombo; combo++)
    {
      long rem = combo;
      MYREAL logq;
      for (i = 0; i < window_n; i++)
        {
          long digit = rem % d;
          rem /= d;
          window[i]->pop = digit;
          window[i]->actualpop = digit;
        }
      logq = sample_or_score_labels (world, ctx, window, window_n, FALSE,
                                      new_pop);
      if (logq > -HUGE / 2.0)
        total += EXP (logq);
    }
  myfree (new_pop);

  for (i = 0; i < window_n; i++)
    {
      window[i]->pop = saved_pop[i];
      window[i]->actualpop = saved_actualpop[i];
    }
  myfree (saved_pop);
  myfree (saved_actualpop);

  ok = fabs ((double) (total - 1.0)) < 1e-6;
  fprintf (stderr, "[window_ffbs_selftest] window_n=%li d=%li ncombo=%li "
                    "total=%.8f  %s\n",
            window_n, d, ncombo, (double) total, ok ? "OK" : "FAIL");
}

static void
windowupdate_selftest (void)
{
  const MYREAL M = 137.0;
  const MYREAL T = 0.02;
  const long d = 2;
  MYREAL Q[4];
  unif_fmt u;
  MYREAL mu;
  MYREAL P[4];
  MYREAL want_P00, got_P00;
  MYREAL want_lbd, got_lbd;
  long trial, ntrials = 20000, zero_count = 0;
  MYREAL empirical_zero_frac, want_zero_frac;
  boolean ok = TRUE;

  Q[0 * d + 0] = -M; Q[0 * d + 1] = M;
  Q[1 * d + 0] = M;  Q[1 * d + 1] = -M;

  mu = uniformization_rate (Q, d);
  unif_build (Q, d, T, mu, &u);
  ctmc_transition_matrix (&u, mu * T, P);
  unif_free (&u);
  got_P00 = P[0 * d + 0];
  want_P00 = 0.5 * (1.0 + EXP (-2.0 * M * T));
  fprintf (stderr, "[window_selftest] check 1 (transition matrix vs closed form): "
                    "got=%.8f want=%.8f  %s\n",
            (double) got_P00, (double) want_P00,
            (fabs ((double) (got_P00 - want_P00)) < 1e-8) ? "OK" : "FAIL");
  ok = ok && (fabs ((double) (got_P00 - want_P00)) < 1e-8);

  got_lbd = log_bridge_density (Q, d, T, 0, 0, NULL, 0, 0.0);
  want_lbd = (-M * T) - LOG (want_P00);
  fprintf (stderr, "[window_selftest] check 2 (bridge density, 0-event path): "
                    "got=%.8f want=%.8f  %s\n",
            (double) got_lbd, (double) want_lbd,
            (fabs ((double) (got_lbd - want_lbd)) < 1e-6) ? "OK" : "FAIL");
  ok = ok && (fabs ((double) (got_lbd - want_lbd)) < 1e-6);

  for (trial = 0; trial < ntrials; trial++)
    {
      long nout;
      migr_table_fmt *path = ctmc_bridge (Q, d, T, 0, 0, 0.0, &nout);
      if (nout == 0)
        zero_count++;
      if (path != NULL)
        myfree (path);
    }
  empirical_zero_frac = (MYREAL) zero_count / (MYREAL) ntrials;
  want_zero_frac = EXP (got_lbd);
  fprintf (stderr, "[window_selftest] check 3 (sampler vs density, Monte Carlo n=%li): "
                    "empirical=%.4f want=%.4f  %s\n",
            ntrials, (double) empirical_zero_frac, (double) want_zero_frac,
            (fabs ((double) (empirical_zero_frac - want_zero_frac)) < 0.02) ? "OK" : "FAIL");
  ok = ok && (fabs ((double) (empirical_zero_frac - want_zero_frac)) < 0.02);

  /* check 4: a genuine >0-event path against an independent closed form --
     none of checks 1-3 exercise the migr_table_fmt.from/.to field
     convention at all (all use the empty path), which is exactly how a
     real from/to swap bug slipped past them earlier in this file's
     history and was only caught by a live NODATA run hitting
     insert_migr_node()'s own internal consistency check. A single event
     0->1 at any time t1 in (0,T) has log_path_density = log(M) - M*T
     (the two diagonal-survival terms sum to -M*T regardless of where t1
     falls, by direct calculation) -- deliberately checked at a
     non-symmetric t1 (0.3*T) so an accidental T/2-only correctness
     coincidence can't hide a bug. */
  {
    migr_table_fmt path1[1];
    MYREAL want_lpd1, got_lpd1;
    path1[0].from = 1;      /* older/root-ward side = 1 (b) */
    path1[0].to = 0;        /* younger/tip-ward side = 0 (a) */
    path1[0].time = 0.3 * T;
    path1[0].event = 'm';
    got_lpd1 = log_path_density (Q, d, T, 0, path1, 1, 0.0);
    want_lpd1 = LOG (M) - M * T;
    fprintf (stderr, "[window_selftest] check 4 (single-event path density vs closed form): "
                      "got=%.8f want=%.8f  %s\n",
              (double) got_lpd1, (double) want_lpd1,
              (fabs ((double) (got_lpd1 - want_lpd1)) < 1e-6) ? "OK" : "FAIL");
    ok = ok && (fabs ((double) (got_lpd1 - want_lpd1)) < 1e-6);
  }

  /* check 5: parity invariant. In a 2-state chain, reaching b==a requires
     an EVEN number of real transitions and b!=a requires an ODD number,
     for EVERY sample, not just on average -- a hard structural check that
     the sampler's endpoint bookkeeping (and the from/to field convention
     used to collapse virtual self-jumps into real events) is self
     -consistent, independent of the density-formula cross-check above. */
  {
    long trial2, bad_parity = 0;
    long ntrials2 = 2000;
    for (trial2 = 0; trial2 < ntrials2; trial2++)
      {
        long nout;
        migr_table_fmt *path = ctmc_bridge (Q, d, T, 0, 1, 0.0, &nout);
        if ((nout % 2) == 0)
          bad_parity++;
        if (path != NULL)
          myfree (path);
      }
    for (trial2 = 0; trial2 < ntrials2; trial2++)
      {
        long nout;
        migr_table_fmt *path = ctmc_bridge (Q, d, T, 0, 0, 0.0, &nout);
        if ((nout % 2) != 0)
          bad_parity++;
        if (path != NULL)
          myfree (path);
      }
    fprintf (stderr, "[window_selftest] check 5 (endpoint-parity invariant, n=%li): "
                      "violations=%li  %s\n",
              2 * ntrials2, bad_parity, (bad_parity == 0) ? "OK" : "FAIL");
    ok = ok && (bad_parity == 0);
  }

  fprintf (stderr, "[window_selftest] overall: %s\n", ok ? "PASS" : "FAIL");
}

/*========================================================================
  The move itself: propose M' by a local step, jointly redraw migration
  history on a small window of branches, accept/reject via the exact
  Hastings ratio.
  ========================================================================*/

/* Existing migration/divergence events currently on branch `up` (walking
   toward the ancestor via showtop(->back), exactly like erase_migr_nodes()
   -- see that function's own comment for why `up->back`/`->next->back`
   already encode "toward parent"/"toward child"). Returns a mycalloc'd
   array of length *out_n (possibly 0, in which case the returned pointer
   is still valid but should not be dereferenced -- caller checks *out_n).
   Needed both to score the REVERSE-direction bridge density (what this
   move's Hastings ratio calls logq_bridge_old) and to restore a branch's
   history exactly if the move is rejected after already committing a
   trial redraw. */
static migr_table_fmt *
extract_migr_events (node * up, long *out_n)
{
  long n = 0;
  node *p = showtop (up->back);
  migr_table_fmt *tab;
  long i;

  while (p->type == 'm' || p->type == 'd')
    {
      n++;
      p = showtop (p->back);
    }
  tab = (migr_table_fmt *) mycalloc (n > 0 ? n : 1, sizeof (migr_table_fmt));
  p = showtop (up->back);
  i = 0;
  while (p->type == 'm' || p->type == 'd')
    {
      tab[i].from = p->pop;
      tab[i].to = p->actualpop;
      tab[i].time = p->tyme;
      tab[i].event = p->type;
      i++;
      p = showtop (p->back);
    }
  *out_n = n;
  return tab;
}

// Joint local-M + windowed migration-history move.
//
// Propose M' via local_M_propose() (a local, log-uniform step -- not an
// independence sampler, see that function's own comment for why); pick a
// window via window_pick()/window_touched_branches(); run the mean-field
// forward pass under M' (forward_bystander()) and draw the window nodes'
// new populations via windowed FFBS (sample_or_score_labels()); redraw
// each touched branch's migration history via the CTMC bridge sampler
// (ctmc_bridge()) conditioned on its (possibly freshly drawn) endpoint
// populations; accept/reject via the exact Hastings ratio, reusing
// migrate's own probg_treetimes()/calculate_prior()/bayes_accept() for the
// target-density and prior terms (no need to reimplement migrate's own
// coalescent/migration likelihood) and the bridge/label densities computed
// above for the proposal (q) terms. On rejection, every mutation made
// while trialling the move (parameter, window populations, touched
// branches' migration history) is explicitly undone, mirroring
// scaler_update()'s own reject-path pattern -- there is no journal/copy-on
// -write mechanism in migrate's tree representation to fall back on.
//
// One structural note not in the Python prototype this was ported from:
// a touched branch can be "the branch above a window node", and if that
// window node is the MRCA, its logical parent is world->root -- a pure
// bookkeeping sentinel placed 10000 time units above the MRCA (tree.c,
// e.g. ~line 1794), not a real branch with any likelihood content. Such
// branches are skipped for bridge redrawing (and excluded from every
// Hastings-ratio term below) -- there is nothing there to redraw.
long
windowed_joint_update (world_fmt * world)
{
  const long window_size = world->options->window_size;
  const long d = world->numpop;
  long which;
  MYREAL Mold, Mnew, logpriorjac;
  node **window, **touched;
  long window_n, touched_n;
  MYREAL *Qp, *Qm;
  forward_ctx ctxp, ctxo;
  long *new_pop, *dummy_pop;
  MYREAL logq_labels_new, logq_labels_old;
  MYREAL oldfull, oldlogprior;
  MYREAL hastings;
  boolean success;
  long i;

  /* per-touched-branch working arrays (only entries [0,nb) are real
     branches; branches skipped for the root-sentinel reason above are
     simply not appended, so nb <= touched_n) */
  node **b_up, **b_down;
  MYREAL *b_T;
  long *b_a_old, *b_b_old, *b_a_new, *b_b_new;
  migr_table_fmt **b_old_events, **b_new_events;
  long *b_old_n, *b_new_n;
  long nb;
  MYREAL logq_bridge_old = 0.0, logq_bridge_new = 0.0;
  boolean bridge_failed = FALSE;
  static boolean selftest_done = FALSE;
  boolean debug = (getenv ("MIGRATE_WINDOW_DEBUG") != NULL);

  if (debug && !selftest_done)
    {
      selftest_done = TRUE;
      windowupdate_selftest ();
    }

  world->window_trials += 1;

  if (!local_M_propose (world, &which, &Mold, &Mnew, &logpriorjac))
    return 0;

  window = (node **) mycalloc (window_size, sizeof (node *));
  touched = (node **) mycalloc (3 * window_size, sizeof (node *));
  window_n = window_pick (world, window, window_size);
  touched_n = window_touched_branches (window, window_n, touched,
                                        3 * window_size);

  if (debug)
    {
      long ii, jj;
      fprintf (stderr, "[window_debug] which=%li Mold=%f Mnew=%f window_n=%li touched_n=%li\n",
                which, (double) Mold, (double) Mnew, window_n, touched_n);
      for (ii = 0; ii < touched_n; ii++)
        for (jj = ii + 1; jj < touched_n; jj++)
          if (touched[ii] == touched[jj])
            fprintf (stderr, "  !! DUPLICATE touched entry at %li,%li\n", ii, jj);
      for (ii = 0; ii < window_n; ii++)
        for (jj = ii + 1; jj < window_n; jj++)
          if (window[ii] == window[jj])
            fprintf (stderr, "  !! DUPLICATE window entry at %li,%li\n", ii, jj);
      {
        MYREAL *Qdbg = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
        forward_ctx ctxdbg;
        build_Q (world, -1, 0.0, Qdbg);
        forward_bystander (world, Qdbg, &ctxdbg);
        window_ffbs_selftest (world, &ctxdbg, window, window_n);
        forward_ctx_free (&ctxdbg);
        myfree (Qdbg);
      }
    }

  /* forward pass + FFBS under M' (fresh draw for window nodes) */
  Qp = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  build_Q (world, which, Mnew, Qp);
  forward_bystander (world, Qp, &ctxp);
  new_pop = (long *) mycalloc (ctxp.nnodes, sizeof (long));
  logq_labels_new = sample_or_score_labels (world, &ctxp, window, window_n,
                                             TRUE, new_pop);
  if (logq_labels_new <= -HUGE / 2.0)
    {
      forward_ctx_free (&ctxp);
      myfree (Qp);
      myfree (new_pop);
      myfree (window);
      myfree (touched);
      return 0;
    }

  /* forward pass + FFBS under the CURRENT M (scoring the tree's own
     current window populations -- the reverse-direction logq term) */
  Qm = (MYREAL *) mycalloc (d * d, sizeof (MYREAL));
  build_Q (world, -1, 0.0, Qm);
  forward_bystander (world, Qm, &ctxo);
  dummy_pop = (long *) mycalloc (ctxo.nnodes, sizeof (long));
  logq_labels_old = sample_or_score_labels (world, &ctxo, window, window_n,
                                             FALSE, dummy_pop);
  myfree (dummy_pop);

  /* per-branch bookkeeping, and trial bridge draws -- all before any
     mutation, so a failure here (or above) costs nothing to undo */
  b_up = (node **) mycalloc (touched_n, sizeof (node *));
  b_down = (node **) mycalloc (touched_n, sizeof (node *));
  b_T = (MYREAL *) mycalloc (touched_n, sizeof (MYREAL));
  b_a_old = (long *) mycalloc (touched_n, sizeof (long));
  b_b_old = (long *) mycalloc (touched_n, sizeof (long));
  b_a_new = (long *) mycalloc (touched_n, sizeof (long));
  b_b_new = (long *) mycalloc (touched_n, sizeof (long));
  b_old_events = (migr_table_fmt **) mycalloc (touched_n, sizeof (migr_table_fmt *));
  b_new_events = (migr_table_fmt **) mycalloc (touched_n, sizeof (migr_table_fmt *));
  b_old_n = (long *) mycalloc (touched_n, sizeof (long));
  b_new_n = (long *) mycalloc (touched_n, sizeof (long));
  nb = 0;

  for (i = 0; i < touched_n; i++)
    {
      node *up = touched[i];
      node *down = showtop (crawlback (up));
      long aold, bold, anew, bnew;
      migr_table_fmt *oldtab;
      long oldn;

      if (down->type == 'r')
        continue;               /* fictitious branch above the MRCA: skip */

      aold = up->pop;
      bold = down->pop;
      anew = (up->type == 't') ? up->pop : new_pop[node_index (world, up)];
      bnew = (down->type == 't') ? down->pop : new_pop[node_index (world, down)];

      oldtab = extract_migr_events (up, &oldn);
      logq_bridge_old += log_bridge_density (Qm, d, down->tyme - up->tyme,
                                              aold, bold, oldtab, oldn,
                                              up->tyme);

      b_up[nb] = up;
      b_down[nb] = down;
      b_T[nb] = down->tyme - up->tyme;
      b_a_old[nb] = aold;
      b_b_old[nb] = bold;
      b_a_new[nb] = anew;
      b_b_new[nb] = bnew;
      b_old_events[nb] = oldtab;
      b_old_n[nb] = oldn;

      {
        long nout;
        migr_table_fmt *seg = ctmc_bridge (Qp, d, b_T[nb], anew, bnew,
                                            up->tyme, &nout);
        if (nout < 0)
          {
            bridge_failed = TRUE;
            b_new_events[nb] = NULL;
            b_new_n[nb] = 0;
            nb++;
            break;
          }
        b_new_events[nb] = seg;
        b_new_n[nb] = nout;
        logq_bridge_new += log_bridge_density (Qp, d, b_T[nb], anew, bnew,
                                                seg, nout, up->tyme);
      }
      nb++;
    }

  if (bridge_failed)
    {
      for (i = 0; i < nb; i++)
        {
          if (b_old_events[i] != NULL)
            myfree (b_old_events[i]);
          if (b_new_events[i] != NULL)
            myfree (b_new_events[i]);
        }
      myfree (b_up); myfree (b_down); myfree (b_T);
      myfree (b_a_old); myfree (b_b_old); myfree (b_a_new); myfree (b_b_new);
      myfree (b_old_events); myfree (b_new_events);
      myfree (b_old_n); myfree (b_new_n);
      forward_ctx_free (&ctxp); forward_ctx_free (&ctxo);
      myfree (Qp); myfree (Qm); myfree (new_pop);
      myfree (window); myfree (touched);
      return 0;
    }

  /* --- everything above is read-only; commit the trial proposal now --- */

  oldfull = probg_treetimes (world);
  oldlogprior = calculate_prior (world);

  world->param0[which] = Mnew;
  precalc_world (world);

  {
    long *win_old_pop = (long *) mycalloc (window_n, sizeof (long));
    MYREAL newfull, newlogprior, newval, oldval;

    for (i = 0; i < window_n; i++)
      {
        long idx = node_index (world, window[i]);
        win_old_pop[i] = window[i]->pop;
        /* set_pop() (tree.c) propagates pop/actualpop to EVERY nodelet in
           the ring (top, ->next, ->next->next for an 'i' node) -- directly
           assigning window[i]->pop/->actualpop only touches the top
           nodelet, leaving the two child-facing nodelets holding stale
           values that insert_migr_node()'s consistency check (via
           crawlback(up), which lands on exactly one of those child-facing
           nodelets, not the top one) then reads. Missing this was a real,
           corrupting bug caught by a live NODATA run -- existing migrate
           code (mcmc2.c's own insert_migr_node() callers) always calls
           set_pop() before insert_migr_node() for exactly this reason. */
        set_pop (window[i], new_pop[idx], new_pop[idx]);
      }

    for (i = 0; i < nb; i++)
      {
        erase_migr_nodes (world, b_up[i]);
        /* insert_migr_node() indexes migr_table[*counter-1] unconditionally
           (tools.c ~line 1964, no counter>0 guard on that particular
           access), so a genuine zero-event bridge -- a normal, common
           outcome, not a failure -- must not be passed through it at all.
           erase_migr_nodes() alone already leaves up->back == b_down[i]
           directly, which IS the correct zero-event state. */
        if (b_new_n[i] > 0)
          {
            if (getenv ("MIGRATE_WINDOW_DEBUG") != NULL)
              fprintf (stderr, "  [commit] branch %li: up.tyme=%f down.tyme=%f "
                                "down.type=%c down.actualpop=%li anew=%li bnew=%li "
                                "nevents=%li last.from=%li last.to=%li last.time=%f\n",
                        i, (double) b_up[i]->tyme, (double) b_down[i]->tyme,
                        b_down[i]->type, b_down[i]->actualpop,
                        b_a_new[i], b_b_new[i], b_new_n[i],
                        b_new_events[i][b_new_n[i]-1].from,
                        b_new_events[i][b_new_n[i]-1].to,
                        (double) b_new_events[i][b_new_n[i]-1].time);
            /* insert_migr_node()'s `down` argument must be the SPECIFIC
               child-facing nodelet (crawlback(up)'s raw, un-normalized
               result), not the canonical top nodelet b_down[i] holds for
               field reads (type/tyme/pop/actualpop are shared across a
               coalescent node's ring, but ->back is per-nodelet, one for
               each of the up-to-3 logical connections). Passing the
               showtop'd nodelet here overwrites the PARENT-facing link of
               b_down[i]'s logical node instead of the correct
               child-specific one -- a real, corrupting bug caught by a
               live NODATA run failing much later in add_partlineages(),
               well downstream of insert_migr_node()'s own (unaffected,
               since .actualpop is shared) consistency check. */
            {
              node *rawdown = crawlback (b_up[i]);
              node *topdown = showtop (rawdown);
              /* Defensive sync, matching the established idiom in
                 mcmc2.c's own insert_migr_node() callers (which call
                 set_pop() immediately before insert_migr_node() every
                 time): a coalescent node's non-top nodelets are not
                 necessarily kept in sync with the top one until something
                 actually needs them (confirmed live -- a node never
                 previously touched by any migration-event insertion had
                 its non-top nodelets still at allocate_nodelet()'s -1
                 default despite a perfectly valid top-nodelet pop). This
                 is idempotent: it only re-asserts the already-authoritative
                 top nodelet's own value onto its siblings, never changes
                 what population `down` actually represents. */
              set_pop (topdown, topdown->pop, topdown->actualpop);
              insert_migr_node (world, b_up[i], rawdown,
                                 b_new_events[i], &b_new_n[i]);
            }
          }
      }

    construct_tymelist (world, &world->treetimes[0]);

    newfull = probg_treetimes (world);
    newlogprior = calculate_prior (world);

    /* logpriorjac already carries the FULL prior-ratio + Jacobian for
       parameter `which` (via migrate's own log_prior_ratio[which], see
       local_M_propose()); newlogprior/oldlogprior are recomputed here only
       to keep world->logprior's bookkeeping current for other code that
       reads it (reporting, other moves) -- they are NOT added again into
       the Hastings ratio, which would double-count the same prior ratio. */
    hastings = (newfull - oldfull)
               + (logq_labels_old + logq_bridge_old)
               - (logq_labels_new + logq_bridge_new)
               + logpriorjac;

    if (world->options->prioralone)
      {
        newval = 0.0;
        oldval = 0.0;
      }
    else
      {
        /* unchanged either way: this move never touches branch lengths or
           topology, only migration history and node labels, so the data
           likelihood is identical before and after */
        newval = world->likelihood[world->G];
        oldval = world->likelihood[world->G];
      }

    success = bayes_accept (newval, oldval, world->heat, hastings);

    if (debug)
      fprintf (stderr, "  [decision] hastings=%f oldfull=%f newfull=%f -> %s\n",
                (double) hastings, (double) oldfull, (double) newfull,
                success ? "ACCEPT" : "REJECT");

    if (success)
      {
        world->window_accept += 1;
        world->logprior = newlogprior;
        world->bayes->oldval = newfull;
        world->param_like = world->bayes->oldval;
      }
    else
      {
        /* --- reject: undo every mutation, exactly --- */
        /* Window populations MUST be restored before the branches: a
           touched branch's `down` end can itself be a window node, and
           insert_migr_node() checks the LAST restored event's `.from`
           against down->actualpop's CURRENT value -- if that is still the
           NEW (not-yet-reverted) population when the old events are
           reinserted, the check fails even though nothing is actually
           wrong (this was a real bug caught by a live NODATA run: the
           commit-phase order, window pops then branches, is correct
           there, but the reject path needs the SAME order for the SAME
           reason, and originally had it backwards). */
        for (i = 0; i < window_n; i++)
          {
            set_pop (window[i], win_old_pop[i], win_old_pop[i]);
          }
        for (i = 0; i < nb; i++)
          {
            if (debug)
              fprintf (stderr, "  [rollback] branch %li: up.tyme=%f down.tyme=%f "
                                "down.type=%c down.actualpop=%li aold=%li bold=%li "
                                "nevents=%li%s\n",
                        i, (double) b_up[i]->tyme, (double) b_down[i]->tyme,
                        b_down[i]->type, b_down[i]->actualpop,
                        b_a_old[i], b_b_old[i], b_old_n[i],
                        b_old_n[i] > 0 ? "" : " (erase only)");
            erase_migr_nodes (world, b_up[i]);
            if (b_old_n[i] > 0)
              {
                node *rawdown = crawlback (b_up[i]);
                node *topdown = showtop (rawdown);
                set_pop (topdown, topdown->pop, topdown->actualpop);
                insert_migr_node (world, b_up[i], rawdown,
                                   b_old_events[i], &b_old_n[i]);
              }
            if (debug)
              {
                node *chk = showtop (crawlback (b_up[i]));
                long verify_n;
                migr_table_fmt *verify_tab = extract_migr_events (b_up[i], &verify_n);
                fprintf (stderr, "    [verify] branch %li: crawlback(up)==down? %s "
                                  "up.pop=%li(want %li) re-extracted n=%li(want %li)\n",
                          i, (chk == b_down[i]) ? "yes" : "NO",
                          b_up[i]->pop, b_a_old[i], verify_n, b_old_n[i]);
                myfree (verify_tab);
              }
          }
        world->param0[which] = Mold;
        precalc_world (world);
        construct_tymelist (world, &world->treetimes[0]);
        world->logprior = oldlogprior;
        world->bayes->oldval = probg_treetimes (world);
        world->param_like = world->bayes->oldval;
      }

    myfree (win_old_pop);
  }

  for (i = 0; i < nb; i++)
    {
      myfree (b_old_events[i]);
      myfree (b_new_events[i]);
    }
  myfree (b_up); myfree (b_down); myfree (b_T);
  myfree (b_a_old); myfree (b_b_old); myfree (b_a_new); myfree (b_b_new);
  myfree (b_old_events); myfree (b_new_events);
  myfree (b_old_n); myfree (b_new_n);
  forward_ctx_free (&ctxp);
  forward_ctx_free (&ctxo);
  myfree (Qp); myfree (Qm); myfree (new_pop);
  myfree (window); myfree (touched);

  return success ? 1 : 0;
}
