/*------------------------------------------------------
 Mittag-Leffler Hastings correction for genealogy moves

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
/*! \file mlf_hastings.c

 The genealogy move (newtree_update() in speciate.c) removes the branch above
 the origin and re-simulates it on the residual tree R with independent
 competing clocks (coalescence, migration) that restart at every event of R
 and at every event of the moving lineage. For exponential waiting times this
 is the conditional prior and the acceptance ratio is the data likelihood
 ratio alone. For Mittag-Leffler waiting times (alpha < 1) the proposal
 differs from the tree density used by probg_treetimes(), so the move is
 accepted with
     L(G')/L(G) * p(G')/p(G) * q(G | G') / q(G' | G).
 For alpha = 1 the p and q factors cancel exactly.

 Two details of the move matter for q(G | G'):
 - above the top of R only the root lineage remains. The moving lineage runs
   alone (the root lineage frozen) up to the horizon max(t_old, top of R), and
   jointly with the root lineage above it; the root lineage's migrations there
   become events of G'. The reverse move therefore works on
   R' = R - D + P2, where P2 are those new root-lineage events and D the
   events of R that the move removes (migration nodes above a new root);
 - the reverse move's horizon is max(t_new, top of R').
 See docs/mittag_leffler_in_migrate.tex.
 */
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "mittag_leffler.h"
#include "mlf_hastings.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

/* one event of a genealogy above the origin: the interval ending at the
   event has k[] lineages per population */
typedef struct
{
  double age;
  char type;      /* 'c' coalescence, 'm' migration, 't' other */
  long below;     /* population below the event */
  long above;     /* population above the event */
} mlh_event;

typedef struct
{
  long n;
  long alloc;
  long numpop;
  long rootpop;   /* population of the root lineage above the last event */
  mlh_event *ev;
  long *k;
} mlh_list;

/* a lineage path: start, migrations (time, from = population above,
   to = population below), end */
typedef struct
{
  double start;
  long startpop;
  long n;
  const migr_table_fmt *mig;
  double end;
} mlh_path;

static mlh_list mlh_R2;      /* R' */
static mlh_list mlh_Gold;    /* predicted G above the origin */
static mlh_list mlh_Gnew;    /* predicted G' above the origin */
static migr_table_fmt *mlh_oldmig = NULL;
static long mlh_oldmig_alloc = 0;
static migr_table_fmt *mlh_dropped = NULL;
static long mlh_dropped_alloc = 0;

static void mlh_reserve(mlh_list *L, long n, long numpop)
{
  if (L->alloc < n || L->numpop != numpop)
    {
      long a = n + 64;
      L->ev = (mlh_event *) myrealloc(L->ev, (size_t) a * sizeof(mlh_event));
      L->k = (long *) myrealloc(L->k, (size_t) (a * numpop) * sizeof(long));
      L->alloc = a;
      L->numpop = numpop;
    }
}

static void mlh_push(mlh_list *L, double age, char type, long below, long above,
                     const long *k)
{
  mlh_reserve(L, L->n + 1, L->numpop);
  L->ev[L->n].age = age;
  L->ev[L->n].type = type;
  L->ev[L->n].below = below;
  L->ev[L->n].above = above;
  memcpy(L->k + L->n * L->numpop, k, (size_t) L->numpop * sizeof(long));
  L->n++;
}

static char mlh_type(char nodetype)
{
  switch (nodetype)
    {
    case 'i':
      return 'c';
    case 'm':
    case 'd':
      return 'm';
    default:
      return 't';
    }
}

/* the events of a timelist strictly above age s */
static void mlh_from_tl(mlh_list *L, timelist_fmt *tv, double s, long numpop)
{
  long i;
  L->n = 0;
  mlh_reserve(L, tv->T, numpop);
  L->rootpop = tv->tl[tv->T - 2].eventnode->pop;
  for (i = 0; i <= tv->T - 2; i++)
    {
      vtlist *t = &tv->tl[i];
      if (t->age <= s)
        continue;
      mlh_push(L, t->age, mlh_type(t->eventnode->type),
               t->eventnode->actualpop, t->eventnode->pop, t->lineages);
    }
}

/* ---- waiting-time densities ---- */
static double mlh_logS(double lambda, double alpha, double t)
{
  if (t <= 0.0 || lambda <= 0.0)
    return 0.0;
  if (alpha >= 1.0)
    return -lambda * t;
  return creal(mittag_leffler(alpha, 1.0, -lambda * pow(t, alpha)));
}

static double mlh_logf(double lambda, double alpha, double t)
{
  if (alpha >= 1.0)
    return log(lambda) - lambda * t;
  return log(lambda) + (alpha - 1.0) * log(t)
    + creal(mittag_leffler(alpha, alpha, -lambda * pow(t, alpha)));
}

static double mlh_alpha(world_fmt *world, long pop)
{
  if (world->has_mlalpha && world->options->mlalphapops[pop] != 0)
    return world->mlalpha[world->options->mlalphapops[pop] - 1];
  return 1.0;
}

/* pair coalescence rate in pop: 2/(mu theta) */
static double mlh_pairrate(world_fmt *world, long pop)
{
  return 2.0 / (world->options->mu_rates[world->locus]
                * world->param0[pop] * world->timek[pop]);
}

/* migration rate of one lineage in pop 'below' moving to 'above' (backwards
   in time), as in time_to_migration() */
static double mlh_migrate(world_fmt *world, long above, long below)
{
  long i = m2mmm(above, below, world->numpop);
  if (world->options->custm2[i] == '0')
    return 0.0;
  return world->data->geo[i] * world->param0[i] * world->timek[i]
    / world->options->mu_rates[world->locus];
}

/* total migration rate out of pop for one lineage */
static double mlh_migtotal(world_fmt *world, long pop)
{
  long i;
  double sum = 0.0;
  for (i = world->mstart[pop]; i < world->mend[pop]; i++)
    {
      if (world->options->custm2[i] == '0')
        continue;
      sum += world->data->geo[i] * world->param0[i] * world->timek[i];
    }
  return sum / world->options->mu_rates[world->locus];
}

/* log survival of all migration clocks of a lineage in pop over dt */
static double mlh_logS_mig(world_fmt *world, long pop, double dt)
{
  long i;
  double alpha = mlh_alpha(world, pop), logs = 0.0;
  const double mu = world->options->mu_rates[world->locus];
  if (alpha >= 1.0)
    return -mlh_migtotal(world, pop) * dt;
  for (i = world->mstart[pop]; i < world->mend[pop]; i++)
    {
      if (world->options->custm2[i] == '0')
        continue;
      logs += mlh_logS(world->data->geo[i] * world->param0[i] * world->timek[i] / mu,
                       alpha, dt);
    }
  return logs;
}

#ifdef MIGRATE_MLH_VERIFY
/* ---- the interval density of probg_treetimes() before the per-population
   model (one alpha per whole-tree interval); kept as a reference: equal to
   mlh_logp_perpop() for one population or alpha = 1 ---- */
static double mlh_interval(world_fmt *world, double dt, const long *k,
                           const mlh_event *e)
{
  const long numpop = world->numpop;
  const double alpha = mlh_alpha(world, e->below);
  double lam = 0.0, rate;
  long pop;
  if (e->type == 't')
    return 0.0;
  for (pop = 0; pop < numpop; pop++)
    {
      if (k[pop] > 1)
        lam += k[pop] * (k[pop] - 1) * mlh_pairrate(world, pop) / 2.0;
      if (k[pop] > 0)
        lam += k[pop] * mlh_migtotal(world, pop);
    }
  rate = (e->type == 'c') ? mlh_pairrate(world, e->below)
    : mlh_migrate(world, e->above, e->below);
  if (alpha >= 1.0)
    return log(rate) - lam * dt;
  return (alpha - 1.0) * log(dt)
    + creal(mittag_leffler(alpha, alpha, -lam * pow(dt, alpha))) + log(rate);
}

static double mlh_logp_list(world_fmt *world, const mlh_list *L, double s)
{
  long i;
  double age = s, logp = 0.0;
  for (i = 0; i < L->n; i++)
    {
      logp += mlh_interval(world, L->ev[i].age - age, L->k + i * L->numpop, &L->ev[i]);
      age = L->ev[i].age;
    }
  return logp;
}
#endif

/* ---- per-population tree density (docs/mittag_leffler_in_migrate.tex,
   section "Per-population alpha") ----
   Every population runs its own ML clock with its own alpha. The clock of
   pop carries its coalescences and the migrations of its lineages
   (backwards in time; forwards these are immigrations into pop) and runs
   over a spell, the stretch during which k[pop] does not change. A spell
   ends with an event of pop, or is censored when k[pop] changes for another
   reason (a lineage arriving, a sample); either way the clock restarts. */

/* total rate of pop's clock with kp lineages */
static double mlh_lambda(world_fmt *world, long pop, long kp)
{
  double lam = 0.0;
  if (kp > 1)
    lam += kp * (kp - 1) * mlh_pairrate(world, pop) / 2.0;
  if (kp > 0)
    lam += kp * mlh_migtotal(world, pop);
  return lam;
}

/* log density of the events of L (and of everything L's populations wait
   through) given the start a[pop] of each population's spell open below the
   first event of L; a[] is updated. The process stops at the last event. */
static double mlh_logp_perpop(world_fmt *world, const mlh_list *L, double *a)
{
  const long numpop = L->numpop;
  long i, pop, x;
  double logp = 0.0;
  for (i = 0; i < L->n; i++)
    {
      const mlh_event *e = &L->ev[i];
      const long *kb = L->k + i * numpop;
      const long *ka = (i + 1 < L->n) ? L->k + (i + 1) * numpop : NULL;
      x = -1;
      if (e->type == 'c' || e->type == 'm')
        {
          double alpha, dt, lam, rate;
          x = e->below;
          alpha = mlh_alpha(world, x);
          dt = e->age - a[x];
          lam = mlh_lambda(world, x, kb[x]);
          rate = (e->type == 'c') ? mlh_pairrate(world, x)
            : mlh_migrate(world, e->above, x);
          /* f(dt) * rate / lam */
          if (alpha >= 1.0)
            logp += log(rate) - lam * dt;
          else
            logp += (alpha - 1.0) * log(dt)
              + creal(mittag_leffler(alpha, alpha, -lam * pow(dt, alpha))) + log(rate);
          a[x] = e->age;
        }
      for (pop = 0; pop < numpop; pop++)
        {
          if (pop == x || (ka != NULL && ka[pop] == kb[pop]) || (ka == NULL && kb[pop] == 0))
            continue;
          if (ka != NULL)   /* censored: k[pop] changed by another population's event */
            logp += mlh_logS(mlh_lambda(world, pop, kb[pop]), mlh_alpha(world, pop),
                             e->age - a[pop]);
          a[pop] = e->age;
        }
    }
  return logp;
}

/* spell starts of all populations at time s, from a timelist whose events
   at or below s are those of the tree (R or G); pop 'atpop' has an event at
   s itself (the origin) */
static void mlh_spellstarts(timelist_fmt *tv, double s, long numpop, long atpop,
                            double *a)
{
  long i, pop;
  for (pop = 0; pop < numpop; pop++)
    a[pop] = 0.0;
  for (i = 0; i + 1 <= tv->T - 2 && tv->tl[i].age < s; i++)
    for (pop = 0; pop < numpop; pop++)
      if (tv->tl[i + 1].lineages[pop] != tv->tl[i].lineages[pop])
        a[pop] = tv->tl[i].age;
  if (atpop >= 0)
    a[atpop] = s;
}

/* probg_treetimes() for Mittag-Leffler runs within mlh_supported() */
double mlh_probg_treetimes(world_fmt *world)
{
  static mlh_list L;
  static double *a = NULL;
  static long aalloc = 0;
  const long numpop = world->numpop;
  long pop;
  if (aalloc < numpop)
    {
      a = (double *) myrealloc(a, (size_t) numpop * sizeof(double));
      aalloc = numpop;
    }
  mlh_from_tl(&L, world->treetimes, -1.0, numpop);
  for (pop = 0; pop < numpop; pop++)
    a[pop] = (L.n > 0) ? L.ev[0].age : 0.0;
#ifdef MIGRATE_MLH_VERIFY
  {
    static double worst = 0.0;
    boolean alpha1 = TRUE;
    const double a0 = a[0];
    double pp = mlh_logp_perpop(world, &L, a);
    double pi = mlh_logp_list(world, &L, a0);
    for (pop = 0; pop < numpop; pop++)
      if (mlh_alpha(world, pop) < 1.0)
        alpha1 = FALSE;
    if ((numpop == 1 || alpha1) && fabs(pp - pi) > worst)
      {
        worst = fabs(pp - pi);
        fprintf(stderr, "MLHVERIFY probg per-pop vs interval density: worst |diff| = %g (logp %g)\n",
                worst, pp);
      }
    return pp;
  }
#else
  return mlh_logp_perpop(world, &L, a);
#endif
}

/* population of a path just below time t */
static long mlh_pathpop(const mlh_path *P, double t)
{
  long j, pop = P->startpop;
  for (j = 0; j < P->n && P->mig[j].time < t; j++)
    pop = P->mig[j].from;
  return pop;
}

/* G = background B (above s) + moving path P coalescing at P->end;
   B already holds the root lineage's events (drops and additions applied);
   beyond the last event of B only the root lineage remains */
static void mlh_merge(mlh_list *G, const mlh_list *B, const mlh_path *P,
                      long numpop)
{
  long ib = 0, im = 0, pop;
  long *k = (long *) mycalloc(numpop, sizeof(long));
  const long rootpop = B->rootpop;
  boolean done = FALSE;
  G->n = 0;
  mlh_reserve(G, B->n + P->n + 2, numpop);
  for (;;)
    {
      const double tb = (ib < B->n) ? B->ev[ib].age : (double) HUGE;
      const double tm = done ? (double) HUGE : ((im < P->n) ? P->mig[im].time : P->end);
      if (ib >= B->n && done)
        break;
      if (tm < tb)
        {                               /* event of the moving lineage */
          if (ib < B->n)
            memcpy(k, B->k + ib * numpop, (size_t) numpop * sizeof(long));
          else
            {
              memset(k, 0, (size_t) numpop * sizeof(long));
              k[rootpop] = 1;
            }
          pop = mlh_pathpop(P, tm);
          k[pop] += 1;
          if (im < P->n)
            {
              mlh_push(G, tm, 'm', P->mig[im].to, P->mig[im].from, k);
              im++;
            }
          else
            {
              mlh_push(G, tm, 'c', pop, pop, k);
              done = TRUE;
            }
          continue;
        }
      memcpy(k, B->k + ib * numpop, (size_t) numpop * sizeof(long));
      if (!done)
        k[mlh_pathpop(P, tb)] += 1;
      mlh_push(G, tb, B->ev[ib].type, B->ev[ib].below, B->ev[ib].above, k);
      ib++;
    }
  myfree(k);
}

/* ---- proposal density of a path ---- */
/* first entry of B above time t, with newtree_update()'s SMALLEPSILON rule */
static long mlh_first(const mlh_list *B, double t)
{
  long i = 0;
  while (i < B->n && (B->ev[i].age < t || B->ev[i].age - t < SMALLEPSILON))
    i++;
  return i;
}

/* log density of all clocks of a lineage in pop firing 'which' at dt:
   which = -1 coalescence, -2 none (survival only), else a migration to
   population 'which' (backwards) */
static double mlh_clocks(world_fmt *world, long pop, long kcoal, long which,
                         double dt)
{
  const double alpha = mlh_alpha(world, pop);
  const double lc = (kcoal > 0) ? kcoal * mlh_pairrate(world, pop) : 0.0;
  double logq;
  if (which == -1)
    logq = mlh_logf(lc, alpha, dt);
  else
    logq = mlh_logS(lc, alpha, dt);
  if (which >= 0)
    {
      double lm = mlh_migrate(world, which, pop);
      if (alpha >= 1.0)
        logq += log(lm) - mlh_migtotal(world, pop) * dt;
      else
        logq += mlh_logS_mig(world, pop, dt) - mlh_logS(lm, alpha, dt)
          + mlh_logf(lm, alpha, dt);
    }
  else
    logq += mlh_logS_mig(world, pop, dt);
  return logq;
}

/* q(P, Q | B, h): moving path P on background B (the residual tree), root
   lineage events Q above the horizon h */
static double mlh_path_logq(world_fmt *world, const mlh_list *B,
                            const mlh_path *P, const migr_table_fmt *Q, long nq,
                            double h)
{
  double age = P->start, e;
  long pop = P->startpop, im = 0, iq = 0, i, kp;
  long p2 = B->rootpop;
  double logq = 0.0;
  /* inside the residual tree: one slice per entry of B */
  for (i = mlh_first(B, P->start); i < B->n; i++)
    {
      const double b = B->ev[i].age;
      for (;;)
        {
          kp = B->k[i * B->numpop + pop];
          e = (im < P->n) ? P->mig[im].time : P->end;
          if (e >= b)
            break;
          if (im < P->n)
            {
              logq += mlh_clocks(world, pop, kp, P->mig[im].from, e - age);
              pop = P->mig[im].from;
              im++;
              age = e;
              continue;
            }
          if (kp < 1)
            return (double) -HUGE;
          return logq + mlh_clocks(world, pop, kp, -1, e - age) - log((double) kp);
        }
      logq += mlh_clocks(world, pop, kp, -2, b - age);
      age = b;
    }
  /* phase 1: the root lineage waits in p2 up to the horizon */
  while (age < h)
    {
      e = (im < P->n) ? P->mig[im].time : P->end;
      if (e >= h)
        {
          logq += mlh_clocks(world, pop, pop == p2 ? 1 : 0, -2, h - age);
          age = h;
          break;
        }
      if (im < P->n)
        {
          logq += mlh_clocks(world, pop, pop == p2 ? 1 : 0, P->mig[im].from, e - age);
          pop = P->mig[im].from;
          im++;
          age = e;
          continue;
        }
      if (pop != p2)
        return (double) -HUGE;
      return logq + mlh_clocks(world, pop, 1, -1, e - age);
    }
  if (nq > 0 && Q[0].time < age)
    return (double) -HUGE;   /* the root lineage cannot move below h */
  /* phase 2: both lineages, all clocks restart at every event */
  for (;;)
    {
      double e1 = (im < P->n) ? P->mig[im].time : P->end;
      double e2 = (iq < nq) ? Q[iq].time : (double) HUGE;
      if (e1 < e2)
        {
          logq += mlh_clocks(world, p2, 0, -2, e1 - age);
          if (im < P->n)
            {
              logq += mlh_clocks(world, pop, pop == p2 ? 1 : 0, P->mig[im].from, e1 - age);
              pop = P->mig[im].from;
              im++;
              age = e1;
              continue;
            }
          if (pop != p2)
            return (double) -HUGE;
          return logq + mlh_clocks(world, pop, 1, -1, e1 - age);
        }
      logq += mlh_clocks(world, pop, pop == p2 ? 1 : 0, -2, e2 - age);
      logq += mlh_clocks(world, p2, 0, Q[iq].from, e2 - age);
      p2 = Q[iq].from;
      iq++;
      age = e2;
    }
}

/* ---- the move ---- */
static void mlh_push_mig(migr_table_fmt **tab, long *alloc, long n,
                         double time, long from, long to)
{
  if (n >= *alloc)
    {
      *alloc = n + 32;
      *tab = (migr_table_fmt *) myrealloc(*tab, (size_t) *alloc * sizeof(migr_table_fmt));
    }
  (*tab)[n].time = time;
  (*tab)[n].from = from;
  (*tab)[n].to = to;
  (*tab)[n].event = 'm';
}

/* timeelements == 2 is the no-skyline default (one open segment) */
boolean mlh_supported(world_fmt *world)
{
  return world->timeelements <= 2 && !world->has_speciation
    && !world->has_growth && !world->options->has_datefile;
}

/* population of the reassigned tip before an assignment move */
static long mlh_assign_oldpop = -1;

void mlh_set_assignment(long oldpop)
{
  mlh_assign_oldpop = oldpop;
}

/* log of p(G')/p(G) * q(G | G')/q(G' | G) for the proposal that is about to
   be decided in acceptlike(); R is the residual timelist. For an assignment
   move the tip's lineage is not part of R, so the same formula applies with
   the old path starting in the tip's old population (the choice of the
   individual and of the new population is symmetric) */
double mlh_log_correction(world_fmt *world, proposal_fmt *proposal,
                          timelist_fmt *R, boolean assign)
{
  static mlh_list B;
  const long numpop = world->numpop;
  const double s = proposal->origin->tyme;
  const double t_old = proposal->oback->tyme;
  const double t_new = proposal->time;
  long nold = 0, ndrop = 0, i, ctop;
  double top, top2, h_f, h_r;
  mlh_path Pnew, Pold;
  node *p;
  double logp_new, logp_old, logq_new, logq_old;

  mlh_from_tl(&B, R, s, numpop);
  /* the old path: migration nodes on the branch from the origin to oback */
  for (p = proposal->origin->back; p->type == 'm' || p->type == 'd'; p = p->next->back)
    {
      node *pt = showtop(p);
      mlh_push_mig(&mlh_oldmig, &mlh_oldmig_alloc, nold++, pt->tyme, pt->pop, pt->actualpop);
    }
  Pold.start = Pnew.start = s;
  Pnew.startpop = proposal->origin->pop;
  Pold.startpop = assign ? mlh_assign_oldpop : proposal->origin->pop;
  Pold.n = nold;
  Pold.mig = mlh_oldmig;
  Pold.end = t_old;
  Pnew.n = proposal->migr_table_counter;
  Pnew.mig = proposal->migr_table;
  Pnew.end = t_new;

  top = (B.n > 0) ? B.ev[B.n - 1].age : s;
  h_f = (t_old > top) ? t_old : top;

  /* G: R plus the old path */
  mlh_merge(&mlh_Gold, &B, &Pold, numpop);

  /* D: a tree has no migration nodes above its top coalescence, so the
     root lineage's events of R above max(t_new, top coalescence of R) are
     removed */
  ctop = B.n - 1;
  while (ctop >= 0 && B.ev[ctop].type == 'm')
    ctop--;
  mlh_R2.n = 0;
  mlh_reserve(&mlh_R2, B.n + proposal->migr_table_counter2 + 1, numpop);
  mlh_R2.rootpop = B.rootpop;
  for (i = 0; i < B.n; i++)
    {
      if (i > ctop && B.ev[i].age > t_new)
        {
          if (ndrop == 0)
            mlh_R2.rootpop = B.ev[i].below;
          mlh_push_mig(&mlh_dropped, &mlh_dropped_alloc, ndrop++,
                       B.ev[i].age, B.ev[i].above, B.ev[i].below);
        }
      else
        mlh_push(&mlh_R2, B.ev[i].age, B.ev[i].type, B.ev[i].below,
                 B.ev[i].above, B.k + i * numpop);
    }
  /* P2: the root lineage's new migrations above the horizon */
  {
    long *k = (long *) mycalloc(numpop, sizeof(long));
    for (i = 0; i < proposal->migr_table_counter2; i++)
      {
        memset(k, 0, (size_t) numpop * sizeof(long));
        k[proposal->migr_table2[i].to] = 1;
        mlh_push(&mlh_R2, proposal->migr_table2[i].time, 'm',
                 proposal->migr_table2[i].to, proposal->migr_table2[i].from, k);
        mlh_R2.rootpop = proposal->migr_table2[i].from;
      }
    myfree(k);
  }
  /* G': R' plus the new path */
  mlh_merge(&mlh_Gnew, &mlh_R2, &Pnew, numpop);

  top2 = (mlh_R2.n > 0) ? mlh_R2.ev[mlh_R2.n - 1].age : s;
  h_r = (t_new > top2) ? t_new : top2;

  /* spells open at s started below s and are the same in G and G' (the
     tip of an assignment move lies at time 0, where all spells start) */
  {
    static double *a = NULL;
    static long aalloc = 0;
    if (aalloc < numpop)
      {
        a = (double *) myrealloc(a, (size_t) numpop * sizeof(double));
        aalloc = numpop;
      }
    mlh_spellstarts(R, s, numpop, Pold.startpop, a);
    logp_old = mlh_logp_perpop(world, &mlh_Gold, a);
    mlh_spellstarts(R, s, numpop, Pnew.startpop, a);
    logp_new = mlh_logp_perpop(world, &mlh_Gnew, a);
  }
  logq_new = mlh_path_logq(world, &B, &Pnew, proposal->migr_table2,
                           proposal->migr_table_counter2, h_f);
  logq_old = mlh_path_logq(world, &mlh_R2, &Pold, mlh_dropped, ndrop, h_r);
#ifdef MIGRATE_MLH_VERIFY
  mlh_verify_record(world, logp_new - logp_old, logq_old - logq_new);
#endif
  if (logq_old <= -HUGE / 2.0)
    return (double) -HUGE;
  return (logp_new - logp_old) + (logq_old - logq_new);
}

#ifdef MIGRATE_MLH_VERIFY
/* ---- verification harness: compare the predicted lists with the trees
   the move actually builds, and check that the correction is 1 for
   alpha = 1 ---- */
static double mlh_worst_alpha1 = 0.0;
static long mlh_nrecord = 0, mlh_nbad_g = 0, mlh_nbad_gold = 0, mlh_nbad_r = 0;

void mlh_verify_record(world_fmt *world, double dlogp, double dlogq)
{
  long pop;
  boolean alpha1 = TRUE;
  for (pop = 0; pop < world->numpop; pop++)
    if (mlh_alpha(world, pop) < 1.0)
      alpha1 = FALSE;
  mlh_nrecord++;
  if (alpha1 && fabs(dlogp + dlogq) > mlh_worst_alpha1)
    {
      mlh_worst_alpha1 = fabs(dlogp + dlogq);
      fprintf(stderr, "MLHVERIFY alpha=1 worst |dlogp+dlogq| = %g (dlogp %g) after %li moves\n",
              mlh_worst_alpha1, dlogp, mlh_nrecord);
    }
}

static boolean mlh_same(const mlh_list *A, const mlh_list *B, const char *what,
                        long *nbad)
{
  long i, pop;
  boolean ok = (A->n == B->n);
  for (i = 0; ok && i < A->n; i++)
    {
      if (fabs(A->ev[i].age - B->ev[i].age) > 1e-12 * (1.0 + fabs(A->ev[i].age))
          || A->ev[i].type != B->ev[i].type
          || (A->ev[i].type != 't' && (A->ev[i].below != B->ev[i].below
                                       || A->ev[i].above != B->ev[i].above)))
        ok = FALSE;
      for (pop = 0; ok && pop < A->numpop; pop++)
        if (A->k[i * A->numpop + pop] != B->k[i * B->numpop + pop])
          ok = FALSE;
    }
  if (!ok && (*nbad)++ < 20)
    {
      fprintf(stderr, "MLHVERIFY %s mismatch (predicted %li events, actual %li)\n",
              what, A->n, B->n);
      for (i = 0; i < A->n || i < B->n; i++)
        {
          if (i < A->n)
            fprintf(stderr, "  P %.10f %c %li->%li k0=%li", A->ev[i].age, A->ev[i].type,
                    A->ev[i].below, A->ev[i].above, A->k[i * A->numpop]);
          else
            fprintf(stderr, "  P -");
          if (i < B->n)
            fprintf(stderr, "   A %.10f %c %li->%li k0=%li\n", B->ev[i].age, B->ev[i].type,
                    B->ev[i].below, B->ev[i].above, B->k[i * B->numpop]);
          else
            fprintf(stderr, "   A -\n");
        }
    }
  return ok;
}

/* before the move: the predicted G must be the current tree */
void mlh_verify_before(world_fmt *world, double s)
{
  static mlh_list A;
  mlh_from_tl(&A, &world->treetimes[0], s, world->numpop);
  mlh_same(&mlh_Gold, &A, "G", &mlh_nbad_gold);
}

/* after an accepted move: the predicted G' and R' must be the new tree and
   its residual for the same origin */
void mlh_verify_after(world_fmt *world, proposal_fmt *R2prop, timelist_fmt *R2, double s)
{
  static mlh_list A;
  (void) R2prop;
  mlh_from_tl(&A, &world->treetimes[0], s, world->numpop);
  mlh_same(&mlh_Gnew, &A, "G'", &mlh_nbad_g);
  mlh_from_tl(&A, R2, s, world->numpop);
  mlh_same(&mlh_R2, &A, "R'", &mlh_nbad_r);
  if (mlh_nrecord % 100000 == 0)
    fprintf(stderr, "MLHVERIFY %li moves: G mismatches %li, G' %li, R' %li\n",
            mlh_nrecord, mlh_nbad_gold, mlh_nbad_g, mlh_nbad_r);
}
#endif
