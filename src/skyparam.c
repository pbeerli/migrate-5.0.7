// skyline plot parametrized
#include "migration.h"
#include "world.h"
#include "tree.h"
#include "slice.h"
#include "tools.h"
#include "skyparam.h"
#include "random.h"
#include "sighandler.h"
#include "bayes.h"
#include "priors.h"
#include "pretty.h"
#include "joint_combine.h"
#include <math.h>

/* skyline=PARAM multipliers world->timek[s * numpop2 + k] scale Theta_k
   (k < numpop) or M (custm2 order) in segment s; timek of the first segment
   is 1, the others have the prior of bayes-priors= SKYLINE (default: random
   walk on log tau, sigma 1), see docs/mittag_leffler_in_migrate.tex */
/* half-width of the symmetric step on log tau */
#define SKYLINE_TIMEPARAM_STEP 0.5

#ifdef MPI
void      skyline_param_record(float *temp, long *z, world_fmt * world);
#else
void      skyline_param_record(char *temp, long *c, world_fmt * world);
#endif

//##
/* log prior of the log multipliers of parameter k (segments 1.., the first
   segment's multiplier is 1): random walk on log tau, or log-uniform */
static MYREAL skyline_logprior(world_fmt *world, long k)
{
  const long nseg = world->timeelements - 1;
  const long numpop2 = world->numpop2;
  MYREAL lp = 0.0, prev = 0.0, x;
  long s;
  for (s = 1; s < nseg; s++)
    {
      x = LOG(world->timek[s * numpop2 + k]);
      if (world->skyptype[k] == SKYPRIOR_LOGUNIFORM)
        {
          if (x < LOG(world->skypa[k]) || x > LOG(world->skypb[k]))
            return (MYREAL) -HUGE;
        }
      else
        {
          const MYREAL sd = world->skypa[k];
          lp += -0.5 * (x - prev) * (x - prev) / (sd * sd);
        }
      prev = x;
    }
  return lp;
}

/// Metropolis update of one skyline=PARAM multiplier: a free parameter
/// ('*' in skyline-model) and a segment other than the first (anchored at 1)
/// are chosen uniformly; the step is symmetric in log tau and the prior is
/// on log tau, so the acceptance is the p(G) ratio times the prior ratio.
MYREAL bayes_update_timeparam(world_fmt * world, boolean *success)
{
  const long nseg = world->timeelements - 1;
  const long numpop2 = world->numpop2;
  long nfree = 0, k, pick, s, idx;
  MYREAL oldval, newval, oldtau, oldprior, newprior;
  *success = FALSE;
  for (k = 0; k < numpop2; k++)
    if (world->skymodel2[k] == '*')
      nfree++;
  if (nseg < 2 || nfree == 0)
    return world->param_like;
  pick = RANDINT(0, nfree - 1);
  for (k = 0; k < numpop2; k++)
    if (world->skymodel2[k] == '*' && pick-- == 0)
      break;
  s = RANDINT(1, nseg - 1);
  idx = s * numpop2 + k;
  world->sky_trials++;
  oldval = probg_treetimes(world);
  oldtau = world->timek[idx];
  oldprior = skyline_logprior(world, k);
  world->timek[idx] = EXP(LOG(oldtau) + SKYLINE_TIMEPARAM_STEP * (2.0 * UNIF_RANDUM() - 1.0));
  /* hard limits exp(+-SKYLINE_HISTMAX): the random-walk prior is unbounded,
     and heated chains (which heat the prior too) otherwise let a multiplier
     drift until the coalescence rates vanish and the times overflow */
  if (fabs(LOG(world->timek[idx])) > SKYLINE_HISTMAX)
    {
      world->timek[idx] = oldtau;
      skyline_param_count(world, FALSE);
      return oldval;
    }
  newprior = skyline_logprior(world, k);
  if (newprior <= (MYREAL) -HUGE / 2.0)
    {
      world->timek[idx] = oldtau;
      skyline_param_count(world, FALSE);
      return oldval;
    }
  newval = probg_treetimes(world);
  /* p(G|params) and prior: not heated under the power posterior */
  *success = bayes_accept(newval, oldval, 1.0, newprior - oldprior);
  skyline_param_count(world, *success);
  if (*success)
    {
      world->sky_accept++;
      world->bayes->oldval = newval;
      return newval;
    }
  world->timek[idx] = oldtau;
  return oldval;
}


/* slice_bayes_update_timeparam() -- the slice-sampling counterpart to
   bayes_update_timeparam() above -- was removed 2026-09-08: confirmed
   zero callers anywhere in the tree (skyline=PARAM's proposal always
   goes through the Metropolis path above), and it had its own,
   independent, genuine out-of-bounds bug: it indexed
   bayes->minparam[]/maxparam[] (both sized to world->numparam, the
   real Bayesian-parameter count) with `wb = whichblock*npx + w`, a
   combined time-slice*parameter index that can exceed that allocation
   for any blocks/timeelements count above 1 -- a real out-of-bounds
   read, not just the wrong-bounds *value* bug fixed above, and with no
   caller to ever trigger or verify a fix against. Removed outright
   rather than left as an unreachable landmine, matching the mixfile/
   write-summary/long-chains/mu-rates precedent for confirmed-dead
   code. */

void insert_time_boundaries(timelist_fmt *timevector, world_fmt *world)
{
  MYREAL *times = world->times;
  long timeelements = world->timeelements;
  long start = timevector->T;
  long z;
  long i;
  if(timeelements>1)
    {
      if (timevector->allocT <= timevector->T + timeelements-1)
	{
	  increase_timelist2 (&timevector, timeelements); 
	}
      timevector->T += timeelements;
      z = 1;
      for(i=start; i < timevector->T; i++)
	{
	  timevector->tl[i].age = times[z];
	  timevector->tl[i].eventnode = NULL;
	  timevector->tl[i].slice = i;
	  timevector->tl[i].from = -1;
	  timevector->tl[i].to = -1;
	} 
    }
}

#ifdef MPI
void      skyline_param_record(float *temp, long *z, world_fmt * world)
{
  long i;
  long addition=0;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  for (i=0; i < end; i++)
    {
      temp[*z] = (float) world->times[i]; //times and timek 
      *z += 1;
    }
}
#else
void      skyline_param_record(char *temp, long *c, world_fmt * world)
{
  long i;
  long addition=0;
  long end = world->timeelements +  world->timeelements * (world->numpop2 + addition);
  for (i=0; i < end; i++)
    *c += mysnprintf(temp+ *c,LINESIZE,"\t%f", world->times[i]);//times and timek 
}
#endif

/* ---- recording and reporting of the skyline=PARAM multipliers ---- */

extern int myID;

/* entries: every free parameter ('*') in every segment but the first */
void skyline_param_alloc(world_fmt *world, long loci)
{
  const long nseg = world->timeelements - 1;
  long k, s, e = 0, l;
  world->sky_nent = 0;
  world->sky_ng = 0;
  world->sky_loci = loci;
  world->skyloc = NULL;
  world->sky_k = world->sky_s = NULL;
  world->sky2sum = world->sky2xlo = NULL;
  world->sky2loci = NULL;
  world->sky_gk = NULL;
  world->sky2gsum = world->sky2gxlo = NULL;
  world->sky2gloci = NULL;
  world->sky_lastlocus = -1;
  for (k = 0; k < world->numpop2; k++)
    if (world->skymodel2[k] == '*')
      world->sky_nent += nseg - 1;
  for (k = 0; k < world->numpop; k++)
    if (world->skymodel2[k] == 'g')
      world->sky_ng++;
  if (world->sky_nent == 0 && world->sky_ng == 0)
    return;
  if (world->sky_ng > 0)
    {
      long j = 0;
      world->sky_gk = (long *) mycalloc(world->sky_ng, sizeof(long));
      for (k = 0; k < world->numpop; k++)
        if (world->skymodel2[k] == 'g')
          world->sky_gk[j++] = k;
      world->sky2gsum = (double *) mycalloc(world->sky_ng * SKY2_NX * SKY2_NG, sizeof(double));
      world->sky2gxlo = (double *) mycalloc(world->sky_ng, sizeof(double));
      world->sky2gloci = (long *) mycalloc(world->sky_ng, sizeof(long));
      for (j = 0; j < world->sky_ng; j++)
        world->sky2gxlo[j] = HUGE_VAL;   /* set from the prior at first use */
    }
  world->skyloc = (skyloc_fmt *) mycalloc(loci, sizeof(skyloc_fmt));
  for (l = 0; l < loci; l++)
    if (world->sky_ng > 0)
      {
        world->skyloc[l].gess = (double *) mycalloc(world->sky_ng * nseg, sizeof(double));
        world->skyloc[l].gq = (double *) mycalloc(6 * world->sky_ng * nseg, sizeof(double));
      }
  if (world->sky_nent == 0)
    return;
  world->sky_k = (long *) mycalloc(world->sky_nent, sizeof(long));
  world->sky_s = (long *) mycalloc(world->sky_nent, sizeof(long));
  for (k = 0; k < world->numpop2; k++)
    if (world->skymodel2[k] == '*')
      for (s = 1; s < nseg; s++)
        {
          world->sky_k[e] = k;
          world->sky_s[e] = s;
          e++;
        }
  world->sky2sum = (double *) mycalloc(world->sky_nent * SKY2_NX * SKY2_NY, sizeof(double));
  world->sky2xlo = (double *) mycalloc(world->sky_nent, sizeof(double));
  world->sky2loci = (long *) mycalloc(world->sky_nent, sizeof(long));
  for (e = 0; e < world->sky_nent; e++)
    world->sky2xlo[e] = HUGE_VAL;   /* set from the prior at first use */
  for (l = 0; l < loci; l++)
    {
      world->skyloc[l].hist = (double *) mycalloc(world->sky_nent * SKYLINE_HISTBINS, sizeof(double));
      world->skyloc[l].ess = (double *) mycalloc(world->sky_nent, sizeof(double));
      world->skyloc[l].vq = (double *) mycalloc(3 * world->sky_nent, sizeof(double));
    }
}

/* index of parameter k in the prior arrays (-1 if not estimated); an
   entry that shares a parameter (symmetric, mean) has its prior too */
static long skyline_prior_index(world_fmt *world, long k)
{
  const long i = world->bayes->map[k][1];
  return (i == INVALID) ? -1 : i;
}

/* the log-parameter axis spans SKY2_XRANGE below the log of the prior maximum */
static double skyline_xlo_k(world_fmt *world, long k)
{
  long i = skyline_prior_index(world, k);
  double mx = (i >= 0) ? world->bayes->maxparam[i] : 1.0;
  return log(mx > 0.0 ? mx : 1.0) - SKY2_XRANGE;
}

static double skyline_xlo(world_fmt *world, long e)
{
  if (world->sky2xlo[e] == HUGE_VAL)
    world->sky2xlo[e] = skyline_xlo_k(world, world->sky_k[e]);
  return world->sky2xlo[e];
}

/* growth-linked Theta_p: the growth parameter of p (prior index and value);
   the growth axis of its grid is the growth prior range */
static long skyline_growth_index(world_fmt *world, long p)
{
  return world->numparamcumvec[SPLITSTDPRIOR] + world->options->growpops[p] - 1;
}

static double skyline_growth(world_fmt *world, long p)
{
  return world->growth[world->options->growpops[p] - 1];
}

static double skyline_gxlo(world_fmt *world, long j)
{
  if (world->sky2gxlo[j] == HUGE_VAL)
    world->sky2gxlo[j] = skyline_xlo_k(world, world->sky_gk[j]);
  return world->sky2gxlo[j];
}

static void skyline_grange(world_fmt *world, long j, double *lo, double *width)
{
  const long gi = skyline_growth_index(world, world->sky_gk[j]);
  *lo = world->bayes->minparam[gi];
  *width = world->bayes->maxparam[gi] - *lo;
  if (!(*width > 0.0))
    *width = 1.0;
}

static long skyline_bin(double x)
{
  long b = (long) floor((x + SKYLINE_HISTMAX) / (2.0 * SKYLINE_HISTMAX) * SKYLINE_HISTBINS);
  if (b < 0)
    b = 0;
  if (b >= SKYLINE_HISTBINS)
    b = SKYLINE_HISTBINS - 1;
  return b;
}

static long skyline_bin2(double v, double lo, double width, long n)
{
  long b = (long) floor((v - lo) / width * n);
  if (b < 0)
    b = 0;
  if (b >= n)
    b = n - 1;
  return b;
}

static void skyline_joint_add(world_fmt *world, long locus);

/* one posterior sample of the cold chain */
void skyline_param_sample(world_fmt *world)
{
  skyloc_fmt *L;
  long e, j;
  if ((world->sky_nent == 0 && world->sky_ng == 0) || !world->cold || world->locus >= world->sky_loci)
    return;
  /* serial runs: a locus is complete once sampling moves on */
  if (myID == MASTER && world->sky_lastlocus >= 0 && world->sky_lastlocus != world->locus)
    skyline_joint_add(world, world->sky_lastlocus);
  world->sky_lastlocus = world->locus;
  L = &world->skyloc[world->locus];
  if (L->n + 1 > L->alloc)
    {
      L->alloc = 2 * L->alloc + 1024;
      if (world->sky_nent > 0)
        L->samples = (float *) myrealloc(L->samples, (size_t) (L->alloc * world->sky_nent) * sizeof(float));
      if (world->sky_ng > 0)
        L->gsamples = (float *) myrealloc(L->gsamples, (size_t) (L->alloc * world->sky_ng * 2) * sizeof(float));
    }
  if (L->hist2 == NULL && world->sky_nent > 0)
    L->hist2 = (float *) mycalloc(world->sky_nent * SKY2_NX * SKY2_NY, sizeof(float));
  if (L->ghist2 == NULL && world->sky_ng > 0)
    L->ghist2 = (float *) mycalloc(world->sky_ng * SKY2_NX * SKY2_NG, sizeof(float));
  for (j = 0; j < world->sky_ng; j++)
    {
      const long k = world->sky_gk[j];
      double x = log(world->param0[k]), g = skyline_growth(world, k), glo, gw;
      skyline_grange(world, j, &glo, &gw);
      L->gsamples[(L->n * world->sky_ng + j) * 2] = (float) x;
      L->gsamples[(L->n * world->sky_ng + j) * 2 + 1] = (float) g;
      L->ghist2[(j * SKY2_NX + skyline_bin2(x, skyline_gxlo(world, j), SKY2_XRANGE, SKY2_NX)) * SKY2_NG
                + skyline_bin2(g, glo, gw, SKY2_NG)] += 1.0f;
    }
  for (e = 0; e < world->sky_nent; e++)
    {
      const long k = world->sky_k[e];
      double y = log(world->timek[world->sky_s[e] * world->numpop2 + k]);
      double x = log(world->param0[k]);
      L->samples[L->n * world->sky_nent + e] = (float) y;
      L->hist[e * SKYLINE_HISTBINS + skyline_bin(y)] += 1.0;
      L->hist2[(e * SKY2_NX + skyline_bin2(x, skyline_xlo(world, e), SKY2_XRANGE, SKY2_NX)) * SKY2_NY
               + skyline_bin2(y, -SKY2_YMAX, 2.0 * SKY2_YMAX, SKY2_NY)] += 1.0f;
    }
  L->n++;
  L->nsamples += 1.0;
}

void skyline_param_count(world_fmt *world, boolean accepted)
{
  if ((world->sky_nent == 0 && world->sky_ng == 0) || !world->cold || world->locus >= world->sky_loci)
    return;
  world->skyloc[world->locus].trials += 1.0;
  if (accepted)
    world->skyloc[world->locus].accept += 1.0;
}

/* effective sample size, initial positive sequence estimator (Geyer 1992) */
static double skyline_ess(const float *x, long n, long stride)
{
  double mean = 0.0, var = 0.0, sum = 0.0, rho0, rho1;
  long i, k;
  if (n < 4)
    return (double) n;
  for (i = 0; i < n; i++)
    mean += x[i * stride];
  mean /= n;
  for (i = 0; i < n; i++)
    var += (x[i * stride] - mean) * (x[i * stride] - mean);
  var /= n;
  if (var <= 0.0)
    return (double) n;
  for (k = 0; k + 1 < n; k += 2)
    {
      double c0 = 0.0, c1 = 0.0;
      for (i = 0; i + k < n; i++)
        c0 += (x[i * stride] - mean) * (x[(i + k) * stride] - mean);
      for (i = 0; i + k + 1 < n; i++)
        c1 += (x[i * stride] - mean) * (x[(i + k + 1) * stride] - mean);
      rho0 = c0 / (n * var);
      rho1 = c1 / (n * var);
      if (rho0 + rho1 <= 0.0)
        break;
      sum += rho0 + rho1;
    }
  sum = 2.0 * sum - 1.0;
  return (sum > 1.0) ? n / sum : (double) n;
}

void skyline_locus_finish(world_fmt *world, long locus)
{
  skyloc_fmt *L;
  long e, j, s, i;
  const long nseg = world->timeelements - 1;
  if ((world->sky_nent == 0 && world->sky_ng == 0) || locus >= world->sky_loci)
    return;
  L = &world->skyloc[locus];
  if (L->finished || L->n == 0)
    return;
  for (e = 0; e < world->sky_nent; e++)
    L->ess[e] += skyline_ess(L->samples + e, L->n, world->sky_nent);
  if (world->sky_ng > 0)
    {   /* ESS of the log segment value log Theta - g t */
      float *v = (float *) mycalloc(L->n, sizeof(float));
      for (j = 0; j < world->sky_ng; j++)
        for (s = 1; s < nseg; s++)
          {
            for (i = 0; i < L->n; i++)
              v[i] = L->gsamples[(i * world->sky_ng + j) * 2]
                - L->gsamples[(i * world->sky_ng + j) * 2 + 1] * (float) world->times[s];
            L->gess[j * nseg + s] += skyline_ess(v, L->n, 1);
          }
      myfree(v);
      myfree(L->gsamples);
      L->gsamples = NULL;
    }
  if (L->samples != NULL)
    myfree(L->samples);
  L->samples = NULL;
  L->n = L->alloc = 0;
  L->finished = TRUE;
}

/* quantiles p[0..2] of a distribution given on a regular grid of n cells
   over [lo, lo + width), returned on the natural scale */
static void skyline_grid_quantiles(const double *m, long n, double lo, double width, double *q)
{
  const double p[3] = {0.025, 0.5, 0.975};
  double tot = 0.0, c = 0.0;
  long b, j = 0;
  for (b = 0; b < n; b++)
    tot += m[b];
  q[0] = q[1] = q[2] = 1.0;
  if (tot <= 0.0)
    return;
  for (b = 0; b < n && j < 3; b++)
    {
      while (j < 3 && c + m[b] >= p[j] * tot)
        {
          double f = (m[b] > 0.0) ? (p[j] * tot - c) / m[b] : 0.5;
          q[j] = exp(lo + (b + f) * width / n);
          j++;
        }
      c += m[b];
    }
}

/* distribution of x + y (the log of the segment value) from a joint grid */
static void skyline_value_marginal(const double *P, double xlo, double *m, long nm, double *mlo, double *mwidth)
{
  const double dx = SKY2_XRANGE / SKY2_NX, dy = 2.0 * SKY2_YMAX / SKY2_NY;
  long a, b;
  *mlo = xlo - SKY2_YMAX;
  *mwidth = SKY2_XRANGE + 2.0 * SKY2_YMAX;
  memset(m, 0, (size_t) nm * sizeof(double));
  for (a = 0; a < SKY2_NX; a++)
    for (b = 0; b < SKY2_NY; b++)
      {
        double v = (xlo + (a + 0.5) * dx) + (-SKY2_YMAX + (b + 0.5) * dy);
        m[skyline_bin2(v, *mlo, *mwidth, nm)] += P[a * SKY2_NY + b];
      }
}

/* Gaussian smoothing of a count grid (separable, bandwidths by Silverman's
   rule per axis), normalised; afterwards each axis is contracted about its
   mean by sqrt(var / (var + bandwidth^2)) so that the smoothed density keeps
   the sample variance (otherwise every locus is widened, and dividing by the
   prior loci-1 times turns that into a density that grows towards the grid
   edges) */
static void skyline_smooth2g(const float *H, double *P, long nx, long ny, double dx, double dy)
{
  double n = 0.0, mx = 0.0, my = 0.0, vx = 0.0, vy = 0.0, bx, by, tot = 0.0;
  double *T = (double *) mycalloc(nx * ny, sizeof(double));
  double *kx = (double *) mycalloc(2 * nx + 1, sizeof(double));
  double *ky = (double *) mycalloc(2 * ny + 1, sizeof(double));
  long a, b, c, rx, ry;
  for (a = 0; a < nx; a++)
    for (b = 0; b < ny; b++)
      {
        double h = H[a * ny + b];
        n += h;
        mx += h * a * dx;
        my += h * b * dy;
      }
  mx /= n;
  my /= n;
  for (a = 0; a < nx; a++)
    for (b = 0; b < ny; b++)
      {
        double h = H[a * ny + b];
        vx += h * (a * dx - mx) * (a * dx - mx);
        vy += h * (b * dy - my) * (b * dy - my);
      }
  bx = 1.06 * sqrt(vx / n) * pow(n, -0.2);
  by = 1.06 * sqrt(vy / n) * pow(n, -0.2);
  if (bx < dx) bx = dx;
  if (by < dy) by = dy;
  rx = (long) ceil(4.0 * bx / dx);
  ry = (long) ceil(4.0 * by / dy);
  if (rx > nx) rx = nx;
  if (ry > ny) ry = ny;
  for (c = -rx; c <= rx; c++)
    kx[c + rx] = exp(-0.5 * (c * dx / bx) * (c * dx / bx));
  for (c = -ry; c <= ry; c++)
    ky[c + ry] = exp(-0.5 * (c * dy / by) * (c * dy / by));
  for (a = 0; a < nx; a++)
    for (b = 0; b < ny; b++)
      {
        double sum = 0.0;
        for (c = -rx; c <= rx; c++)
          if (a + c >= 0 && a + c < nx)
            sum += kx[c + rx] * H[(a + c) * ny + b];
        T[a * ny + b] = sum;
      }
  for (a = 0; a < nx; a++)
    for (b = 0; b < ny; b++)
      {
        double sum = 0.0;
        for (c = -ry; c <= ry; c++)
          if (b + c >= 0 && b + c < ny)
            sum += ky[c + ry] * T[a * ny + b + c];
        P[a * ny + b] = sum;
        tot += sum;
      }
  for (a = 0; a < nx * ny; a++)
    P[a] /= tot;
  {   /* variance correction: P'(u) = P(m + (u - m) / c) / c per axis */
    const double cx = (vx > 0.0) ? sqrt((vx / n) / (vx / n + bx * bx)) : 1.0;
    const double cy = (vy > 0.0) ? sqrt((vy / n) / (vy / n + by * by)) : 1.0;
    tot = 0.0;
    for (a = 0; a < nx; a++)
      for (b = 0; b < ny; b++)
        {
          double u = (mx + (a * dx - mx) / cx) / dx, w = (my + (b * dy - my) / cy) / dy, fu, fw, v = 0.0;
          long i0 = (long) floor(u), j0 = (long) floor(w);
          fu = u - i0;
          fw = w - j0;
          if (i0 >= 0 && i0 + 1 < nx && j0 >= 0 && j0 + 1 < ny)
            v = (1 - fu) * (1 - fw) * P[i0 * ny + j0] + fu * (1 - fw) * P[(i0 + 1) * ny + j0]
              + (1 - fu) * fw * P[i0 * ny + j0 + 1] + fu * fw * P[(i0 + 1) * ny + j0 + 1];
          T[a * ny + b] = v;
          tot += v;
        }
    if (tot > 0.0)
      for (a = 0; a < nx * ny; a++)
        P[a] = T[a] / tot;
  }
  myfree(T);
  myfree(kx);
  myfree(ky);
}

static void skyline_smooth2(const float *H, double *P)
{
  skyline_smooth2g(H, P, SKY2_NX, SKY2_NY, SKY2_XRANGE / SKY2_NX, 2.0 * SKY2_YMAX / SKY2_NY);
}

/* distribution of cx x + cy y on a grid of nx x ny cells over
   [xlo, xlo + xw) x [ylo, ylo + yw) */
static void skyline_lin_marginal(const double *P, long nx, double xlo, double xw, long ny, double ylo, double yw,
                                 double cx, double cy, double *m, long nm, double *mlo, double *mwidth)
{
  const double lo = fmin(cx * xlo, cx * (xlo + xw)) + fmin(cy * ylo, cy * (ylo + yw));
  const double hi = fmax(cx * xlo, cx * (xlo + xw)) + fmax(cy * ylo, cy * (ylo + yw));
  long a, b;
  *mlo = lo;
  *mwidth = (hi > lo) ? hi - lo : 1e-12;
  memset(m, 0, (size_t) nm * sizeof(double));
  for (a = 0; a < nx; a++)
    for (b = 0; b < ny; b++)
      {
        double v = cx * (xlo + (a + 0.5) * xw / nx) + cy * (ylo + (b + 0.5) * yw / ny);
        m[skyline_bin2(v, *mlo, *mwidth, nm)] += P[a * ny + b];
      }
}

/* growth-linked Theta: quantiles of the multiplier exp(-g t) and of the value
   Theta exp(-g t) at the start t of every segment s >= 1, 6 per segment */
static void skyline_g_quantiles(world_fmt *world, long j, const double *P, double *q)
{
  const long nseg = world->timeelements - 1, nm = 400;
  const double xlo = skyline_gxlo(world, j);
  double *m = (double *) mycalloc(nm, sizeof(double));
  double glo, gw, mlo, mw;
  long sg;
  skyline_grange(world, j, &glo, &gw);
  for (sg = 1; sg < nseg; sg++)
    {
      const double t = world->times[sg];
      skyline_lin_marginal(P, SKY2_NX, xlo, SKY2_XRANGE, SKY2_NG, glo, gw, 0.0, -t, m, nm, &mlo, &mw);
      skyline_grid_quantiles(m, nm, mlo, mw, q + 6 * sg);
      skyline_lin_marginal(P, SKY2_NX, xlo, SKY2_XRANGE, SKY2_NG, glo, gw, 1.0, -t, m, nm, &mlo, &mw);
      skyline_grid_quantiles(m, nm, mlo, mw, q + 6 * sg + 3);
    }
  myfree(m);
}

/* add a complete locus to the running sums (and record its value quantiles) */
static void skyline_joint_add(world_fmt *world, long locus)
{
  skyloc_fmt *L;
  double *P, *m, mlo, mwidth;
  long e, j;
  const long nm = 400;
  if ((world->sky_nent == 0 && world->sky_ng == 0) || locus < 0 || locus >= world->sky_loci)
    return;
  L = &world->skyloc[locus];
  if (L->combined || (L->hist2 == NULL && L->ghist2 == NULL))
    return;
  skyline_locus_finish(world, locus);
  if (L->ghist2 != NULL)
    {
      const long nseg = world->timeelements - 1;
      double *PG = (double *) mycalloc(SKY2_NX * SKY2_NG, sizeof(double));
      long g;
      for (g = 0; g < world->sky_ng; g++)
        {
          const float *H = L->ghist2 + g * SKY2_NX * SKY2_NG;
          double n = 0.0, glo, gw;
          for (j = 0; j < SKY2_NX * SKY2_NG; j++)
            n += H[j];
          if (n <= 0.0)
            continue;
          skyline_grange(world, g, &glo, &gw);
          skyline_smooth2g(H, PG, SKY2_NX, SKY2_NG, SKY2_XRANGE / SKY2_NX, gw / SKY2_NG);
          for (j = 0; j < SKY2_NX * SKY2_NG; j++)
            world->sky2gsum[g * SKY2_NX * SKY2_NG + j] += log(PG[j] + 1e-12);
          world->sky2gloci[g]++;
          skyline_g_quantiles(world, g, PG, L->gq + 6 * g * nseg);
        }
      myfree(PG);
      myfree(L->ghist2);
      L->ghist2 = NULL;
    }
  if (L->hist2 == NULL)
    {
      L->combined = TRUE;
      return;
    }
  P = (double *) mycalloc(SKY2_NX * SKY2_NY, sizeof(double));
  m = (double *) mycalloc(nm, sizeof(double));
  for (e = 0; e < world->sky_nent; e++)
    {
      const float *H = L->hist2 + e * SKY2_NX * SKY2_NY;
      double n = 0.0;
      for (j = 0; j < SKY2_NX * SKY2_NY; j++)
        n += H[j];
      if (n <= 0.0)
        continue;
      skyline_smooth2(H, P);
      for (j = 0; j < SKY2_NX * SKY2_NY; j++)
        world->sky2sum[e * SKY2_NX * SKY2_NY + j] += log(P[j] + 1e-12);
      world->sky2loci[e]++;
      skyline_value_marginal(P, skyline_xlo(world, e), m, nm, &mlo, &mwidth);
      skyline_grid_quantiles(m, nm, mlo, mwidth, L->vq + 3 * e);
    }
  myfree(P);
  myfree(m);
  myfree(L->hist2);
  L->hist2 = NULL;
  L->combined = TRUE;
}

/* MPI: appended to the skyline buffer of a locus (see pack_skyline_buffer) */
long skyline_param_bufsize(world_fmt *world)
{
  return 2 + world->sky_nent * (SKYLINE_HISTBINS + 1 + SKY2_NX * SKY2_NY)
    + world->sky_ng * ((world->timeelements - 1) + SKY2_NX * SKY2_NG) + 3;
}

long skyline_param_pack(world_fmt *world, long locus, MYREAL *buf)
{
  long z = 0, e, b;
  const long nseg = world->timeelements - 1;
  skyloc_fmt *L;
  buf[z++] = (MYREAL) world->sky_nent;
  buf[z++] = (MYREAL) world->sky_ng;
  if ((world->sky_nent == 0 && world->sky_ng == 0) || locus >= world->sky_loci)
    {
      for (e = 2; e < skyline_param_bufsize(world); e++)
        buf[z++] = 0.0;
      return z;
    }
  skyline_locus_finish(world, locus);
  L = &world->skyloc[locus];
  for (e = 0; e < world->sky_nent; e++)
    {
      for (b = 0; b < SKYLINE_HISTBINS; b++)
        buf[z++] = (MYREAL) L->hist[e * SKYLINE_HISTBINS + b];
      buf[z++] = (MYREAL) L->ess[e];
      for (b = 0; b < SKY2_NX * SKY2_NY; b++)
        buf[z++] = (L->hist2 != NULL) ? (MYREAL) L->hist2[e * SKY2_NX * SKY2_NY + b] : 0.0;
    }
  for (e = 0; e < world->sky_ng; e++)
    {
      for (b = 0; b < nseg; b++)
        buf[z++] = (MYREAL) L->gess[e * nseg + b];
      for (b = 0; b < SKY2_NX * SKY2_NG; b++)
        buf[z++] = (L->ghist2 != NULL) ? (MYREAL) L->ghist2[e * SKY2_NX * SKY2_NG + b] : 0.0;
    }
  buf[z++] = (MYREAL) L->nsamples;
  buf[z++] = (MYREAL) L->trials;
  buf[z++] = (MYREAL) L->accept;
  return z;
}

long skyline_param_unpack(world_fmt *world, long locus, MYREAL *buf)
{
  long z = 0, e, b;
  const long nseg = world->timeelements - 1;
  skyloc_fmt *L;
  long nent = (long) buf[z++];
  long ng = (long) buf[z++];
  if (nent != world->sky_nent || ng != world->sky_ng)
    error("skyline=PARAM: worker and master disagree on the number of multipliers");
  if ((nent == 0 && ng == 0) || locus >= world->sky_loci)
    return skyline_param_bufsize(world);
  L = &world->skyloc[locus];
  if (L->hist2 == NULL && !L->combined && nent > 0)
    L->hist2 = (float *) mycalloc(nent * SKY2_NX * SKY2_NY, sizeof(float));
  if (L->ghist2 == NULL && !L->combined && ng > 0)
    L->ghist2 = (float *) mycalloc(ng * SKY2_NX * SKY2_NG, sizeof(float));
  for (e = 0; e < nent; e++)
    {
      for (b = 0; b < SKYLINE_HISTBINS; b++)
        L->hist[e * SKYLINE_HISTBINS + b] += buf[z++];
      L->ess[e] += buf[z++];
      for (b = 0; b < SKY2_NX * SKY2_NY; b++)
        {
          if (L->hist2 != NULL)
            L->hist2[e * SKY2_NX * SKY2_NY + b] += (float) buf[z];
          z++;
        }
    }
  for (e = 0; e < ng; e++)
    {
      for (b = 0; b < nseg; b++)
        L->gess[e * nseg + b] += buf[z++];
      for (b = 0; b < SKY2_NX * SKY2_NG; b++)
        {
          if (L->ghist2 != NULL)
            L->ghist2[e * SKY2_NX * SKY2_NG + b] += (float) buf[z];
          z++;
        }
    }
  L->nsamples += buf[z++];
  L->trials += buf[z++];
  L->accept += buf[z++];
  L->finished = TRUE;
  /* a locus arrives in one piece unless its replicates ran on several
     workers; then it is combined in skyline_param_summary() */
  if (!(world->options->replicate && world->options->replicatenum > 1))
    skyline_joint_add(world, locus);
  return z;
}

/* quantiles of a histogram of log tau, returned as tau */
static void skyline_quantiles(const double *h, double *q025, double *q50, double *q975)
{
  double q[3];
  skyline_grid_quantiles(h, SKYLINE_HISTBINS, -SKYLINE_HISTMAX, 2.0 * SKYLINE_HISTMAX, q);
  *q025 = q[0];
  *q50 = q[1];
  *q975 = q[2];
}

static void skyline_label(world_fmt *world, long k, char *label, size_t size)
{
  long from, to;
  if (k < world->numpop)
    snprintf(label, size, "Theta_%li", k + 1);
  else
    {
      m2mm(k, world->numpop, &from, &to);
      if (world->options->usem)
        snprintf(label, size, "M_%li->%li", from + 1, to + 1);
      else
        snprintf(label, size, "Theta_%li*M_%li->%li", to + 1, from + 1, to + 1);
    }
}

/* growth-linked Theta_p: the loci combined on the (log Theta, g) grid,
   product of the per-locus densities divided by the prior (loci-1) times;
   returns the number of loci used */
static long skyline_growth_combined(world_fmt *world, long j, double *P)
{
  const long used = world->sky2gloci[j];
  const long k = world->sky_gk[j];
  const long pi = skyline_prior_index(world, k);
  const long gi = skyline_growth_index(world, k);
  const double dx = SKY2_XRANGE / SKY2_NX, xlo = skyline_gxlo(world, j);
  double glo, gw, mxv = -HUGE_VAL, lpmax = -HUGE_VAL;
  double *lpgrid = (double *) mycalloc(SKY2_NX * SKY2_NG, sizeof(double));
  long a, b;
  if (used == 0)
    {
      myfree(lpgrid);
      return 0;
    }
  skyline_grange(world, j, &glo, &gw);
  for (a = 0; a < SKY2_NX; a++)
    {
      double x = xlo + (a + 0.5) * dx;
      double lpx = (pi >= 0) ? (*log_prior_1[pi])(world, pi, exp(x)) : 0.0;
      lpx = (lpx <= (double) -HUGE / 2.0) ? -HUGE_VAL : lpx + x;
      for (b = 0; b < SKY2_NG; b++)
        {
          double g = glo + (b + 0.5) * gw / SKY2_NG;
          double lpy = (*log_prior_1[gi])(world, gi, g), v;
          lpy = (lpy <= (double) -HUGE / 2.0) ? -HUGE_VAL : lpy;
          lpgrid[a * SKY2_NG + b] = lpx + lpy;
          v = world->sky2gsum[j * SKY2_NX * SKY2_NG + a * SKY2_NG + b];
          if (used > 1)
            v -= (used - 1) * (lpx + lpy);
          if (!(lpx + lpy > -HUGE_VAL) || isnan(v))
            v = -HUGE_VAL;
          P[a * SKY2_NG + b] = v;
          if (v > mxv)
            mxv = v;
          if (lpgrid[a * SKY2_NG + b] > lpmax)
            lpmax = lpgrid[a * SKY2_NG + b];
        }
    }
  for (a = 0; a < SKY2_NX * SKY2_NG; a++)
    {
      P[a] = (mxv > -HUGE_VAL && P[a] > -HUGE_VAL) ? exp(P[a] - mxv) : 0.0;
      if (used > 1 && !(lpgrid[a] > lpmax - SKY2_PRIORSPAN))
        P[a] = 0.0;   /* as for the multipliers: only where the prior has mass */
    }
  myfree(lpgrid);
  return used;
}

/* rows of a growth-linked Theta_p: in every segment the multiplier
   exp(-g t) and the value Theta exp(-g t) at the segment start t */
static long skyline_growth_rows(world_fmt *world, long k, skyline_row_fmt *rows, boolean perlocus)
{
  const long nseg = world->timeelements - 1;
  const long loci = world->sky_loci;
  double *P = (double *) mycalloc(SKY2_NX * SKY2_NG, sizeof(double));
  double *q = (double *) mycalloc(6 * nseg, sizeof(double));
  long j, s, l, n = 0, used;
  for (j = 0; j < world->sky_ng; j++)
    if (world->sky_gk[j] == k)
      break;
  used = skyline_growth_combined(world, j, P);
  if (used > 0)
    skyline_g_quantiles(world, j, P, q);
  for (s = 0; s < nseg; s++)
    {
      skyline_row_fmt *r = &rows[n++];
      snprintf(r->label, sizeof(r->label), "Theta_%li (growth)", k + 1);
      r->seg = s;
      r->t0 = world->times[s];
      r->t1 = (s + 1 < nseg) ? world->times[s + 1] : -1.0;
      r->locus = -1;
      r->v025 = r->v50 = r->v975 = -1.0;
      if (s == 0)
        {   /* anchor: today */
          r->q025 = r->q50 = r->q975 = 1.0;
          r->ess = -1.0;
          continue;
        }
      r->ess = 0.0;
      for (l = 0; l < loci; l++)
        r->ess += world->skyloc[l].gess[j * nseg + s];
      if (used > 0)
        {
          r->q025 = q[6 * s]; r->q50 = q[6 * s + 1]; r->q975 = q[6 * s + 2];
          r->v025 = q[6 * s + 3]; r->v50 = q[6 * s + 4]; r->v975 = q[6 * s + 5];
        }
      else
        r->q025 = r->q50 = r->q975 = 1.0;
      if (perlocus)
        for (l = 0; l < loci; l++)
          {
            const double *lq = world->skyloc[l].gq + 6 * (j * nseg + s);
            skyline_row_fmt *rl = &rows[n++];
            *rl = *r;
            rl->locus = l;
            rl->ess = world->skyloc[l].gess[j * nseg + s];
            rl->q025 = lq[0]; rl->q50 = lq[1]; rl->q975 = lq[2];
            rl->v025 = lq[3]; rl->v50 = lq[4]; rl->v975 = lq[5];
          }
    }
  myfree(P);
  myfree(q);
  return n;
}

/* rows of the skyline table: for every free parameter and segment the
   multiplier and the segment value (parameter times multiplier) over all
   loci from the joint grid, and, with at most SKYLINE_MAXLOCUSROWS loci,
   each locus */
#define SKYLINE_MAXLOCUSROWS 10
long skyline_param_summary(world_fmt *world, skyline_row_fmt **rows, double *trials, double *accept)
{
  const long nseg = world->timeelements - 1;
  const long loci = world->sky_loci;
  const double dx = SKY2_XRANGE / SKY2_NX, dy = 2.0 * SKY2_YMAX / SKY2_NY;
  const long nm = 400;
  long e, l, a, b, n = 0, k, s;
  double qjt[3], qjv[3];
  long perlocus = (loci > 1 && loci <= SKYLINE_MAXLOCUSROWS);
  double *P = (double *) mycalloc(SKY2_NX * SKY2_NY, sizeof(double));
  double *mt = (double *) mycalloc(SKY2_NY, sizeof(double));
  double *mv = (double *) mycalloc(nm, sizeof(double));
  double *lpgrid = (double *) mycalloc(SKY2_NX * SKY2_NY, sizeof(double));
  *trials = *accept = 0.0;
  *rows = (skyline_row_fmt *) mycalloc(world->numpop2 * nseg * (loci + 1) + 1, sizeof(skyline_row_fmt));
  for (l = 0; l < loci; l++)
    {
      skyline_joint_add(world, l);
      *trials += world->skyloc[l].trials;
      *accept += world->skyloc[l].accept;
    }
  for (k = 0; k < world->numpop2; k++)
    {
      if (world->skymodel2[k] == 'g')
        {
          n += skyline_growth_rows(world, k, *rows + n, perlocus);
          continue;
        }
      if (world->skymodel2[k] != '*')
        continue;
      for (s = 0; s < nseg; s++)
        {
          skyline_row_fmt *r = &(*rows)[n++];
          skyline_label(world, k, r->label, sizeof(r->label));
          r->seg = s;
          r->t0 = world->times[s];
          r->t1 = (s + 1 < nseg) ? world->times[s + 1] : -1.0;
          r->locus = -1;
          r->v025 = r->v50 = r->v975 = -1.0;
          if (s == 0)
            {   /* anchor: today */
              r->q025 = r->q50 = r->q975 = 1.0;
              r->ess = -1.0;
              continue;
            }
          for (e = 0; e < world->sky_nent; e++)
            if (world->sky_k[e] == k && world->sky_s[e] == s)
              break;
          r->ess = 0.0;
          for (l = 0; l < loci; l++)
            r->ess += world->skyloc[l].ess[e];
          if (jc_sky_quantiles(world, k, s, qjt, qjv))
            {   /* the joint multi-locus combination (joint_combine.c) */
              r->q025 = qjt[0]; r->q50 = qjt[1]; r->q975 = qjt[2];
              r->v025 = qjv[0]; r->v50 = qjv[1]; r->v975 = qjv[2];
            }
          else if (world->sky2loci[e] > 0)
            {
              const long used = world->sky2loci[e];
              const long pi = skyline_prior_index(world, k);
              const double xlo = skyline_xlo(world, e);
              double mxv = -HUGE_VAL, qv[3], qt[3];
              for (a = 0; a < SKY2_NX; a++)
                {
                  double x = xlo + (a + 0.5) * dx;
                  double lpx = (pi >= 0) ? (*log_prior_1[pi])(world, pi, exp(x)) : 0.0;
                  /* the prior functions return -HUGE (finite) outside their bounds */
                  lpx = (lpx <= (double) -HUGE / 2.0) ? -HUGE_VAL : lpx + x;
                  for (b = 0; b < SKY2_NY; b++)
                    {
                      double y = -SKY2_YMAX + (b + 0.5) * dy, lpy = 0.0, v;
                      if (world->skyptype[k] == SKYPRIOR_RANDOMWALK)
                        lpy = -0.5 * y * y / (s * world->skypa[k] * world->skypa[k]);
                      else if (y < log(world->skypa[k]) || y > log(world->skypb[k]))
                        lpy = -HUGE_VAL;
                      lpgrid[a * SKY2_NY + b] = lpx + lpy;
                      /* the segment value (parameter x multiplier) must lie in
                         the parameter's prior range, as in the sampler */
                      if (pi >= 0 && (x + y > log(world->bayes->maxparam[pi])
                                      || (world->bayes->minparam[pi] > 0.0
                                          && x + y < log(world->bayes->minparam[pi]))))
                        lpy = -HUGE_VAL;
                      v = world->sky2sum[e * SKY2_NX * SKY2_NY + a * SKY2_NY + b];
                      if (used > 1)
                        v -= (used - 1) * (lpx + lpy);
                      if (!(lpx + lpy > -HUGE_VAL) || isnan(v))
                        v = -HUGE_VAL;
                      P[a * SKY2_NY + b] = v;
                      if (v > mxv)
                        mxv = v;
                    }
                }
              for (a = 0; a < SKY2_NX * SKY2_NY; a++)
                P[a] = (mxv > -HUGE_VAL && P[a] > -HUGE_VAL) ? exp(P[a] - mxv) : 0.0;
              if (used > 1)
                {   /* only where the prior has mass: in cells the loci never
                       visited their densities are floored, and dividing by a
                       tiny prior loci-1 times would inflate them */
                  double lpmax = -HUGE_VAL;
                  for (a = 0; a < SKY2_NX * SKY2_NY; a++)
                    if (lpgrid[a] > lpmax)
                      lpmax = lpgrid[a];
                  for (a = 0; a < SKY2_NX * SKY2_NY; a++)
                    if (!(lpgrid[a] > lpmax - SKY2_PRIORSPAN))
                      P[a] = 0.0;
                }
              memset(mt, 0, SKY2_NY * sizeof(double));
              for (a = 0; a < SKY2_NX; a++)
                for (b = 0; b < SKY2_NY; b++)
                  mt[b] += P[a * SKY2_NY + b];
              skyline_grid_quantiles(mt, SKY2_NY, -SKY2_YMAX, 2.0 * SKY2_YMAX, qt);
              {
                double mlo, mwidth;
                skyline_value_marginal(P, xlo, mv, nm, &mlo, &mwidth);
                skyline_grid_quantiles(mv, nm, mlo, mwidth, qv);
              }
              r->q025 = qt[0]; r->q50 = qt[1]; r->q975 = qt[2];
              r->v025 = qv[0]; r->v50 = qv[1]; r->v975 = qv[2];
            }
          else
            r->q025 = r->q50 = r->q975 = 1.0;
          if (perlocus)
            for (l = 0; l < loci; l++)
              {
                skyline_row_fmt *rl = &(*rows)[n++];
                *rl = *r;
                rl->locus = l;
                rl->ess = world->skyloc[l].ess[e];
                skyline_quantiles(world->skyloc[l].hist + e * SKYLINE_HISTBINS,
                                  &rl->q025, &rl->q50, &rl->q975);
                rl->v025 = world->skyloc[l].vq[3 * e];
                rl->v50 = world->skyloc[l].vq[3 * e + 1];
                rl->v975 = world->skyloc[l].vq[3 * e + 2];
              }
        }
    }
  myfree(P);
  myfree(mt);
  myfree(mv);
  myfree(lpgrid);
  return n;
}

static void skyline_cells(const skyline_row_fmt *r, char c[12][64])
{
  snprintf(c[0], 64, "%s", r->label);
  snprintf(c[1], 64, "%li", r->seg + 1);
  snprintf(c[2], 64, "%.5f", r->t0);
  if (r->t1 < 0.0) snprintf(c[3], 64, "inf"); else snprintf(c[3], 64, "%.5f", r->t1);
  if (r->locus < 0) snprintf(c[4], 64, "All"); else snprintf(c[4], 64, "%li", r->locus + 1);
  snprintf(c[5], 64, "%.4f", r->q025);
  snprintf(c[6], 64, "%.4f", r->q50);
  snprintf(c[7], 64, "%.4f", r->q975);
  if (r->v50 < 0.0)
    {
      snprintf(c[8], 64, "-"); snprintf(c[9], 64, "-"); snprintf(c[10], 64, "-");
    }
  else
    {
      snprintf(c[8], 64, "%.5g", r->v025);
      snprintf(c[9], 64, "%.5g", r->v50);
      snprintf(c[10], 64, "%.5g", r->v975);
    }
  if (r->ess < 0.0) snprintf(c[11], 64, "-"); else snprintf(c[11], 64, "%.0f", r->ess);
}

void skyline_param_print(FILE *file, world_fmt *world)
{
  skyline_row_fmt *rows;
  double trials, accept;
  long n, i;
  if (world->sky_nent == 0 && world->sky_ng == 0)
    return;
  n = skyline_param_summary(world, &rows, &trials, &accept);
  fprintf(file, "\n\nParametric skyline: multipliers of Theta and M per time segment\n");
  fprintf(file, "================================================================\n");
  fprintf(file, "Segment 1 (today) is the reference (multiplier 1); the value in a segment\n");
  fprintf(file, "is the parameter times its multiplier. All loci: %s\n\n",
          world->bayes->jointrow
          ? "from the joint multi-locus combination\n(the same samples as the Bayesian estimates of the parameters); a parameter or\nmultiplier flagged Joint* (large genealogy-sampling error) uses the product of\nthe per-locus posteriors instead."
          : "the per-locus joint\nposteriors of parameter and multiplier are multiplied and divided by the\nprior (loci-1) times, as for the other parameters.");
  if (world->sky_ng > 0)
    fprintf(file, "Theta (growth): a growth-linked Theta (skyline-model g) has no free multipliers;\nits multiplier is exp(-g t) at the segment start t, from the joint posterior of\nTheta and its growth rate g.\n\n");
  fprintf(file, "                                                           Multiplier                   Value in segment\n");
  fprintf(file, "Parameter            Seg  From       To         Locus     2.5%%   Median    97.5%%       2.5%%     Median      97.5%%     ESS\n");
  fprintf(file, "-----------------------------------------------------------------------------------------------------------------------\n");
  for (i = 0; i < n; i++)
    {
      char c[12][64];
      skyline_cells(&rows[i], c);
      fprintf(file, "%-20s %3s  %-10s %-10s %-5s %8s %8s %8s  %10s %10s %10s %7s\n",
              c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11]);
    }
  if (trials > 0.0)
    fprintf(file, "\nAcceptance ratio of the multiplier moves: %.0f/%.0f = %.5f\n", accept, trials, accept / trials);
  fflush(file);
  myfree(rows);
}
