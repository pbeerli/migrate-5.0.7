/*------------------------------------------------------------------------
  Mittag-Leffler alpha in the joint multi-locus combination (2026-11-01).

  Waiting-time functions. A population's clock with rate lambda and index
  alpha, with u = lambda^(1/alpha) t, survives with S(u) = E_alpha(-u^alpha)
  and has the density f(u) = u^(alpha-1) E_{alpha,alpha}(-u^alpha) in u.
  Both are Laplace transforms of the spectral density of the Mittag-Leffler
  distribution; with r = e^s
    S(u) = int exp(-u e^s) g(s) ds,   f(u) = int exp(s - u e^s) g(s) ds,
    g(s) = sin(alpha pi) / (pi (4 sinh^2(alpha s / 2) + 4 cos^2(alpha pi / 2))).
  The table holds log S and log f on a grid of log u for rows in
  b = log(1 - alpha) (the alpha -> 1 limit is not uniform in u: a tail of
  weight about (1 - alpha) u^-alpha takes over from exp(-u); in b it is
  smooth), with cubic splines in log u and four-point Lagrange
  interpolation in b. Above JML_AONE alpha is 1 (exponential).

  Spells and the normalizer: docs/joint_mlalpha_design.md, sections 6-7.
  A population's spells depend on its free immigration rates M_j only
  through M_tot (in lambda) and log M_j at a migration to j; with uniform
  priors U(0, mmax) the integral over the M_j reduces to M_tot with the
  factor D(M_tot) = int_simplex prod w_j^m_j 1{w_j <= mmax/M_tot} dw (closed
  form for one and two routes). The inner integral over (Theta, M_tot) at
  fixed alpha is an adaptive Gauss-Hermite rule around the mode in logit
  coordinates; the caller integrates alpha with Gauss-Legendre nodes.
------------------------------------------------------------------------*/
#include "migration.h"
#include "tools.h"
#include "sighandler.h"
#include "joint_ml.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define JML_LU0 (-25.0)
#define JML_LUH 0.02
#define JML_NLU 2501
#define JML_BH 0.05
#define JML_AONE (1.0 - 1e-6)
#define JML_GH 9          /* Gauss-Hermite nodes per dimension */

static struct
{
  long nb;
  double b0;              /* b of row 0; rows increase in b */
  double *S, *S2, *F, *F2;
  double amin, amax;
} T;

/* ------------------------------------------------------- quadrature */

typedef struct
{
  double u, e;            /* u, 1 - alpha */
  double a, sa, ca;       /* alpha, sin(alpha pi), cos(alpha pi / 2) */
  int dens;               /* 0: S, 1: f */
} jml_fx;

static double
jml_g (const jml_fx *p, double s)
{
  const double sh = sinh (0.5 * p->a * s);
  return p->sa / (PI * (4.0 * sh * sh + 4.0 * p->ca * p->ca));
}

static double
jml_integrand (const jml_fx *p, double s)
{
  const double es = exp (s);
  const double x = p->u * es;
  if (x > 745.0)
    return 0.0;
  return (p->dens ? es : 1.0) * exp (-x) * jml_g (p, s);
}

static const double gk_x[8] = {0.991455371120812639, 0.949107912342758525, 0.864864423359769073,
                               0.741531185599394440, 0.586087235467691130, 0.405845151377397167,
                               0.207784955007898468, 0.0};
static const double gk_wk[8] = {0.022935322010529225, 0.063092092629978553, 0.104790010322250184,
                                0.140653259715525919, 0.169004726639267903, 0.190350578064785410,
                                0.204432940075298892, 0.209482141084727828};
static const double gk_wg[4] = {0.129484966168869693, 0.279705391489276668, 0.381830050505118945,
                                0.417959183673469388};

static double
jml_gk15 (const jml_fx *p, double a, double b, double *err)
{
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double fc = jml_integrand (p, c);
  double rk = fc * gk_wk[7], rg = fc * gk_wg[3];
  long j;
  for (j = 0; j < 7; j++)
    {
      const double d = h * gk_x[j];
      const double f = jml_integrand (p, c - d) + jml_integrand (p, c + d);
      rk += gk_wk[j] * f;
      if (j % 2 == 1)
        rg += gk_wg[j / 2] * f;
    }
  *err = fabs ((rk - rg) * h);
  return rk * h;
}

static double
jml_adapt (const jml_fx *p, double a, double b, double tol, long depth)
{
  double err;
  const double v = jml_gk15 (p, a, b, &err);
  if (err <= tol || depth >= 40)
    return v;
  const double m = 0.5 * (a + b);
  return jml_adapt (p, a, m, 0.5 * tol, depth + 1) + jml_adapt (p, m, b, 0.5 * tol, depth + 1);
}

static int
jml_cmp (const void *x, const void *y)
{
  const double a = *(const double *) x, b = *(const double *) y;
  return a < b ? -1 : (a > b ? 1 : 0);
}

/* log S (dens 0) or log f (dens 1) at log u for 1 - alpha = e */
static double
jml_direct (double lu, double e, int dens)
{
  jml_fx p;
  double pts[12], lo, hi, w, v = 0.0, est = 0.0, err;
  long n = 0, i;
  const double kk[8] = {-30.0, -10.0, -3.0, -1.0, 1.0, 3.0, 10.0, 30.0};
  p.u = exp (lu);
  p.e = e;
  p.a = 1.0 - e;
  p.sa = sin (PI * e);
  p.ca = sin (0.5 * PI * e);
  p.dens = dens;
  lo = -45.0 / p.a;
  hi = -lu + log (60.0);
  if (hi < 5.0)
    hi = 5.0;
  w = PI * e / p.a;   /* width of the peak of g at s = 0 */
  pts[n++] = lo;
  if (0.0 > lo && 0.0 < hi)
    pts[n++] = 0.0;
  if (-lu > lo && -lu < hi)
    pts[n++] = -lu;
  for (i = 0; i < 8; i++)
    if (kk[i] * w > lo && kk[i] * w < hi)
      pts[n++] = kk[i] * w;
  pts[n++] = hi;
  qsort (pts, (size_t) n, sizeof (double), jml_cmp);
  for (i = 0; i + 1 < n; i++)
    est += jml_gk15 (&p, pts[i], pts[i + 1], &err);
  for (i = 0; i + 1 < n; i++)
    if (pts[i + 1] > pts[i])
      v += jml_adapt (&p, pts[i], pts[i + 1], 1e-12 * fabs (est) / (double) (n - 1), 0);
  return log (v);
}

/* natural cubic spline second derivatives of y on the uniform log u grid */
static void
jml_spline (const double *y, double *y2)
{
  double *u = (double *) mycalloc (JML_NLU, sizeof (double));
  long i;
  y2[0] = u[0] = 0.0;
  for (i = 1; i < JML_NLU - 1; i++)
    {
      const double p = 0.5 * y2[i - 1] + 2.0;
      y2[i] = -0.5 / p;
      u[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / JML_LUH;
      u[i] = (3.0 * u[i] / JML_LUH - 0.5 * u[i - 1]) / p;
    }
  y2[JML_NLU - 1] = 0.0;
  for (i = JML_NLU - 2; i >= 0; i--)
    y2[i] = y2[i] * y2[i + 1] + u[i];
  myfree (u);
}

/// builds the table for alpha in [amin, amax] (once per run)
void
jml_table_init (double amin, double amax)
{
  long r, j;
  double bhi, blo;
  if (amax > JML_AONE)   /* clamp first: the check below compares clamped values */
    amax = JML_AONE;
  if (amin < 0.05)
    amin = 0.05;
  if (T.S != NULL && amin >= T.amin && amax <= T.amax)
    return;
  myfree (T.S);
  myfree (T.S2);
  myfree (T.F);
  myfree (T.F2);
  blo = log (1.0 - amax) - 2.0 * JML_BH;   /* stencil margins */
  bhi = log (1.0 - amin) + 2.0 * JML_BH;
  if (bhi > log (0.97))
    bhi = log (0.97);
  T.nb = (long) ceil ((bhi - blo) / JML_BH) + 1;
  T.b0 = blo;
  T.S = (double *) mycalloc ((size_t) (T.nb * JML_NLU), sizeof (double));
  T.S2 = (double *) mycalloc ((size_t) (T.nb * JML_NLU), sizeof (double));
  T.F = (double *) mycalloc ((size_t) (T.nb * JML_NLU), sizeof (double));
  T.F2 = (double *) mycalloc ((size_t) (T.nb * JML_NLU), sizeof (double));
  for (r = 0; r < T.nb; r++)
    {
      const double e = exp (T.b0 + JML_BH * (double) r);
      for (j = 0; j < JML_NLU; j++)
        {
          const double lu = JML_LU0 + JML_LUH * (double) j;
          T.S[r * JML_NLU + j] = jml_direct (lu, e, 0);
          T.F[r * JML_NLU + j] = jml_direct (lu, e, 1);
        }
      jml_spline (T.S + r * JML_NLU, T.S2 + r * JML_NLU);
      jml_spline (T.F + r * JML_NLU, T.F2 + r * JML_NLU);
    }
  T.amin = amin;
  T.amax = amax;
}

/// the table at one alpha: Lagrange combination of four rows
void
jml_alpha_set (jml_alpha *A, double alpha)
{
  long r0, j, q;
  double w[4];
  A->alpha = alpha;
  A->one = alpha >= JML_AONE;
  if (A->one)
    return;
  if (A->S == NULL)
    {
      A->S = (double *) mycalloc (JML_NLU, sizeof (double));
      A->S2 = (double *) mycalloc (JML_NLU, sizeof (double));
      A->F = (double *) mycalloc (JML_NLU, sizeof (double));
      A->F2 = (double *) mycalloc (JML_NLU, sizeof (double));
    }
  const double x = (log (1.0 - alpha) - T.b0) / JML_BH;
  r0 = (long) floor (x) - 1;
  if (r0 < 0)
    r0 = 0;
  if (r0 > T.nb - 4)
    r0 = T.nb - 4;
  for (q = 0; q < 4; q++)
    {
      double v = 1.0;
      long m;
      for (m = 0; m < 4; m++)
        if (m != q)
          v *= (x - (double) (r0 + m)) / (double) (q - m);
      w[q] = v;
    }
  for (j = 0; j < JML_NLU; j++)
    {
      A->S[j] = A->S2[j] = A->F[j] = A->F2[j] = 0.0;
      for (q = 0; q < 4; q++)
        {
          const long i = (r0 + q) * JML_NLU + j;
          A->S[j] += w[q] * T.S[i];
          A->S2[j] += w[q] * T.S2[i];
          A->F[j] += w[q] * T.F[i];
          A->F2[j] += w[q] * T.F2[i];
        }
    }
  /* the spline as a cubic in t per interval:
     y = y_j + t (y_j+1 - y_j) + h^2/6 [y2_j (-2t + 3t^2 - t^3) + y2_j+1 (t^3 - t)] */
  if (A->cS == NULL)
    {
      A->cS = (double *) mycalloc (4 * JML_NLU, sizeof (double));
      A->cF = (double *) mycalloc (4 * JML_NLU, sizeof (double));
    }
  {
    const double h6 = JML_LUH * JML_LUH / 6.0;
    for (j = 0; j < JML_NLU - 1; j++)
      {
        A->cS[4 * j] = A->S[j];
        A->cS[4 * j + 1] = A->S[j + 1] - A->S[j] - h6 * (2.0 * A->S2[j] + A->S2[j + 1]);
        A->cS[4 * j + 2] = 3.0 * h6 * A->S2[j];
        A->cS[4 * j + 3] = h6 * (A->S2[j + 1] - A->S2[j]);
        A->cF[4 * j] = A->F[j];
        A->cF[4 * j + 1] = A->F[j + 1] - A->F[j] - h6 * (2.0 * A->F2[j] + A->F2[j + 1]);
        A->cF[4 * j + 2] = 3.0 * h6 * A->F2[j];
        A->cF[4 * j + 3] = h6 * (A->F2[j + 1] - A->F2[j]);
      }
  }
}

void
jml_alpha_free (jml_alpha *A)
{
  myfree (A->S);
  myfree (A->S2);
  myfree (A->F);
  myfree (A->F2);
  myfree (A->cS);
  myfree (A->cF);
  memset (A, 0, sizeof (*A));
}

static inline double
jml_spl (const double *c, double lu)
{
  const double x = (lu - JML_LU0) * (1.0 / JML_LUH);
  long j = (long) x;
  if (j > JML_NLU - 2)
    j = JML_NLU - 2;
  const double t = x - (double) j;
  c += 4 * j;
  return ((c[3] * t + c[2]) * t + c[1]) * t + c[0];
}

static const double jml_luhi = JML_LU0 + JML_LUH * (JML_NLU - 1);

inline double
jml_logS (const jml_alpha *A, double lu)
{
  const double a = A->alpha;
  if (A->one)
    return -exp (lu);
  if (lu < JML_LU0)   /* u -> 0: S = 1 - u^a / Gamma(1 + a) */
    return log1p (-exp (a * lu - lgamma (1.0 + a)));
  if (lu > jml_luhi)  /* u -> inf: S ~ u^-a / Gamma(1 - a) */
    return A->S[JML_NLU - 1] - a * (lu - jml_luhi);
  return jml_spl (A->cS, lu);
}

double
jml_logf (const jml_alpha *A, double lu)
{
  const double a = A->alpha;
  if (A->one)
    return -exp (lu);
  if (lu < JML_LU0)   /* f = u^(a-1) / Gamma(a) */
    return (a - 1.0) * lu - lgamma (a);
  if (lu > jml_luhi)  /* f ~ a u^(-a-1) / Gamma(1 - a) */
    return A->F[JML_NLU - 1] - (a + 1.0) * (lu - jml_luhi);
  return jml_spl (A->cF, lu);
}

/* ------------------------------------------------------- spells */

/// log density of one population's spells; M[j]: rate of a lineage of the
/// population moving to j (0 when there is no route). The spells carry log dt
/// (jml_spells_prepare()) and are sorted by k, so that log lambda is computed
/// once per k.
double
jml_pop_logp (const double *sp, long n, long numpop, double mu, double inh,
              double theta, const double *M, const jml_alpha *A)
{
  double mtot = 0.0, v = 0.0, lrate[3 + 2];
  long i, j, lastk = -1;
  double lam = 0.0, ll = 0.0;
  if (!(theta > 0.0))
    return -HUGE_VAL;
  for (j = 0; j < numpop; j++)
    {
      mtot += M[j];
      lrate[2 + j] = M[j] > 0.0 ? log (M[j] / mu) : -HUGE_VAL;
    }
  const double cpair = 1.0 / (mu * inh * theta);   /* k(k-1) cpair = coalescence rate */
  const double mrate = mtot / mu;
  const double a = A->alpha, ia = 1.0 / a;
  lrate[1] = log (2.0 * cpair);
  for (i = 0; i < n; i++)
    {
      const long kind = JML_PKIND (sp[2 * i]);
      const long kk = JML_PK (sp[2 * i]);
      const double ldt = sp[2 * i + 1];
      if (kk != lastk)
        {
          const double k = (double) kk;
          lam = k * (k - 1.0) * cpair + k * mrate;
          ll = lam > 0.0 ? log (lam) : -HUGE_VAL;
          lastk = kk;
        }
      if (kind >= 2 && lrate[kind] == -HUGE_VAL)
        return -HUGE_VAL;
      if (!(lam > 0.0))
        {
          if (kind > 0)
            return -HUGE_VAL;
          continue;
        }
      if (A->one)
        v += (kind > 0 ? lrate[kind] : 0.0) - lam * exp (ldt);
      else
        {
          const double lu = ll * ia + ldt;
          if (kind == 0)
            v += jml_logS (A, lu);
          else
            v += jml_logf (A, lu) + ll * ia - ll + lrate[kind];
        }
    }
  return v;
}

/* sort key: population, then k */
static int
jml_spcmp (const void *x, const void *y)
{
  const double a = *(const double *) x, b = *(const double *) y;
  return a < b ? -1 : (a > b ? 1 : 0);
}

/// spells as the evaluation needs them: dt -> log dt, sorted by
/// (population, k) (the code orders them so)
void
jml_spells_prepare (double *sp, long n, long numpop)
{
  long i;
  for (i = 0; i < n; i++)   /* a tie (dt = 0) would make an event density infinite for alpha < 1 */
    sp[2 * i + 1] = log (sp[2 * i + 1] > 1e-300 ? sp[2 * i + 1] : 1e-300);
  /* code = ((k numpop + pop)(numpop + 2) + kind): sort by (pop, k) */
  double *key = (double *) mycalloc ((size_t) (3 * (n > 0 ? n : 1)), sizeof (double));
  for (i = 0; i < n; i++)
    {
      const long code = (long) sp[2 * i];
      const long pop = JML_POP (code, numpop), k = JML_K (code, numpop);
      key[3 * i] = (double) pop * 1e9 + (double) k;
      key[3 * i + 1] = sp[2 * i];
      key[3 * i + 2] = sp[2 * i + 1];
    }
  qsort (key, (size_t) n, 3 * sizeof (double), jml_spcmp);
  for (i = 0; i < n; i++)
    {
      const long code = (long) key[3 * i + 1];
      sp[2 * i] = (double) (JML_K (code, numpop) * 8 + JML_KIND (code, numpop));
      sp[2 * i + 1] = key[3 * i + 2];
    }
  myfree (key);
}

/* ------------------------------------------------------- normalizer */

/* regularized incomplete beta I_x(a, b) (continued fraction) */
static double
jml_betacf (double a, double b, double x)
{
  const double tiny = 1e-300;
  double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0), h;
  long m;
  if (fabs (d) < tiny)
    d = tiny;
  d = 1.0 / d;
  h = d;
  for (m = 1; m <= 500; m++)
    {
      const double m2 = 2.0 * m;
      double aa = m * (b - m) * x / ((a - 1.0 + m2) * (a + m2));
      d = 1.0 + aa * d;
      if (fabs (d) < tiny) d = tiny;
      c = 1.0 + aa / c;
      if (fabs (c) < tiny) c = tiny;
      d = 1.0 / d;
      h *= d * c;
      aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1.0 + m2));
      d = 1.0 + aa * d;
      if (fabs (d) < tiny) d = tiny;
      c = 1.0 + aa / c;
      if (fabs (c) < tiny) c = tiny;
      d = 1.0 / d;
      const double del = d * c;
      h *= del;
      if (fabs (del - 1.0) < 1e-15)
        break;
    }
  return h;
}

static double
jml_ibeta (double a, double b, double x)
{
  if (x <= 0.0)
    return 0.0;
  if (x >= 1.0)
    return 1.0;
  const double lbt = lgamma (a + b) - lgamma (a) - lgamma (b) + a * log (x) + b * log1p (-x);
  if (x < (a + 1.0) / (a + b + 2.0))
    return exp (lbt) * jml_betacf (a, b, x) / a;
  return 1.0 - exp (lbt) * jml_betacf (b, a, 1.0 - x) / b;
}

/* log D(M_tot): r = 1 or 2 routes with m[] migrations, c = mmax / M_tot */
static double
jml_logD (long r, const long *m, double c)
{
  if (r <= 1)
    return c >= 1.0 ? 0.0 : -HUGE_VAL;
  const double a = (double) m[0] + 1.0, b = (double) m[1] + 1.0;
  const double lo = c < 1.0 ? 1.0 - c : 0.0, hi = c < 1.0 ? c : 1.0;
  if (hi <= lo)
    return -HUGE_VAL;
  double v;
  if (lo <= 0.0)
    return lgamma (a) + lgamma (b) - lgamma (a + b);
  /* 1 - I_lo - (1 - I_hi): the two tails, accurately for small tails */
  v = 1.0 - jml_ibeta (a, b, lo) - jml_ibeta (b, a, 1.0 - hi);
  if (!(v > 0.0))
    return -HUGE_VAL;
  return lgamma (a) + lgamma (b) - lgamma (a + b) + log (v);
}

static double
jml_logsig (double y)
{
  return y > 0.0 ? -log1p (exp (-y)) : y - log1p (exp (y));
}

static double
jml_sig (double y)
{
  return y > 0.0 ? 1.0 / (1.0 + exp (-y)) : exp (y) / (1.0 + exp (y));
}

typedef struct
{
  const double *sp;
  long n, numpop, pop;
  double mu, inh;
  const jml_alpha *A;
  jml_inner *I;
  long m[2];
  double *M;          /* numpop scratch */
  double *spr;        /* spells with every migration mapped to route 0 */
} jml_ctx;

/* log integrand in y = (logit Theta, logit M_tot / (r mmax)), with
   Jacobians, the Theta prior and the uniform M priors */
static double
jml_h (jml_ctx *c, const double *y)
{
  jml_inner *I = c->I;
  const double th = I->tlo + (I->thi - I->tlo) * jml_sig (y[0]);
  double v = jml_logsig (y[0]) + jml_logsig (-y[0]) + log (I->thi - I->tlo) + I->logprior (I->ctx, th);
  long j;
  for (j = 0; j < c->numpop; j++)
    c->M[j] = 0.0;
  if (I->r > 0)
    {
      const double span = (double) I->r * I->mmax;
      const double mt = span * jml_sig (y[1]);
      if (!(mt > 0.0))
        return -HUGE_VAL;
      c->M[I->mig_to[0]] = mt;   /* every migration at rate M_tot */
      v += jml_logsig (y[1]) + jml_logsig (-y[1]) + log (span)
        + (double) (I->r - 1) * log (mt) + jml_logD (I->r, c->m, I->mmax / mt)
        - (double) I->r * log (I->mmax);
    }
  if (v == -HUGE_VAL)
    return v;
  return v + jml_pop_logp (c->spr, c->n, c->numpop, c->mu, c->inh, th, c->M, c->A);
}

/// log of one population's inner integral at fixed alpha
double
jml_inner_log (const double *sp, long n, long numpop, long pop, double mu,
               double inh, const jml_alpha *A, jml_inner *I)
{
  jml_ctx c;
  const long d = I->r > 0 ? 2 : 1;
  double y[2], g[2], H[4], h0, L[4];
  long i, j, it, q;
  static double ghx[JML_GH], ghw[JML_GH];
  static int ghinit = 0;
  if (!ghinit)
    {   /* Gauss-Hermite nodes and weights (Golub-Welsch would be overkill:
           the 9-point rule) */
      const double x9[5] = {0.0, 0.723551018752837573, 1.468553289216667931,
                            2.266580584531843111, 3.190993201781527830};
      const double w9[5] = {0.720235215606050971, 0.432651559002555750, 0.0884745273943765732,
                            0.00494362427553694721, 0.0000396069772632643819};
      ghx[4] = x9[0];
      ghw[4] = w9[0];
      for (q = 1; q < 5; q++)
        {
          ghx[4 + q] = x9[q];
          ghx[4 - q] = -x9[q];
          ghw[4 + q] = ghw[4 - q] = w9[q];
        }
      ghinit = 1;
    }
  memset (&c, 0, sizeof (c));
  c.sp = sp;
  c.n = n;
  c.numpop = numpop;
  c.pop = pop;
  c.mu = mu;
  c.inh = inh;
  c.A = A;
  c.I = I;
  c.M = (double *) mycalloc ((size_t) numpop, sizeof (double));
  c.spr = (double *) mycalloc ((size_t) (2 * (n > 0 ? n : 1)), sizeof (double));
  for (i = 0; i < n; i++)
    {
      long kind = JML_PKIND (sp[2 * i]);
      c.spr[2 * i] = sp[2 * i];
      c.spr[2 * i + 1] = sp[2 * i + 1];
      if (kind >= 2)
        {
          for (q = 0; q < I->r; q++)
            if (I->mig_to[q] == kind - 2)
              c.m[q < 2 ? q : 1] += 1;
          c.spr[2 * i] += (double) (I->mig_to[0] - (kind - 2));
        }
    }
  /* mode: Newton in y with central differences, from the warm start (or a
     coarse search the first time) */
  if (!I->started)
    {
      double best = -HUGE_VAL;
      for (i = 0; i < 25; i++)
        for (j = 0; j < (d == 2 ? 25 : 1); j++)
          {
            double z[2] = {-12.0 + 1.0 * (double) i, d == 2 ? -14.0 + 1.0 * (double) j : 0.0};
            const double v = jml_h (&c, z);
            if (v > best)
              {
                best = v;
                I->y[0] = z[0];
                I->y[1] = z[1];
              }
          }
      I->started = 1;
    }
  y[0] = I->y[0];
  y[1] = d == 2 ? I->y[1] : 0.0;
  const double hd = 1e-3;
  h0 = jml_h (&c, y);
  for (it = 0; it < 50; it++)
    {
      double yy[2], fpp, fmm, fpm, fmp, step[2], det;
      for (q = 0; q < d; q++)
        {
          yy[0] = y[0];
          yy[1] = y[1];
          yy[q] = y[q] + hd;
          const double fp = jml_h (&c, yy);
          yy[q] = y[q] - hd;
          const double fm = jml_h (&c, yy);
          g[q] = (fp - fm) / (2.0 * hd);
          H[q * 2 + q] = (fp - 2.0 * h0 + fm) / (hd * hd);
        }
      if (d == 2)
        {
          yy[0] = y[0] + hd; yy[1] = y[1] + hd; fpp = jml_h (&c, yy);
          yy[0] = y[0] - hd; yy[1] = y[1] - hd; fmm = jml_h (&c, yy);
          yy[0] = y[0] + hd; yy[1] = y[1] - hd; fpm = jml_h (&c, yy);
          yy[0] = y[0] - hd; yy[1] = y[1] + hd; fmp = jml_h (&c, yy);
          H[1] = H[2] = (fpp - fpm - fmp + fmm) / (4.0 * hd * hd);
        }
      /* Newton on -h with the Hessian made negative definite */
      if (d == 1)
        {
          const double hh = H[0] < -1e-6 ? H[0] : -1e-6;
          step[0] = -g[0] / hh;
          step[1] = 0.0;
        }
      else
        {
          double a = H[0], b = H[1], dd = H[3];
          const double tr = a + dd, dt = a * dd - b * b;
          if (!(a < 0.0 && dt > 0.0))
            {   /* not negative definite: gradient step */
              (void) tr;
              a = dd = -1.0;
              b = 0.0;
            }
          det = a * dd - b * b;
          step[0] = -(dd * g[0] - b * g[1]) / det;
          step[1] = -(-b * g[0] + a * g[1]) / det;
        }
      for (q = 0; q < d; q++)
        {
          if (step[q] > 2.0) step[q] = 2.0;
          if (step[q] < -2.0) step[q] = -2.0;
        }
      double lam = 1.0, hn = -HUGE_VAL;
      for (q = 0; q < 30; q++)
        {
          double yn[2] = {y[0] + lam * step[0], y[1] + lam * step[1]};
          hn = jml_h (&c, yn);
          if (hn >= h0)
            {
              y[0] = yn[0];
              y[1] = yn[1];
              break;
            }
          lam *= 0.5;
        }
      const double gain = hn - h0;
      if (hn >= h0)
        h0 = hn;
      if (fabs (step[0]) + fabs (step[1]) < 1e-6 || (gain >= 0.0 && gain < 1e-10) || q == 30)
        break;
    }
  /* Hessian at the mode, its inverse's Cholesky factor */
  {
    double yy[2];
    for (q = 0; q < d; q++)
      {
        yy[0] = y[0];
        yy[1] = y[1];
        yy[q] = y[q] + hd;
        const double fp = jml_h (&c, yy);
        yy[q] = y[q] - hd;
        const double fm = jml_h (&c, yy);
        H[q * 2 + q] = (fp - 2.0 * h0 + fm) / (hd * hd);
      }
    if (d == 2)
      {
        double fpp, fmm, fpm, fmp;
        yy[0] = y[0] + hd; yy[1] = y[1] + hd; fpp = jml_h (&c, yy);
        yy[0] = y[0] - hd; yy[1] = y[1] - hd; fmm = jml_h (&c, yy);
        yy[0] = y[0] + hd; yy[1] = y[1] - hd; fpm = jml_h (&c, yy);
        yy[0] = y[0] - hd; yy[1] = y[1] + hd; fmp = jml_h (&c, yy);
        H[1] = H[2] = (fpp - fpm - fmp + fmm) / (4.0 * hd * hd);
      }
  }
  I->y[0] = y[0];
  I->y[1] = y[1];
  /* A = -H (positive definite: eigenvalues floored), C = A^-1 = L L^T */
  if (d == 1)
    {
      const double a = -H[0] > 1e-8 ? -H[0] : 1e-8;
      L[0] = 1.0 / sqrt (a);
    }
  else
    {
      double a = -H[0], b = -H[1], dd = -H[3];
      const double m = 0.5 * (a + dd), s = sqrt (0.25 * (a - dd) * (a - dd) + b * b);
      double l1 = m + s, l2 = m - s;
      if (l2 < 1e-8 || l1 < 1e-8)
        {   /* floor the eigenvalues */
          double v1x, v1y, nv;
          if (fabs (b) > 1e-300)
            {
              v1x = l1 - dd;
              v1y = b;
            }
          else
            {
              v1x = a >= dd ? 1.0 : 0.0;
              v1y = a >= dd ? 0.0 : 1.0;
            }
          nv = sqrt (v1x * v1x + v1y * v1y);
          v1x /= nv;
          v1y /= nv;
          if (l1 < 1e-8) l1 = 1e-8;
          if (l2 < 1e-8) l2 = 1e-8;
          a = l1 * v1x * v1x + l2 * v1y * v1y;
          b = (l1 - l2) * v1x * v1y;
          dd = l1 * v1y * v1y + l2 * v1x * v1x;
        }
      const double det = a * dd - b * b;
      const double c00 = dd / det, c01 = -b / det, c11 = a / det;
      L[0] = sqrt (c00);
      L[1] = 0.0;
      L[2] = c01 / L[0];
      L[3] = sqrt (c11 - L[2] * L[2] > 0.0 ? c11 - L[2] * L[2] : 1e-300);
    }
  /* adaptive Gauss-Hermite */
  double mx = -HUGE_VAL, sum = 0.0;
  double vals[JML_GH * JML_GH];
  for (i = 0; i < JML_GH; i++)
    for (j = 0; j < (d == 2 ? JML_GH : 1); j++)
      {
        const double x0 = M_SQRT2 * ghx[i], x1 = d == 2 ? M_SQRT2 * ghx[j] : 0.0;
        double z[2];
        z[0] = y[0] + L[0] * x0;
        z[1] = d == 2 ? y[1] + L[2] * x0 + L[3] * x1 : 0.0;
        double v = jml_h (&c, z);
        v += log (ghw[i]) + ghx[i] * ghx[i];
        if (d == 2)
          v += log (ghw[j]) + ghx[j] * ghx[j];
        vals[i * JML_GH + j] = v;
        if (v > mx)
          mx = v;
      }
  if (mx == -HUGE_VAL)
    {
      myfree (c.M);
      myfree (c.spr);
      return -HUGE_VAL;
    }
  for (i = 0; i < JML_GH; i++)
    for (j = 0; j < (d == 2 ? JML_GH : 1); j++)
      sum += exp (vals[i * JML_GH + j] - mx);
  myfree (c.M);
  myfree (c.spr);
  return mx + log (sum) + log (L[0]) + (d == 2 ? log (L[3]) : 0.0) + 0.5 * (double) d * LOG2;
}

/// Gauss-Legendre: 3 panels of 8 nodes on [amin, amax]
void
jml_alpha_nodes (double amin, double amax, double *x, double *w)
{
  static const double gx[4] = {0.183434642495649805, 0.525532409916328986,
                               0.796666477413626740, 0.960289856497536232};
  static const double gw[4] = {0.362683783378361983, 0.313706645877887287,
                               0.222381034453374471, 0.101228536290376259};
  long p, q, n = 0;
  const double h = (amax - amin) / 3.0;
  for (p = 0; p < 3; p++)
    {
      const double c = amin + h * (p + 0.5);
      for (q = 0; q < 4; q++)
        {
          x[n] = c - 0.5 * h * gx[q];
          w[n++] = 0.5 * h * gw[q];
          x[n] = c + 0.5 * h * gx[q];
          w[n++] = 0.5 * h * gw[q];
        }
    }
}
