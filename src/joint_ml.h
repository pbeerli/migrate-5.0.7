#ifndef JOINT_ML_H
#define JOINT_ML_H
/*------------------------------------------------------------------------
  Mittag-Leffler alpha in the joint multi-locus combination: function
  table, spell densities and per-genealogy normalizers (joint_ml.c,
  docs/joint_mlalpha_design.md sections 6-7)
------------------------------------------------------------------------*/

/* the table collapsed to one alpha: cubic splines in log u */
typedef struct
{
  double alpha;
  int one;          /* alpha treated as 1: exponential waiting times */
  double *S, *S2;   /* log S and its spline second derivatives */
  double *F, *F2;   /* log f */
  double *cS, *cF;  /* per interval: cubic coefficients in t = (lu - lu_j) / h */
} jml_alpha;

extern void jml_table_init (double amin, double amax);
extern void jml_alpha_set (jml_alpha *A, double alpha);
extern void jml_alpha_free (jml_alpha *A);
extern double jml_logS (const jml_alpha *A, double lu);
extern double jml_logf (const jml_alpha *A, double lu);

/* spell code ((k * numpop + pop) * (numpop + 2) + kind), see mlh_spells() */
#define JML_KIND(code, numpop) ((long) (code) % ((numpop) + 2))
#define JML_POP(code, numpop) (((long) (code) / ((numpop) + 2)) % (numpop))
#define JML_K(code, numpop) (((long) (code) / ((numpop) + 2)) / (numpop))

/* prepared spells (jml_spells_prepare()): code k * 8 + kind, log dt */
#define JML_PK(code) (((long) (code)) >> 3)
#define JML_PKIND(code) (((long) (code)) & 7)

/* spells as stored -> as evaluated: log dt, sorted by (population, k) */
extern void jml_spells_prepare (double *sp, long n, long numpop);

/* log p(spells of one population | theta, M[j] (a lineage moving to j),
   alpha), the prepared spells sp[0..2n) of that population */
extern double jml_pop_logp (const double *sp, long n, long numpop, double mu, double inh,
                            double theta, const double *M, const jml_alpha *A);

/* one population's inner integral at fixed alpha over Theta (prior given by
   logprior on [tlo, thi]) and its r free immigration rates, each U(0, mmax)
   (reduced to M_tot); mig_to[q] the target population of route q */
typedef struct
{
  double tlo, thi;
  double (*logprior) (void *ctx, double theta);
  void *ctx;
  long r;
  long mig_to[8];
  double mmax;
  double y[2];      /* warm start (logit coordinates) */
  int started;
} jml_inner;

extern double jml_inner_log (const double *sp, long n, long numpop, long pop, double mu,
                             double inh, const jml_alpha *A, jml_inner *I);

/* Gauss-Legendre nodes for the alpha integral */
#define JML_NALPHA 24
extern void jml_alpha_nodes (double amin, double amax, double *x, double *w);

#endif /*JOINT_ML_H*/
