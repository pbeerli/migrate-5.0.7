// marginal likelihood summaries in migrate
// MIT opensource license
// consolidation from other files to simplify future changes
// (c) Peter Beerli 2021
//

#include "migration.h"
#include "sighandler.h"
#include "bayes.h"
#include "tools.h"
extern int myID;

void calculate_BF(world_fmt **universe, option_fmt *options);
MYREAL combine_scaling_factor(world_fmt *world);
void  print_marginal_order(char *buf, long *bufsize, world_fmt *world);

#if defined(MPI) && !defined(PARALIO)
void      print_marginal_like(float *temp, long *z, world_fmt * world);
#else /*not MPI*/
void      print_marginal_like(char *temp, long *c, world_fmt * world);
#endif

MYREAL sumbezier(long intervals, MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1, MYREAL x2, MYREAL y2, MYREAL *ratio);



MYREAL combine_scaling_factor(world_fmt *world)
{ 
  const long np = world->numparam;//world->numpop2 + world->species_model_size * 2 + world->grownum;
  const long np1 = world->numparamcumvec[SPLITSTDPRIOR];//np -  world->grownum;
  long pop;
  long i;
  MYREAL scaling_factor=0.0;
  bayes_fmt * bayes = world->bayes;
  boolean *visited;
  visited = (boolean *) mycalloc(np,sizeof(boolean));
  for(i=0;i<np;i++)
    {
      if(i<np1)
	{
      if(bayes->map[i][1] == INVALID)
	{
	  /* '0' and constant 'c' entries are fixed: no prior, nothing to
	     integrate, so they add nothing to the scaling factor */
	  continue;
	}
      else
	{
	  pop  = bayes->map[i][1];
	    }
	}
      else
	{
	  pop = i;
	}
      if(visited[pop]==TRUE)
	continue;
      visited[pop] = TRUE;
      //      scaling_factor += exp(bayes->scaling_factors[pop] - bayes->maxmaxvala);
      scaling_factor += bayes->scaling_factors[pop]; //PRODUCT_parameters(scalinfactorcalculation_see_bayes.c)
#ifdef DEBUG
      printf("%i> scaling factor test: %li  %li k=%f k_pop=%f %f\n",myID, i, pop, scaling_factor, bayes->scaling_factors[pop], bayes->maxmaxvala);
#endif
    }
  //  scaling_factor = log(scaling_factor) + bayes->maxmaxvala;
  if(world->options->has_bayesfile)
    {
#ifdef DEBUG
      printf("# Scaling factor %20.20f\n",scaling_factor);
#endif
      fprintf(world->bayesfile, "# Scaling factor %20.20f\n",scaling_factor);
    }
  myfree(visited);
  return scaling_factor;
}



/// integrates over a Bezier curve between two points
/// calculates two handle points that are set to adhoc values
/// so that the x values of the handle are the the x value of the lowest point
/// and the y values are set to about 80% of the min to max interval for the left point
/// and a value that is the the y value from ax + b where a is calculated from a 
/// third point to the right and the second point, the third point is not used for the
/// the Bezier curve otherwise
MYREAL sumbezier(long intervals, MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1, MYREAL x2, MYREAL y2, MYREAL *ratio)
{
  const MYREAL inv_interval = 1./intervals;
  const MYREAL sx0 = x0;
  const MYREAL sx1 = x0;
  const MYREAL sy0 = 0.2 * y0 + 0.8 * y2;
  const MYREAL sy1 = (-x2 * y1 + x1 * y2)/(x1 - x2);
  MYREAL t     = 0.0;
  MYREAL t2    = 0.0;
  MYREAL t3    = 0.0;
  MYREAL onet  = 1.0;
  MYREAL onet2 = 1.0;
  MYREAL onet3 = 1.0;
  MYREAL newx;
  MYREAL newy;
  MYREAL oldx;
  MYREAL oldy;
  MYREAL sum = 0.0;
  // integrate over intervals between x0 and x1 and return sum
  // intialize with t=0.0
  oldx = x0;
  oldy = y0;
  //fprintf(stdout,"\n\n%f %f %f %f %f %f\n",x2,y2,sx0,sy0, x0,y0);
  for(t=inv_interval; t <= 1.0; t += inv_interval)
    {
      onet  = 1.0 - t;
      onet2 = onet * onet;
      onet3 = onet2 * onet;
      onet2 *= 3.0 * t;
      t2 = t * t;
      t3 = t2 * t;
      t2 *=  3.0 * onet;
      //      newx = 3sx0 (1-t)^2 t + 3 sx1 (1-t) t^2 + (1-t)^3 x0 + t^3 x1
      newx = sx0 * onet2 + sx1 * t2 + onet3 * x0 + t3 * x1;
      newy = sy0 * onet2 + sy1 * t2 + onet3 * y0 + t3 * y1;
      //fprintf(stdout,"%f %f\n",newx,newy);
      //printf("\"log mL:\", %i, %f, %f, %f, %f\n", myID, newx, oldx, newy, oldy); 
      sum += (newx - oldx) * (newy + oldy)/2.;
      *ratio += oldy - newy;
      oldx = newx;
      oldy = newy;
    }
#ifdef DEBUG
  //fprintf(stdout,"%f %f %f %f %f %f sum=%f (sum_nobezier %f)\n\n\n",x0,y0,x1,y1,sx1,sy1,sum,(x1-x0)*(y1-y0)/2.0);
#endif
  return sum;
}


/// calculate values for the marginal likelihood using thermodynamic integration
/// based on a method by Friel and Pettitt 2005
/// (http://www.stats.gla.ac.uk/research/TechRep2005/05.10.pdf)
/// this is the same method described in Lartillot and Phillippe 2006 Syst Bio
/// integrate over all temperature using a simple trapezoidal rule
/// prob(D|model) is only accurate with intervals for temperatures from 1/0 to 1/1.
/// reports also the harmonic mean
void calculate_BF(world_fmt **universe, option_fmt *options)
{
  long i;
  world_fmt * world = universe[0];
  MYREAL xx, xx2;
  long locus = universe[0]->locus;
  long hc = options->heated_chains;
  if(world->data->skiploci[locus])
    return;
  if(world->likelihood[world->G] <= (double) -HUGE)
    return;
  //am contains the counter
  world->am[locus] += 1;
  //locus = world->locus;
  xx = world->likelihood[world->G];
  if(xx <= (double) -HUGE)
    {
      warning("%i> l=%li likelihood < -HUGE", myID, locus);
      return;
    }
  if(xx > world->hmscale[locus])
    {
      xx2 = EXP(world->hmscale[locus] - xx);
      world->hm[locus] += (xx2 - world->hm[locus])/ (world->am[locus]);
    }
  else
    {
      world->hm[locus] *= EXP(xx - world->hmscale[locus]);
      world->hmscale[locus] = xx;
      world->hm[locus] += (1. - world->hm[locus])/ (world->am[locus]);
    }
  //thermodynamic section: calculates one-pass averages of the loglike for the different temperatures
  //stored in the cold chain
  if(options->heating)
    {
      //#ifdef DEBUG
      //printf("%i>BF: %li*4*i:",myID, locus);
      //#endif 
      for (i = 0; i < hc; i++)
	{
	  long ii = locus * hc + i;
	  xx = universe[i]->likelihood[universe[i]->G];
	  if (world->am[locus] > 0.0 || xx > (double) -HUGE)
	    world->bf[ii] += (xx - world->bf[ii])/ (world->am[locus]);
	  else
	    {
	      warning("[%i] locus %li chain %li: TI sample skipped (am=%f, log L=%g)\n", myID, locus, i, world->am[locus], xx);
	    }
	  if (options->adaptiveheat != NOTADAPTIVE && xx > (double) -HUGE)
	    {   /* the chain's temperature changes: bin the sample by its beta */
	      const double beta = universe[i]->heat;
	      long bin = (long) (pow (beta > 0.0 ? beta : 0.0, 0.25) * TI_NBINS);
	      double *tb;
	      if (bin >= TI_NBINS)
		bin = TI_NBINS - 1;
	      if (bin < 0)
		bin = 0;
	      tb = world->tibins + 3 * (locus * TI_NBINS + bin);
	      tb[0] += 1.0;
	      tb[1] += xx;
	      tb[2] += beta;
	    }
	  /* stepping stones (Xie et al. 2011): for every chain but the cold
	     one, the running mean of L^(beta_{i-1} - beta_i) from the hotter
	     chain i, as steppingstones[ii] * exp(steppingstone_scalars[ii])
	     (rescaled to the largest term) */
	  if (i > 0)
	    {
	      const double x = (universe[i-1]->heat - universe[i]->heat) * xx;
	      const double n = world->am[locus];
	      if (n <= 1.0)
		{
		  world->steppingstone_scalars[ii] = x;
		  world->steppingstones[ii] = 1.0;
		}
	      else if (x > world->steppingstone_scalars[ii])
		{
		  world->steppingstones[ii] *= EXP(world->steppingstone_scalars[ii] - x);
		  world->steppingstone_scalars[ii] = x;
		  world->steppingstones[ii] += (1.0 - world->steppingstones[ii]) / n;
		}
	      else
		world->steppingstones[ii] += (EXP(x - world->steppingstone_scalars[ii]) - world->steppingstones[ii]) / n;
	    }
	  else
	    {
	      world->steppingstones[ii] = 1.0;   /* log 1 = 0: no ratio for the cold chain */
	      world->steppingstone_scalars[ii] = 0.0;
	    }
#ifdef DEBUG
	  //  printf("%f ",world->bf[locus * hc + i]); 
#endif
	  if (isnan(world->bf[locus * hc + i]))
	    {
	      world->data->skiploci[locus] = TRUE;
	      world->bf[locus * hc + i] = 0.0;
	    }
	}
#ifdef DEBUG
      //printf("\n");
#endif
    }
}

void  print_marginal_order(char *buf, long *bufsize, world_fmt *world)
{
  long i;

  for(i=0;i<world->options->heated_chains;i++)
    *bufsize += mysnprintf(buf+ *bufsize,LINESIZE,"# --  %s = %f\n", "Thermodynamic temperature", world->options->heat[i]);
  *bufsize += mysnprintf(buf+ *bufsize,LINESIZE,"# --  %s\n", "Marginal log(likelihood) [Thermodynamic integration]");
  *bufsize += mysnprintf(buf+ *bufsize,LINESIZE,"# --  %s\n", "Marginal log(likelihood) [Harmonic mean]");
}

#if defined(MPI) && !defined(PARALIO) /* */

void      print_marginal_like(float *temp, long *z, world_fmt * world)
{
  long locus = world->locus;
  long t;
  long hc = world->options->heated_chains; 
  MYREAL lsum; 
  MYREAL heat0, heat1;

  if(world->options->heating)
    {
      lsum = 0.;
      for(t=1; t < hc; t++)
	{
	  heat0 = 1./world->options->heat[t-1];
	  heat1 = 1./world->options->heat[t];
	  // this ignores adaptive heating for MPI!!!!
	  temp[*z] = (float) world->bf[locus * hc + t-1];
	  *z += 1;
	  lsum += (heat0 - heat1) * ((world->bf[locus * hc + t-1] + world->bf[locus * hc + t]) * 0.5);
	}
      temp[(*z)++] =  (float) world->bf[locus * hc + t-1];
      temp[(*z)++] =  (float) lsum;
#ifdef DEBUG
      //printf("@MARGLIKE %f %f\n@",  world->bf[locus * hc + t-1], temp[(*z)-2]);
#endif
    }
  temp[(*z)++] =  (float) (world->hmscale[locus] - log(world->hm[locus]));
#ifdef STEPPINGSTONE
  for(t=0; t < hc; t++)
    {
      temp[(*z)++] = (float) world->steppingstones[locus * hc + t];
      temp[(*z)++] = (float) world->steppingstone_scalars[locus * hc + t];
    }
#endif
}
#else /*not MPI or MPI & PARALIO*/
void      print_marginal_like(char *temp, long *c, world_fmt * world)
{
  long locus = world->locus;
  long t;
  long hc = world->options->heated_chains;  
  MYREAL lsum;
  MYREAL heat0, heat1;
  if(world->options->heating)
    {
      lsum = 0.;
      for(t=1; t < hc; t++)
	{
	  if(world->options->adaptiveheat!=NOTADAPTIVE)
	    {
	      heat0 = world->options->averageheat[t-1];
	      heat1 = world->options->averageheat[t];
	    }
	  else
	    {
	      heat0 = 1./ world->options->heat[t-1];
	      heat1 = 1./ world->options->heat[t];
	    }
	  *c += mysnprintf(temp+ *c,LINESIZE,"\t%f", world->bf[locus * hc + t-1]);
	  lsum += (heat0 - heat1) * ((world->bf[locus * hc + t-1] + world->bf[locus * hc + t]) * 0.5);
	}
      *c += mysnprintf(temp + *c,LINESIZE,"\t%f", world->bf[locus * hc + t-1]);
      *c += mysnprintf(temp + *c,LINESIZE,"\t%f", lsum);
    }
  *c += mysnprintf(temp + *c,LINESIZE,"\t%f", world->hmscale[locus] - log(world->hm[locus]));
#ifdef STEPPINGSTONE
  for(t=0; t < hc; t++)
    {
      *c += mysnprintf(temp+ *c,LINESIZE,"\t%f", world->steppingstones[locus * hc + t]);
      *c += mysnprintf(temp+ *c,LINESIZE,"\t%f",world->steppingstone_scalars[locus * hc + t]);
    }
#endif
}
#endif /*not MPI*/

/// stepping-stone log marginal likelihood of one locus: the sum over the
/// heated chains of log mean L^(beta_{i-1} - beta_i) (samples of the
/// hotter chain i), plus beta_min E[log L] for the step from the hottest
/// chain to beta = 0
double ss_locus_logml (world_fmt *world, long locus)
{
  const long hc = world->options->heated_chains;
  double sum = 0.0;
  long i;
  for (i = 1; i < hc; i++)
    sum += log (world->steppingstones[locus * hc + i]) + world->steppingstone_scalars[locus * hc + i];
  if (hc > 0 && world->options->heat[hc - 1] > 0.0)
    sum += world->bf[locus * hc + hc - 1] / world->options->heat[hc - 1];
  return sum;
}

/// adaptive heating: thermodynamic integration over the samples binned by
/// their beta (the chains' temperatures change during the run, so one mean
/// of log L per chain, placed at its average temperature, is biased): the
/// trapezoid over the bin means from beta 1 down to the hottest bin, and the
/// Bezier version that replaces the hottest segment as sumbezier() does
/// with the chains. FALSE with fewer than two occupied bins.
boolean ti_binned (world_fmt *world, long locus, double *ti, double *bti)
{
  double b[TI_NBINS], l[TI_NBINS], last = 0.0, ratio2 = 0.0;
  long k, n = 0;
  for (k = TI_NBINS - 1; k >= 0; k--)   /* from beta 1 down */
    {
      const double *tb = world->tibins + 3 * (locus * TI_NBINS + k);
      if (tb[0] > 0.0)
        {
          b[n] = tb[2] / tb[0];
          l[n] = tb[1] / tb[0];
          n++;
        }
    }
  if (n < 2)
    return FALSE;
  /* the strips from the top bin to beta 1 and from the bottom bin to
     beta 0, with the end bins' mean log L */
  *ti = (1.0 - b[0]) * l[0] + b[n - 1] * l[n - 1];
  for (k = 0; k + 1 < n; k++)
    {
      last = (b[k] - b[k + 1]) * 0.5 * (l[k] + l[k + 1]);
      *ti += last;
    }
  *bti = *ti;
  if (n >= 3)
    *bti = *ti - last + sumbezier (100L, b[n - 1], l[n - 1], b[n - 2], l[n - 2], b[n - 3], l[n - 3], &ratio2);
  return TRUE;
}

/* ---- locus-level checkpoint of the marginal-likelihood accumulators ----
   With bayes-allfile and recover=YES a restarted run skips the loci whose
   samples are all in the bayes-allfile; their posterior histograms are read
   back from it, and the running TI means are in its records, but not the
   stepping-stone sums or the beta bins of adaptive heating. Whoever finishes
   a locus (serial run, MPI locus worker after merging the replicates)
   therefore writes them to "<bayesallfile>.ckpt.<locus>"; a recovered run
   reads them back for the loci it skipped (ckpt_restore(), called before
   the marginal-likelihood tables, after the MPI results arrived). */
static void ckpt_filename (option_fmt *options, long locus, char *name)
{
  snprintf (name, LINESIZE, "%s.ckpt.%li", options->bayesmdimfilename, locus);
}

/// TRUE when a recovered run skipped every replicate of this locus
boolean ckpt_locus_skipped (option_fmt *options, long locus)
{
  const long repmax = number_replicates2 (options);
  long r;
  if (!options->checkpointing || options->unfinished == NULL)
    return FALSE;
  for (r = 0; r < repmax; r++)
    if (options->unfinished[locus][r] < options->lsteps - 1)
      return FALSE;
  return TRUE;
}

/// writes the marginal-likelihood accumulators of a finished locus
void ckpt_write_locus (world_fmt *world, option_fmt *options, long locus)
{
  const long hc = world->options->heated_chains;
  char name[LINESIZE];
  FILE *f;
  long i;
  if (!options->has_bayesmdimfile || !options->heating || ckpt_locus_skipped (options, locus))
    return;
  ckpt_filename (options, locus, name);
  if ((f = fopen (name, "w")) == NULL)
    return;
  fprintf (f, "migrate-ckpt 1 %li %li %d\n", locus, hc, TI_NBINS);
  fprintf (f, "%.17g %.17g %.17g\n", world->am[locus], world->hmscale[locus], world->hm[locus]);
  for (i = 0; i < hc; i++)
    fprintf (f, "%.17g %.17g %.17g\n", world->bf[locus * hc + i], world->steppingstones[locus * hc + i],
             world->steppingstone_scalars[locus * hc + i]);
  for (i = 0; i < 3 * TI_NBINS; i++)
    fprintf (f, "%.17g\n", world->tibins[3 * TI_NBINS * locus + i]);
  fclose (f);
}

/// a fresh run (no recover) removes checkpoint files of an earlier run
void ckpt_remove_all (option_fmt *options, long loci)
{
  char name[LINESIZE];
  long locus;
  if (!options->has_bayesmdimfile || options->checkpointing)
    return;
  for (locus = 0; locus < loci; locus++)
    {
      ckpt_filename (options, locus, name);
      remove (name);
    }
}

/// a recovered run: the accumulators of the skipped loci from their files
void ckpt_restore (world_fmt *world, option_fmt *options)
{
  const long hc = world->options->heated_chains;
  char name[LINESIZE];
  long locus, i;
  if (!options->checkpointing || !options->heating)
    return;
  for (locus = 0; locus < world->loci; locus++)
    {
      FILE *f;
      long l2, hc2, nb2, version;
      if (!ckpt_locus_skipped (options, locus))
        continue;
      ckpt_filename (options, locus, name);
      if ((f = fopen (name, "r")) == NULL
          || fscanf (f, "migrate-ckpt %li %li %li %li", &version, &l2, &hc2, &nb2) != 4
          || l2 != locus || hc2 != hc || nb2 != TI_NBINS)
        {
          if (f != NULL)
            fclose (f);
          warning ("recover: no checkpoint of the marginal-likelihood sums for locus %li (%s);"
                   " its stepping-stone value is not available\n", locus + 1, name);
          continue;
        }
      if (fscanf (f, "%lf %lf %lf", &world->am[locus], &world->hmscale[locus], &world->hm[locus]) != 3)
        warning ("recover: %s is incomplete\n", name);
      for (i = 0; i < hc; i++)
        if (fscanf (f, "%lf %lf %lf", &world->bf[locus * hc + i], &world->steppingstones[locus * hc + i],
                    &world->steppingstone_scalars[locus * hc + i]) != 3)
          warning ("recover: %s is incomplete\n", name);
      for (i = 0; i < 3 * TI_NBINS; i++)
        if (fscanf (f, "%lf", &world->tibins[3 * TI_NBINS * locus + i]) != 1)
          warning ("recover: %s is incomplete\n", name);
      fclose (f);
    }
}
