// marginal likelihood summaries in migrate
// MIT opensource license
// consolidation from other files to simplify future changes
// (c) Peter Beerli 2021
//

#include "migration.h"
#include "sighandler.h"
#include "bayes.h"
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
  const long np = world->numpop2 + world->species_model_size * 2 + world->grownum;
  const long np1 = np -  world->grownum;
  long pop;
  long i;
  MYREAL scaling_factor=0.0;
  bayes_fmt * bayes = world->bayes;
  boolean *visited;
  double w=1.0;
  double v=0.0;
  double pr=0.0;
  visited = (boolean *) mycalloc(np,sizeof(boolean));
  for(i=0;i<np;i++)
    {
      if(i<np1)
	{
      if(bayes->map[i][1] == INVALID)
	{
	  if (bayes->custm2[i]=='c')
	    {
	      //scaling_factor += logpriors[i][bin_c]
	      w = world->bayes->deltahist[i];
	      //double v = bayes->histogram[0].minima[i] +  + w/2
	      v = ((long) world->param0[i]/w)  + w/2.;
	      pr = scaling_prior(world,i,v);
	      scaling_factor += (1.0-world->loci) * (log(w) + pr);
#ifdef DEBUG
	      printf("%i> scaling factor with 'c': %li k=%f log(w)=%f  w=%f v=%f pr=%f  [%f]\n",myID, i, scaling_factor, log(w), w, v, pr, world->param0[i]);
#endif

	    }
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
	      warning("am or likelihood failed: am=%f",myID, locus, i,world->am[locus]);
	    }
	  world->steppingstones[ii] = universe[i]->steppingstones[ii];
	  world->steppingstone_scalars[ii] = universe[i]->steppingstone_scalars[ii];
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
    *bufsize += sprintf(buf+ *bufsize,"# --  %s = %f\n", "Thermodynamic temperature", world->options->heat[i]);
  *bufsize += sprintf(buf+ *bufsize,"# --  %s\n", "Marginal log(likelihood) [Thermodynamic integration]");
  *bufsize += sprintf(buf+ *bufsize,"# --  %s\n", "Marginal log(likelihood) [Harmonic mean]");
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
      printf("@MARGLIKE %f %f\n@",  world->bf[locus * hc + t-1], temp[(*z)-2]);
#endif
    }
  temp[(*z)++] =  (float) (world->hmscale[locus] - log(world->hm[locus]));
  for(t=0; t < hc; t++)
    {
      temp[(*z)++] = (float) world->steppingstones[locus * hc + t];
      temp[(*z)++] = (float) world->steppingstone_scalars[locus * hc + t];
    }
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
	  *c += sprintf(temp+ *c,"\t%f", world->bf[locus * hc + t-1]);
	  lsum += (heat0 - heat1) * ((world->bf[locus * hc + t-1] + world->bf[locus * hc + t]) * 0.5);
	}
      *c += sprintf(temp + *c,"\t%f", world->bf[locus * hc + t-1]);
      *c += sprintf(temp + *c,"\t%f", lsum);
    }
  *c += sprintf(temp + *c,"\t%f", world->hmscale[locus] - log(world->hm[locus]));
  for(t=0; t < hc; t++)
    {
      *c += sprintf(temp+ *c,"\t%f", world->steppingstones[locus * hc + t]);
      *c += sprintf(temp+ *c,"\t%f",world->steppingstone_scalars[locus * hc + t]);
    }
}
#endif /*not MPI*/
