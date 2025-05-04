//
// mutation models
//
//
//
/*
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
#include "migration.h"
#include "sequence.h"
#include "tree.h"
#include "tools.h"
#include "data.h"

#include "sighandler.h"
#include "watterson.h"
#include "mutationmodel.h"
extern int myID;


void calc_brownian_default(mutationmodel_fmt *s, option_fmt *options, data_fmt *data, long locus, long sublocus);
void calc_seq_basefreq(mutationmodel_fmt *s, option_fmt *options, data_fmt *data, long locus, long sublocus);
void set_subloci_basefrequencies_seq(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long locus, long xs);
void calculate_microsat_steps (mutationmodel_fmt *s);
void adjust_mutationmodel(mutationmodel_fmt *s, option_fmt *options);
void set_subloci_basefrequencies(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long locus, long xs);
void set_tn_model(double a1, double a2, double b, double pia, double pic, double pig, mutationmodel_fmt *mumod);
void set_mutationmodel_eigenmaterial(long z, world_fmt *world);
void copy_micro_steps(mutationmodel_fmt *s, mutationmodel_fmt *old);
void prob_tn93 (MYREAL p[4][4], const MYREAL u, const MYREAL ar, const MYREAL ay, const MYREAL b, const mutationmodel_fmt *s);

long nux(long offspring,long rate, long nuc1, long nuc2, long numsiterates,long numnuc);
#ifndef AVX
void multiply_add(MYREAL p[4][4], const MYREAL *xx1, MYREAL h[4], long numstates);
#else
void multiply_add(double p[4][4], double *xx1, double h[4], long numstates);
void avx_like(double p1[4][4], double p2[4][4], double *xx1, double *xx2, double *xx3);
#endif
void multiply_add_general(MYREAL p[4][4], const MYREAL *xx1, MYREAL h[4], long numstates);
//##
///
/// returns the number of states for a specific mutation model
int get_states(mutationmodel_fmt *s, data_fmt *data, long locus)
{
  long nummaxallele = 4;
  switch(s->datatype)
    {
    case 'a':
    case 'b':
    case 'm':
      //nummaxallele = data->maxalleles[locus];
      nummaxallele = findAllele(data,"\0",locus);
      break;
    case 'n':
    case 'u':
    case 'h':
      //nummaxallele = findAllele(data,"\0",locus);
      //break;
    case 's':
      nummaxallele = 4;
      break;
    }
  return (int) nummaxallele;
}


void calculate_microsat_steps (mutationmodel_fmt *s)
{
  long k;
  long diff;
  const long stepnum = s->micro_threshold;
  MYREAL **steps = s->steps;

  if(stepnum < 1)
    error("micro_threshold is zero or not initialized");

//printf("\n%i> micro-threshold=%li %li\n", myID, s->micro_threshold, stepnum);
  for (diff = 0; diff < stepnum; diff++)
    {
      for (k = 0; k < stepnum; k++)
	{
	  steps[diff][k] = logfac(k + diff) + logfac (k);
	  //printf("%i> steps[%li][%li]=%f\n",myID, diff, k, steps[diff][k]);
	}
    }
}

///
/// calculates the defaults for some datatypes
void set_subloci_basedefaults(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long sublocus)
{
  //long   top;
  //long   pop;
  //long   ind;
  //long   ii;
  //long   btotal=0;
  //MYREAL bsum=0.;
  long locus = world->locus;
  //char a1[DEFAULT_ALLELENMLENGTH];
  //  char a2[DEFAULT_ALLELENMLENGTH];

  //long smallest;
  //long biggest;

  //MYREAL freqa, freqc, freqg, freqt, freqy, freqr;
  switch(s->datatype)
    {
    case 'a':
      // already dealt with in create_alleles
      if (s->maxalleles == 0)
	s->maxalleles = data->maxalleles[sublocus];
      break;
    case 'b':
      calc_brownian_default(s, options, data, locus, sublocus);
      s->maxalleles = XBROWN_SIZE;
      break;
    case 'm':
      calculate_microsat_steps(s);
      break;
    case 'n':
    case 's':
      break;
    default:
      break;
    }
}

void adjust_mutationmodel(mutationmodel_fmt *s, option_fmt *options)
{
  if (s->from_infile == FALSE)
    {
      if(strchr(SEQUENCETYPES,options->datatype) || options->datatype == '@')
	{
	  s->model = options->sequence_model;
	  if(!options->freqsfrom)
	    {
	      if(s->basefreqs==NULL)
		s->basefreqs = (double*) mycalloc((s->numstates+BASEFREQLENGTH-4), sizeof(double));
	      s->basefreqs[NUC_A] = options->freqa;
	      s->basefreqs[NUC_C] = options->freqc;
	      s->basefreqs[NUC_G] = options->freqg;
	      s->basefreqs[NUC_T] = options->freqt;
	    }
	  // this works for HKY, F84 who only need ttratio
	  // and assumes that ttratio is kappa1 for TN, kappa2 is fed into param[2]!!
	  // F81, JC, K2P
	  if (s->parameters == NULL)
	    	s->parameters  = (double *) mycalloc(NUMMUTATIONPARAMETERS,sizeof(double));
	  s->parameters[0] = options->sequence_model_parameters[0]; //ttratio
	  s->ttratio = options->sequence_model_parameters[0]; //ttratio
	  s->parameters[1] = 1.0;
	  if (options->sequence_model_parameters[1] > 0.0)
	    s->parameters[2] = options->sequence_model_parameters[1];
	  else
	    s->parameters[2] = 1.0;
	}
    }
  if(s->numsiterates>1)
    {
      if(options->lambda > 0 && options->lambda != 1.0)
	s->lambda = options->lambda;
      else
	s->lambda = 1.0;
    }
}


///
/// calculates (if necessary) and sets the base frequencies
void set_subloci_basefrequencies(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long locus, long xs)
{
  (void) world;
  //long   top;
  //long   pop;
  //long   ind;
  //long   ii;
  //long   btotal=0;
  //MYREAL bsum=0.;
  //long locus = world->locus;
  //char a1[DEFAULT_ALLELENMLENGTH];
  //char a2[DEFAULT_ALLELENMLENGTH];

  //long smallest;
  //long biggest;

  //MYREAL freqa, freqc, freqg, freqt, freqy, freqr;
  switch(s->datatype)
    {
    case 'a':
      break;
    case 'b':
      calc_brownian_default(s,options, data,locus,xs);
      break;
    case 'm':
      break;
    case 'n':
    case 'h': //inserted 2021 debugging problem with hapmap data
    case 's':
      //set_subloci_basefrequencies_seq(s, world, options, data, xs);
      s->numstates = get_states(s,data,xs);
      if(s->basefreqs==NULL)
	s->basefreqs = (double*) mycalloc((s->numstates+BASEFREQLENGTH-4), sizeof(double));
      if(options->freqsfrom)
	{
	  calc_seq_basefreq(s, options, data, locus, xs);
	}
      else
	{
	  check_basefreq(options);
	  s->basefreqs[NUC_A] = options->freqa;
	  s->basefreqs[NUC_C] = options->freqc;
	  s->basefreqs[NUC_G] = options->freqg;
	  s->basefreqs[NUC_T] = options->freqt;
	}
      break;
    default:
      break;
    }
}


void calc_brownian_default(mutationmodel_fmt *s, option_fmt *options, data_fmt *data, long locus, long sublocus)
{
  long   top;
  long   pop;
  long   ind;
  long   ii;
  long   btotal=0;
  MYREAL bsum=0.;
  char a1[DEFAULT_ALLELENMLENGTH];
  char a2[DEFAULT_ALLELENMLENGTH];
  for (pop = 0; pop < data->numpop; pop++)
    {
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][locus][ii];
	  //error("not fit for subloci");
	  strcpy (a1, data->yy[pop][ind][sublocus][0][0]);
	  if (strcmp (a1, "?"))
	    {
	      btotal++;
	      bsum += atof (a1);
	    }
	  if(!data->oneliner)
	    {
	      strcpy (a2, data->yy[pop][ind][sublocus][1][0]);
	      if (strcmp (a2, "?"))
		{
		  btotal++;
		  bsum += atof (a2);
		}
	    }
	}
    }
  bsum /= btotal;
  s->browniandefault = bsum;
}

void calc_seq_basefreq(mutationmodel_fmt *s, option_fmt *options, data_fmt *data, long locus, long sublocus)
{
  long   top;
  long   pop;
  long   ind;
  long   ii;
  long xsite;
  MYREAL freqa = 0.0;
  MYREAL freqc = 0.0;
  MYREAL freqg = 0.0;
  MYREAL freqt = 0.0;
  MYREAL total = 0.0;
  for (pop = 0; pop < data->numpop; pop++)
    {
#ifdef DEBUG
      printf("%i> Locus %li  Basefrequency calculations:\n",myID,locus);
#endif
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][locus][ii];
	  //printf("\n");
	  for (xsite=0; xsite < s->numsites; xsite++)
	    {
	      //printf("%i> xsite=%li ",myID, xsite);
	      fflush(stdout);
	      char site = data->yy[pop][ind][sublocus][0][xsite][0];
	      total += 1.0;
	      //printf("%c",site);
	      //printf(" +%c+\n",site);
	      switch (site)
		{
		case 'A': freqa += 1.; break;  
		case 'C': freqc += 1.; break;  
		case 'G': freqg += 1.; break;  
		case 'U':
		case 'T': freqt += 1.; break;  
		case 'M': freqa += 0.5; freqc += 0.5; break;  
		case 'R': freqa += 0.5; freqg += 0.5; break;  
		case 'W': freqa += 0.5; freqt += 0.5; break;  
		case 'S': freqc += 0.5; freqg += 0.5; break;  
		case 'Y': freqc += 0.5; freqt += 0.5; break;  
		case 'K': freqg += 0.5; freqt += 0.5; break;  
		case 'B': freqc += 1./3.; freqg += 1./3.; freqt += 1./3.; break;  
		case 'D': freqa += 1./3.; freqg += 1./3.; freqt += 1./3.; break;  
		case 'H': freqa += 1./3.; freqc += 1./3.; freqt += 1./3.; break;  
		case 'V': freqa += 1./3.; freqc += 1./3.; freqg += 1./3.; break;  
		case 'X':  // we count X and N because this means that something is there
		case 'N': freqa += 0.25; freqc += 0.25; freqg += 0.25; freqt += 0.25; break;  
		case '-':
		case '0':
		case '?': total -= 1; //we do not count question mark, which means no data.
		  break;
		default:
		  warning("Wrong Base detected=%c pop=%li, sublocus=%li, total=%f",site,pop,sublocus,total);
		  error("Abort");
		}
	    }
	}
    }
  total = 1./total;
  s->basefreqs[NUC_A] = total*freqa;
  s->basefreqs[NUC_C] = total*freqc;
  s->basefreqs[NUC_G] = total*freqg;
  s->basefreqs[NUC_T] = total*freqt;
  if (total*freqa < EPSILON)
    s->basefreqs[NUC_A] = EPSILON;
  if (total*freqc < EPSILON)
    s->basefreqs[NUC_C] = EPSILON;
  if (total*freqg < EPSILON)
    s->basefreqs[NUC_G] = EPSILON;
  if (total*freqt < EPSILON)
    s->basefreqs[NUC_T] = EPSILON;
  total = s->basefreqs[NUC_A] + s->basefreqs[NUC_C] + s->basefreqs[NUC_G] + s->basefreqs[NUC_T];
  s->basefreqs[NUC_A] /= total;
  s->basefreqs[NUC_C] /= total;
  s->basefreqs[NUC_G] /= total;
  s->basefreqs[NUC_T] /= total;
}

///
/// calculates (if necessary) and sets the base frequencies
void set_subloci_basefrequencies_seq(mutationmodel_fmt *s, world_fmt *world, option_fmt *options, data_fmt *data, long locus, long xs)
{
  (void) world;
  //long   top;
  //long   pop;
  //long   ind;
  //long   ii;
  //long   btotal=0;
  //MYREAL bsum=0.;
  //long locus = world->locus;
  MYREAL freqa, freqc, freqg, freqt, freqy, freqr;
  MYREAL freqgr = 0.0;
  MYREAL freqty = 0.0;
  MYREAL aa = 0.0;
  MYREAL bb = 0.0;

  switch(s->datatype)
    {
    case 'a':
      break;
    case 'b':
      break;
    case 'm':
      break;
    case 'n':
    case 's':
      if(s->basefreqs==NULL)
	s->basefreqs = (double*) mycalloc((s->numstates+BASEFREQLENGTH-4), sizeof(double));
#ifndef MPIDATAONDEMAND      
      if(options->freqsfrom)
	{
	  //empiricalfreqs (world, options, s, xs);
	  calc_seq_basefreq(s, options, data, locus, xs);
	}
      else
	{
	  check_basefreq(options);
	  s->basefreqs[NUC_A] = options->freqa;
	  s->basefreqs[NUC_C] = options->freqc;
	  s->basefreqs[NUC_G] = options->freqg;
	  s->basefreqs[NUC_T] = options->freqt;
	}
#else
      // handled with request_data 
#endif /*MPI_ONDEMAND*/
      freqa = s->basefreqs[NUC_A];
      freqc = s->basefreqs[NUC_C];
      freqg = s->basefreqs[NUC_G];
      freqt = s->basefreqs[NUC_T];
      
      if(s->datatype == 'n' || s->datatype == 'u')
	{
	  MYREAL sum=0.0;
	  if(freqa <= EPSILON)
	    {
	      freqa = BIGEPSILON;
	    }
	  sum += freqa;
	  if(freqc <= EPSILON)
	    {
	      freqc = BIGEPSILON;
	    }
	  sum += freqc;
	  if(freqg <= EPSILON)
	    {
	      freqg = BIGEPSILON;
	    }
	  sum += freqg;
	  if(freqt <= EPSILON)
	    {
	      freqt = BIGEPSILON;
	    }
	  sum += freqt;
	  if(sum>1.0)
	    {
	      freqa = s->basefreqs[NUC_A] = freqa/sum;
	      freqc = s->basefreqs[NUC_C] = freqc/sum;
	      freqg = s->basefreqs[NUC_G] = freqg/sum;
	      freqt = s->basefreqs[NUC_T] = freqt/sum;
	    }
	}
	s->basefreqs[NUC_R] = freqa + freqg;
	s->basefreqs[NUC_Y] = freqc + freqt;
	freqr = s->basefreqs[NUC_R];
	freqy = s->basefreqs[NUC_Y];
	s->basefreqs[NUC_AR]= freqa / freqr;
	s->basefreqs[NUC_CY]= freqc / freqy;
	s->basefreqs[NUC_GR]= freqg / freqr;
	s->basefreqs[NUC_TY]= freqt / freqy;
	//MYREAL freqar = s->basefreqs[NUC_AR];
	//MYREAL freqcy = s->basefreqs[NUC_CY];
	freqgr = s->basefreqs[NUC_GR];
	freqty = s->basefreqs[NUC_TY];
	aa = s->ttratio * (freqr) * (freqy) - freqa * freqg - freqc * freqt;
	bb = freqa * (freqgr) + freqc * (freqty);
	if(s->model==F84)
	  {
	    if(!s->finished)
	      {
		s->finished=TRUE;
		s->parameters[2] = 1.0;
		s->parameters[0] *=  s->basefreqs[NUC_Y]; 
		s->parameters[1] *=  s->basefreqs[NUC_R]; 
	      }
	  }
	if(s->model==HKY)
	  {
	    s->finished=TRUE;
	    //s->parameters[0] *= freqy;
	    //s->parameters[0] *= freqr;
	  }
	if(s->model==TN)
	  {
	    s->finished=TRUE;
	    //s->parameters[0] *= freqy;
	    //s->parameters[0] *= freqr;
	  }
	s->xi = aa / (aa + bb);
	s->xv = 1.0 - s->xi;
	if (s->xi <= 0.0)
	  {
	    //warning ("This transition/transversion ratio (%f)\n",s->ttratio);
	    //warning ("is impossible with these base frequencies (%f, %f, %f, %f)!\n",freqa,freqc,freqg,freqt);
	    s->xi = 0.00001; // do not set this to zero because of the 1/(fracchange=xi*(...))
	    s->xv = 0.99999;
	    s->ttratio = (freqa * freqg + freqc * freqt) / ((freqr) * (freqy));
	    
	    //	    warning (" Transition/transversion parameter reset\n");
	    //warning ("  so transition/transversion ratio is %10.6f\nIF this does not fit, use fixed nucleotide frequencies\n", (s->ttratio));
	  }
	// use 1/frac as precomputation speed up
	s->fracchange = 1. / ((s->xi) * (2. * freqa * (freqgr) + 2. * freqc * (freqty)) + 
			      (s->xv) * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg - freqt * freqt));
	break;
    default:
      break;
    }
}


void set_subloci_frequencies_alleles(world_fmt *world, option_fmt *options, data_fmt *data, long locus)
{
  long   sublocus;
  long   sublocistart = world->sublocistarts[locus];
  long   sublociend = world->sublocistarts[locus+1];
  //long   numpatterns;
  mutationmodel_fmt *s;
 
  for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
    {
      s =&world->mutationmodels[sublocus];
      //s->numstates = get_states(s, world);
      const long xs = sublocus; //- sublocistart;
      set_subloci_basefrequencies(s, world, options, data, locus, xs);
    }
}
void set_subloci_frequencies(world_fmt *world, option_fmt *options, data_fmt *data, long locus)
{
  long   sublocus;
  long   sublocistart = world->sublocistarts[locus];
  long   sublociend = world->sublocistarts[locus+1];
  //long   numpatterns;
  mutationmodel_fmt *s;
 
  for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
    {
      s =&world->mutationmodels[sublocus];
      s->numstates = get_states(s, data,locus);
      const long xs = sublocus;//- sublocistart;
      set_subloci_basefrequencies_seq(s, world, options, data, locus, xs);
    }
}


void set_siterates(long z, world_fmt *world, option_fmt *options)
{
  long i;
  mutationmodel_fmt *s = &world->mutationmodels[z];
  if(s->numsiterates==0)
    {
      s->numsiterates = options->rcategs;
      s->siterates = (double *) mycalloc(s->numsiterates,sizeof(double));
      s->siteprobs = (double *) mycalloc(s->numsiterates,sizeof(double));
    }
  for(i=0; i<s->numsiterates;i++)
    {
      s->siterates[i] = options->rrate[i];
      s->siteprobs[i] = options->probcat[i];
    }
}

void set_tn_model(double a1, double a2, double b, double pia, double pic, double pig, mutationmodel_fmt *mumod)
{
  // TN90 Tamura-Nei model
  //
  //new december 2009 {{0, -b, -ar - b, -ay - b}, {{1, 1, 1, 1}, {-(pY/pR), 1, -(pY/pR), 1}, {-(pG/pA), 0, 1, 0}, {0, -(pT/pC), 0, 1}}}
  // Transpose[evec] = {{1, -(pY/pR), -(pG/pA), 0}, {1, 1, 0, -(pT/pC)}, {1, -(pY/pR), 1, 0}, {1, 1, 0, 1}}
  //
  //
  // older
  //  double tneval[4] = { 0, -b, -a1 + a1*pia - b*pia + a1*pig - b*pig, -b - a2*pia + b*pia - a2*pig + b*pig }; 
  //double tnevec[4 * 4] = {
  // 1.0, 0.5/*1.0 - 1.0/(pia + pig)*/, 0.0,                        -1.0  /*-(pig/pia)*/,
  //1.0, 1.0,               -1.0/*(-1.0 + pia + pic + pig)/pic*/, 0.0,
  //1.0, 0.5 /*1.0 - 1/(pia + pig)*/,   0.0,                          1.0,
  //1.0, 1.0,                   1.0,                          0.0
  //};
  // F84: ar = ay        --> ar = ay = 0..inf,freq variable
  // HKY: ar/ay == pr/py --> ar = pr, ay = py, b=0..inf, freq variable
  // F81:                --> ar = ay = b = 1, freq variable
  // K2P:                --> ar = ay = 1, b= 0..inf, freq are 0.25
  // JC:                 --> ar = ay = b = 1, freq are 0.25
  double pit = 1.0 - pia - pic -pig;
  double pir = pia + pig;
  double piy = pic + pit;

  // eigenvalues
  mumod->eigenvalues[0] = 0.0;
  mumod->eigenvalues[1] = -b;
  mumod->eigenvalues[2] = -a1 - b;
  mumod->eigenvalues[3] = -a2 - b;
  // eigenvector
  mumod->eigenvectormatrix[0]= 1.0;
  mumod->eigenvectormatrix[1]= -piy / pir;
  mumod->eigenvectormatrix[2]= -(pig/pia);
  mumod->eigenvectormatrix[3]= 0.0;

  mumod->eigenvectormatrix[4]= 1.0;
  mumod->eigenvectormatrix[5]= 1.0;
  mumod->eigenvectormatrix[6]= 0.0;
  mumod->eigenvectormatrix[7]= -(pit/pic);

  mumod->eigenvectormatrix[8]= 1.0;
  mumod->eigenvectormatrix[9]= -piy/pir;  
  mumod->eigenvectormatrix[10]= 1.0;
  mumod->eigenvectormatrix[11]= 0.0;

  mumod->eigenvectormatrix[12]= 1.0;
  mumod->eigenvectormatrix[13]= 1.0;
  mumod->eigenvectormatrix[14]= 0.0;
  mumod->eigenvectormatrix[15]= 1.0;

  // inverse(eigenvector)
  //  {{pA, pC, pG, 1 - pA - pC - pG}, {-pA, -((pC (pA + pG))/(-1 + pA +pG)), -pG, ((pA + pG) (-1 + pA + pC + pG))/(-1 + pA + pG)}, 
  //   {-(pA/(pA + pG)), 0, pA/(pA + pG), 0}, {0, pC/(-1 + pA + pG), 0, -(pC/(-1 + pA + pG))}}
  //
  //  {
  //   {pA, pC, pG, pT}, 
  //   {-pA, pC pR/pY, -pG, pR pT/pY}, 
  //   {-pA/pR, 0, pA/pR, 0}, 
  //   {0, -pC/pY, 0, pC/pY}
  //  }
  mumod->inverseeigenvectormatrix[0]= pia;
  mumod->inverseeigenvectormatrix[1]= pic;
  mumod->inverseeigenvectormatrix[2]= pig;
  mumod->inverseeigenvectormatrix[3]= pit;

  mumod->inverseeigenvectormatrix[4]= -pia;
  mumod->inverseeigenvectormatrix[5]=  pic * pir/piy;  
  mumod->inverseeigenvectormatrix[6]= -pig;
  mumod->inverseeigenvectormatrix[7]=  pir * pit/piy;

  mumod->inverseeigenvectormatrix[8]= -pia/pir;
  mumod->inverseeigenvectormatrix[9]=  0.0;
  mumod->inverseeigenvectormatrix[10]= -pia/pir;
  mumod->inverseeigenvectormatrix[11]= 0.0;

  mumod->inverseeigenvectormatrix[12]= 0.0;
  mumod->inverseeigenvectormatrix[13]= -pic/piy;
  mumod->inverseeigenvectormatrix[14]= 0.0;
  mumod->inverseeigenvectormatrix[15]= pic/piy;
}


void set_mutationmodel_eigenmaterial(long z, world_fmt *world) //eigenvectormatrix, inverse eigenvector matrix, eigenvalues
{
  mutationmodel_fmt *mumod = &world->mutationmodels[z];
  double a1,a2,b,pia,pic,pig, piy, pir;
#ifdef DEBUG
  printf("%i> eigenmaterial: model=%i params:%f %f %f\n",myID,mumod->model,mumod->parameters[0],mumod->parameters[1],mumod->parameters[2]);
#endif
  // JC69 model eigenvector matrix
  //double jc69evec[4 * 4] = {
  //  1.0,  2.0,  0.0,  0.5,
  //  1.0,  -2.0,  0.5,  0.0,
  //  1.0,  2.0, 0.0,  -0.5,
  //  1.0,  -2.0,  -0.5,  0.0
  //};
  //double jc69evecPB[4 * 4] = {
  //  -1.0, -1.0, -1.0, 1.0,
  //  0.0, 0.0, 1.0, 1.0,
  //  0.0, 1.0, 0.0, 1.0,
  // 1.0, 0.0, 0.0, 1.0
  //};

  //double k2pevec[4 * 4] = {-1., 1., 0., -1., 1., 1., -1., 0., -1., 1., 0., 1., 1., 1., 1., 0};

  //beagle->inverseeigenvectormatrix = world->mutationmodel->inverseeigenvectormatrix;
  // JC69 model inverse eigenvector matrix
  //  double jc69ivec[4 * 4] = {
  //    0.25,  0.25,  0.25,  0.25,
  //    0.125,  -0.125,  0.125,  -0.125,
  //    0.0,  1.0,  0.0,  -1.0,
  //    1.0,  0.0,  -1.0,  0.0
  //  };
    //double jc69ivecPB[4 * 4] = {
    //  -0.25, -0.25, -0.25, 0.75,
    //  -0.25, -0.25, 0.75, -0.25,
    //  -0.25,  0.75, -0.25, -0.25, 
    //  0.25, 0.25, 0.25, 0.25
    //};

    //double k2pivec[4 * 4] = {-0.25, 0.25, 0., -0.5, 0.25, 0.25, -0.25, 0., -0.25, 0.25, 0., 0.5, 0.25, 0.25, 0.5, 0. };
    
    // JC69 model eigenvalues
    //double jc69eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };
    //double jc69evalPB[4] = {-1.0, -1.0, -1.0, 0.0};
    //K2P needs to be calculated because of kappa
    double kappa = 2.0;

    int named_model = mumod->model;
    const long numstates = mumod->numstates;
    const long numstates2 = numstates * numstates;
    mumod->eigenvalues = (double*) mycalloc(numstates,sizeof(double));
    mumod->eigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));
    mumod->inverseeigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));

    switch ((char) mumod->datatype)
      {
	// Models per datatype
	// Alleles
      case 'A':
      case 'a':
      case 'M':
      case 'm':
      case 'B':
      case 'b':
#ifdef BEAGLE
	error("not implemented in beagle yet!");
#endif
	break;
	// Msats
	// Msats Brownian
	// Sequence
      case 'S':
      case 's':      
      case 'N':
      case 'n':
      case 'u':
      case 'U':
      case 'F':
      case 'f':
	switch(named_model)//this should be through option!
	  {
	  case JC69:
	    //memcpy(mumod->eigenvalues,jc69eval,sizeof(double) * 4);
	    //memcpy(mumod->eigenvectormatrix,jc69evec,sizeof(double) * 16);
	    //memcpy(mumod->inverseeigenvectormatrix,jc69ivec,sizeof(double) * 16);
	    a1 = 1.0;
	    a2 = 1.0;
	    b  = 1.0;
	    pia = 0.25;
	    pic = 0.25;
	    pig = 0.25;
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case K2P:
	    kappa = mumod->ttratio;//mumod->parameters[0];
	    a1 = kappa;
	    a2 = kappa;
	    b  = 1.0;
	    pia = 0.25;
	    pic = 0.25;
	    pig = 0.25;
	    //mumod->eigenvalues[0] = -1.0;
	    //mumod->eigenvalues[1] = 0.0;
	    //mumod->eigenvalues[2] = mumod->eigenvalues[3] = -0.5 * (kappa + 1.0);
	    //memcpy(mumod->eigenvectormatrix,k2pevec,sizeof(double) * 16);
	    //memcpy(mumod->inverseeigenvectormatrix,k2pivec,sizeof(double) * 16);
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case F81:
	    a1 = 1.0;
	    a2 = 1.0;
	    b  = 1.0;
	    pia = mumod->basefreqs[0];
	    pic = mumod->basefreqs[1];
	    pig = mumod->basefreqs[2];
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case F84:
	    //	    kappa = mumod->parameters[0];
	    kappa = mumod->ttratio;
	    pia = mumod->basefreqs[0];
	    pic = mumod->basefreqs[1];
	    pig = mumod->basefreqs[2];
	    piy = pia + pig;//pic
	    pir = 1.0 - piy;
	    b  = 1.0;
	    a1 = (1.0 + kappa) * piy * b;
	    a2 = (1.0 + kappa) * pir * b;
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case HKY:
	    //	    kappa = mumod->parameters[0];
	    kappa = mumod->ttratio;
	    a1 = kappa;
	    a2 = kappa;
	    b  = 1.0;
	    pia = mumod->basefreqs[0];
	    pic = mumod->basefreqs[1];
	    pig = mumod->basefreqs[2];
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case TN:
	    a1  = mumod->parameters[0];
	    a2  = mumod->parameters[1];
	    b   = mumod->parameters[2];
	    pia = mumod->basefreqs[0];
	    pic = mumod->basefreqs[1];
	    pig = mumod->basefreqs[2];
	    set_tn_model(a1,a2,b,pia,pic,pig,mumod);
	    break;
	  case GTR:
	    error("not implemented yet");
	    //    break;
	  default:
	    break;
	}
      break;
      }
  // Sinlge nucleotide is the same as sequence but needs to calculate values for constant site patterns

}

///long read_word_delim(char *input, char *word, char * delim)
///
/// check whether there are #$ comments that contain  mutation model specification
void read_mutationmodel_comments(char *input, data_fmt * data, world_fmt *world)
{
  mutationmodel_fmt *s;
  char * modelname = (char *) mycalloc(LINESIZE,sizeof(char));
  char * word = (char *) mycalloc(LINESIZE,sizeof(char));
  long sublocus = 0;
  long readpos;
  readpos = read_word_delim(input, word, " ", TRUE); // sublocus
  sublocus = atoi(word)-1;
  readpos += read_word_delim(input+readpos, word, " ", TRUE);// scaled or not
  if (sublocus < data->allsubloci)
    {
      s = & world->mutationmodels[sublocus];
      s->from_infile = TRUE;
      if(word[0]=='s')
	{
	  s->scaling=TRUE;
	}
      else
	{
	  s->scaling=FALSE;
	}
      readpos += read_word_delim(input+readpos, modelname, " ", TRUE);//model
      //printf("word |%s|\n",modelname);
      // the placeholder # is giving the number of the sublocus starting with 1,2,3.....
      // labeling all subloci (or oldstyle loci) from first to last
      //$ # s F84: ar = ay        --> ar = ay = 0..inf,freq variable
      //$ # s HKY: ar/ay == pr/py --> ar = pr, ay = py, b=0..inf, freq variable
      //$ # s F81:                --> ar = ay = b = 1, freq variable
      //$ # s K2P:                --> ar = ay = 1, b= 0..inf, freq are 0.25
      //$ # s JC:                 --> ar = ay = b = 1, freq are 0.25
      //$ # s TN:                 --> ar, ay, b
      //$ # s SSM: 
      //$ # s MSM:                --> tune 0 .. 1.0, upchance 0 .. 2/3 
      //$ # s BM:
      //$ # s IAM:
      if(modelname[0]=='F') //F81 and F84
	{
	  if(modelname[2]=='1')
	    {
	      s->finished=TRUE;
	      s->model = F81;
	      s->parameters[0] = s->parameters[1] =  s->parameters[2] = 1.0;
	    }
	  else
	    {
	      s->model = F84;
	      s->finished=FALSE; //still needs empirical frequencies and NUC_Y and NUC_R
	      readpos += read_word_delim(input+readpos, word, " ", TRUE);
	      double kappa = atof(word);
	      s->ttratio = kappa;
	      s->parameters[2] = 1.0;//FIXME 
	      s->parameters[0] = (1.0 + kappa)  /* * s->basefreqs[NUC_Y]*/ * s->parameters[2]; 
	      s->parameters[1] = (1.0 + kappa)  /* * s->basefreqs[NUC_R]*/ * s->parameters[2]; 
	    }
	}
      if(modelname[0]=='K') //
	{
	  s->model= K2P;
	  s->finished=TRUE;
	  readpos += read_word_delim(input+readpos, word, " ", TRUE);
	  double kappa = atof(word);
	  s->ttratio = kappa;
	  s->parameters[2] = 1.0;
	  s->parameters[0] = kappa;
	  s->parameters[1] = 1.0;
	}
      if(modelname[0]=='H') //
	{
	  s->model= HKY;
	  s->finished=FALSE;
	  readpos += read_word_delim(input+readpos, word, " ", TRUE);
	  double kappa = atof(word);
	  s->ttratio = kappa;
	  s->parameters[2] = 1.0;
	  s->parameters[0] = kappa;
	  s->parameters[1] = kappa;
	  //s->parameters[0] = (-1.0 + kappa) * s->basefreqs[NUC_Y] * s->parameters[2]; 
	  //s->parameters[1] = (-1.0 + kappa) * s->basefreqs[NUC_R] * s->parameters[2]; 
	}
      if(modelname[0]=='T')
	{
	  s->finished=FALSE;
	  s->model = TN;
	  readpos += read_word_delim(input+readpos, word, " ", TRUE);
	  double kappa1 = atof(word);
	  readpos += read_word_delim(input+readpos, word, " ", TRUE);
	  double kappa2 = atof(word);
	  s->parameters[2] = (MYREAL) kappa2;
	  s->parameters[0] = (MYREAL) kappa1;
	  s->parameters[1] = 1.0;
	  
	}
      if(modelname[0]=='J')
	{
	  s->finished=TRUE;
	  s->model = JC69;
	  s->parameters[0] = s->parameters[1] =  s->parameters[2] = 1.0;
	}
      
      if(modelname[0]=='S' || modelname[0] == 'M') //Single/Multi-step model [microsatellite]
	{
	  s->finished=TRUE;
	  s->parameters[0] = 0.0;
	  s->parameters[2] = 0.0;
	  if(modelname[0]!='M')
	    {
	      s->model = SSM;
	    }
	  else
	    {
	      s->model = MSM;
	      readpos += read_word_delim(input+readpos, word, " ", TRUE);
	      if(word[0]!='\0')
		{
		  readpos += read_word_delim(input+readpos, word, " ", TRUE);
		  double tune = atof(word);
		  /*readpos +=*/ read_word_delim(input+readpos, word, " ", TRUE);
		  double upchance = atof(word);
		  if(tune >= 0.0 && tune <= 1.0)
		    s->parameters[0] = (MYREAL) tune;
		  
		  if(upchance >= 0.0 && upchance <= 0.666666)
		    s->parameters[2] = (MYREAL) upchance;
		}
	    }
	}
      if(modelname[0]=='B') //Brownian model [microsatellite]
	{
	  s->finished=TRUE;
	  s->model = BM;
	  s->parameters[0] = 0.0;
	  s->parameters[2] = 0.0;
	}
      if(modelname[0]=='I') //infinite allele model [microsatellite and other alleles]
	{
	  s->finished=TRUE;
	  s->model = IAM;
	  s->parameters[0] = 0.0;
	  s->parameters[2] = 0.0;
	}
    }
  myfree(word);
  myfree(modelname);
}


///
/// initialize the mutation model structure
void init_mutationmodel_readsites(mutationmodel_fmt *mumod, char datatype, char * sites) //long sites)
{
  long numsites;
  long startsite;
  char newdatatype;
  char *val;
  // sites must be preprocessed, it could be
  // 10 or  89s10 etc
  // the first one is the standard and has also a valid data type whereas the second
  // gives and explicit start site and hides the datatype, and get fixed here
  val=char_position("abmsnuh", sites);// function any char in sites???????
  if (val!=NULL)
    {
      numsites = atol(val+1);
      newdatatype=val[0];
      val = NULL;
      startsite = atol(sites); 

    }
  else
    {
      numsites = atol(sites);
      startsite = 0;
      newdatatype = datatype;
    }
  mumod->from_infile = FALSE;
  mumod->finished = FALSE;
  mumod->numpatterns = 0L; // number unique patterns
  mumod->numsites    = numsites; // number of sites
  mumod->startsite   = startsite; //begin of sequence in genome 0 is default
  mumod->numstates   = 0L; // number of states in model: DNA=4, DNA+gap=5, msat>2
  mumod->datatype    = newdatatype; //specifices the model
  mumod->dataclass   = (strchr("snuh",datatype)) ? SITECHARACTER : SITEWORD;
  mumod->scaling=TRUE;
  if(strchr(SEQUENCETYPES, mumod->datatype))
    mumod->model     = JC69; 
  else
    {
      switch(datatype)
	{
	case 'b':
	  mumod->model     = BM;
	  break;
	case 'm':
	  mumod->model     = SSM;
	  break;
	case 'a':
	  mumod->model     = IAM;
	  break;
	default:
	  mumod->model     = OTHER;
	}  
    }
  mumod->parameters  = (double *) mycalloc(NUMMUTATIONPARAMETERS,sizeof(double));
  mumod->lambda = 1.0;
  // init of site pattern material
  if (strchr("nuh",mumod->datatype))
    {
      mumod->addon = 4;
      numsites += mumod->addon;
    }
  if (mumod->datatype == 'u')
    numsites *= 5;
  mumod->alias = (long *) mycalloc (numsites, sizeof (long));
  mumod->ally = (long *) mycalloc (numsites, sizeof (long));
  mumod->savealiasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->aliasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->location = (long *) mycalloc (numsites, sizeof (long));
  mumod->category = (long *) mycalloc (numsites, sizeof (long));
  mumod->weight = (short *) mycalloc (numsites, sizeof (short));
  mumod->contribution = NULL;
  mumod->basefreqs = NULL;
  mumod->ttratio = 2.0;
}

void init_mutationmodel_readsites2(mutationmodel_fmt *mumod, char datatype, long numsites) //long sites)
{
  mumod->numpatterns = 0L; // number unique patterns
  mumod->numsites    = numsites; // number of sites
  mumod->numstates   = 0L; // number of states in model: DNA=4, DNA+gap=5, msat>2
  mumod->datatype    = datatype; //specifices the model
  mumod->dataclass   = (strchr("snuh",datatype)) ? SITECHARACTER : SITEWORD;
  mumod->scaling=TRUE;
  if(strchr(SEQUENCETYPES, mumod->datatype))
    mumod->model     = JC69; 
  else
    mumod->model     = OTHER;  
  mumod->parameters  = (double *) mycalloc(NUMMUTATIONPARAMETERS,sizeof(double));
  mumod->lambda = 1.0;
  // init of site pattern material
  if (strchr("nuh",mumod->datatype))
    {
      mumod->addon = 4;
      numsites += mumod->addon;
    }
  if (mumod->datatype == 'u')
    numsites *= 5;
  mumod->alias = (long *) mycalloc (numsites, sizeof (long));
  mumod->ally = (long *) mycalloc (numsites, sizeof (long));
  mumod->savealiasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->aliasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->location = (long *) mycalloc (numsites, sizeof (long));
  mumod->category = (long *) mycalloc (numsites, sizeof (long));
  mumod->weight = (short *) mycalloc (numsites, sizeof (short));
  mumod->contribution = NULL;
  mumod->basefreqs = NULL;
  mumod->ttratio = 2.0;
}

// used for data on demand
void init_mutationmodel_readsites3(mutationmodel_fmt *mumod, char datatype, long numsites) //long sites)
{
  (void) datatype;
  if (mumod->parameters ==NULL)
    mumod->parameters  = (double *) mycalloc(NUMMUTATIONPARAMETERS,sizeof(double));
  mumod->alias = (long *) mycalloc (numsites, sizeof (long));
  mumod->ally = (long *) mycalloc (numsites, sizeof (long));
  mumod->savealiasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->aliasweight = (double *) mycalloc (numsites, sizeof (double));
  mumod->location = (long *) mycalloc (numsites, sizeof (long));
  mumod->category = (long *) mycalloc (numsites, sizeof (long));
  mumod->weight = (short *) mycalloc (numsites, sizeof (short));
  mumod->contribution = NULL;
  mumod->basefreqs = NULL;
}


///
/// initialize the mutation model structure
void init_mutationmodel_first(world_fmt *world, data_fmt *data)
{
  const long numloci   = data->allsubloci;
  world->sublocistarts = (long *) mycalloc(numloci+1, sizeof(long));
  world->numsubloci = (long *) mycalloc(numloci, sizeof(long));
  world->maxnumpattern = (long *) mycalloc(numloci, sizeof(long));
#ifdef DEBUG
  printf("number of mutationmodels initialized: %li\n",numloci);
#endif
  world->mutationmodels = (mutationmodel_fmt *) mycalloc(numloci, sizeof(mutationmodel_fmt));
#ifdef MPI
  data->sublocistarts = world->sublocistarts;
  data->mutationmodels = world->mutationmodels;
#endif
}

///
/// initialize the mutation model structure
void init_mutationmodel_second(world_fmt *world, data_fmt *data)
{
  long totalmodels = 0L;
  long nummodels;
  long numloci   = data->loci;
  long i;
  for(i=0;i<numloci-1;i++)
    {
      nummodels = data->subloci[i];
      totalmodels += nummodels;
      world->sublocistarts[i+1] = totalmodels;
      world->numsubloci[i] = nummodels;
    } 
  nummodels = data->subloci[i];
  world->numsubloci[i] = nummodels;
  world->sublocistarts[i] = totalmodels;
  world->sublocistarts[i+1] = totalmodels + nummodels;


}

///
/// finish the mutation model structure
void finish_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options, long locus)
{
  //  const unsigned long locus = world->locus;
  long j;
  long z = 0L;
  long loci =  data->subloci[locus];
  //long nummaxalleles;
  long smallest, biggest;
  z = world->sublocistarts[locus];
  // for each sublocus
  //weird  for(j=0; j< loci; j++)
  for(j=0; j< loci; j++)
    {
	  mutationmodel_fmt *s = &world->mutationmodels[z];
	  //s->numpatterns = world->data->seq[0]->endsite; // number unique patterns
	  //s->numsites    = world->data->seq[0]->endsite; // number of sites
	  s->estimateseqerror = options->has_estimateseqerror;
	  s->seqerrorcombined = options->seqerrorcombined;
	  if (s->estimateseqerror)
	    {
	      s->scoring_error[0]=options->seqerror[0];
	      s->scoring_error[1]=options->seqerror[1];
	      s->scoring_error[2]=options->seqerror[2];
	      s->scoring_error[3]=options->seqerror[3];
	    }
	  if(world->maxnumpattern[locus] < s->numpatterns)
	    {
	      world->maxnumpattern[locus] = s->numpatterns;
	    }
	  if(data->oldsyntax)
	    {
	      s->datatype  = data->datatype[z]; //specifices the model
	    }
	  s->numstates = get_states(s, data, locus); // number of states in model: DNA=4, DNA+gap=5, msat>2
	  adjust_mutationmodel(s,options); // incase there is no infile definition then set options 
	  set_subloci_basefrequencies(s,world, options, data, locus, z); 
#ifdef DEBUG
	  printf("2@@@@@@@@@@@@@@@@@@@@@@@ %i> %p %i: 0:%f 1:%f 2:%f \n",myID, (void *) s, s->model, s->parameters[0],s->parameters[1],s->parameters[2]);
#endif	      

	  set_mutationmodel_eigenmaterial(z,world); //eigenvectormatrix, inverse eigenvector matrix, eigenvalues
	  set_siterates(z,world,options);

	  switch(s->datatype)
	    {
	    case 'a':
	    case 'b':
	    case 'm':
	      //nummaxalleles = findAllele(data,"\0",z);
	      //if(nummaxalleles == 0)
	      //nummaxalleles = get_states(s, data, locus); 
	      //s->numstates = nummaxalleles;
	      //s->maxalleles = nummaxalleles+1;
	      //s->freq = 1.0 / s->maxalleles;
	      break;
	    case 'x': // this needs to be fixed 
	      find_minmax_msat_allele (world, data, locus, &smallest, &biggest);
	      s->numstates  = get_states(s,data, locus); //biggest - smallest;
	      s->micro_threshold = options->micro_threshold;
	      if(s->micro_threshold > biggest-smallest)
		s->micro_threshold = biggest-smallest;//s->numstates;
	      if(s->micro_threshold % 2 != 0)
		s->micro_threshold += 1; // to make it even
	      s->microstart = smallest - 1 - s->micro_threshold;
	      if (s->microstart < 0)
		s->microstart = 0;
	      s->maxalleles = biggest - smallest + 2 * s->micro_threshold; //this is the scaled max
	      s->microrange = s->maxalleles - s->microstart;
	      s->basefreqs  = (double *) mycalloc((s->numstates+BASEFREQLENGTH-4),sizeof(double));
	      //printf("%i> micro-threshold=%li\n", myID, s->micro_threshold);
	      doublevec2d(&s->steps,s->micro_threshold,s->micro_threshold);
	      break;
	    default:
	      //s->numstates = get_states(s,data,locus);
	      if(s->basefreqs==NULL)
		s->basefreqs = (double *) mycalloc((s->numstates+BASEFREQLENGTH-4),sizeof(double));

	      break;
	    }
	  if(s->basefreqs==NULL)
	    s->basefreqs = (double *) mycalloc((s->numstates+BASEFREQLENGTH-4),sizeof(double));


	  z++;
    }
}

void destroy_mutationmodel(world_fmt* world)
{
  long sublocus;
  long locus;
  long i,j;
  if(world->sublocistarts != NULL)
    {      
      for(locus=0;locus<world->loci;locus++)
	{
	  long   sublocistart = world->sublocistarts[locus];
	  long   sublociend = world->sublocistarts[locus+1];
	  mutationmodel_fmt *s=NULL;      
	  for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
	    {
	      s =&world->mutationmodels[sublocus];
	      if (s!=NULL)
		{
		  myfree(s->parameters);
		  myfree(s->basefreqs);
		  myfree(s->siterates);
		  myfree(s->siteprobs);
		  myfree(s->rate);
		  myfree(s->eigenvectormatrix);
		  myfree(s->inverseeigenvectormatrix);
		  myfree(s->eigenvalues);
		  myfree(s->alias);
		  myfree(s->ally);
		  myfree(s->weight);
		  myfree(s->aliasweight);
		  myfree(s->location);
		  myfree(s->savealiasweight);
		  myfree(s->category);
		  myfree(s->contribution);
		  if(s->tbl!=NULL)
		    {
		      for(i = 0; i < s->numsiterates; i++)
			{
			  for (j = 0; j < world->options->categs; j++)	      
			    myfree(s->tbl[i][j]);
			}
		    }
		  if(s->tbl!=NULL)
		    myfree(s->tbl);
		}
	    }
	}
      myfree(world->mutationmodels);
    }
}


void print_mutationrate_weights(FILE *file, MYREAL *murates, long *segregs, MYREAL *wattersons, long loci)
{
  long locus;
  //  char ***buffer;
  //long rows = loci+2;
  double mumean=0.0;
  double segregmean=0.0;
  double wamean=0.0;
  fprintf(file,"Relative mutation rate among loci estimated from the data \n");
  fprintf(file,"==========================================================\n");
  if(wattersons==NULL)
    {
      fprintf(file,"Locus       Relative        Number of alleles\n");
      fprintf(file,"            mutation rate    \n");
      fprintf(file,"--------    -------------   -----------------\n");
      for (locus=0; locus < loci; locus++)
	{
	  fprintf(file,"%5li          % 5.5f         % 6li\n",locus+1, murates[locus], segregs[locus]);
	  mumean += (murates[locus] - mumean)/(locus+1);
	  segregmean += (segregs[locus] - segregmean)/(locus+1);
	}
      fprintf(file,"%5s          % 5.5f         % 6.1f\n","All", mumean, segregmean);
    }
  else
    {
      fprintf(file,"Locus       Relative        Watterson's Theta  Segregating\n");
      fprintf(file,"            mutation rate   (per site)         sites      \n");
      fprintf(file,"--------    -------------   -----------------  -----------\n");
      for (locus=0; locus < loci; locus++)
	{
	  fprintf(file,"%5li          % 5.5f         %2.8f      % 6li\n",locus+1, murates[locus], wattersons[locus], segregs[locus]);	
	  mumean +=  (murates[locus] - mumean)/(locus+1);
	  wamean +=  (wattersons[locus] - wamean)/(locus+1);
	  segregmean += (segregs[locus] - segregmean)/(locus+1);
	}
      fprintf(file,"%5s          % 5.5f         %2.8f      % 6.1f\n","All", mumean, wamean, segregmean);	
      fprintf(file,"\n\n\n");
    }
}



void copy_micro_steps(mutationmodel_fmt *s, mutationmodel_fmt *old)
{
  //printf("old->micro_threshold=%li\n",old->micro_threshold);
  if(s->steps==NULL)
    {
      doublevec2d(&s->steps,old->micro_threshold+1,old->micro_threshold+1);
      s->micro_threshold = old->micro_threshold;
    }
  else
    {
      free_doublevec2d(s->steps);
      doublevec2d(&s->steps,old->micro_threshold+1,old->micro_threshold+1);
      s->micro_threshold = old->micro_threshold;
    }
  long i;
  for(i=0;i<s->micro_threshold;i++)
    memcpy(s->steps[i],old->steps[i],sizeof(MYREAL) * (size_t) old->micro_threshold);
}

void klone_mutationmodel(world_fmt *newcopy, world_fmt *original, data_fmt *data, long locus)
{
  long z, mi;
  long numstates;
  long numstates2;
  size_t sites;

  if(newcopy->mutationmodels == NULL)
    init_mutationmodel_first(newcopy, data);

  memcpy(newcopy->sublocistarts, original->sublocistarts, sizeof(long) * (size_t) (data->loci+1));
  memcpy(newcopy->numsubloci, original->numsubloci, sizeof(long) * (size_t) data->loci);
  memcpy(newcopy->maxnumpattern, original->maxnumpattern, sizeof(long) * (size_t) data->loci);
  //z=newcopy->sublocistarts[locus];
  //i = locus;
  //  for(j=0; j< data->subloci[i]; j++)
  for(z=original->sublocistarts[locus]; z<original->sublocistarts[locus+1];z++)
    {
      mutationmodel_fmt *muold = &original->mutationmodels[z];
      mutationmodel_fmt *munew = &newcopy->mutationmodels[z];

      munew->estimateseqerror = muold->estimateseqerror;
      munew->seqerrorcombined = muold->seqerrorcombined;
      if (munew->estimateseqerror)
	{
	  munew->scoring_error[0]=muold->scoring_error[0];
	  munew->scoring_error[1]=muold->scoring_error[1];
	  munew->scoring_error[2]=muold->scoring_error[2];
	  munew->scoring_error[3]=muold->scoring_error[3];
	}
      
      
      munew->numpatterns = muold->numpatterns; // number unique patterns
      munew->numsites    = muold->numsites; // number of sites
      sites              = (size_t) muold->numsites + (strchr(SNPTYPES, muold->datatype) ? 4 : 0); //takes care of snp + totalsites
      munew->numstates   = muold->numstates; // number of states in model: DNA=4, DNA+gap=5, msat>2
      numstates          = munew->numstates;
      numstates2         = numstates * numstates;
      munew->numsiterates= muold->numsiterates;
      munew->lambda      = muold->lambda;

      


      munew->siterates   = (double *) mycalloc(munew->numsiterates,sizeof(double));
      munew->siteprobs   = (double *) mycalloc(munew->numsiterates,sizeof(double));
      memcpy(munew->siterates, muold->siterates, (size_t) muold->numsiterates * sizeof(double));
      memcpy(munew->siteprobs, muold->siteprobs, (size_t) muold->numsiterates * sizeof(double));

      munew->datatype        = muold->datatype; //specifices the model
      munew->dataclass       = muold->dataclass;
      munew->model           = muold->model;

      if(munew->parameters == NULL)
	munew->parameters  = (double *) mycalloc(NUMMUTATIONPARAMETERS,sizeof(double));
      memcpy(munew->parameters, muold->parameters,NUMMUTATIONPARAMETERS * sizeof(double));
      if(strchr(SEQUENCETYPES, muold->datatype))
	{
	  if(muold->basefreqs != NULL)
	    {
	      if(munew->basefreqs == NULL)
		munew->basefreqs = (double*) mycalloc((numstates+BASEFREQLENGTH-4), sizeof(double));
	      else
		munew->basefreqs = (double*) myrealloc(munew->basefreqs, (numstates+BASEFREQLENGTH-4) * sizeof(double));
	      memcpy(munew->basefreqs, muold->basefreqs, (size_t) (numstates+BASEFREQLENGTH-4) * sizeof(double));
	    }
	  if(muold->eigenvalues!=NULL)
	    {
	      if(munew->eigenvalues == NULL)
		{
		  munew->eigenvalues = (double*) mycalloc(numstates,sizeof(double));
		  munew->eigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));
		  munew->inverseeigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));
		}
	      else
		{
		  munew->eigenvalues = (double*) myrealloc(munew->eigenvalues, numstates * sizeof(double));
		  munew->eigenvectormatrix = (double*) myrealloc(munew->eigenvectormatrix, numstates2 * sizeof(double));
		  munew->inverseeigenvectormatrix = (double*) myrealloc(munew->inverseeigenvectormatrix, numstates2 * sizeof(double));
		}
	      memcpy(munew->eigenvalues, muold->eigenvalues, (size_t) numstates * sizeof(double));
	      memcpy(munew->eigenvectormatrix, muold->eigenvectormatrix, (size_t) numstates2 * sizeof(double));
	      memcpy(munew->inverseeigenvectormatrix, muold->inverseeigenvectormatrix, (size_t) numstates2 * sizeof(double));
	    }
	}
      if(muold->alias!=NULL)
	{
	  if(munew->alias==NULL)
	    {
	      munew->alias = (long *) mycalloc (sites, sizeof (long));
	      munew->ally = (long *) mycalloc (sites, sizeof (long));
	      munew->location = (long *) mycalloc (sites, sizeof (long));
	      munew->category = (long *) mycalloc (sites, sizeof (long));
	      munew->weight = (short *) mycalloc (sites, sizeof (short));
	      munew->savealiasweight = (double *) mycalloc (sites, sizeof (double));
	      munew->aliasweight = (double *) mycalloc (sites, sizeof (double));
	    }
	  memcpy(munew->alias, muold->alias, sizeof(long) * sites);
	  memcpy(munew->ally, muold->ally, sizeof(long) * sites);
	  memcpy(munew->location, muold->location, sizeof(long) * sites);
	  memcpy(munew->category, muold->category, sizeof(long) * sites);
	  memcpy(munew->weight, muold->weight, sizeof(short) * sites);
	  memcpy(munew->savealiasweight, muold->savealiasweight, sizeof(double) * sites);
	  memcpy(munew->aliasweight, muold->aliasweight, sizeof(double) * sites);
	  munew->rate = (MYREAL *) mycalloc (1, sizeof (MYREAL) * sites);
	}
      munew->contribution =  (contribarr *) mycalloc ((muold->numpatterns + muold->addon), sizeof (contribarr));
      munew->numcategs = muold->numcategs;
      munew->ttratio = muold->ttratio;
      munew->addon = muold->addon;
      munew->xi = muold->xi;
      munew->xv = muold->xv;
      munew->fracchange = muold->fracchange;

      munew->freq = muold->freq;
      munew->freqlast = muold->freqlast;
      munew->maxalleles = muold->maxalleles;

      munew->browniandefault = muold->browniandefault;

      munew->micro_threshold = muold->micro_threshold;
      munew->microrange = muold->microrange;
      munew->microstart = muold->microstart;
      if(!strchr(SEQUENCETYPES, muold->datatype))
	copy_micro_steps(munew,muold);

      for (mi=0;mi<4;mi++)
	munew->scoring_error[mi] = muold->scoring_error[mi];
      munew->scaling = muold->scaling;
      //z++;
    }
}

int get_mutationmodel(char x)
{
  switch(x)
    {
    case 'J' : return JC69;
      //break;
    case 'K' : return K2P;
      //break;
    case 'F' : return F84;
      //break;
    case 'D' : return F81;
      //break;
    case 'H' : return HKY;
      //break;
    case 'T' : return TN;
      //break;
    case 'M': return MSM;
      //break;
    case 'S': return SSM;
      //break;
    case 'B': return BM;
      //break;
    default:  
      return OTHER; 
      //break;
    }
}

void get_mutationmodel_nameparam(char *modelname, char *modelparams, mutationmodel_fmt *s)
{
  switch(s->model)
    {
    case JC69: strcpy(modelname,"Jukes-Cantor"); 
      sprintf(modelparams,"[Basefreq: =0.25]");
      break;
    case K2P: strcpy(modelname,"Kimura"); 
      sprintf(modelparams,"[Basefreq: =0.25, kappa=%.4f]",s->ttratio);
      break;
    case F81: strcpy(modelname,"Felsenstein 81"); 
      sprintf(modelparams,"[Basefreq:%.2f %.2f %.2f %.2f]",s->basefreqs[0],s->basefreqs[1],s->basefreqs[2], 1.0 - s->basefreqs[0] - s->basefreqs[1] - s->basefreqs[2]);
      break;
    case F84: strcpy(modelname,"Felsenstein 84"); 
      sprintf(modelparams,"[Bf:%.2f %.2f %.2f %.2f, t/t ratio=%.3f]",s->basefreqs[0],s->basefreqs[1],s->basefreqs[2], 1.0 - s->basefreqs[0] - s->basefreqs[1] - s->basefreqs[2], s->ttratio);
      break;
    case HKY: strcpy(modelname,"HKY"); 
      sprintf(modelparams,"[Bf:%.2f %.2f %.2f %.2f, kappa=%.3f]",s->basefreqs[0],s->basefreqs[1],s->basefreqs[2], 1.0 - s->basefreqs[0] - s->basefreqs[1] - s->basefreqs[2], s->ttratio);
      break;
    case TN: strcpy(modelname,"Tamura-Nei"); 
      sprintf(modelparams,"[Bf:%.2f %.2f %.2f %.2f, k1=%.3f, k2=%.3f]",s->basefreqs[0],s->basefreqs[1],s->basefreqs[2], 1.0 - s->basefreqs[0] - s->basefreqs[1] - s->basefreqs[2], s->parameters[0],s->parameters[2]);
      break;
    case GTR: strcpy(modelname,"GTR"); 
      sprintf(modelparams,"[not implemented]");
      break;
    case MSM: strcpy(modelname,"Multistep Model"); 
      sprintf(modelparams,"[tune=%.3f, upchance=%.3f]", s->parameters[0],s->parameters[2]);
      break;
    case SSM: strcpy(modelname,"Stepwise Model"); 
      sprintf(modelparams,"[none]");
      break;
    case BM: strcpy(modelname, "Brownian Motion"); 
      sprintf(modelparams,"[none]");
      break;
    case IAM: strcpy(modelname,"Infinite Allele"); 
      sprintf(modelparams,"[none]");
      break;
    default:  
      strcpy(modelname,"UNKNOWN"); sprintf(modelparams,"[none]"); break;
    }
}


void prob_tn93 (MYREAL p[4][4], const MYREAL u, const MYREAL ar, const MYREAL ay, const MYREAL b, const mutationmodel_fmt *s)
{
  /* mathematica derived transition probabilities for TN93 model*/
  const MYREAL pA = s->basefreqs[NUC_A];
  const MYREAL pC = s->basefreqs[NUC_C];
  const MYREAL pG = s->basefreqs[NUC_G];
  const MYREAL pT = s->basefreqs[NUC_T];
  const MYREAL pR = s->basefreqs[NUC_R];
  const MYREAL pY = s->basefreqs[NUC_Y];
  const MYREAL invpR = 1./pR;
  const MYREAL invpY = 1./pY;
  MYREAL x1 = ar * u;
  if (x1 > 100.0)
    x1 = 100.0;
  const MYREAL earu = exp(x1);
  const MYREAL invearu = 1./earu;
  MYREAL x2 = ay * u;
  if (x2 > 100.0)
    x2 = 100.0;
  const MYREAL eayu = exp(x2);
  const MYREAL inveayu = 1./eayu;
  MYREAL x3 = (ar+b) * u;
  if (x3 > 100.0)
    x3 = 100.0;
  const MYREAL earbu = exp(x3);
  const MYREAL invearbu = 1./earbu;
  MYREAL x4 = (ay+b) * u;
  if (x4 > 100.0)
    x4 = 100.0;
  const MYREAL eaybu = exp(x4);
  const MYREAL inveaybu = 1./eaybu;
  MYREAL x5 = b * u;
  if (x5 > 100.0)
    x5 = 100.0;
  const MYREAL ebu = exp(x5);
  const MYREAL invebu = 1./ebu;

  /* A to A  paa = E^(-(ar + b) * u) * (pA + pG)^(-1) (pG - E^(ar u) pA (-1 + pA + pG) + E^((ar + b) u) pA (pA + pG));*/
  p[NUC_A][NUC_A] =  invearbu * invpR * (pG + earu * pA * pY + earbu * pA * pR);
  /* A to C  pac = pC - E^(-b u) pC;*/
  p[NUC_A][NUC_C] = pC - invebu * pC;
  /* A to G  pag = E^(-b u) pG (pA + pG)^(-1) ( 1 - E^(-ar u) - pA - pG + E^(b u) (pA + pG)); */
  p[NUC_A][NUC_G] = invebu * pG * invpR * (pY - invearu + ebu * pR); 
  /* A to T  pat = -E^(-b u) (-1 + E^(b u)) (-1 + pA + pC + pG);*/
  p[NUC_A][NUC_T] = (1.0 - invebu)  * pT;
  /* C to A  pca = pA - E^(-b u) pA;*/
  p[NUC_C][NUC_A] = pA - invebu * pA;
  /* C to C  pcc = pC + E^(-(ay + b) u) (-1 + pA + pG)^(-1) (-1 + pA + pC + pG - E^(ay u) pC (pA + pG)); */
  p[NUC_C][NUC_C] = pC + inveaybu * invpY * (pT + eayu * pC * pR); 
  /* C to G  pcg = pG - E^(-b u) pG; */
  p[NUC_C][NUC_G] = pG - invebu* pG; 
  /* C to T  pct = E^(-(ay + b) u) (-1 + pA + pG)^(-1) (-1 + pA + pC + pG) (-1 - E^((ay + b) u) (-1 + pA + pG) + E^(ay u) (pA + pG));*/
  p[NUC_C][NUC_T] = inveaybu * invpY * pT * (eaybu * pY + eayu * pR - 1.0);
  /* G to A  pga = E^(-b u) pA (pA + pG)^(-1) (1 - E^(-ar u) - pA - pG + E^(b u) (pA + pG)); */
  p[NUC_G][NUC_A] = invebu * pA * invpR * (pY - invearu + ebu * pR); 
  /* G to C  pgc = pC - E^(-b u) pC; */
  p[NUC_G][NUC_C] = pC - invebu * pC; 
  /* G to G  pgg = E^(-b u) (pA + pG)^(-1) (E^(-ar u) pA + pG - pG (pA + pG) + E^(b u) pG (pA + pG));*/
  p[NUC_G][NUC_G] = invebu * invpR * (invearu * pA + pG - pG * pR + ebu * pG * pR); // * (1.0 - ebu);
  /* G to T  pgt = -E^(-b u) (-1 + E^(b u)) (-1 + pA + pC + pG);*/
  p[NUC_G][NUC_T] = (1.0 - invebu) * pT;
  /* T to A  pta = pA - E^(-b u) pA; */
  p[NUC_T][NUC_A] = pA - invebu * pA; 
  /* T to C  ptc = pC + E^(-(ay + b) u) (-1 + pA + pG)^(-1)  (pC - E^(ay u) pC (pA + pG)); */
  p[NUC_T][NUC_C] = pC + inveaybu * invpY * (eayu * pC * pR - pC); 
  /* T to G  ptg = pG - E^(-b u) pG; */
  p[NUC_T][NUC_G] = pG - invebu * pG; 
  /* T to T  ptt = E^(-b u) (-1 + pA + pG)^(-1) (-E^(-ay u) pC - E^( b u) (-1 + pA + pG) (-1 + pA + pC + pG) + (pA + pG) (-1 + pA + pC + pG))*/
  p[NUC_T][NUC_T] = invebu * invpY *(inveayu * pC + pT * (ebu * pY + pR));
}
	  
#ifndef AVX
// not using the AVX SIMD material
void multiply_add(MYREAL p[4][4], const MYREAL *xx1, MYREAL h[4], long numstates)
{
  (void) numstates;
  long i;
  MYREAL sum;
  MYREAL sum2;
  for(i=0;i<FOUR; i++)
    {
      MYREAL *pp = p[i];
      sum = pp[0] * xx1[0];
      sum2 = pp[1] * xx1[1];
      sum += pp[2] * xx1[2];
      sum2 += pp[3] * xx1[3];
      h[i] = sum + sum2;
      //  printf("%f ",h[i]);
    }
  //printf("\n");
}

#else

// get AVX intrinsics
#include <immintrin.h>
void multiply_add(double p[4][4], double *xx1, double h[4], long numstates)
{
  (void) numstates;
  int i;
  __m256d x;
  x = _mm256_load_pd(xx1);
  for(i=0;i<FOUR; i++)
    {
      __m256d pp;
      double * pm = (double *) (&p[i]);
      pp = _mm256_loadu_pd(pm);
      __m256d ht;
      ht = _mm256_mul_pd(x,pp);
      ht = _mm256_hadd_pd(ht,ht);
      h[i] = ((double*)&ht)[0] + ((double*)&ht)[2];
    }
}



///
/// calculates the conditional likelihood for one site using two daughters 
/// p1 and p2 are 4x4 matrices with transition probabilities between the 
/// nucleotides A C G T x A C G T
/// the pointer xx1 and xx2 contain the conditional likelihood for A C G T
/// for example a leave/tip in a tree with 1.0 0.0 0.0 0.0 has A as data
/// a G would be 0.0 0.0 1.0 0.0, 
/// xx3 is to pointer for the result (the mother node) 
void avx_like(double p1[4][4], double p2[4][4], double *xx1, double *xx2, 
	      double *xx3)
{
  // without the post of Norbert P. on
  // http://stackoverflow.com/questions/10833234/4-horizontal-double-precision-sums-in-one-go-with-avx?newreg=360a2ab30c06410199614718ceb9c969
  // I would have not been able to do this.
  __m256d x1,x2,x3;
  x1 = _mm256_loadu_pd(xx1);
  x2 = _mm256_loadu_pd(xx2);
  __m256d pp11,pp12,pp13,pp14,pp21,pp22,pp23,pp24;;
  pp11 = _mm256_loadu_pd((double *) (&p1[0])); 
  __m256d y0 = _mm256_mul_pd(x1,pp11);
  pp12 = _mm256_loadu_pd((double *) (&p1[1])); 
  __m256d y1 = _mm256_mul_pd(x1,pp12);
  __m256d sy01 = _mm256_hadd_pd(y0,y1);//a1,b1,a2,b2
  pp13 = _mm256_loadu_pd((double *) (&p1[2])); 
  __m256d y2 = _mm256_mul_pd(x1,pp13);
  pp14 = _mm256_loadu_pd((double *) (&p1[3])); 
  __m256d y3 = _mm256_mul_pd(x1,pp14);
  __m256d sy02 = _mm256_hadd_pd(y2,y3);//c1,d1,c2,d2
  //
  pp21 = _mm256_loadu_pd((double *) (&p2[0])); 
  __m256d z0 = _mm256_mul_pd(x2,pp21);
  pp22 = _mm256_loadu_pd((double *) (&p2[1])); 
  __m256d z1 = _mm256_mul_pd(x2,pp22);
  pp23 = _mm256_loadu_pd((double *) (&p2[2])); 
  __m256d z2 = _mm256_mul_pd(x2,pp23);
  pp24 = _mm256_loadu_pd((double *) (&p2[3])); 
  __m256d z3 = _mm256_mul_pd(x2,pp24);

  __m256d yblend = _mm256_blend_pd(sy01, sy02, 0b1100);//a1,b1,c2,d2
  __m256d yperm = _mm256_permute2f128_pd(sy01, sy02, 0x21);//a2,b2,c1,d1
  __m256d h1 =  _mm256_add_pd(yperm, yblend);//a,b,c,d

  __m256d sz01 = _mm256_hadd_pd(z0,z1);//a1,b1,a2,b2
  __m256d sz02 = _mm256_hadd_pd(z2,z3);//c1,d1,c2,d2
  __m256d zblend = _mm256_blend_pd(sz01, sz02, 0b1100);//a1,b1,c2,d2
  __m256d zperm = _mm256_permute2f128_pd(sz01, sz02, 0x21);//a2,b2,c1,d1
  __m256d h2 =  _mm256_add_pd(zperm, zblend);//a,b,c,d

  x3 = _mm256_mul_pd(h1,h2);
  // store into the standard pointer location
  _mm256_storeu_pd(xx3,x3);
}


#endif /*AVX multiply add and like calculators*/
	  
void multiply_add_general(MYREAL p[4][4], const MYREAL *xx1, MYREAL h[4], long numstates)
{
  long i;
  long z;
  for(i=0;i<numstates; i++)
    {
      MYREAL sum=0.0;
      MYREAL *pp = p[i];
      for(z=0; z < numstates; z++)
	{
	  sum += pp[z] * xx1[z];
	}
      h[i] = sum;
    }
}

/*
def nux(offspring, rate,  nuc1,  nuc2, numoffspring, numsiterates, numnuc):
  return (offspring * numsiterates * numnuc * numnuc + numnuc * numnuc * rate  + nuc1*numnuc+nuc2)
# test nux
numoffspring=2
numsiterates=3
numnuc=4
print(f"allocated: {numoffspring*numsiterates*numnuc*numnuc}")
for offspring in range(numoffspring):
    print(f"line {offspring}:",end='\n')
    for rate in range(numsiterates):
        print(f"rate {rate}:",end=' ')
        for nuc1 in range(numnuc):
            for nuc2 in range(numnuc):
                print(nux(offspring,rate,nuc1,nuc2,numoffspring, numsiterates,numnuc),end=' ')
        print()
    print()
print("end")
 */

long nux(long offspring,long rate, long nuc1, long nuc2, long numsiterates,long numnuc)
{
  //nuc2 + nuc1 numnuc + numnuc^2 numsiterates offspring + numnuc^2 rate
  long numnuc2 = numnuc * numnuc;
  return (offspring * numsiterates * numnuc2 + numnuc2 * rate  + nuc1*numnuc+nuc2);
}


void nuview_tn93 (mutationmodel_fmt *s, long sublocus, long xs, node * mother, world_fmt * world, const long locus)
{
  (void) sublocus;
  (void) world;
  (void) locus;
  //  static long count=0;
  //MYREAL p1[s->numsiterates][4][4];
  //MYREAL p2[s->numsiterates][4][4];
  MYREAL p1[9][4][4];
  MYREAL p2[9][4][4];
  node *d1;
  node *d2;
  children (mother, &d1, &d2);
  const MYREAL u1 = d1->v;
  const MYREAL u2 = d2->v;
  const long numpatterns = s->numpatterns + s->addon;
  MYREAL *s1 = d1->scale[xs];
  MYREAL *s2 = d2->scale[xs];
  MYREAL * xx3;
  MYREAL sxx3m = -MYREAL_MAX;
  MYREAL invsxx3m;
  //MYREAL *s3 = mother->scale[xs];

  long site;
  long rate;
  sitelike * xx1;
  sitelike * xx2;
  long nuc;
  MYREAL tempsxx3m;
  const long numsiterates = s->numsiterates;
  for(rate=0; rate < numsiterates; rate++)
    {
      prob_tn93(p1[rate], s->siterates[rate] * u1, s->parameters[0], s->parameters[1], s->parameters[2], s);
      prob_tn93(p2[rate], s->siterates[rate] * u2, s->parameters[0], s->parameters[1], s->parameters[2], s);
    }
#ifndef AVX
  MYREAL h1[4][4];
  MYREAL h2[4][4];
  const long numstates    = s->numstates;
#else
  const long numstates    = s->numstates;
#endif
  const phenotype * x1 = &(d1->x[xs].s);
  const phenotype * x2 = &(d2->x[xs].s);
  //phenotype * x3 = &(mother->x[xs].s);
  for(site=0; site < numpatterns; site++)
    {
      for(rate=0; rate < numsiterates; rate++)
	{
	  xx1 = &(*x1)[site][rate];
	  xx2 = &(*x2)[site][rate];
	  //sitelike * xx3 = &(*x3)[site][rate];
#ifdef AVX
	  avx_like(p1[rate],p2[rate],*xx1,*xx2,mother->x[xs].s[site][rate]);
#else
	  multiply_add(p1[rate], *xx1, h1[rate], numstates); 
	  multiply_add(p2[rate], *xx2, h2[rate], numstates);
	  for(nuc=0; nuc < numstates; nuc++)
	    {	     
	      mother->x[xs].s[site][rate][nuc] = h1[rate][nuc] * h2[rate][nuc];
	    }
#endif
	  MYREAL ss=0.0;
	  for(nuc=0; nuc < numstates; nuc++)
	    {	     
	      ss += mother->x[xs].s[site][rate][nuc];
	    }
	  if (ss<=0.0)
	    {
	      for(nuc=0; nuc < numstates; nuc++)
		{	     
		  mother->x[xs].s[site][rate][nuc] = VERYSMALL_VALUE;
		}	      
	    }
	}
      mother->scale[xs][site] = s1[site] + s2[site];	
    }
  if(s->scaling)
    {
	  for (site = 0; site < numpatterns; site++)
	    {
	      sxx3m = -MYREAL_MAX;
	      for (rate = 0; rate < numsiterates; rate++)
		{
		  xx3 = mother->x[xs].s[site][rate];
		  //printf("@%@ site=%li rate=%li: %f %f %f %f\n",site, rate, xx3[0],xx3[1],xx3[2],xx3[3]);
		  tempsxx3m = MAX ((xx3)[0], (xx3)[1]);
		  tempsxx3m = MAX (tempsxx3m, (xx3)[2]);
		  tempsxx3m = MAX (tempsxx3m, (xx3)[3]);
		  if (tempsxx3m > sxx3m)
		    sxx3m = tempsxx3m;
		}
	      if(sxx3m < SMALLEPSILON)
		sxx3m = SMALLEPSILON;
	      invsxx3m = 1. / sxx3m;
	      //MYREAL *mmms = mother->x[xs].s[site][rate];
	      for (rate = 0; rate < numsiterates; rate++)
		{
#ifdef AVX
		  __m256d mom = _mm256_loadu_pd(mother->x[xs].s[site][rate]);
		  __m256d iii = _mm256_broadcast_sd(&invsxx3m);
		  _mm256_storeu_pd(mother->x[xs].s[site][rate],_mm256_mul_pd(iii,mom));
#else
		  mother->x[xs].s[site][rate][0] *=  invsxx3m;
		  mother->x[xs].s[site][rate][1] *=  invsxx3m;
		  mother->x[xs].s[site][rate][2] *=  invsxx3m;
		  mother->x[xs].s[site][rate][3] *=  invsxx3m;
#endif
		}
	      //	      printf("scaler:%f ",log(sxx3m));
	      mother->scale[xs][site] += LOG (sxx3m);
	      //printf("%f\n", mother->scale[xs][site]);
	    }
	  //    }
    }
}


void force_basefreqs(MYREAL ** basefreqs, MYREAL pA, MYREAL pC, MYREAL pG)
{

  MYREAL pT = 1.0 - pA - pC - pG;
  MYREAL qr;
  MYREAL qy;
  if(*basefreqs==NULL)
    *basefreqs = (double*) calloc(BASEFREQLENGTH, sizeof(double));

  (*basefreqs)[NUC_A] = pA;
  (*basefreqs)[NUC_C] = pC;
  (*basefreqs)[NUC_G] = pG;
  (*basefreqs)[NUC_T] = pT;
  (*basefreqs)[NUC_R] = pA + pG;
  (*basefreqs)[NUC_Y] = pC + pT;
  qr = (*basefreqs)[NUC_R];
  qy = (*basefreqs)[NUC_Y];
  (*basefreqs)[NUC_AR]= pA / qr;
  (*basefreqs)[NUC_CY]= pC / qy;
  (*basefreqs)[NUC_GR]= pG / qr;
  (*basefreqs)[NUC_TY]= pT / qy;
  //basefreq is used somewhere else
}

void pseudonu_tn93 (mutationmodel_fmt *s, proposal_fmt *proposal, xarray_fmt *xxx1, MYREAL *lx1, MYREAL v1, xarray_fmt *xxx2, MYREAL *lx2, MYREAL v2, long xs)
{
  (void) proposal;
  //static long count=0;
#ifdef WINDOWS
  MYREAL p1[10][4][4];
  MYREAL p2[10][4][4];
#else
  // clang with -Weverything complains about vla variable length array
  //MYREAL p1[s->numsiterates][4][4];
  //MYREAL p2[s->numsiterates][4][4];
  MYREAL p1[10][4][4];
  MYREAL p2[10][4][4];
#endif
  const MYREAL u1 = v1;
  const MYREAL u2 = v2;
  const long numpatterns = s->numpatterns + s->addon;
  MYREAL *s1 = lx1;
  MYREAL *s2 = lx2;
  MYREAL *s3 = lx1;

  long site;
  long rate;
  const long numsiterates = s->numsiterates;
  const double ar = s->parameters[0];
  const double ay = s->parameters[1];
  const double b = s->parameters[2];
  for(rate=0; rate < numsiterates; rate++)
    {
      prob_tn93(p1[rate], s->siterates[rate] * u1, ar, ay, b, s);
      prob_tn93(p2[rate], s->siterates[rate] * u2, ar, ay, b, s);
    }
#ifndef AVX
  MYREAL h1[4];
  MYREAL h2[4];
  const long numstates    = s->numstates;
#endif
  phenotype * x1 = &xxx1[xs].s;
  phenotype * x2 = &xxx2[xs].s;
  phenotype * x3 = &xxx1[xs].s;
  for(site=0; site < numpatterns; site++)
    {
      for(rate=0; rate < numsiterates; rate++)
	{
	  sitelike * xx1 = &(*x1)[site][rate];
	  sitelike * xx3 = &(*x3)[site][rate];
	  sitelike * xx2 = &(*x2)[site][rate];
#ifdef AVX
	  //__assume_aligned(*xx3, 32); 
	  avx_like(p1[rate],p2[rate],*xx1,*xx2, *xx3);
#else
	  multiply_add(p1[rate], *xx1, h1, numstates); 
	  multiply_add(p2[rate], *xx2, h2, numstates);
	  long nuc;
	  for(nuc=0; nuc < numstates; nuc++)
	    {
	      //	      if(h1[nuc] <= 0.0)
	      //error("f.... in [pseudo]nuview_tn93");
	      //if(h2[nuc] <= 0.0)
	      //	error("f.... in [pseudo]nuview_tn93");
	      (*xx3)[nuc] = h1[nuc] * h2[nuc];
#ifdef DEBUGXXXX
	      if((*xx3)[nuc] - 1.0 > EPSILON)
		{
		  printf("cond like strangely high in %s line  %i (value=%f) \n",__FILE__,__LINE__, (*xx3)[nuc]);
		  //sleep(10);
		}
#endif
	    }
#endif /*AVX*/
	}
      s3[site] = s1[site] + s2[site];
    }
  if(s->scaling)
    {
      //count++;
      //if (count == SCALEINTERVAL)
      //if(count>0)
      //	{
      //	  count = 0;
	  for (site = 0; site < numpatterns; site++)
	    {
	      MYREAL sxx3m = -MYREAL_MAX;
	      for (rate = 0; rate < numsiterates; rate++)
		{
		  sitelike * xx1 = &(*x1)[site][rate];
		  //xx1 = &(xxxx1[i][j]);
		  MYREAL tempsxx3m = MAX ((*xx1)[0], (*xx1)[1]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[2]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[3]);
		  if (tempsxx3m > sxx3m)
		    sxx3m = tempsxx3m;
		}
	      MYREAL invsxx3m = 1. / sxx3m;
	      for (rate = 0; rate < numsiterates; rate++)
		{
		  sitelike * xx1 = &(*x1)[site][rate];
#ifdef AVX
		  __m256d mom = _mm256_loadu_pd(*xx1);
		  __m256d iii = _mm256_broadcast_sd(&invsxx3m);
		  _mm256_storeu_pd(*xx1,_mm256_mul_pd(iii,mom));
#else
		  (*xx1)[0] *= invsxx3m;
		  (*xx1)[3] *= invsxx3m;
		  (*xx1)[1] *= invsxx3m;
		  (*xx1)[2] *= invsxx3m;

#endif
		}
	      s3[site] += LOG (sxx3m);
	    }
    }
}

