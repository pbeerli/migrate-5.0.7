/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 H A P L O T Y P E   R O U T I N E S 
 
 creates and manipulates data structures for haplotyping.
 
 Peter Beerli
 beerli@fsu.edu

the main function in this file is:

 replaced with the function data.c: 
 read_comments(): read_individuals_dictionary(data, options)
 insert_individual_inlist(...)


Copyright 2010 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
*/
#include "migration.h"
#include "random.h"
#include "world.h"
#include "tree.h"
#include "mcmc.h"
#include "tools.h"
#include "sighandler.h"
#include "hash.h"
#include "pretty.h"
#ifdef MPI
#include "migrate_mpi.h"
#endif
#include "haplotype.h"
extern int myID;

//long calc_individual_checksum(int *difference, world_fmt * world, long locus);
long find_inDB(char * thisname, long region, individualDB_fmt *individuals, long numindividuals);
MYREAL extract_showkeys(char **showkeys, int *showvalues, int leaves, int numstates, int i);
long find_in_haplotypingorderDB(char * thisname, haplotyping_orderform_fmt *individuals, long numindividuals);
boolean is_same_sequence(phenotype * xx1, phenotype * xx2, long numpatterns, long numsiterates, long numstates);
boolean is_same_allele(MYREAL ** xx1, MYREAL ** xx2, long numpatterns, long numsiterates, long numstates);
void check_individual_node(char *rawname, node *thenode, long region, world_fmt *world);
void debug_states(individualDB_fmt *w);

// reporting the state changes for each indiviudals
// by counting the number of different states explored for each individual, 
// this could become way too large and may need to be turned off otherwise the system will run out memory
// in these cases perhaps simply printing the string out to a file may help, although that may result in huge files.

void insert_individual_inlist(char *input, data_fmt *data, option_fmt *options)
{
  // read name
  // read subloci model: use 0 for no haplotyping and * for haplotyping
  // an example  looks like this: #% rid1324: *0*00*** 0 0 
  // the name is delimited with a colon, there are 8 subloci in locus 1
  // the second, fourth and fifth are not haplotyped and 2 other loci that are not
  // haplotyped, regions are separated by a space, if there is no data for a
  // a particular locus then use a - (minus sign)
  // for example:  #% rid1324: *0*00-** - 0
  // important: #% name forces that there are 3 characters before the name starts.
  if(data->haplotyping_list == NULL)
    {
      data->haplotyping_list = (haplotyping_orderform_fmt *) mycalloc((data->numhaplotyping+1), sizeof(haplotyping_orderform_fmt));
      data->haplotyping_list[0].key  = (char *) calloc(LINESIZE, sizeof(char));
      charvec2d(&data->haplotyping_list[0].pick, data->loci, LINESIZE);
    }
  else
    {
      data->haplotyping_list = (haplotyping_orderform_fmt *) realloc(data->haplotyping_list,sizeof(haplotyping_orderform_fmt)* (size_t) (data->numhaplotyping+1));
      data->haplotyping_list[data->numhaplotyping].key  = (char *) calloc(LINESIZE, sizeof(char));
      charvec2d(&data->haplotyping_list[data->numhaplotyping].pick, data->loci, LINESIZE);
    }
  char *tmp;
  char *input2 = (char *) strdup(input);
  tmp = strtok (input2, ":\n");
  trim(&tmp);
  if (tmp != NULL)
    {
      int i;
      size_t tmplen = (size_t) strlen(tmp);
      for (i = 0; i < data->numhaplotyping; i++)
	{
	  if(tmplen == (size_t) strlen(data->haplotyping_list[i].key))
	    {
	      if(strstr(tmp,data->haplotyping_list[i].key)!=NULL)
		break;
	    }
	    else
	    {
	      if(strstr(data->haplotyping_list[i].key, tmp)!=NULL)
		break;
	    }
	}
      if(i<data->numhaplotyping)
	error("Identical individual is twice in the haplotyping list");
      strncpy (data->haplotyping_list[data->numhaplotyping].key,tmp, (size_t) options->nmlength);
      for(i=0;i<data->loci;i++)
	{
	  tmp = strtok (NULL, " ");
	  if(tmp == NULL)
	    continue;
	  tmplen = (size_t) strlen(tmp);
	  strncpy (data->haplotyping_list[data->numhaplotyping].pick[i],tmp, tmplen);
	}
    }
  else
    {
      strncpy (data->haplotyping_list[data->numhaplotyping].key,"NONAME", (size_t) options->nmlength);
    }
  //printf("%s %s\n",data->haplotyping_list[data->numhaplotyping].key,data->haplotyping_list[data->numhaplotyping].pick[0]);
  data->numhaplotyping += 1;
  myfree(input2);
}

/*
void read_individuals_dictionary(char **input, data_fmt * data, option_fmt *options)
{
  while(input[0]=='#')
    {
      switch(input[1])
	{
	case '%':
	  FGETS(*input, LINESIZE, data->infile);
	  insert_individual_inlist(*input, data, options);
	  break;
	case '$': // mutation model -- jump out
	  unread_word(data->infile, *input);
	  break;
	default:      // this is a plain comment
	  FGETS(*input, LINESIZE, data->infile);
	  break;
	}
      read_word(data->infile, *input, NULL);
    }
  unread_word(data->infile, input);
  myfree(input);
}
*/


void  set_individuals_request_haplotyping(world_fmt * world, data_fmt * data, long locus)
{
  individualDB_fmt *ind = world->data->individuals;
  const long last = data->numindividuals[locus];
  long id;
  long i;
  if(data->numhaplotyping>0)
    {
      for(i=0; i<data->numhaplotyping; i++)
	{      
	  id = find_inDB(data->haplotyping_list[i].key, locus, ind, last);
	  if(id != -1)
	    {
	      world->data->individuals[id].targeted=i+1; // plus one because zero is valid
#ifdef DEBUG
	      printf("%s === %s\n",data->haplotyping_list[i].key,  world->data->individuals[id].name);
#endif
	    }
	}
      for(id=0; id<last; id++) 
	{
	  if(world->data->individuals[id].targeted>0)
	    {
	      long l;
	      for(l=0;l<data->subloci[locus];l++)
		{
		  if(data->haplotyping_list[world->data->individuals[id].targeted-1].pick[locus][l]==0)
		    world->data->individuals[id].difference[l] = 0;
		  else
		    world->data->individuals[id].difference[l] = 1;
		}
	      world->data->individuals[id].checksum = calc_individual_checksum(world->data->individuals[id].difference, world, locus);
	    }
	  else
	    {
	      world->data->individuals[id].checksum=0;
	    }
	}
    }
  else
    {
      if(world->data->haplotyping)
	{
	  for(id=0;id<last;id++)
	    {
	      world->data->individuals[id].targeted=world->data->individuals[id].checksum;
	    }
	}
    }
}

/// add a count to a particular hapltype configuration
/// assumes that a key os a string of 01011010 etc
void add_state_counter(hash_fmt **hash, size_t *numhash, float *total, int *states, long numstates, int count)
{
  char *key;
  long i;

  key = (char *) mycalloc((numstates + 1),sizeof(char));

  for(i=0; i<numstates; i++)
    {
      key[i] =  (char) (48+states[i]);
    } 
  key[numstates]='\0';
  int offset = 0;
  insert_hash(key, &offset, (size_t) numstates, (size_t) numstates, count, hash, numhash);
  (*total) += (float) count;
  myfree(key);
}

///
/// combines v alues from complement and destroys content of showkeys_i and _j
/// !!!!!!!!!! pay attention: this destroys showkeys !!!!!!!!!!!!!!!!!!!!!!!!!
MYREAL extract_showkeys(char **showkeys, int *showvalues, int leaves, int numstates, int i)
{
  char *key = showkeys[i];
  MYREAL value = (MYREAL) showvalues[i];
  char * transposed_key = (char *) mycalloc((numstates+1),sizeof(char));
  long j;
  for(j=0;j<numstates;j++)
    {
      transposed_key[j] = key[j]=='0' ? '1' : '0';
    }

  for(j=0; j<leaves; j++)
    {
      if(strncmp(transposed_key, showkeys[j], (size_t) numstates)==0)
	{
	  value += (MYREAL) showvalues[j];
	  showkeys[j][0]='\0';
	  myfree(transposed_key);
	  return value;
	}
    }
  myfree(transposed_key);
  return value;
}


//
// print haplotype states
void print_haplotypes2(char **buffer, long *allocbufsize, long *bufsize, world_fmt *world, long locus, long pop, long ind, boolean fraction)
{
  //long numreplicates = (world->options->replicatenum > 0 ? world->options->replicatenum : 1);
  individualDB_fmt * ii = &world->haplotypes[locus][ind];
  //  printf("%li>    %li  %li:%-10.10s\n", locus, ind, pop+1, ii->name);
  long namelength = *bufsize;
  int poplen = (world->numpop > 8 ? ((world->numpop > 98) ? 3 : 2) : 1);
  increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
  *bufsize += sprintf(*buffer + *bufsize,"%*li:%-10.10s\t ", poplen,pop+1, ii->name);
  namelength = *bufsize - namelength + 1;
  if(ii->total1 != 0.0f)
    {
      char **showkeys, **showkeys2;
      long leaves=0;
      hash_leaves(ii->hash,ii->numhash, &leaves);
      charvec2d(&showkeys,leaves,ii->numstates+1);
      charvec2d(&showkeys2,leaves,ii->numstates+1);
      int *showvalues = (int *) mycalloc(leaves,sizeof(int));
      long z=0;
      int i;
      long lines = *bufsize + HUNDRED - namelength;
      show_hash(ii->hash,ii->numhash, (size_t) ii->numstates, &showkeys,&showvalues,&z);
      long count = 0;
      MYREAL * values = (MYREAL*) mycalloc(leaves,sizeof(MYREAL));
      MYREAL total=0.0;
      for (i=0; i < leaves; i++)
	{
	  if(showkeys[i][0]!='\0')
	    {
	      MYREAL value = extract_showkeys(showkeys,showvalues, (int) leaves, (int) ii->numstates, i);///(MYREAL) numreplicates;
	      strncpy(showkeys2[i],showkeys[i],(size_t) ii->numstates);
	      values[count] = value;
	      count++;
	      total += value;
	    }
	}
      count=0;
      for (i=0; i < leaves; i++)
	{
	  if(showkeys[i][0]!='\0')
	    {
	      MYREAL value = values[count];
	      count++;
	      if(fraction)
		{
		  if((ii->numstates+7+ *bufsize) > lines)
		    {
		      lines += HUNDRED - namelength;;
		      increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
		      *bufsize += sprintf(*buffer + *bufsize,"\n\t ");
		      //		  *bufsize += sprintf(*buffer + *bufsize,"\n%-*.*s",namelength, namelength," ");
		    }
		}
	      increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
	      *bufsize += sprintf(*buffer + *bufsize,"%s:%1.3f ",showkeys2[i], fraction ? value/total : value);
	    }
	}
      if(count==0 || !fraction)
	{
	  // no haplotype tried and recorded
	  increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
	  *bufsize += sprintf(*buffer + *bufsize,"\n");
	}
      else
	{
	  increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
	  *bufsize += sprintf(*buffer + *bufsize,"   (%li)\n", (long) total);
	}
      free_charvec2d(showkeys);
      free_charvec2d(showkeys2);
      myfree(showvalues);
    }
  else
    {
      increase_buffer(buffer, allocbufsize, bufsize,LINESIZE);
      *bufsize += sprintf(*buffer + *bufsize,"\n");
    }
}

///
/// print haplotypes() this is a rewrite of the report_states to take into account, replicates distributed loci etc.
/// that will be aggregated to the master,cold node through some other functions ath the end of a replicate or locus.
/// this transfers will be modeled similarly to the histogram transfers [methinks?!not done yet]
void print_haplotypes(char **buffer, long *allocbufsize, long *bufsize, world_fmt* world, data_fmt *data)
{
  //worlddata_fmt *wdata = world->data;
  long ind;
  long locus;
  long pop;
  char *oldname;
  char *newname;
  oldname = (char *) calloc(SMALLBUFSIZE,sizeof(char));//this is 255 char
  newname = (char *) calloc(SMALLBUFSIZE,sizeof(char));//this is 255 char
  *bufsize += sprintf(*buffer+ *bufsize,  "Haplotype states and probabilities\n----------------------------------\n");
  pdf_print_haplotype_title();
  long pdfstartbuf=*bufsize;
  for(locus=0;locus<world->loci;locus++)
    {
      *bufsize += sprintf(*buffer+ *bufsize,  "Locus: %li\n",locus+1);
      for(pop=0;pop<world->numpop;pop++)
	{
	  for(ind=0; ind<world->data->numind[pop][locus]; ind++)	   
	    {
	      // assumes order of names, and second occurrence in tree will
	      // not be used, this seems like a hack but otherwise I nee to extract
	      // or store the pop assignment in the name.
	      size_t numchar = (size_t) strcspn(data->indnames[pop][ind][locus],":");
	      strncpy(newname,data->indnames[pop][ind][locus], numchar);
	      newname[numchar]='\0';
	      size_t ll = (size_t) strlen(newname);
	      if(!strncmp(oldname,newname,  ll))
		continue;
	      else
		strncpy(oldname,newname,ll);
		       
	      long id = find_inDB(data->indnames[pop][ind][locus], locus, world->haplotypes[locus], world->numhaplotypes[locus]);
	      if( id != -1)
		{
		  print_haplotypes2(buffer, allocbufsize, bufsize, world, locus, pop, id,TRUE);
		}
	    }
	}
    }
  pdf_print_haplotypes2(*buffer+pdfstartbuf, *bufsize - pdfstartbuf);
  myfree(oldname);
  myfree(newname);
}



void save_haplotypes(world_fmt* world, long locus)
{
  if(world->cold)// && myID==MASTER)
    {
      world->haplotypes[locus] = world->data->individuals;
      world->numhaplotypes[locus] = world->data->numindividuals;
    }
}

void reset_haplotypes(world_fmt* world, long locus)
{
  (void) locus;
  world->data->individuals = NULL;
  world->data->numindividuals=0;
}

/// looks up the individuals in a simple list for the haplotype_order_list
long find_in_haplotypingorderDB(char * thisname, haplotyping_orderform_fmt *individuals, long numindividuals)
{
  // individuals points to (world->)data->individuals[region]
  long i;
  if(thisname[0]=='\0')
    return numindividuals;

  for (i=0; i< numindividuals; i++)
    {
      haplotyping_orderform_fmt *ii = &individuals[i];
      if(strstr(thisname,ii->key)!=NULL)
	{
	  return i;
	}
    }
  return -1;
}
/// looks up the individuals in a simple list
long find_inDB(char * thisname, long region, individualDB_fmt *individuals, long numindividuals)
{
  // individuals points to (world->)data->individuals[region]
  long i;
  long len = (long) strlen(thisname);

  if(thisname[0]=='\0')
    return -1; //numindividuals;

  for (i=0; i< numindividuals; i++)
    {
      individualDB_fmt *ii = &individuals[i];
      if(len > (long) strlen(ii->name))
	{
	  if(strstr(thisname,ii->name)!=NULL)
	    {
	      if(ii->region==region)
		return i;
	    }
	}
      else
	{
	  if(strstr(ii->name,thisname)!=NULL)
	    {
	      if(ii->region==region)
		return i;
	    }
	}
    }
  return -1;
  //  return numindividuals;
}


void print_haplotype_stat(world_fmt *world, data_fmt *data)
{
    if(world->data->haplotyping_report)
      {
	long allocbufsize = LINESIZE;
	char *buffer = (char *) mycalloc(allocbufsize,sizeof(char));
	long bufsize = 0;
//	long repmax = 0;
//	if (world->options->replicate)
//	  {
//	    if (world->options->replicatenum == 0)
//	      repmax = world->options->lchains;
//	    else
//	      repmax = world->options->replicatenum;
//	  }
//	else
//	  {
//	    repmax = 1;
//	  }
	print_haplotypes(&buffer, &allocbufsize, &bufsize, world, data);
	fprintf(world->outfile, "%s\n",buffer);
	if(world->options->progress)
	  {
	    fprintf(stdout, "%s\n",buffer);
	  }
	myfree(buffer);
      }
}


///
/// save all haplotype assignments from the worker nodes
void
get_haplotypes (world_fmt * world, option_fmt * options)
{
  (void) world;
  (void) options;
#ifdef MPI
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
    
    if (myID == MASTER && world->data->haplotyping_report)
    {
        mpi_results_master (MIGMPI_HAPLOTYPING, world, maxreplicate,
                            unpack_haplotypes_buffer);
    }
#endif
}



/// stores the individuals in a database using the name as indicator
/// all characters before : will be used to find the haplotype-inputs
/// a typical name would look like this
/// ARG34:1
/// ARG34:2
///
void store_indname(char *name, long pop, long ind, long region, data_fmt *data)
{
  char thisname[LINESIZE];
  long id;
  size_t numchar = (size_t) strcspn(name,":");
  strncpy(thisname,name,numchar);
  thisname[numchar]='\0';
  if (-1 == find_in_haplotypingorderDB(thisname, data->haplotyping_list, data->numhaplotyping))
      return;
  id = find_inDB(thisname, region, data->individuals[region], data->numindividuals[region]);
  if (id == -1)
    {
      id = data->numindividuals[region];
      data->numindividuals[region] += 1;
      data->individuals[region] = (individualDB_fmt *) myrealloc(data->individuals[region],sizeof(individualDB_fmt)* (size_t) data->numindividuals[region]);
      if (data->individuals[region][id].ind == NULL)
	data->individuals[region][id].ind = (long *) mycalloc(TWO,sizeof(long));
      data->individuals[region][id].id = id;
      data->individuals[region][id].ind[0] = ind;
      data->individuals[region][id].region = region;
      data->individuals[region][id].nodenum = 0;
      numchar = (size_t) strcspn(name,":");
      strncpy(data->individuals[region][id].name,name,numchar);
      data->individuals[region][id].name[numchar]='\0';
#ifdef DEBUG
      printf("Added new individual to database: %li %s %li %li\n",region, data->individuals[region][id].name, id, data->numindividuals[region]);
#endif
    }
  else
    {
      data->individuals[region][id].ind[1] = ind;
#ifdef DEBUG
      printf("Found individual in database: %s %li %li\n",data->individuals[region][id].name, id, data->numindividuals[region]);
#endif
      if (data->individuals[region][id].region != region)
	warning("Region mismatch for individual %s (%li) in pop %li \n", region, name, ind, pop);
    }
}

void init_individuals(data_fmt **data, long loci)
{
  long locus;
  (*data)->individuals = (individualDB_fmt **) mycalloc(loci, sizeof(individualDB_fmt *));
  
  for(locus=0;locus<loci;locus++)
    (*data)->individuals[locus] = (individualDB_fmt *) calloc(1, sizeof(individualDB_fmt));
  (*data)->numindividuals = (long *) mycalloc(loci, sizeof(long));
}

// for multiple replication we need to reset
void reset_ID_nodelist(long locus, world_fmt * world)
{
  (void) locus;
  long i;
  for(i=0;i<world->data->numindividuals;i++)
    {
      world->data->individuals[i].nodenum = 0;
      if(world->data->individuals[i].nodep != NULL)
	memset(world->data->individuals[i].nodep,0,sizeof(node*) * (size_t) world->data->numindividuals);
    }
}

void link_individual_node(char *rawname, node *thenode, long region, world_fmt *world)
{
  worlddata_fmt *wdata = world->data;
  long subloci = world->sublocistarts[region+1] - world->sublocistarts[region];
  char thisname[LINESIZE];
  size_t numchar = (size_t) strcspn(rawname,":");
  strncpy(thisname,rawname,numchar);
  thisname[numchar]='\0';
  long id = find_inDB(thisname, region, wdata->individuals, wdata->numindividuals);
  long nnum=0;
  if(id == -1)
    {
      return;
      //      usererror("Individual with name \"%s\" (%s)\nwas not found in the haploptype states database\n[File:%s, Line: %s]\n",thisname, rawname,__FILE__,__LINE__);
    }
#ifdef DEBUG
  printf("%i> wdata->individuals[%li].nodep[]=[%s %s %c]\n",myID,id,thenode->truename,thenode->nayme,thenode->tip);
#endif
  if(wdata->individuals[id].nodep == NULL)
    wdata->individuals[id].nodep = (node **) mycalloc(wdata->numindividuals, sizeof(node*));
  nnum = wdata->individuals[id].nodenum;
  wdata->individuals[id].nodep[nnum] = thenode;
  wdata->individuals[id].nodenum += 1;
  if(wdata->individuals[id].nodenum==1) // first filling in of the nodenum is used to establish all arrays
    {
      if(wdata->individuals[id].difference == NULL)
	{
	  wdata->individuals[id].difference = (int *) mycalloc(subloci,sizeof(int));
	}
      else
	{
	  wdata->individuals[id].difference = (int *) myrealloc(wdata->individuals[id].difference, 
								(size_t) subloci * sizeof(int));
	  memset(wdata->individuals[id].difference,0, sizeof (int) * (size_t) subloci); 
	}
      if(wdata->individuals[id].states1==NULL)
	{
	  wdata->individuals[id].states1 = (int *) mycalloc((subloci*2),sizeof(int));
	}
      else
	{
	  wdata->individuals[id].states1 = (int *) myrealloc(wdata->individuals[id].states1, 
							     (size_t) (subloci*2) *sizeof(int));
	  memset(wdata->individuals[id].states1,0, sizeof (int) * (size_t) (2*subloci)); 
	}
      int i;
      for(i=0;i< subloci; i++)
	wdata->individuals[id].states1[i+subloci] = 1;
      
      wdata->individuals[id].total1 = 0.0;
      wdata->individuals[id].total2 = 0.0;
      wdata->individuals[id].numstates = subloci;
    }
}

/// compares haplotype fragments and establishes sameness
boolean is_same_sequence(phenotype * xx1, phenotype * xx2, long numpatterns, long numsiterates, long numstates)
{
  int i, j;
  long site;
  for(site=0; site < numpatterns; site++)
    {
      ratelike x1 = (*xx1)[site];
      ratelike x2 = (*xx2)[site];

      for(i=0;i<numsiterates;i++)
	{
	  for(j=0;j<numstates;j++)
	    {
	      if (fabs(x1[i][j] - x2[i][j])>EPSILON)
		return FALSE;
	    }
	}
    }
  return TRUE;
}

/// compares haplotype fragments and establishes sameness
boolean is_same_allele(MYREAL ** xx1, MYREAL ** xx2, long numpatterns, long numsiterates, long numstates)
{
  (void) numsiterates;
  (void) numstates;
  //int i,j;
  long site;
  for(site=0; site < numpatterns; site++)
    {
      MYREAL x1 =  (*xx1)[site] - (*xx2)[site];
      if (fabs(x1)>EPSILON)
	return FALSE;
    }
  return TRUE;
}



/// checksum to quickly spot individuals that are monomorphic
long calc_individual_checksum(int *difference, world_fmt *world, long locus)
{
  long sublocus;
  long sum=0;

  const long sublocistart = world->sublocistarts[locus];
  const long sublociend   = world->sublocistarts[locus+1];
  const long maxsubloci   = sublociend - sublocistart;
  for(sublocus=0;sublocus<maxsubloci;sublocus++)
    {
      sum +=  (long) difference[sublocus];
    }
  return sum;
}

/// checks whether haplotype fragments are the same and marks the fragments and individuals
/// monomoprhic individuals will be then excluded using the differences
void check_individual_node(char *rawname, node *thenode, long region, world_fmt *world)
{
  (void) thenode;
  worlddata_fmt *wdata = world->data;
  long s1, s2;
  //long site;
  char thisname[LINESIZE];
  size_t numchar = (size_t) strcspn(rawname,":");
  strncpy(thisname,rawname,numchar);
  thisname[numchar]='\0';
  long id = find_inDB(thisname, region, wdata->individuals, wdata->numindividuals);
  if(id == -1)
    {
      return;
      //      usererror("Individual with name \"%s\" (%s)\nwas not found in the haploptype states database\n",thisname, rawname);
    }
  const long sublocistart = world->sublocistarts[region];
  const long sublociend   = world->sublocistarts[region+1];
  const long maxsubloci   = sublociend - sublocistart;

  for(s1=sublocistart, s2=0;s2<maxsubloci; s1++,s2++)
    {
      mutationmodel_fmt * mumod = &world->mutationmodels[s1];
      const long numpatterns = mumod->numpatterns;
      const long numsiterates = mumod->numsiterates;
      const long numstates = mumod->numstates;
      boolean decision;
      if (mumod->dataclass == SITEWORD)
	decision = !(is_same_allele(&(wdata->individuals[id].nodep[0]->x[s2].a), &(wdata->individuals[id].nodep[1]->x[s2].a), numpatterns, numsiterates,numstates));
      else
	decision = !(is_same_sequence(&(wdata->individuals[id].nodep[0]->x[s2].s), &(wdata->individuals[id].nodep[1]->x[s2].s), numpatterns, numsiterates,numstates));
      if(decision)
	{
	  wdata->individuals[id].difference[s2] = 1;    
	}
    }
  wdata->individuals[id].checksum = calc_individual_checksum(wdata->individuals[id].difference, world, region);
}

/// driver to check the individuals nodes for haplotype differences
/// and removing individuals from the database that are monomorphic
void check_individual_nodes(long region, world_fmt *world)
{
  long i;
  node * thenode;
#ifdef DEBUG
  printf("%i> (Heat: %f) executing check_individual_nodes(), world->sumtips=%li\n", myID, world->heat, world->sumtips);
#endif
  for(i=0; i < world->sumtips; i++)
    {
      thenode = world->nodep[i];
#ifdef DEBUG
	  if (thenode->type != 't')
	    {
	      error("failed to show a top in ind assign for heating");
	    }
#endif 	  
      check_individual_node(thenode->truename, thenode, region, world);
    }
}

/// and removing individuals from the database that are monomorphic
void cleanup_individual_nodes(long region, world_fmt *world)
{
  (void) region;
  int i = 0;
  long locus = world->locus;
  individualDB_fmt *ind = world->data->individuals;
  const long last = world->data->numindividuals;
  while (i < last)
    {
      if (ind[i].checksum == 0 || ind[i].numstates<2 || ind[i].targeted == 0)
	{
	  individualDB_fmt tmp = ind[i];
	  world->data->numindividuals -= 1;
	  long end = world->data->numindividuals;
	  if(i==end)
            break;
	  ind[i] = ind[end];
	  ind[end] = tmp;
	  //free(ind[end].difference);
	  //free(ind[end].states1);
	  // this will test the individual we just swapped in
#ifdef DEBUG
	  printf("%i> (Heat: %f) %s removed from total %li in region %li\n",myID, world->heat, 
		 ind[end].name, world->data->numindividuals, region);
#endif
	  i--;
	}
      i++;
      if(i>=world->data->numindividuals)
	break;
    }
  if(world->cold)
    {
      world->numhaplotypes[locus] = world->data->numindividuals;
    }
  if(world->data->numindividuals==0)
    {
      world->data->haplotyping=FALSE;
    }
  else
    {
      //world->data->individuals = (individualDB_fmt *) realloc(world->data->individuals,sizeof(individualDB_fmt)*world->data->numindividuals);
    }
}

/// copies the main haplotype structure from data to world->data
void copy_individuals_from_data(data_fmt * data, worlddata_fmt * wdata, long locus)
{
  long id;
  if(wdata->individuals==NULL)
    wdata->individuals = (individualDB_fmt *) mycalloc(data->numindividuals[locus],sizeof(individualDB_fmt));
  else
    {
      //wdata->individuals = (individualDB_fmt *) realloc(wdata->individuals, data->numindividuals[locus] * sizeof(individualDB_fmt));
      //return;
    }
  wdata->numindividuals = data->numindividuals[locus];
  for (id=0; id < data->numindividuals[locus]; id++)
    {
      wdata->individuals[id].id = id;
      if (wdata->individuals[id].ind == NULL)
	wdata->individuals[id].ind = (long *) mycalloc(TWO,sizeof(long));
      memcpy(wdata->individuals[id].ind, data->individuals[locus][id].ind, sizeof(long)*TWO);
      wdata->individuals[id].region = data->individuals[locus][id].region;
      wdata->individuals[id].nodenum = 0;
      strcpy(wdata->individuals[id].name,data->individuals[locus][id].name);
    }
}


long  set_haplotype_states(individualDB_fmt *winner, char *word)
{
  if(winner->states1==NULL)
    {
      winner->states1 = (int *) calloc(strlen(word),sizeof(int));
    }
  else
    {
      winner->states1 = (int *) realloc(winner->states1,strlen(word)*sizeof(int));
    }
  long ii=0;
  while(word[ii]!='\0')
    {
      if(word[ii]=='1')
	winner->states1[ii]=1;
      else
	winner->states1[ii]=0;
      ii++;
    }
  return ii;
}

void debug_states(individualDB_fmt *w)
{
  long i;
  long sum1=0;
  long sum2=0;
  for(i=0; i<w->numstates;i++)
    {
      sum1 += w->states1[i];
      sum2 += w->states1[i+w->numstates];
    }
  if(sum1+sum2<w->numstates)
    warning("DEBUG total number of states: %li\n",sum1+sum2);
}


/// update haplotypes
long swap_haplotypes(world_fmt *world)
{
  long locus = world->locus;
  long sublocus;
  long sublocus2;
  const long sublociend   = world->sublocistarts[locus+1]-world->sublocistarts[locus];
  if (sublociend==1)
    return -1L; // no swapping of haplotype needed
  long winnerid = RANDINT(0,world->data->numindividuals-1);
  long count=0;
  individualDB_fmt *winner = &(world->data->individuals[winnerid]);
  long rh = RANDINT(0, winner->checksum-1);
  count = 0;
  long i=0;
#ifdef DEBUG
  debug_states(winner);
#endif
  while(1)
    {
      if(winner->difference[i] != 0 )
	count++;
      if(count>rh)
	break;
      i++;
    }
  sublocus = i;
  sublocus2 = i + winner->numstates;
  xarray_fmt one;
  one = winner->nodep[0]->x[sublocus];
  winner->nodep[0]->x[sublocus] = winner->nodep[1]->x[sublocus];
  winner->nodep[1]->x[sublocus] = one;
  // accept or reject move
  MYREAL oldlike = world->likelihood[world->G];
  set_tree_down_dirty (winner->nodep[0]);      
  set_tree_down_dirty (winner->nodep[1]);
  //      set_tree_dirty(world->root);
  smooth (world->root->next, crawlback (world->root->next), world, locus);
  MYREAL newlike = treelikelihood(world);
  if (newlike==0.0)
    error("likelihood zeroed");
  if (oldlike < newlike)
    {
      //accept
      world->likelihood[world->G] = newlike;
      // original state:  0 0 0
      //               :  1 1 1
      // 000,001,010,100,011,101,110,111
      ///111,110,101,011,100,010,001,000
      // these are 2^n choices: say n=10 -> 1024,
      // with many more it will be difficult to trace them all 
      // tracing should be able to be turned off and on depending on need also should be able
      // to write to  file
      long tmp = winner->states1[sublocus];
      winner->states1[sublocus] = winner->states1[sublocus2];
      winner->states1[sublocus2] = (int) tmp;
#ifdef DEBUG
      // printf("S %f %f %i %i | %i %i | %li %li | %f\n",newlike,oldlike,
      //	     winner->states1[0], winner->states1[1],
      //	     winner->states1[2], winner->states1[2], winnerid, sublocus, world->heat);
#endif      
    }
  else
    {
      if(LOG(RANDUM()) < (newlike - oldlike))
	{
	  //accept
	  world->likelihood[world->G] = newlike;
	  long tmp = winner->states1[sublocus];
	  winner->states1[sublocus] = winner->states1[sublocus2];
	  winner->states1[sublocus2] = (int) tmp;
#ifdef DEBUG
	  //    printf("S %f %f %i %i | %i %i | %li %li | %f\n",newlike,oldlike,
	  // winner->states1[0], winner->states1[1],
	  // winner->states1[2], winner->states1[3], winnerid, sublocus, world->heat);
#endif
	}
      else
	{
	  one = winner->nodep[0]->x[sublocus];
	  winner->nodep[0]->x[sublocus] = winner->nodep[1]->x[sublocus];
	  winner->nodep[1]->x[sublocus] = one;
	  set_tree_down_dirty (winner->nodep[0]);      
	  set_tree_down_dirty (winner->nodep[1]);
#ifdef DEBUG
	  //	  	  printf("F %f %f %i %i | %i %i | %li %li | %f\n",newlike,oldlike,
	  // winner->states1[0], winner->states1[1],
	  // winner->states1[2], winner->states1[3], winnerid, sublocus, world->heat);
#endif
	  smooth (world->root->next, crawlback (world->root->next), world, locus);
	  newlike = treelikelihood(world);
#ifdef DEBUGXX
	  if(fabs(newlike-oldlike)> (double) FLT_EPSILON)
	    {
	      set_tree_dirty(world->root);
	      smooth (world->root->next, crawlback (world->root->next), world, locus);
	      newlike = treelikelihood(world);
	      if(fabs(newlike-oldlike)> (double) FLT_EPSILON)
		{
		  warning("Haplotyping: newL=%f != oldL=%f\n", newlike, oldlike);
		}
	    }
	  //	  else
	  //  {
	  //    printf("OK ");
	  //  }
#endif
	} 
    }
  if(world->data->haplotyping_report)
    {
      add_state_counter(&winner->hash, &winner->numhash, &winner->total1, 
			winner->states1, winner->numstates, 1);
    }
  return -1L;
}
  
void free_haplotypes(world_fmt *world, long locus)
{
  individualDB_fmt *ii = world->haplotypes[locus];
  if(ii->difference!=NULL)
    myfree(ii->difference);
  if(ii->states1!=NULL)
    myfree(ii->states1);
  if(ii->states2!=NULL)  
    myfree(ii->states2);
  free_hash(ii->hash,ii->numhash,(size_t) ii->numstates);
  ii = NULL;
}
