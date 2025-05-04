/// hash functionality
/// for migrate, currently uses string with 0 1 characters of 'any' length
/// and adding a count to the values, this is geared towards counting patterns
/// but may be easily expanded to more complex storage devices
/// currently the keys are partitioned into 1<<9 blocks and subkeys are link hashtables
/// I guess they call this a perfect hash, but this one is sparse and not all subtrees are completely
/// allocated.
///
/*
  Copyright (c) 2015 Peter Beerli

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/
#include "migration.h"
#include "sighandler.h"
#include <stdio.h>
#include <stdlib.h>
#include "hash.h"

void create_hash(hash_fmt **hash, size_t number_of_binary_states);
void reset_hash(hash_fmt *hash, size_t numhash, size_t numstates);

void create_hash(hash_fmt **hash, size_t number_of_binary_states)//columns
{
  const int firstset = PARTITIONING;
  if(number_of_binary_states <= firstset)
    {
      (*hash) = (hash_fmt *) calloc(1<<number_of_binary_states,sizeof(hash_fmt));
    }
  else
    {
      printf("in create_hash wrong else condition found");
      exit(-1);
      /*
      int i;
      int first = 1<<firstset;
      (*hash) = (hash_fmt *) mycalloc(first,sizeof(hash_fmt));
      for(i = 0; i < first; i++) 
	{
	  hash_fmt * hh = &((*hash)[i]);
	  create_hash(&hh->next, number_of_binary_states - firstset);
	}
      */
    }
}

boolean insert_hash(char *key, int *offset, size_t numstates, size_t realnumstates, long ptr, hash_fmt **hash, size_t *numhash)
{
  //boolean done=FALSE;
  size_t i,j;
  size_t first=PARTITIONING;
  size_t index = 0;
  if(numstates<=first)
    {
      char *nkey = key + realnumstates - numstates;;
      for(i=0, j=numstates-1; i<numstates; i++, j--)
	{
	  index +=  ((size_t) nkey[i] - 48) * (1 << j);
	}
      if(index < *numhash)
	{
	  hash_fmt * hh = &((*hash)[index]);
	  hh->next = NULL;
	  hh->value += ptr;
	  if(hh->key==NULL)
	    hh->key = (char*) mycalloc((realnumstates+1),sizeof(char));
	  strncpy(hh->key,key,realnumstates);
	  hh->key[realnumstates]='\0';
	}
      else
	{
	  create_hash(hash,numstates);
	  *numhash += 1<<numstates;
	  hash_fmt * hh = &((*hash)[index]);
	  hh->next = NULL;
	  hh->value += ptr;
	  if(hh->key==NULL)
	    hh->key = (char*) mycalloc((realnumstates+1),sizeof(char));
	  strncpy(hh->key,key,realnumstates);
	  hh->key[realnumstates]='\0';
	}
      return TRUE;
    }
  else
    {
      for(i=0, j=first-1; i<first; i++, j--)
	{
	  char *nkey = key + realnumstates - numstates;
	  index +=  ((size_t) nkey[i] - 48) * (1 << j);
	}
      *offset += first;

      if(index < *numhash)
	{
	  hash_fmt * hh = &((*hash)[index]);
	  insert_hash(key, offset, numstates-first, realnumstates, ptr, &(hh->next), &(hh->numhash));
	}
      else
	{
	  create_hash(hash,first);
	  *numhash = 1<<first;
	  hash_fmt * hh = &((*hash)[index]);
	  insert_hash(key, offset, numstates-first, realnumstates, ptr, &(hh->next), &(hh->numhash));
	}
    }
  return FALSE;
}

/*
boolean insert_hash_old(char *key, int offset, int numstates, int realnumstates, long ptr, hash_fmt **hash, size_t *numhash)
{
  boolean done=FALSE;
  int i,j;
  int first=PARTITIONING;
  int index = 0;
  if(numstates<=first)
    {
      // only one level
      for(i=0, j=numstates-1; i<numstates; i++, j--)
	{
	  index +=  ((int) (key + offset)[i] - 48) * (1 << j);
	}
      if(index < *numhash)
	{
	  (*hash)[index].value += ptr;
	  if((*hash)[index].key[0] == '\0')
	    {
	      strncpy((*hash)[index].key,key,realnumstates);
	    }
	}
      else
	{
	  create_hash(hash,numstates);
	  *numhash += 1<<numstates;
	  (*hash)[index].value += ptr;
	  if((*hash)[index].key[0] == '\0')
	    {
	      strncpy((*hash)[index].key,key,realnumstates);
	    }
	}
      return TRUE;
    }
  else
    {
      offset=0;
      while(!done)
	{
	  for(i=0, j=first-1; i<first; i++, j--)
	    {
	      index +=  ((int) (key+offset)[i] - 48) * (1 << j);
	    }
	  offset += first;
	  if (numstates-offset <=0)
	    done=TRUE;
	  if(index < *numhash)
	    {
	      done = insert_hash_old(key, offset, numstates-offset, realnumstates, ptr, &((*hash)[index].next), &((*hash)[index].numhash));
	    }
	  else
	    {
	      create_hash(hash,first);
	      *numhash = 1<<first;
	      done = insert_hash_old(key, offset,numstates-offset, realnumstates, ptr, &(*hash)[index].next, &(*hash)[index].numhash);
	    }
	}
      return FALSE;
    }
}
*/

long find_hash(char *key, size_t numstates, hash_fmt *hash, size_t *numhash)
{
  size_t index = 0;
  const int first = PARTITIONING;
  size_t i, j;
  if(numstates <= first)
    {
      for(i=0, j=numstates-1; i<numstates; i++, j--)
	{
	  index +=  ((size_t) key[i] - 48) * (1 << j);
	}
      if(index < *numhash)
	return hash[index].value;
      else
	{
	  return 0;
	}
    }
  else
    {
      for(i=0, j=first-1; i<first; i++, j--)
	{
	  index +=  ((size_t) key[i] - 48) * (1 << j);
	}
      return find_hash(key+first,numstates-first, hash[index].next, &hash[index].numhash);
    }
}

void hash_elements(hash_fmt *hash, size_t numhash, long *total)
{
  size_t index;
  for(index = 0; index < numhash; index++)
    {
      hash_fmt * hh = &(hash[index]);
      if(hh->next!=NULL)
	{
	  hash_elements(hh->next,hh->numhash, total);
	}
    }
  *total += numhash;
}

void hash_leaves(hash_fmt *hash, size_t numhash, long *leaves)
{
  size_t index;
  for(index = 0; index < numhash; index++)
    {
      hash_fmt * hh = &(hash[index]);
      if(hh->next!=NULL)
	{
	  hash_leaves(hh->next,hh->numhash, leaves);
	}
      else
	{
	  *leaves += 1;
	}
    }
}

void show_hash(hash_fmt *hash, size_t numhash, size_t numstates, char ***showkeys, int **showvalues, long *z)
{
  size_t index;

  for(index = 0; index < numhash; index++)
    {
      hash_fmt * hh = &(hash[index]);
      if(hh->next!=NULL)
	{
	  show_hash(hh->next, hh->numhash, numstates, showkeys,showvalues, z);
	}
      else
	{
	  if(hh->key!=NULL)
	    {
	      strncat((*showkeys)[*z], hh->key, numstates);
	      //printf("retrieve---> %s %s %i\n", (*showkeys)[*z], hh->key, numstates);
	      (*showvalues)[*z] = (int) hh->value;
	      (*z) += 1;
	    }
	}
    }
}

// resetting all values in the hash table to zero
// this is necessary to reuse the table without reallocating the structure
void reset_hash(hash_fmt *hash, size_t numhash, size_t numstates)
{
  size_t index;

  for(index = 0; index < numhash; index++)
    {
      hash_fmt * hh = &(hash[index]);
      if(hh->next!=NULL)
	{
	  reset_hash(hh->next, hh->numhash, numstates);
	}
      else
	{
	  if(hh->key!=NULL)
	    {
	      hh->value=0;
	    }
	}
    }
}
void free_hash(hash_fmt *hash, size_t numhash, size_t numstates)
{
  size_t index;

  for(index = 0; index < numhash; index++)
    {
      hash_fmt * hh = &(hash[index]);
      if(hh->next!=NULL)
	{
	  free_hash(hh->next, hh->numhash, numstates);
	}
    }
  free(hash);
}

// autonomous TEST facility
#ifdef HASHTEST
int main(long argc, char **argv)
{
  int v;
  // 0000 
  // 0001 0010 0100 1000 
  // 0011 0110 1100 0101 1010 1001
  // 0111 1110 1101 1011
  // 1111 
  
  //char mKEYS[][12]={"01010000011","00000000001","01000000000","00010000000","01110000000","11110000000"};
  char mKEYS[][22]={"0101000001100000000011","0101000001100000000001","0101000001101000000000","0101000001100010000000","0101000001101110000000","0101000001111110000000"};
  char **keys;
  hash_fmt *hash;
  size_t numhash=0;
  int key;
  size_t numstates=22;
  keys = (char **) calloc(6,sizeof(char *));
  for(key=0; key<5; key++)
    {
      keys[key] = calloc(numstates+1,sizeof(char));
      strncpy(keys[key],mKEYS[key],numstates);
    }
  keys[key] = calloc(numstates+1,sizeof(char));
  int value=1;
  int offset = 0;
  for(key=0; key<5; key++)
    {
      insert_hash(keys[key], &offset, numstates, numstates, value, &hash, &numhash);
    }
  for(key=0; key<3; key++)
    {
      insert_hash(keys[key], &offset, numstates, numstates, value, &hash, &numhash);
      insert_hash(keys[4], &offset, numstates, numstates, value, &hash, &numhash);
    }
  strncpy(keys[5],mKEYS[5],numstates);
  for(key=0; key<5; key++)
    {
      v = find_hash(keys[key],numstates, hash, &numhash);
      printf("%s: %i\n",keys[key], v);
    }
  long z=0;
  long total=0;
  long leaves=0;
  hash_elements(hash,numhash, &total);
  hash_leaves(hash,numhash, &leaves);
  printf("\nTotal elements in hash: %li\n",total);
  printf("Total leaves in hash  : %li\n\n",leaves);
  char ** showkeys = (char **) calloc(leaves,sizeof(char*));
  int i;
  for (i=0; i < leaves; i++)
    showkeys[i] = (char *) calloc(numstates,sizeof(char));
  int *showvalues = (int *) calloc(leaves,sizeof(int));
  show_hash(hash,numhash,numstates, &showkeys,&showvalues,&z);
  for (i=0; i < leaves; i++)
    {
      if(showkeys[i][0]!='\0')
	printf("Key: %s  Value: %i\n", showkeys[i], showvalues[i]);
    }
  return 0;
}

#endif
