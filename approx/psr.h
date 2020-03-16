/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __PREFIX_SUFFIX_REVERSAL_H
#define __PREFIX_SUFFIX_REVERSAL_H

#include "../util.h"
#include "../permutations.h"

void psr_star_perm(permutation_t *p, int verbose, info_t *info, int nofbps);
int psr_edge_type1p(int *edge, permutation_t *p);
int psr_edge_type2p(int *edge, permutation_t *p);
int psr_edge_type3Ap(int *edge, permutation_t *p);
int psr_edge_type3Bp(int *edge, permutation_t *p);
int psr_edge_type1s(int *edge, permutation_t *p);
int psr_edge_type2s(int *edge, permutation_t *p);
int psr_edge_type3As(int *edge, permutation_t *p);
int psr_edge_type3Bs(int *edge, permutation_t *p);
void alg_2PSR(permutation_t *p, int verbose, info_t* info);



void spsr_star_perm(permutation_t *p, int verbose, info_t *info, int nofbps);
void alg_2SPSR(permutation_t *p, int verbose, info_t *info);

 
#endif 

