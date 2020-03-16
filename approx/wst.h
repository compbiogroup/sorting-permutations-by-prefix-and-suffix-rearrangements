/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __WEIGHTED_SUFFIX_TRANSPOSITION_H
#define __WEIGHTED_SUFFIX_TRANSPOSITION_H

#include "../util.h"
#include "../permutations.h"

long double WST(permutation_t *p, int ini, int verbose, info_t *info);
long double partitionWST(permutation_t *p, int ini, int divisor, int verbose, info_t *info);

void alg_WST(permutation_t *p, int verbose, info_t* info);

#endif /* __WEIGHTED_SUFFIX_TRANSPOSITION_H */

