/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
 
#ifndef __PRPTSRST_H
#define __PRPTSRST_H

#include "../util.h"
#include "../permutations.h"


void alg_2PSRT(permutation_t *p, int verbose, info_t* info);
void alg_2SPSRT(permutation_t *p, int verbose, info_t *info);
 
#endif /* __PRPTSRST_H */

