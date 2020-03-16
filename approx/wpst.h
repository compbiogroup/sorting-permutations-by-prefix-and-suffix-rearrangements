/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WPTST_H
#define __WPTST_H
 
#include "../util.h"
#include "../permutations.h"

void alg_WPSTg(permutation_t *p, int verbose, info_t* info);

int partitionWPST(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info);
int WPST(permutation_t *p, int ini, int type, int verbose, info_t *info);
void alg_WPST(permutation_t *p, int verbose, info_t* info);


#endif /* __WPRSR_H */

