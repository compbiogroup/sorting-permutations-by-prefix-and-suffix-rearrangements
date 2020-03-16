/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WSR_H
#define __WSR_H
 
#include "../util.h"
#include "../permutations.h"

long double partitionWSR(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info);
long double WSR(permutation_t *p, int ini, int type, int verbose, info_t *info);

void alg_WSR(permutation_t *p, int verbose, info_t* info);



long double partitionWSSR(permutation_t *p, permutation_t *sp, int ini, int type, int divisor, int verbose, info_t *info);
long double WSSR(permutation_t *p, permutation_t *sp, int ini, int type, int verbose, info_t *info);

void alg_WSSR(permutation_t *p, int verbose, info_t *info);

#endif /* __WSR_H */

