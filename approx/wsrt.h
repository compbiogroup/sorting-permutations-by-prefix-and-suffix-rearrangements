/**********************************************************
 * Created: Seg 17 Mar 2014 02:52:06 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __WSRST_H
#define __WSRST_H
 
#include "../util.h"
#include "../permutations.h"


long double partitionWSRT(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info);
long double WSRT(permutation_t *p, int ini, int type, int verbose, info_t *info);
void alg_WSRT(permutation_t *p, int verbose, info_t *info);



long double partitionWSSRT(permutation_t *p, int ini, int type, int divisor, int verbose, info_t *info);
long double WSSRT(permutation_t *p, int ini, int type, int verbose, info_t *info);
void alg_WSSRT(permutation_t *p, int verbose, info_t *info);


#endif /* __WSRST_H */


