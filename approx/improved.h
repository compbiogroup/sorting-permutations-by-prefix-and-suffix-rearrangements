/**********************************************************
 * Created: Sáb 14 Dez 2013 17:58:52 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __GREEDY_APPROX_H
#define __GREEDY_APPROX_H

#include "../util.h"
#include "../permutations.h"

void alg_2PRx(permutation_t *p, int verbose, info_t* info);
void alg_2PSRx(permutation_t *p, int verbose, info_t* info);
void alg_2PTx(permutation_t *p, int verbose, info_t* info);
void alg_2PSTx(permutation_t *p, int verbose, info_t* info);
void alg_2PRTx(permutation_t *p, int verbose, info_t* info);
void alg_2PSRTx(permutation_t *p, int verbose, info_t* info);
 
void alg_2SPRx(permutation_t *p, int verbose, info_t* info);
void alg_2SPSRx(permutation_t *p, int verbose, info_t* info);
void alg_2SPRTx(permutation_t *p, int verbose, info_t* info);
void alg_2SPSRTx(permutation_t *p, int verbose, info_t* info);
 
#endif /* __GREEDY_APPROX_H */

