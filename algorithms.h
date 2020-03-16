/**********************************************************
 * Created: Sáb 14 Dez 2013 17:45:49 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __ALGORITHMS_H
#define __ALGORITHMS_H
 
#include "util.h"

#define PR2 1
#define PRx2 2
#define PSR2 3
#define PSRx2 4

#define PT2 5
#define PTx2 6
#define PST2 7
#define PSTx2 8

#define PRT2 9
#define PRTx2 10
#define PSRT2 11
#define PSRTx2 12

#define SPR2 13
#define SPRx2 14
#define SPSR2 15
#define SPSRx2 16

#define SPRT2 17
#define SPRTx2 18
#define SPSRT2 19
#define SPSRTx2 20

#define cWPRm 21
#define cWPRg 22
#define cWPR 23
#define cWSR 24
#define cWPSRg 25
#define cWPSR 26

#define cWPTg 27
#define cWPT 28
#define cWST 29
#define cWPSTg 30
#define cWPST 31

#define cWPRTg 32
#define cWPRT 33
#define cWSRT 34
#define cWPSRTg 35
#define cWPSRT 36

#define cWSPRg 37
#define cWSPR 38
#define cWSSR 39
#define cWSPSRg 40
#define cWSPSR 41

#define cWSPRTg 42
#define cWSPRT 43
#define cWSSRT 44
#define cWSPSRTg 45
#define cWSPSRT 46


void deal_alg(char *alg_name, int *alg);
void run(int alg, permutation_t *p, int verbose, info_t *info);
void print_algs(void);
void alg_approx(char prob[], permutation_t *p, int verbose, info_t *info);

#endif /* __ALGORITHMS_H */

