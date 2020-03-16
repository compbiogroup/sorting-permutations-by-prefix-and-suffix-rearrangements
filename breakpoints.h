/**********************************************************
 * Created: Seg 14 Out 2013 15:14:48 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __BREAKPOINTS_H
#define __BREAKPOINTS_H

#include "permutations.h"
 
int is_UPR_breakpoint(permutation_t *p, int i, int j);
int is_UPSR_breakpoint(permutation_t *p, int i, int j);
int is_UPSRT_breakpoint(permutation_t *p, int i, int j);
int is_P_breakpoint(permutation_t *p, int i, int j);
int is_PS_breakpoint(permutation_t *p, int i, int j);
int is_PSRT_breakpoint(permutation_t *p, int i, int j);
int is_breakpoint(char name[], permutation_t *p, int i, int j);

int nof_UPR_breakpoints(permutation_t *p);
int nof_UPSR_breakpoints(permutation_t *p);
int nof_UPSRT_breakpoints(permutation_t *p);
int nof_P_breakpoints(permutation_t *p);
int nof_PS_breakpoints(permutation_t *p);
int nof_PSRT_breakpoints(permutation_t *p);
int nof_breakpoints(char name[], permutation_t *p);

int elem_is_UPR_breakpoint(permutation_t *p, int i, int j);
int elem_is_UPSR_breakpoint(permutation_t *p, int i, int j);
int elem_is_UPSRT_breakpoint(permutation_t *p, int i, int j);
int elem_is_P_breakpoint(permutation_t *p, int i, int j);
int elem_is_PS_breakpoint(permutation_t *p, int i, int j);
int elem_is_PSRT_breakpoint(permutation_t *p, int i, int j);

int nofbps_after_pr(permutation_t *p, int i, int nofbp);
int nofbps_after_psr(permutation_t *p, char type, int i, int nofbp);
int nofbps_after_pt(permutation_t *p, int i, int j, int nofbp);
int nofbps_after_pst(permutation_t *p, char type, int i, int j, int nofbp);
int nofbps_after_prt(permutation_t *p, char rearr, int i, int j, int nofbp);
int nofbps_after_psrt(permutation_t *p, char rearr, char type, int i, int j, int nofbp);
int nofbps_after_spr(permutation_t *p, int i, int nofbp);
int nofbps_after_spsr(permutation_t *p, char type, int i, int nofbp);
int nofbps_after_sprt(permutation_t *p, char rearr, int i, int j, int nofbp);
int nofbps_after_spsrt(permutation_t *p, char rearr, char type, int i, int j, int nofbp);
int nofbps_after(char prob[], char type, char rearr, permutation_t *p, int i, int j, int k, int nofbp);


#endif /* __BREAKPOINTS_H */

