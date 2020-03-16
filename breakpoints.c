/***********************************************************
 * Created: Seg 14 Out 2013 15:14:52 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "permutations.h"
#include "breakpoints.h"


int is_UPR_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_UPR_breakpoint(%d, %d)\n", i, j);
        exit(0);
    }

    if (i == 0)
        return 0;
    return (abs(p->pi[j] - p->pi[i]) != 1);
}/*}}}*/

int is_UPSR_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_UPSR_breakpoint(%d, %d)\n", i, j);
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    return (abs(p->pi[j] - p->pi[i]) != 1);
}/*}}}*/

int is_UPSRT_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_UPSRT_breakpoint(%d, %d)\n", i, j);
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    if (p->pi[i] == 1 && p->pi[j] == p->size)
        return 0;
    if (p->pi[i] == p->size && p->pi[j] == 1)
        return 0;
    return (abs(p->pi[j] - p->pi[i]) != 1);
}/*}}}*/

int is_P_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_P_breakpoint(%d, %d)\n", i, j);
        exit(0);
    }

    if (i == 0)
        return 0;
    return (p->pi[j] - p->pi[i] != 1);
}/*}}}*/

int is_PS_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_PS_breakpoint(%d, %d)\n", i, j);
        print_permutation(p);
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    return (p->pi[j] - p->pi[i] != 1);
}/*}}}*/

int is_PSRT_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0 || j != i+1) {
        printf("ERRO! is_PSRT_breakpoint(%d, %d)\n", i, j);
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    if (p->pi[i] == -1 && p->pi[j] == -p->size)
        return 0;
    if (p->pi[i] == p->size && p->pi[j] == 1)
        return 0;
    return (p->pi[j] - p->pi[i] != 1);
}/*}}}*/


int is_breakpoint(char name[], permutation_t *p, int i, int j) {/*{{{*/
    if (!strcmp(name, "pr") || !strcmp(name, "prt") ||
                    !strcmp(name, "wpr") || !strcmp(name, "WPRT"))
        return is_UPR_breakpoint(p, i, j);

    else if (!strcmp(name, "psr") || !strcmp(name, "wpsr"))
        return is_UPSR_breakpoint(p, i, j);

    else if (!strcmp(name, "psrt") || !strcmp(name, "wpsrt"))
        return is_UPSRT_breakpoint(p, i, j);

    else if (!strcmp(name, "pt") || !strcmp(name, "spr") || !strcmp(name, "sprt") ||
                    !strcmp(name, "wspr") || !strcmp(name, "wsprt"))
        return is_P_breakpoint(p, i, j);

    else if (!strcmp(name, "pst") || !strcmp(name, "spsr") ||
                    !strcmp(name, "wpst") || !strcmp(name, "wspsr"))
        return is_PS_breakpoint(p, i, j);

    else if (!strcmp(name, "spsrt") || !strcmp(name, "wspsrt"))
        return is_PSRT_breakpoint(p, i, j);

    else {
        printf("ERROR: is_breakpoint(%s, p)\n", name);
        exit(0);
    }
}/*}}}*/




int nof_UPR_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_UPR_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/

int nof_UPSR_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_UPSR_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/

int nof_UPSRT_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_UPSRT_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/

int nof_P_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_P_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/

int nof_PS_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_PS_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/

int nof_PSRT_breakpoints(permutation_t *p) {/*{{{*/
    int i, bp = 0;
    for (i = 0; i <= p->size; i++)
        if (is_PSRT_breakpoint(p, i, i+1))
            bp++;
    return bp;
}/*}}}*/


int nof_breakpoints(char name[], permutation_t *p) {/*{{{*/
    if (!strcmp(name, "pr") || !strcmp(name, "prt") ||
                    !strcmp(name, "wpr") || !strcmp(name, "WPRT"))
        return nof_UPR_breakpoints(p);

    else if (!strcmp(name, "psr") || !strcmp(name, "wpsr"))
        return nof_UPSR_breakpoints(p);

    else if (!strcmp(name, "psrt") || !strcmp(name, "wpsrt"))
        return nof_UPSRT_breakpoints(p);

    else if (!strcmp(name, "pt") || !strcmp(name, "spr") || !strcmp(name, "sprt") ||
                    !strcmp(name, "wspr") || !strcmp(name, "wsprt"))
        return nof_P_breakpoints(p);

    else if (!strcmp(name, "pst") || !strcmp(name, "spsr") ||
                    !strcmp(name, "wpst") || !strcmp(name, "wspsr"))
        return nof_PS_breakpoints(p);

    else if (!strcmp(name, "spsrt") || !strcmp(name, "wspsrt"))
        return nof_PSRT_breakpoints(p);

    else {
        printf("ERROR: nof_breakpoints(%s, p)\n", name);
        exit(0);
    }
}/*}}}*/




int elem_is_UPR_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_UPR_breakpoint()\n");
        exit(0);
    }

    if (i == 0)
        return 0;
    return (abs(j - i) != 1);
}/*}}}*/

int elem_is_UPSR_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_UPSR_breakpoint()\n");
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    return (abs(j - i) != 1);
}/*}}}*/

int elem_is_UPSRT_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_UPSRT_breakpoint()\n");
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    if (i == 1 && j == p->size)
        return 0;
    if (i == p->size && j == 1)
        return 0;
    return (abs(j - i) != 1);
}/*}}}*/

int elem_is_P_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_P_breakpoint()\n");
        exit(0);
    }

    if (i == 0)
        return 0;
    return (j - i != 1);
}/*}}}*/

int elem_is_PS_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_PS_breakpoint()\n");
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    return (j - i != 1);
}/*}}}*/

int elem_is_PSRT_breakpoint(permutation_t *p, int i, int j) {/*{{{*/
    if (i == p->size+1 || j == 0) {
        printf("ERRO! elem_is_PSRT_breakpoint()\n");
        exit(0);
    }

    if (j == p->size+1 || i == 0)
        return 0;
    if (i == -1 && j == -p->size)
        return 0;
    if (i == p->size && j == 1)
        return 0;
    return (j - i != 1);
}/*}}}*/




int nofbps_after_pr(permutation_t *p, int i, int nofbp) {/*{{{*/
    if (elem_is_UPR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
        !elem_is_UPR_breakpoint(p, p->pi[1], p->pi[i+1]))
        nofbp--;

    else if (!elem_is_UPR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
             elem_is_UPR_breakpoint(p, p->pi[1], p->pi[i+1]))
        nofbp++;
    
    return nofbp;
}/*}}}*/

int nofbps_after_psr(permutation_t *p, char type, int i, int nofbp) {/*{{{*/
    if (type == 'p') {
        if (elem_is_UPSR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_UPSR_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_UPSR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_UPSR_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_UPSR_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_UPSR_breakpoint(p, p->pi[i-1], p->pi[p->size]))
            nofbp--;

        else if (!elem_is_UPSR_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_UPSR_breakpoint(p, p->pi[i-1], p->pi[p->size]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_pt(permutation_t *p, int i, int j, int nofbp) {/*{{{*/
    if (elem_is_P_breakpoint(p, p->pi[j-1], p->pi[j]) &&
        !elem_is_P_breakpoint(p, p->pi[i-1], p->pi[j]))
        nofbp--;
    else if (!elem_is_P_breakpoint(p, p->pi[j-1], p->pi[j]) &&
             elem_is_P_breakpoint(p, p->pi[i-1], p->pi[j]))
        nofbp++;

    if (elem_is_P_breakpoint(p, p->pi[i-1], p->pi[i]) &&
        !elem_is_P_breakpoint(p, p->pi[j-1], p->pi[1]))
        nofbp--;
    else if (!elem_is_P_breakpoint(p, p->pi[i-1], p->pi[i]) &&
             elem_is_P_breakpoint(p, p->pi[j-1], p->pi[1]))
        nofbp++;

    return nofbp;
}/*}}}*/

int nofbps_after_pst(permutation_t *p, char type, int i, int j, int nofbp) {/*{{{*/
    if (type == 'p') {
        if (elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (!elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp--;
        else if (!elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp++;
    } else {
        if (elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 !elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_PS_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp--;
        else if (!elem_is_PS_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_PS_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_prt(permutation_t *p, char rearr, int i, int j, int nofbp) {/*{{{*/
    if (rearr == 'r') {
        if (elem_is_UPR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_UPR_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_UPR_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_UPR_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_UPR_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_UPR_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (!elem_is_UPR_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_UPR_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_UPR_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_UPR_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp--;
        else if (!elem_is_UPR_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_UPR_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_psrt(permutation_t *p, char rearr, char type, int i, int j, int nofbp) {/*{{{*/
    if (rearr == 't' && type == 'p') {
        if (elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (!elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp--;
        else if (!elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp++;
    } else if (rearr == 't' && type == 's') {
        if (elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 !elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp--;
        else if (!elem_is_UPSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_UPSRT_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp++;
    } else if (rearr == 'r' && type == 'p') {
        if (elem_is_UPSRT_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_UPSRT_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_UPSRT_breakpoint(p, p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[p->size]))
            nofbp--;

        else if (!elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_UPSRT_breakpoint(p, p->pi[i-1], p->pi[p->size]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_spr(permutation_t *p, int i, int nofbp) {/*{{{*/
    if (elem_is_P_breakpoint(p, p->pi[i], p->pi[i+1]) &&
        !elem_is_P_breakpoint(p, -p->pi[1], p->pi[i+1]))
        nofbp--;

    else if (!elem_is_P_breakpoint(p, p->pi[i], p->pi[i+1]) &&
             elem_is_P_breakpoint(p, -p->pi[1], p->pi[i+1]))
        nofbp++;
    
    return nofbp;
}/*}}}*/

int nofbps_after_spsr(permutation_t *p, char type, int i, int nofbp) {/*{{{*/
    if (type == 'p') {
        if (elem_is_PS_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_PS_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_PS_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_PS_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PS_breakpoint(p, p->pi[i-1], -p->pi[p->size]))
            nofbp--;

        else if (!elem_is_PS_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_PS_breakpoint(p, p->pi[i-1], -p->pi[p->size]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_sprt(permutation_t *p, char rearr, int i, int j, int nofbp) {/*{{{*/
    if (rearr == 'r') {
        if (elem_is_P_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_P_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_P_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_P_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_P_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_P_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (!elem_is_P_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_P_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_P_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_P_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp--;
        else if (!elem_is_P_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_P_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/

int nofbps_after_spsrt(permutation_t *p, char rearr, char type, int i, int j, int nofbp) {/*{{{*/
    if (rearr == 't' && type == 'p') {
        if (elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (!elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp--;
        else if (!elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[1]))
            nofbp++;
    } else if (rearr == 't' && type == 's') {
        if (elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp--;
        else if (elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 !elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[j]))
            nofbp++;

        if (elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
            !elem_is_PSRT_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp--;
        else if (!elem_is_PSRT_breakpoint(p, p->pi[j-1], p->pi[j]) &&
                 elem_is_PSRT_breakpoint(p, p->pi[p->size], p->pi[i]))
            nofbp++;
    } else if (rearr == 'r' && type == 'p') {
        if (elem_is_PSRT_breakpoint(p, p->pi[i], p->pi[i+1]) &&
            !elem_is_PSRT_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp--;

        else if (!elem_is_PSRT_breakpoint(p, p->pi[i], p->pi[i+1]) &&
                 elem_is_PSRT_breakpoint(p, -p->pi[1], p->pi[i+1]))
            nofbp++;
    } else {
        if (elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
            !elem_is_PSRT_breakpoint(p, p->pi[i-1], -p->pi[p->size]))
            nofbp--;

        else if (!elem_is_PSRT_breakpoint(p, p->pi[i-1], p->pi[i]) &&
                 elem_is_PSRT_breakpoint(p, p->pi[i-1], -p->pi[p->size]))
            nofbp++;
    }
    return nofbp;
}/*}}}*/


int nofbps_after(char prob[], char type, char rearr, 
                permutation_t *p, int i, int j, int k, int nofbp) {/*{{{*/
    if (!strcmp(prob, "pr")) {
        return nofbps_after_pr(p, i, nofbp);
    } else if (!strcmp(prob, "psr")) {
        return nofbps_after_psr(p, type, i, nofbp);
    } else if (!strcmp(prob, "pt")) {
        return nofbps_after_pt(p, i, j, nofbp);
    } else if (!strcmp(prob, "pst")) {
        return nofbps_after_pst(p, type, i, j, nofbp);
    } else if (!strcmp(prob, "prt")) {
        return nofbps_after_prt(p, rearr, i, j, nofbp);
    } else if (!strcmp(prob, "psrt")) {
        return nofbps_after_psrt(p, rearr, type, i, j, nofbp);
    } else if (!strcmp(prob, "spr")) {
        return nofbps_after_spr(p, i, nofbp);
    } else if (!strcmp(prob, "spsr")) {
        return nofbps_after_spsr(p, type, i, nofbp);
    } else if (!strcmp(prob, "sprt")) {
        return nofbps_after_sprt(p, rearr, i, j, nofbp);
    } else if (!strcmp(prob, "spsrt")) {
        return nofbps_after_spsrt(p, rearr, type, i, j, nofbp);
    } else {
        printf("ERRO! nofbps_after()\n");
        exit(0);
    }
}/*}}}*/


