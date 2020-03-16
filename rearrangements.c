/***********************************************************
 * Created: Sáb 19 Jan 2013 11:39:13 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "permutations.h"
#include "breakpoints.h"
#include "util.h"



void op_reversal(permutation_t *p, int i, int j) {/*{{{*/
    int k, aux;

    if (p->sig == UNSIGNED) {
        for (k = i; k <= (j-i)/2 + i; k++) {
            aux = p->pi[k];
            p->pi[k] = p->pi[j - k + i];
            p->pi[j - k + i] = aux;
            p->inv_pi[p->pi[k]] = k;
            p->inv_pi[p->pi[j - k + i]] = j - k + i;
        }
    } else {
        for (k = i; k <= (j-i)/2 + i; k++) {
            prm_insert(p, -p->pi[k], k);
            prm_insert(p, -p->pi[j - k + i], j-k+i);
            prm_swap(p, k, j-k+i);
        }
        if ((j-i+1)%2)
            prm_insert(p, -p->pi[(j-i)/2+i], (j-i)/2+i);
    }
}/*}}}*/

void op_prefix_reversal(permutation_t *p, int j) {/*{{{*/
    int i, aux;

    if (p->sig == UNSIGNED) {
        for (i = 1; i <= j/2; i++) {
            aux = p->pi[i];
            p->pi[i] = p->pi[j - i + 1];
            p->pi[j - i + 1] = aux;
            p->inv_pi[p->pi[i]] = i;
            p->inv_pi[p->pi[j - i + 1]] = j - i + 1;
        }
    } else {
        for (i = 1; i <= j/2; i++) {
            prm_insert(p, -p->pi[i], i);
            prm_insert(p, -p->pi[j - i + 1], j-i+1);
            prm_swap(p, i, j-i+1);
        }
        if (j%2)
            prm_insert(p, -p->pi[j/2+1], j/2+1);
    }
}/*}}}*/

void op_suffix_reversal(permutation_t *p, int j) {/*{{{*/
    int i, aux;
    
    if (p->sig == UNSIGNED) {
        for (i = j; i <= (p->size - j)/2 + j; i++) {
            aux = p->pi[i];
            p->pi[i] = p->pi[p->size - i + j];
            p->pi[p->size - i + j] = aux;
            p->inv_pi[p->pi[i]] = i;
            p->inv_pi[p->pi[p->size - i + j]] = p->size - i + j;
        }
    } else {
        for (i = j; i <= (p->size - j)/2 + j; i++) {
            prm_insert(p, -p->pi[i], i);
            prm_insert(p, -p->pi[p->size - i + j], p->size-i+j);
            prm_swap(p, i, p->size-i+j);
        }
        if ((p->size-j+1)%2)
            prm_insert(p, -p->pi[(p->size-j)/2+j], (p->size-j)/2+j);
    }
}/*}}}*/

void op_transposition(permutation_t *p, int i, int j, int k) {/*{{{*/
    int t, l, seg1, seg2, *aux;

    aux = Malloc(sizeof(int) * (j - i + 1));
    seg1 = j - i;
    seg2 = k - j;

    if (seg2 < seg1) {
        for (t = i; t <= seg2 + i - 1; t++) {
            aux[t-i] = p->pi[t];
            p->pi[t] = p->pi[j + t - i];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
        l = 0;
        for (t = i + seg2; t < j; t++) {
            aux[t-i] = p->pi[t];
            p->pi[t] = aux[l++];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
        for (t = j; t < k; t++) {
            p->pi[t] = aux[l++];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
    } else {
        for (t = i; t < j; t++) {
            aux[t-i] = p->pi[t];
            p->pi[t] = p->pi[j + t - i];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
        for (t = j; t <= i + seg2 - 1; t++) {
            p->pi[t] = p->pi[t + seg1];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
        l = 0;
        for (t = i + seg2; t < k; t++) {
            p->pi[t] = aux[l++];
            p->inv_pi[abs(p->pi[t])] = (p->pi[t] < 0) ? -1*t : t;
        }
    }

    free(aux);
}/*}}}*/

void op_prefix_transposition(permutation_t *p, int j, int k) {/*{{{*/
    op_transposition(p, 1, j, k);
}/*}}}*/

void op_suffix_transposition(permutation_t *p, int i, int j) {/*{{{*/
    op_transposition(p, i, j, p->size + 1);
}/*}}}*/




double lower_bound_pr(permutation_t *p) {/*{{{*/
    return nof_UPR_breakpoints(p);
}/*}}}*/

double lower_bound_psr(permutation_t *p) {/*{{{*/
    return nof_UPSR_breakpoints(p);
}/*}}}*/

double lower_bound_pt(permutation_t *p) {/*{{{*/
    return nof_P_breakpoints(p) / 2.0;
}/*}}}*/

double lower_bound_pst(permutation_t *p) {/*{{{*/
    return nof_PS_breakpoints(p) / 2.0;
}/*}}}*/

double lower_bound_prt(permutation_t *p) {/*{{{*/
    return nof_UPR_breakpoints(p) / 2.0;
}/*}}}*/

double lower_bound_psrt(permutation_t *p) {/*{{{*/
    return nof_UPSRT_breakpoints(p) / 2.0;
}/*}}}*/


double lower_bound_spr(permutation_t *p) {/*{{{*/
    return nof_P_breakpoints(p);
}/*}}}*/

double lower_bound_spsr(permutation_t *p) {/*{{{*/
    return nof_PS_breakpoints(p);
}/*}}}*/

double lower_bound_sprt(permutation_t *p) {/*{{{*/
    return nof_P_breakpoints(p) / 2.0;
}/*}}}*/

double lower_bound_spsrt(permutation_t *p) {/*{{{*/
    return nof_PSRT_breakpoints(p) / 2.0;
}/*}}}*/


double lower_bound_wpr(permutation_t *p) {/*{{{*/
    int i = p->size;
    while (p->pi[i] == i && i > 1)
        i--;
    return pot(i, p->alpha);
}/*}}}*/

double lower_bound_wpt(permutation_t *p) {/*{{{*/
    int i = p->size;
    while (p->pi[i] == i && i > 1)
        i--;
    return pot(i, p->alpha);
}/*}}}*/

double lower_bound_wprt(permutation_t *p) {/*{{{*/
    int i = p->size;
    while (p->pi[i] == i && i > 1)
        i--;
    return pot(i, p->alpha);
}/*}}}*/

/* It only works for alpha = 1*/
double lower_bound_wpsr(permutation_t *p) {/*{{{*/
    int i, j, l, k, n = p->size, found;

    for (i = 1; i <= p->size; i++) {
        for (j = i; j <= p->size; j++) {
            l = i;
            while (p->pi[l] == l && l <= j) l++;
            if (l == j+1) {
                found = 1;
                for (k = 1; i <= i-1 && found; k++)
                    if (p->pi[k] >= i) found = 0;
                for (k = j+1; k <= p->size && found; k++)
                    if (p->pi[k] <= j) found = 0;
                if (found) {
                    if (p->size + i - j - 1 < n) {
                        n = p->size + i - j - 1;
                    }
                }
            }
        }
    }

    return pot(n, p->alpha);
}/*}}}*/

double lower_bound_wpst(permutation_t *p) {/*{{{*/
    int i, j, l, k, n = p->size, found;
    /* It only works for alpha = 1*/

    for (i = 1; i <= p->size; i++) {
        for (j = i; j <= p->size; j++) {
            l = i;
            while (p->pi[l] == l && l <= j) l++;
            if (l == j+1) {
                found = 1;
                for (k = 1; i <= i-1 && found; k++)
                    if (p->pi[k] >= i) found = 0;
                for (k = j+1; k <= p->size && found; k++)
                    if (p->pi[k] <= j) found = 0;
                if (found) {
                    if (p->size + i - j - 1 < n) {
                        n = p->size + i - j - 1;
                    }
                }
            }
        }
    }

    return pot(n, p->alpha);
}/*}}}*/

double lower_bound_wpsrt(permutation_t *p) {/*{{{*/
    int i, j, l, k, n = p->size, found;
    /* It only works for alpha = 1*/

    for (i = 1; i <= p->size; i++) {
        for (j = i; j <= p->size; j++) {
            l = i;
            while (p->pi[l] == l && l <= j) l++;
            if (l == j+1) {
                found = 1;
                for (k = 1; i <= i-1 && found; k++)
                    if (p->pi[k] >= i) found = 0;
                for (k = j+1; k <= p->size && found; k++)
                    if (p->pi[k] <= j) found = 0;
                if (found) {
                    if (p->size + i - j - 1 < n) {
                        n = p->size + i - j - 1;
                    }
                }
            }
        }
    }

    return pot(n, p->alpha);
}/*}}}*/


double lower_bound_wspr(permutation_t *p) {/*{{{*/
    int i = p->size;
    while (p->pi[i] == i && i > 1)
        i--;
    return pot(i, p->alpha);
}/*}}}*/

double lower_bound_wsprt(permutation_t *p) {/*{{{*/
    int i = p->size;
    while (p->pi[i] == i && i > 1)
        i--;
    return pot(i, p->alpha);
}/*}}}*/

double lower_bound_wspsr(permutation_t *p) {/*{{{*/
    int i, j, l, k, n = p->size, found;
    /* It only works for alpha = 1*/

    for (i = 1; i <= p->size; i++) {
        for (j = i; j <= p->size; j++) {
            l = i;
            while (p->pi[l] == l && l <= j) l++;
            if (l == j+1) {
                found = 1;
                for (k = 1; i <= i-1 && found; k++)
                    if (abs(p->pi[k]) >= i) found = 0;
                for (k = j+1; k <= p->size && found; k++)
                    if (abs(p->pi[k]) <= j) found = 0;
                if (found) {
                    if (p->size + i - j - 1 < n) {
                        n = p->size + i - j - 1;
                    }
                }
            }
        }
    }

    return pot(n, p->alpha);
}/*}}}*/

double lower_bound_wspsrt(permutation_t *p) {/*{{{*/
    int i, j, l, k, n = p->size, found;
    /* It only works for alpha = 1*/

    for (i = 1; i <= p->size; i++) {
        for (j = i; j <= p->size; j++) {
            l = i;
            while (p->pi[l] == l && l <= j) l++;
            if (l == j+1) {
                found = 1;
                for (k = 1; i <= i-1 && found; k++)
                    if (abs(p->pi[k]) >= i) found = 0;
                for (k = j+1; k <= p->size && found; k++)
                    if (abs(p->pi[k]) <= j) found = 0;
                if (found) {
                    if (p->size + i - j - 1 < n) {
                        n = p->size + i - j - 1;
                    }
                }
            }
        }
    }

    return pot(n, p->alpha);
}/*}}}*/


double lower_bound(char prob_name[], permutation_t *p) {/*{{{*/
    if (!strcmp(prob_name, "pr"))
        return lower_bound_pr(p);
    else if (!strcmp(prob_name, "psr"))
        return lower_bound_psr(p);

    else if (!strcmp(prob_name, "pt"))
        return lower_bound_pt(p);
    else if (!strcmp(prob_name, "pst"))
        return lower_bound_pst(p);

    else if (!strcmp(prob_name, "prt"))
        return lower_bound_prt(p);
    else if (!strcmp(prob_name, "psrt"))
        return lower_bound_psrt(p);

    else if (!strcmp(prob_name, "spr"))
        return lower_bound_spr(p);
    else if (!strcmp(prob_name, "spsr"))
        return lower_bound_spsr(p);

    else if (!strcmp(prob_name, "sprt"))
        return lower_bound_sprt(p);
    else if (!strcmp(prob_name, "spsrt"))
        return lower_bound_spsrt(p);

    else if (!strcmp(prob_name, "wpr"))
        return lower_bound_wpr(p);
    else if (!strcmp(prob_name, "wpt"))
        return lower_bound_wpt(p);
    else if (!strcmp(prob_name, "wpsr"))
        return lower_bound_wpsr(p);
    else if (!strcmp(prob_name, "wpst"))
        return lower_bound_wpst(p);
    else if (!strcmp(prob_name, "wprt"))
        return lower_bound_wprt(p);
    else if (!strcmp(prob_name, "wpsrt"))
        return lower_bound_wpsrt(p);

    else if (!strcmp(prob_name, "wspr"))
        return lower_bound_wspr(p);
    else if (!strcmp(prob_name, "wspsr"))
        return lower_bound_wspsr(p);
    else if (!strcmp(prob_name, "wsprt"))
        return lower_bound_wsprt(p);
    else if (!strcmp(prob_name, "wspsrt"))
        return lower_bound_wspsrt(p);

    else {
        printf("ERROR! lower_bound(): invalid prob_name\n");
        exit(0);
    }
}/*}}}*/




int prt_prohibited(permutation_t *p, char type, int i, int j) {/*{{{*/
    int ini_1, end_1, ini_pi1, end_pi1;

    find_extremes_of_strip(p, 1, &ini_1, &end_1);
    /* if the last strip does not start with 1, then it is not a problem */
    if (end_1 != p->size)
        return 0;

    /* so, the last strip starts at 1 */
    if (type == 't') {
        find_extremes_of_strip(p, p->pi[1], &ini_pi1, &end_pi1);

        if (p->pi[i-1] == p->size && p->pi[i] == 1 && end_pi1 != i-1)
            return 1;

        if (p->pi[j-1] == p->size && p->pi[j] == 1)
            return 1;

        if (p->pi[i] == 1 && p->pi[i-1] != p->size)
            return 1;

        if (j == p->size+1 && end_pi1 != i-1)
            return 1;

    } else {
        if (p->pi[i+1] == 1)
            return 1;
    }

    return 0;
}/*}}}*/

int sprt_prohibited(permutation_t *p, char type, int i, int j) {/*{{{*/
    int ini_1, end_1, ini_pi1, end_pi1;

    if (prm_exists(p, 1) != -1) {
        find_extremes_of_signed_strip(p, 1, &ini_1, &end_1);
        /* if the last strip does not start with 1, then it is not a problem */
        if (end_1 != p->size)
            return 0;
    } else {
        return 0;
    }

    /* so, the last strip starts at 1 */
    if (type == 't') {
        find_extremes_of_signed_strip(p, p->pi[1], &ini_pi1, &end_pi1);

        if (p->pi[i-1] == p->size && p->pi[i] == 1 && end_pi1 != i-1)
            return 1;

        if (p->pi[j-1] == p->size && p->pi[j] == 1)
            return 1;

        if (p->pi[i] == 1 && p->pi[i-1] != p->size)
            return 1;

        if (j == p->size+1 && end_pi1 != i-1)
            return 1;

    } else {
        if (p->pi[i+1] == 1)
            return 1;
    }

    return 0;
}/*}}}*/

