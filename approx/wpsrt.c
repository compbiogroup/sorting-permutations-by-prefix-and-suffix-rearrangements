/***********************************************************
 * Created: Seg 17 Mar 2014 01:40:14 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "../util.h"
#include "../permutations.h"
#include "../rearrangements.h"
#include "../breakpoints.h"
#include "wprt.h"
#include "wsrt.h"

void alg_WPSRTg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight;

    while (!is_identity(p)) {
        type = 0;
        weight = pot(p->size * p->size, info->alpha);

        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            info->weight += pot(p->size, info->alpha);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
            continue;
        }

/* PREFIX TRANSPOSITIONS {{{*/
        j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            weight = pot(j-1, info->alpha);
            type = 1;
            p1 = i;
            p2 = j;
        }

        j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
            weight = pot(j-1, info->alpha);
            type = 1;
            p1 = i;
            p2 = j;
        }

        j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) &&
                is_UPSRT_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
            weight = pot(j-1, info->alpha);
            type = 1;
            p1 = i;
            p2 = j;
        }

        j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) &&
                is_UPSRT_breakpoint(p, j-1, j) && pot(j-1, info->alpha) < weight) {
            weight = pot(j-1, info->alpha);
            type = 1;
            p1 = i;
            p2 = j;
        }
/*}}}*/

/* SUFFIX TRANSPOSITIONS {{{*/
        i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
               is_UPSRT_breakpoint(p, i-1, i) && 
               is_UPSRT_breakpoint(p, j-1, j) && pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 2;
            p1 = i;
            p2 = j;
        }

        i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
               is_UPSRT_breakpoint(p, i-1, i) && 
               is_UPSRT_breakpoint(p, j-1, j) && pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 2;
            p1 = i;
            p2 = j;
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
               is_UPSRT_breakpoint(p, i-1, i) && 
               is_UPSRT_breakpoint(p, j-1, j) && pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 2;
            p1 = i;
            p2 = j;
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
               is_UPSRT_breakpoint(p, i-1, i) && 
               is_UPSRT_breakpoint(p, j-1, j) && pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 2;
            p1 = i;
            p2 = j;
        }
/*}}}*/

        i = p->inv_pi[inc_mod(p->pi[1], p->size)] - 1;
        if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i, i+1) &&
                        pot(i, info->alpha) < weight) {
            weight = pot(i, info->alpha);
            type = 3;
            p1 = i;
        }

        i = p->inv_pi[dec_mod(p->pi[1], p->size)] - 1;
        if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i, i+1) &&
                        pot(i, info->alpha) < weight) {
            weight = pot(i, info->alpha);
            type = 3;
            p1 = i;
        }

        for (i = 2; i <= p->size-1; i++) {
            if (!is_UPSRT_breakpoint(p, i-1, i))
                continue;

            j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
            if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                type = 1;
                p1 = i;
                p2 = j;
            }

            j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
            if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j) &&
                        pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                type = 1;
                p1 = i;
                p2 = j;
            }
        }

        j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
        if ((j == p->size+1 || (j > 2 && j < p->size+1 && is_UPSRT_breakpoint(p, j-1, j)))
                        && pot(j-1, info->alpha) < weight) {
            for (i = 2; i < j; i++) {
                if (is_UPSRT_breakpoint(p, i-1, i)) {
                    weight = pot(j-1, info->alpha);
                    type = 1;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
        if ((j == p->size+1 || (j > 2 && j < p->size+1 && is_UPSRT_breakpoint(p, j-1, j)))
                        && pot(j-1, info->alpha) < weight) {
            for (i = 2; i < j; i++) {
                if (is_UPSRT_breakpoint(p, i-1, i)) {
                    weight = pot(j-1, info->alpha);
                    type = 1;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)] + 1;
        if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i-1, i) &&
                        pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 4;
            p1 = i;
        }

        i = p->inv_pi[dec_mod(p->pi[p->size], p->size)] + 1;
        if (2 <= i && i < p->size && is_UPSRT_breakpoint(p, i-1, i) &&
                        pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size-i+1, info->alpha);
            type = 4;
            p1 = i;
        }

        for (i = 2; i <= p->size-1 && pot(p->size-i+1, info->alpha) < weight; i++) {
            if (!is_UPSRT_breakpoint(p, i-1, i))
                continue;

            j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
            if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
                weight = pot(p->size-i+1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j;
            }

            j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
            if (i < j && j <= p->size && is_UPSRT_breakpoint(p, j-1, j)) {
                weight = pot(p->size-i+1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j;
            }
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
        if ((i == 1 || (i > 1 && i <= p->size-1 && is_UPSRT_breakpoint(p, i-1, i)))
                        && pot(p->size-i+1, info->alpha) < weight) {
            for (j = i+1; j <= p->size; j++) {
                if (is_UPSRT_breakpoint(p, j-1, j)) {
                    weight = pot(p->size-i+1, info->alpha);
                    type = 2;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        i = p->inv_pi[p->pi[p->size] - 1];
        if ((i == 1 || (i > 0 && i <= p->size-1 && is_UPSRT_breakpoint(p, i-1, i)))
                        && pot(p->size-i+1, info->alpha) < weight) {
            for (j = i+1; j <= p->size; j++) {
                if (is_UPSRT_breakpoint(p, j-1, j)) {
                    weight = pot(p->size-i+1, info->alpha);
                    type = 2;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        if (type == 1) {
            op_prefix_transposition(p, p1, p2);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 't', 1, p1, p2);
        } else if (type == 2) {
            op_suffix_transposition(p, p1, p2);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 't', p1, p2, p->size+1);
        } else if (type == 3) {
            op_prefix_reversal(p, p1);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 'r', 1, p1, 0);
        } else if (type == 4) {
            op_suffix_reversal(p, p1);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 'r', p1, p->size, 0);
        } else {
            i = p->inv_pi[1];
            if (p->size == p->pi[i+1]) {
                op_prefix_transposition(p, i+1, p->size+1);
                print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            } else {
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            }
            info->weight += pot(p->size, info->alpha);
            info->nof_pt++;
            info->nmoves++;
        }
    }
}/*}}}*/





long double partitionWPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, k, median, x, y, z;
    long double weight = 0;
        
    median = prm_median(p, 1, p->size);

    x = 0;
    y = p->size + 1;
    while (x < p->size && p->pi[x+1] <= median) x++;
    while (y > 1 && p->pi[y-1] > median) y--;

    if (x == y-1)
        return 0;

    i = x+1;
    while (i < y && p->pi[i] > median) i++;
    z = i-1;
    while (i < y && p->pi[i] <= median) i++;

    if (i == y) {
        if (y-1 <= p->size-x) {
            op_prefix_transposition(p, z+1, y);
            print_verbose(p, verbose, 't', 1, z+1, y);
            info->nmoves++;
            return pot(y-1, info->alpha);
        } else {
            op_suffix_transposition(p, x+1, z+1);
            print_verbose(p, verbose, 't', x+1, z+1, p->size+1);
            info->nmoves++;
            return pot(p->size-x, info->alpha);
        }
    }

    weight += partitionWPRT(p, median, INC, median, verbose, info);
    weight += partitionWSRT(p, median+1, DEC, median, verbose, info);

    i = 1;
    while (p->pi[i] <= median && i <= p->size) i++;
    j = i;
    while (p->pi[i] > median && i <= p->size) i++;
    k = i-1;

    if (k < p->size) {
        op_suffix_reversal(p, j);
        print_verbose(p, verbose, 'r', j, p->size, 0);
        weight += pot(p->size-j+1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int median, n = p->size, i, j;
    long double weight = 0;

    if (is_sorted_interval(p, 1, n, INC))
        return 0;

    if (n <= 1)
        return 0;

    if (n == 2) {
        if (p->pi[1] > p->pi[2]) {
            op_prefix_reversal(p, 2);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 'r', 1, 2, 0);
            info->nmoves++;
        }
        return weight;
    }

    prm_separated(p, &i, &j, &n);

    if (n < p->size || (n == p->size && !(i == 0 && j == 1))) {
        weight += WPRT(p, i, INC, verbose, info);
        weight += WSRT(p, j, INC, verbose, info);
    } else {
        weight += partitionWPSRT(p, verbose, info);
        median = prm_median(p, 1, p->size);
        weight += WPRT(p, median, INC, verbose, info);
        weight += WSRT(p, median+1, INC, verbose, info);
    }

    return weight;
}/*}}}*/

void alg_WPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    if (!is_identity(p))
        info->weight = WPSRT(p, verbose, info);

    if (!is_identity(p))
        printf("ERROR! alg_WPSRT()\n");
}/*}}}*/







void alg_WSPSRTg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight = 0;

    while (!is_identity(p)) {
        type = 0;
        weight = pot(p->size * p->size, info->alpha);

        if (is_signed_reverse(p)) {
            op_prefix_reversal(p, p->size);
            info->weight += pot(p->size, info->alpha);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            continue;
        }

        j = prm_exists(p, dec_mod(p->pi[1], p->size));
        if (j > 0) {
            j++;
            i = prm_exists(p, dec_mod(p->pi[j], p->size));
            if (2 <= i+1 && i+1 < j && j <= p->size &&
                            is_PSRT_breakpoint(p, i, i+1) &&
                            is_PSRT_breakpoint(p, j-1, j)) {
                weight = pot(j-1, info->alpha);
                type = 1;
                p1 = i+1;
                p2 = j;
            }
        }

        i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
        if (i > 0) {
            j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
            if (j > 0 && i < j && j <= p->size &&
                            is_PSRT_breakpoint(p, i-1, i) &&
                            is_PSRT_breakpoint(p, j-1, j) &&
                            pot(p->size-i+1, info->alpha)) {
                weight = pot(p->size-i+1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j;
            }
        }

        i = prm_exists(p, -dec_mod(p->pi[1], p->size));
        if (i >= 2 && i <= p->size && is_PSRT_breakpoint(p, i-1, i) &&
                        pot(i-1, info->alpha) < weight) {
            weight = pot(i-1, info->alpha);
            type = 3;
            p1 = i-1;
        }

        for (i = 2; i <= p->size-1; i++) {
            if (!is_PSRT_breakpoint(p, i-1, i))
                continue;

            j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
            if (i < j && j <= p->size && is_PSRT_breakpoint(p, j-1, j) &&
                            pot(j-1, info->alpha) < weight) {
                weight = pot(j-1, info->alpha);
                type = 1;
                p1 = i;
                p2 = j;
            }
        }

        j = prm_exists(p, dec_mod(p->pi[1], p->size));
        if ((j == p->size || (j > 1 && j < p->size && is_PSRT_breakpoint(p, j, j+1)))
                        && pot(j-1, info->alpha) < weight) {
            j++;
            for (i = 2; i < j; i++) {
                if (is_PSRT_breakpoint(p, i-1, i)) {
                    weight = pot(j-1, info->alpha);
                    type = 1;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        i = prm_exists(p, -inc_mod(p->pi[p->size], p->size));
        if (i >= 1 && i <= p->size-1 && is_PSRT_breakpoint(p, i, i+1) &&
                pot(p->size-i, info->alpha) < weight) {
            weight = pot(p->size-i, info->alpha);
            type = 4;
            p1 = i+1;
        }

        for (i = 2; i <= p->size-1 && pot(p->size-i+1, info->alpha) < weight; i++) {
            if (!is_PSRT_breakpoint(p, i-1, i))
                continue;

            j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
            if (j > 2 && i < j && j <= p->size && is_PSRT_breakpoint(p, j-1, j)) {
                weight = pot(p->size-i+1, info->alpha);
                type = 2;
                p1 = i;
                p2 = j;
            }
        }

        i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
        if ((i == 1 || (i > 0 && i < p->size && is_PSRT_breakpoint(p, i-1, i)))
                        && pot(p->size-i+1, info->alpha) < weight) {
            for (j = i+1; j <= p->size; j++) {
                if (is_PSRT_breakpoint(p, j-1, j)) {
                    weight = pot(p->size-i+1, info->alpha);
                    type = 2;
                    p1 = i;
                    p2 = j;
                }
            }
        }

        if (type == 1) {
            op_prefix_transposition(p, p1, p2);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 't', 1, p1, p2);
        } else if (type == 2) {
            op_suffix_transposition(p, p1, p2);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 't', p1, p2, p->size+1);
        } else if (type == 3) {
            op_prefix_reversal(p, p1);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 'r', 1, p1, 0);
        } else if (type == 4) {
            op_suffix_reversal(p, p1);
            info->weight += weight;
            info->nmoves++;
            print_verbose(p, verbose, 'r', p1, p->size, 0);
        } else {
                if (nof_PSRT_breakpoints(p) != 0) printf("erro\n");
            i = prm_exists(p, 1);
            if (i > 0) {
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            } else {
                i = prm_exists(p, -p->size);
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            }
            info->nof_pt++;
            info->weight += pot(p->size, info->alpha);
            info->nmoves++;
        }

    }
}/*}}}*/




long double partitionWSPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, k, median, x, y, z;
    long double weight = 0;
        
    median = prm_median(p, 1, p->size);
    if (median % 2)
        median++;

    x = 0;
    y = p->size + 1;
    while (x < p->size && p->pi[x+1] <= median) x++;
    while (y > 1 && p->pi[y-1] > median) y--;

    if (x == y-1)
        return 0;

    i = x+1;
    while (i < y && p->pi[i] > median) i++;
    z = i-1;
    while (i < y && p->pi[i] <= median) i++;

    if (i == y) {
        if (y-1 <= p->size-x) {
            op_prefix_transposition(p, z+1, y);
            info->nmoves++;
            return pot(y-1, info->alpha);
        } else {
            op_suffix_transposition(p, x+1, z+1);
            info->nmoves++;
            return pot(p->size-x, info->alpha);
        }
    }

    weight += partitionWSPRT(p, median, INC, median, verbose, info);
    weight += partitionWSSRT(p, median+1, DEC, median, verbose, info);

    i = 1;
    while (p->pi[i] <= median && i <= p->size) i++;
    j = i;
    while (p->pi[i] > median && i <= p->size) i++;
    k = i-1;

    if (k < p->size) {
        op_suffix_reversal(p, j);
        weight += pot(p->size-j+1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int median, n = p->size, i, j;
    long double weight = 0;

    if (is_sorted_interval(p, 1, n, INC))
        return 0;

    if (n <= 1)
        return 0;

    if (n == 2) {
        if (p->pi[1] > p->pi[2]) {
            op_prefix_reversal(p, 2);
            weight += pot(2, info->alpha);
            info->nmoves++;
        }
        return weight;
    }

    prm_separated(p, &i, &j, &n);

    if (n < p->size || (n == p->size && !(i == 0 && j == 1))) {
        weight += WSPRT(p, i, INC, verbose, info);
        weight += WSSRT(p, j, INC, verbose, info);
    } else {
        weight += partitionWSPSRT(p, verbose, info);
        median = prm_median(p, 1, p->size);
        if (median % 2)
            median++;
        weight += WSPRT(p, median, INC, verbose, info);
        weight += WSSRT(p, median+1, INC, verbose, info);
    }

    return weight;
}/*}}}*/

void alg_WSPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    permutation_t up;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    if (!is_identity(p)) {
        info->weight = WSPSRT(&up, verbose, info) / 2;
        prm_back_from_image(p, &up);
    }

    destroy_permutation(&up);

    if (!is_identity(p))
        printf("ERROR! alg_WSPSRT()\n");
}/*}}}*/





