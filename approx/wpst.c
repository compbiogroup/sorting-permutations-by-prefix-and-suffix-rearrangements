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
#include "wpt.h"
#include "wst.h"


void alg_WPSTg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, p1, p2, type;
    long double weight;

    while (!is_identity(p)) {
        weight = pot(p->size * p->size, info->alpha);

        j = p->inv_pi[p->pi[1] - 1] + 1;
        i = p->inv_pi[p->pi[j] - 1] + 1;
        if (2 <= i && i < j && j <= p->size) {
            p1 = i;
            p2 = j;
            weight = pot(j-1, info->alpha);
            type = 1;
        }

        i = p->inv_pi[p->pi[p->size] + 1];
        j = p->inv_pi[p->pi[i - 1] + 1];
        if (2 <= i && i < j && j <= p->size && pot(p->size-i+1, info->alpha) < weight) {
            p1 = i;
            p2 = j;
            weight = pot(p->size - i + 1, info->alpha);
            type = 2;
        }

        for (i = 2; i <= p->size-1; i++) {
            if (!is_PS_breakpoint(p, i-1, i))
                continue;

            j = p->inv_pi[p->pi[i-1] + 1];
            if (i < j && j <= p->size) {
                if (pot(j-1, info->alpha) < weight) {
                    p1 = i;
                    p2 = j;
                    weight = pot(j-1, info->alpha);
                    type = 1;
                }
                if (pot(p->size-i+1, info->alpha) < weight) {
                    p1 = i;
                    p2 = j;
                    weight = pot(p->size-i+1, info->alpha);
                    type = 2;
                }
            }
        }

        j = p->inv_pi[p->pi[1] - 1] + 1;
        for (i = 2; i < j; i++) {
            if (is_PS_breakpoint(p, i-1, i) && pot(j-1, info->alpha) < weight) {
                p1 = i;
                p2 = j;
                weight = pot(j-1, info->alpha);
                type = 1;
            }
        }

        i = p->inv_pi[p->pi[p->size] + 1];
        for (j = i+1; j <= p->size; j++) {
            if (is_PS_breakpoint(p, j-1, j) && pot(p->size-i+1, info->alpha) < weight) {
                p1 = i;
                p2 = j;
                weight = pot(p->size-i+1, info->alpha);
                type = 2;
            }
        }

        if (type == 1) {
            op_prefix_transposition(p, p1, p2);
            print_verbose(p, verbose, 't', 1, p1, p2);
        } else {
            op_suffix_transposition(p, p1, p2);
            print_verbose(p, verbose, 't', p1, p2, p->size+1);
        }
        info->nmoves++;
        info->weight += weight;
    }
}/*}}}*/




long double partitionWPST(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, k, l, median, x, y, z;
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

    weight += partitionWPT(p, median, median, verbose, info);
    weight += partitionWST(p, median+1, median, verbose, info);

    i = 1;
    while (p->pi[i] > median && i <= p->size) i++;
    j = i;
    while (p->pi[i] <= median && i <= p->size) i++;
    k = i-1;
    while (p->pi[i] > median && i <= p->size) i++;
    l = i-1;

    if (j > 1 && j <= k) {
        op_prefix_transposition(p, j, k+1);
        print_verbose(p, verbose, 't', 1, j, k+1);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }
    if (l <= p->size && k <= l) {
        op_suffix_transposition(p, k-j+2, l+1);
        print_verbose(p, verbose, 't', k-j+2, l+1, p->size+1);
        weight += pot(p->size-k+j-1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

int WPST(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    long double weight = 0;
    int median, n = p->size, i, j;

    n = p->size;

    if (is_sorted_interval(p, 1, n, INC))
        return 0;

    if (n <= 1)
        return 0;

    if (n == 2) {
        if (p->pi[1] > p->pi[2]) {
            op_prefix_transposition(p, 2, 3);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 't', 1, 2, 3);
            info->nmoves++;
        }
        return weight;
    }

    prm_separated(p, &i, &j, &n);

    if (n < p->size || (n == p->size && !(i == 0 && j == 1))) {
        weight += WPT(p, i, verbose, info);
        weight += WST(p, j, verbose, info);
    } else {
        weight += partitionWPST(p, verbose, info);
        median = prm_median(p, 1, p->size);
        weight += WPT(p, median, verbose, info);
        weight += WST(p, median+1, verbose, info);
    }

    return weight;
}/*}}}*/

void alg_WPST(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    if (!is_identity(p))
        info->weight = WPST(p, verbose, info);

    if (!is_identity(p))
        printf("ERROR! alg_WPST()\n");
}/*}}}*/

