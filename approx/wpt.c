/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
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


void alg_WPTg(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int i, j, p1, p2, k = p->size;
    long double weight = 0;

    while (k == p->pi[k])
        k--;

    while (!is_identity(p)) {
        weight = pot(p->size * p->size, info->alpha);

        j = p->inv_pi[p->pi[1] - 1] + 1;
        i = p->inv_pi[p->pi[j] - 1] + 1;
        if (1 < i && i < j) {
            p1 = i;
            p2 = j;
            weight = pot(j-1, info->alpha);
        }
        
        for (i = 2; i <= p->size; i++) {
            if (!is_P_breakpoint(p, i-1, i))
                continue;

            j = p->inv_pi[p->pi[i-1]+1];
            if (i < j && pot(j-1, info->alpha) < weight) {
                p1 = i;
                p2 = j;
                weight = pot(j-1, info->alpha);
            }
        }

        j = p->inv_pi[p->pi[1] - 1] + 1;
        for (i = 2; i < j; i++) {
            if (is_P_breakpoint(p, i-1, i) && pot(j-1, info->alpha) < weight) {
                p1 = i;
                p2 = j;
                weight = pot(j-1, info->alpha);
            }
        }

        op_prefix_transposition(p, p1, p2);
        print_verbose(p, verbose, 't', 1, p1, p2);
        info->nmoves++;
        info->weight += weight;
    }
}/*}}}*/



long double partitionWPT(permutation_t *p, int end, int divisor, int verbose, info_t *info) {/*{{{*/
    /* Deixar elementos maiores que median no comeco da permutacao */
    int i, j, k, x, y, z, half;
    long double weight = 0;
        
    if (end <= 1)
        return 0;

    if (end == 2) {
        if (p->pi[1] <= divisor && p->pi[2] > divisor) {
            op_prefix_transposition(p, 2, 3);
            print_verbose(p, verbose, 't', 1, 2, 3);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = 0;
    y = end + 1;
    while (x < end && p->pi[x+1] > divisor) x++;
    while (y > 1 && p->pi[y-1] <= divisor) y--;

    if (x == y-1)
        return 0;

    i = x+1;
    while (i < y && p->pi[i] <= divisor) i++;
    z = i-1;
    while (i < y && p->pi[i] > divisor) i++;
    if (i == y) {
        op_prefix_transposition(p, z+1, y);
        print_verbose(p, verbose, 't', 1, z+1, y);
        info->nmoves++;
        return pot(y-1, info->alpha);
    }

    end = y-1;

    half = end / 2.0;
    if (end % 2)
        half++;

    weight = partitionWPT(p, half, divisor, verbose, info);
    op_prefix_transposition(p, half+1, end+1);
    print_verbose(p, verbose, 't', 1, half+1, end+1);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += partitionWPT(p, end - half, divisor, verbose, info);

    i = 1;
    while (i <= end && p->pi[i] > divisor) i++;
    while (i <= end && p->pi[i] <= divisor) i++;
    j = i;
    while (i <= end && p->pi[i] > divisor) i++;
    k = i-1;

    if (j <= k) {
        op_prefix_transposition(p, j, k+1);
        print_verbose(p, verbose, 't', 1, j, k+1);
        weight += pot(k, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WPT(permutation_t *p, int end, int verbose, info_t *info) {/*{{{*/
    int part, median;
    long double weight = 0;

    if (is_sorted_interval(p, 1, end, INC)) 
        return 0;

    if (end <= 1)
        return 0;

    if (end == 2) {
        if (p->pi[1] > p->pi[2]) {
            op_prefix_transposition(p, 2, 3);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 't', 1, 2, 3);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, 1, end);

    weight += partitionWPT(p, end, median, verbose, info);

    part = 1;
    while (part <= end && p->pi[part] > median)
        part++;
    part--;

    weight += WPT(p, part, verbose, info);
    op_prefix_transposition(p, part+1, end+1);
    print_verbose(p, verbose, 't', 1, part+1, end+1);
    weight += pot(end, info->alpha);
    info->nmoves++;
    weight += WPT(p, end - part, verbose, info);

    return weight;
}/*}}}*/

void alg_WPT(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int k = p->size;

    while (k == p->pi[k] && k >= 1)
        k--;

    if (k > 0) {
        info->weight = WPT(p, k, verbose, info);
        info->nof_pt = info->nmoves;
    }

    if (!is_identity(p)) 
        printf("ERROR! alg_WPT()\n");
}/*}}}*/

