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


long double partitionWST(permutation_t *p, int ini, int divisor, int verbose, info_t *info) {/*{{{*/
    /* Deixar elementos maiores que median no comeco da permutacao */
    int i, j, k, x, y, z, half, n = p->size;
    long double weight = 0;
        
    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if (p->pi[n-1] <= divisor && p->pi[n] > divisor) {
            op_suffix_transposition(p, n-1, n);
            print_verbose(p, verbose, 't', n-1, n, n+1);
            info->nmoves++;
            return pot(2, info->alpha);
        }
        return 0;
    }

    x = ini - 1;
    y = n + 1;
    while (x < n && p->pi[x+1] > divisor) x++;
    while (y > ini && p->pi[y-1] <= divisor) y--;

    if (x == y-1)
        return 0;

    i = x+1;
    while (i < y && p->pi[i] <= divisor) i++;
    z = i-1;
    while (i < y && p->pi[i] > divisor) i++;
    if (i == y) {
        op_suffix_transposition(p, x+1, z+1);
        print_verbose(p, verbose, 't', x+1, z+1, n+1);
        info->nmoves++;
        return pot(n-x, info->alpha);
    }

    ini = x+1;

    half = (n-ini+1) / 2.0;
    if ((n-ini+1) % 2)
        half++;
    half += ini-1;

    weight = partitionWST(p, half+1, divisor, verbose, info);
    op_suffix_transposition(p, ini, half+1);
    print_verbose(p, verbose, 't', ini, half+1, n+1);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += partitionWST(p, n-half+ini, divisor, verbose, info);

    i = ini;
    while (i <= n && p->pi[i] > divisor) i++;
    j = i;
    while (i <= n && p->pi[i] <= divisor) i++;
    k = i-1;

    if (j <= k && k < n) {
        op_suffix_transposition(p, j, k+1);
        print_verbose(p, verbose, 't', j, k+1, n+1);
        weight += pot(n-j+1, info->alpha);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WST(permutation_t *p, int ini, int verbose, info_t *info) {/*{{{*/
    int part, median, n = p->size;
    long double weight = 0;

    if (is_sorted_interval(p, ini, n, INC))
        return 0;

    if (ini >= n)
        return 0;

    if (ini == n-1) {
        if (p->pi[n-1] > p->pi[n]) {
            op_suffix_transposition(p, n-1, n);
            weight += pot(2, info->alpha);
            print_verbose(p, verbose, 't', n-1, n, n+1);
            info->nmoves++;
        }
        return weight;
    }

    median = prm_median(p, ini, n);

    weight += partitionWST(p, ini, median, verbose, info);

    part = ini;
    while (p->pi[part] > median)
        part++;
    part--;

    weight += WST(p, part+1, verbose, info);
    op_suffix_transposition(p, ini, part+1);
    print_verbose(p, verbose, 't', ini, part+1, n+1);
    weight += pot(n-ini+1, info->alpha);
    info->nmoves++;
    weight += WST(p, n - part + ini, verbose, info);

    return weight;
}/*}}}*/

void alg_WST(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int k = 1;

    while (k == p->pi[k] && k <= p->size)
        k++;

    if (k != p->size+1) {
        info->weight = WST(p, k, verbose, info);
        info->nof_st = info->nmoves;
    }

    if (!is_identity(p))
        printf("ERROR! alg_WST()\n");
}/*}}}*/

