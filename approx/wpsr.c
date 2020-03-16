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
#include "wpr.h"
#include "wsr.h"


void wpsr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, l, *lengths, new_n, nofbps;

    nofbps = nof_UPSR_breakpoints(p);
    lengths = malloc(sizeof(int) * (nofbps + 3));
    lengths[0] = 0;
    lengths[nofbps + 2] = 0;
    new_n = 0;
    l = 2;

    if (nofbps % 2) {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]-1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = 1; i <= nofbps+1; i += 2) {
            op_suffix_reversal(p, lengths[i]+1);
            print_verbose(p, verbose, 'r', lengths[i]+1, p->size, 0);
            op_prefix_reversal(p, p->size - lengths[i+1]);
            print_verbose(p, verbose, 'r', 1, p->size - lengths[i+1], 0);
            info->nmoves += 2;
            info->nof_sr++;
            info->nof_pr++;
            info->weight += pot(p->size-lengths[i], info->alpha) +
                    pot(p->size-lengths[i+1], info->alpha);
        }
    } else {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = nofbps+1; i >= 2; i -= 2) {
            op_prefix_reversal(p, p->size - lengths[i]);
            print_verbose(p, verbose, 'r', 1, p->size-lengths[i], 0);
            op_suffix_reversal(p, lengths[i-1]+1);
            print_verbose(p, verbose, 'r', lengths[i-1]+1, p->size, 0);
            info->nmoves += 2;
            info->nof_pr++;
            info->nof_sr++;
            info->weight += pot(p->size-lengths[i], info->alpha) +
                    pot(p->size-lengths[i-1], info->alpha);
        }
        op_prefix_reversal(p, p->size - lengths[1]);
        print_verbose(p, verbose, 'r', 1, p->size-lengths[1], 0);
        info->nmoves++;
        info->nof_pr++;
        info->weight += pot(p->size-lengths[1], info->alpha);
    }

    free(lengths);
}/*}}}*/

void alg_WPSRg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2;
    long double weight = 0;

    while (!is_identity(p)) {
        type = 0;
        weight = pot(p->size * p->size, info->alpha);

        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            info->weight += pot(p->size, info->alpha);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            continue;
        }

        /* Calcula o menor weight de remover um breakpoint com uma operação */
        i = p->inv_pi[p->pi[1]+1] - 1;
        if (i > 1 && i < p->size && is_UPSR_breakpoint(p, i, i+1)) {
            weight = pot(i, info->alpha);
            p1 = i;
            type = 1;
        }

        i = p->inv_pi[p->pi[1]-1] - 1;
        if (i > 1 && i < p->size && is_UPSR_breakpoint(p, i, i+1) && 
                        pot(i, info->alpha) < weight) {
            weight = pot(i, info->alpha);
            p1 = i;
            type = 1;
        }

        i = p->inv_pi[p->pi[p->size]+1] + 1;
        if (i > 1 && i < p->size && is_UPSR_breakpoint(p, i-1, i) &&
                        pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size - i + 1, info->alpha);
            p1 = i;
            type = 2;
        }

        i = p->inv_pi[p->pi[p->size]-1] + 1;
        if (i > 1 && i < p->size && is_UPSR_breakpoint(p, i-1, i) &&
                        pot(p->size-i+1, info->alpha) < weight) {
            weight = pot(p->size - i + 1, info->alpha);
            p1 = i;
            type = 2;
        }

        /* Calcula o menor weight de remover um breakpoint com duas operações */
        for (i = 1; i <= p->size-1; i++) {
            if (!is_UPSR_breakpoint(p, i, i+1))
                continue;

            j = p->inv_pi[p->pi[i] + 1];
            if (j > i+1 && j <= p->size && is_UPSR_breakpoint(p, j, j+1) &&
                            pot(j, info->alpha)+pot(j-i, info->alpha) < weight) {
                weight = pot(j, info->alpha) + pot(j - i, info->alpha);
                type = 3;
                p1 = j;
                p2 = j - i;
            }

            j = p->inv_pi[p->pi[i] - 1];
            if (j > i+1 && j <= p->size && is_UPSR_breakpoint(p, j, j+1) &&
                            pot(j, info->alpha)+pot(j-i, info->alpha) < weight) {
                weight = pot(j, info->alpha) + pot(j - i, info->alpha);
                type = 3;
                p1 = j;
                p2 = j - i;
            }

            if (i > 1) {
                j = p->inv_pi[p->pi[i] + 1];
                if (j > i+1 && j <= p->size && is_UPSR_breakpoint(p, j-1, j) &&
                                pot(i, info->alpha)+pot(j-1, info->alpha) < weight) {
                    weight = pot(i, info->alpha) + pot(j-1, info->alpha);
                    type = 3;
                    p1 = i;
                    p2 = j - 1;
                }

                j = p->inv_pi[p->pi[i] - 1];
                if (j > i+1 && j <= p->size && is_UPSR_breakpoint(p, j-1, j) &&
                                pot(i, info->alpha)+pot(j-1, info->alpha) < weight) {
                    weight = pot(i, info->alpha) + pot(j-1, info->alpha);
                    type = 3;
                    p1 = i;
                    p2 = j - 1;
                }
            }
        }

        for (j = p->size; j >= 2; j--) {
            if (!is_UPSR_breakpoint(p, j-1, j))
                continue;

            i = p->inv_pi[p->pi[j] + 1];
            if (j > i+1 && i > 0 && is_UPSR_breakpoint(p, i-1, i) &&
                      pot((j-i)+1, info->alpha)+pot(p->size-i+1, info->alpha) < weight) {
                weight = pot((j-i)+1, info->alpha) + pot(p->size-i+1, info->alpha);
                type = 4;
                p1 = i;
                p2 = p->size-(j-i)+1;
            }

            i = p->inv_pi[p->pi[j] - 1];
            if (j > i+1 && i > 0 && is_UPSR_breakpoint(p, i-1, i) &&
                      pot((j-i)+1, info->alpha)+pot(p->size-i+1, info->alpha) < weight) {
                weight = pot((j-i)+1, info->alpha) + pot(p->size-i+1, info->alpha);
                type = 4;
                p1 = i;
                p2 = p->size-(j-i)+1;
            }

            if (j < p->size) {
                i = p->inv_pi[p->pi[j] + 1];
                if (j > i+1 && i > 0 && is_UPSR_breakpoint(p, i, i+1) &&
                        pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha) < weight) {
                    weight = pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha);
                    type = 4;
                    p1 = j;
                    p2 = i+1;
                }

                i = p->inv_pi[p->pi[j] - 1];
                if (j > i+1 && i > 0 && is_UPSR_breakpoint(p, i, i+1) &&
                        pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha) < weight) {
                    weight = pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha);
                    type = 4;
                    p1 = j;
                    p2 = i+1;
                }
            }
        }

        if (type == 1) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            info->nmoves++;
            info->nof_pr++;
            info->weight += weight;
        } else if (type == 2) {
            op_suffix_reversal(p, p1);
            print_verbose(p, verbose, 'r', p1, p->size, 0);
            info->nmoves++;
            info->nof_sr++;
            info->weight += weight;
        } else if (type == 3) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            op_prefix_reversal(p, p2);
            print_verbose(p, verbose, 'r', 1, p2, 0);
            info->nmoves += 2;
            info->nof_pr += 2;
            info->weight += weight;
        } else if (type == 4) {
            op_suffix_reversal(p, p1);
            print_verbose(p, verbose, 'r', p1, p->size, 0);
            op_suffix_reversal(p, p2);
            print_verbose(p, verbose, 'r', p2, p->size, 0);
            info->nmoves += 2;
            info->nof_sr += 2;
            info->weight += weight;
        } else {
            wpsr_star_perm(p, verbose, info);
        }
    }
}/*}}}*/



long double partitionWPSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, x, y, z, k, median, cp, cs;
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
        cp = z + y - 1;
        if (x == 0) cp -= z;
        cs = p->size - z + p->size - x;
        if (y == p->size + 1) cs -= p->size-z;

        if (cp <= cs) {
            if (x > 0) {
                op_prefix_reversal(p, z);
                print_verbose(p, verbose, 'r', 1, z, 0);
                weight += pot(z, info->alpha);
                info->nmoves++;
            }
            op_prefix_reversal(p, y-1);
            print_verbose(p, verbose, 'r', 1, y-1, 0);
            weight += pot(y-1, info->alpha);
            info->nmoves++;
        } else {
            if (y < p->size+1) {
                op_suffix_reversal(p, z+1);
                print_verbose(p, verbose, 'r', z+1, p->size, 0);
                weight += pot(p->size-z, info->alpha);
                info->nmoves++;
            }
            op_suffix_reversal(p, x+1);
            print_verbose(p, verbose, 'r', x+1, p->size, 0);
            weight += pot(p->size-x, info->alpha);
            info->nmoves++;
        }
        return weight;
    }

    weight += partitionWPR(p, median, INC, median, verbose, info);
    weight += partitionWSR(p, median+1, DEC, median, verbose, info);

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

long double WPSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
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
        weight += WPR(p, i, INC, verbose, info);
        weight += WSR(p, j, INC, verbose, info);
    } else {
        weight += partitionWPSR(p, verbose, info);
        median = prm_median(p, 1, p->size);
        weight += WPR(p, median, INC, verbose, info);
        weight += WSR(p, median+1, INC, verbose, info);
    }

    return weight;
}/*}}}*/

void alg_WPSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    if (!is_identity(p))
        info->weight = WPSR(p, verbose, info);

    if (!is_identity(p))
        printf("ERROR! alg_WPSR()\n");
}/*}}}*/








void wspsr_star_perm(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, l, *lengths, new_n;
    int nofbps = nof_PS_breakpoints(p);

    lengths = malloc(sizeof(int) * (nofbps + 3));
    lengths[0] = 0;
    lengths[nofbps + 2] = 0;
    new_n = 0;
    l = 2;

    if (nofbps % 2) {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = 1; i <= nofbps+1; i += 2) {
            op_suffix_reversal(p, lengths[i]+1);
            print_verbose(p, verbose, 'r', lengths[i]+1, p->size, 0);
            op_prefix_reversal(p, p->size - lengths[i+1]);
            print_verbose(p, verbose, 'r', 1, p->size - lengths[i+1], 0);
            info->nmoves += 2;
            info->nof_sr++;
            info->nof_pr++;
            info->weight += pot(p->size-lengths[i], info->alpha) +
                    pot(p->size-lengths[i+1], info->alpha);
        }
    } else {
        for (i = 1; i <= nofbps+1; i++) {
            for (; p->pi[l] == p->pi[l-1]+1; l++) ;
            lengths[i] = l-1 - new_n;
            new_n += lengths[i];
            l++;
        }

        for (i = nofbps+1; i >= 2; i -= 2) {
            op_prefix_reversal(p, p->size - lengths[i]);
            print_verbose(p, verbose, 'r', 1, p->size-lengths[i], 0);
            op_suffix_reversal(p, lengths[i-1]+1);
            print_verbose(p, verbose, 'r', lengths[i-1]+1, p->size, 0);
            info->nmoves += 2;
            info->nof_pr++;
            info->nof_sr++;
            info->weight += pot(p->size-lengths[i], info->alpha) +
                    pot(p->size-lengths[i-1], info->alpha);
        }
        op_prefix_reversal(p, p->size - lengths[1]);
        print_verbose(p, verbose, 'r', 1, p->size-lengths[1], 0);
        info->nmoves++;
        info->nof_pr++;
        info->weight += pot(p->size-lengths[1], info->alpha);
    }

    free(lengths);
}/*}}}*/

void alg_WSPSRg(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, type, p1, p2, bps;
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

        i = prm_exists(p, -p->pi[1]+1);
        if (i >= 2 && i <= p->size) {
            weight = pot(i-1, info->alpha);
            type = 1;
            p1 = i-1;
        }

        i = prm_exists(p, -p->pi[p->size]-1);
        if (i >= 1 && pot(p->size-i, info->alpha) < weight) {
            weight = pot(p->size-i, info->alpha);
            p1 = i+1;
            type = 2;
        }

        for (i = 1; i <= p->size-1; i++) {
            if (!is_PS_breakpoint(p, i, i+1))
                continue;

            j = prm_exists(p, -p->pi[i]-1);
            if (j > 0 && i < j && j <= p->size && 
                        pot(j, info->alpha)+pot(j-i, info->alpha) < weight) {
                type = 3;
                weight = pot(j, info->alpha)+pot(j-i, info->alpha);
                p1 = j;
                p2 = j-i;
            }

            j = prm_exists(p, p->pi[i]+1);
            if (j > 0 && i+1 < j && j <= p->size && 
                        pot(i, info->alpha)+pot(j-1, info->alpha) < weight) {
                type = 3;
                weight = pot(i, info->alpha)+pot(j-1, info->alpha);
                p1 = i;
                p2 = j-1;
            }
        }

        for (j = p->size; j >= 3; j--) {
            if (!is_PS_breakpoint(p, j-1, j))
                continue;

            i = prm_exists(p, -p->pi[j]+1);
            if (i >= 1 && i < j && 
                    pot(p->size-i+1, info->alpha)+pot(j-i+1, info->alpha) < weight) {
                type = 4;
                weight = pot(p->size-i+1, info->alpha)+pot(j-i+1, info->alpha);
                p1 = i;
                p2 = p->size+1-(j-i);
            }

            i = prm_exists(p, p->pi[j]-1);
            if (i >= 1 && i+1 < j &&
               pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha) < weight) {
                type = 4;
                weight = pot(p->size-j+1, info->alpha)+pot(p->size-(i+1)+1, info->alpha);
                p1 = j;
                p2 = i+1;
            }
        }

        if (type == 1) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            info->nmoves++;
            info->nof_pr++;
            info->weight += weight;
        } else if (type == 2) {
            op_suffix_reversal(p, p1);
            print_verbose(p, verbose, 'r', p1, p->size, 0);
            info->nmoves++;
            info->nof_sr++;
            info->weight += weight;
        } else if (type == 3) {
            op_prefix_reversal(p, p1);
            print_verbose(p, verbose, 'r', 1, p1, 0);
            op_prefix_reversal(p, p2);
            print_verbose(p, verbose, 'r', 1, p2, 0);
            info->nmoves += 2;
            info->nof_pr += 2;
            info->weight += weight;
        } else if (type == 4) {
            op_suffix_reversal(p, p1);
            print_verbose(p, verbose, 'r', p1, p->size, 0);
            op_suffix_reversal(p, p2);
            print_verbose(p, verbose, 'r', p2, p->size, 0);
            info->nmoves += 2;
            info->nof_sr += 2;
            info->weight += weight;
        } else {
            bps = nof_PS_breakpoints(p);
            if ((bps % 2 == 1 && prm_exists(p, 1) &&
                 prm_exists(p, 1) > prm_exists(p, p->size)) || 
                        (bps % 2 == 0 && prm_exists(p, -1) &&
                         prm_exists(p, -1) < prm_exists(p, -p->size))) {
                op_prefix_reversal(p, p->size);
                print_verbose(p, verbose, 'r', 1, p->size, 0);
                info->nmoves++;
                info->nof_pr++;
                info->weight += pot(p->size, info->alpha);
            }
            wspsr_star_perm(p, verbose, info);
        }
    }
}/*}}}*/



long double partitionWSPSR(permutation_t *p, permutation_t *sp, int verbose, info_t *info) {/*{{{*/
    int i, j, k, median, x, y, z, cp, cs;
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
        cp = z/2 + (y-1)/2;
        if (x == 0) cp -= z/2;
        cs = (p->size-z)/2 + (p->size-x)/2;
        if (y == p->size + 1) cs -= (p->size-z)/2;

        if (cp <= cs) {
            if (x > 0) {
                op_prefix_reversal(p, z);
                op_prefix_reversal(sp, z/2);
                info->weight += pot(z/2, info->alpha);
                print_verbose(sp, verbose, 'r', 1, z/2, 0);
                info->nmoves++;
            }
            op_prefix_reversal(p, y-1);
            op_prefix_reversal(sp, (y-1)/2);
            info->weight += pot((y-1)/2, info->alpha);
            print_verbose(sp, verbose, 'r', 1, (y-1)/2, 0);
            info->nmoves++;
        } else {
            if (y < p->size+1) {
                op_suffix_reversal(p, z+1);
                i = (z+1)/2 + 1;
                op_suffix_reversal(sp, i);
                info->weight += pot((p->size-z)/2, info->alpha);
                print_verbose(sp, verbose, 'r', i, sp->size, 0);
                info->nmoves++;
            }
            op_suffix_reversal(p, x+1);
            i = (x+1)/2 + 1;
            op_suffix_reversal(sp, i);
            info->weight += pot((p->size-x)/2, info->alpha);
            print_verbose(sp, verbose, 'r', i, sp->size, 0);
            info->nmoves++;
        }
        return weight;
    }

    info->weight += partitionWSPR(p, sp, median, INC, median, verbose, info);
    info->weight += partitionWSSR(p, sp, median+1, DEC, median, verbose, info);

    i = 1;
    while (p->pi[i] <= median && i <= p->size) i++;
    j = i;
    while (p->pi[i] > median && i <= p->size) i++;
    k = i-1;

    if (k < p->size) {
        op_suffix_reversal(p, j);
        k = j/2 + 1;
        op_suffix_reversal(sp, k);
        info->weight += pot((p->size-j+1)/2, info->alpha);
        print_verbose(sp, verbose, 'r', k, p->size, 0);
        info->nmoves++;
    }

    return weight;
}/*}}}*/

long double WSPSR(permutation_t *p, permutation_t *sp, int verbose, info_t *info) {/*{{{*/
    int median, n = p->size, i, j;
    long double weight = 0;

    if (is_sorted_interval(p, 1, n, INC))
        return 0;

    if (n <= 1)
        return 0;

    if (n == 2) {
        if (p->pi[1] > p->pi[2]) {
            op_prefix_reversal(p, 2);
            op_prefix_reversal(sp, 1);
            info->weight += pot(1, info->alpha);
            print_verbose(sp, verbose, 'r', 1, 1, 0);
            info->nmoves++;
        }
        return weight;
    }

    prm_separated(p, &i, &j, &n);

    if (n < p->size || (n == p->size && !(i == 0 && j == 1))) {
        info->weight += WSPR(p, sp, i, INC, verbose, info);
        info->weight += WSSR(p, sp, j, INC, verbose, info);
    } else {
        info->weight += partitionWSPSR(p, sp, verbose, info);
        median = prm_median(p, 1, p->size);
        if (median % 2)
            median++;
        info->weight += WSPR(p, sp, median, INC, verbose, info);
        info->weight += WSSR(p, sp, median+1, INC, verbose, info);
    }

    return weight;
}/*}}}*/

void alg_WSPSR(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    permutation_t up;
    long double a;

    create_permutation(&up, 2*p->size, UNSIGNED);
    prm_image(p, &up);

    if (!is_identity(p))
        a = WSPSR(&up, p, verbose, info);
    printf("a:%Lf, %Lf\n", a, info->weight);

    if (!is_identity(p))
        printf("ERROR! alg_WSPSR()\n");

    destroy_permutation(&up);
}/*}}}*/



