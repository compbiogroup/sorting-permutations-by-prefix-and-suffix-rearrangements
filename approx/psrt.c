/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "../permutations.h"
#include "../rearrangements.h"
#include "../util.h"
#include "../breakpoints.h"


void alg_2PSRT(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int i, j, elementsP, elementsS;
    int ip, jp, is, js;

    while (!is_identity(p)) {

        if (is_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
            continue;
        }

        /* try to remove two breakpoints with prefix transposition *//*{{{*/
        j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nmoves++;
            info->nof_pt++;
            continue;
        }

        j = p->inv_pi[dec_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nmoves++;
            info->nof_pt++;
            continue;
        }

        j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[dec_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) &&
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nmoves++;
            info->nof_pt++;
            continue;
        }

        j = p->inv_pi[inc_mod(p->pi[1], p->size)] + 1;
        i = p->inv_pi[inc_mod(p->pi[j], p->size)] + 1;
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) &&
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_prefix_transposition(p, i, j);
            print_verbose(p, verbose, 't', 1, i, j);
            info->nmoves++;
            info->nof_pt++;
            continue;
        }
        /*}}}*/

        /* try to remove two breakpoints with suffix transposition *//*{{{*/
        i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_suffix_transposition(p, i, j);
            print_verbose(p, verbose, 't', i, j, p->size+1);
            info->nmoves++;
            info->nof_st++;
            continue;
        }

        i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_suffix_transposition(p, i, j);
            print_verbose(p, verbose, 't', i, j, p->size+1);
            info->nmoves++;
            info->nof_st++;
            continue;
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[dec_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_suffix_transposition(p, i, j);
            print_verbose(p, verbose, 't', i, j, p->size+1);
            info->nmoves++;
            info->nof_st++;
            continue;
        }

        i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
        j = p->inv_pi[inc_mod(p->pi[i-1], p->size)];
        if (2 <= i && i < j && j <= p->size &&
                is_UPSRT_breakpoint(p, i-1, i) && 
                is_UPSRT_breakpoint(p, j-1, j)) {
            op_suffix_transposition(p, i, j);
            print_verbose(p, verbose, 't', i, j, p->size+1);
            info->nmoves++;
            info->nof_st++;
            continue;
        }
        /*}}}*/

        /* remove only one breakpoint */
        elementsP = elementsS = 2*p->size;
        ip = jp = is = js = 0;

        if (p->pi[1] == dec_mod(p->pi[2], p->size)) {/*{{{*/
            /* first element is in an increasing strip */
            j = p->inv_pi[dec_mod(p->pi[1], p->size)];
            if (!is_UPSRT_breakpoint(p, j-1, j)) {
                /* element in pi_j is in an increasing strip */
                i = j-1;
                while (p->pi[i] == dec_mod(p->pi[i+1], p->size)) i--;
                if (2 <= i && i < j && j <= p->size) {
                    elementsP = j;
                    ip = i+1;
                    jp = j+1;
                }
            } else {
                /* element in pi_j is a singleton or decreasing strip */
                elementsP = j-1;
                ip = j-1;
            }
        } else if (p->pi[1] == inc_mod(p->pi[2], p->size)) {
            /* first element is in a decreasing strip */
            j = p->inv_pi[inc_mod(p->pi[1], p->size)];
            if (!is_UPSRT_breakpoint(p, j-1, j)) {
                /* element in pi_j is in a decreasing strip */
                i = j-1;
                while (p->pi[i] == inc_mod(p->pi[i+1], p->size)) i--;
                if (2 <= i && i < j && j <= p->size) {
                    elementsP = j;
                    ip = i+1;
                    jp = j+1;
                }
            } else {
                /* element in pi_j is a singleton or increasing strip */
                elementsP = j-1;
                ip = j-1;
            }
        } else {
            /* first element is a singleton */
            j = p->inv_pi[dec_mod(p->pi[1], p->size)];
            if (!is_UPSRT_breakpoint(p, j-1, j)) {
                /* element in pi_j is in an increasing strip */
                if (j+1 > 2) {
                    elementsP = j;
                    ip = 2;
                    jp = j+1;
                }
            } else {
                /* element in pi_j is a singleton or begin of a strip */
                elementsP = j-1;
                ip = j-1;
            }
        }/*}}}*/

        if (p->pi[p->size] == inc_mod(p->pi[p->size-1], p->size)) {/*{{{*/
            /* last element is in an increasing strip */
            i = p->inv_pi[inc_mod(p->pi[p->size], p->size)];
            if (!is_UPSRT_breakpoint(p, i, i+1)) {
                /* element in pi_i is in an increasing strip */
                j = i+1;
                while (p->pi[j] == inc_mod(p->pi[j-1], p->size)) j++;
                if (2 <= i && i < j && j <= p->size) {
                    elementsS = p->size-i+1;
                    is = i;
                    js = j;
                }
            } else {
                /* element in pi_i is a singleton or decreasing strip */
                elementsS = p->size-i;
                is = i+1;
            }
        } else if (p->pi[p->size] == dec_mod(p->pi[p->size-1], p->size)) {
            /* last element is in a decreasing strip */
            i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
            if (!is_UPSRT_breakpoint(p, i, i+1)) {
                /* element in pi_i is in a decreasing strip */
                j = i+1;
                while (p->pi[j] == dec_mod(p->pi[j-1], p->size)) j++;
                if (2 <= i && i < j && j <= p->size) {
                    elementsS = p->size-i+1;
                    is = i;
                    js = j;
                }
            } else {
                /* element in pi_i is a singleton or increasing strip */
                elementsS = p->size-i;
                is = i+1;
            }
        } else {
            /* last element is singleton */
            i = p->inv_pi[dec_mod(p->pi[p->size], p->size)];
            if (!is_UPSRT_breakpoint(p, i, i+1)) {
                /* element in pi_i is in a decreasing strip */
                if (i < p->size) {
                    elementsS = p->size-i+1;
                    is = i;
                    js = p->size;
                }
            } else {
                /* element in pi_i is a singleton or begin of a strip */
                elementsS = p->size-i;
                is = i+1;
            }
        }/*}}}*/

        if (elementsP < elementsS || (elementsP == elementsS && elementsP != 2*p->size)) {/*{{{*/
            if (jp == 0) {
                op_prefix_reversal(p, ip);
                print_verbose(p, verbose, 'r', 1, ip, 0);
                info->nof_pr++;
            } else {
                op_prefix_transposition(p, ip, jp);
                print_verbose(p, verbose, 't', 1, ip, jp);
                info->nof_pt++;
            }
        } else if (elementsS < elementsP) {
            if (js == 0) {
                op_suffix_reversal(p, is);
                print_verbose(p, verbose, 'r', is, p->size, 0);
                info->nof_sr++;
            } else {
                op_suffix_transposition(p, is, js);
                print_verbose(p, verbose, 't', is, js, p->size+1);
                info->nof_st++;
            }
        } else {
            /* permutation has two strips */
            i = p->inv_pi[1];
            if (p->size == p->pi[i+1]) {
                op_prefix_transposition(p, i+1, p->size+1);
                print_verbose(p, verbose, 't', 1, i+1, p->size+1);
            } else {
                op_prefix_transposition(p, i, p->size+1);
                print_verbose(p, verbose, 't', 1, i, p->size+1);
            }
            info->nof_pt++;
        }/*}}}*/

        info->nmoves++;
    }
}/*}}}*/



void alg_2SPSRT(permutation_t *p, int verbose, info_t *info) {/*{{{*/
    int i, j, ip, jp, is, js;
    int elementsP, elementsS;

    while (!is_identity(p)) {

        if (is_signed_reverse(p)) {
            op_prefix_reversal(p, p->size);
            print_verbose(p, verbose, 'r', 1, p->size, 0);
            info->nmoves++;
            info->nof_pr++;
            continue;
        }

        /* tries to remove 2 breakpoints with prefix */
        j = prm_exists(p, dec_mod(p->pi[1], p->size));
        if (j > 0) {
            j++;
            i = prm_exists(p, dec_mod(p->pi[j], p->size));
            if (2 <= i+1 && i+1 < j && j <= p->size &&
                            is_PSRT_breakpoint(p, i, i+1) &&
                            is_PSRT_breakpoint(p, j-1, j)) {
                op_prefix_transposition(p, i+1, j);
                print_verbose(p, verbose, 't', 1, i+1, j);
                info->nmoves++;
                info->nof_pt++;
                continue;
            }
        }

        /* tries to remove 2 breakpoints with suffix */
        i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
        if (i > 0) {
            j = prm_exists(p, inc_mod(p->pi[i-1], p->size));
            if (j > 0 && i < j && j <= p->size &&
                            is_PSRT_breakpoint(p, i-1, i) &&
                            is_PSRT_breakpoint(p, j-1, j)) {
                op_suffix_transposition(p, i, j);
                print_verbose(p, verbose, 't', i, j, p->size+1);
                info->nmoves++;
                info->nof_st++;
                continue;
            }
        }

        elementsP = elementsS = 2*p->size;
        ip = jp = is = js = 0;

        j = prm_exists(p, dec_mod(p->pi[1], p->size));
        if (j > 0) {
            i = 2;
            while (p->pi[i] == inc_mod(p->pi[i-1], p->size)) i++;
            if (2 <= i && i < j+1 && j+1 <= p->size) {
                ip = i;
                jp = j+1;
                elementsP = j;
            }
        }

        i = prm_exists(p, -dec_mod(p->pi[1], p->size));
        if (i > 0 && is_PSRT_breakpoint(p, i-1, i)) {
            ip = i-1;
            elementsP = i-1;
        }

        i = prm_exists(p, inc_mod(p->pi[p->size], p->size));
        if (i > 0) {
            j = p->size-1;
            while (p->pi[j] == dec_mod(p->pi[j+1], p->size)) j--;
            if (i <= 2 && i < j+1 && j+1 <= p->size) {
                is = i;
                js = j+1;
                elementsS = p->size-i+1;
            }
        }

        i = prm_exists(p, -inc_mod(p->pi[p->size], p->size));
        if (i > 0 && is_PSRT_breakpoint(p, i, i+1)) {
            is = i+1;
            elementsS = p->size-i;
        }


        if (elementsP < elementsS || (elementsP == elementsS && elementsP != 2*p->size)) {
            if (jp == 0) {
                op_prefix_reversal(p, ip);
                print_verbose(p, verbose, 'r', 1, ip, 0);
                info->nof_pr++;
            } else {
                op_prefix_transposition(p, ip, jp);
                print_verbose(p, verbose, 't', 1, ip, jp);
                info->nof_pt++;
            }
        } else if (elementsS < elementsP) {
            if (js == 0) {
                op_suffix_reversal(p, is);
                print_verbose(p, verbose, 'r', is, p->size, 0);
                info->nof_sr++;
            } else {
                op_suffix_transposition(p, is, js);
                print_verbose(p, verbose, 't', is, js, p->size+1);
                info->nof_st++;
            }
        } else {
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
        }

        info->nmoves++;
    }
}/*}}}*/

