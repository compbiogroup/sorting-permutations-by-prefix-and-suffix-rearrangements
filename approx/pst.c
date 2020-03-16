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



void alg_2PST(permutation_t *p, int verbose, info_t* info) {/*{{{*/
    int ip, jp, is, js;

    while (!is_identity(p)) {
        jp = p->inv_pi[p->pi[1] - 1] + 1;
        ip = p->inv_pi[p->pi[jp] - 1] + 1;

        is = p->inv_pi[p->pi[p->size] + 1];
        js = p->inv_pi[p->pi[is-1] + 1];

        if (2 <= ip && ip < jp && jp <= p->size) {
            if (2 <= is && is < js && js <= p->size) {
                if (jp-1 <= p->size-is+1) {
                    op_prefix_transposition(p, ip, jp);
                    print_verbose(p, verbose, 't', 1, ip, jp); 
                    info->nof_pt++;
                } else {
                    op_suffix_transposition(p, is, js);
                    print_verbose(p, verbose, 't', is, js, p->size+1);
                    info->nof_st++;
                }
            } else {
                op_prefix_transposition(p, ip, jp);
                print_verbose(p, verbose, 't', 1, ip, jp);
                info->nof_pt++;
            }
        } else if (2 <= is && is < js && js <= p->size) {
            op_suffix_transposition(p, is, js);
            print_verbose(p, verbose, 't', is, js, p->size+1);
            info->nof_st++;
        } else {
            /* remove only one */
            ip = 1;
            while (p->pi[ip] == p->pi[ip+1] - 1)
                ip++;

            if (p->pi[ip] == p->size)
                jp = p->inv_pi[p->pi[1] - 1] + 1;
            else
                jp = p->inv_pi[p->pi[ip] + 1];

            js = p->size;
            while (p->pi[js] == p->pi[js-1]+1)
                js--;

            if (p->pi[js] == 1)
                is = p->inv_pi[p->pi[p->size] + 1] - 1;
            else
                is = p->inv_pi[p->pi[js] - 1];

            if (jp-1 <= p->size-is) {
                op_prefix_transposition(p, ip+1, jp);
                print_verbose(p, verbose, 't', 1, ip+1, jp);
                info->nof_pt++;
            } else {
                op_suffix_transposition(p, is+1, js);
                print_verbose(p, verbose, 't', is+1, js, p->size+1);
                info->nof_st++;
            }
        }
        info->nmoves++;
    }
}/*}}}*/

