/***********************************************************
 * Created: Sáb 14 Dez 2013 17:51:12 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/
 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "permutations.h"
#include "rearrangements.h"
#include "algorithms.h"
#include "util.h"

#include "approx/pr.h"
#include "approx/psr.h"
#include "approx/pt.h"
#include "approx/pst.h"
#include "approx/prt.h"
#include "approx/psrt.h"
#include "approx/improved.h"

#include "approx/wpr.h"
#include "approx/wsr.h"
#include "approx/wpsr.h"
#include "approx/wpt.h"
#include "approx/wst.h"
#include "approx/wpst.h"
#include "approx/wprt.h"
#include "approx/wsrt.h"
#include "approx/wpsrt.h"

void deal_alg(char alg_name[], int *alg) {/*{{{*/
    if (!strcmp(alg_name, "2-PR"))
        (*alg) = PR2;
    else if (!strcmp(alg_name, "2-PRx"))
        (*alg) = PRx2;
    
    else if (!strcmp(alg_name, "2-PSR"))
        (*alg) = PSR2;
    else if (!strcmp(alg_name, "2-PSRx"))
        (*alg) = PSRx2;

    else if (!strcmp(alg_name, "2-PT"))
        (*alg) = PT2;
    else if (!strcmp(alg_name, "2-PTx"))
        (*alg) = PTx2;

    else if (!strcmp(alg_name, "2-PST"))
        (*alg) = PST2;
    else if (!strcmp(alg_name, "2-PSTx"))
        (*alg) = PSTx2;

    else if (!strcmp(alg_name, "2-PRT"))
        (*alg) = PRT2;
    else if (!strcmp(alg_name, "2-PRTx"))
        (*alg) = PRTx2;

    else if (!strcmp(alg_name, "2-PSRT"))
        (*alg) = PSRT2;
    else if (!strcmp(alg_name, "2-PSRTx"))
        (*alg) = PSRTx2;

    else if (!strcmp(alg_name, "2-SPR"))
        (*alg) = SPR2;
    else if (!strcmp(alg_name, "2-SPRx"))
        (*alg) = SPRx2;
    
    else if (!strcmp(alg_name, "2-SPSR"))
        (*alg) = SPSR2;
    else if (!strcmp(alg_name, "2-SPSRx"))
        (*alg) = SPSRx2;

    else if (!strcmp(alg_name, "2-SPRT"))
        (*alg) = SPRT2;
    else if (!strcmp(alg_name, "2-SPRTx"))
        (*alg) = SPRTx2;

    else if (!strcmp(alg_name, "2-SPSRT"))
        (*alg) = SPSRT2;
    else if (!strcmp(alg_name, "2-SPSRTx"))
        (*alg) = SPSRTx2;

    else if (!strcmp(alg_name, "WPRm"))
        (*alg) = cWPRm;
    else if (!strcmp(alg_name, "WPRg"))
        (*alg) = cWPRg;
    else if (!strcmp(alg_name, "WPR"))
        (*alg) = cWPR;
    else if (!strcmp(alg_name, "WSR"))
        (*alg) = cWSR;

    else if (!strcmp(alg_name, "WPSRg"))
        (*alg) = cWPSRg;
    else if (!strcmp(alg_name, "WPSR"))
        (*alg) = cWPSR;

    else if (!strcmp(alg_name, "WPTg"))
        (*alg) = cWPTg;
    else if (!strcmp(alg_name, "WPT"))
        (*alg) = cWPT;
    else if (!strcmp(alg_name, "WST"))
        (*alg) = cWST;

    else if (!strcmp(alg_name, "WPSTg"))
        (*alg) = cWPSTg;
    else if (!strcmp(alg_name, "WPST"))
        (*alg) = cWPST;

    else if (!strcmp(alg_name, "WPRTg"))
        (*alg) = cWPRTg;
    else if (!strcmp(alg_name, "WPRT"))
        (*alg) = cWPRT;
    else if (!strcmp(alg_name, "WSRT"))
        (*alg) = cWSRT;

    else if (!strcmp(alg_name, "WPSRTg"))
        (*alg) = cWPSRTg;
    else if (!strcmp(alg_name, "WPSRT"))
        (*alg) = cWPSRT;

    else if (!strcmp(alg_name, "WSPRg"))
        (*alg) = cWSPRg;
    else if (!strcmp(alg_name, "WSPR"))
        (*alg) = cWSPR;
    else if (!strcmp(alg_name, "WSSR"))
        (*alg) = cWSSR;

    else if (!strcmp(alg_name, "WSPSRg"))
        (*alg) = cWSPSRg;
    else if (!strcmp(alg_name, "WSPSR"))
        (*alg) = cWSPSR;

    else if (!strcmp(alg_name, "WSPRTg"))
        (*alg) = cWSPRTg;
    else if (!strcmp(alg_name, "WSPRT"))
        (*alg) = cWSPRT;
    else if (!strcmp(alg_name, "WSSRT"))
        (*alg) = cWSSRT;

    else if (!strcmp(alg_name, "WSPSRTg"))
        (*alg) = cWSPSRTg;
    else if (!strcmp(alg_name, "WSPSRT"))
        (*alg) = cWSPSRT;

    else {
        printf("ERROR! deal_alg(%s): invalid algorithm.\n", alg_name);
        exit(0);
    }
}/*}}}*/

void run(int alg, permutation_t *p, int verbose, info_t *info) {/*{{{*/
    switch (alg) {
        case PR2:
            alg_2PR(p, verbose, info);
            break;
        case PRx2:
            alg_2PRx(p, verbose, info);
            break;
        case PSR2:
            alg_2PSR(p, verbose, info);
            break;
        case PSRx2:
            alg_2PSRx(p, verbose, info);
            break;
        case PT2:
            alg_2PT(p, verbose, info);
            break;
        case PTx2:
            alg_2PTx(p, verbose, info);
            break;
        case PST2:
            alg_2PST(p, verbose, info);
            break;
        case PSTx2:
            alg_2PSTx(p, verbose, info);
            break;
        case PRT2:
            alg_2PRT(p, verbose, info);
            break;
        case PRTx2:
            alg_2PRTx(p, verbose, info);
            break;
        case PSRT2:
            alg_2PSRT(p, verbose, info);
            break;
        case PSRTx2:
            alg_2PSRTx(p, verbose, info);
            break;
        case SPR2:
            alg_2SPR(p, verbose, info);
            break;
        case SPRx2:
            alg_2SPRx(p, verbose, info);
            break;
        case SPSR2:
            alg_2SPSR(p, verbose, info);
            break;
        case SPSRx2:
            alg_2SPSRx(p, verbose, info);
            break;
        case SPRT2:
            alg_2SPRT(p, verbose, info);
            break;
        case SPRTx2:
            alg_2SPRTx(p, verbose, info);
            break;
        case SPSRT2:
            alg_2SPSRT(p, verbose, info);
            break;
        case SPSRTx2:
            alg_2SPSRTx(p, verbose, info);
            break;
        case cWPRm:
            alg_WPRm(p, verbose, info);
            break;
        case cWPRg:
            alg_WPRg(p, verbose, info);
            break;
        case cWPR:
            alg_WPR(p, verbose, info);
            break;
        case cWSR:
            alg_WSR(p, verbose, info);
            break;
        case cWPSRg:
            alg_WPSRg(p, verbose, info);
            break;
        case cWPSR:
            alg_WPSR(p, verbose, info);
            break;
        case cWPTg:
            alg_WPT(p, verbose, info);
            break;
        case cWPT:
            alg_WPT(p, verbose, info);
            break;
        case cWST:
            alg_WST(p, verbose, info);
            break;
        case cWPSTg:
            alg_WPSTg(p, verbose, info);
            break;
        case cWPST:
            alg_WPST(p, verbose, info);
            break;
        case cWPRTg:
            alg_WPRTg(p, verbose, info);
            break;
        case cWPRT:
            alg_WPRT(p, verbose, info);
            break;
        case cWSRT:
            alg_WSRT(p, verbose, info);
            break;
        case cWPSRTg:
            alg_WPSRTg(p, verbose, info);
            break;
        case cWPSRT:
            alg_WPSRT(p, verbose, info);
            break;
        case cWSPRg:
            alg_WSPRg(p, verbose, info);
            break;
        case cWSPR:
            alg_WSPR(p, verbose, info);
            break;
        case cWSSR:
            alg_WSSR(p, verbose, info);
            break;
        case cWSPSRg:
            alg_WSPSRg(p, verbose, info);
            break;
        case cWSPSR:
            alg_WSPSR(p, verbose, info);
            break;
        case cWSPRTg:
            alg_WSPRTg(p, verbose, info);
            break;
        case cWSPRT:
            alg_WSPRT(p, verbose, info);
            break;
        case cWSSRT:
            alg_WSSRT(p, verbose, info);
            break;
        case cWSPSRTg:
            alg_WSPSRTg(p, verbose, info);
            break;
        case cWSPSRT:
            alg_WSPSRT(p, verbose, info);
            break;

        default:
            printf("ERROR! run(): invalid algorithm\n");
            exit(0);
    }
}/*}}}*/

void print_algs(void) {/*{{{*/
    printf("\t\t2-PR: 2-approximation algorithm for SbPR\n");
    printf("\t\t2-PRx: improved 2-approximation algorithm for SbPR\n");
    printf("\t\t2-PSR: 2-approximation algorithm for SbPSR\n");
    printf("\t\t2-PSRx: improved 2-approximation algorithm for SbPSR\n");
    printf("\t\t2-PT: 2-approximation algorithm for SbPT\n");
    printf("\t\t2-PTx: improved 2-approximation algorithm for SbPT\n");
    printf("\t\t2-PST: 2-approximation algorithm for SbPST\n");
    printf("\t\t2-PSTx: improved 2-approximation algorithm for SbPST\n");
    printf("\t\t2-PRT: 2-approximation algorithm for SbPRT\n");
    printf("\t\t2-PRTx: improved 2-approximation algorithm for SbPRT\n");
    printf("\t\t2-PSRT: 2-approximation algorithm for SbPSRT\n");
    printf("\t\t2-PSRTx: improved 2-approximation algorithm for SbPSRT\n");

    printf("\t\t2-SPR: 2-approximation algorithm for SbSPR\n");
    printf("\t\t2-SPRx: improved 2-approximation algorithm for SbSPR\n");
    printf("\t\t2-SPSR: 2-approximation algorithm for SbSPSR\n");
    printf("\t\t2-SPSRx: improved 2-approximation algorithm for SbSPSR\n");
    printf("\t\t2-SPRT: 2-approximation algorithm for SbSPRT\n");
    printf("\t\t2-SPRTx: improved 2-approximation algorithm for SbSPRT\n");
    printf("\t\t2-SPSRT: 2-approximation algorithm for SbSPSRT\n");
    printf("\t\t2-SPSRTx: improved 2-approximation algorithm for SbSPSRT\n");

    printf("\t\tWPRm: algorithm for SbWPR\n");
    printf("\t\tWPTg: algorithm for SbWPR\n");
    printf("\t\tWPR: algorithm for SbWPR\n");
    printf("\t\tWSR: algorithm for SbWSR\n");
    printf("\t\tWPSRg: algorithm for SbWPSR\n");
    printf("\t\tWPSR: algorithm for SbWPSR\n");
    printf("\t\tWPTg: algorithm for SbWPT\n");
    printf("\t\tWPT: algorithm for SbWPT\n");
    printf("\t\tWST: algorithm for SbWST\n");
    printf("\t\tWPSTg: algorithm for SbWPST\n");
    printf("\t\tWPST: algorithm for SbWPST\n");
    printf("\t\tWPRTg: algorithm for SbWPRT\n");
    printf("\t\tWPRT: algorithm for SbWPRT\n");
    printf("\t\tWSRT: algorithm for SbWSRT\n");
    printf("\t\tWPSRTg: algorithm for SbWPSRT\n");
    printf("\t\tWPSRT: algorithm for SbWPSRT\n");

    printf("\t\tWSPTg: algorithm for SbWSPR\n");
    printf("\t\tWSPR: algorithm for SbWSPR\n");
    printf("\t\tWSSR: algorithm for SbWSSR\n");
    printf("\t\tWSPSRg: algorithm for SbWSPSR\n");
    printf("\t\tWSPSR: algorithm for SbWSPSR\n");
    printf("\t\tWSPRTg: algorithm for SbWSPRT\n");
    printf("\t\tWSPRT: algorithm for SbWSPRT\n");
    printf("\t\tWSSRT: algorithm for SbWSSRT\n");
    printf("\t\tWSPSRTg: algorithm for SbWSPSRT\n");
    printf("\t\tWSPSRT: algorithm for SbWSPSRT\n");
}/*}}}*/

void alg_approx(char prob[], permutation_t *p, int verbose, info_t *info) {/*{{{*/
    if (!strcmp(prob, "pr"))
        alg_2PR(p, verbose, info);
    
    else if (!strcmp(prob, "psr"))
        alg_2PSR(p, verbose, info);

    else if (!strcmp(prob, "pt"))
        alg_2PT(p, verbose, info);

    else if (!strcmp(prob, "pst"))
        alg_2PST(p, verbose, info);

    else if (!strcmp(prob, "prt"))
        alg_2PRT(p, verbose, info);

    else if (!strcmp(prob, "psrt"))
        alg_2PSRT(p, verbose, info);

    else if (!strcmp(prob, "spr"))
        alg_2SPR(p, verbose, info);
    
    else if (!strcmp(prob, "spsr"))
        alg_2SPSR(p, verbose, info);

    else if (!strcmp(prob, "sprt"))
        alg_2SPRT(p, verbose, info);

    else if (!strcmp(prob, "spsrt"))
        alg_2SPSRT(p, verbose, info);

    else {
        printf("ERROR! alg_approx(): invalid algorithm.\n");
        exit(0);
    }
}/*}}}*/

