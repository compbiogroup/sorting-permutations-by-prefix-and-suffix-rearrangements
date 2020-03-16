/***********************************************************
 * Created: Qui 17 Jan 2013 15:39:35 BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <math.h>

#include "util.h"
#include "algorithms.h"
#include "permutations.h"


void print_help(char progname[]);
void test_inputs(int size, char *strprm, char *alg_name);

int main(int argc, char *argv[]) {

    extern char *optarg;
    char op, *alg_name = NULL, *strprm = NULL, signal = 0;
    int alg, verbose = 1, size = 0;
    double alpha = 0;
    permutation_t p;
    info_t info;

    struct option longopts[] = {
        {"algorithm", 1, NULL, 'a'},
        {"size", 1, NULL, 'n'},
        {"signal", 1, NULL, 's'},
        {"perm", 1, NULL, 'x'},
        {"alpha", 1, NULL, 'p'},
        {"verbose", 1, NULL, 'v'},
        {"help", 0, NULL, 'h'},
    };

    while ((op = getopt_long(argc, argv, "a:n:s:p:x:v:h", longopts, NULL)) != -1) {
        switch (op) {
            case 'a':
                alg_name = optarg;
                deal_alg(alg_name, &alg);
                break;
            case 'n':
                size = atoi(optarg);
                break;
            case 's':
                signal = atoi(optarg);
                break;
            case 'p':
                strprm = optarg;
                break;
            case 'x':
                alpha = atof(optarg);
                break;
            case 'v':
                verbose = atoi(optarg);
                break;
            case 'h':
                print_help(argv[0]);
                exit(0);
        }
    }

    test_inputs(size, strprm, alg_name);

    create_permutation(&p, size, signal);
    p.alpha = alpha;

    prmfill_fromstr(strprm, &p);
    init_info(&info);
    info.alpha = alpha;

    if (verbose == 3 || verbose == 4)
        print_permutation(&p);

    run(alg, &p, verbose, &info);
    printf("%d %Lf %d %d %d %d %d %d\n", info.nmoves, info.weight, info.nof_r,
                    info.nof_pr, info.nof_sr, info.nof_t, info.nof_pt,
                    info.nof_st);

    if (verbose == 2 || verbose == 4) 
        printf("%d %Lf\n", info.nmoves, info.weight);

    destroy_permutation(&p);

    return 0;
}

void test_inputs(int size, char *strprm, char *alg_name) {/*{{{*/
    if (size == 0) {
        printf("A size of permutation must be specified.\n");
        exit(0);
    } 
    
    if (!strprm) {
        printf("A permutation must be provided.\n");
        exit(0);
    }

    if (!alg_name) {
        printf("A name of algorithm must be provided.\n");
        exit(0);
    }
}/*}}}*/

void print_help(char progname[]) {/*{{{*/
    printf("usage: %s [options] \n", progname);
    printf("options: \n");
    printf("\t--algorithm <name> or -a <name>: name of the algorithm to be tested. They can be:\n");
    print_algs();
    printf("\t--size <size> or -n <size>: size of the permutations to be tested\n");
    printf("\t--perm <p> or -p <p>: permutation to be tested. Example: 1,2,3,4,5,6\n");
    printf("\t--alpha <value> or -x <value>: alpha for length-weighted problems\n");
    printf("\t--verbose <n> or -v <n>: define the level of information to be printed\n");
    printf("\t\tuse 2 for distance\n");
    printf("\t\tuse 3 for sorting sequence\n");
    printf("\t\tuse 4 for both\n");
    printf("\t\tuse 5 for rearrangement events sequence\n");
    printf("\n");
}/*}}}*/

