/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/


#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <float.h>
#include <math.h>

#include "util.h"
#include "algorithms.h"
#include "permutations.h"


void test_inputs(int size, int nperms, FILE *fin, char *alg_name);
void print_help(char progname[]);

int main(int argc, char *argv[]) {

    extern char *optarg;
    char op, strout[500], strprm[6500], *alg_name = NULL, signal = 0;
    int alg, verbose = 1, size = 0, nperms = 0, n;
    double avgnmoves, alpha = 0;
    long double avgweight;
    double avgnof_pr, avgnof_sr, avgnof_pt, avgnof_st, avgnof_t, avgnof_r;
    FILE *fin = NULL, *fout;
    permutation_t p;
    info_t info;

    struct option longopts[] = {
        {"algorithm", 1, NULL, 'a'},
        {"size", 1, NULL, 'n'},
        {"signal", 1, NULL, 's'},
        {"nperms", 1, NULL, 'q'},
        {"input", 1, NULL, 'i'},
        {"alpha", 1, NULL, 'x'},
        {"verbose", 1, NULL, 'v'},
        {"help", 0, NULL, 'h'},
    };

    while ((op = getopt_long(argc, argv, "a:n:s:q:i:x:v:h", longopts, NULL)) != -1) {
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
            case 'q':
                nperms = atoi(optarg);
                break;
            case 'i':
                fin = fopen(optarg, "r");
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

    test_inputs(size, nperms, fin, alg_name);

    create_permutation(&p, size, signal);
    p.alpha = alpha;

    if (alpha != 0)
        sprintf(strout, "%s_%.3f.%d", alg_name, alpha, size);
    else
        sprintf(strout, "%s.%d", alg_name, size);
    fout = fopen(strout, "w");

    avgnmoves = avgweight = 0;
    avgnof_pr = avgnof_sr = avgnof_r = 0;
    avgnof_pt = avgnof_st = avgnof_t = 0;

    for (n = 1; n <= nperms; n++) {
        fscanf(fin, "%s", strprm);
        fprintf(fout, "%s ", strprm);
        prmfill_fromstr(strprm, &p);
        init_info(&info);
        info.alpha = alpha;

        if (verbose == 3 || verbose == 4) {
            printf("\n");
            print_permutation(&p);
        }
        
        run(alg, &p, verbose, &info);

        fprintf(fout, "%d %Lf %d %d %d %d %d %d\n", info.nmoves, info.weight,
                        info.nof_r, info.nof_pr, info.nof_sr, info.nof_t,
                        info.nof_pt, info.nof_st);

        avgnmoves += info.nmoves;
        avgweight += info.weight/nperms;
        avgnof_pr += info.nof_pr;
        avgnof_sr += info.nof_sr;
        avgnof_pt += info.nof_pt;
        avgnof_st += info.nof_st;
        avgnof_r += info.nof_r;
        avgnof_t += info.nof_t;

        if (verbose == 2 || verbose == 4) 
            printf("%d %Lf\n", info.nmoves, info.weight);
    }

    avgnof_pr = avgnof_pr / nperms;
    avgnof_sr = avgnof_sr / nperms;
    avgnof_pt = avgnof_pt / nperms;
    avgnof_st = avgnof_st / nperms;
    avgnof_r = avgnof_r / nperms;
    avgnof_t = avgnof_t / nperms;

    avgnmoves = avgnmoves/nperms;

    printf("%d %f %Lf %.3f %.3f %.3f %.3f %.3f %.3f\n", size, avgnmoves,
                    avgweight, avgnof_r, avgnof_pr, avgnof_sr, avgnof_t,
                    avgnof_pt, avgnof_st);

    fclose(fin);
    fclose(fout);
    destroy_permutation(&p);

    return 0;
}

void test_inputs(int size, int nperms, FILE *fin, char *alg_name) {/*{{{*/
    if (size == 0) {
        printf("A size of permutation must be specified.\n");
        exit(0);
    } 
    
    if (nperms == 0) {
        printf("An amount of permutations must be specified.\n");
        exit(0);
    } 

    if (!fin) {
        printf("A valid input file must be provided.\n");
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
    printf("\t--nperms <n> or -q <n>: number of permutations to be tested\n");
    printf("\t--signal <s> or -s <s>: signal (0=unsigned, 1=signed)\n");
    printf("\t--input <file> or -i <file>: gets the input from <file>\n");
    printf("\t It must contain <nperms> of size <size>, in which the elements are separated by comma\n");
    printf("\t--alpha <value> or -x <value>: alpha for length-weighted problems (default = 0)\n");
    printf("\t--verbose <n> or -v <n>: define the level of information to be printed\n");
    printf("\t\tuse 2 for distance\n");
    printf("\t\tuse 3 for sorting sequence\n");
    printf("\t\tuse 4 for both\n");
    printf("\t\tuse 5 for rearrangement events sequence\n");
    printf("\n");
}/*}}}*/


