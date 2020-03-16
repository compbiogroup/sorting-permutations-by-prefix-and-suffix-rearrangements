/***********************************************************
 * Created: Qui 27 Ago 2015 15:57:33 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************
 *
 * Compare algorithm1 against algorithm2 for each permutation
 *
 ***********************************************************/
 
#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]) {
    FILE *f1, *f2;
    char perm[5000];
    int nmoves1, nmoves2, null;
    long double weight1, weight2;
    unsigned int nperms, alg1_greater_alg2, alg2_greater_alg1;
    unsigned int alg1_less_op_than_alg2, alg2_less_op_than_alg1;

    if (argc != 4) {
        printf("Usage: %s size file_results_alg1 file_results_alg2\n", argv[0]);
        return 0;
    }

    f1 = fopen64(argv[2], "r");
    f2 = fopen64(argv[3], "r");

    nperms = alg1_greater_alg2 = alg2_greater_alg1 = 0;
    alg1_less_op_than_alg2 = alg2_less_op_than_alg1 = 0;

    while (fscanf(f1, "%s %d %Lf %d %d %d %d %d %d", perm, &nmoves1,
                            &weight1, &null, &null, &null, &null, &null, &null) != EOF) {
        fscanf(f2, "%s %d %Lf %d %d %d %d %d %d", perm, &nmoves2,
                            &weight2, &null, &null, &null, &null, &null, &null);

        nperms++;
        if (nmoves1 > nmoves2) {
            alg1_greater_alg2++;
            alg2_less_op_than_alg1 += nmoves1 - nmoves2;
        } else if (nmoves2 > nmoves1) {
            alg2_greater_alg1++;
            alg1_less_op_than_alg2 += nmoves2 - nmoves1;
        }
    }

    printf("%d %d %Lf %d %Lf %d %Lf %Lf %Lf\n", atoi(argv[1]), alg1_greater_alg2,
                    (long double) alg1_greater_alg2/nperms*100,
                    alg2_greater_alg1, (long double) alg2_greater_alg1/nperms*100,
                    (nperms-alg1_greater_alg2-alg2_greater_alg1),
                    (long double) (nperms-alg1_greater_alg2-alg2_greater_alg1)/nperms*100,
                    (long double) alg1_less_op_than_alg2/nperms,
                    (long double) alg2_less_op_than_alg1/nperms);

    fclose(f1);
    fclose(f2);
    
    return 0;
}
