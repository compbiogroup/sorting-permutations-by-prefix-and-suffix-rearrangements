/***********************************************************
 * Created: Sex 28 Ago 2015 09:56:43 BRT
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************
 *
 * Check if approximation factors of algorithm are greater than some value 
 *
 ***********************************************************/
 
#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]) {
    FILE *f;
    int n;
    long double approx;
    double factor;
    unsigned int greater, nperms;

    if (argc != 4) {
        printf("Usage: %s size file_results_approx factor\n", argv[0]);
    }

    n = atoi(argv[1]);
    f = fopen64(argv[2], "r");
    factor = atof(argv[3]);

    nperms = greater = 0;
    while (fscanf(f, "%Lf", &approx) != EOF) {
        if (approx > factor)
            greater++;
        nperms++;
    }

    if (greater > 0)
        printf("%d %d %d %LF\n", n, greater, nperms, (long double) greater/nperms*100);
    
    fclose(f);

    return 0;
}
