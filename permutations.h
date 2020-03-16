/**********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 **********************************************************/
 
 
#ifndef __PERMUTATIONS_H
#define __PERMUTATIONS_H

#define SIGNED 1
#define UNSIGNED 0

struct permutation {
    int *pi;
    int *inv_pi;
    int size;
    int ini;
    double alpha;
    char sig;
};
typedef struct permutation permutation_t;

void create_permutation(permutation_t *p, int n, char sig);
void destroy_permutation(permutation_t *p);
int is_identity(permutation_t *p);
int is_reverse(permutation_t *p);
int is_signed_reverse(permutation_t *p);
void print_permutation(permutation_t *p);
void next_permutation(permutation_t *p);
int nof_inversions(permutation_t *p);
void prmcpy(permutation_t *dst, permutation_t *src);
void uns_prmcpy(permutation_t *dst, permutation_t *src);
void prm_image(permutation_t *p, permutation_t *up);
void prm_back_from_image(permutation_t *p, permutation_t *up);
void prm2str(permutation_t *prm, char str[]);
int prmcmp(permutation_t *p1, permutation_t *p2);
void make_identity(permutation_t *p);
void prm_insert(permutation_t *p, int elem, int pos);
void prmfill_fromstr(char strprm[], permutation_t *p);
void create_random_perm(permutation_t *p);
int prm_median(permutation_t *p, int inicio, int fim);
void prm_swap(permutation_t *p, int pos1, int pos2);
int prm_exists(permutation_t *p, int elem);
void find_extremes_of_strip(permutation_t *p, int elem, int *inicio, int *fim);
void find_extremes_of_signed_strip(permutation_t *p, int elem, int *inicio, int *fim);
void prm_separated(permutation_t *p, int *i, int *j, int *n);

unsigned long rank(permutation_t *p);
void unrank(unsigned long rank, permutation_t *p);

int is_sorted_interval(permutation_t *p, int ini, int fim, int type);

#endif /* __PERMUTATIONS_H */

