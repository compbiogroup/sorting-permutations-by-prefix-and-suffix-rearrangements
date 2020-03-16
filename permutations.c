/***********************************************************
 * Created: Mon 04 Feb 2013 03:35:22 PM BRST
 *
 * Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "permutations.h"
#include "rearrangements.h"
#include "breakpoints.h"
#include "util.h"



static unsigned long fat(int size) {
    if (size == 2 || size == 1)
        return size;
    return size * fat(size - 1);
}

static unsigned long expo(int x, int y) {
    int i;
    unsigned long ret = x;

    for (i = 1; i < y; i++) {
        ret = ret * x;
    }
    return ret;
}
                
static unsigned long total(int size) {
    return expo(2, size) * fat(size);
}

void create_permutation(permutation_t *p, int n, char sig) {/*{{{*/
    p->ini = 1;
    p->size = n;
    p->pi = Malloc(sizeof(int) * (n+2));
    p->inv_pi = Malloc(sizeof(int) * (n+2));
    p->pi[0] = 0;
    p->pi[n+1] = n+1;
    p->inv_pi[0] = 0;
    p->inv_pi[n+1] = n+1;
    p->sig = sig;
}/*}}}*/

void destroy_permutation(permutation_t *p) {/*{{{*/
    free(p->pi);
    free(p->inv_pi);
}/*}}}*/

int is_identity(permutation_t *p) {/*{{{*/
    int i;
    for (i = 1; i <= p->size; i++)
        if (p->pi[i] != i)
            return 0;
    return 1;
}/*}}}*/

int is_reverse(permutation_t *p) {/*{{{*/
    int i;
    for (i = 1; i <= p->size; i++)
        if (p->pi[i] != p->size-i+1)
            return 0;
    return 1;
}/*}}}*/

int is_signed_reverse(permutation_t *p) {/*{{{*/
    int i;
    for (i = 1; i <= p->size; i++)
        if (p->pi[i] != -(p->size-i+1))
            return 0;
    return 1;
}/*}}}*/

void print_permutation(permutation_t *p) {/*{{{*/
    int i;
    printf("%d", p->pi[1]);
    for (i = 2; i <= p->size; i++)
        printf(",%d", p->pi[i]);
    printf("\n");
}/*}}}*/

void next_permutation(permutation_t *p) {/*{{{*/
    int i = p->size;
    int j, prox_i, aux;
    unsigned long r;

    if (p->sig == UNSIGNED) {
        /* Find i such that [i, n] is a decreasing strip */
        while (p->pi[i] < p->pi[i-1])
            i--;

        /* pi[prox_i] must be the smallest element greater than pi[i-1],
         * prox_i \in {i, ..., n} */
        prox_i = i;
        for (j = i+1; j <= p->size; j++) {
            if (p->pi[j] > p->pi[i-1] && 
                    abs(p->pi[j] - p->pi[i-1]) < abs(p->pi[prox_i] - p->pi[i-1]))
                prox_i = j;
        }
        
        if (i != 1) {
            aux = p->pi[i-1];
            p->pi[i-1] = p->pi[prox_i];
            p->pi[prox_i] = aux;
            p->inv_pi[p->pi[i-1]] = i-1;
            p->inv_pi[p->pi[prox_i]] = prox_i;
        }

        /* Invert [i, n] */
        for (j = i; j <= (p->size-i)/2 + i; j++) {
            aux = p->pi[j];
            p->pi[j] = p->pi[p->size - j + i];
            p->pi[p->size - j + i] = aux;
            p->inv_pi[p->pi[j]] = j;
            p->inv_pi[p->pi[p->size - j + i]] = p->size - j + i;
        }
    } else {
        r = rank(p);
        r = (r + 1) % total(p->size) + 1;
        unrank(r, p);
    }
}/*}}}*/

int nof_inversions(permutation_t *p) {/*{{{*/
    int i, j, number = 0;
    for (i = p->ini; i < p->size; i++)
        for (j = i+1; j < p->size+1; j++)
            if (abs(p->pi[i]) > abs(p->pi[j]))
                number++;
    return number;
}/*}}}*/

void prmcpy(permutation_t *dst, permutation_t *src) {/*{{{*/
    int i;
    for (i = 1; i <= src->size; i++) {
        dst->pi[i] = src->pi[i];
        dst->inv_pi[i] = src->inv_pi[i];
    }

    dst->pi[src->size + 1] = src->size + 1;
    dst->pi[0] = 0;
    dst->inv_pi[src->size + 1] = src->size + 1;
    dst->inv_pi[0] = 0;

    dst->ini = src->ini;
    dst->size = src->size;
    dst->sig = src->sig;
}/*}}}*/

void uns_prmcpy(permutation_t *dst, permutation_t *src) {/*{{{*/
    int i;
    for (i = 1; i <= src->size; i++) {
        dst->pi[i] = abs(src->pi[i]);
        dst->inv_pi[dst->pi[i]] = i;
    }

    dst->pi[src->size + 1] = src->size + 1;
    dst->pi[0] = 0;
    dst->inv_pi[src->size + 1] = src->size + 1;
    dst->inv_pi[0] = 0;

    dst->size = src->size;
    dst->sig = UNSIGNED;
}/*}}}*/

void prm_image(permutation_t *p, permutation_t *up) {/*{{{*/
    int i;
   
    up->pi[0] = 0;
    up->pi[2*p->size+1] = 2*p->size+1;

    for (i = 1; i <= p->size; i++) {
        if (p->pi[i] > 0) {
            up->pi[2*i-1] = 2 * p->pi[i] - 1;
            up->inv_pi[2 * p->pi[i] - 1] = 2*i-1;
            up->pi[2*i] = 2 * p->pi[i];
            up->inv_pi[2 * p->pi[i]] = 2*i;
        } else {
            up->pi[2*i-1] = -2 * p->pi[i];
            up->inv_pi[-2 * p->pi[i]] = 2*i-1;
            up->pi[2*i] = -2 * p->pi[i] - 1;
            up->inv_pi[-2 * p->pi[i] - 1] = 2*i;
        }
    }

    up->sig = UNSIGNED;
}/*}}}*/

void prm_back_from_image(permutation_t *p, permutation_t *up) {/*{{{*/
    int i;
   
    p->pi[0] = 0;
    p->pi[up->size/2 + 1] = up->size/2 + 1;

    for (i = 1; i <= up->size-1; i+=2) {
        if (up->pi[i] > up->pi[i+1]) {
            prm_insert(p, -(up->pi[i]/2), (i+1)/2);
        } else {
            prm_insert(p, up->pi[i+1]/2, (i+1)/2);
        }
    }
}/*}}}*/

void prm2str(permutation_t *prm, char str[]) {/*{{{*/
    int i;
    char num[5];
    sprintf(str, "%d", prm->pi[1]);
    for (i = 2; i <= prm->size; i++) {
        sprintf(num, ",%d", prm->pi[i]);
        strcat(str, num);
    }
}/*}}}*/

int prmcmp(permutation_t *p1, permutation_t *p2) {/*{{{*/
    int i;

    for (i = 0; i < p1->size; i++)
        if (p1->pi[i] != p2->pi[i])
            return i;

    return 0;
}/*}}}*/

void make_identity(permutation_t *p) {/*{{{*/
    int i;
    for (i = 0; i <= p->size+1; i++) {
        p->pi[i] = i;
        p->inv_pi[i] = i;
    }
}/*}}}*/

void prm_insert(permutation_t *p, int elem, int pos) {/*{{{*/
    p->pi[pos] = elem;
    if (elem < 0)
        p->inv_pi[-1*elem] = -1*pos;
    else
        p->inv_pi[elem] = pos;
}/*}}}*/

void prmfill_fromstr(char strprm[], permutation_t *p) {/*{{{*/
    char *r = NULL;
    int i;
    int size = 0;

    r = strtok(strprm, ",");
    i = 1;
    while (r != NULL) {
        prm_insert(p, atoi(r), i++);
        r = strtok(NULL, ",");
        if (abs(p->pi[i-1]) > size)
            size = abs(p->pi[i-1]);
        if (p->pi[i-1] < 0)
            p->sig = SIGNED;
    }
    prm_insert(p, 0, 0);
    prm_insert(p, size+1, size+1);
    p->size = size;
}/*}}}*/

void create_random_perm(permutation_t *p) {/*{{{*/
    int i, t;

    for (i = 1; i <= p->size; i++)
        p->pi[i] = 0;

    for (i = 0; i < p->size; i++) {
        t = rand() % (i + 1);
        p->pi[i+1] = p->pi[t+1];
        p->pi[t+1] = i+1;
        p->inv_pi[p->pi[i+1]] = i+1;
        p->inv_pi[p->pi[t+1]] = t+1;
    }

    if (p->sig == SIGNED) {
        for (i = 1; i <= p->size; i++) {
            if (rand() % 2)
                p->pi[i] *= -1;
        }
    }
}/*}}}*/

int prm_median(permutation_t *p, int inicio, int fim) {/*{{{*/
    int median, min, m, i;

    median = (fim - inicio + 1) / 2.0;
    if ((fim-inicio+1) % 2)
        median++;

    min = p->size+1;
    for (i = inicio; i <= fim; i++) {
        if (abs(p->pi[i]) < min)
            min = abs(p->pi[i]);
    }

    m = min - 1 + median;

    return m;
}/*}}}*/

void prm_swap(permutation_t *p, int pos1, int pos2) {/*{{{*/
    int aux;
    aux = p->pi[pos1];
    p->pi[pos1] = p->pi[pos2];
    p->pi[pos2] = aux;

    aux = p->pi[pos1];
    if (aux < 0)
        p->inv_pi[-1*aux] = -1*pos1;
    else
        p->inv_pi[aux] = pos1;
    
    aux = p->pi[pos2];
    if (aux < 0)
        p->inv_pi[-1*aux] = -1*pos2;
    else
        p->inv_pi[aux] = pos2;
}/*}}}*/

int prm_exists(permutation_t *p, int elem) {/*{{{*/
    int pos;
    if (elem < 0)
        pos = p->inv_pi[-1*elem];
    else
        pos = p->inv_pi[elem];

    if (pos < 0) {
        if (p->pi[-1*pos] == elem)
            return -1*pos;
    } else {
        if (p->pi[pos] == elem)
            return pos;
    }

    return -1;
}/*}}}*/

void find_extremes_of_strip(permutation_t *p, int elem, int *inicio, int *fim) {/*{{{*/
    int aux;

    /* find beginning and end of the strip that contains elem */
    aux = p->inv_pi[elem];
    while (aux > 1 && abs(p->pi[aux] - p->pi[aux-1]) == 1)
        aux--;
    (*inicio) = aux;

    while (aux < p->size && abs(p->pi[aux] - p->pi[aux+1]) == 1)
        aux++;
    (*fim) = aux;

    /* inicio is the position of the first element of the strip that contains elem */
    if ((*inicio) > (*fim)) {
        aux = (*inicio);
        (*inicio) = (*fim);
        (*fim) = aux;
    }
}/*}}}*/

void find_extremes_of_signed_strip(permutation_t *p, int elem, int *inicio, int *fim) {/*{{{*/
    int aux;

    /* find beginning and end of the strip that contains elem */
    aux = prm_exists(p, elem);
    while (aux > 1 && p->pi[aux] - p->pi[aux-1] == 1)
        aux--;
    (*inicio) = aux;

    while (aux < p->size && p->pi[aux+1] - p->pi[aux] == 1)
        aux++;
    (*fim) = aux;

    /* inicio is the position of the first element of the strip that contains elem */
    if ((*inicio) > (*fim)) {
        aux = (*inicio);
        (*inicio) = (*fim);
        (*fim) = aux;
    }
}/*}}}*/

static void swap(int pi[], int pos1, int pos2) {/*{{{*/
    int tmp;
    tmp = pi[pos1];
    pi[pos1] = pi[pos2];
    pi[pos2] = tmp;
}/*}}}*/

static unsigned long unsigned_rank(int n, permutation_t *p) {/*{{{*/
    int s;
    if (n == 1)
        return 0;

    s = p->pi[n]-1;
    swap(p->pi, n, p->inv_pi[n]);
    swap(p->inv_pi, s+1, n);
    return (s + n * unsigned_rank(n-1, p));
}/*}}}*/

static void unsigned_unrank(int n, unsigned long r, permutation_t *p) {/*{{{*/
    if (n > 0) {
        swap(p->pi, n, r%n+1);
        unsigned_unrank(n-1, r/n, p);
    }
}/*}}}*/

static void get_base2(unsigned long *base2, int n) {/*{{{*/
	int i;

	base2[0] = 1;
	for(i = 1; i < n; i++){
		base2[i] = 2 * base2[i-1];
	}
}/*}}}*/

unsigned long rank(permutation_t *p) {/*{{{*/
    permutation_t pcpy;
    int i;
    unsigned long value, *base2;

    if (p->sig == UNSIGNED) {
        create_permutation(&pcpy, p->size, UNSIGNED);
        prmcpy(&pcpy, p);
        value = unsigned_rank(p->size, &pcpy);
        destroy_permutation(&pcpy);
    } else {
        base2 = Malloc(sizeof(unsigned long) * (p->size+1));
        get_base2(base2, p->size+1);

        create_permutation(&pcpy, p->size, UNSIGNED);
        uns_prmcpy(&pcpy, p);
	    
        value = unsigned_rank(p->size, &pcpy) * base2[p->size];

        for (i = 1; i <= p->size; i++) {
            if (p->pi[i] < 0) {
                value += base2[p->size - i];
            }
        }

        free(base2);
    }
    return value;
}/*}}}*/

void unrank(unsigned long rank, permutation_t *p) {/*{{{*/
    int i;
    unsigned long r, *base2;

    if (p->sig == UNSIGNED) {
        make_identity(p);
        unsigned_unrank(p->size, rank, p);
    } else {
        base2 = Malloc(sizeof(unsigned long) * (p->size+1));
        get_base2(base2, p->size+1);

        r = rank % base2[p->size];
        rank -= r;
        rank = rank / base2[p->size];

        make_identity(p);
        unsigned_unrank(p->size, rank, p);

        for (i = p->size; i >= 1; i--) {
            if (r % 2 == 1) {
                p->pi[i] = (-1) * p->pi[i];
            }
            r = r / 2;
        }

        free(base2);
    }
}/*}}}*/

void prm_separated(permutation_t *p, int *i, int *j, int *n) {/*{{{*/
    int a, b, found, k;

    (*n) = p->size+1;
    (*i) = (*j) = -1;
    for (a = 0; a <= p->size; a++) {
        for (b = a+1; b <= p->size+1; b++) {
            found = 1;
            for (k = 1; k <= a && found; k++)
                if (abs(p->pi[k]) > a) found = 0;
            for (k = a+1; k <= b-1 && found; k++)
                if (abs(p->pi[k]) != k) found = 0;
            for (k = b; k <= p->size && found; k++)
                if (abs(p->pi[k]) < b) found = 0;
            if (found) {
                if (a + p->size-b+1 < (*n) || (a + p->size-b+1 == (*n) &&
                         (*n) == p->size && abs(a-p->size/2) < abs((*i)-p->size/2))) {
                    (*n) = a + p->size-b+1;
                    (*i) = a;
                    (*j) = b;
                }
            }
        }
    }
}/*}}}*/

int is_sorted_interval(permutation_t *p, int ini, int fim, int type) {/*{{{*/
    int i = ini;

    if (type == INC) {
        while (i < fim && p->pi[i+1] == p->pi[i]+1) i++;
        if (i == fim) return 1;
        return 0;
    } else {
        while (i < fim && p->pi[i+1] == p->pi[i]-1) i++;
        if (i == fim) return 1;
        return 0;
    }
}/*}}}*/




