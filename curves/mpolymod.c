#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpolymod.h"


static int tcmp (struct mpoly_mod_term *s, struct mpoly_mod_term *t, int n)
{ 
    for ( register int i = 0 ; i < n ; i++ ) { if ( s->d[i] < t->d[i] ) return -1; if ( s->d[i] > t->d[i] ) return 1; }
    return 0;
}

static int tcmp_all (const void *a, const void *b)
    { return tcmp ((struct mpoly_mod_term *)a, (struct mpoly_mod_term *)b, MPOLY_MAX_VARS); }


void mpoly_mod_init (mpoly_mod_t mpoly, struct mpoly_term *t, int tcnt, int n, mm_t m)
{
    struct mpoly_mod_term *s;
    struct mpoly_mod_node *p;
    int o[MPOLY_MAX_VARS];
    register int i, j, k;
    
    assert ( n > 1 && n <= MPOLY_MAX_VARS && tcnt > 1 );
    assert ( (m&0x1) );
    memset (mpoly, 0, sizeof(*mpoly));
    mpoly->m = m;
    mpoly->mi = mm_pinv(m);
    mpoly->R = mm_R (m);
    mpoly->R2 = mm_R2 (mpoly->R, m, mpoly->mi);
    mpoly->s4 = mm_from_ui (24, mpoly->R2, m, mpoly->mi);
    mpoly->s5 = mm_from_ui (120, mpoly->R2, m, mpoly->mi);
    s = malloc (tcnt*sizeof(*s));
    // reduce coefficients modulo m
    for ( i = 0 ; i < tcnt ; i++ ) {
        s[i].c = mm_from_si (t[i].c, mpoly->R2, m, mpoly->mi);
        memcpy (s[i].d, t[i].d, MPOLY_MAX_VARS*sizeof(unsigned char));
    }
    for ( i = 0 ; i < tcnt ; i++ ) {
        for ( j = 0 ; j < n ; j++ ) if ( s[i].d[j] > mpoly->d[j] ) mpoly->d[j] = s[i].d[j];     // update max degree in each variable
        for ( ; j < MPOLY_MAX_VARS ; j++ ) s[i].d[j] = 0;                                       // clear out unused degrees
    }
    qsort (s, tcnt, sizeof(*s), tcmp_all);
    mpoly->t[n-1] = malloc (tcnt*sizeof(*mpoly->t[n-1]));
    mpoly->e[n-1] = mpoly->t[n-1]+tcnt;
    for ( j = 0 ; j < tcnt ; j++ ) { mpoly->t[n-1][j].c = s[j].c; mpoly->t[n-1][j].d = s[j].d[n-1]; mpoly->t[n-1][j].s = 0; }
    for ( i = n-2 ; i  >= 0 ; i-- ) {
        for ( j = k = 1 ; j < tcnt ; j++ ) if ( tcmp (s+j, s+j-1, i+1) ) k++;
        mpoly->t[i] =calloc (k,sizeof(*mpoly->t[i]));
        mpoly->e[i] = mpoly->t[i]+k;
    }
    for ( i = 0 ; i < n ; i++ ) o[i] = 0;
    for ( j = 0 ; ; j++ ) {
        for ( i = n-1 ; i > 0; i-- ) {
            mpoly->t[i][o[i]].d = s[j].d[i];
            mpoly->t[i][o[i]].s = mpoly->t[i-1]+o[i-1];
        }
        mpoly->t[0][o[0]].d = s[j].d[0];
        if ( j+1 == tcnt ) break;
        for ( k = 0 ; k < n ; k++ ) if ( s[j+1].d[k] != s[j].d[k] ) break;
        assert ( k < n );
        o[k++]++;
        while ( k < n ) o[k++]++;
    }
    free (s);
    // the array c effectively replaces the array mpoly->t[0], the difference is that c is dense
    mpoly->c = calloc (mpoly->d[0]+1, sizeof(*mpoly->c));
    mpoly->lc = mpoly->c+mpoly->d[0];
    for ( p = mpoly->t[1] ; p < mpoly->e[1] ; p++ ) p->s = (struct mpoly_mod_node *) (mpoly->c + p->s->d);      // c is first element of struct mpoly_node, so this is ok
    mpoly->n = n;
}

void mpoly_mod_clear (mpoly_mod_t mpoly)
{
    int i;
    
    free (mpoly->c);
    for ( i = 0 ; i < mpoly->n ; i++ ) free (mpoly->t[i]);
}

