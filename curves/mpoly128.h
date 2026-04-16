#ifndef _MPOLY128_INCLUDE_
#define _MPOLY128_INCLUDE_

#include <stdint.h>

#define MPOLY_MAX_VARS      8

struct mpoly128_term{
    __int128_t c;
    unsigned char d[MPOLY_MAX_VARS];
};

struct mpoly128_node {
    __int128_t c;                               // c MUST be the first element of mpoly_node
    struct mpoly128_node *s;
    int d;
};


struct mpoly128 {
    struct mpoly128_node *t[MPOLY_MAX_VARS];    // array of terms for each level
    struct mpoly128_node *e[MPOLY_MAX_VARS];    // end pointers for each array (length is e[i]-t[i])
    unsigned char d[MPOLY_MAX_VARS];            // max degree in each variable
    __int128_t *c;                              // points to an array of length d[n]+1 coefficients representing level 0 (last level to be evaluated)
    __int128_t *lc;                             // equal to leading coefficient c+d[n]
    int n;                                      // number of variables, must be > 1
};

typedef struct mpoly128 mpoly128_t[1];

void mpoly128_init (mpoly128_t mpoly, struct mpoly128_term *t, int tcnt, int n);
void mpoly128_clear (mpoly128_t mpoly);
void mpoly128_discriminant_init (mpoly128_t mpoly, int d);
void mpoly128_modified_discriminant_init (mpoly128_t mpoly, int d, int vars[]);

static inline void mpoly128_update (mpoly128_t mpoly, int v, __int128_t x)
{
    __int128_t xp[MPOLY_MAX_VARS+1];
    register struct mpoly128_node *t;
    register int i;

    xp[0] = 1; xp[1] = x;
    for ( i = 2 ; i <= mpoly->d[v] ; i++ ) xp[i] = x*xp[i-1];
    if ( v > 1 ) for ( t = mpoly->t[v-1] ; t < mpoly->e[v-1] ; t++ ) t->c = 0;
    else for ( i = 0 ; i <= mpoly->d[0] ; i++ ) mpoly->c[i] = 0;
    for ( t = mpoly->t[v] ; t < mpoly->e[v] ; t++ )  t->s->c += t->c*xp[t->d];
}

// evaluate univariate poly at top level using Horner's method
static inline __int128_t mpoly128_eval (mpoly128_t mpoly, __int128_t x)
{
    register __int128_t *c, *c0, y;
    
    c0 = mpoly->c;  c = mpoly->lc;  y = *c;
    while ( c > c0 ) y = x*y+*(--c);
    return y;
}

// returns current univariate polynomial at level 0
static inline int mpoly128_upoly (mpoly128_t mpoly, __int128_t f[])
{
    register __int128_t *g = mpoly->c ;
    while ( g <= mpoly->lc ) *f++ = *g++;
    return mpoly->d[0];
}

// give direct access to current univariate polynomial at level 0 without copying (caller must not modify!)
static inline __int128_t *mpoly128_upoly_ptr (mpoly128_t mpoly) { return mpoly->c; }
static inline int mpoly128_upoly_deg (mpoly128_t mpoly) { return mpoly->d[0]; }

#endif
