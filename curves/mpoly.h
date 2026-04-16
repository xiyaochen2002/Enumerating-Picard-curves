#ifndef _MPOLY_INCLUDE_
#define _MPOLY_INCLUDE_

#define MPOLY_MAX_VARS      8

// each term has a long coefficient and an array of variable degrees
struct mpoly_term{
    long c;
    unsigned char d[MPOLY_MAX_VARS];
};

struct mpoly_node {
    long c;                                 // c MUST be the first element of mpoly_node
    struct mpoly_node *s;
    int d;
};


struct mpoly {
    struct mpoly_node *t[MPOLY_MAX_VARS];       // array of terms for each level
    struct mpoly_node *e[MPOLY_MAX_VARS];       // end pointers for each array (length is e[i]-t[i])
    unsigned char d[MPOLY_MAX_VARS];            // max degree in each variable
    long *c;                                    // points to an array of length d[n]+1 coefficients representing level 0 (last level to be evaluated)
    long *lc;                                   // equal to leading coefficient c+d[n]
    int n;                                      // number of variables, must be > 1
};

typedef struct mpoly mpoly_t[1];

void mpoly_init (mpoly_t mpoly, struct mpoly_term *t, int tcnt, int n);
void mpoly_clear (mpoly_t mpoly);

// functions to init mpolys from discriminant polys are defined in mpolydisc.c
// mpoly_modified_discriminant allows caller to zero some variables and choose ordering of nonzero variables
void mpoly_discriminant_init (mpoly_t mpoly, int d);
void mpoly_modified_discriminant_init (mpoly_t mpoly, int d, int vars[]);

static inline void mpoly_update (mpoly_t mpoly, int v, long x)
{
    long xp[MPOLY_MAX_VARS+1];
    register struct mpoly_node *t;
    register int i;

    xp[0] = 1; xp[1] = x;
    for ( i = 2 ; i <= mpoly->d[v] ; i++ ) xp[i] = x*xp[i-1];
    if ( v > 1 ) for ( t = mpoly->t[v-1] ; t < mpoly->e[v-1] ; t++ ) t->c = 0;
    else for ( i = 0 ; i <= mpoly->d[0] ; i++ ) mpoly->c[i] = 0;
    for ( t = mpoly->t[v] ; t < mpoly->e[v] ; t++ )  t->s->c += t->c*xp[t->d];
}

// evaluate univariate poly at top level using Horner's method
static inline long mpoly_eval (mpoly_t mpoly, long x)
{
    register long *c, *c0, y;
    
    c0 = mpoly->c;  c = mpoly->lc;  y = *c;
    while ( c > c0 ) y = x*y+*(--c);
    return y;
}

// copies current univariate polynomial at level 0
static inline int mpoly_upoly (mpoly_t mpoly, long f[])
{
    register long *g = mpoly->c ;
    while ( g <= mpoly->lc ) *f++ = *g++;
    return mpoly->d[0];
}

// give direct access to current univariate polynomial at level 0 without copying (caller must not modify!)
static inline long *mpoly_upoly_ptr (mpoly_t mpoly) { return mpoly->c; }
static inline int mpoly_upoly_deg (mpoly_t mpoly) { return mpoly->d[0]; }


#endif
