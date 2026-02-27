#pragma once

#include "mpoly.h"
#include "cstd.h"

// Impelmentation of monomial tree for efficient amortized evaluation of multi-variate polynomials, as described in section 3.2 of https://arxiv.org/abs/1602.03715

// mpoly_node represents a node in the tree
struct mpoly128_node {
    int128_t c;                                 // integer coefficient
    int d;                                      // degree of this node in the variable of the current level
    int p;                                      // offset of parent
};


// tree datastructure created from an array of mpoly_terms
struct mpoly128 {
    struct mpoly128_node *b;                    // base pointer for all node storage (one block of memory)
    struct mpoly128_node *t[MPOLY_MAX_VARS+1];  // pointers to nodes at each level, t[i+1] points to the end of level t[i] (so t[n] is end of b array)
    unsigned char d[MPOLY_MAX_VARS+1];            // degree of poly in each variable
    int128_t c[MPOLY_MAX_DEG+1];                // array of d[0]+1 coefficients representing level 0 (last level to be evaluated)
    int128_t *lc;                               // points to leading coefficient c+d[0]
    long x[MPOLY_MAX_DEG+1];                    // instantiated values x[n-1], x[n-2], ..., x[m]
    int m;                                      // next level to update (can increase in jumps but decrease only in increments)
    int n;                                      // number of variables/levels
    int o;                                      // level where fan-in starts
    int L;                                      // number of terms = number of leaves
    int T;                                      // number of nodes
};

typedef struct mpoly128 mpoly128_t[1];

// degs is an optional paramter, if specified it is a list of the maximum degrees in each variable and indicates that terms is sorted
void mpoly128_init (mpoly128_t mpoly, struct mpoly_term *terms, int num_terms, int vars, unsigned char *degs);
void mpoly128_clear (mpoly128_t mpoly);

// create a new monomial tree that is a subtree of one that has been partially instantiated -- vars is the number of variables/levels in the new tree
void mpoly128_copy_subtree (mpoly128_t new_mpoly, mpoly128_t mpoly);
// same as copy, but in-place
void mpoly128_prune (mpoly128_t mpoly);

// functions to init mpolys from discriminant polys are defined in mpolydisc.c
// mpoly_modified_discriminant allows caller to zero some variables and choose ordering of nonzero variables
void mpoly128_discriminant_init (mpoly128_t mpoly, int d);
void mpoly128_modified_discriminant_init (mpoly128_t mpoly, int d, int vars[]);

// functions to init mpolys form plane quartic discriminant polys defined in mpolydisc.c
void mpoly128_plane_quartic_discriminant_init (mpoly128_t mpoly);
void mpoly128_plane_quartic_3to1_discriminant_init (mpoly128_t mpoly);
void mpoly128_geometric_hyperelliptic_discriminant_init (mpoly128_t mpoly, int n);

// Update the level of the tree corresponding to the variable v > 0
static inline void _mpoly128_update (mpoly128_t mpoly, int v, int128_t xp[MPOLY_MAX_DEG+1])
{
    struct mpoly128_node *t;

    if ( v == 1 ) {
        memset(mpoly->c,0,(mpoly->d[0]+1)*sizeof(*mpoly->c));
        for ( t = mpoly->t[1] ; t < mpoly->t[2] ; t++ ) mpoly->c[t->p] += t->c*xp[t->d];
    } else if ( v > mpoly->o ) {                                // no fan-in, so no need to sum
        for ( t = mpoly->t[v] ; t < mpoly->t[v+1] ; t++ ) mpoly->b[t->p].c = t->c*xp[t->d];
    } else {
        for ( t = mpoly->t[v-1] ; t < mpoly->t[v] ; t++ ) t->c = 0;
        for ( t = mpoly->t[v] ; t < mpoly->t[v+1] ; t++ ) mpoly->b[t->p].c += t->c*xp[t->d];
    }
    mpoly->x[v] = xp[1];
    mpoly->m = v;
}

static inline void mpoly128_update (mpoly128_t mpoly, int v, int128_t x)
{
    if ( v >= mpoly->m && mpoly->x[v] == x ) return;            // nothing to do

    int128_t xp[MPOLY_MAX_DEG+1];
    xp[0] = 1; xp[1] = x;
    for ( int i = 2 ; i <= mpoly->d[v] ; i++ ) xp[i] = x*xp[i-1];
    _mpoly128_update (mpoly, v, xp);
}

// evaluate univariate poly at top level using Horner's method
static inline int128_t mpoly128_eval (mpoly128_t mpoly, int128_t x)
{
    int128_t *c0 = mpoly->c, *c = mpoly->lc, y = *c--;
    while ( c >= c0 ) y = x*y + *c--;
    return y;
}

// copies current univariate polynomial at level 0
static inline int mpoly128_upoly (mpoly128_t mpoly, int128_t f[])
{
    int128_t *g = mpoly->c;
    while ( g <= mpoly->lc ) *f++ = *g++;
    return mpoly->d[0];
}

// give direct access to current univariate polynomial at level 0 without copying (caller must not modify!)
static inline int128_t *mpoly128_upoly_ptr (mpoly128_t mpoly) { return mpoly->c; }
static inline int mpoly128_upoly_deg (mpoly128_t mpoly) { return mpoly->d[0]; }

void mpoly128_print (mpoly128_t mpoly, int v);
