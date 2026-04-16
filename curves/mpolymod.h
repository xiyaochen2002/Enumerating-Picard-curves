#ifndef _MPOLY_MOD_INCLUDE_
#define _MPOLY_MOD_INCLUDE_

#include "mpoly.h"

// #define MM_BITS             32 // this gives about a 10 percent speedup
#include "mm.h"

// each term has a long coefficient and an array of variable degrees
struct mpoly_mod_term{
    mm_t c;
    unsigned char d[MPOLY_MAX_VARS];
};

struct mpoly_mod_node {
    mm_t c;                                 // c MUST be the first element of mpoly_node
    int d;
    struct mpoly_mod_node *s;
};


struct mpoly_mod {
    struct mpoly_mod_node *t[MPOLY_MAX_VARS];   // array of terms for each level
    struct mpoly_mod_node *e[MPOLY_MAX_VARS];   // end pointers for each array (length is e[i]-t[i])
    unsigned char d[MPOLY_MAX_VARS];            // max degree in each variable
    mm_t *c;                                    // points to an array of length d[n]+1 coefficients representing level 0 (last level to be evaluated)
    mm_t *lc;                                   // equal to leading coefficient c+d[n]
    int n;                                      // number of variables, must be > 1
    mm_t m, mi, R, R2;                          // modulus parameters
    mm_t s4, s5;                                // 4! and 5! in mm
};

typedef struct mpoly_mod mpoly_mod_t[1];

// note that we take in mpoly_term's and convert them to mpoly_mod_term's
void mpoly_mod_init (mpoly_mod_t mpoly, struct mpoly_term *t, int tcnt, int n, mm_t m);
void mpoly_mod_clear (mpoly_mod_t mpoly);

// functions to init mpoly_mods from discriminant polys are defined in mpolydiscmod.c
void mpoly_mod_discriminant_init (mpoly_mod_t mpoly, int d, mm_t m);

static inline void mpoly_mod_update (mpoly_mod_t mpoly, int v, long x)
{
    mm_t xp[MPOLY_MAX_VARS+1];
    register struct mpoly_mod_node *t;
    register int i;

    xp[0] = mpoly->R; xp[1] = mm_from_si(x, mpoly->R2, mpoly->m, mpoly->mi);
    for ( i = 2 ; i <= mpoly->d[v] ; i++ ) xp[i] = mm_mul(xp[1], xp[i-1], mpoly->m, mpoly->mi);
    if ( v > 1 ) for ( t = mpoly->t[v-1] ; t < mpoly->e[v-1] ; t++ ) t->c = 0;
    else for ( i = 0 ; i <= mpoly->d[0] ; i++ ) mpoly->c[i] = 0;
    for ( t = mpoly->t[v] ; t < mpoly->e[v] ; t++ ) {
        t->s->c = mm_addmul(t->s->c, t->c, xp[t->d], mpoly->R, mpoly->m, mpoly->mi);
    }
}

// copies current univariate polynomial at level 0
static inline int mpoly_mod_upoly (mpoly_mod_t mpoly, mm_t f[])
{
    register mm_t *g = mpoly->c;
    while ( g <= mpoly->lc ) *f++ = *g++;
    return mpoly->d[0];
}

// note that input z is not an mm_t
static inline mm_t mod_quintic_eval (mpoly_mod_t mpoly, mm_t f[6], mm_t x)
{
    register mm_t m = mpoly->m;
    register mm_t mi = mpoly->mi;
    register mm_t R = mpoly->R;
    return mm_addmul(f[0],x,mm_addmul(f[1],x,mm_addmul(f[2],x,mm_addmul(f[3],x,mm_addmul(f[4],x,f[5],R,m,mi),R,m,mi),R,m,mi),R,m,mi),R,m,mi);
}

// given fpts[i] = f(x+i*a) for i=0,..,3 computes Df[i] = (D^i f)(x+3*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s4=4!*a^4, and s5=5!*a^5 encode this information
// note that caller needs to convert x to an mm_t!
static inline void mod_quintic_enum_setup (mpoly_mod_t mpoly, mm_t Df[6], mm_t f[6], mm_t x, mm_t fpts[4], mm_t s4, mm_t s5)
{
    register mm_t r;
    register mm_t m = mpoly->m;
    register mm_t mi = mpoly->mi;
    register mm_t R = mpoly->R;

    Df[0] = fpts[3];                        // f(x+3a)
    Df[1] = mm_sub (fpts[3],fpts[2],m);     // (Df)(x+2a)
    r = mm_sub (fpts[2],fpts[1],m);         // (Df)(x+a)
    Df[2] = mm_sub(Df[1],r,m);              // (D^2f)(x+a)
    r = mm_sub(r,mm_sub(fpts[1],fpts[0],m),m); // (D^2f)(x)
    Df[3] = mm_sub(Df[2],r,m);              // (D^3f)(x)
    Df[5] = mm_mul(s5,f[5],m,mi);    // (D^5f)(x) (constant)
    Df[4] = mm_addmul (2*Df[5],s4,mm_addmul(f[4],5*f[5],x,R,m,mi),R,m,mi);     // (D^4f)(x)
    Df[3] = mm_add(Df[3],Df[4],m);          // (D^3f)(x+a)
    Df[2] = mm_add(Df[2],Df[3],m);          // (D^2f)(x+2a)
    Df[1] = mm_add(Df[1],Df[2],m);          // (Df)(x+3a)
    Df[4] = mm_add(Df[4],Df[5],m);          // (D^4f)(x+a)
    Df[3] = mm_add(Df[3],Df[4],m);          // (D^3f)(x+2a)
    Df[2] = mm_add(Df[2],Df[3],m);          // (D^2f)(x+3a)
    Df[4] = mm_add(Df[4],Df[5],m);          // (D^4f)(x+2a)
    Df[3] = mm_add(Df[3],Df[4],m);          // (D^3f)(x+3a)
    Df[4] = mm_add(Df[4],Df[5],m);          // (D^4f)(x+3a)
}

// given Df[i] = (D^i f) (x) for i=0,..,5 computes (D^i f) (x+1) and returns f(x+1)
static inline mm_t mod_quintic_enum (mpoly_mod_t mpoly, mm_t Df[6])
{
    register mm_t m = mpoly->m;
    Df[0] = mm_add(Df[0],Df[1],m);
    Df[1] = mm_add(Df[1],Df[2],m);
    Df[2] = mm_add(Df[2],Df[3],m);
    Df[3] = mm_add(Df[3],Df[4],m);
    Df[4] = mm_add(Df[4],Df[5],m);
    return Df[0];
}


#endif
