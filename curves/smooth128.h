#ifndef _SMOOTH128_INCLUDE_
#define _SMOOTH128_INCLUDE_

#include <stdint.h>
#include <gmp.h>
#include "cstd.h"

void seven_smooth_setup_128 ();

static inline int shash (__int128_t ndiv2) // low bit of odd n has no information in it
{
    register int m = ndiv2;
    return (m&0x1FFFF) | ((m&0x80000)>>2) | ((m&0x800000) >> 5);  // 19 bit concatenation of bits 0..16,19,23 of (n-1)/2
}

// input should be odd and positive, this is not checked
static inline int is_odd_seven_smooth_candidate (__int64_t n)
{
    register int h;
    register __int64_t m;
    extern __uint64_t s128map[8192];
    extern __uint32_t *s128list;

    m = n >> 1;
    h = shash(m);
    if ( ! (s128map[h>>6] & (1UL<<(h&0x3F))) ) return 0;
    return s128list[h] == (__uint32_t)((m&0xFFFFFFFF000000L)>>24);
}

// Input must be positive, dbnd applies to odd-part only
static inline int is_seven_smooth_or_small (__int128_t d, long dbnd)
{
    register __uint64_t x,y;

    if ( ! (__int64_t)d ) d >>= 64;
    d >>= __builtin_ctzll(d);
    if ( d <= dbnd ) return 1;
    if ( ! is_odd_seven_smooth_candidate (d) ) return 0;
    if ( ! (d>>64) ) { // often true and much faster
        x = d;
        for(;;) { y=x/3; if ( 3*y != x ) break; x = y; }
        for(;;) { y=x/5; if ( 5*y != x ) break; x = y; }
        for(;;) { y=x/7; if ( 7*y != x ) break; x = y; }
        return ( x == 1 );
    } else {
        register __int128_t e;
        for(;;) { e=d/3; if ( 3*e != d ) break; d = e; }
        for(;;) { e=d/5; if ( 5*e != d ) break; d = e; }
        for(;;) { e=d/7; if ( 7*e != d ) break; d = e; }
        return ( d == 1 );
    }
}

static inline int is_seven_smooth (__int128_t d) { return is_seven_smooth_or_small (d,9); }

// Input must be positive, dbnd applies to odd-part only
static inline int mpz_is_seven_smooth_or_small (mpz_t D, long Dbnd, mpz_t w[2])
{
    int i = mpz_scan1 (D, 0); mpz_div_2exp (w[0], D, i);
    if ( mpz_cmpabs_ui (w[0], Dbnd) <= 0 ) return 1;
    mpz_set_ui (w[1], 3); mpz_remove (w[0], w[0], w[1]);
    mpz_set_ui (w[1], 5); mpz_remove (w[0], w[0], w[1]);
    mpz_set_ui (w[1], 7); mpz_remove (w[0], w[0], w[1]);
    return ( mpz_cmpabs_ui (w[0], 1) == 0 );
}

static inline int mpz_is_seven_smooth (mpz_t D, long Dbnd, mpz_t w[2]) { return mpz_is_seven_smooth_or_small(D,9,w); }

int ui64_small_radical (uint64_t n);
int ui128_small_radical (uint128_t n);

void smooth20_setup (); // must be called before either ui92_smooth20 or ui123_smooth20 is called

// the _smooth20 functions below return 1 if the input n is 2^20-smooth and has radical <= 2^20 after excluding primes p for which v_p(n) = 10*a+12*b some a,b >= 0
int ui92_smooth20 (uint128_t n);    // requires n < 2^92 (not checked!)
int ui123_smooth20 (uint128_t n);   // requires n < 2^123 (not checked!)

#endif
