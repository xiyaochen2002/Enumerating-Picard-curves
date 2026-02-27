#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntutil.h"
#include "b32.h"
#include "b64.h"
#include "m64.h"
#include "cstd.h"
#include "smooth128.h"

#define LOG3    1.58496250072115618145373894394781650876L
#define LOG5    2.32192809488736234787031942948939017587L
#define LOG7    2.80735492205760410744196931723183080864L

// build a lookup table for the 35727 7-smooth odd numbers less than 2^127
// The 19-bit hash  2^15*(clzll(n) >> 3) + ((n + n>>1) & 0x7FFF) (note that the second summand only depends on 16 low bits)
// We thus obtain a 64K bitmap that should fit in L1-cache that indexes a 2MB array of 32-bit ints (taken from bits 17-48) that should mostly fit in L2
// A match does not guarantee smoothness, but it reduces the probability of non-smoothness by a factor of 2^32
// false positives will still dominate of course, but it makes the remaining cost of smoothness testing negligible.
__uint64_t s128map[8192];
__uint32_t *s128list;   // allocated dynamically

void seven_smooth_setup_128 ()
{
    register __int128_t x, x3, x5, x7;
    register __int64_t m;
    register int e3, e5, e7, h, n;
    char buf[64];

    n = 0;
    memset(s128map,0,sizeof(s128map));
    s128list = malloc ((1<<19)*sizeof(*s128list));
    for ( x3 = 1, e3 = 0 ; e3*LOG3 < 127.0 ; e3++, x3 *= 3 ) {
        for ( x5 = 1, e5 = 0 ; e3*LOG3+e5*LOG5 < 127.0 ; e5++, x5 *=5 ) {
            for ( x7 = 1, e7 = 0 ; e3*LOG3+e5*LOG5+e7*LOG7 < 127.0 ; e7++, x7 *= 7 ) {
                x = x3*x5*x7;
                h = shash(x>>1);  m = (1UL<<(h&0x3F));
                if ( (s128map[h>>6] & m) ) { fprintf(stderr,"collision at h=%d, x=%s, e3=%d, e5=%d, e7=%d\n", h, itoa128(buf, x), e3, e5, e7); exit(0); }
                s128map[h>>6] |= m;
                s128list[h] = (__uint32_t)(((x>>1)&0xFFFFFFFF000000L)>>24);
                n++;
            }
        }
    }
    assert (n == 35727);
}
