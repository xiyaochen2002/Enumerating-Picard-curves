#ifndef _POLYENUM_INCLUDE_
#define _POLYENUM_INCLUDE_

// given fpts[i] = f(x+i*a) for i=0,1 computes Df[i] = (D^i f)(x+a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], and s2=2!*a^2 encode this information
static inline void quadratic_enum_setup (__int128_t Df[3], __int128_t f[3], __int128_t x, __int128_t fpts[2], __int128_t s2)
{
    Df[0] = fpts[1];                // f(x+a)
    Df[1] = fpts[1]-fpts[0];        // (Df)(x)
    Df[2] = s2*f[2];                // (D^2f)
    Df[1] += Df[2];                 // (Df)(x+a)
}

static inline __int128_t quadratic_eval (__int128_t f[3], __int128_t x)
    { return f[0] + x*(f[1] + x*f[2]); }

// given Df[i] = (D^i f) (x) for i=0,..,3 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t quadratic_enum (__int128_t Df[3])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    return Df[0];
}

// given fpts[i] = f(x+i*a) for i=0,..,2 computes Df[i] = (D^i f)(x+2a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], and s3=3!*a^3 encode this information
static inline void cubic_enum_setup (__int128_t Df[4], __int128_t f[4], __int128_t x, __int128_t fpts[3], __int128_t s3)
{
    Df[0] = fpts[2];                    // f(x+2a)
    Df[1] = fpts[2]-fpts[1];            // (Df)(x+a)
    Df[2] = fpts[2]-2*fpts[1]+fpts[0];  // (D^2f)(x)
    Df[3] = s3*f[3];                    // (D^3f)
    Df[2] += Df[3];                     // (D^2f)(x+a)
    Df[1] += Df[2];                     // (Df)(x+2a)
    Df[2] += Df[3];                     // (D^2f)(x+2a)
}

static inline __int128_t cubic_eval (__int128_t f[4], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*f[3])); }

// given Df[i] = (D^i f) (x) for i=0,..,3 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t cubic_enum (__int128_t Df[4])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    return Df[0];
}

// given fpts[i] = f(x+i*a) for i=0,..,2 computes Df[i] = (D^i f)(x+2*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s3=3!*a^3, and s4=4!*a^4 encode this information
static inline void quartic_enum_setup (__int128_t Df[5], __int128_t f[5], __int128_t x, __int128_t fpts[3], __int128_t s3, __int128_t s4)
{
    Df[0] = fpts[2];                        // f(x+2a)
    Df[1] = fpts[2]-fpts[1];                // (Df)(x+a)
    Df[2] = fpts[2]-2*fpts[1]+fpts[0];      // (D^2f)(x)
    Df[4] = s4*f[4];                        // D^4f
    Df[3] = s3*(4*f[4]*x+f[3])+3*Df[4]/2;   // (D^3f)(x)
    Df[2] += Df[3];                         // (D^2f)(x+a)
    Df[1] += Df[2];                         // (Df)(x+2a)
    Df[3] += Df[4];                         // (D^3f)(x+a)
    Df[2] += Df[3];                         // (D^2f)(x+2a)
    Df[3] += Df[4];                         // (D^3f)(x+2a)
}

static inline __int128_t quartic_eval (__int128_t f[5], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*f[4]))); }

// given Df[i] = (D^i f) (x) for i=0,..,4 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t quartic_enum (__int128_t Df[5])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    return Df[0];
}

static inline __int128_t quintic_eval (__int128_t f[6], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4]+x*f[5])))); }

// given fpts[i] = f(x+i*a) for i=0,..,3 computes Df[i] = (D^i f)(x+3*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s4=4!*a^4, and s5=5!*a^5 encode this information
static inline void quintic_enum_setup (__int128_t Df[6], __int128_t f[6], __int128_t x, __int128_t fpts[4], __int128_t s4, __int128_t s5)
{
    register __int128_t r;
    Df[0] = fpts[3];                    // f(x+3a)
    Df[1] = fpts[3]-fpts[2];            // (Df)(x+2a)
    r = fpts[2]-fpts[1];
    Df[2] = Df[1]-r;                    // (D^2f)(x+a)
    r -= fpts[1]-fpts[0];
    Df[3] = Df[2]-r;                    // (D^3f)(x)
    Df[5] = s5*f[5];                    // D^5f
    Df[4] = s4*(5*f[5]*x+f[4])+2*Df[5]; // (D^4f)(x)
    Df[3] += Df[4];                     // (D^3f)(x+a)
    Df[2] += Df[3];                     // (D^2f)(x+2a)
    Df[1] += Df[2];                     // (Df)(x+3a)
    Df[4] += Df[5];                     // (D^4f)(x+a)
    Df[3] += Df[4];                     // (D^3f)(x+2a)
    Df[2] += Df[3];                     // (D^2f)(x+3a)
    Df[4] += Df[5];                     // (D^4f)(x+2a)
    Df[3] += Df[4];                     // (D^3f)(x+3a)
    Df[4] += Df[5];                     // (D^4f)(x+3a)
}

// given Df[i] = (D^i f) (x) for i=0,..,5 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t quintic_enum (__int128_t Df[6])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    return Df[0];
}


static inline __int128_t sextic_eval (__int128_t f[6], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4]+x*(f[5]+x*f[6]))))); }

// given fpts[i] = f(x+i*a) for i=0,..,4 computes Df[i] = (D^i f)(x+4*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s4=4!*a^5, and s6=6!*a^6 encode this information
static inline void sextic_enum_setup (__int128_t Df[7], __int128_t f[7], __int128_t x, __int128_t fpts[5], __int128_t s5, __int128_t s6)
{
    register __int128_t r,rr;
    Df[0] = fpts[4];                        // f(x+4a)
    Df[1] = fpts[4]-fpts[3];                // (D^1f)(x+3a)
    r = fpts[3]-fpts[2];                    // (D^1f)(x+2a)
    Df[2] = Df[1]-r;                        // (D^2f)(x+2a)
    rr = fpts[2]-fpts[1];                   // (D^1f)(x+a)
    r -= rr;                                // (D^2f)(x+a)
    Df[3] = Df[2]-r;                        // (D^3f)(x+a)
    rr -= fpts[1]-fpts[0];                  // (D^2f)(x)
    r -= rr;                                // (D^3f)(x)
    Df[4] = Df[3]-r;                        // (D^4f)(x)
    Df[6] = s6*f[6];                        // (D^6f)(x) (constant)
    Df[5] = s5*(6*f[6]*x+f[5])+5*Df[6]/2;   // (D^5f)(x)
    Df[4] += Df[5];                         // (D^4f)(x+a)
    Df[3] += Df[4];                         // (D^3f)(x+2a)
    Df[2] += Df[3];                         // (D^2f)(x+3a)
    Df[1] += Df[2];                         // (D^1f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+a)
    Df[4] += Df[5];                         // (D^4f)(x+2a
    Df[3] += Df[4];                         // (D^3f)(x+3a)
    Df[2] += Df[3];                         // (D^2f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+2a)
    Df[4] += Df[5];                         // (D^4f)(x+3a)
    Df[3] += Df[4];                         // (D^3f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+3a)
    Df[4] += Df[5];                         // (D^4f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+4a)
}

// given Df[i] = (D^i f) (x) for i=0,..,5 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t sextic_enum (__int128_t Df[7])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    Df[5] += Df[6];
    return Df[0];
}


static inline __int128_t septic_eval (__int128_t f[8], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4] + x*(f[5] + x*(f[6] +x*f[7])))))); }

// given fpts[i] = f(x+i*a) for i=0,..,5 computes Df[i] = (D^i f)(x+5*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s6=6!*a^6, and s7=7!*a^7 encode this information
static inline void septic_enum_setup (__int128_t Df[8], __int128_t f[8], __int128_t x, __int128_t fpts[6], __int128_t s6, __int128_t s7)
{
    register __int128_t r,rr,rrr;
    Df[0] = fpts[5];                        // f(x+5a)
    Df[1] = fpts[5]-fpts[4];                // (D^1f)(x+4a)
    r = fpts[4]-fpts[3];                    // (D^1f)(x+3a)
    Df[2] = Df[1]-r;                        // (D^2f)(x+3a)
    rr = fpts[3]-fpts[2];                   // (D^1f)(x+2a)
    r -= rr;                                // (D^2f)(x+2a)
    Df[3] = Df[2]-r;                        // (D^3f)(x+2a)
    rrr = fpts[2]-fpts[1];                  // (D^1f)(x+a)
    rr -= rrr;                              // (D^2f)(x+a)
    r -= rr;                                // (D^3f)(x+a)
    Df[4] = Df[3]-r;                        // (D^4f)(x+a)
    rrr -= fpts[1]-fpts[0];                 // (D^2f)(x)
    rr -= rrr;                              // (D^3f)(x)
    r -= rr;                                // (D^4f)(x)
    Df[5] = Df[4]-r;                        // (D^5f)(x)
    Df[7] = s7*f[7];                        // (D^7f)(x) (constant)
    Df[6] = s6*(7*f[7]*x+f[6])+3*Df[7];     // (D^6f)(x)
    Df[5] += Df[6];                         // (D^5f)(x+a)
    Df[4] += Df[5];                         // (D^4f)(x+2a)
    Df[3] += Df[4];                         // (D^3f)(x+3a)
    Df[2] += Df[3];                         // (D^2f)(x+4a)
    Df[1] += Df[2];                         // (D^1f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+a)
    Df[5] += Df[6];                         // (D^5f)(x+2a)
    Df[4] += Df[5];                         // (D^4f)(x+3a
    Df[3] += Df[4];                         // (D^3f)(x+4a)
    Df[2] += Df[3];                         // (D^2f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+2a)
    Df[5] += Df[6];                         // (D^5f)(x+3a)
    Df[4] += Df[5];                         // (D^4f)(x+4a)
    Df[3] += Df[4];                         // (D^3f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+3a)
    Df[5] += Df[6];                         // (D^5f)(x+4a)
    Df[4] += Df[5];                         // (D^4f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+5a)
}

// given Df[i] = (D^i f) (x) for i=0,..,7 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t septic_enum (__int128_t Df[8])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    Df[5] += Df[6];
    Df[6] += Df[7];
    return Df[0];
}


static inline __int128_t octic_eval (__int128_t f[9], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4] + x*(f[5] + x*(f[6] +x*(f[7]+x*f[8]))))))); }

// given fpts[i] = f(x+i*a) for i=0,..,6 computes Df[i] = (D^i f)(x+6*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s7=7!*a^7, and s8=8!*a^8 encode this information
static inline void octic_enum_setup (__int128_t Df[9], __int128_t f[9], __int128_t x, __int128_t fpts[7], __int128_t s7, __int128_t s8)
{
    register __int128_t r1,r2,r3,r4;
    Df[0] = fpts[6];                        // f(x+6a)
    Df[1] = fpts[6]-fpts[5];                // (D^1f)(x+5a)
    r1 = fpts[5]-fpts[4];                   // (D^1f)(x+4a)
    Df[2] = Df[1]-r1;                       // (D^2f)(x+4a)
    r2 = fpts[4]-fpts[3];                   // (D^1f)(x+3a)
    r1 -= r2;                               // (D^2f)(x+3a)
    Df[3] = Df[2]-r1;                       // (D^3f)(x+3a)
    r3 = fpts[3]-fpts[2];                   // (D^1f)(x+2a)
    r2 -= r3;                               // (D^2f)(x+2a)
    r1 -= r2;                               // (D^3f)(x+2a)
    Df[4] = Df[3]-r1;                       // (D^4f)(x+2a)
    r4 = fpts[2]-fpts[1];                   // (D^1f)(x+a)
    r3 -= r4;                               // (D^2f)(x+a)
    r2 -= r3;                               // (D^3f)(x+a)
    r1 -= r2;                               // (D^4f)(x+a)
    Df[5] = Df[4]-r1;                       // (D^5f)(x+a)
    r4 -= fpts[1]-fpts[0];                  // (D^2f)(x)
    r3 -= r4;                               // (D^3f)(x)
    r2 -= r3;                               // (D^4f)(x)
    r1 -= r2;                               // (D^5f)(x)
    Df[6] = Df[5]-r1;                       // (D^6f)(x)
    Df[8] = s8*f[8];                        // (D^8f)(x) (constant)
    Df[7] = s7*(8*f[8]*x+f[7])+7*Df[8]/2;   // (D^7f)(x)
    Df[6] += Df[7];                         // (D^6f)(x+a)
    Df[5] += Df[6];                         // (D^5f)(x+2a)
    Df[4] += Df[5];                         // (D^4f)(x+3a)
    Df[3] += Df[4];                         // (D^3f)(x+4a)
    Df[2] += Df[3];                         // (D^2f)(x+5a)
    Df[1] += Df[2];                         // (D^1f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+a)
    Df[6] += Df[7];                         // (D^6f)(x+2a)
    Df[5] += Df[6];                         // (D^5f)(x+3a)
    Df[4] += Df[5];                         // (D^4f)(x+4a
    Df[3] += Df[4];                         // (D^3f)(x+5a)
    Df[2] += Df[3];                         // (D^2f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+2a)
    Df[6] += Df[7];                         // (D^6f)(x+3a)
    Df[5] += Df[6];                         // (D^5f)(x+4a)
    Df[4] += Df[5];                         // (D^4f)(x+5a)
    Df[3] += Df[4];                         // (D^3f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+3a)
    Df[6] += Df[7];                         // (D^6f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+5a)
    Df[4] += Df[5];                         // (D^4f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+4a)
    Df[6] += Df[7];                         // (D^6f)(x+5a)
    Df[5] += Df[6];                         // (D^5f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+6a)
}

// given Df[i] = (D^i f) (x) for i=0,..,8 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t octic_enum (__int128_t Df[9])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    Df[5] += Df[6];
    Df[6] += Df[7];
    Df[7] += Df[8];
    return Df[0];
}

static inline __int128_t nonic_eval (__int128_t f[10], __int128_t x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4] + x*(f[5] + x*(f[6] +x*(f[7]+x*(f[8]+x*f[9])))))))); }

// given fpts[i] = f(x+i*a) for i=0,..,7 computes Df[i] = (D^i f)(x+7*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s8=8!*a^8, and s9=9!*a^9 encode this information
static inline void nonic_enum_setup (__int128_t Df[10], __int128_t f[10], __int128_t x, __int128_t fpts[8], __int128_t s8, __int128_t s9)
{
    register __int128_t r1,r2,r3,r4,r5;
    Df[0] = fpts[7];                        // f(x+7a)
    Df[1] = fpts[7]-fpts[6];                // (D^1f)(x+6a)
    r1 = fpts[6]-fpts[5];                   // (D^1f)(x+5a)
    Df[2] = Df[1]-r1;                       // (D^2f)(x+5a)
    r2 = fpts[5]-fpts[4];                   // (D^1f)(x+4a)
    r1 -= r2;                               // (D^2f)(x+4a)
    Df[3] = Df[2]-r1;                       // (D^3f)(x+4a)
    r3 = fpts[4]-fpts[3];                   // (D^1f)(x+3a)
    r2 -= r3;                               // (D^2f)(x+3a)
    r1 -= r2;                               // (D^3f)(x+3a)
    Df[4] = Df[3]-r1;                       // (D^4f)(x+3a)
    r4 = fpts[3]-fpts[2];                   // (D^1f)(x+2a)
    r3 -= r4;                               // (D^2f)(x+2a)
    r2 -= r3;                               // (D^3f)(x+2a)
    r1 -= r2;                               // (D^4f)(x+2a)
    Df[5] = Df[4]-r1;                       // (D^5f)(x+2a)
    r5 = fpts[2]-fpts[1];                   // (D^1f)(x+a)
    r4 -= r5;                               // (D^2f)(x+a)
    r3 -= r4;                               // (D^3f)(x+a)
    r2 -= r3;                               // (D^4f)(x+a)
    r1 -= r2;                               // (D^5f)(x+a)
    Df[6] = Df[5]-r1;                       // (D^6f)(x+a)
    r5 -= fpts[1]-fpts[0];                  // (D^2f)(x)
    r4 -= r5;                               // (D^3f)(x)
    r3 -= r4;                               // (D^4f)(x)
    r2 -= r3;                               // (D^5f)(x)
    r1 -= r2;                               // (D^6f)(x)
    Df[7] = Df[6]-r1;                       // (D^7f)(x)
    Df[9] = s9*f[9];                        // (D^9f)(x) (constant)
    Df[8] = s8*(9*f[9]*x+f[8])+4*Df[9];     // (D^8f)(x)
    Df[7] += Df[8];                         // (D^7f)(x+a)
    Df[6] += Df[7];                         // (D^6f)(x+2a)
    Df[5] += Df[6];                         // (D^5f)(x+3a)
    Df[4] += Df[5];                         // (D^4f)(x+4a)
    Df[3] += Df[4];                         // (D^3f)(x+5a)
    Df[2] += Df[3];                         // (D^2f)(x+6a)
    Df[1] += Df[2];                         // (D^1f)(x+7a)
    Df[8] += Df[9];                         // (D^8f)(x+a)
    Df[7] += Df[8];                         // (D^7f)(x+2a)
    Df[6] += Df[7];                         // (D^6f)(x+3a)
    Df[5] += Df[6];                         // (D^5f)(x+4a)
    Df[4] += Df[5];                         // (D^4f)(x+5a
    Df[3] += Df[4];                         // (D^3f)(x+6a)
    Df[2] += Df[3];                         // (D^2f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+2a)
    Df[7] += Df[8];                         // (D^7f)(x+3a)
    Df[6] += Df[7];                         // (D^6f)(x+4a)
    Df[5] += Df[6];                         // (D^5f)(x+5a)
    Df[4] += Df[5];                         // (D^4f)(x+6a)
    Df[3] += Df[4];                         // (D^3f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+3a)
    Df[7] += Df[8];                         // (D^7f)(x+4a)
    Df[6] += Df[7];                         // (D^6f)(x+5a)
    Df[5] += Df[6];                         // (D^5f)(x+6a)
    Df[4] += Df[5];                         // (D^4f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+4a)
    Df[7] += Df[8];                         // (D^7f)(x+5a)
    Df[6] += Df[7];                         // (D^6f)(x+6a)
    Df[5] += Df[6];                         // (D^5f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+5a)
    Df[7] += Df[8];                         // (D^7f)(x+6a)
    Df[6] += Df[7];                         // (D^6f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+6a)
    Df[7] += Df[8];                         // (D^7f)(x+7a)
    Df[8] += Df[9];                         // (D^f8)(x+7a)
}

// given Df[i] = (D^i f) (x) for i=0,..,8 computes (D^i f) (x+1) and returns f(x+1)
static inline __int128_t nonic_enum (__int128_t Df[9])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    Df[5] += Df[6];
    Df[6] += Df[7];
    Df[7] += Df[8];
    Df[8] += Df[9];
    return Df[0];
}

static inline __int128_t poly_eval (int d, __int128_t f[d+1], __int128_t x)
    { __int128_t y =f[d];  for ( register int i = d-1 ; i >= 0 ; i-- ) y = x*y+f[i]; return y; }

// given fpts[i] = f(x+i*a) for i=0,..,d-2 computes Df[i] = (D^i f)(x+(d-2)*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], sdm1=(d-1)!*a^(d-1), and sd=d!*a^d encode this information
static inline void poly_enum_setup (int d, __int128_t Df[d+1], __int128_t f[d+1], __int128_t x, __int128_t fpts[d-1], __int128_t sdm1, __int128_t sd)
{
    __int128_t r[d-4];
    register int i,j;

    assert (d >= 5);
    Df[0] = fpts[d-2];
    Df[1] = fpts[d-2]-fpts[d-3];
    for ( i = 0 ; i < d-4 ; i++ ) {
        r[i] = fpts[d-3-i] - fpts[d-4-i];
        for ( j = i-1 ; j >= 0 ; j-- ) r[j] -= r[j+1];
        Df[i+2] = Df[i+1]-r[0];
    }
    r[d-5] -= fpts[1] - fpts[0];
    for ( j = d-6 ; j >= 0 ; j-- ) r[j] -= r[j+1];
    Df[d-2] = Df[d-3]-r[0];
    Df[d] = sd*f[d];
    Df[d-1] = sdm1*(d*f[d]*x+f[d-1])+(d-1)*Df[d]/2;
    for ( i = d-2 ; i > 0 ; i-- ) Df[i] += Df[i+1];
    for ( i = 2 ; i < d ; i++ ) for ( j = d-1 ; j >= i ; j-- ) Df[j] += Df[j+1];
}

static inline __int128_t poly_enum (int d, __int128_t Df[d+1])
{
    for ( register int i = 0 ; i < d ; i++ ) Df[i] += Df[i+1];
    return Df[0];
}

#endif
