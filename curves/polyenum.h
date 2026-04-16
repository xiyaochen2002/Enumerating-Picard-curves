#ifndef _POLYENUM_INCLUDE_
#define _POLYENUM_INCLUDE_

// given fpts[i] = f(x+i*a) for i=0,1 computes Df[i] = (D^i f)(x+a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], and s2=2!*a^2 encode this information
static inline void quadratic_enum_setup (long Df[3], long f[3], long x, long fpts[2], long s2)
{
    Df[0] = fpts[1];                    // f(x+a)
    Df[1] = fpts[1]-fpts[0];                // (Df)(x)
    Df[2] = s2*f[2];                    // (D^2f)
    Df[1] += Df[2];                 // (Df)(x+a)
}

static inline long quadratic_eval (long f[3], long x)
    { return f[0] + x*(f[1] + x*f[2]); }

// given Df[i] = (D^i f) (x) for i=0,..,3 computes (D^i f) (x+1) and returns f(x+1)
static inline long quadratic_enum (long Df[3])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    return Df[0];
}

// given fpts[i] = f(x+i*a) for i=0,..,2 computes Df[i] = (D^i f)(x+2a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], and s3=3!*a^3 encode this information
static inline void cubic_enum_setup (long Df[4], long f[4], long x, long fpts[3], long s3)
{
    Df[0] = fpts[2];                    // f(x+2a)
    Df[1] = fpts[2]-fpts[1];            // (Df)(x+a)
    Df[2] = fpts[2]-2*fpts[1]+fpts[0];  // (D^2f)(x)
    Df[3] = s3*f[3];                    // (D^3f)
    Df[2] += Df[3];                     // (D^2f)(x+a)
    Df[1] += Df[2];                     // (Df)(x+2a)
    Df[2] += Df[3];                     // (D^2f)(x+2a)
}

static inline long cubic_eval (long f[4], long x)
    { return f[0] + x*(f[1] + x*(f[2] + x*f[3])); }

// given Df[i] = (D^i f) (x) for i=0,..,3 computes (D^i f) (x+1) and returns f(x+1)
static inline long cubic_enum (long Df[4])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    return Df[0];
}

// given fpts[i] = f(x+i*a) for i=0,..,2 computes Df[i] = (D^i f)(x+2*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s3=3!*a^3, and s4=4!*a^4 encode this information
static inline void quartic_enum_setup (long Df[5], long f[5], long x, long fpts[3], long s3, long s4)
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

static inline long quartic_eval (long f[5], long x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*f[4]))); }

// given Df[i] = (D^i f) (x) for i=0,..,4 computes (D^i f) (x+1) and returns f(x+1)
static inline long quartic_enum (long Df[5])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    return Df[0];
}

static inline long quintic_eval (long f[6], long x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4]+x*f[5])))); }

// given fpts[i] = f(x+i*a) for i=0,..,3 computes Df[i] = (D^i f)(x+3*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s4=4!*a^4, and s5=5!*a^5 encode this information
static inline void quintic_enum_setup (long Df[6], long f[6], long x, long fpts[4], long s4, long s5)
{
    register long r;
    Df[0] = fpts[3];                    // f(x+3a)
    Df[1] = fpts[3]-fpts[2];            // (Df)(x+2a)
    r = fpts[2]-fpts[1];                // (Df)(x+a)
    Df[2] = Df[1]-r;                    // (D^2f)(x+a)
    r -= fpts[1]-fpts[0];               // (D^2f)(x)
    Df[3] = Df[2]-r;                    // (D^3f)(x)
    Df[5] = s5*f[5];                    // (D^5f)(x) (constant)
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
static inline long quintic_enum (long Df[6])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    return Df[0];
}

static inline long sextic_eval (long f[6], long x)
    { return f[0] + x*(f[1] + x*(f[2] + x*(f[3] + x*(f[4]+x*(f[5]+x*f[6]))))); }

// given fpts[i] = f(x+i*a) for i=0,..,4 computes Df[i] = (D^i f)(x+4*a), where (D^0 f) = f and (D^i f) (x) = (D^(i-1) f) (x+a) - (D^(i-1) f) (x) for i > 0.
// note that this code is independent of a, the inputs fpts[], s4=4!*a^5, and s6=6!*a^6 encode this information
static inline void sextic_enum_setup (long Df[7], long f[7], long x, long fpts[5], long s5, long s6)
{
    register long r,rr;
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
static inline long sextic_enum (long Df[7])
{
    Df[0] += Df[1];
    Df[1] += Df[2];
    Df[2] += Df[3];
    Df[3] += Df[4];
    Df[4] += Df[5];
    Df[5] += Df[6];
    return Df[0];
}

#endif
