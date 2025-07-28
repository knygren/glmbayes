// fmax2.cl - OpenCL Adaptation of fmax2.c
// @provides: expm1
// @depends: nmath
//@includes: nmath


// fmax2.cl – OpenCL port of R's fmax2.c

inline double fmax2(double x, double y)
{
    if (ISNAN(x)) return y;
    if (ISNAN(y)) return x;
    return (x > y) ? x : y;
}