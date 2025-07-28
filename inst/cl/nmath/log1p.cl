// log1p.cl - OpenCL Adaptation of log1p.c
//@provides: Rlog1p
//@depends: nmath, chebyshev
//@includes: nmath



inline double log1p(double x) {
    if (fabs(x) < 1e-4) {
        // series expansion for small x
        return x - 0.5*x*x + x*x*x/3.0;
    }
    return log(1.0 + x);
}