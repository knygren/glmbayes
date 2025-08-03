// chebyshev.cl - OpenCL Adaptation of chebyshev.c
//@provides: chebyshev_eval, chebyshev_init
//@depends: nmath
//@includes: nmath


inline int chebyshev_init(double *dos, int nos, double eta)
{
    int i, ii;
    double err;

    if (nos < 1)
	return 0;

    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
	i = nos - ii;
	err += fabs(dos[i]);
	if (err > eta) {
	    return i;
	}
    }
    return i;
}


inline double chebyshev_eval(double x, const double *a, const int n)
{
static const char fname[] = "chebyshev_eval";
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_WARN_return_NAN;

    if (x < -1.1 || x > 1.1) ML_WARN_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}
