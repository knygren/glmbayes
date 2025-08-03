// gammalims.cl - OpenCL Adaptation of gammalims.c
//@provides: gammalims
//@depends: nmath, d1mach
//@includes: nmath


// gammalims.cl – OpenCL port of R's gammalims.c

inline void gammalims(double *xmin, double *xmax)
{
#ifdef IEEE_754
    *xmin = -170.5674972726612;
    *xmax = 171.61447887182298;
#else
    double alnbig, alnsml, xln, xold;
    int i;

    //In OpenCL, these must be pointers
    int idx = 1;
    int idx2 = 2;

    alnsml = log(d1mach(&idx));

    //alnsml = log(d1mach(1));
    *xmin = -alnsml;

    for (i = 1; i <= 10; ++i) {
        xold = *xmin;
        xln = log(*xmin);
        *xmin -= *xmin * ((*xmin + 0.5) * xln - *xmin - 0.2258 + alnsml) / (*xmin * xln + 0.5);
        if (fabs(*xmin - xold) < 0.005) {
            *xmin = -(*xmin) + 0.01;
            goto find_xmax;
        }
    }

    ML_ERROR(ME_NOCONV, "gammalims");
    *xmin = *xmax = ML_NAN;

find_xmax:
    //In OpenCL, this must be a pointer
//    int idx2 = 2;
    alnbig = log(d1mach(&idx2));
//    alnbig = log(d1mach(2));
    *xmax = alnbig;

    for (i = 1; i <= 10; ++i) {
        xold = *xmax;
        xln = log(*xmax);
        *xmax -= *xmax * ((*xmax - 0.5) * xln - *xmax + 0.9189 - alnbig) / (*xmax * xln - 0.5);
        if (fabs(*xmax - xold) < 0.005) {
            *xmax += -0.01;
            goto done;
        }
    }

    ML_ERROR(ME_NOCONV, "gammalims");
    *xmin = *xmax = ML_NAN;

done:
    *xmin = fmax2(*xmin, -(*xmax) + 1);
#endif
}