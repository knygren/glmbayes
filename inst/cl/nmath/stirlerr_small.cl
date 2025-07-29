// stirlerr_small.cl - OpenCL Adaptation of stirlerr.c (small values only)
//@provides: stirlerr_small
//@depends: nmath, lgamma1p
//@includes: nmath


inline double stirlerr_small(double n) {
    double nn;

        if (n <= 5.25) {
        
        if (n >= 1.0) {
                double l_n = log(n);
                return
//                     lgamma(n)
                     //+ 
                     n * (1.0 - l_n)
                     + ldexp(l_n , -1) /// Temporarily adjusted to exclude -M_LN_2PI

                     + ldexp(l_n - M_LN_2PI, -1)
                     ;
            } else {
                return 
             //lgamma1p(n) 
             - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI
                     ;
            }

        }
        else
        {
        return ML_NAN;
        }

} 















