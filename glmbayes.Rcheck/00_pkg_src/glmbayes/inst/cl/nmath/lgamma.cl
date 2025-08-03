// lgamma.cl - OpenCL Adaptation of lgamma.c
//@provides: lgammafn_sign, lgammafn
//@depends: nmath, lgammacor, gamma
//@includes: nmath


inline double lgammafn_sign(double x, int *sgn)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax_lgamma = 0.;
    static double dxrel_lgamma = 0.;

    if (xmax_lgamma == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax_lgamma = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel_lgamma = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax_lgamma  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel_lgamma = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax_lgamma  2.5327372760800758e+305
#define dxrel_lgamma 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
	*sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	// No warning: this is the best answer; was  ML_WARNING(ME_RANGE, "lgamma");
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax_lgamma) {
	// No warning: +Inf is the best answer
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_WARN_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel_lgamma) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer; e.g. for  lgamma(1e-7 - 11) */

	ML_WARNING(ME_PRECISION, "lgamma");
    }

    return ans;
}

inline double lgammafn(double x)
{
    return lgammafn_sign(x, NULL);
}
