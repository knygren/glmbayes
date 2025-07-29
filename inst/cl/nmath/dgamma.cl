// dgamma.cl - OpenCL Adaptation of dgamma.c
//@provides: dgamma
//@depends: nmath, lgamma, log1p, dpois
//@includes: nmath, dpq


double dgamma(double x, double shape, double scale, int give_log)
{
    double pr;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
        return x + shape + scale;
#endif
    if (shape < 0 || scale <= 0) ML_WARN_return_NAN;
    if (x < 0)
	return R_D__0;
    if (shape == 0) /* point mass at 0 */
	return (x == 0)? ML_POSINF : R_D__0;
    if (x == 0) {
	if (shape < 1) return ML_POSINF;
	if (shape > 1) return R_D__0;
	/* else */
	return give_log ? -log(scale) : 1 / scale;
    }

    if (shape < 1) {
	pr = dpois_raw(shape, x/scale, give_log);
	return (
	    give_log/* NB: currently *always*  shape/x > 0  if shape < 1:
		     * -- overflow to Inf happens, but underflow to 0 does NOT : */
	    ? pr + (R_FINITE(shape/x)
		    ? log(shape/x)
		    : /* shape/x overflows to +Inf */ log(shape) - log(x))
	    : pr*shape / x);
    }
    /* else  shape >= 1 */
    pr = dpois_raw(shape-1, x/scale, give_log);
    return give_log ? pr - log(scale) : pr/scale;
}
