// expm1.cl - OpenCL Adaptation of expm1.c
// @provides: expm1
// @depends: nmath, log1p
//@includes: nmath




#ifndef HAVE_EXPM1
double expm1(double x)
{
    double y, a = fabs(x);

    if (a < DBL_EPSILON) return x;
    if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

    if (a > 1e-8)
	y = exp(x) - 1;
    else /* Taylor expansion, more accurate in this range */
	y = (x / 2 + 1) * x;

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (log1p (y) - x);
    return y;
}
#endif