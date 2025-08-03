// d1mach.cl - OpenCL Adaptation of d1mach.c
// @provides: d1mach, Rf_d1mach
// @depends: nmath
// @include: nmath, Rmath

// @provides: d1mach
// @depends: nmath
// @include: nmath, Rmath

inline double Rf_d1mach(int i)
{
    switch(i) {
        case 1: return DBL_MIN;
        case 2: return DBL_MAX;

        case 3: // = FLT_RADIX ^ -DBL_MANT_DIG
                // for IEEE: = 2^-53 = 1.110223e-16 = 0.5 * DBL_EPSILON
            return 0.5 * DBL_EPSILON;

        case 4: // = FLT_RADIX ^ (1 - DBL_MANT_DIG)
                // for IEEE: = 2^-52 = DBL_EPSILON
            return DBL_EPSILON;

        case 5: return M_LOG10_2;

        default: return 0.0;
    }
}

inline double d1mach(int *i)
{
    return Rf_d1mach(*i);
}