// expm1.cl - OpenCL Adaptation of expm1.c
// @provides: expm1
// @depends: nmath, log1p
//@includes: nmath


#pragma once

inline double expm1(double x) {
    if (fabs(x) < 1e-5) {
        // series expansion for small x
        return x + 0.5*x*x + x*x*x/6.0;
    }
    return exp(x) - 1.0;
}