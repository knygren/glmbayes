// fmax2.cl - OpenCL Adaptation of fmax2.c
// @provides: expm1
// @depends: nmath
//@includes: nmath

#pragma once

inline double fmax2(double x, double y) {
    if (isnan(x)) return y;
    if (isnan(y)) return x;
    return x > y ? x : y;
}