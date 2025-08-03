// @test: rmath
// @provides: rmath_test_kernel

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void rmath_test_kernel(__global double* output) {
    // ── Mathematical constants ──
    output[ 0] = M_E;
    output[ 1] = M_LOG2E;
    output[ 2] = M_LOG10E;
    output[ 3] = M_LN2;
    output[ 4] = M_LN10;
    output[ 5] = M_PI;
    output[ 6] = M_2PI;
    output[ 7] = M_PI_2;
    output[ 8] = M_PI_4;
    output[ 9] = M_1_PI;
    output[10] = M_2_PI;
    output[11] = M_2_SQRTPI;
    output[12] = M_SQRT2;
    output[13] = M_SQRT1_2;
    output[14] = M_SQRT_3;
    output[15] = M_SQRT_32;
    output[16] = M_LOG10_2;
    output[17] = M_SQRT_PI;
    output[18] = M_1_SQRT_2PI;
    output[19] = M_SQRT_2dPI;
    output[20] = M_LN_2PI;
    output[21] = M_LN_SQRT_PI;
    output[22] = M_LN_SQRT_2PI;
    output[23] = M_LN_SQRT_PId2;

    // ── Validation macros ──
    const double finite_x    = 1.23;
    const double nan_x       = NAN;
    const double inf_x       = INFINITY;

    output[24] = ISNAN(finite_x)    ? 1.0 : 0.0;  // should be 0
    output[25] = ISNAN(nan_x)       ? 1.0 : 0.0;  // should be 1
    output[26] = R_FINITE(finite_x) ? 1.0 : 0.0;  // should be 1
    output[27] = R_FINITE(inf_x)    ? 1.0 : 0.0;  // should be 0
    output[28] = ML_VALID(finite_x) ? 1.0 : 0.0;  // should be 1
    output[29] = ML_VALID(nan_x)    ? 1.0 : 0.0;  // should be 0

}