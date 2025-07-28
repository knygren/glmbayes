// @test: nmath
// @provides: nmath_test_kernel

__kernel void nmath_test_kernel(__global float* output) {
    // Input sample values for testing ML_VALID and related checks
    const float valid_x     = 0.42f;
    const float nan_x       = ML_NAN;
    const float inf_x       = ML_POSINF;
    const float neg_inf_x   = ML_NEGINF;

    // ── Test constant macros ──
    output[0] = ML_POSINF;       // +∞
    output[1] = ML_NEGINF;       // −∞
    output[2] = ML_NAN;          // NaN

    // ── Test validation macros ──
    output[3] = ISNAN(valid_x)    ? 1.0f : 0.0f;  // should be 0
    output[4] = ISNAN(nan_x)      ? 1.0f : 0.0f;  // should be 1
    output[5] = R_FINITE(valid_x) ? 1.0f : 0.0f;  // should be 1
    output[6] = R_FINITE(inf_x)   ? 1.0f : 0.0f;  // should be 0
    output[7] = ML_VALID(valid_x) ? 1.0f : 0.0f;  // should be 1
    output[8] = ML_VALID(nan_x)   ? 1.0f : 0.0f;  // should be 0
}

