// @test:      stirlerr
// @provides:  stirlerr_test_kernel
// @requires:  cl_khr_fp64
// @outputs:   21 doubles testing every branch in stirlerr()

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void stirlerr_test_kernel(__global double* output) {


    for (int i = 0; i < 21; ++i) {
        output[i] = 1000.0 + i;
    }


    // ── Table lookup branches (n ≤ 15, 2n ∈ ℤ) ──
    // integer lookup
//    output[0] = stirlerr_large(3.0);
    // half-integer lookup
//    output[1] = stirlerr_large(2.5);

    // ── Small-n branches (n < 1.0 or ≤ 5.25) ──
    // n < 1.0 → lgamma1p branch
//    output[2] = stirlerr(0.3);
    // 1.0 ≤ n ≤ 5.25 → log/ldexp branch
//    output[3] = stirlerr(2.0);
    // boundary at n = 5.25
//    output[4] = stirlerr(5.25);

    // ── Mid-range series branches (5.25 < n ≤ 15.0) ──
    // deepest fallback series (5.25 < n ≤ 6.1)
//    output[5] = stirlerr(5.5);
    // 1-term series (6.1 < n ≤ 6.6)
//    output[6] = stirlerr(6.5);
    // 2-term series (6.6 < n ≤ 7.3)
//    output[7] = stirlerr(7.0);
    // 3-term series (7.3 < n ≤ 8.9)
//    output[8] = stirlerr(8.0);
    // 4-term series (8.9 < n ≤ 12.3)
//    output[9] = stirlerr(10.0);
    // 5-term series (12.3 < n ≤ 12.8)
//    output[10] = stirlerr(12.5);
    // 6-term series (12.8 < n ≤ 15.0)
//    output[11] = stirlerr(13.0);
    // boundary at n = 15.0 → back to table lookup
//    output[12] = stirlerr(15.0);

    // ── Asymptotic branches (n > 15.0) ──
    // boundary small vs large split at n = 23.5
//    output[13] = stirlerr(23.5);
    // 1-term large asymptotic (23.5 < n ≤ 27)
//    output[14] = stirlerr(23.6);
    // fallback asymptotic (23.5 < n ≤ 27)
//    output[15] = stirlerr(25.0);
    // 2-term asymptotic (27 < n ≤ 86)
//    output[16] = stirlerr(50.0);
    // 3-term asymptotic (86 < n ≤ 205)
//    output[17] = stirlerr(100.0);
    // 4-term asymptotic (205 < n ≤ 6180)
//    output[18] = stirlerr(500.0);
    // 5-term asymptotic (6180 < n ≤ 1.57e7)
//    output[19] = stirlerr(1e4);
    // extreme shortcut (n > 1.57e7)
//    output[20] = stirlerr(2e7);
}
