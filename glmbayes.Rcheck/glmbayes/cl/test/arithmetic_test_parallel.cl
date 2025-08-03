__kernel void arithmetic_test_parallel(__global float* a_vec,
                                       __global float* b_vec,
                                       __global float* output) {
    int gid = get_global_id(0);  // One work-item per column

    float a = a_vec[gid];
    float b = b_vec[gid];

    int base = gid * 6;  // Each row is 6 output values

    output[base + 0] = gpu_add(a, b);
    output[base + 1] = subtract(a, b);
    output[base + 2] = gpu_multiply(a, b);
    output[base + 3] = divide(a, b);
    output[base + 4] = mod(a, b);
    output[base + 5] = alpha_form(a, b, a, b);
}