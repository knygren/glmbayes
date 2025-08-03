__kernel void arithmetic_test_v2(const float a,
                                 const float b,
                                 __global float* output) {
    int id = get_global_id(0);
    if (id == 0) {
        output[0] = gpu_add(a, b);              // a + b
        output[1] = subtract(a, b);             // a - b
        output[2] = gpu_multiply(a, b);         // a * b
        output[3] = divide(a, b);               // a / b
        output[4] = mod(a, b);                  // a % b
        output[5] = alpha_form(a, b, a, b);     // Custom operation
    }
}