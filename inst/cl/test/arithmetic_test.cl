__kernel void arithmetic_test(__global float* output) {
    int id = get_global_id(0);
    if (id == 0) {
        float a = 10.0f;
        float b = 3.0f;
        output[0] = gpu_add(a, b);         // 10 + 3 = 13
        output[1] = subtract(a, b);    // 10 - 3 = 7
        output[2] = gpu_multiply(a, b);    // 10 * 3 = 30
        output[3] = divide(a, b);      // 10 / 3 ≈ 3.333
        output[4] = mod(a, b);         // 10 % 3 = 1
        output[5] = xalpha_form(a, b, a, b);        // 10 % 3 = 1
        
        //output[0] = 111.0f;
        //output[1] = 222.0f;
        //output[2] = 333.0f;
        //output[3] = 444.0f;
        //output[4] = 555.0f;
       // output[5] = 666.0f;
        
    }
}