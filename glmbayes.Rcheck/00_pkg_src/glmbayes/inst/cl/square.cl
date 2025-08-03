

__kernel void square(
    __global const float* input,
    __global float* output,
    const int n
) {
    int id = get_global_id(0);
    if (id < n) {
        output[id] = input[id] * input[id];
    }
}