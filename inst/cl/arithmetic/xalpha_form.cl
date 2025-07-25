inline float xalpha_form(float x, float a, float b, float c) {
  float ax2 = gpu_multiply(a, gpu_multiply(x, x));
  float bx  = gpu_multiply(b, x);
  float sum = gpu_add(ax2, bx);
  return gpu_add(sum, c);
}