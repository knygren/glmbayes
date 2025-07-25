inline float divide(float x, float y) {
    return y != 0.0f ? x / y : 0.0f;  // Simple protection against division by zero
}