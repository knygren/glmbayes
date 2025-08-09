// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define MAX_L2 64   // upper bound on l2; tune as needed

#pragma OPENCL EXTENSION cl_khr_fp64   : enable   // for double
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf



__kernel void f2_binomial_logit_prep(
    __global const double* X,      // design matrix, size = l1 × l2, col-major: X[k*l1 + i]
    __global const double* B,      // grid of draws,   size = m1 × l2:       B[j*l2 + k]
    __global const double* mu,     // prior mean,       length = l2
    __global const double* P,      // prior precision,  size = l2 × l2, row-major: P[k*l2 + l]
    __global const double* alpha,  // offset vector,    length = l1
    __global double*       qf,     // output quadratics, length = m1
    __global double*       xb,     // output logistic prep, size = m1 × l1
    const int l1,
    const int l2,
    const int m1
) {
    // one work-item per grid point
    int j = get_global_id(0);
    if (j >= m1) return;




    // temporary storage for P*(B_j - mu)
    double tmp[MAX_L2];

    // Compute B_j - mu  and then tmp[k] = ∑_ℓ P[k,ℓ] * (B_j[ℓ] - mu[ℓ])
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    // Quadratic form: 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
    qf[j] = 0.5 * qsum;

    // Logistic “prep”: xb[j,l] = 1 / (1 + exp(alpha[l] + X[*,l] · B_j))
    // note: X is stored so that X[k*l1 + i] = X[k,i]
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        double dot =- alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot -= X[k*l1 + i] * B[j*l2 + k];
        }
        xb[base + i] = 1.0 / (1.0 + exp(dot));
    }
}