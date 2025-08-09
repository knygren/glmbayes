#pragma once

#include <string>
#include <vector>

// CPU‐fallback / OpenCL proxy runner for the binomial‐logit prep kernel.
// Mirrors the signature you’ll keep when swapping in the real OpenCL launch.
// 
// kernel_source : full OpenCL program text (helpers + kernel)
// kernel_name   : e.g. "f2_binomial_logit_prep_parallel"
// l1, l2, m1    : dims (observations, coeffs, grid size)
// X_flat        : length = l1 * l2, column‐major design matrix
// B_flat        : length = l2 * m1, column‐major grid of β
// mu_flat       : length = l2, prior means
// P_flat        : length = l2 * l2, prior precision
// alpha_flat    : length = l1, offsets
// qf_flat       : OUT, length = m1, quadratic forms
// xb_flat       : OUT, length = l1 * m1, logistic preps
// progbar       : 0 = no progress, >0 = text bar
void f2_binomial_logit_prep_kernel_runner(
    const std::string&            kernel_source,
    const char*                   kernel_name,
    int                           l1,
    int                           l2,
    int                           m1,
    const std::vector<double>&    X_flat,
    const std::vector<double>&    B_flat,
    const std::vector<double>&    mu_flat,
    const std::vector<double>&    P_flat,
    const std::vector<double>&    alpha_flat,
    std::vector<double>&          qf_flat,
    std::vector<double>&          xb_flat,
    int                           progbar = 0
);