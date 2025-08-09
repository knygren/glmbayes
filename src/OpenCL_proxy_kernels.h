#pragma once

#include <vector>

namespace OpenCLProxy {

// CPU‐fallback proxy for your OpenCL prep kernel.
// Mirrors signature you’ll use for the real kernel launch.
//
// Inputs:
//   X_flat     – length = l1 * l2
//   B_flat     – length = l2 * m1
//   mu_flat    – length = l2
//   P_flat     – length = l2 * l2
//   alpha_flat – length = l1
//   l1, l2, m1 – dimensions
//
// Outputs:
//   qf_flat – length = m1
//   xb_flat – length = l1 * m1
//
// progbar = 0 (no progress) or 1 (show simple text bar).
void f2_binomial_logit_prep_kernel_proxy(
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    int l1,
    int l2,
    int m1,
    std::vector<double>& qf_flat,
    std::vector<double>& xb_flat,
    int progbar = 0
);

} // namespace OpenCLProxy