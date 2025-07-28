// @test: dpq
// @provides: dpq_test_kernel

__kernel void dpq_test_kernel(__global float* output) {
    
    const bool log_p = true;
    const bool lower_tail = true;
    
        const float p_valid = 0.3f;
    const float x_valid = 0.7f;
    const float x_lo = -INFINITY;
    const float x_hi = INFINITY;
    const float x_min = 0.2f;
    const float x_max = 0.8f;
    const float f_val  = 1.45f; // example log-density or scaled density value
    const float rf_val = 0.32f; // reference density (might represent fallback or ratio basis)

    
    output[0]  = R_D__0;
    output[1]  = R_D__1;
    output[2]  = R_DT_0;
    output[3]  = R_DT_1;
    output[4]  = R_D_half;

    
    output[5]  = R_D_Lval(p_valid);
    output[6]  = R_D_Cval(p_valid);

    output[7]  = R_D_val(x_valid);
    output[8]  = R_D_qIv(p_valid);
    output[9]  = R_D_exp(x_valid);
    output[10] = R_D_log(p_valid);
    output[11] = R_D_Clog(p_valid);

    output[12] = R_Log1_Exp(x_valid);
    output[13] = R_D_LExp(x_valid);

    output[14] = R_DT_val(x_valid);
    output[15] = R_DT_Cval(x_valid);

    output[16] = R_DT_qIv(p_valid);
    output[17] = R_DT_CIv(p_valid);

    output[18] = R_DT_exp(x_valid);
    output[19] = R_DT_Cexp(x_valid);

    output[20] = R_DT_log(p_valid);
    output[21] = R_DT_Clog(p_valid);
    output[22] = R_DT_Log(p_valid);

    output[23] = R_D_fexp(f_val, x_valid);  // uses x_valid as the second arg
    output[24] = R_D_rtxp(rf_val, x_valid);


}