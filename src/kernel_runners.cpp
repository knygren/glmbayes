

#ifdef USE_OPENCL
#include "kernel_loader.h"

#ifdef USE_DIRECT_CLH
// we passed “-I…/include/CL -DUSE_DIRECT_CLH”
#include <CL/cl.h>
#else
// normal case on Linux/macOS/Windows
#include <CL/cl.h>
#endif
#endif


//#include <Rcpp.h>
# include <RcppArmadillo.h>
#include <vector>
#include <string>


#ifdef USE_OPENCL

void f2_binomial_logit_prep_kernel_runner(
    const std::string& kernel_source,
    const char*        kernel_name,
    int                l1,
    int                l2,
    int                m1,
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    std::vector<double>&       qf_flat,   // output, length = m1
    std::vector<double>&       xb_flat,   // output, length = l1*m1
    int progbar 
) {
  if ((int)X_flat.size() != l1*l2 ||
      (int)B_flat.size() != l2*m1 ||
      (int)mu_flat.size() != l2 ||
      (int)P_flat.size() != l2*l2 ||
      (int)alpha_flat.size() != l1) {
    throw std::runtime_error("Input flat-vector sizes mismatch dimensions.");
  }
  
  // Initialize outputs
  qf_flat.assign(m1, 0.0);
  xb_flat.assign((size_t)l1*m1, 0.0);
  
  cl_int status;
  // 🔍 Platform & Device
  cl_platform_id platform;
  cl_device_id   device;
  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 🌐 Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  cl_queue_properties props[] = {0};
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  // 📦 Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t src_len = kernel_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem bufX = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*X_flat.size(), (void*)X_flat.data(), &status);
  cl_mem bufB = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*B_flat.size(), (void*)B_flat.data(), &status);
  cl_mem bufMu= clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*mu_flat.size(), (void*)mu_flat.data(), &status);
  cl_mem bufP = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*P_flat.size(), (void*)P_flat.data(), &status);
  cl_mem bufA = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*alpha_flat.size(), (void*)alpha_flat.data(), &status);
  
  cl_mem bufQF= clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                               sizeof(double)*qf_flat.size(), nullptr, &status);
  cl_mem bufXB= clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                               sizeof(double)*xb_flat.size(), nullptr, &status);
  
  // ⏫ Transfer Data to Device (already via COPY_HOST_PTR for inputs)
  
  // 🗳️ Set Kernel Args
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufXB);
  clSetKernelArg(kernel, arg++, sizeof(int), &l1);
  clSetKernelArg(kernel, arg++, sizeof(int), &l2);
  clSetKernelArg(kernel, arg++, sizeof(int), &m1);
  
  // 🚀 Launch Parallel Kernel
  size_t global = m1;
  
//  size_t global = m1;
  
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Output
  status = clEnqueueReadBuffer(queue, bufQF, CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(), qf_flat.data(),
                               0, nullptr, nullptr);
  status = clEnqueueReadBuffer(queue, bufXB, CL_TRUE, 0,
                               sizeof(double)*xb_flat.size(), xb_flat.data(),
                               0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(bufX);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufXB);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}
#endif



#ifdef USE_OPENCL

void f2_binomial_logit_prep_grad_kernel_runner(
    const std::string&        kernel_source,  // your .cl contents
    const char*               kernel_name,    // "f2_binomial_logit_prep_grad"
    int                       l1,             // nobs
    int                       l2,             // ncoef
    int                       m1,             // ngrids
    const std::vector<double>& X_flat,        // length = l1*l2
    const std::vector<double>& B_flat,        // length = m1*l2
    const std::vector<double>& mu_flat,       // length = l2
    const std::vector<double>& P_flat,        // length = l2*l2
    const std::vector<double>& alpha_flat,    // length = l1
    const std::vector<double>& y_flat,        // length = l1
    const std::vector<double>& wt_flat,       // length = l1
    std::vector<double>&       qf_flat,       // OUT: length = m1
    std::vector<double>&       xb_flat,       // OUT: length = m1*l1
    std::vector<double>&       grad_flat,     // OUT: length = m1*l2
    int                       progbar 
) {
  // 0) Sanity‐check sizes
  if ((int)X_flat.size()   != l1*l2 ||
      (int)B_flat.size()   != m1*l2 ||
      (int)mu_flat.size()  != l2 ||
      (int)P_flat.size()   != l2*l2 ||
      (int)alpha_flat.size()!= l1 ||
      (int)y_flat.size()   != l1 ||
      (int)wt_flat.size()  != l1) {
    throw std::runtime_error("Input flat‐vector sizes mismatch dimensions.");
  }
  
  // 1) Initialize output buffers
  qf_flat  .assign(m1,           0.0);
  xb_flat  .assign((size_t)l1*m1, 0.0);
  grad_flat.assign((size_t)l2*m1, 0.0);
  
  cl_int status;
  
  // 2) Platform & Device
  cl_platform_id platform;
  cl_device_id   device;
  status  = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs  (platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 3) Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  cl_queue_properties props[] = {0};
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  // 4) Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t      src_len = kernel_source.size();
  cl_program  program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram   (program, 0, nullptr, nullptr, nullptr, nullptr);
  
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 5) Device Buffers
  cl_mem bufX    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*X_flat.size(),   (void*)X_flat.data(),   &status);
  cl_mem bufB    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*B_flat.size(),   (void*)B_flat.data(),   &status);
  cl_mem bufMu   = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*mu_flat.size(),  (void*)mu_flat.data(),  &status);
  cl_mem bufP    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*P_flat.size(),   (void*)P_flat.data(),   &status);
  cl_mem bufA    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*alpha_flat.size(), (void*)alpha_flat.data(), &status);
  cl_mem bufY    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*y_flat.size(),   (void*)y_flat.data(),   &status);
  cl_mem bufW    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*wt_flat.size(),  (void*)wt_flat.data(),  &status);
  
  cl_mem bufQF   = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*qf_flat.size(),   nullptr, &status);
  cl_mem bufXB   = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*xb_flat.size(),   nullptr, &status);
  cl_mem bufGrad = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*grad_flat.size(), nullptr, &status);
  
  // 6) Set Kernel Args (must match the .cl signature exactly)
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufY);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufW);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufXB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufGrad);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l1);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l2);
  clSetKernelArg(kernel, arg++, sizeof(int),    &m1);
  
  // 7) Launch
  size_t global = (size_t)m1;
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  
  // 8) Read back outputs
  status = clEnqueueReadBuffer(queue, bufQF,   CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(),   qf_flat.data(),
                               0, nullptr, nullptr);
  
  status = clEnqueueReadBuffer(queue, bufXB,   CL_TRUE, 0,
                               sizeof(double)*xb_flat.size(),   xb_flat.data(),
                               0, nullptr, nullptr);
  
  status = clEnqueueReadBuffer(queue, bufGrad, CL_TRUE, 0,
                               sizeof(double)*grad_flat.size(), grad_flat.data(),
                               0, nullptr, nullptr);
  
  // 9) Cleanup
  clReleaseMemObject(bufX);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufY);
  clReleaseMemObject(bufW);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufXB);
  clReleaseMemObject(bufGrad);
  
  clReleaseKernel       (kernel);
  clReleaseProgram      (program);
  clReleaseCommandQueue (queue);
  clReleaseContext      (context);
}
#endif
