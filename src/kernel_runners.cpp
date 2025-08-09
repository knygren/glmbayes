


#include <CL/cl.h>
#include <Rcpp.h>
#include <vector>
#include <string>
#include "kernel_runners.h"



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