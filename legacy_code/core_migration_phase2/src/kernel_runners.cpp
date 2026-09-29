

#ifdef USE_OPENCL

#ifdef USE_DIRECT_CLH
#define CL_TARGET_OPENCL_VERSION 300

// we passed “-I…/include/CL -DUSE_DIRECT_CLH”
#include <CL/cl.h>
#else
#define CL_TARGET_OPENCL_VERSION 300

// normal case on Linux/macOS/Windows
#include <CL/cl.h>
#endif
#endif


//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "openclPort.h"
#include "opencl.h"
#include <vector>
#include <string>

using namespace openclPort;
using namespace glmbayes::opencl;

#ifdef USE_OPENCL

namespace glmbayes {

namespace opencl {

void f2_f3_kernel_runner(
    const std::string&        kernel_source,
    const char*               kernel_name,
    int                       l1,
    int                       l2,
    int                       m1,
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    const std::vector<double>& y_flat,
    const std::vector<double>& wt_flat,
    std::vector<double>&       qf_flat,
    std::vector<double>&       grad_flat,
    int                        progbar
) {
  // 0) Sanity-check sizes
  if ((int)X_flat.size()    != l1*l2 ||
      (int)B_flat.size()    != m1*l2 ||
      (int)mu_flat.size()   != l2    ||
      (int)P_flat.size()    != l2*l2 ||
      (int)alpha_flat.size()!= l1    ||
      (int)y_flat.size()    != l1    ||
      (int)wt_flat.size()   != l1) {
    throw std::runtime_error("Input flat-vector sizes mismatch dimensions.");
  }
  
  // 1) Initialize output buffers
  qf_flat.assign(m1, 0.0);
  grad_flat.assign((size_t)l2*m1, 0.0);
  
  cl_int status = 0;
  
  // 2) Platform & Device
  // Rcpp::Rcout << "[runner] P0: before clGetPlatformIDs\n";
  
  cl_platform_id platform = nullptr;
  cl_device_id   device   = nullptr;
  
  // ---- Platform ----
  status = clGetPlatformIDs(1, &platform, nullptr);

  // -1001 = CL_PLATFORM_NOT_FOUND_KHR (not always defined in headers)
  if (status == -1001) {
    throw std::runtime_error(
        "OpenCL error: no OpenCL platforms found (clGetPlatformIDs returned -1001). "
        "Your system does not expose an OpenCL platform."
    );
  }
  if (status != CL_SUCCESS) {
    std::ostringstream msg;
    msg << "OpenCL error: clGetPlatformIDs failed with status " << status << ".";
    throw std::runtime_error(msg.str());
  }
  
  // ---- Device ----
  status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);

  // -9 = CL_DEVICE_NOT_FOUND
  if (status == -9) {
    throw std::runtime_error(
        "OpenCL error: no suitable OpenCL GPU devices found "
        "(clGetDeviceIDs returned -9)."
    );
  }
  if (status != CL_SUCCESS) {
    std::ostringstream msg;
    msg << "OpenCL error: clGetDeviceIDs failed with status " << status << ".";
    throw std::runtime_error(msg.str());
  }
  
  // 3) Context & Queue
  // Rcpp::Rcout << "[runner] P3: before clCreateContext\n";
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  // Rcpp::Rcout << "[runner] P4: after clCreateContext, status=" << status << "\n";
  
  cl_queue_properties props[] = {0};
  // Rcpp::Rcout << "[runner] P5: before clCreateCommandQueueWithProperties\n";
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  // Rcpp::Rcout << "[runner] P6: after clCreateCommandQueueWithProperties, status=" << status << "\n";
  
  // 4) Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t      src_len = kernel_source.size();
  // Rcpp::Rcout << "[runner] P7: before clCreateProgramWithSource, src_len=" << src_len << "\n";
  cl_program  program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  // Rcpp::Rcout << "[runner] P8: after clCreateProgramWithSource, status=" << status << "\n";
  
  // Rcpp::Rcout << "[runner] P9: before clBuildProgram\n";
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  // Rcpp::Rcout << "[runner] P10: after clBuildProgram, status=" << status << "\n";
  
  // Rcpp::Rcout << "[runner] P11: before clCreateKernel\n";
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  // Rcpp::Rcout << "[runner] P12: after clCreateKernel, status=" << status << "\n";
  
  // 5) Device Buffers
  // Rcpp::Rcout << "[runner] A: before buffer creation\n";
  
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
  cl_mem bufGrad = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*grad_flat.size(), nullptr, &status);
  
  // Rcpp::Rcout << "[runner] B: after buffer creation\n";
  
  // 6) Set Kernel Args
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufY);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufW);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufGrad);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l1);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l2);
  clSetKernelArg(kernel, arg++, sizeof(int),    &m1);
  
  // 7) Launch
  size_t global = (size_t)m1;
  // Rcpp::Rcout << "[runner] C: before enqueue\n";
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  // Rcpp::Rcout << "[runner] D: after enqueue\n";
  
  // 8) Read back outputs
  // Rcpp::Rcout << "[runner] E: before read qf\n";
  status = clEnqueueReadBuffer(queue, bufQF,   CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(),   qf_flat.data(),
                               0, nullptr, nullptr);
  // Rcpp::Rcout << "[runner] F: after read qf\n";
  
  // Rcpp::Rcout << "[runner] G: before read grad\n";
  status = clEnqueueReadBuffer(queue, bufGrad, CL_TRUE, 0,
                               sizeof(double)*grad_flat.size(), grad_flat.data(),
                               0, nullptr, nullptr);
  // Rcpp::Rcout << "[runner] H: after read grad\n";
  
  // 8a) Sanity-check: error out if both outputs are all zeros
  {
    auto all_zero = [](auto& vec){
      return std::all_of(vec.begin(), vec.end(),
                         [](double x){ return x == 0.0; });
    };
    
    bool qf_is_zero   = all_zero(qf_flat);
    bool grad_is_zero = all_zero(grad_flat);
    
    if (qf_is_zero || grad_is_zero) {
      std::ostringstream msg;
      msg << "OpenCL kernel returned "
          << (qf_is_zero   ? "qf_flat all zeros "   : "")
          << (grad_is_zero ? "grad_flat all zeros." : "");
      throw std::runtime_error(msg.str());
    }
  }
  
  clFlush(queue);
  clFinish(queue);
  
  clReleaseMemObject(bufGrad);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufW);
  clReleaseMemObject(bufY);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufX);
  
  clReleaseKernel       (kernel);
  clReleaseProgram      (program);
  clReleaseCommandQueue (queue);
  clReleaseContext      (context);
}

}
}

#endif

