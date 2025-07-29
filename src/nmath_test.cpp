#include "kernel_loader.h"
#include <CL/cl.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>
#include "configure_OpenCL.h"

// 🚀 Runner for nmath test kernel
void nmath_test_runner(const std::string& source,
                       const char* kernel_name,
                       std::vector<float>& output) {
  if (output.size() != 9) {
    throw std::runtime_error("Output vector must be preallocated with size 9.");
  }
  
  cl_int status;
  
  // 🔍 Platform & Device
  cl_platform_id platform;
  cl_device_id device;
  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 🌐 Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);

  
  
    cl_queue_properties props[] = {0};
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  // 📦 Program & Kernel
  const char* src_ptr = source.c_str();
  size_t src_len = source.size();

  // Get Configuration (Based on context and device)
  OpenCLConfig cfg = configureOpenCL(context, device);
  const char *opts = cfg.buildOptions.empty()
    ? nullptr
    : cfg.buildOptions.c_str();
  
  
  
  
    cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 
                           0,    // build for all devices
                           nullptr, //device list
                           opts, // Configuration options
//                           nullptr, // Configuration options
                           nullptr,
                           nullptr);
  
  
  // 📣 Retrieve build log
  size_t log_size;
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
  char *log = (char *)malloc(log_size);
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
  printf("Build Log:\n%s\n", log);
  free(log);
  
  // 📦  Kernel
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     sizeof(float) * output.size(), nullptr, &status);
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  
  // 🚀 Launch Kernel (Single Work-Item)
  size_t global_size = 1;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Output
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0,
                      sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}

// [[Rcpp::export]]
Rcpp::NumericVector nmath_test_wrapper() {
  const size_t stride = 9;  // Number of values expected from kernel
  std::vector<float> output(stride);
  
  std::string OPENCL_source     = load_kernel_source("OPENCL.CL");
  std::string rmath_source     = load_kernel_library("rmath");
  std::string nmath_source     = load_kernel_library("nmath");
  std::string dpq_source     = load_kernel_library("dpq");
  std::string nmath_source2     = load_kernel_source("nmath/nmath.cl");
  
  
  std::string chebyshev_source     = load_kernel_source("nmath/chebyshev.cl");
  std::string d1mach_source     = load_kernel_source("nmath/d1mach.cl");
  std::string dnorm_source     = load_kernel_source("nmath/dnorm.cl");
  std::string fmax2_source     = load_kernel_source("nmath/fmax2.cl");
  std::string gammalims_source     = load_kernel_source("nmath/gammalims.cl");
  std::string lgammacor_source     = load_kernel_source("nmath/lgammacor.cl");
  std::string log1p_source     = load_kernel_source("nmath/log1p.cl");
  std::string pnorm_source     = load_kernel_source("nmath/pnorm.cl");
  std::string stirlerr_large_source     = load_kernel_source("nmath/stirlerr_large.cl");
  std::string expm1_source     = load_kernel_source("nmath/expm1.cl");
  std::string gamma_source     = load_kernel_source("nmath/gamma.cl");
  std::string lgamma_source     = load_kernel_source("nmath/lgamma.cl");
  std::string lgamma1p_source     = load_kernel_source("nmath/lgamma1p.cl");
  std::string stirlerr_small_source     = load_kernel_source("nmath/stirlerr_small.cl");
  std::string dgamma_source     = load_kernel_source("nmath/dgamma.cl");
  std::string stirlerr_source     = load_kernel_source("nmath/stirlerr.cl");
  std::string bd0_source     = load_kernel_source("nmath/bd0.cl");
  std::string dbinom_source     = load_kernel_source("nmath/dbinom.cl");
  std::string dpois_source     = load_kernel_source("nmath/dpois.cl");
  std::string pgamma_source     = load_kernel_source("nmath/pgamma.cl");
  
  std::string test_kernel_code = load_kernel_source("test/nmath_test_kernel.cl");
  
//  std::string kernel_code = rmath_source+ nmath_source + test_kernel_code;
  
  std::string kernel_code =
    OPENCL_source +
    rmath_source + 
    "\n" + dpq_source
  +  "\n" +nmath_source2 
  + "\n" + chebyshev_source
  + "\n" + d1mach_source
  + "\n" + dnorm_source
  + "\n" + fmax2_source
  + "\n" + gammalims_source
  + "\n" + lgammacor_source
  + "\n" + log1p_source
  + "\n" + pnorm_source
  + "\n" + stirlerr_large_source
  + "\n" + expm1_source
  + "\n" + gamma_source
  + "\n" + lgamma_source
  + "\n" + lgamma1p_source
  + "\n" + stirlerr_small_source
//  + "\n" + dgamma_source
  + "\n" + stirlerr_source
//  + "\n" + bd0_source
//  + "\n" + dbinom_source
//  + "\n" + dpois_source
//  + "\n" + pgamma_source
  + "\n" + test_kernel_code;
  std::cout << kernel_code << std::endl;
  
  
  // 👾 Dispatch minimal kernel
  nmath_test_runner(kernel_code, "nmath_test_kernel", output);
  
  // 📤 Return result to R
  return Rcpp::NumericVector(output.begin(), output.end());
}