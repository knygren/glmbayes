// rmath_test.cpp

#include "kernel_loader.h"
#include <CL/cl.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>

// 🚀 Runner for rmath test kernel
void rmath_test_runner(const std::string& source,
                       const char* kernel_name,
                       std::vector<double>& output) {
  if (output.size() != 30) {
    throw std::runtime_error("Output vector must be preallocated with size 30.");
  }
  
  cl_int status;
  
  // 🔍 Platform & Device
  cl_platform_id platform;
  cl_device_id device;
  status  = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 🌐 Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  cl_queue_properties props[] = { 0 };
  cl_command_queue queue = clCreateCommandQueueWithProperties(
    context, device, props, &status
  );
  
  // 📦 Program & Kernel
  const char* src_ptr = source.c_str();
  size_t src_len     = source.size();
  cl_program program = clCreateProgramWithSource(
    context, 1, &src_ptr, &src_len, &status
  );
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  // 📣 Retrieve build log
  size_t log_size;
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
  char *log = (char *)malloc(log_size);
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
  printf("Build Log:\n%s\n", log);
  free(log);
  
  // 📦 Kernel
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem output_buf = clCreateBuffer(
    context,
    CL_MEM_WRITE_ONLY,
    sizeof(double) * output.size(),
    nullptr,
    &status
  );
  status |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  
  // 🚀 Launch Kernel (Single Work-Item)
  size_t global_size = 1;
  status = clEnqueueNDRangeKernel(
    queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr
  );
  
  // 📥 Retrieve Output
  status = clEnqueueReadBuffer(
    queue, output_buf, CL_TRUE,
    0,
    sizeof(double) * output.size(),
    output.data(),
    0, nullptr, nullptr
  );
  
  // 🧹 Cleanup
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}

// [[Rcpp::export]]
Rcpp::NumericVector rmath_test_wrapper() {
  const size_t stride = 30;  // Number of values expected from kernel
  std::vector<double> output(stride);
  
  // Load rmath core and test kernel
  std::string OPENCL_source     = load_kernel_source("OPENCL.CL");
  std::string rmath_source      = load_kernel_library("rmath");
  std::string test_kernel_code  = load_kernel_source("test/rmath_test_kernel.cl");
  std::string kernel_code       = OPENCL_source + rmath_source + test_kernel_code;
  
  // Optional: dump the combined kernel to stdout for debugging
  std::cout << kernel_code << std::endl;
  
  // Dispatch minimal kernel
  rmath_test_runner(kernel_code, "rmath_test_kernel", output);
  
  // Return results back to R
  return Rcpp::NumericVector(output.begin(), output.end());
}