#include "kernel_loader.h"
#include <CL/cl.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>

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
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
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
  
  std::string nmath_source     = load_kernel_library("nmath");
  std::string test_kernel_code = load_kernel_source("test/nmath_test_kernel.cl");
  
  std::string kernel_code = nmath_source + test_kernel_code;
  std::cout << kernel_code << std::endl;
  
  // 👾 Dispatch minimal kernel
  nmath_test_runner(kernel_code, "nmath_test_kernel", output);
  
  // 📤 Return result to R
  return Rcpp::NumericVector(output.begin(), output.end());
}