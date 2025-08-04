#ifdef USE_OPENCL
#include "kernel_loader.h"
#include <CL/cl.h>
#endif
#include <iostream>
#include <vector>
#include <Rcpp.h>





//std::string nmath_source = load_kernel_library("nmath");
//std::string arithmetic_source = load_kernel_library("arithmetic");
//std::string test_kernel = load_kernel_source("test/arithmetic_test.cl");

//std::string test_program_source = arithmetic_source + test_kernel;

#ifdef USE_OPENCL

void arithmetic_test_kernel_runner(const std::string& test_program_source,
                                   std::vector<float>& output) {
  cl_int status;
  
  cl_platform_id platform;
  cl_device_id device;
  cl_context context;
  cl_command_queue queue;
  
  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  
  cl_queue_properties props[] = {0};
  queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  const char* src_ptr = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  cl_kernel kernel = clCreateKernel(program, "arithmetic_test", &status);
  
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     sizeof(float) * output.size(), nullptr, &status);
  
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  
  size_t global_size = 1;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0,
                      sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}
#endif




// [[Rcpp::export]]
Rcpp::NumericVector arithmetic_test_wrapper() {

  const int n = 6;
  std::vector<float> output(n);

  #ifdef USE_OPENCL
  
    std::string nmath_source      = load_kernel_library("nmath");
  std::string arithmetic_source = load_kernel_library("arithmetic");
  std::string test_kernel       = load_kernel_source("test/arithmetic_test.cl");
  
  std::string test_program_source = nmath_source + arithmetic_source + test_kernel;
  
  
  arithmetic_test_kernel_runner(test_program_source, output);
  #else
  Rcpp::Rcout << "[INFO] OpenCL not available — returning zero vector.\n";
     #endif
  
  return Rcpp::NumericVector(output.begin(), output.end());
}

