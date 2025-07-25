#include "kernel_loader.h"
#include <CL/cl.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>



//std::string nmath_source = load_kernel_library("nmath");
//std::string arithmetic_source = load_kernel_library("arithmetic");
//std::string test_kernel = load_kernel_source("test/arithmetic_test.cl");

//std::string test_program_source = arithmetic_source + test_kernel;


// [[Rcpp::export]]
Rcpp::NumericVector run_test_kernel() {
  const int n = 6;  // Expecting 5 arithmetic outputs
  std::vector<float> output(n);
  
  cl_platform_id platform;
  cl_device_id device;
  cl_context context;
  cl_command_queue queue;
  cl_program program;
  cl_kernel kernel;
  cl_int status;
  
  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  
  cl_queue_properties props[] = {0};
  queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  std::string nmath_source      = load_kernel_library("nmath");
  std::string arithmetic_source = load_kernel_library("arithmetic");
  std::string test_kernel       = load_kernel_source("test/arithmetic_test.cl");
  
  std::string test_program_source = nmath_source + arithmetic_source + test_kernel;
  
  std::cout << test_program_source << std::endl;
  
  const char* src = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  program = clCreateProgramWithSource(context, 1, &src, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  kernel = clCreateKernel(program, "arithmetic_test", &status);
  
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     sizeof(float) * n, nullptr, &status);
  
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  
  size_t global_size = 1;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0,
                      sizeof(float) * n, output.data(), 0, nullptr, nullptr);
  
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
  
  return Rcpp::NumericVector(output.begin(), output.end());
}



