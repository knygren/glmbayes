#include "kernel_loader.h"
#include <CL/cl.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>



//std::string nmath_source = load_kernel_library("nmath");
//std::string arithmetic_source = load_kernel_library("arithmetic");
//std::string test_kernel = load_kernel_source("test/arithmetic_test.cl");

//std::string test_program_source = arithmetic_source + test_kernel;


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





// [[Rcpp::export]]
Rcpp::NumericVector arithmetic_test_wrapper() {
  std::string nmath_source      = load_kernel_library("nmath");
  std::string arithmetic_source = load_kernel_library("arithmetic");
  std::string test_kernel       = load_kernel_source("test/arithmetic_test.cl");
  
  std::string test_program_source = nmath_source + arithmetic_source + test_kernel;
  
  const int n = 6;
  std::vector<float> output(n);
  
  arithmetic_test_kernel_runner(test_program_source, output);
  
  return Rcpp::NumericVector(output.begin(), output.end());
}


void arithmetic_test_v2_kernel_runner(const std::string& test_program_source,const char *kernel_name,
std::vector<float>& output,
float a,
float b) {
  cl_int status;

  
  // 🔍 Select OpenCL platform and device
  
  cl_platform_id platform;
  cl_device_id device;

  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
    
  // 🌐 Create execution context
  
  
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  
  // 🧵 Create command queue for issuing instructions
  
  cl_queue_properties props[] = {0};
  cl_command_queue queue =clCreateCommandQueueWithProperties(context, device, props, &status);
  

  // 📦 Compile OpenCL program from source string
  
  const char* src_ptr = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  // 🎯 Extract kernel from compiled program
  
  
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Allocate device buffer for kernel output
  
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     sizeof(float) * output.size(), nullptr, &status);
  
  // 🗳️ Set kernel arguments
  
    
  clSetKernelArg(kernel, 0, sizeof(float), &a);
  clSetKernelArg(kernel, 1, sizeof(float), &b);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), &output_buf);
  
  // 🚀 Launch kernel with one work-item (can be expanded later for parallelism)
  
  size_t global_size = 1;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0,
                      sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  // 🧹 Cleanup OpenCL resources
  
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
  }



// [[Rcpp::export]]
Rcpp::NumericVector arithmetic_test_v2_wrapper(float a, float b) {
  // 📦 Load core arithmetic functions (addition, multiplication, etc.)
  std::string arithmetic_source = load_kernel_library("arithmetic");
  
  // 🧪 Load the specific test kernel that exercises arithmetic behavior
  std::string test_kernel = load_kernel_source("test/arithmetic_test_v2.cl");
  
  // 🔗 Combine all source components into one complete OpenCL program string
  std::string test_program_source = arithmetic_source + test_kernel;
  
  // 🎯 Set up output buffer
  const int n = 6;
  std::vector<float> output(n);
  
  // 🚀 Run test kernel with provided inputs
  arithmetic_test_v2_kernel_runner(test_program_source, "arithmetic_test_v2",output, a, b);
  
  // 📤 Return results to R as a numeric vector
  return Rcpp::NumericVector(output.begin(), output.end());
}


//////////////////////////////////////////////////////////////////////


void arithmetic_test_parallel_kernel_runner(const std::string& test_program_source,
                                            const char *kernel_name,
                                            std::vector<float>& output,
                                            const std::vector<float>& a_vec,
                                            const std::vector<float>& b_vec) {
  
  
  if (a_vec.size() != b_vec.size()) {
    throw std::runtime_error("Input vectors 'a' and 'b' must have the same size.");
  }
  
  if (output.size() != a_vec.size() * 6) {
    throw std::runtime_error("Output vector size must be 6 times the input vector size.");
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
  const char* src_ptr = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem a_buf = clCreateBuffer(context, CL_MEM_READ_ONLY,sizeof(float) * a_vec.size(), nullptr, &status);
  cl_mem b_buf = clCreateBuffer(context, CL_MEM_READ_ONLY,sizeof(float) * b_vec.size(), nullptr, &status);
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,sizeof(float) * output.size(), nullptr, &status);  
  
  // ⏫ Transfer Data to Device  - Needed  for each input
  clEnqueueWriteBuffer(queue, a_buf, CL_TRUE, 0,sizeof(float) * a_vec.size(), a_vec.data(), 0, nullptr, nullptr);
  clEnqueueWriteBuffer(queue, b_buf, CL_TRUE, 0, sizeof(float) * b_vec.size(), b_vec.data(), 0, nullptr, nullptr);
  
  // 🗳️ Set Kernel Args
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &a_buf);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), &b_buf);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), &output_buf);

    
  // 🚀 Launch Parallel Kernel 
  size_t global_size = a_vec.size(); //  Number of work items...
  
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Output - if more than one output, this will need to have one for each output
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0, sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(a_buf);
  clReleaseMemObject(b_buf);
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}

/////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::NumericMatrix arithmetic_test_parallel_wrapper(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  // 🧮 Ensure input alignment
  if (a.size() != b.size()) {
    Rcpp::stop("Input vectors 'a' and 'b' must have the same length.");
  }
  
  const size_t n = a.size();
  const size_t stride = 6;
  
  // 🎯 Initialize output buffer
  std::vector<float> output(n * stride);
  
  // 📦 Load core and test kernels
  std::string arithmetic_source = load_kernel_library("arithmetic");
  std::string test_kernel = load_kernel_source("test/arithmetic_test_parallel.cl");
  std::string test_program_source = arithmetic_source + test_kernel;
  
  // 🚀 Launch parallel kernel runner
  std::vector<float> a_vec(a.begin(), a.end());
  std::vector<float> b_vec(b.begin(), b.end());
  arithmetic_test_parallel_kernel_runner(test_program_source, "arithmetic_test_parallel", output, a_vec, b_vec);
  
  // 📤 Reshape and return as a matrix
  Rcpp::NumericMatrix result(n, stride);
  std::copy(output.begin(), output.end(), result.begin());
  
  // 🌟 Optional: Add column names
  Rcpp::CharacterVector colnames = {"add", "subtract", "multiply", "divide", "mod", "alpha_form"};
  result.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  
  // 📤 Return output to R
  return result;
  
}



//////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////


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
  
  /// Build Kernel C++ program (also print for now)
  
  std::string nmath_source      = load_kernel_library("nmath");
  std::string arithmetic_source = load_kernel_library("arithmetic");
  std::string test_kernel       = load_kernel_source("test/arithmetic_test.cl");
  
  
//  std::cout << nmath_source << std::endl;

    std::string test_program_source = nmath_source + arithmetic_source + test_kernel;
  
  std::cout << test_program_source << std::endl;
  
  
  //// Kernel Call
  
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



void dpq_macro_kernel_runner(const std::string& test_program_source,
                                  const char *kernel_name,
                                  std::vector<float>& output,
                                  float p_in,
                                  bool logp_in,
                                  bool lt_in) {
  if (output.empty()) {
    throw std::runtime_error("Output vector must be preallocated and non-empty.");
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
  const char* src_ptr = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * output.size(), nullptr, &status);
  
  // 🗳️ Set Kernel Args
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  clSetKernelArg(kernel, 1, sizeof(float), &p_in);
  clSetKernelArg(kernel, 2, sizeof(cl_bool), &logp_in);
  clSetKernelArg(kernel, 3, sizeof(cl_bool), &lt_in);
  
  // 🚀 Launch Parallel Kernel
//  size_t global_size = output.size();  // Each work item writes one output entry
  
  size_t global_size = 4;
  output.resize(4);  // Only one thread writing 4 floats
//  size_t global_size = output.size();  // Each work item writes one output entry
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Output
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0, sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}




// [[Rcpp::export]]
Rcpp::NumericVector dpq_macro_usage_wrapper(double p, bool log_p, bool lower_tail) {
  const size_t stride = 4;  // Number of macro probes we'll test
  std::vector<float> output(stride);
  
  // 🚀 Load test kernel and required macro layers
//    std::string nmath_source      = load_kernel_library("nmath");
//  std::string dpq_macros = load_kernel_source("dpq");        // Contains the full macro set
  //std::string dpq_prelude = load_kernel_source("nmath/dpq_prelude.cl"); 
  std::string test_kernel = load_kernel_source("test/dpq_macro_usage.cl");
  
//  std::string kernel_code = nmath_source + test_kernel;
//  std::string kernel_code = dpq_prelude + test_kernel;
  std::string kernel_code = test_kernel;
  
  std::cout << kernel_code << std::endl;
  
  
  // 👾 Dispatch kernel
  dpq_macro_kernel_runner(kernel_code, "dpq_macro_usage", output, p, log_p, lower_tail);
  
  // 📤 Return result to R for inspection
  return Rcpp::NumericVector(output.begin(), output.end());
}





void basic_kernel_runner(const std::string& test_program_source,
                         const char* kernel_name,
                         std::vector<float>& output) {
  if (output.size() < 4) {
    output.resize(4);  // Ensure space for 4 results
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
  const char* src_ptr = test_program_source.c_str();
  size_t src_len = test_program_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  // 🐛 Optional Build Log Dump
  char build_log[4096];
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(build_log), build_log, nullptr);
  std::cerr << "Build Log:\n" << build_log << std::endl;
  
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * output.size(), nullptr, &status);
  
  // 🗳️ Set Kernel Arg
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buf);
  
  // 🚀 Launch Kernel with One Work-Item
  size_t global_size = 1;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Results
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0, sizeof(float) * output.size(), output.data(), 0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}



// [[Rcpp::export]]
Rcpp::NumericVector basic_kernel_wrapper() {
  const size_t stride = 5;  // Number of values we expect back
  std::vector<float> output(stride);
  
  //std::string dpq_prelude = load_kernel_source("nmath/dpq_prelude.cl"); 
 //     std::string dpq_source      = load_kernel_library("dpq");
  
  
    // 🚀 Load standalone test kernel
  std::string test_kernel = load_kernel_source("test/test_dpq_kernel.cl");  // Should contain test_kernel()
 
//   std::string kernel_code = dpq_source + test_kernel;
 std::string kernel_code = test_kernel;
  std::cout << kernel_code << std::endl;
 
  // 👾 Dispatch minimal kernel
  basic_kernel_runner(kernel_code, "test_dpq_kernel", output);
  
  // 📤 Return result to R
  return Rcpp::NumericVector(output.begin(), output.end());
}


void basic_kernel_runner_v2(const std::string& source,
                         const char* kernel_name,
                         std::vector<float>& output) {
  if (output.size() != 25) {
    throw std::runtime_error("Output vector must be preallocated with size 4.");
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
Rcpp::NumericVector basic_kernel_wrapper_v2() {
  const size_t stride = 25;  // Number of values we expect back
  std::vector<float> output(stride);
  
  //std::string dpq_prelude = load_kernel_source("nmath/dpq_prelude.cl"); 
       std::string dpq_source      = load_kernel_library("dpq");
  
  
  // 🚀 Load standalone test kernel
  std::string test_kernel = load_kernel_source("test/basic_kernel.cl");  // Should contain test_kernel()
  
     std::string kernel_code = dpq_source + test_kernel;
 // std::string kernel_code = test_kernel;
  std::cout << kernel_code << std::endl;
  
  // 👾 Dispatch minimal kernel
  basic_kernel_runner_v2(kernel_code, "basic_kernel", output);
  
  // 📤 Return result to R
  return Rcpp::NumericVector(output.begin(), output.end());
}
