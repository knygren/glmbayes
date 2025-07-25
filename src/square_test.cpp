#include <Rcpp.h>
#include <CL/cl.h>
#include <fstream>
#include <sstream>

// Helper to load kernel source from inst/cl
std::string load_kernel_source(const std::string& kernel_name) {
  std::string path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", kernel_name, Rcpp::Named("package") = "glmbayes")  // replace with your package name
  );
  std::ifstream file(path);
  std::ostringstream oss;
  oss << file.rdbuf();
  return oss.str();
}

// [[Rcpp::export]]
Rcpp::NumericVector test_square_kernel(Rcpp::NumericVector inputR) {
  const int n = inputR.size();
  std::vector<float> input(inputR.begin(), inputR.end());
  std::vector<float> output(n);

  cl_platform_id platform;
  cl_device_id device;
  cl_context context;
  cl_command_queue queue;
  cl_program program;
  cl_kernel kernel;

  cl_int status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  
  cl_queue_properties props[] = {0};  // default empty property list
  queue = clCreateCommandQueueWithProperties(context, device, props, &status);  // ✅ preferred
  
  // Load and compile kernel
  std::string kernelSource = load_kernel_source("square.cl");
  const char* src = kernelSource.c_str();
  size_t src_len = kernelSource.size();
  program = clCreateProgramWithSource(context, 1, &src, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  
  kernel = clCreateKernel(program, "square", &status);

  // Set up buffers
  cl_mem input_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float) * n, input.data(), &status);
  cl_mem output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     sizeof(float) * n, nullptr, &status);

  // Set kernel args
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_buf);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buf);
  clSetKernelArg(kernel, 2, sizeof(int), &n);

  // Execute kernel
  size_t global_size = n;
  clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_size, nullptr, 0, nullptr, nullptr);

  // Read result
  clEnqueueReadBuffer(queue, output_buf, CL_TRUE, 0, sizeof(float) * n, output.data(), 0, nullptr, nullptr);

  // Cleanup
  clReleaseMemObject(input_buf);
  clReleaseMemObject(output_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);

  return Rcpp::NumericVector(output.begin(), output.end());
}