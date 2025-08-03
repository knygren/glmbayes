// configure_OpenCL.h
#ifdef USE_OPENCL
#ifndef CONFIGURE_OPENCL_H
#define CONFIGURE_OPENCL_H

#include <CL/cl.h>
#include <string>

struct OpenCLConfig {
  bool have_expm1;
  bool have_log1p;
  std::string buildOptions;
};

// Probes the given context/device for built-in functions.
// Returns a filled OpenCLConfig struct.
OpenCLConfig configureOpenCL(cl_context context, cl_device_id device);

#endif // CONFIGURE_OPENCL_H
#endif // USE_OPENCL
