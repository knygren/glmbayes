#include <vector>

#ifdef USE_OPENCL
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#endif

#include "openclPort.h"

using namespace openclPort;

#ifdef USE_OPENCL

namespace openclPort{
int detect_num_gpus_internal() {
  cl_uint num_platforms = 0;
  cl_int status = clGetPlatformIDs(0, nullptr, &num_platforms);
  if (status != CL_SUCCESS || num_platforms == 0) return 0;
  
  std::vector<cl_platform_id> platforms(num_platforms);
  clGetPlatformIDs(num_platforms, platforms.data(), nullptr);
  
  int total_compute_units = 0;
  for (auto& platform : platforms) {
    cl_uint num_devices = 0;
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &num_devices);
    if (status != CL_SUCCESS || num_devices == 0) continue;
    
    std::vector<cl_device_id> devices(num_devices);
    clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, num_devices, devices.data(), nullptr);
    
    for (auto& device : devices) {
      cl_uint units = 0;
      status = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, nullptr);
      if (status == CL_SUCCESS) total_compute_units += units;
    }
  }
  return total_compute_units;
}

}
#endif