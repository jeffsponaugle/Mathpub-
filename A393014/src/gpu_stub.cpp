// gpu_stub.cpp - used when building without Metal.
#include "gpu.h"

GpuEngine* gpu_create(bool) { return nullptr; }
const char* gpu_name(GpuEngine*) { return "none"; }
void gpu_destroy(GpuEngine*) {}
RunStats run_gpu(GpuEngine*, const Plan&, int, std::vector<NumQ>&, double*, bool) { return RunStats{}; }
RunStats run_gpu2(GpuEngine*, const Plan2&, int, std::vector<NumQ>&, bool) { return RunStats{}; }
