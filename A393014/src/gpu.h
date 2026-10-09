// gpu.h - Metal engine interface (implemented in gpu_metal.mm).
#pragma once
#include <vector>

#include "engine.h"
#include "plan2.h"

struct GpuEngine;

GpuEngine* gpu_create(bool verbose);  // nullptr if no usable GPU / kernels
const char* gpu_name(GpuEngine* g);
void gpu_destroy(GpuEngine* g);

// Build the table (on the CPU, straight into GPU-visible memory) and search length pl.L.
RunStats run_gpu(GpuEngine* g, const Plan& pl, int threads, std::vector<NumQ>& found, double* tableSeconds,
                 bool verbose);

// v2 (class-partitioned) search on the GPU.
RunStats run_gpu2(GpuEngine* g, const Plan2& p2, int threads, std::vector<NumQ>& found, bool verbose);
