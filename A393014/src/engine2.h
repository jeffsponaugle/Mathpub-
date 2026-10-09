// engine2.h - v2 engines.
#pragma once
#include "engine.h"
#include "plan2.h"

RunStats run_cpu2(const Plan2& p2, int threads, std::vector<NumQ>& found);
