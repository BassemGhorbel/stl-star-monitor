#pragma once

#include "stl_star/formula.h"
#include "stl_star/robustness_monitor_one.h"   // RobustnessResult, compute_robustness_range
#include "stl_star/trace.h"

namespace stl_star {

// ---------------------------------------------------------------------------
// Robustness of a general STL* formula: same bisection as
// robustness_monitor_one (see that header), using
// boolean_monitor_general for every decision step.
// Precondition: negation-free formula.
// ---------------------------------------------------------------------------
RobustnessResult robustness_monitor_general(const Formula& f, const Trace& trace,
                                            double epsilon);

}  // namespace stl_star
