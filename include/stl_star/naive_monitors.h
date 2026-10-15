#pragma once

#include "stl_star/formula.h"
#include "stl_star/interval_operations.h"
#include "stl_star/trace.h"

namespace stl_star {

// ---------------------------------------------------------------------------
// Reference (naive) monitors.
//
// They evaluate the STL* semantics point by point, directly from the
// definitions, and are used to cross-check the fast monitors:
//   * the Boolean reference monitors return the same results as
//     boolean_monitor_one / boolean_monitor_general;
//   * the robustness reference monitors return the exact robustness, which
//     lies in the bracket [lo, hi] of robustness_monitor_one /
//     robustness_monitor_general.
// Intended for testing and for small traces.
//
// How they work.  A FREEZE node binds its variable to s_x(i) for every
// sample i and re-evaluates its scope on all n samples.  The STL*_1 monitors
// compute a nested FREEZE once (its scope does not read the outer variable);
// the general monitors re-evaluate it at every step of the enclosing FREEZE.
// Timed F and G use prefix counts (Boolean) or Lemire's sliding max / min
// filter (robustness); until scans the window of every sample.
//
// The STL*_1 monitors throw std::invalid_argument for formulas outside
// STL*_1, the robustness monitors for formulas with negation.
// ---------------------------------------------------------------------------

// ---- Boolean ----------------------------------------------------------------
Intervals boolean_monitor_one_naive(const Formula& f, const Trace& trace);
bool      boolean_monitor_one_naive_holds_at_start(const Formula& f, const Trace& trace);

Intervals boolean_monitor_general_naive(const Formula& f, const Trace& trace);
bool      boolean_monitor_general_naive_holds_at_start(const Formula& f, const Trace& trace);

// ---- Robustness (exact value at the first sample; +-inf for empty windows) --
double robustness_monitor_one_naive(const Formula& f, const Trace& trace);
double robustness_monitor_general_naive(const Formula& f, const Trace& trace);

}  // namespace stl_star
