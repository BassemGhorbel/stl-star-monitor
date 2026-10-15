#pragma once

// ---------------------------------------------------------------------------
// stl_star_monitor.h -- one-call API for STL* monitoring.
//
// Every function takes
//   formula     STL* formula as text (see formula_parser.h for the syntax)
//   signal      signal[k][i] = value of signal k at sample i, one signal per
//               letter of the formula, in alphabetical order of the letters:
//               for a formula over x and z, signal[0] is x and signal[1] is z.
//   dt          sampling period of a uniform trace (tau_i = i * dt), or
//   timestamps  the sample times of a non-uniform trace.
//
// Example:
//   auto r = monitor_general("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)",
//                            {{0.0, 1.0, 2.0, 3.0, 4.0},     // signal a
//                             {0.5, 0.5, 0.5, 0.5, 0.5}},    // signal b
//                            1.0);                           // dt
//   r.holds_at_start   // does the formula hold at time tau_0?
//   r.intervals        // all times at which it holds
//
// Which function?
//   monitor_one / robustness_one           STL*_1 formulas: every constraint
//                                          reads only the variable of its
//                                          nearest freeze operator, e.g.
//                                          "G a*. F[0,3] (a - a* > 1)".
//   monitor_general / robustness_general   any STL* formula, including nested
//                                          freeze operators whose constraints
//                                          read outer frozen values.
//
// All functions throw std::runtime_error on syntax errors and
// std::invalid_argument on invalid traces or unsupported formulas.
// ---------------------------------------------------------------------------

#include "stl_star/boolean_monitor_general.h"
#include "stl_star/boolean_monitor_one.h"
#include "stl_star/formula_parser.h"
#include "stl_star/interval_operations.h"
#include "stl_star/robustness_monitor_general.h"
#include "stl_star/robustness_monitor_one.h"
#include "stl_star/trace.h"

#include <string>
#include <vector>

namespace stl_star {

using Signal = std::vector<std::vector<double>>;   // signal[k][i]

struct MonitorResult {
    bool      holds_at_start;  // does the formula hold at the first sample?
    Intervals intervals;       // maximal closed time intervals where it holds
};

namespace detail {
inline MonitorResult to_result(Intervals ivs, const Trace& trace) {
    const double t0 = trace.timestamp(0);
    const bool start = !ivs.empty() && ivs.front().lo == t0;   // endpoints are exact timestamps
    return {start, std::move(ivs)};
}
}  // namespace detail

// ---- Boolean monitoring ----------------------------------------------------
inline MonitorResult monitor_one(const std::string& formula, const Signal& signal,
                                 double dt) {
    const Trace trace = Trace::make_uniform(dt, signal);
    return detail::to_result(boolean_monitor_one(parse_formula(formula), trace), trace);
}

inline MonitorResult monitor_one(const std::string& formula, const Signal& signal,
                                 const std::vector<double>& timestamps) {
    const Trace trace = Trace::make_non_uniform(timestamps, signal);
    return detail::to_result(boolean_monitor_one(parse_formula(formula), trace), trace);
}

inline MonitorResult monitor_general(const std::string& formula, const Signal& signal,
                                     double dt) {
    const Trace trace = Trace::make_uniform(dt, signal);
    return detail::to_result(boolean_monitor_general(parse_formula(formula), trace), trace);
}

inline MonitorResult monitor_general(const std::string& formula, const Signal& signal,
                                     const std::vector<double>& timestamps) {
    const Trace trace = Trace::make_non_uniform(timestamps, signal);
    return detail::to_result(boolean_monitor_general(parse_formula(formula), trace), trace);
}

// ---- Robustness (negation-free formulas) -----------------------------------
// Returns {lo, hi, monitor_calls} with the robustness at tau_0 in [lo, hi];
// epsilon is the relative error (0.01 = 1%).  lo = hi = +-inf when the
// robustness is infinite (e.g. F[5,10] on a trace shorter than 5 time units).
inline RobustnessResult robustness_one(const std::string& formula, const Signal& signal,
                                       double dt, double epsilon = 0.01) {
    return robustness_monitor_one(parse_formula(formula),
                                  Trace::make_uniform(dt, signal), epsilon);
}

inline RobustnessResult robustness_one(const std::string& formula, const Signal& signal,
                                       const std::vector<double>& timestamps,
                                       double epsilon = 0.01) {
    return robustness_monitor_one(parse_formula(formula),
                                  Trace::make_non_uniform(timestamps, signal), epsilon);
}

inline RobustnessResult robustness_general(const std::string& formula, const Signal& signal,
                                           double dt, double epsilon = 0.01) {
    return robustness_monitor_general(parse_formula(formula),
                                      Trace::make_uniform(dt, signal), epsilon);
}

inline RobustnessResult robustness_general(const std::string& formula, const Signal& signal,
                                           const std::vector<double>& timestamps,
                                           double epsilon = 0.01) {
    return robustness_monitor_general(parse_formula(formula),
                                      Trace::make_non_uniform(timestamps, signal), epsilon);
}

}  // namespace stl_star
