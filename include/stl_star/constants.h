#pragma once

namespace stl_star {

// Tolerance used when deciding whether a time offset lies inside a temporal
// window [a, b]: sample j is in the window of sample i iff
//     a - tol <= tau_j - tau_i <= b + tol,
//     tol = max(TIMESTAMP_EPS, TIMESTAMP_ULPS * DBL_EPSILON * max(|tau_0|, |tau_{n-1}|))
// (see Trace).  The absolute part makes e.g. G[0,0.3] on a dt = 0.1 trace
// include tau_3 = 0.30000000000000004; the relative part absorbs the rounding
// of timestamps far from zero (about 1.4e-5 for times around 1e9).
static constexpr double TIMESTAMP_EPS  = 1e-9;
static constexpr double TIMESTAMP_ULPS = 64.0;

// Robustness stopping criterion: the bisection stops when
//   hi - lo <= epsilon * max(|mid|, ROBUSTNESS_RANGE_FLOOR * initial_width).
// The floor keeps the loop finite when the true robustness is (close to) zero.
static constexpr double ROBUSTNESS_RANGE_FLOOR = 0.01;

}  // namespace stl_star
