#pragma once

#include "stl_star/formula.h"
#include "stl_star/trace.h"

#include <utility>

namespace stl_star {

// ---------------------------------------------------------------------------
// Quantitative robustness monitoring by bisection.
//
// Robustness rho(phi) at the first sample (negation-free formulas):
//   rho(h ~ g)      = h(s) - g(x*)  for > and >=,   g(x*) - h(s)  for < and <=
//                     (a signal predicate is the case g = 0)
//   rho(phi && psi) = min,  rho(phi || psi) = max
//   rho(F[a,b] phi) = max over the window (-inf if the window is empty)
//   rho(G[a,b] phi) = min over the window (+inf if the window is empty)
//   rho(phi U psi)  = max over j in the window of min(rho(psi, j), min of rho(phi, k)
//                     over the samples k that phi must cover)
//   rho(x*. phi)    = rho of phi with x* bound to the current value.
// rho > 0 implies the formula holds at tau_0; rho < 0 implies it does not.
//
// Decision step (robustness theorem): rho(phi) >= r  iff  phi'(r) holds, where
// phi'(r) replaces every atom "h > g" / "h >= g" by "h >= g + r" and every
// atom "h < g" / "h <= g" by "h <= g - r".  Each bisection step is one call
// of the Boolean monitor.
//
// Algorithm:
//   1. [lo, hi] = compute_robustness_range: bounds on every finite atomic
//      robustness value (O(n) for STL*_1, see below).
//   2. If phi'(lo) is false then rho = -inf; if phi'(+inf) holds then
//      rho = +inf (both are returned as lo = hi = +-inf).
//   3. Otherwise rho is in [lo, hi]; bisect until
//        hi - lo <= epsilon * max(|mid|, ROBUSTNESS_RANGE_FLOOR * initial width),
//      i.e. epsilon is the RELATIVE error with respect to the robustness value
//      itself (with an absolute floor for values near zero).
//
// Precondition: the formula has no negation (throws std::invalid_argument).
// ---------------------------------------------------------------------------
struct RobustnessResult {
    double lo;             // rho is in [lo, hi]
    double hi;
    int    monitor_calls;  // number of Boolean monitor calls made
};

// Bounds [lo, hi] containing every finite atomic robustness value over the
// trace, hence rho whenever rho is finite.  For a signal constraint every
// combination of trace values of the freeze variables it reads is tried, so
// the bounds are valid for any constraint expression: O(n) for constraints
// over one freeze variable, O(n^k) for constraints that read k variables.
std::pair<double, double> compute_robustness_range(const Formula& f, const Trace& trace);

// Robustness of an STL*_1 formula (uses boolean_monitor_one).
// epsilon in (0, 1]: relative error, e.g. 0.01 = 1%.
RobustnessResult robustness_monitor_one(const Formula& f, const Trace& trace,
                                        double epsilon);

}  // namespace stl_star
