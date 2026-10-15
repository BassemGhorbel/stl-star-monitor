#pragma once
// Bisection shared by the STL*_1 and general robustness monitors (internal).

#include "stl_star/robustness_monitor_one.h"
#include "stl_star/constants.h"
#include "monitor_common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace stl_star {

namespace robustness_detail {

// Run the bisection of robustness_monitor_one.h with `holds_at_start`
// (one of the Boolean monitors) as the decision procedure.
template <class HoldsAtStart>
RobustnessResult bisect(const Formula& f, const Trace& trace, double epsilon,
                        HoldsAtStart holds_at_start) {
    if (!(epsilon > 0.0 && epsilon <= 1.0))
        throw std::invalid_argument("epsilon must be in (0, 1]");
    f.validate();   // the trace is checked by compute_robustness_range below
    monitor_detail::check_negation_free(f);

    // phi'(r): every atom compares non-strictly against its shifted threshold,
    // so that phi'(r) holds exactly when rho >= r.
    Formula shifted = f;
    for (FormulaNode& nd : shifted.nodes) {
        if (nd.type != FormulaType::SIGNAL_PRED && nd.type != FormulaType::SIGNAL_CONSTRAINT)
            continue;
        if (nd.comp_op == CompOp::GT) nd.comp_op = CompOp::GE;
        if (nd.comp_op == CompOp::LT) nd.comp_op = CompOp::LE;
    }
    int calls = 0;
    auto holds = [&](double r) {
        for (size_t k = 0; k < shifted.nodes.size(); ++k) {
            FormulaNode& nd = shifted.nodes[k];
            if (nd.type != FormulaType::SIGNAL_PRED && nd.type != FormulaType::SIGNAL_CONSTRAINT)
                continue;
            const double base = f.nodes[k].threshold;
            nd.threshold = (nd.comp_op == CompOp::GE) ? base + r : base - r;
        }
        ++calls;
        return holds_at_start(shifted, trace);
    };

    // With threshold 0 an atom of phi'(r) is "d >= r" (or "d <= -r") on the
    // very value d = h - g whose robustness is d (or -d), so each decision is
    // exact.  A small margin keeps the bracket valid when an atom carries a
    // non-zero base threshold (then d - base is rounded).
    const double inf = std::numeric_limits<double>::infinity();
    const double big = std::numeric_limits<double>::max();
    auto [lo, hi] = compute_robustness_range(f, trace);
    const double margin = 1e-9 * std::max({1.0, std::abs(lo), std::abs(hi)});
    lo = std::max(lo - margin, -big);
    hi = std::min(hi + margin, big);
    if (!holds(lo)) return {-inf, -inf, calls};   // rho < every finite atom value
    if (holds(inf)) return {inf, inf, calls};     // rho = +inf
    // Now rho is finite, hence an atomic value: lo <= rho <= hi.  Midpoints
    // and widths are computed from halves so that nothing overflows.
    auto half_width = [&] { return 0.5 * hi - 0.5 * lo; };
    const double floor_width = 2.0 * ROBUSTNESS_RANGE_FLOOR * half_width();
    while (hi > lo) {
        const double mid = 0.5 * lo + 0.5 * hi;
        if (2.0 * half_width() <= epsilon * std::max(std::abs(mid), floor_width)) break;
        if (mid <= lo || mid >= hi) break;         // no representable midpoint left
        if (holds(mid)) lo = mid; else hi = mid;
    }
    return {lo, hi, calls};
}

}  // namespace robustness_detail

}  // namespace stl_star
