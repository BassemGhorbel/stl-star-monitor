#include "stl_star/robustness_monitor_one.h"
#include "stl_star/boolean_monitor_one.h"
#include "monitor_common.h"
#include "robustness_common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace stl_star {

// ---------------------------------------------------------------------------
// Bounds on all finite atomic robustness values.
//
// An atom  h(s) - g(x*) ~ t  has robustness  d - t  (>, >=)  or  t - d  (<, <=)
// with d = h - g, h ranging over the samples and g over the reachable frozen
// values (a signal predicate is the case g = 0).  Every reachable vector of
// frozen values takes each variable from a trace value of its component, so
// enumerating those combinations bounds g for ANY expression g (no
// monotonicity assumption).  Floating-point subtraction is monotone, so
// h_min - g_max and h_max - g_min bound every computed difference h - g exactly.
// Bounds that overflow are clamped to the largest finite double.
// ---------------------------------------------------------------------------
std::pair<double, double> compute_robustness_range(const Formula& f, const Trace& trace) {
    monitor_detail::check_input(f, trace, /*require_fragment_one=*/false);
    const double inf = std::numeric_limits<double>::infinity();
    const double big = std::numeric_limits<double>::max();
    double lo = inf, hi = -inf;
    auto take = [&](double v) {
        if (std::isnan(v)) return;
        v = std::max(-big, std::min(big, v));
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    };

    // Trace values of the component frozen by each variable.  Duplicates are
    // removed (O(n log n)) only for constraints over several variables, where
    // they would multiply the number of combinations; a single-variable
    // constraint just evaluates g at all n values, keeping the range O(n).
    std::vector<const std::vector<double>*> raw_values(f.num_freeze_vars, nullptr);
    std::vector<std::vector<double>> distinct_values(f.num_freeze_vars);
    for (const FormulaNode& nd : f.nodes)
        if (nd.type == FormulaType::FREEZE) raw_values[nd.freeze_var] = &trace.signal[nd.freeze_component];
    auto values_of = [&](int var, bool distinct) -> const std::vector<double>& {
        if (!distinct) return *raw_values[var];
        std::vector<double>& d = distinct_values[var];
        if (d.empty()) {
            d = *raw_values[var];
            std::sort(d.begin(), d.end());
            d.erase(std::unique(d.begin(), d.end()), d.end());
        }
        return d;
    };

    std::vector<double> s(trace.num_components), frozen(f.num_freeze_vars, 0.0);
    for (const FormulaNode& nd : f.nodes) {
        if (nd.type != FormulaType::SIGNAL_PRED && nd.type != FormulaType::SIGNAL_CONSTRAINT)
            continue;
        const bool greater = nd.comp_op == CompOp::GT || nd.comp_op == CompOp::GE;

        double h_min = inf, h_max = -inf;
        for (int i = 0; i < trace.length; ++i) {
            trace.sample_at(i, s);
            const double h = nd.signal_fn(s);
            if (std::isfinite(h)) { h_min = std::min(h_min, h); h_max = std::max(h_max, h); }
        }
        if (h_min > h_max) continue;                       // no finite value at all

        double g_min = 0.0, g_max = 0.0;                   // signal predicate: g = 0
        if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
            g_min = inf; g_max = -inf;
            const std::vector<int>& vars = nd.constraint_vars;
            const bool distinct = vars.size() > 1;
            std::vector<const std::vector<double>*> vals;
            for (int v : vars) vals.push_back(&values_of(v, distinct));
            std::vector<size_t> pick(vars.size(), 0);    // odometer over value combinations
            for (;;) {
                for (size_t v = 0; v < vars.size(); ++v) frozen[vars[v]] = (*vals[v])[pick[v]];
                const double g = nd.frozen_fn(frozen);
                if (std::isfinite(g)) { g_min = std::min(g_min, g); g_max = std::max(g_max, g); }
                size_t v = 0;
                while (v < vars.size() && ++pick[v] == vals[v]->size()) pick[v++] = 0;
                if (v == vars.size()) break;
            }
            if (g_min > g_max) continue;
        }
        const double d_min = h_min - g_max, d_max = h_max - g_min, t = nd.threshold;
        if (greater) { take(d_min - t); take(d_max - t); }
        else         { take(t - d_max); take(t - d_min); }
    }
    if (lo > hi) { lo = 0.0; hi = 0.0; }               // no finite atomic value
    return {lo, hi};
}

RobustnessResult robustness_monitor_one(const Formula& f, const Trace& trace,
                                        double epsilon) {
    return robustness_detail::bisect(f, trace, epsilon, boolean_monitor_one_holds_at_start);
}

}  // namespace stl_star
