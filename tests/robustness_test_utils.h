#pragma once
// Robustness checks shared by the two robustness test files.
//
// The robustness monitors are approximate: they return a bracket [lo, hi]
// that must contain the exact robustness rho and be narrower than the
// epsilon stopping rule.  The exact rho comes from the reference monitors;
// the reference monitors are themselves checked against the independent
// oracle on random formulas.
#include "random_generators.h"

namespace test_util {

using RobustnessFn = RobustnessResult (*)(const Formula&, const Trace&, double);
using ReferenceFn  = double (*)(const Formula&, const Trace&);

// Two exact robustness values agree (same infinity, or equal up to rounding
// of hand-written decimal expectations).
inline bool same_robustness(double x, double y) {
    if (std::isinf(x) || std::isinf(y)) return x == y;
    return std::abs(x - y) <= 1e-9 * std::max(1.0, std::abs(y));
}

// Is `rho` (exact value) inside the returned bracket, and does the bracket meet
// the stopping criterion  hi - lo <= eps * max(|mid|, 1% of the initial width)?
inline bool robustness_ok(const RobustnessResult& r, double rho, double eps,
                          const Formula& f, const Trace& tr, std::string& why) {
    if (std::isinf(rho)) {
        if (r.lo == rho && r.hi == rho) return true;
        why = "expected " + num(rho);
        return false;
    }
    const double tol = 1e-9 * std::max(1.0, std::abs(rho));
    if (!(r.lo - tol <= rho && rho <= r.hi + tol)) { why = "[OUT OF RANGE]"; return false; }
    auto [lo0, hi0] = compute_robustness_range(f, tr);
    const double margin = 1e-9 * std::max({1.0, std::abs(lo0), std::abs(hi0)});
    const double floor_width = ROBUSTNESS_RANGE_FLOOR * (hi0 - lo0 + 2 * margin);
    const double mid = 0.5 * (r.lo + r.hi);
    const double allowed = eps * std::max(std::abs(mid), floor_width);
    if (r.hi - r.lo > allowed * (1 + 1e-12) &&
        !(0.5 * (r.lo + r.hi) <= r.lo || 0.5 * (r.lo + r.hi) >= r.hi)) {
        why = "[NOT TIGHT] width " + num(r.hi - r.lo) + " > " + num(allowed);
        return false;
    }
    return true;
}

// The bisection monitor `fn` must bracket the exact robustness `rho` within
// epsilon, and the reference monitors must return `rho` itself (the STL*_1
// reference only for STL*_1 formulas).
inline void check_robustness(const std::string& label, RobustnessFn fn, const Formula& f,
                             const Trace& tr, double rho, double eps = 0.01,
                             const std::string& extra = "") {
    const RobustnessResult r = fn(f, tr, eps);
    std::string why;
    const bool ok = robustness_ok(r, rho, eps, f, tr, why);
    check(ok, label, extra + "\n    exact " + num(rho) + "  got [" + num(r.lo) + ", " +
                         num(r.hi) + "] " + why + "\n" + fmt_trace(tr));

    const double ref_general = robustness_monitor_general_naive(f, tr);
    check(same_robustness(ref_general, rho), label + " (general reference monitor)",
          extra + "\n    exact " + num(rho) + "  reference " + num(ref_general) + "\n" + fmt_trace(tr));
    if (is_fragment_one(f)) {
        const double ref_one = robustness_monitor_one_naive(f, tr);
        check(same_robustness(ref_one, rho), label + " (STL*_1 reference monitor)",
              extra + "\n    exact " + num(rho) + "  reference " + num(ref_one) + "\n" + fmt_trace(tr));
    }
}

// The bisection monitor `fn` must bracket the robustness computed by the
// reference monitor `reference` within epsilon.
inline void check_vs_reference(const std::string& label, RobustnessFn fn, ReferenceFn reference,
                               const Formula& f, const Trace& tr, double eps,
                               const std::string& extra = "") {
    const double rho = reference(f, tr);
    const RobustnessResult r = fn(f, tr, eps);
    std::string why;
    check(robustness_ok(r, rho, eps, f, tr, why), label,
          extra + "\n    reference " + num(rho) + "  got [" + num(r.lo) + ", " + num(r.hi) + "] " +
              why + "\n" + fmt_trace(tr));
}

// Random negation-free formulas on small traces: the exact robustness of the
// oracle must be returned exactly by the reference monitors and bracketed by
// the bisection monitor(s); its sign must agree with the Boolean semantics.
inline void random_robustness(RobustnessFn fn, bool fragment_one, int count, unsigned seed,
                              RobustnessFn other = nullptr) {
    Rng R(seed);
    long infinite = 0;
    Coverage coverage;
    for (int it = 0; it < count; ++it) {
        GenOptions o;
        o.comps = R.uni(1, 2);
        o.fragment_one = fragment_one;
        o.negation_free = true;
        o.max_depth = R.uni(1, 4);
        o.max_freeze_depth = fragment_one ? R.uni(0, 2) : 2;
        o.freeze_bias = fragment_one ? 0.0 : 0.4;
        const OracleFormula F = FormulaGen(R, o).make();
        coverage.add(F);
        const std::string text = to_text(F);
        const Trace tr = random_trace(R, R.uni(1, 8), F.num_signals());
        const Formula f = parse_formula(text);
        const double rho = oracle_robustness(F, tr);
        const double eps = R.pick(std::vector<double>{0.1, 0.01, 0.001});
        const std::string extra = "    formula: " + text + "  eps=" + num(eps);
        check_robustness("random robustness #" + std::to_string(it), fn, f, tr, rho, eps, extra);
        if (other) check_robustness("random robustness (cross) #" + std::to_string(it), other, f, tr, rho, eps, extra);
        // Dyadic data: the reference monitor must agree with the oracle exactly.
        const double ref = fragment_one ? robustness_monitor_one_naive(f, tr)
                                        : robustness_monitor_general_naive(f, tr);
        check(ref == rho, "reference == oracle exactly #" + std::to_string(it),
              extra + "\n    oracle " + num(rho) + "  reference " + num(ref));
        if (std::isinf(rho)) ++infinite;
        // Sign of the robustness vs the Boolean semantics.
        std::vector<double> env(F.num_binders(), 0.0);
        const bool holds = sat(F.root, tr, 0, env);
        const double tol = 1e-9 * std::max(1.0, std::abs(rho));
        if (rho > tol)  check(holds,  "rho > 0 implies the formula holds: " + text);
        if (rho < -tol) check(!holds, "rho < 0 implies the formula fails: " + text);
    }
    std::printf("  %d random formulas (%ld with infinite robustness)\n", count, infinite);
    coverage.require(Coverage::negation_free_operators());
}

// Random negation-free formulas on traces of 1 to 5000 samples, noisy and
// smooth: the bisection monitor must bracket the reference monitor's
// robustness within epsilon.  Formulas with nested dependent freezes use at
// most 100 samples (the general reference re-evaluates them from scratch).
inline void bisection_vs_reference(RobustnessFn fn, ReferenceFn reference, bool fragment_one,
                                   int count, unsigned seed) {
    Rng R(seed);
    int largest = 0;
    for (int it = 0; it < count; ++it) {
        int n = random_trace_size(R);
        GenOptions o;
        o.comps = R.uni(1, 2);
        o.fragment_one = fragment_one;
        o.negation_free = true;
        o.max_depth = R.uni(2, 4);
        o.max_freeze_depth = R.uni(1, 2);
        o.freeze_bias = fragment_one ? 0.2 : 0.4;
        o.allow_until = n <= 500;     // the reference until scans every window
        const std::string text = to_text(FormulaGen(R, o).make());
        const Formula f = parse_formula(text);
        if (!fragment_one && freeze_depth(f) >= 2) n = std::min(n, 100);
        largest = std::max(largest, n);
        const Trace tr = varied_trace(R, n, f.num_signals);
        const double eps = R.pick(std::vector<double>{0.1, 0.01, 0.001, 1e-6});
        check_vs_reference("bisection vs reference #" + std::to_string(it), fn, reference, f, tr, eps,
                           "    formula: " + text + "  eps=" + num(eps) + "  n=" + std::to_string(n));
    }
    std::printf("  %d random formulas on traces of 1 to %d samples\n", count, largest);

    // Chosen formulas on 5000 samples: a smooth signal on a uniform grid and
    // a noisy signal on an irregular grid.
    int fixed = 0;
    for (const char* fs : {"a*. F[0,10] (a - a* > 1)", "G[0,20] (a*. F[0,5] (a > a* + 0.5))"})
        for (bool smooth : {true, false}) {
            check_vs_reference(std::string("n=5000: ") + fs, fn, reference, parse_formula(fs),
                               trace_5000_style(5000, smooth, /*uniform=*/smooth, 88), 0.01);
            ++fixed;
        }
    std::printf("  %d chosen formula/trace pairs with 5000 samples\n", fixed);
}

}  // namespace test_util
