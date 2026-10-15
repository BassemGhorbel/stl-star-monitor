// examples.cpp -- how to use the STL* monitors.
//
//   1. Boolean monitoring of an STL*_1 formula
//   2. Boolean monitoring of a general STL* formula (nested dependent freezes)
//   3. Non-uniform traces (explicit timestamps)
//   4. Robustness
//   5. Errors
//   6. Low-level API: parsed formulas, traces, formulas built in code
//
// Build with the library (target `examples`) and run it.
#include "stl_star/stl_star_monitor.h"
#include "stl_star/naive_monitors.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

// Every name of the library lives in namespace stl_star, e.g.
// stl_star::monitor_one(...).  This example imports them all for brevity.
using namespace stl_star;

static void print(const char* formula, const MonitorResult& r) {
    std::printf("  %-45s %s at t0, holds on ", formula, r.holds_at_start ? "holds" : "fails");
    if (r.intervals.empty()) std::printf("(nowhere)");
    for (const Interval& iv : r.intervals) std::printf("[%g, %g] ", iv.lo, iv.hi);
    std::printf("\n");
}

int main() {
    // -----------------------------------------------------------------------
    std::printf("1. STL*_1 formulas: monitor_one\n");
    // Signal 'a' sampled every dt = 1 time unit.
    const Signal rising = {{0.0, 0.5, 1.2, 2.0, 2.8, 3.5}};

    // "During [0,3], whenever we look, a rises by at least 0.5 within 2 time units."
    print("G[0,3] a*. F[0,2] (a - a* >= 0.5)",
          monitor_one("G[0,3] a*. F[0,2] (a - a* >= 0.5)", rising, 1.0));
    // "a never drops below its current value."
    print("a*. G (a - a* >= 0)", monitor_one("a*. G (a - a* >= 0)", rising, 1.0));
    // Plain STL (no freeze operator) works too.
    print("F[0,2] a > 2", monitor_one("F[0,2] a > 2", rising, 1.0));

    // -----------------------------------------------------------------------
    std::printf("\n2. General STL* formulas: monitor_general\n");
    // One signal per letter of the formula, in alphabetical order of the
    // letters: a = signal[0], b = signal[1].
    const Signal two = {{0.0, 1.0, 2.0, 3.0, 4.0},
                        {0.5, 0.5, 0.5, 0.5, 0.5}};
    // The constraint reads a* (outer freeze) and b* (inner freeze): not STL*_1.
    print("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)",
          monitor_general("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", two, 1.0));
    // Two frozen values of the same signal: some later value exceeds their sum.
    print("a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))",
          monitor_general("a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))", {{1, 2, 5, 3, 4}}, 1.0));
    // Any lowercase letters can name signals: here x = signal[0] (the speed)
    // and z = signal[1] (the limit), since x comes before z.
    const Signal speed_and_limit = {{1.0, 2.0, 4.0, 3.0, 5.0},
                                    {3.0, 3.0, 4.5, 4.5, 4.5}};
    print("x*. F[0,2] (x - x* >= 1 && x <= z)",
          monitor_general("x*. F[0,2] (x - x* >= 1 && x <= z)", speed_and_limit, 1.0));

    // -----------------------------------------------------------------------
    std::printf("\n3. Non-uniform traces: pass timestamps instead of dt\n");
    const std::vector<double> times = {0.0, 0.5, 1.2, 2.8, 4.0};
    print("a*. G[0,4] (a - a* >= 0)",
          monitor_one("a*. G[0,4] (a - a* >= 0)", {{1.0, 1.5, 2.0, 2.5, 3.0}}, times));

    // -----------------------------------------------------------------------
    std::printf("\n4. Robustness: how strongly the formula holds at t0\n");
    // The result is a bracket [lo, hi] around the robustness; epsilon is the
    // relative error (here 0.1%).
    RobustnessResult r = robustness_one("G[0,3] a*. F[0,2] (a - a* >= 0)", rising, 1.0, 0.001);
    std::printf("  %-45s rho in [%g, %g] (%d monitor calls)\n",
                "G[0,3] a*. F[0,2] (a - a* >= 0)", r.lo, r.hi, r.monitor_calls);
    r = robustness_general("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", two, 1.0, 0.001);
    std::printf("  %-45s rho in [%g, %g] (%d monitor calls)\n",
                "a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", r.lo, r.hi, r.monitor_calls);
    // An empty window gives infinite robustness: lo = hi = -inf.
    r = robustness_one("F[10,20] a > 0", rising, 1.0);
    std::printf("  %-45s rho in [%g, %g]\n", "F[10,20] a > 0", r.lo, r.hi);

    // -----------------------------------------------------------------------
    std::printf("\n5. Errors\n");
    try {   // a general formula given to the STL*_1 monitor
        monitor_one("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", two, 1.0);
    } catch (const std::invalid_argument& e) { std::printf("  invalid_argument: %s\n", e.what()); }
    try {   // syntax error: binary operators need parentheses
        monitor_one("a > 0 && a < 5", rising, 1.0);
    } catch (const std::runtime_error& e) { std::printf("  runtime_error:    %s\n", e.what()); }
    try {   // the formula reads two signals (a and b) but the trace has one
        monitor_one("a > b", rising, 1.0);
    } catch (const std::invalid_argument& e) { std::printf("  invalid_argument: %s\n", e.what()); }

    // -----------------------------------------------------------------------
    std::printf("\n6. Low-level API\n");
    // Parse once and reuse the formula on several traces.
    const Formula f = parse_formula("a*. F[0,2] (a - a* > 1)");
    const Trace t1 = Trace::make_uniform(0.5, {{0, 1, 3, 2, 4}});
    const Trace t2 = Trace::make_non_uniform({0, 1, 3, 4, 7}, {{0, 1, 3, 2, 4}});
    for (const Trace* t : {&t1, &t2}) {
        const Intervals fast = boolean_monitor_one(f, *t);
        const Intervals ref  = boolean_monitor_one_naive(f, *t);   // reference monitor
        std::printf("  %-45s", "a*. F[0,2] (a - a* > 1)");
        for (const Interval& iv : fast) std::printf("[%g, %g] ", iv.lo, iv.hi);
        const bool same = fast.size() == ref.size() &&
                          std::equal(fast.begin(), fast.end(), ref.begin(), [](const Interval& x, const Interval& y) {
                              return x.lo == y.lo && x.hi == y.hi;
                          });
        std::printf("(%s the reference monitor)\n", same ? "same as" : "DIFFERENT from");
    }

    // Build a formula in code, e.g. with a function the parser does not know:
    //   G[0,0.5] a*. (sin(a) >= a* - 1)
    Formula g;
    g.num_freeze_vars = 1;
    const int c = g.add_node(FormulaNode::make_signal_constraint(
        [](const std::vector<double>& s) { return std::sin(s[0]); },   // h(s), signal part
        [](const std::vector<double>& x) { return x[0] - 1.0; },       // g(x*), frozen part
        CompOp::GE, /*freeze variables read by g=*/{0}));
    const int fr = g.add_node(FormulaNode::make_freeze(/*variable=*/0, /*signal=*/0, c));
    g.add_node(FormulaNode::make_always(0.0, 0.5, fr));   // the last node is the root
    g.validate();
    std::printf("  %-45s %s at t0\n", "G[0,0.5] a*. (sin(a) >= a* - 1)",
                boolean_monitor_one_holds_at_start(g, Trace::make_uniform(0.1, {{0, 0.2, 0.4, 0.6, 0.8, 1.0}}))
                    ? "holds" : "fails");
    return 0;
}
