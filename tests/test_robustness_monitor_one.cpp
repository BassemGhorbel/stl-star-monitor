// test_robustness_monitor_one.cpp
//
// STL*_1 robustness monitor (approximate, by bisection), checked against the
// STL*_1 reference robustness monitor (exact): the bracket [lo, hi] must
// contain the reference value and meet the epsilon stopping rule.
//   A. Every operator: predicates, Boolean and temporal operators, until,
//      freeze, nesting, non-uniform traces, windows between samples.
//   B. Infinite robustness, non-monotone constraints, extreme values, negation.
//   C. Epsilon: relative stopping criterion, nested brackets, call counts.
//   D. Random STL*_1 formulas: the reference monitor equals the independent
//      oracle exactly, and the bisection brackets it (also cross-checked with
//      the general robustness monitor).
//   E. Random formulas on traces of 1 to 5000 samples vs the reference monitor.
// Exit code is non-zero on any failure.
#include "robustness_test_utils.h"

using namespace test_util;

static void rob(const char* formula, const Trace& tr, double eps = 0.01) {
    check_vs_reference(std::string("robustness of ") + formula, robustness_monitor_one,
                       robustness_monitor_one_naive, parse_formula(formula), tr, eps);
}
static Trace U(std::vector<double> a) { return Trace::make_uniform(1.0, {std::move(a)}); }

// ===========================================================================
static void operators() {
    section("A. Every operator vs the reference monitor");
    // Predicates.
    rob("a > 0",       U({7, 7, 7}));
    rob("a > 0",       U({-3, -3}));
    rob("a - 5 > 0",   U({8, 8}));
    rob("a < 10",      U({6, 6}));
    rob("a <= 0",      U({-4}));
    // Timed F and G.
    rob("F[0,3] a > 0",      U({1, 5, 3, 2, 4}));
    rob("F[0,3] a > 0",      U({-1, -3, -2, -4}));
    rob("F[1,4] a > 0",      U({2, 5, 3, 7, 4}));
    rob("F[2,4] a > 0",      U({9, 1, 3, 5, 2, 8}));
    rob("G[0,3] a > 0",      U({3, 1, -2, 4, 2}));
    rob("G[1,3] a > 0",      U({8, 2, 4, 3, 6}));
    rob("G[2,5] a > 0",      U({9, 8, 3, 5, 2, 4, 1}));
    // Untimed F and G.
    rob("G a > 0",           U({3, 5, -2, 4, 1}));
    rob("F a - 7 > 0",       U({2, 5, 9, 3, 4}));
    // Boolean.
    rob("(a > 2 && a < 8)",  U({3, 3}));
    rob("(a > 5 || a < 2)",  U({3, 3}));
    // Until, strict and non-strict, timed and untimed.
    rob("(a > 0 U[0,2] a > 3)",  U({3, 1, 4, 2}));
    rob("(a > 0 UNS[0,2] a > 3)", U({3, 1, 4, 2}));
    rob("(a > 0 U a > 3)",       U({3, 1, -1, 4}));
    rob("(a > 0 UNS a > 3)",      U({3, 1, 2, 4}));
    // Freeze.
    rob("a*. F[0,4] (a - a* > 0)", U({2, 5, 3, 8, 4, 1}));
    rob("a*. F[1,3] (a - a* > 0)", U({7, 2, 3, 4, 9}));
    rob("a*. G[1,3] (a - a* > 0)", U({2, 5, 4, 6, 3}));
    rob("a*. G (a - a* > 0)",      U({5, 3, 4, 6, 2}));
    rob("a*. F (a - a* > 0)",      U({3, 7, 2, 5, 8, 1}));
    // Nesting.
    rob("F[0,2] G[0,1] (a > 0)",               U({3, 5, 2, 4, 1}));
    rob("G[0,2] F[0,1] (a > 0)",               U({2, 5, 1, 4, 3}));
    rob("a*. F[0,3] G[0,1] (a - a* > 0)",      U({3, 5, 4, 6, 2, 1}));
    rob("G[0,2] (a*. F[0,1] (a - a* > 0))",    U({4, 7, 2, 5, 3}));
    rob("a*. G[0,3] F[0,1] (a - a* > 0)",      U({2, 5, 3, 7, 4, 6}));
    // Non-uniform traces, and a window that falls between samples.
    rob("F[0,3] a > 0",  Trace::make_non_uniform({0, 1, 3, 5}, {{2, 7, 4, 1}}));
    rob("G[1,3] a > 0",  Trace::make_non_uniform({0, 0.5, 1.5, 3, 4}, {{8, 3, 5, 2, 6}}));
    rob("a*. F[1,4] (a - a* > 0)", Trace::make_non_uniform({0, 1, 2, 4, 5}, {{3, 5, 2, 6, 4}}));
    rob("F[0,2] F[1.2,1.8] a > 0", Trace::make_non_uniform({0, 1, 2.5}, {{4, 6, 9}}));
    // Two signals.
    rob("b*. F (a - b* > 0)", Trace::make_uniform(1.0, {{1, 3, 5, 2}, {2, 1, 3, 4}}));
}

// ===========================================================================
static void infinite_non_monotone_extreme() {
    section("B. Infinite robustness, non-monotone constraints, extreme values, negation");
    // Empty windows: -inf for F, +inf for G, and their combinations.
    rob("F[5,10] a > 0",             U({1, 2, 3}));
    rob("G[5,10] a > 0",             U({1, 2, 3}));
    rob("(G[5,10] a > 0 && a > 1)",  U({4, 2, 3}));
    rob("(F[5,10] a > 0 || a > 1)",  U({4, 2, 3}));
    rob("a*. G[3,4] (a - a* > 0)",   U({1, 2}));
    // a* * a* is not monotone in a*: the initial bounds must still hold.
    rob("a*. F (a - a* * a* > -3)",  U({0, 5, -2, 3}));
    // Large magnitudes: differences of nearby huge values.
    check_vs_reference("large frozen values (1e16)", robustness_monitor_one,
                       robustness_monitor_one_naive,
                       parse_formula("(a*. F[0,1] (a - a* >= 0) && b >= -1000)"),
                       Trace::make_uniform(1.0, {{1e16, 1e16 + 2}, {1000, 1000}}), 0.01);
    {
        const double t0 = 1.7e9, t1 = std::nextafter(std::nextafter(1.7e9, 2e9), 2e9);
        check_vs_reference("Unix-time-like values", robustness_monitor_one,
                           robustness_monitor_one_naive,
                           parse_formula("a*. F[0,1] (a - a* >= 0)"), U({t0, t1}), 1e-6);
    }
    // Near the largest double: bounds are clamped, and an overflowing
    // difference gives +inf.
    rob("a*. F (a - a* > 0)", U({-6e307, 9e307}));
    rob("a*. F (a - a* > 0)", U({-1e308, 1e308}));
    // Negation is rejected by both monitors.
    bool threw = false;
    try { robustness_monitor_one(parse_formula("! (a > 0)"), U({1}), 0.01); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "negation is rejected");
    threw = false;
    try { robustness_monitor_one_naive(parse_formula("! (a > 0)"), U({1})); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "negation is rejected by the reference monitor");
}

// ===========================================================================
static void epsilon() {
    section("C. Epsilon");
    const Formula f = parse_formula("a*. F[0,5] (a - a* > 0)");
    const Trace tr = U({2, 5, 8, 3, 6, 1, 9});
    RobustnessResult prev{-1e300, 1e300, 0};
    for (double eps : {0.1, 0.01, 0.001, 1e-6}) {
        const RobustnessResult r = robustness_monitor_one(f, tr, eps);
        check_vs_reference("eps=" + num(eps), robustness_monitor_one, robustness_monitor_one_naive,
                           f, tr, eps);
        check(r.lo >= prev.lo - 1e-12 && r.hi <= prev.hi + 1e-12 && r.monitor_calls >= prev.monitor_calls,
              "smaller epsilon gives a nested bracket with more calls (eps=" + num(eps) + ")");
        check(r.monitor_calls <= 2 + 64, "bisection terminates quickly");
        std::printf("  eps=%-6g [%.9f, %.9f] %d calls\n", eps, r.lo, r.hi, r.monitor_calls);
        prev = r;
    }
    bool threw = false;
    try { robustness_monitor_one(f, tr, 0.0); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "epsilon <= 0 is rejected");
}

// ===========================================================================
static void random_formulas_vs_oracle() {
    section("D. Random STL*_1 formulas: reference == oracle, bisection brackets it");
    random_robustness(robustness_monitor_one, /*fragment_one=*/true, 6000, 8008,
                      robustness_monitor_general);
}

// ===========================================================================
static void traces_of_different_sizes() {
    section("E. Traces of 1 to 5000 samples, noisy and smooth: bisection vs reference");
    bisection_vs_reference(robustness_monitor_one, robustness_monitor_one_naive,
                           /*fragment_one=*/true, 150, 12012);
}

int main() {
    std::printf("=== STL*_1 robustness monitor ===\n");
    operators();
    infinite_non_monotone_extreme();
    epsilon();
    random_formulas_vs_oracle();
    traces_of_different_sizes();
    return finish("test_robustness_monitor_one");
}
