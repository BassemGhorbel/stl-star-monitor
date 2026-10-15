// test_robustness_monitor_general.cpp
//
// General STL* robustness monitor (approximate, by bisection), checked
// against the general reference robustness monitor (exact): the bracket
// [lo, hi] must contain the reference value and meet the epsilon stopping rule.
//   A. Constraints over several frozen values, nested dependent freezes.
//   B. On STL*_1 formulas it agrees with the STL*_1 robustness monitor.
//   C. Infinite robustness, non-monotone constraints, negation.
//   D. Random general formulas: the reference monitor equals the independent
//      oracle exactly, and the bisection brackets it.
//   E. Random formulas on traces of 1 to 5000 samples vs the reference monitor.
// Exit code is non-zero on any failure.
#include "robustness_test_utils.h"

using namespace test_util;

static void rob(const char* formula, const Trace& tr, double eps = 0.001) {
    check_vs_reference(std::string("robustness of ") + formula, robustness_monitor_general,
                       robustness_monitor_general_naive, parse_formula(formula), tr, eps);
}
static Trace U(std::vector<std::vector<double>> s) { return Trace::make_uniform(1.0, std::move(s)); }

// ===========================================================================
static void several_frozen_values() {
    section("A. Several frozen values and nested dependent freezes vs the reference monitor");
    // Two binders at the same time.
    rob("a*_1. a*_2. F (a - a*_1 - a*_2 > 0)", U({{1, 5, 2}}));
    rob("a*_1. a*_2. F (a - a*_1 - a*_2 > 0)", U({{5, 5, 5}}));
    rob("a*_1. a*_2. F (a - a*_1 - a*_2 > 0)", U({{3}}));
    rob("a*_1. a*_2. F (a - a*_1 - a*_2 > 0)", U({{1, 10, 2, 10, 3}}));
    rob("a*_1. a*_2. G (a - a*_1 - a*_2 >= 0)", U({{0, 0, 0}}));
    rob("a*_1. a*_2. G (a - a*_1 - a*_2 >= 0)", U({{3, 1, 2}}));
    rob("a*_1. a*_2. F[0,1] (a - a*_1 - a*_2 > 0)", U({{1, 5, 2, 3}}));
    rob("a*_1. a*_2. F (a - a*_1 + a*_2 > 0)", U({{1, 2, 5}}));
    rob("a*. b*. F (a + b - a* - b* > 0)", U({{1, 3}, {2, 4}}));
    // Nested dependent freezes at different times.
    rob("a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))", U({{1, 2, 5, 3}}));
    rob("a*_1. G[0,1] (a*_2. F (a - a*_1 - a*_2 > 0))", U({{1, 4, 2}}));
    rob("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", U({{0, 1, 2, 3, 4}, {0.5, 0.5, 0.5, 0.5, 0.5}}));
    rob("a*_1. F[0,2] (a*_2. (a >= a*_2 U[0,2] a - a*_1 > 1))", U({{1, 2, 0, 3, 4}}));
}

// ===========================================================================
static void agreement_with_robustness_one() {
    section("B. Agreement with the STL*_1 robustness monitor");
    Rng R(9009);
    for (const char* text : {"a*. F (a - a* > 0)", "a*. G (a - a* >= 0)", "a*. F[0,1] (a - a* > 0)",
                             "a*. G[0,1] (a - a* >= 0)", "a*. F (a - a* < 0)",
                             "((a - 2 > 0) && a*. F (a - a* > 0))", "G[0,2] (a*. F[0,1] (a > a* + 1))"}) {
        const Formula f = parse_formula(text);
        for (int t = 0; t < 20; ++t) {
            const Trace tr = random_trace(R, R.uni(1, 10), 1);
            const RobustnessResult r1 = robustness_monitor_one(f, tr, 0.001);
            const RobustnessResult rg = robustness_monitor_general(f, tr, 0.001);
            // Same decision procedure results, so the brackets must be identical.
            check(r1.lo == rg.lo && r1.hi == rg.hi && r1.monitor_calls == rg.monitor_calls,
                  std::string("STL*_1 vs general robustness: ") + text,
                  "    stl1 [" + num(r1.lo) + "," + num(r1.hi) + "]  gen [" + num(rg.lo) + "," +
                      num(rg.hi) + "]\n" + fmt_trace(tr));
        }
    }
}

// ===========================================================================
static void infinite_and_non_monotone() {
    section("C. Infinite robustness, non-monotone constraints, negation");
    rob("a*_1. F[5,6] (a*_2. F (a - a*_1 - a*_2 > 0))", U({{1, 2}}));
    rob("a*_1. G[5,6] (a*_2. F (a - a*_1 - a*_2 > 0))", U({{1, 2}}));
    // a*_1 * a*_2 is not monotone: the initial bounds must still hold.
    rob("a*_1. a*_2. F (a - a*_1 * a*_2 > 0)", U({{0, -2, 4}}));
    bool threw = false;
    try { robustness_monitor_general(parse_formula("a*. ! (a > a*)"), U({{1}}), 0.01); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "negation is rejected");
    threw = false;
    try { robustness_monitor_general_naive(parse_formula("a*. ! (a > a*)"), U({{1}})); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "negation is rejected by the reference monitor");
}

// ===========================================================================
static void random_formulas_vs_oracle() {
    section("D. Random general formulas: reference == oracle, bisection brackets it");
    random_robustness(robustness_monitor_general, /*fragment_one=*/false, 5000, 10010);
}

// ===========================================================================
static void traces_of_different_sizes() {
    section("E. Traces of 1 to 5000 samples, noisy and smooth: bisection vs reference");
    bisection_vs_reference(robustness_monitor_general, robustness_monitor_general_naive,
                           /*fragment_one=*/false, 120, 13013);
}

int main() {
    std::printf("=== General STL* robustness monitor ===\n");
    several_frozen_values();
    agreement_with_robustness_one();
    infinite_and_non_monotone();
    random_formulas_vs_oracle();
    traces_of_different_sizes();
    return finish("test_robustness_monitor_general");
}
