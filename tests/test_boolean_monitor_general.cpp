// test_boolean_monitor_general.cpp
//
// Correctness and complexity of the general STL* monitors (fast and naive),
// i.e. formulas with nested dependent freeze operators:
//   A. Corner cases: fast vs the reference monitor.
//   B. Random general formulas (up to 3 nested freeze operators, constraints
//      reading outer variables): oracle == naive == fast.
//   C. The STL*_1 monitors reject formulas outside STL*_1.
//   D. Traces of 1 to 5000 samples, noisy and smooth: naive == fast.
//   E. Scaling with two nested dependent freeze variables (~n^2 log n).
// Exit code is non-zero on any failure.
#include "random_generators.h"

using namespace test_util;

static void check_general(const std::string& label, const Formula& f, const Trace& tr,
                          const Intervals& expected, const std::string& extra = "");

// The fast general monitor must return the result of the reference monitor.
static void check_vs_reference(const std::string& label, const Formula& f, const Trace& tr,
                               const std::string& extra = "") {
    check_general(label, f, tr, boolean_monitor_general_naive(f, tr), extra);
}

static void check_general(const std::string& label, const Formula& f, const Trace& tr,
                          const Intervals& expected, const std::string& extra) {
    const Intervals ng = boolean_monitor_general_naive(f, tr);
    const Intervals fg = boolean_monitor_general(f, tr);
    const bool start = !expected.empty() && expected.front().lo == tr.timestamp(0);
    const bool ok = same(ng, expected) && same(fg, expected) &&
                    boolean_monitor_general_holds_at_start(f, tr) == start &&
                    boolean_monitor_general_naive_holds_at_start(f, tr) == start;
    check(ok, label,
          extra + "\n    expected   " + fmt(expected) + "\n    gen naive  " + fmt(ng) +
          "\n    gen fast   " + fmt(fg) + "\n" + fmt_trace(tr));
}

static Trace U(std::vector<std::vector<double>> s) { return Trace::make_uniform(1.0, std::move(s)); }

// ===========================================================================
static void corner_cases() {
    section("A. Corner cases: fast vs the reference monitor");
    struct Case { const char* formula; Trace trace; };
    const std::vector<Case> cases = {
        // the README example: the constraint reads a* and b*
        {"a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", U({{0, 1, 2, 3, 4}, {0.5, 0.5, 0.5, 0.5, 0.5}})},
        // nested dependent freezes of the same signal
        {"a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))", U({{1, 2, 5, 3, 4}})},
        // both binders at the same position
        {"a*. b*. F (a + b - a* - b* > 0)", U({{2, 1, 4}, {1, 3, 2}})},
        // sibling freeze operators
        {"(a*. F (a - a* > 0) && b*. G (b - b* >= 0))", U({{1, 3, 2, 4}, {2, 1, 3, 3}})},
        // the inner scope reads the outer variable a*_1 as well as its own
        {"a*_1. F (a*_2. F[0,1] ((a - a*_1 >= 0) && (a - a*_2 > 0)))", U({{1, 3, 2, 4}})},
        // an inner scope that reads only its own variable, inside an outer scope
        {"a*_1. F ((a > a*_1 && b*_2. G[0,1] (b - b*_2 >= 0)))", U({{1, 3, 2, 4}, {2, 1, 3, 3}})},
        // non-uniform, two nested timed windows
        {"a*_1. F[0,2] (a*_2. F[0,1] (a - a*_1 - a*_2 >= 0))",
         Trace::make_non_uniform({0, 0.5, 1.5, 3, 3.5}, {{1, 0, 2, 4, 1}})},
        // the same letter frozen twice; the outer value is unused
        {"a*_1. G (a*_2. G (a - a*_2 >= 0))", U({{3, 1, 2, 2}})},
        // until inside nested freezes
        {"a*_1. F[0,2] (a*_2. (a >= a*_2 U[0,2] a - a*_1 > 1))", U({{1, 2, 0, 3, 4}})},
        // single sample
        {"a*_1. F (a*_2. G (a - a*_1 - a*_2 >= 0))", U({{0}})},
        {"a*_1. F (a*_2. G (a - a*_1 - a*_2 > 0))", U({{0}})},
    };
    for (const Case& c : cases)
        check_vs_reference(std::string("corner case: ") + c.formula, parse_formula(c.formula), c.trace);
    std::printf("  %zu corner cases\n", cases.size());
}

// ===========================================================================
static void random_formulas_vs_oracle() {
    section("B. Random general formulas vs the oracle");
    Rng R(6006);
    long cases = 0, outside_fragment_one = 0;
    Coverage coverage;
    for (int it = 0; it < 20000; ++it) {
        GenOptions o;
        o.comps = R.uni(1, 3);
        o.fragment_one = false;
        o.max_depth = R.uni(2, 5);
        o.max_freeze_depth = R.coin(0.15) ? 3 : 2;
        o.freeze_bias = R.coin(0.7) ? 0.4 : 0.0;
        const OracleFormula F = FormulaGen(R, o).make();
        coverage.add(F);
        const std::string text = to_text(F);
        const int n = o.max_freeze_depth == 3 ? R.uni(1, 7) : R.uni(1, 10);
        const Trace tr = random_trace(R, n, F.num_signals());
        Formula f;
        try { f = parse_formula(text); }
        catch (const std::exception& e) { check(false, "parse: " + text, e.what()); continue; }
        check_general("random general #" + std::to_string(it), f, tr, oracle_intervals(F, tr),
                      "    formula: " + text);
        if (!is_fragment_one(f)) ++outside_fragment_one;
        ++cases;
    }
    std::printf("  %ld random formulas x traces (%ld outside STL*_1)\n", cases, outside_fragment_one);
    coverage.require(Coverage::all_operators());
}

// ===========================================================================
static void monitor_one_rejects_general_formulas() {
    section("C. STL*_1 monitors reject general formulas");
    for (const char* text : {"a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))",
                             "a*. G (b*. F (a - a* > 0))",            // reads the outer variable only
                             "a*. b*. F (a + b - a* - b* > 0)"}) {
        const Formula f = parse_formula(text);
        const Trace tr = f.num_signals == 1 ? U({{1, 2, 3}}) : U({{1, 2, 3}, {3, 2, 1}});
        bool threw1 = false, threw2 = false;
        try { boolean_monitor_one(f, tr); } catch (const std::invalid_argument&) { threw1 = true; }
        try { boolean_monitor_one_naive(f, tr); } catch (const std::invalid_argument&) { threw2 = true; }
        check(threw1 && threw2 && !is_fragment_one(f), std::string("rejected: ") + text);
        // ... while the general monitors accept them.
        check(same(boolean_monitor_general(f, tr), boolean_monitor_general_naive(f, tr)),
              std::string("general accepts: ") + text);
    }
}

// ===========================================================================
static void traces_of_different_sizes() {
    section("D. Traces of 1 to 5000 samples, noisy and smooth: naive vs fast");
    Rng R(7007);
    long cases = 0;
    int largest_nested = 0;
    for (int it = 0; it < 120; ++it) {
        int n = random_trace_size(R);
        GenOptions o;
        o.comps = R.uni(1, 2);
        o.fragment_one = false;
        o.max_depth = R.uni(2, 4);
        o.max_freeze_depth = R.uni(1, 2);
        o.freeze_bias = 0.4;
        o.allow_until = n <= 500;     // the reference until scans every window
        const std::string text = to_text(FormulaGen(R, o).make());
        const Formula f = parse_formula(text);
        // The reference monitor re-evaluates a nested freeze from scratch at
        // every outer step, so nested formulas use at most 100 samples.
        if (freeze_depth(f) >= 2) { n = std::min(n, 100); largest_nested = std::max(largest_nested, n); }
        const Trace tr = varied_trace(R, n, f.num_signals);
        check_vs_reference("size test #" + std::to_string(it), f, tr,
                           "    formula: " + text + "  n=" + std::to_string(n));
        ++cases;
    }
    // Chosen formulas on 5000 samples (freeze operators not nested in each
    // other): a smooth signal on a uniform grid, a noisy one on an irregular grid.
    for (const char* fs : {"G[0,20] (a*. F[0,5] (a > a* + 0.5))",
                           "(a*. F[0,10] (a - a* > 1) && a*. G[0,3] (a - a* >= -2))"})
        for (bool smooth : {true, false}) {
            check_vs_reference(std::string("n=5000: ") + fs, parse_formula(fs),
                               trace_5000_style(5000, smooth, /*uniform=*/smooth, 77));
            ++cases;
        }
    std::printf("  %ld cases (nested dependent freezes up to %d samples)\n", cases, largest_nested);
}

// ===========================================================================
static void scaling() {
    section("E. Scaling, two nested dependent freezes (n^2 log n gives ~x4.4 per doubling)");
    const Formula f = parse_formula("a*_1. F[0,30] (a*_2. F[0,10] (a - a*_1 - a*_2 > 0.5))");
    std::vector<double> ms;
    const std::vector<int> sizes = {400, 800, 1600};
    for (int n : sizes) {
        std::vector<double> sig(n);
        for (int i = 0; i < n; ++i) sig[i] = std::sin(i * 6.283 * 3 / n) * 2 + 2;
        const Trace tr = Trace::make_uniform(0.5, {sig});
        ms.push_back(best_time_ms([&] { boolean_monitor_general(f, tr); }));
    }
    const double g = std::sqrt(ms[2] / std::max(ms[0], 1e-3));
    std::printf("  n=%d:%.1fms n=%d:%.1fms n=%d:%.1fms  x%.2f per doubling\n",
                sizes[0], ms[0], sizes[1], ms[1], sizes[2], ms[2], g);
    // Three nested freezes (depth 3): fast == naive on a moderate trace.
    const Formula f3 = parse_formula(
        "a*_1. F[0,2] (a*_2. (a >= a*_2 && F[0,20] (a*_3. F[0,20] (a - a*_3 > a*_1 - a*_2 + 1))))");
    {
        std::vector<double> sig(60);
        for (int i = 0; i < 60; ++i) sig[i] = std::sin(i * 0.4) * 3;
        const Trace tr = Trace::make_uniform(1.0, {sig});
        check_vs_reference("three nested freezes", f3, tr);
    }
#ifdef NDEBUG
    check(g < 6.0, "general monitor with 2 nested freezes grows like n^2 log n (" + num(g) + " < 6)");
#else
    std::printf("  (scaling is only asserted in optimized builds)\n");
#endif

    // The scaling formula itself: fast == reference.
    std::vector<double> sig(200);
    for (int i = 0; i < 200; ++i) sig[i] = std::sin(i * 0.1) * 2 + 2;
    check_vs_reference("scaling formula, n=200", f, Trace::make_uniform(0.5, {sig}));
}

int main() {
    std::printf("=== General STL* monitors: correctness and complexity ===\n");
    corner_cases();
    random_formulas_vs_oracle();
    monitor_one_rejects_general_formulas();
    traces_of_different_sizes();
    scaling();
    return finish("test_boolean_monitor_general");
}
