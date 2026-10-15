// test_boolean_monitor_one.cpp
//
// Correctness and complexity of the STL*_1 monitors (fast and naive):
//   A. Corner cases: every monitor vs the reference monitor.
//   B. Interval operators vs their pointwise definitions (random sets/windows).
//   C. SignalConstraintUpdater vs brute force (ties, NaN/inf keys and thresholds, resets).
//   D. Random STL*_1 formulas: oracle == naive == fast (both STL*_1 and the
//      general monitors), on uniform and non-uniform traces.
//   E. Traces of 1 to 5000 samples, noisy and smooth: naive == fast.
//   F. Scaling of the fast monitor (near n log n where the theory says so).
//   G. Timestamps near 1e9, an inexact dt, values near 1e16.
// Exit code is non-zero on any failure.
#include "stl_star/signal_constraint_updater.h"
#include "random_generators.h"

using namespace test_util;

// Run all four Boolean monitors that accept STL*_1 formulas and compare them
// with `expected`.
static void check_all_monitors(const std::string& label, const Formula& f, const Trace& tr,
                               const Intervals& expected, const std::string& extra = "") {
    const Intervals n1 = boolean_monitor_one_naive(f, tr);
    const Intervals f1 = boolean_monitor_one(f, tr);
    const Intervals ng = boolean_monitor_general_naive(f, tr);
    const Intervals fg = boolean_monitor_general(f, tr);
    const bool start = !expected.empty() && expected.front().lo == tr.timestamp(0);
    const bool ok = same(n1, expected) && same(f1, expected) && same(ng, expected) &&
                    same(fg, expected) &&
                    boolean_monitor_one_holds_at_start(f, tr) == start &&
                    boolean_monitor_one_naive_holds_at_start(f, tr) == start &&
                    boolean_monitor_general_holds_at_start(f, tr) == start &&
                    boolean_monitor_general_naive_holds_at_start(f, tr) == start;
    check(ok, label,
          extra + "\n    expected   " + fmt(expected) + "\n    s1 naive   " + fmt(n1) +
          "\n    s1 fast    " + fmt(f1) + "\n    gen naive  " + fmt(ng) +
          "\n    gen fast   " + fmt(fg) + "\n" + fmt_trace(tr));
}

static Trace U(double dt, std::vector<std::vector<double>> s) { return Trace::make_uniform(dt, std::move(s)); }
static Trace NU(std::vector<double> t, std::vector<std::vector<double>> s) {
    return Trace::make_non_uniform(std::move(t), std::move(s));
}

// Every monitor must return the result of the STL*_1 reference monitor.
static void check_vs_reference(const std::string& label, const Formula& f, const Trace& tr,
                               const std::string& extra = "") {
    check_all_monitors(label, f, tr, boolean_monitor_one_naive(f, tr), extra);
}

// ===========================================================================
// A. Corner cases vs the reference monitor
// ===========================================================================
static void corner_cases() {
    section("A. Corner cases: all monitors vs the reference monitor");
    struct Case { const char* formula; Trace trace; };
    const Trace T0 = U(1, {{1, 3, 5, 7, 9}});
    const Trace T1 = U(1, {{3, 7, 2, 8, 5, 1, 9, 4, 6, 2}});
    const Trace T2 = NU({0, 1, 2.5, 4, 5.5, 7, 8, 10}, {{3, 8, 2, 9, 5, 1, 7, 4}});
    const Trace T3 = U(0.1, {{0, 1, 2, 3, 4}});
    const Trace T4 = U(1, {{2, 2, 2, 3, 2}});
    const Trace T5 = U(1, {{5}});
    const Trace T6 = U(1, {{1, 3, 5, 2, 4}, {2, 1, 3, 4, 1}});
    const Trace T7 = U(1, {{0.0, 0.5, 1.2, 2.0, 2.8, 3.5}});

    const std::vector<Case> cases = {
        // --- predicates / Boolean ------------------------------------------------
        {"a >= 5",                          T0}, {"! a >= 5", T0},
        {"(a >= 3 && a <= 7)",              T0}, {"(a < 3 || a > 7)", T0},
        {"a >= 5",                          T1}, {"! (a >= 5)", T1},
        {"(a >= 3 && a <= 7)",              T1}, {"(a < 3 || a > 7)", T1},
        {"(a + 1) * 2 >= 10",               T1},
        // --- timed operators, uniform ------------------------------------------------
        {"F[0,2] a >= 7",                   T1}, {"G[0,2] a >= 3", T1},
        {"F[1,3] a >= 8",                   T1},
        {"G[1,2] a >= 4",                   T1},   // empty window at the last sample
        {"F[0.3,0.6] a >= 0",               T1},   // no sample in any window
        {"G[0.3,0.6] a >= 100",             T1},   // vacuously true
        {"F[0.5,1.5] a >= 8",               T1},
        {"F[2,inf] a >= 9",                 T1},
        // --- until -----------------------------------------------------------------------
        {"(a >= 3 U[0,2] a >= 8)",          T1},
        {"(a < 8 U[0,2] a >= 8)",           T1},   // strict: phi not needed at j
        {"(a < 8 UNS[0,2] a >= 8)",          T1},   // non-strict: phi needed at j
        {"(a > 0 U a >= 9)",                T1},
        // --- untimed -------------------------------------------------------------------
        {"G a >= 2",                        T1}, {"F a >= 9", T1},
        // --- freeze --------------------------------------------------------------------
        {"a*. F (a - a* > 0)",              T1}, {"a*. G (a - a* >= 0)", T1},
        {"a*. F[0,2] (a - a* >= 2)",        T1}, {"a*. F[0,2] (a >= a* + 2)", T1},
        {"F[0,1] (a*. F[0,2] (a - a* >= 2))", T1},
        {"G (a*. F[0,2] (a - a* >= 2))",    T1}, {"F (a*. F[0,2] (a - a* >= 2))", T1},
        {"(a >= 5 && a*. F (a - a* > 0))",  T1}, {"(a*. F (a - a* > 0) && a >= 5)", T1},
        {"(a >= 9 || a*. G (a - a* >= 0))", T1}, {"! (a*. F (a - a* > 0))", T1},
        {"a*. F (a*. G (a - a* >= 0))",     T1},   // inner a* shadows the outer one
        {"a*. (a >= 3 U[0,3] a - a* >= 4)", T1},
        {"G[0,3] a*. F[0,2] (a - a* >= 0.5)", T7},
        // --- non-uniform -----------------------------------------------------------------
        {"F[0,1] a >= 7",                   T2},
        {"F[0.5,1] a <= 5",                 T2},   // windows miss the samples
        {"G[0.5,1] a >= 100",               T2},   // vacuous where the window is empty
        {"(a >= 2 U[1,3] a >= 9)",          T2},
        {"a*. F[0,2.5] (a - a* > 3)",       T2},
        // --- dt = 0.1: tau_3 = 0.30000000000000004 must count as 0.3 -----------------------
        {"F[0.3,0.3] a >= 3",               T3}, {"G[0,0.2] a >= 1", T3},
        // --- ties: >= vs > -------------------------------------------------------------
        {"a*. F[1,1] (a - a* >= 0)",        T4}, {"a*. F[1,1] (a - a* > 0)", T4},
        {"a*. G (a <= a*)",                 T4},
        // --- single sample ----------------------------------------------------------------
        {"F[1,2] a >= 0",                   T5}, {"G[1,2] a < 0", T5},
        {"a*. G (a - a* >= 0)",             T5},
        {"(a > 9 U a >= 5)",                T5}, {"(a > 9 UNS a >= 5)", T5},
        // --- two components --------------------------------------------------------------
        {"b*. F (a - b* > 0)",              T6}, {"a*. G[0,1] (b >= a* - 1)", T6},
    };
    for (const Case& c : cases)
        check_vs_reference(std::string("corner case: ") + c.formula, parse_formula(c.formula), c.trace);
    std::printf("  %zu corner cases\n", cases.size());
}

// ===========================================================================
// B. Interval operators vs pointwise definitions
// ===========================================================================
static void interval_operators_vs_pointwise() {
    section("B. Interval operators vs pointwise definitions");
    Rng R(1001);
    const double inf = std::numeric_limits<double>::infinity();
    long cases = 0;
    for (int it = 0; it < 5000; ++it) {
        // Mostly short traces; some longer ones, with sparse and dense sets
        // (one list much shorter than the other).
        const int n = it % 5 == 4 ? R.uni(26, 300) : R.uni(1, 25);
        const Trace tr  = random_trace(R, n, 1);
        auto random_set = [&]() {
            std::vector<char> p(n);
            const double density = R.coin(0.3) ? R.pick(std::vector<double>{0.01, 0.03, 0.97, 0.99})
                                               : R.real(0.0, 1.0);
            for (auto& x : p) x = R.coin(density);
            return p;
        };
        const std::vector<char> P = random_set(), Q = random_set();
        const IndexIntervals A = points_to_intervals(P), B = points_to_intervals(Q);
        double a = R.pick(std::vector<double>{0, 0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3.75});
        double b = a + R.pick(std::vector<double>{0, 0.05, 0.25, 0.5, 1, 2, 4, inf});
        const IndexIntervals W = window_set(tr, a, b);

        auto expect_pts = [&](const std::function<bool(int)>& pred) {
            std::vector<char> e(n);
            for (int i = 0; i < n; ++i) e[i] = pred(i) ? 1 : 0;
            return points_to_intervals(e);
        };
        auto eq = [](const IndexIntervals& x, const IndexIntervals& y) {
            if (x.size() != y.size()) return false;
            for (size_t k = 0; k < x.size(); ++k)
                if (x[k].lo != y[k].lo || x[k].hi != y[k].hi) return false;
            return true;
        };
        auto in_w = [&](int i, int k) { return in_window(tr, i, k, a, b); };
        const std::string ctx = "\n    window [" + num(a) + "," + num(b) + "]\n" + fmt_trace(tr);

        check(eq(complement_of(A, n), expect_pts([&](int i) { return !P[i]; })), "complement", ctx);
        check(eq(union_of(A, B), expect_pts([&](int i) { return P[i] || Q[i]; })), "union", ctx);
        check(eq(intersection_of(A, B), expect_pts([&](int i) { return P[i] && Q[i]; })), "intersection", ctx);
        check(eq(W, expect_pts([&](int i) {
                  for (int k = i; k < n; ++k) if (in_w(i, k)) return true;
                  return false; })), "window_set", ctx);
        check(eq(eventually_of(A, a, b, W, tr), expect_pts([&](int i) {
                  for (int k = i; k < n; ++k) if (in_w(i, k) && P[k]) return true;
                  return false; })), "eventually", ctx);
        check(eq(always_of(A, a, b, W, tr), expect_pts([&](int i) {
                  for (int k = i; k < n; ++k) if (in_w(i, k) && !P[k]) return false;
                  return true; })), "always", ctx);
        // Until, from the definition: some k in the window with Q[k] and P on
        // [i, k) (strict) or [i, k] (non-strict); "P on [i, l)" by prefix counts.
        std::vector<int> ptrue(n + 1, 0);
        for (int l = 0; l < n; ++l) ptrue[l + 1] = ptrue[l] + P[l];
        auto p_on = [&](int from, int to) { return ptrue[to] - ptrue[from] == to - from; };  // P on [from, to)
        for (bool strict : {true, false}) {
            check(eq(until_of(A, B, a, b, strict, W, tr), expect_pts([&](int i) {
                      for (int k = i; k < n; ++k)
                          if (in_w(i, k) && Q[k] && p_on(i, strict ? k : k + 1)) return true;
                      return false; })), strict ? "until (strict)" : "until (non-strict)", ctx);
        }
        check(eq(eventually_untimed_of(A), expect_pts([&](int i) {
                  for (int k = i; k < n; ++k) if (P[k]) return true;
                  return false; })), "eventually untimed", ctx);
        check(eq(always_untimed_of(A, n), expect_pts([&](int i) {
                  for (int k = i; k < n; ++k) if (!P[k]) return false;
                  return true; })), "always untimed", ctx);
        // Window index queries agree with the predicates.
        for (int j = 0; j < n; ++j) {
            int lo = tr.first_in_window_end(j, b), hi = tr.last_reaching_window_start(j, a);
            int elo = 0;
            while (!tr.within_window_end(elo, j, b)) ++elo;
            int ehi = j;
            while (ehi >= 0 && !tr.reaches_window_start(ehi, j, a)) --ehi;
            check(lo == elo && hi == ehi, "window index queries", ctx);
        }
        ++cases;
    }
    std::printf("  %ld random operator cases\n", cases);
}

// ===========================================================================
// C. SignalConstraintUpdater vs brute force
// ===========================================================================
static void signal_constraint_updater_vs_brute_force() {
    section("C. SignalConstraintUpdater vs brute force");
    Rng R(2002);
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    long steps = 0;
    for (int it = 0; it < 1500; ++it) {
        const int n = it < 1400 ? R.uni(1, 40) : R.uni(200, 1500);   // large n: both rebuild branches
        std::vector<double> sig(n);
        // Keys: 99 -> NaN, 3 -> +inf, -3 -> -inf, otherwise small integers (ties).
        for (auto& v : sig) v = R.coin(0.05) ? 99.0 : R.uni(-3, 3);
        const Trace tr = Trace::make_uniform(1.0, {sig});
        const CompOp op = R.pick(std::vector<CompOp>{CompOp::GT, CompOp::GE, CompOp::LT, CompOp::LE});
        FormulaNode sc = FormulaNode::make_signal_constraint(
            [inf](const std::vector<double>& s) {
                return s[0] == 99.0 ? std::nan("") : s[0] == 3.0 ? inf : s[0] == -3.0 ? -inf : s[0];
            },
            [](const std::vector<double>& x) { return x[0]; }, op, {0});
        // Constraint threshold: 0 (Boolean monitoring) or a robustness-style shift.
        sc.threshold = R.pick(std::vector<double>{0, 0, 0, 1, -1.5, inf, -inf});
        SignalConstraintUpdater up;
        up.init(sc, tr);
        for (int pass = 0; pass < 3; ++pass) {
            up.reset();
            const int t = R.uni(0, n - 1);
            for (int i = n - 1; i >= t; --i) {
                double g;                                    // frozen part
                switch (R.uni(0, 9)) {
                    case 0:  g = nan; break;
                    case 1:  g = R.coin() ? inf : -inf; break;
                    case 2:  g = R.real(-4, 4); break;
                    default: g = R.uni(-4, 4); break;        // ties with keys
                }
                if (n > 100 && R.coin(0.5)) g = sig[i];      // slowly moving frozen value
                up.advance(i, g);
                std::vector<char> e(n, 0);
                for (int l = i; l < n; ++l) {
                    std::vector<double> s = {sig[l]};
                    e[l] = compare(op, sc.signal_fn(s) - g, sc.threshold);
                }
                const IndexIntervals exp = points_to_intervals(e), got = up.intervals();
                bool ok = exp.size() == got.size();
                for (size_t k = 0; ok && k < exp.size(); ++k)
                    ok = exp[k].lo == got[k].lo && exp[k].hi == got[k].hi;
                if (!check(ok, "SignalConstraintUpdater step", "    n=" + std::to_string(n) + " i=" +
                           std::to_string(i) + " op=" + comp_op_to_string(op)))
                    goto next;
                ++steps;
            }
        }
    next:;
    }
    std::printf("  %ld updater steps checked\n", steps);
}

// ===========================================================================
// D. Random STL*_1 formulas: oracle vs all monitors
// ===========================================================================
static void random_formulas_vs_oracle() {
    section("D. Random STL*_1 formulas vs the oracle");
    Rng R(3003);
    long cases = 0, mixed = 0;   // mixed: result neither empty nor the whole trace
    Coverage coverage;
    for (int it = 0; it < 25000; ++it) {
        GenOptions o;
        o.comps = R.uni(1, 3);
        o.fragment_one = true;
        o.max_depth = R.uni(1, 4);
        o.max_freeze_depth    = R.uni(0, 2);
        const OracleFormula F = FormulaGen(R, o).make();
        coverage.add(F);
        const std::string text = to_text(F);
        const Trace tr = random_trace(R, R.uni(1, 11), F.num_signals());
        Formula f;
        try { f = parse_formula(text); }
        catch (const std::exception& e) { check(false, "parse: " + text, e.what()); continue; }
        if (!check(is_fragment_one(f), "generated formula is STL*_1: " + text)) continue;
        const Intervals expected = oracle_intervals(F, tr);
        check_all_monitors("random STL*_1 #" + std::to_string(it), f, tr, expected,
                           "    formula: " + text);
        ++cases;
        if (!expected.empty() && !(expected.size() == 1 && expected[0].lo == tr.timestamp(0) &&
                                   expected[0].hi == tr.timestamp(tr.length - 1)))
            ++mixed;
    }
    std::printf("  %ld random formulas x traces (%ld with a mixed true/false result)\n", cases, mixed);
    coverage.require(Coverage::all_operators());
}

// ===========================================================================
// E. Traces of different sizes and types: reference vs fast
// ===========================================================================
static void traces_of_different_sizes() {
    section("E. Traces of 1 to 5000 samples, noisy and smooth: reference vs fast");
    Rng R(4004);
    long cases = 0;
    for (int it = 0; it < 300; ++it) {
        const int n = random_trace_size(R);
        GenOptions o;
        o.comps = R.uni(1, 2);
        o.max_depth = R.uni(2, 4);
        o.max_freeze_depth = R.uni(1, 2);
        o.allow_until = n <= 500;     // the reference until scans every window
        const std::string text = to_text(FormulaGen(R, o).make());
        const Formula f = parse_formula(text);
        const Trace tr = varied_trace(R, n, f.num_signals);
        const Intervals naive = boolean_monitor_one_naive(f, tr);
        const Intervals fast  = boolean_monitor_one(f, tr);
        // The general monitor re-evaluates nested freezes at every outer step.
        const bool general = n <= 1000;
        const Intervals gfast = general ? boolean_monitor_general(f, tr) : naive;
        check(same(naive, fast) && same(naive, gfast), "size test #" + std::to_string(it),
              "    formula: " + text + "  n=" + std::to_string(n) + "\n    naive " + fmt(naive) +
                  "\n    fast  " + fmt(fast) + "\n    gfast " + fmt(gfast));
        ++cases;
    }
    // Chosen cases on 500 samples (every comparison operator) and 5000 samples.
    for (int n : {500, 5000}) {
        std::vector<double> sig(n);
        for (int i = 0; i < n; ++i) sig[i] = dyadic(std::sin(i * 0.37) * 5 + R.real(-1, 1));
        for (bool uniform : {true, false}) {
            std::vector<double> ts(n);
            for (int i = 0; i < n; ++i) ts[i] = uniform ? i : (i ? ts[i - 1] + R.pick(std::vector<double>{0.25, 0.5, 1, 2}) : 0);
            const Trace tr = uniform ? Trace::make_uniform(1.0, {sig}) : Trace::make_non_uniform(ts, {sig});
            for (const char* op : {">", ">=", "<", "<="}) {
                if (n == 5000 && std::string(op) != ">=") continue;   // one operator at n = 5000
                for (const std::string& fs : {
                         std::string("a*. F[0,40] (a - a* ") + op + " 1)",
                         std::string("G (a*. F[2,30] (a ") + op + " a* + 0.5))",
                         std::string("a*. (a > -4 U[0,5] a ") + op + " a* - 2)"}) {
                    const Formula f = parse_formula(fs);
                    check(same(boolean_monitor_one_naive(f, tr), boolean_monitor_one(f, tr)),
                          "n=" + std::to_string(n) + (uniform ? " uniform " : " non-uniform ") + fs);
                    ++cases;
                }
            }
        }
    }
    std::printf("  %ld cases\n", cases);
}

// ===========================================================================
// F. Scaling
// ===========================================================================
// Average growth factor of the run time when n doubles.  Each measurement
// repeats the run 10 times, since the traces (at most 5000 samples) are fast.
static double growth(const std::function<void(int)>& run, const std::vector<int>& sizes,
                     const char* label) {
    std::vector<double> ms;
    for (int n : sizes) ms.push_back(best_time_ms([&] { for (int r = 0; r < 10; ++r) run(n); }, 5) / 10);
    double prod = 1.0;
    std::printf("  %-48s", label);
    for (size_t k = 0; k < sizes.size(); ++k) std::printf(" n=%d:%.2fms", sizes[k], ms[k]);
    for (size_t k = 1; k < ms.size(); ++k) prod *= ms[k] / std::max(ms[k - 1], 1e-3);
    const double g = std::pow(prod, 1.0 / (ms.size() - 1));
    std::printf("  x%.2f per doubling\n", g);
    return g;
}

static void scaling() {
    section("F. Scaling of the fast STL*_1 monitor (n log n gives ~x2.2 per doubling, n^2 gives x4)");
    const std::vector<int> sizes = {625, 1250, 2500, 5000};

    // Slowly varying signal: few runs per constraint, so O(n log n) overall.
    auto smooth = [](int n, bool uniform) {
        std::vector<double> sig(n), ts(n);
        for (int i = 0; i < n; ++i) {
            sig[i] = std::sin(i * 6.283 * 4 / n) * 10;
            ts[i] = uniform ? i * 0.5 : (i ? ts[i - 1] + ((i % 3) ? 0.25 : 1.0) : 0);
        }
        return uniform ? Trace::make_uniform(0.5, {sig}) : Trace::make_non_uniform(ts, {sig});
    };
    // (No untimed G/F directly above the freeze: those parents let the
    // monitor stop early, which would make the measurement meaningless.)
    const Formula f1 = parse_formula("G[0,10] (a*. F[0,20] (a - a* > 0.5))");
    double g1 = growth([&](int n) { Trace t = smooth(n, true);  boolean_monitor_one(f1, t); },
                       sizes, "G[0,10](a*.F[0,20](a-a*>0.5)), smooth, uniform");
    double g2 = growth([&](int n) { Trace t = smooth(n, false); boolean_monitor_one(f1, t); },
                       sizes, "same, non-uniform");
    const Formula f2 = parse_formula("a*. (a > -11 U[0,15] a - a* >= 3)");
    double g3 = growth([&](int n) { Trace t = smooth(n, true);  boolean_monitor_one(f2, t); },
                       sizes, "a*.(a>-11 U[0,15] a-a*>=3), smooth");

    // Regression test for the flip loop: in the second half of the trace the
    // threshold jumps between 100 and 110 at every step, and all keys of the
    // (still inactive) first half lie in between.  Visiting inactive positions
    // would cost O(n) per step, i.e. O(n^2) in total.
    const Formula f3 = parse_formula("a*. F (a - a* > 100)");
    double g4 = growth([&](int n) {
        std::vector<double> sig(n);
        for (int i = 0; i < n; ++i) sig[i] = i < n / 2 ? 100 + 10.0 * i / n : ((i % 2) ? 10.0 : 0.0);
        Trace t = Trace::make_uniform(1.0, {sig});
        boolean_monitor_one(f3, t);
    }, sizes, "a*.F(a-a*>100), inactive keys between thresholds");

    // A scope part that does not depend on the frozen value (G[0,1] b > 0) is
    // computed once, not at every step.  b alternates, so its input b > 0 has
    // n/2 intervals (computing G[0,1] at every step would cost O(n^2)); the
    // result G[0,1] b > 0 itself is empty.
    const Formula f4 = parse_formula("a*. (G[0,1] b > 0 || F[0,5] (a - a* > 100))");
    double g5 = growth([&](int n) {
        std::vector<double> a(n), b(n);
        for (int i = 0; i < n; ++i) { a[i] = std::sin(i * 0.001); b[i] = (i % 2) ? 1.0 : -1.0; }
        Trace t = Trace::make_uniform(1.0, {a, b});
        boolean_monitor_one(f4, t);
    }, sizes, "a*.(G[0,1] b>0 || F[0,5](a-a*>100)), invariant part");

#ifdef NDEBUG
    for (double g : {g1, g2, g3, g4, g5})
        check(g < 3.3, "fast monitor scales like n log n, not n^2 (growth " + num(g) + " < 3.3)");
#else
    // Debug builds (e.g. _GLIBCXX_DEBUG checks every binary search in linear
    // time) distort the timings: report them only.
    (void)g1; (void)g2; (void)g3; (void)g4; (void)g5;
    std::printf("  (scaling is only asserted in optimized builds)\n");
#endif

    // Noise with long windows: many runs per constraint.
    Rng R(5005);
    const int n = 2000;
    std::vector<double> sig(n);
    for (auto& v : sig) v = dyadic(R.real(0, 100));
    const Trace tr = Trace::make_uniform(1.0, {sig});
    const Formula f = parse_formula("G[0,100] (a*. F[0,400] (a - a* > 0))");
    check(same(boolean_monitor_one_naive(f, tr), boolean_monitor_one(f, tr)),
          "noisy signal with long windows: reference == fast");
}

// ===========================================================================
// G. Large times and magnitudes
// ===========================================================================
// Sample indices covered by a result.
static std::vector<int> indices_of(const Intervals& ivs, const Trace& tr) {
    std::vector<int> out;
    for (int i = 0; i < tr.length; ++i)
        for (const auto& iv : ivs)
            if (iv.lo <= tr.timestamp(i) && tr.timestamp(i) <= iv.hi) { out.push_back(i); break; }
    return out;
}

static void large_times_and_magnitudes() {
    section("G. Large times and magnitudes");
    // The same samples shifted to t ~ 1e9 (Unix-time-like) must give the same
    // samples: window membership may not depend on where the trace sits in time.
    {
        const int n = 60;
        std::vector<double> sig(n), t0(n), t1(n);
        for (int i = 0; i < n; ++i) {
            sig[i] = (i * 7) % 5;
            t0[i] = 0.1 * i;
            t1[i] = 1e9 + 0.1 * i;
        }
        const Trace A = Trace::make_non_uniform(t0, {sig}), B = Trace::make_non_uniform(t1, {sig});
        for (const char* fs : {"F[0.3,0.3] a > 2", "G[0.2,0.5] a > 0", "a*. F[0.3,0.3] (a > a*)",
                               "(a > 0 U[0.3,0.3] a > 2)", "a*. G[0,0.4] (a >= a* - 1)"}) {
            const Formula f = parse_formula(fs);
            const Intervals ra = boolean_monitor_one(f, A), rb = boolean_monitor_one(f, B);
            check(indices_of(ra, A) == indices_of(rb, B) && !indices_of(ra, A).empty() &&
                      same(rb, boolean_monitor_one_naive(f, B)),
                  std::string("time shift to 1e9 does not change the result: ") + fs,
                  "    t=0:   " + fmt(ra) + "\n    t=1e9: " + fmt(rb));
        }
    }
    // A uniform trace with a non-representable dt (offsets are (j - i) * dt)
    // and windows that are exact multiples of it.
    {
        const int n = 5000;
        std::vector<double> sig(n);
        for (int i = 0; i < n; ++i) sig[i] = (i * 7) % 5;
        const Trace tr = Trace::make_uniform(100.1, {sig});
        for (const char* fs : {"F[300.3,300.3] a > 2", "G[100.1,200.2] a > 0", "a*. F[0,300.3] (a > a*)"})
            check_vs_reference(std::string("uniform dt=100.1, n=5000: ") + fs, parse_formula(fs), tr);
    }
    // Large signal values: differences of nearby large numbers.
    {
        const Trace tr = Trace::make_uniform(1.0, {{1e16, 1e16 + 2, 1e16 + 4}});
        for (const char* fs : {"a*. F[1,1] (a - a* >= 2)", "a*. F[1,1] (a - a* > 2)"})
            check_vs_reference(std::string("values near 1e16: ") + fs, parse_formula(fs), tr);
    }
}

int main() {
    std::printf("=== STL*_1 monitors: correctness and complexity ===\n");
    corner_cases();
    interval_operators_vs_pointwise();
    signal_constraint_updater_vs_brute_force();
    random_formulas_vs_oracle();
    traces_of_different_sizes();
    scaling();
    large_times_and_magnitudes();
    return finish("test_boolean_monitor_one");
}
