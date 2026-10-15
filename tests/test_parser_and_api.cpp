// test_parser_and_api.cpp
//
//   A. Parser: accepted syntax (all keyword forms, numbers, leading
//      parentheses, windows, freeze binders) checked by evaluating the result.
//   B. Parser: rejected inputs.
//   C. Signal letters <-> components, trace validation.
//   D. Low-level API: building formulas in code, Formula::validate().
//   E. One-call API (stl_star_monitor.h), including every README example,
//      compared with the reference monitors.
// Exit code is non-zero on any failure.
#include "stl_star/stl_star_monitor.h"
#include "test_utils.h"

#include <cmath>

using namespace test_util;

// A test trace with one signal per letter of the formula.
static Trace test_trace(const Formula& f) {
    const Signal rows = {{3, 7, 2, 8, 5, 1, 9, 4, 6, 2}, {1, 1, 2, 2, 3, 3, 4, 4, 5, 5}};
    return Trace::make_uniform(1.0, Signal(rows.begin(), rows.begin() + f.num_signals));
}

// Two formulas must denote the same set of samples on the test trace.
static void equivalent(const char* text, const char* reference) {
    std::string err;
    try {
        const Formula f = parse_formula(text), g = parse_formula(reference);
        const Intervals got = boolean_monitor_general(f, test_trace(f));
        const Intervals exp = boolean_monitor_general(g, test_trace(f));
        check(same(got, exp), std::string("parse \"") + text + "\" == \"" + reference + "\"",
              "    got " + fmt(got) + "\n    exp " + fmt(exp));
    } catch (const std::exception& e) {
        check(false, std::string("parse \"") + text + "\"", std::string("    threw: ") + e.what());
    }
}

static void rejected(const char* text, const char* fragment) {
    try {
        parse_formula(text);
        check(false, std::string("reject \"") + text + "\"", "    was accepted");
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        check(msg.find(fragment) != std::string::npos,
              std::string("reject \"") + text + "\" with a message containing \"" + fragment + "\"",
              "    message: " + msg);
    }
}

// ===========================================================================
static void accepted_syntax() {
    section("A. Accepted syntax");
    // Keyword forms (case-insensitive) and symbols.
    equivalent("Always[0,2] a >= 3",            "G[0,2] a >= 3");
    equivalent("BOX[0,2] a >= 3",               "G[0,2] a >= 3");
    equivalent("eventually a >= 9",             "F a >= 9");
    equivalent("Diamond[1,2] a >= 9",           "F[1,2] a >= 9");
    equivalent("not a >= 5",                    "! a >= 5");
    equivalent("(a >= 3 and a <= 7)",           "(a >= 3 && a <= 7)");
    equivalent("(a < 3 OR a > 7)",              "(a < 3 || a > 7)");
    equivalent("(a > 2 Until[0,3] a >= 8)",     "(a > 2 U[0,3] a >= 8)");
    equivalent("(a > 2 untilNS[0,3] a >= 8)",   "(a > 2 UNS[0,3] a >= 8)");
    equivalent("(a > 2 uns[0,3] a >= 8)",       "(a > 2 UNS[0,3] a >= 8)");
    // Untimed forms and their alternatives: window [0, inf).
    equivalent("Always a >= 2",                 "G[0,inf] a >= 2");
    equivalent("Box a >= 2",                    "G a >= 2");
    equivalent("Eventually a >= 9",             "F[0,inf] a >= 9");
    equivalent("diamond a >= 9",                "F a >= 9");
    equivalent("(a > 2 U a >= 8)",              "(a > 2 U[0,inf] a >= 8)");
    equivalent("(a > 2 Until a >= 8)",          "(a > 2 U a >= 8)");
    equivalent("(a > 2 UNS a >= 8)",            "(a > 2 UNS[0,inf] a >= 8)");
    equivalent("(a > 2 UntilNS a >= 8)",        "(a > 2 UNS a >= 8)");
    equivalent("F[2,inf] a >= 9",               "F[2,1e9] a >= 9");
    // Numbers.
    equivalent("a >= 1e1",                      "a >= 10");
    equivalent("a >= 2.5E-1",                   "a >= 0.25");
    equivalent("a > 1e-400",                    "a > 0");      // underflows to 0
    // Arithmetic, precedence, unary minus, leading parentheses.
    equivalent("(a + 1) * 2 >= 10",             "a >= 4");
    equivalent("((a)) >= 4",                    "a >= 4");
    equivalent("-a <= -4",                      "a >= 4");
    equivalent("2 * a - b / 2 > 3",             "4 * a - b > 6");
    equivalent("(a - b) * (0 - 1) < 0",         "a > b");
    equivalent("((a >= 3) && (a <= 7))",        "(a >= 3 && a <= 7)");
    // Freeze binders, shadowing, subscripts, constraints written in any order.
    equivalent("a*. F (a* < a)",                "a*. F (a - a* > 0)");
    equivalent("a*. F (2 * (a - a*) > 1)",      "a*. F (a - a* > 0.5)");
    equivalent("a*. F ((a + a*) / 2 > a*)",     "a*. F (a > a*)");
    equivalent("a*. F ((a - a*) / 4 > 0.25)",   "a*. F (a - a* > 1)");
    equivalent("a*. F (a - 3 * a* + 2 * a* > 0)", "a*. F (a - a* > 0)");
    equivalent("b*. G[0,2] (a > b* + b)",       "b*. G[0,2] (a - b > b*)");
    equivalent("a*. G (a*. F (a > a*))",        "G (a*. F (a > a*))");          // inner binder shadows
    equivalent("a*_1. G (a*_2. F (a > a*_2))",  "G (a*. F (a > a*))");
    equivalent("a*. (a*_7. F (a > a*_7) && a >= a*)", "(a*. F (a > a*) && a >= a)");
    // Products / quotients of frozen values only (non-linear in a*) are fine.
    equivalent("a*. F (a > a* * a* / 8)",       "a*. F (8 * a > a* * a*)");
    // "a*" without a space is always the freeze variable: a*-1 is (a*) - 1.
    equivalent("a*. F (a > a*-1)",              "a*. F (a - a* > -1)");
    equivalent("a*. F (a * 2 > a*)",            "a*. F (2 * a - a* > 0)");
    equivalent("a*.F(a>a*)",                    "a*. F (a > a*)");             // spaces optional
    // Negation and nested binary operators.
    equivalent("! ! a >= 5",                    "a >= 5");
    equivalent("((a >= 3 && a <= 7) || a > 8)", "! ((a < 3 || a > 7) && a <= 8)");
    equivalent("a*. ((a > 2 U[0,3] a >= 8) && F a > a*)", "a*. (F a > a* && (a > 2 U[0,3] a >= 8))");

    // README "accepted" table: classification STL*_1 vs general STL*.
    const std::vector<std::pair<const char*, bool>> table = {
        {"G[0,5] a*. F[0,3] (a - a* >= 0.5)",         true},
        {"a*. G (b*. F (b - b* > 0) && a >= a*)",     true},
        {"a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1)",  false},
        {"a*. G (b*. F (b > a*))",                    false},
        {"a*. F (a* > 2)",                            true},
        {"a*. F (a > a*-1)",                          true},
    };
    for (const auto& [text, fragment_one] : table) {
        try {
            const Formula f = parse_formula(text);
            const Trace tr = test_trace(f);
            check(is_fragment_one(f) == fragment_one, std::string("README table: \"") + text + "\" is " +
                                               (fragment_one ? "STL*_1" : "general STL*"));
            bool threw = false;
            try { boolean_monitor_one(f, tr); } catch (const std::invalid_argument&) { threw = true; }
            check(threw == !fragment_one, std::string("README table: monitor_one() ") +
                                         (fragment_one ? "accepts" : "rejects") + " \"" + text + "\"");
            if (fragment_one)
                check(same(boolean_monitor_one(f, tr), boolean_monitor_general(f, tr)),
                      std::string("README table: both monitors agree on \"") + text + "\"");
        } catch (const std::exception& e) {
            check(false, std::string("README table: parse \"") + text + "\"", e.what());
        }
    }
}

// ===========================================================================
static void rejected_input() {
    section("B. Rejected input");
    rejected("a*. a * a* > 1",              "may not combine signals with freeze variables");
    rejected("a*. (a - 1) * a* > 1",        "may not combine signals with freeze variables");
    rejected("a*. F (a / a* > 1)",          "may not combine signals with freeze variables");
    rejected("a*. F (a * (a* + 1) > 0)",    "may not combine signals with freeze variables");
    rejected("a > a*",                      "not bound");
    rejected("a*. F (a - b* > 0)",          "not bound");
    rejected("a*_1. a*_1. F (a > a*_1)",    "already bound");
    rejected("A > 0",                       "Unknown operator");
    rejected("ab > 0",                      "Unknown identifier");
    rejected("(a > 2 UN[0,3] a >= 8)",      "Unknown identifier");   // the non-strict until is UNS
    rejected("F[2,1] a > 0",                "Invalid window");
    rejected("F[-1,2] a > 0",               "window start");
    rejected("F[0,2 a > 0",                 "Expected ']'");
    rejected("(a > 0",                      "Expected ')'");
    rejected("a > 0 && a < 1",              "Unexpected token");
    rejected("a >",                         "Expected a number");
    rejected("a ? 3",                       "Unexpected character");
    rejected("F[0,2]",                      "Expected a number");
    rejected(".5 < a",                      "Expected a number");
    rejected("a*_. F (a > a*_)",            "Expected digits");
    rejected("a > 1e999",                   "Number out of range");
    rejected("F[0,1e999] a > 0",            "Number out of range");
    rejected("a*_12345678901. F (a > a*_12345678901)", "Subscript too long");
    rejected("a*. F (a*2 > 1)",              "to multiply");
    rejected("a*. F (a*b > 1)",              "to multiply");
    rejected("a*. F (a*(b + 1) > 1)",        "to multiply");
    rejected("a*2 > 1",                      "to multiply");
    rejected("a*_1. F (a > a*)",             "not bound");      // a* and a*_1 are different names
    rejected("(a > 0 && b > 0 && a > b)",    "Expected ')'");   // one binary operator per pair
    rejected("a*. F (a > a*) && a > 1",      "Unexpected token");
    rejected("G[0,2) a > 0",                 "Expected ']'");
    rejected("a = 0",                        "Unexpected character");
    rejected("a != 0",                       "Unexpected character");
    rejected("a > 0 && b > 0",               "Unexpected token");
    rejected("G[0,2] A > 0",                 "Unknown operator");
    rejected("a*. F ((a + 1) * (a* - 1) > 0)", "may not combine signals with freeze variables");
}

// ===========================================================================
static void signals_and_traces() {
    section("C. Signals and traces");
    // The letters of the formula, sorted, are the components 0, 1, 2, ...
    auto names = [](const char* text) { return parse_formula(text).signal_names; };
    check(names("a > 0") == "a", "signal_names(a) == \"a\"");
    check(names("(c > 0 && a > 0)") == "ac", "signal_names(c, a) == \"ac\"");
    check(names("(z > 0 && x > 0)") == "xz", "signal_names(z, x) == \"xz\"");
    check(names("b*. F (a > b*)") == "ab", "freeze letters count as signals");
    check(names("(x > a && F[0,2] (x > 1 || a < 0))") == "ax", "repeated letters count once");
    check(parse_formula("(z > 0 && x > 0)").num_signals == 2, "num_signals(z, x) == 2");

    // Any letters: the result equals the same formula written over a, b, ...
    // on the same signals.
    const Signal s2 = {{3, 7, 2, 8, 5, 1, 9, 4}, {1, 1, 2, 2, 3, 3, 4, 4}};
    const Trace t2 = Trace::make_uniform(1.0, s2);
    for (const auto& [text, renamed] : std::vector<std::pair<const char*, const char*>>{
             {"(z - x > 1 && F[0,2] x > 6)",               "(b - a > 1 && F[0,2] a > 6)"},
             {"(x > 2 && a < 3)",                           "(b > 2 && a < 3)"},
             {"x*. F[0,3] (z - x > -x* - 4 && g > x*)",     "b*. F[0,3] (c - b > -b* - 4 && a > b*)"},
             {"u*. G (f*. F (u - u* + f - f* > 2))",        "b*. G (a*. F (b - b* + a - a* > 2))"}}) {
        const Formula f = parse_formula(text), g = parse_formula(renamed);
        const Trace tr = f.num_signals == 2 ? t2 : Trace::make_uniform(1.0, {s2[0], s2[1], {0, 5, 1, 6, 2, 7, 3, 8}});
        const Intervals reference = boolean_monitor_general_naive(g, tr);
        check(same(boolean_monitor_general(f, tr), reference) &&
                  same(monitor_general(text, tr.signal, 1.0).intervals, reference),
              std::string("letters in alphabetical order: \"") + text + "\" == \"" + renamed + "\"");
    }
    {   // all 26 letters
        std::string text = "a";
        for (char c = 'b'; c <= 'z'; ++c) text += std::string(" + ") + c;
        text = "z*. F (" + text + " > z* + 40)";
        Signal s(26);
        for (int k = 0; k < 26; ++k) s[k] = {double(k % 3), double(k % 5), double(k % 2), double(k % 7)};
        const Formula f = parse_formula(text);
        const Trace tr = Trace::make_uniform(1.0, s);
        check(f.num_signals == 26 && same(boolean_monitor_one(f, tr), boolean_monitor_one_naive(f, tr)),
              "26 signals a..z");
    }

    auto throws_invalid = [](const std::function<void()>& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        return false;
    };
    auto message = [](const std::function<void()>& fn) {
        try { fn(); } catch (const std::invalid_argument& e) { return std::string(e.what()); }
        return std::string();
    };
    const Trace one   = Trace::make_uniform(1.0, {{1, 2}});
    const Trace three = Trace::make_uniform(1.0, {{1, 2}, {3, 4}, {5, 6}});
    check(throws_invalid([&] { boolean_monitor_one(parse_formula("(a > 0 && x > 0)"), one); }),
          "fewer signals than letters is rejected");
    check(throws_invalid([&] { boolean_monitor_one(parse_formula("(a > 0 && x > 0)"), three); }),
          "more signals than letters is rejected");
    check(message([&] { boolean_monitor_general(parse_formula("x*. F (a > x*)"), one); }).find("(a, x)") !=
              std::string::npos,
          "the error names the letters in order");
    check(throws_invalid([&] { boolean_monitor_one_naive(parse_formula("(a > 0 && x > 0)"), three); }),
          "reference monitor rejects a signal count mismatch");
    check(throws_invalid([&] { robustness_monitor_one(parse_formula("(a > 0 && x > 0)"), one, 0.01); }),
          "robustness rejects a signal count mismatch");
    check(throws_invalid([&] { robustness_monitor_general(parse_formula("x*. F (a > x*)"), three, 0.01); }),
          "general robustness rejects a signal count mismatch");

    check(throws_invalid([] { Trace::make_uniform(1.0, {}); }), "no components");
    check(throws_invalid([] { Trace::make_uniform(1.0, {{}}); }), "no samples");
    check(throws_invalid([] { Trace::make_uniform(1.0, {{1, 2}, {1}}); }), "components of different lengths");
    check(throws_invalid([] { Trace::make_uniform(0.0, {{1, 2}}); }), "dt = 0");
    check(throws_invalid([] { Trace::make_uniform(-1.0, {{1, 2}}); }), "dt < 0");
    check(throws_invalid([] { Trace::make_uniform(1.0, {{1, std::nan("")}}); }), "NaN value");
    check(throws_invalid([] { Trace::make_uniform(1.0, {{1, INFINITY}}); }), "infinite value");
    check(throws_invalid([] { Trace::make_non_uniform({0, 0}, {{1, 2}}); }), "equal timestamps");
    check(throws_invalid([] { Trace::make_non_uniform({1, 0}, {{1, 2}}); }), "decreasing timestamps");
    check(throws_invalid([] { Trace::make_non_uniform({0, 1e-12}, {{1, 2}}); }), "timestamps closer than the tolerance");
    check(throws_invalid([] { Trace::make_non_uniform({0}, {{1, 2}}); }), "timestamp count mismatch");
    check(throws_invalid([] { Trace::make_uniform(1.0, Signal(27, std::vector<double>{1})); }),
          "more than 26 components");
}

// ===========================================================================
static void low_level_api() {
    section("D. Low-level API");
    // G[0,0.5] (a*. sin(a) >= a* - 1), built in code.
    Formula f;
    f.num_freeze_vars = 1;
    int c = f.add_node(FormulaNode::make_signal_constraint(
        [](const std::vector<double>& s) { return std::sin(s[0]); },
        [](const std::vector<double>& x) { return x[0] - 1.0; }, CompOp::GE, {0}));
    int fr = f.add_node(FormulaNode::make_freeze(/*var=*/0, /*component=*/0, c));
    f.add_node(FormulaNode::make_always(0.0, 0.5, fr));
    const Trace tr = Trace::make_uniform(0.5, {{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}});
    const Intervals reference = boolean_monitor_one_naive(f, tr);
    check(same(boolean_monitor_one(f, tr), reference) && same(boolean_monitor_general(f, tr), reference),
          "programmatic formula with a custom signal function", fmt(reference));

    auto invalid = [](Formula g) {
        try { g.validate(); } catch (const std::invalid_argument&) { return true; }
        return false;
    };
    {   // shared child (a DAG, not a tree)
        Formula g;
        int p = g.add_node(FormulaNode::make_signal_pred([](const std::vector<double>& s) { return s[0]; }, CompOp::GT));
        g.add_node(FormulaNode::make_and(p, p));
        check(invalid(g), "validate rejects shared children");
    }
    {   // constraint outside the scope of its variable
        Formula g;
        g.num_freeze_vars = 1;
        int sc = g.add_node(FormulaNode::make_signal_constraint(
            [](const std::vector<double>& s) { return s[0]; },
            [](const std::vector<double>& x) { return x[0]; }, CompOp::GT, {0}));
        int p  = g.add_node(FormulaNode::make_signal_pred([](const std::vector<double>& s) { return s[0]; }, CompOp::GT));
        int fz = g.add_node(FormulaNode::make_freeze(0, 0, p));
        g.add_node(FormulaNode::make_and(sc, fz));
        check(invalid(g), "validate rejects a constraint outside its freeze scope");
    }
    {   // child index not smaller than the parent
        Formula g;
        g.add_node(FormulaNode::make_not(0));
        check(invalid(g), "validate rejects bad child indices");
    }
    {   // invalid window
        Formula g;
        int p = g.add_node(FormulaNode::make_signal_pred([](const std::vector<double>& s) { return s[0]; }, CompOp::GT));
        g.add_node(FormulaNode::make_eventually(2.0, 1.0, p));
        check(invalid(g), "validate rejects a window with b < a");
    }
}

// ===========================================================================
// A one-call result must equal the reference monitor on the same trace
// (whose first sample is at t = 0), including holds_at_start.
static void check_one_call(const std::string& label, const MonitorResult& r, const Intervals& reference) {
    const bool start = !reference.empty() && reference.front().lo == 0.0;
    check(same(r.intervals, reference) && r.holds_at_start == start,
          label, "    got " + fmt(r.intervals) + "  reference " + fmt(reference));
}

static void one_call_api() {
    section("E. One-call API and README examples vs the reference monitors");
    {   // No signal map: letters follow the order of signal[].
        const char* f = "a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)";
        const Signal s = {{0.0, 1.0, 2.0, 3.0, 4.0}, {0.5, 0.5, 0.5, 0.5, 0.5}};
        check_one_call("monitor_general example", monitor_general(f, s, 1.0),
                       boolean_monitor_general_naive(parse_formula(f), Trace::make_uniform(1.0, s)));
    }
    {
        const char* f = "G[0,3] a*. F[0,2] (a - a* >= 0.5)";
        const Signal s = {{0.0, 0.5, 1.2, 2.0, 2.8, 3.5}};
        check_one_call("monitor_one example", monitor_one(f, s, 1.0),
                       boolean_monitor_one_naive(parse_formula(f), Trace::make_uniform(1.0, s)));
    }
    {
        const char* f = "a*. G[0,4] (a - a* >= 0)";
        const Signal s = {{1.0, 1.5, 2.0, 2.5, 3.0}};
        const std::vector<double> t = {0.0, 0.5, 1.2, 2.8, 4.0};
        check_one_call("monitor_one non-uniform example", monitor_one(f, s, t),
                       boolean_monitor_one_naive(parse_formula(f), Trace::make_non_uniform(t, s)));
    }
    {
        const char* f = "G[0,3] a*. F[0,2] (a - a* >= 0.5)";
        const Signal s = {{0.0, 0.6, 1.3, 2.1, 3.0}};
        const std::vector<double> t = {0.0, 0.4, 1.0, 2.2, 3.5};
        check_one_call("monitor_one non-uniform example 2", monitor_one(f, s, t),
                       boolean_monitor_one_naive(parse_formula(f), Trace::make_non_uniform(t, s)));
    }
    {
        const char* f = "a >= 5";
        const Signal s = {{1, 3, 5, 7, 9}};
        const MonitorResult r = monitor_one(f, s, 1.0);
        check(!r.holds_at_start, "holds_at_start is false when the first sample fails");
        check_one_call("monitor_one, formula false at the start", r,
                       boolean_monitor_one_naive(parse_formula(f), Trace::make_uniform(1.0, s)));
    }
    {   // Robustness: the bracket must contain the reference value.
        const char* f = "G[0,3] a*. F[0,2] (a - a* >= 0)";
        const Signal s = {{0.0, 0.5, 1.2, 2.0, 2.8, 3.5}};
        const double exact = robustness_monitor_one_naive(parse_formula(f), Trace::make_uniform(1.0, s));
        const RobustnessResult r = robustness_one(f, s, 1.0);
        check(r.lo <= exact && exact <= r.hi && r.hi - r.lo <= 0.01 * std::abs(exact) + 1e-9,
              "robustness_one example", "[" + num(r.lo) + ", " + num(r.hi) + "] reference " + num(exact));
        const RobustnessResult g = robustness_general(f, s, 1.0);
        check(g.lo == r.lo && g.hi == r.hi, "robustness_general agrees");
        const RobustnessResult nu = robustness_one(f, s, std::vector<double>{0, 1, 2, 3, 4, 5});
        check(nu.lo == r.lo && nu.hi == r.hi, "non-uniform overload with unit spacing agrees");
    }
    {
        const char* f = "a*. b*. F (a + b - a* - b* > 0)";
        const Signal s = {{1, 2, 0}, {0, 2, 1}};
        const double exact = robustness_monitor_general_naive(parse_formula(f), Trace::make_uniform(1.0, s));
        const RobustnessResult r = robustness_general(f, s, 1.0, 0.001);
        check(r.lo <= exact && exact <= r.hi, "robustness_general example",
              "[" + num(r.lo) + ", " + num(r.hi) + "] reference " + num(exact));
    }
    {   // monitor_one() refuses non-STL*_1 formulas with a helpful message.
        bool ok = false;
        try { monitor_one("a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)", {{0, 1}, {0, 1}}, 1.0); }
        catch (const std::invalid_argument& e) {
            ok = std::string(e.what()).find("general") != std::string::npos;
        }
        check(ok, "monitor_one() points to monitor_general for general formulas");
    }
}

int main() {
    std::printf("=== Parser, traces and API ===\n");
    accepted_syntax();
    rejected_input();
    signals_and_traces();
    low_level_api();
    one_call_api();
    return finish("test_parser_and_api");
}
