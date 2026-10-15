#pragma once
// ---------------------------------------------------------------------------
// Test utilities shared by every test program: pass/fail reporting, printing
// of interval lists and traces, exact comparison of results, timing.
// ---------------------------------------------------------------------------
#include "stl_star/boolean_monitor_general.h"
#include "stl_star/boolean_monitor_one.h"
#include "stl_star/formula_parser.h"
#include "stl_star/interval_operations.h"
#include "stl_star/naive_monitors.h"
#include "stl_star/robustness_monitor_general.h"
#include "stl_star/robustness_monitor_one.h"
#include "stl_star/trace.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace stl_star;   // the tests use the library names unqualified

namespace test_util {


// ===========================================================================
// Reporting
// ===========================================================================
struct Stats { long pass = 0, fail = 0; };
inline Stats& stats() { static Stats s; return s; }

inline bool check(bool ok, const std::string& label, const std::string& detail = "") {
    if (ok) { ++stats().pass; return true; }
    ++stats().fail;
    if (stats().fail <= 40) {
        std::printf("  FAIL %s\n", label.c_str());
        if (!detail.empty()) std::printf("%s\n", detail.c_str());
    }
    return false;
}

inline void section(const char* title) { std::printf("\n--- %s ---\n", title); std::fflush(stdout); }

inline int finish(const char* name) {
    std::printf("\n=== %s: %ld passed, %ld failed ===\n", name, stats().pass, stats().fail);
    return stats().fail == 0 ? 0 : 1;
}

inline std::string num(double v) {
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char b[40]; std::snprintf(b, sizeof b, "%.17g", v); return b;
}
inline std::string short_num(double v) {
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char b[40]; std::snprintf(b, sizeof b, "%g", v); return b;
}

inline std::string fmt(const Intervals& ivs) {
    std::string s = "{";
    for (size_t i = 0; i < ivs.size(); ++i)
        s += (i ? "," : "") + std::string("[") + short_num(ivs[i].lo) + "," + short_num(ivs[i].hi) + "]";
    return s + "}";
}

// Endpoints are exact trace timestamps, so lists can be compared exactly.
inline bool same(const Intervals& a, const Intervals& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].lo != b[i].lo || a[i].hi != b[i].hi) return false;
    return true;
}

inline std::string fmt_trace(const Trace& tr) {
    std::string s = "    trace: n=" + std::to_string(tr.length) +
                    (tr.is_uniform ? " dt=" + short_num(tr.dt) : std::string(" non-uniform")) + "\n    t=[";
    for (int i = 0; i < tr.length; ++i) s += (i ? "," : "") + short_num(tr.timestamp(i));
    s += "]";
    for (int c = 0; c < tr.num_components; ++c) {
        s += "\n    " + std::string(1, char('a' + c)) + "=[";
        for (int i = 0; i < tr.length; ++i) s += (i ? "," : "") + short_num(tr.value(c, i));
        s += "]";
    }
    return s;
}

// Largest number of FREEZE nodes on a root-to-leaf path.
inline int freeze_depth(const Formula& f) {
    std::vector<int> depth(f.nodes.size(), 0);   // children come before parents
    for (size_t j = 0; j < f.nodes.size(); ++j) {
        int d = 0;
        for (int c : f.nodes[j].children) d = std::max(d, depth[c]);
        depth[j] = d + (f.nodes[j].type == FormulaType::FREEZE ? 1 : 0);
    }
    return f.nodes.empty() ? 0 : depth[f.root_index()];
}

// Intervals of the samples where pred(i) holds.
inline Intervals intervals_where(const Trace& tr, const std::function<bool(int)>& pred) {
    std::vector<char> pt(tr.length);
    for (int i = 0; i < tr.length; ++i) pt[i] = pred(i) ? 1 : 0;
    return to_time_intervals(points_to_intervals(pt), tr);
}

// ===========================================================================
// Timing helper: best of `reps` runs, in milliseconds.
// ===========================================================================
template <class Fn>
double best_time_ms(Fn fn, int reps = 3) {
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        fn();
        auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

}  // namespace test_util
