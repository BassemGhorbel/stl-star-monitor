#include "stl_star/interval_operations.h"

#include "stl_star/exponential_search.h"

#include <algorithm>

namespace stl_star {

// Append [lo, hi] to a list whose intervals are sorted by lo and whose new
// element does not start before the previous one; merges overlapping or
// adjacent intervals so the list stays canonical.
static void append_merged(IndexIntervals& out, int lo, int hi) {
    if (!out.empty() && lo <= out.back().hi + 1)
        out.back().hi = std::max(out.back().hi, hi);
    else
        out.push_back({lo, hi});
}

// ---------------------------------------------------------------------------
IndexIntervals points_to_intervals(const std::vector<char>& point) {
    IndexIntervals out;
    const int n = static_cast<int>(point.size());
    for (int i = 0; i < n;) {
        if (!point[i]) { ++i; continue; }
        int start = i;
        while (i < n && point[i]) ++i;
        out.push_back({start, i - 1});
    }
    return out;
}

Intervals to_time_intervals(const IndexIntervals& a, const Trace& trace) {
    Intervals out;
    out.reserve(a.size());
    for (const auto& iv : a)
        out.push_back({trace.timestamp(iv.lo), trace.timestamp(iv.hi)});
    return out;
}

bool contains_index(const IndexIntervals& a, int i) {
    // First interval with hi >= i.
    auto it = std::lower_bound(a.begin(), a.end(), i,
                               [](const IndexInterval& iv, int x) { return iv.hi < x; });
    return it != a.end() && it->lo <= i;
}

// ---------------------------------------------------------------------------
IndexIntervals complement_of(const IndexIntervals& a, int n) {
    IndexIntervals out;
    out.reserve(a.size() + 1);
    int next = 0;                         // first index not yet covered
    for (const auto& iv : a) {
        if (iv.lo > next) out.push_back({next, iv.lo - 1});
        next = iv.hi + 1;
    }
    if (next <= n - 1) out.push_back({next, n - 1});
    return out;
}

IndexIntervals union_of(const IndexIntervals& a, const IndexIntervals& b) {
    IndexIntervals out;
    out.reserve(a.size() + b.size());
    size_t i = 0, j = 0;
    while (i < a.size() || j < b.size()) {
        const IndexInterval& nx =
            (j == b.size() || (i < a.size() && a[i].lo <= b[j].lo)) ? a[i++] : b[j++];
        append_merged(out, nx.lo, nx.hi);
    }
    return out;
}

// Intersection of two canonical lists is canonical: two output pieces can only
// be adjacent if both adjacent samples lie in the same interval of `a` and the
// same interval of `b`, in which case they are emitted as one piece.
IndexIntervals intersection_of(const IndexIntervals& a, const IndexIntervals& b) {
    IndexIntervals out;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        int lo = std::max(a[i].lo, b[j].lo);
        int hi = std::min(a[i].hi, b[j].hi);
        if (lo <= hi) out.push_back({lo, hi});
        if (a[i].hi < b[j].hi) ++i; else ++j;
    }
    return out;
}

// Intersection with a window set, cost O(|a| log(1 + |w|/|a|) + output)
// instead of O(|a| + |w|): the window set is fixed for the whole run and may
// have up to n/2 intervals on irregular traces, while `a` is usually much
// shorter.  The first window interval reaching iv.lo is found by exponential
// search from the one found for the previous interval of `a`.
static IndexIntervals intersect_with_window(const IndexIntervals& a,
                                            const IndexIntervals& w) {
    if (w.size() == 1) {                  // common case: one contiguous block
        IndexIntervals out;
        out.reserve(a.size());
        for (const auto& iv : a) {
            int lo = std::max(iv.lo, w[0].lo), hi = std::min(iv.hi, w[0].hi);
            if (lo <= hi) out.push_back({lo, hi});
        }
        return out;
    }
    IndexIntervals out;
    const int m = static_cast<int>(w.size());
    int first = 0;
    for (const auto& iv : a) {
        first = first_true_from(first, m, [&](int k) { return w[k].hi >= iv.lo; });
        for (int k = first; k < m && w[k].lo <= iv.hi; ++k) {
            int lo = std::max(iv.lo, w[k].lo), hi = std::min(iv.hi, w[k].hi);
            if (lo <= hi) out.push_back({lo, hi});
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// window_set: two-pointer sweep.  L(i) = first j with tau_j - tau_i >= a is
// non-decreasing in i; the window of i is non-empty iff that j exists and
// tau_j - tau_i <= b.
// ---------------------------------------------------------------------------
IndexIntervals window_set(const Trace& trace, double a, double b) {
    const int n = trace.length;
    std::vector<char> ok(n, 0);
    int L = 0;
    for (int i = 0; i < n; ++i) {
        L = std::max(L, i);
        while (L < n && !trace.reaches_window_start(i, L, a)) ++L;
        ok[i] = (L < n && trace.within_window_end(i, L, b)) ? 1 : 0;
    }
    return points_to_intervals(ok);
}

// ---------------------------------------------------------------------------
// F[a,b] phi.
//
// For one run [p, q] of phi, sample i is satisfied iff some j in [p, q] lies in
// its window.  Writing L(i) / R(i) for the first / last sample in the window
// (both non-decreasing in i), this is
//     L(i) <= q   and   R(i) >= p   and   L(i) <= R(i)
// The first two conditions give a contiguous range [lo, hi]:
//     lo = first i with tau_p - tau_i <= b,   hi = last i with tau_q - tau_i >= a
// and the third is membership in the window set.  Ranges of successive runs
// have non-decreasing endpoints, so they are merged on the fly, and each
// window lookup starts from the previous answer (exponential search).
// ---------------------------------------------------------------------------
IndexIntervals eventually_of(const IndexIntervals& phi, double a, double b,
                             const IndexIntervals& window, const Trace& trace) {
    IndexIntervals ranges;
    ranges.reserve(phi.size());
    int lo = 0, hi = -1;                  // previous answers: hints for the next lookups
    for (const auto& run : phi) {
        lo = trace.first_in_window_end(run.lo, b, lo);
        hi = trace.last_reaching_window_start(run.hi, a, hi);
        if (lo <= hi) append_merged(ranges, lo, hi);
    }
    return intersect_with_window(ranges, window);
}

IndexIntervals always_of(const IndexIntervals& phi, double a, double b,
                         const IndexIntervals& window, const Trace& trace) {
    const int n = trace.length;
    return complement_of(eventually_of(complement_of(phi, n), a, b, window, trace), n);
}

// ---------------------------------------------------------------------------
// Until.
//
// Take a run I = [p, q] of phi.  A witness j for i in I needs phi on [i, j)
// (strict) or [i, j] (non-strict), so j may range over [p, q+1] (strict) or
// [p, q] (non-strict); j >= i is implied by tau_j - tau_i >= a >= 0.  For each
// run X = [x1, x2] of psi inside that range the satisfied samples are exactly
// those of F[a,b] X restricted to I.  Pieces are produced in order, so they are
// merged on the fly; every (I, X) pair overlaps, and there are at most
// |phi| + |psi| such pairs.  The x1's and the x2's are non-decreasing over
// the pairs (the runs of phi are separated by at least one sample, so the
// pieces of the next run start after jmax), so each window lookup starts from
// the previous answer (exponential search).
//
// Strict until additionally accepts j = i without any phi requirement when
// 0 lies in [a, b]: those samples are exactly psi.
// ---------------------------------------------------------------------------
IndexIntervals until_of(const IndexIntervals& phi, const IndexIntervals& psi,
                        double a, double b, bool strict,
                        const IndexIntervals& window, const Trace& trace) {
    const int n = trace.length;
    IndexIntervals ranges;
    size_t k = 0;
    int first = 0, last = -1;             // previous answers: hints for the next lookups
    for (const auto& I : phi) {
        const int p = I.lo, q = I.hi;
        const int jmax = (strict && q + 1 < n) ? q + 1 : q;
        while (k < psi.size() && psi[k].hi < p) ++k;
        for (size_t m = k; m < psi.size() && psi[m].lo <= jmax; ++m) {
            int x1 = std::max(psi[m].lo, p);
            int x2 = std::min(psi[m].hi, jmax);
            first = trace.first_in_window_end(x1, b, first);
            last  = trace.last_reaching_window_start(x2, a, last);
            int lo = std::max(first, p);
            int hi = std::min(last, q);
            if (lo <= hi) append_merged(ranges, lo, hi);
        }
    }
    IndexIntervals out = intersect_with_window(ranges, window);
    if (strict && trace.reaches_window_start(0, 0, a))   // j = i is allowed (a ~ 0)
        out = union_of(out, psi);
    return out;
}

// ---------------------------------------------------------------------------
IndexIntervals eventually_untimed_of(const IndexIntervals& phi) {
    if (phi.empty()) return {};
    return {{0, phi.back().hi}};
}

IndexIntervals always_untimed_of(const IndexIntervals& phi, int n) {
    if (phi.empty() || phi.back().hi != n - 1) return {};
    return {{phi.back().lo, n - 1}};
}

}  // namespace stl_star
