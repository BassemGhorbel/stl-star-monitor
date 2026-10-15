#pragma once

#include "stl_star/trace.h"

#include <vector>

namespace stl_star {

// ---------------------------------------------------------------------------
// Result type returned by the monitors: closed intervals of time.
//
// A list of Intervals describes the set of trace samples at which a formula
// holds: every endpoint is an actual trace timestamp, and each interval is a
// maximal run of consecutive satisfying samples.  Lists are sorted.
// ---------------------------------------------------------------------------
struct Interval {
    double lo;
    double hi;
};
using Intervals = std::vector<Interval>;

// ---------------------------------------------------------------------------
// Internal representation: closed intervals of SAMPLE INDICES.
//
// Every interval list used inside the monitors is CANONICAL:
//   sorted, lo <= hi, and separated by at least one missing index
//   (next.lo >= prev.hi + 2), i.e. each interval is a maximal run.
// A canonical list is therefore a unique encoding of a set of samples, which
// makes the fast monitors' output directly comparable with the naive ones.
//
// Complexity notation: |A| is the number of intervals in A, n the trace length.
// ---------------------------------------------------------------------------
struct IndexInterval {
    int lo;
    int hi;
};
using IndexIntervals = std::vector<IndexInterval>;

// Maximal runs of true values in point[0..n-1].                       O(n)
IndexIntervals points_to_intervals(const std::vector<char>& point);

// Convert index intervals to timestamp intervals.                     O(|A|)
Intervals to_time_intervals(const IndexIntervals& a, const Trace& trace);

// Is sample i in A?                                                   O(log |A|)
bool contains_index(const IndexIntervals& a, int i);

// Set operations over the samples {0, ..., n-1}.
IndexIntervals complement_of(const IndexIntervals& a, int n);                 // O(|A|)
IndexIntervals union_of(const IndexIntervals& a, const IndexIntervals& b);    // O(|A|+|B|)
IndexIntervals intersection_of(const IndexIntervals& a, const IndexIntervals& b); // O(|A|+|B|)

// ---------------------------------------------------------------------------
// Window sets.
//
// window_set(trace, a, b) = { i : some sample tau_j satisfies tau_j - tau_i in [a,b] }
// ("the window of i is not empty").  On a non-uniform trace, or when [a,b] is
// narrower than the sampling period, a window can fall entirely between two
// samples; such an i can never satisfy F[a,b] and always satisfies G[a,b].
// The set depends only on the trace and [a,b], so the monitors compute it once
// per temporal operator (O(n)) and reuse it in every evaluation.
// ---------------------------------------------------------------------------
IndexIntervals window_set(const Trace& trace, double a, double b);

// ---------------------------------------------------------------------------
// Temporal operators (discrete-time STL semantics over the trace samples).
//
// `window` must be window_set(trace, a, b).
// k = |output|; the log n factor appears on non-uniform traces with timed
// windows only (the window lookups are O(1) for a = 0 and for b = inf).
// ---------------------------------------------------------------------------

// F[a,b] phi : exists j, tau_j - tau_i in [a,b], phi(j).
//   O(|phi| log n + |phi| log |window| + k)
IndexIntervals eventually_of(const IndexIntervals& phi, double a, double b,
                             const IndexIntervals& window, const Trace& trace);

// G[a,b] phi : for all j with tau_j - tau_i in [a,b], phi(j) (vacuously true
// when the window is empty).  Computed as  not F[a,b] not phi.
IndexIntervals always_of(const IndexIntervals& phi, double a, double b,
                         const IndexIntervals& window, const Trace& trace);

// phi U[a,b] psi (strict):     exists j >= i, tau_j - tau_i in [a,b], psi(j),
//                              phi(k) for all i <= k < j.
// phi UNS[a,b] psi (non-strict): same but phi(k) for all i <= k <= j.
//   O((|phi| + |psi|) log n + (|phi| + |psi|) log |window| + k)
IndexIntervals until_of(const IndexIntervals& phi, const IndexIntervals& psi,
                        double a, double b, bool strict,
                        const IndexIntervals& window, const Trace& trace);

// Untimed F phi and G phi (window [0, inf)).                          O(1)
IndexIntervals eventually_untimed_of(const IndexIntervals& phi);
IndexIntervals always_untimed_of(const IndexIntervals& phi, int n);

}  // namespace stl_star
