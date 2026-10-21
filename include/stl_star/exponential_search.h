#pragma once

namespace stl_star {

// ---------------------------------------------------------------------------
// Exponential (galloping) search.
//
// first_true_from(lo, hi, pred) returns the smallest i in [lo, hi] with
// pred(i), for a predicate that is false...false true...true on [lo, hi] and
// true at hi (pred(hi) itself is never evaluated).  It probes lo, lo+1, lo+3,
// lo+7, ... and then binary-searches the last gap, so it makes
//     at most 2 * ceil(log2(d + 1)) + 1 calls of pred,   d = answer - lo.
//
// The monitors answer a sequence of such queries whose answers never
// decrease (the queries come from sorted interval lists), starting each one
// at the previous answer.  The distances d_1, d_2, ... then add up to at most
// n, so q queries cost O(sum log(d_k + 1)) = O(q log(1 + n/q)) by concavity of
// log, which is at most O(min(n, q log n)).
// ---------------------------------------------------------------------------
template <class Pred>
int first_true_from(int lo, int hi, Pred pred) {
    if (lo >= hi) return hi;
    if (pred(lo)) return lo;
    int prev = lo;                 // pred(prev) is false
    int step = 1;
    int cur;                       // pred(cur) is true, or cur == hi
    for (;;) {
        cur = (hi - prev > step) ? prev + step : hi;
        if (cur == hi || pred(cur)) break;
        prev = cur;
        step *= 2;
    }
    int l = prev + 1, h = cur;     // answer in [l, h]
    while (l < h) {
        int mid = l + (h - l) / 2;
        if (pred(mid)) h = mid; else l = mid + 1;
    }
    return l;
}

}  // namespace stl_star
