#pragma once

#include "stl_star/constants.h"

#include <vector>

namespace stl_star {

// Maximum number of signal components: one per lowercase letter.  A formula
// names its signals by letters; sorted alphabetically, the k-th letter is
// component k (a formula over x and z reads x = signal[0], z = signal[1]).
static constexpr int MAX_SIGNAL_COMPONENTS = 26;

/**
 * A finite, sampled signal trace.
 *
 *   signal[k][i]  = value of component k (the k-th letter of the formula) at sample i
 *   timestamps[i] = tau_i (strictly increasing; uniform traces use tau_i = i * dt)
 *
 * Length and number of components are inferred from `signal`.
 * Construction validates the input and throws std::invalid_argument on:
 *   empty signal, components of different lengths, more than
 *   MAX_SIGNAL_COMPONENTS components, non-finite values or timestamps,
 *   or consecutive timestamps that are not more than window_tol apart.
 *
 * Time windows.  "Sample j lies in the window [a, b] of sample i" is decided
 * by reaches_window_start / within_window_end below, everywhere in the library
 * (naive monitors, fast monitors, tests), so all implementations agree
 * bit-for-bit.  The offset tau_j - tau_i is computed as (j - i) * dt on
 * uniform traces (it then depends on j - i only) and as a difference of
 * timestamps otherwise; it is compared with the tolerance
 *     window_tol = max(1e-9, 64 * DBL_EPSILON * max(|tau_0|, |tau_{n-1}|))
 * (constants.h), so that rounding errors in the timestamps never decide
 * whether a sample is in a window, also for traces at large absolute times.
 */
struct Trace {
    int    length         = 0;     // number of samples n
    int    num_components = 0;     // number of signal components
    bool   is_uniform     = true;
    double dt             = 0.0;   // sampling period (uniform traces only)
    double window_tol     = 0.0;   // tolerance of window comparisons (see above)

    std::vector<double>              timestamps;  // tau_0 .. tau_{n-1}
    std::vector<std::vector<double>> signal;      // signal[component][sample]

    // Uniform trace: tau_i = i * dt.
    static Trace make_uniform(double dt, std::vector<std::vector<double>> signal);

    // Non-uniform trace with explicit, strictly increasing timestamps.
    static Trace make_non_uniform(std::vector<double> timestamps,
                                  std::vector<std::vector<double>> signal);

    // tau_i and signal value (no bounds checks: 0 <= i < length).
    double timestamp(int i) const { return timestamps[i]; }
    double value(int component, int i) const { return signal[component][i]; }

    // Copy sample i of every component into s (s must have num_components slots).
    void sample_at(int i, std::vector<double>& s) const {
        for (int c = 0; c < num_components; ++c) s[c] = signal[c][i];
    }

    // Time offset tau_j - tau_i (see the class comment).
    double offset(int i, int j) const {
        return is_uniform ? (j - i) * dt : timestamps[j] - timestamps[i];
    }
    // tau_j - tau_i >= a   and   tau_j - tau_i <= b   (with tolerance window_tol;
    // the second is always true for b = +inf).
    bool reaches_window_start(int i, int j, double a) const { return offset(i, j) >= a - window_tol; }
    bool within_window_end(int i, int j, double b) const    { return offset(i, j) <= b + window_tol; }

    // Window index queries used by the fast temporal operators.
    //   first_in_window_end(j, b):        smallest i in [0, j] with within_window_end(i, j, b)
    //   last_reaching_window_start(j, a): largest i in [0, j] with reaches_window_start(i, j, a),
    //                                  or -1 if there is none.
    // Both answers are non-decreasing in j.  `hint` must be a lower bound of
    // the answer (typically the answer for a smaller j; 0, resp. -1, if none).
    // O(1) on uniform traces and for the untimed window [0, inf) (b = inf,
    // resp. a = 0).  On non-uniform traces, exponential search from the hint
    // (exponential_search.h): O(log(answer - hint + 1)), so a sequence of q
    // queries with increasing j, each hinted with the previous answer, costs
    // O(q log(1 + n/q)) = O(min(n, q log n)) in total.
    int first_in_window_end(int j, double b, int hint = 0) const;
    int last_reaching_window_start(int j, double a, int hint = -1) const;

private:
    void validate();
};

}  // namespace stl_star
