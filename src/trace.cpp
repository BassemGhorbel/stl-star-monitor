#include "stl_star/trace.h"

#include "stl_star/exponential_search.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace stl_star {

Trace Trace::make_uniform(double dt, std::vector<std::vector<double>> signal) {
    Trace t;
    t.is_uniform     = true;
    t.dt             = dt;
    t.signal         = std::move(signal);
    t.num_components = static_cast<int>(t.signal.size());
    t.length         = t.num_components > 0 ? static_cast<int>(t.signal[0].size()) : 0;
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("Sampling period dt must be finite and positive");
    t.timestamps.resize(t.length);
    for (int i = 0; i < t.length; ++i) t.timestamps[i] = i * dt;
    t.validate();
    return t;
}

Trace Trace::make_non_uniform(std::vector<double> timestamps,
                              std::vector<std::vector<double>> signal) {
    Trace t;
    t.is_uniform     = false;
    t.dt             = 0.0;
    t.timestamps     = std::move(timestamps);
    t.signal         = std::move(signal);
    t.num_components = static_cast<int>(t.signal.size());
    t.length         = t.num_components > 0 ? static_cast<int>(t.signal[0].size()) : 0;
    t.validate();
    return t;
}

void Trace::validate() {
    if (num_components < 1 || num_components > MAX_SIGNAL_COMPONENTS)
        throw std::invalid_argument(
            "The signal must have between 1 and " +
            std::to_string(MAX_SIGNAL_COMPONENTS) + " components (got " +
            std::to_string(num_components) + ")");
    if (length <= 0)
        throw std::invalid_argument("The trace must contain at least one sample");
    for (int k = 0; k < num_components; ++k) {
        if (static_cast<int>(signal[k].size()) != length)
            throw std::invalid_argument("All signal components must have the same length");
        for (double v : signal[k])
            if (!std::isfinite(v)) throw std::invalid_argument("Signal values must be finite");
    }
    if (static_cast<int>(timestamps.size()) != length)
        throw std::invalid_argument("The number of timestamps must equal the signal length");
    for (double t : timestamps)
        if (!std::isfinite(t)) throw std::invalid_argument("Timestamps must be finite");

    const double scale = std::max(std::abs(timestamps.front()), std::abs(timestamps.back()));
    window_tol = std::max(TIMESTAMP_EPS, TIMESTAMP_ULPS * std::numeric_limits<double>::epsilon() * scale);
    // Consecutive samples must be further apart than the tolerance, otherwise
    // a sample could fall in a window that starts after it.
    for (int i = 1; i < length; ++i)
        if (!(offset(i - 1, i) > window_tol))
            throw std::invalid_argument(
                "Timestamps must be strictly increasing, by more than the window "
                "tolerance (" + std::to_string(window_tol) + ")");
}

// ---------------------------------------------------------------------------
// Window index queries.
//
// Both predicates are monotone in i for fixed j (the offset decreases as i
// grows and floating-point arithmetic is monotone), so the answer is a
// boundary, found by exponential search from the hint (a lower bound of the
// answer).  On uniform traces an arithmetic guess is corrected by a constant
// number of predicate checks, which gives the same boundary.
// ---------------------------------------------------------------------------
static int clamp_guess(double g, int lo, int hi) {
    if (!(g > lo)) return lo;     // also catches NaN / -inf
    if (g >= hi)   return hi;
    return static_cast<int>(g);
}

int Trace::first_in_window_end(int j, double b, int hint) const {
    auto ok = [&](int i) { return within_window_end(i, j, b); };
    if (std::isinf(b)) return 0;
    if (is_uniform) {
        int i = clamp_guess(j - std::floor(b / dt), 0, j);
        while (i > 0 && ok(i - 1)) --i;
        while (i < j && !ok(i)) ++i;
        return i;
    }
    // ok is false...false true...true on [hint, j], and ok(j) is always true.
    return first_true_from(std::min(std::max(hint, 0), j), j, ok);
}

int Trace::last_reaching_window_start(int j, double a, int hint) const {
    auto ok = [&](int i) { return reaches_window_start(i, j, a); };
    if (a <= window_tol) return j;  // ok(j): offset 0 reaches the start (untimed windows)
    if (is_uniform) {
        int i = clamp_guess(j - std::ceil(a / dt), -1, j);
        while (i < j && ok(i + 1)) ++i;
        while (i >= 0 && !ok(i)) --i;
        return i;
    }
    // ok is true...true false...false on [0, j] and ok(j) is false (a > tol),
    // so the answer is one before the first i in [hint + 1, j] with !ok(i).
    const int from = std::min(std::max(hint, -1), j - 1) + 1;
    return first_true_from(from, j, [&](int i) { return !ok(i); }) - 1;
}

}  // namespace stl_star
