#pragma once

#include "stl_star/formula.h"
#include "stl_star/interval_operations.h"
#include "stl_star/trace.h"

#include <cstdint>
#include <vector>

namespace stl_star {

// ---------------------------------------------------------------------------
// RankSet: a set of integers in [0, size) supporting insert and "smallest
// element >= r" queries in O(log_64 size) (at most 4 levels for size <= 16M).  Used to iterate only over ACTIVE positions of the sorted
// key order.
// ---------------------------------------------------------------------------
class RankSet {
public:
    void init(int size);
    void clear();                    // O(size / 64)
    void insert(int r);
    int  next(int r) const;          // smallest element >= r, or size() if none
    int  size() const { return size_; }

private:
    int size_ = 0;
    std::vector<std::vector<std::uint64_t>> levels_;  // levels_[0] = one bit per element
};

// ---------------------------------------------------------------------------
// SignalConstraintUpdater: maintains intvl(c) for one signal constraint
//
//     c:  h(s_l) - g ~ threshold      (g = g(x*), see FormulaNode)
//
// during the backward FREEZE loop i = n-1, n-2, ..., t.  At step i the
// positions l >= i are "active"; position i becomes active at step i, and the
// frozen value (hence g) changes.
//
//   * The keys h(s_l) are computed and sorted ONCE (init, O(n log n)).  Since
//     key - g is monotone in the key, for a fixed g the true positions form a
//     suffix (>, >=) or a prefix (<, <=) of the sorted order (also with
//     infinite keys or g), found by binary search with the exact predicate.
//     Moving g therefore flips exactly the active positions whose rank lies
//     between the old and the new boundary.
//   * Those positions are enumerated through a RankSet of active ranks, so
//     inactive positions are never visited: O(log n + k_i) for k_i flips.
//   * Each flip is applied to the truth array in O(1) and recorded as pending
//     additions/removals of run starts and run ends (lazy updates).
//   * The run lists are then rebuilt either by sorting the k_i pending
//     changes and merging them into the sorted lists, O(m_i + k_i log k_i)
//     with m_i the current number of runs, or by rescanning the truth array,
//     O(n - i); the cheaper option is chosen at every step.
//
// Per step: O(log n + k_i + min(k_i log k_i, n - i) + m_i).
//
// Truth values are decided with compare(op, h(s_l) - g, threshold) on exactly
// the same doubles as the naive monitors, so both always agree.  Only finite
// keys enter the sorted order: NaN keys are never true, and infinite keys
// (possible only when h overflows) are re-evaluated individually at each step.
// ---------------------------------------------------------------------------
class SignalConstraintUpdater {
public:
    // Compute and sort the keys h(s_l) of constraint `sc`.         O(n log n)
    void init(const FormulaNode& sc, const Trace& trace);

    // Start a new backward pass: no position is active.            O(n / 64 + previous pass length)
    void reset();

    // Activate position i (n-1 right after reset(), then i-1, i-2, ...) with
    // g = frozen_fn(x*) for the current frozen values.
    void advance(int i, double g);

    // Runs of true positions among the active ones (canonical, indices >= i).
    const IndexIntervals& intervals() const { return ivs_; }

private:
    int  boundary_for(double g) const;     // rank boundary of the true set
    bool true_at_rank(int r) const;        // is rank r true under boundary_?
    void flip(int l);                      // toggle point_[l], record run changes
    void rebuild(int i);                   // refresh starts_/finishes_/ivs_

    CompOp op_        = CompOp::GE;
    double threshold_ = 0.0;
    int    n_         = 0;
    int    m_         = 0;                // number of finite keys

    std::vector<double> key_;              // key_[l] = h(s_l)
    std::vector<int>    order_;            // positions with a finite key, sorted by key
    std::vector<double> sorted_key_;       // sorted_key_[r] = key_[order_[r]]
    std::vector<int>    rank_;             // rank_[l] = index in order_, -1 if not finite
    std::vector<int>    active_inf_;       // active positions with an infinite key

    RankSet           active_;            // ranks of active positions
    int               boundary_ = 0;      // GT/GE: true iff rank >= boundary_
                                          // LT/LE: true iff rank <  boundary_
    int               lowest_   = 0;      // lowest active position (n_ when none)
    std::vector<char> point_;             // point_[l] = truth value of active position l

    std::vector<int> starts_, finishes_;  // sorted run starts / ends
    std::vector<int> add_s_, add_f_, del_s_, del_f_;   // pending changes of this step
    IndexIntervals   ivs_;
};

}  // namespace stl_star
