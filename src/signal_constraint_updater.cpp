#include "stl_star/signal_constraint_updater.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace stl_star {

#if defined(_MSC_VER)
static inline int lowest_bit(std::uint64_t x) {
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return static_cast<int>(idx);
}
#else
static inline int lowest_bit(std::uint64_t x) { return __builtin_ctzll(x); }
#endif

// ===========================================================================
// RankSet
// ===========================================================================
void RankSet::init(int size) {
    size_ = size;
    levels_.clear();
    std::size_t bits = static_cast<std::size_t>(std::max(size, 1));
    do {
        std::size_t words = (bits + 63) / 64;
        levels_.emplace_back(words, 0);
        bits = words;
    } while (bits > 1);
}

void RankSet::clear() {
    for (auto& lv : levels_) std::fill(lv.begin(), lv.end(), 0);
}

void RankSet::insert(int r) {
    std::size_t pos = static_cast<std::size_t>(r);
    for (auto& lv : levels_) {
        lv[pos >> 6] |= std::uint64_t(1) << (pos & 63);
        pos >>= 6;
    }
}

int RankSet::next(int r) const {
    if (r >= size_) return size_;
    // Climb until a set bit at or after the current position is found.
    std::size_t pos = static_cast<std::size_t>(r);
    std::size_t lvl = 0;
    for (;;) {
        if (lvl == levels_.size()) return size_;
        const std::size_t w = pos >> 6;
        if (w >= levels_[lvl].size()) return size_;
        const std::uint64_t bits = levels_[lvl][w] & (~std::uint64_t(0) << (pos & 63));
        if (bits) { pos = (w << 6) + lowest_bit(bits); break; }
        pos = w + 1;   // continue with the next word, one level up
        ++lvl;
    }
    // Descend: bit `pos` at level `lvl` means word `pos` of level lvl-1 is non-zero.
    while (lvl > 0) {
        --lvl;
        pos = (pos << 6) + lowest_bit(levels_[lvl][pos]);
    }
    return static_cast<int>(pos);
}

// ===========================================================================
// SignalConstraintUpdater
// ===========================================================================
void SignalConstraintUpdater::init(const FormulaNode& sc, const Trace& trace) {
    if (sc.type != FormulaType::SIGNAL_CONSTRAINT)
        throw std::invalid_argument("SignalConstraintUpdater::init: not a signal constraint");
    op_        = sc.comp_op;
    threshold_ = sc.threshold;
    n_         = trace.length;

    key_.resize(n_);
    std::vector<double> s(trace.num_components);
    for (int l = 0; l < n_; ++l) {
        trace.sample_at(l, s);
        key_[l] = sc.signal_fn(s);
    }

    // Only finite keys enter the sorted order.  NaN keys are never true;
    // infinite keys (overflow in h) are re-evaluated individually.
    order_.clear();
    for (int l = 0; l < n_; ++l)
        if (std::isfinite(key_[l])) order_.push_back(l);
    std::sort(order_.begin(), order_.end(), [this](int x, int y) {
        return key_[x] < key_[y] || (key_[x] == key_[y] && x < y);
    });
    m_ = static_cast<int>(order_.size());
    sorted_key_.resize(m_);
    rank_.assign(n_, -1);
    for (int r = 0; r < m_; ++r) {
        sorted_key_[r]   = key_[order_[r]];
        rank_[order_[r]] = r;
    }

    active_.init(m_);
    point_.assign(n_, 0);
    lowest_ = n_;
    reset();
}

void SignalConstraintUpdater::reset() {
    for (int l = lowest_; l < n_; ++l) point_[l] = 0;   // only the last pass set anything
    active_.clear();
    active_inf_.clear();
    lowest_   = n_;
    boundary_ = (op_ == CompOp::GT || op_ == CompOp::GE) ? m_ : 0;   // empty true set
    starts_.clear(); finishes_.clear();
    add_s_.clear(); add_f_.clear(); del_s_.clear(); del_f_.clear();
    ivs_.clear();
}

// Boundary rank for frozen part g.  With truth(key) = compare(op, key - g, t):
//   GT/GE: truth is false...false true...true along the sorted keys
//          -> boundary = first true rank
//   LT/LE: truth is true...true false...false -> boundary = first false rank
// The order only holds finite keys, so key - g is non-decreasing in the key
// for every g (all -inf, all +inf or all NaN when g is not finite) and binary
// search with the exact predicate is valid.
int SignalConstraintUpdater::boundary_for(double g) const {
    const CompOp op = op_;
    const double t  = threshold_;
    auto holds = [op, g, t](double key) { return compare(op, key - g, t); };
    auto b = sorted_key_.begin(), e = sorted_key_.end();
    if (op == CompOp::GT || op == CompOp::GE)
        return static_cast<int>(std::partition_point(b, e, [&](double k) { return !holds(k); }) - b);
    return static_cast<int>(std::partition_point(b, e, holds) - b);
}

bool SignalConstraintUpdater::true_at_rank(int r) const {
    return (op_ == CompOp::GT || op_ == CompOp::GE) ? r >= boundary_ : r < boundary_;
}

// Toggle point_[l] and record how the run starts/ends change.  Inactive
// neighbours are false, so they are treated like any false position.
void SignalConstraintUpdater::flip(int l) {
    const bool right = (l + 1 < n_) && point_[l + 1];
    const bool left  = (l > 0) && point_[l - 1];
    if (point_[l]) {                       // true -> false
        if (right) add_s_.push_back(l + 1); else del_f_.push_back(l);
        if (left)  add_f_.push_back(l - 1); else del_s_.push_back(l);
        point_[l] = 0;
    } else {                               // false -> true
        if (right) del_s_.push_back(l + 1); else add_f_.push_back(l);
        if (left)  del_f_.push_back(l - 1); else add_s_.push_back(l);
        point_[l] = 1;
    }
}

void SignalConstraintUpdater::advance(int i, double g) {
    // 1. Flip active positions whose rank lies between the old and new boundary.
    const int nb = boundary_for(g);
    const int lo = std::min(boundary_, nb), hi = std::max(boundary_, nb);
    for (int r = active_.next(lo); r < hi; r = active_.next(r + 1))
        flip(order_[r]);
    boundary_ = nb;
    for (int l : active_inf_)                // infinite keys: direct evaluation
        if (compare(op_, key_[l] - g, threshold_) != (point_[l] != 0)) flip(l);

    // 2. Activate position i with its truth value under the new boundary.
    lowest_ = i;
    if (rank_[i] >= 0) {
        active_.insert(rank_[i]);
        if (true_at_rank(rank_[i])) flip(i);
    } else if (!std::isnan(key_[i])) {
        active_inf_.push_back(i);
        if (compare(op_, key_[i] - g, threshold_)) flip(i);
    }

    // 3. Refresh the run lists.
    rebuild(i);
}

// Merge sorted `base` with sorted additions and removals (multiset semantics;
// every value ends with multiplicity 0 or 1).
static void merge_changes(std::vector<int>& base, std::vector<int>& add,
                          std::vector<int>& del, std::vector<int>& tmp) {
    std::sort(add.begin(), add.end());
    std::sort(del.begin(), del.end());
    tmp.clear();
    tmp.reserve(base.size() + add.size());
    std::size_t a = 0, b = 0, d = 0;
    while (a < base.size() || b < add.size()) {
        int v;
        if (b == add.size() || (a < base.size() && base[a] <= add[b])) v = base[a++];
        else                                                           v = add[b++];
        while (d < del.size() && del[d] < v) ++d;
        if (d < del.size() && del[d] == v) { ++d; continue; }   // cancelled
        tmp.push_back(v);
    }
    base.swap(tmp);
    add.clear();
    del.clear();
}

void SignalConstraintUpdater::rebuild(int i) {
    const std::size_t k = add_s_.size() + add_f_.size() + del_s_.size() + del_f_.size();
    if (k == 0) return;                    // nothing changed: ivs_ is current

    const double sort_cost = static_cast<double>(k) * std::log2(static_cast<double>(k) + 1.0);
    if (sort_cost <= static_cast<double>(n_ - i)) {
        // Sort-and-merge branch: O(m + k log k).
        std::vector<int> tmp;
        merge_changes(starts_,   add_s_, del_s_, tmp);
        merge_changes(finishes_, add_f_, del_f_, tmp);
    } else {
        // Rescan branch: O(n - i).
        add_s_.clear(); add_f_.clear(); del_s_.clear(); del_f_.clear();
        starts_.clear(); finishes_.clear();
        for (int l = i; l < n_;) {
            if (!point_[l]) { ++l; continue; }
            starts_.push_back(l);
            while (l < n_ && point_[l]) ++l;
            finishes_.push_back(l - 1);
        }
    }

    ivs_.resize(starts_.size());
    for (std::size_t r = 0; r < starts_.size(); ++r) ivs_[r] = {starts_[r], finishes_[r]};
}

}  // namespace stl_star
