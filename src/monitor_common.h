#pragma once
// Internal helpers shared by the monitor implementations (not public API).

#include "stl_star/formula.h"
#include "stl_star/interval_operations.h"
#include "stl_star/trace.h"

#include <vector>

namespace stl_star {

namespace monitor_detail {

// parent[i] = index of the parent of node i (-1 for the root).
std::vector<int> parents_of(const Formula& f);

// owner[i] = nearest FREEZE node strictly above node i, or -1.
std::vector<int> freeze_owners(const Formula& f);

// Validate the formula and check it against the trace (a parsed formula needs
// exactly one trace signal per letter; freeze components must exist in the
// trace).  With require_fragment_one, also require the STL*_1 fragment.
// Throws std::invalid_argument.
void check_input(const Formula& f, const Trace& trace, bool require_fragment_one);

// Window sets (see interval_operations.h) for every timed temporal node;
// empty vectors for all other nodes.
std::vector<IndexIntervals> compute_window_sets(const Formula& f, const Trace& trace);

// Truth of a signal predicate at sample s.
inline bool predicate_holds(const FormulaNode& nd, const std::vector<double>& s) {
    return compare(nd.comp_op, nd.signal_fn(s), nd.threshold);
}
// Truth of a signal constraint given h = signal_fn(s) and g = frozen_fn(x*).
// The reference monitors decide constraints with this function, and
// SignalConstraintUpdater evaluates the same expression compare(op, h - g, threshold).
inline bool constraint_holds(const FormulaNode& nd, double h, double g) {
    return compare(nd.comp_op, h - g, nd.threshold);
}

// Throws std::invalid_argument if the formula contains a negation.
void check_negation_free(const Formula& f);

// Intervals of a SIGNAL_PRED node.                                      O(n)
IndexIntervals signal_pred_intervals(const FormulaNode& nd, const Trace& trace);

// intvl of a Boolean / temporal node from its children's intervals.
IndexIntervals compute_node(int j, const Formula& f, const Trace& trace,
                            const std::vector<IndexIntervals>& windows,
                            const std::vector<IndexIntervals>& iv);

// ---------------------------------------------------------------------------
// Acceleration heuristics for the fast monitors.  They exploit the parent of
// a FREEZE node; they are only valid because formulas are trees (the FREEZE
// result is used by its parent alone).
//
//   AND_SIBLING: parent is (phi1 && x*.phi) or (x*.phi && phi1): when phi1 is
//              false at i the parent is false at i, so the instantiation of
//              x* at i is skipped.
//   OR_SIBLING:  parent is (phi1 || x*.phi): when phi1 is true at i, same.
//   ALWAYS_UNTIMED:     parent is G(x*.phi): once x*.phi is false at i,
//              G is false at every j <= i, so the loop can stop.
//   EVENTUALLY_UNTIMED: parent is F(x*.phi): once true at i, F is true at
//              every j <= i.
//
// AND_SIBLING / OR_SIBLING need the intervals of phi1 before the FREEZE loop
// runs; each monitor checks that and otherwise drops the heuristic (NONE).
// ---------------------------------------------------------------------------
enum class FreezeCtxKind { NONE, AND_SIBLING, OR_SIBLING, ALWAYS_UNTIMED, EVENTUALLY_UNTIMED };

struct FreezeCtx {
    FreezeCtxKind kind = FreezeCtxKind::NONE;
    int sibling = -1;   // phi1 for AND_SIBLING / OR_SIBLING
};

FreezeCtx analyze_freeze_ctx(int freeze_idx, const Formula& f, const std::vector<int>& parent);

// Membership test for a decreasing sequence of queries i: `cursor` starts at
// ivs.size()-1 and only moves left, so a whole backward pass costs O(n + |ivs|).
inline bool contains_decreasing(const IndexIntervals& ivs, int i, int& cursor) {
    while (cursor >= 0 && ivs[cursor].lo > i) --cursor;
    return cursor >= 0 && ivs[cursor].hi >= i;
}

}  // namespace monitor_detail

}  // namespace stl_star
