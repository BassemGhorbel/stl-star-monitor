#pragma once

#include "stl_star/formula.h"
#include "stl_star/interval_operations.h"
#include "stl_star/trace.h"

namespace stl_star {

// ---------------------------------------------------------------------------
// Fast Boolean monitor for STL*_1 (every signal constraint reads only the
// variable of its nearest enclosing freeze operator).
//
// Returns intvl(phi): the samples at which the formula holds, as maximal
// closed timestamp intervals.  Throws std::invalid_argument if the formula is
// not STL*_1 or does not match the trace.
//
// Algorithm.  Nodes are processed bottom-up on interval lists.  A FREEZE node
// x*.phi runs a backward loop i = n-1 .. 0: the frozen value becomes s_x(i),
// every signal constraint of its scope is updated incrementally by a
// SignalConstraintUpdater, the rest of the scope is recomputed from intervals, and
// x*.phi holds at i iff phi holds at i.  Inner FREEZE nodes do not depend on
// the outer frozen value (STL*_1), so they are computed once.
//
// Complexity: O(n log n) to sort each constraint's keys, plus per step i the
// SignalConstraintUpdater cost O(log n + k_i + min(k_i log k_i, n - i) + m_i) and the interval
// operations of the scope, which are linear in the interval counts
// (times log n on non-uniform traces).
// ---------------------------------------------------------------------------
Intervals boolean_monitor_one(const Formula& f, const Trace& trace);

// Does the formula hold at the first sample (time tau_0)?
bool boolean_monitor_one_holds_at_start(const Formula& f, const Trace& trace);

}  // namespace stl_star
