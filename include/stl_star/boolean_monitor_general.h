#pragma once

#include "stl_star/formula.h"
#include "stl_star/interval_operations.h"
#include "stl_star/trace.h"

namespace stl_star {

// ---------------------------------------------------------------------------
// Fast Boolean monitor for general STL* (signal constraints may read several
// freeze variables, e.g. "a*. G[0,5] b*. F[0,3] (a - a* >= b*)").
//
// Rec-STL*: a FREEZE node runs the backward loop of the STL*_1 monitor.  Every
// nested FREEZE is assumed to depend on the outer frozen values: it is
// re-evaluated (on the suffix starting at i) at every step i of its nearest
// enclosing FREEZE, so the recursion depth is the freeze nesting depth.  A
// FREEZE not nested in another one is evaluated once, and scope nodes built
// only from signal predicates are computed once per evaluation.
//
// Signal constraints use SignalConstraintUpdater; their keys are sorted once for the
// whole run and every recursive call only resets the updater.
//
// Complexity: for freeze nesting depth d, O(n^(d-1)) evaluations of the
// innermost loops, i.e. O(n^d) loop steps in total, each step costing what one
// STL*_1 step costs (about O(log n) when the intervals are few).
// ---------------------------------------------------------------------------
Intervals boolean_monitor_general(const Formula& f, const Trace& trace);

bool boolean_monitor_general_holds_at_start(const Formula& f, const Trace& trace);

}  // namespace stl_star
