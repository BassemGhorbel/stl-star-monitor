#pragma once

#include "stl_star/formula.h"

#include <string>

namespace stl_star {

// ---------------------------------------------------------------------------
// Parse a text STL* formula.
//
// Signals
//   Any lowercase letter names a signal.  The letters used in the formula,
//   in alphabetical order, are components 0, 1, 2, ...: a formula over x and
//   z reads x from signal[0] and z from signal[1]; one over a and x reads a
//   from signal[0] and x from signal[1].  The trace must have exactly one
//   signal per letter (Formula::signal_names lists the letters in order).
//
// Atomic propositions:  expr  op  expr      with op in  >=  >  <=  <
//   expr: numbers (2, 0.5, 1e-3), signals (a, b), freeze variables (a*, a*_2),
//         + - * / and unary minus, parentheses.
//   A product or quotient may not mix signals with freeze variables
//   ("a - 2*a* > 0" is fine, "a * a* > 1" is rejected).
//
// Operators (keywords are case-insensitive):
//   G[a,b] phi   Always / Box            F[a,b] phi   Eventually / Diamond
//   G phi, F phi (untimed, window [0, inf)), with the same alternatives
//   (phi U[a,b] psi)    Until    strict until:     phi on [i, j), psi at j
//   (phi UNS[a,b] psi)  UntilNS  non-strict until: phi on [i, j], psi at j
//   (phi U psi), (phi UNS psi)   untimed until, window [0, inf)
//   ! phi / not phi,  (phi && psi) / and,  (phi || psi) / or
//   Binary operators must be parenthesised.  The upper bound may be "inf".
//
// Freeze binders:
//   a*. phi     binds a* to the current value of signal 'a'
//   a*_2. phi   subscripted name: use it to keep two values of the same signal
//               (a*_1 and a*_2) alive at once.  Re-binding a plain a* inside
//               another a*. scope shadows the outer binding.
//
// Examples:
//   "a >= 2"
//   "G[0,5] a*. F[0,3] (a - a* >= 0.5)"
//   "a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)"
//
// Throws std::runtime_error (with the position in the input) on any error.
// ---------------------------------------------------------------------------
Formula parse_formula(const std::string& input);

}  // namespace stl_star
