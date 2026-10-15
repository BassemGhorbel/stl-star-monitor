# STL\* Monitor

A C++17 library for monitoring **Signal Temporal Logic with the freeze operator (STL\*)** on finite, sampled traces.

Given a formula and a trace, the library answers two questions:

- **Boolean monitoring:** at which sample times does the formula hold?
- **Robustness monitoring:** by how much does the formula hold or fail at the first sample?

Both questions are answered by our algorithms. Each monitor also has a simple **reference (naive) implementation** written directly from the semantics. The reference monitors are used only for testing: the tests compare our monitors against them. Calling one of our monitors never runs a reference monitor.

- **Boolean:** Our monitors and reference monitors return identical results.
- **Robustness:** Our monitors approximate the robustness to a chosen precision `epsilon`. The reference monitors compute it exactly, and the exact value always lies inside our monitor's bracket.

---

## Contents

1. [Quick start](#quick-start)
2. [Build and test](#build-and-test)
3. [Inputs and outputs](#inputs-and-outputs)
4. [Writing formulas](#writing-formulas)
5. [Semantics](#semantics)
6. [Robustness](#robustness)
7. [Errors](#errors)
8. [Low-level API](#low-level-api)
9. [How the monitors work](#how-the-monitors-work)
10. [Project layout](#project-layout)
11. [How correctness is tested](#how-correctness-is-tested)
12. [License](#license)

---

## Quick start

```cpp
#include "stl_star/stl_star_monitor.h"
#include <iostream>

int main() {
    // Freeze a now (a*). At every time j within the next 3 time units, freeze
    // b (b*); within 2 time units after j, a must have risen above a* by at
    // least b* - 1.
    auto result = stl_star::monitor_general(
        "a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1.0)",
        {
            {0.0, 1.0, 2.0, 3.0, 4.0},   // signal a
            {0.5, 0.5, 0.5, 0.5, 0.5}    // signal b
        },
        1.0);                            // sampling period dt

    std::cout << (result.holds_at_start ? "holds" : "fails") << " at t = 0\n";
    for (const stl_star::Interval& iv : result.intervals)
        std::cout << "[" << iv.lo << ", " << iv.hi << "]\n";   // prints [0, 4]
}
```

Signals are named by lowercase letters. Give one signal per letter of the formula, in alphabetical order of the letters. Here the letters are `a` and `b`, so `signal[0]` is `a` and `signal[1]` is `b`. A formula over `x` and `z` takes `x` as `signal[0]` and `z` as `signal[1]`. There is no signal map to pass.

Every name of the library is in the namespace `stl_star`. The other examples in this README assume `using namespace stl_star;`.

The one-call API in `stl_star/stl_star_monitor.h` has four functions:

| Function | Use it for |
|---|---|
| `monitor_one(formula, signal, dt_or_timestamps)` | Boolean monitoring of **STL\*₁** formulas |
| `monitor_general(formula, signal, dt_or_timestamps)` | Boolean monitoring of **any** STL\* formula |
| `robustness_one(formula, signal, dt_or_timestamps, epsilon = 0.01)` | Robustness of STL\*₁ formulas |
| `robustness_general(formula, signal, dt_or_timestamps, epsilon = 0.01)` | Robustness of any STL\* formula |

**Which monitor to use:**

- **STL\*₁ formula → `monitor_one()` / `robustness_one()`.** In STL\*₁, every constraint reads only the variable of its nearest enclosing freeze operator, as in `G[0,5] a*. F[0,3] (a - a* >= 0.5)`. Freeze operators may be nested, but an inner one never reads an outer frozen value, so the scopes are independent.
- **Any other STL\* formula → `monitor_general()` / `robustness_general()`.** This includes nested freeze operators whose constraints read outer frozen values, as in the example above, where the constraint reads both `a*` and `b*`. The general monitor treats every nested freeze operator as dependent on the ones around it. It re-evaluates the nested operator at every step of the enclosing one.

`monitor_one()` and `robustness_one()` reject a formula outside STL\*₁, and the error message says to use the general version. The general functions accept STL\*₁ formulas too and give the same result. They are slower when STL\*₁ freezes are nested, because they re-evaluate the inner freeze at every outer step.

---

## Build and test

Requirements: CMake 3.14 or newer and a C++17 compiler (GCC, Clang or MSVC).

```bash
mkdir build && cd build
cmake ..                # Release by default
cmake --build .
./examples              # usage examples of every monitor
ctest                   # runs the 5 test programs (about 30 seconds)
```

The examples and tests are optional. Both are built by default. To build only the library, pass these two options to the `cmake ..` command above:

```bash
cmake .. -DSTL_STAR_BUILD_EXAMPLES=OFF -DSTL_STAR_BUILD_TESTS=OFF
```

The options are defined in [`CMakeLists.txt`](CMakeLists.txt).

To use the library, link against the `stl_star` target and include `stl_star/stl_star_monitor.h`. [`examples/examples.cpp`](examples/examples.cpp) shows every function in use: both Boolean monitors, non-uniform traces, robustness, error handling and the low-level API.

---

## Inputs and outputs

**Signal.** `signal[k][i]` is the value of signal `k` at sample `i`. Signal `k` is the `k`-th letter of the formula in alphabetical order, and there must be exactly one signal per letter (up to 26). All signals must have the same length. Values must be finite.

```cpp
// letters a and x: a = signal[0], x = signal[1]
monitor_one("x*. F[0,2] (x - x* > a)",
            {{0.5, 0.5, 1.0, 1.0},     // a
             {0.0, 1.0, 2.0, 2.5}},    // x
            1.0);
```

**Time.** Pass either of these:

- a sampling period `dt`, so sample `i` is at time `i * dt`; or
- a `std::vector<double>` of timestamps, which must be strictly increasing.

```cpp
monitor_one("a*. G[0,4] (a - a* >= 0)",
            {{1.0, 1.5, 2.0, 2.5, 3.0}},
            std::vector<double>{0.0, 0.5, 1.2, 2.8, 4.0});   // non-uniform trace
```

**Boolean result** (`MonitorResult`):

```cpp
struct MonitorResult {
    bool      holds_at_start;  // does the formula hold at the first sample?
    Intervals intervals;       // the sample times at which it holds
};
struct Interval { double lo, hi; };   // Intervals = std::vector<Interval>
```

`intervals` lists maximal runs of consecutive satisfying samples. Every endpoint is an actual timestamp of the trace. For example, `{[0, 1], [4, 4]}` on a `dt = 1` trace means the formula holds at samples 0, 1 and 4.

**Robustness result** (`RobustnessResult`):

```cpp
struct RobustnessResult {
    double lo, hi;       // the robustness at the first sample lies in [lo, hi]
    int monitor_calls;   // Boolean monitor calls used by the bisection
};
```

---

## Writing formulas

### Grammar

```
formula := x*. formula                       freeze (x is a signal letter)
         | G formula | G[lo,hi] formula      always
         | F formula | F[lo,hi] formula      eventually
         | ! formula
         | ( formula && formula )
         | ( formula || formula )
         | ( formula U formula )   | ( formula U[lo,hi] formula )     strict until
         | ( formula UNS formula ) | ( formula UNS[lo,hi] formula )   non-strict until
         | ( formula )
         | expr op expr                      atomic proposition, op is >=  >  <=  <

expr    := number | signal | freeze-var | -expr | expr (+ - * /) expr | ( expr )
```

- **Signals** are single lowercase letters, any of `a` to `z`. The letters used in the formula, in alphabetical order, are `signal[0]`, `signal[1]`, …: in `(x > 2 && a < 1)`, `a` is `signal[0]` and `x` is `signal[1]`. A letter used several times counts once, and a frozen letter (`x` in `x*.`) counts as a signal. The trace must have exactly as many signals as the formula has letters.
- **Freeze variables** are a letter followed by `*`: `a*`, or with a subscript `a*_1`, `a*_2`. `a*.` binds `a*`, and `a*_2.` binds `a*_2`.
- **Numbers** are written as `2`, `0.5`, `1e-3` or `2.5E+2`. A leading dot such as `.5` is not accepted, so write `0.5`.
- **Windows** are `[lo,hi]` with `0 <= lo <= hi`. `hi` may be `inf`.
- `*` and `/` bind tighter than `+` and `-`. Spaces are optional.

### Atomic propositions

```
a >= 2
2 * a - b / 2 > 3
(a + 1) * 2 >= 10
a - a* >= 0.5            constraint: uses a frozen value
a - a* >= b* - 1.0       constraint over two frozen values (general STL*)
```

A proposition that reads a freeze variable is a **constraint**. A constraint must split into a signal part and a frozen part, `f(signals) ~ g(frozen values)`. In practice, a product or quotient may not combine a signal with a freeze variable:

```
a - 2 * a* > 1          accepted
a* * a* < a + 1         accepted (the frozen part may be non-linear)
(a - a*) * 2 > 1        accepted (multiplying by a constant is fine)
(a - a*) / 4 > 1        accepted
a * a* > 1              rejected
a / a* > 2              rejected
(a + 1) * (a* - 1) > 0  rejected
```

The algorithm depends on this split, and the parser reports a clear error when it is not possible. The constraint is evaluated as `f(signals) − g(frozen values) ~ 0`. With extreme magnitudes this can differ from evaluating the literal expression: `(a - a*) * 1e308 > 0` overflows on both sides.

**`a*` always means the freeze variable.** There is no space-free way to write "a times …" after a letter:

- `a*-1` is read as `a* - 1`;
- `a*2`, `a*b` and `a*(b+1)` are rejected with a hint to write `a * 2`.

Put a space between a signal and `*` when you mean multiplication.

### Operators

Each binary operator (`&&`, `||`, `U`, `UNS`) must be wrapped in its own pair of parentheses, so `(p && q && r)` must be written `((p && q) && r)`. Unary operators (`!`, `G`, `F`, freeze) apply to the formula that follows them. Keywords are case-insensitive. The one-letter operators `G`, `F` and `U` must be uppercase, because lowercase `g`, `f` and `u` are signals.

| Syntax | Alternatives | Meaning |
|---|---|---|
| `! phi` | `not phi` | negation |
| `(phi && psi)` | `(phi and psi)` | conjunction |
| `(phi \|\| psi)` | `(phi or psi)` | disjunction |
| `G[a,b] phi` | `Always[a,b]`, `Box[a,b]` | phi at every sample whose time offset is in [a, b] |
| `F[a,b] phi` | `Eventually[a,b]`, `Diamond[a,b]` | phi at some sample whose time offset is in [a, b] |
| `G phi` | `Always phi`, `Box phi` | untimed always, window [0, ∞) |
| `F phi` | `Eventually phi`, `Diamond phi` | untimed eventually, window [0, ∞) |
| `(phi U[a,b] psi)` | `(phi Until[a,b] psi)` | strict until: phi on [i, j), psi at j |
| `(phi UNS[a,b] psi)` | `(phi UntilNS[a,b] psi)` | non-strict until: phi on [i, j], psi at j |
| `(phi U psi)` | `(phi Until psi)` | untimed strict until, window [0, ∞) |
| `(phi UNS psi)` | `(phi UntilNS psi)` | untimed non-strict until, window [0, ∞) |
| `x*. phi` | | freeze: bind `x*` to the current value of signal `x` |

Windows are closed and need `0 <= a <= b`. The upper bound may be `inf`, as in `F[2,inf] phi`.

### Freeze variables

`a*. phi` binds `a*` to the value of `a` at the time where `a*.` is evaluated. Inside `phi`, `a*` refers to that value.

- **Plain names.** A plain name such as `a*` refers to the nearest enclosing `a*.` binder. A nested `a*.` shadows the outer one.
- **Subscripted names.** Use `a*_1`, `a*_2`, … to keep two frozen values of the same signal alive at once, for example `a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))`. `a*` and `a*_1` are different names, so `a*_1. F (a > a*)` is rejected. A subscripted binder cannot be re-bound inside its own scope, so `a*_1. G a*_1. …` is rejected. The same subscript may be reused in a separate, non-nested scope.
- **Different signals.** Different letters can be frozen together, as in `a*. b*. F (a + b - a* - b* > 0)`.
- **Scope.** A freeze variable may only be used inside the scope of its binder. Otherwise the parser reports `not bound`.

### Accepted and rejected formulas at a glance

| Formula | Result |
|---|---|
| `G[0,5] a*. F[0,3] (a - a* >= 0.5)` | accepted, STL\*₁ |
| `a*. G (b*. F (b - b* > 0) && a >= a*)` | accepted, STL\*₁ (nested but independent) |
| `a*. G[0,3] b*. F[0,2] (a - a* >= b* - 1)` | accepted, general STL\* (`monitor_one()` rejects it, use `monitor_general()`) |
| `a*. G (b*. F (b > a*))` | accepted, general STL\* (the constraint reads the outer `a*`) |
| `a*. F (a* > 2)` | accepted (a constraint may read only frozen values) |
| `a*. F (a > a*-1)` | accepted, means `a > a* - 1` |
| `a > a*` | rejected: `a*` is not bound |
| `a*. F (a * a* > 1)` | rejected: product of a signal and a freeze variable |
| `a*. F (a*2 > 1)` | rejected: write `a * 2` |
| `a > 0 && b > 0` | rejected: binary operators need parentheses |
| `(a > 0 && b > 0 && a > b)` | rejected: one binary operator per pair of parentheses |
| `F[2,1] a > 0`, `F[-1,2] a > 0` | rejected: invalid window |
| `a = 0`, `a != 0` | rejected: only `>= > <= <` are supported |
| `.5 < a` | rejected: write `0.5` |
| `a > 1e999` | rejected: number out of range |
| `G[0,2] A > 0`, `ab > 0` | rejected: signals are single lowercase letters |

### More examples

```
G[0,5] a*. F[0,3] (a - a* >= 0.5)              a always rises by 0.5 within 3 time units
a*. G (a - a* >= 0)                            a never drops below its current value
(a >= 5 && a*. F (a - a* > 0))                 a >= 5 now and a rises later
a*. (a >= 3 U[0,3] a - a* >= 4)                a stays >= 3 until it has risen by 4
a*_1. F (a*_2. F (a - a*_1 - a*_2 > 0))        some later value exceeds the sum of two earlier ones
```

---

## Semantics

Formulas are evaluated at the **samples** of the trace (discrete-time semantics). Write `t_i` for the time of sample `i` and `s_i` for its values. The following holds at sample `i`:

| Formula | Holds at sample i when |
|---|---|
| `h ~ g` | `h(s_i) − g(frozen values) ~ 0` |
| `F[a,b] phi` | some `j` with `t_j - t_i` in `[a, b]` satisfies phi |
| `G[a,b] phi` | every `j` with `t_j - t_i` in `[a, b]` satisfies phi (vacuously true if there is none) |
| `(phi U[a,b] psi)` | some `j` with `t_j - t_i` in `[a, b]` satisfies psi, and phi holds at `i … j-1` |
| `(phi UNS[a,b] psi)` | the same, but phi must also hold at `j` |
| `x*. phi` | phi holds at `i` with `x*` = value of `x` at sample `i` |

Some consequences are worth knowing:

- **Values between samples are not interpolated.** A window can fall between two samples, for example `F[0.5,1]` on a trace sampled every 1.5 time units. At such a time, `F` is false and `G` is vacuously true.
- **Windows near the end of the trace are cut off.** A window that extends past the last sample only contains the samples that exist. `F[1,2] phi` is false at the last sample, and `G[1,2] phi` is true there.

---

## Robustness

Robustness `rho` measures how strongly the formula holds at the first sample:

- `rho(h > g)` and `rho(h >= g)` are `h - g`; `rho(h < g)` and `rho(h <= g)` are `g - h`;
- `&&` and `||` take the minimum and maximum;
- `F` and `G` take the maximum and minimum over their window. An empty window gives `-inf` for `F` and `+inf` for `G`;
- until takes the maximum over `j` of the minimum of psi at `j` and phi on the samples it must cover.

If `rho > 0`, the formula holds. If `rho < 0`, it fails.

```cpp
auto r = robustness_one("G[0,3] a*. F[0,2] (a - a* >= 0)",
                        {{0.0, 0.5, 1.2, 2.0, 2.8, 3.5}}, 1.0, 0.01);
// r.lo <= 1.2 <= r.hi
```

**How it is computed.** Robustness is found by bisection. Each step asks the Boolean monitor whether `rho >= r`, using the shifted formula in which `h ~ g` becomes `h − g >= r` or `h − g <= −r`. Each atom compares the same difference that defines its robustness, so every decision is exact, even for values around 1e16.

1. Bounds `[lo, hi]` covering every possible atomic value are computed. For a constraint over several frozen variables, every combination of trace values is tried, so the bounds are valid for any expression. Bounds that overflow are clamped to the largest finite double.
2. Two calls detect infinite robustness. In that case `lo = hi = -inf` or `lo = hi = +inf`. For example, `F[5,10] a > 0` on a 3-second trace has `rho = -inf`.
3. Bisection runs until `hi - lo <= epsilon * max(|mid|, 1% of the initial width)`. So **`epsilon` is the error relative to the robustness value itself**, with an absolute floor for values near 0.

Robustness requires a **negation-free** formula. Formulas containing `!` or `not` are rejected.

---

## Errors

| Problem | Exception |
|---|---|
| Syntax error, unbound or re-bound freeze variable, non-separable constraint, invalid window | `std::runtime_error`, with the position in the formula |
| Invalid trace: empty, different signal lengths, non-finite values, `dt <= 0`, timestamps not increasing | `std::invalid_argument` |
| The trace does not have exactly one signal per letter of the formula | `std::invalid_argument` |
| `monitor_one` or `robustness_one` called with a non-STL\*₁ formula | `std::invalid_argument`, saying to use the general version |
| Robustness of a formula with negation, or `epsilon` not in (0, 1] | `std::invalid_argument` |

---

## Low-level API

The one-call functions are thin wrappers. You can also build the pieces yourself:

```cpp
#include "stl_star/formula_parser.h"              // Formula parse_formula(const std::string&)
#include "stl_star/trace.h"                       // Trace::make_uniform / Trace::make_non_uniform
#include "stl_star/boolean_monitor_one.h"         // our STL*1 monitor
#include "stl_star/boolean_monitor_general.h"     // our general monitor
#include "stl_star/robustness_monitor_one.h"      // RobustnessResult, compute_robustness_range
#include "stl_star/robustness_monitor_general.h"
#include "stl_star/naive_monitors.h"              // the four reference monitors
using namespace stl_star;

Formula   f  = parse_formula("a*. F[0,2] (a - a* > 1)");   // parse once, reuse on many traces
Trace     tr = Trace::make_uniform(0.5, {{0, 1, 3, 2, 4}});
Intervals r  = boolean_monitor_one(f, tr);                  // same as boolean_monitor_one_naive(f, tr)
bool      ok = boolean_monitor_one_holds_at_start(f, tr);
RobustnessResult rho = robustness_monitor_one(f, tr, /*epsilon=*/0.01);   // bracket [lo, hi]
double    exact      = robustness_monitor_one_naive(f, tr);                 // lies in [rho.lo, rho.hi]
```

`f.signal_names` lists the letters of a parsed formula in the order the trace must give them, for example `"ax"` for a formula over `x` and `a`.

Formulas can also be built in code. This is useful for functions the parser does not support, such as `sin`. A constraint is given as a signal part `h(s)` and a frozen part `g(x*)`, compared as `h(s) ~ g(x*)`:

```cpp
Formula f;
f.num_freeze_vars = 1;
int c  = f.add_node(FormulaNode::make_signal_constraint(
             [](const std::vector<double>& s) { return std::sin(s[0]); },   // h(s)
             [](const std::vector<double>& x) { return x[0] - 1.0; },       // g(x*)
             CompOp::GE, /*freeze variables read by g=*/{0}));
int fr = f.add_node(FormulaNode::make_freeze(/*variable=*/0, /*signal=*/0, c));
f.add_node(FormulaNode::make_always(0.0, 0.5, fr));   // root = last node
f.validate();   // checks the tree structure and variable scopes
```

Nodes are stored children-first, and the last node is the root. The formula must be a tree: no node may be shared. In a formula built in code, `s[k]` is `signal[k]`. There are no letters to check, so the trace only needs to contain every signal that the functions read.

---

## How the monitors work

**Our monitors.** A satisfaction set is a list of intervals of sample indices. Boolean and temporal operators combine these lists directly. A freeze operator `x*. phi` scans the trace backwards. At each sample `i` it sets `x*` to the value of `x` at `i` and updates the scope `phi`.

- **Signal constraints.** These are updated incrementally by `SignalConstraintUpdater`. The signal part of each sample is sorted once. A new frozen value only moves a boundary through that order, so only the samples whose truth value actually changes are visited.
- **STL\*₁ monitor.** The scopes of nested freeze operators are independent, so each freeze operator is evaluated once.
- **General monitor.** Nested freeze operators are assumed to depend on the freeze operators around them. Each one is re-evaluated at every step of its enclosing freeze operator.

**Acceleration heuristics.** The freeze loop runs backwards and skips work that cannot change the result. This depends on the operator directly above the freeze:

- **Conjunctive freeze** `(phi1 && x*. phi2)`: where `phi1` is false at sample `i`, the conjunction is false, so `x*` is not instantiated at `i`.
- **Disjunctive freeze** `(phi1 || x*. phi2)`: where `phi1` is true at `i`, the disjunction is true, so `x*` is not instantiated at `i`.
- **Box with freeze** `G x*. phi`: once `x*. phi` is false at `i`, `G` is false at every earlier sample, so the loop stops.
- **Diamond with freeze** `F x*. phi`: once `x*. phi` is true at `i`, `F` is true at every earlier sample, so the loop stops.

`phi1` may be any formula whose value is known before the freeze loop runs.

**Our robustness monitors.** These find the robustness by bisection, calling our Boolean monitor on shifted formulas (see [Robustness](#robustness)).

**Reference monitors.** These evaluate the semantics point by point, with no intervals: truth values for the Boolean monitors and robustness values for the robustness monitors. A freeze operator re-evaluates its scope for every sample.

- **Timed `F` and `G`, Boolean:** prefix counts give the number of true samples in each window.
- **Timed `F` and `G`, robustness:** the sliding max/min of each window is computed with Lemire's streaming filter.
- **Until:** checks the window of every sample.

---

## Project layout

```
include/stl_star/                       public headers, everything in namespace stl_star
  stl_star_monitor.h                    one-call API (start here)
  formula_parser.h                      parse_formula: text -> Formula
  formula.h                             FormulaNode, Formula, is_fragment_one
  trace.h                               Trace and the time-window predicates
  interval_operations.h                 Interval (results) and the interval operators
  signal_constraint_updater.h           SignalConstraintUpdater (incremental constraint update)
  boolean_monitor_one.h                 STL*1 Boolean monitor
  boolean_monitor_general.h             general STL* Boolean monitor
  robustness_monitor_one.h              RobustnessResult, STL*1 robustness monitor
  robustness_monitor_general.h          general STL* robustness monitor
  naive_monitors.h                      the four reference monitors (Boolean and exact robustness)
  constants.h                           time tolerance, robustness stopping rule
src/
  *.cpp                                 implementations of the headers above
  monitor_common.h/.cpp                 internals shared by our monitors
  robustness_common.h                   the bisection shared by both robustness monitors
examples/
  examples.cpp                          how to use every monitor
tests/
  test_utils.h                          reporting, printing and comparing results
  oracle.h                              independent evaluator of the STL* semantics
  random_generators.h                   random traces and formulas
  robustness_test_utils.h               shared robustness checks
  test_parser_and_api.cpp               parser, traces, low-level and one-call API
  test_boolean_monitor_one.cpp          STL*1 Boolean monitors, interval operators, constraint updater
  test_boolean_monitor_general.cpp      general Boolean monitors
  test_robustness_monitor_one.cpp       STL*1 robustness monitor
  test_robustness_monitor_general.cpp   general robustness monitor
```

---

## How correctness is tested

Our monitors were tested against the reference monitors:

- **Boolean monitors:** our results and the reference results **match exactly on every test**.
- **Robustness monitors:** our results are approximations. On every test, the reference monitor's exact robustness lies inside the returned bracket, and the bracket is as narrow as `epsilon` requires.

The reference monitors are in turn checked against an independent evaluator of the semantics in the tests, called the oracle. Every test program returns a non-zero exit code on failure, and `ctest` runs all five.

1. **Corner cases.** About 130 chosen formulas and traces are run through our monitors and the reference monitors. They cover:
   - every operator, including all four kinds of until;
   - ties between `>=` and `>`;
   - empty windows, and windows that fall between samples;
   - `dt = 0.1` rounding;
   - non-uniform traces and single-sample traces;
   - shadowed, sibling and nested dependent freeze operators;
   - for robustness: infinite robustness, non-monotone constraints such as `a* * a*`, values near 1e16 and near the largest double, Unix-time-like values, and `epsilon` from 0.1 down to 1e-6.
2. **Random formulas against the oracle.** Random formulas are generated together with a separate expression tree. The oracle computes that tree straight from the definitions in [Semantics](#semantics) and [Robustness](#robustness), with no library code.
   - **Boolean:** every monitor, ours and the reference ones, STL\*₁ and general, must give exactly the oracle's result. This runs on 25,000 random STL\*₁ formulas and 20,000 random general formulas.
   - **Robustness:** on 11,000 random formulas, the reference monitors must return exactly the oracle's robustness, and our monitors must bracket it within `epsilon`.
   - **Exact ties:** data values are exact binary fractions, so ties are tested deterministically.
   - **Operator coverage:** each test counts the random formulas that contain each operator, and fails if one is too rare. The operators counted are:
     - `!`, `&&`, `||`, timed and untimed `F` and `G`, and the four kinds of until: `U[a,b]`, `UNS[a,b]`, `U` and `UNS`;
     - freeze and nested freeze;
     - a freeze directly under `&&`, `||`, untimed `G` or untimed `F`, which are the cases where the acceleration heuristics apply.
3. **Traces of different sizes and types.** Our monitors and the reference monitors were compared on about 730 formulas. The traces had:
   - **sizes** from 1 to 5,000 samples;
   - **time grids** that were uniform, with several sampling periods, or non-uniform;
   - **values** that were noisy or slowly varying.

   Timestamps near 1e9 and values near 1e16 are also tested.
4. **Building blocks against brute force.** Each interval operator is compared with its pointwise definition on random sets and windows. `SignalConstraintUpdater` is checked after each of about 180,000 steps, including NaN and infinite values.
5. **Mutation checks.** Bugs were injected on purpose to confirm that the suite catches them, and every one was detected. The injected bugs were:
   - breaking conjunction, disjunction, strict and non-strict until, or the window lookups of the interval operators;
   - missing one truth-value flip in the signal constraint updater;
   - dropping the empty-window rule;
   - flipping a sign in the parser;
   - inverting each acceleration heuristic, or using a sibling before it is computed;
   - skipping the re-evaluation of nested freeze operators;
   - in the reference monitors: breaking the prefix counts, Lemire's filter, the empty-window robustness, and until (strict and non-strict, Boolean and robustness);
   - stopping the robustness bisection with a bracket ten times too wide.

The suite also passes in a debug build with libstdc++ bounds checking (`-D_GLIBCXX_DEBUG`) and the undefined-behaviour sanitizer.

---

## License

MIT. See [LICENSE](LICENSE).
