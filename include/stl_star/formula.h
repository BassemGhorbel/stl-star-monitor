#pragma once

#include <functional>
#include <string>
#include <vector>

namespace stl_star {

// ---------------------------------------------------------------------------
// Comparison operators of atomic propositions.
// ---------------------------------------------------------------------------
enum class CompOp { LT, LE, GT, GE };

std::string comp_op_to_string(CompOp op);

// lhs (op) rhs.  Any comparison involving NaN is false.
inline bool compare(CompOp op, double lhs, double rhs) {
    switch (op) {
        case CompOp::LT: return lhs <  rhs;
        case CompOp::LE: return lhs <= rhs;
        case CompOp::GT: return lhs >  rhs;
        case CompOp::GE: return lhs >= rhs;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Formula node types.
// ---------------------------------------------------------------------------
enum class FormulaType {
    // Atomic propositions
    SIGNAL_PRED,        // h(s) ~ 0                       (no freeze variable)
    SIGNAL_CONSTRAINT,  // h(s) ~ g(x*)                   (uses freeze variables)

    // Boolean
    NOT,                // !phi
    AND,                // phi && psi
    OR,                 // phi || psi

    // Timed temporal operators, closed windows [a, b] with 0 <= a <= b <= inf
    EVENTUALLY,         // F[a,b] phi
    ALWAYS,             // G[a,b] phi
    UNTIL,              // phi U[a,b] psi  (strict:     phi on [i, j), psi at j)
    UNTIL_NS,           // phi UNS[a,b] psi (non-strict: phi on [i, j], psi at j)

    // Untimed temporal operators (window [0, inf))
    EVENTUALLY_UNTIMED, // F phi
    ALWAYS_UNTIMED,     // G phi

    // Freeze binder
    FREEZE,             // x*. phi : bind x* to the current value of component x
};

std::string formula_type_to_string(FormulaType t);

// ---------------------------------------------------------------------------
// Atomic propositions are kept in SEPARATED form
//
//     SIGNAL_PRED:        signal_fn(s) ~ threshold
//     SIGNAL_CONSTRAINT:  signal_fn(s) - frozen_fn(x*) ~ threshold
//
// where s is the current sample (s[k] = component k) and x* is the vector of
// frozen values (x*[v] = value of freeze variable v).  `threshold` is 0 for
// Boolean monitoring; the robustness monitors shift it by +-r, so that the
// shifted atom holds exactly when the atom's robustness h - g is >= r.
//
// The separation is what makes the fast signal-constraint update possible:
// the signal side is sorted once, and each new frozen value only moves a
// boundary through that sorted order.  The parser produces this form from
// any expression in which no product/quotient mixes signals and freeze
// variables (e.g. "a - a* > 2*b*" is fine, "a * a* > 1" is rejected).
// Because the two sides are evaluated separately, an intermediate overflow
// can differ from the literal expression: "(a - a*) * 1e308 > 0" is evaluated
// as a*1e308 - (a*)*1e308, which overflows to inf - inf for a ~ 2.
// ---------------------------------------------------------------------------
using SignalFn = std::function<double(const std::vector<double>& s)>;
using FrozenFn = std::function<double(const std::vector<double>& frozen)>;

// ---------------------------------------------------------------------------
// A node of the formula tree.  Nodes live in Formula::nodes in
// reverse-topological order: every child index is smaller than its parent's,
// and the root is the last node.
//
// Fields used per type:
//   SIGNAL_PRED        signal_fn, comp_op, threshold
//   SIGNAL_CONSTRAINT  signal_fn, frozen_fn, constraint_vars, comp_op, threshold
//   NOT                children[0]
//   AND, OR            children[0], children[1]
//   EVENTUALLY, ALWAYS children[0], interval_lo, interval_hi
//   UNTIL, UNTIL_NS    children[0] (phi), children[1] (psi), interval_lo, interval_hi
//   *_UNTIMED          children[0]
//   FREEZE             children[0], freeze_var, freeze_component
// ---------------------------------------------------------------------------
struct FormulaNode {
    FormulaType type = FormulaType::SIGNAL_PRED;
    std::vector<int> children;

    // Atomic propositions
    SignalFn         signal_fn;
    FrozenFn         frozen_fn;
    std::vector<int> constraint_vars;   // freeze variables read by frozen_fn (sorted)
    CompOp           comp_op   = CompOp::GE;
    double           threshold = 0.0;

    // FREEZE
    int freeze_var       = -1;  // index of the bound variable (slot in x*)
    int freeze_component = -1;  // signal component whose current value is frozen

    // Timed temporal operators
    double interval_lo = 0.0;
    double interval_hi = 0.0;

    std::string label;         // optional, for debugging

    // ----- Factories --------------------------------------------------------
    // h(s) ~ 0
    static FormulaNode make_signal_pred(SignalFn h, CompOp op, const std::string& lbl = "");
    // h(s) ~ g(x*), where g reads the freeze variables listed in vars
    static FormulaNode make_signal_constraint(SignalFn h, FrozenFn g, CompOp op,
                                              std::vector<int> vars,
                                              const std::string& lbl = "");
    static FormulaNode make_not(int child, const std::string& lbl = "");
    static FormulaNode make_and(int left, int right, const std::string& lbl = "");
    static FormulaNode make_or (int left, int right, const std::string& lbl = "");
    static FormulaNode make_eventually(double lo, double hi, int child, const std::string& lbl  = "");
    static FormulaNode make_always     (double lo, double hi, int child, const std::string& lbl = "");
    static FormulaNode make_until      (double lo, double hi, int left, int right, const std::string& lbl = "");
    static FormulaNode make_until_ns   (double lo, double hi, int left, int right, const std::string& lbl = "");
    static FormulaNode make_eventually_untimed(int child, const std::string& lbl = "");
    static FormulaNode make_always_untimed    (int child, const std::string& lbl = "");
    static FormulaNode make_freeze(int freeze_var, int freeze_comp, int child, const std::string& lbl = "");
};

// ---------------------------------------------------------------------------
// A complete formula: a tree stored as a flat node array (root = last node).
// ---------------------------------------------------------------------------
struct Formula {
    std::vector<FormulaNode> nodes;
    int num_freeze_vars = 0;  // freeze variables are 0 .. num_freeze_vars-1
    int num_signals     = 0;  // the formula reads components 0 .. num_signals-1
                              // (set by the parser; 0 = unknown / not checked)
    std::string signal_names; // letter of each component, in alphabetical order
                              // (set by the parser: "ax" means a = 0, x = 1)

    int root_index() const { return static_cast<int>(nodes.size()) - 1; }
    int add_node(FormulaNode node);

    // Throws std::invalid_argument unless the formula is well formed:
    //   children precede parents and every node except the root has exactly
    //   one parent (a tree); arities, functions and windows are valid; every
    //   freeze variable is bound by exactly one FREEZE node; and every signal
    //   constraint only reads variables bound by an enclosing FREEZE node.
    void validate() const;
};

// Is the formula in the STL*_1 fragment?  This holds when every signal
// constraint reads exactly one freeze variable: the one bound by its nearest
// enclosing FREEZE node.  If not and `why` is non-null, a reason is stored.
bool is_fragment_one(const Formula& f, std::string* why = nullptr);

}  // namespace stl_star
