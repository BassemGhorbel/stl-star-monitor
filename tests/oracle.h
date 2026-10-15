#pragma once
// ---------------------------------------------------------------------------
// Independent ORACLE for STL* semantics.
//
// An oracle formula is a small tree that is generated together with its text
// form (to_text).  The oracle evaluates the tree directly from the STL*
// definitions: pointwise, recursive, without intervals and without any code
// of the monitors except the two window predicates of trace.h (which define
// the time tolerance).  Comparing  parse_formula(text) + monitor  with the
// oracle therefore checks the parser AND the monitors.
//
//   sat(node, trace, i, env)   Boolean semantics at sample i
//   rob(node, trace, i, env)   robustness at sample i (negation-free formulas)
// ---------------------------------------------------------------------------
#include "test_utils.h"

#include <memory>

namespace test_util {


// ===========================================================================
// Oracle formulas
// ===========================================================================
struct OracleExpr;
using OracleExprPtr = std::shared_ptr<OracleExpr>;
struct OracleExpr {
    enum Kind { Num, Sig, Frz, Add, Sub, Mul, Div, Neg } kind;
    double num = 0.0;
    int    idx = -1;     // Sig: component, Frz: binder id
    OracleExprPtr     l, r;
};

inline OracleExprPtr expr_num(double v)  { auto e = std::make_shared<OracleExpr>(); e->kind = OracleExpr::Num; e->num = v; return e; }
inline OracleExprPtr expr_signal(int c)     { auto e = std::make_shared<OracleExpr>(); e->kind = OracleExpr::Sig; e->idx = c; return e; }
inline OracleExprPtr expr_frozen(int b)     { auto e = std::make_shared<OracleExpr>(); e->kind = OracleExpr::Frz; e->idx = b; return e; }
inline OracleExprPtr expr_binary(OracleExpr::Kind k, OracleExprPtr l, OracleExprPtr r) { auto e = std::make_shared<OracleExpr>(); e->kind = k; e->l = l; e->r = r; return e; }
inline OracleExprPtr expr_neg(OracleExprPtr x)      { auto e = std::make_shared<OracleExpr>(); e->kind = OracleExpr::Neg; e->l = x; return e; }

struct OracleNode;
using OracleNodePtr = std::shared_ptr<OracleNode>;
struct OracleNode {
    enum Kind { Atom, Not, And, Or, Ev, Al, EvU, AlU, Until, UntilNS, Freeze } kind;
    OracleExprPtr     lhs, rhs;            // Atom
    CompOp op = CompOp::GE;     // Atom
    double a = 0.0, b = 0.0;    // Ev, Al, Until, UntilNS (b may be +inf)
    int    binder = -1;         // Freeze: binder id
    int    comp   = -1;         // Freeze: frozen component
    OracleNodePtr     c0, c1;
};

struct OracleFormula {
    OracleNodePtr root;
    std::vector<std::string> binder_names;   // text name of each binder id ("a*", "b*_3")
    std::string letters;                     // letter of each component, sorted ("dqx")
    int num_binders() const { return static_cast<int>(binder_names.size()); }
    int num_signals() const { return static_cast<int>(letters.size()); }
};

// ---- Text form ---------------------------------------------------------------
inline std::string to_text(const OracleExprPtr& e, const OracleFormula& F) {
    switch (e->kind) {
        case OracleExpr::Num: return e->num < 0 ? "(0 - " + num(-e->num) + ")" : num(e->num);
        case OracleExpr::Sig: return std::string(1, F.letters[e->idx]);
        case OracleExpr::Frz: return F.binder_names[e->idx];
        case OracleExpr::Neg: return "-" + to_text(e->l, F);
        default: break;
    }
    const char* op = e->kind == OracleExpr::Add ? " + " : e->kind == OracleExpr::Sub ? " - "
                   : e->kind == OracleExpr::Mul ? " * " : " / ";
    return "(" + to_text(e->l, F) + op + to_text(e->r, F) + ")";
}

inline std::string window(double a, double b) { return "[" + num(a) + "," + num(b) + "]"; }
// Until windows: [0, inf) is written in the untimed form "(phi U psi)".
inline std::string until_window(double a, double b) {
    return a == 0 && std::isinf(b) ? "" : window(a, b);
}

inline std::string to_text(const OracleNodePtr& n, const OracleFormula& F) {
    switch (n->kind) {
        case OracleNode::Atom:    return to_text(n->lhs, F) + " " + comp_op_to_string(n->op) + " " + to_text(n->rhs, F);
        case OracleNode::Not:     return "! (" + to_text(n->c0, F) + ")";
        case OracleNode::And:     return "(" + to_text(n->c0, F) + " && " + to_text(n->c1, F) + ")";
        case OracleNode::Or:      return "(" + to_text(n->c0, F) + " || " + to_text(n->c1, F) + ")";
        case OracleNode::Ev:      return "F" + window(n->a, n->b) + " (" + to_text(n->c0, F) + ")";
        case OracleNode::Al:      return "G" + window(n->a, n->b) + " (" + to_text(n->c0, F) + ")";
        case OracleNode::EvU:     return "F (" + to_text(n->c0, F) + ")";
        case OracleNode::AlU:     return "G (" + to_text(n->c0, F) + ")";
        case OracleNode::Until:   return "((" + to_text(n->c0, F) + ") U" + until_window(n->a, n->b) + " (" + to_text(n->c1, F) + "))";
        case OracleNode::UntilNS: return "((" + to_text(n->c0, F) + ") UNS" + until_window(n->a, n->b) + " (" + to_text(n->c1, F) + "))";
        case OracleNode::Freeze:  return F.binder_names[n->binder] + ". (" + to_text(n->c0, F) + ")";
    }
    return "?";
}
inline std::string to_text(const OracleFormula& F) { return to_text(F.root, F); }

// ---- Evaluation ----------------------------------------------------------------
inline double eval_expr(const OracleExprPtr& e, const Trace& tr, int i, const std::vector<double>& env) {
    switch (e->kind) {
        case OracleExpr::Num: return e->num;
        case OracleExpr::Sig: return tr.value(e->idx, i);
        case OracleExpr::Frz: return env[e->idx];
        case OracleExpr::Neg: return -eval_expr(e->l, tr, i, env);
        case OracleExpr::Add: return eval_expr(e->l, tr, i, env) + eval_expr(e->r, tr, i, env);
        case OracleExpr::Sub: return eval_expr(e->l, tr, i, env) - eval_expr(e->r, tr, i, env);
        case OracleExpr::Mul: return eval_expr(e->l, tr, i, env) * eval_expr(e->r, tr, i, env);
        case OracleExpr::Div: return eval_expr(e->l, tr, i, env) / eval_expr(e->r, tr, i, env);
    }
    return 0.0;
}

// Is tau_k in the window [a, b] of tau_i?  (same tolerance as the library)
inline bool in_window(const Trace& tr, int i, int k, double a, double b) {
    return tr.reaches_window_start(i, k, a) && tr.within_window_end(i, k, b);
}

// Boolean semantics (rho, i, env) |= phi, straight from the definition.
inline bool sat(const OracleNodePtr& n, const Trace& tr, int i, std::vector<double>& env) {
    const int len = tr.length;
    switch (n->kind) {
        case OracleNode::Atom:
            return compare(n->op, eval_expr(n->lhs, tr, i, env), eval_expr(n->rhs, tr, i, env));
        case OracleNode::Not: return !sat(n->c0, tr, i, env);
        case OracleNode::And: return sat(n->c0, tr, i, env) && sat(n->c1, tr, i, env);
        case OracleNode::Or:  return sat(n->c0, tr, i, env) || sat(n->c1, tr, i, env);
        case OracleNode::Ev:  case OracleNode::EvU: {
            const double a = n->kind == OracleNode::Ev ? n->a : 0.0;
            const double b = n->kind == OracleNode::Ev ? n->b : std::numeric_limits<double>::infinity();
            for (int k = i; k < len; ++k)
                if (in_window(tr, i, k, a, b) && sat(n->c0, tr, k, env)) return true;
            return false;
        }
        case OracleNode::Al:  case OracleNode::AlU: {
            const double a = n->kind == OracleNode::Al ? n->a : 0.0;
            const double b = n->kind == OracleNode::Al ? n->b : std::numeric_limits<double>::infinity();
            for (int k = i; k < len; ++k)
                if (in_window(tr, i, k, a, b) && !sat(n->c0, tr, k, env)) return false;
            return true;
        }
        case OracleNode::Until: case OracleNode::UntilNS: {
            // exists k >= i in the window with psi(k), and phi on [i, k) (strict)
            // or on [i, k] (non-strict)
            const bool strict = n->kind == OracleNode::Until;
            for (int k = i; k < len; ++k) {
                if (!in_window(tr, i, k, n->a, n->b)) continue;
                if (!sat(n->c1, tr, k, env)) continue;
                bool phi_ok = true;
                for (int l = i; l < (strict ? k : k + 1) && phi_ok; ++l) phi_ok = sat(n->c0, tr, l, env);
                if (phi_ok) return true;
            }
            return false;
        }
        case OracleNode::Freeze: {
            const double saved = env[n->binder];
            env[n->binder] = tr.value(n->comp, i);
            const bool r = sat(n->c0, tr, i, env);
            env[n->binder] = saved;
            return r;
        }
    }
    return false;
}

// Quantitative semantics (negation-free formulas).
inline double rob(const OracleNodePtr& n, const Trace& tr, int i, std::vector<double>& env) {
    const double inf = std::numeric_limits<double>::infinity();
    const int len = tr.length;
    switch (n->kind) {
        case OracleNode::Atom: {
            const double d = eval_expr(n->lhs, tr, i, env) - eval_expr(n->rhs, tr, i, env);
            return (n->op == CompOp::GT || n->op == CompOp::GE) ? d : -d;
        }
        case OracleNode::And: return std::min(rob(n->c0, tr, i, env), rob(n->c1, tr, i, env));
        case OracleNode::Or:  return std::max(rob(n->c0, tr, i, env), rob(n->c1, tr, i, env));
        case OracleNode::Ev:  case OracleNode::EvU: {
            const double a = n->kind == OracleNode::Ev ? n->a : 0.0, b = n->kind == OracleNode::Ev ? n->b : inf;
            double r = -inf;
            for (int k = i; k < len; ++k) if (in_window(tr, i, k, a, b)) r = std::max(r, rob(n->c0, tr, k, env));
            return r;
        }
        case OracleNode::Al:  case OracleNode::AlU: {
            const double a = n->kind == OracleNode::Al ? n->a : 0.0, b = n->kind == OracleNode::Al ? n->b : inf;
            double r = inf;
            for (int k = i; k < len; ++k) if (in_window(tr, i, k, a, b)) r = std::min(r, rob(n->c0, tr, k, env));
            return r;
        }
        case OracleNode::Until: case OracleNode::UntilNS: {
            const bool strict = n->kind == OracleNode::Until;
            double best = -inf;
            for (int k = i; k < len; ++k) {
                if (!in_window(tr, i, k, n->a, n->b)) continue;
                double v = rob(n->c1, tr, k, env);
                for (int l = i; l < (strict ? k : k + 1); ++l) v = std::min(v, rob(n->c0, tr, l, env));
                best = std::max(best, v);
            }
            return best;
        }
        case OracleNode::Freeze: {
            const double saved = env[n->binder];
            env[n->binder] = tr.value(n->comp, i);
            const double r = rob(n->c0, tr, i, env);
            env[n->binder] = saved;
            return r;
        }
        case OracleNode::Not: break;
    }
    return std::numeric_limits<double>::quiet_NaN();   // negation: not defined here
}

inline Intervals oracle_intervals(const OracleFormula& F, const Trace& tr) {
    std::vector<double> env(F.num_binders(), 0.0);
    return intervals_where(tr, [&](int i) { return sat(F.root, tr, i, env); });
}
inline double oracle_robustness(const OracleFormula& F, const Trace& tr) {
    std::vector<double> env(F.num_binders(), 0.0);
    return rob(F.root, tr, 0, env);
}

}  // namespace test_util
