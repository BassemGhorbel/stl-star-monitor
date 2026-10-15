#include "stl_star/formula.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace stl_star {

// ---------------------------------------------------------------------------
// String helpers
// ---------------------------------------------------------------------------
std::string comp_op_to_string(CompOp op) {
    switch (op) {
        case CompOp::LT: return "<";
        case CompOp::LE: return "<=";
        case CompOp::GT: return ">";
        case CompOp::GE: return ">=";
    }
    return "?";
}

std::string formula_type_to_string(FormulaType t) {
    switch (t) {
        case FormulaType::SIGNAL_PRED:        return "SIGNAL_PRED";
        case FormulaType::SIGNAL_CONSTRAINT:  return "SIGNAL_CONSTRAINT";
        case FormulaType::NOT:                return "NOT";
        case FormulaType::AND:                return "AND";
        case FormulaType::OR:                 return "OR";
        case FormulaType::EVENTUALLY:         return "EVENTUALLY";
        case FormulaType::ALWAYS:             return "ALWAYS";
        case FormulaType::UNTIL:              return "UNTIL";
        case FormulaType::UNTIL_NS:           return "UNTIL_NS";
        case FormulaType::EVENTUALLY_UNTIMED: return "EVENTUALLY_UNTIMED";
        case FormulaType::ALWAYS_UNTIMED:     return "ALWAYS_UNTIMED";
        case FormulaType::FREEZE:             return "FREEZE";
    }
    return "UNKNOWN";
}

// ---------------------------------------------------------------------------
// Factories
// ---------------------------------------------------------------------------
FormulaNode FormulaNode::make_signal_pred(SignalFn h, CompOp op, const std::string& lbl) {
    FormulaNode n;
    n.type      = FormulaType::SIGNAL_PRED;
    n.signal_fn = std::move(h);
    n.comp_op   = op;
    n.label     = lbl;
    return n;
}

FormulaNode FormulaNode::make_signal_constraint(SignalFn h, FrozenFn g, CompOp op,
                                                std::vector<int> vars,
                                                const std::string& lbl) {
    FormulaNode n;
    n.type      = FormulaType::SIGNAL_CONSTRAINT;
    n.signal_fn = std::move(h);
    n.frozen_fn = std::move(g);
    n.comp_op   = op;
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    n.constraint_vars = std::move(vars);
    n.label    = lbl;
    return n;
}

static FormulaNode make_unary(FormulaType t, int child, const std::string& lbl) {
    FormulaNode n;
    n.type     = t;
    n.children = {child};
    n.label    = lbl;
    return n;
}

static FormulaNode make_timed(FormulaType t, double lo, double hi,
                              std::vector<int> children, const std::string& lbl) {
    FormulaNode n;
    n.type        = t;
    n.interval_lo = lo;
    n.interval_hi = hi;
    n.children    = std::move(children);
    n.label       = lbl;
    return n;
}

FormulaNode FormulaNode::make_not(int c, const std::string& l)               { return make_unary(FormulaType::NOT, c, l); }
FormulaNode FormulaNode::make_eventually_untimed(int c, const std::string& l) { return make_unary(FormulaType::EVENTUALLY_UNTIMED, c, l); }
FormulaNode FormulaNode::make_always_untimed(int c, const std::string& l)     { return make_unary(FormulaType::ALWAYS_UNTIMED, c, l); }

FormulaNode FormulaNode::make_and(int a, int b, const std::string& l) {
    FormulaNode n; n.type = FormulaType::AND; n.children = {a, b}; n.label = l; return n;
}
FormulaNode FormulaNode::make_or(int a, int b, const std::string& l) {
    FormulaNode n; n.type = FormulaType::OR; n.children = {a, b}; n.label = l; return n;
}
FormulaNode FormulaNode::make_eventually(double lo, double hi, int c, const std::string& l) {
    return make_timed(FormulaType::EVENTUALLY, lo, hi, {c}, l);
}
FormulaNode FormulaNode::make_always(double lo, double hi, int c, const std::string& l) {
    return make_timed(FormulaType::ALWAYS, lo, hi, {c}, l);
}
FormulaNode FormulaNode::make_until(double lo, double hi, int a, int b, const std::string& l) {
    return make_timed(FormulaType::UNTIL, lo, hi, {a, b}, l);
}
FormulaNode FormulaNode::make_until_ns(double lo, double hi, int a, int b, const std::string& l) {
    return make_timed(FormulaType::UNTIL_NS, lo, hi, {a, b}, l);
}

FormulaNode FormulaNode::make_freeze(int var, int comp, int child, const std::string& lbl) {
    FormulaNode n;
    n.type             = FormulaType::FREEZE;
    n.freeze_var       = var;
    n.freeze_component = comp;
    n.children         = {child};
    n.label            = lbl;
    return n;
}

// ---------------------------------------------------------------------------
// Formula
// ---------------------------------------------------------------------------
int Formula::add_node(FormulaNode node) {
    nodes.push_back(std::move(node));
    return static_cast<int>(nodes.size()) - 1;
}

static int expected_arity(FormulaType t) {
    switch (t) {
        case FormulaType::SIGNAL_PRED:
        case FormulaType::SIGNAL_CONSTRAINT: return 0;
        case FormulaType::AND:
        case FormulaType::OR:
        case FormulaType::UNTIL:
        case FormulaType::UNTIL_NS:          return 2;
        default:                             return 1;
    }
}

void Formula::validate() const {
    auto fail = [](int i, const std::string& msg) {
        throw std::invalid_argument("Formula node " + std::to_string(i) + ": " + msg);
    };
    const int N = static_cast<int>(nodes.size());
    if (N == 0) throw std::invalid_argument("Formula has no nodes");

    // Tree shape: children precede parents, every non-root node has one parent.
    std::vector<int> parent(N, -1);
    for (int i = 0; i < N; ++i) {
        const FormulaNode& nd = nodes[i];
        if (static_cast<int>(nd.children.size()) != expected_arity(nd.type))
            fail(i, formula_type_to_string(nd.type) + " has the wrong number of children");
        for (int c : nd.children) {
            if (c < 0 || c >= i) fail(i, "child index must be in [0, " + std::to_string(i) + ")");
            if (parent[c] != -1)  fail(c, "has more than one parent (formulas must be trees)");
            parent[c] = i;
        }
    }
    for (int i = 0; i + 1 < N; ++i)
        if (parent[i] == -1) fail(i, "is not reachable from the root");

    // Freeze variables: each bound by exactly one FREEZE node.
    std::vector<int> binder(num_freeze_vars, -1);
    for (int i = 0; i < N; ++i) {
        const FormulaNode& nd = nodes[i];
        if (nd.type != FormulaType::FREEZE) continue;
        if (nd.freeze_var < 0 || nd.freeze_var >= num_freeze_vars)
            fail(i, "freeze variable index out of range");
        if (binder[nd.freeze_var] != -1) fail(i, "freeze variable is bound twice");
        if (nd.freeze_component < 0) fail(i, "freeze component not set");
        binder[nd.freeze_var] = i;
    }

    for (int i = 0; i < N; ++i) {
        const FormulaNode& nd = nodes[i];
        switch (nd.type) {
            case FormulaType::SIGNAL_PRED:
                if (!nd.signal_fn) fail(i, "SIGNAL_PRED needs signal_fn");
                break;
            case FormulaType::SIGNAL_CONSTRAINT: {
                if (!nd.signal_fn || !nd.frozen_fn)
                    fail(i, "SIGNAL_CONSTRAINT needs signal_fn and frozen_fn");
                if (nd.constraint_vars.empty())
                    fail(i, "SIGNAL_CONSTRAINT must read at least one freeze variable");
                for (int v : nd.constraint_vars) {
                    if (v < 0 || v >= num_freeze_vars || binder[v] == -1)
                        fail(i, "reads an undefined freeze variable");
                    bool bound = false;   // is binder[v] an ancestor of i?
                    for (int p = parent[i]; p != -1; p = parent[p])
                        if (p == binder[v]) { bound = true; break; }
                    if (!bound) fail(i, "reads a freeze variable outside its scope");
                }
                break;
            }
            case FormulaType::EVENTUALLY:
            case FormulaType::ALWAYS:
            case FormulaType::UNTIL:
            case FormulaType::UNTIL_NS:
                if (!(nd.interval_lo >= 0.0) || !std::isfinite(nd.interval_lo) ||
                    !(nd.interval_hi >= nd.interval_lo))
                    fail(i, "invalid window [" + std::to_string(nd.interval_lo) + ", " +
                            std::to_string(nd.interval_hi) + "]: need 0 <= a <= b");
                break;
            default:
                break;
        }
    }
}

bool is_fragment_one(const Formula& f, std::string* why) {
    const int N = static_cast<int>(f.nodes.size());
    std::vector<int> parent(N, -1);
    for (int i = 0; i < N; ++i)
        for (int c : f.nodes[i].children) parent[c] = i;
    for (int i = 0; i < N; ++i) {
        const FormulaNode& nd = f.nodes[i];
        if (nd.type != FormulaType::SIGNAL_CONSTRAINT) continue;
        int p = parent[i];
        while (p != -1 && f.nodes[p].type != FormulaType::FREEZE) p = parent[p];
        if (p == -1 || nd.constraint_vars.size() != 1 ||
            nd.constraint_vars[0] != f.nodes[p].freeze_var) {
            if (why)
                *why = "a signal constraint reads a freeze variable other than (or in "
                       "addition to) the one bound by its nearest enclosing freeze "
                       "operator";
            return false;
        }
    }
    return true;
}

}  // namespace stl_star
