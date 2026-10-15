#pragma once
// ---------------------------------------------------------------------------
// Random traces and random oracle formulas.
//
// Signal values, constants and coefficients are dyadic rationals (multiples
// of 1/64), so every arithmetic operation in an atom is exact: the oracle's
// "lhs - rhs" and the library's separated "h(s) - g(x*)" agree exactly, and
// ties (where >= and > differ) are exercised deterministically.
// ---------------------------------------------------------------------------
#include "oracle.h"

#include <map>
#include <random>
#include <set>

namespace test_util {


// ===========================================================================
// Random generators
// ===========================================================================
struct Rng {
    std::mt19937 g;
    explicit Rng(unsigned seed) : g(seed) {}
    int    uni(int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(g); }
    double real(double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(g); }
    bool   coin(double p = 0.5) { return real(0.0, 1.0) < p; }
    template <class T> const T& pick(const std::vector<T>& v) { return v[uni(0, (int)v.size() - 1)]; }
};

// Dyadic value: multiple of 1/64, so sums/products stay exact.
inline double dyadic(double x) { return std::round(x * 64.0) / 64.0; }

// Random trace.  Values: small integers (many ties), halves, or fine dyadics.
// Timestamps: uniform with assorted dt (including inexact 0.1 / 0.3), or
// irregular gaps (tiny, regular, and large), sometimes starting at t0 != 0.
inline Trace random_trace(Rng& R, int n, int comps, int uniform_mode = -1) {
    const int value_mode = R.uni(0, 2);
    std::vector<std::vector<double>> sig(comps, std::vector<double>(n));
    for (auto& c : sig)
        for (auto& v : c) {
            if (value_mode == 0)      v = R.uni(-2, 3);
            else if (value_mode == 1) v = 0.5 * R.uni(-6, 6);
            else                     v = dyadic(R.real(-4.0, 4.0));
        }
    const bool uniform = uniform_mode < 0 ? R.coin(0.5) : uniform_mode == 1;
    if (uniform) {
        static const std::vector<double> dts = {1.0, 1.0, 0.5, 0.25, 2.0, 0.1, 0.3, 0.75};
        return Trace::make_uniform(R.pick(dts), sig);
    }
    std::vector<double> ts(n);
    const int gap_mode = R.uni(0, 2);
    ts[0] = R.coin(0.3) ? 0.25 * R.uni(1, 20) : 0.0;
    for (int i = 1; i < n; ++i) {
        double gap;
        if (gap_mode == 0)      gap = 0.25 * R.uni(1, 10);
        else if (gap_mode == 1) gap = dyadic(R.real(0.05, 2.5)) + 1.0 / 64;
        else                   gap = R.pick(std::vector<double>{0.01, 0.1, 0.5, 1.0, 1.0, 2.0, 5.0});
        ts[i] = ts[i - 1] + gap;
    }
    return Trace::make_non_uniform(ts, sig);
}

// Trace sizes of the random size tests: from a single sample up to 1000.
// (Traces of 5000 samples are tested with a few chosen formulas, since the
// reference monitors re-evaluate every freeze scope at every sample.)
inline int random_trace_size(Rng& R) {
    return R.pick(std::vector<int>{1, 2, 5, 10, 50, 100, 200, 500, 1000});
}

// A smooth or noisy signal of n samples on a uniform (dt = 0.5) or
// irregular time grid, for the chosen 5000-sample cases.
inline Trace trace_5000_style(int n, bool smooth, bool uniform, unsigned seed) {
    Rng R(seed);
    std::vector<double> sig(n), ts(n);
    for (int i = 0; i < n; ++i) {
        sig[i] = smooth ? dyadic(3.0 * std::sin(i * 6.283 * 7 / n) + 0.5 * std::sin(i * 0.9))
                        : dyadic(R.real(-4.0, 4.0));
        ts[i] = i == 0 ? 0.0 : ts[i - 1] + R.pick(std::vector<double>{0.25, 0.5, 1.0});
    }
    return uniform ? Trace::make_uniform(0.5, {sig}) : Trace::make_non_uniform(ts, {sig});
}

// Traces of different types: half are the noisy traces of random_trace, half
// are slowly varying signals (two sines, rounded to dyadic values) on a
// uniform or irregular time grid.
inline Trace varied_trace(Rng& R, int n, int comps) {
    if (R.coin(0.5)) return random_trace(R, n, comps);
    std::vector<std::vector<double>> sig(comps, std::vector<double>(n));
    for (auto& c : sig) {
        const double p1 = R.real(10.0, 10.0 + n), p2 = R.real(3.0, 30.0), phase = R.real(0.0, 6.283);
        for (int i = 0; i < n; ++i)
            c[i] = dyadic(3.0 * std::sin(6.283 * i / p1 + phase) + 0.5 * std::sin(6.283 * i / p2));
    }
    if (R.coin(0.5))
        return Trace::make_uniform(R.pick(std::vector<double>{1.0, 0.5, 0.1, 2.0}), sig);
    std::vector<double> ts(n);
    for (int i = 1; i < n; ++i) ts[i] = ts[i - 1] + R.pick(std::vector<double>{0.25, 0.5, 1.0, 1.5});
    return Trace::make_non_uniform(ts, sig);
}

struct GenOptions {
    int    comps            = 1;      // signals the generator may use (at most 26)
    bool   fragment_one     = true;   // STL*_1: constraints read only the nearest binder
    bool   negation_free    = false;
    int    max_depth        = 4;
    int    max_freeze_depth = 2;
    bool   allow_until      = true;
    double freeze_bias      = 0.0;    // extra probability of opening a freeze scope
};

class FormulaGen {
public:
    FormulaGen(Rng& R, GenOptions o) : R_(R), o_(o) {}

    // The components the formula uses are renumbered 0 .. m-1 and named by m
    // letters in alphabetical order (half the time a, b, ..., otherwise random
    // letters such as "dqx"), so the formula needs a trace of exactly
    // F.num_signals() signals.
    OracleFormula make() {
        F_ = OracleFormula{};
        scope_.clear();
        binder_comp_.clear();
        binder_suffix_.clear();
        F_.root = node(o_.max_depth);

        std::vector<int> used(o_.comps, 0);
        visit(F_.root, [&](int& c) { used[c] = 1; });
        std::vector<int> rank(o_.comps, -1);
        int m = 0;
        for (int c = 0; c < o_.comps; ++c) if (used[c]) rank[c] = m++;
        visit(F_.root, [&](int& c) { c = rank[c]; });
        for (int& c : binder_comp_) c = rank[c];

        std::string all = "abcdefghijklmnopqrstuvwxyz";
        if (R_.coin(0.5)) std::shuffle(all.begin(), all.end(), R_.g);
        F_.letters = all.substr(0, m);
        std::sort(F_.letters.begin(), F_.letters.end());
        for (size_t b = 0; b < binder_comp_.size(); ++b)
            F_.binder_names.push_back(std::string(1, F_.letters[binder_comp_[b]]) + "*" + binder_suffix_[b]);
        return F_;
    }

private:
    Rng& R_;
    GenOptions o_;
    OracleFormula F_;
    std::vector<int> scope_;                  // binder ids, innermost last
    std::vector<int> binder_comp_;            // frozen component of each binder id
    std::vector<std::string> binder_suffix_;  // "" or "_k" for each binder id

    // Calls fn on every component index of the formula (signals and freezes).
    template <class Fn> static void visit(const OracleExprPtr& e, Fn& fn) {
        if (!e) return;
        if (e->kind == OracleExpr::Sig) fn(e->idx);
        visit(e->l, fn);
        visit(e->r, fn);
    }
    template <class Fn> static void visit(const OracleNodePtr& n, Fn&& fn) {
        if (!n) return;
        if (n->kind == OracleNode::Freeze) fn(n->comp);
        visit(n->lhs, fn);
        visit(n->rhs, fn);
        visit(n->c0, fn);
        visit(n->c1, fn);
    }

    double constant() { return R_.pick(std::vector<double>{-2, -1, -0.5, 0, 0, 0.5, 1, 1.5, 2, 3}); }
    double coef()     { return R_.pick(std::vector<double>{1, 1, 1, 2, 0.5, -1, 3}); }

    OracleExprPtr scaled(OracleExprPtr x) {
        double c = coef();
        if (c == 1) return x;
        if (c == -1) return expr_neg(x);
        return R_.coin() ? expr_binary(OracleExpr::Mul, expr_num(c), x) : expr_binary(OracleExpr::Mul, x, expr_num(c));
    }
    OracleExprPtr sig_term() { return scaled(expr_signal(R_.uni(0, o_.comps - 1))); }
    OracleExprPtr sig_expr() {
        OracleExprPtr e = sig_term();
        if (R_.coin(0.3)) e = expr_binary(R_.coin() ? OracleExpr::Add : OracleExpr::Sub, e, sig_term());
        if (R_.coin(0.3)) e = expr_binary(OracleExpr::Add, e, expr_num(constant()));
        if (R_.coin(0.1)) e = expr_binary(OracleExpr::Div, e, expr_num(R_.pick(std::vector<double>{2, 4, 0.5})));
        return e;
    }
    OracleExprPtr frz_expr(const std::vector<int>& vars) {
        OracleExprPtr e = scaled(expr_frozen(vars[0]));
        for (size_t k = 1; k < vars.size(); ++k)
            e = expr_binary(R_.coin() ? OracleExpr::Add : OracleExpr::Sub, e, scaled(expr_frozen(vars[k])));
        if (R_.coin(0.3)) e = expr_binary(OracleExpr::Add, e, expr_num(constant()));
        return e;
    }

    OracleNodePtr atom() {
        auto n = std::make_shared<OracleNode>();
        n->kind = OracleNode::Atom;
        n->op = R_.pick(std::vector<CompOp>{CompOp::GT, CompOp::GE, CompOp::LT, CompOp::LE});
        if (scope_.empty() || R_.coin(0.3)) {             // signal predicate
            n->lhs = sig_expr();
            n->rhs = R_.coin(0.7) ? expr_num(constant()) : sig_expr();
            return n;
        }
        std::vector<int> vars;                             // signal constraint
        if (o_.fragment_one) {
            vars.push_back(scope_.back());
        } else {
            // Mostly several frozen values and/or outer ones (outside STL*_1).
            for (int b : scope_) if (R_.coin(0.65)) vars.push_back(b);
            if (vars.empty()) vars.push_back(R_.pick(scope_));
        }
        OracleExprPtr s = sig_expr(), z = frz_expr(vars);
        switch (R_.uni(0, 4)) {
            case 0:  n->lhs = s; n->rhs = z; break;                            // s ~ z
            case 1:  n->lhs = z; n->rhs = s; break;                            // z ~ s
            case 2:  n->lhs = expr_binary(OracleExpr::Sub, s, z); n->rhs = expr_num(constant()); break;
            case 3:  n->lhs = expr_binary(OracleExpr::Mul, expr_num(R_.pick(std::vector<double>{2, -1, 0.5})),
                                   expr_binary(OracleExpr::Sub, s, z));
                     n->rhs = expr_num(constant()); break;                         // c*(s - z) ~ c'
            default: n->lhs = expr_binary(OracleExpr::Div, expr_binary(OracleExpr::Add, s, z), expr_num(2));
                     n->rhs = sig_expr(); break;                                // (s + z)/2 ~ s'
        }
        return n;
    }

    void timed_window(OracleNode& n) {
        n.a = R_.pick(std::vector<double>{0, 0, 0, 0.25, 0.5, 1, 1, 1.5, 2, 3});
        n.b = n.a + R_.pick(std::vector<double>{0, 0.1, 0.25, 0.5, 1, 1, 2, 2, 3.5,
                                                std::numeric_limits<double>::infinity()});
    }

    // Untimed until "(phi U psi)" (window [0, inf)) in a quarter of the cases.
    void until_window(OracleNode& n) {
        if (R_.coin(0.25)) { n.a = 0; n.b = std::numeric_limits<double>::infinity(); }
        else timed_window(n);
    }

    OracleNodePtr node(int depth) {
        if (depth <= 0) return atom();
        std::vector<int> kinds = {0, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 10, 10};
        if (!o_.negation_free) kinds.push_back(9);
        if (o_.allow_until) { kinds.push_back(8); kinds.push_back(11); }
        if ((int)scope_.size() < o_.max_freeze_depth) {
            for (int k = 0; k < (o_.fragment_one ? 3 : 7); ++k) kinds.push_back(12);
        }
        auto n = std::make_shared<OracleNode>();
        const bool force_freeze = (int)scope_.size() < o_.max_freeze_depth && R_.coin(o_.freeze_bias);
        switch (force_freeze ? 12 : R_.pick(kinds)) {
            case 0: case 1: return atom();
            case 2: n->kind = OracleNode::And; break;
            case 3: n->kind = R_.coin() ? OracleNode::And : OracleNode::Or; break;
            case 4: n->kind = OracleNode::Ev; timed_window(*n); break;
            case 5: n->kind = OracleNode::Al; timed_window(*n); break;
            case 6: n->kind = OracleNode::EvU; break;
            case 7: n->kind = OracleNode::AlU; break;
            case 8: n->kind = OracleNode::Until; until_window(*n); break;
            case 9: n->kind = OracleNode::Not; break;
            case 10: n->kind = R_.coin() ? OracleNode::Ev : OracleNode::Al; timed_window(*n); break;
            case 11: n->kind = OracleNode::UntilNS; until_window(*n); break;
            default: {                                   // freeze binder
                n->kind = OracleNode::Freeze;
                n->comp = R_.uni(0, o_.comps - 1);
                n->binder = static_cast<int>(binder_comp_.size());
                // STL*_1 constraints only read the innermost binder, so plain
                // names (which shadow) are fine; general formulas may read
                // outer binders, so every binder gets a unique subscript.
                // The name ("x*_3") is completed by make().
                binder_comp_.push_back(n->comp);
                binder_suffix_.push_back(!o_.fragment_one || R_.coin(0.3)
                                             ? "_" + std::to_string(n->binder + 1) : "");
                scope_.push_back(n->binder);
                n->c0 = node(depth - 1);
                scope_.pop_back();
                return n;
            }
        }
        n->c0 = node(depth - 1);
        if (n->kind == OracleNode::And || n->kind == OracleNode::Or ||
            n->kind == OracleNode::Until || n->kind == OracleNode::UntilNS)
            n->c1 = node(depth - 1);
        return n;
    }
};

// ===========================================================================
// Operator coverage of the random formulas
// ===========================================================================
// Counts how many random formulas contain each operator, each kind of until
// (strict / non-strict, timed / untimed), and each parent of a freeze operator
// that triggers an acceleration heuristic of the fast monitors.  A test calls
// require() at the end, so that a change in the generator cannot silently
// stop testing an operator.
class Coverage {
public:
    void add(const OracleFormula& F) {
        std::set<std::string> seen;
        visit(F.root, nullptr, 0, seen);
        for (const auto& s : seen) ++count_[s];
    }

    void require(const std::vector<std::string>& names, long at_least = 20) {
        std::printf("  operator coverage (formulas containing each):");
        int col = 0;
        for (const auto& name : names) {
            std::printf("%s %s: %ld", col++ % 5 == 0 ? "\n   " : ",", name.c_str(), count_[name]);
            check(count_[name] >= at_least, "random formulas cover \"" + name + "\"",
                  "    only " + std::to_string(count_[name]) + " formulas");
        }
        std::printf("\n");
    }

    static std::vector<std::string> all_operators() {
        return {"!", "&&", "||", "F[a,b]", "G[a,b]", "F", "G", "U[a,b]", "UNS[a,b]", "U", "UNS",
                "freeze", "nested freeze", "freeze in &&", "freeze in ||", "freeze &&/|| freeze",
                "freeze under G", "freeze under F"};
    }
    static std::vector<std::string> negation_free_operators() {
        std::vector<std::string> v = all_operators();
        v.erase(v.begin());   // "!"
        return v;
    }

private:
    std::map<std::string, long> count_;

    void visit(const OracleNodePtr& n, const OracleNode* parent, int freezes_above,
               std::set<std::string>& seen) {
        if (!n) return;
        const bool untimed = n->a == 0 && std::isinf(n->b);
        switch (n->kind) {
            case OracleNode::Not:     seen.insert("!"); break;
            case OracleNode::And:     seen.insert("&&"); break;
            case OracleNode::Or:      seen.insert("||"); break;
            case OracleNode::Ev:      seen.insert("F[a,b]"); break;
            case OracleNode::Al:      seen.insert("G[a,b]"); break;
            case OracleNode::EvU:     seen.insert("F"); break;
            case OracleNode::AlU:     seen.insert("G"); break;
            case OracleNode::Until:   seen.insert(untimed ? "U" : "U[a,b]"); break;
            case OracleNode::UntilNS: seen.insert(untimed ? "UNS" : "UNS[a,b]"); break;
            case OracleNode::Freeze:
                seen.insert("freeze");
                if (freezes_above > 0) seen.insert("nested freeze");
                if (parent && parent->kind == OracleNode::And) seen.insert("freeze in &&");
                if (parent && parent->kind == OracleNode::Or)  seen.insert("freeze in ||");
                if (parent && parent->kind == OracleNode::AlU) seen.insert("freeze under G");
                if (parent && parent->kind == OracleNode::EvU) seen.insert("freeze under F");
                if (parent && (parent->kind == OracleNode::And || parent->kind == OracleNode::Or) &&
                    parent->c0->kind == OracleNode::Freeze && parent->c1->kind == OracleNode::Freeze)
                    seen.insert("freeze &&/|| freeze");
                break;
            default: break;
        }
        const int above = freezes_above + (n->kind == OracleNode::Freeze ? 1 : 0);
        visit(n->c0, n.get(), above, seen);
        visit(n->c1, n.get(), above, seen);
    }
};

}  // namespace test_util
