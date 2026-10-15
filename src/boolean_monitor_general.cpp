#include "stl_star/boolean_monitor_general.h"
#include "monitor_common.h"
#include "stl_star/signal_constraint_updater.h"

#include <algorithm>
#include <vector>

namespace stl_star {

using namespace monitor_detail;

namespace {

// General STL*: a nested FREEZE is assumed to depend on the frozen values of
// the freeze operators around it, so it is re-evaluated at every step of its
// nearest enclosing FREEZE (exactly as in the naive monitor).  A FREEZE that
// is not nested in another one is evaluated once.
struct Scope {
    std::vector<int> constraints;   // signal constraints whose nearest FREEZE is this one
    std::vector<int> invariant;     // scope nodes that do not change with this FREEZE's
                                    // value (computed once per evaluation)
    std::vector<int> steps;         // all scope nodes except atoms (nested FREEZE nodes
                                    // included), in evaluation order; the loop skips
                                    // the invariant ones
    FreezeCtx        ctx;
};

class GeneralMonitor {
public:
    GeneralMonitor(const Formula& f, const Trace& trace) : f_(f), tr_(trace) {}

    IndexIntervals run() {
        check_input(f_, tr_, /*require_fragment_one=*/false);
        const int N = static_cast<int>(f_.nodes.size());
        n_ = tr_.length;
        parent_  = parents_of(f_);
        owner_   = freeze_owners(f_);
        windows_ = compute_window_sets(f_, tr_);
        iv_.assign(N, {});
        frozen_.assign(f_.num_freeze_vars, 0.0);
        scopes_.assign(N, {});
        updaters_.assign(N, {});
        points_.assign(N, {});

        // Scope membership.  A scope node "varies" when its subtree (inside
        // the scope) contains a signal constraint of the scope or a nested
        // FREEZE; all other scope nodes (built from signal predicates only)
        // are constant during one evaluation of the FREEZE.
        varies_.assign(N, 0);
        for (int k = 0; k < N; ++k) {
            const FormulaNode& nd = f_.nodes[k];
            if (nd.type == FormulaType::SIGNAL_PRED)
                iv_[k] = signal_pred_intervals(nd, tr_);       // frozen-value independent
            if (nd.type == FormulaType::SIGNAL_CONSTRAINT)
                updaters_[k].init(nd, tr_);                   // keys sorted once
            if (nd.type == FormulaType::FREEZE) points_[k].assign(n_, 0);
            const int own = owner_[k];
            if (own == -1) continue;
            if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
                scopes_[own].constraints.push_back(k);
                varies_[k] = 1;
            } else if (nd.type == FormulaType::FREEZE) {
                varies_[k] = 1;
            } else if (nd.type != FormulaType::SIGNAL_PRED) {
                for (int c : nd.children) varies_[k] = varies_[k] || varies_[c];
                if (!varies_[k]) scopes_[own].invariant.push_back(k);
            }
        }

        // Evaluation order of every scope, and of the nodes outside all scopes.
        position_.assign(N, -1);
        std::vector<int> top;
        order(f_.root_index(), top);
        for (int j = 0; j < N; ++j)
            if (f_.nodes[j].type == FormulaType::FREEZE) order(f_.nodes[j].children[0], scopes_[j].steps);
        for (int j = 0; j < N; ++j)
            if (f_.nodes[j].type == FormulaType::FREEZE) scopes_[j].ctx = freeze_ctx(j);

        for (int k : top) {
            if (f_.nodes[k].type == FormulaType::FREEZE) eval_freeze(k, 0);
            else iv_[k] = compute_node(k, f_, tr_, windows_, iv_);
        }
        return iv_[f_.root_index()];
    }

private:
    // Post-order of the nodes of one scope (or of the nodes outside all
    // scopes), without signal predicates and constraints, and without the
    // nodes inside nested FREEZE nodes.  For (phi1 && x*.phi) and
    // (phi1 || x*.phi), phi1 comes first, so that the FREEZE can use it (see
    // freeze_ctx).
    void order(int k, std::vector<int>& out) {
        const FormulaNode& nd = f_.nodes[k];
        if (nd.type == FormulaType::SIGNAL_PRED || nd.type == FormulaType::SIGNAL_CONSTRAINT) return;
        if (nd.type != FormulaType::FREEZE) {
            std::vector<int> children = nd.children;
            if (children.size() == 2 && f_.nodes[children[0]].type == FormulaType::FREEZE &&
                f_.nodes[children[1]].type != FormulaType::FREEZE)
                std::swap(children[0], children[1]);
            for (int c : children) order(c, out);
        }
        position_[k] = static_cast<int>(out.size());
        out.push_back(k);
    }

    // Heuristics of FREEZE node j (see monitor_common.h).  The sibling phi1 of
    // (phi1 && x*.phi) / (phi1 || x*.phi) must be computed before x*.phi in
    // the same scope: a signal predicate or constraint, a node that does not
    // vary, or a node placed earlier by order().
    FreezeCtx freeze_ctx(int j) const {
        FreezeCtx ctx = analyze_freeze_ctx(j, f_, parent_);
        if (ctx.sibling < 0) return ctx;
        const int s = ctx.sibling;
        const FormulaType t = f_.nodes[s].type;
        const bool known = t == FormulaType::SIGNAL_PRED || t == FormulaType::SIGNAL_CONSTRAINT ||
                           (owner_[s] != -1 && !varies_[s]) || position_[s] < position_[j];
        return known ? ctx : FreezeCtx{};
    }
    // Evaluate FREEZE node F on the samples t .. n-1 (positions below t are
    // left false), with the current outer frozen values in frozen_.
    void eval_freeze(int F, int t) {
        const FormulaNode& node = f_.nodes[F];
        const Scope& sc = scopes_[F];
        for (int c : sc.constraints) updaters_[c].reset();
        for (int k : sc.invariant) iv_[k] = compute_node(k, f_, tr_, windows_, iv_);
        std::vector<char>& point = points_[F];
        std::fill(point.begin(), point.end(), 0);
        int sib_cursor = sc.ctx.sibling >= 0 ? static_cast<int>(iv_[sc.ctx.sibling].size()) - 1 : -1;

        for (int i = n_ - 1; i >= t; --i) {
            frozen_[node.freeze_var] = tr_.value(node.freeze_component, i);
            for (int c : sc.constraints) {
                updaters_[c].advance(i, f_.nodes[c].frozen_fn(frozen_));
                iv_[c] = updaters_[c].intervals();
            }

            if (sc.ctx.kind == FreezeCtxKind::AND_SIBLING &&
                !contains_decreasing(iv_[sc.ctx.sibling], i, sib_cursor))
                continue;
            if (sc.ctx.kind == FreezeCtxKind::OR_SIBLING &&
                contains_decreasing(iv_[sc.ctx.sibling], i, sib_cursor)) {
                point[i] = 1;
                continue;
            }

            // Nested FREEZE nodes depend on frozen_[var]: re-evaluated on i..n-1.
            for (int k : sc.steps) {
                if (f_.nodes[k].type == FormulaType::FREEZE) eval_freeze(k, i);
                else if (varies_[k]) iv_[k] = compute_node(k, f_, tr_, windows_, iv_);
            }
            point[i] = contains_index(iv_[node.children[0]], i) ? 1 : 0;

            if (sc.ctx.kind == FreezeCtxKind::ALWAYS_UNTIMED && !point[i]) break;
            if (sc.ctx.kind == FreezeCtxKind::EVENTUALLY_UNTIMED && point[i]) {
                for (int r = i - 1; r >= t; --r) point[r] = 1;
                break;
            }
        }
        iv_[F] = points_to_intervals(point);
    }

    const Formula& f_;
    const Trace&   tr_;
    int n_ = 0;
    std::vector<int> parent_, owner_;
    std::vector<IndexIntervals> windows_, iv_;
    std::vector<double> frozen_;
    std::vector<char> varies_;
    std::vector<int> position_;     // index of a node in its evaluation order
    std::vector<Scope> scopes_;
    std::vector<SignalConstraintUpdater> updaters_;
    std::vector<std::vector<char>> points_;
};

}  // namespace

Intervals boolean_monitor_general(const Formula& f, const Trace& trace) {
    return to_time_intervals(GeneralMonitor(f, trace).run(), trace);
}

bool boolean_monitor_general_holds_at_start(const Formula& f, const Trace& trace) {
    return contains_index(GeneralMonitor(f, trace).run(), 0);
}

}  // namespace stl_star
