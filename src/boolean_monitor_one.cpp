#include "stl_star/boolean_monitor_one.h"
#include "monitor_common.h"
#include "stl_star/signal_constraint_updater.h"

#include <vector>

namespace stl_star {

using namespace monitor_detail;

namespace {

// STL*_1: every constraint reads only the variable of its nearest FREEZE, so
// the scopes are independent and every FREEZE is evaluated exactly once.
//
// A node is VARIANT when it lies in the scope of a FREEZE and its value
// depends on that FREEZE's frozen value (its subtree has a constraint of the
// scope); variant nodes are recomputed at every step of the FREEZE loop.  All
// other nodes have one fixed value and are computed once, on demand (eval).
class MonitorOne {
public:
    MonitorOne(const Formula& f, const Trace& trace) : f_(f), tr_(trace) {}

    IndexIntervals run() {
        check_input(f_, tr_, /*require_fragment_one=*/true);
        const int N = static_cast<int>(f_.nodes.size());
        n_       = tr_.length;
        parent_  = parents_of(f_);
        owner_   = freeze_owners(f_);
        windows_ = compute_window_sets(f_, tr_);
        iv_.assign(N, {});
        done_.assign(N, 0);
        busy_.assign(N, 0);
        variant_.assign(N, 0);
        frozen_.assign(f_.num_freeze_vars, 0.0);

        for (int k = 0; k < N; ++k) {
            const FormulaNode& nd = f_.nodes[k];
            if (nd.type == FormulaType::SIGNAL_PRED) {
                iv_[k]   = signal_pred_intervals(nd, tr_);   // never depends on frozen values
                done_[k] = 1;
            } else if (owner_[k] != -1) {
                if (nd.type == FormulaType::SIGNAL_CONSTRAINT) variant_[k] = 1;
                else if (nd.type != FormulaType::FREEZE)         // a nested FREEZE is independent
                    for (int c : nd.children) variant_[k] = variant_[k] || variant_[c];
            }
        }
        eval(f_.root_index());
        return iv_[f_.root_index()];
    }

private:
    // Compute a node that is not variant (and its subtree), once.
    void eval(int k) {
        if (done_[k]) return;
        const FormulaNode& nd = f_.nodes[k];
        if (nd.type == FormulaType::FREEZE) {
            run_freeze(k);
        } else {
            for (int c : nd.children) eval(c);
            iv_[k] = compute_node(k, f_, tr_, windows_, iv_);
        }
        done_[k] = 1;
    }

    // ---- FREEZE node j: x*.phi ----------------------------------------------
    void run_freeze(int j) {
        const FormulaNode& node = f_.nodes[j];

        // Heuristics from the parent (see monitor_common.h).  The sibling phi1
        // of (phi1 && x*.phi) / (phi1 || x*.phi) must be known before the loop:
        // that is the case unless it depends on an enclosing frozen value, or
        // is itself a FREEZE whose evaluation is waiting for this one
        // ((x*.phi && y*.psi): the FREEZE started first evaluates the other one
        // and uses it; the other one cannot use the first).
        busy_[j] = 1;
        FreezeCtx ctx = analyze_freeze_ctx(j, f_, parent_);
        if (ctx.sibling >= 0) {
            if (variant_[ctx.sibling] || busy_[ctx.sibling]) ctx = {};
            else                                             eval(ctx.sibling);
        }

        // Scope = nodes whose nearest enclosing FREEZE is j (ascending order,
        // so children come before parents).  Non-variant scope nodes (built
        // from signal predicates and nested FREEZE nodes) are computed once.
        std::vector<int> constraints, compute;
        for (int k = 0; k < j; ++k) {
            if (owner_[k] != j) continue;
            const FormulaType t = f_.nodes[k].type;
            if (t == FormulaType::SIGNAL_CONSTRAINT) constraints.push_back(k);
            else if (variant_[k])                    compute.push_back(k);
            else                                     eval(k);
        }
        std::vector<SignalConstraintUpdater> updaters(constraints.size());
        for (size_t k = 0; k < constraints.size(); ++k)
            updaters[k].init(f_.nodes[constraints[k]], tr_);

        const int child = node.children[0];
        int sib_cursor = ctx.sibling >= 0 ? static_cast<int>(iv_[ctx.sibling].size()) - 1 : -1;
        std::vector<char> point(n_, 0);

        for (int i = n_ - 1; i >= 0; --i) {
            frozen_[node.freeze_var] = tr_.value(node.freeze_component, i);

            // The updaters must see every step to keep their state consistent,
            // even when the rest of the scope is skipped below.
            for (size_t k = 0; k < constraints.size(); ++k) {
                updaters[k].advance(i, f_.nodes[constraints[k]].frozen_fn(frozen_));
                iv_[constraints[k]] = updaters[k].intervals();
            }

            if (ctx.kind == FreezeCtxKind::AND_SIBLING &&
                !contains_decreasing(iv_[ctx.sibling], i, sib_cursor))
                continue;                       // parent AND is false at i anyway
            if (ctx.kind == FreezeCtxKind::OR_SIBLING &&
                contains_decreasing(iv_[ctx.sibling], i, sib_cursor)) {
                point[i] = 1;                   // parent OR is true at i anyway
                continue;
            }

            for (int k : compute) iv_[k] = compute_node(k, f_, tr_, windows_, iv_);
            point[i] = contains_index(iv_[child], i) ? 1 : 0;

            if (ctx.kind == FreezeCtxKind::ALWAYS_UNTIMED && !point[i])
                break;                          // G(x*.phi) is false at every j <= i
            if (ctx.kind == FreezeCtxKind::EVENTUALLY_UNTIMED && point[i]) {
                for (int r = i - 1; r >= 0; --r) point[r] = 1;   // F(x*.phi) true at j <= i
                break;
            }
        }
        iv_[j] = points_to_intervals(point);
    }

    const Formula& f_;
    const Trace&   tr_;
    int n_ = 0;
    std::vector<int> parent_, owner_;
    std::vector<IndexIntervals> windows_, iv_;
    std::vector<char> done_, busy_, variant_;   // busy_: FREEZE loop being prepared
    std::vector<double> frozen_;
};

}  // namespace

Intervals boolean_monitor_one(const Formula& f, const Trace& trace) {
    return to_time_intervals(MonitorOne(f, trace).run(), trace);
}

bool boolean_monitor_one_holds_at_start(const Formula& f, const Trace& trace) {
    return contains_index(MonitorOne(f, trace).run(), 0);
}

}  // namespace stl_star
