// naive_monitors.cpp -- reference (naive) monitors.
//
// The four monitors evaluate the STL* semantics point by point, with no
// intervals and no incremental updates, so that they can be checked against
// the definitions by reading them.  They are used to cross-check the fast
// monitors.  Boolean monitors work on truth vectors, robustness monitors on
// vectors of robustness values; apart from that the two pairs have the same
// structure.
#include "stl_star/naive_monitors.h"
#include "monitor_common.h"

#include <algorithm>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace stl_star {

using namespace monitor_detail;

namespace {

// ===========================================================================
// Pointwise evaluation of one formula node
// ===========================================================================

// Robustness of an atom given h = signal_fn(s) and g = frozen_fn(x*) (g = 0
// for a signal predicate): (h - g) - threshold for > and >=, and
// threshold - (h - g) for < and <=.  This is exactly the largest r for which
// the shifted atom of the robustness bisection holds.
double atom_robustness(const FormulaNode& nd, double h, double g) {
    const double d = h - g;
    return (nd.comp_op == CompOp::GT || nd.comp_op == CompOp::GE) ? d - nd.threshold
                                                                  : nd.threshold - d;
}

// The window of sample i -- the samples k with t_k - t_i in [a, b] -- is the
// index range [first, last).  As i increases both ends only move forward, so
// one pass over the trace moves each end at most n times.
struct SlidingWindow {
    const Trace& trace;
    double a, b;
    int first = 0, last = 0;

    void move_to(int i) {
        while (last < trace.length && trace.within_window_end(i, last, b)) ++last;
        if (first < i) first = i;
        while (first < last && !trace.reaches_window_start(i, first, a)) ++first;
    }
};

// F[a,b] phi / G[a,b] phi on truth values, with prefix counts: count[k] is
// the number of samples before k where phi holds, so the window [first, last)
// of sample i contains count[last] - count[first] true samples.
void eventually_always_points(const std::vector<char>& phi, const Trace& trace, double a,
                              double b, bool eventually, std::vector<char>& out) {
    const int n = trace.length;
    std::vector<int> count(n + 1, 0);
    for (int k = 0; k < n; ++k) count[k + 1] = count[k] + (phi[k] ? 1 : 0);

    SlidingWindow window{trace, a, b};
    for (int i = 0; i < n; ++i) {
        window.move_to(i);
        const int size  = window.last - window.first;
        const int trues = count[window.last] - count[window.first];
        out[i] = eventually ? trues > 0        // some sample holds (false if empty)
                            : trues == size;   // every sample holds (true if empty)
    }
}

// F[a,b] phi / G[a,b] phi on robustness values: the sliding max / min of rho
// over the windows, with Lemire's streaming filter.
//
// `candidates` holds the indices of the current window that can still be the
// max (min) of this or a later window: their values decrease (increase) from
// front to back, so the front is the answer.  A new value removes the older
// values it beats, since those leave the window first.  Every index is pushed
// and popped at most once: O(n).  An empty window gives -inf (F) / +inf (G).
void eventually_always_robustness(const std::vector<double>& rho, const Trace& trace,
                                  double a, double b, bool eventually, std::vector<double>& out) {
    const int n = trace.length;
    const double inf = std::numeric_limits<double>::infinity();
    auto beats = [eventually](double x, double y) { return eventually ? x >= y : x <= y; };

    std::deque<int> candidates;
    SlidingWindow window{trace, a, b};
    for (int i = 0; i < n; ++i) {
        const int old_last = window.last;
        window.move_to(i);
        for (int k = old_last; k < window.last; ++k) {        // samples entering the window
            while (!candidates.empty() && beats(rho[k], rho[candidates.back()]))
                candidates.pop_back();
            candidates.push_back(k);
        }
        while (!candidates.empty() && candidates.front() < window.first)   // samples leaving it
            candidates.pop_front();
        out[i] = candidates.empty() ? (eventually ? -inf : inf) : rho[candidates.front()];
    }
}

// Truth values of node j (a signal predicate, or a Boolean / temporal operator)
// at every sample, from its children's truth values.
void eval_pointwise(int j, const Formula& f, const Trace& trace,
                    std::vector<std::vector<char>>& pt) {
    const FormulaNode& nd = f.nodes[j];
    const int n = trace.length;
    std::vector<char>& out = pt[j];
    out.assign(n, 0);
    const double a = nd.interval_lo, b = nd.interval_hi;

    switch (nd.type) {
        case FormulaType::SIGNAL_PRED: {
            std::vector<double> s(trace.num_components);
            for (int i = 0; i < n; ++i) { trace.sample_at(i, s); out[i] = predicate_holds(nd, s); }
            return;
        }
        case FormulaType::NOT: {
            const auto& c = pt[nd.children[0]];
            for (int i = 0; i < n; ++i) out[i] = !c[i];
            return;
        }
        case FormulaType::AND: {
            const auto &l = pt[nd.children[0]], &r = pt[nd.children[1]];
            for (int i = 0; i < n; ++i) out[i] = l[i] && r[i];
            return;
        }
        case FormulaType::OR: {
            const auto &l = pt[nd.children[0]], &r = pt[nd.children[1]];
            for (int i = 0; i < n; ++i) out[i] = l[i] || r[i];
            return;
        }
        case FormulaType::EVENTUALLY_UNTIMED: {   // exists j >= i
            const auto& c = pt[nd.children[0]];
            out[n - 1] = c[n - 1];
            for (int i = n - 2; i >= 0; --i) out[i] = c[i] || out[i + 1];
            return;
        }
        case FormulaType::ALWAYS_UNTIMED: {       // for all j >= i
            const auto& c = pt[nd.children[0]];
            out[n - 1] = c[n - 1];
            for (int i = n - 2; i >= 0; --i) out[i] = c[i] && out[i + 1];
            return;
        }
        case FormulaType::EVENTUALLY:             // phi at some sample of the window
            eventually_always_points(pt[nd.children[0]], trace, a, b, /*eventually=*/true, out);
            return;
        case FormulaType::ALWAYS:                 // phi at every sample of the window
            eventually_always_points(pt[nd.children[0]], trace, a, b, /*eventually=*/false, out);
            return;
        case FormulaType::UNTIL:                  // phi on [i, k), psi at k
        case FormulaType::UNTIL_NS: {             // phi on [i, k], psi at k
            const bool strict = nd.type == FormulaType::UNTIL;
            const auto &phi = pt[nd.children[0]], &psi = pt[nd.children[1]];
            for (int i = 0; i < n; ++i) {
                bool sat = false;
                bool phi_before = true;            // phi holds on [i, k)
                for (int k = i; k < n && trace.within_window_end(i, k, b); ++k) {
                    const bool phi_ok = strict ? phi_before : (phi_before && phi[k]);
                    if (!phi_ok) break;
                    if (trace.reaches_window_start(i, k, a) && psi[k]) { sat = true; break; }
                    phi_before = phi_before && phi[k];
                }
                out[i] = sat;
            }
            return;
        }
        default:
            throw std::logic_error("eval_pointwise: unexpected node type " +
                                   formula_type_to_string(nd.type));
    }
}

// Robustness of node j (a signal predicate, or a Boolean / temporal operator
// other than negation) at every sample, from its children's robustness.
void eval_robustness_pointwise(int j, const Formula& f, const Trace& trace,
                               std::vector<std::vector<double>>& rho) {
    const FormulaNode& nd = f.nodes[j];
    const int n = trace.length;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double>& out = rho[j];
    out.assign(n, 0.0);
    const double a = nd.interval_lo, b = nd.interval_hi;

    switch (nd.type) {
        case FormulaType::SIGNAL_PRED: {
            std::vector<double> s(trace.num_components);
            for (int i = 0; i < n; ++i) {
                trace.sample_at(i, s);
                out[i] = atom_robustness(nd, nd.signal_fn(s), 0.0);
            }
            return;
        }
        case FormulaType::AND: {                  // min
            const auto &l = rho[nd.children[0]], &r = rho[nd.children[1]];
            for (int i = 0; i < n; ++i) out[i] = std::min(l[i], r[i]);
            return;
        }
        case FormulaType::OR: {                   // max
            const auto &l = rho[nd.children[0]], &r = rho[nd.children[1]];
            for (int i = 0; i < n; ++i) out[i] = std::max(l[i], r[i]);
            return;
        }
        case FormulaType::EVENTUALLY_UNTIMED: {   // max over j >= i
            const auto& c = rho[nd.children[0]];
            out[n - 1] = c[n - 1];
            for (int i = n - 2; i >= 0; --i) out[i] = std::max(c[i], out[i + 1]);
            return;
        }
        case FormulaType::ALWAYS_UNTIMED: {       // min over j >= i
            const auto& c = rho[nd.children[0]];
            out[n - 1] = c[n - 1];
            for (int i = n - 2; i >= 0; --i) out[i] = std::min(c[i], out[i + 1]);
            return;
        }
        case FormulaType::EVENTUALLY:             // max over the window
            eventually_always_robustness(rho[nd.children[0]], trace, a, b, /*eventually=*/true, out);
            return;
        case FormulaType::ALWAYS:                 // min over the window
            eventually_always_robustness(rho[nd.children[0]], trace, a, b, /*eventually=*/false, out);
            return;
        case FormulaType::UNTIL:                  // max over k of min(psi at k, phi on [i, k))
        case FormulaType::UNTIL_NS: {             // max over k of min(psi at k, phi on [i, k])
            const bool strict = nd.type == FormulaType::UNTIL;
            const auto &phi = rho[nd.children[0]], &psi = rho[nd.children[1]];
            for (int i = 0; i < n; ++i) {
                double best = -inf;
                double phi_min = inf;              // min of phi on [i, k)
                for (int k = i; k < n && trace.within_window_end(i, k, b); ++k) {
                    if (!strict) phi_min = std::min(phi_min, phi[k]);
                    if (trace.reaches_window_start(i, k, a))
                        best = std::max(best, std::min(psi[k], phi_min));
                    if (strict) phi_min = std::min(phi_min, phi[k]);
                }
                out[i] = best;
            }
            return;
        }
        default:
            throw std::logic_error("eval_robustness_pointwise: unexpected node type " +
                                   formula_type_to_string(nd.type));
    }
}

// ===========================================================================
// Boolean, STL*_1
// ===========================================================================

std::vector<char> run_boolean_one(const Formula& f, const Trace& trace) {
    check_input(f, trace, /*require_fragment_one=*/true);
    const int N = static_cast<int>(f.nodes.size());
    const int n = trace.length;
    const std::vector<int> owner = freeze_owners(f);

    std::vector<std::vector<char>> pt(N, std::vector<char>(n, 0));
    std::vector<double> frozen(f.num_freeze_vars, 0.0);
    std::vector<double> s(trace.num_components);

    for (int j = 0; j < N; ++j) {
        const FormulaNode& node = f.nodes[j];
        if (node.type == FormulaType::SIGNAL_PRED) { eval_pointwise(j, f, trace, pt); continue; }
        if (node.type == FormulaType::SIGNAL_CONSTRAINT) continue;   // done by its FREEZE
        if (node.type != FormulaType::FREEZE) {
            if (owner[j] == -1) eval_pointwise(j, f, trace, pt);
            continue;
        }

        // FREEZE node j: for every i, freeze s_x(i) and evaluate the scope.
        std::vector<int> scope;   // ascending: children before parents
        for (int k = 0; k < j; ++k) if (owner[k] == j) scope.push_back(k);

        std::vector<char> result(n, 0);
        for (int i = n - 1; i >= 0; --i) {
            frozen[node.freeze_var] = trace.value(node.freeze_component, i);
            for (int k : scope) {
                const FormulaNode& nd = f.nodes[k];
                if (nd.type == FormulaType::SIGNAL_PRED || nd.type == FormulaType::FREEZE)
                    continue;   // independent of this frozen value: already computed
                if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
                    const double g = nd.frozen_fn(frozen);
                    for (int l = 0; l < n; ++l) {
                        trace.sample_at(l, s);
                        pt[k][l] = constraint_holds(nd, nd.signal_fn(s), g);
                    }
                } else {
                    eval_pointwise(k, f, trace, pt);
                }
            }
            result[i] = pt[node.children[0]][i];
        }
        pt[j] = std::move(result);
    }
    return pt[f.root_index()];
}

// ===========================================================================
// Boolean, general STL*
// ===========================================================================

class BooleanGeneral {
public:
    BooleanGeneral(const Formula& f, const Trace& trace) : f_(f), tr_(trace) {}

    std::vector<char> run() {
        check_input(f_, tr_, /*require_fragment_one=*/false);
        const int N = static_cast<int>(f_.nodes.size());
        owner_ = freeze_owners(f_);
        pt_.assign(N, std::vector<char>(tr_.length, 0));
        frozen_.assign(f_.num_freeze_vars, 0.0);
        s_.assign(tr_.num_components, 0.0);

        for (int j = 0; j < N; ++j) {
            if (owner_[j] != -1) continue;   // evaluated by the enclosing FREEZE
            eval_node(j);
        }
        return pt_[f_.root_index()];
    }

private:
    // Evaluate node j on every sample with the current frozen values.
    void eval_node(int j) {
        const FormulaNode& nd = f_.nodes[j];
        if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
            const double g = nd.frozen_fn(frozen_);
            for (int l = 0; l < tr_.length; ++l) {
                tr_.sample_at(l, s_);
                pt_[j][l] = constraint_holds(nd, nd.signal_fn(s_), g);
            }
        } else if (nd.type == FormulaType::FREEZE) {
            std::vector<int> scope;              // ascending: children first
            for (int k = 0; k < j; ++k) if (owner_[k] == j) scope.push_back(k);
            std::vector<char> result(tr_.length, 0);
            for (int i = 0; i < tr_.length; ++i) {
                frozen_[nd.freeze_var] = tr_.value(nd.freeze_component, i);
                for (int k : scope) eval_node(k);
                result[i] = pt_[nd.children[0]][i];
            }
            pt_[j] = std::move(result);
        } else {
            eval_pointwise(j, f_, tr_, pt_);
        }
    }

    const Formula& f_;
    const Trace&   tr_;
    std::vector<int> owner_;
    std::vector<std::vector<char>> pt_;
    std::vector<double> frozen_, s_;
};

// ===========================================================================
// Robustness, STL*_1
// ===========================================================================

// Same structure as run_boolean_one, on robustness values.
double run_robustness_one(const Formula& f, const Trace& trace) {
    check_input(f, trace, /*require_fragment_one=*/true);
    check_negation_free(f);
    const int N = static_cast<int>(f.nodes.size());
    const int n = trace.length;
    const std::vector<int> owner = freeze_owners(f);

    std::vector<std::vector<double>> rho(N, std::vector<double>(n, 0.0));
    std::vector<double> frozen(f.num_freeze_vars, 0.0);
    std::vector<double> s(trace.num_components);

    for (int j = 0; j < N; ++j) {
        const FormulaNode& node = f.nodes[j];
        if (node.type == FormulaType::SIGNAL_PRED) { eval_robustness_pointwise(j, f, trace, rho); continue; }
        if (node.type == FormulaType::SIGNAL_CONSTRAINT) continue;   // done by its FREEZE
        if (node.type != FormulaType::FREEZE) {
            if (owner[j] == -1) eval_robustness_pointwise(j, f, trace, rho);
            continue;
        }

        // FREEZE node j: for every i, freeze s_x(i) and evaluate the scope.
        std::vector<int> scope;   // ascending: children before parents
        for (int k = 0; k < j; ++k) if (owner[k] == j) scope.push_back(k);

        std::vector<double> result(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            frozen[node.freeze_var] = trace.value(node.freeze_component, i);
            for (int k : scope) {
                const FormulaNode& nd = f.nodes[k];
                if (nd.type == FormulaType::SIGNAL_PRED || nd.type == FormulaType::FREEZE)
                    continue;   // independent of this frozen value: already computed
                if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
                    const double g = nd.frozen_fn(frozen);
                    for (int l = 0; l < n; ++l) {
                        trace.sample_at(l, s);
                        rho[k][l] = atom_robustness(nd, nd.signal_fn(s), g);
                    }
                } else {
                    eval_robustness_pointwise(k, f, trace, rho);
                }
            }
            result[i] = rho[node.children[0]][i];
        }
        rho[j] = std::move(result);
    }
    return rho[f.root_index()][0];
}

// ===========================================================================
// Robustness, general STL*
// ===========================================================================

// Same structure as BooleanGeneral, on robustness values.
class RobustnessGeneral {
public:
    RobustnessGeneral(const Formula& f, const Trace& trace) : f_(f), tr_(trace) {}

    double run() {
        check_input(f_, tr_, /*require_fragment_one=*/false);
        check_negation_free(f_);
        const int N = static_cast<int>(f_.nodes.size());
        owner_ = freeze_owners(f_);
        rho_.assign(N, std::vector<double>(tr_.length, 0.0));
        frozen_.assign(f_.num_freeze_vars, 0.0);
        s_.assign(tr_.num_components, 0.0);

        for (int j = 0; j < N; ++j) {
            if (owner_[j] != -1) continue;   // evaluated by the enclosing FREEZE
            eval_node(j);
        }
        return rho_[f_.root_index()][0];
    }

private:
    // Evaluate node j on every sample with the current frozen values.
    void eval_node(int j) {
        const FormulaNode& nd = f_.nodes[j];
        if (nd.type == FormulaType::SIGNAL_CONSTRAINT) {
            const double g = nd.frozen_fn(frozen_);
            for (int l = 0; l < tr_.length; ++l) {
                tr_.sample_at(l, s_);
                rho_[j][l] = atom_robustness(nd, nd.signal_fn(s_), g);
            }
        } else if (nd.type == FormulaType::FREEZE) {
            std::vector<int> scope;              // ascending: children first
            for (int k = 0; k < j; ++k) if (owner_[k] == j) scope.push_back(k);
            std::vector<double> result(tr_.length, 0.0);
            for (int i = 0; i < tr_.length; ++i) {
                frozen_[nd.freeze_var] = tr_.value(nd.freeze_component, i);
                for (int k : scope) eval_node(k);
                result[i] = rho_[nd.children[0]][i];
            }
            rho_[j] = std::move(result);
        } else {
            eval_robustness_pointwise(j, f_, tr_, rho_);
        }
    }

    const Formula& f_;
    const Trace&   tr_;
    std::vector<int> owner_;
    std::vector<std::vector<double>> rho_;
    std::vector<double> frozen_, s_;
};

}  // namespace

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------
Intervals boolean_monitor_one_naive(const Formula& f, const Trace& trace) {
    return to_time_intervals(points_to_intervals(run_boolean_one(f, trace)), trace);
}

bool boolean_monitor_one_naive_holds_at_start(const Formula& f, const Trace& trace) {
    return run_boolean_one(f, trace)[0] != 0;
}

Intervals boolean_monitor_general_naive(const Formula& f, const Trace& trace) {
    return to_time_intervals(points_to_intervals(BooleanGeneral(f, trace).run()), trace);
}

bool boolean_monitor_general_naive_holds_at_start(const Formula& f, const Trace& trace) {
    return BooleanGeneral(f, trace).run()[0] != 0;
}

double robustness_monitor_one_naive(const Formula& f, const Trace& trace) {
    return run_robustness_one(f, trace);
}

double robustness_monitor_general_naive(const Formula& f, const Trace& trace) {
    return RobustnessGeneral(f, trace).run();
}

}  // namespace stl_star
