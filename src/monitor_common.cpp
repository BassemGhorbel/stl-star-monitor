#include "monitor_common.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace stl_star {

namespace monitor_detail {

std::vector<int> parents_of(const Formula& f) {
    std::vector<int> parent(f.nodes.size(), -1);
    for (int i = 0; i < static_cast<int>(f.nodes.size()); ++i)
        for (int c : f.nodes[i].children) parent[c] = i;
    return parent;
}

std::vector<int> freeze_owners(const Formula& f) {
    const int N = static_cast<int>(f.nodes.size());
    std::vector<int> owner(N, -1);
    // Parents have larger indices, so a top-down pass (high to low) sees the
    // owner of the parent before the child.
    for (int i = N - 1; i >= 0; --i) {
        const FormulaNode& nd = f.nodes[i];
        int own = (nd.type == FormulaType::FREEZE) ? i : owner[i];
        for (int c : nd.children) owner[c] = own;
    }
    return owner;
}

void check_input(const Formula& f, const Trace& trace, bool require_fragment_one) {
    f.validate();
    // Parsed formulas: one signal per letter, in alphabetical order.
    if (f.num_signals > 0 && f.num_signals != trace.num_components) {
        std::string letters;
        for (char c : f.signal_names) letters += std::string(letters.empty() ? "" : ", ") + c;
        throw std::invalid_argument(
            "The formula uses " + std::to_string(f.num_signals) + " signal(s)" +
            (letters.empty() ? "" : " (" + letters + ")") + " but the trace has " +
            std::to_string(trace.num_components) +
            "; give one signal per letter, in alphabetical order of the letters");
    }
    for (const FormulaNode& nd : f.nodes)
        if (nd.type == FormulaType::FREEZE && nd.freeze_component >= trace.num_components)
            throw std::invalid_argument("A freeze operator reads signal component " +
                                        std::to_string(nd.freeze_component) +
                                        ", which the trace does not have");
    std::string why;
    if (require_fragment_one && !is_fragment_one(f, &why))
        throw std::invalid_argument("The formula is not in STL*_1 (" + why +
                                    "); use the general STL* monitor instead "
                                    "(monitor_general / robustness_general)");
}

std::vector<IndexIntervals> compute_window_sets(const Formula& f, const Trace& trace) {
    std::vector<IndexIntervals> w(f.nodes.size());
    for (size_t j = 0; j < f.nodes.size(); ++j) {
        const FormulaNode& nd = f.nodes[j];
        switch (nd.type) {
            case FormulaType::EVENTUALLY:
            case FormulaType::ALWAYS:
            case FormulaType::UNTIL:
            case FormulaType::UNTIL_NS:
                w[j] = window_set(trace, nd.interval_lo, nd.interval_hi);
                break;
            default:
                break;
        }
    }
    return w;
}

IndexIntervals signal_pred_intervals(const FormulaNode& nd, const Trace& trace) {
    std::vector<char> pt(trace.length);
    std::vector<double> s(trace.num_components);
    for (int i = 0; i < trace.length; ++i) {
        trace.sample_at(i, s);
        pt[i] = predicate_holds(nd, s) ? 1 : 0;
    }
    return points_to_intervals(pt);
}

IndexIntervals compute_node(int j, const Formula& f, const Trace& trace,
                           const std::vector<IndexIntervals>& windows,
                           const std::vector<IndexIntervals>& iv) {
    const FormulaNode& nd = f.nodes[j];
    const int n = trace.length;
    switch (nd.type) {
        case FormulaType::NOT:
            return complement_of(iv[nd.children[0]], n);
        case FormulaType::AND:
            return intersection_of(iv[nd.children[0]], iv[nd.children[1]]);
        case FormulaType::OR:
            return union_of(iv[nd.children[0]], iv[nd.children[1]]);
        case FormulaType::EVENTUALLY:
            return eventually_of(iv[nd.children[0]], nd.interval_lo, nd.interval_hi,
                                windows[j], trace);
        case FormulaType::ALWAYS:
            return always_of(iv[nd.children[0]], nd.interval_lo, nd.interval_hi,
                            windows[j], trace);
        case FormulaType::UNTIL:
        case FormulaType::UNTIL_NS:
            return until_of(iv[nd.children[0]], iv[nd.children[1]], nd.interval_lo,
                           nd.interval_hi, nd.type == FormulaType::UNTIL, windows[j], trace);
        case FormulaType::EVENTUALLY_UNTIMED:
            return eventually_untimed_of(iv[nd.children[0]]);
        case FormulaType::ALWAYS_UNTIMED:
            return always_untimed_of(iv[nd.children[0]], n);
        default:
            // Atomic propositions and FREEZE nodes are produced by the monitors.
            throw std::logic_error("compute_node: unexpected node type " +
                                   formula_type_to_string(nd.type));
    }
}

void check_negation_free(const Formula& f) {
    for (const FormulaNode& nd : f.nodes)
        if (nd.type == FormulaType::NOT)
            throw std::invalid_argument("Robustness monitoring requires a negation-free formula");
}

FreezeCtx analyze_freeze_ctx(int freeze_idx, const Formula& f, const std::vector<int>& parent) {
    const int p = parent[freeze_idx];
    if (p < 0) return {};
    const FormulaNode& pn = f.nodes[p];
    if (pn.type == FormulaType::ALWAYS_UNTIMED)     return {FreezeCtxKind::ALWAYS_UNTIMED, -1};
    if (pn.type == FormulaType::EVENTUALLY_UNTIMED) return {FreezeCtxKind::EVENTUALLY_UNTIMED, -1};
    if (pn.type == FormulaType::AND || pn.type == FormulaType::OR) {
        const int sib = pn.children[0] == freeze_idx ? pn.children[1] : pn.children[0];
        return {pn.type == FormulaType::AND ? FreezeCtxKind::AND_SIBLING
                                            : FreezeCtxKind::OR_SIBLING, sib};
    }
    return {};
}

}  // namespace monitor_detail

}  // namespace stl_star
