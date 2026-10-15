#include "stl_star/formula_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace stl_star {

namespace {

// A parse error that remembers where it happened, so that when two parse
// alternatives fail we can report the one that got further.  Semantic errors
// (unbound or re-bound variables, unsupported constraints, bad windows) are
// final: they are reported even if another alternative could be tried.
struct ParseError : std::runtime_error {
    int  pos;
    bool final;
    ParseError(const std::string& msg, int p, bool is_final = false)
        : std::runtime_error(msg + " at position " + std::to_string(p)), pos(p), final(is_final) {}
};

// ============================================================================
// Lexer
// ============================================================================
enum class TK {
    ALWAYS, EVENTUALLY, UNTIL, UNTIL_NS, NOT, AND, OR, INF,
    PLUS, MINUS, STAR, SLASH, GEQ, LEQ, GT, LT,
    LPAREN, RPAREN, LBRACKET, RBRACKET, COMMA, DOT,
    NUMBER, SIGNAL, FREEZE_VAR, END
};

struct Token {
    TK     type;
    int    pos           = 0;
    int    end           = 0;      // one past the last character
    double num_val       = 0.0;
    char   letter        = '\0';   // SIGNAL / FREEZE_VAR
    int    subscript     = 0;      // FREEZE_VAR
    bool   has_subscript = false;  // FREEZE_VAR
};

const std::map<std::string, TK> KEYWORDS = {
    {"always", TK::ALWAYS},   {"box", TK::ALWAYS},
    {"eventually", TK::EVENTUALLY}, {"diamond", TK::EVENTUALLY},
    {"until", TK::UNTIL},     {"uns", TK::UNTIL_NS}, {"untilns", TK::UNTIL_NS},
    {"not", TK::NOT}, {"and", TK::AND}, {"or", TK::OR}, {"inf", TK::INF},
};

class Lexer {
public:
    explicit Lexer(const std::string& src) : src_(src) {}

    std::vector<Token> tokenize() {
        std::vector<Token> out;
        for (;;) {
            while (pos_ < size() && std::isspace(static_cast<unsigned char>(src_[pos_]))) ++pos_;
            if (pos_ >= size()) { out.push_back({TK::END, pos_, pos_}); return out; }
            out.push_back(read_one());
            out.back().end = pos_;
        }
    }

private:
    const std::string& src_;
    int pos_ = 0;

    int  size() const { return static_cast<int>(src_.size()); }
    char at(int p) const { return p < size() ? src_[p] : '\0'; }
    static bool is_digit(char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; }

    Token read_one() {
        const int start = pos_;
        const char c = src_[pos_];
        if (is_digit(c)) return read_number();
        if (std::isalpha(static_cast<unsigned char>(c))) return read_word();

        const char c2 = at(pos_ + 1);
        if (c == '>' && c2 == '=') { pos_ += 2; return {TK::GEQ, start}; }
        if (c == '<' && c2 == '=') { pos_ += 2; return {TK::LEQ, start}; }
        if (c == '&' && c2 == '&') { pos_ += 2; return {TK::AND, start}; }
        if (c == '|' && c2 == '|') { pos_ += 2; return {TK::OR,  start}; }
        ++pos_;
        switch (c) {
            case '>': return {TK::GT, start};       case '<': return {TK::LT, start};
            case '!': return {TK::NOT, start};      case '+': return {TK::PLUS, start};
            case '-': return {TK::MINUS, start};    case '*': return {TK::STAR, start};
            case '/': return {TK::SLASH, start};    case '(': return {TK::LPAREN, start};
            case ')': return {TK::RPAREN, start};   case '[': return {TK::LBRACKET, start};
            case ']': return {TK::RBRACKET, start}; case ',': return {TK::COMMA, start};
            case '.': return {TK::DOT, start};
            default: break;
        }
        throw ParseError(std::string("Unexpected character '") + c + "'", start);
    }

    // digits [. digits] [(e|E) [+|-] digits]
    Token read_number() {
        const int start = pos_;
        while (is_digit(at(pos_))) ++pos_;
        if (at(pos_) == '.') { ++pos_; while (is_digit(at(pos_))) ++pos_; }
        if ((at(pos_) == 'e' || at(pos_) == 'E') &&
            (is_digit(at(pos_ + 1)) ||
             ((at(pos_ + 1) == '+' || at(pos_ + 1) == '-') && is_digit(at(pos_ + 2))))) {
            pos_ += 2;
            while (is_digit(at(pos_))) ++pos_;
        }
        Token t{TK::NUMBER, start};
        // strtod: an underflow gives 0 or a subnormal value (accepted), an
        // overflow gives inf (rejected).
        t.num_val = std::strtod(src_.substr(start, pos_ - start).c_str(), nullptr);
        if (!std::isfinite(t.num_val)) throw ParseError("Number out of range", start);
        return t;
    }

    Token read_word() {
        const int start = pos_;
        while (std::isalpha(static_cast<unsigned char>(at(pos_)))) ++pos_;
        std::string word = src_.substr(start, pos_ - start);

        if (word.size() == 1 && std::isupper(static_cast<unsigned char>(word[0]))) {
            switch (word[0]) {
                case 'G': return {TK::ALWAYS, start};
                case 'F': return {TK::EVENTUALLY, start};
                case 'U': return {TK::UNTIL, start};
                default: throw ParseError("Unknown operator '" + word + "'", start);
            }
        }
        if (word.size() > 1) {
            std::string lower = word;
            for (char& ch : lower) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
            auto it = KEYWORDS.find(lower);
            if (it == KEYWORDS.end())
                throw ParseError("Unknown identifier '" + word +
                                 "' (signals are single lowercase letters)", start);
            return {it->second, start};
        }

        // Single lowercase letter: a signal, or a freeze variable if followed by '*'.
        Token t{TK::SIGNAL, start};
        t.letter = word[0];
        if (at(pos_) == '*') {
            ++pos_;
            t.type = TK::FREEZE_VAR;
            if (at(pos_) == '_') {
                if (!is_digit(at(pos_ + 1)))
                    throw ParseError("Expected digits after '_'", pos_);
                ++pos_;
                const int ds = pos_;
                while (is_digit(at(pos_))) ++pos_;
                if (pos_ - ds > 9) throw ParseError("Subscript too long", ds);
                t.subscript    = std::stoi(src_.substr(ds, pos_ - ds));
                t.has_subscript = true;
            }
        }
        return t;
    }
};

// ============================================================================
// Arithmetic expressions
// ============================================================================
enum class ExprKind { Num, Signal, FreezeVar, Neg, Add, Sub, Mul, Div };

struct Expr;
using ExprP = std::shared_ptr<const Expr>;

struct Expr {
    ExprKind kind;
    double   num   = 0.0;  // Num
    int      index = -1;   // Signal: component, FreezeVar: variable
    ExprP    l, r;         // Neg uses l only
};

ExprP mk(ExprKind k, ExprP l = nullptr, ExprP r = nullptr) {
    auto e = std::make_shared<Expr>();
    e->kind = k; e->l = std::move(l); e->r = std::move(r);
    return e;
}
ExprP mk_num(double v)              { auto e = std::make_shared<Expr>(); e->kind = ExprKind::Num; e->num = v; return e; }
ExprP mk_leaf(ExprKind k, int idx)  { auto e = std::make_shared<Expr>(); e->kind = k; e->index = idx; return e; }

bool has_kind(const ExprP& e, ExprKind k) {
    if (!e) return false;
    if (e->kind == k) return true;
    return has_kind(e->l, k) || has_kind(e->r, k);
}
bool has_signal(const ExprP& e) { return has_kind(e, ExprKind::Signal); }
bool has_frozen(const ExprP& e) { return has_kind(e, ExprKind::FreezeVar); }

void collect_vars(const ExprP& e, std::set<int>& out) {
    if (!e) return;
    if (e->kind == ExprKind::FreezeVar) out.insert(e->index);
    collect_vars(e->l, out);
    collect_vars(e->r, out);
}

// Builders that treat nullptr as the constant 0.
ExprP add(ExprP a, ExprP b) { if (!a) return b; if (!b) return a; return mk(ExprKind::Add, a, b); }
ExprP sub(ExprP a, ExprP b) { if (!b) return a; if (!a) return mk(ExprKind::Neg, b); return mk(ExprKind::Sub, a, b); }
ExprP neg(ExprP a)          { return a ? mk(ExprKind::Neg, a) : nullptr; }
ExprP mul(ExprP a, ExprP c) { return a ? mk(ExprKind::Mul, a, c) : nullptr; }
ExprP divide(ExprP a, ExprP c) { return a ? mk(ExprKind::Div, a, c) : nullptr; }

// Split e into  e = sig + frz  where sig reads only signals and constants and
// frz reads only freeze variables and constants (nullptr = 0).
struct Parts { ExprP sig, frz; };

Parts separate(const ExprP& e, int pos) {
    if (!has_frozen(e)) return {e, nullptr};
    if (!has_signal(e)) return {nullptr, e};
    // e mixes both kinds of leaves.
    switch (e->kind) {
        case ExprKind::Neg: { Parts p = separate(e->l, pos); return {neg(p.sig), neg(p.frz)}; }
        case ExprKind::Add: { Parts a = separate(e->l, pos), b = separate(e->r, pos);
                              return {add(a.sig, b.sig), add(a.frz, b.frz)}; }
        case ExprKind::Sub: { Parts a = separate(e->l, pos), b = separate(e->r, pos);
                              return {sub(a.sig, b.sig), sub(a.frz, b.frz)}; }
        case ExprKind::Mul: {
            const bool l_const = !has_signal(e->l) && !has_frozen(e->l);
            const bool r_const = !has_signal(e->r) && !has_frozen(e->r);
            if (l_const) { Parts p = separate(e->r, pos); return {mul(p.sig, e->l), mul(p.frz, e->l)}; }
            if (r_const) { Parts p = separate(e->l, pos); return {mul(p.sig, e->r), mul(p.frz, e->r)}; }
            break;
        }
        case ExprKind::Div: {
            const bool r_const = !has_signal(e->r) && !has_frozen(e->r);
            if (r_const) { Parts p = separate(e->l, pos); return {divide(p.sig, e->r), divide(p.frz, e->r)}; }
            break;
        }
        default: break;
    }
    throw ParseError("Unsupported constraint: a product or quotient may not combine "
                     "signals with freeze variables (the constraint must have the "
                     "form f(signals) ~ g(freeze variables))", pos, /*final=*/true);
}

// ---- Compilation to closures ----------------------------------------------
// `leaf_kind` is the only kind of leaf e may contain (guaranteed by separate());
// every leaf reads x[index] of the vector passed to the closure.
using VecFn = std::function<double(const std::vector<double>&)>;

VecFn compile(const ExprP& e, ExprKind leaf_kind) {
    if (!e) return [](const std::vector<double>&) { return 0.0; };
    switch (e->kind) {
        case ExprKind::Num: { double v = e->num; return [v](const std::vector<double>&) { return v; }; }
        case ExprKind::Signal:
        case ExprKind::FreezeVar: {
            int idx = e->index;   // the separation guarantees e->kind == leaf_kind
            (void)leaf_kind;
            return [idx](const std::vector<double>& x) { return x[idx]; };
        }
        case ExprKind::Neg: { auto f = compile(e->l, leaf_kind); return [f](const std::vector<double>& x) { return -f(x); }; }
        default: break;
    }
    auto f = compile(e->l, leaf_kind), g = compile(e->r, leaf_kind);
    switch (e->kind) {
        case ExprKind::Add: return [f, g](const std::vector<double>& x) { return f(x) + g(x); };
        case ExprKind::Sub: return [f, g](const std::vector<double>& x) { return f(x) - g(x); };
        case ExprKind::Mul: return [f, g](const std::vector<double>& x) { return f(x) * g(x); };
        default:            return [f, g](const std::vector<double>& x) { return f(x) / g(x); };
    }
}

// ============================================================================
// Parser
// ============================================================================
struct Binding { char letter; int subscript; bool has_subscript; int var; };

class Parser {
public:
    // The signal letters of the formula, sorted, name the components:
    // the k-th letter in alphabetical order is component k.
    explicit Parser(const std::vector<Token>& tokens) : toks_(tokens) {
        std::set<char> letters;
        for (const Token& t : tokens)
            if (t.type == TK::SIGNAL || t.type == TK::FREEZE_VAR) letters.insert(t.letter);
        signal_names_.assign(letters.begin(), letters.end());    }

    Formula parse() {
        parse_subformula();
        if (cur().type == TK::AND || cur().type == TK::OR || cur().type == TK::UNTIL ||
            cur().type == TK::UNTIL_NS)
            throw ParseError("Unexpected token after the end of the formula (binary operators "
                             "must be enclosed in parentheses, e.g. \"(phi && psi)\")", cur().pos);
        if (cur().type != TK::END)
            throw ParseError("Unexpected token after the end of the formula", cur().pos);
        f_.num_freeze_vars = next_var_;
        f_.num_signals     = static_cast<int>(signal_names_.size());
        f_.signal_names    = signal_names_;
        f_.validate();
        return std::move(f_);
    }

private:
    const std::vector<Token>& toks_;
    int pos_ = 0;
    Formula f_;
    std::vector<Binding> scope_;     // active freeze bindings, innermost last
    int next_var_ = 0;
    std::string signal_names_;       // the signal letters, sorted

    const Token& cur() const { return toks_[pos_]; }
    const Token& peek(int k) const {
        int i = pos_ + k;
        return toks_[std::min(i, static_cast<int>(toks_.size()) - 1)];
    }
    bool at(TK t) const { return cur().type == t; }
    Token take() { return toks_[pos_++]; }
    void expect(TK t, const std::string& what) {
        if (!at(t)) throw ParseError("Expected " + what, cur().pos);
        ++pos_;
    }

    int signal_index(char letter, int pos) {
        const size_t idx = signal_names_.find(letter);
        if (letter < 'a' || letter > 'z' || idx == std::string::npos)
            throw ParseError("Invalid signal name", pos);
        return static_cast<int>(idx);
    }

    // ---- formula := freeze | G | F | ! | ( ... ) | atom ------------------
    int parse_subformula() {
        if (at(TK::FREEZE_VAR) && peek(1).type == TK::DOT) return parse_freeze();
        if (at(TK::ALWAYS))     return parse_unary_temporal(true);
        if (at(TK::EVENTUALLY)) return parse_unary_temporal(false);
        if (at(TK::NOT)) { take(); int c = parse_subformula(); return f_.add_node(FormulaNode::make_not(c)); }
        if (at(TK::LPAREN)) return parse_paren_or_atom();
        return parse_atom();
    }

    int parse_freeze() {
        Token v = take();
        expect(TK::DOT, "'.' after the freeze variable");
        if (v.has_subscript)
            for (const Binding& b : scope_)
                if (b.letter == v.letter && b.has_subscript && b.subscript == v.subscript)
                    throw ParseError(std::string("Freeze variable '") + v.letter + "*_" +
                                     std::to_string(v.subscript) + "' is already bound", v.pos,
                                     /*final=*/true);
        const int comp = signal_index(v.letter, v.pos);
        const int var  = next_var_++;
        scope_.push_back({v.letter, v.subscript, v.has_subscript, var});
        const int child = parse_subformula();
        scope_.pop_back();
        return f_.add_node(FormulaNode::make_freeze(var, comp, child));
    }

    int parse_unary_temporal(bool always) {
        take();
        if (at(TK::LBRACKET)) {
            auto [lo, hi] = parse_window();
            int c = parse_subformula();
            return f_.add_node(always ? FormulaNode::make_always(lo, hi, c)
                                     : FormulaNode::make_eventually(lo, hi, c));
        }
        int c = parse_subformula();
        return f_.add_node(always ? FormulaNode::make_always_untimed(c)
                                 : FormulaNode::make_eventually_untimed(c));
    }

    // A leading '(' is either a formula parenthesis "(phi && psi)" or the start
    // of an arithmetic expression "(a + 1) * 2 >= 0".  Try the formula reading
    // first; if it fails with a syntax error, rewind and try an atom.  Report
    // the error that got further into the input (final errors immediately).
    int parse_paren_or_atom() {
        const int save_pos = pos_, save_nodes = static_cast<int>(f_.nodes.size());
        const int save_var = next_var_;
        const size_t save_scope = scope_.size();
        auto restore = [&] {
            pos_ = save_pos; f_.nodes.resize(save_nodes);
            next_var_ = save_var; scope_.resize(save_scope);
        };
        try {
            return parse_paren_formula();
        } catch (const ParseError& e1) {
            if (e1.final) throw;
            restore();
            try {
                return parse_atom();
            } catch (const ParseError& e2) {
                if (e2.final || e2.pos > e1.pos) throw;
                throw e1;
            }
        }
    }

    int parse_paren_formula() {
        expect(TK::LPAREN, "'('");
        int left = parse_subformula();
        int node;
        if (at(TK::AND) || at(TK::OR)) {
            bool is_and = at(TK::AND);
            take();
            int right = parse_subformula();
            node = f_.add_node(is_and ? FormulaNode::make_and(left, right)
                                    : FormulaNode::make_or(left, right));
        } else if (at(TK::UNTIL) || at(TK::UNTIL_NS)) {
            bool strict = at(TK::UNTIL);
            take();
            double lo = 0.0, hi = std::numeric_limits<double>::infinity();
            if (at(TK::LBRACKET)) std::tie(lo, hi) = parse_window();
            int right = parse_subformula();
            node = f_.add_node(strict ? FormulaNode::make_until(lo, hi, left, right)
                                     : FormulaNode::make_until_ns(lo, hi, left, right));
        } else {
            node = left;   // plain parentheses
        }
        expect(TK::RPAREN, "')'");
        return node;
    }

    // '[' number ',' (number | inf) ']'
    std::pair<double, double> parse_window() {
        expect(TK::LBRACKET, "'['");
        if (!at(TK::NUMBER)) throw ParseError("Expected a number for the window start", cur().pos);
        double lo = take().num_val;
        expect(TK::COMMA, "','");
        double hi;
        if (at(TK::INF)) { take(); hi = std::numeric_limits<double>::infinity(); }
        else if (at(TK::NUMBER)) hi = take().num_val;
        else throw ParseError("Expected a number or 'inf' for the window end", cur().pos);
        if (hi < lo)
            throw ParseError("Invalid window: the end is smaller than the start", cur().pos,
                             /*final=*/true);
        expect(TK::RBRACKET, "']'");
        return {lo, hi};
    }

    // atom := expr op expr.  "lhs op rhs" is turned into "sig(s) op frz(x*)"
    // via lhs - rhs = S + F  <=>  S op -F.
    int parse_atom() {
        const int start = cur().pos;
        ExprP lhs = parse_expr();
        CompOp op = parse_comp_op();
        ExprP rhs = parse_expr();
        ExprP diff = mk(ExprKind::Sub, lhs, rhs);

        if (!has_frozen(diff)) {
            auto h = compile(diff, ExprKind::Signal);
            return f_.add_node(FormulaNode::make_signal_pred(std::move(h), op));
        }
        Parts parts = separate(diff, start);
        std::set<int> vars;
        collect_vars(diff, vars);
        auto h = compile(parts.sig, ExprKind::Signal);
        auto g = compile(neg(parts.frz), ExprKind::FreezeVar);
        return f_.add_node(FormulaNode::make_signal_constraint(
            std::move(h), std::move(g), op, std::vector<int>(vars.begin(), vars.end())));
    }

    CompOp parse_comp_op() {
        switch (cur().type) {
            case TK::GEQ: take(); return CompOp::GE;
            case TK::LEQ: take(); return CompOp::LE;
            case TK::GT:  take(); return CompOp::GT;
            case TK::LT:  take(); return CompOp::LT;
            default: throw ParseError("Expected a comparison operator (>=, >, <=, <)", cur().pos);
        }
    }

    // expr := term (('+'|'-') term)*
    ExprP parse_expr() {
        ExprP e = parse_term();
        while (at(TK::PLUS) || at(TK::MINUS)) {
            ExprKind k = at(TK::PLUS) ? ExprKind::Add : ExprKind::Sub;
            take();
            e = mk(k, e, parse_term());
        }
        return e;
    }

    // term := factor (('*'|'/') factor)*
    ExprP parse_term() {
        ExprP e = parse_factor();
        while (at(TK::STAR) || at(TK::SLASH)) {
            ExprKind k = at(TK::STAR) ? ExprKind::Mul : ExprKind::Div;
            take();
            e = mk(k, e, parse_factor());
        }
        return e;
    }

    // factor := number | signal | freeze_var | '-' factor | '(' expr ')'
    ExprP parse_factor() {
        if (at(TK::NUMBER)) return mk_num(take().num_val);
        if (at(TK::SIGNAL)) {
            Token t = take();
            return mk_leaf(ExprKind::Signal, signal_index(t.letter, t.pos));
        }
        if (at(TK::FREEZE_VAR)) {
            Token t = take();
            // "a*2", "a*b", "a*(..)": the user most likely meant a product.
            const TK nx = cur().type;
            if (cur().pos == t.end && (nx == TK::NUMBER || nx == TK::SIGNAL ||
                                       nx == TK::FREEZE_VAR || nx == TK::LPAREN))
                throw ParseError(std::string("'") + t.letter + "*' directly followed by a factor: '" +
                                 t.letter + "*' is a freeze variable; to multiply, write '" +
                                 t.letter + " * ...' with a space", t.pos, /*final=*/true);
            for (int i = static_cast<int>(scope_.size()) - 1; i >= 0; --i) {
                const Binding& b = scope_[i];
                if (b.letter == t.letter && b.has_subscript == t.has_subscript &&
                    (!t.has_subscript || b.subscript == t.subscript))
                    return mk_leaf(ExprKind::FreezeVar, b.var);
            }
            std::string name = std::string(1, t.letter) + "*";
            if (t.has_subscript) name += "_" + std::to_string(t.subscript);
            throw ParseError("Freeze variable '" + name + "' is not bound by an enclosing '" +
                             name + ".'", t.pos, /*final=*/true);
        }
        if (at(TK::MINUS)) { take(); return mk(ExprKind::Neg, parse_factor()); }
        if (at(TK::LPAREN)) {
            take();
            ExprP e = parse_expr();
            expect(TK::RPAREN, "')'");
            return e;
        }
        throw ParseError("Expected a number, signal, freeze variable or '('", cur().pos);
    }
};

}  // namespace

Formula parse_formula(const std::string& input) {
    std::vector<Token> tokens = Lexer(input).tokenize();
    return Parser(tokens).parse();
}

}  // namespace stl_star
