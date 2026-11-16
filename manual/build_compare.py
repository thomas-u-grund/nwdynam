#!/usr/bin/env python3
"""Reads ex/r/goldfish_results.csv (written by ex/r/compare.R) and ex/r/nwdynam_results.csv (written by ex/ex91_compare_stata.do)
and writes gen_compare.tex (the table of chapter 11) and gen_compare_nums.tex (the maxima quoted in the text)."""
import csv, os
here = os.path.dirname(os.path.abspath(__file__))
def read(f):
    return {(r["model"], r["term"].lower()): r for r in csv.DictReader(open(os.path.join(here, "ex", "r", f)))}
g = read("goldfish_results.csv"); s = read("nwdynam_results.csv")
labels = {
 "soc_choice": ("Social Evolution", "choice: inertia, recip, indeg"),
 "soc_rate": ("Social Evolution", "rate: indeg, outdeg"),
 "soc_rate_int": ("Social Evolution", "rate with intercept: indeg, outdeg"),
 "em_choice": ("Emails", "choice: inertia, recip, indeg"),
 "em_choice_trans": ("Emails", "choice: inertia, recip, indeg, trans"),
 "em_choice_win": ("Emails", "choice: inertia (window 0.5), recip, indeg"),
 "em_rate": ("Emails", "rate: indeg, outdeg"),
 "em_rate_int": ("Emails", "rate with intercept: indeg, outdeg"),
}
def rel(a, b): return abs(a - b) / max(abs(b), 1e-300)
out = []
mx_b = mx_se = mx_ll = 0.0
where = None
done = set()
for (m, t), r in s.items():
    gr = g[(m, "intercept" if t == "intercept" else t)]
    b1, b2 = float(r["estimate"]), float(gr["estimate"])
    s1, s2 = float(r["se"]), float(gr["se"])
    if rel(b1, b2) > mx_b: mx_b, where = rel(b1, b2), (m, t)
    mx_se = max(mx_se, rel(s1, s2))
    ll1, ll2 = float(r["loglik"]), float(gr["loglik"])
    if m not in done:
        mx_ll = max(mx_ll, rel(ll1, ll2))
for m in labels:
    rows = [(t, r) for (mm, t), r in s.items() if mm == m]
    ll1 = float(rows[0][1]["loglik"]); ll2 = float(g[(m, rows[0][0])]["loglik"])
    out.append(r"\midrule")
    out.append(r"\multicolumn{6}{@{}p{14.6cm}@{}}{\textbf{%s, %s}; log likelihood %.4f (nwdynam), %.4f (goldfish)}\\" % (labels[m][0], labels[m][1], ll1, ll2))
    for t, r in rows:
        gr = g[(m, t)]
        b1, b2 = float(r["estimate"]), float(gr["estimate"]); s1, s2 = float(r["se"]), float(gr["se"])
        out.append(r"\quad \texttt{%s} & %.5f & %.5f & %.1e & %.5f & %.5f \\" % (t.replace("_", r"\_").capitalize() if t == "intercept" else t, b1, b2, rel(b1, b2), s1, s2))
tab = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
       r"Effect & nwdynam & goldfish & rel.\ diff. & s.e.\ nwdynam & s.e.\ goldfish \\"] + out + [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(here, "gen_compare.tex"), "w").write("\n".join(tab) + "\n")
nums = [r"\newcommand{\cmpMaxB}{%.1e}" % mx_b, r"\newcommand{\cmpMaxSE}{%.1e}" % mx_se, r"\newcommand{\cmpMaxLL}{%.1e}" % mx_ll,
        r"\newcommand{\cmpNModels}{%d}" % len(labels),
        r"\newcommand{\cmpMaxWhere}{%s, \texttt{%s}}" % (labels[where[0]][0] + ", " + labels[where[0]][1].split(":")[0], where[1])]
open(os.path.join(here, "gen_compare_nums.tex"), "w").write("\n".join(nums) + "\n")
print("max rel diff: estimates %.2e, s.e. %.2e, ll %.2e" % (mx_b, mx_se, mx_ll))
