"""TTree::Draw-equivalent histogram filling from a single pass over a chain.

chain.Draw("expr>>h", "(weight) * (cut)", "goff") evaluates the TTreeFormula strings per
entry, skips entries where (weight) * (cut) is 0 and adds the others in entry order through
TH1::FillN.  evaluate_draw_requests() evaluates every distinct expr / weight / cut of many
such requests with the same TTreeFormula objects in ONE loop over the chain (compiled C++),
and fill_like_draw() fills each histogram with TH1::FillN over the same entries in the same
order, so contents, errors, entries and statistics are identical to one Draw per histogram.

TTreeFormula (not a C++ translation) is used on purpose: it is not plain C++ arithmetic,
e.g. "(w) * (s) / (c)" is evaluated as w * (s / c), and Entry$ is the chain entry.
"""

from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from root_compat import import_pyroot

ROOT = import_pyroot()

_CPP = r"""
#include "TTree.h"
#include "TTreeFormula.h"
#include <string>
#include <vector>

namespace single_pass_draw {
// values: n_values x n (row major); bits: ceil(n_cuts / 64) x n, bit k%64 of word k/64 = cut k != 0.
// Returns n on success, the number of entries read if the loop stopped early, or
// -(index + 1) of the first formula that does not compile.
Long64_t Evaluate(TTree *tree, const std::vector<std::string> &valueExprs, const std::vector<std::string> &cutExprs,
                  double *values, long *signedBits, Long64_t n)
{
   ULong64_t *bits = reinterpret_cast<ULong64_t *>(signedBits);  // numpy int64 binds to long*, uint64 to nothing
   std::vector<TTreeFormula *> vf, cf;
   Long64_t status = n;
   if (n > 0) tree->LoadTree(0);  // formulas are compiled against the first tree, as in TTree::Draw
   auto build = [&](const std::vector<std::string> &exprs, std::vector<TTreeFormula *> &out, Long64_t offset) {
      for (size_t k = 0; k < exprs.size(); ++k) {
         auto *f = new TTreeFormula("single_pass_draw", exprs[k].c_str(), tree);
         out.push_back(f);
         if (f->GetNdim() == 0 && status == n) status = -(offset + Long64_t(k) + 1);
      }
   };
   build(valueExprs, vf, 0);
   build(cutExprs, cf, Long64_t(valueExprs.size()));
   if (status == n) {
      Int_t currentTree = -1;
      for (Long64_t i = 0; i < n; ++i) {
         if (tree->LoadTree(i) < 0) {
            status = i;
            break;
         }
         if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            for (auto *f : vf) f->UpdateFormulaLeaves();
            for (auto *f : cf) f->UpdateFormulaLeaves();
         }
         for (size_t k = 0; k < vf.size(); ++k) {
            vf[k]->GetNdata();
            values[k * n + i] = vf[k]->EvalInstance(0);
         }
         for (size_t k = 0; k < cf.size(); ++k) {
            cf[k]->GetNdata();
            if (cf[k]->EvalInstance(0) != 0) bits[(k / 64) * n + i] |= (ULong64_t(1) << (k % 64));
         }
      }
   }
   for (auto *f : vf) delete f;
   for (auto *f : cf) delete f;
   return status;
}
}
"""

_DECLARED = False


def _declare() -> None:
    global _DECLARED
    if not _DECLARED:
        ROOT.gInterpreter.Declare(_CPP)
        _DECLARED = True


def evaluate_draw_requests(chain, requests: Iterable[Tuple[str, str, str]]) -> Dict[Tuple[str, str], np.ndarray]:
    """{("x", expr) | ("w", weight) | ("c", cut): per-entry float64 array} for all (expr, weight, cut).

    Cuts come back as 1.0 / 0.0 (the value of a TTreeFormula boolean); in Draw the selection
    "(weight) * (cut)" is exactly weight * that value.
    """
    _declare()
    value_keys, cut_keys = [], []
    for expr, weight, cut in requests:
        for key in (("x", expr), ("w", weight)):
            if key not in value_keys:
                value_keys.append(key)
        if cut not in cut_keys:
            cut_keys.append(cut)
    n_entries = int(chain.GetEntries())
    values = np.zeros(len(value_keys) * n_entries, dtype=np.float64)
    bits = np.zeros(((len(cut_keys) + 63) // 64) * n_entries, dtype=np.int64)
    value_exprs = ROOT.std.vector("std::string")([key[1] for key in value_keys])
    cut_exprs = ROOT.std.vector("std::string")(cut_keys)
    status = ROOT.single_pass_draw.Evaluate(chain, value_exprs, cut_exprs, values, bits, n_entries)
    if status < 0:
        bad = -status - 1
        bad_expr = value_keys[bad][1] if bad < len(value_keys) else cut_keys[bad - len(value_keys)]
        raise RuntimeError(f"TTreeFormula cannot compile '{bad_expr}'")
    if status != n_entries:
        raise RuntimeError(f"read {status} of {n_entries} entries")
    values = values.reshape(len(value_keys), n_entries)
    bits = bits.view(np.uint64).reshape(-1, n_entries)
    out = {key: values[k] for k, key in enumerate(value_keys)}
    for k, cut in enumerate(cut_keys):
        out[("c", cut)] = ((bits[k // 64] >> np.uint64(k % 64)) & np.uint64(1)).astype(np.float64)
    return out


def fill_like_draw(hist, values, expr: str, weight: str, cut: str, rows: Optional[np.ndarray] = None) -> None:
    """Fill `hist` as chain.Draw(expr>>hist, "(weight) * (cut)") would; `rows` selects/orders the entries."""
    x = values[("x", expr)]
    fill_weight = values[("w", weight)] * values[("c", cut)]
    if rows is not None:
        x, fill_weight = x[rows], fill_weight[rows]
    selected = fill_weight != 0
    n_fill = int(np.count_nonzero(selected))
    if n_fill:
        hist.FillN(n_fill, np.ascontiguousarray(x[selected]), np.ascontiguousarray(fill_weight[selected]))
