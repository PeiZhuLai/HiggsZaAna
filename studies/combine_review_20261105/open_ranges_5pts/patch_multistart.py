#!/usr/bin/env python3
"""Add a multi-start data fit to fTest for the masses with opened step/smearing ranges (mA10, 15, 21-29).
With the opened turn-on window (95-125 GeV) a single start can stop in a local minimum: at mA10 the single-start
fits sat at turn-on 95 / width 50 with NLL 2.0-3.2 above the best of 16 starts, at mA15 Pow1 +0.55.
The starts are the ones of range_full_refit/refit_limit.py: the first-fit result, then turn-on parameters
(name contains 'turnon', or Bern_step) at 97, 101, 105, 110, 115 GeV times width parameters (name contains
'width', or Bern_stepWidth) at 0.5, 5, 20 GeV, other parameters at the first-fit values. Kept: the lowest NLL
among converged fits (a failed fit only if none converged). Only data fits (MaxTries >= 3) are affected; the
toy fits (MaxTries 1, 2) keep a single start. The grid starts are set after sanitizePdfParams, which would
otherwise push width 0.5 up to 15% of the range.
Input : <revert.cpp>   Output: <patched.cpp>
Usage : patch_multistart.py <in.cpp> <out.cpp>"""
import sys

src = open(sys.argv[1]).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, 'anchor not found exactly once: %r' % old[:60]
    src = src.replace(old, new)

sub('bool PLOT_ONLY = false;\n', '''bool PLOT_ONLY = false;
// Multi-start data fits for the masses whose step/smearing ranges are opened in PdfModelBuilder.cc
// (openStepGausRanges), without the merged massALP == 0. Set in main from --mass_ALP.
bool MULTI_START = false;
static bool multiStartMass(int m) { return m == 10 || m == 15 || (m >= 21 && m <= 29); }
''')

sub('void runFit(RooAbsPdf *pdf, RooAbsData *data, double *NLL, int *stat_t, int MaxTries){\n', r'''// Refit from the turn-on x width grid and keep the lowest-NLL converged solution. The selection uses an
// un-offset NLL (each offset NLL has its own offset); *nllOut keeps the runFit convention of the chosen fit.
static void multiStartRefit(RooAbsPdf* pdf, RooAbsData* data, int* status, double* nllOut, bool fastMode) {
  std::unique_ptr<RooArgSet> params(pdf->getParameters(*data));
  if (!params) return;
  std::vector<RooRealVar*> pars, step, soft;
  for (RooAbsArg* a : *params) {
    RooRealVar* v = dynamic_cast<RooRealVar*>(a);
    if (!v || v->isConstant()) continue;
    pars.push_back(v);
    const std::string n = v->GetName();
    if (n.find("turnon") != std::string::npos || n == "Bern_step") step.push_back(v);
    else if (n.find("width") != std::string::npos || n == "Bern_stepWidth") soft.push_back(v);
  }
  if (step.empty() && soft.empty()) return;
  auto plainNll = [&]() { std::unique_ptr<RooAbsReal> n(pdf->createNLL(*data)); return n->getVal(); };
  auto snap = [&]() { std::vector<double> s; for (auto* v : pars) s.push_back(v->getVal()); return s; };
  auto load = [&](const std::vector<double>& s) { for (size_t i = 0; i < pars.size(); ++i) pars[i]->setVal(s[i]); };
  auto clampSet = [](RooRealVar* v, double x) { v->setVal(std::min(std::max(x, v->getMin()), v->getMax())); };
  struct Sol { bool failed; double plain; double nll; int st; std::vector<double> vals; int start; };
  Sol best{*status != 0, plainNll(), *nllOut, *status, snap(), 0};
  const Sol first = best;
  const std::vector<double> vp = best.vals;
  int k = 0;
  for (double t : {97.0, 101.0, 105.0, 110.0, 115.0}) {
    for (double w : {0.5, 5.0, 20.0}) {
      ++k;
      load(vp);
      for (auto* v : step) clampSet(v, t);
      for (auto* v : soft) clampSet(v, w);
      std::unique_ptr<RooAbsReal> nll(pdf->createNLL(*data, RooFit::Offset(true), RooFit::Optimize(2)));
      RooMinimizer minim(*nll);
      minim.setPrintLevel(-1);
      minim.setStrategy(fastMode ? 0 : 1);
      minim.setOffsetting(true);
      minim.optimizeConst(1);
      minim.setEps(fastMode ? 1e-4 : 1e-6);
      minim.setMaxFunctionCalls(fastMode ? 2000 : 20000);
      minim.setMaxIterations(fastMode ? 2000 : 20000);
      int st = minim.minimize("Minuit2", "migrad");
      if (!fastMode && st != 0) { minim.setStrategy(2); st = minim.minimize("Minuit2", "migrad"); }
      const double nv = nll->getVal();
      Sol s{st != 0, plainNll(), nv, st, snap(), k};
      if (!std::isfinite(s.plain)) continue;
      if ((!s.failed && best.failed) || (s.failed == best.failed && s.plain < best.plain)) best = s;
    }
  }
  load(best.vals);
  *status = best.st;
  *nllOut = best.nll;
  std::cout << "[INFO] multi-start " << pdf->GetName() << ": NLL first " << first.plain << " (status " << first.st
            << "), best " << best.plain << " (status " << best.st << ", start " << best.start << "/15)";
  for (auto* v : step) std::cout << " " << v->GetName() << "=" << v->getVal();
  for (auto* v : soft) std::cout << " " << v->GetName() << "=" << v->getVal();
  std::cout << std::endl;
}

void runFit(RooAbsPdf *pdf, RooAbsData *data, double *NLL, int *stat_t, int MaxTries){
''')

sub('  if (stat_t) *stat_t = status;\n  if (NLL) *NLL = bestNll;\n}\n',
    '  if (MULTI_START && MaxTries >= 3) multiStartRefit(pdf, data, &status, &bestNll, fastMode);\n\n'
    '  if (stat_t) *stat_t = status;\n  if (NLL) *NLL = bestNll;\n}\n')

sub('    std::cout << "[INFO] Processing for mA = " << mass_ALP << std::endl;\n',
    '    std::cout << "[INFO] Processing for mA = " << mass_ALP << std::endl;\n'
    '    MULTI_START = multiStartMass(mass_ALP);\n'
    '    std::cout << "[INFO] multi-start data fits: " << (MULTI_START ? "on" : "off") << std::endl;\n')

open(sys.argv[2], 'w').write(src)
