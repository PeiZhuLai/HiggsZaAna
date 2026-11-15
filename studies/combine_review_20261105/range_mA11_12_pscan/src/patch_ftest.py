#!/usr/bin/env python3
"""Adds the EXPSCAN_MASSES study switch to the fTest source of the production binary. At the listed mass points the
multi-start refit of the data fit also starts the exponent of the one-term shapes (Pow_p1_pow1, Lau_p1_lau1,
Exp_p1_exp1) from 15%, 50% and 85% of its range, on top of the unchanged 15 turn-on x width starts (60 starts).
The lowest-NLL converged solution is kept, so the result can only equal or improve on the nominal one.
Reason: at mA11 and mA12 the stored Pow1 is a local minimum; the Lau1 solution lies inside the Pow1 ranges with a
lower NLL (mA11 478.975 vs 479.81, mA12 525.625 vs 527.117), the multi-start only scans turn-on x width.
Unset by default, so the nominal configuration is unchanged. Every edit must match its anchor exactly once.
Usage : patch_ftest.py <input fTest source> <output fTest source>"""
import sys
s = open(sys.argv[1]).read()


def rep(old, new):
    global s
    assert s.count(old) == 1, ('anchor not found exactly once', old[:80], s.count(old))
    s = s.replace(old, new)


rep('''static void multiStartRefit(RooAbsPdf* pdf, RooAbsData* data, int* status, double* nllOut, bool fastMode) {''',
    '''// Study switch (Combine review): EXPSCAN_MASSES="11,12" (comma separated) also scans the exponent of the one-term
// shapes in the multi-start refit at those mass points. Unset by default, so the nominal configuration is untouched.
static bool isExpScanMass(int m)
{
  const char* env = getenv("EXPSCAN_MASSES");
  if (!env || !*env) return false;
  std::stringstream ss{std::string(env)};
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    try { if (std::stoi(tok) == m) return true; } catch (...) {}
  }
  return false;
}
static bool EXP_SCAN = false;

static void multiStartRefit(RooAbsPdf* pdf, RooAbsData* data, int* status, double* nllOut, bool fastMode) {''')
rep('''  std::vector<RooRealVar*> pars, step, soft;
  for (RooAbsArg* a : *params) {
    RooRealVar* v = dynamic_cast<RooRealVar*>(a);
    if (!v || v->isConstant()) continue;
    pars.push_back(v);
    const std::string n = v->GetName();
    if (n.find("turnon") != std::string::npos || n == "Bern_step") step.push_back(v);
    else if (n.find("width") != std::string::npos || n == "Bern_stepWidth") soft.push_back(v);
  }''', '''  std::vector<RooRealVar*> pars, step, soft, expo;
  for (RooAbsArg* a : *params) {
    RooRealVar* v = dynamic_cast<RooRealVar*>(a);
    if (!v || v->isConstant()) continue;
    pars.push_back(v);
    const std::string n = v->GetName();
    if (n.find("turnon") != std::string::npos || n == "Bern_step") step.push_back(v);
    else if (n.find("width") != std::string::npos || n == "Bern_stepWidth") soft.push_back(v);
    else if (EXP_SCAN && (n == "Pow_p1_pow1" || n == "Lau_p1_lau1" || n == "Exp_p1_exp1")) expo.push_back(v);
  }
  // f < 0: keep the exponent at its first-fit value (the nominal 15 starts); otherwise start it at min + f (max - min).
  std::vector<double> expoFrac = {-1.0};
  if (!expo.empty()) { expoFrac.push_back(0.15); expoFrac.push_back(0.50); expoFrac.push_back(0.85); }
  const int nStarts = 15 * (int)expoFrac.size();''')
rep('''  int k = 0;
  for (double t : {97.0, 101.0, 105.0, 110.0, 115.0}) {
    for (double w : {0.5, 5.0, 20.0}) {
      ++k;
      load(vp);
      for (auto* v : step) clampSet(v, t);
      for (auto* v : soft) clampSet(v, w);''', '''  int k = 0;
  for (double f : expoFrac) {
  for (double t : {97.0, 101.0, 105.0, 110.0, 115.0}) {
    for (double w : {0.5, 5.0, 20.0}) {
      ++k;
      load(vp);
      for (auto* v : step) clampSet(v, t);
      for (auto* v : soft) clampSet(v, w);
      if (f >= 0.) for (auto* v : expo) clampSet(v, v->getMin() + f * (v->getMax() - v->getMin()));''')
rep('''      if ((!s.failed && best.failed) || (s.failed == best.failed && s.plain < best.plain)) best = s;
    }
  }''', '''      if ((!s.failed && best.failed) || (s.failed == best.failed && s.plain < best.plain)) best = s;
    }
  }
  }''')
rep('''<< "), best " << best.plain << " (status " << best.st << ", start " << best.start << "/15)";''',
    '''<< "), best " << best.plain << " (status " << best.st << ", start " << best.start << "/" << nStarts << ")";''')
rep('''  for (auto* v : soft) std::cout << " " << v->GetName() << "=" << v->getVal();
  std::cout << std::endl;''', '''  for (auto* v : soft) std::cout << " " << v->GetName() << "=" << v->getVal();
  for (auto* v : expo) std::cout << " " << v->GetName() << "=" << v->getVal();
  std::cout << std::endl;''')
rep('''    MULTI_START = multiStartMass(mass_ALP);''', '''    MULTI_START = multiStartMass(mass_ALP);
    EXP_SCAN = isExpScanMass(mass_ALP);
    std::cout << "[INFO] exponent scan in the multi-start: " << (EXP_SCAN ? "on" : "off") << std::endl;''')
open(sys.argv[2], 'w').write(s)
print('patched')
