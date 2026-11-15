#!/usr/bin/env python3
"""Adds the ORDER2_MASSES study switch to the fTest source of the production binary. At the listed mass points the
next order of the Exponential, PowerLaw and Laurent families (the order at which the envelope F-test stops) is also put
into the envelope when its sideband GOF passes (> minEnvelopeGof); the members chosen by the nominal rule are kept and
the Bernstein member is not touched. Unset by default, so the nominal configuration is unchanged.
Every edit must match its anchor exactly once.
Usage : patch_ftest.py <input fTest source> <output fTest source>"""
import sys
s = open(sys.argv[1]).read()


def rep(old, new):
    global s
    assert s.count(old) == 1, ('anchor not found exactly once', old[:80], s.count(old))
    s = s.replace(old, new)


rep('''static int envMaxOrder(const std::string& funcType, int mass_ALP)
{''', '''// Study switch (Combine review): ORDER2_MASSES="11,12" (comma separated) also puts the next order of the
// Exponential, PowerLaw and Laurent families into the envelope at those mass points, when its sideband GOF passes,
// although the F-test does not ask for it. Unset by default, so the nominal configuration is untouched.
static bool isOrder2Mass(int m)
{
  const char* env = getenv("ORDER2_MASSES");
  if (!env || !*env) return false;
  std::stringstream ss{std::string(env)};
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    try { if (std::stoi(tok) == m) return true; } catch (...) {}
  }
  return false;
}

static int envMaxOrder(const std::string& funcType, int mass_ALP)
{''')
rep('''				std::vector<EnvelopeCandidate> acceptedEnvelopeCandidates;
''', '''				std::vector<EnvelopeCandidate> acceptedEnvelopeCandidates;
				std::vector<EnvelopeCandidate> order2Extras;
''')
rep('''									if ((prob < upperEnvThreshold) ) {
''', '''									const bool order2Extra = !(prob < upperEnvThreshold) && order > 1
									                         && *funcType != "Bernstein" && isOrder2Mass(mass_ALP);
									if (order2Extra && fitStatus != 5 && gofProb > minEnvelopeGof) {
										EnvelopeCandidate extra;
										extra.pdf = bkgPdf;
										extra.family = *funcType;
										extra.name = bkgPdf->GetName();
										extra.order = order;
										extra.gof = gofProb;
										extra.score = myNll + bkgPdf->getVariables()->getSize();
										extra.isTruth = (order == truthOrder);
										order2Extras.push_back(extra);
										std::cout << "[INFO] ORDER2_MASSES: next order " << extra.name << " gof " << gofProb
										          << " 2xNLL + c is " << extra.score << std::endl;
									} else if (order2Extra) {
										std::cout << "[INFO] ORDER2_MASSES: next order " << bkgPdf->GetName()
										          << " not added, fit status " << fitStatus << " gof " << gofProb << std::endl;
									}
									if ((prob < upperEnvThreshold) ) {
''')
rep('''		      if (finalEnvelopeCandidates.empty()) {
''', '''		      for (const auto& cand : order2Extras) {
		        std::cout << "[INFO] ORDER2_MASSES: adding " << cand.name << " score " << cand.score << " gof " << cand.gof << std::endl;
		        addUniqueEnvelopeCandidate(cand);
		      }

		      if (finalEnvelopeCandidates.empty()) {
''')
open(sys.argv[2], 'w').write(s)
print('patched')
