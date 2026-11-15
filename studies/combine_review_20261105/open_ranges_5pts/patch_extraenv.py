#!/usr/bin/env python3
"""Study option for fTest: EXTRA_ENV_PDFS=<name>[,<name>...] adds the named functions to the envelope on top of the
one-per-family envelope, if they were accepted by the envelope scan (F-test threshold and GOF > minEnvelopeGof).
Unset, the envelope is unchanged. The accepted candidates are always printed (name, score 2NLL+npar, GOF) so the
choice can be made from the log. A second Bernstein order is refused: all Bernstein orders share Bern_step,
Bern_stepWidth and Bern_gsigma, so two of them in one envelope cannot keep their own fitted values.
Applied after patch_multistart.py.
Input : <in.cpp>   Output: <out.cpp>
Usage : patch_extraenv.py <in.cpp> <out.cpp>"""
import sys

src = open(sys.argv[1]).read()
old = '\t\t      if (finalEnvelopeCandidates.empty()) {\n'
assert src.count(old) == 1, 'anchor not found exactly once'
new = r'''		      for (const auto& cand : acceptedEnvelopeCandidates) {
		        std::cout << "[INFO] accepted candidate " << cand.name << " family " << cand.family
		                  << " score " << cand.score << " gof " << cand.gof << std::endl;
		      }
		      if (const char* extraEnv = getenv("EXTRA_ENV_PDFS")) {
		        std::stringstream ss(extraEnv);
		        std::string want;
		        while (std::getline(ss, want, ',')) {
		          if (want.empty()) continue;
		          const EnvelopeCandidate* hit = nullptr;
		          for (const auto& cand : acceptedEnvelopeCandidates) if (cand.name == want) hit = &cand;
		          if (!hit) { std::cerr << "[WARN] EXTRA_ENV_PDFS: " << want << " is not an accepted candidate, not added" << std::endl; continue; }
		          bool bernTaken = false;
		          for (const auto& c : finalEnvelopeCandidates) if (c.family == "Bernstein" && hit->family == "Bernstein") bernTaken = true;
		          if (bernTaken) { std::cerr << "[WARN] EXTRA_ENV_PDFS: " << want << " refused, a Bernstein order is already in the envelope (shared parameters)" << std::endl; continue; }
		          std::cout << "[INFO] EXTRA_ENV_PDFS: adding " << want << " score " << hit->score << " gof " << hit->gof << std::endl;
		          addUniqueEnvelopeCandidate(*hit);
		        }
		      }

		      if (finalEnvelopeCandidates.empty()) {
'''
src = src.replace(old, new)
if '#include <sstream>' not in src:
    src = src.replace('#include <cstdlib>\n', '#include <cstdlib>\n#include <sstream>\n', 1)
open(sys.argv[2], 'w').write(src)
