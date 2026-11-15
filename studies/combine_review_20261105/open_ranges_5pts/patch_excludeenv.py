#!/usr/bin/env python3
"""Study option for fTest: EXCLUDE_ENV_PDFS=<name>[,<name>...] removes the named functions from the final envelope,
after the one-per-family selection, the fill up to the target size and EXTRA_ENV_PDFS. The envelope is not refilled
and no family-fallback member takes the place of a removed one, so the envelope is the unmodified one minus the
named functions, with every remaining member at the same fitted values. Unset, the envelope is unchanged.
Used at mA10, where the bias screen with truth Pow1 fails the 0.2 criterion.
Applied after patch_extraenv.py.
Input : <in.cpp>   Output: <out.cpp>
Usage : patch_excludeenv.py <in.cpp> <out.cpp>"""
import sys

src = open(sys.argv[1]).read()
old = '\t\t      if (finalEnvelopeCandidates.empty()) {\n'
assert src.count(old) == 1, 'anchor not found exactly once'
assert 'EXTRA_ENV_PDFS' in src, 'apply patch_extraenv.py first'
new = r'''		      if (const char* exclEnv = getenv("EXCLUDE_ENV_PDFS")) {
		        std::stringstream ss(exclEnv);
		        std::string drop;
		        while (std::getline(ss, drop, ',')) {
		          if (drop.empty()) continue;
		          const size_t before = finalEnvelopeCandidates.size();
		          finalEnvelopeCandidates.erase(std::remove_if(finalEnvelopeCandidates.begin(), finalEnvelopeCandidates.end(),
		                                        [&](const EnvelopeCandidate& c) { return c.name == drop; }),
		                                        finalEnvelopeCandidates.end());
		          if (finalEnvelopeCandidates.size() == before)
		            std::cerr << "[WARN] EXCLUDE_ENV_PDFS: " << drop << " is not in the envelope, nothing removed" << std::endl;
		          else
		            std::cout << "[INFO] EXCLUDE_ENV_PDFS: removing " << drop << " from the envelope" << std::endl;
		        }
		      }

		      if (finalEnvelopeCandidates.empty()) {
'''
src = src.replace(old, new)
if '#include <algorithm>' not in src:
    src = src.replace('#include <cstdlib>\n', '#include <cstdlib>\n#include <algorithm>\n', 1)
open(sys.argv[2], 'w').write(src)
