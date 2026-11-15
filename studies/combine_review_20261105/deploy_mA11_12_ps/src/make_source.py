#!/usr/bin/env python3
"""Production fTest source of the mA11/mA12 deployment: the side-study source (EXPSCAN_MASSES switch) with the mass
list written into the code, so the production command line of shellScripts/bkg/fit_bkg.sh reproduces it without
any environment variable. Everything else is the side-study source unchanged.
Input : ../range_mA11_12_pscan/src/fTest_ALP_turnOn.pscan.cpp
Output: <output fTest source>
Usage : make_source.py <input fTest source> <output fTest source>"""
import sys
s = open(sys.argv[1]).read()
start = s.index('// Study switch (Combine review): EXPSCAN_MASSES')
end = s.index('static bool EXP_SCAN = false;')
assert s.count('// Study switch (Combine review): EXPSCAN_MASSES') == 1 and s.count('static bool EXP_SCAN = false;') == 1
old = s[start:end]
assert 'getenv("EXPSCAN_MASSES")' in old and old.count('static bool isExpScanMass(int m)') == 1, old
new = ('// Combine review: at mA11 and mA12 the multi-start refit of the data fit also starts the exponent of the one-term\n'
       '// shapes (Pow1, Lau1, Exp1) from 15/50/85% of its range. With the turn-on x width starts alone the stored Pow1 was\n'
       '// a local minimum (the Lau1 solution lies inside the Pow1 ranges with a lower NLL).\n'
       'static bool isExpScanMass(int m) { return m == 11 || m == 12; }\n')
s = s[:start] + new + s[end:]
assert 'EXPSCAN_MASSES' not in s
open(sys.argv[2], 'w').write(s)
print('written')
