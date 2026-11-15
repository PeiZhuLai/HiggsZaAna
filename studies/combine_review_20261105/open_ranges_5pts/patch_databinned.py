#!/usr/bin/env python3
"""fTest fix: the binned data saved to the multipdf workspace (roohist_data_mass_<cat>, the data_obs of the datacards)
was built from the full workspace dataset (95-180 GeV) instead of the dataset reduced to the fit window. With a
window narrower than the dataset (--mhLow 100 on a 95-180 GeV dataset) the events below the window were piled into
the first bin (mA10 test: first bin 27 instead of 8, total 120 instead of 101). The fits always used the reduced
dataset. When the dataset already lies inside the window (every production point, mA1 included) data == dataFull
and the output is unchanged.
Input : <in.cpp>   Output: <out.cpp>
Usage : patch_databinned.py <in.cpp> <out.cpp>"""
import sys
src = open(sys.argv[1]).read()
old = '"data",*mass,*dataFull);'
assert src.count(old) == 2, 'expected 2 dataBinned lines, found %d' % src.count(old)
open(sys.argv[2], 'w').write(src.replace(old, '"data",*mass,*data);'))
