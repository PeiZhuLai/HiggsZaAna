import os, re, sys
new, v1, out = sys.argv[1:4]; bad = 0; rows = []
E = {'2022preEE': '2022', '2022postEE': '2022EE', '2023preBPix': '2023', '2023postBPix': '2023BPix', '2024': '2024'}
for m in range(1, 31):
    for s in ('', '_xsec'):
        n = '%d_pruned_datacard_leptons%s.txt' % (m, s)
        a = open(os.path.join(v1, n)).read(); b = open(os.path.join(new, n)).read()
        k = len(re.findall(r'\binterpolate_yield_', a))
        a2 = re.sub(r'\binterpolate_yield_(2022preEE|2022postEE|2023preBPix|2023postBPix|2024)_(ele|mu)\b',
                    lambda x: 'CMS_NPS25014_interp_yield_%s_%s' % (E[x.group(1)], x.group(2)), a)
        same = [l.split() for l in a2.splitlines()] == [l.split() for l in b.splitlines()]
        bad += 0 if same and 'interpolate_yield' not in b else 1
        rows.append('%s %-34s old-name occurrences renamed: %d' % ('[OK]   ' if same else '[DIFF] ', n, k))
open(out, 'w').write('\n'.join(rows) + '\n%s\n' % ('ALL OK' if not bad else 'PROBLEMS %d' % bad)); sys.exit(1 if bad else 0)
