#!/bin/bash
# Run combine AsymptoticLimits on each merged low-mA datacard, collect results.
# MUST run inside the flashgg combine cmsenv (text2workspace.py / combine on PATH).
set -euo pipefail

OUTDIR="${1:-$(cd "$(dirname "$0")" && pwd)/output}"
CARD_DIR="${OUTDIR}/datacards"
LIM_DIR="${OUTDIR}/limits"
mkdir -p "${LIM_DIR}"

if ! command -v combine >/dev/null 2>&1; then
    echo "[limits] ERROR: 'combine' not on PATH. Did you cmsenv the flashgg area?" >&2
    exit 1
fi

shopt -s nullglob
cards=("${CARD_DIR}"/datacard_*.txt)
if [[ ${#cards[@]} -eq 0 ]]; then
    echo "[limits] no datacards in ${CARD_DIR}" >&2
    exit 1
fi

echo "tag,exp_m2,exp_m1,exp_med,exp_p1,exp_p2" > "${LIM_DIR}/limits.csv"

for card in "${cards[@]}"; do
    tag="$(basename "${card}" .txt | sed 's/^datacard_//')"
    echo "==== combine ${tag} ===="
    ws="${LIM_DIR}/ws_${tag}.root"
    text2workspace.py "${card}" -o "${ws}"
    # parametric bkg -> need RooMultiPdf; AsymptoticLimits, blinded (expected only)
    # merged limits are O(100) in r (huge bkg) -> default rMax=20 fails to
    # bracket; widen the scan ([[ref_hza_limit_rmax_mh_bug]]).
    combine -M AsymptoticLimits "${ws}" \
        --run blind --noFitAsimov \
        --rMin 0 --rMax "${RMAX:-1000}" --rRelAcc 0.01 \
        -n "_merged_${tag}" --mass 125 \
        > "${LIM_DIR}/log_${tag}.txt" 2>&1 || {
            echo "[limits] combine failed for ${tag}; see ${LIM_DIR}/log_${tag}.txt" >&2
            continue
        }
    res="higgsCombine_merged_${tag}.AsymptoticLimits.mH125.root"
    [[ -f "${res}" ]] && mv -f "${res}" "${LIM_DIR}/"
    python3 - "${LIM_DIR}/${res}" "${tag}" "${LIM_DIR}/limits.csv" <<'PY'
import sys, ROOT
res, tag, csv = sys.argv[1], sys.argv[2], sys.argv[3]
f = ROOT.TFile.Open(res)
t = f.Get("limit") if f and not f.IsZombie() else None
# map combine quantileExpected -> our 5 bands; tolerate missing quantiles
want = {0.025: "m2", 0.16: "m1", 0.5: "med", 0.84: "p1", 0.975: "p2"}
got = {}
if t:
    for e in t:
        q = round(e.quantileExpected, 3)
        for qe, key in want.items():
            if abs(q - qe) < 0.01:
                got[key] = e.limit
if "med" in got:
    order = ["m2", "m1", "med", "p1", "p2"]
    row = [f"{got.get(k, float('nan')):.4g}" for k in order]
    with open(csv, "a") as fh:
        fh.write(f"{tag}," + ",".join(row) + "\n")
    print(f"[limits] {tag}: median exp = {got['med']:.4g}  "
          f"(bands: {sum(k in got for k in order)}/5)")
else:
    print(f"[limits] {tag}: combine did not produce a median (channel too weak / "
          f"no crossing in [0,rMax]); skipped")
PY
done

echo "[limits] summary -> ${LIM_DIR}/limits.csv"
cat "${LIM_DIR}/limits.csv"
