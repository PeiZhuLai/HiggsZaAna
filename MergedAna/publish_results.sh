#!/bin/bash
# Publish the plots and text tables of one HZa_merged chain run to EOS home:
#   ${EOS_RESULTS}/output_plots/<tag>/   limits, bkg model, signal model, bias plots
#   ${EOS_RESULTS}/output_tables/<tag>/  expected limits, envelope/fTest tables,
#                                        datacards, bias summary, MANIFEST.txt
# Logs and the ROOT files the next stage reads (multipdf, t2w workspaces, combine
# roots) stay where the chain wrote them; logs belong on AFS.
#
# Called at the end of 1_grand_merged.sh (Stage 8), or by hand after the fact:
#   PLOTDIR=... BKG_WSDIR=... LIMDIR=... bash publish_results.sh <tag>
#   BIAS_DIR=$FF/Combine/Checks/Bias_merged_phosel/envtest2_s2 bash publish_results.sh envtest2_s2
# Each source is optional: an unset or missing source is reported and skipped.
#
# EOS_RESULTS uses the real namespace /eos/user/p/pelai: /eos/home-p is only a FUSE
# alias, and xrootd paths built from it fail with "No such file" on condor.
# Every copied file is checked by size on EOS; any mismatch, or nothing copied at
# all, exits non-zero so a silent failed write on a batch node cannot pass.
set -u -o pipefail

TAG="${1:?usage: publish_results.sh <tag>}"
EOS_RESULTS="${EOS_RESULTS:-/eos/user/p/pelai/HZa_merged}"
MASSES="${MASSES:-0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9}"
LEPS="${LEPS:-ele mu}"
YEAR="${YEAR:-2024}"
PLOTDIR="${PLOTDIR:-}"
BKG_WSDIR="${BKG_WSDIR:-}"
LIMDIR="${LIMDIR:-}"
DC_OUT="${DC_OUT:-}"
SIG_DIR="${SIG_DIR:-}"
BIAS_DIR="${BIAS_DIR:-}"

P="${EOS_RESULTS}/output_plots/${TAG}"
T="${EOS_RESULTS}/output_tables/${TAG}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/publish_results.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
MANIFEST_TMP="$WORK/manifest.lines"; : > "$MANIFEST_TMP"
n_ok=0; n_bad=0

# copy <src> <dest dir>: FUSE cp first, xrdcp as fallback, then verify the size.
put() {
    local src="$1" dst_dir="$2" dst
    dst="${dst_dir}/$(basename "$src")"
    mkdir -p "$dst_dir" 2>/dev/null
    cp -f "$src" "$dst" 2>/dev/null || \
        xrdcp -f -p -s "$src" "root://eosuser.cern.ch/${dst}" 2>/dev/null
    local s_src s_dst
    s_src=$(stat -c %s "$src" 2>/dev/null || echo -1)
    # FUSE may be missing on a batch node: fall back to asking the MGM directly
    s_dst=$(stat -c %s "$dst" 2>/dev/null || \
            xrdfs root://eosuser.cern.ch stat "$dst" 2>/dev/null | awk '/^Size:/{print $2}')
    s_dst="${s_dst:--2}"
    if [[ "$s_src" == "$s_dst" ]]; then
        n_ok=$((n_ok+1)); echo "$s_src $src -> $dst" >> "$MANIFEST_TMP"
    else
        n_bad=$((n_bad+1)); echo "[ERROR] copy failed: $src -> $dst (src $s_src B, eos $s_dst B)" >&2
    fi
}

# put_glob <dest dir> <glob...>: copy every existing match.
put_glob() {
    local dst_dir="$1"; shift
    local f
    for f in "$@"; do [[ -f "$f" ]] && put "$f" "$dst_dir"; done
}

have() { [[ -n "$1" && -d "$1" ]] || { echo "  [SKIP] $2: ${1:-unset}"; return 1; }; }

echo "publish ${TAG} -> ${P} , ${T}"

# limits: plots + a plain-text table read from the combine logs
if have "$PLOTDIR" "limit plots (PLOTDIR)"; then
    put_glob "$P/limits" "$PLOTDIR"/*.pdf "$PLOTDIR"/*.png
fi
if have "$LIMDIR" "limit tables (LIMDIR)"; then
    tab="$WORK/limits_expected.txt"
    {
        echo "# median expected limit on r per m_a; source ${LIMDIR}/combine_<m>.log"
        printf "%-6s %10s %10s %10s %10s %10s\n" mA exp2p5 exp16 exp50 exp84 exp97p5
        for m in $MASSES; do
            lg="${LIMDIR}/combine_${m}.log"
            [[ -f "$lg" ]] || continue
            vals=$(for q in " 2.5%" "16.0%" "50.0%" "84.0%" "97.5%"; do
                       v=$(grep "Expected ${q}" "$lg" | sed 's/.*< //' | head -1); echo "${v:--}"
                   done | tr '\n' ' ')
            printf "%-6s %10s %10s %10s %10s %10s\n" "$m" $vals
        done
    } > "$tab"
    put "$tab" "$T"
    # run_variants.sh tars the variant's datacards after the chain ends, so inside the
    # chain this tgz is the previous run's; use it only for an after-the-fact publish
    [[ -z "$DC_OUT" ]] && put_glob "$T/limits" "$LIMDIR"/sigmodel_datacards.tgz
fi

# background model: banded Data+envelope plot, per-function fTest/gof plots and tables
if have "$BKG_WSDIR" "bkg model (BKG_WSDIR)"; then
    for m in $MASSES; do
        b="${BKG_WSDIR}/${m}"; [[ -d "$b" ]] || continue
        put_glob "$P/bkg_model/M${m}" "$b"/AllFitResults/*.pdf "$b"/BkgPlots/*.pdf \
            "$b"/HZAmassInde_fTest/*.pdf "$b"/HZAmassInde_fTest/*.png
        put_glob "$T/bkg_envelope/M${m}" "$b"/HZAmassInde_fTest/*.txt
    done
fi

# datacards (fixed path, overwritten by the next run: publish right after the run)
if have "$DC_OUT" "datacards (DC_OUT)"; then
    for m in $MASSES; do put_glob "$T/datacards" "${DC_OUT}/${m}_pruned_datacard_leptons.txt"; done
fi

# signal model plots (fixed path per lep, also overwritten by the next run)
if have "$SIG_DIR" "signal model (SIG_DIR)"; then
    for lep in $LEPS; do
        for m in $MASSES; do
            put_glob "$P/signal_model/${lep}" "${SIG_DIR}/outdir_${lep}/signalFit/Plots"/*_${m}_${YEAR}_${lep}*.pdf
        done
    done
fi

# bias study: pull plots, gauss-fit json and the summary
if have "$BIAS_DIR" "bias study (BIAS_DIR)"; then
    for m in $MASSES; do
        d="${BIAS_DIR}/mA_${m}/merged"; [[ -d "$d" ]] || continue
        put_glob "$P/bias/M${m}" "$d"/BiasPlots/*.pdf "$d"/BiasPlots/*.png
        put_glob "$T/bias/M${m}" "$d"/BiasJson/*.json
    done
    put_glob "$T/bias" "$BIAS_DIR"/bias_summary.txt
fi

if [[ "$n_ok" -gt 0 ]]; then
    { echo "# publish ${TAG} $(date '+%F %T') host=$(hostname) ok=${n_ok} failed=${n_bad}"
      echo "# <bytes> <source> -> <eos path>"
      cat "$MANIFEST_TMP"; } > "$WORK/MANIFEST.txt"
    put "$WORK/MANIFEST.txt" "$T"
fi
echo "publish ${TAG}: copied ${n_ok}, failed ${n_bad}"
[[ "$n_ok" -gt 0 && "$n_bad" -eq 0 ]] || exit 1
