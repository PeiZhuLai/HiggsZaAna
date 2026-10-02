"""MC overlap roles for MCOverlapTagger and for the post-production check
(scripts/check_mc_overlap_veto.py). Kept dependency-free so the check can import it without HiggsDNA."""
import logging
import os

logger = logging.getLogger(__name__)

# [PZ 2026-09-29] Fail-closed sample classification.
# Until 2026-09-27 the roles were a chain of `"X" in file` tests and anything that matched none of
# them was passed through with NO overlap removal and no message. The 2024 DY+jets samples are
# named DYto2E/DYto2Mu/DYto2Tau (not DYto2L), so their iso-photon veto was silently skipped and
# ~21-23% of their weight double-counted DY+gamma, found only months later in the BDT inputs.
# Now every MC file must match exactly one of the three lists below; an unmatched MC file raises.
# Patterns are substrings of the input URL (the tagger sees the original /store/... URL).
#   veto : "+jets" side, keep events with NO isolated gen photon
#   keep : "+gamma" side, keep events WITH an isolated gen photon
#   none : no overlapping partner sample in this analysis (signal, Higgs samples, ...)
# Adding a new MC sample = adding its name here on purpose. Override for quick studies only:
#   HZA_MC_OVERLAP_ALLOW_UNCLASSIFIED=1 (logs an error and applies no removal).
OVERLAP_VETO = ("DYto2L", "DYto2E", "DYto2Mu", "DYto2Tau", "DYJetsToLL", "EWKZ2Jets",
                "TTTo2L2Nu", "TTto2L2Nu", "WJets", "WZ_", "WW_", "ZZ_")
OVERLAP_KEEP = ("DYGto2LG", "ZGToLLG", "ZGamma2J", "ZG2J", "TTGJets", "TTG-1Jets",
                "WGTo", "WZG_", "WWG_", "ZZG_")
# TTGG (tt+gamma gamma) and ZGG are deliberately "none": their partner would be TTG/DYG with a
# second photon, which this veto does not model.
OVERLAP_NONE = ("HZa-Zto2L-ato2G", "HZa_Zto2L_ato2G", "HtoZG", "HToZG", "HtoZGto2LG", "Hto2Mu",
                "HToMuMu", "TTGG_", "ZGG_", "TTJets_", "TTtoLNu2Q", "TGJets", "ttZJets", "LLAJJ_EWK",
                "HZaTo2l2g", "/MLNanoAOD/", "ZZGTo4L")   # MLNanoAOD = private HZa_merged signal


def classify_overlap_role(file):
    """Return 'veto', 'keep' or 'none' for an MC input file; raise if it is in no list.
    Same precedence as the pre-2026-09-29 if/elif chain: veto patterns are tested first."""
    if any(p in file for p in OVERLAP_VETO):
        return "veto"
    if any(p in file for p in OVERLAP_KEEP):
        return "keep"
    if any(p in file for p in OVERLAP_NONE):
        return "none"
    msg = ("[MCOverlapTagger] MC file matches no overlap role (veto/keep/none): %s -- add its "
           "dataset name to OVERLAP_VETO/OVERLAP_KEEP/OVERLAP_NONE in mc_overlap_tagger.py" % file)
    if os.environ.get("HZA_MC_OVERLAP_ALLOW_UNCLASSIFIED") == "1":
        logger.error(msg + " (HZA_MC_OVERLAP_ALLOW_UNCLASSIFIED=1: no overlap removal applied)")
        return "none"
    raise RuntimeError(msg)
