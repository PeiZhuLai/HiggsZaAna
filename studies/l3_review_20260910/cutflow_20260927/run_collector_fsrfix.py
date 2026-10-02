"""Run HiggsDNA/scripts/5_collect_cutflow.py UNCHANGED on the provenance-checked FSR-fix .out mirror.
main() reads the module globals DEFAULT_BASE_DIR / outputDir at call time, so they are redirected
here instead of editing the script (its CUTFLOW_BASEDIR env var is not used by main())."""
import importlib.util, os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
spec=importlib.util.spec_from_file_location("cc","/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/5_collect_cutflow.py")
cc=importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
cc.DEFAULT_BASE_DIR=os.path.join(HERE,"logs_mirror")
cc.outputDir=os.path.join(HERE,"cutflow_list_fsrfix")
sys.exit(cc.main())
