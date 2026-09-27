import json, sys
from higgs_dna.utils.logger_utils import setup_logger
cfg = json.load(open(sys.argv[1]))
setup_logger('INFO')
from higgs_dna.analysis import run_analysis
run_analysis(cfg)
