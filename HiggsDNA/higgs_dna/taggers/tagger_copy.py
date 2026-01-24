import awkward
import vector
import json
import logging
import numpy
from higgs_dna.utils.logger_utils import simple_logger
# logger = logging.getLogger(__name__)
logger = simple_logger(__name__)
import json
import re

from higgs_dna.utils import misc_utils

from higgs_dna.constants import NOMINAL_TAG

class Tagger():
    """
    Abstract base class for taggers
    :param name: Name to identify this tagger
    :type name: str
    :param options: options with selection info
    :type options: str, dict
    :param is_data: whether this tagger is being run on data
    :type is_data: bool
    :param year: which year this tagger is being run on
    :type year: str
    """
    def __init__(self, name = "tagger", options = {}, is_data = None, year = None):
        self.name = name
        self.options = options
        self.is_data = is_data
        self.year = year

        self.selection = {}
        self.events = {}
        self.cut_summary = {}

        self.options = misc_utils.load_config(options)
        self.cutflow_order = {}   # cut_type -> [cut1, cut2, ...]
        self.cutflow_yields = {}  # cut_type -> {cut: yields}

        # NEW: cutflow 輸出控制（只影響印出，不影響計數/儲存）
        cfopt = (self.options or {}).get("cutflow", {})
        self.cutflow_dump_enable = bool(cfopt.get("dump_enable", True))
        # allow list：若非空，只 dump 這些 cut_type
        self.cutflow_dump_allow = set(cfopt.get("dump_allow", []))
        # deny list：永遠不 dump 這些 cut_type
        self.cutflow_dump_deny = set(cfopt.get("dump_deny", []))
        # prefix deny：cut_type 以這些前綴開頭就不 dump（適合 trigeff_*）
        self.cutflow_dump_deny_prefix = tuple(cfopt.get("dump_deny_prefix", []))

        # NEW: dump 風格控制
        # style: "table"（逐行）或 "one_line"（單行摘要）
        self.cutflow_dump_style = str(cfopt.get("dump_style", "table"))
        # 這些前綴的 cut_type 強制用 one_line（預設 trigeff_）
        self.cutflow_dump_one_line_prefix = tuple(cfopt.get("dump_one_line_prefix", ["trigeff_"]))

        # NEW: one_line JSON (pt-binned trigeff aggregation) 去重：
        # 以前用 prefix 去重會導致「只有第一個 syst 印得出來」，其他 syst 被擋掉
        # 改成以 (syst, prefix) 去重：每個 syst 各印一次
        self._trigeff_json_dumped = set()

    def select(self, events):
        """
        Convenience function for running tagger in standalone contexts.

        :param events: awkward array of events
        :type events: awkward.Array
        :returns: awkward array of selected events
        :rtype: awkward.Array
        """
        self.current_syst = NOMINAL_TAG
        selection, events_updated = self.get_selection(NOMINAL_TAG, events)
        return events_updated[selection]


    def run(self, events, syst_tag = NOMINAL_TAG):
        """
        Return dictionary of boolean arrays of events to be selected
        by this tagger for each systematic variation,
        along with updated dictionary of events containing any additional fields
        added by this tagger.
        :param events: sets of events to perform selection for (nominal + any independent collections for syst variations
        :type events: dict
        :return: boolean arrays of events to be selected by this tagger for each systematic variation, events dict with added fields computed by this tagger
        :rtype: dict, dict 
        """

        self.current_syst = syst_tag
        if not len(events) >= 1:
            logger.debug("[Tagger] %s : event set : %s : 0 events passed to tagger, skipping running this tagger." % (self.name, syst_tag))
            selection = awkward.ones_like(events, dtype=bool)
            self.selection[syst_tag] = selection

        else:
            selection, events = self.calculate_selection(events)
            logger.debug("after selection",events)
            self.selection[syst_tag] = selection
            logger.debug("[Tagger] %s : event set : %s : %d (%d) events before (after) selection" % (self.name, syst_tag, len(selection), awkward.sum(selection)))
        
        
        # NEW: 依照配置決定哪些 cut_type 要 dump，避免 trigger study 洗版
        if self.cutflow_dump_enable:
            for cut_type in self.cutflow_order:
                if not self._should_dump_cut_type(cut_type):
                    continue
                style = self.cutflow_dump_style
                if self.cutflow_dump_one_line_prefix and cut_type.startswith(self.cutflow_dump_one_line_prefix):
                    style = "one_line"
                self.dump_cutflow(cut_type, style=style)

        return selection, events

    def _should_dump_cut_type(self, cut_type: str) -> bool:
        if cut_type in self.cutflow_dump_deny:
            return False
        if self.cutflow_dump_deny_prefix and cut_type.startswith(self.cutflow_dump_deny_prefix):
            return False
        if self.cutflow_dump_allow:
            return cut_type in self.cutflow_dump_allow
        return True

    def _trigeff_parse_cut_type(self, cut_type: str):
        """
        trigeff_{lep}_{ord}_{trig}_{ptXtoY}
        trig 允許包含底線，例如 double_ele / OR_ele
        """
        m = re.match(
            r"^trigeff_(?P<lep>[^_]+)_(?P<ord>[^_]+)_(?P<trig>.+?)_(?P<ptbin>pt\d+to(?:\d+|Inf))$",
            cut_type,
        )
        return m.groupdict() if m else None

    def _trigeff_build_summary_json(self, prefix_key: str):
        """
        prefix_key: "trigeff_ele_lead_double_ele" 之類
        回傳 dict：包含每個 pt bin 的 in_bin/pass/eff，與 overall。
        """
        bins = []
        overall_in = 0.0
        overall_pass = 0.0

        for ct in self.cutflow_order.keys():
            if not ct.startswith(prefix_key + "_"):
                continue
            meta = self._trigeff_parse_cut_type(ct)
            if not meta:
                continue

            ymap = self.cutflow_yields.get(ct, {})
            in_y = float(ymap.get("in_bin", 0.0))
            pass_y = float(ymap.get("pass_trigger", 0.0))
            eff = (pass_y / in_y) if in_y > 0 else 0.0

            bins.append({
                "ptbin": meta["ptbin"],
                "in_bin": in_y,
                "pass_trigger": pass_y,
                "eff": round(eff, 3),
            })
            overall_in += in_y
            overall_pass += pass_y

        # sort by low edge number
        def _low_edge(b):
            m = re.match(r"^pt(?P<lo>\d+)to", b["ptbin"])
            return int(m.group("lo")) if m else 10**9

        bins.sort(key=_low_edge)

        out = {
            "tagger": self.name,
            "syst": self.current_syst,
            "cut_type_prefix": prefix_key,
            "bins": bins,
            "overall": {
                "in_bin": overall_in,
                "pass_trigger": overall_pass,
                "eff": round((overall_pass / overall_in), 3) if overall_in > 0 else 0.0,
            }
        }
        return out

    def set_selection(self, selection):
        """
        Update selection (TagSequence will do this to remove overlap
        with other tags)
        :param selection: boolean array for each set of events
        :type selection: dict
        """
        for name in selection.keys():
            self.selection[name] = selection[name]


    def get_selection(self, syst_tag, syst_events): 
        """
        Return boolean array of events to be selected for a given event set.
        Checks if the selection has already been calculated and
        calculates it if not.
        :param syst_tag: name of current systematic variation with independent collection
        :type syst_tag: str
        :param syst_events: events for current systematic variation with independent collection
        :type syst_events: awkward.Array
        """

        if syst_tag in self.selection.keys():
            return self.selection[syst_tag], self.events[syst_tag] 

        else:
            return self.calculate_selection(syst_events) 

        
    def calculate_selection(self, events): 
        """
        Abstract function that should be reimplemented for each tagger.
        """

        raise NotImplementedError()


    @staticmethod
    def get_range_cut(array, ranges):
        """
        Return mask corresponding to selecting elements
        in <array> that are within the ranges (inclusive)
        specified in <ranges>.
        :param array: array of quantity to be cut on
        :type array: awkward.Array
        :param ranges: list of 2-element lists corresponding to allowed ranges of variable in <array>
        :type ranges: list
        """
        cuts = []

        for range in ranges:
            if not len(range) == 2:
                message = "Range of allowed values should be 2 numbers ([lower bound, upper bound]), but you gave %s" % str(range)
                logger.exception(message)
                raise AssertionError(message)
            cut_low = array >= range[0]
            cut_high = array <= range[1]
            cuts.append(cut_low & cut_high)

        for idx, cut in enumerate(cuts):
            if idx == 0:
                final_cut = cut
            else:
                final_cut = final_cut | cut

        return final_cut
    
    def register_event_cuts(self, names, results, events, cut_type = "event", weighted = False):
        """
        Record a given cut in the tagger instance.
        cut_type could be event-level (default) or object-level ("photon", "muon", etc)
        :param names: names to identify cuts
        :type names: list of str or str
        :param results: boolean arrays with results of applying given cut
        :type results: list of awkward.Array
        :param cut_type: flag to indicate the type of cut
        :type cut_type: str, optional
        """
        if cut_type not in self.cut_summary.keys():
            self.cut_summary[cut_type] = {}

        if not isinstance(names, list):
            names = [names]
        if not isinstance(results, list):
            results = [results]

        for name, result in zip(names, results):
            if awkward.count(result) > 0:
                yields = 0
                individual_eff = 0
                if weighted:
                    # use SF-corrected central event weight
                    # (assumes weight_central exists for MC; for safety, fall back to Generator_weight)
                    wbranch = "weight_central" if hasattr(events, "weight_central") else "Generator_weight"
                    w_all = getattr(events, wbranch)

                    # sum of weights after this cut
                    yields = float(numpy.sum(w_all[result].to_numpy().astype("float64")))

                    # standard weighted selection efficiency = sumw(pass) / sumw(total)
                    sumw_total = float(numpy.sum(w_all.to_numpy().astype("float64")))
                    individual_eff = (yields / sumw_total) if sumw_total != 0.0 else 0.0
                else:
                    yields = numpy.sum(result)
                    if float(awkward.count(result)) > 0:
                        individual_eff = awkward.sum(result) / awkward.count(result)
            self.cut_summary[cut_type][name] = {
                    "individual_eff" : float(individual_eff) 
                    #TODO: add eff as N-1 cut
            }

            # logger.debug("[Tagger] : %s, syst variation : %s, cut type : %s, cut : %s, yields : %.3f"
            #         % (self.name, self.current_syst, cut_type, name, yields))

            if cut_type not in self.cutflow_order:
                self.cutflow_order[cut_type] = []
                self.cutflow_yields[cut_type] = {}

            if name not in self.cutflow_yields[cut_type]:
                self.cutflow_order[cut_type].append(name)

            self.cutflow_yields[cut_type][name] = float(yields)


    def register_cuts(self, names, results, cut_type = "event"):
        """
        Record a given cut in the tagger instance.
        cut_type could be event-level (default) or object-level ("photon", "muon", etc)
        :param names: names to identify cuts
        :type names: list of str or str
        :param results: boolean arrays with results of applying given cut
        :type results: list of awkward.Array
        :param cut_type: flag to indicate the type of cut
        :type cut_type: str, optional
        """
        if cut_type not in self.cut_summary.keys():
            self.cut_summary[cut_type] = {}

        if not isinstance(names, list):
            names = [names]
        if not isinstance(results, list):
            results = [results]

        for name, result in zip(names, results):
            if awkward.count(result) > 0:
                individual_eff = float(awkward.sum(result)) / float(awkward.count(result))
                yields = float(awkward.sum(result))
            else:
                individual_eff = 0.
                yields = 0
            self.cut_summary[cut_type][name] = {
                    "individual_eff" : float(individual_eff) 
                    #TODO: add eff as N-1 cut
            }

            # logger.debug("[Tagger] : %s, syst variation : %s, cut type : %s, cut : %s, yields : %d"
                    # % (self.name, self.current_syst, cut_type, name, yields))
            if cut_type not in self.cutflow_order:
                self.cutflow_order[cut_type] = []
                self.cutflow_yields[cut_type] = {}

            if name not in self.cutflow_yields[cut_type]:
                self.cutflow_order[cut_type].append(name)

            self.cutflow_yields[cut_type][name] = float(yields)




    def dump_cutflow(self, cut_type, style="table"):
        if cut_type not in self.cutflow_order:
            return
        # NEW: 保險起見，即使外部直接呼叫 dump_cutflow 也會套用過濾
        if hasattr(self, "_should_dump_cut_type") and not self._should_dump_cut_type(cut_type):
            return

        order = self.cutflow_order[cut_type]
        ymap  = self.cutflow_yields[cut_type]

        if not order:
            return

        first_cut = order[0]
        last_cut = "all cuts" if "all cuts" in ymap else order[-1]
        all_y = float(ymap.get(first_cut, 0.0))
        final_y = float(ymap.get(last_cut, 0.0))
        eff = (final_y / all_y) if all_y > 0 else 0.0

        # --- one-line summary (for trigger studies etc.) ---
        if style == "one_line":
            logger.debug(
                "[CutFlow] tagger=%s syst=%s cut_type=%s final=%s %.2f / all %.2f (=%.4f)",
                self.name, self.current_syst, cut_type, last_cut, final_y, all_y, eff
            )

            # NEW: one_line 也輸出 JSON（針對 trigeff_* 匯總所有 pt bins；每個 prefix 印一次）
            if isinstance(cut_type, str) and cut_type.startswith("trigeff_"):
                meta = self._trigeff_parse_cut_type(cut_type)
                if meta:
                    prefix = f"trigeff_{meta['lep']}_{meta['ord']}_{meta['trig']}"
                    dump_key = (str(self.current_syst), prefix)
                    if dump_key not in self._trigeff_json_dumped:
                        payload = self._trigeff_build_summary_json(prefix)
                        logger.info("CutFlow JSON %s", json.dumps(payload, separators=(",", ":")))
                        logger.debug("-" * 80)
                        self._trigeff_json_dumped.add(dump_key)

            return

        # --- table style (legacy) ---
        logger.debug("-" * 80)
        logger.debug("[CutFlow] tagger=%s syst=%s cut_type=%s", self.name, self.current_syst, cut_type)
        prev = None
        for cut in order:
            y = float(ymap.get(cut, 0.0))
            cum = (y / all_y * 100.0) if all_y > 0 else 0.0
            if prev is None:
                step = None
            else:
                step = (y / prev * 100.0) if prev > 0 else 0.0

            if step is None:
                logger.debug(f"{cut:<22} : {y:>10.2f}   (all: {cum:>6.1f}%)   (step: {'-':>6})")
            else:
                logger.debug(f"{cut:<22} : {y:>10.2f}   (all: {cum:>6.1f}%)   (step: {step:>6.1f}%)")
            prev = y

        # JSON output（保留；改成兩位小數）
        payload = {
            "tagger": self.name,
            "syst": self.current_syst,
            "cut_type": cut_type,
            "cuts": {k: round(float(ymap.get(k, 0.0)), 2) for k in order},
        }
        logger.info("CutFlow JSON %s", json.dumps(payload, separators=(",", ":")))
        logger.debug("-" * 80)


    def get_summary(self):
        """
        Get dictionary of tag diagnostic info.
        :return: A dictionary of tag diagnostic info
        :rtype: dict
        """
        summary = {
                "options" : self.options,
                "cut_summary" : self.cut_summary
        }
        return summary
