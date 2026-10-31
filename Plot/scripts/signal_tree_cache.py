"""
Signal ROOT 讀檔快取：每個 (檔案, tree) 的每個 branch 在同一個 process 內只真正讀一次。

signal_eff_sumw / plot_MVASigEffVmA / plot_mAmigrated{Hist,Bar,Matrix} 的各個 accumulate 函式
原本各自 uproot.open 同一批 70 個 signal 檔（逐年、逐 channel、逐 mA、Run3 再重讀一遍），
訊號 reweight 還用 t.arrays(library="pd") 讀整棵 tree 的 379 個 branch。

open_cached(fp) 回傳與 uproot file 同介面的替身（with / in / [tree] / keys / arrays / iterate），
呼叫端只要把 uproot.open(fp) 換成 open_cached(fp)，其餘程式不必改：
  - 第一次要某些 branch 時才讀，讀進來的 awkward array 留在記憶體，之後直接取用；
    同時把 signal reweight 會用到的欄位（rw_columns）一起讀，避免之後再開一次檔。
  - iterate(step_size=...) 先確認 uproot 會不會切成多塊（見 _single_chunk）；
    只有一塊時從快取給（與 uproot 給的同一塊）；超過一塊就交回真正的 uproot iterate，切法不變。
  - 開檔失敗不快取，下次照樣重試、照樣在 open_cached 丟例外（與每次 uproot.open 的行為相同）。
  - 檔案在 process 結束前保持開啟；輸入檔若在執行中被改寫，快取不會察覺（與一次讀完再算的做法相同）。
"""
import os

import awkward as ak
import uproot

_FILES = {}


def rw_names(rw):
    """
    sideband reweight（SidebandReweighter.weights_for_dataframe 與各腳本的別名補欄）可能查的所有欄位名稱。
    weights_for_dataframe 只按名稱查欄位（param 由呼叫端先算好放進 frame），所以 frame 裡只要有這些欄位，
    權重就與讀整棵 tree 的 frame 相同。遇到不認得的衍生變數回傳 None（= 只能讀整棵 tree）。
    """
    from sideband_reweight import VAR_ALIASES, DERIVED_VARS
    names = {"n_electrons", "n_muons",
             "pho1ECALIso", "pho1PIso_noCorr", "pho2ECALIso", "pho2PIso_noCorr",
             "H_m", "H_mass", "ALP_m", "ALP_mass"}
    for iteration in rw.iterations:
        for step in iteration.get("steps", []):
            var = step["var"]
            names.add(var)
            names.update(VAR_ALIASES.get(var, ()))
            if var in DERIVED_VARS:
                if var != "pho_pt_asym":
                    return None
                for base in ("pho1Pt_oHm", "pho2Pt_oHm"):
                    names.add(base)
                    names.update(VAR_ALIASES.get(base, ()))
    return names


def rw_columns(rw, keys, extra=()):
    """rw_names 中 tree 裡有的欄位（加上 extra，例如 weight）；None = 讀整棵 tree。"""
    names = rw_names(rw)
    if names is None:
        return None
    keys = set(map(str, keys))
    return sorted(n for n in names | set(extra) if n in keys)


# 第一次讀某棵 tree 時一併讀進來的欄位名稱（實際只讀 tree 裡有的）；腳本用 set_prefetch 設定
PREFETCH = set()


def set_prefetch(rw, extra=()):
    names = rw_names(rw)
    if names is not None:
        PREFETCH.update(names | set(extra))


class _CachedTree:
    def __init__(self, cfile, name):
        self._cfile = cfile
        self._name = name
        self._keys = None
        self._arrs = {}

    def _real(self):
        return self._cfile._real()[self._name]

    @property
    def num_entries(self):
        return self._real().num_entries

    def keys(self):
        if self._keys is None:
            self._keys = list(self._real().keys())
        return list(self._keys)

    def _load(self, branches):
        branches = list(dict.fromkeys(branches))
        missing = [b for b in branches if b not in self._arrs]
        if missing:
            keys = set(self.keys())
            extra = [b for b in sorted(PREFETCH) if b in keys and b not in self._arrs and b not in missing]
            try:
                rec = self._real().arrays(missing + extra, library="ak")
            except Exception:
                if not extra:
                    raise
                rec = self._real().arrays(missing, library="ak")
            for b in rec.fields:
                self._arrs[b] = rec[b]
        return {b: self._arrs[b] for b in branches}

    def arrays(self, branches, library="ak"):
        d = self._load(branches)
        if library == "ak":
            return ak.Array(d)
        if library == "np":
            return {b: ak.to_numpy(v) for b, v in d.items()}
        if library == "pd":
            import pandas as pd
            return pd.DataFrame({b: ak.to_numpy(v) for b, v in d.items()})
        raise ValueError(f"unsupported library {library!r}")

    def _single_chunk(self, branches, step_size, n):
        """uproot iterate 是否只切一塊。uproot 用「所選 branch 的 compressed basket bytes」換算每塊的 entry 數，
        這個總數不會超過檔案大小，所以檔案不大於 step_size 時必定一塊（免算；num_entries_for 要掃 379 個 branch）。"""
        try:
            target = uproot._util.memory_size(step_size)
        except Exception:
            target = None
        if target is not None and os.path.getsize(self._cfile.fp) <= target:
            return True
        return self._real().num_entries_for(step_size, expressions=list(branches)) >= n

    def iterate(self, branches, library="ak", step_size="100 MB"):
        real = self._real()
        n = real.num_entries
        if n == 0:
            return
        if self._single_chunk(branches, step_size, n):
            yield self.arrays(branches, library=library)
        else:
            yield from real.iterate(branches, library=library, step_size=step_size)


class _CachedFile:
    def __init__(self, fp):
        self.fp = fp
        self._file = None
        self._has = {}
        self._trees = {}

    def _real(self):
        if self._file is None:
            self._file = uproot.open(self.fp)
        return self._file

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __contains__(self, name):
        if name not in self._has:
            self._has[name] = name in self._real()
        return self._has[name]

    def __getitem__(self, name):
        if name not in self._trees:
            self._real()[name]  # 不存在時與 uproot 一樣丟 KeyError
            self._trees[name] = _CachedTree(self, name)
        return self._trees[name]


def open_cached(fp):
    fp = str(fp)
    if fp not in _FILES:
        cf = _CachedFile(fp)
        cf._real()  # 開檔失敗在這裡丟例外（與 uproot.open 同一個位置），且不快取
        _FILES[fp] = cf
    return _FILES[fp]
