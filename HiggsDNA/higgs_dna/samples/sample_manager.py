import os
import glob
import json
import copy
from tqdm import tqdm 

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.samples.sample import Sample
from higgs_dna.samples.file import File
from higgs_dna.utils.misc_utils import load_config 
from higgs_dna.utils import metis_utils, misc_utils

class SampleManager():
    """

    """
    def __init__(self, sample_list, years, catalog = "metadata/samples_catalog.json"):
        self.sample_list = sample_list
        print("catalog: ", catalog)
        if not isinstance(self.sample_list, list):
            self.sample_list = [self.sample_list]
        self.years = years
        if not isinstance(self.years, list):
            self.years = [self.years]

        self.catalog_name = misc_utils.expand_path(catalog)
        self.catalog = load_config(catalog)

        self.samples = {}
        self.loaded_samples = False
        self.data = None
        self.process_id_map = {}


    def get_samples(self):
        if self.loaded_samples:
            return self.data
        samples = []
        # Loop through each sample

        for s_idx, sample in enumerate(tqdm(self.sample_list)):
            if sample not in self.catalog.keys():
                logger.exception("[SampleManager : get_samples] Could not find sample '%s' in samples catalog." % (sample))
                raise ValueError()

            self.samples[sample] = {}
            info = self.catalog[sample]
            # Loop through years and create a separate Sample object for each year
            for year in self.years:
                # Skip unsupported sample/year pairs before looking up xs/bf.
                if year not in info["files"].keys():
                    logger.warning("[SampleManager : get_samples] Could not find any information about 'files' in samples catalog for sample '%s', year '%s'." % (sample, year))
                    continue

                # Get xs and bf info
                if "xs" in info.keys():
                    if isinstance(info["xs"], dict): # different xs for different years
                        if year not in info["xs"].keys():
                            logger.warning("[SampleManager : get_samples] Could not find any information about 'xs' in samples catalog for sample '%s', year '%s'. Skipping this sample/year pair." % (sample, year))
                            continue
                        xs = info["xs"][year]
                    else:
                        xs = info["xs"]

                    if "bf" in info.keys():
                        if isinstance(info["bf"], dict): # different bf for different years
                            if year not in info["bf"].keys():
                                logger.warning("[SampleManager : get_samples] Could not find any information about 'bf' in samples catalog for sample '%s', year '%s'. Skipping this sample/year pair." % (sample, year))
                                continue
                            bf = info["bf"][year]
                        else:
                            bf = info["bf"]
                    else:
                        bf = 1.

                else:
                    xs = None
                    bf = None

                if xs is not None:
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', found xs of %.6f pb and bf of %.10f" % (sample, year, xs, bf))
                    is_data = False
                else:
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s' no 'xs' info was found, treating this as data." % (sample, year))
                    is_data = True
                
                # Get input files
                files = []

                self.samples[sample][year] = {}

                grabbed_files = False

                if isinstance(info["files"][year], str) and "*" in info["files"][year]: # directory with wildcards
                    files = self.get_files_from_wildcard(info["files"][year], is_data)
                    grabbed_files = True

                # Is this a list? Could be a list of hard-coded files (Option 1) or list of dirs (Option 2b) 
                elif isinstance(info["files"][year], list):
                    if info["files"][year][0].endswith(".root"): # Option 1
                        logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', getting files from hard-coded list." % (sample, year))
                        files = [File(name = x, is_data = is_data) for x in info["files"][year]]
                        grabbed_files = True


                if not grabbed_files:
                    if isinstance(info["files"][year], str):  # Option 2a/3a
                        info["files"][year] = [info["files"][year]]  # recast this as Option 2b/3b format

                    for path in info["files"][year]:
                        files_dir = []  # ensure defined for all branches

                        if path.startswith("root://"):  # access via xrd
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a directory to be accessed via xrd" % (sample, year, path))

                            proxy = misc_utils.check_proxy()
                            if proxy is None:
                                logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired. Since you are accessing files through xrd, a valid proxy is necessary.")
                                raise RuntimeError()
                            files_dir = self.get_files_from_xrd(path, is_data)

                        elif os.path.exists(path):  # Option 1: local
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a local directory, whose files will be grabbed with <glob>." % (sample, year, path))
                            files_dir = self.get_files_from_local_dir(path, is_data)

                        elif path.endswith(("NANOAOD", "NANOAODSIM", "USER")):  # Option 3: DAS
                            logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', we interpreted the specified files '%s' as a DAS dataset, whose files will be grabbed with <dasgoclient>." % (sample, year, path))

                            proxy = misc_utils.check_proxy()
                            if proxy is None:
                                logger.exception("[CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found to be expired. Since you are accessing files through dasgoclient, a valid proxy is necessary.")
                                raise RuntimeError()
                            files_dir = self.get_files_from_dasgoclient(path, is_data)

                        else:
                            logger.warning(
                                "[SampleManager : get_samples] Unrecognized files path '%s' for sample '%s', year '%s'. "
                                "Expected: xrootd path (root://...), existing local dir, or DAS dataset ending with NANOAOD/NANOAODSIM/USER. Skipping.",
                                path, sample, year
                            )
                            continue

                        files += files_dir
                    grabbed_files = True

                files = sorted(files, key = lambda x : x.name)
                # files = files[1:] # remove first file (usually the smallest one)
                logger.debug("[SampleManager : get_samples] For sample '%s', year '%s', found %d input files:" % (sample, year, len(files)))
                if len(files) < 50: # don't print out if more than 50
                    for file in files:
                        logger.debug("\t %s" % file.name)

                # Check for sample-specific systematics
                if "systematics" in info.keys():
                    if year in info["systematics"].keys(): # year-specific syst
                        systematics = info["systematics"][year]
                    else: # same systs for all years
                        systematics = info["systematics"]
                    logger.debug("[SampleManager : get_samples] For sample '%s', year '%s' adding sample-specific systematics: %s" % (sample, year, str(systematics)))
                else:
                    systematics = None

                self.samples[sample][year] = {
                        "files" : files,
                        "xs" : xs,
                        "bf" : bf,
                }
                if systematics is not None:
                    self.samples[sample][year]["systematics"] = systematics

                # Check if fpo is specified
                if "fpo" in info.keys():
                    fpo = info["fpo"]
                else:
                    fpo = None

                samples.append(
                        Sample(
                            process = sample,
                            year = year,
                            files = files,
                            xs = xs,
                            bf = bf,
                            process_id = s_idx,
                            fpo = fpo,
                            systematics = systematics
                        )
                )
                self.process_id_map[sample] = s_idx
        self.data = samples
        self.loaded_samples = True

        self.update_catalog(samples)

        return samples

    def update_catalog(self, samples):
        """
        Rewrite the input catalog with the full list of files written out explicitly
        This way it is not necessary to rerun the glob/xrootd/dasgoclient commands and we save some time.
        """

        if "_sample_manager_full.json" in self.catalog_name: # we already made the full sample list on a previous run
            return

        catalog_full = copy.deepcopy(self.catalog)

        for sample in samples:
            catalog_full[sample.process]["files"][sample.year] = [x.name for x in sample.files]

        with open(self.catalog_name.replace(".json", "_sample_manager_full.json"), "w") as f_out:
            json.dump(catalog_full, f_out, indent = 4)

    # -----------Read T2_CN_Beijing private mc-------------------
    def _guess_dbs_instance(self, dataset: str) -> str:
        """
        Heuristic:
          - USER datasets live in prod/phys03
          - Most centrally-produced datasets live in prod/global
        """
        # dataset is "three-slash" format: /A/B/C
        if dataset.rstrip("/").endswith("/USER"):
            return "prod/phys03"
        return "prod/global"

    def _run_das_query_json(self, query_str: str):
        """
        Run dasgoclient with -json and return parsed JSON.
        """
        cmd = f"/cvmfs/cms.cern.ch/common/dasgoclient -query '{query_str}' -json"
        out = metis_utils.do_cmd(cmd)
        if not out:
            return []

        payload = self._extract_json_payload(out)
        try:
            return json.loads(payload) if payload else []
        except json.JSONDecodeError:
            logger.error(
                "[SampleManager] Failed to parse DAS JSON for query '%s'. Raw output was:\n%s",
                query_str,
                out,
            )
            raise

    def _extract_json_payload(self, raw_output: str) -> str:
        """
        dasgoclient output is not always a clean standalone JSON document. In
        practice it may contain leading warnings, a preceding "null", or even
        multiple concatenated JSON documents. Find the first decodable JSON
        object/array and return exactly that slice. If DAS only emitted scalar
        JSON values such as "null", treat it as an empty result.
        """
        raw_output = raw_output.strip()
        if not raw_output:
            return raw_output

        decoder = json.JSONDecoder()
        start = 0

        while start < len(raw_output):
            try:
                value, end = decoder.raw_decode(raw_output, start)
            except json.JSONDecodeError:
                start += 1
                continue

            if isinstance(value, (list, dict)):
                return raw_output[start:end]

            start = max(end, start + 1)

        return ""
    # -----------Read T2_CN_Beijing private mc-------------------

    def get_files_from_dasgoclient(self, sample, is_data, instance=None, redirector="root://xrootd-cms.infn.it"):
        """
        Get list of files, along with number of events and size in GB for each.

        Parameters
        ----------
        sample : str
            DBS dataset in three-slash format, e.g.
            "/EGamma/Run2018A.../NANOAOD" or "/.../.../USER"
        is_data : bool
        instance : str or None
            "prod/global" or "prod/phys03". If None, auto-guess.
        redirector : str
            xrootd redirector prefix, e.g. "root://xrootd-cms.infn.it".
            Overridable via env HDNA_REDIRECTOR (e.g. root://cmsxrootd.fnal.gov) to
            steer reads away from a flaky site (RAL) during backfill; dormant unless set.
        """
        redirector = os.environ.get("HDNA_REDIRECTOR", redirector)
        # 1) decide instance
        guessed = self._guess_dbs_instance(sample)
        inst = instance or guessed

        # 2) build query (NOTE: instance must be inside query string for your dasgoclient)
        #    Use "system=dbs" if you want DBS-only; but keep default to let DAS aggregate.
        query_str = f"file dataset={sample} instance={inst}"

        # 3) run
        query = self._run_das_query_json(query_str)

        # 4) fallback: if user didn't specify instance and nothing returned, try the other instance
        if (not query) and (instance is None):
            alt = "prod/phys03" if inst == "prod/global" else "prod/global"
            logger.warning(f"[SampleManager] No files with instance={inst}, retry with instance={alt} for sample={sample}")
            query_str_alt = f"file dataset={sample} instance={alt}"
            query = self._run_das_query_json(query_str_alt)
            inst = alt  # keep for logging

        files = []
        for j in query:
            # Expected schema: j["file"][0] is a dict with keys: name, nevents, size
            f = j.get("file", [None])[0]
            if not f:
                continue

            lfn = f.get("name")
            nevt = f.get("nevents", 0)
            size = f.get("size", 0)

            if not lfn:
                continue

            # LFN -> PFN via redirector (works for /store/... including /store/user/...)
            files.append(
                File(
                    name=(redirector.rstrip("/") + "/" + (lfn if lfn.startswith("/") else "/" + lfn)),
                    is_data=is_data,
                    n_events=int(nevt) if nevt is not None else 0,
                    size_gb=round(float(size) * 1e-9, 2) if size is not None else 0.0,
                )
            )

        for das_file in files:
            logger.debug(
                "[SampleManager:get_files_from_dasgoclient] instance=%s file=%s size=%.2f GB events=%d",
                inst, das_file.name, das_file.size_gb, das_file.n_events
            )

        return files

    def get_files_from_wildcard(self, path, is_data):
        """

        """
        files_dir = glob.glob(path)
        logger.debug("[SampleManager : get_files_from_wildcard] Found %d files from specified wildcard '%s'." % (len(files_dir), path))
        files = []
        for f in files_dir:
            if ".root" not in f:
                logger.debug("[SampleManager : get_files_from_wildcard] File '%s' was grabbed under your specified wildcard '%s' but is not a .root file, so we are skipping it." % (f, path))
                continue

            files.append(
                    File(
                        name = f,
                        is_data = is_data,
                        size_gb = round(os.path.getsize(f)*1e-9,2)
                    )
            )

        return files


    def get_files_from_local_dir(self, directory, is_data):
        """

        """
        #files_dir = glob.glob(info["files"][year] + "/*.root")
        files_dir = glob.glob(directory + "/*.root")
        logger.debug("[SampleManager : get_files_from_local_dir] Found %d files in dir '%s'." % (len(files_dir), directory))
        files = []
        for f in files_dir:
            files.append(
                    File(
                        name = f,
                        is_data = is_data,
                        size_gb = round(os.path.getsize(f)*1e-9,2)
                    )
            )

        return files


    def get_files_from_xrd(self, directory, is_data):
        """
        Get list of files from an xrootd directory with xrdfs ls command"

        :param directory: directory (in xrootd format) to glob files from
        :type directory: str
        :return: list of all root files from directory
        :rtype: list of str
        """
        files = []

        idx = directory.find("//store") + 1
        redirector = directory[:idx]
        dir = directory[idx:]

        command = "xrdfs %s ls %s" % (redirector, dir)

        logger.debug("[SampleManager : get_files_from_xrd] We will find files for dir '%s' with the command: \n %s" % (directory, command))

        contents = os.popen(command).read().split("\n")
        for x in contents:
            if x.endswith(".root"):
                files.append(File(name = redirector + x, is_data = is_data))

        logger.debug("[SampleManager : get_files_from_xrd] Found %d files in dir '%s'." % (len(files), directory))
            
        return files
