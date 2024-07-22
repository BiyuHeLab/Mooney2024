import pickle
from datetime import datetime
from functools import partial
from pathlib import Path
from time import time
from sys import getsizeof
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

from .analysis import PWMVAnalysis
from ..paths import FilePath
from ..plotting import Subjects
from ..ROIs import pwmvROIs
from ..stats import Array, ClusterCorrectedMultipleDof, MAXWORKERS, SEM_within_subject
from ..data import load_mne_subject


class MeanROIAnalysis(PWMVAnalysis):
    """
    Mean ROI analysis class.
    """

    def __init__(
        self,
        subjects: Subjects,
        toi=slice(-0.2, 1.5),
        source: Literal["HGP", "ERP"] = "HGP",
        smoothWindowSecs=0.1,
        roiGroupName: str = "network",
        outputPath: Path = Path("..", "Analysis", "MeanROI"),
        rngSeed=20240522,
        maxWorkers=MAXWORKERS,
    ) -> None:
        super().__init__(subjects, outputPath, smoothWindowSecs=smoothWindowSecs)
        self.toi = toi
        self.source = source
        self.roiGroupName = roiGroupName
        self.saveFolder.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.default_rng(rngSeed)
        self.maxWorkers = maxWorkers

        self.rois = pwmvROIs(Path("Electrodes", "ROI_Masks", roiGroupName))

    def __sizeof__(self) -> int:
        sizeSum = (
            getsizeof(self.toi)
            + getsizeof(self.source)
            + getsizeof(self.roiGroupName)
            + getsizeof(self.rng)
        )
        if hasattr(self, "sigBars"):
            sizeSum += getsizeof(self.sigBars)
        if hasattr(self, "clusterStats"):
            sizeSum += sum(
                [getsizeof(cs) for vals in self.clusterStats.values() for cs in vals.values()]
            )
        return super().__sizeof__() + sizeSum

    @property
    def saveFolder(self) -> Path:
        # Figure output folder
        if not hasattr(self, "_saveFolder"):
            self._saveFolder = self.rootFolder / self.roiGroupName
            self._saveFolder.mkdir(parents=True, exist_ok=True)
        return self._saveFolder

    def save(self, path: Path = None):
        if hasattr(self, "data"):
            delattr(self, "data")
        super().save(path)

    def run(self, loadFromDisk=True, rng: np.random.Generator = None):
        self.load_results(loadFromDisk)
        self.load_significance_bars(loadFromDisk, rng)
        print("Results loaded.")

    def load(self, pklPath: Path = None, runIfNotSaved=True):
        self = super().load(pklPath, runIfNotSaved)
        # Reset saved paths
        if hasattr(self, "_saveFolder"):
            delattr(self, "_saveFolder")
        return self

    def load_data(self, loadFromDisk=True) -> None:
        # Verify datafile is newer than most recent patient data
        dataFilename = Path(self.dataFolder, f"MeanROI_allData_{self.source}_{self.subjects}.pkl")
        if (
            loadFromDisk
            and dataFilename.is_file()
            and (datetime.fromtimestamp(dataFilename.stat().st_mtime) > self._update_before_date)
        ):
            print(f"Loading from disk: {dataFilename}.")
            self.data: pd.DataFrame = pd.read_pickle(dataFilename)
        else:
            print(f"Generating data for {self.make_detail_str()}.")
            self.data = self._load_data()
            print(f"Saving data: {dataFilename.resolve()}.")
            self.data.to_pickle(dataFilename)
        self.data = self.data.droplevel(["TrialNum"])
        self.data.columns = self.data.columns.astype("float")  # Shouldn't be necessary
        self.data = self.data.loc[:, self.toi]
        self.elecCondMeans = self.data.groupby(["Condition", "SubjectID", "Electrode"]).mean()

    def _load_data(self) -> pd.DataFrame:
        allSubjDatas = []
        for pNum in self.subjects:
            subj = load_mne_subject(pNum, loadGamma=self.source == "HGP")
            # behavs[pNum] = subj.behav
            if self.source == "HGP" and self._smoothingWindow is not None:
                assert subj.fsample == self.dataSampleRate, (
                    "Window calculation expects HGP data to be sampled at"
                    f" {self.dataSampleRate} Hz."
                )
                subjData = subj.data.rolling(
                    self._smoothingWindow, center=True, min_periods=1, axis=1
                ).mean()
            else:
                subjData = subj.data
            trialConditions = (
                pd.from_dummies(subj.behav[["pre", "post", "grey"]]).astype("string").squeeze()
            )
            condIndex = trialConditions.loc[subjData.index.get_level_values("TrialNum")]
            condIndex.loc[~subj.behav.disambig] = "noDisambig"
            subjData = (
                subjData.assign(Condition=condIndex.values)
                .set_index("Condition", append=True)
                .drop(index="noDisambig", level="Condition")
            )
            # Fix columns.dtype back to float after momentarily adding the 'Condition' column
            subjData.columns = subjData.columns.astype("float")
            allSubjDatas.append(subjData)
        allData = pd.concat(
            allSubjDatas,
            keys=self.subjects,
            names=["SubjectID", "Electrode", "TrialNum", "Condition"],
        )
        return allData

    def load_results(self, loadFromDisk=True):
        resultsFile = self.dataFolder / f"MeanROI_results_{self.make_detail_str()}.pkl"
        if (
            loadFromDisk
            and resultsFile.is_file()
            and (datetime.fromtimestamp(resultsFile.stat().st_mtime) > self._update_before_date)
        ):
            print(f"Loading from disk: {resultsFile}.")
            self.results = pd.read_pickle(resultsFile)
        else:
            print(f"Generating results for {self.make_detail_str()}.")
            self.results = self._load_results(loadFromDisk=loadFromDisk)
            print(f"Saving results: {resultsFile.resolve()}.")
            self.results.to_pickle(resultsFile)

    def _load_results(self, loadFromDisk=True) -> pd.DataFrame:
        if not hasattr(self, "data"):
            self.load_data(loadFromDisk=loadFromDisk)
        presentROIs = self.memberships.columns[self.memberships.any()]
        roiResults = []
        for roiName in presentROIs:
            roiElecMeans = self._roi_elec_means(roiName)

            # Unpaired analysis
            roiCondData = pd.concat(
                [
                    roiElecMeans.groupby("Condition").mean(),
                    roiElecMeans.groupby("Condition").sem(),
                    roiElecMeans.groupby("Condition").count(),
                ],
                keys=["mean", "sem", "numElec"],
                names=["Stat", "Condition"],
            ).reorder_levels([1, 0])

            # Paired analysis
            roiPreElecMeans, roiPostElecMeans = self._pair_pre_post(roiElecMeans)
            # Update roiCondData with updated paired means & numElec
            roiCondData.loc[("pre", "mean"), :] = roiPreElecMeans.mean().values
            roiCondData.loc[("pre", "numElec"), :] = roiPreElecMeans.count().values
            roiCondData.loc[("post", "mean"), :] = roiPostElecMeans.mean().values
            roiCondData.loc[("post", "numElec"), :] = roiPostElecMeans.count().values
            # Calc paired SEM for Pre/Post
            prePostSEM = pd.DataFrame(
                {
                    t: SEM_within_subject(
                        np.vstack(
                            [roiPreElecMeans.loc[:, t].values, roiPostElecMeans.loc[:, t].values]
                        ).T
                    )
                    for t in roiElecMeans.columns
                },
                index=["pre", "post"],
                columns=roiElecMeans.columns,
            )
            roiCondData.loc[(["pre", "post"], "sem"), :] = prePostSEM.values
            roiResults.append(roiCondData)
        return pd.concat(roiResults, keys=presentROIs, names=["ROI", "Condition", "Stat"])

    @property
    def memberships(self) -> pd.DataFrame:
        if not hasattr(self, "_memberships"):
            self._memberships = self.rois.get_patient_memberships(self.subjects)
            # Unload Niftis to save memory
            self.rois = pwmvROIs(Path("Electrodes", "ROI_Masks", self.roiGroupName))
        return self._memberships

    def load_significance_bars(
        self,
        loadFromDisk=True,
        rng: np.random.Generator = None,
        tests=["Post>Pre", "Pre>0", "Post>0"],
    ):
        if rng is None:
            rng = self.rng
        roiMemberships = self.memberships
        presentROIs = roiMemberships.columns[roiMemberships.any()]
        with tqdm(presentROIs, unit="ROI") as roiBar:
            roiTests = []
            roiClusterStats: dict[str, dict[str, ClusterCorrectedMultipleDof]] = {}
            for roiName in roiBar:
                roiBar.set_description(roiName)
                sigPoints = []
                testClusterStats: dict[str, ClusterCorrectedMultipleDof] = {}
                for sigTest in tests:
                    sigFname = self.dataFolder / (
                        f"{roiName}_electrode_mean_{self.source}_{self.subjects}_"
                        + f"{sigTest.replace(' ','_')}.pkl"
                    )
                    if (
                        loadFromDisk
                        and sigFname.is_file()
                        and (
                            datetime.fromtimestamp(sigFname.stat().st_mtime)
                            > self._update_before_date
                        )
                    ):
                        print(f"Loading from disk: {sigFname}.")
                        with open(sigFname, "rb") as f:
                            clusterStats = pickle.load(f)
                    else:
                        print(f"Generating {sigTest} significance for {roiName}.")
                        clusterStats = self._load_significance_stats(sigTest, roiName, rng)
                        print(f"Saving clusterStats: {sigFname.resolve()}.")
                        with open(sigFname, "wb") as f:
                            pickle.dump(clusterStats, f)
                    sigPoints.append(pd.Series(clusterStats.sigClusterPlotPoints))
                    testClusterStats[sigTest] = clusterStats
                roiTests.append(pd.concat(sigPoints, keys=tests, names=["Test"], axis=1).T)
                roiClusterStats[roiName] = testClusterStats
        sigBarsData = pd.concat(roiTests, keys=presentROIs, names=["ROI"])
        sigBarsData.columns = clusterStats.refIndex
        self.sigBars = sigBarsData
        self.clusterStats = roiClusterStats

    def _roi_elec_means(self, roiName: str) -> pd.DataFrame:
        roiElecsMask = self.elecCondMeans.index.droplevel(["Condition"]).isin(
            self.memberships[self.memberships[roiName]].index
        )
        return self.elecCondMeans[roiElecsMask]

    def _pair_pre_post(self, roiElecMeans: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Create paired pre/post dataframes.
        """
        roiPreElecMeans = roiElecMeans.loc["pre"].sort_index()
        roiPostElecMeans = roiElecMeans.loc["post"].sort_index()
        # Filter out unpaired electrodes
        roiPreElecMeans = roiPreElecMeans[roiPreElecMeans.index.isin(roiPostElecMeans.index)]
        roiPostElecMeans = roiPostElecMeans[roiPostElecMeans.index.isin(roiPreElecMeans.index)]
        # Set unpaired timepoints to NaN (should not include a condition when paired is gone)
        roiPreElecMeans = roiPreElecMeans.mask(roiPostElecMeans.isnull())
        roiPostElecMeans = roiPostElecMeans.mask(roiPreElecMeans.isnull())
        return roiPreElecMeans, roiPostElecMeans

    def _load_significance_stats(
        self,
        sigTest: Literal["Post>Pre", "Pre>0", "Post>0"],
        roiName: str,
        rng: np.random.Generator,
    ) -> ClusterCorrectedMultipleDof:
        if not hasattr(self, "data"):
            self.load_data()
        # Establish testSlice so pre-stim data is not tested against post-stim data
        if self.source == "HGP":
            zeroIdx = self.data.columns.to_list().index(0)
            self.testSlice = slice(
                (
                    None
                    if self._smoothingWindow is None
                    else zeroIdx + (self._smoothingWindow - 1) // 2
                ),
                None,
            )
        else:
            self.testSlice = slice(None)
        roiElecMeans = self._roi_elec_means(roiName)
        ttest = partial(stats.ttest_rel, nan_policy="omit", alternative="greater")

        if sigTest == "Gray>0":
            roiGreyElecMeans = roiElecMeans.loc["grey"].sort_index()
            clusterStats = ClusterCorrectedMultipleDof(
                roiGreyElecMeans,
                np.zeros_like(roiGreyElecMeans),
                statisticCalculator=ttest,
                test="greater",
                testName="Grey>0",
                rng=rng,
                testWindow=self.testSlice,
                maxWorkers=self.maxWorkers,
            )
            return clusterStats

        roiPreElecMeans, roiPostElecMeans = self._pair_pre_post(roiElecMeans)
        # Update numElecs for paired electrodes only
        # numElecs = roiPreElecMeans.index.nunique()

        if sigTest == "Post>Pre":
            A, B = roiPostElecMeans, roiPreElecMeans
        elif sigTest == "Pre>0":
            A, B = roiPreElecMeans, np.zeros_like(roiPreElecMeans)
        elif sigTest == "Post>0":
            A, B = roiPostElecMeans, np.zeros_like(roiPostElecMeans)
        clusterStats = ClusterCorrectedMultipleDof(
            A,
            B,
            statisticCalculator=ttest,
            test="greater",
            testName=sigTest,
            rng=rng,
            testWindow=self.testSlice,
            maxWorkers=self.maxWorkers,
        )
        return clusterStats

    def make_detail_str(self) -> str:
        return (
            super()._basic_detail_str(roiGroup=self.roiGroupName, source=self.source)
            + f"_{self.subjects}"
        )

    def clusters_to_df(self, conds=["Pre>0", "Post>0", "Post>Pre"]) -> pd.DataFrame:
        """
        Convert clusterStats to a DataFrame.

        Args:
            conds (list[str], optional): List of conditions to include. Will be in same order in
                the DataFrame. Defaults to ["Pre>0", "Post>0", "Post>Pre"].

        Returns:
            pd.DataFrame: MultiIndex with levels: [ROI, Condition, ClusterNum]. Columns are:
                [pVal, clusterStat, startCluster, endCluster, ClusterDoF, NullSize].
        """
        if hasattr(self, "clusterStats"):
            roiClusterDFs = []
            rois = self.clusterStats.keys()
            for roi in rois:
                condClusterDFs = []
                for condition in conds:
                    clusterStats = self.clusterStats[roi][condition]
                    condClusterDFs.append(
                        pd.DataFrame({
                            "pVal": clusterStats.clusterPVals,
                            "clusterStat": [
                                clusterStats.stats[start:end].sum()
                                for start, end in clusterStats.clusterEndpoints
                            ],
                            "startCluster": [
                                clusterStats.refIndex[points[0]]
                                for points in clusterStats.clusterEndpoints
                            ],
                            "endCluster": [
                                clusterStats.refIndex[points[1] - 1]
                                for points in clusterStats.clusterEndpoints
                            ],
                            "ClusterDoF": clusterStats.clusterDoFs,
                            "NullSize": [
                                clusterStats.nullSizes[bn] for bn in clusterStats.clusterBins
                            ],
                        })
                    )
                roiClusterDFs.append(
                    pd.concat(condClusterDFs, keys=list(conds), names=["Condition", "ClusterNum"])
                )
            return pd.concat(
                roiClusterDFs, keys=list(rois), names=["ROI", "Condition", "ClusterNum"]
            )
        else:
            raise AttributeError("No clusterStats attribute found. Run load_significance_bars().")


def bootstrap_results(
    mra: MeanROIAnalysis,
    saveDir: FilePath,
    rngSeed: int | None = None,
    roiTests: dict[str, list[str]] | None = None,
    numBootstraps=1000,
    saveTime=1200,
    maxWorkers=MAXWORKERS,
) -> None:
    """
    Generate bootstrapped results of a MeanROIAnalysis. Results saved to disk.

    Allows for effect onset time estimation.

    Args:
        mra (MeanROIAnalysis): MeanROIAnalysis object.
        saveDir (FilePath): Directory to save results.
        rngSeed (int, optional): Random number generator seed. Defaults to a random seed.
        roiTests (dict[str, list[str]] | None, optional): The ROIs and tests to bootstrap. Keys are
            ROI names, values are a list of tests to bootstrap in that ROI.Defaults to all ROIs in
            the mra.memberships attribute and all tests, i.e. ["Post>Pre", "Pre>0", "Post>0"].
        numBootstraps (int, optional): Number of replications. Defaults to 1000.
        saveTime (int, optional): Results will be saved to disk no more frequently than this many
            seconds. Defaults to 1200.
    """
    rng = np.random.default_rng(rngSeed)
    Path(saveDir).mkdir(parents=True, exist_ok=True)
    startTime = time()
    if roiTests is None:
        roiTests = {roi: ["Post>Pre", "Pre>0", "Post>0"] for roi in mra.memberships.columns}
    mra.load_data()
    if mra.source == "HGP":
        zeroIdx = mra.data.columns.to_list().index(0)
        testSlice = slice(
            (None if mra._smoothingWindow is None else zeroIdx + (mra._smoothingWindow - 1) // 2),
            None,
        )
    else:
        testSlice = slice(None)
    ttest = partial(stats.ttest_rel, nan_policy="omit", alternative="greater")
    print(f"Generating bootstrap results for {roiTests}.")
    for roiName, tests in roiTests.items():
        roiElecMeans = mra._roi_elec_means(roiName)
        roiPreElecMeans, roiPostElecMeans = mra._pair_pre_post(roiElecMeans)
        for test in tests:
            savePath = Path(
                saveDir,
                f"{roiName}_electrode_mean_{mra.source}_{mra.subjects}_{test.replace(' ','_')}"
                "_bootstrap_results.pkl",
            )
            bootResults = []
            for _ in tqdm(
                range(numBootstraps),
                desc="-".join([roiName, test]),
                unit="replication",
                smoothing=0,
            ):
                bootPre = roiPreElecMeans.sample(frac=1, replace=True, random_state=rng)
                bootPost = roiPostElecMeans.loc[bootPre.index]
                if test == "Post>Pre":
                    A: Array = bootPost
                    B: Array = bootPre
                elif test == "Pre>0":
                    A, B = bootPre, np.zeros_like(bootPre)
                elif test == "Post>0":
                    A, B = bootPost, np.zeros_like(bootPost)
                clusterStats = ClusterCorrectedMultipleDof(
                    A,
                    B,
                    statisticCalculator=ttest,
                    test="paired greater" if (test == "Post>Pre") else "greater",
                    testName=test,
                    rng=rng,
                    testWindow=testSlice,
                    maxWorkers=maxWorkers,
                )
                bootResults.append(clusterStats)
                if time() - startTime > saveTime:
                    with open(savePath, "wb") as f:
                        pickle.dump(bootResults, f)
            with open(savePath, "wb") as f:
                pickle.dump(bootResults, f)
