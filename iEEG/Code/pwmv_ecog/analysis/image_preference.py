import logging
import pickle
from datetime import datetime
from functools import partial
from pathlib import Path
from sys import getsizeof
from time import time
from typing import Literal, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from tqdm import tqdm

from .analysis import PWMVAnalysis
from ..logger import get_dated_FileHandler, set_MNE_logfile
from ..paths import FilePath
from ..plotting import plot_sig_bar, Subjects
from ..process import patient_last_updated
from ..ROIs import pwmvROIs
from ..stats import (
    ClusterCorrectedMultipleDof,
    InsufficientPointGenerationError,
    MAXWORKERS,
    SEM_within_subject,
)
from ..data import ECoGDataMNE, load_mne_subject

polyfit = np.polynomial.Polynomial.fit

IMAGEPREFLOGS = Path("../Analysis/ImagePref/logs")

# Setup logging
logger = logging.getLogger(__name__)
logger.propagate = False
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
file_handler = get_dated_FileHandler(IMAGEPREFLOGS / "imagePref.log", debug_filelogging=True)
console_handler.setLevel(logging.WARNING)
logger.addHandler(console_handler)
logger.addHandler(file_handler)
set_MNE_logfile(file_handler.baseFilename)


def image_pref_wilcoxon_perm_test(
    data: pd.DataFrame,
    oneSideConds: list[str] = ["Pre", "Post"],
    pairedConds: list[tuple[str, str]] = [("Post", "Pre")],
    numRepeats=1000,
    cdt=0.05,
    rngSeed: int | np.random.Generator = 20220506,
    chunksize: int | None = None,
    testWindow: slice = slice(None),
    maxWorkers: int = MAXWORKERS,
) -> tuple[pd.DataFrame, dict[str, ClusterCorrectedMultipleDof]]:
    """
    Uses the wilcoxon test at each timepoint and across the electrodes within the provided
    conditions -- Pre, Post, and paired Post-Pre by default. Cluster based permutation test for
    multiple comparisons correction.

    Expects data to have MultiIndex with a 'Condition' level, with the timeseries along the
    columns. Having extra MultiIndex levels is fine, and unused.

    Returns (fullResults, allClusterStats)
    fullresults is a DataFrame with MultiIndex['Condition', 'Stats']. Condition is Pre and Post for
    the single condition compared to zero, and Post - Pre is the between condition test. Stats is
    ['statistic', 'significant', 'slope mean', 'slope SEM', 'numElec']. The columns are the
    timepoints.
    allClusterStats is a dict mapping each Condition name to the ClusterCorrectedMultipleDof object
    used to produce the stats forth condition.
    """
    rng = np.random.default_rng(seed=rngSeed)
    results = []
    allClusterStats: dict[str, ClusterCorrectedMultipleDof] = {}
    condNames = []
    for cond in oneSideConds:
        condData = data.xs(cond, level="Condition")
        condNames.append(cond)
        logger.info(f"Running cluster corrected one-sided wilcoxon on condition {cond}.")
        # Have to use zsplit to ensure all permutations use the same rank range
        clusterStats = ClusterCorrectedMultipleDof(
            condData,
            np.zeros_like(condData),
            partial(always_signedT_wilcoxon, test="greater", zero_method="zsplit"),
            "greater",
            condNames[-1],
            rng,
            numRepeats=numRepeats,
            cdt=cdt,
            chunksize=chunksize,
            testWindow=testWindow,
            maxWorkers=maxWorkers,
        )
        allClusterStats[condNames[-1]] = clusterStats
        sigPoints, wStats = clusterStats.sigClusterPlotPoints, clusterStats.stats
        results.append(
            pd.DataFrame(
                [wStats, sigPoints, condData.mean(), condData.sem(), condData.count()],
                index=pd.Index(
                    ["statistic", "significant", "slope mean", "slope SEM", "numElec"],
                    name="Stats",
                ),
                columns=condData.columns,
            )
        )
    for condA, condB in pairedConds:
        condAData = data.xs(condA, level="Condition")
        condBData = data.xs(condB, level="Condition")
        condNames.append(f"{condA} - {condB}")
        logger.info(
            f"Running cluster corrected one-sided paired wilcoxon on condition {condNames[-1]}."
        )
        # Have to use zsplit to ensure all permutations use the same rank range
        clusterStats = ClusterCorrectedMultipleDof(
            condAData,
            condBData,
            partial(always_signedT_wilcoxon, test="greater", zero_method="zsplit"),
            "paired greater",
            condNames[-1],
            rng,
            numRepeats=numRepeats,
            cdt=cdt,
            chunksize=chunksize,
            testWindow=testWindow,
            maxWorkers=maxWorkers,
        )
        allClusterStats[condNames[-1]] = clusterStats
        sigPoints, wStats = clusterStats.sigClusterPlotPoints, clusterStats.stats
        # Calculate paired SEM
        semData = {}
        for t in condAData.columns:
            tData = np.vstack([condAData[t].values, condBData[t].values]).T
            tDataNoNaNs = tData[~np.isnan(tData).any(axis=1), :]
            semData[t] = SEM_within_subject(tDataNoNaNs)
        pairedSEM = pd.DataFrame(semData, columns=condAData.columns)
        results.append(
            pd.DataFrame(
                [
                    wStats,
                    sigPoints,
                    (condAData - condBData).mean(),
                    pairedSEM.iloc[0, :].values,
                    (condAData - condBData).count(),
                ],
                index=pd.Index(
                    ["statistic", "significant", "slope mean", "slope SEM", "numElec"],
                    name="Stats",
                ),
                columns=condAData.columns,
            )
        )
    fullResults = pd.concat(
        results, axis=0, keys=condNames, names=["Condition"] + results[0].index.names
    )
    return fullResults, allClusterStats


def always_signedT_wilcoxon(
    A: np.ndarray,
    B: np.ndarray,
    test: str,
    zero_method: Literal["pratt", "wilcox", "zsplit"] = "zsplit",
    nan_policy: Literal["omit", "propagate", "raise"] = "omit",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the wilcoxon but ensure the returned statistic is always summable.

    This means the stat will be the sum of the ranks of the differences above zero for a 'greater'
    test; the negated sum of the ranks of the differences below zero for a 'lesser' test; and the
    sum of the ranks of the differences above or below zero, whichever is bigger, signed according
    to which was bigger, for the 'two-sided' test.

    Args:
        A (np.ndarray): First dataset. Should be same shape as B.
        B (np.ndarray): Second dataset. Should be same shape as A.
        test (str): Type of wilcoxon being performed. Can be 'greater', 'lesser, or 'two-sided'.
        zero_method (str, optional): {“pratt”, “wilcox”, “zsplit”},
            The following options are available (default is “zsplit”):
            - “pratt”: Includes zero-differences in the ranking process, but drops the ranks of the
                zeros, (more conservative).
            - “wilcox”: Discards all zero-differences.
            - “zsplit”: Includes zero-differences in the ranking process and split the zero rank
                between positive and negative ones. Defaults to "zsplit".
        nan_policy (Literal["omit", "propagate", "raise"], optional): Passed to scipy. Default is
            "omit".

    Returns:
        tuple[np.ndarray, np.ndarray]: _description_
    """
    assert A.shape == B.shape
    import scipy

    if scipy.__version__ == "1.13.0":
        if nan_policy != "omit":
            raise ValueError("Wilcoxon can only use the omit nan_policy in scipy version 1.13.0.")
        # Bug with scipy.stats.wilcoxon. This is the work around
        wstats, pVals = stats.wilcoxon(
            A, B, alternative=test, zero_method=zero_method, _no_deco=True
        )
    else:
        wstats, pVals = stats.wilcoxon(
            A, B, alternative=test, zero_method=zero_method, nan_policy=nan_policy
        )
    numRanks = A.shape[0] * (A.shape[0] + 1) / 2
    if test == "lesser":
        wstats = numRanks - wstats
    elif test == "two-sided":
        # scipy returns min(T-,T+) in the two-sided test, which cannot be summed to form cluster
        # stats. By using zsplit, we can rely on max(T-,T+) = numRanks - min(T-,T+)
        minT = wstats.copy()
        wstats = numRanks - minT
        signedRanks = pd.DataFrame(A - B).abs().rank() * np.sign(A - B)
        numZeros = np.sum(signedRanks.values == 0)
        assert numZeros == 0
        tPlus = signedRanks[signedRanks > 0].sum() + numZeros * (numZeros + 1) / 4
        # if T+ = min(T+,T-), then there were more negative ranks, and the stat should be negative
        wstats[minT == tPlus] *= -1
    return wstats, pVals


def tune_electrodes(
    data: ECoGDataMNE,
    imageNums: list[int],
    smoothWindowDuration: float,
    greyActivationWindow=None,
    requiredDuration=1.5,
) -> pd.DataFrame:
    # Given ECoGDataMNE and image numbers, find the best fit lines for the pre and post conditions
    # according to the greyscale activation ordering of the images.
    # smoothWindowDuration is the width of a rolling mean smoothing window in SECONDS.
    # requiredDuration is a detection threshold: if an image has no trials for some condition
    #   that last at least requiredDuration seconds, a warning is raised.
    # Returned DataFrame has MultiIndex with levels=['Electrode', 'Condition', 'Best Fit Component']
    # and columns=data.data.columns
    electrodes = data.data.index.get_level_values("Electrode").unique()
    collectResults = []
    assert data.fsample == 512, "Significance testing testWindow assumes 512Hz."
    window = int(np.floor(smoothWindowDuration * data.fsample))
    # TODO: Install numba to improve performance of this smoothing. Following line needs numba
    # smoothed = ecogData.data.rolling(window, center=True, min_periods=6, axis=1, method='table').mean(engine='numba')
    smoothed = data.data.rolling(window, center=True, min_periods=1, axis=1).mean()
    barPosition = (int(data.subjectID) % 1000) % 15
    for elec in tqdm(
        electrodes,
        unit="electrode",
        position=barPosition,
        leave=False,
        desc=f"Patient_{data.subjectID}",
    ):
        logger.info(f"Requesting image preference for electrode {elec}, Subject {data.subjectID}")
        elecData = smoothed.loc[elec]
        collectResults.append(
            electrode_tuning(elecData, data.behav, imageNums, activationTOI=greyActivationWindow)
        )

    bestFits = pd.concat(collectResults, axis=0, keys=electrodes)
    # Report anomolous data
    # Warn if any images have insufficient trial lengths
    for condition in ["Pre", "Post"]:
        condFits = bestFits.xs((condition, "slope"), level=["Condition", "Best Fit Component"])
        # Assuming each timeseries has one end-point after which all NaN values are found
        imageDurations = condFits.apply(
            lambda x: condFits.columns[x.dropna().shape[0] - 1], axis=1
        )
        if (imageDurations <= requiredDuration).any():
            tooShort = imageDurations[imageDurations < requiredDuration]
            # If this warning ever actually fires on real data, we will need to make an intentional
            # decision on how to calculate the slope at a timepoint after an image has dropped out
            logger.warning(
                f"Following image(s) and their max {condition} trial duration(s) fail to"
                f" exceed {requiredDuration}s:\n{tooShort}"
            )
    return bestFits


def get_activations_at_each_timepoint(
    dataView: pd.DataFrame,
    behav: pd.DataFrame,
    imageNums: Sequence[int],
    condition: str,
    disambigOnly=True,
) -> pd.DataFrame:
    """
    Get mean activations for each image in imageNums during given condition independently at
    each timepoint.

    Args:
        dataView (pd.DataFrame): Activation DataFrame. Index expected to be TrialNums. Columns
            expected to be timepoints.
        behav (pd.DataFrame): Behavioral DataFrame. Index expected to be TrialNums. Needs
            columns "ImageNum", "disambig", and one that matches passed condition.
        imageNums (Sequence[int]): Only these images will be included in output.
        condition (str): Only trials for this condition will be included in the output.
        disambigOnly (bool, optional): Only disambiguated trials will be included in the
            output. Defaults to True.

    Returns:
        pd.DataFrame: Mean activation for each ImageNums (index) at each timepoint (columns).
    """
    activations = np.zeros((len(imageNums), dataView.shape[1]))
    for i_imageNum, imageNum in enumerate(imageNums):
        if disambigOnly:
            trialsMask = behav[condition] & (behav.ImageNum == imageNum) & behav["disambig"]
        else:
            trialsMask = behav[condition] & (behav.ImageNum == imageNum)
        trialNums = behav.index[trialsMask]
        imageTrials = dataView.loc[trialNums, :]
        if ((condition == "grey") and (imageTrials.shape[0] != 2)) or (
            (condition != "grey") and (imageTrials.shape[0] < 3)
        ):
            logger.info(
                f"Summing {dataView.loc[trialNums, :].shape[0]} {condition} activations for"
                f" ImageNum {imageNum}"
            )
        # This is pd.DataFrame.mean -- all NaN input => NaN, otherwise drops NaNs from calc
        activations[i_imageNum, :] = dataView.loc[trialNums, :].mean()
    return pd.DataFrame(data=activations, index=imageNums, columns=dataView.columns)


def get_activations(
    dataView: pd.DataFrame,
    behav: pd.DataFrame,
    imageNums: Sequence[int],
    condition: str,
    disambigOnly=True,
    activationDuration=0.65,
) -> pd.Series:
    """
    Get mean activations for each image in imageNums during given condition as the mean from 0s
    to activeDuration.

    Args:
        dataView (pd.DataFrame): Activation DataFrame. Index expected to be TrialNums. Columns
            expected to be timepoints.
        behav (pd.DataFrame): Behavioral DataFrame. Index expected to be TrialNums. Needs
            columns "ImageNum", "disambig", and one that matches passed condition.
        imageNums (Sequence[int]): Only these images will be included in output.
        condition (str): Only trials for this condition will be included in the output.
        disambigOnly (bool, optional): Only disambiguated trials will be included in the
            output. Defaults to True.
        activationDuration (float, optional): Activation will be calculated for each
            image as the mean from 0s to this parameter. Defaults to 0.65.

    Returns:
        pd.Series[index=ImageNum]: Activation for each image.
    """
    activations = np.zeros(len(imageNums))
    for i_imageNum, imageNum in enumerate(imageNums):
        if disambigOnly:
            trialsMask = behav[condition] & (behav.ImageNum == imageNum) & behav["disambig"]
        else:
            trialsMask = behav[condition] & (behav.ImageNum == imageNum)
        trialNums = behav.index[trialsMask]
        imageTrials = dataView.loc[trialNums, 0:activationDuration]
        if ((condition == "grey") and (imageTrials.shape[0] != 2)) or (
            (condition != "grey") and (imageTrials.shape[0] < 3)
        ):
            logger.info(
                f"Summing {imageTrials.shape[0]} {condition} activations for ImageNum {imageNum}"
            )
        # Take the mean of every non-NaN trial timepoint within the TOI
        activations[i_imageNum] = np.nanmean(imageTrials)
    return pd.Series(data=activations, index=imageNums)


def get_best_fits(A: np.ndarray) -> dict[Literal["slope", "intercept"], np.ndarray]:
    # Return the slopes and intercepts for each column vector's best fit line
    # NaN entries are excluded.
    # Slope for a column vector with <2 non-NaN entries is also NaN.
    # Returned dict with keys=['slope','intercept']
    _, numTimepoints = A.shape
    results = np.full((2, numTimepoints), np.nan)
    for i_time in range(numTimepoints):
        elementMask = ~np.isnan(A[:, i_time])
        if np.count_nonzero(elementMask) < 2:
            continue
        maskedY = A[elementMask, i_time]
        maskedX = np.linspace(0, 1, len(maskedY), endpoint=True)
        bounds = min(maskedX.min(), maskedY.min()), max(maskedX.max(), maskedY.max())
        polyResult = polyfit(maskedX, maskedY, 1, window=bounds, domain=bounds)
        results[0, i_time] = polyResult.coef[-1]
        results[1, i_time] = polyResult.coef[0]
    return {"slope": results[0, :], "intercept": results[1, :]}


def electrode_tuning(
    elecData: pd.DataFrame,
    behav: pd.DataFrame,
    imageNums: list[int],
    activationTOI: Optional[float] = None,
) -> pd.DataFrame:
    """
    Calc best fit lines of pre and post condition from grey ordering of imageNums for given
    electrode.
    Returns DF with index=MultiIndex.from_product([['Pre','Post'],['slope','intercept']]),
    columns=elecData.columns
    activationTOI: Time in seconds after image onset to mean for grey activation calculation.
        All images must have at least one trial that lasts this long, unless 0 or None is provided
        in which case each timepoint will calculate and use its own grey activation.
    """
    # Since grey accuracy is almost always near ceiling, always use all grey trials
    if activationTOI:
        greyActivations = pd.DataFrame(
            get_activations(elecData, behav, imageNums, "grey", activationDuration=activationTOI)
        )
        # Unify greyActivations to the same shape as the _at_each_timepoint version
        greyActivations = pd.DataFrame(
            np.broadcast_to(greyActivations.values, (len(greyActivations), elecData.shape[1])),
            columns=elecData.columns,
            index=greyActivations.index,
        )
    else:
        greyActivations = get_activations_at_each_timepoint(elecData, behav, imageNums, "grey")
    if greyActivations.isna().all(axis=1).any():
        logger.warning(
            f"No grey activations available for Patient_{behav.Subject.unique()[0]} images"
            f" {greyActivations.index[greyActivations.isna().all(axis=1)]}. Omitting from"
            " analysis."
        )
        greyActivations = greyActivations.loc[~greyActivations.isna().all(axis=1)]
        imageNums = greyActivations.index.tolist()

    # Get pre & post activations. Remove values where there is no grey activation
    preActivations = get_activations_at_each_timepoint(elecData, behav, imageNums, "pre")
    preActivations = preActivations.mask(greyActivations.isnull())
    postActivations = get_activations_at_each_timepoint(elecData, behav, imageNums, "post")
    postActivations = postActivations.mask(greyActivations.isnull())

    # Sort by grey activations, get best fits for other conditions
    greySortIdx = np.argsort(greyActivations.values, axis=0)
    sortedPre = np.take_along_axis(preActivations.values, greySortIdx, axis=0)
    preFits = pd.DataFrame(get_best_fits(sortedPre))
    preFits.index = elecData.columns
    logger.debug(f"Pre fits: {preFits}")
    logger.info(f"Found {preFits['slope'].count()} pre slopes.")
    sortedPost = np.take_along_axis(postActivations.values, greySortIdx, axis=0)
    postFits = pd.DataFrame(get_best_fits(sortedPost))
    postFits.index = elecData.columns
    logger.debug(f"Post fits: {postFits}")
    logger.info(f"Found {postFits['slope'].count()} post slopes.")

    bestFits = pd.concat([preFits, postFits], axis=1, keys=["Pre", "Post"]).T
    bestFits.rename_axis(index=["Condition", "Best Fit Component"], inplace=True)
    return bestFits


def make_detail_str(subjectID: str, **kwargs) -> str:
    """
    Create a string that can be used to uniquely identify a set of image preference fits.

    Args:
        subjectID (str): Patient ID number.
        kwargs: Additional image preference parameters to include in the string.

    Returns:
        str: Detail string.
    """
    detail = f"{subjectID}"
    if "smoothWindow" in kwargs:
        detail = f"{kwargs.pop('smoothWindow')}s-" + detail
    detail = "-".join([f"{val}" for val in kwargs.values()] + [detail])
    return detail


def load_subject_imagePref_fits(
    subjectID: str, imagePrefDataPath=Path("../Analysis/ImagePref/data"), **kwargs
) -> pd.DataFrame:
    """
    Load image preference fits from file if they exist and are up to date, otherwise generate.

    Args:
        subjectID (str): Patient ID number.
        smoothWindow (float): Image preference parameter.
        dataFileName (str): Image preference parameter

    Returns:
        pd.DataFrame: Best fits dataframe with MultiIndex with levels ["Electrode", "Condition",
            "Best Fit Component"]. Columns are timepoints. Loads from disk if it exists and is up
            to date.
    """
    detailsStr = make_detail_str(subjectID, **kwargs)
    ptCSVPath = imagePrefDataPath / f"imagePref_Fits_{detailsStr}.csv"
    if ptCSVPath.is_file():
        patientUpdated = patient_last_updated(int(subjectID))
        patientUpdated = datetime.min if patientUpdated is None else patientUpdated
        fitsUpdated = datetime.fromtimestamp(ptCSVPath.stat().st_mtime)
        if patientUpdated < fitsUpdated:
            logger.info(f"Loading existing image pref fits: {ptCSVPath}.")
            fits = pd.read_csv(
                ptCSVPath,
                index_col=["Electrode", "Condition", "Best Fit Component"],
            )
            fits.columns = fits.columns.astype("float")
            return fits
        else:
            logger.debug(
                f"Found outdated image pref fits for {subjectID} last updated {fitsUpdated}."
                f" Patient last updated {patientUpdated}."
            )
    fits = get_subject_image_pref_fits(subjectID, **kwargs)
    Path(ptCSVPath.parent).mkdir(parents=True, exist_ok=True)
    fits.to_csv(ptCSVPath)
    return fits


def get_subject_image_pref_fits(
    patientID: str,
    smoothWindow: float,
    disambigOnly=True,
    dataFileName="HGP",
    **kwargs,
) -> pd.DataFrame:
    data = load_mne_subject(patientID, loadGamma=dataFileName == "HGP")
    # Select images
    if disambigOnly:
        imageNums = data.behav.loc[data.behav.disambig, "ImageNum"].unique().tolist()
    else:
        imageNums = data.behav.ImageNum.unique().tolist()
    # Patient_18 errata
    if str(patientID) == "18":
        imageNums = [img for img in imageNums if img != 27]
    results = tune_electrodes(data, imageNums, smoothWindow, **kwargs)
    return results


def parallelize_get_subject_image_prefs(
    patientID: str, smoothWindow: float, **kwargs
) -> pd.DataFrame:
    logger.info(f"Requesting image pref results for patientID {patientID}.")
    results = load_subject_imagePref_fits(patientID, smoothWindow=smoothWindow, **kwargs)
    logger.info(f"Subject {patientID} results.shape = {results.shape}.")
    return results


def get_subjects_image_prefs(subjectIDs: list[int], smoothWindow: float, **kwargs) -> pd.DataFrame:
    # Outputs the best fit lines for each subject's disambiguated images in the pre and post
    # condition as each electrode is tuned by the grey condition.
    # smoothWindow is the width of the mean smoothing window to apply to the data in SECONDS.
    # Output df has a MultiIndex with levels=['SubjectID', 'Electrode', 'Condition', 'Best Fit
    # Component'], columns=timeseries
    fixed_para_get_subject_image_prefs = partial(
        parallelize_get_subject_image_prefs, smoothWindow=smoothWindow, **kwargs
    )
    subjectResults = Parallel(n_jobs=min(8, len(subjectIDs)))(
        delayed(fixed_para_get_subject_image_prefs)(subjectID) for subjectID in subjectIDs
    )
    return pd.concat(
        subjectResults,
        axis=0,
        keys=subjectIDs,
        names=["SubjectID"] + subjectResults[0].index.names,
    )


def test_subject_ROI_image_prefs(
    rois: pwmvROIs, data: pd.DataFrame, **kwargs
) -> tuple[pd.DataFrame, dict[str, dict[str, ClusterCorrectedMultipleDof]]]:
    """
    Run wilcoxon test on data for each ROI in rois.

    Looks up each patient's electrode's coordinates and determines ROI membership in rois.

    Args:
        rois (pwmvROIs): ROIs to test.
        data (pd.DataFrame): Expected to match the output from get_subjects_image_prefs:
            MultiIndex with levels=['SubjectID', 'Electrode', 'Condition', 'Best Fit
            Component'], columns=timeseries

    Returns:
        pd.DataFrame: Stats per timeseries. MultiIndex.levels=['ROI', 'Condition', 'Stats']. See
            image_pref_wilcoxon_perm_test for details on Condition & Stats. Columns are the
            timepoints.
        dict[str, dict[str, ClusterCorrectedMultipleDof]]: Dict mapping ROI to Condition to
            ClusterCorrectedMultipleDof object used to produce corresponding stats.
    """
    roiNames = sorted(rois.keys())
    roiResults = []
    memberships = rois.get_patient_memberships(
        list(data.index.get_level_values("SubjectID").unique())
    )
    coveredROIs = []
    allClusterStats: dict[str, dict[str, ClusterCorrectedMultipleDof]] = {}
    with tqdm(roiNames, unit="ROI") as roiBar:
        for roiName in roiBar:
            roiBar.set_description(roiName)
            roiElecs = memberships[memberships[roiName]].index
            numElectrodes = len(roiElecs)
            if numElectrodes:
                roiData = data[
                    data.index.droplevel(["Condition", "Best Fit Component"]).isin(roiElecs)
                ].xs("slope", level="Best Fit Component")
                logger.info(f"Analyzing ROI {roiName} consisting of {numElectrodes} electrodes.")
                coveredROIs.append(roiName)
                results, clusterStats = image_pref_wilcoxon_perm_test(roiData, **kwargs)
                roiResults.append(results)
                allClusterStats[roiName] = clusterStats
            else:
                logger.warning(f"No eletrodes found for ROI {roiName}. Skipping.")
    return (
        pd.concat(roiResults, axis=0, keys=coveredROIs, names=["ROI"] + roiResults[0].index.names),
        allClusterStats,
    )


def imagePrefs_survival_plot(
    preChannels: pd.Series,
    postChannels: pd.Series,
    ax: plt.Axes,
    roiName: str,
    toi: slice = None,
    ls=["-", "-"],
    lw=[1, 1],
):
    """
    Plot channel survival for Pre and Post conditions.

    Args:
        preChannels (pd.Series): Surviving channels at each timepoint in Pre.
        postChannels (pd.Series): Surviving channels at each timepoint in Post.
        ax (plt.Axes): Axes to plot on.
        roiName (str): ROI name to plot in title.
        toi (slice, optional): Timepoints to plot. Defaults to None, which plots all.
        ls (list, optional): Linestyles as [Pre, Post]. Defaults to ["-", "-"].
        lw (list, optional): Linewidths as [Pre, Post]. Defaults to [1, 1].
    """
    if not toi:
        toi = slice(None)
    ax.plot(preChannels.index, preChannels, label="Pre", ls=ls[0], lw=lw[0])
    ax.plot(postChannels.index, postChannels, label="Post", ls=ls[1], lw=lw[1])
    ax.legend(loc="lower left")
    ax.set_xlabel("seconds")
    ax.set_ylabel("surviving channels")
    ax.set_title(f"{roiName} Channel Survival", fontsize="large")
    ax.grid()


def imagePrefs_figure(
    roiData: pd.DataFrame,
    ax: plt.Axes,
    roiName: str,
    yUnit: str,
    toi: slice = None,
    colors: Mapping[str, tuple[float, float, float] | str] = {
        "Pre": (0, 0, 1),
        "Post": (1, 0, 0),
        "Post - Pre": "k",
    },
    alpha=0.2,
    semToUse: str | None = "Post - Pre",
):
    """
    Print Image Preference Analysis plot for one ROI to an axes.

    Args:
        roiData (pd.DataFrame): MultiIndex.levels=['Condition','Stats'], where Stats includes
            'slope mean', 'slope SEM', 'significant', and 'numElec'. Columns are the
            timepoints to plot.
        ax (plt.Axes): Axes on which to plot.
        roiName (str): Used for the title.
        yUnit (str): Y-axis units label.
        toi (slice, optional): Columns will be sliced to this, if given. Defaults to None.
        colors (dict, optional): Colors to use for each condition. Defaults to {"Pre": (0, 0, 1),
            "Post": (1, 0, 0), "Post - Pre": "k"}.
        alpha (float, optional): Alpha for the SEM fill_between bands. Defaults to 0.2.
        semToUse (str | None, optional): Which condition's SEM to use for the fill_between bands.
            Set to None for each condition to use its own SEM. Defaults to "Post - Pre".
    """
    if toi:
        roiData = roiData.loc[:, toi]
    columns = roiData.columns

    ax.axvline(0, c="k", linestyle=":")
    ax.axhline(0, c="k", linestyle=":")

    means = roiData.loc[("Pre", "slope mean"), :]
    sems = roiData.loc[("Pre" if semToUse is None else semToUse, "slope SEM"), :]
    ax.fill_between(
        list(columns), list(means - sems), list(means + sems), color=colors["Pre"], alpha=alpha
    )
    ax.plot(columns, means, color=colors["Pre"], label="Pre", linewidth=2.0)

    means = roiData.loc[("Post", "slope mean"), :]
    sems = roiData.loc[("Post" if semToUse is None else semToUse, "slope SEM"), :]
    ax.fill_between(
        list(columns), list(means - sems), list(means + sems), color=colors["Post"], alpha=alpha
    )
    ax.plot(columns, means, color=colors["Post"], label="Post", linewidth=2.0)

    plot_sig_bar(
        ax, roiData.loc[("Post - Pre", "significant"), :].values, colors["Post - Pre"], columns
    )  # , label="Sig Post-Pre")
    plot_sig_bar(ax, roiData.loc[("Pre", "significant"), :].values, colors["Pre"], columns)
    plot_sig_bar(ax, roiData.loc[("Post", "significant"), :].values, colors["Post"], columns)

    if roiData.shape[1]:
        numElecs = int(roiData.loc[("Pre", "numElec"), 0])
        assert (
            roiData.xs("numElec", level="Stats")[0] == numElecs
        ).all(), f"Number of electrodes varies across conditions for {roiName}"
    else:
        numElecs = 0

    ax.set_title(f"{roiName} ({numElecs} elecs)", fontsize="large")
    ax.set_ylabel(yUnit)
    ax.set_xlabel("seconds")
    ax.grid(axis="x")
    # return fig,ax


def imagePrefs_figures(
    data: Union[str, pd.DataFrame], toi: slice = None, yUnit="log[µV²]", source="HGP"
) -> plt.Figure:
    """
    Call imagePrefs_figure on all ROIs in data. Produces one figure with a subplot for each ROI.

    Args:
        data (Union[str, pd.DataFrame]): Requires a MultiIndex with the top level being ROI.
        toi (slice, optional): See imagePrefs_figure. Defaults to None.
        yUnit (str, optional): imagePrefs_figure. Defaults to "log[µV²]".
        source (str, optional): Used for whole figure suptitle. Defaults to "HGP".

    Returns:
        plt.Figure: _description_
    """
    if isinstance(data, str):
        data = pd.read_csv(data, index_col=[0, 1, 2])
        data.columns = data.columns.astype("float")
    if toi:
        data = data.loc[:, toi]

    alpha = 0.2
    numROI = data.index.get_level_values("ROI").nunique()
    plotCols = int(np.ceil(np.sqrt(numROI)))
    plotRows = int(np.ceil(numROI / plotCols))
    fig, axarr = plt.subplots(
        plotRows, plotCols, figsize=(4 * plotCols, 3 * plotRows), sharex=True
    )
    for i_roi, roiName in enumerate(data.index.get_level_values("ROI").unique()):
        ax = axarr.flatten()[i_roi]
        roiData = data.loc[roiName]
        imagePrefs_figure(
            roiData, ax, roiName, yUnit if i_roi % plotCols == 0 else "", toi=toi, alpha=alpha
        )

    # axarr[-1, -1].set_visible(False)
    if i_roi < plotRows * plotCols:
        # Put a blank axes in last subplot so .legend() there will populate correctly
        imagePrefs_figure(pd.DataFrame(index=roiData.index), axarr[-1, -1], "", "")
    axarr[-1, -1].legend(loc="lower left")
    fig.suptitle(f"Image Preference in Pre, Post Conditions, {source}", fontsize=22)
    return fig


def generate_and_plot(
    subjects: Subjects,
    smoothWindowSecs: float,
    toi: slice,
    dataFileName: Literal["HGP", "ERP"] = "HGP",
    imagePrefPath="../Analysis/ImagePref",
    roiGroupName="network",
    numRepeats=1000,
    loadStatsFromDisk=False,
    greyActivationWindow=None,
    plotLegend=True,
    roiFilenamePrefixes: dict[str, str] = None,
):
    """
    Generate image preference fits for subjects and all ROIs in roiGroupName and plot them.

    Args:
        subjects (Subjects): Subject numbers.
        smoothWindowSecs (float): Width in seconds to center smooth data.
        toi (slice): Time Of Interest on which to run stats. This does not affect how the fits are
            calculated.
        dataFileName (Literal["HGP", "ERP"], optional): Whether to load the data
            as the high gamma power transform or the raw ERP. Defaults to "HGP".
        imagePrefPath (str, optional): Path to analysis resutls folder. Defaults to
            "../Analysis/ImagePref".
        roiGroupName (str, optional): ROI group name. Must refer to a folder in
            Electrodes/ROI_Masks/. Defaults to "network".
        numRepeats (int, optional): Target null distribution size for permutation testing. Defaults
            to 1000.
        loadStatsFromDisk (bool, optional): Whether to try and load previously calculated stats
            data from disk. Defaults to False, which will always generate the stats from the fits.
        plotLegend (bool, optional): Whether to plot legends on the stats figures. Defaults to
            True.
        roiFilenamePrefixes (dict[str,str], optional): Prefixes to add to saved figure filenames,
            keyed by ROI. Can be used to induce an alternate alphabetical ordering of filenames
            produced. Defaults to no prefixes.
    """
    if roiFilenamePrefixes is None:
        roiFilenamePrefixes = {}
    savePath = Path(imagePrefPath, f"{roiGroupName}_ROIs")
    savePath.mkdir(parents=True, exist_ok=True)
    rois = pwmvROIs(f"Electrodes/ROI_Masks/{roiGroupName}")
    detailsStr = make_detail_str(
        str(subjects),
        dataFileName=dataFileName,
        smoothWindow=smoothWindowSecs,
        roiGroupName=roiGroupName,
        greyActivationWindow=greyActivationWindow,
    )
    print(f"Processing {detailsStr}")
    statsPath = Path(savePath, "CSV")
    statsPath.mkdir(parents=True, exist_ok=True)
    csvPath = statsPath / f"imagePrefs_Stats_{detailsStr}.csv"
    if loadStatsFromDisk and csvPath.is_file():
        print(f"LOADING FROM DISK:{csvPath}")
        roiStats = pd.read_table(csvPath, sep=",", index_col=["ROI", "Condition", "Stats"])
        roiStats.columns = roiStats.columns.astype("float")
    else:
        fits = get_subjects_image_prefs(
            subjects,
            smoothWindowSecs,
            dataFileName=dataFileName,
            greyActivationWindow=greyActivationWindow,
        )
        if not fits.columns.is_monotonic_increasing:
            fits = fits.sort_index(axis=1)
        trimmedFits = fits.loc[:, (fits.columns >= toi.start) & (fits.columns <= toi.stop)]
        # Window assumes 512Hz sample rate
        window = int(np.floor(smoothWindowSecs * 512)) if dataFileName == "HGP" else None
        zeroIdx = trimmedFits.columns.to_list().index(0)
        testSlice = slice(None if window is None else zeroIdx + (window - 1) // 2, None)
        roiStats, clusterStats = test_subject_ROI_image_prefs(
            rois, trimmedFits, numRepeats=numRepeats, testWindow=testSlice
        )
        roiStats.to_csv(csvPath)
        print(f"Saving {csvPath}.")
        with open(csvPath.with_suffix(".pkl"), "wb") as f:
            pickle.dump(clusterStats, f)
        print(f"Saving stats as {csvPath.with_suffix('.pkl')}.")
    if dataFileName == "HGP":
        titleDetails = f"HGP N={len(subjects)} Window={smoothWindowSecs*1000:.0f}ms"
        yUnit = "log[µV²]"
    else:
        titleDetails = f"ERP N={len(subjects)} Window={smoothWindowSecs*1000:.0f}ms"
        yUnit = "µV"

    # Print group plot
    fig = imagePrefs_figures(roiStats, yUnit=yUnit, source=titleDetails, toi=toi)
    if not plotLegend:
        fig.axes[-1, -1].get_legend().remove()
    fig.tight_layout(pad=0.5)
    fname = Path(savePath, f"imagePrefs_cluster-corrected_{detailsStr}.png")
    print(f"Saving {fname}.")
    fig.savefig(fname)
    plt.close(fig)

    # Print individual plots
    plt.rcParams.update({
        "font.size": 18,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "ytick.right": False,
    })
    fullNameSavePath = Path(savePath, "fullROINames")
    fullNameSavePath.mkdir(parents=True, exist_ok=True)
    for roiName in roiStats.index.get_level_values("ROI").unique():
        fig, ax = plt.subplots(figsize=(12, 8))
        roiData = roiStats.loc[roiName]
        imagePrefs_figure(roiData, ax, rois.fullNames[roiName], yUnit, toi=toi)
        if plotLegend:
            ax.legend(loc="lower left", labelspacing=0.3)
        else:
            ax.get_legend().remove()
        fig.tight_layout(pad=0.5)
        fname = Path(
            fullNameSavePath,
            f"{roiFilenamePrefixes.get(roiName,'')}imagePrefs_cc_{detailsStr}_{roiName}.png",
        )
        print(f"Saving {fname}.")
        fig.savefig(fname)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(9, 6))
        imagePrefs_survival_plot(
            roiData.loc[("Pre", "numElec")],
            roiData.loc[("Post", "numElec")],
            ax,
            roiName,
            toi=toi,
        )
        fname = Path(
            fullNameSavePath,
            f"{roiFilenamePrefixes.get(roiName,'')}imagePrefsSurvival_cc_{detailsStr}_{roiName}.png",
        )
        print(f"Saving {fname}.")
        fig.savefig(fname)
        plt.close(fig)


class ImagePreferenceAnalysis(PWMVAnalysis):
    def __init__(
        self,
        subjects: Subjects,
        toi=slice(-0.2, 1.5),
        source: Literal["HGP", "ERP"] = "HGP",
        smoothWindowSecs=0.1,
        roiGroupName: str = "network",
        numRepeats=1000,
        greyActivationWindow=None,
        outputPath: Path = Path("..", "Analysis", "ImagePref"),
        rngSeed=20220506,
    ) -> None:
        super().__init__(subjects, outputPath, smoothWindowSecs=smoothWindowSecs)
        self.toi = toi
        self.source = source
        self.roiGroupName = roiGroupName
        self.numRepeats = numRepeats
        self.greyActivationWindow = greyActivationWindow
        self.saveFolder.mkdir(parents=True, exist_ok=True)
        self.rng = np.random.default_rng(rngSeed)

        self.rois = pwmvROIs(Path("Electrodes", "ROI_Masks", roiGroupName))

    def __sizeof__(self) -> int:
        sizeSum = (
            getsizeof(self.toi)
            + getsizeof(self.source)
            + getsizeof(self.roiGroupName)
            + getsizeof(self.numRepeats)
            + getsizeof(self.greyActivationWindow)
            + getsizeof(self.rng)
            + getsizeof(self.rois)
        )
        if hasattr(self, "testSlice"):
            sizeSum += getsizeof(self.testSlice)
        if hasattr(self, "clusterStats"):
            sizeSum += getsizeof(self.clusterStats)
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
        if hasattr(self, "fits"):
            delattr(self, "fits")
        super().save(path)

    def run(self, loadFromDisk=True, rng: np.random.Generator = None):
        """
        Loads or generates image preference fits and stats.

        Sets .results and .clusterStats attributes.

        Args:
            loadFromDisk (bool, optional): If False, will never try to load saved results from
                disk. Does not affect whether generated fits load from disk. Defaults to True.
            rng (np.random.Generator, optional): Random number generator to use for permutation
                tests. Defaults to a generator using this object's rngSeed.
        """
        self.load_results(loadFromDisk, rng)
        print("Results loaded.")

    def load_results(self, loadFromDisk=True, rng: np.random.Generator = None):
        """
        Sets .results and .clusterStats attributes from disk or by generating them.

        Args:
            loadFromDisk (bool, optional): If False, will never try to load saved results from
                disk. Does not affect whether generated fits load from disk. Defaults to True.
            rng (np.random.Generator, optional): Random number generator to use for permutation
                tests. Defaults to a generator using this object's rngSeed.
        """
        resultsFile = self.dataFolder / f"imagePrefs_Stats_{self.make_detail_str()}.csv"
        # Assumes the cluster stats are always saved with the same name as the stat results
        if (
            loadFromDisk
            and resultsFile.is_file()
            and (datetime.fromtimestamp(resultsFile.stat().st_mtime) > self._update_before_date)
        ):
            print(f"Loading from disk: {resultsFile}.")
            self.results = pd.read_table(
                resultsFile, sep=",", index_col=["ROI", "Condition", "Stats"]
            )
            self.results.columns = self.results.columns.astype("float")
            with open(resultsFile.with_suffix(".pkl"), "rb") as f:
                self.clusterStats = pickle.load(f)
        else:
            print(f"Generating results for {self.make_detail_str()}.")
            self.results, self.clusterStats = self._load_results(rng)
            print(f"Saving {resultsFile}.")
            self.results.to_csv(resultsFile)
            with open(resultsFile.with_suffix(".pkl"), "wb") as f:
                pickle.dump(self.clusterStats, f)

    def _load_results(
        self, rng: np.random.Generator = None
    ) -> tuple[pd.DataFrame, dict[str, dict[str, ClusterCorrectedMultipleDof]]]:
        if rng is None:
            rng = self.rng
        if not hasattr(self, "fits"):
            self.load_data()
        assert self.fits.columns.is_monotonic_increasing
        trimmedFits = self.fits.loc[
            :, (self.fits.columns >= self.toi.start) & (self.fits.columns <= self.toi.stop)
        ]
        if (self.source == "HGP") and (self._smoothingWindow is not None):
            zeroIdx = trimmedFits.columns.to_list().index(0)
            self.testSlice: slice = slice(zeroIdx + (self._smoothingWindow - 1) // 2, None)
        else:
            self.testSlice = slice(None)
        results = test_subject_ROI_image_prefs(
            self.rois,
            trimmedFits,
            numRepeats=self.numRepeats,
            testWindow=self.testSlice,
            rngSeed=rng,
            maxWorkers=self.maxWorkers,
        )
        # Unload Niftis to save memory
        self.rois = pwmvROIs(Path("Electrodes", "ROI_Masks", self.roiGroupName))
        return results

    def load_data(self):
        print("Getting fits.")
        self.fits = get_subjects_image_prefs(
            self.subjects,
            self.smoothWindowSecs,
            dataFileName=self.source,
            greyActivationWindow=self.greyActivationWindow,
            imagePrefDataPath=self.dataFolder,
        )
        self.data = self.fits

    def make_detail_str(self) -> str:
        return (
            super()._basic_detail_str(
                source=self.source,
                roiGroup=self.roiGroupName,
                greyWindow=self.greyActivationWindow,
                smoothing=self.smoothWindowSecs,
            )
            + f"_{self.subjects}"
        )


def bootstrap_results(
    analysis: ImagePreferenceAnalysis,
    saveDir: FilePath,
    rngSeed: int | None = None,
    roiTests: dict[str, list[str]] | None = None,
    numBootstraps=1000,
    saveTime=1200,
    maxWorkers=MAXWORKERS,
) -> None:
    """
    Bootstrap the results of an ImagePreferenceAnalysis and saves the ruselts to disk.

    Results are saved as a pickled list of ClusterCorrectedMultipleDof objects.

    Args:
        analysis (ImagePreferenceAnalysis): Analysis object to bootstrap.
        saveDir (FilePath): Path to save the bootstrapped results.
        rngSeed (int | None, optional): Seed for the random number generator. Defaults to a random
            seed.
        roiTests (dict[str, list[str]] | None, optional): Dict of ROIs and which tests to
            bootstrap. Keys are ROI names, values are a list of tests to bootstrap in that ROI.
            Defaults to all ROI in analysis.rois and ["Post>Pre", "Pre>0", "Post>0"] for each.
        numBootstraps (int, optional): Number of bootstraps to run. Defaults to 1000.
        saveTime (int, optional): Bootstrapped results generated so far are saved to disk if at
            least this many seconds have passed since the last save to disk. Defaults to 1200.
        maxWorkers (_type_, optional): Maximum number of parallel processes to generate. Defaults
            to MAXWORKERS.

    Raises:
        ValueError: For very patholigical data that fails to generate points repeatedly.
    """
    rng = np.random.default_rng(rngSeed)
    Path(saveDir).mkdir(parents=True, exist_ok=True)
    if roiTests is None:
        roiTests = {roi: ["Post>Pre", "Pre>0", "Post>0"] for roi in sorted(analysis.rois.keys())}
    startTime = time()
    analysis.load_data()
    trimmedFits = analysis.fits.loc[
        :,
        (analysis.fits.columns >= analysis.toi.start)
        & (analysis.fits.columns <= analysis.toi.stop),
    ]
    if (analysis.source == "HGP") and (analysis._smoothingWindow is not None):
        zeroIdx = trimmedFits.columns.to_list().index(0)
        testSlice = slice(zeroIdx + (analysis._smoothingWindow - 1) // 2, None)
    else:
        testSlice = slice(None)
    memberships = analysis.rois.get_patient_memberships(
        list(analysis.data.index.get_level_values("SubjectID").unique())
    )
    statTest = partial(always_signedT_wilcoxon, test="greater", zero_method="zsplit")
    print(f"Generating bootstrap results for {roiTests}.")
    for roiName, tests in roiTests.items():
        # Select ROI data
        roiElecs = memberships[memberships[roiName]].index
        if not roiElecs.size:
            raise ValueError(f"No electrodes found for ROI {roiName}.")
        roiData = trimmedFits[
            trimmedFits.index.droplevel(["Condition", "Best Fit Component"]).isin(roiElecs)
        ].xs("slope", level="Best Fit Component")
        for test in tests:
            savePath = Path(
                saveDir,
                f"{roiName}_image_preference_{analysis.source}_{analysis.subjects}_"
                f"{test.replace(' ','_')}_bootstrap_results.pkl",
            )
            bootResults = []
            for _ in tqdm(
                range(numBootstraps),
                desc="-".join([roiName, test]),
                unit="replication",
                smoothing=0,
            ):
                # Guard against very problematic data failing to generate points repeatedly
                maxRetries = 10
                for i_retries in range(maxRetries):
                    bootstrapPre = roiData.xs("Pre", level="Condition").sample(
                        frac=1, replace=True
                    )
                    bootstrapPost = roiData.xs("Post", level="Condition").loc[bootstrapPre.index]
                    if test == "Post>Pre":
                        condA = bootstrapPost
                        condB = bootstrapPre
                    else:
                        condA = bootstrapPre if test == "Pre>0" else bootstrapPost
                        condB = np.zeros_like(condA)
                    try:
                        bootResults.append(
                            ClusterCorrectedMultipleDof(
                                condA,
                                condB,
                                statTest,
                                "paired greater" if (test == "Post>Pre") else "greater",
                                test,
                                rng,
                                numRepeats=analysis.numRepeats,
                                cdt=0.05,
                                testWindow=testSlice,
                                maxWorkers=maxWorkers,
                            )
                        )
                        break
                    except InsufficientPointGenerationError as e:
                        if i_retries == maxRetries - 1:
                            e.add_note(
                                f"Failed to generate points {maxRetries} censecutive times for"
                                f" {roiName} {test}."
                            )
                            raise
                if time() - startTime > saveTime:
                    with open(savePath, "wb") as f:
                        pickle.dump(bootResults, f)
            with open(savePath, "wb") as f:
                pickle.dump(bootResults, f)
