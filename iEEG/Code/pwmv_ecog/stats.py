import logging
from functools import partial
from sys import getsizeof
from typing import Any, Callable, Generator, Literal, Sequence, TypeAlias


import numpy as np
import pandas as pd
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

from pwmv_ecog import MAXWORKERS

logger = logging.getLogger(__name__)
Array: TypeAlias = np.ndarray | pd.DataFrame
Vector: TypeAlias = np.ndarray | pd.Series
Endpoints: TypeAlias = list[tuple[int, int]]


def get_cluster_SSs(
    statValues: Sequence[float], pvalues: Sequence[float], cdt=0.05, twoHanded=False
) -> tuple[list[float], Endpoints]:
    """
    Calc cluster summary statistics and boundaries of the clusters. Assumes 1d time series data.

    Clusters are defined as adjacent significant (pvalue < cdt) statistics with matching sign.

    Args:
        statValues (Sequence[float]): Test statistics. Should be same length as pvalues.
        pvalues (Sequence[float]): P-value for each test statistic. Should be same length as
            statValues.
        cdt (float, optional): Cluster Determined Threshold. Defaults to 0.05.
        twoHanded (bool, optional): Whether to treat as two-handed stats: the value of cdt will be
            halved.

    Returns:
        Tuple[List[float], Endpoints]: Returns a list of cluster summary stats, and the endpoints
            for each cluster. First endpoint inclusive, second exclusive. If no points are
            significant, the first list will contain a single element of zero, and second list will
            be empty.
    """
    if twoHanded:
        cdt /= 2
    results: list[float] = []
    clusterEndpoints = []
    inCluster = False
    start = 0
    for i_val, (stat, p) in enumerate(zip(statValues, pvalues)):
        if inCluster and ((p > cdt) or (np.sign(stat) != np.sign(results[-1]))):
            clusterEndpoints.append((start, i_val))
            inCluster = False
        if p <= cdt:
            if inCluster:
                results[-1] += stat
            else:
                results.append(stat)
                inCluster = True
                start = i_val
    if inCluster:
        clusterEndpoints.append((start, i_val + 1))
    if len(results) == 0:
        results = [0]
    return results, clusterEndpoints


def paired_permute_and_cluster(
    A: np.ndarray,
    B: np.ndarray,
    statisticCalculator: Callable[[Array, Array], tuple[Vector, Vector]],
    twoSided: bool,
    rng: np.random.Generator,
    numRepeats=1000,
    cdt=0.05,
    maxWorkers: int = MAXWORKERS,
    **tqdmKWargs,
) -> tuple[list[list[float]], list[Endpoints]]:
    """
    Generate cluster stats for random permutations of the rows of A and B.

    If 2^(A.shape[0]) < numRepeats, all possible permutations will be generated.

    Args:
        A (np.ndarray): First dataset. Should have same shape as B.
        B (np.ndarray): First dataset. Should have same shape as A.
        statisticCalculator (Callable[[Array, Array], tuple[Vector,Vector]]): Function to calculate
            the test statistic. Should take two arrays of the same shape and return two vectors:
            the first of test statistics, the second of p-values. Must be a symmetric test and done
            on the columns of the arrays.
        twoSided (bool): Whether the cluster stats should follow a two-sided significance test.
        rng (np.random.Generator): Initialized generator to perform the permutations.
        numRepeats (int, optional): Target number of permutations to generate. Defaults to 1000.
        cdt (float, optional): Cluster Determined Threshold sets eligibility for cluster inclusion
            of each test point. Defaults to 0.05. Halved if twoSided is True.
        maxWorkers (int, optional): Maximum number of workers to use for parallel processing.

    Returns:
        tuple[list[list[float]], list[Endpoints]]: (list of cluster summary stats, list
            of cluster endpoints). Both outer lists' lengths equal to the number of permutations
            generated which will be no greater than numRepeats. Each pair of inner
            lists are of equal length corresponding to the number of clusters found for the
            corresponding permutation; unless no points are significant for a given permutation,
            then the stats list will contain a single element of 0, and Endpoints will be
            empty. Each cluster endpoint tuple contains the start and end indices of the cluster;
            first endpoint inclusive, second exclusive.
    """
    # Make permutations
    allStats, allPVals = create_permutations(
        A,
        statisticCalculator,
        rng,
        altData=B,
        numRepeats=numRepeats,
        maxWorkers=maxWorkers,
        **tqdmKWargs,
    )
    numPerms = len(allStats)
    # Calc null distribution from permutations' cluster stats
    clusterStats = []
    clusterEndpoints = []
    for i_perm in tqdm(range(numPerms), desc="Cluster Stats", unit="permutation", leave=False):
        rawClusterStats, endpoints = get_cluster_SSs(
            allStats[i_perm], allPVals[i_perm], cdt=cdt, twoHanded=twoSided
        )
        clusterStats.append(rawClusterStats)
        clusterEndpoints.append(endpoints)
    return clusterStats, clusterEndpoints


def _test_adjust(
    data: Sequence[float], test: Literal["greater", "lesser", "two-sided", "paired greater"]
) -> np.ndarray[float]:
    """
    Adjust statistics to allow always testing for extreme positive values.
    """
    if test == "two-sided":
        adjustedData = np.abs(np.array(data))
    elif test == "lesser":
        adjustedData = -1 * np.array(data)
    elif test in ["greater", "paired greater"]:
        adjustedData = np.array(data)
    else:
        raise ValueError(f"Unknown test type: {test}")
    return adjustedData


def get_cluster_DoF(endpoints: tuple[int, int], numSamples: np.ndarray) -> float:
    """
    Return degrees of freedom for cluster with endpoints in numSamples.

    Args:
        endpoints (tuple[int, int]): Start and end index of the cluster.
        numSamples (np.ndarray): Number of samples in each time point.

    Returns:
        float: Mean number of samples in the cluster.
    """
    return numSamples[endpoints[0] : endpoints[1]].mean()


def bin_DoF(DoF: float, bins: np.ndarray[Any, np.dtype[np.signedinteger[Any]]]) -> int:
    """
    Assign DoF to one of bins.
    """
    return int(np.round(DoF))


def define_DoF_bins(
    dofsOfInterest: Sequence[float], numSamples: np.ndarray
) -> np.ndarray[Any, np.dtype[np.signedinteger[Any]]]:
    """
    Set up integer bins for degrees of freedom null distributions.

    Args:
        dofsOfInterest (Sequence[float]): DoFs of interest, i.e. of the real data. Currently not
            used.
        numSamples (np.ndarray): Number of samples (channels) at each column of real data.

    Returns:
        np.ndarray[Any, np.dtype[np.signedinteger[Any]]]: Bin names for the DoFs. DoFs will be
            rounded to the nearest integer to determine bin membership.
    """
    return np.arange(int(np.floor(numSamples.min())), int(np.ceil(numSamples.max())) + 1)


class InsufficientPointGenerationError(ValueError):
    pass


class ClusterCorrectedMultipleDof:
    hardCapPerms: int = 1_000_000

    def __init__(
        self,
        A: Array,
        B: Array,
        statisticCalculator: Callable[[Array, Array], tuple[Vector, Vector]],
        test: Literal["greater", "lesser", "two-sided", "paired greater"],
        testName: str,
        rng: np.random.Generator,
        numRepeats=1000,
        cdt=0.05,
        clusterAlpha=0.05,
        chunksize: int | None = None,
        testWindow: slice = slice(None),
        maxWorkers: int = MAXWORKERS,
        errorIfLessThan: int = 800,
    ):
        """
        Paired test of the rows of A - B at each column, with multiple comparison correction by
        cluster-based permutation testing.

        Handles clusters having varying degrees of freedom (DoF) due to sample drop out as
        indicated by NaNs in A and B. Assuming NaNs are expected to be omitted in
        statisticCalculator, the DoF of the produced test statistics will differ between columns
        with different numbers of NaNs.

        A null distribution is created for each DoF calculated from a significant cluster
        in A-B. The cluster from each permutation with the most extreme cluster statistic according
        to test is added to the null distrution matching the cluster's DoF. This is
        repeated until all required null distributions have at least numRepeats data points, or all
        possible permutations have been generated. Cluster DoF is calculated as the mean
        number of samples (non-NaN entries in a column) across all columns in the cluster.

        The permutation cap is the lesser of 2^(A.shape[0]) (all possible permutations) and 1e6
        (arbitrary hard cap).

        Args:
            A (np.ndarray): First dataset. Should have same shape as B.
            B (np.ndarray): Second dataset. Should have same shape as A.
            statisticCalculator (Callable[[Array, Array], tuple[Vector,Vector]]): Function to
                calculate the test statistic. Should take two arrays of the same shape and return
                two vectors: the first of test statistics, the second of p-values. Must be a
                symmetric test.
            test (Literal['greater', 'lesser', 'two-sided']): Type of statistical test being
                performed.
            testName (str): Display name of current test. For logging purposes only.
            rng (np.random.Generator): Initialized generator to perform the permutations.
            numRepeats (int, optional): Number of permutations to generate. Defaults to 1000.
            cdt (float, optional): Cluster Determined Threshold sets eligibility for cluster
                inclusion of each test point. Defaults to 0.05.
            clusterAlpha (float, optional): Cluster significance alpha level. A cluster with summed
                statistic greater than the [(1 - clusterAlpha) multiplied by the number of
                points in the appropriate null distirbution] percentile in the null is reported as
                significant.
            chunksize (float, optional): Starting number of permutations to generate in a batch
                until all necessary null distributions have at least numRepeats data points. The
                batch size used will automatically increase if the hunt requires more than a few
                batches. Defaults to numRepeats.
            testWindow (slice, optional): Slice of columns to test and produce stats for. Defaults
                to slice(None), all.
            maxWorkers (int, optional): Maximum number of workers to use for parallel processing.
            errorIfLessThan (int, optional): Once the rate of null distribution points generated
                per permutation has stablized, if the rate indicates less than this value will be
                produced within self.hardCapPerms permutations, raise an
                InsufficientPointGenerationError. This prevents unproductive permutations with
                problematic data. Defaults to 800, per Pernet, et al. (2015)
                [https://doi.org/10.1016/j.jneumeth.2014.08.003].

        Sets attributes:
            sigClusterPlotPoints (np.ndarray): 1d array of length A.shape[1]. Contains a 1 at
                statistically significant points, a NaN otherwise.
            clusterPVals (np.ndarray): P-value for each cluster.
            clusterEndpoints (list[tuple[int, int]]): List of endpoints for each cluster. First
                endpoint inclusive, second exclusive. In same order as clusterPVals.
            clusterDoFs (list[float]): List of degrees of freedom for each cluster. In same order
                as clusterPVals.
            clusterBins (list[int]): List of bins for each cluster, i.e. which null distribution
                was used for the cluster's significance testing. In same order as clusterPVals.
            stats (np.ndarray): Contains the test stat for each time point.
            numSamples (np.ndarray): Number of samples (channels) at each time point, i.e. for each
                column, how many rows are not NaN in either A or B.
            nullSizes (dict[int,int]): Mapping from bin label (DoF) to final size of corresponding
                null distribution.
            maxPerms (int): Cap on permutations to be generated.
            permsCreated (int): Number permutations actually generated.
            tested: (np.ndarray): 1d array of length A.shape[1]. Contains a 1 at points tested, a
                NaN otherwise. Validates testWindow.
            refIndex: (pd.Index): A.columns if A is a DataFrame, None otherwise. Convenience for
                interpreting sigClusterPlotPoints, clusterEndpoints, stats, numSamples and tested
                without A.
        """
        self.testWindow = testWindow
        self.tested = np.full(A.shape[1], np.nan)
        self.tested[testWindow] = 1
        if chunksize is None:
            chunksize = numRepeats
        self.refIndex = A.columns if isinstance(A, pd.DataFrame) else None
        if isinstance(A, pd.DataFrame):
            A = A.values
        if isinstance(B, pd.DataFrame):
            B = B.values
        # Trim data to testWindow
        A = A[:, testWindow]
        B = B[:, testWindow]
        # Calc real data's stats
        realStats, realPVals = statisticCalculator(A, B)
        self.stats = np.array(realStats)
        # Calc real clusters and cluster stats
        rawClusterStats, clusterEndPoints = get_cluster_SSs(
            realStats, realPVals, cdt=cdt, twoHanded=test == "two-sided"
        )
        self.clusterEndpoints = clusterEndPoints
        self.sigClusterPlotPoints = np.full(A.shape[1], np.nan)
        self.numSamples = (~np.isnan(A - B)).sum(axis=0)
        if max(self.numSamples) <= 1:
            # There is no test for which a single sample can be significant.
            return
        self.clusterDoFs = [get_cluster_DoF(eps, self.numSamples) for eps in clusterEndPoints]
        nullDistBins = define_DoF_bins(self.clusterDoFs, self.numSamples)
        nullDistributions: dict[int, list[float]] = {dfbin: [] for dfbin in nullDistBins}
        self.clusterBins = [bin_DoF(df, nullDistBins) for df in self.clusterDoFs]
        # Hard cap permutations
        self.maxPerms = min(self.hardCapPerms, 2 ** int(max(self.numSamples)))
        self.permsCreated = 0
        numNoClusterPerms = 0  # Track number of permutations with no clusters
        needMorePerms = len(clusterEndPoints) > 0
        with (
            tqdm(
                total=self.maxPerms,
                desc="Max Perms",
                position=1,
                unit_scale=True,
                leave=False,
                unit="perms",
                smoothing=0,
            ) as maxPermBar,
            tqdm(
                total=numRepeats,
                desc="Null Size",
                position=2,
                leave=False,
                unit="point",
                smoothing=0,
            ) as minDistBar,
        ):
            # Create permutations until all required null distributions have at least numRepeats
            # points or we've created the maximum number of permutations
            tqdmParams = {
                "position": 3,
                "leave": False,
                "desc": testName,
                "unit": "permutation",
                "smoothing": 0,
            }
            nextChunksize = chunksize
            while needMorePerms:
                # Request chunksize permutations
                permutationStats, permEndpoints = paired_permute_and_cluster(
                    A,
                    B,
                    statisticCalculator,
                    test == "two-sided",
                    rng,
                    nextChunksize,
                    cdt,
                    maxWorkers=maxWorkers,
                    **tqdmParams,
                )
                curNumPerms = len(permutationStats)
                self.permsCreated += curNumPerms
                maxPermBar.update(curNumPerms)
                # Assign permutations' cluster stats to appropriate null distributions
                for i_perm, permStat in enumerate(
                    tqdm(permutationStats, leave=False, desc="Assigning")
                ):
                    if permStat == [0]:
                        numNoClusterPerms += 1
                        continue
                    adjustedStats = _test_adjust(permStat, test)
                    permClusterIdx = np.argmax(adjustedStats)
                    permDoF = get_cluster_DoF(
                        permEndpoints[i_perm][permClusterIdx], self.numSamples
                    )
                    nullDistributions[bin_DoF(permDoF, nullDistBins)].append(
                        adjustedStats[permClusterIdx]
                    )
                curSizes = np.array([len(v) for v in nullDistributions.values()])
                # Calculate number of zeros to assign to each null distribution bin
                zeroDist = rng.choice(
                    list(nullDistributions.keys()), numNoClusterPerms, p=curSizes / sum(curSizes)
                )
                estTotalTargetSizes = np.array([
                    size + (zeroDist == bn).sum()
                    for bn, size in zip(nullDistributions.keys(), curSizes)
                    if bn in self.clusterBins
                ])
                if len(estTotalTargetSizes) > 0 and min(estTotalTargetSizes) != minDistBar.n:
                    minDistBar.update(min(estTotalTargetSizes) - minDistBar.n)

                # Calc how many more permutations are needed to reach numRepeats as a mult of perms
                # generated so far
                if min(estTotalTargetSizes) <= 10 and self.permsCreated < self.hardCapPerms * 0.05:
                    # Use fixed growth rate until enough iterations to estimate true rate
                    morePermsFactor = 2
                elif (
                    self.maxPerms * (min(estTotalTargetSizes) / self.permsCreated)
                    < errorIfLessThan
                ):
                    raise InsufficientPointGenerationError(
                        f"{self.permsCreated} permutations and only"
                        f" {min(estTotalTargetSizes)} points indicates {errorIfLessThan} null"
                        " distribution points will not be generated within"
                        f" {self.maxPerms} permutations."
                    )
                else:
                    morePermsFactor = numRepeats / min(estTotalTargetSizes) - 1
                nextChunksize = max(chunksize, int(self.permsCreated * morePermsFactor))
                nextChunksize = min(nextChunksize, self.maxPerms - self.permsCreated)
                needMorePerms = any(estTotalTargetSizes < numRepeats) and (
                    self.permsCreated < self.maxPerms
                )
                assert (
                    curNumPerms < self.maxPerms
                ), f"curNumPerms = {curNumPerms}, self.maxPerms = {self.maxPerms}"
        if self.permsCreated >= self.maxPerms:
            logger.warning(
                f"Produced {self.permsCreated} permutations with max {self.maxPerms}. A.shape ="
                f" {A.shape}, testName = {testName}"
            )
        # Assign no-cluster-permutation zeros proportionately to each null distribution bin
        if numNoClusterPerms > 0:
            for bn in nullDistributions.keys():
                nullDistributions[bn].extend([0] * (zeroDist == bn).sum())
        self.nullSizes = {bn: len(dist) for bn, dist in nullDistributions.items()}
        # Asign a p-value to each real cluster
        adjustedClusterStats = _test_adjust(rawClusterStats, test)
        assert rawClusterStats == [0] or len(adjustedClusterStats) == len(self.clusterBins), (
            f"rawClusterStats = {rawClusterStats}, adjustedClusterStats = {adjustedClusterStats},"
            f" realClusterBins = {self.clusterBins}"
        )
        # P-values cannot be zero, just less than 1 / size of null distribution
        self.clusterPVals = np.array([
            max(np.mean(nullDistributions[bn] > stat), 1 / len(nullDistributions[bn]))
            for stat, bn in zip(adjustedClusterStats, self.clusterBins)
        ])
        if clusterEndPoints:
            for start, end in np.array(clusterEndPoints)[self.clusterPVals <= clusterAlpha]:
                self.sigClusterPlotPoints[start:end] = 1
        if np.isnan(self.tested).any():
            self._fillBeyondTestWindow()

    def __sizeof__(self) -> int:
        return (
            getsizeof(self.sigClusterPlotPoints)
            + getsizeof(self.clusterPVals)
            + getsizeof(self.clusterEndpoints)
            + getsizeof(self.clusterDoFs)
            + getsizeof(self.stats)
            + getsizeof(self.numSamples)
            + getsizeof(self.nullSizes)
            + getsizeof(self.maxPerms)
            + getsizeof(self.permsCreated)
            + getsizeof(self.tested)
            + getsizeof(self.refIndex)
        )

    def _fillBeyondTestWindow(self):
        """
        Adjust sigClusterPlotPoints and stats to A.shape[1] size and fill in beyond testWindow with
        NaNs. Correct clusterEndpoints to reflect true size of A.shape[1].
        """
        sigClusterPlotPoints = np.full(self.tested.shape, np.nan)
        sigClusterPlotPoints[self.testWindow] = self.sigClusterPlotPoints
        self.sigClusterPlotPoints = sigClusterPlotPoints
        stts = np.full(self.tested.shape, np.nan)
        stts[self.testWindow] = self.stats
        self.stats = stts
        idxAdjustment = 0 if self.testWindow.start is None else self.testWindow.start
        self.clusterEndpoints = [
            (start + idxAdjustment, end + idxAdjustment) for start, end in self.clusterEndpoints
        ]

    def clusterDF(self) -> pd.DataFrame:
        """
        Return a DataFrame of the cluster stats, p-values, endpoints, and size of null used.

        Returns:
            pd.DataFrame: DataFrame with columns "pVal", "clusterStat", "startCluster",
                "endCluster", "ClusterDoF", "NullSize".
        """
        return pd.DataFrame({
            "pVal": self.clusterPVals,
            "clusterStat": [self.stats[start:end].sum() for start, end in self.clusterEndpoints],
            "startCluster": [x[0] for x in self.clusterEndpoints],
            "endCluster": [x[1] for x in self.clusterEndpoints],
            "ClusterDoF": self.clusterDoFs,
            "NullSize": [self.nullSizes[bn] for bn in self.clusterBins],
        })


def create_permutations(
    data: np.ndarray,
    statisticCalculator: Callable[[Array, Array], tuple[Vector, Vector]],
    rng: np.random.Generator,
    altData: np.ndarray = None,
    numRepeats=1000,
    maxWorkers: int = MAXWORKERS,
    **tqdmKwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Create permutations for paired statistical testing.

    If 2^(data.shape[0]) < numRepeats, all possible permutations will be generated.

    Args:
        data (np.ndarray): numElec by numTime 2D data array to test. If altData is not provided,
            permutated against zero (likely should be one-sided test). Otherwise, permutated
            against altData (likely should be two-sided test).
        statisticCalculator (Callable[[Array, Array], tuple[Vector,Vector]]): Function to calculate
            the test statistic. Should take two arrays of the same shape and return two vectors:
            the first of test statistics, the second of p-values.
        test (Literal['greater', 'lesser', 'two-sided']): Type of statistical test being performed.
        rng (np.random.Generator): An initialized generator to perform the permutations.
        altData (np.ndarray, optional): Optional second dataset for paired sample testing. If
            present, must be same shape as data. Defaults to None.
        numRepeats (int, optional): Target number of permutations to generate. Defaults to 1000.
        maxWorkers (int, optional): Maximum number of workers to use for parallel processing.

    Returns:
        tuple[np.ndarray, np.ndarray]: (allStats, allPVals), both number of permutations generated
            by numTime 2D arrays.
    """

    def indices_to_boolean(n: int, inds: list[int]) -> list[bool]:
        bools = [False] * n
        for ind in inds:
            bools[ind] = True
        return bools

    def all_boolean_masks(n: int) -> Generator[list[bool], Any, None]:
        for i in range(2**n):
            yield indices_to_boolean(n, [j for j in range(n) if (i & (1 << j))])

    if altData is None:
        altData = np.zeros_like(data)
    else:
        assert data.shape == altData.shape
    numElec = data.shape[0]
    if np.log2(numRepeats) > numElec:
        logger.info(
            f"Too many permutations requested, generating all possible {2**numElec} permutations."
        )
        numRepeats = 2**numElec
        permMask = np.array(list(all_boolean_masks(numElec)))
    else:
        # Produce numRepeats permutations: each electrode is flipped (entirely) with .5 chance
        permMask = rng.integers(2, size=(numRepeats, numElec), dtype="bool")
    logger.debug(
        f"Created permutation masks {permMask.shape} with grand mean"
        f" {permMask.mean():.4}, between {permMask.mean(axis=0).min():.4} and"
        f" {permMask.mean(axis=0).max()} across electrodes."
    )

    fixed_parallelize_create_permutations = partial(
        _parallelize_create_permutations, data, altData, statisticCalculator
    )
    # Optimize chunksize for large numRepeats (from https://stackoverflow.com/a/42096963/5722359)
    results = process_map(
        fixed_parallelize_create_permutations,
        permMask,
        max_workers=maxWorkers,
        chunksize=1 + (int(numRepeats / maxWorkers / 14) if numRepeats > 1000 else 0),
        **tqdmKwargs,
    )
    allStats = np.array([x[0] for x in results])
    allPVals = np.array([x[1] for x in results])
    return allStats, allPVals


def _parallelize_create_permutations(
    data: np.ndarray,
    altData: np.ndarray,
    statisticCalculator: Callable[[Array, Array], tuple[Vector, Vector]],
    mask: np.ndarray,
) -> tuple[Vector, Vector]:
    """
    Apply a shuffle mask to the data, calculate the test statistic, and return the results.
    """
    shuffleCopy = data - altData
    shuffleCopy[mask, :] *= -1
    stats, pVals = statisticCalculator(shuffleCopy, np.zeros_like(shuffleCopy))
    return stats, pVals


def SEM_within_subject(data: np.ndarray) -> np.ndarray:
    """
    Calculate the standard error of the mean within subjects.

    All rows with any NaN are silently dropped from the calculation.

    Adapted from Baria, A. T., Maniscalco, B., & He, B. J. (2017). Initial-state-dependent, robust,
    transient neural dynamics encode conscious visual perception. PLoS computational biology,
    13(11), e1005806.

    Args:
        data (np.ndarray): numSubjects by numConditions 2D array of data.

    Returns:
        np.ndarray: numConditions 1D array of SEM values.
    """
    cleanData = data[~np.isnan(data).any(axis=1), :]
    N, M = cleanData.shape
    normedData = cleanData - cleanData.mean(axis=1, keepdims=True) + np.mean(cleanData)
    varNormed = np.var(normedData, axis=0, ddof=1) * M / (M - 1)
    stdNormed = np.sqrt(varNormed)
    return stdNormed / np.sqrt(N)
