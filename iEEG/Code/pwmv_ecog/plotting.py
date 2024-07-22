import logging
import os
from typing import List, Sequence, TypeVar, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.markers as pltMarkers
import mne
import numpy as np
import pandas as pd
from tqdm import tqdm

from pwmv_ecog.channels import segment_ch_name
from pwmv_ecog.process import ECoGProcessing
from pwmv_ecog.triggers import EVENT_DICT, extract_trigger

logger = logging.getLogger(__name__)

Axes = mpl.axes.Axes
_SubjectIDs = TypeVar("_SubjectIDs", bound=str | int)

EPS = 1e-6  # Used for extracting tstats
MIN_TSTAT, MAX_TSTAT = 2.0, 5.0  # TODO: move to arguments

CARETLEFT = pltMarkers.MarkerStyle(pltMarkers.CARETLEFT)
CARETRIGHT = pltMarkers.MarkerStyle(pltMarkers.CARETRIGHT)
DEFAULT_TRIGGER_PLOT_SIZE = (12, 5)
# Change to list of named tuples
DEFAULT_PLOT_TYPES = [
    {"type": "raw", "vals": "raw", "desc": "Before Prep"},
    {"type": "raw", "vals": "raw_notched", "desc": "After Notch"},
    {"type": "raw", "vals": "rawFilteredRuns", "desc": "After Filtering"},
    {"type": "raw", "vals": "rawFiltDetrendRuns", "desc": "After Detrending"},
    {"type": "raw", "vals": "rawEKGCorrectedRuns", "desc": "After EKG Correction"},
    {"type": "raw", "vals": "preprocessedRuns", "desc": "After CAR"},
    {"type": "psd", "vals": "raw", "desc": "Before Prep"},
    {"type": "psd", "vals": "rawFilteredRuns", "desc": "After Filtering"},
    {"type": "psd", "vals": "rawEKGCorrectedRuns", "desc": "After EKG Correction"},
]


class Subjects(list[_SubjectIDs]):
    _rangeDelimiter = "-"
    _breakDelimiter = ","

    # Managing list of integer subject IDs with consistent shorthand printing
    def __init__(self, iterable, autoSort=True):
        super().__init__(iterable)
        self.autoSort = autoSort
        if self.autoSort:
            self.sort()

    def __iter__(self):
        copy = self.copy()
        if self.autoSort:
            copy.sort(key=lambda e: 1 if int(e) == 1001 else int(e))
        yield from copy

    def sort(self, **_):
        super().sort(key=lambda e: 1 if int(e) == 1001 else int(e))

    def __str__(self):
        # str always sorted. Use repr for actual current order
        copy = [str(e) for e in self.copy()]
        copy.sort(key=lambda e: 1 if int(e) == 1001 else int(e))
        base = self._breakDelimiter.join(copy)
        i_subj = 0
        while i_subj < len(copy):
            lastPrintedIdx = i_subj
            # Find and remove runs of consecutive integers
            while i_subj + 1 < len(copy) and int(copy[i_subj]) + 1 == int(copy[i_subj + 1]):
                i_subj += 1
            if i_subj - lastPrintedIdx > 1:
                toRemove = self._breakDelimiter.join(copy[lastPrintedIdx : i_subj + 1])[
                    len(copy[lastPrintedIdx]) : -len(copy[i_subj])
                ]
                base = base.replace(toRemove, self._rangeDelimiter)
            i_subj += 1

        return base


def plot_runstart_triggers(
    triggerdata: np.ndarray,
    timedata: np.ndarray,
    runstart_indices: List[int],
    pos_threshold=None,
    neg_threshold=None,
    normalize_on_y=False,
    show=True,
    figsize=DEFAULT_TRIGGER_PLOT_SIZE,
):
    """
    Plot the trigger data for each run with the aligned behavioral events

    :param triggerdata: trigger data
    :param timedata: time data
    :param runstart_indices: indices of run_start events in triggerdata/timedata frame (used for
        creating events)
    :param pos_threshold: threshold for positive value to get triggers
    :param neg_threshold: threshold for negative value to get triggers
    :param normalize_on_y: normalize the y axis
    :param show: show the plot
    """

    trigdata = triggerdata - np.mean(triggerdata)
    if normalize_on_y:
        trigdata = trigdata / np.max(trigdata)
    yval = np.max(trigdata) * 0.8

    fig = plt.figure(figsize=figsize)
    plt.plot(timedata, trigdata, alpha=0.9)
    plt.scatter(
        timedata[runstart_indices],
        [yval] * len(runstart_indices),
        c="c",
        marker="v",
        edgecolors="black",
    )
    if pos_threshold:
        pos_th = np.max(trigdata) * pos_threshold
        plt.axhline(pos_th, c="k", linestyle="--")
    if neg_threshold:
        neg_th = np.min(trigdata) * neg_threshold
        plt.axhline(neg_th, c="k", linestyle="--")
    plt.title(f"Trigger data found {len(runstart_indices)} run_start events")
    plt.xlabel("time (s)")
    plt.tight_layout()
    if show:
        plt.show()
    return fig


def plot_runs(
    timedata: np.ndarray,
    triggerdata: np.ndarray,
    runstart_indices: List[int],
    events: List[List[int]],
    sfreq: int,
    output_dir: str,
    runBuffer=10.0,
    prefix: str = None,
    drawTrialLines=False,
    figsize=DEFAULT_TRIGGER_PLOT_SIZE,
):
    """
    Plot the raw trigger data for each run with the aligned behavioral events

    :param timedata: time data
    :param triggerdata: trigger data
    :param runstart_indices: indices of run_start events (used for creating events)
    :param events: the generated MNE events. Indices must be in timedata/triggerdata frame.
    :param sfreq: sampling frequency
    :param output_dir: directory to save plots to
    :param prefix: prefix for the filename (could patient name)
    :param drawTrialLines: draw vertical lines at each trial start and end to confirm alignment
    """
    trial_start = ["Pre", "Post", "Grey", "mooney3"]

    for runid, (run_idx1, run_idx2) in enumerate(
        zip(runstart_indices, runstart_indices[1:] + [len(triggerdata)])
    ):
        runSlice = slice(int(run_idx1 - sfreq), run_idx2 - 1)
        yval = np.max(triggerdata[runSlice]) * 1.0

        trials_start_indices = [
            idx for idx, _, eid in events if sum([eid == EVENT_DICT[ts] for ts in trial_start]) > 0
        ]
        trials_end_indices = [idx for idx, _, eid in events if eid == EVENT_DICT["TrialEnd"]]

        fig = plt.figure(figsize=figsize)
        plt.plot(timedata, triggerdata)
        plt.scatter(timedata[runstart_indices], [yval] * len(runstart_indices), c="c", marker="v")
        plt.scatter(
            timedata[trials_start_indices],
            [yval] * len(trials_start_indices),
            c="r",
            marker=CARETRIGHT,
        )
        plt.scatter(
            timedata[trials_end_indices],
            [yval] * len(trials_end_indices),
            c="g",
            marker=CARETLEFT,
        )
        plt.ylim(np.min(triggerdata[runSlice]) * 1.1, yval * 1.1)
        if drawTrialLines:
            ymin, ymax = plt.ylim()
            trialLineStyle = {"lines.linestyle": "dotted"}
            with plt.rc_context(trialLineStyle):
                plt.vlines(timedata[trials_start_indices], ymin, ymax, colors="r", alpha=0.65)
                plt.vlines(timedata[trials_end_indices], ymin, ymax, colors="g", alpha=0.65)
            plt.ylim(ymin, ymax)
            del trialLineStyle

        trial_starts_cnt = len(
            [idx for idx in trials_start_indices if run_idx1 < idx and idx < run_idx2]
        )
        trial_ends_cnt = len(
            [idx for idx in trials_end_indices if run_idx1 < idx and idx < run_idx2]
        )
        assert trial_starts_cnt == trial_ends_cnt

        runstart_ind = [
            idx
            for idx, _, ev in events
            if run_idx1 < (idx + runBuffer * sfreq)
            and (idx + runBuffer * sfreq) < run_idx2
            and ev == EVENT_DICT["RunStart"]
        ][0]
        runend_ind = [
            idx
            for idx, _, ev in events
            if run_idx1 < (idx - runBuffer * sfreq)
            and (idx - runBuffer * sfreq) < run_idx2
            and ev == EVENT_DICT["RunEnd"]
        ][-1]
        plt.axvline(x=timedata[runstart_ind], color="r")
        plt.axvline(x=timedata[runend_ind], color="g")

        plt.xlim(timedata[runSlice.start], timedata[runSlice.stop])
        plt.title(
            "{} Trigger data for run {} (out of {}), {} trials".format(
                prefix.replace("_", " ") if prefix else "",
                runid + 1,
                len(runstart_indices),
                trial_starts_cnt,
            )
        )
        plt.xlabel("time (s)")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_dir,
                "{}run_{}-events.png".format(prefix + "_" if prefix else "", runid + 1),
            )
        )
        plt.close(fig)


def plot_triggers_outputdir_suffix(proc: ECoGProcessing, suffix: str, **kwargs):
    """
    Wrapper around plot_triggers to tag a suffix onto the default output_dir.
    """
    output_dir = os.path.join(
        proc.proj_vars["path_base"],
        proc.proj_vars["subjects"][proc.patient_name]["path_analysis_figures"],
        "TriggerAlignment" + suffix,
    )
    plot_triggers(proc, output_dir, **kwargs)


def plot_triggers(proc: ECoGProcessing, output_dir: Union[str, None] = None, drawTrialLines=True):
    """
    Plot the trigger data and save to output_dir

    :param proc: ECoGProcessing object
    :param output_dir: directory to save plots to
    """
    if not output_dir:
        plot_triggers_outputdir_suffix(proc, "", drawTrialLines=drawTrialLines)
        return
    os.makedirs(output_dir, exist_ok=True)

    fig = plot_runstart_triggers(
        proc.triggerdata, proc.timedata, proc.runstart_indices, show=False
    )
    plt.savefig(os.path.join(output_dir, "trigger-runs.png"), dpi=150)
    plt.close(fig)

    # Adjust event indices for trigger plotting
    adjustedEvents = np.array(proc.events) - [proc.si["index_TriggersStart"], 0, 0]
    plot_runs(
        proc.timedata,
        proc.triggerdata,
        proc.runstart_indices,
        adjustedEvents,
        proc.raw.info["sfreq"],
        output_dir,
        prefix=proc.patient_name,
        drawTrialLines=drawTrialLines,
    )

    if "Photodiode_chanlabel" in proc.si and proc.si["Photodiode_chanlabel"]:
        # Only if the channel exists
        pd_data, pd_timedata = extract_trigger(proc.raw, proc.si, proc.si["Photodiode_chanlabel"])
        plot_runs(
            pd_timedata,
            pd_data,
            proc.runstart_indices,
            adjustedEvents,
            proc.raw.info["sfreq"],
            output_dir,
            prefix=None if proc.patient_name is None else proc.patient_name + "_Photodiode",
            drawTrialLines=drawTrialLines,
        )

    # plot MNE events figure
    presentEvents = {
        k: v for k, v in EVENT_DICT.items() if v in pd.DataFrame(proc.events)[2].unique()
    }
    fig = mne.viz.plot_events(
        proc.events, event_id=presentEvents, sfreq=proc.raw.info["sfreq"], show=False
    )
    plt.savefig(os.path.join(output_dir, "MNE_Events.png"), dpi=150)
    plt.close(fig)

    logger.info("Finished plotting at:")
    logger.info(output_dir)


def _plot_run_raw(ax: Axes, rawRun: mne.io.Raw, ch: str) -> Axes:
    """
    Plot raw data for a single channel
    :param ax: matplotlib axes to plot on
    :param rawRun: raw data for a single run
    :param ch: channel name
    """
    ax.axhline(0, color="black", linestyle="dotted")
    data = rawRun.get_data(ch)[0]
    for badAnnot in filter(lambda a: a["description"] == "bad_acq_skip", rawRun.annotations):
        startBad = rawRun.time_as_index(badAnnot["onset"])[0] - rawRun.first_samp
        endBad = (
            rawRun.time_as_index(badAnnot["onset"] + badAnnot["duration"])[0] - rawRun.first_samp
        )
        data[startBad:endBad] = np.nan
    ax.plot(data)
    ax.tick_params(axis="x", labelbottom=False, top=False, bottom=False)
    ticks = [np.nanmin(data), np.nanmax(data)]
    color = "red" if (np.abs(ticks) > 1e-3).any() else "black"
    ax.set_yticks(
        ticks,
        [f"{tick:.1e}" for tick in ticks],
        rotation=90,
        va="center",
        fontsize="x-large",
        c=color,
    )
    return ax


def _plot_psd(ax: Axes, psd: mne.time_frequency.Spectrum, ch: str) -> Axes:
    """
    Plot the PSD of the raw data for a given channel and annotation
    :param ax: matplotlib axis
    :param raw: mne spectrum object
    :param ch: channel name to plot
    """
    psd.plot(axes=ax, picks=[ch], show=False, average=True)
    return ax


def plot_processing_channel(
    proc: ECoGProcessing,
    channel: str,
    output_path: str,
    plot_types: list[dict] = DEFAULT_PLOT_TYPES,
    rawYLimAnchorPlot: str = "preprocessedRuns",
) -> None:
    """
    Plot the raw data and PSD for a given channel, for all runs in the raw data.

    :param proc: ECoGProcessing object. All plot_types['vals'] expected to refer to attributes of
        proc that are each a list of raws, each raw corresponding to a run.

    :param channel: channel to plot
    :param output_path: path to save the plot
    :param plot_types: list of dict of plot types to plot, with keys 'type' ('raw' or 'psd'),
        'vals' (attribute in proc containing runs to plot), and 'desc' (title for subplot row).
    :param rawYLimAnchorPlot: Use this plot's ylim for all raw plots in each column.
    """
    ncols = len(proc.run_annotations)
    nrows = len(plot_types)
    fig, axr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
    rawYLims = np.zeros((ncols, 2))
    for row_id, (axc, pt) in enumerate(zip(axr, plot_types)):
        data: list = getattr(proc, pt["vals"])

        for ax, run in zip(axc, data):
            ann = run.annotations[0]  # Assume first annotation has run info
            col_id = int(run.annotations[0]["description"].removeprefix("Run ")) - 1
            if channel in run.ch_names:
                if pt["type"] == "raw":
                    ax = _plot_run_raw(ax, run, channel)
                    if pt["vals"] == rawYLimAnchorPlot:
                        rawYLims[col_id] = ax.get_ylim()
                elif (pt["type"] == "psd") and (
                    run.get_channel_types([channel])[0] in ["ecog", "seeg"]
                ):
                    ax = _plot_psd(ax, run, channel)

            title = ann["description"] if row_id == 0 else None
            left_ticks = col_id == 0 or pt["type"] == "raw"
            ax.tick_params(
                axis="y", labelleft=left_ticks, left=left_ticks, right=False, labelright=False
            )

            if title:
                ax.set_title(title, fontsize="xx-large")
            else:
                ax.set_title("")
            if col_id == 0:
                ax.set_ylabel(pt["desc"], rotation=85, fontsize="xx-large")
            else:
                ax.set_ylabel("")
    for row_id, pt in enumerate(plot_types):
        if pt["type"] == "raw":
            for col_id in range(ncols):
                axr[row_id, col_id].set_ylim(rawYLims[col_id])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_processing(
    proc: ECoGProcessing,
    output_dir: str,
    plot_types: list[dict] = DEFAULT_PLOT_TYPES,
    prefix="",
    skip_cached=False,
) -> None:
    """
    Plot the raw data for each preprocessed channel in the ECoGProcessing object.

    :param proc: ECoGProcessing object
    :param output_dir: Output directory
    :param plot_types: Dictionary of plot types to plot. Defaults to DEFAULT_PLOT_TYPES
    :param prefix: Prefix to add to the output filename
    :param skip_cached: Skip plotting if the output file already exists

    :return: None
    """
    if not proc.si["ApplyPulseArtifactRejection"]:
        # Skip EKG plots if no EKG correction done
        plot_types = [ptype for ptype in plot_types if ptype["vals"] != "rawEKGCorrectedRuns"]
    # Chop whole experiment raws into runs
    for i_p, ptype in enumerate(plot_types):
        if ptype["vals"].endswith("Runs"):
            continue
        setattr(
            proc,
            f"{ptype['vals']}Runs",
            getattr(proc, ptype["vals"]).crop_by_annotations(proc.run_annotations),
        )
        plot_types[i_p]["vals"] = f"{ptype['vals']}Runs"
    # Run PSD once per run, if needed
    for i_p, ptype in enumerate(tqdm(plot_types, desc="PSD Creation")):
        if ptype["type"] == "psd":
            raws: list[mne.io.Raw] = getattr(proc, ptype["vals"])
            psds = []
            for raw in raws:
                psds.append(raw.compute_psd(fmax=200, verbose="error"))
                # Need annotations to track run number. Probably lighterweight way to do this...
                setattr(psds[-1], "annotations", raw.annotations.copy())
            setattr(proc, f"PSD{ptype['vals']}", psds)
            plot_types[i_p]["vals"] = f"PSD{ptype['vals']}"
    for ch in tqdm(proc.ch_names, desc="Plot Generation", unit="Channel"):
        output_path = os.path.join(output_dir, f"{prefix}Processing-{ch}.png")
        if skip_cached and os.path.isfile(output_path):
            continue
        plot_processing_channel(proc, ch, output_path, plot_types=plot_types)


def generate_video(ddir, outpath):
    import cv2

    video_name = os.path.join(ddir, outpath)
    images = sorted([img for img in os.listdir(ddir) if img.endswith(".png")])

    fourcc = cv2.VideoWriter_fourcc(*"MP4V")

    # Array images should only consider
    # the image files ignoring others if any
    frame = cv2.imread(os.path.join(ddir, images[0]))

    # setting the frame width, height width
    # the width, height of first image
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name, fourcc, 1, (width, height))

    # Appending the images to the video one by one
    for image in images:
        video.write(cv2.imread(os.path.join(ddir, image)))

    # Deallocating memories taken for window creation
    cv2.destroyAllWindows()
    video.release()  # releasing the video generated


def middleLabelOnly(labels: list[str]) -> list[str]:
    """
    Produce label list with one label per electrode placed in the middle.

    Args:
        labels (list[str]): List of channel labels. Assumes channel digits are always
            contiguous and starting from 1.

    Returns:
        list[str]: Label list for use in plotting the channels.
    """
    labelsDF = pd.DataFrame(
        (pd.Series(labels).apply(lambda label: segment_ch_name(label, verbose=False)).to_list()),
        columns=["elec", "chNum", "suffix"],
    ).drop(columns="suffix")
    labelsDF["chNum"] = labelsDF.chNum.astype(int)
    numChannels = labelsDF.groupby("elec").count()
    labelPos = {label: (num + 1) // 2 for label, num in numChannels.itertuples()}
    return [
        (elec if chNum == labelPos[elec] else "")
        for elec, chNum in labelsDF.itertuples(index=False)
    ]


def meanROIs_figure(
    roiData: pd.DataFrame,
    ax: plt.Axes,
    roiName: str,
    numElecs: int,
    yUnit: str,
    condsToPlot=["grey", "pre", "post"],
    meanLineColors={"pre": (0, 0, 1), "post": (1, 0, 0), "grey": (0.5, 0.5, 0.5)},
    semAlpha: Union[float, dict[str, float], Sequence[float]] = 0.2,
):
    """
    Paint mean ROI figure on ax.

    Args:
        roiData (pd.DataFrame): Data to plot. Expects MultiIndex with (cond, "mean") and (cond,
            "sem") to be valid for each cond in condsToPlot. Columns are the timepoints to plot.
        ax (plt.Axes): Axes on which to paint.
        roiName (str): ROI name to use in plot title.
        numElecs (int): Number of electrodes to list in plot title.
        yUnit (str): Y-axis label.
        condsToPlot (list, optional): Conditions to plot. Defaults to ["pre", "post", "grey"].
        meanLineColors (dict, optional): Colors to use for each condition. Can specify alpha
            individually by specifying a fourth value in each color tuple. Defaults to {"pre":
            (0, 0, 1), "post": (1, 0, 0), "grey": (0.5, 0.5, 0.5)}.
        semAlpha (Union[float, dict[str, float], Sequence[float]], optional): Transparency for the
            SEM clouds. Can specify a scalar to apply to all clouds, a list in the same order
            as condsToPlot, or a dict explicitly mapping a value to each cond. Defaults to 0.2.
    """
    if isinstance(semAlpha, float):
        semAlpha = {cond: semAlpha for cond in condsToPlot}
    elif isinstance(semAlpha, Sequence):
        semAlpha = {cond: semAlpha[i] for i, cond in enumerate(condsToPlot)}

    columns = roiData.columns
    ax.axvline(0, c="k", linestyle=":")
    ax.axhline(0, c="k", linestyle=":")

    for cond in condsToPlot:
        means = roiData.loc[(cond, "mean"), :]
        sems = roiData.loc[(cond, "sem"), :]
        meanColor = meanLineColors[cond]
        ax.fill_between(
            list(columns),
            list(means - sems),
            list(means + sems),
            color=meanColor[:3] if type(meanColor) is tuple else meanColor,
            alpha=semAlpha[cond],
        )
        ax.plot(columns, means, color=meanColor, label=cond.capitalize(), linewidth=2.0)

    title = roiName + f" ({numElecs} elecs)" if numElecs > 1 else roiName
    ax.set_title(title, fontsize="large")
    ax.set_ylabel(yUnit)
    ax.set_xlabel("seconds")
    ax.grid(axis="x")

    ax.legend(labelspacing=0.3, loc="lower left")


def plot_sig_bar(
    ax: plt.Axes,
    boolData: np.ndarray,
    color: str | tuple[float, float, float] | tuple[float, float, float, float],
    xData: Sequence[float] = None,
    label: str = "",
):
    """
    Plot new significance bar across top of ax.

    Args:
        ax (plt.Axes): Axes to plot on.
        boolData (np.ndarray): Positions to plot the bar. Should be boolean, or only 1 and 0/NaN.
        color (str | tuple[float, float, float, float]): Any valid matplotlib color code.
        xData (Sequence[float], optional): X-axis values matching boolData if needed. Defaults is
            to plot boolData without any x-axis specification.
        label (str, optional): Label for the bar. Defaults to no label.
    """
    boolData[boolData == False] = np.nan
    boolData[boolData == 0] = np.nan
    boolData[~np.isnan(boolData)] = 1
    if any(boolData):
        _, ymax = ax.get_ylim()
        sigPoints = boolData * ymax
        if xData is not None:
            ax.plot(xData, sigPoints, ".", color=color, label=label if label else None)
        else:
            ax.plot(sigPoints, ".", color=color, label=label if label else None)


def scatterBox(
    data: pd.DataFrame,
    markers: Sequence[str] = None,
    columnColors: Sequence = None,
    linestyle="-",
    linecolor=(0.0, 0.0, 0.0, 0.3),
):
    """
    Boxplot + scatter of each column, with lines connecting the row-points in each column plot.

    Args:
        data (pd.DataFrame): Rows are observations, columns are conditions to plot.
        markers (Sequence[str], optional): Markers to use for each row. If provided, must be the
            same length as rows of data. Defaults to None.
        columnColors (Sequence, optional): Colors to use in the boxplot lines and scatter markers
            for each column. If provided, must be same length as columns of data. Defaults to using
            matplotlib's default colors.
        linestyle (str, optional): Linestyle for connecting lines. Defaults to "-".
        linecolor (tuple, optional): Color for connecting lines, including alpha. Defaults to
            (0.0, 0.0, 0.0, 0.3).
    """
    if not columnColors:
        columnColors = [None] * data.shape[1]
    boxplot = data.boxplot(figsize=(12, 9), showfliers=False, fontsize=22, return_type="dict")
    for lineGroup in boxplot.values():
        # Some groups have more than one line per boxplot. All lines for same boxplot are together,
        # so need to grab runs of lines by stride
        stride = len(lineGroup) // len(columnColors)
        for i_color, color in enumerate(columnColors):
            for line in lineGroup[i_color * stride : (i_color + 1) * stride]:
                line.set_color(color)
    plt.ylabel("Accuracy", fontsize=18)
    Xs = np.zeros(data.shape)
    for i, d in enumerate(data):
        y = data[d]
        Xs[:, i] = np.random.normal(i + 1, 0.01, len(y))
        if markers is None:
            plt.scatter(Xs[:, i], y)
        else:
            assert data.shape[0] == len(markers), "Must be same number of markers as subjects."
            scatter = plt.scatter(Xs[0, i], y.iloc[0], marker=markers[0], color=columnColors[i])
            for i_subj, subjID in enumerate(y.index[1:]):
                i_subj += 1
                plt.scatter(
                    Xs[i_subj, i], y[subjID], marker=markers[i_subj], color=scatter.get_fc()
                )
    for i, row in enumerate(Xs):
        plt.plot(row, data.iloc[i, :], ls=linestyle, c=linecolor)
    autoYLims = plt.ylim()
    plt.ylim([min(0, autoYLims[0]), max(1, autoYLims[1])])
