from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pylab import Rectangle

from bjhlab_ecog_tools.plotting.coverage import merged_hemispheres_coverage_plot, VALID_SIDES
from pwmv_ecog.analysis.mean_roi import MeanROIAnalysis
from pwmv_ecog.analysis.image_preference import (
    ImagePreferenceAnalysis,
    imagePrefs_figure,
    imagePrefs_survival_plot,
)
from pwmv_ecog.plotting import meanROIs_figure, plot_sig_bar, scatterBox, Subjects

SUBJECTS: Subjects = Subjects(
    [1001, 2, 4, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] + [3, 5, 9]
)
subjectCMap = plt.colormaps["tab20"]
SUBJECTCOLORS = [subjectCMap(i) for i in np.linspace(0, 1, len(SUBJECTS))]

# CONDITIONCOLORS = {"pre": "#1f77b4", "post": "#ff7f0e", "grey": "#2ca02c"}
CONDITIONCOLORS = {"pre": "#1f77b4", "post": "#ff7f0e", "grey": "0.5"}
SIGBARCOLORS = {"Pre>0": "#1f77b4", "Post>0": "#ff7f0e", "Post>Pre": "k"}

ROICOLORS = {
    "EVC": "#00FF00",
    "HLVC": "#FF0000",
    "Dorsal": "#0000FF",
    "FPN": "#FF00FF",
    "DMN": "#00FFFF",
    "Limbic": "#862D2D",
    "other": "#000000",
}

# These channels appear outside the brain and should be omitted from the coverage plot
# Already not in any ROI
OMITCHANNELS: list[tuple[int, str]] = [
    (3, "O1"),
    (1001, "FP3"),
    (1001, "FP2"),
    (3, "SF1"),
    (3, "SF2"),
    (5, "RAF2"),
    (5, "RAF3"),
    (5, "RPF2"),
    (5, "RPF3"),
    (3, "SF3"),
    (3, "SF4"),
    (3, "PT3"),
    (17, "STP1.5"),
    (3, "PT1"),
    (10, "SFOF1.5"),
]
# These channels have more than two ROIs and will only be shown as being Dorsal-FPN for coverage
# plot since the DMN ROI analysis is only shown in the supplementary
TRIPLECHANS: dict[str, list[tuple[int, str]]] = {
    "DMN-Dorsal-FPN": [
        (1001, "PO3"),
        (1001, "PO4"),
        (2, "RC3"),
        (2, "RC4"),
        (13, "LPCPP7.5"),
        (13, "LPCPP8.5"),
        (16, "LSPPC9.5"),
        (16, "RPIOC11.5"),
        (16, "RSPPC8.5"),
        (16, "RSPPC9.5"),
    ]
}

FIGSPATH = Path("..", "Analysis", "Figures", "Manuscript2024")
BEHAVPATH = Path("pwmv_ecog", "manuscript", "manuscript_behav.pd.pkl")


def supp_table_3(subjects=SUBJECTS, outPath=Path(FIGSPATH, "SuppTable3")) -> pd.Series:
    """Electrode details table."""
    meanROI: MeanROIAnalysis = MeanROIAnalysis(subjects).load()
    ipa: ImagePreferenceAnalysis = ImagePreferenceAnalysis(subjects).load()
    nonExcludedBySubject = pd.concat(
        [
            analysis.data.loc[:, 0]
            .groupby(["SubjectID", "Electrode"])
            .max()
            .groupby("SubjectID")
            .count()
            for analysis in [meanROI, ipa]
        ],
        axis=1,
    )  # type: ignore
    assert (
        nonExcludedBySubject.iloc[:, 0] - nonExcludedBySubject.iloc[:, 1] == 0
    ).all(), "Mean ROI & ImagePref analyses have differing number of included electrodes."
    nonExcluded = meanROI.data.loc[:, 0].groupby(["SubjectID", "Electrode"]).max().count()
    includedChannels = meanROI.memberships.any(axis=1)
    analyzedChannels = (
        meanROI.data.loc[includedChannels]
        .loc[:, 0]
        .groupby(["SubjectID", "Electrode"])
        .max()
        .count()
    )
    noNA = meanROI.data.loc[:, 0].dropna()
    roiElecIndexs = {}
    for roi in meanROI.memberships.columns:
        roiElecIndexs[roi] = (
            noNA[
                noNA.index.droplevel(["Condition"]).isin(
                    meanROI.memberships[meanROI.memberships[roi]].index
                )
            ]
            .index.droplevel(["Condition"])
            .drop_duplicates()
        )
    meanROICounts = pd.Series({roi: len(roiIdx) for roi, roiIdx in roiElecIndexs.items()})
    twoROIIntersects = {}
    threeROIIntersects = {}
    for i_roi, roi1 in enumerate(meanROI.memberships.columns):
        for i_roi2, roi2 in enumerate(meanROI.memberships.columns[i_roi + 1 :]):
            twoIntersect = roiElecIndexs[roi1].intersection(roiElecIndexs[roi2])
            if len(twoIntersect):
                twoROIIntersects[f"{roi1} and {roi2}"] = len(twoIntersect)
            for roi3 in meanROI.memberships.columns[i_roi + 1 + i_roi2 + 1 :]:
                threeIntersect = twoIntersect.intersection(roiElecIndexs[roi3])
                if len(threeIntersect):
                    threeROIIntersects[" and ".join([roi1, roi2, roi3])] = len(threeIntersect)
    twoROICounts = pd.Series(twoROIIntersects).sort_values(ascending=False)
    threeROICounts = pd.Series(threeROIIntersects).sort_values(ascending=False)
    assert (
        meanROI.memberships.sum(axis=1) > 3
    ).sum() == 0, "More than 3 ROI memberships are not counted"
    assert meanROICounts.sum() - twoROICounts.sum() + threeROICounts.sum() == analyzedChannels
    suppTable = pd.concat([
        pd.Series([nonExcluded, analyzedChannels], index=["All recorded", "All analyzed"]),
        pd.Series(meanROICounts).sort_values(ascending=False),
        twoROICounts,
        threeROICounts,
    ]).rename("Electrode Counts")
    suppTable.index = suppTable.index.str.replace("Dorsal", "Dorsal Stream")
    outPath.mkdir(parents=True, exist_ok=True)
    suppTable.to_csv(outPath / "SuppTable3.csv")
    print("SuppTable3.csv saved.")
    return suppTable


def supp_fig_6(
    subjects=SUBJECTS,
    figPath=Path(FIGSPATH, "SuppFig6"),
    coords_path=Path("Electrodes", "ROI_Masks", "network", "networkCoords.pd.pkl"),
) -> plt.Figure:
    def axes_letter(ax: plt.Axes, letter: str, xy=(0, 1.1), xycoords="axes fraction"):
        ax.annotate(letter, xy=xy, xycoords=xycoords, fontsize=28, weight="bold")

    figPath.mkdir(parents=True, exist_ok=True)
    figWidth = 8.5 * 2
    nrows = int(np.ceil(len(subjects) / 3)) + 1
    ncols = 3
    fig = plt.figure(layout="constrained", figsize=(figWidth, figWidth / (ncols * 2) * nrows))
    gs = GridSpec(nrows, ncols, figure=fig, height_ratios=[1] * (nrows - 1) + [1.5])
    # Trial timeline plot
    trialTimeline = fig.add_subplot(gs[-2, 1:])
    trialTimeline.set_anchor("E")
    l, b, w, h = trialTimeline.get_position().bounds
    trialTimeline.set_position((l + 0.09, b - 0.12, w, h + 0.12))
    trialTimeline.imshow(plt.imread(figPath / "trialTimeline.png"))
    axes_letter(trialTimeline, "b", xy=(-0.1, 0.85))
    trialTimeline.axis("off")
    # Don't have a title to accomodate, so increase the height used for the image
    # Task performance plots
    with plt.rc_context({
        "font.size": 16,
        "axes.titlesize": 28,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "ytick.right": False,
    }):
        taskPerfFig = fig.add_subfigure(gs[-1, :])
        originalPerfAx, catchPerfAx, blank = taskPerfFig.subplots(1, 3)
        blank.axis("off")
        task_performance_plot(subjects, originalPerfAx, catchPerfAx)
        originalPerfAx.set_title("Original Trials")
        axes_letter(originalPerfAx, "c", xy=(-0.195, 1.1))
        catchPerfAx.set_title("Catch Trials")
    # Subject coverage plots
    coords: pd.DataFrame = pd.read_pickle(coords_path)
    # Clear ROI data since plotting single color per patient
    coords.loc[:, "DMN":] = False
    coords.loc[:, "EVC"] = True
    with plt.rc_context({"font.size": 16, "axes.titlesize": 22}):
        for i_subj, subjID in enumerate(subjects):
            subjCoords = coords[coords.index.get_level_values(0) == subjID]
            coverageFig = fig.add_subfigure(gs[i_subj // 3, i_subj % 3])
            coverageAxs: np.ndarray = coverageFig.subplots(1, 2)
            if i_subj == 0:
                axes_letter(coverageAxs[0], "a", xy=(0, 1))
            subjColor = {roi: SUBJECTCOLORS[i_subj] for roi in ROICOLORS.keys()}
            manuscript_coverage_plot(coverageAxs, subjCoords, subjColor, axesWidth=figWidth / 6)
            coverageFig.text(0.5, 0.85, str(i_subj + 1), horizontalalignment="center", fontsize=28)
    fig.savefig(figPath / "SuppFig6.png")
    print("SuppFig6.png saved.")
    return fig


def supp_fig_7(
    subjects=SUBJECTS,
    figPath=Path(FIGSPATH, "SuppFig7"),
    roisToPlot=["EVC", "HLVC", "Dorsal", "FPN", "Limbic", "DMN"],
) -> plt.Figure:
    figPath.mkdir(parents=True, exist_ok=True)
    figWidth = 8.5 * 2
    nrows, ncols = 3, 6
    fig = plt.figure(layout="constrained", figsize=(figWidth, figWidth / ncols * nrows))
    gs = GridSpec(nrows, ncols, figure=fig)
    plotSettings = {
        "font.size": 12,
        "axes.titlesize": 20,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "ytick.right": False,
    }

    # RT histograms
    behav: pd.DataFrame = pd.read_pickle(BEHAVPATH)
    with plt.rc_context(plotSettings):
        toBoxAxs = []
        for i_ax, cond in enumerate([
            ("pre", "rec"),
            ("pre", "un"),
            ("post", "rec"),
            ("post", "un"),
            ("grey", "rec"),
            ("grey", "un"),
        ]):
            # Share axes
            if i_ax == 0:
                rtHist = fig.add_subplot(gs[0, i_ax])
                rtRefAx = rtHist
            else:
                rtHist = fig.add_subplot(gs[0, i_ax], sharey=rtRefAx)
            (
                behav.loc[
                    behav[cond[0]] & (behav.OverrideACC == int(cond[1] == "rec")), "ReactionTime"
                ]
                / 1000  # Convert from ms to s
            ).plot.hist(ax=rtHist, color=CONDITIONCOLORS[cond[0]], xlabel="seconds", bins=10)
            perfStr = "Recog." if cond[1] == "rec" else "Unrec."
            condStr = "Gray" if cond[0] == "grey" else cond[0].capitalize()
            rtHist.set_title(" ".join([perfStr, condStr]))
            if cond in [("pre", "un"), ("post", "rec"), ("grey", "rec")]:
                toBoxAxs.append(rtHist)
        # Draw green box around the conditions used in analyses
        for rtHist in toBoxAxs:
            axBounds = rtHist.axis()
            rec = Rectangle(
                (axBounds[0] - 0.2, axBounds[2] - 2.2),
                (axBounds[1] - axBounds[0]) + 0.3,
                (axBounds[3] - axBounds[2]) + 4,
                fill=False,
                lw=2,
                color="g",
            )
            rec = rtHist.add_patch(rec)
            rec.set_clip_on(False)
    fig.text(0.005, 0.95, "a", fontsize=28, weight="bold")

    # Mean ROI survival plots
    meanROIResults: MeanROIAnalysis = MeanROIAnalysis(subjects).load()
    with plt.rc_context(plotSettings):
        for i_ax, roi in enumerate(roisToPlot):
            meanROIAx = fig.add_subplot(gs[1, i_ax])
            manuscript_meanROI_survival_plot(meanROIResults, roi, meanROIAx)
            if i_ax == 0:
                meanROIAx.legend(loc="lower left")
                meanROIAx.set_ylabel("# channels", fontsize=16)
            else:
                meanROIAx.set_ylabel("")
            meanROIAx.tick_params(labelbottom=False)
            meanROIAx.set_xlabel("")
    fig.text(0.005, 0.6, "b", fontsize=28, weight="bold")

    # Image preference survival plots
    imagePrefResults: ImagePreferenceAnalysis = ImagePreferenceAnalysis(subjects).load()
    with plt.rc_context(plotSettings):
        for i_ax, roi in enumerate(roisToPlot):
            imagePrefAx = fig.add_subplot(gs[2, i_ax])
            manuscript_imagePref_survival_plot(imagePrefResults, roi, imagePrefAx)
            if i_ax == 0:
                imagePrefAx.legend(loc="lower left")
                imagePrefAx.set_ylabel("# channels", fontsize=16)
            else:
                imagePrefAx.set_ylabel("")
            imagePrefAx.set_title("")
    fig.text(0.005, 0.3, "c", fontsize=28, weight="bold")
    fig.savefig(figPath / "SuppFig7.png")
    print("SuppFig7.png saved.")
    return fig


def supp_fig_8(subjects=SUBJECTS, figPath=Path(FIGSPATH, "SuppFig8")) -> plt.Figure:
    figPath.mkdir(parents=True, exist_ok=True)
    roisToPlot = ["Limbic", "DMN"]
    figWidth = 8.5 * 2 / 2
    with plt.rc_context({
        "font.size": 16,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "xtick.labelsize": "small",
        "ytick.right": False,
        "ytick.labelsize": "small",
    }):
        # Figure height is the width of one results subplot multiplied by number of rows
        fig = plt.Figure((figWidth, figWidth / 2 * 2), constrained_layout=True)
        manuscript_results_rows(fig, roisToPlot, subjects=subjects)
    fig.savefig(figPath / "SuppFig8.png")
    print("SuppFig8.png saved.")
    return fig


def fig_4(
    subjects=SUBJECTS,
    roisToPlot=["EVC", "HLVC", "Dorsal", "FPN"],
    coords_path=Path("networkCoords.pd.pkl"),
    figPath=Path(FIGSPATH, "Fig4"),
):
    assert len(roisToPlot) == 4, f"Only the 4-wide ROI design is possible. Got {roisToPlot}."
    figPath.mkdir(parents=True, exist_ok=True)
    figWidth = 8.5 * 2
    fig = plt.figure(layout="constrained", figsize=(figWidth, figWidth / 3 * 2 + 1))
    gs = GridSpec(3, 4, figure=fig, height_ratios=[4, 3, 3])
    topRow = fig.add_subfigure(gs[0, :])
    topRowGS = GridSpec(1, 3, figure=topRow, width_ratios=[2, 2, 7])
    imagePrefSchematic = fig.add_subplot(topRowGS[0, -1])
    imagePrefSchematic.set_anchor("W")
    imagePrefSchematic.imshow(plt.imread(Path(figPath, "ImagePrefSchematic.png")))
    imagePrefSchematic.axis("off")
    fig.text(4 / 11, 0.95, "b", fontsize=28, weight="bold")

    coords: pd.DataFrame = pd.read_pickle(coords_path)
    coverageFig = fig.add_subfigure(topRowGS[0, :2])
    coverageAxs = coverageFig.subplots(1, 2)
    with plt.rc_context({"font.size": 16, "axes.titlesize": 22}):
        manuscript_coverage_plot(coverageAxs, coords, axesWidth=figWidth / 4)
        handles, labels = coverageAxs[1].get_legend_handles_labels()
        roiSortOrder = ["EVC", "HLVC", "Dorsal", "FPN", "Limbic", "DMN"]
        sortedHandles = [handles[labels.index(roiName)] for roiName in roiSortOrder]
        sortedLabels = [labels[labels.index(roiName)] for roiName in roiSortOrder]
        coverageFig.legend(
            sortedHandles,
            sortedLabels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.27),
            ncol=2,
            markerscale=1.5,
            fontsize=12,
            handlelength=0.5,
            handletextpad=0.4,
            columnspacing=1,
        )
        fig.text(0.005, 0.95, "a", fontsize=28, weight="bold")

    with plt.rc_context({
        "font.size": 16,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "xtick.labelsize": "x-small",
        "ytick.right": False,
        "ytick.labelsize": "x-small",
    }):
        resultsRows = fig.add_subfigure(gs[1:3, :])
        manuscript_results_rows(resultsRows, roisToPlot, subjects=subjects)
        fig.text(0.005, 0.57, "c", fontsize=28, weight="bold")
    fig.savefig(figPath / "Fig4.png")
    print("Fig4.png saved.")
    return fig


def manuscript_results_rows(fig: plt.Figure, roisToPlot: list[str], subjects=SUBJECTS):
    meanROIResults: MeanROIAnalysis = MeanROIAnalysis(subjects).load()
    imagePrefResults: ImagePreferenceAnalysis = ImagePreferenceAnalysis(subjects).load()
    meanROIAxs = []
    imagePrefAxs = []
    gs = GridSpec(3, len(roisToPlot), figure=fig, height_ratios=[4, 4, 1])
    for i_ax, roi in enumerate(roisToPlot):
        # Mean ROI plot
        roiCondData = meanROIResults.results.loc[roi]
        roiSigPoints = meanROIResults.sigBars.loc[roi]
        meanROIAxs.append(fig.add_subplot(gs[0, i_ax]))
        manuscript_meanROI_plot(roiCondData, roiSigPoints, meanROIAxs[-1], roi)
        # Image Pref plot
        roiData = imagePrefResults.results.loc[roi]
        imagePrefAxs.append(fig.add_subplot(gs[1, i_ax]))
        manuscript_imagePref_plot(roiData, imagePrefAxs[-1], roi, imagePrefResults.toi)
    imagePrefAxs[0].set_ylabel("Grayscale Similarity (log[µV²])", fontsize=14)
    for ax in meanROIAxs[1:] + imagePrefAxs[1:]:
        ax.set_ylabel("")
    for ax in meanROIAxs:
        ax.tick_params(labelbottom=False)
        ax.set_xlabel("")
        if ax.get_title() == "Dorsal":
            ax.set_title("Dorsal Stream")
    for ax in imagePrefAxs:
        ax.set_title("")
    for ax in meanROIAxs + imagePrefAxs:
        ax.locator_params(nbins=5)
        # ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, -2), useMathText=True)
        ax.set_xticks([0.25, 0.75, 1.25], minor=True)

        ax.grid(True, which="minor", axis="x")
    handles, labels = meanROIAxs[0].get_legend_handles_labels()
    handles[labels.index("Post>Pre")] = plt.Line2D([], [], lw=4, c=SIGBARCOLORS["Post>Pre"])
    handles[labels.index("Pre>0")] = plt.Line2D([], [], lw=4, color=SIGBARCOLORS["Pre>0"])
    handles[labels.index("Post>0")] = plt.Line2D([], [], lw=4, color=SIGBARCOLORS["Post>0"])
    handles[labels.index("Pre")] = (
        meanROIAxs[0].fill(np.nan, np.nan, CONDITIONCOLORS["pre"], alpha=0.2)[0],
        handles[labels.index("Pre")],
    )
    handles[labels.index("Post")] = (
        meanROIAxs[0].fill(np.nan, np.nan, CONDITIONCOLORS["post"], alpha=0.2)[0],
        handles[labels.index("Post")],
    )
    legendFig = fig.add_subfigure(gs[2, :])
    legendAxTrace, legendAxSig = legendFig.subplots(1, 2)
    for legendAx in [legendAxTrace, legendAxSig]:
        legendAx.axis("off")
    legendParams = {"ncol": 3, "fontsize": "small", "title_fontsize": "small"}
    if len(roisToPlot) <= 2:
        legendParams["handlelength"] = 1.0
        legendParams["columnspacing"] = 1.0
    else:
        legendParams["handlelength"] = 2.0
        legendParams["columnspacing"] = 2.0
    legendAxTrace.legend(
        handles[:3], labels[:3], loc="center right", title="mean, paired SEM", **legendParams
    )
    legendAxSig.legend(
        handles[3:], labels[3:], loc="center left", title="p<0.05, cluster-corr.", **legendParams
    )


def manuscript_imagePref_plot(
    roiData: pd.DataFrame, ax: plt.Axes, roiFullName: str, toi
) -> plt.Axes:
    imagePrefColors = {
        "Pre": SIGBARCOLORS["Pre>0"],
        "Post": SIGBARCOLORS["Post>0"],
        "Post - Pre": SIGBARCOLORS["Post>Pre"],
    }
    imagePrefs_figure(
        roiData,
        ax,
        roiFullName,
        "Grayscale Similarity (log[µV²])",
        toi,
        colors=imagePrefColors,
        alpha=0.2,
        semToUse="Post - Pre",
    )
    if (legend := ax.get_legend()) is not None:
        legend.remove()
    return ax


def manuscript_imagePref_survival_plot(
    imagePrefResults: ImagePreferenceAnalysis, roiName: str, ax: plt.Axes
):
    roiImgPrefData = imagePrefResults.results.loc[roiName]
    imagePrefs_survival_plot(
        roiImgPrefData.loc[("Pre", "numElec")],
        roiImgPrefData.loc[("Post", "numElec")],
        ax,
        roiName,
        toi=imagePrefResults.toi,
        ls=["-", "--"],
        lw=[3, 3],
    )
    if (legend := ax.get_legend()) is not None:
        legend.remove()
    _, top = ax.get_ylim()
    ax.set_ylim(bottom=0, top=top * 1.01)


def manuscript_meanROI_survival_plot(
    meanROIResults: MeanROIAnalysis, roiName: str, ax: plt.Axes, colors=CONDITIONCOLORS
):
    roiCounts = meanROIResults.results.loc[roiName].xs("numElec", level="Stat")
    pd.DataFrame(
        {"Gray": roiCounts.loc["grey"], "Pre": roiCounts.loc["pre"], "Post": roiCounts.loc["post"]}
    ).plot(
        ax=ax,
        color=[colors["grey"], colors["pre"], colors["post"]],
        style=["-", "--", ":"],
        linewidth=3,
        grid=True,
        xlabel="seconds",
        ylabel="surviving channels",
        title=roiName,
    )
    if (legend := ax.get_legend()) is not None:
        legend.remove()
    _, top = ax.get_ylim()
    ax.set_ylim(bottom=0, top=top * 1.01)


def manuscript_meanROI_plot(
    roiCondData: pd.DataFrame, condSigPoints: pd.DataFrame, ax: plt.Axes, roiFullName: str
) -> plt.Axes:
    """
    Mean ROI plot for paper. Paints post, pre, grey mean trace with post and pre SEM clouds.

    Args:
        roiCondData (pd.DataFrame): Data to plot. Expects MultiIndex with (cond, "mean") and (cond,
            "sem") for each cond in ["grey", "pre", "post"]. Columns are the timepoints to plot.
        condSigPoints (pd.DataFrame): Positions of significant points to plot. Each row should be
            boolean, or only 1 and 0/NaN. Index values must all be in SIGBARCOLORS.keys(). Columns
            are the timepoints to plot.
        ax (plt.Axes): Axes on which to paint.
        roiFullName (str): ROI name to use in plot title.

    Returns:
        plt.Axes: Drawn axes.
    """
    # Ensure no grey SEM cloud is drawn
    noGreySEM = roiCondData.copy()
    noGreySEM.loc[("grey", "sem"), :] = 0
    meanROIs_figure(
        noGreySEM,
        ax,
        roiFullName,
        numElecs=0,  # No NumElecs in title
        yUnit="HGP (log[µV²])",
        condsToPlot=["grey", "pre", "post"],
        meanLineColors=CONDITIONCOLORS,
        semAlpha=0.2,
    )
    for cond, sigPoints in condSigPoints.iterrows():
        plot_sig_bar(
            ax,
            sigPoints.to_numpy(),
            color=SIGBARCOLORS[str(cond)],
            label=str(cond),
            xData=sigPoints.index.astype(float).to_list(),
        )
    # legend handled at overall figure level
    ax.get_legend().remove()
    return ax


def manuscript_coverage_plot(
    axs: np.ndarray,
    coords: pd.DataFrame,
    roiColors=ROICOLORS,
    axesWidth=8.5 / 4,
    noROIColor="black",
    displaySide: Literal["r", "l"] = "r",
) -> np.ndarray:
    f"""
    Prep and generate coverage plot for paper.

    Args:
        axs (np.ndarray): Array of axes to plot on. Should be 1D with 2 axes.
        coords (pd.DataFrame): Expects ['x','y','z'] coordinates columns, followed by a sensortype
            column, and bool ROI columns starting at coords.columns[4].
        roiColors (dict[str, str], optional): Colors to plot ROI electrodes. Must have an entry for
            each ROI in roiNames. Defaults to {ROICOLORS}.
        axesWidth (_type_, optional): Width of one plot. Used to properly scale plot markers.
            Defaults to 8.5/4.
        noROIColor (str, optional): Color for electrodes with no ROI assignment. Defaults to black.
        displaySide (Literal["r", "l"], optional): Hemisphere side to unify to and display.
            Defaults to right side.

    Returns:
        np.ndarray: Array of axes containing the plots, lateral then medial.
    """
    roiNames = coords.columns[4:].to_list()
    # Define sides
    coords = coords.assign(side="")
    coords.loc[coords.x < -20, "side"] = "l"
    coords.loc[coords.x >= -20, "side"] = "lmid"
    coords.loc[coords.x >= 0, "side"] = "rmid"
    coords.loc[coords.x >= 20, "side"] = "r"
    assert coords.side.isin(VALID_SIDES).all()
    # Remove out of brain channels
    cleanedDF = coords[~coords.index.isin(OMITCHANNELS)]

    # Clean up channels with more than two ROIs
    cleanedDF.loc[cleanedDF.index.isin(TRIPLECHANS["DMN-Dorsal-FPN"]), "DMN"] = False

    axs = merged_hemispheres_coverage_plot(
        axs, cleanedDF, roiNames, roiColors, axesWidth, noROIColor, displaySide
    )
    axs[1].get_legend().remove()
    axs[0].set_xlim(-0.5 + 50, 799.5 - 50)
    axs[1].set_xlim(-0.5 + 50 * 1.35, 799.5 - 50)
    axs[0].set_ylim(799.5 - 100, -0.5 + 100)
    axs[1].set_ylim(799.5 - 100 * 1.35, -0.5 + 100)

    return axs


def task_performance_plot(
    subjects: Subjects,
    originalAx: plt.Axes,
    catchAx: plt.Axes,
    conditions=["pre", "post", "grey"],
    performancePath=Path("../Analysis/TaskPerformance"),
):
    """
    Generate task performance plots for paper onto provided axs.

    Args:
        subjects (Subjects): Subjects to include.
        originalAx (plt.Axes): Axes for original trials performance plot.
        catchAx (plt.Axes): Axes for catch trials performance plot.
        conditions (list, optional): Conditions to include. Defaults to ["pre", "post", "grey"].
        performancePath (Path, optional): Path to save performance data to disk. Defaults to
            Path("../Analysis/TaskPerformance").
    """
    performancePath.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.top": False,
        "ytick.right": False,
    })
    jmb = {}
    accData = pd.DataFrame(index=subjects, columns=conditions, dtype="float")
    catchData = pd.DataFrame(index=subjects, columns=conditions, dtype="float")
    allBehav: pd.DataFrame = pd.read_pickle(BEHAVPATH)
    for subjID in subjects:
        jmb[subjID] = allBehav.loc[subjID]
        for cond in conditions:
            accData.loc[subjID, cond] = np.sum(
                ~jmb[subjID].Catch & jmb[subjID][cond] & jmb[subjID].OverrideACC
            ) / np.sum(~jmb[subjID].Catch & jmb[subjID][cond])
            catchData.loc[subjID, cond] = np.sum(
                jmb[subjID].Catch & jmb[subjID][cond] & jmb[subjID].OverrideACC
            ) / np.sum(jmb[subjID].Catch & jmb[subjID][cond])
    numPatients = np.sum(~np.isnan(accData[conditions[0]]))
    markers = ["o" for _ in subjects]
    colors = [CONDITIONCOLORS[cond] for cond in conditions]

    fname = f"Behav_OriginalAcc_{subjects}"
    accData.to_csv(performancePath / f"{fname}.csv")
    plt.sca(originalAx)
    scatterBox(accData, markers, columnColors=colors, linestyle=":")
    originalAx.plot([1, 2, 3], accData.mean(), "k", linewidth=3, marker="o")
    originalAx.grid(False)
    originalAx.set_title(f"Original Trials, N={numPatients} Patients", fontsize=26)

    fname = f"Behav_CatchAcc_{subjects}"
    catchData.to_csv(performancePath / f"{fname}.csv")
    plt.sca(catchAx)
    scatterBox(catchData, markers, columnColors=colors, linestyle=":")
    catchAx.plot([1, 2, 3], catchData.mean(), "k", linewidth=3, marker=None)
    catchAx.grid(False)
    catchAx.set_title(f"Catch Trials, N={numPatients} Patients", fontsize=26)
