import logging
from pathlib import Path

import mne
import numpy as np
import pandas as pd

from pwmv_ecog.process import ECoGHGB, ECoGProcessing

# Setup logging
logger = logging.getLogger(__name__)
logger.propagate = False
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler("PWMV_ECoG_ECoGData.log")
console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.INFO)
# Will use the default formatter for the console, so only define one for file
file_format = logging.Formatter(
    "{asctime}\t{name}\t{levelname}: {message}", style="{", datefmt="%Y%m%d %H:%M:%S"
)
file_handler.setFormatter(file_format)
logger.addHandler(console_handler)
logger.addHandler(file_handler)


class ECoGDataMNE:
    """
    Canonical data structure for PWMV analysis generated from mne.Epochs.

    Attributes:
        data (pd.DataFrame): index=pd.MultiIndex(['Electrode','TrialNum']), columns=epochs.times.
            Taken from passed epochs, so subjected to whatever automatic modifications set there,
            such as baseline correction. Also, each trial's data is set to NaN beyond image offset.
        behav (pd.DataFrame): Passed behavioral data.
    """

    def __init__(self, subjectID: str, epochs: mne.Epochs, behav: pd.DataFrame):
        """
        Args:
            subjectID (str): Patient number.
            epochs (mne.Epochs): Trial data.
            behav (pd.DataFrame): should be like the output from
                pwmv_ecog.behavior.BehavioralAnalysis.get_mooney.
        """
        self.subjectID = subjectID
        self._add_behav(behav)
        self.numChannels = len(epochs.ch_names)
        self.numTimepoints = len(epochs.times)
        self.fsample = epochs.info["sfreq"]
        self._init_data(epochs)
        _, self.numTrials = self.data.index.levshape
        self._init_edf2mni(epochs)

    def _init_data(self, epochs: mne.Epochs) -> None:
        if epochs.metadata is None:
            raise ValueError("Passed epochs must have metadata.")
        epochs.load_data()
        # ImagePresentedDuration in ms, but time_as_index expects seconds
        epochs.metadata["ImageOffsetIdx"] = epochs.time_as_index(
            epochs.metadata["ImagePresentedDuration"] / 1000
        )
        allData = epochs.get_data()
        for i_epoch, trialID in enumerate(epochs.metadata.index):
            # Remove data after image offset (set to NaN)
            allData[i_epoch, :, epochs.metadata.loc[trialID, "ImageOffsetIdx"] + 1 :] = np.nan
        self._allData = pd.DataFrame(
            allData.reshape(allData.shape[0] * allData.shape[1], allData.shape[2]),
            index=pd.MultiIndex.from_product(
                [epochs.metadata.MooneyTrialNum, epochs.ch_names], names=["TrialNum", "Electrode"]
            ),
            columns=epochs.times,
        ).reorder_levels([1, 0])
        selectedChannels = epochs.copy().pick(["seeg", "ecog"]).ch_names
        data = self._allData.loc[selectedChannels, :]
        # Ensure removed channels are not lingering in the index
        self.data = data.set_index(data.index.remove_unused_levels())

    def _init_edf2mni(self, epochs: mne.Epochs) -> None:
        self.edf2mni = {ch: ch for ch in epochs.ch_names}

    def _add_behav(self, behav: pd.DataFrame) -> None:
        # Validate behavioral DataFrame and set self.behav
        BEHAVREQS = [
            "Subject",
            "ImageNum",
            "ImageType",
            "Responded",
            "ACC",
            "OverrideACC",
            "Catch",
            "pre",
            "post",
            "grey",
            "disambig",
        ]
        memberships = np.isin(BEHAVREQS, list(behav))
        if all(memberships):
            self.behav = behav
        else:
            errormsg = f"behav missing mandatory columns: {np.array(BEHAVREQS)[~memberships]}"
            raise ValueError(errormsg)


def load_mne_subject(subjectID: str, loadGamma=False) -> ECoGDataMNE:
    """
    Load ECoGData for one subject from the standard ECoGProcessing object.

    Args:
        subjectID (str): Integer identifier for the subject.

    Returns:
        ECoGData: Data & behavior for subjectID.
    """
    if loadGamma:
        proc: ECoGProcessing = ECoGHGB(local=False)
    else:
        proc = ECoGProcessing(local=False)
    proc.load(f"Patient_{subjectID}")
    allBehav = pd.read_pickle(Path("pwmv_ecog", "manuscript", "manuscript_behav.pd.pkl"))
    behav = allBehav.loc[int(subjectID)]
    # Trim behav of any trials eliminated in epochs processing
    behav = behav.loc[behav.index.isin(proc.epochs.metadata["MooneyTrialNum"])]
    return ECoGDataMNE(subjectID, proc.epochs, behav)
