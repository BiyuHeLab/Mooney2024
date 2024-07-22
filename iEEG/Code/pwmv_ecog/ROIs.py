import logging
import os
from pathlib import Path
from typing import Any, Iterator

import nibabel as nib
import numpy as np
import numpy.linalg as npl
import pandas as pd
from nibabel.affines import apply_affine
from nibabel.nifti1 import Nifti1Image
from tqdm import tqdm

from pwmv_ecog.channels import (
    get_exclude_channels_from_list,
    get_patient_invalid_channels,
    get_patient_rejected_channels,
    standardize_ch_name,
)
from pwmv_ecog.paths import ECoGPaths, FilePath

logger = logging.getLogger(__name__)


def get_nifti_masks(maskFolder: Path) -> dict[str, nib.nifti1.Nifti1Image]:
    """
    Load all niftis in maskFolder.

    Args:
        maskFolder (Path): Path to the folder.

    Returns:
        dict[str, nib.nifti1.Nifti1Image]: {base filename: nifti image object}
    """
    maskFiles = [
        os.path.join(maskFolder, f)
        for f in os.listdir(maskFolder)
        if (f.endswith(".nii")) or (f.endswith(".nii.gz"))
    ]
    loaded = {os.path.basename(fname).split(".")[0]: nib.load(fname) for fname in maskFiles}
    niftis = {name: img for name, img in loaded.items() if isinstance(img, nib.nifti1.Nifti1Image)}
    nonNiftis = set(loaded.keys()) - set(niftis.keys())
    if nonNiftis:
        logger.error(f"Excluding non-nifti files found in mask folder: {nonNiftis}")
    return niftis


def bipolar_remap(coords: pd.DataFrame, patientNum: int, drop=False, local=False) -> pd.DataFrame:
    """
    Update the MNI coords dataframe with any bipolar referencing done for the patient.

    Always standardizes coords.index.
    There is no check that coords is for the same patient.
    If drop=True, removes the channels used to produce the new virtual channels. Channels that
    are not in coords are silently ignored, but the new virtual channels made from them are
    still added to the updated coords.

    Args:
        coords (pd.DataFrame): Dataframe read in from f"Patient_{patientNum}_MNI_coords.txt").
            Expected index=channel names, and columns=['x','y','z','chanType']
        patientNum (int): Patient ID number.
        drop (bool): Whether to remove the depth electrodes used to create the virtual
            electrodes from the returned coords.

    Returns:
        pd.DataFrame: Updated coords.
    """
    # Remove the dependency on ECoGProcessing as it creates a circular loop
    proc = ECoGPaths(local=local)
    remapPath = proc.bipolarRemapPath(f"Patient_{patientNum}")
    coords = coords.set_index(map(standardize_ch_name, coords.index))
    if os.path.isfile(remapPath):
        logger.info(f"Reading bipolar remapping details from {remapPath}.")
        remap = pd.read_table(remapPath, sep=",", index_col=0)
        meanCoords = lambda row: coords.loc[row[["anode", "cathode"]], ["x", "y", "z"]].mean()
        # remapCoords will have remap's index and columns='x','y','z'
        remapCoords = remap.apply(meanCoords, axis=1)
        remapCoords["sensortype"] = "D"
        if drop:
            # Drop the rereferenced channels, append the new virtual channels
            remappedChans = set(remap.anode).union(remap.cathode).intersection(coords.index)
            return pd.concat([coords.drop(list(remappedChans)), remapCoords])
        else:
            return pd.concat([coords, remapCoords])
        # return pd.concat([coords.drop(remappedChans), remapCoords])
    else:
        logger.warning("No bipolar mapping CSV found. Nothing to do.")
        return coords


def get_patient_coords(patientDir: Path, dropBad=True, dropNonBipolar=True) -> pd.DataFrame:
    """
    Get MNI coordinates from a patient folder. Optionally restrict only to good channels.

    Args:
        patientDir (Path): Folder where Patient_#_MNI_coords.txt is. Final component of path
            expected to be of the form Patient_# so the patient number can be extracted.
        dropBad (bool, optional): If True, bad channels are dropped, per processed subject yaml.
            Defaults to True.
        dropNonBipolar (bool, optional): When bipolar referencing is used to generate a virtual
            electrode, drop the input electrodes from the final results. Defaults to True.
    Returns:
        pd.DataFrame: Index=channel names. Columns=['x','y','z','sensortype']
    """
    patientNum = int(patientDir.name.split("_")[-1])
    mniCoordsFile = Path(patientDir, f"Patient_{patientNum}_MNI_coords.txt")
    coords = (
        pd.read_table(mniCoordsFile, sep=" ", header=None, index_col=0, usecols=range(5))
        .rename_axis(index={0: "Channel"})
        .rename(columns={1: "x", 2: "y", 3: "z", 4: "sensortype"})
    )
    coords = bipolar_remap(coords, patientNum, drop=dropNonBipolar)
    coords.drop(
        index=get_exclude_channels_from_list(list(coords.index), patientNum=patientNum),
        inplace=True,
    )
    if dropBad:
        coords.drop(
            index=get_patient_invalid_channels(patientNum)
            + get_patient_rejected_channels(patientNum),
            inplace=True,
        )
    return coords


class pwmvROIs:
    """
    Manages access to a group of nifti mask defined ROIs.

    Will read any full_names.csv in the roiMasksFolder for full print names for each ROI. Create
    such a CSV from a dict like consolidatedROINames with the following code:
    pd.DataFrame(consolidatedROINames, index=[0]).T.to_csv('full_names.csv', header=False)
    """

    def __init__(self, roiMasksFolder: FilePath = "Electrodes/ROI_Masks/network"):
        self.roiMasksFolder = Path(roiMasksFolder)
        self.allROIs: dict[str, Nifti1Image] = get_nifti_masks(Path(self.roiMasksFolder))
        if (fullNamesCSV := self.roiMasksFolder / "full_names.csv").is_file():
            self.fullNames: dict[str, str] = (
                pd.read_csv(fullNamesCSV, header=None, index_col=0).squeeze().to_dict()
            )
        else:
            self.fullNames = {k: k for k in self.keys()}

    @property
    def invAffines(self) -> dict[str, np.ndarray[Any, np.dtype[np.float64]] | None]:
        if not hasattr(self, "_invAffines"):
            self._invAffines: dict[str, np.ndarray[Any, np.dtype[np.float64]] | None] = {
                name: None if img.affine is None else npl.inv(img.affine)
                for name, img in self.allROIs.items()
            }
        return self._invAffines

    def get_patient_coords_and_memberships(
        self, patientNums: list[int], patientsFolder=Path("../rawData")
    ) -> pd.DataFrame:
        """
        Get MNI coordinates and ROI memberships for a group of patients.

        Args:
            patientNums (list[int]): Patient numbers to include.
            patientsFolder (_type_, optional): Folder containing Patient folders. Defaults to
                Path("../rawData").

        Returns:
            pd.DataFrame: Columns=['x','y','z','sensortype', *ROI names], where ROI names are
                self.keys().
        """
        # Get coords
        channelCoords: list[pd.DataFrame] = []
        loadedPNums = []
        for pNum in patientNums:
            try:
                channelCoords.append(get_patient_coords(Path(patientsFolder, f"Patient_{pNum}")))
                loadedPNums.append(pNum)
            except Exception as e:
                print(f"FAILED: Could not get Patient_{pNum} coordinates.")
                raise e
        coords = pd.concat(channelCoords, keys=loadedPNums, names=["PatientNum", "Channel"])
        # Get memberships
        memberships = self.get_patient_memberships(patientNums)
        coordsAndMemberships = pd.concat((coords, memberships), axis=1)
        return coordsAndMemberships

    def get_patient_memberships(
        self, patientNums: list[int], patientsFolder=Path("../rawData")
    ) -> pd.DataFrame:
        """
        Convenience method to get ROI memberships for a group of patients. Returned DataFrame's
        columns are the ROI names with bool dtype.
        """
        channelCoords: list[pd.DataFrame] = []
        for pNum in patientNums:
            channelCoords.append(get_patient_coords(Path(patientsFolder, f"Patient_{pNum}")))
        coords = pd.concat(channelCoords, keys=patientNums, names=["PatientNum", "Channel"])
        return self.memberships(coords, standardize=False).sort_index(axis=1)

    def memberships(self, coords: pd.DataFrame, standardize=True) -> pd.DataFrame:
        """
        Test each coordinates in coords for membership in each ROI.

        Args:
            coords (pd.DataFrame): Index should be channel names if standardize is True. Columns
                must include ['x','y','z'].
            standardize (bool, optional): Standardize channel names in returned DataFrame index.
                Defaults to True.

        Returns:
            pd.DataFrame: Index == coords.index. One column for each ROI.
        """
        channels = list(map(standardize_ch_name, coords.index)) if standardize else coords.index
        rois = {}
        for roiName, image in tqdm(self, unit="Mask", leave=False):
            maskMembership: list[bool] = []
            image_data = image.get_fdata()
            translatedCoords = (
                apply_affine(self.invAffines[roiName], coords[["x", "y", "z"]].values)
                .round()
                .astype("int")
            )
            for elecCoords in translatedCoords:
                try:
                    # Any non-zero value considered in the ROI
                    maskMembership.append(bool(image_data[tuple(elecCoords)]))
                except IndexError:
                    maskMembership.append(False)
            rois[roiName] = pd.Series(maskMembership, index=channels, name=roiName)
        return pd.DataFrame(rois)

    def __str__(self) -> str:
        return str(self.allROIs.keys())

    def __repr__(self) -> str:
        return f'pwmv_ecog.ROIs.pwmvROIs("{self.roiMasksFolder}")'

    def __getitem__(self, key) -> Nifti1Image:
        return self.allROIs[key]

    def __len__(self) -> int:
        return len(self.allROIs)

    def __iter__(self) -> Iterator[tuple[str, Nifti1Image]]:
        return iter(self.allROIs.items())

    def __contains__(self, item) -> bool:
        return item in self.allROIs

    def keys(self) -> list[str]:
        return list(self.allROIs.keys())

    def values(self) -> list[Nifti1Image]:
        return list(self.allROIs.values())

    def items(self) -> list[tuple[str, Nifti1Image]]:
        return list(self.allROIs.items())
