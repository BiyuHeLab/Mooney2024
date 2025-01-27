# This script is for calculating within-image distances (Supplementary Fig. 5)

import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from nilearn import plotting
from nilearn.maskers import NiftiMasker
from nilearn.image import get_data, load_img, new_img_like

import rsatoolbox
from rsatoolbox.inference import eval_fixed
from rsatoolbox.model import ModelFixed
from rsatoolbox.rdm import RDMs
import rsatoolbox.data as rsd  # abbreviation to deal with dataset
import rsatoolbox.rdm as rsr

from glob import glob
from rsatoolbox.util.searchlight import (
    get_volume_searchlight,
    get_searchlight_RDMs,
    evaluate_models_searchlight,
)
import matplotlib.colors
import os
from glob import glob
from nilearn import image
from nilearn.plotting import plot_img, plot_stat_map, show
import pickle
import json
from scipy.io import savemat
import re
import sys

# source Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/rsatoolbox-main/.venv/bin/activate

# Load subject:
if len(sys.argv) != 2:
    print("Usage: python your_script.py <subjnum>")
    sys.exit(1)

subjnum = sys.argv[1]
print(f"Processing subject: {subjnum}")

# Regressor titles:
evTitles = [
    "Original1",
    "CW1",
    "CCW1",
    "Rshift1",
    "Lshift1",
    "LRinvert1",
    "Size1",
    "Original5",
    "CW5",
    "CCW5",
    "Rshift5",
    "Lshift5",
    "LRinvert5",
    "Size5",
    "Original6",
    "CW6",
    "CCW6",
    "Rshift6",
    "Lshift6",
    "LRinvert6",
    "Size6",
    "Original7",
    "CW7",
    "CCW7",
    "Rshift7",
    "Lshift7",
    "LRinvert7",
    "Size7",
    "Original9",
    "CW9",
    "CCW9",
    "Rshift9",
    "Lshift9",
    "LRinvert9",
    "Size9",
    "Original2",
    "CW2",
    "CCW2",
    "Rshift2",
    "Lshift2",
    "LRinvert2",
    "Size2",
    "Original3",
    "CW3",
    "CCW3",
    "Rshift3",
    "Lshift3",
    "LRinvert3",
    "Size3",
    "Original4",
    "CW4",
    "CCW4",
    "Rshift4",
    "Lshift4",
    "LRinvert4",
    "Size4",
    "Original8",
    "CW8",
    "CCW8",
    "Rshift8",
    "Lshift8",
    "LRinvert8",
    "Size8",
    "Original10",
    "CW10",
    "CCW10",
    "Rshift10",
    "Lshift10",
    "LRinvert10",
    "Size10",
]

#Subject 005 is missing the last 4 runs = 2 pseudoruns
if subjnum == "005": 
    runname = [
        "PseudoMain1",
        "PseudoMain2",
        "PseudoMain3",
        "PseudoMain4",
        "PseudoMain5",
        "PseudoMain6",
    ]
else:
    runname = [
        "PseudoMain1",
        "PseudoMain2",
        "PseudoMain3",
        "PseudoMain4",
        "PseudoMain5",
        "PseudoMain6",
        "PseudoMain7",
        "PseudoMain8",
    ]

roinames = [
    "new_v1_2highres",
    "new_v2_2highres",
    "new_v3_2highres",
    "new_hv4_2highres",
    "new_LO1_2highres",
    "new_LO2_2highres",
    "FG_3percent_2highres",
    "new_IPS0_2highres",
    "new_IPS1_2highres",
    "new_IPS2_2highres",
    "new_IPS3_2highres",
    "new_IPS4_2highres",
    "new_IPS5_2highres",
    "new_SPL1_2highres",
    "FPN_F_2highres",
    "FPN_P_2highres",
    "DMN_P_2highres",
    "DMN_PCC_2highres",
    "DMN_MPFC_2highres"
]
roi_niftis = [roi + ".nii.gz" for roi in roinames]

def upper_tri(RDM):
    """upper_tri returns the upper triangular index of an RDM

    Args:
        RDM 2Darray: squareform RDM

    Returns:
        1D array: upper triangular vector of the RDM
    """
    # returns the upper triangle
    m = RDM.shape[0]
    r, c = np.triu_indices(m, 1)
    return RDM[r, c]


def RDMcolormapObject(direction=1):
    """
    Returns a matplotlib color map object for RSA and brain plotting
    """
    if direction == 0:
        cs = ["yellow", "red", "gray", "turquoise", "blue"]
    elif direction == 1:
        cs = ["blue", "turquoise", "gray", "red", "yellow"]
    else:
        raise ValueError("Direction needs to be 0 or 1")
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", cs)
    return cmap


#Loading ROIs:
roiFilePath = (
    "Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/" + subjnum + "/T1_rois/"
)
all_roi_paths = glob(f"{roiFilePath}/*.nii.gz")  # all nifti files
roi_paths = [path for path in all_roi_paths if path.split("/")[-1] in roi_niftis]  # selected ROIs
roi_paths.sort()

#Load one image to get the dimensions and make the mask
tmp_roi = nib.load(roi_paths[0])
x, y, z = tmp_roi.get_fdata().shape

#Loop over all roi masks
roimasks = np.zeros((len(roi_paths), x, y, z))
for i in range(len(roinames)):
    roiname2find = roinames[i]
    pattern = r"\b{}\b".format(re.escape(roiname2find))
    thisROI_path = [path for path in roi_paths if re.search(pattern, path)]
    print(thisROI_path)
    roimasks[i] = nib.load(thisROI_path[0]).get_fdata()

rownum = len(evTitles) * len(runname)
totalcope = len(evTitles)
totalrun = len(runname)

#Rearranging by image, not manipulation: use [allevbyimg] as indices.
allev = np.array(range(70))
allevbymanip = allev.reshape(7, 10)
allevbyimg_temp = []
for r in range(10):
    allevbyimg_temp.append(allevbymanip[:, r])
allevbyimg = np.array(np.concatenate(allevbyimg_temp))

#Load data
noise_loc = []
data = np.zeros((rownum, x, y, z))
counter = -1
for c in range(totalcope):
    copeind = int(c)
    copenum = allevbyimg[copeind]
    # copenum = copeind #if analyzing by manipulation
    for r in range(totalrun):
        counter = counter + 1
        currun = runname[r]
        file_loc = os.path.join(
            "/fMRI_Data/"
            + subjnum
            + currun
            + "/tstat"
            + str(copenum + 1)
            + "_2highres.nii.gz"
        )
        data[counter] = nib.load(file_loc).get_fdata()

betamap = data.reshape(rownum, x * y * z)

n_rep = totalrun
nCond = totalcope

conds = np.array(["cond_%02d" % x for x in np.arange(nCond)])
sessions = np.tile(np.arange(n_rep), totalcope)
print(sessions)

conds = np.repeat(conds, n_rep)
obs_des = {"conds": conds, "sessions": sessions}
des = {"subj": 1}

# Calculate RDM for each ROI.
# Do not pass a noise covariance matrix; calculating c.v. Euclidean distance
subset_data_list = []
for r in range(len(roinames)):
    chn_des = roimasks[r].reshape(-1)  # stretch out into vector form
    roifilename = roinames[r]
    chn_des[chn_des > 0] = 1
    chn_des = {"ROIs": np.array(chn_des)}
    dataset = rsd.Dataset(
        measurements=betamap, obs_descriptors=obs_des, channel_descriptors=chn_des
    )
    subset_temp = dataset.subset_channel(by="ROIs", value=[1])
    
    rdm_cv = rsatoolbox.rdm.calc_rdm(
        subset_temp, method="crossnobis", descriptor="conds", cv_descriptor="sessions"
    )
    rdm2save = rdm_cv.get_matrices()
    file_path = f"fMRI_Data/Results/{subjnum}/rdms/{subjnum}_{roifilename}_crossnobisBYIMG_pseudo_rdm.mat"
    savemat(file_path, {"rdm": rdm2save})
