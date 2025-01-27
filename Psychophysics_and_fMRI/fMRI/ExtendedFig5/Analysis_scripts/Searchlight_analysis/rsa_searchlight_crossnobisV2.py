import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import nibabel as nib
import seaborn as sns
from nilearn import plotting
from nilearn.maskers import NiftiMasker
from nilearn.image import get_data, load_img, new_img_like
from nltools.stats import fdr, threshold, fisher_r_to_z, one_sample_permutation

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
import re
import sys

if len(sys.argv) != 2:
    print("Usage: python your_script.py <subjnum>")
    sys.exit(1)

subjnum = sys.argv[1]
print(f"Processing subject: {subjnum}")

# updated for crossnobis implementation
# plotted value is averaged lower triangle for every RDM; no model comparison

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

totalimages = 10

# for EV list, pick every 10 images.
# Generating a sequence
manipNumber = 7
evlist = list(range(1, 71))  # 70 total images

# Creating ev-specific sequences:
imagenum = []
for i in range(totalimages):
    imagenum.append(evlist[i : manipNumber * totalimages : totalimages])


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


# Loading a mask:
maskFile = nib.load(
    "fMRI_data/" + subjnum + "/highresMask.nii.gz"
)

x, y, z = maskFile.get_fdata().shape
affine2save = maskFile.affine

# loop over all images
rownum = manipNumber * len(runname)
totalcope = manipNumber
totalrun = len(runname)

# Get searchlight centers and neighbors:

mask = maskFile.get_fdata()
centers, neighbors = get_volume_searchlight(mask, radius=6, threshold=0.1)
byimg_avgtril = []

for i in range(totalimages):
    imgList = imagenum[i]
    data = np.zeros((rownum, x, y, z))
    counter = -1
    for c in range(manipNumber):
        copenum = int(imgList[c])
        for r in range(totalrun):
            counter = counter + 1
            currun = runname[r]
            # file_loc = os.path.join('Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/' + subjnum + '/Processed/Session1/' + currun + '/1stLevel_GLM_' + currun + '.feat/stats/tstat' + str(copenum) + '_2highres.nii.gz')
            file_loc = os.path.join(
                "/fMRI_Data/"
                + subjnum
                + currun
                + "/tstat"
                + str(copenum)
                + "_2highres.nii.gz"
            )
            data[counter] = nib.load(file_loc).get_fdata()

    betamap = data.reshape(rownum, x * y * z)
    # Get RDMs
    n_rep = totalrun
    nCond = totalcope

    conds = np.array(["cond_%02d" % x for x in np.arange(nCond)])
    sessions = np.tile(np.arange(n_rep), totalcope)

    SL_RDM = get_searchlight_RDMs(
        betamap, centers, neighbors, sessions, method="crossnobis", n_jobs=8
    )

    allrdm = SL_RDM.dissimilarities  # dissimilarities is just the lower triangle of the matrix.
    avg_tril = np.mean(allrdm, axis=1, keepdims=True)
    byimg_avgtril.append(avg_tril)

stacked_arrays = np.stack(byimg_avgtril, axis=-1)
avg_tril_final = np.mean(stacked_arrays, axis=-1)

# Create an 3D array, with the size of mask, and
x, y, z = mask.shape
RDM_brain = np.zeros([x * y * z])

# For lower triangle average
RDM_brain[list(SL_RDM.rdm_descriptors["voxel_index"])] = avg_tril_final.flatten()
RDM_brain2 = RDM_brain.reshape([x, y, z])

#Plot the voxels above the 99th percentile
# threshold = np.percentile(eval_score, 95)
#os.chdir("invar_bp_fmri/")
os.chdir("/fMRI_Data/")

#Save output:
tmp_img = maskFile
plot_img = new_img_like(tmp_img, RDM_brain2)
sl_outputfile = f"{subjnum}_wholebrain.nii.gz"
plot_img.to_filename(sl_outputfile)
nib.save(plot_img, sl_outputfile)  # save nifti file.
