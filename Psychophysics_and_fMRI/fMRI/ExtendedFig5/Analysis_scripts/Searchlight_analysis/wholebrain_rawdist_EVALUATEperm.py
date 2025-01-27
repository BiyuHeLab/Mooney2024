import pickle
import numpy as np
import random
import nibabel as nib
from rsatoolbox.util.rdm_utils import batch_to_matrices
import os
from joblib import Parallel, delayed
from statistics import NormalDist
import sys
from scipy.stats import norm
import scipy.stats as st

#Permutation test for neural distance analysis in searchlight

allsubs = ["001", "004", "005", "006", "007", "008", "009", "010", "012", "013"]
nReps = 1000  # nperm is 1000, as run in ~Scripts/wholebrain_rawdist_perm.py


for s in allsubs:
    maskFile = nib.load(
        "fMRI_Data/"
        + str(s)
        + "/highresMask.nii.gz"
    )
    x, y, z = maskFile.get_fdata().shape
    affine2save = maskFile.affine
    mask = maskFile.get_fdata()
    print(s)
    file_loc = os.path.join(
        "fMRI_Data/Searchlight_T1/"
        + str(s)
        + "_SL2.nii.gz"
    )
    real = nib.load(file_loc).get_fdata()
    perm = np.zeros((nReps, x, y, z))

    for p in range(nReps):
        perm_loc = os.path.join(
            "/All_searchlight/"
            + str(s)
            + "/perm/"
            + str(s)
            + "_perm"
            + str(p)
            + ".nii.gz"
        )
        permdata = nib.load(perm_loc).get_fdata()
        perm[p, :, :, :] = permdata
    # print(perm.shape)

    totalvox = len(real)
    zmap = np.zeros((x, y, z))
    pmap = np.zeros((x, y, z))
    empirical_zmap = np.zeros((x, y, z))
    for i in range(x):
        for j in range(y):
            for k in range(z):
                if mask[i, j, k] > 0:
                    permvector = perm[:, i, j, k]
                    # print(len(permvector))
                    zmap[i, j, k] = (real[i, j, k] - np.nanmedian(permvector)) / np.nanstd(
                        permvector
                    )
                    pmap[i, j, k] = (np.sum(permvector <= real[i, j, k]) + 1) / (nReps + 1)
                    empirical_zmap[i, j, k] = st.norm.ppf(pmap[i, j, k] + 1e-10)
                else:
                    zmap[i, j, k] = 0
                    pmap[i, j, k] = 0
                    empirical_zmap[i, j, k] = 0

    zmap[np.isnan(zmap)] = 0
    empirical_zmap = np.nan_to_num(empirical_zmap)
    empz_img = nib.Nifti1Image(empirical_zmap, affine2save)
    ni_img = nib.Nifti1Image(zmap, affine2save)
    pmap_img = nib.Nifti1Image(pmap, affine2save)

    pathdir = (
        "fMRI_Data/Searchlight_rawdist_zmaps"
    )
    os.chdir(pathdir)
    nib.save(empz_img, f"{s}_EMPzmap.nii.gz")
    nib.save(ni_img, f"{s}_zmap.nii.gz")
    nib.save(pmap_img, f"{s}_pmap.nii.gz")
