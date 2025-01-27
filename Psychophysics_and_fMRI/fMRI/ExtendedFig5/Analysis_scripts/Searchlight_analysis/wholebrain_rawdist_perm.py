import pickle
import numpy as np
import random
import nibabel as nib
from rsatoolbox.util.rdm_utils import batch_to_matrices
import os
from joblib import Parallel, delayed
import sys

if len(sys.argv) != 2:
    print("Usage: python your_script.py <subjnum>")
    sys.exit(1)

subjnum = sys.argv[1]
print(f"Processing subject: {subjnum}")

#Load mask:
maskFile = nib.load(
    "/fMRI_Data/"
    + subjnum
    + "/highresMask.nii.gz"
)
x, y, z = maskFile.get_fdata().shape
affine2save = maskFile.affine
mask = maskFile.get_fdata()

wholebrainRDM_path = (
    "~/Desktop/Ayaka/BigPurpleCluster/70by70RDM_wholebrain/SL_rdms12mm_subject_"
    + subjnum
    + ".pickle"
)

with open(wholebrainRDM_path, "rb") as file:
    wholebrain = pickle.load(file)

numManip = 7
ind = numManip
numlowertri = np.array(numManip * (numManip - 1) / 2)
imgbyimgMAT = []
counter = -1
submatrix = []
totalvox = len(wholebrain)
permutedMatrixall = np.zeros((totalvox, 45, 21))

for r in range(totalvox):
    rdmbyvoxel = wholebrain[r].get_matrices()
    rdmbyvoxel = np.array(rdmbyvoxel[0])
    rdmbyvoxel = rdmbyvoxel.reshape(70, 70)
    permutedMatrix = []  # created a new cross-image matrix for every voxel cluster
    for i in range(1, 11):
        for j in range(1, 11):
            submatrix = rdmbyvoxel[(i - 1) * ind : i * ind, (j - 1) * ind : j * ind]
            if i > j:
                submatrix = submatrix - np.diag(np.diag(submatrix))
                temp = np.tril(submatrix, -1)
                temp2 = temp[temp != 0]
                if len(temp2) != 0:
                    temp2 = temp2.reshape(1, 21)
                    permutedMatrix.append(np.array(temp2))
    if len(permutedMatrix) != 0:
        permutedMatrix = np.vstack(permutedMatrix)
        permutedMatrixall[r, :, :] = permutedMatrix


def compute_perm(subjnum, permutedMatrixall, totalvox, p):
    perm_all = np.zeros((totalvox, 1))
    for r in range(totalvox):
        perm2 = permutedMatrixall[r, :, :]
        pairing = []
        for d in range(21):
            RANDimgselection = random.sample(range(45), 10)
            pairing.append(perm2[np.array(RANDimgselection), d])
        perm_all[r] = np.mean(pairing)
        # 10 images x 21 lower triangle vals
    x, y, z = mask.shape
    RDM_brain = np.zeros([x * y * z])
    perm_all = perm_all.reshape((totalvox,))
    RDM_brain[list(wholebrain.rdm_descriptors["voxel_index"])] = perm_all
    RDM_brain = RDM_brain.reshape([x, y, z])
    ni_img = nib.Nifti1Image(RDM_brain, affine2save)
    pathdir = (
        "All_searchlight/"
        + subjnum
        + "/perm/"
    )
    if not os.path.exists(pathdir):
        os.makedirs(pathdir)
    os.chdir(pathdir)
    nib.save(ni_img, f"{subjnum}_perm{p}.nii.gz")


num_iterations = 1000
results = Parallel(n_jobs=8)(
    delayed(compute_perm)(subjnum, permutedMatrixall, totalvox, p) for p in range(num_iterations)
)
