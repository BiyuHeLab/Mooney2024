# Model-based RSA searchlight analysis (Fig. 3d)

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
import rsatoolbox.data as rsd # abbreviation to deal with dataset
import rsatoolbox.rdm as rsr

from glob import glob
from rsatoolbox.util.searchlight import get_volume_searchlight, get_searchlight_RDMs, evaluate_models_searchlight
import matplotlib.colors
import os
from glob import glob
from nilearn import image
from nilearn.plotting import plot_img, plot_stat_map, show
import re
import sys
import pickle
from joblib import Parallel, delayed

#Since whole brain searchlight maps have already been generated, we don't need to re-run this.
#Just correlating to the model RDM.

if len(sys.argv) != 2:
    print("Usage: python your_script.py <subjnum>")
    sys.exit(1)

subjnum = sys.argv[1]
print(f"Processing subject: {subjnum}")

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

maskFile = nib.load('/data/' + subjnum + '/Processed/highresMask.nii.gz')
x, y, z = maskFile.get_fdata().shape
affine2save = maskFile.affine

#SL_RDM = get_searchlight_RDMs(betamap, centers, neighbors, sessions,method='crossnobis',n_jobs=1)
filename = f'data/SL_rdms12mm_subject_{subjnum}.pickle'
with open(filename, 'rb') as file:
    SL_RDM = pickle.load(file)

allrdm = SL_RDM.dissimilarities #dissimilarities is just the lower triangle of the matrix.

#Load model:
ind = 10  # number of images
numManip = 7
modelRDM = np.ones((numManip*ind, numManip*ind))

combined_matrix = [
    [0, 0, 0.5, 0.5, 0.5, 0.5, 0.5],
    [0, 0, 0.5, 0.5, 0.5, 0.5, 0.5],
    [0.5, 0.5, 0, 0.5, 0.5, 0.5, 0.5],
    [0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5],
    [0.5, 0.5, 0.5, 0.5, 0, 0.5, 0.5],
    [0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5],
    [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0]
]

for i in range(1, ind+1):
    for j in range(1, ind+1):
        if i == j:
            modelRDM[(i-1)*numManip:i*numManip, (j-1)*numManip:j*numManip] = combined_matrix

an_model = ModelFixed('model RDM', modelRDM)

eval_results = evaluate_models_searchlight(SL_RDM, an_model, eval_fixed, method='rho-a', n_jobs=32)

# get the evaulation score for each voxel
# We only have one model, but evaluations returns a list. By using float we just grab the value within that list
eval_score = [float(e.evaluations) for e in eval_results]

# Create an 3D array, with the size of mask, and 
RDM_brain = np.zeros([x*y*z])
RDM_brain[list(SL_RDM.rdm_descriptors['voxel_index'])] = eval_score
RDM_brain = RDM_brain.reshape([x, y, z])

#Convert eval scores to nifti for later use
#eval_array = np.array(eval_score, dtype = np.float32)
# to save this 3D (ndarry) numpy use this
ni_img = nib.Nifti1Image(RDM_brain, affine2save)

#This is not included in the uploaded code/data, since the files are too big.
os.chdir('/Analysis/fMRI_analysis/Fig3_neural_dist_analysis/Searchlight_analysis_scripts/RDM_data/')
nib.save(ni_img, f'subject12mm_{subjnum}_output.nii.gz')
