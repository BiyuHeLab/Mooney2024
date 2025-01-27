#Model-based RSA script, using c.v. euclidean distance RDM outputs (Fig. 3)
#Model is based on behavior results:
#Full invariance for size manipulation, within the same image
#Parital invariance for position & rotation, within the same image
#No invariance (specificity) for cross-image distances

import numpy as np
from scipy import io
import matplotlib.pyplot as plt
import rsatoolbox
from itertools import chain, repeat
import os
from scipy.io import savemat
import pickle
from scipy.io import savemat
import os

directory_path = f'~/results/'
if not os.path.exists(directory_path):
    os.makedirs(directory_path)

allsubs=["001","004","005","006","007","008","009","010","012","013"]
subjnum = len(allsubs)

roiall = [
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

roinum = len(roiall)

#Load model:
ind = 10  # number of images
numManip = 7 # number of manipulations

combined_matrix = [
        [0, 0, 0.5, 0.5, 0.5, 0.5, 0.5],
        [0, 0, 0.5, 0.5, 0.5, 0.5, 0.5],
        [0.5, 0.5, 0, 0.5, 0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5, 0.5, 0, 0.5, 0.5],
        [0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5],
        [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0]

np_matrix = np.array(combined_matrix)
result_matrix = np_matrix 
modelRDM = np.ones((numManip*ind, numManip*ind)) #specificity: distance = 1

for i in range(1, ind+1):
    for j in range(1, ind+1):
        if i == j:
            modelRDM[(i-1)*numManip:i*numManip, (j-1)*numManip:j*numManip] = result_matrix

lower_triangular_indices = np.tril_indices(numManip*ind, k=-1)

modelrdm_array = modelRDM[lower_triangular_indices]

#Load subject data:

for r in range(roinum):
    roi=[roiall[r]]
    #print(roi)
    model = rsatoolbox.model.ModelFixed(roi[0], modelrdm_array)
    fmri_rdms=[]
    rdms_array = [] #create empty dictionary to store matrices
    for s in range(subjnum):
        m = allsubs[s];
        filepath = '/Analysis/fMRI_analysis/Fig3_neural_dist_analysis/ROI_analysis_scripts/data/{}/rdms/{}_{}_crossnobisBYIMG_pseudo_rdm.mat'.format(m,m,roi[0])
        rdm_matrix = io.matlab.loadmat(filepath)['rdm'] #rdm is variable name in the .mat file
        rdm_matrix = np.squeeze(rdm_matrix)  # Remove singleton dimensions
        rdm_matrix = rdm_matrix.reshape((numManip*ind, numManip*ind))
        rdms_array.append(rdm_matrix[lower_triangular_indices])

    n_models = subjnum
    rdmarray = np.array([rdms_array[i] for i in range(n_models)]) #RDM data

    #n_models = len(rdmarray)
    roi_labels = list(repeat(roi[0], subjnum))
    subject_label = allsubs

    fmri_rdms = rsatoolbox.rdm.RDMs(rdmarray,
                                rdm_descriptors={'ROIs':roi_labels,
                                                 'SubjectNum':subject_label},
                                dissimilarity_measure='crossnobis'
                               )

    results_1 = rsatoolbox.inference.eval_fixed(model, fmri_rdms, method='rho-a')
    print(results_1)

    os.chdir(directory_path)

    with open(f'modelcomparison_{roi[0]}.pickle', 'wb') as file:
        pickle.dump(results_1, file)

os.chdir(directory_path)
for r in range(roinum):
    roi=[roiall[r]]    
    with open(f'modelcomparison_{roi[0]}.pickle', 'rb') as file:
        file_name = os.path.join(directory_path,f'{roi[0]}_modelcomparison.mat')
        result = pickle.load(file)
        mat_data = {
    'evaluations': result.evaluations,
    'noise_ceiling': result.noise_ceiling,
    'variances': result.variances if result.variances is not None else np.array([]),
    'degrees_of_freedom': np.array([result.dof]),
    'n_rdm': np.array([result.n_rdm]),
    'n_pattern': np.array([result.n_pattern if result.n_pattern is not None else 0])
}
        savemat(file_name, mat_data)


