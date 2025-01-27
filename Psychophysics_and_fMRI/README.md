## **Psychophysics**


### *Behavioral Experiment 1 (Fig. 2a-2e)*

Source data:  
`~/Psychophysics/Fig2/source_data/meanbysubj_exp1.mat`  
Data averaged across trials, per subject. Separated by condition (pre/grayscale/post) and manipulation.
meanbysubj_exp1.mat is a 10x30 matrix where each row is the manipulation and each column is a subject.
Note that only 7 out of 10 manipulations are used in the final analysis for Experiment 1.

Analysis:  
`~/Psychophysics/Fig2/Analysis_scripts/mooney_behavior_analysis_experiment1.m`

Plotting:  
`~/Psychophysics/Fig2/plot_Fig2_exp1.m`

Statistics:  
`~/Psychophysics/Fig2/Statistics/Mooney_invar_Exp1_alltrials.jasp`

-----

### *Behavioral Experiment 2 (Fig. 2f-2j)*

Source data:  
`~/Psychophysics/Fig2/source_data/meanbysubj_exp2.mat`  
Data averaged across trials, per subject. Separated by condition (pre/grayscale/post) and manipulation.

Analysis:  
 `~/Psychophysics/Fig2/Analysis_scripts/mooney_behavior_analysis_experiment2.m`

Plotting:  
`~/Psychophysics/Fig2/plot_Fig2_exp2.m`

Statistics:  
`~/Psychophysics/Fig2/Statistics/Mooney_invar_Exp2_alltrials.jasp`

-----

### *Population receptive field (pRF) Analysis (Supplementary Fig. 1)*

Source data:  
`~/Psychophysics/ExtendedFig1/source_data/featureRF.mat`

Analysis:  
`~/Psychophysics/ExtendedFig1/Analysis_scripts/feature_analysis_final.m`

Plotting:  
`~/Psychophysics/ExtendedFig1/plot_ExtendedFig1a_image_prf_overlay.m`  
`~/Psychophysics/ExtendedFig1/plot_ExtendedFig1c.m`

-----

### *Behavioral Control Analysis 1--Excluding gaze-shifted trials for Position manipulation (Supplementary Fig. 2):*

Source data:  
`~/Psychophysics/ExtendedFig2/source_data/meanbysubjVF.mat`  
Data averaged across trials, per subject. Separated by condition (pre/grayscale/post) and manipulation.
meabysubjVF.mat is a 10x30 matrix where each row is the manipulation and each column is a subject.
Note that only 7 out of 10 manipulations are used in the final analysis for Experiment 1.

Analysis:  
`~/Psychophysics/ExtendedFig2/Analysis_scripts/ExtendedFig2_mooney_analysis_eyetrackingVF.m`

Plotting:  
`~/Psychophysics/ExtendedFig2/plot_ExtendedFig2.m`

-----

### *Behavioral Control Analysis 2—Excluding unrecognized grayscale (lapse) trials for Experiments 1&2 (Supplementary Fig. 3):*

Source data:  
`~/Psychophysics/ExtendedFig3/source_data/Gyes_bySubj_pre_exp1.mat`  
`~/Psychophysics/ExtendedFig3/source_data/Gyes_bySubj_post_exp1.mat`  
These are 10x30 matrices where each row is the manipulation and each column is a subject.  
Note that only 7 out of 10 manipulations are used in the final analysis for Experiment 1.

`~/Psychophysics/ExtendedFig3/source_data/Gyes_bySubj_pre_exp2.mat`  
`~/Psychophysics/ExtendedFig3/source_data/Gyes_bySubj_post_exp2.mat`  
These are 6x12 matrices where each row is the manipulation and each column is a subject.

Analysis scripts:  
`~/Psychophysics/ExtendedFig3/Analysis_scripts/Psychophysics_analysis/ExtendedFig3_Gsyesnoanalysis/mooney_analysis_experiment1_GSyesno.m`  
`~/Psychophysics/ExtendedFig3/Analysis_scripts/Psychophysics_analysis/ExtendedFig3_Gsyesnoanalysis/mooney_analysis_experiment2_Gsyesno.m`  

Statistics:  
`~/Psychophysics/ExtendedFig3/Mooney_invar_Exp1_GSyesonly.jasp`  
`~/Psychophysics/ExtendedFig3/Mooney_invar_Exp2_GSyesonly.jasp`

-----

## **fMRI:**

### *Neural Distance ROI Analysis (Supplementary Fig. 5b):*

ROI analysis: First, randomly bin GLM tstat outputs from 16 runs into 2 bins per “resampled” (or, “pseudo”) run  
`~/fMRI/ExtendedFig5/Analysis_scripts/create_resamplerun.sh`

Subject-level ROI analysis:  
`~/fMRI/ExtendedFig5/Analysis_scripts/ROI_analysis/rsa_plots_crossnobis_resamplerun.py`  
(Alternatively, to run all subjects, use: multi_subject_crossnobis.sh and make sure the correct subjects are listed in `subjects2.txt`)

Group-level ROI analysis (Supplementary Fig. 5b):  
`~/fMRI/ExtendedFig5/Analysis_scripts/ROI_analysis/create_rdm_fromrsatoolboxGENERALv2.m`

Plotting:  
`~/fMRI/ExtendedFig5/plot_FigExtendedFig5b_roianalysis.m`

Source data:  
`~/ExtendedFig5/ROIanalysis_source_data/*`

Includes:  
`meanwithinimg_dist.mat` for averaged within-image, cross-subject distances.  
`SEMwithinimg_dist.mat` for SEM across within-image, cross-subject distances.  
`perm_avgsubj.mat` for subsampled & averaged cross-image, cross-subject distances.

-----

### *Neural Distance Analysis - Searchlight Part 1, invariance trend across HLVC (Supplementary Fig. 5c):*

Subject-level analysis: Main searchlight script-run searchlight analysis in each subject’s T1 space:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/rsa_searchlight_crossnobisV2.py`

Convert all resulting searchlight maps (in T1 space) to standard space:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/t12standard.sh`

Spatially constrain the results using a LO1, LO2, and Fusiform Cortex (“HLVC”) mask:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/constrain_searchlight.sh`  
LO1 and LO2 are taken from the Wang, et. al (2015) probabilistic atlas.  
FC is taken from the Harvard-Oxford (“Occipital temporal fusiform cortex”).

Group-level analysis—to characterize trend of distances across HLVC:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/searchlight_invariance_trend.m`

Source data:  
`~/fMRI/ExtendedFig5/Searchlight_source_data/*`  
`ap_vals.mat` contains averaged within-image distances, for every plane along the posterior-anterior axis in the HLVC.  
`ap_vals_SEM.mat` contains SEM of ap_vals.mat  
`ap_vals_SHUFFLE.mat` contains distances where voxels are spatially shuffled along the posterior-anterior axis.  
`ap_vals_SEM_SHUFFLE.mat` contains SEM of ap_vals_SHUFFLE.mat

-----

### *Neural Distance Analysis - Searchlight Part 2, statistical testing:*

Take resulting searchlight maps from above:  
To generate maps with distributions of cross-image distances per voxel cluster (within each subject):  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/run_wholebrain_rawdist_perms.sh`  
(contains `wholebrain_rawdist_perm.py`)   
This will run 1000 permutations per subject.

To generate permutated values:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/wholebrain_rawdist_EVALUATEperm.py`  
This will give you z-score maps, per subject (both the standard z-map and empirically derived z map. The empirical method first calculates p values as `(sum(perm data < real data) / nReps)`, where nReps = 1000, and converts them to z maps by taking the inverse CDF.

Convert all z maps to standard space:  
`~/fMRI/ExtendedFig5/Analysis_scripts/Searchlight_analysis/t1standard_wholebrain_rawdist.sh`

Convert all z maps to a negative version of it, for a left-tailed one sampled t-test.  
Type in terminal:  
`fslmaths allEMPzmap -mul -1 neg_allEMPzmap`

Group-level analysis: Run the Winkler 2014 method for permutation testing with `fsl randomise`. Type in terminal:  
`randomise -i neg_allEMPzmap -m standardMask -o onesamp_t12mmEMP/less_zero -1 -T -n 5000 -v 12`  
(For spatial smoothing with a 12mm Gaussian kernel (FWHM) with TFCE)

Convert to negative log P and plot in connectome workbench:  
`~/fMRI/ExtendedFig5/General_searchlight_scripts/convert2neglogP.py`

-----

### *Model-based RSA ROI analysis (Fig. 3c):*

Subject level model-based RSA ROI analysis:  
`~/fMRI/Fig3/Analysis_scripts/ROI_analysis/rsa_modelbasedrsa.py`

Group level model-based RSA ROI analysis:  
`~/fMRI/Fig3/Analysis_scripts/ROI_analysis/modelbasedrsa_group.m`

-----

### *Model-based RSA Searchlight analysis for every subject in T1 space (Fig. 3d):*

`~/fMRI/Fig3/Analysis_scripts/Searchlight_analysis/rsa_searchlight_crossnobis_modelbased.py`

Convert to zmap:  
`~/fMRI/Fig3/Analysis_scripts/Searchlight_analysis/wholebrain_modelbased_rtoz.py`

To run a group-level permutation test:  
Convert all searchlight maps to standard & aggregate all maps into a single nifti file:  
`~/fMRI/Fig3/Analysis_scripts/Searchlight_analysis/t12standard_modelbased.sh`

Type this into the terminal (using FSL randomise):   
`randomise -i subjavg_data.nii.gz -m Standard.nii.gz -o onesamp_t/greater_zero -1 -T -n 5000 -v 12`  
(For spatial smoothing with a 12mm Gaussian kernel (FWHM) with TFCE)

Convert to -log(pvalue) and plot in connectome workbench:  
`~/fMRI/Fig3/General_searchlight_scripts/convert2neglogP.py`

-----

### ~/custom_toolbox_python/:

To ensure all dependencies are correctly installed:
1) `git clone git@github.com:BiyuHeLab/Mooney2024.git` to clone this github repository.
2) cd to `custom_toolboxes_python`
3) Run `pipenv install`
4) Run `pipenv shell`
5) Copy over `my_rsatoolbox/rsatoolbox` and overwrite original rsatoolbox files.

-----

### ~/toolboxes_matlab/:

Contains toolboxes used in Matlab analysis:
For plotting: `Violinplot-Matlab-master`, `subaxis`  
For eyetracking analysis: `edf-converter-master`  
For reading nifti files: `vistasoft-master`  
Custom-modified rsatoolbox: `rsatoolbox-main`  
