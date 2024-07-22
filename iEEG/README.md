# **iEEG**

`~/Code/genManuscriptFigures.py` is the entry point for the code used for the analysis and figures in
this manuscript.

### *Mean ROI Analysis (Fig. 4c, Extended Data Fig. 8)*

Source data:

`~/Analysis/MeanROI/network/MeanROIAnalysis_network_HGP_1001,2-6,9-21.pkl`
Pickled object of custom Python class `MeanROIAnalysis`, defined at `~/Code/pwmv_ecog/analysis/mean_roi.py`.  
Attribute `results` is the mean, SEM, and electrode survival time courses for each ROI and condition.  
Attribute `clusterStats` are the results of the permutation-based cluster correction testing for each ROI and test in the form of an object of class `ClusterCorrectedMultipleDof`, defined at `~/Code/pwmv_ecog/stats.py`.

Analysis:

`~/Code/pwmv_ecog/analysis/mean_roi.py`  
`~/Code/pwmv_ecog/stats.py`, class `ClusterCorrectedMultipleDof` and function `SEM_within_subject`

Plotting:

`~/Code/pwmv_ecog/manuscript/manuscript_figures.py`  
`~/Code/pwmv_ecog/plotting.py`

-----

### *Image Preference Analysis (Fig. 4c, Extended Data Fig. 8)*

Source data:

`~/Analysis/ImagePref/network/ImagePreferenceAnalysis_HGP_network_None_0.1_1001,2-6,9-21.pkl`
Pickled object of custom Python class `ImagePreferenceAnalysis`, defined at `~/Code/pwmv_ecog/analysis/image_preference.py`.  
Attribute `results` is the mean, SEM, test statistic, significance test outcome, and electrode survival time courses for each ROI and test.  
Attribute `clusterStats` are the results of the permutation-based cluster correction testing for each ROI and test in the form of an object of class `ClusterCorrectedMultipleDof`, defined at `~/Code/pwmv_ecog/stats.py`.

Analysis:

`~/Code/pwmv_ecog/analysis/image_preference.py`  
`~/Code/pwmv_ecog/stats.py`, class `ClusterCorrectedMultipleDof` and function `SEM_within_subject`

Plotting:

 `~/Code/pwmv_ecog/manuscript/manuscript_figures.py`  
`~/Code/pwmv_ecog/plotting.py`

-----

### *Behavioral Analysis 1 -- Recognition Accuracy (Extended Data Fig. 6c)*

Analysis:

 `~/Code/pwmv_ecog/manuscript/manuscript_figures.py`, function `task_performance_plot`

Plotting:

`~/Code/pwmv_ecog/manuscript/manuscript_figures.py`, function `task_performance_plot`  
`~/Code/pwmv_ecog/plotting.py`

Statistics:

`~/Analysis/Figures/Manuscript2024/SuppFig6/PerformanceANOVA.jasp`

-----

### *Behavioral Analysis 2 -- Reaction Time (Extended Data Fig. 7a):*

Analysis:

 `~/Code/pwmv_ecog/manuscript/manuscript_figures.py`, function `supp_fig_7`

Plotting:

 `~/Code/pwmv_ecog/manuscript/manuscript_figures.py`, function `supp_fig_7`

-----

### *Electrode Survival Analysis (Extended Data Fig. 7b-c):*

Source data:

`~/Analysis/MeanROI/network/MeanROIAnalysis_network_HGP_1001,2-6,9-21.pkl` and
`~/Analysis/ImagePref/network/ImagePreferenceAnalysis_HGP_network_None_0.1_1001,2-6,9-21.pkl` as
explained above.

Analysis scripts:

`~/Code/pwmv_ecog/analysis/mean_roi.py`  
`~/Code/pwmv_ecog/analysis/image_preference.py`  
`~/Code/pwmv_ecog/stats.py`, class `ClusterCorrectedMultipleDof`

Plotting:

`~/Code/pwmv_ecog/manuscript/manuscript_figures.py`, function `manuscript_meanROI_survival_plot`  
`~/Code/pwmv_ecog/analysis/image_preference.py`, function `imagePrefs_survival_plot`

-----

### Other Included Files

`~/Code/pwmv_ecog/data.py`: data organization for iEEG data.  
`~/Code/pwmv_ecog/ROIs.py`: managing iEEG electrode memberships in ROIs defined by niftis in `iEEG/Code/Electrodes/ROI_Masks`.
