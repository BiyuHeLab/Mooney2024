import nibabel as nib
import numpy as np

# Load the target image and the mask
target_img = nib.load(
    "Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/Searchlight_standard/All_searchlight/VTC_noEVC/VTC_noEVCstandard_allsubjAVG.nii.gz"
)
mask_img = nib.load(
    "Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_brain_rois_combined/VTC_noEVC.nii.gz"
)

# Get the data from both images
target_data = target_img.get_fdata()
mask_data = mask_img.get_fdata()

# Set voxels in the target image to NaN where the mask is 0
target_data[mask_data == 0] = np.nan

# Create a new NIfTI image for the modified data
new_img = nib.Nifti1Image(target_data, affine=target_img.affine)

# Save the new image
nib.save(
    new_img,
    "Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/Searchlight_standard/All_searchlight/VTC_noEVC/VTC_SL_final.nii.gz",
)
