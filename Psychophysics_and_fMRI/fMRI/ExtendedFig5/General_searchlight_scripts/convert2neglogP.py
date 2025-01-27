import numpy as np
import nibabel as nib

map_affine = nib.load("onesamp_t12mm/greater_zero_tfce_corrp_tstat1.nii.gz")
stat_map = nib.load(f"onesamp_t12mm/greater_zero_tfce_corrp_tstat1.nii.gz").get_fdata()
stat_map[stat_map < 0.95] = np.nan  # p<0.05 only

neglogP = -1 * np.log10(1 - stat_map)  # convert to negative log P
output = nib.Nifti1Image(neglogP, map_affine.affine)
nib.save(output, "neg_logP9mm.nii")

# toolboxes/workbench/bin_rh_linux64/wb_command -volume-to-surface-mapping ~/Desktop/Ayaka/BigPurpleCluster/transfer_files0306202024/Searchlight_VFis0point5/standard_outputs/neg_logP9mm.nii toolboxes/workbench/HCP_S1200_GroupAvg_v1/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii toolboxes/workbench/SL_neglogP9mm_left.shape.gii -trilinear

# toolboxes/workbench/bin_rh_linux64/wb_command -volume-to-surface-mapping ~/Desktop/Ayaka/BigPurpleCluster/transfer_files0306202024/Searchlight_VFis0point5/standard_outputs/neg_logP9mm.nii toolboxes/workbench/HCP_S1200_GroupAvg_v1/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii toolboxes/workbench/SL_neglogP9mm_right.shape.gii -trilinear
