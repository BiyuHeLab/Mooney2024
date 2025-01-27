import numpy as np
import nibabel as nib
import os


allsubs = ["001", "004", "005", "006", "007", "008", "009", "010", "012", "013"]

for s in allsubs:
    rmap_file = os.path.join(
        "BigPurpleCluster/Searchlight_Spearman_H2d/subject12mm_" + str(s) + "_output.nii.gz"
    )
    rmap_nifti = nib.load(rmap_file)
    rmap = rmap_nifti.get_fdata()

    # fisher's r-to-z transformation
    zmap = 0.5 * np.log((1 + rmap) / (1 - rmap))
    zmap_nifti = nib.Nifti1Image(zmap, affine=rmap_nifti.affine)
    os.chdir("BigPurpleCluster/Searchlight_Spearman_H2d/")
    nib.save(zmap_nifti, f"{s}_modelcompare_zmap.nii.gz")
