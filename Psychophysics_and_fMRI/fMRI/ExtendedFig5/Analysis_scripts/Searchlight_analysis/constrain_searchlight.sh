#Constrain searchlight voxels to LO1/2 + FC only.
subjects="001 004 005 006 007 008 009 010 012 013"

for subj in $subjects; do
    echo $subj
    if [ $subj == "006" ]; then runname='05'
        else runname='Main1'
    fi

    loROI=Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_winnertakeall_rois
    subj_dir=fMRI_Data/$subj
    searchlightOrigin=All_searchlight
    outputdir=/All_searchlight/LOandFG
    fslmaths $searchlightOrigin/Subj${subj}_SL_2standard.nii.gz -mas $loROI/LOandFG_3per_2mm.nii.gz $outputdir/${subj}_LOandFG_3per_2mm.nii.gz
done

flirt -in Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_prob_rois/LOandFG.nii.gz -ref Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/005/Processed/Session1/Main1/1stLevel_GLM_Main1.feat/reg/standard.nii.gz -out Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_prob_rois/LOandFG_2mm.nii.gz -applyisoxfm 2

cd Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/Searchlight_standard/All_searchlight/LOandFG/

fslmerge -t LOandFG_3per2mm_searchlight_allsubs 001_LOandFG_3per_2mm 004_LOandFG_3per_2mm 005_LOandFG_3per_2mm 006_LOandFG_3per_2mm 007_LOandFG_3per_2mm 008_LOandFG_3per_2mm 009_LOandFG_3per_2mm 010_LOandFG_3per_2mm 012_LOandFG_3per_2mm 013_LOandFG_3per_2mm
