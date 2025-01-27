#Constrain searchlight voxels to LO1/2 + FC only for permutation outputs.
#These outputs are used in searchlight_invariance_trend.m
subjects="001 004 005 006 007 008 009 010 012 013"

for subj in $subjects; do
    for p in {1..1000}; do
    echo $subj
    if [ $subj == "006" ]; then runname='05'
        else runname='Main1'
    fi
    loROI=/isilon/LFMI/VMdrive/Ayaka/Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_winnertakeall_rois
    subj_dir=/isilon/LFMI/VMdrive/Ayaka/Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/$subj
    searchlightOrigin=/isilon/LFMI/VMdrive/Ayaka/Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/Searchlight_standard/All_searchlight/$subj/perm
    outputdir=/isilon/LFMI/VMdrive/Ayaka/Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/Searchlight_standard/All_searchlight/LOandFG/perm
    fslmaths $searchlightOrigin/Subj${subj}_SL_2standard.nii.gz -mas $loROI/LOandFG_3per_2mm.nii.gz $outputdir/${subj}/${p-1}_LOandFG_3per_2mm.nii.gz    
    done
done
