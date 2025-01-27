#Convert searchlight results in T1 space to standard for averaging:

subjects="001 004 005 006 007 008 009 010 012 013"

subject_array=($subjects)
data_dir=BigPurpleCluster/Searchlight_Spearman_H2d/
root_dir=Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/
std=Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_brain_rois_combined
prob=Mooney_PriorInvariance/InvarianceMap_fmri/Data/AnalysisPipeline/standard_prob_rois/
#Convert VTC_noEVC to each subject's T1:
for subj in $subjects; do
    echo $subj
    if [[ $subj == "001" || $subj == "006" ]]; then
        runname='05'
    else
        runname='Main1'
    fi
    subj_dir=$root_dir/$subj
    reg_dir=$subj_dir/Processed/Session1/${runname}/1stLevel_GLM_${runname}.feat/reg
    standardbrain=$reg_dir/standard.nii.gz
    #flirt -in $data_dir/${subj}/subject_${subj}_output.nii.gz -ref $reg_dir/standard.nii.gz -init $reg_dir/highres2standard.mat -out $data_dir/${subj}_standard.nii.gz -applyxfm
    flirt -in $data_dir/${subj}_modelcompare_zmap.nii.gz -ref $reg_dir/standard.nii.gz -init $reg_dir/highres2standard.mat -out $data_dir/${subj}_standard.nii.gz -applyxfm
done

cd $data_dir

fslmerge -t standard_allsubj.nii.gz 001_standard.nii.gz 004_standard.nii.gz 006_standard.nii.gz 007_standard.nii.gz  009_standard.nii.gz 010_standard.nii.gz 012_standard.nii.gz 013_standard.nii.gz 008_standard.nii.gz 005_standard.nii.gz

#fslmaths standard_allsubj.nii.gz -div ${#subject_array[@]} modelcompare_allsubjAVG.nii.gz

