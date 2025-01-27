#Convert searchlight results in T1 space to standard for averaging:

subjects="001 004 005 006 007 008 009 010 012 013"
mask="LOandFG_2highresProb"

subject_array=($subjects)
data_dir=fMRI_Data/Searchlight_standard/All_searchlight/${mask}
root_dir=Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData/

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
    flirt -in $data_dir/${subj}_searchlight_${mask}.nii.gz -ref $reg_dir/standard.nii.gz -init $reg_dir/highres2standard.mat -out $data_dir/${subj}_${mask}standard.nii.gz -applyxfm
done

cd $data_dir

fslmerge -t ${mask}standard_allsubj.nii.gz 001_${mask}standard.nii.gz 004_${mask}standard.nii.gz 005_${mask}standard.nii.gz 006_${mask}standard.nii.gz 007_${mask}standard.nii.gz 008_${mask}standard.nii.gz 009_${mask}standard.nii.gz 010_${mask}standard.nii.gz 012_${mask}standard.nii.gz 013_${mask}standard.nii.gz

fslmaths ${mask}standard_allsubj.nii.gz -div ${#subject_array[@]} ${mask}standard_allsubjAVG.nii.gz
