#Create "pseudoruns" by averaging 2 runs per pseudorun, randomly chosen for every t value.

subjects="001"
data_dir=Mooney_PriorInvariance/InvarianceMap_fmri/Data/fMRISubjectData

for subj in $subjects; do
    echo $subj
    #Assign run names separately because of mishap in my coding history:
    if [ $subj == "001" ]; then runname="04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19"
        else runname="Main1 Main2 Main3 Main4 Main5 Main6 Main7 Main8 Main9 Main10 Main11 Main12 Main13 Main14 Main15 Main16"
    fi

    #Define directories
    subj_dir=$data_dir/$subj

    #Randomize runs
    string_runname=($runname) #convert array to string
    shuffled_runname=($(shuf -e "${string_runname[@]}")) #shuffle runname array
    echo "${shuffled_runname[@]}" #display for sanity check

    counter=0
    for pseudo in {1..8}; do
        echo $counter
        echo ${shuffled_runname[$counter]}
        echo ${shuffled_runname[$counter + 1]}
        outputdir=$subj_dir/Processed/Session1/PseudoMain${pseudo}
        if [ ! -d "$outputdir" ]; then
            mkdir "$outputdir"
        fi

        runname1=${shuffled_runname[$counter]}
        runname2=${shuffled_runname[$counter + 1]}
        echo $runname1
        echo $runname2
        reg_dir1=$subj_dir/Processed/Session1/${runname1}/1stLevel_GLM_${runname1}.feat/stats
        reg_dir2=$subj_dir/Processed/Session1/${runname2}/1stLevel_GLM_${runname2}.feat/stats
        for t in {1..80}; do
	    fslmaths $reg_dir1/tstat${t}_2highres.nii.gz -add $reg_dir2/tstat${t}_2highres.nii.gz -div 2 $outputdir/tstat${t}_2highres.nii.gz
        done
        counter=$((counter + 2))
    done
done
