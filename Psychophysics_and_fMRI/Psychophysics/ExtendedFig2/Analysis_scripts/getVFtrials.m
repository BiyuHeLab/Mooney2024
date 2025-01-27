%% This is to get a list of VF trials for eyetrack analysis:

function [VFtrials] = getVFtrials(sortedManipList,sortedImgsList)
%Find which trials had manipulation # 4 (which is VF trials)
VFtrials = zeros(length(sortedManipList),1);
for trial = 1:length(sortedManipList)
    if sortedManipList(trial) == 4
        VFtrials (trial) = sortedImgsList(trial);
    end
end

%get rid of zeros
VFtrials = VFtrials(VFtrials~=0);

end

