function save_info(FolderName,info)
% save data_info without overwriting .mat (to save ROIs etc)

if ~exist(FolderName,'dir')
    mkdir(FolderName)
end

if exist([FolderName,'/exp_info.mat'], 'file')
    save([FolderName,'/exp_info.mat'],'info','-append');
else
    save([FolderName,'/exp_info.mat'],'info');
end
end