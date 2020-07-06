root='/Users/Desktop/Study/cogsub/Modules/';
copyroot='/Users/Desktop/Study/cogsub/';
mkdir(root);
time=3; %number of timepoints
groups={'exp1'; 'exp2'; 'exp3'};

for i=1: time
    for j=1: length(groups)
    
        mkdir(fullfile(root,sprintf('T%d',i),sprintf('%s',groups{j})))
        cd(fullfile(copyroot,sprintf('tp%d',i),'analysis'))
        tmp=dir('conn*');
        cd(tmp(1).name);
        cd(fullfile('modularanalysis',sprintf('%s',groups{j})))
        tmpmod=dir('module*');
        for k=1:length(tmpmod)
            if isdir(fullfile(sprintf('module%d',k),'Regional','Hubs_AUC'))
                copyfile(fullfile(sprintf('module%d',k),'Regional','Hubs_AUC','HubSummary.csv'),fullfile(root,sprintf('T%d',i),sprintf('%s',groups{j}),sprintf('Mod%d_HubSummary.csv',k)))
                copyfile(fullfile(sprintf('module%d',k),'AUC_Results','AUC_Mean.mat'),fullfile(root,sprintf('T%d',i),sprintf('%s',groups{j}),sprintf('Mod%d_AUC_Mean.mat',k)))
            end
        end
    
    end
end

