% Vikram Rao
% Shelli Kesler
% Dalu Yang

% Computes individual network contribution from group level structural correlation networks
% Requires bramilla_mantel.m (https://users.aalto.fi/~eglerean/permutations.html) 
% Reference:  Saggar, M, et al. 2015. NeuroImage, 120, pp. 274-284.

clc
clear all
close all

%% Load the Structural correlation network
[FileName,root,FilterIndex] = uigetfile('*.mat','Select the file that contains the Network Connectivity Data');
if FilterIndex==0 error('File not selected'); end % Errors if the file is not selected
NCdata=load(fullfile(root,FileName)); % Loads the data
if length(fieldnames(NCdata))==1 NCdata=NCdata.(char(fieldnames(NCdata))); end
ngroups=length(cell2mat(strfind(fieldnames(NCdata),'data')));
ngroups=input('Enter the number of groups: ');
for igroup=1:ngroups
    
    ntp.(sprintf('group%d',igroup))=input(sprintf('Enter the number of time points for group%d: ',igroup));
end
%% Routine if there is only 1 timepoint for groups
if ntp.group1==1
    %% Loops across groups
    for igroup=1:ngroups
        disp(sprintf('working on group%d',igroup))
        
        %% Computing Correlation matrix with every sample in the group
        loo.(sprintf('group%d',igroup)).data=NCdata.(sprintf('data%d',igroup));
        loo.(sprintf('group%d',igroup)).groupname=NCdata.(sprintf('g%d',igroup));
        loo.(sprintf('group%d',igroup)).SCN=corr(NCdata.(sprintf('data%d',igroup)));
        
        %% Computing Correlation matrix without one sample at a time in the group
        for i = 1:size(NCdata.(sprintf('data%d',igroup)),1)
            DI = NCdata.(sprintf('data%d',igroup)); DI(i,:) = [];
            CI = corr(DI);
            [loo.(sprintf('group%d',igroup)).mantelscore(i,1),loo.(sprintf('group%d',igroup)).mantelscorepval(i,1)] = bramila_mantel(CI,loo.(sprintf('group%d',igroup)).SCN,5000,'pearson');
        end
        %% Computing the contribution scores defined as 1-(mantel test score) and the average contribution in a group
        loo.(sprintf('group%d',igroup)).contribution = 1-loo.(sprintf('group%d',igroup)).mantelscore(:,1);
        loo.(sprintf('group%d',igroup)).avgcontribution=mean(loo.(sprintf('group%d',igroup)).contribution);
    end
    save(fullfile(root,'loomantelcontrib.mat'),'loo') % saves all the computed variables in a structure 'loo'
    
    %% Routine if there are more than 1 timepoint for groups
else
    %% Loops across groups
    for igroup=1:ngroups
        for itp=1:ntp.(sprintf('group%d',igroup))
            disp(sprintf('working on group%d timepoint%d',igroup,itp))
            
            %% Computing Correlation matrix with every sample in the group
            loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).data=NCdata.(sprintf('data%dt%d',igroup,itp));
            loo.(sprintf('group%d',igroup)).groupname=NCdata.(sprintf('g%d',igroup));%Storing group names in the new structure
            loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).SCN=corr(NCdata.(sprintf('data%dt%d',igroup,itp))); 
            
            %% Computing Correlation matrix without one sample at a time in the group
            for i = 1:size(NCdata.(sprintf('data%dt%d',igroup,itp)),1) 
                DI = NCdata.(sprintf('data%dt%d',igroup,itp)); DI(i,:) = []; 
                CI = corr(DI); %compute the correlation matrix without subject i
                [loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).mantelscore(i,1),loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).mantelscorepval(i,1)] = bramila_mantel(CI,loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).SCN,5000,'pearson');
            end
            %% Computing the contribution scores defined as 1-(mantel test score) and the average contribution in a group
            loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).contribution = 1-loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).mantelscore(:,1); %contribution score of 1 indicates highest contribution, 0 indicates lowest contribution
            loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).avgcontribution=mean(loo.(sprintf('group%d',igroup)).(sprintf('tp%d',itp)).contribution);
        end
        
    end
    save(fullfile(root,'loomantelcontrib.mat'),'loo') % saves all the computed variables in a structure 'loo'
end




