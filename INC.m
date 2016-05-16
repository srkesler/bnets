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
ngroups=length(cell2mat(strfind(fieldnames(NCdata),'G')));
    
%% Loops across groups
for igroup=1:ngroups
    disp(sprintf('working on group%d',igroup)) % Displays the group the script is working on
    
    %% Computing Correlation matrix with every sample in the group
    loo.(sprintf('group%d',igroup)).data=NCdata.(sprintf('G%d',igroup));%Storing Original data in the new structure
    loo.(sprintf('group%d',igroup)).SCN=corr(NCdata.(sprintf('G%d',igroup))); % Auto Correlate data to find SCN with all subjects in a group
    
    %% Computing Correlation matrix without one sample at a time in the group
    for i = 1:size(NCdata.(sprintf('G%d',igroup)),1) % loops across subjects for Group1
        DI = NCdata.(sprintf('G%d',igroup)); DI(i,:) = []; %create a data matrix without subject i
        CI = corr(DI); %compute the correlation matrix without subject i
        %Mantel test to get the statistics and pvalue between the correlation matrices with and without subject i
        [loo.(sprintf('group%d',igroup)).mantelscore(i,1),loo.(sprintf('group%d',igroup)).mantelscorepval(i,1)] = bramila_mantel(CI,loo.(sprintf('group%d',igroup)).SCN,5000,'pearson');
    end
    %% Computing the contribution scores defined as 1-(mantel test score) and the average contribution in a group
    loo.(sprintf('group%d',igroup)).contribution = 1-loo.(sprintf('group%d',igroup)).mantelscore(:,1); %contribution score of 1 indicates highest contribution, 0 indicates lowest contribution
    loo.(sprintf('group%d',igroup)).avgcontribution=mean(loo.(sprintf('group%d',igroup)).contribution);
end
save(fullfile(root,'loomantelcontrib.mat'),'loo') % saves all the computed variables in a structure 'loo'


