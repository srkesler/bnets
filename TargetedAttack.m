% Shelli Kesler 8/5/16
% Conducts targeted network attack by iteratively removing nodes in
% descending order of degree.  Global efficiency is calculated after each
% attack iteration and then area under the curve is computed across all
% iterations

clear all;
clc;

[subjs]=textread('/StudyFolder/paths.txt','%s');
for j = 1:length(subjs)
    tic
    load(fullfile(subjs{j,1},'rotcorr_th.mat'));
    % calculate nodal degrees and sort them
    [degs,DI] = sort((degrees_und(rotcorr_th)),'descend');
    %[degs,DI] = sort(unique((degrees_und(rotcorr_th))),'descend');
    % targeted attack
    G = rotcorr_th; % make a copy of the network
    for i = 1:length(DI)
        G(DI(i),:) = 0; G(:,DI(i)) = 0; % remove the node
        TAgeff(:,i) = efficiency_bin(G,0); % calculate global efficiency after each attack
    end
    TAgeffAUC(j,:) = trapz(1:size(G,2),TAgeff); % compute AUC for each subject
    toc
end
save TAresults.mat TAgeff DI TAgeffAUC;
%plot(1:size(G,2),TAgeff)