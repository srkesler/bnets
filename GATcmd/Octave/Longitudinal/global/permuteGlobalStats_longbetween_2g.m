function permuteGlobalStats_longbetween_2g (nperm)
% Shelli Kesler May 2020
% longitudinal permutation testing for global properties
% run GAT 4 groups first with G1T1, G1T2, G2T1, G2T2
% run this script from AUC_Results folder


%-----Enter info here-------------
FileSaveSuffix = 'long_global';
Meas = {'MClust','MDeg','Trans','Assort','GEff','MLocEff','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBetw','MNodeBetw'}; % Measures of interest
%---------------------------------

load('AUC_Results_G1_vs_G2/AUC_NetMesPerDens_AveAcrossDens.mat')
load('AUC_Results_G3_vs_G4/AUC_NetMesPerDens_AveAcrossDens.mat')

for i = 1:length(Meas)
    slopeG1(i,:) = nanmean(eval(['auc_' Meas{i} '_2']) - eval(['auc_' Meas{i} '_1']));
    slopeG2(i,:) = nanmean(eval(['auc_' Meas{i} '_4']) - eval(['auc_' Meas{i} '_3']));
end

trueDiff = slopeG2-slopeG1;
trueDiffabs= abs(trueDiff);

for i = 1:nperm
    for j = 1:length(Meas)
        randSlope1(j,i) = nanmean(eval(['auc_' Meas{j} '_2_rand' '(:,:,i)']) - eval(['auc_' Meas{j} '_1_rand' '(:,:,i)']));
        randSlope2(j,i) = nanmean(eval(['auc_' Meas{j} '_4_rand' '(:,:,i)']) - eval(['auc_' Meas{j} '_3_rand' '(:,:,i)']));
    end
end

randDiff = abs(randSlope2-randSlope1);

for i = 1:length(Meas)
    pVal(i,:) = mean(randDiff(i,:) > trueDiffabs(i,:));
end

[~, FDRcritP, pValFDR]=fdr_bh(pVal);

% calculate the 95% confidence intervals
lolim = .025*nperm;
hilim = nperm-lolim;
for i = 1:length(Meas)
    confint(i,:) = sort(randDiff(i,:));
    CI(i,1) = confint(i,lolim);
    CI(i,2) = confint(i,hilim);
end

rownames = {'MClust','MDegree','Trans','Assort','Eglob','MEloc','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBtwn','MNodeBtwn'}'; %names of measures
colnames = {'VariableNames','meanDiff','LL','UL','p','pFDR'};
temp = cat(2,trueDiff,CI,pVal,pValFDR);

resultsTable(1,:)=colnames;
resultsTable(2:14,1)=rownames;
resultsTable(2:end,2:end)=num2cell(temp);
xlswrite(['resultsTable_' FileSaveSuffix '.xlsx'],resultsTable)

save(['PermStats_' FileSaveSuffix '.mat'], 'trueDiff', 'trueDiffabs','randDiff',...
    'pVal','pValFDR','CI','FDRcritP', 'resultsTable','-v7');