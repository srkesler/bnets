function permuteMatrixStats_sparselongwithin(groupwithin,FileSaveSuffix_groupname, nperm)
% Shelli Kesler 
% longitudinal permutation testing for edge properties e.g. edge betweenness

%-----Enter info here-------------
load('RegMesAUC.mat'); %nxnxm matrices of edge metrics
if groupwithin==1
    T1 = EdgeAUCG1T1; % matrix for group 1 time 1
    T2 = EdgeAUCG1T2; % matrix for group 1 time 2
else
    T1 = EdgeAUCG2T1; % matrix for group 1 time 1
    T2 = EdgeAUCG2T2; % matrix for group 1 time 2
end
%---------------------------------

trueDiff = mean(T2,3)-mean(T1,3); %diff in slope
trueDiffabs = abs(mean(T2,3) - mean(T1,3)); 

for k = 1:nperm
    rng(1000+k)
    randT1 = datasample(T1,size(T1,3),3,'Replace',true); 
    rng(2000+k)
    randT2 = datasample(T2,size(T2,3),3,'Replace',true); 
    
    randDiff(:,:,k) = abs(mean(randT2,3)-mean(randT1,3)); 
end    


pVal = mean(randDiff > trueDiffabs,3);
pVal=triu(pVal);
pVal0=pVal;
pVal0(trueDiffabs==0) = 999; %Remove instances where edge auc is 0 for all subjects. later we ll remove these instances of 999

%% Extract just the upper triangular matrix of pVal0 as a vector (removing the diagonal as well)
pVal0=pVal0';
m=tril(true(size(pVal0)));
m=logical(m-diag(diag(m)));
true_pVal0  = pVal0(m).';
true_pVal0(true_pVal0==999)=[];

%%
pValnan = pVal;
pValnan(trueDiffabs==0) = NaN;
pValnan=triu(pValnan);

[signEdgesFDR, FDRcritP, pValFDR]=fdr_bh(true_pVal0);
%indx = find(signEdgesFDR'==1);

pValSurv = pValnan;
pValSurv(pValSurv <= FDRcritP & pValSurv>0) = 555;
pValSurv(pValSurv ~= 555) = 0;
pValSurv(pValSurv == 555) = 1;
trueDiffFDR = trueDiff.*pValSurv;
trueDiffFDR =trueDiffFDR'; % to convert final fdr result to lower triangle

indx = find(pValSurv'==1);

for j = 1:length(indx)
    [row(j,1),row(j,2)] = ind2sub(size(pValSurv),indx(j));
end

if ~exist('row') row=[];end

save(['PermStats_longwithin_' FileSaveSuffix_groupname],'trueDiff','trueDiffabs','randDiff', 'pVal', 'pVal0', 'true_pVal0',...
    'pValnan','FDRcritP','pValFDR','pValSurv','indx','trueDiffFDR','row');

save(['TrueDiffFDR_longwithin_' FileSaveSuffix_groupname],'trueDiffFDR')
writetable(array2table(trueDiffFDR),['TrueDiffFDR_longwithin_' FileSaveSuffix_groupname '.xlsx'])