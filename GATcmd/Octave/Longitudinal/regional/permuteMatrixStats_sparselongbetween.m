function permuteMatrixStats_sparselongbetween (nperm)
% Shelli Kesler 
% longitudinal permutation testing for edge properties e.g. edge betweenness

%-----Enter info here-------------
load('RegMesAUC.mat'); %nxnxm matrices of edge metrics
G1T1 = EdgeAUCG1T1; % matrix for group 1 time 1
G1T2 = EdgeAUCG1T2; % matrix for group 1 time 2
G2T1 = EdgeAUCG2T1; % matrix for group 2 time 1
G2T2 = EdgeAUCG2T2; % matrix for group 2 time 2
FileSaveSuffix = 'long_between';
%---------------------------------

trueDiff = (mean(G2T2,3)-mean(G2T1,3)) - (mean(G1T2,3)-mean(G1T1,3)); %diff in slope
trueDiffabs = abs((mean(G2T2,3)-mean(G2T1,3)) - (mean(G1T2,3)-mean(G1T1,3))); %diff in slope

for k = 1:nperm
  sprintf('%d',k)
    rand("seed",1000+k);
    rndIdxG1T1 = randsample(size(G1T1,3),size(G1T1,3),true);
    randG1T1=G1T1(:,:,rndIdxG1T1);
    
    rand("seed",2000+k);
    rndIdxG1T2 = randsample(size(G1T2,3),size(G1T2,3),true);
    randG1T2=G1T2(:,:,rndIdxG1T2);
    
    rand("seed",3000+k);
    rndIdxG2T1 = randsample(size(G2T1,3),size(G2T1,3),true);
    randG2T1=G2T1(:,:,rndIdxG2T1);
    
    rand("seed",4000+k);
    rndIdxG2T2 = randsample(size(G2T2,3),size(G2T2,3),true);
    randG2T2=G2T2(:,:,rndIdxG2T2);
    
    randDiff(:,:,k) = abs((mean(randG2T2,3)-mean(randG2T1,3)) - (mean(randG1T2,3)-mean(randG1T1,3)));   
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

save(['PermStats_' FileSaveSuffix '.mat'], 'trueDiff', 'trueDiffabs', 'randDiff', 'pVal', ...
    'pVal0', 'true_pVal0', 'pValnan', 'FDRcritP', 'pValFDR', 'pValSurv', 'indx', ...
    'trueDiffFDR', 'row','-v7');
save(['TrueDiffFDR_' FileSaveSuffix '.mat'],'trueDiffFDR','-v7');
xlswrite(['TrueDiffFDR_' FileSaveSuffix '.xlsx'],trueDiffFDR)