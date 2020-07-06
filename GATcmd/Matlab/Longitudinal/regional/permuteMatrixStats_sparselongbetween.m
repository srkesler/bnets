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
    rng(1000+k)
    randG1T1 = datasample(G1T1,size(G1T1,3),3,'Replace',true);
    rng(2000+k)
    randG1T2 = datasample(G1T2,size(G1T2,3),3,'Replace',true); 
    rng(3000+k)
    randG2T1 = datasample(G2T1,size(G2T1,3),3,'Replace',true);
    rng(4000+k)
    randG2T2 = datasample(G2T2,size(G2T2,3),3,'Replace',true);
    
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

save(['PermStats_' FileSaveSuffix],'trueDiff', 'trueDiffabs', 'randDiff', 'pVal', 'pVal0', 'true_pVal0',...
    'pValnan','FDRcritP','pValFDR','pValSurv','indx','trueDiffFDR','row');


save(['TrueDiffFDR_' FileSaveSuffix],'trueDiffFDR')
writetable(array2table(trueDiffFDR),['TrueDiffFDR_' FileSaveSuffix '.xlsx'])