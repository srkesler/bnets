function permuteMatrixStats_sparse(inputname,FileSaveSuffix_groupname, nperm)
% Shelli Kesler 
% permutation testing for edge properties e.g. edge betweenness

%-----Enter info here-------------
load(inputname); %nxnxm matrices of edge metrics
G1 = EdgeAUCG1; % matrix for group 1
G2 = EdgeAUCG2; % matrix for group 2
%---------------------------------

trueDiff = mean(G2,3) - mean(G1,3); % true mean difference
trueDiffabs = abs(mean(G2,3) - mean(G1,3)); % true mean difference

for k = 1:nperm
    rand("seed",1000+k);
    rndIdxG1 = randsample(size(G1,3),size(G1,3),true);
    randG1=G1(:,:,rndIdxG1);
    
    
    rand("seed",2000+k);
    rndIdxG2 = randsample(size(G2,3),size(G2,3),true);
    randG2=G2(:,:,rndIdxG2);
    
    randDiff(:,:,k) = abs(mean(randG2,3) - mean(randG1,3));
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

save(['PermStats_' FileSaveSuffix_groupname '.mat'],'trueDiff', 'trueDiffabs', 'randDiff', 'pVal', 'pVal0', 'true_pVal0',...
    'pValnan','FDRcritP','pValFDR','pValSurv','indx','trueDiffFDR','row','-v7');

save('TrueDiffFDR.mat','trueDiffFDR','-v7');
xlswrite('TrueDiffFDR.xlsx',trueDiffFDR)