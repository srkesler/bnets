% Shelli Kesler
% Vikram Rao
% 2/29/16

clear all;

% Threshold connectivity matrices at a common minimum density (data driven)
% Requires BCT version 2016_01_16 (https://sites.google.com/site/bctnet/)

% Input:  unthresholded connectivity matrices
% Output: minD.mat which contains:
%            iMD: minimum density for each subject
%            cMD: common minimum density across subjects
%            sMD: correlation coefficient or other stat/value corresponding to common minimum density
%            GT:  binary matrices thresholded at a common minimum density

tic
% paths to subject data
subjs = {...
'/path1/'
'/path2/'
};

M = 'brain.mat'; % matrix file name

% concatenate matrices
for i = 1:length(subjs)
   G{i} = load(fullfile(subjs{i},M));
end;
% calculate densities across all possible thresholds in each subject's
% matrix and identify minimum connection density for each subject
for isubj=1:length(subjs)
    R1=G{1,isubj}.(char(fields(G{1,1})));
    V = R1(tril(true(length(R1)),-1)); 
    V = round(V*10^2)/10^2;
    V1{isubj} = sort(unique(abs(V)),'descend'); 
        for j = 1:length(V1{1,isubj})
            G1 = R1; 
            G1(logical(eye(size(G1)))) = 0; 
            inx = find(R1<V1{isubj}(j));
            G1(inx) = 0; G1(G1>0)=1; 
            D{1,isubj}(j,1)=density_und(G1); 
            degs = degrees_und(G1); 
            Dsc=0;
                if any(degs<1) % determine if any nodes are isolated
                    Dsc = 1;
                end
                disconnected(j,:)=Dsc;
        end   
    thresh1(isubj)=V1{1,isubj}(find(disconnected,1,'last')+1); 
    % binarize the matrix at minimum density and record the density
    G2=abs(R1); G2(logical(eye(size(G2)))) = 0;
    inx2 = find(G2<thresh1(isubj));
    G2(inx2) = 0; G2(G2>0)=1;
    iMD(isubj)=density_und(G2);
    clear disconnected;
end

% Select thresholds based on highest min density across subjects and save thresholded matrices
[~, ix]=max(iMD); % Find the index with highest min density
cMD = round(iMD(ix)*10^2)/10^2; 
for nsubj=1:size(G,2)
    V2=V1{1,nsubj}; D2=D{1,nsubj}; D2=round(D2*10^2)/10^2;
    [~,id] = min(abs(D2-cMD)); sMD(nsubj)=V2(id); 
    G3=abs(G{1,nsubj}.(char(fields(G{1,1}))));
    G3(logical(eye(size(G3)))) = 0; 
    inx3 = find(G3<sMD(:,nsubj)); G3(inx3) = 0; G3(G3>0)=1; 
    GT{nsubj}=G3;  
end
save('minD.mat','iMD', 'cMD', 'sMD', 'GT', '-v7.3');
toc
