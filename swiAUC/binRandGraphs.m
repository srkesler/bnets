function ZTrand = binRandGraphs(Z,V1,D,S,Perm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

  for subj = 1:size(Z,3)
    V2=V1{1,Perm}{1,subj}; 
    D2=D{1,Perm}{1,subj}; 
        [~,id] = min(abs(D2-S)); sMD(subj)=V2(id); % find best threshold for the subject at S density
        Z3=Z(:,:,subj);
        Z3(logical(eye(size(Z3)))) = 0; 
        inx3 = find(Z3<sMD(:,subj)); Z3(inx3) = 0; Z3(Z3>0)=1; % binarize
        ZTrand{subj}=Z3; % save the thresholded matrices
   end
end

