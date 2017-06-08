function ZT = binGraph(Z,V1,D,S)
% binarize graphs at specific densities
    for i = 1:size(Z,3)
    V2=V1{1,i}; D2=D{1,i}; 
        [~,id] = min(abs(D2-S)); sMD(i)=V2(id); % find best threshold for the subject at S density
        Z3=Z(:,:,i);
        Z3(logical(eye(size(Z3)))) = 0; 
        inx3 = find(Z3<sMD(:,i)); Z3(inx3) = 0; Z3(Z3>0)=1; % binarize
        ZT{i}=Z3; % save the thresholded matrices
    end
end
