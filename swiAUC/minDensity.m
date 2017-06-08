function [minD,V1,D] = minDensity(Z)
%Shelli Kesler

%find minimum connection density across subjects
for isubj=1:size(Z,3)
    R1 = Z(:,:,isubj);
    R1(R1<0)=0;  % set negative edges to zero
    V = R1(tril(true(length(R1)),-1)); % extract the lower triangle
    V1{isubj} = sort(abs(V),'descend'); % remove duplicate values and sort
        for j = 1:length(V1{1,isubj})
            G1 = R1; % make a copy of matrix
            G1(logical(eye(size(G1)))) = 0; % zero the diagonal
            inx = find(R1<V1{isubj}(j));% find values lower than the jth highest value
            G1(inx) = 0; G1(G1>0)=1; % binarize based on j value
            D{1,isubj}(j,1)=density_und(G1); % record density of each graph
            degs = degrees_und(G1); % calculate nodal degrees
            Dsc=0;
                if any(degs<1) % determine if any nodes are isolated
                    Dsc = 1;
                end
                disconnected(j,:)=Dsc;
        end  
        if all(disconnected==1)
            LD = length(disconnected);
            thresh1(isubj)=V1{1,isubj}(LD); % all disconnected, record lowest threshold
        elseif all(disconnected==0)  %no disconnection, record highest threshold
            thresh1(isubj)=V1{1,isubj}(1);
        else
            thresh1(isubj)=V1{1,isubj}(find(disconnected,1,'last')+1);
        end
              
    % binarize the matrix at minimum density and record the density
    G2=abs(R1); G2(logical(eye(size(G2)))) = 0;
    inx2 = find(G2<thresh1(isubj));
    G2(inx2) = 0; G2(G2>0)=1;
    iMD(isubj)=density_und(G2);
    clear disconnected;
end

[~, ix]=max(iMD); 
minD = iMD(ix); % common minimum density
end

