function [grpDens,V1rand] = randDensity(Zrand,NSUB)
%Shelli Kesler
%get densities of random bootstrapped networks
% for randperm = 1:size(randG,2)
%     Zrand = randG{1,randperm}; %nxnxm Z matrix from randperm permutation
    for rsub = 1:NSUB
        Grand = Zrand(:,:,rsub); % one subject from Zrand matrix
        Grand(Grand<0)=0;  % set negative edges to zero
        Vrand = Grand(tril(true(length(Grand)),-1)); % extract the lower triangle
        V1rand{rsub} = sort(abs(Vrand),'descend'); % remove duplicate values and sort
            for j = 1:length(V1rand{1,rsub})
                Grand2=Grand;
                Grand2(logical(eye(size(Grand2)))) = 0; % zero the diagonal
                inx = find(Grand<V1rand{rsub}(j));% find values lower than the jth highest value
                Grand2(inx) = 0; Grand2(Grand2>0)=1; % binarize based on j value
                grpDens{1,rsub}(j,1)=density_und(Grand2); % density at each jth value for all subjects in Zrand
            end
     end
    
end