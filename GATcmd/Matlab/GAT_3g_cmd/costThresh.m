function ThreshCorr = costThresh(Corr,Cost)

%Input:
%   Corr: a 1xN cell array containing N correlation matrices
%   Corr{i}: R(nroi,nroi) 

%Output:
%   ThreshCorr: a cell array containing N binary correaltion matrices
%   thrresholded to have a specific density (=Cost)

%% -- Hadi Hosseini, March 2014


Corr = cellfun(@zero_diag,Corr,'UniformOutput', false);
Corr = cellfun(@neg_zero,Corr,'UniformOutput', false);

nnet = size(Corr,2);
nroi = size(Corr{1},1);
ind_triu = find(triu(ones(nroi,nroi),1));

parfor i = 1:nnet
    
    R = Corr{i};
    MaxD = density_und(R);
    
    R_ = zeros(nroi,nroi);
    [srt,ind_srt] = sort(R(ind_triu),'descend');
    
    if Cost <= MaxD 
        
        R_(ind_triu(ind_srt(1:round(Cost*(nroi*(nroi-1)/2))))) = 1;
        R_ = R_ + R_';

        ThreshCorr{i} = R_;
        
    else
        
        R_(ind_triu(ind_srt(1:floor(MaxD*(nroi*(nroi-1)/2))))) = 1;
        R_ = R_ + R_';

        ThreshCorr{i} = R_;
    
    end
end

