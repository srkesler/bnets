% Shelli Kesler
clear all;
% Calculate small-worldness properties for individual level connectivity matrices
% Reference: Bassett, DS, et al. 2006. The neuroscientist, 12, pp. 512-23.

% Input:  binary or weighted, undirected and thresholded matrices
% Output: swiData.mat which contains:
%           1: smallworldness index
%           2: normalized clustering coefficient
%           3: normalized characteristic path length
%           4: global efficiency
%           5: mean local efficiency
%           6: network size
          
% Requires BCT version 2016_01_16 (https://sites.google.com/site/bctnet/)

% start parallel pool if not already running
isempty(gcp);

[dirs] = textread('/Users/skesler/Desktop/tracks/Paths.txt' , '%s');
% Txt file contains full paths to each subj e.g. /Volumes/STUDY/SUB001

for j=1:length(dirs)
    cd(dirs{j});
    load('G.mat') % matrix file name

    % Generate random graphs and calculate clustering coefficient and
    % characteristic path length for each random graph
    tic
    parfor i = 1:20 
        rgraph = randmio_und_connected(G,1);
        
        % for binary graphs
        ccrand(i,1) = mean(clustering_coef_bu(rgraph));
        plrand(i,1) = charpath(distance_bin(rgraph));  
        
        % for weighted graphs
        % ccrand(i,1) = mean(clustering_coef_wu(rgraph));
        % plrand(i,1) = charpath(distance_wei(rgraph));
        
    end

    % normalize clustering coefficient and characteristic path length by
    % these same values obtained from random graphs to compute
    % small-worldness index
    
    % for binary graphs
    CC = mean(clustering_coef_bu(G))/mean(ccrand);
    PL = charpath(distance_bin(G))/mean(plrand);
    
    %for weighted graphs
    %CC = mean(clustering_coef_wu(G))/mean(ccrand);
    %PL = charpath(distance_wei(G))/mean(plrand);
    
    SWI = CC/PL;
    
    swiData(j,1) = SWI;
    swiData(j,2) = CC;
    swiData(j,3) = PL;
    
    % for binary graphs
    swiData(j,4) = efficiency_bin(G,0); %global efficiency
    swiData(j,5) = mean(efficiency_bin(G,1)); % local efficiency
    
    % for weighted graphs
    %swiData(j,4) = efficiency_wei(G,0); %global efficiency
    %swiData(j,5) = mean(efficiency_wei(G,1)); % local efficiency
    
    swiData(j,6) = size(G,2); % size

    toc
end
save('swiData.mat', 'swiData');
