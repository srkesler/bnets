function swiData = swiFx(ZT)
% Shelli Kesler
% Measure small-world connectome properties
for j=1:size(ZT,2)
    G = ZT{1,j};
    % Generate random graphs and calculate clustering coefficient and
    % characteristic path length for each random graph
    for i = 1:20 
        rgraph = randmio_und_connected(G,1);
        ccrand(i,1) = mean(clustering_coef_bu(rgraph));
        plrand(i,1) = charpath(distance_bin(rgraph));  
    end
    CC = mean(clustering_coef_bu(G))/mean(ccrand);
    PL = charpath(distance_bin(G))/mean(plrand);
    SWI = CC/PL;
    
    swiData(j,1) = SWI;
    swiData(j,2) = CC;
    swiData(j,3) = PL;
    swiData(j,4) = efficiency_bin(G,0); %global efficiency
    swiData(j,5) = mean(efficiency_bin(G,1)); % local efficiency
   % swiData(j,6) = transitivity_bu(G);
    [~,swiData(j,6)] = modularity_und(G);
    % EDIT: add other measures here as needed
end

end

