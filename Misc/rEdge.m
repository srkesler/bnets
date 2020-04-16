% Dalu Yang
% Shelli Kesler
% 2/2/15
% updated 2/5/16

%% Random edge attack with efficiencies as attack response
% Requires BCT version 2016_01_16 (https://sites.google.com/site/bctnet/)
% Reference: Kesler, SR, et al. 2015. Neurobiol Aging, 36, pp. 2429-42.

load('datafile.mat','datavar')

rng(100); %set seed

for i = 1:size(datavar,3)
    Network = datavar(:,:,i);           
    Netv = squareform(Network);         
    V_edges_all = find(Netv);           
    n_edges_all = length(V_edges_all);  
    p = randperm(n_edges_all);  %randomize the order of indices to attack
    results.attack_edges{i} = p';       
    GE = zeros(n_edges_all,1);          
    LE = zeros(n_edges_all,1);          
    
    for j = 1:n_edges_all   
    %set one edge to 0 (all the edges between these two nodes are disconnected)            			
        Netv(V_edges_all(p(j))) = 0;   				
        Network = squareform(Netv);     			
        GE(j,1) = efficiency_bin(Network);  		%compute global efficiency
        LE(j,1) = mean(efficiency_bin(Network,1));  %compute mean local efficiency
    end
    results.global_efficiency{i} = GE;  
    results.mean_local_efficiency{i} = LE;
    results.AUCge(i) = trapz(GE)./length(GE);   %compute AUCs
    results.AUCle(i) = trapz(LE)./length(LE);
        
end

    
    
