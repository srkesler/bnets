function  NetMes_Bin = NetMeasures_Binary(Graph,Null,Iter,varargin)


%%
%-- Hadi Hosseini, March, 2011
%-- updated on Apr 20,2011 for GUI
%-- updated on Mar 2014 (HH method)

%Ref (Hosseini et al.,Plos One 2012)


NetMes_Bin=cell(22,2);

Graph(Graph>0)=1;


NetMes_Bin{1,1}=degrees_und(Graph);
NetMes_Bin{1,2}='degrees';
NetMes_Bin{2,1}=mean(NetMes_Bin{1,1});
NetMes_Bin{2,2}='mean-degree';
NetMes_Bin{1,3}=NetMes_Bin{1,1}/NetMes_Bin{2,1};



NetMes_Bin{3,1} = sum(Graph);
NetMes_Bin{3,2} = 'strengths';
NetMes_Bin{4,1} = mean(NetMes_Bin{3,1});
NetMes_Bin{4,2} = 'mean-strength';
NetMes_Bin{3,3}=NetMes_Bin{3,1}/NetMes_Bin{4,1};


NetMes_Bin{5,1}=assortativity_bin(Graph,0); %changed from assortativity(Graph,0); as the new BCT toolbox doesnt have an assortativity function 
NetMes_Bin{5,2}='assortativity';

NetMes_Bin{6,1}=density_und(Graph);
NetMes_Bin{6,2}='density';


NetMes_Bin{7,1}=clustering_coef_bu(Graph);
NetMes_Bin{7,2}='clustering-coefs'; 
NetMes_Bin{8,1}=mean(NetMes_Bin{7,1});
NetMes_Bin{8,2}='mean-clustering-coef';
NetMes_Bin{7,3}=NetMes_Bin{7,1}/NetMes_Bin{8,1};


NetMes_Bin{9,1}=transitivity_bu(Graph);
NetMes_Bin{9,2}='transitivity';


NetMes_Bin{10,1}=efficiency_bin(Graph); % Changed from efficiency(Graph);
NetMes_Bin{10,2}='global-efficiency';


NetMes_Bin{11,1}=efficiency_bin(Graph,1); %Changed from efficiency(Graph);
NetMes_Bin{11,2}='local-efficiency';
NetMes_Bin{12,1}=mean(NetMes_Bin{11,1});
NetMes_Bin{12,2}='mean-local-efficiency';
NetMes_Bin{11,3}=NetMes_Bin{11,1}/NetMes_Bin{12,1};


[Mi,NetMes_Bin{13,1}]=modularity_und(Graph);
NetMes_Bin{13,2}='modularity';
[Mi_L,NetMes_Bin{14,1}]=community_louvain(Graph);%Changed from modularity_louvain_und(Graph);
NetMes_Bin{14,2}='modularity-louvain';

NetMes_Bin{25,1}=module_degree_zscore(Graph,Mi,0);
NetMes_Bin{25,2}='module_degree_zscore';
NetMes_Bin{26,1}=mean(NetMes_Bin{25,1});
NetMes_Bin{26,2}='mean-moduledegreezscore';
NetMes_Bin{25,3}=NetMes_Bin{25,1}/NetMes_Bin{26,1};



NetMes_Bin{27,1}=participation_coef(Graph,Mi);
NetMes_Bin{27,2}='participation_coefficient';
NetMes_Bin{28,1}=mean(NetMes_Bin{27,1});
NetMes_Bin{28,2}='mean-participation-coefficient';
NetMes_Bin{27,3}=NetMes_Bin{27,1}/NetMes_Bin{28,1};



NetMes_Bin{15,1}=charpath(distance_bin(Graph),0,0);% Changed from charpath(distance_bin(Graph));
NetMes_Bin{15,2}='characteristic-path-length';


NetMes_Bin{16,1}=betweenness_bin(Graph);
NetMes_Bin{16,2}='node-betweenness';
NetMes_Bin{17,1}=mean(NetMes_Bin{16,1});
NetMes_Bin{17,2}='mean-node-betweenness';
NetMes_Bin{16,3}=NetMes_Bin{16,1}/NetMes_Bin{17,1};


NetMes_Bin{18,1}=edge_betweenness_bin(Graph);
NetMes_Bin{18,2}='edge-betweenness';
NetMes_Bin{19,1}=mean(mean(NetMes_Bin{18,1}));
NetMes_Bin{19,2}='mean-edge-betweenness';
NetMes_Bin{18,3}=NetMes_Bin{18,1}/NetMes_Bin{19,1};





%Small worldness

switch Null
    
   
    case 2
        
        [NetMes_Bin{21,1},NetMes_Bin{20,1},NetMes_Bin{22,1}] = RandomNetGen_SameDeg(Graph,Iter);
        NetMes_Bin{20,2}='lambda_N';
        NetMes_Bin{21,2}='gamma_N';
        NetMes_Bin{22,2}='sigma_N';
    
    case 3
        
        if ~isempty(varargin)
            
            Dat = varargin{1};
            
            cor_null_unthr = RandomNetGen_HQS (cov(Dat),Iter);
            
            TargetDens = nnz(triu(Graph));
            
            clust = []; pathl = []; sw= [];
            
            for ee = 1:size(cor_null_unthr,3)
            
                
                cor_null_thr = match_density(cor_null_unthr(:,:,ee),TargetDens);
                
                cor_null_thr(cor_null_thr>0) = 1;
                
                clust(ee) = mean(clustering_coef_bu(cor_null_thr));
                pathl(ee) = charpath(distance_bin(cor_null_thr));
                sw(ee) = clust(ee)./pathl(ee);
                
            end
            
            NetMes_Bin{20,1} = charpath(distance_bin(Graph))./mean(pathl);
            NetMes_Bin{20,2}='lambda_hqs';
            NetMes_Bin{21,1} = mean(clustering_coef_bu(Graph))./mean(clust);
            NetMes_Bin{21,2}='gamma_hqs';
            NetMes_Bin{22,1} = NetMes_Bin{21,1}./NetMes_Bin{20,1};
            NetMes_Bin{22,2}='sigma_hqs';
            
        
        else
            error('provide covariance matrix')
            
        end

end




function G2thr = match_density(G2,TD)
    
    Step = 0.0001;
    
    
    G2max = max(max(G2));
    ss = -1;
    ii = 0;
    
    while ii <= G2max 
        
        ss = ss + 1;
        ii = Step + ss.*Step; %thr
        
        G2(G2<=ii) = 0;
        
        
        NN = nnz(triu(G2));%number of edges    
        
        
        if ss ~= 0
            
            delta_o = delta;
        
        else
            
            delta_o = NN - TD;
            
        end
        
        
        delta = NN - TD;
        
        
        
        if delta <= 0
            
            if abs(delta) <= delta_o
                
                Thresh = ii;
                
                
                ii = G2max + .1;
                
            elseif abs(delta) > delta_o
                
                Thresh = ii- Step;
                
                
                ii = G2max + .1;
                
            end
            
        end
           
    end
    
    G2thr = G2; G2thr(G2thr <= Thresh) = 0;
    
    clear G2 Thresh_Dens Ind
        
end           
    
    
    
    
end
        







