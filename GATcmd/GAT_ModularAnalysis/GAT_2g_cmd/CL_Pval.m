function p = CL_Pval(Mes_perm,Mes_net,MesName,Tail)

%% 
%  accepts the bootstrap distribution (for difference between random 
%  networks' net measures) and the observed value (difference between
%  actual groups' net measures) and returns the p-value for the observed
%  statistics

% Input:
% Mes_perm:m*n with n permutations for a specified variable at m network densities (e.g. 41*1000)
% Mes_net:m*1 with a specified net measure at m network densities (e.g. 41*1)
% Tail: One- vs. two-tail (1 or 2)

% Output:
% p: p-value of the observed measure (if the difference is a smaller than
% the mean of the permutation, a negative p-value! is generated meaning that the difference is negative)  

%% Reference
%  Bootstrap methods and permutation tests (Ch. 14) 
%  (Tim Hesterberg, David. S. Moore, Shaum Monaghan, 
%   Ashley Clipson, Rachel Epestein) 

%%

%-- Hadi Hosseini: May 3,2011


%%
[m,n] = size(Mes_perm);
p=[];

for k=1:m
    
    Mes_rand=Mes_perm(k,:);
    Mes=Mes_net(k);
    
    if isequal(Tail,1)
        
        as=sort(Mes_rand);
        med=mean(Mes_rand);
        
        if isnan(med)
            
            p=[p;NaN];
            
        end
        
        if Mes>=med
            
            [i,j]=find(as>=Mes);
            L=size(i,2);
            
            if isequal(L,0)
            
                L = 1;%0.000000001*n;
                
            end
            
            p=[p;L./n];
            
        end
        
        if Mes<med
        
            [i,j]=find(as<=Mes);
            L=size(i,2);
            
            if isequal(L,0)
            
                L = 1;%0.000000001*n;
                
            end
            
            p=[p;-L./n];
            
        end
        
    elseif isequal(Tail,2)
        
        Mes_rand(isnan(Mes_rand)) = nanmean(Mes_rand);
        MM = mean(Mes_rand);
        Mes_rand = Mes_rand - MM;
        
        as_abs=sort(abs(Mes_rand));
        [ii,jj]=find(as_abs>=abs(Mes - MM));
        LL=size(ii,2);
        
        if isequal(LL,0)
            
            LL = 1;%=0.000000001*n;
            
        end
        
        p=[p;LL./n];
        
    end

    
end
% % % 
% % % if ~isempty(p)
% % %     save([MesName '_' num2str(Tail) 'Tail_pval'],'p')
% % % end


