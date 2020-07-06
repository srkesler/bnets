function maxdisc = Test_Discon_NetMesBin(Data,data_log)
%%
%finds the minimum threshold at which the network becomes disconnected
%from NetMes_Bin

%%
%-- updated on August 21,2011 for functionals
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)

%%  Main Body

e1 = load(Data,'NetMes_Bin');NetMes_Bin = e1.NetMes_Bin;
ff = load(data_log,'mat4GAT');
MinMes = ff.mat4GAT.MinThr;MaxMes = ff.mat4GAT.MaxThr;MesStep = ff.mat4GAT.MesStep;
DensInterval = MinMes:MesStep:MaxMes;

[nSbj, nThr] = size(NetMes_Bin);
Disc = zeros(1,nSbj);

for i = 1:nSbj    
    for j = nThr:-1:1        
        deg=NetMes_Bin{i,j}{1,1};        
        if any(deg<1)            
            Disc(i) = DensInterval(j);            
            break            
        end        
    end    
end

fprintf('%-4s\n',['Fragmentation occurs at network density of ' num2str(max(Disc))]);
maxdisc=max(Disc);

if isequal(single(max(Disc)), single(DensInterval(end)))    
    warning('some of the networks remain fragmented at the highest designated density');
end

% % save DiscResults Disc            
fprintf('%-4s\n',' ..........Done '); 