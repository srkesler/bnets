function Disc = Test_Discon_NetMesBin(NetMes_Bin,DensInterval)

%%
%finds the minimum threshold at which the network becomes disconnected
%from NetMes_Bin

%%
%-- updated on August 21,2011 for functionals

%Ref (Hosseini et al.,Plos One 2012)


%%
GUIinput = 0;

if nargin < 1
    
    Data = spm_select(1,'NetMesBin_f','Select NetMesBin_f_*.mat file that you want to examine');
    e1 = load(Data,'NetMes_Bin');NetMes_Bin = e1.NetMes_Bin;
    
    %DensInterval = input('type the density range and steps associated with the selected NetMesBin_f file (e.g. .1:.01:.2): ');
    data_log = spm_select(1,'mat4GATf*','Select mat4GATf .mat (output from "Network Measures")');
    ff = load(data_log,'mat4GATf');
    MinMes = ff.mat4GATf.MinThr;MaxMes = ff.mat4GATf.MaxThr;MesStep = ff.mat4GATf.MesStep;
    DensInterval = MinMes:MesStep:MaxMes;
    
    GUIinput = 1;
    
end

if isequal(GUIinput,0)
    
    if nargin > 2 
    
        error('too many inputs')
        
    elseif nargin > 1
        
        error('missing inputs')
    
    end
    
end

[nSbj, nThr] = size(NetMes_Bin);
Disc = zeros(1,nSbj);

for i = 1:nSbj
    
    for j = nThr:-1:1
        
        deg=NetMes_Bin{i,j}{1,1};
        
        if any(deg<1)
            
            %Disc(i) = j;
            Disc(i) = DensInterval(j);
            
            break
            
        end
        
    end
    
end

fprintf('%-4s\n',['Fragmentation occurs at network density of ' num2str(max(Disc))]);

if isequal(single(max(Disc)), single(DensInterval(end)))
    
    warning('some of the networks remain fragmented at the highest designated density');

end


save DiscResults Disc
            
fprintf('%-4s\n',' ..........Done ');
            
       
        