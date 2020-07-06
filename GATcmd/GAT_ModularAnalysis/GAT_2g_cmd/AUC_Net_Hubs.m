function AUC_Net_Hubs(data_log,Data1,Data2,SD)

%%
% This function accept the
% NetMes: the output (NetMes_B1 or NetMes_B2) network measures from NetMEs_vs_Density...
% Dens: The NetMes_B1/B2 index corresponding to minimum density for fully connected graph (output from KAN_FixedDensity)
% ROIs: name of the ROIs associated with each node (output regionName from rex)
% Param: The criterion used for defining hubs. available options are:
%       -'deg' for Degree
%       -'betw' for Node Betweenness
%SD: the threshold for accepting a node as hub. e.g SD=2 means that nodes
%    whose degrees (or betweenness) are at least 2SD greater than mean will be considered as hub

%Output:
%HubNames: Name of the hubs
%HubInd: index of the hubs (node number)

%%
%-- Hadi Hosseini: Created April 12, 2011
%-- Hadi Hosseini: Updated for GUI April 21,2011
%-- Hadi Hosseini: Updated for functionals Sept. 13,2011
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)


%% Main Body
ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;
Group2 = ff.mat4GAT.g2;
load(Data1);load(Data2);
fprintf('%-4s\n',['loading inputs...']);
ROIs = ff.mat4GAT.roi1;

dd=pwd;
mkdir('Hubs_AUC');
cd([dd '/Hubs_AUC']);

Mdeg1norm=mean(AUC_MDeg_1norm);
SDdeg1norm=std(AUC_MDeg_1norm);
iDeg1norm=find(AUC_MDeg_1norm >= (Mdeg1norm+SDdeg1norm*SD));
HubInd1norm=iDeg1norm;

if ~isempty(HubInd1norm)
    count=0;
    for i=HubInd1norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Degnorm_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1norm');
clear HubNames

Mbetw1norm=mean(AUC_MNodeBetw_1norm);
SDbetw1norm=std(AUC_MNodeBetw_1norm);
iBetw1norm=find(AUC_MNodeBetw_1norm >= (Mbetw1norm+SDbetw1norm*SD));
HubInd1norm=iBetw1norm;

if ~isempty(HubInd1norm)
    count=0;
    for i=HubInd1norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Betwnorm_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1norm');

clear HubNames
Mclust1norm=mean(AUC_MClust_1norm);
SDclust1norm=std(AUC_MClust_1norm);
iDeg1norm=find(AUC_MClust_1norm >= (Mclust1norm+SDclust1norm*SD));
HubInd1norm=iDeg1norm;

if ~isempty(HubInd1norm)
    count=0;
    for i=HubInd1norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Clustnorm_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1norm');

clear HubNames

Mdeg1=mean(AUC_MDeg_1);
SDdeg1=std(AUC_MDeg_1);
iDeg1=find(AUC_MDeg_1 >= (Mdeg1+SDdeg1*SD));
HubInd1=iDeg1;

if ~isempty(HubInd1)
    count=0;
    for i=HubInd1
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Deg_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1');

clear HubNames
Mbetw1=mean(AUC_MNodeBetw_1);
SDbetw1=std(AUC_MNodeBetw_1);
iBetw1=find(AUC_MNodeBetw_1 >= (Mbetw1+SDbetw1*SD));
HubInd1=iBetw1;

if ~isempty(HubInd1)
    count=0;
    for i=HubInd1
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Betw_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1');

clear HubNames
Mclust1=mean(AUC_MClust_1);
SDclust1=std(AUC_MClust_1);
iDeg1=find(AUC_MClust_1 >= (Mclust1+SDclust1*SD));
HubInd1=iDeg1;

if ~isempty(HubInd1)
    count=0;
    for i=HubInd1
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Clust_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1');

%%
clear HubNames
Mdeg2norm=mean(AUC_MDeg_2norm);
SDdeg2norm=std(AUC_MDeg_2norm);
iDeg2norm=find(AUC_MDeg_2norm >= (Mdeg2norm+SDdeg2norm*SD));
HubInd2norm=iDeg2norm;

if ~isempty(HubInd2norm)
    count=0;
    for i=HubInd2norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Degnorm_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2norm');
clear HubNames

Mbetw2norm=mean(AUC_MNodeBetw_2norm);
SDbetw2norm=std(AUC_MNodeBetw_2norm);
iBetw2norm=find(AUC_MNodeBetw_2norm >= (Mbetw2norm+SDbetw2norm*SD));
HubInd2norm=iBetw2norm;

if ~isempty(HubInd2norm)
    count=0;
    for i=HubInd2norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Betwnorm_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2norm');
clear HubNames
Mclust2norm=mean(AUC_MClust_2norm);
SDclust2norm=std(AUC_MClust_2norm);
iDeg2norm=find(AUC_MClust_2norm >= (Mclust2norm+SDclust2norm*SD));
HubInd2norm=iDeg2norm;

if ~isempty(HubInd2norm)
    count=0;
    for i=HubInd2norm
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Clustnorm_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2norm');
clear HubNames

Mdeg2=mean(AUC_MDeg_2);
SDdeg2=std(AUC_MDeg_2);
iDeg2=find(AUC_MDeg_2 >= (Mdeg2+SDdeg2*SD));
HubInd2=iDeg2;

if ~isempty(HubInd2)
    count=0;
    for i=HubInd2
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Deg_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2');

clear HubNames

Mbetw2=mean(AUC_MNodeBetw_2);
SDbetw2=std(AUC_MNodeBetw_2);
iBetw2=find(AUC_MNodeBetw_2 >= (Mbetw2+SDbetw2*SD));
HubInd2=iBetw2;

if ~isempty(HubInd2)
    count=0;
    for i=HubInd2
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end

save(['Net_Hubs_Betw_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2');
clear HubNames
Mclust2=mean(AUC_MClust_2);
SDclust2=std(AUC_MClust_2);
iDeg2=find(AUC_MClust_2 >= (Mclust2+SDclust2*SD));
HubInd2=iDeg2;

if ~isempty(HubInd2)
    count=0;
    for i=HubInd2
        count=count+1;
        HubNames{count}=ROIs{i};
    end
else
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Clust_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2');
clear HubNames

fprintf('%-4s\n','.......... done');