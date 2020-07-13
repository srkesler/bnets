function Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)

% data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
% MinMax: mninmum, maximum thresholds and the thresholded step
% nperm: number of permutation (sampling) (default 100)

%-- Hadi Hosseini (Created Apr 12,2011)
%-- Hadi Hosseini (Updated Apr 21,2011 for GUI)
%-- Hadi Hosseini (Updated Sept 14,2011 for functionals)
%-- updated on Mar 2014 (parallel processing)
%-- Shelli Kesler 10/30/15 corrected parallel processing
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao
%Ref (Hosseini et al.,Plos One 2012)

%% Main body common to all group combinations

if isempty(nperm)
    nperm=100;
end
ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;
ROI = ff.mat4GAT.roi1;

if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha')
    Alpha = ff.mat4GAT.alpha;
    Tail = ff.mat4GAT.tail;
else
    Alpha = .05;
    Tail = 2;
end

MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;

fprintf('%-4s\n',['loading inputs...']);

%% Group1 Vs Group2

f1=load(Data1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
f2=load(Data2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;

Sz1=size(NetMes1,1);Sz2=size(NetMes2,1);
data=NetMes1;data(Sz1+1:Sz1+Sz2,:)=NetMes2;
rand("seed",1001); % Added by Vikram Rao on 03/22/2020 in order to fix the seed and have consistent results
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

NetMes1_rand=cell(1,nperm);NetMes2_rand=cell(1,nperm);

for i=1:nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rand("seed",1000+i); % Added by Vikram Rao on 03/22/2020 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    rand("seed",2000+i); % Added by Vikram Rao on 03/22/2020 in order to fix the seed and have consistent results
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    NetMes1_rand{i}=Randata(Samp1,:);
    NetMes2_rand{i}=Randata(Samp2,:);
end

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];

MinThr=ff.mat4GAT.MinThr;
MinIdx=find(single(xxx)==single(MinThr));
MaxThr=ff.mat4GAT.MaxThr;
MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(NetMes1,2)
    MaxIdx = size(NetMes1,2);
end

if isempty(MaxIdx) || isempty(MinIdx)
    errordlg('the selected density range should correspond to the density range specified in the original analysis!',...
        'Error', 'modal');
    return
end

Xax=MinMesPlot:MesStepPlot:MaxMesPlot;
Xax=Xax(MinIdx:MaxIdx);

dd=pwd;
mkdir('Regional');
cd([dd '/Regional']);

fprintf('%-4s\n',' calculating regional network measures....');

NetMes1=NetMes1(:,MinIdx:MaxIdx);
NetMes2=NetMes2(:,MinIdx:MaxIdx);

MClust_1norm=[];MDeg_1norm=[];MNodeBetw_1norm=[];MLocEff_1norm=[];
MClust_2norm=[];MDeg_2norm=[];MNodeBetw_2norm=[];MLocEff_2norm=[];
AUC_MClust_1norm=[];AUC_MDeg_1norm=[];AUC_MNodeBetw_1norm=[];AUC_MLocEff_1norm=[];
AUC_MClust_2norm=[];AUC_MDeg_2norm=[];AUC_MNodeBetw_2norm=[];AUC_MLocEff_2norm=[];

fda_MClust_1norm=[];fda_MDeg_1norm=[];fda_MNodeBetw_1norm=[];fda_MLocEff_1norm=[];
fda_MClust_2norm=[];fda_MDeg_2norm=[];fda_MNodeBetw_2norm=[];fda_MLocEff_2norm=[];

for i=1:size(NetMes1,1)
    fprintf('%-4s\n',['calculating group1 subject ' num2str(i) ' regional network measures....']);
    temp_clust1=[];temp_deg1=[];temp_nodeb1=[];temp_leff1=[];
    
    for j=1:size(NetMes1,2)
        temp_clust1=[temp_clust1;NetMes1{i,j}{7,3}'];
        temp_deg1=[temp_deg1;NetMes1{i,j}{1,3}];
        temp_nodeb1=[temp_nodeb1;NetMes1{i,j}{16,3}];
        temp_leff1=[temp_leff1;NetMes1{i,j}{11,3}'];
    end
    
    MClust_1norm=[MClust_1norm;mean(temp_clust1)];
    MDeg_1norm=[MDeg_1norm;mean(temp_deg1)];
    MNodeBetw_1norm=[MNodeBetw_1norm;mean(temp_nodeb1)];
    MLocEff_1norm=[MLocEff_1norm;mean(temp_leff1)];
    
    AUC_MClust_1norm=[AUC_MClust_1norm;trapz(Xax,temp_clust1)];
    AUC_MDeg_1norm=[AUC_MDeg_1norm;trapz(Xax,temp_deg1)];
    AUC_MNodeBetw_1norm=[AUC_MNodeBetw_1norm;trapz(Xax,temp_nodeb1)];
    AUC_MLocEff_1norm=[AUC_MLocEff_1norm;trapz(Xax,temp_leff1)];
    
    fda_MClust_1norm=[fda_MClust_1norm;sum(temp_clust1)];
    fda_MDeg_1norm=[fda_MDeg_1norm;sum(temp_deg1)];
    fda_MNodeBetw_1norm=[fda_MNodeBetw_1norm;sum(temp_nodeb1)];
    fda_MLocEff_1norm=[fda_MLocEff_1norm;sum(temp_leff1)];
end

save(['Indiv_NetMesReg_' Group1 '.mat'],'MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm','-v7');
save(['Indiv_AUC_NetMesReg_' Group1 '.mat'],'AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm','-v7');
save(['Indiv_fda_NetMesReg_' Group1 '.mat'],'fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm','-v7');

MClust_1norm=mean(MClust_1norm);
MDeg_1norm=mean(MDeg_1norm);
MNodeBetw_1norm=mean(MNodeBetw_1norm);
MLocEff_1norm=mean(MLocEff_1norm);

AUC_MClust_1norm=mean(AUC_MClust_1norm);
AUC_MDeg_1norm=mean(AUC_MDeg_1norm);
AUC_MNodeBetw_1norm=mean(AUC_MNodeBetw_1norm);
AUC_MLocEff_1norm=mean(AUC_MLocEff_1norm);

fda_MClust_1norm=mean(fda_MClust_1norm);
fda_MDeg_1norm=mean(fda_MDeg_1norm);
fda_MNodeBetw_1norm=mean(fda_MNodeBetw_1norm);
fda_MLocEff_1norm=mean(fda_MLocEff_1norm);

for i=1:size(NetMes2,1)
    fprintf('%-4s\n',['calculating group2 subject ' num2str(i) ' regional network measures....']);
    temp_clust2=[];temp_deg2=[];temp_nodeb2=[];temp_leff2=[];
    
    for j=1:size(NetMes2,2)
        temp_clust2=[temp_clust2;NetMes2{i,j}{7,3}'];
        temp_deg2=[temp_deg2;NetMes2{i,j}{1,3}];
        temp_nodeb2=[temp_nodeb2;NetMes2{i,j}{16,3}];
        temp_leff2=[temp_leff2;NetMes2{i,j}{11,3}'];
    end
    
    MClust_2norm=[MClust_2norm;mean(temp_clust2)];
    MDeg_2norm=[MDeg_2norm;mean(temp_deg2)];
    MNodeBetw_2norm=[MNodeBetw_2norm;mean(temp_nodeb2)];
    MLocEff_2norm=[MLocEff_2norm;mean(temp_leff2)];
    
    AUC_MClust_2norm=[AUC_MClust_2norm;trapz(Xax,temp_clust2)];
    AUC_MDeg_2norm=[AUC_MDeg_2norm;trapz(Xax,temp_deg2)];
    AUC_MNodeBetw_2norm=[AUC_MNodeBetw_2norm;trapz(Xax,temp_nodeb2)];
    AUC_MLocEff_2norm=[AUC_MLocEff_2norm;trapz(Xax,temp_leff2)];
    
    fda_MClust_2norm=[fda_MClust_2norm;sum(temp_clust2)];
    fda_MDeg_2norm=[fda_MDeg_2norm;sum(temp_deg2)];
    fda_MNodeBetw_2norm=[fda_MNodeBetw_2norm;sum(temp_nodeb2)];
    fda_MLocEff_2norm=[fda_MLocEff_2norm;sum(temp_leff2)];
end

save(['Indiv_NetMesReg_' Group2 '.mat'],'MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm','-v7');
save(['Indiv_AUC_NetMesReg_' Group2 '.mat'],'AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm', '-v7');
save(['Indiv_fda_NetMesReg_' Group2 '.mat'],'fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm', '-v7');

MClust_2norm=mean(MClust_2norm);
MDeg_2norm=mean(MDeg_2norm);
MNodeBetw_2norm=mean(MNodeBetw_2norm);
MLocEff_2norm=mean(MLocEff_2norm);

AUC_MClust_2norm=mean(AUC_MClust_2norm);
AUC_MDeg_2norm=mean(AUC_MDeg_2norm);
AUC_MNodeBetw_2norm=mean(AUC_MNodeBetw_2norm);
AUC_MLocEff_2norm=mean(AUC_MLocEff_2norm);

fda_MClust_2norm=mean(fda_MClust_2norm);
fda_MDeg_2norm=mean(fda_MDeg_2norm);
fda_MNodeBetw_2norm=mean(fda_MNodeBetw_2norm);
fda_MLocEff_2norm=mean(fda_MLocEff_2norm);

MClust_1_randnorm=[];MDeg_1_randnorm=[];MNodeBetw_1_randnorm=[];MLocEff_1_randnorm=[];
MClust_2_randnorm=[];MDeg_2_randnorm=[];MNodeBetw_2_randnorm=[];MLocEff_2_randnorm=[];

AUC_MClust_1_randnorm=[];AUC_MDeg_1_randnorm=[];AUC_MNodeBetw_1_randnorm=[];AUC_MLocEff_1_randnorm=[];
AUC_MClust_2_randnorm=[];AUC_MDeg_2_randnorm=[];AUC_MNodeBetw_2_randnorm=[];AUC_MLocEff_2_randnorm=[];

fda_MClust_1_randnorm=[];fda_MDeg_1_randnorm=[];fda_MNodeBetw_1_randnorm=[];fda_MLocEff_1_randnorm=[];
fda_MClust_2_randnorm=[];fda_MDeg_2_randnorm=[];fda_MNodeBetw_2_randnorm=[];fda_MLocEff_2_randnorm=[];

for k=1:size(NetMes1_rand,2)
    fprintf('%-4s\n',['calculating random net ' num2str(k) ' regional network measures....']);
    temp_rand1=NetMes1_rand{1,k};
    temp_rand2=NetMes2_rand{1,k};
    temp_rand1=temp_rand1(:,MinIdx:MaxIdx);
    temp_rand2=temp_rand2(:,MinIdx:MaxIdx);
    
    MClust_1norm_rand=[];MDeg_1norm_rand=[];MNodeBetw_1norm_rand=[];MLocEff_1norm_rand=[];
    MClust_2norm_rand=[];MDeg_2norm_rand=[];MNodeBetw_2norm_rand=[];MLocEff_2norm_rand=[];
    AUC_MClust_1norm_rand=[];AUC_MDeg_1norm_rand=[];AUC_MNodeBetw_1norm_rand=[];AUC_MLocEff_1norm_rand=[];
    AUC_MClust_2norm_rand=[];AUC_MDeg_2norm_rand=[];AUC_MNodeBetw_2norm_rand=[];AUC_MLocEff_2norm_rand=[];
    fda_MClust_1norm_rand=[];fda_MDeg_1norm_rand=[];fda_MNodeBetw_1norm_rand=[];fda_MLocEff_1norm_rand=[];
    fda_MClust_2norm_rand=[];fda_MDeg_2norm_rand=[];fda_MNodeBetw_2norm_rand=[];fda_MLocEff_2norm_rand=[];
    
    for i=1:size(temp_rand1,1)
        temp_rand_clust1=[];temp_rand_deg1=[];temp_rand_nodeb1=[];temp_rand_leff1=[];
        
        for j=1:size(temp_rand1,2)
            temp_rand_clust1=[temp_rand_clust1;temp_rand1{i,j}{7,3}'];
            temp_rand_deg1=[temp_rand_deg1;temp_rand1{i,j}{1,3}];
            temp_rand_nodeb1=[temp_rand_nodeb1;temp_rand1{i,j}{16,3}];
            temp_rand_leff1=[temp_rand_leff1;temp_rand1{i,j}{11,3}'];
        end
        
        MClust_1norm_rand=[MClust_1norm_rand;mean(temp_rand_clust1)];
        MDeg_1norm_rand=[MDeg_1norm_rand;mean(temp_rand_deg1)];
        MNodeBetw_1norm_rand=[MNodeBetw_1norm_rand;mean(temp_rand_nodeb1)];
        MLocEff_1norm_rand=[MLocEff_1norm_rand;mean(temp_rand_leff1)];
        
        AUC_MClust_1norm_rand=[AUC_MClust_1norm_rand;trapz(Xax,temp_rand_clust1)];
        AUC_MDeg_1norm_rand=[AUC_MDeg_1norm_rand;trapz(Xax,temp_rand_deg1)];
        AUC_MNodeBetw_1norm_rand=[AUC_MNodeBetw_1norm_rand;trapz(Xax,temp_rand_nodeb1)];
        AUC_MLocEff_1norm_rand=[AUC_MLocEff_1norm_rand;trapz(Xax,temp_rand_leff1)];
        
        fda_MClust_1norm_rand=[fda_MClust_1norm_rand;sum(temp_rand_clust1)];
        fda_MDeg_1norm_rand=[fda_MDeg_1norm_rand;sum(temp_rand_deg1)];
        fda_MNodeBetw_1norm_rand=[fda_MNodeBetw_1norm_rand;sum(temp_rand_nodeb1)];
        fda_MLocEff_1norm_rand=[fda_MLocEff_1norm_rand;sum(temp_rand_leff1)];
    end
    
    MClust_1_randnorm=[MClust_1_randnorm;mean(MClust_1norm_rand)];
    MDeg_1_randnorm=[MDeg_1_randnorm;mean(MDeg_1norm_rand)];
    MNodeBetw_1_randnorm=[MNodeBetw_1_randnorm;mean(MNodeBetw_1norm_rand)];
    MLocEff_1_randnorm=[MLocEff_1_randnorm;mean(MLocEff_1norm_rand)];
    
    AUC_MClust_1_randnorm=[AUC_MClust_1_randnorm;mean(AUC_MClust_1norm_rand)];
    AUC_MDeg_1_randnorm=[AUC_MDeg_1_randnorm;mean(AUC_MDeg_1norm_rand)];
    AUC_MNodeBetw_1_randnorm=[AUC_MNodeBetw_1_randnorm;mean(AUC_MNodeBetw_1norm_rand)];
    AUC_MLocEff_1_randnorm=[AUC_MLocEff_1_randnorm;mean(AUC_MLocEff_1norm_rand)];
    
    fda_MClust_1_randnorm=[fda_MClust_1_randnorm;mean(fda_MClust_1norm_rand)];
    fda_MDeg_1_randnorm=[fda_MDeg_1_randnorm;mean(fda_MDeg_1norm_rand)];
    fda_MNodeBetw_1_randnorm=[fda_MNodeBetw_1_randnorm;mean(fda_MNodeBetw_1norm_rand)];
    fda_MLocEff_1_randnorm=[fda_MLocEff_1_randnorm;mean(fda_MLocEff_1norm_rand)];
    
    for i=1:size(temp_rand2,1)
        temp_rand_clust2=[];temp_rand_deg2=[];temp_rand_nodeb2=[];temp_rand_leff2=[];
        for j=1:size(temp_rand2,2)
            temp_rand_clust2=[temp_rand_clust2;temp_rand2{i,j}{7,3}'];
            temp_rand_deg2=[temp_rand_deg2;temp_rand2{i,j}{1,3}];
            temp_rand_nodeb2=[temp_rand_nodeb2;temp_rand2{i,j}{16,3}];
            temp_rand_leff2=[temp_rand_leff2;temp_rand2{i,j}{11,3}'];
        end
        
        MClust_2norm_rand=[MClust_2norm_rand;mean(temp_rand_clust2)];
        MDeg_2norm_rand=[MDeg_2norm_rand;mean(temp_rand_deg2)];
        MNodeBetw_2norm_rand=[MNodeBetw_2norm_rand;mean(temp_rand_nodeb2)];
        MLocEff_2norm_rand=[MLocEff_2norm_rand;mean(temp_rand_leff2)];
        
        AUC_MClust_2norm_rand=[AUC_MClust_2norm_rand;trapz(Xax,temp_rand_clust2)];
        AUC_MDeg_2norm_rand=[AUC_MDeg_2norm_rand;trapz(Xax,temp_rand_deg2)];
        AUC_MNodeBetw_2norm_rand=[AUC_MNodeBetw_2norm_rand;trapz(Xax,temp_rand_nodeb2)];
        AUC_MLocEff_2norm_rand=[AUC_MLocEff_2norm_rand;trapz(Xax,temp_rand_leff2)];
        
        fda_MClust_2norm_rand=[fda_MClust_2norm_rand;sum(temp_rand_clust2)];
        fda_MDeg_2norm_rand=[fda_MDeg_2norm_rand;sum(temp_rand_deg2)];
        fda_MNodeBetw_2norm_rand=[fda_MNodeBetw_2norm_rand;sum(temp_rand_nodeb2)];
        fda_MLocEff_2norm_rand=[fda_MLocEff_2norm_rand;sum(temp_rand_leff2)];
    end
    
    MClust_2_randnorm=[MClust_2_randnorm;mean(MClust_2norm_rand)];
    MDeg_2_randnorm=[MDeg_2_randnorm;mean(MDeg_2norm_rand)];
    MNodeBetw_2_randnorm=[MNodeBetw_2_randnorm;mean(MNodeBetw_2norm_rand)];
    MLocEff_2_randnorm=[MLocEff_2_randnorm;mean(MLocEff_2norm_rand)];
    
    AUC_MClust_2_randnorm=[AUC_MClust_2_randnorm;mean(AUC_MClust_2norm_rand)];
    AUC_MDeg_2_randnorm=[AUC_MDeg_2_randnorm;mean(AUC_MDeg_2norm_rand)];
    AUC_MNodeBetw_2_randnorm=[AUC_MNodeBetw_2_randnorm;mean(AUC_MNodeBetw_2norm_rand)];
    AUC_MLocEff_2_randnorm=[AUC_MLocEff_2_randnorm;mean(AUC_MLocEff_2norm_rand)];
    
    fda_MClust_2_randnorm=[fda_MClust_2_randnorm;mean(fda_MClust_2norm_rand)];
    fda_MDeg_2_randnorm=[fda_MDeg_2_randnorm;mean(fda_MDeg_2norm_rand)];
    fda_MNodeBetw_2_randnorm=[fda_MNodeBetw_2_randnorm;mean(fda_MNodeBetw_2norm_rand)];
    fda_MLocEff_2_randnorm=[fda_MLocEff_2_randnorm;mean(fda_MLocEff_2norm_rand)];
end

nROI=size(MClust_1norm,2);

save(['NetMesReg_rand_' Group1 '.mat'],'MClust_1_randnorm','MDeg_1_randnorm','MNodeBetw_1_randnorm','MLocEff_1_randnorm','-v7');
save(['NetMesReg_rand_' Group2 '.mat'],'MClust_2_randnorm','MDeg_2_randnorm','MNodeBetw_2_randnorm','MLocEff_2_randnorm','-v7');
save(['AUC_NetMesReg_rand_' Group1 '.mat'],'AUC_MClust_1_randnorm','AUC_MDeg_1_randnorm','AUC_MNodeBetw_1_randnorm','AUC_MLocEff_1_randnorm' ,'-v7');
save(['AUC_NetMesReg_rand_' Group2 '.mat'],'AUC_MClust_2_randnorm','AUC_MDeg_2_randnorm','AUC_MNodeBetw_2_randnorm','AUC_MLocEff_2_randnorm','-v7');
save(['fda_NetMesReg_rand_' Group1 '.mat'],'fda_MClust_1_randnorm','fda_MDeg_1_randnorm','fda_MNodeBetw_1_randnorm','fda_MLocEff_1_randnorm','-v7');
save(['fda_NetMesReg_rand_' Group2 '.mat'],'fda_MClust_2_randnorm','fda_MDeg_2_randnorm','fda_MNodeBetw_2_randnorm','fda_MLocEff_2_randnorm','-v7');

save(['NetMesReg_' Group1 '.mat'],'MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm','-v7');
save(['NetMesReg_' Group2 '.mat'],'MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm','-v7');
save(['AUC_NetMesReg_' Group1 '.mat'],'AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm','-v7');
save(['AUC_NetMesReg_' Group2 '.mat'],'AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm','-v7');
save(['fda_NetMesReg_' Group1 '.mat'],'fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm','-v7');
save(['fda_NetMesReg_' Group2 '.mat'],'fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm','-v7');

[mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(MClust_1_randnorm,0.05);
[mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(MDeg_1_randnorm,0.05);
[mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(MNodeBetw_1_randnorm,0.05);
[mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(MLocEff_1_randnorm,0.05);

[mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(MClust_2_randnorm,0.05);
[mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(MDeg_2_randnorm,0.05);
[mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(MNodeBetw_2_randnorm,0.05);
[mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(MLocEff_2_randnorm,0.05);

[AUC_mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(AUC_MClust_1_randnorm,0.05);
[AUC_mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(AUC_MDeg_1_randnorm,0.05);
[AUC_mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(AUC_MNodeBetw_1_randnorm,0.05);
[AUC_mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(AUC_MLocEff_1_randnorm,0.05);

[AUC_mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(AUC_MClust_2_randnorm,0.05);
[AUC_mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(AUC_MDeg_2_randnorm,0.05);
[AUC_mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(AUC_MNodeBetw_2_randnorm,0.05);
[AUC_mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(AUC_MLocEff_2_randnorm,0.05);

[fda_mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(fda_MClust_1_randnorm,0.05);
[fda_mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(fda_MDeg_1_randnorm,0.05);
[fda_mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(fda_MNodeBetw_1_randnorm,0.05);
[fda_mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(fda_MLocEff_1_randnorm,0.05);

[fda_mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(fda_MClust_2_randnorm,0.05);
[fda_mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(fda_MDeg_2_randnorm,0.05);
[fda_mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(fda_MNodeBetw_2_randnorm,0.05);
[fda_mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(fda_MLocEff_2_randnorm,0.05);

trueDiff_MClust = MClust_2norm-MClust_1norm;
trueDiff_MDeg = MDeg_2norm-MDeg_1norm;
trueDiff_MNodeBetw = MNodeBetw_2norm-MNodeBetw_1norm;
trueDiff_MLocEff = MLocEff_2norm-MLocEff_1norm;

trueDiff_AUC_MClust = AUC_MClust_2norm-AUC_MClust_1norm;
trueDiff_AUC_MDeg = AUC_MDeg_2norm-AUC_MDeg_1norm;
trueDiff_AUC_MNodeBetw = AUC_MNodeBetw_2norm-AUC_MNodeBetw_1norm;
trueDiff_AUC_MLocEff = AUC_MLocEff_2norm-AUC_MLocEff_1norm;

trueDiff_fda_MClust = fda_MClust_2norm-fda_MClust_1norm;
trueDiff_fda_MDeg = fda_MDeg_2norm-fda_MDeg_1norm;
trueDiff_fda_MNodeBetw = fda_MNodeBetw_2norm-fda_MNodeBetw_1norm;
trueDiff_fda_MLocEff = fda_MLocEff_2norm-fda_MLocEff_1norm;


save truediff.mat -v7 trueDiff_MClust trueDiff_MDeg trueDiff_MNodeBetw trueDiff_MLocEff ...
    trueDiff_AUC_MClust trueDiff_AUC_MDeg trueDiff_AUC_MNodeBetw trueDiff_AUC_MLocEff ....
    trueDiff_fda_MClust trueDiff_fda_MDeg trueDiff_fda_MNodeBetw trueDiff_fda_MLocEff
    
N_rand=size(NetMes1_rand,2);

Pvalue = Alpha;

Ci_MClustnorm=CL_per(MClust_2_randnorm-MClust_1_randnorm,Pvalue);
Ci_MDegnorm=CL_per(MDeg_2_randnorm-MDeg_1_randnorm,Pvalue);
Ci_MNodeBetwnorm=CL_per(MNodeBetw_2_randnorm-MNodeBetw_1_randnorm,Pvalue);
Ci_MLocEffnorm=CL_per(MLocEff_2_randnorm-MLocEff_1_randnorm,Pvalue);

AUC_Ci_MClustnorm=CL_per(AUC_MClust_2_randnorm-AUC_MClust_1_randnorm,Pvalue);
AUC_Ci_MDegnorm=CL_per(AUC_MDeg_2_randnorm-AUC_MDeg_1_randnorm,Pvalue);
AUC_Ci_MNodeBetwnorm=CL_per(AUC_MNodeBetw_2_randnorm-AUC_MNodeBetw_1_randnorm,Pvalue);
AUC_Ci_MLocEffnorm=CL_per(AUC_MLocEff_2_randnorm-AUC_MLocEff_1_randnorm,Pvalue);

fda_Ci_MClustnorm=CL_per(fda_MClust_2_randnorm-fda_MClust_1_randnorm,Pvalue);
fda_Ci_MDegnorm=CL_per(fda_MDeg_2_randnorm-fda_MDeg_1_randnorm,Pvalue);
fda_Ci_MNodeBetwnorm=CL_per(fda_MNodeBetw_2_randnorm-fda_MNodeBetw_1_randnorm,Pvalue);
fda_Ci_MLocEffnorm=CL_per(fda_MLocEff_2_randnorm-fda_MLocEff_1_randnorm,Pvalue);

p_RegClust_norm = CL_Pval((MClust_2_randnorm-MClust_1_randnorm)',(MClust_2norm'-MClust_1norm'),'RegClustNorm',Tail);
p_RegDeg_norm = CL_Pval((MDeg_2_randnorm-MDeg_1_randnorm)',(MDeg_2norm'-MDeg_1norm'),'RegDegNorm',Tail);
p_RegNodeBetw_norm = CL_Pval((MNodeBetw_2_randnorm-MNodeBetw_1_randnorm)',(MNodeBetw_2norm'-MNodeBetw_1norm'),'RegNodeBetwNorm',Tail);
p_RegLocEff_norm = CL_Pval((MLocEff_2_randnorm-MLocEff_1_randnorm)',(MLocEff_2norm'-MLocEff_1norm'),'RegLocEffNorm',Tail);

p_AUC_RegClust_norm = CL_Pval((AUC_MClust_2_randnorm-AUC_MClust_1_randnorm)',(AUC_MClust_2norm'-AUC_MClust_1norm'),'AUC_RegClustNorm',Tail);
p_AUC_RegDeg_norm = CL_Pval((AUC_MDeg_2_randnorm-AUC_MDeg_1_randnorm)',(AUC_MDeg_2norm'-AUC_MDeg_1norm'),'AUC_RegDegNorm',Tail);
p_AUC_RegNodeBetw_norm = CL_Pval((AUC_MNodeBetw_2_randnorm-AUC_MNodeBetw_1_randnorm)',(AUC_MNodeBetw_2norm'-AUC_MNodeBetw_1norm'),'AUC_RegNodeBetwNorm',Tail);
p_AUC_RegLocEff_norm = CL_Pval((AUC_MLocEff_2_randnorm-AUC_MLocEff_1_randnorm)',(AUC_MLocEff_2norm'-AUC_MLocEff_1norm'),'AUC_RegLocEffNorm',Tail);

p_fda_RegClust_norm = CL_Pval((fda_MClust_2_randnorm-fda_MClust_1_randnorm)',(fda_MClust_2norm'-fda_MClust_1norm'),'fda_RegClustNorm',Tail);
p_fda_RegDeg_norm = CL_Pval((fda_MDeg_2_randnorm-fda_MDeg_1_randnorm)',(fda_MDeg_2norm'-fda_MDeg_1norm'),'fda_RegDegNorm',Tail);
p_fda_RegNodeBetw_norm = CL_Pval((fda_MNodeBetw_2_randnorm-fda_MNodeBetw_1_randnorm)',(fda_MNodeBetw_2norm'-fda_MNodeBetw_1norm'),'fda_RegNodeBetwNorm',Tail);
p_fda_RegLocEff_norm = CL_Pval((fda_MLocEff_2_randnorm-fda_MLocEff_1_randnorm)',(fda_MLocEff_2norm'-fda_MLocEff_1norm'),'fda_RegLocEffNorm',Tail);

save pvals.mat -v7 p_RegClust_norm p_RegDeg_norm p_RegNodeBetw_norm p_RegLocEff_norm p_AUC_RegClust_norm p_AUC_RegDeg_norm p_AUC_RegNodeBetw_norm ...
    p_AUC_RegLocEff_norm p_fda_RegClust_norm p_fda_RegDeg_norm p_fda_RegNodeBetw_norm p_fda_RegLocEff_norm

[~,~,trueDiff_MClust_FDR_pval]=fdr_bh(abs(p_RegClust_norm),Alpha);
[~,~,trueDiff_MDeg_FDR_pval]=fdr_bh(abs(p_RegDeg_norm),Alpha);
[~,~,trueDiff_MNodeBetw_FDR_pval]=fdr_bh(abs(p_RegNodeBetw_norm),Alpha);
[~,~,trueDiff_MLocEff_FDR_pval]=fdr_bh(abs(p_RegLocEff_norm),Alpha);

[~,~,trueDiff_AUC_MClust_FDR_pval]=fdr_bh(abs(p_AUC_RegClust_norm),Alpha);
[~,~,trueDiff_AUC_MDeg_FDR_pval]=fdr_bh(abs(p_AUC_RegDeg_norm),Alpha);
[~,~,trueDiff_AUC_MNodeBetw_FDR_pval]=fdr_bh(abs(p_AUC_RegNodeBetw_norm),Alpha);
[~,~,trueDiff_AUC_MLocEff_FDR_pval]=fdr_bh(abs(p_AUC_RegLocEff_norm),Alpha);


[~,~,trueDiff_fda_MClust_FDR_pval]=fdr_bh(abs(p_fda_RegClust_norm),Alpha);
[~,~,trueDiff_fda_MDeg_FDR_pval]=fdr_bh(abs(p_fda_RegDeg_norm),Alpha);
[~,~,trueDiff_fda_MNodeBetw_FDR_pval]=fdr_bh(abs(p_fda_RegNodeBetw_norm),Alpha);
[~,~,trueDiff_fda_MLocEff_FDR_pval]=fdr_bh(abs(p_fda_RegLocEff_norm),Alpha);

save pvals_FDR.mat -v7 trueDiff_MClust_FDR_pval trueDiff_MDeg_FDR_pval trueDiff_MNodeBetw_FDR_pval trueDiff_MLocEff_FDR_pval...
    trueDiff_AUC_MClust_FDR_pval trueDiff_AUC_MDeg_FDR_pval trueDiff_AUC_MNodeBetw_FDR_pval trueDiff_AUC_MLocEff_FDR_pval ...
    trueDiff_fda_MClust_FDR_pval trueDiff_fda_MDeg_FDR_pval trueDiff_fda_MNodeBetw_FDR_pval trueDiff_fda_MLocEff_FDR_pval

regmes_auc =num2cell([trueDiff_AUC_MClust' trueDiff_AUC_MClust_FDR_pval trueDiff_AUC_MDeg' trueDiff_AUC_MDeg_FDR_pval trueDiff_AUC_MNodeBetw' trueDiff_AUC_MNodeBetw_FDR_pval trueDiff_AUC_MLocEff' trueDiff_AUC_MLocEff_FDR_pval]);
regmes_auc=horzcat(ff.mat4GAT.roi1',regmes_auc);
regmes_auc(2:length(ROI)+1,1:9)=regmes_auc(1:length(ROI),1:9);
regmes_auc(1,1:9) = {'ROI_name','Diff_AUC_MClust', 'FDR_pval_AUC_MClust', 'Diff_AUC_MDegree', 'FDR_pval_AUC_MDegree', 'Diff_AUC_MNodeBtwn', 'FDR_pval_AUC_MNodeBtwn', 'Diff_AUC_MEloc', 'FDR_pval_AUC_MEloc'};
save RegMes.mat -v7 regmes_auc
xlswrite('RegMes.xlsx',regmes_auc);

cd ..

fprintf('%-4s\n','.... done ....');