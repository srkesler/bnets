function Net_RegDiff_withBootstrap(data_log,data1,data2,Nperm)

% data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
% MinMax: mninmum, maximum thresholds and the thresholded step
% Nperm: number of permutation (sampling) (default 100)

%%

%-- Hadi Hosseini (Created Apr 12,2011)
%-- Hadi Hosseini (Updated Apr 21,2011 for GUI)
%-- Hadi Hosseini (Updated Sept 14,2011 for functionals)
%-- updated on Mar 2014 (parallel processing)
%-- Shelli Kesler 10/30/15 corrected parallel processing
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao
%Ref (Hosseini et al.,Plos One 2012)

%% Main Body

ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;
ROI = ff.mat4GAT.roi1;
Alpha = ff.mat4GAT.alpha;
Tail = ff.mat4GAT.tail;

MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;
fprintf('%-4s\n',['loading inputs...']);
f1=load(data1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
f2=load(data2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;

Sz1=size(NetMes1,1);Sz2=size(NetMes2,1);
data=NetMes1;data(Sz1+1:Sz1+Sz2,:)=NetMes2;
rand("seed",1001); % Added by Vikram Rao on 02/05/2019 in order to fix the seed and have consistent results
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);
NetMes1_rand=cell(1,Nperm);NetMes2_rand=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rand("seed",1000+i); % Added by Vikram Rao on 02/05/2019 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    rand("seed",2000+i); % Added by Vikram Rao on 02/05/2019 in order to fix the seed and have consistent results
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    NetMes1_rand{i}=Randata(Samp1,:);
    NetMes2_rand{i}=Randata(Samp2,:);
end

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=MinMesPlot;MinIdx=find(single(xxx)==single(MinThr));
MaxThr=MaxMesPlot;MaxIdx=find(single(xxx)==single(MaxThr));
Xax = [MinThr:MesStepPlot:MaxThr];

dd=pwd;
mkdir('Regional');
cd([dd '/Regional']);

fprintf('%-4s\n',' calculating regional network measures....');

NetMes1=NetMes1(:,MinIdx:MaxIdx);
NetMes2=NetMes2(:,MinIdx:MaxIdx);

MClust_1norm=[];MDeg_1norm=[];MNodeBetw_1norm=[];MLocEff_1norm=[]; MPartCoeff_1norm=[];
MClust_2norm=[];MDeg_2norm=[];MNodeBetw_2norm=[];MLocEff_2norm=[]; MPartCoeff_2norm=[];
MClust_1=[];MDeg_1=[];MNodeBetw_1=[];MLocEff_1=[]; MPartCoeff_1=[];
MClust_2=[];MDeg_2=[];MNodeBetw_2=[];MLocEff_2=[]; MPartCoeff_2=[];

AUC_MClust_1norm=[];AUC_MDeg_1norm=[];AUC_MNodeBetw_1norm=[];AUC_MLocEff_1norm=[];
AUC_MClust_2norm=[];AUC_MDeg_2norm=[];AUC_MNodeBetw_2norm=[];AUC_MLocEff_2norm=[];

AUC_MClust_1=[];AUC_MDeg_1=[];AUC_MNodeBetw_1=[];AUC_MLocEff_1=[];
AUC_MClust_2=[];AUC_MDeg_2=[];AUC_MNodeBetw_2=[];AUC_MLocEff_2=[];

fda_MClust_1norm=[];fda_MDeg_1norm=[];fda_MNodeBetw_1norm=[];fda_MLocEff_1norm=[];
fda_MClust_2norm=[];fda_MDeg_2norm=[];fda_MNodeBetw_2norm=[];fda_MLocEff_2norm=[];

fda_MClust_1=[];fda_MDeg_1=[];fda_MNodeBetw_1=[];fda_MLocEff_1=[];
fda_MClust_2=[];fda_MDeg_2=[];fda_MNodeBetw_2=[];fda_MLocEff_2=[];

for i=1:size(NetMes1,1)
    fprintf('%-4s\n',['calculating group1 subject ' num2str(i) ' regional network measures....']);
    temp_clust1norm=[];temp_deg1norm=[];temp_nodeb1norm=[];temp_leff1norm=[]; temp_partcoeff1norm=[];
    temp_clust1=[];temp_deg1=[];temp_nodeb1=[];temp_leff1=[]; temp_partcoeff1=[];
    
    for j=1:size(NetMes1,2)
        temp_clust1norm=[temp_clust1norm;NetMes1{i,j}{7,3}'];
        temp_deg1norm=[temp_deg1norm;NetMes1{i,j}{1,3}];
        temp_nodeb1norm=[temp_nodeb1norm;NetMes1{i,j}{16,3}];
        temp_leff1norm=[temp_leff1norm;NetMes1{i,j}{11,3}'];
        temp_partcoeff1norm=[temp_partcoeff1norm;NetMes1{i,j}{27,3}'];
        
        temp_clust1=[temp_clust1;NetMes1{i,j}{7,1}'];
        temp_deg1=[temp_deg1;NetMes1{i,j}{1,1}];
        temp_nodeb1=[temp_nodeb1;NetMes1{i,j}{16,1}];
        temp_leff1=[temp_leff1;NetMes1{i,j}{11,1}'];
        temp_partcoeff1=[temp_partcoeff1;NetMes1{i,j}{27,1}'];
    end
    
    MClust_1norm=[MClust_1norm;mean(temp_clust1norm)];
    MDeg_1norm=[MDeg_1norm;mean(temp_deg1norm)];
    MNodeBetw_1norm=[MNodeBetw_1norm;mean(temp_nodeb1norm)];
    MLocEff_1norm=[MLocEff_1norm;mean(temp_leff1norm)];
    MPartCoeff_1norm=[MPartCoeff_1norm;mean(temp_partcoeff1norm)];
    
    MClust_1=[MClust_1;mean(temp_clust1)];
    MDeg_1=[MDeg_1;mean(temp_deg1)];
    MNodeBetw_1=[MNodeBetw_1;mean(temp_nodeb1)];
    MLocEff_1=[MLocEff_1;mean(temp_leff1)];
    MPartCoeff_1=[MPartCoeff_1;mean(temp_partcoeff1)];
    
    AUC_MClust_1norm=[AUC_MClust_1norm;trapz(Xax,temp_clust1norm)];
    AUC_MDeg_1norm=[AUC_MDeg_1norm;trapz(Xax,temp_deg1norm)];
    AUC_MNodeBetw_1norm=[AUC_MNodeBetw_1norm;trapz(Xax,temp_nodeb1norm)];
    AUC_MLocEff_1norm=[AUC_MLocEff_1norm;trapz(Xax,temp_leff1norm)];
    
    AUC_MClust_1=[AUC_MClust_1;trapz(Xax,temp_clust1)];
    AUC_MDeg_1=[AUC_MDeg_1;trapz(Xax,temp_deg1)];
    AUC_MNodeBetw_1=[AUC_MNodeBetw_1;trapz(Xax,temp_nodeb1)];
    AUC_MLocEff_1=[AUC_MLocEff_1;trapz(Xax,temp_leff1)];
    
    fda_MClust_1norm=[fda_MClust_1norm;sum(temp_clust1norm)];
    fda_MDeg_1norm=[fda_MDeg_1norm;sum(temp_deg1norm)];
    fda_MNodeBetw_1norm=[fda_MNodeBetw_1norm;sum(temp_nodeb1norm)];
    fda_MLocEff_1norm=[fda_MLocEff_1norm;sum(temp_leff1norm)];
    
    fda_MClust_1=[fda_MClust_1;sum(temp_clust1)];
    fda_MDeg_1=[fda_MDeg_1;sum(temp_deg1)];
    fda_MNodeBetw_1=[fda_MNodeBetw_1;sum(temp_nodeb1)];
    fda_MLocEff_1=[fda_MLocEff_1;sum(temp_leff1)];
end

save(['Indiv_NetMesReg_' Group1 '.mat'],'MClust_1','MDeg_1','MNodeBetw_1','MLocEff_1','MPartCoeff_1','MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm','MPartCoeff_1norm','-v7');
save(['Indiv_AUC_NetMesReg_' Group1 '.mat'],'AUC_MClust_1','AUC_MDeg_1','AUC_MNodeBetw_1','AUC_MLocEff_1','AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm','-v7');
save(['Indiv_fda_NetMesReg_' Group1 '.mat'],'fda_MClust_1','fda_MDeg_1','fda_MNodeBetw_1','fda_MLocEff_1','fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm','-v7');

MClust_1norm=mean(MClust_1norm);
MDeg_1norm=mean(MDeg_1norm);
MNodeBetw_1norm=mean(MNodeBetw_1norm);
MLocEff_1norm=mean(MLocEff_1norm);
MPartCoeff_1norm=mean(MPartCoeff_1norm);

AUC_MClust_1norm=mean(AUC_MClust_1norm);
AUC_MDeg_1norm=mean(AUC_MDeg_1norm);
AUC_MNodeBetw_1norm=mean(AUC_MNodeBetw_1norm);
AUC_MLocEff_1norm=mean(AUC_MLocEff_1norm);

fda_MClust_1norm=mean(fda_MClust_1norm);
fda_MDeg_1norm=mean(fda_MDeg_1norm);
fda_MNodeBetw_1norm=mean(fda_MNodeBetw_1norm);
fda_MLocEff_1norm=mean(fda_MLocEff_1norm);

MClust_1=mean(MClust_1);
MDeg_1=mean(MDeg_1);
MNodeBetw_1=mean(MNodeBetw_1);
MLocEff_1=mean(MLocEff_1);
MPartCoeff_1=mean(MPartCoeff_1);

AUC_MClust_1=mean(AUC_MClust_1);
AUC_MDeg_1=mean(AUC_MDeg_1);
AUC_MNodeBetw_1=mean(AUC_MNodeBetw_1);
AUC_MLocEff_1=mean(AUC_MLocEff_1);

fda_MClust_1=mean(fda_MClust_1);
fda_MDeg_1=mean(fda_MDeg_1);
fda_MNodeBetw_1=mean(fda_MNodeBetw_1);
fda_MLocEff_1=mean(fda_MLocEff_1);

for i=1:size(NetMes2,1)
    fprintf('%-4s\n',['calculating group2 subject ' num2str(i) ' regional network measures....']);
    temp_clust2norm=[];temp_deg2norm=[];temp_nodeb2norm=[];temp_leff2norm=[]; temp_partcoeff2norm=[];
    temp_clust2=[];temp_deg2=[];temp_nodeb2=[];temp_leff2=[]; temp_partcoeff2=[];
    
    for j=1:size(NetMes2,2)
        temp_clust2norm=[temp_clust2norm;NetMes2{i,j}{7,3}'];
        temp_deg2norm=[temp_deg2norm;NetMes2{i,j}{1,3}];
        temp_nodeb2norm=[temp_nodeb2norm;NetMes2{i,j}{16,3}];
        temp_leff2norm=[temp_leff2norm;NetMes2{i,j}{11,3}'];
        temp_partcoeff2norm=[temp_partcoeff2norm;NetMes2{i,j}{27,3}'];
        
        temp_clust2=[temp_clust2;NetMes2{i,j}{7,1}'];
        temp_deg2=[temp_deg2;NetMes2{i,j}{1,1}];
        temp_nodeb2=[temp_nodeb2;NetMes2{i,j}{16,1}];
        temp_leff2=[temp_leff2;NetMes2{i,j}{11,1}'];
        temp_partcoeff2=[temp_partcoeff2;NetMes2{i,j}{27,1}'];
    end
    
    MClust_2norm=[MClust_2norm;mean(temp_clust2norm)];
    MDeg_2norm=[MDeg_2norm;mean(temp_deg2norm)];
    MNodeBetw_2norm=[MNodeBetw_2norm;mean(temp_nodeb2norm)];
    MLocEff_2norm=[MLocEff_2norm;mean(temp_leff2norm)];
    MPartCoeff_2norm=[MPartCoeff_2norm;mean(temp_partcoeff2norm)];
    
    MClust_2=[MClust_2;mean(temp_clust2)];
    MDeg_2=[MDeg_2;mean(temp_deg2)];
    MNodeBetw_2=[MNodeBetw_2;mean(temp_nodeb2)];
    MLocEff_2=[MLocEff_2;mean(temp_leff2)];
    MPartCoeff_2=[MPartCoeff_2;mean(temp_partcoeff2)];
    
    AUC_MClust_2norm=[AUC_MClust_2norm;trapz(Xax,temp_clust2norm)];
    AUC_MDeg_2norm=[AUC_MDeg_2norm;trapz(Xax,temp_deg2norm)];
    AUC_MNodeBetw_2norm=[AUC_MNodeBetw_2norm;trapz(Xax,temp_nodeb2norm)];
    AUC_MLocEff_2norm=[AUC_MLocEff_2norm;trapz(Xax,temp_leff2norm)];
    
    AUC_MClust_2=[AUC_MClust_2;trapz(Xax,temp_clust2)];
    AUC_MDeg_2=[AUC_MDeg_2;trapz(Xax,temp_deg2)];
    AUC_MNodeBetw_2=[AUC_MNodeBetw_2;trapz(Xax,temp_nodeb2)];
    AUC_MLocEff_2=[AUC_MLocEff_2;trapz(Xax,temp_leff2)];
    
    fda_MClust_2norm=[fda_MClust_2norm;sum(temp_clust2norm)];
    fda_MDeg_2norm=[fda_MDeg_2norm;sum(temp_deg2norm)];
    fda_MNodeBetw_2norm=[fda_MNodeBetw_2norm;sum(temp_nodeb2norm)];
    fda_MLocEff_2norm=[fda_MLocEff_2norm;sum(temp_leff2norm)];
    
    fda_MClust_2=[fda_MClust_2;sum(temp_clust2)];
    fda_MDeg_2=[fda_MDeg_2;sum(temp_deg2)];
    fda_MNodeBetw_2=[fda_MNodeBetw_2;sum(temp_nodeb2)];
    fda_MLocEff_2=[fda_MLocEff_2;sum(temp_leff2)];
end

save(['Indiv_NetMesReg_' Group2 '.mat'],'MClust_2','MDeg_2','MNodeBetw_2','MLocEff_2','MPartCoeff_2','MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm','MPartCoeff_2norm','-v7');
save(['Indiv_AUC_NetMesReg_' Group2 '.mat'],'AUC_MClust_2','AUC_MDeg_2','AUC_MNodeBetw_2','AUC_MLocEff_2','AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm','-v7');
save(['Indiv_fda_NetMesReg_' Group2 '.mat'],'fda_MClust_2','fda_MDeg_2','fda_MNodeBetw_2','fda_MLocEff_2','fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm','-v7');

MClust_2norm=mean(MClust_2norm);
MDeg_2norm=mean(MDeg_2norm);
MNodeBetw_2norm=mean(MNodeBetw_2norm);
MLocEff_2norm=mean(MLocEff_2norm);
MPartCoeff_2norm=mean(MPartCoeff_2norm);

AUC_MClust_2norm=mean(AUC_MClust_2norm);
AUC_MDeg_2norm=mean(AUC_MDeg_2norm);
AUC_MNodeBetw_2norm=mean(AUC_MNodeBetw_2norm);
AUC_MLocEff_2norm=mean(AUC_MLocEff_2norm);

fda_MClust_2norm=mean(fda_MClust_2norm);
fda_MDeg_2norm=mean(fda_MDeg_2norm);
fda_MNodeBetw_2norm=mean(fda_MNodeBetw_2norm);
fda_MLocEff_2norm=mean(fda_MLocEff_2norm);

MClust_2=mean(MClust_2);
MDeg_2=mean(MDeg_2);
MNodeBetw_2=mean(MNodeBetw_2);
MLocEff_2=mean(MLocEff_2);
MPartCoeff_2=mean(MPartCoeff_2);

AUC_MClust_2=mean(AUC_MClust_2);
AUC_MDeg_2=mean(AUC_MDeg_2);
AUC_MNodeBetw_2=mean(AUC_MNodeBetw_2);
AUC_MLocEff_2=mean(AUC_MLocEff_2);

fda_MClust_2=mean(fda_MClust_2);
fda_MDeg_2=mean(fda_MDeg_2);
fda_MNodeBetw_2=mean(fda_MNodeBetw_2);
fda_MLocEff_2=mean(fda_MLocEff_2);

MClust_1_randnorm=[];MDeg_1_randnorm=[];MNodeBetw_1_randnorm=[];MLocEff_1_randnorm=[];
MClust_2_randnorm=[];MDeg_2_randnorm=[];MNodeBetw_2_randnorm=[];MLocEff_2_randnorm=[];

AUC_MClust_1_randnorm=[];AUC_MDeg_1_randnorm=[];AUC_MNodeBetw_1_randnorm=[];AUC_MLocEff_1_randnorm=[];
AUC_MClust_2_randnorm=[];AUC_MDeg_2_randnorm=[];AUC_MNodeBetw_2_randnorm=[];AUC_MLocEff_2_randnorm=[];

fda_MClust_1_randnorm=[];fda_MDeg_1_randnorm=[];fda_MNodeBetw_1_randnorm=[];fda_MLocEff_1_randnorm=[];
fda_MClust_2_randnorm=[];fda_MDeg_2_randnorm=[];fda_MNodeBetw_2_randnorm=[];fda_MLocEff_2_randnorm=[];

MClust_1_rand=[];MDeg_1_rand=[];MNodeBetw_1_rand=[];MLocEff_1_rand=[];
MClust_2_rand=[];MDeg_2_rand=[];MNodeBetw_2_rand=[];MLocEff_2_rand=[];

AUC_MClust_1_rand=[];AUC_MDeg_1_rand=[];AUC_MNodeBetw_1_rand=[];AUC_MLocEff_1_rand=[];
AUC_MClust_2_rand=[];AUC_MDeg_2_rand=[];AUC_MNodeBetw_2_rand=[];AUC_MLocEff_2_rand=[];

fda_MClust_1_rand=[];fda_MDeg_1_rand=[];fda_MNodeBetw_1_rand=[];fda_MLocEff_1_rand=[];
fda_MClust_2_rand=[];fda_MDeg_2_rand=[];fda_MNodeBetw_2_rand=[];fda_MLocEff_2_rand=[];

for k=1:size(NetMes1_rand,2)
    fprintf('%-4s\n',['calculating random net ' num2str(k) ' regional network measures....']);
    temp_rand1norm=NetMes1_rand{1,k};
    temp_rand2norm=NetMes2_rand{1,k};
    
    temp_rand1norm=temp_rand1norm(:,MinIdx:MaxIdx);
    temp_rand2norm=temp_rand2norm(:,MinIdx:MaxIdx);
    
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
    
    MClust_1_rand=[];MDeg_1_rand=[];MNodeBetw_1_rand=[];MLocEff_1_rand=[];
    MClust_2_rand=[];MDeg_2_rand=[];MNodeBetw_2_rand=[];MLocEff_2_rand=[];
    AUC_MClust_1_rand=[];AUC_MDeg_1_rand=[];AUC_MNodeBetw_1_rand=[];AUC_MLocEff_1_rand=[];
    AUC_MClust_2_rand=[];AUC_MDeg_2_rand=[];AUC_MNodeBetw_2_rand=[];AUC_MLocEff_2_rand=[];
    fda_MClust_1_rand=[];fda_MDeg_1_rand=[];fda_MNodeBetw_1_rand=[];fda_MLocEff_1_rand=[];
    fda_MClust_2_rand=[];fda_MDeg_2_rand=[];fda_MNodeBetw_2_rand=[];fda_MLocEff_2_rand=[];
    
    for i=1:size(temp_rand1norm,1)
        temp_rand_clust1norm=[];temp_rand_deg1norm=[];temp_rand_nodeb1norm=[];temp_rand_leff1norm=[];
        
        for j=1:size(temp_rand1norm,2)
            temp_rand_clust1norm=[temp_rand_clust1norm;temp_rand1norm{i,j}{7,3}'];
            temp_rand_deg1norm=[temp_rand_deg1norm;temp_rand1norm{i,j}{1,3}];
            temp_rand_nodeb1norm=[temp_rand_nodeb1norm;temp_rand1norm{i,j}{16,3}];
            temp_rand_leff1norm=[temp_rand_leff1norm;temp_rand1norm{i,j}{11,3}'];
        end
        
        MClust_1norm_rand=[MClust_1norm_rand;mean(temp_rand_clust1norm)];
        MDeg_1norm_rand=[MDeg_1norm_rand;mean(temp_rand_deg1norm)];
        MNodeBetw_1norm_rand=[MNodeBetw_1norm_rand;mean(temp_rand_nodeb1norm)];
        MLocEff_1norm_rand=[MLocEff_1norm_rand;mean(temp_rand_leff1norm)];
        
        AUC_MClust_1norm_rand=[AUC_MClust_1norm_rand;trapz(Xax,temp_rand_clust1norm)];
        AUC_MDeg_1norm_rand=[AUC_MDeg_1norm_rand;trapz(Xax,temp_rand_deg1norm)];
        AUC_MNodeBetw_1norm_rand=[AUC_MNodeBetw_1norm_rand;trapz(Xax,temp_rand_nodeb1norm)];
        AUC_MLocEff_1norm_rand=[AUC_MLocEff_1norm_rand;trapz(Xax,temp_rand_leff1norm)];
        
        fda_MClust_1norm_rand=[fda_MClust_1norm_rand;sum(temp_rand_clust1norm)];
        fda_MDeg_1norm_rand=[fda_MDeg_1norm_rand;sum(temp_rand_deg1norm)];
        fda_MNodeBetw_1norm_rand=[fda_MNodeBetw_1norm_rand;sum(temp_rand_nodeb1norm)];
        fda_MLocEff_1norm_rand=[fda_MLocEff_1norm_rand;sum(temp_rand_leff1norm)];
    end
    
    for i=1:size(temp_rand1,1)
        temp_rand_clust1=[];temp_rand_deg1=[];temp_rand_nodeb1=[];temp_rand_leff1=[];
        
        for j=1:size(temp_rand1,2)
            temp_rand_clust1=[temp_rand_clust1;temp_rand1{i,j}{7,3}'];
            temp_rand_deg1=[temp_rand_deg1;temp_rand1{i,j}{1,3}];
            temp_rand_nodeb1=[temp_rand_nodeb1;temp_rand1{i,j}{16,3}];
            temp_rand_leff1=[temp_rand_leff1;temp_rand1{i,j}{11,3}'];
        end
        
        MClust_1_rand=[MClust_1_rand;mean(temp_rand_clust1)];
        MDeg_1_rand=[MDeg_1_rand;mean(temp_rand_deg1)];
        MNodeBetw_1_rand=[MNodeBetw_1_rand;mean(temp_rand_nodeb1)];
        MLocEff_1_rand=[MLocEff_1_rand;mean(temp_rand_leff1)];
        
        AUC_MClust_1_rand=[AUC_MClust_1_rand;trapz(Xax,temp_rand_clust1)];
        AUC_MDeg_1_rand=[AUC_MDeg_1_rand;trapz(Xax,temp_rand_deg1)];
        AUC_MNodeBetw_1_rand=[AUC_MNodeBetw_1_rand;trapz(Xax,temp_rand_nodeb1)];
        AUC_MLocEff_1_rand=[AUC_MLocEff_1_rand;trapz(Xax,temp_rand_leff1)];
        
        fda_MClust_1_rand=[fda_MClust_1_rand;sum(temp_rand_clust1)];
        fda_MDeg_1_rand=[fda_MDeg_1_rand;sum(temp_rand_deg1)];
        fda_MNodeBetw_1_rand=[fda_MNodeBetw_1_rand;sum(temp_rand_nodeb1)];
        fda_MLocEff_1_rand=[fda_MLocEff_1_rand;sum(temp_rand_leff1)];
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
    
    MClust_1_rand=[MClust_1_rand;mean(MClust_1_rand)];
    MDeg_1_rand=[MDeg_1_rand;mean(MDeg_1_rand)];
    MNodeBetw_1_rand=[MNodeBetw_1_rand;mean(MNodeBetw_1_rand)];
    MLocEff_1_rand=[MLocEff_1_rand;mean(MLocEff_1_rand)];
    
    AUC_MClust_1_rand=[AUC_MClust_1_rand;mean(AUC_MClust_1_rand)];
    AUC_MDeg_1_rand=[AUC_MDeg_1_rand;mean(AUC_MDeg_1_rand)];
    AUC_MNodeBetw_1_rand=[AUC_MNodeBetw_1_rand;mean(AUC_MNodeBetw_1_rand)];
    AUC_MLocEff_1_rand=[AUC_MLocEff_1_rand;mean(AUC_MLocEff_1_rand)];
    
    fda_MClust_1_rand=[fda_MClust_1_rand;mean(fda_MClust_1_rand)];
    fda_MDeg_1_rand=[fda_MDeg_1_rand;mean(fda_MDeg_1_rand)];
    fda_MNodeBetw_1_rand=[fda_MNodeBetw_1_rand;mean(fda_MNodeBetw_1_rand)];
    fda_MLocEff_1_rand=[fda_MLocEff_1_rand;mean(fda_MLocEff_1_rand)];
    
    for i=1:size(temp_rand2norm,1)
        temp_rand_clust2norm=[];temp_rand_deg2norm=[];temp_rand_nodeb2norm=[];temp_rand_leff2norm=[];
        
        for j=1:size(temp_rand2norm,2)
            temp_rand_clust2norm=[temp_rand_clust2norm;temp_rand2norm{i,j}{7,3}'];
            temp_rand_deg2norm=[temp_rand_deg2norm;temp_rand2norm{i,j}{1,3}];
            temp_rand_nodeb2norm=[temp_rand_nodeb2norm;temp_rand2norm{i,j}{16,3}];
            temp_rand_leff2norm=[temp_rand_leff2norm;temp_rand2norm{i,j}{11,3}'];
        end
        
        MClust_2norm_rand=[MClust_2norm_rand;mean(temp_rand_clust2norm)];
        MDeg_2norm_rand=[MDeg_2norm_rand;mean(temp_rand_deg2norm)];
        MNodeBetw_2norm_rand=[MNodeBetw_2norm_rand;mean(temp_rand_nodeb2norm)];
        MLocEff_2norm_rand=[MLocEff_2norm_rand;mean(temp_rand_leff2norm)];
        
        AUC_MClust_2norm_rand=[AUC_MClust_2norm_rand;trapz(Xax,temp_rand_clust2norm)];
        AUC_MDeg_2norm_rand=[AUC_MDeg_2norm_rand;trapz(Xax,temp_rand_deg2norm)];
        AUC_MNodeBetw_2norm_rand=[AUC_MNodeBetw_2norm_rand;trapz(Xax,temp_rand_nodeb2norm)];
        AUC_MLocEff_2norm_rand=[AUC_MLocEff_2norm_rand;trapz(Xax,temp_rand_leff2norm)];
        
        fda_MClust_2norm_rand=[fda_MClust_2norm_rand;sum(temp_rand_clust2norm)];
        fda_MDeg_2norm_rand=[fda_MDeg_2norm_rand;sum(temp_rand_deg2norm)];
        fda_MNodeBetw_2norm_rand=[fda_MNodeBetw_2norm_rand;sum(temp_rand_nodeb2norm)];
        fda_MLocEff_2norm_rand=[fda_MLocEff_2norm_rand;sum(temp_rand_leff2norm)];
    end
    
    for i=1:size(temp_rand2,1)
        temp_rand_clust2=[];temp_rand_deg2=[];temp_rand_nodeb2=[];temp_rand_leff2=[];
        
        for j=1:size(temp_rand2,2)
            temp_rand_clust2=[temp_rand_clust2;temp_rand2{i,j}{7,3}'];
            temp_rand_deg2=[temp_rand_deg2;temp_rand2{i,j}{1,3}];
            temp_rand_nodeb2=[temp_rand_nodeb2;temp_rand2{i,j}{16,3}];
            temp_rand_leff2=[temp_rand_leff2;temp_rand2{i,j}{11,3}'];
        end
        
        MClust_2_rand=[MClust_2_rand;mean(temp_rand_clust2)];
        MDeg_2_rand=[MDeg_2_rand;mean(temp_rand_deg2)];
        MNodeBetw_2_rand=[MNodeBetw_2_rand;mean(temp_rand_nodeb2)];
        MLocEff_2_rand=[MLocEff_2_rand;mean(temp_rand_leff2)];
        
        AUC_MClust_2_rand=[AUC_MClust_2_rand;trapz(Xax,temp_rand_clust2)];
        AUC_MDeg_2_rand=[AUC_MDeg_2_rand;trapz(Xax,temp_rand_deg2)];
        AUC_MNodeBetw_2_rand=[AUC_MNodeBetw_2_rand;trapz(Xax,temp_rand_nodeb2)];
        AUC_MLocEff_2_rand=[AUC_MLocEff_2_rand;trapz(Xax,temp_rand_leff2)];
        
        fda_MClust_2_rand=[fda_MClust_2_rand;sum(temp_rand_clust2)];
        fda_MDeg_2_rand=[fda_MDeg_2_rand;sum(temp_rand_deg2)];
        fda_MNodeBetw_2_rand=[fda_MNodeBetw_2_rand;sum(temp_rand_nodeb2)];
        fda_MLocEff_2_rand=[fda_MLocEff_2_rand;sum(temp_rand_leff2)];
    end
    
    MClust_2_rand=[MClust_2_rand;mean(MClust_2_rand)];
    MDeg_2_rand=[MDeg_2_rand;mean(MDeg_2_rand)];
    MNodeBetw_2_rand=[MNodeBetw_2_rand;mean(MNodeBetw_2_rand)];
    MLocEff_2_rand=[MLocEff_2_rand;mean(MLocEff_2_rand)];
    
    AUC_MClust_2_rand=[AUC_MClust_2_rand;mean(AUC_MClust_2_rand)];
    AUC_MDeg_2_rand=[AUC_MDeg_2_rand;mean(AUC_MDeg_2_rand)];
    AUC_MNodeBetw_2_rand=[AUC_MNodeBetw_2_rand;mean(AUC_MNodeBetw_2_rand)];
    AUC_MLocEff_2_rand=[AUC_MLocEff_2_rand;mean(AUC_MLocEff_2_rand)];
    
    fda_MClust_2_rand=[fda_MClust_2_rand;mean(fda_MClust_2_rand)];
    fda_MDeg_2_rand=[fda_MDeg_2_rand;mean(fda_MDeg_2_rand)];
    fda_MNodeBetw_2_rand=[fda_MNodeBetw_2_rand;mean(fda_MNodeBetw_2_rand)];
    fda_MLocEff_2_rand=[fda_MLocEff_2_rand;mean(fda_MLocEff_2_rand)];
end

nROI=size(MClust_1norm,2);

save(['NetMesReg_rand_' Group1 '.mat'],'MClust_1_rand','MDeg_1_rand','MNodeBetw_1_rand','MLocEff_1_rand','MClust_1_randnorm','MDeg_1_randnorm','MNodeBetw_1_randnorm','MLocEff_1_randnorm','-v7');
save(['NetMesReg_rand_' Group2 '.mat'],'MClust_2_rand','MDeg_2_rand','MNodeBetw_2_rand','MLocEff_2_rand','MClust_2_randnorm','MDeg_2_randnorm','MNodeBetw_2_randnorm','MLocEff_2_randnorm','-v7');
save(['AUC_NetMesReg_rand_' Group1 '.mat'],'AUC_MClust_1_rand','AUC_MDeg_1_rand','AUC_MNodeBetw_1_rand','AUC_MLocEff_1_rand','AUC_MClust_1_randnorm','AUC_MDeg_1_randnorm','AUC_MNodeBetw_1_randnorm','AUC_MLocEff_1_randnorm','-v7');
save(['AUC_NetMesReg_rand_' Group2 '.mat'],'AUC_MClust_2_rand','AUC_MDeg_2_rand','AUC_MNodeBetw_2_rand','AUC_MLocEff_2_rand','AUC_MClust_2_randnorm','AUC_MDeg_2_randnorm','AUC_MNodeBetw_2_randnorm','AUC_MLocEff_2_randnorm','-v7');
save(['fda_NetMesReg_rand_' Group1 '.mat'],'fda_MClust_1_rand','fda_MDeg_1_rand','fda_MNodeBetw_1_rand','fda_MLocEff_1_rand','fda_MClust_1_randnorm','fda_MDeg_1_randnorm','fda_MNodeBetw_1_randnorm','fda_MLocEff_1_randnorm','-v7');
save(['fda_NetMesReg_rand_' Group2 '.mat'],'fda_MClust_2_rand','fda_MDeg_2_rand','fda_MNodeBetw_2_rand','fda_MLocEff_2_rand','fda_MClust_2_randnorm','fda_MDeg_2_randnorm','fda_MNodeBetw_2_randnorm','fda_MLocEff_2_randnorm','-v7');

save(['NetMesReg_' Group1 '.mat'],'MClust_1','MDeg_1','MNodeBetw_1','MLocEff_1','MPartCoeff_1','MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm','MPartCoeff_1norm','-v7');
save(['NetMesReg_' Group2 '.mat'],'MClust_2','MDeg_2','MNodeBetw_2','MLocEff_2','MPartCoeff_2','MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm','MPartCoeff_2norm','-v7');
save(['AUC_NetMesReg_' Group1 '.mat'],'AUC_MClust_1','AUC_MDeg_1','AUC_MNodeBetw_1','AUC_MLocEff_1','AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm','-v7');
save(['AUC_NetMesReg_' Group2 '.mat'],'AUC_MClust_2','AUC_MDeg_2','AUC_MNodeBetw_2','AUC_MLocEff_2','AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm','-v7');
save(['fda_NetMesReg_' Group1 '.mat'],'fda_MClust_1','fda_MDeg_1','fda_MNodeBetw_1','fda_MLocEff_1','fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm','-v7');
save(['fda_NetMesReg_' Group2 '.mat'],'fda_MClust_2','fda_MDeg_2','fda_MNodeBetw_2','fda_MLocEff_2','fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm','-v7');

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

fprintf('%-4s\n','.... done ....');