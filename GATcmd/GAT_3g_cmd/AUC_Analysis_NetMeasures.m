function AUC_Analysis_NetMeasures(data_log,NetM_1_2,NetM_2_3,NetM_1_3,MinThr,MaxThr)

%%
%-- updated on August 21,2011 for functionals
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)

%% Main Body of the script common to all group combinations

ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2; Group3 = ff.mat4GAT.g3;

if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha')
    Alpha = ff.mat4GAT.alpha;
    Tail = ff.mat4GAT.tail;
else
    Alpha = .05;
    Tail = 2;
end

fprintf('%-4s\n',['loading inputs...']);

load(NetM_1_2);
MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinIdx=find(single(xxx)==single(MinThr));
MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(Dens_1_rand,2)
    MaxIdx = size(Dens_1_rand,2);
end

%% Group1 vs Group2

load(NetM_1_2);
Dens_1_rand=Dens_1_rand(:,MinIdx:MaxIdx,:);
MClust_1_rand=MClust_1_rand(:,MinIdx:MaxIdx,:);
MDeg_1_rand=MDeg_1_rand(:,MinIdx:MaxIdx,:);
Trans_1_rand=Trans_1_rand(:,MinIdx:MaxIdx,:);
Assort_1_rand=Assort_1_rand(:,MinIdx:MaxIdx,:);
GEff_1_rand=GEff_1_rand(:,MinIdx:MaxIdx,:);
MLocEff_1_rand=MLocEff_1_rand(:,MinIdx:MaxIdx,:);
Mod_1_rand=Mod_1_rand(:,MinIdx:MaxIdx,:);
ModL_1_rand=ModL_1_rand(:,MinIdx:MaxIdx,:);
PathL_1_rand=PathL_1_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_1_rand=MNodeBetw_1_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_1_rand=MEdgeBetw_1_rand(:,MinIdx:MaxIdx,:);
Lambda_1_rand=Lambda_1_rand(:,MinIdx:MaxIdx,:);
Gamma_1_rand=Gamma_1_rand(:,MinIdx:MaxIdx,:);
Sigma_1_rand=Sigma_1_rand(:,MinIdx:MaxIdx,:);

Dens_2_rand=Dens_2_rand(:,MinIdx:MaxIdx,:);
MClust_2_rand=MClust_2_rand(:,MinIdx:MaxIdx,:);
MDeg_2_rand=MDeg_2_rand(:,MinIdx:MaxIdx,:);
Trans_2_rand=Trans_2_rand(:,MinIdx:MaxIdx,:);
Assort_2_rand=Assort_2_rand(:,MinIdx:MaxIdx,:);
GEff_2_rand=GEff_2_rand(:,MinIdx:MaxIdx,:);
MLocEff_2_rand=MLocEff_2_rand(:,MinIdx:MaxIdx,:);
Mod_2_rand=Mod_2_rand(:,MinIdx:MaxIdx,:);
ModL_2_rand=ModL_2_rand(:,MinIdx:MaxIdx,:);
PathL_2_rand=PathL_2_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_2_rand=MNodeBetw_2_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_2_rand=MEdgeBetw_2_rand(:,MinIdx:MaxIdx,:);
Lambda_2_rand=Lambda_2_rand(:,MinIdx:MaxIdx,:);
Gamma_2_rand=Gamma_2_rand(:,MinIdx:MaxIdx,:);
Sigma_2_rand=Sigma_2_rand(:,MinIdx:MaxIdx,:);

Dens_1=Dens_1(:,MinIdx:MaxIdx);
MClust_1=MClust_1(:,MinIdx:MaxIdx);
MDeg_1=MDeg_1(:,MinIdx:MaxIdx);
Trans_1=Trans_1(:,MinIdx:MaxIdx);
Assort_1=Assort_1(:,MinIdx:MaxIdx);
GEff_1=GEff_1(:,MinIdx:MaxIdx);
MLocEff_1=MLocEff_1(:,MinIdx:MaxIdx);
Mod_1=Mod_1(:,MinIdx:MaxIdx);
ModL_1=ModL_1(:,MinIdx:MaxIdx);
PathL_1=PathL_1(:,MinIdx:MaxIdx);
MNodeBetw_1=MNodeBetw_1(:,MinIdx:MaxIdx);
MEdgeBetw_1=MEdgeBetw_1(:,MinIdx:MaxIdx);
Lambda_1=Lambda_1(:,MinIdx:MaxIdx);
Gamma_1=Gamma_1(:,MinIdx:MaxIdx);
Sigma_1=Sigma_1(:,MinIdx:MaxIdx);

Dens_2=Dens_2(:,MinIdx:MaxIdx);
MClust_2=MClust_2(:,MinIdx:MaxIdx);
MDeg_2=MDeg_2(:,MinIdx:MaxIdx);
Trans_2=Trans_2(:,MinIdx:MaxIdx);
Assort_2=Assort_2(:,MinIdx:MaxIdx);
GEff_2=GEff_2(:,MinIdx:MaxIdx);
MLocEff_2=MLocEff_2(:,MinIdx:MaxIdx);
Mod_2=Mod_2(:,MinIdx:MaxIdx);
ModL_2=ModL_2(:,MinIdx:MaxIdx);
PathL_2=PathL_2(:,MinIdx:MaxIdx);
MNodeBetw_2=MNodeBetw_2(:,MinIdx:MaxIdx);
MEdgeBetw_2=MEdgeBetw_2(:,MinIdx:MaxIdx);
Lambda_2=Lambda_2(:,MinIdx:MaxIdx);
Gamma_2=Gamma_2(:,MinIdx:MaxIdx);
Sigma_2=Sigma_2(:,MinIdx:MaxIdx);

fprintf('%-4s\n',['saving net measures across new threshold range...']);

dd=pwd;
mkdir('AUC_Results/AUC_Results_G1_vs_G2');
cd([dd '/AUC_Results/AUC_Results_G1_vs_G2']);

save NetMesPerDens_AveAcrossDens Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
    Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
    PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
    Sigma_2_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
    Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ModL_2 PathL_2 MNodeBetw_2 ...
    MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 ...
    
Xax=MinMesPlot:MesStepPlot:MaxMesPlot;
iXax=Xax(MinIdx:MaxIdx);

auc_MClust_1 = trapz(iXax,MClust_1,2);
auc_MDeg_1 = trapz(iXax,MDeg_1,2);
auc_Trans_1 = trapz(iXax,Trans_1,2);
auc_Assort_1 = trapz(iXax,Assort_1,2);
auc_GEff_1 = trapz(iXax,GEff_1,2);
auc_MLocEff_1 = trapz(iXax,MLocEff_1,2);
auc_Mod_1 = trapz(iXax,Mod_1,2);
auc_ModL_1 = trapz(iXax,ModL_1,2);
auc_PathL_1 = trapz(iXax,PathL_1,2);
auc_MNodeBetw_1 = trapz(iXax,MNodeBetw_1,2);
auc_MEdgeBetw_1 = trapz(iXax,MEdgeBetw_1,2);
auc_Lambda_1 = trapz(iXax,Lambda_1,2);
auc_Gamma_1 = trapz(iXax,Gamma_1,2);
auc_Sigma_1 = trapz(iXax,Sigma_1,2);

auc_MClust_2 = trapz(iXax,MClust_2,2);
auc_MDeg_2 = trapz(iXax,MDeg_2,2);
auc_Trans_2 = trapz(iXax,Trans_2,2);
auc_Assort_2 = trapz(iXax,Assort_2,2);
auc_GEff_2 = trapz(iXax,GEff_2,2);
auc_MLocEff_2 = trapz(iXax,MLocEff_2,2);
auc_Mod_2 = trapz(iXax,Mod_2,2);
auc_ModL_2 = trapz(iXax,ModL_2,2);
auc_PathL_2 = trapz(iXax,PathL_2,2);
auc_MNodeBetw_2 = trapz(iXax,MNodeBetw_2,2);
auc_MEdgeBetw_2 = trapz(iXax,MEdgeBetw_2,2);
auc_Lambda_2 = trapz(iXax,Lambda_2,2);
auc_Gamma_2 = trapz(iXax,Gamma_2,2);
auc_Sigma_2 = trapz(iXax,Sigma_2,2);

auc_MClust_1_rand = trapz(iXax,MClust_1_rand,2);
auc_MDeg_1_rand = trapz(iXax,MDeg_1_rand,2);
auc_Trans_1_rand = trapz(iXax,Trans_1_rand,2);
auc_Assort_1_rand = trapz(iXax,Assort_1_rand,2);
auc_GEff_1_rand = trapz(iXax,GEff_1_rand,2);
auc_MLocEff_1_rand = trapz(iXax,MLocEff_1_rand,2);
auc_Mod_1_rand = trapz(iXax,Mod_1_rand,2);
auc_ModL_1_rand = trapz(iXax,ModL_1_rand,2);
auc_PathL_1_rand = trapz(iXax,PathL_1_rand,2);
auc_MNodeBetw_1_rand = trapz(iXax,MNodeBetw_1_rand,2);
auc_MEdgeBetw_1_rand = trapz(iXax,MEdgeBetw_1_rand,2);
auc_Lambda_1_rand = trapz(iXax,Lambda_1_rand,2);
auc_Gamma_1_rand = trapz(iXax,Gamma_1_rand,2);
auc_Sigma_1_rand = trapz(iXax,Sigma_1_rand,2);

auc_MClust_2_rand = trapz(iXax,MClust_2_rand,2);
auc_MDeg_2_rand = trapz(iXax,MDeg_2_rand,2);
auc_Trans_2_rand = trapz(iXax,Trans_2_rand,2);
auc_Assort_2_rand = trapz(iXax,Assort_2_rand,2);
auc_GEff_2_rand = trapz(iXax,GEff_2_rand,2);
auc_MLocEff_2_rand = trapz(iXax,MLocEff_2_rand,2);
auc_Mod_2_rand = trapz(iXax,Mod_2_rand,2);
auc_ModL_2_rand = trapz(iXax,ModL_2_rand,2);
auc_PathL_2_rand = trapz(iXax,PathL_2_rand,2);
auc_MNodeBetw_2_rand = trapz(iXax,MNodeBetw_2_rand,2);
auc_MEdgeBetw_2_rand = trapz(iXax,MEdgeBetw_2_rand,2);
auc_Lambda_2_rand = trapz(iXax,Lambda_2_rand,2);
auc_Gamma_2_rand = trapz(iXax,Gamma_2_rand,2);
auc_Sigma_2_rand = trapz(iXax,Sigma_2_rand,2);

save AUC_NetMesPerDens_AveAcrossDens  auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_MClust_2_rand auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
    auc_MClust_2 auc_MDeg_2 auc_Trans_2 auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2


AUC_indiv_mat= num2cell([[1:ff.mat4GAT.n1 1:ff.mat4GAT.n2]' [ones(ff.mat4GAT.n1,1)*1; ones(ff.mat4GAT.n2,1)*2]  [auc_MClust_1; auc_MClust_2] ...
    [auc_MDeg_1; auc_MDeg_2] [auc_Trans_1; auc_Trans_2] [auc_Assort_1; auc_Assort_2] [auc_GEff_1; auc_GEff_2] [auc_MLocEff_1; auc_MLocEff_2] ...
    [auc_Mod_1; auc_Mod_2] [auc_PathL_1; auc_PathL_2] [auc_Lambda_1; auc_Lambda_2] [auc_Gamma_1; auc_Gamma_2] [auc_Sigma_1; auc_Sigma_2]...
    [auc_MEdgeBetw_1; auc_MEdgeBetw_2] [auc_MNodeBetw_1; auc_MNodeBetw_2]]);
colnames = {'Subj','Group','MClust','MDegree','Trans','Assort','Eglob','MEloc','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBtwn','MNodeBtwn'};
AUCdata = cell2table(AUC_indiv_mat,'VariableNames',colnames);
save AUCdata AUCdata
writetable(AUCdata,'AUCdata.xlsx')


nRand=size(auc_MClust_1_rand,3);

auc_MClust_1 = mean(auc_MClust_1);
auc_MDeg_1 = mean(auc_MDeg_1);
auc_Trans_1 = mean(auc_Trans_1);
auc_Assort_1 = mean(auc_Assort_1);
auc_GEff_1 = mean(auc_GEff_1);
auc_MLocEff_1 = mean(auc_MLocEff_1);
auc_Mod_1 = mean(auc_Mod_1);
auc_ModL_1 = mean(auc_ModL_1);
auc_PathL_1 = mean(auc_PathL_1);
auc_MNodeBetw_1 = mean(auc_MNodeBetw_1);
auc_MEdgeBetw_1 = mean(auc_MEdgeBetw_1);
auc_Lambda_1 = mean(auc_Lambda_1);
auc_Gamma_1 = mean(auc_Gamma_1);
auc_Sigma_1 = mean(auc_Sigma_1);

auc_MClust_2 = mean(auc_MClust_2);
auc_MDeg_2 = mean(auc_MDeg_2);
auc_Trans_2 = mean(auc_Trans_2);
auc_Assort_2 = mean(auc_Assort_2);
auc_GEff_2 = mean(auc_GEff_2);
auc_MLocEff_2 = mean(auc_MLocEff_2);
auc_Mod_2 = mean(auc_Mod_2);
auc_ModL_2 = mean(auc_ModL_2);
auc_PathL_2 = mean(auc_PathL_2);
auc_MNodeBetw_2 = mean(auc_MNodeBetw_2);
auc_MEdgeBetw_2 = mean(auc_MEdgeBetw_2);
auc_Lambda_2 = mean(auc_Lambda_2);
auc_Gamma_2 = mean(auc_Gamma_2);
auc_Sigma_2 = mean(auc_Sigma_2);

auc_MClust_1_rand = mean(auc_MClust_1_rand,1);auc_MClust_1_rand = reshape(auc_MClust_1_rand,1,nRand);
auc_MDeg_1_rand = mean(auc_MDeg_1_rand,1);auc_MDeg_1_rand = reshape(auc_MDeg_1_rand,1,nRand);
auc_Trans_1_rand = mean(auc_Trans_1_rand,1);auc_Trans_1_rand = reshape(auc_Trans_1_rand ,1,nRand);
auc_Assort_1_rand = mean(auc_Assort_1_rand,1);auc_Assort_1_rand = reshape(auc_Assort_1_rand ,1,nRand);
auc_GEff_1_rand = mean(auc_GEff_1_rand,1);auc_GEff_1_rand = reshape(auc_GEff_1_rand ,1,nRand);
auc_MLocEff_1_rand = mean(auc_MLocEff_1_rand,1);auc_MLocEff_1_rand = reshape(auc_MLocEff_1_rand ,1,nRand);
auc_Mod_1_rand = mean(auc_Mod_1_rand,1);auc_Mod_1_rand = reshape(auc_Mod_1_rand ,1,nRand);
auc_ModL_1_rand = mean(auc_ModL_1_rand,1);auc_ModL_1_rand = reshape(auc_ModL_1_rand ,1,nRand);
auc_PathL_1_rand = mean(auc_PathL_1_rand,1);auc_PathL_1_rand = reshape(auc_PathL_1_rand ,1,nRand);
auc_MNodeBetw_1_rand = mean(auc_MNodeBetw_1_rand,1);auc_MNodeBetw_1_rand = reshape(auc_MNodeBetw_1_rand ,1,nRand);
auc_MEdgeBetw_1_rand = mean(auc_MEdgeBetw_1_rand,1);auc_MEdgeBetw_1_rand = reshape(auc_MEdgeBetw_1_rand ,1,nRand);
auc_Lambda_1_rand = mean(auc_Lambda_1_rand,1);auc_Lambda_1_rand = reshape( auc_Lambda_1_rand,1,nRand);
auc_Gamma_1_rand = mean(auc_Gamma_1_rand,1);auc_Gamma_1_rand = reshape(auc_Gamma_1_rand ,1,nRand);
auc_Sigma_1_rand = mean(auc_Sigma_1_rand,1);auc_Sigma_1_rand = reshape(auc_Sigma_1_rand ,1,nRand);

auc_MClust_2_rand = mean(auc_MClust_2_rand,1);auc_MClust_2_rand = reshape(auc_MClust_2_rand,1,nRand);
auc_MDeg_2_rand = mean(auc_MDeg_2_rand,1);auc_MDeg_2_rand = reshape(auc_MDeg_2_rand,1,nRand);
auc_Trans_2_rand = mean(auc_Trans_2_rand,1);auc_Trans_2_rand = reshape(auc_Trans_2_rand ,1,nRand);
auc_Assort_2_rand = mean(auc_Assort_2_rand,1);auc_Assort_2_rand = reshape(auc_Assort_2_rand ,1,nRand);
auc_GEff_2_rand = mean(auc_GEff_2_rand,1);auc_GEff_2_rand = reshape(auc_GEff_2_rand ,1,nRand);
auc_MLocEff_2_rand = mean(auc_MLocEff_2_rand,1);auc_MLocEff_2_rand = reshape(auc_MLocEff_2_rand ,1,nRand);
auc_Mod_2_rand = mean(auc_Mod_2_rand,1);auc_Mod_2_rand = reshape(auc_Mod_2_rand ,1,nRand);
auc_ModL_2_rand = mean(auc_ModL_2_rand,1);auc_ModL_2_rand = reshape(auc_ModL_2_rand ,1,nRand);
auc_PathL_2_rand = mean(auc_PathL_2_rand,1);auc_PathL_2_rand = reshape(auc_PathL_2_rand ,1,nRand);
auc_MNodeBetw_2_rand = mean(auc_MNodeBetw_2_rand,1);auc_MNodeBetw_2_rand = reshape(auc_MNodeBetw_2_rand ,1,nRand);
auc_MEdgeBetw_2_rand = mean(auc_MEdgeBetw_2_rand,1);auc_MEdgeBetw_2_rand = reshape(auc_MEdgeBetw_2_rand ,1,nRand);
auc_Lambda_2_rand = mean(auc_Lambda_2_rand,1);auc_Lambda_2_rand = reshape( auc_Lambda_2_rand,1,nRand);
auc_Gamma_2_rand = mean(auc_Gamma_2_rand,1);auc_Gamma_2_rand = reshape(auc_Gamma_2_rand ,1,nRand);
auc_Sigma_2_rand = mean(auc_Sigma_2_rand,1);auc_Sigma_2_rand = reshape(auc_Sigma_2_rand ,1,nRand);

save AUC_Mean auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_MClust_2_rand  auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
    auc_MClust_2 auc_MDeg_2 auc_Trans_2 ...
    auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2

%Difference between Group2 and Group1
p = Alpha;

CL_auc_21_Geff = CL_per(auc_GEff_2_rand' - auc_GEff_1_rand', p);
CL_auc_21_Leff = CL_per(auc_MLocEff_2_rand' - auc_MLocEff_1_rand', p);
CL_auc_21_Trans = CL_per(auc_Trans_2_rand' - auc_Trans_1_rand', p);
CL_auc_21_Mod = CL_per(auc_Mod_2_rand' - auc_Mod_1_rand', p);
CL_auc_21_Lambda = CL_per(auc_Lambda_2_rand' - auc_Lambda_1_rand', p);
CL_auc_21_Gamma = CL_per(auc_Gamma_2_rand' - auc_Gamma_1_rand', p);
CL_auc_21_Sigma = CL_per(auc_Sigma_2_rand' - auc_Sigma_1_rand', p);
CL_auc_21_Assort = CL_per(auc_Assort_2_rand' - auc_Assort_1_rand', p);
CL_auc_21_PathL = CL_per(auc_PathL_2_rand' - auc_PathL_1_rand', p);
CL_auc_21_MClust = CL_per(auc_MClust_2_rand' - auc_MClust_1_rand', p);
CL_auc_21_MDeg = CL_per(auc_MDeg_2_rand' - auc_MDeg_1_rand', p);
CL_auc_21_MEdgeBetw = CL_per(auc_MEdgeBetw_2_rand' - auc_MEdgeBetw_1_rand', p);
CL_auc_21_MNodeBetw = CL_per(auc_MNodeBetw_2_rand' - auc_MNodeBetw_1_rand', p);

AUC_CLs_21{1,1} = 'MClust';
AUC_CLs_21{2,1} = 'MDegree';
AUC_CLs_21{3,1} = 'Trans';
AUC_CLs_21{4,1} = 'Assort';
AUC_CLs_21{5,1} = 'Eglob';
AUC_CLs_21{6,1} = 'MEloc';
AUC_CLs_21{7,1} = 'Mod';
AUC_CLs_21{8,1} = 'PathL';
AUC_CLs_21{9,1} = 'Lambda';
AUC_CLs_21{10,1} = 'Gamma';
AUC_CLs_21{11,1} = 'Sigma';
AUC_CLs_21{12,1} = 'MEdgeBtwn';
AUC_CLs_21{13,1} = 'MNodeBtwn';

% AUC diff
AUC_CLs_21{1,2} = auc_MClust_2 - auc_MClust_1;
AUC_CLs_21{2,2} = auc_MDeg_2 - auc_MDeg_1;
AUC_CLs_21{3,2} = auc_Trans_2 - auc_Trans_1;
AUC_CLs_21{4,2} = auc_Assort_2 - auc_Assort_1;
AUC_CLs_21{5,2} = auc_GEff_2 - auc_GEff_1;
AUC_CLs_21{6,2} = auc_MLocEff_2 - auc_MLocEff_1;
AUC_CLs_21{7,2} = auc_Mod_2 - auc_Mod_1;
AUC_CLs_21{8,2} = auc_PathL_2 - auc_PathL_1;
AUC_CLs_21{9,2} = auc_Lambda_2 - auc_Lambda_1;
AUC_CLs_21{10,2} = auc_Gamma_2 - auc_Gamma_1;
AUC_CLs_21{11,2} = auc_Sigma_2 - auc_Sigma_1;
AUC_CLs_21{12,2} = auc_MEdgeBetw_2 - auc_MEdgeBetw_1;
AUC_CLs_21{13,2} = auc_MNodeBetw_2 - auc_MNodeBetw_1;

% Lower limit CL
AUC_CLs_21{1,3} = CL_auc_21_MClust(2,:);
AUC_CLs_21{2,3} = CL_auc_21_MDeg(2,:);
AUC_CLs_21{3,3} = CL_auc_21_Trans(2,:);
AUC_CLs_21{4,3} = CL_auc_21_Assort(2,:);
AUC_CLs_21{5,3} = CL_auc_21_Geff(2,:);
AUC_CLs_21{6,3} = CL_auc_21_Leff(2,:);
AUC_CLs_21{7,3} = CL_auc_21_Mod(2,:);
AUC_CLs_21{8,3} = CL_auc_21_PathL(2,:);
AUC_CLs_21{9,3} = CL_auc_21_Lambda(2,:);
AUC_CLs_21{10,3} = CL_auc_21_Gamma(2,:);
AUC_CLs_21{11,3} = CL_auc_21_Sigma(2,:);
AUC_CLs_21{12,3} = CL_auc_21_MEdgeBetw(2,:);
AUC_CLs_21{13,3} = CL_auc_21_MNodeBetw(2,:);

% Upper limit CL
AUC_CLs_21{1,4} = CL_auc_21_MClust(1,:);
AUC_CLs_21{2,4} = CL_auc_21_MDeg(1,:);
AUC_CLs_21{3,4} = CL_auc_21_Trans(1,:);
AUC_CLs_21{4,4} = CL_auc_21_Assort(1,:);
AUC_CLs_21{5,4} = CL_auc_21_Geff(1,:);
AUC_CLs_21{6,4} = CL_auc_21_Leff(1,:);
AUC_CLs_21{7,4} = CL_auc_21_Mod(1,:);
AUC_CLs_21{8,4} = CL_auc_21_PathL(1,:);
AUC_CLs_21{9,4} = CL_auc_21_Lambda(1,:);
AUC_CLs_21{10,4} = CL_auc_21_Gamma(1,:);
AUC_CLs_21{11,4} = CL_auc_21_Sigma(1,:);
AUC_CLs_21{12,4} = CL_auc_21_MEdgeBetw(1,:);
AUC_CLs_21{13,4} = CL_auc_21_MNodeBetw(1,:);

% p values
AUC_CLs_21{1,5} = CL_Pval(auc_MClust_2_rand - auc_MClust_1_rand, auc_MClust_2 - auc_MClust_1,'AUC_MClust',Tail);
AUC_CLs_21{2,5} = CL_Pval(auc_MDeg_2_rand - auc_MDeg_1_rand, auc_MDeg_2 - auc_MDeg_1,'AUC_MDeg',Tail);
AUC_CLs_21{3,5} = CL_Pval(auc_Trans_2_rand - auc_Trans_1_rand, auc_Trans_2 - auc_Trans_1,'AUC_Trans',Tail);
AUC_CLs_21{4,5} = CL_Pval(auc_Assort_2_rand - auc_Assort_1_rand, auc_Assort_2 - auc_Assort_1,'AUC_Assort',Tail);
AUC_CLs_21{5,5} = CL_Pval(auc_GEff_2_rand - auc_GEff_1_rand, auc_GEff_2 - auc_GEff_1,'AUC_GEff',Tail);
AUC_CLs_21{6,5} = CL_Pval(auc_MLocEff_2_rand - auc_MLocEff_1_rand, auc_MLocEff_2 - auc_MLocEff_1,'AUC_MLocEff',Tail);
AUC_CLs_21{7,5} = CL_Pval(auc_Mod_2_rand - auc_Mod_1_rand, auc_Mod_2 - auc_Mod_1,'AUC_Mod',Tail);
AUC_CLs_21{8,5} = CL_Pval(auc_PathL_2_rand - auc_PathL_1_rand, auc_PathL_2 - auc_PathL_1,'AUC_PathL',Tail);
AUC_CLs_21{9,5} = CL_Pval(auc_Lambda_2_rand - auc_Lambda_1_rand, auc_Lambda_2 - auc_Lambda_1,'AUC_Lambda',Tail);
AUC_CLs_21{10,5} = CL_Pval(auc_Gamma_2_rand - auc_Gamma_1_rand, auc_Gamma_2 - auc_Gamma_1,'AUC_Gamma',Tail);
AUC_CLs_21{11,5} = CL_Pval(auc_Sigma_2_rand - auc_Sigma_1_rand, auc_Sigma_2 - auc_Sigma_1,'AUC_Sigma',Tail);
AUC_CLs_21{12,5} = CL_Pval(auc_MEdgeBetw_2_rand - auc_MEdgeBetw_1_rand, auc_MEdgeBetw_2 - auc_MEdgeBetw_1,'AUC_MEdgeBetw',Tail);
AUC_CLs_21{13,5} = CL_Pval(auc_MNodeBetw_2_rand - auc_MNodeBetw_1_rand, auc_MNodeBetw_2 - auc_MNodeBetw_1,'AUC_MNodeBetw',Tail);


[h crit_p adj_p]=fdr_bh(abs(cell2mat(AUC_CLs_21(:,5))));
AUC_CLs_21(:,6)=num2cell(adj_p);

colnames = {'measure','meandiff','LL','UL','pval','pFDR'};

resultsTable_21 = cell2table(AUC_CLs_21,'VariableNames',colnames);
save AUC_ResultsTable_21 resultsTable_21
writetable(resultsTable_21,'AUCresultsTable_21.xlsx')

%% Group2 vs Group3
load(NetM_2_3);

Dens_2_rand=Dens_2_rand(:,MinIdx:MaxIdx,:);
MClust_2_rand=MClust_2_rand(:,MinIdx:MaxIdx,:);
MDeg_2_rand=MDeg_2_rand(:,MinIdx:MaxIdx,:);
Trans_2_rand=Trans_2_rand(:,MinIdx:MaxIdx,:);
Assort_2_rand=Assort_2_rand(:,MinIdx:MaxIdx,:);
GEff_2_rand=GEff_2_rand(:,MinIdx:MaxIdx,:);
MLocEff_2_rand=MLocEff_2_rand(:,MinIdx:MaxIdx,:);
Mod_2_rand=Mod_2_rand(:,MinIdx:MaxIdx,:);
ModL_2_rand=ModL_2_rand(:,MinIdx:MaxIdx,:);
PathL_2_rand=PathL_2_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_2_rand=MNodeBetw_2_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_2_rand=MEdgeBetw_2_rand(:,MinIdx:MaxIdx,:);
Lambda_2_rand=Lambda_2_rand(:,MinIdx:MaxIdx,:);
Gamma_2_rand=Gamma_2_rand(:,MinIdx:MaxIdx,:);
Sigma_2_rand=Sigma_2_rand(:,MinIdx:MaxIdx,:);

Dens_3_rand=Dens_3_rand(:,MinIdx:MaxIdx,:);
MClust_3_rand=MClust_3_rand(:,MinIdx:MaxIdx,:);
MDeg_3_rand=MDeg_3_rand(:,MinIdx:MaxIdx,:);
Trans_3_rand=Trans_3_rand(:,MinIdx:MaxIdx,:);
Assort_3_rand=Assort_3_rand(:,MinIdx:MaxIdx,:);
GEff_3_rand=GEff_3_rand(:,MinIdx:MaxIdx,:);
MLocEff_3_rand=MLocEff_3_rand(:,MinIdx:MaxIdx,:);
Mod_3_rand=Mod_3_rand(:,MinIdx:MaxIdx,:);
ModL_3_rand=ModL_3_rand(:,MinIdx:MaxIdx,:);
PathL_3_rand=PathL_3_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_3_rand=MNodeBetw_3_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_3_rand=MEdgeBetw_3_rand(:,MinIdx:MaxIdx,:);
Lambda_3_rand=Lambda_3_rand(:,MinIdx:MaxIdx,:);
Gamma_3_rand=Gamma_3_rand(:,MinIdx:MaxIdx,:);
Sigma_3_rand=Sigma_3_rand(:,MinIdx:MaxIdx,:);

Dens_2=Dens_2(:,MinIdx:MaxIdx);
MClust_2=MClust_2(:,MinIdx:MaxIdx);
MDeg_2=MDeg_2(:,MinIdx:MaxIdx);
Trans_2=Trans_2(:,MinIdx:MaxIdx);
Assort_2=Assort_2(:,MinIdx:MaxIdx);
GEff_2=GEff_2(:,MinIdx:MaxIdx);
MLocEff_2=MLocEff_2(:,MinIdx:MaxIdx);
Mod_2=Mod_2(:,MinIdx:MaxIdx);
ModL_2=ModL_2(:,MinIdx:MaxIdx);
PathL_2=PathL_2(:,MinIdx:MaxIdx);
MNodeBetw_2=MNodeBetw_2(:,MinIdx:MaxIdx);
MEdgeBetw_2=MEdgeBetw_2(:,MinIdx:MaxIdx);
Lambda_2=Lambda_2(:,MinIdx:MaxIdx);
Gamma_2=Gamma_2(:,MinIdx:MaxIdx);
Sigma_2=Sigma_2(:,MinIdx:MaxIdx);

Dens_3=Dens_3(:,MinIdx:MaxIdx);
MClust_3=MClust_3(:,MinIdx:MaxIdx);
MDeg_3=MDeg_3(:,MinIdx:MaxIdx);
Trans_3=Trans_3(:,MinIdx:MaxIdx);
Assort_3=Assort_3(:,MinIdx:MaxIdx);
GEff_3=GEff_3(:,MinIdx:MaxIdx);
MLocEff_3=MLocEff_3(:,MinIdx:MaxIdx);
Mod_3=Mod_3(:,MinIdx:MaxIdx);
ModL_3=ModL_3(:,MinIdx:MaxIdx);
PathL_3=PathL_3(:,MinIdx:MaxIdx);
MNodeBetw_3=MNodeBetw_3(:,MinIdx:MaxIdx);
MEdgeBetw_3=MEdgeBetw_3(:,MinIdx:MaxIdx);
Lambda_3=Lambda_3(:,MinIdx:MaxIdx);
Gamma_3=Gamma_3(:,MinIdx:MaxIdx);
Sigma_3=Sigma_3(:,MinIdx:MaxIdx);

fprintf('%-4s\n',['saving net measures across new threshold range...']);
cd ../..
dd=pwd;
mkdir('AUC_Results/AUC_Results_G2_vs_G3');
cd([dd '/AUC_Results/AUC_Results_G2_vs_G3']);

save NetMesPerDens_AveAcrossDens Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
    Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
    PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
    Sigma_2_rand Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ModL_2 PathL_2 MNodeBetw_2 ...
    MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 ...
    Dens_3_rand MClust_3_rand MDeg_3_rand Trans_3_rand ...
    Assort_3_rand GEff_3_rand MLocEff_3_rand Mod_3_rand ModL_3_rand ...
    PathL_3_rand MNodeBetw_3_rand MEdgeBetw_3_rand Lambda_3_rand Gamma_3_rand ...
    Sigma_3_rand Dens_3 MClust_3 MDeg_3 Trans_3 Assort_3 GEff_3 MLocEff_3 Mod_3 ModL_3 PathL_3 MNodeBetw_3 ...
    MEdgeBetw_3 Lambda_3 Gamma_3 Sigma_3

Xax=MinMesPlot:MesStepPlot:MaxMesPlot;
iXax=Xax(MinIdx:MaxIdx);

auc_MClust_2 = trapz(iXax,MClust_2,2);
auc_MDeg_2 = trapz(iXax,MDeg_2,2);
auc_Trans_2 = trapz(iXax,Trans_2,2);
auc_Assort_2 = trapz(iXax,Assort_2,2);
auc_GEff_2 = trapz(iXax,GEff_2,2);
auc_MLocEff_2 = trapz(iXax,MLocEff_2,2);
auc_Mod_2 = trapz(iXax,Mod_2,2);
auc_ModL_2 = trapz(iXax,ModL_2,2);
auc_PathL_2 = trapz(iXax,PathL_2,2);
auc_MNodeBetw_2 = trapz(iXax,MNodeBetw_2,2);
auc_MEdgeBetw_2 = trapz(iXax,MEdgeBetw_2,2);
auc_Lambda_2 = trapz(iXax,Lambda_2,2);
auc_Gamma_2 = trapz(iXax,Gamma_2,2);
auc_Sigma_2 = trapz(iXax,Sigma_2,2);

auc_MClust_3 = trapz(iXax,MClust_3,2);
auc_MDeg_3 = trapz(iXax,MDeg_3,2);
auc_Trans_3 = trapz(iXax,Trans_3,2);
auc_Assort_3 = trapz(iXax,Assort_3,2);
auc_GEff_3 = trapz(iXax,GEff_3,2);
auc_MLocEff_3 = trapz(iXax,MLocEff_3,2);
auc_Mod_3 = trapz(iXax,Mod_3,2);
auc_ModL_3 = trapz(iXax,ModL_3,2);
auc_PathL_3 = trapz(iXax,PathL_3,2);
auc_MNodeBetw_3 = trapz(iXax,MNodeBetw_3,2);
auc_MEdgeBetw_3 = trapz(iXax,MEdgeBetw_3,2);
auc_Lambda_3 = trapz(iXax,Lambda_3,2);
auc_Gamma_3 = trapz(iXax,Gamma_3,2);
auc_Sigma_3 = trapz(iXax,Sigma_3,2);

auc_MClust_2_rand = trapz(iXax,MClust_2_rand,2);
auc_MDeg_2_rand = trapz(iXax,MDeg_2_rand,2);
auc_Trans_2_rand = trapz(iXax,Trans_2_rand,2);
auc_Assort_2_rand = trapz(iXax,Assort_2_rand,2);
auc_GEff_2_rand = trapz(iXax,GEff_2_rand,2);
auc_MLocEff_2_rand = trapz(iXax,MLocEff_2_rand,2);
auc_Mod_2_rand = trapz(iXax,Mod_2_rand,2);
auc_ModL_2_rand = trapz(iXax,ModL_2_rand,2);
auc_PathL_2_rand = trapz(iXax,PathL_2_rand,2);
auc_MNodeBetw_2_rand = trapz(iXax,MNodeBetw_2_rand,2);
auc_MEdgeBetw_2_rand = trapz(iXax,MEdgeBetw_2_rand,2);
auc_Lambda_2_rand = trapz(iXax,Lambda_2_rand,2);
auc_Gamma_2_rand = trapz(iXax,Gamma_2_rand,2);
auc_Sigma_2_rand = trapz(iXax,Sigma_2_rand,2);

auc_MClust_3_rand = trapz(iXax,MClust_3_rand,2);
auc_MDeg_3_rand = trapz(iXax,MDeg_3_rand,2);
auc_Trans_3_rand = trapz(iXax,Trans_3_rand,2);
auc_Assort_3_rand = trapz(iXax,Assort_3_rand,2);
auc_GEff_3_rand = trapz(iXax,GEff_3_rand,2);
auc_MLocEff_3_rand = trapz(iXax,MLocEff_3_rand,2);
auc_Mod_3_rand = trapz(iXax,Mod_3_rand,2);
auc_ModL_3_rand = trapz(iXax,ModL_3_rand,2);
auc_PathL_3_rand = trapz(iXax,PathL_3_rand,2);
auc_MNodeBetw_3_rand = trapz(iXax,MNodeBetw_3_rand,2);
auc_MEdgeBetw_3_rand = trapz(iXax,MEdgeBetw_3_rand,2);
auc_Lambda_3_rand = trapz(iXax,Lambda_3_rand,2);
auc_Gamma_3_rand = trapz(iXax,Gamma_3_rand,2);
auc_Sigma_3_rand = trapz(iXax,Sigma_3_rand,2);

save AUC_NetMesPerDens_AveAcrossDens auc_MClust_2_rand  auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_MClust_2 auc_MDeg_2 auc_Trans_2 auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 auc_MClust_3_rand auc_MDeg_3_rand auc_Trans_3_rand ...
    auc_Assort_3_rand auc_GEff_3_rand auc_MLocEff_3_rand auc_Mod_3_rand auc_ModL_3_rand ...
    auc_PathL_3_rand auc_MNodeBetw_3_rand auc_MEdgeBetw_3_rand auc_Lambda_3_rand auc_Gamma_3_rand ...
    auc_Sigma_3_rand auc_MClust_3 auc_MDeg_3 auc_Trans_3 auc_Assort_3 auc_GEff_3 auc_MLocEff_3 auc_Mod_3 auc_ModL_3 auc_PathL_3 ...
    auc_MNodeBetw_3 auc_MEdgeBetw_3 auc_Lambda_3 auc_Gamma_3 auc_Sigma_3


AUC_indiv_mat= num2cell([[1:ff.mat4GAT.n2 1:ff.mat4GAT.n3]' [ones(ff.mat4GAT.n2,1)*2; ones(ff.mat4GAT.n3,1)*3]  [auc_MClust_2; auc_MClust_3] ...
    [auc_MDeg_2; auc_MDeg_3] [auc_Trans_2; auc_Trans_3] [auc_Assort_2; auc_Assort_3] [auc_GEff_2; auc_GEff_3] [auc_MLocEff_2; auc_MLocEff_3] ...
    [auc_Mod_2; auc_Mod_3] [auc_PathL_2; auc_PathL_3] [auc_Lambda_2; auc_Lambda_3] [auc_Gamma_2; auc_Gamma_3] [auc_Sigma_2; auc_Sigma_3]...
    [auc_MEdgeBetw_2; auc_MEdgeBetw_3] [auc_MNodeBetw_2; auc_MNodeBetw_3]]);
colnames = {'Subj','Group','MClust','MDegree','Trans','Assort','Eglob','MEloc','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBtwn','MNodeBtwn'};
AUCdata = cell2table(AUC_indiv_mat,'VariableNames',colnames);
save AUCdata AUCdata
writetable(AUCdata,'AUCdata.xlsx')


nRand=size(auc_MClust_2_rand,3);

auc_MClust_2 = mean(auc_MClust_2);
auc_MDeg_2 = mean(auc_MDeg_2);
auc_Trans_2 = mean(auc_Trans_2);
auc_Assort_2 = mean(auc_Assort_2);
auc_GEff_2 = mean(auc_GEff_2);
auc_MLocEff_2 = mean(auc_MLocEff_2);
auc_Mod_2 = mean(auc_Mod_2);
auc_ModL_2 = mean(auc_ModL_2);
auc_PathL_2 = mean(auc_PathL_2);
auc_MNodeBetw_2 = mean(auc_MNodeBetw_2);
auc_MEdgeBetw_2 = mean(auc_MEdgeBetw_2);
auc_Lambda_2 = mean(auc_Lambda_2);
auc_Gamma_2 = mean(auc_Gamma_2);
auc_Sigma_2 = mean(auc_Sigma_2);

auc_MClust_3 = mean(auc_MClust_3);
auc_MDeg_3 = mean(auc_MDeg_3);
auc_Trans_3 = mean(auc_Trans_3);
auc_Assort_3 = mean(auc_Assort_3);
auc_GEff_3 = mean(auc_GEff_3);
auc_MLocEff_3 = mean(auc_MLocEff_3);
auc_Mod_3 = mean(auc_Mod_3);
auc_ModL_3 = mean(auc_ModL_3);
auc_PathL_3 = mean(auc_PathL_3);
auc_MNodeBetw_3 = mean(auc_MNodeBetw_3);
auc_MEdgeBetw_3 = mean(auc_MEdgeBetw_3);
auc_Lambda_3 = mean(auc_Lambda_3);
auc_Gamma_3 = mean(auc_Gamma_3);
auc_Sigma_3 = mean(auc_Sigma_3);

auc_MClust_2_rand = mean(auc_MClust_2_rand,1);auc_MClust_2_rand = reshape(auc_MClust_2_rand,1,nRand);
auc_MDeg_2_rand = mean(auc_MDeg_2_rand,1);auc_MDeg_2_rand = reshape(auc_MDeg_2_rand,1,nRand);
auc_Trans_2_rand = mean(auc_Trans_2_rand,1);auc_Trans_2_rand = reshape(auc_Trans_2_rand ,1,nRand);
auc_Assort_2_rand = mean(auc_Assort_2_rand,1);auc_Assort_2_rand = reshape(auc_Assort_2_rand ,1,nRand);
auc_GEff_2_rand = mean(auc_GEff_2_rand,1);auc_GEff_2_rand = reshape(auc_GEff_2_rand ,1,nRand);
auc_MLocEff_2_rand = mean(auc_MLocEff_2_rand,1);auc_MLocEff_2_rand = reshape(auc_MLocEff_2_rand ,1,nRand);
auc_Mod_2_rand = mean(auc_Mod_2_rand,1);auc_Mod_2_rand = reshape(auc_Mod_2_rand ,1,nRand);
auc_ModL_2_rand = mean(auc_ModL_2_rand,1);auc_ModL_2_rand = reshape(auc_ModL_2_rand ,1,nRand);
auc_PathL_2_rand = mean(auc_PathL_2_rand,1);auc_PathL_2_rand = reshape(auc_PathL_2_rand ,1,nRand);
auc_MNodeBetw_2_rand = mean(auc_MNodeBetw_2_rand,1);auc_MNodeBetw_2_rand = reshape(auc_MNodeBetw_2_rand ,1,nRand);
auc_MEdgeBetw_2_rand = mean(auc_MEdgeBetw_2_rand,1);auc_MEdgeBetw_2_rand = reshape(auc_MEdgeBetw_2_rand ,1,nRand);
auc_Lambda_2_rand = mean(auc_Lambda_2_rand,1);auc_Lambda_2_rand = reshape( auc_Lambda_2_rand,1,nRand);
auc_Gamma_2_rand = mean(auc_Gamma_2_rand,1);auc_Gamma_2_rand = reshape(auc_Gamma_2_rand ,1,nRand);
auc_Sigma_2_rand = mean(auc_Sigma_2_rand,1);auc_Sigma_2_rand = reshape(auc_Sigma_2_rand ,1,nRand);

auc_MClust_3_rand = mean(auc_MClust_3_rand,1);auc_MClust_3_rand = reshape(auc_MClust_3_rand,1,nRand);
auc_MDeg_3_rand = mean(auc_MDeg_3_rand,1);auc_MDeg_3_rand = reshape(auc_MDeg_3_rand,1,nRand);
auc_Trans_3_rand = mean(auc_Trans_3_rand,1);auc_Trans_3_rand = reshape(auc_Trans_3_rand ,1,nRand);
auc_Assort_3_rand = mean(auc_Assort_3_rand,1);auc_Assort_3_rand = reshape(auc_Assort_3_rand ,1,nRand);
auc_GEff_3_rand = mean(auc_GEff_3_rand,1);auc_GEff_3_rand = reshape(auc_GEff_3_rand ,1,nRand);
auc_MLocEff_3_rand = mean(auc_MLocEff_3_rand,1);auc_MLocEff_3_rand = reshape(auc_MLocEff_3_rand ,1,nRand);
auc_Mod_3_rand = mean(auc_Mod_3_rand,1);auc_Mod_3_rand = reshape(auc_Mod_3_rand ,1,nRand);
auc_ModL_3_rand = mean(auc_ModL_3_rand,1);auc_ModL_3_rand = reshape(auc_ModL_3_rand ,1,nRand);
auc_PathL_3_rand = mean(auc_PathL_3_rand,1);auc_PathL_3_rand = reshape(auc_PathL_3_rand ,1,nRand);
auc_MNodeBetw_3_rand = mean(auc_MNodeBetw_3_rand,1);auc_MNodeBetw_3_rand = reshape(auc_MNodeBetw_3_rand ,1,nRand);
auc_MEdgeBetw_3_rand = mean(auc_MEdgeBetw_3_rand,1);auc_MEdgeBetw_3_rand = reshape(auc_MEdgeBetw_3_rand ,1,nRand);
auc_Lambda_3_rand = mean(auc_Lambda_3_rand,1);auc_Lambda_3_rand = reshape( auc_Lambda_3_rand,1,nRand);
auc_Gamma_3_rand = mean(auc_Gamma_3_rand,1);auc_Gamma_3_rand = reshape(auc_Gamma_3_rand ,1,nRand);
auc_Sigma_3_rand = mean(auc_Sigma_3_rand,1);auc_Sigma_3_rand = reshape(auc_Sigma_3_rand ,1,nRand);

save AUC_Mean auc_MClust_2_rand auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_MClust_2 auc_MDeg_2 auc_Trans_2 ...
    auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 auc_MClust_3_rand auc_MDeg_3_rand auc_Trans_3_rand ...
    auc_Assort_3_rand auc_GEff_3_rand auc_MLocEff_3_rand auc_Mod_3_rand auc_ModL_3_rand ...
    auc_PathL_3_rand auc_MNodeBetw_3_rand auc_MEdgeBetw_3_rand auc_Lambda_3_rand auc_Gamma_3_rand ...
    auc_Sigma_3_rand auc_MClust_3 auc_MDeg_3 auc_Trans_3 ...
    auc_Assort_3 auc_GEff_3 auc_MLocEff_3 auc_Mod_3 auc_ModL_3 auc_PathL_3 auc_MNodeBetw_3 ...
    auc_MEdgeBetw_3 auc_Lambda_3 auc_Gamma_3 auc_Sigma_3

CL_auc_32_Geff = CL_per(auc_GEff_3_rand' - auc_GEff_2_rand', p);
CL_auc_32_Leff = CL_per(auc_MLocEff_3_rand' - auc_MLocEff_2_rand', p);
CL_auc_32_Trans = CL_per(auc_Trans_3_rand' - auc_Trans_2_rand', p);
CL_auc_32_Mod = CL_per(auc_Mod_3_rand' - auc_Mod_2_rand', p);
CL_auc_32_Lambda = CL_per(auc_Lambda_3_rand' - auc_Lambda_2_rand', p);
CL_auc_32_Gamma = CL_per(auc_Gamma_3_rand' - auc_Gamma_2_rand', p);
CL_auc_32_Sigma = CL_per(auc_Sigma_3_rand' - auc_Sigma_2_rand', p);
CL_auc_32_Assort = CL_per(auc_Assort_3_rand' - auc_Assort_2_rand', p);
CL_auc_32_PathL = CL_per(auc_PathL_3_rand' - auc_PathL_2_rand', p);
CL_auc_32_MClust = CL_per(auc_MClust_3_rand' - auc_MClust_2_rand', p);
CL_auc_32_MDeg = CL_per(auc_MDeg_3_rand' - auc_MDeg_2_rand', p);
CL_auc_32_MEdgeBetw = CL_per(auc_MEdgeBetw_3_rand' - auc_MEdgeBetw_2_rand', p);
CL_auc_32_MNodeBetw = CL_per(auc_MNodeBetw_3_rand' - auc_MNodeBetw_2_rand', p);

AUC_CLs_32{1,1} = 'MClust';
AUC_CLs_32{2,1} = 'MDegree';
AUC_CLs_32{3,1} = 'Trans';
AUC_CLs_32{4,1} = 'Assort';
AUC_CLs_32{5,1} = 'Eglob';
AUC_CLs_32{6,1} = 'MEloc';
AUC_CLs_32{7,1} = 'Mod';
AUC_CLs_32{8,1} = 'PathL';
AUC_CLs_32{9,1} = 'Lambda';
AUC_CLs_32{10,1} = 'Gamma';
AUC_CLs_32{11,1} = 'Sigma';
AUC_CLs_32{12,1} = 'MEdgeBtwn';
AUC_CLs_32{13,1} = 'MNodeBtwn';

% AUC diff
AUC_CLs_32{1,2} = auc_MClust_3 - auc_MClust_2;
AUC_CLs_32{2,2} = auc_MDeg_3 - auc_MDeg_2;
AUC_CLs_32{3,2} = auc_Trans_3 - auc_Trans_2;
AUC_CLs_32{4,2} = auc_Assort_3 - auc_Assort_2;
AUC_CLs_32{5,2} = auc_GEff_3 - auc_GEff_2;
AUC_CLs_32{6,2} = auc_MLocEff_3 - auc_MLocEff_2;
AUC_CLs_32{7,2} = auc_Mod_3 - auc_Mod_2;
AUC_CLs_32{8,2} = auc_PathL_3 - auc_PathL_2;
AUC_CLs_32{9,2} = auc_Lambda_3 - auc_Lambda_2;
AUC_CLs_32{10,2} = auc_Gamma_3 - auc_Gamma_2;
AUC_CLs_32{11,2} = auc_Sigma_3 - auc_Sigma_2;
AUC_CLs_32{12,2} = auc_MEdgeBetw_3 - auc_MEdgeBetw_2;
AUC_CLs_32{13,2} = auc_MNodeBetw_3 - auc_MNodeBetw_2;

% Lower limit CL
AUC_CLs_32{1,3} = CL_auc_32_MClust(2,:);
AUC_CLs_32{2,3} = CL_auc_32_MDeg(2,:);
AUC_CLs_32{3,3} = CL_auc_32_Trans(2,:);
AUC_CLs_32{4,3} = CL_auc_32_Assort(2,:);
AUC_CLs_32{5,3} = CL_auc_32_Geff(2,:);
AUC_CLs_32{6,3} = CL_auc_32_Leff(2,:);
AUC_CLs_32{7,3} = CL_auc_32_Mod(2,:);
AUC_CLs_32{8,3} = CL_auc_32_PathL(2,:);
AUC_CLs_32{9,3} = CL_auc_32_Lambda(2,:);
AUC_CLs_32{10,3} = CL_auc_32_Gamma(2,:);
AUC_CLs_32{11,3} = CL_auc_32_Sigma(2,:);
AUC_CLs_32{12,3} = CL_auc_32_MEdgeBetw(2,:);
AUC_CLs_32{13,3} = CL_auc_32_MNodeBetw(2,:);


% Upper limit CL
AUC_CLs_32{1,4} = CL_auc_32_MClust(1,:);
AUC_CLs_32{2,4} = CL_auc_32_MDeg(1,:);
AUC_CLs_32{3,4} = CL_auc_32_Trans(1,:);
AUC_CLs_32{4,4} = CL_auc_32_Assort(1,:);
AUC_CLs_32{5,4} = CL_auc_32_Geff(1,:);
AUC_CLs_32{6,4} = CL_auc_32_Leff(1,:);
AUC_CLs_32{7,4} = CL_auc_32_Mod(1,:);
AUC_CLs_32{8,4} = CL_auc_32_PathL(1,:);
AUC_CLs_32{9,4} = CL_auc_32_Lambda(1,:);
AUC_CLs_32{10,4} = CL_auc_32_Gamma(1,:);
AUC_CLs_32{11,4} = CL_auc_32_Sigma(1,:);
AUC_CLs_32{12,4} = CL_auc_32_MEdgeBetw(1,:);
AUC_CLs_32{13,4} = CL_auc_32_MNodeBetw(1,:);

% p values
AUC_CLs_32{1,5} = CL_Pval(auc_MClust_3_rand - auc_MClust_2_rand, auc_MClust_3 - auc_MClust_2,'AUC_MClust',Tail);
AUC_CLs_32{2,5} = CL_Pval(auc_MDeg_3_rand - auc_MDeg_2_rand, auc_MDeg_3 - auc_MDeg_2,'AUC_MClust',Tail);
AUC_CLs_32{3,5} = CL_Pval(auc_Trans_3_rand - auc_Trans_2_rand, auc_Trans_3 - auc_Trans_2,'AUC_Trans',Tail);
AUC_CLs_32{4,5} = CL_Pval(auc_Assort_3_rand - auc_Assort_2_rand, auc_Assort_3 - auc_Assort_2,'AUC_Assort',Tail);
AUC_CLs_32{5,5} = CL_Pval(auc_GEff_3_rand - auc_GEff_2_rand, auc_GEff_3 - auc_GEff_2,'AUC_GEff',Tail);
AUC_CLs_32{6,5} = CL_Pval(auc_MLocEff_3_rand - auc_MLocEff_2_rand, auc_MLocEff_3 - auc_MLocEff_2,'AUC_MLocEff',Tail);
AUC_CLs_32{7,5} = CL_Pval(auc_Mod_3_rand - auc_Mod_2_rand, auc_Mod_3 - auc_Mod_2,'AUC_Mod',Tail);
AUC_CLs_32{8,5} = CL_Pval(auc_PathL_3_rand - auc_PathL_2_rand, auc_PathL_3 - auc_PathL_2,'AUC_PathL',Tail);
AUC_CLs_32{9,5} = CL_Pval(auc_Lambda_3_rand - auc_Lambda_2_rand, auc_Lambda_3 - auc_Lambda_2,'AUC_Lambda',Tail);
AUC_CLs_32{10,5} = CL_Pval(auc_Gamma_3_rand - auc_Gamma_2_rand, auc_Gamma_3 - auc_Gamma_2,'AUC_Gamma',Tail);
AUC_CLs_32{11,5} = CL_Pval(auc_Sigma_3_rand - auc_Sigma_2_rand, auc_Sigma_3 - auc_Sigma_2,'AUC_Sigma',Tail);
AUC_CLs_32{12,5} = CL_Pval(auc_MEdgeBetw_3_rand - auc_MEdgeBetw_2_rand, auc_MEdgeBetw_3 - auc_MEdgeBetw_2,'AUC_MEdgeBetw',Tail);
AUC_CLs_32{13,5} = CL_Pval(auc_MNodeBetw_3_rand - auc_MNodeBetw_2_rand, auc_MNodeBetw_3 - auc_MNodeBetw_2,'AUC_MNodeBetw',Tail);

[h crit_p adj_p]=fdr_bh(abs(cell2mat(AUC_CLs_32(:,5))));
AUC_CLs_32(:,6)=num2cell(adj_p);

colnames = {'measure','meandiff','LL','UL','pval','pFDR'};
resultsTable_32 = cell2table(AUC_CLs_32,'VariableNames',colnames);
save AUC_ResultsTable_32 resultsTable_32
writetable(resultsTable_32,'AUCresultsTable_32.xlsx')

%% Group1 vs Group3
load(NetM_1_3);

Dens_1_rand=Dens_1_rand(:,MinIdx:MaxIdx,:);
MClust_1_rand=MClust_1_rand(:,MinIdx:MaxIdx,:);
MDeg_1_rand=MDeg_1_rand(:,MinIdx:MaxIdx,:);
Trans_1_rand=Trans_1_rand(:,MinIdx:MaxIdx,:);
Assort_1_rand=Assort_1_rand(:,MinIdx:MaxIdx,:);
GEff_1_rand=GEff_1_rand(:,MinIdx:MaxIdx,:);
MLocEff_1_rand=MLocEff_1_rand(:,MinIdx:MaxIdx,:);
Mod_1_rand=Mod_1_rand(:,MinIdx:MaxIdx,:);
ModL_1_rand=ModL_1_rand(:,MinIdx:MaxIdx,:);
PathL_1_rand=PathL_1_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_1_rand=MNodeBetw_1_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_1_rand=MEdgeBetw_1_rand(:,MinIdx:MaxIdx,:);
Lambda_1_rand=Lambda_1_rand(:,MinIdx:MaxIdx,:);
Gamma_1_rand=Gamma_1_rand(:,MinIdx:MaxIdx,:);
Sigma_1_rand=Sigma_1_rand(:,MinIdx:MaxIdx,:);

Dens_3_rand=Dens_3_rand(:,MinIdx:MaxIdx,:);
MClust_3_rand=MClust_3_rand(:,MinIdx:MaxIdx,:);
MDeg_3_rand=MDeg_3_rand(:,MinIdx:MaxIdx,:);
Trans_3_rand=Trans_3_rand(:,MinIdx:MaxIdx,:);
Assort_3_rand=Assort_3_rand(:,MinIdx:MaxIdx,:);
GEff_3_rand=GEff_3_rand(:,MinIdx:MaxIdx,:);
MLocEff_3_rand=MLocEff_3_rand(:,MinIdx:MaxIdx,:);
Mod_3_rand=Mod_3_rand(:,MinIdx:MaxIdx,:);
ModL_3_rand=ModL_3_rand(:,MinIdx:MaxIdx,:);
PathL_3_rand=PathL_3_rand(:,MinIdx:MaxIdx,:);
MNodeBetw_3_rand=MNodeBetw_3_rand(:,MinIdx:MaxIdx,:);
MEdgeBetw_3_rand=MEdgeBetw_3_rand(:,MinIdx:MaxIdx,:);
Lambda_3_rand=Lambda_3_rand(:,MinIdx:MaxIdx,:);
Gamma_3_rand=Gamma_3_rand(:,MinIdx:MaxIdx,:);
Sigma_3_rand=Sigma_3_rand(:,MinIdx:MaxIdx,:);

Dens_1=Dens_1(:,MinIdx:MaxIdx);
MClust_1=MClust_1(:,MinIdx:MaxIdx);
MDeg_1=MDeg_1(:,MinIdx:MaxIdx);
Trans_1=Trans_1(:,MinIdx:MaxIdx);
Assort_1=Assort_1(:,MinIdx:MaxIdx);
GEff_1=GEff_1(:,MinIdx:MaxIdx);
MLocEff_1=MLocEff_1(:,MinIdx:MaxIdx);
Mod_1=Mod_1(:,MinIdx:MaxIdx);
ModL_1=ModL_1(:,MinIdx:MaxIdx);
PathL_1=PathL_1(:,MinIdx:MaxIdx);
MNodeBetw_1=MNodeBetw_1(:,MinIdx:MaxIdx);
MEdgeBetw_1=MEdgeBetw_1(:,MinIdx:MaxIdx);
Lambda_1=Lambda_1(:,MinIdx:MaxIdx);
Gamma_1=Gamma_1(:,MinIdx:MaxIdx);
Sigma_1=Sigma_1(:,MinIdx:MaxIdx);

Dens_3=Dens_3(:,MinIdx:MaxIdx);
MClust_3=MClust_3(:,MinIdx:MaxIdx);
MDeg_3=MDeg_3(:,MinIdx:MaxIdx);
Trans_3=Trans_3(:,MinIdx:MaxIdx);
Assort_3=Assort_3(:,MinIdx:MaxIdx);
GEff_3=GEff_3(:,MinIdx:MaxIdx);
MLocEff_3=MLocEff_3(:,MinIdx:MaxIdx);
Mod_3=Mod_3(:,MinIdx:MaxIdx);
ModL_3=ModL_3(:,MinIdx:MaxIdx);
PathL_3=PathL_3(:,MinIdx:MaxIdx);
MNodeBetw_3=MNodeBetw_3(:,MinIdx:MaxIdx);
MEdgeBetw_3=MEdgeBetw_3(:,MinIdx:MaxIdx);
Lambda_3=Lambda_3(:,MinIdx:MaxIdx);
Gamma_3=Gamma_3(:,MinIdx:MaxIdx);
Sigma_3=Sigma_3(:,MinIdx:MaxIdx);

fprintf('%-4s\n',['saving net measures across new threshold range...']);

cd ../..
dd=pwd;
mkdir('AUC_Results/AUC_Results_G1_vs_G3');
cd([dd '/AUC_Results/AUC_Results_G1_vs_G3']);

save NetMesPerDens_AveAcrossDens Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
    Dens_3_rand MClust_3_rand MDeg_3_rand Trans_3_rand ...
    Assort_3_rand GEff_3_rand MLocEff_3_rand Mod_3_rand ModL_3_rand ...
    PathL_3_rand MNodeBetw_3_rand MEdgeBetw_3_rand Lambda_3_rand Gamma_3_rand ...
    Sigma_3_rand Dens_3 MClust_3 MDeg_3 Trans_3 Assort_3 GEff_3 MLocEff_3 Mod_3 ModL_3 PathL_3 MNodeBetw_3 ...
    MEdgeBetw_3 Lambda_3 Gamma_3 Sigma_3

Xax=MinMesPlot:MesStepPlot:MaxMesPlot;
iXax=Xax(MinIdx:MaxIdx);

auc_MClust_1 = trapz(iXax,MClust_1,2);
auc_MDeg_1 = trapz(iXax,MDeg_1,2);
auc_Trans_1 = trapz(iXax,Trans_1,2);
auc_Assort_1 = trapz(iXax,Assort_1,2);
auc_GEff_1 = trapz(iXax,GEff_1,2);
auc_MLocEff_1 = trapz(iXax,MLocEff_1,2);
auc_Mod_1 = trapz(iXax,Mod_1,2);
auc_ModL_1 = trapz(iXax,ModL_1,2);
auc_PathL_1 = trapz(iXax,PathL_1,2);
auc_MNodeBetw_1 = trapz(iXax,MNodeBetw_1,2);
auc_MEdgeBetw_1 = trapz(iXax,MEdgeBetw_1,2);
auc_Lambda_1 = trapz(iXax,Lambda_1,2);
auc_Gamma_1 = trapz(iXax,Gamma_1,2);
auc_Sigma_1 = trapz(iXax,Sigma_1,2);

auc_MClust_3 = trapz(iXax,MClust_3,2);
auc_MDeg_3 = trapz(iXax,MDeg_3,2);
auc_Trans_3 = trapz(iXax,Trans_3,2);
auc_Assort_3 = trapz(iXax,Assort_3,2);
auc_GEff_3 = trapz(iXax,GEff_3,2);
auc_MLocEff_3 = trapz(iXax,MLocEff_3,2);
auc_Mod_3 = trapz(iXax,Mod_3,2);
auc_ModL_3 = trapz(iXax,ModL_3,2);
auc_PathL_3 = trapz(iXax,PathL_3,2);
auc_MNodeBetw_3 = trapz(iXax,MNodeBetw_3,2);
auc_MEdgeBetw_3 = trapz(iXax,MEdgeBetw_3,2);
auc_Lambda_3 = trapz(iXax,Lambda_3,2);
auc_Gamma_3 = trapz(iXax,Gamma_3,2);
auc_Sigma_3 = trapz(iXax,Sigma_3,2);

auc_MClust_1_rand = trapz(iXax,MClust_1_rand,2);
auc_MDeg_1_rand = trapz(iXax,MDeg_1_rand,2);
auc_Trans_1_rand = trapz(iXax,Trans_1_rand,2);
auc_Assort_1_rand = trapz(iXax,Assort_1_rand,2);
auc_GEff_1_rand = trapz(iXax,GEff_1_rand,2);
auc_MLocEff_1_rand = trapz(iXax,MLocEff_1_rand,2);
auc_Mod_1_rand = trapz(iXax,Mod_1_rand,2);
auc_ModL_1_rand = trapz(iXax,ModL_1_rand,2);
auc_PathL_1_rand = trapz(iXax,PathL_1_rand,2);
auc_MNodeBetw_1_rand = trapz(iXax,MNodeBetw_1_rand,2);
auc_MEdgeBetw_1_rand = trapz(iXax,MEdgeBetw_1_rand,2);
auc_Lambda_1_rand = trapz(iXax,Lambda_1_rand,2);
auc_Gamma_1_rand = trapz(iXax,Gamma_1_rand,2);
auc_Sigma_1_rand = trapz(iXax,Sigma_1_rand,2);

auc_MClust_3_rand = trapz(iXax,MClust_3_rand,2);
auc_MDeg_3_rand = trapz(iXax,MDeg_3_rand,2);
auc_Trans_3_rand = trapz(iXax,Trans_3_rand,2);
auc_Assort_3_rand = trapz(iXax,Assort_3_rand,2);
auc_GEff_3_rand = trapz(iXax,GEff_3_rand,2);
auc_MLocEff_3_rand = trapz(iXax,MLocEff_3_rand,2);
auc_Mod_3_rand = trapz(iXax,Mod_3_rand,2);
auc_ModL_3_rand = trapz(iXax,ModL_3_rand,2);
auc_PathL_3_rand = trapz(iXax,PathL_3_rand,2);
auc_MNodeBetw_3_rand = trapz(iXax,MNodeBetw_3_rand,2);
auc_MEdgeBetw_3_rand = trapz(iXax,MEdgeBetw_3_rand,2);
auc_Lambda_3_rand = trapz(iXax,Lambda_3_rand,2);
auc_Gamma_3_rand = trapz(iXax,Gamma_3_rand,2);
auc_Sigma_3_rand = trapz(iXax,Sigma_3_rand,2);

save AUC_NetMesPerDens_AveAcrossDens  auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
    auc_MClust_3_rand auc_MDeg_3_rand auc_Trans_3_rand ...
    auc_Assort_3_rand auc_GEff_3_rand auc_MLocEff_3_rand auc_Mod_3_rand auc_ModL_3_rand ...
    auc_PathL_3_rand auc_MNodeBetw_3_rand auc_MEdgeBetw_3_rand auc_Lambda_3_rand auc_Gamma_3_rand ...
    auc_Sigma_3_rand auc_MClust_3 auc_MDeg_3 auc_Trans_3 auc_Assort_3 auc_GEff_3 auc_MLocEff_3 auc_Mod_3 auc_ModL_3 auc_PathL_3 ...
    auc_MNodeBetw_3 auc_MEdgeBetw_3 auc_Lambda_3 auc_Gamma_3 auc_Sigma_3


AUC_indiv_mat= num2cell([[1:ff.mat4GAT.n1 1:ff.mat4GAT.n3]' [ones(ff.mat4GAT.n1,1)*1; ones(ff.mat4GAT.n3,1)*3]  [auc_MClust_1; auc_MClust_3] ...
    [auc_MDeg_1; auc_MDeg_3] [auc_Trans_1; auc_Trans_3] [auc_Assort_1; auc_Assort_3] [auc_GEff_1; auc_GEff_3] [auc_MLocEff_1; auc_MLocEff_3] ...
    [auc_Mod_1; auc_Mod_3] [auc_PathL_1; auc_PathL_3] [auc_Lambda_1; auc_Lambda_3] [auc_Gamma_1; auc_Gamma_3] [auc_Sigma_1; auc_Sigma_3]...
    [auc_MEdgeBetw_1; auc_MEdgeBetw_3] [auc_MNodeBetw_1; auc_MNodeBetw_3]]);
colnames = {'Subj','Group','MClust','MDegree','Trans','Assort','Eglob','MEloc','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBtwn','MNodeBtwn'};
AUCdata = cell2table(AUC_indiv_mat,'VariableNames',colnames);
save AUCdata AUCdata
writetable(AUCdata,'AUCdata.xlsx')


nRand=size(auc_MClust_1_rand,3);

auc_MClust_1 = mean(auc_MClust_1);
auc_MDeg_1 = mean(auc_MDeg_1);
auc_Trans_1 = mean(auc_Trans_1);
auc_Assort_1 = mean(auc_Assort_1);
auc_GEff_1 = mean(auc_GEff_1);
auc_MLocEff_1 = mean(auc_MLocEff_1);
auc_Mod_1 = mean(auc_Mod_1);
auc_ModL_1 = mean(auc_ModL_1);
auc_PathL_1 = mean(auc_PathL_1);
auc_MNodeBetw_1 = mean(auc_MNodeBetw_1);
auc_MEdgeBetw_1 = mean(auc_MEdgeBetw_1);
auc_Lambda_1 = mean(auc_Lambda_1);
auc_Gamma_1 = mean(auc_Gamma_1);
auc_Sigma_1 = mean(auc_Sigma_1);

auc_MClust_3 = mean(auc_MClust_3);
auc_MDeg_3 = mean(auc_MDeg_3);
auc_Trans_3 = mean(auc_Trans_3);
auc_Assort_3 = mean(auc_Assort_3);
auc_GEff_3 = mean(auc_GEff_3);
auc_MLocEff_3 = mean(auc_MLocEff_3);
auc_Mod_3 = mean(auc_Mod_3);
auc_ModL_3 = mean(auc_ModL_3);
auc_PathL_3 = mean(auc_PathL_3);
auc_MNodeBetw_3 = mean(auc_MNodeBetw_3);
auc_MEdgeBetw_3 = mean(auc_MEdgeBetw_3);
auc_Lambda_3 = mean(auc_Lambda_3);
auc_Gamma_3 = mean(auc_Gamma_3);
auc_Sigma_3 = mean(auc_Sigma_3);

auc_MClust_1_rand = mean(auc_MClust_1_rand,1);auc_MClust_1_rand = reshape(auc_MClust_1_rand,1,nRand);
auc_MDeg_1_rand = mean(auc_MDeg_1_rand,1);auc_MDeg_1_rand = reshape(auc_MDeg_1_rand,1,nRand);
auc_Trans_1_rand = mean(auc_Trans_1_rand,1);auc_Trans_1_rand = reshape(auc_Trans_1_rand ,1,nRand);
auc_Assort_1_rand = mean(auc_Assort_1_rand,1);auc_Assort_1_rand = reshape(auc_Assort_1_rand ,1,nRand);
auc_GEff_1_rand = mean(auc_GEff_1_rand,1);auc_GEff_1_rand = reshape(auc_GEff_1_rand ,1,nRand);
auc_MLocEff_1_rand = mean(auc_MLocEff_1_rand,1);auc_MLocEff_1_rand = reshape(auc_MLocEff_1_rand ,1,nRand);
auc_Mod_1_rand = mean(auc_Mod_1_rand,1);auc_Mod_1_rand = reshape(auc_Mod_1_rand ,1,nRand);
auc_ModL_1_rand = mean(auc_ModL_1_rand,1);auc_ModL_1_rand = reshape(auc_ModL_1_rand ,1,nRand);
auc_PathL_1_rand = mean(auc_PathL_1_rand,1);auc_PathL_1_rand = reshape(auc_PathL_1_rand ,1,nRand);
auc_MNodeBetw_1_rand = mean(auc_MNodeBetw_1_rand,1);auc_MNodeBetw_1_rand = reshape(auc_MNodeBetw_1_rand ,1,nRand);
auc_MEdgeBetw_1_rand = mean(auc_MEdgeBetw_1_rand,1);auc_MEdgeBetw_1_rand = reshape(auc_MEdgeBetw_1_rand ,1,nRand);
auc_Lambda_1_rand = mean(auc_Lambda_1_rand,1);auc_Lambda_1_rand = reshape( auc_Lambda_1_rand,1,nRand);
auc_Gamma_1_rand = mean(auc_Gamma_1_rand,1);auc_Gamma_1_rand = reshape(auc_Gamma_1_rand ,1,nRand);
auc_Sigma_1_rand = mean(auc_Sigma_1_rand,1);auc_Sigma_1_rand = reshape(auc_Sigma_1_rand ,1,nRand);

auc_MClust_3_rand = mean(auc_MClust_3_rand,1);auc_MClust_3_rand = reshape(auc_MClust_3_rand,1,nRand);
auc_MDeg_3_rand = mean(auc_MDeg_3_rand,1);auc_MDeg_3_rand = reshape(auc_MDeg_3_rand,1,nRand);
auc_Trans_3_rand = mean(auc_Trans_3_rand,1);auc_Trans_3_rand = reshape(auc_Trans_3_rand ,1,nRand);
auc_Assort_3_rand = mean(auc_Assort_3_rand,1);auc_Assort_3_rand = reshape(auc_Assort_3_rand ,1,nRand);
auc_GEff_3_rand = mean(auc_GEff_3_rand,1);auc_GEff_3_rand = reshape(auc_GEff_3_rand ,1,nRand);
auc_MLocEff_3_rand = mean(auc_MLocEff_3_rand,1);auc_MLocEff_3_rand = reshape(auc_MLocEff_3_rand ,1,nRand);
auc_Mod_3_rand = mean(auc_Mod_3_rand,1);auc_Mod_3_rand = reshape(auc_Mod_3_rand ,1,nRand);
auc_ModL_3_rand = mean(auc_ModL_3_rand,1);auc_ModL_3_rand = reshape(auc_ModL_3_rand ,1,nRand);
auc_PathL_3_rand = mean(auc_PathL_3_rand,1);auc_PathL_3_rand = reshape(auc_PathL_3_rand ,1,nRand);
auc_MNodeBetw_3_rand = mean(auc_MNodeBetw_3_rand,1);auc_MNodeBetw_3_rand = reshape(auc_MNodeBetw_3_rand ,1,nRand);
auc_MEdgeBetw_3_rand = mean(auc_MEdgeBetw_3_rand,1);auc_MEdgeBetw_3_rand = reshape(auc_MEdgeBetw_3_rand ,1,nRand);
auc_Lambda_3_rand = mean(auc_Lambda_3_rand,1);auc_Lambda_3_rand = reshape( auc_Lambda_3_rand,1,nRand);
auc_Gamma_3_rand = mean(auc_Gamma_3_rand,1);auc_Gamma_3_rand = reshape(auc_Gamma_3_rand ,1,nRand);
auc_Sigma_3_rand = mean(auc_Sigma_3_rand,1);auc_Sigma_3_rand = reshape(auc_Sigma_3_rand ,1,nRand);

save AUC_Mean auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
    auc_MClust_3_rand auc_MDeg_3_rand auc_Trans_3_rand ...
    auc_Assort_3_rand auc_GEff_3_rand auc_MLocEff_3_rand auc_Mod_3_rand auc_ModL_3_rand ...
    auc_PathL_3_rand auc_MNodeBetw_3_rand auc_MEdgeBetw_3_rand auc_Lambda_3_rand auc_Gamma_3_rand ...
    auc_Sigma_3_rand auc_MClust_3 auc_MDeg_3 auc_Trans_3 ...
    auc_Assort_3 auc_GEff_3 auc_MLocEff_3 auc_Mod_3 auc_ModL_3 auc_PathL_3 auc_MNodeBetw_3 ...
    auc_MEdgeBetw_3 auc_Lambda_3 auc_Gamma_3 auc_Sigma_3

p = Alpha;

CL_auc_31_Geff = CL_per(auc_GEff_3_rand' - auc_GEff_1_rand', p);
CL_auc_31_Leff = CL_per(auc_MLocEff_3_rand' - auc_MLocEff_1_rand', p);
CL_auc_31_Trans = CL_per(auc_Trans_3_rand' - auc_Trans_1_rand', p);
CL_auc_31_Mod = CL_per(auc_Mod_3_rand' - auc_Mod_1_rand', p);
CL_auc_31_Lambda = CL_per(auc_Lambda_3_rand' - auc_Lambda_1_rand', p);
CL_auc_31_Gamma = CL_per(auc_Gamma_3_rand' - auc_Gamma_1_rand', p);
CL_auc_31_Sigma = CL_per(auc_Sigma_3_rand' - auc_Sigma_1_rand', p);
CL_auc_31_Assort = CL_per(auc_Assort_3_rand' - auc_Assort_1_rand', p);
CL_auc_31_PathL = CL_per(auc_PathL_3_rand' - auc_PathL_1_rand', p);
CL_auc_31_MClust = CL_per(auc_MClust_3_rand' - auc_MClust_1_rand', p);
CL_auc_31_MDeg = CL_per(auc_MDeg_3_rand' - auc_MDeg_1_rand', p);
CL_auc_31_MEdgeBetw = CL_per(auc_MEdgeBetw_3_rand' - auc_MEdgeBetw_1_rand', p);
CL_auc_31_MNodeBetw = CL_per(auc_MNodeBetw_3_rand' - auc_MNodeBetw_1_rand', p);

AUC_CLs_31{1,1} = 'MClust';
AUC_CLs_31{2,1} = 'MDegree';
AUC_CLs_31{3,1} = 'Trans';
AUC_CLs_31{4,1} = 'Assort';
AUC_CLs_31{5,1} = 'Eglob';
AUC_CLs_31{6,1} = 'MEloc';
AUC_CLs_31{7,1} = 'Mod';
AUC_CLs_31{8,1} = 'PathL';
AUC_CLs_31{9,1} = 'Lambda';
AUC_CLs_31{10,1} = 'Gamma';
AUC_CLs_31{11,1} = 'Sigma';
AUC_CLs_31{12,1} = 'MEdgeBtwn';
AUC_CLs_31{13,1} = 'MNodeBtwn';

% AUC diff
AUC_CLs_31{1,2} = auc_MClust_3 - auc_MClust_1;
AUC_CLs_31{2,2} = auc_MDeg_3 - auc_MDeg_1;
AUC_CLs_31{3,2} = auc_Trans_3 - auc_Trans_1;
AUC_CLs_31{4,2} = auc_Assort_3 - auc_Assort_1;
AUC_CLs_31{5,2} = auc_GEff_3 - auc_GEff_1;
AUC_CLs_31{6,2} = auc_MLocEff_3 - auc_MLocEff_1;
AUC_CLs_31{7,2} = auc_Mod_3 - auc_Mod_1;
AUC_CLs_31{8,2} = auc_PathL_3 - auc_PathL_1;
AUC_CLs_31{9,2} = auc_Lambda_3 - auc_Lambda_1;
AUC_CLs_31{10,2} = auc_Gamma_3 - auc_Gamma_1;
AUC_CLs_31{11,2} = auc_Sigma_3 - auc_Sigma_1;
AUC_CLs_31{12,2} = auc_MEdgeBetw_3 - auc_MEdgeBetw_1;
AUC_CLs_31{13,2} = auc_MNodeBetw_3 - auc_MNodeBetw_1;

% Lower limit CL
AUC_CLs_31{1,3} = CL_auc_31_MClust(2,:);
AUC_CLs_31{2,3} = CL_auc_31_MDeg(2,:);
AUC_CLs_31{3,3} = CL_auc_31_Trans(2,:);
AUC_CLs_31{4,3} = CL_auc_31_Assort(2,:);
AUC_CLs_31{5,3} = CL_auc_31_Geff(2,:);
AUC_CLs_31{6,3} = CL_auc_31_Leff(2,:);
AUC_CLs_31{7,3} = CL_auc_31_Mod(2,:);
AUC_CLs_31{8,3} = CL_auc_31_PathL(2,:);
AUC_CLs_31{9,3} = CL_auc_31_Lambda(2,:);
AUC_CLs_31{10,3} = CL_auc_31_Gamma(2,:);
AUC_CLs_31{11,3} = CL_auc_31_Sigma(2,:);
AUC_CLs_31{12,3} = CL_auc_31_MEdgeBetw(2,:);
AUC_CLs_31{13,3} = CL_auc_31_MNodeBetw(2,:);

% Upper limit CL
AUC_CLs_31{1,4} = CL_auc_31_MClust(1,:);
AUC_CLs_31{2,4} = CL_auc_31_MDeg(1,:);
AUC_CLs_31{3,4} = CL_auc_31_Trans(1,:);
AUC_CLs_31{4,4} = CL_auc_31_Assort(1,:);
AUC_CLs_31{5,4} = CL_auc_31_Geff(1,:);
AUC_CLs_31{6,4} = CL_auc_31_Leff(1,:);
AUC_CLs_31{7,4} = CL_auc_31_Mod(1,:);
AUC_CLs_31{8,4} = CL_auc_31_PathL(1,:);
AUC_CLs_31{9,4} = CL_auc_31_Lambda(1,:);
AUC_CLs_31{10,4} = CL_auc_31_Gamma(1,:);
AUC_CLs_31{11,4} = CL_auc_31_Sigma(1,:);
AUC_CLs_31{12,4} = CL_auc_31_MEdgeBetw(1,:);
AUC_CLs_31{13,4} = CL_auc_31_MNodeBetw(1,:);

% p values
AUC_CLs_31{1,5} = CL_Pval(auc_MClust_3_rand - auc_MClust_1_rand, auc_MClust_3 - auc_MClust_1,'AUC_MClust',Tail);
AUC_CLs_31{2,5} = CL_Pval(auc_MDeg_3_rand - auc_MDeg_1_rand, auc_MDeg_3 - auc_MDeg_1,'AUC_MDeg',Tail);
AUC_CLs_31{3,5} = CL_Pval(auc_Trans_3_rand - auc_Trans_1_rand, auc_Trans_3 - auc_Trans_1,'AUC_Trans',Tail);
AUC_CLs_31{4,5} = CL_Pval(auc_Assort_3_rand - auc_Assort_1_rand, auc_Assort_3 - auc_Assort_1,'AUC_Assort',Tail);
AUC_CLs_31{5,5} = CL_Pval(auc_GEff_3_rand - auc_GEff_1_rand, auc_GEff_3 - auc_GEff_1,'AUC_GEff',Tail);
AUC_CLs_31{6,5} = CL_Pval(auc_MLocEff_3_rand - auc_MLocEff_1_rand, auc_MLocEff_3 - auc_MLocEff_1,'AUC_MLocEff',Tail);
AUC_CLs_31{7,5} = CL_Pval(auc_Mod_3_rand - auc_Mod_1_rand, auc_Mod_3 - auc_Mod_1,'AUC_Mod',Tail);
AUC_CLs_31{8,5} = CL_Pval(auc_PathL_3_rand - auc_PathL_1_rand, auc_PathL_3 - auc_PathL_1,'AUC_PathL',Tail);
AUC_CLs_31{9,5} = CL_Pval(auc_Lambda_3_rand - auc_Lambda_1_rand, auc_Lambda_3 - auc_Lambda_1,'AUC_Lambda',Tail);
AUC_CLs_31{10,5} = CL_Pval(auc_Gamma_3_rand - auc_Gamma_1_rand, auc_Gamma_3 - auc_Gamma_1,'AUC_Gamma',Tail);
AUC_CLs_31{11,5} = CL_Pval(auc_Sigma_3_rand - auc_Sigma_1_rand, auc_Sigma_3 - auc_Sigma_1,'AUC_Sigma',Tail);
AUC_CLs_31{12,5} = CL_Pval(auc_MEdgeBetw_3_rand - auc_MEdgeBetw_1_rand, auc_MEdgeBetw_3 - auc_MEdgeBetw_1,'AUC_MEdgeBetw',Tail);
AUC_CLs_31{13,5} = CL_Pval(auc_MNodeBetw_3_rand - auc_MNodeBetw_1_rand, auc_MNodeBetw_3 - auc_MNodeBetw_1,'AUC_MEdgeBetw',Tail);

[h crit_p adj_p]=fdr_bh(abs(cell2mat(AUC_CLs_31(:,5))));
AUC_CLs_31(:,6)=num2cell(adj_p);

colnames = {'measure','meandiff','LL','UL','pval','pFDR'};
resultsTable_31 = cell2table(AUC_CLs_31,'VariableNames',colnames);
save AUC_ResultsTable_31 resultsTable_31
writetable(resultsTable_31,'AUCresultsTable_31.xlsx')

cd ..
save AUC_CLs AUC_CLs_21  AUC_CLs_32 AUC_CLs_31
cd ..