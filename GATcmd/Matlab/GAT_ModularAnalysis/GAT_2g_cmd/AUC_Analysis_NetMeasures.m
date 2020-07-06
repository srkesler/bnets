function AUC_Analysis_NetMeasures(data_log,NetM,MinThr,MaxThr)

%%
%-- updated on August 21,2011 for functionals
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)
%% Main body

ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;

if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha')
    Alpha = ff.mat4GAT.alpha;
    Tail = ff.mat4GAT.tail;
else
    Alpha = .05;
    Tail = 2;
end

fprintf('%-4s\n',['loading inputs...']);
load(NetM);
MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinIdx=find(single(xxx)==single(MinThr));
MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(Dens_1_rand,2)
    MaxIdx = size(Dens_1_rand,2);
end

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
Mmoddegzscore_1_rand=Mmoddegzscore_1_rand(:,MinIdx:MaxIdx,:);
Mpartcoeff_1_rand=Mpartcoeff_1_rand(:,MinIdx:MaxIdx,:);

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
Mmoddegzscore_2_rand=Mmoddegzscore_2_rand(:,MinIdx:MaxIdx,:);
Mpartcoeff_2_rand=Mpartcoeff_2_rand(:,MinIdx:MaxIdx,:);

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
Mmoddegzscore_1=Mmoddegzscore_1(:,MinIdx:MaxIdx);
Mpartcoeff_1=Mpartcoeff_1(:,MinIdx:MaxIdx);

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
Mmoddegzscore_2=Mmoddegzscore_2(:,MinIdx:MaxIdx);
Mpartcoeff_2=Mpartcoeff_2(:,MinIdx:MaxIdx);

fprintf('%-4s\n',['saving net measures across new threshold range...']);

dd=pwd;
mkdir('AUC_Results');
cd([dd '/AUC_Results']);

save NetMesPerDens_AveAcrossDens Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand Mmoddegzscore_1_rand Mpartcoeff_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
    Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
    PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
    Sigma_2_rand Mmoddegzscore_2_rand Mpartcoeff_2_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 Mmoddegzscore_1 Mpartcoeff_1...
    Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ModL_2 PathL_2 MNodeBetw_2 ...
    MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 Mmoddegzscore_2 Mpartcoeff_2

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
auc_Mmoddegzscore_1 = trapz(iXax,Mmoddegzscore_1,2);
auc_Mpartcoeff_1 = trapz(iXax,Mpartcoeff_1,2);

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
auc_Mmoddegzscore_2 = trapz(iXax,Mmoddegzscore_2,2);
auc_Mpartcoeff_2 = trapz(iXax,Mpartcoeff_2,2);

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
auc_Mmoddegzscore_1_rand = trapz(iXax,Mmoddegzscore_1_rand,2);
auc_Mpartcoeff_1_rand = trapz(iXax,Mpartcoeff_1_rand,2);

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
auc_Mmoddegzscore_2_rand = trapz(iXax,Mmoddegzscore_2_rand,2);
auc_Mpartcoeff_2_rand = trapz(iXax,Mpartcoeff_2_rand,2);

save AUC_NetMesPerDens_AveAcrossDens  auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_Mmoddegzscore_1_rand auc_Mpartcoeff_1_rand auc_MClust_2_rand auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_Mmoddegzscore_2_rand auc_Mpartcoeff_2_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 auc_Mmoddegzscore_1 auc_Mpartcoeff_1...
    auc_MClust_2 auc_MDeg_2 auc_Trans_2 auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 auc_Mmoddegzscore_2 auc_Mpartcoeff_2

AUC_indiv_mat= num2cell([[1:ff.mat4GAT.n1]' [ones(ff.mat4GAT.n1,1)*1]  [auc_MClust_1] ...
    [auc_MDeg_1] [auc_Trans_1] [auc_Assort_1] [auc_GEff_1] [auc_MLocEff_1] ...
    [auc_Mod_1] [auc_PathL_1] [auc_Lambda_1] [auc_Gamma_1] [auc_Sigma_1]...
    [auc_MEdgeBetw_1] [auc_MNodeBetw_1] [auc_Mmoddegzscore_1] [auc_Mpartcoeff_1] ]);

colnames = {'Subj','Group','MClust','MDegree','Trans','Assort','Eglob','MEloc','Mod','PathL','Lambda','Gamma','Sigma','MEdgeBtwn','MNodeBtwn','Mmoddegzscore','Mpartcoeff'};
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
auc_Mmoddegzscore_1 =mean(auc_Mmoddegzscore_1);
auc_Mpartcoeff_1=mean(auc_Mpartcoeff_1);

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
auc_Mmoddegzscore_2 =mean(auc_Mmoddegzscore_2);
auc_Mpartcoeff_2=mean(auc_Mpartcoeff_2);

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
auc_Mmoddegzscore_1_rand = mean(auc_Mmoddegzscore_1_rand,1);auc_Mmoddegzscore_1_rand = reshape(auc_Mmoddegzscore_1_rand ,1,nRand);
auc_Mpartcoeff_1_rand = mean(auc_Mpartcoeff_1_rand,1);auc_Mpartcoeff_1_rand = reshape(auc_Mpartcoeff_1_rand ,1,nRand);

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
auc_Mmoddegzscore_2_rand = mean(auc_Mmoddegzscore_2_rand,1);auc_Mmoddegzscore_2_rand = reshape(auc_Mmoddegzscore_2_rand ,1,nRand);
auc_Mpartcoeff_2_rand = mean(auc_Mpartcoeff_2_rand,1);auc_Mpartcoeff_2_rand = reshape(auc_Mpartcoeff_2_rand ,1,nRand);

save AUC_Mean auc_MClust_1_rand auc_MDeg_1_rand auc_Trans_1_rand ...
    auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
    auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
    auc_Sigma_1_rand auc_Mmoddegzscore_1_rand auc_Mpartcoeff_1_rand auc_MClust_2_rand auc_MDeg_2_rand auc_Trans_2_rand ...
    auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
    auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
    auc_Sigma_2_rand auc_Mmoddegzscore_2_rand auc_Mpartcoeff_2_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
    auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 auc_Mmoddegzscore_1 auc_Mpartcoeff_1...
    auc_MClust_2 auc_MDeg_2 auc_Trans_2 ...
    auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
    auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 auc_Mmoddegzscore_2 auc_Mpartcoeff_2

fprintf('%-4s\n',' ..........Done ');