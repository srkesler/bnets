function AUC_Analysis_NetMeasures(data_log,NetM,MinThr,MaxThr)

%%
%-- updated on August 21,2011 for functionals
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)

%% Main Body of the script common to all group combinations

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

%%  Group1 Vs Group2

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
mkdir('AUC_Results');
cd([dd '/AUC_Results']);

save NetMesPerDens_AveAcrossDens Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
     Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
     PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
     Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
     Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
     PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
     Sigma_2_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
     Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
     Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ModL_2 PathL_2 MNodeBetw_2 ...
     MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 


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
save AUCdata AUCdata ;
writetable(AUCdata,'AUCdata.xlsx');

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
     auc_Sigma_1_rand auc_MClust_2_rand auc_MDeg_2_rand auc_Trans_2_rand ...
     auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
     auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
     auc_Sigma_2_rand auc_MClust_1 auc_MDeg_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
     auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
     auc_MClust_2 auc_MDeg_2 auc_Trans_2 ...
     auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
     auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 

p = Alpha;

CL_auc_Geff = CL_per(auc_GEff_2_rand' - auc_GEff_1_rand', p);
CL_auc_Leff = CL_per(auc_MLocEff_2_rand' - auc_MLocEff_1_rand', p);
CL_auc_Trans = CL_per(auc_Trans_2_rand' - auc_Trans_1_rand', p);
CL_auc_Mod = CL_per(auc_Mod_2_rand' - auc_Mod_1_rand', p);
CL_auc_Lambda = CL_per(auc_Lambda_2_rand' - auc_Lambda_1_rand', p);
CL_auc_Gamma = CL_per(auc_Gamma_2_rand' - auc_Gamma_1_rand', p);
CL_auc_Sigma = CL_per(auc_Sigma_2_rand' - auc_Sigma_1_rand', p);
CL_auc_Assort = CL_per(auc_Assort_2_rand' - auc_Assort_1_rand', p);
CL_auc_PathL = CL_per(auc_PathL_2_rand' - auc_PathL_1_rand', p);
CL_auc_MClust = CL_per(auc_MClust_2_rand' - auc_MClust_1_rand', p);
CL_auc_MDeg = CL_per(auc_MDeg_2_rand' - auc_MDeg_1_rand', p);
CL_auc_MNodeBetw = CL_per(auc_MNodeBetw_2_rand' - auc_MNodeBetw_1_rand', p);
CL_auc_MEdgeBetw = CL_per(auc_MEdgeBetw_2_rand' - auc_MEdgeBetw_1_rand', p);

AUC_CLs{1,1} = 'MClust';
AUC_CLs{2,1} = 'MDegree';
AUC_CLs{3,1} = 'Trans';
AUC_CLs{4,1} = 'Assort';
AUC_CLs{5,1} = 'Eglob';
AUC_CLs{6,1} = 'MEloc';
AUC_CLs{7,1} = 'Mod';
AUC_CLs{8,1} = 'PathL';
AUC_CLs{9,1} = 'Lambda';
AUC_CLs{10,1} = 'Gamma';
AUC_CLs{11,1} = 'Sigma';
AUC_CLs{12,1} = 'MEdgeBtwn';
AUC_CLs{13,1} = 'MNodeBtwn';

% AUC diff
AUC_CLs{1,2} = auc_MClust_2 - auc_MClust_1;
AUC_CLs{2,2} = auc_MDeg_2 - auc_MDeg_1;
AUC_CLs{3,2} = auc_Trans_2 - auc_Trans_1;
AUC_CLs{4,2} = auc_Assort_2 - auc_Assort_1;
AUC_CLs{5,2} = auc_GEff_2 - auc_GEff_1;
AUC_CLs{6,2} = auc_MLocEff_2 - auc_MLocEff_1;
AUC_CLs{7,2} = auc_Mod_2 - auc_Mod_1;
AUC_CLs{8,2} = auc_PathL_2 - auc_PathL_1;
AUC_CLs{9,2} = auc_Lambda_2 - auc_Lambda_1;
AUC_CLs{10,2} = auc_Gamma_2 - auc_Gamma_1;
AUC_CLs{11,2} = auc_Sigma_2 - auc_Sigma_1;
AUC_CLs{12,2} = auc_MEdgeBetw_2 - auc_MEdgeBetw_1;
AUC_CLs{13,2} = auc_MNodeBetw_2 - auc_MNodeBetw_1;

% Lower limit CL
AUC_CLs{1,3} = CL_auc_MClust(2,:);
AUC_CLs{2,3} = CL_auc_MDeg(2,:);
AUC_CLs{3,3} = CL_auc_Trans(2,:);
AUC_CLs{4,3} = CL_auc_Assort(2,:);
AUC_CLs{5,3} = CL_auc_Geff(2,:);
AUC_CLs{6,3} = CL_auc_Leff(2,:);
AUC_CLs{7,3} = CL_auc_Mod(2,:);
AUC_CLs{8,3} = CL_auc_PathL(2,:);
AUC_CLs{9,3} = CL_auc_Lambda(2,:);
AUC_CLs{10,3} = CL_auc_Gamma(2,:);
AUC_CLs{11,3} = CL_auc_Sigma(2,:);
AUC_CLs{12,3} = CL_auc_MEdgeBetw(2,:);
AUC_CLs{13,3} = CL_auc_MNodeBetw(2,:);

% Upper limit CL
AUC_CLs{1,4} = CL_auc_MClust(1,:);
AUC_CLs{2,4} = CL_auc_MDeg(1,:);
AUC_CLs{3,4} = CL_auc_Trans(1,:);
AUC_CLs{4,4} = CL_auc_Assort(1,:);
AUC_CLs{5,4} = CL_auc_Geff(1,:);
AUC_CLs{6,4} = CL_auc_Leff(1,:);
AUC_CLs{7,4} = CL_auc_Mod(1,:);
AUC_CLs{8,4} = CL_auc_PathL(1,:);
AUC_CLs{9,4} = CL_auc_Lambda(1,:);
AUC_CLs{10,4} = CL_auc_Gamma(1,:);
AUC_CLs{11,4} = CL_auc_Sigma(1,:);
AUC_CLs{12,4} = CL_auc_MEdgeBetw(1,:);
AUC_CLs{13,4} = CL_auc_MNodeBetw(1,:);

% p values
AUC_CLs{1,5} = CL_Pval(auc_MClust_2_rand - auc_MClust_1_rand, auc_MClust_2 - auc_MClust_1,'AUC_MClust',Tail);
AUC_CLs{2,5} = CL_Pval(auc_MDeg_2_rand - auc_MDeg_1_rand, auc_MDeg_2 - auc_MDeg_1,'AUC_MDeg',Tail);
AUC_CLs{3,5} = CL_Pval(auc_Trans_2_rand - auc_Trans_1_rand, auc_Trans_2 - auc_Trans_1,'AUC_Trans',Tail);
AUC_CLs{4,5} = CL_Pval(auc_Assort_2_rand - auc_Assort_1_rand, auc_Assort_2 - auc_Assort_1,'AUC_Assort',Tail);
AUC_CLs{5,5} = CL_Pval(auc_GEff_2_rand - auc_GEff_1_rand, auc_GEff_2 - auc_GEff_1,'AUC_GEff',Tail);
AUC_CLs{6,5} = CL_Pval(auc_MLocEff_2_rand - auc_MLocEff_1_rand, auc_MLocEff_2 - auc_MLocEff_1,'AUC_MLocEff',Tail);
AUC_CLs{7,5} = CL_Pval(auc_Mod_2_rand - auc_Mod_1_rand, auc_Mod_2 - auc_Mod_1,'AUC_Mod',Tail);
AUC_CLs{8,5} = CL_Pval(auc_PathL_2_rand - auc_PathL_1_rand, auc_PathL_2 - auc_PathL_1,'AUC_PathL',Tail);
AUC_CLs{9,5} = CL_Pval(auc_Lambda_2_rand - auc_Lambda_1_rand, auc_Lambda_2 - auc_Lambda_1,'AUC_Lambda',Tail);
AUC_CLs{10,5} = CL_Pval(auc_Gamma_2_rand - auc_Gamma_1_rand, auc_Gamma_2 - auc_Gamma_1,'AUC_Gamma',Tail);
AUC_CLs{11,5} = CL_Pval(auc_Sigma_2_rand - auc_Sigma_1_rand, auc_Sigma_2 - auc_Sigma_1,'AUC_Sigma',Tail);
AUC_CLs{12,5} = CL_Pval(auc_MEdgeBetw_2_rand - auc_MEdgeBetw_1_rand, auc_MEdgeBetw_2 - auc_MEdgeBetw_1,'AUC_MEdgeBetw',Tail);
AUC_CLs{13,5} = CL_Pval(auc_MNodeBetw_2_rand - auc_MNodeBetw_1_rand, auc_MNodeBetw_2 - auc_MNodeBetw_1,'AUC_MNodeBetw',Tail);

[h crit_p adj_p]=fdr_bh(abs(cell2mat(AUC_CLs(:,5))));
AUC_CLs(:,6)=num2cell(adj_p);

colnames = {'measure','meandiff','LL','UL','pval','pFDR'};
resultsTable = cell2table(AUC_CLs,'VariableNames',colnames);
save AUC_ResultsTable resultsTable
writetable(resultsTable,'AUCresultsTable.xlsx')

cd ..
fprintf('%-4s\n',' ..........Done ');