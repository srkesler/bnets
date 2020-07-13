function [N1,N2,N_G1,N_G2]=Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,Nperm)

%data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
%MinMax: mninmum, maximum thresholds and the thresholded step
%Nperm: number of permutation (sampling) (default 100)

%output
%'N1': Original unthresholded correlation matrix for group1
%'N2': Original unthresholded correlation matrix for group2
%'N_G1': Cell array containing Randomly generated unthresholded correlation matrix for group1
%'N_G2': Cell array containing Randomly generated unthresholded correlation matrix for group2

%%

%-- Hadi Hosseini, Mar 21, 2011
%-- updated on Apr 20,2011 for GUI
%-- updated on August 18,2011 for functionals
%-- updated on Apr 2014 (parallel processing)
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

% Shelli Kesler 10/30/15 corrected parallel processing


%% Group1 vs Group2

if isempty(Nperm)
    Nperm=100;
end

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
f1=load(Data1,'NetMes_Bin');data1=f1.NetMes_Bin;
f2=load(Data2,'NetMes_Bin');data2=f2.NetMes_Bin;

MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;
N1=data1;N2=data2;

%% 

Sz1=size(data1,1);Sz2=size(data2,1);

data=data1;data(Sz1+1:Sz1+Sz2,:)=data2;
rand("seed",1001); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

N_G1=cell(1,Nperm);N_G2=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rand("seed",1000+i); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    rand("seed",2000+i); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    N_G1{i}=Randata(Samp1,:);
    N_G2{i}=Randata(Samp2,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
MClust_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MDeg_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MStr_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Trans_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Assort_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
GEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MLocEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Mod_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));ModL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
PathL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MNodeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MEdgeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
Lambda_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Gamma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Sigma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));

Dens_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
MClust_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MDeg_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MStr_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Trans_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Assort_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
GEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MLocEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Mod_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));ModL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
PathL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MNodeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MEdgeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
Lambda_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Gamma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Sigma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));

Dens_1=zeros(size(N1));
MClust_1=zeros(size(N1));MDeg_1=zeros(size(N1));MStr_1=zeros(size(N1));Trans_1=zeros(size(N1));Assort_1=zeros(size(N1));
GEff_1=zeros(size(N1));MLocEff_1=zeros(size(N1));Mod_1=zeros(size(N1));ModL_1=zeros(size(N1));
PathL_1=zeros(size(N1));MNodeBetw_1=zeros(size(N1));MEdgeBetw_1=zeros(size(N1));
Lambda_1=zeros(size(N1));Gamma_1=zeros(size(N1));Sigma_1=zeros(size(N1));

Dens_2=zeros(size(N2));
MClust_2=zeros(size(N2));MDeg_2=zeros(size(N2));MStr_2=zeros(size(N2));Trans_2=zeros(size(N2));Assort_2=zeros(size(N2));
GEff_2=zeros(size(N2));MLocEff_2=zeros(size(N2));Mod_2=zeros(size(N2));ModL_2=zeros(size(N2));
PathL_2=zeros(size(N2));MNodeBetw_2=zeros(size(N2));MEdgeBetw_2=zeros(size(N2));
Lambda_2=zeros(size(N2));Gamma_2=zeros(size(N2));Sigma_2=zeros(size(N2));

for i=1:size(N_G1,2)
    for j=1:size(N1,1)
        for k=1:size(N1,2)
            Dens_1_rand(j,k,i)=N_G1{i}{j,k}{6};
            MClust_1_rand(j,k,i)=N_G1{i}{j,k}{8};
            MDeg_1_rand(j,k,i)=N_G1{i}{j,k}{2};
            Trans_1_rand(j,k,i)=N_G1{i}{j,k}{9};
            Assort_1_rand(j,k,i)=N_G1{i}{j,k}{5};
            GEff_1_rand(j,k,i)=N_G1{i}{j,k}{10};
            MLocEff_1_rand(j,k,i)=N_G1{i}{j,k}{12};
            Mod_1_rand(j,k,i)=N_G1{i}{j,k}{13};
            ModL_1_rand(j,k,i)=N_G1{i}{j,k}{14};
            PathL_1_rand(j,k,i)=N_G1{i}{j,k}{15};
            MNodeBetw_1_rand(j,k,i)=N_G1{i}{j,k}{17};
            MEdgeBetw_1_rand(j,k,i)=N_G1{i}{j,k}{19};
            Lambda_1_rand(j,k,i)=N_G1{i}{j,k}{20};
            Gamma_1_rand(j,k,i)=N_G1{i}{j,k}{21};
            Sigma_1_rand(j,k,i)=N_G1{i}{j,k}{22};
        end
    end
    
    for j=1:size(N2,1)
        for k=1:size(N2,2)
            Dens_2_rand(j,k,i)=N_G2{i}{j,k}{6};
            MClust_2_rand(j,k,i)=N_G2{i}{j,k}{8};
            MDeg_2_rand(j,k,i)=N_G2{i}{j,k}{2};
            Trans_2_rand(j,k,i)=N_G2{i}{j,k}{9};
            Assort_2_rand(j,k,i)=N_G2{i}{j,k}{5};
            GEff_2_rand(j,k,i)=N_G2{i}{j,k}{10};
            MLocEff_2_rand(j,k,i)=N_G2{i}{j,k}{12};
            Mod_2_rand(j,k,i)=N_G2{i}{j,k}{13};
            ModL_2_rand(j,k,i)=N_G2{i}{j,k}{14};
            PathL_2_rand(j,k,i)=N_G2{i}{j,k}{15};
            MNodeBetw_2_rand(j,k,i)=N_G2{i}{j,k}{17};
            MEdgeBetw_2_rand(j,k,i)=N_G2{i}{j,k}{19};
            Lambda_2_rand(j,k,i)=N_G2{i}{j,k}{20};
            Gamma_2_rand(j,k,i)=N_G2{i}{j,k}{21};
            Sigma_2_rand(j,k,i)=N_G2{i}{j,k}{22}; 
        end
    end
end

for j=1:size(N1,1)
    for k=1:size(N1,2)
        Dens_1(j,k)=N1{j,k}{6};
        MClust_1(j,k)=N1{j,k}{8};
        MDeg_1(j,k)=N1{j,k}{2};
        Trans_1(j,k)=N1{j,k}{9};
        Assort_1(j,k)=N1{j,k}{5};
        GEff_1(j,k)=N1{j,k}{10};
        MLocEff_1(j,k)=N1{j,k}{12};
        Mod_1(j,k)=N1{j,k}{13};
        ModL_1(j,k)=N1{j,k}{14};
        PathL_1(j,k)=N1{j,k}{15};
        MNodeBetw_1(j,k)=N1{j,k}{17};
        MEdgeBetw_1(j,k)=N1{j,k}{19};
        Lambda_1(j,k)=N1{j,k}{20};
        Gamma_1(j,k)=N1{j,k}{21};
        Sigma_1(j,k)=N1{j,k}{22};
    end
end

for j=1:size(N2,1)
    for k=1:size(N2,2)
        Dens_2(j,k)=N2{j,k}{6};
        MClust_2(j,k)=N2{j,k}{8};
        MDeg_2(j,k)=N2{j,k}{2};
        Trans_2(j,k)=N2{j,k}{9};
        Assort_2(j,k)=N2{j,k}{5};
        GEff_2(j,k)=N2{j,k}{10};
        MLocEff_2(j,k)=N2{j,k}{12};
        Mod_2(j,k)=N2{j,k}{13};
        ModL_2(j,k)=N2{j,k}{14};
        PathL_2(j,k)=N2{j,k}{15};
        MNodeBetw_2(j,k)=N2{j,k}{17};
        MEdgeBetw_2(j,k)=N2{j,k}{19};
        Lambda_2(j,k)=N2{j,k}{20};
        Gamma_2(j,k)=N2{j,k}{21};
        Sigma_2(j,k)=N2{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens.mat -v7 Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
     Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
     PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
     Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand Assort_2_rand GEff_2_rand ...
     MLocEff_2_rand Mod_2_rand ModL_2_rand PathL_2_rand MNodeBetw_2_rand ...
     MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand Sigma_2_rand ...
     Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
     Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
     Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ...
     ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2;

fprintf('%-4s\n',' Done ...');