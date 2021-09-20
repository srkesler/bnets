function [N1,N2,N3,N4,N_G1,N_G2,N_G3,N_G4]=Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1_g1g2,Data2_g1g2,Data2_g2g3,Data3_g2g3, Data1_g1g3, Data3_g1g3,Data1_g1g4,Data4_g1g4,Data2_g2g4,Data4_g2g4,Data3_g3g4,Data4_g3g4,Nperm)

%data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
%MinMax: mninmum, maximum thresholds and the thresholded step
%Nperm: number of permutation (sampling) (default 100)

%output
%'N1': Original unthresholded correlation matrix for group1
%'N2': Original unthresholded correlation matrix for group2
%'N3': Original unthresholded correlation matrix for group3
%'N4': Original unthresholded correlation matrix for group4
%'N_G1': Cell array containing Randomly generated unthresholded correlation matrix for group1
%'N_G2': Cell array containing Randomly generated unthresholded correlation matrix for group2
%'N_G3': Cell array containing Randomly generated unthresholded correlation matrix for group3
%'N_G4': Cell array containing Randomly generated unthresholded correlation matrix for group4

%%

%-- Hadi Hosseini, Mar 21, 2011
%-- updated on Apr 20,2011 for GUI
%-- updated on August 18,2011 for functionals
%-- Shelli Kesler 10/30/15 corrected parallel processing
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%% Main body common to all group combinations

if isempty(Nperm)
    Nperm=100;
end

ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2; Group3 = ff.mat4GAT.g3; Group4 = ff.mat4GAT.g4;

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

%% Group1 vs Group2

f1=load(Data1_g1g2,'NetMes_Bin');data1=f1.NetMes_Bin;
f2=load(Data2_g1g2,'NetMes_Bin');data2=f2.NetMes_Bin;
N1=data1;N2=data2;

Sz1=size(data1,1);Sz2=size(data2,1);

data=data1;data(Sz1+1:Sz1+Sz2,:)=data2;
rng(1001,'twister')% Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results % edited 9/20/21 by 
%Shelli Kesler to address legacy random number generator override
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

N_G1=cell(1,Nperm);N_G2=cell(1,Nperm);


for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
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

save NetMesPerDens_G1_vs_G2 Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand Assort_2_rand GEff_2_rand ...
    MLocEff_2_rand Mod_2_rand ModL_2_rand PathL_2_rand MNodeBetw_2_rand ...
    MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand Sigma_2_rand ...
    Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
    Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ...
    ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2...
    
%% Group2 vs Group3

f2=load(Data2_g2g3,'NetMes_Bin');data2=f2.NetMes_Bin;
f3=load(Data3_g2g3,'NetMes_Bin');data3=f3.NetMes_Bin;
N2=data2; N3=data3;

Sz2=size(data2,1); Sz3=size(data3,1);

data=data2;data(Sz2+1:Sz2+Sz3,:)=data3;
rng(1001,'twister') % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz2+Sz3);
Randata(1:Sz2+Sz3,:)=data(RandIndex(1:Sz2+Sz3),:);

N_G2=cell(1,Nperm); N_G3=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp2=randsample(Sz2+Sz3,Sz2,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp3=randsample(Sz2+Sz3,Sz3,'true');
    N_G2{i}=Randata(Samp2,:);
    N_G3{i}=Randata(Samp3,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
MClust_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MDeg_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MStr_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Trans_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Assort_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
GEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MLocEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Mod_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));ModL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
PathL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MNodeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MEdgeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
Lambda_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Gamma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Sigma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));

Dens_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
MClust_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MDeg_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MStr_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Trans_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Assort_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
GEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MLocEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Mod_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));ModL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
PathL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MNodeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MEdgeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
Lambda_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Gamma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Sigma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));

Dens_2=zeros(size(N2));
MClust_2=zeros(size(N2));MDeg_2=zeros(size(N2));MStr_2=zeros(size(N2));Trans_2=zeros(size(N2));Assort_2=zeros(size(N2));
GEff_2=zeros(size(N2));MLocEff_2=zeros(size(N2));Mod_2=zeros(size(N2));ModL_2=zeros(size(N2));
PathL_2=zeros(size(N2));MNodeBetw_2=zeros(size(N2));MEdgeBetw_2=zeros(size(N2));
Lambda_2=zeros(size(N2));Gamma_2=zeros(size(N2));Sigma_2=zeros(size(N2));

Dens_3=zeros(size(N3));
MClust_3=zeros(size(N3));MDeg_3=zeros(size(N3));MStr_3=zeros(size(N3));Trans_3=zeros(size(N3));Assort_3=zeros(size(N3));
GEff_3=zeros(size(N3));MLocEff_3=zeros(size(N3));Mod_3=zeros(size(N3));ModL_3=zeros(size(N3));
PathL_3=zeros(size(N3));MNodeBetw_3=zeros(size(N3));MEdgeBetw_3=zeros(size(N3));
Lambda_3=zeros(size(N3));Gamma_3=zeros(size(N3));Sigma_3=zeros(size(N3));

for i=1:size(N_G2,2)
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
    
    for j=1:size(N3,1)
        for k=1:size(N3,2)
            Dens_3_rand(j,k,i)=N_G3{i}{j,k}{6};
            MClust_3_rand(j,k,i)=N_G3{i}{j,k}{8};
            MDeg_3_rand(j,k,i)=N_G3{i}{j,k}{2};
            Trans_3_rand(j,k,i)=N_G3{i}{j,k}{9};
            Assort_3_rand(j,k,i)=N_G3{i}{j,k}{5};
            GEff_3_rand(j,k,i)=N_G3{i}{j,k}{10};
            MLocEff_3_rand(j,k,i)=N_G3{i}{j,k}{12};
            Mod_3_rand(j,k,i)=N_G3{i}{j,k}{13};
            ModL_3_rand(j,k,i)=N_G3{i}{j,k}{14};
            PathL_3_rand(j,k,i)=N_G3{i}{j,k}{15};
            MNodeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{17};
            MEdgeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{19};
            Lambda_3_rand(j,k,i)=N_G3{i}{j,k}{20};
            Gamma_3_rand(j,k,i)=N_G3{i}{j,k}{21};
            Sigma_3_rand(j,k,i)=N_G3{i}{j,k}{22};
        end
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


for j=1:size(N3,1)
    for k=1:size(N3,2)
        Dens_3(j,k)=N3{j,k}{6};
        MClust_3(j,k)=N3{j,k}{8};
        MDeg_3(j,k)=N3{j,k}{2};
        Trans_3(j,k)=N3{j,k}{9};
        Assort_3(j,k)=N3{j,k}{5};
        GEff_3(j,k)=N3{j,k}{10};
        MLocEff_3(j,k)=N3{j,k}{12};
        Mod_3(j,k)=N3{j,k}{13};
        ModL_3(j,k)=N3{j,k}{14};
        PathL_3(j,k)=N3{j,k}{15};
        MNodeBetw_3(j,k)=N3{j,k}{17};
        MEdgeBetw_3(j,k)=N3{j,k}{19};
        Lambda_3(j,k)=N3{j,k}{20};
        Gamma_3(j,k)=N3{j,k}{21};
        Sigma_3(j,k)=N3{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens_G2_vs_G3 Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
    Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
    PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
    Sigma_2_rand ...
    Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 ...
    Mod_2 ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 ...
    Dens_3_rand MClust_3_rand MDeg_3_rand Trans_3_rand ...
    Assort_3_rand GEff_3_rand MLocEff_3_rand Mod_3_rand ModL_3_rand ...
    PathL_3_rand MNodeBetw_3_rand MEdgeBetw_3_rand Lambda_3_rand Gamma_3_rand ...
    Sigma_3_rand Dens_3 MClust_3 MDeg_3 Trans_3 Assort_3 GEff_3 MLocEff_3 ...
    Mod_3 ModL_3 PathL_3 MNodeBetw_3 MEdgeBetw_3 Lambda_3 Gamma_3 Sigma_3


%% Group1 vs Group3
f1=load(Data1_g1g3,'NetMes_Bin');data1=f1.NetMes_Bin;
f3=load(Data3_g1g3,'NetMes_Bin');data3=f3.NetMes_Bin;
N1=data1; N3=data3;

Sz1=size(data1,1);Sz3=size(data3,1);

data=data1;data(Sz1+1:Sz1+Sz3,:)=data3;
rng(1001,'twister') % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz1+Sz3);
Randata(1:Sz1+Sz3,:)=data(RandIndex(1:Sz1+Sz3),:);

N_G1=cell(1,Nperm); N_G3=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz3,Sz1,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp3=randsample(Sz1+Sz3,Sz3,'true');
    N_G1{i}=Randata(Samp1,:);
    N_G3{i}=Randata(Samp3,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));%[];
MClust_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MDeg_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MStr_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Trans_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Assort_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
GEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MLocEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Mod_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));ModL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
PathL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MNodeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MEdgeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
Lambda_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Gamma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Sigma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));

Dens_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
MClust_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MDeg_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MStr_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Trans_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Assort_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
GEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MLocEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Mod_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));ModL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
PathL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MNodeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));MEdgeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));
Lambda_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Gamma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));Sigma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G3,2));

Dens_1=zeros(size(N1));
MClust_1=zeros(size(N1));MDeg_1=zeros(size(N1));MStr_1=zeros(size(N1));Trans_1=zeros(size(N1));Assort_1=zeros(size(N1));
GEff_1=zeros(size(N1));MLocEff_1=zeros(size(N1));Mod_1=zeros(size(N1));ModL_1=zeros(size(N1));
PathL_1=zeros(size(N1));MNodeBetw_1=zeros(size(N1));MEdgeBetw_1=zeros(size(N1));
Lambda_1=zeros(size(N1));Gamma_1=zeros(size(N1));Sigma_1=zeros(size(N1));

Dens_3=zeros(size(N3));
MClust_3=zeros(size(N3));MDeg_3=zeros(size(N3));MStr_3=zeros(size(N3));Trans_3=zeros(size(N3));Assort_3=zeros(size(N3));
GEff_3=zeros(size(N3));MLocEff_3=zeros(size(N3));Mod_3=zeros(size(N3));ModL_3=zeros(size(N3));
PathL_3=zeros(size(N3));MNodeBetw_3=zeros(size(N3));MEdgeBetw_3=zeros(size(N3));
Lambda_3=zeros(size(N3));Gamma_3=zeros(size(N3));Sigma_3=zeros(size(N3));

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
    
    for j=1:size(N3,1)
        for k=1:size(N3,2)
            Dens_3_rand(j,k,i)=N_G3{i}{j,k}{6};
            MClust_3_rand(j,k,i)=N_G3{i}{j,k}{8};
            MDeg_3_rand(j,k,i)=N_G3{i}{j,k}{2};
            Trans_3_rand(j,k,i)=N_G3{i}{j,k}{9};
            Assort_3_rand(j,k,i)=N_G3{i}{j,k}{5};
            GEff_3_rand(j,k,i)=N_G3{i}{j,k}{10};
            MLocEff_3_rand(j,k,i)=N_G3{i}{j,k}{12};
            Mod_3_rand(j,k,i)=N_G3{i}{j,k}{13};
            ModL_3_rand(j,k,i)=N_G3{i}{j,k}{14};
            PathL_3_rand(j,k,i)=N_G3{i}{j,k}{15};
            MNodeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{17};
            MEdgeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{19};
            Lambda_3_rand(j,k,i)=N_G3{i}{j,k}{20};
            Gamma_3_rand(j,k,i)=N_G3{i}{j,k}{21};
            Sigma_3_rand(j,k,i)=N_G3{i}{j,k}{22};
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

for j=1:size(N3,1)
    for k=1:size(N3,2)
        Dens_3(j,k)=N3{j,k}{6};
        MClust_3(j,k)=N3{j,k}{8};
        MDeg_3(j,k)=N3{j,k}{2};
        Trans_3(j,k)=N3{j,k}{9};
        Assort_3(j,k)=N3{j,k}{5};
        GEff_3(j,k)=N3{j,k}{10};
        MLocEff_3(j,k)=N3{j,k}{12};
        Mod_3(j,k)=N3{j,k}{13};
        ModL_3(j,k)=N3{j,k}{14};
        PathL_3(j,k)=N3{j,k}{15};
        MNodeBetw_3(j,k)=N3{j,k}{17};
        MEdgeBetw_3(j,k)=N3{j,k}{19};
        Lambda_3(j,k)=N3{j,k}{20};
        Gamma_3(j,k)=N3{j,k}{21};
        Sigma_3(j,k)=N3{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens_G1_vs_G3 Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand ...
    Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
    Dens_3_rand MClust_3_rand MDeg_3_rand Trans_3_rand ...
    Assort_3_rand GEff_3_rand MLocEff_3_rand Mod_3_rand ModL_3_rand ...
    PathL_3_rand MNodeBetw_3_rand MEdgeBetw_3_rand Lambda_3_rand Gamma_3_rand ...
    Sigma_3_rand Dens_3 MClust_3 MDeg_3 Trans_3 Assort_3 GEff_3 MLocEff_3 ...
    Mod_3 ModL_3 PathL_3 MNodeBetw_3 MEdgeBetw_3 Lambda_3 Gamma_3 Sigma_3


%% Group1 vs Group4
f1=load(Data1_g1g4,'NetMes_Bin');data1=f1.NetMes_Bin;
f4=load(Data4_g1g4,'NetMes_Bin');data4=f4.NetMes_Bin;
N1=data1; N4=data4;

Sz1=size(data1,1);Sz4=size(data4,1);

data=data1;data(Sz1+1:Sz1+Sz4,:)=data4;
rng(1001,'twister') % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz1+Sz4);
Randata(1:Sz1+Sz4,:)=data(RandIndex(1:Sz1+Sz4),:);

N_G1=cell(1,Nperm); N_G4=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp1=randsample(Sz1+Sz4,Sz1,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/17/2018 in order to fix the seed and have consistent results
    Samp4=randsample(Sz1+Sz4,Sz4,'true');
    N_G1{i}=Randata(Samp1,:);
    N_G4{i}=Randata(Samp4,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));%[];
MClust_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MDeg_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MStr_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Trans_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Assort_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
GEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MLocEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Mod_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));ModL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
PathL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MNodeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MEdgeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
Lambda_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Gamma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Sigma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));

Dens_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
MClust_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MDeg_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MStr_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Trans_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Assort_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
GEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MLocEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Mod_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));ModL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
PathL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MNodeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MEdgeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
Lambda_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Gamma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Sigma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));

Dens_1=zeros(size(N1));
MClust_1=zeros(size(N1));MDeg_1=zeros(size(N1));MStr_1=zeros(size(N1));Trans_1=zeros(size(N1));Assort_1=zeros(size(N1));
GEff_1=zeros(size(N1));MLocEff_1=zeros(size(N1));Mod_1=zeros(size(N1));ModL_1=zeros(size(N1));
PathL_1=zeros(size(N1));MNodeBetw_1=zeros(size(N1));MEdgeBetw_1=zeros(size(N1));
Lambda_1=zeros(size(N1));Gamma_1=zeros(size(N1));Sigma_1=zeros(size(N1));

Dens_4=zeros(size(N4));
MClust_4=zeros(size(N4));MDeg_4=zeros(size(N4));MStr_4=zeros(size(N4));Trans_4=zeros(size(N4));Assort_4=zeros(size(N4));
GEff_4=zeros(size(N4));MLocEff_4=zeros(size(N4));Mod_4=zeros(size(N4));ModL_4=zeros(size(N4));
PathL_4=zeros(size(N4));MNodeBetw_4=zeros(size(N4));MEdgeBetw_4=zeros(size(N4));
Lambda_4=zeros(size(N4));Gamma_4=zeros(size(N4));Sigma_4=zeros(size(N4));


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
    
    for j=1:size(N4,1)
        for k=1:size(N4,2)
            Dens_4_rand(j,k,i)=N_G4{i}{j,k}{6};
            MClust_4_rand(j,k,i)=N_G4{i}{j,k}{8};
            MDeg_4_rand(j,k,i)=N_G4{i}{j,k}{2};
            Trans_4_rand(j,k,i)=N_G4{i}{j,k}{9};
            Assort_4_rand(j,k,i)=N_G4{i}{j,k}{5};
            GEff_4_rand(j,k,i)=N_G4{i}{j,k}{10};
            MLocEff_4_rand(j,k,i)=N_G4{i}{j,k}{12};
            Mod_4_rand(j,k,i)=N_G4{i}{j,k}{13};
            ModL_4_rand(j,k,i)=N_G4{i}{j,k}{14};
            PathL_4_rand(j,k,i)=N_G4{i}{j,k}{15};
            MNodeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{17};
            MEdgeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{19};
            Lambda_4_rand(j,k,i)=N_G4{i}{j,k}{20};
            Gamma_4_rand(j,k,i)=N_G4{i}{j,k}{21};
            Sigma_4_rand(j,k,i)=N_G4{i}{j,k}{22};
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

for j=1:size(N4,1)
    for k=1:size(N4,2)
        Dens_4(j,k)=N4{j,k}{6};
        MClust_4(j,k)=N4{j,k}{8};
        MDeg_4(j,k)=N4{j,k}{2};
        Trans_4(j,k)=N4{j,k}{9};
        Assort_4(j,k)=N4{j,k}{5};
        GEff_4(j,k)=N4{j,k}{10};
        MLocEff_4(j,k)=N4{j,k}{12};
        Mod_4(j,k)=N4{j,k}{13};
        ModL_4(j,k)=N4{j,k}{14};
        PathL_4(j,k)=N4{j,k}{15};
        MNodeBetw_4(j,k)=N4{j,k}{17};
        MEdgeBetw_4(j,k)=N4{j,k}{19};
        Lambda_4(j,k)=N4{j,k}{20};
        Gamma_4(j,k)=N4{j,k}{21};
        Sigma_4(j,k)=N4{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens_G1_vs_G4 Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
    Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
    PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
    Sigma_1_rand ...
    Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
    Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
    Dens_4_rand MClust_4_rand MDeg_4_rand Trans_4_rand ...
    Assort_4_rand GEff_4_rand MLocEff_4_rand Mod_4_rand ModL_4_rand ...
    PathL_4_rand MNodeBetw_4_rand MEdgeBetw_4_rand Lambda_4_rand Gamma_4_rand ...
    Sigma_4_rand Dens_4 MClust_4 MDeg_4 Trans_4 Assort_4 GEff_4 MLocEff_4 ...
    Mod_4 ModL_4 PathL_4 MNodeBetw_4 MEdgeBetw_4 Lambda_4 Gamma_4 Sigma_4


%% Group2 vs Group4

f2=load(Data2_g2g4,'NetMes_Bin');data2=f2.NetMes_Bin;
f4=load(Data4_g2g4,'NetMes_Bin');data4=f4.NetMes_Bin;
N2=data2; N4=data4;

Sz2=size(data2,1); Sz4=size(data4,1);

data=data2;data(Sz2+1:Sz2+Sz4,:)=data4;
rng(1001,'twister') % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz2+Sz4);
Randata(1:Sz2+Sz4,:)=data(RandIndex(1:Sz2+Sz4),:);

N_G2=cell(1,Nperm); N_G4=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
    Samp2=randsample(Sz2+Sz4,Sz2,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
    Samp4=randsample(Sz2+Sz4,Sz4,'true');
    N_G2{i}=Randata(Samp2,:);
    N_G4{i}=Randata(Samp4,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
MClust_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MDeg_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MStr_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Trans_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Assort_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
GEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MLocEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Mod_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));ModL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
PathL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MNodeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MEdgeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
Lambda_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Gamma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Sigma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));

Dens_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
MClust_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MDeg_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MStr_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Trans_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Assort_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
GEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MLocEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Mod_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));ModL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
PathL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MNodeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MEdgeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
Lambda_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Gamma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Sigma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));

Dens_2=zeros(size(N2));
MClust_2=zeros(size(N2));MDeg_2=zeros(size(N2));MStr_2=zeros(size(N2));Trans_2=zeros(size(N2));Assort_2=zeros(size(N2));
GEff_2=zeros(size(N2));MLocEff_2=zeros(size(N2));Mod_2=zeros(size(N2));ModL_2=zeros(size(N2));
PathL_2=zeros(size(N2));MNodeBetw_2=zeros(size(N2));MEdgeBetw_2=zeros(size(N2));
Lambda_2=zeros(size(N2));Gamma_2=zeros(size(N2));Sigma_2=zeros(size(N2));

Dens_4=zeros(size(N4));
MClust_4=zeros(size(N4));MDeg_4=zeros(size(N4));MStr_4=zeros(size(N4));Trans_4=zeros(size(N4));Assort_4=zeros(size(N4));
GEff_4=zeros(size(N4));MLocEff_4=zeros(size(N4));Mod_4=zeros(size(N4));ModL_4=zeros(size(N4));
PathL_4=zeros(size(N4));MNodeBetw_4=zeros(size(N4));MEdgeBetw_4=zeros(size(N4));
Lambda_4=zeros(size(N4));Gamma_4=zeros(size(N4));Sigma_4=zeros(size(N4));


for i=1:size(N_G2,2)
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
    
    for j=1:size(N4,1)
        for k=1:size(N4,2)
            Dens_4_rand(j,k,i)=N_G4{i}{j,k}{6};
            MClust_4_rand(j,k,i)=N_G4{i}{j,k}{8};
            MDeg_4_rand(j,k,i)=N_G4{i}{j,k}{2};
            Trans_4_rand(j,k,i)=N_G4{i}{j,k}{9};
            Assort_4_rand(j,k,i)=N_G4{i}{j,k}{5};
            GEff_4_rand(j,k,i)=N_G4{i}{j,k}{10};
            MLocEff_4_rand(j,k,i)=N_G4{i}{j,k}{12};
            Mod_4_rand(j,k,i)=N_G4{i}{j,k}{13};
            ModL_4_rand(j,k,i)=N_G4{i}{j,k}{14};
            PathL_4_rand(j,k,i)=N_G4{i}{j,k}{15};
            MNodeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{17};
            MEdgeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{19};
            Lambda_4_rand(j,k,i)=N_G4{i}{j,k}{20};
            Gamma_4_rand(j,k,i)=N_G4{i}{j,k}{21};
            Sigma_4_rand(j,k,i)=N_G4{i}{j,k}{22};
        end
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

for j=1:size(N4,1)
    for k=1:size(N4,2)
        Dens_4(j,k)=N4{j,k}{6};
        MClust_4(j,k)=N4{j,k}{8};
        MDeg_4(j,k)=N4{j,k}{2};
        Trans_4(j,k)=N4{j,k}{9};
        Assort_4(j,k)=N4{j,k}{5};
        GEff_4(j,k)=N4{j,k}{10};
        MLocEff_4(j,k)=N4{j,k}{12};
        Mod_4(j,k)=N4{j,k}{13};
        ModL_4(j,k)=N4{j,k}{14};
        PathL_4(j,k)=N4{j,k}{15};
        MNodeBetw_4(j,k)=N4{j,k}{17};
        MEdgeBetw_4(j,k)=N4{j,k}{19};
        Lambda_4(j,k)=N4{j,k}{20};
        Gamma_4(j,k)=N4{j,k}{21};
        Sigma_4(j,k)=N4{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens_G2_vs_G4 Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
    Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
    PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
    Sigma_2_rand ...
    Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 ...
    Mod_2 ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 ...
    Dens_4_rand MClust_4_rand MDeg_4_rand Trans_4_rand ...
    Assort_4_rand GEff_4_rand MLocEff_4_rand Mod_4_rand ModL_4_rand ...
    PathL_4_rand MNodeBetw_4_rand MEdgeBetw_4_rand Lambda_4_rand Gamma_4_rand ...
    Sigma_4_rand Dens_4 MClust_4 MDeg_4 Trans_4 Assort_4 GEff_4 MLocEff_4 ...
    Mod_4 ModL_4 PathL_4 MNodeBetw_4 MEdgeBetw_4 Lambda_4 Gamma_4 Sigma_4


%% Group3 vs Group4

f3=load(Data3_g3g4,'NetMes_Bin');data3=f3.NetMes_Bin;
f4=load(Data4_g3g4,'NetMes_Bin');data4=f4.NetMes_Bin;
N3=data3; N4=data4;

Sz3=size(data3,1); Sz4=size(data4,1);

data=data3;data(Sz3+1:Sz3+Sz4,:)=data4;
rng(1001,'twister') % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
RandIndex=randperm(Sz3+Sz4);
Randata(1:Sz3+Sz4,:)=data(RandIndex(1:Sz3+Sz4),:);

N_G3=cell(1,Nperm); N_G4=cell(1,Nperm);

for i=1:Nperm
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    rng(1000+i,'twister'); % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
    Samp3=randsample(Sz3+Sz4,Sz3,'true');
    rng(2000+i,'twister'); % Added by Vikram Rao on 09/18/2018 in order to fix the seed and have consistent results
    Samp4=randsample(Sz3+Sz4,Sz4,'true');
    N_G3{i}=Randata(Samp3,:);
    N_G4{i}=Randata(Samp4,:);
end

fprintf('%-4s\n',' extracting network measures...');

Dens_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));
MClust_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));MDeg_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));MStr_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));Trans_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));Assort_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));
GEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));MLocEff_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));Mod_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));ModL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));
PathL_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));MNodeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));MEdgeBetw_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));
Lambda_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));Gamma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));Sigma_3_rand=zeros(size(N3,1),size(N3,2),size(N_G2,2));

Dens_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
MClust_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MDeg_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MStr_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Trans_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Assort_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
GEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MLocEff_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Mod_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));ModL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
PathL_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MNodeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));MEdgeBetw_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));
Lambda_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Gamma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));Sigma_4_rand=zeros(size(N4,1),size(N4,2),size(N_G4,2));

Dens_3=zeros(size(N3));
MClust_3=zeros(size(N3));MDeg_3=zeros(size(N3));MStr_3=zeros(size(N3));Trans_3=zeros(size(N3));Assort_3=zeros(size(N3));
GEff_3=zeros(size(N3));MLocEff_3=zeros(size(N3));Mod_3=zeros(size(N3));ModL_3=zeros(size(N3));
PathL_3=zeros(size(N3));MNodeBetw_3=zeros(size(N3));MEdgeBetw_3=zeros(size(N3));
Lambda_3=zeros(size(N3));Gamma_3=zeros(size(N3));Sigma_3=zeros(size(N3));

Dens_4=zeros(size(N4));
MClust_4=zeros(size(N4));MDeg_4=zeros(size(N4));MStr_4=zeros(size(N4));Trans_4=zeros(size(N4));Assort_4=zeros(size(N4));
GEff_4=zeros(size(N4));MLocEff_4=zeros(size(N4));Mod_4=zeros(size(N4));ModL_4=zeros(size(N4));
PathL_4=zeros(size(N4));MNodeBetw_4=zeros(size(N4));MEdgeBetw_4=zeros(size(N4));
Lambda_4=zeros(size(N4));Gamma_4=zeros(size(N4));Sigma_4=zeros(size(N4));


for i=1:size(N_G3,2)
    for j=1:size(N3,1)
        for k=1:size(N3,2)
            Dens_3_rand(j,k,i)=N_G3{i}{j,k}{6};
            MClust_3_rand(j,k,i)=N_G3{i}{j,k}{8};
            MDeg_3_rand(j,k,i)=N_G3{i}{j,k}{2};
            Trans_3_rand(j,k,i)=N_G3{i}{j,k}{9};
            Assort_3_rand(j,k,i)=N_G3{i}{j,k}{5};
            GEff_3_rand(j,k,i)=N_G3{i}{j,k}{10};
            MLocEff_3_rand(j,k,i)=N_G3{i}{j,k}{12};
            Mod_3_rand(j,k,i)=N_G3{i}{j,k}{13};
            ModL_3_rand(j,k,i)=N_G3{i}{j,k}{14};
            PathL_3_rand(j,k,i)=N_G3{i}{j,k}{15};
            MNodeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{17};
            MEdgeBetw_3_rand(j,k,i)=N_G3{i}{j,k}{19};
            Lambda_3_rand(j,k,i)=N_G3{i}{j,k}{20};
            Gamma_3_rand(j,k,i)=N_G3{i}{j,k}{21};
            Sigma_3_rand(j,k,i)=N_G3{i}{j,k}{22};
        end
    end
    
    for j=1:size(N4,1)
        for k=1:size(N4,2)
            Dens_4_rand(j,k,i)=N_G4{i}{j,k}{6};
            MClust_4_rand(j,k,i)=N_G4{i}{j,k}{8};
            MDeg_4_rand(j,k,i)=N_G4{i}{j,k}{2};
            Trans_4_rand(j,k,i)=N_G4{i}{j,k}{9};
            Assort_4_rand(j,k,i)=N_G4{i}{j,k}{5};
            GEff_4_rand(j,k,i)=N_G4{i}{j,k}{10};
            MLocEff_4_rand(j,k,i)=N_G4{i}{j,k}{12};
            Mod_4_rand(j,k,i)=N_G4{i}{j,k}{13};
            ModL_4_rand(j,k,i)=N_G4{i}{j,k}{14};
            PathL_4_rand(j,k,i)=N_G4{i}{j,k}{15};
            MNodeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{17};
            MEdgeBetw_4_rand(j,k,i)=N_G4{i}{j,k}{19};
            Lambda_4_rand(j,k,i)=N_G4{i}{j,k}{20};
            Gamma_4_rand(j,k,i)=N_G4{i}{j,k}{21};
            Sigma_4_rand(j,k,i)=N_G4{i}{j,k}{22};
        end
    end
end

for j=1:size(N3,1)
    for k=1:size(N3,2)
        Dens_3(j,k)=N3{j,k}{6};
        MClust_3(j,k)=N3{j,k}{8};
        MDeg_3(j,k)=N3{j,k}{2};
        Trans_3(j,k)=N3{j,k}{9};
        Assort_3(j,k)=N3{j,k}{5};
        GEff_3(j,k)=N3{j,k}{10};
        MLocEff_3(j,k)=N3{j,k}{12};
        Mod_3(j,k)=N3{j,k}{13};
        ModL_3(j,k)=N3{j,k}{14};
        PathL_3(j,k)=N3{j,k}{15};
        MNodeBetw_3(j,k)=N3{j,k}{17};
        MEdgeBetw_3(j,k)=N3{j,k}{19};
        Lambda_3(j,k)=N3{j,k}{20};
        Gamma_3(j,k)=N3{j,k}{21};
        Sigma_3(j,k)=N3{j,k}{22};
    end
end

for j=1:size(N4,1)
    for k=1:size(N4,2)
        Dens_4(j,k)=N4{j,k}{6};
        MClust_4(j,k)=N4{j,k}{8};
        MDeg_4(j,k)=N4{j,k}{2};
        Trans_4(j,k)=N4{j,k}{9};
        Assort_4(j,k)=N4{j,k}{5};
        GEff_4(j,k)=N4{j,k}{10};
        MLocEff_4(j,k)=N4{j,k}{12};
        Mod_4(j,k)=N4{j,k}{13};
        ModL_4(j,k)=N4{j,k}{14};
        PathL_4(j,k)=N4{j,k}{15};
        MNodeBetw_4(j,k)=N4{j,k}{17};
        MEdgeBetw_4(j,k)=N4{j,k}{19};
        Lambda_4(j,k)=N4{j,k}{20};
        Gamma_4(j,k)=N4{j,k}{21};
        Sigma_4(j,k)=N4{j,k}{22};
    end
end

fprintf('%-4s\n',' saving output...');

save NetMesPerDens_G3_vs_G4 Dens_3_rand MClust_3_rand MDeg_3_rand Trans_3_rand ...
    Assort_3_rand GEff_3_rand MLocEff_3_rand Mod_3_rand ModL_3_rand ...
    PathL_3_rand MNodeBetw_3_rand MEdgeBetw_3_rand Lambda_3_rand Gamma_3_rand ...
    Sigma_3_rand ...
    Dens_3 MClust_3 MDeg_3 Trans_3 Assort_3 GEff_3 MLocEff_3 ...
    Mod_3 ModL_3 PathL_3 MNodeBetw_3 MEdgeBetw_3 Lambda_3 Gamma_3 Sigma_3 ...
    Dens_4_rand MClust_4_rand MDeg_4_rand Trans_4_rand ...
    Assort_4_rand GEff_4_rand MLocEff_4_rand Mod_4_rand ModL_4_rand ...
    PathL_4_rand MNodeBetw_4_rand MEdgeBetw_4_rand Lambda_4_rand Gamma_4_rand ...
    Sigma_4_rand Dens_4 MClust_4 MDeg_4 Trans_4 Assort_4 GEff_4 MLocEff_4 ...
    Mod_4 ModL_4 PathL_4 MNodeBetw_4 MEdgeBetw_4 Lambda_4 Gamma_4 Sigma_4

fprintf('%-4s\n','......Done ...');