function [Com_g1,Coml_g1,Com_g2,Coml_g2] = Module_Analysis_all(data_log,Data1,Data2,nIter,MidThr)
%%%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao
%Ref (Hosseini et al.,Plos One 2012)


ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;ROI = ff.mat4GAT.roi1;
MinMesPlot=ff.mat4GAT.MinThr;
MaxMesPlot=ff.mat4GAT.MaxThr;
MesStepPlot=ff.mat4GAT.MesStep;

%% Group1 Vs Group2

f1=load(Data1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
f2=load(Data2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;

Xax = [MinMesPlot:MesStepPlot:MaxMesPlot];
MidIdx = find(single(Xax) == single(MidThr));

if isempty(MidIdx)
    error('no such target density');
end

j = MidIdx;
count =0;Q1=[];
for i = 1:size(NetMes1,1)
    count = count +1;
    Q1(:,:,count) = NetMes1{i,j}{24,1};
end

Q1=mean(Q1,3);
O1 = Q1;

count =0;Q2=[];
for i = 1:size(NetMes2,1)
    count = count +1;
    Q2(:,:,count) = NetMes2{i,j}{24,1};
end

Q2=mean(Q2,3);
O2 = Q2;

for i = 1:nIter
    [Cl_g1(i,:),Ql_g1(i)] = community_louvain(O1);
    Nmodl_g1(i) = max(Cl_g1(i,:));
    
    [Cl_g2(i,:),Ql_g2(i)] = community_louvain(O2);
    Nmodl_g2(i) = max(Cl_g2(i,:));
end

imaxl_g1 = find(Ql_g1==max(Ql_g1));
[All_Nmodl_g1, Idxl_g1] = max(Nmodl_g1(imaxl_g1));
Cl_g1 = Cl_g1(imaxl_g1(Idxl_g1),:);

imaxl_g2 = find(Ql_g2==max(Ql_g2));
[CommonModl_g2, Idxl_g2] = max(Nmodl_g2(imaxl_g2));
Cl_g2 = Cl_g2(imaxl_g2(Idxl_g2),:);

for i=1:max(Cl_g1)
    Coml_g1{i}=ROI(Cl_g1==i);
end

for i=1:max(Cl_g2)
    Coml_g2{i}=ROI(Cl_g2==i);
end

for i = 1:nIter
    [C_g1(i,:),Q_g1(i)] = modularity_und(O1);
    Nmod_g1(i) = max(C_g1(i,:));
    
    [C_g2(i,:),Q_g2(i)] = modularity_und(O2);
    Nmod_g2(i) = max(C_g2(i,:));
end

imax_g1 = find(Q_g1==max(Q_g1));
[All_Nmod_g1, Idx_g1] = max(Nmod_g1(imax_g1));
C_g1 = C_g1(imax_g1(Idx_g1),:);

imax_g2 = find(Q_g2==max(Q_g2));
[CommonMod_g2, Idx_g2] = max(Nmod_g2(imax_g2));
C_g2 = C_g2(imax_g2(Idx_g2),:);

for i=1:max(C_g1)
    Com_g1{i}=ROI(C_g1==i);
end

for i=1:max(C_g2)
    Com_g2{i}=ROI(C_g2==i);
end

save(['ModularCommunity_' Group1],'Com_g1','Coml_g1');
save(['ModularCommunity_' Group2],'Com_g2','Coml_g2');

fprintf('%-4s\n','.................done');