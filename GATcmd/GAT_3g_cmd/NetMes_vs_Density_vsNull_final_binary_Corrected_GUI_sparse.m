function NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(Null_Iter, Tail, Alpha, Precision, connfile,MinMes,MaxMes,MesStep,Group1,Group2, Group3, n1,n2,n3)  

% The program accepts the Results_ROI (output from conn) of all groups combined
%'MinMes': The minimum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.05 for density)
%'MaxMes': The maximum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.6 for density)


%%
%-- Hadi Hosseini March 22, 2011
%-- Hadi Hosseini August 18, 2011 (updated for functional connectivity)
%-- Hadi Hosseini August 18, 2011 (updated for covariates of no interest)
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao


%% Main body common to all group combinations

%Default Precision=10;
Step=0.0001;

filename = connfile;
fprintf('%-4s\n',' loading the input data....');
load(filename,'Z','names')
[nROI,nROI2]=size(Z);

MinMesPlot=MinMes;
MaxMesPlot=MaxMes;
MesStepPlot=MesStep;
MinThr = MinMesPlot;
MaxThr = MaxMesPlot;

mat4GAT.g1 = Group1;mat4GAT.g2 = Group2; mat4GAT.g3 = Group3;
mat4GAT.n1 = n1;mat4GAT.n2 = n2; mat4GAT.n3 = n3;
mat4GAT.MinThr = MinThr;
mat4GAT.MaxThr = MaxThr;
mat4GAT.MesStep = MesStepPlot;
mat4GAT.roi1 = names;
mat4GAT.roi2 = names;
mat4GAT.roi3 = names;

%null_flag = input('type 2 (same deg dist), 3 (H-Q-S), 4 (RandCorr): '); 
null_flag = 2;%at this moment, only upports same deg dist
mat4GAT.null = null_flag;
NullType = null_flag;

% % % Null_Iter = input('number of null networks (e.g. 20):  ');
mat4GAT.nullIter = Null_Iter;

% % % Tail = input('one (type 1) or two (type 2) -tailed test:  ');
mat4GAT.tail = Tail;

% % % Alpha = input('significance threshold level (e.g. 0.05):  ');
mat4GAT.alpha = Alpha;


save(['mat4GAT_' Group1 '_' Group2 '_' Group3],'mat4GAT');

%%

N = size(Z,1);
MinEdge=MinMes*((N^2-N)/2);MinMes=floor(MinEdge);
MaxEdge=MaxMes*((N^2-N)/2);MaxMes=ceil(MaxEdge);
EdgeStep=MesStep*((N^2-N)/2);MesStep=round(EdgeStep);


%% Group1 vs Group2
ZZ=Z;
Sz=size(ZZ,1);

for i=1:size(ZZ,3)
    z=ZZ(:,:,i);
    z(1:Sz+1:end)=0;z(z<0)=0;
    ZZ(:,:,i)=z;
end

ZZ=ZZ(:,:,[1:n1 n1+1:n1+n2]);

Step = max(max(max(ZZ)))./10000;
fprintf('%-4s\n',' calculating threshold level. it might take tens of minutes depending on the number of networks....');
Thresh_Density=zeros(size([Step:Step:max(max(max(ZZ)))],2),2,size(ZZ,3));

for n=1:size(ZZ,3)
    C=ZZ(:,:,n);
    cn=0;
    for i=Step:Step:max(max(C))
        cn=cn+1;
        C(C<=i)=0;
        Thresh_Density(cn,:,n)=[i nnz(triu(C))];
    end
end

count=0;

for u=MinMes:MesStep:MaxMes
    count=count+1;
    fprintf('%-4s\n',[' calculating network measures at each density: ' num2str(count) ' out of ' num2str(ceil((MaxMes-MinMes)./MesStep)) '...']);
    
    for n=1:size(Thresh_Density,3)
        
        T_Density=Thresh_Density(:,:,n);
        Ind=find(T_Density(:,2) >= u-Precision & T_Density(:,2) <= u+Precision );
        Dind=T_Density(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
            fprintf('%-4s\n',[' No same density within the default precision: decreasing the precision for subject#' num2str(n)]);
            Ind=find(T_Density(:,2) >= u-100*Precision & T_Density(:,2) <= u+100*Precision );
            Dind=Thresh_Density(Ind,2)-u;
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
        end
        
        Thresh=T_Density(Ind,1);
        R=ZZ(:,:,n);R(R<=Thresh)=0;
        NetMes = NetMeasures_Binary(R,NullType,Null_Iter);
        NetMes{23,1}=Thresh;
        NetMes{24,1}=R;
        NetMes_Bin{n,count}=NetMes;
    end
end
MinThr = MinMesPlot;MinIdx = 1;
MaxThr = MaxMesPlot;MaxIdx =  size(NetMes_Bin,2);
NetMes_Bin=NetMes_Bin(1:n1+n2,MinIdx:MaxIdx);

Spare=NetMes_Bin;
NetMes_Bin=Spare(1:n1,:);save(['NetMesBin_G1_vs_G2_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare(n1+1:n1+n2,:);save(['NetMesBin_G1_vs_G2_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare;

%% Group2 vs Group3

ZZ=Z;
Sz=size(ZZ,1);

for i=1:size(ZZ,3)
    z=ZZ(:,:,i);
    z(1:Sz+1:end)=0;z(z<0)=0;
    ZZ(:,:,i)=z;
end

ZZ=ZZ(:,:,[n1+1:n1+n2 n1+n2+1:n1+n2+n3]);
Step = max(max(max(ZZ)))./10000;
fprintf('%-4s\n',' calculating threshold level. it might take tens of minutes depending on the number of networks....');
Thresh_Density=zeros(size([Step:Step:max(max(max(ZZ)))],2),2,size(ZZ,3));

for n=1:size(ZZ,3)
    C=ZZ(:,:,n);
    cn=0;
    for i=Step:Step:max(max(C))
        cn=cn+1;
        C(C<=i)=0;
        Thresh_Density(cn,:,n)=[i nnz(triu(C))];
    end
end

count=0;

for u=MinMes:MesStep:MaxMes
    count=count+1;
    fprintf('%-4s\n',[' calculating network measures at each density: ' num2str(count) ' out of ' num2str(ceil((MaxMes-MinMes)./MesStep)) '...']);
    
    for n=1:size(Thresh_Density,3)
        
        T_Density=Thresh_Density(:,:,n);
        Ind=find(T_Density(:,2) >= u-Precision & T_Density(:,2) <= u+Precision );
        Dind=T_Density(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
            fprintf('%-4s\n',[' No same density within the default precision: decreasing the precision for subject#' num2str(n)]);
            Ind=find(T_Density(:,2) >= u-100*Precision & T_Density(:,2) <= u+100*Precision );
            Dind=Thresh_Density(Ind,2)-u;
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
        end
        
        Thresh=T_Density(Ind,1);
        R=ZZ(:,:,n);R(R<=Thresh)=0;
        NetMes = NetMeasures_Binary(R,NullType,Null_Iter);
        NetMes{23,1}=Thresh;
        NetMes{24,1}=R;
        NetMes_Bin{n,count}=NetMes;
    end
    
end

MinThr = MinMesPlot;MinIdx = 1;
MaxThr = MaxMesPlot;MaxIdx =  size(NetMes_Bin,2);
NetMes_Bin=NetMes_Bin(1:n2+n3,MinIdx:MaxIdx);
Spare=NetMes_Bin;

NetMes_Bin=Spare(1:n2,:);save(['NetMesBin_G2_vs_G3_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin = Spare(n2+1:n2+n3,:);save(['NetMesBin_G2_vs_G3_' Group3 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare;


%% Group1 vs Group3
ZZ=Z;
Sz=size(ZZ,1);

for i=1:size(ZZ,3)
    z=ZZ(:,:,i);
    z(1:Sz+1:end)=0;z(z<0)=0;
    ZZ(:,:,i)=z;
end

ZZ=ZZ(:,:,[1:n1 n1+n2+1:n1+n2+n3]); 
Step = max(max(max(ZZ)))./10000;
fprintf('%-4s\n',' calculating threshold level. it might take tens of minutes depending on the number of networks....');
Thresh_Density=zeros(size([Step:Step:max(max(max(ZZ)))],2),2,size(ZZ,3));

for n=1:size(ZZ,3)
    C=ZZ(:,:,n);
    cn=0;
    for i=Step:Step:max(max(C))
        cn=cn+1;
        C(C<=i)=0;
        Thresh_Density(cn,:,n)=[i nnz(triu(C))];
    end
end

count=0;

for u=MinMes:MesStep:MaxMes
    
    count=count+1;
    fprintf('%-4s\n',[' calculating network measures at each density: ' num2str(count) ' out of ' num2str(ceil((MaxMes-MinMes)./MesStep)) '...']);
    
    for n=1:size(Thresh_Density,3)
        T_Density=Thresh_Density(:,:,n);
        Ind=find(T_Density(:,2) >= u-Precision & T_Density(:,2) <= u+Precision );
        Dind=T_Density(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
            fprintf('%-4s\n',[' No same density within the default precision: decreasing the precision for subject#' num2str(n)]);
            Ind=find(T_Density(:,2) >= u-100*Precision & T_Density(:,2) <= u+100*Precision );
            Dind=Thresh_Density(Ind,2)-u;
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
        end
        
        Thresh=T_Density(Ind,1);
        R=ZZ(:,:,n);R(R<=Thresh)=0;
        NetMes = NetMeasures_Binary(R,NullType,Null_Iter);
        NetMes{23,1}=Thresh;
        NetMes{24,1}=R;
        NetMes_Bin{n,count}=NetMes;
    end
end

MinThr = MinMesPlot;MinIdx = 1;
MaxThr = MaxMesPlot;MaxIdx =  size(NetMes_Bin,2);
NetMes_Bin=NetMes_Bin(1:n1+n3,MinIdx:MaxIdx);

Spare=NetMes_Bin;
NetMes_Bin=Spare(1:n1,:);save(['NetMesBin_G1_vs_G3_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin = Spare(n1+1:n1+n3,:);save(['NetMesBin_G1_vs_G3_' Group3 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare;

fprintf('%-4s\n',' ..........Done ');

