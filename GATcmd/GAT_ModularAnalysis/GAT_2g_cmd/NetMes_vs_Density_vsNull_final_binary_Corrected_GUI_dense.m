function NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(IfPara, NumCores, Null_Iter, Tail, Alpha,connfile,MinMes,MaxMes,MesStep,Group1,Group2,n1,n2)

% The program accepts the Results_ROI (output from conn) of both groups
% combined (note that the groups should be in order e.g. 1:34 BC; 35:61 CON
%'MinMes': The minimum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.05 for density)
%'MaxMes': The maximum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.6 for density)


%%
%-- Hadi Hosseini March 22, 2011
%-- Hadi Hosseini August 18, 2011 (updated for functional connectivity)
%-- Hadi Hosseini August 18, 2011 (updated for covariates of no interest)
%-- updated on Mar 2014 (parallel processing)
%-- updated on Apr 2020 by Shelli Kesler and Vikram Rao

%Ref (Hosseini et al.,Plos One 2012)


%% Main Body

if isequal(IfPara,1)
    MaxCores = feature('numCores');
    if NumCores > MaxCores
        NumCores = 0;
        error('too many cores!')
    end
else
    NumCores = 0;
end

mat4GAT.para = IfPara;
mat4GAT.ncore = NumCores;

%null_flag = input('type 2 (same deg dist), 3 (H-Q-S), 4 (RandCorr): ');
null_flag = 2;%at this moment, only supports same deg dist
mat4GAT.null = null_flag;
NullType = null_flag;

% % % Null_Iter = input('number of null networks (e.g. 20):  ');
mat4GAT.nullIter = Null_Iter;

% % % Tail = input('one (type 1) or two (type 2) -tailed test:  ');
mat4GAT.tail = Tail;

% % % Alpha = input('significance threshold level (e.g. 0.05):  ');
mat4GAT.alpha = Alpha;

ff = connfile; %%%spm_select(1,'resultsROI_*','Select resultsROI_Condition*.mat file (output from Conn toolbox)');
fprintf('%-4s\n',' loading the input data....');
load(ff,'Z','names');

[nROI,nROI2]=size(Z);

ZZ=tanh(Z);
MaxD = [];

for i=1:size(ZZ,3)
    z=ZZ(:,:,i);
    z(1:nROI+1:end)=0;z(z<0)=0;
    ZZ(:,:,i)=z;
    MaxD(i) = density_und(z);
end

MaxMaxD = max(MaxD);

MinMesPlot = MinMes;
MaxMesPlot = MaxMes;
MesStepPlot = MesStep;

Steps = [MinMes:MesStep:MaxMes];
NStep = length(Steps);
nSbj = size(ZZ,3);
NetMes_Bin = cell(nSbj,NStep);

fprintf('%-4s\n',[' calculating network measures ...']);

if isequal(IfPara,1)
    
    parfor_progress(nSbj);
    parfor n = 1:nSbj
        parfor_progress;
        R = ZZ(:,:,n);
        for count = 1:NStep
            TDens = Steps(count);
            R_ = costThresh({R},TDens);
            R_ = R_{1};
            try
                NetMes = NetMeasures_Binary(R_,NullType,Null_Iter);
            catch
                error('Minimum threshold density too low, try increasing the minimum threshold to a higher value')
            end
            NetMes{24,1} = R_;
            NetMes_Bin{n,count} = NetMes;
        end
    end
    parfor_progress(0);
else
    for n = 1:nSbj
        R = ZZ(:,:,n);
        for count = 1:NStep
            TDens = Steps(count);
            R_ = costThresh({R},TDens);
            R_ = R_{1};
            try
                NetMes = NetMeasures_Binary(R_,NullType,Null_Iter);
            catch
                error('Minimum threshold density too low, try increasing the minimum threshold to a higher value')
            end
            NetMes{24,1} = R_;
            NetMes_Bin{n,count} = NetMes;
        end
    end
end

if ~isequal(n1+n2,nSbj)
    warning('the number of subjects is not consistent with the input data');
end

mat4GAT.g1 = Group1;mat4GAT.g2 = Group2;
mat4GAT.n1 = n1;mat4GAT.n2 = n2;
mat4GAT.MinThr = MinMesPlot;
mat4GAT.MaxThr = MaxMesPlot;
mat4GAT.MesStep = MesStepPlot;
mat4GAT.roi1 = names;
mat4GAT.roi2 = names;
save(['mat4GAT_' Group1 '_' Group2],'mat4GAT');

Spare = NetMes_Bin;
NetMes_Bin = Spare(1:n1,:);save(['NetMesBin_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin = Spare(n1+1:n1+n2,:);save(['NetMesBin_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin = Spare;

fprintf('%-4s\n',' ..........Done ');