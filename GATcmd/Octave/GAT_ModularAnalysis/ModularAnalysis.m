function  ModularAnalysis(root, roifilename,mat4GATfilename,sparse,ngroups,nperm,Precision)
%% Inputs: root: Folder containing connectivity data for which you want to run modular analysis
%%       : roifilename: this is the name of the mat file output from conn analysis (ex: resultsROI_Condition.mat)
%%       : mat4GATfilename: this is the matgat file that contains details about groups and number of subjects, saved by gatcmd
%%       : sparse % Whether data is sparse or dense. Use 0 for dense (e.g. fMRI data) and 1 for sparse data (e.g. DTI data)
%%       : ngroups % Number of groups. Options are 2, 3 or 4
%%       : nperm % Number of random permutations (typically 5000)
%%       : Precision;% Adjusts density threshold for sparse data; increase for increasingly sparse data (typically 10)
%% This function does the following things:
% Measure properties (Eglob, Eloc, trans, path length,etc.) for each module
% Next, obtain the above measures for each module for each group as the AUC across densities from minthr:stepthr:maxthr
% Measure ROI specific properties with addition of: Within module degree z score: module_degree_zscore.m
% Next Find network hubs
% Once those are calculated, you determine which ROIs have one or more of those values that are 1 standard deviation or more above the group network mean for that metric
% (i.e. the mean degree for a specific ROI across densities, across subjects is 1 SD or more above the mean degree across densities, across subjects for all ROIs).
% Classify hub type according to the following rule:
% If the mean (not AUC) participation coefficient of the hub across densities is greater than 0.3,  it is a connector hub.  If it is less than 0.3, it is a provincial hub
% The script will output hubsummary.csv
% Note: Network Modules must be run to determine the module assignment (with minimum connection density) for each ROI for each group
% (use regular modularity) before running this script

%% Load required files
load(fullfile(root,roifilename))
ff=load(fullfile(root,mat4GATfilename));

%% Load Modular Community for each group
load(fullfile(root,sprintf('ModularCommunity_%s.mat',ff.mat4GAT.g1)))
load(fullfile(root,sprintf('ModularCommunity_%s.mat',ff.mat4GAT.g2)))
if ngroups>2
    load(fullfile(root,sprintf('ModularCommunity_%s.mat',ff.mat4GAT.g3)))
end
if ngroups>3
    load(fullfile(root,sprintf('ModularCommunity_%s.mat',ff.mat4GAT.g4)))
end

%% Make modularanalysis folder for each group
mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g1))
mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g2))
if ngroups>2
    mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g3))
end
if ngroups>3
    mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g4))
end

%% Load original Z from conn which has data for all three groups
Zorig=Z;
namesorig=names;


%% Group 1 Analysis

%% Group 1  analysis with all rois included regardless of the module membership
cd(fullfile(root,'modularanalysis',ff.mat4GAT.g1))
Z=Zorig(:,:,1:ff.mat4GAT.n1);
names=namesorig;
connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('%s_%s_allrois.mat',roifilename(1:end-4),ff.mat4GAT.g1));
save(connfile,'Z','names','-v7')
if sparse==0
    NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g1,'none',ff.mat4GAT.n1,0);
elseif sparse==1
    NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g1,'none',ff.mat4GAT.n1,0);
end
data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g1,['mat4GAT_' ff.mat4GAT.g1 '_none.mat']);
Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g1,['NetMesBin_' ff.mat4GAT.g1 '_FinalThrRange.mat']);
Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g1,['NetMesBin_none_FinalThrRange.mat']);
Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g1 '.mat']);
Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional',['AUC_NetMesReg_none.mat']);
AUC_Net_Hubs(data_log,Data1,Data2,1)


for i = 2:length(Com_g1) % Loop across each module
    mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i))) % make directory for each module
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i)))% and cd into it
    %% From the original Z from conn for all three groups, extract only subjects for the group and also the rois in a specific module and save it
    [~,~,ib]=intersect(Com_g1{i},namesorig);
    names=namesorig(ib);
    Z=Zorig(ib,ib,1:ff.mat4GAT.n1);
    connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),sprintf('%s_%s_m%d.mat',roifilename(1:end-4),ff.mat4GAT.g1,i));
    save(connfile,'Z','names', '-v7')
    %% Run network measures with all subjects in 1st group and 0/none in the 2nd group
    if sparse==0
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g1,'none',ff.mat4GAT.n1,0);
    elseif sparse==1
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g1,'none',ff.mat4GAT.n1,0);
    end
    
    %% Run Compare Networks (Across Densities) for the two groups (group1=allsubjects vs group2=none)
    % added module degree zscore and Participation coefficient as measures
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i)))
    data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),['mat4GAT_' ff.mat4GAT.g1 '_none.mat']);
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),['NetMesBin_' ff.mat4GAT.g1 '_FinalThrRange.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),['NetMesBin_none_FinalThrRange.mat']);
    
    Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,nperm);
    close all;
    
    %% Run Compare Networks (AUC) for both groups
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i)))
    NetM=fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),'NetMesPerDens.mat');
    
    AUC_Analysis_NetMeasures(data_log,NetM,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr);
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i)))
    
    %% Run Regional Measures (AUC) with nperm repetitions
    % have saved both the normalized measures and non normalized measures.
    % We will be using non normalized measures for the final summary
    
    Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm);
    
    % % % % %     %% Run Networks Hubs (AUC) with 1 standard deviation
    % % % % %     Data1=fullfile(root,'modularanalysis',mat4GAT.g1,sprintf('module%d',i),'Regional',['AUC_NetMesReg_' mat4GAT.g1 '.mat']);
    % % % % %     Data2=fullfile(root,'modularanalysis',mat4GAT.g1,sprintf('module%d',i),'Regional',['AUC_NetMesReg_none.mat']);
    % % % % %     AUC_Net_Hubs(data_log,Data1,Data2,1)
    % % % % %
    %% Creating CSV
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),'Regional',['NetMesReg_' ff.mat4GAT.g1 '.mat']),'MPartCoeff_1');
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional','Hubs_AUC', ['Net_Hubs_Deg_' ff.mat4GAT.g1 '_1SD.mat']),'HubInd1');
    [deghub_all,~, deghub] = intersect(HubInd1,ib);%indices for degree that are more than 1std deviation more than mean
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional','Hubs_AUC', ['Net_Hubs_Betw_' ff.mat4GAT.g1 '_1SD.mat']),'HubInd1');
    [betwhub_all,~, betwhub] = intersect(HubInd1,ib); %indices for node betweenness that are more than 1std deviation more than mean
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional','Hubs_AUC', ['Net_Hubs_Clust_' ff.mat4GAT.g1 '_1SD.mat']),'HubInd1');
    [clusthub_all,~, clusthub] = intersect(HubInd1,ib);  %indices for clustering coefficient that are more than 1std deviation more than mean
    hubs=unique([deghub'; betwhub'; clusthub']); %% Find all hubs that are common (i.e. or-ing) between the three measures
    hubs_all=unique([deghub_all'; betwhub_all'; clusthub_all']);
    [~,~,ideghubs]=intersect(deghub,hubs);
    [~,~,ibetwhubs]=intersect(betwhub,hubs);
    [~,~,iclusthubs]=intersect(clusthub,hubs);
    gg=load(data_log);
    ROI = gg.mat4GAT.roi1;
    results={};
    results(1,1:6)={'Region' 'Degree' 'Betweenness' 'Clustering' 'Participation Coefficient' 'Hub Type'}; % headers for csv
    results(2:length(hubs)+1,1)=ROI(hubs);
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g1,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g1 '.mat']));
    results(2:length(hubs)+1,2)=num2cell(AUC_MDeg_1(hubs_all)'); %Mean AUC of Degree hubs
    for i1=1:length(ideghubs) results{ideghubs(i1)+1,2}=[num2str(results{ideghubs(i1)+1,2}) '*'];end %degree measures that are significant and add a * to those instances
    results(2:length(hubs)+1,3)=num2cell(AUC_MNodeBetw_1(hubs_all)');
    for i2=1:length(ibetwhubs) results{ibetwhubs(i2)+1,3}=[num2str(results{ibetwhubs(i2)+1,3}) '*'];end
    
    results(2:length(hubs)+1,4)=num2cell(AUC_MClust_1(hubs_all)');
    for i3=1:length(iclusthubs) results{iclusthubs(i3)+1,4}=[num2str(results{iclusthubs(i3)+1,4}) '*'];end
    % Find the mean (not auc) participation coefficient of the hub across
    % densities and if it is greater than 0.3,  it is a connector hub
    % otherwise provincial
    results(2:length(hubs)+1,5)=num2cell(MPartCoeff_1(hubs)');
    for i4=1:length(hubs)
        if results{i4+1,5}>=0.3 results{i4+1,6}='Connector'; else results{i4+1,6}='Provincial';end
    end
    
    cell2csv(fullfile(root,'modularanalysis',ff.mat4GAT.g1,sprintf('module%d',i),'Regional','HubSummary.csv'),results)
end
%%Start editing from here
%% Group 2 Analysis



%% Group 2  analysis with all rois included regardless of the module membership
cd(fullfile(root,'modularanalysis',ff.mat4GAT.g2))
Z=Zorig(:,:,ff.mat4GAT.n1+1:ff.mat4GAT.n1+ff.mat4GAT.n2);
connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('%s_%s_allrois.mat',roifilename(1:end-4),ff.mat4GAT.g2));
names=namesorig;
save(connfile,'Z','names','-v7')
if sparse==0
    NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g2,'none',ff.mat4GAT.n2,0);
elseif sparse==1
    NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g2,'none',ff.mat4GAT.n2,0);
end
data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g2,['mat4GAT_' ff.mat4GAT.g2 '_none.mat']);
Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g2,['NetMesBin_' ff.mat4GAT.g2 '_FinalThrRange.mat']);
Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g2,['NetMesBin_none_FinalThrRange.mat']);
Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g2 '.mat']);
Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional',['AUC_NetMesReg_none.mat']);
AUC_Net_Hubs(data_log,Data1,Data2,1)


for i = 1:length(Com_g2)
    mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i)))
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i)))
    %% From the original Z from conn for all three groups, extract only subjects for the group and also the rois in a specific module and save it
    [~,~,ib]=intersect(Com_g2{i},namesorig);
    names=namesorig(ib);
    Z=Zorig(ib,ib,ff.mat4GAT.n1+1:ff.mat4GAT.n1+ff.mat4GAT.n2);
    connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),sprintf('%s_%s_m%d.mat',roifilename(1:end-4),ff.mat4GAT.g2,i));
    save(connfile,'Z','names','-v7')
    
    %% Run network measures with all subjects in 1st group and 0/none in the 2nd group
    if sparse==0
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g2,'none',ff.mat4GAT.n2,0);
    elseif sparse==1
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g2,'none',ff.mat4GAT.n2,0);
    end
    
    %% Run Compare Networks (Across Densities) for the two groups (group1=allsubjects vs group2=none)
    % added module degree zscore and Participation coefficient as measures
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i)))
    data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),['mat4GAT_' ff.mat4GAT.g2 '_none.mat']);
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),['NetMesBin_' ff.mat4GAT.g2 '_FinalThrRange.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),['NetMesBin_none_FinalThrRange.mat']);
    
    Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,nperm)
    close all;
    
    %% Run Compare Networks (AUC) for both groups
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i)))
    NetM=fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),'NetMesPerDens.mat');
    
    AUC_Analysis_NetMeasures(data_log,NetM,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr);
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i)))
    %% Run Regional Measures (AUC) with nperm repetitions
    % have saved both the normalized measures and non normalized measures.
    % We will be using non normalized measures for the final summary
    
    Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm);
    % % % %     %% Run Networks Hubs (AUC) with 1 standard deviation
    % % % %     Data1=fullfile(root,'modularanalysis',mat4GAT.g2,sprintf('module%d',i),'Regional',['AUC_NetMesReg_' mat4GAT.g2 '.mat']);
    % % % %     Data2=fullfile(root,'modularanalysis',mat4GAT.g2,sprintf('module%d',i),'Regional',['AUC_NetMesReg_none.mat']);
    % % % %     AUC_Net_Hubs(data_log,Data1,Data2,1)
    
    %% Creating CSV
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),'Regional',['NetMesReg_' ff.mat4GAT.g2 '.mat']),'MPartCoeff_1');
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional','Hubs_AUC', ['Net_Hubs_Deg_' ff.mat4GAT.g2 '_1SD.mat']),'HubInd1');
    [deghub_all,~, deghub] = intersect(HubInd1,ib);%indices for degree that are more than 1std deviation more than mean
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional','Hubs_AUC', ['Net_Hubs_Betw_' ff.mat4GAT.g2 '_1SD.mat']),'HubInd1');
    [betwhub_all,~, betwhub] = intersect(HubInd1,ib); %indices for node betweenness that are more than 1std deviation more than mean
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional','Hubs_AUC', ['Net_Hubs_Clust_' ff.mat4GAT.g2 '_1SD.mat']),'HubInd1');
    [clusthub_all,~, clusthub] = intersect(HubInd1,ib);  %indices for clustering coefficient that are more than 1std deviation more than mean
    hubs=unique([deghub'; betwhub'; clusthub']); %% Find all hubs that are common (i.e. or-ing) between the three measures
    hubs_all=unique([deghub_all'; betwhub_all'; clusthub_all']);
    [~,~,ideghubs]=intersect(deghub,hubs);
    [~,~,ibetwhubs]=intersect(betwhub,hubs);
    [~,~,iclusthubs]=intersect(clusthub,hubs);
    gg=load(data_log);
    ROI = gg.mat4GAT.roi1;
    results={};
    results(1,1:6)={'Region' 'Degree' 'Betweenness' 'Clustering' 'Participation Coefficient' 'Hub Type'}; % headers for csv
    results(2:length(hubs)+1,1)=ROI(hubs);
    load(fullfile(root,'modularanalysis',ff.mat4GAT.g2,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g2 '.mat']));
    results(2:length(hubs)+1,2)=num2cell(AUC_MDeg_1(hubs_all)'); %Mean AUC of Degree hubs
    for i1=1:length(ideghubs) results{ideghubs(i1)+1,2}=[num2str(results{ideghubs(i1)+1,2}) '*'];end %degree measures that are significant and add a * to those instances
    results(2:length(hubs)+1,3)=num2cell(AUC_MNodeBetw_1(hubs_all)');
    for i2=1:length(ibetwhubs) results{ibetwhubs(i2)+1,3}=[num2str(results{ibetwhubs(i2)+1,3}) '*'];end
    
    results(2:length(hubs)+1,4)=num2cell(AUC_MClust_1(hubs_all)');
    for i3=1:length(iclusthubs) results{iclusthubs(i3)+1,4}=[num2str(results{iclusthubs(i3)+1,4}) '*'];end
    % Find the mean (not auc) participation coefficient of the hub across
    % densities and if it is greater than 0.3,  it is a connector hub
    % otherwise provincial
    results(2:length(hubs)+1,5)=num2cell(MPartCoeff_1(hubs)');
    for i4=1:length(hubs)
        if results{i4+1,5}>=0.3 results{i4+1,6}='Connector'; else results{i4+1,6}='Provincial';end
    end
    
    cell2csv(fullfile(root,'modularanalysis',ff.mat4GAT.g2,sprintf('module%d',i),'Regional','HubSummary.csv'),results)
end

%% Group 3 Analysis
if ngroups>2
    
    %% Group 3  analysis with all rois included regardless of the module membership
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g3))
    Z=Zorig(:,:,ff.mat4GAT.n1+ff.mat4GAT.n2+1:ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3);
    connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('%s_%s_allrois.mat',roifilename(1:end-4),ff.mat4GAT.g3));
    names=namesorig;
    save(connfile,'Z','names','-v7')
    if sparse==0
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g3,'none',ff.mat4GAT.n3,0);
    elseif sparse==1
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g3,'none',ff.mat4GAT.n3,0);
    end
    data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g3,['mat4GAT_' ff.mat4GAT.g3 '_none.mat']);
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g3,['NetMesBin_' ff.mat4GAT.g3 '_FinalThrRange.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g3,['NetMesBin_none_FinalThrRange.mat']);
    Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g3 '.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional',['AUC_NetMesReg_none.mat']);
    AUC_Net_Hubs(data_log,Data1,Data2,1)
    
    %% Group 3  modular analysis
    for i = 1:length(Com_g3)
        mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i)))
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i)))
        %% From the original Z from conn for all three groups, extract only subjects for the group and also the rois in a specific module and save it
        [~,~,ib]=intersect(Com_g3{i},namesorig);
        names=namesorig(ib);
        Z=Zorig(ib,ib,ff.mat4GAT.n1+ff.mat4GAT.n2+1:ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3);
        connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),sprintf('%s_%s_m%d.mat',roifilename(1:end-4),ff.mat4GAT.g3,i));
        save(connfile,'Z','names','-v7')
        %% Run network measures with all subjects in 1st group and 0/none in the 2nd group
        if sparse==0
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g3,'none',ff.mat4GAT.n3,0);
        elseif sparse==1
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g3,'none',ff.mat4GAT.n3,0);
        end
        
        %% Run Compare Networks (Across Densities) for the two groups (group1=allsubjects vs group2=none)
        % added module degree zscore and Participation coefficient as measures
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i)))
        data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),['mat4GAT_' ff.mat4GAT.g3 '_none.mat']);
        Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),['NetMesBin_' ff.mat4GAT.g3 '_FinalThrRange.mat']);
        Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),['NetMesBin_none_FinalThrRange.mat']);
        %Rand_Net_Bootstrap_AcrossSbjComp_f(data_log,data1,data2,Nperm)
        Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,nperm)
        close all;
        
        %% Run Compare Networks (AUC) for both groups
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i)))
        NetM=fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),'NetMesPerDens.mat');
        
        AUC_Analysis_NetMeasures(data_log,NetM,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr)
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i)))
        
        %% Run Regional Measures (AUC) with nperm repetitions
        % have saved both the normalized measures and non normalized measures.
        % We will be using non normalized measures for the final summary
        
        Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
        
        % % %         %% Run Networks Hubs (AUC) with 1 standard deviation
        % % %         Data1=fullfile(root,'modularanalysis',mat4GAT.g3,sprintf('module%d',i),'Regional',['AUC_NetMesReg_' mat4GAT.g3 '.mat']);
        % % %         Data2=fullfile(root,'modularanalysis',mat4GAT.g3,sprintf('module%d',i),'Regional',['AUC_NetMesReg_none.mat']);
        % % %         AUC_Net_Hubs(data_log,Data1,Data2,1)
        
        %% Creating CSV
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),'Regional',['NetMesReg_' ff.mat4GAT.g3 '.mat']),'MPartCoeff_1');
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional','Hubs_AUC', ['Net_Hubs_Deg_' ff.mat4GAT.g3 '_1SD.mat']),'HubInd1');
        [deghub_all,~, deghub] = intersect(HubInd1,ib);%indices for degree that are more than 1std deviation more than mean
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional','Hubs_AUC', ['Net_Hubs_Betw_' ff.mat4GAT.g3 '_1SD.mat']),'HubInd1');
        [betwhub_all,~, betwhub] = intersect(HubInd1,ib); %indices for node betweenness that are more than 1std deviation more than mean
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional','Hubs_AUC', ['Net_Hubs_Clust_' ff.mat4GAT.g3 '_1SD.mat']),'HubInd1');
        [clusthub_all,~, clusthub] = intersect(HubInd1,ib);  %indices for clustering coefficient that are more than 1std deviation more than mean
        hubs=unique([deghub'; betwhub'; clusthub']); %% Find all hubs that are common (i.e. or-ing) between the three measures
        hubs_all=unique([deghub_all'; betwhub_all'; clusthub_all']);
        [~,~,ideghubs]=intersect(deghub,hubs);
        [~,~,ibetwhubs]=intersect(betwhub,hubs);
        [~,~,iclusthubs]=intersect(clusthub,hubs);
        gg=load(data_log);
        ROI = gg.mat4GAT.roi1;
        results={};
        results(1,1:6)={'Region' 'Degree' 'Betweenness' 'Clustering' 'Participation Coefficient' 'Hub Type'}; % headers for csv
        results(2:length(hubs)+1,1)=ROI(hubs);
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g3,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g3 '.mat']));
        results(2:length(hubs)+1,2)=num2cell(AUC_MDeg_1(hubs_all)'); %Mean AUC of Degree hubs
        for i1=1:length(ideghubs) results{ideghubs(i1)+1,2}=[num2str(results{ideghubs(i1)+1,2}) '*'];end %degree measures that are significant and add a * to those instances
        results(2:length(hubs)+1,3)=num2cell(AUC_MNodeBetw_1(hubs_all)');
        for i2=1:length(ibetwhubs) results{ibetwhubs(i2)+1,3}=[num2str(results{ibetwhubs(i2)+1,3}) '*'];end
        
        results(2:length(hubs)+1,4)=num2cell(AUC_MClust_1(hubs_all)');
        for i3=1:length(iclusthubs) results{iclusthubs(i3)+1,4}=[num2str(results{iclusthubs(i3)+1,4}) '*'];end
        % Find the mean (not auc) participation coefficient of the hub across
        % densities and if it is greater than 0.3,  it is a connector hub
        % otherwise provincial
        results(2:length(hubs)+1,5)=num2cell(MPartCoeff_1(hubs)');
        for i4=1:length(hubs)
            if results{i4+1,5}>=0.3 results{i4+1,6}='Connector'; else results{i4+1,6}='Provincial';end
        end
        
        cell2csv(fullfile(root,'modularanalysis',ff.mat4GAT.g3,sprintf('module%d',i),'Regional','HubSummary.csv'),results)
    end
end


%% Group 4 Analysis
if ngroups>3
    
    %% Group 4  analysis with all rois included regardless of the module membership
    cd(fullfile(root,'modularanalysis',ff.mat4GAT.g4))
    Z=Zorig(:,:,ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3+1:ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3+ff.mat4GAT.n4);
    connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('%s_%s_allrois.mat',roifilename(1:end-4),ff.mat4GAT.g4));
    names=namesorig;
    save(connfile,'Z','names','-v7')
    if sparse==0
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g4,'none',ff.mat4GAT.n4,0);
    elseif sparse==1
        NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g4,'none',ff.mat4GAT.n4,0);
    end
    data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g4,['mat4GAT_' ff.mat4GAT.g4 '_none.mat']);
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g4,['NetMesBin_' ff.mat4GAT.g4 '_FinalThrRange.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g4,['NetMesBin_none_FinalThrRange.mat']);
    Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
    Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g4 '.mat']);
    Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional',['AUC_NetMesReg_none.mat']);
    AUC_Net_Hubs(data_log,Data1,Data2,1)
    
    for i = 1:length(Com_g4)
        mkdir(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i)))
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i)))
        %% From the original Z from conn for all three groups, extract only subjects for the group and also the rois in a specific module and save it
        [~,~,ib]=intersect(Com_g4{i},namesorig);
        names=namesorig(ib);
        Z=Zorig(ib,ib,ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3+1:ff.mat4GAT.n1+ff.mat4GAT.n2+ff.mat4GAT.n3+ff.mat4GAT.n4);
        connfile=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),sprintf('%s_%s_m%d.mat',roifilename(1:end-4),ff.mat4GAT.g4,i));
        save(connfile,'Z','names','-v7')
        %% Run network measures with all subjects in 1st group and 0/none in the 2nd group
        if sparse==0
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(0,8,100, 2, 0.05,connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g4,'none',ff.mat4GAT.n4,0);
        elseif sparse==1
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, connfile,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr,ff.mat4GAT.MesStep,ff.mat4GAT.g4,'none',ff.mat4GAT.n4,0);
        end
        
        %% Run Compare Networks (Across Densities) for the two groups (group1=allsubjects vs group2=none)
        % added module degree zscore and Participation coefficient as measures
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i)))
        data_log=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),['mat4GAT_' ff.mat4GAT.g4 '_none.mat']);
        Data1=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),['NetMesBin_' ff.mat4GAT.g4 '_FinalThrRange.mat']);
        Data2=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),['NetMesBin_none_FinalThrRange.mat']);
        
        Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,nperm)
        close all;
        
        %% Run Compare Networks (AUC) for both groups
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i)))
        NetM=fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),'NetMesPerDens.mat');
        
        AUC_Analysis_NetMeasures(data_log,NetM,ff.mat4GAT.MinThr,ff.mat4GAT.MaxThr)
        cd(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i)))
        
        %% Run Regional Measures (AUC) with nperm repetitions
        % have saved both the normalized measures and non normalized measures.
        % We will be using non normalized measures for the final summary
        
        Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm)
        
        % % % %         %% Run Networks Hubs (AUC) with 1 standard deviation
        % % % %
        % % % %         Data1=fullfile(root,'modularanalysis',mat4GAT.g4,sprintf('module%d',i),'Regional',['AUC_NetMesReg_' mat4GAT.g4 '.mat']);
        % % % %         Data2=fullfile(root,'modularanalysis',mat4GAT.g4,sprintf('module%d',i),'Regional',['AUC_NetMesReg_none.mat']);
        % % % %         AUC_Net_Hubs(data_log,Data1,Data2,1)
        
        %% Creating CSV
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),'Regional',['NetMesReg_' ff.mat4GAT.g4 '.mat']),'MPartCoeff_1');
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional','Hubs_AUC', ['Net_Hubs_Deg_' ff.mat4GAT.g4 '_1SD.mat']),'HubInd1');
        [deghub_all,~, deghub] = intersect(HubInd1,ib);%indices for degree that are more than 1std deviation more than mean
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional','Hubs_AUC', ['Net_Hubs_Betw_' ff.mat4GAT.g4 '_1SD.mat']),'HubInd1');
        [betwhub_all,~, betwhub] = intersect(HubInd1,ib); %indices for node betweenness that are more than 1std deviation more than mean
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional','Hubs_AUC', ['Net_Hubs_Clust_' ff.mat4GAT.g4 '_1SD.mat']),'HubInd1');
        [clusthub_all,~, clusthub] = intersect(HubInd1,ib);  %indices for clustering coefficient that are more than 1std deviation more than mean
        hubs=unique([deghub'; betwhub'; clusthub']); %% Find all hubs that are common (i.e. or-ing) between the three measures
        hubs_all=unique([deghub_all'; betwhub_all'; clusthub_all']);
        [~,~,ideghubs]=intersect(deghub,hubs);
        [~,~,ibetwhubs]=intersect(betwhub,hubs);
        [~,~,iclusthubs]=intersect(clusthub,hubs);
        gg=load(data_log);
        ROI = gg.mat4GAT.roi1;
        results={};
        results(1,1:6)={'Region' 'Degree' 'Betweenness' 'Clustering' 'Participation Coefficient' 'Hub Type'}; % headers for csv
        results(2:length(hubs)+1,1)=ROI(hubs);
        load(fullfile(root,'modularanalysis',ff.mat4GAT.g4,'Regional',['AUC_NetMesReg_' ff.mat4GAT.g4 '.mat']));
        results(2:length(hubs)+1,2)=num2cell(AUC_MDeg_1(hubs_all)'); %Mean AUC of Degree hubs
        for i1=1:length(ideghubs) results{ideghubs(i1)+1,2}=[num2str(results{ideghubs(i1)+1,2}) '*'];end %degree measures that are significant and add a * to those instances
        results(2:length(hubs)+1,3)=num2cell(AUC_MNodeBetw_1(hubs_all)');
        for i2=1:length(ibetwhubs) results{ibetwhubs(i2)+1,3}=[num2str(results{ibetwhubs(i2)+1,3}) '*'];end
        
        results(2:length(hubs)+1,4)=num2cell(AUC_MClust_1(hubs_all)');
        for i3=1:length(iclusthubs) results{iclusthubs(i3)+1,4}=[num2str(results{iclusthubs(i3)+1,4}) '*'];end
        % Find the mean (not auc) participation coefficient of the hub across
        % densities and if it is greater than 0.3,  it is a connector hub
        % otherwise provincial
        results(2:length(hubs)+1,5)=num2cell(MPartCoeff_1(hubs)');
        for i4=1:length(hubs)
            if results{i4+1,5}>=0.3 results{i4+1,6}='Connector'; else results{i4+1,6}='Provincial';end
        end
        
        cell2csv(fullfile(root,'modularanalysis',ff.mat4GAT.g4,sprintf('module%d',i),'Regional','HubSummary.csv'),results)
    end
end