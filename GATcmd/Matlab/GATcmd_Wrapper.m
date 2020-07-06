% GAT-CMD relies on Brain Connectivity Toolbox
% https://www.nitrc.org/projects/bct/
% This needs to be in your Matlab path

% GAT-CMD requires Matlab with Parallel Computing, Statistics and Machine Learning,
% and Image Processing Toolboxes

% GAT-CMD is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes.

% The connfile should be a nxnxm matrix of weights representing nxn symmetric
% connectivity matrices for m individuals

% USER EDIT HERE-----------------------------------------------------------------------------
GAT='/Users/Desktop/GATcmd'; %edit path to GAT-CMD folder
BCTpath='/Users/Desktop/GATcmd/2019_03_03_BCT'; %edit path to Brain Connectivity Toolbox folder

root='/Users/Desktop/Study/'; %% Folder containing connectivity data. This is also where the results will be stored
cd(root)
connfile='MyMatrix.mat'; % Name of connectivity file

sparse=0;    % Use 0 for dense (e.g. fMRI data) and 1 for sparse data (e.g. DTI data)
nperm=5000;  % Number of random permutations
Precision=10;% Adjusts density threshold for sparse data; increase for increasingly sparse data

testdisc=1; % Minimum density = minimum connection density; otherwise set to 0 and change values below
minthr=0.01; % Minimum density threshold
maxthr=0.5; % Maximum density threshold
stepthr=0.01; %Density Interval

ngroups=2; % Number of groups. Options are 2, 3 or 4

g1='exp1'; % Name of Group 1
g2='exp2';  % Name of Group 2
g3='exp3'; % Name of Group 3, otherwise, change this to ''
g4='exp4'; % Name of Group 4, otherwise, change this to ''

n1=125; % Number of subjects in Group 1
n2=85; % Number of subjects in Group 2
n3=100; % Number of subjects in Group 3, otherwise, change this to 0
n4=171; % Number of subjects in Group 4, otherwise, change this to 0

modular=0; % Set to 1 if you want modular analysis and hub characterization
regional=0; % Set to 1 if you want regional analysis
long=0; % Set to 1 if you want longitudinal analysis


%% No further edits needed below -------------------------------------------

GAT2=fullfile(GAT,'GAT_2g_cmd');
GAT3=fullfile(GAT,'GAT_3g_cmd');
GAT4=fullfile(GAT,'GAT_4g_cmd');
MODpath=fullfile(GAT,'GAT_ModularAnalysis');
longpath=fullfile(GAT,'Longitudinal');

addpath(genpath(BCTpath))

switch ngroups
    
    case 2
        
        addpath(genpath(GAT2))
        if sparse==0 % Functional or Structural data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(1, feature('numcores'), 100, 2, 0.05,fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,n1,n2);
        elseif sparse==1 % Diffusion Data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05,Precision,fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,n1,n2);
        end
        data_log=fullfile(root,sprintf('mat4GAT_%s_%s.mat',g1,g2));
        load(data_log);
        Data1=sprintf('NetMesBin_%s_FinalThrRange.mat',mat4GAT.g1);
        Data2=sprintf('NetMesBin_%s_FinalThrRange.mat',mat4GAT.g2);
        Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1,Data2,nperm);
        NetM=fullfile(root,'NetMesPerDens.mat');
        
        if regional
            Net_RegDiff_withBootstrap(data_log,Data1,Data2,nperm);
            regMes_2g(data_log,Data1,Data2,[g1 '_' g2]);
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g2 '.mat'] ,[g1 '_' g2], nperm);
            load (['PermStats_' g1 '_' g2 '.mat'], 'trueDiffFDR');
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g2]);
            cd ..
        end
        
        maxdisc_data1 = Test_Discon_NetMesBin(Data1,data_log);
        maxdisc_data2 = Test_Discon_NetMesBin(Data2,data_log);
        disc=max(maxdisc_data1,maxdisc_data2);
        
        if testdisc
            AUC_Analysis_NetMeasures(data_log,NetM,disc,maxthr);
        else
            AUC_Analysis_NetMeasures(data_log,NetM,minthr,maxthr);
        end
        
        if long
            addpath(genpath(longpath))
            cd AUC_Results
            permuteGlobalStats_longbetween_1g (nperm)
            cd ..
            
            regMes_long_1g(data_log,Data1,Data2);
            permuteMatrixStats_sparselongwithin(1,[g1 '_' g2],nperm)
            load (['PermStats_longwithin_' g1 '_' g2 '.mat'], 'trueDiffFDR');
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_long_' g1 '_' g2]);
            rmpath(genpath(longpath))
            cd ..
        end
        
        if modular
            cd (root)
            Module_Analysis_all(data_log,Data1,Data2,100,disc)
            rmpath(genpath(GAT2))
            addpath(genpath(MODpath))
            ModularAnalysis(root,connfile,sprintf('mat4GAT_%s_%s.mat',g1,g2),sparse,ngroups,nperm,Precision)
            rmpath(genpath(MODpath))
        end
        rmpath(genpath(GAT2))
        
    case 3
        
        addpath(genpath(GAT3))
        if sparse==0 % Functional or Structural data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(1, feature('numcores'), 100, 2, 0.05,fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,g3,n1,n2,n3);
        elseif sparse==1 % Diffusion Data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,g3,n1,n2,n3);
        end
        data_log=fullfile(root,sprintf('mat4GAT_%s_%s_%s.mat',g1,g2,g3));
        load(data_log);
        Data1_g1g2=sprintf('NetMesBin_G1_vs_G2_%s_FinalThrRange.mat',mat4GAT.g1);
        Data2_g1g2=sprintf('NetMesBin_G1_vs_G2_%s_FinalThrRange.mat',mat4GAT.g2);
        Data2_g2g3=sprintf('NetMesBin_G2_vs_G3_%s_FinalThrRange.mat',mat4GAT.g2);
        Data3_g2g3=sprintf('NetMesBin_G2_vs_G3_%s_FinalThrRange.mat',mat4GAT.g3);
        Data1_g1g3=sprintf('NetMesBin_G1_vs_G3_%s_FinalThrRange.mat',mat4GAT.g1);
        Data3_g1g3=sprintf('NetMesBin_G1_vs_G3_%s_FinalThrRange.mat',mat4GAT.g3);
        Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1_g1g2,Data2_g1g2,Data2_g2g3,Data3_g2g3, Data1_g1g3, Data3_g1g3,nperm);
        NetM_1_2=fullfile(root,'NetMesPerDens_G1_vs_G2.mat');
        NetM_2_3=fullfile(root,'NetMesPerDens_G2_vs_G3.mat');
        NetM_1_3=fullfile(root,'NetMesPerDens_G1_vs_G3.mat');
        
        if regional
            Net_RegDiff_withBootstrap(data_log,Data1_g1g2,Data2_g1g2,Data2_g2g3,Data3_g2g3, Data1_g1g3, Data3_g1g3,nperm);
            
            regMes_2g(data_log,Data1_g1g2,Data2_g1g2,[g1 '_' g2]);
            movefile(['RegMesAUC_' g1 '_' g2 '.mat'],[root '/' 'Regional' '/' 'Regional_G1_vs_G2' '/' 'RegMesAUC_' g1 '_' g2 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G1_vs_G2'])
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g2 '.mat'] ,[g1 '_' g2], nperm)
            load (['PermStats_' g1 '_' g2 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g2])
            cd ../..
            
            regMes_2g(data_log,Data2_g2g3,Data3_g2g3,[g2 '_' g3]);
            movefile(['RegMesAUC_' g2 '_' g3 '.mat'],[root '/' 'Regional' '/' 'Regional_G2_vs_G3' '/' 'RegMesAUC_' g2 '_' g3 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G2_vs_G3'])
            permuteMatrixStats_sparse(['RegMesAUC_' g2 '_' g3 '.mat'] ,[g2 '_' g3], nperm)
            load (['PermStats_' g2 '_' g3 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g2 '_' g3])
            cd ../..
            
            regMes_2g(data_log,Data1_g1g3,Data3_g1g3,[g1 '_' g3]);
            movefile(['RegMesAUC_' g1 '_' g3 '.mat'],[root '/' 'Regional' '/' 'Regional_G1_vs_G3' '/' 'RegMesAUC_' g1 '_' g3 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G1_vs_G3'])
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g3 '.mat'] ,[g1 '_' g3], nperm)
            load (['PermStats_' g1 '_' g3 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g3])
            cd ../..
            
        end
        
        maxdisc_Data1_g1g2 = Test_Discon_NetMesBin(Data1_g1g2,data_log);
        maxdisc_Data2_g1g2 = Test_Discon_NetMesBin(Data2_g1g2,data_log);
        maxdisc_Data2_g2g3 = Test_Discon_NetMesBin(Data2_g2g3,data_log);
        maxdisc_Data3_g2g3 = Test_Discon_NetMesBin(Data3_g2g3,data_log);
        maxdisc_Data1_g1g3 = Test_Discon_NetMesBin(Data1_g1g3,data_log);
        maxdisc_Data3_g1g3 = Test_Discon_NetMesBin(Data3_g1g3,data_log);
        disc=max([ maxdisc_Data1_g1g2, maxdisc_Data2_g1g2,maxdisc_Data2_g2g3, maxdisc_Data3_g2g3, maxdisc_Data1_g1g3, maxdisc_Data3_g1g3]);
        
        if testdisc
            AUC_Analysis_NetMeasures(data_log,NetM_1_2,NetM_2_3,NetM_1_3,disc,maxthr);
        else
            AUC_Analysis_NetMeasures(data_log,NetM_1_2,NetM_2_3,NetM_1_3,minthr,maxthr);
        end
        
        if modular
            cd (root)
            Module_Analysis_all(data_log,Data1_g1g2,Data2_g1g2,Data3_g2g3, 100,disc)
            rmpath(genpath(GAT3))
            addpath(genpath(MODpath))
            ModularAnalysis(root,connfile,sprintf('mat4GAT_%s_%s_%s.mat',g1,g2,g3),sparse,ngroups,nperm,Precision)
            rmpath(genpath(MODpath))
        end
        rmpath(genpath(GAT3))
        
    case 4
        
        addpath(genpath(GAT4))
        if sparse==0 % Functional or Structural data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_dense(1, feature('numcores'), 100, 2, 0.05,fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,g3,g4,n1,n2,n3,n4);
        elseif sparse==1 % Diffusion Data
            NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_sparse(100, 2, 0.05, Precision, fullfile(root,connfile),minthr,maxthr,stepthr,g1,g2,g3,g4,n1,n2,n3,n4);
        end
        data_log=fullfile(root,sprintf('mat4GAT_%s_%s_%s_%s.mat',g1,g2,g3,g4));
        load(data_log);
        Data1_g1g2=sprintf('NetMesBin_G1_vs_G2_%s_FinalThrRange.mat',mat4GAT.g1);
        Data2_g1g2=sprintf('NetMesBin_G1_vs_G2_%s_FinalThrRange.mat',mat4GAT.g2);
        Data2_g2g3=sprintf('NetMesBin_G2_vs_G3_%s_FinalThrRange.mat',mat4GAT.g2);
        Data3_g2g3=sprintf('NetMesBin_G2_vs_G3_%s_FinalThrRange.mat',mat4GAT.g3);
        Data1_g1g3=sprintf('NetMesBin_G1_vs_G3_%s_FinalThrRange.mat',mat4GAT.g1);
        Data3_g1g3=sprintf('NetMesBin_G1_vs_G3_%s_FinalThrRange.mat',mat4GAT.g3);
        Data1_g1g4=sprintf('NetMesBin_G1_vs_G4_%s_FinalThrRange.mat',mat4GAT.g1);
        Data4_g1g4=sprintf('NetMesBin_G1_vs_G4_%s_FinalThrRange.mat',mat4GAT.g4);
        Data2_g2g4=sprintf('NetMesBin_G2_vs_G4_%s_FinalThrRange.mat',mat4GAT.g2);
        Data4_g2g4=sprintf('NetMesBin_G2_vs_G4_%s_FinalThrRange.mat',mat4GAT.g4);
        Data3_g3g4=sprintf('NetMesBin_G3_vs_G4_%s_FinalThrRange.mat',mat4GAT.g3);
        Data4_g3g4=sprintf('NetMesBin_G3_vs_G4_%s_FinalThrRange.mat',mat4GAT.g4);
        Rand_Net_Bootstrap_AcrossSbjComp(data_log,Data1_g1g2,Data2_g1g2,Data2_g2g3,Data3_g2g3, Data1_g1g3, Data3_g1g3,Data1_g1g4,Data4_g1g4,Data2_g2g4,Data4_g2g4,Data3_g3g4,Data4_g3g4,nperm);
        NetM_1_2=fullfile(root,'NetMesPerDens_G1_vs_G2.mat');
        NetM_2_3=fullfile(root,'NetMesPerDens_G2_vs_G3.mat');
        NetM_1_3=fullfile(root,'NetMesPerDens_G1_vs_G3.mat');
        NetM_1_4=fullfile(root,'NetMesPerDens_G1_vs_G4.mat');
        NetM_2_4=fullfile(root,'NetMesPerDens_G2_vs_G4.mat');
        NetM_3_4=fullfile(root,'NetMesPerDens_G3_vs_G4.mat');
        
        if regional
            Net_RegDiff_withBootstrap(data_log,Data1_g1g2,Data2_g1g2,Data2_g2g3,Data3_g2g3, Data1_g1g3, Data3_g1g3,Data1_g1g4,Data4_g1g4,Data2_g2g4,Data4_g2g4,Data3_g3g4,Data4_g3g4,nperm);
            
            regMes_2g(data_log,Data1_g1g2,Data2_g1g2,[g1 '_' g2]);
            movefile(['RegMesAUC_' g1 '_' g2 '.mat'],[root '/' 'Regional' '/' 'Regional_G1_vs_G2' '/' 'RegMesAUC_' g1 '_' g2 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G1_vs_G2'])
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g2 '.mat'] ,[g1 '_' g2], nperm)
            load (['PermStats_' g1 '_' g2 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g2])
            cd ../..
            
            regMes_2g(data_log,Data2_g2g3,Data3_g2g3,[g2 '_' g3]);
            movefile(['RegMesAUC_' g2 '_' g3 '.mat'],[root '/' 'Regional' '/' 'Regional_G2_vs_G3' '/' 'RegMesAUC_' g2 '_' g3 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G2_vs_G3'])
            permuteMatrixStats_sparse(['RegMesAUC_' g2 '_' g3 '.mat'] ,[g2 '_' g3], nperm)
            load (['PermStats_' g2 '_' g3 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g2 '_' g3])
            cd ../..
            
            regMes_2g(data_log,Data1_g1g3,Data3_g1g3,[g1 '_' g3]);
            movefile(['RegMesAUC_' g1 '_' g3 '.mat'],[root '/' 'Regional' '/' 'Regional_G1_vs_G3' '/' 'RegMesAUC_' g1 '_' g3 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G1_vs_G3'])
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g3 '.mat'] ,[g1 '_' g3], nperm)
            load (['PermStats_' g1 '_' g3 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g3])
            cd ../..
            
            regMes_2g(data_log,Data1_g1g4,Data4_g1g4,[g1 '_' g4]);
            movefile(['RegMesAUC_' g1 '_' g4 '.mat'],[root '/' 'Regional' '/' 'Regional_G1_vs_G4' '/' 'RegMesAUC_' g1 '_' g4 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G1_vs_G4'])
            permuteMatrixStats_sparse(['RegMesAUC_' g1 '_' g4 '.mat'] ,[g1 '_' g4], nperm)
            load (['PermStats_' g1 '_' g4 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g1 '_' g4])
            cd ../..
            
            regMes_2g(data_log,Data2_g2g4,Data4_g2g4,[g2 '_' g4]);
            movefile(['RegMesAUC_' g2 '_' g4 '.mat'],[root '/' 'Regional' '/' 'Regional_G2_vs_G4' '/' 'RegMesAUC_' g2 '_' g4 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G2_vs_G4'])
            permuteMatrixStats_sparse(['RegMesAUC_' g2 '_' g4 '.mat'] ,[g2 '_' g4], nperm)
            load (['PermStats_' g2 '_' g4 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g2 '_' g4])
            cd ../..
            
            regMes_2g(data_log,Data3_g3g4,Data4_g3g4,[g3 '_' g4]);
            movefile(['RegMesAUC_' g3 '_' g4 '.mat'],[root '/' 'Regional' '/' 'Regional_G3_vs_G4' '/' 'RegMesAUC_' g3 '_' g4 '.mat'])
            cd([root '/' 'Regional' '/' 'Regional_G3_vs_G4'])
            permuteMatrixStats_sparse(['RegMesAUC_' g3 '_' g4 '.mat'] ,[g3 '_' g4], nperm)
            load (['PermStats_' g3 '_' g4 '.mat'], 'trueDiffFDR')
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_' g3 '_' g4])
            cd ../..
            
        end
        
        maxdisc_Data1_g1g2 = Test_Discon_NetMesBin(Data1_g1g2,data_log);
        maxdisc_Data2_g1g2 = Test_Discon_NetMesBin(Data2_g1g2,data_log);
        maxdisc_Data2_g2g3 = Test_Discon_NetMesBin(Data2_g2g3,data_log);
        maxdisc_Data3_g2g3 = Test_Discon_NetMesBin(Data3_g2g3,data_log);
        maxdisc_Data1_g1g3 = Test_Discon_NetMesBin(Data1_g1g3,data_log);
        maxdisc_Data3_g1g3 = Test_Discon_NetMesBin(Data3_g1g3,data_log);
        maxdisc_Data1_g1g4 = Test_Discon_NetMesBin(Data1_g1g4,data_log);
        maxdisc_Data4_g1g4 = Test_Discon_NetMesBin(Data4_g1g4,data_log);
        maxdisc_Data2_g2g4 = Test_Discon_NetMesBin(Data2_g2g4,data_log);
        maxdisc_Data4_g2g4 = Test_Discon_NetMesBin(Data4_g2g4,data_log);
        maxdisc_Data3_g3g4 = Test_Discon_NetMesBin(Data3_g3g4,data_log);
        maxdisc_Data4_g3g4 = Test_Discon_NetMesBin(Data4_g3g4,data_log);
        disc=max([ maxdisc_Data1_g1g2, maxdisc_Data2_g1g2,maxdisc_Data2_g2g3, maxdisc_Data3_g2g3, maxdisc_Data1_g1g3, maxdisc_Data3_g1g3, maxdisc_Data1_g1g4, maxdisc_Data4_g1g4, maxdisc_Data2_g2g4, maxdisc_Data4_g2g4, maxdisc_Data3_g3g4, maxdisc_Data4_g3g4]);
        
        if testdisc
            AUC_Analysis_NetMeasures(data_log,NetM_1_2,NetM_2_3,NetM_1_3,NetM_1_4,NetM_2_4,NetM_3_4,disc,maxthr);
        else
            AUC_Analysis_NetMeasures(data_log,NetM_1_2,NetM_2_3,NetM_1_3,NetM_1_4,NetM_2_4,NetM_3_4,minthr,maxthr);
        end
        
        if long
            addpath(genpath(longpath))
            cd AUC_Results
            permuteGlobalStats_longbetween_2g (nperm)
            cd ..
            
            regMes_long_2g(data_log,Data1_g1g2,Data2_g1g2,Data3_g3g4,Data4_g3g4);
            permuteMatrixStats_sparselongbetween (nperm)
            permuteMatrixStats_sparselongwithin(1,[g1 '_' g2], nperm)
            load (['PermStats_longwithin_' g1 '_' g2 '.mat'], 'trueDiffFDR');
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_long_' g1 '_' g2]);
            
            
            permuteMatrixStats_sparselongwithin(2,[g3 '_' g4], nperm)
            load (['PermStats_longwithin_' g3 '_' g4 '.mat'], 'trueDiffFDR');
            plotHeatMap(mat4GAT.roi1,trueDiffFDR,['Heatmap_long_' g3 '_' g4]);
            rmpath(genpath(longpath))
            cd ..
        end
        
        if modular
            cd (root)
            Module_Analysis_all(data_log,Data1_g1g2,Data2_g1g2,Data3_g2g3,Data4_g1g4, 100,disc)
            rmpath(genpath(GAT4))
            addpath(genpath(MODpath))
            ModularAnalysis(root,connfile,sprintf('mat4GAT_%s_%s_%s_%s.mat',g1,g2,g3,g4),sparse,ngroups,nperm, Precision)
            rmpath(genpath(MODpath))
        end
        rmpath(genpath(GAT4))
end

rmpath(genpath(BCTpath))
