% Shelli Kesler
% NOTE: adjust y axis if necessary using e.g. ylim([.9 1.5])

%Just have to make sure when you execute this function you are in the main
%folder of the analysis and not in a particular AUC folder, for example if
%you are in this folder it will work perfectly: /Volumes/wdext/scripts/GAT_2015b_3g/exampledata/Finaltesting_diffusion/3grps_using3grps_iter1

%clear all;


function aucFigs(Group1,Group2,Group3,Group4,minD,dens1,densInt,dens2)

cd ('AUC_Results/AUC_Results_G1_vs_G2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Sigma_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Sigma_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Gamma_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Gamma_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Lambda_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Lambda_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Mod_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Mod_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Assort_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Assort_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'MLocEff_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'MLocEff_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'GEff_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'GEff_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Trans_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Trans_2')
load('NetMesPerDens_AveAcrossDens_D.mat', 'PathL_1')
load('NetMesPerDens_AveAcrossDens_D.mat', 'PathL_2')

cd ../..
cd ('AUC_Results/AUC_Results_G2_vs_G3')

load('NetMesPerDens_AveAcrossDens_D.mat', 'Sigma_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Gamma_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Lambda_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Mod_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Assort_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'MLocEff_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'GEff_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Trans_3')
load('NetMesPerDens_AveAcrossDens_D.mat', 'PathL_3')


cd ../..
cd ('AUC_Results/AUC_Results_G2_vs_G4')

load('NetMesPerDens_AveAcrossDens_D.mat', 'Sigma_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Gamma_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Lambda_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Mod_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Assort_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'MLocEff_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'GEff_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'Trans_4')
load('NetMesPerDens_AveAcrossDens_D.mat', 'PathL_4')

cd ..

% Xvar = .01:.01:.6;
% if length(Xvar)~= size(Sigma_1,2);
%     Xvar = Xvar(1:size(Sigma_1,2));
% end

Xvar = dens1:densInt:dens2;

% Xvar = min:interval:max;
%%
figure
errorbar(Xvar,mean(Sigma_1,1),std(Sigma_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Sigma_2,1),std(Sigma_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Sigma_3,1),std(Sigma_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Sigma_4,1),std(Sigma_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Small-worldness Index  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figSWI.fig')

figure
errorbar(Xvar,mean(Gamma_1,1),std(Gamma_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Gamma_2,1),std(Gamma_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Gamma_3,1),std(Gamma_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Gamma_4,1),std(Gamma_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Normalized Clustering Coefficient  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figGamma.fig')

figure
errorbar(Xvar,mean(Lambda_1,1),std(Lambda_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Lambda_2,1),std(Lambda_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Lambda_3,1),std(Lambda_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Lambda_4,1),std(Lambda_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Characteristic Path Length  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figLambda.fig')

figure
errorbar(Xvar,mean(Mod_1,1),std(Mod_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Mod_2,1),std(Mod_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Mod_3,1),std(Mod_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Mod_4,1),std(Mod_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Modularity  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figMod.fig')

figure
errorbar(Xvar,mean(Assort_1,1),std(Assort_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Assort_2,1),std(Assort_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Assort_3,1),std(Assort_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Assort_4,1),std(Assort_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Assortativity  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figAssort.fig')

figure
errorbar(Xvar,mean(MLocEff_1,1),std(MLocEff_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(MLocEff_2,1),std(MLocEff_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(MLocEff_3,1),std(MLocEff_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(MLocEff_4,1),std(MLocEff_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Mean Local Efficiency  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figMLEFF.fig')

figure
errorbar(Xvar,mean(GEff_1,1),std(GEff_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(GEff_2,1),std(GEff_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(GEff_3,1),std(GEff_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(GEff_4,1),std(GEff_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Global Efficiency  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figGEFF.fig')

figure
errorbar(Xvar,mean(Trans_1,1),std(Trans_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Trans_2,1),std(Trans_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(Trans_3,1),std(Trans_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(Trans_4,1),std(Trans_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Transitivity  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figTrans.fig')

figure
errorbar(Xvar,mean(PathL_1,1),std(PathL_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(PathL_2,1),std(PathL_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
hold on
errorbar(Xvar,mean(PathL_3,1),std(PathL_3,1),'Color',[0 1 0],'LineWidth',3,'DisplayName',Group3)
hold on
errorbar(Xvar,mean(PathL_4,1),std(PathL_4,1),'Color',[0 1 1],'LineWidth',3,'DisplayName',Group4)
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Path Length  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figPL.fig')

%%
cd ..
cd ('AUC_Results/AUC_Results_G1_vs_G2')
load('AUC_NetMesPerDens_AveAcrossDens_D.mat')
G1AUC = cat(2,auc_MClust_1,auc_Trans_1,auc_Assort_1,auc_GEff_1,auc_MLocEff_1,auc_Mod_1,auc_PathL_1,auc_Lambda_1,auc_Gamma_1,auc_Sigma_1);
G2AUC = cat(2,auc_MClust_2,auc_Trans_2,auc_Assort_2,auc_GEff_2,auc_MLocEff_2,auc_Mod_2,auc_PathL_2,auc_Lambda_2,auc_Gamma_2,auc_Sigma_2);
cd ../..
cd ('AUC_Results/AUC_Results_G2_vs_G3')
load('AUC_NetMesPerDens_AveAcrossDens_D.mat')
G3AUC = cat(2,auc_MClust_3,auc_Trans_3,auc_Assort_3,auc_GEff_3,auc_MLocEff_3,auc_Mod_3,auc_PathL_3,auc_Lambda_3,auc_Gamma_3,auc_Sigma_3);

cd ../..
cd ('AUC_Results/AUC_Results_G2_vs_G4')
load('AUC_NetMesPerDens_AveAcrossDens_D.mat')
G4AUC = cat(2,auc_MClust_4,auc_Trans_4,auc_Assort_4,auc_GEff_4,auc_MLocEff_4,auc_Mod_4,auc_PathL_4,auc_Lambda_4,auc_Gamma_4,auc_Sigma_4);

AUCmat = cat(1,G1AUC,G2AUC,G3AUC,G4AUC);
cd ..
save AUCmat AUCmat
csvwrite('AUCmat.csv',AUCmat);

end

