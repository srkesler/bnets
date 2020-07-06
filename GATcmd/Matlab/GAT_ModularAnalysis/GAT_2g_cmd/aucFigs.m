% Shelli Kesler
% NOTE: adjust y axis if necessary using e.g. ylim([.9 1.5])

%clear all;


function aucFigs(Group1,Group2,minD)
load('NetMesPerDens_AveAcrossDens_f.mat', 'Sigma_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Sigma_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Gamma_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Gamma_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Lambda_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Lambda_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Mod_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Mod_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Assort_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Assort_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'MLocEff_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'MLocEff_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'GEff_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'GEff_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Trans_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'Trans_2')
load('NetMesPerDens_AveAcrossDens_f.mat', 'PathL_1')
load('NetMesPerDens_AveAcrossDens_f.mat', 'PathL_2')

% Xvar = .1:.01:.6;
% if length(Xvar)~= size(Sigma_1,2);
%     Xvar = Xvar(1:size(Sigma_1,2));
% end

Xvar = .01:.01:.6;

% Xvar = min:interval:max;

figure
errorbar(Xvar,mean(Sigma_1,1),std(Sigma_1,1),'Color',[0 0.447058826684952 0.74117648601532],'LineWidth',3,'DisplayName',Group1)
hold on
errorbar(Xvar,mean(Sigma_2,1),std(Sigma_2,1),'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],'LineWidth',3,'DisplayName',Group2)
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
hx = graph2d.constantline(minD, 'LineStyle',':', 'Color',[0 0 0],'LineWidth',3,'DisplayName','MinDens');
changedependvar(hx,'x');
title('Path Length  ','FontWeight','bold','FontSize',16);
xlabel('Density','FontWeight','bold','FontSize',16);
legend('show')
set(gca,'fontsize',16,'FontWeight','bold')
savefig('figPL.fig')
end

