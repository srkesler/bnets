function plotHeatMap(labs,R,outname)
% plot heatmaps

%Edit here-------------------------------------------
% labs = names; %element labels
% R = trueDiffFDR; %what to plot
ftsz = 14; %font size
colmap = 'hot'; % gray, hot, jet, cool, winter, etc.
%----------------------------------------------------
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
figure; imagesc(R);
set(gca, 'XTick', 1:length(labs), 'XTickLabel', labs, 'FontSize',ftsz, 'TickDir', 'out'); 
xticklabel_rotate([],90,[],'interpreter','none');
set(gca, 'YTick', 1:length(labs), 'YTickLabel', labs, 'FontSize',ftsz, 'TickDir', 'out');
colorbar()
colormap(eval(colmap))
Pix_SS = get(0,'screensize');
SS=min(Pix_SS(3),Pix_SS(4));
figure(1, 'position',[1,1,SS,SS]);
saveas(gcf,sprintf('%s.png',outname))
saveas(gcf,sprintf('%s.fig',outname))
close all
end