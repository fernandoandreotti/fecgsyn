
clearvars -except stat stats_struct
stat = stats_struct;
% 
% for kk=1:length(stat_nopca)
%     mean_FASTICA_DEF(kk) = 100.*mean(stat_nopca{kk}.stats_FASTICA_DEF(:,1));
%     median_FASTICA_DEF(kk) = 100.*median(stat_nopca{kk}.stats_FASTICA_DEF(:,1));
%     mean_FASTICA_SYM(kk) = 100.*mean(stat_nopca{kk}.stats_FASTICA_SYM(:,1));
%     median_FASTICA_SYM(kk) = 100.*median(stat_nopca{kk}.stats_FASTICA_SYM(:,1));
%     mean_JADEICA(kk) = 100.*mean(stat_nopca{kk}.stats_JADEICA(:,1));
%     median_JADEICA(kk) = 100.*median(stat_nopca{kk}.stats_JADEICA(:,1));
%     mean_pca(kk) = 100.*mean(stat_nopca{kk}.stats_pca(:,1));
%     median_pca(kk) = 100.*median(stat_nopca{kk}.stats_pca(:,1));
% end


cc = linspecer(5);
ch = [2,4,6,8,12,16];
FONT_SIZE = 14;
LWIDTH = 2;

% h2 = subplot(1,2,1);
% hold on
% plot(ch,mean_JADEICA,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:))
% plot(ch,median_JADEICA,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:))
% plot(ch,mean_FASTICA_SYM,'d--','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:))
% plot(ch,median_FASTICA_SYM,'d-','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:))
% plot(ch,mean_FASTICA_DEF,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:))
% plot(ch,median_FASTICA_DEF,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:))
% plot(ch,mean_pca,'o--','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:))
% plot(ch,median_pca,'o-','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:))
% xlim([2,32])
% legend('mean ICA (JADE)','median ICA (JADE)','mean ICA (sym. FAST-ICA)','median ICA (sym. FAST-ICA)',...
%     'mean ICA (defl. FAST-ICA)','median ICA (defl. FAST-ICA)','mean PCA','median PCA')
% set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
% set(gca,'FontSize',FONT_SIZE)
% set(gca,'XTick',ch)
% ylabel('F_1 (in %)')
% xlabel('                                                                                                 Number of Channels')
% box on
% title('No dimension reduction')
% hold off


for kk=1:length(stat)
    mean_FASTICA_DEF(kk) = 100.*mean(stat{kk}.stats_FASTICA_DEF(:,1));
    median_FASTICA_DEF(kk) = 100.*median(stat{kk}.stats_FASTICA_DEF(:,1));
    mean_FASTICA_SYM(kk) = 100.*mean(stat{kk}.stats_FASTICA_SYM(:,1));
    median_FASTICA_SYM(kk) = 100.*median(stat{kk}.stats_FASTICA_SYM(:,1));
    mean_JADEICA(kk) = 100.*mean(stat{kk}.stats_JADEICA(:,1));
    median_JADEICA(kk) = 100.*median(stat{kk}.stats_JADEICA(:,1));
    mean_pca(kk) = 100.*mean(stat{kk}.stats_pca(:,1));
    median_pca(kk) = 100.*median(stat{kk}.stats_pca(:,1));
end

for kk=1:length(stat)
    FASTICA_DEF(:,kk) = 100.*stat{kk}.stats_FASTICA_DEF(:,1);
    FASTICA_SYM(:,kk) = 100.*stat{kk}.stats_FASTICA_SYM(:,1);
    JADEICA(:,kk) = 100.*stat{kk}.stats_JADEICA(:,1);
    pca(:,kk) = 100.*stat{kk}.stats_pca(:,1);
end

% for kk=1:length(stat_nopca)
% FASTICA_DEF(1:size(stat_nopca{kk}.stats_FASTICA_DEF,1),kk) = 100.*stat_nopca{kk}.stats_FASTICA_DEF(:,1);
% FASTICA_SYM(1:size(stat_nopca{kk}.stats_FASTICA_SYM,1),kk) = 100.*stat_nopca{kk}.stats_FASTICA_SYM(:,1);
% JADEICA(1:size(stat_nopca{kk}.stats_JADEICA,1),kk) = 100.*stat_nopca{kk}.stats_JADEICA(:,1);
% pca(1:size(stat_nopca{kk}.stats_pca,1),kk) = 100.*stat_nopca{kk}.stats_pca(:,1);
% end


ch = [2,4,6,8,12,16,32];

close
figure(1)
% h1 = subplot(1,2,2);
hold on
x = ch;
h(1)=shadedErrorBar(x, FASTICA_DEF, {@nanmean, @(x) nanstd(x)  },{'d--','markerfacecolor',[0.7 0.7 0.7],'color',[0.7 0.7 0.7],'linewidth',LWIDTH},0);
h(2)=shadedErrorBar(x, FASTICA_SYM, {@nanmean, @(x) nanstd(x)  },{'^--','markerfacecolor',cc(1,:),'color',cc(1,:),'linewidth',LWIDTH},1);
h(3)=shadedErrorBar(x, JADEICA, {@nanmean, @(x) nanstd(x)  },{'s--','markerfacecolor',cc(2,:),'color',cc(2,:),'linewidth',LWIDTH},1);
h(4)=shadedErrorBar(x, pca, {@nanmean, @(x) nanstd(x)  },{'o--','markerfacecolor',cc(4,:),'color',cc(4,:),'linewidth',LWIDTH},1);
han = arrayfun(@(x) x.mainLine,h,'UniformOutput',0);
legend([han{:}],'defl. FAST-ICA','sym. FAST-ICA','JADE','PCA')
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
set(gca,'FontSize',FONT_SIZE)
set(gca,'XTick',ch)
ylabel('F_1 (in %)')
xlabel('Number of Channels')
box on
% title('PCA reduction')
xlim([2 32])
ylim([86 100])
hold off

figure(2)
hold on
plot(ch,mean_JADEICA,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:),'linewidth',LWIDTH)
plot(ch,median_JADEICA,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:),'linewidth',LWIDTH)
plot(ch,mean_FASTICA_SYM,'d--','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:),'linewidth',LWIDTH)
plot(ch,median_FASTICA_SYM,'d-','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:),'linewidth',LWIDTH)
plot(ch,mean_FASTICA_DEF,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:),'linewidth',LWIDTH)
plot(ch,median_FASTICA_DEF,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:),'linewidth',LWIDTH)
plot(ch,mean_pca,'o--','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:),'linewidth',LWIDTH)
plot(ch,median_pca,'o-','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:),'linewidth',LWIDTH)
xlim([2,32])
legend('mean ICA (JADE)','median ICA (JADE)','mean ICA (sym. FAST-ICA)','median ICA (sym. FAST-ICA)',...
    'mean ICA (defl. FAST-ICA)','median ICA (defl. FAST-ICA)','mean PCA','median PCA')
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
set(gca,'FontSize',FONT_SIZE)
set(gca,'XTick',ch)
ylabel('F_1 (in %)')
xlabel('Number of Channels')
box on
% title('PCA reduction')
hold off




