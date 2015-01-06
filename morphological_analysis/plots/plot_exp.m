
load('stats_ch_5')
FONT_SIZE = 14;
LWIDTH = 2;
figure(1)
ch = [2,4,6,8,16];
h1 = subplot(1,2,1);
hold on
plot(ch,mean_JADEICA,'sr--','linewidth',LWIDTH,'MarkerFaceColor','r')
plot(ch,median_JADEICA,'sr-','linewidth',LWIDTH,'MarkerFaceColor','r')
plot(ch,mean_FASTICA_SYM,'dk--','linewidth',LWIDTH,'MarkerFaceColor','k')
plot(ch,median_FASTICA_SYM,'dk-','linewidth',LWIDTH,'MarkerFaceColor','k')
plot(ch,mean_pca,'ob--','linewidth',LWIDTH,'MarkerFaceColor','b')
plot(ch,median_pca,'ob-','linewidth',LWIDTH,'MarkerFaceColor','b')
xlim([2,16])
legend('mean ICA (JADE)','median ICA (JADE)','mean ICA (FAST-ICA)','median ICA (FAST-ICA)',...
    'mean PCA','median PCA')
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
set(gca,'FontSize',FONT_SIZE)
set(gca,'XTick',ch)
ylabel('F_1 (in %)')
xlabel('                                                                                                 Number of Channels')
box on
hold off
h2 = subplot(1,2,2);
hold on
plot(ch,mean_FASTICA_SYM,'dk--','linewidth',LWIDTH,'MarkerFaceColor','k')
plot(ch,median_FASTICA_SYM,'dk-','linewidth',LWIDTH,'MarkerFaceColor','k')
plot(ch,mean_FASTICA_DEF,'sm--','linewidth',LWIDTH,'MarkerFaceColor','m')
plot(ch,median_FASTICA_DEF,'sm-','linewidth',LWIDTH,'MarkerFaceColor','m')

xlim([2,16])

legend('mean (FAST-ICA SYMMETRIC)','median ICA (FAST-ICA SYMMETRIC)',...
    'mean (FAST-ICA DEFLATIONARY)','median (FAST-ICA DEFLATIONARY)')
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
set(gca,'FontSize',FONT_SIZE)
set(gca,'XTick',ch)
box on
hold off

