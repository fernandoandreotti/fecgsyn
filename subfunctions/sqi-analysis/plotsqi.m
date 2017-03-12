ch = 4;
collist=[31;26;22;9;25;17;14;12;21;32];
collist=[8;12;21;26;31;36;40;42;44;46]-4;
features = feats;



% plot definitions
markers = {'-+','-o','-*','-.','-x','-s','-d','-^','-v','>-','-<','-p','-h'};
markerorder = cell(length(collist),1);
for m = 1:length(collist); markerorder{m} = markers{randi(length(markers),1)};end
markerorder = cellstr(markerorder);

% plotting
% set(gca(), 'LineStyleOrder',markerorder, 'NextPlot','replacechildren')
% plot(Data.preprocessed(ch,:),'-','Color',[0.7 0.7 0.7])
% hold on
% plot(Data.fECGestim(ch,:),'-r')
% plot(Data.mRef.samplestamp,-1.*ones(length(Data.mRef.samplestamp),1),'dg')
% plot(table2array(Data.sqi(Data.sqi.channel==ch,4)),features)
% legend(['rawdata','extracted','mRef',Data.sqi.Properties.VariableNames(collist)])

close all
ax(1)=subplot(3,1,1);
tm = 1:5:Data.length;
prep = resample(Data.preprocessed(ch,:),1,5);
plot(tm,prep,'-','Color',[0.7 0.7 0.7])
hold on
plot(Data.mRef.samplestamp,1.*ones(length(Data.mRef.samplestamp),1),'dg')
legend('preprocessed','MQRS')
ax(2)=subplot(3,1,2);
fecg = resample(Data.fECGestim(ch,:),1,5);
plot(tm,fecg,'-b')
hold on
plot(Data.fRef.samplestamp,0.2.*ones(length(Data.fRef.samplestamp),1),'or')
legend('FECG','FQRS')
ax(3)=subplot(3,1,3);
% set(gca(), 'LineStyleOrder',markerorder, 'NextPlot','replacechildren')
plot(table2array(Data.sqi(Data.sqi.channel==ch,4)),features(:,collist,ch))
hold on
plot(table2array(Data.sqi(Data.sqi.channel==ch,4)),sqinorm(:,ch),'--b','LineWidth',2)
legend([Data.sqi.Properties.VariableNames(collist+4) 'final'])

linkaxes(ax,'x')

 xlim([[950000,958000]])
% xlim([707500 718500])



% % Selecting parts
% xlim([2.92 3].*10e4)
% cleanfigure
% matlab2tikz('standalone',true,'seg3.tikz')