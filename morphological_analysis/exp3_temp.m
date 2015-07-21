% This script runs through comparison between QT and Th values for applying
% or not the 
cd('/home/andreotti/Desktop/Simulator Paper/bssdist')
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

pca = 1:8:14400; %fix
ica = 2:8:14400; %fix
srcnan = 0;
srctot = 0;
timenan = 0;
timetot = 0;
exp3tab = [];
exp3ica = [];
for i = 1:length(fls)
    load(fls{i})
    
     if ismember(str2double(strtok(fls{i},'_')),pca)
%         color = 'b';
%         figure(1)
%         hold on
        exp3ica(5*(i-1)+1:5*(i-1)+5,1) = false;
    elseif ismember(str2double(strtok(fls{i},'_')),ica)
        exp3ica(5*(i-1)+1:5*(i-1)+5,1) = true;
%         color = 'r';
%         figure(2)
%         hold on
    else
        color = 'k';
     end
    %% QT
    qt_ref(cellfun(@isempty, qt_ref)) = {NaN};
    qt_ref2(cellfun(@isempty, qt_ref2)) = {NaN};
    th_ref(cellfun(@isempty, th_ref)) = {NaN};
    th_ref2(cellfun(@isempty, th_ref2)) = {NaN};  
%     subplot(2,1,1)
    % Counting NaNs
    srcnan = srcnan + sum(sum(isnan(cell2mat(qt_ref))));
    srctot = srctot + numel(cell2mat(qt_ref));
    timenan = timenan + sum(sum(isnan(cell2mat(qt_ref2))));
    timetot = timetot + numel(cell2mat(qt_ref2));
    
    % Max FQT per segment
    tmpres = min(cell2mat(qt_ref)-cell2mat(qt_ref2));
    exp3tab = [exp3tab; tmpres'];
    
    clear tmpresßß
    
%     hold on
%     
%     plot(nanmean(),nanmean(cell2mat(qt_ref2)),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
%     hold off

%     
%     %% TH
%     subplot(2,1,2)
%     hold on
%     plot(nanmean(cell2mat(th_ref)),nanmean(cell2mat(th_ref2)),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
%     hold off

    
    
   clear qt_ref qt_ref2 th_ref th_ref2
end

[h,p]=ttest(exp3tab(:,1),exp3tab(:,2));
if h
    fprintf('Null hypothesis can be rejected with p= %d \n',p);
else
    fprintf('Null hypothesis CANNOT be rejected');
end
fprintf('NaNs on source domain:  %3.2f percent \n',srcnan/srctot*100);
fprintf('NaNs on time domain:  %3.2f percent \n',timenan/timetot*100);

figure(1)
% subplot(2,1,1)
color = 'b';
plot(exp3tab(find(exp3ica),1),exp3tab(find(exp3ica),2),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
xlim([100 260]),ylim([100 260])
xlabel('FQT (time domain)')
ylabel('FQT (source domain)')  
title('BSS_{ica}')
matlab2tikz(['fqtica.tikz'], 'height', '\figureheight', 'width', '\figurewidth');

figure(2)
color = 'r';
plot(exp3tab(find(~exp3ica),1),exp3tab(find(~exp3ica),2),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
xlim([100 260]),ylim([100 260])
xlabel('FQT (time domain)')
ylabel('FQT (source domain)')  
title('BSS_{pca}')
matlab2tikz(['fqtpca.tikz'], 'height', '\figureheight', 'width', '\figurewidth');

% %     subplot(2,1,2)
% %     xlabel('T_h (time domain)')
% %     ylabel('T_h (spectral domain)')
% %     xlim([0 50]),ylim([0 50])
%     
%     figure(2)
% %     subplot(2,1,1)
%     xlim([100 300]),ylim([100 300])
%     xlabel('FQT (time domain)')
%     ylabel('FQT (spectral domain)')    
% %     subplot(2,1,2)
%     xlabel('T_h (time domain)')
%     ylabel('T_h (spectral domain)')
%     xlim([0 50]),ylim([0 50])