% This script runs through comparison between QT and Th values for applying
% or not the 
cd('/home/andreotti/Desktop/Simulator Paper/bssdist')
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

pca = 1:8:14400; %fix
ica = 2:8:14400; %fix

for i = 1:length(fls)
    load(fls{i})
    
    if ismember(str2double(strtok(fls{i},'_')),pca)
        color = 'b';
        figure(1)
        hold on
    elseif ismember(str2double(strtok(fls{i},'_')),ica)
        color = 'r';
        figure(2)
        hold on
    else
        color = 'k';
    end
    %% QT
    qt_ref(cellfun(@isempty, qt_ref)) = {NaN};
    qt_ref2(cellfun(@isempty, qt_ref2)) = {NaN};
    th_ref(cellfun(@isempty, th_ref)) = {NaN};
    th_ref2(cellfun(@isempty, th_ref2)) = {NaN};  
    subplot(2,1,1)
    hold on
    plot(nanmean(cell2mat(qt_ref)),nanmean(cell2mat(qt_ref2)),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
    hold off
    xlabel('FQT (time domain)')
    ylabel('FQT (spectral domain)')
    
    %% TH
    subplot(2,1,2)
    hold on
    plot(nanmean(cell2mat(th_ref)),nanmean(cell2mat(th_ref2)),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',3)
    hold off
    xlabel('T_h (time domain)')
    ylabel('T_h (spectral domain)')
    
    
   clear qt_ref qt_ref2 th_ref th_ref2
end
    