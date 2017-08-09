function FECGSYN_benchFQRS_plot(stats,fls_orig)
% function FECGSYN_benchFQRS_plot(stats,fls_orig)
% This function produces plots for Experiment 2 of Andreotti et al 2016. 
% 
% Input:
%  stats            Structure containing statistic results for FQRS
%                   detections
%  fls_orig         Files containing original data (prior to extraction).
%                   See FECGSYN_benchFQRS.m
% 
% See also:
% FECGSYN_benchFQRS
% 
% --
% fecgsyn toolbox, version 1.2, Jan 2017
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% University of Oxford, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
% 
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. 
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. 
%   Physiological Measurement.35 1537-1550. 2014.
%
% Last updated : 10-03-2016
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Plots and statistics generation

%== Hardcoded bits..

% methods used
met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }; 
Nmet = 8; %number of methods
% cases and snrs used..
c0 = cellfun(@(x) ~isempty(regexp(x,'.c0','ONCE')),fls_orig);
c1 = cellfun(@(x) ~isempty(regexp(x,'.c1','ONCE')),fls_orig);
c2 = cellfun(@(x) ~isempty(regexp(x,'.c2','ONCE')),fls_orig);
c3 = cellfun(@(x) ~isempty(regexp(x,'.c3','ONCE')),fls_orig);
c4 = cellfun(@(x) ~isempty(regexp(x,'.c4','ONCE')),fls_orig);
c5 = cellfun(@(x) ~isempty(regexp(x,'.c5','ONCE')),fls_orig);
base = ~(c0|c1|c2|c3|c4|c5);
casename = {'base','c0','c1','c2','c3','c4','c5'};
casename(sum([base c0 c1 c2 c3 c4 c5]) == 0) = [];
Ncases = length(casename); % number of cases
snr00 = cellfun(@(x) ~isempty(regexp(x,'.snr00dB','ONCE')),fls_orig);
snr03 = cellfun(@(x) ~isempty(regexp(x,'.snr03dB','ONCE')),fls_orig);
snr06 = cellfun(@(x) ~isempty(regexp(x,'.snr06dB','ONCE')),fls_orig);
snr09 = cellfun(@(x) ~isempty(regexp(x,'.snr09dB','ONCE')),fls_orig);
snr12 = cellfun(@(x) ~isempty(regexp(x,'.snr12dB','ONCE')),fls_orig);
snrname =  {'snr00' 'snr03' 'snr06' 'snr09' 'snr12'};
snrname(sum([snr00 snr03 snr06 snr09 snr12]) == 0) = [];
Nsnr = length(snrname); % number of SNRs

if (Nsnr==0)||(Ncases==0)
    error('File name not following the standard FECGSYN nomenclature. Please fix dataset.')    
end


%% Generating Tables
% Generate Table for cases (different SNRs are average)
% table format is simplified for publishing, -1 represents $\pm$ (plus/minus symbol)
% while -2 can be substituted by "&" in latex (column separation)
counter1 = 1;
table = zeros(2*Nmet,4*Ncases); % 8 methods (2 rows each), 7 cases (4 columns each)
for k = 1:Nmet    
    eval(['stat = stats.' met{k} ';']);
    % F1
    statscase = 100*[stat(base,1) stat(c0,1) stat(c1,1) stat(c2,1) stat(c3,1) stat(c4,1) stat(c5,1)];
    auxtab = [mean(statscase)',-1.*ones(Ncases,1),std(statscase)',-2.*ones(Ncases,1)];
    auxtab2(counter1,:) = median(statscase)';
    table(counter1,:) = reshape(auxtab',1,Ncases*4);
    counter1 = counter1 + 1;
    
    % MAE
    statscase = [stat(base,2) stat(c0,2) stat(c1,2) stat(c2,2) stat(c3,2) stat(c4,2) stat(c5,2)];
    auxtab = [mean(statscase)',-1.*ones(Ncases,1),std(statscase)',-2.*ones(Ncases,1)];
    table(counter1,:) = reshape(auxtab',1,Ncases*4);
    counter1 = counter1 + 1;
end
table = round(table.*10)./10; %rounding numbers

%= Boxplot for cases and SNR
figure
count2 = 1;
colors = rand(35,3); 
for k = 1:Nmet
    eval(['stattmp = stats.'  met{k} ';']);
    statscasesnr = NaN(size(stattmp,1)/Ncases/Nsnr,Ncases*Nsnr); 

    count1 = 1;
    for cs = 1:length(casename)
        cases = casename(cs);
        for snrs = 1:length(snrname)
            snr = snrname(snrs);
            try
                statscasesnr(:,count1) = 100*stattmp(eval(cases{:})&eval(snr{:}),1);
            catch
                disp(['No available case ' cases{:} ' or SNR ' snr{:} '. Skipping'])
            end
            count1 = count1 + 1;
        end
    end
    labelssnr = repmat(1:Nsnr,1,Ncases);
    labelscase = reshape(repmat(casename,Nsnr,1),1,Ncases*Nsnr);
    subplot(2,4,count2)
    boxplot(statscasesnr,{labelscase labelssnr},'factorgap',3,'color',colors,...
        'medianstyle','target','plotstyle','compact','boxstyle','filled') 
    h=findobj(gca,'tag','Outliers'); % not ploting outliers
    delete(h) 
    for i = 1:Ncases % go through each case and do a Kruskal-Wallis test to check significance
        p(count2,i)=kruskalwallis(statscasesnr(:,(i-1)*Nsnr+1:Nsnr*i),[],'off');
    end
    
    count2 = count2 +1;
end   
    
%% Generate Boxplots
% Loading data

%=F1
% re-calculating to include baseline
statsf1 = zeros(Nsnr,Ncases,8); statsmae = statsf1;
count1 = 1;
for k = 1:Nmet
    eval(['stat = stats.' met{k} ';']);
    count2 = 1;
     for snrs = 1:length(snrname)
         snr = snrname(snrs);
         snrloop = eval(snr{:});
         statsf1(count2,1:Ncases,count1) = median(100*[stat(base&snrloop,1) stat(c0&snrloop,1)...
             stat(c1&snrloop,1) stat(c2&snrloop,1) stat(c3&snrloop,1) ...
             stat(c4&snrloop,1) stat(c5&snrloop,1)]); % F1
          statsmae(count2,1:Ncases,count1) = median([stat(base&snrloop,2) stat(c0&snrloop,2)...
             stat(c1&snrloop,2) stat(c2&snrloop,2) stat(c3&snrloop,2) ...
             stat(c4&snrloop,2) stat(c5&snrloop,2)]);% MAE
         count2 = count2 + 1; 
     end
     count1 = count1 + 1;
end

% Evaluating across different SNRs
statsf1(:,1,:) = []; statsmae(:,1,:) = []; % removing baseline since it does not make sense on significance analysis
count1 = 1;
for var = {'statsf1' 'statsmae'}
    stastuse = eval(var{:});    
    figure
    if any(any(any(isnan(stastuse))))
        continue; % loop if there are NaNs
    end
    for snr = 1:size(stastuse,1)
        tempstat = reshape(stastuse(snr,:,:),[],8);
        p(snr,1) = friedman(tempstat,1,'off');
        p(snr,2) = friedman(tempstat',1,'off');
        for i = 1:8
            for j = 1:8
                [psig(snr,i,j),hsig(snr,i,j)] = signtest(tempstat(:,i),tempstat(:,j));
            end
        end
        gridp = reshape(psig(snr,:,:),8,8);
        gridp(gridp >= 0.05) = 0;
        gridp(gridp < 0.01&gridp>0) =  2;
        gridp(gridp < 0.05&gridp>0) =  1;
        subplot(1,Nsnr,snr)
        pcolor([gridp NaN(8,1);NaN(1,9)])%,[0 2])
        colormap(flipud(gray))
        %xlim([0 8]),ylim([0 8])
        set(gca,'xtick', linspace(0.5,8.5,8), 'ytick', linspace(0.5,8.5,8));
        set(gca,'xticklabel', {[1:8]}, 'yticklabel', {[1:8]});
        axis square
        %     set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
        grid on
    end
end
if exist('psig','var'), disp(psig),end
clear p statuse gridp statsf1 statsmae psig hsig


% Evaluating across different cases (baseline not included)
statsf1 = zeros(Ncases-1,Nsnr,Nmet); statsmae = statsf1;
count1 = 1;
for k=1:Nmet
    eval(['stat = stats.' met{k} ';']);
    count2 = 1;
     for cs = 1:length(casename)
         cases = casename(cs);
         caseloop = eval(cases{:});
         statsf1(count2,1:Nsnr,count1) = median(100*[stat(snr00&caseloop,1)...
             stat(snr03&caseloop,1) stat(snr06&caseloop,1) stat(snr09&caseloop,1) ...
             stat(snr12&caseloop,1)]); % F1
          statsmae(count2,1:Nsnr,count1) = median([stat(snr00&caseloop,2)...
             stat(snr03&caseloop,2) stat(snr06&caseloop,2) stat(snr09&caseloop,2) ...
             stat(snr12&caseloop,2)]);% MAE
         count2 = count2 + 1; 
     end
     count1 = count1 + 1;
end

count1 = 1;
for var = {'statsf1' 'statsmae'}
    stastuse = eval(var{:});
    figure
    for snr = 1:size(stastuse,1)
        tempstat = reshape(stastuse(snr,:,:),[],Nmet);
        try % sometimes NaNs occur
            p(snr,1) = friedman(tempstat,1,'off');
        catch
            p(snr,1) = 1;
        end
        try
            p(snr,2) = friedman(tempstat',1,'off');
        catch
            p(snr,2) = 1;
        end
        for i = 1:Nmet
            for j = 1:Nmet
                [psig(snr,i,j),hsig(snr,i,j)] = signtest(tempstat(:,i),tempstat(:,j));
            end
        end
        gridp = reshape(psig(snr,:,:),Nmet,Nmet);
        gridp(gridp >= 0.05) = 0;
        gridp(gridp < 0.01&gridp>0) =  2;
        gridp(gridp < 0.05&gridp>0) =  1;
        subplot(1,size(stastuse,1),snr)
        pcolor([gridp NaN(Nmet,1);NaN(1,9)])%,[0 2])
        colormap(flipud(gray))
        %xlim([0 Nmet]),ylim([0 Nmet])
        set(gca,'xtick', linspace(0.5,Nmet + 0.5,Nmet), 'ytick', linspace(0.5,Nmet + 0.5,Nmet));
        set(gca,'xticklabel', {[1:Nmet]}, 'yticklabel', {[1:Nmet]});
        axis square
        %     set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
        grid on
    end
end
disp(psig)
