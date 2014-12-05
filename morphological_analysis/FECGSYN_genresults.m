function FECGSYN_genresults(path_orig,path_ext,fs,ch)
% Input:
% path_orig: Path for original dataset
% path_ext: Path for extracted dataset
% fs: Sampling frequency
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing
% effects, foetal movement etc.).
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 20-11-2014
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

global debug

%% == Parameters
INTERV = round(0.05*fs); % BxB acceptance interval
TEMP_SAMPS = round(60*fs); % samples used for building templates
morph = 1; % turn on/off morphological analysis
%% Run through extracted datasets
cd(path_orig)
slashchar = char('/'*isunix + '\'*(~isunix));
fls_orig = dir('*.mat'); % looking for .mat (creating index)
fls_orig = arrayfun(@(x)x.name,fls_orig,'UniformOutput',false);
cd(path_ext)
fls_ext = dir('*.mat'); % looking for .mat (creating index)
fls_ext = arrayfun(@(x)x.name,fls_ext,'UniformOutput',false);

stats_ica = zeros(length(fls_orig),4);
stats_pca = zeros(length(fls_orig),4);
stats_tsc = zeros(length(fls_orig),4);
stats_tspca = zeros(length(fls_orig),4);
stats_tsekf = zeros(length(fls_orig),4);
stats_alms = zeros(length(fls_orig),4);
stats_arls = zeros(length(fls_orig),4);
stats_aesn = zeros(length(fls_orig),4);

morph_ica = zeros(length(fls_orig),2);
morph_tspca = zeros(length(fls_orig),2);
morph_tskf = zeros(length(fls_orig),2);
morph_aesn = zeros(length(fls_orig),2);
idx = cellfun(@(x) strcmp(x(1:3),'rec'),fls_ext);
fls_ext(~idx) = [];

for i = 1:length(fls_ext)
    disp(fls_ext{i})
    %= loading extracted file
    [rec,met] = strtok(fls_ext(i),'_');
    file = strcat(path_ext,fls_ext(i));
    load(file{:})
    %= loading original file
    origrec = str2double(rec{:}(4:end));
    file = strcat(path_orig,fls_orig(origrec));
    load(file{:});
    fecg = double(out.fecg{1}(ch,:)); % selecting channels
    if exist('outdata','var') % uniform naming residuals
        residual = outdata;
    end
    %= Resampling original data to match extracted (fs - if necessary)
    if size(out.mecg,2) ~= size(residual,2)
        % fref:          fetal QRS reference
        % fecgref:       fetal ECG signals (before mixture)
        fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
        fecgref = zeros(size(residual));
        for k = 1:size(fecgref,1)
            fecgref(k,:) = resample(fecg(k,:),fs,out.param.fs);
        end
    else
        fecgref = fecg;
        fref = out.fqrs{1};
    end
    [elif,~]=strtok(file{:}(end:-1:1),slashchar);
    disp(elif(end:-1:1))
    clear fecg outdata rec file elif k
    switch met{:}(2:end-4)
        case 'JADEICA'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_ica(origrec,:) = [F1,MAD,PPV,SE];
            if morph
                [qt_err,theight_err]=morpho(fecgref,residual,fref,maxch,fs,TEMP_SAMPS);
                morph_ica(origrec,:) = [mean(qt_err); mean(theight_err)];
%                 print('-dpng','-r72',[cd '/plots/' fls_ext{i} '.png'])
            end
            clear fqrs F1 MAD PPV SE
        case 'PCA'
            %= generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_pca(origrec,:) = [F1,MAD,PPV,SE];
            clear fqrs F1 MAD PPV SE
        case 'tsc'
            %= generating QRS detection statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_tsc(origrec,:) = [F1,MAD,PPV,SE];
            clear fqrs F1 MAD PPV SE
        case 'tspca'
            %= generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_tspca(origrec,:) = [F1,MAD,PPV,SE];
            %= morphological statistics
            if morph
                [qt_err,theight_err]=morpho(fecgref,residual,fref,maxch,fs,TEMP_SAMPS);
                morph_tspca(origrec,:) = [mean(qt_err); mean(theight_err)];
%                 print('-dpng','-r72',[cd '/plots/' fls_ext{i} '.png'])
            end
            clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
        case 'tsekf'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_tsekf(origrec,:) = [F1,MAD,PPV,SE];
            %= morphological statistics
            if morph
                [qt_err,theight_err]=morpho(fecgref,residual,fref,maxch,fs,TEMP_SAMPS);
                morph_tskf(origrec,:) = [mean(qt_err); mean(theight_err)];
%                 print('-dpng','-r72',[cd '/plots/' fls_ext{i} '.png'])
            end
            clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
        case 'alms'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_alms(origrec,:) = [F1,MAD,PPV,SE];
            clear fecg residual fqrs F1 MAD PPV SE
        case 'arls'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_arls(origrec,:) = [F1,MAD,PPV,SE];
            clear fecg residual fqrs F1 MAD PPV SE
        case 'aesn'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
            stats_aesn(origrec,:) = [F1,MAD,PPV,SE];
            % morphological statistics
            if morph
                [qt_err,theight_err]=morpho(fecgref,residual,fref,maxch,fs,TEMP_SAMPS);
                morph_aesn(origrec,:) = [mean(qt_err); mean(theight_err)];
%                 print('-dpng','-r72',[cd '/plots/' fls_ext{i} '.png'])
            end
            clear fecg residual fqrs F1 MAD PPV SE
    end
    clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
end

%% Statistics Generation

LWIDTH = 1.5;
FSIZE = 15;
% F1
stats_f1 = 100*[stats_ica(:,1) stats_pca(:,1) stats_tsc(:,1) stats_tspca(:,1) ...
    stats_tsekf(:,1) stats_alms(:,1) stats_arls(:,1) stats_aesn(:,1)];
h = boxplot(stats_f1,{[repmat({'BSSica'},1,length(fls_orig)) repmat({'BSSpca'},1,length(fls_orig)) ...
    repmat({'TSc'},1,length(fls_orig)) repmat({'TSpca'},1,length(fls_orig)) ...
    repmat({'TSekf'},1,length(fls_orig)) repmat({'Alms'},1,length(fls_orig)) ...
    repmat({'Arls'},1,length(fls_orig)) repmat({'Aesn'},1,length(fls_orig))]});
set(h, 'LineWidth',LWIDTH)
ylabel('F1 (%)','FontSize',FSIZE)
% MAD
figure
stats_MAD = [stats_ica(:,2) stats_pca(:,2) stats_tsc(:,2) stats_tspca(:,2) ...
    stats_tsekf(:,2) stats_alms(:,2) stats_arls(:,2) stats_aesn(:,2)];
h = boxplot(stats_MAD,{[repmat({'ICA'},1,length(fls_orig)) repmat({'PCA'},1,length(fls_orig)) ...
    repmat({'TSc'},1,length(fls_orig)) repmat({'TSpca'},1,length(fls_orig)) ...
    repmat({'EKF'},1,length(fls_orig)) repmat({'LMS'},1,length(fls_orig)) ...
    repmat({'RLS'},1,length(fls_orig)) repmat({'ESN'},1,length(fls_orig))]});
set(h, 'LineWidth',LWIDTH)
ylabel('MAD (ms)','FontSize',FSIZE)

% generate case plots
c0 = cellfun(@(x) ~isempty(regexp(x,'.c0','ONCE')),fls_orig);
c1 = cellfun(@(x) ~isempty(regexp(x,'.c1','ONCE')),fls_orig);
c2 = cellfun(@(x) ~isempty(regexp(x,'.c2','ONCE')),fls_orig);
c3 = cellfun(@(x) ~isempty(regexp(x,'.c3','ONCE')),fls_orig);
c4 = cellfun(@(x) ~isempty(regexp(x,'.c4','ONCE')),fls_orig);
c5 = cellfun(@(x) ~isempty(regexp(x,'.c5','ONCE')),fls_orig);

for met = {'tsekf' 'tspca' 'aesn' 'ica'}
    figure
   
    eval(['stat = stats_' met{:} ';']);
    statscase = [stat(c0,1) stat(c1,1) stat(c2,1) stat(c3,1) stat(c4,1) stat(c5,1)];
    h = boxplot(statscase,{[repmat({'Case 0'},1,250) repmat({'Case 1'},1,250)...
       repmat({'Case 2'},1,250) repmat({'Case 3'},1,250) repmat({'Case 4'},1,250) ...
       repmat({'Case 5'},1,250)]});
set(h, 'LineWidth',LWIDTH)
title(met)   
end

end

function [qt_err,theight_err]=morpho(fecg,residual,fqrs,maxch,fs,SAMPS)
%% Function to perform morphological analysis for TS/BSS extracted data
%
% >Inputs
% fecg: Propagated fetal signal before mixture with noise sources
% residual: Result of fetal extraction from abdominal signals
% fqrs: Reference fetal QRS samplestamps
% maxch: Channel with highest F1-accuracy for FQRS detections
% SAMPS: Number of samples used for generating templates
%
% > Outputs
% qt_err: Array containing QT error for each template
% theight_err: Array containing T-height error for each template
%

% generating reference template
W = fecg*pinv(residual);
srcfecg = W*fecg;
qt_err = zeros(length(residual)/SAMPS,1); theight_err = zeros(length(residual)/SAMPS,1);
block = 1;
for j = 1:SAMPS:length(residual)
    % checking borders
    if j+SAMPS > length(residual)
        endsamp = length(residual);
    else
        endsamp = j + SAMPS -1;
    end
    % qrs complexes in interval
    qrstmp = fqrs(fqrs>j&fqrs<endsamp)-j;
    % abdominal signal template
    temp_abdm = FECGSYN_tgen(residual(maxch(block),j:endsamp),qrstmp);
    % reference template
    temp_ref = FECGSYN_tgen(srcfecg(maxch(block),j:endsamp),qrstmp);
    temp_abdm = temp_abdm.avg; temp_ref = temp_ref.avg;
    % evaluating morphological features
    [qt_err(block),theight_err(block)] = FECGSYN_manalysis(temp_abdm,temp_ref,fs);
    block = block+1;
end
end




%xlim([1500 2000])
%print('-dpng','-r72',['/media/fernando/FetalEKG/2014.10_fecgsyn_simulations(5.0)/extracted3Hz/plots/' fls_ext{i} '.png'])

% % % Plot about cases
% % if debug
% % % This script plots boxplots with 2 groups
% % % Boxplot multicolor
% % stats_ica(:,[1 3 4]) = 100.*stats_ica(:,[1 3 4]);
% % stats_tsc(:,[1 3 4]) = 100.*stats_tsc(:,[1 3 4]);
% % LINE_WIDTH = 1.2;
% % FONT_SIZE = 12;
% % MARKER_SIZE = 7;
% % fig1=figure(1)
% % base = cellfun(@(x) ~isempty(regexp(x, '_l\d.mat$', 'match')), fls_orig);
% % c0 = cellfun(@(x) ~isempty(regexp(x, '_c0.mat$', 'match')), fls_orig);
% % c1 = cellfun(@(x) ~isempty(regexp(x, '_c1.mat$', 'match')), fls_orig);
% % c2 = cellfun(@(x) ~isempty(regexp(x, '_c2.mat$', 'match')), fls_orig);
% % c3 = cellfun(@(x) ~isempty(regexp(x, '_c3.mat$', 'match')), fls_orig);
% % c4 = cellfun(@(x) ~isempty(regexp(x, '_c4.mat$', 'match')), fls_orig);
% % c5 = cellfun(@(x) ~isempty(regexp(x, '_c5.mat$', 'match')), fls_orig);
% %
% % stats_ica = stats_ica(:,1);
% % stats_tsc = stats_tsc(:,1);
% %
% % N = sum(base)+sum(c0)+sum(c1)+sum(c2)+sum(c3)+sum(c4)+sum(c5);
% % bh=boxplot([stats_ica(base); stats_ica(c0); stats_ica(c1); stats_ica(c2); stats_ica(c3); stats_ica(c4); stats_ica(c5);...
% % stats_tsc(base); stats_ica(c0); stats_tsc(c1); stats_tsc(c2); stats_tsc(c3); stats_tsc(c4); stats_tsc(c5)], ...
% % {[repmat({'Base'},1,sum(base)) repmat({'Case 0'},1,sum(c0)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ...
% % repmat({'Base'},1,sum(base)) repmat({'Case 0'},1,sum(c0)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ] ...
% % [repmat({'ICA'},N,1);repmat({'TS'},N,1)]},'colors',repmat('rk',1,6),'factorgap',10,'labelverbosity','minor','labelorientation','inline')
% % set(bh,'linewidth',LINE_WIDTH);
% % ylabel('F1 (in %)','FontSize',FONT_SIZE)
% % xlabel('Recording','FontSize',FONT_SIZE)
% % set(gca,'FontSize',FONT_SIZE)
% % set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
% % % ylim([0 105])
% % xl=xlabel('Recording Number','FontSize',FONT_SIZE);
% % pos=get(xl,'Pos');
% % set(xl,'Pos',[pos(1) pos(2)-30 pos(3)])
% % save2pdf('boxplot',fig1,600)
% %
% % % Plot about SNRmn
% % a = cellfun(@(x) strsplit(x,'_snr'), fls_orig,'UniformOutput',0);
% %
% % b = cellfun(@(x) length(x)>1,a);
% % end