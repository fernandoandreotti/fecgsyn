function FECGSYN_genresults_exp1(path_orig,path_ext,fs)
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

for i = 1:length(fls_ext)
    disp(fls_ext{i})
    %= loading extracted file
    rec = regexp(fls_ext(i),'rec+[0-1]?[0-9]?[0-9]?[0-9]','match');
    file = strcat(path_ext,fls_ext(i));
    met = strcmp(char('DEF', 'PCA', 'ICA'),file(end-7:end-5));
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
    % getting statistics
    [F1,MAD,PPV,SE] = getStats(out,residual,fs_new);
    switch met{:}(2:end-4)
        case 'JADEICA'
            % generating statistics
            stats_ica(origrec,:) = [F1,MAD,PPV,SE];
            if morph
                [qt_err,theight_err]=morpho(fecgref,residual,fref,maxch,fs,TEMP_SAMPS);
                morph_ica(origrec,:) = [mean(qt_err); mean(theight_err)];
                %                 print('-dpng','-r72',[cd '/plots/' fls_ext{i} '.png'])
            end
            clear fqrs F1 MAD PPV SE
        case 'PCA'
            %= generating statistics
            stats_pca(origrec,:) = [F1,MAD,PPV,SE];
            clear fqrs F1 MAD PPV SE
        case 'tsc'
            %= generating QRS detection statistics
            stats_tsc(origrec,:) = [F1,MAD,PPV,SE];
            clear fqrs F1 MAD PPV SE
        case 'tspca'
            %= generating statistics
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
            stats_alms(origrec,:) = [F1,MAD,PPV,SE];
            clear fecg residual fqrs F1 MAD PPV SE
        case 'arls'
            % generating statistics
            stats_arls(origrec,:) = [F1,MAD,PPV,SE];
            clear fecg residual fqrs F1 MAD PPV SE
        case 'aesn'
            % generating statistics
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
save ['wkspace_' date]

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

for met = {'pca' 'FASTICA_DEF' 'JADEICA' 'FASTICA_SYM'}
    figure   
    eval(['stat = stats_struct{4,1}.stats_' met{:} ';']);
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

function [F1,RMS,PPV,SE] = getStats(out,residual,fs_new)
    INTERV = round(0.05*fs_new);    % BxB acceptance interval

    % detect QRS on each channel
    fqrs = cell(1,size(mixture,1));
    for j = 1:length(ch)
        fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
    end
    % creating statistics in 1-min blocks

    min = 1;
    maxch = zeros(1,length(residual)/fs_new/60);
    fqrs_temp = cell(1,length(residual)/fs_new/60);
    while min <= length(residual)/fs_new/60;
        F1max = 0;
        idxref = (out.fqrs{1}>=(min-1)*fs_new*60+1)&(out.fqrs{1}<=min*fs_new*60);
        for j = 1:length(ch)
            idx = (fqrs{j}>=(min-1)*fs_new*60+1)&(fqrs{j}<=min*fs_new*60);
            [F1,~,~,~] = Bxb_compare(out.fqrs{1}(idxref),fqrs{j}(idx),INTERV);
            if F1 > F1max    % compare and see if this channel provides max F1
                maxch(min) = j;
                F1max = F1;
                fqrs_temp{min} = fqrs{j}(idx);%+ (min-1)*fs_new*60;    % adding fqrs detections to temporary cell
            end
        end
        min = min+1;
    end
    fqrs = cell2mat(fqrs_temp);
    [F1,RMS,PPV,SE] = Bxb_compare(out.fqrs{1}(idxref),fqrs{j}(idx),INTERV);
    
    
end

