function FECGSYN_genresults(path_orig,fs,ch,exp3)
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing
% effects, foetal movement etc).
% Experiment 2 - Statistics (F1,MAE,PPV,SE)
% Experiment 3 - Morphologycal Analysis
%
% Input:
% path_orig:        Path for original dataset
% fs:               Sampling frequency
% ch:               Channels to be used
% exp3:            Boolean, if 0 runs exp2 and 1 exp3
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 30-05-2014
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

global debug filesproc
slashchar = char('/'*isunix + '\'*(~isunix));
if debug,  mkdir([path_orig 'plots' slashchar]), end

%% == Parameters
INTERV = round(0.05*fs); % BxB acceptance interval
TEMP_SAMPS = round(60*fs); % samples used for building templates
fs_new = 250;

%% Run through extracted datasets
cd(path_orig)

% abdominal mixtures
fls_orig = dir([path_orig '*.mat']); % looking for .mat (creating index)
fls_orig = arrayfun(@(x)x.name,fls_orig,'UniformOutput',false);
% extracted files
if ~exp3
    path_ext = [path_orig 'exp2' slashchar];
    fls_ext = dir([path_ext '*.mat']); % looking for .mat (creating index)
else
    path_ext = [path_orig 'exp3' slashchar];
    fls_ext = dir([path_ext '*.mat']); % looking for .mat (creating index)
    % Preprocessing Filter coefficiens
    % high-pass filter
    Fstop = 0.5;  % Stopband Frequency
    Fpass = 1;    % Passband Frequency
    Astop = 20;   % Stopband Attenuation (dB)
    Apass = 0.1;  % Passband Ripple (dB)
    h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, fs_new);
    Hhp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf', ...
        'SystemObject', true);
    [b_hp,a_hp] = tf(Hhp);
    % low-pass filter
    Fpass = 100;   % Passband Frequency
    Fstop = 110;  % Stopband Frequency
    Astop = 20;   % Stopband Attenuation (dB)
    Apass = 0.1;    % Passband Ripple (dB)
    h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, fs_new);
    Hlp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf');
    [b_lp,a_lp] = tf(Hlp);
    clear Fstop Fpass Astop Apass h Hhp Hlp
    
end
fls_ext = arrayfun(@(x)x.name,fls_ext,'UniformOutput',false);
idx = cellfun(@(x) strcmp(x(1:3),'rec'),fls_ext); % ignoring files not begining with 'rec'
fls_ext(~idx) = [];

% pre-allocation
stats.JADEICA = zeros(length(fls_orig),4);
stats.PCA = zeros(length(fls_orig),4);
stats.tsc = zeros(length(fls_orig),4);
stats.tspca = zeros(length(fls_orig),4);
stats.tsekf = zeros(length(fls_orig),4);
stats.alms = zeros(length(fls_orig),4);
stats.arls = zeros(length(fls_orig),4);
stats.aesn = zeros(length(fls_orig),4);
morph.JADEICA = cell(length(fls_orig),7);
morph.PCA = cell(length(fls_orig),7);
morph.tsc = cell(length(fls_orig),7);
morph.tspca = cell(length(fls_orig),7);
morph.tsekf = cell(length(fls_orig),7);
morph.alms = cell(length(fls_orig),7);
morph.arls = cell(length(fls_orig),7);
morph.aesn = cell(length(fls_orig),7);

% = Runs through list of extracted files
for i = filesproc%length(fls_ext)
% for i = randperm(length(fls_ext))
    disp(fls_ext{i})
    fprintf('Data %d out of %d \n',i,length(fls_ext));

    %= loading extracted file
    [rec,met] = strtok(fls_ext(i),'_');
    file = strcat(path_ext,fls_ext(i));
    load(file{:})
    %= loading original file
    origrec = str2double(rec{:}(4:end));
    file = strcat(path_orig,fls_orig(origrec));
    load(file{:});
    fecg = double(out.fecg{1}(ch,:)); % selecting channels
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
    
    % Figuring out which extraction method was used, possibilities are:
    % (JADEICA,PCA,tsc,tspca,tsekf,alms,arls,aesn)
    method = met{:}(2:end-4);
    
    %= Getting statistics (exp 2)
    if ~exp3
        [F1,MAE,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
        MAE = MAE*1000/fs_new;
        stats.(method)(origrec,:) = [F1,MAE,PPV,SE]; % dynamic naming
        %= Getting statistics (exp 3)
    else
        if ~exist([path_orig 'wfdb'],'dir')
            mkdir([path_orig 'wfdb'])
        end
        cd([path_orig 'wfdb'])
        bss = strcmp(method,'JADEICA')|strcmp(method,'PCA'); % apply coordinate transformation or not
        fname = [path_orig 'plots' slashchar fls_ext{i}(1:end-4)];
        [outputs{1:7}]= morpho(fecgref,residual,fref,fs,TEMP_SAMPS,bss,fname,[b_hp,a_hp,b_lp,a_lp],i);
        morph.(method)(origrec,:) = outputs;
    end
    clear fecg residual fqrs F1 MAE PPV SE qt_err theight_err
end

save([path_orig 'wksp' num2str(filesproc(end))])

%% Plots and statistics generation
if debug
    LWIDTH = 1.5;
    FSIZE = 15;
    %== F1 plot
    figure
    stats_f1 = 100*[stats_ica(:,1) stats_pca(:,1) stats_tsc(:,1) stats_tspca(:,1) ...
        stats_tsekf(:,1) stats_alms(:,1) stats_arls(:,1) stats_aesn(:,1)];
    h = boxplot(stats_f1);
    set(gca,'XTick',[1:8])  % This automatically sets
    set(gca,'XTickLabel',{'BSSica';'BSSpca';'TSc';'TSpca';'TSekf';'Alms';'Arls';'Aesn'})
    set(h, 'LineWidth',LWIDTH)
    ylabel('F_1 (%)','FontSize',FSIZE)
    xlabel('Method','FontSize',FSIZE)
    h=gca;
    rotateticklabel(h,45);
    set(gca,'FontSize',FSIZE)
    set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE)
    
    % MAE
    figure
    stats_MAE = [stats_ica(:,2) stats_pca(:,2) stats_tsc(:,2) stats_tspca(:,2) ...
        stats_tsekf(:,2) stats_alms(:,2) stats_arls(:,2) stats_aesn(:,2)];
    h = boxplot(stats_MAE);
    set(gca,'XTick',1:8)  % This automatically sets
    set(gca,'XTickLabel',{'BSSica';'BSSpca';'TSc';'TSpca';'TSekf';'Alms';'Arls';'Aesn'})
    set(h, 'LineWidth',LWIDTH)
    ylabel('MAE (ms)','FontSize',FSIZE)
    xlabel('Method','FontSize',FSIZE)
    h=gca;
    rotateticklabel(h,45);
    set(gca,'FontSize',FSIZE)
    set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE)
    
    % generate case plots
    c0 = cellfun(@(x) ~isempty(regexp(x,'.c0','ONCE')),fls_orig);
    c1 = cellfun(@(x) ~isempty(regexp(x,'.c1','ONCE')),fls_orig);
    c2 = cellfun(@(x) ~isempty(regexp(x,'.c2','ONCE')),fls_orig);
    c3 = cellfun(@(x) ~isempty(regexp(x,'.c3','ONCE')),fls_orig);
    c4 = cellfun(@(x) ~isempty(regexp(x,'.c4','ONCE')),fls_orig);
    c5 = cellfun(@(x) ~isempty(regexp(x,'.c5','ONCE')),fls_orig);
    base = ~(c0|c1|c2|c3|c4|c5);
    
    % Generate Table
    counter1 = 1;
    table = zeros(16,28);
    for met = {'ica' 'pca' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
        eval(['stat = stats_' met{:} ';']);
        % F1
        statscase = 100*[stat(base,1) stat(c0,1) stat(c1,1) stat(c2,1) stat(c3,1) stat(c4,1) stat(c5,1)];
        auxtab = [median(statscase)',-1.*ones(7,1),std(statscase)',-2.*ones(7,1)];
        table(counter1,:) = reshape(auxtab',1,7*4);
        counter1 = counter1 + 1;
        
        % MAE
        statscase = [stat(base,2) stat(c0,2) stat(c1,2) stat(c2,2) stat(c3,2) stat(c4,2) stat(c5,2)];
        auxtab = [median(statscase)',-1.*ones(7,1),std(statscase)',-2.*ones(7,1)];
        table(counter1,:) = reshape(auxtab',1,7*4);
        counter1 = counter1 + 1;
    end
    table = round(table.*10)./10;
    % F1
    c=1;
    for met = {'ica' 'aesn'}%{'ica' 'pca' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
        figure(c)
        c = c+1;
        eval(['stat = stats_' met{:} ';']);
        statscase = 100*[stat(base,1) stat(c0,1) stat(c1,1) stat(c2,1) stat(c3,1) stat(c4,1) stat(c5,1)];
        h = boxplot(statscase,{});
        set(gca,'XTick',1:7)  % This automatically sets
        set(gca,'XTickLabel',{'Baseline','Case 0','Case 1','Case 2','Case 3','Case 4','Case 5'})
        set(h, 'LineWidth',LWIDTH)
        ylabel('F_1 (%)','FontSize',FSIZE)
        title(met)
        h=gca;
        rotateticklabel(h,45);
        set(gca,'FontSize',FSIZE);
        set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE);
    end
    
    % MAE
    
    for met = {'ica' 'aesn'}%{'tsekf' 'tspca' 'aesn' 'ica'}
        figure(c)
        c= c+1;
        eval(['stat = stats_' met{:} ';']);
        statscase = [stat(base,2) stat(c0,2) stat(c1,2) stat(c2,2) stat(c3,2) stat(c4,2) stat(c5,2)];
        h = boxplot(statscase,{});
        set(gca,'XTick',1:7)  % This automatically sets
        set(gca,'XTickLabel',{'Baseline','Case 0','Case 1','Case 2','Case 3','Case 4','Case 5'})
        set(h, 'LineWidth',LWIDTH)
        h = findobj('Tag','Box');
        set(h,'Color',([187 81 112]./255));
        ylabel('MAE (ms)','FontSize',FSIZE)
        h=gca;
        set(gca,'FontSize',FSIZE);
        set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE);
    end
end

end

function [qt_test,qt_ref,th_test,th_ref,qt_err,theight_err,numbNaN]=...
    morpho(fecg,residual,fqrs,fs,SAMPS,bss,fname,filterc,filen)
%% Function to perform morphological analysis for TS/BSS extracted data
%
% >Inputs
% fecg:         Propagated fetal signal before mixture with noise sources
% residual:     Result of fetal extraction from abdominal signals
% fqrs:         Reference fetal QRS samplestamps
% SAMPS:        Number of samples used for generating templates
% bss:          Boolean true if using BSS technique
% fname:        Filename to be used in saving plots
% filterc:      Filter coefficients [b_hp,a_hp,b_lp,a_lp] being
%               highpass (hp) and lowpass (lp)
%
% > Outputs
% qt_err: Array containing QT error for each template
% theight_err: Array containing T-height error for each template
%
global debug
% if bss, propagating reference to source domain
if bss
    W = fecg*pinv(residual);
    srcfecg = W*fecg;
else
    srcfecg = fecg;
end

% Allocatting
qt_test = cell(size(residual,1),length(residual)/SAMPS,1);
qt_ref = qt_test;
th_test = qt_test;
th_ref = qt_test;
qt_err = qt_test;
theight_err = qt_test;
%= Block-wise calculation and template generation
for ch = 1:size(residual,1)
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
        [temp_abdm,qrs_abdm,status1] = FECGSYN_tgen(residual(ch,j:endsamp),qrstmp,fs);
        % reference template
        [temp_ref,qrs_ref,status2] = FECGSYN_tgen(srcfecg(ch,j:endsamp),qrstmp,fs);
        temp_abdm = temp_abdm.avg; temp_ref = temp_ref.avg;
        
        if (~status1||~status2)
            qt_test{ch,block} = NaN;
            qt_ref{ch,block} = NaN;
            th_test{ch,block} = NaN;
            th_ref{ch,block} = NaN;
            qt_err{ch,block} = NaN;
            theight_err{ch,block} = NaN;
        else
            % evaluating morphological features
            [qt_test{ch,block},qt_ref{ch,block},th_test{ch,block},th_ref{ch,block},...
                qt_err{ch,block},theight_err{ch,block}] = FECGSYN_manalysis(temp_abdm,temp_ref,qrs_abdm,qrs_ref,fs,filterc,filen);
        end
        
        if debug && ~isnan(qt_test{ch,block}) && ~isnan(qt_ref{ch,block})
            try
                drawnow
                print('-dpng','-r72',[fname '_ch' num2str(ch) '_s' num2str(block) '.png'])
            catch
                warning('Failed to save plot')
            end
            
        end
        block = block+1;
    end
end
% Figuring out how many NaNs were output per channel
try
    id1 = cellfun(@(x) isnan(x),qt_test);
    id2 = cellfun(@(x) isnan(x),qt_ref);
    numbNaN=sum(sum(id1|id2));
catch
    disp
end
end
