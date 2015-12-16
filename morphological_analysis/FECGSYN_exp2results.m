% function FECGSYN_exp2results(path_orig,fs,ch,exp3)
%
% this script generates statistics for the Experiment 2 from Andreotti2016,
% namely F1,MAE,PPV and SE.
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 16-12-2014
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


% == Parameters
path = 'D:\Users\Andreotti\Dropbox\sharelatex\[Andreotti et al 2016] Standardising FECGSYN (in review)\fecgsyn_morphology\mats\';
INTERV = round(0.05*fs); % BxB acceptance interval
TEMP_SAMPS = round(60*fs); % samples used for building templates
fs_new = 250;

% Run through extracted datasets
cd(path)
load('exp2.mat')

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
    Fpass = 90;   % Passband Frequency
    Fstop = 100;  % Stopband Frequency
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

% = Runs through list of extracted files
for i = filesproc%length(fls_ext)
    % for i = randperm(length(fls_ext))
    disp(fls_ext{i})
    fprintf('Data %d out of %d \n',i,length(fls_ext));
    
    %= loading extracted file
    [rec,met] = strtok(fls_ext(i),'_');
    % Figuring out which extraction method was used, possibilities are:
    % (JADEICA,PCA,tsc,tspca,tsekf,alms,arls,aesn)
    method = met{:}(2:end-4);
    if ~strcmp(method,'JADEICA')
        continue
    end
    file = strcat(path_ext,fls_ext(i));
    load(file{:})
    %= loading original file
    origrec = str2double(rec{:}(4:end));
    file = strcat(path_orig,fls_orig(origrec));
    cas = regexp(file{:},'_c[0-7]','match'); % find out which case it depicts
    if isempty(cas)
        cas = {'bas'};
    end
    
    load(file{:});
    fecg = double(out.fecg{1}(ch,:)); % selecting channels
    %= Resampling original data to match extracted (fs - if necessary)
    if size(out.mecg,2) ~= size(residual,2)
        % fref:          fetal QRS reference
        % fecgref:       fetal ECG signals (before mixture)
        fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
        fecgref = zeros(length(ch),size(residual,2));
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
    
    %= Getting statistics (exp 2)
    if ~exp3
        [F1,MAE,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
        MAE = MAE*1000/fs_new;
        stats.(method)(origrec,:) = [F1,MAE,PPV,SE]; % dynamic naming
    else %= Getting statistics (exp 3)
        if ~exist([path_orig 'wfdb'],'dir')
            mkdir([path_orig 'wfdb'])
        end
        cd([path_orig 'wfdb'])
        mkdir(num2str(i))
        cd(num2str(i))
        fname = [path_orig 'plots' slashchar fls_ext{i}(1:end-4) cas];
        fname = strcat(fname{:});
        [outputs{1:7}]= morpho_loop(fecgref,residual,fref,fs,TEMP_SAMPS,fname,[b_hp,a_hp,b_lp,a_lp]);
        
        %         end
        morph.(method)(origrec,:) = outputs;
    end
    clear fecg residual fqrs F1 MAE PPV SE qt_err theight_err outputs
    cd ..
end

save([path_orig 'wksp' num2str(filesproc(end))])

%= Plots and statistics generation
LWIDTH = 1.5;
FSIZE = 15;
%== F1 plot
figure
stats_f1 = 100*[stats.JADEICA(:,1) stats.PCA(:,1) stats.tsc(:,1) stats.tspca(:,1) ...
    stats.tsekf(:,1) stats.alms(:,1) stats.arls(:,1) stats.aesn(:,1)];
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
stats_MAE = [stats.JADEICA(:,2) stats.PCA(:,2) stats.tsc(:,2) stats.tspca(:,2) ...
    stats.tsekf(:,2) stats.alms(:,2) stats.arls(:,2) stats.aesn(:,2)];
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
snr00 = cellfun(@(x) ~isempty(regexp(x,'.snr00dB','ONCE')),fls_orig);
snr03 = cellfun(@(x) ~isempty(regexp(x,'.snr03dB','ONCE')),fls_orig);
snr06 = cellfun(@(x) ~isempty(regexp(x,'.snr06dB','ONCE')),fls_orig);
snr09 = cellfun(@(x) ~isempty(regexp(x,'.snr09dB','ONCE')),fls_orig);
snr12 = cellfun(@(x) ~isempty(regexp(x,'.snr12dB','ONCE')),fls_orig);

% Generate Table
counter1 = 1;
table = zeros(16,28);
for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
    eval(['stat = stats.' met{:} ';']);
    % F1
    statscase = 100*[stat(base,1) stat(c0,1) stat(c1,1) stat(c2,1) stat(c3,1) stat(c4,1) stat(c5,1)];
    auxtab = [mean(statscase)',-1.*ones(7,1),std(statscase)',-2.*ones(7,1)];
    auxtab2(counter1,:) = median(statscase)';
    table(counter1,:) = reshape(auxtab',1,7*4);
    counter1 = counter1 + 1;
    
    % MAE
    statscase = [stat(base,2) stat(c0,2) stat(c1,2) stat(c2,2) stat(c3,2) stat(c4,2) stat(c5,2)];
    auxtab = [mean(statscase)',-1.*ones(7,1),std(statscase)',-2.*ones(7,1)];
    table(counter1,:) = reshape(auxtab',1,7*4);
    counter1 = counter1 + 1;
end
table = round(table.*10)./10;
% F1
c=3;
for met = {'JADEICA' 'aesn'}%{'ica' 'pca' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
    figure(c)
    c = c+1;
    eval(['stat = stats.' met{:} ';']);
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

for met = {'JADEICA' 'aesn'}%{'tsekf' 'tspca' 'aesn' 'ica'}
    figure(c)
    c= c+1;
    eval(['stat = stats.' met{:} ';']);
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





