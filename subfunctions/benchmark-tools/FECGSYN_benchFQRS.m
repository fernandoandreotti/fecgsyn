function FECGSYN_benchFQRS(path_orig)
% function FECGSYN_benchFQRS(path_orig,fs,ch,exp3)
%
% this script generates statistics as in Experiment 2 by Andreotti et al 2016,
% namely F1,MAE,PPV and SE of fetal QRS detections
%
%
%  Examples:
% TODO
%
% See also:
% TODO
%
% 
% 
% --
% fecgsyn toolbox, version 1.1, March 2016
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
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


slashchar = char('/'*isunix + '\'*(~isunix));
if debug,  mkdir([path_orig 'plots' slashchar]), end


% == Parameters
fs = 1000;
INTERV = round(0.05*fs); % BxB acceptance interval
TEMP_SAMPS = round(60*fs); % samples used for building templates
fs_new = 250;


% Run through extracted datasets
cd(path_orig)
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
for i = length(fls_ext)
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

save([path_orig 'wksp_exp2'])


%


