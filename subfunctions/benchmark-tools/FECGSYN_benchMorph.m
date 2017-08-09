function morph=FECGSYN_benchMorph(pathdir,ch,debug)
% function FECGSYN_benchFQRS(pathdir,debug)
%
% this script generates statistics as in Experiment 3 by Andreotti et al 2016,
% namely fetal QT and T/QRS morphological features. Methods available are
% hardcoded in the function, make sure to alter it in case you include your
% own method.
%
% Input:
%  pathdir             Root directory where original data is saved. It is
%                   expected that the pathdir contains a subfolder "ext" with
%                   extracted data AND and /ext/index.mat file containing
%                   the mapping between original and extracted data.
%
%  ch              Channels used (WFDB only)
%
% debug             If debug is active (=true) will save additional plots
%                   into a 'pathdir/ext/plots' folder
%
% Output:
%  morph          Structure containing benchmark results for all files
%                 contained in pathdir.
%
%
% Examples:
% TODO
%
% See also:
% FECGSYN_benchMorph_plot
% FECGSYN_benchFQRS
% FEGSYN_main_extract
% FEGSYN_main_extract
%
%
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


slashchar = char('/'*isunix + '\'*(~isunix));
if ~strcmp(pathdir(end),slashchar), pathdir = [pathdir slashchar];end

if debug,  mkdir([pathdir 'ext' slashchar 'plots' slashchar]), end
% == Parameters
fs_new = 250;              % ECGPUWAVE function works at 250 Hz
    
% Find out if *mat or wfdb
cd(pathdir)
fls = dir([pathdir '*.mat']); % looking for .mat (creating index)
fls = arrayfun(@(x)x.name,fls,'UniformOutput',false);
if isempty(fls)
    fls = dir([pathdir '*.hea']);
    fls =  arrayfun(@(x) x.name,fls,'UniformOutput',false);
    remfls = cellfun(@(x) length(strtok(x(end:-1:1),'_')),fls);
    for i = 1:length(fls), fls{i} = fls{i}(1:end-remfls(i)-1);end
    fls = unique(fls); % had to remove duplicates
    if isempty(fls)
        error('No file found')
    else
        wfdb = 1;
    end
else
    wfdb = 0;   % working with math files
end


% Reading extracted information
path_ext = [pathdir 'ext' slashchar];
try
    load([path_ext 'index.mat']);
catch
    error('/ext/index.mat not found!')
end
fls_ext = dir([path_ext '*.mat']); % looking for .mat (creating index) NOT INCLUDING WFDB SUPPPORT AT THIS POINT
fls_ext = arrayfun(@(x)x.name,fls_ext,'UniformOutput',false);
idx = cellfun(@(x) strcmp(x(1:3),'rec'),fls_ext); % ignoring files not begining with 'rec'
fls_ext(~idx) = [];

disp('Designing filters ... may take a little while..')
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

% pre-allocation
stats.JADEICA = zeros(length(fls),4);
stats.PCA = zeros(length(fls),4);
stats.tsc = zeros(length(fls),4);
stats.tspca = zeros(length(fls),4);
stats.tsekf = zeros(length(fls),4);
stats.alms = zeros(length(fls),4);
stats.arls = zeros(length(fls),4);
stats.aesn = zeros(length(fls),4);

% = Runs through list of extracted files
for i = 1:length(fls_ext)
    % for i = randperm(length(fls_ext))
    disp(fls_ext{i})
    fprintf('Data %d out of %d \n',i,length(fls_ext));
    
    %= loading extracted file
    [rec,met] = strtok(fls_ext(i),'_');
    method = met{:}(2:end-4);     % Figuring out which extraction method was used, possibilities are:
    % (JADEICA,PCA,tsc,tspca,tsekf,alms,arls,aesn)
    load([path_ext fls_ext{i}])
    %= loading original signal
    file = index(cellfun(@(x) strcmp(rec,x),index(:,2)),1);
    if wfdb
        out = wfdb2fecgsyn([pathdir file{:}],ch);
    else
        load([pathdir file{:}])     %= loading original file
    end
    cas = regexp(file{:},'_c[0-7]','match'); % find out which case is depicted
    if isempty(cas)
        cas = {'bas'};
    else
        cas = cas{1}(2:end);
    end
    clear file
    
    %= re-mixing original signal (necessary after compression)
    fecg = double(out.fecg{1}); % selecting channels
    %= Resampling original data to match extracted (fs - if necessary)
    if size(out.mecg,2) ~= size(residual,2)
        % fref:          fetal QRS reference
        % fecgref:       fetal ECG signals (before mixture)
        fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
        fecgref = zeros(length(ch),size(residual,2));
        for k = 1:size(fecgref,1)
            fecgref(k,:) = resample(fecg(k,:),fs_new,out.param.fs);
        end
    else
        fecgref = fecg;
        fref = out.fqrs{1};
    end
    %= Getting statistics (exp 3)_orig
    TEMP_SAMPS = min(round(60*fs_new),length(residual)/5); % samples used for building templates
                                                           % either 60sec or the whole dataset's length 
                                                           
        
    fname = [pathdir 'plots' slashchar fls_ext{i}(1:end-4) cellstr(cas)];
    fname = strcat(fname{:});
    [outputs{1:7}]= FECGSYN_morpho_loop(fecgref,residual,fref,fs_new,TEMP_SAMPS,fname,[b_hp,a_hp,b_lp,a_lp],debug);
    
    morph.(method)(str2double(rec{1}(4:end)),:) = outputs;
    clear fecg residual fqrs F1 MAE PPV SE  fecg outdata rec  k
    
end

if debug
    save('stats_log','stats','fls')
    FECGSYN_benchMorph_plot(stats,fls) % CHECKOUT this function for some more information!
end


%


