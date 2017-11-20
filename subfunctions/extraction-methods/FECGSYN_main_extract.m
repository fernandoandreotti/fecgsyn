function FECGSYN_main_extract(path,narrowband,wfdb,ch,refchs,debug)
% function FECGSYN_main_extract(path,narrowband,wfdb)
% Extraction script for FECG morphological analysis
%
% This script performs NIFECG extractin on fecgsyn files in a given path using pre-defined
% bandwidths (defined in Experiment 2 and 3 of Andreotti2016). Results are
% saved in Matlab format.
%
% Input:
% path                  Path where datasets are saved
% narrowband            Bandpass to be used [boolean]. 1 (3-100 Hz) or 0 (0.5-100 Hz)
% wfdb                  Load data in WFDB format? [boolean]
% ch                    Channels to extract. In Andreotti et al 2016 ch =  [1 8 11 14 19 22 25 32]
% refchs                Reference maternal channels for adaptive methods. 
%                       In Andreotti et al 2016 refchs = 33:34; 
% debug                 Boolean to define if plots are shown
%
% Examples:
% exp_datagen1
% main_extract_data(cd,1,1)
%
% See also:
% exp_datagen1
% FECGSYN_kf_extraction
% FECGSYN_adaptfilt_extraction
% FECGSYN_bss_extraction
% FECGSYN_ts_extraction
% 
% --
% fecgsyn toolbox, version 1.2, March 2017
% Released under the GNU General Public License
%
% Copyright (C) 2017  Joachim Behar & Fernando Andreotti
% Department of Engineering Science, University of Oxford
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Behar, J., Andreotti, F., Zaunseder, S., Li, Q., Oster, J., & Clifford, G. D. (2014). An ECG Model for Simulating 
% Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings. Physiol. Meas., 35(8), 1537–1550.
% 
%
% Last updated : 15-03-2017
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
%


%% Parameters
% Channels to be used
fs_new = 250;           % extraction occurs at 250 Hz, data will be resampled, if necessary. 
                        % Only relevant for ECG segmantation using ECGPUWAVE
gain = 3000;            % gain apploed saving files

%% Preprocessing
slashchar = char('/'*isunix + '\'*(~isunix));
spath = [path 'ext' slashchar]; % saving folder
if ~exist(spath,'dir'), mkdir(spath);end
if ~strcmp(path(end),slashchar), path = [path slashchar];end

% = Defining preprocessing bands (narrow/wide)
if narrowband
    HF_CUT = 100;% high cut frequency
    LF_CUT = 3; % low cut frequency
    [b_lp,a_lp] = butter(5,HF_CUT/(fs_new/2),'low');
    [b_hp,a_hp] = butter(3,LF_CUT/(fs_new/2),'high');
    clear HF_CUT LF_CUT
else
    % Preprocessing more carefully
    % high-pass filter
    Fstop = 0.5;  % Stopband Frequency
    Fpass = 1;    % Passband Frequency
    Apass = 0.1;  % Passband Ripple (dB)
    Astop = 20;   % Stopband Attenuation (dB)
    h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, fs_new);
    Hhp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf', ...
        'SystemObject', true);
    [b_hp,a_hp] = tf(Hhp);
    % low-pass filter
    Fpass = 100;   % Passband Frequency
    Fstop = 110;  % Stopband Frequency
    Apass = 0.1;    % Passband Ripple (dB)
    Astop = 20;   % Stopband Attenuation (dB)
    h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, fs_new);
    Hlp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf');
    [b_lp,a_lp] = tf(Hlp);
    clear Fstop Fpass Astop Apass h Hhp Hlp
end


%= Read files
if wfdb
    fls = dir([path '*.hea']);
    fls =  arrayfun(@(x) x.name,fls,'UniformOutput',false);
    remfls = cellfun(@(x) length(strtok(x(end:-1:1),'_')),fls);
    for i = 1:length(fls), fls{i} = fls{i}(1:end-remfls(i)-1);end
    fls = unique(fls); % had to remove duplicates
else    
    fls = dir([path '*.mat']);     % looking for .mat (creating index)
    fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
end
% cd(path)
index = cell(length(fls),2); index(:,1) = fls;
for i = 1:length(fls)
    disp(['Extracting file ' fls{i} '..'])
    filename = [spath 'rec' num2str(i)];
    index(i,2) = {['rec' num2str(i)]};
    disp(num2str(i))
    % = loading data (wfdb or mat)
    if wfdb
        out = wfdb2fecgsyn([path fls{i}],[ch refchs]);
        ch = 1:length(ch);
        refchs = length(ch)+[1:length(refchs)];
        warning('For batch processing, avoid converting data to WFDB format, since it may cause slow downs.')
    else    
        load([path fls{i}]) % out structure
    end
    %% Mixing separate sources
    if isempty(out.noise)
        noise = zeros(size(out.mecg));
    else
        noise = sum(cat(3,out.noise{:}),3);
    end
    fs = out.param.fs;
    
    mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
        + noise;     % re-creating abdominal mixture
    refs = zeros(length(refchs),length(mixture)/(fs/fs_new));
    for j = 1:length(refchs)
        refs(j,:) = resample(mixture(refchs(j),:),fs_new,fs);   % reference maternal channels
    end
    mixture = mixture(ch,:)./gain;
    out.fecg{1} = out.fecg{1}(ch,:)./gain;
    fref = round(out.fqrs{1}/(fs/fs_new));
    mref = round(out.mqrs/(fs/fs_new));
    
    %% Preprocessing channels
    ppmixture = zeros(size(mixture,1),size(mixture,2)/(fs/fs_new));
    fecg = ppmixture;
    for j=1:length(ch)
        ppmixture(j,:) = resample(mixture(j,:),fs_new,fs);    % reducing number of channels
        fecg(j,:) = resample(double(out.fecg{1}(j,:)),fs_new,fs);    % propagated FECG signal (unmixed)
        lpmix = filtfilt(b_lp,a_lp,ppmixture(j,:));
        ppmixture(j,:) = filtfilt(b_hp,a_hp,lpmix);
    end
    mixture = ppmixture;
    
    
    %% Extraction Methods
    %-------------------
    %ICA Independent Component Analysis
    %-------------------
    disp('ICA extracthion ..')
    loopsec = 60;   % in seconds
    [residual,~,A] = FECGSYN_bss_extraction(mixture,'JADEICA',fs_new,loopsec,debug);     %#ok<*ASGLU> % extract using IC
    [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);     % detect QRS and channel with highest F1
    %== saving results
    save([filename '_JADEICA'],'maxch','residual','fqrs','fref','A')
    clear fqrs residual maxch
    
    % -------------------
    % PCA Principal Component Analysis
    % -------------------
    disp('PCA extraction ..')
    [residual,~,A] = FECGSYN_bss_extraction(mixture,'PCA',fs_new,loopsec,debug);     % extract using IC
    [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
    % == saving results
    save([filename '_PCA'],'maxch','residual','fqrs','fref','A')
    clear fqrs residual maxch loopsec
    
    % -------------------
    % TS-CERUTTI
    % -------------------
    disp('TS-CERUTTI extraction ..')
    % parameters
    NbCycles = 20;
    residual = zeros(size(mixture));
    for j = 1:length(ch)
        residual(j,:) = FECGSYN_ts_extraction(mref,mixture(j,:),'TS-CERUTTI',debug,...
            NbCycles,'',fs_new);
    end
    [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
    % == saving results
    save([filename '_tsc'],'residual','maxch','fqrs','fref');
    clear maxch residual fqrs
    
    % -------------------
    % TS-PCA
    % -------------------
    disp('TS-PCA extraction ..')
    % parameters
    NbPC = 2;
    residual = zeros(size(mixture));
    for j = 1:length(ch)
        residual(j,:) = FECGSYN_ts_extraction(mref,mixture(j,:),'TS-PCA',debug,...
            NbCycles,NbPC,fs_new);
    end
    [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
    % == saving results
    save([filename '_tspca'],'residual','maxch','fqrs','fref');
    
    clear maxch residual fqrs NbCycles NbPC
     
     % ----------------------------
     % EKF Extended Kalman Filter
     % ----------------------------
     disp('EKF extraction ..')
     NbCycles = 30; % first 30 cycles will be used for template generation
     residual = zeros(size(mixture));
     for j = 1:length(ch)
       residual(j,:) = FECGSYN_kf_extraction(mref,mixture(j,:),debug,NbCycles,fs_new);
     end
     [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
     save([filename '_tsekf'],'residual','maxch','fqrs','fref');
     % == saving results
     clear maxch residual fqrs NbCyclesmat
     % ----------------------
     % LMS Least Mean Square
     % ----------------------
     disp('LMS extraction ..')
     %parameters
     refch = 1;      % pick reference channel
     mirrow = min(30*fs_new,size(mixture,2));    % mirrow 30 seconds or whole signal (depending on size)
     % channel loop
     residual = zeros(size(mixture));
     for j = 1:length(ch)
         res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)], ...
             [refs(refch,mirrow:-1:1) refs(refch,:)],'LMS',debug,fs_new);
         residual(j,:) = res(mirrow+1:end);
     end
     [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
     % == saving results
     save([filename '_alms'],'residual','maxch','fqrs','fref');
     clear maxch residual fqrs lmsStruct
   
     % ----------------------
     % RLS Recursive Least Square
     % ----------------------
     disp('RLS extraction ..')
     % channel loop
     residual = zeros(size(mixture));
     for j = 1:length(ch)
         res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)],...
             [refs(refch,mirrow:-1:1) refs(refch,:)],'RLS',debug,fs_new);
         residual(j,:) = res(mirrow+1:end);
     end
     [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
     % == saving results
     save([filename '_arls'],'residual','maxch','fqrs','fref');
     clear maxch residual fqrs rlsStruct
     
     % ----------------------
     % ESN Echo State Neural Network
     % ----------------------
     disp('ESN extraction ..')
     % channel loop
     residual = zeros(size(mixture));
     for j = 1:length(ch)
         res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)]...
             ,[refs(refch,mirrow:-1:1) refs(refch,:)],'ESN',debug,fs_new);
         residual(j,:) = res(mirrow+1:end);
     end
     [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);    % detect QRS and channel with highest F1
     % == saving results
     save([filename '_aesn'],'residual','maxch','fqrs','fref');
     clear maxch residual fqrs ESNparam
end
save([spath 'index'],'index')
