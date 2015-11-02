%% Extraction script for FECG morphological analysis
%
% This script performs FECG extractin on the given path using pre-defined
% bandwidths (defined in Experiment 2 and 3 of Andreotti2016)
%
% Input:
% path                  Path where datasets are saved
% narrowband            Bandpass to be used [boolean]. 1 (3-100 Hz) or 0 (0.5-100 Hz)
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-11-2015
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free So                                                           it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function FECGSYN_extract_data(path,narrowband)

%% Parameters
% Channels to be used
ch = [1 8 11 14 19 22 25 32]; % using 8 channels (decided considering Exp. 1)
refchs = 33:34;               % reference channels

% = Defining preprocessing bands (narrow/wide)
if narrowband
    HF_CUT = 100;% high cut frequency
    LF_CUT = 3; % low cut frequency
    [b_lp,a_lp] = butter(5,HF_CUT/(fs_new/2),'low');
    [b_hp,a_hp] = butter(3,LF_CUT/(fs_new/2),'high');
    clear HF_CUT LF_CUT
    path2save = [path 'exp2/'] ;
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
    path2save = [path 'exp3/'] ;
end


%= Read files
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

for i = 1:length(fls)
    disp(['Extracting file ' fls{i} '..'])
    filename = [path2save 'rec' num2str(i)];
    % = loading data
    load(fls{i})
    disp(num2str(i))
    if isempty(out.noise)
        noise = zeros(size(out.mecg));
    else
        noise = sum(cat(3,out.noise{:}),3);
    end
    fs = out.param.fs;
    
    mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
        + noise;     % re-creating abdominal mixture
    mixture = mixture./3000;    % removing gain given during int conversion
    refs = zeros(length(refchs),length(mixture)/(fs/fs_new));
    for j = 1:length(refchs)
        refs(j,:) = resample(mixture(refchs(j),:),fs_new,fs);   % reference maternal channels
    end
    mixture = mixture(ch,:);
    fref = round(out.fqrs{1}/(fs/fs_new));
    mref = round(out.mqrs/(fs/fs_new));
    
    % = preprocessing channels
    ppmixture = zeros(size(mixture,1),size(mixture,2)/(fs/fs_new));
    fecg = ppmixture;
    for j=1:length(ch)
        ppmixture(j,:) = resample(mixture(j,:),fs_new,fs);    % reducing number of channels
        fecg(j,:) = resample(double(out.fecg{1}(j,:))./3000,fs_new,fs);    % propagated FECG signal (unmixed)
        lpmix = filtfilt(b_lp,a_lp,ppmixture(j,:));
        ppmixture(j,:) = filtfilt(b_hp,a_hp,lpmix);
    end
    mixture = ppmixture;
    
    %             % == Extraction
    %-------------------
    %ICA Independent Component Analysis
    %-------------------
    disp('ICA extracthion ..')
    loopsec = 60;   % in seconds
    [residual,~,A] = FECGSYN_bss_extraction(mixture,'JADEICA',fs_new,loopsec,1);     % extract using IC
    [fqrs,maxch] = FECGSYN_QRSmincompare(residual,fref,fs_new);     % detect QRS and channel with highest F1
    %== saving results
    save([filename '_JADEICA'],'maxch','residual','fqrs','fref','A')
    clear fqrs residual maxch
    
    % -------------------
    % PCA Principal Component Analysis
    % -------------------
    disp('PCA extraction ..')
    [residual,~,A] = FECGSYN_bss_extraction(mixture,'PCA',fs_new,loopsec,0);     % extract using IC
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
        residual(j,:) = FECGSYN_ts_extraction(mref,mixture(j,:),'TS-CERUTTI',0,...
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
        residual(j,:) = FECGSYN_ts_extraction(mref,mixture(j,:),'TS-PCA',0,...
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
        residual(j,:) = FECGSYN_kf_extraction(mref,mixture(j,:),NbCycles,fs_new);
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
    mirrow = 30*fs_new;    % mirrow 30 seconds of signal to train method
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
