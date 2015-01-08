%% Main script for FECG morphological analysis
%
% This is the mains script for testing the morphological consistency of
% extracting the foetal signal using various methods.
% Used extraction methods:
%  - ICA
%  - PCA
%  - TS ...
%
% Used morphological measures:
% - T/QRS ratio
% - ST segment
% - QT interval
%
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
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

global debug fref filename channel
debug = 1;
%% Input parameters
% importing path
% clear all; close all; clc;
slashchar = char('/'*isunix + '\'*(~isunix));
[status,result] = system('hostname');
if status~=0
    error(result);
end
result = strtrim(result);

path = '/media/fernando/FetalEKG/tuning stuff/training/';
path2save = '/media/fernando/FetalEKG/tuning stuff/training/extracted3Hz/';

cd(path)

fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

%% Set-up parameters
generate = 0;   % boolean, data should be generated?
extract = 1;
fs_new = 250;       % signals will be resample to 250 Hz

%% Experiment 2 (later 2 and 3)
% Channels to be used
ch = [1 8 11 22 25 32]; % using 6 channels (decided considering Exp. 1)
refchs = 33:34;

if extract
    for i = 1:10%length(fls)
        tic
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
        INTERV = round(0.05*fs_new);    % BxB acceptance interval
        TH = 0.3;                   % detector threshold
        REFRAC = .15;               % detector refractory period (in s)
        mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
            + noise;     % re-creating abdominal mixture
        mixture = mixture./3000;    % removing gain given during int conversion
        refs = zeros(length(refchs),length(mixture)/(fs/fs_new));
        for j = 1:length(refchs)
            refs(j,:) = resample(mixture(refchs(j),:),fs_new,fs);   % reference maternal channels
        end
        mixture = mixture(ch,:);
        fref = round(out.fqrs{1}/(fs/fs_new));        
        out.mqrs = round(out.mqrs/(fs/fs_new));
        
        % = preprocessing channels
        HF_CUT = 100; % high cut frequency
        LF_CUT = 3; % low cut frequency
        [b_lp,a_lp] = butter(5,HF_CUT/(fs_new/2),'low');
        [b_bas,a_bas] = butter(3,LF_CUT/(fs_new/2),'high');
        ppmixture = zeros(size(mixture,1),size(mixture,2)/(fs/fs_new));
        for j=1:length(ch)
            ppmixture(j,:) = resample(mixture(j,:),fs_new,fs);    % reducing number of channels
            lpmix = filtfilt(b_lp,a_lp,ppmixture(j,:));
            ppmixture(j,:) = filtfilt(b_bas,a_bas,lpmix);
        end
        mixture = ppmixture;
        clear HF_CUT LF_CUT a_bas a_lp b_bas b_lp bw wo lpmix ppmixture
        % == Extraction
        
        % ----------------------------
        % EKF Extended Kalman Filter
        % ----------------------------
        disp('EKF extraction ..')
        NbCycles = 30; % first 30 cycles will be used for template generation
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for channel = 1:length(ch)            
            residual(channel,:) = FECGx_kf_extraction(out.mqrs,mixture(channel,:),NbCycles,fs_new);
            fqrs{channel} = qrs_detect(residual(channel,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
            F1max = 0;
            idxref = (fref>=(min-1)*fs_new*60+1)&(fref<=min*fs_new*60);
            for j = 1:length(ch)
                idx = (fqrs{j}>=(min-1)*fs_new*60+1)&(fqrs{j}<=min*fs_new*60);
                [F1,~,~,~] = Bxb_compare(fref(idxref),fqrs{j}(idx),INTERV);
                if F1 > F1max    % compare and see if this channel provides max F1
                    maxch(min) = j;
                    F1max = F1;
                    fqrs_temp{min} = fqrs{j}(idx);%+ (min-1)*fs_new*60;    % adding fqrs detections to temporary cell
                end
            end
            min = min+1;
        end
        fqrs = cell2mat(fqrs_temp);
        % == saving results
        save([filename '_tsekf'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs NbCycles
        
        
        % ----------------------------
        % TS-PCA
        % ----------------------------
        disp('TS-PCA extraction ..')
        % parameters
        NbCycles = 20;
        NbPC = 2;
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            residual(j,:) = FECGSYN_ts_extraction(out.mqrs,mixture(j,:),'TS-PCA',0,...
                NbCycles,NbPC,fs_new);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
            F1max = 0;
            idxref = (fref>=(min-1)*fs_new*60+1)&(fref<=min*fs_new*60);
            for j = 1:length(ch)
                idx = (fqrs{j}>=(min-1)*fs_new*60+1)&(fqrs{j}<=min*fs_new*60);
                [F1,~,~,~] = Bxb_compare(fref(idxref),fqrs{j}(idx),INTERV);
                if F1 > F1max    % compare and see if this channel provides max F1
                    maxch(min) = j;
                    F1max = F1;
                    fqrs_temp{min} = fqrs{j}(idx);%+ (min-1)*fs_new*60;    % adding fqrs detections to temporary cell
                end
            end
            min = min+1;
        end
        fqrs = cell2mat(fqrs_temp);
        % == saving results
        save([filename '_tspca'],'residual','maxch','fqrs');
        
        clear F1 RMS PPV SE maxch residual fqrs NbCycles NbPC
        
    end
end

