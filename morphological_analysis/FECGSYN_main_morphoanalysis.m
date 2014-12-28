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

global debug
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
switch result
    case 'ENGS-19715' % Joachim
        path = '/netshares/ipmprojects3/JB_Experimental_Data/2014.07_fecgsyn_simulations(3.0)/';
        path2functions = '/local/shil3432/Dropbox/DPhil_My_reading_list/fecgsyn-Morph_Analysis/';
        path2save = '/netshares/ipmprojects3/JB_Experimental_Data/out_ICA/';
    otherwise % Fernando loves having to flip slashes
        if isunix
            path = '/home/fernando/tmp/2014.12_fecgsyn_simulations(5.0)/';
            path2save = '/home/fernando/tmp/2014.12_fecgsyn_simulations(5.0)/extracted3Hz/';
        else
            path = 'D:\2014.10_fecgsyn_simulations(5.0)\';
            path2save = 'D:\2014.10_fecgsyn_simulations(5.0)\extracted3Hz\';
        end
end


%% Set-up parameters
generate = 0;   % boolean, data should be generated?
extract = 1;
fs_new = 250;       % signals will be resample to 250 Hz

%% Data Generation
if generate
    mkdir(path)
    FECGSYN_generate_data(path)  % generates set of unique simulated data for testing
else
    cd(path)
end

%% Experiments
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

% % %% Experiment 1
% % % PCA&ICA and number of channels. Testing how the number of channels
% % % available is affecting the results of PCA and ICA.
% % %
% % % == experiments param
% % 
% % ch = {[11 22],[1 11 22 32],[1 8 11 22 25 32],[1 8 11 14 19 22 25 32], ...
% %     [1 3 6 8 9 11 14 16 17 19 22 24 25  27 30 32], ...
% %     1:32}; % trying with 4, 6, 8, 16 and 32 channels
% % 
% % 
% % == core function
% % %NB_REC = 100; % for testing on a few records
% % NB_REC = length(fls);
% % NB_RUN = length(ch)-1; % !! FIXME. NOT RUNNING 32 Channels for now!!
% % stats_struct = cell(NB_RUN,1);
% % cd([path slashchar 'exp1' slashchar])
% % for k = 1:NB_RUN
% %     disp('>>>>>>>>>>>>>>>>>>>>>')
% %     fprintf('processing case with %f channels \n',length(ch{k}));
% %     
% %     for i = 1:NB_REC
% %         for icamethod = {'FASTICA_DEF','FASTICA_SYM','JADEICA'}
% %             disp('==============================');
% %             disp(['Extracting file ' fls{i} '..']);
% %             disp(['The ICA method used is ' icamethod{:}]);
% %             fprintf(' Processing record %f / %f /n (%f abdominal channels) \n',i,NB_REC,length(ch{k}));
% %             
% %             % = loading data
% %             load([path fls{i}])
% %             disp(num2str(i))
% %             if isempty(out.noise)
% %                 noise = zeros(size(out.mecg));
% %             else
% %                 noise = sum(cat(3,out.noise{:}),3);
% %             end
% %             fs = out.param.fs;
% %             INTERV = round(0.05*fs);     % BxB acceptance interval
% %             TH = 0.3;                    % detector threshold
% %             REFRAC = round(.15*fs)/1000; % detector refractory period
% %             mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
% %                 + noise;                 % re-creating abdominal mixture
% %             mixture = mixture(ch{k},:);  % reducing number of channels
% %             
% %             % = preprocessing channels
% %             HF_CUT = 100; % high cut frequency
% %             LF_CUT = 3; % low cut frequency
% %             [b_lp,a_lp] = butter(5,HF_CUT/(fs_new/2),'low');
% %             [b_bas,a_bas] = butter(3,LF_CUT/(fs_new/2),'high');
% %             ppmixture = zeros(size(mixture,1),size(mixture,2)/(fs/fs_new));
% %             for j=1:length(ch{k})
% %                 ppmixture(j,:) = resample(mixture(j,:),fs_new,fs);    % reducing number of channels
% %                 lpmix = filtfilt(b_lp,a_lp,ppmixture(j,:));
% %                 ppmixture(j,:) = filtfilt(b_bas,a_bas,lpmix);
% %                 fref = round(out.fqrs{1}./(fs/fs_new));
% %             end
% %             
% %             % == extraction
% %             % = using ICA (FASTICA or JADE)
% %             disp('ICA extraction ..')
% %             loopsec = 60;   % in seconds
% %             filename = [path2save icamethod{:} '_nbch' num2str(length(ch{k})) '_rec' num2str(i)];
% %             icasig = FECGSYN_bss_extraction(ppmixture,icamethod{:},fs_new,fref,loopsec,filename);     % extract using IC
% %             % Calculate quality measures
% %             qrs = qrs_detect(icasig,TH,REFRAC,fs_new);
% %             if isempty(qrs)
% %                 F1= 0;
% %                 RMS = NaN;
% %                 PPV = 0;
% %                 SE = 0;
% %             else
% %                 [F1,RMS,PPV,SE] = Bxb_compare(fref,qrs,INTERV);
% %             end
% %             eval(['stats_' icamethod{:} '(i,:) = [F1,RMS,PPV,SE];'])            
% %         end
% %         % = using PCA
% %         disp('PCA extraction ..')
% %         loopsec = 60;   % in seconds
% %         filename = [path2save 'PCA_nbch' num2str(length(ch{k})) '_rec' num2str(i)];
% %         icasig = FECGSYN_bss_extraction(ppmixture,'PCA',fs_new,fref,loopsec,filename);     % extract using IC
% %         % Calculate quality measures
% %         qrs = qrs_detect(icasig,TH,REFRAC,fs_new);
% %         if isempty(qrs)
% %             F1= 0;
% %             RMS = NaN;
% %             PPV = 0;
% %             SE = 0;
% %         else
% %             [F1,RMS,PPV,SE] = Bxb_compare(fref,qrs,INTERV);
% %         end
% %         stats_pca(i,:) = [F1,RMS,PPV,SE];
% %     end
% %     stats_struct{k}.stats_pca = stats_pca;
% %     stats_struct{k}.stats_FASTICA_DEF = stats_FASTICA_DEF;
% %     stats_struct{k}.stats_FASTICA_SYM = stats_FASTICA_SYM;
% %     stats_struct{k}.stats_JADEICA = stats_JADEICA;
% %     
% %     save(['stats_ch_' num2str(k)],'stats_struct');
% % end
% % 
% % 
% % % == statistics
% % mean_ica = zeros(NB_RUN,1);
% % median_ica = zeros(NB_RUN,1);
% % mean_pca = zeros(NB_RUN,1);
% % median_pca = zeros(NB_RUN,1);
% % for kk=1:NB_RUN
% %     mean_FASTICA_DEF(kk) = mean(stats_struct{kk}.stats_FASTICA_DEF(1:NB_REC,1));
% %     median_FASTICA_DEF(kk) = median(stats_struct{kk}.stats_FASTICA_DEF(1:NB_REC,1));
% %         mean_FASTICA_SYM(kk) = mean(stats_struct{kk}.stats_FASTICA_SYM(1:NB_REC,1));
% %     median_FASTICA_SYM(kk) = median(stats_struct{kk}.stats_FASTICA_SYM(1:NB_REC,1));
% %         mean_JADEICA(kk) = mean(stats_struct{kk}.stats_JADEICA(1:NB_REC,1));
% %     median_JADEICA(kk) = median(stats_struct{kk}.stats_JADEICA(1:NB_REC,1));
% %     mean_pca(kk) = mean(stats_struct{kk}.stats_pca(1:NB_REC,1));
% %     median_pca(kk) = median(stats_struct{kk}.stats_pca(1:NB_REC,1));
% % end
% 
% % % save(['workspace_exp1_', icamethod]); % save the workspace for history

%% Experiment 2 (later 2 and 3)
% Channels to be used
ch = [1 8 11 22 25 32]; % using 6 channels (decided considering Exp. 1)
refchs = 33:34;

if extract
    for i = 1:length(fls)
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
        out.fqrs{1} = round(out.fqrs{1}/(fs/fs_new));
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
        
        %-------------------
        %ICA Independent Component Analysis
        %-------------------
        disp('ICA extraction ..')
        loopsec = 60;   % in seconds
        icasig = FECGSYN_bss_extraction(mixture,'JADEICA',fs_new,out.fqrs{1},loopsec,filename);     % extract using IC
        % Calculate quality measures
        fqrs = qrs_detect(icasig,TH,REFRAC,fs_new);
        %== saving results
        load([filename '_JADEICA'])
        save([filename '_JADEICA'],'maxch','outdata','fqrs')
        clear fqrs icasig F1 RMS PPV SE outdata
        
        % -------------------
        % PCA Principal Component Analysis
        % -------------------
        disp('PCA extraction ..')
        pcasig = FECGSYN_bss_extraction(mixture,'PCA',fs_new,out.fqrs{1},loopsec,filename);     % extract using IC
        % Calculate quality measures
        fqrs = qrs_detect(pcasig,TH,REFRAC,fs_new);
        % == saving results
        load([filename '_PCA'])
        save([filename '_PCA'],'maxch','outdata','fqrs')
        clear fqrs pcasig qrs F1 RMS PPV SE loopsec outdata
        
        
        % -------------------
        % TS-CERUTTI
        % -------------------
        disp('TS-CERUTTI extraction ..')
        % parameters
        NbCycles = 20;
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            residual(j,:) = FECGSYN_ts_extraction(out.mqrs,mixture(j,:),'TS-CERUTTI',0,...
                NbCycles,'',fs_new);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
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
        % == saving results
        save([filename '_tsc'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs
        
        % -------------------
        % TS-PCA
        % -------------------
        disp('TS-PCA extraction ..')
        % parameters
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
        % == saving results
        save([filename '_tspca'],'residual','maxch','fqrs');
        
        clear F1 RMS PPV SE maxch residual fqrs NbCycles NbPC
        
        % ----------------------------
        % EKF Extended Kalman Filter
        % ----------------------------
        disp('EKF extraction ..')
        NbCycles = 30; % first 30 cycles will be used for template generation
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            residual(j,:) = FECGx_kf_extraction(out.mqrs,mixture(j,:),NbCycles,fs_new);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
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
        % == saving results
        save([filename '_tsekf'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs NbCycles
        
        % ----------------------
        % LMS Least Mean Square
        % ----------------------
        disp('LMS extraction ..')
        %parameters
        refch = 1;      % pick reference channel
        mirrow = 30*fs_new;    % mirrow 30 seconds of signal to train method
        % channel loop
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)], ...
                [refs(refch,mirrow:-1:1) refs(refch,:)],'LMS',debug,fs_new);
            residual(j,:) = res(mirrow+1:end);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
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
        % == saving results
        save([filename '_alms'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs lmsStruct
        
        % ----------------------
        % RLS Recursive Least Square
        % ----------------------
        disp('RLS extraction ..')
        % channel loop
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)],...
                [refs(refch,mirrow:-1:1) refs(refch,:)],'RLS',debug,fs_new);
            residual(j,:) = res(mirrow+1:end);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
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
        % == saving results
        save([filename '_arls'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs rlsStruct
        
        % ----------------------
        % ESN Echo State Neural Network
        % ----------------------
        disp('ESN extraction ..')
        % channel loop
        residual = zeros(size(mixture));
        fqrs = cell(1,size(mixture,1));
        for j = 1:length(ch)
            res = FECGSYN_adaptfilt_extraction([mixture(j,mirrow:-1:1) mixture(j,:)]...
                ,[refs(refch,mirrow:-1:1) refs(refch,:)],'ESN',debug,fs_new);
            residual(j,:) = res(mirrow+1:end);
            fqrs{j} = qrs_detect(residual(j,:),TH,REFRAC,fs_new);
        end
        
        % creating statistics in 1-min blocks
        min = 1;
        maxch = zeros(1,length(mixture)/fs_new/60);
        fqrs_temp = cell(1,length(mixture)/fs_new/60);
        while min <= length(mixture)/fs_new/60;
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
        % == saving results
        save([filename '_aesn'],'residual','maxch','fqrs');
        clear F1 RMS PPV SE maxch residual fqrs ESNparam
        toc
    end
end

%% Generate Results
FECGSYN_genresults(path,path2save,fs_new,ch)

