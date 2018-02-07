function sqi = main_sqi_analysis(mref_all,mecg,mixture,residual,fs,varargin)
% Obtain SQI values for window segment
%
% Input:
%  mref_all     - maternal reference QRS locations (in samples)
%  mecg         - maternal reference ECG lead (1xN vecotr)
%  mixture      - abdominal mixture (prior to extraction) (MxN matrix)
%  residual     - extracted fetal ECG signals (MxN matrix)
%  fs           - sampling rate (in Hz)
%  varargin     - optional name given to file, in case several recordings
%                 are being evaluated. Introduces column in sqi table.
%
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
% (SQI indices)
% Andreotti, F., Gräßer, F., Malberg, H., & Zaunseder, S. (2017). Non-Invasive Fetal ECG Signal Quality Assessment for
% Multichannel Heart Rate Estimation. IEEE Trans. Biomed. Eng., (in press).
%
% (FECGSYN Toolbox)
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


if isempty(varargin)
    name = 'unnamed_file';
else
    name = varargin{1};
end


% Parameters
win0 = 5;   % window size for SQI calculation [s]
olap0 = 1;  % overalaping [s]
NDET = 5;   % number of detection algorithms tested
win_alg = 0.05; % window to use around FQRS (in s)
NCHAN = size(mixture,1);
LEN = size(mixture,2);


sqi_table = table;
tline = 1; % table line
win = win0*fs;
olap = olap0*fs;
fqrs_sav = cell(NDET*NCHAN,1);

% Loop through segments
k = 1;
seg = 1;
while k <= LEN-win+1
    fqrs_ann = cell(NCHAN,NDET); % Local detections
    for ch = 1:NCHAN
        % Shift window through data
        res = residual(ch,k:k+win-1); % current segment
        rawsignal = mixture(ch,k:k+win-1);
        sqi_table.filename(tline,1) = {name};
        sqi_table.channel(tline,1) = ch;
        sqi_table.segment(tline,1) = seg;
        sqi_table.tstamp(tline,1) = k+win/2-1; % window center
        if all(res==0) % check if all zeros
            tline = tline+1;
            continue; % skip == 0
        end
        
        %% FQRS detections and references
        % Algorithms denoted by:
        % 1- MaxSearch
        % 2- jqrs
        % 3- PanTompkins
        % 4- gqrs
        % 5- wqrs
        % Performing FQRS detections
        FHR = 2.3; % expected frquency (in Hz)
        fqrs_ann{ch,1} =  OSET_MaxSearch(res,FHR/fs);
        
        fqrs_ann{ch,2} = jqrs_mod(res,0.25,0.3,fs,''); % added refractory period and threshold
        recordName = [name '_seg' num2str(seg) '_ch' num2str(ch)]; % just a temp file
        % WFDB for further detection
        fs_wfdb = 250; gain = 1000;
        res_wfdb = resample(res,2*fs_wfdb,fs);              % resampling to fit adult frequencies at 250 Hz
        [~,I] = max(abs(res_wfdb));  wsign = sign(res_wfdb(I)); % looking for signal sign
        res_wfdb = wsign(1)*gain*res_wfdb/max(abs(res_wfdb)); % fit 2mV
        tm = 1/fs_wfdb:1/fs_wfdb:length(res_wfdb)/fs_wfdb;
        wrsamp(tm',res_wfdb',recordName,fs_wfdb,gain)
        % DanTompkins
        fqrs_ann{ch,3} = pantompkins_qrs(res_wfdb,fs_wfdb);        % making suitable for FECG faster heart rate
        % gqrs from WFDB
        %         if isunix
        %             system(['gqrs -r ' recordName ' -f 0 -s 0']);
        %         else
        gqrs(recordName)
        %         end
        try
            fqrs_ann{ch,4} = rdann(recordName,'qrs');
            fqrs_ann{ch,4} = 2.*fqrs_ann{ch,4}';
        catch
            fqrs_ann{ch,4} = NaN;
        end
        
        % wqrs
        %         if isunix
        %             system(['wqrs -r ' recordName ' -p 50 -R']);
        %         else
        wqrs(recordName)
        %         end
        try
            fqrs_ann{ch,5} = rdann(recordName,'wqrs');
            fqrs_ann{ch,5} = 2.*fqrs_ann{ch,4}';
        catch
            fqrs_ann{ch,5} = NaN;
        end
        fqrs_ann{ch,5} = 2.*fqrs_ann{ch,5}';
        delete([recordName '.*'])
        clear tm res_wfdb I fs_wfdb gain
        %         markers = {'o','x','s','d','^','v','>','<','p','h','+','*','.'};
        %         plot(res)
        %         for c = 1:NDET
        %             hold on
        %             plot(fqrs_ann{ch,c},0.3.*ones(1,length(fqrs_ann{ch,c}))-0.02*c,markers{c})
        %         end
        mref = mref_all(mref_all>k&mref_all<k+win)-k;
        %% Running SQIs
        %= temporal SQI
        sqi_table.ssqi(tline,1) = ssqi(res);                               % Skewness
        sqi_table.ksqi(tline,1) = ksqi(res);                               % Kurtosis
        sqi_table.stdsqi(tline,1) = stdsqi(res);                           % Standard deviation
        %= spectral SQI
        sqi_table.psqi(tline,1) = psqi(res,fs);                            % Peak band power
        sqi_table.bassqi(tline,1) = bassqi(res,fs);                        % Baseline power
        %= detection based SQI
        combs = nchoosek(1:NDET,2);
        for c = 1:size(combs); eval(['sqi_table.bsqi' strcat(num2str(combs(c,:)'))' '(tline,1) = bsqi(fqrs_ann{ch,' ... % Calculate bSQI
                num2str(combs(c,1)) '},fqrs_ann{ch,' num2str(combs(c,2)) '},0.05,fs);']); end
        for c = 1:NDET;eval(['sqi_table.rsqi' num2str(c) '(tline,1) = rsqi(fqrs_ann{ch,' num2str(c) '},fs,0.96);']);end % rSQI metric
        for c = 1:NDET;eval(['sqi_table.csqi' num2str(c) '(tline,1) = csqi(res,fqrs_ann{ch,' num2str(c) '},fs,win_alg);']);end % cSQI
        for c = 1:NDET; eval(['sqi_table.xsqi' num2str(c) '(tline,1) = xsqi(res,fqrs_ann{ch,' num2str(c) '},fs,win_alg);']); end % xSQI
        % Maternal SQIs
        sqi_table.mxsqi(tline,1) = 1- xsqi(res,mref,fs,2*win_alg);         % mxSQI (same as xSQI, broader window, inverse)
        for c = 1:NDET; eval(['sqi_table.misqi' num2str(c) '(tline,1) = 1 - bsqi(mref,fqrs_ann{ch,' num2str(c) '},0.05,fs);']); end; % miSQI metric
        sqi_table.mpsqi(tline,1) = mpsqi1(res,mref,fs);
        sqi_table.mpsqi2(tline,1) = mpsqi2(res,mref,fs);
        sqi_table.mcsqi(tline,1) = mcsqi1(res,mecg(k:k+win-1),fs);          % Spectral Coherence
        sqi_table.mcsqi2(tline,1) = mcsqi2(rawsignal,res,fs);          % Spectral Coherence
        
        tline = tline+1;
        try
            [a, MSGID] = lastwarn();
            warning('off', MSGID)
        catch
        end
    end
    for ch = 1:NCHAN % separated loop because detections are only now available
        for d = 1:NDET
            isqi_temp = isqi(fqrs_ann(:,d),0.05,fs);
            eval(['sqi_table.isqi' num2str(d) '(tline-NCHAN+ch-1,1) = isqi_temp(ch);'])
            % Save detections
            fqrs_sav{ch*d} = [fqrs_sav{ch*d} (fqrs_ann{ch,d}+k-1)];
        end
    end
    k = k+olap;  % move window
    seg = seg+1; % loop table indexes
    disp(['Sample: ' num2str(k)])
end

% saving detected FQRS to data
sqi = sqi_table(:,[1:4 7 5:6 8:19 45:49 20:35 41:44 36:40]); % order according to Andreotti 2017 (Fig. 5)
%  This step is important for using the Naive Bayes classfifier correctly
sqi = sortrows(sqi,'channel','ascend');

