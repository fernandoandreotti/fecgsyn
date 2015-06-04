function [qt_test,qt_ref,th_test,th_ref,qt_err,th_err] = FECGSYN_manalysis(abdm_temp,ref_temp,fs)
% This function calculates morphological features form signals given two
% templates (reference and abdm). Statistics are give as %.
%
% Input:
%  abdm_temp:       Template to be tested
%  path_ext:        Path for extracted dataset
%  fs:              Sampling frequency
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 26-09-2014
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

% resampling and repeating templates
FS_ECGPU = 250;     % default sampling frequency for ECGPUWAVE
gain = 200;        % saving gain for WFDB format

%% Preprocessing

% resample case input data not compatible with ECGPUWAVE
% upsampling to 500Hz so that foetal heart looks like an adult heart
abdm_temp = resample(abdm_temp,2*FS_ECGPU,fs);
ref_temp = resample(ref_temp,2*FS_ECGPU,fs);
T_LEN = length(abdm_temp);  % template length
    
wsign = abs(max(abdm_temp))>abs(min(abdm_temp));
wsign = 2*wsign - 1;
abdm_temp = gain*wsign*abdm_temp/max(abs(abdm_temp)); % normalizing for
ref_temp = gain*wsign*ref_temp/max(abs(ref_temp));
abdm_sig = repmat(abdm_temp,1,20)';
ref_sig = repmat(ref_temp,1,20)';

% Preprocessing reference channel
% high-pass filter
Fstop = 0.5;  % Stopband Frequency
Fpass = 1;    % Passband Frequency
Astop = 20;   % Stopband Attenuation (dB)
Apass = 0.1;  % Passband Ripple (dB)
h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, FS_ECGPU);
Hhp = design(h, 'butter', ...
    'MatchExactly', 'stopband', ...
    'SOSScaleNorm', 'Linf', ...
    'SystemObject', true);
[b_hp,a_hp] = tf(Hhp);
% low-pass filter
Fpass = 100;   % Passband Frequency
Fstop = 110;  % Stopband Frequency
Apass = 1;    % Passband Ripple (dB)
Astop = 20;   % Stopband Attenuation (dB)
h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, FS_ECGPU);
Hlp = design(h, 'butter', ...
    'MatchExactly', 'stopband', ...
    'SOSScaleNorm', 'Linf');
[b_lp,a_lp] = tf(Hlp);
clear Fstop Fpass Astop Apass h Hhp Hlp
ref_sig = filtfilt(b_lp,a_lp,ref_sig);
ref_sig = filtfilt(b_hp,a_hp,ref_sig);

%% Saving data as WFDB
% looking for peaks in temporary signal
[~,qrsref] = findpeaks(ref_sig,'MinPeakDistance',T_LEN-30);
[~,qrsabdm] = findpeaks(abdm_sig,'MinPeakDistance',T_LEN-30);

% writting to WFDB
tm1 = 1:length(abdm_sig); tm1 = tm1'-1;
tm2 = 1:length(ref_sig); tm2 = tm2'-1;
wrsamp(tm1,abdm_sig,'absig',FS_ECGPU,gain,'')
wrsamp(tm2,ref_sig,'refsig',FS_ECGPU,gain,'')
wrann('absig','qrs',qrsabdm,repmat('N',20,1));
wrann('refsig','qrs',qrsref,repmat('N',20,1));

%% Segmentation using ECGPUWAVE
% ref signal
ecgpuwave('refsig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allref,alltypes_r] = rdann('refsig','edr');
% if debug   
%     figure(1)
%     ax(1)=subplot(2,1,1);
%     plot(ref_sig./gain)
%     hold on
%     plot(allref,ref_sig(allref)./gain,'or')
%     text(allref,ref_sig(allref)./gain+0.1,alltypes_r)
%     title('Reference Signal')
% end
% test signal
ecgpuwave('absig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[alltest,alltypes_t] = rdann('absig','edr');
% if debug
%     figure(1)
%     ax(2)=subplot(2,1,2);
%     plot(abdm_sig./gain)
%     hold on
%     plot(alltest,abdm_sig(alltest)./gain,'or')
%     text(alltest,abdm_sig(alltest)./gain+0.2,alltypes_t)
%     linkaxes(ax,'x')
%     title('Test Signal')
% end

% == Calculate error on morphological analysis made by extracted data

%% QT-intervals from ref

[qt_ref,th_ref,qs,tends,tpeak] = QTcalc(alltypes_r,allref,ref_sig,T_LEN);
% test if QT analysis feasible
if isnan(qt_ref)||isnan(th_ref)
    qt_test = NaN;
    qt_ref = NaN;
    th_test = NaN;
    th_ref = NaN;
    qt_err = NaN;
    th_err = NaN;   
    disp('manalysis: Could not encounter QT wave for the template.')
    return
end
isoel = median(ref_temp(round(qrsref+0.185*fs):end));
qt_ref = qt_ref*1000/(2*FS_ECGPU);          % in ms
th_ref = (th_ref-isoel)./gain;              % in mV (or not)

if debug
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])   
    pskip = sum(qrsref<qs(1)); % peaks to be skipped
    offset = pskip*T_LEN;
    ax(1)=subplot(2,1,1);
    cla
    plot(ref_temp./gain,'k','LineWidth',2)
    hold on
    plot(qs(pskip)-offset,ref_temp(qs(pskip)-offset)./gain,'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(tends(pskip)-offset,ref_temp(tends(pskip)-offset)./gain,'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(tpeak(pskip)-offset,ref_temp(tpeak(pskip)-offset)./gain,'go','MarkerSize',10,'MarkerFaceColor','g')
    title('Reference Signal')   
    clear qs tends twave offset
end

%% QT-intervals from test

[qt_test,th_test,qs,tends,tpeak] = QTcalc(alltypes_t,alltest,abdm_sig,T_LEN);
% test if QT analysis feasible
if isnan(qt_test)||isnan(th_test)
    qt_test = NaN;
    qt_ref = NaN;
    th_test = NaN;
    th_ref = NaN;
    qt_err = NaN;
    th_err = NaN;   
    return
end

isoel = median(abdm_temp(round(qrsabdm+0.185*fs):end));
qt_test = qt_test*1000/(2*FS_ECGPU);          % in ms
th_test = (th_test-isoel)./gain;                  % in mV (or not)

if debug   
    pskip = sum(qrsref<qs(1));
    offset = pskip*T_LEN;
    figure(1)
    ax(2)=subplot(2,1,2);
    cla
    plot(abdm_temp./gain,'k','LineWidth',2)
    hold on
    plot(qs(pskip)-offset,abdm_temp(qs(pskip)-offset)./gain,'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(tends(pskip)-offset,abdm_temp(tends(pskip)-offset)./gain,'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(tpeak(pskip)-offset,abdm_temp(tpeak(pskip)-offset)./gain,'go','MarkerSize',10,'MarkerFaceColor','g')
    title('Test Signal')
    linkaxes(ax,'x')
    clear qs tend twave offset
end

%% Final results
%== QT error
qt_err = qt_test - qt_ref;        % absolute error in ms
%== T-height estimation
th_err = abs(th_test/th_ref);     % only considering abs value

end


function [qtint,th,qs,tends,tpeak] = QTcalc(ann_types,ann_stamp,signal,T_LEN)
%% Function that contains heuristics behind QT interval calculation
% Based on assumption that ECGPUWAVE only outputs a wave (p,N,t) if it can 
% detect its begin and end. Only highest peak of T-waves marked as biphasic 
% are considered for further analysis.
% 
% 
% Inputs
% ann_types:          Type of ALL annotations obtained from ECGPUWAVE
% ann_stamp:          Samplestamp of ALL annotations obtained from ECGPUWAVE
% T_LEN:              Length of template
% 
% Outputs
% qtint:              Length of QT (samples)
% th:                 Height T-wave (no unit)
% qs:                 Q onset locations
% tends:              Locations of T-wave (end)
% twave:              Locations of T-waves (peak)
%
%
%

temp_types = ann_types;     % allows exclusion of unsuitable annotations
temp_stamp = ann_stamp;

%== Treat biphasic T-waves
annstr = strcat({temp_types'});
idxbi=cell2mat(regexp(annstr,'tt')); % biphasic
nonbi=cell2mat(regexp(annstr,'\(t\)')) +1; % regular
temp_types(idxbi) = [];    % temporarilly clearing first T in biphasic cases
temp_stamp(idxbi) = [];
clear biphasic tees

%== Disregard R-peaks not followed by T-waves
obrackts = arrayfun(@(x) strcmp(x,'('),temp_types);      % '('
cbrackts = arrayfun(@(x) strcmp(x,')'),temp_types);      % ')'
pees = arrayfun(@(x) strcmp(x,'p'),temp_types);      % 'p'
temp_types2 = temp_types; temp_stamp2 = temp_stamp;
temp_types2(obrackts|cbrackts|pees) = [];
temp_stamp2(obrackts|cbrackts|pees) = [];
annstr = strcat({temp_types2'});
idxR = cell2mat(regexp(annstr,'Nt'));  % looking for 'N's followed by 't's
validR = temp_stamp2(idxR);           % valid R-peak sample-stamps
validT = temp_stamp2(idxR+1);           % valid first T-peak sample-stamps
clear idxR annstr pees obrackts cbrackts temp_stamp2 temp_types2
% == Q wave (start)
% is defined as an open bracket before the R-peak (no annotation between)
rees = arrayfun(@(x) strcmp(x,'N'),temp_types);          % 'R'
goodR = ismember(temp_stamp(rees),validR);
Rpeaks = find(rees);   % annotation number
Rpeaks(~goodR) = [];   % eliminating R's without T
qs = temp_stamp(Rpeaks-1);  % Q locations
clear Rpeaks goodR rees

% == T-wave (end)
% Defined as closing parenthesis after T-wave peak
tees = arrayfun(@(x) strcmp(x,'t'),temp_types);
goodT = ismember(temp_stamp(tees),validT);
Tpeaks = find(tees);   % annotation number
Tpeaks(~goodT) = [];   % eliminating R's without T
tends = temp_stamp(Tpeaks+1); % T ends
% assure that T-ends are within template limits
if isempty(tends)||tends(1)>2*T_LEN, tends = NaN; end
%clear Rpeaks goodR
if isempty(qs), qs = NaN; end
qtint = mean(tends-qs); % qt interval in samples

% == T-peak
if sum(idxbi) > 0   % case biphasic waves occured
    posmax = [idxbi' idxbi'+1];
    [valbi,bindx]=max(abs(signal(ann_stamp(posmax)))'); % max abs value between tt
    valnonbi = abs(signal(ann_stamp(nonbi))');
    th = mean([valbi valnonbi]); % theight with gain
    for i = 1:length(bindx); tpeak(i)=ann_stamp(posmax(i,bindx(i)));end
    tpeak = sort([tpeak ann_stamp(nonbi)'])';
else
    th = mean(abs(signal(ann_stamp(Tpeaks))));   
    tpeak = ann_stamp(Tpeaks);
end
if isempty(tpeak), tpeak = NaN; end

% == isoeletric line
% % waves = find(ann_stamp<twave+length(ref_temp));
% % tbeg = ann_stamp(waves(end)) -length(ref_temp);
% % speak = ann_stamp(waves(end-1))-length(ref_temp);
end
