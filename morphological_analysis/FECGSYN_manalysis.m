function [qt_test,qt_ref,th_test,th_ref,qt_err,th_err] = FECGSYN_manalysis(abdm_temp,ref_temp,qrsabdm,qrsref,fs)
% This function calculates morphological features form signals given two
% templates (reference and abdm). Statistics are give as %.
%
% Input:
%  abdm_temp:               Template to be tested
%  ref_temp:                Reference template
%  qrs_abdm/qrs_ref:        Location of qrs in template
%  fs:                      Sampling frequency
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
qrsref = cumsum([0 repmat(T_LEN,1,19)])+2*qrsref;
qrsabdm = cumsum([0 repmat(T_LEN,1,19)])+2*qrsabdm;
qrsref([1,20]) = []; qrsabdm([1,20]) = [];
% writting to WFDB
tm1 = 1:length(abdm_sig); tm1 = tm1'-1;
tm2 = 1:length(ref_sig); tm2 = tm2'-1;
wrsamp(tm1,abdm_sig,'absig',FS_ECGPU,gain,'')
wrsamp(tm2,ref_sig,'refsig',FS_ECGPU,gain,'')
wrann('absig','qrs',qrsabdm,repmat('N',20,1));
wrann('refsig','qrs',qrsref,repmat('N',20,1));

%% Segmentation using ECGPUWAVE
% ref signal
ecgpuwave('refsig','edr',[],[],'qrsref'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allref,alltypes_r] = rdann('refsig','edr');
% if debug   
%     figure(2)
%     ax(1)=subplot(2,1,1);
%     cla
%     plot(ref_sig./gain)
%     hold on
%     plot(allref,ref_sig(allref)./gain,'or')
%     text(allref,ref_sig(allref)./gain+0.1,alltypes_r)
%     title('Reference Signal')
% end
% test signal
ecgpuwave('absig','edr',[],[],'qrsabdm'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[alltest,alltypes_t] = rdann('absig','edr');
% if debug
%     figure(2)
%     ax(2)=subplot(2,1,2);
%     cla
%     plot(abdm_sig./gain)
%     hold on
%     plot(alltest,abdm_sig(alltest)./gain,'or')
%     text(alltest,abdm_sig(alltest)./gain+0.2,alltypes_t)
%     linkaxes(ax,'x')
%     title('Test Signal')
% end

% == Calculate error on morphological analysis made by extracted data

%% QT-intervals from ref

[qt_ref,th_ref,qs,tends,tpeak] = QTcalc(alltypes_r,allref,ref_sig,fs);
% test if QT analysis feasible
if isnan(qt_ref)||isnan(th_ref)
    qt_test = NaN;
    qt_ref = NaN;
    th_test = NaN;
    th_ref = NaN;
    qt_err = NaN;
    th_err = NaN;   
    disp('manalysis: Could not encounter QT wave for REFERENCE.')
    return
end
if isnan(tpeak)
    disp
end

isoel = median(ref_temp(round(qrsref+0.185*fs):end));
qt_ref = qt_ref*1000/(2*FS_ECGPU);          % in ms
th_ref = (th_ref-isoel)./gain;              % in mV (or not)

if debug
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])   
    ax(1)=subplot(2,1,1);
    cla
    plot(ref_temp./gain,'k','LineWidth',2)
    hold on
    rpeak = qrsref(1)-T_LEN;
    plot(rpeak+qs,ref_temp(rpeak+qs)./gain,'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(rpeak+tends,ref_temp(rpeak+tends)./gain,'ms','MarkerSize',10,'MarkerFaceColor','m')
    try
        plot(rpeak+tpeak,ref_temp(rpeak+tpeak)./gain,'go','MarkerSize',10,'MarkerFaceColor','g')
    catch
    disp    
    end
    
    title('Reference Signal')   
    clear qs tends twave offset rpeak
end

%% QT-intervals from test

[qt_test,th_test,qs,tends,tpeak] = QTcalc(alltypes_t,alltest,abdm_sig,fs);
% test if QT analysis feasible
if isnan(qt_test)||isnan(th_test)
    qt_test = NaN;
    qt_ref = NaN;
    th_test = NaN;
    th_ref = NaN;
    qt_err = NaN;
    th_err = NaN;   
    disp('manalysis: Could not encounter QT wave for TEST.')
    return
end

isoel = median(abdm_temp(round(qrsabdm+0.185*fs):end));
qt_test = qt_test*1000/(2*FS_ECGPU);          % in ms
th_test = (th_test-isoel)./gain;                  % in mV (or not)

if debug   
    figure(1)
    ax(2)=subplot(2,1,2);
    cla
    plot(abdm_temp./gain,'k','LineWidth',2)
    hold on
    try
    rpeak = qrsref(1)-T_LEN;
    plot(rpeak+qs,abdm_temp(rpeak+qs)./gain,'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(rpeak+tends,abdm_temp(rpeak+tends)./gain,'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(rpeak+tpeak,abdm_temp(rpeak+tpeak)./gain,'go','MarkerSize',10,'MarkerFaceColor','g')
    catch
        disp('Failed to plot, results do not make sense!!!')
    end
    clear qs tend twave offset rpeak
end

%% Final results
%== QT error
qt_err = qt_test - qt_ref;        % absolute error in ms
%== T-height estimation
th_err = abs(th_test/th_ref);     % only considering abs value

end


function [qtint,th,qs,tends,tpeak] = QTcalc(ann_types,ann_stamp,signal,fs)
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
rees = arrayfun(@(x) strcmp(x,'N'),temp_types);          % 'R'
goodR = ismember(temp_stamp(rees),validR);
Rpeaks = find(rees);   % annotation number
Rpeaks(~goodR) = [];   % eliminating R's without T
qs = round(mean(temp_stamp(Rpeaks-1)-temp_stamp(Rpeaks)));  % Q locations

clear goodR rees

% == T-wave (end)
% Defined as closing parenthesis after T-wave peak
tees = arrayfun(@(x) strcmp(x,'t'),temp_types);
goodT = ismember(temp_stamp(tees),validT);
Tpeaks = find(tees);   % annotation number
Tpeaks(~goodT) = [];   % eliminating R's without T
tends = round(mean(temp_stamp(Tpeaks+1) - temp_stamp(Rpeaks)));

qtint = tends-qs; % qt interval in samples

% == T-peak
if sum(idxbi) > 0   % case biphasic waves occured
    posmax = [idxbi' idxbi'+1];
    [valbi,bindx]=max(abs(signal(ann_stamp(posmax)))'); % max abs value between tt
    valnonbi = abs(signal(ann_stamp(nonbi))');
    th = mean([valbi valnonbi]); % theight with gain
    for i = 1:length(bindx); tpeak(i)=ann_stamp(posmax(i,bindx(i)));end
    tpeak = sort([tpeak ann_stamp(nonbi)'])';
    tpeak = round(mean(tpeak-temp_stamp(Rpeaks)));
else
    th = mean(abs(signal(ann_stamp(Tpeaks))));   
    tpeak = ann_stamp(Tpeaks);
    tpeak = round(mean(tpeak-temp_stamp(Rpeaks)));
end


if isempty(tends), tends = NaN; end
if qtint>0.5*fs, qtint = NaN; end   % limit QT interval length: MIGHT VARY IN OTHER APPLICATIONS!
if isempty(qs), qs = NaN; end
if isempty(tpeak), tpeak = NaN; end

% == isoeletric line
% % waves = find(ann_stamp<twave+length(ref_temp));
% % tbeg = ann_stamp(waves(end)) -length(ref_temp);
% % speak = ann_stamp(waves(end-1))-length(ref_temp);
end
