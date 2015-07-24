function [qt_ref,qt_test,th_ref,th_test,qt_err,th_err] = FECGSYN_manalysis(abdm_temp,ref_temp,qrsabdm,qrsref,fs,filterc,filen)
% This function calculates morphological features form signals given two
% templates (reference and abdm). Statistics are give as %.
%
% Input:
%  abdm_temp:               Template to be tested
%  ref_temp:                Reference template
%  qrs_abdm/qrs_ref:        Location of qrs in template
%  fs:                      Sampling frequency
%  filterc:                 Filter coefficients [b_hp,a_hp,b_lp,a_lp] being
%                           highpass (hp) and lowpass (lp)
%  filen:                   number to be added to ecgpuwaves outputs
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
gain = 1000;        % saving gain for WFDB format

%% Preprocessing
b_hp = filterc(1); a_hp = filterc(2); b_lp = filterc(3);a_lp= filterc(4);
% resample case input data not compatible with ECGPUWAVE
% upsampling to 500Hz so that foetal heart looks like an adult heart
if isnan(qrsabdm)||isnan(qrsref)
    qt_test = NaN;
    qt_ref = NaN;
    th_test = NaN;
    th_ref = NaN;
    qt_err = NaN;
    th_err = NaN;   
    disp('manalysis: could not generate a TEMPLATE.')
    return % does not extract from test
end

abdm_temp = resample(abdm_temp,2*FS_ECGPU,fs);
ref_temp = resample(ref_temp,2*FS_ECGPU,fs);
qrsref = round(qrsref*2*FS_ECGPU/fs);
qrsabdm = round(qrsabdm*2*FS_ECGPU/fs);
T_LENa = length(abdm_temp);  % template length
T_LENr = length(ref_temp);  % template length

abdm_sig = repmat(abdm_temp',1,20);
ref_sig = repmat(ref_temp',1,20);

% Preprocessing reference channel
ref_sig = filtfilt(b_lp,a_lp,ref_sig);
ref_sig = filtfilt(b_hp,a_hp,ref_sig);
abdm_sig = filtfilt(b_lp,a_lp,abdm_sig); 
abdm_sig = filtfilt(b_hp,a_hp,abdm_sig); 

% Normalizing signal
wsign = sign(max(abdm_sig)); % looking for signal sign
abdm_sig = 2*gain*wsign(1)*abdm_sig/max(abs(abdm_sig)); % normalizing in 2 mV
wsign = sign(max(ref_temp)); % looking for signal sign
ref_sig = 2*gain*wsign(1)*ref_sig/max(abs(ref_sig));


%% Saving data as WFDB
% looking for peaks in temporary signal
qrsref = cumsum([0 repmat(T_LENr,1,19)])+qrsref;
qrsabdm = cumsum([0 repmat(T_LENa,1,19)])+qrsabdm;
qrsref([1,20]) = []; qrsabdm([1,20]) = [];
% writting to WFDB
tm1 = 1:length(abdm_sig); tm1 = tm1'-1;
tm2 = 1:length(ref_sig); tm2 = tm2'-1;
filen = filen(regexp(filen,'rec'):end);
counter = 1; % avoind rewriting file
while exist(['absig_' filen '_' num2str(counter) '.hea'],'file')
    counter = counter + 1;
end
filen = [filen '_' num2str(counter)];
wrsamp(tm1,abdm_sig',['absig_' filen],FS_ECGPU,gain,'')
wrsamp(tm2,ref_sig',['refsig_' filen],FS_ECGPU,gain,'')
wrann(['absig_' filen],'qrs',qrsabdm',repmat('N',20,1));
wrann(['refsig_' filen],'qrs',qrsref',repmat('N',20,1));

%% Segmentation using ECGPUWAVE
% ref signal
ecgpuwave(['refsig_' filen],'ecgpuwave',[],[],'qrsref'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allref,alltypes_r] = rdann(['refsig_' filen],'ecgpuwave');
% if debug   
%     figure(2)
%     ax(1)=subplot(2,1,1);
%     cla
%     plot(ref_sig./gain)
%     hold on
%     plot(allref,ref_sig(allref)./gain,'or')
%     plot(qrsref,0,'sg')
%     text(allref,ref_sig(allref)./gain+0.1,alltypes_r)
%     title('Reference Signal')
% end
% test signal
ecgpuwave(['absig_' filen],'ecgpuwave',[],[],'qrsabdm'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[alltest,alltypes_t] = rdann(['absig_' filen],'ecgpuwave');
% if debug
%     figure(2)
%     ax(2)=subplot(2,1,2);
%     cla
%     plot(abdm_sig./gain)
%     hold on
%     plot(qrsabdm,0,'sg')
%     plot(alltest,abdm_sig(alltest)./gain,'or')
%     text(alltest,abdm_sig(alltest)./gain+0.1,alltypes_t)
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
    return % does not extract from test
end

isoel = median(ref_temp(round(qrsref(1)-T_LENr+0.185*fs):end)./gain);
qt_ref = qt_ref*1000/(2*FS_ECGPU);          % in ms
th_ref = th_ref./gain-isoel;              % in mV (or not)




if debug
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])   
    ax(1)=subplot(2,1,1);
    cla
    ref_temp = ref_sig(1:length(ref_temp));
    plot(ref_temp./gain,'k','LineWidth',2)
    hold on
    rpeak = qrsref(1)-T_LENr;
    plot(rpeak+qs,ref_temp(rpeak+qs)./gain,'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(rpeak+tends,ref_temp(rpeak+tends)./gain,'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(rpeak+tpeak,ref_temp(rpeak+tpeak)./gain,'go','MarkerSize',10,'MarkerFaceColor','g')
   
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
isoel = median(abdm_temp(round(qrsabdm(1)-T_LENa+0.185*fs):end)./gain);
qt_test = qt_test*1000/(2*FS_ECGPU);          % in ms
th_test = th_test./gain-isoel;                  % in mV (or not)

if debug   
    figure(1)
    ax(2)=subplot(2,1,2);
    cla
    abdm_temp = abdm_sig(1:length(abdm_temp));
    plot(abdm_temp./gain,'k','LineWidth',2)
    hold on
    try
    rpeak = qrsref(1)-T_LENa;
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


if ~isnan(qt_err)&&isnan(th_err)
    disp('Hold your horses!')
end

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
% fs:                 Sampling frequency
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
QT_MAX = 0.5; % Maximal QT length (in s)  MAY VARY DEPENDING ON APPLICATION!
QT_MIN = 0.1; % Minimal QT length (in s)  MAY VARY DEPENDING ON APPLICATION!
temp_types = ann_types;     % allows exclusion of unsuitable annotations
temp_stamp = ann_stamp;

%== Treat biphasic T-waves
annstr = strcat({temp_types'});
idxbi=cell2mat(regexp(annstr,'tt')); % biphasic
nonbi1=cell2mat(regexp(annstr,'\(t\)')) +1; % regular
nonbi2= cell2mat(regexp(annstr,'\)t\(')) +1; % weird t
nonbi3= cell2mat(regexp(annstr,'\)t\)')) +1; % weird t2
if ~isempty(nonbi2)||~isempty(nonbi3)
    disp('ecpuwave: might be returning weird waves.')
end
nonbi = sort([nonbi1 nonbi2 nonbi3]);
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

if any(Tpeaks >= length(temp_stamp)) % avoid bugs
    Tpeaks(Tpeaks >= length(temp_stamp)) = [];
    Rpeaks = Rpeaks(1:length(Tpeaks));
end
tends = round(mean(temp_stamp(Tpeaks+1) - temp_stamp(Rpeaks)));

qtint = tends-qs; % qt interval in samples

% == T-peak
if sum(idxbi) > 0   % case biphasic waves occured
    posmax = [idxbi' idxbi'+1];
    [valbi,bindx]=max(abs(signal(ann_stamp(posmax)))'); % max abs value between tt
    valnonbi = abs(signal(ann_stamp(nonbi))');
    th = mean([valbi valnonbi']); % theight with gain
    for i = 1:length(bindx); tpeak(i)=ann_stamp(posmax(i,bindx(i)));end
    tpeak = sort([tpeak ann_stamp(nonbi)'])';
    tpeak = round(mean(tpeak-temp_stamp(Rpeaks)));
else
    Tpeaks(ann_stamp(Tpeaks)>length(signal)) = [];
    th = mean(abs(signal(ann_stamp(Tpeaks))));   
    tpeak = ann_stamp(Tpeaks);
    temp = temp_stamp(Rpeaks);
    tpeak = round(mean(tpeak-temp(1:length(tpeak))));
end


if isempty(tends), tends = NaN; end
if qtint>QT_MAX*fs, qtint = NaN; end  
if qtint<QT_MIN*fs, qtint = NaN; end  

if isempty(qs), qs = NaN; end
if isempty(tpeak), tpeak = NaN; end

% == isoeletric line
% % waves = find(ann_stamp<twave+length(ref_temp));
% % tbeg = ann_stamp(waves(end)) -length(ref_temp);
% % speak = ann_stamp(waves(end-1))-length(ref_temp);
end
