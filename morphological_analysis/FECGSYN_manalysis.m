function [qt,theight] = FECGSYN_manalysis(abdm_temp,ref_temp)
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

% == generate signal from template (necessary for ecgpuwave)
T_LEN = length(abdm_temp);  % template length

% resampling and repeating templates
fs = 250;           % default sampling frequency for ECGPUWAVE
fsnew = 500;        % upsampling to 500Hz so that foetal
% heart looks like adult
gain = 200;        % saving gain for WFDB format
wsign1 = abs(max(abdm_temp))>abs(min(abdm_temp));
wsign1 = 2*wsign1 - 1;
abdm_temp = 1000*wsign1*abdm_temp/max(abs(abdm_temp)); % normalizing for
wsign2 = abs(max(ref_temp))>abs(min(ref_temp));      % comparing T-height
wsign2 = 2*wsign2 - 1;
ref_temp = 1000*wsign2*ref_temp/max(abs(ref_temp));
abdm_sig = repmat(abdm_temp,1,20)';
ref_sig = repmat(ref_temp,1,20)';

% high-passing reference signal
LF_CUT = 0.7;
[b_bas,a_bas] = butter(3,LF_CUT/(fs/2),'high');
ref_sig = filtfilt(b_bas,a_bas,ref_sig);

% == getting annotations right
qrsref = round((0.5 - 1/6)*T_LEN);
qrsabdm = round((0.5 - 1/6)*T_LEN);
qrsref = arrayfun(@(x) qrsref + x*T_LEN,0:19)';
qrsabdm = arrayfun(@(x) qrsabdm + x*T_LEN,0:19)';

% writting to WFDB

tm1 = 1:length(abdm_sig); tm1 = tm1'-1;
tm2 = 1:length(ref_sig); tm2 = tm2'-1;
wrsamp(tm1,abdm_sig,'absig',fs,gain,'')
wrsamp(tm2,ref_sig,'refsig',fs,gain,'')
wrann('absig','qrs',qrsabdm,repmat('N',20,1));
wrann('refsig','qrs',qrsref,repmat('N',20,1));

% == Segmentation using ECGPUWAVE
% ref signal
ecgpuwave('refsig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allref,alltypes_r] = rdann('refsig','edr');
% if debug
%     close all
%     figure(1)
%     ax(1)=subplot(2,1,1);
%     plot(ref_sig)
%     hold on
%     plot(allref,ref_sig(allref),'or')
%     text(allref,ref_sig(allref)+10,alltypes_r)
%     title('Reference Signal')
% end
% test signal
ecgpuwave('absig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[alltest,alltypes_t] = rdann('absig','edr');
% if debug
%     figure(1)
%     ax(2)=subplot(2,1,2);
%     plot(abdm_sig)
%     hold on
%     plot(alltest,abdm_sig(alltest),'or')
%     text(alltest,abdm_sig(alltest)+50,alltypes_t)
%     linkaxes(ax,'x')
%     title('Test Signal')
% end

% == Calculate error on morphological analysis made by extracted data

%% QT-intervals from ref

[qs,tends,twave] = QTcalc(alltypes_r,allref,ref_sig,T_LEN);
% test if QT analysis feasible
if isempty(tends)
    theight = NaN;
    qt = NaN;
    disp('manalysis: Could not encounter QT wave for the template.')
    return
end

offset = sum(qrsref<qs(1))*T_LEN;
thref = abs(ref_temp(twave-offset));



qt_ref = mean(tends-qs)*1000/fsnew;    % in ms

if debug
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    ax(1)=subplot(2,1,1);
    plot(ref_temp,'k','LineWidth',2)
    hold on
    plot(qs(1)-offset,ref_temp(qs(1)-offset),'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(tends(1)-offset,ref_temp(tends(1)-offset),'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(twave-offset,ref_temp(twave-offset),'go','MarkerSize',10,'MarkerFaceColor','g')
    title('Reference Signal')   
end
clear qs tends twave
%% QT-intervals from test

[qs,tends,twave] = QTcalc(alltypes_t,alltest,abdm_sig,T_LEN);
% test if QT analysis feasible
if isempty(tends)
    theight = NaN;
    qt = NaN;
    return
end

offset = sum(qrsref<qs(1))*T_LEN;
thtest = abs(abdm_temp(twave-offset));

if debug   
    figure(1)
    ax(2)=subplot(2,1,2);
    plot(abdm_temp,'k','LineWidth',2)
    hold on
    plot(qs(1)-offset,abdm_temp(qs(1)-offset),'rv','MarkerSize',10,'MarkerFaceColor','r')
    plot(tends(1)-offset,abdm_temp(tends(1)-offset),'ms','MarkerSize',10,'MarkerFaceColor','m')
    plot(twave-offset,abdm_temp(twave-offset),'go','MarkerSize',10,'MarkerFaceColor','g')
    title('Test Signal')
    linkaxes(ax,'x')

end

qt_test = mean(tends-qs)*1000/fsnew;   % in ms
clear qs tends twave
%% QT error
qt = qt_test - qt_ref;        % absolute error in ms

%% T-height estimation
theight = thtest/thref;

end


function [qs,tends,twave] = QTcalc(ann_types,ann_stamp,signal,T_LEN)
%% Function that contains heuristics behind QT interval calculation
% > Inputs
% ann_types:          Type of ALL annotations obtained from ECGPUWAVE
% ann_stamp:          Samplestamp of ALL annotations obtained from ECGPUWAVE
% T_LEN:              Length of template
% 
% > Outputs
% qs:                 Q onset locations
% tends:              Locations of T-wave (end)
% twave:              Locations of T-waves (peak)

% == Q wave
% is defined as an open bracket before the R-peak (no annotation between)
rees = arrayfun(@(x) strcmp(x,'N'),ann_types);
obrackts = arrayfun(@(x) strcmp(x,'('),ann_types);
idxr = find(rees);
idxqomplete = obrackts(idxr-1);     % finding QRS complexes with begin/end
idxincomp = idxr(~idxqomplete);     % R-peak location of incomplete complexes
qs = ann_stamp(idxr(idxqomplete)-1);

% throw some beats away
% throw T-waves away if there is no Q
cleanqs = ones(size(rees));
if ~isempty(idxincomp)
    for i = 1:length(idxincomp)
        idx = find(idxr == idxincomp(i));
        cleanqs(idxr(idx):idxr(idx+1)) = 0;
    end
end

% == T-wave (end)
% Defined as closing parenthesis after T-wave peak
tees = arrayfun(@(x) strcmp(x,'t'),ann_types);
cbrackts = arrayfun(@(x) strcmp(x,')'),ann_types);
tees = tees&cleanqs;            % ignoring T's without Q's

% treating T-waves detected as biphasic
biphasic = filter([1 1],1,tees);
idxbi = biphasic==2; idxbi = circshift(idxbi,-1);
tees_all = tees;    % saving for theight analysis
tees(idxbi) = 0;    % only considering latter T annotation
% no2tees = tees_all;
% no2tees(idxbi|circshift(idxbi,1)) = 0;

% looking for T ends
idxcbrackt = find(tees)+1;
idxcbrackt = idxcbrackt(cbrackts(idxcbrackt)); % which c-brackts come right after T's
tends = ann_stamp(idxcbrackt); % T-end locations

% == T-height
if sum(idxbi) > 0
    twave = find(idxbi&tees_all,1);
    twave = [twave twave+1];
    twave = ann_stamp(twave);
    [~,idx] = max(abs(signal(twave)));
    twave = twave(idx);
else
    twave = ann_stamp(tees_all);
    csum = cumsum([0 ; T_LEN*ones(length(twave)-1,1)],1); % removing shift between beats
    twave = mean(twave-csum);
end

% % % isoeletric line
% % waves = find(ann_stamp<twave+length(ref_temp));
% % tbeg = ann_stamp(waves(end)) -length(ref_temp);
% % speak = ann_stamp(waves(end-1))-length(ref_temp);
end
