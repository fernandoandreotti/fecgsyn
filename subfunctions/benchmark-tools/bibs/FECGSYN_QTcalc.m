function [qtint,th,qs,tends,tpeak,qrs] = FECGSYN_QTcalc(ann_types,ann_stamp,signal,fs)
% function [qtint,th,qs,tends,tpeak,qrs] = FECGSYN_QTcalc(ann_types,ann_stamp,signal,fs)
% Function that contains heuristics behind QT interval calculation
% Based on assumption that ECGPUWAVE only outputs a wave (p,N,t) if it can
% detect its begin and end. Only highest peak of T-waves marked as biphasic
% are considered for further analysis.
%
%
% Inputs:
% ann_types:          Type of ALL annotations obtained from ECGPUWAVE
% ann_stamp:          Samplestamp of ALL annotations obtained from ECGPUWAVE
% fs:                 Sampling frequency
%
% Outputs:
% qtint:              Length of QT (samples)
% th:                 Height T-wave (no unit)
% qs:                 Q onset locations
% tends:              Locations of T-wave (end)
% twave:              Locations of T-waves (peak)
%
% Examples:
% TODO
%
% See also:
% FECGSYN_benchMorph
% FECGSYN_manalysis
% FECGSYN_morpho_loop
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
ss = round(mean(temp_stamp(Rpeaks+1)-temp_stamp(Rpeaks)));  % S locations
if isempty(Rpeaks) % case signal is really bad
    qs = NaN; 
    tpeak = NaN; 
    qtint = NaN;
    th = NaN;
    tends=NaN;
    qrs = NaN;
    return; 
end
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
    try
    tpeak = round(mean(tpeak-temp(1:length(tpeak))));
    catch
        warning('No peaks found')
    end
end

if numel(Rpeaks) >= 1
    qrs = [qs ss]+temp_stamp(Rpeaks(1));
    qrs = max(signal(qrs(1):qrs(2)))-min(signal(qrs(1):qrs(2)));
else
    qrs = NaN;
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
