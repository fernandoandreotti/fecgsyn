function [qt,theight] = FECGSYN_manalysis(abdm_temp,ref_temp,fs,debug)
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

% == generate signal from template (necessary for ecgpuwave)

% resampling and repeating templates
fsnew = 500;        % upsampling to 500Hz so that foetal 
                    % heart looks like adult
wsign = abs(max(abdm_temp))>abs(min(abdm_temp));
wsign = 2*wsign - 1;
abdm_temp = 1000*wsign*abdm_temp/max(abs(abdm_temp)); % normalizing for 
ref_temp = 1000*wsign*ref_temp/max(abs(ref_temp));    % comparing T-height
abdm_temp = resample(abdm_temp,fsnew,fs);
ref_temp = resample(ref_temp,fsnew,fs);
abdm_sig = repmat(abdm_temp,1,20)';
ref_sig = repmat(ref_temp,1,20)';

% == getting annotations right
qrsref = round((0.5 - 1/6)*length(ref_temp));
qrsabdm = round((0.5 - 1/6)*length(abdm_temp));
qrsref = arrayfun(@(x) qrsref + x*length(ref_temp),0:19)';
qrsabdm = arrayfun(@(x) qrsabdm + x*length(abdm_temp),0:19)';

% writting to WFDB
tm1 = 1:length(abdm_sig); tm1 = tm1'-1;
tm2 = 1:length(ref_sig); tm2 = tm2'-1;
wrsamp(tm1,abdm_sig,'absig','','','')
wrsamp(tm2,ref_sig,'refsig','','','')
wrann('absig','qrs',qrsabdm,repmat('N',20,1));
wrann('refsig','qrs',qrsref,repmat('N',20,1));

% == Segmentation using ECGPUWAVE
% ref signal
ecgpuwave('refsig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allref,alltypes_r] = rdann('refsig','edr');
if debug
    close all
    figure(1)
    ax(1)=subplot(2,1,1);
    plot(ref_sig)
    hold on
    plot(allref,ref_sig(allref),'or')
    text(allref,ref_sig(allref)+10,alltypes_r)
end
% test signal
ecgpuwave('absig','edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[alltest,alltypes_t] = rdann('absig','edr');
if debug
    figure(1)
    ax(2)=subplot(2,1,2);
    plot(abdm_sig)
    hold on
    plot(alltest,abdm_sig(alltest),'or')
    text(alltest,abdm_sig(alltest)+50,alltypes_t)
    linkaxes(ax,'x')
end

% == Calculate error on morphological analysis made by extracted data

%% == QT-intervals from ref
% Q
rees = arrayfun(@(x) strcmp(x,'N'),alltypes_r);
obrackts = arrayfun(@(x) strcmp(x,'('),alltypes_r);
idxr = find(rees);
idxqomplete = obrackts(idxr-1);     % finding QRS complexes with begin/end
idxincomp = idxr(~idxqomplete);     % R-peak location of incomplete complexes
quus = allref(idxr(idxqomplete)-1);

cleants = ones(size(rees)); 
if ~isempty(idxincomp)  % throw some beats away
    % throw T-waves away if there is no Q
    for i = 1:length(idxincomp)
       idx = find(idxr == idxincomp(i));
       cleants(idxr(idx):idxr(idx+1)) = 0;        
    end
end

% T
tees = arrayfun(@(x) strcmp(x,'t'),alltypes_r);
cbrackts = arrayfun(@(x) strcmp(x,')'),alltypes_r);
tees = tees&cleants;            % removing T's without Q's

% treating T-waves detected as biphasic
biphasic = filter([1 1],1,tees);    
idxbi = biphasic==2; idxbi = circshift(idxbi,1);
tees_all = tees;    % saving for theight analysis
tees(idxbi) = 0;    % only considering second T
% looking for T ends
idxcbrackt = find(tees)+1;
idxcbrackt = idxcbrackt(cbrackts(idxcbrackt)); % which c-brackts come right after T's
tends = allref(idxcbrackt); % T-end locations

% test if QT analysis feasible
if isempty(tends)
    theight = NaN;
    qt = NaN;
    return
end
qtref = mean(tends-quus)*1000/fsnew/2;    % in ms
% looking for T height
if sum(idxbi) > 0
    twave = allref(find(tees_all,2))-length(ref_temp);
    [~,idx] = max(abs(ref_temp(twave)));
    twave = twave(idx);
else
    twave = allref(find(tees_all,1))-length(ref_temp);
end
thref = abs(ref_temp(twave));
% % % isoeletric line
% % waves = find(allref<twave+length(ref_temp));
% % tbeg = allref(waves(end)) -length(ref_temp);
% % speak = allref(waves(end-1))-length(ref_temp);
%% == QT-intervals from test
% Q
rees = arrayfun(@(x) strcmp(x,'N'),alltypes_t);
obrackts = arrayfun(@(x) strcmp(x,'('),alltypes_t);
idxr = find(rees);
idxqomplete = obrackts(idxr-1);     % finding QRS complexes with begin/end
idxincomp = idxr(~idxqomplete);     % R-peak location of incomplete complexes
quus = alltest(idxr(idxqomplete)-1);

cleants = ones(size(rees)); 
if ~isempty(idxincomp)  % throw some beats away
    % throw T-waves away if there is no Q
    for i = 1:length(idxincomp)
       idx = find(idxr == idxincomp(i));
       cleants(idxr(idx):idxr(idx+1)) = 0;        
    end
end

% T
tees = arrayfun(@(x) strcmp(x,'t'),alltypes_t);
cbrackts = arrayfun(@(x) strcmp(x,')'),alltypes_t);
tees = tees&cleants;            % removing T's without Q's

% treating T-waves detected as biphasic
biphasic = filter([1 1],1,tees);    
idxbi = biphasic==2; idxbi = circshift(idxbi,1);
tees_all = tees;    % saving for theight analysis
tees(idxbi) = 0;    % only considering second T
% looking for T ends
idxcbrackt = find(tees)+1;
idxcbrackt = idxcbrackt(cbrackts(idxcbrackt)); % which c-brackts come right after T's
tends = alltest(idxcbrackt); % T-end locations

% test if QT analysis feasible
if isempty(tends)
    theight = NaN;
    qt = NaN;
    return
end

qttest = mean(tends-quus)*1000/fsnew/2;   % in ms
% looking for T height
if sum(idxbi) > 0
    twave = alltest(find(tees_all,2))-length(abdm_temp);
    [~,idx] = max(abs(abdm_temp(twave)));
    twave = twave(idx);
else
    twave = alltest(find(tees_all,1))-length(abdm_temp);
end
thtest = abs(abdm_temp(twave));
% % % isoeletric line
% % waves = find(allref<twave+length(ref_temp));
% % tbeg = allref(waves(end)) -length(ref_temp);
% % speak = allref(waves(end-1))-length(ref_temp);

%% == QT error
qt = qttest - qtref;        % absolute error in ms

%% == T-height estimation
theight = thtest/thref;

