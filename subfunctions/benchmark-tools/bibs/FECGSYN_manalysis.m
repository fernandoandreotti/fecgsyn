function [qt_ref,qt_test,tqrs_ref,tqrs_test] = FECGSYN_manalysis(abdm_temp,ref_temp,qrsabdm,qrsref,fs,filterc,filen,debug)
% function [qt_ref,qt_test,tqrs_ref,tqrs_test] = FECGSYN_manalysis(abdm_temp,ref_temp,qrsabdm,qrsref,fs,filterc,filen,debug)
% This function calculates morphological features form signals given two
% templates (reference and abdominal). Statistics are given in percentage.
% This functin makes use of the ECGPUWAVE script (Jane et al 1996) and
% required that the WFDB-App is installed and in Matlab's path.
%
% Input:
%  abdm_temp:               Template to be tested
%  ref_temp:                Reference template
%  qrs_abdm/qrs_ref:        Location of qrs in each template
%  fs:                      Sampling frequency
%  filterc:                 Filter coefficients [b_hp,a_hp,b_lp,a_lp] being
%                           highpass (hp) and lowpass (lp)
%  filen:                   number to be added to ecgpuwaves outputs
%
%
%
% Reference to functions:
%
%  ECGPUWAVE: Jane, R., Blasi, A., Garcia, J., & Laguna, P. (1997). Evaluation of an automatic
%  threshold based detector of waveform limits in Holter ECG with the QT database. In Computers in
%  Cardiology 1997 (pp. 295–298). IEEE. http://doi.org/10.1109/CIC.1997.647889
%  This script is available at http://www.physionet.org/physiotools/ecgpuwave/
%  and
%  in the WFDB-App Toolbox: Silva, I, Moody, G. "An Open-source Toolbox for Analysing and Processing
%  PhysioNet Databases in MATLAB and Octave." Journal of Open Research Software 2(1):e27
%  [http://dx.doi.org/10.5334/jors.bi] ; 2014.
%
% Examples:
% TODO
%
% Examples:
% TODO
%
% See also:
% FECGSYN_benchMorph
% FECGSYN_morpho_loop
% FECGSYN_QTcalc
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

% Case second run of ICA tests (speed up a bit)
if isequal(abdm_temp,ref_temp); ident = 1; else ident = 0; end

% resampling and repeating templates
FS_ECGPU = 250;     % default sampling frequency for ECGPUWAVE
gain = 3000;        % saving gain for WFDB format

%% Preprocessing
b_hp = filterc(1); a_hp = filterc(2); b_lp = filterc(3);a_lp= filterc(4);
% resample case input data not compatible with ECGPUWAVE
% upsampling to 500Hz so that foetal heart looks like an adult heart
if isnan(qrsabdm)||isnan(qrsref)
    qt_test = NaN;
    qt_ref = NaN;
    tqrs_test = NaN;
    tqrs_ref = NaN;
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
[~,I] = max(abs(abdm_sig));
wsign = sign(abdm_sig(I)); % looking for signal sign
abdm_sig = 2*gain*wsign(1)*abdm_sig/max(abs(abdm_sig)); % normalizing in 2 mV
[~,I] = max(abs(ref_temp));
wsign = sign(ref_temp(I)); % looking for signal sign
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


%% Segmentation using ECGPUWAVE
% ref signal
counter = 1; % avoind rewriting file
while exist(['refsig_' filen '_' num2str(counter) '.hea'],'file')
    counter = counter + 1;
end
filen = [filen '_' num2str(counter)];
recordName = ['refsig_' filen];


wrsamp(tm2,ref_sig',recordName,FS_ECGPU,gain,'')
wrann(recordName,'qrs',qrsref',repmat('N',20,1));

% if isunix
%     system(['ecgpuwave -r ' recordName '.hea' ' -a ecgpu -i qrsref | awk ''{print $1}'' > tmpout']);
%     fp = fopen('tmpout','r');
%     tmp = textscan(fp,'%s','delimiter','\n');
%     tmp = tmp{1};
%     tmp = str2double(tmp);
%     tmp = round(tmp * fs)+1; % convert 0-1 indexing
% else

ecgpuwave(recordName,'ecgpu',[],[],'qrsref'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
% end
[allref,alltypes_r] = rdann(recordName,'ecgpu');
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

if ~ident
    wrsamp(tm1,abdm_sig',['absig_' filen],FS_ECGPU,gain,'')
    wrann(['absig_' filen],'qrs',qrsabdm',repmat('N',20,1));
    ecgpuwave(['absig_' filen],'ecgpu',[],[],'qrsabdm'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
    [alltest,alltypes_t] = rdann(['absig_' filen],'ecgpu');
    if isempty(alltest)
        qt_test = NaN;
        qt_ref = NaN;
        tqrs_test = NaN;
        tqrs_ref = NaN;
        tqrs_test = NaN;
        tqrs_ref = NaN;
        disp('manalysis: Could not encounter QT wave for TEST signal.')
        return
    end
end
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

[qt_ref,th_ref,qs,tends,tpeak,qrs_ref] = FECGSYN_QTcalc(alltypes_r,allref,ref_sig,fs);
% test if QT analysis feasible
if isnan(qt_ref)||isnan(th_ref)
    qt_test = NaN;
    qt_ref = NaN;
    tqrs_test = NaN;
    tqrs_ref = NaN;
    tqrs_test = NaN;
    tqrs_ref = NaN;
    disp('manalysis: Could not encounter QT wave for REFERENCE signal.')
    return % does not extract from test
end

isoel = median(ref_temp(round(qrsref(1)-T_LENr+0.185*fs):end)./gain);
qt_ref = qt_ref*1000/(2*FS_ECGPU);          % in ms
th_ref = th_ref./gain-isoel;              % in mV (or not)
qrs_ref = qrs_ref./gain;
tqrs_ref = 100*abs(th_ref/qrs_ref);

if debug
    figure(1)
    clf,cla
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    ax(1)=subplot(2,1,1);
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

if ident
    qt_test = '';
    tqrs_test = '';
    return
end

%% QT-intervals from test
[qt_test,th_test,qs,tends,tpeak,qrs_test] = FECGSYN_QTcalc(alltypes_t,alltest,abdm_sig,fs);
% test if QT analysis feasible
if isnan(qt_test)||isnan(th_test)
    qt_test = NaN;
    qt_ref = NaN;
    tqrs_test = NaN;
    tqrs_ref = NaN;
    tqrs_test = NaN;
    tqrs_ref = NaN;
    disp('manalysis: Could not encounter QT wave for TEST.')
    return
end
isoel = median(abdm_temp(round(qrsabdm(1)-T_LENa+0.185*fs):end)./gain);
qt_test = qt_test*1000/(2*FS_ECGPU);          % in ms
th_test = th_test./gain-isoel;                  % in mV (or not)
qrs_test = qrs_test./gain;
tqrs_test = 100*abs(th_test/qrs_test);

if debug&&~ident
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
        warning('Failed to plot, results do not make sense!!!')
    end
    clear qs tend twave offset rpeak
end
%% Final results
%== QT error
qt_err = qt_test - qt_ref;        % absolute error in ms
%== T-height estimation
th_err = abs(th_test/th_ref);     % only considering abs value
end


