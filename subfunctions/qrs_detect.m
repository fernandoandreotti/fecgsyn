function [qrs_pos,sign,en_thres] = qrs_detect(ecg,varargin)
% QRS detector based on the P&T method. This is an offline implementation
% of the detector.
%
% inputs
%   ecg:            one ecg channel on which to run the detector (required)
%   varargin
%       THRES:      energy threshold of the detector (default: 0.6)
%       REF_PERIOD: refractory period in sec between two R-peaks (default: 0.250)
%       fs:         sampling frequency (default: 1KHz)
%       fid_vec:    if some subsegments should not be used for finding the
%                   optimal threshold of the P&Tthen input the indices of
%                   the corresponding points here
%       SIGN_FORCE: force sign of peaks (positive value/negative value).
%                   Particularly usefull if we do window by window detection and want to
%                   unsure the sign of the peaks to be the same accross
%                   windows (which is necessary to build an FECG template)
%       debug:      1: plot to bebug, 0: do not plot
%
% outputs
%   qrs_pos:        indexes of detected peak (in samples)
%   sign:           sign of the peaks (a pos or neg number)
%   en_thres:       energy threshold used
%
%
%
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 28-01-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
REF_PERIOD = 0.250; 
THRES = 0.6; 
fs = 1000; 
fid_vec = [];
SIGN_FORCE = [];
debug = 0;

switch nargin
    case 2
        REF_PERIOD=varargin{1};
    case 3
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2};
    case 4
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2};
        fs=varargin{3};  
    case 5
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3}; 
        fid_vec=varargin{4};
    case 6
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3};
        fid_vec=varargin{4}; 
        SIGN_FORCE=varargin{5};
    case 7
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3};
        fid_vec=varargin{4};
        SIGN_FORCE=varargin{5};          
        debug=varargin{6};
    otherwise
        error('qrs_detect: wrong number of input arguments \n');
end

[a b] = size(ecg);
if(a>b); NB_SAMP=a; elseif(b>a); NB_SAMP=b; ecg=ecg'; end;
tm = 1/fs:1/fs:ceil(NB_SAMP/fs);

% == constants
MED_SMOOTH_NB_COEFF = round(fs/100);
INT_NB_COEFF = round(7*fs/256); % length is 7 for fs=256Hz
LOW_CUT_FREQ = 5;
HIGH_CUT_FREQ = 45;
SEARCH_BACK = 1; % perform search back (FIXME: should be in function param)
MAX_FORCE = []; % if you want to force the energy threshold value (FIXME: should be in function param)

try
    % == prefiltering
    [b_lp,a_lp] = butter(5,HIGH_CUT_FREQ/(fs/2),'high');
    [b_bas,a_bas] = butter(2,LOW_CUT_FREQ/(fs/2),'high'); % FIXME: use more coefficients?
    ecg = ecg-mean(ecg);                    % (1) centre
    bpfecg = ecg'-filtfilt(b_lp,a_lp,ecg'); % (2) remove higher freq (zero phase)
    bpfecg = filtfilt(b_bas,a_bas,bpfecg);  % (3) remove baseline (zero phase)

    % == P&T operations
    dffecg = diff(bpfecg');  % (4) differentiate (one datum shorter)
    sqrecg = dffecg.*dffecg; % (5) square ecg
    intecg = filter(ones(1,INT_NB_COEFF),1,sqrecg); % (6) integrate
    mdfint = medfilt1(intecg,MED_SMOOTH_NB_COEFF);  % (7) smooth
    delay  = ceil(INT_NB_COEFF/2); 
    mdfint = circshift(mdfint,-delay); % remove filter delay for scanning back through ECG

    % look for some measure of signal quality with signal fid_vec? (FIXME)
    if isempty(fid_vec); mdfintFidel = mdfint; else mdfintFidel(fid_vec>2) = 0; end;

    % == P&T threshold
    if NB_SAMP/fs>90; xs=sort(mdfintFidel(fs:fs*90)); else xs = sort(mdfintFidel(fs:end)); end;

    if isempty(MAX_FORCE)
       if NB_SAMP/fs>10
            ind_xs = ceil(98/100*length(xs)); 
            en_thres = xs(ind_xs); % if more than ten seconds of ecg then 98% CI
        else
            ind_xs = ceil(99/100*length(xs)); 
            en_thres = xs(ind_xs); % else 99% CI  
        end 
    else
       en_thres = MAX_FORCE;
    end

    % build an array of segments to look into
    poss_reg = mdfint>(THRES*en_thres); 

    % in case empty because force threshold and crap in the signal
    if isempty(poss_reg); poss_reg(10) = 1; end;

    % == P&T QRS detection & search back
    if SEARCH_BACK
        indAboveThreshold = find(poss_reg); % ind of samples above threshold
        RRv = diff(tm(indAboveThreshold));  % compute RRv
        medRRv = median(RRv(RRv>0.01));
        indMissedBeat = find(RRv>1.5*medRRv); % missed a peak?
        % find interval onto which a beat might have been missed
        indStart = indAboveThreshold(indMissedBeat);
        indEnd = indAboveThreshold(indMissedBeat+1);

        for i=1:length(indStart)
            % look for a peak on this interval by lowering the energy threshold
            poss_reg(indStart(i):indEnd(i)) = mdfint(indStart(i):indEnd(i))>(0.5*THRES*en_thres);
        end
    end

    % find indices into boudaries of each segment
    left  = find(diff([0 poss_reg'])==1);  % remember to zero pad at start
    right = find(diff([poss_reg' 0])==-1); % remember to zero pad at end

    % looking for max/min?
    if SIGN_FORCE
        sign = SIGN_FORCE;
    else
        nb_s = length(left<30*fs);
        loc  = zeros(1,nb_s);
        for j=1:nb_s
            [~,loc(j)] = max(abs(bpfecg(left(j):right(j))));
            loc(j) = loc(j)-1+left(j);
        end
        sign = mean(ecg(loc));  % FIXME: change to median?  
    end

    % loop through all possibilities 
    compt=1;
    NB_PEAKS = length(left);
    maxval = zeros(1,NB_PEAKS);
    maxloc = zeros(1,NB_PEAKS);
    for i=1:NB_PEAKS
        if sign>0
            % if sign is positive then look for positive peaks
            [maxval(compt) maxloc(compt)] = max(ecg(left(i):right(i)));
        else
            % if sign is negative then look for negative peaks
            [maxval(compt) maxloc(compt)] = min(ecg(left(i):right(i)));
        end
        maxloc(compt) = maxloc(compt)-1+left(i); % add offset of present location

        % refractory period - has proved to improve results
        if compt>1
            if maxloc(compt)-maxloc(compt-1)<fs*REF_PERIOD && abs(maxval(compt))<abs(maxval(compt-1))
                maxloc(compt)=[]; maxval(compt)=[];
            elseif maxloc(compt)-maxloc(compt-1)<fs*REF_PERIOD && abs(maxval(compt))>abs(maxval(compt-1))
                maxloc(compt-1)=[]; maxval(compt-1)=[];
            else
                compt=compt+1;
            end
        else
            % if first peak then increment
            compt=compt+1;
        end
    end

    qrs_pos = maxloc; % datapoints QRS positions 
    R_t = tm(maxloc); % timestamps QRS positions
    R_amp = maxval; % amplitude at QRS positions
    hrv = 60./diff(R_t); % heart rate

catch ME
    rethrow(ME);
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    qrs_pos = [1 10 20]; sign = 1; en_thres = 0.5; 
end

% == plots
if debug
    figure(1);
    ax(1) = subplot(4,1,1); plot(tm,ecg); hold on;plot(tm,bpfecg,'r')
        title('raw ECG (blue) and zero-pahse FIR filtered ECG (red)'); ylabel('ECG');
        xlim([0 tm(end)]);  hold off;
    ax(2) = subplot(4,1,2); plot(tm(1:length(mdfint)),mdfint);hold on;
        plot(tm,max(mdfint)*bpfecg/(2*max(bpfecg)),'r',tm(left),mdfint(left),'og',tm(right),mdfint(right),'om'); 
        title('Integrated ecg with scan boundaries over scaled ECG');
        ylabel('Int ECG'); xlim([0 tm(end)]); hold off;
    ax(3) = subplot(4,1,3); plot(tm,bpfecg,'r');hold on;
        plot(R_t,R_amp,'+k');
        title('ECG with R-peaks (black) and S-points (green) over ECG')
        ylabel('ECG+R+S'); xlim([0 tm(end)]); hold off;
    ax(4) = subplot(4,1,4); hold off; plot(R_t(1:length(hrv)),hrv,'r+')
        hold on, title('HR')
        ylabel('RR (s)'); xlim([0 tm(end)]);
    
    linkaxes(ax,'x');
    set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end


% NOTES
%   Finding the P&T energy threshold: in order to avoid crash due to local 
%   huge bumps, threshold is choosen at 98-99% of amplitude distribution. 
%   first sec removed for choosing the thres because of filter init lag.
%   
%   Search back: look for missed peaks by lowering the threshold in area where the 
%   RR interval variability (RRv) is higher than 1.5*medianRRv
% 
%   Sign of the QRS (signForce): look for the mean sign of the R-peak over the
%   first 30sec when looking for max of abs value. Then look for the
%   R-peaks over the whole record that have this given sign. This allows to
%   not alternate between positive and negative detections which might
%   happen in some occasion depending on the ECG morphology. It is also
%   better than forcing to look for a max or min systematically.








