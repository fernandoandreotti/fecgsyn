function sqi = isqi(qrs,win,fs)
%BSQI_MATLAB Calculate iSQI of all channels
%
% Input:
%  qrs:          Cell containing annotations for all channels
%  win:           Acceptance interval for TP detection (in s)
%  fs:            Sampling frequency (in Hz)
%
% Output:
%  sqi           Array containing for each channel iSQI metric (percentage
%                of detections on each channel, present in all channels)
% References
%
% [1] Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2015).
% Multimodal heart beat detection using signal quality indices, Physiological Measurement
% 36 (2015): 1665-1677.
% [2] Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2014).
% R-peak estimation using multimodal lead switching, Computing in Cardiology Conference
% (CinC), 2014, Vol. 41, pp. 281-284.
%
% [3] Behar, Joachim, et al. "ECG signal quality during arrhythmia and its
% application to false alarm reduction." Biomedical Engineering, IEEE Transactions on 60.6 (2013): 1660-1666.
%
% [4] Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation from multiple asynchronous noisy
% sources using signal quality indices and a Kalman filter." Physiological measurement 29.1 (2008): 15.
%
% Multimodal peak detection using ECG, ABP, PPG or SV
% Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J.
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%

% == managing inputs
if nargin<2; win=0.05; end;
if nargin<3; fs=250; end;

win = win * fs;
NCHAN = length(qrs);
% == core function
sqi = zeros(NCHAN,1);
for ch = 1:NCHAN
    qrstmp = qrs;
    qrstmp(ch) = [];
    ref = qrs{ch};
    if size(ref,2) > size(ref,1); ref = ref';end
    refbin = ones(size(ref));
    for ch2 = 1:length(qrstmp)
        qrstest = qrs{ch2};
        if size(qrstest,2) > size(qrstest,1); qrstest = qrstest';end
        if isempty(qrstest)||isempty(ref); refbin = zeros(size(refbin));continue;end;
        [IndMatch,Dist] = dsearchn(ref,qrstest);         % closest ref for each point in test qrs        
        refbin(IndMatch(Dist<win)) = refbin(IndMatch(Dist<win))+1;         % number of channels where this annotation is present
    end    
    sqi(ch) = sum(refbin==NCHAN)/length(refbin);
end

sqi(isnan(sqi)) = 0;
end
