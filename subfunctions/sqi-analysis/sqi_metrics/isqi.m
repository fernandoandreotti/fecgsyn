function sqi = isqi(qrs,win,fs)
% iSQI Function
% Calculate iSQI of all channels, i.e. for each detection in which
% percentage of channels it was detected e.g. on a multichannel ECG.
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
% [1] Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation from multiple asynchronous noisy
% sources using signal quality indices and a Kalman filter." Physiological measurement 29.1 (2008): 15.
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
