function sqi = xsqi(signal,qrs,fs,win)
%xSQI Extravagance SQI
% 
% Despite the fancy name, this function does some pretty boring
% calculations, while attempting to describe how different a QRS complex is
% from its surrounding signal. This is particularly important for FECG
% signals, which are often buried into noise.
% 
% Input:
%   signal:         single channel (F)ECG [1xN double]
%   qrs:            list with (F)QRS locations [1xNp double]
%   fs:             signal sampling frequency [Hz]
%   win:            half of the window length around (F)QRS complex [ms],
%                   same window is used to noise area
% 
% Output:
%   sqi:            resulting xSQI for segment
%
%
% References:
% [1] Andreotti, F., Riedl, M., Himmelsbach, T., Wedekind, D., Wessel, N., 
% Stepan, H., … Zaunseder, S. (2014). Robust fetal ECG extraction and detection
% from abdominal leads. Physiol. Meas., 35(8), 1551–1567. 
% https://doi.org/10.1088/0967-3334/35/8/1551
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


%% Generate a template
% removing extremities detections
if size(qrs,2) > size(qrs,1); qrs = qrs';end
win = win*fs;
hwin = round(win/2); %half window
extremities = (qrs <= win+hwin | qrs >= length(signal)-(win+hwin));        % test if there are peaks on the border that may lead to error
qrs = round(qrs(~extremities));                                    % remove extremity peaks
if length(qrs) <3
    disp('xsqi: skipping due to low number of beats available')
    sqi = 0;
    return
end
% Stacking cycles
M = arrayfun(@(x) signal(1,x-hwin:x+hwin)'.^2,qrs,'UniformOutput',false);    % creates beat matrix
M = cell2mat(M);                                                        % converting cell output to array form (matrix is 2*width+1 x
N = arrayfun(@(x) signal(1,[(x-win-hwin):(x-hwin-1) (x+hwin):(x+win+hwin)])'.^2,qrs,'UniformOutput',false);    % creates a surrounding matrix
N = cell2mat(N);                                                        % converting cell output to array form (matrix is 2*width+1 x

%% Check power
Psrd = sum(N); % Power of sorroundings
Ppeak=sum(M); % Power of peaks 2 times because interval is half as long
sqi=mean(Ppeak./(Psrd+Ppeak)); % percentage of power that the peaks represent



end

