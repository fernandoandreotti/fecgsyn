function sqi = csqi(signal,qrs,fs,win)
%cSQI Conformity SQI
% 
% Takes the median correlation coefficient between a template beat and each
% individual beat on a chunk of signal
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
% References:
%  [1] Andreotti, F., Riedl, M., Himmelsbach, T., Wedekind, D., Wessel, N., 
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
if size(qrs,1) > size(qrs,2); qrs = qrs';end
win = win*fs; % convert to samples
extremities = (qrs <= win | qrs >= length(signal)-win);        % test if there are peaks on the border that may lead to error
qrs = round(qrs(~extremities));                                    % remove extremity peaks
if length(qrs) <5
    disp('csqi: skipping due to low number of beats available')
    sqi = 0;
    return
end
% Stacking cycles
M = arrayfun(@(x) signal(1,x-win:x+win)',qrs,'UniformOutput',false);    % creates a maternal beat matrix
M = cell2mat(M);                                                        % converting cell output to array form (matrix is 2*width+1 x
avgbeat = median(M,2)';                                                 % generates template from detected QRS

%% Calculate individual corrcoef
corr = zeros(1,size(M,2));
for i=1:size(M,2)
    c = corrcoef(avgbeat,M(:,i));
    corr(i) = c(1,2);
end

%% Median correlation
sqi = mean(corr);
sqi(sqi<0) = 0; % negative correlations are zero

end

