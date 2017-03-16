function sqi = psqi(signal,fs)
%pSQI Power of QRS SQI
% 
% Returns the relative power on band P(5-20Hz)/P(5-45Hz)
% 
% References
%
% [1] Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation from multiple asynchronous noisy
% sources using signal quality indices and a Kalman filter." Physiological measurement 29.1 (2008): 15.
%
% 
% Input:
%   signal:         single channel (F)ECG [1xN double]
% 
% Output:
%   sqi:            resulting sSQI for segment
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

    

pband = bandpower(signal,fs,[5 15]);
ptot = bandpower(signal,fs,[5 40]);
sqi = 1-(pband/ptot);

end



