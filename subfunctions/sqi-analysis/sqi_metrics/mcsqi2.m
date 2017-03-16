function sqi = mcsqi2(raw,fecg,fs)
%mcSQIb Spectral Coherence metric
% 
% Returns the coherence SQI between raw and extracted FECG signals.
% 
% 
% Input:
%   raw:         Raw data segment
%   fecg:        Residual data segment
%   mecg:        Chest lead segment
% 
% Output:
%   sqi:            resulting sSQI for segment
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

    

%[Prm,f]=mscohere(raw,mecg,[],[],1024,fs);
[Prf,f]=mscohere(raw,fecg,[],[],1024,fs);
%[Pmf,f]=mscohere(mecg,fecg,[],[],1024,fs);


sqi1 = mean(Prf(f>=60&f<=100)); % average coherence for 60-100Hz (no high frequent noise added)
%sqi2 = 1-mean(Pmf(f<=100));     % 1- average cross-spectrum MECG-FECG for f<100Hz

sqi = sqi1;


end



