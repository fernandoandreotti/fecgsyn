function sqi = rsqi(qrs,fs,CI)
% rSQI Function
% 
% Compute smoothness of HRV or of HR time series given a confidence interval (CI).
% The underlying assumption for using this function is that the more smooth
% the QRS time series is the most likely it is to be a meaningful qrs time series.
%
% Inputs
%   qrs:    qrs fiducials (required, in data points number)
%   fs:     sampling frequency (default: 1kHz)
%   CI:     confidence interval (default: 0.96)
%
% Outputs
%   sqi     percentage of intervals inside CI
%
% References:
% [1] Behar, J., Oster, J., & Clifford, G. D. (2014). Combining and Benchmarking 
% Methods of Foetal ECG Extraction Without Maternal or Scalp Electrode Data. 
% Physiol. Meas., 35(8), 1569–1589.
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


if size(qrs,2) > size(qrs,1); qrs = qrs';end

% == manage inputs
if nargin<1; error('assess_regularity: wrong number of input arguments \n'); end;
if nargin<2; secDer = 1; end;
if nargin<3; fs = 1000; end;
if nargin<4; CI = 1; end;
if nargin<5; segL = 60; end;
if nargin<6; debug = 0; end;
if length(qrs) < 5;
    sqi = 0;
    disp('rSQI: too few qrs detection points')
    return
end

% == core function
%try
% == compute variability
hr = 60./(diff(qrs)/fs);


% rather than looking at the distribution of hr this option is looking at
% the distribution of hrv. This makes more sense because this way
% we are looking into high variation in deltaHR from a measure to the
% following one rather than the variability of absolute value of HR (which
% might be high if the foetus HR is changing significantly)



% now taking the derivative of smoothed hear rate. This will give the
% hrv
hrv = sort(diff(hr));
hrv_N = length(hrv);
% plot(hr); hold on, plot(yi,'r');
% we tolerate some mistakes using a confidence interval
if CI~=1
    CI_sup = ceil(hrv_N*(CI+(1-CI)/2));
    CI_inf = ceil(hrv_N*((1-CI)/2));
    hrv_CI = hrv(CI_inf:CI_sup);
else
    hrv_CI = hrv;
end
% output the std of the hrv
% SMI = std(hrv_CI); % OLD version

% new version (25-08-2013)
SMI = length(find(abs(hrv_CI)>30));
% this looks at the absolute number of outliers in hrv_CI with >
% 30bpm drop or increase from a point to the next.



sqi = 1 - SMI/hrv_N;
if sqi <0
    disp('What')
end

% == plots
%     hist(hrv_CI,40); xlabel('hr or hrv histogram');
%     title(['assess regularity in term of hr or hrv, REGULARITY:' SMI]);
%     set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');

end
