function sqi = mpsqi2(signal,qrs,fs)
%mpSQIb Spectral Power around MQRS complexes on extracted FECG
% 
% Returns spectral power around MQRS and its first harmonics in frequency domain
% up until 10 Hz.
% 
% 
% Input:
%   signal:         single channel (F)ECG [1xN double]
%   qrs:            maternal QRS locations (for calculating fundamental
%                   frequency)
%   fs:             sampling frequency
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


% Parameters
win = 0.3; % plus minus 0.3 Hz window
lband = 0.5;  % lowest band considered [Hz]
hband = 10;  % highest band considered [Hz]
Nest = 8;    % number of spected visible peaks on spektrum
% Calculate magnitude spectrum
N = length(signal);
freq = linspace(-fs/2, fs/2, N+1); freq(end) = [];
ydft = fft(signal);
shiftSpectrum = abs(fftshift(ydft));
% plot(freq,shiftSpectrum,'r');
% xlabel ('Frequency (Hz)');
% ylabel('Magnitude');

% Inteval of interest
fundfreq = fs/median(diff(qrs));        % finding median fundamental frequency
Nh = floor(hband/fundfreq);             % number of harmonics used
if (fundfreq > 3)||(fundfreq < 0.7) % skipping too high/low heart rates
    sqi = NaN; 
    warning('mpsqi:Maternal heart rate too high, skipping') 
    return; 
elseif isnan(fundfreq)||isinf(fundfreq)
    sqi = NaN; 
    warning('mpsqi:Maternal heart rate not detected, skipping') 
    return; 
end


freqcut = freq(freq<hband&freq>lband);
speccut = shiftSpectrum(freq<hband&freq>lband);  % chopping band of interest
harm = [1:Nh].*fundfreq;               % harmonics centers
harmlow = arrayfun(@(x) x-win,harm);  % lower bounds for harmonics
harmhigh = arrayfun(@(x) x+win,harm); % higher bounds for harmonics
idxpks = false(1,length(freqcut));
for p = 1:Nh    
    idxpks = idxpks|(freqcut<harmhigh(p)&freqcut>harmlow(p));
end

nhr = log(Nest/length(harm))+1; % normalizing factor to deal with different heart rates

sqi = 1-nhr*sum(speccut(idxpks).^2)/sum(speccut.^2);
sqi(sqi>1)=1;
sqi(sqi<0)= 0;
% Plot bands use
% figure
% fill(freqcut,max(speccut).*double(idxpks),[0.9 0.9 0.9])
% hold on
% plot(freqcut,speccut)
% hold off
% xlabel(num2str(harm))
% title(num2str(sqi))

end



