function [ECGmean,ECGsd,meanPhase] = FECGSYN_kf_templategen(ecg,phase,NB_BINS)
% function [ECGmean,ECGsd,meanPhase] = FECGSYN_kf_templategen(ecg,phase,NB_BINS)
% Calculation of the mean and SD of ECG waveforms in different beats.
% Although this function structure is based on OSET's toolbox, the
% averaging procedure itself is based on Dr. Oster's approach for stacking
% and averaging beats.
%
%  Inputs
%       ecg:            input ECG signal
%       phase:          ECG phase
%       NB_BINS:        number of desired phase bins
%
%  Outputs
%       ECGmean:        mean ECG beat
%       ECGsd:          standard deviation of ECG beats
%       meanPhase:      the corresponding phase for one ECG beat
%
% References
% (Andreotti 2014) Andreotti, F., Riedl, M., Himmelsbach, T., Wedekind, D., 
% Wessel, N., Stepan, H., … Zaunseder, S. (2014). Robust fetal ECG extraction and 
% detection from abdominal leads. Physiol. Meas., 35(8), 1551–1567. 
% 
% (OSET) Sameni, R. (2010). The Open-Source Electrophysiological Toolbox (OSET). 
% Retrieved from http://www.oset.ir
% 
% (Oster 2015) Oster, J., Behar, J., Sayadi, O., Nemati, S., Johnson, A., & Clifford, G. (2015). 
% Semi-supervised ECG Ventricular Beat Classification with Novelty Detection Based on Switching 
% Kalman Filters. IEEE Trans. Biomed. Eng., 62(9), 2125–2134. http://doi.org/10.1109/TBME.2015.2402236
% 
% 
%
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% TODO
%
% See also:
% FECGSYN_kf_modelling
% FECGSYN_kf_extraction
% FECGSYN_kf_linearization
% FECGSYN_kf_EKFilter
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
%
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

ini_cycles = find(phase(2:end)<0&phase(1:end-1)>0)+1; % start of cycles
cycle_len = diff(ini_cycles); % distance between cycles
end_cycles = ini_cycles(1:end-1)+cycle_len-1; % start of cycles
meanPhase = linspace(-pi,pi,NB_BINS);
% stacking cycles
cycle = arrayfun(@(x) interp1(phase(ini_cycles(x):end_cycles(x)),...
    ecg(1,ini_cycles(x):end_cycles(x)),meanPhase,'spline'),...
    1:length(ini_cycles)-1,'UniformOutput',0);
cycle = cell2mat(cycle');
ECGmean = mean(cycle);
ECGsd = std(cycle);

end

