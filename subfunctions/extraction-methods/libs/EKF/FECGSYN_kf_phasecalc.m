function phase = FECGSYN_kf_phasecalc(peaks,NbSamples)
% function phase = FECGSYN_kf_phasecalc(peaks,NbSamples)
% Phase calculation, this function generates a saw-tooth phase signal, based on QRS locations
%
%  Inputs
%       peaks:          fidutials location, in these points phase is zero
%       NbSamples:      number of samples on the signal
%
%  Output
%       phase:          generated phase signal for NbSamples
% 
% Reference
% (Andreotti 2014) Andreotti, F., Riedl, M., Himmelsbach, T., Wedekind, D., 
% Wessel, N., Stepan, H., … Zaunseder, S. (2014). Robust fetal ECG extraction and 
% detection from abdominal leads. Physiol. Meas., 35(8), 1551–1567. 
% 
% (OSET) Sameni, R. (2010). The Open-Source Electrophysiological Toolbox (OSET). 
% Retrieved from http://www.oset.ir
% 
% 
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% TODO
%
% See also:
% FECGSYN_kf_extraction
% FECGSYN_kf_linearization
% FECGSYN_kf_EKFilter
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
% FECGSYN_kf_ECGmodelling
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
phase = zeros(1,NbSamples);
m = diff(peaks);            % gets distance between peaks
% = dealing with borders (first and last peaks may not be full waves)

% first interval uses second interval as reference
L = peaks(1);       %length of first interval
if isempty(m)       % only ONE peak was detected
    phase(1:NbSamples) = linspace(-2*pi,2*pi,NbSamples);
else
    phase(1:L) = linspace(2*pi-L*2*pi/m(1),2*pi,L);
    % beats in the middle
    for i = 1:length(peaks)-1;      % generate phases between 0 and 2pi for almos all peaks
        phase(peaks(i):peaks(i+1)) = linspace(0,2*pi,m(i)+1);
    end                             % 2pi is overlapped by 0 on every loop
    % last interval
    % uses second last interval as reference
    L = length(phase)-peaks(end);   %length of last interval
    phase(peaks(end):end) = linspace(0,L*2*pi/m(end),L+1);
end
phase = mod(phase,2*pi);
phase(phase>pi) = phase(phase>pi)- 2*pi;
end