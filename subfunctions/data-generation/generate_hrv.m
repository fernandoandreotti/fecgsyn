function [theta,w] = generate_hrv(strhrv,n,fs,theta0)
% function [theta,w] = generate_hrv(strhrv,n,fs,theta0)
% generate variable heart rate (HR). Add suddent change of HR 
% in the middle of a time interval. This is meant, as an example, 
% to model high HR variation to test the robustness of a NI-FECG 
% extraction method to adapt to change of ECG morphology due to 
% HR variability (HRV).
%
% Input:
%   strhrv:  
%       - strhrv.hr:        mean heart rate [bpm]
%       - strhrv.lfhfr:     low to high frequency ratio of the two Gaussians
%                           hf->simulates respiratory sinus arrythmia
%                           lf->simulates Mayer
%       - strhrv.hrstd:     standard deviation of heart rate [bpm]
%       - strhrv.flo:       center freqency of low frequency Gaussian (Mayer) [Hz]
%       - strhrv.flhi:      center freqency of high frequency Gaussian (RSA) [Hz]
%       - strhrv.acc:       amplitude of acceleration(positive) or deceleration
%                           (negative) [bpm]
%       - strhrv.typeacc:   type of curve used, e.g. 'none', 'tanh', 'gauss', 
%                           'mexhat'
%       - strhrv.accmean    point of inflexion point or mean in percent [-1,1]
%       - strhrv.accstd     standard deviation (case necessary) 
%   n:                      number of samples   
%   fs:                     sampling frequency    
%   theta0:                  initial phase of the synthetic dipole
%
% Output:
%   theta:                   generated phase signal theta(t)
%   w:                      angular frequency
%
%
% Examples:
% TODO
%
% See also:
% run_ecg_generator
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

% == checking inputs
if nargin<=2; error('generate_hrv: not enough input arguments'); end;
if nargin<3; fs=250; end;

% == constants
NB_SUB = ceil(strhrv.hr*(n/fs)/60); % number of subdivisions
RRtemp = rrprocess(NB_SUB, strhrv.hr, strhrv.lfhfr, strhrv.hrstd,[], strhrv.flo, strhrv.fhi,[],[]); % generate RR
NB_SUB = length(RRtemp); % NB_SUB can change from +/- 1 so adjust
HRV = 60./RRtemp;

% == core function to modulate heart rate
% the example functions always perform the modulation around the half of the measurement
switch strhrv.typeacc
    %   strhrv.typeacc:    acceleration of  heart rate with some pre-defined functions around middle point of RR series:
    %          'none'            No acceleration/deceleration
    %          'tanh'            Hyperbolical tangent
    %          'mexhat'          Mexican hat
    %          'gauss'           Gaussian
    case 'none'
        rr = 60/strhrv.hr;
        RR = repmat(rr,NB_SUB,1);
    case 'nsr'
        % = normal sinus rythm
        RR = 60./HRV;
    case 'tanh'
        % = hyperbolic tangent (tanh)
        tmp = linspace(-20,20,NB_SUB);
        strhrv.accmean = strhrv.accmean*20;   % normalizing mean
        tmp = ((tanh(tmp-strhrv.accmean)+1)/2)*strhrv.acc;    % make HR suddenly ramp up in the middle, using tanh

        HRV = HRV+tmp';
        RR = 60./HRV;
     case 'mexhat'
        % = sombrero modulation (for foetal response to contraction)
        x = linspace(-3,3,NB_SUB);
        strhrv.accmean = strhrv.accmean*3;   % normalizing mean
        c = 2/(sqrt(3*strhrv.accstd)*pi^(1/4));
        sombrero = (c*(1-(x-strhrv.accmean).^2/strhrv.accstd^2).*exp(-(x-strhrv.accmean).^2/(2*strhrv.accstd)))';
        sombrero = strhrv.acc*sombrero/max(sombrero);        % amplitude correction    

        HRV = HRV + sombrero; % making an acceleration/decceleration HR
        RR = 60./HRV;
    case 'gauss'
        % = gaussian modulation
        x = linspace(-3,3,NB_SUB);
        strhrv.accmean = strhrv.accmean*3;   % normalizing mean
        gauss = 1/(sqrt(2*pi)*strhrv.accstd)*exp(-(x-strhrv.accmean).^2/(2*strhrv.accstd^2));
        gauss = strhrv.acc*gauss/max(gauss);        % amplitude correction    

        HRV = HRV + gauss';  % generating an acceleration/decceleration HR
        RR = 60./HRV;
end

% == Generating a phase trend w(t)
csum = cumsum(RR); RR(csum>n/fs)=[]; 
csum(round(csum*fs)>n) = [];
RR_rs = interp1(cumsum(RR),RR,RR(1):1/fs:sum(RR));
nbm_str = ceil(RR(1)*fs);               % number of points missed at the begining
nbm_end = ceil(n-sum(RR)*fs);           % number of points missed at the end
RR_rs = [repmat(RR_rs(1),1,nbm_str) RR_rs repmat(RR_rs(end),1,nbm_end)];
hr = 1./RR_rs(1:n);             % heart rate in Hz
w = 2*pi*hr;                    % angular frequency

% == Generating theta
theta = FECGSYN_kf_phasecalc(round(csum*fs),n);
nshift = find(theta>theta0,1,'first'); % considering theta0
theta = FECGSYN_kf_phasecalc(round(csum*fs),n+nshift-1);
theta(1:nshift-1) = [];
end