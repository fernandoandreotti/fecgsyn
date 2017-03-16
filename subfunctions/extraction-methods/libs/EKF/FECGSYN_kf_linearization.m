function [M,N] = FECGSYN_kf_linearization(x,alphai,bi,tetai,w,fs,flag)
%LINEARIZATION linearization of the non linear system
% The code is based on the code provided in OSET Toolbox (http://www.oset.ir/)
% by Dr. Reza Sameni and in (Andreotti 2014)
%
% inputs:
%   x: non linear state
%   alphai:     value for alpha parameters (Gaussian amplitude)
%   bi:         value for beta parameters (Gaussian width)
%   tetai:      value for theta parameters (Gaussian position - mean)
%   w:          average heart-rate in rads.
%   fs:         sampling frequency
%   flag: If flag = 0 use EKF, if flag=1 use EKS
%
% outputs:
% M: H Jacobian
% N: V Jacobian
%
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
% 
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% TODO
%
% See also:
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

L=length(alphai);
dt = 1/fs;
    

% Linearize state equation
if flag==0,
    dtetai = rem(x(1) - tetai,2*pi);

    M(1,1) = 1;                                                                     % dF1/dteta
    M(1,2) = 0;                                                                     % dF1/dz

    M(2,1) = -dt*sum( w*alphai./(bi.^2).*(1 - dtetai.^2./bi.^2).*exp(-dtetai.^2./(2*bi.^2)) ) ;    % dF2/dteta
    M(2,2) = 1 ;                                                                    % dF2/dz

    % W = [alpha1, ..., alphan, b1, ..., bn, teta1, ..., tetan, omega, N]
    N(1,1:3*L) = 0;
    N(1,3*L+1) = dt;
    N(1,3*L+2) = 0;

    N(2,1:L) = -dt*w./(bi.^2).*dtetai .* exp(-dtetai.^2./(2*bi.^2));
    N(2,L+1:2*L) = 2*dt.*alphai.*w.*dtetai./bi.^3.*(1 - dtetai.^2./(2*bi.^2)).*exp(-dtetai.^2./(2*bi.^2));
    N(2,2*L+1:3*L) = dt*w*alphai./(bi.^2).*exp(-dtetai.^2./(2*bi.^2)) .* (1 - dtetai.^2./bi.^2);
    N(2,3*L+1) = -sum(dt*alphai.*dtetai./(bi.^2).*exp(-dtetai.^2./(2*bi.^2)));
    N(2,3*L+2) = 1;

    % Linearize output equation
elseif flag==1,
    M=eye(2);
    N=eye(2);
end
end

