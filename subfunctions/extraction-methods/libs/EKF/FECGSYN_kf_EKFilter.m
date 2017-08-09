function [Xhat,X0,P0,e] = FECGSYN_kf_EKFilter(Y,X0,P0,Q0,R0,Wmean,Vmean,ModelParam,w,fs,flag)
% Extended kalman filter or smoothing
% The code is based on the code provided in OSET Toolbox (http://www.oset.ir/)
% by Dr. Reza Sameni and in (Andreotti 2014)
%
% Inputs
%       Z:        matrix of observation signals (samples x 2). First column 
%                 corresponds to the phase observations and the second column 
%                 corresponds to the noisy ECG.
%       X0:       initial estimate for state vector
%       P0:       covariance matrix of the initial state vector
%       Q:        covariance matrix of the process noise vector
%       R:        covariance matrix of the observation noise vector
%       Wmean:    mean process noise vector
%       Vmean:    mean observation noise vector
%       ModelParam:     Values of the parameters for the KF (values for
%                       gaussians)
%       w:        average heart-rate in rads.
%       fs:       sampling frequency
%       flag:     if flag = 0 use EKF, if flag=1 use EKS
%       u:        control input (not used)
% 
% Output
%       Xhat: state vectors estimated by the EKF or EKS (depending on the flag).
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

%% Parameters Initialization

% Define the model parameters
N = (length(ModelParam)/3);
alphai = ModelParam(1:N);
bi = ModelParam(N+1:2*N);
tetai = ModelParam(2*N+1:3*N);
Samples = length(Y);
L = length(X0);

% Initalizing filter
R = R0;
Q = Q0;
X = X0;
P = P0;
H = [eye(2) zeros(2,L-2)];   % measuring matrix,

% Allocating variables
Xhat = zeros(L,Samples);
Phat = zeros(L,L,Samples);


%% Filtering
e = zeros(2,Samples);
for k = 1 : Samples
    %= Solution for some numerical problems of singularities
    idx = find(bi<0.005);   % look for narrow Gaussians
    if ~isempty(idx)
        alphai(idx) = 0;     % if too narrow, alpha is zeroed so it has no contribution
    end
    
    % Prevention of 'Xminus' mis-calculations on phase jumps
    if(abs(X(1)-Y(1,k))>pi/5)
        X(1) = Y(1,k);
    end
    
    % FORECAST STNdrateEP (predict)
    X = StateProp(X,alphai,bi,tetai,w,fs);
    [A,dQ] = FECGSYN_kf_linearization(X,alphai,bi,tetai,w,fs,0);   % linearizing model(x,,w,fs,flag)
    P = A*P*A' + dQ*Q*dQ';                          % a priori error covariance matrix
    
    % DATA ASSIMILATION STEP (update)
    e(:,k) = Y(:,k)-H*X;                        % innovation
    S = H*P*H'+R;                               % 'innovation covariance'
    
    %= usual data assimilation
    K = (P*H')/S;                                 % compute Kalman gain
    X = X + K*e(:,k);                           % a posteriori state update
    P = (eye(L)-K*H)*P;                         % a posteriori state error covariance update   
    
    % Store results
    Xhat(:,k) = X; %Store the state (it's important that the state is a vector of 2 dimensions, phase and real state)
    Phat(:,:,k) = P; %Store error covaraiance
end

%% Smoothing: (Only enabled if flag)
if(flag==1)
    disp('Smoothing Estimation ..')
    Xhat = FECGSYN_kf_EKsmoothing(Phat,Xhat,Pbar,Xbar,Wmean,alphai,bi,tetai,w,fs);
end

end
%% Auxiliar Functions
function xout = StateProp(x,alphai,bi,tetai,w,fs)
% State propagation (prediction step) for state x
% The code is based on the code provided in OSET Toolbox (http://www.oset.ir/)
% by Dr. Reza Sameni and in (Andreotti 2014)
% 
% Input:
%    x:          current state x
%    alphai:     value for alpha parameters (Gaussian amplitude)
%    bi:         value for beta parameters (Gaussian width)
%    tetai:      value for theta parameters (Gaussian position - mean)
%    w:          average heart-rate in rads.
%    fs:         sampling frequency
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
% 
% FECGSYN_kf_extraction
% FECGSYN_kf_linearization
% FECGSYN_kf_ECGmodelling
% 
% fecgsyn toolbox, version 1.2, Jan 2017
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% University of Oxford, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
% 
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. 
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. 
%   Physiological Measurement.35 1537-1550. 2014.
% 
% 
% Last updated : 10-03-2016
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
dt = 1/fs;
xout(1,1) = x(1) + w*dt;            % teta state variable
if(xout(1,1)>pi),
    xout(1,1) = xout(1,1) - 2*pi;
end
dtetai = rem(xout(1,1) - tetai,2*pi);
xout(2,1) = x(2) - dt*sum(w*alphai./(bi.^2).*dtetai.*exp(-dtetai.^2./(2*bi.^2))); % z state variable
end
