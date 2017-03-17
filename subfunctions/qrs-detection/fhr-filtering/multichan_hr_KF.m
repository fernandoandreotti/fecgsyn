function [fhr_final,fhr_estim,v_estim] = multichan_hr_KF(fhr,sqi,p,Q0,R0)
% Kalman filter for Multichannel Heart Rate estimation
%
% This function runs KF for a single camera PPG (cPPG) data, attempting to
% make a better HR estimate, based on the signal quality index (SQI) given
% as input, together with the previously estimated HR.
%
% Input
%       fhr:          HR signal (NxM) where M are diferent channels
%       sqi:        signal quality index for heart rate estimates in, same
%                   dimensions as y. NaNs are converted into zeros (NxM)
%       p:          Order of autogressive (AR) model, i.e. number of past
%                   points used
%       Q0:         Model error covariance matrix
%       R0:         Observation error covariance matrix
% Output
%       fhr_final:   Multichannel FHR estimate
%       fhr_estim:   Singlechannel FHR estimates
%       v_estim:     Estimated innovation signal
%
% Reference:
% Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation from multiple asynchronous noisy
% sources using signal quality indices and a Kalman filter." Physiological measurement 29.1 (2008): 15.
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


sqi(isnan(sqi)) = 10^-5; % replace NaNs for zeros
switch nargin
    case 1
        error('multichan_hr_KF: too few input arguments')
    case 2
        p = 1;
        Q0 = 1;
        R0 = 1e-3;
    case 3
        Q0 = 1;
        R0 = 1e-3;
    case 4
        R0 = 1e-3;
end
if size(fhr)~=size(sqi); error('multichan_hr_KF: SQI and FHR require same dimensions');end;


fhr_estim = zeros(size(fhr));
v_estim = fhr_estim;
for ch = 1:size(fhr,2)
    y = fhr(:,ch)';
    tsqi = sqi(:,ch)';
    %% Treat missing observations
    % Removing NaNs, keep last value
    idx = (~isnan(y)); %non nans
    vr = y(idx); %v non nan
    keepidx = cumsum(idx);
    keepidx(keepidx == 0) = [];
    vr = vr(keepidx); %use cumsum to build index into vr
    % Mirrowing begin of data
%     y = [150.*ones(1,length(y)-length(vr)) vr];
    
    %% KF parameters may change/be calibrated
    % other parameters
    y = [fliplr(y(1:p-1)) y];     % mirrowing beginning of dataset to filter first samples
    tsqi = [fliplr(tsqi(1:p)) tsqi];
    
    X0 = 150; % estimated bpm
    P0 = diag(100*ones(p,1));    
    [xhat,v,~] = runKF_SQI(y,X0,P0,R0,Q0,tsqi,p);
    xhat(:,1:p-1) = []; % removing mirrowed section
    v(:,1:p-1) = []; % removing mirrowed section
    fhr_estim(:,ch) = xhat(1,:)'; % saving for latter
    v_estim(:,ch) = v(1,:)';      % saving for latter
end
NSAMP = size(fhr,1);
NCHAN = size(fhr,2);
fhr_final  = zeros(NCHAN,NSAMP);

%% Merging multiple channels as in Oster 2015 (Signal quality indices for state space ...)
sqi(sqi==0) = 10e-10;
v_estim(v_estim==0) = 1e-5; % zeros generate singularities (null divisions)
m = (v_estim./sqi).^2; % weighting factors from innovation and SQI
weights = zeros(size(fhr_estim));
for k = 1:NSAMP % samples    
    for n = 1:NCHAN
        % numerator
        num = m(k,:)';
        num(n,:) = [];
        num = prod(num); % multiplication
        
        % denominator
        combs = nchoosek(1:NCHAN,NCHAN-1);
        cden = zeros(size(combs,1),1);
        for r = 1:size(combs,1); 
            cden(r) = prod(m(k,combs(r,:))); 
        end
        den = sum(cden);                            
        fhr_final(n,k) = (num/den)*fhr_estim(k,n); 
        weights(k,n) = (num/den);
    end
end
fhr_final = sum(fhr_final)'; % final merge
end

function [Xhat,v,Perror] = runKF_SQI(y,X0,P0,R0,Q0,SQI,p)
%Initial estimates
P = P0;
x = zeros(p,1);
x(1) = X0;
Q = zeros(p);Q(1) = Q0;
Samples = length(y);
Xhat = zeros(p,Samples);
Phat = zeros(p,p,Samples);
% H = eye(p);
H = zeros(1,p); H(1) = 1;

aa=flip(linspace(0.2,1,p)); % linear regression for filter coefficients
sumaa=sum(aa);
A = [aa/sumaa.*ones(1,p); eye(p-1), zeros(p-1,1)];
v = zeros(p,Samples); % innovation

% KF loop
for k = 1:Samples
    
    % Prediction for state vector and covariance:
    x = A*x;                   % a priori state estimate
    P = A*P*A' + Q;            % a priori error cov. matrix update
%     R = R0*(1-SQI(k));  % observational cov. matrix SQI-dependent
     R = R0*exp(1/SQI(k)^2-1); % from Li et al 2008
     if isinf(R); R = 1e6*R0;  end % substituting infinity by large number
    %% Suggestion from Sebastian
    S = (H*P*H'+R);
    K = H*P/S;           % Compute Kalman gain
    v(:,k) = y(k)'-H*x;   % innovation
    x = x + K*v(:,k);              % state update
    P = P - K*H'*P;                 % covariance update
    Perror(:,:,k) = P;
    % Store results
    Xhat(:,k) = x;
    Phat(:,:,k) = P'; %EKS filter use this information to calculate backwards
end

v = v(1,:); % checking just first state
end