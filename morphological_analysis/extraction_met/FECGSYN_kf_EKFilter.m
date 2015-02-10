function [Xhat,X0,P0,e] = FECGSYN_kf_EKFilter(Z,X0,P0,Q,R0,Wmean,Vmean,ModelParam,w,fs,flag,u)
%% EXTENDED KALMAN FILTER/SMOOTHER Extended kalman filter or smoothing
%
% > Inputs
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

% > Output
%       Xhat: state vectors estimated by the EKF or EKS (depending on the flag).
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
%
%
% Current version:
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
% Last updated : 24-07-2014
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%
%% Equation Initialization
Inits2 = [ModelParam w fs];
StateProp(Inits2);           % Initialize state equation

%% Parameters Initialization

% Define the model parameters exported
L = (length(ModelParam)/3);
alphai = ModelParam(1:L);
bi = ModelParam(L+1:2*L);
tetai = ModelParam(2*L+1:3*L);

%lenghts calculation
Samples = length(Z);
L = length(X0);

%Initial estimates
Pminus = P0;
Xminus = X0;

%Covariance matrix of the observation noise vector initialization
R = R0;

%Storing initialization
Xbar = zeros(L,Samples);
Pbar = zeros(L,L,Samples);
Xhat = zeros(L,Samples);
Phat = zeros(L,L,Samples);


%% Filtering
e = zeros(2,Samples);
B = [0;1];  % control matrix
for k = 1 : Samples
    
    % Prevention of 'Xminus' mis-calculations on phase jumps
    if(abs(Xminus(1)-Z(1,k))>pi)
        Xminus(1) = Z(1,k);
    end
    
    % Store results
    Xbar(:,k) = Xminus'; %Store the state (it's important that the state is a vector of 2 dimensions, phase and real state)
    Pbar(:,:,k) = Pminus'; %Store error covaraiance
    
    XX = Xminus; %Store a priori state to don't overlap it
    PP = Pminus; %Store a priori error covariance to don't overlap it
    
    for jj = 1:size(Z,1);
        %MEASUREMENT UPDATE (CORRECT)
        %A posteriori updates)
        Yminus = ObservationProp(XX,Vmean); %Calculate output estimate (both state and theta)
        YY = Yminus(jj);
        [HH,VV] = FECGSYN_kf_linearization(XX,alphai,bi,tetai,w,fs,1);  % Linearized observation eq.
        H = HH(jj,:); % Get row jj   -> H Jacobian
        V = VV(jj,:); % Get  jj   -> V
        
        K = PP*H'/(H*PP*H' + V*R(jj,jj)*V');                  % Compute Kalman Gain
        XX = XX + K*(Z(jj,k)-YY);                             % Update estimate with measurement
        PP = (eye(L)-K*H)*PP;                               % Update the error covariance
    end
    
    % TIME UPDATE (PREDICT)
    e(:,k) = XX-Xminus;
    Xminus = StateProp(XX,B,u(k));                             % Project the state ahead
    [A,F] = FECGSYN_kf_linearization(XX,alphai,bi,tetai,w,fs,0);   % Linearized equations
    Pminus = A*PP*A' + F*Q*F';                             % Project the error covariance ahead
    
    % Store results
    Xhat(:,k) = XX';
    Phat(:,:,k) = PP'; %EKS filter use this information to calculate backwards
    
    
end

%Create initial values, in case of entering again on KF they will be
%used
X0=XX;
P0=PP;
%% Smoothing: (Only enabled if flag)
if(flag==1)
    disp('Smoothing Estimation ..')
    Xhat = FECGSYN_kf_EKsmoothing(Phat,Xhat,Pbar,Xbar,Wmean,alphai,bi,tetai,w,fs);
end


%% Auxiliar Functions
function xout = StateProp(x,B,u)

% Make variables static
persistent tetai alphai bi fs w dt;

% Check if variables should be initialized
if nargin==1,
    % mean of the noise parameters
    % Inits = [alphai bi tetai w fs];
    L = (length(x)-2)/3;
    alphai = x(1:L);
    bi = x(L+1:2*L);
    tetai = x(2*L+1:3*L);
    w = x(3*L+1);
    fs = x(3*L+2);   
    dt = 1/fs;
    return
end

xout(1,1) = x(1) + w*dt;            % teta state variable
if(xout(1,1)>pi),
    xout(1,1) = xout(1,1) - 2*pi;
end

dtetai = rem(xout(1,1) - tetai,2*pi);
xout(2,1) = x(2) - dt*sum(w*alphai./(bi.^2).*dtetai.*exp(-dtetai.^2./(2*bi.^2))) + B(2)*u; % z state variable

function y = ObservationProp(x,v)

% Calculate output estimate
y = zeros(2,1);
y(1) = x(1) + v(1);   % teta observation
y(2) = x(2) + v(2);   % amplidute observation
