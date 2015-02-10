function [Xeks] = FECGSYN_kf_EKsmoothing(Phat,Xhat,Pbar,Xbar,Wmean,alphai,bi,tetai,w,fs)
%EKS Smoothing add-on for Extended Kalman Filter
%   Called inside EKF filter (just before filtering)
%
% inputs:
% Phat: State covariance matrix
% Xhat: State vector estimation
% Pbar: Stored covariance matrix
% Xbar: Stor state vector estimation
% Wmean: mean process noise vector
% alphai: Gaussian amplitude of the model (vector)
% betai: Gaussian width of the model (vector)
% tetai: Gaussian angular position of the model (vector)
% w: average heart-rate in rads.
% fs: sampling frequency

% outputs:
% Xeks: State vector estimation of EKS
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

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
    %Storing initialization
    PSmoothed = zeros(size(Phat));
    Xeks = zeros(size(Xhat));
    %Define number of samples
    Samples = length(Xhat);

    PSmoothed(:,:,Samples) = Phat(:,:,Samples);
    Xeks(:,Samples) = Xhat(:,Samples);
    for k = Samples-1 : -1 : 1,
        [A] = Linearization(Xhat(:,k),alphai,bi,tetai,w,fs,0);
        S = Phat(:,:,k) * A' * inv(Pbar(:,:,k+1));
        Xeks(:,k) = Xhat(:,k) + S * (Xeks(:,k+1) - Xbar(:,k+1));
        PSmoothed(:,:,k) = Phat(:,:,k) - S * (Pbar(:,:,k+1) - PSmoothed(:,:,k+1)) * S';

    end
end
