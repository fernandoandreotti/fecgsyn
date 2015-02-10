function [M,N] = FECGSYN_kf_Linearization(x,alphai,bi,tetai,w,fs,flag)
%LINEARIZATION linearization of the non linear system
%   linearization of the non linear system
%
% inputs:
% x: non linear
% Wmean: mean process noise vector
% Vmean: mean observation noise vector
% Modelparameters: Values of the parameters for the KF (values for
% gaussians)
% w: average heart-rate in rads.
% fs: sampling frequency
% flag: If flag = 0 use EKF, if flag=1 use EKS
%
% outputs:
% M: H Jacobian
% N: V Jacobian
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

