function [y,e,W] = FECGESN_lms(u,tar,mu,nC)
% least mean square adaptive filter.
%
% inputs
%   u:   input signal
%   tar: target signal
%   mu:  step size
%   nC:  number of filter coefficients
%   
% outputs
%   y: predicted signal
%   e: history error vector
%   W: history of the weights
%
%
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% University of Oxford, Intelligent Patient Monitoring Group
% joachim.behar@oxfordalumni.org
%
% Last updated : 28-01-2014
%
% Adapted from:
%       http://www.mathworks.co.uk/matlabcentral/answers/6524
%       (no copyright)
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

    N = length(u); % number of data samples
    y = zeros(N,1); % initialize filter output vector
    w = zeros(nC,1); % initialize filter coefficient vector
    e = zeros(N,1); % initialize error vector
    W = zeros(nC,N); % filter coefficient matrix for coeff. history
    for n = 1:N
      if n <= nC
          k = n:-1:1;
          x1 = [u(k); zeros(nC-numel(k),1)];
      else
          x1 = u(n:-1:n-nC+1); % nC samples of u in reverse order
      end
      y(n) = w'*x1; % filter output
      e(n) = tar(n) - y(n); % error
      w = w + mu*e(n)'*x1; % update filter coefficients
      W(:,n) = w; % store current filter coefficients
    end
end
