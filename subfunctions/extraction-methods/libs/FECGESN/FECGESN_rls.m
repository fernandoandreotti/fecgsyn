function [out,P,w,e] = FECGESN_rls(state,tar,lda,P,w)
% run the recursive least square algorithm. This function consists of one
% step of the RLS i.e. given the covariance matrix and weight at a given
% timestep the filtered new datapoint is computed and the covariance matrix 
% and the weights are updated. In other words you have to put a loop around
% it to filter a set of datapoits. This is made on purpose in case the
% entry datapoints are projected in a higher dimentional space/go through 
% some transformation (as it is the case for the ESN) before applying the
% RLS step.
%
% inputs
%   state: state of the ESN (or any state of a network neurons) at a given 
%          timestep. In the more classical case this corresponds to the last 
%          N datapoints of the input signal
%   tar:   target point
%   lda:   'lambda' = forgetting factor (in [0 1], 0-> mostly current 
%          datapoint error is considered, 1-> past estimation errors 
%          are also taken into account)
%   P:     covariance matrix
%   w:     filter weights
%
% outputs
%   out:   output of the RLS
%   P:     updated covariance matrix
%   w:     updated filter weights
%   e:     estimation error (absolute error)
%
%
% Reference:
% RLS was implemented as described in
% [1] Petrenas, Andrius, et al. "An echo state neural network for qrst 
% cancellation during atrial fibrillation." (2012): 1-1.
%
% The original paper describing this algorithm is
%
% [2] Douglas, S. C. "Numerically-robust O (N< sup> 2</sup>) RLS algorithms 
% using least-squares prewhitening." Acoustics, Speech, and Signal Processing, 
% 2000. ICASSP'00. Proceedings. 2000 IEEE International Conference on. Vol. 1. IEEE, 2000.
%
%
% --
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% University of Oxford, Intelligent Patient Monitoring Group
% joachim.behar@oxfordalumni.org
%
% Last updated : 28-01-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == check inputs
if nargin<5; error('run_rls: wrong number of input arguments'); end;
if lda>1 || lda<0; error('run_rls: lambda MUST be in range [0 1]'); end;

% == compute output given past estimated w and P matrices
out = w*state;  

% == RLS
v = P*state;
u = P'*v;
k = 1/(lda+v'*v+sqrt(lda*(lda+v'*v)));
e = tar-out; % estimate error
w = w'+(e*u)/(lda+v'*v); % update the weights 
P = (P-k*v*u')/sqrt(lda); % estimate the covariance matrix

w = w';
end

