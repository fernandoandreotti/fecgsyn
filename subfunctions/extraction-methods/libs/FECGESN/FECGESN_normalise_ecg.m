function [necg,scalings] = FECGESN_normalise_ecg(ecg,start,stop)
% normalisation function. See STEP 1-3 in the code below to understand 
% what it does.
% 
% inputs:
%   ecg: ecg before normalization
%   start/stop:  starting and ending points (in samples) defining
%   representative ecg segments on which to compute scalings and shifts.
%   The whole ecg is not used for that otherwise the scaling would be
%   dependant on the local artefacts etc.
%   
% outputs:
%   necg: the normalized ecg which values will be within [0 1].
%   scalings: scaling factor by which the ecg(s) have been multiplied 
%   (NOTE: that it does not take the detrand or tanh operations into account).
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
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == manage inputs
if size(ecg,2)>size(ecg,1)
    ecg = ecg';
end

% == core function
try
    % STEP 1: Normalizes a multivariate dataset data (where each data vector is a row in data) 
    % by scaling such that in each column the min/max becomes 0/1. No
    % shifting is performed (it is assumed to be zero DC after
    % prefiltering).
    scalings=1./(max(ecg(start:stop,:))-min(ecg(start:stop,:)));
    
    % STEP 2: detrend the ecg
    rescaled_ecg = bsxfun(@times,ecg,scalings);
    detrendedInput = bsxfun(@minus,rescaled_ecg,mean(rescaled_ecg(start:stop,:))); 
    
    % STEP 3: take the tanh of the sigal . This step is applied in order to
    % avoid outliers whih could result in the reservoir state or LMS
    % weights to take unexpected values that are not covered by the normal
    % trajectory given the optimized global parameters or output learned.
    necg = tanh(detrendedInput);
catch ME
    rethrow(ME);
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    necg = ecg;
end

end

 