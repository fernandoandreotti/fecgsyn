function residual = FECGESN_lmsrls_canceller(ref_ecg,tar_ecg,fs,metStruct,debug)
% Adaptive Noise cancelling filter for FECG source separation.
%
% This function uses the least mean square (LMS) or the recursive least square
% (RLS) algorithm in order to remove the MECG contribution on the abdominal
% mixture. The LMS can be run 'online' (i.e the weights are left evolving)
% or 'offline' (the weights are optimised on a training set and then constant
% for the test set). The RLS is run online.
%
% inputs
%   ref_ecg:            reference chest channel
%   tar_ecg:            abdominal target channel
%   fs:                 sampling frequency
%   metStruct:          structure with LMS or RLS parameters
%       metStruct.mu
%       metStruct.N
%       metStruct.learningMode ('online'/'offline')
%       metStruct.method ('LMS'/'RLS')
%
% tar_ecg
%   residual:           FECG residual
%
% NOTE: time used for initialising the weights = 30sec (see TIME_INIT constant)
%
%
% References
%   1.  B. Widrow, J. Glover, J. McCool, J. Kaunitz, C. Williams, R. Hearn,
%       J. Zeidler, E. Dong, and R. Goodlin, Adaptive noise cancelling:
%       Principles and applications, Proceedings of the IEEE, vol. 63, no. 12,
%       pp. 1692?1716, 1975.
%   2.  lms implementation: http://www.mathworks.co.uk/matlabcentral/answers/6524
%       (no copyright)
%
%
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% University of Oxford, Intelligent Patient Monitoring Group
% joachim.behar@oxfordalumni.org
%
% Last updated : 22-07-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
if nargin<4; error('lms_canceller: wrong number of input arguments \n'); end;
if size(ref_ecg,1)>size(ref_ecg,2); ref_ecg = ref_ecg'; end;
if size(tar_ecg,1)>size(tar_ecg,2); tar_ecg = tar_ecg'; end;

% == general
TIME_INIT = 30;

% == core function
% == normalize the data
% assumes that the first 5sec are representative (no huge noise) for normalizaiton purposes.
[input_norm,~] = FECGESN_normalise_ecg(ref_ecg,1,5*fs);
[output_norm,~] = FECGESN_normalise_ecg(tar_ecg,1,5*fs);

if strcmp(metStruct.method,'LMS')
    % == least mean square
    if strcmp(metStruct.learningMode,'online')
        [y,residual,~] = FECGESN_lms(input_norm,output_norm,metStruct.mu,metStruct.Nunits);
    else
        % training readout on first 30sec and test on ALL signal using
        % constant weights
        
        % split into train and test
        train_fraction = TIME_INIT/(length(input_norm)/fs);
        [train_in,~] = split_train_test(input_norm,train_fraction);
        [train_out,~] = split_train_test(output_norm,train_fraction);
        
        % find lms param using training sequence
        [y,~,W] = FECGESN_lms(train_in,train_out,metStruct.mu,metStruct.Nunits);
        w = W(:,end);
        
        % now apply trained lms to all the signal
        residual = output_norm-filter(w,1,input_norm);
        
        % [~,w] = build_design_matrix_and_whopf(train_in,train_out,metStruct.Nunits);
        % chestf = conv(input_norm,w,'valid'); % idem as filtering
        % plot(output_norm(1:end)); hold on, plot(chestf(metStruct.Nunits+1:end),'r');
        % residual = output_norm-chestf(metStruct.Nunits+1:end);
    end
elseif strcmp(metStruct.method,'RLS')
    % == recursive least square
    p0 = 2*eye(metStruct.Nunits);
    RLSFilt = dsp.RLSFilter(metStruct.Nunits,'ForgettingFactor',metStruct.mu,'InitialInverseCovariance',p0);
    [y,residual] = RLSFilt(input_norm,output_norm);
end

end









