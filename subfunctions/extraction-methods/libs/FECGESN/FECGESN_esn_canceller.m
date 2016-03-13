function residual = FECGESN_esn_canceller(ref_ecg,tar_ecg,fs,esn,debug)
% ESN network using modified Frank Jaeger ESN toolbox and for its application
% to source separation for NI-FECG extraction. This function runs the ESN
% on the inputed signal ref_ecg to match the tar_ecg signal.
%
% inputs
%   ref_ecg:            reference ECG channels
%   tar_ecg:            target ECG channels
%   fs:                 sampling frequency
%   esn:                Echo State Network structure. The same structure should be used for all
%                       channels (by opposition to regenerating the network systematically for 
%                       each channel). See run_extraction for the structure
%                       parameters.
% ouput
%   residual:           return the residual containing the FECG on both
%                       train and test in a concatenated vector.
%
% NOTE: time used for initialising the weights = 30sec (see TIME_INIT constant)
%
% References:
%   1. ESN toolbox: http://reservoir-computing.org/node/129
%   2. Practical understanding of ESN:
%   http://www.scholarpedia.org/article/Echo_State_Network
%
% --
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
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

close all;

% == managing inputs
if nargin<4; error('esn_canceller: wrong number of input arguments \n'); end;
if size(ref_ecg,1)>size(ref_ecg,2); ref_ecg = ref_ecg'; end;
if size(tar_ecg,1)>size(tar_ecg,2); tar_ecg = tar_ecg'; end;

% == general
TIME_INIT = 30; % time used for initialising the weights

% == core function
try
    % == normalize the data
    % assumes that the first 1 to 5sec are representative (no huge noise) for normalizaiton purposes.   
    [input_norm,~]  = FECGESN_normalise_ecg(ref_ecg,1,5*fs);
    [output_norm,~] = FECGESN_normalise_ecg(tar_ecg,1,5*fs);
    
    train_frac = TIME_INIT/(length(input_norm)/fs);
    
    % == apply esn
    if strcmp(esn.learningMode,'online')
        % training readout using RLS (init but only 'training set')
        nForgetPoints = 50; % discard the first ... points
        [~,~,esn_out] = ESNTOOL_train_esn(input_norm,output_norm,esn,nForgetPoints,debug);
        residual = output_norm - esn_out;
    else
        % training readout using ridge regression or pseudo inverse (train/test sets)

        % split into train and test   
        [train_in,~] = ESNTOOL_split_train_test(input_norm,train_frac);
        [train_out,test_out] = ESNTOOL_split_train_test(output_norm,train_frac);

        % find esn param using training sequence
        nForgetPoints = 50 ; % discard the first ... points
        [trainedEsn, ~] = ESNTOOL_train_esn(train_in,train_out,esn,nForgetPoints,debug); 

        % now aplly trained esn to all the signal
        nForgetPoints = 0;
        esn_out = ESNTOOL_test_esn(input_norm,trainedEsn,nForgetPoints);
        residual = output_norm-esn_out;
    end
catch ME
    rethrow(ME);
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    residual = [];
end

% == display train and test results
if debug
    nb_of_points = length(ref_ecg);
    tm = 1/fs:1/fs:nb_of_points/fs;

    close all;
    % = plot normalized ref_ecg signals
    ax(1) = subplot(311); 
        plot(tm,input_norm(:,1),tm,output_norm,'r');
        legend('MECG chest','ABD');
    % = plot the ESN tar_ecg and the abdominal channel
    if strcmp(esn.learningMode,'online')
    ax(2) = subplot(312); 
        hold on, plot(tm,tar_ecg,tm,esn_out,'r');
        legend('ABD','ESN out'); ylabel('Amplitude');        
    else
    ax(2) = subplot(312); 
        hold on, plot(tm(train_frac*nb_of_points+1:end),test_out,tm,esn_out,'r');
        legend('ABD','ESN out'); ylabel('Amplitude');
    end
    ax(3) = subplot(313);
    % = plot residual
        plot(tm,residual);
        ylabel('Residual amplitude'); xlabel('Time [s]');    
    linkaxes(ax,'x'); 
    xlim([1 nb_of_points/fs]); % remove first sec beause init of ESN = big artefact
end

end













