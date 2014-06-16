%% Script for testing ICA & PCA simulations for morphological analysis
% 
% this script generates a series of abdominal mixtures and apply PCA and ICA on the
% array of generated channels. The purpose is to illustrate how well PCA
% and ICA are performing in i) the stationary case (i.e stationary mixing matrix)
% and ii) non-stationary case (when adding breathing effects, foetal movement
% etc.). As one of the main assumption behind the blind souce separation
% method (in their classical forms) is the stationarity of the mixing
% matrix it is expected that they will fail in separating the FECG from the
% MECG in the second case ii). The analysis should be performed on both
% spatial and temporal (back propagated) signals.
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
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

% clear all; close all;
debug = 0;

%% Generating simulated data with various SNR
param.n = 60000; % number of data points to generate
param.fs = 1000;
param.SNRfm = -20;

for SNRmn = -20:0
    param.SNRmn = SNRmn;    % varying SNRmn
    param.fhr = randi([110,160],1,1);   % choosing fhr in normal interval
    param.mhr = randi([60,100],1,1);    % choosing maternal heart rate
    % stationary mixture
    param.mtypeacc = 'none'; % force constant mother heart rate
    param.ftypeacc = {'none'}; % force constant foetal heart rate
    param.posdev = 0; % no random deviation from default hearts and electrodes positions
    out_st = run_ecg_generator(param,debug); % stationary output
    plotmix(out_st)
    
    % non-stationary mixture
    param.macc = 20; % maternal acceleration in HR [bpm]
    param.mtypeacc = 'tanh'; % hyperbolic tangent acceleration
    param.facc = -40; % foetal decceleration in HR [bpm]
    param.ftypeacc = {'mexhat'}; % gaussian drop and recovery
    param.mres = 0.25; % mother respiration frequency
    param.fres = 0.8; % foetus respiration frequency
    param.ftraj{1} = 'spline'; % giving spiral-like movement to fetus
    param.ntype = {'MA','EM','BW'}; % noise types
    param.noise_fct = {1,1,1}; % constant SNR (each noise may be modulated by a function)
    
    out_nst = run_ecg_generator(param,debug); % non-stationary output
    plotmix(out_nst)
    
end
