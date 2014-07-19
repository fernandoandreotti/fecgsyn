% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

% Testing that run_scenario is working as desired
clear all; close all; clc;

choice = 10;

THR = 0.2; % threshold of QRS detector
mVCG = 5; % choose mother VCG (if empty then the simulator randomly choose one within the set of available VCGs)
fVCG = 4; % choose foetus VCG (ibid)
debug = 0; % debug level
CH_CANC = 5; % channel onto which to perform MECG cancellation
POS_DEV = 0; % slight deviation from default hearts and electrodes positions 
             % (0: hard coded values, 1: random deviation and phase initialisation)

[out, param] = eval_scenario( choice, THR, mVCG, fVCG, debug, CH_CANC, POS_DEV);
%[out, param] = eval_scenario( choice );

cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);