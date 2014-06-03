%% main function for the fecg model
% Joachim Behar & Fernando Andreotti - December 2013/ January 2014

% this script presents the various physiological events modelled by the
% mecg-fecg model. For each example only one of these events is used. 
% This is in order to better highlight the effect of each individual event
% on the ECG morphology.

% TODO: Joachim - change the 'isempty' to 'any'
% FIXME: in load_param NEED TO MAKE SENSE OF THE VCG 2-4 for the ECTOPIC
% BEATS

%% == GLOBAL
clear all; close all; clc;
THR = 0.2; % threshold of QRS detector
mVCG = 5; % choose mother VCG (if empty then the simulator randomly choose one within the set of available VCGs)
fVCG = 4; % choose foetus VCG (ibid)
debug = 5; % debug level
CH_CANC = 5; % channel onto which to perform MECG cancellation
POS_DEV = 0; % slight deviation from default hearts and electrodes positions 
             % (0: hard coded values, 1: random deviation and phase initialisation)

%% == (1) SIMPLE RUN
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (1): SIMPLE RUN ----');
param.fs = 1000; % sampling frequency [Hz]

if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;
    
out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% == (2) ADDING NOISE
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (2): ADDING NOISE ----');
param.fs = 1000;
param.ntype = {'MA'}; % noise types
param.noise_fct = {1}; % constant SNR (each noise may be modulated by a function)
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

% % multiple noise sources
% param.ntype = {'MA' 'EM' 'BW'};     % noise types
% param.noise_fct = {1 1 1};          % constant SNR
 
% % multiple noise sources (varying noise power along measurement)
% % in order to do so, a modulating function should be given
% param.n = 20000;      % then number of samples should be defined here
% param.noise_fct{1} = 1+0.2*tanh(linspace(-pi,pi,param.n));  % tanh function
% param.noise_fct{2} = 1+sin(linspace(-2*pi,2*pi,param.n));   % sinus function
% param.noise_fct{3} = 0.1;                                   % e.g. of weighting noise power 

%% == (3) ADDING RESPIRATION
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (3): ADDING RESPIRATION ----');
param.fs = 1000;
param.mres = 0.25; % mother respiration frequency
param.fres = 0.8; % foetus respiration frequency
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% == (4) ADDING FOETAL MOVEMENT
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (4): ADDING FOETAL MOVEMENT ----');
param.fs = 1000;
debug=4;
param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% == (5) ADDING HEART RATE VARIABILITY
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (5): ADDING HEART RATE VARIABILITY ----');
param.fs = 1000;

% Case 4
param.fhr = 130; param.mhr = 130;

% Case 5
% param.macc = 40; % maternal acceleration in HR [bpm]
% param.mtypeacc = 'tanh';
% param.mhr = 110;
% param.fhr = 140;

if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% == (6) ADDING UTERINE CONTRACTION
% simulating uterus activity (through muscular noise) and heart rate changes for both 
% fetus (umbilical cord compression) and mother (acceleration followed by decelleration)
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (6): ADDING UTERINE CONTRACTION ----');
param.fs = 1000;
param.n = 60000;
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

% Case 6 (early deceleration)
x = linspace(-param.n/10,param.n/10,param.n);
mu = 0;

% Case 6b (late deceleration)
 mu = 0.5; % distance from center-beginning [%]
param.maccmean = -mu;

gauss = (100/(param.n*sqrt(2*pi)))*exp(-(x-(x(1)*mu)).^2/(2*(param.n/50)^2)); % approximating
gauss = gauss/max(gauss);                      % uterine contraction by gaussian modulated MA    

param.ntype = {'MA'};   
param.noise_fct = {gauss};
param.SNRmn = -10;

% heart rate changes
param.macc = 40;
param.mtypeacc = 'gauss';
param.facc = -30;
param.fhr = 130;
param.ftypeacc = {'mexhat'};
param.faccstd{1} = 0.5;
out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);


%% == (7) ADDING ECTOPIC BEATS
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (7): ADDING ECTOPIC BEATS ----');

param.fs = 1000; % sampling frequency [Hz]
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;
param.mectb = 1; param.fectb = 1; 

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% == (8) MULTIPLE PREGNANCIES (e.g twins)
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (8): MULTIPLE PREGNANCIES ----');

param.fs = 1000;
param.fhr = [120 150];
param.fres = [0 0];
param.facc = [0 0];
param.ftypeacc = {'none' 'none'};
if ~isempty(POS_DEV); param.posdev = 0; end;
if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = [fVCG (fVCG+3)]; end;

param.fheart{1} = [-pi/10 0.35 -0.1];
param.fheart{2} = [pi/10 0.4 -0.2];

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

%% random
% save2pdf('illustrateVolConductor');
% save2pdf('illustrateVolConductorUpperView');
% save2pdf('nstdb_noise_eg')
% save2pdf('generated_noise_eg')
% save2pdf('poles_ar_model')
% save2pdf('psd_nstdb_and_ar_generated')
% save2pdf('illustration_hrv')
% save2pdf('umbilicalcordcompression')
% save2pdf('mhr_tanh_increase')
% save2pdf('plot_respiration_waveform')
% save2pdf('ica_stationary')
% save2pdf('ica_nonstationary')
% save2pdf('multiple_preg_ica_8chan')
% save2pdf('multiple_preg_ica_5chan')
% save2pdf('TRratio')












