%% script_qt.m
% 
%
% This script was used to produce the example in the paper related to morphological analysis. 
% The T/R ratio was computed while considering various amount of MA noise
% added to the FECG mixture.
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

%% == GLOBAL
clear all; close all; clc;
THR = 0.2; % threshold of QRS detector
mVCG = 5; % choose mother VCG (if empty then the simulator randomly choose one within the set of available VCGs)
fVCG = 4; % choose foetus VCG (ibid)
debug = 0; % debug level
CH_CANC = 5; % channel onto which to perform MECG cancellation
POS_DEV = 0; % slight deviation from default hearts and electrodes positions 
             % (0: hard coded values, 1: random deviation and phase initialisation)
             
%% == (1) SIMPLE RUN
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (1): SIMPLE RUN ----');
param.fs = 1000; % sampling frequency [Hz]
param.n = 60000;

if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;
    
out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

out = qt_detect(res,qrs_det,param.fs,debug);

NB_MEAS = length(out.RTr);
PACE = 5;
NB_POINTS = floor(NB_MEAS/PACE);
bRTr = zeros(NB_POINTS,1); bqrs = zeros(NB_POINTS,1);

for ii=1:NB_POINTS
   bRTr(ii) = median(out.RTr(PACE*(ii-1)+1:PACE*ii));
   bqrs(ii) = median(out.qrs(PACE*(ii-1)+1:PACE*ii));
end
res1 = res;

%% == (2) ADDING NOISE
close all; clear param; clear res; clear out; clear cmqrs; clear qrs_det;
disp('---- Example (2): ADDING NOISE ----');

param.fs = 1000;
param.SNRmn = 20; % signal to noise ratio between mother ECG power and background noise
param.ntype = {'MA'}; % noise types
param.noise_fct = {1}; % constant SNR (each noise may be modulated by a function)
param.n = 60000;

if ~isempty(mVCG); param.mvcg = mVCG; end;
if ~isempty(fVCG); param.fvcg = fVCG; end;
if ~isempty(POS_DEV); param.posdev = 0; end;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

out = qt_detect(res,qrs_det,param.fs,debug);

NB_MEAS = length(out.RTr);
PACE = 5;
NB_POINTS = floor(NB_MEAS/PACE);
nRTr = zeros(NB_POINTS,1); nqrs = zeros(NB_POINTS,1);

for ii=1:NB_POINTS
   nRTr(ii) = median(out.RTr(PACE*(ii-1)+1:PACE*ii));
   nqrs(ii) = median(out.qrs(PACE*(ii-1)+1:PACE*ii));
end

res2 = res;

% == more noise
param.SNRmn = 15;

out = run_ecg_generator(param,debug);
cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
res = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI');
[qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

out = qt_detect(res,qrs_det,param.fs,debug);

NB_MEAS = length(out.RTr);
PACE = 5;
NB_POINTS = floor(NB_MEAS/PACE);
nnRTr = zeros(NB_POINTS,1); nnqrs = zeros(NB_POINTS,1);

for ii=1:NB_POINTS
   nnRTr(ii) = median(out.RTr(PACE*(ii-1)+1:PACE*ii));
   nnqrs(ii) = median(out.qrs(PACE*(ii-1)+1:PACE*ii));
end

res3 = res;

% == plots
close all;
FONT_SIZE = 20;
LINE_WIDTH = 3;
plot(bqrs/param.fs,bRTr,'b','LineWidth',LINE_WIDTH); 
hold on, plot(nqrs/param.fs,nRTr,'r','LineWidth',LINE_WIDTH);
hold on, plot(nnqrs/param.fs,nnRTr,'k','LineWidth',LINE_WIDTH);
xlabel('Time [sec]'); ylabel('T/R [NU]');
legend('No noise','with noise SNRmn=20','with noise SNRmn=15');
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
set(gca,'FontSize',FONT_SIZE);

tm = 1/param.fs:1/param.fs:length(res1)/param.fs;
figure; plot(tm,res1,'LineWidth',LINE_WIDTH); 
hold on, plot(tm,res2,'r','LineWidth',LINE_WIDTH);
hold on, plot(tm,res3,'k','LineWidth',LINE_WIDTH);
legend('residual with no noise','residual with SNRmn=20','residual with SNRmn=15');
set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
set(gca,'FontSize',FONT_SIZE);


