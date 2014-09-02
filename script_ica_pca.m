%% script for ICA & PCA simulations
%
% this script generates two abdominal mixture and apply PCA and ICA on the
% array of generated channels. The purpose is to illustrate how well PCA
% and ICA are performing in i) the stationary case (i.e stationary mixing matrix)
% and ii) non-stationary case (when adding breathing effects, foetal movement
% etc.). A second simulation shows the results of ICA when considering
% single and twin pregnancy.
%
%
% fecgsyn toolbox, version 1.0, July 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 31-08-2014
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

%% == (1) BSS examples (single pregnancy)

% = general
clear all; close all;
param.fvcg = 3; % force to keep same vcgs for the stationary 
                % and non-stationary cases in order to make a direct
                % comparison
param.mvcg = 5;
param.fs = 1000;
param.n = 60000; % number of data points to generate
param.SNRfm = -30; % signal to noise ratio between foetus and mother ECG 
debug = 0;
tm = 1/param.fs:1/param.fs:param.n/param.fs;

% = (a) stationary case
param.mtypeacc = 'none'; % force constant mother heart rate
param.ftypeacc = {'none'}; % force constant foetal heart rate
param.posdev = 0; % no random deviation from default hearts and electrodes positions

out_st = run_ecg_generator(param,debug); % stationary output
out_st.mixture = out_st.mixture(1:6,:); % just select the first 8 channels because ICA computation takes too long with the default 32

% = (b) non-stationary case
param.macc = 35; % maternal acceleration in HR [bpm]
param.mtypeacc = 'tanh'; % hyperbolic tangent acceleration
param.facc = -34; % foetal decceleration in HR [bpm]
param.ftypeacc = {'mexhat'}; % gaussian drop and recovery
param.mres = 0.25; % mother respiration frequency
param.fres = 0.8; % foetus respiration frequency
param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
param.SNRmn = -5; % signal to noise ratio between mother ECG power and added noise power
param.ntype = {'MA','BW'}; % noise types
param.noise_fct = {1,1}; % constant SNR (each noise may be modulated by a function)

out_nst = run_ecg_generator(param,debug); % non-stationary output
out_nst.mixture = out_nst.mixture(1:6,:);

% = normalise data
stmix_trans = bsxfun(@minus,out_st.mixture,mean(out_st.mixture,2)); % remove mean
stmix = bsxfun(@rdivide,stmix_trans,std(out_st.mixture,0,2)); % divide by standard deviation

nstmix_trans = bsxfun(@minus,out_nst.mixture,mean(out_nst.mixture,2)); % remove mean
nstmix = bsxfun(@rdivide,nstmix_trans,std(out_nst.mixture,0,2)); % divide by standard deviation

LINE_WIDTH = 2;
FONT_SIZE = 24;

% = applying PCA
close all;
[~,src_st] = princomp(stmix');
[~,src_nst] = princomp(nstmix');
NB_CH = size(src_st,2); axst = zeros(8,1); axnst = zeros(8,1);
if NB_CH>6; NB_CH = 6; end; % because otherwise too many channels on one single plot
figure('name','PCA stationary case');
figure('name','PCA non-stationary case');
for cc=1:NB_CH
    figure(1); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
    figure(2); axnst(cc) = subplot(NB_CH,1,cc); plot(tm,src_nst(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);    
end
linkaxes(axst,'x');
linkaxes(axnst,'x');
xlabel('Time [sec]');
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);
figure(1); 
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);
figure(2); 
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);

% = applying ICA
close all;
[~,src_st] = jade(ceil(1000*stmix)); % FIXME: NOT SURE WHY NEED TO MULTIPLY BY 1000
[~,src_nst] = jade(ceil(1000*nstmix)); % FIXME: NOT SURE WHY NEED TO MULTIPLY BY 1000
NB_CH = size(src_st,2); axst = zeros(8,1); axnst = zeros(8,1);
figure('name','ICA stationary case');
figure('name','ICA non-stationary case');
for cc=1:NB_CH
    figure(1); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    figure(2); axnst(cc) = subplot(NB_CH,1,cc); plot(tm,src_nst(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end
linkaxes(axst,'x');
linkaxes(axnst,'x');
figure(1); 
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);
figure(2); 
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);


%% == (2) BSS examples (multiple pregnancy)
% this script is meant to illustrate the limitation of ICA when the number
% of components representing the sources is higher than the number of
% channels that are recorded. A twin pregnancy is simulated and ICA run on
% 8 abdominal channels and then 5 abdominal channels.

% = general
clear all; close all;
param.fvcg = [3,4]; % force to keep same vcgs for the stationary 
                % and non-stationary cases in order to make a direct
                % comparison
param.mvcg = 5;
param.fheart{1} = [-pi/10 0.4 -0.3];
param.fheart{2} = [-pi/1.1 0.3 -0.1];
param.fhr = [120 150]; % foetuses mean heart rate
param.posdev = 0; % no random deviation from default hearts and electrodes positions
param.fs = 1000;
debug = 0;
param.n = 10000;
tm = 1/param.fs:1/param.fs:param.n/param.fs;

% = stationary case
param.mtypeacc = 'none'; % force constant mother heart rate
param.ftypeacc = {'none','none'}; % force constant foetal heart rate
param.facc = [0 0]; % foetuses accelerations on mean HR (no acceleration)

out_st = run_ecg_generator(param,debug); % stationary output
out_st.mixture = out_st.mixture(1:6,:); % just select the first 8 channels because ICA computation takes too long with the default 32

% = normalise data
stmix_trans = bsxfun(@minus,out_st.mixture,mean(out_st.mixture,2)); % remove mean
stmix = bsxfun(@rdivide,stmix_trans,std(out_st.mixture,0,2)); % divide by standard deviation

LINE_WIDTH = 2;
FONT_SIZE = 15;

% = applying ICA
close all;
% = 6 channels
[~,src_st] = jade(ceil(1000*stmix));
NB_CH = size(src_st,2); axst = zeros(8,1);
figure('name','ICA stationary case - 6 channels');
for cc=1:NB_CH
    figure(1); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end
linkaxes(axst,'x');
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);

% = 4 channels
[~,src_st_red] = jade(ceil(1000*stmix(1:4,:)));
NB_CH = size(src_st_red,2); axst = zeros(8,1);
figure('name','ICA stationary case - 4 channels');
for cc=1:NB_CH
    figure(2); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st_red(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end
linkaxes(axst,'x');
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);






