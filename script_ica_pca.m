%% script for ICA & PCA simulations
% Joachim Behar & Fernando Andreotti - December 2013/ January 2014

% this script generates two abdominal mixture and apply PCA and ICA on the
% array of generated channels. The purpose is to illustrate how well PCA
% and ICA are performing in i) the stationary case (i.e stationary mixing matrix)
% and ii) non-stationary case (when adding breathing effects, foetal movement
% etc.). As one of the main assumption behind the blind souce separation
% method (in their classical forms) is the stationarity of the mixing
% matrix it is expected that they will fail in separating the FECG from the
% MECG in the second case ii).

%% == (1) BSS examples (single pregnancy)

% = general
clear all; close all;
param.fvcg = 3; % force to keep same vcgs for the stationary 
                % and non-stationary cases in order to make a direct
                % comparison
param.mvcg = 5;
param.fs = 1000;
param.n = 60000; % number of data points to generate
param.SNRfm = -20; % signal to noise ratio between foetus and mother ECG 
debug = 0;
tm = 1/param.fs:1/param.fs:param.n/param.fs;

% = (a) stationary case
param.mtypeacc = 'none'; % force constant mother heart rate
param.ftypeacc = {'none'}; % force constant foetal heart rate
param.posdev = 0; % no random deviation from default hearts and electrodes positions

out_st = run_ecg_generator(param,debug); % stationary output

% = (b) non-stationary case
param.macc = 20; % maternal acceleration in HR [bpm]
param.mtypeacc = 'tanh'; % hyperbolic tangent acceleration
param.facc = -40; % foetal decceleration in HR [bpm]
param.ftypeacc = {'mexhat'}; % gaussian drop and recovery
param.mres = 0.25; % mother respiration frequency
param.fres = 0.8; % foetus respiration frequency
param.ftraj{1} = 'spline'; % giving spiral-like movement to fetus
param.SNRmn = -2; % signal to noise ratio between mother ECG power and added noise power
param.ntype = {'MA','EM','BW'}; % noise types
param.noise_fct = {1,1,1}; % constant SNR (each noise may be modulated by a function)

out_nst = run_ecg_generator(param,debug); % non-stationary output

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

% = stats for each channel
% in the stationary case
THR = 0.2;
F1 = zeros(NB_CH,1); Se = zeros(NB_CH,1); PPV = zeros(NB_CH,1);
for cc=1:NB_CH
    [qrs_det,~,~] = qrs_detect(src_st(:,cc),THR,0.150,param.fs,[],[],debug);
    [F1(cc),Se(cc),PPV(cc),~] =stats(out_st.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out_st.param.n/param.fs,param.fs);
end
[~,ind] = max(F1);
fprintf('Best Score in the stationary case is: Se=%f, PPV=%f, F1=%f \n', Se(ind), PPV(ind), F1(ind));

% in the non-stationary case
THR = 0.2;
F1 = zeros(NB_CH,1); Se = zeros(NB_CH,1); PPV = zeros(NB_CH,1);
for cc=1:NB_CH
    [qrs_det,~,~] = qrs_detect(src_nst(:,cc),THR,0.150,param.fs,[],[],debug);
    [F1(cc),Se(cc),PPV(cc),~] =stats(out_nst.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out_nst.param.n/param.fs,param.fs);
end
[~,ind] = max(F1);
fprintf('Best Score in the stationary case is: Se=%f, PPV=%f, F1=%f \n', Se(ind), PPV(ind), F1(ind));

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
debug = 5;
param.n = 10000;
tm = 1/param.fs:1/param.fs:param.n/param.fs;

% = stationary case
param.mtypeacc = 'none'; % force constant mother heart rate
param.ftypeacc = {'none','none'}; % force constant foetal heart rate
param.facc = [0 0]; % foetuses accelerations on mean HR (no acceleration)

out_st = run_ecg_generator(param,debug); % stationary output

% = normalise data
stmix_trans = bsxfun(@minus,out_st.mixture,mean(out_st.mixture,2)); % remove mean
stmix = bsxfun(@rdivide,stmix_trans,std(out_st.mixture,0,2)); % divide by standard deviation

LINE_WIDTH = 2;
FONT_SIZE = 15;

% = applying ICA
close all;
% = 8 channels
[~,src_st] = jade(ceil(1000*stmix));
NB_CH = size(src_st,2); axst = zeros(8,1);
figure('name','ICA stationary case - 8 channels');
for cc=1:NB_CH
    figure(1); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end
linkaxes(axst,'x');
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);

% = 5 channels
[~,src_st_red] = jade(ceil(1000*stmix(1:5,:)));
NB_CH = size(src_st_red,2); axst = zeros(8,1);
figure('name','ICA stationary case - 5 channels');
for cc=1:NB_CH
    figure(2); axst(cc) = subplot(NB_CH,1,cc); plot(tm,src_st_red(:,cc),'LineWidth',LINE_WIDTH); xlim([5 10]);
    set(gca,'FontSize',FONT_SIZE); set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end
linkaxes(axst,'x');
xlabel('Time [sec]','fontSize',FONT_SIZE);
ylabel('Amplitude [NU]','fontSize',FONT_SIZE);






