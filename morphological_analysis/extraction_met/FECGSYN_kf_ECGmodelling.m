function [OptimumParams,phase,ECGsd,w,wsd] = FECGSYN_kf_ECGmodelling(x,peaksidx,NbCycles,fs)
%ecg_filt = FECGSYN_kf_ECGfiltering(ecg,peaks,nbCycles,fs,debug);
%ECG FILTERING BLOCK Generates a model and call EKF/EKS
%   > Inputs
%      x:            data to be preprocessed
%      peaksidx:     maternal peak location
%      NbCycles:     number of cycles used to initialize template
%      fs:           sampling frequency [Hz]
%      flag:         flag = 0 to use EKF, if flag=1 to use EKS
%
%  > Output
%       Xhat:           filtered signal
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
% Last updated : 24-07-2014
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
global debug
%% Parameters
Nkernels = 7;          % number of kernels for Gaussian modelling
NB_BINS = 500;          % number of phase bins (SWT depends on this number!!)
scala = 6;              % number from scales used by SWT approach

% = parametrization
L = length(x);
Li = length(peaksidx);
peaks = zeros(1,L);
peaks(peaksidx(1:Li)) = 1;
fm = fs./diff(peaksidx);          % heart-rate
w = mean(2*pi*fm);          % average heart-rate in rads.
wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.


%% Generating average ECG beat
phase = PhaseCalc(find(peaks),length(x)); % phase calculation
phase_tmp = PhaseCalc(peaksidx(1:NbCycles),peaksidx(NbCycles)+300); % phase calculation
[ECGmean,ECGsd,meanphase] = ECG_tgen(x(1:peaksidx(NbCycles)+300),phase_tmp,NB_BINS); % mean ECG extraction

if debug>1
    % = phase calculation figure
    I = find(peaks);
    t = (0:length(x)-1)/fs;
    figure;
    plot(t,x*2*pi/max(x),'b');
    hold on
    plot(t(I),peaks(I)*2,'ro');
    plot(t,phase,'g','linewidth',2);
    grid;
    title('Phase Calculation');
    xlabel('time (sec.)');
    ylabel('Arbitrary units');
    legend('Scaled ECG','mECG references','Asigned phase');
    % phase Wrapping figure
    [X,Y,Z] = pol2cart(phase,1,x); %NOTE: pol2cart -> Transform polar or cylindrical coordinates to Cartesian.
    %pol2cart(THETA,RHO,Z) transforms the cylindrical coordinate data stored in corresponding
    %elements of THETA, RHO, and Z to three-dimensional Cartesian, or xyz coordinates.
    %The arrays THETA, RHO, and Z must be the same size (or any can be scalar).
    %The values in THETA must be in radians.
    figure;
    plot3(X,Y,Z);
    grid;
    title('Phase-wrapped ECG');
    xlabel('X (Arbitrary units)');
    ylabel('Y (Arbitrary units)');
    zlabel('ECG (Arbitrary units)');
end


%% Approximating average beat using Gaussians
% using Stationary Wavelet Transform approach (as in Andreotti CINC2014)
clear bi s gp scalex

% optimization procedure options
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'Display','off');

% loading quadratic spline coefficients
load('FECGSYN_kf_wtFilterCoefs.mat');

GaussPos = zeros(1,Nkernels);
Optpre = zeros(Nkernels,3);
ECGmean_aux = ECGmean;
samptorads = linspace(-pi,pi,length(ECGmean));
amax = max(abs(ECGmean));
amin = amax/100;
gp = zeros(1,Nkernels);
scalex = zeros(scala,NB_BINS);
bi = zeros(1,Nkernels);

if debug
    plot(meanphase,ECGmean,'Color',[0.7 0.7 0.7],'LineWidth',2)
    xlabel('Phase (rads.)');
    ylabel('Arbitrary units');
    hold on
end
% = iteratively obtaining candidate Gaussians
for i = 1:Nkernels
    scl = cell(1,scala); % initializing scales
    for k = 1:scala
        scalex(k,:) = conv(ECGmean_aux,wtLow{k},'same');
        lenc = length(wtLow{k});
        scl{k}=xcov(ECGmean_aux,[wtLow{k} zeros(1,NB_BINS-lenc)],'coeff');
        scl{k} = fliplr(scl{k});
        scl{k} = abs(scl{k}(ceil(lenc/2):NB_BINS+ceil(lenc/2)-1));
        scalex(k,:) = scalex(k,:);
    end
    scl = cell2mat(scl');
    [~,s]=max(max(scl')); % picking scale with highest cross-covariance
    [~,gp(i)] = max(scalex(s,:).^2);    % defining maximum
    
    %== Gaussian positioning constraint
    % avoids getting stuck on discontinuities
    if (i>1)&&(abs(gp(i)-gp(i-1))<= 5)
        [~,gp(i)] = max(scalex(randi([1 scala],1,1),:).^2);    % defining maximum
    end
    
    GaussPos(i) = samptorads(gp(i));  % converting to interval [-pi,pi]
    
    %== Gaussian standard deviation constraint
    % calculating expected standard deviation for Gaussians
    reswt = resample(wtLow{s},5,1); % just for finding stds
    hwhh = fwhm(1:length(reswt),reswt)/2/5;    % half width at half height
    bi(i) = abs(samptorads(1)-samptorads(ceil(hwhh)));
    
    % = small adjusts for fitting
    % extremities
    [~,idx] = min(abs(meanphase-GaussPos(i)));
    if idx <= 0
        idx = 1;
    elseif idx >= 500
        idx = 499;
    end
    
    % initial tetai in bins
    tetai = meanphase(idx);
    
    %== Amplitude directional constraint
    % initial alphai
    alphai = scalex(s,idx); % proposed initial point for gaussian amplitude
    if abs(amin) > abs(alphai)
        alphai = sign(alphai)*amin;
    end
    if alphai == 0
        alphai = amax;
    end
    
    
    %== Initial parameters for optimization
    InitParams = [alphai bi(i) tetai];
    
    %== Setting up bounds for optimization
    % alphai: between +/-amin and +/-amax
    % bi:     0.001 and estimated value + 1
    % tetai:  current value +/- pi/5
    LowBound = [min(amin*sign(alphai),amax*sign(alphai)),0.001,tetai-pi/5];
    UpBound = [max(amin*sign(alphai),amax*sign(alphai)),bi(i)+1,tetai+pi/5];
    OptimPar = lsqnonlin(@(InitParams) FECGSYN_kf_ECGModelError(InitParams,ECGmean_aux,meanphase),InitParams,LowBound,UpBound,options);
    
    %== Optimization procedure
    % Plot and Calculate average error in template
    [~,Model] = FECGSYN_kf_ECGModelError(OptimPar,ECGmean_aux,meanphase);
    Optpre(i,:) = OptimPar;
    ECGmean_aux = ECGmean_aux - Model;
    if debug
        % Plot intermediate results
        plot(meanphase,Model,'LineWidth',2)
    end
    
end

% = Re-run optimization procedure for definite solution
% optimization procedure options (more consuming than first)
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',Nkernels*100,'MaxFunEval',Nkernels*1000,'Display','off');
% adapting Gaussians coordinates to bins coordinates
[~,idx] = arrayfun(@(x)(min(abs(meanphase-GaussPos(x)))),1:length(GaussPos),'UniformOutput', false);
idx = cell2mat(idx);
tetai = meanphase(idx);     %closest values for teta (based on GaussPos that belongs to meanphase)

% purposed initial point for gaussian amplitude
alphai = 1.2*ECGmean(idx);

% initial bi's are sugested by SWT procedure
InitParams = [alphai bi tetai];

% setting bounds
% alphai: 50% +/-
% bi:     0.001 and estimated value + 1
% tetai:  current value +/- pi/2
LowBound = zeros(size(InitParams));
UpBound = zeros(size(InitParams));
LowBound(1:Nkernels) = alphai - amax/4;
LowBound(Nkernels+1:2*Nkernels) = repmat(0.001,1,Nkernels);
LowBound(2*Nkernels+1:end) = tetai-pi/2;
idx2low = (LowBound(2*Nkernels+1:end)<-pi); % checking if not too low
if any(idx2low)
    LowBound([false(1,2*Nkernels) idx2low]) = -pi;
end
UpBound(1:Nkernels) = alphai + amax/4;
UpBound(Nkernels+1:2*Nkernels) = bi+1;
UpBound(2*Nkernels+1:end) = tetai +pi/2;
idx2high = (UpBound(2*Nkernels+1:end)>pi); % checking if not too low
if any(idx2high)
    UpBound([false(1,2*Nkernels) idx2high]) = pi;
end


% = Final Gaussian fitting
OptimumParams = lsqnonlin(@(InitParams) FECGSYN_kf_ECGModelError(InitParams,ECGmean,meanphase),InitParams,LowBound,UpBound,options);

clear L Li LowBound Model OptimPar ECGmean_aux


%% Remove invalid Gaussians
% Look if there is any invalid Gaussian and remove it
% If gaussian width or height is less than 0.001 it will excluded.
N = Nkernels;
rmidx = (abs(OptimumParams(1:N))<0.001.*amax)|(abs(OptimumParams(N+1:2*N))<0.001);
if any(rmidx)
    OptimumParams([rmidx rmidx rmidx]) = [];
    disp('Throwing Gaussians away!')
end
N = length(OptimumParams)/3;     %new number of Gaussian kernels


% Plot final resultsNew Folder
if debug && ~isempty(OptimumParams)
    figure('units','normalized','outerposition',[0 0 1 1])
    [Error,Model] = FECGSYN_kf_ECGModelError(OptimumParams,ECGmean,meanphase);
    errorbar(meanphase,ECGmean,ECGsd/2);
    hold on;
    plot(meanphase,ECGmean,'r');
    plot(meanphase,Model,'m','linewidth',2)
    plot(meanphase,Error,'-g','linewidth',2);
    plot(OptimumParams(2*N:end),1,'xr','LineWidth',2)
    legend('SD bar','Mean ECG','Gaussian Approx.','Error','Kernel position');%,'Gaussian');
    title(['N_k = ' num2str(N)]);
    xlabel('Phase (rads.)');
    ylabel('Arbitrary units');
    grid
    hold off;
end

end

%% Phase Calculation
%
% This function generates a saw-tooth phase signal, based on QRS locations
%
% > Inputs
%       peaks:          fidutials location, in these points phase is zero
%       NbSamples:      number of samples on the signal
%
% This function is based on Dr. Sameni's OSET
%
%
function phase = PhaseCalc(peaks,NbSamples)
phase = zeros(1,NbSamples);
m = diff(peaks);            % gets distance between peaks
% = dealing with borders (first and last peaks may not be full waves)

% first interval uses second interval as reference
L = peaks(1);       %length of first interval
if isempty(m)       % only ONE peak was detected
    phase(1:NbSamples) = linspace(-2*pi,2*pi,NbSamples);
else
    phase(1:L) = linspace(2*pi-L*2*pi/m(1),2*pi,L);
    % beats in the middle
    for i = 1:length(peaks)-1;      % generate phases between 0 and 2pi for almos all peaks
        phase(peaks(i):peaks(i+1)) = linspace(0,2*pi,m(i)+1);
    end                             % 2pi is overlapped by 0 on every loop
    % last interval
    % uses second last interval as reference
    L = length(phase)-peaks(end);   %length of last interval
    phase(peaks(end):end) = linspace(0,L*2*pi/m(end),L+1);
end
phase = mod(phase,2*pi);
phase(phase>pi) = phase(phase>pi)- 2*pi;
end


%% Calculation of the mean and SD of ECG waveforms in different beats
%
% > Inputs
%       ecg:            input ECG signal
%       phase:          ECG phase
%       NB_BINS:        number of desired phase bins
%
% > Outputs
%       ECGmean:        mean ECG beat
%       ECGsd:          standard deviation of ECG beats
%       meanPhase:      the corresponding phase for one ECG beat
%
%
% Although this function structure is based on OSET's toolbox, the
% averaging procedure itself is based on Dr. Oster's approach for stacking
% and averaging beats.
function [ECGmean,ECGsd,meanPhase] = ECG_tgen(ecg,phase,NB_BINS)

ini_cycles = find(phase(2:end)<0&phase(1:end-1)>0)+1; % start of cycles
cycle_len = diff(ini_cycles); % distance between cycles
end_cycles = ini_cycles(1:end-1)+cycle_len-1; % start of cycles
meanPhase = linspace(-pi,pi,NB_BINS);
% stacking cycles
cycle = arrayfun(@(x) interp1(phase(ini_cycles(x):end_cycles(x)),...
    ecg(1,ini_cycles(x):end_cycles(x)),meanPhase,'spline'),...
    1:length(ini_cycles)-1,'UniformOutput',0);
cycle = cell2mat(cycle');
ECGmean = mean(cycle);
ECGsd = std(cycle);

end

