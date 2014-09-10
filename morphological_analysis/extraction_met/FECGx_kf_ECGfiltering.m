function [ output,Data] = FECGx_kf_ECGfiltering(x,fs,flag,indicepeaks,debug)
%ECG FILTERING BLOCK Generates a model and call EKF/EKS
%   INPUTS:
%
%   x: data with preprocessing
%   fs: Sampling frequency
%   flag: If flag = 0 use EKF, if flag=1 use EKS
%   X0: Initial estimate por x
%   P0: Initial estimate por templatelen
%   indicepeaks: Maternal references (Optional)
%
%   OUTPUTS:
%   output: Processed signal
%   X0: Initial estimate por x (used for next segment processed)
%   P0: Initial estimate por templatelen (used for next segment processed)

%//////////////////////////////////////////////////////////////////////////
% PEAK DETECTION BLOCK
%//////////////////////////////////////////////////////////////////////////
Nkernels = 10;

% If detection provided is provided load annotations
L = length(x);
Li = length(indicepeaks);
peaks = zeros(1,L);
peaks(indicepeaks(1:Li)) = 1;


%//////////////////////////////////////////////////////////////////////////
% MEAN PHASE EXTRACTION BLOCK
%//////////////////////////////////////////////////////////////////////////

phase = PhaseCalc(find(peaks),length(x)); % phase calculation
%Data.phase{end+1} = phase;
NB_BINS = 300;          % number of phase bins
[ECGmean,ECGsd,meanphase] = meanbeat(x,phase,NB_BINS); % mean ECG extraction

%//////////////////////////////////////////////////////////////////////////
% PLOTS AND RESULTS
%//////////////////////////////////////////////////////////////////////////
if debug
    %% Phase calculation figure
    
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
    
    
    %% Phase Wrapping figure
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

%//////////////////////////////////////////////////////////////////////////
%% PARAMETER EXTRACTION
%//////////////////////////////////////////////////////////////////////////

fm = fs./diff(indicepeaks);          % heart-rate
w = mean(2*pi*fm);          % average heart-rate in rads.
wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.
scala = 6;  % number from scales used

GaussPos = zeros(1,Nkernels);
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'Display','off');  %Optimization options

Optpre = zeros(Nkernels,3);
ECGmean_aux = ECGmean;
samptorads = linspace(-pi,pi,length(ECGmean));

amax = max(abs(ECGmean));
amin = amax/100;

% Stationary Wavelet Transform (as in Andreotti CINC2014)
load('FECGx_kf_wtFilterCoefs.mat'); % loading quadratic spline coeff
clear bi s gp scalex
gp = zeros(1,Nkernels);
scalex = zeros(scala,NB_BINS);
timeoptim = [];

for i = 1:Nkernels
    c1 = cell(1,scala);
    c2 = cell(1,scala);
    for k = 1:scala
        scalex(k,:) = conv(ECGmean_aux,wtLow{k},'same');
        lenc = length(wtLow{k});
        c2{k}=xcov(ECGmean_aux,[wtLow{k} zeros(1,NB_BINS-lenc)],'coeff');
        c2{k} = fliplr(c2{k});
        c2{k} = abs(c2{k}(ceil(lenc/2):NB_BINS+ceil(lenc/2)-1));
        scalex(k,:) = scalex(k,:);%./mean(scalex(k,:).^2);
    end
        
    c2 = cell2mat(c2');
%     figure
%     surf(c2','EdgeColor','none')
%     axis square
%     xlabel('Scale (j)')
%     ylabel('Samples [n]')
%     zlabel('Normalized Cross-Covariance (NU)')
%     
    [~,s]=max(max(c2')); % picking scale with biggest xcov
    [~,gp(i)] = max(scalex(s,:).^2);
    GaussPos(i) = samptorads(gp(i));  % converting to interval [-pi,pi]
    
    % calculating expected width
    reswt = resample(wtLow{s},5,1); % just for finding std
    hwhh = fwhm(1:length(reswt),reswt)/2/5;    % half width at half height
    bi(i) = abs(samptorads(1)-samptorads(ceil(hwhh)));
    
    [~,idx] = min(abs(meanphase-GaussPos(i)));
    % extremities gausspos
    if idx <= 0
        idx = 1;
    elseif idx >= 500
        idx = 499;
    end
    tetai = meanphase(idx);
    alphai = mean(ECGmean_aux(idx)); % proposed initial point for gaussian amplitude
    if abs(amin) > abs(alphai)
        alphai = sign(alphai)*amin;
    end
    if alphai == 0
        alphai = amax;
    end
    tstart = tic;
    InitParams = [alphai bi(i) tetai]; % initial parameters for optimization
    LowBound = [min(amin*sign(alphai),amax*sign(alphai)),0.001,tetai-pi/5];
    UpBound = [max(amin*sign(alphai),amax*sign(alphai)),bi(i)+1,tetai+pi/5];
    OptimPar = lsqnonlin(@(InitParams) FECGx_ECGModelError(InitParams,ECGmean_aux,meanphase),InitParams,LowBound,UpBound,options); 
%Optimization
    timeoptim(end+1)=toc(tstart);
    % Plot and Calculate average error in template
    [Error,Model] = FECGx_ECGModelError(OptimPar,ECGmean_aux,meanphase);
   
    gspot(i,:) = Model;
    Optpre(i,:) = OptimPar;
    ECGmean_aux = ECGmean_aux - Model;
end
OptimumParams = reshape(Optpre,1,3*Nkernels);
[Error,Model] = FECGx_ECGModelError(OptimumParams,ECGmean,meanphase);

%%%%%%%%%%%%%%%%
%% SWT 4 %%%%%%%
%%%%%%%%%%%%%%%%
% Wavelet Positions and allow Optimization
% Given Wavelet Positions
tswt6=tic;
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',Nkernels*100,'MaxFunEval',Nkernels*1000,'Display','off');  %Optimization options
[~,idx] = arrayfun(@(x)(min(abs(meanphase-GaussPos(x)))),1:length(GaussPos),'UniformOutput', false); 
   % looks for minimal distance from GaussPos and meanphase
idx = cell2mat(idx);        %converts cell to array (both are numeric)

tetai = meanphase(idx); %closest values for teta (based on GaussPos that belongs to meanphase)
alphai = 1.2*ECGmean(idx); % purposed initial point for gaussian amplitude
% bi = .04*ones(size(alphai)); % purposed initial point for gaussian width
InitParams = [alphai bi tetai];
OptimumParams = lsqnonlin(@(InitParams) FECGx_ECGModelError(InitParams,ECGmean,meanphase),InitParams,InitParams-2,InitParams+2,options); 
%Optimization
[Error,Model] = FECGx_ECGModelError(OptimumParams,ECGmean,meanphase);

clear L Li LowBound Model OptimPar ECGmean_aux


%% Remove invalid Gaussians
% Look if there is any invalid Gaussian and remove it
% If gaussian width or height is less than 0.001 it will excluded.
N = Nkernels;
yy = 1;
while(yy<N+1)
    if((abs(OptimumParams(N+yy))<0.001)||(abs(OptimumParams(yy))<0.001))
        OptimumParams(2*N+yy)=[];OptimumParams(N+yy)=[];OptimumParams(yy)=[];
        N = N-1;
        disp('Throwing Gaussians away!')
    else
        yy=yy+1;
    end
end
N = length(OptimumParams)/3;     %new number of Gaussian kernels

% Plot final results
% if debug && ~isempty(Optpos)
%     figure(2)
%     %     errorbar(meanphase,ECGmean,ECGsd/2);
%     hold on;
%     plot(meanphase,ECGmean,'r');
%     plot(meanphase,Model,'m','linewidth',2)
%     plot(meanphase,Error,'-g','linewidth',2);
%     plot(Optpos(:,end),1,'or','MarkerSize',7)
%     legend('SD Bar','Mean ECG','Model','Error');%,'Gaussian');
%     legend('Mean ECG','Gaussian Approx.','Error');%,'Gaussian');
%     title('Mean, SD extraction and Model');
%     xlabel('Phase (rads.)');
%     ylabel('Arbitrary units');
%     grid
%     hold off;
% end

%% Kalman Filter Parametrization
% [~,Model] = FECGx_ECGModelError(OptimumParams,ECGmean,meanphase);
% Error2 = sum((ECGmean-Model).^2)/sum((ECGmean-mean(ECGmean)).^2); % Normalized Mean Square Error (%)
% matrix of observation signals (samples x 2). First column corresponds
% to the phase observations and the second column corresponds to the noisy
% ECG
y = [phase ; x];

%covariance matrix of the process noise vector
% should train these gains (grid search)
GQ = 0.1;
% Q = diag( [(.1*OptimumParams(1:N)).^2 (.1*ones(1,N)).^2 (.1*ones(1,N)).^2 (0.1*wsd)^2 , GQ*var(Error)]);
Q = diag( [(.1*OptimumParams(1:N)).^2 (.1*ones(1,N)).^2 (.1*ones(1,N)).^2 (0.1*wsd)^2 , (GQ*mean(ECGsd))^2]);

%covariance matrix of the observation noise vector
GR = 1;
R = [(w/fs).^2 0 ;0 GR*mean(ECGsd(1:round(length(ECGsd)))).^2];
Wmean = [OptimumParams w 0]';

%mean observation noise vector
Vmean = [0 0]';
X0 = [-pi 0]';
P0 = [(2*pi)^2 0 ;0 (10*max(abs(x))).^2];
u = zeros(1,length(x));

%//////////////////////////////////////////////////////////////////////////
%% ECG PROCESSING
%//////////////////////////////////////////////////////////////////////////
%Use EKF or EKS
% disp('Parameters estimated. Filtering...')
Xfiltered = FECGx_EKFilter(y,X0,P0,Q,R,Wmean,Vmean,OptimumParams,w,fs,flag,u);                          
output = Xfiltered(2,:);
% plot(output,'k')
end

%% Phase Calculation
function phase = PhaseCalc(peaks,lengthx)
% Based on Sameni's method for phase calculation
phase = zeros(1,lengthx);
m = diff(peaks);            % gets distance between peaks
% first interval
% dealing with borders (first and last peaks may not be full waves)
% uses second interval as reference
L = peaks(1);   %length of first interval
if isempty(m) % only ONE peak was detected
    phase(1:lengthx) = linspace(-2*pi,2*pi,lengthx);
else
    phase(1:L) = linspace(2*pi-L*2*pi/m(1),2*pi,L);
    % beats in the middle
    for i = 1:length(peaks)-1;  % generate phases between 0 and 2pi for almos all peaks
        phase(peaks(i):peaks(i+1)) = linspace(0,2*pi,m(i)+1);
    end                                         % 2pi is overlapped by 0 on every loop
    % last interval
    % uses second last interval as reference
    L = length(phase)-peaks(end);   %length of last interval
    phase(peaks(end):end) = linspace(0,L*2*pi/m(end),L+1);
end
phase = mod(phase,2*pi);
phase(find(phase>pi)) = phase(find(phase>pi))- 2*pi;
end

function [ECGmean,ECGsd,meanPhase] = meanbeat(ecg,phase,NB_BINS)
%
% [ECGmean,ECGsd,meanPhase] = MeanECGExtraction(x,phase,bins,flag)
% Calculation of the mean and SD of ECG waveforms in different beats
%
% inputs:
% x: input ECG signal
% phase: ECG phase
% NB_BINS: number of desired phase bins
%
% outputs:
% ECGmean: mean ECG beat
% ECGsd: standard deviation of ECG beats
% meanPhase: the corresponding phase for one ECG beat
%
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
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
%
%
% == using Julien's approach
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

