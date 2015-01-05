function Xhat = FECGx_kf_ECGfiltering(x,peaksidx,NbCycles,fs)
%ecg_filt = FECGx_kf_ECGfiltering(ecg,peaks,nbCycles,fs,debug);
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
global debug GR GQ

%% Parameters
Nkernels = 10;          % number of kernels for Gaussian modelling
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

if debug
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
load('FECGx_kf_wtFilterCoefs.mat');

GaussPos = zeros(1,Nkernels);
Optpre = zeros(Nkernels,3);
ECGmean_aux = ECGmean;
samptorads = linspace(-pi,pi,length(ECGmean));
amax = max(abs(ECGmean));
amin = amax/100;
gp = zeros(1,Nkernels);
scalex = zeros(scala,NB_BINS);
bi = zeros(1,Nkernels);

% = iteratively obtaining candidate Gaussians
for i = 1:Nkernels
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
    [~,s]=max(max(c2')); % picking scale with highest cross-covariance
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
    InitParams = [alphai bi(i) tetai]; % initial parameters for optimization
    LowBound = [min(amin*sign(alphai),amax*sign(alphai)),0.001,tetai-pi/5];
    UpBound = [max(amin*sign(alphai),amax*sign(alphai)),bi(i)+1,tetai+pi/5];
    OptimPar = lsqnonlin(@(InitParams) FECGx_kf_ECGModelError(InitParams,ECGmean_aux,meanphase),InitParams,LowBound,UpBound,options);
    %Optimization
    % Plot and Calculate average error in template
    [~,Model] = FECGx_kf_ECGModelError(OptimPar,ECGmean_aux,meanphase);
    Optpre(i,:) = OptimPar;
    ECGmean_aux = ECGmean_aux - Model;
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

% = Final Gaussian fitting
OptimumParams = lsqnonlin(@(InitParams) FECGx_kf_ECGModelError(InitParams,ECGmean,meanphase),InitParams,InitParams-2,InitParams+2,options);

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

% Plot final resultsNew Folder
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

% = ECG
y = [phase ; x];

% = covariance matrix of the process noise vector
Q = diag( [(0.1*OptimumParams(1:N)).^2 (0.5*ones(1,N)).^2 (0.5*ones(1,N)).^2 wsd^2 , GQ*(0.5*mean(ECGsd))^2]);

% = covariance matrix of the observation noise vector
R = diag([(w/fs).^2/12      GR*mean(ECGsd).^2]);

% = covariance matrix for error
P0 = diag([(2*pi)^2,(10*max(abs(x))).^2]); % error covariance matrix

% = noises
Wmean = [OptimumParams w 0]';
Vmean = [0 0]'; % mean observation noise vector

% = initialize state
X0 = [-pi 0]';  % state initialization

% = control input
u = zeros(1,length(x));


%% Filtering
Xhat = FECGx_kf_EKFilter(y,X0,P0,Q,R,Wmean,Vmean,OptimumParams,w,fs,flag,u);

end

function phase = PhaseCalc(peaks,NbSamples)
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

function [ECGmean,ECGsd,meanPhase] = ECG_tgen(ecg,phase,NB_BINS)
% Calculation of the mean and SD of ECG waveforms in different beats
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

