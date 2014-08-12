% Data = FECGx_kf_ECGfiltering(x,fs,flag,X0,P0,peaks,debug)
%
% This function calls the Kalman Filter/Smoother to extract the fetal ECG.
% It is an adaptation of the work developed by Reza Sameni named OSET Toolbox
% available at: http://www.oset.ir/
%
% Inputs:
% x:        Data to be filtered
% peaks: Maternal references (Optional)
% NbCycles: how many cycles (on the begining) will be used to initialize
%           model
% fs:       Sampling frequency [in Hz]
% flag:     If flag = 0 use EKF, if flag=1 use EKS [bool]
% debug:    Debugging function
% 
% 
% Output:
%   xfiltered: Processed signal
%   X0: Initial estimate por x (used for next segment processed)
%   P0: Initial estimate por templatelen (used for next segment processed)
%
% Fetal Extraction Toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Fernando Andreotti
% Dresden University of Technology, Institute of Biomedical Engineering
% fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 24-07-2014
%
% Based on: Synthetic ECG model error
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
function [xfiltered,X0,P0] = FECGx_kf_ECGfiltering(x,peaks,NbCycles,fs,flag,debug)

%//////////////////////////////////////////////////////////////////////////
% MEAN PHASE EXTRACTION BLOCK
%//////////////////////////////////////////////////////////////////////////
      
% == finding maternal peaks within NbCycles first cycles
peak_tmp=peaks(peaks<=peaks(NbCycles));      
bins = fs/2;  % number of phase bins    
phase_tmp = FECGx_kf_PhaseCalc(peak_tmp,peak_tmp(end)); % phase calculation temporary
[ECGmean,ECGsd,meanphase] = FECGx_kf_MeanECGExtraction(x,phase_tmp,bins,1); % mean ECG extraction
phase = FECGx_kf_PhaseCalc(peaks,length(x)); % phase calculation

%//////////////////////////////////////////////////////////////////////////
% PLOTS AND RESULTS
%//////////////////////////////////////////////////////////////////////////
if debug
    % Phase calculation figure
    
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
    
    
    % Phase Wrapping figure
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
% PARAMETER EXTRACTION
%//////////////////////////////////////////////////////////////////////////

fm = fs./diff(peaks);          % heart-rate
w = mean(2*pi*fm);          % average heart-rate in rads.
wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.


%% Optimizing Gaussian Kernels to mean ECG (as in Andreotti CINC2014)

%parameters
Nkernels = 10; % number of Gaussian kernels to be used
scala = 7;  % number from scales used


GaussPos = zeros(1,Nkernels);
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'Display','off');  %Optimization options
Optpre = zeros(Nkernels,3);
ECGmean_aux = ECGmean;
samptorads = linspace(-pi,pi,length(ECGmean));
amax = max(abs(ECGmean));
amin = amax/100;
rms = zeros(1,Nkernels);

load('FECGx_kf_wtFilterCoefs.mat'); % loading quadratic spline coeff
for i = 1:Nkernels
    %     %+++ using wavelets
    scalex(1,:) = conv(ECGmean_aux,wtLow{1},'same');
    for k = 2:scala
        scalex(k,:) = conv(scalex(k-1,:),wtLow{k},'same');
    end

    % optimization settings
    [s,GaussPos(i)] = find(scalex.^2 == max(max(scalex.^2)));
    GaussPos(i) = samptorads(GaussPos(i));  % converting to interval [-pi,pi]
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
        alphai = amin;
    end
    bi = 0.04;  % proposed initial point for gaussian width

    InitParams = [alphai bi tetai]; % initial parameters for optimization
    LowBound = [min(amin*sign(alphai),amax*sign(alphai)),0.001,tetai-pi/5];
    UpBound = [max(amin*sign(alphai),amax*sign(alphai)),bi+1,tetai+pi/5];
    OptimPar = lsqnonlin(@(InitParams) FECGx_kf_ECGModelError(InitParams,ECGmean_aux,meanphase),InitParams,LowBound,UpBound,options); %Optimization
        
    % Plot and Calculate average error in template
    [Error,Model] = FECGx_kf_ECGModelError(OptimPar,ECGmean_aux,meanphase);
    rms(i) = sqrt(sum(Error.^2)/500);
    % plot
    if debug
        figure(1)
        subplot(1,2,1)
        h = findobj('type','line');
        delete(h)
        waterfall(scalex)
        hold on
        plot3(GaussPos(i),s,scalex(s,idx),'rv')
        hold off
        subplot(1,2,2)
        hold on;
        plot(meanphase,ECGmean,'r');
        plot(meanphase,Model,'m','linewidth',2)
        plot(meanphase,Error,'-g','linewidth',2);
        legend('Mean ECG','Model','Error');%,'Gaussian');
        legend('Mean ECG','Gaussian Approx.','Error');%,'Gaussian');
        title('Mean, SD extraction and Model');
        xlabel('Phase (rads.)');
        ylabel('Arbitrary units');
        %         frame = getframe(figure(1));
        %         writeVideo(writerObj,frame);
        hold off;
    end
    Optpre(i,:) = OptimPar;
    ECGmean_aux = ECGmean_aux - Model;
    
    
end
OptimParam = reshape(Optpre,1,3*Nkernels);
[modelerror,model] = FECGx_kf_ECGModelError(OptimParam,ECGmean,meanphase);

% Remove invalid Gaussians
% Look if there is any invalid Gaussian and remove it
% If gaussian width or height is less than 0.001 it will excluded.
N = Nkernels;
yy = 1;
while(yy<N+1)
    if((abs(OptimParam(N+yy))<0.001)||(abs(OptimParam(yy))<0.001))
        OptimParam(2*N+yy)=[];OptimParam(N+yy)=[];OptimParam(yy)=[];
        N = N-1;
    else
        yy=yy+1;
    end
end
Nkernels = length(OptimParam)/3;     %new number of Gaussian kernels

% Plot final results
if debug && ~isempty(Optpos)
    figure(2)
    [Error,Model] = ECGModelError(OptimParam,ECGmean,meanphase);
    %     errorbar(meanphase,ECGmean,ECGsd/2);
    hold on;
    plot(meanphase,ECGmean,'r');
    plot(meanphase,Model,'m','linewidth',2)
    plot(meanphase,Error,'-g','linewidth',2);
    plot(Optpos(:,end),1,'or','MarkerSize',7)
    legend('SD Bar','Mean ECG','Model','Error');%,'Gaussian');
    legend('Mean ECG','Gaussian Approx.','Error');%,'Gaussian');
    title('Mean, SD extraction and Model');
    xlabel('Phase (rads.)');
    ylabel('Arbitrary units');
    grid
    hold off;
end

%% Optimizing Gaussian Kernels to mean ECG (Random Search)
% % nIter = 100;
% % Nkernels = 10; % number of Gaussian kernels to be used
% % 
% % [OptimParam,model,modelerror] = FECGSYN_kf_gaussfit(ECGmean,nIter,Nkernels,debug);
% % Nkernels = length(OptimParam)/3;
% % 
% % 

%% Kalman Filter Parametrization

% matrix of observation signals (samples x 2). First column corresponds
% to the phase observations and the second column corresponds to the noisy
% ECG
y = [phase ; x];
X0 = [-pi;0];       % initialize phase
P0 = [(2*pi)^2 0 ;0 (10*max(abs(x))).^2]; % initialize error cov. mat.
%covariance matrix of the process noise vector
% should train these gains (grid search)
GQ = 0.01;
Q = diag( [(.1*OptimParam(1:Nkernels)).^2 (.1*ones(1,Nkernels)).^2 (.1*ones(1,Nkernels)).^2 (0.1*wsd)^2 , GQ*var(modelerror)]);

%covariance matrix of the observation noise vector
GR = 1;
R = [(w/fs).^2 0 ;0 GR*mean(ECGsd(1:round(length(ECGsd)))).^2];

%mean process noise vector
Wmean = [OptimParam w 0]';

%mean observation noise vector
Vmean = [0 0]';

u = zeros(1,length(x));

%//////////////////////////////////////////////////////////////////////////
% ECG PROCESSING
%//////////////////////////////////////////////////////////////////////////
%Use EKF or EKS
[Xfiltered,X0,P0,e] = FECGx_kf_EKFilter(y,X0,P0,Q,R,Wmean,Vmean,OptimParam,w,fs,flag,u);
xfiltered = Xfiltered(2,:);
end


%% Mean ECG Extraction
% This function calculates the average ECG beat
%
% [ECGmean,ECGsd,meanPhase] = MeanECGExtraction(x,phase,bins,flag)
% Calculation of the mean and SD of ECG waveforms in different beats
%
% inputs:
% x: input ECG signal
% phase: ECG phase
% bins: number of desired phase bins
% flag
%     1: aligns the baseline on zero, by using the mean of the first 10%
%     segment of the calculated mean ECG beat
%     0: no baseline alignment
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

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
function [ECGmean,ECGsd,meanPhase] = FECGx_kf_MeanECGExtraction(x,phase,bins,flag)

meanPhase = zeros(1,bins);
ECGmean = zeros(1,bins);
ECGsd = zeros(1,bins);

I = find( phase>=(pi-pi/bins)  | phase<(-pi+pi/bins) );
if length(I)>1
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x(I));
    ECGsd(1) = std(x(I));
else
    meanPhase(1) = 0;
    ECGmean(1) =0;
    ECGsd(1) = -1;
end
for i = 1 : bins-1;
    I = find( phase >= 2*pi*(i-0.5)/bins-pi & phase < 2*pi*(i+0.5)/bins-pi );
    if(~isempty(I))
        meanPhase(i + 1) = mean(phase(I));
        ECGmean(i + 1) = mean(x(I));
        ECGsd(i + 1) = std(x(I));
    else
        meanPhase(i + 1) = 0;
        ECGmean(i + 1) = 0;
        ECGsd(i + 1) = -1;
    end
end
K = find(ECGsd==-1);
for i = 1:length(K),
    switch K(i)
        case 1
            meanPhase(K(i)) = -pi;
            ECGmean(K(i)) = ECGmean(K(i)+1);
            ECGsd(K(i)) = ECGsd(K(i)+1);
        case bins
            meanPhase(K(i)) = pi;
            ECGmean(K(i)) = ECGmean(K(i)-1);
            ECGsd(K(i)) = ECGsd(K(i)-1);
        otherwise
            meanPhase(K(i)) = mean([meanPhase(K(i)-1),meanPhase(K(i)+1)]);
            ECGmean(K(i)) = mean([ECGmean(K(i)-1),ECGmean(K(i)+1)]);
            ECGsd(K(i)) = mean([ECGsd(K(i)-1),ECGsd(K(i)+1)]);
    end
end

if(flag==1)
    ECGmean = ECGmean - mean(ECGmean(1:ceil(length(ECGmean)/10)));
end

end