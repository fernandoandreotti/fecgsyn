function [ecgme,ecgsd,phaseme,nbcyc] = FECGSYN_tgen(ecg,qrs,phase,bins,fs,flag,debug)
% this function is used to contruct a template ECG. The steps for building.
% Different methods were tested.
% The procedure for building the template is:
% 1. identify the beats that have similar RR
% 2. apply Reza's phase wrapping between [-pi pi]
% Note: the median is used to build the template accross the cyles. This
% naturally remove the outliers noisy cycles and turned to work better than
% using corcoeff.
%
%
% inputs
%   ecg:     ecg signal
%   qrs:     qrs fiducials
%   phase:   phase of the EGC
%   bins:    number of bins for the phase wrapping
%   fs:      sampling frequency
%   flag:    adjust 'iso' to zero
%   debug:   debug level (1/2)
%
% outputs
%   ecgme:   template ECG
%   ecgsd:   ECG standard deviation used as a proxy of signal quality 
%            (i.e. trust of the Kalman in observations)
%   phaseme: mean phase of the template ECG
%   nbcyc:   number of cycles used to build the template
%
%
% Dual EKF, version 1.0, March 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 07-07-2014
% This function is adapted from the OSET toolbox of Dr Reza Sameni
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

% == constants
qrs_cut = qrs(2:end-2); % remove first and last QRS to account for border effects
RR = diff(qrs_cut)/fs;
RRm = median(RR);
TOL_RR = 0.015; % cycles that are TOL_RR different from the min RR are excluded to prevent innacurate linear phase stretching
NB_CYC = length(qrs_cut); % number of cycles
WINDOW = round((RRm/2)*fs); % window either way of the R-peak

% == baseline (make sure baseline wander has been cancelled)
LF_CUT = 0.7;
[bs,as] = butter(3,LF_CUT/(fs/2),'high');
ecg = filtfilt(bs,as,ecg);

% == only select beat with similar RR interval (i.e. HR)
indrr = 0;
while sum(indrr)<ceil(NB_CYC/1.5)
    indrr = abs(RR-RRm)<TOL_RR;
    TOL_RR = TOL_RR + 0.005; % be more tolerant if not enough beats
end
indrr = ~indrr; % the one to remove
rposrem = qrs_cut(indrr); % position of the R-peaks of the cycles to remove

% == now build the template
phaseme = zeros(1,bins);
ecgme = zeros(1,bins);
ecgsd = zeros(1,bins);

I = find( phase>=(pi-pi/bins)  | phase<(-pi+pi/bins) );
if ~isempty(rposrem)
    [~,D] = dsearchn(rposrem',I');
    I = I(D>1.2*WINDOW);
end

if(~isempty(I))
    phaseme(1) = -pi;
    ecgme(1) = median(ecg(I));
    ecgsd(1) = std(ecg(I));
else
    phaseme(1) = 0;
    ecgme(1) = 0;
    ecgsd(1) = -1;
end

for i = 1 : bins-1;
    I = find( phase >= 2*pi*(i-0.5)/bins-pi & phase < 2*pi*(i+0.5)/bins-pi );
    if ~isempty(rposrem)
        [~,D] = dsearchn(rposrem',I');
        I = I(D>1.2*WINDOW);
    end
   
    if(~isempty(I))
        phaseme(i + 1) = median(phase(I));
        ecgme(i + 1) = median(ecg(I));
        ecgsd(i + 1) = std(ecg(I)); 
    else
        phaseme(i + 1) = 0;
        ecgme(i + 1) = 0;
        ecgsd(i + 1) = -1;
    end
end

K = find(ecgsd==-1);
for i = 1:length(K),
    switch K(i)
        case 1
            phaseme(K(i)) = -pi;
            ecgme(K(i)) = ecgme(K(i)+1);
            ecgsd(K(i)) = ecgsd(K(i)+1);
        case bins
            phaseme(K(i)) = pi;
            ecgme(K(i)) = ecgme(K(i)-1);
            ecgsd(K(i)) = ecgsd(K(i)-1);
        otherwise
            phaseme(K(i)) = mean([phaseme(K(i)-1),phaseme(K(i)+1)]);
            ecgme(K(i)) = mean([ecgme(K(i)-1),ecgme(K(i)+1)]);
            ecgsd(K(i)) = mean([ecgsd(K(i)-1),ecgsd(K(i)+1)]);
    end
end

if(flag==1)
    ecgme = ecgme - mean(ecgme(1:ceil(length(ecgme)/10)));
end

nbcyc = NB_CYC-sum(indrr);

% == debug
if debug>=1
    fprintf('Number of cycles detected: %f \n',length(qrs_cut));
    fprintf('Number of cycles selected: %f \n',NB_CYC-sum(indrr));
    if isnan(ecgme); disp('WARNING: Not able to produce ecgme, low correlation'); end  
elseif debug>=2
    LINEWIDTH = 2;
    FONTSIZE = 20;
    ph = linspace(-pi,pi,bins);
    plot(ph,ecgme,'LineWidth',LINEWIDTH);
    xlim([-pi pi]);
    xlabel('Phase [rad]');
    ylabel('Amplitude [NU]');
    set(findall(gcf,'type','text'),'fontSize',FONTSIZE);
    set(gca,'FontSize',FONTSIZE);
end





