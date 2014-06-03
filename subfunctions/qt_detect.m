function out = qt_detect(ecg,qrs,fs,debug)
% this code is used to return the ECG QT interval given an input ecg with
% QRS locations. NOTE: this code was produced for DEMO purposes. It WAS NOT
% validated on a large database for the purpose of QT measurement.
%
% inputs
%   ecg:    one ecg channel
%   qrs:    qrs position [numbder of samples]
%   fs:     sampling frequency [Hz]
%   type:   'adult'; 'foetal' [string]
%   debug:  debug plots [bool]
%
% outputs
%   out: output structure containing
%       out.Qs:     Q-start
%       out.Te:     T-end
%       out.QT:     QT length
%       out.RTr:    T/R peak amplitude ratio
%
%
% reference
% [1] Vazquez-Seisdedos, Carlos R., et al. "New approach for T-wave end 
% detection on electrocardiogram: Performance in noisy conditions." 
% Biomedical engineering online 10.1 (2011): 1-11.
%
%
% Copyright (C) 2013  Joachim Behar & David Springer
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, david.springer@eng.ox.ac.uk
%
% Last updated : 05-02-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


% == check inputs
if nargin<3; fs = 1000; end; 
if nargin<4; debug = 0; end;
if size(ecg,2)<size(ecg,1); ecg=ecg'; end;

% == defines constants
HIGH_CUT_FREQ = 30;
LOW_CUT_FREQ = 0.5;
HR_MED = median(60./diff(qrs/fs));
C = 85/HR_MED; % will adjust the constant with respect to heart rate (85pbm is the reference for adults )
GAP100 = round(fs*0.1*C); 
GAP200 = round(fs*0.2*C);
GAP400 = round(fs*0.4*C); 
NB_QRS = length(qrs);

% == Pre-processing
[b_lp,a_lp] = butter(4,HIGH_CUT_FREQ/(fs/2),'high');
[b_bas,a_bas] = butter(4,LOW_CUT_FREQ/(fs/2),'high');
bpfecg_lp = ecg-filtfilt(b_lp,a_lp,ecg); % remove higher freq (zero phase)
bpfecg = filtfilt(b_bas,a_bas,bpfecg_lp); % remove baseline (zero phase)
ecgd = diff(bpfecg); % First derivative of input signal
SIGN_PEAKS = median(bpfecg(qrs));
abpfecg = abs(bpfecg);

% == core function
I = 2:NB_QRS-1; cpt=0; % have to account for borders
Te = zeros(NB_QRS-2,1,1);
Qs = zeros(NB_QRS-2,1,1);
if debug; Xms=zeros(NB_QRS-2,1,1); Xrs=zeros(NB_QRS-2,1,1); Tps=zeros(NB_QRS-2,1,1); end;
% == for each r-peak look for T-end
for rp = I
    try
    cpt=cpt+1;
    % == search for Q-start
    windq = bpfecg(qrs(rp)-GAP100:qrs(rp));
    if SIGN_PEAKS>0
        [~,ind_q] = min(windq);
    else
        [~,ind_q] = max(windq);
    end
    Qs(cpt) = ind_q+qrs(rp)-GAP100-1;
    
    % == search for T-end
    % = select working segment
    start = qrs(rp)+GAP100; % start of the windtow of interest
    stop = qrs(rp+1)-GAP200; % stop of the windtow of interest
    windt = bpfecg(start:stop); % select windtow to search cycle T-end
    lg_windt = length(windt);
    
    % = (1) Determination of the point identified as xm
    [~,Tp] = max(abs(windt(1:round(lg_windt/2)))); % top of T-wave indice       % - FIXME: where to look for top of T-wave?

    Tp = Tp-1;
    if Tp+GAP200>lg_windt; 
        stopXp = start+lg_windt-1; 
    else
        stopXp = start+Tp+GAP200; 
    end;
    [~,Xp] = max(abs(ecgd(start+Tp:stopXp))); % look for max of derivative
    Xp = Xp-1;
    xm = Tp+Xp;
    if debug; Xms(cpt) = start+xm; end;
    Tps(cpt) = start+Tp;
    
    % = (2) Determination of the point identified as xr
    str_iso = Tp+GAP100; % change from intial paper
    end_iso = Tp+GAP400;
    
    if end_iso>lg_windt; end_iso=lg_windt; end;
    
    [~,Xr] = min(abs(ecgd(str_iso:end_iso)));
    Xr = Xr-1;
    xr = Tp+GAP100+Xr; % change from intial paper
    if debug; Xrs(cpt) = start+xr; end; 
    
    % = (3) Calculation of the trapeziums areas
    xi = xm:xr;
    [~,old_index] = max(abs(0.5.*(windt(xm)-windt(xi)).*(2.*xr-xi-xm)));

    Te(cpt) = qrs(rp)+GAP100+xi(old_index);
    catch ME
        Qs(cpt) = 1;
        Te(cpt) = 1;
        if debug; Tps(cpt) = 1; end;
        if debug; Xms(cpt) = 1; end;
        if debug; Xrs(cpt) = 1; end; 
    end
end

% == format output
out.Qs = Qs;
out.Te = Te;
out.QT = Te-Qs;
out.RTr = ecg(Tps)./ecg(qrs(I)); % T/R peak amplitude ratio
out.qrs = qrs(I); % the QRS that were used

% == debug
if debug 
    LINE_WIDTH = 2;
    FONT_SIZE = 15;
    %tm = 1/fs:1/fs:NB_SAMP/fs;
    %close all; figure; 
    ax(1) = subplot(311); plot(ecg,'LineWidth',LINE_WIDTH);
    hold on, plot(Qs,ecg(Qs),'+k','LineWidth',LINE_WIDTH+3);
    hold on, plot(Te,ecg(Te),'+r','LineWidth',LINE_WIDTH+3);
    hold on, plot(Tps,ecg(Tps),'+m','LineWidth',LINE_WIDTH+3);
    title('ECG with T-end'); xlabel('Time [sec]');
    legend('ecg','Q','T','Tmax');
    linkaxes(ax,'x');
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);    
    
    ax(2) = subplot(312); plot(abpfecg,'b','LineWidth',LINE_WIDTH);
    hold on, plot(Qs,abpfecg(Qs),'+k','LineWidth',LINE_WIDTH+3);
    hold on, plot(Te,abpfecg(Te),'+r','LineWidth',LINE_WIDTH+3);
    hold on, plot(Tps,abpfecg(Tps),'+m','LineWidth',LINE_WIDTH+3);
    hold on, plot(qrs,abpfecg(qrs),'+g','LineWidth',LINE_WIDTH+3);
    title('abs of pre-filtered ECG'); xlabel('Time [sec]');
    legend('ecgfiltered','Q','T','Tmax','R-peak');
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    
    ax(3) = subplot(313); plot(ecgd,'b','LineWidth',LINE_WIDTH);
    hold on, plot(Xms,ecgd(Xms),'+r','LineWidth',LINE_WIDTH+3);
    hold on, plot(Xrs,ecgd(Xrs),'+k','LineWidth',LINE_WIDTH+3);
    title('first derivarive'); xlabel('Time [sec]');
    legend('ecg','xm','xr');
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    
    linkaxes(ax,'x');
end








