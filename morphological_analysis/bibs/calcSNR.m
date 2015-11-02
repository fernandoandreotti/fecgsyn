%  (mecg, fecg and noise).
%
% inputs
%       signal          signal to be used for energy calculation
%       mref            maternal QRS reference for SNR calculus
%       fref            fetal QRS reference for SNR calculus
%       fs              sampling frequency
% 
% output
%       SNRfm:      fetal-maternal SNR
%       SNRnf:      noise to fetal signal SNR
%
% FECG extraction toolbox, version 1.0, December 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Fernando Andreotti
% Dresden University of Technology, Institute of Biomedical Engineering
% fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 06-01-2015
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
function [SNRfm,SNRnf,Pmat,Pfet,Pn] = calcSNR(signal,mref,fref,fs)
%% Parameters
MWIN = 50;     % half of window around maternal complexes considered for signal energy (in ms)
FWIN = 25;     % half of window length around fetal complexes (in ms)
MHR = 60; %     [in bpm] constants to normalize SNR
FHR = 120; %    [in bpm]


%% Processing
MWIN = round(MWIN/(1000/fs));   % converting to samples
FWIN = round(FWIN/(1000/fs));   % at fs=1000Hz does nothing

% adjusting extremities
extremities = (mref <= MWIN | mref >= length(signal)-MWIN);        % if there are peaks on the border that may lead to error
mref = mref(~extremities);                                          % remove extremity peaks
extremities = (fref <= FWIN | fref >= length(signal)-FWIN);        % if there are peaks on the border that may lead to error
fref = fref(~extremities);                                          % remove extremity peaks

% selecting overlapping peaks
% fetal side
f_over = false(1,length(fref));
idxover = arrayfun(@(x) find(x - MWIN < fref & fref < x + MWIN),mref,'UniformOutput',false);  % looking for peaks overlapped peaks (within MWIN)
f_over(cell2mat(idxover)) = true;          % fet. peaks (indices) which are overlapped with maternal peaks (completely overlapped only)
f_nover = ~f_over;                      % fet. peaks indices where no overlap occur
% maternal side
m_over = false(1,length(mref));
idxover = arrayfun(@(x) find(x-MWIN < mref & mref<x+MWIN),fref,'UniformOutput',false);  % corresponding maternal peaks
m_over(cell2mat(idxover)) = true;           % mat. indices with overlap
m_nover = ~m_over;                       % mat. indices with no overlap

% defining limits for energy regions
sigtime = false(1,length(signal));
mat_sig = sigtime; mat_sig(mref(m_nover)) = true;          % creates a maternal allowed parts
mat_sig = logical(filter(ones(1,2*MWIN),1,mat_sig));       % maternal useful signal
mat_sig = [mat_sig(MWIN+1:end) false(1,MWIN)];               % taking moving average window delay

mat_nosig = sigtime; mat_nosig(mref(m_over)) = true;       % creates a maternal NOT allowed parts
mat_nosig = logical(filter(ones(1,2*MWIN),1,mat_nosig));   % maternal useless signal
mat_nosig = [mat_nosig(MWIN+1:end) false(1,MWIN)];               % taking moving average window delay

fet_sig = sigtime; fet_sig(fref(f_nover)) = true;          % creates a fetal allowed parts
fet_sig = logical(filter(ones(1,2*FWIN),1,fet_sig));       % fetal useful signal
fet_sig = [fet_sig(FWIN+1:end) false(1,FWIN)];               % taking moving average window delay

fet_nosig = sigtime; fet_nosig(fref(f_over)) = true;       % creates a maternal NOT allowed parts
fet_nosig = logical(filter(ones(1,2*FWIN),1,fet_nosig));   % maternal useless signal
fet_nosig = [fet_nosig(FWIN+1:end) false(1,FWIN)];           % taking moving average window delay

noise_sig = ~(mat_sig + mat_nosig + fet_sig + fet_nosig);

%% SNR Calculation

mbeats = 60*fs*length(mref)/length(signal); % now im bpm
fbeats = 60*fs*length(fref)/length(signal); % now im bpm

Pmat = sum(signal(:,mat_sig).^2,2)*(length(mref)/sum(m_nover)); % maternal power (in each channel)
Pmat = Pmat.*(MHR/mbeats);                          % normalized power
Pfet = sum(signal(:,fet_sig).^2,2)*(length(fref)/sum(f_nover)); % fetal power (in each channel)
Pfet = Pfet.*(FHR/fbeats);                          % normalized power
Pn   = sum(signal(:,noise_sig).^2,2); % noise power (in each channel)

SNRfm = 10*log10(Pfet./Pmat);
SNRnf= 10*log10(Pn./Pfet);

%% Plotting
figure
hold on
chan = 1; % channel to plot
maxes = max(abs(signal)');
t = 1:length(signal);
marea = area(t',maxes(chan)*[-mat_sig',2*mat_sig'],'FaceColor',[156 177 219]./255,'EdgeColor','none');
farea = area(t',maxes(chan)*[-fet_sig',2*fet_sig'],'FaceColor',[187 142 187]./255,'EdgeColor','none');
narea = area(t',maxes(chan)*[-noise_sig',2*noise_sig'],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
alpha = 0.1;
plot(signal(chan,:),'Color',[0.7, 0.7, 0.7])

xlabel('Time(s)','FontSize',14,'FontWeight','bold'), ylabel('Amplitude (mV)','FontSize',14,'FontWeight','bold')
set(gca,'XTick',1:fs:length(signal))  % This automatically sets
set(gca,'XTickLabel',num2cell(0:1:length(signal)/fs))
set(gca,'FontSize',12)

hsAnno = get(marea(1), 'Annotation');
hsLegend = get(hsAnno, 'LegendInformation');
set(hsLegend, 'IconDisplayStyle', 'off');
hsAnno = get(farea(1), 'Annotation');
hsLegend = get(hsAnno, 'LegendInformation');
set(hsLegend, 'IconDisplayStyle', 'off');
hsAnno = get(narea(1), 'Annotation');
hsLegend = get(hsAnno, 'LegendInformation');
set(hsLegend, 'IconDisplayStyle', 'off');

legend('mat.','fet.','noise')

% % for fast plot substituting logical to NaN
% mat_plot = NaN(1,length(mat_sig));    
% mat_plot(mat_sig) =  4*median(signal(chan,mat_sig));
% plot(mat_plot,'Color',[0 0.8 0]);plot(-mat_plot,'Color',[0 0.8 0]);       % maternal blocks
% 
% fet_plot = NaN(1,length(fet_sig));    
% fet_plot(fet_sig) =  4*median(signal(chan,fet_sig));
% plot(fet_plot,'r');plot(-fet_plot,'r');          % fetal blocks
% 
% noise_plot = NaN(1,length(noise_sig));    
% noise_plot(noise_sig) =  4*median(signal(chan,noise_sig));
% plot(noise_plot,'b');plot(-noise_plot,'b');          % noise blocks
% 
% legend('preproc. Signal', 'mat. Reference','fet. Reference')

