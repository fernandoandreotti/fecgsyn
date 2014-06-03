function bwa = generate_breathing(fs,N,fres,debug)
% variable sawtooth shape for modelling respiration. The intention is to
% model the respiration pattern in order to modulate the orientation 
% of the cardiac dipole with breathing.
%
% reference: http://mathworld.wolfram.com/FourierSeriesSawtoothWave.html
% Only the first three coefficients of the Fourier Transform of the
% sawtooth function are kept to make the bwa smooth.
%
% inputs
%   fs:     sampling frequency [Hz]
%   N:      number of datapoints
%   fres:   respiratory frequency [Hz]
%   debug:  [bool]
%
% output
%   bwa:    normalised breathing waveform (range [-0.5: 0.5])
% 
%
% reference:
% [1] Petrenas et al. "An Echo State Neural Network for QRST Cancellation During Atrial
% Fibrillation". IEEE Trans Biomed. Eng , VOL. 59, NO. 10, OCTOBER 2012.
% [description for generating the breathing waveform for rotation matrix modulation]
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 19-01-2014
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

deltaA = 0.3;   % modulation amplitude
fa = 0.1;       % amplitude modulation frequency
deltaF = 0.05;  % maximum frequency deviation
ff = 0.1;       % modulation frequency

a = @(n) 2./pi*(1+deltaA*sin(2*pi*fa.*n/fs));
b = @(n) 2./(2*pi)*(1+deltaA*sin(2*pi*fa.*n/fs));
c = @(n) 2./(3*pi)*(1+deltaA*sin(2*pi*fa.*n/fs));

s = @(n) 0.5 - (a(n).*sin(2*pi*fres.*n/fs + deltaF/ff*sin(2*pi*ff*n/fs))...
    + b(n).*sin(4*pi*fres.*n/fs + deltaF/ff*sin(2*pi*ff*n/fs))...
    + c(n).*sin(6*pi*fres.*n/fs + deltaF/ff*sin(2*pi*ff*n/fs)));

% == add random noise to vary the bwa
signal_noisy = s(1:N) + 0.1*rand(1,N);

% == then cubic spline to resmooth
signal_noisy_interp = interp1(1:1:N,signal_noisy,1:100:N);

% == reinterpolate to go back to initial fs
bwa = resample(signal_noisy_interp, 100,1);

% == normalise in [-0.5 0.5]
bwa = bwa/(abs(max(bwa))+abs(min(bwa)));
bwa = bwa + abs(min(bwa)) - 0.5;

% == debug
if debug
    FONT_SIZE = 15;
    figure('name','respiration waveform');
    tm = 1/fs:1/fs:N/fs;
    plot(tm,bwa,'LineWidth',2);
    xlabel('Time [sec]'); ylabel('Amplitude [Normalised]');
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
end

end

