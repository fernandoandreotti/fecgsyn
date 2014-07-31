function [P1dB,P2dB] = plot_psd(y_1,y_2,fs,method,leg)
% plots two normalised power spectral density (PSD) on the same plot for 
% analysis. Two methods can be specified for computing the PSD
% (welch or burg). This function was designed for allowing direct
% comparison of two PSD (e.g. y_1 as the raw ecg and y_2 as the filtered ecg)  
% and also in order to look into the difference in PSD evaluation using 
% both the welch and the burg methods.
%
% inputs
%   y_1: signal 1 for which you want to compute the PSD 
%   y_2: signal 2 for which you want to compute the PSD 
%   fs:  sampling frequency
%   method: 'welch' or 'burg'
%
% outputs
%   P1dB: normalised PSD of signal 1
%   P1dB: normalised PSD of signal 2
% 
% IMPORTANT NOTE: the parameters of the welch and the burg functions are
% hard-coded below. You need to play with them for your application.
%
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
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

% == use welch of burg to compute PSD
if strcmp(method,'welch')
    [P1 F1] = pwelch(y_1,5*fs,fs,10000,fs,'onesided'); % 10000 - 5000
    [P2 F2] = pwelch(y_2,5*fs,fs,10000,fs,'onesided');
elseif strcmp(method,'burg')
    [P1 F1] = pburg(y_1,20,[],fs);
    [P2 F2] = pburg(y_2,20,[],fs);
end

% == convert to db
P1dB = 10*log10(P1/(mean(P1)));
P2dB = 10*log10(P2/(mean(P2)));

% == plots everything
hold all;
plot(F1,P1dB,'--r','LineWidth',2); 
plot(F2,P2dB,'LineWidth',2);
allAxesInFigure = findall(gcf,'type','axes'); set(allAxesInFigure,'fontsize',16);

xlabel('Frequency [Hz]','fontsize',16,'color','k');
ylabel('Normalized power [db]','fontsize',16,'color','k');
legend(leg{1},leg{2});
hold off; grid on;

end

