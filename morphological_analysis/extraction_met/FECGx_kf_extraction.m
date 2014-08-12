function residual = FECGx_kf_extraction(peaks,ecg,method,debug,varargin)
% MECG cancellation algorithms using the Extended Kalman Filter/Smoother.
% The code is based on the PhD from Reza Sameni and the code provided in
% OSET Toolbox (http://www.oset.ir/).
%
% Inputs
%   peaks:      MQRS markers in ms. Each marker corresponds to the
%               position of a MQRS
%   ecg:        matrix of abdominal ecg channels
%   method:     method to use (TS,TS-CERUTTI,TS-SUZANNA,TS-LP,TS-PCA)
%   varargin:
%       nbCycles:   number of cycles to use in order to build the mean MECG template
%       fs:         sampling frequency (NOTE: this code is meant to work at 1kHz)
%
% output
%   residual:   residual containing the FECG
%
%
% Adapted from:
% Fetal Extraction Toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Fernando Andreotti
% Dresden University of Technology, Institute of Biomedical Engineering
% fernando.andreotti@mailbox.tu-dresden.de
% Available at: http://fernando.planetarium.com.br
%
% Current version:
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

% == manage inputs
optargs = {20 1000};  % default values for [nbCycles fs]
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);
[nbCycles,fs] = optargs{:};

if nargin > 6
    error('kf_extraction: too many input arguments \n');
end

% check that we have more peaks than nbCycles
if nbCycles>length(peaks)
    error('MECGcancellation Error: more peaks than number of cycles for average ecg');
end

% check which method to use
switch method
    case 'EKF'
        flag = 0;
    case 'EKS'
        flag = 1;
    otherwise
        error('kf_extraction: method not implemented.')
end

% == MECG estimation using KF
% try   
    ecg_filt = FECGx_kf_ECGfiltering(ecg,peaks,nbCycles,fs,flag,debug);
    % == compute residual
    residual = ecg - ecg_filt;
    
% catch ME
%     for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
%     residual = ecg;
% end

% == debug
if debug
    FONT_SIZE = 15;
    tm = 1/fs:1/fs:length(residual)/fs;
    figure('name','MECG cancellation');
    plot(tm,ecg,'LineWidth',3);
    hold on, plot(tm,ecg-residual,'--k','LineWidth',3);
    hold on, plot(tm,residual-1.5,'--r','LineWidth',3);
    hold on, plot(tm(peaks),ecg(peaks),'+r','LineWidth',2);
    
    legend('mixture','template','residual','MQRS');
    title('Template subtraction for extracting the FECG');
    xlabel('Time [sec]'); ylabel('Amplitude [NU]')
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end

end






















