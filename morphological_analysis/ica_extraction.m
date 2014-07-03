function [F1,RMS,PPV,SE] = ica_extraction(data,FS,chan,refqrs,varargin)
% Uses Independent Component Analysis (ICA) over selected channels input 
% data and choses best channel based on th F1 measure.
% 
% Input
% data:      Matrix containing signals to serve as input to ICA.
% FS:        Sampling frequency [Hz]
% chan:      Channels from data to be used in ICA
% refqrs:    Array containing reference QRS detections for F1 measure
% (optional)
% blength    Iterates ICA every blength (in seconds)
% 
% Output
% bestsig:   Selected unmixed signal
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 30-06-2014
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

%% Input test
if size(data,1) > size(data,2)
    data = data';
end

data = data(chan,:);    % only using selected channels

switch length(varargin)
    case 0
        blength = size(data,2);
    case 1
        blength = varargin{1};
    otherwise
        error('ica_extraction: too many inputs given to function')
end
% Parameters
INTERV = round(0.05*FS);    % BxB acceptance interval
TH = 0.3;   % detector threshold
REFRAC = round(.15*FS)/1000; % detector refractory period

%% Run ICA on every block
blength = blength * FS;
ssamp = 1;          % starting sample to filter (offset)
endsamp = ssamp + blength - 1;      % ending sample to filter
loop = 1;           % allows iterations
icamorph = zeros(1,length(data));      % best produced ICA channel
qrsica = icamorph;                  % resulting qrs detections
if size(data,2)<endsamp
    blength = Data.length;
    endsamp = size(data,2);
end

% = normalise data
trans = bsxfun(@minus,data,mean(data,2)); % remove mean
data = bsxfun(@rdivide,trans,std(data,0,2)); % divide by standard deviation
icasignal = zeros(size(data));

while (loop)  %quit will be given as soon as complete signal is filtered
    if (size(data,2) - ssamp) < 1.5*blength        % if there is less than 1.5x wrapping
        endsamp = size(data,2);              % interval, it should filter all
        loop = 0;                           % and be the last loop iteration
    end
    samp2filt = ssamp:endsamp;              % creating a list with samples to filter
     idx = refqrs>=ssamp & refqrs <= endsamp;
     refint = refqrs(idx) - ssamp;  % within interval
    % = apply ICA
    % not re-using matrices since stationarity is not a possibility
    [dataICA,~,~] = fastica(data(:,samp2filt),'approach','symm','g','tanh','verbose','off');
    dataICA = diag(1./max(dataICA'))*dataICA; % may not have the same size as "data"
    % = QRS detect each component and take F1, RMS measure
    qrsdet = cell(1,size(dataICA,1));
    F1 = zeros(1,size(dataICA,1));
    RMS = zeros(1,size(dataICA,1));
    for ch = 1:size(dataICA,1)
        qrsdet{ch} = qrs_detect(dataICA(ch,:),TH,REFRAC,FS);
        if ~isempty(qrsdet{ch})
            [F1(ch),RMS(ch)] = Bxb_compare(refint,qrsdet{ch},INTERV);
        else
            F1(ch) = 0;
            RMS(ch) = NaN;
        end
    end
    
    % = decide for one component
    %Saving segment
    [maxF1,maxch] = max(F1);
    if maxF1 > .8
        icamorph(samp2filt) = dataICA(maxch,:);
        qrsica(qrsdet{maxch}+ssamp-1) = 1;
    end
    
    %Augment offsets
    ssamp = endsamp+1;
    endsamp = endsamp+blength;
end
qrsica = find(qrsica);

%% Calculate F1 measure
[F1,RMS,PPV,SE] = Bxb_compare(refqrs,qrsica,INTERV);
