function [fqrs,maxch] = FECGSYN_QRSmincompare(data,fref,fs,varargin)
% function [fqrs,maxch] = FECGSYN_QRSmincompare(data,fref,fs,varargin)
% BxB compare on minute basis
%
% This is function detects and performs comparison between reference and
% detection on a minute basis or pre-defined interval.
%
% data                         Extracted signals
% fref                         Fetal referenec signal (samplestamps)
% fs                           Sampling frequency (in Hz)
% window (optional input)      Window to compare beats (in seconds) 
% 
% 
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 17-12-2015
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

% Input check
switch length(varargin)
    case 0
        window = 60; % default minute window
    case 1
        window = varargin{1};
    otherwise
        error('FECGSYN_QRSmincompare: Too many inputs')
end



% Parameters
INTERV = round(0.05*fs);    % BxB acceptance interval
TH = 0.3;                   % detector threshold
REFRAC = .15;               % detector refractory period (in s)

% Detect QRS complexes
fqrs = cell(1,size(data,1));
for j = 1:size(data,1)
    fqrs{j} = qrs_detect(data(j,:),TH,REFRAC,fs);
end

% creating statistics in 1-min blocks
min = 1;
numblock = length(data)/fs/window;
if rem(numblock,1) ~= 0
    warning('FECGSYN_QRSmincompare: non-integer division of data length and block window. Dataset will be only partially evaluated!')
    numblock = floor(numblock);
end
maxch = zeros(1,numblock);
fqrs_temp = cell(1,numblock);
while min <= length(data)/fs/window;
    F1max = 0;
    idxref = (fref>=(min-1)*fs*window+1)&(fref<=min*fs*window);
    for j = 1:size(data,1)
        idx = (fqrs{j}>=(min-1)*fs*window+1)&(fqrs{j}<=min*fs*window);
        [F1,~,~,~] = Bxb_compare(fref(idxref),fqrs{j}(idx),INTERV);
        if F1 > F1max    % compare and see if this channel provides max F1
            maxch(min) = j;
            F1max = F1;
            fqrs_temp{min} = fqrs{j}(idx);%+ (min-1)*fs*60;    % adding fqrs detections to temporary cell
        end
    end
    min = min+1;
end
fqrs = cell2mat(fqrs_temp);
end