function bestsig = ica_extraction(data,chan,refqrs,varargin)
% Uses Independent Component Analysis (ICA) over selected channels input 
% data and choses best channel based on th F1 measure.
% 
% Input
% data:      Matrix containing signals to serve as input to ICA.
% chan:      Channels from data to be used in ICA
% refqrs:    Array containing reference QRS detections for F1 measure
% (optional)
% blength    Iterates ICA every blength (in samples)
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

%% Run ICA on every block
ssamp = 1;          % starting sample to filter (offset)
endsamp = ssamp + blength - 1;      % ending sample to filter
loop = 1;           % allows iterations

if size(data,2)<endsamp
    blength = Data.length;
    endsamp = size(data,2);
end

% = normalise data
trans = bsxfun(@minus,data,mean(data,2)); % remove mean
X = bsxfun(@rdivide,trans,std(data,0,2)); % divide by standard deviation
icasignal = zeros(size(X));

while (loop)  %quit will be given as soon as complete signal is filtered
    if (size(data,2) - ssamp) < 1.5*blength        % if there is less than 1.5x wrapping
        endsamp = size(data,2);              % interval, it should filter all
        loop = 0;                           % and be the last loop iteration
    end
    samp2filt = ssamp:endsamp;              % creating a list with samples to filter
    
    % = apply JADE
    % X = A.S + N
    tic
    %     [A,S] = jade(X(:,samp2filt));    % FIXME: NOT SURE WHY NEED TO MULTIPLY BY 1000
    [dataICA,A,W] = fastica(X(:,samp2filt),'approach','symm','g','tanh','verbose','off');
    toc
    
    %Saving this part
    %     icasignal(:,samp2filt)= S';
    icasignal(:,samp2filt)= dataICA;
    %Augment offsets
    ssamp = endsamp+1;
    endsamp = endsamp+blength;
end

%% Calculate F1 measure
bestsig = ref;

% BxB function by Behar (CinC2013)
 % == core function
 [IndMatch,Dist] = dsearchn(refqrs,testqrs);         % closest ref for each point in test qrs
 IndMatchInWindow = IndMatch(Dist<thres*fs);         % keep only the ones within a certain window
 NB_MATCH_UNIQUE = length(unique(IndMatchInWindow)); % how many unique matching
 TP = NB_MATCH_UNIQUE;                               % number of identified ref QRS
 FN = NB_REF-TP;                                     % number of missed ref QRS
 FP = NB_TEST-TP;                                    % how many extra detection?
 Se  = TP/(TP+FN);
 PPV = TP/(FP+TP);
 F1 = 2*Se*PPV/(Se+PPV);                             % accuracy measure
