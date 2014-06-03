function [F1,Se,PPV,Nb] = stats(refqrs,testqrs,thres,margin,windowlen,fs)
% sqi = bsqi(refqrs,testqrs,thres,margin,windowlen,fs)
% compare two sets of annotation with one as the reference (refqrs) and one
% as the test (testqrs)
%
% inputs
%     refqrs:       reference qrs annotation (in sec)
%     testqrs:      test qrs annotations (in sec)
%     thres:        threshold (in sec,default 0.05s)
%     margin:       margin time not include in comparison (in sec,default 2s)
%     windowlen:    length of the comparison window (in sec,default 60s)
%     fs:           sampling frequency
%
% output
%     sqi: match proportion according to some criteria you can change
%     depending on what you are looking for (can be Se, PPV or F1 measure).
%     See at the end of the function.
%
% Reference: inspired from
%     Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation 
%     from multiple asynchronous noisy sources using signal quality indices and 
%     a Kalman filter." Physiological measurement 29.1 (2008): 15.
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

% == managing inputs
if nargin<2; error('bsqi: wrong number of input arguments \n'); end;
if nargin<3; thres=0.05; end;
if nargin<4; margin=2; end;
if nargin<5; windowlen=60; end;
if nargin<6; fs=1000; end;

if size(refqrs,1)>size(refqrs,2); refqrs=refqrs';end
if size(testqrs,1)>size(testqrs,2); testqrs=testqrs';end

start = margin*fs;
stop = (windowlen-margin)*fs;
refqrs = refqrs*fs;
testqrs = testqrs*fs;

try
    refqrs  = refqrs(refqrs>start & refqrs<stop)'; % reference annotations
    testqrs = testqrs(testqrs>start & testqrs<stop)'; % test annotations
    
    NB_REF  = length(refqrs);
    NB_TEST = length(testqrs);
    
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

    Nb.TP = TP;
    Nb.FN = FN;
    Nb.FP = FP;
    
    fprintf('Se: %f \n',Se);
    fprintf('PPV: %f \n',PPV);
    fprintf('F1: %f \n',F1);    
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    rethrow(ME);
end

end































