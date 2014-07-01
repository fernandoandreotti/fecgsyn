function [F1,RMS,ACC,PPV,SE,TP,FN,FP] = Bxb_compare(refqrs, testqrs, acceptint)
% This function is similar to the function bxb.exe from Physionet. It
% compares in a beat-by-beat basis if the detections match the reference.
% The algorithm is based on the entry by Joachim Behar on the Physionet / 
% Computing in Cardiology Challenge 2013 and on ANSI/AAMI EC57 Norm 1998
% 
% Input
% refqrs:        reference QRS detections
% testqrs:       detections to be tested against 
% acceptint:     acceptance interval (left and right) in samples
% 
% 
% Output
% F1:            F1-measure (Joachim Behar - Computing in Cardiology 2013)
% ACC:           accuracy (by Karvounis 2007)
% PPV:           positive predictive value
% SE:            sensitivity
% TP:            number of true positives
% FN:            number of false negatives
% FP:            number of false positives
%
% References
%
% TODO
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

% == input test

if size(refqrs,2) > size(refqrs,1)
    refqrs = refqrs';
end

if size(testqrs,2) > size(testqrs,1)
    testqrs = testqrs';
end

NB_REF = length(refqrs);
NB_TEST = length(testqrs);

% == core function
 [idxmatch,dist] = dsearchn(refqrs,testqrs);         % closest ref for each point in test qrs
 idxmatchwin = idxmatch(dist<acceptint);         % keep only the ones within a certain window
 RMS = mean(dist(dist<acceptint));        % RMS value
 NB_MATCH_UNIQUE = length(unique(idxmatchwin)); % how many unique matching
 
 
% == generating output
 TP = NB_MATCH_UNIQUE;             % number of identified ref QRS
 FN = NB_REF-TP;                   % number of missed ref QRS
 FP = NB_TEST-TP;                  % how many extra detection?
 SE  = TP/(TP+FN);
 PPV = TP/(FP+TP);
 ACC=TP/(TP+FP+FN);         
 F1 = 2*SE*PPV/(SE+PPV);           % accuracy measure
