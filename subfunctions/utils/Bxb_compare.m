function [F1,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(refqrs, testqrs, acceptint)
% This function is similar to the function bxb.exe from Physionet's
% WFDB-App toolbox. It compares in a beat-by-beat basis if the detections 
% match the reference. The algorithm is based on the entry by Joachim Behar 
% on the Physionet / Computing in Cardiology Challenge 2013 and on ANSI/AAMI 
% EC57 Norm 1998
% 
% Input
% refqrs:        reference QRS detections
% testqrs:       detections to be tested against 
% acceptint:     acceptance interval (left and right) in samples
% 
% 
% Output
% F1:            F1-measure (Joachim Behar - Computing in Cardiology 2013)
% ACC:           accuracy (by Karvounis 2007) - alternative to F1
% PPV:           positive predictive value
% SE:            sensitivity
% TP:            number of true positives
% FN:            number of false negatives
% FP:            number of false positives
%
% References
% [ANSI/AAMI Norm]  American National Standard ANSI/AAMI EC57:1998, Testing and Reporting Performance 
% Results of Cardiac Rhythm and ST Segment Measurement Algorithms 
% 
% [WFDB-APP] Silva, I, Moody, G. "An Open-source Toolbox for Analysing and Processing PhysioNet Databases 
% in MATLAB and Octave." Journal of Open Research Software 2(1):e27 [http://dx.doi.org/10.5334/jors.bi]; 
% 2014 (September 24). 
% 
% [Behar2014] Behar, J., Oster, J., & Clifford, G. D. (2014). Combining and Benchmarking Methods of Foetal
% ECG Extraction Without Maternal or Scalp Electrode Data. Physiological Measurement, 35(8), 1569–1589.
% 
% 
% --
% fecgsyn toolbox, version 1.2, Jan 2017
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% University of Oxford, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
% 
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. 
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. 
%   Physiological Measurement.35 1537-1550. 2014.
%
% Last updated : 10-03-2016
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
try
    [idxmatch,dist] = dsearchn(refqrs,testqrs);     % closest ref for each point in test qrs
    idxmatchwin = idxmatch(dist<acceptint);         % keep only the ones within a certain window
    %  RMS = sqrt(mean(dist(dist<acceptint).^2));      % RMS value
    MAE = sum(abs(dist(dist<acceptint)));      % RMS value
catch % case testqrs is empty
    idxmatchwin = [];
    MAE = NaN;
end
 
 NB_MATCH_UNIQUE = length(unique(idxmatchwin)); % how many unique matching
 
 
% == generating output
 TP = NB_MATCH_UNIQUE;             % number of identified ref QRS
 MAE = MAE/TP;                     % mean absolute distance
 FN = NB_REF-TP;                   % number of missed ref QRS
 FP = NB_TEST-TP;                  % how many extra detection?
 SE  = TP/(TP+FN);
 
 if isnan(MAE)                     % define PPV for zero detections
     PPV = 0;
     F1 = 0;                       % F1 definition
 else
     PPV = TP/(FP+TP);
     F1 = 2*SE*PPV/(SE+PPV);           % accuracy measure
 end
%  ACC=TP/(TP+FP+FN);         

