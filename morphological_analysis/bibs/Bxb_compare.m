function [TP,FN,FP,SE,PPV,ACC,F1] = Bxb_compare(refqrs, testqrs, acceptint)
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
% TP:            number of true positives
% FN:            number of false negatives
% FP:            number of false positives
% SE:            sensitivity
% PPV:           positive predictive value
% ACC:           accuracy (by Karvounis 2007)
% F1:            F1-measure (Joachim Behar - Computing in Cardiology 2013)

NB_REF = length(refqrs);
NB_TEST = length(testqrs);

% == core function
 [idxmatch,dist] = dsearchn(refqrs,testqrs);         % closest ref for each point in test qrs
 idxmatchwin = idxmatch(dist<acceptint);         % keep only the ones within a certain window
 NB_MATCH_UNIQUE = length(unique(idxmatchwin)); % how many unique matching
 
% == 
 TP = NB_MATCH_UNIQUE;             % number of identified ref QRS
 FN = NB_REF-TP;                   % number of missed ref QRS
 FP = NB_TEST-TP;                  % how many extra detection?
 SE  = TP/(TP+FN);
 PPV = TP/(FP+TP);
 ACC=TP/(TP+FP+FN);         
 F1 = 2*SE*PPV/(SE+PPV);           % accuracy measure
