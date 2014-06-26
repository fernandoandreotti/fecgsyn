function [TP,FN,FP,Se,pP,Acc,F1,RMS] = Bxb_compare(refqrs, testqrs, acceptint)
% This function calculates the 

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
