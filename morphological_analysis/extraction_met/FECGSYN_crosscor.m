function [match,coeff] = FECGSYN_crosscor(cycleA,cycleB,thres)
% cross correlation between two cycles to identify wether or not they match
% with each other. cycleA and cycleB are the two cycles to compare
% thres is a decimal number in [0 1] what degree of similarity is requested
% to consider that the two cycles are matching.
%
% inputs
%   cycleA & cycleB: cycles to compare with each other
%   thres:           threshold at which to decide whether the cycles match or not
%                    (range 0-1)
%
% outputs
%   match:           are the two cycles matching? [bool]
%   coeff:           correlation coefficient
%
% Copyright (C) 2014  Julien Oster, Joachim Behar, Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 18-08-2014, Joachim

% == format the intput
[la,ca] = size(cycleA);
[lb,cb] = size(cycleB);
if la>ca; cycleA = cycleA'; end;
if lb>cb; cycleB = cycleB'; end;

% == init variables
match = 0;
nbLeads = size(cycleA,1);
NB_BINS = size(cycleA,2);
coeff = zeros(1,nbLeads);

% == compute correlation coeff
for i=1:nbLeads
    C = corrcoef(cycleA,cycleB);    
    coeff(i) = C(1,2);
end

if mean(coeff)>thres; 
    match = 1; 
end

end