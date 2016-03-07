function internalWeights = ESNTOOL_generate_internal_weights(nInternalUnits, ...
                                                  connectivity)
% GENERATE_INTERNAL_WEIGHTS creates a random reservoir for an ESN
%  
% inputs:
% nInternalUnits = the number of internal units in the ESN
% connectivity \in [0,1], says how many weights should be non-zero 
%
% output:
% internalWeights = matrix of size nInternalUnits x nInternalUnits
% internalWeights(i,j) = value of weight(synapse) from unit i to unit j
% internalWeights(i,j) might be different from internalWeights(j,i)
% 
%
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, Feb 23, 2007, H. Jaeger
% Revision 2, March 10, 2007, H. Jaeger (replaced eigs by myeigs)
% Reference
% Jaeger, H. (2001). The “echo state” approach to analysing and training recurrent neural networks.
% 
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 30, 2006, H. Jaeger
% Revision 2, Feb 23, 2007, H. Jaeger
% 
% ============================================================
% This software is free for non-commercial use. If commercial use is
% intended, contact Fraunhofer IAIS (www.iais.fraunhofer.de) who have
% claimed international patents on ESN algorithms (pending).
% 
% This software is intended for research use by experienced Matlab users
% and includes no warranties or services.  Bug reports are gratefully
% received by Herbert Jaeger (h.jaeger [at] jacobs-university.de)
% 
% ============================================================

success = 0 ;                                               
while success == 0
    % following block might fail, thus we repeat until we obtain a valid
    % internalWeights matrix
    try
        internalWeights = sprand(nInternalUnits, nInternalUnits, connectivity);
        internalWeights(internalWeights ~= 0) = internalWeights(internalWeights ~= 0)  - 0.5;
        %maxVal = max(abs(myeigs(internalWeights,1)));
        maxVal = max(abs(eigs(internalWeights,1)));
        internalWeights = internalWeights/maxVal;
        success = 1 ;
    catch
    %    success = 0 ; 
    end
end

% % Prof Peter Tino type of connectivity
% internalWeights = zeros(nInternalUnits,nInternalUnits);
% internalWeights(1:nInternalUnits,1)     = 1;
% internalWeights(1,end)                  = 1;
% v = 0.5;
% internalWeights = v*internalWeights;



