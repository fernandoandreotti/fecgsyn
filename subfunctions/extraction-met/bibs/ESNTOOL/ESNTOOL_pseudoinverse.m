function outputWeights = ESNTOOL_pseudoinverse(stateCollectMat, teachCollectMat)

% PSEUDOINVERSE computes the outputWeights using the standard pseudoinverse
% operation
%
% inputs:
% stateCollectMat = matrix of size (nTrainingPoints-nForgetPoints) x
% nInputUnits + nInternalUnits 
% stateCollectMat(i,j) = internal activation of unit j after the 
% (i + nForgetPoints)th training point has been presented to the network
% teacherCollectMat is a nSamplePoints * nOuputUnits matrix that keeps
% the expected output of the ESN
% teacherCollectMat is the transformed(scaled, shifted etc) output see
% compute_teacher for more documentation
%
% output:
% outputWeights = vector of size (nInputUnits + nInternalUnits) x 1
% containing the learnt weights
%
% 
% 
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% 
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

  
outputWeights = (pinv(stateCollectMat)*teachCollectMat)' ; 
