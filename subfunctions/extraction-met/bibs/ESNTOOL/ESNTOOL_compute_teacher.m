function teachCollectMat = ESNTOOL_compute_teacher(outputSequence, esn, ...
    nForgetPoints)
% COMPUTE_TEACHER scales, shifts and applies the inverse output
% activation function on the exepcted teacher. 
% the first nForgetPoints are being disregarded
%
% inputs:
% outputSequence = teacher vector of size nTrainingPoints x nOutputDimension
% esn = an ESN structure, which contains the information about the
% transformation we need to apply to the teacher 
% nForgetPoints: an integer, may be negative, positive or zero.
%    If positive: the first nForgetPoints will be disregarded (washing out
%    initial reservoir transient)
%    If negative: the network will be initially driven from zero state with
%    the first input repeated |nForgetPoints| times; size(inputSequence,1)
%    many states will be sorted into state matrix
%    If zero: no washout accounted for, all states except the zero starting
%    state will be sorted into state matrix
%
% outputs:
% teachCollectMat = matrix of size (nOutputPoints - nForgetPoints) x
% nOutputUnits
% teachCollectMat contains the shifted and scaled output
% 
%
% Version 1.0, April 30, 2006
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 7, 2006, H.Jaeger
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


nOutputPoints  = length(outputSequence(:,1)) ; 
teachCollectMat = zeros(nOutputPoints - max([0, nForgetPoints]), esn.nOutputUnits) ;

% delete the first nForgetPoints elements from outputSequence
if nForgetPoints >= 0
    outputSequence = outputSequence(nForgetPoints+1:end,:) ; 
end

% update the size of outputSequence
nOutputPoints  = length(outputSequence(:,1)) ; 

teachCollectMat = [(diag(esn.teacherScaling) * outputSequence')' + ...
        repmat(esn.teacherShift',[nOutputPoints 1])];

teachCollectMat = feval(esn.inverseOutputActivationFunction, teachCollectMat);




