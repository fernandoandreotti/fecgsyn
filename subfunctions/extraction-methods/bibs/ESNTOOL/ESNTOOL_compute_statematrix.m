function stateCollectMat = ...
    ESNTOOL_compute_statematrix(inputSequence, outputSequence, esn, nForgetPoints, varargin)
% compute_statematrix  runs the input through the ESN and writes the
% obtained input+reservoir states into stateCollectMat.
% The first nForgetPoints will be deleted, as the first few states could be
% not reliable due to initial transients  
%
% inputs:
% inputSequence = input time series of size nTrainingPoints x nInputDimension
% outputSequence = output time series of size nTrainingPoints x nOutputDimension
% esn = an ESN structure, through which we run our input sequence
% nForgetPoints: an integer, may be negative, positive or zero.
%    If positive: the first nForgetPoints will be disregarded (washing out
%    initial reservoir transient)
%    If negative: the network will be initially driven from zero state with
%    the first input repeated |nForgetPoints| times; size(inputSequence,1)
%    many states will be sorted into state matrix
%    If zero: no washout accounted for, all states except the zero starting
%    state will be sorted into state matrix
%
% Note: one of inputSequence and outputSequence may be the empty list [],
% but not both. If the inputSequence is empty, we are dealing with a purely
% generative task; states are then computed by teacher-forcing
% outputSequence. If outputSequence is empty, we are using this function to
% test a trained ESN; network output is then computed from network dynamics
% via output weights. If both are non-empty, states are computed by
% teacher-forcing outputSequence.
%
% optional input argument:
% there may be one optional input, the starting vector by which the esn is
% started. The starting vector must be given as a column vector of
% dimension esn.nInternalUnits + esn.nOutputUnits + esn.nInputUnits  (that
% is, it is a total state, not an internal reservoir state). If this input
% is desired, call test_esn with fourth input 'startingState' and fifth
% input the starting vector.
%
% output:
% stateCollectMat = matrix of size (nTrainingPoints-nForgetPoints) x
% nInputUnits + nInternalUnits 
% stateCollectMat(i,j) = internal activation of unit j after the 
% (i + nForgetPoints)th training point has been presented to the network
%
%
% Version 1.0, April 30, 2006
% Copyright: Fraunhofer IAIS 2006 / Patents pending
% Revision 1, June 6, 2006, H. Jaeger
% Revision 2, June 23, 2007, H. Jaeger (added optional starting state
% input)
% Revision 3, July 1, 2007, H. Jaeger (added leaky1_esn update option)
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


if isempty(inputSequence) && isempty(outputSequence)
    error('error in compute_statematrix: two empty input args');
end

if isempty(outputSequence)
    teacherForcing = 0;
    nDataPoints = length(inputSequence(:,1));
else
    teacherForcing = 1;
    nDataPoints = length(outputSequence(:,1));
end

if nForgetPoints >= 0
    stateCollectMat = ...
        zeros(nDataPoints - nForgetPoints, esn.nInputUnits + esn.nInternalUnits) ; 
else
    stateCollectMat = ...
        zeros(nDataPoints, esn.nInputUnits + esn.nInternalUnits) ; 
end

%% set starting state

externalStartStateFlag = 0;
args = varargin; 
nargs= length(args);
for i=1:2:nargs
    switch args{i},
        case 'startingState', 
            totalstate = args{i+1} ; 
            internalState = totalstate(1:esn.nInternalUnits,1) ; 
            externalStartStateFlag = 1;
        otherwise error('the option does not exist'); 
    end      
end

if externalStartStateFlag == 0
    totalstate = zeros(esn.nInputUnits + esn.nInternalUnits + esn.nOutputUnits, 1);
    internalState = zeros(esn.nInternalUnits, 1);
end

%%%% if nForgetPoints is negative, ramp up ESN by feeding first input
%%%% |nForgetPoints| many times

if nForgetPoints < 0
    for i = 1:-nForgetPoints
        if esn.nInputUnits > 0
            in = esn.inputScaling .* inputSequence(1,:)' + esn.inputShift;  % in is column vector
        else in = [];
        end
        if esn.nInputUnits > 0
            totalstate(esn.nInternalUnits+1:esn.nInternalUnits + esn.nInputUnits) = in;
        end 
        % the internal state is computed based on the type of the network
        switch esn.type
            case 'ESNTOOL_leaky_esn'
                typeSpecificArg = [];
            otherwise
                error('Not implemented in FECGSYN')
        end
        internalState = feval(esn.type, totalstate, esn, typeSpecificArg) ; 
        
        if teacherForcing
            netOut = esn.teacherScaling .* outputSequence(1,:)' + esn.teacherShift;
        else
            netOut = feval(esn.outputActivationFunction, esn.outputWeights * [internalState; in]);
        end
        
        totalstate = [internalState; in; netOut];
    end
end

collectIndex = 0;
for i = 1:nDataPoints
    
    % scale and shift the value of the inputSequence
    if esn.nInputUnits > 0
        in = esn.inputScaling .* inputSequence(i,:)' + esn.inputShift;  % in is column vector
    else in = [];
    end
    
    % write input into totalstate
    if esn.nInputUnits > 0
        totalstate(esn.nInternalUnits+1:esn.nInternalUnits + esn.nInputUnits) = in;
    end    
    
    % the internal state is computed based on the type of the network
    switch esn.type
        case 'ESNTOOL_leaky_esn'
            typeSpecificArg = [];
        otherwise
            error('Not implemented within FECGSYN!')      
    end
    internalState = feval(esn.type, totalstate, esn, typeSpecificArg) ; 
    
    if teacherForcing
        netOut = esn.teacherScaling .* outputSequence(i,:)' + esn.teacherShift;
    else        
        netOut = feval(esn.outputActivationFunction, esn.outputWeights * [internalState; in]);
    end
    
    
    totalstate = [internalState; in; netOut];
    
    %collect state
    if nForgetPoints >= 0 &&  i > nForgetPoints
        collectIndex = collectIndex + 1;
        stateCollectMat(collectIndex,:) = [internalState' in']; 
    elseif nForgetPoints < 0
        collectIndex = collectIndex + 1;
        stateCollectMat(collectIndex,:) = [internalState' in']; 
    end
    
end
