function esn = ESNTOOL_generate_esn(nInputUnits, nInternalUnits, nOutputUnits, varargin)
% Creates an ESN set up for use in multiple-channel output association tasks. 
% The number of input, internal, and output 
% units have to be set. Any other option is set using the format 
% 'name_of_options1',value1,'name_of_option2',value2, etc.
% 
%%%%% input arguments:
% nInputUnits: the dimension of the input 
% nInternalUnits: size of the Esn
% nOutputUnits: the dimension of the output
%
%%%%% optional arguments:
% 'inputScaling': a nInputUnits x 1 vector
%
% 'inputShift': a nInputUnits x 1 vector. 
%
% 'teacherScaling': a nOutputUnits x 1 vector
%
% 'teacherShift': a nOutputUnits x 1 vector. 
%
% 'noiseLevel': a small number containing the amount of uniform noise to be
%  added when computing the internal states
%
% 'learningMode': a string ('offline_singleTimeSeries', 'offline_multipleTimeSeries' or 'online')
%     1. Case 'offline_singleTimeSeries': trainInput and trainOutput each represent a 
%        single time series in an array of size sequenceLength x sequenceDimension
%     2. Case 'offline_multipleTimeSeries': trainInput and trainOutput each represent a 
%        collection of K time series, given in cell arrays of size K x 1, where each cell is an
%        array of size individualSequenceLength x sequenceDimension
%     3. Case 'online': trainInput and trainOutput are a single time
%        series, output weights are adapted online
%
% 'reservoirActivationFunction': a string ("tanh", "identity", "sigmoid01") ,
%
% 'outputActivationFunction': a string("tanh", "identity", "sigmoid01") ,
%
% 'inverseOutputActivationFunction': the inverse to
%    outputActivationFunction, one of 'atanh', 'identity', 'sigmoid01_inv'.
%    When choosing the activation function, make sure the inverse
%    activation function is corectly set.
%
% 'methodWeightCompute': a string ('pseudoinverse', 'wiener_hopf'). It  
%    specifies which method to use to compute the output weights given the
%    state collection matrix and the teacher
%
% 'spectralRadius': a positive number less than 1. 
%
% 'feedbackScaling': a nOutputUnits x 1 vector, indicating the scaling
%     factor to be applied on the output before it is fed back into the network
%
% 'type': a string ('plain_esn', 'leaky_esn' or 'twi_esn')
% 'trained': a flag indicating whether the network has been trained already
% 'timeConstants': option used in networks with type == "leaky_esn", "leaky1_esn" and "twi_esn".
%                      Is given as column vector of size esn.nInternalUnitsm, where each entry 
%                      signifies a time constant for a reservoir neuron.
% 'leakage': option used in networks with type == "leaky_esn" or "twi_esn"
% 'RLS_lambda': option used in online training(learningMode == "online") 
% 'RLS_delta': option used in online training(learningMode == "online")
%
% for more information on the Echo State network approach take a look at
% the following tutorial : 
% http://www.faculty.iu-bremen.de/hjaeger/pubs/ESNTutorialRev.pdf

% Version 1.0, April 30, 2006
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 6, 2006, H. Jaeger
% Revision 2, Feb 23, 2007, H. Jaeger
% Revision 3, June 27, 2007, H. Jaeger 
% Revision 4, July 1, 2007, H. Jaeger: 
%    - changed esn.timeConstant to esn.timeConstants
%    - deleted esn.retainment
%    - deleted esn.internalWeights (to enforce that this is set outside
%    this generation script)
% Revision 5, July 29, 2007, H. Jaeger: bugfix (for cases of zero input or
%    output length)
% Revision 6, Jan 28, 2009, H. Jaeger: bugfix (deleted defunct error
%                             catching routine for leaky / twiesn
%                             reservoirs)
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


%%%% set the number of units
esn.nInternalUnits = nInternalUnits; 
esn.nInputUnits = nInputUnits; 
esn.nOutputUnits = nOutputUnits; 
  
connectivity = min([20/nInternalUnits 1]);
nTotalUnits = nInternalUnits + nInputUnits + nOutputUnits; 

esn.internalWeights_UnitSR = ESNTOOL_generate_internal_weights(nInternalUnits, ...
                                                connectivity);
esn.nTotalUnits = nTotalUnits; 

% input weight matrix has weight vectors per input unit in colums
esn.inputWeights = 2.0 * rand(nInternalUnits, nInputUnits)- 1.0;

% output weight matrix has weights for output units in rows
% includes weights for input-to-output connections
esn.outputWeights = zeros(nOutputUnits, nInternalUnits + nInputUnits);

%output feedback weight matrix has weights in columns
esn.feedbackWeights = (2.0 * rand(nInternalUnits, nOutputUnits)- 1.0);

%init default parameters
if nInputUnits > 0
    esn.inputScaling  = ones(nInputUnits, 1);  esn.inputShift    = zeros(nInputUnits, 1);
else
    esn.inputScaling = []; esn.inputShift = [];
end
if nOutputUnits > 0
    esn.teacherScaling = ones(nOutputUnits, 1); esn.teacherShift  = zeros(nOutputUnits, 1);
else
    esn.teacherScaling = []; esn.teacherShift  =  [];
end
esn.noiseLevel = 0.0 ; 
esn.reservoirActivationFunction = 'tanh';
esn.outputActivationFunction = 'ESNTOOL_identity';       % options: identity or tanh or sigmoid01
esn.methodWeightCompute = 'ESNTOOL_pseudoinverse' ; % options: pseudoinverse and wiener_hopf and ridge_reg
esn.inverseOutputActivationFunction = 'ESNTOOL_identity' ;  % options:
                                                    % identity or
                                                    % atanh or sigmoid01_inv
esn.spectralRadius = 1 ; 
esn.feedbackScaling = zeros(nOutputUnits, 1); 
esn.trained = 0 ; 
esn.type = 'plain_esn' ; 
esn.timeConstants = ones(esn.nInternalUnits,1); 
esn.leakage = 0.5;  
esn.learningMode = 'offline_singleTimeSeries' ; 
esn.RLS_lambda = 1 ; 

args = varargin; 
nargs= length(args);
for i=1:2:nargs
  switch args{i},
   case 'inputScaling', esn.inputScaling = args{i+1} ; 
   case 'inputShift', esn.inputShift= args{i+1} ; 
   case 'teacherScaling', esn.teacherScaling = args{i+1} ; 
   case 'teacherShift', esn.teacherShift = args{i+1} ;     
   case 'noiseLevel', esn.noiseLevel = args{i+1} ; 
   case 'learningMode', esn.learningMode = args{i+1} ; 
   case 'reservoirActivationFunction',esn.reservoirActivationFunction=args{i+1};
   case 'outputActivationFunction',esn.outputActivationFunction=  ...
                        args{i+1};        
   case 'inverseOutputActivationFunction', esn.inverseOutputActivationFunction=args{i+1}; 
   case 'methodWeightCompute', esn.methodWeightCompute = args{i+1} ; 
   case 'spectralRadius', esn.spectralRadius = args{i+1} ;  
   case 'feedbackScaling',  esn.feedbackScaling = args{i+1} ; 
   case 'type' , esn.type = args{i+1} ; 
   case 'timeConstants' , esn.timeConstants = args{i+1} ; 
   case 'leakage' , esn.leakage = args{i+1} ; 
   case 'RLS_lambda' , esn.RLS_lambda = args{i+1};
   case 'RLS_delta' , esn.RLS_delta = args{i+1};
       
   otherwise error('the option does not exist'); 
  end      
end

%%%% error checking
% check that inputScaling has correct format
if size(esn.inputScaling,1) ~= esn.nInputUnits
    error('the size of the inputScaling does not match the number of input units'); 
end
if size(esn.inputShift,1) ~= esn.nInputUnits
    error('the size of the inputScaling does not match the number of input units'); 
end
if size(esn.teacherScaling,1) ~= esn.nOutputUnits
    error('the size of the teacherScaling does not match the number of output units'); 
end
if size(esn.teacherShift,1) ~= esn.nOutputUnits
    error('the size of the teacherShift does not match the number of output units'); 
end
if length(esn.timeConstants) ~= esn.nInternalUnits
    error('timeConstants must be given as column vector of length esn.nInternalUnits'); 
end
if ~strcmp(esn.learningMode,'offline_singleTimeSeries') &&...
        ~strcmp(esn.learningMode,'offline_multipleTimeSeries') && ...
        ~strcmp(esn.learningMode,'online')
    error('learningMode should be either "offline_singleTimeSeries", "offline_multipleTimeSeries" or "online" ') ; 
end
if ~((strcmp(esn.outputActivationFunction,'ESNTOOL_identity') && ...
        strcmp(esn.inverseOutputActivationFunction,'ESNTOOL_identity')) || ...
        (strcmp(esn.outputActivationFunction,'tanh') && ...
        strcmp(esn.inverseOutputActivationFunction,'atanh')) || ...
        (strcmp(esn.outputActivationFunction,'sigmoid01') && ...
        strcmp(esn.inverseOutputActivationFunction,'sigmoid01_inv')))  ...
error('outputActivationFunction and inverseOutputActivationFunction do not match'); 
end


