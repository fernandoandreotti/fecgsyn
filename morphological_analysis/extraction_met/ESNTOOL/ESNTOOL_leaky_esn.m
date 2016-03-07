function internalState = ESNTOOL_leaky_esn(totalstate , esn , varargin )
% Update internal state using the leaky integrator neuron model.
%
% input arguments:
% totalstate: the previous totalstate, vector of size 
%     (esn.nInternalUnits + esn.nInputUnits + esn.nOutputUnits) x 1
% esn: the ESN structure
%
% output: 
% internalState: the updated internal state, size esn.nInternalUnits x 1
%
% Created June 7, 2006, H.Jaeger
% Revision 1, June 23, 2007, H.Jaeger  (include esn.feedbackScaling)
% Revision 2, July 1, 2007, H.Jaeger (change from uniform timeConstant to
%                                     neuron-specific timeConstants

    previousInternalState = totalstate(1:esn.nInternalUnits, 1);
        %% WAS - changed Joachim Behar 07/01/2013
        internalState = (1 -  esn.leakage * esn.timeConstants) .* previousInternalState + esn.timeConstants .* ...
        feval(esn.reservoirActivationFunction ,...
        [ esn.internalWeights, esn.inputWeights, esn.feedbackWeights * diag(esn.feedbackScaling )] * totalstate) ; 
        %% now is.. according to Petrenas paper although my alpha is
        % petrenas 1-alpha
%         internalState = (1 -  esn.leakage * esn.timeConstants) .* previousInternalState + esn.leakage*esn.timeConstants .* ...
%         feval(esn.reservoirActivationFunction ,...
%         [ esn.internalWeights, esn.inputWeights, esn.feedbackWeights * diag(esn.feedbackScaling )] * totalstate) ; 
%     
    %%%% add noise to the Esn 
internalState = internalState + esn.noiseLevel * (rand(esn.nInternalUnits,1) - 0.5) ; 
