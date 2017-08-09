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
% Revision 3, January 07, 2013 J. Behar (alpha and leakage adapted acording to Petrenas 2012)
% 
% References
% A. Petrenas, V. Marozas, L. Sornmo and A. Lukosevi¿ius, "An Echo State Neural Network for QRST 
% Cancellation During Atrial Fibrillation," in IEEE Transactions on Biomedical Engineering, vol. 59, 
% no. 10, pp. 2950-2957, Oct. 2012. doi: 10.1109/TBME.2012.2212895
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
% This function, originally from the ESN toolbox by H. Jaeger was modified:
% Updated, Jan 2014 - Joachim Behar (changes RLS implementation and other 
% minor syntax changes)
%
% Copyright (C) 2014  Joachim Behar
% University of Oxford, Intelligent Patient Monitoring Group
% joachim.behar@oxfordalumni.org
%
% Last updated : 24-02-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


previousInternalState = totalstate(1:esn.nInternalUnits, 1);
% Modified by Joachim Behar 07/01/2013
internalState = (1 -  esn.leakage * esn.timeConstants) .* previousInternalState + esn.timeConstants .* ...
    feval(esn.reservoirActivationFunction ,...
    [ esn.internalWeights, esn.inputWeights, esn.feedbackWeights * diag(esn.feedbackScaling )] * totalstate) ;
% Adapted alpha to match Petrenas
% petrenas 1-alpha
%         internalState = (1 -  esn.leakage * esn.timeConstants) .* previousInternalState + esn.leakage*esn.timeConstants .* ...
%         feval(esn.reservoirActivationFunction ,...
%         [ esn.internalWeights, esn.inputWeights, esn.feedbackWeights * diag(esn.feedbackScaling )] * totalstate) ;
%
% add noise to the Esn
internalState = internalState + esn.noiseLevel * (rand(esn.nInternalUnits,1) - 0.5) ;
