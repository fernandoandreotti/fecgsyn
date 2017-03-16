function [template,qrsloc,status] = FECGSYN_tgen(ecg,qrs,fs,debug)
% function [template,qrsloc,status] = FECGSYN_tgen(ecg,qrs,fs,debug)
% This function contructs a template ECG based on the location of the R-peaks, 
% based on (Oster 2015). A series of peaks that match with each other are stacked to
% build a template. This template can then be used for ecg morphological
% analysis or further processing. Note that the qrs location inputed must
% be as precise as possible. This approach for building the template ECG
% seems to be the best of the alternatives that were tested and leaves the
% freedom of having more than one mode (i.e. multiple ECG template can be
% built if there are different cycle morphology such as PVC)  but it is not
% particularly fast.
% The procedure for building the template is:
% 1. create average wrapped template
% 2. identify different modes present
%
% inputs
%   ecg:            the ecg channel(s)
%   qrs:            qrs location [number of samples]
%    fs:            sampling frequency
%
% outputs
% 
%   status:         bool, success or failed to extract a dominant mode
%
%
% Dependencies: FECGSYN_phase_wrap.m
%
%
% inputs
%   ecg:     ecg signal
%   qrs:     qrs fiducials
%   phase:   phase of the EGC
%   bins:    number of bins for the phase wrapping
%   fs:      sampling frequency
%   flag:    adjust 'iso' to zero
%   debug:   debug level (1/2)
%
% outputs
%   ecgme:   template ECG
%   ecgsd:   ECG standard deviation used as a proxy of signal quality
%            (i.e. trust of the Kalman in observations)
%   phaseme: mean phase of the template ECG
%   nbcyc:   number of cycles used to build the template
%
%
% Reference
% Oster, J., Behar, J., Sayadi, O., Nemati, S., Johnson, A., & Clifford, G. (2015). Semi-supervised ECG 
% Ventricular Beat Classification with Novelty Detection Based on Switching Kalman Filters. IEEE Trans. 
% Biomed. Eng., 62(9), 2125–2134. http://doi.org/10.1109/TBME.2015.2402236
% 
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% TODO
%
% See also:
% FECGSYN_kf_extraction
% FECGSYN_kf_linearization
% FECGSYN_kf_EKFilter
% 
% 
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% TODO
%
% See also:
% FECGSYN_kf_modelling
% FECGSYN_kf_extraction
% FECGSYN_kf_linearization
% FECGSYN_kf_ECGmodelling
% 
% --
% fecgsyn toolbox, version 1.2, March 2017
% Released under the GNU General Public License
%
% Copyright (C) 2017  Joachim Behar & Fernando Andreotti
% Department of Engineering Science, University of Oxford
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Behar, J., Andreotti, F., Zaunseder, S., Li, Q., Oster, J., & Clifford, G. D. (2014). An ECG Model for Simulating 
% Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings. Physiol. Meas., 35(8), 1537–1550.
% 
%
% Last updated : 15-03-2017
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

% == manage inputs
if nargin<2; error('ecg_template_build: wrong number of input arguments \n'); end;


% == constants
NB_BINS = 250; % number of beans onto which to wrap a cycle
NB_LEADS = size(ecg,1);
NB_SAMPLES = size(ecg,2);
NB_REL = 10; % number relevance. How many cycles minimum to consider that a mode is relevant? - UPDATE ME DEPENDING ON APPLICATION
MIN_NB_CYC = 30; % mininum number of cycles (will decrease THRES until this number of cycles is achieved) - UPDATE ME DEPENDING ON APPLICATION
THRES = 0.9; % threshold at which to decide whether the cycle match or not - UPDATE ME DEPENDING ON APPLICATION
MIN_TLEN = 0.35;     % minimum template length in ms - UPDATE ME DEPENDING ON APPLICATION
PACE = 0.1;
MIN_THRES = 0.6;
cycle = zeros(NB_LEADS,NB_BINS);
startCycle = 2;
NbModes = 1; % initialisation variable
relevantModeInd = []; % - UPDATE ME DEPENDING ON APPLICATION
relevantMode.NbCycles = 0;

% == linear phase wrapping (shift in -pi/6)
phase = FECGSYN_kf_phasecalc(qrs,NB_SAMPLES);
PhaseChangePoints = find(phase(2:end)<0&phase(1:end-1)>0);
NB_CYCLES = length(PhaseChangePoints);

% ==== core function ====
% == creating the different modes
cycleIndex = (PhaseChangePoints(startCycle)+1:PhaseChangePoints(startCycle+1));
for j=1:NB_LEADS
    cycle(j,:) = interp1(phase(cycleIndex),ecg(j,cycleIndex),linspace(-pi,pi,NB_BINS),'spline');
end
Mode{NbModes}.cycles = zeros(NB_LEADS,NB_BINS,1);       % cycles included
Mode{NbModes}.cycles(:,:,1) = cycle;
Mode{NbModes}.cycleMean = cycle;                        % average cycle
Mode{NbModes}.cycleStd = zeros(3,NB_BINS);              % standard deviation from cycles
Mode{NbModes}.NbCycles = 1;                             % number of cycles present
Mode{NbModes}.cycleLen = [];                            % length of cycles (in samples)

while relevantMode.NbCycles<MIN_NB_CYC && THRES>MIN_THRES
    % the THRES is lowered until a mode with mode than MIN_NB_CYC cycles is
    % detected or the THREShold is too low (which means no mode can be identified)
    for i=startCycle+1:NB_CYCLES-2
        cycleIndex = (PhaseChangePoints(i)+1:PhaseChangePoints(i+1));
        for j=1:NB_LEADS
            cycle(j,:) = interp1(phase(cycleIndex),ecg(j,cycleIndex),linspace(-pi,pi,NB_BINS),'spline');
        end
        match = 0;
        indMode = 1;
        while (~match&&indMode<=NbModes)
%             [r,~] = corrcoef(cycle,Mode{indMode}.cycleMean);
%             coeff = abs(r(1,2));
            % 50% computation time 
            coeff = prcorr2(cycle,Mode{indMode}.cycleMean);
            match = coeff>THRES;
            indMode = indMode+1;
        end
        if ~match  % if the new cycle does not match with the average cycle of any mode
            % then create a new mode
            NbModes=NbModes+1;
            Mode{NbModes}.cycles = zeros(NB_LEADS,NB_BINS,1);
            Mode{NbModes}.cycles(:,:,1) = cycle;
            Mode{NbModes}.cycleMean = cycle;
            Mode{NbModes}.cycleStd = zeros(3,NB_BINS);
            Mode{NbModes}.NbCycles = 1;
            Mode{NbModes}.cycleLen = length(cycleIndex);
        else % it it correlates then integrate it to the corresponding mode
            Mode{indMode-1}.NbCycles = Mode{indMode-1}.NbCycles+1;
            temp = Mode{indMode-1}.cycles ;
            Mode{indMode-1}.cycles = zeros(NB_LEADS,NB_BINS,Mode{indMode-1}.NbCycles);
            Mode{indMode-1}.cycles(:,:,1:end-1)= temp;
            Mode{indMode-1}.cycles(:,:,end)= cycle;
            Mode{indMode-1}.cycleMean = mean(Mode{indMode-1}.cycles,3);
            Mode{indMode-1}.cycleStd = std(Mode{indMode-1}.cycles,0,3);
            Mode{indMode-1}.cycleLen = [Mode{indMode-1}.cycleLen length(cycleIndex)];
        end
        
    end
    
    % == detecting what mode is relevant
    %   relevantMode:   structure containing cycle, cycleMean and cycleStd
    %                   representing how many cycles have been selected to build the stack, the
    %                   mean ecg cycle that is built upon these selected cycles and the
    %                   standard deviation for each point of the template cycle as an indicator
    %                   of the precision of the estimation. *Only the dominant mode is outputted
    %                   for this application.*
    for i=1:length(Mode)
        % minimum amount of cycles and length
        if Mode{i}.NbCycles>NB_REL && mean(Mode{i}.cycleLen) >= MIN_TLEN*fs
            relevantModeInd = [relevantModeInd i];
        end
    end
    relevantMode = Mode(relevantModeInd);
    
    if isempty(relevantMode)
        % if we did not catch anything then output rubbish
        relevantMode.cycleMean = ones(NB_BINS,1);
        relevantMode.cycleStd = ones(NB_BINS,1);
        relevantMode.NbCycles  = 0;
        relevantMode.cycleLen = 0;
        status = 0;
        template.avg = NaN;
        template.stdev = NaN;
        qrsloc = NaN;
    else
        relevantModeStruc = [relevantMode{:}];
        [~,pos] = max([relevantModeStruc.NbCycles]); % look for dominant mode
        relevantMode = relevantMode{pos}; % only return the dominant mode for this application
        status = 1;
        
        % == Converting template from bins back to samples
        desl = round(NB_BINS/6);
        template.avg = circshift(relevantMode.cycleMean',-desl);
        template.stdev = circshift(relevantMode.cycleStd',-desl); 
        template.avg = resample(template.avg,round(mean(relevantMode.cycleLen)),NB_BINS);
        template.stdev = resample(template.avg,round(mean(relevantMode.cycleLen)),NB_BINS)';
        qrsloc = round(round((NB_BINS/2 - desl)*(length(template.avg)/NB_BINS)));
        [~,delay]=max(abs(template.avg));
        if abs(qrsloc-delay)<3, qrsloc = delay; end  % allow some alignment
    end
    THRES = THRES-PACE;
end

% if debug
%     fprintf('The number of cycles constituting the dominant mode was %f \n',relevantMode.NbCycles);
%     fprintf('The correlation threshold at which this happened was    %f \n',THRES+PACE);
%     fprintf('========================================================== \n');
% end


end