function [qt_test,qt_ref,tqrs_test,tqrs_ref,qt_err,theight_err,numNaN]=...
    FECGSYN_morpho_loop(fecg,residual,fqrs,fs,SAMPS,fname,filterc,debug)
% function [qt_test,qt_ref,tqrs_test,tqrs_ref,qt_err,theight_err,numNaN]= FECGSYN_morpho_loop(fecg,residual,fqrs,fs,SAMPS,fname,filterc,debug)
% Function to perform morphological analysis for TS/BSS extracted data
% This function loops through available channels on reference and test
% signals and generates morphological features (if possible) using the
% ECGPUWAVE algorithm.
%
% Inputs:
% fecg         Propagated fetal signal before mixture with noise sources
% residual     Result of fetal extraction from abdominal signals
% fqrs         Reference fetal QRS samplestamps
% SAMPS        Number of samples used for generating templates
% fname        Filename to be used in saving plots
% filterc      Filter coefficients [b_hp,a_hp,b_lp,a_lp] being
%               highpass (hp) and lowpass (lp)
% debug        toggle for debugging
%
% Outputs:
% qt_err        Array containing QT error for each template
% theight_err   Array containing T-height error for each template
%
% Reference:
% Jane, R et al. 1997 Evaluation of An Automatic Threshold Based Detector of Waveform Limits in Holter ECG
% with the QT Database. Computers in Cardiology. Moody, G B, Mark, R G, Zoccola, A and Mantero, S
%
% Examples:
% TODO
%
% See also:
% FECGSYN_benchMorph
% FECGSYN_manalysis
% FECGSYN_QTcalc
%
% --
% fecgsyn toolbox, version 1.2, Jan 2017
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% University of Oxford, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
%
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
%
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D.
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings.
%   Physiological Measurement.35 1537-1550. 2014.
%
% Last updated : 10-03-2016
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

numNaN = 0;
% Allocatting
qt_test = cell(size(residual,1),length(residual)/SAMPS,1);
qt_ref = qt_test;
tqrs_test = qt_test;
tqrs_ref = qt_test;
qt_err = qt_test;
theight_err = qt_test;
%= Block-wise calculation and template generation
block = 1;
for j = 1:SAMPS:length(residual)
    for ch = 1:size(fecg,1)
        % checking borders
        if j+SAMPS > length(residual)
            endsamp = length(residual);
        else
            endsamp = j + SAMPS -1;
        end
        % qrs complexes in interval
        qrstmp = fqrs(fqrs>j&fqrs<endsamp)-j;
        %% Template Generation
        % reference template
        [temp_ref,qrs_ref,status2] = FECGSYN_tgen(fecg(ch,j:endsamp),qrstmp,fs,debug);
        % abdominal signal template
        if ch <= size(residual,1)
            [temp_abdm,qrs_abdm,status1] = FECGSYN_tgen(residual(ch,j:endsamp),qrstmp,fs,debug);
        else % usually relevant for ICA cases, where number of components is smaller than the number of input channels
            temp_abdm = temp_ref;
            qrs_abdm = qrs_ref;
            status2 = status1;
        end
        
        temp_abdm = temp_abdm.avg; temp_ref = temp_ref.avg;
        % crop end of templates which have steps on them
        try
            per80 = round(0.8*length(temp_abdm));
            [~,idx]=findpeaks(abs(diff(temp_abdm(per80:end))),'Threshold',10*median(abs(diff(temp_abdm(per80:end)))));
            if ~isempty(idx)
                idx = idx-1;
                med1 = median(temp_abdm(per80:per80+idx)); med2 = median(temp_abdm(per80+idx:end));
                temp_abdm(per80+idx:end) = temp_abdm(per80+idx:end)+(med1-med2); % removing step in signals
            end
            clear idx med1 med2 per 80
            per80 = round(0.8*length(temp_ref));
            [~,idx]=findpeaks(abs(diff(temp_ref(per80:end))),'Threshold',10*median(abs(diff(temp_ref(per80:end)))));
            if ~isempty(idx)
                idx = idx-1;
                med1 = median(temp_ref(per80:per80+idx)); med2 = median(temp_ref(per80+idx:end));
                temp_ref(per80+idx:end) = temp_ref(per80+idx:end)+(med1-med2); % removing step in signals
            end
            clear idx med1 med2 per 80
        catch
            disp('Template could not be generated for segment. Skipping..')
        end
        
        if (~status1||~status2)
            qt_test{ch,block} = NaN;
            qt_ref{ch,block} = NaN;
            tqrs_test{ch,block} = NaN;
            tqrs_ref{ch,block} = NaN;
            qt_err{ch,block} = NaN;
            theight_err{ch,block} = NaN;
        else
            %% Performs morphological analysis
            [qt_ref{ch,block},qt_test{ch,block},tqrs_ref{ch,block},tqrs_test{ch,block}] = FECGSYN_manalysis(temp_abdm,temp_ref,qrs_abdm,qrs_ref,fs,filterc,fname,debug);
        end
        % Saves generated plots
        try
            if debug && ~isnan(qt_test{ch,block}) && ~isnan(qt_ref{ch,block})
                
                drawnow
                subplot(2,1,1)
                hold on
                text(0,0,['QT = ' strcat(num2str(qt_ref{ch,block}))])
                
                subplot(2,1,2)
                hold on
                text(0,0,['QT = ' strcat(num2str(qt_test{ch,block}))])
                print('-dpng','-r72',[fname '_ch' num2str(ch) '_s' num2str(block) '.png'])
            end
        catch
            warning('Failed to save plot')
        end
    end
    block = block+1;
end
end
