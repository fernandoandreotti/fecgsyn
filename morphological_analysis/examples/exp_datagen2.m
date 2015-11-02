%% Another example of data generation
%  Here are the specifics
% - 5 feto-maternal combinations
% - 3 SNR levels (3,6,9 dB)
% - Each dataset with 1,5min duration
% - 5x repetition for statistics
%
% * Cases/events:
% - Case 0 - Baseline
% - Case 1 - HR abrupt change (by 1/3 using tanh() normally distributed)
% - Case 2 - SNR abrupt change (by 1/3 using tanh() modulation, amplitude and direction normally distributed)
% - Case 3 - overall ECG amplitude change (sinusoidal 1-10 cycles/recording)
% - Case 4 - overall ECG amplitude change (skewed with Gamma distribution, 1-10 cycles/recording, amplitude, width and direction normally dist.)
% - Case 5 - T+P waves + amplitude change (sinusoidal)
% - Case 6 - T+P waves + amplitude change (skewed)
% - Case 7 -T+P waves + amplitude change (skewed) + ECG amplitude change (sinusoidal)
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-11-2015
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
function exp_datagen2()

path = cd;
debug = 5;

% Generating simulated data with various SNR for morphological analysis
% global parameters
paramorig.fs = 1000;            % sampling frequency [Hz]
paramorig.n = 90*paramorig.fs;  % number of data points to generate (5 min)

% electrode positions
x = pi/12*[3 4 5 6 7 8 9 10]' -pi/2;     % 32 abdominal channels
y = .5*ones(8,1);
xy = repmat([x y],4,1);
z = repmat([-.1 -.2 -.3 -.4],8,1); z = reshape(z,32,1);
abdmleads = [xy z];
refs = [-pi/4 0.5 0.4;(5/6-.5)*pi 0.5 0.4];  % + 2 reference leads
paramorig.elpos = [abdmleads;refs];
cd(path)
for i = 1:5           % generate 5 cases of each
    close all
    paramst = paramorig;
    paramst.fhr = 135+25*randn;   % choosing foetal heart rate
    % mean=135, std= 25 [bpm]
    paramst.mhr = 80+20*randn;    % choosing maternal heart rate
    % mean = 80, std = 20 [bpm]
    
    % setting up stationary mixture
    paramst.mtypeacc = 'nsr';      % force constant mother heart rate
    paramst.ftypeacc = {'nsr'};    % force constant foetal heart rate
    paramst.SNRfm = -7 + 2*randn;
    out = run_ecg_generator(paramst,debug);  % stationary output
    %plotmix(out)
    out = clean_compress(out);
    paramst = out.param;                    % keeping same parameters
    clear out
    % adding some noise
    for SNRmn = 3:3:9 % five noise levels
        for loop = 1:5 % repeat same setup
            % just recalculating noise five times
            % reseting config    outst = out;
            %% Case 0: Baseline (noise and hearts, no event)
            disp(['Generating for SNRmn=' num2str(SNRmn) ' simulation number ' num2str(i) '.'])
            param = paramst;
            param.SNRmn = SNRmn;    % varying SNRmn
            param.ntype = {'MA','MA'}; % noise types
            disp(['Generating for SNRmn=' num2str(SNRmn) ' simulation number ' num2str(i) '.'])
            param = paramst;
            param.SNRfm = -9 + 2*randn;
            param.SNRmn = SNRmn;    % varying SNRmn
            param.ntype = {'MA','MA'}; % noise types
            param.fheart{1} = [pi*(2*rand-1)*sign(randn)/10 (0.1*rand*sign(randn)+0.15) -0.3*rand]; % define first foetus position
            param.noise_fct = {1+.5*randn,1+.5*randn}; % constant SNR (each noise may be modulated by a function)
            param.mres = 0.25 + 0.05*randn; % mother respiration frequency
            param.fres = 0.9 + 0.05*randn; % foetus respiration frequency
            parambase = param;
            out = run_ecg_generator(param,debug);  % stationary output
            out = clean_compress(out);
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c0',i,SNRmn,loop)],'out')
            clear out
% %             
% %             %% Case 1: rate rate accelerations
% %             param = parambase;
% %             param.macc = (20+10*abs(randn))*sign(randn); % maternal acceleration in HR [bpm]
% %             param.mtypeacc = 'tanh';                % hyperbolic tangent acceleration
% %             out = run_ecg_generator(param,debug);   % stationary output
% %             out = clean_compress(out);
% %             save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c1',i,SNRmn,loop)],'out')
% %             clear out
% %             
% %             %% Case 2: SNR abrupt change
% %             param = parambase;
% %             param.noise_fct{1} = 1+sign(randn)*(rand+0.3)*tanh(linspace(-pi,2*pi,param.n));  % tanh function           
% %             param.noise_fct{2} = param.noise_fct{1};  % tanh function           
% %             param.ntype = {'MA' 'MA'};
% %             param.SNRmn = -6;         % put additional contraction with strong power
% %             out = run_ecg_generator(param,debug);  % stationary output
% %             out = clean_compress(out);
% %             save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c2',i,SNRmn,loop)],'out')
% %             clear out
            %% Case 3: overall ECG amplitude change (sinusoidal 1-10 cycles/recording)
            param = parambase;
            out = run_ecg_generator(param,debug);  % stationary output
            cyccount = randi([1,10],1,1);
            piinit = (2*rand-1)*pi; % [-pi,pi]
            modfun = (1+sin(linspace(piinit,cyccount*2*pi+piinit,param.n))*(0.2*rand+0.001));
            out.mecg = repmat(modfun,34,1).*out.mecg;
            out = clean_compress(out);
            out.modfun = modfun;
            out.cycles = cyccount;
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c3',i,SNRmn,loop)],'out')
            clear out
            %% Case 4: previous case skewed      
            param = parambase;
            samps = param.n/cyccount;
            epsilon = round(samps/2);
            w = 0.3*rand*samps;
            alpha = 5*randn;
            skewgauss = @(x) (1/(2*w*pi))*exp(-(x-epsilon)^2/(2*w)^2)*(1+erf((alpha*(x-epsilon)/w)/sqrt(2)));
            modfun = arrayfun(skewgauss,1:samps);
            modfun = repmat(modfun,1,cyccount);
            modfun = 1+(rand-rand)*(modfun./max(modfun));
            out = run_ecg_generator(param,debug);  % stationary output           
            out.mecg = repmat(modfun,34,1).*out.mecg;
            out = clean_compress(out);
            out.cycles = cyccount;
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c4',i,SNRmn,loop)],'out')
            clear out
            %% Case 5 - T+P waves + amplitude change (sinusoidal)
            param = parambase;
            out = run_ecg_generator(param,debug);  % stationary output           
            for b = 1:length(out.mqrs)
               size1 = randn*0.1*out.fs;
               size2 = randn*0.2*out.fs;
               modfun1 = sin(linspace(0,pi,size1));
               modfun2 = sin(linspace(0,pi,size2));
                
            end
            
            
        end
    end
end
end
% this function eliminates some of the substructures from "out" and
% compresses the variables to int16 for saving disk space
function out=clean_compress(out)
gain = 3000;
out_tmp=rmfield(out,{'f_model' 'm_model' 'vols' 'selvcgm' 'selvcgf'});
out = struct();
out.mecg = int16(round(3000*out_tmp.mecg));
if ~isempty(out_tmp.fecg)
    for i = 1:length(out_tmp.fecg)
        out.fecg{i} = int16(round(3000*out_tmp.fecg{i}));
    end
else
    out.fecg = {};
end
if ~isempty(out_tmp.noise)
    for i = 1:length(out_tmp.noise)
        out.noise{i} = int16(round(gain*out_tmp.noise{i}));
    end
else
    out.noise = {};
end
out.mqrs = out_tmp.mqrs;
out.fqrs = out_tmp.fqrs;
out.param = out_tmp.param;
end