%% Script for generating data for stress testing spatial filtering techniques
% 
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing 
% effects, foetal movement etc.).
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
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
function generate_data(path)

% clear all; close all;
debug = 0;

%% Generating simulated data with various SNR for morphological analysis
% global parameters
paramorig.fs = 1000;            % sampling frequency [Hz]
paramorig.n = 60*paramorig.fs;  % number of data points to generate (5 min)
paramorig.SNRfm = -6;          % fetal SNR [dB]    

% electrode positions
x = pi/12*[3 4 5 6 7 8 9 10]' -pi/2;     % 32 abdominal channels 
y = .5*ones(8,1);
xy = repmat([x y],4,1);
z = repmat([-.1 -.2 -.3 -.4],8,1); z = reshape(z,32,1);
abdmleads = [xy z];
refs = [-pi/4 0.5 0.4;(5/6-.5)*pi 0.5 0.4];  % + 2 reference leads
paramorig.elpos = [abdmleads;refs];
cd(path)
for i = 1:10            % generate 5 cases of each
        close all
        paramst = paramorig;
        paramst.fhr = 135+25*randn;   % choosing foetal heart rate
                                      % mean=135, std= 25 [bpm]
        paramst.mhr = 80+20*randn;    % choosing maternal heart rate
                                      % mean = 80, std = 20 [bpm]
        
        %% stationary mixture
        paramst.mtypeacc = 'none';      % force constant mother heart rate
        paramst.ftypeacc = {'none'};    % force constant foetal heart rate
        out = run_ecg_generator(paramst,debug);  % stationary output
        plotmix(out)
        save([path 'fecgsyn' sprintf('%2.2d',i)],'out')
        paramst = out.param;                     % keeping same parameters
            
    %% adding some noise    
    for SNRmn = -3:3:12 % six noise levels
        for loop = 1:5 % repeat same setup 
                       % just recalculating noise five times
            % reseting config
            disp(['Generating for SNRmn=' num2str(SNRmn) ' simulation number ' num2str(i) '.'])
            param = paramst;
            param.SNRmn = SNRmn;    % varying SNRmn
            param.ntype = {'MA','EM'}; % noise types
            param.noise_fct = {1+.5*randn,1+.5*randn,1+.5*randn}; % constant SNR (each noise may be modulated by a function)
            param.mres = 0.25 + 0.05*randn; % mother respiration frequency
            param.fres = 0.9 + 0.05*randn; % foetus respiration frequency
            out = run_ecg_generator(param,debug);  % stationary output
            plotmix(out)
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d',i,SNRmn,loop)],'out')
            paramnst = out.param;
            %% non-stationary mixture
            % Case 1: foetal movement
            param = paramnst;
            param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
            
            out = run_ecg_generator(param,debug);  % stationary output
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c1',i,SNRmn,loop)],'out')
            % Case 2: rate rate accelerations
            param = paramnst;
            param.macc = 20; % maternal acceleration in HR [bpm]
            param.mtypeacc = 'tanh'; % hyperbolic tangent acceleration
            param.facc = -40; % foetal decceleration in HR [bpm]
            param.ftypeacc = {'mexhat'}; % gaussian drop and recovery
            
            out = run_ecg_generator(param,debug);  % stationary output
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c2',i,SNRmn,loop)],'out')
            % Case 3: contraction
            param = paramnst;
            x = linspace(-param.n/10,param.n/10,param.n);
            mu = 0;
            gauss = (100/(param.n*sqrt(2*pi)))*exp(-(x-(x(1)*mu)).^2/(2*(param.n/50)^2)); % approximating
            gauss = gauss/max(gauss);                      % uterine contraction by gaussian modulated MA
            param.noise_fct{1} = gauss;
            param.macc = 40;
            param.mtypeacc = 'gauss';
            param.facc = -30;
            param.ftypeacc = {'mexhat'};
            param.faccstd{1} = 0.5;
            
            out = run_ecg_generator(param,debug);  % stationary output
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c3',i,SNRmn,loop)],'out')
            % Case 4: ectopic beats
            param = paramnst;
            param.mectb = 1; param.fectb = 1;
            
            out = run_ecg_generator(param,debug);  % stationary output
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c4',i,SNRmn,loop)],'out')
            % Case 5: twins
            param = paramnst;
            param.fhr(2) = 135+25*randn;
            param.fres(2) = 0.9 + 0.05*randn;
            param.facc = [0 0];
            param.ftypeacc = {'none' 'none'};
            param.fheart{2} = [pi/10 0.4 -0.2];
            param.fvcg(2) = randi([1,9]);
            
            out = run_ecg_generator(param,debug);  % stationary output
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c5',i,SNRmn,loop)],'out')
        end
    end
end
