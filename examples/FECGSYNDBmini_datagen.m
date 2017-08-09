function FECGSYNDBmini_datagen(path,debug,varargin)
%function FECGSYN_generate_data(path,debug)
% Script used for generating data in FECGSYNDB (see https://www.physionet.org/physiotools/ipmcode/fecgsyn/)
% 
%
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing
% effects, foetal movement etc.).
%
% Input:
%   path        saving path (default pwd)
%   debug       toggle debug (default true)
% 
%
% More detailed help is in the <a href="https://fernandoandreotti.github.io/fecgsyn/">FECGSYN website</a>.
%
% Examples:
% FECGSYNDB_datagen(pwd,5) % generate data and plots
%
% See also:
% exp_datagen1
% exp_datagen2 
% FECGSYNDB_datagen
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

%% == check inputs
if nargin >2, error('Too many inputs to data generation function'),end
slashchar = char('/'*isunix + '\'*(~isunix));
optargs = {[pwd slashchar] 5};  % default values for input arguments
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);
[path,debug] = optargs{:};
if ~strcmp(path(end),slashchar), path = [path slashchar];end

% Generating simulated data with various SNR for morphological analysis
% global parameters
paramorig.fs = 1000;            % sampling frequency [Hz]
paramorig.n = 30*paramorig.fs;  % number of data points to generate (20 sec)

% electrode positions
x = pi/12*[3 4 5 6 7 8 9 10]' -pi/2;     % 32 abdominal channels
y = .5*ones(8,1);
xy = repmat([x y],4,1);
z = repmat([-.1 -.2 -.3 -.4],8,1); z = reshape(z,32,1);
abdmleads = [xy z];
refs = [-pi/4 0.5 0.4;(5/6-.5)*pi 0.5 0.4];  % + 2 reference leads
paramorig.elpos = [abdmleads;refs];
cd(path)
    close all
    paramst = paramorig;
    paramst.fhr = 135+25*randn;   % choosing foetal heart rate
    % mean=135, std= 25 [bpm]
    paramst.mhr = 80+20*randn;    % choosing maternal heart rate
    % mean = 80, std = 20 [bpm]
    
    %% stationary mixture
    paramst.mtypeacc = 'nsr';      % force constant mother heart rate
    paramst.ftypeacc = {'nsr'};    % force constant foetal heart rate
    paramst.SNRfm = -9 + 2*randn;
    out = run_ecg_generator(paramst,debug);  % stationary output
    %plotmix(out)
    out = clean_compress(out);
    paramst = out.param;                    % keeping same parameters
    clear out
    %% adding some noise
    i = 1;
    for SNRmn = [6 12] % five noise levels
        for loop = 1:3 % repeat same setup
            % just recalculating noise five times
            % reseting config    outst = out;

            % Case 0: FECG + MECG + noise
            disp(['Generating for SNRmn=' num2str(SNRmn) ' simulation number ' num2str(i) '.'])
            param = paramst;
            param.SNRmn = SNRmn;    % varying SNRmn
            param.ntype = {'MA','MA'}; % noise types
            disp(['Generating for SNRmn=' num2str(SNRmn) ' simulation number ' num2str(i) '.'])
            param = paramst;
            param.SNRfm = -9 + 2*randn;
            param.SNRmn = SNRmn;    % varying SNRmn
            param.ntype = {'MA','MA'}; % noise types
            param.fheart{1} = [pi*(2*rand-1)/10 (0.1*rand+0.25) -.2*(1+rand)]; % define first foetus position
            param.noise_fct = {1+.5*randn,1+.5*randn}; % constant SNR (each noise may be modulated by a function)
            param.mres = 0.25 + 0.05*randn; % mother respiration frequency
            param.fres = 0.9 + 0.05*randn; % foetus respiration frequency
            out = run_ecg_generator(param,debug);  % stationary output
            %plotmix(out)
            out = clean_compress(out);                       
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c0',i,SNRmn,loop)],'out')

            out_noise = out;
            param_noise = out.param;
            param_noise = rmfield(param_noise,{'ntype' 'noise_fct'});         % will not be re-simulated
            param_noise.posdev = 0;    % maternal and fetal hearts fixed
            
            % Case 0: baseline, no noise added
            out.noise = [];
            out.param = rmfield(out.param,{'ntype' 'noise_fct'});
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d',i,SNRmn,loop)],'out')
            clear out
            
            % Case 1: foetal movement
            % no need to simulate again
            param = param_noise;
            param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
            out = run_ecg_generator(param,debug);  % stationary output
            out = clean_compress(out);
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c1',i,SNRmn,loop)],'out')
            clear out
            
            % Case 2: rate rate accelerations
            param = param_noise;
            param.macc = (20+10*randn)*sign(randn); % maternal acceleration in HR [bpm]
            param.mtypeacc = 'tanh';                % hyperbolic tangent acceleration
            param.facc = (30 + 10*randn)*sign(randn); % foetal decceleration in HR [bpm]
            param.ftypeacc = {'mexhat'};            % gaussian drop and recovery
            out = run_ecg_generator(param,debug);   % stationary output
            out = clean_compress(out);
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c2',i,SNRmn,loop)],'out')
            clear out
            
            % Case 3: contraction
            param = param_noise;
            x = linspace(-param.n/10,param.n/10,param.n);
            mu = 0;
            gauss = (100/(param.n*sqrt(2*pi)))*exp(-(x-(x(1)*mu)).^2/(2*(param.n/50)^2)); % approximating
            gauss = gauss/max(gauss);                      % uterine contraction by gaussian modulated MA
            param.noise_fct{1} = gauss;
            param.ntype = {'MA'};
            param.SNRmn = -6;         % put additional contraction with strong power
            param.macc = 40;
            param.mtypeacc = 'gauss';
            param.facc = -30;
            param.ftypeacc = {'mexhat'};
            param.faccstd{1} = 0.5;
            out = run_ecg_generator(param,debug);  % stationary output
            out = clean_compress(out);
            out.noise = [out.noise out_noise.noise];
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c3',i,SNRmn,loop)],'out')
            clear out
            
            % Case 4: ectopic beats
            param = param_noise;
            param.mectb = 1; param.fectb = 1;
            out = run_ecg_generator(param,debug);  % stationary output
            out = clean_compress(out);
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c4',i,SNRmn,loop)],'out')
            clear out
            
            % Case 5: twins
            param = param_noise;
            param.fhr(2) = 135+25*randn;
            param.fres(2) = 0.9 + 0.05*randn;
            param.fheart{2} = [pi*(2*rand-1)/10 (0.1*rand+0.25) -(0.2+.2*rand)];
            param.fvcg(2) = randi([1,9]);
            param=rmfield(param,{'faccmean' 'facc' 'ftypeacc' 'faccstd' 'ftraj'});
            out = run_ecg_generator(param,debug);  % stationary output
            out = clean_compress(out);
            out.noise = out_noise.noise;    % re-inserting noise
            save([path 'fecgsyn' sprintf('%2.2d_snr%2.2ddB_l%d_c5',i,SNRmn,loop)],'out')
            clear out
        end
    end
end

