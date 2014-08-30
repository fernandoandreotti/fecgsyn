% FECGSYN is fruit of the collaboration between the Department of Engineering 
% Science, University of Oxford (DES-OX) and the Institute of Biomedical Engineering, 
% TU Dresden (IBMT-TUD). The authors are Joachim Behar (DES-OX), Fernando Andreotti 
% (IBMT-TUD), Julien Oster (DES-OX), Sebastian Zaunseder (IBMT-TUD) and 
% Gari Clifford (DES-OX). 
%
% The present user interface was contributed by Mohsan Alvi (DES-OX) under
% the supervision of Joachim Behar (DES-OX) and Fernando Andreotti (IBMT-TUD).
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

function [ param_struct ] = FECGSYN_UI_create_param_structs( varargin )
%CREATE_PARAM_STRUCTS Creates a all param structs for default scenarios
%   
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

disp('Creating default scenarios...');

% deal with optional inputs
if nargin < 2;    THR = 0.2;    else    THR = varargin{1};      end
if nargin < 3;    mVCG = 5;     else    mVCG = varargin{2};     end
if nargin < 4;    fVCG = 4;     else    fVCG = varargin{3};     end
if nargin < 5;    debug = 11;   else    debug = varargin{4};    end
if nargin < 6;    CH_CANC = 5;  else    CH_CANC = varargin{5};  end
if nargin < 7;    POS_DEV = 0;  else    POS_DEV = varargin{6};  end


%% Create default struct (copied from run_ecg_generator.m)
default_param = struct;


%% Initialise structures for default scenarios
% 8 default scenarios and one entry for the user-defined scenario
param_struct = cell(8,1);

for i=1:length(param_struct)
    param_struct{i} = default_param;
end

%% DS 1 - Simple run
param = param_struct{1};

    param.title = 'Simple';

    param.fs = 1000; % sampling frequency [Hz]
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

param_struct{1} = param;


%% DS 2 - Adding noise
param = param_struct{2};

    param.title = 'Noise';

    param.fs = 1000;
    param.ntype = {'MA'}; % noise types
    param.noise_fct = {1}; % constant SNR (each noise may be modulated by a function)
    param.noise_fct_str = {'1'};
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

param_struct{2} = param;

%% DS 3 - Respiration

param = param_struct{3};

    param.title = 'Respiration';

    param.fs = 1000;
    param.mres = 0.25; % mother respiration frequency
    param.fres = 0.8; % foetus respiration frequency
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

param_struct{3} = param;


%% DS 4 - Adding foetal movement

param = param_struct{4};

    param.title = 'Foetal movement';

    param.fs = 1000;
    param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

param_struct{4} = param;

%% DS 5 - Adding heart rate variability

param = param_struct{5};

    param.title = 'Heart rate variability';
    
    param.fhr = 135; param.mhr = 130;
    param.mtypeacc = 'nsr';
    param.ftypeacc = {'nsr'};

    param.fs = 1000;
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

param_struct{5} = param;

% %% DS 6 - Adding uterine contraction
% % simulating uterus activity (through muscular noise) and heart rate changes for both 
% % fetus (umbilical cord compression) and mother (acceleration followed by decelleration)
% param = param_struct{6};
% 
%     param.title = 'Uterine contraction';
% 
%     param.fs = 1000;
%     param.n = 60000;
%     if ~isempty(mVCG); param.mvcg = mVCG; end;
%     if ~isempty(fVCG); param.fvcg = fVCG; end;
%     if ~isempty(POS_DEV); param.posdev = 0; end;
% 
%     % Case 6 (early deceleration)
%     x = linspace(-param.n/10,param.n/10,param.n);
%     mu = 0;
% 
%     % Case 6b (late deceleration)
%      mu = 0.5; % distance from center-beginning [%]
%     param.maccmean = -mu;
% 
%     gauss = (100/(param.n*sqrt(2*pi)))*exp(-(x-(x(1)*mu)).^2/(2*(param.n/50)^2)); % approximating
%     gauss = gauss/max(gauss);                      % uterine contraction by gaussian modulated MA    
% 
%     param.ntype = {'MA'};   
%     param.noise_fct = {gauss};
%     param.SNRmn = -10;
% 
%     % heart rate changes
%     param.macc = 40;
%     param.mtypeacc = 'gauss';
%     param.facc = -30;
%     param.fhr = 130;
%     param.ftypeacc = {'mexhat'};
%     param.faccstd{1} = 0.5;
% 
% param_struct{6} = param;

%% DS 7 - Adding ectopic beats

param = param_struct{6};

    param.title = 'Ectopic beats';

    param.fs = 1000; % sampling frequency [Hz]
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;
    param.mectb = 1; param.fectb = 1; 

param_struct{6} = param;

%% DS 8 - Multiple pregnancies

param = param_struct{7};

    param.title = 'Multiple pregnancies';

    param.fs = 1000;
    param.fhr = [120 150];
    param.fres = [0 0];
    param.facc = [0 0];
    param.ftypeacc = {'none' 'none'};
    if ~isempty(POS_DEV); param.posdev = 0; end;
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = [fVCG (fVCG+3)]; end;

    param.fheart{1} = [-pi/10 0.35 -0.1];
    param.fheart{2} = [pi/10 0.4 -0.2];

param_struct{7} = param;

%% Custom Scenario
param_struct{8} = param_struct{1};
param_struct{8}.title = 'Custom';

%% Add missing values
for i = 1:length(param_struct)
    param_struct{i} = FECGSYN_UI_add_default_params(param_struct{i}, varargin);
end

disp('Default scenarios ready.');


end

