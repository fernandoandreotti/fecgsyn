function [ recg_out, param ] = eval_scenario( choice , varargin)
%CHOOSE_SCENARIO Runs scenarios of run_ecg_generator 
%   INPUTS:
%       choice (integer) - corresponds to the scenario choice
%       (optional) THR - threshold of QRS detector
%       (optional) mVCG - mother VCG (if empty then the simulator randomly choose one within the set of available VCGs)
%       (optional) fVCG - foetus VCG (ibid)
%       (optional) debug - debug level
%       (optional) CH_CANC - channel onto which to perform MECG cancellation
%       (optional) POS_DEV - slight deviation from default hearts and electrodes positions 
%                            (0: hard coded values, 1: random deviation and phase initialisation)
%
%   OUTPUTS:
%       recg_out (struct) - the output of run_ecg_generator()
%
%
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

% deal with optional inputs
if nargin < 2;    THR = 0.2;    else    THR = varargin{1};      end
if nargin < 3;    mVCG = 5;     else    mVCG = varargin{2};     end
if nargin < 4;    fVCG = 4;     else    fVCG = varargin{3};     end
if nargin < 5;    debug = 11;   else    debug = varargin{4};    end
if nargin < 6;    CH_CANC = 5;  else    CH_CANC = varargin{5};  end
if nargin < 7;    POS_DEV = 0;  else    POS_DEV = varargin{6};  end

% This function will always be called by the gui, so the gui flag is true
param.gui = 1;

if choice == 10
    % Simple run
    disp('---- Example (1): SIMPLE RUN ----');
    param.fs = 1000; % sampling frequency [Hz]

    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    recg_out = run_ecg_generator(param,debug);

elseif choice == 20
    % Adding noise
    disp('---- Example (2): ADDING NOISE ----');
    param.fs = 1000;
    param.ntype = {'MA'}; % noise types
    param.noise_fct = {1}; % constant SNR (each noise may be modulated by a function)
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    recg_out = run_ecg_generator(param,debug);
        
elseif choice == 30
    % Noise + Respiration
    disp('---- Example (3): ADDING RESPIRATION ----');
    param.fs = 1000;
    param.mres = 0.25; % mother respiration frequency
    param.fres = 0.8; % foetus respiration frequency
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    recg_out = run_ecg_generator(param,debug);

elseif choice == 40
    % Adding foetal movement
    disp('---- Example (4): ADDING FOETAL MOVEMENT ----');
    param.fs = 1000;
    param.ftraj{1} = 'helix'; % giving spiral-like movement to fetus
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    recg_out = run_ecg_generator(param,debug);

elseif choice == 50
    % Adding heart rate variability
    disp('---- Example (5): ADDING HEART RATE VARIABILITY ----');
    param.fs = 1000;

    % Case 4
    param.fhr = 130; param.mhr = 130;

    % Case 5
    % param.macc = 40; % maternal acceleration in HR [bpm]
    % param.mtypeacc = 'tanh';
    % param.mhr = 110;
    % param.fhr = 140;

    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    recg_out = run_ecg_generator(param,debug);
    
elseif choice == 60
    % Adding uterine contraction
    %
    % simulating uterus activity (through muscular noise) and heart rate changes for both 
    % fetus (umbilical cord compression) and mother (acceleration followed by decelleration)
    disp('---- Example (6): ADDING UTERINE CONTRACTION ----');
    param.fs = 1000;
    param.n = 60000;
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;

    % Case 6 (early deceleration)
    x = linspace(-param.n/10,param.n/10,param.n);
    mu = 0;

    % Case 6b (late deceleration)
     mu = 0.5; % distance from center-beginning [%]
    param.maccmean = -mu;

    gauss = (100/(param.n*sqrt(2*pi)))*exp(-(x-(x(1)*mu)).^2/(2*(param.n/50)^2)); % approximating
    gauss = gauss/max(gauss);                      % uterine contraction by gaussian modulated MA    

    param.ntype = {'MA'};   
    param.noise_fct = {gauss};
    param.SNRmn = -10;

    % heart rate changes
    param.macc = 40;
    param.mtypeacc = 'gauss';
    param.facc = -30;
    param.fhr = 130;
    param.ftypeacc = {'mexhat'};
    param.faccstd{1} = 0.5;
    
    recg_out = run_ecg_generator(param,debug);
    
elseif choice == 70
    % Adding ectopic beats
    disp('---- Example (7): ADDING ECTOPIC BEATS ----');

    param.fs = 1000; % sampling frequency [Hz]
    if ~isempty(mVCG); param.mvcg = mVCG; end;
    if ~isempty(fVCG); param.fvcg = fVCG; end;
    if ~isempty(POS_DEV); param.posdev = 0; end;
    param.mectb = 1; param.fectb = 1; 

    recg_out = run_ecg_generator(param,debug);
    
elseif choice == 80
    % Multiple pregnancies
    disp('---- Example (8): MULTIPLE PREGNANCIES ----');

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

    recg_out = run_ecg_generator(param,debug);
    
else
    % Defaults to the simplest example
    [recg_out, param] = choose_scenario( 10 , varargin);
end

