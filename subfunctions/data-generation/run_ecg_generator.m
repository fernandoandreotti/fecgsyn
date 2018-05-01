function out = run_ecg_generator(param,debug)
% function out = run_ecg_generator(param,debug)
% generate a realistic FECG-MECG mixture. The physiological phenomenon that
% are modelled include: heart rate variability and sudden change in hr,
% rotation of the heart axis with respect to breathing rate, translation
% of the foetal heart modelling foetal movement and realistic noise.
% Note that in order to make repeated simulation more variable (and so more
% representative of the reality) the positions of the maternal and foetal
% hearts are randomly placed around the default coordinates. Note that the
% coordinate system is normalised and centred on the cylinder (which is modeling
% the volume conductor). As a consequence, electodes location, hearts
% location are defined relative to each other (relative coordinate system).
%
%
% list of abbreviation used in the toolbox:
%   ECG:        electrocardiogram
%   MECG:       maternal ECG
%   FECG:       foetal ECG
%   NI-FECG:    non invasive FECG
%   HR:         heart rate
%   HRV:        heart rate variability
%   FHR:        foetal heart rate
%   MHR:        maternal heart rate
%   RSA:        respiratory sinus arrhythmia
%
% inputs:
%   param: structure containing all the information needed for running the
%   simulation.
%       param.mheart    maternal heart origin - actual location will be
%                       picked randomly around it* -(default [2*pi/3 0.4
%                       0.4]) [normalised]
%       param.fheart    foetus heart origin - actual location will be
%                       picked randomly around it* - (default [-pi/10 0.4
%                       -0.3]) [normalised]
%       param.elpos:    electrode pair locations in polar coordinate [normalised]
%       param.n:        number of datapoint to generate (default 15000)
%       param.fs:       sampling frequency (default 1000) [Hz]
%       param.ntype:    noise type (default MA) [string]
%       param.noise_fct function of modulating noise (each noise may be modulated by a
%                       function, e.g param.noise_fct=sin(linspace(-2*pi,2*pi,param.n)))
%       param.SNRn:     SNR (MECG+FECG)/Noise (default 6)
%       param.SNRm:     SNR FECG/MECG (default -10)
%       param.mhr:      mother reference HR (default 110) [bpm]
%       param.fhr:      foetus reference HR (default 150) [bpm]
%       param.macc:     maternal acceleration in HR
%       param.facc:     foetus acceleration in HR
%       param.mtypeacc: maternal acceleration type (chosen from switch inside function,
%                       e.g. 'none', 'mexhat', 'gauss' or 'flattop')
%       param.ftypeacc: foetus acceleration type (chosen from switch inside function,
%                       e.g. 'none', 'mexhat', 'gauss' or 'flattop')
%       param.ftraj:    trajectory given to fetus heart (e.g. 'none', 
%                       'step','linear', 'spline' or 'spiral')
%       param.mtraj     movement trajectory to maternal heart (e.g. case of
%                       maternal movement). Analogous to fetal trajectory
%       param.fname:    record name for saving output (default 'aecg') if
%                       empty then no file is saved
%       param.mres:     respiratory frequency of mother (default 0) [Hz]
%       param.fres:     respiratory frequency of foetus (default 0) [Hz]
%       param.mvcg:     mother vcg chosen (1-9)
%       param.fvcg:     foetus vcg chosen (1-9)
%       param.evcg:     ectopic beat params (1-4)
%       param.posdev:   position deviation (bool). 1: the position of the
%                       electrodes and hearts are slightly varied around
%                       their specified or default positions. 0: no variation.
%       param.mectb:    add ectopic beats to maternal ECG (bool)
%       param.fectb:    add ectopic beats to foetal ECG (bool)
%   debug:          debug level (1-5), (default: 1)
%       1:              Some generated AECG
%       2:              VCG plots
%       3:              Projected FECG and MECG before being mixed
%       4:              Volume conductor
%       5:              Heart rate
% 
% * This is in order to be able to produce many simulations with the heart
% position changing at every iteration without having to respecify a
% location. Locations are specified in polar coordinate because it is
% easier to visualise.
%
% [normalised]: belly represented as a cylinder of radius 0.5 (so diameter
% 1) and height 1. The cylinder is centered on zero.
%
% NOTE: in the case only one input is entered, the function considers it to
% be an 'out' structure containing all the saved parameters from a previous
% simulation. The function only plots them.
%
% outputs:
%   out: outputs of the run_ecg_generator function containing the ecg
%   mixture and all the important model information that would allow
%   reproducing the results. 
%       out.mixture: generated ecg mixture [NB_EL x param.n matrix]
%       out.m_model: structure contaning dipole model for the foetus [struct]
%               m_model.H: Dower-like matrix for dipole either 2D (time invariant) or 3D
%                          (variant case).
%               m_model.VCG: VCG for dipole
%               m_model.type: maternal (1) or foetal (2) dipole
%       out.f_model: structure contaning dipole model for the foetus [struct]
%               ibid m_model
%
%       out.mecg:   mecg projected signal
%       out.fecg:   fecg projected signal
%       out.vols:   contains volume conductor information (electrodes and heart
%                   position)
%       out.mqrs:   maternal reference QRS
%       out.fqrs:   foetal reference QRS
%       out.param:  parameters used in the simulation [struct]
%       selvcgm:    selected maternal vcg [cell]
%       selvcgf:    selected foetal vcg [cell]
%
%
% Examples:
% TODO
%
% See also:
% exp_datagen1
% exp_datagen2
% exp_datagen3
% FECGSYNDB_datagen
% 
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

% == the two subfolders are assumed to exist
addpath(genpath('subfunctions'));
addpath(genpath('noise_sources'));

% == if we just want to plot using the inputed parameters
if nargin == 1
    debug_plots(param,5);
    return;
end

f_handles = [];
noise_f_handles = [];

% == default parameters
if ~any(strcmp('mheart',fieldnames(param))); param.mheart = [2*pi/3 0.2 0.4]; end;
if ~any(strcmp('fheart',fieldnames(param))); param.fheart{1} = [-pi/10 0.35 -0.3]; end;
if ~any(strcmp('elpos',fieldnames(param))); x = pi/12*[3 4 5 6 7 8 9 10]' -pi/2;     % 32 abdominal channels 
    y = .5*ones(8,1); xy = repmat([x y],4,1); z = repmat([-.1 -.2 -.3 -.4],8,1); z = reshape(z,32,1);
    abdmleads = [xy z]; refs = [-pi/4 0.5 0.4;(5/6-.5)*pi 0.5 0.4];  % + 2 reference leads
    param.elpos = [abdmleads;refs]; 
end   
if ~any(strcmp('refpos',fieldnames(param))); param.refpos = [pi 0.5 -0.3];end;
NB_FOETUSES = size(param.fheart,2); % number of foetuses figured out from the number of foetal heart locations entered
if ~any(strcmp('n',fieldnames(param))); param.n = 60000; end;
if ~any(strcmp('fs',fieldnames(param))); param.fs = 1000; end;
if ~any(strcmp('ntype',fieldnames(param))); param.ntype = ''; end;
if ~any(strcmp('noise_fct',fieldnames(param))); param.noise_fct(1:length(param.ntype)) = {1}; end;
if ~any(strcmp('SNRfm',fieldnames(param))); param.SNRfm = -9; end;
if ~any(strcmp('SNRmn',fieldnames(param))); param.SNRmn = 10; end;
if ~any(strcmp('mhr',fieldnames(param))); param.mhr = 90; end;
if ~any(strcmp('fhr',fieldnames(param))); param.fhr = repmat(150,NB_FOETUSES,1); end;
if ~any(strcmp('macc',fieldnames(param))); param.macc = 0; end;
if ~any(strcmp('facc',fieldnames(param))); param.facc = zeros(1,NB_FOETUSES); end;
if ~any(strcmp('mtypeacc',fieldnames(param))); param.mtypeacc = 'nsr'; end;
if ~any(strcmp('maccmean',fieldnames(param))); param.maccmean = 0; end;
if ~any(strcmp('maccstd',fieldnames(param))); param.maccstd = 1; end;
if ~any(strcmp('ftypeacc',fieldnames(param))); param.ftypeacc = arrayfun(@(x){sprintf('nsr',x)},1:NB_FOETUSES); end;
if ~any(strcmp('faccmean',fieldnames(param))); param.faccmean = repmat({0},1,NB_FOETUSES); end;
if ~any(strcmp('faccstd',fieldnames(param))); param.faccstd = repmat({1},1,NB_FOETUSES); end;
if ~any(strcmp('ftraj',fieldnames(param))); for cc=1:length(param.fhr); param.ftraj{cc} = 'none'; end; end;
if ~any(strcmp('mtraj',fieldnames(param))); param.mtraj = 'none'; end;
if ~any(strcmp('fname',fieldnames(param))); param.fname = 'aecg'; end;
if ~any(strcmp('mres',fieldnames(param))); param.mres = 0; end;
if ~any(strcmp('fres',fieldnames(param))); param.fres = zeros(1,NB_FOETUSES); end;
if ~any(strcmp('mvcg',fieldnames(param))); param.mvcg = randi([1,9]); end;
if ~any(strcmp('fvcg',fieldnames(param))); param.fvcg = randi([1,9],NB_FOETUSES,1); end;
if ~any(strcmp('evcg',fieldnames(param))); param.evcg = randi([1,4]); end;
if ~any(strcmp('posdev',fieldnames(param))); param.posdev = 1; end;
if ~any(strcmp('mectb',fieldnames(param))); param.mectb = 0; end;
if ~any(strcmp('fectb',fieldnames(param))); param.fectb = 0; end;
if isempty(debug); debug = 1; end;

% == check that parameters make sense
if size(param.elpos,2)~=3; error('run_ecg_generator: an electrode location must be defined using 3D polar coordinates \n'); end;
if param.fs>param.n; error('run_ecg_generator: the number of data point requested is inferior to the sampling frequency \n'); end;

% == constants
if any(strcmp('ntype',fieldnames(param)))
    NB_NOISES = length(param.ntype);
else
    NB_NOISES = 0;
end

%% == MATERNAL heart dipole generation
param.elpos = [param.elpos; param.refpos]; % calculating reference in same manner as other electrodes
[gp_m.norm,selvcgm] = load_gparam(param.mvcg,'normal'); % randomly pick VCG model for mother
if param.mectb % add ectopic beats?
    [gp_m.ecto,~] = load_gparam(param.evcg,'ectopic');
    % normalising the alphai to unsure same as normal beats
    VCGect = ecg_model([gp_m.ecto{2}.x gp_m.ecto{3}.x gp_m.ecto{1}.x] ,linspace(-pi,pi,250));
    VCGnorm = ecg_model([gp_m.norm{2}.x gp_m.norm{3}.x gp_m.norm{1}.x ] ,linspace(-pi,pi,250));
    gp_m.ecto{2}.x = gp_m.ecto{2}.x.*(max(abs(VCGnorm))/max(abs(VCGect)));
    VCGect = ecg_model([gp_m.ecto{2}.y gp_m.ecto{3}.y gp_m.ecto{1}.y ] ,linspace(-pi,pi,250));
    VCGnorm = ecg_model([gp_m.norm{2}.y gp_m.norm{3}.y gp_m.norm{1}.y ] ,linspace(-pi,pi,250));
    gp_m.ecto{2}.y = gp_m.ecto{2}.y.*(max(abs(VCGnorm))/max(abs(VCGect)));
    VCGect = ecg_model([gp_m.ecto{2}.z gp_m.ecto{3}.z gp_m.ecto{1}.z ] ,linspace(-pi,pi,250));
    VCGnorm = ecg_model([gp_m.norm{2}.z gp_m.norm{3}.z gp_m.norm{1}.z ] ,linspace(-pi,pi,250));
    gp_m.ecto{2}.z = gp_m.ecto{2}.z.*(max(abs(VCGnorm))/max(abs(VCGect)));
end
rm = 0.05; % radius around origin allowed for maternal heart to be
L_m = eye(3); % scaling of dipole in each direction
theta0_m = pi/3; % inital phase of the model for the MECG
vols.Rm = struct('x', 0, 'y', 0, 'z', 0); % rotation matrix

mh_cart = param.mheart;
[mh_cart(1),mh_cart(2)] = pol2cart(param.mheart(1),param.mheart(2));

if param.posdev
    % = if we do not force the location of the heart we shift it slightly
    [xp,yp,zp] = sph2cart(2*pi*rand-1,asin(2*rand-1),rm*rand); % random deviation from mat. heart origin
    mh_cart = mh_cart + [xp,yp,zp]; % final location for mat. heart (cartesian)
end

[vols.mheart{1}(1),vols.mheart{1}(2),vols.mheart{1}(3)] = cart2pol(mh_cart(1),mh_cart(2),mh_cart(3));  % new location (cyl. coord.)
param.mheart = vols.mheart{1};

% == maternal heart rate variability
strhrv.hr = param.mhr;
strhrv.lfhfr = 0.6;
strhrv.hrstd = 2;
strhrv.flo = param.mres;
strhrv.fhi = 0.25;
strhrv.acc = param.macc;
strhrv.typeacc = param.mtypeacc;
strhrv.accmean = param.maccmean;
strhrv.accstd = param.maccstd;
[theta_m,w_m] = generate_hrv(strhrv,param.n,param.fs,theta0_m);

% == electrodes position and reference electrode
vols.elpos = param.elpos;
[Xc,Yc] = pol2cart(vols.elpos(:,1),vols.elpos(:,2)); % converting from polar to cartesian coordinate system
epos = [Xc,Yc,vols.elpos(:,3)]; % electrodes position

%== adding maternal heart trajectory (e.g. position change)
if strcmp(param.mtraj,'none')
   mtraj = mh_cart;
else
   xl=linspace(0,mh_cart(1));yl=linspace(0,mh_cart(2));zl=linspace(0,mh_cart(3)); % line towards origin
   idx=randi([50,100],1,3); % picking three coordinates on second half
   mh_cart2 = [xl(idx(1)) yl(idx(2)) zl(idx(3))];
   mtraj = traject_generator(param.n,mh_cart,mh_cart2,param.mtraj); % defining a trajectory to foetal movement
end

%== generate MATERNAL heart dipole
disp('Generating maternal model...')
m_model = add_cardiacdipole(param.n,param.fs,gp_m,L_m,theta_m,w_m,param.mres,vols.Rm,epos,mtraj,0);
m_model.type = 1; % maternal ecg is type 1

%% == Foetal heart(s) generation
L_f = eye(3); % scaling of dipole in each direction
Rfh = 0.1; % radius allowed for foetal heart to appear

% = foetal dipole generation
f_model = cell(NB_FOETUSES,1); % allocating memory
w_f = f_model; theta_f = f_model; gp_f = f_model; vols.param.fheart = f_model;
selvcgf = cell(NB_FOETUSES,1);
for fet=1:NB_FOETUSES
    disp(['Generating model for fetus ' num2str(fet) ' ..'])
    
    % = foetal dipole parameters
    fh_cart = param.fheart{fet};
    [fh_cart(1),fh_cart(2)] = pol2cart(param.fheart{fet}(1),param.fheart{fet}(2)); % convert in cartesian coordinates
    
    if param.posdev
        % picking random location around origin for fetus to be
        [xp,yp,zp] = sph2cart(2*pi*rand, ...   % new random location 
                              asin(2*rand-1),...   % for fet. heart within sphere
                              Rfh*rand); 
        posf_start = [xp+fh_cart(1),yp+fh_cart(2),zp+fh_cart(3)]; % new foetal start position
    else
        posf_start = [fh_cart(1),fh_cart(2),fh_cart(3)];
    end
    
    % picking location for final position (translation)
    xl=linspace(0,posf_start(1));yl=linspace(0,posf_start(2));zl=linspace(0,posf_start(3)); % line towards origin
    idx=randi([50,100],1,3); % picking three coordinates on second half
    posf_end = [xl(idx(1)) yl(idx(2)) zl(idx(3))];
    
    [vols.fheart{fet}(1), vols.fheart{fet}(2), vols.fheart{fet}(3)] = cart2pol(posf_start(1),posf_start(2),posf_start(3));
    param.fheart{fet} = vols.fheart{fet};
    % = randomly pick VCG model for fetus (load Gaussian parameters)
    [gp_f{fet}.norm,selvcgf{fet}] = load_gparam(param.fvcg(fet),'normal');
    if param.fectb;
        [gp_f{fet}.ecto,~] = load_gparam(param.evcg,'ectopic');
        % normalising the alphai to unsure same as normal beats
        VCGect = ecg_model([gp_f{fet}.ecto{2}.x gp_f{fet}.ecto{3}.x gp_f{fet}.ecto{1}.x] ,linspace(-pi,pi,250));
        VCGnorm = ecg_model([gp_f{fet}.norm{2}.x gp_f{fet}.norm{3}.x gp_f{fet}.norm{1}.x ] ,linspace(-pi,pi,250));
        gp_f{fet}.ecto{2}.x = gp_f{fet}.ecto{2}.x.*(max(abs(VCGnorm))/max(abs(VCGect)));
        VCGect = ecg_model([gp_f{fet}.ecto{2}.y gp_f{fet}.ecto{3}.y gp_f{fet}.ecto{1}.y ] ,linspace(-pi,pi,250));
        VCGnorm = ecg_model([gp_f{fet}.norm{2}.y gp_f{fet}.norm{3}.y gp_f{fet}.norm{1}.y ] ,linspace(-pi,pi,250));
        gp_f{fet}.ecto{2}.y = gp_f{fet}.ecto{2}.y.*(max(abs(VCGnorm))/max(abs(VCGect)));
        VCGect = ecg_model([gp_f{fet}.ecto{2}.z gp_f{fet}.ecto{3}.z gp_f{fet}.ecto{1}.z ] ,linspace(-pi,pi,250));
        VCGnorm = ecg_model([gp_f{fet}.norm{2}.z gp_f{fet}.norm{3}.z gp_f{fet}.norm{1}.z ] ,linspace(-pi,pi,250));
        gp_f{fet}.ecto{2}.z = gp_f{fet}.ecto{2}.z.*(max(abs(VCGnorm))/max(abs(VCGect)));           
    end;
    % == rotation
    if param.posdev
        theta0_f = (2*rand-1)*pi; % inital phase of the model for the foetus ECG (random [-pi,pi])
        r0 = (2*rand(1,3)-[1 1 1]).*pi; % initial rotation angles
        vols.Rf{fet} = struct('x', r0(1), 'y', r0(2), 'z', r0(3));
    else
        theta0_f = -pi/2;
        vols.Rf{fet} = struct('x', -3*pi/4, 'y', 0, 'z', -pi/2);
    end
    
    %== Heart cycle parameters
    strhrv.hr = param.fhr(fet);
    strhrv.lfhf = 0.8;
    strhrv.hrstd = 3;
    strhrv.flo = param.fres(fet);
    strhrv.flhi = 0.25;
    strhrv.acc = param.facc(fet);
    strhrv.typeacc = param.ftypeacc{fet};
    strhrv.accmean = param.faccmean{fet};
    strhrv.accstd = param.faccstd{fet};
    
    [theta_f{fet},w_f{fet}] = generate_hrv(strhrv,param.n,param.fs,theta0_f);
    
    %== Translation
    traj = traject_generator(param.n,posf_start,posf_end,param.ftraj{fet}); % defining a trajectory to foetal movement
    
    % = Generating foetal dipole
    f_model{fet} = add_cardiacdipole(param.n,param.fs,gp_f{fet},L_f,theta_f{fet},w_f{fet},...
        param.fres(fet),vols.Rf{fet},epos,traj,0);
    f_model{fet}.type = 2; % foetal ecg is type 2
end

%% == NOISE DIPOLE(s)
% considering that noise sources are stationary (no rotation nor translation). Case noise
% should follow a cardiac dipole, just use the H matrix from the heart dipole for
% propagating the noise source.

n_model = cell(NB_NOISES,1);
for n = 1:NB_NOISES
    disp(['Generating model for noise source ' num2str(n) ' ..'])
    [xn,yn] = pol2cart(2*pi*rand,0.1*rand); % random location for noise
    % each even noise source is going to be around center of cylinder,
    % while odd sources are around base
    pos_noise = [xn,yn,0.1*rand-(0.5*mod(n,2))]; % inside small semi-sphere
    % generating dipole
    [n_model{n}, tmp_handle, noise_misc] = add_noisedipole(param.n,param.fs,param.ntype{n},...
        epos,pos_noise,debug);
    if ~isempty(tmp_handle)
        noise_f_handles = [noise_f_handles, tmp_handle];
    end
    n_model{n}.SNRfct = param.noise_fct{n};
    n_model{n}.pos = pos_noise;
    n_model{n}.type = 3;
end

% == Obtaining information about peak locations
fqrs = cell(NB_FOETUSES,1);
mqrs = phase2qrs(m_model.theta);
for ff=1:NB_FOETUSES
    fqrs{ff} = phase2qrs(f_model{ff}.theta);
end

% == PROPAGATION onto ELECTRODES
disp('Projecting dipoles...')
[mixture,mecg,fecg,noise, ecg_f_handle] = generate_ecg_mixture(debug,param.SNRfm,...
    param.SNRmn,mqrs,fqrs,param.fs,m_model,f_model{:},n_model{:});

% == an electrode in the back of the cylinder is used as the reference 
%    electrode (i.e. ground, location: [pi 0.5 -0.3])
ground = mixture(end,:);
mixture = mixture(1:end-1,:)-repmat(ground,size(mixture,1)-1,1);
ground = mecg(end,:);
mecg = mecg(1:end-1,:)-repmat(ground,size(mecg,1)-1,1);
if ~isempty(fecg)
    fecg = cellfun(@(x) x(1:end-1,:)-repmat(x(end,:),size(x,1)-1,1),fecg,'UniformOutput',0);
end
if ~isempty(noise)
    noise = cellfun(@(x) x(1:end-1,:)-repmat(x(end,:),size(x,1)-1,1),noise,'UniformOutput',0);
end
vols.refpos = param.refpos;
vols.elpos = vols.elpos(1:end-1,:); % removing ground electrode

%% == FORMATING OUTPUT ARGUMENTS
out.mixture = mixture;
out.mecg = mecg;
out.fecg = fecg;
out.noise = noise;
out.m_model = m_model;
out.f_model = f_model;
out.vols = vols;
out.mqrs = mqrs;
out.fqrs = fqrs;
out.param = param;
out.selvcgm = selvcgm;
out.selvcgf = selvcgf;
out.param.elpos = out.param.elpos(1:end-1,:);

noise_misc.n_model = n_model;
out.noise_misc = noise_misc;

% == Plotting results
if debug
    % the gui needs the figure handles
    f_handles = debug_plots(out,debug); 
end

out.f_handles = [f_handles, noise_f_handles, ecg_f_handle];

end


%% This function generates plots for the fecgsyn model
function f_handles = debug_plots(out,debug)
f_handles = [];

% == debug
NB_FOETUSES = size(out.f_model,1); % number of foetuses figured out from the number of foetal hea
disp(['Selected VCG for the mother: ' out.selvcgm]);
for vv=1:NB_FOETUSES
    disp(['Selected VCG for the foetus: ' out.selvcgf{vv}]);
end

disp('Generating plots ..')
% == plots a few final AECG channels
tmp_handle = figure('name','Some generated AECG');
if debug == 11 % corresponds to running the code from the gui
    set(tmp_handle, 'Visible', 'off');
end
f_handles = [f_handles, tmp_handle];

col = [1,0,0; % red
    0,0,1; % blue
    0,0.8,0; % green
    0.4,0.4,0; % dark yellow
    0,0.8,0.8; % cyan
    0.4,0,0.8; % dark magenta
    0.8,0.4,1; % light magenta
    0.4,0.4,1]; % lilac

LINE_WIDTH = 3;
FONT_SIZE = 15;
NB_EL2PLOT = min(3,size(out.mixture,1)); % number of electodes to plot
PACE = 1;
if NB_EL2PLOT<=NB_EL2PLOT*PACE
    compt = 0;
    tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
    ax = zeros(NB_EL2PLOT,1);
    for ee=1:PACE:PACE*NB_EL2PLOT
        compt = compt+1;
        ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mixture(ee,:),...
            'color',col(compt,:),'LineWidth',LINE_WIDTH);
        xlabel('Time [sec]'); ylabel('Amplitude [NU]');
        set(gca,'FontSize',FONT_SIZE);
    end
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    linkaxes(ax,'x');
else
    warning('not enough input channel for plotting with default configuration \n');
end

if debug>1
    if debug ~= 11   % only close all figures if not running in gui mode
        close all;
    end
    % == plot mother and foetuse(s) VCGs
    tmp_handle = figure('name','VCG plots');
    if debug == 11 % corresponds to running the code from the gui
        set(tmp_handle, 'Visible', 'off');
    end
    f_handles = [f_handles, tmp_handle];
    
    LegCell = cell(NB_FOETUSES*2,1);
    for vv=1:3
        % == plot
        ax(2*vv-1)=subplot(3,2,2*vv-1); plot(tm,out.m_model.VCG(vv,:),'color',col(2*vv-1,:),'LineWidth',LINE_WIDTH);
        hold on, plot(tm(out.mqrs),out.m_model.VCG(vv,out.mqrs),'+k','LineWidth',LINE_WIDTH);
        legend(['mother VCG channel: ' int2str(vv)],'MQRS');
        set(gca,'FontSize',FONT_SIZE);
        xlabel('Time [sec]'); ylabel('Amplitude [NU]');
        for fet=1:NB_FOETUSES
            ax(2*vv)=subplot(3,2,2*vv); plot(tm,out.f_model{fet}.VCG(vv,:),'color',col(2*vv+fet,:),'LineWidth',LINE_WIDTH);
            hold on, plot(tm(out.fqrs{fet}),out.f_model{fet}.VCG(vv,out.fqrs{fet}),'+k','LineWidth',LINE_WIDTH);
            LegCell(2*(fet-1)+1) = {['foetus ' int2str(fet) ' VCG channel: ' int2str(vv)]};
            LegCell(2*fet) = {['FQRS ' int2str(fet)]};
        end
        legend(LegCell);
        set(gca,'FontSize',FONT_SIZE);
        xlabel('Time [sec]'); ylabel('Amplitude [NU]');
    end
    linkaxes(ax,'x'); xlim([0 tm(end)]);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end

if debug>2
    % == plot the projection of mother and foetuse(s) VCGs
    tmp_handle = figure('name','Projected FECG and MECG before being mixed');
    if debug == 11 % corresponds to running the code from the gui
        set(tmp_handle, 'Visible', 'off');
    end
    f_handles = [f_handles, tmp_handle];
    
    LegCell = cell(NB_FOETUSES+1,1);
    GAIN_F = 1;
    if NB_EL2PLOT<=NB_EL2PLOT*PACE
        compt = 0;
        tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
        ax = zeros(NB_EL2PLOT,1);
        for ee=1:PACE:PACE*NB_EL2PLOT
            compt = compt+1;
            ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mecg(ee,:),...
                'color','b','LineWidth',LINE_WIDTH);
            LegCell(1) = {'MECG'};
            for fet=1:NB_FOETUSES
                hold on, plot(tm,GAIN_F*out.fecg{fet}(ee,:),'color',col(fet+3,:),'LineWidth',LINE_WIDTH);
                LegCell(1+fet) = {['FECG ' int2str(fet) ', Gain: ' int2str(GAIN_F)]};
            end
            xlabel('Time [sec]'); ylabel('Amplitude [NU]');
            set(gca,'FontSize',FONT_SIZE);
            legend(LegCell);
        end
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
        linkaxes(ax,'x');
    else
        warning('not enough input channel for plotting with default configuration \n');
    end
end

if debug>3
    % == plot the volume conductor
    tmp_handle = figure('name','Volume conductor');
    if debug == 11 % corresponds to running the code from the gui
        set(tmp_handle, 'Visible', 'off');
    end
    f_handles = [f_handles, tmp_handle];
    
    plot3_volume(out.vols);
    hold on
    for i=1:NB_FOETUSES % plotting each foetuses trajectory
        plot3(out.f_model{i}.traj(:,1),out.f_model{i}.traj(:,2),out.f_model{i}.traj(:,3),'--g','LineWidth',2);
    end
    axis square
end

if debug>4
    % == plot mother and foetus heart rates
    tmp_handle = figure('name','Heart rate');
    if debug == 11 % corresponds to running the code from the gui
        set(tmp_handle, 'Visible', 'off');
    end
    f_handles = [f_handles, tmp_handle];
    
    fhr = cell(NB_FOETUSES,1);
    legstr = cell(NB_FOETUSES+1,1);
    mhr = 60./diff(out.mqrs/1000);
    plot(tm(out.mqrs(1:end-1)),mhr,'color','r','LineWidth',LINE_WIDTH,'LineStyle','--');
    legstr{1} = 'MQRS';
    for ff=1:NB_FOETUSES
        fhr{ff} = 60./diff(out.fqrs{ff}/1000);
        hold on, plot(tm(out.fqrs{ff}(1:end-1)),fhr{ff},'color',col(ff+1,:),'LineWidth',LINE_WIDTH,'LineStyle','-');
        legstr{ff+1} = ['FQRS' int2str(ff)];
    end
    legend(legstr);
    xlabel('Time [sec]'); ylabel('FHR [bpm]')
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    legend boxoff;
end

end
