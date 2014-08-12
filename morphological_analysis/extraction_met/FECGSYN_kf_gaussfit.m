function [ker,model,abserr,status] = FECGSYN_kf_gaussfit(tecg,maxIter,nbgauss,debug)
% this function fits seven Gaussians to a template ecg and return the optimum
% Gaussian parameters [alphai bi tetai]. If the fitting residual e is too high then the fiducial
% location are slightly shifted (randomly) and another fitting is
% performed. This operation is done iteratively until the Gaussian are well
% fitted or that the maximum number of iterations has been reached.
%
% inputs
%   tecg:       Template ECG
%   maxIter:    Max number of iterations
%   nbgauss:    Number of Gaussians
%   debug:      Enter debug mode? 0/1
%
% outputs
%   ker:        Gaussian ker [alphai bi tetai]
%   model:      Mapping of tecg using the Gaussians
%   e:          Measure of how well the fitting was performed (RMS)
%   status:     Success/failure (failure if failed to each e < TOL)
% 
%
% Dual EKF Toolbox, version 1.0, March 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 08-04-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == manage inputs
if nargin<3; disp('gauss_fit: wrong number of input arguments \n'); return; end;
if nargin<4; debug=0; end;

% == constants
NB_BINS = length(tecg);
MAX_NB_ITER = maxIter;
NB_GAUSS = nbgauss; % number of Gaussians to use
phase = linspace(-pi,pi,NB_BINS);
status = 1;
TOL = 5; % 0.5 percent tolerance

% == core function
try
    [ker,model] = FECGSYN_kf_findgaussparam(phase,tecg,NB_GAUSS);
    eo = 100*mean((model-tecg).^2)/mean(tecg.^2); % e old
    compt = 0;
    %udt = 0;

    while (eo>TOL && compt<=MAX_NB_ITER) %|| udt==0 
        % while we have not reached a satisfactory e try to optimise
        [kerTmp,modTemp] = FECGSYN_kf_findgaussparam(phase,tecg,NB_GAUSS); % find Gaussian param that best fits
        en = 100*mean((modTemp-tecg).^2)/mean(tecg.^2); % compute e associated with approximation-> e new

        bi = kerTmp(NB_GAUSS+1:2*NB_GAUSS);
        
        if en<eo && (isempty(find(abs(bi)<0.03,1)))
            % if improvement then store the corresponding results
            model = modTemp;
            ker = kerTmp;
            eo = en;
            %udt = 1;
        end

        compt = compt+1;
    end
    
    e = eo;
    if e<TOL; status=0; end;
    fprintf('Number of iterations for Gaussian fitting was %f \n',compt);
    abserr = model-tecg;

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(end)); end;
    status=0;
end

% == plots
if debug || ~nargout
    col = [1,0,0; % red
    0,0,1; % blue
    0,0.8,0; % green
    0.4,0.4,0; % dark yellow
    0,0.8,0.8; % cyan
    0.4,0,0.8; % dark magenta
    0.8,0.4,1; % light magenta
    0.4,0.4,1]; % lilac
    
    hold on, plot(tecg,'--r','LineWidth',2);
    fprintf('Init e for data:  %f \n',e);    
    hold on, plot(model,'k','LineWidth',2);
    xlabel('bin number'); ylabel('amplitude (mV)');
    xlim([0 250]);
    
    alphai = ker(1:NB_GAUSS);
    bi = ker(NB_GAUSS+1:2*NB_GAUSS);
    tetai = ker(2*NB_GAUSS+1:3*NB_GAUSS);
    for jj=1:NB_GAUSS
        y = alphai(jj)*gaussmf(phase,[bi(jj);tetai(jj)]);
        hold on, plot(y,'color',col(jj,:));
    end
end
end


function [goptparam,model] = FECGSYN_kf_findgaussparam(phase,tecg,nbgauss,debug)
% find Gaussians optimal paramters given an estimate of the Gaussian
% centres. The non linear fitting function is constrained to force
% associating two Gaussians with the P-wave, three with the QRS complex and
% two with the T-wave.
%
% inputs
%   phase:  tecg phase
%   tecg:   template ecg
%   debug:  enter debug mode?
% 
% outputs 
%   goptparam: Gaussians optim parameters [alphaix7 bix7 tetaix7]
%   model:     model ecg cycle
% 
% Dual EKF, version 1.0, March 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 24-03-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == manage inputs
if nargin<3; error('findgaussparam: wrong number of input arguments \n'); end;
if nargin<4; debug=0; end;

% == general
NB_GAUSS = nbgauss;
I = randi([1,round(length(tecg))],1,NB_GAUSS);
P = length(I);


% == core function
try
    % initialize tetai, alphai and bi
    tetai = phase(I(1:P));
    alphai = 1.2*tecg(I(1:P));
    bi = .04*ones(size(alphai));
    
    % param for nonlinear regression
    options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'Display','off');
    InitParams = [alphai bi tetai];
    % constraints on P1s,P2s,Qs,Rs,Ss,T1s,T2s for the theta parameter
    if NB_GAUSS==7
        lb = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf...
            -pi/1.2 -pi/1.2 -pi/3 -pi/3 -pi/2 pi/5 pi/5]; % lower bound
        ub = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf ...
            -pi/10 -pi/10 pi/3 pi/3 pi/2 pi/1.2 pi/1.2]; % upper bound
    else
        lb = [-inf(1,2*NB_GAUSS), -pi*ones(1,NB_GAUSS)];
        ub = [inf(1,2*NB_GAUSS), pi*ones(1,NB_GAUSS)];
    end
    
    % perform nonlinear regression
    OptParams = lsqcurvefit(@ecg_model,InitParams,phase,tecg,lb,ub,options);
    
    % == remove problematic Gaussians
    % Look if there is any invalid Gaussian and remove it
    % If gaussian width or height is less than 0.001 it will excluded.
    N = NB_GAUSS;
    yy = 1;
    while(yy<N+1)
        if((abs(OptParams(N+yy))<0.001)||(abs(OptParams(yy))<0.001))
            OptParams(2*N+yy)=[];OptParams(N+yy)=[];OptParams(yy)=[];
            N = N-1;
        else
            yy=yy+1;
        end
    end
    
    % == generating output
    model = ecg_model(OptParams,phase);
    goptparam = OptParams;

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
end

if debug
   plot(tecg,'LineWidth',2); hold on, plot(model,'r','LineWidth',2);
   xlabel('Samp Nb'); ylabel('Amp [NU]'); legend('tecg','Gaussian reconstructed tecg');
end

end






