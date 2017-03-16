function dmodel = add_cardiacdipole(N,fs,gp_all,L,...
    theta,w,fres,R0,epos,traj,debug)
% function dmodel = add_cardiacdipole(N,fs,gp_all,L,theta,w,fres,R0,epos,traj,debug)
% synthetic cardiac dipole generator using the 'direct form' of the cardiac
% dipole equation. This function generates the vectocardiogram (VCG) of 
% the mother or foetus as well as the Dower-like matrix that allows the 
% projection of the VCG onto the electrode locations as specified by (elec.ep). 
% In this function the rotation of the heart axis and movement of the heart 
% are modelled. Rotation of the heart is modelled by dynamically updating 
% the rotation matrix R with respect to breathing. Movement of the foetus 
% heart is modelled by updating the dower-like projection matrix H. These
% components of the model intend to reflect the non-stationary aspect of 
% the FECG signal. The function returns the VCG signal, the projection
% matrix H (which is 2D in the case there is foetal movement and 3D in the
% case the foetus is moving).
% Input:
%        N:     signal length [number of points]
%       fs:     sampling rate [Hz]
%       gp_all: Gaussian parameters [cell of cells]
%           gp{1}{:}: Gaussian parameters of mother/foetus ecg - normal beats
%               gp{i}{1}:  structure contaning the phase of Gaussian functions used for
%                   modeling the x, y, and z coordinates of the cardiac
%                   dipole
%               gp{i}{2}:  structure contaning the amplitudes of Gaussian functions used for
%                   modeling the x, y, and z coordinates of the cardiac dipole
%               gp{i}{3}:  structure contaning the widths of Gaussian functions used for
%                   modeling the x, y, and z coordinates of the cardiac dipole
%           gp{2}{:}: Gaussian parameters of moether/foetus ecg - ectopic beat
%               
%       L:      scaling of dipole in each direction [3x3 matrix]
%     theta:     phase for heart dipole model
%        w:     angular frequency
%     fres:     respiration frequency (for heart dipole rotation) [Hz]
%       R0:     initial rotation angles
%     epos:     position of electrodes [normalised units]
%     traj:     [Nx3 matrix] representing trajectory taken by the dipole, if it has some 
%               trainslation other than respiration. If no translation is given, traj    
%               is the initial position of dipole [1x3 matrix]. Columns are trajectory in
%               x,y and z direction.
%    debug:     debug [bool]
%
% Output:
%     dmodel   structure contaning the dipole model i.e.:
%        dmodel.H:      Dower-like matrix for dipole either 2D (time invariant) or 3D 
%                           (variant case).
%        dmodel.VCG:    VCG for dipole
%        dmodel.type:   maternal (1) or foetal (2) dipole
%        dmodel.traj:   trajectory taken by dipole
%        dmodel.stm: state transition matrix
%        dmodel.rax: max respiration angle in radian around X 
%        dmodel.ray: max respiration angle in radian around Y
%        dmodel.raz: max respiration angle in radian around Z
%        dmodel.rht: volumes height allowed for heart translation due to respiration
%
%
% References
% [1] Leanderson et al. "Estimation of the respiratory frequency using spatial information in
% the VCG". Medical Engineering & Physics 25 (2003) 501-507. 
% this paper gives an idea of the Rx, Ry, Rz angles variation with breathing for the VCG.
%
% [2] Sameni, Reza, et al. "Multichannel ECG and noise modeling: application to 
% maternal and foetal ECG signals." EURASIP Journal on Advances in Signal Processing 
% 2007 (2007).
%
% [3] Oster, Julien, and Gari D. Clifford. "An Artificial Model of the 
% Electrocardiogram during Paroxysmal Atrial Fibrillation." Computing in
% cardiology 2013.
%
%
% Examples:
% TODO
%
% See also:
% run_ecg_generator
% add_noisedipole
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

% == checking inputs
if nargin<10; error('add_cardiacdipole: not enough input argments \n'); end;
if any(strcmp('ecto',fieldnames(gp_all)));
    ect = 1; % we are including ectopic beats into the model
    gp{2} = gp_all.ecto;
    rn = 0.7 + 0.1*rand; % a different transition matrix is generated at each iteration
    re = 0.2 + 0.1*rand;
    STM = [rn (1-rn); (1-re) re]; % state transition matrix HMM model
    S = [1;0]; % initial state (normal beat)
else
    ect = 0;
end;
gp{1} = gp_all.norm;

% trajectory check
if any(sqrt(traj(:,1).^2+traj(:,2).^2) > 0.5) || any(abs(traj(:,3))>0.5) % test if trajectories make sense
    error('add_cardiacdipole: Trajectory extrapolates the volume conductor. Limiting trajectory to volume conductor');
end

% == constants
RESP_ANG_X = 2*0.1; % max respiration angle in radian around X
RESP_ANG_Y = 2*0.08; % max respiration angle in radian around Y
RESP_ANG_Z = 2*0.07; % max respiration angle in radian around Z
HEART_T_RESP = 0.05; % allowing 5% of volumes height for heart translation due to respiration
NB_EL = size(epos,1); % number of electrodes

dt = 1/fs; % time pace
VCG = zeros(3,N);
ncy = find(diff(theta)<0); % new cycle's locations

% == generates breathing waveform for rotation matrix modulation
if fres==0
    brwave = zeros(N,1)';
else
    brwave = generate_breathing(fs,N,fres,0);
end
% == Computing VCG and Dower H matrix
if (fres==0 && size(traj,1)==1) % respiration (FALSE), trajectory (FALSE)
    den_norm = diag(1./sqrt(sum((epos-repmat(traj,NB_EL,1)).^2,2)).^3); % denominator's norm (from H equation)
    h_1 = den_norm*(epos-repmat(traj,NB_EL,1));
    H = h_1; % H for invariant case
    % trajectory is only a vector
elseif  (fres~=0 && size(traj,1)==1) % respiration (TRUE), trajectory (FALSE)
    H = zeros(NB_EL,3,N);
    traj = repmat(traj,N,1);
else % respiration (FALSE) AND trajectory (TRUE) OR respiration (TRUE) AND trajectory (TRUE)
    H = zeros(NB_EL,3,N);    
end

% == heart translation with respiration
if fres==0
    dpos = traj; % no translation of the heart
else    
    dpos = traj + [zeros(N,2) HEART_T_RESP*(brwave+0.5)']; % dipole absolute position for each time
end

crst = 1; X = 0; Y = 0; Z = 0;
% == loop every sample
for i=1:N
    if ect && sum(i==ncy)
       % = if considering ectopic beats
       rd_nb = rand;
       crst = find(S==1); % current state
       if rd_nb>STM(crst,1) % state 1 (normal) or 2 (ectopic)
            S = [0 1];
       else
            S = [1 0];
       end 
    end
    
    % = guarantees that dtheta stay in [-pi,pi]
    dthetaix = mod(theta(ones(length(gp{crst}{1}.x),1),i)' - gp{crst}{1}.x + pi , 2*pi) - pi;
    dthetaiy = mod(theta(ones(length(gp{crst}{1}.y),1),i)' - gp{crst}{1}.y + pi , 2*pi) - pi;
    dthetaiz = mod(theta(ones(length(gp{crst}{1}.z),1),i)' - gp{crst}{1}.z + pi , 2*pi) - pi;  

    % = differential expression
    X = X - dt*sum(w(i)*gp{crst}{2}.x ./ (gp{crst}{3}.x .^ 2) .* dthetaix .* exp(-dthetaix .^2 ./ (2*gp{crst}{3}.x .^ 2)),2);
    Y = Y - dt*sum(w(i)*gp{crst}{2}.y ./ (gp{crst}{3}.y .^ 2) .* dthetaiy .* exp(-dthetaiy .^2 ./ (2*gp{crst}{3}.y .^ 2)),2);
    Z = Z - dt*sum(w(i)*gp{crst}{2}.z ./ (gp{crst}{3}.z .^ 2) .* dthetaiz .* exp(-dthetaiz .^2 ./ (2*gp{crst}{3}.z .^ 2)),2);
    
    % = analytical expression
    % X = sum(gp{crst}{2}.x .* exp(-dthetaix .^2 ./ (2*gp{crst}{3}.x .^ 2)),2);
    % Y = sum(gp{crst}{2}.y .* exp(-dthetaiy .^2 ./ (2*gp{crst}{3}.y .^ 2)),2);
    % Z = sum(gp{crst}{2}.z .* exp(-dthetaiz .^2 ./ (2*gp{crst}{3}.z .^ 2)),2);
    
    % rotation due to respiration
    thetax = R0.x + RESP_ANG_X*brwave(i);
    thetay = R0.y + RESP_ANG_Y*brwave(i);
    thetaz = R0.z + RESP_ANG_Z*brwave(i);
    
    R = rotatexyz(thetax,thetay,thetaz); % rotation matrix
    
    VCG(:,i) = R*L*[X; Y; Z]; % VCG with rotation
    
    % generating the projection matrix (variant case)
    if size(traj,1) > 1 % which occurs case respiration or translation is added
        dr = repmat(dpos(i,:),NB_EL,1); % changing position of dipole     
        den_norm = 1./sqrt(sum((epos-dr).^2,2)).^3; % denominator's norm (from H equation)
        h_1 = diag(den_norm)*(epos-dr);
        H(:,:,i) = h_1; 
    end
end

% == format outputs
dmodel.H = H;
[B,A] = butter(5,.7*2/fs,'low'); % high-pass VCG filter with .5 Hz
opol = 6;  % polynome order

% == avoid trends in the VCG (due to non-zero differentials  
for i = 1:3                      
    % gross trends
    [p,~,mu] = polyfit(1:N,VCG(i,:),opol);   % polynomial fit
    f_y = polyval(p,1:N,[],mu);       
    VCG(i,:) = VCG(i,:) - f_y; % remove rough trends
    % fine trends
    VCG(i,:) = VCG(i,:) - filtfilt(B,A,VCG(i,:));  % at the end/beggining of cycle)
    % normalizing VCG for propagating
    VCG(i,:) = VCG(i,:)./max(abs(VCG(i,:)));   
end
dmodel.VCG = VCG;
dmodel.theta = theta;
dmodel.traj = traj;
if ect; dmodel.stm = STM; end;
dmodel.rax = RESP_ANG_X;
dmodel.ray = RESP_ANG_Y;
dmodel.raz = RESP_ANG_Z;
dmodel.rht = HEART_T_RESP;

% == debug
if debug
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
   tm = 1/fs:1/fs:N/fs;
   ax = zeros(3,1);
   for vv=1:3
        ax(vv)=subplot(3,1,vv); plot(tm,VCG(vv,:),'color',col(vv,:),'LineWidth',LINE_WIDTH);
        set(gca,'FontSize',FONT_SIZE);
   end
   xlabel('Time [sec]'); ylabel('Amplitude [NU]');   
   set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
   linkaxes(ax,'x');
end

end






