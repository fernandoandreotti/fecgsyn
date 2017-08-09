function traj =  traject_generator(N,pos_i,pos_f,type)
% function traj =  traject_generator(N,pos_i,pos_f,type)
% this function is used to generate trajectories within volume conductor. These
% trajectories may be applied to heart dipole in order to generate a more realistic
% modelling, e.g. respiration or fetal movements.
%
% inputs:
%   N      size of the noise to generate at fs (sampling frequency) [datapoint number]
%   pos_i  initial position for trajectory [1x3 number array]
%   pos_f  final position for trajectory [1x3 number array]
%   type   type of trajectory to be build, e.g. 'none','linear', 'spline' or 'helix' [string]% fecgsyn toolbox, version 1.2, Jan 2017
%
% output:
%   traj   generated trajectory (Nx3 number matrix)
%
%
%
% Examples:
% TODO
%
% See also:
% run_ecg_generator
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

% == check inputs
if nargin<4; error('traject_generator: not enough input arguments'); end;
    
% == constants
NB_CIRC = 3.5; % number of cycles for helix trajectory

% == main
traj = zeros(N,3);
switch type
    case 'none'
        traj = pos_i;
    case 'step' % step change in trajectory
        traj = repmat(pos_i,N,1);
        traj(round(N/2)+1:end,:) = repmat(pos_f,round(N/2),1);
    case 'linear' % linear trajectory
        trajx = linspace(pos_i(1),pos_f(1),N)'; % linear trajectory for X
        trajy = linspace(pos_i(2),pos_f(2),N)'; % linear trajectory for Y
        trajz = linspace(pos_i(3),pos_f(3),N)'; % linear trajectory for Z
        traj = [trajx trajy trajz];
    case 'spline' % spline trajectory for one coordinate
        idx = randperm(3);
        pos_mid = (pos_f-pos_i)./2; % inserting middle point for spline
        for i = idx(1:2)
           traj(:,i) = linspace(pos_i(i),pos_f(i),N)';            
        end
        traj(:,idx(3))=spline([pos_i(idx(2)) pos_mid(idx(2)) pos_f(idx(2))],[pos_i(idx(3)) 0.5*pos_mid(idx(3)) pos_f(idx(3))],traj(:,idx(2))); % spline trajectory for third coordinate
    case 'helix' % circular trajectory in xy
        traj(:,3) = linspace(pos_i(3),pos_f(3),N)'; % linear trajectory for Z
        center = (pos_i(1:2)+pos_f(1:2))/2;
        r = sqrt(sum((pos_i(1:2)-center).^2)); % radius half the distance between pos_i and pos_f
        w = linspace(0,2*pi*NB_CIRC,N)'; % linear in Z
        phi = atan((pos_i(2)-center(2))/(pos_i(1)-center(1)));
        traj(:,1) = r*cos(w+phi)+center(1); % circle in X Y
        traj(:,2) = r*sin(w+phi)+center(2);
    otherwise
        error('TrajectGenerator: Unknown trajectory type')        
end
    