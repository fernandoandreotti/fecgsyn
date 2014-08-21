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

function f_handle = FECGSYN_UI_create_electrode_plot(param, out)
% CREATE_ELECTRODE_PLOT Generates the plot for visualizing the electrodes

% give priority to the out structure, if it is given. otherwise use the
% param structure.
%
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

if nargin < 2
    NB_FOETUSES = size(param.fheart,1); % number of foetuses figured out from the number of foetal hearts
else
    NB_FOETUSES = size(out.f_model,1); % number of foetuses figured out from the number of foetal hea
    param = out.param;
end



% prep: copied out of run_ecg_generator() so that a set of trial values for
% electrodes can be plotted

rm = 0.01; % radius around origin allowed for maternal heart to be
% L_m = eye(3); % scaling of dipole in each direction
% teta0_m = pi/3; % inital phase of the model for the MECG
vols.Rm = struct('x', 0, 'y', 0, 'z', 0); % rotation matrix

mh_cart = param.mheart;
[mh_cart(1),mh_cart(2)] = pol2cart(param.mheart(1),param.mheart(2));

if param.posdev
    % = if we do not force the location of the heart we shift it slightly
    [xp,yp,zp] = sph2cart(2*pi*rand-1,asin(2*rand-1),rm*rand); % random deviation from mat. heart origin
    mh_cart = mh_cart + [xp,yp,zp]; % final location for mat. heart (cartesian)
end

[vols.mheart{1}(1),vols.mheart{1}(2),vols.mheart{1}(3)] = cart2pol(mh_cart(1),mh_cart(2),mh_cart(3));  % new location (cyl. coord.)

% electrodes position and reference electrode
vols.elpos = param.elpos;
[Xc,Yc] = pol2cart(vols.elpos(:,1),vols.elpos(:,2)); % converting from polar to cartesian coordinate system
epos = [Xc,Yc,vols.elpos(:,3)]; % electrodes position

f_model = cell(NB_FOETUSES,1); % allocating memory
% w_f = f_model; teta_f = f_model; gp_f = f_model; 
vols.param.fheart = f_model;
% selvcgf = cell(NB_FOETUSES,1);
for fet=1:NB_FOETUSES
%     disp(['Generating model for fetus ' num2str(fet) ' ..'])
    
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
    idx=randi([50,100],1,3);     % picking three coordinates on second half
%     posf_end = [xl(idx(1)) yl(idx(2)) zl(idx(3))];
    
    [vols.fheart{fet}(1), vols.fheart{fet}(2), vols.fheart{fet}(3)] = cart2pol(posf_start(1),posf_start(2),posf_start(3));

    % == rotation
    if param.posdev
%         teta0_f = (2*rand-1)*pi; % inital phase of the model for the foetus ECG (random [-pi,pi])
        r0 = (2*rand(1,3)-[1 1 1]).*pi; % initial rotation angles
        vols.Rf{fet} = struct('x', r0(1), 'y', r0(2), 'z', r0(3));
    else
%         teta0_f = -pi/2;
        vols.Rf{fet} = struct('x', -3*pi/4, 'y', 0, 'z', -pi/2);
    end

end

vols.refpos = param.refpos;
% vols.elpos = vols.elpos(1:end-1,:); % removing ground electrode


% == plot the volume conductor
f_handle = figure('name','Volume conductor');
set(f_handle, 'Visible', 'off');

plot3_volume(vols);
axis square