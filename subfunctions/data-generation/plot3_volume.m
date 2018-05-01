function plot3_volume(vols)
% function plot3_volume(vols)
% plot volume conductor together with foetuses and mother heart positions 
% as well as electrodes positions. This allows for visual representation 
% of simulation configuration. The volume conductor is represented by a
% cylinder centred on zero and normalised alond each direction (in [-0.5
% 0.5]). Remember that the hearts are represented as dipoles which position 
% in space is determined by vols.fheart and vols.mheart with a certain
% orientation specified by the rotation matrices (vols.Rf and vols.Rm).
%
% inputs:
%   vols: volume structure containing the different objects to position.
%       vols.fheart:   foetal heart position in polar coordinates, 
%                      (cell of 3x1 vectors, one cell per foetus thus 
%                      supporting multiple pregnancies representation)
%       vols.mheart:   mother heart position in polar coordinates, 
%                      (one cell of 3x1 vector)
%       vols.elpos:    electrodes position in polar coordinates
%                      (3 x NB_ELECTRODES vector)
%       vols.Rf:       rotation parameters for the foetus heart (thetaX,
%                      thetaY,thetaZ - in radian)
%       vols.Rm:       rotation parameters for the mother heart (thetaX,
%                      thetaY,thetaZ - in radian)
%
%
%
% Examples:
% TODO
%
% See also:
% TODO
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
if nargin<1; error('plot3_volume: wrong number of input argument \n'); end
if isempty(vols.mheart); error('plot3_volume: vols.mheart is empty \n'); end
if isempty(vols.elpos); error('plot3_volume: vols.elpos is empty \n'); end
if isempty(vols.Rm); error('plot3_volume: vols.Rm is empty \n'); end

% == general
HEART_SIZE = 0.07; % normalised heart size
FONT_SIZE  = 20;
if exist('vols.fheart','var')
    NB_FOETUSES = length(vols.fheart);
else
    NB_FOETUSES = 0;
end
NB_ELECTRODES = size(vols.elpos,1);
CYLINDER_RADIUS = 0.5;
BASE_ANGLE = 90;
ARROWHEAD_LENGTH = 15;

% == convert polar to cartesian coordinate system
for ff=1:NB_FOETUSES
    [xf,yf] = pol2cart(vols.fheart{ff}(:,1),vols.fheart{ff}(:,2));
    vols.fheart{ff}(:,1) = xf;
    vols.fheart{ff}(:,2) = yf;
end

[xm,ym] = pol2cart(vols.mheart{1}(:,1),vols.mheart{1}(:,2));
vols.mheart{1}(:,1) = xm;
vols.mheart{1}(:,2) = ym;    
    
% == plot conductor volume (Xc,Yc,Zc)
[Xc,Yc,Zc] = cylinder(CYLINDER_RADIUS,30);
surf(Xc,Yc,Zc-CYLINDER_RADIUS,'FaceColor','none');

% == plot foetal and mother hearts: (Xmh,Ymh,Zmh) and (Xfh,Yfh,Zfh)
hold on,
[Xmh,Ymh,Zmh] = sphere(20);
surf(vols.mheart{1}(1)+HEART_SIZE*Xmh,vols.mheart{1}(2)+...
    HEART_SIZE*Ymh,vols.mheart{1}(3)+HEART_SIZE*Zmh);

for ff=1:NB_FOETUSES
    hold on,
    [Xfh,Yfh,Zfh] = sphere(20);
    surf(vols.fheart{ff}(1)+0.5*HEART_SIZE*Xfh,vols.fheart{ff}(2)+...
        0.5*HEART_SIZE*Yfh,vols.fheart{ff}(3)+0.5*HEART_SIZE*Zfh); % foetus heart size = 0.5 the one of the mother for display purposes
end

% == plot axis of the heart
% = for the mother
hold on,
R = rotatexyz(vols.Rm.x,vols.Rm.y,vols.Rm.z);
M = [vols.mheart{1};vols.mheart{1};vols.mheart{1}] + R*0.3*eye(3,3);
MX = [vols.mheart{1};M(1,:)];
MY = [vols.mheart{1};M(2,:)];
MZ = [vols.mheart{1};M(3,:)];
if isempty(regexp(version('-release'),'2014'))
    arrow(MX(1,:),MX(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','r'); text(MX(2,1),MX(2,2),MX(2,3)+0.05,'x','FontSize',FONT_SIZE+5);
    arrow(MY(1,:),MY(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','r'); text(MY(2,1),MY(2,2),MY(2,3)+0.05,'y','FontSize',FONT_SIZE+5);
    arrow(MZ(1,:),MZ(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','r'); text(MZ(2,1),MZ(2,2),MZ(2,3)+0.05,'z','FontSize',FONT_SIZE+5);
else
    disp('Skipping drawing arrows on heart dipoles, only compatible with Matlab 2014. See plot3_volume.m')
end

% = for the foetuses
for ff=1:NB_FOETUSES
    R = rotatexyz(vols.Rf{ff}.x,vols.Rf{ff}.y,vols.Rf{ff}.z);
    %F = [vols.fheart{ff};vols.fheart{ff};vols.fheart{ff}] + R*0.3*eye(3,3);
    F = [vols.fheart{ff};vols.fheart{ff};vols.fheart{ff}] + R*0.3*eye(3,3);
    FX = [vols.fheart{ff};F(1,:)];
    FY = [vols.fheart{ff};F(2,:)];
    FZ = [vols.fheart{ff};F(3,:)];
    if isempty(regexp(version('-release'),'2014'))
        arrow(FX(1,:),FX(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','b'); text(FX(2,1),FX(2,2),FX(2,3)-0.05,'x','FontSize',FONT_SIZE+5);
        arrow(FY(1,:),FY(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','b'); text(FY(2,1),FY(2,2),FY(2,3)-0.05,'y','FontSize',FONT_SIZE+5);
        arrow(FZ(1,:),FZ(2,:), ARROWHEAD_LENGTH, 'BaseAngle',BASE_ANGLE,'Width',2,'FaceColor','b'); text(FZ(2,1),FZ(2,2),FZ(2,3)-0.05,'z','FontSize',FONT_SIZE+5);
    end
end

% == plot electrodes position (Xe,Ye,Ze)
[Xe, Ye] = pol2cart(vols.elpos(:,1),vols.elpos(:,2));
Ze = vols.elpos(:,3);

hold on,
plot3(Xe,Ye,Ze,'s','MarkerSize',30,'LineWidth',4)
for i=1:NB_ELECTRODES
    text(Xe(i),Ye(i),Ze(i),int2str(i),'Color','r','FontSize',20,'FontWeight','bold');
end

% == plot reference electrode
[Xe, Ye] = pol2cart(vols.refpos(:,1),vols.refpos(:,2));
Ze = vols.refpos(:,3);
plot3(Xe,Ye,Ze,'s','MarkerSize',30,'LineWidth',4,'Color','g')
text(Xe,Ye,Ze,'GD','Color','g','FontSize',20,'FontWeight','bold');

% == set display properties
xlabel('X','fontsize',FONT_SIZE); 
ylabel('Y','fontsize',FONT_SIZE); 
zlabel('Z','fontsize',FONT_SIZE);
set(gca,'FontSize',FONT_SIZE);
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'fontsize',FONT_SIZE);
xlim([-0.6 0.6]);
ylim([-0.6 0.6]);
zlim([-0.6 0.6]);
end



