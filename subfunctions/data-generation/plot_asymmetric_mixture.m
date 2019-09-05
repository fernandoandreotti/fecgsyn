function plot_asymmetric_mixture(out, model)
% function plot_asymmetric_mixture(out, model)
% Plots the source and sensor positions for the asymmetric volume 
% conductor model and the generated maternal and fetal signals at each 
% sensor position according to the vectorcardiogram models in 'out'
%
% inputs:
%   out: outputs of the run_ecg_generator function containing the ecg
%        mixture and all the important model information that would allow
%        reproducing the results. 
%   model: Struct of anatomic model
%       model.header = information on anatomic model type
%       model.compartments = struct array of model compartments
%       model.variants =  struct array of model variants
%       model.folder = folder of anatomic model
%       model.name = name of anatomic model
%       model.sources = name of sources file
%       model.sensors = name of sensors file
%       model.fem = finite element mesh ready for solver
%       model.leadfield = struct array of maternal and fetal lead field matrices
%
%
% --
% fecgsyn toolbox, version 1.3-alpha, August 2019
% Released under the GNU General Public License
%
% Copyright (C) 2019  Emerson Keenan
% Department of Electrical and Electronic Engineering, University of
% Melbourne
% e.keenan@ieee.org
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Keenan E., Karmakar C K. and Palaniswami M., The effects of asymmetric volume conductor modeling on non-invasive fetal ECG extraction. 
% Physiol Meas 39(10), pp. 105013, 2018.
% 
%
% Last updated : 05-09-2019
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

conf = FECGSYN_check_platform();

%% Read in models, sources and sensors
[fetusPos, fetusTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.fetus.mesh.tri, '.off'));
[amnioticPos, amnioticTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.amniotic_fluid.mesh.tri, '.off'));
[abdoPos, abdoTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.maternal_abdomen.mesh.tri, '.off'));

[maternalSource, fetalSource] = FECGSYN_read_sources(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.sources, '.sources'));
sensors = FECGSYN_read_sensors(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.sensors, '.sensors'));
sources = struct('fetal',fetalSource,'maternal',maternalSource);
plotConfig = struct('sourcesOn',true,'sensorsOn',true);

%% Check if vernix is present
if isempty(model.compartments.vernix.mesh.tri) % No vernix present
    vertices = struct('fetus',fetusPos,'amnioticFluid',amnioticPos,'maternalAbdomen',abdoPos);
    faces = struct('fetus',fetusTri,'amnioticFluid',amnioticTri,'maternalAbdomen',abdoTri);
else % Vernix present
    [vernixPos, vernixTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_vernix_', model.compartments.vernix.mesh.tri, '.off'));
    vertices = struct('fetus',fetusPos,'amnioticFluid',amnioticPos,'maternalAbdomen',abdoPos,'vernix',vernixPos);
    faces = struct('fetus',fetusTri,'amnioticFluid',amnioticTri,'maternalAbdomen',abdoTri,'vernix',vernixTri);
end

if ~isstruct(vertices) || ~isfield(vertices,'fetus') || ...
   ~isfield(vertices,'amnioticFluid') || ...
   ~isfield(vertices, 'maternalAbdomen')
   error('Vertices struct is in incorrect format');
end

if ~isstruct(faces) || ~isfield(faces,'fetus') || ...
   ~isfield(faces,'amnioticFluid') || ...
   ~isfield(faces, 'maternalAbdomen')
   error('Faces struct is in incorrect format');
end

if ~isstruct(sources) || ~isfield(sources,'maternal') || ...
   ~isfield(sources,'fetal')
   error('Sources struct is in incorrect format');
end

fetalAlpha = 0.8;
fetalFaceColor = [0.28 0.31 1];
amnioticFluidAlpha = 0.3;
amnioticFluidFaceColor = [1 0.31 0.28];
vernixAlpha = 0.1;
vernixFaceColor = [1 0.95 0.1];
maternalAbdomenAlpha = 0.3;
maternalAbdomenFaceColor = [0.28 1 0.31];

f = figure('Name','Maternal-Fetal Model');
hold on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
lightangle(90,30);
lightangle(0,30);
view(0,0);
set(f, 'Position', [100, 100, 500, 400])
rotate3d;

fetusPatch = patch('Vertices', vertices.fetus, 'Faces', faces.fetus,... 
    'FaceColor', fetalFaceColor, 'facealpha', fetalAlpha,...
    'EdgeColor', 'none','SpecularColorReflectance', 0, ...
    'SpecularStrength', 0, 'AmbientStrength', 0.4);

amnioticFluidPatch = patch('Vertices', vertices.amnioticFluid, 'Faces', faces.amnioticFluid,... 
    'FaceColor', amnioticFluidFaceColor, 'facealpha', amnioticFluidAlpha,...
    'EdgeColor', 'none','SpecularColorReflectance', 0, ...
    'SpecularStrength', 0, 'AmbientStrength', 0.4);

if (~isempty(vertices.vernix) && ~isempty(faces.vernix))
    vernixPatch = patch('Vertices', vertices.vernix, 'Faces', faces.vernix,... 
        'FaceColor', vernixFaceColor, 'facealpha',vernixAlpha,...
        'EdgeColor', 'none','SpecularColorReflectance', 0, ...
        'SpecularStrength', 0, 'AmbientStrength', 0.4);
end

maternalPatch = patch('Vertices', vertices.maternalAbdomen, 'Faces', faces.maternalAbdomen,... 
    'FaceColor', maternalAbdomenFaceColor, 'facealpha', maternalAbdomenAlpha,...
    'EdgeColor', 'none','SpecularColorReflectance', 0, ...
    'SpecularStrength', 0, 'AmbientStrength', 0.4);

if plotConfig.sourcesOn && ~isempty(sources.maternal) && ~isempty(sources.fetal)
    maternalSourcePlot = scatter3([sources.maternal(1), sources.fetal(1)],...
        [sources.maternal(2), sources.fetal(2)],...
        [sources.maternal(3) sources.fetal(3)],...
        'filled','MarkerFaceColor',[0 0 0]);
    fetalSourcePlot = scatter3([sources.fetal(1)],...
        [sources.fetal(2)],...
        [sources.fetal(3)],...
        'filled','MarkerFaceColor',[1 0 0]);
end

if plotConfig.sensorsOn
    sensorPositions = vertices.maternalAbdomen(sensors(:,2),:);
    sensorsPlot = scatter3(sensorPositions(:,1),...
        sensorPositions(:,2), sensorPositions(:,3),...
        'MarkerEdgeColor',[0 0 0]);
    labels = 1:size(model.fem.sensors,1);
    labels = strtrim(cellstr([num2str(labels')])');
    text(sensorPositions(:,1),...
        sensorPositions(:,2), sensorPositions(:,3),...
        labels);
end

f_model = out.f_model{1,1};
m_model = out.m_model(1);
fVCG = f_model.VCG;
mVCG = m_model.VCG;

fetal_potentials = model.leadfield.fetal*fVCG;
maternal_potentials = model.leadfield.maternal*mVCG;

figure('Name','Fetal Source Potentials');
hold on;
title('Fetal Source Potentials (unscaled)');
plot(fetal_potentials');
legend(labels);

figure('Name','Maternal Source Potentials');
hold on;
title('Maternal Source Potentials (unscaled)');
plot(maternal_potentials');
legend(labels);

end
