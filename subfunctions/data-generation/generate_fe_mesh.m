function model = generate_fe_mesh(model)
% function model = generate_fe_mesh(model)
% Generate finite element mesh from input model
%
% inputs:
%   model: Struct of anatomic model
%       model.header = information on anatomic model type
%       model.compartments = struct array of model compartments
%       model.variants =  struct array of model variants
%       model.folder = folder of anatomic model
%       model.name = name of anatomic model
%       model.sources = name of sources file
%       model.sensors = name of sensors file
%
% outputs:
%   model: Struct of anatomic model
%       model.header = information on anatomic model type
%       model.compartments = struct array of model compartments
%       model.variants =  struct array of model variants
%       model.folder = folder of anatomic model
%       model.name = name of anatomic model
%       model.sources = name of sources file
%       model.sensors = name of sensors file
%       model.fem = finite element mesh ready for solver
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

if ispc
   clear tempdir
   setenv('TEMP',char(strcat(conf.fecgsynpath, 'data', conf.slashchar, 'temp', conf.slashchar, model.name, conf.slashchar))); 
elseif isunix
   clear tempdir
   setenv('TMP',char(strcat(conf.fecgsynpath, 'data', conf.slashchar, 'temp', conf.slashchar, model.name, conf.slashchar))); 
else
    error('Platform not supported');
end

%% Read in models, sources and sensors
[fetusPos, fetusTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.fetus.mesh.tri, '.off'));
[amnioticPos, amnioticTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.amniotic_fluid.mesh.tri, '.off'));
[abdoPos, abdoTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.compartments.maternal_abdomen.mesh.tri, '.off'));

[maternalSource, fetalSource] = FECGSYN_read_sources(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.sources, '.sources'));
sensors = FECGSYN_read_sensors(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_', model.sensors, '.sensors'));
sources = struct('fetal',fetalSource,'maternal',maternalSource);

%% Define conductivities for each region
% c(1) = Fetus
% c(2) = Vernix Caseosa
% c(3) = Amniotic Fluid
% c(4) = Maternal abdomen
c = [];
c(1) = model.compartments.fetus.conductivity;
c(2) = model.compartments.vernix.conductivity;
c(3) = model.compartments.amniotic_fluid.conductivity;
c(4) = model.compartments.maternal_abdomen.conductivity;

%% Define internal point for each region (columns 1-3) and maxvol (column 4)
pFetus = [model.compartments.fetus.mesh.innerpoint model.compartments.fetus.meshsize];
pAmniotic = [model.compartments.amniotic_fluid.mesh.innerpoint model.compartments.amniotic_fluid.meshsize];
pAbdo = [model.compartments.maternal_abdomen.mesh.innerpoint model.compartments.maternal_abdomen.meshsize];

%% Check if vernix is present
if isempty(model.compartments.vernix.mesh.tri) % No vernix present
    c = [c(1) c(3) c(3) c(4)];
    regions = [pFetus; pAmniotic; pAbdo];
    [fullV,fullF] = mergemesh(fetusPos,fetusTri,amnioticPos,amnioticTri,abdoPos,abdoTri); 
else % Vernix present
    [vernixPos, vernixTri] = FECGSYN_read_off(strcat(model.folder, conf.slashchar, model.name, conf.slashchar, model.name, '_vernix_', model.compartments.vernix.mesh.tri, '.off'));   
    pVernix = [model.compartments.vernix.mesh.innerpoint model.compartments.vernix.meshsize];
    regions = [pFetus; pVernix; pAmniotic; pAbdo];
    [fetusAndVernixPos,fetusAndVernixTri]= FECGSYN_join_nodes(fetusPos,fetusTri,vernixPos,vernixTri);                             
    [fullV,fullF]=mergemesh(fetusAndVernixPos,fetusAndVernixTri,amnioticPos,amnioticTri,abdoPos,abdoTri); 
end

%% Compute tetrahedralization using iso2mesh
[node,elem, ~]=surf2mesh(fullV,fullF,[],[],1,[],regions,[]);
faceRegions = elem(:,5);

%% Fix region labels for no vernix case
if isempty(model.compartments.vernix.mesh.tri) % No vernix present
    for i=1:length(faceRegions)
        if faceRegions(i) > 1
        faceRegions(i) = faceRegions(i)+1;
        end
    end
end

elem=meshreorient(node,elem(:,1:4)); % Ensure all elements are oriented
                                     % consistently   
elem = [elem faceRegions];         
abdoIndex = sensors(:,2);
abdoSensors = abdoPos(abdoIndex,:);

%% Create output finite element mesh structure
model.fem = struct();
model.fem.pos = node;
model.fem.tet = elem(:,1:4);
model.fem.tissuelabel = {'fetus','vernix','amniotic', 'abdomen'};
model.fem.tissue = faceRegions;
model.fem.sensors = abdoSensors;
model.fem.sources = sources;
model.fem.conductivity = c;

end