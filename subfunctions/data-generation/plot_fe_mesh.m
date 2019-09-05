function plot_fe_mesh(model)
% function plot_fe_mesh(model)
% Plots the finite element mesh for the selected anatomic model
%
% inputs:
%   model: Struct of anatomic model
%       model.header = information on anatomic model type
%       model.compartments = struct array of model compartments
%       model.variants =  struct array of model variants
%       model.folder = folder of anatomic model
%       model.name = name of anatomic model
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

no_fem_str = ['There is no finite element mesh in the input model.' newline ...
        'Try running the following code before calling plot_fe_mesh' newline ...
        'model = generate_fe_mesh(model)'];

%% Print error if model does not contain a finite element mesh
if ~isfield(model,'fem')
    error(no_fem_str);
end

figure;
axis equal;
hold on;
title('Finite Element Mesh')
map = [0.28 0.31 1
    1 1 0
    1 0.31 0.28
    0.28 1 0.31];
colormap(map);
h=plotmesh(model.fem.pos,[model.fem.tet model.fem.tissue],'x >600 & x < 1000 & z > 800 & z < 1300 & y > 125 & y < 150','LineWidth',0.5);
view([0 0]);

end