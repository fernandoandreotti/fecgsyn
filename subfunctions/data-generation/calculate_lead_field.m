function model = calculate_lead_field(model)
% function model = calculate_lead_field(model)
% Calculate lead field matrix from input model
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
%       model.fem = finite element mesh ready for solver
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

% Write variables to temporary storage
sensors = model.fem.sensors;
sources = model.fem.sources;
conductivity = model.fem.conductivity;
ft_sources = [sources.maternal; sources.fetal];

% Create a set of electrodes at chosen sensor positions
sens.elecpos = sensors;
sens.label = {};
nsens = size(sens.elecpos,1);
sens.unit = 'mm';

for ii=1:nsens
    sens.label{ii} = sprintf('vertex%03d', ii);
end

% Remove fields from model before passing to solver
model.fem = rmfield(model.fem,'sensors');
model.fem  = rmfield(model.fem,'sources');
model.fem  = rmfield(model.fem,'conductivity');

% Pass model to FieldTrip-Simbio solver
headmodel = ft_headmodel_simbio(model.fem,'conductivity',conductivity);
headmodel.unit = 'mm';
[headmodel, ~] = ft_prepare_vol_sens(headmodel,sens);
lf = ft_compute_leadfield(ft_sources,sens,headmodel);

% Add fields to model after returning from solver
model.fem.sensors = sensors;
model.fem.sources = sources;
model.fem.conductivity = conductivity;
model.leadfield = struct('maternal',lf(:,1:3),'fetal',lf(:,4:6));

end