function model = read_anatomic_model(modelInfo)
% function model = read_anatomic_model(modelInfo)
% Reads anatomic model data from the input .mat file
%
% inputs:
%   modelInfo: Path to the anatomic model
%       modelInfo.folder = folder of anatomic model
%       modelInfo.name = name of anatomic model
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
% Last updated : 30-08-2019
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

filename = char(fullfile(modelInfo.folder, char(strcat(modelInfo.name, '.mat'))));
model = load(filename);
model.folder = modelInfo.folder;
model.name = modelInfo.name;
end