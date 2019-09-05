function  sensors = FECGSYN_read_sensors(filename) 
% sensors = FECGSYN_read_sensors(filename) 
% Reads in the selected .sensors file
%
% inputs:
%   filename: string - filename of input file
%
% outputs:
%   sensors: n x 2 - Sensor number and node index in maternal abdomen mesh
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

fid = fopen(filename,'r');

if (fid == -1)
    error('Unable to open the selected file.');
end

filetype = fgets(fid); 
if ~strcmpi(filetype(1:7), 'sensors')
    error('This filetype is not supported. Please select a .sources file.');    
end

numSensors = str2double(fgets(fid));

if  mod(numSensors,1) ~= 0
    error('Invalid number of sensors.');  
end

% Read in sensors data
% Sensor number and node index in maternal abdomen mesh
[sensorsData,readCount] = fscanf(fid,'%f %f\n', 2*numSensors);
if readCount~=2*numSensors
     error('Sensors file contains invalid data.');
end

sensors = reshape(sensorsData, 2, readCount/2)';
fclose(fid);

end